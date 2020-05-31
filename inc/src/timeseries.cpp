#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <numeric>
#include <iomanip>




#include "../hdr/timeseries.hpp"

using namespace std;

pulse::pulse()
{
	trig_pos = 0;
	right_pos = 0;
	left_pos = 0;
	baseline = 0;
	max = 0;
	pos = 0;
	width = 0;
	del = false;
}

void pulse::reset()
{
	trig_pos = 0;
	right_pos = 0;
	left_pos = 0;
	baseline = 0;
	max = 0;
	pos = 0;
	width = 0;
	del = false;
}

//###########################################

timeseries::timeseries()
{
}

timeseries::timeseries(allpar_set *init_AP)
{
	AP = init_AP;
	ipar_dbl_set ip(*AP);

	len_dbl = ip.out_time;	
	len = len_dbl/ip.dt;

	t.resize(len);
	I.resize(len);
	X.resize(len);
}

timeseries::timeseries(double init_len, allpar_set *init_AP)
{
	*AP = *init_AP;
	ipar_dbl_set ip(*AP);

	len_dbl = init_len;	
	len = len_dbl/ip.dt;
	
	t.resize(len);
	I.resize(len);
	X.resize(len);
}

void timeseries::write_file(string file_name)
{
	stringstream full_name;
	full_name.precision(7);
	
	
	srand (time(NULL));
	int batch_no = rand()%10000;

	full_name << "data/";
	full_name << file_name;
	full_name << "_D";
	full_name << AP->IP.D.par_dbl;
//	full_name << "_ba";
//	full_name << batch_no;
	full_name << ".ts.dat";


	ofstream data;
	data.open(full_name.str().c_str(),ios::trunc);


	data << "#t" << '\t';
	data << "ER" << '\t';
	data << "EI" << '\t';
	data << "G" << '\t';
	data << "Q" << '\t';
	data << "J" << '\t';
	data << "I" << '\t';
	data << endl;
	
	
	
	for(long i = 0; i < len; i++)
	{
		data << setprecision(15);

		data << t[i] << '\t';
		data << X[i].ER << '\t';
		data << X[i].EI << '\t';
		data << X[i].G << '\t';
		data << X[i].Q << '\t';
		data << X[i].J << '\t';
		data << I[i] << '\t';

		data << endl;
	}
	
	data.close();
}




vector<pulse> timeseries::pulse_analysis()
{
	vector<pulse> pulse_list;
	
	auto max_it = max_element(I.begin(), I.end());
	I_max = I[distance(I.begin(), max_it)];

	double thres = 0.75 * AP->IP.D.par_dbl * I_max;
	
	double thres_left = 0.75 * AP->IP.D.par_dbl * I_max;
	double thres_right = 0.75 * AP->IP.D.par_dbl * I_max;
	
	bool pulse_detected = false;
	long unit = ceil(AP->LP.T.par_dbl/(1000.0*AP->IP.dt.par_dbl));
	long whs = 5*unit;
	
	long baseline_gap = 2 * unit;
	long baseline_interval = 4 * unit;
	long baseline_space = baseline_interval + baseline_gap;
	
	double Iw_mean = 0.0;
	
	bool skip_pulse = false;

	pulse curr_pulse;
	
	cout << I_max << endl;
	cout << unit << endl;
	cout << whs << endl;
	cout << thres << endl;

	for(long i = whs; i < len - whs; i++)
	{
	
		//sliding window average
		for(unsigned int j = i - whs; j < i + whs; j++)
		{
			Iw_mean += I[j];
		}
		
		Iw_mean /= (double)(2*whs);
		
		//trigger event
		if(Iw_mean >= thres && pulse_detected==false && i - baseline_space > 0)
		{
//			cout << "trigger at: " << t[i] << " / " << i << endl;
			pulse_detected = true;
			curr_pulse.trig_pos = i;
			

			for(int j = i - baseline_space; j < i - baseline_gap; j++)
			{
				//alert if j_start < 0
				curr_pulse.baseline += I[j];
			}
			
			curr_pulse.baseline /= (double)(baseline_interval);
		}
		
		//find right edge and expand left edge
		if(Iw_mean <= curr_pulse.baseline && pulse_detected==true)
		{
			curr_pulse.right_pos = i;
			
			long k = curr_pulse.trig_pos;
			while(k >= 0)
			{
				if(k == 0)
				{
					skip_pulse = true;
					break;
				}
				else if(I[k] <= curr_pulse.baseline)
				{
					curr_pulse.left_pos = k;
					break;
				}
				
				k--;
			}
			
			if(skip_pulse == false)
			{
				curr_pulse.width = t[curr_pulse.right_pos]-t[curr_pulse.left_pos];

				double I_pulse_sum = 0.0;
				
				for(long l = curr_pulse.left_pos; l < curr_pulse.right_pos; l++)
				{
					if(I[l] > curr_pulse.max)
					{
						curr_pulse.max = I[l];
					}
					
					I_pulse_sum += I[l];
					curr_pulse.pos += I[l] * t[l];
				}
				
				curr_pulse.pos /= I_pulse_sum;
				
				pulse_list.push_back(curr_pulse);
				curr_pulse.reset();
				
				

				pulse_detected = false;
			}
			else
			{
				i = curr_pulse.right_pos;
				curr_pulse.reset();
				pulse_detected = false;
			}
		}
		
		
	}
	
	//erase all nan valued pulses -- this doesnt work
	for(unsigned int i = 0; i < pulse_list.size(); i++)
	{
		if(isnan(pulse_list[i].pos))
		{
			pulse_list.erase(pulse_list.begin()+i);
		}
	}
	
	
	
	//scan for small pulse distances -- but what is small?
	if(pulse_list.size() > 1)
	{
		double dist_mean = 0.0;
		for(unsigned int i = 1; i < pulse_list.size(); i++)
		{
			dist_mean += pulse_list[i].pos - pulse_list[i-1].pos;
			cout << pulse_list[i].pos << endl;
//			cout << dist_mean << endl;
		}
		cout << (double)(pulse_list.size()-1) << endl;
		dist_mean /= (double)(pulse_list.size()-1);
		
		cout << dist_mean << endl;

		
		for(unsigned int i = 1; i < pulse_list.size(); i++)
		{
			double dist = pulse_list[i].pos - pulse_list[i-1].pos;
			if(dist < 0.9 * dist_mean)
			{
				cout << "small" << endl;
				int index = pulse_list[i].max < pulse_list[i-1].max ? i : i-1;
				cout << pulse_list[index].pos << endl;
			}
		}
	}
	
	//erase all marked pulses
	for(unsigned int i = 0; i < pulse_list.size(); i++)
	{
		if(pulse_list[i].del)
		{
			pulse_list.erase(pulse_list.begin()+i);
		}
	}

	return pulse_list;
}


void timeseries::cout_pulse_data(vector<pulse> pulse_list)
{

	cout << "counted: " << pulse_list.size() << endl;

	cout << "trig_pos" << '\t';
	cout << "right_pos" << '\t';
	cout << "left_pos" << '\t';
	cout << "baseline" << '\t';
	cout << "max" << '\t';
	cout << "pos" << '\t';
	cout << "width" << '\t';
	cout << "del" << '\t';
	cout << endl;


	for(unsigned int i = 0; i < pulse_list.size(); i++)
	{
		cout << pulse_list[i].trig_pos << '\t';
		cout << pulse_list[i].right_pos << '\t';
		cout << pulse_list[i].left_pos << '\t';
		cout << pulse_list[i].baseline << '\t';
		cout << pulse_list[i].max << '\t';
		cout << pulse_list[i].pos << '\t';
		cout << pulse_list[i].width << '\t';
		cout << pulse_list[i].del << '\t';
		cout << endl;
	}
}

vector<double> timeseries::get_pulse_dist(vector<pulse> pulse_list)
{
	vector<double> pulse_dist;
	
	for(unsigned i = 1; i < pulse_list.size(); i++)
	{
		pulse_dist.push_back(pulse_list[i].pos - pulse_list[i-1].pos);
	}
	
	return pulse_dist;
}













