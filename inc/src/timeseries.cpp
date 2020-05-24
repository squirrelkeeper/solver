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
	thres_pos = 0;
	right_pos = 0;
	left_pos = 0;
	baseline = 0;
	max = 0;
	pos = 0;
}

void pulse::reset()
{
	thres_pos = 0;
	right_pos = 0;
	left_pos = 0;
	baseline = 0;
	max = 0;
	pos = 0;
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
	bool pulse_detected = false;
	long unit = ceil(AP->LP.T.par_dbl/(1000.0*AP->IP.dt.par_dbl));
	long whs = 5*unit;
	
	double Iw_mean = 0.0;

	pulse curr_pulse;
	
	for(long i = whs; i < len - whs; i++)
	{
		
		for(unsigned int j = i - whs; j < i + whs; j++)
		{
			Iw_mean += I[j];
		}
		
		Iw_mean /= (double)(2*whs);
		
		
		if(Iw_mean >= thres && pulse_detected==false && i - 6 * unit > 0)
		{
			pulse_detected = true;
			
//			cout << "pulse detected at: " << t[i] << " / " << i << endl;
			
			curr_pulse.thres_pos = i;
			

			for(int j = i - 6 * unit; j < i - 2*unit; j++)
			{
				//alert if j_start < 0
				curr_pulse.baseline += I[j];
			}
			
			curr_pulse.baseline /= (double)(4*unit);
		}
		else if(Iw_mean >= thres && pulse_detected==false && i - 6 * unit < 0)
		{
			i += 1000*unit;
		}
		
		if(Iw_mean <= curr_pulse.baseline && pulse_detected==true)
		{
			curr_pulse.right_pos = i;
			
			
			long j = curr_pulse.thres_pos;
			while(j >= 0)
			{
				if(j == 0)
				{
					pulse_detected = false;
					curr_pulse.reset();
					i += 1000*unit;
					break;
				}
				else if(I[j] <= curr_pulse.baseline)
				{
					curr_pulse.left_pos = j;
					break;
				}
				
				j--;
			}
			
			double I_pulse_sum = 0.0;
			
			for(long j = curr_pulse.left_pos; j < curr_pulse.right_pos; j++)
			{
				if(I[j] > curr_pulse.max)
				{
					curr_pulse.max = I[j];
				}
				
				I_pulse_sum += I[j];
				curr_pulse.pos += I[j] * t[i];
			}
			
			curr_pulse.pos /= I_pulse_sum;
			
			pulse_list.push_back(curr_pulse);
			
			i = curr_pulse.left_pos;

			curr_pulse.reset();
			pulse_detected = false;
			cout << "done" << endl;
		}
		
		if(i > 7000) break;
		
	}
	

	return pulse_list;
}


void timeseries::cout_pulse_data(vector<pulse> pulse_list)
{
	cout << "counted: " << pulse_list.size() << endl;
	
	cout << "thres_pos" << '\t';
	cout << "right_pos" << '\t';
	cout << "left_pos" << '\t';
	cout << "baseline" << '\t';
	cout << "max" << '\t';
	cout << "pos" << '\t';
	cout << endl;
	
	for(unsigned int i = 0; i < pulse_list.size(); i++)
	{
		cout << pulse_list[i].thres_pos << '\t';
		cout << pulse_list[i].right_pos << '\t';
		cout << pulse_list[i].left_pos << '\t';
		cout << pulse_list[i].baseline << '\t';
		cout << pulse_list[i].max << '\t';
		cout << pulse_list[i].pos << '\t';
		cout << endl;
	}
}



/*

for(unsigned int i=1; i<res_tt.size();i++)
{
    if(res_Eabs[i-1]<res_Eabs[i]&&res_Eabs[i+1]<res_Eabs[i])
    {
        vector<double> temp;
        temp.push_back(res_tt[i]);
        temp.push_back(res_Eabs[i]);
        res_maxima.push_back(temp);
    }
}

for(unsigned int i=0; i<=res_Eabs.size();i++)
{
    if(i==0)
    {
        absolute_max = res_Eabs[0];
        absolute_min = res_Eabs[0];
        absolute_mean = res_Eabs[0];
    }
    if(absolute_max < res_Eabs[i]) absolute_max = res_Eabs[i];
    if(absolute_min > res_Eabs[i]) absolute_min = res_Eabs[i];
    absolute_mean += res_Eabs[i];
    if(i==res_Eabs.size()) absolute_mean/=i;
}


for(unsigned int i=0; i<res_tt.size();i++)
{
    if(res_Eabs[i]>0.1*absolute_max+Lpara->D*1.1 && pulse_detected==false)
    {
        pulse_detected = true;
        pulse_left_index.push_back(i);
    }
    if(res_Eabs[i]<0.1*absolute_max && pulse_detected==true)
    {
        pulse_counter++;
        pulse_detected = false;
        pulse_right_index.push_back(i);
    }
}


if(pulse_counter != 0)

{    

if(pulse_left_index.size()!=pulse_right_index.size())
{
    pulse_left_index.pop_back();
}

for(unsigned int i=0;i<pulse_left_index.size();i++)
{
	if(pulse_left_index[i] != 0)
	{
		for(unsigned int j=pulse_left_index[i];res_Eabs[j]>=res_Eabs[pulse_right_index[i]];j--)
		{
			if(pulse_left_index[i] == 0)
			{
				break;
			}
		}
	}
}



for(unsigned int i=0;i<pulse_left_index.size();i++)
{
    double weighted_amp_sum = 0;
    double amp_sum = 0;
    for(unsigned int j=pulse_left_index[i]; j<pulse_right_index[i]; j++)
    {
        weighted_amp_sum += res_Eabs[j] * res_tt[j];
        amp_sum += res_Eabs[j];
    }
    pulse_postions.push_back(weighted_amp_sum/amp_sum);

    out_pulse_postions.push_back(Lpara->Jg);
    out_pulse_postions.push_back(Lpara->gg);

    out_pulse_postions.push_back(Lpara->gq);
    out_pulse_postions.push_back(Lpara->q0);
    out_pulse_postions.push_back(FBpara->K);
    out_pulse_postions.push_back(FBpara->tau);
    out_pulse_postions.push_back(FBpara->tauLP);
    out_pulse_postions.push_back(r);
    out_pulse_postions.push_back(weighted_amp_sum/amp_sum);

}

for(unsigned int i=0;i<pulse_postions.size()-1;i++)
{
    t_isi.push_back(pulse_postions[i+1]-pulse_postions[i]);
//    cout << "Ppos" << '\t' << pulse_postions[i] << endl;
//    cout << "Ppos+1" << '\t' << pulse_postions[i+1] << endl;
//    cout << "t_isi[i]" << '\t' << t_isi[i] << endl;
}



*/





