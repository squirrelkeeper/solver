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
	full_name << batch_no;
	full_name << ".dat";


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



vector<double> timeseries::pulse_positions()
{
	vector<double> pp_out;
		
	auto max_iterator = max_element(I.begin(), I.end());
	I_max = I[distance(I.begin(), max_iterator)];
	
	auto min_iterator = min_element(I.begin(), I.end());
	I_min = I[distance(I.begin(), min_iterator)];
	
	I_mid = (I_max - I_min)*0.5;
	
	I_mean = accumulate(I.begin(), I.end(), (double)(0.0));
	I_mean/= (double)(len);

	double thres_front = (I_mid - I_mean)*0.5;
	double thres_tail = (I_mid - I_mean)*0.25;
	
	double sqrtD = sqrt(AP->IP.D.par_dbl);
	bool pulse_detected = false;
	int pulse_counter = 0;	
	
	vector<long> pulse_left_index;
	vector<long> pulse_right_index;
	
	double front_t;
	double tail_t;
	
	
	
	for(long i = 0; i < len; i++)
	{
		if(I[i] > thres_front && pulse_detected==false)
		{
			pulse_detected = true;
			pulse_left_index.push_back(i);
			front_t = t[i];
			cout << t[i] << '\t';
		}
		if(I[i] < thres_tail && pulse_detected==true)
		{
			pulse_counter++;
			pulse_detected = false;
			pulse_right_index.push_back(i);
			tail_t = t[i];
			cout << t[i] << '\t';
			cout << tail_t-front_t << endl;
		}
	}
	
	cout << pulse_counter << endl;
	
	cout << I_max << '\t';
	cout << I_min << '\t';
	cout << I_mid << '\t';
	cout << I_mean << endl;
	

		

	
	return pp_out;
}



vector<double> timeseries::pulse_positions2()
{
	vector<double> pp_out;
		
	auto max_iterator = max_element(I.begin(), I.end());
	I_max = I[distance(I.begin(), max_iterator)];
	
	auto min_iterator = min_element(I.begin(), I.end());
	I_min = I[distance(I.begin(), min_iterator)];
	
	I_mid = (I_max - I_min)*0.5;
	
	I_mean = accumulate(I.begin(), I.end(), (double)(0.0));
	I_mean/= (double)(len);

	double thres_front = (I_mid - I_mean)*0.5;
	double thres_tail = (I_mid - I_mean)*0.25;
	
	double sqrtD = sqrt(AP->IP.D.par_dbl);
	bool pulse_detected = false;
	int pulse_counter = 0;	
	
	vector<long> pulse_left_index;
	vector<long> pulse_right_index;
	
	double front_t;
	double tail_t;
	
	
	
	for(long i = 0; i < len; i++)
	{
		if(I[i] > thres_front && pulse_detected==false)
		{
			pulse_detected = true;
			pulse_left_index.push_back(i);
			front_t = t[i];
			cout << t[i] << '\t';
		}
		if(I[i] < thres_tail && pulse_detected==true)
		{
			pulse_counter++;
			pulse_detected = false;
			pulse_right_index.push_back(i);
			tail_t = t[i];
			cout << t[i] << '\t';
			cout << tail_t-front_t << endl;
		}
	}
	
	cout << pulse_counter << endl;
	
	cout << I_max << '\t';
	cout << I_min << '\t';
	cout << I_mid << '\t';
	cout << I_mean << endl;
	

		

	
	return pp_out;
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





