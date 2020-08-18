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
	max_pos = 0;

	right_val = 0.0;
	left_val = 0.0;
	max_val = 0.0;


	baseline = 0.0;

	pos = 0.0;
	width = 0.0;
	fwhm = 0.0;
	area = 0.0;


	del = false;
}

pulse::pulse(long init_trig_pos)
{
	trig_pos = init_trig_pos;
	right_pos = 0;
	left_pos = 0;
	max_pos = 0;

	right_val = 0.0;
	left_val = 0.0;
	max_val = 0.0;


	baseline = 0.0;

	pos = 0.0;
	width = 0.0;
	fwhm = 0.0;
	area = 0.0;


	del = false;
}

void pulse::reset()
{
	trig_pos = 0;
	right_pos = 0;
	left_pos = 0;
	max_pos = 0;

	right_val = 0.0;
	left_val = 0.0;
	max_val = 0.0;


	baseline = 0.0;

	pos = 0.0;
	width = 0.0;
	fwhm = 0.0;
	area = 0.0;


	del = false;
}

//###########################################

timeseries::timeseries()
{
}

timeseries::timeseries(allpar_set init_AP)
{
	AP = init_AP;
	ipar_dbl_set ip(AP);

	len_dbl = ip.out_time;
	len = len_dbl/ip.dt;

	t.resize(len);
	I.resize(len);
	X.resize(len);
}

timeseries::timeseries(double init_len, allpar_set init_AP)
{
	AP = init_AP;
	ipar_dbl_set ip(AP);

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
//	int batch_no = rand()%10000;

	full_name << "data/";
	full_name << file_name;
	full_name << "_D";
	full_name << AP.IP.D.par_dbl;
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
	data << "I";
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
		data << I[i];

		data << endl;
	}
	
	data.close();
}




void timeseries::cut_series(string opt, double time_interval)
{
	long new_len = time_interval/AP.IP.dt.par_dbl;

	vector<double> new_t(new_len);
	vector<double> new_I(new_len);
	vector<var> new_X(new_len);
	
	long start;
	
	if(opt=="last")
	{
		start = len - new_len;
	}
	else if(opt=="first")
	{
		start = 0;
	}
	
	for(long i = 0; i < new_len; i++)
	{
		new_t[i] = t[i+start]; 
		new_I[i] = I[i+start];
		new_X[i] = X[i+start];
	}
	
	len = new_len;
	len_dbl = time_interval;
	t = new_t;
	I = new_I;
	X = new_X;
}

void timeseries::reverse_series()
{
	vector<double> new_t(len);
	vector<double> new_I(len);
	vector<var> new_X(len);
	
	for(long i = 0; i < len; i++)
	{
		new_t[i] = t[len-1 - i]; 
		new_I[i] = I[len-1 - i];
		new_X[i] = X[len-1 - i];
	}
	
	t = new_t;
	I = new_I;
	X = new_X;
}

void timeseries::cc_series()
{
	for(long i = 0; i < len; i++)
	{
		X[i].EI = -X[i].EI;
	}
}

void timeseries::reset_time()
{
	double dt = AP.IP.dt.par_dbl;
	
	for(long i = 0; i < len; i++)
	{
		t[i] = (double)(i) * dt;
	}
}



//###########################################



ts_evaluation::ts_evaluation(vector<double> *I_ptr_init, vector<double> *t_ptr_init, allpar_set AP_init)
{
	I_ptr = I_ptr_init;
	t_ptr = t_ptr_init;
	len = (long)(*I_ptr).size();
	t_offset = t_ptr->at(0);
	
	AP = AP_init;
	dt = AP.IP.dt.par_dbl;
	
	MaxPos.resize(0);
	MaxVal.resize(0);
	MinPos.resize(0);
	MinVal.resize(0);
	UniqMaxPos.resize(0);
	UniqMaxVal.resize(0);
	UniqMinPos.resize(0);
	UniqMinVal.resize(0);
	
	average = 0;
	period = 0;

	GlobalSupr = I_ptr->at(0);
	GlobalInfi = I_ptr->at(0);
}


ts_evaluation::ts_evaluation(timeseries *TS_init, allpar_set AP_init)
{
	I_ptr = &TS_init->I;
	t_ptr = &TS_init->t;
	len = TS_init->len;
	t_offset = t_ptr->at(0);
	
	AP = AP_init;
	dt = AP.IP.dt.par_dbl;
	
	MaxPos.resize(0);
	MaxVal.resize(0);
	MinPos.resize(0);
	MinVal.resize(0);
	UniqMaxPos.resize(0);
	UniqMaxVal.resize(0);
	UniqMinPos.resize(0);
	UniqMinVal.resize(0);
	
	average = 0;
	period = 0;

	GlobalSupr = I_ptr->at(0);
	GlobalInfi = I_ptr->at(0);
}





void ts_evaluation::FindUniqMax()
{
	int max_cnt;
	
	if(GlobalSupr/GlobalInfi < (1.0 + MaxMinDscr) || average < AverageThres)
	{
		max_cnt = 0;
	}
	else
	{
		max_cnt = 0;
		bool doublecount;
		
		for(unsigned i = 0; i < MaxVal.size(); i++)
		{
			doublecount = false;
			
			for(int j = 0; j < max_cnt; j++)
			{
				if(abs(MaxVal[i]-average) > (1.0-DblCountTol) * abs(UniqMaxVal[j]-average) && abs(MaxVal[i]-average) < (1.0+DblCountTol) * abs(UniqMaxVal[j]-average))
				{
					doublecount = true;
				}
			}
			if(doublecount == false)
			{
				UniqMaxPos.push_back(MaxPos[i]);
				UniqMaxVal.push_back(MaxVal[i]);
				max_cnt += 1;
			}
		}
	}
	
	if(max_cnt > 1)
	{
		av_max_dist = (MaxPos[0]-MaxPos[MaxPos.size()-1])/((double)(max_cnt-1));
	}
	
	uniq_max_num = max_cnt;
}

void ts_evaluation::FindPeriod()
{
	vector<double> MaxPattern = {};
	int seq_num = 0;
	int per_num = 0;

	
	if(MaxVal.size()==0)
	{
		per_num = -1.0;
	}
	else
	{
		for(unsigned i = 0; i < MaxVal.size(); i++)
		{
			MaxPattern.push_back(MaxVal[i]);
			bool con;
			
			for(unsigned j = 0; j < MaxVal.size(); j++)
			{
				con = true;
				if( !is_approx(MaxPattern[j%MaxPattern.size()], MaxVal[j], "rel", DblCountTol) )
				{
					seq_num = 0;
					con = false;
					break;
				}
				if(j%MaxPattern.size() == MaxPattern.size()-1)
				{
					seq_num++;
				}
			}
			if(con)
			{
				per_num = i+1;
				break;
			}
		}
	}
	
	if(seq_num > 1)
	{
		period = (MaxPos[seq_num*per_num-1]-MaxPos[0])/((double)seq_num);
	}
	else
	{
		period = -1.0;
	}
}


void ts_evaluation::FindState()
{
	double T0 = AP.LP.T.par_dbl; 
	
	//test for transients
	if(uniq_max_num == 0 && average < AverageThres)
	{
		//off
		state = 0;
	}
	else if(uniq_max_num == 0 && average > AverageThres)
	{
		//cw
		state = 1000;
	}
	else if(uniq_max_num == 1 && T0/av_max_dist > 0.5)
	{
		//FML + HML
		state = (int)(round(T0/period));
	}
	else if(uniq_max_num > 1 && uniq_max_num <= 15 && T0/av_max_dist > 1.5 && T0/period < 1.4)
	{
		//SP - satelite pulses
		state = 99 + uniq_max_num;
	}
	else if(uniq_max_num > 1 && T0/av_max_dist <= 1.4 && T0/av_max_dist > 0.05)
	{
		//QS ML
		state = 200 + uniq_max_num;
	}
	else if(uniq_max_num >= 1 && T0/av_max_dist < 0.05)
	{
		//QS
		state = 1001;
	}
	else if(uniq_max_num >= 15 && T0/av_max_dist > 1.5)
	{
		//QP - chaotic
		state = 400 + uniq_max_num;
	}
	else
	{
		//error
		state =  -1;
	}
}




bool ts_evaluation::is_approx(double x, double y, string opt, double val)
{
	bool con = false;

	if(opt == "eq")
	{
		con = (x == y);
	}
	else if(opt == "abs")
	{
		con = (x >= y-val && x <= y+val);
	}
	else if(opt == "rel")
	{
		con = (x >= y*(1-val) && x <= y*(1+val));
	}

	return con;
}






double ts_evaluation::InterpolQuadExtrPos(double a, double b, double c, double fa, double fb, double fc)
{
	return b + (-0.5*(-a + b)*(-a + b)*(-fb + fc) + 0.5*(-b + c)*(-b + c)*(fa - fb))/((-a + b)*(-fb + fc) + (-b + c)*(fa - fb));
}


double ts_evaluation::InterpolQuadExtrVal(double a, double b, double c, double fa, double fb, double fc, double tmax)
{
	return (fa*(b - c)*(b - tmax)*(c - tmax) - fb*(a - c)*(a - tmax)*(c - tmax) + fc*(a - b)*(a - tmax)*(b - tmax))/((a - b)*(a - c)*(b - c));
}





//###########################################





pp_evaluation::pp_evaluation(vector<double> *I_ptr_init, vector<double> *t_ptr_init, allpar_set AP_init)
{
	I_ptr = I_ptr_init;
	t_ptr = t_ptr_init;
	len = (long)(*I_ptr).size();
	t_offset = t_ptr->at(0);
	
	AP = AP_init;
	dt = AP.IP.dt.par_dbl;

	average = 0;

	GlobalSupr = I_ptr->at(0);
	GlobalInfi = I_ptr->at(0);

	pulse_list.resize(0);
	pulse_list_len = 0;

}


pp_evaluation::pp_evaluation(timeseries *TS_init, allpar_set AP_init)
{
	I_ptr = &TS_init->I;
	t_ptr = &TS_init->t;
	len = TS_init->len;
	t_offset = t_ptr->at(0);
	
	AP = AP_init;
	dt = AP.IP.dt.par_dbl;
	
	average = 0;

	GlobalSupr = I_ptr->at(0);
	GlobalInfi = I_ptr->at(0);

	pulse_list.resize(0);
	pulse_list_len = 0;
}

void pp_evaluation::FindGlobalExtrAvrg()
{
	average = 0.0;
	
	for(long i = 0; i < len; i++)
	{
		average += (*I_ptr)[i];
		
		if((*I_ptr)[i] > GlobalSupr)
		{
			GlobalSupr = (*I_ptr)[i];
		}
		if((*I_ptr)[i] < GlobalInfi)
		{
			GlobalInfi = (*I_ptr)[i];
		}
	}

	average /= (double)(len);
}

void pp_evaluation::DetectPulses_simple()
{
	double RelTrigThres =  0.2;
	double AbsTrigThres = RelTrigThres * GlobalSupr + 1.1 * AP.IP.D.par_dbl;
	
	double AbsResetThres;
	
	if(RelTrigThres * GlobalSupr < average)
	{
		AbsResetThres = RelTrigThres * GlobalSupr;
	}
	else if(RelTrigThres * GlobalSupr >= average)
	{
		AbsResetThres = average;
	}
	
	
	bool triggered = false;
	int curr_pulse = -1;
		
	for(long i = 0; i < len; i++)
	{
		if(!triggered && (*I_ptr)[i] >= AbsTrigThres)
		{
			triggered = true;

			pulse P(i);
			pulse_list.push_back(P);
			pulse_list_len++;
			curr_pulse++;
		}

		if(triggered)
		{
			if((*I_ptr)[i] >= pulse_list[curr_pulse].max_val)
			{
				pulse_list[curr_pulse].max_pos = i;
				pulse_list[curr_pulse].max_val = (*I_ptr)[i];
			}
		}
		
		if(triggered && (*I_ptr)[i] < AbsResetThres)
		{
			triggered = false;
			
			pulse_list[curr_pulse].right_pos = i;
			pulse_list[curr_pulse].right_val = (*I_ptr)[i];
			
			
			for(long j = pulse_list[curr_pulse].trig_pos; j >= 0; j--)
			{	if((*I_ptr)[j] <= AbsResetThres)
				{
					pulse_list[curr_pulse].left_pos = j;
					pulse_list[curr_pulse].left_val = (*I_ptr)[j];
					break;
				}
				else if(j == 0)
				{
					pulse_list[curr_pulse].del = true;
					break;
				}
				
				
			}
		}
	}
	
	if(triggered)
	{
		pulse_list[curr_pulse].del = true;
	}
	
	pulse_list = PurgeList(pulse_list);
	pulse_list_len = pulse_list.size();
	
}



void pp_evaluation::DetectPulses_MWA()
{
	double RelTrigThres =  0.2;
	double AbsTrigThres = RelTrigThres * GlobalSupr + 1.1 * AP.IP.D.par_dbl;
	
	double AbsResetThres;
	
	if(RelTrigThres * GlobalSupr < average)
	{
		AbsResetThres = RelTrigThres * GlobalSupr;
	}
	else if(RelTrigThres * GlobalSupr >= average)
	{
		AbsResetThres = average;
	}
	
	
	bool triggered = false;
	int curr_pulse = -1;
	long u = (long)(1.0/(100.0 * dt));	
	
	for(long i = u; i < len; i++)
	{
		double mwa = 0.0;
		
		for(long j = i - u; j < i + u; j++)
		{
			mwa += (*I_ptr)[j];
		}
		
		mwa /= 2.0*u;
		
		if(!triggered && mwa >= AbsTrigThres)
		{
			triggered = true;

			pulse P(i);
			pulse_list.push_back(P);
			pulse_list_len++;
			curr_pulse++;
 		}

		if(triggered)
		{			
			if((*I_ptr)[i] >= pulse_list[curr_pulse].max_val)
			{
				pulse_list[curr_pulse].max_pos = i;
				pulse_list[curr_pulse].max_val = (*I_ptr)[i];
			}
		}
		
		if(triggered && (*I_ptr)[i] < AbsResetThres)
		{
			triggered = false;
			
			pulse_list[curr_pulse].right_pos = i;
			pulse_list[curr_pulse].right_val = (*I_ptr)[i];
			
			
			for(long j = pulse_list[curr_pulse].trig_pos; j >= 0; j--)
			{	if((*I_ptr)[j] <= AbsResetThres)
				{
					pulse_list[curr_pulse].left_pos = j;
					pulse_list[curr_pulse].left_val = (*I_ptr)[j];
					break;
				}
				else if(j == 0)
				{
					pulse_list[curr_pulse].del = true;
					break;
				}
				
				
			}
		}
	}
	
	if(triggered)
	{
		pulse_list[curr_pulse].del = true;
	}
	
	pulse_list = PurgeList(pulse_list);
	pulse_list_len = pulse_list.size();
}


void pp_evaluation::DetectPulses_MWC()
{
	double RelTrigThres =  0.2;
	double AbsTrigThres = RelTrigThres * GlobalSupr + 1.1 * AP.IP.D.par_dbl;
	
	double AbsResetThres;
	
	if(RelTrigThres * GlobalSupr < average)
	{
		AbsResetThres = RelTrigThres * GlobalSupr;
	}
	else if(RelTrigThres * GlobalSupr >= average)
	{
		AbsResetThres = average;
	}
	
	
	bool triggered = false;
	int curr_pulse = -1;
	long u = (long)(1.0/(100.0 * dt));	
	
	for(long i = u; i < len; i++)
	{
		double mwc1 = 0.0;
		double mwc2 = 0.0;
		
		
		for(long j = i - u; j < i + u; j++)
		{
			mwc1 += (*I_ptr)[j];
		}
		
		mwc1 /= 2.0*u;
		
		for(long j = i - u; j < i + u; j++)
		{
			mwc2 += ( (*I_ptr)[j] - mwc1 ) * ( (*I_ptr)[j] - mwc1 );
		}
		
		mwc2 /= 2.0*u;
		
		if(!triggered && abs(mwc1 * mwc2) >= AbsTrigThres)
		{
			triggered = true;

			pulse P(i);
			pulse_list.push_back(P);
			pulse_list_len++;
			curr_pulse++;
		}

		if(triggered)
		{			
			if((*I_ptr)[i] >= pulse_list[curr_pulse].max_val)
			{
				pulse_list[curr_pulse].max_pos = i;
				pulse_list[curr_pulse].max_val = (*I_ptr)[i];
			}
		}
		
		if(triggered && (*I_ptr)[i] < AbsResetThres)
		{
			triggered = false;
			
			pulse_list[curr_pulse].right_pos = i;
			pulse_list[curr_pulse].right_val = (*I_ptr)[i];
			
			
			for(long j = pulse_list[curr_pulse].trig_pos; j >= 0; j--)
			{	if((*I_ptr)[j] <= AbsResetThres)
				{
					pulse_list[curr_pulse].left_pos = j;
					pulse_list[curr_pulse].left_val = (*I_ptr)[j];
					break;
				}
				else if(j == 0)
				{
					pulse_list[curr_pulse].del = true;
					break;
				}
				
				
			}
		}
	}
	
	if(triggered)
	{
		pulse_list[curr_pulse].del = true;
	}
	
	pulse_list = PurgeList(pulse_list);
	pulse_list_len = pulse_list.size();
}



void pp_evaluation::FilterSamePos()
{
	for(int i = 1; i < pulse_list_len; i++)
	{
		bool con1 = (pulse_list[i].trig_pos == pulse_list[i-1].trig_pos);
		bool con2 = (pulse_list[i].left_pos == pulse_list[i-1].left_pos);
		bool con3 = (pulse_list[i].right_pos == pulse_list[i-1].right_pos);
		bool con4 = (pulse_list[i].max_pos == pulse_list[i-1].max_pos);
		
		
		if(con1 || con2 || con3 || con4)
		{
			pulse_list[i].del = true;
		}
	}
	
	pulse_list = PurgeList(pulse_list);
}


void pp_evaluation::FilterTooClose()
{
	for(int i = 1; i < pulse_list_len; i++)
	{
		bool con1 = ((pulse_list[i].max_pos - pulse_list[i-1].max_pos) < 0.2);
		
		if(con1)
		{
			pulse_list[i].del = true;
		}
	}
	
	pulse_list = PurgeList(pulse_list);
	pulse_list_len = pulse_list.size();
}




void pp_evaluation::FindPos()
{
	for(int i = 0; i < pulse_list_len; i++)
	{
		double p_sum = 0.0;
		double I_sum = 0.0;
		
		for(long j = pulse_list[i].left_pos; j < pulse_list[i].right_pos; j++)
		{
			p_sum += (*I_ptr)[j] * (*t_ptr)[j];
			I_sum += (*I_ptr)[j];
		}
		
		pulse_list[i].pos = p_sum / I_sum;
	}
}

void pp_evaluation::FindProps()
{
	for(int i = 0; i < pulse_list_len; i++)
	{
		double p_sum = 0.0;
		double I_sum = 0.0;
		
		for(long j = pulse_list[i].left_pos; j < pulse_list[i].right_pos; j++)
		{
			p_sum += (*I_ptr)[j] * (*t_ptr)[j];
			I_sum += (*I_ptr)[j];
		}
		
		pulse_list[i].pos = p_sum / I_sum;
		pulse_list[i].area = I_sum * dt;
		pulse_list[i].width = (*t_ptr)[pulse_list[i].right_pos-1] - (*t_ptr)[pulse_list[i].left_pos];
	}
}


void pp_evaluation::FindProps(string opt)
{
	for(int i = 0; i < pulse_list_len; i++)
	{
		double p_sum = 0.0;
		double I_sum = 0.0;
		
		for(long j = pulse_list[i].left_pos; j < pulse_list[i].right_pos; j++)
		{
			p_sum += (*I_ptr)[j] * (*t_ptr)[j];
			I_sum += (*I_ptr)[j];
		}
		
		pulse_list[i].pos = p_sum / I_sum;
		pulse_list[i].area = I_sum * dt;
		pulse_list[i].width = (*t_ptr)[pulse_list[i].right_pos-1] - (*t_ptr)[pulse_list[i].left_pos];
		
		if(opt=="fwhm")
		{
			long u = (long)(1.0/(100.0 * dt));
			
			double fwhm_left = 0.0;
			double fwhm_right = 0.0;
			
			
			for(long j = pulse_list[i].max_pos; j < pulse_list[i].right_pos - u; j++)
			{
				double mwa = 0.0;
				
				for(long k = j-u; j < j+u; k++)
				{
					mwa += (*I_ptr)[k];
				}
				
				mwa /= 2.0*u;
				
				if(mwa < 0.5 * pulse_list[i].max_val)
				{
					fwhm_right = (*t_ptr)[j];
				}
			}
			
			for(long j = pulse_list[i].max_pos; j > pulse_list[i].left_pos + u; j--)
			{				
				double mwa = 0.0;
				
				for(long k = j-u; j < j+u; k++)
				{
					mwa += (*I_ptr)[k];
				}
				
				mwa /= 2.0*u;
				
				if(mwa < 0.5 * pulse_list[i].max_val)
				{
					fwhm_left = (*t_ptr)[j];
				}
			}
			
			if(fwhm_right != 0 && fwhm_left != 0)
			{
				pulse_list[i].fwhm = fwhm_right - fwhm_left;
			}
			
		}
	}
}








vector<pulse> pp_evaluation::PurgeList(vector<pulse> init_pulse_list)
{
	vector<pulse> new_pulse_list;
	new_pulse_list.resize(0);
	
	for(unsigned i = 0; i < init_pulse_list.size(); i++)
	{
		if(!init_pulse_list[i].del)
		{
			new_pulse_list.push_back(init_pulse_list[i]);
		}
	}
	
	return new_pulse_list;
}







vector<double> pp_evaluation::GetPulseDist()
{
	vector<double> out(pulse_list_len-1, 0.0);
	
	for(int i = 0; i < pulse_list_len-1; i++)
	{
		out[i] = pulse_list[i+1].pos - pulse_list[i].pos;
	}
	
	return out;
}










/*
vector<pulse> timeseries::pulse_analysis()
{
	vector<pulse> pulse_list;
	
	auto max_it = max_element(I.begin(), I.end());
	I_max = I[distance(I.begin(), max_it)];

	double thres = 0.75 * AP.IP.D.par_dbl * I_max;
	
//	double thres_left = 0.75 * AP.IP.D.par_dbl * I_max;
//	double thres_right = 0.75 * AP.IP.D.par_dbl * I_max;
	
	bool pulse_detected = false;
	long unit = ceil(AP.LP.T.par_dbl/(1000.0*AP.IP.dt.par_dbl));
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
*/



/*

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


*/





/*

vector<double> timeseries::get_pulse_dist(vector<pulse> pulse_list)
{
	vector<double> pulse_dist;
	
	for(unsigned i = 1; i < pulse_list.size(); i++)
	{
		pulse_dist.push_back(pulse_list[i].pos - pulse_list[i-1].pos);
	}
	
	return pulse_dist;
}

*/

