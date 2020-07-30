#ifndef __TIMESERIES_HPP_INCLUDED__
#define __TIMESERIES_HPP_INCLUDED__

#include "../hdr/parameter.hpp"
#include "../hdr/variable.hpp"

class pulse
{
public:
	long trig_pos;
	long right_pos;
	long left_pos;
	long max_pos;
	
	double right_val;
	double left_val;
	double max_val;
	
	
	double baseline;

	double pos;
	double width;
	double fwhm;
	double area;

	
	bool del;
	
	
	pulse();
	pulse(long);
	
	void reset();
};

class timeseries
{
public:
	long len;
	double len_dbl;
	std::vector<double> t;
	std::vector<double> I;
	std::vector<var> X;
	
	double I_max;
	double I_min;
	double I_mid;
	double I_mean;
	
	allpar_set AP;
	double seed;

	timeseries();
	timeseries(allpar_set);
	timeseries(double, allpar_set);
	
	std::vector<pulse> pulse_analysis();
	std::vector<double> get_pulse_dist(std::vector<pulse>);
	void cout_pulse_data(std::vector<pulse>);
	
	void cut_series(std::string, double);
	void reverse_series();
	void cc_series();
	void reset_time();
	
	void write_file(std::string);
};

class ts_evaluation
{
public:
// Tolerances
	double ExtrTol = 1e-16;
	double MaxMinDscr = 1e-5;
	double AverageThres = 1e-6;
	double DblCountTol = 0.05;



	allpar_set AP;
	std::vector<double>* I_ptr;
	std::vector<double>* t_ptr;
	long len;
	
	int state, uniq_max_num;
	double av_max_dist, average, period;
	double GlobalSupr, GlobalInfi;
	double t_offset, dt;
	std::vector<double> MaxPos;
	std::vector<double> MaxVal;
	std::vector<double> MinPos;
	std::vector<double> MinVal;
	std::vector<double> UniqMaxPos;
	std::vector<double> UniqMaxVal;
	std::vector<double> UniqMinPos;
	std::vector<double> UniqMinVal;


	ts_evaluation(std::vector<double>*, std::vector<double>*, allpar_set);
	ts_evaluation(timeseries*, allpar_set);
	
	void FindUniqMax();
	void FindPeriod();
	void FindState();
	
	double InterpolQuadExtrPos(double, double, double, double, double, double);
	double InterpolQuadExtrVal(double, double, double, double, double, double, double);
	
	bool is_approx(double, double, std::string, double);
};

class pp_evaluation
{
public:
// Tolerances



	allpar_set AP;
	std::vector<double>* I_ptr;
	std::vector<double>* t_ptr;
	long len;
	
	double average;
	double GlobalSupr, GlobalInfi;
	double t_offset, dt;

	std::vector<pulse> pulse_list;
	int pulse_list_len;


	pp_evaluation(std::vector<double>*, std::vector<double>*, allpar_set);
	pp_evaluation(timeseries*, allpar_set);
	
	
	void FindGlobalExtrAvrg();
	
	void DetectPulses_simple();
	void DetectPulses_MWA();
	void DetectPulses_MWC();
	
	void FilterSamePos();
	void FilterTooClose();
	void FilterTooSmallArea();
	
	void FindPos();
	void FindProps();
	void FindProps(std::string);
	
	
	std::vector<pulse> PurgeList(std::vector<pulse>);
	
	
	std::vector<double> GetPulseDist();
	
//	void FindUniqMax();
//	void FindPeriod();
//	void FindState();
};



#endif
