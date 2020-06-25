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
	
	double baseline;
	double max;
	double pos;
	double width;
	
	bool del;
	
	pulse();
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

#endif
