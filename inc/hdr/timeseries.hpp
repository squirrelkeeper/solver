#ifndef __TIMESERIES_HPP_INCLUDED__
#define __TIMESERIES_HPP_INCLUDED__

#include "../hdr/parameter.hpp"
#include "../hdr/variable.hpp"

class pulse
{
public:
	long thres_pos;
	long right_pos;
	long left_pos;
	double baseline;
	double max;
	double pos;

	
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
	
	allpar_set *AP;
	double seed;

	timeseries();
	timeseries(allpar_set*);
	timeseries(double, allpar_set*);
	
	std::vector<pulse> pulse_analysis();
	
	void cout_pulse_data(std::vector<pulse>);

	
	void write_file(std::string);
};

#endif
