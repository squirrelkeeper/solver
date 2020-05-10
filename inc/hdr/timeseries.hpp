#ifndef __TIMESERIES_HPP_INCLUDED__
#define __TIMESERIES_HPP_INCLUDED__

#include "../hdr/parameter.hpp"
#include "../hdr/variable.hpp"

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
	
	allpar_set *AP;
	double seed;

	timeseries();
	timeseries(allpar_set*);
	timeseries(double, allpar_set*);
	
	std::vector<double> pulse_positions();
	
	void write_file(std::string);
};

#endif
