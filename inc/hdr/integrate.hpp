#ifndef __INTEGRATE_HPP_INCLUDED__
#define __INTEGRATE_HPP_INCLUDED__

#include "../hdr/parameter.hpp"
#include "../hdr/variable.hpp"

class initial_condition
{
public:
	std::vector<var> X;
	double t;
	
	initial_condition();
};

//###########################################

class timeseries
{
public:
	std::vector<var> X;
	std::vector<double> I;
	std::vector<double> t;
	
	timeseries();
};

//###########################################

class integrator
{
public:
	std::vector<varC> C;
	std::vector<var> R;
	allpar_set *AP;
	long it;
	long dim1;
	long dim2;
	long pos0;
	long pos1;
	long pos2;
	double t;
	
	integrator(allpar_set *AP);
	
	timeseries integrate(std::string opt);
	void initialize(std::vector<varC> &C, std::string opt);
	void initialize(std::vector<var> &R, std::string opt);
	void initialize(std::vector<varC> &C, initial_condition ic);
	void initialize(std::vector<var> &R, initial_condition ic);
};

//###########################################



#endif
