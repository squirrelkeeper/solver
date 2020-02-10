#ifndef __INTEGRATE_HPP_INCLUDED__
#define __INTEGRATE_HPP_INCLUDED__

#include "../hdr/parameter.hpp"
#include "../hdr/variable.hpp"

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
	long it;
	long dim1;
	long dim2;
	long pos0;
	long pos1;
	long pos2;
	double t;
	
	integrator(allpar_set *AP);
};

//###########################################



#endif
