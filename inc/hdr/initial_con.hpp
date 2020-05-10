#ifndef __INITIAL_CON_HPP_INCLUDED__
#define __INITIAL_CON_HPP_INCLUDED__

#include "../hdr/parameter.hpp"
#include "../hdr/variable.hpp"
#include "../hdr/timeseries.hpp"

class initial_con
{
public:
	std::vector<var> hist;

	initial_con(allpar_set*);
	initial_con(std::string, allpar_set*);
	initial_con(std::string, double, double, allpar_set*);
	initial_con(std::string, var, allpar_set*);
	initial_con(timeseries);
};


#endif
