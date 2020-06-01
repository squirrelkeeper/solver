#ifndef __SWEEP_HPP_INCLUDED__
#define __SWEEP_HPP_INCLUDED__

#include "../hdr/parameter.hpp"
#include "../hdr/variable.hpp"
#include "../hdr/timeseries.hpp"
#include "../hdr/initial_con.hpp"
#include "../hdr/integrate.hpp"

class sweep
{
public:
	par *ptr_sweep_par1;
	par *ptr_sweep_par2;
	
	allpar_set AP;
	initial_con IC;
	integrator IN;
	timeseries TS;
};



#endif
