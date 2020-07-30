#ifndef __SCAN_HPP_INCLUDED__
#define __SCAN_HPP_INCLUDED__

#include "../hdr/parameter.hpp"
#include "../hdr/variable.hpp"
#include "../hdr/timeseries.hpp"
#include "../hdr/initial_con.hpp"
#include "../hdr/integrate.hpp"

class scan
{
public:
	par *ptr_scan_par1;
	par *ptr_scan_par2;
	
	allpar_set AP;
	initial_con IC;
	integrator IN;
	timeseries TS;
	
	scan(int argc, char* argv[]);
};



#endif
