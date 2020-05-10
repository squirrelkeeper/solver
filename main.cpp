#include <iostream>
#include <iomanip>
#include <fstream>

#include <sstream>
#include <string>

#include <cmath>
#include <complex>
#include <vector>
#include <random>

#include <algorithm>
#include <chrono>
#include <ctime>

#include "timer.hpp"
#include "parameter.hpp"
#include "variable.hpp"

#include "initial_con.hpp"
#include "timeseries.hpp"

#include "integrate.hpp"

typedef unsigned int uint;
typedef unsigned long int ulint;

using namespace std;



int main(int argc, char* argv[])
{
	timer time_total;

	allpar_set *AP = new allpar_set("JAU15", "TAU1", "quick");
	AP->check_cmd_line(argc, argv);


	initial_con IC(AP);
	integrator IN(AP);
 	timeseries TS(AP);

	IN.initialize(IC);
/*	
	IN.integrate_noise();

	
	TS = IN.integrate_simple_TS_noise("g");
	
	
	initial_con IC2(TS);

	
	TS.write_file("template_ts_ba");

*/
	
	time_total.stop();
	time_total.print_elaps();
	
	
	return 0;
}

