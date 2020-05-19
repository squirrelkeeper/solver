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


	var X_IC;
	X_IC.ER = 0.4;
	X_IC.EI = 0.0;
	X_IC.G = 4.0;
	X_IC.Q = 1.0;
	X_IC.J = 0.4;
	
	initial_con IC("const", X_IC, AP);
	integrator IN(AP);
 	timeseries TS(AP);

	IN.initialize(IC);
	
	TS = IN.integrate_noise();
	
	TS.pulse_positions();
	
	TS.write_file("template_ts_ba");

	
	time_total.stop();
	time_total.print_elaps();
	
	
	return 0;
}

