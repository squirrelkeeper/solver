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
#include "sweep.hpp"



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
	
	TS = IN.integrate();
	
	TS.cut_series("last", 10);
	TS.reverse_series();
	TS.cc_series();
	TS.reset_time();
	
	
//	vector<pulse> pulse_list = TS.pulse_analysis();
	
//	TS.cout_pulse_data(pulse_list);
	
//	vector<double> pulse_dist = TS.get_pulse_dist(pulse_list);
	
	
	
	TS.write_file("test_ts");

	
	time_total.stop();
	time_total.print_elaps();
	
	
	return 0;
}

