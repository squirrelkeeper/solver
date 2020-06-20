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

/*
	allpar_set *AP = new allpar_set("JAU15", "TAU1", "quick");
	AP->check_cmd_line(argc, argv);

	
	var X_IC;
	X_IC.ER = 0.4;
	X_IC.EI = 0.0;
	X_IC.G = 4.0;
	X_IC.Q = 1.0;
	X_IC.J = 0.4;
	initial_con hom_IC("const", X_IC, AP);
	

	var Y_IC;
	Y_IC.ER = 2.0;
	Y_IC.EI = 0.0;
	Y_IC.G = 1.0;
	Y_IC.Q = 1.0;
	Y_IC.J = 0.0;
	initial_con ret_IC("const", Y_IC, AP);

	
	var Z_IC;
	Z_IC.ER = 0.0;
	Z_IC.EI = 1.0;
	Z_IC.G = 0.0;
	Z_IC.Q = 0.0;
	Z_IC.J = 0.0;
	initial_con adj_IC("const", Z_IC, AP);
	
	
	AP->IP.int_time.par_dbl = 1000;
	AP->IP.out_time.par_dbl = 400;
	
	integrator hom_IN(AP);
	timeseries hom_TS(AP);
	hom_IN.initialize(hom_IC);
	hom_TS = hom_IN.integrate();

	
	

		
	
	
	
	AP->IP.int_time.par_dbl = 200;
	AP->IP.out_time.par_dbl = 50;

	integrator ret_IN(AP);
	timeseries ret_TS(AP);
	ret_IN.initialize(ret_IC);
	ret_TS = ret_IN.integrate_ret(hom_TS);

	integrator adj_IN(AP);
	timeseries adj_TS(AP);
	adj_IN.initialize(ret_IC);
	adj_TS = adj_IN.integrate_adj(hom_TS);

	


	adj_TS.cut_series("last", 2);


//	adj_TS.reverse_series();
//	adj_TS.cc_series();
	adj_TS.reset_time();


	/*	
	ret_TS.cut_series("last", 10);
	ret_TS.reset_time();
	
	hom_TS.cut_series("last", 10);
	hom_TS.reset_time();
*/
	

	
/*	
	
	
	double period = hom_IN.get_period();
	
	cout << period << endl;
	
	
	
	
	
	
	
	
	
	
	
	
	
//	vector<double> bil_prod = adj_IN.bilinear_prod(hom_TS, ret_TS, adj_TS);
	
//	vector<pulse> pulse_list = TS.pulse_analysis();
	
//	TS.cout_pulse_data(pulse_list);
	
//	vector<double> pulse_dist = TS.get_pulse_dist(pulse_list);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
//	hom_TS.write_file("test_hom_ts");
	ret_TS.write_file("test_ret");
	adj_TS.write_file("test_adj");
	
	time_total.stop();
	time_total.print_elaps();
*/	
	
	return 0;
}

