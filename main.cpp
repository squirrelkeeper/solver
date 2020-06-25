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


	allpar_set AP("JAU15", "TAU1", "quick");
	AP.check_cmd_line(argc, argv);

	vector<double> hom_const_IC = {0.4, 0.0, 4.0, 1.0, 1.0};
	vector<double> adj1_const_IC = {2.0, 0.0, 1.0, 1.0, 1.0};
	vector<double> adj2_const_IC = {0.0, 1.0, 0.0, 0.0, 0.0};
	
	
	
	initial_con hom_IC("const", hom_const_IC, AP);
	initial_con adj1_IC("const", adj1_const_IC, AP);
	initial_con adj2_IC("const", adj2_const_IC, AP);

	
	AP.IP.int_time.par_dbl = 1000;
	AP.IP.out_time.par_dbl = 200;
	
	integrator hom_IN(AP);
	timeseries hom_TS(AP);
	hom_IN.initialize(hom_IC);
	hom_TS = hom_IN.integrate();

	AP.IP.int_time.par_dbl = 100;
	AP.IP.out_time.par_dbl = 10;
	
/*	
	integrator temp_IN(AP);
	temp_IN = hom_IN;
	double period = temp_IN.get_period();
	
	cout << period << endl;
// */	
	integrator adj1_IN(AP);
	timeseries adj1_TS(AP);
	adj1_IN.initialize(adj1_IC);
	adj1_TS = adj1_IN.integrate_adj(hom_TS);

	cout << 1 << endl;
	
	integrator adj2_IN(AP);
	timeseries adj2_TS(AP);
	adj2_IN.initialize(adj2_IC);
	adj2_TS = adj2_IN.integrate_adj(hom_TS);
	
	
	vector<timeseries> neutral_modes = hom_IN.integrate_get_neutral_modes();

	neutral_modes[0].cut_series("last", 10);
	neutral_modes[1].cut_series("last", 10);

	
	
	
	
	double b11 = adj1_IN.bilinear_one_step(hom_TS, neutral_modes[0], adj1_TS);
	double b12 = adj1_IN.bilinear_one_step(hom_TS, neutral_modes[1], adj1_TS);
	double b21 = adj2_IN.bilinear_one_step(hom_TS, neutral_modes[0], adj1_TS);
	double b22 = adj2_IN.bilinear_one_step(hom_TS, neutral_modes[1], adj1_TS);
	
	cout << b11 << endl;
	cout << b12 << endl;
	cout << b21 << endl;
	cout << b22 << endl;
	
/*	
	hom_TS.write_file("test_hom");
	adj2_TS.write_file("test_adj1");
	adj1_TS.write_file("test_adj2");
	neutral_modes[0].write_file("test_neut1");
	neutral_modes[1].write_file("test_neut2");
*/	
	
	
	
	
	
	
	
	
	
	
	time_total.stop();
	time_total.print_elaps();
	
	
	return 0;
}

