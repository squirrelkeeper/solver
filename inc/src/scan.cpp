#include <algorithm>
#include <chrono>
#include <random>

#include <complex>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <string>
#include <vector>

#include "../hdr/scan.hpp"


using namespace std;


scan::scan(
	string opt,
	string ini_ptr_scan_par1,
	string ini_ptr_scan_par2,
	vector<double> par_scope,
	allpar_set init_AP
)
{

	AP = init_AP;
	
	double start = par_scope[0];
	double stop  = par_scope[1];
	double step  = par_scope[2];
	
	vector<par*> par_coll_ptr = AP.collect_ptr();

	int check = 0;
	
	for(unsigned i = 0; i < par_coll_ptr.size(); i++)
	{
		if((*par_coll_ptr[i]).par_str == ini_ptr_scan_par1)
		{
			ptr_scan_par1 = par_coll_ptr[i];
			check++;
		}
		if((*par_coll_ptr[i]).par_str == ini_ptr_scan_par2)
		{
			ptr_scan_par2 = par_coll_ptr[i];
			check++;
		}
	}
	
	if(check!=2)
	{
		cout << "err06" << endl;
	}
	
	(*ptr_scan_par1).par_dbl = start;
	
	
	if(opt == "semianalytic_tj")
	{
		while( (*ptr_scan_par1).par_dbl < stop )
		{
			
			allpar_set hom_AP = AP;
			hom_AP.IP.int_time.par_dbl = 1000.0 + AP.larger_delay();
			hom_AP.IP.out_time.par_dbl = AP.larger_delay();
		
			vector<double> hom_const_IC = {0.4, 0.0, 4.0, 1.0, 1.0};
			initial_con hom_IC("const", hom_const_IC, AP);
	
			
			integrator hom_IN(hom_AP);
			hom_IN.initialize(hom_IC);
			
			tuple<
				timeseries,
				timeseries,
				timeseries,
				ts_evaluation
				> hom_TS_NM1_NM2_EV = hom_IN.integrate_neutral_modes_analysis("period");

			timeseries hom_TS = get<0>(hom_TS_NM1_NM2_EV).reverse_series();
			timeseries hom_NM1 = get<1>(hom_TS_NM1_NM2_EV).reverse_series();
			timeseries hom_NM2 = get<2>(hom_TS_NM1_NM2_EV).reverse_series();
			double period = get<3>(hom_TS_NM1_NM2_EV).period;
			

			
			(*ptr_scan_par1).par_dbl += step;
		}
	}

}
