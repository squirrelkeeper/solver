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
			double dt = AP.IP.dt.par_dbl;
			double tau1 = AP.LP.T.par_dbl;
			double tau2 = AP.FP.tau.par_dbl;
			double tauL = AP.larger_delay();
			
			
			allpar_set hom_AP = AP;
			hom_AP.IP.int_time.par_dbl = 1000.0 + 200.0 * hom_AP.LP.T.par_dbl;
			hom_AP.IP.out_time.par_dbl = 100.0;
			
			cout << "check1" << endl;
		
			vector<double> hom_const_IC = {0.4, 0.0, 4.0, 1.0, 1.0};
			initial_con hom_IC("const", hom_const_IC, AP);
	
			
			integrator hom_IN(hom_AP);
			hom_IN.initialize(hom_IC);
			
			tuple<
				timeseries,
				timeseries,
				timeseries,
				ts_evaluation
				> hom_TS_NM1_NM2_EV = hom_IN.integrate_neutral_modes_analysis("full");

			cout << "check2" << endl;
			timeseries hom_TS = get<0>(hom_TS_NM1_NM2_EV);
			timeseries NM1 = get<1>(hom_TS_NM1_NM2_EV);
			timeseries NM2 = get<2>(hom_TS_NM1_NM2_EV);

			hom_TS.write_file("test_sc_hom");
			
			double period = get<3>(hom_TS_NM1_NM2_EV).period;
			long period_int = (long)(period / dt);

			cout << period << endl;
			cout << get<3>(hom_TS_NM1_NM2_EV).state << endl;
			
			cout << "check3" << endl;			
			allpar_set adj_AP = AP;
			
			adj_AP.IP.int_time.par_dbl = 100.0;
			adj_AP.IP.out_time.par_dbl = 50.0;
			
			cout << "check4" << endl;			
			vector<double> adj1_const_IC = {2.0, 0.0, 1.0, 1.0, 1.0};
			initial_con adj1_IC("const", adj1_const_IC, adj_AP);
			
			integrator adj1_IN(adj_AP);
			adj1_IN.initialize(adj1_IC);
			
			timeseries adj1_TS = adj1_IN.integrate_adj(hom_TS);
		
			adj1_TS.write_file("test_sc_adj1");
			
			
			cout << "check5" << endl;			
			vector<double> adj2_const_IC = {0.0, 2.0, 0.0, 0.0, 0.0};
			initial_con adj2_IC("const", adj2_const_IC, AP);
			
			integrator adj2_IN(adj_AP);
			adj2_IN.initialize(adj2_IC);
			
			timeseries adj2_TS = adj2_IN.integrate_adj(hom_TS);

			adj2_TS.write_file("test_sc_adj2");
			
/*			
			hom_TS.cut_series("last", 2*tauL);
			NM1.cut_series("last", 2*tauL);
			NM2.cut_series("last", 2*tauL);
			adj1_TS.cut_series("last", 2*tauL);
			adj2_TS.cut_series("last", 2*tauL);
*/
			
			integrator bil_IN(adj_AP);

			
			bil_IN.initialize("bil");
			double b_NM1_adj1 = bil_IN.bilinear_prod(hom_TS, NM1, adj1_TS);
			
			bil_IN.initialize("bil");
			double b_NM1_adj2 = bil_IN.bilinear_prod(hom_TS, NM1, adj2_TS);
			
			bil_IN.initialize("bil");
			double b_NM2_adj1 = bil_IN.bilinear_prod(hom_TS, NM2, adj1_TS);
			
			bil_IN.initialize("bil");
			double b_NM2_adj2 = bil_IN.bilinear_prod(hom_TS, NM2, adj2_TS);
			
			
			double c_adj2 = 1/(b_NM1_adj2 - (b_NM1_adj1 * b_NM2_adj2) / b_NM2_adj1);
			double c_adj1 = 1/(b_NM1_adj1 - (b_NM1_adj2 * b_NM2_adj1) / b_NM2_adj2);
			
			cout << c_adj2 << '\t';
			cout << c_adj1 << '\t';
			cout << endl;
			
			
	//		while( < period)

//			tj = sum( c_adj1 * adj1 + c_adj2 * adj2 )
			
			break;
			(*ptr_scan_par1).par_dbl += step;
		}
	}

}
