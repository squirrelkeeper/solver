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
#include "scan.hpp"



using namespace std;



int main(int argc, char* argv[])
{
	timer time_total;

	allpar_set AP("JAU15", "TAU1", "quick");
	AP.check_cmd_line(argc, argv);

	vector<double> hom_const_IC = {0.4, 0.0, 4.0, 1.0, 1.0};

	
	initial_con hom_IC("const", hom_const_IC, AP);
//	initial_con adj1_IC("const", adj1_const_IC, AP);
//	initial_con adj2_IC("const", adj2_const_IC, AP);

	
//	AP.IP.D.par_dbl = 0.2;

	
//	string mode = "ts";
//	string mode = "sweep";
//	string mode = "ts_noise";
	string mode = "pp_noise";
	
	
//	scan semianalytic_timing_jitter();
	
	
	if(mode == "pp_noise")
	{
		AP.IP.int_time.par_dbl = 5000;
		AP.IP.out_time.par_dbl = AP.larger_delay();
		
		integrator init_IN(AP);
		
		init_IN.initialize(hom_IC);
		
		timeseries init_TS = init_IN.integrate_noise();
		
		
		AP.IP.int_time.par_dbl = 500;
		AP.IP.out_time.par_dbl = 100;
		
		int r = (int)(AP.IP.rea.par_dbl);

		string ofile_name = "testscan_g";
		ofile_name += to_string(AP.LP.g.par_dbl);
		ofile_name += "_Jg";
		ofile_name += to_string(AP.LP.Jg.par_dbl);	
		ofile_name += ".sc.dat";

		
		
		ofstream ofile;
		ofile.open(ofile_name);
		
		ofile << "#g=" << AP.LP.g.par_dbl << endl;
		ofile << "#Jg=" << AP.LP.Jg.par_dbl << endl;
		ofile << "#r" << '\t' << "pp" << endl;
		
		
		for(int i = 0; i < r; i++)
		{
			integrator IN(AP);
			IN.initialize(init_TS);
			
			tuple<timeseries, pp_evaluation> TS_PP = IN.integrate_noise_analysis("simple");
			
//			get<0>(TS_PP).write_file("test_hom");
			
			vector<double> ppos = get<1>(TS_PP).GetPulseDist();
			
			for(unsigned j = 0; j < ppos.size(); j++)
			{
				ofile << i << '\t';
				ofile << ppos[j] << '\t';
				ofile << endl;
			}
		}
		
		ofile.close();
		
//		init_TS.write_file("test_hom");
	}
	
	
	
	
	if(mode == "ts")
	{
		integrator hom_IN(AP);
		hom_IN.initialize(hom_IC);
		tuple<timeseries, ts_evaluation> hom_TS_EV = hom_IN.integrate_analysis("full");
		
		get<0>(hom_TS_EV).write_file("test_hom");
	}

	if(mode == "ts_noise")
	{
		integrator hom_IN(AP);
		hom_IN.initialize(hom_IC);
		
 		tuple<timeseries, pp_evaluation> hom_TS_PP = hom_IN.integrate_noise_analysis("simple");

		cout << get<1>(hom_TS_PP).GlobalSupr << endl;
		cout << get<1>(hom_TS_PP).GlobalInfi << endl;
		cout << get<1>(hom_TS_PP).average << endl;
		cout << get<1>(hom_TS_PP).pulse_list_len << endl;
		
/*
		for(long i = 0; i < get<1>(hom_TS_PP).pulse_list_len; i++)
		{
			
			cout << get<1>(hom_TS_PP).pulse_list[i].max_pos << '\t';
			cout << get<1>(hom_TS_PP).pulse_list[i].max_val << '\t';
			cout << endl;
		}
*/


/*
		for(long i = 1; i < get<1>(hom_TS_PP).pulse_list_len; i++)
		{
			
			cout << get<1>(hom_TS_PP).pulse_list[i].max_pos-get<1>(hom_TS_PP).pulse_list[i-1].max_pos << '\t';
			cout << endl;
		}
*/


		for(long i = 0; i < get<1>(hom_TS_PP).pulse_list_len; i++)
		{
			
			cout << get<1>(hom_TS_PP).pulse_list[i].left_pos << '\t';
			cout << get<1>(hom_TS_PP).pulse_list[i].right_pos << '\t';
			cout << get<1>(hom_TS_PP).pulse_list[i].del << '\t';
			cout << endl;
		}


		get<0>(hom_TS_PP).write_file("test_hom");
	}
	
	
	
	
	if(mode == "sweep")
	{
		int pts = 100;
		
		double Jg_start = 0.1;
		double Jg_stop = 12.0;
		double Jg_step = (Jg_stop-Jg_start)/pts;

		
		double q0_start = 0.1;
		double q0_stop = 8.0;
		double q0_step = (q0_stop-q0_start)/pts;
		
		AP.LP.Jg.par_dbl = Jg_start;
		AP.LP.q0.par_dbl = q0_start;
		
		

		ofstream out_sweep;
		out_sweep.open("data/out_sweep.swp.dat");
		
		for(int i = 0; i < pts; i++)
		{
			
			AP.LP.Jg.par_dbl+= Jg_step;

			for(int j = 0; j < pts; j++)
			{
				
				AP.LP.q0.par_dbl+= q0_step;
				
				
				integrator hom_IN(AP);
				hom_IN.initialize(hom_IC);

				tuple<timeseries, ts_evaluation> hom_TS_EV = hom_IN.integrate_analysis("full");
				
				
				cout << AP.LP.Jg.par_dbl << '\t';
				cout << AP.LP.q0.par_dbl << '\t';
				cout << get<1>(hom_TS_EV).state << '\t';
				cout << '\n';
				
				
				
				out_sweep << AP.LP.Jg.par_dbl << '\t';
				out_sweep << AP.LP.q0.par_dbl << '\t';
				out_sweep << get<1>(hom_TS_EV).state << '\t';
				out_sweep << '\n';
				
				break;
			}
		}
		
		out_sweep.close();
	}
	

	
	

//	AP.IP.int_time.par_dbl = 100;
//	AP.IP.out_time.par_dbl = 10;
	
/*	
	integrator temp_IN(AP);
	temp_IN = hom_IN;
	double period = temp_IN.get_period();
	
	cout << period << endl;
// */	
//	integrator adj1_IN(AP);
//	timeseries adj1_TS(AP);
//	adj1_IN.initialize(adj1_IC);
//	adj1_TS = adj1_IN.integrate_adj(hom_TS);

//	cout << 1 << endl;
	
//	integrator adj2_IN(AP);
//	timeseries adj2_TS(AP);
//	adj2_IN.initialize(adj2_IC);
//	adj2_TS = adj2_IN.integrate_adj(hom_TS);
	
	
//	vector<timeseries> neutral_modes = hom_IN.integrate_get_neutral_modes();

//	neutral_modes[0].cut_series("last", 10);
//	neutral_modes[1].cut_series("last", 10);

	
	
	
/*	
	double b11 = adj1_IN.bilinear_one_step(hom_TS, neutral_modes[0], adj1_TS);
	double b12 = adj1_IN.bilinear_one_step(hom_TS, neutral_modes[1], adj1_TS);
	double b21 = adj2_IN.bilinear_one_step(hom_TS, neutral_modes[0], adj1_TS);
	double b22 = adj2_IN.bilinear_one_step(hom_TS, neutral_modes[1], adj1_TS);
	
	cout << b11 << endl;
	cout << b12 << endl;
	cout << b21 << endl;
	cout << b22 << endl;
*/
	
/*	
	adj2_TS.write_file("test_adj1");
	adj1_TS.write_file("test_adj2");
	neutral_modes[0].write_file("test_neut1");
	neutral_modes[1].write_file("test_neut2");
*/	
	
//	cout << get<1>(hom_TS_EV).UniqMaxVal.size() << endl;
//	cout << get<1>(hom_TS_EV).UniqMinVal.size() << endl;
	
//	get<0>(hom_TS_EV).write_file("test_hom");
	
	
	
	
	
	
	
	
	
	time_total.stop();
	time_total.print_elaps();
	
	
	return 0;
}

