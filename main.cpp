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




using namespace std;



int main(int argc, char* argv[])
{
	timer time_total;

	allpar_set AP("JAU15", "NOFB", "quick");
	AP.check_cmd_line(argc, argv);

	icpar_set IC("std");
	IC.check_cmd_line(argc, argv);
	

	mode_cmd mode(argc, argv);

	string file_name = "test.dat";
	string first_line = "#\n";
	
	par *par1_ptr;
	par *par2_ptr;
	int pts = 1;
	double incr = 0.0;
	
	
	
	if(mode.mode_str == "TS_pp")
	{
		file_name = "num_D"
		+ to_string(AP.IP.D.par_dbl)
		+ ".ts.dat";
	
		first_line = "#start\tend\tpos\tdt\n";
	}
	
	
	if(mode.mode_str == "long")
	{
		file_name = "num_D"
		+ to_string(AP.IP.D.par_dbl)
		+ ".ts.dat";
	
		first_line = "#rea\ttisi\n";
	}
	
	
	
	if(mode.mode_str == "lscan")
	{
		par1_ptr = AP.get_par_ptr(mode.par1_str);
		par2_ptr = AP.get_par_ptr(mode.par2_str);
	
	
		(*par1_ptr).par_dbl = mode.par1_start;
		
		pts = mode.par1_steps;
		incr = (mode.par1_stop-mode.par1_start)/mode.par1_steps;
		
		file_name = "num_pp_"
		+ (*par2_ptr).par_str
		+ "_"
		+ to_string((*par2_ptr).par_dbl)
		+ "_"
		+(*par1_ptr).par_str
		+ "_"
		+ to_string(mode.par1_start)
		+ "_to_"
		+ to_string(mode.par1_stop)
		+ ".ls.dat";
	
		first_line = "#" 
		+ (*par1_ptr).par_str 
		+ '\t' 
		+ (*par2_ptr).par_str 
		+ '\t' 
		+ "rea" 
		+ '\t' 
		+ "pnum" 
		+ '\t' 
		+ "pp" 
		+ '\n';
	}


	
	
	
	
	ofstream out;
	out.open(file_name);
	out << first_line;
	
	

	
	
		
	if(mode.mode_str == "TS_pp")
	{
		if(AP.FP.K.par_dbl == 0)
		{
			IC.j_ic.par_dbl = 0.0;
		}
			
			
		vector<double> hom_const_IC = {
			IC.er_ic.par_dbl,
			IC.ei_ic.par_dbl,
			IC.g_ic.par_dbl, 
			IC.q_ic.par_dbl, 
			IC.j_ic.par_dbl
		};
			
		initial_con hom_IC("const", hom_const_IC, AP);
		
		integrator IN(AP);
		IN.initialize(hom_IC);
	
		tuple<timeseries, pp_evaluation> hom_TS_PP = IN.integrate_noise_analysis("mwa");
		

		for(long i2 = 0; i2 < get<1>(hom_TS_PP).pulse_list_len; i2++)
		{
			out << get<1>(hom_TS_PP).pulse_list[i2].left_pos  << '\t';
			out << get<1>(hom_TS_PP).pulse_list[i2].right_pos << '\t';
			out << setprecision(15);
			out << get<1>(hom_TS_PP).pulse_list[i2].pos << '\t';
			if(i2 > 0)
			{
				out << get<1>(hom_TS_PP).pulse_list[i2].pos - get<1>(hom_TS_PP).pulse_list[i2-1].pos;
			}
			else
			{
				out << 0.0;
			}
			
			
			out << endl;
		}
		
		get<0>(hom_TS_PP).write_file("TS_"+file_name);
		
	}
	
	
	
	
	
	if(mode.mode_str == "lscan")
	{
		for(int i1 = 0; i1 < pts; i1++)
		{
			int realisations = (int)(AP.IP.rea.par_dbl);
	
			if(AP.FP.K.par_dbl == 0)
			{
				IC.j_ic.par_dbl = 0.0;
			}
			
			
			vector<double> hom_const_IC = {
				IC.er_ic.par_dbl,
				IC.ei_ic.par_dbl,
				IC.g_ic.par_dbl, 
				IC.q_ic.par_dbl, 
				IC.j_ic.par_dbl
			};
			
			initial_con hom_IC("const", hom_const_IC, AP);
				
			for(int i2 = 0; i2 < realisations; i2++)
			{
				integrator hom_IN(AP);
				hom_IN.initialize(hom_IC);
				
				vector<double> PP = hom_IN.integrate_noise_conc_analysis();
				
				for(unsigned i3 = 1; i3 < PP.size(); i3++)
				{
					
					out << (*par1_ptr).par_dbl << '\t';
					out << (*par2_ptr).par_dbl << '\t';
					out << i2 << '\t';
					out << setprecision(15);
					out << PP[i3];
					out << endl;
				}
				
				
				
			}

			(*par1_ptr).par_dbl += incr;
			
		}
	
	}

	if(mode.mode_str == "long")
	{
		int realisations = (int)(AP.IP.rea.par_dbl);

		if(AP.FP.K.par_dbl == 0)
		{
			IC.j_ic.par_dbl = 0.0;
		}
			
			
		vector<double> hom_const_IC = {
			IC.er_ic.par_dbl,
			IC.ei_ic.par_dbl,
			IC.g_ic.par_dbl, 
			IC.q_ic.par_dbl, 
			IC.j_ic.par_dbl
		};
			
		initial_con hom_IC("const", hom_const_IC, AP);

		for(int i1 = 0; i1 < realisations; i1++)
		{
			integrator hom_IN(AP);
			hom_IN.initialize(hom_IC);
			
			vector<double> PP = hom_IN.integrate_noise_conc_analysis();
			
			for(unsigned i2 = 0; i2 < PP.size(); i2++)
			{
				
				out << i1 << '\t';
				out << setprecision(15);
				out << PP[i2];
				out << endl;


				cout << i1 << '\t';
				cout << setprecision(15);
				cout << PP[i2];
				cout << endl;

				
			}
			
			
			
		}

	}


	
	out.close();
	
	time_total.stop();
	time_total.print_elaps();
	
	
	return 0;
}





/*
 * 
 * 	

	
	
	
	
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
		

		for(long i = 0; i < get<1>(hom_TS_PP).pulse_list_len; i++)
		{
			
			cout << get<1>(hom_TS_PP).pulse_list[i].max_pos << '\t';
			cout << get<1>(hom_TS_PP).pulse_list[i].max_val << '\t';
			cout << endl;
		}




		for(long i = 1; i < get<1>(hom_TS_PP).pulse_list_len; i++)
		{
			
			cout << get<1>(hom_TS_PP).pulse_list[i].max_pos-get<1>(hom_TS_PP).pulse_list[i-1].max_pos << '\t';
			cout << endl;



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
	

	integrator temp_IN(AP);
	temp_IN = hom_IN;
	double period = temp_IN.get_period();
	
	cout << period << endl;
	
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

	
	
	
	double b11 = adj1_IN.bilinear_one_step(hom_TS, neutral_modes[0], adj1_TS);
	double b12 = adj1_IN.bilinear_one_step(hom_TS, neutral_modes[1], adj1_TS);
	double b21 = adj2_IN.bilinear_one_step(hom_TS, neutral_modes[0], adj1_TS);
	double b22 = adj2_IN.bilinear_one_step(hom_TS, neutral_modes[1], adj1_TS);
	
	cout << b11 << endl;
	cout << b12 << endl;
	cout << b21 << endl;
	cout << b22 << endl;

	
	
	adj2_TS.write_file("test_adj1");
	adj1_TS.write_file("test_adj2");
	neutral_modes[0].write_file("test_neut1");
	neutral_modes[1].write_file("test_neut2");
	
	
//	cout << get<1>(hom_TS_EV).UniqMaxVal.size() << endl;
//	cout << get<1>(hom_TS_EV).UniqMinVal.size() << endl;
	
//	get<0>(hom_TS_EV).write_file("test_hom");
	
*/
