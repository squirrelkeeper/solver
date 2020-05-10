#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>
#include <random>


#include "../hdr/initial_con.hpp"


using namespace std;


initial_con::initial_con(allpar_set* init_AP)
{
	double tau = init_AP->larger_delay();
	double dt = init_AP->IP.dt.par_dbl;
	long len = (long)(floor(tau / dt)+1);
	
	hist.resize(len);
	
	for(long i = 0; i < len; i++)
	{
		hist[i].ER = 0.0;
		hist[i].EI = 0.0;
		hist[i].G = 0.0;
		hist[i].Q = 0.0;
		hist[i].J = 0.0;
	}
}

initial_con::initial_con(string opt, allpar_set* init_AP)
{
	double tau = init_AP->larger_delay();
	double dt = init_AP->IP.dt.par_dbl;
	long len = (long)(floor(tau / dt)+1);
	
	hist.resize(len);
	
	if(opt == "ones")
	{
		for(long i = 0; i < len; i++)
		{
			hist[i].ER = 1.0;
			hist[i].EI = 1.0;
			hist[i].G = 1.0;
			hist[i].Q = 1.0;
			hist[i].J = 1.0;
		}
	}
	else
	{
		cout << "err003" << endl;
	}
}



initial_con::initial_con(string opt, double first, double sec, allpar_set* init_AP)
{
	double tau = init_AP->larger_delay();
	double dt = init_AP->IP.dt.par_dbl;
	long len = (long)(floor(tau / dt)+1);
	
	hist.resize(len);
	
	if(opt == "uniform")
	{

		random_device true_rnd; 
		mt19937 rnd_gen(true_rnd()); 
		uniform_real_distribution<double> rnd_distr(first,sec); 
		
		for(long i = 0; i < len; i++)
		{
			hist[i].ER = rnd_distr(rnd_gen);
			hist[i].EI = rnd_distr(rnd_gen);
			hist[i].G = rnd_distr(rnd_gen);
			hist[i].Q = rnd_distr(rnd_gen);
			hist[i].J = rnd_distr(rnd_gen);
		}

		
	}
	else if(opt == "normal")
	{
		random_device true_rnd; 
		mt19937 rnd_gen(true_rnd()); 
		normal_distribution<double> rnd_distr(first,sec); 
		
		for(long i = 0; i < len; i++)
		{
			hist[i].ER = rnd_distr(rnd_gen);
			hist[i].EI = rnd_distr(rnd_gen);
			hist[i].G = rnd_distr(rnd_gen);
			hist[i].Q = rnd_distr(rnd_gen);
			hist[i].J = rnd_distr(rnd_gen);
		}
	}
	
	else
	{
		cout << "err003" << endl;
	}
}


initial_con::initial_con(string opt, var X_IC, allpar_set* init_AP)
{
	double tau = init_AP->larger_delay();
	double dt = init_AP->IP.dt.par_dbl;
	long len = (long)(floor(tau / dt)+1);
	
	hist.resize(len);
	
	if(opt == "const")
	{
		for(long i = 0; i < len; i++)
		{
			hist[i] = X_IC;
		}
	}	
	else
	{
		cout << "err003" << endl;
	}
}

initial_con::initial_con(timeseries TS)
{
	double tau = TS.AP->larger_delay();
	double dt = TS.AP->IP.dt.par_dbl;
	long len = (long)(floor(tau / dt)+1);

	hist.resize(len);

	if(len < TS.len)
	{
		cout << "err004" << endl;
	}
	for(long i = 0; i < len; i++)
	{
		hist[i] = TS.X[i+TS.len-len];
	}
}