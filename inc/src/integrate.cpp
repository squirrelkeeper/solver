#include <algorithm>
#include <chrono>
#include <random>

#include <complex>
#include <cmath>
#include <numeric>


#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <string>
#include <vector>
#include <tuple>


#include "../hdr/integrate.hpp"


using namespace std;




static void LOOKUP_EXP_INIT()
{
	for(int i=0; i < LOOKUP_SAMPLES; i++)
	{
		LOOKUP_EXP_Y0[i] = exp(i*LOOKUP_STEP);
    }

	for(int i=0; i < LOOKUP_SAMPLES - 1; i++)
	{
		LOOKUP_EXP_S[i] = LOOKUP_EXP_Y0[i+1] - LOOKUP_EXP_Y0[i];
		LOOKUP_EXP_Y0[i] = LOOKUP_EXP_Y0[i] - i * LOOKUP_EXP_S[i];
    }
}

static void LOOKUP_SIN_INIT()
{
	for(int i=0; i < LOOKUP_SAMPLES; i++)
	{
		LOOKUP_SIN_Y0[i] = sin(i*LOOKUP_STEP);
    }

	for(int i=0; i < LOOKUP_SAMPLES - 1; i++)
	{
		LOOKUP_SIN_S[i] = LOOKUP_SIN_Y0[i+1] - LOOKUP_SIN_Y0[i];
		LOOKUP_SIN_Y0[i] = LOOKUP_SIN_Y0[i] - i * LOOKUP_SIN_S[i];
    }
}

static void LOOKUP_COS_INIT()
{
	for(int i=0; i < LOOKUP_SAMPLES; i++)
	{
		LOOKUP_COS_Y0[i] = cos(i*LOOKUP_STEP);
    }

	for(int i=0; i < LOOKUP_SAMPLES - 1; i++)
	{
		LOOKUP_COS_S[i] = LOOKUP_COS_Y0[i+1] - LOOKUP_COS_Y0[i];
		LOOKUP_COS_Y0[i] = LOOKUP_COS_Y0[i] - i * LOOKUP_COS_S[i];
    }
}




//###########################################

integrator::integrator(allpar_set ap_init)
{
	AP = ap_init;
	
	lpar_dbl_set lp(AP);
	fpar_dbl_set fp(AP);
	ipar_dbl_set ip(AP);

	max_rea = (int)(ip.rea);
	
	it = (long)(ip.int_time / ip.dt);

	double tau1 = lp.T;
	double tau2 = lp.T + fp.tau;

	dim1 = (long)(floor(tau1 / ip.dt));
	dim2 = (long)(floor(tau2 / ip.dt)+1);
	
	pos0 = dim2 - 1;
	pos1 = pos0 - dim1;
	pos2 = 0;

	Time = 0.0;
	rea = 0;

	noise = (AP.IP.D.par_dbl == 0.0) ? false : true; 
	
	X.resize(dim2);
}


void integrator::initialize(initial_con IC)
{
	noise = (AP.IP.D.par_dbl == 0.0) ? false : true; 

	if(dim2 != ((long)IC.hist.size()))
	{
		cout << "err005" << endl;
	}
	
	for(long i = 0; i < dim2; i++)
	{
		X[i] = IC.hist[i];
	}	
}





timeseries integrator::integrate()
{
	timeseries TS(AP);
// 	
	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();


	lpar_dbl_set *lp = new lpar_dbl_set(AP);
	fpar_dbl_set *fp = new fpar_dbl_set(AP);
	ipar_dbl_set *ip = new ipar_dbl_set(AP);
	
	
	for(long i=0; i < it-TS.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	}
	
	for(long i=0; i < TS.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	
		TS.X[i] = Xnew;
		TS.t[i] = Time;
		TS.I[i] = Xnew.ER*Xnew.ER+Xnew.EI*Xnew.EI;
	}

	
	return TS;
}



tuple<timeseries,ts_evaluation> integrator::integrate_analysis(string opt)
{
	timeseries TS(AP);

	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();


	lpar_dbl_set *lp = new lpar_dbl_set(AP);
	fpar_dbl_set *fp = new fpar_dbl_set(AP);
	ipar_dbl_set *ip = new ipar_dbl_set(AP);
	
	double I3[3] = {0,0,0}; 
	double sum = 0;
	
	
	
	for(long i=0; i < it-TS.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	}
	
	ts_evaluation EV(&TS, AP);

	
	for(long i=0; i < TS.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	
		I3[i%3] = Xnew.ER*Xnew.ER+Xnew.EI*Xnew.EI;
		
		TS.X[i] = Xnew;
		TS.t[i] = Time;
		TS.I[i] = I3[i%3];
		
		sum += I3[i%3];
		
		if(EV.GlobalSupr < I3[i%3])
		{
			EV.GlobalSupr = I3[i%3];
		}
		
		if(EV.GlobalInfi > I3[i%3])
		{
			EV.GlobalInfi = I3[i%3];
		}
		
		if(i>1 && I3[(i-2)%3] < (1.0-EV.ExtrTol) * I3[(i-1)%3] && (1.0-EV.ExtrTol) * I3[(i-1)%3] > I3[i%3])
		{
			double tmax = InterpolQuadExtrPos(
				Time - 2*ip->dt,
				Time - ip->dt,
				Time,
				I3[(i-2)%3],
				I3[(i-1)%3],
				I3[i%3]
			);
			double Imax = InterpolQuadExtrVal(
				Time - 2*ip->dt,
				Time - ip->dt,
				Time,
				I3[(i-2)%3],
				I3[(i-1)%3],
				I3[i%3],
				tmax
			);
			
			EV.MaxPos.push_back(tmax);
			EV.MaxVal.push_back(Imax);
		}
		
	}
	
	EV.average = sum / (TS.len * 1.0);
	
	if(opt == "full")
	{
		EV.FindUniqMax();
		EV.FindPeriod();
		EV.FindState();
	}
	else if(opt == "period")
	{
		EV.FindPeriod();
	}
	
	
	tuple<timeseries, ts_evaluation> OUT = make_tuple(TS, EV);
	return OUT;
}



double integrator::InterpolQuadExtrPos(double a, double b, double c, double fa, double fb, double fc)
{
	return b + (-0.5*(-a + b)*(-a + b)*(-fb + fc) + 0.5*(-b + c)*(-b + c)*(fa - fb))/((-a + b)*(-fb + fc) + (-b + c)*(fa - fb));
}


double integrator::InterpolQuadExtrVal(double a, double b, double c, double fa, double fb, double fc, double tmax)
{
	return (fa*(b - c)*(b - tmax)*(c - tmax) - fb*(a - c)*(a - tmax)*(c - tmax) + fc*(a - b)*(a - tmax)*(b - tmax))/((a - b)*(a - c)*(b - c));
}



tuple<
	timeseries,
	timeseries,
	timeseries,
	ts_evaluation
> integrator::integrate_neutral_modes_analysis(string opt)
{
	timeseries TS(AP);
	timeseries NM1(AP);
	timeseries NM2(AP);


	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();


	lpar_dbl_set *lp = new lpar_dbl_set(AP);
	fpar_dbl_set *fp = new fpar_dbl_set(AP);
	ipar_dbl_set *ip = new ipar_dbl_set(AP);
	
	double I3[3] = {0,0,0}; 
	double sum = 0;
	
	
	
	for(long i=0; i < it-TS.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	}
	
	ts_evaluation EV(&TS, AP);

	
	for(long i=0; i < TS.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	
		I3[i%3] = Xnew.ER*Xnew.ER+Xnew.EI*Xnew.EI;
		
		TS.X[i] = Xnew;
		TS.t[i] = Time;
		TS.I[i] = I3[i%3];
		
		NM1.X[i] = dX;
		NM1.t[i] = Time;
		NM1.I[i] = dX.ER * dX.ER + dX.EI * dX.EI;
		
		NM2.X[i].ER = -Xnew.EI;
		NM2.X[i].EI =  Xnew.ER;
		NM2.X[i].G =  0.0;
		NM2.X[i].Q =  0.0;
		NM2.X[i].J =  0.0;
		
		NM2.t[i] = Time;
		NM2.I[i] = Xnew.ER*Xnew.ER+Xnew.EI*Xnew.EI;
		
		
		
		sum += I3[i%3];
		
		if(EV.GlobalSupr < I3[i%3])
		{
			EV.GlobalSupr = I3[i%3];
		}
		
		if(EV.GlobalInfi > I3[i%3])
		{
			EV.GlobalInfi = I3[i%3];
		}
		
		if(i>1 && I3[(i-2)%3] < (1.0-EV.ExtrTol) * I3[(i-1)%3] && (1.0-EV.ExtrTol) * I3[(i-1)%3] > I3[i%3])
		{
			double tmax = InterpolQuadExtrPos(
				Time - 2*ip->dt,
				Time - ip->dt,
				Time,
				I3[(i-2)%3],
				I3[(i-1)%3],
				I3[i%3]
			);
			double Imax = InterpolQuadExtrVal(
				Time - 2*ip->dt,
				Time - ip->dt,
				Time,
				I3[(i-2)%3],
				I3[(i-1)%3],
				I3[i%3],
				tmax
			);
			
			EV.MaxPos.push_back(tmax);
			EV.MaxVal.push_back(Imax);
		}
		
	}
	
	EV.average = sum / (TS.len * 1.0);
	
	if(opt == "full")
	{
		EV.FindUniqMax();
		EV.FindPeriod();
		EV.FindState();
	}
	else if(opt == "period")
	{
		EV.FindPeriod();
	}
	
	
	tuple<
		timeseries,
		timeseries,
		timeseries,
		ts_evaluation
		> OUT = make_tuple(TS, NM1, NM2, EV);
		
	return OUT;
}








timeseries integrator::integrate_noise()
{
	timeseries TS(AP);

	random_device true_rnd; 
	TS.seed = true_rnd();
	mt19937 rnd_gen(TS.seed); 
	normal_distribution<double> rnd_distr(0,1); 
	
	
	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();


	lpar_dbl_set *lp = new lpar_dbl_set(AP);
	fpar_dbl_set *fp = new fpar_dbl_set(AP);
	ipar_dbl_set *ip = new ipar_dbl_set(AP);
	

	for(long i=0; i < it-TS.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		Xnew.ER+= ip->D * ip->sqrtdt * rnd_distr(rnd_gen);
		Xnew.EI+= ip->D * ip->sqrtdt * rnd_distr(rnd_gen);
		
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	}
	
	for(long i=0; i < TS.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		Xnew.ER+= ip->D * ip->sqrtdt * rnd_distr(rnd_gen);
		Xnew.EI+= ip->D * ip->sqrtdt * rnd_distr(rnd_gen);
		
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	
		TS.X[i] = Xnew;
		TS.t[i] = Time;
		TS.I[i] = Xnew.ER*Xnew.ER+Xnew.EI*Xnew.EI;
	}

	return TS;
}




tuple<timeseries, pp_evaluation> integrator::integrate_noise_analysis(string opt)
{
	timeseries TS(AP);
	
	random_device true_rnd; 
	TS.seed = true_rnd();
	mt19937 rnd_gen(TS.seed); 
	normal_distribution<double> rnd_distr(0,1); 
	
	
	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();


	lpar_dbl_set *lp = new lpar_dbl_set(AP);
	fpar_dbl_set *fp = new fpar_dbl_set(AP);
	ipar_dbl_set *ip = new ipar_dbl_set(AP);
	
	double sum = 0.0;
	double I;
	

	for(long i=0; i < it-TS.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		Xnew.ER+= ip->D * ip->sqrtdt * rnd_distr(rnd_gen);
		Xnew.EI+= ip->D * ip->sqrtdt * rnd_distr(rnd_gen);
		
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	}
	
	pp_evaluation PP(&TS, AP);

	cout << lp->T << endl;
	cout << fp->tau << endl;
	cout << dim1 << endl;
	cout << dim2 << endl;

	
	for(long i=0; i < TS.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		Xnew.ER+= ip->D * ip->sqrtdt * rnd_distr(rnd_gen);
		Xnew.EI+= ip->D * ip->sqrtdt * rnd_distr(rnd_gen);
		
		
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		
		Time += ip->dt;
		I = Xnew.ER*Xnew.ER+Xnew.EI*Xnew.EI;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	
		TS.X[i] = Xnew;
		TS.t[i] = Time;
		TS.I[i] = I;
		
		sum += I;
		
		if(I > PP.GlobalSupr)
		{
			PP.GlobalSupr = I;
		}
		
		if(I < PP.GlobalInfi)
		{
			PP.GlobalInfi = I;
		}
	}
	
	PP.average = sum/(TS.len * 1.0);
	
		
	if(opt == "simple")
	{
		PP.DetectPulses_simple();
		
		PP.FilterSamePos();
		PP.FilterTooClose();
		
		PP.FindPos();
	}
	else if(opt == "mwa")
	{
		PP.DetectPulses_MWA();
		
		PP.FilterSamePos();
		PP.FilterTooClose();
		
		PP.FindPos();
	}
	
	
	tuple<timeseries, pp_evaluation> OUT = make_tuple(TS, PP);
	return OUT;
}



vector<timeseries> integrator::integrate_get_neutral_modes()
{
	timeseries mode1(AP);
	timeseries mode2(AP);

	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();


	lpar_dbl_set *lp = new lpar_dbl_set(AP);
	fpar_dbl_set *fp = new fpar_dbl_set(AP);
	ipar_dbl_set *ip = new ipar_dbl_set(AP);
	
	
	for(long i=0; i < it-mode1.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	}
	
	for(long i=0; i < mode1.len; i++)
	{
		
		dX = derive_real(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
	
		
		mode1.X[i] = dX;
		mode1.t[i] = Time;
		mode1.I[i] = dX.ER*dX.ER+dX.EI*dX.EI;
		
		
		mode2.X[i].ER = -Xnew.EI;
		mode2.X[i].EI =  Xnew.ER;
		mode2.X[i].G =  0;
		mode2.X[i].Q =  0;
		mode2.X[i].J =  0;
		
		mode2.t[i] = Time;
		mode2.I[i] = Xnew.ER*Xnew.ER+Xnew.EI*Xnew.EI;
	}
	
	vector<timeseries> modes = {mode1,mode2};

	
 	return modes;
}











timeseries integrator::integrate_ret(timeseries hom)
{
	timeseries TS(AP);

	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();
	
	lpar_dbl_set *lp = new lpar_dbl_set(AP);
	fpar_dbl_set *fp = new fpar_dbl_set(AP);
	ipar_dbl_set *ip = new ipar_dbl_set(AP);

	double tau1 = AP.smaller_delay();
	double tau2 = AP.larger_delay();
	
	long hpos0 = (long)((double)(hom.X.size())*0.5);
	long hpos1 = hpos0 - (long)(floor(tau1 / ip->dt));
	long hpos2 = hpos0 - (long)(floor(tau2 / ip->dt)+1);
	

	for(long i=0; i < it-TS.len; i++)
	{
		
		dX = derive_ret(
			hom.X[hpos0], hom.X[hpos1], hom.X[hpos2],
			X[pos0], X[pos1], X[pos2],
			lp, fp
		);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
		
		hpos0++;
		hpos1++;
		hpos2++;
	}
	
	for(long i=0; i < TS.len; i++)
	{
		
		dX = derive_ret(
			hom.X[hpos0], hom.X[hpos1], hom.X[hpos2],
			X[pos0], X[pos1], X[pos2],
			lp, fp
		);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
		
		hpos0++;
		hpos1++;
		hpos2++;
	
		TS.X[i] = Xnew;
		TS.t[i] = Time;
		TS.I[i] = Xnew.ER*Xnew.ER+Xnew.EI*Xnew.EI;
	}

	
	return TS;
}


timeseries integrator::integrate_adj(timeseries hom)
{
	timeseries TS(AP);

	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();
	
	lpar_dbl_set *lp = new lpar_dbl_set(AP);
	fpar_dbl_set *fp = new fpar_dbl_set(AP);
	ipar_dbl_set *ip = new ipar_dbl_set(AP);

	long hpos0 = hom.X.size()-1;
	
	if(hpos0+1 < (long)(AP.IP.int_time.par_dbl*AP.IP.dt.par_dbl))
	{
		cout << "err07" << endl;
	}
	
	
	for(long i=0; i < it-TS.len; i++)
	{
		
		dX = derive_adj(
			hom.X[hpos0], X[pos0], X[pos1], X[pos2],
			lp, fp
		);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
		
		hpos0--;
	}
	
	for(long i=0; i < TS.len; i++)
	{
		
		dX = derive_adj(
			hom.X[hpos0], X[pos0], X[pos1], X[pos2],
			lp, fp
		);
		
		Xnew.ER = X[pos0].ER + ip->dt * dX.ER;
		Xnew.EI = X[pos0].EI + ip->dt * dX.EI;
		Xnew.G  = X[pos0].G  + ip->dt * dX.G;
		Xnew.Q  = X[pos0].Q  + ip->dt * dX.Q;
		Xnew.J  = X[pos0].J  + ip->dt * dX.J;
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += ip->dt;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;
		
		hpos0--;
	
		TS.X[i] = Xnew;
		TS.t[i] = Time;
		TS.I[i] = Xnew.ER*Xnew.ER+Xnew.EI*Xnew.EI;
	}

	
	return TS;
}


	




double integrator::bilinear_one_step(timeseries hom, timeseries ret, timeseries adj)
{
	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();
	
	lpar_dbl_set *lp = new lpar_dbl_set(AP);
	fpar_dbl_set *fp = new fpar_dbl_set(AP);
//	ipar_dbl_set *ip = new ipar_dbl_set(AP);
	
	pos0 = (long)(hom.X.size() * 0.5);
	
	double bil;
	
	
	bil = bilinear_step(hom.X, ret.X, adj.X, lp, fp);

	
	return bil;
}





double integrator::bilinear_prod(timeseries hom, timeseries ret, timeseries adj)
{
	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();
	
	lpar_dbl_set *lp = new lpar_dbl_set(AP);
	fpar_dbl_set *fp = new fpar_dbl_set(AP);
//	ipar_dbl_set *ip = new ipar_dbl_set(AP);

//	pos0?
//	dim1 / dim2?
	
	
	double bil = bilinear_step(hom.X, ret.X, adj.X, lp, fp);

	
	return bil;
}






/*REPLACE START*/
var integrator::derive_real(var &X, var &XT, var &Xtau, lpar_dbl_set *l, fpar_dbl_set *f)
{
	var d;

	d.ER = -X.ER*l->g + XT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + XT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw);

	d.EI = -X.EI*l->g + XT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - XT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw);

	d.G = -X.EI*X.EI*expf(X.G - X.Q) + X.EI*X.EI*expf(-X.Q) - X.ER*X.ER*expf(X.G - X.Q) + X.ER*X.ER*expf(-X.Q) - X.G*l->gg + l->Jg;

	d.Q = -X.EI*X.EI*l->rs + X.EI*X.EI*l->rs*expf(-X.Q) - X.ER*X.ER*l->rs + X.ER*X.ER*l->rs*expf(-X.Q) - X.J*X.Q + X.J*l->q0 - X.Q*l->gq + l->gq*l->q0;

	d.J = -X.J*f->wLP + Xtau.EI*Xtau.EI*f->K*f->wLP + Xtau.ER*Xtau.ER*f->K*f->wLP;


	return d;
}


var integrator::derive_ret(var &X, var &XT, var &Xtau, var &Y, var &YT, var &Ytau, lpar_dbl_set *l, fpar_dbl_set *f)
{
	//X is the homogenous solution, Y the pertubation
	var d;

	d.ER = 0.5*XT.EI*YT.G*l->ag*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.EI*YT.G*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.EI*YT.Q*l->aq*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.EI*YT.Q*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.ER*YT.G*l->ag*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.ER*YT.G*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.ER*YT.Q*l->aq*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.ER*YT.Q*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - Y.ER*l->g + YT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + YT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw);

	d.EI = -0.5*XT.EI*YT.G*l->ag*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.EI*YT.G*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.EI*YT.Q*l->aq*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.EI*YT.Q*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.ER*YT.G*l->ag*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.ER*YT.G*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.ER*YT.Q*l->aq*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.ER*YT.Q*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - Y.EI*l->g + YT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - YT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw);

	d.G = -X.EI*X.EI*Y.G*expf(X.G - X.Q) + X.EI*X.EI*Y.Q*expf(X.G - X.Q) - X.EI*X.EI*Y.Q*expf(-X.Q) - 2*X.EI*Y.EI*expf(X.G - X.Q) + 2*X.EI*Y.EI*expf(-X.Q) - X.ER*X.ER*Y.G*expf(X.G - X.Q) + X.ER*X.ER*Y.Q*expf(X.G - X.Q) - X.ER*X.ER*Y.Q*expf(-X.Q) - 2*X.ER*Y.ER*expf(X.G - X.Q) + 2*X.ER*Y.ER*expf(-X.Q) - Y.G*l->gg;

	d.Q = -X.EI*X.EI*Y.Q*l->rs*expf(-X.Q) - 2*X.EI*Y.EI*l->rs + 2*X.EI*Y.EI*l->rs*expf(-X.Q) - X.ER*X.ER*Y.Q*l->rs*expf(-X.Q) - 2*X.ER*Y.ER*l->rs + 2*X.ER*Y.ER*l->rs*expf(-X.Q) - X.J*Y.Q - X.Q*Y.J + Y.J*l->q0 - Y.Q*l->gq;

	d.J = 2*Xtau.EI*Ytau.EI*f->K*f->wLP + 2*Xtau.ER*Ytau.ER*f->K*f->wLP - Y.J*f->wLP;


	return d;
}


var integrator::derive_adj(var &X, var &Z, var &ZT, var &Ztau, lpar_dbl_set *l, fpar_dbl_set *f)
{
	//X is the homogenous solution, Z the adjoint pertubation
	var d;

	d.ER = -2*X.ER*Z.G*expf(X.G - X.Q) + 2*X.ER*Z.G*expf(-X.Q) - 2*X.ER*Z.Q*l->rs + 2*X.ER*Z.Q*l->rs*expf(-X.Q) + 2*X.ER*Ztau.J*f->K*f->wLP - Z.ER*l->g - ZT.EI*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*sinf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) + ZT.ER*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*cosf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw);

	d.EI = -2*X.EI*Z.G*expf(X.G - X.Q) + 2*X.EI*Z.G*expf(-X.Q) - 2*X.EI*Z.Q*l->rs + 2*X.EI*Z.Q*l->rs*expf(-X.Q) + 2*X.EI*Ztau.J*f->K*f->wLP - Z.EI*l->g + ZT.EI*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*cosf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) + ZT.ER*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*sinf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw);

	d.G = -X.EI*X.EI*Z.G*expf(X.G - X.Q) - X.ER*X.ER*Z.G*expf(X.G - X.Q) - 0.5*X.EI*ZT.EI*l->ag*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*sinf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) + 0.5*X.EI*ZT.EI*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*cosf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) + 0.5*X.EI*ZT.ER*l->ag*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*cosf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) + 0.5*X.EI*ZT.ER*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*sinf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) - 0.5*X.ER*ZT.EI*l->ag*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*cosf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) - 0.5*X.ER*ZT.EI*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*sinf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) - 0.5*X.ER*ZT.ER*l->ag*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*sinf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) + 0.5*X.ER*ZT.ER*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*cosf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) - Z.G*l->gg;

	d.Q = X.EI*X.EI*Z.G*expf(X.G - X.Q) - X.EI*X.EI*Z.G*expf(-X.Q) - X.EI*X.EI*Z.Q*l->rs*expf(-X.Q) + X.ER*X.ER*Z.G*expf(X.G - X.Q) - X.ER*X.ER*Z.G*expf(-X.Q) - X.ER*X.ER*Z.Q*l->rs*expf(-X.Q) - X.J*Z.Q + 0.5*X.EI*ZT.EI*l->aq*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*sinf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) - 0.5*X.EI*ZT.EI*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*cosf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) - 0.5*X.EI*ZT.ER*l->aq*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*cosf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) - 0.5*X.EI*ZT.ER*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*sinf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) + 0.5*X.ER*ZT.EI*l->aq*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*cosf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) + 0.5*X.ER*ZT.EI*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*sinf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) + 0.5*X.ER*ZT.ER*l->aq*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*sinf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) - 0.5*X.ER*ZT.ER*l->g*l->sqrtkap*expf(0.5*X.G - 0.5*X.Q)*cosf(0.5*X.G*l->ag - 0.5*X.Q*l->aq + l->T*l->dw) - Z.Q*l->gq;

	d.J = -X.Q*Z.Q - Z.J*f->wLP + Z.Q*l->q0;


	return d;
}


double integrator::bilinear_step(vector<var> &X, vector<var> &Y, vector<var> &Z, lpar_dbl_set *l, fpar_dbl_set *f)
{
	//X is the homogenous solution, Y the pertubation, Z the adjoint pertubation
	

	double b = Y[pos0].EI*Z[pos0].EI + Y[pos0].ER*Z[pos0].ER + Y[pos0].G*Z[pos0].G + Y[pos0].J*Z[pos0].J + Y[pos0].Q*Z[pos0].Q;

	for(long r = -dim1; r < 0; r++)
	{
		b+=Y[pos0+r].EI*(Z[pos0+r+dim1].EI*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*cosf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) + Z[pos0+r+dim1].ER*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*sinf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw)) + Y[pos0+r].ER*(-Z[pos0+r+dim1].EI*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*sinf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) + Z[pos0+r+dim1].ER*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*cosf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw)) + Y[pos0+r].G*(Z[pos0+r+dim1].EI*(-0.5*X[pos0+r+dim1].EI*l->ag*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*sinf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) + 0.5*X[pos0+r+dim1].EI*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*cosf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) - 0.5*X[pos0+r+dim1].ER*l->ag*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*cosf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) - 0.5*X[pos0+r+dim1].ER*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*sinf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw)) + Z[pos0+r+dim1].ER*(0.5*X[pos0+r+dim1].EI*l->ag*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*cosf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) + 0.5*X[pos0+r+dim1].EI*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*sinf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) - 0.5*X[pos0+r+dim1].ER*l->ag*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*sinf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) + 0.5*X[pos0+r+dim1].ER*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*cosf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw))) + Y[pos0+r].Q*(Z[pos0+r+dim1].EI*(0.5*X[pos0+r+dim1].EI*l->aq*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*sinf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) - 0.5*X[pos0+r+dim1].EI*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*cosf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) + 0.5*X[pos0+r+dim1].ER*l->aq*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*cosf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) + 0.5*X[pos0+r+dim1].ER*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*sinf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw)) + Z[pos0+r+dim1].ER*(-0.5*X[pos0+r+dim1].EI*l->aq*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*cosf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) - 0.5*X[pos0+r+dim1].EI*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*sinf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) + 0.5*X[pos0+r+dim1].ER*l->aq*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*sinf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw) - 0.5*X[pos0+r+dim1].ER*l->g*l->sqrtkap*expf(0.5*X[pos0+r+dim1].G - 0.5*X[pos0+r+dim1].Q)*cosf(0.5*X[pos0+r+dim1].G*l->ag - 0.5*X[pos0+r+dim1].Q*l->aq + l->T*l->dw)));
	}

	for(long r = -dim2; r < 0; r++)
	{
		b+=2*X[pos0+r+dim2].EI*Y[pos0+r].EI*Z[pos0+r+dim2].J*f->K*f->wLP + 2*X[pos0+r+dim2].ER*Y[pos0+r].ER*Z[pos0+r+dim2].J*f->K*f->wLP;
	}

	return b;
}/*REPLACE END*/
