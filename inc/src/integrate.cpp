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

integrator::integrator(allpar_set *ap_init)
{
	AP = ap_init;
	
	lpar_dbl_set lp(*this->AP);
	fpar_dbl_set fp(*this->AP);
	ipar_dbl_set ip(*this->AP);


	it = (long)(ip.int_time / ip.dt);

	double tau1 = AP->smaller_delay();
	double tau2 = AP->larger_delay();

	dim1 = (long)(floor(tau1 / ip.dt));
	dim2 = (long)(floor(tau2 / ip.dt)+1);
	
	pos0 = dim2 - 1;
	pos1 = pos0 - dim1;
	pos2 = 0;

	Time = 0.0;
	
	X.resize(dim2);
}


void integrator::integrate(vector<double> &result_Time, vector<var> &result_X, string opt)
{
	
	double icer = 0.4;
	double icei = 0.0;
	double icg  = 4.0;
	double icq  = 1.0;
	double icj  = 0.0;
	
	
	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();


	lpar_dbl_set *lp = new lpar_dbl_set(*AP);
	fpar_dbl_set *fp = new fpar_dbl_set(*AP);
	ipar_dbl_set *ip = new ipar_dbl_set(*AP);
	
	for(long i = 0; i < dim2; i++)
	{
		X[i].ER = icer;
		X[i].EI = icei;
		X[i].G = icg;
		X[i].Q = icq;
		X[i].J = icj;
	}
	
	double h = ip->dt;
	long it = ip->out_time/ip->dt;
	long rtout = ip->out_time/ip->dt;

	for(long i=1; i <= it; i++)
	{
		
		dX = derive(X[pos0], X[pos1], X[pos2], lp, fp);
		
		Xnew.ER = X[pos0].ER + h * dX.ER;
		Xnew.EI = X[pos0].EI + h * dX.EI;
		Xnew.G  = X[pos0].G  + h * dX.G;
		Xnew.Q  = X[pos0].Q  + h * dX.Q;
		Xnew.J  = X[pos0].J  + h * dX.J;
		
		
		X[pos2].ER = Xnew.ER;
		X[pos2].EI = Xnew.EI;
		X[pos2].G = Xnew.G;
		X[pos2].Q = Xnew.Q;
		X[pos2].J = Xnew.J;
		Time += h;
		
		pos0 = pos2;
		pos2 = (pos2+1) % dim2;
		pos1 = (pos1+1) % dim2;

		
		

		if(i > it - rtout)
		{
			result_Time.push_back(Time);
			result_X.push_back(Xnew);
		}
	}
	
	
	
	
	
}

var integrator::derive(var &X, var &XT, var &Xtau, lpar_dbl_set *l, fpar_dbl_set *f)
{
	var d;

	d.ER = -X.ER*l->g + XT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + XT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw);

	d.EI = -X.EI*l->g + XT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - XT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw);

	d.G = -X.EI*X.EI*expf(X.G - X.Q) + X.EI*X.EI*expf(-X.Q) - X.ER*X.ER*expf(X.G - X.Q) + X.ER*X.ER*expf(-X.Q) - X.G*l->gg + l->Jg;

	d.Q = -X.EI*X.EI*l->rs + X.EI*X.EI*l->rs*expf(-X.Q) - X.ER*X.ER*l->rs + X.ER*X.ER*l->rs*expf(-X.Q) - X.J*X.Q + X.J*l->q0 - X.Q*l->gq + l->gq*l->q0;

	d.J = -X.J*f->wLP + Xtau.EI*Xtau.EI*f->K*f->wLP + Xtau.ER*Xtau.ER*f->K*f->wLP;


	return d;
}


