#include <algorithm>
#include <complex>
#include <cmath>
#include <iostream>
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



initial_condition::initial_condition()
{
	
}

//###########################################

timeseries::timeseries()
{
	
}

//###########################################

integrator::integrator(allpar_set *AP)
{
	this->AP = AP;
	
	lpar_dbl_set lp(*this->AP);
	fpar_dbl_set fp(*this->AP);
	ipar_dbl_set ip(*this->AP);

	it = (long)(ip.int_time / ip.dt);

	if(lp.T < fp.tau)
	{
		dim1 = (long)(floor(lp.T / ip.dt));
		dim2 = (long)(floor(fp.tau / ip.dt));
		
		pos0 = dim2-1;
		pos1 = pos0 - dim1;
		pos2 = 0;
	}
	else
	{
		dim1 = (long)(floor(fp.tau / ip.dt));
		dim2 = (long)(floor(lp.T / ip.dt));

		pos0 = dim2-1;
		pos2 = pos0 - dim1;
		pos1 = 0;
	}
	dim2++;


	
	time = 0.0;
}

void integrator::initialize(vector<varC> &C, string opt)
{
	if(opt == "ones")
	{
		for(unsigned long i = 0; i < C.size(); i++)
		{
			C[i].E = 1.0 + 1.0i;
			C[i].G = 1.0;
			C[i].Q = 1.0;
			C[i].J = 1.0;
		}
	}
}


timeseries integrator::integrate(string opt)
{
	timeseries TS;
	
	if(opt == "complex")
	{
		LOOKUP_EXP_INIT();
		LOOKUP_COS_INIT();
		LOOKUP_SIN_INIT();
		
		lpar_dbl_set lp(*this->AP);
		fpar_dbl_set fp(*this->AP);
		ipar_dbl_set ip(*this->AP);

		this->C.resize(dim2);
	
	
		initialize(this->C, "ones");
		
		long res_size = long(ip.out_time / ip.dt + 0.5);
		
		TS.X.resize(res_size);
		TS.t.resize(res_size);
		TS.I.resize(res_size);
		
	
		for(long i = 0; i < (it - res_size); i++)
		{
			varC dC = derive_full(C[pos0], C[pos1], C[pos2], &lp, &fp);

			C[pos2] = C[pos0] + dC * ip.dt;
			
		
			time += ip.dt;
			
			pos0 = pos2;
			pos1 = (pos1 + 1) % dim2;
			pos2 = (pos2 + 1) % dim2;
			
		}
	
		for(long i = 0; i < res_size; i++)
		{
			varC dC = derive_full(this->C[pos0], this->C[pos1], this->C[pos2], &lp, &fp);

			this->C[pos2] = this->C[pos0] + (dC * ip.dt);

			TS.X[i] = C[pos2];			
			TS.t[i] = time;
			TS.I[i] = norm(C[pos2].E);

			
			
			time += ip.dt;
			
			pos0 = pos2;
			pos1 = (pos1 + 1) % dim1;
			pos2 = (pos2 + 1) % dim2;
			
			 
		}
	
	}
	else if(opt == "linearized")
	{
		
	}
	else if(opt == "adjoint")
	{
		
	}
	else
	{
		cout << "err003" << endl;
	}
	
	
	return TS;
}



double integrator::test(double x)
{
	LOOKUP_EXP_INIT();
	LOOKUP_SIN_INIT();
	LOOKUP_COS_INIT();

	return LOOKUP_EXP_S[0];
}

/*REPLACE START*/
varC integrator::derive_full(varC &X, varC &XT, varC &Xtau, lpar_dbl_set *l, fpar_dbl_set *f)
{
	varC d;

	d.E = -1.0i*l->g*l->sqrtkap*expf(0.5*XT.G)*expf(-0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + 1.0*l->T*l->dw)*XT.E.real() + l->g*l->sqrtkap*expf(0.5*XT.G)*expf(-0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + 1.0*l->T*l->dw)*XT.E.imag() + l->g*l->sqrtkap*expf(0.5*XT.G)*expf(-0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + 1.0*l->T*l->dw)*XT.E.real() + 1.0i*l->g*l->sqrtkap*expf(0.5*XT.G)*expf(-0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + 1.0*l->T*l->dw)*XT.E.imag() - l->g*X.E.real() - 1.0i*l->g*X.E.imag();

	d.G = -X.G*l->gg + l->Jg - expf(X.G - X.Q)*norm(X.E) + expf(-X.Q)*norm(X.E);

	d.Q = -X.J*X.Q + X.J*l->q0 - X.Q*l->gq + l->gq*l->q0 - l->rs*norm(X.E) + l->rs*expf(-X.Q)*norm(X.E);

	d.J = -X.J*f->wLP + f->K*f->wLP*norm(Xtau.E);


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


var integrator::derive_adj(var &X, var &XT, var &Xtau, var &Z, var &ZT, var &Ztau, lpar_dbl_set *l, fpar_dbl_set *f)
{
	//X is the homogenous solution, Z the adjoint pertubation
	var d;

	d.ER = -2*X.ER*Z.G*expf(X.G - X.Q) + 2*X.ER*Z.G*expf(-X.Q) - 2*X.ER*Z.Q*l->rs + 2*X.ER*Z.Q*l->rs*expf(-X.Q) + 2*Xtau.ER*Ztau.J*f->K*f->wLP - Z.ER*l->g - ZT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + ZT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw);

	d.EI = -2*X.EI*Z.G*expf(X.G - X.Q) + 2*X.EI*Z.G*expf(-X.Q) - 2*X.EI*Z.Q*l->rs + 2*X.EI*Z.Q*l->rs*expf(-X.Q) + 2*Xtau.EI*Ztau.J*f->K*f->wLP - Z.EI*l->g + ZT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + ZT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw);

	d.G = -X.EI*X.EI*Z.G*expf(X.G - X.Q) - X.ER*X.ER*Z.G*expf(X.G - X.Q) - 0.5*XT.EI*ZT.EI*l->ag*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.EI*ZT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.EI*ZT.ER*l->ag*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.EI*ZT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.ER*ZT.EI*l->ag*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.ER*ZT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.ER*ZT.ER*l->ag*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.ER*ZT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - Z.G*l->gg;

	d.Q = X.EI*X.EI*Z.G*expf(X.G - X.Q) - X.EI*X.EI*Z.G*expf(-X.Q) - X.EI*X.EI*Z.Q*l->rs*expf(-X.Q) + X.ER*X.ER*Z.G*expf(X.G - X.Q) - X.ER*X.ER*Z.G*expf(-X.Q) - X.ER*X.ER*Z.Q*l->rs*expf(-X.Q) - X.J*Z.Q + 0.5*XT.EI*ZT.EI*l->aq*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.EI*ZT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.EI*ZT.ER*l->aq*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.EI*ZT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.ER*ZT.EI*l->aq*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.ER*ZT.EI*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) + 0.5*XT.ER*ZT.ER*l->aq*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*sinf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - 0.5*XT.ER*ZT.ER*l->g*l->sqrtkap*expf(0.5*XT.G - 0.5*XT.Q)*cosf(0.5*XT.G*l->ag - 0.5*XT.Q*l->aq + l->T*l->dw) - Z.Q*l->gq;

	d.J = -X.Q*Z.Q - Z.J*f->wLP + Z.Q*l->q0;


	return d;
}


double integrator::bilinear_step(vector<var> &X, vector<var> &Y, vector<var> &Z, lpar_dbl_set *l, fpar_dbl_set *f)
{
	//X is the homogenous solution, Y the pertubation, Z the adjoint pertubation
	double b = 0.0;
	long t = this->pos0;
	long T = this->dim1;
	long tau = this->dim2;

	b+=Y[t].EI*Z[t].EI + Y[t].ER*Z[t].ER + Y[t].G*Z[t].G + Y[t].J*Z[t].J + Y[t].Q*Z[t].Q;

	for(long r = -T; r < 0; r++)
	{
		b+=Y[t+r].EI*(Z[t+r+T].EI*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*cosf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) + Z[t+r+T].ER*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*sinf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw)) + Y[t+r].ER*(-Z[t+r+T].EI*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*sinf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) + Z[t+r+T].ER*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*cosf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw)) + Y[t+r].G*(Z[t+r+T].EI*(-0.5*X[t+r].EI*l->ag*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*sinf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) + 0.5*X[t+r].EI*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*cosf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) - 0.5*X[t+r].ER*l->ag*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*cosf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) - 0.5*X[t+r].ER*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*sinf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw)) + Z[t+r+T].ER*(0.5*X[t+r].EI*l->ag*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*cosf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) + 0.5*X[t+r].EI*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*sinf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) - 0.5*X[t+r].ER*l->ag*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*sinf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) + 0.5*X[t+r].ER*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*cosf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw))) + Y[t+r].Q*(Z[t+r+T].EI*(0.5*X[t+r].EI*l->aq*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*sinf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) - 0.5*X[t+r].EI*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*cosf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) + 0.5*X[t+r].ER*l->aq*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*cosf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) + 0.5*X[t+r].ER*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*sinf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw)) + Z[t+r+T].ER*(-0.5*X[t+r].EI*l->aq*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*cosf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) - 0.5*X[t+r].EI*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*sinf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) + 0.5*X[t+r].ER*l->aq*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*sinf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw) - 0.5*X[t+r].ER*l->g*l->sqrtkap*expf(0.5*X[t+r].G - 0.5*X[t+r].Q)*cosf(0.5*X[t+r].G*l->ag - 0.5*X[t+r].Q*l->aq + l->T*l->dw)));
	}

	for(long r = -tau; r < 0; r++)
	{
		b+=2*X[t+r].EI*Y[t+r].EI*Z[t+r+tau].J*f->K*f->wLP + 2*X[t+r].ER*Y[t+r].ER*Z[t+r+tau].J*f->K*f->wLP;
	}

	return b;
}
/*REPLACE END*/
