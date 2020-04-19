#include <algorithm>
#include <complex>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../hdr/integrate.hpp"


using namespace std;


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

	dim1 = (long)(lp.T / ip.dt  + 0.5);
	dim2 = (long)(fp.tau / ip.dt + 0.5);

	pos0 = dim1 - 1;
	pos1 = pos0 - dim1;
	pos2 = 0;
	
	t = 0.0;
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
			varC dC = derive_full(this->C[pos0], this->C[pos1], this->C[pos2], &lp, &fp);

			this->C[pos2] = this->C[pos0] + (dC * ip.dt);
			
			pos0 = pos2;
			pos1 = (pos1 + 1) % dim1;
			pos2 = (pos2 + 1) % dim2;
		}

		for(long i = 0; i < res_size; i++)
		{
			
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


varC integrator::derive_full(varC &X, varC &XT, varC &Xtau, lpar_dbl_set *l, fpar_dbl_set *f)
{
	varC dX;
	
	dX.E = 1.0i;
	
	
	
	
	return dX;
}

double integrator::test(double x)
{
	return 5.0;
}

