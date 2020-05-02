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


/*

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



*/



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
}
