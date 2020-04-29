#ifndef __INTEGRATE_HPP_INCLUDED__
#define __INTEGRATE_HPP_INCLUDED__

#include "../hdr/parameter.hpp"
#include "../hdr/variable.hpp"

//###########################################

#define LOOKUP_SAMPLES 1000000
#define LOOKUP_STEP 0.0001
#define LOOKUP_INVERS_STEP 1/0.0001
#define LOOKUP_MAX LOOKUP_STEP*LOOKUP_SAMPLES


static double LOOKUP_EXP_Y0[LOOKUP_SAMPLES];
static double LOOKUP_EXP_S[LOOKUP_SAMPLES];

static double LOOKUP_SIN_Y0[LOOKUP_SAMPLES];
static double LOOKUP_SIN_S[LOOKUP_SAMPLES];

static double LOOKUP_COS_Y0[LOOKUP_SAMPLES];
static double LOOKUP_COS_S[LOOKUP_SAMPLES];

inline double expf(double x)
{
	if(x >= LOOKUP_MAX)
	{
		std::cout << "exp out of bounds" << std::endl;
		exit(1);
	}
	
	if(x < 0)
	{
		return 1.0/expf(-x);
	}
	return LOOKUP_EXP_S[(int)(x*LOOKUP_INVERS_STEP)]*x*LOOKUP_INVERS_STEP + LOOKUP_EXP_Y0[(int)(x*LOOKUP_INVERS_STEP)];
}

inline double sinf(double x)
{
	if(x >= LOOKUP_MAX)
	{
		std::cout << "sin out of bounds" << std::endl;
		exit(1);
	}
	
	if(x < 0)
	{
		return -sinf(-x);
	}
	return LOOKUP_SIN_S[(int)(x*LOOKUP_INVERS_STEP)]*x*LOOKUP_INVERS_STEP + LOOKUP_SIN_Y0[(int)(x*LOOKUP_INVERS_STEP)];
}

inline double cosf(double x)
{
	if(x >= LOOKUP_MAX)
	{
		std::cout << "cos out of bounds" << std::endl;
		exit(1);
	}
	
	if(x < 0)
	{
		return cosf(-x);
	}
	return LOOKUP_COS_S[(int)(x*LOOKUP_INVERS_STEP)]*x*LOOKUP_INVERS_STEP + LOOKUP_COS_Y0[(int)(x*LOOKUP_INVERS_STEP)];
}



//###########################################


class initial_condition
{
public:
	std::vector<var> X;
	double t;
	
	initial_condition();
};

//###########################################

class timeseries
{
public:
	std::vector<var> X;
	std::vector<double> I;
	std::vector<double> t;
	
	timeseries();
};

//###########################################

class integrator
{
public:
	std::vector<varC> C;
	std::vector<var> R;
	allpar_set *AP;
	long it;
	long dim1;
	long dim2;
	long pos0;
	long pos1;
	long pos2;
	double time;
	
	integrator(allpar_set *AP);
	
	timeseries integrate(std::string opt);
	
	void initialize(std::vector<varC> &C, std::string opt);
	void initialize(std::vector<var> &R, std::string opt);
	void initialize(std::vector<varC> &C, initial_condition ic);
	void initialize(std::vector<var> &R, initial_condition ic);
	
	varC derive_full(varC&, varC&, varC&, lpar_dbl_set*, fpar_dbl_set*);
	var derive_ret(var&, var&, var&, var&, var&, var&, lpar_dbl_set*, fpar_dbl_set*);
	var derive_adj(var&, var&, var&, var&, var&, var&, lpar_dbl_set*, fpar_dbl_set*);
	double bilinear_step(std::vector<var>&, std::vector<var>&, std::vector<var>&, lpar_dbl_set*, fpar_dbl_set*);


	
	double test(double);
};

//###########################################


//varC derive_full(varC &X, varC &XT, varC &Xtau, lpar_dbl_set &l, fpar_dbl_set &f);
//var derive_full(var &X, var &XT, var &Xtau, lpar_dbl_set &l, fpar_dbl_set &f);

#endif
