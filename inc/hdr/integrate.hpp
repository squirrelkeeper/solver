#ifndef __INTEGRATE_HPP_INCLUDED__
#define __INTEGRATE_HPP_INCLUDED__

#include "../hdr/parameter.hpp"
#include "../hdr/variable.hpp"
#include "../hdr/timeseries.hpp"
#include "../hdr/initial_con.hpp"



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

class integrator
{
public:
	std::vector<var> X;
	var dX;
	var Xnew;
	allpar_set AP;
	long it;
	long dim1;
	long dim2;
	long pos0;
	long pos1;
	long pos2;
	double Time;
	int rea;
	int max_rea;
	bool noise;
	
	integrator(allpar_set);
	
	void initialize(initial_con IC);
	void initialize(std::string opt);

	
	timeseries integrate();
	std::tuple<timeseries,ts_evaluation> integrate_analysis(std::string);
	std::tuple<timeseries,timeseries,timeseries,ts_evaluation> integrate_neutral_modes_analysis(std::string);
	
	
	timeseries integrate_noise();
	std::tuple<timeseries, pp_evaluation> integrate_noise_analysis(std::string);
	std::vector<double> integrate_noise_conc_analysis();

	timeseries integrate_adj(timeseries);
	double bilinear_prod(timeseries, timeseries, timeseries);
	
	
	
	timeseries integrate_ret(timeseries);
	
	std::vector<timeseries> integrate_get_neutral_modes();

	double bilinear_one_step(timeseries, timeseries, timeseries);

	var derive_real(var&, var&, var&, lpar_dbl_set*, fpar_dbl_set*);
	var derive_ret(var&, var&, var&, var&, var&, var&, lpar_dbl_set*, fpar_dbl_set*);
	var derive_adj(var&, var&, var&, var&, lpar_dbl_set*, fpar_dbl_set*);

	double bilinear_step(std::vector<var>&, std::vector<var>&, std::vector<var>&, lpar_dbl_set*, fpar_dbl_set*);

	double InterpolQuadExtrPos(double, double, double, double, double, double);
	double InterpolQuadExtrVal(double, double, double, double, double, double, double);	
};

#endif
