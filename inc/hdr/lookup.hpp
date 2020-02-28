#ifndef __LOOKUP_HPP_INCLUDED__
#define __LOOKUP_HPP_INCLUDED__

#define LOOKUP_SAMPLES 1000000

class lookup_exp
{
public:
	int samples = LOOKUP_SAMPLES;
	double step = 1e-4;
	double step_inv = 1.0/step;
	double y[LOOKUP_SAMPLES];
	double slope[LOOKUP_SAMPLES];
	
	lookup_exp();
	lookup_exp(int init_samples, double init_step);

	inline double f(double &x);
};

//###########################################

class lookup_cos
{
public:
	int samples = LOOKUP_SAMPLES;
	double step = 1e-4;
	double step_inv = 1.0/step;
	double y[LOOKUP_SAMPLES];
	double slope[LOOKUP_SAMPLES];
	
	lookup_cos();
	lookup_cos(int init_samples, double init_step);
	inline double f(double &x);
};

//###########################################

class lookup_sin
{
public:
	int samples = LOOKUP_SAMPLES;
	double step = 1e-4;
	double step_inv = 1.0/step;
	double y[LOOKUP_SAMPLES];
	double slope[LOOKUP_SAMPLES];
	
	lookup_sin();
	lookup_sin(int init_samples, double init_step);
	inline double f(double &x);
};

#endif
