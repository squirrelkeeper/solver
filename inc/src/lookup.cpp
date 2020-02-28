#include <cmath>
#include <complex>
#include <iostream>

#include "../hdr/lookup.hpp"


using namespace std;


lookup_exp::lookup_exp()
{
	for(int i = 0; i < samples; i++)
	{
		y[i] = exp(i*step);
	}
	for(int i = 0; i < (samples-1); i++)
	{
		slope[i] = y[i+1] - y[i];
		y[i] = y[i] - i * (y[i+1] - y[i]);
    }
}

lookup_exp::lookup_exp(int init_samples, double init_step)
{
	for(int i = 0; i < samples; i++)
	{
		y[i] = exp(i*step);
	}
	for(int i = 0; i < (samples-1); i++)
	{
		slope[i] = y[i+1] - y[i];
		y[i] = y[i] - i * (y[i+1] - y[i]);
    }
}

inline double lookup_exp::f(double &x)
{
	if(x >= samples * step)
	{
		cout << "exp out of bounds" << endl;
		exit(1);
	}
	
    if(x < 0.0)
	{
		x = -x;
		return 1.0/f(x);
	}
	
    double EE = x * step_inv;
    int EEi = (int)(EE);
	
    return slope[EEi] * EE + y[EEi];
}

//###########################################

lookup_cos::lookup_cos()
{
	for(int i = 0; i < samples; i++)
	{
		y[i] = cos(i*step);
	}
	for(int i = 0; i < (samples-1); i++)
	{
		slope[i] = y[i+1] - y[i];
		y[i] = y[i] - i * (y[i+1] - y[i]);
    }
}

inline double lookup_cos::f(double &x)
{
	if(x >= samples * step)
	{
		cout << "cos out of bounds" << endl;
		exit(1);
	}
	
    if(x < 0.0)
	{
		x = -x;
		return f(x);
	}
	
    double EE = x * step_inv;
    int EEi = (int)(EE);
	
    return slope[EEi] * EE + y[EEi];
}

//###########################################

lookup_sin::lookup_sin()
{
	for(int i = 0; i < samples; i++)
	{
		y[i] = sin(i*step);
	}
	for(int i = 0; i < (samples-1); i++)
	{
		slope[i] = y[i+1] - y[i];
		y[i] = y[i] - i * (y[i+1] - y[i]);
    }
}

inline double lookup_sin::f(double &x)
{
	if(x >= samples * step)
	{
		cout << "sin out of bounds" << endl;
		exit(1);
	}
	
    if(x < 0.0)
	{
		x = -x;
		return -f(x);
	}
	
    double EE = x * step_inv;
    int EEi = (int)(EE);
	
    return slope[EEi] * EE + y[EEi];
}
