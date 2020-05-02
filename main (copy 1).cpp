#include <iostream>
#include <iomanip>
#include <fstream>

#include <sstream>
#include <string>

#include <cmath>
#include <complex>
#include <vector>

#include <algorithm>
#include <chrono>
#include <ctime>

#include "timer.hpp"
#include "parameter.hpp"
#include "variable.hpp"







typedef unsigned int uint;
typedef unsigned long int ulint;

using namespace std;


#define EXP_SAMPLES 1000000
#define EXP_STEP 0.0001
#define EXP_STEP_invers 10000.0

#define SIN_SAMPLES 1000000
#define SIN_STEP 0.0001
#define SIN_STEP_invers 10000.0


double exp_y[EXP_SAMPLES], exp_slope[EXP_SAMPLES];
void init_exp() {
    for(int i=0; i<EXP_SAMPLES; i++) {
        exp_y[i] = exp(i*EXP_STEP);
    }
    for(int i=0; i<EXP_SAMPLES-1; i++) {
        exp_slope[i] = exp_y[i+1]-exp_y[i];
        exp_y[i] = exp_y[i]-i*(exp_y[i+1]-exp_y[i]);
    }

}


double sin_y[SIN_SAMPLES], sin_slope[SIN_SAMPLES];
double cos_y[SIN_SAMPLES], cos_slope[SIN_SAMPLES];
void init_sin() {
    for(int i=0; i<SIN_SAMPLES; i++) {
        sin_y[i] = sin(i*SIN_STEP);
    }
    for(int i=0; i<SIN_SAMPLES-1; i++) {
        sin_slope[i] = sin_y[i+1]-sin_y[i];
        sin_y[i] = sin_y[i]-i*(sin_y[i+1]-sin_y[i]);
    }

}
void init_cos() {
    for(int i=0; i<SIN_SAMPLES; i++) {
        cos_y[i] = cos(i*SIN_STEP);
    }
    for(int i=0; i<SIN_SAMPLES-1; i++) {
        cos_slope[i] = cos_y[i+1]-cos_y[i];
        cos_y[i] = cos_y[i]-i*(cos_y[i+1]-cos_y[i]);
    }

}


	                   
                    

double expf(double );
double sinf(double );
double cosf(double );

void deriv(
	var &d,
	var &X,
	var &XT,
	var &Xtau,
	lpar_dbl_set *LP,
	fpar_dbl_set *FP
);





int main(int argc, char* argv[])
{
	timer time_total;

	
	
	init_exp();
	init_sin();
	init_cos();

		
	allpar_set *AP = new allpar_set("JAU15", "TAU1", "quick");
	AP->check_cmd_line(argc, argv);
	
	lpar_dbl_set *lp = new lpar_dbl_set(*AP);
	fpar_dbl_set *fp = new fpar_dbl_set(*AP);
	ipar_dbl_set *ip = new ipar_dbl_set(*AP);

//######################### initial values for history arrays   ###########################

	double icer = 0.4;
	double icei = 0.0;
	double icg  = 4.0;
	double icq  = 1.0;
	double icj  = 0.0;

//############################################ integration parameters ###########################


	double h = ip->dt;                  // time step size
	ulint it = ip->out_time/ip->dt; // must be larger than it_used (lttj), stepnpforfft (fft) and 200000 (timeseries)
	uint rtout = ip->out_time; // for time series, if turnon is false, then this gives the number of round trips that are read out


//############################## reading in of parameter values and flags #########################################
	

	

	
	if(argc > 0) 
	{
		for(int i=0; i < argc; i++)
		{
			if(string(argv[i]) == string("-icer") && argc>i+1)
			{
				icer = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-icei") && argc>i+1)
			{
				icei = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-icg") && argc>i+1)
			{
				icg = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-icq") && argc>i+1)
			{
				icq = atof(argv[i+1]);
			}

		}
	}


//########################### defining history vectors #####################################


	uint dim2 =  floor((lp->T+fp->tau)/h+1);   // dimension of history vector
	uint dim1 = floor(lp->T/h);          // internal delay - used later to specify position of the internal delay

	
	uint pos0 = dim2 - 1;         // indices of the delay terms
	uint posT = pos0 - dim1;
	uint postau = 0;

	

//####################### defining time for the integration ###########################################

	long double Time = 0.0; //Time


/**************** from here on are the various calculations that are carried out if the respective flags are set to true ***********************/

//############################### calculation of a timeseries ######################################//


	vector< double > result_Time;	
	vector<var> result_X;
	
	vector<var> X(dim2);
	var dX;
	var Xnew;
	
	for(ulint i = 0; i < dim2; i++)
	{
		X[i].ER = icer;
		X[i].EI = icei;
		X[i].G = icg;
		X[i].Q = icq;
		X[i].J = icj;
	}
		


	for(ulint i=1; i <= it; i++)
	{
		deriv(dX, X[pos0], X[posT], X[postau], lp, fp);
		
		Xnew.ER = X[pos0].ER + h * dX.ER;
		Xnew.EI = X[pos0].EI + h * dX.EI;
		Xnew.G  = X[pos0].G  + h * dX.G;
		Xnew.Q  = X[pos0].Q  + h * dX.Q;
		Xnew.J  = X[pos0].J  + h * dX.J;
		
		
		X[postau].ER = Xnew.ER;
		X[postau].EI = Xnew.EI;
		X[postau].G = Xnew.G;
		X[postau].Q = Xnew.Q;
		X[postau].J = Xnew.J;
		Time += h;
		
		pos0 = postau;
		postau = (postau+1) % dim2;
		posT = (posT+1) % dim2;

		
		

		if(i > it - rtout * lp->T / h)
		{
			result_Time.push_back(Time);
			result_X.push_back(Xnew);
		}
	}

	
	
	
	
	stringstream paramFileName;
	paramFileName.precision(7);
	
	
	srand (time(NULL));
	int batch_no = rand()%10000;

	paramFileName << "data/template_ts_ba" << batch_no <<".dat";


	ofstream data;
	data.open(paramFileName.str().c_str(),ios::trunc);


	data << "#t" << '\t';
	data << "ER" << '\t';
	data << "EI" << '\t';
	data << "G" << '\t';
	data << "Q" << '\t';
	data << "J" << '\t';
	data << "I" << '\t';
	data << endl;
	
	
	
	for(ulint i = 0; i < result_Time.size(); i++)
	{
		data << setprecision(15);
		data << result_Time[i] << '\t';
		data << result_X[i].ER << '\t';
		data << result_X[i].EI << '\t';
		data << result_X[i].G << '\t';
		data << result_X[i].Q << '\t';
		data << result_X[i].J << '\t';
		data << result_X[i].ER*result_X[i].ER+result_X[i].EI*result_X[i].EI << '\t';
		data << endl;
	}
	
	data.close();
         
        
	time_total.stop();
	time_total.print_elaps();
	
	
	return 0;
}

//######################################################################################################################################################

inline double expf( double ex )
{
	if(ex>=EXP_SAMPLES*EXP_STEP){cout <<" exp out of bounds"<<endl;
	exit(1);}
	if(ex < 0) return 1.0/expf(-ex);
	double EE = ex*EXP_STEP_invers;
	int EEi = (int)EE;
	return exp_slope[EEi]*EE + exp_y[EEi];
}

inline double sinf( double ex )
{
	if(ex>=SIN_SAMPLES*SIN_STEP){ cout <<" sin out of bounds"<<endl;
	exit(1);}
	if(ex < 0) return -sinf(-ex);
	double EE = ex*SIN_STEP_invers;
	int EEi = (int)EE;
	return sin_slope[EEi]*EE + sin_y[EEi];
}

inline double cosf( double ex )
{
	if(ex>=SIN_SAMPLES*SIN_STEP) {cout <<" cos out of bounds"<<endl;
	exit(1);}
	if(ex < 0) return cosf(-ex);
	double EE = ex*SIN_STEP_invers;
	int EEi = (int)EE;
	return cos_slope[EEi]*EE + cos_y[EEi];
}


void deriv(
	var &d,
	var &X,
	var &XT,
	var &Xtau,
	lpar_dbl_set *LP,
	fpar_dbl_set *FP
)
{
	d.ER  = -LP->g * X.ER;
	d.ER +=  LP->g * LP->sqrtkap * expf(0.5*(XT.G-XT.Q)) * cosf(0.5*(-LP->ag*XT.G+LP->aq*XT.Q)+LP->dw*LP->T) * XT.ER;
	d.ER += -LP->g * LP->sqrtkap * expf(0.5*(XT.G-XT.Q)) * sinf(0.5*(-LP->ag*XT.G+LP->aq*XT.Q)+LP->dw*LP->T) * XT.EI;

	d.EI  = -LP->g * X.EI;
	d.EI +=  LP->g * LP->sqrtkap * expf(0.5*(XT.G-XT.Q)) * cosf(0.5*(-LP->ag*XT.G+LP->aq*XT.Q)+LP->dw*LP->T) * XT.EI;
	d.EI +=  LP->g * LP->sqrtkap * expf(0.5*(XT.G-XT.Q)) * sinf(0.5*(-LP->ag*XT.G+LP->aq*XT.Q)+LP->dw*LP->T) * XT.ER;
	
	d.G = LP->Jg - LP->gg * X.G - expf(-X.Q) * (expf(X.G)-1) * (X.ER*X.ER+X.EI*X.EI);
	
	d.Q = (LP->gq + X.J) * (LP->q0 - X.Q) - LP->rs * expf(-X.Q) * (expf(X.Q)-1) * (X.ER*X.ER+X.EI*X.EI);

	d.J = FP->wLP * (FP->K * (Xtau.ER*Xtau.ER + Xtau.EI*Xtau.EI) - X.J);
}
