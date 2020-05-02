#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <algorithm>
#include <chrono>








typedef unsigned int uint;
typedef unsigned long int ulint;

using namespace std;


#define EXP_SAMPLES 1000000
#define EXP_STEP 0.0001
#define EXP_STEP_invers 10000.0 // == 1.0/EXP_STEP
// exp_y = Schnittpunkt mit y-Achse
// exp_slope = Steigung
// Exponentialfunktionergibt sich aus exp(x) = exp_y[(int)xx] +exp_slope[(int)xx] * xx, mit xx := x / EXP_STEP
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

#define SIN_SAMPLES 1000000
#define SIN_STEP 0.0001
#define SIN_STEP_invers 10000.0
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


struct LP
{
	double a;
	double ag;
	double aq;
	double Jg;
	double q0;
	double rs;
	double g;
	double gg;
	double gq;
	double sqrtkap;
	double dw;
	double T;
};

struct FP
{
	double K;
	double tau;
	double sinC;
	double cosC;
	double omegaLP;
};

#define  intstep    double dERnew; double dEInew; double dGnew; double dQnew; double dJnew;\
                    deriv(ER[pos0], ER[posT], ER[postau], EI[pos0], EI[posT], EI[postau], G[pos0], G[posT], G[postau], Q[pos0], Q[posT], Q[postau], J[pos0], J[posT], J[postau], lp, fp, dERnew, dEInew, dGnew, dQnew, dJnew);\
                    double ERnew = ER[pos0] + h * dERnew;\
                    double EInew = EI[pos0] + h * dEInew;\
                    double Gnew  = G[pos0]  + h * dGnew;\
                    double Qnew  = Q[pos0]  + h * dQnew;\
                    double Jnew  = J[pos0]  + h * dJnew;\
                    ER[postau] = ERnew; EI[postau] = EInew; G[postau] = Gnew; Q[postau] = Qnew; J[postau] = Jnew; Time += h;\
                    pos0 = postau;\
                    postau = (postau+1) % dim;\
                    posT = (posT+1) % dim;\



double Expf(double );
double Sinf(double );
double Cosf(double );

//declare functions
void deriv(
	double & ER,
	double & ERintdel,
	double & ERdelay,
	double & EI,
	double & EIintdel,
	double & EIdelay,
	double & G,
	double & Gintdel,
	double & Gdelay,
	double & Q,
	double & Qintdel,
	double & Qdelay,
	double & J,
	double & Jintdel,
	double & Jdelay,
	LP *LP,
	FP *FP, 
	double & dERnew,
	double & dEInew,
	double & dGnew,
	double & dQnew,
	double & dJnew
);




int main(int argc, const char* argv[])
{
	init_exp();  // initial declaration of sin, cos, exp lookup functions
	init_sin();
	init_cos();

//#######################   laser and feedback parameters for JAU15 ##########################

	LP *lp = new LP();

	lp->a  = 0.0;
	lp->ag = 0.0;
	lp->aq = 0.0;
	lp->Jg = 3.0;
	lp->q0 = 4.0;
	lp->rs = 25.0;
	lp->g  = 66.5;
	lp->gg = 0.025;
	lp->gq = 1.875;
	lp->dw = 0.0;
	lp->T  = 1.0;
	lp->sqrtkap = sqrt(0.1);

	FP *fp = new FP();

	fp->tau = 0.1;
	fp->K = 0.0;
	fp->sinC = sin(0.0);
	fp->cosC = cos(0.0);
	fp->omegaLP = 2*M_PI / 1.0;
	

//######################### initial values for history arrays   ###########################

	double icer = 0.4;
	double icei = 0.0;
	double icg  = 4.0;
	double icq  = 1.0;
	double icj  = 0.0;

//############################################ integration parameters ###########################

	ulint it = 1E6; // must be larger than it_used (lttj), stepnpforfft (fft) and 200000 (timeseries)
	double h = 1E-4;                  // time step size

//########################################### calculation specific variables ########################

	uint rtout = ((double)(it))*h; // for time series, if turnon is false, then this gives the number of round trips that are read out


//############################## reading in of parameter values and flags #########################################

	if(argc > 0) 
	{
		for(int i=0; i < argc; i++)
		{
			if((string(argv[i])) == string("-Jg") && argc>i+1)
			{
				lp->Jg = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-g") && argc>i+1)
			{
				lp->g = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-gg") && argc>i+1)
			{
				lp->gg = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-gq") && argc>i+1)
			{
				lp->gq = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-q0") && argc>i+1)
			{
				lp->q0 = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-kap") && argc>i+1)
			{
				lp->sqrtkap = sqrt(atof(argv[i+1]));
			}
			if(string(argv[i]) == string("-a") && argc>i+1)
			{
				lp->a = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-ag") && argc>i+1)
			{
				lp->ag = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-aq") && argc>i+1)
			{
				lp->aq = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-rs") && argc>i+1)
			{
				lp->rs = atof(argv[i+1]);
			}
			if((string(argv[i])) == string("-dw") && argc>i+1)
			{
				lp->dw = atof(argv[i+1]);
			}

			if(string(argv[i]) == string("-K") && argc>i+1)
			{
				fp->K = atof(argv[i+1]);
			}
			if((string(argv[i])) == string("-tau") && argc>i+1)
			{
				fp->tau = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-C") && argc>i+1)
			{
				fp->sinC = sin(atof(argv[i+1]));
				fp->cosC = cos(atof(argv[i+1]));
			}
			if(string(argv[i]) == string("-tauLP") && argc>i+1)
			{
				fp->omegaLP = 2*M_PI/atof(argv[i+1]);
			}

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

			if((string(argv[i])) == string("-rtout") && argc>i+1)
			{
				rtout = atof(argv[i+1]);
			}
			if((string(argv[i])) == string("-it") && argc>i+1)
			{
				it = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-h") && argc>i+1)
			{
				h = atof(argv[i+1]);
			}
		}
	}


//########################### defining history vectors #####################################


	uint dim =  floor((lp->T+fp->tau)/h+1);   // dimension of history vector
	uint dimT = floor(lp->T/h);          // internal delay - used later to specify position of the internal delay

	vector< double > ER(dim, icer);
	vector< double > EI(dim, icei);
	vector< double >  G(dim, icg);
	vector< double >  Q(dim, icq);
	vector< double >  J(dim, icj);
	
	uint pos0 = dim - 1;         // indices of the delay terms
	uint posT = pos0 - dimT;
	uint postau = 0;

	

//####################### defining time for the integration ###########################################

	long double Time = 0.0; //Time


/**************** from here on are the various calculations that are carried out if the respective flags are set to true ***********************/

//############################### calculation of a timeseries ######################################//


	vector< double > result_Time;
	vector< double > result_ER;
	vector< double > result_EI;
	vector< double > result_G;
	vector< double > result_Q;
	vector< double > result_J;

	for(ulint i=1; i <= it; i++)
	{

		
			intstep;
		
		if(i > it - rtout * lp->T / h)
		{
			result_Time.push_back(Time);
			result_ER.push_back(ERnew);
			result_EI.push_back(EInew);
			result_G.push_back(Gnew);
			result_Q.push_back(Qnew);
			result_J.push_back(Jnew);
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
		data << result_ER[i] << '\t';
		data << result_EI[i] << '\t';
		data << result_G[i] << '\t';
		data << result_Q[i] << '\t';
		data << result_J[i] << '\t';
		data << result_ER[i]*result_ER[i]+result_EI[i]*result_EI[i] << '\t';
		data << endl;
	}
	
	data.close();
         
        


	return 0;
}

//######################################################################################################################################################

inline double Expf( double ex )
{
	if(ex>=EXP_SAMPLES*EXP_STEP){cout <<" exp out of bounds"<<endl;
	exit(1);}
	if(ex < 0) return 1.0/Expf(-ex);
	double EE = ex*EXP_STEP_invers;
	int EEi = (int)EE;
	return exp_slope[EEi]*EE + exp_y[EEi];
}

inline double Sinf( double ex )
{
	if(ex>=SIN_SAMPLES*SIN_STEP){ cout <<" sin out of bounds"<<endl;
	exit(1);}
	if(ex < 0) return -Sinf(-ex);
	double EE = ex*SIN_STEP_invers;
	int EEi = (int)EE;
	return sin_slope[EEi]*EE + sin_y[EEi];
}

inline double Cosf( double ex )
{
	if(ex>=SIN_SAMPLES*SIN_STEP) {cout <<" cos out of bounds"<<endl;
	exit(1);}
	if(ex < 0) return Cosf(-ex);
	double EE = ex*SIN_STEP_invers;
	int EEi = (int)EE;
	return cos_slope[EEi]*EE + cos_y[EEi];
}

void deriv(
	double & ER,
	double & ERintdel,
	double & ERdelay,
	double & EI,
	double & EIintdel,
	double & EIdelay,
	double & G,
	double & Gintdel,
	double & Gdelay,
	double & Q,
	double & Qintdel,
	double & Qdelay,
	double & J,
	double & Jintdel,
	double & Jdelay,
	LP *LP,
	FP *FP, 
	double & dERnew,
	double & dEInew,
	double & dGnew,
	double & dQnew,
	double & dJnew	
)
{
	dERnew  = -LP->g * ER + LP->a * EI;
	dERnew +=  LP->g * LP->sqrtkap * Expf(0.5*(Gintdel-Qintdel)) * Cosf(0.5*(-LP->ag*Gintdel+LP->aq*Qintdel)+LP->dw*LP->T) * ERintdel;
	dERnew += -LP->g * LP->sqrtkap * Expf(0.5*(Gintdel-Qintdel)) * Sinf(0.5*(-LP->ag*Gintdel+LP->aq*Qintdel)+LP->dw*LP->T) * EIintdel;

	dEInew  = -LP->g * EI - LP->a * ER;
	dEInew +=  LP->g * LP->sqrtkap * Expf(0.5*(Gintdel-Qintdel)) * Cosf(0.5*(-LP->ag*Gintdel+LP->aq*Qintdel)+LP->dw*LP->T) * EIintdel;
	dEInew +=  LP->g * LP->sqrtkap * Expf(0.5*(Gintdel-Qintdel)) * Sinf(0.5*(-LP->ag*Gintdel+LP->aq*Qintdel)+LP->dw*LP->T) * ERintdel;
	
	dGnew = LP->Jg - LP->gg * G-Expf(-Q) * (Expf(G)-1) * (ER*ER+EI*EI);
	
	dQnew = (LP->gq + J) * (LP->q0 - Q) - LP->rs * Expf(-Q) * (Expf(Q)-1) * (ER*ER+EI*EI);

	dJnew = FP->omegaLP * (FP->K * (ERdelay*ERdelay + EIdelay*EIdelay) - J);
}


