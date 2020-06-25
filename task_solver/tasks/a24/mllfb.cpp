#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include <algorithm>
#include <complex>
#include <iomanip>

#include<boost/random.hpp>



using namespace std;


#define EXP_SAMPLES 100000
#define EXP_STEP 0.001
#define EXP_STEP_invers 1000.0 // == 1.0/EXP_STEP
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

#define SIN_SAMPLES 100000
#define SIN_STEP 0.001
#define SIN_STEP_invers 1000.0
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



double Expf(double );
double Sinf(double );
double Cosf(double );

//declare function intensity
double intens(double &,double &);

double derivER(double &,double &, double &,double &, double &,double &, double &,double &, double &,double &, double &,double &, double &,double &, double &,double &, double &,double &, double &, double &);
double derivEI(double &,double &, double &,double &, double &,double &, double &,double &, double &,double &, double &,double &, double &,double &, double &,double &, double &,double &, double &, double &);

double derivG(double &,double &, double &,double &, double &);
double derivQ(double &,double &, double &,double &, double &,double &);



int main(int argc, const char* argv[] )
{



	init_exp();
	init_sin();
	init_cos();

	


	unsigned iterations  = 3e7;
	unsigned iterations2 = 3e7;
  

	double h = 1e-4;
	
	double ag = 0;
	double aq = 0;
	
	double g = 66.5;
	
	double gg = 0.025;
	double Jg = 3.0;
	
	double gq = 1.875;
	double q0 = 4.0;
	double Jq = gq*q0;

    double T = 1.0;
	double sqrtkap = sqrt(0.1);
	double rs = 25.0;
	
	double wLP = 16.0*M_PI;
	double tau = 1.0;
	double K = 0.0;

	double D = 0.2;


	double A1ic = 2.0;
	double A2ic = 0.0;
	double A3ic = 1.0;
	double A4ic = 1.0;
	double A5ic = 1.0;

	double B1ic = 0.0;
	double B2ic = 1.0;
	double B3ic = 0.0;
	double B4ic = 0.0;
	double B5ic = 0.0;

	double ERic = 0.4;
	double EIic = 0.0;
	double Gic  = 4.0;
	double Qic  = 1.0;
	double Jic  = 1.0;



	if(argc > 0)
	{
		for(int i=0; i < argc; i++)
		{
			if(string(argv[i]) == string("-K") && argc>i+1)
			{
				K = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-tau") && argc>i+1)
			{
				tau = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-D") && argc>i+1)
			{
				D = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-Jg") && argc>i+1)
			{
				Jg = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-it") && argc>i+1)
			{
				iterations = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-it2") && argc>i+1)
			{
				iterations2 = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-g") && argc>i+1)
			{
				g = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-gg") && argc>i+1)
			{
				gg = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-gq") && argc>i+1)
			{
				gq = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-q0") && argc>i+1)
			{
				q0 = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-kap") && argc>i+1)
			{
				sqrtkap = sqrt(atof(argv[i+1]));
			}
			if(string(argv[i]) == string("-h") && argc>i+1)
			{
				h = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-ag") && argc>i+1)
			{
				ag = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-aq") && argc>i+1)
			{
				aq = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-a1ic") && argc>i+1)
			{
				A1ic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-a2ic") && argc>i+1)
			{
				A2ic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-a3ic") && argc>i+1)
			{
				A3ic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-a4ic") && argc>i+1)
			{
				A4ic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-a5ic") && argc>i+1)
			{
				A5ic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-b1ic") && argc>i+1)
			{
				B1ic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-b2ic") && argc>i+1)
			{
				B2ic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-b3ic") && argc>i+1)
			{
				B3ic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-b4ic") && argc>i+1)
			{
				B4ic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-b5ic") && argc>i+1)
			{
				B5ic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-eric") && argc>i+1)
			{
				ERic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-eiic") && argc>i+1)
			{
				EIic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-gic") && argc>i+1)
			{
				Gic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-qic") && argc>i+1)
			{
				Qic = atof(argv[i+1]);
			}
			if(string(argv[i]) == string("-jic") && argc>i+1)
			{
				Jic = atof(argv[i+1]);
			}
		}
	}


	int i1_max = 50;
	int i2_max = 50;
	
	ofstream out1;
	out1.open("a24.dat");
	
	out1 << "#tau" << '\t' << "K" << '\t' << "jit" << endl;

	for(int i1 = 0; i1 < i1_max; i1++)
	{
		tau = 0 + (double)(i1)*1.25/((double)i1_max);
	
		for(int i2 = 0; i2 < i2_max; i2++)
		{
			K = 0 + (double)(i2)*4.0/((double)i2_max);
			
			double jitter1 = 0.0;

			
			
			long double Time = 0.0;

			// the following part of the program calculates and out puts time series of E, G and Q
		
			unsigned int dim=floor((T+tau)/h+1);
			unsigned int dimT=floor(T/h);

			vector<double> A1(dim,A1ic);
			vector<double> A2(dim,A2ic);
			vector<double> A3(dim,A3ic);
			vector<double> A4(dim,A4ic);
			vector<double> A5(dim,A5ic);

			vector<double> B1(dim,B1ic);
			vector<double> B2(dim,B2ic);
			vector<double> B3(dim,B3ic);
			vector<double> B4(dim,B4ic);
			vector<double> B5(dim,B5ic);

			vector<double> ER(dim,ERic);
			vector<double> EI(dim,EIic);
			vector<double> G(dim,Gic);
			vector<double> Q(dim,Qic);
			vector<double> J(dim,Jic);

			unsigned int currPos=dim-1;
			unsigned int intDel=currPos-dimT;
			unsigned int delay=0;

			double intensity=intens(ER[currPos],EI[currPos]);

			vector<double> peakfinding(3,0.0);
			vector<double> peaktimefinding(3,0.0);

			bool peakrepeat=false;


			vector<double> maxima(0);
			maxima.reserve(40);
			vector<double> maxtime(0);
			maxtime.reserve(40);


			
			for(unsigned long int i = 1; i <= iterations; i++)
			{
				double ERnew = ER[currPos] + 
				h*(
					- g * ER[currPos]
					+ g * sqrtkap * Expf(0.5*(G[intDel]-Q[intDel]))
					* (
						ER[intDel]*Cosf(0.5*ag*G[intDel]-0.5*aq*Q[intDel])
						+EI[intDel]*Sinf(0.5*ag*G[intDel]-0.5*aq*Q[intDel])
					)
				);
				
				double EInew = EI[currPos] +
				h*(
					- g * EI[currPos]
					+ g * sqrtkap * Expf(0.5*(G[intDel]-Q[intDel]))
					* (
						EI[intDel]*Cosf(0.5*ag*G[intDel]-0.5*aq*Q[intDel])
						-ER[intDel]*Sinf(0.5*ag*G[intDel]-0.5*aq*Q[intDel])
					)
				);
				
				double Gnew = G[currPos] +
				h*(
					Jg - gg * G[currPos] - Expf(-Q[currPos]) * (Expf(G[currPos])-1.0) * intensity
				);
				
				double Qnew = Q[currPos] +
				h*(
					(q0 - Q[currPos])*(gq + J[currPos]) - rs*Expf(-Q[currPos]) * (Expf(Q[currPos])-1) * intensity
				);
				
				double Jnew = J[currPos] +
				h*(
					wLP * K * (ER[delay]*ER[delay] + EI[delay]*EI[delay]) - wLP * J[currPos]
				);
				
				
				ER[delay] = ERnew;
				EI[delay] = EInew;
				G[delay]  = Gnew;
				Q[delay]  = Qnew;
				J[delay]  = Jnew;
				
				Time += h;
				
				intensity = intens(ERnew,EInew);
				
				currPos=delay;                
				delay=(delay+1) % dim;
				intDel=(intDel+1) % dim;

				double Qtotal = Qnew+abs(log(sqrtkap*sqrtkap*(1+K)*(1+K)));
				double gtotal = Gnew - Qtotal;
				double absE = sqrt(intensity);
				

				if(i>iterations-50*T/h-2)
				{
					peakfinding.erase (peakfinding.begin());
					peaktimefinding.erase (peaktimefinding.begin());
					peakfinding.push_back (gtotal);
					peaktimefinding.push_back (Time);
				}		
				
				
				if(i>(iterations-50*T/h) and gtotal>0)
				{
					if(peakfinding[1]>peakfinding[2] and peakfinding[1]>peakfinding[0])
					{
						for(unsigned int j=0; j< maxima.size();j++ )
						{
							if (peakfinding[1] >= maxima[j]*0.99 and peakfinding[1] <= maxima[j]*1.01)
							{
								peakrepeat=true;
							}
						}

						if(peakrepeat==false)
						{
							// Fitting three data points nearest to the maximum to a parabola to obtain a more accurate peak height and time  f(t)=a(t-tmax)^2+c
							double t1=-h;
							double t2=0;
							double t3=h;

							double p1=peakfinding[0];
							double p2=peakfinding[1];
							double p3=peakfinding[2];


							double c=((p3*(t1 - t2)*(t1 - t2) - p1*(t2 - t3)*(t2 - t3))*(p3*(t1 - t2)*(t1 - t2) - p1*(t2 - t3)*(t2 - t3)) - 2*p2*(p3*(t1 - t2)*(t1 - t2) + p1*(t2 - t3)*(t2 - t3))*(t1 - t3)*(t1 - t3) + p2*p2*(t1 - t3)*(t1 - t3)*(t1 - t3)*(t1 - t3))/(4*(t1 - t2)*(t1 - t3)*(t2 - t3)*(p3*(-t1 + t2) + p2*(t1 - t3) + p1*(-t2 + t3)));

							double tmax=peaktimefinding[1]+(p2*t1*t1 - p3*t1*t1 - p1*t2*t2 + p3*t2*t2+ p1*t3*t3 - p2*t3*t3)/(2*(p2*t1 - p3*t1 - p1*t2 + p3*t2 + p1*t3 - p2*t3));
							maxima.push_back (c);
							maxtime.push_back (tmax);
						}

						if(peakrepeat==true)
						{  
							// Fitting three data points nearest to the maximum to a parabola to obtain a more accurate peak height and time  f(t)=a(t-tmax)^2+c
							double t1=-h;
							double t2=0;
							double t3=h;

							double p1=peakfinding[0];
							double p2=peakfinding[1];
							double p3=peakfinding[2];


							double c=((p3*(t1 - t2)*(t1 - t2) - p1*(t2 - t3)*(t2 - t3))*(p3*(t1 - t2)*(t1 - t2) - p1*(t2 - t3)*(t2 - t3)) - 2*p2*(p3*(t1 - t2)*(t1 - t2) + p1*(t2 - t3)*(t2 - t3))*(t1 - t3)*(t1 - t3) + p2*p2*(t1 - t3)*(t1 - t3)*(t1 - t3)*(t1 - t3))/(4*(t1 - t2)*(t1 - t3)*(t2 - t3)*(p3*(-t1 + t2) + p2*(t1 - t3) + p1*(-t2 + t3)));

							double tmax=peaktimefinding[1]+(p2*t1*t1 - p3*t1*t1 - p1*t2*t2 + p3*t2*t2+ p1*t3*t3 - p2*t3*t3)/(2*(p2*t1 - p3*t1 - p1*t2 + p3*t2 + p1*t3 - p2*t3));
							maxtime.push_back (tmax);
						}
						peakrepeat=false;
					}
				}
			}
			

			
			
			
			if(maxtime.size()>0)
			{
				double dT=(maxtime[maxtime.size()-1]-maxtime[0])/(maxtime.size()-1);
				double Tisi=dT;
				double variance=0;

				for(int i=0;i<(maxtime.size()-1);i++)
				{
					variance+=(maxtime[i+1]-maxtime[i]-Tisi)*(maxtime[i+1]-maxtime[i]-Tisi)/(maxtime.size()-1);
				}

				unsigned int roundtripdim=floor(Tisi/h);
				
				if(Tisi<1.02)
				{
					vector<double> ER0(0);
					// vectors containging periodic solution for homogenous equations. increasing index is going forward in time
					ER0.reserve(iterations2+4*dim);
					vector<double> EI0(0);
					EI0.reserve(iterations2+4*dim);
					vector<double> G0(0);
					G0.reserve(iterations2+4*dim);
					vector<double> Q0(0);
					Q0.reserve(iterations2+4*dim);
					vector<double> J0(0);
					J0.reserve(iterations2+4*dim);
					
				
		//			ofstream out1;
		//			out1.open("Psi0.ts.dat");


					for(unsigned long int i = 1; i <= iterations2+4*dim; i++)
					{
						double ERnew = ER[currPos] + 
						h*(
							- g * ER[currPos]
							+ g * sqrtkap * Expf(0.5*(G[intDel]-Q[intDel]))
							* (
								ER[intDel]*Cosf(0.5*ag*G[intDel]-0.5*aq*Q[intDel])
								+EI[intDel]*Sinf(0.5*ag*G[intDel]-0.5*aq*Q[intDel])
							)
						);
						
						double EInew = EI[currPos] +
						h*(
							- g * EI[currPos]
							+ g * sqrtkap * Expf(0.5*(G[intDel]-Q[intDel]))
							* (
								EI[intDel]*Cosf(0.5*ag*G[intDel]-0.5*aq*Q[intDel])
								-ER[intDel]*Sinf(0.5*ag*G[intDel]-0.5*aq*Q[intDel])
							)
						);
						
						double Gnew = G[currPos] +
						h*(
							Jg - gg * G[currPos] - Expf(-Q[currPos]) * (Expf(G[currPos])-1.0) * intensity
						);
						
						double Qnew = Q[currPos] +
						h*(
							(q0 - Q[currPos])*(gq + J[currPos]) - rs*Expf(-Q[currPos]) * (Expf(Q[currPos])-1) * intensity
						);
						
						double Jnew = J[currPos] +
						h*(
							wLP * K * (ER[delay]*ER[delay] + EI[delay]*EI[delay]) - wLP * J[currPos]
						);

						ER[delay] = ERnew;
						EI[delay] = EInew;
						G[delay]  = Gnew;
						Q[delay]  = Qnew;
						J[delay]  = Jnew;
						
						Time += h;

						intensity = intens(ERnew,EInew);
						currPos = delay;
						delay = (delay+1) % dim;
						intDel = (intDel+1) % dim;

						ER0.push_back(ERnew);
						EI0.push_back(EInew);
						G0.push_back(Gnew);
						Q0.push_back(Qnew);
						J0.push_back(Jnew);
					
		/*
						if(i > iterations2+4*dim - 40000)
						{						
							out1 << ERnew << '\t';
							out1 << EInew << '\t';
							out1 << Gnew << '\t';
							out1 << Qnew << '\t';
							out1 << Jnew << '\t';
							out1 << endl;
						}
		*/

					}
					
		//			out1.close();



					currPos=dim-1;
					intDel=currPos-dimT;
					delay=0;
						
					intensity=intens(ER[currPos],EI[currPos]);

					unsigned int homoPos=iterations2+4*dim-1;
					double intensity0=intens(ER0[homoPos],EI0[homoPos]);

		/*					
					ofstream out2;
					out2.open("PsiA.ts.dat");
					
					ofstream out3;
					out3.open("PsiB.ts.dat");
		*/
							
							
							
					for(unsigned long int i=1;i<=iterations2;i++)
					{
					
						intensity0=intens(ER0[homoPos],EI0[homoPos]);		
					
					
						double A1new = A1[currPos];
						A1new += h * (
							- g * A1[currPos]
							- Expf(-Q0[homoPos])
								* (Expf(G0[homoPos])-1.0) 
								* 2.0 * ER0[homoPos] * A3[currPos]
							- rs * Expf(-Q0[homoPos])
								* (Expf(Q0[homoPos])-1.0)
								* 2.0 * ER0[homoPos] * A4[currPos]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * A1[intDel]
							- g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * A2[intDel]
							+ 2.0 * wLP * K * ER0[homoPos] * A5[delay]
						);
						
						double A2new = A2[currPos];
						A2new += h * (
							- g * A2[currPos]
							- Expf(-Q0[homoPos])
								* (Expf(G0[homoPos])-1.0)
								* 2.0 * EI0[homoPos] * A3[currPos]
							- rs * Expf(-Q0[homoPos])
								* (Expf(Q0[homoPos])-1.0)
								* 2.0 * EI0[homoPos] * A4[currPos]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * A1[intDel]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos]) 
								* Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * A2[intDel]
							+ 2.0 * wLP * K * EI0[homoPos] * A5[delay]
						);
						
						
						double A3new = A3[currPos];
						A3new += h * (
							- gg * A3[currPos]
							- Expf(-Q0[homoPos])*Expf(G0[homoPos])*intensity0*A3[currPos]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)
								+(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)*ag
							) * A1[intDel]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								-1.0*(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)*ag
								+(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)
							) * A2[intDel]
						);
						
						double A4new = A4[currPos];
						A4[currPos] += h * (
							- Expf(-Q0[homoPos])*(1-Expf(G0[homoPos]))*intensity0 * A3[currPos]
							- (gq + J0[homoPos]) * A4[currPos]
							- rs * Expf(-Q0[homoPos]) * intensity0 * A4[currPos]
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								-1.0*(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)
								-(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos])*aq
							) * A1[intDel]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)*aq
								-(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)
							) * A2[intDel]
						);
									
						
						double A5new = A5[currPos];
						A5[currPos] += h * (
							A4[currPos] * (q0 - Q0[homoPos])- wLP * A5[currPos]
						);


						
						
						double B1new = B1[currPos];
						B1new += h * (
							- g * B1[currPos]
							- Expf(-Q0[homoPos])
								* (Expf(G0[homoPos])-1.0) 
								* 2.0 * ER0[homoPos] * B3[currPos]
							- rs * Expf(-Q0[homoPos])
								* (Expf(Q0[homoPos])-1.0)
								* 2.0 * ER0[homoPos] * B4[currPos]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * B1[intDel]
							- g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * B2[intDel]
							+ 2.0 * wLP * K * ER0[homoPos] * B5[delay]
						);


						
						double B2new = B2[currPos];
						B2new += h * (
							- g * B2[currPos]
							- Expf(-Q0[homoPos])
								* (Expf(G0[homoPos])-1.0)
								* 2.0 * EI0[homoPos] * B3[currPos]
							- rs * Expf(-Q0[homoPos])
								* (Expf(Q0[homoPos])-1.0)
								* 2.0 * EI0[homoPos] * B4[currPos]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * B1[intDel]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos]) 
								* Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * B2[intDel]
							+ 2.0 * wLP * K * EI0[homoPos] * B5[delay]
						);
						
						double B3new = B3[currPos];
						B3new += h * (
							- gg * B3[currPos]
							- Expf(-Q0[homoPos])*Expf(G0[homoPos])*intensity0*B3[currPos]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)
								+(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)*ag
							) * B1[intDel]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								-1.0*(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)*ag
								+(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)
							) * B2[intDel]
						);

						double B4new = B4[currPos];
						B4[currPos] += h * (
							- Expf(-Q0[homoPos])*(1-Expf(G0[homoPos]))*intensity0 * B3[currPos]
							- (gq + J0[homoPos]) * B4[currPos]
							- rs * Expf(-Q0[homoPos]) * intensity0 * B4[currPos]
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								-1.0*(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)
								-(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos])*aq
							) * B1[intDel]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)*aq
								-(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)
							) * B2[intDel]
						);

						double B5new = B5[currPos];
						B5[currPos] += h * (
							B4[currPos] * (q0 - Q0[homoPos])- wLP * B5[currPos]
						);




						

						
						A1[delay] = A1new;
						A2[delay] = A2new;
						A3[delay] = A3new;
						A4[delay] = A4new;
						A5[delay] = A5new;
						
						B1[delay] = B1new;
						B2[delay] = B2new;
						B3[delay] = B3new;
						B4[delay] = B4new;
						B5[delay] = B5new;


						Time += h;
				
						currPos=delay;
						delay=(delay+1) % dim;
						intDel=(intDel+1) % dim;

						homoPos--;
						

								
		/*						
						if(i > iterations2 - 40000)
						{
							out2 << A1new << '\t';
							out2 << A2new << '\t';
							out2 << A3new << '\t';
							out2 << A4new << '\t';
							out2 << A5new << '\t';
							
							out2 << endl;
							
							
							out3 << B1new << '\t';
							out3 << B2new << '\t';
							out3 << B3new << '\t';
							out3 << B4new << '\t';
							out3 << B5new << '\t';
							
							out3 << endl;
						}
		*/						
								
					}

		//			out2.close();
		//			out3.close();
					

					// vectors containging periodic solution for homogenous equations. increasing index is going forward in time
					vector<double> dER0(0);   
					dER0.reserve(dim*2+100);
					vector<double> dEI0(0);
					dEI0.reserve(dim*2+100);
					vector<double> dG0(0);
					dG0.reserve(dim*2+100);
					vector<double> dQ0(0);
					dQ0.reserve(dim*2+100);
					vector<double> dJ0(0);
					dJ0.reserve(dim*2+100);

					// vectors containging periodic solution for homogenous equations. increasing index is going forward in time
					vector<double> negEI0(0);   
					negEI0.reserve(dim*2+100);
					vector<double> ER0final(0);
					ER0final.reserve(dim*2+100);


					// vectors containging periodic solution for homogenous equations. increasing index is going forward in time
					vector<double> A1final(0);
					A1final.reserve(dim*2+100);
					vector<double> A2final(0);
					A2final.reserve(dim*2+100);
					vector<double> A3final(0);
					A3final.reserve(dim*2+100);
					vector<double> A4final(0);
					A4final.reserve(dim*2+100);
					vector<double> A5final(0);
					A5final.reserve(dim*2+100);

						
					// vectors containging periodic solution for homogenous equations. increasing index is going forward in time
					vector<double> B1final(0);
					B1final.reserve(dim*2+100);
					vector<double> B2final(0);
					B2final.reserve(dim*2+100);
					vector<double> B3final(0);
					B3final.reserve(dim*2+100);
					vector<double> B4final(0);
					B4final.reserve(dim*2+100);
					vector<double> B5final(0);
					B5final.reserve(dim*2+100);

					
					// vectors containging periodic variation of the non-zero entries of the C matrix - needed to the bilinear form
					vector<double> C11(0);
					C11.reserve(dim*2+100);
					vector<double> C12(0);
					C12.reserve(dim*2+100);
					vector<double> C13(0);
					C13.reserve(dim*2+100);
					vector<double> C14(0);
					C14.reserve(dim*2+100);
					vector<double> C21(0);
					C21.reserve(dim*2+100);
					vector<double> C22(0);
					C22.reserve(dim*2+100);
					vector<double> C23(0);
					C23.reserve(dim*2+100);
					vector<double> C24(0);
					C24.reserve(dim*2+100);

					// vectors containging periodic variation of the non-zero entries of the C matrix - needed to the bilinear form
					vector<double> D51(0);
					D51.reserve(dim*2+100);
					vector<double> D52(0);
					D52.reserve(dim*2+100);

					double homoPosintdelay = homoPos - dimT;
					double homoPosdelay = homoPos - (dim - 1);

					double Yhomo1=0;
					double Yhomo2=0;
					double Zhomo1=0;
					double Zhomo2=0;
					



		//			ofstream out4;
		//			out4.open("dPsi.ts.dat");
							
							
					for(unsigned long int i=1;i<=dim*2+100;i++)
					{
						double A1new = A1[currPos];
						A1new += h * (
							- g * A1[currPos]
							- Expf(-Q0[homoPos])
								* (Expf(G0[homoPos])-1.0) 
								* 2.0 * ER0[homoPos] * A3[currPos]
							- rs * Expf(-Q0[homoPos])
								* (Expf(Q0[homoPos])-1.0)
								* 2.0 * ER0[homoPos] * A4[currPos]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * A1[intDel]
							- g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * A2[intDel]
							+ 2.0 * wLP * K * ER0[homoPos] * A5[delay]
						);

						double A2new = A2[currPos];
						A2new += h * (
							- g * A2[currPos]
							- Expf(-Q0[homoPos])
								* (Expf(G0[homoPos])-1.0)
								* 2.0 * EI0[homoPos] * A3[currPos]
							- rs * Expf(-Q0[homoPos])
								* (Expf(Q0[homoPos])-1.0)
								* 2.0 * EI0[homoPos] * A4[currPos]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * A1[intDel]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos]) 
								* Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * A2[intDel]
							+ 2.0 * wLP * K * EI0[homoPos] * A5[delay]
						);
						
						double A3new = A3[currPos];
						A3new += h * (
							- gg * A3[currPos]
							- Expf(-Q0[homoPos])*Expf(G0[homoPos])*intensity0*A3[currPos]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)
								+(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)*ag
							) * A1[intDel]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								-1.0*(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)*ag
								+(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)
							) * A2[intDel]
						);

						double A4new = A4[currPos];
						A4[currPos] += h * (
							- Expf(-Q0[homoPos])*(1-Expf(G0[homoPos]))*intensity0 * A3[currPos]
							- (gq + J0[homoPos]) * A4[currPos]
							- rs * Expf(-Q0[homoPos]) * intensity0 * A4[currPos]
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								-1.0*(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)
								-(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos])*aq
							) * A1[intDel]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)*aq
								-(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)
							) * A2[intDel]
						);

						double A5new = A5[currPos];
						A5[currPos] += h * (
							A4[currPos] * (q0 - Q0[homoPos])- wLP * A5[currPos]
						);


						double B1new = B1[currPos];
						B1new += h * (
							- g * B1[currPos]
							- Expf(-Q0[homoPos])
								* (Expf(G0[homoPos])-1.0) 
								* 2.0 * ER0[homoPos] * B3[currPos]
							- rs * Expf(-Q0[homoPos])
								* (Expf(Q0[homoPos])-1.0)
								* 2.0 * ER0[homoPos] * B4[currPos]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * B1[intDel]
							- g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * B2[intDel]
							+ 2.0 * wLP * K * ER0[homoPos] * B5[delay]
						);

						double B2new = B2[currPos];
						B2new += h * (
							- g * B2[currPos]
							- Expf(-Q0[homoPos])
								* (Expf(G0[homoPos])-1.0)
								* 2.0 * EI0[homoPos] * B3[currPos]
							- rs * Expf(-Q0[homoPos])
								* (Expf(Q0[homoPos])-1.0)
								* 2.0 * EI0[homoPos] * B4[currPos]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])
								* Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * B1[intDel]
							+ g * sqrtkap
								* Expf(0.5*G0[homoPos]-0.5*Q0[homoPos]) 
								* Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]) * B2[intDel]
							+ 2.0 * wLP * K * EI0[homoPos] * B5[delay]
						);
						
						double B3new = B3[currPos];
						B3new += h * (
							- gg * B3[currPos]
							- Expf(-Q0[homoPos])*Expf(G0[homoPos])*intensity0*B3[currPos]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)
								+(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)*ag
							) * B1[intDel]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								-1.0*(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)*ag
								+(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)
							) * B2[intDel]
						);

						double B4new = B4[currPos];
						B4[currPos] += h * (
							- Expf(-Q0[homoPos])*(1-Expf(G0[homoPos]))*intensity0 * B3[currPos]
							- (gq + J0[homoPos]) * B4[currPos]
							- rs * Expf(-Q0[homoPos]) * intensity0 * B4[currPos]
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								-1.0*(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)
								-(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos])*aq
							) * B1[intDel]
							
							+ g * sqrtkap * 0.5 * Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(
								(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
									+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
								)*aq
								-(
									Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]
									-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]
								)
							) * B2[intDel]
						);

						double B5new = B5[currPos];
						B5[currPos] += h * (
							B4[currPos] * (q0 - Q0[homoPos])- wLP * B5[currPos]
						);

						
						double dERnew =  
						(
							- g * ER0[homoPos]
							+ g * sqrtkap * Expf(0.5*(G0[homoPosintdelay]-Q0[homoPosintdelay]))
							* (
								ER0[homoPosintdelay]*Cosf(0.5*ag*G0[homoPosintdelay]-0.5*aq*Q0[homoPosintdelay])
								+EI0[homoPosintdelay]*Sinf(0.5*ag*G0[homoPosintdelay]-0.5*aq*Q0[homoPosintdelay])
							)
						);
						
						double dEInew = 
						(
							- g * EI0[homoPos]
							+ g * sqrtkap * Expf(0.5*(G0[homoPosintdelay]-Q0[homoPosintdelay]))
							* (
								EI0[homoPosintdelay]*Cosf(0.5*ag*G0[homoPosintdelay]-0.5*aq*Q0[homoPosintdelay])
								-ER0[homoPosintdelay]*Sinf(0.5*ag*G0[homoPosintdelay]-0.5*aq*Q0[homoPosintdelay])
							)
						);
						
						double dGnew = 
						(
							Jg - gg * G0[homoPos] - Expf(-Q0[homoPos]) * (Expf(G0[homoPos])-1.0) * intensity0
						);
						
						double dQnew = 
						(
							(q0 - Q0[homoPos])*(gq + J0[homoPos]) - rs*Expf(-Q0[homoPos]) * (Expf(Q0[homoPos])-1) * intensity0
						);
						
						double dJnew = 
						(
							wLP * K * (ER0[homoPosdelay]*ER0[homoPosdelay] + EI0[homoPosdelay]*EI0[homoPosdelay]) - wLP * J0[homoPos]
						);
		

		/*						
						out4 << dERnew << '\t';
						out4 << dEInew << '\t';
						out4 << dGnew << '\t';
						out4 << dQnew << '\t';
						out4 << dJnew << '\t';
						out4 << endl;
		*/						
								

						double C11new=g*sqrtkap*Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]);
						double C12new=g*sqrtkap*Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]);
						double C22new=g*sqrtkap*Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]);
						double C21new=-g*sqrtkap*Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos]);

						double C13new=g*sqrtkap*0.5*Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*((Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos])+(Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos])*ag);
						
						double C14new=g*sqrtkap*0.5*Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(-1*(Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos])-(Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos])*aq);
						
						double C23new=g*sqrtkap*0.5*Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*(-1*(Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos])*ag+(Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]));
						
						double C24new=g*sqrtkap*0.5*Expf(0.5*G0[homoPos]-0.5*Q0[homoPos])*((Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]+Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos])*aq-(Cosf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*EI0[homoPos]-Sinf(0.5*ag*G0[homoPos]-0.5*aq*Q0[homoPos])*ER0[homoPos]));

						double D51new = 2.0 * K * wLP * ER0[homoPos];
						double D52new = 2.0 * K * wLP * EI0[homoPos];

						negEI0.push_back(-EI0[homoPos]);
						ER0final.push_back(ER0[homoPos]);  
						// these are in reverse time


						A1[delay]=A1new;
						A2[delay]=A2new;
						A3[delay]=A3new;
						A4[delay]=A4new;
						A5[delay]=A5new;
						
						B1[delay]=B1new;
						B2[delay]=B2new;
						B3[delay]=B3new;
						B4[delay]=B4new;
						B5[delay]=B5new;
						
						Time+=h;
						
						
						currPos=delay;
						delay=(delay+1) % dim;
						intDel=(intDel+1) % dim;

						homoPos=homoPos-1;
						homoPosintdelay=homoPos-dimT;
						homoPosdelay=homoPos-(dim-1);
						intensity0=intens(ER0[homoPos],EI0[homoPos]);

						// these are in reverse time
						dER0.push_back(dERnew);
						dEI0.push_back(dEInew);
						dG0.push_back(dGnew);
						dQ0.push_back(dQnew);
						dJ0.push_back(dJnew);

						// these are in reverse time
						A1final.push_back(A1new);
						A2final.push_back(A2new);
						A3final.push_back(A3new);
						A4final.push_back(A4new);
						A5final.push_back(A5new);

						// these are in reverse time
						B1final.push_back(B1new); 
						B2final.push_back(B2new);
						B3final.push_back(B3new);
						B4final.push_back(B4new);
						B5final.push_back(B5new);

						// these are in reverse time
						C11.push_back(C11new);
						C12.push_back(C12new);
						C13.push_back(C13new);
						C14.push_back(C14new);
						C21.push_back(C21new);
						C22.push_back(C22new);
						C23.push_back(C23new);
						C24.push_back(C24new);

						// these are in reverse time
						D51.push_back(D51new);
						D52.push_back(D52new);
					}
							
							
		//			out4.close();
							



					for(long int i=dimT;i>=0;i--)
					{

						if(i==0 or i==dimT)
						{
							//cout << dim+50-dimT+i<< '\t'<< dim+50+i<< '\n';

							Zhomo1 += 0.5*h * (
								B1final[dim+50-dimT+i] * (
									C11[dim+50+i] * negEI0[dim+50+i]
									+ C12[dim+50+i] * ER0final[dim+50+i]
								)
								+ B2final[dim+50-dimT+i] * (
									C21[dim+50+i]*negEI0[dim+50+i]
									+ C22[dim+50+i]*ER0final[dim+50+i]
								)
							);

							Zhomo2 += 0.5*h * (
								B1final[dim+50-dimT+i] * (
									C11[dim+50+i]*dER0[dim+50+i]
									+C12[dim+50+i]*dEI0[dim+50+i]
									+C13[dim+50+i]*dG0[dim+50+i]
									+C14[dim+50+i]*dQ0[dim+50+i]
								)
								+B2final[dim+50-dimT+i] * (
									C21[dim+50+i]*dER0[dim+50+i]
									+C22[dim+50+i]*dEI0[dim+50+i]
									+C23[dim+50+i]*dG0[dim+50+i]
									+C24[dim+50+i]*dQ0[dim+50+i]
								)
							);

							Yhomo1 += 0.5* h * (
								A1final[dim+50-dimT+i] * (
									C11[dim+50+i]*negEI0[dim+50+i]
									+C12[dim+50+i]*ER0final[dim+50+i]
								)
								+A2final[dim+50-dimT+i] * (
									C21[dim+50+i]*negEI0[dim+50+i]
									+C22[dim+50+i]*ER0final[dim+50+i]
								)
							);

							Yhomo2 += 0.5* h * (
								A1final[dim+50-dimT+i]*(C11[dim+50+i]*dER0[dim+50+i]+C12[dim+50+i]*dEI0[dim+50+i]+C13[dim+50+i]*dG0[dim+50+i]+C14[dim+50+i]*dQ0[dim+50+i])+A2final[dim+50-dimT+i]*(C21[dim+50+i]*dER0[dim+50+i]+C22[dim+50+i]*dEI0[dim+50+i]+C23[dim+50+i]*dG0[dim+50+i]+C24[dim+50+i]*dQ0[dim+50+i]));

						}
						
						else
						{

							Zhomo1+=h*(B1final[dim+50-dimT+i]*(C11[dim+50+i]*negEI0[dim+50+i]+C12[dim+50+i]*ER0final[dim+50+i])+B2final[dim+50-dimT+i]*(C21[dim+50+i]*negEI0[dim+50+i]+C22[dim+50+i]*ER0final[dim+50+i]));

							Zhomo2+=h*(B1final[dim+50-dimT+i]*(C11[dim+50+i]*dER0[dim+50+i]+C12[dim+50+i]*dEI0[dim+50+i]+C13[dim+50+i]*dG0[dim+50+i]+C14[dim+50+i]*dQ0[dim+50+i])+B2final[dim+50-dimT+i]*(C21[dim+50+i]*dER0[dim+50+i]+C22[dim+50+i]*dEI0[dim+50+i]+C23[dim+50+i]*dG0[dim+50+i]+C24[dim+50+i]*dQ0[dim+50+i]));

							Yhomo1+=h*(A1final[dim+50-dimT+i]*(C11[dim+50+i]*negEI0[dim+50+i]+C12[dim+50+i]*ER0final[dim+50+i])+A2final[dim+50-dimT+i]*(C21[dim+50+i]*negEI0[dim+50+i]+C22[dim+50+i]*ER0final[dim+50+i]));

							Yhomo2+=h*(A1final[dim+50-dimT+i]*(C11[dim+50+i]*dER0[dim+50+i]+C12[dim+50+i]*dEI0[dim+50+i]+C13[dim+50+i]*dG0[dim+50+i]+C14[dim+50+i]*dQ0[dim+50+i])+A2final[dim+50-dimT+i]*(C21[dim+50+i]*dER0[dim+50+i]+C22[dim+50+i]*dEI0[dim+50+i]+C23[dim+50+i]*dG0[dim+50+i]+C24[dim+50+i]*dQ0[dim+50+i]));
						}
					}
					
					for(long int i=(dim-1);i>=0;i--)
					{

						if(i==0 or i==(dim-1))
						{
							Zhomo1 += 0.5*h * (
								B5final[dim+50-(dim-1)+i] * (
									D51[dim+50+i]*ER0final[dim+50+i]
									+D52[dim+50+i]*negEI0[dim+50+i]
								)
							);					
							
							Zhomo2 += 0.5*h * (
								B5final[dim+50-(dim-1)+i] * (
									D51[dim+50+i]*dER0[dim+50+i]
									+D52[dim+50+i]*dEI0[dim+50+i]
								)
							);

							Yhomo1 += 0.5*h * (
								A5final[dim+50-(dim-1)+i] * (
									D51[dim+50+i]*ER0final[dim+50+i]
									+D52[dim+50+i]*negEI0[dim+50+i]
								)
							);
							
							Yhomo2 += 0.5*h * (
								A5final[dim+50-(dim-1)+i] * (
									D51[dim+50+i]*dER0[dim+50+i]
									+D52[dim+50+i]*dEI0[dim+50+i]
								)
							);
							
						}
						else
						{
							Zhomo1 += h * (
								B5final[dim+50-(dim-1)+i] * (
									D51[dim+50+i]*ER0final[dim+50+i]
									+D52[dim+50+i]*negEI0[dim+50+i]
								)
							);					
							
							Zhomo2 += h * (
								B5final[dim+50-(dim-1)+i] * (
									D51[dim+50+i]*dER0[dim+50+i]
									+D52[dim+50+i]*dEI0[dim+50+i]
								)
							);

							Yhomo1 += h * (
								A5final[dim+50-(dim-1)+i] * (
									D51[dim+50+i]*ER0final[dim+50+i]
									+D52[dim+50+i]*negEI0[dim+50+i]
								)
							);
							
							Yhomo2 += h * (
								A5final[dim+50-(dim-1)+i] * (
									D51[dim+50+i]*dER0[dim+50+i]
									+D52[dim+50+i]*dEI0[dim+50+i]
								)
							);
						}
					}
						
					Zhomo1=Zhomo1+B1final[dim+50]*negEI0[dim+50]+B2final[dim+50]*ER0final[dim+50];

					Zhomo2=Zhomo2+B1final[dim+50]*dER0[dim+50]+B2final[dim+50]*dEI0[dim+50]+B3final[dim+50]*dG0[dim+50]+B4final[dim+50]*dQ0[dim+50]+B5final[dim+50]*dJ0[dim+50];
					
					Yhomo1=Yhomo1+A1final[dim+50]*negEI0[dim+50]+A2final[dim+50]*ER0final[dim+50];

					Yhomo2=Yhomo2+A1final[dim+50]*dER0[dim+50]+A2final[dim+50]*dEI0[dim+50]+A3final[dim+50]*dG0[dim+50]+A4final[dim+50]*dQ0[dim+50]+A5final[dim+50]*dJ0[dim+50];


					double c1ad1=1/(Zhomo1-(Yhomo1*Zhomo2)/Yhomo2);
					double c2ad1=1/(Yhomo1-(Zhomo1*Yhomo2)/Zhomo2);

					double c1ad2=1/(Zhomo2-(Yhomo2*Zhomo1)/Yhomo1);
					double c2ad2=1/(Yhomo2-(Zhomo2*Yhomo1)/Zhomo1);

					
					vector<double> ad11(0);
					ad11.reserve(roundtripdim);
					vector<double> ad12(0);
					ad12.reserve(roundtripdim);
					vector<double> ad13(0);
					ad13.reserve(roundtripdim);
					vector<double> ad14(0);
					ad14.reserve(roundtripdim);
					vector<double> ad21(0);
					ad21.reserve(roundtripdim);
					vector<double> ad22(0);
					ad22.reserve(roundtripdim);
					vector<double> ad23(0);
					ad23.reserve(roundtripdim);
					vector<double> ad24(0);
					ad24.reserve(roundtripdim);

					
					Time=0;

					for(unsigned long int i=1;i<=roundtripdim;i++)
					{

						ad11.push_back(c1ad1*B1final[i-1]+c2ad1*A1final[i-1]);
						ad12.push_back(c1ad1*B2final[i-1]+c2ad1*A2final[i-1]);
						ad13.push_back(c1ad1*B3final[i-1]+c2ad1*A3final[i-1]);
						ad14.push_back(c1ad1*B4final[i-1]+c2ad1*A4final[i-1]);

						ad21.push_back(c1ad2*B1final[i-1]+c2ad2*A1final[i-1]);
						ad22.push_back(c1ad2*B2final[i-1]+c2ad2*A2final[i-1]);
						ad23.push_back(c1ad2*B3final[i-1]+c2ad2*A3final[i-1]);
						ad24.push_back(c1ad2*B4final[i-1]+c2ad2*A4final[i-1]);

						Time+=h;
					}

					jitter1=0;
						
						


					for(unsigned long int i=1;i<=roundtripdim;i++)
					{
						if(i==1 or i==roundtripdim)
						{
							jitter1+=0.5*h*(ad21[i-1]*ad21[i-1]+ad22[i-1]*ad22[i-1]);
						}
						else
						{
							jitter1+=h*(ad21[i-1]*ad21[i-1]+ad22[i-1]*ad22[i-1]);
						}
					}
						
//					cout << "tj: " << sqrt(jitter1)*0.2*25000.0 << endl;

					
					peakfinding.clear();
					peaktimefinding.clear();

					maxima.clear();
					maxtime.clear();
					ad11.clear();

					ad12.clear();
					ad13.clear();
					ad14.clear();
					ad21.clear();
					ad22.clear();
					ad23.clear();
					ad24.clear();

					negEI0.clear();
					ER0final.clear();


					dER0.clear();
					dEI0.clear();
					dG0.clear();
					dQ0.clear();
					dJ0.clear();


					A1final.clear();
					A2final.clear();
					A3final.clear();
					A4final.clear();
					A5final.clear();

					B1final.clear();
					B2final.clear();
					B3final.clear();
					B4final.clear();
					B5final.clear();

					C11.clear();
					C12.clear();
					C13.clear();
					C14.clear();
					C21.clear();
					C22.clear();
					C23.clear();
					C24.clear();

					D51.clear();
					D52.clear();

					ER0.clear();
					EI0.clear();
					G0.clear();
					Q0.clear();
					J0.clear();

					ER.clear();
					EI.clear();
					G.clear();
					Q.clear();
					J.clear();

					A1.clear();
					A2.clear();
					A3.clear();
					A4.clear();
					A5.clear();

					B1.clear();
					B2.clear();
					B3.clear();
					B4.clear();
					B5.clear();
					
					
				}
					
				else
				{
					double jitter1=0;
					double variance=0;		
				}
			}
			else
			{
				double jitter1=0;
				double Tisi=0;
				double variance=0;
			}
	
			out1 << tau << '\t';
			out1 << K << '\t';
			out1 << sqrt(jitter1)*0.2*25000.0 << '\t';
			out1 << endl;
			
		}
	}
	
	out1.close();

	
	return 0;
}



inline double Expf( double ex ) {
    if(ex < 0) return 1.0/Expf(-ex);
    double EE = ex*EXP_STEP_invers;
    int EEi = (int)EE;
    return exp_slope[EEi]*EE + exp_y[EEi];
}

inline double Sinf( double ex ) {
    if(ex < 0) return -Sinf(-ex);
    double EE = ex*SIN_STEP_invers;
    int EEi = (int)EE;
    return sin_slope[EEi]*EE + sin_y[EEi];
}

inline double Cosf( double ex ) {
    if(ex < 0) return Cosf(-ex);
    double EE = ex*SIN_STEP_invers;
    int EEi = (int)EE;
    return cos_slope[EEi]*EE + cos_y[EEi];
}


inline double intens(double & Er, double & Ei)
{
    return Er*Er+Ei*Ei;
}
