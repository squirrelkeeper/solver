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

#include "integrate.hpp"

typedef unsigned int uint;
typedef unsigned long int ulint;

using namespace std;



int main(int argc, char* argv[])
{
	timer time_total;

	allpar_set *AP = new allpar_set("JAU15", "TAU1", "quick");
	AP->check_cmd_line(argc, argv);


//######################### initial values for history arrays   ###########################

	double icer = 0.4;
	double icei = 0.0;
	double icg  = 4.0;
	double icq  = 1.0;
	double icj  = 0.0;


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





	vector<double> result_Time;	
	vector<var> result_X;
	
	
	integrator IN(AP);
	IN.integrate(result_Time, result_X, "27");
	
	
	
	
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

