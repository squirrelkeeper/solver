#include <iostream>
#include <iomanip>
#include <fstream>

#include <sstream>
#include <string>

#include <cmath>
#include <complex>
#include <vector>
#include <random>

#include <algorithm>
#include <chrono>
#include <ctime>

#include "timer.hpp"
#include "parameter.hpp"
#include "variable.hpp"
#include "timeseries.hpp"


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

	


	
	integrator IN(AP);	
	
 	timeseries TS(AP);

	TS = IN.integrate_simple_TS_noise("g");
	
	TS.write_file("template_ts_ba");
        
	time_total.stop();
	time_total.print_elaps();
	
	
	return 0;
}

