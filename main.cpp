#include <cmath>
#include <vector>
#include <array>
#include <complex>

#include <chrono>

#include <iostream>
#include <iomanip>

#include <fstream>
#include <sstream>
#include <string>

#include <algorithm>
#include <thread>

#include "timer.hpp"
#include "parameter.hpp"
#include "variable.hpp"
#include "integrate.hpp"


using namespace std;


int main(int argc, char* argv[])
{
	timer total;
	
	
	allpar_set *AP = new allpar_set("JAU15", "TAU1", "quick");
	
	AP->check_cmd_line(argc, argv);
	AP->cout_pars(AP->collect());
	


	
	integrator intor(AP);
	
	timeseries TS;
	TS = intor.integrate("complex");

	
	total.stop();
	total.print_elaps("s");
	
	return 0;
}
