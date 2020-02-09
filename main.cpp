#include <cmath>
#include <vector>
#include <array>

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


using namespace std;


int main(int argc, char* argv[])
{
	timer total;
	
	this_thread::sleep_for(std::chrono::milliseconds(1000));

	allpar_set AP;
	
	


	
	
	
	
	
	
	
	total.stop();
	total.print_elaps("s");
	
	return 0;
}
