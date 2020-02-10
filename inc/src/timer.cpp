#include <cmath>

#include <chrono>

#include <iostream>

#include <string>


#include "../hdr/timer.hpp"


using namespace std;
using namespace chrono;


timer::timer()
{
	this->start = system_clock::now();
}

void timer::stop()
{
	this->end = system_clock::now();
	this->elaps_duration = this->end - this->start;
	this->elaps = this->elaps_duration.count();
}

double timer::get_elaps()
{
	return this->elaps;
}

void timer::print_elaps()
{
	cout << "com time: " << this->elaps << "s" << endl;
}

void timer::print_elaps(string opt)
{
	if(opt == "s" or opt == "sec" or opt == "SEC")
	{
		cout << "com time: " << this->elaps << "s" << endl;
	}
	else if(opt == "ms")
	{
		cout << "com time: " << this->elaps * 1e3 << "ms" << endl;
	}
	else
	{
		print_elaps();
	}
}
