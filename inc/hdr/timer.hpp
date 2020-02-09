#ifndef __TIMER_HPP_INCLUDED__
#define __TIMER_HPP_INCLUDED__

class timer
{
	std::chrono::time_point<std::chrono::system_clock> start;
	std::chrono::time_point<std::chrono::system_clock> end;
	std::chrono::duration<double> elaps_duration;
	double elaps;

public:	
	timer();
	
	void stop();
	
	double get_elaps();
	double get_elaps(std::string);

	void print_elaps();
	void print_elaps(std::string);
};

#endif
