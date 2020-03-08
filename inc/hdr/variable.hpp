#ifndef __VARIABLE_HPP_INCLUDED__
#define __VARIABLE_HPP_INCLUDED__

class varC
{
public:
	std::complex<double> E;
	double G;
	double Q;
	double J;

	varC();
	varC(std::complex<double> iE, double iG, double iQ, double iJ);
	
	varC operator+(const varC&);
	varC operator*(const double&);
};

//###########################################

class var
{
public:
	double ER;
	double EI;
	double G;
	double Q;
	double J;

	var();
	var(double iER, double iEI, double iG, double iQ, double iJ);
	
	var operator+(const var&);
	var operator*(const double&);
};

#endif
