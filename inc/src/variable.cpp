#include <complex>

#include "../hdr/variable.hpp"


using namespace std;

varC::varC()
{
	E = 0.0;
	G = 0.0;
	Q = 0.0;
	J = 0.0;
}

varC::varC(complex<double> iE, double iG, double iQ, double iJ)
{
	E = iE;
	G = iG;
	Q = iQ;
	J = iJ;
}

//###########################################

var::var()
{
	ER = 0.0;
	EI = 0.0;
	G = 0.0;
	Q = 0.0;
	J = 0.0;
}

var::var(double iER, double iEI, double iG, double iQ, double iJ)
{
	ER = iER;
	EI = iEI;
	G = iG;
	Q = iQ;
	J = iJ;
}
