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

varC varC::operator+(const varC& a)
{
	varC sum;
	
	sum.E = a.E + this->E;
	sum.G = a.G + this->G;
	sum.Q = a.Q + this->Q;
	sum.J = a.J + this->J;
	
	return sum;
}

varC varC::operator*(const double& a)
{
	varC prod;
	
	prod.E = a * this->E;
	prod.G = a * this->G;
	prod.Q = a * this->Q;
	prod.J = a * this->J;
	
	return prod;
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

var var::operator+(const var& a)
{
	var sum;
	
	sum.ER = a.ER + this->ER;
	sum.EI = a.EI + this->EI;
	sum.G = a.G + this->G;
	sum.Q = a.Q + this->Q;
	sum.J = a.J + this->J;

	return sum;
}

var var::operator*(const double& a)
{
	var prod;
	
	prod.ER = a * this->ER;
	prod.EI = a * this->EI;
	prod.G = a * this->G;
	prod.Q = a * this->Q;
	prod.J = a * this->J;

	return prod;
}
/* DOESNT WORK
var var::operator=(const varC &C)
{
	var R;
	R.ER = C.E.real();
	R.EI = C.E.imag();
	R.G = C.G;
	R.Q = C.Q;
	R.J = R.J;
	
	return R;
}
*/
