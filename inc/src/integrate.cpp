#include <algorithm>
#include <complex>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../hdr/integrate.hpp"


using namespace std;

integrator::integrator(allpar_set *AP)
{
	lpar_dbl_set lp(*AP);
	fpar_dbl_set fp(*AP);
	ipar_dbl_set ip(*AP);

	it = (long)(ip.int_time / ip.dt);

	dim1 = (long)floor(lp.T / ip.dt );
	dim2 = (long)floor(fp.tau / ip.dt + 1.0);

	pos0 = dim1 - 1;
	pos1 = pos0 - dim1;
	pos2 = 0;
	
	t = 0.0;
	
	C.resize(dim2);
	

	
	
	
	cout << "it" << '\t' << "dim1" << '\t' << "dim2" << endl;
	cout << it << '\t' << dim1 << '\t' << dim2 << endl;
}
