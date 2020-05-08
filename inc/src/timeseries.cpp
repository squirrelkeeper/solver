#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>


#include "../hdr/timeseries.hpp"

using namespace std;

timeseries::timeseries()
{
}

timeseries::timeseries(allpar_set *init_AP)
{
	AP = init_AP;
	ipar_dbl_set ip(*AP);

	len_dbl = ip.out_time;	
	len = len_dbl/ip.dt;

	t.resize(len);
	I.resize(len);
	X.resize(len);
}

timeseries::timeseries(double init_len, allpar_set *init_AP)
{
	*AP = *init_AP;
	ipar_dbl_set ip(*AP);

	len_dbl = init_len;	
	len = len_dbl/ip.dt;
	
	t.resize(len);
	I.resize(len);
	X.resize(len);
}

void timeseries::write_file(string file_name)
{
	stringstream full_name;
	full_name.precision(7);
	
	
	srand (time(NULL));
	int batch_no = rand()%10000;

	full_name << "data/";
	full_name << file_name;
	full_name << batch_no;
	full_name << ".dat";


	ofstream data;
	data.open(full_name.str().c_str(),ios::trunc);


	data << "#t" << '\t';
	data << "ER" << '\t';
	data << "EI" << '\t';
	data << "G" << '\t';
	data << "Q" << '\t';
	data << "J" << '\t';
	data << "I" << '\t';
	data << endl;
	
	
	
	for(long i = 0; i < len; i++)
	{
		data << setprecision(15);

		data << t[i] << '\t';
		data << X[i].ER << '\t';
		data << X[i].EI << '\t';
		data << X[i].G << '\t';
		data << X[i].Q << '\t';
		data << X[i].J << '\t';
		data << I[i] << '\t';

		data << endl;
	}
	
	data.close();
}
