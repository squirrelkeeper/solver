#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>


#include "../hdr/parameter.hpp"


using namespace std;

par::par(double init_dbl, 
		 string init_str,
		 string init_tex
		)
{
	par_dbl = init_dbl;
	par_str = init_str;
	par_tex = init_tex;
}

//###########################################

par_cmd::par_cmd(char **be, char **en)
	: begin(be), end(en)
{}

double par_cmd::get_dbl(const string &par_str, double dbl_default)
{
	string new_str = "-" + par_str;
	char **itr = find(begin, end, new_str);
	if (itr != end && ++itr != end)
	{
		char *b;
		double d;
		d = strtod(*itr, &b);
		if (0 == d && *itr == b)
		{
			cout << "Input of "<< par_str << " was wrong. ";
			cout << "Setting to default: "<< dbl_default << endl;
			return dbl_default;
		}
		return d;
	}
	return dbl_default;
}

//###########################################


lpar_set::lpar_set(string option)
{
	if(option == "JAU15")
	{
		this->Jg.par_dbl = 3.0;
		this->q0.par_dbl = 4.0;
		this->rs.par_dbl = 25.0;
		this->g.par_dbl = 66.5;
		this->gg.par_dbl = 0.025;
		this->gq.par_dbl = 1.875;
		this->sqrtkap.par_dbl = sqrt(0.1);
		this->T.par_dbl = 1.0;
	}
	else if(option == "2MM")
	{
		this->Jg.par_dbl = 3.0;
		this->q0.par_dbl = 4.0;
		this->rs.par_dbl = 25.0;
		this->g.par_dbl = 66.5;
		this->gg.par_dbl = 0.025;
		this->gq.par_dbl = 1.875;
		this->sqrtkap.par_dbl = sqrt(0.1);
		this->T.par_dbl = 2.0;
	}
	else if(option == "def")
	{
	}
	else
	{
		cout << "err000" << endl;
	}
}

std::vector<par> lpar_set::collect()
{
	vector<par> collection;
	
	collection.push_back(this->ag);
	collection.push_back(this->aq);
	collection.push_back(this->Jg);
	collection.push_back(this->q0);
	collection.push_back(this->rs);
	collection.push_back(this->g);
	collection.push_back(this->gg);
	collection.push_back(this->gq);
	collection.push_back(this->sqrtkap);
	collection.push_back(this->dw);
	collection.push_back(this->T);
	
	return collection;
}

void lpar_set::cout_pars(vector<par> collection)
{
	for(unsigned int i = 0; i < collection.size(); i++)
	{
		cout << collection[i].par_str;
		cout << " = ";
		cout << collection[i].par_dbl;
		cout << endl;
	}
}

//###########################################

fpar_set::fpar_set(string option)
{
	if(option == "TAU1")
	{
		this->K.par_dbl = 1.0;
		this->tau.par_dbl = 2.0;
		this->wLP.par_dbl = 2.0 * M_PI / 1.0;
	}
	else if(option == "TAU2")
	{
		this->K.par_dbl = 1.0;
		this->tau.par_dbl = 3.0;
		this->wLP.par_dbl = 2.0 * M_PI / 3.0;
	}
	else if(option == "TAU5")
	{
		this->K.par_dbl = 1.0;
		this->tau.par_dbl = 6.0;
		this->wLP.par_dbl = 2.0 * M_PI / 6.0;
	}
	else if(option == "def")
	{
	}
	else
	{
		cout << "err001" << endl;
	}
}

std::vector<par> fpar_set::collect()
{
		vector<par> collection;
		
		collection.push_back(this->K);
		collection.push_back(this->tau);
		collection.push_back(this->wLP);

		return collection;
}

void fpar_set::cout_pars(vector<par> collection)
{
	for(unsigned int i = 0; i < collection.size(); i++)
	{
		cout << collection[i].par_str;
		cout << " = ";
		cout << collection[i].par_dbl;
		cout << endl;
	}
}

//###########################################

ipar_set::ipar_set(string option)
{
	if(option == "quick")
	{
		this->int_time.par_dbl = 300.0;
		this->out_time.par_dbl = 100.0;
		this->dt.par_dbl = 1e-4;
		this->sqrtdt.par_dbl = sqrt(1e-4);
		this->D.par_dbl = 0.0;
	}
	else if(option == "def")
	{
		this->int_time.par_dbl = 100.0;
		this->out_time.par_dbl = 50.0;
		this->dt.par_dbl = 1e-3;
		this->sqrtdt.par_dbl = sqrt(1e-3);
		this->D.par_dbl = 0.0;
	}
	else if(option == "noise")
	{
		this->int_time.par_dbl = 300.0;
		this->out_time.par_dbl = 100.0;
		this->dt.par_dbl = 1e-4;
		this->sqrtdt.par_dbl = sqrt(1e-4);
		this->D.par_dbl = 0.2;
		this->rea.par_dbl = 1.0;
	}
	else
	{
		cout << "err002" << endl;
	}
}

vector<par> ipar_set::collect()
{
	vector<par> collection;
	
	collection.push_back(this->int_time);
	collection.push_back(this->out_time);
	collection.push_back(this->dt);
	collection.push_back(this->sqrtdt);
	collection.push_back(this->D);

	
	return collection;
}

void ipar_set::cout_pars(vector<par> collection)
{
	for(unsigned int i = 0; i < collection.size(); i++)
	{
		cout << collection[i].par_str;
		cout << " = ";
		cout << collection[i].par_dbl;
		cout << endl;
	}
}

//###########################################

allpar_set::allpar_set()
{

}

allpar_set::allpar_set(string optL, string optF, string optI)
{
	lpar_set LP_init{optL};
	fpar_set FP_init{optF};
	ipar_set IP_init{optI};

	this->LP = LP_init;
	this->FP = FP_init;
	this->IP = IP_init;
}

allpar_set::allpar_set(lpar_set LP_init, fpar_set FP_init, ipar_set IP_init)
{
	this->LP = LP_init;
	this->FP = FP_init;
	this->IP = IP_init;
}

vector<par> allpar_set::collect()
{
	vector<par> vecLP;
	vector<par> vecFP;
	vector<par> vecIP;
	
	vecLP = this->LP.lpar_set::collect();
	vecFP = this->FP.fpar_set::collect();
	vecIP = this->IP.ipar_set::collect();
	
	vecLP.insert(vecLP.end(), vecFP.begin(), vecFP.end());
	vecLP.insert(vecLP.end(), vecIP.begin(), vecIP.end());

	return vecLP;
}

void allpar_set::cout_pars(vector<par> collection)
{
	for(unsigned int i = 0; i < collection.size(); i++)
	{
		cout << collection[i].par_str;
		cout << " = ";
		cout << collection[i].par_dbl;
		cout << endl;
	}
}

void allpar_set::check_cmd_line(int argc, char* argv[])
{
	par_cmd cmd(argv, argv+argc);
	
	this->LP.ag.par_dbl = cmd.get_dbl(this->LP.ag.par_str, this->LP.ag.par_dbl);
	this->LP.aq.par_dbl = cmd.get_dbl(this->LP.aq.par_str, this->LP.aq.par_dbl);
	this->LP.Jg.par_dbl = cmd.get_dbl(this->LP.Jg.par_str, this->LP.Jg.par_dbl);
	this->LP.q0.par_dbl = cmd.get_dbl(this->LP.q0.par_str, this->LP.q0.par_dbl);
	this->LP.rs.par_dbl = cmd.get_dbl(this->LP.rs.par_str, this->LP.rs.par_dbl);
	this->LP.g.par_dbl = cmd.get_dbl(this->LP.g.par_str, this->LP.g.par_dbl);
	this->LP.gg.par_dbl = cmd.get_dbl(this->LP.gg.par_str, this->LP.gg.par_dbl);
	this->LP.gq.par_dbl = cmd.get_dbl(this->LP.gq.par_str, this->LP.gq.par_dbl);
	this->LP.sqrtkap.par_dbl = cmd.get_dbl(this->LP.sqrtkap.par_str, this->LP.sqrtkap.par_dbl);
	this->LP.dw.par_dbl = cmd.get_dbl(this->LP.dw.par_str, this->LP.dw.par_dbl);
	this->LP.T.par_dbl = cmd.get_dbl(this->LP.T.par_str, this->LP.T.par_dbl);

	this->FP.K.par_dbl = cmd.get_dbl(this->FP.K.par_str, this->FP.K.par_dbl);
	this->FP.tau.par_dbl = cmd.get_dbl(this->FP.tau.par_str, this->FP.tau.par_dbl);
	this->FP.wLP.par_dbl = cmd.get_dbl(this->FP.wLP.par_str, this->FP.wLP.par_dbl);
	
	this->IP.int_time.par_dbl = cmd.get_dbl(this->IP.int_time.par_str, this->IP.int_time.par_dbl);
	this->IP.out_time.par_dbl = cmd.get_dbl(this->IP.out_time.par_str, this->IP.out_time.par_dbl);
	this->IP.dt.par_dbl = cmd.get_dbl(this->IP.dt.par_str, this->IP.dt.par_dbl);
	this->IP.sqrtdt.par_dbl = sqrt(cmd.get_dbl(this->IP.dt.par_str, sqrt(this->IP.dt.par_dbl)));
	this->IP.D.par_dbl = cmd.get_dbl(this->IP.D.par_str, this->IP.D.par_dbl);
	this->IP.rea.par_dbl = cmd.get_dbl(this->IP.rea.par_str, this->IP.rea.par_dbl);
}

double allpar_set::larger_delay()
{
	double T1 = this->LP.T.par_dbl;
	double T2 = this->FP.tau.par_dbl;
	double max_tau = T1 < T2 ? T2 : T1; 
	return max_tau;
}

double allpar_set::smaller_delay()
{
	double T1 = this->LP.T.par_dbl;
	double T2 = this->FP.tau.par_dbl;
	double min_tau = T1 > T2 ? T2 : T1; 
	return min_tau;
}

//###########################################

lpar_dbl_set::lpar_dbl_set(lpar_set LP)
{
	this->aq = LP.aq.par_dbl;
	this->ag = LP.ag.par_dbl;
	this->Jg = LP.Jg.par_dbl;
	this->q0 = LP.q0.par_dbl;
	this->rs = LP.rs.par_dbl;
	this->g = LP.g.par_dbl;
	this->gg = LP.gg.par_dbl;
	this->gq = LP.gq.par_dbl;
	this->sqrtkap = LP.sqrtkap.par_dbl;
	this->dw = LP.dw.par_dbl;
	this->T = LP.T.par_dbl;
}

lpar_dbl_set::lpar_dbl_set(allpar_set AP)
{
	this->aq = AP.LP.aq.par_dbl;
	this->ag = AP.LP.ag.par_dbl;
	this->Jg = AP.LP.Jg.par_dbl;
	this->q0 = AP.LP.q0.par_dbl;
	this->rs = AP.LP.rs.par_dbl;
	this->g = AP.LP.g.par_dbl;
	this->gg = AP.LP.gg.par_dbl;
	this->gq = AP.LP.gq.par_dbl;
	this->sqrtkap = AP.LP.sqrtkap.par_dbl;
	this->dw = AP.LP.dw.par_dbl;
	this->T = AP.LP.T.par_dbl;
}

//###########################################

fpar_dbl_set::fpar_dbl_set(fpar_set FP)
{
	this->K = FP.K.par_dbl;
	this->tau = FP.tau.par_dbl;
	this->wLP = FP.wLP.par_dbl;
}

fpar_dbl_set::fpar_dbl_set(allpar_set AP)
{
	this->K = AP.FP.K.par_dbl;
	this->tau = AP.FP.tau.par_dbl;
	this->wLP = AP.FP.wLP.par_dbl;
}

//###########################################

ipar_dbl_set::ipar_dbl_set(ipar_set IP)
{
	this->int_time = IP.int_time.par_dbl;
	this->out_time = IP.out_time.par_dbl;
	this->dt = IP.dt.par_dbl;
	this->sqrtdt = IP.sqrtdt.par_dbl;
	this->D = IP.D.par_dbl;
	this->rea = IP.rea.par_dbl;
}

ipar_dbl_set::ipar_dbl_set(allpar_set AP)
{
	this->int_time = AP.IP.int_time.par_dbl;
	this->out_time = AP.IP.out_time.par_dbl;	
	this->dt = AP.IP.dt.par_dbl;
	this->sqrtdt = AP.IP.sqrtdt.par_dbl;
	this->D = AP.IP.D.par_dbl;
	this->rea = AP.IP.rea.par_dbl;
}
