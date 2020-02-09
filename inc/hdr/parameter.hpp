#ifndef __PARAMETER_HPP_INCLUDED__
#define __PARAMETER_HPP_INCLUDED__

class par
{
public:
	double par_dbl;
	std::string par_str;
	std::string par_tex;

	par(double init_dbl, 
	    std::string init_str,
	    std::string init_tex
	   );
};

//###########################################

class par_cmd
{
	char **begin;
	char **end;
public:
	par_cmd(char **be,char **en);
	double get_dbl(const std::string &par_str, double dbl_default);
};

//###########################################

class lpar_set
{
public:
	par a{0.0, "a", "\\alpha"};
	par ag{0.0, "ag", "\\alpha_g"};
	par aq{0.0, "aq", "\\alpha_q"};
	par Jg{0.0, "Jg", "J_g"};
	par q0{0.0, "q0", "q_0"};
	par rs{0.0, "rs", "r_s"};
	par g{0.0, "g", "\\gamma"};
	par gg{0.0, "gg", "\\gamma_g"};
	par gq{0.0, "gq", "\\gamma_q"};
	par sqrtkap{0.0, "sqrtkap", "\\sqrt{\\kappa}"};
	par dw{0.0, "dw", "\\Delta\\Omega"};
	par T{0.0, "T", "T"};

	lpar_set(std::string option);

	std::vector<par> collect();
	void cout_pars(std::vector<par> collection);

};

//###########################################

class fpar_set
{
public:class cmd_line
{
	char **begin;
	char **end;
public:
	cmd_line(char **be,char **en);
	double get_dbl(const std::string &par_str, double dbl_default);
	bool get_bool(const std::string &par_str, bool bool_default);
	bool get_swp(const std::string &par_str, bool swp_default);

};
	par K{0.0, "K", "K"};
	par tau{0.0, "tau", "\\tau"};
	par omegaLP{0.0, "omegaLP", "\\omega_{\\text{LP}}"};

	fpar_set(std::string option);

	std::vector<par> collect();
	void cout_pars(std::vector<par> collection);
};

//###########################################

class ipar_set
{
public:
	par int_time{0.0, "int_time", "t_{\\text{int}}"};
	par out_time{0.0, "out_time", "t_{\\text{out}}"};
	par dt{0.0, "dt", "\\mathrm{d}t"};
	par sqrtdt{sqrt(0.0), "sqrtdt", "\\sqrt{\\mathrm{d}t}"};

	par t{0.0, "t", "t"};
	
	ipar_set(std::string option);
	
	std::vector<par> collect();
	void cout_pars(std::vector<par> collection);
};

//###########################################

class allpar_set
{
public:
	lpar_set LP{"def"};
	fpar_set FP{"def"};
	ipar_set IP{"def"};
	
	allpar_set(std::string optL, std::string optF, std::string optI);
	allpar_set(lpar_set LP_init, fpar_set FP_init, ipar_set IP_init);

	
	std::vector<par> collect();
	void cout_pars(std::vector<par> collection);
	void check_cmd_line(int argc, char* argv[]);
	double larger_delay();
};

//###########################################

class lpar_dbl_set
{
public:
	double a;
	double ag;
	double aq;
	double Jg;
	double q0;
	double rs;
	double g;
	double gg;
	double gq;
	double sqrtkap;
	double dw;
	double T;

	lpar_dbl_set(lpar_set LP);
	lpar_dbl_set(allpar_set AP);
};

//###########################################

class fpar_dbl_set
{
public:
	double K;
	double tau;
	double omegaLP;
  
	fpar_dbl_set(fpar_set FP);
	fpar_dbl_set(allpar_set AP);
};

//###########################################

class ipar_dbl_set
{
public:
	double int_time;
	double out_time;
	double t;
	double dt;
	double sqrtdt;
  
	ipar_dbl_set(ipar_set IP);
	ipar_dbl_set(allpar_set AP);
};

#endif
