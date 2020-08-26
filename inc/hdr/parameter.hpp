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
	double get_dbl(const std::string&, double);
};

//###########################################

class mode_cmd
{
	char **argv;
	int argc;
	
public:
	std::string cmd_line;
	std::string mode_str;

	std::string par1_str;
	std::string par2_str;
	std::string par3_str;
	
	double par1_start;
	double par1_stop;
	int par1_steps;
	
	std::string up_down;
	
	mode_cmd(int, char**);
};

//###########################################

class lpar_set
{
public:
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
public:
	par K{0.0, "K", "K"};
	par tau{0.0, "tau", "\\tau"};
	par wLP{0.0, "wLP", "\\omega_{\\text{LP}}"};

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
	par D{0.0, "D", "D"};
	par rea{0.0, "rea", "rea"};
	
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
	
	allpar_set();
	allpar_set(std::string optL, std::string optF, std::string optI);
	allpar_set(lpar_set LP_init, fpar_set FP_init, ipar_set IP_init);

	
	std::vector<par> collect();
	std::vector<par*> collect_ptr();
	void cout_pars(std::vector<par> collection);
	
	void check_cmd_line(int argc, char* argv[]);
	par* get_par_ptr(std::string);
	
	double larger_delay();
	double smaller_delay();
};

//###########################################

class icpar_set
{
public:
	par er_ic{0.0, "er_ic", "er_ic"};
	par ei_ic{0.0, "ei_ic", "ei_ic"};
	par g_ic{0.0, "g_ic"  , "g_ic"};
	par q_ic{0.0, "q_ic"  , "q_ic"};
	par j_ic{0.0, "j_ic"  , "j_ic"};
	
	par a1_ic{0.0, "a1_ic", "a1_ic"};
	par a2_ic{0.0, "a2_ic", "a2_ic"};
	par a3_ic{0.0, "a3_ic", "a3_ic"};
	par a4_ic{0.0, "a4_ic", "a4_ic"};
	par a5_ic{0.0, "a5_ic", "a5_ic"};
	
	par b1_ic{0.0, "b1_ic", "b1_ic"};
	par b2_ic{0.0, "b2_ic", "b2_ic"};
	par b3_ic{0.0, "b3_ic", "b3_ic"};
	par b4_ic{0.0, "b4_ic", "b4_ic"};
	par b5_ic{0.0, "b5_ic", "b5_ic"};
		
	
	icpar_set(std::string option);
	
	std::vector<par> collect();
	void cout_pars(std::vector<par> collection);
	
	void check_cmd_line(int argc, char* argv[]);
};


//###########################################

class lpar_dbl_set
{
public:
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
	double wLP;
  
	fpar_dbl_set(fpar_set FP);
	fpar_dbl_set(allpar_set AP);
};

//###########################################

class ipar_dbl_set
{
public:
	double int_time;
	double out_time;
	double dt;
	double sqrtdt;
	double D;
	double rea;

  
	ipar_dbl_set(ipar_set IP);
	ipar_dbl_set(allpar_set AP);
};

#endif
