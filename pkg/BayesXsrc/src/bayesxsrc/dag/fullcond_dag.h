/* BayesX - Software for Bayesian Inference in
Structured Additive Regression Models.
Copyright (C) 2011  Christiane Belitz, Andreas Brezger,
Thomas Kneib, Stefan Lang, Nikolaus Umlauf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */



#include"../export_type.h"

#if !defined (FULLCOND_DAG_INCLUDED)

#define FULLCOND_DAG_INCLUDED

#include"mcmc.h"
#include"distribution.h"
#include"fullcond.h"
#include"adjacency.h"

//#include <FULLCOND_DAG_D_IA.h>


namespace MCMC
{

class FULLCOND_dag_ia;



class __EXPORT_TYPE FULLCOND_dag : public FULLCOND
  {

  protected:

	  datamatrix beta_help;  // matrix which is corresponding to beta, but without zeros

	  datamatrix y;
	  datamatrix lin;
	  datamatrix lin_prop;

	  Matrix<unsigned> adcol;		// responding column of adja
	  datamatrix x;					// matrix of covariates
	  datamatrix xx;				// x'x
	  datamatrix Sigma;				// covariance matrix of the full conditional distribution of beta

	  // matrices for the proposed paramters in a death or birth step
	  // (avoids to compute always a new matrix in eaxh rj step
	  datamatrix b_new_b;
	  datamatrix b_new_d;
	  datamatrix x_new_b;
	  datamatrix x_new_d;
	  datamatrix xx_new_b;
	  datamatrix xx_new_d;

	  datamatrix beta_mean;			// mean of the full conditional distribution of beta

	  double sigma_i;
	  double alpha;					//scaling paramater

	  double sigma_prop;			// variance of the proposal distribution

	  double SQT_x;		//sum of squares in total:  sum_1^n (x-µ)^2
	  double SQT_b;		//sum of squares in total:  sum_1^n (b_i-0)^2

	  double SQT_x_n;		// the new SQT_x for the proposed model
	  double SQT_b_n;		// the new SQT_b for the proposed model

	  double a_invg;	// 1. parameter of the inv-gamma distribution of sigma_i
	  double b_invg;	// 2. parameter of the inv-gamma distribution of sigma_i

	  unsigned int ncoef;		// ncoef = ncoef_m+ncoef_ia
	  unsigned int ncoef_m;		// number of coefficients of main effects;
								// = number of parents + 1
	  unsigned int ncoef_ia;	// number of interactions (=0, if continous)

	  unsigned int nvar; //number of variables
	  unsigned int nobs; //number of observations
	  unsigned int self; //number of regarded regression model

	  bool print_dags; //if true, all estimations of the dags are listend in the outpu
	  ST::string priori_sigma;
	  ST::string priori_beta;

	  char var_type;		//  'c' for continuous, 'd' for discrete

  // add here additional variables needed for the new class
  // add here additional functions/methods needed for the inherited class
  // (e.g. a method for centering beta)

  public:

  // DEFAULT CONSTRUCTOR:

  FULLCOND_dag(void) : FULLCOND()
    {
    }

  // CONSTRUCTOR


  // CONSTRUCTOR
  // o    : pointer to MCMCoptions object
  // dp   : pointer to DISTRIBUTION object
  // t    : title of the full conditional (for example "fixed effects")
  // rows : number of rows of the beta matrix (i.e. number of parameters)
  // cols : number of columns of the beta matrix
  //        (i.e. number of categories of the response variable)
  // fp   : file path for storing sampled parameters






  // CONSTRUCTOR 1
  FULLCOND_dag (double s_i, unsigned int num,
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp);

  // CONSTRUCTOR 2
  FULLCOND_dag (double value_a, double value_b, ST::string prio_sig, bool dags_all,
               const datamatrix & res, double s_i, unsigned int num,
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp);


  // COPY CONSTRUCTOR

  FULLCOND_dag(const FULLCOND_dag & fc) : FULLCOND(FULLCOND(fc))
  {
	  beta_help=fc.beta_help;
	  y = fc.y;
	  lin = fc.lin;
	  lin_prop = fc.lin_prop;
	  adcol= fc.adcol;
	  x = fc.x;
	  xx = fc.xx;
	  Sigma = fc.Sigma;
	  b_new_b = fc.b_new_b;
	  b_new_d = fc.b_new_d;
	  x_new_b = fc.x_new_b;
	  x_new_d = fc.x_new_d;
	  xx_new_b = fc.xx_new_b;
	  xx_new_d = fc.xx_new_d;
	  beta_mean = fc.beta_mean;
	  sigma_i = fc.sigma_i;
	  alpha= fc.alpha;
	  sigma_prop = fc.sigma_prop;
	  SQT_x = fc.SQT_x;
	  SQT_b = fc.SQT_b;
	  SQT_x_n = fc.SQT_x_n;
	  SQT_b_n = fc.SQT_b_n;
	  ncoef = fc.ncoef;
	  ncoef_m = fc.ncoef_m;
	  ncoef_ia = fc.ncoef_ia;
	  nvar = fc.nvar;
	  nobs = fc.nobs;
	  self = fc.self;
	  a_invg = fc.a_invg;
	  b_invg = fc.b_invg;
	  print_dags = fc.print_dags;
	  priori_sigma = fc.priori_sigma;
	  priori_beta = fc.priori_beta;
	  var_type = fc.var_type;

	  // assign here ONLY additional variables of the inherited class
    // e.g. y = fc.y
    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_dag & operator=(const FULLCOND_dag & fc)
  {
	  if (this==&fc)
		  return *this;

	  beta_help=fc.beta_help;
	  y = fc.y;
	  lin = fc.lin;
	  lin_prop = fc.lin_prop;
	  adcol= fc.adcol;
	  x= fc.x;
	  xx= fc.xx;
	  Sigma = fc.Sigma;
	  b_new_b = fc.b_new_b;
	  b_new_d = fc.b_new_d;
	  x_new_b = fc.x_new_b;
	  x_new_d = fc.x_new_d;
	  xx_new_b = fc.xx_new_b;
	  xx_new_d = fc.xx_new_d;
	  beta_mean= fc.beta_mean;
	  sigma_i= fc.sigma_i;
	  alpha= fc.alpha;
	  sigma_prop = fc.sigma_prop;
	  SQT_x = fc.SQT_x;
	  SQT_b = fc.SQT_b;
	  SQT_x_n = fc.SQT_x_n;
	  SQT_b_n = fc.SQT_b_n;
	  ncoef= fc.ncoef;
      ncoef_m = fc.ncoef_m;
	  ncoef_ia = fc.ncoef_ia;
	  nvar= fc.nvar;
	  nobs = fc.nobs;
	  self = fc.self;
	  a_invg = fc.a_invg;
	  b_invg = fc.b_invg;
	  print_dags = fc.print_dags;
	  priori_sigma = fc.priori_sigma;
	  priori_beta = fc.priori_beta;
	  var_type = fc.var_type;


	  // assign here ONLY additional variables of the inherited class
    // e.g. y = fc.y

    return *this;

    }


  // OVERLOADED ASSIGNMENT OPERATOR
  //FULLCOND_dag   operator= (const FULLCOND_dag_d_ia & dag_ia);






  // DESTRUCTOR
  ~FULLCOND_dag() {}

  virtual ST::string get_family(void);



	// FUNCTION: get_current_ia
	// TASK: returns current_ia
	virtual double get_ln_prop_beta (void)
	{
		return 0;
	}


  // FUNCTION get_nobs
  // TASK: gives nobs
  unsigned get_nobs(void)
  {
	  return nobs;
  }



  // FUNCTION update
  // TASK: updates parameters (i.e. matrix beta and beta_help)
  void update(void);


  // FUNCTION initialize
  // TASK: initializes x and xx  and y for pred_mod[i] (regression model i)
  void initialize (const adja & zeta, unsigned int j);



  // FUNCTION initialize_ia
  // TASK: initializes
  virtual void initialize_ia (const adja & zeta, unsigned int j)
  {
  }

  // FUNCTION: write_ia_to_x
  // TASK: adds the interaction to x when x is initialized
  // is called in FULLCOND_dag::initialize
	virtual void write_ia_to_x(void)
	{
	}





  // FUNCTION create_matrices();
  // TASK: creates matrices which are needed for proposals of an rj-step
  virtual void create_matrices (void);


  // FUNCTION create_matrices();
  // TASK: creates matrices which are needed for proposals of an rj-step
  void create_matrices (ST::string step, unsigned ncoef_new);


  // FUNCTION: ia_of_i
  // TASK: counts the number of interactions containing variable i
  virtual unsigned ia_of_i(unsigned)
  {
	  return 0;
  }

   // FUNCTION: ia_of_i
  // TASK:adds interactions containing variable i to v
  virtual void ia_of_i( unsigned i, vector <vector <unsigned > > & v)
  {
  }


  // FUNCTION: ia_of_i
  // TASK:adds interactions containing variable i to v
  // when i IS NOT already a main effect
  virtual void new_ia_of_i( unsigned i, vector <vector <unsigned > > & v)
  {
  }



  // TASK: returns the number of interactions of the main effect i
  // if an interaction between all pairs of main effects is assumed
  virtual unsigned ia_of_i( void)
  {
	  return 0;
  }






  // FUNCTION: get_pos
	// TASK: gives back the position of main effect i and
	//		 the corresponding ia in the regression model
	//		 starting with 0 for the intercept
	virtual void get_pos(unsigned i, vector <unsigned> & pos)
	{
	}



	// FUNCTION: get_current_ia
	// TASK: returns current_ia
	virtual vector <vector <unsigned > > & get_current_ia (void)
	{
		vector <vector <unsigned > > test;
		vector <vector <unsigned > >  & reference = test;
		return reference;
	}

	// FUNCTION: get_current_ia
	// TASK: returns k-th element of current_ia
	virtual vector<unsigned> get_current_ia (unsigned k)
	{
		vector<unsigned> test;
		return test;
	}

	// FUNCTION: change_occur()
	// TARGET: adds or delets ia in occur after birth or death step
	virtual void change_occur(char step, vector <vector <unsigned > >  ia_new)
	{
	}

	// FUNCTION: change_current()
	// TARGET: adds or delets ia.term to/from current_ia
	virtual void change_current(char step, vector <vector <unsigned > > term)
	{
	}

	virtual void make_new_b (unsigned i, vector <vector <unsigned> > ia_new,
							datamatrix & beta_new, datamatrix & xx_new,
						datamatrix & b_new, datamatrix & x_new)
	{
	}




	// FUNCTION: make_new_b_int
	// TASK: computes the new values for a birth-step when one new main effect
	// and ALL corresponding interactions are added
	// is called from birth_step in rj_int
	virtual void make_new_b_int (ST::string step, unsigned i, vector <vector <unsigned> > ia_new,
				datamatrix & beta_new, datamatrix & xx_new, datamatrix & b_new,
				datamatrix & x_new)
	{
	}


   // FUNCTION: make_new_d_int
   // TASK: computes the new values for a death-step
   virtual void make_new_d_int (ST::string step, unsigned i, unsigned j, unsigned ia_del,
						datamatrix & beta_old, vector <vector <unsigned > > & current_ia_n,
						datamatrix & xx_new, datamatrix & b_new, datamatrix & x_new)
   {

   }




  // FUNCTION calc_kq_est
  // TASK: sets beta_n equal to the kq-estimator for given x, xx
  // where beta_n denotes ALL  regression coefficions of the models (old ones and added)
  void calc_kq_est (datamatrix & beta_n, const datamatrix & x, const datamatrix & xx);


  // FUNCTION calc_kq_est
  // TASK: sets beta_n equal to the kq-estimator for given x_ia, xx_ia
  // where beta_n denotes ONLY the regression coefficions of the added covariables
  void calc_kq_est (datamatrix & beta_n, const datamatrix & x_ia,
							const datamatrix & xx_ia, const datamatrix & y_ia);



  // FUNCTION calc_beta_mean
  // TASK: calculates beta_mean
  void calc_beta_mean (void);



  // FUNCTION calc_Sigma
  // TASK: calculates Sigma
  void calc_Sigma (void);



   // FUNCTION get_y
  // TASK: returns y
  datamatrix  & get_y (void)
  {
	  datamatrix & reference = y;
	  return reference;
  }

   // FUNCTION get_y(i)
  // TASK: returns y(i,0)
  double get_y (unsigned int i)
  {
	  return y(i,0);
  }

   // FUNCTION put_y(value,i)
  // TASK:  sets y(i,0) equal to value (needed after updating utilities)
  void put_y (double value, unsigned int i)
  {
	  y(i,0) = value;
  }




  // FUNCTION get_sigma_i
  // TASK: returns sigma_i, the conditioned variance of the regression model
  double get_sigma_i (void)
  {
	  return sigma_i;
  }



  // FUNCTION get_x
  // TASK: gives x
  datamatrix get_x (void)
  {
	  datamatrix & reference = x;
	  return reference;
  }

  // FUNCTION getV_x
  // TASK: gives x.getV()
  double *  getV_x (void)
  {
	  return x.getV();
  }



  // FUNCTION get_x
  // TASK: gives i-th row of x
  datamatrix get_x (unsigned i)
  {
	  return x.getRow(i);
  }

  // FUNCTION get_x(i,j)
  // TASK: gives ij-th coefficient of xx
  double get_x (unsigned int i, unsigned int j)
  {
	  assert(i<x.rows() && j<x.cols());
	  return x(i,j);
  }


  // FUNCTION get_xx(i,j)
  // TASK: gives xx(i,j)
  double get_xx (unsigned int i, unsigned int j)
  {
	  assert(i<xx.rows() && j<xx.cols());
	  return xx(i,j);
  }



  // FUNCTION get_xx
  // TASK: gives xx
  datamatrix & get_xx (void)
  {
	  datamatrix & reference = xx;
	  return reference;
  }

  // FUNCTION getV_xx
  // TASK: gives pointer to xx
  double *  getV_xx (void)
  {
	  return xx.getV();
  }



  // FUNCTION getV_beta,
  // TASK: gives pointer to beta_help
  double *  getV_beta_help(void)
  {
	  return beta_help.getV();
  }


  // FUNCTION get_beta_help,
  // TASK: gives i-th component beta
  double get_beta(unsigned  i)
  {
	  return beta(i,0);
  }



  // FUNCTION get_beta_help,
  // TASK: gives beta_help
  datamatrix & get_beta_help (void)
  {
	  datamatrix & reference = beta_help;
	  return reference;
  }




  // FUNCTION get_beta_help(i,j)
  // TASK: gives beta(i,j)
  double get_beta_help (unsigned int i, unsigned int j)
  {
	  return beta_help(i,j);
  }



  // FUNCTION get_beta_mean
  // TASK: gives beta_mean
  datamatrix & get_beta_mean (void)
  {
	  datamatrix & reference = beta_mean;
	  return reference;
  }

  // FUNCTION get_beta_mean
  // TASK: gives i-th element of beta_mean
  double get_beta_mean (unsigned int i)
  {
	  return beta_mean(i,0);
  }



  // FUNCTION get_ncoef()
  // TASK: gives ncoef()
  double get_ncoef (void)
  {
	  return ncoef;
  }


  // FUNCTION get_ncoef()
  // TASK: gives ncoef_ia
  unsigned get_ncoef_ia (void)
  {
	  return ncoef_ia;
  }


  // FUNCTION get_ncoef_m()
  // TASK: gives ncoef_m()
  unsigned get_ncoef_m (void)
  {
	  return ncoef_m;
  }


  // FUNCTION get_b_new_b
  // TASK: returns b_new_b as reference
  datamatrix get_b_new_b (void)
  {
	  return  b_new_b;
  }


  // FUNCTION get_x_new_b
  // TASK: returns x_new_b as reference
  datamatrix   get_x_new_b (void)
  {
	  return  x_new_b;
  }


  // FUNCTION get_xx_new_b
  // TASK: returns xx_new_b as reference
  datamatrix   get_xx_new_b (void)
  {
	  return xx_new_b;
  }


  // FUNCTION get_b_new_d
  // TASK: returns b_new_d as reference
  datamatrix  & get_b_new_d (void)
  {
	  datamatrix & reference = b_new_d;
	  return reference;
  }


  // FUNCTION get_x_new_d
  // TASK: returns x_new_d as reference
  datamatrix  & get_x_new_d (void)
  {
	  datamatrix & reference = x_new_d;
	  return reference;
  }


  // FUNCTION get_xx_new_d
  // TASK: returns xx_new_d as reference
  datamatrix  & get_xx_new_d (void)
  {
	  datamatrix & reference = xx_new_d;
	  return reference;
  }




  // FUNCTION: num_pa
  // TASK: returns the number of parents of variable i
  unsigned int  num_pa(void) const
  {
	  unsigned int num =0;
	  for(unsigned int k=0; k<nvar; k++)
	  {
		  num = num + adcol(k,0);
	  }

    return num;
  }


  // FUNCTION: log_p_x
  // TASK: returns the log-likelihood of the x
  double log_p_x();



  // FUNCTION: log_p_x
  // TASK: returns the log-likelihood of the x with the proposed beta_help
  double log_p_x(const datamatrix & b_new, const datamatrix & x_new);


  // FUNCTION: log_p_x
  // TASK: returns the log-likelihood of the x with the proposed beta_help and sig_i_new
  double log_p_x(const datamatrix & b_new,const datamatrix & x_new, double sig_i_new);


  // FUNCTION: log_p_b1
  // TASK: returns the log-probability of beta_help with variance regarded as constant
  double log_p_b1(void);


  // FUNCTION: log_p_b1
  // TASK: returns the log-probability of the proposed b_n with variance regarded as constant
  double log_p_b1(const datamatrix & b_n);


  // FUNCTION: log_p_b2
  // TASK: returns the log-probability of beta_help with variance NOT regarded as constant
  double log_p_b2(void);




  // FUNCTION: log_p_b
  // TASK: returns the log-probability of the proposed beta_help
  double log_p_b(const datamatrix & b_new);


   // FUNCTION: log_p_b
  // TASK: returns the log-probability of the proposed beta_help and the new sigma
  double log_p_b(const datamatrix & b_new, double sigma);


  // FUNCTION: log_p_b
  // TASK: returns the log-probability of the proposed beta_help
//  double FULLCOND_dag::log_p_b(const datamatrix & b_new, double sig_i_new);




	// FUNCTION: p_prop()
	// TARGET: returns the density of the proposal u, which is normaldistributed
	double p_prop(double prop);


		// FUNCTION: p_prop()
	// TARGET: returns the log-density of the proposal u, which is normaldistributed
	double p_prop(double prop, double mu, double sigma);


	// FUNCTION: p_prop()
	// TARGET: returns the log-density of the proposal u, which is normaldistributed
	double p_prop(const datamatrix & prop, const datamatrix & mu, double sigma);


	// FUNCTION: p_prop()
	// TARGET: returns the log-density of the proposal u, which is normaldistributed
	double p_prop(const datamatrix & prop, const datamatrix & mu, const datamatrix & concentration);


	// FUNCTION: b_distr()
	// TARGET:
	double b_distr(void);


	double p_prop(void) ;



  // FUNCTION: change
  // TASK: changes the values after step has been accepted
  virtual void change(unsigned i, const datamatrix & beta_help_new, const datamatrix & x_new,
								const datamatrix & xx_new, unsigned int ncoef_new);

  // FUNCTION: change_adcol
  // TASK: changes adcol after step has been accepted
  void change_adcol(unsigned int i, unsigned int value)
  {
	   adcol(i,0) = value;
  }

  Matrix<unsigned> get_adcol(void)
  {
	   return adcol;
  }


  // FUNCTION: new_comp_xx (m,n,i)
  // TASK: calculates the mn-th element of the new column or row of x'x
  // when i is added as a new covariate (i is the t-th coefficient of j)
	double new_comp_xx(unsigned int m, unsigned int n, unsigned int i,
						unsigned int t);



	// FUNCTION: calc_lin()
	// TASK: calculates linear predictor
	void calc_lin(void);


	void calc_lin_o(void);



	// FUNCTION: get_lin()
	// TASK: get linear predictor
	 datamatrix get_lin(void)
  {
	  return lin;
  }

	// FUNCTION: calc_lin(x_prop,  b_prop)
	// TASK: calculates the proposed linear predictor
	void calc_lin_prop(const datamatrix & x_prop, const datamatrix & b_prop);




	// FUNCTION: calc_xx()
	// TASK: calculates xx
	void calc_xx(void);


	// FUNCTION: calc_SQT_x()
	// TASK: calculates the SQT_x
	double calc_SQT_x(void);


	// FUNCTION: calc_SQT_x(X_n, b_n)
	// TASK: calculates the SQT_x for the proposed X_n and b_n
	double calc_SQT_x(const datamatrix & X_n, const datamatrix & b_n);



	// FUNCTION: calc_yXb(yy, XX, bb)
	// TASK: calculates the (yy-Xxbb)'(yy-XXbb)
	double calc_yXb(const datamatrix & yy, const datamatrix & XX,
							  const datamatrix & bb) ;


	// FUNCTION: calc_yXb(y, XX, bb)
	// TASK: calculates the (y-Xxbb)'(y-XXbb)
	double calc_yXb( const datamatrix & XX, const datamatrix & bb) ;


	// FUNCTION: get_SQT_x()
	// TASK: returns the SQT_x
	double get_SQT_x(void)
	{
		return SQT_x;
	}


	// FUNCTION: calc_SQT_b()
	// TASK: calculates the SQT_b
	double calc_SQT_b(void);


	// FUNCTION: calc_SQT_b(b_n)
	// TASK: calculates the SQT_b for the proposed b_n
	double calc_SQT_b(const datamatrix & b_n);


	// FUNCTION: get_SQT_b()
	// TASK: returns the SQT_b
	double get_SQT_b(void)
	{
		return SQT_b;
	}


	// FUNCTION: log_u()
	// TASK: returns the log-value of the N(beta_mean, Sigma)
	double log_u(void);

	double log_u(const datamatrix & mean, const datamatrix & beta,
								const datamatrix & Sigma, unsigned int ncoef);


	// FUNCTION: update_sigma_i
	// TASK: draws new sigma_i from posterior distribution
	void update_sigma_i(void);


	// FUNCTION: write_to_beta
	// TASK: writes beta_help to the corresponding beta
	void write_to_beta(void);

	// FUNCTION: write_to_x
	// TASK: writes data to the corresponding x
	void write_to_x(const adja & zeta);


	// FUNCTION: write_to_x
	// TASK: writes utilities to the corresponding x
	void write_to_x(const adja & zeta, const datamatrix & uti);

	// FUNCTION: write_to_y
	// TASK: writes utilities to the corresponding y
	void write_to_y(const datamatrix & uti);


  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file (after estimation)
  void outresults(void);


    // FUNCTION: tell_var_type()
	// TARGET: tells the var_type of the response
	char tell_var_type(void)
	{
		return var_type;
	}


		// FUNCTION: num_ia_new
   // TASK: returns the number of allowed new interactions of the new main effect i
   virtual unsigned num_ia_new(unsigned i)
   {
	   return 0;
   }




   void tell_flags(void)
   {
	   cout<<flags[0]<<flags[1]<<flags[2]<<endl;
	   cout<<endl;

   }

   // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream
  virtual void outoptions(void);



  // FUNCTION: reset
  // TASK: resets all parameters for a new simulation

  void reset(void)
    {
    FULLCOND::reset();   // include this command, because the reset function
                         // of the base class automatically resets all beta
                         // matrices to their starting values (beta = 0)
                         // Sets in addition acceptance =  0;
                         //                  nrtrials = 0;

    // reset here additional variables of the inherited class
    // (e.g. additional scale parameters, auxiliary variables etc.)

    }




  // FUNCTION: predict
  // TASK: used to predict mu for a new observation Xnew
  //       computes the part of the linear predictor belonging
  //       to the full conditional and adds this part to 'linpred'

 // void predict(const datamatrix & newX, datamatrix & linpred)
 //   {

 //   assert(Xnew.cols() == data.cols());

    // add the part of the linear predictor that belongs to this full
    // conditional

 //   }


  };

  } //namespace

#endif
