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



#include "../export_type.h"
#if !defined (FULLCOND_RJ_INCLUDED)

#define FULLCOND_RJ_INCLUDED

#include"mcmc.h"
#include"distribution.h"
#include"fullcond.h"
#include"fullcond_dag.h"
#include"adjacency.h"
#include"func_dag.h"


namespace MCMC
{




struct modfreq
{

	ST::string model;
	unsigned nedges;
	unsigned freq;    // absolute frequency


	modfreq(void)
	{
		model="";
		nedges=0;
		freq=0;
	}


	modfreq(const ST::string &  m, const unsigned & n, const unsigned & f)
	{
		model=m;
		nedges=n;
		freq=f;
	}

	modfreq(const modfreq & m)
	{
		model = m.model;
		nedges = m.nedges;
		freq = m.freq;
	}


    const modfreq & operator=(const modfreq & m)
	{
		if (this == &m)
			return *this;

		model = m.model;
		nedges = m.nedges;
		freq=m.freq;

		return *this;
	}



	// "<" (smaller) in the sense of smaller frequency
	int operator<(const modfreq & m1) const
	{
		if (freq < m1.freq)
			return 1;
		return 0;
	}


	friend int operator>(const modfreq & m1, const modfreq & m2)
	{
		return m2 < m1;
	}




	// "<<" (smaller) in the sense of lexicographical ordering of the models
	int operator<<(const modfreq & m1) const
	{
		if (model < m1.model)
			return 1;
		return 0;
	}


	friend int operator>>(const modfreq & m1, const modfreq & m2)
	{
		return m2 << m1;
	}

	};


class __EXPORT_TYPE FULLCOND_rj : public FULLCOND
  {

  protected:

  // add here additional variables needed for the new class
  // add here additional functions/methods needed for the inherited class
  // (e.g. a method for centering beta)

	  unsigned int nvar;			// number of variables
	  unsigned int nobs;			// number of observations
	  unsigned int type;			// denotes with which type of graph to start

	  adja zeta;					// represents the graph
	  adja zeta_fix;				// represents conditions that have been made on the graph

	  vector <FULLCOND_dag *> preg_mods;	// vector of the regression models
	  vector < modfreq > freq;				// vector, that stores the chosen models and their frequencies
	  vector < essfreq > list_ess;			// vector, with the essential models of freq


	  double sigma_prop; //Varianz of the normaldistribution from which the proposal u is chosen

	  unsigned long acceptance_b;		// number of accepted birth-steps
	  unsigned long acceptance_d;		// number of accepted death-steps
	  unsigned long acceptance_s;		// number of accepted switch-steps

	  unsigned long nrtrials_b;			// number of trials for birth-step
	  unsigned long nrtrials_d;			// number of trialsf for death-step
	  unsigned long nrtrials_s;			// number of trials for switch-step

	  double alpha_sig_i;		// parameter alpha of the IG-distribution from
								// which the new variance sigma_new_i is sampled in the switch step
	  double beta_sig_i;		// parameter beta of the IG-distribution from
								// which the proposed variance sigma_new_i is sampled in the switch step

	  double alpha_sig_j;		// parameter alpha of the IG-distribution from
								// which the new variance sigma_new_j is sampled in the switch step
	  double beta_sig_j;		// parameter beta of the IG-distribution from
								// which the proposed variance sigma_new_j is sampled in the switch step

	  double gamma_nobs1;		// log(gamma(0.5*nobs))
	  double gamma_nobs2;		// log(gamma(0.5*(nobs+1)))

	  bool step_aborted;		// equal to TRUE if cycles detected

	  datamatrix zetamean;
	  datamatrix zetameanold;

	  unsigned int limit_number;// maximal number of models which are written out
								// when cum. rel. frequencies of the first "number" models is less than 1-alpha

	  double alpha;
	  ST::string switch_type;	// switch_type=normal  ->  normal switch step
                                // switch_type=equi  ->  only switching to equivalent classes
							    // switch_type=mix  ->  mixture of nomal and equi

	  ST::string print_models;  // determines how many dags are printed in the output
								// possible options: "normal", "limit", "prob", "all"

	  bool mixed_case;
	  bool file_of_results;
	  bool conditions;

	  ST::string path_res;		// name of path where results are stored

      ST::string family;       // continuous, binary with(out) interaction or mixed






  public:


  // DEFAULT CONSTRUCTOR:
  FULLCOND_rj(void) : FULLCOND() {}


  // CONSTRUCTOR_1
  FULLCOND_rj( MCMCoptions * o, const datamatrix & d,
               const ST::string & t, const unsigned & r, const unsigned & c,
               const ST::string & fp);


  // CONSTRUCTOR_2
  FULLCOND_rj( vector < FULLCOND_dag * > dagp,
			   MCMCoptions * o, const datamatrix & d,
               const ST::string & t, const unsigned & r, const unsigned & c,
               const ST::string & fp);


  // CONSTRUCTOR_3
  FULLCOND_rj (ST::string fix, const ST::string & rp, unsigned int lim, double alph,
				ST::string switch_t, ST::string print_modc, unsigned & type,
				vector < FULLCOND_dag * > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp);


  // COPY CONSTRUCTOR
  FULLCOND_rj(const FULLCOND_rj & fc);


  // OVERLOADED ASSIGNMENT OPERATOR
  const FULLCOND_rj & operator=(const FULLCOND_rj & fc);


  // DESTRUCTOR
  ~FULLCOND_rj() {};


  // FUNCTION: rj_step
  // chooses randomly new edge, decides which kind of step and makes it
  void rj_step(void);


  // FUNCTION: birth_step
  // makes birth step
  virtual void birth_step(unsigned int v_i, unsigned int v_j);


  // FUNCTION: death_step
  // makes death step, tries to delete edge i->j
  //virtual void death_step(unsigned int i, unsigned int j)
  //{
  //}

  // FUNCTION: death_step
  // makes death step, tries to delete edge i->j
  virtual void death_step(unsigned int i, unsigned int j);


  // FUNCTION: switch_step
  // makes switch step from j->i to i->j
  virtual void switch_step(unsigned int i, unsigned int j);
//  {
//  }



  // FUNCTION: switch_version_1
  // makes switch step without regarding equivalences
  virtual void switch_version_1(unsigned i, unsigned j) ;


  // FUNCTION: switch_version_2
  // makes switch step and considers equivalences
 virtual void switch_version_2(unsigned i, unsigned j) ;





	 // FUNCTION: ratio_s
	// TARGET: computes acceptance ratio in the switch step
	 double ratio_s(unsigned int i,unsigned int j,
							const datamatrix & b_new_i, const datamatrix & b_new_j,
							const datamatrix & x_new_i, const datamatrix & x_new_j,
							const datamatrix & mean_i, const datamatrix & mean_j,
							const datamatrix & sig_mean_i, const datamatrix & sig_mean_j,
							const datamatrix & xx_new_i, const datamatrix & xx_new_j,
							double sigma_new_i, double sigma_new_j);




	// FUNCTION: ratio_s
	// TARGET: computes acceptance ratio in the switch step with interactions ;
	 virtual double ratio_s_int (unsigned int i,unsigned int j,
							const datamatrix & b_new_i, const datamatrix & b_new_j,
							const datamatrix & x_new_i, const datamatrix & x_new_j)
	 {
		 return 0;
	 };



  // FUNCTION: make_new_b
  // TASK: computes the new values for a birth-step
  void make_new_b (ST::string step, unsigned int i, unsigned int j, double beta_new,
				datamatrix & xx_new, datamatrix & b_new, datamatrix & x_new);




  // FUNCTION: make_new_d
  // TASK: computes the new values for a death-step
  void make_new_d (ST::string step, unsigned int i, unsigned int j,datamatrix & xx_new,
				double & beta_old, datamatrix & b_new, datamatrix & x_new);



  // FUNCTION: sample_sigma
  // TARGET: samples the new variance of the regression model i in the switch step
  double sample_sigma(char vertex, unsigned int i, unsigned int ncoef_new_i,
							const datamatrix & mean_i, const datamatrix & x_new_i);



  // FUNCTION: ratio_b()
  // TARGET: calculates the ratio of a birth-step
  double ratio_b(unsigned int j,  double u,
							const datamatrix & b_new, const datamatrix & x_new);


  // FUNCTION: ratio_d()
  // TARGET: calculates the ratio of a death-step
  double ratio_d(unsigned int j, double u,
					const datamatrix & b_new, const datamatrix & x_new);



  // FUNCTION: ratio_s
  // TARGET: computes acceptance ratio in the switch step
  //double ratio_s(unsigned int i,unsigned int j,
//							const datamatrix & b_new_i, const datamatrix & b_new_j,
//							const datamatrix & x_new_i, const datamatrix & x_new_j,
//							const datamatrix & mean_i, const datamatrix & mean_j,
//							double sigma_new_i, double sigma_new_j);





  // FUNCTION: log_gamma
  // TASK: returns the logarithm of the gammafunction when value=0.5*(unsigned)
  double log_gamma(double value) const;


  // FUNCTION: log_gamma1
  // TASK: returns the logarithm of the gammafunction when value=0.5*(nobs-k)
  double log_gamma1(double x) const;



  // FUNCTION: p_prop()
  // TARGET: returns the density of the proposal u, which is normaldistributed
	double p_prop(double prop);




	// FUNCTION: p_prop()
  // TARGET: returns the density of the proposal vector u, which is normal and independent
	double p_prop(const datamatrix & prop);



  // FUNCTION: accept
  // TASK: returns true with probability ratio
  //bool accept (double ratio);


  // FUNCTION update
  // TASK: updates parameters (i.e. matrix beta)
  void update(void);


  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file (after estimation)
  void outresults(void);


  // FUNCTION:  outres_dags
  // TASK: writes out the different dags
  void outres_dags(void);



  // FUNCTION: outres_essentials
  // TASK: writes out the different essential graphs
  void outres_essentials(void);

  // FUNCTION: write_out_essential
  // TASK: writes out the essential graph m and its properties
  void write_out_essential(essfreq & m, unsigned total_number);


  // FUNCTION: write_out_ess_short
  // TASK: writes out the essential graph ess into a separate file
  // void write_out_ess_short(essfreq & ess, unsigned total_number);


  // FUNCTION: write_out_resfile
  // TASK: writes out results (= 10 most important essential graphs and the adjacency matrix)
  // into a separate file which has to be named before; good for simulation studies
   void write_out_resfile(void);


  // FUNCTION: make_list_essential
 // TASK: computes the list of the essential graphs from the list of all dags
 void make_list_essential(void);

  // FUNCTION: add_ess_to_list
  // TASK: adds ess_new to ther lists and updates frequencies
  //void add_ess_to_list(vector <essfreq> & list_ess, essfreq & ess_new);

  // FUNCTION: store_model
 // TASK: stores the models in a vector
 void store_model(void);


  // FUNCTION: update_zeta
  // TASK: updates zetamean and the auxiliary variables like zeta_max etc.
  void update_zeta(void);


  // FUNCTION: setzeta
  // TASK: initializes zeatmean etc.
  void setzeta(const adja & zetanew);


  // FUNCTION: ini_structure
   // TASK: initializes structure
	void ini_structure(void);

	// FUNCTION: ini_structure
   // TASK: initializes structure
	void ini_structure(unsigned t);


  // FUNCTION: ini_ratio
  // TASK: initializes ratios.
  void ini_ratio(void);


  // FUNCTION: ini_hyperpar
  // TASK: initializes hyperparamaters etc
  void ini_hyperpar(void);


  // FUNCTION: set_options
  // TASK: sets the options to default parameter
  void set_options(void);


  // FUNCTION: conditions_okay
  // TASK: returns true if conditions are fullfilled
  bool conditions_okay(unsigned int i, unsigned int j);


	// FUNCTION: conditions_okay_d
	// TASK: returns true if conditions are fullfilled
	bool conditions_okay_d (unsigned int i, unsigned int j);


	// FUNCTION: conditions_okay_b
	// TASK: returns true if conditions are fullfilled
	bool conditions_okay_b (unsigned int i, unsigned int j);


	// FUNCTION: conditions_okay_s
	// TASK: returns true if conditions are fullfilled
	bool conditions_okay_s (unsigned int i, unsigned int j);


  // FUNCTION: set_options
  // TASK: sets the options
  void set_options( unsigned lim, double alph, ST::string switch_t,
									ST::string print_mod, ST::string fix_path);


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

     void outoptions(void) ;




 /**************


	 void predict(const datamatrix & newX, datamatrix & linpred)
    {

    assert(newX.cols() == data.cols());

    // add the part of the linear predictor that belongs to this full
    // conditional

    }
******/

  };

  } //namespace

#endif

