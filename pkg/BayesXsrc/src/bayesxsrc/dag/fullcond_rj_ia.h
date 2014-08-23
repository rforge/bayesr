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




/****************** DESCRIPTION *******************************************

extension of fullcond_rj_c when there are interactions;

has to be used together with
	* fullcond_rj
	* fullcond_rj_c
	* fullcond_dag
	* fullcond_dag_d
	* fullcond_dag_d_ia

fullcond_rj_ia allows a death-step where all corresponding
ia-terms are deleted, too. The birth-step is only for main
effects. There is no switch step (to complicate).

The intake of interactions (and also the deletion) takes
place in fullcond_dag_d_ia.

The corresponding main-file is test_rj_ia.cpp

IMPORTANT: This procedure is NOT reversible.

****************************************************************************/

#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FULLCOND_RJ_IA_INCLUDED)

#define FULLCOND_RJ_IA_INCLUDED

#include"mcmc.h"
#include"distribution.h"
#include"fullcond.h"
#include"fullcond_rj.h"
#include"fullcond_dag.h"
#include"fullcond_dag_ia_mixed.h"
#include"adjacency.h"
#include"functions.h"


namespace MCMC
{




class __EXPORT_TYPE FULLCOND_rj_ia : public FULLCOND_rj
  {

  protected:

	  //vector <FULLCOND_dag_d_ia *> preg_mods;




  public:


  // DEFAULT CONSTRUCTOR:
  FULLCOND_rj_ia(void) : FULLCOND_rj() {}


  // CONSTRUCTOR_1
  FULLCOND_rj_ia( vector < FULLCOND_dag_d_ia * > dagp,
			   MCMCoptions * o, const datamatrix & d,
               const ST::string & t, const unsigned & r, const unsigned & c,
               const ST::string & fp);


  // CONSTRUCTOR_2
  FULLCOND_rj_ia (unsigned int lim, double alph, ST::string swi,
					ST::string print_mod, unsigned & type,
					vector < FULLCOND_dag * > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp);


  // COPY CONSTRUCTOR
  FULLCOND_rj_ia(const FULLCOND_rj_ia & fc);


  // OVERLOADED ASSIGNMENT OPERATOR
  const FULLCOND_rj_ia & operator=(const FULLCOND_rj_ia & fc);


  // DESTRUCTOR
  ~FULLCOND_rj_ia() {}




  // FUNCTION: change_preg_mods
  // TASK: changes preg_mods
  void change_preg_mods( vector <FULLCOND_dag_d_ia *> );


  // FUNCTION: death_step
  // makes death step, tries to delete edge i->j
  void death_step(unsigned int i, unsigned int j);



   // FUNCTION: make_new_d_ia
   // TASK: computes the new values for a death-step
   void make_new_d_ia (ST::string step, unsigned i, unsigned j, unsigned ia_del,
						datamatrix & beta_old, vector <vector <unsigned > > & current_ia_n,
						datamatrix & xx_new, datamatrix & b_new, datamatrix & x_new);


/*
  // FUNCTION: rj_step
  // chooses randomly new edge, decides which kind of step and makes it
  void rj_step(void);


  // FUNCTION: birth_step
  // makes birth step
  void birth_step(unsigned int i, unsigned int j);


  // FUNCTION: death_step
  // makes death step, tries to delete edge i->j
  void death_step(unsigned int i, unsigned int j);


  // FUNCTION: switch_step
  // makes switch step from j->i to i->j
  void switch_step(unsigned int i, unsigned int j);



  // FUNCTION: make_new_b
  // TASK: computes the new values for a birth-step
  void FULLCOND_rj::make_new_b (ST::string step, unsigned int i, unsigned int j, double beta_new,
				datamatrix & xx_new, datamatrix & b_new, datamatrix & x_new);


  // FUNCTION: make_new_d
  // TASK: computes the new values for a death-step
  void make_new_d (ST::string step, unsigned int i, unsigned int j,datamatrix & xx_new,
				double & beta_old, datamatrix & b_new, datamatrix & x_new);


  // FUNCTION: sample_sigma
  // TARGET: samples the new variance of the regression model i in the switch step
  double FULLCOND_rj::sample_sigma(char vertex, unsigned int i, unsigned int ncoef_new_i,
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
  double ratio_s(unsigned int i,unsigned int j,
							const datamatrix & b_new_i, const datamatrix & b_new_j,
							const datamatrix & x_new_i, const datamatrix & x_new_j,
							const datamatrix & mean_i, const datamatrix & mean_j,
							double sigma_new_i, double sigma_new_j);





  // FUNCTION: log_gamma
  // TASK: returns the logarithm of the gammafunction when value=0.5*(unsigned)
  double log_gamma(double value) const;


  // FUNCTION: log_gamma1
  // TASK: returns the logarithm of the gammafunction when value=0.5*(nobs-k)
  double FULLCOND_rj::log_gamma1(double x) const;



  // FUNCTION: p_prop()
  // TARGET: returns the density of the proposal u, which is normaldistributed
	double p_prop(double prop);


  // FUNCTION: accept
  // TASK: returns true with probability ratio
  //bool accept (double ratio);


  // FUNCTION update
  // TASK: updates parameters (i.e. matrix beta)
  void update(void);


  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file (after estimation)
  void outresults(void);

  // FUNCTION: store_model
 // TASK: stores the models in a vector
 void FULLCOND_rj::store_model(void);


  // FUNCTION: update_zeta
  // TASK: updates zetamean and the auxiliary variables like zeta_max etc.
  void update_zeta(void);


  // FUNCTION: setzeta
  // TASK: initializes zeatmean etc.
  void setzeta(const adja & zetanew);


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



*/


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

