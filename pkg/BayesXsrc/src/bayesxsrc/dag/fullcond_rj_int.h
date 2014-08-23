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

extension of fullcond_rj when there are interactions;

has to be used together with
	* fullcond_rj
	* fullcond_dag
	* fullcond_dag_d
	* fullcond_dag_ia

fullcond_rj_int always incorporates all possible interaction of the main effects.

The updating of the ia_terms takes place in fullcond_dag_d_ia.

The corresponding main-file is test_int.cpp

IMPORTANT: This procedure is HOPEFULLY reversible.

****************************************************************************/


#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (FULLCOND_RJ_INT_INCLUDED)

#define FULLCOND_RJ_INT_INCLUDED

#include"mcmc.h"
#include"distribution.h"
#include"fullcond.h"
#include"fullcond_rj.h"
#include"fullcond_dag.h"
#include"fullcond_dag_ia.h"
#include"fullcond_dag_ia_mixed.h"
#include"adjacency.h"
#include"functions.h"


namespace MCMC
{




class __EXPORT_TYPE FULLCOND_rj_int : public FULLCOND_rj
{


  protected:

	  //vector <FULLCOND_dag_d_ia *> preg_mods;




  public:


  // DEFAULT CONSTRUCTOR:
  FULLCOND_rj_int(void) : FULLCOND_rj() {}


  // CONSTRUCTOR_1
  FULLCOND_rj_int( vector < FULLCOND_dag_ia * > dagp,
			   MCMCoptions * o, const datamatrix & d,
               const ST::string & t, const unsigned & r, const unsigned & c,
               const ST::string & fp);


   // CONSTRUCTOR_1a
  FULLCOND_rj_int( vector < FULLCOND_dag_ia_mixed * > dagp,
			   MCMCoptions * o, const datamatrix & d,
               const ST::string & t, const unsigned & r, const unsigned & c,
               const ST::string & fp);


  // CONSTRUCTOR_2
  FULLCOND_rj_int (ST::string fix, const ST::string & rp,unsigned int lim, double alph,
					ST::string switch_t, ST::string print_mod, unsigned & type,
					vector < FULLCOND_dag_ia* > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp);

  // CONSTRUCTOR_2a
  FULLCOND_rj_int (ST::string fix, const ST::string & rp, unsigned int lim,
									double alph , ST::string switch_t, ST::string print_mod,
								  unsigned & ty, vector < FULLCOND_dag_ia_mixed* > dagp,
								  MCMCoptions * o, const datamatrix & d, const ST::string & t,
								  const unsigned & r, const unsigned & c, const ST::string & fp);



  // COPY CONSTRUCTOR
  FULLCOND_rj_int(const FULLCOND_rj_int & fc);


  // OVERLOADED ASSIGNMENT OPERATOR
  const FULLCOND_rj_int & operator=(const FULLCOND_rj_int & fc);


  // DESTRUCTOR
  ~FULLCOND_rj_int() {}




  // FUNCTION: change_preg_mods
  // TASK: changes preg_mods
  void change_preg_mods( vector <FULLCOND_dag_ia *> );


  // FUNCTION: change_preg_mods
  // TASK: changes preg_mods
  void change_preg_mods( vector <FULLCOND_dag_ia_mixed* > dagp);


  // FUNCTION: death_step
  // makes death step, tries to delete edge i->j
  // and all interaction terms that contain i
  void death_step(unsigned int i, unsigned int j);



  // FUNCTION: birth_step
  // makes birth-step, tries to add edge i->j
  // and all interaction terms between i and the already existing main effects
  void birth_step(unsigned int v_i, unsigned int v_j);




   // FUNCTION: make_new_b
  // TASK: computes the new values for a birth-step
  void make_new_b (ST::string step, unsigned int i, unsigned int j, unsigned ia_new,
						datamatrix & beta_new, datamatrix & xx_new, datamatrix & b_new,
						datamatrix & x_new);


  // FUNCTION: switch_version_1
  // makes switch step and DOES NOT consider equivalences
 void switch_version_1(unsigned i, unsigned j);


 // FUNCTION: switch_version_2
  // makes switch step and considers equivalences
  void switch_version_2 (unsigned i, unsigned j);


 // FUNCTION: ratio_s
  // TARGET: computes acceptance ratio in the switch step
  double ratio_s_int (unsigned int i,unsigned int j,
							const datamatrix & b_new_i, const datamatrix & b_new_j,
							const datamatrix & x_new_i, const datamatrix & x_new_j);

   // FUNCTION: switch_step
// makes switch step from j->i to i->j
void switch_step(unsigned int i, unsigned int j);


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

