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

extension of fullcond_dag_d when there are interactions;

has to be used together with
	* fullcond_rj
	* fullcond_rj_c
	* fullcond_rj_ia
	* fullcond_dag
	* fullcond_dag_d

Interactions are added or deleted while the main effects stay fix.
(The latter are changed in fullcond_rj_ia.)

The corresponding main-file is test_rj_ia.cpp

IMPORTANT: This procedure is NOT reversible.

For the future: The program could be also used, when the general model (=the edges)
is already given and one is only interested in the interactions.

****************************************************************************/

#include"../export_type.h"

#if !defined (FULLCOND_DAG_IA_MIXED_INCLUDED)

#define FULLCOND_DAG_IA_MIXED_INCLUDED

#include"mcmc.h"
#include"distribution.h"
#include"fullcond.h"
#include"fullcond_dag.h"
#include"fullcond_dag_d.h"
#include"fullcond_dag_ia.h"
#include"func_dag.h"


namespace MCMC
{


	class __EXPORT_TYPE FULLCOND_dag_ia_mixed : public FULLCOND_dag_ia
{

	protected:

		unsigned num_continuous_pa;				// number of continuous parents
		unsigned num_discrete_pa;				// number of discrete parents
												// true, when mixed (= discrete and continuous variable)




  public:


  // DEFAULT CONSTRUCTOR:

  FULLCOND_dag_ia_mixed(void) : FULLCOND_dag_ia()
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
  FULLCOND_dag_ia_mixed  ( IA * iap , double s_i, unsigned int number,
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp);

  // CONSTRUCTOR 2
  FULLCOND_dag_ia_mixed  (bool detail_ia, IA * iap, double value_a, double value_b, ST::string prio_sig, bool dags_all,
                  const datamatrix & res, double s_i, unsigned int number,
				  MCMCoptions * o, const datamatrix & d, const ST::string & t,
				  const unsigned & r, const unsigned & c,  const ST::string & fp);



  // COPY CONSTRUCTOR

  FULLCOND_dag_ia_mixed(const FULLCOND_dag_ia_mixed & fc) : FULLCOND_dag_ia(FULLCOND_dag_ia(fc))
  {
	  num_continuous_pa = fc.num_continuous_pa;
	  num_discrete_pa = fc.num_discrete_pa;
  }



  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_dag_ia_mixed & operator=(const FULLCOND_dag_ia_mixed & fc)
  {
	  if (this==&fc)
		  return *this;

	  num_continuous_pa = fc.num_continuous_pa;
	  num_discrete_pa = fc.num_discrete_pa;

	  return *this;
  }


// DESTRUCTOR

  ~FULLCOND_dag_ia_mixed(){}




    ST::string get_family(void);



	// FUNCTION initialize
	// TASK: initializes x and xx  and y for pred_mod[i] (regression model i)
	void initialize (const adja & zeta, unsigned int j);


	// FUNCTION: ini_ia
	// TASK: creates and adds the interactions of the given main effects
	void initialize_ia (const adja & zeta, unsigned int j);


	// FUNCTION: ia_of_i
    // TASK: counts the number of allowed interactions containing variable i
	unsigned ia_of_i(unsigned i);


	// FUNCTION: new_ia_of_i
	// TASK:adds ALLOWED interactions containing variable i to v
	// when i IS NOT already a main effect
	void new_ia_of_i( unsigned i, vector <vector <unsigned > > & v);


	// FUNCTION: num_ia_of_i
	// TASK: returns the number of allowed interactions of the existing main effect i
	unsigned num_ia_of_i(unsigned i);


	// FUNCTION: num_ia_new
   // TASK: returns the number of allowed new interactions of the new main effect i
   unsigned num_ia_new(unsigned i);



   // FUNCTION: change
   // TASK: changes the values after step has been accepted
   void change(unsigned i, const datamatrix & beta_help_new, const datamatrix & x_new,
								const datamatrix & xx_new, unsigned int ncoef_new);


/*

	// FUNCTION: create_matrices
	// TASK: makes nothing, is called by fullcond_dag
	void create_matrices (void);









  // FUNCTION update
  // TASK: updates parameters (i.e. matrix beta and beta_help)
  void update(void);
*/

  };

  } //namespace

#endif
