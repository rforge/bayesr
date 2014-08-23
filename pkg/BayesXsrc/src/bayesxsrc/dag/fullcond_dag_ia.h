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

#if !defined (FULLCOND_DAG_IA_INCLUDED)

#define FULLCOND_DAG_IA_INCLUDED

#include"mcmc.h"
#include"distribution.h"
#include"fullcond.h"
#include"fullcond_dag.h"
#include"fullcond_dag_d.h"
#include"ia.h"
#include"func_dag.h"


namespace MCMC
{


	class __EXPORT_TYPE FULLCOND_dag_ia : public FULLCOND_dag_d
{

	protected:

		IA * pia;								// pointer to a IA-object
		vector <vector <unsigned> >	current_ia;	// vector wih the current interactions of the model
		vector <vector <unsigned> > all_ia;		// contains the vectors of ALL possible interactions
												// needed for outresults

		vector<int > occurrence;			// elements are 1 if corresponding ia occurs
												// otherwise elements are zero

		datamatrix occmean;						// stores mean of occurrence


		datamatrix x_ia_d;						// help matrices for the death_step
		datamatrix xx_ia_d;
		datamatrix x_ia_b;						// help matrices for the birth_step
		datamatrix xx_ia_b;
		datamatrix y_ia;						// help matrix for death and birth-step

		bool ia_d_there;						// indicates if x_ia_d and xx_ia_d have already been created
		bool ia_b_there;						// indicates if x_ia_b and xx_ia_b have already been created


		unsigned max_ia_order;					// maximal order of interaction terms
		unsigned all_possible_ia;				// maximal number of coefficients (main + ia);
		double ln_prop_beta;					// logaritmized density of proposal

		bool detail;
		bool mixed_case;						// true, when mixed (= discrete and continuous variable)


		unsigned proposal_version;


  public:


  // DEFAULT CONSTRUCTOR:

  FULLCOND_dag_ia(void) : FULLCOND_dag_d()
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
  FULLCOND_dag_ia (IA * iap , double s_i, unsigned int number,
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp);

  // CONSTRUCTOR 2
  FULLCOND_dag_ia (bool detail_ia, IA * iap, double value_a, double value_b, ST::string prio_sig, bool dags_all,
                  const datamatrix & res, double s_i, unsigned int number,
				  MCMCoptions * o, const datamatrix & d, const ST::string & t,
				  const unsigned & r, const unsigned & c,  const ST::string & fp);

// CONSTRUCTOR 3
  FULLCOND_dag_ia (bool detail_ia, char typ, IA * iap , double s_i, unsigned int number,
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp);

  // CONSTRUCTOR 3
  FULLCOND_dag_ia (bool detail_ia, char typ, IA * iap, double value_a, double value_b, ST::string prio_sig, bool dags_all,
                  const datamatrix & res, double s_i, unsigned int number,
				  MCMCoptions * o, const datamatrix & d, const ST::string & t,
				  const unsigned & r, const unsigned & c,  const ST::string & fp);



  // COPY CONSTRUCTOR

  FULLCOND_dag_ia(const FULLCOND_dag_ia & fc) : FULLCOND_dag_d(FULLCOND_dag_d(fc))
  {
	  pia = fc.pia;
	  current_ia = fc.current_ia;
	  all_ia = fc.all_ia;
	  occurrence = fc.occurrence;
	  y_ia = fc.y_ia;
	  x_ia_b = fc.x_ia_b;
	  xx_ia_b = fc.xx_ia_b;
	  x_ia_d = fc.x_ia_d;
	  xx_ia_d = fc.xx_ia_d;
	  ia_d_there = fc.ia_d_there;
	  ia_b_there = fc.ia_b_there;
	  occmean = fc.occmean;
	  max_ia_order = fc.max_ia_order;
	  all_possible_ia = fc.all_possible_ia;
	  ln_prop_beta = fc.ln_prop_beta;
	  detail = fc.detail;
	  mixed_case = fc.mixed_case;
	  proposal_version = fc.proposal_version;
  }



  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_dag_ia & operator=(const FULLCOND_dag_ia & fc)
  {
	  if (this==&fc)
		  return *this;

	  pia = fc.pia;
	  current_ia = fc.current_ia;
	  all_ia = fc.all_ia;
	  occurrence = fc.occurrence;
	  y_ia = fc.y_ia;
	  x_ia_b = fc.x_ia_b;
	  xx_ia_b = fc.xx_ia_b;
	  x_ia_d = fc.x_ia_d;
	  xx_ia_d = fc.xx_ia_d;
	  ia_d_there = fc.ia_d_there;
	  ia_b_there = fc.ia_b_there;
	  occmean = fc.occmean;
	  max_ia_order = fc.max_ia_order;
	  all_possible_ia = fc.all_possible_ia;
	  ln_prop_beta = fc.ln_prop_beta;
	  detail = fc.detail;
	  mixed_case = fc.mixed_case;
	  proposal_version = fc.proposal_version;

	  return *this;
  }


// DESTRUCTOR

  ~FULLCOND_dag_ia(){}




  virtual ST::string get_family(void);


  // FUNCTION get_nvar
  // TASK: gives nvar
  unsigned get_nvar(void);



	// FUNCTION: get_current_ia
	// TASK: returns current_ia
	double get_ln_prop_beta (void)
	{
		  return ln_prop_beta;
	}


    // FUNCTION: get_current_ia
	// TASK: returns current_ia
	 vector <vector <unsigned > > & get_current_ia (void)
	{
		vector <vector <unsigned > > & reference = current_ia;
		return reference;
	}


	// FUNCTION: get_current_ia
	// TASK: returns k-th element of current_ia
	vector<unsigned> get_current_ia (unsigned k)
	{
		return current_ia[k];
	}



    // FUNCTION: ini_ia
	// TASK: creates and adds the interactions of the given main effects
	virtual void initialize_ia (const adja & zeta, unsigned int j);


	// FUNCTION: create_matrices
	// TASK: create_matrices for the next birth or death step
	void create_matrices (void);


	// FUNCTION: write_ia_to_x
	// TASK: adds the interaction to x when x is initialized
	// is called in FULLCOND_dag::initialize
	void write_ia_to_x(void);




  // FUNCTION: ia_of_i
  // TASK:adds interactions containing variable i to v
  void ia_of_i(unsigned i,  vector <vector <unsigned > > & v);



	// FUNCTION: ia_of_i
	// TASK:adds interactions containing variable i to v
	// when i IS NOT already a main effect
  virtual void new_ia_of_i( unsigned i, vector <vector <unsigned > > & v);


  // FUNCTION: ia_of_i
  // TASK: counts the number of interactions containing variable i
  unsigned ia_of_i(unsigned);


  // FUNCTION: ia_of_i
  // TASK: returns the number of interactions of the main effect i
  // if an interaction between all pairs of main effects is assumed
  unsigned ia_of_i(void);



  // FUNCTION update
  // TASK: updates parameters (i.e. matrix beta and beta_help)
  void update(void);


  // FUNCTION update_occ
  // TASK: updates occurrence
   void update_occ(void);


  // FUNCTION birth_step
  // TASK: tries to add an additional interaction term
  void birth_step (vector<unsigned> new_ia);


  // FUNCTION death_step
  // TASK: tries to delete an interaction term
  void death_step (vector<unsigned> old_ia);

    // FUNCTION: make_new_b
	// TASK: computes the new values for a birth-step
	void make_new_b (vector<unsigned> ia_new, double beta_new,
				datamatrix & xx_new, datamatrix & b_new, datamatrix & x_new);

	// FUNCTION: make_new_b_int
	// TASK: computes the new values for a birth-step when one new main effect
	// and ALL corresponding interactions are added
	// is called from birth_step in rj_int
	void make_new_b_int (ST::string step, unsigned i, vector <vector <unsigned> > ia_new,
					 datamatrix & beta_new, datamatrix & xx_new,
					 datamatrix & b_new, datamatrix & x_new);



   // FUNCTION: make_new_d_int
   // TASK: computes the new values for a death-step
   void make_new_d_int (ST::string step, unsigned i, unsigned j, unsigned ia_del,
						datamatrix & beta_old, vector <vector <unsigned > > & current_ia_n,
						datamatrix & xx_new, datamatrix & b_new, datamatrix & x_new);





	// FUNCTION: make_prop_beta
	// TASK: computes the proposed regression coefficient beta_new
	// and stores its new components in u
	void make_prop_beta (char step, datamatrix & beta_new,
			 datamatrix & u, const datamatrix & x_new, const datamatrix & xx_new,
			unsigned ncoef_new, const vector <unsigned> vec_t);


	// FUNCTION: make_new_d
   // TASK: computes the new values for a death-step
   void make_new_d ( vector<unsigned> ia_old, datamatrix & xx_new,
					double & beta_old, datamatrix & b_new, datamatrix & x_new);



   // FUNCTION: get_pos
	// TASK: gives back the position of main effect i and
	//		 the corresponding ia in the regression model
	//		 starting with 0 for the intercept
	void get_pos(unsigned i, vector <unsigned> & pos);


	// FUNCTION: get_pos_cur
	// TASK: gives back the position of ia in the vector current_ia
	unsigned get_pos_cur(vector <unsigned> ia);




	// FUNCTION: write_to_beta_ia()
	// TARGET: writes the coefficients of the interactions-variables to the vector beta
	void write_to_beta_ia(void);


	// FUNCTION: change_occur()
	// TARGET: adds or delets ia in occur after birth or death step
	void change_occur(char step, vector<unsigned> ia_new);


	// FUNCTION: change_occur()
	// TARGET: adds or delets ia in occur after birth or death step
	void change_occur(char step, vector <vector <unsigned > > ia_vec);



	// FUNCTION: change_current()
	// TARGET: adds or delets ia.term to/from current_ia
	void change_current(char step, vector<unsigned> term);


	// FUNCTION: change_current()
	// TARGET: adds or delets ia.term to/from current_ia
	void change_current(char step, vector <vector <unsigned > > term);


  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream
  void outoptions(void);

  void outresults(void);


	// FUNCTION: num_ia_new
    // TASK: returns the number of allowed new interactions of the new main effect i
    unsigned num_ia_new(unsigned i);






/*
  // FUNCTION normal_step
  // TASK: updates main effects and interaction term by drawing from full conditional
  void normal_step (void);

  */


  };

  } //namespace

#endif
