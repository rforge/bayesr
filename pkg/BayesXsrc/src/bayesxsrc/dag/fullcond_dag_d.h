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

#if !defined (FULLCOND_DAG_D_INCLUDED)

#define FULLCOND_DAG_D_INCLUDED

#include"mcmc.h"
#include"distribution.h"
#include"fullcond.h"
#include"fullcond_dag.h"



namespace MCMC
{


class __EXPORT_TYPE FULLCOND_dag_d : public FULLCOND_dag
{

	protected:

		datamatrix y_true;			// the real response
		//vector <interact> * pia;	// pointer to the vector which stores all interaction terms


  public:


  // DEFAULT CONSTRUCTOR:

  FULLCOND_dag_d(void) : FULLCOND_dag()
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
  FULLCOND_dag_d (double s_i, unsigned int num,
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp);

  // CONSTRUCTOR 2
  FULLCOND_dag_d (double value_a, double value_b, ST::string prio_sig, bool dags_all,
                  const datamatrix & res, double s_i, unsigned int num,
				  MCMCoptions * o, const datamatrix & d, const ST::string & t,
				  const unsigned & r, const unsigned & c,  const ST::string & fp);



  // COPY CONSTRUCTOR

  FULLCOND_dag_d(const FULLCOND_dag_d & fc) : FULLCOND_dag(FULLCOND_dag(fc))
  {
	  y_true = fc.y_true;
  }



  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_dag_d & operator=(const FULLCOND_dag_d & fc)
  {
	  if (this==&fc)
		  return *this;

	  y_true = fc.y_true;

	  return *this;
  }


// DESTRUCTOR

  ~FULLCOND_dag_d(){}





  virtual ST::string get_family(void);


  // FUNCTION ini_dag_d
  // TASK: initializes y_true
  void ini_dag_d(void);


  // FUNCTION draw_utilities
  // TASK: draws the utilities and sets them equal to y
  void draw_utilities(void);


  // FUNCTION update
  // TASK: updates parameters (i.e. matrix beta and beta_help)
  void update(void);





  };

  } //namespace

#endif
