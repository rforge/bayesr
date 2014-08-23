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

extension of fullcond_rj_int when there are interactions;

has to be used together with
	* fullcond_rj
	* fullcond_rj_int
	* fullcond_dag
	* fullcond_dag_d
	* fullcond_dag_ia_mixed

fullcond_rj_mix always incorporates all possible interaction of the main effects.

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

#if !defined (FULLCOND_RJ_MIX_INCLUDED)

#define FULLCOND_RJ_MIX_INCLUDED

#include"mcmc.h"
#include"distribution.h"
#include"fullcond.h"
#include"fullcond_rj.h"
#include"fullcond_rj_int.h"
#include"fullcond_dag.h"
#include"fullcond_dag_ia.h"
#include"fullcond_dag_ia_mixed.h"
#include"adjacency.h"
#include"functions.h"


namespace MCMC
{




class __EXPORT_TYPE FULLCOND_rj_mix : public FULLCOND_rj_int
{


  protected:

	  //vector <FULLCOND_dag_d_ia *> preg_mods;




  public:


  // DEFAULT CONSTRUCTOR:
  FULLCOND_rj_mix(void) : FULLCOND_rj_int() {}


  // CONSTRUCTOR_1
  FULLCOND_rj_mix( vector < FULLCOND_dag_ia * > dagp,
			   MCMCoptions * o, const datamatrix & d,
               const ST::string & t, const unsigned & r, const unsigned & c,
               const ST::string & fp);

   // CONSTRUCTOR_1
  FULLCOND_rj_mix( vector < FULLCOND_dag_ia_mixed * > dagp,
			   MCMCoptions * o, const datamatrix & d,
               const ST::string & t, const unsigned & r, const unsigned & c,
               const ST::string & fp);


  // CONSTRUCTOR_2
  FULLCOND_rj_mix (ST::string fix, const ST::string & rp, unsigned int lim, double alph, ST::string switch_t,
					ST::string print_mod, unsigned & type,
					vector < FULLCOND_dag_ia * > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp);

   // CONSTRUCTOR_2
  FULLCOND_rj_mix (ST::string fix, const ST::string & rp, unsigned int lim, double alph, ST::string switch_t,
					ST::string print_mod, unsigned & type,
					vector < FULLCOND_dag_ia_mixed * > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp);



  // COPY CONSTRUCTOR
  FULLCOND_rj_mix(const FULLCOND_rj_mix & fc);


  // OVERLOADED ASSIGNMENT OPERATOR
  const FULLCOND_rj_mix & operator=(const FULLCOND_rj_mix & fc);


  // DESTRUCTOR
  ~FULLCOND_rj_mix() {}



  void update(void);





  };

  } //namespace

#endif

