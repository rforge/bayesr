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



#ifndef randomeffect_stepwiseH
#define randomeffect_stepwiseH

#include"../export_type.h"
#include"randomeffect.h"


namespace MCMC
{


//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_random_stepwise ----------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_random_stepwise : public FULLCOND_random
  {


  protected:

  FULLCOND_nonp_basis * fbasisp;

  double intercept;

  datamatrix data_varcoeff_fix;
  datamatrix effmodi;
  double df_unstruct;

  FULLCOND fc_df;
  ST::string utype;   // gaussian || iwlsmode

  public:

  // DEFAULT CONSTRUCTOR:

  FULLCOND_random_stepwise(void) : FULLCOND_random()
    {
    }

  // CONSTRUCTOR1
  // random intercept

  FULLCOND_random_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & d, const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const double & la, const unsigned & c=0);

  // CONSTRUCTOR2
  // random slope

  FULLCOND_random_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & intvar,const datamatrix & effmod,
                  const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const ST::string & prf,
                  const double & la, const bool & inclfix,
                  const unsigned & c=0);


  // COPY CONSTRUCTOR

  FULLCOND_random_stepwise(const FULLCOND_random_stepwise & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_random_stepwise & operator=(
                        const FULLCOND_random_stepwise & fc);

  // DESTRUCTOR

  ~FULLCOND_random_stepwise() {}


  void set_nofixed(bool fix);

  bool posteriormode(void);

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect
  ST::string get_effect(void);

  void init_names(const vector<ST::string> & na);

  // FUNCTION: reset_effect
  // TASK: resets the effect, subtracts the current effect from linearpred
  void reset_effect(const unsigned & pos);

  void compute_lambdavec(vector<double> & lvec, int & number);

  void create_weight(datamatrix & w);

  void hierarchie_fix(vector<double> & untervector, int dfo);

  void update_fix_effect(double & intercept);

  void const_varcoeff(void);

  void hierarchical(ST::string & possible);

  void set_pointer_to_interaction(FULLCOND * inter);

  void get_interactionspointer(vector<FULLCOND*> & inter);

  const datamatrix & get_data_forfixedeffects(void);

  // FUNCTION: compute_df
  // TASK: returns the approximate degrees of freedom of a smoother
  double compute_df(void);

//  double compute_df_andererteil(void);

//  void set_gleichwertig(const bool & gleich, bool weiter);

  void set_dfunstruct(const double & df_unstr);

  // FUNCTION: update_stepwise
  // TASK: returns (usually) the current smoothing parameter
  void update_stepwise(double la);

  double get_lambda(void);

  void update(void);

  void set_utype(void)
    {
    utype = "iwlsmode";
    }

  void update_gauss(void);

  void update_nongauss(void);

  void update_spatialtotal(void);

  void update_bootstrap(const bool & uncond=false);

  void save_betamean(void);

  void update_bootstrap_betamean(void);

  void update_bootstrap_df(void);

  void outresults_df(unsigned & size);

  void change_Korder(double lamb);

  void undo_Korder(void);

  void init_spatialtotal(FULLCOND_nonp_basis * sp,const ST::string & pnt,
                         const ST::string & prt);

  };     // end: class FULLCOND_random_stepwise

}   // end: namespace MCMC


#endif

