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



#ifndef fullcond_nonp_gaussian_stepwiseH
#define fullcond_nonp_gaussian_stepwiseH

#include"../export_type.h"
#include"fullcond_nonp_gaussian.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_nonp_gaussian_stepwise ---------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_nonp_gaussian_stepwise : public FULLCOND_nonp_gaussian
  {


  protected:

  double intercept;
  double df_lambdaold;
  double lambdaold;
  double lambdaold_unstr;
  double df_lambdaold_unstr;

  datamatrix data_varcoeff_fix;
  datamatrix effmodi;
  datamatrix XVX;

  FULLCOND * fcunstruct;
  bool spatialtotal;

  vector<envmatdouble> all_precenv;      // vector of all possible (X'X + lambda_i P)
  vector<double> lambdavec;

  FULLCOND fc_df;
  bool isbootstrap;

  envmatdouble Kenv2;
  envmatdouble Kenv3;
  vector<double> kappa;
  vector<double> kappaold;
  vector<double> kappa_prec;
  vector<FULLCOND *> otherfullcond;

  void init_priorassumptions(const ST::string & na);

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_nonp_gaussian_stepwise(void) : FULLCOND_nonp_gaussian()
    {
    }


  // additive Effekte, RW1 RW2 und season

  FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                      DISTRIBUTION * dp,
                      const datamatrix & d,
                      FULLCOND_const * fcc,
                      const unsigned & maxint,const fieldtype & ft,
                      const ST::string & ti,
                      const ST::string & fp, const ST::string & pres,
                      const unsigned & c,const double & l,
                      const unsigned & per=12);

  // varying coefficients , RW1 RW2 und season

  FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                       const datamatrix & d,
                       const datamatrix & intvar,
                       FULLCOND_const * fcc,
                       const unsigned & maxint,
                       const fieldtype & ft,const ST::string & ti,
                       const ST::string & fp, const ST::string & pres,
                       const unsigned & c,const double & l, const bool & vccent,
                       const unsigned & per=12);

  // spatial covariates

  FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                        DISTRIBUTION * dp,const datamatrix & d,
                        FULLCOND_const * fcc,
                        const MAP::map & m, const ST::string & mn,
                        const ST::string & ti,
                        const ST::string & fp, const ST::string & pres,
                        const unsigned & c,const double & l,
                        const fieldtype & ft, const MAP::map & m2 = MAP::map());

  // varying coefficients , spatial covariates as effect modifier

  FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                        DISTRIBUTION * dp,
                        FULLCOND_const * fcc,
                        const MAP::map & m,
                        const ST::string & mn,
                        const datamatrix & d,
                        const datamatrix & d2,
                        const ST::string & ti,
                        const ST::string & fp, const ST::string & pres,
                        const unsigned & c, const double & l, const bool & vccent,
                        const fieldtype & ft);

  // COPY CONSTRUCTOR

  FULLCOND_nonp_gaussian_stepwise(const FULLCOND_nonp_gaussian_stepwise & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_nonp_gaussian_stepwise & operator=(const FULLCOND_nonp_gaussian_stepwise & fc);


// --------------------------- FOR STEPWISE ------------------------------------

  double compute_df(void);

  void update_stepwise(double la);

  double get_lambda(void)
    {
    return lambda;
    }

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect

  ST::string get_effect(void);

  void init_names(const vector<ST::string> & na);

  void reset_effect(const unsigned & pos);

  void hierarchie_rw1(vector<double> & untervector, int dfo);

  void compute_lambdavec(vector<double> & lvec, int & number);

  void update_fix_effect(double & intercept);

  void const_varcoeff(void);

  void hierarchical(ST::string & possible);

  void set_pointer_to_interaction(FULLCOND * inter);

  void get_interactionspointer(vector<FULLCOND*> & inter);

  void init_spatialtotal(FULLCOND * unstructp);

  const datamatrix & get_data_forfixedeffects(void);

  void create_weight(datamatrix & w);

  void update_bootstrap(const bool & uncond=false);

  void save_betamean(void);

  void update_bootstrap_betamean(void);

  void update(void);

  void update_gauss(void);

  void update_IWLS(void);

  void update_bootstrap_df(void);

  void outresults_df(unsigned & size);

  void change_Korder(double lam);

  void undo_Korder(void);

  void outresults(void);

// ------------------------- END: FOR STEPWISE ---------------------------------

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  bool posteriormode_kombi(void);

  double compute_df_kombi(void);

  void set_otherfullcond(FULLCOND * ofullc)
    {
    otherfullcond.push_back(ofullc);
    }

  //void search_for_interaction(void);

  //void hierarchical(ST::string & possible);


  // DESTRUCTOR

  ~FULLCOND_nonp_gaussian_stepwise() {}

  };



} // end: namespace MCMC

#endif
