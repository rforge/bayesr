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



#ifndef fullcond_pspline_stepwiseH
#define fullcond_pspline_stepwiseH

#include"../export_type.h"
#include"fullcond_pspline_gaussian.h"
#include "fullcond_nonp_gaussian.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_pspline_stepwise ---------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FULLCOND_pspline_stepwise : public FULLCOND_pspline_gaussian
  {


  protected:

  datamatrix data_varcoeff_fix;
  datamatrix effmodi;
  datamatrix XVX;

  double df_lambdaold;
  double lambdaold;

  vector<envmatdouble> all_precenv;      // vector of all possible (X'X + lambda_i P)
  vector<double> lambdavec;

  envmatdouble Menv;
  bool concave;
  bool convex;
  double lambdamono;

  envmatdouble Kenv2;
  vector<double> kappa;
  vector<double> kappaold;
  vector<double> kappa_prec;
  vector<FULLCOND *> otherfullcond;
  envmatdouble Kenv3;

  FULLCOND fc_df;
  bool isbootstrap;
  updatetype utype;   // gaussian || iwls


  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_pspline_stepwise(void) : FULLCOND_pspline_gaussian()
    {
    }

  // CONSTRUCTOR 1  (for additive models)
  // o    : pointer to MCMCoptions object
  // dp   : pointer to DISTRIBUTION object
  // fcc  : pointer to FULLCOND_const object
  // d    : data
  // nrk  : number of knots
  // degr : degree of splines
  // kp   : position of knots (equidistant or quantiles)
  // ft   : field type (RW1, RW2)
  // monotone: increasing || decreasing || unrestricted
  // ti   : title of the object
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored
  // deriv: should the first derivative be computed?
  // l    : starting value for lambda
  // gs   : gridsize
  // diag : should the diagonal transformation be performed?
  // c    : column of the linear predictor (ususally 0)

  FULLCOND_pspline_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                          FULLCOND_const * fcc, const datamatrix & d,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                         const ST::string & fp, const ST::string & pres, const bool & deriv,
                         const double & l, const int & gs, const bool & diag, const unsigned & c=0);

  // CONSTRUCTOR 2  (for  varying coefficients term)
  // effmod: values of the effect modifier
  // intact: values of the interaction variable

  FULLCOND_pspline_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,
                         const datamatrix & effmod, const datamatrix & intact,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                         const ST::string & fp, const ST::string & pres, const bool & deriv,
                         const double & l, const int & gs, const bool & vccent,
                         const unsigned & c=0);


  // COPY CONSTRUCTOR

  FULLCOND_pspline_stepwise(const FULLCOND_pspline_stepwise & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_pspline_stepwise & operator=(const FULLCOND_pspline_stepwise & fc);


  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  bool changeposterior3(const datamatrix & betamain, const datamatrix & main, const double & inter);

  bool changeposterior_varcoeff(const datamatrix & betamain, const datamatrix & main, const double & inter);

  void reset_effect(const unsigned & pos);

  void reset(void);

  void hierarchie_rw1(vector<double> & untervector, int dfo);

  void compute_lambdavec(vector<double> & lvec, int & number);

  double compute_df(void);

  void update_stepwise(double la);

  double get_lambda(void)
    {
    return lambda;
    }

  void create_weight(datamatrix & w);

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect

  ST::string get_effect(void);

  void init_names(const vector<ST::string> & na);

  const datamatrix & get_data_forfixedeffects(void);

  void update_fix_effect(void);

  void const_varcoeff(void);

  void set_pointer_to_interaction(FULLCOND * inter);

  void get_interactionspointer(vector<FULLCOND*> & inter);

  bool search_for_interaction(void);

  void hierarchical(ST::string & possible);

  void updateMenv(void);

  void set_spmonotone(double & spmono)
    {
    lambdamono = spmono;
    }

  void update_bootstrap(const bool & uncond=false);

  void update_beta_average(unsigned & samplesize);

  void save_betamean(void);

  void update_bootstrap_betamean(void);

  void update(void);

  void update_bootstrap_df(void);

  void outresults_df(unsigned & size);

  void change_Korder(double lam);

  void undo_Korder(void);

  void get_samples(const ST::string & filename,const unsigned & step) const;

  void change_varcoeff(const datamatrix & betamain,const datamatrix & main,const double & inter);

  void update_gauss(void);

  void update_IWLS(void);

  void set_utype(void)
    {
    utype = iwls;
    }

  void outresults(void);

  vector<int>::iterator get_freqoutputit(void)
    {
    return freqoutput.begin();
    }

  void set_spline(datamatrix & sp)
    {
    spline.assign(sp);
    }

  bool posteriormode_kombi(void);

  double compute_df_kombi(void);

  void set_otherfullcond(FULLCOND * ofullc)
    {
    otherfullcond.push_back(ofullc);
    }

  // DESTRUCTOR

  ~FULLCOND_pspline_stepwise() {}

  };


} // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
