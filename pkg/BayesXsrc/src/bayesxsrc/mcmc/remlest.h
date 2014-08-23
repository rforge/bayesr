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


#if !defined (remlest_INCLUDED)

#define remlest_INCLUDED

#include "statmat.h"
#include "fullcond.h"
#include "distribution.h"

#if defined(JAVA_OUTPUT_WINDOW)
#include<adminparse_basic.h>
#endif

using std::ofstream;

//------------------------------------------------------------------------------
//------------------------------ CLASS: remlest --------------------------------
//------------------------------------------------------------------------------

class remlest
  {

//------------------------------------------------------------------------------
//----------------------------- Variables --------------------------------------
//------------------------------------------------------------------------------

  private:

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adminb_p;
  #endif

  ostream * logout;                         // Pointer to filestream for writing
                                            // output

  vector<MCMC::FULLCOND*> fullcond;

  ST::string respfamily;                         // distribution of the response

  ST::string outfile;
  int maxit;
  double lowerlim;
  double eps;
  double maxchange;
  double maxvar;

  bool fisher;                                  // store final fisher infor-matrix?

  bool constlambda;        // fix smoothing parameters at their starting values?
  bool constscale;         // fix scale parameter at its starting value?

  statmatrix<double> X;                         // fixed effects
  statmatrix<double> Z;                         // random effects

  vector<unsigned> xcut;                        // partition of X
  vector<unsigned> zcut;                        // partition of Z

  statmatrix<double> beta;                          // regression coefficients
  statmatrix<double> theta;                         // variance parameters

  int nrint;                               // Anzahl intervallzensierter,
  int nrright;                             // rechtszensierter,
  int nrlefttrunc;                         // linkstrunkierter,
  int nruncens;                            // unzensierter Beobachtungen

  int nrobspos;                           // Beobachtungen mit positivem Gewicht

  double loglike;
  double df;
  double aic;
  double bic;
  double gcv;

  public:

//------------------------------------------------------------------------------
//----------------------------- Constructor ------------------------------------
//------------------------------------------------------------------------------


  // DEFAULT CONSTRUCTOR

  remlest(void) {}

  // CONSTRUCTOR 1: Initialize from fullcond-objects

  remlest(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adb,
  #endif
  vector<MCMC::FULLCOND*> & fc,datamatrix & re,bool dispers,
          const ST::string & family, const ST::string & ofile,
          const int & maxiter, const double & lowerlimit, const double & epsi,
          const double & maxch, const double & maxv, const bool & fi,
          const bool & cl, const bool & cs,
          ostream * lo=&cout);

//------------------------------------------------------------------------------
//----------------------------- REML estimation --------------------------------
//------------------------------------------------------------------------------


  // Function: estimate
  // Task: Perform REML-estimation without dispersion parameter
  //       returns true if an error or user break occured

  bool estimate(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);

  // Function: estimate_dispers
  // Task: Perform REML-estimation with dispersion parameter
  //       returns true if an error or user break occured

  bool estimate_dispers(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);

  // Function: estimate_glm
  // Task: compute estimates if only fix effects are present (without dispersion
  //       parameter)

  bool estimate_glm(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);

  // Function: estimate_glm_dispers
  // Task: compute estimates if only fix effects are present (with dispersion
  //       parameter)

  bool estimate_glm_dispers(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);

  // Function: estimate_survival
  // Task: compute estimates for survival data

  bool estimate_survival(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);

  // Function: estimate_survival_interval
  // Task: compute estimates for survival data in the presence of interval
  //       censoring

  bool estimate_survival_interval(datamatrix resp,
                const datamatrix & offset, const datamatrix & weight);

  // Function: estimate_survival_interval2
  // Task: compute estimates for survival data in the presence of interval
  //       censoring and left truncation

  bool estimate_survival_interval2(datamatrix resp,
                const datamatrix & offset, const datamatrix & weight,
                const bool & aiccontrol);

  // Function: estimate_aft
  // Task: compute estimates for survival data based on an AFT model with
  // smoothed error distribution

  bool estimate_aft(datamatrix resp, const datamatrix & offset,
                    const datamatrix & weight);

  // Function: estimate_aft
  // Task: compute estimates for survival data based on an AFT model with
  // smoothed error distribution if only fixed effects are present

  bool estimate_aft_glm(datamatrix resp, const datamatrix & offset,
                    const datamatrix & weight);

//------------------------------------------------------------------------------
//----------------------------- Object description -----------------------------
//------------------------------------------------------------------------------

  void outoptions();

//------------------------------------------------------------------------------
//----------------------------- Writing results --------------------------------
//------------------------------------------------------------------------------

  void make_plots(ofstream & outtex,ST::string path_batch,
                  ST::string path_splus);

  void make_model(ofstream & outtex, const ST::string & rname);

  void make_predictor(ofstream & outtex);

  void make_prior(ofstream & outtex);

  void make_options(ofstream & outtex);

  void make_fixed_table(ofstream & outtex);

  void make_graphics(const ST::string & title,
                     const ST::string & path_batch,
                     const ST::string & path_tex,
                     const ST::string & path_splus,
                     const ST::string & rname,
                     const bool & dispers);

  bool check_pause();

  // FUNCTION: out
  // TASK: writes results to outputstream or
  //       in Output window if BORLAND_OUTPUT_WINDOW is defined

  void out(const ST::string & s,bool thick=false,bool italic = false,
           unsigned size = 12,int r=0,int g=0, int b=0);

  void outerror(const ST::string & s);

//------------------------------------------------------------------------------
//------------------------------- Miscellanea ----------------------------------
//------------------------------------------------------------------------------

  };

#endif





