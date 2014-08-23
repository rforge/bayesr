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


#if !defined (remlest_multi2_INCLUDED)

#define remlest_multi2_INCLUDED

#include "statmat.h"
#include "fullcond.h"
#include "mcmc_const.h"

#if defined(JAVA_OUTPUT_WINDOW)
#include"adminparse_basic.h"
#endif

using std::ofstream;

//------------------------------------------------------------------------------
//-------------------------- CLASS: remlest_ordinal ----------------------------
//------------------------------------------------------------------------------

class remlest_ordinal
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

  double refcat;                                   // Referenz-Kategorie
  unsigned nrcat;                                    // Anzahl Kategorien
  unsigned nrcat2;                                   // Anzahl nicht redundanter Kategorien
  statmatrix<double> cats;                         // Mögliche Kategorien
  bool catspec;                                 // are there category-specific effects
                                                // (apart from thresholds)

  unsigned totalnrfixed;
  unsigned totalnrpar;

  unsigned nrobs;
  unsigned nrobspos;

  double loglike;
  double df;
  double aic;
  double bic;
  double gcv;

  statmatrix<double> X;                         // fixed effects
  statmatrix<double> Z;                         // random effects

  vector<unsigned> xcut;                        // partition of X
  vector<unsigned> zcut;                        // partition of Z

  vector<unsigned> xcutbeta;                    // partition of fixed part in beta
  vector<unsigned> zcutbeta;                    // partition of random part in beta

  vector<bool> catspecific;                     // true if corresponding effect is
                                                // assumed to be category-specific
  vector<bool> catspecific_fixed;               // similar vector for fixed effects


  statmatrix<double> beta;                          // regression coefficients
  statmatrix<double> theta;                         // variance parameters

  public:

//------------------------------------------------------------------------------
//----------------------------- Constructor ------------------------------------
//------------------------------------------------------------------------------


  // DEFAULT CONSTRUCTOR

  remlest_ordinal(void) {}

  // Initialize from fullcond-objects

  remlest_ordinal(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adb,
  #endif
  vector<MCMC::FULLCOND*> & fc,datamatrix & re,
          const ST::string & family, const ST::string & ofile,
          const int & maxiter, const double & lowerlimit, const double & epsi,
          const double & maxch, const double & maxv,
          const datamatrix & categories, const datamatrix & weight,
          const bool & fi, ostream * lo=&cout);
//------------------------------------------------------------------------------
//--------------------------- get characteristics ------------------------------
//------------------------------------------------------------------------------

  bool get_catspec()
    {
    return catspec;
    }

//------------------------------------------------------------------------------
//----------------------------- REML estimation --------------------------------
//------------------------------------------------------------------------------


  // Function: estimate
  // Task: Perform REML-estimation with nonparametric terms

  bool estimate(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);

  // Function: estimate2
  // Task: Perform REML-estimation with nonparametric terms and
  //       category-specific effects

  bool estimate2(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);

  // Function: estimate_glm
  // Task: Perform REML-estimation without nonparametric terms

  bool estimate_glm(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);

  // Function: estimate_glm
  // Task: Perform REML-estimation without nonparametric terms and
  //       category-specific effects

  bool estimate_glm2(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight);

//------------------------------------------------------------------------------
//------------- Weights, expectation, linear predictor, etc --------------------
//------------------------------------------------------------------------------

  // FUNCTION: compute_respind
  // TASK: Computes the indicator version of the response

  void compute_respind(const datamatrix & re, datamatrix & respind);

  // FUNCTION: compute_weights
  // TASK: Computes mu, weights and the working observations

  void compute_weights(datamatrix & mu, datamatrix & workweights,
                       datamatrix & worky, datamatrix & eta,
                       datamatrix & respind, const datamatrix & weights);

  // FUNCTION: compute_eta
  // TASK: Computes the linear predictor X*beta

  void compute_eta(datamatrix & eta);

  // FUNCTION: compute_eta
  // TASK: Computes the linear predictor X*beta+Z*b

  void compute_eta2(datamatrix & eta);

  // FUNCTION: compute_sscp
  // TASK: Computes the SSCP matrices X'WX and X'W worky

  void compute_sscp(datamatrix & H, datamatrix & H1, datamatrix & workweight,
                    datamatrix & worky);

  // FUNCTION: compute_sscp2
  // TASK: Computes the weighted SSCP matrices (X Z)'W(X Z) and (X Z)W worky

  void compute_sscp2(datamatrix & H, datamatrix & H1, datamatrix & workweight,
                    datamatrix & worky);

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
                     const ST::string & rname);

  bool check_pause();

  // FUNCTION: out
  // TASK: writes results to outputstream or
  //       in Output window if BORLAND_OUTPUT_WINDOW is defined

  void out(const ST::string & s,bool thick=false,bool italic = false,
           unsigned size = 12,int r=0,int g=0, int b=0);

  void outerror(const ST::string & s);

  };

#endif

//---------------------------------------------------------------------------
#ifndef remlest_multi2H
#define remlest_multi2H
//---------------------------------------------------------------------------
#endif

