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



#ifndef mcmc_psplineH
#define mcmc_psplineH

#include"../export_type.h"
#include "mcmc.h"
#include "fullcond.h"
#include "time.h"
#include<deque>
#include "sparsemat.h"
#include "mcmc_nonpbasis.h"
#include "spline_basis.h"


namespace MCMC
{


//------------------------------------------------------------------------------
//---------------------------- class: FULLCOND_pspline -------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_pspline : public spline_basis
  {


  protected:

  bool maxtoobig;
  bool mintoobig;
  unsigned min;                      // Minimum Blocksize
  unsigned max;                      // Maximum Blocksize

  unsigned minauto;
  unsigned maxauto;
  bool automatic;
  unsigned oldacceptance;
  unsigned oldnrtrials;

  SparseMatrix K;

  vector<datamatrix> KAB;            // vector of all possible K_ab^-1
  vector<datamatrix> KABroot;        // vector of all possible sqrt(K_ab^-1)

  vector<SparseMatrix> KABr_sp;      // vector of all possible K_left matrices
                                     // (sparse matrix version)
  vector<SparseMatrix> KABl_sp;      // vector of all possible K_right matrices
                                     // (sparse matrix version)

  vector<unsigned> begin;            // begin[size-min] gives the Position of
                                     // the first KAB,KABl,KABr,KABroot Matrices
                                     // for blocksize 'size'
  vector<unsigned> matquant;         // matquant[size-min] gives the number of
                                     // blocks for blocksize 'size'

  vector<datamatrix> fc_random;
  vector<datamatrix> randnorm;


  // FUNCTION: make_Kab_list
  // TASK: creates KAB, KAbroot,KABr,KABl matrices
  // NOTE: functions depends on:
  //       - Penalty matrix 'K'
  //       - Minimum and Maximum Block size 'min' and 'max'
  //       - Type of MRF 'type'

  void make_Kab_list(void);

  void adjust_blocksize(const unsigned & alphamin,const unsigned & alphamax);



  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_pspline(void) : spline_basis()
    {
    }

  // CONSTRUCTOR 1

  FULLCOND_pspline(MCMCoptions * o,DISTRIBUTION * dp,
                    FULLCOND_const * fcc,const datamatrix & d,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const double & lk, const double & uk,
                    const double & lg, const double & ug,
                    const unsigned & c);

  // CONSTRUCTOR 2 für variierende Koeffizienten

  FULLCOND_pspline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                   const datamatrix & effmod,const datamatrix & intact,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const unsigned & c);

  // CONSTRUCTOR 3 für Cox Modell (nur um spline_basis aufzurufen)

  FULLCOND_pspline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const fieldtype & ft,const ST::string & ti,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const unsigned & c);


  // COPY CONSTRUCTOR

  FULLCOND_pspline(const FULLCOND_pspline & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_pspline & operator=(const FULLCOND_pspline & fc);

  void update(void);

  void outresults(void);

  void outoptions(void);

  bool posteriormode(void)
    {
    return true;
    }

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND::reset();
    }

  // FUNCTION: predict
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred);

  // FUNCTION: confidencebands
  // TASK: computes confidence bands for the estimated values and writes them into
  //       the file 'filename'

  void confidencebands(const ST::string & filename);

  // FUNCTION: compute_mse
  // TASK: computes the MSE for a model with known true function

  double compute_mse(const datamatrix & m);

  // FUNCTION: compute mu
  // TASK: computes the mean of the conditional prior f[a,b]
  //       'beta' is the current state of the Markov chain
  //        bs = blocksize
  //        v = column of beta

  void compute_mu(const datamatrix & beta, const unsigned & bs,
                  const unsigned & a,const unsigned & b,
                  const unsigned & v);

  // FUNCTION: compute_fc
  // TASK:
  // ADDITIONAL NOTE: implementation is independent of the type of MRF

  void compute_fc(const datamatrix & beta, const unsigned & bs,
                  const unsigned & a,const unsigned & b,
                  const double & Q,const unsigned & v=0);

  // FUNCTION: compute_quadform
  // TASK: returns beta(.,v)' K beta(.,v) where K is the penalty matrix

  double compute_quadform(void);

  // DESTRUCTOR

  ~FULLCOND_pspline() {}

  };



}   // end: namespace MCMC

#endif
