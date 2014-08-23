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



#ifndef mcmc_pspline_surfH
#define mcmc_pspline_surfH

#include"../export_type.h"
#include "mcmc.h"
#include "fullcond.h"
#include "mcmc_nonp.h"
#include "time.h"
#include<deque>
#include "sparsemat.h"
#include "mcmc_pspline.h"
#include "bandmat.h"
#include "bandmat_penalty.h"
#include "spline_basis_surf.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//---------------------------- class: FULLCOND_pspline_surf --------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_pspline_surf : public spline_basis_surf
  {


  protected:

  unsigned min;
  unsigned max;
  bool mintoobig;
  bool maxtoobig;

  unsigned minauto;
  unsigned maxauto;
  bool automatic;
  unsigned oldacceptance;
  unsigned oldnrtrials;

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


  // FUNCTION: add_linearpred_multBS_Block
  // TASK: updates the linear predictor if beta is replaced by the values in
  //       'fcrand' between position 'begupdate' and 'endupdate'

  void add_linearpred_multBS_Block(const unsigned a,const unsigned e,const unsigned beg,const unsigned end);

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

  void make_Kab_list(void);

  void adjust_blocksize(const unsigned & alphamin,const unsigned & alphamax);


  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_pspline_surf(void) : spline_basis_surf()
    {
    }

  // CONSTRUCTOR 1

  FULLCOND_pspline_surf(MCMCoptions * o,DISTRIBUTION * dp,
                        FULLCOND_const * fcc,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const unsigned & minb, const unsigned & maxb,
                         const double & l, const int & gs, const fieldtype & ft,
                         const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const unsigned & c=0);


  // CONSTRUCTOR 2: geosplines

  FULLCOND_pspline_surf(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                         const datamatrix & region, const MAP::map & mp,
                         const ST::string & mn, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const unsigned & minb, const unsigned & maxb,
                         const double & l, const int & gs, const fieldtype & ft,
                         const ST::string & fp, const ST::string & pres, const unsigned & c=0);


  // CONSTRUCTOR 3 varying coefficients

  FULLCOND_pspline_surf(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,const datamatrix & intact,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const unsigned & minb, const unsigned & maxb,
                         const double & l, const int & gs, const fieldtype & ft,
                         const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const unsigned & c=0);


  // CONSTRUCTOR 4: geosplines varying coefficients

  FULLCOND_pspline_surf(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,const datamatrix & intact,
                         const datamatrix & region, const MAP::map & mp, const ST::string & mn, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const unsigned & minb, const unsigned & maxb,
                         const double & l, const int & gs, const fieldtype & ft,
                         const ST::string & fp, const ST::string & pres, const unsigned & c=0);


  void init_maineffects(spline_basis * mp1,spline_basis * mp2,
                        const ST::string & pnt,const ST::string & prt);


  // COPY CONSTRUCTOR

  FULLCOND_pspline_surf(const FULLCOND_pspline_surf & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_pspline_surf & operator=(const FULLCOND_pspline_surf & fc);

  void update(void);

  void outresults(void);

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND::reset();
    }

  bool posteriormode(void)
    {
    return true;
    }

  double compute_quadform(void)
    {
    if(centertotal)
      return Ksp.compute_quadform(beta,0);
    else
      return Ksp.compute_quadform(beta_uncentered,0);
    }


  // DESTRUCTOR

  ~FULLCOND_pspline_surf() {}

  };



}   // end: namespace MCMC

#endif
