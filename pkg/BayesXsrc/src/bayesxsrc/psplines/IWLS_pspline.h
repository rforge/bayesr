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



#ifndef IWLS_psplineH
#define IWLS_psplineH

#include"../export_type.h"
#include "mcmc.h"
#include "fullcond.h"
#include "time.h"
#include <deque>
#include "sparsemat.h"
#include "mcmc_nonpbasis.h"
#include "spline_basis.h"
#include "fullcond_nonp_gaussian.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//---------------------------- class: IWLS_pspline -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE IWLS_pspline : public spline_basis
  {

  protected:

  updatetype utype;                         // iwls || iwlsmode || hyperblock || hyperblockmode
                                            // increasing || decreasing || diagtransform

// für gemeinsames updaten von beta und sigma2 (hyperblock)

  double a_invgamma;                        // Parameter a der IG(a,b) für sigma2
  double b_invgamma;                        // Parameter b der IG(a,b) für sigma2
  double kappa;                             // 1/sigma2
  double kappaprop;                         // vorgeschlagenes kappa
  double kappamode;                         // kappa für hyperblock
  double kappamean;                         // Hilfvariable für Startwert von kappamode

  bool diagtransform;                       // Tranformation, so dass 'prec_env' eine Diagonalmatrix ist

  unsigned updateW;                         // jede wievielte Iteration soll IWLS-Gewicht W neu berechnet werden?

  void create_iwls(void);

  void update_IWLS(void);                   // update nach IWLS

  void update_IWLS_mode(void);              // update nach IWLS basierend auf dem posterior mode

  void update_IWLS_hyperblock(void);        // gemeinsames updaten (IWLS) von beta und kappa

  void update_IWLS_hyperblock_mode(void);   // gemeinsames updaten (IWLS basierend auf dem posterior mode) von beta und kappa

  void update_isotonic(void);               // update bei monotoner Regression

  void update_diagtransform(void);          // update bei diagtransform == true


  public:

  // DEFAULT CONSTRUCTOR

  IWLS_pspline(void) : spline_basis()
    {
    }

  // CONSTRUCTOR 1

  IWLS_pspline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const fieldtype & ft, const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const double & lk,
                    const double & uk, const double & lg, const double & ug,
                    const unsigned & c);

  // CONSTRUCTOR 2 (für variierende Koeffizienten)

  IWLS_pspline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & effmod,const datamatrix & intact,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const fieldtype & ft, const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c);

  // CONSTRUCTOR 3 (für Cox)

  IWLS_pspline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const fieldtype & ft, const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c);

  // COPY CONSTRUCTOR

  IWLS_pspline(const IWLS_pspline & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const IWLS_pspline & operator=(const IWLS_pspline & fc);

  void update(void);

  void outresults(void);

  void outoptions(void);

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    f = 100.0;
    oldacceptance = 0;
    oldnrtrials = 0;
    FULLCOND_nonp_basis::reset();
    }

  // FUNCTION: predict
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred);

  // FUNCTION: compute_quadform
  // TASK: returns beta(.,v)' K beta(.,v) where K is the penalty matrix

  double compute_quadform(void);


  // DESTRUCTOR

  ~IWLS_pspline() {}

  };


}   // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
