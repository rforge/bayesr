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



#ifndef IWLS_baselineH
#define IWLS_baselineH

#include"../export_type.h"
#include"mcmc.h"
#include"fullcond.h"
#include "time.h"
#include <deque>
#include"sparsemat.h"
#include"mcmc_nonpbasis.h"
#include"IWLS_pspline.h"
#include"fullcond_nonp_gaussian.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//---------------------------- class: IWLS_baseline -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE IWLS_baseline : public IWLS_pspline
  {

  protected:

  bool begin0;
  datamatrix int_knots;
  datamatrix int_D;
  MCMC::bsplinemat testmat;
  vector<MCMC::bsplinemat> gaussmat;
  vector<IWLS_baseline*> baselinep;
  datamatrix zi;
  unsigned gauss_n;
  datamatrix coeff;
  datamatrix z_vc;
  datamatrix zi_ges;
  datamatrix beg_i;
  datamatrix A;
  datamatrix distance;
  datamatrix interval;
  datamatrix AWA;
  datamatrix Wbase;
  statmatrix<int> zi_index;
  statmatrix<int> ges_index;
  datamatrix spline_ges;
  datamatrix spline_ges2;
  datamatrix spline_zi;
  datamatrix spline_zi2;
  datamatrix gaussspline;
  datamatrix int_ti_help;
  datamatrix int_deriv;
  datamatrix int_H;
  datamatrix response_help;
  datamatrix Xdelta;
  datamatrix Adelta;
  datamatrix Eins;
  datamatrix score;
  datamatrix deltaexact;
  datamatrix An;
  datamatrix DeltaN;
  datamatrix cov_cp;
  bool vc_dummy1;

  void update_IWLS(void);

  void update_IWLS_mode(void);


  public:

  // DEFAULT CONSTRUCTOR

  IWLS_baseline(void) : IWLS_pspline()
    {
    }

  // CONSTRUCTOR 1

  IWLS_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const fieldtype & ft, const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c, const datamatrix & anfang);

  // CONSTRUCTOR 2 (für variierende Koeffizienten)

  IWLS_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & effmod,const datamatrix & intact,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const fieldtype & ft, const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c, const datamatrix & anfang);

  // COPY CONSTRUCTOR

  IWLS_baseline(const IWLS_baseline & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const IWLS_baseline & operator=(const IWLS_baseline & fc);

  void update(void);

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

  void set_fcconst(FULLCOND_const * fcc)
    {
    fcconst = fcc;
    }

//-------- für baseline --------------------------------------------------------

  void compute_int_ti(const datamatrix & b);
  void compute_int_ti_linear(const double & b);
  void compute_int_ti(unsigned beg);
  void compute_int_ti_vc_di0(const vector<double *>,const vector<double *>);
  void compute_int_ti_vc_di(const int,const vector<double *>,const vector<double *>);
  void compute_int_gauss(void);
  void compute_int_gauss_DIC(void);
  void update_baseline(void);
  void compute_int_ti_mean(void);
  void compute_int_deriv(const datamatrix & b);
  void compute_int_H(const datamatrix & b);

  void set_baselinep(vector<IWLS_baseline*> bp)
    {
    baselinep=bp;
    }

  double * get_int_D(void)
    {
    return int_D.getV();
    }

  double * get_z_vc(void)
    {
    return z_vc.getV();
    }

  datamatrix get_z_vc_np(void)
    {
    return z_vc;
    }

  double * get_spline_zi(void)
    {
    multBS(spline_zi,beta);
    return spline_zi.getV();
    }

  double * get_spline_ges(void)
    {
    testmat.mult_index(spline_ges,beta);
    return spline_ges.getV();
    }

  double * get_spline_ges2(void)
    {
    testmat.mult_index(spline_ges2,beta);
    return spline_ges2.getV();
    }

  double * get_spline_ges_mean(void)
    {
    testmat.mult_index(spline_ges,betamean);
    return spline_ges.getV();
    }

  double * get_gaussspline(void);


  double * get_gaussspline_mean(void);


  double * get_betamean(void)
    {
    return betamean.getV();
    }

  double * get_beta(void)
    {
    return beta.getV();
    }

  double * get_spline_zi_mean(void)
    {
    multBS(spline_zi,betamean);
    return spline_zi.getV();
    }

  void compute_AWA(void);
  void compute_Wbase(void);
  void compute_score(void);

//---------- ENDE: für baseline ------------------------------------------------

  // DESTRUCTOR

  ~IWLS_baseline() {}

  };


}   // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
