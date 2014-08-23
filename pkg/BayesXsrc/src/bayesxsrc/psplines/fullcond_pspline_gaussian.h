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



#if !defined (MCMCpsplinegaussian_INCLUDED)
#define MCMCpsplinegaussian_INCLUDED

#include"../export_type.h"
#include<deque>
#include "mcmc_nonpbasis.h"
#include "spline_basis.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_pspline_gaussian ---------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_pspline_gaussian : public spline_basis
  {


  protected:

  bool samplecentered;           // Samplen unter der Nebenbedingung, dass die Zentrierungskonstante = 0 ist.
  bool diagtransform;            // Tranformation, so dass 'prec_env' eine Diagonalmatrix ist
  bool hierarchical;             // hierarchical centering

  // für hierachical centering

  double lineff;                 // linearer Anteil
  double lineffsum;              // Summe der samples für den linearen Anteil
  vector<double> lineffsamples;  // Samples für linearen Anteil
  datamatrix gamma;              // Hilfsmatrix: lineff*X*gamma ergibt Gerade mit Steigung lineff


  // FUNCTION: update_isotonic
  // TASK: updates beta (monotonic regression)

  void update_isotonic(void);

  // FUNCTION: update_diagtransform
  // TASK: updates beta (if diagtransform == true)

  void update_diagtransform(void);


  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_pspline_gaussian(void) : spline_basis()
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

  FULLCOND_pspline_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                          FULLCOND_const * fcc, const datamatrix & d,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                         const ST::string & fp, const ST::string & pres, const bool & deriv,
                         const double & l, const int & gs, const bool & diag,
                         const double & lk, const double & uk, const double & lg,
                         const double & ug, const unsigned & c=0);

  // CONSTRUCTOR 2  (for  varying coefficients term)
  // effmod: values of the effect modifier
  // intact: values of the interaction variable

  FULLCOND_pspline_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,
                         const datamatrix & effmod, const datamatrix & intact,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                         const ST::string & fp, const ST::string & pres, const bool & deriv,
                         const double & l, const int & gs, bool  ce, const unsigned & c=0);


  // COPY CONSTRUCTOR

  FULLCOND_pspline_gaussian(const FULLCOND_pspline_gaussian & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_pspline_gaussian & operator=(const FULLCOND_pspline_gaussian & fc);

  void update(void);

  void outresults(void);

  // FUNCTION: updateK
  // TASK: updates the penalty matrix K with q (K = D'diag(q)D)

  void updateK(const datamatrix & q)
    {
    FULLCOND_nonp_basis::updateK(q);
    }

  void compute_u(datamatrix & u)
    {
    FULLCOND_nonp_basis::compute_u(u);
    }

  double compute_quadform(void);

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND_nonp_basis::reset();
    }

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  // FUNCTION: predict
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred);

  // FUNCTION: write_contour
  // TASK: writes the mean of the full conditional, 1/scale, 1/sigma2 etc. to a
  // temporary ASCII file. (Necessary to compute countour probabilities)

  void write_contour();

  // FUNCTION: compute_contourprob
  // TASK: computes the contour probabilities for beta

  void compute_contourprob(void);

  // FUNCTION: compute_contourprob
  // TASK: computes the contour probabilities for differences of order 'diff'

  void compute_contourprob(const int & diff);

  // FUNCTION: compute_pseudocontourprob
  // TASK: computes the pseudo contour probabilities for differences of order 'diff'

  void compute_pseudocontourprob(const int & diff);

  // DESTRUCTOR

  ~FULLCOND_pspline_gaussian() {}

  };


} // end: namespace MCMC

#endif

