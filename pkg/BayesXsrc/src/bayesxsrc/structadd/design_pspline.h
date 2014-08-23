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



#if !defined (DESIGNpsplineINCLUDED)

#define DESIGNpsplineINCLUDED

#include"../export_type.h"
#include<deque>
#include<design.h>

using std::deque;

namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_pspline ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_pspline : public DESIGN
  {

  protected:

  bool ccov;                    // true if covariate should be centered

  bool multeffect;

  deque<double> knot;          // Vektor der Knoten (sichtbare und unsichtbare)

  vector<double> weightK;      // weights to compute penalty matrix K

  // FUNCTION: bspline
  // TASK: computes B-splines at position x

  datamatrix bspline(const double & x);

  // FUNCTION: bspline
  // TASK: computes B-splines first derivative at position x

  datamatrix bspline_derivative(const double & x);

  // FUNCTION: make_Bspline
  // TASK: computes knot, Zout and index_Zout

  void make_Bspline(void);

  // FUNCTION: make_Bspline_derivative
  // TASK: computes Bsplines derivative design matrix Zout_derivative

  void make_Bspline_derivative(void);

  void compute_betaweight(datamatrix & betaweight);

  public:

  long nrknots;                     // Anzahl der (sichtbaren) Knoten
  long degree;                      // Grad des Splines
  long difforder;                   // Differenzenordnung (1,2,3)
  double round;
  double binning;
  double minBS;
  double maxBS;


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_pspline(void);

  // CONSTRUCTOR

  DESIGN_pspline(datamatrix & dm, datamatrix & iv,
             GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl,
             vector<ST::string> & op,
             vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  DESIGN_pspline(const DESIGN_pspline & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_pspline & operator=(const DESIGN_pspline & m);

  // virtual functions

  void compute_penalty(void);

  void compute_penalty2(const datamatrix & pen);

  void compute_basisNull(void);

  void compute_precision(double l);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  void outoptions(GENERAL_OPTIONS * op);

  void outbasis_R(ofstream & out);


// DESTRUCTOR

  ~DESIGN_pspline() {}

  };


} // end: namespace MCMC

#endif


