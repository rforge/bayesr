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



#if defined (__BUILDING_THE_DLL)
#define __EXPORT_TYPE __export
#elif defined (__BUILDING_GNU)
#define __EXPORT_TYPE  __attribute__((dllexport))
#else
#define __EXPORT_TYPE __import
#endif

#if !defined (DESIGNpsplinegroup_INCLUDED)

#define DESIGNpsplinegroup_INCLUDED

#include<deque>
#include<design.h>

using std::deque;

namespace MCMC
{


//------------------------------------------------------------------------------
//------------------------- CLASS: DESIGN_pspline_group ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_pspline_group  : public DESIGN
  {

  protected:

  public:

  //----------------------- CONSTRUCTORS, DESTRUCTOR ---------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_pspline_group(void);

  // CONSTRUCTOR

  DESIGN_pspline_group(DISTR * dp,FC_linear * fcl);

  // COPY CONSTRUCTOR

  DESIGN_pspline_group(const DESIGN_pspline_group & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_pspline_group & operator=(const DESIGN_pspline_group & m);

  // ------------------------- VIRTUAL FUNCTIONS -------------------------------

  // FUNCTION: init_data
  // TASK: sorts the data,
  //       creates data, intvar, intvar2
  //       initializes index_data, data, intvar, intvar2, ind
  //       posbeg, posend, effectvalues
  //       meaneffectnr,
  //       effectvalues

                                                                 virtual void init_data(const datamatrix & dm, const datamatrix & iv);

  // FUNCTION: compute_penalty
  // TASK: computes the penalty matrix and determines rankK

  virtual void compute_penalty(void);

  virtual double penalty_compute_quadform(datamatrix & beta);

  // FUNCTION: compute_basisNull
  // TASK: computes the basis of the null space of the penalty matrix

  virtual void compute_basisNull(void);

  // FUNCTION: computes XWres
  // TASK: computes XWres, res is the partial residual

  virtual void compute_XtransposedWres(datamatrix & partres, double l);

  // FUNCTION: compute_XtransposedWX_XtransposedWres
  // TASK: computes XWX and XWres, res is the partial residual

  virtual void compute_XtransposedWX(void);

  // FUNCTION: compute_meaneffect
  // TASK: computes the meaneffect of a particular model term

  virtual void compute_meaneffect(DISTR * level1_likep,double & meaneffect,
                                datamatrix & beta,datamatrix & meaneffectbeta,
                                bool computemeaneffect);


  virtual void compute_precision(double l);

  // FUNCTION: read_options
  // TASK: reads options and initializes varnames stored in datanames

  virtual void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  virtual void outoptions(GENERAL_OPTIONS * op);

  virtual void compute_orthogonaldecomp(void);

  // --------------------- END: VIRTUAL FUNCTIONS ------------------------------

  // DESTRUCTOR

  ~DESIGN_pspline_group() {}


  };






} // end: namespace MCMC

#endif


