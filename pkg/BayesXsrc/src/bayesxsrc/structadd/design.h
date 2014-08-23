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



#if !defined (DESIGNINCLUDED)

#define DESIGNINCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"sparsemat.h"
#include "statmat_penalty.h"
#include"Random.h"
#include"envmatrix_penalty.h"
#include"../values.h"
#include<fstream>
#include<vector>
#include<bitset>
#include"GENERAL_OPTIONS.h"
#include"FC_linear.h"
#include"clstring.h"
#include"distr.h"
#include<cmath>


namespace MCMC
{


enum ttype2 {   Rw1,
                Rw2,
                Rw3,
                Mrf,
                Hrandom,
                Grf
                };


enum effecttype2 {
                Function,
                Varcoefftotal
   };


// Centermethods
// 1. Without sampling
// a) meansimple
//    center parameters around zero
// b) integralsimple
//    integral f = 0
// 2. Within sampling
// a) meancoeff
//    center parameters around zero
// b) meanf
//    sum_i f(x_i) = 0 (over all observations)
// c) meanfd
//    sum_j f(x_(j)) = 0 (over the different observations)

enum centerm {meancoeff,nullspace,meansimple,integralsimple,
              cmeanintegral,cmeaninvvar,meanf,meanfd,meansum2,meansimplevar};

//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN ------------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN
  {

  protected:

  unsigned compute_modecategorie(void);

  void compute_meaneffectintvar(void);

  GENERAL_OPTIONS * optionsp;
  DISTR * likep;                             // Pointer to DISTR obejct

  //----------------------------------------------------------------------------

  vector< vector<double> > ZoutT;            // Nonzero Elements of Z'
  vector< vector<int> > index_ZoutT;         // Columns of nonzero elements of
                                             // Z'
  int consecutive_ZoutT;                     // -1 = not tested
                                             //  0 = not consecutive
                                             //  1 = consecutive

  void compute_Zout_transposed(void);        // Computes Z', i.e. ZoutT from
                                             // Zout

  bool check_ZoutT_consecutive(void);        // checks if non zero elements are
                                             // consecutive

  //----------------------------------------------------------------------------

  public:

  bool changingdesign;

  // Variables determined by function init_data

  bool discrete;                             // data discrete or continuous
  datamatrix data;                           // data matrix
  datamatrix intvar;                         // interaction variable for
                                             // varying coefficients
  datamatrix intvar2;                        // intvar^2 for varying coefficients

  statmatrix<unsigned> ind;                  // the category in the sorted data
                                             // vector

  statmatrix<int> index_data;                // index for sort of data
  vector<ST::string> datanames;              // names of covariates


  vector<ST::string> effectvalues;           // values of the different
                                             // covariates

  unsigned meaneffectnr;                    // position of meaneffect value
  unsigned meaneffectnr_intvar;             // position in intvar for meaneffect
                                            // value of intvar
  double meaneffectintvar;                  // value of the interaction variable
                                            // at which the meaneffect will be
                                            // computed

  //------------------------------- Zout ---------------------------------------

  datamatrix Zout;                           // Design matrix (only non null
                                             // elements for output of results
  statmatrix<int> index_Zout;                // stores the columns of the
                                             // non null elements of Zout

  vector<int> posbeg;                        // begin and end of equal covariate
                                             // values in data
  vector<int> posend;

  int consecutive;                           // -1 = not tested
                                             //  0 = not consecutive
                                             //  1 = consecutive

  bool identity;                             // true if Zout identity matrix
  bool full;                                 // true if Zout full matrix

  bool check_Zout_consecutive(void);

  // --------------------------- check intvar ----------------------------------

  // FUNCTION: compute_kernel_intvar
  // TASK: computes integral |intvar| g(intvar) dintvar where g is the epanechnikov
  //       kernel density estimator with AMISE optimal bandwidth

  double compute_kernel_intvar(bool absolute);

  double compute_epanechnikov(double & x, statmatrix<int> & intindex, vector<int> & posb, vector<int> & pose);

  // ---------------------------------------------------------------------------

  //---------------------------- Zout_derivative -------------------------------

  bool derivative;                         // true if derivatives should be
                                           // computed/estimated

  datamatrix Zout_derivative;              // Design matrix of first derivative

  bool derivative_full;                    // true if derivative design marix
                                           // is full


  void compute_f_derivative(datamatrix & beta,datamatrix & betalin,
                                  datamatrix & f_der, datamatrix & ftot);


  //----------------------------------------------------------------------------

  vector<double> ZoutTZout;
  vector<int> beg_ZoutTZout;
  vector<int> Wsump;

  vector<double> ZoutTZout_d;

  vector<int> Wsump_d;

  void compute_ZoutTZout(unsigned & i, unsigned & j);

  void compute_ZoutTZout(void);

  //----------------------------------------------------------------------------

  unsigned nrpar;                            // number of parameters

  // --------------------------- for center ------------------------------------

  double compute_sumBk(unsigned & k);
  double compute_sumBk_different(unsigned & k);

  bool center;
  centerm centermethod;


  // for nullspace centering
  datamatrix basisNull;                     // contains a basis of the null
                                            // space of the penalty K
  vector<datamatrix> basisNullt;            // contains the transposed of
                                            // basisNull

  FC_linear * FClinearp;                    // Pointer to linear effects
  int position_lin;                         // position in the designmatrix
                                            // of linear effects
  datamatrix designlinear;                  // designmatrix linear effects
  // end for nullpsace centering


  // ---------------------------------------------------------------------------

  // Variables determined by function compute_penalty

  envmatdouble K;                            // Penalty Matrix
  double rankK;

  // ---------------------------------------------------------------------------

  // Variables determined by function compute_precision

  envmatdouble precision;                    // precision matrix
  bool precisiondeclared;                    // true if precision is already
                                             // defined

  // ---------------------------------------------------------------------------

  // Variables determined by function  compute_XtransposedWX_XtransposedWres
  // and compute_XtransposedWres

  datamatrix Wsum;

  envmatdouble XWX;                          // X'WX

  datamatrix XWres;                          // X'W(y-eta)
  datamatrix * XWres_p;                      // Pointer to the current
                                             // XWres object

  // ---------------------------------------------------------------------------

  ttype2 type;                                // Term type


  //--------------------- for orthogonal transformation ------------------------

  datamatrix s;                             // contains eigenvalues of
                                             // R^-1 K R-T
  datamatrix QtRinv;
  datamatrix RtinvQ;

  datamatrix u;

  //----------------------------------------------------------------------------

  bool errors;                // true if an error has occured during
                              // initialization

  //----------------------- CONSTRUCTORS, DESTRUCTOR ---------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN(void);

  // CONSTRUCTOR

  DESIGN(GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl);

  // COPY CONSTRUCTOR

  DESIGN(const DESIGN & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN & operator=(const DESIGN & m);

  //----------------------------------------------------------------------------

  // FUNCTION: compute_f
  // TASK: compute Zout*beta, i.e. the estimated/current function evaluated at
  //       the different observations in data

  void compute_f(datamatrix & beta,datamatrix & betalin,
                       datamatrix & f, datamatrix & ftot);

  // FUNCTION: compute_effect
  // TASK: computes the effect vector

  void compute_effect(datamatrix & effect,datamatrix & f,
                      effecttype2 et = Function);

  void set_intvar(datamatrix & iv, double add=0);

  // FUNCTION: update_linpred
  // TASK: updates the predictor based on the current function f

  void update_linpred(datamatrix & f);

  bool update_linpred_save(datamatrix & f);

  // FUNCTION: compute_partres
  // TASK: computes

  void compute_partres(datamatrix & res,datamatrix & f,bool cwsum=false);

  double compute_ZtZ(unsigned & i, unsigned & j);


  void outdesign(ofstream & out);

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
  virtual void compute_penalty2(const datamatrix & par);

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
                                bool computemeaneffect,double meaneffectconstant);


  virtual void compute_precision(double l);

  // FUNCTION: read_options
  // TASK: reads options and initializes varnames stored in datanames

  virtual void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  virtual void outoptions(GENERAL_OPTIONS * op);

  virtual void outbasis_R(ofstream & out);

  virtual void compute_orthogonaldecomp(void);

  // --------------------- END: VIRTUAL FUNCTIONS ------------------------------

  // DESTRUCTOR

  ~DESIGN() {}


  };






} // end: namespace MCMC

#endif


