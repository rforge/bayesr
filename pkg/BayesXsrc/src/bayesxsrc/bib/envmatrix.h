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


#if !defined(ENVMATRIX_INCLUDED)

#define ENVMATRIX_INCLUDED

#include "statmat.h"
#include "map.h"
#include "bandmat.h"

using std::cout;
using std::endl;

//------------------------------------------------------------------------------
//------------------------------ CLASS: envmatrix ------------------------------
//------------------------------------------------------------------------------

const double sqrtmin = 0.0000000000000000001;
const double sqrtmax = 100000000000000000000.0;
const double logmin =  0.0000000000000000001;
const double logmax =  100000000000000000000.0;

template<class T>
class envmatrix
  {

  private:

  vector<T> diag;                    //vector containing diagonal elements

  vector<T> env;                     //vector containing nonzero-entries

  vector<T> ldiag;                   //vector containing the diagonal elements
                                     //of the lower triangular factor of the
                                     //matrix if it is cholesky decomposed and
                                     //the elements of the matrix D if it is
                                     //rational cholesky decomposed.

  vector<T> lenv;                    //vector containing the envelope elements
                                     //of the lower triangular factor of the
                                     //matrix if it is cholesky decomposed and
                                     //the elements of the envelope of U without
                                     //if it is rational cholesky decomposed.

  vector<unsigned> xenv;             //vector of length dim+1 containing the
                                     //envelope-structure and the size of the
                                     //enevelope stored in xenv[dim]

  unsigned dim;                      //dimension of the matrix

  bool decomposed;                   //wether the matrix is cholesky
                                     //decomposed or not

  bool rational_decomposed;          //wether the matrix is rational cholesky
                                     //decomposed or not

  int bandwidth;                     //the bandwidth of the matrix. 0 indicates
                                     //a diagonal matrix, -1 indicates a matrix
                                     //with a more general envelope structure

  public:


//-----------------------------------------------------------------------------
//----------------------------- Constructor------------------------------------
//-----------------------------------------------------------------------------


  // DEFAULT CONSTRUCTOR

  envmatrix(void) { dim = 0; }


  // CONSTRUCTOR1
  // Initializes an enevelope-matrix with dimension d, envelope-structure xenv
  // and value v for all nonzero entries


  envmatrix(const vector<unsigned> & xe, const T v, const unsigned d);

  // CONSTRUCTOR2
  // Initializes an envelope-matrix and stores the data provided in v in the
  // envelope of the matrix according to the storing scheme provided in xe
  // and the data provided in d in the diag-vector of the matrix
  // Therefore d.size()+1==xe.size() and v.size()==xe[d.size()]


  envmatrix(const vector<T> & v, const vector<T> & d,
            const vector<unsigned> & xe);


  // CONSTRUCTOR3
  // Initializes an envelope-matrix with the values given in the statmatrix<T> X
  // It is assumed that the matrix is symmetric positive definite and the
  // envelope-structure is exploited in the storing scheme

  envmatrix(const statmatrix<T> & X, const double epsilon=0.0);


  // CONSTRUCTOR4
  // Initializes a diagonal matrix with value v for all diagonal elements
  // and dimension d

  envmatrix(const T &v, const unsigned &d);


  // CONSTRUCTOR5
  // Initializes a diagonal matrix with values v and dimension d

  envmatrix(const vector<T> &v, const unsigned &d);


  // CONSTRUCTOR6
  // Initializes a band matrix with bandwidth bw, the values in d as diagonal
  // elements and the values in v as envelope elements. The envelope elements
  // in v are assumed to have the structure specified in xe. Therefore
  // v.size()==bw*dim-(bw+1)*bw/2 where dim is the dimension of the matrix.

  envmatrix(const vector<T> & v, const vector<T> & d,
            const vector<unsigned> & xe, const int & bw);


  // CONSTRUCTOR7
  // Initializes an envelope-matrix with the values given in the symbandmatrix X

  envmatrix(const symbandmatrix<T> & X);

  // CONSTRUCTOR8
  // Initializes an envelope-matrix with bandwidth bw, value v for all
  // nonzero elements and dimension d.

  envmatrix(const T &v, const unsigned &d, const unsigned bw);


  // Copy CONSTRUCTOR

  envmatrix(const envmatrix & em);

  // OVERLOADED ASSIGNMENT OPERATOR

  const envmatrix & operator=(const envmatrix & em);

  //

  T operator()(const unsigned & i, const unsigned & j) const;

  T get(const unsigned & i, const unsigned & j) const;

  statmatrix<T> get(void) const;

  // DESTRUCTOR

  ~envmatrix() {}

//-----------------------------------------------------------------------------
//--- Functions that decompose a matrix, solve systems of linear equations ----
//--------------- or compute the envelope of the inverse ----------------------
//-----------------------------------------------------------------------------


  // FUNCTION: decomp
  // TASK: Computes the cholesky decomposition and stores it in lenv and ldiag

  void decomp();

  bool decomp_save();

  // FUNCTION: decomp2
  // TASK: Computes the cholesky decomposition and stores it in lenv and ldiag
  // NOTE: The function starts in row 'start' of the calling matrix and assumes
  //       that all element up to row 'start-1' are already available

  void decomp2(unsigned start=0);

  // FUNCTION: decomp_rational
  // TASK: Computes the rational cholesky decomposition U'D^-1U of the calling
  //       matrix and stores U in lenv and D in ldiag

  void decomp_rational();

  // FUNCTION: solve
  // TASK: solves the equation Ax=b where A is the calling matrix
  //       and stores the result in the b.
  //       b may be either a vector or a datamatrix

  void solve(vector<T> &b);

  void solve(datamatrix &b);


  // FUNCTION: solve
  // TASK: solves the equation Ax=b where A is the calling matrix
  //       and stores the result in the res.

  void solve(const datamatrix &b, datamatrix &res);


  // FUNCTION: solve
  // TASK:     solves the equation Ax=b where A is the calling matrix adds
  //           bhelp and stores the result in res.

  void solve(const datamatrix &b, const datamatrix &bhelp, datamatrix &res);


  // FUNCTION: solveL
  // TASK: solves the equation Lx = b, where L is a lower triangular matrix
  //       stored in the envelope-format and stores the result in b.
  //       b may be either a vector or a datamatrix.
  //       If res is specified in addition, the result is stored in res.

  void solveL(vector<T> & b);

  void solveL(datamatrix & b);

  void solveL(const datamatrix & b, datamatrix & res);


  // FUNCTION: solveU
  // TASK: solves the equation Ux = b, where U is an upper triangular matrix
  //       and stores the result in b.
  //       b may be either a vector or a datamatrix.
  //       U must be stored in a cloumn-enevlope-format.
  //       Equivalently the transpose of U may be stored in the usual
  //       row-enevlope-format.

  void solveU(vector<T> &b);

  void solveU(datamatrix &b);


  // FUNCTION: solveU
  // TASK: solves the equation Ux = b, where U is an upper triangular matrix,
  //       adds bhelp and stores the result in b.

  void solveU(datamatrix &b, const datamatrix &bhelp);

  // FUNCTION: inverse_envelope
  // TASK: computes the envelope of the inverse of the calling matrix and stores
  //       it in inv. inv is assumed to have the same envelope structure as the
  //       calling matrix.

  void inverse_envelope(envmatrix<T> & inv);

//------------------------------------------------------------------------------
//---- Functions to assess elements or characteristics of the matrix------------
//------------------------------------------------------------------------------


  // FUNCTION: getL
  // TASK:     returns the (i,j) element of the cholesky-factor of the calling
  //           matrix

  T getL(const unsigned & i, const unsigned & j) const;


  // FUNCTION: getL
  // TASK:     returns the cholesky-factor of the calling
  //           matrix  in statmatrix format

  statmatrix<T> getL(void) const;


  // FUNCTION: getBandwidth
  // TASK:     returns the bandwidth of the matrix. -1 indicates a matrix
  //           without banded structure

  int getBandwidth(void) const;


  // FUNCTION: getDim
  // TASK: returns the dimension of the matrix

  unsigned getDim(void) const;


  //FUNCTION: getXenv
  //TASK: returns the value of xenv(i)

  unsigned getXenv(const unsigned & i) const;

  // FUNCTION: getEnv
  // TASK: returns the value of env(i)

  T & getEnv(const unsigned & i);

  // FUNCTION: getDiag
  // TASK: returns diag, the vector with the diagonal elements
  //       of the calling matrix

  vector<T> getDiag();


  // FUNCTION: getDiag
  // TASK: returns the ith element of diag, the vector with the diagonal
  //       elements of the calling matrix

  T getDiag(unsigned i);


  // FUNCTION: getEnv
  // TASK: returns env, the envelope of the calling matrix

  vector<T> getEnv();

  // FUNCTION: getXenv
  // TASK: returns xenv, the vector containing the indices of the envelope
  //       of the calling matrix

  vector<unsigned> getXenv() const;

  // FUNCTION: getLogDet
  // TASK: returns the logarithm of the determinant of the calling matrix

  T getLogDet();

  T getLogDet_save(bool error);

  // FUNCTION: traceOfProduct
  // TASK: returns the trace of the product A*B, where A is the calling matrix.

  T traceOfProduct(envmatrix<T> & B);

  // FUNCTION: computeMaxXenv
  // TASK: Computes an envelope structure that combines the envelopes of the
  //       calling matrix and B

  vector<unsigned> computeMaxXenv(const envmatrix<T> & B);


//------------------------------------------------------------------------------
//---------- Functions to get pointers to elements of the matrix----------------
//------------------------------------------------------------------------------

  #if defined(__BUILDING_GNU)
  typename
  #endif
  vector<T>::iterator getDiagIterator();

  #if defined(__BUILDING_GNU)
  typename
  #endif
  vector<T>::iterator getEnvIterator();

  vector<unsigned>::iterator getXenvIterator();


//------------------------------------------------------------------------------
//------------- Functions for changing elements of the matrix-------------------
//------------------------------------------------------------------------------


  // FUNCTION setDiag
  // TASK: sets the ith diagonal-element to t

  void setDiag(const unsigned & i, const T & t);

  // FUNCTION set
  // TASK: sets the element (i,j) to t
  // NOTE: (i,j) must be an element of the envelope

  void set(const unsigned & i, const unsigned & j, const T & t);

  // FUNCTION setDecomposed
  // TASK: changes decomposed to the specified value. For decomposed=true the
  //       elements of ldiag and lenv are treated as the cholesky factor of the
  //       matrix

  void setDecomposed(const bool &t);

  // FUNCTION setRational_decomposed
  // TASK: changes rational_decomposed to the specified value. For
  //       rational_decomposed=true the elements of ldiag and lenv are treated
  //       as the rational cholesky factorisation of the matrix

  void setRational_decomposed(const bool &t);

//------------------------------------------------------------------------------
//-------------------- Functions for computing new matrices---------------------
//------------------------------------------------------------------------------


  // FUNCTION: addtodiag
  // TASK: computes f1*X + f2*K and assigns it to the calling matrix. X is
  //       assumed to be a diagonal matrix and K must have the same
  //       envelope-structure as the calling matrix

  void addtodiag(envmatrix &X, envmatrix &K, const T &f1,
                 const T &f2);

  // FUNCTION: addto
  // TASK: computes f1*X + f2*K and assigns it to the calling matrix.
  //       The calling matrix is assumed to have the maximum envelope of X and K
  //       which may be computed with the function getMaxXenv()

  void addto(envmatrix &X, envmatrix &K, const T &f1,
                 const T &f2);

  // FUNCTION compute_quadform
  // TASK: computes x'Kx
  //       column 'c' of x ist used to compute the quadratic form

  T compute_quadform(const statmatrix<T> & x,const unsigned & c);

  // FUNCTION compute_sumfabsdiff
  // TASK: computes sum_{i~j} -m_ij*|x_i-x_j|, where m_ij is the elemen (i,j) of the calling matrix
  //       column 'c' of x is used

  T compute_sumfabsdiff(const statmatrix<T> & x,const unsigned & c);

  // FUNCTION compute_quadformblock
  // TASK: computes x[a:b]'K[a:b,a:b]x[a:b]
  //       column 'c' of x ist used to compute the quadratic form.
  //       The calling matrix must be a band matrix.

  T compute_quadformblock(const statmatrix<T> & x,const unsigned & c,
                          const unsigned & a, const unsigned &b);

//------------------------------------------------------------------------------
//----------------- Functions for printing matrices-----------------------------
//------------------------------------------------------------------------------


  // FUNCTION: print1
  // TASK: writes the elements of the matrix to out
  //       writes only the diagonal elements and the nonzero
  //       elements under the diagonal row by row

  void print1(ostream & out = cout);


  // FUNCTION: print2
  // TASK: writes the elements of the matrix to out
  //       zero elements are considered in addition

  void print2(ostream & out = cout);


  // FUNCTION: print3
  // TASK: writes out the elements of the envelope (1st line),
  //       the elements of the envelope`s index-vector (2nd line)
  //       and the elements of the diagonal (3rd line).
  //       If the matrix is decomposed:
  //       elements of the envelope of the lower triangular factor (4th line)
  //       elements of the diagonal of the lower triangular factor (5th line)

  void print3(ostream & out = cout);


  // FUNCTION: print4
  // TASK: writes the complete matrix to out.

  void print4(ostream & out = cout);


  // FUNCTION: print4L
  // TASK: writes the complete lower triangular factor to out (if the matrix
  //       is decomposed)

  void print4L(ostream & out = cout);


  };

  typedef envmatrix<double> envmatdouble;
  typedef long double ldouble;

#if !defined(ENVMATRIX_CPP_INCLUDED)
#include"envmatrix.cpp"
#endif

#endif





