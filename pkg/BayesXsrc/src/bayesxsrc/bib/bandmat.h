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



#if !defined(SYMBANDMATRIX_INCLUDED)

#define SYMBANDMATRIX_INCLUDED

#include "statmat.h"
#include "map.h"

//------------------------------------------------------------------------------
//-------------------------- CLASS: symbandmatrix ------------------------------
//------------------------------------------------------------------------------

// USAGE OF CLASS symbandmatrix

// class 'symbandmatrix' is used to store symmetric positive definite band
// matrices and to solve linear equation systems.
// Since the matrices are symmetric only the main diagonal and the non zero
// elemens above the main diagonal are stored.


// CREATING MATRICES:

// Currently there are two constructors for creating matrices.

// CONSTRUCTOR 1:

// symbandmatrix<double> m(diagelem,upperdiagelem);

// 'diagelem' and 'upperelem'  are matrices of type statmatrix, see below.
// 'diagelem' must be a  n x 1 matrix where n is the number of rows of 'm'.
// Its elements store the diagonal elements of 'm'.
// 'upperelem' stores the non zero entries above the main diagonal of 'm'.
// It must be of size  n x b, where b is the bandsize of the bandmatrix 'm'

// EXAMPLE:

// statmatrix<double> u(100,1);
// statmatrix<double> ud(100,3);

// Insert now the elements of u and ud. This achieved for example by
// u(0,0) = 1.2;
// u(1,0) = 1.5;
// etc.
// or in a for loop.

// Create now the symbandmatrix m
// symbandmatrix<double> m(u,ud);
// This creates a 100 x 100 bandmatrix with bandsize 3


// CONSTRUCTOR 2: (for diagonal matrices)


//  symbandmatrix<double> m(de);

// de must be a n x 1 matrix of type statmatrix storing the diagonal elements
// de is constructed in the same way as for constructor 1


// SOLVING LINEAR EQUATIONS:

// Suppose we want to solve the system A x = b where A is a bandmatrix

// We first have to create the bandmatrix as shown above.

// symbandmatrix<double> A(de,ud);

// Matrices x (n x 1) and b (n x 1) must be created as well.
// They have to be of type statmatrix<double>

//  statmatrix<double> x(n,1);
//  statmatrix<double> b(n,1);

// Fill now the elements of b using for loops etc.

// We can now solve the system by:

//   A.solve(b,x,0,0);

// The solution will now be stored in x. The last two arguments denote the
// column of b and x, respectively, to use for solving the systems.
// This means that x as well as b may consist of more than one column.


// USING THE CLASS IN A GIBBS SAMPLER:

// The problem in Gaussian linear or nonparametric models is, that we
// usually have a posterior precision matrix P, say, with bandmatrix structure.
// The full conditionals are Gaussian with mean

// E(beta) = P^{-1} X' 1/sigma^2 (y-offset)

// and Covariance matrix Sigma = P^{-1}
// Here X is the design matrix. For instance, X may consist of the Basis functions
// of a B-spline (smoothing splines) evaluated at the observations.

// Additional problem: In every iteration of the sampler the precision matrix
// usually changes, because of a smoothing parameter.


// The objective is now to draw random numbers from the full conditional.

// 1. step:

// recompute the precision matrix
// To save computing time it is better not to reconstruct the matrix but
// just to overwrite the existing matrices
// Usually the precision matrix is composed of the sum of two matrices
// P = (X'X/sigma^2 + K/tau^2)
// where sigma^2 is the variance component of the observation model and
// tau^2 is the variance component of the smoothing spline or what ever.

// to update the precision we can use the functions addto and addto2,
// see below for a description.

// 2. step

// Compute the Cholesky factorisation P = L'L and solve the System
// L beta = a, where a is N(0,I) (standard normal).

// This achieved using the function solveL (see below)

// P.solveL(a,beta);

// beta is now a random number from a N(0,P^{-1}) distribution

// 3. Step

// Compute b = 1/sigma^2 X'(y-offset)

// 4. step

// Solve the system P u = b using the function

// P.solve(b,betahelp,0,0);

// Here betahelp is an auxiliary matrix with the same saize as beta

// 5. step

// set beta = beta + betahelp using

// beta.plus(beta,betahelp);

// This function adds betahelp to beta without reallocation of memory (faster)

// Then beta is distributed according to N(E(beta),P^{-1})

using std::endl;

#if defined(__BUILDING_GNU)
template <class T>
class symbandmatrix;

template <class T>
symbandmatrix<T> operator*(const T & v,const symbandmatrix<T> & bm);
#endif

template <class T>
class symbandmatrix
  {

  private:

  unsigned bands;                    // bandisze of the matrix
  unsigned dim;                      // dimension of the matrix

  bool decomposedonly;               // decomposedonly=true indicates that
                                     // the matrix is available only in
                                     // decomposed form, i.e. diagelem and
                                     // upperelem are NOT accessible

  statmatrix<T> diagelem;
  statmatrix<T> upperelem;

  bool decomposed;

  double det;                         // the derminant

  statmatrix<T> D;
  statmatrix<T> R;

  statmatrix<T> r;
  statmatrix<T> z;

  public:

  // DEFAULT CONSTRUCTOR

  symbandmatrix(void) {}

  // CONSTRUCTOR 0
  // Initializes a symmetric band matrix with dim 'd' and bandsize 'bs'

  symbandmatrix(const unsigned & d,const unsigned & bs,const double & v=0);

  // CONSTRUCTOR 1

  symbandmatrix(const statmatrix<T> & de,const statmatrix<T> & ud,
                bool decomp = false);

  // CONSTRUCTOR 2 (for diagonalmatrices)

  symbandmatrix(const statmatrix<T> & de);

  // COPY CONSTRUCTOR

  symbandmatrix(const symbandmatrix & bm);

  // OVERLOADED ASSIGNMENT OPERATORS

  const symbandmatrix & operator=(const symbandmatrix & bm);

  symbandmatrix & operator=(const T & v)
    {
    return symbandmatrix(1,0,v);
    }

  void set_decomposed(void)
    {
    decomposed = false;
    decomposedonly = false;
    }

  // FUNCTION: get_det
  // TASK: returns the logarithm of the determinant

  T get_det(void);

  // Multiplication operator
  // TASK: computes v*bm, where v is a scalar. The resulting bandmatrix is
  //       assigned to the calling matrix

  #if defined(__BUILDING_GNU)
  friend symbandmatrix operator*<>(const T & v,const symbandmatrix & bm);
  #else
  friend symbandmatrix operator*(const T & v,const symbandmatrix & bm);
  #endif

  // ACCESS TO THE i,j th element

  T operator()(const unsigned & i, const unsigned & j) const;

  // TASK: sets the element (i,j) to t
  // NOTE: (i,j) must be an element within the bands

  void set(const unsigned & i, const unsigned & j, const T & t);

  // FUNCTION: assign
  // TASK: reassigns the values of the calling matrix without realocation of
  //       memory
  // NOTE: 'de' and 'ud' must have proper size

  void assign(const statmatrix<T> & de,const statmatrix<T> & ud,
              bool decomp = false);

  // FUNCTION: dim
  // TASK returns the number of rows (and/or columns)

  const unsigned & getdim(void) const
    {
    return dim;
    }

  // FUNCTION: bandsize
  // TASK: returns the bandsize
  // NOTE: Diagonal matrices have bandsize 0, triangular matrices have size 1
  //       etc.

  const unsigned & bandsize(void) const
    {
    return bands;
    }

  // FUNCTION: print
  // TASK: writes the elements of the matrix to out
  //       writes only the diagonal elements (first column) and the nonzero
  //       elements above the diagonal (subsequent columns)

  void print(ostream & out);

  // FUNCTION: print2
  // TASK: writes the elements of the matrix to out
  //       zero elements are considered in addition

  void print2(ostream & out);

  // FUNCTION compute_quadform
  // TASK: computes x'Kx
  //       column 'c' of x ist used to compute the quadratic form

  T compute_quadform(const statmatrix<T> & x,const unsigned & c) const;

  // FUNCTION compute_quadform
  // TASK: computes x'K_[a,b]x
  //       column 'c' of x ist used to compute the quadratic form

  T compute_quadformblock(const statmatrix<T> & x,const unsigned & c,
                          const unsigned & a,const unsigned & b);

  // FUNCTION: decomp
  // TASK: computes the R'DR decomposition of the calling matrix

  void decomp(void);

  // FUNCTION: decomp_block
  // TASK: computes the R'DR decomposition of the calling matrix
  //        thecalling matrix is assumed to be a block bandmatrix

  // FUNCTION: printdecomp
  // TASK: writes the R'D'R decomposition, first column diagonal elements of D,
  //       subsequent columns: Elements of R above the diagonal (that is non-
  //       zero elements of R

  void printdecomp(ostream & out);

  // FUNCTION: solve
  // TASK: solves the equation Ax = a(,colres) and stores the result x in
  //       'res(,colres)'.
  // NOTE: 'res' must have proper size

  void solve(const statmatrix<T> & a, statmatrix<T> & res,const unsigned & cola,
             const unsigned & colres);

  void solveL(const datamatrix & z,datamatrix & res);

  // FUNCTION: inverse

  void inverse(statmatrix<T> & res);

  // FUNCTION: addto
  // TASK: computes Q = f1*X+f2*K, where f1 and f2 are scalars and assigns the
  //       result to the calling matrix
  //       X is assumed to be a diagonal matrix
  //       All matrices must have proper size

  void addto(const symbandmatrix<T> & X,const symbandmatrix<T> & K,
             const T & f1,const T & f2);

  // FUNCTION: addto
  // TASK: changes ONLY the diagonal elements of the calling matrix
  //       comutes for the diagonal: f*X+K where f is a scaler
  //       X is assumed to be a diagonal matrix
  //       All matrices must have proper size

  void addtodiag(const symbandmatrix<T> & X,const symbandmatrix<T> & K,
             const T & f);

  void addtoblock(const symbandmatrix<T> & X,const symbandmatrix<T> & K,
                  const T & f1,const T & f2,const unsigned & a,
                  const unsigned & b);

  // FUNCTION: addto2
  // TASK: computes Q = f1*X+f2*K, where f1 and f2 are scalars and assigns the
  //       result to the calling matrix
  //       X is assumed to be nondiagonal bandmatrix
  //       All matrices must have proper size
  //       IMPORTANT: bandsize of X must be LOWER or equal to the size of K

  void addto2(const symbandmatrix<T> & X,const symbandmatrix<T> & K,
              const T & f1,const T & f2);

  void addtoblock2(const symbandmatrix<T> & X1,const symbandmatrix<T> & X2,
                  const T & f1,const T & f2,const unsigned & a,
                  const unsigned & b);

  // FUNCTION: mult
  // TASK: computes A*X where A is the calling matrix and assigns the result to
  //       'res'
  //       X and res must have proper size that is they must be declared before

  void mult(const statmatrix<T> & X,statmatrix<T> & res) const;

  // FUNCTION: multBlock
  // TASK: taks the submatrix determined by rowfirst-rowlast
  //       and colfirst-collast and multiplies it with the column vektor
  //       'x'. The multiplication starts in row 'startx' of 'x'.
  //       The result is stored in 'res'

  void multBlock(const statmatrix<T> & x,statmatrix<T> & res,
                 unsigned rowfirst,unsigned colfirst, unsigned rowlast,
                 unsigned collast,unsigned startx);

  // FUNCTION: getBlock
  // TASK: returns the block from
  //       row    'rowfirst' - 'rowlast' and
  //       column 'colfirst' - 'collast'

  statmatrix<T> getBlock(const unsigned & rowfirst,const unsigned & colfirst,
                      const unsigned & rowlast,const unsigned & collast);

  // FUNCTION: getBlockasband
  // TASK: extracts the block between row 'first' - 'last' and column
  //       'first' - 'last' and stores the result in the bandmatrix 'res'
  //       If 'res' has already proper size, the elements of 'res' will not be
  //       recreated. Otherwise, res will be recreated.

  void getBlockasband(const unsigned & first,const unsigned & last,symbandmatrix<T> & res);


  // FUNCTION: getL
  // TASK: assigns the cholesky root of the calling matrix to L

  void getL(datamatrix & L);


  // FUNCTION: getdiagpointer
  // TASK: returns a pointer to the first element of diagelem

  T * getdiagpointer(void) const
    {
    return diagelem.getV();
    }

  // FUNCTION: getupperpointer
  // TASK: returns a pointer to the first element of upperelem

  T * getupperpointer(void) const
    {
    return upperelem.getV();
    }

  // DESTRUCTOR

  ~symbandmatrix() {}

  };

typedef symbandmatrix<double> bandmatdouble;


#if !defined(BANDMAT_CPP_INCLUDED)
#include"bandmat.cpp"
#endif

#endif

