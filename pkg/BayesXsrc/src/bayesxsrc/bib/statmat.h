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



#if !defined(STATMATRIX_INCLUDED)

#define STATMATRIX_INCLUDED

#include "tmatrix.h"
#include <math.h>
#include <vector>
#include "clstring.h"

using std::vector;

//------------------------------------------------------------------------------
//---------------------------- CLASS: statmatrix -------------------------------
//------------------------------------------------------------------------------

#if !defined(__BUILDING_GNU)
class SparseMatrix;

// AENDERUNG (Eva)
class adja;
#endif

template <class T>
class statmatrix : public Matrix<T>

  {

#if defined(__BUILDING_GNU)
class SparseMatrix;

// AENDERUNG (Eva)
class adja;
#endif

  public:

  // ---------------------------- CONSTRUCTORS ---------------------------------

  // DEFAULT CONSTRUCTOR

  statmatrix(void) : Matrix<T> () {}

  // CONSTRUCTOR 1

  statmatrix(unsigned rows, unsigned cols = 1) : Matrix<T> (rows,cols) {}

  // CONSTRUCTOR 2

  statmatrix(unsigned rows, unsigned cols, const T & init)
				 : Matrix<T> (rows, cols, init) {}

  // CONSTRUCTOR 3

  statmatrix(const SparseMatrix & m);

  // CONSTRUCTOR 4

  statmatrix(const vector<T> & v);

  // COPY CONSTRUCTORS

  statmatrix(const Matrix<T> & m) : Matrix<T> (m) {}
  statmatrix(const statmatrix & s) : Matrix<T> (s) {}

  // OVERLOADED ASSIGNMENT OPERATORS

  const statmatrix & operator=(const SparseMatrix & m);

  // AENDERUNG (Eva)
  // OVERLOADED ASSIGNMENT OPERATOR
  statmatrix   operator= (const adja & a);


  // ------------------------- PUBLIC FUNCTIONS --------------------------------

  // FUNCTION: solveroot
  // TASK: solves Ax=b where A is the calling matrix, A is supposed to be a
  //       choleksy root

  void solveroot_t(const statmatrix & b,statmatrix & x);


  void solveroot(const statmatrix & b,statmatrix & help,
                              statmatrix & x);

  // FUNCTION: assign
  // TASK: assigns the elements of A to the elements of the calling matrix
  //        faster than B = A or B.putBlock(A,0,0,B.rows,B.cols())

  void assign(const statmatrix & A);

  // FUNCTION: plus
  // TASK: assigns A+B to the calling matrix
  //       faster than C = A+B

  void plus(const statmatrix & A,const statmatrix & B);

  // FUNCTION: plus
  // TASK: adds A to the calling matrix

  void plus(const statmatrix & A);

  // FUNCTION: minus
  // TASK: assigns A-B to the calling matrix
  //       faster than C = A-B

  void minus(const statmatrix & A,const statmatrix & B);

  // FUNCTION: minus
  // TASK: substracts the colA th column of A from the colB th column B and
  //       assigns the result to the calling matrix (must be column vector)

  void minus(const statmatrix & A,const statmatrix & B,const unsigned & colA,
             const unsigned & colB);

  // FUNCTION: mult
  // TASK: assigns A*B to the calling matrix

  void mult(const statmatrix & A,const statmatrix & B);

  // FUNCTION: mult_scalar
  // TASK: multiplies b*A and assigns the result to the calling matrix

  void mult_scalar(const statmatrix & A, const T & b);

  // FUNCTION: addmult
  // TASK: computes A*B and adds the result to the calling matrix

  void addmult(const statmatrix & A,const statmatrix & B);

  // FUNCTION: addmult
  // TASK: computes A*B and adds the result to the calling matrix
  //       A is assumed to be symmetric with elements stored in the lower
  //       triangular

  void addmultsym(const statmatrix & A,const statmatrix & B);

  // FUNCTION: inverse
  // TASK: computes the inverse of the calling matrix

  statmatrix<T> inverse(void);

  // FUNCTION: multdiagback
  // TASK: multiplies the calling matrix with the diagonal matrix D (from the
  //       right). The diagonal elements of D are stored in d

  void multdiagback(const statmatrix & d);

  // FUNCTION: multdiagfront
  // TASK: multiplies the calling matrix with the diagonal matrix D (from the
  //       left). The diagonal elements of D are stored in d

  void multdiagfront(const statmatrix & d);

  // FUNCTION: multdiagfront
  // TASK: multiplies A with the diagonal matrix D (from the left) and assigns
  //       the result to the calling matrix. The diagonal elements of D are
  //       stored in d

  void multdiagfront(const statmatrix & A, const statmatrix & d);

  // FUNCTION: addtodiag
  // TASK: Adds the elements in d to the diagonal elements of the calling matrix
  //       The first element of d is added to A(first,first), the last element
  //       of d is added to A(last-1,last-1)
  //       The calling matrix is assumed to have rows==cols

  void addtodiag(const statmatrix & d, unsigned first, unsigned last);

  // FUNCTION: subfromdiag
  // TASK: Subcracts the elements in d to the diagonal elements of the calling
  //       matrix. The first element of d is subtracted from A(first,first),
  //       the last element of d is subtracted from A(last-1,last-1)
  //       The calling matrix is assumed to have rows==cols

  void subfromdiag(const statmatrix & d, unsigned first, unsigned last);

  // FUNCTION: elemmult
  // TASK: computes the elementwise product of the calling matrix and A

  void elemmult(const statmatrix & A);

  // FUNCTION: elemquot
  // TASK: computes the elementwise ratio B/A and assigns it to the calling
  //       matrix B

  void elemquot(const statmatrix & A);

  // FUNCTION: weightedsscp
  // TASK: Computes the weighted sums of squares and crossproducts
  //       matrix X'WX

  void weightedsscp(const statmatrix & X, const statmatrix & w);

  // FUNCTION: weightedsscp2
  // TASK: Computes the weighted sums of squares and crossproducts
  //       matrix (X Z)'W(X Z)

  void weightedsscp2(const statmatrix & X, const statmatrix & Z,
                     const statmatrix & w);

  // FUNCTION: weightedsscp_resp
  // TASK: Computes X'W y

  void weightedsscp_resp(const statmatrix & X, const statmatrix & y,
                          const statmatrix & w);

  // FUNCTION: weightedsscp_resp2
  // TASK: Computes (X Z)'W y

  void weightedsscp_resp2(const statmatrix & X, const statmatrix & Z,
                     const statmatrix & y, const statmatrix & w);

  // ------------------ functions for sorting a column -------------------------

  // FUNCTION: sort
  // TASK: sorts the statmatrix between the 'start' th and the 'ende' th row
  //       according to its 'col' th column

  void sort (int start ,int ende,int col);

  // FUNCTION: sort
  // TASK: sorts the 'col' th column of the statmatrix between the 'start' th
  //       and the 'ende' th row

  void sortcol (int start ,int ende,int col);

  // FUNCTION: indexinit
  // TASK: initializes a column vector with elements (i,0) = i

  void indexinit (void);

  // FUNCTION: indexsort
  // TASK: index sort of the 'col' th column within the 'start' th and 'ende'
  //       th row. after sorting 'indexcol' of index contains the rank vector
  //       of the 'col' th column of the calling matrix

  void indexsort (statmatrix<int> & index,int start,int ende,
                  int col,int indexcol) const;

  // FUNCTION: indexsort2d
  // TASK: index sort of the 'col' th column within the 'start' th and 'ende'
  //       th row followed by sorting accoring to the 'col2' th row.
  //       after sorting 'indexcol' of index contains the rank vector
  //       of the 'col' th column of the calling matrix

  void indexsort2d (statmatrix<int> & index,int start,int ende,
                    int col,int col2, int indexcol) const;

  // FUNCTION: rank
  // TASK: computes the ranks of the elements within the 'start'th and 'ende'th
  //       row of the 'col'th column.
  // Remark: The elements within the 'start'th and 'ende'th row of the 'col'th
  //         column must be sorted with "indexsort" beforehand.

  void rank (statmatrix<double> & rang,statmatrix<int> & index,
             int start,int ende,int col) const;

//------------------------- end: sorting a column ------------------------------


  // FUNCTION: sum
  // TASK: returns the sum of the elements of the col. column of the calling
  //       matrix

  T  sum (const unsigned & col) const;

  // FUNCTION: sum2
  // TASK: returns the sum of the squared elements of the col th column of the
  //       calling matrix

  T  sum2 (const unsigned & col) const;

  T  sum2 (const unsigned & col,const statmatrix<T> & weight) const;

  // FUNCTION: sum2
  // TASK: returns a column vector whose elements are the squared column sums of
  //       the calling matrix

  statmatrix<T> sum2();

  // FUNCTION: sumcomplete
  // TASK: returns the sum of the elements of the calling
  //       matrix

  T  sumcomplete(void) const;

  // FUNCTION: sum
  // TASK: returns a column vector whose elements are the column sums of the
  //       calling matrix

  statmatrix<T> sum() const;

  // FUNCTION: mean
  // TASK: computes the mean of column 'col' of the calling matrix

  T mean(const unsigned & col) const
    {
    return sum(col)/T(this->rows());
    }

  T mean(const unsigned & col,const statmatrix<T> & weight) const;

  // FUNCTION: meancomplete
  // TASK:  returns the mean of the elements of the calling matrix

  T meancomplete(void) const
    {
    return sumcomplete()/T(this->rows());
    }

  // FUNCTION: norm
  // TASK: computes the euclidean norm of column col

  T norm(unsigned col) const;

  // FUNCTION: norm
  // TASK: returns a column vector of the norms of the columns of the calling
  //       matrix

  statmatrix<T> norm();


  // FUNCTION
  // TASK: returns the euclidean distance between the col-th column and the
  //       colA-th column of A

  T euclidean_dist(unsigned col, const statmatrix<T> & A, const unsigned & colA) const;



  // FUNCTION: mean
  // TASK: computes the mean of the columns of the calling matrix
  //       returns a column vector of means

  statmatrix<T> mean() const;

  // FUNCTION: var
  // TASK: computes the variance of column 'col' of the calling matrix

  T var(const unsigned & col) const;

  T var(const unsigned & col,const statmatrix<double> & weight) const;

  // FUNCTION: min
  // TASK: returns the minimum of column 'col'

  T min(const unsigned & col) const;

  // FUNCTION: max
  // TASK: returns the maximum of column 'col'

  T max(const unsigned & col) const;


  // FUNCTION: quantile
  // TASK: returns the 'percent' (0 < percent < 100) percent quantile
  //       of the 'col' th column

  T quantile  (const T & percent,const unsigned & col) const;

  // FUNCTION: quantile
  // TASK: returns the 'percent' (0 < percent < 100) percent quantile
  //       of the 'col' th column
  //       index contains the index sort of the col-th column, i.e. it is
  //       assumed that the indexsort is already done

  T quantile  (const T & percent,const unsigned & col,
  statmatrix<int> & index) const;



  // FUNCTION: quantile
  // TASK: computes the 'percent' percent quantile of the columns of the
  //       calling matrix, returns a column vektor with the quantiles

  statmatrix<T> quantile  (T percent);

  // FUNCTION: autocorr
  // TASK: computes the autocorrelation function for the 'col'. column and
  //			  lag 'lag'

  T             autocorr  (const unsigned & lag,const unsigned & col) const;

  // FUNCTION: autocorr
  // TASK: computes the autocorrelation function for lags 1 to 'lag' for
  //       all columns
  //       returns a  lag x column matrix of autocorrelations

  statmatrix<T> autocorr  (const unsigned & lag) const;

  // FUNCTION: autocorr
  // TASK: computes the autocorrelation function for lags 'beginlag' to 'endlag'
  //       for column 'col'
  //       returns a column vector of autocorrelations

  statmatrix<T> autocorr (const unsigned & beginlag,const unsigned & endlag,
                          const unsigned & col) const;


  // FUNCTION: cov
  // TASK: computes the covariance matrix of the calling matrix

  statmatrix<T> cov();

  statmatrix<T> corr();

  T compute_quadform(const statmatrix<T> & x,const unsigned & c=0);

  	//	Datenzeiger zugreifbar

  T *getV() const { return this->m_v; }


   statmatrix<T> strike (unsigned int k);

   statmatrix<T> get_cov_iX (int i, int j);

   statmatrix<T> partial_var(void);

   // FUNCTION: round
   // TASK: rounds the elements of the calling matrix between min and max to the
   //       precision specified in digits
   void round(const int digits, const unsigned mincol, const unsigned maxcol,
              const unsigned minrow, const unsigned maxrow);

   // FUNCTION: round
   // TASK: rounds the elements of the calling matrix to the precision specified
   //       in digits
   void round(const int digits);

   // FUNCTION: check_sorted
   // TASK: checks whether the col-th row of the calling matrix is sorted

   bool check_ascending(unsigned & col);

  };

typedef statmatrix<double> datamatrix;

// FUNCTION: multdiagback
// TASK: Computes X*D, where D is a diagonal matrix whose elements are stored
//       in d

//statmatrix<double> multdiagback(datamatrix X, const datamatrix & d);
template <class T>
statmatrix<T> multdiagback(statmatrix<T> X, const statmatrix<T> & d);

// FUNCTION: multdiagback
// TASK: Computes D*X, where D is a diagonal matrix whose elements are stored
//       in d

//statmatrix<double> multdiagfront(datamatrix X, const datamatrix & d);
template <class T>
statmatrix<T> multdiagfront(statmatrix<T> X, const statmatrix<T> & d);

#if !defined(STATMAT_CPP_INCLUDED)
#include "statmat.cpp"
#endif

#endif
