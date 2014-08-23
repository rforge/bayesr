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



#if !defined(SPARSEMATRIX_INCLUDED)

#define SPARSEMATRIX_INCLUDED

#include"../export_type.h"
#include"statmat.h"

#include"map.h"

using namespace std;

double __EXPORT_TYPE norm(const datamatrix & v);

// FUNCTION: norm
// TASK: computes the euclidian norm of the col th column

double __EXPORT_TYPE norm (const datamatrix & v, const unsigned col);


//------------------------------------------------------------------------------
//------------------------- class: SparseMatrix -------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE SparseMatrix
  {

  protected:

  unsigned cols;
  unsigned rows;
  vector< vector<double> > values;
  vector< vector<unsigned> > nonseros;


  public:

  // DEFAULT CONSTRUCTOR

  SparseMatrix(void)
    {
    cols = 1;
    rows = 1;
    values = vector< vector<double> >(1,vector<double>());
    nonseros = vector< vector<unsigned> >(1,vector<unsigned>());
    }

  // CONSTRUCTOR 0
  // TASK: creates a sparse matrix with size row x col
  //       maximum number of nonsero elements before memory reallocation
  //       = 'maxnonseros'

  SparseMatrix(const unsigned & row,const unsigned & col,
               const unsigned & maxnonseros = 0 );

  // CONSTRUCTOR 1
  // TASK: creates a sparse matrix containing the elements of the datamarix 'm'
  //       if opimize = true, storage size will be minimized

  SparseMatrix(const datamatrix & m,const bool optimize = false);

  // COPY CONSTRUCTOR

  SparseMatrix(const SparseMatrix & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const SparseMatrix & operator=(const SparseMatrix & m);

  // OVERLOADED () OPERATOR
  // TASK: returns the i,j th element

  double operator()(const unsigned & row,const unsigned & col) const;

  // FUNCTION: put
  // TASK: sets K(row,col) = v

  void put(const unsigned & row,const unsigned & col,const double & v);


  double compute_quadform(const datamatrix & x,const unsigned & col);

  double compute_condmean(const unsigned & i,const datamatrix & beta);

  // FUNCTION: get_rows
  // TASK: returns the number of rows of the matrix

  const unsigned & get_rows(void) const
    {
    return rows;
    }

  // FUNCTION: get_cols
  // TASK: returns the number of columns of the matrix

  const unsigned & get_cols(void) const
    {
    return cols;
    }


  unsigned getbandsize(void) const;

  // FUNCTION: kronecker
  // TASK: returns the kronecker product between the calling matrix and 'm'

  SparseMatrix kronecker(const SparseMatrix & m) const;

  // FUNCTION: getBlock
  // TASK: returns the block from
  //       row    'rowfirst' - 'rowlast'-1 and
  //       column 'colfirst' - 'collast'-1

  datamatrix getBlock(const unsigned & rowfirst,const unsigned & colfirst,
                      const unsigned & rowlast,const unsigned & collast);

  SparseMatrix getBlockasSparse(const unsigned & rowfirst,
                                const unsigned & colfirst,
                                const unsigned & rowlast,
                                const unsigned & collast);


  // FUNCTION: reorder

  SparseMatrix reorder(const statmatrix<int> & index);

  // FUNCTION: mult
  // TASK: multiplies the calling sparse matrix with the col-th column
  //       of 'vec' (starting row = 'a') and stores the result in 'res'
  // IMPORTANT:
  // - 'vec' must have proper size, i.e. 'vec' must have at least
  //   a + this.rows() and col+1-columns
  // - 'vec' must have proper size,
  //   i.e. rows (of the calling sparse matrix) = res.rows()

  void mult(const datamatrix & vec,const unsigned & a,const unsigned & col,
            datamatrix & res);

  // FUNCTION: add_mult
  // TASK: multiplies the calling sparse matrix with the 'col' th column
  //       of 'vec' (starting row = a)
  //       and ADDS the result to 'res'
  // IMPORTANT:
  // - 'vec' must have proper size, i.e. 'vec' must have at least
  //   a + this.rows() and col+1-columns
  // - 'vec' must have proper size,
  //   i.e. rows (of the calling sparse matrix) = res.rows()

  void add_mult(const datamatrix & vec,const unsigned & a,const unsigned & col,
                datamatrix & res);

  // FUNCTION: substr_mult
  // TASK: multiplies the calling sparse matrix with the 'col'-th column of
  //      'vec' (starting row = a)
  //       and substracts the result from 'res'
  // IMPORTANT:
  // - 'vec' must have proper size, i.e. 'vec' must have at least
  //   a + this.rows() and col+1-columns
  // - 'vec' must have proper size,
  //   i.e. rows (of the calling sparse matrix) = res.rows()

  void substr_mult(const datamatrix & vec,const unsigned & a,
                   const unsigned & col, datamatrix & res,
                   const unsigned & resrow=0);


  void print(ostream & o);

  void print2(ostream & o);


  // DESTRUCTOR

  ~SparseMatrix(void) {};

  };


  //----------------------------------------------------------------------------
  //---------------- Functions for computing penalty matrices ------------------
  //----------------------------------------------------------------------------

  // FUNCTION: Kmrf
  // TASK: returns the penalty matrix for MRF with characteristics stored in map

  SparseMatrix __EXPORT_TYPE Kmrf(const MAP::map & m);

  // FUNCTION: Kmrflinear
  // TASK: returns the penalty matrix for a 2 dimensional first order random
  //       walk


  SparseMatrix __EXPORT_TYPE Kmrflinear(const unsigned & nr1,const unsigned & nr2);

  // FUNCTION: Krw1
  // TASK: returns the penalty matrix for first order random walk

  SparseMatrix __EXPORT_TYPE Krw1(const vector<double> & weight);

  // FUNCTION: Krw2
  // TASK: returns the penalty matrix for second order random walk

  SparseMatrix __EXPORT_TYPE Krw2(const vector<double> & weight);

  // FUNCTION: Kseason
  // TASK: returns the penalty matrix for a sesonal component with period 'per'

  SparseMatrix __EXPORT_TYPE Kseason(unsigned per,unsigned s);

#endif
