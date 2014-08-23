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



#if !defined(STATMATRIXPENALTY_INCLUDED)

#define STATMATRIXPENALTY_INCLUDED

#include"../export_type.h"
#include "statmat.h"
#include "map.h"

namespace STATMAT_PENALTY
{

statmatrix<double> __EXPORT_TYPE Kmrf(const MAP::map & m);

statmatrix<double> __EXPORT_TYPE K2dim_pspline(const unsigned & nknots);

statmatrix<double> __EXPORT_TYPE K2dim_pspline_rw2(const unsigned & nknots, const unsigned & ox, const unsigned & oy);

statmatrix<double> __EXPORT_TYPE K2dim_pspline_biharmonic(const unsigned & nknots);
}

// FUNCTION: diffmat
// TASK: Computes a difference matrix with dimension (d-k x d), where k is the
//       difference order (k=1,2)

statmatrix<double> diffmat(const int k, const int d);

// FUNCTION: diffmat
// TASK: Computes a difference matrix with dimension (d-k x d), where k is the
//       difference order (k=0,1,2,...,d-1)

statmatrix<double> diffmat_k(const int k, const int d);

// FUNCTION: weighteddiffmat
// TASK: Computes a weighted difference matrix with dimension (d-k x d), where
//       k is the difference order and d is given by weight.size()

statmatrix<double> weighteddiffmat(const int k, const vector<double> & weight);

// FUNCTION: seasonalfactor
// TASK:     computes the factor of the penalty matrix for a seasonal effect
//           with period per and s time periods.

statmatrix<double> seasonalfactor(const unsigned & per, const unsigned & s);

// FUNCTION: seasonalX
// TASK:     computes the deterministic part of a seasonal effect
//           with period per and s time periods.

statmatrix<double> seasonalX(const unsigned & per, const unsigned & s);

// FUNCTION: rotate
// TASK:     help-function for the computation of eigenvalues and -vectors

void rotate(statmatrix<double> & a, const double & s, const double & tau,
            const int & i, const int & j, const int & k, const int & l);

// FUNCTION: tridiag
// TASK:     reduces the (symmetric) matrix a to tridiagonal form. On output,
//           the diagonal elements are stored in d and the subdiagonal elements
//           are stored in e. d and e are assumed to be n x 1 matrices.
//           The first element of e is zero on output.
// NOTE:     a is replaced by the orthogonal matrix effecting its transformation !!

void tridiag(statmatrix<double> & a, statmatrix<double> & d,
             statmatrix<double> & e);

// FUNCTION: eigentridiag
// TASK:     computes the eigenvalues and -vectors of a tridiagonal matrix. d
//           and e are assumed to be n x 1 matrices containing the diagonal
//           elements and the subdiagonal elements (with e(0,0) arbitrary)
//           respectively. If the eigenvectors of an originally tridiagonal
//           matrix are desired, z is input as the identity matrix, otherwise
//           z contains the output a from tridiag. On output d contains the
//           eigenvalues and z contains the eigenvectors. e ist destroyed on
//           output.

bool eigentridiag(statmatrix<double> & d, statmatrix<double> & e,
                  statmatrix<double> & z);

// FUNCTION: pythag
// TASK:     computes a^2+b^2 without destructive underflow or overflow

double pythag(const double & a, const double & b);

double sqr(const double & a);

double SIGN(const double & a, const double & b);


// FUNCTION: eigen
// TASK:     computes eigenvalues and eigenvectors of a and
//           stores them in values or vectors respectively.
//           The return value gives the number of iterations, that where needed
//           in the computation. If the return value equals 50, no convergence
//           could be achieved.
// NOTE:     a is modified in the computation of the eigenvalues and -vectors!!

int eigen(statmatrix<double> & a, statmatrix<double> & values,
           statmatrix<double> & vectors);

// FUNCTION: eigen2
// TASK:     computes eigenvalues and eigenvectors of a. On output, eigenvectors
//           are stored in a and eigenvalues are stored in d. If the return
//           value is false, no convergence could be achieved

bool eigen2(datamatrix & a, datamatrix & d);

// FUNCTION: eigensort
// TASK:     sorts the eigenvalues in values into descending order and
//           rearranges the columns of vectors correspondingly

void eigensort(datamatrix & values, datamatrix & vectors);

// FUNCTION: kronecker
// TASK:     computes the kronecker product of the matrices Aand B

statmatrix<double> kronecker(const statmatrix<double> & A, const statmatrix<double> & B);

void compare(const datamatrix & ref, const datamatrix & neu, double limit, unsigned col, const ST::string & colname, vector<ST::string> & out);
void compare_nonp(const ST::string & ref, const ST::string & neu, double limit, vector<ST::string> & out);



#endif
