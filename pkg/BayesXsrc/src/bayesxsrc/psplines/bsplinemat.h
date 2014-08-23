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



#ifndef bsplinematH
#define bsplinematH

#include"../export_type.h"
#include<deque>
#include "statmat.h"
#include "mcmc_nonpbasis.h"

using std::deque;

namespace MCMC
{


//---------------------------------------------------------------------------
//----------------------- class: bsplinemat -------------------------------
//---------------------------------------------------------------------------


class __EXPORT_TYPE bsplinemat
  {

  protected:

  datamatrix B;
  datamatrix BS;

  unsigned nrknots;
  unsigned degree;
  unsigned nrdiffobs;
  unsigned nrpar;

  knotpos knpos;

  vector<int> freq;
  vector<int> freqoutput;
  vector<int> index2;
  vector<int> begcol;

  deque<int> firstnonzero;
  deque<int> lastnonzero;
  deque<double> knot;

  statmatrix<int> index;

  datamatrix Bcolmean;


  void make_index(const datamatrix & md);

  void make_Bspline(const bool & deriv, const datamatrix & md, const bool & minnull = false);

  datamatrix bspline(const double & x);
  datamatrix bspline_derivative(const double & x);


  public:

  // DEFAULT CONSTRUCTOR

  bsplinemat(void)
    {
    }

  // CONSTRUCTOR

  bsplinemat(const datamatrix & data, const unsigned & nrk, const unsigned & degr,
                const MCMC::knotpos & kp, const bool & minnull = false, const deque<double> & k = deque<double>());

  // CONSTRUCTOR 2 (for derivatives)

  bsplinemat(const bool & deriv, const datamatrix & data, const unsigned & nrk, const unsigned & degr,
                const MCMC::knotpos & kp, const bool & minnull = false, const deque<double> & k = deque<double>());


  // COPY CONSTRUCTOR

  bsplinemat(const bsplinemat & bmat);

  // OVERLOADED ASSIGNMENT OPERATOR

  const bsplinemat & operator=(const bsplinemat & bmat);

  void mult(datamatrix & res, const datamatrix & beta);

  void mult_index(datamatrix & res, const datamatrix & beta);

  // DESTRUCTOR

  ~bsplinemat(){}

  };


}   // END: namespace MCMC

//---------------------------------------------------------------------------
#endif
