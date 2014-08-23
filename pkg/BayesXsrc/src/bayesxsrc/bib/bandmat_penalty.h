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



#if !defined(BANDMATRIX_PENALTY)

#define BANDMATRIX_PENALTY

#include"../export_type.h"
#include "bandmat.h"


  //----------------------------------------------------------------------------
  //----------------- Functions for computing penalty matrices -----------------
  //----------------------------------------------------------------------------

  // FUNCTION: Krw1
  // TASK: returns the penalty matrix for a first order random walk

  bandmatdouble Krw1band(const vector<double> & weight);

  // FUNCTION: Krw2
  // TASK: returns the penalty matrix for a second order random walk

  bandmatdouble Krw2band(const vector<double> & weight);

  // FUNCTION: Kseason
  // TASK: returns the penalty matrix for a sesonal component with period 'per'

  bandmatdouble Kseasonband(const unsigned & per,const unsigned & s);

  // FUNCTION Kmrfband
  // TASK: returns the penalty matrix for MRF with characteristics stored in map

  bandmatdouble Kmrfband(const MAP::map & m);

  // FUNCTION: Kmrflinear
  // TASK: returns the penalty matrix for a columnwise first order random
  //       walk

  bandmatdouble Kmrflinearband(const unsigned & nr1,const unsigned & nr2);

#endif
