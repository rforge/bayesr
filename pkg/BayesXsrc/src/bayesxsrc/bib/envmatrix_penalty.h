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



#if !defined(ENVMATRIX_PENALTY)

#define ENVMATRIX_PENALTY

#include"../export_type.h"
#include "envmatrix.h"

  //----------------------------------------------------------------------------
  //----------------- Functions for computing penalty matrices -----------------
  //----------------------------------------------------------------------------

  // FUNCTION Kmrfenv
  // TASK: returns the penalty matrix for MRF with characteristics stored in map

  envmatrix<double> Kmrfenv(const MAP::map & m);

  // FUNCTION: Krw1env
  // TASK: returns the penalty matrix for a first order random walk

  envmatrix<double> Krw1env(const vector<double> & weight);
  envmatrix<double> Krw1env(const datamatrix & weight);

  // FUNCTION: Krw2env
  // TASK: returns the penalty matrix for a second order random walk

  envmatrix<double> Krw2env(const vector<double> & weight);

  envmatrix<double> Krw2env(const unsigned & nrpar);

  // FUNCTION: Krw3env
  // TASK: returns the penalty matrix for a third order random walk

  envmatrix<double> Krw3env(const unsigned & nrpar);

  // FUNCTION: Kseasonenv
  // TASK: returns the penalty matrix for a sesonal component with period 'per'

  envmatrix<double> Kseasonenv(const unsigned & per,const unsigned & s);

  // FUNCTION: Krw0env
  // TASK: returns the identity matrix (penalty for a random effect)

  envmatrix<double> Krw0env(const unsigned & nrpar);


#endif
