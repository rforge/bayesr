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



#include "design_pspline_group.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//---------- CLASS: DESIGN_pspline_group implementation of member functions ----
//------------------------------------------------------------------------------


void DESIGN_pspline_group::compute_orthogonaldecomp(void)
  {

  }


void DESIGN_pspline_group::init_data(const datamatrix & dm,const datamatrix & iv)
  {


  }


// DEFAULT CONSTRUCTOR

DESIGN_pspline_group::DESIGN_pspline_group(void)
  {

  }

// CONSTRUCTOR

DESIGN_pspline_group::DESIGN_pspline_group(DISTR * lp,FC_linear * fcp)
//: DESIGN(dp,fcl)
  {

  }


// COPY CONSTRUCTOR

DESIGN_pspline_group::DESIGN_pspline_group(const DESIGN_pspline_group & m)
: DESIGN(DESIGN(m))
  {

  }


// OVERLOADED ASSIGNMENT OPERATOR

const DESIGN_pspline_group & DESIGN_pspline_group::operator=(const DESIGN_pspline_group & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));

  return *this;
  }


void DESIGN_pspline_group::compute_penalty(void)
  {

  }

double  DESIGN_pspline_group::penalty_compute_quadform(datamatrix & beta)
  {

  }

void DESIGN_pspline_group::compute_basisNull(void)
  {


  }


void DESIGN_pspline_group::compute_XtransposedWX(void)
  {


  }


void DESIGN_pspline_group::compute_XtransposedWres(datamatrix & partres, double l)
  {


  }


void DESIGN_pspline_group::compute_precision(double l)
  {

  }


void DESIGN_pspline_group::compute_meaneffect(DISTR * level1_likep,double & meaneffect,
                                datamatrix & beta,datamatrix & meaneffectbeta,
                                bool computemeaneffect)

  {


  }


void DESIGN_pspline_group::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {

  }


void DESIGN_pspline_group::outoptions(GENERAL_OPTIONS * op)
  {

  }


} // end: namespace MCMC




