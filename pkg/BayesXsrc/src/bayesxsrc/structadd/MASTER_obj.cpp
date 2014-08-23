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


#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"

#endif


#include "MASTER_obj.h"



namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------------- CLASS: MASTER_OBJ ------------------------------
//------------------------------------------------------------------------------


MASTER_OBJ::MASTER_OBJ(void)
  {
  }


MASTER_OBJ::MASTER_OBJ(const MASTER_OBJ & o)
  {
  level1_likep = o.level1_likep;
  }


const MASTER_OBJ & MASTER_OBJ::operator=(const MASTER_OBJ & o)
  {
  if (this == &o)
    return *this;
  level1_likep = o.level1_likep;
  return *this;
  }



} // end: namespace MCMC



