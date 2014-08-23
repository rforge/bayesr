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



#if !defined (MASTEROBJ)

#define MASTEROBJ

#include"../export_type.h"
#include"distr.h"

#if defined(JAVA_OUTPUT_WINDOW)
#include"adminparse_basic.h"
#endif


namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: MASTER_OBJ --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE MASTER_OBJ
  {


  protected:


  public:

  vector<DISTR *> level1_likep;

  // DEFAULT CONSTRUCTOR

  MASTER_OBJ(void);

  // COPY CONSTRUCTOR

  MASTER_OBJ(const MASTER_OBJ & o);

  // OVERLOADED ASSIGNMENT OPERATOR

  const MASTER_OBJ & operator=(const MASTER_OBJ & o);


  // DESTRUCTOR

  ~MASTER_OBJ() {}


  };

//------------------------------------------------------------------------------
//------------------------ End: CLASS MASTER_OBJ -------------------------------
//------------------------------------------------------------------------------

} // end: namespace MCMC

#endif
