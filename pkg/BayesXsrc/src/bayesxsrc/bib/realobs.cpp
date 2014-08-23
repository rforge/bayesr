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





#include "realobs.h"

namespace realob
{

//------------------------------------------------------------------------------
//------------ CLASS realobs: implementation of member functions ---------------
//------------------------------------------------------------------------------


realobs sqrt(realobs & o)
  {
  if ( (o.value < 0) || (o.value == NA) )
	 return NA;
  else
	 return std::sqrt(o.value);
  }


realobs abs(realobs & o)
  {
  if (o.value == NA)
	 return NA;
  else
	 return fabs(o.value);
  }


realobs exp(realobs & o)
  {
  if (o.value == NA)
	 return NA;
  else
	 return std::exp(o.value);
  }


realobs cos(realobs & o)
  {
  if (o.value == NA)
	 return NA;
  else
	 return std::cos(o.value);
  }


realobs sin(realobs & o)
  {
  if (o.value == NA)
	 return NA;
  else
	 return std::sin(o.value);
  }


realobs log(realobs & o)
  {
  if  ( (o.value <= 0) || (o.value == NA) )
	 return NA;
  else
	 return std::log(o.value);
  }


realobs log10(realobs & o)
  {
  if  ( (o.value <= 0) || (o.value == NA) )
	 return NA;
  else
	 return std::log10(o.value);
  }


realobs pow(const realobs & o,const realobs & p)
  {
  if ((o.value == NA) || (p.value == NA))
	 return NA;
  else
	 return std::pow(o.value,p.value);
  }


realobs pow(realobs & o, double & p)
  {
  if ((o.value == NA) || (p == NA))
	 return NA;
  else
	 return std::pow(o.value,p);
  }


realobs pow(double o,realobs & p)
  {
  if ((o == NA) || (p.value == NA))
	 return NA;
  else
	 return std::pow(o,p.value);
  }

realobs floor(realobs & o)
  {
  if (o.value == NA)
	 return NA;
  else
	 return std::floor(o.value);
  }

// --------------------- forward friends decls ------------------

#if defined (__BUILDING_GNU)
__EXPORT_TYPE realobs _uniform(void);

realobs __EXPORT_TYPE sqrt(realobs & o);
realobs __EXPORT_TYPE abs(realobs & o);
realobs __EXPORT_TYPE exp(realobs & o);
realobs __EXPORT_TYPE cos(realobs & o);
realobs __EXPORT_TYPE sin(realobs & o);
realobs __EXPORT_TYPE log(realobs & o);
realobs __EXPORT_TYPE log10(realobs & o);
realobs __EXPORT_TYPE pow(const realobs & o,const realobs & p);
realobs __EXPORT_TYPE pow(realobs & o, double & p);
realobs __EXPORT_TYPE pow(double o,realobs & p);
realobs __EXPORT_TYPE floor(realobs & o);
#endif

} // end: namespace realobs

