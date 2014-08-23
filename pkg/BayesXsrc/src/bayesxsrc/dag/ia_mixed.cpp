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





#include "ia_mixed.h"
#include <algorithm>
#include <iterator>



namespace MCMC
{

	// DEFAULT CONSTRUCTOR:
	IA_MIXED::IA_MIXED(void): IA()
	{
	}
	


	// CONSTRUCTOR_1
	IA_MIXED::IA_MIXED(const datamatrix & d):IA(d)
	{
		
	}




	// CONSTRUCTOR_2
	// for interactions of order>2 (some day in future....) 
	IA_MIXED::IA_MIXED(unsigned order, const datamatrix & d):IA(order, d)
	{
	}


	// COPY CONSTRUCTOR
	IA_MIXED::IA_MIXED(const IA_MIXED & a) : IA (IA(a))
	{	
	}


	 // OVERLOADED ASSIGNMENT OPERATOR
	const IA_MIXED & IA_MIXED::operator=(const IA_MIXED & a)
	{
		if (this==&a)
		  return *this;

		// x = a.x;

		return *this;
	}





	


} //namespace MCMC
