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


// DATE: 14.11.97


#if !defined (ERRORM_INCLUDED)

#define ERRORM_INCLUDED


#include<iostream>
#include"clstring.h"
#include<vector>


namespace errorm
{

	class messages : public vector<ST::string>
  {


  public:


  //--------------------------- PUBLIC FUNCTIONS -------------------------------

  // DEFAULT CONSTRUCTOR

	  messages(void) : vector<ST::string>() {}

  // COPY CONSTRUCTOR

	  messages(const messages & e) : vector<ST::string>(vector<ST::string>(e)) {}

  // OVERLOADED ASSIGNMENT OPERATORS

  const messages & operator=(const messages & e);

  // FUNKTION: insert_back

  void insert_back(const messages & e)
	 {
	 if (! e.empty())
		insert(end(),e.begin(),e.end());
	 }

  // OVERLOADED << OPERATOR

  friend ostream & operator<<(ostream & c, const messages & em);

  // FUNCTIONS: clear
  // TASK: deletes all errormessages

  void clear(void)
	 {
	 if (! empty())
		erase(begin(),end());
	 }


  };

}

#endif
