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



#if !defined (REALOBS_INCLUDED)

#define REALOBS_INCLUDED

#include"../export_type.h"
#if defined(MICROSOFT_VISUAL)
#include<limits>
#else
#include"../values.h"
#endif

#include<iostream>
#include<math.h>
#include<cmath>
#include"Random.h"

namespace realob
{

// missing value


#if defined(MICROSOFT_VISUAL)
  const double NA = DBL_MAX;

#else
  const double NA = MAXDOUBLE;
#endif
//------------------------------------------------------------------------------
//--------------------------- CLASS realobs ------------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE realobs
  {

  protected:


  //----------------------- PROTECTED VARIABLES --------------------------------

  // value of the real observation (can be NA = missing)

  double value;


  public :


  //------------------------ PUBLIC FUNCTIONS ----------------------------------

  // DEFAULT CONSTRUCTOR
  // TASK: value is set to NA (missing)

  realobs()
	 {
	 value = NA;
	 }

  // CONSTRUCTOR
  // TASK: value =v

  realobs(const double & v)
	 {
	 value = v;
	 }

  // COPY CONSTRUCTOR

  realobs(const realobs & o)
	 {
	 value = o.value;
	 }

  // DESTRUCTOR

  ~realobs() {}

  // OVERLOADED ASSIGNMENT OPERATOR
  // TASK: assigns realobs 'o' to the realobs

  const realobs & operator=(const realobs & o)
	 {
	 value = o.value;
	 return *this;
	 }

  // OVERLOADED ASSIGNMENT OPERATOR
  // TASK: assigns double value 'o' to the realobs

  const realobs & operator=(const double & v)
	 {
	 value = v;
	 return *this;
	 }

  // FUNCTION: getvalue

  const double & getvalue() const
	 {
	 return value;
	 }

  // BINARY ADDITION OPERATOR
  // ADDITIONAL INFORMATION:
  // NA + o = NA
  // o + NA = NA

  realobs operator+(const realobs & o) const
	 {
	 if ((value == NA) || (o.value == NA))
		return NA;
	 else
		return value+o.value;
	 }

  // BINARY ADDITION OPERATOR
  // ADDITIONAL INFORMATION:
  // v + NA = NA

  friend realobs __EXPORT_TYPE operator+(const double v,const realobs & o)
	 {
	 if (o.value==NA)
		return NA;
	 else
		return v+o.value;
	 }

  // BINARY ADDITION OPERATOR
  // ADDITIONAL INFORMATION:
  // NA + v = NA

  realobs operator+(const double v) const
	 {
	 return v + *this;
	 }

  // UNARY + OPERATOR

  realobs operator+()
	 {
	 return *this;
	 }

  // ++ OPERATOR

  realobs operator++()
	 {
	 if (value == NA)
		return NA;
	 else
		{
		value++;
		return *this;
		}
	 }

  // BINARY SUBSTRACTION OPERATOR
  // ADDITIONAL INFORMATION:
  // NA - o = NA

  realobs operator-(const realobs & o) const
	 {
	 if ((value == NA) || (o.value == NA))
		return NA;
	 else
		return value - o.value;
	 }

  // BINARY SUBSTRACTION OPERATOR
  // ADDITIONAL INFORMATION:
  // v - NA = NA

  friend realobs __EXPORT_TYPE operator-(const double v,const realobs & o)
	 {
	 if (o.value == NA)
		return NA;
	 else
		return v - o.value;
	 }

  // BINARY SUBSTRACTION OPERATOR
  // ADDITIONAL INFORMATION:
  // NA - v = NA

  realobs operator-(const double v) const
	 {
	 if (value == NA)
		return NA;
	 else
		return value - v;
	 }

  // UNARY - OPERATOR

  realobs operator-()
	 {
	 return -value;
	 }

  // -- OPERATOR

  realobs operator--()
	 {
	 if (value == NA)
		return NA;
	 else
		{
		value--;
		return *this;
		}
	 }

  // BINARY MULTIPLICATION OPERATOR

  realobs operator*(const realobs & o) const
	 {
	 if ((value == NA) || (o.value == NA))
		return NA;
	 else
		return value * o.value;
	 }

  // BINARY MULTIPLICATION OPERATOR

  realobs operator*(const double & o) const
	 {
	 if (value == NA)
		return NA;
	 else
		return value * o;
	 }

  // BINARY MULTIPLICATION OPERATOR

  friend realobs __EXPORT_TYPE operator*(const double v,const realobs & o)
	 {
	 if (o.value == NA)
		return NA;
	 else
		return v * o.value;
	 }

  // BINARY DIVISION OPERATOR

  realobs operator/(const realobs & o) const
	 {
	 if ((value == NA) || (o.value == NA) || (o.value == 0))
		return NA;
	 else
		return value/o.value;
	 }

  // BINARY DIVISION OPERATOR

  realobs operator/(const double v) const
	 {
	 if ((value == NA) || (v == 0))
		return NA;
	 else
		return value / v;
	 }

  // BINARY DIVISION OPERATOR

  friend realobs __EXPORT_TYPE operator/(const double v,const realobs & o)
	 {
	 if ((o.value == NA) || (o.value == 0))
		return NA;
	 else
		return v / o.value;
	 }

  // OVERLOADED COMPARISON OPERATORS

  // OVERLOADED == OPERATORS

  int operator==(realobs & o2) const
	 {
	 return value == o2.value;
	 }

  int operator==(const double & o2) const
	 {
	 return value == o2;
	 }

  friend int __EXPORT_TYPE operator==(double & o1, const realobs & o2)
	 {
	 return o1 == o2.value;
	 }

	 // OVERLOADED != OPERATORS

  int operator!=(const realobs & o2) const
	 {
	 return value != o2.value;
	 }

  int operator!=(const double  o2) const
	 {
	 return value != o2;
	 }

  friend int __EXPORT_TYPE operator!=(const double & o1,const realobs & o2)
	 {
	 return o1 != o2.value;
	 }

  // OVERLOADED > OPERATORS

  friend int __EXPORT_TYPE operator>(const realobs & o1,const realobs & o2)
	 {
	 return o1.value > o2.value;
	 }

  friend int __EXPORT_TYPE operator>(const realobs & o1,const double & o2)
	 {
	 return o1.value > o2;
	 }

  friend int __EXPORT_TYPE operator>(const double & o1,const realobs & o2)
	 {
	 return o1 > o2.value;
	 }

  // OVERLOADED < OPERATORS

  friend int __EXPORT_TYPE operator<(const realobs & o1,const realobs & o2)
	 {
	 return o1.value < o2.value;
	 }

  friend int __EXPORT_TYPE operator<(const realobs & o1,const double & o2)
	 {
	 return o1.value < o2;
	 }

  friend int __EXPORT_TYPE operator<(const double & o1,const realobs & o2)
	 {
	 return o1 < o2.value;
	 }

  // OVERLOADED >= OPERATORS

  friend int __EXPORT_TYPE operator>=(const realobs & o1,const realobs & o2)
	 {
	 return o1.value >= o2.value;
	 }

  friend int __EXPORT_TYPE operator>=(const realobs & o1,const double & o2)
	 {
	 return o1.value >= o2;
	 }

  friend int __EXPORT_TYPE operator>=(const double & o1,const realobs & o2)
	 {
	 return o1 >= o2.value;
	 }

  // OVERLOADED <= OPERATORS

  friend int __EXPORT_TYPE operator<=(const realobs & o1, const realobs & o2)
	 {
	 return o1.value <= o2.value;
	 }

  int operator<=(const double & o2) const
    {
    return value <= o2;
    }

  friend int __EXPORT_TYPE operator<=(const double & o1,const realobs & o2)
	 {
	 return o1 <= o2.value;
	 }

  // OVERLOADED << OPERATOR

  friend ostream & __EXPORT_TYPE operator<<(ostream & out,const realobs & o)
	 {
	 if (o.value==NA)
		return out << "NA";
	 else
		return out << o.value;
	 }

  // OVERLOADED >> OPERATOR

  friend istream & __EXPORT_TYPE operator>>(istream & in, realobs & o)
	 {
	 return in;
	 }

  // FUNCTION: applied
  // TASK: applies function func (which must be implemented for double values)
  //       to the realobs

  realobs applied(double (*func)(double))
	 {
	 if (value == NA)
		return NA;
	 else
		return func(value);
	 }


  // MATHEMATICAL FUNCTIONS

  friend realobs __EXPORT_TYPE sqrt(realobs & o);

  friend realobs __EXPORT_TYPE abs(realobs & o);

  friend realobs __EXPORT_TYPE exp(realobs & o);

  friend realobs __EXPORT_TYPE cos(realobs & o);

  friend realobs __EXPORT_TYPE sin(realobs & o);

  friend realobs __EXPORT_TYPE log(realobs & o);

  friend realobs __EXPORT_TYPE log10(realobs & o);

  friend realobs __EXPORT_TYPE pow(const realobs & o,const realobs & p);

  friend realobs __EXPORT_TYPE pow(realobs & o, double & p);

  friend realobs __EXPORT_TYPE pow(double o,realobs & p);

  friend realobs __EXPORT_TYPE floor(realobs & o);


//--------------------- statistical functions ----------------------------------


// FUNCTION: _uniform
// TASK: returns a in(0,1) unformly distributed random number

friend __EXPORT_TYPE realobs _uniform(void)
  {
  return randnumbers::uniform();
  }


// FUNCTION: _normal
// TASK: returns a standard normal distributed random number

friend __EXPORT_TYPE realobs _normal(void)
  {
  return randnumbers::rand_normal();
  }


// FUNCTION: _exponetial
// TASK: returns a random number, which is eponentially distributed with
//       parameter 'lambda'

friend __EXPORT_TYPE realobs _exponential(realobs lambda)
  {
  return randnumbers::rand_expo(lambda.value);
  }


  };

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

}  // end: namespace realob

#endif

