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


#if !defined(INCL_ARITHM_H)
#define INCL_ARITHM_H

class arithm
{
public:
   arithm() : m_v(0.0) {};
   arithm(int v) : m_v(v) {}
   arithm(const arithm &from) { m_v = from.m_v; }

   arithm operator+(const arithm &ob) const
      { return arithm(m_v + ob.m_v); }
   arithm operator-(const arithm &ob) const
      { return arithm(m_v - ob.m_v); }
   arithm operator*(const arithm &ob) const
      { return arithm(m_v * ob.m_v); }
   arithm operator/(const arithm &ob) const
      { return arithm(m_v / ob.m_v); }

   arithm operator-()
     { return arithm(-m_v); }

   int operator==(const arithm &other) const
      { return m_v == other.m_v; } 
   int operator<(const arithm &other)
     { return m_v < other.m_v; } 
   int operator<=(const arithm &other)
     { return m_v <= other.m_v; } 
   arithm operator=(const arithm &from)
      { m_v = from.m_v; return *this; }
   arithm operator+=(const arithm &from)
      { m_v += from.m_v; return *this; }
   arithm operator-=(const arithm &from)
      { m_v -= from.m_v; return *this; }
   arithm operator*=(const arithm &from)
      { m_v *= from.m_v; return *this; }
   friend ostream &operator << (ostream &os, const arithm &ob)
      { os.precision(3); os << ob.m_v; return os; } 
   friend istream &operator >> (istream &is, arithm &ob)
      { is >> ob.m_v; return is; } 

 private:
   arithm(double v) : m_v(v) {}
   double m_v;
};

#endif
