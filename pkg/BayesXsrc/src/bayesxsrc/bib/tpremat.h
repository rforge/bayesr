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



// Deklaration einer template-Klasse PreMatrix
//
// Der Typ T muss ein "numerischer" Datentyp sein, der neben
// Vergleichen und einfachen arithmetischen Operationen auch
// bestimmte Funktionen kennt. 

#if !defined(TPREMAT_H_INCLUDED)
#define TPREMAT_H_INCLUDED

#include "tarray2d.h"

template <class T>
class PreMatrix : public Array2D<T>
{
public :

   // 1. Konstruktoren

   // Default-Konstruktor

   PreMatrix(void) : Array2D<T>() {}

   // Konstruktor fuer eine uninitalisierte PreMatrix

   PreMatrix(unsigned rows, unsigned cols = 1) : 
   Array2D<T>(rows, cols) { }

   // Konstruktor fuer eine mit einem festen Wert initialisierte
   // PreMatrix

   PreMatrix(unsigned rows, unsigned cols, const T init) : 
   Array2D<T>(rows, cols, init) { }

   // Kopierkonstruktor

   PreMatrix(const PreMatrix &init) : 
   Array2D<T>(Array2D<T>(init)) { }

   // Destruktor

   virtual ~PreMatrix() {}

   // Zweistellige Operatoren

   PreMatrix operator+(const PreMatrix &m) const;
   PreMatrix operator-(const PreMatrix &m) const;
   PreMatrix operator*(const PreMatrix &m) const;
   PreMatrix operator*(const T v) const;
   PreMatrix operator/(const T v) const;

   // Concatenation
		
   PreMatrix operator&(const PreMatrix &m) const;
   PreMatrix operator|(const PreMatrix &m) const;

   PreMatrix kronecker(const PreMatrix &m) const;

   PreMatrix vec() const;
   PreMatrix vech() const;

   // unaere Operatoren
		
   PreMatrix operator+();
   PreMatrix operator-();

   // Skalarmultiplikation von links

   friend PreMatrix<T> operator *(const T v, const PreMatrix<T> &m)
      { return m.operator*(v); }

   const PreMatrix &operator+=(const PreMatrix &m);
   const PreMatrix &operator-=(const PreMatrix &m);
   const PreMatrix &operator*=(const PreMatrix &m)
      { *this = operator*(m); return *this; }
   const PreMatrix &operator*=(const T v);
   const PreMatrix &operator/=(const T v);
   const PreMatrix &operator|=(const PreMatrix &m);
   const PreMatrix &operator&=(const PreMatrix &m);

   // Funktionen von Matrizen
   //
   // transposed - Transponieren
   // sscp - SSCP-PreMatrix

   PreMatrix transposed(void) const;
   PreMatrix sscp(void) const;
	
   // Lesbare Ausgabe

   void prettyPrint(ostream &out);
   int prettyScan(istream &in);

   // Ausgabe mit einem bestimmten Delimiter
	
   void print(ostream& out, char delimiter) const;
   void print(ostream& out, char* delimiter) const;


   // Abfrage von Eigenschaften

   bool symmetric(const T epsilon = T(0)) const;
   bool zero(const T epsilon = T(0)) const;

   // Elementtypwertige Funktionen
		
   T det() const;
   T trace() const;

   PreMatrix getRow (unsigned i) const
      { PreMatrix<T> res; Array2D<T>::getRow(i).purge(res); return res; }
   PreMatrix getCol (unsigned j) const
      { PreMatrix<T> res; Array2D<T>::getCol(j).purge(res); return res; }
   PreMatrix getBlock(unsigned int rl, unsigned int ru, unsigned int cl, 
			 unsigned int cu) const;

   PreMatrix getRowBlock(unsigned int rl, unsigned int ru) const;
   PreMatrix getColBlock(unsigned int cl, unsigned int cu) const;
									 
   PreMatrix strikedOut(unsigned int row, unsigned int col)	const
      { PreMatrix<T> res; Array2D<T>::strikedOut(row, col).purge(res); return res; }
   PreMatrix strikedOutRow(unsigned int row)	const
      { PreMatrix<T> res; Array2D<T>::strikedOutRow(row).purge(res); return res; }
   PreMatrix strikedOutCol(unsigned int col)	const
      { PreMatrix<T> res; Array2D<T>::strikedOutCol(col).purge(res); return res; }

   // Lesen der Hauptdiagonalen.

   PreMatrix diag() const;

   // Besetzten der Hauptdiagonalen einer _neuen_ Matrix

   static PreMatrix diag(unsigned dim, const T v);
   static PreMatrix diag(const PreMatrix<T> &m);

   // Erzeugen einer Tridiagonalmatrix

   static PreMatrix tridiag(const PreMatrix &m, const PreMatrix &lm, 
			    const PreMatrix &um);

   PreMatrix blockdiag(const PreMatrix &m);
   
   const PreMatrix &operator=(const PreMatrix &from);

   // Anwenden einer Funktion
		
   PreMatrix applied(const PreMatrix &m, T (* f)(T, T)) const;
   PreMatrix applied(T (* f)(T)) const;

   PreMatrix vcat(const PreMatrix &bottom) const
     { PreMatrix <T> res; Array2D<T>::vcat(bottom).purge(res); return res; }
   PreMatrix hcat(const PreMatrix &right) const
     { PreMatrix<T> res; Array2D<T>::hcat(right).purge(res); return res; }

   PreMatrix luinverse() const;

   PreMatrix inverse() const
     { return luinverse(); }

   // L"ost ein lineares Gleichungsystem: *this * x = b mit LU-Zerlegung

   PreMatrix solve (const PreMatrix &bIn) const;

protected :

   PreMatrix decompLU(int *Index, int *IsEven = 0, int unique = 0) const;
   static PreMatrix backsubstLU(const PreMatrix &LU, const PreMatrix &bIn, 
                                int *Index );
};

template <class T>
inline
PreMatrix<T>
PreMatrix<T>::
getBlock(unsigned int rl, unsigned int ru, unsigned int cl, unsigned int cu) const
{ 
   PreMatrix<T> res; 

   Array2D<T>::getBlock(rl, ru, cl, cu).purge(res); 
   return res;
}

template <class T>
inline
PreMatrix<T>
PreMatrix<T>::
getColBlock(unsigned int cl, unsigned int cu) const
{
	PreMatrix<T> res;

	Array2D<T>::getColBlock(cl, cu).purge(res);
	return res;
}

template <class T>
inline
PreMatrix<T>
PreMatrix<T>::
getRowBlock(unsigned int rl, unsigned int ru) const
{
	PreMatrix<T> res;

	Array2D<T>::getRowBlock(rl, ru).purge(res);
	return res;
}

#if defined(TEMPL_INCL_DEF)
#	if defined(CC_SOURCE)
#		include <tpremat.cc>
#	else
#		include "tpremat.cpp"
#	endif
#endif

#endif
