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


// Deklaration einer template-Klasse Matrix
//
// Der Typ T muss ein "arithmetischer" Datentyp sein, der neben
// Vergleichen und einfachen arithmetischen Operationen auch
// bestimmte Funktionen kennt. Fuer T muss zumindest definiert
// sein:
//
// sqrt()   Quadratwurzel

#if !defined( TMATRIX_H_INCLUDED )

#define TMATRIX_H_INCLUDED

#include "tpremat.h"

template <class T>
class Matrix : public PreMatrix<T>
{
public :

   // 1. Konstruktoren

   // Default-Konstruktor

   Matrix() : PreMatrix<T>() { }

   // Konstruktor fuer eine uninitalisierte Matrix

   Matrix(unsigned rows, unsigned cols = 1) : PreMatrix<T>(rows, cols) { }
  
   // Konstruktor fuer eine mit einem festen Wert initialisierte Matrix

   Matrix(unsigned rows, unsigned cols, const T init) : 
   PreMatrix<T>(rows, cols, init) { }

   // Kopierkonstruktor

   Matrix(const Matrix &init) : PreMatrix<T>(PreMatrix<T>(init)) { }

   // Destruktor

   virtual ~Matrix() { }

   // Zweistellige Operatoren

   Matrix operator+(const Matrix &m) const
      { Matrix<T> res;  PreMatrix<T>::operator+(m).purge(res); return res; }
   Matrix operator-(const Matrix &m) const
      { Matrix<T> res; PreMatrix<T>::operator-(m).purge(res); return res; }
   Matrix operator*(const Matrix &m) const
      { Matrix<T> res; PreMatrix<T>::operator*(m).purge(res); return res; }
   Matrix operator*(const T v) const
      { Matrix<T> res; PreMatrix<T>::operator*(v).purge(res); return res; }
   Matrix operator/(const T v) const
      { Matrix<T> res; PreMatrix<T>::operator/(v).purge( res ); return res; }
   Matrix operator|(const Matrix<T> &m) const
      { Matrix<T> res; PreMatrix<T>::operator|(m).purge(res); return res; }
   Matrix operator&(const Matrix<T> &m) const
      { Matrix<T> res; PreMatrix<T>::operator&(m).purge(res); return res; }

   // Einstellige Operatoren
	
   Matrix operator+()
      { Matrix<T> res; PreMatrix<T>::operator+().purge(res); return res; }
   Matrix operator-()
      { Matrix<T> res; PreMatrix<T>::operator-().purge(res); return res; }

   // Skalarmultiplikation von links

   friend Matrix<T> operator *(const T v, const Matrix<T> &m)
      { return m.operator*(v); }

   // Zuweisungsoperatoren

   const Matrix &operator=(const Matrix &m)
      { this->PreMatrix<T>::operator=(m); return *this; }
   const Matrix &operator+=(const Matrix &m)
      { this->PreMatrix<T>::operator+=(m); return *this; }
   const Matrix &operator-=(const Matrix &m)
      { this->PreMatrix<T>::operator-=(m); return *this; }
   const Matrix &operator*=(const Matrix &m)
      { *this = operator*(m); return *this;}
   const Matrix &operator*=(const T v)
      { this->PreMatrix<T>::operator*=(v); return *this; }
   const Matrix &operator/=(const T v)
      { this->PreMatrix<T>::operator/=(v); return *this; }
   const Matrix &operator&=(const Matrix &m)
      { this->PreMatrix<T>::operator&=(m); return *this; }
   const Matrix &operator|=(const Matrix &m)
      { this->PreMatrix<T>::operator|=(m); return *this; }

   // Matrixwertige Funktionen von Matrizen
   //
   // Zweistellige:
   // kronecker - Kronecker-Produkt
   //
   // Einstellige:
   // vec - Vektorisierung
   // vech - horizontale Vektorisierung
   // inverse - Invertieren mit geeignetem Algorithmus
   // cinverse - Mit Choleskyzerlegung (positiv definit)
   // luinverse - Mit LU-Zerlegung (invertierbar)
   // transposed - Transponieren
   // sscp - SSCP-Matrix
   // root - Choleskyzerlegung

   Matrix kronecker(const Matrix &m) const
      { Matrix<T> res; PreMatrix<T>::kronecker(m).purge(res); return res; }

   Matrix vec() const
      { Matrix<T> res; PreMatrix<T>::vec().purge(res); return res; }
   Matrix vech() const
      { Matrix<T> res; PreMatrix<T>::vech().purge(res); return res; }
   Matrix transposed() const
      { Matrix<T> res; PreMatrix<T>::transposed().purge(res); return res; }
   Matrix sscp() const
      { Matrix<T> res; PreMatrix<T>::sscp().purge(res);  return res; }

   // Anwendung von Funktionen auf Matrizen
		
   Matrix applied(T (*f)(T)) const
		{ Matrix<T> res; Array2D<T>::applied(f).purge(res); return res; }
   Matrix applied(const Matrix &m, T (* f)(T, T)) const
	   { Matrix<T> res; Array2D<T>::applied(m, f).purge(res); return res; }

   // Weitere Funktionen

   Matrix getRow (unsigned i) const
      { Matrix<T> res; Array2D<T>::getRow(i).purge(res); return res; }
   Matrix getCol (unsigned j) const
      { Matrix<T> res; Array2D<T>::getCol(j).purge(res); return res; }
   Matrix getBlock(unsigned int rl, unsigned int ru,
		   unsigned int cl, unsigned int cu) const;

   Matrix strikedOut(unsigned int row, unsigned int col) const
      { Matrix<T> res; Array2D<T>::strikedOut(row, col).purge(res); return res;}
   Matrix strikedOutRow(unsigned int row) const;
   Matrix strikedOutCol(unsigned int col) const;

   Matrix getRowBlock(unsigned int rl, unsigned int ru) const
      { Matrix<T> res; Array2D<T>::getRowBlock(rl, ru).purge(res); return res; }
   Matrix getColBlock(unsigned int cl, unsigned int cu) const
      { Matrix<T> res; Array2D<T>::getColBlock(cl, cu).purge(res); return res; }
   Matrix diag() const
      { Matrix<T> res; PreMatrix<T>::diag().purge(res); return res; }
   Matrix<T> blockdiag(const Matrix &m)
      { Matrix<T> res; PreMatrix<T>::blockdiag(m).purge(res); return res; }

   static Matrix diag(unsigned int dim, const T v)
	   { Matrix<T> res; PreMatrix<T>::diag(dim, v).purge(res); return res; } 
   static Matrix diag(const Matrix<T> &m)
	   { Matrix<T> res; PreMatrix<T>::diag(m).purge(res); return res; }
   static Matrix tridiag(const Matrix<T> &m, const Matrix<T> &lm, 
			 const Matrix<T> &um)
	   { Matrix<T> res; PreMatrix<T>::tridiag(m,lm,um).purge(res); return res; }

   Matrix luinverse() const;
		
   Matrix inverse(const T epsilon = T(0)) const;

   Matrix cinverse() const;
   Matrix root () const; 

   Matrix vcat(const Matrix &bottom) const
	   { Matrix<T> res; PreMatrix<T>::vcat(bottom).purge(res); return res; }
   Matrix hcat(const Matrix &right) const
	   { Matrix<T> res; PreMatrix<T>::hcat(right).purge(res); return res; }

   // L"ost ein lineares Gleichungsystem: *this * x = b mit LU-Zerlegung

   Matrix solve (const Matrix &bIn) const;
  //    { Matrix<T> res; PreMatrix<T>::solve(bIn).purge(res); return res; }

   // L"ost ein lineares Gleichungssystem: *this * x = b mit Cholesky-
   // Zerlegung
   
   Matrix csolve(const Matrix &bIn) const;
protected:

   static Matrix solveCholesky(const Matrix &CH, unsigned index);
   Matrix decompCholesky() const;
};

template <class T>
inline 
Matrix<T>
Matrix<T>::
getBlock(unsigned int rl, unsigned int ru, unsigned int cl, unsigned int cu) 
const
{ 
   Matrix<T> res; 

   Array2D<T>::getBlock(rl, ru, cl, cu).purge(res);
   return res;
}

template <class T>
inline
Matrix<T>
Matrix<T>::
strikedOutRow(unsigned int row) const
{ 
   Matrix<T> res;
		
   Array2D<T>::strikedOutRow(row).purge(res);
   return res;
}

template <class T>
inline
Matrix<T>
Matrix<T>::
strikedOutCol(unsigned int col) const
{
   Matrix<T> res;

   Array2D<T>::strikedOutCol(col).purge(res);
   return res;
}

#if defined( TEMPL_INCL_DEF )
#	if defined( CC_SOURCE )
#		include <tmatrix.cc>
#	else
#		include "tmatrix.cpp"
#	endif
#endif

#endif
