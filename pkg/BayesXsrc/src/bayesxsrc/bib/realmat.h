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


// realmat.h 1.3 97/08/07 17:01:30
//
// Deklaration der Spezialisierung von Matrix auf einen Typ real, der vom
// Anwender in einer Header-Datei "real.h" bereitgestellt wird.

#include <tmatrix.h>

#if !defined(REALMAT_H_INCLUDED)
#define REALMAT_H_INCLUDED

class realMatrix : public Matrix<real>
{
public :

   // 1. Konstruktoren

   // Default-Konstruktor

   realMatrix() : Matrix<real>() { }

   // Konstruktor fuer eine uninitalisierte Matrix

   realMatrix(unsigned rows, unsigned cols = 1) : Matrix<real>(rows, cols) { }
  
   // Konstruktor fuer eine mit einem festen Wert initialisierte Matrix

   realMatrix(unsigned rows, unsigned cols, const real init) : 
   Matrix<real>(rows, cols, init) { }

  // Kopierkonstruktor

  realMatrix(const realMatrix &init) : Matrix<real>(Matrix<real>(init)) { }


  // Destruktor

  virtual ~realMatrix() { }



   // Zweistellige Operatoren

   realMatrix operator+(const realMatrix &m) const
      { realMatrix res;  Matrix<real>::operator+(m).purge(res); return res; }
   realMatrix operator-(const realMatrix &m) const
      { realMatrix res; Matrix<real>::operator-(m).purge(res); return res; }
   realMatrix operator*(const realMatrix &m) const
      { realMatrix res; Matrix<real>::operator*(m).purge(res); return res; }
   realMatrix operator*(const real v) const
      { realMatrix res; Matrix<real>::operator*(v).purge(res); return res; }
   realMatrix operator/(const real v) const
      { realMatrix res; Matrix<real>::operator/(v).purge( res ); return res; }
   realMatrix operator|(const realMatrix &m) const
      { realMatrix res; Matrix<real>::operator|(m).purge(res); return res; }
   realMatrix operator&(const realMatrix &m) const
      { realMatrix res; Matrix<real>::operator&(m).purge(res); return res; }

   // Einstellige Operatoren
	
   realMatrix operator+()
      { realMatrix res; Matrix<real>::operator+().purge(res); return res; }
   realMatrix operator-()
      { realMatrix res; Matrix<real>::operator-().purge(res); return res; }

   // Skalarmultiplikation von links

   friend realMatrix operator *(const real v, const realMatrix &m)
      { return m.operator*(v); }

   // Zuweisungsoperatoren

   const realMatrix &operator=(const realMatrix &m)
      { this->Matrix<real>::operator=(m); return *this; }
   const realMatrix &operator+=(const realMatrix &m)
      { this->Matrix<real>::operator+=(m); return *this; }
   const realMatrix &operator-=(const realMatrix &m)
	   { this->Matrix<real>::operator-=(m); return *this; }
   const realMatrix &operator*=(const realMatrix &m)
      { *this = operator*(m); return *this;}
   const realMatrix &operator*=(const real v)
	   { this->Matrix<real>::operator*=(v); return *this; }
   const realMatrix &operator/=(const real v)
	   { this->Matrix<real>::operator/=(v); return *this; }
   const realMatrix &operator&=(const realMatrix &m)
	   { this->Matrix<real>::operator&=(m); return *this; }
   const realMatrix &operator|=(const realMatrix &m)
	   { this->Matrix<real>::operator&=(m); return *this; }

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

   realMatrix kronecker(const realMatrix &m) const
      { realMatrix res; Matrix<real>::kronecker(m).purge(res); return res; }

   realMatrix vec() const
      { realMatrix res; Matrix<real>::vec().purge(res); return res; }
   realMatrix vech() const
      { realMatrix res; Matrix<real>::vech().purge(res); return res; }
   realMatrix transposed() const
      { realMatrix res; Matrix<real>::transposed().purge(res); return res; }
   realMatrix sscp() const
      { realMatrix res; Matrix<real>::sscp().purge(res);  return res; }

   // Anwendung von Funktionen auf Matrizen
		
   realMatrix applied(real (*f)(real)) const
      { realMatrix res; Array2D<real>::applied( f ).purge(res);  return res; }
      	
   realMatrix applied(const realMatrix &m, real (* f)(real, real)) const
      { realMatrix res; Array2D<real>::applied( m, f ).purge(res);  return res; }

   // Weitere Funktionen

   realMatrix getRow (unsigned i) const
      { realMatrix res; Array2D<real>::getRow(i).purge(res); return res; }
   realMatrix getCol (unsigned j) const
      { realMatrix res; Array2D<real>::getCol(j).purge(res); return res; }
   realMatrix getBlock(unsigned int rl, unsigned int ru,
		   unsigned int cl, unsigned int cu) const;

   realMatrix strikedOut(unsigned int row, unsigned int col) const
      { realMatrix res; Array2D<real>::strikedOut(row, col).purge(res); return res;}
   realMatrix strikedOutRow(unsigned int row) const;
   realMatrix strikedOutCol(unsigned int col) const;

   realMatrix getRowBlock(unsigned int rl, unsigned int ru) const
      { realMatrix res; Array2D<real>::getRowBlock(rl, ru).purge(res); return res; }
   realMatrix getColBlock(unsigned int cl, unsigned int cu) const
      { realMatrix res; Array2D<real>::getColBlock(cl, cu).purge(res); return res; }
   realMatrix diag() const
      { realMatrix res; Matrix<real>::diag().purge(res); return res; }
   realMatrix blockdiag(const realMatrix &m)
      { realMatrix res; Matrix<real>::blockdiag(m).purge(res); return res; }

   static realMatrix diag(unsigned int dim, const real v)
      { realMatrix res; Matrix<real>::diag(dim, v).purge(res); return res; }
   static realMatrix diag(const realMatrix &m)
      { realMatrix res; Matrix<real>::diag(m).purge(res); return res; }
   static realMatrix tridiag(const realMatrix &m, const realMatrix &lm, 
			 const realMatrix &um)
      { realMatrix res; Matrix<real>::tridiag(m,lm,um).purge(res); return res; }

   realMatrix luinverse() const
      { realMatrix res; Matrix<real>::luinverse().purge(res); return res; }
   realMatrix inverse(const real epsilon = 0.0) const
      { realMatrix res; Matrix<real>::inverse(epsilon).purge(res); return res; }

   realMatrix cinverse() const
      { realMatrix res; Matrix<real>::cinverse().purge(res); return res; }
   realMatrix root () const
      { realMatrix res; Matrix<real>::root().purge(res); return res; } 

   realMatrix vcat(const realMatrix &bottom) const
      { realMatrix res; Matrix<real>::vcat(bottom).purge(res); return res; }
   realMatrix hcat(const realMatrix &right) const
      { realMatrix res; Matrix<real>::hcat(right).purge(res); return res; }

   // L"ost ein lineares Gleichungsystem: *this * x = b mit LU-Zerlegung

   realMatrix solve (const realMatrix &bIn) const
      { realMatrix res; Matrix<real>::solve(bIn).purge(res); return res; }

   // L"ost ein lineares Gleichungssystem: *this * x = b mit Cholesky-
   // Zerlegung
   
   realMatrix csolve(const realMatrix &bIn) const
      { realMatrix res; Matrix<real>::csolve(bIn).purge(res); return res; }
protected:

   static realMatrix solveCholesky(const realMatrix &CH, unsigned index)
      { realMatrix res; Matrix<real>::solveCholesky(CH,index).purge(res); return res; }
   realMatrix decompCholesky() const
      { realMatrix res; Matrix<real>::decompCholesky().purge(res); return res; }


  // nur fuer realMatrix
public:
  
  bool containsNaNs() const;

  realMatrix eigenvalSym() const;
  realMatrix eigenvecSym() const;
		
};


inline 
realMatrix
realMatrix::
getBlock(unsigned int rl, unsigned int ru, unsigned int cl, unsigned int cu) 
const
{ 
   realMatrix res; 

   Array2D<real>::getBlock(rl, ru, cl, cu).purge(res);
   return res;
}


inline
realMatrix
realMatrix::
strikedOutRow(unsigned int row) const
{ 
   realMatrix res;
		
   Array2D<real>::strikedOutRow(row).purge(res);
   return res;
}


inline
realMatrix
realMatrix::
strikedOutCol(unsigned int col) const
{
   realMatrix res;

   Array2D<real>::strikedOutCol(col).purge(res);
   return res;
}



#endif



