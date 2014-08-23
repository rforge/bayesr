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


// Muster fuer die outline-Funktionen der template-Klasse Matrix
//
// SCHALTER :
//
// SAVE_ALGORITHMS
//
// Mit dem Schalter kann zwischen einer "theorienahen" und einer
// Speicherdarstellungsnahen Variante der Algorithmen gewaehlt 
// werden. Die an der Speicherdarstellung orientierten Algorithmen
// sind fuer Fliesskommadatentypen den "theorienahen" Algorithmen 
// nicht deutlich ueberlegen (Sie sparen einige Prozent Rechenzeit)
// Interessanter sind diese Algorithmen bei Festkommaarithmetik,
// wenn der Overhead der Algorithmen gegenueber der Rechenzeit
// ins Gewicht faellt.

#include "tmatrix.h"
#include "tlinklst.h"
#include "tarray.h"

//#include <strstream>
#include <string.h>


template <class T>
Matrix<T> 
Matrix<T>::
root() const
{
   if (this->operator!() || this->rows() != this->cols())
      return Matrix<T>(0);

   unsigned n = this->rows( );

   Matrix<T> result = *this;
   if (!result)
      return result;

   if (n == 1)
      {
      T x = this->get(0, 0);

      if (T(0) < x)
	 result(0, 0) =  T(sqrt(x));
      else if (x == T(0))
	 result(0, 0) = T(0);
      else
	 return Matrix<T>( 0 );
      }
   else
      {
      unsigned i, j;

      for (i = 0; i < n; ++i)
	 for (j = i + 1; j < n; ++j)
	    result(i, j) = T(0);

      for (i = 0; i < n; ++i)
	 {
	 T sum = result(i, i);
	 for (j = 0; j < i; ++j)
	    {
	    T r = result(i, j);
	    sum -= r * r;
	    }
	 if (sum <= T(0))
	    return Matrix<T>(0);
	 result(i, i) =  T(sqrt(sum));
	 for (j = i + 1; j < n; j++)
	    {
	    unsigned k;

	    sum = result(j, i);
	    for (k = 0; k < i; ++k)
	       sum -= result(i, k) * result(j, k);
	    result(j, i) = sum / result(i, i);
	    }
	 }
      }
   return result;
}



template <class T>
Matrix<T> Matrix<T>::decompCholesky( void ) const
{
	Matrix<T> mat( *this );
	unsigned i, j, k;
	unsigned n = this->rows( );

	for ( j = 0; j < n; j++ )
	{
		for ( i = 0; i < j; i++ )
		{
			T h = mat( i, j );
			if ( mat( i, i ) != T( 0 ) )
				mat( i, j ) = h / mat( i, i );
			else
				return Matrix<T>( 0 );
			for ( k = i + 1; k <= j; k++ )
				mat( k, j ) = mat( k, j ) - h * mat( i, k );
		}
	}
	return mat;
}

template <class T>
Matrix<T> Matrix<T>::solveCholesky( const Matrix<T> &CH, unsigned index )
{
	Matrix<T> xvec( CH.rows( ), 1 );
	unsigned i, j;
	unsigned n = CH.rows( );
	T sum;

	if ( CH.get( 0, 0 ) <= T( 0 ) )
		return Matrix<T>( 0 );

	for( j = 0; j < n; j++ )
	{
		if ( j == index )
			sum = T( 1 );
		else
			sum = T( 0 );
		for ( i = 0; i < j; i++ )
			sum -= CH.get( i, j ) * xvec( i, 0 );
		xvec( j, 0 ) = sum;
		if ( CH.get( j, j ) <= T( 0 ) )
			return Matrix<T>( 0 );
	}
	for ( j = 0; j < n; j++ )
		xvec( j, 0 ) /= CH.get( j, j );

	assert( n );

	for ( j = n - 1; 1; j-- )
	{
		sum = xvec( j, 0 );
		for ( i = j + 1; i < n; i++ )
			sum -= CH.get( j, i ) * xvec( i, 0 );
		xvec( j, 0 ) = sum;
		if ( !j )
			break;
	}
	return xvec;
}


template <class T>
Matrix<T> 
Matrix<T>::
cinverse() const
{
   assert(!(this->operator!()));
   assert(this->rows() == this->cols());
   assert(this->symmetric());
   assert(this->rows() > 0);

   if (this->rows() == 1)
      {
      T v = this->get(0, 0);
      if (v == T(0))
	 return Matrix<T>(0);

      return Matrix<T>( 1, 1, T(1) / v );
      }

   Matrix<T> CH = decompCholesky();
   if (!CH)
      return Matrix<T>(0);

   Matrix<T> Inverse(this->rows(), this->cols());
   if (!Inverse)
      return Matrix<T>(0);

   unsigned j;

   for(j = 0; j < this->cols(); ++j)
      {
      Matrix<T> xvec = solveCholesky(CH, j);
      if (!xvec)
	 return Matrix<T>(0);

      Inverse.putCol(j, xvec);
      }
   return Inverse;
}


template <class T>
Matrix<T> Matrix<T>::solve(const Matrix &bIn) const
  {
  Matrix<T> res;
  PreMatrix<T>::solve(bIn).purge(res);
  return res;
  }



template <class T>
Matrix<T>
Matrix<T>::
inverse(const T) const
{
   if (this->operator!())
      return Matrix<T>(0);
   if (this->rows() != this->cols())
      return Matrix<T>(0);

   Matrix<T> res;

   PreMatrix<T>::inverse().purge(res);
   return res;
}



