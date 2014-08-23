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



// Deklaration der template-Klasse Array2D

#if !defined(TARRAY2D_H_INCLUDED)
#define TARRAY2D_H_INCLUDED

#include <iostream>
#include <limits.h>
#include <assert.h>
// #include <bool.h>

using std::istream;
using std::ostream;

//	Array2D - 2-dimensionales Feld fester Groesse

template <class T> class Array2D
{
   public :

	//	1. Konstruktoren und Destruktoren

	//	Default-Konstruktor

	Array2D()
		{ m_rows = 1;	m_cols = 1; create(); }

	//	Konstruktor fuer ein uninitalisiertes 2-dimensionales
	//	Feld

	Array2D(unsigned int rows, unsigned int cols)
		{ m_rows = rows; m_cols = cols; create(); }

	//	Konstruktor fuer ein mit einem festen Wert initialisiertes
	//	2-dimensionales Feld

	Array2D(unsigned int rows, unsigned int cols, const T init);

	//	Kopierkonstruktor

	Array2D(const Array2D &init);

	//	Destruktor

	virtual ~Array2D()
		{ discard(); }

	//	2. oeffentliche Elementfunktionen

	//	Zeilendimension

	const unsigned int rows() const
		{ return m_rows; }

	//   Spaltendimension

	const unsigned int cols() const
		{ return m_cols; }

	//	Lesender Zugriff auf ein Element des 2-dimensionalen
	//	Feldes. Dieser Aufruf ist auch dann zugelassen, wenn der
	//	Aufrufer nur ein const Array2D hat

	const T &get(unsigned int i, unsigned int j) const
	{
		assert(i < m_rows);
		assert(j < m_cols);
		return m_row[ i ][ j ];
	}

	//   Schreibender Zugriff auf ein Element des 2-dimensionalen
	//	Feldes

	void put(unsigned int i, unsigned int j, const T &v)
	{
		assert(i < m_rows);
		assert(j < m_cols);
		m_row[ i ][ j ] = v;
	}

	Array2D vcat(const Array2D &bottom) const;
	Array2D hcat(const Array2D &right) const;
	Array2D join(const Array2D &right,
                     bool (* pred)(const Array2D &,
                                   const Array2D &,
                                   unsigned int,
                                    unsigned int)) const;

	Array2D proj(bool (* pred)(const Array2D<T> &, unsigned int)) const;
	Array2D sel(bool (* pred)(const Array2D<T> &, unsigned int)) const;

	//	Blockweiser Ausschnitt aus einem zweidimensionalen Feld

	Array2D
	getBlock(unsigned int rl, unsigned int cl, unsigned int ru = UINT_MAX,
				unsigned int cu = UINT_MAX) const;

	//	Blockweises Einf"ugen in ein zweidimensionales Feld

	void putBlock(const Array2D<T> & m, unsigned int rl, unsigned int cl,
		          unsigned int ru = UINT_MAX, unsigned int cu = UINT_MAX);

	// getRow - i-te Zeile als Zeilenvektor zurueckliefern
	// putRow - i-te Zeile mit Zeilenvektor besetzen

	Array2D<T> getRow (unsigned i) const
	{
		assert(!(operator!()));
		assert(i < rows());
		return getBlock(i, 0, i + 1, cols());
	}

	void putRow (unsigned i, const Array2D<T> &from)
	{
		assert(!(operator!()));
		assert(i < rows());
		assert(from.rows() == 1);
		assert(from.cols() == cols() );

		putBlock(from, i, 0, i + 1, cols());
	}

   // getCol - j-te Spalte als Spaltenvektor zurueckliefern
   // putCol - j-te Spalte mit Spaltenvektor besetzen

   Array2D<T> getCol (unsigned j) const
   {
		assert(!(operator!()));
		assert(j < cols());

		return getBlock(0, j, rows(), j + 1);
   }

   void putCol (unsigned j, const Array2D<T> &from)
   {
	 	assert(!(operator!()));
		assert(j < cols());
		assert(from.cols() == 1);
		assert(from.rows() == rows() );

		putBlock(from, 0, j, rows(), j + 1);
   }


	Array2D<T> getRowBlock(unsigned int rl, unsigned int ru) const
	{
		assert(!(operator!()));
		assert(rl < ru);
		assert(ru <= rows());

		return getBlock(rl, 0, ru, cols());
	}

	void putRowBlock(unsigned int rl, unsigned int ru, const Array2D<T> &from)
	{
		assert(!(operator!()));
		assert(!from.operator!());
		assert(rl < ru);
		assert(ru <= rows());
		assert(from.cols() == cols());

		putBlock(from, rl, 0, ru, cols());
	}

	Array2D<T> getColBlock(unsigned int cl, unsigned int cu) const
	{
  		assert(!(operator!()));
		assert(cl < cu);
		assert(cu <= cols());

		return getBlock(0, cl, rows(), cu);
	}

	void putColBlock(unsigned int cl, unsigned int cu, const Array2D<T> &from)
	{
	  	assert(!(operator!()));
		assert(cl < cu);
		assert(cu <= cols());
		assert(from.rows() == rows());

		putBlock(from, 0, cl, rows(), cu);
	}

	Array2D<T> strikedOut(unsigned int row, unsigned int col) const;

	Array2D<T> strikedOutRow(unsigned int i) const;

	Array2D<T> strikedOutCol(unsigned int j) const;

	//	3. Operatoren

	//	Zugriff auf ein Element als lvalue ueber X(i, j)

	T &operator()(unsigned i, unsigned j)
	{
		assert(i < m_rows);
		assert(j < m_cols);
		return m_row[ i ][ j ];
	}

	//	Zugriff auf ein Element als rvalue ueber X(i, j)

	const T &operator()(unsigned i, unsigned j) const
	{
		assert(i < m_rows);
		assert(j < m_cols);
		return m_row[ i ][ j ];
	}

	//	Vergleichsoperator

	bool operator==(const Array2D<T> &to);

	//	Zuweisungsoperator

	const Array2D<T> &operator=(const Array2D<T> &from);

	//	negative Validitaetspruefung

	bool operator!() const
		{ return m_v ? false : true; }

	//	positive Validitaetspruefung

	operator bool() const
		{ return m_v ? true : false; }

	//	Anwendung einer Funktion auf eine Matrix

	Array2D<T> applied(T (* f)(T)) const;

	//	Anwendung einer Funktion auf zwei Matrizen

	Array2D<T> applied(const Array2D<T> &m, T (* f)(T, T)) const;


	//	4. Befreundete Klassen und Funktionen

	friend ostream &operator << (ostream &out, const Array2D<T> &fld)
		{ fld.writeOn(out); return out; }

	friend istream &operator >> (istream &in, Array2D<T> &fld)
		{ fld.readFrom(in); return in; }

	void purge(Array2D<T> &into);

protected :

	//	Klartextausgabe

	void writeOn(ostream &out) const;

	//	Klartexteingabe

	void readFrom(istream &in);


	//	4. Instanzvariablen

	//	Daten - Ein Datenbereich fuer das ganze 2-dimensionale
	//	Feld

	T *m_v;

	//	Zeiger auf Zeilenanfang

	T **m_row;

	//	Zeilendimension

	unsigned m_rows;

	//	Spaltendimension

	unsigned m_cols;

	//	5. Hilfsfunktionen

	//	Speicher freigeben

	void discard()
		{ if (m_v) { delete [] m_v; delete [] m_row; } }

	//	Speicher bereitstellen

	void create();

	//	Inhalte kopieren

	void copyContents(const Array2D &from);

	//	Datenzeiger zugreifbar

 	T *getV() const
   		{ return m_v; }
};

#if defined(TEMPL_INCL_DEF)
#	if defined(CC_SOURCE)
#		include <tarray2d.cc>
#	else
#		include "tarray2d.cpp"
#	endif
#endif

#endif
