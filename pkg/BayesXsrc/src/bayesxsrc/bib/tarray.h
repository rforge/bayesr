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


using namespace std;

//	Deklaration der template-Klasse Array

#if !defined(TARRAY_H_INCLUDED)
#define TARRAY_H_INCLUDED

#include <iostream>
#include <assert.h>

using std::istream;
using std::ostream;

//	Array - Feld fester Groesse

template <class T> class Array
{
   public :

	//	1. Konstruktoren und Destruktor

	//	Default-Konstruktor

	Array()
		{ m_size = 1; create(); }

	//	Konstruktor fuer uninitalisiertes Array

	Array(unsigned int size)
		{ m_size = size; create(); }

	//	Konstruktor fuer mit festem Wert initialisiertes Array

	Array(unsigned int size, const T init);

	//	Konstruktor fuer mit mehreren Werten initialisiertes Array.
	//	Der Vektor init muss ein Zeiger auf size Objekte der Klasse
	//	T sein.
	//
	//	Dieser Konstruktor erm"oglicht die Umwandlung eines built-in
	//	Arrays in ein Array der Bibliothek.

	Array(unsigned int size, const T *init);

	//	Kopierkonstruktor

	Array(const Array &init);

	//	Destruktor

	virtual ~Array()
		{ discard(); }

	//	2. "Offentliche Elementfunktionen

	//	Groesse - Anzahl der Elemente im Array

	const unsigned int size() const
		{ return m_size; }

	//	Lesender Zugriff auf ein Element des Arrays. Dieser ist
	//	auch dann zugelassen, wenn der Aufrufer nur ein const Array
	//	hat.

	const T &get(unsigned int at) const
		{ assert(at < m_size); return m_v[ at ]; }

	//   schreibender Zugriff auf ein Array-Element

	void put(unsigned int at, const T &v)
		{ assert(at < m_size); m_v[ at ] = v; }

	//	3. Operatoren

	//	Zugriff auf eine Array-Element als lvalue

	T &operator()(unsigned int at)
		{ assert(at < m_size); return m_v[at]; }

	const T &operator()(unsigned int at) const
		{ assert(at < m_size); return m_v[at]; }

	//	Zuweisungsoperator fuer Arrays

	const Array &operator=(const Array &from);

	//	negative Validitaetspruefung

	int operator!() const
		{ return m_v ? 0 : 1; }

	//	 positive Validitaetspruefung

	operator int() const
		{ return m_v ? 1 : 0; }

	//	4. Befreundete Klassen und Funktionen

	friend ostream &operator << (ostream &out, const Array<T> &fld)
		{ fld.writeOn(out); return out; }

	friend istream &operator >> (istream &in, Array<T> &fld)
		{ fld.readFrom(in); return in; }

	int operator==(const Array<T> &other) const;

  private :

	//	5. Instanzvariablen

	//	Daten - Zeiger auf die Elemente des Arrays

	T *m_v;

	//	Groesse - Elemente des Arrays

	unsigned int m_size;

	//	6. Hilfsfunktionen

	//	Klartextausgabe

	void writeOn(ostream &out) const;

	//	Klartexteingabe

	void readFrom(istream &in);

	//	Speicher freigeben

	void discard()
		{ if (m_v) delete [] m_v; m_v = 0; }

	//	Speicher bereitstellen

	void create();


	//	Inhalte kopieren

	void copyContents(const Array &from);
};

#if defined(TEMPL_INCL_DEF)
#	if defined(CC_SOURCE)
#		include <tarray.cc>
#	else
#		include "tarray.cpp"
#	endif
#endif

#endif
