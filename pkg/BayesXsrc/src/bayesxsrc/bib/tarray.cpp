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



//	Muster fuer die outline-Funktionen der template-Klasse Array

#if !defined(TARRAY_H_INCLUDED)
#include "tarray.h"
#endif

using std::ws;
using std::istream;
using std::ostream;

//	Array( unsigned size, const T init )
//
//	Konstruktor fuer mit festem Wert initialisiertes Array
//
//	Parameter :
//
//	size - Anzahl der Elemente des Arrays
//	init - Fester Wert, auf den alle Elemente des Arrays gesetzt werden
//
//	Bemerkung :
//
//	init ist ein konstanter Wert, keine Referenz, weil dadurch zur
//	Initialisierung auch ein T-Literal (d.h. eine echte Konstante)
//	zugelassen ist. Bei einer Referenz muss der Wert eine Adresse haben,
//	also eine echte Variable sein.
//
//	(Diese Einschraenkung gilt f"ur "altere Implementatinen)

template <class T>
Array<T>::Array( unsigned size, const T init )
{
	m_size = size;
	create( );
	if ( m_v )
	{
		register T *work;
		register unsigned i;

		for( work = m_v, i = 0; i < m_size; i++, work++ )
			*work = init;
	}
}

// Array( unsigned size, const T *init )
//
//	Konstruktor fuer mit mehreren Werten initialisiertes Array.
//
//	Parameter :
//
//	size - Anzahl der Elemente des Arrays
//	init - Zeiger auf Initialisierungswerte
//
//	Bemerkung :
//
//	Der Vektor init muss ein Zeiger auf mindestens size Objekte der Klasse
//	T sein.

template <class T>
Array<T>::Array( unsigned size, const T *init )
{
	m_size = size;
	create( );
	if ( m_v )
	{
		register T *workTo;
		const T *workFrom;
		register unsigned i;

		for ( workTo = m_v, workFrom = init, i = 0; i < m_size;
			 i++, workTo++, workFrom++ )
			*workTo = *workFrom;
	}
}

//	Array( const Array & init )
//
//	Kopierkonstruktor
//
//	Parameter :
//
//	init - Initialisierer

template <class T>
Array<T>::Array( const Array<T> &init )
{
	m_size = init.m_size;
	create( );
	if ( m_v )
		copyContents( init );
}

//	const Array &operator=( const Array &from)
//
//	Zuweisungsoperator
//
//	Parameter :
//
//	from - Zugewiesener Wert
//
//	Ergebnis :
//
//	Wert des Arrays nach der Zuweisung

template <class T>
const Array<T> &Array<T>::operator=( const Array<T> &from )
{
	discard( );
	m_size = from.m_size;
	create( );
	copyContents( from );
	return *this;
}

//	void writeOn( ostream &out ) const
//
//	Klartextausgabe
//
//	Parameter :
//
// out - Stream, in den das Objekt geschrieben werden soll


template <class T>
void Array<T>::writeOn( ostream &out ) const
{
	out << '[' << m_size << ']' << endl << '{' << endl;

	register T *work;
	register unsigned i;

	for ( work = m_v, i = 0; i < m_size; i++, work++ )
		out << *work << endl;
	out << '}';
}

template <class T>
void Array<T>::readFrom( istream &in )
{
	discard( );
	in >> ws;
	if ( in.get( ) != '[' )
		return;

	in >> m_size >> ws;
	if ( in.get( ) != ']' )
		return;
	in >> ws;
	if ( in.get( ) != '{' )
		return;

	create( );

	register T *work;
	register unsigned i;

	for ( work = m_v, i = 0; i < m_size; i++, work++ )
		in >> ws >> *work;
	in >> ws;
	if ( in.get( ) != '}' )
		discard( );
}


//	void copyContents( const Array &from )
//
// Inhalte kopieren
//
// Parameter :
//
// from - Array, dessen Inhalte kopiert werden sollen
//
//	Bemerkung :
//
//	Die Kopie erfolgt ueber den Zuweisungsoperator des Elementtyps.
//	Diese Vorgehensweise ist bei "flachen" Datentypen" u.U.
//	weniger effizient.

template <class T> void Array<T>::copyContents( const Array<T> &from )
{
	assert( m_v != 0 );
	assert( m_size > 0 );
	assert( from != 0 );
	assert( from.m_size == m_size );

	register unsigned int i;
	register T *workTo;
	register const T *workFrom;

	for( workTo = m_v, workFrom = from.m_v, i = 0; i < m_size;
		i++, workTo++, workFrom++ )
		*workTo = *workFrom;
}

//	operator==
//
//	Vergleichsoperator
//
//	Parameter:
//	other - zweiter Operand des Vergleichs
//
//	Ergebnis:
//	true (1) f"ur gleiche Arrays, false (0) f"ur Arrays, die sich
//	in der Gr"o"se oder in mindestens einem Element unterscheiden.

template <class T> int Array<T>::operator==(const Array<T> &other) const
{
	if (m_size != other.m_size)
		return 0;

	if (m_size)
	{
		assert(m_v != 0);
		assert(other.m_v != 0);
	}

	register T *work;
	register const T *workOther;

	for (work = m_v, workOther = other.m_v;
	     work < m_v + m_size;
		 ++work, ++workOther)
		if (*work != *workOther)
			return 0;

	return 1;
}

//	create
//
//	Array im Speicher anlegen. Wenn die Arraygr"o"se nicht
//	gr"o"ser als 0 ist, wird kein Array angelegt.
//
//	Wenn de Speicher f"ur das Array nicht angelegt werden
//	kann, wird das Array als ung"ultig erkl"art (d.h.
//	opeator int lifert 0 (false) und operator! liefert
//	1 (true)

template <class T>
void Array<T>::
create()
{
	assert( m_size > 0 );

	m_v = new T[ m_size ];
	if (m_v == 0)
	{
		m_size = 0;
	}
}
