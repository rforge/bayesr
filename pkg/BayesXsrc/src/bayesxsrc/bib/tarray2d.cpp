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



//	Muster fuer die outline-Funktionen der template-Klasse
//	Array2D

#if !defined(TARRAY2D_H_INCLUDED)
#include "tarray2d.h"
#endif

#include "tlinklst.h"
#include "tarray.h"
#include <assert.h>

//	Array2D(unsigned rows, unsigned cols, const T init)
//
//	Konstruktor fuer ein mit einem festen Wert initialisiertes
//	2-dimensionales Feld
//
//	Parameter :
//
//	rows - Zeilen des Feldes
//	cols - Spalten des Feldes
//	init - Fester Wert, auf den alle Elemente des Feldes gesetzt werden

template <class T>
Array2D<T>::
Array2D(unsigned rows, unsigned cols, const T init)
{
	m_rows = rows;
	m_cols = cols;
	create();
	if (m_v)
	{
		register T *work;
		unsigned size = m_rows * m_cols;
		register unsigned i;

		for(work = m_v, i = 0; i < size; i++, work++)
			*work = init;
	}
}

//	Array2D(const Array2D &init)
//
//   Kopierkonstruktor
//
//	Parameter :
//
//	init - Initialisierer

template <class T> 
Array2D<T>::
Array2D(const Array2D<T> &init)
{
	m_rows = init.m_rows;
	m_cols = init.m_cols;
	create();
	if (m_v)
		copyContents(init);
}


//	const Array2D &operator=(const Array2D &from)
//
//   Zuweisungsoperator
//
//	Parameter :
//
//	from - Zugewiesener Wert
//
//	Ergebnis :
//
//	Wert des Objekts nach der Zuweisung

template <class T>
const Array2D<T> &
Array2D<T>::
operator=(const Array2D<T> &from)
{
	discard();
	m_rows = from.m_rows;
	m_cols = from.m_cols;
	create();
	if (m_v)
		copyContents(from);
	return *this;
}

// Array2D vcat(const Array2D &bottom) const
//
// Untereinanderh"angen von 2Dimensionalen Arrays
//
// Parameter:
//
// bottom - Array, das unter den Empf"anger der Botschaft geh"angt wird
//
// Ergebnis:
// Neues Array2D geeigneter Dimension, in dem der Inhalt des Empf"aengers
// dieser Botschaft und das Argument der Botschaft untereinander angeordnet
// sind.

template <class T>
Array2D<T>
Array2D<T>::
vcat(const Array2D<T> &bottom) const
{
  assert(cols() == bottom.cols());

  Array2D<T> res(rows() + bottom.rows(), cols());

  res.putBlock(*this, 0, 0);
  res.putBlock(bottom, rows(), 0);
  return res;
}

// Array2D hcat(const Array2D &right) const
//
// Nebeneinanderhaengen vin 2-dimensionalen Arrays
//
// Parameter:
// right - Array, das rechts neben den Empf"anger geh"angt wird
//
// Ergebnis:
// Neues Array2D geeigneter Dimension, in dem der Inhalt des Empf"aengers
// dieser Botschaft und das Argument der Botschaft nebeneinander angeordnet
// sind.

template <class T>
Array2D<T>
Array2D<T>::
hcat(const Array2D<T> &right) const
{
  assert(rows() == right.rows());

  Array2D<T> res(rows(), cols() + right.cols());

  res.putBlock(*this, 0, 0);
  res.putBlock(right, 0, cols());
  return res;
}

// Array2D join(const Array2D<T> right, Predicate pred)
//
// Einen join "uber das Predikat pred durchf"uhren
//
// Parameter:
// right - rechter Operand des join
// pred - Kontrollierendes Predikat des joins
//
// Ergebnis: 

template <class T>	
Array2D<T> 
Array2D<T>::
join(const Array2D<T> &right, bool (* pred)(const Array2D<T> &, 
                                           const Array2D<T>&,
                                           unsigned int il,
                                           unsigned int ir)) const
{
  Stack< Array<T> > found;
  unsigned int width = cols() + right.cols(), i, j, k;

  for (i = 0; i < rows(); ++i)
    {
    for (j = 0; j < right.rows(); ++j)
      {
      if (pred(*this, right, i, j))
	{
	  Array<T> line(width);

	  for (k = 0; k < cols(); ++k)
	    line(k) = get(i, k);
	  for (; k < width; ++k)
	    line(k) = right.get(j, k - cols());
	  found.insert(line);
	}
      }
    }
  if (found.empty())
    return Array2D<T>(0, 0);

  Array2D<T> res(found.len(), width);
  for (i = 0; i < res.rows(); ++i)
    {
      Array<T> l = found.top();
      found.remove();
      for (j = 0; j < res.cols(); ++j)
	res(res.rows() - i - 1, j) = l(j);
    }
  assert(found.empty());
  return res;
}

//  Array2D proj(Predicate pred)
//
//  Projektion (Auswahl von Spalten nach einem Predikat
//
//  Parameter:
//  pred - liefert true f"ur Spalten, die ausgewaehlt werden sollen
//
//  Ergebnis:
//  Auswahlfeld

template <class T>
Array2D<T>
Array2D<T>::
proj(bool (* pred)(const Array2D<T> &, unsigned int j)) const
{
  Stack< Array<T> > found;
  unsigned int i, j;

  for (j = 0; j < cols(); ++j)
    {
    if (pred(*this, j))
      {
      Array<T> column(rows());

      for (i = 0; i < rows(); ++i)
	 column(i) = get(i, j);
      found.insert(column);
      }
    }
   
  if (found.empty())
    return Array2D<T>(0, 0);

  Array2D<T> res(rows(), found.len());
  for (j = 0; j < res.cols(); ++j)
    {
      Array<T> l = found.top();
      found.remove();
      for (i = 0; i < res.rows(); ++i)
	res(i, res.cols() - j - 1) = l(i);
    }
  assert(found.empty());
  return res;
}

// Array2D sel(Predicate pred)
//
// Selektion in einer Tafel durchf"uhren

template <class T>
Array2D<T>
Array2D<T>::
sel(bool (*pred)(const Array2D<T> &, unsigned int i)) const 
{
  Stack< Array<T> > found;
  unsigned int i, j;

  for (i = 0; i < rows(); ++i)
    {
    if (pred(*this, i))
      {
      Array<T> row(cols());

      for (j = 0; j < cols(); ++j)
	 row(j) = get(i, j);
      found.insert(row);
      }
    }
   
  if (found.empty())
    return Array2D<T>(0, 0);

  Array2D<T> res(found.len(), cols());
  for (i = 0; i < res.rows(); ++i)
    {
      Array<T> l = found.top();
      found.remove();
      for (j = 0; j < res.cols(); ++j)
	res(res.rows() - i - 1, j) = l(j);
    }
  assert(found.empty());
  return res;
}

//	void writeOn(ostream &out) const
//
//	Klartextausgabe
//
//	Parameter :
//
//	out - Stream, in den das Objekt geschrieben werden soll

template <class T>
void 
Array2D<T>::
writeOn(ostream &out) const
{
	out << '[' << m_rows << ',' << m_cols << ']' << endl << '{' << endl;

	unsigned i;
	register T * work = m_v;

	for (i = 0; i < m_rows; i++)
	{
		register unsigned j;

		out << '{' << endl;
		for (j = 0; j < m_cols; j++, work++)
			out << *work << endl;
		out << '}' << endl;
	}
	out << '}';
}

// void readFrom(istream &in)
//
// Klartexteingabe aus einer Datei

template <class T>
void Array2D<T>::readFrom(istream &in)
{
	//	Aktuellen Inhalt der Tafel loeschen

	discard();
	m_rows = 1;
	m_cols = 1;
	create();

	//	Dimension der Tafel aus der Datei lesen

	in >> ws;
	if (in.get() != '[')
		return;

	unsigned rows, cols;
	in >> rows >> ws;
	if (in.get() != ',')
		return;
	in >> cols >> ws;
	if (in.get() != ']')
		return;
	in >> ws;
	if (in.get() != '{')
		return;

	//	Tafel entsprechend der eingelesenen Dimension
	//	anlegen

	discard();
	m_rows = rows;
	m_cols = cols;
	create();

	//	Datenzeilen der Tafel einlesen

	register T * work = m_v;
	register unsigned i, j;
	for (i = 0; i < m_rows; i++)
	{
		in >> ws;
		if (in.get() != '{')
		{
			discard();
			m_rows = 1;
			m_cols = 1;
			create();
			return;
		 }

		//	Daten einlesen

		for (j = 0; j < m_cols; j++, work++)
			in >> *work;
		in >> ws;
		if (in.get() != '}')
		{
			discard();
			m_rows = 1;
			m_cols = 1;
			create();
			return;
		 }
	}
	in >> ws;
	if (in.get() != '}')
	{
		discard();
		m_rows = 1;
		m_cols = 1;
		create();
	}
}

//	void create()
//
//   Speicher bereitstellen
//
//	Bemerkung :
//
//	create hinterlaesst entweder ein vollstaendig initialisiertes Objekt
//	oder ein als ungueltig markiertes Objekt.

template <class T>
void Array2D<T>::create()
{
	if (m_rows == 0 || m_cols == 0)
	{
		m_v = 0;
		m_rows = 0;
		m_cols = 0;
	}
	else
	{
		unsigned size = m_rows * m_cols;
		m_v = new T[ size ];
		if (m_v)
		{
			m_row = new T*[ m_rows ];
			if (m_row)
			{
				register T* work;
				T **workRow;
				register unsigned i;

				//	Zeiger auf die Zeilen des 2D-Feldes besetzen:
				//	die Elemente von m_row zeigen auf die Zeilen
				//	des Feldes, die im Speicher "hintereinander"
				//	angeordnet sind.
				//
				//	m_row[ 0 ] = m_v;
				//	m_row[ 1 ] = m_v + (1 * m_cols);
				//	m_row[ 2 ] = m_v + (2 * m_cols);
				//		.
				//		.

				//	i = 0; work = m_v; workRow = m_row;
				//	while(i < m_rows)
				//	{
				//        *workRow = work;
				//		i++; work += m_cols; workRow++;
				//	}

				for(i = 0, work = m_v, workRow = m_row;
					i < m_rows;
					i++, work += m_cols, workRow++)
					*workRow = work;
			}
			else
			{
				//	Feld als ungueltig markieren

				delete [] m_v;
				m_v = 0;
				m_rows = 0;
				m_cols = 0;
			}
		}
	}
}

//	void copyContents(const Array2D &from)
//
//   Inhalte kopieren
//
//	Parameter :
//
//	from - Objekt, dessen Inhalt kopiert werden soll
//
//	Bemerkung :
//
//	Die Kopie erfolgt nicht fuer Speicherinhalte, sondern Elementweise

template <class T>
void Array2D<T>::copyContents(const Array2D<T> &from)
{
	assert(m_v != 0);
	assert(m_rows > 0);
	assert(m_cols > 0);
	assert(from != 0);
	assert(from.m_rows == m_rows);
	assert(from.m_cols == m_cols);

	register unsigned i;
	register T *workTo;
	register const T *workFrom;
	unsigned size = m_rows * m_cols;

	for(workTo = m_v, workFrom = from.m_v, i = 0;
		i < size;
		i++, workTo++, workFrom++)
		*workTo = *workFrom;
}


#if defined ( BODYHERE )

template <class T>
Array2D<T> Array2D<T>::getCol (unsigned j) const
{
   assert(!(operator!()));
   if (operator!())
   	return Array2D<T>(0, 0);

	assert(j < cols());

   Array2D<T> result(rows(), 1);
   assert (result);
	if (!result)
   	return Array2D<T>(0, 0);

   for (register unsigned i = 0; i < rows(); i++)
         result (i, 0) = get(i, j);

   return result;
}

 
template <class T>
void Array2D<T>::putCol(unsigned j, const Array2D<T> &from)
{
	assert(!(operator!()));
	assert(j < cols());
	assert(from);
	assert(from.cols() == 1);
	assert(from.rows() == rows() );

   for (register unsigned i = 0; i < rows(); i++)
         put(i,j,  from.get(i, 0));
}

#endif

template <class T>
Array2D<T> Array2D<T>::applied(T (* f)(T)) const
{
	assert(!(operator!()));
	
	Array2D<T> result(rows(), cols());

	for (register unsigned int i = 0; i < rows(); ++i)
		for (register unsigned int j = 0; j < cols(); ++j)
			result(i,j) = (*f)(get(i, j));
	return result;
}

//	Anwendung einer Funktion auf zwei Matrizen

template <class T>
Array2D<T> Array2D<T>::applied(const Array2D<T> &m, T (*f)(T, T)) const
{
	assert(!(operator!()));
	assert(!(!m));
	assert(rows() == m.rows());
	assert(cols() == m.cols());

	Array2D<T> result(rows(), cols());

	for (register unsigned int i = 0; i < rows(); ++i)
		for (register unsigned int j = 0; j < cols(); ++j)
			result(i,j) = (*f)(get(i, j), m.get(i, j));
	return result;
}

//	Blockweiser Ausschnitt aus einem zweidimensionalen Feld
//
//	Die Indizierung erfogt im "C"-Stil.

template <class T>	
Array2D<T> 
Array2D<T>::
getBlock(unsigned int rl, unsigned int cl, unsigned int ru, unsigned int cu) const
{
	assert(!(operator!()));
	assert(rl < rows());
	assert(cl < cols());
	assert(rl < ru);
	assert(cl < cu);
	assert(ru == UINT_MAX || ru <= rows());
	assert(cu == UINT_MAX || cu <= cols());

	if (operator!())
   		return Array2D<T>(0, 0);
	if (ru == UINT_MAX)
		ru = rows();
	if (cu == UINT_MAX)
		cu = cols();

	Array2D<T> result(ru - rl, cu - cl);
	assert(result);
	if (!result)
	{
		return Array2D<T>(0,0);
	}

	for (register unsigned int i = rl; i < ru; ++i)
		for (register unsigned int j = cl; j < cu; ++j)
			result.put(i - rl, j - cl, get(i, j));
	return result;
}
	
//	Blockweises Einf"ugen in ein zweidimensionales Feld

template <class T>
void 
Array2D<T>::
putBlock(const Array2D<T> & m, unsigned int rl, unsigned int cl, 
		 unsigned int ru, unsigned int cu)
{
	assert(!(operator!()));
	assert(rl < rows());
	assert(cl < cols());
	assert(rl < ru);
	assert(cl < cu);
	assert(ru == UINT_MAX || ru <= rows());
	assert(cu == UINT_MAX || cu <= cols());
	assert(rl + m.rows() <= rows() || ru <= rows());
	assert(cl + m.cols() <= cols() || cu <= cols());

	if (ru == UINT_MAX)
	{
		ru = rows();
		if (ru - rl > m.rows())
			ru = rl + m.rows();
	}
	if (cu == UINT_MAX)	
	{
		cu = cols();
		if (cu - cl > m.cols())
			cu = cl + m.cols();
	}

	assert(ru <= rows());
	assert(cu <= cols());
	assert(ru - rl <= m.rows());
	assert(cu - cl <= m.cols());


	for (register unsigned int i = rl; i < ru; ++i)
		for (register unsigned int j = cl; j < cu; ++j)
			put(i, j, m.get(i - rl, j - cl));
}

template <class T>
bool 
Array2D<T>::
operator==(const Array2D<T> &to)
{
   T *v, *w;

   assert(!operator!());
   assert(!to.operator!());

   if (rows() != to.rows() || cols() != to.cols())
      return false;
	
   for (v = m_v, w = to.m_v; v < m_v + (m_rows * m_cols); ++v, ++w)
      if (*v != *w)
         return false;
   return true;
}

template <class T>
Array2D<T> 
Array2D<T>::
strikedOut(unsigned int row, unsigned int col) const
{
	assert(!operator!());
	assert(row < rows());
	assert(col < cols());
	assert(rows() > 1);
	assert(cols() > 1);

	Array2D<T> result(rows() - 1, cols() - 1);

	register unsigned int i, j, k, l;

	for (i = 0, k = 0; k < rows(); ++i, ++k )
	{
		if (k == row)
		{
			--i;
			continue;
		}

		for (j = 0, l = 0; l < cols(); ++j, ++l)
		{
			if (l == col)
			{
				--j;
				continue;
			}
			result.put(i, j, get(k, l));
		}
	}
	return result;
}

template <class T>
Array2D<T>
Array2D<T>::
strikedOutRow(unsigned int row) const
{
	assert(!operator!());
	assert(row < rows());
	assert(rows() > 1);

	Array2D<T> result(rows() - 1, cols());

	register unsigned int i, j, k;

	for (i = 0, k = 0; k < rows(); ++i, ++k )
	{
		if (k == row)
		{
			--i;
			continue;
		}
		for (j = 0; j < cols(); ++j)
			result.put(i, j, get(k, j));
	}
	return result;
}

template <class T>
Array2D<T>
Array2D<T>::
strikedOutCol(unsigned int col) const
{
	assert(!operator!());
	assert(col < cols());
	assert(cols() > 1);

	Array2D<T> result(rows(), cols() - 1);

	register unsigned int i, j, l;

	for (i = 0; i < rows(); ++i)
	{
		for (j = 0, l = 0; l < cols(); ++j, ++l)
		{
			if (l == col)
			{
				--j;
				continue;
			}
			result.put(i, j, get(i, l));
		}
	}
	return result;
}

//	void purge(Array2D &)
//
//	Uebertragen des Speicheraufbaus ein ein andere Array2D.
//	Diese Operation entspricht konzeptionell einem "move"-
//	Konstruktor.

template <class T>
void
Array2D<T>::
purge(Array2D<T> &into)
{
	into.m_rows = m_rows;
	into.m_cols = m_cols;
	if (into.m_v) 
		delete [] into.m_v;
	into.m_v = m_v;
	if (into.m_row) 
		delete [] into.m_row;
	into.m_row = m_row;
	m_rows = m_cols = 0;
	m_row = 0;
	m_v = 0;
}



