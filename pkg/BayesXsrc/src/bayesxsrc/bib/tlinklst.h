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


// Deklaration der template-Klassen LinkedList und
// LinkedListIterator

#if !defined(TLINKLST_H_INCLUDED)
#define TLINKLST_H_INCLUDED

#include <assert.h>

template <class T>
class LinkedListIterator;

//	ListNode
//
//	Nutzdatenknoten einer einfach verzeigerten Liste

template <class T>
class ListNode
{
public:
	//	Nutzdaten des Knotens

	T v;

	//	Verweis auf das naechste Listenelement (Offset vom
	//	Anfang des Speichers fuer Knoten plus 1). Werte mit
	//	besonderer Bedeutung sind:
	//
	//	0 - Kein gueltiger Verweis
	//	UINT_MAX - nicht belegter Knoten

	unsigned next;

	//	Verweis auf das vorangehende Listenelement 
	//
	//	Werte mit besonderer Bedeutung: wie oben

	unsigned prev;
};

//	LinkedList - einfach verzeigerte Liste

template <class T>
class LinkedList
{
public :
        typedef ListNode<T> Node;

	//	1. Konstruktoren und Destruktoren -----------------------

	//	Default-Konstruktor

	LinkedList()
		{ create(); }

	//	Destruktor

	virtual ~LinkedList()
		{ discard(); }

protected:

	//	2. Elementfunktionen -----------------------

	//	Ein Element an Anfang der Liste einfuegen

	void addHead(const T &v);

	//	Ein Element am Ende der Liste einfuegen

	void addTail(const T &v);

	//	Ein Element hinter dem Element einfuegen, dessen Datenanteil
	//	an der Position node steht

	int addAfter(const T &node, const T &v);

	//	Ein Element vor dem Element einfuegen, dessen Datenanteil
	//	an der Position node steht

	int addBefore(const T &node, const T &v);

	//	Ein Element am Anfang der Liste entfernen

	void removeHead();

	//	Ein Element am Ende der Liste entfernen

	void removeTail();

	//	Das Element entfernen, dessen Datenanteil an der Position
	//	Node steht

	int remove(const T &node);

	//	Das Element am Anfang der Liste zurueckliefern

	const T& head() const
	{
		assert(m_memSize);
		return m_nodeMemory[ m_head - 1 ].v;
	}

	//	Das Element am Ende der Liste zurueckliefern

	const T& tail() const
	{
		assert(m_memSize);
		return m_nodeMemory[m_tail - 1].v;
	}

	//	Oeffentliche Elementfunktionen

public:

	//	Feststellen, ob die Liste leer ist

	int empty() const
		{ assert(m_valid); return m_head ? 0 : 1; }

	//	Die Laenge der Liste ermitteln

	unsigned len() const
		{ assert(m_valid); return m_len; }

	//	3. Operatoren -------------------------------------------

	//	negative Validitaetspruefung

	int operator!() const
		{ return m_valid ? 0 : 1; }

	//	 positive Validitaetspruefung

	operator int() const
		{ return m_valid ? 1 : 0; }

private :

	//	5. Instanzvariablen -------------------------------------

	//	Inhalt der doppelt verzeigerten Liste

	ListNode<T> *m_nodeMemory;

	//	Verfuegbarer Speicherplatz der verzeigerten Liste

	unsigned m_memSize;

protected:

	//	Index des ersten Elements in der Liste

	unsigned m_head;

	//	Index des letzten Elements in der Liste

	unsigned m_tail;

	//	Anzahl der Elemente in der Liste

	unsigned m_len;

	//	Gueltigkeitsanzeiger der Liste

	unsigned m_valid;

	//	6. Hilfsfunktionen --------------------------------------

	//	Einen freien Knoten der Liste ermitteln und seinen Index
	//	zurueckgeben

	unsigned newNode();

	//	Einen Knoten als geloescht markieren und den dynamisch
	//	fuer die Daten bereitgestellten Speicher freigeben

	void deleteNode(unsigned at);

	//	Einen Knoten an einer bestimmten Position (Index) als
	// Zeiger ermitteln

#if defined(NO_TEMPLATE_FRIENDS)
public:
#endif
	ListNode<T> *getNode(unsigned at) const
	{
		assert(at);
		assert(at <= m_memSize);
		return m_nodeMemory + (at - 1);
	}

// GNU:
    unsigned get_mhead(void)
	{
		return(m_head);
	}

#if defined(NO_TEMPLATE_FRIENDS)
protected:
#endif

	//	Einen Knoten finden, dessen Nutzdaten an der Position
	//	node stehen

	unsigned findNode(const T *node) const;

	//	Eine leere Liste anlegen

	void create()
	{
		m_nodeMemory = 0;
		m_memSize = 0;
		m_head = 0;
		m_tail = 0;
		m_valid = 1;
		m_len = 0;
	}

	//	Den Speicherbereich einer Liste neuen Anforderungen anpassen

	void resize();

	//	Den Speicherbereich einer Liste loeschen

	void discard()
		{ if (m_nodeMemory) delete [] m_nodeMemory; }

	//	Eine Liste als ungueltig markieren

	void fail()
		{ discard(); m_valid = 0; }

	//	7.	Klassenvariablen ------------------------------------

#if !defined(NO_STATIC_DATA)
	static int growSize;
#endif
};

template <class T>
class Stack : public LinkedList<T>
{
public:

	//	1. Konstruktoren und Destruktoren -----------------------

	//	Default-Konstruktor

	Stack()	: LinkedList<T>() { }

	//	Destruktor

	virtual ~Stack() { }

	//	2. Schnittstelle eines Stacks

	void insert(const T &v)
	{
		this->addHead(v);
	}

	void remove()
	{
// GNU:
		this->removeHead();
//		removeHead();
	}

	const T &top() const
	{
// GNU:
		return this->head();
//		return head();
	}
};

template <class T>
class Heap : public LinkedList<T>
{
public:
	//	1. Konstruktoren und Destruktoren -----------------------

	//	Default-Konstruktor

	Heap()	: LinkedList<T>() { }

	//	Destruktor

	virtual ~Heap() { }

	//	2. Schnittstelle eines Stacks

	void insert(T v); 

	void remove()
	{
// GNU:
		this->removeHead();
//		removeHead();
	}

	const T &top() 
	{
// GNU:
	  return this->head();
//	  return head();
	}

};

template <class T>
class Queue : public LinkedList<T>
{
public:

	//	1. Konstruktoren und Destruktoren -----------------------

	//	Default-Konstruktor

	Queue()	: LinkedList<T>() { }

	//	Destruktor

	virtual ~Queue() { }

	//	2. Schnittstelle einer Queue

	void insert(T v)
	{
		addTail(v);
	}

	void remove()
	{
// GNU:
		this->removeHead();
//		removeHead();
	}

	const T& first() const
	{
// GNU:
	  return this->head();
//	  return head();
	}
};

//	class List

//	Simuliert eine einfach verzeigerte Liste mit Einfügen am Ende

template <class T>
class List : public LinkedList<T>
{
public:

	//	1. Konstruktoren und Destruktoren -----------------------

	//	Default-Konstruktor

	List()	: LinkedList<T>() { }

	//	Destruktor

	virtual ~List() { }

	//	2. Schnittstelle einer Liste

	void insert(T v) 
	{
		addTail(v);
	}

	void insert(const T &node, T v)
	{
		addAfter(node, v);
	}

	void remove()
	{
// GNU:
		this->removeHead();
//		removeHead();
	}


	const T& head() const
	{
		return head();
	}

	const T& tail() const
	{
		return tail();
	}

	//	4. befreundete Klassen und Funktionen -------------------

#if !defined(NO_TEMPLATE_FRIENDS)
// GNU:
	friend class ListIterator;
//	friend class ListIterator<T>;

private:
// GNU:
	unsigned int ihead() { return this->get_mhead(); }
//	unsigned int ihead() { return m_head; }
#else
// GNU:
	unsigned int ihead() { return this->get_mhead(); }
//	unsigned int ihead() { return m_head; }
#endif
};

//	ListIterator - Durchlaufen einer verzeigerten Liste

template <class T>
class ListIterator
{
public :

	//	1. Konstruktoren und Destruktoren -----------------------

	//	Konstruktor, der die Verbindung des Iterators zu einer
	//	Liste herstellt

	ListIterator(List<T> &list) : 
        m_current(list.ihead()), m_list(list) { }

	//	Destruktor
		
	virtual ~ListIterator() { }

   //	2. Oeffentliche Elementfunktionen -----------------------

   //	Zugriff auf das aktuelle Element bei einer konstanten
   //	Liste

	const T &current() const
	{ 
		assert(m_current); 
		return m_list.getNode(m_current)->v;
	}

   //	Auf den Listenanfang zurueckstellen

	void reset()
   	{ m_current = m_list.ihead(); }

	//	3. Operatoren -------------------------------------------

	//	Gewoehnlicher Elementzugriff

   T &operator()()
	{
		assert(m_current);
		return m_list.getNode(m_current)->v;
	}

	//	Zum naechsten Listenelement uebergehen

	void operator++()
	{
		assert(m_current);
		m_current = m_list.getNode(m_current)->next;
	}

	//	Feststellen, ob die Liste bereits vollstaendig durchlaufen
	//	ist. Durch diesen Operator sind for-, while- und do-Schleifen
	//	ueber die Liste moeglich.

	operator int()
		{ return m_current ? 1 : 0; }

private :

	//	5. Instanzvariablen

	//	Aktuelles Element der Liste

	unsigned m_current;

	//	Liste, die der Iterator durchlaeuft

	List<T> &m_list;
};


#if defined(TEMPL_INCL_DEF)
#	if defined (CC_SOURCE)
#		include <tlinklst.cc>
#	else
#		include "tlinklst.cpp"
#	endif
#endif

#endif
