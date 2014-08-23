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


// Schablonen fuer die Klassen LinkedList und LinkedListIterator

#if !defined(TLINKLST_H_INCLUDED)
#include "tlinklst.h"
#endif

#include <limits.h>

#if !defined(NDEBUG)
#include <iostream>
#endif

using std::cerr;
using std::endl;

//	Vorwaertssuche in der Speicherliste kann mit FIND_FORWARD
//	eingeschaltet werden.

#define TLINKLST_FIND_FORWARD  1
#define TLINKLST_LOG_FILE      0

//	int growSize
//
//	Betrag, um den eine Liste waechst, wenn eine Speicheranforderung
//	nicht erfuellt werden kann. Der tatsaechliche Betrag, um den die
//	Speicheranforderung der Liste vergroessert wird, ergibt sich aus
//	dieser Klassenvariablen und der Groesse eines Knotens der Liste.

#if defined(NO_STATIC_DATA)
const int growSize = 8192;
#else
template <class T>
int LinkedList<T>::growSize = 8192;
#endif

//	void insert(T v)
//
//	Ein Datenlement am Kopf der Liste einfuegen
//
// Parameter :
//
//	v - Ein Datenelement
//
//	Bemerkung :
//
//	Durch die Einfuegeoperation kann die Liste ungueltig
//	werden.

template <class T>
void LinkedList<T>::
addHead(const T &v)
{
	assert(m_valid);

	//	Neuen Knoten anlegen

	unsigned at = newNode();
	if (!m_valid)
	{
#if !defined(NDEBUG)
		cerr << "LinkedList::addHead(T) : Speichermangel!" << endl;
#endif
		return;
	}

	ListNode<T> *aNode = getNode(at);
	assert(aNode);

	//	Datenkomponente des Knotens besetzen

	aNode->v = v;

	//	Knoten am Kopf der Liste einhaengen

	if (m_head)
	{
		ListNode<T> *oldHead = getNode(m_head);
		oldHead->prev = at;
	}
	aNode->next = m_head;
	aNode->prev = 0;
	m_head = at;
	if (!m_tail)
		m_tail = at;

	//	Laenge der Liste hat sich erhoeht

	m_len++;
}

template <class T>
void
LinkedList<T>::
addTail(const T &v)
{
	assert(m_valid);

	//	Neuen Knoten anlegen

	unsigned at = newNode();
	if (!m_valid)
	{
#if !defined(NDEBUG)
		cerr << "LinkedList::addHead(T) : Speichermangel!" << endl;
#endif
		return;
	}

	ListNode<T> *aNode = getNode(at);
	assert(aNode);

	//	Datenkomponente des Knotens besetzen

	aNode->v = v;

	//	Knoten am Ende der Liste einhaengen

	if (m_tail)
	{
		ListNode<T> *oldTail = getNode(m_tail);
		oldTail->next = at;
	}
	aNode->prev = m_tail;
	aNode->next = 0;
	m_tail = at;
	if (!m_head)
		m_head = at;

	//	Laenge der Liste hat sich erhoeht

	m_len++;

}

//	int addAfter(const T &node, const T &v)
//
//	Ein Datenlement hinter dem Element mit der Addresse
//	node einfuegen
//
//  Parameter :
//
//	node - Die Addresse eines existierenden Elements der
//	Liste
//  v - Ein neues Datenelement
//
//	Ergebnis :
//
//	0 - Der Knoten "node" wurde nicht gefunden oder die Liste
//	    wurde wegen Speicherplatzmangels ungueltig
//	1 - Die Einfuegeoperation war erfolgreich
//
//	Bemerkung :
//
//	Durch die Einfuegeoperation kann die Liste ungueltig
//	werden.
//
//	Voraussetzung :
//
//	Die Liste darf nicht ungueltig sein

template <class T>
int
LinkedList<T>::
addAfter(const T &node, const T &v)
{
   assert(m_valid);

   //	Durch "node" bezeichneten Vorgaenger finden

    unsigned prev = findNode(&node);
	if (!prev)
	{
#if !defined(NDEBUG)
		cerr << "LinkedList::addAfter(const T&, T)"
			  << " : Knoten nicht gefunden : "
			  << node << endl;
#endif
		return 0;
	}


	ListNode<T> *prevNode = getNode(prev);
	ListNode<T> *nextNode = prevNode->next ? getNode(prevNode->next) : 0;

	//	Neuen Knoten anlegen

	unsigned at = newNode();
	if (!m_valid)
	{
#if !defined(NDEBUG)
		cerr << "LinkedList::addAfter(const T&, T) "
			  << " : Speichermangel!" << endl;
#endif
		return 0;
	}

	ListNode<T> *aNode = getNode(at);

	//	Datenkomponente des Knotens besetzen

	aNode->v = v;

	//	Den neuen Knoten in die Liste einhaengen


	aNode->next = prevNode->next;
	aNode->prev = prev;

	prevNode->next = at;

	//	Der Knoten ist entweder der Vorgaenger des Nachfolgers seines
	//	Vorgaengers, oder das neue Ende der Liste

	if (nextNode != 0)
		nextNode->prev = at;
	else
		m_tail = at;

	//	Laenge der Liste hat sich erhoeht

	m_len++;

	//	Einfuegeoperation war erfolgreich

	return 1;
}

//	int addBefore(const T &node, const T &v)
//
//	Ein Datenlement vor dem Element mit der Addresse
//	node einfuegen
//
//  Parameter :
//
//	node - Die Addresse eines existierenden Elements der
//	       Liste
//  v - Ein neues Datenelement
//
//	Ergebnis :
//
//	0 - Der Knoten "node" wurde nicht gefunden oder die Liste
//	    wurde wegen Speicherplatzmangels ungueltig
//	1 - Die Einfuegeoperation war erfolgreich
//
//	Bemerkung :
//
//	Durch die Einfuegeoperation kann die Liste ungueltig
//	werden.
//
//	Voraussetzung :
//
//	Die Liste darf nicht ungueltig sein

template <class T>
int
LinkedList<T>::
addBefore(const T &node, const T &v)
{
   assert(m_valid);

   //	Durch "node" bezeichneten Vorgaenger finden

   unsigned next = findNode(&node);
	if (!next)
	{
#if !defined(NDEBUG)
		cerr << "LinkedList::addBefore(const T&, T)"
			  << " : Knoten nicht gefunden : "
			  << node << endl;
#endif
		return 0;
	}


	ListNode<T> *nextNode = getNode(next);
	ListNode<T> *prevNode = nextNode->prev ? getNode(nextNode->prev) : 0;

	//	Neuen Knoten anlegen

	unsigned at = newNode();
	if (!m_valid)
	{
#if !defined(NDEBUG)
		cerr << "LinkedList::insert(const T&, T) "
			  << " : Speichermangel!" << endl;
#endif
		return 0;
	}

	ListNode<T> *aNode = getNode(at);

	//	Datenkomponente des Knotens besetzen

	aNode->v = v;

	//	Den neuen Knoten in die Liste einhaengen


	aNode->next = next;
	aNode->prev = nextNode->prev;

	//	Der Neue Knoten ist entweder der Nachfolger des Vorgaengers
	//	seines Nachfolgers, oder der neue Kopf der Liste

	if (prevNode != 0)
		prevNode->next = at;
	else
		m_head = at;

	nextNode->prev = at;

	//	Laenge der Liste hat sich erhoeht

	m_len++;

	//	Einfuegeoperation war erfolgreich

	return 1;
}

// void removeHead()
//
//	Das erste Element der Liste entfernen
//
//	Voraussetzung :
//
//	Die Liste darf nicht leer sein

template <class T>
void LinkedList<T>::removeHead()
{
	assert(m_head);

	//	Nachfolger festhalten

	unsigned nextHead = getNode(m_head)->next;

	//	Knoten loeschen

	deleteNode(m_head);

	//	Listenanfang auf den Nachfolger setzen

    m_head = nextHead;
	if (m_head)
	{
		getNode(m_head)->prev = 0;
	}

   //	Laenge der Liste hat sich vermindert

	m_len--;
}


// void removeTail()
//
//	Das letzte Element der Liste entfernen
//
//	Voraussetzung :
//
//	Die Liste darf nicht leer sein

template <class T>
void
LinkedList<T>::
removeTail()
{
	assert(m_tail);

	//	Vorgaenger festhalten

	unsigned prevTail = getNode(m_tail)->prev;

	//	Knoten loeschen

	deleteNode(m_tail);

	//	Listenende auf den Vorgaenger setzen

    m_tail = prevTail;
	if (m_tail)
	{
		getNode(m_tail)->next = 0;
	}

   //	Laenge der Liste hat sich vermindert

	m_len--;
}

// int remove(const T &node)
//
//	Einen bestimmten Knoten aus der Liste entfernen
//
//	Parameter :
//
//	node - Adresse des Datenanteils des Knotens, der geloescht
//	werden soll
//
//	Ergebnis :
//
//	0 - Es gibt keinen Knoten mit Inhalt an der angegebenen
//	Adresse, oder dieser Knoten ist nicht in die Liste
//	eingehaengt.
//	1 - Der Knoten wurde erfolgreich entfernt

template <class T>
int
LinkedList<T>::
remove(const T &node)
{
   assert(m_head);

   // Den durch "node" bezeichneten Knoten finden

   unsigned curr = findNode(&node);
   if (!curr)
   		return 0;

   //	Nachfolger und Vorgaenger des Knotens festhalten

   ListNode<T> *currNode = getNode(curr);
   unsigned int next = currNode->next;
   unsigned int prev = currNode->prev;
   ListNode<T> *nextNode = next ? getNode(next) : 0;
   ListNode<T> *prevNode = prev ? getNode(prev) : 0;

	//	Knoten loeschen

	deleteNode(curr);

	//	Geloeschten Knoten aus der Liste entfernen.

	if (prevNode != 0)
	{
		prevNode->next = next;
	}
	else
	{
		m_head = next;
		if (m_head)
			getNode(m_head)->prev = 0;
	}

	if (nextNode != 0)
	{
		nextNode->prev = prev;
	}
	else
	{
		m_tail = prev;
		if (m_tail)
			getNode(m_tail)->next = 0;
	}

	return 1;
}

//	unsigned newNode()
//
//	Einen freien Knoten im Listenspeicher bereitstellen
//
//	Ergebnis :
//
//	Index (um 1 verschoben) eines freien Knotens im (evtl.
//	vergroesserten Listenspeicher, wenn der Knoten
//	bereitgestellt werden konnte,
//	0 - kein Knoten wurde bereitgestellt.


template <class T>
unsigned LinkedList<T>::newNode()
{
	assert(m_valid);

	//	Wenn die Liste voll ist, muss der Listenspeicher vergroessert
	//	werden, wenn ein neuer Knoten bereitgestellt werden muss.

	if (m_len == m_memSize)
	{
   		resize();
		if (m_valid)
	   		return newNode();
		else
      		return 0;
   }

   //	Nach einem unbelegten Knoten in der Liste suchen, diese Suche
   //	muss vom Ende der Liste her erfolgen, da dort deutlich eher
   //	mit einem unbelegten Knoten zu rechnen ist.


	register unsigned i;
	register ListNode<T> *work;

	//	Jetzt muß wirklich Speicher bereitgestellt sein, damit in
	//	diesem Speicher nach einem freien Knoten gesucht werden
	//	kann.

	assert(m_memSize);

#if (TLINKLST_FIND_FORWARD)
	for (i = 0, work = m_nodeMemory; i < m_memSize; i++, work++)
#else
	i = m_memSize - 1;
	work = &m_nodeMemory[ i ];
	for (; 1; i--, work--)
#endif
	{
#if (TLINKLST_LOG_FILE)
		cerr << "i = " << i
			  << ", work->next = " << work->next
			  << endl;
#endif
		if (work->next == UINT_MAX)
			return i + 1;

		//	Kein freier Knoten in der Liste der Knoten gefunden.
		//	Dieser Zustand sollte eigentlich nie erreicht werden.

#if (TLINKLST_FIND_FORWARD == 0)
		if (!i)
			return 0;
#endif
	}

	//	Bei Vorwaertssuche : Kein freier Knoten in der Liste gefunden

#if (TLINKLST_FIND_FORWARD)
	return 0;
#endif
}

//	void deleteNode(unsigned at)
//
//	Einen Knoten mit einem bestimmten Index loeschen (Speicherfreigabe,
//	nicht logisches entfernen aus der Liste)
//
//	Parameter :
//
//	at - Position (um 1 verschoben) an der ein Knoten geloescht wird

template <class T>
void LinkedList<T>::deleteNode(unsigned at)
{
	ListNode<T> *aNode = getNode(at);

	//	Listenelement mit einem durch den Default-Konsturktor der
	//	Klasse T erzeugten Objekt ueberschreiben

	T v = T();
	aNode->v = v;

	//	Knoten fuer die Speicherverwaltung der Liste als ungueltig
	//	markieren

	aNode->next = UINT_MAX;
	aNode->prev = UINT_MAX;
}

//	unsigned findNode(const T * node) const
//
//	Einen Knoten finden, dessen Nutzdaten an der Position
//	node stehen
//
//	Parameter :
//
//	node - Addresse der Nutzdaten eines Knotens
//
//	Ergebnis :
//
//	Knotenindex (Zeigeräquivalent) des Knoten oder 0, wenn
//	kein passender Knoten gefunden wurde.

template <class T>
unsigned LinkedList<T>::findNode(const T *node) const
{
	assert(m_len);

	//	Suche durch Bisektion zwischen dem ersten und dem letzten Knoten
	//	in der Speicherliste. Dieser Ansatz macht es erforderlich,
	//	dass Zeigervergleiche durchgefuehrt werden duerfen.

	unsigned first = 1;
	unsigned last  = m_memSize;

	while(last - first >  1 )
	{
		unsigned center = (first + last) / 2;
#if (TLINKLST_LOG_FILE)
		cerr << "first = " << first << '\t'
			  << "last = " << last << '\t'
			  << "center = " << center << endl;
#endif
		ListNode<T> *aNode = getNode(center);
		if (&(aNode->v) < node)
			first = center;
		else if (&(aNode->v) > node)
			last = center;
		else if (&(aNode->v) == node)
			return center;
	}

	//	Jetzt unterscheiden sich first und last hoechstens noch
	//	um eins. Der gesuchte Knoten ist entweder first, oder last
	//	oder keines von beiden.

	if (&(getNode(first)->v) == node)
		return first;
	else if (&(getNode(last)->v) == node)
		return last;
	else
		return 0;
}

//	void resize()
//
//	Nach einem Zugriffsversuch jenseits der Grenzen des bisherigen
//	Speichers den verfügbaren Speicher anpassen.
//
//	Bemerkung:
//
//	Durch diese Operation kann die Liste ungültig werden.

template <class T>
void LinkedList<T>::
resize()
{
  // Anzahl der Elemente festlegen, um die die Liste wachsen soll

  unsigned grow = growSize / sizeof(ListNode<T>);
  if (grow < 4)
    grow = 4;

  // Neuen Speicher fuer die vergroesserte Liste anlegen

  ListNode<T> *newMemory = new ListNode<T>[ m_memSize + grow ];
  if (!newMemory)
    {
    // Bei Misserfolg ist die Liste jetzt ungueltig

    fail();
    return;
    }

  // Bisherigen Inhalt der LIste in den neuen Speicherbereich kopieren
  // (durch Anwendung des assignment-Operators, da moeglicherweise eine
  // tiefe Kopie notwendig ist.

  ListNode<T> *workOld;
  ListNode<T> *workNew;
  register unsigned i;

  for (i = 0, workOld = m_nodeMemory, workNew = newMemory;
       i < m_memSize;
       i++, workOld++, workNew++)
    {
    workNew->v    = workOld->v;
    workNew->next = workOld->next;
    }

  // Restliche (neue) Knoten als frei markieren

  m_memSize += grow;
  for (; i < m_memSize; i++, workNew++)
    workNew->next = UINT_MAX;

  // evtl. alten Arbeitsspeicher freigeben

  if (m_nodeMemory)
    delete [] m_nodeMemory;

  // neu bereitgestellten Speicherbereich zum Speicherbereich der
  // Liste machen

  m_nodeMemory = newMemory;
}

template <class T>
void
Heap<T>::
insert(T v)
{
// GNU:
  if (this->empty())
//  if (empty())
    addTail(v);
  else
    {
  // GNU:
    ListNode<T> *curr = getNode(this->get_mhead());
//    ListNode<T> *curr = getNode(m_head);

    while(curr && (curr->v < v))
       curr = curr->next ? getNode(curr->next) : 0;
    if (curr)
      addBefore(curr->v, v);
    else
      addTail(v);
    }
}



