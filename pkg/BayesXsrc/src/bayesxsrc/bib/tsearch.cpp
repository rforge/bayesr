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


// tsearcht.cc 3.1 97/08/05 01:57:40 
//
// AVL-Suchbaum

#include "tsearcht.h"

template <class T>
TreeNode<T>::
TreeNode( T init ) :
nodeCount(0),
nodeData(init) 
{
   leftTree  = 0;
   rightTree = 0;
}

template <class T>
void TreeNode<T>::
print( unsigned h, ostream &out )
{
   if ( leftTree )
      leftTree->print( h + 1, out );
   for ( unsigned i = 0; i < h; ++i )
      out << "  ";
   out << nodeData << endl;
   if ( rightTree )
      rightTree->print( h + 1, out );
}

template <class T>
void 
TreeNode<T>::
deleteTree( )
{
   if ( leftTree  ) 
   { 
      leftTree->deleteTree( );
      delete leftTree;
   }
   if ( rightTree ) 
   {
      rightTree->deleteTree( );
      delete rightTree;
   }
}

template <class T>
unsigned
TreeNode<T>::
height( ) const
{
   unsigned leftHeight  = leftTree ? leftTree->height( ) : 0;
   unsigned rightHeight = rightTree ? rightTree->height( ) : 0;

   return 1 + (leftHeight > rightHeight ? leftHeight : rightHeight);
}

template <class T>
int 
TreeNode<T>::
compare( const T &x )
{
   if ( x < nodeData )
      return -1;
   else if ( nodeData < x )
      return 1;
   else
      return 0;
}

template <class T>
void 
TreeNode<T>::
set( const TreeNode<T> &to )
{
   nodeData  = to.nodeData;
   nodeCount = to.nodeCount;
}

template <class T>
AVLNode<T>::
AVLNode( T Init )
:TreeNode<T>( Init )
{
   nodeBalance = 0;
}


// Der rechte Unterbaum des linken Unterbaus wird zum linken
// Unterbaum der Wurzel.
//   
// Die Wurzel wird zum rechten Unterbaum ihres bisherigen
// linken Unterbaumes.

template <class T>
void 
TreeNode<T>::
rotateLL( )
{
   TreeNode<T> *temp = leftTree;
   leftTree = temp->rightTree;
   temp->rightTree = this;
}

// Der linke Unterbaum des rechten Unterbaus wird zum rechten
// Unterbaum der Wurzel.
//
// Die Wurzel wird zum linken Unterbaum ihres bisherigen
// rechten Unterbaumes.

template <class T>
void
TreeNode<T>::
rotateRR( )
{
   TreeNode<T> *temp = rightTree;
   rightTree = temp->leftTree;
   temp->leftTree = this;
}

// Der rechte Unterbaum des linken Unterbaumes wird durch
// seinen linken Unterbaum ersetzt.
//
// Der linke Unterbaum wird zum linken Unterbaum seines
// bisherigen rechtenen Unterbaumes. (1. Rotation)
//
// Der rechte Unterbaum der linken Unterbaumes wird zum
// neuen linken Unterbaum der Wurzel.
//
// Die Wurzel wird zum rechten Unterbaum des bisherigen
// rechten Unterbaues des linken Unterbaumes, der die
// neue Wurzel bildet. (2. Rotation)

template <class T>
void
TreeNode<T>::
rotateLRLR( )
{
   TreeNode<T> *temp = leftTree;
   TreeNode<T> *sub = leftTree->rightTree;
   temp->rightTree = sub->leftTree;
   sub->leftTree = temp;
   leftTree = sub->rightTree;
   sub->rightTree = this;
}

// Der linke Unterbaum des rechten Unterbaumes wird durch
// seinen rechten Unterbaum ersetzt.
//
// Der rechte Unterbaum wird zum rechten Unterbaum seines
// bisherigen linken Unterbaumes. (1. Rotation)
//
// Der linke Unterbaum der rechten Unterbaumes wird zum
// neuen rechten Unterbaum der Wurzel.
//
// Die Wurzel wird zum linken Unterbaum des bisherigen
// linken Unterbaumes des rechten Unterbaumes, der die
// neue Wurzel bildet. (2. Rotation)

template <class T>
void
TreeNode<T>::
rotateRLRL( )
{
   TreeNode<T> *temp = rightTree;
   TreeNode<T> *sub = rightTree->leftTree;
   temp->leftTree = sub->rightTree;
   sub->rightTree = temp;
   rightTree = sub->leftTree;
   sub->leftTree  = this;
}

template <class T>
void 
TreeNode<T>::foreach(void (* f)(void *, void *), void *data)
{
	if (leftTree != 0)
		leftTree->foreach(f, data);
	f(&nodeData, data);
	if (rightTree != 0)
		rightTree->foreach(f, data);
}

template <class T>
SearchTree<T>::
SearchTree( )
{
   tree = 0;
}

template <class T>
SearchTree<T>::
~SearchTree( )
{
   if (tree)
      tree->deleteTree( );
}


template <class T>
void 
SearchTree<T>::
foreach(void (* f)(void *, void *), void *)
{
	if (tree != 0) tree->foreach(f);
}

template <class T>
void 
SearchTree<T>::
ins( const T &x, TreeNode<T> *&node )
{
   if ( node == 0 )
      node = new TreeNode<T>( x );
   else
   {
      int cmp = node->compare( x );

      if ( cmp < 0 )
         ins( x, node->left( ) );
      else if ( cmp > 0 )
         ins( x, node->right( ) );
      else
         node->incCount( );
   }
}


template <class T>
void 
SearchTree<T>::
del( TreeNode<T> *&node, TreeNode<T> *&branch )
{
   if ( branch->right( ) != 0 )
      del( node, branch->right( ) );
   else
   {
      node->set( *branch );
      node  = branch;
      branch = branch->left( );
   }
}

template <class T>
void 
SearchTree<T>::
rm( const T &x,TreeNode<T> *&node )
{
   if ( node == 0 )
   {
      return;
   }
   else
   {
      int cmp = node->compare( x );

      if ( cmp < 0 )
         rm( x, node->left( ) );
      else if ( cmp > 0 )
         rm( x, node->right( ) );
      else
      {
         TreeNode<T> *temp = node;
      
         if ( temp->count( ) > 0 )
         {
            temp->decCount( );
            return;
         }

         if ( temp->right( ) == 0 )
            node = temp->left( );
         else if ( temp->left( ) == 0 )
            node = temp->right( );
         else
            del( temp, temp->left( ) );   
         delete temp;
      }
   }
}

// AVL-ausgeglichener sortierter Binaerbaum

template <class T>
void 
AVLTree<T>::
rebalanceLeft( AVLNode<T> *&node )
{
   AVLNode<T> *leftTree = node->left( );

   if ( leftTree->balance( ) == -1 )
   {
      // Einfache LL - Rotation
      
      node->rotateLL( );
      node->setBalance( );

      // Der bisherige linke Unterbaum wird zur neuen Wurzel des
      // Baumes.

      node = leftTree;
   }
   else
   {
      // Zweifache LR - Rotation

      AVLNode<T> *subTree = leftTree->right( );
      node->rotateLRLR( );
      node->setBalance( subTree->balance( ) == -1 ?  1 : 0 );
      leftTree->setBalance( subTree->balance( ) ==  1 ? -1 : 0 );
      
      // Der bisherige rechte Unterbaum des linken Unterbaumes
      // wird zur neuen Wurzel.

      node = subTree;
   }
}

template <class T>
void 
AVLTree<T>::
rebalanceRight( AVLNode<T> *&node )
{
   AVLNode<T> *rightTree = node->right( );

   if ( rightTree->balance( ) == 1 )
   {
      // Einfache RR - Rotation :
      
      node->rotateRR( );
      node->setBalance( );

      // Der bisherige rechte Unterbaum wird zur neuen Wurzel des
      // Baumes.

      node = rightTree;
   }
   else
   {
      // Zweifache RL - Rotation

      AVLNode<T> *subTree = rightTree->left( );
      node->rotateRLRL( );
      node->setBalance( subTree->balance( )  ==  1 ? -1 : 0 );
      rightTree->setBalance( subTree->balance( ) == -1 ?  1 : 0 );
      
      // Der bisherige linke Unterbaum des rechten Unterbaumes
      // wird zur neuen Wurzel.
     
      node = subTree;
   }
} 

template <class T>
void 
AVLTree<T>::
rebalance( AVLNode<T> *&node, int &change, int sign )
{
   if ( node->balance( ) == sign )
   {
      node->setBalance( );
      change = 0;
   }
   else if ( node->balance( ) == 0 )
   {
      node->setBalance( -sign );
   }
   else
   {
      if ( sign == 1 )
         rebalanceLeft( node );
      else 
         rebalanceRight( node );
      node->setBalance( );
      change = 0;
   }
}

#if !defined( NDEBUG )

template <class T>
int 
AVLTree<T>::
isAVL( AVLNode<T> *&node )
{
   if ( node == 0 )
      return 1;
   else
   {
      if ( isAVL( node->left( ) ) && isAVL( node->right( ) ) )
      {
    int hLeft = node->left( ) ? node->left( )->height( ) : 0;
    int hRight = node->right( ) ? node->right( )->height( ) : 0;
         int diff = hLeft - hRight; 
         return int( diff > -2 && diff < 2 );
      }
      else
         return 0;
   }
}

#endif
  
template <class T>
void 
AVLTree<T>::
ins( const T &x, AVLNode<T> *&node, int &change )
{
   if ( node == 0 )
   {
      node = new AVLNode<T>( x );
      change = 1;
   }
   else
   {
      int cmp = node->compare( x );

      if ( cmp < 0 )
      {
         ins( x, node->left( ), change );
         if ( change )
            rebalance( node, change , 1 );
      }
      else if ( cmp > 0 )
      {
         ins( x, node->right( ), change );
         if ( change )
            rebalance( node, change, -1 );
      }
      else
         node->incCount( );
   }
} 

/************************************************************************/
/*                                                                      */
/* void insertAVL( const T & )                                          */
/*                                                                      */
/* Einen Knotem mit einem bestimmten Schluessel in einen AVL-Baum       */
/* einfuegen.                                                           */
/*                                                                      */
/* Parameter:                                                           */
/* x     - Schluessel, der eingefuegt werden soll                       */
/*                                                                      */
/************************************************************************/
 
template <class T>
void 
AVLTree<T>::
insert( const T &x )
{
   int change = 0;

   assert( isAVL( (AVLNode<T>*&)( tree ) ) );
   ins( x, (AVLNode<T>*&)( tree ), change );
}

/************************************************************************/
/*                                                                      */
/* void balanceLeft( AVLNode *&, int & )                               */
/*                                                                      */
/* Die Hoehe des linken Unterbaumes ist um eins geschrumpft. Der Baum   */
/* wird neu ausgeglichen                                                */
/*                                                                      */
/* node   - Verweis auf den Baum, der ausgegleichen werden muss         */
/* change - Verweis uaf das Flag fuer Hoehenaenderungen                 */
/*                                                                      */
/************************************************************************/

template <class T> 
void 
AVLTree<T>::
balanceLeft( AVLNode<T> *&node, int &change )
{
   if ( node->balance( ) == -1 )
      node->setBalance( );
   else if ( node->balance( ) == 0 )
   {
      node->setBalance( 1 );
      change = 0;
   }
   else
   {
      AVLNode<T> *rightTree = node->right( );
      int rightBalance = rightTree->balance( );
      if ( rightBalance >= 0 )
      {
         // RR - Einfachrotation

         node->rotateRR( );

         if ( rightBalance == 0 )
         {
            node->setBalance( 1 );
            rightTree->setBalance( -1 );
            change = 0;
         }
         else
         {
            node->setBalance( );
            rightTree->setBalance( );
         }

         node = rightTree;
      }
      else
      {
         // RL - Doppelrotation

         AVLNode<T> *sub = rightTree->left( );
         int subBalance = sub->balance( );
         
         node->rotateRLRL( );

         node->setBalance( subBalance ==  1 ? -1 : 0 );
         rightTree->setBalance( subBalance == -1 ?  1 : 0 );

         node = sub;
         node->setBalance( );
      }
   }
}

/************************************************************************/
/*                                                                      */
/* Die Hoehe des rechten Unterbaumes ist um eins geschrumpft. Der Baum  */
/* wird neu ausgeglichen                                                */
/*                                                                      */
/************************************************************************/

template <class T>
void
AVLTree<T>::
balanceRight( AVLNode<T> *&node, int &change )
{
   if ( node->balance( ) == 1 )
      node->setBalance( );
   else if ( node->balance( ) == 0 )
   {
      node->setBalance( -1 );
      change = 0;
   }
   else
   {
      AVLNode<T> *leftTree = node->left( );
      int leftBalance = leftTree->balance( );
      if ( leftBalance <= 0 )
      {
          // LL - Einfachrotation

         node->rotateLL( );

         if ( leftBalance == 0 )
         {
            node->setBalance( -1 );
            leftTree->setBalance( 1 );
            change = 0;
         }
         else
         {
            node->setBalance( );
            leftTree->setBalance( );
         }

         node = leftTree;
      }
      else
      {
         // LR - Doppelrotation

         AVLNode<T> *sub = leftTree->right( );
         int subBalance = sub->balance( );

         node->rotateLRLR( );

         node->setBalance( subBalance == -1 ?  1 : 0 );
         leftTree->setBalance( subBalance ==  1 ? -1 : 0 );

         node = sub;
         node->setBalance( );
      }
   }
}

template <class T>
void
AVLTree<T>::
del( AVLNode<T> *&node, AVLNode<T> *&branch, int &change )
{
   if ( branch->right( ) != 0 )
   {
      del( node, branch->right( ), change );
      if ( change )
         balanceRight( branch, change );
   }
   else
   {
      node->set( *branch );
      node = branch;
      branch = branch->left( );
      change = 1;
   }
}

template <class T>
void
AVLTree<T>::
rmNode( AVLNode<T> *&node, int &change )
{
   AVLNode<T> *temp = node;

   // Wenn der Knoten den Schluesselwert mehrfach repraesentiert,
   // Wird nur der Referenzzaehler der Knoten vermindert.

   if ( temp->count( ) > 0 )
   {
      temp->decCount( );
      return;
   }

   // Andernfalls wird der Knoten echt geloescht. Dabei wird er
   // durch seinen linken Teilbaum ersetzt, wenn der rechte 
   // Teilbaum leer ist, bzw. durch seinen rechten Teilbaum,
   // wenn sein linker Teilbaum leer ist.

   if ( temp->right( ) == 0 )
   {
      node = temp->left( );
      change = 1;
   }
   else if ( temp->left( ) == 0 )
   {
      node = temp->right( );
      change = 1;
   }
   else
   {
      del( temp, temp->left( ), change );
      if ( change )
         balanceLeft( node, change );
   }
   delete temp;
}

template <class T>
void 
AVLTree<T>::
rm( const T &x, AVLNode<T> *&node, int &change )
{
   if  ( node == 0 )
   { 
      return;
   }
   else
   {
      // Durch Vergleichen des zu loeschenden Schluessels mit
      // dem Schluessel des aktuellen Knotens den Teilbaum
      // identifizieren, in dem der Schluessel geloescht wird.

      int cmp = node->compare( x );

      if ( cmp < 0 )
      {
         rm( x, node->left( ), change );
         if ( change )
            balanceLeft( node, change );
      }
      else if ( cmp > 0 )
      {
         rm( x, node->right( ), change );
         if ( change )
            balanceRight( node, change );
      }
      else
      {
         // Der Schluesselwert stimmt mit dem aktuellen
         // Knoten ueberein. Der aktuelle Knoten wird entfernt.

         rmNode( node, change );
      }
   }
}

template <class T>
void 
AVLTree<T>::
remove( const T &x )
{
   int change = 0;
   assert( isAVL( (AVLNode<T>*&)( tree ) ) );
   rm( x, (AVLNode<T>*&)( tree ), change );
}

template <class T>
T *
SearchTree<T>::
find(const T &like) const
{
   TreeNode<T> *ptr = tree;
   int dir;

   while (ptr && (dir = ptr->compare(like)) != 0)
      ptr = (dir < 0) ? ptr->left() : ptr->right();
   if (dir == 0 && ptr)
     return &ptr->data();
   else
     return 0;
}







