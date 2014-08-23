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





#if !defined (VECTORN_CPP_INCLUDED)
#define VECTORN_CPP_INCLUDED
#endif

#if !defined (VECTORN_INCLUDED)
#include "vectorn.h"
#endif

//------------------------------------------------------------------------------
//------------- CLASS vectornum: implementation of member functions ------------
//------------------------------------------------------------------------------


template<class T>
const vectornum<T> & vectornum<T>::operator=(const T v)
  {
  class vectornum<T>::iterator pos;
  for (pos=this->begin(); pos != this->end();++pos)
	 *pos = v;
  return *this;
  }


template<class T>
vectornum<T> operator+(const vectornum<T> & v1,const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2,pos3;
  vectornum<T> h(size());
  for (pos1=begin(),pos2=v.begin(),pos3=h.begin();pos1 != end();
		 ++pos1,++pos2,++pos3)
	 *pos3 = (*pos1)+(*pos2);
*/
  vectornum<T> h(v1.size());
  unsigned i;
  for(i=0;i<v1.size();i++)
    h[i] = v1[i]+v2[i];
  return h;
  }


template<class T>
vectornum<T> operator+(const T & v1,const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v2.size());
  for (pos1=v2.begin(),pos2=h.begin(); pos1 != v2.end();++pos1,++pos2)
	 *pos2 = v1+(*pos1);
*/
  int i;
  vectornum<T> h(v2.size());
  for(i=0;i<v2.size();i++)
    h[i] = v1+v2[i];
  return h;
  }


template<class T>
vectornum<T> operator-(const vectornum<T> & v1,const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2,pos3;
  vectornum<T> h(size());
  for (pos1=begin(),pos2=v.begin(),pos3=h.begin();pos1 != end();
		 ++pos1,++pos2,++pos3)
	 *pos3 = (*pos1)-(*pos2);
*/
  unsigned i;
  vectornum<T> h(v2.size());
  for(i=0;i<v2.size();i++)
    h[i] = v1[i]-v2[i];

  return h;
  }


template<class T>
vectornum<T> operator-(const T & v1,const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v2.size());
  for (pos1=v2.begin(),pos2=h.begin(); pos1 != v2.end();++pos1,++pos2)
	 *pos2 = v1-(*pos1);
*/
  unsigned i;
  vectornum<T> h(v2.size());
  for(i=0;i<v2.size();i++)
    h[i] = v1-v2[i];
  return h;
  }


template<class T>
vectornum<T> operator-(const vectornum<T> & v1,const T & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v1.size());
  for (pos1=v1.begin(),pos2=h.begin(); pos1 != v1.end();++pos1,++pos2)
	 *pos2 = (*pos1)-v2;
 */
  int i;
  vectornum<T> h(v1.size());
  for(i=0;i<v1.size();i++)
    h[i] = v1[i]-v2;
  return h;
  }


template<class T>
vectornum<T> operator*(const vectornum<T> & v1,const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2,pos3;
  vectornum<T> h(v1.size());
  for (pos1=v1.begin(),pos2=v2.begin(),pos3=h.begin();
		 pos1 != v1.end();++pos1,++pos2,++pos3)
	 *pos3 = (*pos1) * (*pos2);
*/
  unsigned i;
  vectornum<T> h(v2.size());
  for(i=0;i<v2.size();i++)
    h[i] = v1[i]*v2[i];
  return h;
  }


template<class T>
vectornum<T> operator*(const T & v1,const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v2.size());
  for (pos1=v2.begin(),pos2=h.begin(); pos1 != v2.end();++pos1,++pos2)
	 *pos2 = v1*(*pos1);
*/
  unsigned i;
  vectornum<T> h(v2.size());
  for(i=0;i<v2.size();i++)
    h[i] = v1*v2[i];
  return h;
  }


template<class T>
vectornum<T> operator/(const vectornum<T> & v1,const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2,pos3;
  vectornum<T> h(size());
  for (pos1=begin(),pos2=v.begin(),pos3=h.begin();pos1 != end();
		  ++pos1,++pos2,++pos3)
	 *pos3 = (*pos1)/(*pos2);
*/
  unsigned  i;
  vectornum<T> h(v2.size());
  for(i=0;i<v2.size();i++)
    h[i] = v1[i]/v2[i];
  return h;
  }


template<class T>
vectornum<T> operator/(const T & v1,const vectornum<T> & v2)
  {
  unsigned i;
  vectornum<T> h(v2.size());
  for(i=0;i<v2.size();i++)
    h[i] = v1/v2[i];
/*
  vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v2.size());
  for (pos1=v2.begin(),pos2=h.begin(); pos1 != v2.end();++pos1,++pos2)
	 *pos2 = v1/(*pos1);
*/
  return h;
  }


template<class T>
vectornum<T> operator>(const vectornum<T> & v1,const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2,pos3;
  vectornum<T> h(v2.size());
  for (pos1=v1.begin(),pos2=v2.begin(),pos3=h.begin();
		 pos1 != v1.end();++pos1,++pos2,++pos3)
	 *pos3 = (*pos1 > *pos2);
*/

  vectornum<T> h(v2.size());
  unsigned long i;
  for(i=0;i<v2.size();i++)
    h[i] = (v1[i] > v2[i]);


  return h;
  }


template<class T>
vectornum<T> operator>(const T & v1,const vectornum<T> & v2)
  {
  class vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v2.size());
  for (pos1=v2.begin(),pos2=h.begin(); pos1 != v2.end();++pos1,++pos2)
	 *pos2 = (v1 > (*pos1));
  return h;
  }


template<class T>
vectornum<T> operator>(const vectornum<T> & v1,const T & v2)
  {
  class vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v1.size());
  for (pos1=v1.begin(),pos2=h.begin(); pos1 != v1.end();++pos1,++pos2)
	 *pos2 = (*pos1 > v2);
  return h;
  }


template<class T>
vectornum<T> operator<(const vectornum<T> & v1,const vectornum<T> & v2)
  {
  return v2 > v1;
  }


template<class T>
vectornum<T> operator<(const vectornum<T> & v1,const T & v2)
  {
  return isgreater(v2,v1);
  }


template<class T>
vectornum<T> operator<(const T & v2,const vectornum<T> & v1)
  {
  return isgreater(v2,v1);
  }


template<class T>
vectornum<T>  vectornum<T>::isequal(vectornum<T> & v2)
  {
  VEC_ITER_TYPE vectornum<T>::iterator pos1,pos2,pos3;
  vectornum<T> h(v2.size());
  for (pos1=this->begin(),pos2=v2.begin(),pos3=h.begin();
		 pos1 != this->end();++pos1,++pos2,++pos3)
	 *pos3 = (*pos1 == *pos2);
  return h;
  }


template<class T>
vectornum<T> isequal(T v1, vectornum<T> & v2)
  {
  class vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v2.size());
  for (pos1=v2.begin(),pos2=h.begin(); pos1 != v2.end();++pos1,++pos2)
	 *pos2 = (v1 == (*pos1));
  return h;
  }


template<class T>
vectornum<T> isequal(vectornum<T> & v1,const T v2)
  {
  class vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v1.size());
  for (pos1=v1.begin(),pos2=h.begin(); pos1 != v1.end();++pos1,++pos2)
	 *pos2 = (*pos1 == v2);
  return h;
  }


template<class T>
vectornum<T> vectornum<T>::isnotequal(vectornum<T> & v2)
  {
  #if defined(__BUILDING_GNU)
  typename vectornum<T>::iterator pos1,pos2,pos3;
  #else
  vectornum<T>::iterator pos1,pos2,pos3;
  #endif
  vectornum<T> h(v2.size());
  for (pos1=this->begin(),pos2=v2.begin(),pos3=h.begin();
		 pos1 != this->end();++pos1,++pos2,++pos3)
	 *pos3 = (*pos1 != *pos2);
  return h;
  }


template<class T>
vectornum<T> isnotequal(T v1, vectornum<T> & v2)
  {
  class vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v2.size());
  for (pos1=v2.begin(),pos2=h.begin(); pos1 != v2.end();++pos1,++pos2)
	 *pos2 = (v1 != (*pos1));
  return h;
  }


template<class T>
vectornum<T> isnotequal(vectornum<T> & v1,const T v2)
  {
  class vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v1.size());
  for (pos1=v1.begin(),pos2=h.begin(); pos1 != v1.end();++pos1,++pos2)
	 *pos2 = (*pos1 != v2);
  return h;
  }


template<class T>
vectornum<T> operator>=(const vectornum<T> & v1,const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2,pos3;
  vectornum<T> h(v2.size());
  for (pos1=v1.begin(),pos2=v2.begin(),pos3=h.begin();
		 pos1 != v1.end();++pos1,++pos2,++pos3)
	 *pos3 = (*pos1 >= *pos2);
*/
  unsigned i;
  vectornum<T> h(v2.size());
  for(i=0;i<v2.size();i++)
    h[i] = (v1[i]>=v2[i]);
  return h;
  }


template<class T>
vectornum<T> operator>=(const T v1, const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v2.size());
  for (pos1=v2.begin(),pos2=h.begin(); pos1 != v2.end();++pos1,++pos2)
	 *pos2 = (v1 >= (*pos1));
*/
  int i;
  vectornum<T> h(v2.size());
  for(i=0;i<v2.size();i++)
    h[i] = (v1>=v2[i]);

  return h;
  }


template<class T>
vectornum<T> operator>=(vectornum<T> & v1,const T v2)
  {

/*
  vectornum<T>::iterator pos1,pos2;
  vectornum<T> h(v1.size());
  for (pos1=v1.begin(),pos2=h.begin(); pos1 != v1.end();++pos1,++pos2)
	 *pos2 = (*pos1 >= v2);
*/
  int i;
  vectornum<T> h(v1.size());
  for(i=0;i<v1.size();i++)
    h[i] = (v1[i]>=v2);
  return h;

  }


template<class T>
vectornum<T> operator<=(const vectornum<T> & v1,const vectornum<T> & v2)
  {
  return v2 >= v1;
  }


template<class T>
vectornum<T> operator<=(vectornum<T> & v1,const T v2)
  {
  return (v2 >= v1);
  }


template<class T>
vectornum<T> operator<=(const T v2, const vectornum<T> & v1)
  {
  return (v2 >= v1);
  }


template<class T>
vectornum<T> operator||(const vectornum<T> & v1,const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2,pos3;
  vectornum<T> h(v2.size());
  for (pos1=begin(),pos2=v2.begin(),pos3=h.begin();
		 pos1 != end();++pos1,++pos2,++pos3)
	 *pos3 = ((*pos1 != 0.0) || (*pos2 != 0.0));
*/
  unsigned i;
  vectornum<T> h(v2.size());
  for(i=0;i<v2.size();i++)
    h[i] = ((v1[i] != 0.0) || (v2[i] != 0.0));
  return h;
  }


template<class T>
vectornum<T> operator&&(const vectornum<T> & v1,const vectornum<T> & v2)
  {
/*
  vectornum<T>::iterator pos1,pos2,pos3;
  vectornum<T> h(v2.size());
  for (pos1=begin(),pos2=v2.begin(),pos3=h.begin();
		 pos1 != end();++pos1,++pos2,++pos3)
	 *pos3 = ((*pos1 != 0.0) && (*pos2 != 0.0));
 */
  unsigned i;
  vectornum<T> h(v2.size());
  for(i=0;i<v2.size();i++)
    h[i] = ((v1[i] != 0.0) && (v2[i] != 0.0));;

  return h;
  }


template<class T>
ostream & operator<<(ostream & out,vectornum<T> & v)
  {
  class vectornum<T>::iterator pos;
  for (pos=v.begin(); pos != v.end();++pos)
	 out << (*pos) << endl;
  return out;
  }


template<class T>
vectornum<T> vectornum<T>::applied(T (*func)(const T &))
  {
  class vectornum<T>::iterator pos,pos2;
  vectornum<T> h(this->size());
  for (pos=this->begin(),pos2=h.begin(); pos != this->end();++pos,++pos2)
	 *pos2 = (func(*pos));
  return h;
  }


template<class T>
vectornum<T> vectornum<T>::applied(T (*func)(T))
  {
  class vectornum<T>::iterator pos;
  class vectornum<T>::iterator pos2;
  vectornum<T> h(this->size());
  for (pos=this->begin(),pos2=h.begin(); pos != this->end();++pos,++pos2)
	 *pos2 = (func(*pos));
  return h;
  }


template<class T>
void vectornum<T>::sort(vector<int> & index,int left,int right)
  {
  int help;
  int i = left;
  int j = right;
  T x = this->operator[](index[(left+right)/2]);
  do
	 {
	 while (this->operator[](index[i]) < x)
		i++;
	 while (x < this->operator[](index[j]))
		j--;
	 if (i <= j)
		{
		help = index[i];
		index[i] = index[j];
		index[j] = help;
		i++;
		j--;
		}
	 }
  while ( i <= j);

  if (left < j)
	 sort(index,left,j);
  if (i < right)
	 sort(index,i,right);
  }


