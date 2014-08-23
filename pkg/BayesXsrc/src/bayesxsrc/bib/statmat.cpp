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





#if !defined(STATMAT_CPP_INCLUDED)
#define STATMAT_CPP_INCLUDED
#endif

#include "statmat.h"

#include <fstream>

//------------------------------------------------------------------------------
//----------- CLASS statmatrix: implementation of member functions -------------
//------------------------------------------------------------------------------

template<class T>
statmatrix<T>::statmatrix(const SparseMatrix & m)
                            : Matrix<T>(m.get_rows(),m.get_cols())
  {
  register unsigned i,j;
  double * work = getV();
  for(i=0;i<this->rows();i++)
    for(j=0;j<this->cols();j++,work++)
      *work = m(i,j);
  }


template<class T>
statmatrix<T>::statmatrix(const vector<T> & v)
                            : Matrix<T>(v.size(),1)
  {
  register unsigned i;
  double * work = this->getV();
  for(i=0;i<this->rows();i++,work++)
    *work = v[i];
  }


template<class T>
const statmatrix<T> & statmatrix<T>::operator=(const SparseMatrix & m)
  {
  statmatrix<T> res(m.get_rows(),m.get_cols());
  double * work = res.getV();
  register unsigned i,j;
  for(i=0;i<res.rows();i++)
    for(j=0;j<res.cols();j++,work++)
      *work = T(m(i,j));
  return res;
  }


// OVERLOADED ASSIGNMENT OPERATOR
	template<class T>
   statmatrix<T>  statmatrix<T>::operator= (const adja& a)
	{
		datamatrix res(a.rows(),a.cols());

		double * work = res.getV();
		register unsigned i,j;
		for(i=0;i<res.rows();i++)
			for(j=0;j<res.cols();j++,work++)
				 *work = a(i,j);
		return res;
	}



template<class T>
void statmatrix<T>::solveroot_t(const statmatrix & b,statmatrix & x)
  {
  int i,j;
  T h;
  T * xp;
  T * xrp = &x(this->rows()-1,0);
  T * bp = b.getV()+this->rows()-1;

  for (i=this->rows()-1;i>=0;i--,xrp--,bp--)
    {
    h=0;
    if (i < this->rows()-1)
      {
      xp = &x(i+1,0);
      for (j=i+1;j<this->rows();j++,xp++)
        h+= (*this)(j,i)* (*xp);
      }
    *xrp = (*bp-h)/((*this)(i,i));
    }
  }


template<class T>
void statmatrix<T>::solveroot(const statmatrix & b,statmatrix & help,
                              statmatrix & x)
  {
  int i,j;
  T h;
  T * mr;
  T * hp;
  T * hrp = help.getV();
  T * bp = b.getV();

  for (i=0;i<this->rows();i++,hrp++,bp++)
    {
    h=0;
    mr = &(*this)(i,0);
    hp = help.getV();
    for (j=0;j<i;j++,mr++,hp++)
      h+= (*mr) * (*hp);
    *hrp = (*bp-h)/((*this)(i,i));

/*
    for (j=0;j<i;j++)
      h+= (*this)(i,j)* help(j,0);
    help(i,0) = (b(i,0)-h)/((*this)(i,i));
*/
    }

  solveroot_t(help,x);

  }


template<class T>
void statmatrix<T>::assign(const statmatrix & A)
  {
  assert(this->rows()==A.rows());
  unsigned size = this->rows( ) * this->cols( );

  T *workA = A.getV( );
  T *workR = this->getV( );
  register unsigned i;
  for ( i = 0;i < size;i++, workA++,workR++ )
    *workR = *workA;
  }


template<class T>
void statmatrix<T>::plus(const statmatrix & A,const statmatrix & B)
  {

  unsigned size = this->rows( ) * this->cols( );

  T *workA = A.getV( );
  T *workB = B.getV( );
  T *workR = this->getV( );
  register unsigned i;

  for ( i = 0;i < size;i++, workA++, workB++, workR++ )
        *workR = *workA + *workB;
  }

template<class T>
void statmatrix<T>::plus(const statmatrix & A)
  {

  unsigned size = this->rows( ) * this->cols( );

  T *workA = A.getV( );
  T *workR = this->getV( );
  register unsigned i;

  for ( i = 0;i < size;i++, workA++, workR++ )
        *workR += *workA;
  }


template<class T>
void statmatrix<T>::minus(const statmatrix & A,const statmatrix & B)
  {
  unsigned size = this->rows( ) * this->cols( );

  T *workA = A.getV( );
  T *workB = B.getV( );
  T *workR = this->getV( );
  register unsigned i;

  for ( i = 0;i < size;i++, workA++, workB++, workR++ )
    *workR = *workA - *workB;
  }


template<class T>
void statmatrix<T>::minus(const statmatrix & A,const statmatrix & B,
                         const unsigned & colA,const unsigned & colB)
  {

  unsigned size = this->rows();
  unsigned sizeA = A.cols();
  unsigned sizeB = B.cols();
  register unsigned i;
  T * workA = A.getV()+colA;
  T * workB = B.getV()+colB;
  T* workR = this->getV();
  for (i=0;i<size;i++,workA+=sizeA,workB+=sizeB,workR++)
    *workR = *workA- *workB;

  }


template<class T>
void statmatrix<T>::mult(const statmatrix & A,const statmatrix & B)
  {

  assert(A.cols() == B.rows());
  assert(this->rows() == A.rows());
  assert(this->cols() == B.cols());

  T * workA;
  T * workB;
  T * workR = this->getV();
  unsigned n = this->cols();
  unsigned size = this->rows()*n;
  register unsigned i, k;

  for (i=0;i <size;i++,workR++)
    {

    *workR=T(0);
    workA = A.getV( ) +  (i / n) * A.cols();
    workB = B.getV( ) +  i % n;

    for (k = 0; k < A.cols(); ++k, ++workA, workB += n )
      if (!(*workA == T(0) || *workB == T(0)))
        *workR += *workA * *workB;

    }
  }


template<class T>
void statmatrix<T>::mult_scalar(const statmatrix & A, const T & b)
  {

  assert(this->rows() == A.rows());
  assert(this->cols() == A.cols());

  unsigned size = A.rows()*A.cols();

  register unsigned i;

  T * workA = A.getV();
  T * workR = this->getV();

  for (i=0;i<size;i++,workR++,workA++)
    *workR = b*(*workA);

  }



template<class T>
void statmatrix<T>::addmult(const statmatrix & A,const statmatrix & B)
  {

  assert(A.cols() == B.rows());
  assert(this->rows() == A.rows());
  assert(this->cols() == B.cols());

  T * workA;
  T * workB;
  T * workR = this->getV();
  unsigned n = this->cols();
  unsigned size = this->rows()*n;
  register unsigned i, k;

  for (i=0;i <size;i++,workR++)
    {

    workA = A.getV( ) +  (i / n) * A.cols();
    workB = B.getV( ) +  i % n;

    for (k = 0; k < A.cols(); ++k, ++workA, workB += n )
      if (!(*workA == T(0) || *workB == T(0)))
        *workR += *workA * *workB;

    }
  }


template<class T>
void statmatrix<T>::addmultsym(const statmatrix & A,const statmatrix & B)
  {

  assert(A.cols() == B.rows());
  assert(this->rows() == A.rows());
  assert(this->cols() == B.cols());

  T * workA;
  T * workB;
  T * workR = this->getV();
  unsigned n = this->cols();
  unsigned m = A.cols();

  register unsigned i,j,k;

  for (i=0;i<this->rows();i++)
    for(j=0;j<n;j++,workR++)
      {

      workA = A.getV( ) +  i*m;
      workB = B.getV( ) +  j;

      for(k=0;k<=i;k++,workB+=n,workA++)
        *workR += * workA * *workB;
//        *workR += A(i,k)*B(j,k);

      workA = A.getV() + (i+1)*m + i;

      for(k=i+1;k<this->rows();k++,workB+=n,workA+=m)
        *workR += *workA * *workB;
//        *workR += A(k,i)*B(j,k);

      }

  }

template<class T>
statmatrix<T> statmatrix<T>::inverse(void)
  {
  assert(this->rows()==this->cols());
  if (this->rows() == 1)
    {
    assert( *getV() != T(0) );
    return statmatrix<T>(1,1,T(1)/(*getV()));
    }
  else if (this->rows()==2)
    {
    T det = this->get(0,0)*this->get(1,1)-this->get(0,1)*this->get(1,0);
    assert(det !=  T(0));
    statmatrix<T> result(2,2);
    T* work = result.getV();
    *work =  this->get(1,1)/det;                 // result(0,0)
    work++;
    *work =  -this->get(0,1)/det;                // result(0,1)
    work++;
    *work =  -this->get(1,0)/det;                // result(1,0)
    work++;
    *work =  this->get(0,0)/det;                 // result(0,0)
    return result;
    }
  else
    return Matrix<T>::inverse();
  }

template<class T>
void statmatrix<T>::multdiagback(const statmatrix & d)
  {
  T* dpoint;
  T* thispoint=getV();
  unsigned i, j;
  for(i=0; i<this->rows(); i++)
    {
    for(j=0, dpoint=d.getV(); j<this->cols(); j++, dpoint++, thispoint++)
      {
      *thispoint *= *dpoint;
      }
    }
  }

template<class T>
void statmatrix<T>::multdiagfront(const statmatrix & d)
  {
  T* dpoint=d.getV();
  T* thispoint=this->getV();
  unsigned i, j;
  for(i=0; i<this->rows(); i++, dpoint++)
    {
    for(j=0; j<this->cols(); j++, thispoint++)
      {
      *thispoint *= *dpoint;
      }
    }
  }

template<class T>
void statmatrix<T>::multdiagfront(const statmatrix & A, const statmatrix & d)
  {
  assert(A.rows()==this->rows());
  assert(A.cols()==this->cols());
  assert(A.rows()==d.rows());
  T* dpoint=d.getV();
  T* apoint=A.getV();
  T* thispoint=this->getV();
  unsigned i, j;
  for(i=0; i<this->rows(); i++, dpoint++)
    {
    for(j=0; j<this->cols(); j++, thispoint++, apoint++)
      {
      *thispoint = *apoint * *dpoint;
      }
    }
  }

template<class T>
void statmatrix<T>::addtodiag(const statmatrix & d, unsigned first,
                              unsigned last)
  {
  assert(this->rows()==this->cols());
  T* dpoint=d.getV();
  T* thispoint=getV()+first*this->cols()+first;
  unsigned i;
  for(i=first; i<last; i++, dpoint++, thispoint+=this->cols()+1)
    {
    *thispoint += *dpoint;
    }
  }

template<class T>
void statmatrix<T>::subfromdiag(const statmatrix & d, unsigned first,
                              unsigned last)
  {
  assert(this->rows()==this->cols());
  T* dpoint=d.getV();
  T* thispoint=this->getV()+first*this->cols()+first;
  unsigned i;
  for(i=first; i<last; i++, dpoint++, thispoint+=this->cols()+1)
    {
    *thispoint -= *dpoint;
    }
  }

template<class T>
void statmatrix<T>::elemmult(const statmatrix<T> & A)
  {
  assert(A.cols()==this->cols());
  assert(A.rows()==this->rows());
  T* Apoint = A.getV();
  T* thispoint = this->getV();
  unsigned i;
  for(i=0; i<this->cols()*this->rows(); i++, thispoint++, Apoint++)
    {
    *thispoint *= *Apoint;
    }
  }

template<class T>
void statmatrix<T>::elemquot(const statmatrix<T> & A)
  {
  assert(A.cols()==this->cols());
  assert(A.rows()==this->rows());
  T* Apoint = A.getV();
  T* thispoint = this->getV();
  unsigned i;
  for(i=0; i<this->cols()*this->rows(); i++, thispoint++, Apoint++)
    {
    *thispoint /= *Apoint;
    }
  }

template<class T>
void statmatrix<T>::weightedsscp(const statmatrix<T> & X,
                                 const statmatrix<T> & w)
  {
  unsigned xcols=X.cols();
  unsigned n=X.rows();

  assert(this->cols()==xcols);
  assert(this->rows()==xcols);
  assert(w.rows()==n);

  T* xpointi;
  T* xpointj;
  T* wpoint;

  register unsigned i,j,k;
  double sum;

  for(i=0; i<xcols; i++)
    {
    for(j=i; j<xcols; j++)
      {
      sum=0;
      xpointi=X.getV()+i;
      xpointj=X.getV()+j;
      wpoint=w.getV();

//      for(k=0; k<n; k++)
      for(k=0; k<n; k++, wpoint++, xpointi+=xcols, xpointj+=xcols)
        {
//        sum += X(k,i)*X(k,j)*w(k,0);
        if(!(*xpointi==T(0)||*xpointj==T(0)))
          {
          sum += *xpointi * *xpointj * *wpoint;
          }
        }
      this->put(i,j,sum);
      if(i!=j)
        {
        this->put(j,i,sum);
        }
      }
    }
  }

template<class T>
void statmatrix<T>::weightedsscp2(const statmatrix<T> & X, const statmatrix<T> & Z,
                                  const statmatrix<T> & w)
  {
  unsigned xcols=X.cols();
  unsigned zcols=Z.cols();
  unsigned n=Z.rows();

  assert(this->cols()==xcols+zcols);
  assert(this->rows()==xcols+zcols);
  assert(w.rows()==n);
  assert(X.rows()==n);

  T* xpointi;
  T* xpointj;
  T* zpointi;
  T* zpointj;
  T* wpoint;

  register unsigned i,i1,j,j1,k;
  double sum;

// compute X'WX, X'WZ and Z'WX
  for(i=0; i<xcols; i++)
    {
    for(j=i; j<xcols; j++)
      {
      sum=0;
      xpointi=X.getV()+i;
      xpointj=X.getV()+j;
      wpoint=w.getV();

//      for(k=0; k<n; k++)
      for(k=0; k<n; k++, wpoint++, xpointi+=xcols, xpointj+=xcols)
        {
//        sum += X(k,i)*X(k,j)*w(k,0);
        if(!(*xpointi==T(0)||*xpointj==T(0)))
          {
          sum += *xpointi * *xpointj * *wpoint;
          }
        }
      this->put(i,j,sum);
      if(i!=j)
        {
        this->put(j,i,sum);
        }
      }
    for(j=0, j1=xcols; j<zcols; j++, j1++)
      {
      sum=0;
      xpointi=X.getV()+i;
      zpointj=Z.getV()+j;
      wpoint=w.getV();

//      for(k=0; k<n; k++)
      for(k=0; k<n; k++, wpoint++, xpointi+=xcols, zpointj+=zcols)
        {
//        sum += X(k,i)*Z(k,j)*w(k,0);
        if(!(*xpointi==T(0)||*zpointj==T(0)))
          {
          sum += *xpointi * *zpointj * *wpoint;
          }
        }
      this->put(i,j1,sum);
      this->put(j1,i,sum);
      }
    }

// compute Z'WZ
  for(i=0, i1=xcols; i<zcols; i++, i1++)
    {
    for(j=i, j1=i+xcols; j<zcols; j++, j1++)
      {
      sum=0;
      zpointi=Z.getV()+i;
      zpointj=Z.getV()+j;
      wpoint=w.getV();

//      for(k=0; k<n; k++)
      for(k=0; k<n; k++, wpoint++, zpointi+=zcols, zpointj+=zcols)
        {
//        sum += Z(k,i)*Z(k,j)*w(k,0);
        if(!(*zpointi==T(0)||*zpointj==T(0)))
          {
          sum += *zpointi * *zpointj * *wpoint;
          }
        }
      this->put(i1,j1,sum);
      if(i!=j)
        {
        this->put(j1,i1,sum);
        }
      }
    }
  }

template<class T>
void statmatrix<T>::weightedsscp_resp(const statmatrix & X, const statmatrix & y,
                          const statmatrix & w)
  {
  unsigned xcols=X.cols();
  unsigned n=X.rows();

  assert(this->rows()==xcols);
  assert(w.rows()==n);
  assert(y.rows()==n);

  register unsigned i,k;
  double sum;

  T* wpoint=w.getV();
  T* ypoint=y.getV();

  statmatrix<T>wy(n,1);
  T* wypoint=wy.getV();

  for(i=0; i<n; i++, ++wpoint, ++ypoint, ++wypoint)
    {
    *wypoint = *ypoint * *wpoint;
    }

  T* xpointi;
  T* thispoint=getV();

  for(i=0; i<xcols; i++, ++thispoint)
    {
    sum=0;
    wypoint=wy.getV();
    xpointi=X.getV()+i;
    for(k=0; k<n; k++, xpointi+=xcols, ++wypoint)
      {
      if(*xpointi!=T(0))
        {
        sum += *xpointi * *wypoint;
        }
      }
    *thispoint = sum;
    }
  }

template<class T>
void statmatrix<T>::weightedsscp_resp2(const statmatrix<T> & X,
                                       const statmatrix<T> & Z,
                                       const statmatrix<T> & y,
                                       const statmatrix<T> & w)
  {
  unsigned xcols=X.cols();
  unsigned zcols=Z.cols();
  unsigned n=Z.rows();

  assert(this->rows()==xcols+zcols);
  assert(w.rows()==n);
  assert(y.rows()==n);
  assert(X.rows()==n);

  register unsigned i,k;
  double sum;

  T* wpoint=w.getV();
  T* ypoint=y.getV();

  statmatrix<T>wy(n,1);
  T* wypoint=wy.getV();

  for(i=0; i<n; i++, ++wpoint, ++ypoint, ++wypoint)
    {
    *wypoint = *ypoint * *wpoint;
    }

  T* xpointi;
  T* zpointi;
  T* thispoint=getV();

// compute X'Wy
  for(i=0; i<xcols; i++, ++thispoint)
    {
    sum=0;
    wypoint=wy.getV();
    xpointi=X.getV()+i;
    for(k=0; k<n; k++, xpointi+=xcols, ++wypoint)
      {
      if(*xpointi!=T(0))
        {
        sum += *xpointi * *wypoint;
        }
      }
    *thispoint = sum;
    }

// compute Z'Wy
  for(i=0; i<zcols; i++, ++thispoint)
    {
    sum=0;
    zpointi=Z.getV()+i;
    wypoint=wy.getV();

    for(k=0; k<n; k++, wypoint++, zpointi+=zcols)
      {
      if(*zpointi!=T(0))
        {
        sum += *zpointi * *wypoint;
        }
      }
    *thispoint = sum;
    }
  }

template<class T>
void statmatrix<T>::sort(int start,int ende,int col)
  {
  int i = start;
  int j = ende;
  T x = this->get((start+ende)/2,col);
  statmatrix<T> hilfe;
  do
	 {
	 while (this->get(i,col) < x)
		i++;
	 while (x < this->get(j,col))
		j--;
	 if (i <= j)
		{
		hilfe = this->getRow(i);
		this->putRow(i,this->getRow(j));
		this->putRow(j,hilfe);
		i++;
		j--;
		}
	 }
  while ( i <= j );
	 if (start < j)
		sort(start,j,col);
	 if (i < ende)
		sort(i,ende,col);
  }


template<class T>
void statmatrix<T>::sortcol(int start,int ende,int col)
  {
  int i = start;
  int j = ende;
  T x = this->get((start+ende)/2,col);
  T hilfe;
  do
	 {
	 while (this->get(i,col) < x)
		i++;
	 while (x < this->get(j,col))
		j--;
	 if (i <= j)
		{
		hilfe = this->get(i,col);
		this->put(i,col,this->get(j,col));
		this->put(j,col,hilfe);
		i++;
		j--;
		}
	 }
  while ( i <= j );
	 if (start < j)
		sortcol(start,j,col);
	 if (i < ende)
		sortcol(i,ende,col);
  }


template<class T>
void statmatrix<T>::indexinit(void)
  {
  unsigned i;
  unsigned j;
  for (i=0;i<this->cols();i++)
	 for (j=0;j<this->rows();j++)
		this->put(j,i,j);
  }


template<class T>
void statmatrix<T>::indexsort(statmatrix<int> & index,int start,int ende,
										int col,int indexcol) const
  {
  int i = start;
  int j = ende;
  T x = this->get(index((start+ende)/2,indexcol),col);
  int hilfe;
  do
	 {
	 while (this->get(index(i,indexcol),col) < x)
		i++;
	 while (x < this->get(index(j,indexcol),col))
		j--;
	 if (i <= j)
		{
		hilfe = index(i,indexcol);
		index(i,indexcol) = index(j,indexcol);
		index(j,indexcol) = hilfe;
		i++;
		j--;
		}
	 }
  while ( i <= j );
	 if (start < j)
		indexsort(index,start,j,col,indexcol);
	 if (i < ende)
		indexsort(index,i,ende,col,indexcol);
  }


template<class T>
void  statmatrix<T>::indexsort2d(statmatrix<int> & index,int start,int ende,
                                 int col,int col2, int indexcol) const
  {
  indexsort(index,start,ende,col,indexcol);

  unsigned j;
  vector<unsigned> posbeg;
  vector<unsigned> posend;
  posbeg.push_back(0);

  double help = this->get(index(0,indexcol),col);
  for(j=1;j<this->rows();j++)
    {
    if (this->get(index(j,indexcol),col) != help)
      {
      posend.push_back(j-1);
      if (j < this->rows())
        posbeg.push_back(j);
      }

    help = this->get(index(j,indexcol),col);

    }

  if (posend.size() < posbeg.size())
    posend.push_back(this->rows()-1);


  for (j=0;j<posbeg.size();j++)
    indexsort(index,posbeg[j],posend[j],col2,indexcol);

  }


template<class T>
void statmatrix<T>::rank(statmatrix<double> & rang,statmatrix<int> & index,
                        int start,int ende,int col) const
  {
  assert(index.rows()==ende+1-start);
  assert(index.cols()==1);
  unsigned j;
  for(j=0;j<rang.rows();j++)
    rang.put(j,0,j+1);

  unsigned i = 1;
  unsigned unten;
  unsigned anzahl = 0;
  double neurang;

  while(i<=ende-start)
    {
    unten = i-1;
    while( (i<=ende-start) && (this->get(index(i,0),col)-this->get(index(i-1,0),col))<pow(10.0,-10.0) )
      {
      anzahl++;
      i++;
      }
    if(anzahl!=0)
      {
      neurang = (rang(unten,0) + rang(unten+anzahl,0)) / 2;
      for(j=unten;j<=unten+anzahl;j++)
      rang.put(j,0,neurang);
      }
    anzahl = 0;
    i++;
    }
  }


template<class T>
statmatrix<T> statmatrix<T>::sum (void) const
  {
  statmatrix<T> s(this->cols(),1,0);
  unsigned col;
  for (col=0;col<this->cols();col++)
    s(col,0) = sum(col);
  return s;
  }


template<class T>
T statmatrix<T>::sum (const unsigned & col) const
  {

  assert(col < this->cols());

  T sum = 0;
  register unsigned i;
  T* work = this->getV()+col;
  for (i=0;i<this->rows();i++,work+=this->cols())
    sum += *work;
  return sum;
  }


template<class T>
T statmatrix<T>::sum2 (const unsigned & col) const
  {

  assert(col < this->cols());

  T sum = 0;
  register unsigned i;
  T* work = this->getV()+col;
  for (i=0;i<this->rows();i++,work+=this->cols())
    sum += *work * *work;
  return sum;
  }


template<class T>
T  statmatrix<T>::sum2(const unsigned & col,const statmatrix<T> & weight) const
  {

  assert(col < this->cols());

  T sum = 0;
  T* work = this->getV()+col;
  T* workweight = weight.getV();
  register unsigned i;
  for (i=0;i<this->rows();i++,work+=this->cols(),workweight++)
    {
    sum += *workweight* (*work) * (*work);
    }

  return sum;

  }

template<class T>
statmatrix<T> statmatrix<T>::sum2()
  {
  statmatrix<T>res(this->cols(),1,0);
  unsigned i;
  for(i=0; i<this->cols(); i++)
    {
    res(i,0)=sum2(i);
    }
  return res;
  }

template<class T>
T statmatrix<T>::mean(const unsigned & col,
                      const statmatrix<T> & weight) const
  {
  assert(col < this->cols());

  T sum = 0;
  T sumweight = 0;
  register unsigned i;
  T* work = this->getV()+col;
  T* workweight = weight.getV();
  for (i=0;i<this->rows();i++,work+=this->cols(),workweight++)
    {
    sumweight+= *workweight;
    sum += *workweight * *work;
    }
  return sum/sumweight;
  }


template<class T>
T statmatrix<T>::min (const unsigned & c) const
  {
  T* work = this->getV()+c;
  T minv = *work;
  work+=this->cols();
  unsigned i;
  for (i=1;i<this->rows();i++,work+=this->cols())
    {
    if ((*work) < minv)
      minv = *work;
    }
  return minv;
  }


template<class T>
T statmatrix<T>::max (const unsigned & c) const
  {
  T* work = this->getV()+c;
  T maxv = *work;
  work+=this->cols();
  unsigned i;
  for (i=1;i<this->rows();i++,work+=this->cols())
    {
    if ((*work) > maxv)
      maxv = *work;
    }
  return maxv;
  }


template<class T>
T statmatrix<T>::sumcomplete(void) const
  {
  register unsigned i;
  unsigned size = this->rows()*this->cols();
  T* work = this->getV();
  T sum = 0;
  for (i=0;i<size;i++,work++)
	 sum += *work;
  return sum;
  }

template<class T>
T statmatrix<T>::norm(unsigned col) const
  {
  T norm=0;
  norm = sqrt(sum2(col));
  return norm;
  }


template<class T>
T statmatrix<T>::euclidean_dist(unsigned col, const statmatrix<T> & A,
const unsigned & colA) const

  {


  assert(col < this->cols());
  assert(colA < A.cols());
  assert(this->rows() == A.rows());

  T sum = 0;
  register unsigned i;
  T* work = this->getV()+col;
  T* workA = A.getV()+colA;
  for (i=0;i<this->rows();i++,work+=this->cols(),workA+=A.cols())
    sum += pow((*work- (*workA)),2)  ;
  return sqrt(sum);
  }


template<class T>
statmatrix<T> statmatrix<T>::norm()
  {
  statmatrix<T>res(this->cols(),1,0);
  unsigned i;
  for(i=0; i<this->cols(); i++)
    {
    res(i,0)=norm(i);
    }
  return res;
  }

template<class T>
statmatrix<T> statmatrix<T>::mean() const
  {
  statmatrix<T> m(this->cols(),1);
  for (unsigned col=0;col<this->cols();col++)
	 m(col,0) = mean(col);
  return m;
  }

template<class T>
T statmatrix<T>::var(const unsigned & col) const
  {
  T m = mean(col);
  return T(1)/T(this->rows())*sum2(col)-m*m;
  }

template<class T>
T statmatrix<T>::var(const unsigned & col,
                     const statmatrix<double> & weight) const
  {
  T m = mean(col,weight);
  T ws = weight.sum(0);
  T s2 = sum2(col,weight);
  return T(1)/ws*s2-m*m;
  }

template<class T>
T statmatrix<T>::quantile(const T & percent,const unsigned & col) const
  {

  T k = this->rows()*(percent/100.0);           // (alpha * Anzahl der Beobachtungen)
  unsigned kganz = unsigned(k);

  statmatrix<int> index(this->rows(),1);
  index.indexinit();
  indexsort(index,0,this->rows()-1,col,0);

  if(kganz == 0)                             // T==0 => Minimum
     return this->get(index(0,0),col);
  else if(kganz == this->rows())                   // T==100 => Maximum
     return this->get(index(this->rows()-1,0),col);
  else if (k == kganz)                             // Falls k ganzzahlig ist
	 return (this->get(index(kganz-1,0),col) + this->get(index(kganz,0),col))/2.0;
  else                                       // Falls k nicht ganzzahlig ist
	 return this->get(index(kganz,0),col);

  }


template<class T>
T statmatrix<T>::quantile(const T & percent,const unsigned & col, statmatrix<int> & index) const
  {

  T k = this->rows()*(percent/100.0);           // (alpha * Anzahl der Beobachtungen)
  unsigned kganz = unsigned(k);

  if(kganz == 0)                             // T==0 => Minimum
     return this->get(index(0,0),col);
  else if(kganz == this->rows())                   // T==100 => Maximum
     return this->get(index(this->rows()-1,0),col);
  else if (k == kganz)                             // Falls k ganzzahlig ist
	 return (this->get(index(kganz-1,0),col) + this->get(index(kganz,0),col))/2.0;
  else                                       // Falls k nicht ganzzahlig ist
	 return this->get(index(kganz,0),col);

  }



template<class T>
statmatrix<T> statmatrix<T>::quantile(T percent)
  {
  statmatrix<T> quant(this->cols(),1);
  for (int col=0;col<this->cols();col++)
	 quant(col,0) = quantile(percent,col);
  return quant;
  }


template<class T>
T statmatrix<T>::autocorr (const unsigned & lag,const unsigned & col) const
  {

  T sum = 0;                               // Summe der Werte
  T sum2 = 0;                              // Quadratsumme der Werte
  T sum_lag = 0;                           // Summe der verzögerten Werte
  T sum_lag2 = 0;                          // Quadratsumme der verzögerten Werte
  T sum_wertlag = 0;                       // Summe Wert * verzögerter Wert
  T mean,mean_lag;                         // Mittelwert, verzögerter Mittelwert
  T anz = this->rows()-lag;                      // Anzahl Beobachtungen

  for (unsigned k=lag;k<this->rows();k++)
	 {
	 sum = sum + this->get(k,col);
	 sum2 = sum2 + this->get(k,col)*this->get(k,col);
	 sum_lag = sum_lag + this->get(k-lag,col);
	 sum_lag2 = sum_lag2 + this->get(k-lag,col)*this->get(k-lag,col);
	 sum_wertlag = sum_wertlag + this->get(k,col)*this->get(k-lag,col);
	 }

  mean = (1.0/anz)*sum;
  mean_lag = (1.0/anz)*sum_lag;

  return	(sum_wertlag - anz*mean*mean_lag)/
			 sqrt( (sum2-anz*mean*mean)*(sum_lag2-anz*mean_lag*mean_lag) );

  }

//------------------------------------------------------------------------------

template<class T>
statmatrix<T> statmatrix<T>::autocorr (const unsigned & lag) const
  {

  statmatrix<T> corr(lag,this->cols());

  for (unsigned i=1;i<=lag;i++)
	 for (unsigned j=0;j<this->cols();j++)
		corr(i-1,j) = autocorr(i,j);

  return corr;
  }

//------------------------------------------------------------------------------


template<class T>
statmatrix<T> statmatrix<T>::autocorr(const unsigned & beginlag,
                                      const unsigned & endlag,
                                      const unsigned & col) const
  {



  unsigned rowstot = endlag-beginlag+1;
  unsigned i;
  statmatrix corr(rowstot,1);


  T sum = 0;                               // Summe der Werte
  T sum2 = 0;                              // Quadratsumme der Werte
  T sum_lag = 0;                           // Summe der verzögerten Werte
  T sum_lag2 = 0;                          // Quadratsumme der verzögerten Werte
  T sum_wertlag = 0;                       // Summe Wert * verzögerter Wert
  T mean,mean_lag;                         // Mittelwert, verzögerter Mittelwert
  T anz = this->rows() - beginlag;                   // Anzahl Beobachtungen

  for (unsigned k=beginlag;k<this->rows();k++)
    {
	 sum = sum + this->get(k,col);
	 sum2 = sum2 + this->get(k,col)*this->get(k,col);
	 sum_lag = sum_lag + this->get(k-beginlag,col);
	 sum_lag2 = sum_lag2 + this->get(k-beginlag,col)*this->get(k-beginlag,col);
	 sum_wertlag = sum_wertlag + this->get(k,col)*this->get(k-beginlag,col);
	 }

  mean = (1.0/anz)*sum;
  mean_lag = (1.0/anz)*sum_lag;

  if ((sum2-anz*mean*mean)*(sum_lag2-anz*mean_lag*mean_lag) <= 0)
    corr(0,0) = 2;
  else
    corr(0,0) = (sum_wertlag - anz*mean*mean_lag)/
			     sqrt( (sum2-anz*mean*mean)*(sum_lag2-anz*mean_lag*mean_lag) );


  for(i=beginlag+1;i<=endlag;i++)
    {

    sum -= this->get(i-1,col);
    sum2 -= this->get(i-1,col)*this->get(i-1,col);
    sum_lag -= this->get(this->rows()-1-(i-1),col);
    sum_lag2 -=  this->get(this->rows()-1-(i-1),col)*this->get(this->rows()-1-(i-1),col);

    sum_wertlag = 0;
    for (unsigned k=i;k<this->rows();k++)
      sum_wertlag += this->get(k,col)*this->get(k-i,col);

    anz--;
    mean = (1.0/anz)*sum;
    mean_lag = (1.0/anz)*sum_lag;


    if ((sum2-anz*mean*mean)*(sum_lag2-anz*mean_lag*mean_lag) <= 0)
      corr(i-beginlag,0) = 2;
    else
      corr(i-beginlag,0) =  (sum_wertlag - anz*mean*mean_lag)/
			     sqrt( (sum2-anz*mean*mean)*(sum_lag2-anz*mean_lag*mean_lag) );

//    corr(i-beginlag,0) = autocorr(i,0);

    }


  return corr;

  }


template<class T>
statmatrix<T> statmatrix<T>::cov()
  {
  statmatrix<T> one(this->rows(),1,1);

  return (1.0/(this->rows()-1))*( (*this).transposed()*(*this) -
			(1.0/this->rows())*(*this).transposed()*one*one.transposed()*(*this) );
  }


template<class T>
statmatrix<T> statmatrix<T>::corr()
  {
  int i,j;
  statmatrix<T> c = cov();
  statmatrix<T> co(this->cols(),this->cols());
  for (i=0;i<c.rows();i++)
	 for(j=0;j<c.cols();j++)
		co(i,j) = c(i,j)/sqrt(c(i,i)*c(j,j));
  return co;
  }


template<class T>
T statmatrix<T>::compute_quadform(const statmatrix<T> & x,const unsigned & c)
  {
  unsigned i,j;
  T res=0;
  T * xi=x.getV()+c;
  T * xj;
  T * workm=this->getV();
  unsigned d = x.cols();
  for (i=0;i<this->rows();i++,xi+=d)
    {
    workm+=i;
//    res+= x(i,0)*x(i,0)*get(i,i);
    res+= *xi * *xi * *workm;
    xj = xi+d;
    workm++;
    for(j=i+1;j<this->cols();j++,xj+=d,workm++)
      {
//      res+=2*x(i,0)*x(j,0)*get(i,j);
      res+=2* *xi * *xj * *workm;

      }

    }

  return res;

  }


template<class T>
statmatrix<T> statmatrix<T>::strike (unsigned int k)
{
	unsigned int i,j;
	unsigned rows_new = this->rows()-1;

	statmatrix<T> matrix_new (rows_new,rows_new);

	if(k==0)
	{
		for(i=0; i<rows_new; i++)
			for(j=0; j<rows_new; j++)
				matrix_new(i,j)=this->get(i+1,j+1);
	}
	else if(k==rows_new+1)
	{
		for(i=0; i<rows_new; i++)
			for(j=0; j<rows_new; j++)
				matrix_new(i,j)=this->get(i,j);
	}

	else
	{
		for(i=0; i<rows_new; i++)
		{
			for(j=0; j<rows_new; j++)
			{
				if(i<k && j<k)
					 matrix_new(i,j)=this->get(i,j);
				else if(i<k && j>k-1)
					matrix_new(i,j)=this->get(i,j+1);
				else if(i>k-1 && j<k)
					matrix_new(i,j)=this->get(i+1,j);
				else if (i>k-1 && j>k-1)
					matrix_new(i,j)=this->get(i+1,j+1);
			}
		}
	}

	return matrix_new;
}

template<class T>
statmatrix<T> statmatrix<T>::get_cov_iX (int i, int j)
{
	assert(this->rows()==this->cols());
	int k,l;
	datamatrix res (1,this->rows()-2);

	l=0;
	for(k=0; k<this->rows(); k++)
	{
		if(k==i)
			l--;
		else if(k==j)
			l--;
		else
			res(0,l) = this->get(i,k);
		l++;
	}




/*	if(i<j)
	{
		for(k=0; k<rows()-2; k++)
		{
			if(k<i)
				res(0,k) = get(i,k);

			else if(k>i-1 && k<j && i!=j-1)
				res(0,k) = get(i,k+1);

			else if( k>i-1  && k<j && i==j-1)
				res(0,k) = get(i,k+2);

			else if(k>j-1)
				res(0,k) = get(i,k+2);

			else
				res(0,k) = get(i,k+2);
		}
	}
	else if (j<i)
	{
		for(k=0; k<rows()-2; k++)
		{
			if(k<j)
				res(0,k) = get(i,k);

			else if(k>j-1 && k<i && j!=i-1)
				res(0,k) = get(i,k+1);

			else if(k>j-1 && k<i && j==i-1)
				res(0,k) = get(i,k+2);

			else if(k>i-1)
				res(0,k) = get(i,k+2);

			else
				res(0,k) = get(i,k+2);
		}
	} */


	return res;
}



template<class T>
statmatrix<T> statmatrix<T>::partial_var(void)
{
	unsigned i,j;
	unsigned nvar = this->cols();


	double numerator, denominator;

	datamatrix cov_all (nvar,nvar);
	datamatrix par_var (nvar,nvar,-999);

	cov_all.assign(cov());

	datamatrix var_iX (nvar-1,nvar-1);
	datamatrix var_jX (nvar-1,nvar-1);
	datamatrix cov_iX (1,nvar-2);
	datamatrix cov_jX (1,nvar-2);
	datamatrix var_X (nvar-2, nvar-2);

	datamatrix help1 (nvar-2,1);
	datamatrix help2 (1, 1);

	datamatrix test1 (1,nvar-2);
	datamatrix test2 (1,1);


	double var_i_X;
	double var_j_X;

	for(i=0; i<nvar; i++)
	{
		for(j=0; j<nvar; j++)
		{
			if(i<j)
			{
				cov_iX.assign(cov_all.get_cov_iX(i,j));
				cov_jX.assign(cov_all.get_cov_iX(j,i));

			//	cout<<cov_iX<<endl;
			//	cout<<cov_jX<<endl;

				datamatrix help (nvar-1, nvar-1);

				var_X.assign((cov_all.strike(i)).strike(j-1));

			//	cout<<var_X<<endl;

				test1.mult(cov_iX,var_X.inverse());
				test2.mult(test1,cov_iX.transposed());

				var_i_X = cov_all(i,i) - test2(0,0);

				test1.mult(cov_jX,var_X.inverse());
				test2.mult(test1,cov_jX.transposed());

				var_j_X = cov_all(j,j) - test2(0,0);

				help1.mult(var_X.inverse(), cov_jX.transposed());
				help2.mult(cov_iX,help1);

				numerator = cov_all(i,j) - help2(0,0);

				denominator = sqrt(var_i_X *var_j_X );

				par_var(i,j)=numerator/denominator;
				par_var(j,i)=numerator/denominator;
			}

			else if (i==j)
				par_var(i,j)=1;
		}
	}

	return par_var;
}

template<class T>
void statmatrix<T>::round(const int digits, const unsigned mincol,
                          const unsigned maxcol, const unsigned minrow,
                          const unsigned maxrow)
{
unsigned cols = maxcol-mincol;
unsigned colinc = this->cols()-cols;
unsigned rows = maxrow-minrow;
unsigned i,j;
double mult = pow(10.0,static_cast<double>(digits));
T* p = this->getV() + mincol + minrow*this->cols();

for(i=0; i<rows; i++, p+=colinc)
  {
  for(j=0; j<cols; j++, p++)
    {
    *p = floor(*p * mult + 0.5) / mult;
    }
  }
}

template<class T>
void statmatrix<T>::round(const int digits)
{
/*unsigned cols = this->cols();
unsigned rows = this->rows();
unsigned i,j;
double mult = pow(10,digits);
T* p = this->getV();
for(i=0; i<rows; i++)
  {
  for(j=0; j<cols; j++, p++)
    {
    *p = floor(*p * mult + 0.5) / mult;
    }
  }*/
this->round(digits, 0, this->cols(), 0, this->rows());
}


template<class T>
bool statmatrix<T>::check_ascending(unsigned & col)
  {
  unsigned i;
  bool asc=true;
  T * p = getV()+col;
  T last = *p;
  i=1;
  while (i<this->rows() && asc==true)
    {
    if (*p < last)
      asc = false;

    last = (*p);
    i++;
    if (i < this->rows()-1)
      p+=this->cols();
    }

  return asc;
  }

/*
template<class T>
statmatrix<T> statmatrix<T>::diag_one(void)
{
	assert(rows()==cols());
	unsigned i,j;
	statmatrix<T> matrix_new (rows(),cols(),0);

	for(i=0; i<rows();i++)
		for(j=0; j<rows();j++)
			matrix_new(i,j) = get(i,j)/ get(i,i);

	return matrix_new;
}


*/

/*statmatrix<double> multdiagback(datamatrix X, const datamatrix & d)
  {
  double* dpoint;
  double* Xpoint=X.getV();
  unsigned i, j;
  for(i=0; i<X.rows(); i++)
    {
    for(j=0, dpoint=d.getV(); j<X.cols(); j++, dpoint++, Xpoint++)
      {
      *Xpoint *= *dpoint;
      }
    }
  return X;
  }*/
template <class T>
statmatrix<T> multdiagback(statmatrix<T> X, const statmatrix<T> & d)
  {
  T* dpoint;
  T* Xpoint=X.getV();
  unsigned i, j;
  for(i=0; i<X.rows(); i++)
    {
    for(j=0, dpoint=d.getV(); j<X.cols(); j++, dpoint++, Xpoint++)
      {
      *Xpoint *= *dpoint;
      }
    }
  return X;
  }

/*statmatrix<double> multdiagfront(datamatrix X, const datamatrix & d)
  {
  double* dpoint=d.getV();
  double* Xpoint=X.getV();
  unsigned i, j;
  for(i=0; i<X.rows(); i++, dpoint++)
    {
    for(j=0; j<X.cols(); j++, Xpoint++)
      {
      *Xpoint *= *dpoint;
      }
    }
  return X;
  }*/
template <class T>
statmatrix<T> multdiagfront(statmatrix<T> X, const statmatrix<T> & d)
  {
  T* dpoint=d.getV();
  T* Xpoint=X.getV();
  unsigned i, j;
  for(i=0; i<X.rows(); i++, dpoint++)
    {
    for(j=0; j<X.cols(); j++, Xpoint++)
      {
      *Xpoint *= *dpoint;
      }
    }
  return X;
  }



