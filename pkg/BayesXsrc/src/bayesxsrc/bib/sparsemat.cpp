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





#include "sparsemat.h"

double norm(const datamatrix & v)
  {
  double* work = v.getV();
  register unsigned i;
  unsigned size = v.rows()*v.cols();
  double sum2=0;
  for(i=0;i<size;i++,work++)
    sum2+= *work * *work;

  return sqrt(sum2);
  }


double norm (const datamatrix & v, const unsigned col)
  {
  register unsigned i;
  double sum=0;
  for(i=0;i<v.rows();i++)
    sum+= v(i,col)*v(i,col);
  return sqrt(sum);
  }


//------------------------------------------------------------------------------
//------------- class SparseMatrix: implementation of member functions ---------
//------------------------------------------------------------------------------


SparseMatrix::SparseMatrix(const unsigned & row,const unsigned & col,
                           const unsigned & maxnonseros)
  {
  rows=row;
  cols=col;
  values = vector< vector<double> >(rows,vector<double>());
  nonseros = vector< vector<unsigned> >(rows,vector<unsigned>());

  if (maxnonseros > 0)
    {
    unsigned i;
    for(i=0;i<rows;i++)
      {
      values[i].reserve(maxnonseros);
      nonseros[i].reserve(maxnonseros);
      }
    }

  }


SparseMatrix::SparseMatrix(const datamatrix & m,const bool optimize)
  {
  cols = m.cols();
  rows= m.rows();
  values = vector< vector<double> >(rows,vector<double>());
  nonseros = vector< vector<unsigned> >(rows,vector<unsigned>());
  unsigned i,j;

  if (optimize == false)
    {
    unsigned div;
    if (cols <= 20)
      div = 1;
    else
      div = 20;
    for (i=0;i<rows;i++)
      for(j=0;j<cols;j++)
        if (m(i,j) != 0)
          {
          if (values[i].capacity() == 0)
            {
            values[i].reserve(cols/div);
            nonseros[i].reserve(cols/div);
            }
          values[i].push_back(m(i,j));
          nonseros[i].push_back(j);
          }
    } // end: if (optimize == false)
  else
    {    // optimize == true

    unsigned count;
    for (i=0;i<rows;i++)
      {
      count = 0;
      for(j=0;j<cols;j++)
        if (m(i,j) != 0)
          count++;
      values[i].reserve(count);
      nonseros[i].reserve(count);
      }

    for (i=0;i<rows;i++)
      for(j=0;j<cols;j++)
        if (m(i,j) != 0)
          {
          values[i].push_back(m(i,j));
          nonseros[i].push_back(j);
          }

    }   // end: optimize == true

  }   // end: constructor


SparseMatrix::SparseMatrix(const SparseMatrix & m)
  {
  rows = m.rows;
  cols = m.cols;
  values = m.values;
  nonseros = m.nonseros;
  }


const SparseMatrix & SparseMatrix::operator=(const SparseMatrix & m)
  {
  if (this == &m)
    return *this;
  rows = m.rows;
  cols = m.cols;
  values = m.values;
  nonseros = m.nonseros;
  return *this;
  }



double SparseMatrix::operator()(const unsigned & row,const unsigned & col) const
  {
  unsigned i=0;
  while (i<nonseros[row].size())
    {
    if (nonseros[row][i] == col)
      return values[row][i];

    i++;
    }
  return 0.0;
  }


void SparseMatrix::put(const unsigned & row,const unsigned & col,
                       const double & v)
  {

  if (v != 0)
    {
    unsigned i=0;
    bool end=false;
    while ( (end == false) && (i<nonseros[row].size()) )
      {
      if (nonseros[row][i] == col)
        {
        values[row][col] = v;
        end = true;
        }
      else if (nonseros[row][i] > col)
        {
        nonseros[row].insert(nonseros[row].begin()+i,col);
        values[row].insert(values[row].begin()+i,v);
        end = true;
        }

      i++;
      }


    if (end==false)
      {
      values[row].push_back(v);
      nonseros[row].push_back(col);
      }


    }  // end: if v != 0


  }



SparseMatrix SparseMatrix::reorder(const statmatrix<int> & index)
  {
  unsigned i,j,c,k;
  SparseMatrix S(rows,cols);
  double value;
  for (i=0;i<rows;i++)
    {
    for(j=0;j<nonseros[index(i,0)].size();j++)
       {
       value = values[index(i,0)][j];
       c = nonseros[index(i,0)][j];
       for(k=0;k<index.rows();k++)
         if (index(k,0) == c)
           S.put(i,k  ,value);
//       S.put(i,index(nonseros[index(i,0)][j],0),value);
       }

    }

  return S;  
  }


double SparseMatrix::compute_quadform(const datamatrix & x,const unsigned & col)
  {

  double sum=0;
  register unsigned i,j;
  double * workx = x.getV()+col;
  unsigned v = x.cols();
  std::vector< vector<double> >::iterator valuesit = values.begin();
  std::vector< vector<unsigned> >::iterator nonserosit = nonseros.begin();
  for(i=0;i<rows;i++,workx+=v,++nonserosit,++valuesit)
    {
    j=0;
//    while ( (j<nonseros[i].size()) && (nonseros[i][j] <= i) )
      while ( (j< nonserosit->size()) && ((*nonserosit)[j] <= i) )
      {
      if ((*nonserosit)[j] == i)
        sum+= (*valuesit)[j]* *workx * *workx;       //x(i,col)*x(i,col);
      else
        sum+= 2* (*valuesit)[j] *  *workx *x((*nonserosit)[j],col);
      j++;
      }

    }

  return sum;
  }


double SparseMatrix::compute_condmean(const unsigned & i,
                                      const datamatrix & beta)
  {
  double m=0;
  double sumweight=0.0;
  unsigned j;

  std::vector<unsigned>::iterator nonserosit = nonseros[i].begin();
  std::vector<double>::iterator valuesit = values[i].begin();

  unsigned size = nonseros[i].size();
  for (j=0;j<size;j++,++nonserosit,++valuesit)
    {
    if (*nonserosit == i)
      sumweight=*valuesit;
    else
      m-= (*valuesit) * beta(*nonserosit,0);
    }

  return m/sumweight;
  }


SparseMatrix SparseMatrix::kronecker(const SparseMatrix & m) const
  {
  SparseMatrix res(rows*m.rows,cols*m.cols);
  unsigned i,j;
  unsigned k,l;
  unsigned c,r;
  double value;
  for(i=0;i<rows;i++)
    for(j=0;j<nonseros[i].size();j++)
      {
      value = values[i][j];

      for(k=0;k<m.rows;k++)
        for(l=0;l<m.nonseros[k].size();l++)
          {
          r = i*m.rows+k;
          c = nonseros[i][j]*m.cols+m.nonseros[k][l];
          res.values[r].push_back(value*m.values[k][l]);
          res.nonseros[r].push_back(c);

          }

      }
  return res;
  }


datamatrix SparseMatrix::getBlock(const unsigned & rowfirst,
                                  const unsigned & colfirst, const unsigned & rowlast,
                                  const unsigned & collast)

  {
  datamatrix res(rowlast-rowfirst,collast-colfirst);
  unsigned i,j;
  for(i=0;i<res.rows();i++)
    for(j=0;j<res.cols();j++)
      res(i,j) = operator()(rowfirst+i,colfirst+j);

  return res;
  }


SparseMatrix SparseMatrix::getBlockasSparse(const unsigned & rowfirst,
                              const unsigned & colfirst,
                              const unsigned & rowlast,
                              const unsigned & collast)
  {
  SparseMatrix res(rowlast-rowfirst,collast-colfirst);
  unsigned i,j;
  bool stop;

  for(i=rowfirst;i<rowlast;i++)
    {
    j=0;
    stop = false;
    unsigned size;
    while (j<nonseros[i].size())
      {
      if ( (nonseros[i][j] >= colfirst) && (nonseros[i][j] < collast) )
        {
        if (stop==false)
          {
          size = nonseros[i].size()-j+1;
          res.nonseros[i-rowfirst].reserve(size);
          res.values[i-rowfirst].reserve(size);
          stop = true;
          } // end: if (stop==false)

        res.nonseros[i-rowfirst].push_back(nonseros[i][j]-colfirst);
        res.values[i-rowfirst].push_back(values[i][j]);

        } // end: if (nonseros[i][j] >= colfirst)

      j++;
      } // end: while (j<nonseros[i].size())

    }


  return res;
  }



unsigned SparseMatrix::getbandsize(void) const
  {
  unsigned bandsize=0;
  unsigned i;
  unsigned lastcol;
  for(i=0;i<rows;i++)
    {
    lastcol = nonseros[i][nonseros[i].size()-1];
    if (lastcol-i > bandsize)
      bandsize = lastcol-i;
    }

  return bandsize;
  }

void SparseMatrix::print(ostream & o)
  {
  unsigned i,j,k;
  unsigned before;
  for(i=0;i<rows;i++)
    {
    before = 0;
    for(j=0;j<nonseros[i].size();j++)
      {
      for(k=before;k<nonseros[i][j];k++)
        o << 0 << "    ";
      o << values[i][j] << "   ";
      before = nonseros[i][j]+1;
      }

    for(k=before;k<cols;k++)
      o << 0 << "   ";
    o << endl;
    } // end: for(i=0;i<rows;i++)

  }


void SparseMatrix::print2(ostream & o)
  {
  unsigned i,j;
  for(i=0;i<rows;i++)
    {
    for (j=0;j<nonseros[i].size();j++)
      o << values[i][j] << " ";
    o << endl;
    }

  o << endl;

  for(i=0;i<rows;i++)
    {
    for (j=0;j<nonseros[i].size();j++)
      o << nonseros[i][j] << " ";
    o << endl;
    }



  }


void SparseMatrix::mult(const datamatrix & vec,const unsigned & a,
                        const unsigned & col,datamatrix & res)
  {

  assert(col < vec.cols());
  assert(rows == res.rows());

  unsigned i;
  unsigned j;
  std::vector< vector<unsigned> >::iterator nonserosit = nonseros.begin();
  std::vector< vector<double> >::iterator valuesit = values.begin();
  double* reswork = res.getV();

  for (i=0;i<rows;i++,++nonserosit,++valuesit,reswork++)
    {
    *reswork=0;
    for (j=0;j<nonserosit->size();j++)
      *reswork += (*valuesit)[j]*vec(a+(*nonserosit)[j],col);
    }

  }


void SparseMatrix::add_mult(const datamatrix & vec,const unsigned & a,
                            const unsigned & col,datamatrix & res)
  {

  unsigned i;
  unsigned j;
  double* reswork = res.getV();
  std::vector< vector<unsigned> >::iterator nonserosit = nonseros.begin();
  std::vector< vector<double> >::iterator valuesit = values.begin();

  for (i=0;i<rows;i++,reswork++,++nonserosit,++valuesit)
    {
    for (j=0;j<nonserosit->size();j++)
      *reswork += (*valuesit)[j]*vec(a+(*nonserosit)[j],col);

    }

  }


void SparseMatrix::substr_mult(const datamatrix & vec,const unsigned & a,
                               const unsigned & col,datamatrix & res,
                               const unsigned & resrow)
  {

  unsigned i;
  unsigned j;
  std::vector< vector<unsigned> >::iterator nonserosit = nonseros.begin();
  std::vector< vector<double> >::iterator valuesit = values.begin();
  double* reswork = res.getV()+resrow;

  for (i=0;i<rows;i++,++nonserosit,++valuesit,reswork++)
    {

    for (j=0;j<nonserosit->size();j++)
      *reswork -= (*valuesit)[j]*vec(a+(*nonserosit)[j],col);

    }
  }




SparseMatrix Kmrf(const MAP::map & m)
  {

  unsigned S = m.get_nrregions();
  SparseMatrix K(S,S,m.get_maxn());
  unsigned i,j;
  for(i=0;i<S;i++)
    {
    K.put(i,i,m.get_weightssum(i));
    for(j=0;j<m.get_neighbors()[i].size();j++)
      K.put(i,m.get_neighbors()[i][j],-m.get_weights()[i][j]);
    }

  return K;
  }


SparseMatrix Kmrflinear(const unsigned & nr1, const unsigned & nr2)
  {

  SparseMatrix K(nr1*nr2,nr1*nr2,4);

  unsigned i,j;
  unsigned row = 0;
  for(i=0;i<nr1;i++)
    {

    for (j=0;j<nr2;j++)
      {

      if ( (i==0) && (j==0) )
        {
        K.put(0,0,2);
        K.put(0,1,-1);
        K.put(0,nr2,-1);
        }
      else if ( (i==0) && (j==nr2-1) )
        {
        K.put(row,row,2);
        K.put(row,row-1,-1);
        K.put(row,row+nr2,-1);
        }
      else if ( (i==nr1-1) && (j==0) )
        {
        K.put(row,row,2);
        K.put(row,row+1,-1);
        K.put(row,row-nr2,-1);
        }
      else if ( (i==nr1-1) && (j==nr2-1) )
        {
        K.put(row,row,2);
        K.put(row,row-1,-1);
        K.put(row,row-nr2,-1);
        }
      else if ( (i==0) && (j>0) && (j < nr2-1) )
        {
        K.put(row,row,3);
        K.put(row,row-1,-1);
        K.put(row,row+1,-1);
        K.put(row,row+nr2,-1);
        }
      else if ( (i==0) && (j>0) && (j < nr2-1) )
        {
        K.put(row,row,3);
        K.put(row,row-1,-1);
        K.put(row,row+1,-1);
        K.put(row,row+nr2,-1);
        }
      else if ( (i==nr1-1) && (j>0) && (j < nr2-1) )
        {
        K.put(row,row,3);
        K.put(row,row-1,-1);
        K.put(row,row+1,-1);
        K.put(row,row-nr2,-1);
        }
      else if ( (j==0) && (i>0) && (i < nr1-1) )
        {
        K.put(row,row,3);
        K.put(row,row+1,-1);
        K.put(row,row-nr2,-1);
        K.put(row,row+nr2,-1);
        }
      else if ( (j==nr2-1) && (i>0) && (i < nr1-1) )
        {
        K.put(row,row,3);
        K.put(row,row-1,-1);
        K.put(row,row-nr2,-1);
        K.put(row,row+nr2,-1);
        }
      else
        {
        K.put(row,row,4);
        K.put(row,row+1,-1);
        K.put(row,row-1,-1);
        K.put(row,row-nr2,-1);
        K.put(row,row+nr2,-1);
        }


      row++;
      }

    }

  return K;

  }


SparseMatrix Krw1(const vector<double> & weight)
  {

  unsigned S = weight.size();
  datamatrix K(S,S,0);
  unsigned i;
  for (i=1;i<S-1;i++)
    {
    K(i,i) = 1.0/weight[i]+1.0/weight[i+1];
    K(i,i-1) = -1.0/(weight[i]);
    K(i,i+1) = -1.0/(weight[i+1]);
    }
  K(0,0) = 1.0/weight[1];
  K(0,1) = -1.0/weight[1];
  K(S-1,S-1) = 1.0/weight[S-1];
  K(S-1,S-2) = -1.0/weight[S-1];

  return SparseMatrix(K,true);

  }


SparseMatrix Krw2(const vector<double> & weight)
  {
  unsigned i;
  int S = weight.size();

  datamatrix F(S-2,S,0);
  for (i=0;i<F.rows();i++)
    {
    F(i,i)   = weight[2+i]/weight[1+i];
    F(i,i+1) = -(1+weight[2+i]/weight[1+i]);
    F(i,i+2) = 1;
    }
  datamatrix Q(S-2,S-2,0);
  for(i=0;i<Q.rows();i++)
//     Q(i,i) = 1;
	 Q(i,i) =   weight[2+i]*(1+weight[2+i]/weight[1+i]);

  datamatrix K = F.transposed()*Q.inverse()*F;

//  ofstream out("c:\\temp\\K.raw");
//  K.prettyPrint(out);
//  out.close();

  return SparseMatrix(K,true);
  }


SparseMatrix Kseason(unsigned per,unsigned s)
  {

  unsigned k = per-1;
  datamatrix F(s-k,s,0);
  unsigned i,j;
  for(i=0;i<F.rows();i++)
    for(j=i;j<i+per;j++)
      F(i,j) = 1;

  datamatrix Q(s-k,s-k,0);
  for(i=0;i<Q.rows();i++)
    Q(i,i) = 1;

  datamatrix FT = F.transposed();

  datamatrix K;

  K = FT*Q*F;

  return SparseMatrix(K,true);

  }




