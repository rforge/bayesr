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





#if !defined(BANDMAT_CPP_INCLUDED)
#define BANDMAT_CPP_INCLUDED
#endif

#include "bandmat.h"

template<class T>
symbandmatrix<T>::symbandmatrix(const unsigned & d,const unsigned & bs,
                                const double & v)
  {
  diagelem = statmatrix<T>(d,1,v);
  D = diagelem;
  decomposed = true;
  if (bs > 0)
    {
    r = statmatrix<T>(d,1);
    z = r;
    upperelem = statmatrix<T>(d,bs,v);
    R = upperelem;
    decomposed = false;
    }
  dim = d;
  bands = bs;
  decomposedonly = false;
  }


template<class T>
symbandmatrix<T>::symbandmatrix(const statmatrix<T> & de,
                                const statmatrix<T> & ud, bool decomp)
  {
  r = statmatrix<T>(de.rows(),1);
  z = r;
  if (!decomp)
    {
    diagelem = de;
    upperelem = ud;
    D = de;
    R = upperelem;
    dim = de.rows();
    bands = ud.cols();
    decomposedonly = false;
    decomposed = false;
    }
  else
    {
    D = de;
    R = ud;
    diagelem = D;
    upperelem = R;
    dim = de.rows();
    bands = ud.cols();
    decomposed = true;
    decomposedonly = true;
    }

  }


template<class T>
symbandmatrix<T>::symbandmatrix(const statmatrix<T> & de)
  {
  diagelem = de;
  dim = de.rows();
  bands = 0;
  decomposedonly = false;
  }


template<class T>
symbandmatrix<T> operator*(const T & v,const symbandmatrix<T> & bm)
  {
  if (bm.decomposedonly == false)
    {
    statmatrix<T> de(bm.dim,1);
    statmatrix<T> ud(bm.dim,bm.bands);
    unsigned i,j;
    for(i=0;i<bm.dim;i++)
      {
      de(i,0) = v*bm.diagelem(i,0);
      for (j=0;j<bm.bands;j++)
        ud(i,j) = v*bm.upperelem(i,j);
      }

    return symbandmatrix<T>(de,ud,false);
    }
  else
    return symbandmatrix<T>();

  }


template<class T>
void symbandmatrix<T>::assign(const statmatrix<T> & de,const statmatrix<T> & ud,
                              bool decomp)
  {

  T * workdiag;
  T * workupper;

  T * workde = de.getV();
  T * workud = ud.getV();

  if (!decomp)
    {
    workdiag = diagelem.getV();
    workupper = upperelem.getV();
    decomposed = false;
    decomposedonly = false;
    }
  else
    {
    workdiag  = D.getV();
    workupper = R.getV();
    decomposed = true;
    decomposedonly = true;
    }

  register unsigned i;
  for(i=0;i<dim;i++,workde++,workdiag++)
    *workdiag = *workde;

  unsigned size = dim*bands;

  for(i=0;i<size;i++,workupper++,workud++)
    *workupper = *workud;

  }


template<class T>
symbandmatrix<T>::symbandmatrix(const symbandmatrix<T> & bm)
  {
  bands = bm.bands;
  dim = bm.dim;
  decomposedonly = bm.decomposedonly;
  decomposed = bm.decomposed;
  diagelem = bm.diagelem;
  upperelem = bm.upperelem;
  D = bm.D;
  R = bm.R;
  r = bm.r;
  z = bm.z;
  det = bm.det;
  }


template<class T>
const symbandmatrix<T> & symbandmatrix<T>::operator=(const symbandmatrix<T> & bm)
  {
  if (this == &bm)
	 return *this;
  bands = bm.bands;
  dim = bm.dim;
  decomposedonly = bm.decomposedonly;
  decomposed = bm.decomposed;
  diagelem = bm.diagelem;
  upperelem = bm.upperelem;
  D = bm.D;
  R = bm.R;
  r = bm.r;
  z = bm.z;
  det = bm.det;
  return * this;
  }


template<class T>
T symbandmatrix<T>::operator()(const unsigned & i, const unsigned & j) const
  {

  assert(i < dim);
  assert(j < dim);
  assert(!decomposedonly);

  if (i==j)
    return diagelem(i,0);
  else if ( (j > i) && (j <= i+bands) )
    return upperelem(i,j-i-1);
  else if ( (j < i) && (i <= j+bands) )
    return upperelem(j,i-j-1);
  else
    return T(0);
  }


template<class T>
void symbandmatrix<T>::set(const unsigned & i, const unsigned & j, const T & t)
  {

  assert(i < dim);
  assert(j < dim);
  assert(!decomposedonly);

  if (i==j)
    {
    diagelem(i,0) = t;
    decomposed = false;
    }
  else if ( (j > i) && (j <= i+bands) )
    {
    upperelem(i,j-i-1) = t;
    decomposed = false;
    }
  else if ( (j < i) && (i <= j+bands) )
    {
    upperelem(j,i-j-1) = t;
    decomposed = false;
    }

  }


template<class T>
void symbandmatrix<T>::print(ostream & out)
  {
  unsigned i,j;

  double * workd;
  double * workup;
  if (!decomposedonly)
    {
    workd = diagelem.getV();
    workup = upperelem.getV();
    }
  else
    {
    workd = D.getV();
    workup = R.getV();
    }

  for(i=0;i<dim;i++,workd++)
    {
    out << *workd << " ";
      for (j=0;j<bands;j++,workup++)
        {
        if (i < dim-j-1)
          out << *workup << " ";
        }

    out << endl;
    }

  }


template<class T>
void symbandmatrix<T>::print2(ostream & out)
  {

  if(!decomposedonly)
    {
    int i,j;
    for(i=0;i<dim;i++)
      {
      for(j=0;j<dim;j++)
        {
        if (j == i)
          out << diagelem(i,0) << " ";
        else if ( (j > i) && (j <= i+bands) )
          out << upperelem(i,j-i-1) << " ";
        else if ( (j < i) &&  (j >= i-int(bands)) )
          out << upperelem(j,i-j-1)  << " ";
        else
          out << 0 << " ";

        }
      out << endl;
      }

    }


  }


template<class T>
T symbandmatrix<T>::compute_quadform(const statmatrix<T> & x,
                                     const unsigned & c) const
  {
  unsigned i;
  unsigned j;
  unsigned pos;
  T sum = T(0);

  assert(decomposedonly==false);

  for(i=0;i<dim;i++)
    {
    sum+= x(i,c)*x(i,c)*diagelem(i,0);
    if (dim-1-i >= bands)
      pos = bands;
    else
      pos = dim-1-i;
    for(j=0;j<pos;j++)
      sum+= 2*x(i,c)*x(i+j+1,c)*upperelem(i,j);
    }

  return sum;
  }


template<class T>
T symbandmatrix<T>::compute_quadformblock(const statmatrix<T> & x,
                            const unsigned & c,const unsigned & a,
                            const unsigned & b)
  {
  unsigned i;
  unsigned j;
  unsigned pos;
  T sum = T(0);
  unsigned size = b-a+1;
  if (!decomposedonly)
    {

    for(i=0;i<size;i++)
      {
      sum+= x(i,c)*x(i,c)*diagelem(a+i,0);
      if (b-a-i >= bands)
        pos = bands;
      else
        pos = b-a-i;
      for(j=0;j<pos;j++)
        sum+= 2*x(i,c)*x(i+j+1,c)*upperelem(a+i,j);
      }

    }

  return sum;
  }


template<class T>
T symbandmatrix<T>::get_det(void)
  {
  unsigned i;
  T logdet = T(0);
  if (!decomposed)
    decomp();

  T * Dp = D.getV();


  if (bands==1 || bands == 2)
    {
    for(i=0;i<dim;i++,Dp++)
      logdet += log(*Dp);
    }
  else
    {
    for(i=0;i<dim;i++,Dp++)
      logdet += log(*Dp);
    logdet *= 2.0;
    }

  return logdet;

  }

template<class T>
void symbandmatrix<T>::decomp(void)
  {
  if (!decomposed)
    {
//    D = statmatrix<T>(dim,1);
//    R = statmatrix<T>(dim,bands);

    if (bands == 1)                                                // triangular
      {
      D(0,0) = diagelem(0,0);                                        // 1.1
      R(0,0) = upperelem(0,0)/D(0,0);                                // 1.2
      register unsigned i;
      for(i=1;i<dim-1;i++)                                           // 1.3
        {
        D(i,0) = diagelem(i,0)-upperelem(i-1,0)*R(i-1,0);            // 1.3.1
        R(i,0) = upperelem(i,0)/D(i,0);                              // 1.3.2
        }
      D(dim-1,0) = diagelem(dim-1,0)-upperelem(dim-2,0)*R(dim-2,0);  // 1.4

/*
      T * Dp = D.getV();
      det = 1;
      for(i=0;i<dim;i++,Dp++)
        {
        det *= *Dp;
        }
*/

      }

    else if (bands == 2)
      {


      T * workD = D.getV();
      T * workdiag = diagelem.getV();
//      D(0,0) = diagelem(0,0);                                        // 1.1
//      R(0,0) = upperelem(0,0)/D(0,0);                                // 1.2
//      R(0,1) = upperelem(0,1)/D(0,0);                                // 1.3

//      D(1,0) = diagelem(1,0)-upperelem(0,0)*R(0,0);                  // 1.4
//      R(1,0) = (upperelem(1,0)-upperelem(0,1)*R(0,0))/D(1,0);        // 1.5
//      R(1,1) = upperelem(1,1)/D(1,0);                                // 1.6

      *workD = *workdiag;                                              // 1.1
      R(0,0) = upperelem(0,0)/(*workD);                                // 1.2
      R(0,1) = upperelem(0,1)/(*workD);                                // 1.3

      workD++;
      workdiag++;

      *workD = *workdiag -upperelem(0,0)*R(0,0);                       // 1.4
      R(1,0) = (upperelem(1,0)-upperelem(0,1)*R(0,0))/(*workD);        // 1.5
      R(1,1) = upperelem(1,1)/(*workD);                                // 1.6

      workD++;
      workdiag++;

      register unsigned i;
      for(i=2;i<dim-2;i++,workD++,workdiag++)
        {
//        D(i,0) = diagelem(i,0) - upperelem(i-2,1)*R(i-2,1) -
//                 D(i-1,0)*R(i-1,0)*R(i-1,0);                         // 1.7.1
//        R(i,0) = (upperelem(i,0)-upperelem(i-1,1)*R(i-1,0))/D(i,0);  // 1.7.2
//        R(i,1) = upperelem(i,1)/D(i,0);                              // 1.7.3

        *workD = *workdiag - upperelem(i-2,1)*R(i-2,1) -
                  (*(workD-1))  *R(i-1,0)*R(i-1,0);                    // 1.7.1
        R(i,0) = (upperelem(i,0)-upperelem(i-1,1)*R(i-1,0))/(*workD);  // 1.7.2
        R(i,1) = upperelem(i,1)/(*workD);                               // 1.7.3

        }

//      D(dim-2,0) = diagelem(dim-2,0) - upperelem(dim-4,1)*R(dim-4,1)
//                   - D(dim-3,0)*R(dim-3,0)*R(dim-3,0);               // 1.8
//      R(dim-2,0) = (upperelem(dim-2,0) -
//                    upperelem(dim-3,1)*R(dim-3,0))/                  // 1.9
//                    D(dim-2,0);
//      D(dim-1,0) = diagelem(dim-1,0) - upperelem(dim-3,1)*R(dim-3,1)
//                   - D(dim-2,0)*R(dim-2,0)*R(dim-2,0);               // 1.10

      *workD = *workdiag  - upperelem(dim-4,1)*R(dim-4,1)
                   - (*(workD-1)) *R(dim-3,0)*R(dim-3,0);               // 1.8
      R(dim-2,0) = (upperelem(dim-2,0) -
                    upperelem(dim-3,1)*R(dim-3,0))/                  // 1.9
                    (*workD);

      workD++;
      workdiag++;

      *workD = *workdiag - upperelem(dim-3,1)*R(dim-3,1)
                   - (*(workD-1)) *R(dim-2,0)*R(dim-2,0);               // 1.10



/*
      T * workD = D.getV();
      T * workdiag = diagelem.getV();
      T * workR0 = R.getV();
      T * workR1 = R.getV()+1;
      T * workupper0 = upperelem.getV();
      T * workupper1 = upperelem.getV()+1;

//      D(0,0) = diagelem(0,0);                                        // 1.1
//      R(0,0) = upperelem(0,0)/D(0,0);                                // 1.2
//      R(0,1) = upperelem(0,1)/D(0,0);                                // 1.3


      *workD = *workdiag;                                              // 1.1
      *workR0 = upperelem(0,0)/(*workD);                               // 1.2
      *workR1 = upperelem(0,1)/(*workD);                               // 1.3

      workD++;
      workdiag++;
      workR0+=2;
      workR1+=2;

//      D(1,0) = diagelem(1,0)-upperelem(0,0)*R(0,0);                  // 1.4
//      R(1,0) = (upperelem(1,0)-upperelem(0,1)*R(0,0))/D(1,0);        // 1.5
//      R(1,1) = upperelem(1,1)/D(1,0);                                // 1.6


      *workD  = *workdiag - (upperelem(0,0) *R(0,0));                  // 1.4

      workupper0+=2;

      *workR0 = upperelem(1,0) - (upperelem(0,1) * R(0,0))/(*workD);       // 1.5

      workupper1+=2;

      *workR1 = upperelem(1,1)/(*workD);                                  // 1.6

      workD++;
      workdiag++;
      workR0+=2;
      workR1+=2;
      workupper0+=2;
      workupper1+=2;

      register unsigned i;
      for(i=2;i<dim-2;i++,workD++,workdiag++,workR0+=2,workR1+=2)
        {
//        D(i,0) = diagelem(i,0) - upperelem(i-2,1)*R(i-2,1) -
//                 D(i-1,0)*R(i-1,0)*R(i-1,0);                         // 1.7.1
//        R(i,0) = (upperelem(i,0)-upperelem(i-1,1)*R(i-1,0))/D(i,0);  // 1.7.2
//        R(i,1) = upperelem(i,1)/D(i,0);                              // 1.7.3

        *workD = *workdiag - upperelem(i-2,1)*R(i-2,1) -
                  (*(workD-1))  * R(i-1,0)*R(i-1,0);                      // 1.7.1
        *workR0 = upperelem(i,0) - (upperelem(i-1,1)*R(i-1,0))/(*workD);  // 1.7.2
        *workR1 = upperelem(i,1)/(*workD);                                // 1.7.3

        }

//      D(dim-2,0) = diagelem(dim-2,0) - upperelem(dim-4,1)*R(dim-4,1)
//                   - D(dim-3,0)*R(dim-3,0)*R(dim-3,0);               // 1.8
//      R(dim-2,0) = (upperelem(dim-2,0) -
//                    upperelem(dim-3,1)*R(dim-3,0))/                  // 1.9
//                    D(dim-2,0);
//      D(dim-1,0) = diagelem(dim-1,0) - upperelem(dim-3,1)*R(dim-3,1)
//                   - D(dim-2,0)*R(dim-2,0)*R(dim-2,0);               // 1.10

      *workD = *workdiag  - upperelem(dim-4,1)*R(dim-4,1)
                   - (*(workD-1)) *R(dim-3,0)*R(dim-3,0);               // 1.8
      *workR0 = upperelem(dim-2,0) -
                    (upperelem(dim-3,1)*R(dim-3,0))/                    // 1.9
                    (*workD);

      workD++;
      workdiag++;

      *workD = *workdiag - upperelem(dim-3,1)*R(dim-3,1)
                   - (*(workD-1)) *R(dim-2,0)*R(dim-2,0);               // 1.10

*/

/*
      T * Dp = D.getV();
      det = 1;
      for(i=0;i<dim;i++,Dp++)
        {
        det *= *Dp;
        }
*/

      } // end: bands == 2

    else
      {      // bands > 2

      register int k,i,j;
      unsigned p;

      double * workD = D.getV();
      double * workdiag = diagelem.getV();
      double * workR = R.getV();
      double * workupper = upperelem.getV();

//      for (k=0;k<dim;k++)
      for (k=0;k<dim;k++,workD++,workdiag++,workR+=bands,workupper+=bands)
        {

//        D(k,0) = sqrt(diagelem(k,0));
        *workD = sqrt(*workdiag);

        if (k+bands < dim)
          p = bands;
        else
          p = dim-1-k;

        for (i=0;i<p;i++)
          {

//          R(k,i) = upperelem(k,i)/D(k,0);
          *(workR+i) = *(workupper+i)/(*workD);

//          if(R(k,i)!=0.0)
          if(*(workR+i)!=0.0)
            {
            for (j = 1 ; j <= i ;j++)
//              upperelem(k+j,i-j) = upperelem(k+j,i-j) - R(k,i) * R(k,j-1);
              upperelem(k+j,i-j) -= *(workR+i) * *(workR+j-1);

//            diagelem(k+i+1,0) = diagelem(k+i+1,0) - R(k,i) * R(k,i);
            *(workdiag+i+1) -= *(workR+i) * *(workR+i);
            }

          }


        }

/*
      T * Dp = D.getV();
      det = 1;
      for(i=0;i<dim;i++,Dp++)
        {
        det *= *Dp * *Dp;
        }
*/

      }  // end: bands > 2

    decomposed = true;
    }

  }


template<class T>
void symbandmatrix<T>::inverse(statmatrix<T> & res)
  {

  if (!decomposed)
    decomp();

   register int i;
   register int j;

   T * workres;

   if (bands==0)
     {
     for (i=0;i<dim;i++)
       for (j=0;j<dim;j++)
         {
         if (i==j)
           res(i,j) = 1/diagelem(i,0);
         else
           res(i,j) = 0;
         }
     }
   else if (bands==1)
     {
     T * workz;
     T * workr;
     T * workD;
     T * workR;

//     T*workres = res.getV() + dim*dim - 1;

     for (i=dim-1;i>=0;i--)
       {

       workz = z.getV()+i;
       workr = r.getV()+i;
       workD = D.getV()+i;
       workR = R.getV()+i;
//       z(i,0) = 1;
//       r(i,0) = 1/D(i,0);
       *workz = 1;
       *workr = 1/(*workD);
       workz++;
       workr++;
       workD++;
       for(j=i+1;j<dim;j++,workz++,workr++,workD++,workR++)
         {
         *workz = - (*workR) * *(workz-1);
         *workr = *workz/(*workD);
//         z(j,0) = -R(j-1,0)*z(j-1,0);
//         r(j,0) = z(j,0)/D(j,0);
         }

//       T * workres = res.getV()+dim-1;

       workr = r.getV()+dim-1;


//       *workres = *workr;
//       res(dim-1,i) = *workres;
//       res(i,dim-1) = *workres;
       res(dim-1,i) = *workr;
//       res(i,dim-1) = res(dim-1,i);
//       res(dim-1,i) = r(dim-1,0);

       workr--;
       workR = R.getV()+dim-2;
//       workres--;
       for(j=dim-2;j>=i;j--,workr--,workR--)
         {
//         *workres = *workr - *workR * *(workres+1);
         res(j,i) = *workr - *workR * res(j+1,i);
//         res(j,i) = r(j,0)-R(j,0)*res(j+1,i);
//         res(i,j) = res(j,i);
//         res(i,j) = *workres;
         }

       }  // end: for (i=dim-1;i>=0;i--)


     }   // end: else if (bands==1)

   else if (bands == 2)
     {

     unsigned beg;

     T * workz;
     T * workr;
     T* workD;
     T * workR1;
     T * workR2;

     for (i=dim-1;i>=0;i--)
       {

       workz = z.getV()+i;
       workr = r.getV()+i;
       workD = D.getV()+i;
       workR1 = R.getV()+2*i;

       // forward

//       z(i,0) = 1;
//       r(i,0) = 1/D(i,0);

       *workz = 1;
       *workr = 1/(*workD);

       if (i==0)
         {

//         z(1,0) = - R(0,0)*z(0,0);                          // 2.2
//         r(1,0) = z(1,0)/D(1,0);                            // 2.4

         workz++;
         workr++;
         workD++;

         *workz = - *workR1  * (*(workz-1));                     // 2.2
         *workr =   (*workz)/(*workD);                        // 2.4
         workz++;
         workr++;
         workD++;
         workR1+=2;
         workR2 = R.getV()+1;

         beg = 2;

         }
       else
         {
         beg = i+1;

//       z(i-1,0) = 0;
//       r(i-1,0) = 0;

         *(workz-1) = 0;
         *(workr-1) = 0;

         workr++;
         workz++;
         workD++;
         workR2 = R.getV()+2*(i-1)+1;
         }

       for (j=beg;j<dim;j++,workr++,workz++,workD++,workR1+=2,workR2+=2)
         {
//         z(j,0) = - R(j-1,0)*z(j-1,0) -
//                  R(j-2,1)*z(j-2,0);                       // 2.3

//         r(j,0) = z(j,0)/D(j,0);                           // 2.4

         *workz = -  *(workR1) * (*(workz-1)) -
                     *(workR2)   * (*(workz-2));                   // 2.3

         *workr = (*workz)/(*workD);                           // 2.4

         }

       // backward

       workr = r.getV()+dim-1;



//       res(dim-1,i) = r(dim-1,0);

       res(dim-1,i) = *workr;

       if (i < dim-1)
         {
         workr--;
         workR1 = R.getV()+2*(dim-2);
//         res(dim-2,i) =r(dim-2,0) - R(dim-2,0)*res(dim-1,i);
         res(dim-2,i) = *workr - *workR1 * res(dim-1,i);
         workr--;

         workR1-=2;
         workR2 = R.getV()+2*(dim-3)+1;

         for (j=dim-3;j>=i;j--,workr--,workR1-=2,workR2-=2)
           {
//           res(j,i) = r(j,0)-R(j,0)*res(j+1,i)-R(j,1)*res(j+2,i);
           res(j,i) = *workr - *workR1 * res(j+1,i)- *workR2 * res(j+2,i);

           }
         }

       } // end: for (i=dim-1;i>=0;i--)

     }  // end: bands == 2

   else
     {

     int p;

     int k;

     double * workr;
     double * workD;

     for (k=dim-1;k>=0;k--)
       {

       // forward

       if (k < bands)
         p = 0;
       else
         p = k-bands;

//       for (j=k-1;j>=p;j--)
//         r(j,0) = 0;

       workr = r.getV() + p;

       for (j=p;j<=k-1;j++,workr++)
         *workr = 0;


//       r(k,0) = 1.0/D(k,0);
       workD = D.getV()+k;
       *workr = 1.0/(*workD);

       workr++;
       workD++;

       for (i=k+1;i<dim;i++,workr++,workD++)
         {
//         r(i,0) = 0;
         *workr = 0;

         if (i < bands)
           p = 0;
         else
           p = i-bands;

         for(j=i-1;j>=p;j--)
//           r(i,0)-= R(j,i-1-j)*r(j,0);
         *workr -= R(j,i-1-j)*r(j,0);

//         r(i,0) = r(i,0)/D(i,0);

         *workr = *workr/(*workD);

         }

       // backward

//       workr = r.getV()+dim-1;
//       workD = D.getV()+dim-1;

       workr--;
       workD--;

       workres= res.getV()+(dim-1)*dim+k;

       for(i=dim-1;i>=k;i--,workr--,workD--,workres-=dim)
         {
//         res(i,k) = r(i,0);
          *workres = *workr;

         if (i+bands > dim-1)
           p = dim-1;
         else
           p = i+bands;

         for (j=i+1;j<=p;j++)
           *workres-= R(i,j-i-1)*res(j,k);
//           res(i,k)-= R(i,j-i-1)*res(j,k);

//         res(i,k) = res(i,k)/D(i,0);

         *workres = *workres/(*workD);

         }


       }  // end: for (k=dim-1;k>=0;k--)

     }  // end: bands > 2

/*
   workres = res.getV()+1;
   for (i=0;i<dim-1;i++,workres+=i+1)
     {
     for(j=i+1;j<dim;j++,workres++)
       *workres = res(j,i);
     }
 */

   }


template<class T>
void symbandmatrix<T>::solveL(const datamatrix & z,datamatrix & res)
  {
  if (!decomposed)
    decomp();

  register int i;

  if (bands==1)
    {
    res(dim-1,0) = z(dim-1,0)/sqrt(D(dim-1,0));

    for(i=dim-2;i>=0;i--)
      res(i,0) = z(i,0)/sqrt(D(i,0)) - R(i,0)*res(i+1,0);
    }

  else if (bands==2)
    {
    res(dim-1,0) = z(dim-1,0)/sqrt(D(dim-1,0));
    res(dim-2,0) = z(dim-2,0)/sqrt(D(dim-2,0)) - R(dim-2,0)*res(dim-1,0);
    for(i=dim-3;i>=0;i--)
      res(i,0) = z(i,0)/sqrt(D(i,0)) - R(i,0)*res(i+1,0)-R(i,1)*res(i+2,0);
    }

  else
    {

    register int j;
    unsigned p;

    double * workres = res.getV()+dim-1;
    double * workz = z.getV()+dim-1;
    double * workD = D.getV()+dim-1;
    double * workR = R.getV()+R.cols()*(dim-1);
    unsigned Rcols = R.cols();

//    for(i=dim-1;i>=0;i--)
    for(i=dim-1;i>=0;i--,--workres,--workz,--workD,workR-=Rcols)
      {
//      res(i,0) = z(i,0);
      *workres = *workz;

      if (i+bands > dim-1)
        p = dim-1;
      else
        p = i+bands;

//      for (j=i+1;j<=p;j++)
      for (j=0;j+i+1<=p;j++)
//        res(i,0)-= R(i,j-i-1)*res(j,0);
        *workres -= *(workR+j) * *(workres+j+1);

//      res(i,0) = res(i,0)/D(i,0);
      *workres /= *workD;

      }

    }

  }


template<class T>
void symbandmatrix<T>::printdecomp(ostream & out)
  {
  if (decomposed)
    {
    unsigned i,j;

    T  * workd;
    T * workup;

    workd = D.getV();
    workup = R.getV();

    for(i=0;i<dim;i++,workd++)
      {
      out << *workd << " ";
        for (j=0;j<bands;j++,workup++)
          {
          if (i < dim-j-1)
            out << *workup << " ";
          }

      out << endl;
      }
    }
  }


template<class T>
void symbandmatrix<T>::solve(const statmatrix<T> & a, statmatrix<T> & res,
                             const unsigned & cola, const unsigned & colres)
  {

  register int i;

  if (!decomposed)
    decomp();

  if (bands == 1)
    {
    // forward

    z(0,0) = a(0,cola);                                // 2.1
    r(0,0) = z(0,0)/D(0,0);                            // 2.3


    for(i=1;i<dim;i++)                                 // 2.2
      {
      z(i,0) = a(i,cola) - R(i-1,0)*z(i-1,0);
      r(i,0) = z(i,0)/D(i,0);                          // 2.3
      }

    // backward

    res(dim-1,colres) = r(dim-1,0);                    // 3.1
    for(i=dim-2;i>=0;i--)                              // 3.2
      res(i,colres) = r(i,0) - R(i,0)*res(i+1,colres);
    }  // end: bands == 1

  else if (bands == 2)
    {

    // forward

    T * workz = z.getV();
    T * workr = r.getV();
    T * workD = D.getV();
    T * workR0 = R.getV();
    T * workR1 = R.getV()+1;
    unsigned asize = a.cols();
    T * worka = a.getV()+cola;


//    z(0,0) = a(0,cola);                                // 2.1
    *workz = *worka;

//    r(0,0) = z(0,0)/D(0,0);                            // 2.4
    *workr = *workz/(*workD);                            // 2.4

    workz++;
    workr++;
    workD++;
    worka+=asize;

//    z(1,0) = a(1,cola)- R(0,0)*z(0,0);                 // 2.2
    *workz = *worka - *workR0  * *(workz-1);             // 2.2


//    r(1,0) = z(1,0)/D(1,0);                            // 2.4
    *workr = *workz/(*workD);                            // 2.4

    workz++;
    workr++;
    workD++;
    workR0+=2;
    worka+=asize;

    for(i=2;i<dim;i++,workz++,workr++,workD++,workR0+=2,workR1+=2,worka+=asize)
      {

//      z(i,0) = a(i,cola) - R(i-1,0)*z(i-1,0) -
//               R(i-2,1)*z(i-2,0);                        // 2.3
      *workz = *worka - *workR0 * *(workz-1) -
                *workR1 * *(workz-2);                      // 2.3

//      r(i,0) = z(i,0)/D(i,0);                            // 2.4
      *workr = *workz/(*workD);                            // 2.4
      }

    // backward

    workr--;
    unsigned sizeres = res.cols();
    T * workres = res.getV()+(dim-1)*sizeres+colres;
    T * workresold = workres;

//    res(dim-1,colres) = r(dim-1,0);                    // 3.1
    *workres = *workr;                                   // 3.1

    workres-=sizeres;
    workr--;
    workR0-=2;

//    res(dim-2,colres) = *workr
//                        - R(dim-2,0)*
//                        res(dim-1,colres);             // 3.2

    *workres = *workr - *workR0 * *workresold ;          // 3.2


    workr--;
    workR0-=2;
    workR1-=2;
    T * workresolder = workresold;
    workresold = workres;
    workres-=sizeres;
    for (i=dim-3;i>=0;i--,workr--,workR0-=2,workR1-=2,workresolder=workresold,
         workresold=workres,workres-=sizeres)
      {
//      res(i,colres) = r(i,0) - R(i,0)*res(i+1,colres)
//                      - R(i,1)*res(i+2,colres);        // 3.3

      *workres = *workr - *workR0 * *workresold
                      - *workR1 * *workresolder;        // 3.3

      }

    } // end: bands == 2

  else
    {

    register int j;
    int p;

    double * workr;
    double * worka;
    double * workD;
    double * workR;
    double * workres;

    unsigned acols = a.cols();
    unsigned rescols = res.cols();
    unsigned Rcols = R.cols();

    // forward

    workr = r.getV();
    worka = a.getV()+cola;
    workD = D.getV();

//    for (i=0;i<dim;i++)
    for (i=0;i<dim;i++,workr++,worka+=acols,workD++)
      {
//      r(i,0) = a(i,cola);
      *workr = *worka;

      if (i < bands)
        p = 0;
      else
        p = i-bands;

      for(j=i-1;j>=p;j--)
//        r(i,0)-= R(j,i-1-j)*r(j,0);
        *workr -= R(j,i-1-j)*r(j,0);

//      r(i,0) = r(i,0)/D(i,0);
      *workr /= *workD;

      }

    // backward

    workr = r.getV()+dim-1;
    workD = D.getV()+dim-1;
    workR = R.getV()+Rcols*(dim-1);
    workres = res.getV()+rescols*(dim-1)+colres;

//    for(i=dim-1;i>=0;i--)
    for(i=dim-1;i>=0;i--,--workr,--workD,workres-=rescols,workR-=Rcols)
      {
//      res(i,colres) = r(i,0);
      *workres = *workr;

      if (i+bands > dim-1)
        p = dim-1;
      else
        p = i+bands;

//      for (j=i+1;j<=p;j++)
      for (j=0;j+i+1<=p;j++)
//        res(i,colres)-= R(i,j-i-1)*res(j,colres);
        *workres -= *(workR+j) * *(workres+(j+1)*rescols);

//      res(i,colres) = res(i,colres)/D(i,0);
      *workres /= *workD;

      }

    } // end: bands > 2

  }


template<class T>
void symbandmatrix<T>::addto(const symbandmatrix<T> & X,
                             const symbandmatrix<T> & K,
                             const T & f1,const T & f2)
  {
  unsigned register i,j;
  T * workK = K.diagelem.getV();
  T * workX = X.diagelem.getV();
  T * workres = diagelem.getV();
  for(i=0;i<dim;i++,workK++,workX++,workres++)
    *workres = f1 * *workX + f2 * *workK;

  workK = K.upperelem.getV();
  workres = upperelem.getV();

  for(i=0;i<dim;i++)
    {
    for (j=0;j<bands;j++,workres++,workK++)
      *workres = f2 * *workK;
    }

  decomposed = false;
  decomposedonly = false;

  }


template<class T>
void symbandmatrix<T>::addtodiag(const symbandmatrix<T> & X,
                                 const symbandmatrix<T> & K, const T & f)
  {
  unsigned register i;
  T * workK = K.diagelem.getV();
  T * workX = X.diagelem.getV();
  T * workres = diagelem.getV();
  for(i=0;i<dim;i++,workK++,workX++,workres++)
    *workres = f * *workX +  *workK;

  decomposed = false;
  decomposedonly = false;

  }



template<class T>
void symbandmatrix<T>::addtoblock(const symbandmatrix<T> & X,
                                  const symbandmatrix<T> & K,
                                  const T & f1,const T & f2,
                                  const unsigned & a,const unsigned & b)
  {
  unsigned register i,j;
  T * workK = K.diagelem.getV()+a;
  T * workX = X.diagelem.getV()+a;
  T * workres = diagelem.getV();
  unsigned size = b-a+1;
  for(i=0;i<size;i++,workK++,workX++,workres++)
    *workres = f1 * *workX + f2 * *workK;

  workres = upperelem.getV();

  for(i=0;i<size;i++)
    {
    for (j=0;j<bands;j++,workres++)
      *workres = f2 * K.upperelem(a+i,j);
    }

  decomposed = false;
  decomposedonly = false;

  }



template<class T>
void symbandmatrix<T>::addtoblock2(const symbandmatrix<T> & X1,
                                  const symbandmatrix<T> & X2,
                                  const T & f1,const T & f2,
                                  const unsigned & a,const unsigned & b)
  {
  unsigned register i,j;
  T * workX1 = X1.diagelem.getV()+a;
  T * workX2 = X2.diagelem.getV()+a;
  T * workres = diagelem.getV();
  unsigned size = b-a+1;
  for(i=0;i<size;i++,workX1++,workX2++,workres++)
    *workres = f1 * *workX1 + f2 * *workX2;

  workX1 = X1.upperelem.getV()+a*X1.bands;
  workX2 = X2.upperelem.getV()+a*X2.bands;
  workres = upperelem.getV();

  unsigned bands1,bands2;

  if(X1.bands > X2.bands)
    {
    bands1 = X2.bands;
    bands2 = X1.bands;

    if(bands1>bands)
      bands1 = bands;
    if(bands2>bands)
      bands2 = bands;

    for(i=0;i<size;i++,workX1+=X1.bands,workX2+=X2.bands)
      {
      for (j=0;j<bands1;j++,workres++)
//        *workres = f1 * X1.upperelem(a+i,j) + f2 * X2.upperelem(a+i,j);
        *workres = f1 * *(workX1+j) + f2 * *(workX2+j);
      for(   ;j<bands2;j++,workres++)
//        *workres = f1 * X1.upperelem(a+i,j);
        *workres = f1 * *(workX1+j);
      }
    }
  else
    {
    bands1 = X1.bands;
    bands2 = X2.bands;

    if(bands1>bands)
      bands1 = bands;
    if(bands2>bands)
      bands2 = bands;

    for(i=0;i<size;i++,workX1+=X1.bands,workX2+=X2.bands)
      {
      for (j=0;j<bands1;j++,workres++)
//        *workres = f2 * X2.upperelem(a+i,j) + f1 * X1.upperelem(a+i,j);
        *workres = f1 * *(workX1+j) + f2 * *(workX2+j);
      for(   ;j<bands2;j++,workres++)
//        *workres = f2 * X2.upperelem(a+i,j);
        *workres = f2 * *(workX2+j);
      }
    }

  decomposed = false;
  decomposedonly = false;

  }



template<class T>
void symbandmatrix<T>::addto2(const symbandmatrix<T> & X,
                              const symbandmatrix<T> & K,
                              const T & f1,const T & f2)
  {
  unsigned register i,j;
  T * workK = K.diagelem.getV();
  T * workX = X.diagelem.getV();
  T * workres = diagelem.getV();
  for(i=0;i<dim;i++,workK++,workX++,workres++)
    *workres = f1 * *workX + f2 * *workK;

  workK = K.upperelem.getV();
  workX = X.upperelem.getV();
  workres = upperelem.getV();

  for(i=0;i<dim;i++)
    {

    for (j=0;j<X.bandsize();j++,workres++,workX++,workK++)
      *workres = f1* *workX + f2 * *workK;

    for (j=X.bandsize();j<K.bandsize();j++,workres++,workK++)
      *workres = f2 * *workK;

    }

  decomposed = false;
  decomposedonly = false;

  }



template<class T>
void symbandmatrix<T>::mult(const statmatrix<T> & X,statmatrix<T> & res) const
  {

  T * workR = res.getV();
  int beg;
  int end;

  register unsigned i,j,k;

  for (i=0;i<dim;i++)
    for (j=0;j<X.cols();j++,workR++)
      {

      *workR=T(0);
      beg = i - bands;
      if (beg < 0)
        beg=0;
      end = i+bands;
      if (end >= dim)
        end = dim-1;

      for (k = beg; k <= end; k++)
        {
        if (k < i)
          *workR += upperelem(k,i-k-1) * X(k,j);
        else if (k == i)
          *workR += diagelem(i,0)* X(k,j);
        else
          *workR += upperelem(i,k-i-1) * X(k,j);
        }

      }

  }


template<class T>
void symbandmatrix<T>::multBlock(const statmatrix<T> & x,statmatrix<T> & res,
                         unsigned rowfirst,unsigned colfirst, unsigned rowlast,
                         unsigned collast,unsigned startx)
  {

  assert(decomposedonly == false);
  assert(rowlast >= rowfirst);
  assert(rowlast < dim);
  assert(collast >= colfirst);
  assert(collast < dim);
  assert(res.rows() == rowlast-rowfirst+1);
  assert(x.rows() >= collast-colfirst+1);

  unsigned i,j;
  T * resp = res.getV();
  for (i=rowfirst;i<=rowlast;i++,resp++)
    {
    if (colfirst > i+bands)
      *resp = 0;
    else if (collast+bands < i)
      *resp = 0;
    else
      {
      unsigned start;
      unsigned end;
      unsigned h;
      if (i >= bands)
        h = i-bands;
      else
        h= 0;


      if (colfirst > h)
        start = colfirst;
      else
        start = h;

      if (collast < i+bands)
        end = collast;
      else
        end = i+bands;

      *resp = 0;
      T * xp = x.getV()+startx+(start-colfirst);
      for(j=start;j<=end;j++,xp++)
        {
        if (j==i)
          *resp += diagelem(j,0) * *xp;
        else if (j > i)
          *resp += upperelem(i,j-i-1) * *xp;
        else
          *resp += upperelem(j,i-j-1) * *xp;
        } // end: for(j=start;j<=end;j++,xp++)


      } // end: else


    }


  }


template<class T>
statmatrix<T> symbandmatrix<T>::getBlock(const unsigned & rowfirst,const unsigned & colfirst,
                    const unsigned & rowlast,const unsigned & collast)
  {

  assert(rowlast >= rowfirst);
  assert(collast >= colfirst);

  int i,j;
  statmatrix<T> block(rowlast-rowfirst+1,collast-colfirst+1);

  double * workblock = block.getV();
/*
  for(i=rowfirst;i<=rowlast;i++)
    for(j=colfirst;j<=collast;j++,workblock++)
      {
      if (i==j)
        *workblock = diagelem(i,0);
      else if ( (j > i) && (j <= i+bands) )
        *workblock = upperelem(i,j-i-1);
      else if ( (j < i) && (i <= j+bands) )
        *workblock = upperelem(j,i-j-1);
      else
        *workblock = T(0);
      }
*/
  for(i=rowfirst;i<=rowlast;i++)
    for(j=colfirst;j<=collast;j++,workblock++)
      *workblock = operator()(i,j);

  return block;
  }


template<class T>
void symbandmatrix<T>::getBlockasband(const unsigned & first,const unsigned & last,symbandmatrix<T> & res)
  {

  assert(decomposedonly==false);
  assert(last >= first);
  assert(last < dim);

  int dimneu = last-first+1;
  int bandsneu;
  if (last-first < bands)
    bandsneu = last - first;
  else
    bandsneu = bands;

  if ( (res.getdim() == dimneu) && (res.bandsize() == bandsneu) )
    { // no new creation

    unsigned i,j;
    T * resdiag = res.getdiagpointer();
    T * diagp = diagelem.getV()+first;
    T * resup = res.getupperpointer();
    T * upperp = upperelem.getV()+bands*first;
    for (i=first;i<=last;i++,resdiag++,diagp++)
      {
      *resdiag = *diagp;
      for(j=0;j<bandsneu;j++,resup++,upperp++)
        {
        *resup = *upperp;
        }

      }

    res.set_decomposed();

    }
  else   // new creation
    {

    datamatrix resd(dimneu,1,0);
    datamatrix resu(dimneu,bandsneu,0);

    unsigned i,j;
    T * resdiag = resd.getV();
    T * diagp = diagelem.getV()+first;
    T * resup = resu.getV();
    T * upperp = upperelem.getV()+bands*first;
    for (i=first;i<=last;i++,resdiag++,diagp++)
      {
      *resdiag = *diagp;
      for(j=0;j<bandsneu;j++,resup++,upperp++)
        {
        *resup = *upperp;
        }

      }

    res = symbandmatrix<T>(resd,resu);

    }

  }


template<class T>
void symbandmatrix<T>::getL(datamatrix & L)
  {
  assert(dim == L.rows());
  assert(dim == L.cols());

  unsigned i,j;

  if (!decomposed)
    decomp();

  for(i=0;i<dim;i++)
    L(i,i) = D(i,0);

  for(i=0;i<bands;i++)
    for(j=0;j<dim-i-1;j++)
      L(j,j+1+i) = R(j,i);

  }


