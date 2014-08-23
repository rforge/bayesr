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





#if !defined(ENVMATRIX_CPP_INCLUDED)
#define ENVMATRIX_CPP_INCLUDED
#endif

#include "envmatrix.h"

//-----------------------------------------------------------------------------
//----------------------------- Constructor------------------------------------
//-----------------------------------------------------------------------------

template<class T>
envmatrix<T>::envmatrix(const vector<unsigned> &xe, const T v, const unsigned d)
  {
  xenv=xe;
  dim=d;
  diag=vector<T>(dim,v);
  ldiag=vector<T>(dim,0);
  if(xenv[dim]>0)
    {
    env=vector<T>(xenv[dim],v);
    lenv=vector<T>(xenv[dim],0);
    bandwidth=-1;
    }
  else
    {
    bandwidth=0;
    }
  decomposed=false;
  rational_decomposed=false;
  }


template<class T>
envmatrix<T>::envmatrix(const vector<T> & v, const vector<T> & d,
                        const vector<unsigned> & xe)
  {
  assert(d.size()+1==xe.size());
  assert(v.size()==xe[d.size()]);
  xenv=xe;
  diag=d;
  dim=diag.size();
  env=v;
  ldiag=vector<T>(dim,0);
  lenv=vector<T>(env.size(),0);
  decomposed=false;
  rational_decomposed=false;
  bandwidth=-1;
  }

template<class T>
envmatrix<T>::envmatrix(const statmatrix<T> & X, const double epsilon)
  {
  assert(X.symmetric(epsilon));
  dim=X.cols();
  diag=vector<T>(1,X(0,0));
  xenv=vector<unsigned>(2,0);
  env=vector<T>();
  unsigned c, i, k;
  c=0;
  for(i=1; i<dim; i++)
    {
    diag.push_back(X(i,i));
    k=0;
    for(k; (k<i) && (fabs(X(i,k))<=epsilon); k++)
      {
      }
    for(k; k<i; k++)
      {
      env.push_back(X(i,k));
      c++;
      }
    xenv.push_back(c);
    }
  if(xenv[dim]>0)
    {
    bandwidth=-1;
    lenv=vector<T>(env.size(),0);
    }
  else
    {
    bandwidth=0;
    }
  ldiag=vector<T>(dim,0);
  decomposed=false;
  rational_decomposed=false;
  }

template<class T>
envmatrix<T>::envmatrix(const T &v, const unsigned &d)
  {
  diag=vector<T>(d,v);
  ldiag=vector<T>(d,0);
  dim=d;
  xenv=vector<unsigned>(d+1,0);
  decomposed=false;
  rational_decomposed=false;
  bandwidth=0;
  }

template<class T>
envmatrix<T>::envmatrix(const vector<T> &v, const unsigned &d)
  {
  assert(v.size()==d);
  diag=v;
  ldiag=vector<T>(d,0);
  dim=d;
  xenv=vector<unsigned>(d+1,0);
  decomposed=false;
  rational_decomposed=false;
  bandwidth=0;
  }

template<class T>
envmatrix<T>::envmatrix(const vector<T> & v, const vector<T> & d,
            const vector<unsigned> & xe, const int & bw)
  {
  assert(d.size()+1==xe.size());
  assert(v.size()==xe[d.size()]);
  assert(bw<d.size());
  xenv=xe;
  diag=d;
  dim=diag.size();
  assert(v.size()==bw*dim-(bw+1)*bw/2);
  env=v;
  ldiag=vector<T>(dim,0);
  lenv=vector<T>(env.size(),0);
  decomposed=false;
  rational_decomposed=false;
  bandwidth=bw;
  }

template<class T>
envmatrix<T>::envmatrix(const symbandmatrix<T> & X)
  {
  dim=X.getdim();
  bandwidth=X.bandsize();
  diag=vector<T>(dim,0);
  env=vector<T>(bandwidth*dim-(bandwidth+1)*bandwidth/2,0);
  xenv=vector<unsigned>(dim+1,0);

  T* dx=X.getdiagpointer();
  T* ex=X.getupperpointer();
  VEC_ITER_TYPE vector<T>::iterator di = diag.begin();
  VEC_ITER_TYPE vector<T>::iterator e = env.begin();

  unsigned i, l, k;

//  for(i=0, h=0; i<bandwidth; i++, dx++)
  for(i=0; i<(unsigned)bandwidth; i++, dx++, ++di)
    {
    xenv[i+1]=xenv[i]+i;
//    diag[i]=*dx;
    *di=*dx;
    l=bandwidth-1;
//    for(k=0, ex = X.getupperpointer()+i-1; k<i; k++, ex+=l, h++)
    for(k=0, ex = X.getupperpointer()+i-1; k<i; k++, ex+=l, ++e)
      {
//      env[h]=*ex;
      *e = *ex;
      }
    }

//  for(i=bandwidth; i<dim; i++, dx++)
  for(i=bandwidth; i<dim; i++, dx++, ++di)
    {
    xenv[i+1]=xenv[i]+bandwidth;
//    diag[i]=*dx;
    *di=*dx;
    }

  for(k=0; k<(unsigned)bandwidth; k++)
    {
//    for(i=0, l=bandwidth*(bandwidth-1)/2+k, ex=X.getupperpointer()+(k+1)*bandwidth-k-1;
//     i<dim-bandwidth; l+=bandwidth, ex+=bandwidth, i++)
    for(i=0, e=env.begin()+bandwidth*(bandwidth-1)/2+k,
        ex=X.getupperpointer()+(k+1)*bandwidth-k-1;
        i<dim-bandwidth; ex+=bandwidth, i++, e+=bandwidth)
      {
//      env[l]=*ex;
      *e=*ex;
      }
    }

  ldiag=vector<T>(dim,0);
  lenv=vector<T>(env.size(),0);
  decomposed=false;
  rational_decomposed=false;
  }

template<class T>
envmatrix<T>::envmatrix(const T &v, const unsigned &d, const unsigned bw)
  {
  assert(bw<d);
  dim=d;
  bandwidth=bw;
  diag=vector<T>(dim,v);
  env=vector<T>(bandwidth*dim-(bandwidth+1)*bandwidth/2,v);
  xenv=vector<unsigned>(dim+1,0);
  unsigned i;

  for(i=0; i<bw; i++)
    {
    xenv[i+1]=xenv[i]+i;
    }
  for(i=bandwidth; i<dim; i++)
    {
    xenv[i+1]=xenv[i]+bandwidth;
    }

  ldiag=vector<T>(dim,0);
  lenv=vector<T>(env.size(),0);
  decomposed=false;
  rational_decomposed=false;
  }

template<class T>
T envmatrix<T>::operator()(const unsigned & i, const unsigned & j) const
  {
  assert(i<dim);
  assert(j<dim);

  unsigned ih, jh;
  int kl, ku, zeroes;

  if(i>j)
    {
    ih = i;
    jh = j;
    }
  else if(i<j)
    {
    ih = j;
    jh = i;
    }
  else
    return diag[i];

  kl=xenv[ih];
  ku=xenv[ih+1];
  zeroes=ih-ku+kl;

  if(jh<zeroes)
    return T(0);
  else
    return env[kl+jh-zeroes];

  }


template<class T>
T envmatrix<T>::get(const unsigned & i, const unsigned & j) const
  {
  assert(i<dim);
  assert(j<dim);

  unsigned ih, jh;
  int kl, ku, zeroes;

  if(i>j)
    {
    ih = i;
    jh = j;
    }
  else if(i<j)
    {
    ih = j;
    jh = i;
    }
  else
    return diag[i];

  kl=xenv[ih];
  ku=xenv[ih+1];
  zeroes=ih-ku+kl;

  if(jh<zeroes)
    return T(0);
  else
    return env[kl+jh-zeroes];
  }


template<class T>
envmatrix<T>::envmatrix(const envmatrix & em)
  {
  diag = em.diag;
  env = em.env;
  ldiag = em.ldiag;
  lenv = em.lenv;
  xenv = em.xenv;
  dim = em.dim;
  decomposed = em.decomposed;
  rational_decomposed = em.rational_decomposed;
  bandwidth = em.bandwidth;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

template<class T>
const envmatrix<T> & envmatrix<T>::operator=(const envmatrix<T> & em)
  {
  if (this == &em)
    return *this;

  diag = em.diag;
  env = em.env;
  ldiag = em.ldiag;
  lenv = em.lenv;
  xenv = em.xenv;
  dim = em.dim;
  decomposed = em.decomposed;
  rational_decomposed = em.rational_decomposed;
  bandwidth = em.bandwidth;

  return *this;
  }


//-----------------------------------------------------------------------------
//--- Functions that decompose a matrix or solve systems of linear equations---
//-----------------------------------------------------------------------------

template<class T>
void envmatrix<T>::decomp()
  {
  if(!decomposed)
    {
    if(bandwidth==0)
      {
//      unsigned i;
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
//      for(i=0; i<dim; i++)
      for(; d!=diag.end();++ld, ++d)
        {
//        ldiag[i]=sqrt(diag[i]);
        *ld = sqrt(*d);
        }
      }
    else if(bandwidth==1)
      {
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();

      unsigned i;
//      ldiag[0]=sqrt(diag[0]);
      *ld=sqrt(*d);
      ++d;
//      lenv[0]=env[0]/ldiag[0];
      if(*e!=0)
        {
        *le = *e/ *ld;
        }
      else
        {
        *le=0;
        }
      ++e;
      ++ld;
      for(i=1; i<dim-1; i++, ++e, ++ld, ++d)
        {
//        ldiag[i]=sqrt(diag[i]-lenv[i-1]*lenv[i-1]);
        *ld=sqrt(*d-(*le* *le));
        ++le;
//        lenv[i]=env[i]/ldiag[i];
        if(*e!=0)
          {
          *le=*e/ *ld;
          }
        else
          {
          *le=0;
          }
        }
//      ldiag[dim-1]=sqrt(diag[dim-1]-lenv[dim-2]*lenv[dim-2]);
      *ld=sqrt(*d-*le* *le);
      }
    else if(bandwidth==2)
      {
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();

      unsigned i, h;
//      ldiag[0]=sqrt(diag[0]);
      *ld=sqrt(*d);
//      lenv[0]=env[0]/ldiag[0];
      if(*e!=0)
        {
        *le=*e/ *ld;
        }
      else
        {
        *le=0;
        }
      ++d; ++e, ++ld;
//      ldiag[1]=sqrt(diag[1]-lenv[0]*lenv[0]);
      *ld=sqrt(*d-*le* *le);
      ++le;
      for(i=2, h=1; i<dim; i++)
        {
//        lenv[h]=env[h]/ldiag[i-2];
        if(*e!=0)
          {
          *le=*e/ *(ld-1);
          }
        else
          {
          *le=0;
          }
        h++;
        ++e; ++le;
//        lenv[h]=(env[h]-lenv[h-1]*lenv[h-2])/ldiag[i-1];
        *le=(*e-*(le-1)* *(le-2))/ *ld;
        h++;
        ++e; ++ld; ++d;
//        ldiag[i]=sqrt(diag[i]-lenv[h-1]*lenv[h-1]-lenv[h-2]*lenv[h-2]);
        *ld=sqrt(*d-*le* *le-*(le-1)* *(le-1));
        ++le;
        }
      }
    else if(bandwidth>2)
      {
//      unsigned i, j, k, h, l;
      unsigned i, k, h, l;
      h=0;
      VEC_ITER_TYPE vector<T>::iterator di = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator ldk;
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      VEC_ITER_TYPE vector<T>::iterator lej;
      VEC_ITER_TYPE vector<T>::iterator lel;
      vector<unsigned>::iterator xe;

//      for(i=0; i<bandwidth; i++)
      for(i=0; i<(unsigned)bandwidth; i++, ++di, ++ldi)
        {
//        ldiag[i]=diag[i];
        *ldi=*di;
//        for(k=0; k<i; k++, h++)
        for(k=0, xe = xenv.begin(), ldk = ldiag.begin(); k<i;
            k++, h++, ++e, ++le, ++xe, ++ldk)
          {
//          j=xenv[k];
          lej=lenv.begin()+*xe;
//          lenv[h]=env[h];
          *le=*e;
//          for(l=h-k; l<h; l++, j++)
          for(l=h-k, lel=lenv.begin()+h-k; l<h; l++, ++lel, ++lej)
            {
//            lenv[h]-=lenv[j]*lenv[l];
            *le-=*lej* *lel;
            }
//          lenv[h]/=ldiag[k];
          *le/=*ldk;
//          ldiag[i]-=lenv[h]*lenv[h];
          *ldi-=*le* *le;
          }
//        ldiag[i]=sqrt(ldiag[i]);
        *ldi=sqrt(*ldi);
        }
//      for(i=bandwidth; i<dim; i++)
      for(i=bandwidth; i<dim; i++, ++di, ++ldi)
        {
//        ldiag[i]=diag[i];
        *ldi=*di;
//        for(k=0; k<bandwidth; k++, h++)
        for(k=0, xe = xenv.begin()+i-bandwidth+1, ldk=ldiag.begin()+i-bandwidth;
            k<(unsigned)bandwidth; k++, h++, ++e, ++le, ++xe, ++ldk)
          {
//          lenv[h]=env[h];
          *le=*e;
//          j=xenv[i-bandwidth+k+1]-k;
          lej=lenv.begin()+*xe-k;
//          for(l=h-k; l<h; l++, j++)
          for(l=h-k, lel=lenv.begin()+h-k; l<h; l++, ++lel, ++lej)
            {
//            lenv[h]-=lenv[l]*lenv[j];
            *le-=*lel* *lej;
            }
//          lenv[h]/=ldiag[i-bandwidth+k];
          *le/=*ldk;
//          ldiag[i]-=lenv[h]*lenv[h];
          *ldi-=*le* *le;
          }
//        ldiag[i]=sqrt(ldiag[i]);
        *ldi=sqrt(*ldi);
        }
      }
    else
      {
      unsigned i, ixenv, iband, ifirst, last, k, l, mstop, m, jstop, j;
//        unsigned i, ixenv, iband, ifirst, last, k, mstop, m, jstop, j;
      int kband;
      T temp = T(0);
      T s = T(0);
      T s1 = T(0);

      VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator di = diag.begin();
      vector<unsigned>::iterator xei = xenv.begin()+1;
      vector<unsigned>::iterator xek;
      VEC_ITER_TYPE vector<T>::iterator ek;
      VEC_ITER_TYPE vector<T>::iterator lem;
      VEC_ITER_TYPE vector<T>::iterator lel;
      VEC_ITER_TYPE vector<T>::iterator lek;
      VEC_ITER_TYPE vector<T>::iterator lej;
      VEC_ITER_TYPE vector<T>::iterator ldk;

      assert(*di>0);
      *ldi=sqrt(*di);
      ++di; ++ldi;
//        for(i=1; i<dim; i++)
      for(i=1; i<dim; i++, ++ldi, ++di)
        {
//          ixenv = xenv[i];
//          iband = xenv[i+1]-ixenv;
//          temp = diag[i];

        ixenv = *xei;
        ++xei;
        iband = *xei-ixenv;
        temp = *di;

        if(iband>0)
          {
          ifirst = i-iband;
          last = ixenv;
//            for(k=0; k<iband; k++)
          for(k=0, xek=xenv.begin()+ifirst, ek=env.begin()+ixenv,
              lek=lenv.begin()+ixenv, ldk=ldiag.begin()+ifirst;
              k<iband;
              k++, ++ek, ++lek, ++ldk)
            {
//              kband = xenv[k+ifirst+1]-xenv[k+ifirst];

            kband = -(int)*xek;
            ++xek;
            kband += *xek;

            if(kband > k)
              {
              kband = k;
              }

//              s=env[k+ixenv];
            s = *ek;

            l = k+ixenv-kband;
            lel = lenv.begin()+k+ixenv-kband;
 //           int h1 = *lel;

            if((kband>0) && (last>=l))
//              if((kband>0) && (last>=*lel))
              {
//                mstop=xenv[k+ifirst+1]-1;
//                for(m=xenv[k+ifirst+1]-kband; m<=mstop; m++)
              mstop = *xek - 1;
              for(lem=lenv.begin() + *xek - kband, m=*xek-kband; m<=mstop;
                                 ++lem, ++lel, m++)
                  {
//                  s=s-lenv[m]*lenv[l];
                l++;
                if(*lel!=0 && *lem!=0)
                  {
                  s -= *lel**lem;
                  }
                }
              }
            if(s!=0)
              {
//                lenv[k+ixenv]=s/ldiag[k+ifirst];
              *lek = s/ *ldk;
              last=k+ixenv;
              }
            else
              {
//                lenv[k+ixenv]=0;
              *lek=0;
              }
            }
//            jstop = xenv[i+1]-1;
          jstop = *xei-1;
//           for(j=ixenv; j<=jstop; j++)
         for(j=ixenv, lej=lenv.begin()+ixenv; j<=jstop; j++, ++lej)
            {
//              s1 = lenv[j];
            s1 = *lej;
            if(s1!=0)
              {
              temp = temp - s1*s1;
              }
            }
          }
        assert(temp>0);
//          ldiag[i]=sqrt(temp);
        *ldi = sqrt(temp);
        }
      }
    }
  decomposed=true;
  rational_decomposed=false;
  }



//-----------------------------------------------------------------------------
//--- Functions that decompose a matrix or solve systems of linear equations---
//-----------------------------------------------------------------------------

template<class T>
bool envmatrix<T>::decomp_save()
  {
  bool error = false;
  T help;
  if(!decomposed)
    {
    if(bandwidth==0)
      {
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
      for(; d!=diag.end();++ld, ++d)
        {
        if ((*d) > sqrtmin && (*d) < sqrtmax)
          *ld = sqrt(*d);
        else
          {
          error=true;
          *ld = sqrt(sqrtmin);
          }
        }
      }
    else if(bandwidth==1)
      {
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();

      unsigned i;

      if ((*d) > sqrtmin &&  (*d) < sqrtmax)
        *ld=sqrt(*d);
      else
        {
        error=true;
        *ld = sqrt(sqrtmin);
        }

      ++d;

      if(*e!=0)
        {
        *le = *e/ *ld;
        }
      else
        {
        *le=0;
        }
      ++e;
      ++ld;

      for(i=1; i<dim-1; i++, ++e, ++ld, ++d)
        {
        help = *d-(*le* *le);
        if (help > sqrtmin && help < sqrtmax)
          *ld=sqrt(help);
        else
          {
          error=true;
          *ld = sqrt(sqrtmin);
          }
        ++le;

        if(*e!=0)
          {
          *le=*e/ *ld;
          }
        else
          {
          *le=0;
          }
        }

      help = *d-*le* *le;
      if (help > sqrtmin && help < sqrtmax)
        *ld=sqrt(help);
      else
        {
        error=true;
        *ld = sqrt(sqrtmin);
        }

      }
    else if(bandwidth==2)
      {
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();

      unsigned i, h;

      if ((*d) > sqrtmin && (*d) < sqrtmax)
        *ld=sqrt(*d);
      else
        {
        error=true;
        *ld = sqrt(sqrtmin);
        }

      if(*e!=0)
        {
        *le=*e/ *ld;
        }
      else
        {
        *le=0;
        }
      ++d; ++e, ++ld;

      help = *d-*le* *le;
      if (help > sqrtmin && help < sqrtmax)
        *ld=sqrt(help);
      else
        {
        error=true;
        *ld = sqrt(sqrtmin);
        }

      ++le;
      for(i=2, h=1; i<dim; i++)
        {
        if(*e!=0)
          {
          *le=*e/ *(ld-1);
          }
        else
          {
          *le=0;
          }
        h++;
        ++e; ++le;

        *le=(*e-*(le-1)* *(le-2))/ *ld;
        h++;
        ++e; ++ld; ++d;

        help = *d-*le* *le-*(le-1)* *(le-1);
        if (help > sqrtmin && help < sqrtmax)
          *ld=sqrt(help);
        else
          {
          error=true;
          *ld = sqrt(sqrtmin);
          }

        ++le;
        }
      }
    else if(bandwidth>2)
      {
      unsigned i, k, h, l;
      h=0;
      VEC_ITER_TYPE vector<T>::iterator di = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator ldk;
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      VEC_ITER_TYPE vector<T>::iterator lej;
      VEC_ITER_TYPE vector<T>::iterator lel;
      vector<unsigned>::iterator xe;

      for(i=0; i<(unsigned)bandwidth; i++, ++di, ++ldi)
        {
        *ldi=*di;

        for(k=0, xe = xenv.begin(), ldk = ldiag.begin(); k<i;
            k++, h++, ++e, ++le, ++xe, ++ldk)
          {

          lej=lenv.begin()+*xe;
          *le=*e;
          for(l=h-k, lel=lenv.begin()+h-k; l<h; l++, ++lel, ++lej)
            {
            *le-=*lej* *lel;
            }
          *le/=*ldk;
          *ldi-=*le* *le;
          }

        if ((*ldi) > sqrtmin && (*ldi) < sqrtmax)
          *ldi=sqrt(*ldi);
        else
          {
          error=true;
          *ldi = sqrt(sqrtmin);
          }

        }

      for(i=bandwidth; i<dim; i++, ++di, ++ldi)
        {

        *ldi=*di;

        for(k=0, xe = xenv.begin()+i-bandwidth+1, ldk=ldiag.begin()+i-bandwidth;
            k<(unsigned)bandwidth; k++, h++, ++e, ++le, ++xe, ++ldk)
          {

          *le=*e;
          lej=lenv.begin()+*xe-k;
          for(l=h-k, lel=lenv.begin()+h-k; l<h; l++, ++lel, ++lej)
            {
            *le-=*lel* *lej;
            }

          *le/=*ldk;

          *ldi-=*le* *le;
          }


        if ((*ldi) > sqrtmin && (*ldi) < sqrtmax)
          *ldi=sqrt(*ldi);
        else
          {
          error=true;
          *ldi = sqrt(sqrtmin);
          }

        }
      }
    else
      {
      unsigned i, ixenv, iband, ifirst, last, k, l, mstop, m, jstop, j;

      int kband;
      T temp = T(0);
      T s = T(0);
      T s1 = T(0);

      VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator di = diag.begin();
      vector<unsigned>::iterator xei = xenv.begin()+1;
      vector<unsigned>::iterator xek;
      VEC_ITER_TYPE vector<T>::iterator ek;
      VEC_ITER_TYPE vector<T>::iterator lem;
      VEC_ITER_TYPE vector<T>::iterator lel;
      VEC_ITER_TYPE vector<T>::iterator lek;
      VEC_ITER_TYPE vector<T>::iterator lej;
      VEC_ITER_TYPE vector<T>::iterator ldk;


      if ((*di) > sqrtmin && (*di) < sqrtmax)
          *ldi=sqrt(*di);
      else
        {
        error=true;
        *ldi = sqrt(sqrtmin);
        }

      ++di; ++ldi;

      for(i=1; i<dim; i++, ++ldi, ++di)
        {

        ixenv = *xei;
        ++xei;
        iband = *xei-ixenv;
        temp = *di;

        if(iband>0)
          {
          ifirst = i-iband;
          last = ixenv;

          for(k=0, xek=xenv.begin()+ifirst, ek=env.begin()+ixenv,
              lek=lenv.begin()+ixenv, ldk=ldiag.begin()+ifirst;
              k<iband;
              k++, ++ek, ++lek, ++ldk)
            {

            kband = -(int)*xek;
            ++xek;
            kband += *xek;

            if(kband > k)
              {
              kband = k;
              }

            s = *ek;

            l = k+ixenv-kband;
            lel = lenv.begin()+k+ixenv-kband;

            if((kband>0) && (last>=l))
              {
              mstop = *xek - 1;
              for(lem=lenv.begin() + *xek - kband, m=*xek-kband; m<=mstop;
                                 ++lem, ++lel, m++)
                  {
                l++;
                if(*lel!=0 && *lem!=0)
                  {
                  s -= *lel**lem;
                  }
                }
              }
            if(s!=0)
              {

              *lek = s/ *ldk;
              last=k+ixenv;
              }
            else
              {

              *lek=0;
              }
            }

          jstop = *xei-1;

         for(j=ixenv, lej=lenv.begin()+ixenv; j<=jstop; j++, ++lej)
            {
            s1 = *lej;
            if(s1!=0)
              {
              temp = temp - s1*s1;
              }
            }
          }

        if (temp > sqrtmin && temp < sqrtmax)
          *ldi = sqrt(temp);
        else
          {
          error=true;
          *ldi = sqrt(sqrtmin);
          }

        }
      }
    }
  if (error == false)
    {
    decomposed=true;
    rational_decomposed=false;
    }
  return error;

  }


template<class T>
void envmatrix<T>::decomp2(unsigned start)
  {
  if(!decomposed)
    {
    if(bandwidth==0)
      {
//      unsigned i;
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
//      for(i=0; i<dim; i++)
      for(; d!=diag.end();++ld, ++d)
        {
//        ldiag[i]=sqrt(diag[i]);
        *ld = sqrt(*d);
        }
      }
    else if(bandwidth==1)
      {
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();

      unsigned i;
//      ldiag[0]=sqrt(diag[0]);
      *ld=sqrt(*d);
      ++d;
//      lenv[0]=env[0]/ldiag[0];
      if(*e!=0)
        {
        *le = *e/ *ld;
        }
      else
        {
        *le=0;
        }
      ++e;
      ++ld;
      for(i=1; i<dim-1; i++, ++e, ++ld, ++d)
        {
//        ldiag[i]=sqrt(diag[i]-lenv[i-1]*lenv[i-1]);
        *ld=sqrt(*d-(*le* *le));
        ++le;
//        lenv[i]=env[i]/ldiag[i];
        if(*e!=0)
          {
          *le=*e/ *ld;
          }
        else
          {
          *le=0;
          }
        }
//      ldiag[dim-1]=sqrt(diag[dim-1]-lenv[dim-2]*lenv[dim-2]);
      *ld=sqrt(*d-*le* *le);
      }
    else if(bandwidth==2)
      {
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();

      unsigned i, h;
//      ldiag[0]=sqrt(diag[0]);
      *ld=sqrt(*d);
//      lenv[0]=env[0]/ldiag[0];
      if(*e!=0)
        {
        *le=*e/ *ld;
        }
      else
        {
        *le=0;
        }
      ++d; ++e, ++ld;
//      ldiag[1]=sqrt(diag[1]-lenv[0]*lenv[0]);
      *ld=sqrt(*d-*le* *le);
      ++le;
      for(i=2, h=1; i<dim; i++)
        {
//        lenv[h]=env[h]/ldiag[i-2];
        if(*e!=0)
          {
          *le=*e/ *(ld-1);
          }
        else
          {
          *le=0;
          }
        h++;
        ++e; ++le;
//        lenv[h]=(env[h]-lenv[h-1]*lenv[h-2])/ldiag[i-1];
        *le=(*e-*(le-1)* *(le-2))/ *ld;
        h++;
        ++e; ++ld; ++d;
//        ldiag[i]=sqrt(diag[i]-lenv[h-1]*lenv[h-1]-lenv[h-2]*lenv[h-2]);
        *ld=sqrt(*d-*le* *le-*(le-1)* *(le-1));
        ++le;
        }
      }
    else if(bandwidth>2)
      {
//      unsigned i, j, k, h, l;
      unsigned i, k, h, l;

      h = unsigned(start*(start-1)/2);
      if(start > (unsigned)bandwidth)
        h = unsigned(bandwidth*(bandwidth-1)/2 + bandwidth*(start-bandwidth));

      VEC_ITER_TYPE vector<T>::iterator di = diag.begin()+start;
      VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin()+start;
      VEC_ITER_TYPE vector<T>::iterator ldk;
      VEC_ITER_TYPE vector<T>::iterator e = env.begin()+h;
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin()+h;
      VEC_ITER_TYPE vector<T>::iterator lej;
      VEC_ITER_TYPE vector<T>::iterator lel;
      vector<unsigned>::iterator xe;

//      for(i=0; i<bandwidth; i++)
      for(i=start; i<(unsigned)bandwidth; i++, ++di, ++ldi)
        {
//        ldiag[i]=diag[i];
        *ldi=*di;
//        for(k=0; k<i; k++, h++)
        for(k=0, xe = xenv.begin(), ldk = ldiag.begin(); k<i;
            k++, h++, ++e, ++le, ++xe, ++ldk)
          {
//          j=xenv[k];
          lej=lenv.begin()+*xe;
//          lenv[h]=env[h];
          *le=*e;
//          for(l=h-k; l<h; l++, j++)
          for(l=h-k, lel=lenv.begin()+h-k; l<h; l++, ++lel, ++lej)
            {
//            lenv[h]-=lenv[j]*lenv[l];
            *le-=*lej* *lel;
            }
//          lenv[h]/=ldiag[k];
          *le/=*ldk;
//          ldiag[i]-=lenv[h]*lenv[h];
          *ldi-=*le* *le;
          }
//        ldiag[i]=sqrt(ldiag[i]);
        *ldi=sqrt(*ldi);
        }

      if(start<=bandwidth)
        start=bandwidth;

//      for(i=bandwidth; i<dim; i++)
      for(i=start; i<dim; i++, ++di, ++ldi)
        {
//        ldiag[i]=diag[i];
        *ldi=*di;
//        for(k=0; k<bandwidth; k++, h++)
        for(k=0, xe = xenv.begin()+i-bandwidth+1, ldk=ldiag.begin()+i-bandwidth;
            k<(unsigned)bandwidth; k++, h++, ++e, ++le, ++xe, ++ldk)
          {
//          lenv[h]=env[h];
          *le=*e;
//          j=xenv[i-bandwidth+k+1]-k;
          lej=lenv.begin()+*xe-k;
//          for(l=h-k; l<h; l++, j++)
          for(l=h-k, lel=lenv.begin()+h-k; l<h; l++, ++lel, ++lej)
            {
//            lenv[h]-=lenv[l]*lenv[j];
            *le-=*lel* *lej;
            }
//          lenv[h]/=ldiag[i-bandwidth+k];
          *le/=*ldk;
//          ldiag[i]-=lenv[h]*lenv[h];
          *ldi-=*le* *le;
          }
//        ldiag[i]=sqrt(ldiag[i]);
        *ldi=sqrt(*ldi);
        }
      }
    else
      {
      unsigned i, ixenv, iband, ifirst, last, k, l, mstop, m, jstop, j;
//        unsigned i, ixenv, iband, ifirst, last, k, mstop, m, jstop, j;
      int kband;
      T temp = T(0);
      T s = T(0);
      T s1 = T(0);

      VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator di = diag.begin();
      vector<unsigned>::iterator xei = xenv.begin()+1;
      vector<unsigned>::iterator xek;
      VEC_ITER_TYPE vector<T>::iterator ek;
      VEC_ITER_TYPE vector<T>::iterator lem;
      VEC_ITER_TYPE vector<T>::iterator lel;
      VEC_ITER_TYPE vector<T>::iterator lek;
      VEC_ITER_TYPE vector<T>::iterator lej;
      VEC_ITER_TYPE vector<T>::iterator ldk;

      assert(*di>0);
      *ldi=sqrt(*di);
      ++di; ++ldi;
//        for(i=1; i<dim; i++)
      for(i=1; i<dim; i++, ++ldi, ++di)
        {
//          ixenv = xenv[i];
//          iband = xenv[i+1]-ixenv;
//          temp = diag[i];

        ixenv = *xei;
        ++xei;
        iband = *xei-ixenv;
        temp = *di;

        if(iband>0)
          {
          ifirst = i-iband;
          last = ixenv;
//            for(k=0; k<iband; k++)
          for(k=0, xek=xenv.begin()+ifirst, ek=env.begin()+ixenv,
              lek=lenv.begin()+ixenv, ldk=ldiag.begin()+ifirst;
              k<iband;
              k++, ++ek, ++lek, ++ldk)
            {
//              kband = xenv[k+ifirst+1]-xenv[k+ifirst];

            kband = -(int)*xek;
            ++xek;
            kband += *xek;

            if(kband > k)
              {
              kband = k;
              }

//              s=env[k+ixenv];
            s = *ek;

            l = k+ixenv-kband;
            lel = lenv.begin()+k+ixenv-kband;
 //           int h1 = *lel;

            if((kband>0) && (last>=l))
//              if((kband>0) && (last>=*lel))
              {
//                mstop=xenv[k+ifirst+1]-1;
//                for(m=xenv[k+ifirst+1]-kband; m<=mstop; m++)
              mstop = *xek - 1;
              for(lem=lenv.begin() + *xek - kband, m=*xek-kband; m<=mstop;
                                 ++lem, ++lel, m++)
                  {
//                  s=s-lenv[m]*lenv[l];
                l++;
                if(*lel!=0 && *lem!=0)
                  {
                  s -= *lel**lem;
                  }
                }
              }
            if(s!=0)
              {
//                lenv[k+ixenv]=s/ldiag[k+ifirst];
              *lek = s/ *ldk;
              last=k+ixenv;
              }
            else
              {
//                lenv[k+ixenv]=0;
              *lek=0;
              }
            }
//            jstop = xenv[i+1]-1;
          jstop = *xei-1;
//           for(j=ixenv; j<=jstop; j++)
         for(j=ixenv, lej=lenv.begin()+ixenv; j<=jstop; j++, ++lej)
            {
//              s1 = lenv[j];
            s1 = *lej;
            if(s1!=0)
              {
              temp = temp - s1*s1;
              }
            }
          }
        assert(temp>0);
//          ldiag[i]=sqrt(temp);
        *ldi = sqrt(temp);
        }
      }
    }
  decomposed=true;
  rational_decomposed=false;
  }


template<class T>
void envmatrix<T>::decomp_rational()
  {
  if(!rational_decomposed)
    {
    if(bandwidth==0)
      {
//      unsigned i;
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
//      for(i=0; i<dim; i++)
      for(; d!=diag.end();++ld, ++d)
        {
//        ldiag[i]=1/diag[i];
        *ld = 1/ *d;
        }
      }
    else if(bandwidth==1)
      {
      unsigned i;
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();

//      ldiag[0]=1/diag[0];
      *ld=1/ *d;
//      lenv[0]=env[0]*ldiag[0];
      *le=*e * *ld;
      ++ld; ++d;
      for(i=1; i<dim-1; i++)
        {
//        ldiag[i]=1/(diag[i]-env[i-1]*lenv[i-1]);
        *ld=1/(*d- *e* *le);
        ++e; ++le;
//        lenv[i]=env[i]*ldiag[i];
        *le=*e* *ld;
        ++d; ++ld;
        }
//      ldiag[dim-1]=1/(diag[dim-1]-env[dim-2]*lenv[dim-2]);
      *ld=1/(*d-*e* *le);
      }
    else if(bandwidth==2)
      {
//      unsigned i, h;
      unsigned i;
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();

//      ldiag[0]=1/diag[0];
      *ld=1/ *d;
//      lenv[0]=env[0]*ldiag[0];
      *le=*e* *ld;
      ++ld; ++d;
//      ldiag[1]=1/(diag[1]-env[0]*lenv[0]);
      *ld=1/(*d-*e* *le);
      ++e; ++le;
//      for(i=2, h=1; i<dim; i++)
      for(i=2; i<dim; i++)
        {
//      lenv[h]=env[h]*ldiag[i-2];
        *le=*e* *(ld-1);
        ++le; ++e;
//        h++;
//      lenv[h]=(env[h]-lenv[h-1]*lenv[h-2]/ldiag[i-2])*ldiag[i-1];
        *le=(*e-*(le-1)**(le-2)/ *(ld-1))**ld;
        ++ld; ++d;
//        ldiag[i]=1/(diag[i]-lenv[h]*lenv[h]/ldiag[i-1]
//               -lenv[h-1]*lenv[h-1]/ldiag[i-2]);
        *ld=1/(*d-*le**le/ *(ld-1)
                 -*(le-1)**(le-1)/ *(ld-2));
        ++le; ++e;
//        h++;
        }
      }
    else if(bandwidth>2)
      {
//      unsigned i, j, k, h, l, m;
      unsigned i, k, h, l;
      h=0;
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator ldk;
      VEC_ITER_TYPE vector<T>::iterator ldm;
      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      VEC_ITER_TYPE vector<T>::iterator lej;
      VEC_ITER_TYPE vector<T>::iterator lel;
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();
      vector<unsigned>::iterator xe;

      for(i=0; i<(unsigned)bandwidth; i++)
        {
//        ldiag[i]=diag[i];
        *ld=*d;
        xe=xenv.begin();
        ldk=ldiag.begin();
        for(k=0; k<i; k++, h++)
          {
//          j=xenv[k];
//          lenv[h]=env[h];
          *le=*e;
//          m=0;
          ldm=ldiag.begin();
          lej=lenv.begin()+*xe;
          ++xe;
          lel=lenv.begin()+h-k;
//          for(l=h-k; l<h; l++, j++,m++)
          for(l=h-k; l<h; l++)
            {
//            lenv[h] -= lenv[j]*lenv[l]/ldiag[m];
            *le -= *lej**lel/ *ldm;
            ldm++; ++lej; ++lel;
            }
//          lenv[h] *= ldiag[k];
          *le *= *ldk;
//          ldiag[i] -= lenv[h]*lenv[h]/ldiag[k];
          *ld -= *le* *le/ *ldk;
          ++le; ++e; ++ldk;
          }
//        ldiag[i] = 1/ldiag[i];
        *ld = 1/ *ld;
        ++ld; ++d;
        }

      for(i=bandwidth; i<dim; i++)
        {
//        ldiag[i]=diag[i];
        *ld=*d;
        xe=xenv.begin()+i-bandwidth+1;
        ldk=ldiag.begin()+i-bandwidth;
        for(k=0; k<(unsigned)bandwidth; k++, h++)
          {
//          lenv[h]=env[h];
          *le=*e;
//          j=xenv[i-bandwidth+k+1]-k;
//          m=i-bandwidth;
          ldm=ldiag.begin()+i-bandwidth;
          lel=lenv.begin()+h-k;
          lej=lenv.begin()+*xe-k;
          ++xe;
//          for(l=h-k; l<h; l++, j++, m++)
          for(l=h-k; l<h; l++)
            {
//            lenv[h] -= lenv[j]*lenv[l]/ldiag[m];
            *le -= *lej* *lel/ *ldm;
            ++ldm; ++lel; ++lej;
            }
//          lenv[h]*=ldiag[i-bandwidth+k];
          *le*=*ldk;
//          ldiag[i] -= lenv[h]*lenv[h]/ldiag[m];
          *ld -= *le**le/ *ldm;
          ++le; ++e; ++ldk;
          }
//        ldiag[i] = 1/ldiag[i];
        *ld = 1/ *ld;
        ++ld; ++d;
        }
      }
    else
      {
      int i, k, j, h, jstart;
      h=0;

      VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
      VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
      VEC_ITER_TYPE vector<T>::iterator ldk;
      VEC_ITER_TYPE vector<T>::iterator ldj;
      VEC_ITER_TYPE vector<T>::iterator e = env.begin();
      VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
      vector<unsigned>::iterator xe = xenv.begin();
      vector<unsigned>::iterator xek;

      for(i=0; i<(int)dim; i++)
        {
//        ldiag[i]=diag[i];
        *ld=*d;
//        for(k=i-xenv[i+1]+xenv[i]; k<i; k++, h++)
        xek=xenv.begin()+i-*(xe+1)+*xe;
        ldk=ldiag.begin()+i-*(xe+1)+*xe;
        for(k=i-*(xe+1)+*xe; k<i; k++, h++)
          {
//          lenv[h] = env[h];
          *le = *e;
//          if((xenv[i+1]-xenv[i])>(xenv[k+1]-xenv[k]))
          if((*(xe+1)-*xe)>(*(xek+1)-*xek))
            {
            jstart = i-*(xe+1)+*xe;
            }
          else
            {
            jstart = i-*(xek+1)+*xek;
            }
          ldj=ldiag.begin()+jstart;
          for(j=jstart; j<k; j++)
            {
//            lenv[h] -= getL(k,j)*getL(i,j)/getL(j,j);
            *le -= getL(k,j)*getL(i,j)/ *ldj;
            ++ldj;
            }
//          lenv[h] *= getL(k,k);
          *le *= *ldk;
//          ldiag[i] -= lenv[h]*lenv[h]/ldiag[k];
          *ld -= *le**le/ *ldk;
          ++le; ++e; ++xek; ++ldk;
          }
//        ldiag[i] = 1/ldiag[i];
        *ld = 1/ *ld;
        ++d; ++ld; ++xe;
        }
      }
    }
  rational_decomposed=true;
  decomposed=false;
  }


template<class T>
void envmatrix<T>::solve(vector<T> &b)
  {
  solveL(b);
  solveU(b);
  }

template<class T>
void envmatrix<T>::solve(datamatrix &b)
  {
  solveL(b);
  solveU(b);
  }


template<class T>
void envmatrix<T>::solve(const datamatrix &b, datamatrix &res)
  {
  assert(b.rows()==res.rows());
  assert(b.cols()==res.cols()==1);
  solveL(b, res);
  solveU(res);
  }

template<class T>
void envmatrix<T>::solve(const datamatrix &b,const datamatrix &bhelp,
                         datamatrix &res)
  {
  assert(b.rows()==res.rows());
  assert(b.cols()==res.cols()==1);
  solveL(b, res);
  solveU(res, bhelp);
  }

template<class T>
void envmatrix<T>::solveL(vector<T> & b)
  {
  if(!decomposed)
    {
    decomp();
    }
  if(bandwidth!=0)
    {
    unsigned ifirst = 0;
    while((ifirst<dim) && (b[ifirst]==0))
      {
      ifirst++;
      }
    unsigned last = 0;
    unsigned i, l, kstop, k;
    int iband;
    T s = T(0);

    vector<unsigned>::iterator xei = xenv.begin()+ifirst;

    for(i=ifirst; i<dim; i++)
      {
      iband = xenv[i+1]-xenv[i];
      s=b[i];
      l = i-iband;
      if((iband>0) && (last>=l))
        {
        kstop=xenv[i+1]-1;
        k;
         for(k=xenv[i+1]-iband; k<=kstop; k++)
          {
          s=s-lenv[k]*b[l];
          l++;
          }
        }
      if(s!=0)
        {
        b[i]=s/ldiag[i];
        last=i;
        }
      }
    }
  else
    {
    unsigned i;
    for(i=0; i<dim; i++)
      {
      b[i] /= ldiag[i];
      }
    }
  }

template<class T>
void envmatrix<T>::solveL(datamatrix & b)
  {
  if(!decomposed)
    {
    decomp();
    }
  if(bandwidth==0)
    {
//    unsigned i;
    VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin();
    T* bi = b.getV();

//    for(i=0; i<dim; i++)
    for(; ldi!=ldiag.end(); bi++, ++ldi)
      {
//      b(i,0) /= ldiag[i];
      *bi /= *ldi;
      }
    }
  else if(bandwidth==1)
    {
    unsigned i;
    VEC_ITER_TYPE vector<T>::iterator le=lenv.begin();
    VEC_ITER_TYPE vector<T>::iterator ld=ldiag.begin();
    T* bi = b.getV();

//    b(0,0) = b(0,0)/ldiag[0];
    *bi /= *ld;
    bi++; ++ld;
//    for(i=1; i<dim; i++)
    for(i=1; i<dim; i++, bi++, ++ld, ++le)
      {
//      b(i,0)=(b(i,0)-b(i-1,0)*lenv[i-1])/ldiag[i];
      *bi=(*bi-*(bi-1)* *le)/ *ld;
      }
    }
  else if(bandwidth==2)
    {
//    unsigned i, k;
    unsigned i;
    VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
    VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
    T* bi = b.getV();

//    b(0,0) = b(0,0)/ldiag[0];
    *bi /= *ld;
    bi++; ++ld;
//    b(1,0) = (b(1,0)-b(0,0)*lenv[0])/ldiag[1];
    *bi = (*bi-*(bi-1)* *le)/ *ld;
    bi++; ld++; le++;
//    for(i=2, k=1; i<dim; i++, k+=2)
    for(i=2; i<dim; i++, bi++, ++ld, le+=2)
      {
//      b(i,0) = (b(i,0)-b(i-2,0)*lenv[k]-b(i-1,0)*lenv[k+1])/ldiag[i];
      *bi = (*bi-*(bi-2)* *le-*(bi-1)* *(le+1))/ *ld;
      }
    }
  else if(bandwidth>2)
    {
//    unsigned i, k, l;
    unsigned i, k;
    VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
    VEC_ITER_TYPE vector<T>::iterator le;
    vector<unsigned>::iterator xe=xenv.begin();
    T* bi = b.getV();
    T* bk;

//    for(i=0; i<bandwidth; i++)
    for(i=0; i<(unsigned)bandwidth; i++, ++xe, bi++, ++ld)
      {
//      for(k=0, l=xenv[i]; k<i; k++, l++)
      for(k=0, le=lenv.begin()+*xe, bk = b.getV(); k<i; k++, bk++, ++le)
        {
//          b(i,0)-=b(k,0)*lenv[l];
        *bi -= *bk* *le;
        }
//      b(i,0)/=ldiag[i];
      *bi /= *ld;
      }
//    for(i=bandwidth; i<dim; i++)
    for(i=bandwidth; i<dim; i++, ++xe, bi++, ++ld)
      {
//      for(k=i-bandwidth, l=xenv[i]; k<i; k++, l++)
      for(k=i-bandwidth, le=lenv.begin()+*xe, bk=b.getV()+i-bandwidth; k<i;
          k++, bk++, ++le)
        {
//          b(i,0)-=b(k,0)*lenv[l];
        *bi -= *bk* *le;
        }
//      b(i,0)/=ldiag[i];
      *bi /= *ld;
      }
    }
  else
    {
    unsigned ifirst = 0;
    T* bi = b.getV();
//    while((ifirst<dim) && (b(ifirst,0)==0))
    while((ifirst<dim) && (*bi==0))
      {
//      ifirst++;
      ifirst++, bi++;
      }
    unsigned last = 0;
    unsigned i, l, kstop, k;
    int iband;
    T s = T(0);
    vector<unsigned>::iterator xei = xenv.begin()+ifirst;
    VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin()+ifirst;
    VEC_ITER_TYPE vector<T>::iterator lek;
    T* bl;

//    for(i=ifirst; i<dim; i++)
      for(i=ifirst; i<dim; i++, bi++, ++ldi)
      {
//      iband = xenv[i+1]-xenv[i];
      iband = -(int)*xei;
      ++xei;
      iband += *xei;

//      s=b(i,0);
      s=*bi;

      l = i-iband;
      bl = b.getV()+l;
//      if((iband>0) && (last>=l))
      if((iband>0) && (last>=l))
        {
//        kstop=xenv[i+1]-1;
//        for(k=xenv[i+1]-iband; k<=kstop; k++)
        kstop = *xei - 1;
        for(k=*xei-iband, lek=lenv.begin()+*xei-iband; k<=kstop;
                          k++, ++lek, bl++, l++)
          {
//          s=s-lenv[k]*b(l,0);
//          l++;
          s-=*lek**bl;
          }
        }
      if(s!=0)
        {
//        b(i,0)=s/ldiag[i];
        *bi= s/ *ldi;
        last=i;
        }
      }
    }
  }

template<class T>
void envmatrix<T>::solveL(const datamatrix & b, datamatrix & res)
  {
  if(!decomposed)
    {
    decomp();
    }
  if(bandwidth==0)
    {
//    unsigned i;
    T* ri = res.getV();
    T* bi = b.getV();
    VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin();
//    for(i=0; i<dim; i++)
    for(; ldi!=ldiag.end(); ri++, ++ldi, bi++)
      {
//      res(i,0) /= ldiag[i];
      *ri = *bi/ *ldi;
      }
    }
  else if(bandwidth==1)
    {
    unsigned i;
    VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
    VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
    T* bi = b.getV();
    T* resi = res.getV();

//    res(0,0) = b(0,0)/ldiag[0];
    *resi = *bi/ *ld;
    resi++; bi++; ++ld;
//    for(i=1; i<dim; i++)
    for(i=1; i<dim; i++, resi++, bi++, ++ld, ++le)
      {
//      res(i,0)=(b(i,0)-res(i-1,0)*lenv[i-1])/ldiag[i];
      *resi=(*bi-*(resi-1)* *le)/ *ld;
      }
    }
  else if(bandwidth==2)
    {
//    unsigned i, k;
    unsigned i;
    VEC_ITER_TYPE  vector<T>::iterator ld = ldiag.begin();
    VEC_ITER_TYPE vector<T>::iterator le = lenv.begin();
    T* bi = b.getV();
    T* resi = res.getV();

//    res(0,0) = b(0,0)/ldiag[0];
    *resi = *bi/ *ld;
    resi++; bi++; ++ld;
//    res(1,0) = (b(1,0)-res(0,0)*lenv[0])/ldiag[1];
    *resi = (*bi-*(resi-1)* *le)/ *ld;
    resi++; bi++; ++le; ++ld;
//    for(i=2, k=1; i<dim; i++, k+=2)
    for(i=2; i<dim; i++, resi++, bi++, ++ld, le+=2)
      {
//      res(i,0) = (b(i,0)-res(i-2,0)*lenv[k]-res(i-1,0)*lenv[k+1])/ldiag[i];
      *resi = (*bi-*(resi-2)* *le-*(resi-1)* *(le+1))/ *ld;
      }
    }
  else if(bandwidth>2)
    {
    unsigned i, k;
    VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
    VEC_ITER_TYPE vector<T>::iterator le;
    vector<unsigned>::iterator xe = xenv.begin();
    T* bi = b.getV();
    T* resi = res.getV();
    T* resk;

//    for(i=0; i<bandwidth; i++)
    for(i=0; i<(unsigned)bandwidth; i++, resi++, bi++, ++ld, ++xe)
      {
//      res(i,0)=b(i,0);
      *resi=*bi;
//      for(k=0, l=xenv[i]; k<i; k++, l++)
      for(k=0, resk = res.getV(), le=lenv.begin()+*xe; k<i; k++, resk++, ++le)
        {
//        res(i,0) -= res(k,0)*lenv[l];
        *resi -= *resk* *le;
        }
//      res(i,0) /= ldiag[i];
      *resi /= *ld;
      }
//    for(i=bandwidth; i<dim; i++)
    for(i=bandwidth; i<dim; i++, ++ld, bi++, resi++, ++xe)
      {
//      res(i,0)=b(i,0);
      *resi=*bi;
//      for(k=i-bandwidth, l=xenv[i]; k<i; k++, l++)
      for(k=i-bandwidth, resk=res.getV()+i-bandwidth , le=lenv.begin()+*xe;
          k<i; k++, ++le, resk++)
        {
//        res(i,0) -= res(k,0)*lenv[l];
        *resi -= *resk* *le;
        }
//      res(i,0) /= ldiag[i];
      *resi /= *ld;
      }
    }
  else
    {
    unsigned ifirst = 0;
    T* bi = b.getV();
    T* ri = res.getV();

//    while((ifirst<dim) && (b(ifirst,0)==0))
    while((ifirst<dim) && (*bi==0))
      {
//      res(ifirst,0)=0;
      *ri=0;
      ri++;
      bi++;
      ifirst++;
      }
    unsigned last = 0;
    unsigned i,l, kstop, k;
    int iband;
    T s = T(0);
    vector<unsigned>::iterator xei = xenv.begin()+ifirst;
    VEC_ITER_TYPE vector<T>::iterator lek;
    VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin()+ifirst;
    T* rl;

    for(i=ifirst; i<dim; i++, bi++, ri++, ++ldi)
      {
//      iband = xenv[i+1]-xenv[i];
      iband = -(int)*xei;
      ++xei;
      iband += *xei;

//    s=b(i,0);
      s=*bi;
      l = i-iband;
      if((iband>0) && (last>=l))
        {
//       kstop=xenv[i+1]-1;
//        for(k=xenv[i+1]-iband; k<=kstop; k++)
        kstop=*xei-1;
        for(k=*xei-iband, lek = lenv.begin()+*xei-iband, rl = res.getV()+l;
                          k<=kstop; k++, ++lek, rl++, l++)
          {
//          s=s-lenv[k]*res(l,0);
//          l++;
          s-= *lek**rl;
          }
        }
      if(s!=0)
        {
//        res(i,0)=s/ldiag[i];
        *ri = s/ *ldi;
        last=i;
        }
      else
        {
//        res(i,0)=0;
        *ri=0;
        }
      }
    }
  }

template<class T>
void envmatrix<T>::solveU(vector<T> &b)
  {
  if(!decomposed)
    {
    decomp();
    }
  if(bandwidth!=0)
    {
    int i;
    unsigned iband, l, k;
    T s = T(0);
    for(i=dim-1; i>=0; i--)
      {
      if(b[i]!=0)
        {
        s = b[i]/ldiag[i];
        b[i] = s;
        iband = xenv[i+1]-xenv[i];
        if(iband>0)
          {
          l = xenv[i+1]-iband;
          k;
          for(k=i-iband; k<i; k++)
            {
            b[k] = b[k] - s*lenv[l];
            l++;
            }
          }
        }
      }
    }
  else
    {
    unsigned i;
    for(i=0; i<dim; i++)
      {
      b[i] /= ldiag[i];
      }
    }
  }

template<class T>
void envmatrix<T>::solveU(datamatrix &b)
  {
  if(!decomposed)
    {
    decomp();
    }
  if(bandwidth==0)
    {
//    unsigned i;
    T* bi = b.getV();
    VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin();

//    for(i=0; i<dim; i++)
    for(; ldi!=ldiag.end();bi++, ++ldi)
      {
//      b(i,0) /= ldiag[i];
      *bi /= *ldi;
      }
    }
  else if(bandwidth==1)
    {
    VEC_ITER_TYPE vector<T>::iterator ld = ldiag.end()-1;
    VEC_ITER_TYPE vector<T>::iterator le = lenv.end()-1;
    T* bi = b.getV()+dim-1;
    int i;

//    b(dim-1,0)=b(dim-1,0)/ldiag[dim-1];
    *bi /= *ld;
    bi--; --ld;
//    for(i=dim-2; i>=0; i--)
    for(i=dim-2; i>=0; i--, bi--, --le, --ld)
      {
//      b(i,0)=(b(i,0)-b(i+1,0)*lenv[i])/ldiag[i];
      *bi=(*bi-*(bi+1)* *le)/ *ld;
      }
    }
  else if(bandwidth==2)
    {
    VEC_ITER_TYPE vector<T>::iterator ld = ldiag.end()-1;
    VEC_ITER_TYPE vector<T>::iterator le = lenv.end()-1;
    vector<unsigned>::iterator xe = xenv.end()-2;
    T* bi = b.getV()+dim-1;
//    int i, k;
    int i;

//    b(dim-1,0)=b(dim-1,0)/ldiag[dim-1];
    *bi/=*ld;
    bi--; --ld;
//    b(dim-2,0)=(b(dim-2,0)-b(dim-1,0)*lenv[xenv[dim]-1])/ldiag[dim-2];
    *bi=(*bi-*(bi+1)* *le)/ *ld;
    bi--; --ld;
//    for(i=dim-3, k=xenv[dim-1]; i>=0; i--, k-=2)
    for(i=dim-3, le=lenv.begin()+*xe; i>=0; i--, --xe, --ld, bi--, le-=2)
      {
//      b(i,0)=(b(i,0)-b(i+1,0)*lenv[k-1]-b(i+2,0)*lenv[k])/ldiag[i];
      *bi=(*bi-*(bi+1)* *(le-1)-*(bi+2)* *le)/ *ld;
      }
    }
  else if(bandwidth>2)
    {
    VEC_ITER_TYPE vector<T>::iterator ld = ldiag.end()-1;
    vector<unsigned>::iterator xe = xenv.end()-1;
    VEC_ITER_TYPE vector<T>::iterator le;
    T* bi = b.getV()+dim-1;
    T* bk;

//    int i, k, l, iband;
    int i, k, iband;
    T s = T(0);
//    for(i=dim-1; i>=bandwidth; i--)
    for(i=dim-1; i>=bandwidth; i--, bi--, --ld, --xe)
      {
//      s=b(i,0)/ldiag[i];
      s=*bi / *ld;
//      b(i,0)=s;
      *bi=s;
//      for(k=i-bandwidth, l=xenv[i+1]-bandwidth; k<i; k++, l++)
      for(k=i-bandwidth, bk = b.getV()+i-bandwidth,
          le=lenv.begin()+*xe-bandwidth; k<i; k++, bk++, ++le)
        {
//        b(k,0)=b(k,0)-s*lenv[l];
        *bk-=s* *le;
        }
      }
//    for(i=bandwidth-1; i>=0; i--)
    for(i=bandwidth-1; i>=0; i--, bi--, --ld, --xe)
      {
//      s=b(i,0)/ldiag[i];
      s=*bi/ *ld;
//      b(i,0)=s;
      *bi=s;
//      iband = xenv[i+1]-xenv[i];
      iband = *xe-*(xe-1);
//      for(k=i-iband, l=xenv[i+1]-iband; k<i; k++, l++)
      for(k=i-iband, bk=b.getV()+i-iband, le=lenv.begin()+*xe-iband; k<i;
          k++, ++le, ++bk)
        {
//        b(k,0)=b(k,0)-s*lenv[l];
        *bk-=s* *le;
        }
      }
    }
  else
    {
    int i;
    unsigned iband, l, k, h;
    T s = T(0);

    T* bi = b.getV() + dim -1;
    T* bk;
//    vector<unsigned>::iterator xei = xenv.begin()+dim;
    vector<unsigned>::iterator xei = xenv.end()-1;
//    vector<T>::iterator ldi = ldiag.begin()+dim-1;
    VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.end()-1;
    VEC_ITER_TYPE vector<T>::iterator lel;

//    for(i=dim-1; i>=0; i--)
    for(i=dim-1; i>=0; i--, bi--, --ldi)
      {
//      if(b(i,0)!=0)
      if(*bi!=0)
        {
//        s = b(i,0)/ldiag[i];
        s = *bi/ *ldi;
//        b(i,0) = s;
        *bi = s;
//      iband = xenv[i+1]-xenv[i];
        h=*xei;
        --xei;
        iband=h-*xei;
        if(iband>0)
          {
//          l = xenv[i+1]-iband;
          l = h-iband;
//          for(k=i-iband; k<i; k++, l++)
          for(k=i-iband, bk = b.getV()+i-iband, lel = lenv.begin()+h-iband;
                         k<i; k++, l++, bk++, ++lel)
            {
//            b(k,0) = b(k,0) - s*lenv[l];
            *bk -= s**lel;
            }
          }
        }
      }
    }
  }

template<class T>
void envmatrix<T>::solveU(datamatrix &b, const datamatrix &bhelp)
  {
  if(!decomposed)
    {
    decomp();
    }
  if(bandwidth==0)
    {
//    unsigned i;
    T* bi = b.getV();
    T* bhelpi = bhelp.getV();
    VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.begin();

//    for(i=0; i<dim; i++)
    for(; ldi!=ldiag.end();bi++, ++ldi, bhelpi++)
      {
//      b(i,0) /= ldiag[i];
//      *bi =  *bi/ *ldi + bhelp(i,0);
      *bi =  *bi/ *ldi + *bhelpi;
      }
    }
  else if(bandwidth==1)
    {
    VEC_ITER_TYPE vector<T>::iterator ld = ldiag.end()-1;
    VEC_ITER_TYPE vector<T>::iterator le = lenv.end()-1;
    T* bi = b.getV()+dim-1;
    T* bhelpi = bhelp.getV()+dim-1;

    int i;
    T s = T(0);
//    for(i=dim-1; i>0; i--)
    for(i=dim-1; i>0; i--, bhelpi--, --ld, --le)
      {
//      s=b(i,0)/ldiag[i];
      s=*bi/ *ld;
//      b(i,0)=s+bhelp(i,0);
      *bi=s+*bhelpi;
      bi--;
//      b(i-1,0)-=s*lenv[i-1];
      *bi-=s* *le;
      }
//    b(0,0)=b(0,0)/ldiag[0]+bhelp(0,0);
    *bi=*bi/ *ld+*bhelpi;
    }
  else if(bandwidth==2)
    {
    VEC_ITER_TYPE vector<T>::iterator ld = ldiag.end()-1;
    VEC_ITER_TYPE vector<T>::iterator le = lenv.end()-1;
    T* bi = b.getV()+dim-1;
    T* bhelpi = bhelp.getV()+dim-1;

    int i;
    T s = T(0);
//    for(i=dim-1, k=xenv[dim]; i>1; i--, k-=2)
    for(i=dim-1; i>1; i--, --le, --ld, bhelpi--)
      {
//      s=b(i,0)/ldiag[i];
      s=*bi/ *ld;
//      b(i,0)=s+bhelp(i,0);
      *bi=s+*bhelpi;
      bi--;
//      b(i-1,0)-=s*lenv[k-1];
      *bi-=s* *le;
      --le;
//      b(i-2,0)-=s*lenv[k-2];
      *(bi-1)-=s* *le;
      }
//    s=b(1,0)/ldiag[1];
    s=*bi/ *ld;
//    b(1,0)=s+bhelp(1,0);
    *bi=s+*bhelpi;
    bi--; bhelpi--; --ld;
//    b(0,0)-=s*lenv[0];
    *bi -= s* *le;
//    b(0,0)=b(0,0)/ldiag[0]+bhelp(0,0);
    *bi=*bi/ *ld+*bhelpi;
    }
  else if(bandwidth>2)
    {
    VEC_ITER_TYPE vector<T>::iterator ld = ldiag.end()-1;
    VEC_ITER_TYPE vector<T>::iterator le;
    vector<unsigned>::iterator xe = xenv.end()-1;
    T* bi = b.getV()+dim-1;
    T* bhelpi = bhelp.getV()+dim-1;
    T* bk;

    int i, k, iband;
    T s = T(0);
//    for(i=dim-1; i>=bandwidth; i--)
    for(i=dim-1; i>=bandwidth; i--, bi--, --ld, --xe, bhelpi--)
      {
//      s=b(i,0)/ldiag[i];
      s=*bi/ *ld;
//      b(i,0)=s+bhelp(i,0);
      *bi=s+*bhelpi;
//      for(k=i-bandwidth, l=xenv[i+1]-bandwidth; k<i; k++, l++)
      for(k=i-bandwidth, le=lenv.begin()+*xe-bandwidth, bk=b.getV()+i-bandwidth;
          k<i; k++, ++le, bk++)
        {
//        b(k,0)=b(k,0)-s*lenv[l];
        *bk-=s* *le;
        }
      }
//    for(i=bandwidth-1; i>=0; i--)
    for(i=bandwidth-1; i>=0; i--, bi--, bhelpi--, --ld, --xe)
      {
//      s=b(i,0)/ldiag[i];
      s=*bi/ *ld;
//      b(i,0)=s+bhelp(i,0);
      *bi=s+*bhelpi;
//      iband = xenv[i+1]-xenv[i];
      iband = *xe-*(xe-1);
//      for(k=i-iband, l=xenv[i+1]-iband; k<i; k++, l++)
      for(k=i-iband, le=lenv.begin()+*xe-iband, bk=b.getV()+i-iband;
          k<i; k++, ++le, bk++)
        {
//        b(k,0)=b(k,0)-s*lenv[l];
        *bk-=s* *le;
        }
      }
    }
  else
    {
    int i;
    unsigned iband, l, k, h;
    T s = T(0);

    T* bi = b.getV() + dim -1;
    T* bhelpi = bhelp.getV() + dim -1;
    T* bk;
//    vector<unsigned>::iterator xei = xenv.begin()+dim;
    vector<unsigned>::iterator xei = xenv.end()-1;
//    vector<T>::iterator ldi = ldiag.begin()+dim-1;
    VEC_ITER_TYPE vector<T>::iterator ldi = ldiag.end()-1;
    VEC_ITER_TYPE vector<T>::iterator lel;

//    for(i=dim-1; i>=0; i--)
    for(i=dim-1; i>=0; i--, bi--, bhelpi--, --ldi)
      {
//      if(b(i,0)!=0)
      if(*bi!=0)
        {
//        s = b(i,0)/ldiag[i];
        s = *bi/ *ldi;
//        b(i,0) = s;
//        *bi = s + bhelp(i,0);
        *bi = s + *bhelpi;
//      iband = xenv[i+1]-xenv[i];
        h=*xei;
        --xei;
        iband=h-*xei;
        if(iband>0)
          {
//          l = xenv[i+1]-iband;
          l = h-iband;
//          for(k=i-iband; k<i; k++, l++)
          for(k=i-iband, bk = b.getV()+i-iband, lel = lenv.begin()+h-iband;
                         k<i; k++, l++, bk++, ++lel)
            {
//            b(k,0) = b(k,0) - s*lenv[l];
            *bk -= s**lel;
            }
          }
        }
      }
    }
  }

template <class T>
void envmatrix<T>::inverse_envelope(envmatrix<T> & inv)
  {
  assert(xenv==inv.getXenv());
  assert(bandwidth==inv.getBandwidth());
  assert(dim==inv.getDim());
  if(bandwidth==0)
    {
    unsigned i;
    VEC_ITER_TYPE vector<T>::iterator invdiag=inv.getDiagIterator();
    VEC_ITER_TYPE vector<T>::iterator d=diag.begin();
    for(i=0; i<dim; i++, invdiag++, ++d)
      {
      *invdiag = 1/ *d;
      }
    }
  else if(bandwidth==1)
    {
    int i;
    VEC_ITER_TYPE vector<T>::iterator invdiag=inv.getDiagIterator()+dim-1;
    VEC_ITER_TYPE vector<T>::iterator ld=ldiag.end()-1;
    VEC_ITER_TYPE vector<T>::iterator invenv=inv.getEnvIterator()+env.size()-1;
    VEC_ITER_TYPE vector<T>::iterator le=lenv.end()-1;

    decomp_rational();

//    *invdiag = ldiag[dim-1];
    *invdiag = *ld;

    for(i=dim-2; i>=0; i--)
      {
//      *invenv = -lenv[i]* *invdiag;
      *invenv = -*le* *invdiag;
      invdiag--; ld--;
//      *invdiag = ldiag[i]-lenv[i]* *invenv;
      *invdiag = *ld-*le* *invenv;
      invenv--; le--;
      }
    }
  else if(bandwidth==2)
    {
    int i;
    VEC_ITER_TYPE vector<T>::iterator invdiag=inv.getDiagIterator()+dim-1;
    VEC_ITER_TYPE vector<T>::iterator ld=ldiag.end()-1;
    VEC_ITER_TYPE vector<T>::iterator invenv=inv.getEnvIterator()+env.size()-1;
    VEC_ITER_TYPE vector<T>::iterator lenvit=lenv.end()-1;

    decomp_rational();

//    *invdiag = ldiag[dim-1];
    *invdiag = *ld;
    *invenv = -*invdiag* *lenvit;
    invenv--; lenvit--;
//    *invenv = -inv(dim-1,dim-2)*getL(dim-2,dim-3)
//              -inv(dim-1,dim-1)*getL(dim-1,dim-3);
    *invenv = -*(invenv+1)**(lenvit-1)
              -*invdiag**lenvit;

    invdiag--; ld--;
//    *invdiag = ldiag[dim-2]-inv(dim-1,dim-2)*getL(dim-1,dim-2);
    *invdiag = *ld-*(invenv+1)**(lenvit+1);
    invenv--; lenvit--;

    for(i=dim-2; i>1; i--)
      {
//      *invenv = -inv(i,i)*getL(i,i-1)-inv(i,i+1)*getL(i+1,i-1);
      *invenv = -*invdiag**lenvit-*(invenv+2)**(lenvit+1);
      invenv--; lenvit--;
//      *invenv = -inv(i,i-1)*getL(i-1,i-2)-inv(i,i)*getL(i,i-2);
      *invenv = -*(invenv+1)**(lenvit-1)-*invdiag**lenvit;

      invdiag--; ld--;
//      *invdiag = ldiag[i-1]-inv(i,i-1)*getL(i,i-1)-inv(i+1,i-1)*getL(i+1,i-1);
      *invdiag = *ld-*(invenv+1)**(lenvit+1)-*(invenv+2)**(lenvit+2);
      invenv--; lenvit--;
      }

//    *invenv = -inv(1,1)*getL(1,0)-inv(1,2)*getL(2,0);
    *invenv = -*invdiag**lenvit-*(invenv+2)**(lenvit+1);
    invdiag--; ld--;
//    *invdiag = ldiag[0]-inv(1,0)*getL(1,0)-inv(2,0)*getL(2,0);
    *invdiag = *ld-*invenv**lenvit-*(invenv+1)**(lenvit+1);

    }
  else if(bandwidth>2)
    {
    int i, k, l, upperk, upperl;
    VEC_ITER_TYPE vector<T>::iterator invdiag=inv.getDiagIterator()+dim-1;
    VEC_ITER_TYPE vector<T>::iterator ld=ldiag.end()-1;
    VEC_ITER_TYPE vector<T>::iterator invenv=inv.getEnvIterator()+env.size()-1;

    decomp_rational();

    for(i=dim-1; i>0; i--)
      {
      if(dim-i>(unsigned)bandwidth)
        {
        upperk=bandwidth;
        }
      else
        {
        upperk=dim-i-1;
        }
//      *invdiag = ldiag[i];
      *invdiag = *ld;
      for(k=0; k<upperk; k++)
        {
        *invdiag -= inv(i+k+1,i)*getL(i+k+1,i);
        }
      invdiag--; ld--;

      if(i>bandwidth)
        {
        upperk=bandwidth;
        }
      else
        {
        upperk=i;
        }
      for(k=0; k<upperk; k++)
        {
//        if(i==(int)dim-1)
        if(dim-i+k<bandwidth)
          {
//          upperl = k+1;
          upperl = dim-i+k;
          }
        else
          {
          upperl = bandwidth;
          }
        *invenv=0;
        for(l=0; l<upperl; l++)
          {
          *invenv -= inv(i,i-k+l)*getL(i-k+l,i-k-1);
          }
        invenv--;
        }
      }
//    *invdiag = ldiag[0];
    *invdiag = *ld;
    for(k=0; k<(int)bandwidth; k++)
      {
      *invdiag -= inv(k+1,0)*getL(k+1,0);
      }
    }
  else
    {
    int i, k, l, upperk, upperl;
    VEC_ITER_TYPE vector<T>::iterator invdiag=inv.getDiagIterator()+dim-1;
    VEC_ITER_TYPE vector<T>::iterator invenv=inv.getEnvIterator()+env.size()-1;

    VEC_ITER_TYPE vector<T>::iterator invenv1;
    for(invenv1=inv.getEnvIterator(); invenv1<invenv; invenv1++)
      {
      *invenv1=0;
      }

    int maxbw=0;
    for(i=0; i<(int)dim; i++)
      {
      if((int)xenv[i+1]-(int)xenv[i]>maxbw)
        {
        maxbw=xenv[i+1]-xenv[i];
        }
      }

    decomp_rational();

    for(i=dim-1; i>=0; i--)
      {
      *invdiag = ldiag[i];

      if(i+maxbw>=(int)dim)
        {
        upperk=dim;
        }
      else
        {
        upperk=i+maxbw+1;
        }

      for(k=i+1; k<upperk; k++)
        {
        *invdiag -= inv(k,i)*getL(k,i);
        }
      invdiag--;

      for(k=i-1; k>=(i-(int)xenv[i+1]+(int)xenv[i]); k--)
        {

        if(k+maxbw>=(int)dim)
          {
          upperl=dim;
          }
        else
          {
          upperl=k+maxbw+1;
          }

        for(l=k+1; l<upperl; l++)
          {
          *invenv -= inv(i,l)*getL(l,k);
          }
        invenv--;
        }
      }
    }
  inv.setDecomposed(false);
  inv.setRational_decomposed(false);
  }

//------------------------------------------------------------------------------
//------ Functions to assess elements or characteristics of the matrix----------
//------------------------------------------------------------------------------

template<class T>
T envmatrix<T>::getL(const unsigned & i, const unsigned & j) const
  {
  int kl, ku, zeroes;
  unsigned ih, jh;

  assert(i<dim);
  assert(j<dim);

  if(i>j)
    {
    ih = i;
    jh = j;
    }
  else if(i<j)
    {
    ih = j;
    jh = i;
    }
  else
    return ldiag[i];

  kl=xenv[ih];
  ku=xenv[ih+1];
  zeroes=ih-ku+kl;

  if(jh<zeroes)
    return T(0);
  else
    return lenv[kl+jh-zeroes];
  }

template<class T>
statmatrix<T> envmatrix<T>::getL(void) const
  {
  statmatrix<T> R(dim,dim,0);
  unsigned i,j;
  for (i=0;i<dim;i++)
    for (j=0;j<=i;j++)
      R(i,j) = getL(i,j);

  return R;
  }


template<class T>
statmatrix<T> envmatrix<T>::get(void) const
  {
  statmatrix<T> S(dim,dim,0);
  unsigned i,j;
  for (i=0;i<dim;i++)
    for (j=0;j<dim;j++)
      S(i,j) = get(i,j);

  return S;
  }


template<class T>
int envmatrix<T>::getBandwidth(void) const
  {
  return bandwidth;
  }

template<class T>
unsigned envmatrix<T>::getDim(void) const
  {
  return dim;
  }

template<class T>
unsigned envmatrix<T>::getXenv(const unsigned & i) const
  {
  return xenv[i];
  }

template<class T>
T & envmatrix<T>::getEnv(const unsigned & i)
  {
  return env[i];
  }

template<class T>
vector<T> envmatrix<T>::getDiag()
  {
  return diag;
  }

template<class T>
T envmatrix<T>::getDiag(unsigned i)
  {
  return diag[i];
  }

template<class T>
vector<T> envmatrix<T>::getEnv()
  {
  return env;
  }

template<class T>
vector<unsigned> envmatrix<T>::getXenv() const
  {
  return xenv;
  }

template<class T>
T envmatrix<T>::getLogDet()
  {
  if(!decomposed)
    {
    decomp();
    }
  T logdet=0;
  VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
  for(; ld!=ldiag.end(); ++ld)
    {
    logdet+=log(*ld);
    }
  return 2*logdet;
  }


template<class T>
T envmatrix<T>::getLogDet_save(bool error)
  {
  error=false;
  if(!decomposed)
    {
    error = decomp_save();
    }

  T logdet=0;
  if (error == false)
    {

    VEC_ITER_TYPE vector<T>::iterator ld = ldiag.begin();
    for(; ld!=ldiag.end(); ++ld)
      {
      if ((*ld) > logmin && (*ld) < logmax)
        logdet+=log(*ld);
      else
        {
        error = true;
        logdet += logmin;
        }
      }
    }

  return 2*logdet;
  }


template<class T>
T envmatrix<T>::traceOfProduct(envmatrix<T> & B)
  {
  T trace=0;
  if((bandwidth==0)||(B.getBandwidth()==0))
    {
    VEC_ITER_TYPE vector<T>::iterator d1=diag.begin();
    VEC_ITER_TYPE vector<T>::iterator d2=B.getDiagIterator();
    unsigned i;

    for(i=0; i<dim; i++, d1++, d2++)
      {
      trace += *d1* *d2;
      }
    }
  else if(bandwidth>0 && B.getBandwidth()>0)
    {
    if(bandwidth==B.getBandwidth())
      {
      VEC_ITER_TYPE vector<T>::iterator d1=diag.begin();
      VEC_ITER_TYPE vector<T>::iterator d2=B.getDiagIterator();
      VEC_ITER_TYPE vector<T>::iterator env1=env.begin();
      VEC_ITER_TYPE vector<T>::iterator env2=B.getEnvIterator();
      for(; d1<diag.end(); d1++, d2++)
        {
        trace += *d1* *d2;
        }
      for(; env1<env.end(); env1++, env2++)
        {
        trace += 2* *env1* *env2;
        }
      }
    else if(bandwidth<B.getBandwidth())
      {
      VEC_ITER_TYPE vector<T>::iterator d1=diag.begin();
      VEC_ITER_TYPE vector<T>::iterator d2=B.getDiagIterator();
      VEC_ITER_TYPE vector<T>::iterator env1=env.begin();
      VEC_ITER_TYPE vector<T>::iterator env2=B.getEnvIterator();
      int i, k, diff,bbw;
      bbw=B.getBandwidth();
      diff=bbw-bandwidth;

      for(; d1<diag.end(); d1++, d2++)
        {
        trace += *d1* *d2;
        }
      for(i=0; i<bandwidth; i++)
        {
        for(k=0; (k<i) && (env1<env.end()) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      for(i=bandwidth; i<bbw; i++, env2+=(i-bandwidth))
        {
        for(k=0; (k<bandwidth) && (env1<env.end()) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      for(i=bandwidth; i<(int)dim; i++, env2+=diff)
        {
        for(k=0; (k<bandwidth) && (env1<env.end()) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      }
    else
      {
      VEC_ITER_TYPE vector<T>::iterator d1=diag.begin();
      VEC_ITER_TYPE vector<T>::iterator d2=B.getDiagIterator();
      VEC_ITER_TYPE vector<T>::iterator env1=env.begin();
      VEC_ITER_TYPE vector<T>::iterator env2=B.getEnvIterator();
      VEC_ITER_TYPE vector<T>::iterator env2end=B.getEnvIterator()+B.getXenv(dim);

      int i, k, diff, bbw;
      bbw=B.getBandwidth();
      diff=bandwidth-bbw;

      for(; d1<diag.end(); d1++, d2++)
        {
        trace += *d1* *d2;
        }
      for(i=0; i<bbw; i++)
        {
        for(k=0; (k<i) && (env2<env2end) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      for(i=bbw; i<bandwidth; i++, env1+=(i-bbw))
        {
        for(k=0; (k<bbw) && (env2<env2end) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      for(i=bandwidth; i<(int)dim; i++, env1+=diff)
        {
        for(k=0; (k<bbw) && (env2<env2end) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      }
    }
  else
    {
    VEC_ITER_TYPE vector<T>::iterator d1=diag.begin();
    VEC_ITER_TYPE vector<T>::iterator d2=B.getDiagIterator();
    VEC_ITER_TYPE vector<T>::iterator env1=env.begin();
    VEC_ITER_TYPE vector<T>::iterator env2=B.getEnvIterator();
    vector<unsigned>::iterator xenv1=xenv.begin();
    vector<unsigned>::iterator xenv2=B.getXenvIterator();

    unsigned i, k;

    for(; d1<diag.end(); d1++, d2++)
      {
      trace += *d1* *d2;
      }

    for(i=0; i<dim; i++, ++xenv1, ++xenv2)
      {
      if(*(xenv1+1)-*xenv1>=*(xenv2+1)-*xenv2)
        {
        for(k=0; k<*(xenv1+1)-*xenv1-*(xenv2+1)+*xenv2; k++)
          {
          ++env1;
          }
        for(k=0; k<*(xenv2+1)-*xenv2; k++)
          {
          trace += 2* *env1 * *env2;
          ++env1; ++env2;
          }
        }
      else
        {
        for(k=0; k<*(xenv2+1)-*xenv2-*(xenv1+1)+*xenv1; k++)
          {
          ++env2;
          }
        for(k=0; k<*(xenv1+1)-*xenv1; k++)
          {
          trace += 2* *env1 * *env2;
          ++env1; ++env2;
          }
        }
      }
    }
  return trace;
  }

template<class T>
vector<unsigned> envmatrix<T>::computeMaxXenv(const envmatrix<T> & B)
  {
  assert(dim==B.getDim());
  vector<unsigned> maxxenv(xenv.size(),0);
  unsigned i;

  vector<unsigned>::iterator bxe = B.getXenvIterator()+1;
  vector<unsigned>::iterator xe = xenv.begin()+1;
  vector<unsigned>::iterator maxxe = maxxenv.begin()+1;

  for(i=1; i<xenv.size(); i++, ++bxe, ++xe, ++maxxe)
    {
//    if(xenv[i]-xenv[i-1]>=B.getXenv(i)-B.getXenv(i-1))
    if(*xe-*(xe-1)>=*bxe-*(bxe-1))
      {
//      maxxenv[i]=maxxenv[i-1]+xenv[i]-xenv[i-1];
      *maxxe=*(maxxe-1)+*xe-*(xe-1);
      }
    else
      {
//      maxxenv[i]=maxxenv[i-1]+B.getXenv(i)-B.getXenv(i-1);
      *maxxe=*(maxxe-1)+*bxe-*(bxe-1);
      }
    }
  return maxxenv;
  }


//------------------------------------------------------------------------------
//---------- Functions to get pointers to elements of the matrix----------------
//------------------------------------------------------------------------------


template<class T>
  VEC_ITER_TYPE vector<T>::iterator envmatrix<T>::getDiagIterator()
  {
  return diag.begin();
  }

template<class T>
  VEC_ITER_TYPE vector<T>::iterator envmatrix<T>::getEnvIterator()
  {
  return env.begin();
  }

template<class T>
vector<unsigned>::iterator envmatrix<T>::getXenvIterator()
  {
  return xenv.begin();
  }


//------------------------------------------------------------------------------
//------------- Functions for changing elements of the matrix-------------------
//------------------------------------------------------------------------------


template<class T>
void envmatrix<T>::setDiag(const unsigned & i, const T & t)
  {
  diag[i]=t;
  decomposed=false;
  rational_decomposed=false;
  }

template<class T>
void envmatrix<T>::set(const unsigned & i, const unsigned & j, const T & t)
  {

  unsigned ih, jh;
  int kl, ku, zeroes;

  if(i>j)
    {
    ih = i;
    jh = j;
    }
  else if(i<j)
    {
    ih = j;
    jh = i;
    }
  else
    {
    diag[i]=t;
    decomposed=false;
    rational_decomposed=false;
    return;
    }

  kl=xenv[ih];
  ku=xenv[ih+1];
  zeroes=ih-ku+kl;

  if(jh>=zeroes)
    {
    env[kl+jh-zeroes]=t;
    decomposed=false;
    rational_decomposed=false;
    }

  }

template<class T>
void envmatrix<T>::setDecomposed(const bool &t)
  {
  decomposed=t;
  rational_decomposed=false;
  }

template<class T>
void envmatrix<T>::setRational_decomposed(const bool &t)
  {
  decomposed=false;
  rational_decomposed=t;
  }


//------------------------------------------------------------------------------
//-------------------- Functions for computing new matrices---------------------
//------------------------------------------------------------------------------


template<class T>
void envmatrix<T>::addtodiag(envmatrix &X, envmatrix &K,
                             const T &f1, const T &f2)
  {
  assert(X.getDim()==K.getDim());
  assert(dim==K.getDim());
  assert(xenv==K.getXenv());
//  unsigned i;
  VEC_ITER_TYPE vector<T>::iterator xd = X.getDiagIterator();
  VEC_ITER_TYPE vector<T>::iterator kd = K.getDiagIterator();
  VEC_ITER_TYPE vector<T>::iterator ke = K.getEnvIterator();
  VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
  VEC_ITER_TYPE vector<T>::iterator e = env.begin();

//  for(i=0; i<dim; i++)
  for(; d!=diag.end();++d, ++kd, ++xd)
    {
//    diag[i]=f1* X.getDiag(i) + f2* K.getDiag(i);
      *d = f1* *xd + f2* *kd;
    }
//  unsigned n = K.getXenv(dim);
//  for(i=0; i<n; i++)
  for(; e!=env.end(); ++ke, ++e)
    {
//    env[i]=f2* K.getEnv(i);
    *e = f2* *ke;
    }

  decomposed = false;
  rational_decomposed=false;
  }

template<class T>
void envmatrix<T>::addto(envmatrix &X, envmatrix &K,
                             const T &f1, const T &f2)
  {
  assert(X.getDim()==K.getDim());
  assert(dim==K.getDim());

  VEC_ITER_TYPE vector<T>::iterator xd = X.getDiagIterator();
  VEC_ITER_TYPE vector<T>::iterator kd = K.getDiagIterator();
  VEC_ITER_TYPE vector<T>::iterator d = diag.begin();
  for(; d!=diag.end(); ++d, ++kd, ++xd)
    {
    *d = f1* *xd + f2* *kd;
    }

  if(X.getBandwidth()>=0 && K.getBandwidth()>=0)
    {
    unsigned bwx = X.getBandwidth();
    unsigned bwk = K.getBandwidth();

    VEC_ITER_TYPE vector<T>::iterator ke = K.getEnvIterator();
    VEC_ITER_TYPE vector<T>::iterator xe = X.getEnvIterator();
    VEC_ITER_TYPE vector<T>::iterator e = env.begin();

    if(bwx>bwk)
      {
      assert(bandwidth==X.getBandwidth());
      unsigned i, k;
      for(i=0; i<bwk; i++)
        {
        for(k=0; k<i; k++, ++e, ++ke, ++xe)
          {
          *e = f1* *xe + f2* *ke;
          }
        }
      for(i=bwk; i<bwx; i++)
        {
        for(k=0; k<i-bwk; k++, ++e, ++xe)
          {
          *e = f1* *xe;
          }
        for(k=bwx-bwk; k<(unsigned)bandwidth; k++, ++e, ++xe, ++ke)
          {
          *e = f1* *xe + f2* *ke;
          }
        }
      for(i=bwx; i<dim; i++)
        {
        for(k=0; k<bwx-bwk; k++, ++e, ++xe)
          {
          *e = f1* *xe;
          }
        for(k=bwx-bwk; k<(unsigned)bandwidth; k++, ++e, ++xe, ++ke)
          {
          *e = f1* *xe + f2* *ke;
          }
        }
      }
    else if(bwx<bwk)
      {
      assert(bandwidth==K.getBandwidth());
      unsigned i, k;
      for(i=0; i<bwx; i++)
        {
        for(k=0; k<i; k++, ++e, ++ke, ++xe)
          {
          *e = f1* *xe + f2* *ke;
          }
        }
      for(i=bwx; i<bwk; i++)
        {
        for(k=0; k<i-bwx; k++, ++e, ++ke)
          {
          *e = f2* *ke;
          }
        for(k=bwk-bwx; k<(unsigned)bandwidth; k++, ++e, ++xe, ++ke)
          {
          *e = f1* *xe + f2* *ke;
          }
        }
      for(i=bwk; i<dim; i++)
        {
        for(k=0; k<bwk-bwx; k++, ++e, ++ke)
          {
          *e = f2* *ke;
          }
        for(k=bwk-bwx; k<(unsigned)bandwidth; k++, ++e, ++xe, ++ke)
          {
          *e = f1* *xe + f2* *ke;
          }
        }
      }
    else
      {
      assert(bandwidth==X.getBandwidth());
      for(; e!=env.end(); ++ke, ++e, ++xe)
        {
        *e = f1* *xe + f2* *ke;
        }
      }
    }
  else
    {
    VEC_ITER_TYPE vector<T>::iterator ke = K.getEnvIterator();
    VEC_ITER_TYPE vector<T>::iterator xe = X.getEnvIterator();
    VEC_ITER_TYPE vector<T>::iterator e = env.begin();
    vector<unsigned>::iterator kxe = K.getXenvIterator();
    vector<unsigned>::iterator xxe = X.getXenvIterator();

    unsigned i, k;
    for(i=0; i<dim; i++, ++kxe, ++xxe)
      {
//      if(K.getXenv(i+1)-K.getXenv(i)>=X.getXenv(i+1)-X.getXenv(i))
      if(*(kxe+1)-*kxe>=*(xxe+1)-*xxe)
        {
//        for(k=0; k<K.getXenv(i+1)-K.getXenv(i)-X.getXenv(i+1)+X.getXenv(i); k++)
        for(k=0; k<*(kxe+1)-*kxe-*(xxe+1)+*xxe; k++)
          {
          *e = f2* *ke;
          ++e; ++ke;
          }
//        for(k=0; k<X.getXenv(i+1)-X.getXenv(i); k++)
        for(k=0; k<*(xxe+1)-*xxe; k++)
          {
          *e = f1* *xe + f2* *ke;
          ++e; ++ke; ++xe;
          }
        }
      else
        {
//        for(k=0; k<X.getXenv(i+1)-X.getXenv(i)-K.getXenv(i+1)+K.getXenv(i); k++)
        for(k=0; k<*(xxe+1)-*xxe-*(kxe+1)+*kxe; k++)
          {
          *e = f1* *xe;
          ++e; ++xe;
          }
//        for(k=0; k<K.getXenv(i+1)-K.getXenv(i); k++)
        for(k=0; k<*(kxe+1)-*kxe; k++)
          {
          *e = f1* *xe + f2* *ke;
          ++e; ++ke; ++xe;
          }
        }
      }
    }
  decomposed=false;
  rational_decomposed=false;
  }

template<class T>
T envmatrix<T>::compute_quadform(const statmatrix<T> & x,
                                 const unsigned & c)
  {
  T sum = T(0);
  if(bandwidth==0)
    {
    unsigned i;
    VEC_ITER_TYPE vector<T>::iterator di = diag.begin();
    T* xi = x.getV()+c;
    unsigned d = x.cols();
//    for(i=0; i<dim; i++)
    for(i=0; i<dim; i++, ++di, xi+=d)
      {
//      sum += diag[i]*x(i,c)*x(i,c);
      sum += *di* *xi* *xi;
      }
    }
  else if(bandwidth==1)
    {
    unsigned i;
    unsigned d = x.cols();
    VEC_ITER_TYPE vector<T>::iterator di = diag.begin();
    VEC_ITER_TYPE vector<T>::iterator e = env.begin();
    T* xi = x.getV()+c;
//    sum += diag[0]*x(0,c)*x(0,c);
    sum += *di**xi**xi;
    ++di; xi+=d;
//    for(i=1; i<dim; i++)
    for(i=1; i<dim; i++, ++e, ++di, xi+=d)
      {
      if(*e!=0)
        {
//        sum += 2*env[i-1]*x(i,c)*x(i-1,c);
        sum += 2**e**xi**(xi-d);
        }
//      sum += diag[i]*x(i,c)*x(i,c);
      sum += *di**xi**xi;
      }
    }
  else if(bandwidth==2)
    {
//    unsigned i, k;
    unsigned i;
    unsigned d=x.cols();
    VEC_ITER_TYPE vector<T>::iterator di = diag.begin();
    VEC_ITER_TYPE vector<T>::iterator e = env.begin();
    T* xi = x.getV()+c;

//    sum += diag[0]*x(0,c)*x(0,c);
    sum += *di**xi**xi;
    ++di; xi+=d;
//    sum += 2*env[0]*x(1,c)*x(0,c);
    sum += 2**e**xi**(xi-d);
    ++e;
//    sum += diag[1]*x(1,c)*x(1,c);
    sum += *di**xi**xi;
    xi+=d; ++di;
//    for(i=2, k=1; i<dim; i++, k+=2)
    for(i=2; i<dim; i++, ++e, xi+=d, ++di)
      {
//      sum += 2*env[k]*x(i,c)*x(i-2,c);
      sum += 2**e**xi**(xi-2*d);
      e++;
//      sum += 2*env[k+1]*x(i,c)*x(i-1,c);
      sum += 2**e**xi**(xi-d);
//      sum += diag[i]*x(i,c)*x(i,c);
      sum += *di**xi**xi;
      }
    }
  else if(bandwidth>2)
    {
//    unsigned i, k, j;
    unsigned i, k;
    unsigned d = x.cols();
    VEC_ITER_TYPE vector<T>::iterator di = diag.begin();
    VEC_ITER_TYPE vector<T>::iterator e = env.begin();
    T* xi = x.getV()+c;
    T* xk;

//    for(i=0; i<bandwidth; i++)
    for(i=0; i<(unsigned)bandwidth; i++, ++di, xi+=d)
      {
//      for(k=0, j=xenv[i]; k<i; k++, j++)
      for(k=0, xk = x.getV()+c; k<i; k++, ++e, xk+=d)
        {
//        sum+=2*env[j]*x(i,c)*x(k,c);
        sum+=2**e**xi**xk;
        }
//      sum += diag[i]*x(i,c)*x(i,c);
      sum += *di**xi**xi;
      }
//    for(i=bandwidth; i<dim; i++)
    for(i=bandwidth; i<dim; i++, ++di, xi+=d)
      {
//      for(k=i-bandwidth, j=xenv[i]; k<i; k++, j++)
      for(k=i-bandwidth, xk=x.getV()+c+k*d; k<i; k++, ++e, xk+=d)
        {
//        sum+=2*env[j]*x(i,c)*x(k,c);
        sum+=2**e**xi**xk;
        }
//      sum += diag[i]*x(i,c)*x(i,c);
      sum += *di**xi**xi;
     }
    }
  else
    {
//    unsigned i,j,jstop,k;
    unsigned i,j,jstop;
    unsigned d = x.cols();
    VEC_ITER_TYPE vector<T>::iterator di = diag.begin();
    VEC_ITER_TYPE vector<T>::iterator ej;
    vector<unsigned>::iterator xei1 = xenv.begin();
    vector<unsigned>::iterator xei2 = xenv.begin()+1;
    T* xi = x.getV()+c;
    T* xk;

//    for(i=0; i<dim; i++)
    for(i=0; i<dim; i++, ++xei1, ++xei2, xi+=d, ++di)
      {
//      if(xenv[i+1]-xenv[i]>0)
     if(*xei2-*xei1>0)
        {
//        jstop = xenv[i+1];
        jstop = *xei2;
//        k = i-*xei2+*xei1;
        xk = x.getV()+c+(i-*xei2+*xei1)*d;
//        for(j=xenv[i]; j<jstop; j++)
        for(j=*xei1, ej= env.begin()+*xei1; j<jstop; j++, xk+=d, ++ej)
          {
//          sum += 2*env[j]*x(i,c)*x(k,c);
//          k++;
          if(*ej!=0)
            {
            sum += 2* *ej* *xi* *xk;
            }
          }
        }
//      sum += diag[i]*x(i,c)*x(i,c);
      sum += *di* *xi* *xi;
      }
    }
  return sum;
  }

template<class T>
T envmatrix<T>::compute_sumfabsdiff(const statmatrix<T> & x,
                                 const unsigned & c)
  {
  T sum = T(0);
  unsigned i,j,jstop;
  unsigned d = x.cols();
  VEC_ITER_TYPE vector<T>::iterator ej;
  vector<unsigned>::iterator xei1 = xenv.begin();
  vector<unsigned>::iterator xei2 = xenv.begin()+1;
  T* xi = x.getV()+c;
  T* xk;

  for(i=0; i<dim; i++, ++xei1, ++xei2, xi+=d)
    {
   if(*xei2-*xei1>0)
      {
      jstop = *xei2;
      xk = x.getV()+c+(i-*xei2+*xei1)*d;
      for(j=*xei1, ej= env.begin()+*xei1; j<jstop; j++, xk+=d, ++ej)
        {
        if(*ej!=0)
          {
          sum -= *ej * fabs(*xi - *xk);
          }
        }
      }
    }
  return sum;
  }

template<class T>
T envmatrix<T>::compute_quadformblock(const statmatrix<T> & x,
             const unsigned & c, const unsigned & a, const unsigned & b)
  {
  assert(bandwidth>=0);

  T sum = T(0);
  if(bandwidth==0)
    {
    unsigned i;
//    unsigned d=x.cols();
    for(i=a; i<b+1; i++)
      {
      sum += diag[i]*x(i,c)*x(i,c);
      }
    }
  else if(bandwidth==1)
    {
    unsigned i;
//    unsigned d=x.cols();
    sum += diag[a]*x(a,c)*x(a,c);
    for(i=a+1; i<b+1; i++)
      {
      if(env[i-1]!=0)
        {
        sum += 2*env[i-1]*x(i,c)*x(i-1,c);
        }
      sum +=  diag[i]*x(i,c)*x(i,c);
      }
    }
  else if(bandwidth>1)
    {
    unsigned i,k,j;
//    unsigned d=x.cols();
    unsigned helpbw;
    if(b-a+1<(unsigned)bandwidth)
      {
      helpbw=b-a+1;
      }
    else
      {
      helpbw=bandwidth;
      }
    for(i=0; i<helpbw; i++)
      {
      sum += diag[i+a]*x(i+a,c)*x(i+a,c);
      for(k=0, j=xenv[i+a+1]-i; k<i; k++, j++)
        {
        sum += 2*env[j]*x(i+a,c)*x(k+a,c);
        }
      }
    if(helpbw==(unsigned)bandwidth)
      {
      for(i=bandwidth; i<b-a+1; i++)
        {
        sum += diag[i+a]*x(i+a,c)*x(i+a,c);
        for(k=i-bandwidth, j=xenv[i+a]; k<i; k++, j++)
          {
          sum += 2*env[j]*x(i+a,c)*x(k+a,c);
          }
        }
      }
    }
  return sum;
  }


//------------------------------------------------------------------------------
//----------------- Functions for printing matrices-----------------------------
//------------------------------------------------------------------------------

template<class T>
void envmatrix<T>::print1(ostream & out)
  {
  unsigned i;
  int k, kl, ku;
  for(i=0; i<dim; i++)
    {
    kl=xenv[i];
    ku=xenv[i+1];
    for(k=kl; k<ku; k++)
      {
      out << env[k] <<" \t";
      }
    out << diag[i] << "\n";
    }
  out << endl;
  }


template<class T>
void envmatrix<T>::print2(ostream & out)
  {
  unsigned i;
  int kl, ku, zeroes, k;
  for(i=0; i<dim; i++)
    {
    kl=xenv[i];
    ku=xenv[i+1];
    zeroes=i-ku+kl;
    for(k=0; k<zeroes; k++)
      {
      out << "0 \t";
      }
    for(k=kl; k<ku; k++)
      {
      out << env[k] <<" \t";
      }
    out << diag[i] << "\n";
    }
  out << endl;
  }

template<class T>
void envmatrix<T>::print3(ostream & out)
  {
  int i;
  int n=xenv[dim];
  out << "env: ";
  for(i=0; i<n; i++)
    {
    out << env[i] << " ";
    }
  out << "\nxenv: ";
  for(i=0; i<dim+1; i++)
    {
    out << xenv[i] << " ";
    }
  out << "\ndiag: ";
  for(i=0; i<dim; i++)
    {
    out << diag[i] << " ";
    }
  out << "\n";
  out << "bandwidth: " << bandwidth <<endl;
  out << "decomposed: " << decomposed << endl;
  if(decomposed)
    {
    out << "lenv: ";
    for(i=0; i<n; i++)
      {
      out << lenv[i] << " ";
      }
    out << "\nldiag: ";
    for(i=0; i<dim; i++)
      {
      out << ldiag[i] << " ";
      }
    }
  out << endl;
  }

template<class T>
void envmatrix<T>::print4(ostream & out)
  {
  unsigned i;
  int kl, ku, k;
  datamatrix help = datamatrix(dim,dim,0);
  for(i=0; i<dim; i++)
    {
    kl=xenv[i];
    ku=xenv[i+1];
    for(k=kl; k<ku; k++)
      {
      help(i,i-ku+k) = env[k];
      help(i-ku+k,i) = env[k];
      }
    help(i,i) = diag[i];
    }
  help.prettyPrint(out);
  }

template<class T>
void envmatrix<T>::print4L(ostream & out)
  {
  unsigned i;
  int kl, ku, k;
  datamatrix help = datamatrix(dim,dim,0);
  for(i=0; i<dim; i++)
    {
    kl=xenv[i];
    ku=xenv[i+1];
    for(k=kl; k<ku; k++)
      {
      help(i,i-ku+k) = lenv[k];
//      help(i-ku+k,i) = lenv[k];
      }
    help(i,i) = ldiag[i];
    }
  help.prettyPrint(out);
  }

//------------------------------------------------------------------------------
//------------------------------- old functions --------------------------------
//------------------------------------------------------------------------------

/*template<class T>
void envmatrix<T>::decomp(envmatrix<T> & L)
  {
  if(!decomposed)
    {
    if(!isDiag)
      {
      assert(diag[0]>0);
      L.setDiag(0,sqrt(diag[0]));
      if(dim>1)
        {
        unsigned i;
        for(i=1; i<dim; i++)
          {
          unsigned ixenv = xenv[i];
          unsigned iband = xenv[i+1]-ixenv;
          T temp = diag[i];
          if(iband>0)
            {
            unsigned ifirst = i-iband;
            L.solveL2(env, xenv[i], ifirst, iband);
            unsigned jstop = xenv[i+1]-1;
            unsigned j;
            for(j=ixenv; j<=jstop; j++)
              {
              T s = L.getEnv(j);
              temp = temp - s*s;
              }
            }
          assert(temp>0);
          L.setDiag(i,sqrt(temp));
          }
        }
      }
    else
      {
      unsigned i;
      for(i=0; i<dim; i++)
        {
        L.setDiag(i,sqrt(diag[i]));
        }
      }
    }
  L.setDecomposed(true);
  }

template<class T>
void envmatrix<T>::solveL2(const vector<T> & aenv, const unsigned & aoffset,
                           const unsigned & offset, const unsigned & neqns)
  {
  unsigned last = aoffset;
  unsigned i;
  for(i=0; i<neqns; i++)
    {
    unsigned iband = xenv[i+offset+1]-xenv[i+offset];
    if(iband > i)
      {
      iband = i;
      }
    T s=aenv[i+aoffset];
    unsigned l = i+aoffset-iband;
    if((iband>0) && (last>=l))
      {
      unsigned kstop=xenv[i+offset+1]-1;
      unsigned k;
      for(k=xenv[i+offset+1]-iband; k<=kstop; k++)
        {
        s=s-env[k]*env[l];
        l++;
        }
      }
    if(s!=0)
      {
      env[i+aoffset]=s/diag[i+offset];
      last=i+aoffset;
      }
    else
      {
      env[i+aoffset]=0;
      }
    }
  }

  template<class T>
void envmatrix<T>::solveL1(const vector<T> & b,vector<T> & res,
                           const unsigned & offset, const unsigned & neqns)
  {
  unsigned ifirst = 0;
  while((ifirst<neqns) && (b[ifirst]==0))
    {
    res[ifirst]=0;
    ifirst++;
    }
  unsigned last = 0;
  unsigned i;
  for(i=ifirst; i<neqns; i++)
    {
    unsigned iband = xenv[i+offset+1]-xenv[i+offset];
    if(iband > i)
      {
      iband = i;
      }
    T s=b[i];
    unsigned l = i-iband;
    if((iband>0) && (last>=l))
      {
      unsigned kstop=xenv[i+offset+1]-1;
      unsigned k;
      for(k=xenv[i+offset+1]-iband; k<=kstop; k++)
        {
        s=s-env[k]*res[l];
        l++;
        }
      }
    if(s!=0)
      {
      res[i]=s/diag[i+offset];
      last=i;
      }
    else
      {
      res[i]=0;
      }
    }
  }*/






