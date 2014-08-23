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





#include "envmatrix_penalty.h"

//------------------------------------------------------------------------------
//----------------- Functions for computing penalty matrices -------------------
//------------------------------------------------------------------------------

envmatrix<double> Kmrfenv(const MAP::map & m)
  {
  unsigned n = m.get_nrregions();
  vector<double> diag(n,0);
  vector<double> env;
  vector<unsigned> xenv(1,0);
  vector<unsigned> ne;
  unsigned i,k,l,minne,c;
  vector< vector<double> > w = m.get_weights();
  for(i=0; i<n; i++)
    {
    diag[i]=m.get_weightssum(i);
    ne = m.get_neighbors()[i];
    c = ne.size();
    sort(ne.begin(), ne.end());
    minne = *ne.begin();
    if(minne<i)
      {
      l=env.size();
      env.resize(l+i-minne,0);
      for(k=0; (ne[k]<i) && (k<c); k++)
        {
        env[l+ne[k]-minne]=-m.get_weights()[i][k];
        }
      xenv.push_back(xenv[i]+i-minne);
      }
    else
      {
      xenv.push_back(xenv[i]);
      }
    }
  return envmatrix<double>(env, diag, xenv);
  }

envmatrix<double> Krw1env(const vector<double> & weight)
  {
  unsigned d = weight.size();
  vector<double> diag(d,0);
  vector<double> env(d-1,0);
  vector<unsigned> xenv(d+1,0);

  vector<double>::iterator e = env.begin();
  vector<double>::iterator di = diag.begin();
  vector<unsigned>::iterator x = xenv.begin()+2;

  unsigned i;

//  diag[0] = 1.0/weight[1];
  *di = 1.0/weight[1];
  ++di;
//  env[0]=-1.0/weight[1];
  *e=-1.0/weight[1];
  ++e;
  for(i=1; i<d-1; i++, ++di, ++e, ++x)
    {
//    diag[i] = 1.0/weight[i]+1.0/weight[i+1];
//    env[i] = -1.0/weight[i+1];
//    xenv[i+1] = i;
    *di = 1.0/weight[i]+1.0/weight[i+1];
    *e = -1.0/weight[i+1];
    *x = i;
    }
//  diag[d-1] = 1.0/weight[d-1];
  *di = 1.0/weight[d-1];
//  xenv[d] = d-1;
  *x = d-1;

  return envmatrix<double>(env, diag, xenv, 1);
  }


envmatrix<double> Krw1env(const datamatrix & weight)
  {
  unsigned d = weight.rows();
  vector<double> diag(d,0);
  vector<double> env(d-1,0);
  vector<unsigned> xenv(d+1,0);

  vector<double>::iterator e = env.begin();
  vector<double>::iterator di = diag.begin();
  vector<unsigned>::iterator x = xenv.begin()+2;

  unsigned i;

  *di = 1.0/weight(1,0);
  ++di;
  *e=-1.0/weight(1,0);
  ++e;
  for(i=1; i<d-1; i++, ++di, ++e, ++x)
    {
    *di = 1.0/weight(i,0)+1.0/weight(i+1,0);
    *e = -1.0/weight(i+1,0);
    *x = i;
    }

  *di = 1.0/weight(d-1,0);
  *x = d-1;

  return envmatrix<double>(env, diag, xenv, 1);
  }


envmatrix<double> Krw2env(const vector<double> & weight)
  {
  unsigned d = weight.size();
  vector<double> diag(d,0);
  vector<double> env(2*d-3,0);
  vector<unsigned> xenv(d+1,0);

  vector<double>::iterator e = env.begin();
  vector<double>::iterator di = diag.begin();
  vector<unsigned>::iterator x = xenv.begin()+2;

//  diag[0] = weight[2]*weight[2]/
//            (weight[1]*weight[1]*weight[2]*(1+weight[2]/weight[1]));
//  diag[0] = weight[2]/
  *di = weight[2]/
            (weight[1]*(weight[1]+weight[2]));
  ++di;
//  diag[1] = (1+weight[2]/weight[1])*(1+weight[2]/weight[1])/
  *di = (1+weight[2]/weight[1])*(1+weight[2]/weight[1])/
            (weight[2]*(1+weight[2]/weight[1]))+
//            weight[3]*weight[3]/
//            (weight[2]*weight[2]*weight[3]*(1+weight[3]/weight[2]));
            weight[3]/
            (weight[2]*(weight[2]+weight[3]));
  ++di;
//  env[0] = -(1+weight[2]/weight[1])*weight[2]/
//            (weight[1]*weight[2]*(1+weight[2]/weight[1]));
//  env[0] = -(1+weight[2]/weight[1])/
  *e = -(1+weight[2]/weight[1])/
            (weight[1]+weight[2]);
  ++e;
//  env[1] = weight[2]/(weight[1]*weight[2]*(1+weight[2]/weight[1]));
//  env[1] = 1/(weight[1]+weight[2]);
  *e = 1/(weight[1]+weight[2]);
  ++e;
//  env[2] = -(1+weight[2]/weight[1])/(weight[2]*(1+weight[2]/weight[1]))+
  *e = -(1+weight[2]/weight[1])/(weight[2]*(1+weight[2]/weight[1]))+
//           -(1+weight[3]/weight[2])*weight[3]/(weight[2]*weight[3]*(1+weight[3]/weight[2]));
           -(1+weight[3]/weight[2])/(weight[2]+weight[3]);
  ++e;
//  env[3]= weight[3]/(weight[2]*weight[3]*(1+weight[3]/weight[2]));
//  env[3]= 1/(weight[2]+weight[3]);
  *e= 1/(weight[2]+weight[3]);
  ++e;
//  xenv[2]=1;
  *x = 1;
  ++x;
  unsigned i;
//  for(i=2; i<d-2; i++)
  for(i=2; i<d-2; i++, ++x, ++di, ++e)
    {
//    xenv[i+1] = 2*i-1;
    *x = 2*i-1;
//    diag[i] = (1+weight[i+1]/weight[i])*(1+weight[i+1]/weight[i])/
    *di = (1+weight[i+1]/weight[i])*(1+weight[i+1]/weight[i])/
              (weight[i+1]*(1+weight[i+1]/weight[i]))+
//              weight[i+2]*weight[i+2]/
//              (weight[i+1]*weight[i+1]*weight[i+2]*(1+weight[i+2]/weight[i+1]))+
              weight[i+2]/
              (weight[i+1]*(weight[i+1]+weight[i+2]))+
              1/(weight[i]*(1+weight[i]/weight[i-1]));
//    env[2*i+1] = weight[i+2]/
//                 (weight[i+1]*weight[i+2]*(1+weight[i+2]/weight[i+1]));
//    env[2*i] = -(1+weight[i+1]/weight[i])/
    *e = -(1+weight[i+1]/weight[i])/
               (weight[i+1]*(1+weight[i+1]/weight[i]))+
//               -(1+weight[i+2]/weight[i+1])*weight[i+2]/
//               (weight[i+1]*weight[i+2]*(1+weight[i+2]/weight[i+1]));
               -(1+weight[i+2]/weight[i+1])/
               (weight[i+1]+weight[i+2]);
    ++e;
//    env[2*i+1] = 1/(weight[i+1]+weight[i+2]);
    *e = 1/(weight[i+1]+weight[i+2]);
    }
//  diag[d-2] = (1+weight[d-1]/weight[d-2])*(1+weight[d-1]/weight[d-2])/
  *di = (1+weight[d-1]/weight[d-2])*(1+weight[d-1]/weight[d-2])/
              (weight[d-1]*(1+weight[d-1]/weight[d-2]))+
              1/(weight[d-2]*(1+weight[d-2]/weight[d-3]));
  ++di;
//  diag[d-1] = 1/(weight[d-1]*(1+weight[d-1]/weight[d-2]));
  *di = 1/(weight[d-1]*(1+weight[d-1]/weight[d-2]));
//  env[2*d-4] = -(1+weight[d-1]/weight[d-2])/
  *e = -(1+weight[d-1]/weight[d-2])/
               (weight[d-1]*(1+weight[d-1]/weight[d-2]));
//  xenv[d-1] = 2*d-5;
  *x = 2*d-5;
  ++x;
//  xenv[d] = 2*d-3;
  *x = 2*d-3;

  return envmatrix<double>(env, diag, xenv, 2);
  }


envmatrix<double> Krw2env(const unsigned & nrpar)
  {

  unsigned i;


  datamatrix D2(nrpar-2,nrpar,0);

  for (i=0;i<D2.rows();i++)
    {
    D2(i,i) = 1;
    D2(i,i+1) = -2;
    D2(i,i+2) = 1;
    }


  datamatrix K2 = D2.transposed()*D2;

  return envmatrix<double>(K2);

  }


envmatrix<double> Krw3env(const unsigned & nrpar)
  {

  unsigned i;

  datamatrix D1(nrpar-3,nrpar-2,0);

  for (i=0;i<D1.rows();i++)
    {
    D1(i,i) = -1;
    D1(i,i+1) = 1;
    }

  datamatrix D2(nrpar-2,nrpar,0);

  for (i=0;i<D2.rows();i++)
    {
    D2(i,i) = 1;
    D2(i,i+1) = -2;
    D2(i,i+2) = 1;
    }

  datamatrix D3;

  D3 = D1*D2;

  datamatrix K3 = D3.transposed()*D3;

  return envmatrix<double>(K3);

  }


envmatrix<double> Kseasonenv(const unsigned & per,const unsigned & s)
  {
  assert(s > 2*(per-1));
  vector<double> diag(s,per);
  vector<double> env((per-1)*s-(per-1)*per/2,1);
  vector<unsigned> xenv(s+1,0);

  vector<double>::iterator e = env.begin();
  vector<double>::iterator d = diag.begin();

//  int i,k;,l;
  unsigned i,k;

//  for(i=0,l=0; i<per; i++)
  for(i=0; i<per; i++, ++d)
    {
    xenv[i+1]=xenv[i]+i;
//    diag[i]=i+1;
    *d = i+1;
    diag[s-i-1]=i+1;
//    for(k=0;k<i;k++,l++)
    for(k=0;k<i;k++,++e)
      {
//      env[l]=k+1;
      *e = k+1;
      }
    }
  for(i=per; i<s-per+1; i++)
    {
    xenv[i+1]=xenv[i]+per-1;
//    diag[i]=per;
//    *d=per;
//    for(k=1;k<per;k++,l++)
    for(k=1;k<per;k++,++e)
      {
//      env[l]=k;
      *e = k;
      }
    }
  for(i=s-per+1;i<s;i++)
    {
    xenv[i+1]=xenv[i]+per-1;
//    for(k=1;k<s-i;k++,l++)
    for(k=1;k<s-i;k++,++e)
      {
//      env[l]=k;
      *e = k;
      }
//    for(k=s-i;k<per;k++,l++)
    for(k=s-i;k<per;k++,++e)
      {
//      env[l]=s-i;
      *e = s-i;
      }
    }
  return envmatrix<double>(env, diag, xenv, per-1);
  }


envmatrix<double> Krw0env(const unsigned & nrpar)
  {
  vector<double> diag(nrpar,1);
  vector<double> env;
  vector<unsigned> xenv(nrpar+1,0);

  return envmatrix<double>(env, diag, xenv, 0);
  }
