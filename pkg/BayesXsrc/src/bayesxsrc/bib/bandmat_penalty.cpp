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





#include "bandmat_penalty.h"

//------------------------------------------------------------------------------
//----------------- Functions for computing penalty matrices -------------------
//------------------------------------------------------------------------------


bandmatdouble Krw1band(const vector<double> & weight)
  {

  unsigned S = weight.size();

  datamatrix diagelem(S,1);
  datamatrix upperdiagelem(S,1);

  unsigned i;
  for (i=1;i<S-1;i++)
    {
    diagelem(i,0) = 1.0/weight[i]+1.0/weight[i+1];
    upperdiagelem(i,0) = -1.0/(weight[i+1]);
    }

  diagelem(0,0) = 1.0/weight[1];
  diagelem(S-1,0) = 1.0/weight[S-1];
  upperdiagelem(0,0) = -1.0/weight[1];


  return bandmatdouble(diagelem,upperdiagelem);

  }


bandmatdouble Krw2band(const vector<double> & weight)
  {
  unsigned i;
  unsigned S = weight.size();

  datamatrix F(S-2,S,0);
  for (i=0;i<F.rows();i++)
    {
    F(i,i)   = weight[2+i]/weight[1+i];
    F(i,i+1) = -(1+weight[2+i]/weight[1+i]);
    F(i,i+2) = 1;
    }
  datamatrix Q(S-2,S-2,0);
  for(i=0;i<Q.rows();i++)
	 Q(i,i) =   weight[2+i]*(1+weight[2+i]/weight[1+i]);
//     Q(i,i) = 1;

  datamatrix K = F.transposed()*Q.inverse()*F;

  datamatrix diagelem(S,1);

  for (i=0;i<S;i++)
    diagelem(i,0) = K(i,i);

  datamatrix upperdiagelem(S,2);

  for (i=0;i<S-2;i++)
    {
    upperdiagelem(i,0) = K(i,i+1);
    upperdiagelem(i,1) = K(i,i+2);
    }
  upperdiagelem(S-2,0) = K(S-2,S-1);

  return bandmatdouble(diagelem,upperdiagelem);

  }


bandmatdouble Kseasonband(const unsigned & per,const unsigned & s)
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

  datamatrix diagelem(s,1);

  for (i=0;i<s;i++)
    diagelem(i,0) = K(i,i);

  datamatrix upperdiagelem(s,k);

  int p;

  for (i=0;i<s;i++)
    {
    if (i+k > s-1)
      p = s-1;
    else
      p = i+k;
    for (j=i+1;j<=p;j++)
      upperdiagelem(i,j-i-1) = K(i,j);
    }

  return bandmatdouble(diagelem,upperdiagelem);

  }


bandmatdouble Kmrfband(const MAP::map & m)
  {
  unsigned bs = m.get_bandsize();
  unsigned r = m.get_nrregions();

  datamatrix de(r,1);
  datamatrix ud(r,bs,0);
  unsigned i,j;
  unsigned ind;
  for (i=0;i<r;i++)
    {

    de(i,0) = m.get_weightssum(i);
    for(j=0;j<m.get_neighbors()[i].size();j++)
      {
      ind = m.get_neighbors()[i][j];
      if (ind > i)
      ud(i,ind-i-1) = -m.get_weights()[i][j];
      }

    }

  return bandmatdouble(de,ud);
  }


bandmatdouble Kmrflinearband(const unsigned & nr1,const unsigned & nr2)
  {
  datamatrix de(nr1*nr2,1,2);
  datamatrix ud(nr1*nr2,1,-1);

  unsigned i;
  for(i=0;i<nr1;i++)
    {
    de(i*nr2,0) = 1;
    de(i*nr2+nr2-1,0) = 1;
    ud((i+1)*nr2-1,0) = 0;
    }

  return bandmatdouble(de,ud);
  }

