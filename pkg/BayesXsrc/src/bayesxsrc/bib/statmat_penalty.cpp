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



#include "statmat_penalty.h"

using std::ifstream;

namespace STATMAT_PENALTY
{

statmatrix<double> Kmrf(const MAP::map & m)
  {
  unsigned r,i,j,n;
  r = m.get_nrregions();

  statmatrix<double>res(r,r,0);

  for(i=0; i<r; i++)
    {
    res(i,i) = m.get_weightssum(i);
    for(j=0; j<m.get_neighbors()[i].size(); j++)
      {
      n=m.get_neighbors()[i][j];
      res(i,n) = -m.get_weights()[i][j];
      res(n,i) = res(i,n);
      }
    }

  return res;
  }

statmatrix<double> K2dim_pspline(const unsigned & nknots)
  {
  statmatrix<double> res(nknots*nknots,nknots*nknots,0);
  unsigned i,j;

// 4 corners:
  res(0,0)=2;
  res(0,1)=-1;
  res(0,nknots)=-1;

  res(nknots-1,nknots-1)=2;
  res(nknots-1,nknots-2)=-1;
  res(nknots-1,2*nknots-1)=-1;

  res((nknots-1)*nknots,(nknots-1)*nknots)=2;
  res((nknots-1)*nknots,(nknots-2)*nknots)=-1;
  res((nknots-1)*nknots,(nknots-1)*nknots+1)=-1;

  res(nknots*nknots-1,nknots*nknots-1)=2;
  res(nknots*nknots-1,(nknots-1)*nknots-1)=-1;
  res(nknots*nknots-1,nknots*nknots-2)=-1;

// upper and lower edge

  for(i=1; i<nknots-1; i++)
    {
    res(i,i)=3;
    res(i,i-1)=-1;
    res(i,i+1)=-1;
    res(i,i+nknots)=-1;
    }

  for(i=(nknots-1)*nknots+1; i<nknots*nknots-1; i++)
    {
    res(i,i)=3;
    res(i,i-1)=-1;
    res(i,i+1)=-1;
    res(i,i-nknots)=-1;
    }

// rest

  for(i=2; i<nknots; i++)
    {
    res((i-1)*nknots,(i-1)*nknots)=3;
    res((i-1)*nknots,(i-2)*nknots)=-1;
    res((i-1)*nknots,(i-1)*nknots+1)=-1;
    res((i-1)*nknots,i*nknots)=-1;

    for(j=1; j<nknots-1; j++)
      {
      res((i-1)*nknots+j,(i-1)*nknots+j)=4;
      res((i-1)*nknots+j,(i-2)*nknots+j)=-1;
      res((i-1)*nknots+j,(i-1)*nknots+j-1)=-1;
      res((i-1)*nknots+j,(i-1)*nknots+j+1)=-1;
      res((i-1)*nknots+j,i*nknots+j)=-1;
      }

    res(i*nknots-1,i*nknots-1)=3;
    res(i*nknots-1,(i-1)*nknots-1)=-1;
    res(i*nknots-1,i*nknots-2)=-1;
    res(i*nknots-1,(i+1)*nknots-1)=-1;
    }
  return res;
  }

//}


statmatrix<double> K2dim_pspline_rw2(const unsigned & nknots, const unsigned & ox, const unsigned & oy)
  {
  statmatrix<double> res(nknots*nknots,nknots*nknots,0);
//  unsigned i,j;

  statmatrix<double>  Dx = diffmat(ox,nknots);       // difference matrix of order ox
  statmatrix<double> DDx = Dx.transposed() * Dx;     // D'D

  statmatrix<double>  Dy = diffmat(oy,nknots);       // difference matrix of order oy
  statmatrix<double> DDy = Dy.transposed() * Dy;

  statmatrix<double> I = statmatrix<double>::diag(nknots,1);  // identity matrix

  statmatrix<double> Px = kronecker(I,DDx);      // penalty matrix of x-direction
  statmatrix<double> Py = kronecker(DDy,I);      // penalty matrix of y-direction
  res = Px + Py;
  return res;
  }

statmatrix<double> K2dim_pspline_biharmonic(const unsigned & nknots)
  {
  unsigned i,j;
  statmatrix<double> res(nknots*nknots,nknots*nknots,0);

// 4 corners

// upper left
res(0,0) = 4;
res(0,1) = res(0,nknots) = -4;
res(0,2) = res(0,2*nknots) = 1;
res(0,nknots+1) = 2;

// upper right
res(nknots-1,nknots-1) = 4;
res(nknots-1,nknots-2) = res(nknots-1,2*nknots-1) = -4;
res(nknots-1,nknots-3) = res(nknots-1,3*nknots-1) = 1;
res(nknots-1,2*nknots-2) = 2;

// lower left
res(nknots*(nknots-1),nknots*(nknots-1)) = 4;
res(nknots*(nknots-1),nknots*(nknots-1)+1) = res(nknots*(nknots-1),nknots*(nknots-2)) = -4;
res(nknots*(nknots-1),nknots*(nknots-1)+2) = res(nknots*(nknots-1),nknots*(nknots-3)) = 1;
res(nknots*(nknots-1),nknots*(nknots-2)+1) = 2;

// lower right
res(nknots*nknots-1,nknots*nknots-1) = 4;
res(nknots*nknots-1,nknots*nknots-2) = res(nknots*nknots-1,nknots*(nknots-1)-1) = -4;
res(nknots*nknots-1,nknots*nknots-3) = res(nknots*nknots-1,nknots*(nknots-2)-1) = 1;
res(nknots*nknots-1,nknots*(nknots-1)-2) = 2;

// 4 second order corners

// upper left
res(nknots+1,nknots+1) = 18;
res(nknots+1,nknots+2) = res(nknots+1,2*nknots+1) = -8;
res(nknots+1,nknots) = res(nknots+1,1) = -6;
res(nknots+1,0) = res(nknots+1,2) = res(nknots+1,2*nknots) = res(nknots+1,2*nknots+2) = 2;
res(nknots+1,nknots+3) = res(nknots+1,3*nknots+1) = 1;

// upper right
res(2*nknots-2,2*nknots-2) = 18;
res(2*nknots-2,2*nknots-3) = res(2*nknots-2,3*nknots-2) = -8;
res(2*nknots-2,2*nknots-1) = res(2*nknots-2,nknots-2) = -6;
res(2*nknots-2,nknots-3) = res(2*nknots-2,nknots-1) = res(2*nknots-2,3*nknots-3) = res(2*nknots-2,3*nknots-1) = 2;
res(2*nknots-2,4*nknots-2) = res(2*nknots-2,2*nknots-4) = 1;

// lower left
res(nknots*(nknots-2)+1,nknots*(nknots-2)+1) = 18;
res(nknots*(nknots-2)+1,nknots*(nknots-2)+2) = res(nknots*(nknots-2)+1,nknots*(nknots-3)+1) = -8;
res(nknots*(nknots-2)+1,nknots*(nknots-2)) = res(nknots*(nknots-2)+1,nknots*(nknots-1)+1) = -6;
res(nknots*(nknots-2)+1,nknots*(nknots-3)) = res(nknots*(nknots-2)+1,nknots*(nknots-3)+2) = res(nknots*(nknots-2)+1,nknots*(nknots-1)) = res(nknots*(nknots-2)+1,nknots*(nknots-1)+2) = 2;
res(nknots*(nknots-2)+1,nknots*(nknots-2)+3) = res(nknots*(nknots-2)+1,nknots*(nknots-4)+1) = 1;

// lower right
res(nknots*(nknots-1)-2,nknots*(nknots-1)-2) = 18;
res(nknots*(nknots-1)-2,nknots*(nknots-1)-3) = res(nknots*(nknots-1)-2,nknots*(nknots-2)-2) = -8;
res(nknots*(nknots-1)-2,nknots*(nknots-1)-1) = res(nknots*(nknots-1)-2,nknots*nknots-2) = -6;
res(nknots*(nknots-1)-2,nknots*(nknots-2)-3) = res(nknots*(nknots-1)-2,nknots*(nknots-2)-1) = res(nknots*(nknots-1)-2,nknots*nknots-3) = res(nknots*(nknots-1)-2,nknots*nknots-1) = 2;
res(nknots*(nknots-1)-2,nknots*(nknots-3)-2) = res(nknots*(nknots-1)-2,nknots*(nknots-1)-4) = 1;

// 8 next to corner edges

// upper left
res(1,1) = 10;
res(1,2) = -6;
res(1,nknots+1) = -6;
res(1,0) = -4;
res(1,nknots) = res(1,nknots+2) = 2;
res(1,3) = res(1,2*nknots+1) = 1;

res(nknots,nknots) = 10;
res(nknots,2*nknots) = -6;
res(nknots,nknots+1) = -6;
res(nknots,0) = -4;
res(nknots,1) = res(nknots,2*nknots+1) = 2;
res(nknots,nknots+2) = res(nknots,3*nknots) = 1;

// upper right
res(nknots-2,nknots-2) = 10;
res(nknots-2,nknots-3) = -6;
res(nknots-2,2*nknots-2) = -6;
res(nknots-2,nknots-1) = -4;
res(nknots-2,2*nknots-3) = res(nknots-2,2*nknots-1) = 2;
res(nknots-2,nknots-4) = res(nknots-2,3*nknots-2) = 1;

res(2*nknots-1,2*nknots-1) = 10;
res(2*nknots-1,3*nknots-1) = -6;
res(2*nknots-1,2*nknots-2) = -6;
res(2*nknots-1,nknots-1) = -4;
res(2*nknots-1,nknots-2) = res(2*nknots-1,3*nknots-2) = 2;
res(2*nknots-1,2*nknots-3) = res(2*nknots-1,4*nknots-1) = 1;

// lower left
res(nknots*(nknots-2),nknots*(nknots-2)) = 10;
res(nknots*(nknots-2),nknots*(nknots-3)) = -6;
res(nknots*(nknots-2),nknots*(nknots-2)+1) = -6;
res(nknots*(nknots-2),nknots*(nknots-1)) = -4;
res(nknots*(nknots-2),nknots*(nknots-3)+1) = res(nknots*(nknots-2),nknots*(nknots-1)+1) = 2;
res(nknots*(nknots-2),nknots*(nknots-2)+2) = res(nknots*(nknots-2),nknots*(nknots-4)) = 1;

res(nknots*(nknots-1)+1,nknots*(nknots-1)+1) = 10;
res(nknots*(nknots-1)+1,nknots*(nknots-1)+2) = -6;
res(nknots*(nknots-1)+1,nknots*(nknots-2)+1) = -6;
res(nknots*(nknots-1)+1,nknots*(nknots-1)) = -4;
res(nknots*(nknots-1)+1,nknots*(nknots-2)) = res(nknots*(nknots-1)+1,nknots*(nknots-2)+2) = 2;
res(nknots*(nknots-1)+1,nknots*(nknots-1)+3) = res(nknots*(nknots-1)+1,nknots*(nknots-3)+1) = 1;

// lower right
res(nknots*(nknots-1)-1,nknots*(nknots-1)-1) = 10;
res(nknots*(nknots-1)-1,nknots*(nknots-2)-1) = -6;
res(nknots*(nknots-1)-1,nknots*(nknots-1)-2) = -6;
res(nknots*(nknots-1)-1,nknots*nknots-1) = -4;
res(nknots*(nknots-1)-1,nknots*nknots-2) = res(nknots*(nknots-1)-1,nknots*(nknots-2)-2) = 2;
res(nknots*(nknots-1)-1,nknots*(nknots-1)-3) = res(nknots*(nknots-1)-1,nknots*(nknots-3)-1) = 1;

res(nknots*nknots-2,nknots*nknots-2) = 10;
res(nknots*nknots-2,nknots*nknots-3) = -6;
res(nknots*nknots-2,nknots*(nknots-1)-2) = -6;
res(nknots*nknots-2,nknots*nknots-1) = -4;
res(nknots*nknots-2,nknots*(nknots-1)-3) = res(nknots*nknots-2,nknots*(nknots-1)-1) = 2;
res(nknots*nknots-2,nknots*nknots-4) = res(nknots*nknots-2,nknots*(nknots-2)-2) = 1;

for(i=2; i<nknots-2; i++)
  {
  // real edges

  // upper edge
  res(i,i) = 11;
  res(i,i-1) = res(i,i+1) = res(i,nknots+i) = -6;
  res(i,nknots+i-1) = res(i,nknots+i+1) = 2;
  res(i,i-2) = res(i,i+2) = res(i,2*nknots+i) = 1;

  // left edge
  res(i*nknots,i*nknots) = 11;
  res(i*nknots,i*nknots+1) = res(i*nknots,(i-1)*nknots) = res(i*nknots,(i+1)*nknots) = -6;
  res(i*nknots,(i-1)*nknots+1) = res(i*nknots,(i+1)*nknots+1) = 2;
  res(i*nknots,i*nknots+2) = res(i*nknots,(i-2)*nknots) = res(i*nknots,(i+2)*nknots) = 1;

  // right edge
  res((i+1)*nknots-1,(i+1)*nknots-1) = 11;
  res((i+1)*nknots-1,i*nknots-1) = res((i+1)*nknots-1,(i+1)*nknots-2) = res((i+1)*nknots-1,(i+2)*nknots-1) = -6;
  res((i+1)*nknots-1,i*nknots-2) = res((i+1)*nknots-1,(i+2)*nknots-2) = 2;
  res((i+1)*nknots-1,(i+1)*nknots-3) = res((i+1)*nknots-1,(i-1)*nknots-1) = res((i+1)*nknots-1,(i+3)*nknots-1) = 1;

  // lower edge
  res(nknots*(nknots-1)+i,nknots*(nknots-1)+i) = 11;
  res(nknots*(nknots-1)+i,nknots*(nknots-1)+i-1) = res(nknots*(nknots-1)+i,nknots*(nknots-1)+i+1) = res(nknots*(nknots-1)+i,nknots*(nknots-2)+i) = -6;
  res(nknots*(nknots-1)+i,nknots*(nknots-2)+i-1) = res(nknots*(nknots-1)+i,nknots*(nknots-2)+i+1) = 2;
  res(nknots*(nknots-1)+i,nknots*(nknots-1)+i-2) = res(nknots*(nknots-1)+i,nknots*(nknots-1)+i+2) = res(nknots*(nknots-1)+i,nknots*(nknots-3)+i) = 1;

  // second order edges

  // upper edge
  res(nknots+i,nknots+i) = 19;
  res(nknots+i,nknots+i-1) = res(nknots+i,nknots+i+1) = res(nknots+i,2*nknots+i) = -8;
  res(nknots+i,i) = -6;
  res(nknots+i,i-1) = res(nknots+i,i+1) = res(nknots+i,2*nknots+i-1) = res(nknots+i,2*nknots+i+1) = 2;
  res(nknots+i,nknots+i-2) = res(nknots+i,nknots+i+2) = res(nknots+i,3*nknots+i) = 1;

  // left edge
  res(i*nknots+1,i*nknots+1) = 19;
  res(i*nknots+1,i*nknots+2) = res(i*nknots+1,(i-1)*nknots+1) = res(i*nknots+1,(i+1)*nknots+1) = -8;
  res(i*nknots+1,i*nknots) = -6;
  res(i*nknots+1,(i-1)*nknots) = res(i*nknots+1,(i-1)*nknots+2) = res(i*nknots+1,(i+1)*nknots) = res(i*nknots+1,(i+1)*nknots+2) = 2;
  res(i*nknots+1,i*nknots+3) = res(i*nknots+1,(i-2)*nknots+1) = res(i*nknots+1,(i+2)*nknots+1) = 1;

  // right edge
  res((i+1)*nknots-2,(i+1)*nknots-2) = 19;
  res((i+1)*nknots-2,(i+1)*nknots-3) = res((i+1)*nknots-2,i*nknots-2) = res((i+1)*nknots-2,(i+2)*nknots-2) = -8;
  res((i+1)*nknots-2,(i+1)*nknots-1) = -6;
  res((i+1)*nknots-2,i*nknots-3) = res((i+1)*nknots-2,i*nknots-1) = res((i+1)*nknots-2,(i+2)*nknots-3) = res((i+1)*nknots-2,(i+2)*nknots-1) = 2;
  res((i+1)*nknots-2,(i+1)*nknots-4) = res((i+1)*nknots-2,(i-1)*nknots-2) = res((i+1)*nknots-2,(i+3)*nknots-2) = 1;

  // lower edge
  res(nknots*(nknots-2)+i,nknots*(nknots-2)+i) = 19;
  res(nknots*(nknots-2)+i,nknots*(nknots-2)+i-1) = res(nknots*(nknots-2)+i,nknots*(nknots-2)+i+1) = res(nknots*(nknots-2)+i,nknots*(nknots-3)+i) = -8;
  res(nknots*(nknots-2)+i,nknots*(nknots-1)+i) = -6;
  res(nknots*(nknots-2)+i,nknots*(nknots-3)+i-1) = res(nknots*(nknots-2)+i,nknots*(nknots-3)+i+1) = res(nknots*(nknots-2)+i,nknots*(nknots-1)+i-1) = res(nknots*(nknots-2)+i,nknots*(nknots-1)+i+1) = 2;
  res(nknots*(nknots-2)+i,nknots*(nknots-2)+i-2) = res(nknots*(nknots-2)+i,nknots*(nknots-2)+i+2) = res(nknots*(nknots-2)+i,nknots*(nknots-4)+i) = 1;
  }

// interior

for(i=2; i<nknots-2; i++)
  {
  for(j=2; j<nknots-2; j++)
    {
    res(i*nknots+j,i*nknots+j) = 20;
    res(i*nknots+j,i*nknots+j-1) = res(i*nknots+j,i*nknots+j+1) = res(i*nknots+j,(i-1)*nknots+j) = res(i*nknots+j,(i+1)*nknots+j) = -8;
    res(i*nknots+j,(i-1)*nknots+j-1) = res(i*nknots+j,(i-1)*nknots+j+1) = res(i*nknots+j,(i+1)*nknots+j-1) = res(i*nknots+j,(i+1)*nknots+j+1) = 2;
    res(i*nknots+j,i*nknots+j-2) = res(i*nknots+j,i*nknots+j+2) = res(i*nknots+j,(i-2)*nknots+j) = res(i*nknots+j,(i+2)*nknots+j) = 1;
    }
  }

/*
// Corners
  res(0,1) = 1;
  res(0,nknots) = 1;

  res(nknots-1,nknots-2) = 1;
  res(nknots-1,2*nknots-1) = 1;

  res(nknots*(nknots-1),nknots*(nknots-2)) = 1;
  res(nknots*(nknots-1),nknots*(nknots-1)+1) = 1;

  res(nknots*nknots-1,nknots*(nknots-1)-1) = 1;
  res(nknots*nknots-1,nknots*nknots-2) = 1;

// edges

  for(i=1; i<nknots-1; i++)
    {
    res(i,i-1) = 1;
    res(i,i+1) = 1;
    res(i,nknots+i) = 1;

    res(nknots*(nknots-1)+i,nknots*(nknots-1)+i-1) = 1;
    res(nknots*(nknots-1)+i,nknots*(nknots-1)+i+1) = 1;
    res(nknots*(nknots-1)+i,nknots*(nknots-2)+i) = 1;

    res(i*nknots,(i-1)*nknots) = 1;
    res(i*nknots,i*nknots+1) = 1;
    res(i*nknots,(i+1)*nknots) = 1;

    res((i+1)*nknots-1,i*nknots-1) = 1;
    res((i+1)*nknots-1,(i+1)*nknots-2) = 1;
    res((i+1)*nknots-1,(i+2)*nknots-1) = 1;
    }

// interior

  for(i=1; i<nknots-1; i++)
    {
    for(j=1; j<nknots-1; j++)
      {
      res(i*nknots+j ,i*nknots+j-1) = 1;
      res(i*nknots+j ,i*nknots+j+1) = 1;
      res(i*nknots+j ,(i-1)*nknots+j) = 1;
      res(i*nknots+j ,(i+1)*nknots+j) = 1;
      }
    }

// diagonal

  for(i=0; i<nknots*nknots; i++)
    {
    res(i,i) = -res.sum(i);
    }

  res = res.transposed()*res;
  */
  return res;
  }
}

statmatrix<double> diffmat(const int k, const int d)
  {
  assert(k>=1);
  assert(k<=2);
  int i;
  statmatrix<double>res(d-k,d,0);
  for(i=0; i<d-k; i++)
    {
    if(k==1)
      {
      res(i,i)=-1;
      res(i,i+1)=1;
      }
    else if(k==2)
      {
      res(i,i)=1;
      res(i,i+1)=-2;
      res(i,i+2)=1;
      }
    }
  return res;
  }

statmatrix<double> diffmat_k(const int k, const int d)
  {
  assert(k>=0);
  assert(k<d);
  int i,j,m;
  statmatrix<double>res(d-k,d,0);
  for(i=0;i<int(res.rows());i++)
    res(i,i+k) = 1.0;
  for(i=1;i<=k;i++)
    {
    for(j=0;j<int(res.rows());j++)
      for(m=k-i;m<k;m++)
        res(j,j+m) = res(j,j+m)-res(j,j+m+1);
    }
  return res;
  }


statmatrix<double> weighteddiffmat(const int k, const vector<double> & weight)
  {
  assert(k>=1);
  assert(k<=2);
  unsigned i;
  unsigned d = weight.size();
  statmatrix<double>res(d-k,d,0);
  if(k==1)
    {
    for(i=0; i<d-1; i++)
      {
      res(i,i)=1/sqrt(weight[i+1]);
      res(i,i+1)=-res(i,i);
      }
    }
  else
    {
    for(i=0; i<d-2; i++)
      {
      res(i,i+2)=1/sqrt(weight[i+2]*(1+weight[i+2]/weight[i+1]));
      res(i,i+1)=-(1+weight[i+2]/weight[i+1])*res(i,i+2);
      res(i,i)=weight[i+2]/weight[i+1]*res(i,i+2);
      }
    }
  return res;
  }

statmatrix<double> seasonalfactor(const unsigned & per, const unsigned & s)
  {
  unsigned i,j;
  statmatrix<double> res(s,s-per+1,0);
  for(i=0; i<s-per+1; i++)
    {
    for(j=0; j<per; j++)
      {
      res(i+j,i)=1;
      }
    }
  return res;
  }

statmatrix<double> seasonalX(const unsigned & per, const unsigned & s)
  {
  unsigned i,k;
  statmatrix<double> res(s,per-1,0);
  for(i=0; i<s; i++)
    {
    k=(i+1)%per;
    if(k>0)
      {
      res(i,k-1)=1;
      }
    else
      {
      for(k=0; k<per-1; k++)
        {
        res(i,k)=-1;
        }
      }
    }
  return res;
  }

void rotate(statmatrix<double> & a, const double & s, const double & tau,
            const int & i, const int & j, const int & k, const int & l)
  {
  double g,h;

  g=a(i,j);
  h=a(k,l);
  a(i,j)=g-s*(h+g*tau);
  a(k,l)=h+s*(g-h*tau);
  }

void tridiag(statmatrix<double> & a, statmatrix<double> & d,
             statmatrix<double> & e)
  {
  assert(a.rows()==a.cols());
  assert(a.rows()==d.rows());
  assert(d.cols()==1);
  assert(e.rows()==a.rows());
  assert(e.cols()==1);

  int l, k, j, i;
  double scale, hh, h, g, f;

  int n = d.rows();
  for(i=n-1; i>0; i--)
    {
    l=i-1;
    h=scale=0;
    if(l>0)
      {
      for(k=0; k<l+1; k++)
        {
        scale += fabs(a(i,k));
        }
      if(scale==0)
        {
        e(i,0) = a(i,l);
        }
      else
        {
        for(k=0; k<l+1; k++)
          {
          a(i,k) /= scale;
          h += a(i,k)*a(i,k);
          }
        f = a(i,l);
        g = (f >= 0 ? -sqrt(h) : sqrt(h));
        e(i,0) = scale*g;
        h -= f*g;
        a(i,l) = f-g;
        f=0;
        for(j=0; j<l+1; j++)
          {
          a(j,i) = a(i,j)/h;
          g=0;
          for(k=0; k<j+1; k++)
            {
            g += a(j,k)*a(i,k);
            }
          for(k=j+1; k<l+1; k++)
            {
            g += a(k,j)*a(i,k);
            }
          e(j,0) = g/h;
          f += e(j,0)*a(i,j);
          }
        hh = f/(h+h);
        for(j=0; j<l+1; j++)
          {
          f=a(i,j);
          e(j,0)=g=e(j,0)-hh*f;
          for(k=0; k<j+1; k++)
            {
            a(j,k) -= (f*e(k,0) + g*a(i,k));
            }
          }
        }
      }
    else
      {
      e(i,0) = a(i,l);
      }
    d(i,0)=h;
    }
  d(0,0)=0;
  e(0,0)=0;
  for(i=0; i<n; i++)
    {
    l=i;
    if(d(i,0) != 0)
      {
      for(j=0; j<l; j++)
        {
        g=0;
        for(k=0; k<l; k++)
          {
          g += a(i,k)*a(k,j);
          }
        for(k=0; k<l; k++)
          {
          a(k,j) -= g*a(k,i);
          }
        }
      }
    d(i,0) = a(i,i);
    a(i,i) = 1;
    for(j=0; j<l; j++)
      {
      a(j,i) = a(i,j) = 0;
      }
    }
  }

bool eigentridiag(statmatrix<double> & d, statmatrix<double> & e,
                  statmatrix<double> & z)
  {
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;

  int n=d.rows();
  for(i=1; i<n; i++)
    {
    e(i-1,0)=e(i,0);
    }
  e(n-1,0)=0;

  for(l=0; l<n; l++)
    {
    iter=0;
    do
      {
      for(m=l; m<n-1; m++)
        {
        dd=fabs(d(m,0))+fabs(d(m+1,0));
        if(fabs(e(m,0))+dd == dd)
          {
          break;
          }
        }
      if(m != l)
        {
        if(iter++ == 30)
          {
          return false;
          }
        g=(d(l+1,0)-d(l,0))/(2*e(l,0));
        r=pythag(g,1);
        g=d(m,0)-d(l,0)+e(l,0)/(g+SIGN(r,g));
        s=c=1;
        p=0;
        for(i=m-1; i>=l; i--)
          {
          f=s*e(i,0);
          b=c*e(i,0);
          e(i+1,0)=(r=pythag(f,g));
          if(r==0)
            {
            d(i+1,0) -= p;
            e(m,0)=0;
            break;
            }
          s=f/r;
          c=g/r;
          g=d(i+1,0)-p;
          r=(d(i,0)-g)*s+2*c*b;
          d(i+1,0)=g+(p=s*r);
          g=c*r-b;
          for(k=0; k<n; k++)
            {
            f=z(k,i+1);
            z(k,i+1)=s*z(k,i)+c*f;
            z(k,i)=c*z(k,i)-s*f;
            }
          }
        if(r==0 && i>=1)
          {
          continue;
          }
        d(l,0) -= p;
        e(l,0) = g;
        e(m,0) = 0;
        }
      }
    while (m!=l);
    }
  return true;
  }

double pythag(const double & a, const double & b)
  {
  double absa, absb;

  absa=fabs(a);
  absb=fabs(b);
  if(absa > absb)
    {
    return absa*sqrt(1+sqr(absb/absa));
    }
  else
    {
    return (absb == 0 ? 0 : absb*sqrt(1+sqr(absa/absb)));
    }
  }
double sqr(const double & a)
  {
  return a*a;
  }
double SIGN(const double & a, const double & b)
  {
  return b>=0 ? (a>=0 ? a : -a) : (a>=0 ? -a : a);
  }


int eigen(statmatrix<double> & a, statmatrix<double> & values,
           statmatrix<double> & vectors)
  {
  assert(a.cols()==vectors.cols());
  assert(a.rows()==vectors.rows());
  assert(a.rows()==values.rows());
  assert(values.cols()==1);
  assert(a.cols()==a.rows());

  int i,j,ip,iq;
  double tresh,theta,tau,t,sm,s,h,g,c;

  int n=a.cols();

  vectors=statmatrix<double>::diag(n,1);
  statmatrix<double> b = a.diag();
  values = a.diag();
  statmatrix<double>z(n,1,0);

  for(i=1; i<=50; i++)
    {
    sm=0;
    for(ip=0; ip<n-1; ip++)
      {
      for(iq=ip+1; iq<n; iq++)
        {
        sm += fabs(a(ip,iq));
        }
      }

    if(sm==0)
      {
      return i;
      }

    if(i<4)
      {
      tresh=0.2*sm/(n*n);
      }
    else
      {
      tresh=0;
      }

    for(ip=0; ip<n-1; ip++)
      {
      for(iq=ip+1; iq<n; iq++)
        {
        g=100*fabs(a(ip,iq));
        if(i>4 && (fabs(values(ip,0))+g)==fabs(values(ip,0)) &&
                  (fabs(values(iq,0))+g)==fabs(values(iq,0)))
          {
          a(ip,iq)=0;
          }
        else if(fabs(a(ip,iq))>tresh)
          {
          h=values(iq,0)-values(ip,0);
          if((fabs(h)+g)==fabs(h))
            {
            t=a(ip,iq)/h;
            }
          else
            {
            theta=0.5*h/(a(ip,iq));
            t=1/(fabs(theta)+sqrt(1+theta*theta));
            if(theta<0)
              {
              t=-t;
              }
            }
          c=1/sqrt(1+t*t);
          s=t*c;
          tau=s/(1+c);
          h=t*a(ip,iq);
          z(ip,0) -= h;
          z(iq,0) += h;
          values(ip,0) -= h;
          values(iq,0) += h;
          a(ip,iq)=0;
          for(j=0; j<ip; j++)
            {
            rotate(a,s,tau,j,ip,j,iq);
            }
          for(j=ip+1; j<iq; j++)
            {
            rotate(a,s,tau,ip,j,j,iq);
            }
          for(j=iq+1; j<n; j++)
            {
            rotate(a,s,tau,ip,j,iq,j);
            }
          for(j=0; j<n; j++)
            {
            rotate(vectors,s,tau,j,ip,j,iq);
            }
          }
        }
      }
    for(ip=0; ip<n; ip++)
      {
      b(ip,0) += z(ip,0);
      values(ip,0) = b(ip,0);
      z(ip,0) = 0;
      }
    }
  return i;
  }

bool eigen2(datamatrix & a, datamatrix & d)
  {
  datamatrix e(d.rows(),1,0);
  tridiag(a,d,e);
  return eigentridiag(d,e,a);
  }


void eigensort(datamatrix & values, datamatrix & vectors)
  {
  int i,j,k;
  double p;

  int n=values.rows();
  for(i=0; i<n-1; i++)
    {
    p=values(k=i,0);
    for(j=i; j<n; j++)
      {
      if(values(j,0)>=p)
        {
        p=values(k=j,0);
        }
      }
    if(k!=i)
      {
      values(k,0)=values(i,0);
      values(i,0)=p;
      for(j=0; j<n; j++)
        {
        p=vectors(j,i);
        vectors(j,i)=vectors(j,k);
        vectors(j,k)=p;
        }
      }
    }
  }

statmatrix<double> kronecker(const statmatrix<double> & A, const statmatrix<double> & B)
  {
  statmatrix<double> res(A.rows()*B.rows(),A.cols()*B.cols(),0);
  unsigned i,j;
  unsigned brows=B.rows();
  unsigned bcols=B.cols();
  for(i=0; i<A.rows(); i++)
    {
    for(j=0; j<A.cols(); j++)
      {
      res.putBlock(A(i,j)*B,i*brows,j*bcols,(i+1)*brows,(j+1)*bcols);
      }
    }
  return res;
  }

void compare(const datamatrix & refdata, const datamatrix & neudata, double limit,
             unsigned col, const ST::string & colname, vector<ST::string> & out)
  {

  double diff;
  datamatrix help = datamatrix(refdata.rows(),1,0.0);

  help.minus(neudata.getCol(col),refdata.getCol(col));
  diff = help.norm(0)/refdata.norm(col);

// Ausgabe

  if(diff > limit)
    out.push_back("WARNUNG:");

  out.push_back("  '" + colname + "': " + ST::doubletostring(diff,4) );

  if(diff > limit)
    out.push_back("  Toleranzgrenze: " + ST::doubletostring(limit) );

  }


void compare_nonp(const ST::string & ref, const ST::string & neu, double limit,
                  vector<ST::string> & out)
  {

/* pmean   pqu2p5   pqu10   pmed   pqu90   pqu97p5 */
// pmode   ci95lower   ci80lower   std   ci80upper   ci95upper
// pstd
// variance
// sigma2
// linpred   mu   saturated_deviance   leverage
// eta

  datamatrix refdata;
  datamatrix neudata;

//  Dateien lesen

  ST::string header;

  ifstream in(ref.strtochar());
  if(!in.fail())
    {
    ST::getline(in,50000,header);
    header = header.eatallcarriagereturns();
    refdata.prettyScan(in);
    }

  ifstream in2(neu.strtochar());
  if(!in2.fail())
    {
    ST::getline(in2,50000,header);
    header = header.eatallcarriagereturns();
    neudata.prettyScan(in2);
    }

  list<ST::string> colnames = header.strtokenlist(" \t",false);

// diff berechnen

  out.push_back("Relative quadratische Abweichung zur Referenz in der Datei");
  out.push_back("  '" + neu + "':");
  out.push_back("\n");

  unsigned i = 0;
  list<ST::string>::iterator it = colnames.begin();
  while(it != colnames.end())
    {
    if(*it == "pmean")
      compare(refdata,neudata,limit,i,*it,out);
    i++;
    it++;
    }

//  compare(refdata,neudata,limit,2,"pmean",out);
  compare(refdata,neudata,limit,5,"pmed",out);

  compare(refdata,neudata,2.0*limit,4,"pqu10",out);
  compare(refdata,neudata,2.0*limit,6,"pqu90",out);

  compare(refdata,neudata,2.5*limit,3,"pqu2p5",out);
  compare(refdata,neudata,2.5*limit,7,"pqu97p5",out);

  out.push_back("\n");
  }




