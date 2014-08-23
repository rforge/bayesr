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



#ifdef BORLAND_OUTPUT_WINDOW
#include <vcl.h>
#endif
#pragma hdrstop

#include "bsplinemat.h"

namespace MCMC
{

//---------------------------------------------------------------------------
//----------------------- class: bsplinemat -------------------------------
//---------------------------------------------------------------------------

bsplinemat::bsplinemat(const datamatrix & data, const unsigned & nrk, const unsigned & degr,
                const MCMC::knotpos & kp, const bool & minnull, const deque<double> & k)
  {
  nrknots = nrk;
  degree = degr;
  knpos = kp;
  nrpar = nrknots-1+degree;
  knot = k;

  make_index(data);
  make_Bspline(false,data,minnull);
  }

bsplinemat::bsplinemat(const bool & deriv, const datamatrix & data, const unsigned & nrk, const unsigned & degr,
                const MCMC::knotpos & kp, const bool & minnull, const deque<double> & k)
  {
  nrknots = nrk;
  degree = degr;
  knpos = kp;
  nrpar = nrknots-1+degree;
  knot = k;

  make_index(data);
  make_Bspline(deriv,data,minnull);
  }

bsplinemat::bsplinemat(const bsplinemat & bmat)
  {
  B = bmat.B;
  BS = bmat.BS;

  nrpar = bmat.nrpar;
  nrknots = bmat.nrknots;
  degree = bmat.degree;
  nrdiffobs = bmat.nrdiffobs;

  knpos = bmat.knpos;

  freq = bmat.freq;
  freqoutput = bmat.freqoutput;
  index2 = bmat.index2;
  begcol = bmat.begcol;

  firstnonzero = bmat.firstnonzero;
  lastnonzero = bmat.lastnonzero;
  knot = bmat.knot;

  index = bmat.index;

  Bcolmean = bmat.Bcolmean;
  }


const bsplinemat & bsplinemat::operator=(const bsplinemat & bmat)
  {
  if(&bmat==this)
    return *this;

  B = bmat.B;
  BS = bmat.BS;

  nrpar = bmat.nrpar;
  nrknots = bmat.nrknots;
  degree = bmat.degree;
  nrdiffobs = bmat.nrdiffobs;

  knpos = bmat.knpos;

  freq = bmat.freq;
  freqoutput = bmat.freqoutput;
  index2 = bmat.index2;
  begcol = bmat.begcol;

  firstnonzero = bmat.firstnonzero;
  lastnonzero = bmat.lastnonzero;
  knot = bmat.knot;

  index = bmat.index;

  Bcolmean = bmat.Bcolmean;

  return *this;
  }


void bsplinemat::make_index(const datamatrix & moddata)
  {

  index = statmatrix<int>(moddata.rows(),1);
  index.indexinit();
  moddata.indexsort(index,0,moddata.rows()-1,0,0);
  unsigned i,j;

  int *workindex = index.getV();
  freq.reserve(moddata.rows());

  workindex++;
  freq.push_back(0);
  i = 0;
  for(j=1;j<moddata.rows();j++,workindex++)
    {
    if ( moddata(*workindex,0) != moddata(*(workindex-1),0))
      {
      i++;
      }
    freq.push_back(i);
    }

  freqoutput = freq;
  nrdiffobs = i+1;

  index2.push_back(index(0,0));
  for(i=1;i<moddata.rows();i++)
    index2.push_back(index(i,0)-index(i-1,0));

  }


void bsplinemat::make_Bspline(const bool & deriv, const datamatrix & md, const bool & minnull)
  {

  unsigned i=0,j=0,k=0;
  double value=0.0;
  double * work=NULL;

  vector<int>::iterator freqwork;
  datamatrix help = datamatrix(nrpar,1,0.0);
  Bcolmean = datamatrix(nrpar,1,0.0);

  double min = md(index(0,0),0);
  double max = md(index(md.rows()-1,0),0);
  double dist = max-min;

  min -= 0.01*dist;
  max += 0.01*dist;

  if(minnull)
    min = 0.0;

  if(knot.size()==0)
    {
    if(knpos == equidistant)
      {
      dist = (max - min)/(nrknots-1);
      knot.push_back(min - degree*dist);
      for(i=1;i<nrknots+2*degree;i++)
        knot.push_back(knot[i-1] + dist);
      }
    else if(knpos == quantiles)
      {
      double distfirst, distlast;

      knot.push_back(min);
      for(i=1;i<nrknots-1;i++)
        knot.push_back(md.quantile((i*100)/double(nrknots-1),0));
      knot.push_back(max);

      distfirst = knot[1] - knot[0];
      distlast = knot[nrknots-1] - knot[nrknots-2];

      for(i=0;i<degree;i++)
        {
        knot.push_front(min - (i+1)*distfirst);
        knot.push_back(max + (i+1)*distlast);
        }
      }
    }

  for(i=0;i<nrpar;i++)
    {
    lastnonzero.push_back(-1);
    firstnonzero.push_back(0);
    }

  B = datamatrix(*(freq.end()-1)+1,degree+1,0.0);
  work = B.getV();

  freqwork = freq.begin();
  for(i=0;i<md.rows();i++,++freqwork)
//  for(freqwork=freq.begin();freqwork<freq.end();++freqwork)
    {
    value = md(index(i,0),0);
//    value = md(index(*freqwork,0),0);
    if(freqwork == freq.begin() || *freqwork != *(freqwork-1))
      {
      j=0;
      while(knot[degree+j+1] <= value)
        j++;
      begcol.push_back(j);

      if(deriv)
        help.assign(bspline_derivative(value));
      else
        help.assign(bspline(value));

      for(k=0;k<degree+1;k++,work++)
        {
        *work = help(k+j,0);
        Bcolmean(k+j,0) += *work;
        }
      }

    for(k=j;k<nrpar;k++)
      lastnonzero[k] += 1;
    for(k=j+degree+1;k<nrpar;k++)
      firstnonzero[k] += 1;

    }

  for(i=0;i<nrpar;i++)
    Bcolmean(i,0) /= double(nrdiffobs);

  }


datamatrix bsplinemat::bspline(const double & x)
  {

  datamatrix b(nrpar,1,0.0);
  datamatrix help(nrpar+degree,1,0.0);
  unsigned j;
  double * bwork;
  double * helpwork;

// Grad 0

  for(j=0;j<nrpar;j++)
    if( knot[j]<=x && x<knot[j+1])
      b(j,0) = 1.0;

  for(unsigned l=1;l<=degree;l++)
    {
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
//      help(j,0) = b(j,0);
      *helpwork = *bwork;
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
      {
//      b(j,0) = (x-knot[j])*help(j,0)/(knot[j+l]-knot[j])
//                  + (knot[j+l+1]-x)*help(j+1,0)/(knot[j+l+1]-knot[j+1]);
      *bwork = (x-knot[j])**helpwork/(knot[j+l]-knot[j])
                  + (knot[j+l+1]-x)**(helpwork+1)/(knot[j+l+1]-knot[j+1]);

      }
    }

  return b;

  }


datamatrix bsplinemat::bspline_derivative(const double & x)
  {

  datamatrix b(nrpar,1,0.0);
  datamatrix help(nrpar+degree,1,0.0);
  unsigned j;
  double * bwork;
  double * helpwork;

// Grad 0

  for(j=0;j<nrpar;j++)
    if( knot[j]<=x && x<knot[j+1])
      b(j,0) = 1.0;

  for(unsigned l=1;l<=degree-1;l++)
    {
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
      *helpwork = *bwork;
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
      {
      *bwork = (x-knot[j])**helpwork/(knot[j+l]-knot[j])
                  + (knot[j+l+1]-x)**(helpwork+1)/(knot[j+l+1]-knot[j+1]);
      }
    }

// Hämmerlin/Hoffmann Seite 263

  bwork = b.getV();
  helpwork = help.getV();
  for(j=0;j<nrpar;j++,++helpwork,++bwork)
    *helpwork = *bwork;
  bwork = b.getV();
  helpwork = help.getV();
  for(j=0;j<nrpar;j++,++helpwork,++bwork)
    *bwork = degree*( *helpwork/(knot[j+degree]-knot[j]) - *(helpwork+1)/(knot[j+degree+1]-knot[j+1]) );

  return b;

  }


void bsplinemat::mult(datamatrix & res, const datamatrix & beta)
  {

  double *workres;
  double *workbeta;
  double *workBS;

  vector<int>::iterator freqwork = freq.begin();

  workBS = B.getV();

  unsigned col = degree+1;
  unsigned j,k;
  int i,stop;

  workres = res.getV();
  for(j=0;j<res.rows()*res.cols();j++,workres++)
    *workres = 0.0;

  i = 0;
  k = 0;
  workres = res.getV();
  while (k<nrpar)
    {
    stop = lastnonzero[k];
    while (i <= stop)
      {
      workbeta = beta.getV();
      for(j=0;j<col;j++,workBS++,workbeta++)
        *workres += *workBS * *(workbeta + k);
      if((freqwork+1)!=freq.end() && *freqwork==*(freqwork+1))
        {
        workBS -= col;
        workbeta -= col;
        }
      i++;
      workres++;
      freqwork++;
      }
    k++;
    }

  }


void bsplinemat::mult_index(datamatrix & res, const datamatrix & beta)
  {

  double *workres;
  double *workbeta;
  double *workBS;

  int *workindex;

  vector<int>::iterator freqwork = freq.begin();

  workBS = B.getV();

  unsigned col = degree+1;
  unsigned j,k;
  int i,stop;

  workres = res.getV();
  for(j=0;j<res.rows()*res.cols();j++,workres++)
    *workres = 0.0;

  i = 0;
  k = 0;
  workres = res.getV();
  workindex = index.getV();
  while (k<nrpar)
    {
    stop = lastnonzero[k];
    while (i <= stop)
      {
      workbeta = beta.getV();
      for(j=0;j<col;j++,workBS++,workbeta++)
        res(*workindex,0) += *workBS * *(workbeta + k);
      if((freqwork+1)!=freq.end() && *freqwork==*(freqwork+1))
        {
        workBS -= col;
        workbeta -= col;
        }
      i++;
      freqwork++;
      workindex++;
      }
    k++;
    }

  }


}   // END: namespace MCMC


//---------------------------------------------------------------------------
#if !defined(__BUILDING_GNU)
#pragma package(smart_init)
#endif

