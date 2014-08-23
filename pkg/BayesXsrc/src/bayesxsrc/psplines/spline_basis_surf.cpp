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



#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include "StatwinFrame.h"
#endif

#include "spline_basis_surf.h"


namespace MCMC
{


void spline_basis_surf::make_index2(void)
  {
  unsigned i;

  index2.push_back(index(0,0));
  for(i=1;i<index.rows();i++)
    index2.push_back(index(i,0)-index(i-1,0));
  }


spline_basis_surf::spline_basis_surf(MCMCoptions * o, DISTRIBUTION * dp,
                FULLCOND_const * fcc,const fieldtype & ft,
                const ST::string & ti, const unsigned & nrk, const unsigned & degr,
                const MCMC::knotpos & kp, const double & l, const int & gsx, const int & gsy,
                const ST::string & fp, const ST::string & pres, const unsigned & c)
  : FULLCOND_nonp_basis(o,dp,ft,ti,fp,pres,c)
  {

  fcconst = fcc;

  lambdaconst = false;

  mapexisting = false;

  centertotal = true;

  nrknots = nrk;
  degree = degr;
  knpos = kp;

  lambda = l;
  sigma2 = 1.0/l;

  nrpar1dim = nrknots-1+degree;
  setbeta(nrpar1dim*nrpar1dim,1,0.0);
  beta_uncentered = datamatrix(nrpar,1,0.0);

  intercept = 0.0;
  spline = datamatrix(likep->get_nrobs(),1,0.0);

  gridsizex = gsx;
  gridsizey = gsy;
  if(gridsizex < 0 || gridsizey < 0)
    gridsize = -1;
  else
    gridsize = gridsizex*gridsizey;

  }

// REML

spline_basis_surf::spline_basis_surf(MCMCoptions * o, const datamatrix & v1, const datamatrix & v2,
                      const unsigned & nrk, const unsigned & degr,
                      const fieldtype & ft, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres, const double & l,
                      const double & sl, const bool & catsp, const unsigned & grx,
                      const unsigned & gry)
  : FULLCOND_nonp_basis(o,ti)
  {
  catspecific = catsp;

  mapexisting = false;
  varcoeff = false;

//------------------------------------------------------------------------------

  plotstyle=noplot;

  lambdaconst = false;
//  lambda_prec = -1.0;

  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;

  nrknots = nrk;
  degree = degr;
  knpos = MCMC::equidistant;

  gridsizex = grx;
  gridsizey = gry;
  gridsize = grx*gry;

  transformnonlinear = false;
  transformed =  false;

  type = ft;

  nrpar1dim = nrknots-1+degree;
  nrpar = nrpar1dim*nrpar1dim;

//------------------------------------------------------------------------------

//  samplepath = pres;
  samplepath=fp;

  likep = NULL;
  if(type == mrflinear)
    {
    dimX = 0;
    }
  else if(type == mrfquadratic8)
    {
    dimX = 3;
    }
  else if(type == mrfquadratic12)
    {
    dimX=2;
    }
  dimZ = nrpar-dimX-1;

//------------------------------------------------------------------------------
  if(gridsize>0)
    {
    X_grid = datamatrix(gridsize,dimX,1.0);
    Z_grid = datamatrix(gridsize,dimZ,0.0);
    }

  spline = datamatrix(v1.rows(),1,0);

  lambda = l;
  startlambda = sl;

  make_index(v1,v2);
  compute_knots(v1,v2);
  make_B(v1,v2);

  index2.push_back(index(0,0));
  for(unsigned i=1;i<v1.rows();i++)
    index2.push_back(index(i,0)-index(i-1,0));

  make_xy_v(v1,v2);

  if(gridsize>0)
    {
    make_xy_values_REML(v1,v2);
    make_DG_REML();
    }

  }

// REML (varcoeff)

spline_basis_surf::spline_basis_surf(MCMCoptions * o, const datamatrix & intact,
                      const datamatrix & v1, const datamatrix & v2,
                      const unsigned & nrk, const unsigned & degr,
                      const fieldtype & ft, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres,
                      const double & l, const double & sl, const bool & catsp,
                      const bool & ctr)
  : FULLCOND_nonp_basis(o,ti)
  {
  catspecific = catsp;
  centervcm = ctr;

  mapexisting = false;
  varcoeff=true;

  plotstyle=noplot;

  lambdaconst = false;

  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;

  nrknots = nrk;
  degree = degr;
  knpos = MCMC::equidistant;

  gridsizex = -1;
  gridsizey = -1;
  gridsize = -1;

  transformnonlinear = false;
  transformed =  false;

  type = ft;

  nrpar1dim = nrknots-1+degree;
  nrpar = nrpar1dim*nrpar1dim;

  samplepath=fp;

  likep = NULL;
  dimX = 1;
  dimZ = nrpar-1;

  if(centervcm)
    {
    dimX = dimX-1;
    }

  X_VCM = datamatrix(intact.rows(),dimX,1.0);
  Z_VCM = datamatrix(intact.rows(),dimZ,0.0);

  data_forfixed = intact;

  spline = datamatrix(v1.rows(),1,0);

  lambda = l;
  startlambda = sl;

  make_index(v1,v2);
  compute_knots(v1,v2);
  make_B(v1,v2);

  index2.push_back(index(0,0));
  for(unsigned i=1;i<v1.rows();i++)
    index2.push_back(index(i,0)-index(i-1,0));

  make_xy_v(v1,v2);
  }

// REML geospline

spline_basis_surf::spline_basis_surf(MCMCoptions * o, const datamatrix & region,
                      const MAP::map & mp, const ST::string & mn,
                      const unsigned & nrk, const unsigned & degr,
                      const fieldtype & ft, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres, const double & l,
                      const double & sl, const bool & catsp, const unsigned & grx,
                      const unsigned & gry)
  : FULLCOND_nonp_basis(o,ti)
  {
  catspecific = catsp;

  m = mp;
  mapexisting = true;
  mapname = mn;

  gridsizex = grx;
  gridsizey = gry;
  gridsize = grx*gry;

  if(mp.polygones_existing() == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

  datamatrix v1 = datamatrix(region.rows(),1,0.0);
  datamatrix v2 = datamatrix(region.rows(),1,0.0);

  ST::string regname;
  for(unsigned i=0;i<region.rows();i++)
    {
    regname = ST::doubletostring(region(i,0));
    regionnames.push_back(regname);
    v1(i,0) = m.get_region(m.getnr(regname)).get_xcenter();
    v2(i,0) = m.get_region(m.getnr(regname)).get_ycenter();
    }

//------------------------------------------------------------------------------

  lambdaconst = false;
//  lambda_prec = -1.0;

  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;

  nrknots = nrk;
  degree = degr;
  knpos = MCMC::equidistant;

  varcoeff = false;

  transformnonlinear = false;
  transformed =  false;

  type = ft;

  nrpar1dim = nrknots-1+degree;
  nrpar = nrpar1dim*nrpar1dim;

//------------------------------------------------------------------------------

//  samplepath = pres;
  samplepath=fp;

  likep = NULL;
  if(type == mrflinear)
    {
    dimX = 0;
    }
  else if(type == mrfquadratic8)
    {
    dimX = 3;
    }
  else if(type == mrfquadratic12)
    {
    dimX=2;
    }
  dimZ = nrpar-dimX-1;

//------------------------------------------------------------------------------

  if(gridsize>0)
    {
    X_grid = datamatrix(gridsize,dimX,1.0);
    Z_grid = datamatrix(gridsize,dimZ,0.0);
    }

  spline = datamatrix(v1.rows(),1,0);

  lambda = l;
  startlambda = sl;

  make_index(v1,v2);
  compute_knots(v1,v2);
  make_B(v1,v2);

  index2.push_back(index(0,0));
  for(unsigned i=1;i<v1.rows();i++)
    index2.push_back(index(i,0)-index(i-1,0));

  make_xy_v(v1,v2);

  if(gridsize>0)
    {
    make_xy_values_REML(v1,v2);
    make_DG_REML();
    }

  }

// REML geospline varcoeff

spline_basis_surf::spline_basis_surf(MCMCoptions * o, const datamatrix & intact,
               const datamatrix & region, const MAP::map & mp,
               const ST::string & mn,
               const unsigned & nrk, const unsigned & degr,
               const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp, const bool & ctr)
  : FULLCOND_nonp_basis(o,ti)
  {
  catspecific = catsp;
  centervcm = ctr;

  m = mp;
  mapexisting = true;
  mapname = mn;
  if(mp.polygones_existing() == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

  datamatrix v1 = datamatrix(region.rows(),1,0.0);
  datamatrix v2 = datamatrix(region.rows(),1,0.0);
  data_forfixed = intact;

  ST::string regname;
  for(unsigned i=0;i<region.rows();i++)
    {
    regname = ST::doubletostring(region(i,0));
    regionnames.push_back(regname);
    v1(i,0) = m.get_region(m.getnr(regname)).get_xcenter();
    v2(i,0) = m.get_region(m.getnr(regname)).get_ycenter();
    }

  lambdaconst = false;

  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;

  nrknots = nrk;
  degree = degr;
  knpos = MCMC::equidistant;

  gridsizex = -1;
  gridsizey = -1;
  gridsize = -1;

  varcoeff = true;

  transformnonlinear = false;
  transformed =  false;

  type = ft;

  nrpar1dim = nrknots-1+degree;
  nrpar = nrpar1dim*nrpar1dim;

  samplepath=fp;

  likep = NULL;
  dimX = 1;
  dimZ = nrpar-1;

  if(centervcm)
    {
    dimX = dimX-1;
    }

  X_VCM = datamatrix(region.rows(),dimX,1.0);
  Z_VCM = datamatrix(region.rows(),dimZ,0.0);

  spline = datamatrix(v1.rows(),1,0);

  lambda = l;
  startlambda = sl;

  make_index(v1,v2);
  compute_knots(v1,v2);
  make_B(v1,v2);

  index2.push_back(index(0,0));
  for(unsigned i=1;i<v1.rows();i++)
    index2.push_back(index(i,0)-index(i-1,0));

  make_xy_v(v1,v2);
  }

  
spline_basis_surf::spline_basis_surf(const spline_basis_surf & sp)
  : FULLCOND_nonp_basis(FULLCOND_nonp_basis(sp))
  {

  interactvar = sp.interactvar;  

  fcconst = sp.fcconst;

  lambdaconst = sp.lambdaconst;

  index2 = sp.index2;

  m = sp.m;
  mapexisting = sp.mapexisting;
  mapname = sp.mapname;
  regionnames = sp.regionnames;

  fchelp = sp.fchelp;
  fctotal = sp.fctotal;

  nrdiffobs = sp.nrdiffobs;
  nrpar1dim = sp.nrpar1dim;

  fctotalrespath = sp.fctotalrespath;
  outfile = sp.outfile;

  centertotal = sp.centertotal;

  gridsize = sp.gridsize;
  gridsizex = sp.gridsizex;
  gridsizey = sp.gridsizey;

  effectvaluesx = sp.effectvaluesx;
  effectvaluesy = sp.effectvaluesy;
  xv = sp.xv;
  yv = sp.yv;
  xvalues = sp.xvalues;
  yvalues = sp.yvalues;

  beta1 = sp.beta1;
  beta2 = sp.beta2;
  beta_uncentered = sp.beta_uncentered;
  he1 = sp.he1;
  he2 = sp.he2;

  mainp1 = sp.mainp1;
  mainp2 = sp.mainp2;

  betaweightx = sp.betaweightx;
  betaweighty = sp.betaweighty;
  betaweight = sp.betaweight;
  betaweight_main = sp.betaweight_main;

  spline = sp.spline;
  intercept = sp.intercept;
  splinehelp = sp.splinehelp;

  nrknots = sp.nrknots;
  degree = sp.degree;
  knpos = sp.knpos;

  freq = sp.freq;
  freqoutput = sp.freqoutput;
  knot1 = sp.knot1;
  knot2 = sp.knot2;

  B = sp.B;
  Bout = sp.Bout;
  first = sp.first;

  DG = sp.DG;
  DGfirst = sp.DGfirst;

  X_VCM=sp.X_VCM;
  Z_VCM=sp.Z_VCM;

  X_grid=sp.X_grid;
  Z_grid=sp.Z_grid;

  effectvaluesxgrid = sp.effectvaluesxgrid;
  effectvaluesygrid = sp.effectvaluesygrid;
  xvaluesgrid = sp.xvaluesgrid;
  yvaluesgrid = sp.yvaluesgrid;
  }


const spline_basis_surf & spline_basis_surf::operator=(const spline_basis_surf & sp)
  {
  if(&sp==this)
    return *this;
  FULLCOND_nonp_basis::operator=(FULLCOND_nonp_basis(sp));

  interactvar = sp.interactvar;

  fcconst = sp.fcconst;

  lambdaconst = sp.lambdaconst;

  index2 = sp.index2;

  m = sp.m;
  mapexisting = sp.mapexisting;
  mapname = sp.mapname;
  regionnames = sp.regionnames;

  fchelp = sp.fchelp;
  fctotal = sp.fctotal;

  nrdiffobs = sp.nrdiffobs;
  nrpar1dim = sp.nrpar1dim;

  fctotalrespath = sp.fctotalrespath;
  outfile = sp.outfile;

  centertotal = sp.centertotal;

  gridsize = sp.gridsize;
  gridsizex = sp.gridsizex;
  gridsizey = sp.gridsizey;

  effectvaluesx = sp.effectvaluesx;
  effectvaluesy = sp.effectvaluesy;
  xv = sp.xv;
  yv = sp.yv;
  xvalues = sp.xvalues;
  yvalues = sp.yvalues;

  beta1 = sp.beta1;
  beta2 = sp.beta2;
  beta_uncentered = sp.beta_uncentered;  
  he1 = sp.he1;
  he2 = sp.he2;

  mainp1 = sp.mainp1;
  mainp2 = sp.mainp2;

  betaweightx = sp.betaweightx;
  betaweighty = sp.betaweighty;
  betaweight = sp.betaweight;
  betaweight_main = sp.betaweight_main;

  spline = sp.spline;
  intercept = sp.intercept;
  splinehelp = sp.splinehelp;

  nrknots = sp.nrknots;
  degree = sp.degree;
  knpos = sp.knpos;

  freq = sp.freq;
  freqoutput = sp.freqoutput;
  knot1 = sp.knot1;
  knot2 = sp.knot2;

  B = sp.B;
  Bout = sp.Bout;
  first = sp.first;

  DG = sp.DG;
  DGfirst = sp.DGfirst;

  X_VCM=sp.X_VCM;
  Z_VCM=sp.Z_VCM;

  X_grid=sp.X_grid;
  Z_grid=sp.Z_grid;

  effectvaluesxgrid = sp.effectvaluesxgrid;
  effectvaluesygrid = sp.effectvaluesygrid;
  xvaluesgrid = sp.xvaluesgrid;
  yvaluesgrid = sp.yvaluesgrid;

  return *this;
  }


void spline_basis_surf::make_index(const datamatrix & var1,const datamatrix & var2)
  {

  unsigned i,j,k,beg,end;
  unsigned nrobs = var1.rows();

  index = statmatrix<int>(nrobs,1);
  index.indexinit();
  var1.indexsort(index,0,nrobs-1,0,0);

  i = 0;
  j = 1;
//  freq.reserve(likep->get_nrobs());
  freq.push_back(0);
  while(j<nrobs)
    {
    while(j<nrobs && var1(index(j,0),0)!=var1(index(j-1,0),0))
      {
      i++;
      freq.push_back(i);
      j++;
      }
    beg = j-1;
    while(j<nrobs && var1(index(j,0),0) == var1(index(j-1,0),0))
      j++;
    end = j-1;
    if(end!=beg)
      {
      var2.indexsort(index,beg,end,0,0);
      for(k=beg+1;k<=end;k++)
        {
        if(var2(index(k,0),0) != var2(index(k-1,0),0))
          i++;
        freq.push_back(i);
        }
      }
    }

  freqoutput = freq;
  nrdiffobs = i+1;

  if(varcoeff)
    {
    for(i=0;i<freq.size();i++)
      freq[i] = i;
    nrdiffobs = freq.size();
//    freqoutput = freq;
    }

  }


void spline_basis_surf::compute_knots(const datamatrix & x,const datamatrix & y)
  {

  unsigned i;
  double min,max,dist;

// compute knots
  if(knpos == equidistant)
    {
    min = x.min(0);
    max = x.max(0);
    dist = max - min;
    min -= 0.01*dist;
    max += 0.01*dist;
    dist = (max - min)/(nrknots-1);
    knot1.push_back(min - degree*dist);
    for(i=1;i<nrknots+2*degree;i++)
      knot1.push_back(knot1[i-1] + dist);

    min = y.min(0);
    max = y.max(0);
    dist = max - min;
    min -= 0.01*dist;
    max += 0.01*dist;
    dist = (max - min)/(nrknots-1);
    knot2.push_back(min - degree*dist);
    for(i=1;i<nrknots+2*degree;i++)
      knot2.push_back(knot2[i-1] + dist);
    }
  else if(knpos == quantiles)
    {
    double distfirst, distlast;
    min = x.min(0);
    max = x.max(0);
    dist = max - min;
    min -= 0.01*dist;
    max += 0.01*dist;
    knot1.push_back(min);
    for(i=1;i<nrknots-1;i++)
      knot1.push_back(x.quantile((i*100)/double(nrknots-1),0));
    knot1.push_back(max);
    distfirst = knot1[1] - knot1[0];
    distlast = knot1[nrknots-1] - knot1[nrknots-2];
    for(i=0;i<degree;i++)
      {
      knot1.push_front(min - (i+1)*distfirst);
      knot1.push_back(max + (i+1)*distlast);
      }

    min = y.min(0);
    max = y.max(0);
    dist = max - min;
    min -= 0.01*dist;
    max += 0.01*dist;
    knot2.push_back(min);
    for(i=1;i<nrknots-1;i++)
      knot2.push_back(y.quantile((i*100)/double(nrknots-1),0));
    knot2.push_back(max);

    distfirst = knot2[1] - knot2[0];
    distlast = knot2[nrknots-1] - knot2[nrknots-2];

    for(i=0;i<degree;i++)
      {
      knot2.push_front(min - (i+1)*distfirst);
      knot2.push_back(max + (i+1)*distlast);
      }
    }
// ENDE: compute knots

  }


void spline_basis_surf::make_DG(void)
  {

  int i,j;
  unsigned k,l,m,n;
  datamatrix betahelp(nrpar,1,0);

  DG = datamatrix(gridsize,(degree+1)*(degree+1),0);
  DGfirst = vector<int>(gridsize);

  effectvaluesx = vector<double>(gridsize);
  effectvaluesy = vector<double>(gridsize);

  for(i=0;i<gridsizex;i++)
    {
    for(j=0;j<gridsizey;j++)
      {

      betahelp.assign(bspline(xvalues(i,0),yvalues(j,0)));

      k=degree+1;
      while(knot2[k] <= yvalues(j,0) && k<nrknots+degree)
        k++;
      l=degree+1;
      while(knot1[l] <= xvalues(i,0) && l<nrknots+degree)
        l++;

      for(m=0;m<degree+1;m++)
        for(n=0;n<degree+1;n++)
          DG(i*gridsizey + j,n + m*(degree+1)) = betahelp(n+l-(degree+1) + (m+k-(degree+1))*nrpar1dim,0);

      DGfirst[i*gridsizey + j] = l-(degree+1) + (k-(degree+1))*nrpar1dim;
      effectvaluesx[i*gridsizey + j] = xvalues(i,0);
      effectvaluesy[i*gridsizey + j] = yvalues(j,0);

      }
    }

  }

void spline_basis_surf::make_DG_REML(void)
  {

  int i,j;
  unsigned k,l,m,n;
  datamatrix betahelp(nrpar,1,0);

  DG = datamatrix(gridsize,(degree+1)*(degree+1),0);
  DGfirst = vector<int>(gridsize);

  effectvaluesxgrid = vector<double>(gridsize);
  effectvaluesygrid = vector<double>(gridsize);

  for(i=0;i<gridsizex;i++)
    {
    for(j=0;j<gridsizey;j++)
      {

      betahelp.assign(bspline(xvaluesgrid(i,0),yvaluesgrid(j,0)));

      k=degree+1;
      while(knot2[k] <= yvaluesgrid(j,0) && k<nrknots+degree)
        k++;
      l=degree+1;
      while(knot1[l] <= xvaluesgrid(i,0) && l<nrknots+degree)
        l++;

      for(m=0;m<degree+1;m++)
        for(n=0;n<degree+1;n++)
          DG(i*gridsizey + j,n + m*(degree+1)) = betahelp(n+l-(degree+1) + (m+k-(degree+1))*nrpar1dim,0);

      DGfirst[i*gridsizey + j] = l-(degree+1) + (k-(degree+1))*nrpar1dim;
      effectvaluesxgrid[i*gridsizey + j] = xvaluesgrid(i,0);
      effectvaluesygrid[i*gridsizey + j] = yvaluesgrid(j,0);

      }
    }

  }

void spline_basis_surf::make_xy_v(datamatrix var1,datamatrix var2)
  {

  unsigned i;

  var1.sort(0,var1.rows()-1,0);
  var2.sort(0,var2.rows()-1,0);

  xv.push_back(var1(0,0));
  yv.push_back(var2(0,0));

  for(i=1;i<var1.rows();i++)
    {
    if(var1(i,0)!=var1(i-1,0))
      xv.push_back(var1(i,0));
    if(var2(i,0)!=var2(i-1,0))
      yv.push_back(var2(i,0));
    }

  }


void spline_basis_surf::make_xy_values(const datamatrix & var1,const datamatrix & var2)
  {

  int i;

// Design-Daten erzeugen
  double xmin = var1.min(0);
  double xmax = var1.max(0);
  double ymin = var2.min(0);
  double ymax = var2.max(0);
  xvalues = datamatrix(gridsizex,1);
  yvalues = datamatrix(gridsizey,1);

  for(i=0;i<gridsizex;i++)
    xvalues(i,0) = xmin + i*(xmax-xmin)/double(xvalues.rows()-1);
  for(i=0;i<gridsizey;i++)
    yvalues(i,0) = ymin + i*(ymax-ymin)/double(yvalues.rows()-1);

  }

void spline_basis_surf::make_xy_values_REML(const datamatrix & var1,const datamatrix & var2)
  {

  int i;

// Design-Daten erzeugen
  double xmin = var1.min(0);
  double xmax = var1.max(0);
  double ymin = var2.min(0);
  double ymax = var2.max(0);
  xvaluesgrid = datamatrix(gridsizex,1);
  yvaluesgrid = datamatrix(gridsizey,1);

  for(i=0;i<gridsizex;i++)
    xvaluesgrid(i,0) = xmin + i*(xmax-xmin)/double(xvaluesgrid.rows()-1);
  for(i=0;i<gridsizey;i++)
    yvaluesgrid(i,0) = ymin + i*(ymax-ymin)/double(yvaluesgrid.rows()-1);

  }

void spline_basis_surf::make_B(const datamatrix & x,const datamatrix & y)
  {
/*
  Xcol = vector<vector<double>>(nrpar);
  datamatrix helpXcolbeg = datamatrix(nrpar,nrpar1dim,likep->get_nrobs());
  datamatrix helpXcolend = datamatrix(nrpar,nrpar1dim,0);
  Xcolbeg = vector<vector<int>>(nrpar);
  Xcolend = vector<vector<int>>(nrpar);
  Xcolbeg_ind = vector<vector<int>>(nrpar);
*/
  unsigned i,k,l,m,n;
  unsigned position;
//  unsigned position2;
  vector<int>::iterator freqwork;
  vector<int>::iterator freqoutputwork;

  int * workindex;
  int * workhelpindex;
  vector<int>::iterator workfirst;
  vector<int>::iterator helpfreqwork;

  datamatrix helpbspline(nrpar,1,0);
  statmatrix<int> helpindex(x.rows(),1,0);
  vector<int> helpfreq(x.rows());

  B = datamatrix(nrdiffobs,(degree+1)*(degree+1),0);
  first = vector<int>(x.rows());

  workhelpindex = helpindex.getV();
  workfirst = first.begin();
  helpfreqwork = helpfreq.begin();

  vector<ST::string> regionnameshelp = regionnames;
  regionnames = vector<ST::string>(0);

  position = 0;
//  position2 = 0;
  for(k=degree+1;k<nrknots+degree;k++)
    {
    for(l=degree+1;l<nrknots+degree;l++)
      {
      freqwork = freq.begin();
      freqoutputwork = freqoutput.begin();
      workindex = index.getV();
      for(i=0;i<x.rows();i++,freqwork++,workindex++,freqoutputwork++)
        {

        if(knot1[l-1] <= x(*workindex,0) && x(*workindex,0) < knot1[l]
        && knot2[k-1] <= y(*workindex,0) && y(*workindex,0) < knot2[k])
          {
/*
          for(m=0;m<degree+1;m++)
            for(n=0;n<degree+1;n++)
              {
              if(position2 < helpXcolbeg(n+l-(degree+1) + (m+k-(degree+1))*nrpar1dim,k-(degree+1)))
                helpXcolbeg(n+l-(degree+1) + (m+k-(degree+1))*nrpar1dim,k-(degree+1)) = position2;
              if(position2 > helpXcolend(n+l-(degree+1) + (m+k-(degree+1))*nrpar1dim,k-(degree+1)))
                helpXcolend(n+l-(degree+1) + (m+k-(degree+1))*nrpar1dim,k-(degree+1)) = position2;
              }
*/
          if(freqwork == freq.begin() || *freqwork != *(freqwork-1))
            {
            effectvaluesx.push_back(x(*workindex,0));
            effectvaluesy.push_back(y(*workindex,0));
            if(mapexisting)
                {
                regionnames.push_back(regionnameshelp[*workindex]);
                }
            helpbspline.assign(bspline(x(*workindex,0),y(*workindex,0)));
            for(m=0;m<degree+1;m++)
              for(n=0;n<degree+1;n++)
                {
                B(position,n + m*(degree+1)) = helpbspline(n+l-(degree+1) + (m+k-(degree+1))*nrpar1dim,0);
//                Xcol[n+l-(degree+1) + (m+k-(degree+1))*nrpar1dim].push_back(B(position,n + m*(degree+1)));
                }

            position++;
            }
          *workfirst = l-(degree+1) + (k-(degree+1))*nrpar1dim;
          *helpfreqwork = position-1;
          *workhelpindex = *workindex;
          workfirst++;
          helpfreqwork++;
          workhelpindex++;
//          position2++;
          }
        }
      }
    }
/*
  for(i=0;i<helpXcolbeg.rows();i++)
    {
    Xcolbeg_ind[i].push_back(0);
    for(unsigned j=0;j<helpXcolbeg.cols();j++)
      {
      if(helpXcolbeg(i,j)<=helpXcolend(i,j))
        {
        Xcolbeg[i].push_back(helpXcolbeg(i,j));
        Xcolend[i].push_back(helpXcolend(i,j));
        Xcolbeg_ind[i].push_back(helpXcolend(i,j)-helpXcolbeg(i,j)+1);
        }
      }
    }

  for(i=0;i<helpXcolbeg.rows();i++)
    {
    for(unsigned j=1;j<Xcolbeg_ind[i].size();j++)
      Xcolbeg_ind[i][j] += Xcolbeg_ind[i][j-1];
    }
*/
  freq = helpfreq;
  freqoutput = freq;
  index.assign(helpindex);

  make_index2();

  if(varcoeff)
    Bout = B;

  }


void spline_basis_surf::make_BVC(const datamatrix & intact)
  {

  unsigned i,j;
  vector<int>::iterator freqwork = freq.begin();

  for(i=0;i<likep->get_nrobs();i++,freqwork++)
    {
    if(freqwork == freq.begin() || *freqwork != *(freqwork-1))
      {
      for(j=0;j<B.cols();j++)
        {
        B(*freqwork,j) *= intact(index(i,0),0);
        }
      }
    }

  }


datamatrix spline_basis_surf::bspline(const double x, const double y)
  {

  datamatrix b;

  datamatrix b1(nrpar1dim,1,0.0);
  datamatrix b2(nrpar1dim,1,0.0);

  unsigned i,j;

// x
  b = datamatrix(nrpar1dim+degree,1,0.0);
  for(j=0;j<nrpar1dim;j++)
    {
    if( knot1[j]<=x && x<knot1[j+1])
      b1(j,0) = 1.0;
    }
  for(unsigned l=1;l<=degree;l++)
    {
    for(j=0;j<nrpar1dim;j++)
      b(j,0) = b1(j,0);
    for(j=0;j<nrpar1dim;j++)
      b1(j,0) = (x-knot1[j])*b(j,0)/(knot1[j+l]-knot1[j])
                    + (knot1[j+l+1]-x)*b(j+1,0)/(knot1[j+l+1]-knot1[j+1]);
    }

// y
  b = datamatrix(nrpar1dim+degree,1,0.0);
  for(j=0;j<nrpar1dim;j++)
    {
    if( knot2[j]<=y && y<knot2[j+1])
      b2(j,0) = 1.0;
    }
  for(unsigned l=1;l<=degree;l++)
    {
    for(j=0;j<nrpar1dim;j++)
      b(j,0) = b2(j,0);
    for(j=0;j<nrpar1dim;j++)
      b2(j,0) = (y-knot2[j])*b(j,0)/(knot2[j+l]-knot2[j])
                    + (knot2[j+l+1]-y)*b(j+1,0)/(knot2[j+l+1]-knot2[j+1]);
    }

// xy
  b = datamatrix(nrpar,1,0.0);
  for(i=0;i<nrpar1dim;i++)
    for(j=0;j<nrpar1dim;j++)
      b(i+j*nrpar1dim,0) = b1(i,0) * b2(j,0);

  return b;

  }


datamatrix spline_basis_surf::bspline(const double x, const deque<double> knot)
  {

  datamatrix b(nrpar1dim,1,0.0);
  datamatrix help(nrpar1dim+degree,1,0.0);
  unsigned j;
  double * bwork;
  double * helpwork;

// Grad 0

  for(j=0;j<nrpar1dim;j++)
    if( knot[j]<=x && x<knot[j+1])
      b(j,0) = 1.0;

  for(unsigned l=1;l<=degree;l++)
    {
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar1dim;j++,++helpwork,++bwork)
      *helpwork = *bwork;
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar1dim;j++,++helpwork,++bwork)
      {
      *bwork = (x-knot[j])**helpwork/(knot[j+l]-knot[j])
                  + (knot[j+l+1]-x)**(helpwork+1)/(knot[j+l+1]-knot[j+1]);
      }
    }

  return b;

  }


double spline_basis_surf::bspline_rek(unsigned l, unsigned nu, const double value, const bool xvar)
  {

  if(xvar)
    {
    if(l==0)
      {
      if(knot1[nu] <= value && value < knot1[nu+1])
        return 1.0;
      else
        return 0.0;
      }
    else
      {
      return (value-knot1[nu])*bspline_rek(l-1,nu,value,xvar)/(knot1[nu+l]-knot1[nu])
        +(knot1[nu+l+1]-value)*bspline_rek(l-1,nu+1,value,xvar)/(knot1[nu+l+1]-knot1[nu+1]);
      }
    }
  else
    {
    if(l==0)
      {
      if(knot2[nu] <= value && value < knot2[nu+1])
        return 1.0;
      else
        return 0.0;
      }
    else
      {
      return (value-knot2[nu])*bspline_rek(l-1,nu,value,xvar)/(knot2[nu+l]-knot2[nu])
        +(knot2[nu+l+1]-value)*bspline_rek(l-1,nu+1,value,xvar)/(knot2[nu+l+1]-knot2[nu+1]);
      }
    }

  }


void spline_basis_surf::outresults(void)
  {

  optionsp->out("  Results are stored in file\n");
  optionsp->out("  " + pathcurrent + "\n");
  optionsp->out("\n");

  if (transformnonlinear)
    {
    if (transformtype=="elasticity")
      {
      optionsp->out("  Results for elasticities could not be computed because of\n");
      optionsp->out("  missing derivatives\n");
      optionsp->out("\n");
      }
    else
      {
      ST::string suffix = "";
      fchelp.set_transform(suffix,transformtype);
      vector<FULLCOND*> fcvec(1);
      fcvec[0] = &fchelp;
      likep->transform_nonlinear(fcvec,transformtype);
      }
    }


if(mapexisting)
  {
  #if defined(JAVA_OUTPUT_WINDOW)
  optionsp->out("  Postscript files are stored in files\n");
  ST::string psfile;
  psfile = pathcurrent.substr(0,pathcurrent.length()-4) + ".ps";
  optionsp->out("  " + psfile + "\n");
  psfile = pathcurrent.substr(0,pathcurrent.length()-4) + "_pcatbig" + ".ps";
  optionsp->out("  " + psfile + "\n");
  psfile = pathcurrent.substr(0,pathcurrent.length()-4) + "_pcatsmall" + ".ps";
  optionsp->out("  " + psfile + "\n");
  optionsp->out("\n");
  optionsp->out("  Results may be visualized using method 'drawmap'\n");
  optionsp->out("  Type for example:\n");
  optionsp->out("  objectname.drawmap " + ST::inttostring(fcnumber) + "\n");
  #else
  optionsp->out("  Results may be visualized using the R function 'drawmap' \n");
  #endif
  optionsp->out("\n");
  }
else
  {
  optionsp->out(ST::string("  Results may be visualized using the R function 'plotsurf' \n"));
  ST::string doublebackslash = "/";
  ST::string spluspath = pathcurrent.insert_string_char('\\',doublebackslash);
  optionsp->out("  Type for example:\n");
  optionsp->out("  plotsurf(\"" + spluspath + "\")\n");
  optionsp->out("\n");
  }

  unsigned i;

  fchelp.outresults();

  ofstream outres(pathcurrent.strtochar());

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  outres << "intnr" << "   ";
if(mapexisting)
  outres << datanames[0] << "   " << "x   " << "y   ";
else
  outres << datanames[0] << "   " << datanames[1] << "   ";
  outres << "pmean   ";
  outres << "pqu"  << l1  << "   ";
  outres << "pqu"  << l2  << "   ";
  outres << "pmed   ";
  outres << "pqu"  << u1  << "   ";
  outres << "pqu"  << u2  << "   ";
  outres << "pcat" << level1 << "   ";
  outres << "pcat" << level2 << "   ";

  outres << endl;

  double * workmean = fchelp.get_betameanp();
  double * workbetaqu_l1_lower_p = fchelp.get_beta_lower1_p();
  double * workbetaqu_l2_lower_p = fchelp.get_beta_lower2_p();
  double * workbetaqu_l1_upper_p = fchelp.get_beta_upper1_p();
  double * workbetaqu_l2_upper_p = fchelp.get_beta_upper2_p();
  double * workbetaqu50 = fchelp.get_betaqu50p();

  vector<double>::iterator effitx = effectvaluesx.begin();
  vector<double>::iterator effity = effectvaluesy.begin();

  for(i=0;i<effectvaluesx.size();i++,workmean++,
                     workbetaqu_l1_lower_p++,workbetaqu_l2_lower_p++,
                    workbetaqu50++,
                   workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                  effitx++,effity++)
    {
    outres << (i+1) << "   ";
  if(mapexisting)
    outres << regionnames[i] << "   ";
    outres << *effitx << "   " << *effity << "   ";
    outres << *workmean << "   ";
    outres << *workbetaqu_l1_lower_p << "   ";
    outres << *workbetaqu_l2_lower_p << "   ";
    outres << *workbetaqu50 << "   ";
    outres << *workbetaqu_l2_upper_p << "   ";
    outres << *workbetaqu_l1_upper_p << "   ";

    if (*workbetaqu_l1_lower_p > 0)
      outres << 1 << "   ";
    else if (*workbetaqu_l1_upper_p < 0)
      outres << -1 << "   ";
    else
      outres << 0 << "   ";

    if (*workbetaqu_l2_lower_p > 0)
      outres << 1 << "   ";
    else if (*workbetaqu_l2_upper_p < 0)
      outres << -1 << "   ";
    else
      outres << 0 << "   ";

    outres << endl;
    }

if(!mapexisting)
  {
  if(centertotal==false)
    {

    optionsp->out("  Results for total effect are stored in file\n");
    optionsp->out("  " + fctotalrespath + "\n");
    optionsp->out("\n");

    fctotal.outresults();

    ofstream outrestotal(fctotalrespath.strtochar());

    outrestotal << "intnr" << "   ";
    outrestotal << datanames[0] << "   ";
    outrestotal << datanames[1] << "   ";
    outrestotal << "pmean   ";
    outrestotal << "pqu"  << l1  << "   ";
    outrestotal << "pqu"  << l2  << "   ";
    outrestotal << "pmed   ";
    outrestotal << "pqu"  << u1  << "   ";
    outrestotal << "pqu"  << u2  << "   ";
    outrestotal << "pcat" << level1 << "   ";
    outrestotal << "pcat" << level2 << "   ";
    outrestotal << endl;

    workmean = fctotal.get_betameanp();
    workbetaqu_l1_lower_p = fctotal.get_beta_lower1_p();
    workbetaqu_l2_lower_p = fctotal.get_beta_lower2_p();
    workbetaqu_l1_upper_p = fctotal.get_beta_upper1_p();
    workbetaqu_l2_upper_p = fctotal.get_beta_upper2_p();
    workbetaqu50 = fctotal.get_betaqu50p();
    effitx = effectvaluesx.begin();
    effity = effectvaluesy.begin();


    for(i=0;i<effectvaluesx.size();i++,workmean++,
                     workbetaqu_l1_lower_p++,workbetaqu_l2_lower_p++,
                      workbetaqu50++,
                   workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
    effitx++,effity++)
      {
      outrestotal << (i+1) << "   ";
      outrestotal << *effitx << "   ";
      outrestotal << *effity << "   ";
      outrestotal << *workmean << "   ";
      outrestotal << *workbetaqu_l1_lower_p << "   ";
      outrestotal << *workbetaqu_l2_lower_p << "   ";
      outrestotal << *workbetaqu50 << "   ";
      outrestotal << *workbetaqu_l2_upper_p << "   ";
      outrestotal << *workbetaqu_l1_upper_p << "   ";

      if (*workbetaqu_l1_lower_p > 0)
        outrestotal << 1 << "   ";
      else if (*workbetaqu_l1_upper_p < 0)
        outrestotal << -1 << "   ";
      else
        outrestotal << 0 << "   ";

      if (*workbetaqu_l2_lower_p > 0)
        outrestotal << 1 << "   ";
      else if (*workbetaqu_l2_upper_p < 0)
        outrestotal << -1 << "   ";
      else
        outrestotal << 0 << "   ";

      outrestotal << endl;
      }
    }
  else if(centertotal && optionsp->get_samplesize()>0)
    {
    ST::string of;

    optionsp->out("  Extracted results for main effects are stored in files\n");

    if(likep->get_linearpred().cols()==1 || column==0)
      of = outfile + "_f_" + datanames[0] + "_pspline.res";
    else
      of = outfile + "_f_" + ST::inttostring(column+1) + "_" + datanames[0] + "_pspline.res";

    optionsp->out("  " + of + "\n");
    ofstream xout(of.strtochar());

    xout << "intnr" << "   ";
    xout << datanames[0] << "   ";
    xout << "pmean   ";
    xout << "pqu"  << l1  << "   ";
    xout << "pqu"  << l2  << "   ";
    xout << "pmed   ";
    xout << "pqu"  << u1  << "   ";
    xout << "pqu"  << u2  << "   ";
    xout << "pcat" << level1 << "   ";
    xout << "pcat" << level2 << "   ";
    xout << endl;

    if(likep->get_linearpred().cols()==1 || column==0)
      of = outfile + "_f_" + datanames[1] + "_pspline.res";
    else
      of = outfile + "_f_" + ST::inttostring(column+1) + "_" + datanames[1] + "_pspline.res";

    optionsp->out("  " + of + "\n");
    optionsp->out("\n");
    ofstream yout(of.strtochar());

    yout << "intnr" << "   ";
    yout << datanames[1] << "   ";
    yout << "pmean   ";
    yout << "pqu"  << l1  << "   ";
    yout << "pqu"  << l2  << "   ";
    yout << "pmed   ";
    yout << "pqu"  << u1  << "   ";
    yout << "pqu"  << u2  << "   ";
    yout << "pcat" << level1 << "   ";
    yout << "pcat" << level2 << "   ";
    yout << endl;

    if(likep->get_linearpred().cols()==1 || column==0)
      of = outfile + "_f_" + datanames[0] + "_" + datanames[1] + "_pspline_interact.res";
    else
      of = outfile + "_f_" + ST::inttostring(column+1) + "_" + datanames[0] + "_" + datanames[1] + "_pspline_interact.res";

    optionsp->out("  Extracted results for interaction effect are stored in file\n");
    optionsp->out("  " + of + "\n");
    optionsp->out("\n");
    ofstream interactout(of.strtochar());

    interactout << "intnr" << "   ";
    interactout << datanames[0] << "   ";
    interactout << datanames[1] << "   ";
    interactout << "pmean   ";
    interactout << "pqu"  << l1  << "   ";
    interactout << "pqu"  << l2  << "   ";
    interactout << "pmed   ";
    interactout << "pqu"  << u1  << "   ";
    interactout << "pqu"  << u2  << "   ";
    interactout << "pcat" << level1 << "   ";
    interactout << "pcat" << level2 << "   ";
    interactout << endl;

    unsigned j,k,l;

    double l1i,l2i,u1i,u2i;

    datamatrix interact;
    datamatrix x1mean;
    datamatrix x2mean;

    if(gridsize<0)
      {

      datamatrix sample_k(nrpar,1,0);
      double * x1meanp;
      double * x2meanp;
      double * he1p;
      double * he2p;
      double * interactp;

      x1mean = datamatrix(optionsp->get_samplesize(),xv.size(),0);
      x2mean = datamatrix(optionsp->get_samplesize(),yv.size(),0);
      interact = datamatrix(optionsp->get_samplesize(),nrdiffobs,0);

      he1 = datamatrix(xv.size(),1,0);
      he2 = datamatrix(yv.size(),1,0);

      effitx = effectvaluesx.begin();
      effity = effectvaluesy.begin();

      fchelp.readsample3(interact);

      x1meanp = x1mean.getV();
      x2meanp = x2mean.getV();
      interactp = interact.getV();
      for(k=0;k<optionsp->get_samplesize();k++)
        {

        if(breakpause())
          break;

        readsample2(sample_k,k);

        he1.mult(betaweightx,sample_k);
        he2.mult(betaweighty,sample_k);

        he1p = he1.getV();
        for(i=0;i<he1.rows();i++,x1meanp++,he1p++)
          *x1meanp += *he1p;
        he2p = he2.getV();
        for(i=0;i<he2.rows();i++,x2meanp++,he2p++)
          *x2meanp += *he2p;

        }

      for(i=0;i<nrdiffobs;i++,++effitx,++effity)
        {

        if(breakpause())
          break;

        interactp = interact.getV()+i;
        vector<double>::iterator xvit = xv.begin();
        vector<double>::iterator yvit = yv.begin();
        x1meanp = x1mean.getV();
        x2meanp = x2mean.getV();

        while(*effitx!=*xvit)
          {
          xvit++;
          x1meanp++;
          }
        while(*effity!=*yvit)
          {
          yvit++;
          x2meanp++;
          }
        for(k=0;k<optionsp->get_samplesize();k++,interactp+=interact.cols(),
                                        x1meanp+=x1mean.cols(),x2meanp+=x2mean.cols())
          *interactp -= *x1meanp + *x2meanp;

        }

      for(i=0;i<x1mean.cols();i++)
        {

        l1i = x1mean.quantile(lower1,i);
        l2i = x1mean.quantile(lower2,i);
        u1i = x1mean.quantile(upper1,i);
        u2i = x1mean.quantile(upper2,i);

        xout << (i+1) << "   ";
        xout << xv[i] << "   ";
        xout << x1mean.mean(i) << "   ";
        xout << l1i << "   ";
        xout << l2i << "   ";
        xout << x1mean.quantile(50,i) << "   ";
        xout << u1i << "   ";
        xout << u2i << "   ";

        if (l1i > 0)
          xout << 1 << "   ";
        else if (u2i < 0)
          xout << -1 << "   ";
        else
          xout << 0 << "   ";

        if (l2i > 0)
          xout << 1 << "   ";
        else if (u1i < 0)
          xout << -1 << "   ";
        else
          xout << 0 << "   ";

        xout << endl;
        }

      for(i=0;i<x2mean.cols();i++)
        {

        l1i = x2mean.quantile(lower1,i);
        l2i = x2mean.quantile(lower2,i);
        u1i = x2mean.quantile(upper1,i);
        u2i = x2mean.quantile(upper2,i);

        yout << (i+1) << "   ";
        yout << yv[i] << "   ";
        yout << x2mean.mean(i) << "   ";
        yout << l1i << "   ";
        yout << l2i << "   ";
        yout << x2mean.quantile(50,i) << "   ";
        yout << u1i << "   ";
        yout << u2i << "   ";

        if (l1i > 0)
          yout << 1 << "   ";
        else if (u2i < 0)
          yout << -1 << "   ";
        else
          yout << 0 << "   ";

        if (l2i > 0)
          yout << 1 << "   ";
        else if (u1i < 0)
          yout << -1 << "   ";
        else
          yout << 0 << "   ";

        yout << endl;
        }

      effitx = effectvaluesx.begin();
      effity = effectvaluesy.begin();
      for(i=0;i<interact.cols();i++,effitx++,effity++)
        {

        l1i = interact.quantile(lower1,i);
        l2i = interact.quantile(lower2,i);
        u1i = interact.quantile(upper1,i);
        u2i = interact.quantile(upper2,i);

        interactout << (i+1) << "   ";
        interactout << *effitx << "   ";
        interactout << *effity << "   ";
        interactout << interact.mean(i) << "   ";
        interactout << l1i << "   ";
        interactout << l2i << "   ";
        interactout << interact.quantile(50,i) << "   ";
        interactout << u1i << "   ";
        interactout << u2i << "   ";

        if (l1i > 0)
          interactout << 1 << "   ";
        else if (u2i < 0)
          interactout << -1 << "   ";
        else
          interactout << 0 << "   ";

        if (l2i > 0)
          interactout << 1 << "   ";
        else if (u1i < 0)
          interactout << -1 << "   ";
        else
          interactout << 0 << "   ";

        interactout << endl;
        }
      }
    else
      {

      x1mean = datamatrix(optionsp->get_samplesize(),gridsizex,0);
      x2mean = datamatrix(optionsp->get_samplesize(),gridsizey,0);
      interact = datamatrix(optionsp->get_samplesize(),gridsize,0);

      fchelp.readsample3(interact);

      for(i=0;i<gridsize;i++)
        for(j=0;j<xvalues.rows();j++)
          if(effectvaluesx[i]==xvalues(j,0))
            for(k=0;k<optionsp->get_samplesize();k++)
              x1mean(k,j) += interact(k,i)/gridsizey;
      for(i=0;i<gridsize;i++)
        for(j=0;j<yvalues.rows();j++)
          if(effectvaluesy[i]==yvalues(j,0))
            for(k=0;k<optionsp->get_samplesize();k++)
              x2mean(k,j) += interact(k,i)/gridsizex;

      for(i=0;i<gridsize;i++)
        for(j=0;j<xvalues.rows();j++)
          for(l=0;l<yvalues.rows();l++)
            if(effectvaluesx[i]==xvalues(j,0) && effectvaluesy[i]==yvalues(l,0))
              for(k=0;k<optionsp->get_samplesize();k++)
                interact(k,i) -= x1mean(k,j) + x2mean(k,l);

      for(i=0;i<x1mean.cols();i++)
        {

        l1i = x1mean.quantile(lower1,i);
        l2i = x1mean.quantile(lower2,i);
        u1i = x1mean.quantile(upper1,i);
        u2i = x1mean.quantile(upper2,i);

        xout << (i+1) << "   ";
        xout << xvalues(i,0) << "   ";
        xout << x1mean.mean(i) << "   ";
        xout << l1i << "   ";
        xout << l2i << "   ";
        xout << x1mean.quantile(50,i) << "   ";
        xout << u1i << "   ";
        xout << u2i << "   ";

        if (l1i > 0)
          xout << 1 << "   ";
        else if (u2i < 0)
          xout << -1 << "   ";
        else
          xout << 0 << "   ";

        if (l2i > 0)
          xout << 1 << "   ";
        else if (u1i < 0)
          xout << -1 << "   ";
        else
          xout << 0 << "   ";

        xout << endl;
        }

      for(i=0;i<x2mean.cols();i++)
        {

        l1i = x2mean.quantile(lower1,i);
        l2i = x2mean.quantile(lower2,i);
        u1i = x2mean.quantile(upper1,i);
        u2i = x2mean.quantile(upper2,i);

        yout << (i+1) << "   ";
        yout << yvalues(i,0) << "   ";
        yout << x2mean.mean(i) << "   ";
        yout << l1i << "   ";
        yout << l2i << "   ";
        yout << x2mean.quantile(50,i) << "   ";
        yout << u1i << "   ";
        yout << u2i << "   ";

        if (l1i > 0)
          yout << 1 << "   ";
        else if (u2i < 0)
          yout << -1 << "   ";
        else
          yout << 0 << "   ";

        if (l2i > 0)
          yout << 1 << "   ";
        else if (u1i < 0)
          yout << -1 << "   ";
        else
          yout << 0 << "   ";

        yout << endl;
        }

      effitx = effectvaluesx.begin();
      effity = effectvaluesy.begin();
      for(i=0;i<interact.cols();i++,effitx++,effity++)
        {

        l1i = interact.quantile(lower1,i);
        l2i = interact.quantile(lower2,i);
        u1i = interact.quantile(upper1,i);
        u2i = interact.quantile(upper2,i);

        interactout << (i+1) << "   ";
        interactout << *effitx << "   ";
        interactout << *effity << "   ";
        interactout << interact.mean(i) << "   ";
        interactout << l1i << "   ";
        interactout << l2i << "   ";
        interactout << interact.quantile(50,i) << "   ";
        interactout << u1i << "   ";
        interactout << u2i << "   ";

        if (l1i > 0)
          interactout << 1 << "   ";
        else if (u2i < 0)
          interactout << -1 << "   ";
        else
          interactout << 0 << "   ";

        if (l2i > 0)
          interactout << 1 << "   ";
        else if (u1i < 0)
          interactout << -1 << "   ";
        else
          interactout << 0 << "   ";

        interactout << endl;
        }

      } // ENDE: if(gridsize<0)
    } // ENDE: if(centertotal==false) else if(centertotal)
  } // ENDE: if(!mapexisting)

  }


void spline_basis_surf::add_linearpred_multBS_Block(const datamatrix & b, const unsigned a,const unsigned e,
                                                    const unsigned beg,const unsigned end)
  {

  datamatrix *workl = &(likep->get_linearpred(true));
  int *workindex;

  unsigned i,j,k,pos;

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator firstit = first.begin();

  workindex = index.getV() + beg;
  freqwork = freq.begin() + beg;
  firstit = first.begin() + beg;

  for(i=beg;i<=end;i++,freqwork++,workindex++,firstit++)
    {
    for(j=0;j<degree+1;j++)
      {
      for(k=0;k<degree+1;k++)
        {
        pos = *firstit+j+k*nrpar1dim;
        if(a <= pos && pos <= e)
          (*workl)(*workindex,column) += B(*freqwork,j+k*(degree+1)) * b(pos-a,0);
        }
      }
    }

  }


void spline_basis_surf::multBS_index_Block(datamatrix & res,const datamatrix & b, const unsigned a,const unsigned e,
                                                            const unsigned beg,const unsigned end)
  {

  int *workindex;
  double *work;

  unsigned i,j,k,pos;

  vector<int>::iterator freqwork;
  vector<int>::iterator firstit;

  work = res.getV();
  for(i=0;i<res.rows();i++,++work)
    *work = 0.0;

  workindex = index.getV() + beg;
  freqwork = freq.begin() + beg;
  firstit = first.begin() + beg;

  for(i=beg;i<=end;i++,freqwork++,workindex++,firstit++)
    {
    for(j=0;j<degree+1;j++)
      {
      for(k=0;k<degree+1;k++)
        {
        pos = *firstit+j+k*nrpar1dim;
        if(a <= pos && pos <= e)
          res(*workindex,0) += B(*freqwork,j+k*(degree+1)) * b(pos,0);
        }
      }
    }

  }


void spline_basis_surf::multBS(datamatrix & res, const datamatrix & b)
  {

  unsigned i,j,k;
  double *workres;
  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator firstit = first.begin();

  workres = res.getV();
  for(i=0;i<res.rows()*res.cols();i++,workres++)
    *workres = 0.0;

  workres = res.getV();
  for(i=0;i<likep->get_nrobs();i++,firstit++,freqwork++,workres++)
    {
    for(j=0;j<degree+1;j++)
      {
      for(k=0;k<degree+1;k++)
        {
        *workres += B(*freqwork,j+k*(degree+1)) * b(*firstit + j + k*nrpar1dim,0);
        }
      }
    }

  }


void spline_basis_surf::multBout(datamatrix & res, const datamatrix & b)
  {

  unsigned i,j,k;
  double *workres;
  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator firstit = first.begin();

  workres = res.getV();
  for(i=0;i<res.rows()*res.cols();i++,workres++)
    *workres = 0.0;

  workres = res.getV();
  for(i=0;i<likep->get_nrobs();i++,firstit++,freqwork++,workres++)
    {
    for(j=0;j<degree+1;j++)
      {
      for(k=0;k<degree+1;k++)
        {
        *workres += Bout(*freqwork,j+k*(degree+1)) * b(*firstit + j + k*nrpar1dim,0);
        }
      }
    }

  }


void spline_basis_surf::multBS_index(datamatrix & res, const datamatrix & b)
  {
/*
  unsigned i,j,k;
  double *workres;
  int *workindex;
  vector<int>::iterator freqwork = freq.begin();

  workres = res.getV();
  for(i=0;i<res.rows()*res.cols();i++,workres++)
    *workres = 0.0;

  workindex = index.getV();
  for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
    {
    for(j=0;j<degree+1;j++)
      {
      for(k=0;k<degree+1;k++)
        {
        res(*workindex,0) += B(*freqwork,j+k*(degree+1)) * b(first[i] + j + k*nrpar1dim,0);
        }
      }
    }
*/

  int i;
  unsigned j,k;
  double val=0.0;

  double *workbeta=NULL;
  double *workB = B.getV();
  int *workindex = index.getV();

  int maxfirst = first[res.rows()-1];

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator firstit = first.begin();

  i = 0;
  while(i<=maxfirst)
    {

    workbeta = b.getV() + i;
    while(*firstit==i)
      {
      if( freqwork==freq.begin() || *freqwork!=*(freqwork-1) )
        {
        val = 0.0;
        for(j=0;j<degree+1;j++)
          for(k=0;k<degree+1;k++,workB++)
            val += *workB * *(workbeta+k+j*nrpar1dim);
        }
      res(*workindex,0) = val;
      workindex++;
      freqwork++;
      firstit++;
      }
    i++;

    }

  }


void spline_basis_surf::multDG(datamatrix & res, const datamatrix & b)
  {

  int i;
  unsigned j,k;
  double *workres;

  vector<int>::iterator firstit = DGfirst.begin();

  workres = res.getV();
  for(i=0;i<gridsize;i++,workres++)
    *workres = 0.0;

  workres = res.getV();
  for(i=0;i<gridsize;i++,++workres,++firstit)
    {
    for(j=0;j<degree+1;j++)
      {
      for(k=0;k<degree+1;k++)
        {
        *workres += DG(i,j+k*(degree+1)) * b(*firstit + j + k*nrpar1dim,0);
        }
      }
    }

  }


void spline_basis_surf::predict(const datamatrix & newX, datamatrix & linpred)
  {

  assert(newX.rows() == 1 && newX.cols() == 2);
  unsigned i,j;

  datamatrix betac(beta.rows(),beta.cols());
  datamatrix bspline1(1,nrpar1dim,0);
  datamatrix bspline2(1,nrpar1dim,0);
  datamatrix bspline(1,nrpar,0);
  double * worklin = linpred.getV();

  for(i=0;i<nrpar1dim;i++)
    {
    bspline1(0,i) = bspline_rek(degree,i,newX(0,0),true);
    bspline2(0,i) = bspline_rek(degree,i,newX(0,1),false);
    }

  for(i=0;i<nrpar1dim;i++)
    for(j=0;j<nrpar1dim;j++)
      bspline(0,j + i*nrpar1dim) = bspline1(0,j)*bspline2(0,i);

  for(i=0;i<optionsp->get_samplesize();i++,worklin++)
    {
    readsample2(betac,i);
    for(j=0;j<nrpar;j++)
      {
      *worklin += betac(j,0) * bspline(0,j);
      }
    }

  }


void spline_basis_surf::compute_intercept(void)
  {

  double *workbeta = beta.getV();
  double *workbetaweight = betaweight.getV();
  unsigned i;

  intercept = 0.0;
  for(i=0;i<nrpar;i++,workbetaweight++,workbeta++)
    intercept += *workbetaweight * *workbeta;

  }


void spline_basis_surf::compute_betaweight(void)
  {

  unsigned i,j;

  betaweight = datamatrix(nrpar,1,0);
  betaweight_main = datamatrix(nrpar1dim,1,1);

  datamatrix bw = datamatrix(nrpar1dim,1,1);

  if(knpos == equidistant && (degree==1 || degree==2 || degree==3))
    {
    if(degree==1)
      {
      bw(0,0) = 0.5;
      bw(nrpar1dim-1,0) = 0.5;
      }
    else if(degree==2)
      {
      bw(0,0) = 1/6.0;
      bw(nrpar1dim-1,0) = 1/6.0;
      bw(1,0) = 5/6.0;
      bw(nrpar1dim-2,0) = 5/6.0;
      }
    else if(degree==3)
      {
      bw(0,0) = 1/24.0;
      bw(nrpar1dim-1,0) = 1/24.0;
      bw(1,0) = 12/24.0;
      bw(nrpar1dim-2,0) = 12/24.0;
      bw(2,0) = 23/24.0;
      bw(nrpar1dim-3,0) = 23/24.0;
      }

    for(i=0;i<nrpar1dim;i++)
      bw(i,0) = bw(i,0)/double(nrknots-1);

    }
  else
    {
    bw = datamatrix(nrpar1dim,1,1.0/double(nrpar1dim));
    }

  betaweight_main = bw;

  for(i=0;i<nrpar1dim;i++)
    for(j=0;j<nrpar1dim;j++)
      betaweight(i+j*nrpar1dim,0) = bw(i,0) * bw(j,0);

  }


void spline_basis_surf::compute_betaweightxy(void)
  {
/*
  unsigned i,j,k;

  datamatrix help(nrpar,1);

  betaweightx = datamatrix(xv.size(),nrpar,0);
  betaweighty = datamatrix(yv.size(),nrpar,0);

  if(knpos == equidistant && (degree==1 || degree==2 || degree==3))
    {
    datamatrix bw = datamatrix(nrpar1dim,1,1);
    if(degree==1)
      {
      bw(0,0) = 0.5;
      bw(nrpar1dim-1,0) = 0.5;
      }
    else if(degree==2)
      {
      bw(0,0) = 1/6.0;
      bw(nrpar1dim-1,0) = 1/6.0;
      bw(1,0) = 5/6.0;
      bw(nrpar1dim-2,0) = 5/6.0;
      }
    else if(degree==3)
      {
      bw(0,0) = 1/24.0;
      bw(nrpar1dim-1,0) = 1/24.0;
      bw(1,0) = 12/24.0;
      bw(nrpar1dim-2,0) = 12/24.0;
      bw(2,0) = 23/24.0;
      bw(nrpar1dim-3,0) = 23/24.0;
      }

    datamatrix b;

// betaweightx
    for(i=0;i<xv.size();i++)
      {
      b = bspline(xv[i],knot1);
      for(j=0;j<nrpar1dim;j++)
        for(k=0;k<nrpar1dim;k++)
          betaweightx(i,j+k*nrpar1dim) = (b(j,0) * bw(k,0))/double(nrknots-1);
      }
// betaweighty
    for(i=0;i<yv.size();i++)
      {
      b = bspline(yv[i],knot2);
      for(j=0;j<nrpar1dim;j++)
        for(k=0;k<nrpar1dim;k++)
          betaweighty(i,j+k*nrpar1dim) = (bw(j,0) * b(k,0))/double(nrknots-1);
      }

    }   // ENDE: if
  else  // Sehnentrapezregel an 'nr' \E4quidistanten Knoten zwischen x(y)min und x(y)max
    {
    unsigned nr = 30;

    datamatrix kx = datamatrix(nr,1,0);
    datamatrix ky = datamatrix(nr,1,0);

    for(i=0;i<nr;i++)
      {
      kx(i,0) = knot1[degree] + i*(knot1[degree+nrknots-1]-knot1[degree])/(nr-1);
      ky(i,0) = knot2[degree] + i*(knot2[degree+nrknots-1]-knot2[degree])/(nr-1);
      }

// betaweightx
    for(i=0;i<xv.size();i++)
      {
      help = bspline(xv[i],ky(0,0));
      help.plus(help,bspline(xv[i],ky(nr-1,0)));
      for(k=0;k<nrpar;k++)
        betaweightx(i,k) += 0.5*help(k,0);
      for(j=1;j<nr-1;j++)
        {
        help = bspline(xv[i],ky(j,0));
        for(k=0;k<nrpar;k++)
          betaweightx(i,k) += help(k,0);
        }
      }
// betaweighty
    for(i=0;i<yv.size();i++)
      {
      help = bspline(kx(0,0),yv[i]);
      help.plus(help,bspline(kx(nr-1,0),yv[i]));
      for(k=0;k<nrpar;k++)
        betaweighty(i,k) += 0.5*help(k,0);
      for(j=1;j<nr-1;j++)
        {
        help = bspline(kx(j,0),yv[i]);
        for(k=0;k<nrpar;k++)
          betaweighty(i,k) += help(k,0);
        }
      }

    for(i=0;i<betaweightx.rows();i++)
      for(j=0;j<nrpar;j++)
        betaweightx(i,j) /= (nr-1);
    for(i=0;i<betaweighty.rows();i++)
      for(j=0;j<nrpar;j++)
        betaweighty(i,j) /= (nr-1);

    }   // ENDE: else
*/

  unsigned i,j,k;

  datamatrix help(nrpar,1);

  betaweightx = datamatrix(xv.size(),nrpar,0);
  betaweighty = datamatrix(yv.size(),nrpar,0);

  datamatrix bw = datamatrix(nrpar1dim,1,1);

  if(knpos == equidistant && (degree==1 || degree==2 || degree==3))
    {
    if(degree==1)
      {
      bw(0,0) = 0.5;
      bw(nrpar1dim-1,0) = 0.5;
      }
    else if(degree==2)
      {
      bw(0,0) = 1/6.0;
      bw(nrpar1dim-1,0) = 1/6.0;
      bw(1,0) = 5/6.0;
      bw(nrpar1dim-2,0) = 5/6.0;
      }
    else if(degree==3)
      {
      bw(0,0) = 1/24.0;
      bw(nrpar1dim-1,0) = 1/24.0;
      bw(1,0) = 12/24.0;
      bw(nrpar1dim-2,0) = 12/24.0;
      bw(2,0) = 23/24.0;
      bw(nrpar1dim-3,0) = 23/24.0;
      }
    }
  else
    {
    double w = (nrknots-1)/double(nrpar1dim);
    bw = datamatrix(nrpar1dim,1,w);
    }

  datamatrix b;

// betaweightx
  for(i=0;i<xv.size();i++)
    {
    b = bspline(xv[i],knot1);
    for(j=0;j<nrpar1dim;j++)
      for(k=0;k<nrpar1dim;k++)
        betaweightx(i,j+k*nrpar1dim) = (b(j,0) * bw(k,0))/double(nrknots-1);
    }
// betaweighty
  for(i=0;i<yv.size();i++)
    {
    b = bspline(yv[i],knot2);
    for(j=0;j<nrpar1dim;j++)
      for(k=0;k<nrpar1dim;k++)
        betaweighty(i,j+k*nrpar1dim) = (bw(j,0) * b(k,0))/double(nrknots-1);
    }

  }


void spline_basis_surf::compute_beta(void)
  {

  unsigned i,j;

//* ------------- betas zeilen- und spaltenweise zentrieren ---------------------

  for(i=0;i<nrpar1dim;i++)
    {
    beta1(i,0) = 0;
    beta2(i,0) = 0;
    }

// betax1mean, betax2mean
  for(i=0;i<nrpar1dim;i++)
    for(j=0;j<nrpar1dim;j++)
      beta1(i,0) += beta(i+j*nrpar1dim,0)*betaweight_main(j,0);

  for(i=0;i<nrpar1dim;i++)
    for(j=0;j<nrpar1dim;j++)
      beta2(i,0) += beta(j+i*nrpar1dim,0)*betaweight_main(j,0);

// betas \E4ndern
  for(i=0;i<nrpar1dim;i++)
    for(j=0;j<nrpar1dim;j++)
      beta(i+j*nrpar1dim,0) -= beta1(i,0);

  for(i=0;i<nrpar1dim;i++)
    for(j=0;j<nrpar1dim;j++)
      beta(j+i*nrpar1dim,0) -= beta2(i,0);

  for(i=0;i<nrpar;i++)
    beta(i,0) += intercept;

// ------------- ENDE: betas zeilen- und spaltenweise zentrieren -------------*/

  } // end: compute_beta()


void spline_basis_surf::compute_main(void)
  {

  unsigned i;

  double *workspline;
  int *workindex = index.getV();
  vector<int>::iterator freqwork;

// Haupteffekte berechnen
  he1.mult(betaweightx,beta);
  he2.mult(betaweighty,beta);

// 'spline' \E4ndern
  freqwork = mainp1->get_freqit();
  workindex = mainp1->get_indexp();
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    spline(*workindex,0) -= he1(*freqwork,0);
  freqwork = mainp2->get_freqit();
  workindex = mainp2->get_indexp();
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    spline(*workindex,0) -= he2(*freqwork,0);

  workspline = spline.getV();
  for(i=0;i<spline.rows();i++,workspline++)
    *workspline += intercept;

  }


void spline_basis_surf::init_names(const vector<ST::string> & na)
  {
  FULLCOND::init_names(na);

  ST::string underscore = "\\_";

  if(!varcoeff)
    {
    if(mapexisting)
      {
      ST::string helpname0 = na[0].insert_string_char('_',underscore);
      term_symbolic = "f_{" + helpname0 + "}(" + helpname0 + ")";
      }
    else
      {
      ST::string helpname0 = na[0].insert_string_char('_',underscore);
      ST::string helpname1 = na[1].insert_string_char('_',underscore);
      term_symbolic = "f_{" +  helpname0 + "," + helpname1 + "}(" + helpname0 + "," + helpname1 + ")";
      }
    }
  else
    {
    if(mapexisting)
      {
      ST::string helpname0 = na[0].insert_string_char('_',underscore);
      ST::string helpname1 = na[1].insert_string_char('_',underscore);
      term_symbolic = "f_{" + helpname0 + "}(" + helpname0 + ")" + " \\cdot " + helpname1;
      }
    else
      {
      ST::string helpname0 = na[0].insert_string_char('_',underscore);
      ST::string helpname1 = na[2].insert_string_char('_',underscore);
      ST::string helpname2 = na[1].insert_string_char('_',underscore);
      term_symbolic = "f_{" +  helpname0 + "," + helpname1 + "}(" + helpname0 + "," + helpname1 + ")" + " \\cdot " + helpname2;
      }

    }
  if (column > 0)
    priorassumptions.push_back("$" + term_symbolic + "$" +
     " (" + ST::inttostring(column+1) + ". response category)");
  else
    priorassumptions.push_back("$" + term_symbolic + "$");

  if(type==MCMC::mrfkronecker)
    priorassumptions.push_back("P-spline with Kronecker product interaction penalty");
  else if(type==MCMC::mrflinear)
    {
    priorassumptions.push_back("P-spline with 2 dimensional first order random walk penalty");
    priorassumptions.push_back("(Kronecker sum of two first order random walks)");
    }
  else if(type==MCMC::mrfquadratic8)
    {
    priorassumptions.push_back("P-spline with 2 dimensional second order random walk penalty");
    priorassumptions.push_back("(Kronecker sum of two second order random walks)");
    }
  else if(type==MCMC::mrfquadratic12)
    {
    priorassumptions.push_back("P-spline with 2 dimensional second order random walk penalty");
    priorassumptions.push_back("(Approximation to the biharmonic differential operator)");
    }
  else if(type==MCMC::mrfkr1)
    priorassumptions.push_back("P-spline with Kronecker product interaction (RW1*RW1) penalty");
  else if(type==MCMC::mrfkr2)
    priorassumptions.push_back("P-spline with Kronecker product interaction (RW2*RW2) penalty");

  ST::string knotstr;
  if (knpos == equidistant)
    knotstr = "equidistant";
  else if (knpos == quantiles)
    knotstr = "quantiles";
  priorassumptions.push_back("Number of knots: " + ST::inttostring(nrknots));
  priorassumptions.push_back("Knot choice: " + knotstr);
  priorassumptions.push_back("Degree of Splines: " + ST::inttostring(degree));
  }


void spline_basis_surf::set_lambdaconst(double la)
  {
  lambda=la;
  lambdaconst = true;
  }

bool spline_basis_surf::breakpause(void)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
    Application->ProcessMessages();

    if (Frame->pause)
      {
      optionsp->out("\n");
      optionsp->out("SIMULATION PAUSED\n");
      optionsp->out("Click CONTINUE to proceed\n");
      optionsp->out("\n");

      while (Frame->pause)
        {
        Application->ProcessMessages();
        }

      optionsp->out("SIMULATION CONTINUED\n");
      optionsp->out("\n");
      }

    return Frame->stop;
#elif defined(JAVA_OUTPUT_WINDOW)
    return optionsp->adminb_p->breakcommand();
#endif
    return true;
  }


void spline_basis_surf::createreml(datamatrix & X,datamatrix & Z,
                                const unsigned & Xpos, const unsigned & Zpos)
  {
  unsigned i,j;

  double * workdata;
  double * workZ;

  datamatrix Kstat;
  if(type == mrflinear)
    {
    Kstat=STATMAT_PENALTY::K2dim_pspline(nrpar1dim);
    datamatrix vals(Kstat.rows(),1,0);

    bool eigentest=eigen2(Kstat,vals);
    if(eigentest==false)
      {
      errors.push_back("ERROR: Unable to compute eigen decomposition for 2 dimensional P-spline.\n");
      }
    else
      {
      eigensort(vals,Kstat);
      for(i=0; i<vals.rows()-1; i++)
        {
        vals(i,0)=1/sqrt(vals(i,0));
        }
      vals(vals.rows()-1,0)=0;
      Kstat = multdiagback(Kstat,vals).getColBlock(0,Kstat.cols()-1);
      }
    }
  else if(type == mrfquadratic8 || type == mrfquadratic12)
    {
    if(type == mrfquadratic8)
      {
      Kstat=STATMAT_PENALTY::K2dim_pspline_rw2(nrpar1dim,2,2);
      }
    else if(type == mrfquadratic12)
      {
      Kstat=STATMAT_PENALTY::K2dim_pspline_biharmonic(nrpar1dim);
      }
    datamatrix vals(Kstat.rows(),1,0);

    bool eigentest=eigen2(Kstat,vals);
    if(eigentest==false)
      {
      errors.push_back("ERROR: Unable to compute eigen decomposition for 2 dimensional P-spline.\n");
      }
    else
      {
      eigensort(vals,Kstat);
      for(i=0; i < dimZ; i++)
        {
        vals(i,0)=1/sqrt(vals(i,0));
        }
      for(i=0; i < dimX+1; i++)
        {
        vals(dimZ+i,0)=0;
        }
      Kstat = multdiagback(Kstat,vals).getColBlock(0,dimZ);
      }
    }

/*
// BEGIN: Charris
Kstat = statmatrix<double>::diag(nrpar,1.0);
datamatrix design(X.rows(),nrpar,0);
datamatrix h1(nrpar,1,0);
datamatrix h2(X.rows(),1,0);

for(i=0; i<nrpar; i++)
  {
  h1 = Kstat.getCol(i);
  multBS_index(h2,h1);
  for(j=0; j<X.rows(); j++)
    {
    design(j,i) = h2(j,0);
    }
  }

ofstream out1("c:\\temp\\design.raw");
design.prettyPrint(out1);
out1.close();
// END: Charris
*/

/*
// BEGIN: Susanne
  Z_VCM = Kstat;
  datamatrix knoten1(knot1.size(),1,0);
  datamatrix knoten2(knot1.size(),1,0);
  for(i=0; i<knot1.size(); i++)
    {
    knoten1(i,0)=knot1[i];
    knoten2(i,0)=knot2[i];
    }
  ofstream out1((pathcurrent.substr(0,pathcurrent.length()-4) + "_knots1.res").strtochar());
  knoten1.prettyPrint(out1);
  out1.close();
  ofstream out2((pathcurrent.substr(0,pathcurrent.length()-4) + "_knots2.res").strtochar());
  knoten2.prettyPrint(out2);
  out2.close();
//END: Susanne
*/
  datamatrix spline2;
  if(gridsize>0)
    {
    spline2 = datamatrix(gridsize,1,0);
    }

// X berechnen (varcoeff)
  if(varcoeff & !centervcm)
    {
    double * workX;
    unsigned Xcols = X.cols();
    workX = X.getV()+Xpos;
    double * workintact = data_forfixed.getV();
    for (i=0;i<spline.rows();i++,workintact++,workX+=Xcols)
      {
      *workX = *workintact;
      }
    }
// X berechnen (mrfquadratic8)

  if(type==mrfquadratic8)
    {
    datamatrix knoten1(nrpar1dim,1,0);
    datamatrix knoten2(nrpar1dim,1,0);
    for(i=0; i<nrpar1dim; i++)
      {
      knoten1(i,0)=knot1[i];
      knoten2(i,0)=knot2[i];
      }
    statmatrix<double> I = statmatrix<double>(nrpar1dim,1,1);
    knoten1 = kronecker(I,knoten1);
    knoten2 = kronecker(knoten2,I);

    datamatrix knoten3(knoten1.rows(),1,0);
    for(i=0; i<knoten1.rows(); i++)
      {
      knoten3(i,0) = knoten1(i,0)*knoten2(i,0);
      }

    double mean=0;

    multBS_index(spline,knoten1);
    mean=spline.mean(0);
    for(i=0; i<spline.rows(); i++)
      {
      X(i,Xpos) = spline(i,0)-mean;
      }
    if(gridsize>0)
      {
      multDG(spline2,knoten1);
      for(i=0; i<spline2.rows(); i++)
        {
        X_grid(i,0) = spline2(i,0)-mean;
        }
      }

    multBS_index(spline,knoten2);
    mean=spline.mean(0);
    for(i=0; i<spline.rows(); i++)
      {
      X(i,Xpos+1) = spline(i,0)-mean;
      }
    if(gridsize>0)
      {
      multDG(spline2,knoten2);
      for(i=0; i<spline2.rows(); i++)
        {
        X_grid(i,1) = spline2(i,0)-mean;
        }
      }

    multBS_index(spline,knoten3);
    mean=spline.mean(0);
    for(i=0; i<spline.rows(); i++)
      {
      X(i,Xpos+2) = spline(i,0)-mean;
      }
    if(gridsize>0)
      {
      multDG(spline2,knoten3);
      for(i=0; i<spline2.rows(); i++)
        {
        X_grid(i,2) = spline2(i,0)-mean;
        }
      }
    }

  if(type==mrfquadratic12)
    {
    datamatrix knoten1(nrpar1dim,1,0);
    datamatrix knoten2(nrpar1dim,1,0);
    for(i=0; i<nrpar1dim; i++)
      {
      knoten1(i,0)=knot1[i];
      knoten2(i,0)=knot2[i];
      }
    statmatrix<double> I = statmatrix<double>(nrpar1dim,1,1);
    knoten1 = kronecker(I,knoten1);
    knoten2 = kronecker(knoten2,I);

    double mean=0;

    multBS_index(spline,knoten1);
    mean=spline.mean(0);
    for(i=0; i<spline.rows(); i++)
      {
      X(i,Xpos) = spline(i,0)-mean;
      }
    if(gridsize>0)
      {
      multDG(spline2,knoten1);
      for(i=0; i<spline2.rows(); i++)
        {
        X_grid(i,0) = spline2(i,0)-mean;
        }
      }

    multBS_index(spline,knoten2);
    mean=spline.mean(0);
    for(i=0; i<spline.rows(); i++)
      {
      X(i,Xpos+1) = spline(i,0)-mean;
      }
    if(gridsize>0)
      {
      multDG(spline2,knoten2);
      for(i=0; i<spline2.rows(); i++)
        {
        X_grid(i,1) = spline2(i,0)-mean;
        }
      }
    }

// Z berechnen

  unsigned Zcols = Z.cols();
  for(j=0;j<dimZ;j++)
    {
    multBS_index(spline,Kstat.getCol(j));

    workdata = spline.getV();
    workZ = Z.getV()+Zpos+j;

    if(varcoeff)
      {
      double * workintact = data_forfixed.getV();
      double * workZ_VCM = Z_VCM.getV()+j;
      for (i=0;i<spline.rows();i++,workdata++,workintact++,workZ+=Zcols,workZ_VCM+=dimZ)
        {
        *workZ = *workdata**workintact;
        *workZ_VCM = *workdata;
        }
      }
    else
      {
      for (i=0;i<spline.rows();i++,workdata++,workZ+=Zcols)
        {
        *workZ = *workdata;
        }
      }
    }

  if(gridsize>0)
    {
    for(j=0;j<dimZ;j++)
      {
      multDG(spline2,Kstat.getCol(j));

      workdata = spline2.getV();
      workZ = Z_grid.getV()+j;

      for (i=0;i<spline2.rows();i++,workdata++,workZ+=dimZ)
        {
        *workZ = *workdata;
        }
      }
    }
  }


double spline_basis_surf::outresultsreml(datamatrix & X,datamatrix & Z,
                                  datamatrix & betareml,datamatrix & betacov,
                                  datamatrix & thetareml,
                                  const unsigned & Xpos,
                                  const unsigned & Zpos,
                                  const unsigned & thetapos,
                                  const bool & dispers,
                                  const unsigned & betaXpos,
                                  const unsigned & betaZpos,
                                  const double & category,
                                  const bool & ismultinomial,
                                  const unsigned plotpos)
  {
  double mean=0;
  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  unsigned nr = effectvaluesx.size();
  unsigned i,j;

  betamean=datamatrix(nr,1,0);
  datamatrix betastd=datamatrix(nr,1,0);
  betaqu_l1_lower=datamatrix(nr,1,0);
  betaqu_l1_upper=datamatrix(nr,1,0);
  betaqu_l2_lower=datamatrix(nr,1,0);
  betaqu_l2_upper=datamatrix(nr,1,0);

  vector<int>::iterator indexit = index2.begin();
  unsigned k = *indexit;
  vector<int>::iterator freqwork = freqoutput.begin();

  if(varcoeff)
    {
    for(i=0,j=0;i<spline.rows();i++,indexit++,freqwork++,k+=*indexit)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        if(!centervcm)
          {
          betamean(j,0) = betareml(betaXpos,0)*X_VCM(k,0) + (Z_VCM.getRow(k)*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1))(0,0);
          betastd(j,0) = sqrt(
                              (
                               X_VCM(k,0)*betacov(betaXpos,betaXpos)
                               +
                               (Z_VCM.getRow(k)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+1))(0,0)
                              )*X_VCM(k,0)
                              +
                              (
                               (
                                X_VCM(k,0)*betacov.getBlock(betaXpos,betaZpos,betaXpos+1,betaZpos+dimZ)
                                +
                                Z_VCM.getRow(k)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                               )*(Z_VCM.getRow(k).transposed())
                              )(0,0)
                             );
          j++;
          }
        else
          {
          betamean(j,0) = (Z_VCM.getRow(k)*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1))(0,0);
          betastd(j,0) = sqrt((Z_VCM.getRow(k)*
                               betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)*
                               (Z_VCM.getRow(k).transposed()))(0,0));
          j++;
          }
        }
      }
    }
  else
    {
    if(type==mrflinear)
      {
      for(i=0,j=0;i<spline.rows();i++,indexit++,freqwork++,k+=*indexit)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          betamean(j,0) = (Z.getBlock(k,Zpos,k+1,Zpos+nrpar-1)*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1))(0,0);
          betastd(j,0) = sqrt((Z.getBlock(k,Zpos,k+1,Zpos+nrpar-1)*
                   betacov.getBlock(betaZpos,betaZpos,betaZpos+nrpar-1,betaZpos+nrpar-1)*
                   Z.getBlock(k,Zpos,k+1,Zpos+nrpar-1).transposed())(0,0));
          j++;
          }
        }
      }
    else if(type==mrfquadratic8 || type==mrfquadratic12)
      {
      for(i=0,j=0;i<spline.rows();i++,indexit++,freqwork++,k+=*indexit)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          betamean(j,0) = (X.getBlock(k,Xpos,k+1,Xpos+dimX)*betareml.getBlock(betaXpos,0,betaXpos+dimX,1))(0,0) + (Z.getBlock(k,Zpos,k+1,Zpos+dimZ)*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1))(0,0);
          betastd(j,0) = sqrt(
                             ((
                              X.getBlock(k,Xpos,k+1,Xpos+dimX)*betacov.getBlock(betaXpos,betaXpos,betaXpos+dimX,betaXpos+dimX)
                              +
                              (Z.getBlock(k,Zpos,k+1,Zpos+dimZ)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+dimX))
                             )*X.getBlock(k,Xpos,k+1,Xpos+dimX).transposed())(0,0)
                             +
                             (
                              (
                               X.getBlock(k,Xpos,k+1,Xpos+dimX)*betacov.getBlock(betaXpos,betaZpos,betaXpos+dimX,betaZpos+dimZ)
                               +
                               Z.getBlock(k,Zpos,k+1,Zpos+dimZ)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                              )*(Z.getBlock(k,Zpos,k+1,Zpos+dimZ).transposed())
                             )(0,0)
                            );
          j++;
          }
        }
      }
    }

  if(!varcoeff)
    {
    mean = betamean.mean(0);
    for(j=0; j<nr; j++)
      {
      betamean(j,0)=betamean(j,0)-mean;
      betaqu_l1_lower(j,0) = betamean(j,0)+randnumbers::invPhi2(lower1/100)*betastd(j,0);
      betaqu_l1_upper(j,0) = betamean(j,0)+randnumbers::invPhi2(upper2/100)*betastd(j,0);
      betaqu_l2_lower(j,0) = betamean(j,0)+randnumbers::invPhi2(lower2/100)*betastd(j,0);
      betaqu_l2_upper(j,0) = betamean(j,0)+randnumbers::invPhi2(upper1/100)*betastd(j,0);
      }
    }
  else
    {
    for(j=0; j<nr; j++)
      {
      betaqu_l1_lower(j,0) = betamean(j,0)+randnumbers::invPhi2(lower1/100)*betastd(j,0);
      betaqu_l1_upper(j,0) = betamean(j,0)+randnumbers::invPhi2(upper2/100)*betastd(j,0);
      betaqu_l2_lower(j,0) = betamean(j,0)+randnumbers::invPhi2(lower2/100)*betastd(j,0);
      betaqu_l2_upper(j,0) = betamean(j,0)+randnumbers::invPhi2(upper1/100)*betastd(j,0);
      }
    }


  optionsp->out("\n");
  if(ismultinomial)
    {
    optionsp->out("  " + title + " (cat."+ST::doubletostring(category,6)+")\n",true);
    }
  else
    {
    optionsp->out("  " + title + "\n",true);
    }
  optionsp->out("\n");
  optionsp->out("\n");

  optionsp->out("  Estimated variance:  "
                + ST::doubletostring(thetareml(thetapos,0),6) + "\n");
  double smoothpar;
  if(dispers==true)
    {
    smoothpar = thetareml(thetareml.rows()-1,0)/thetareml(thetapos,0);
    optionsp->out("  Inverse variance:    "
                  + ST::doubletostring(1/thetareml(thetapos,0),6) + "\n");
    optionsp->out("  Smoothing parameter: "
                + ST::doubletostring(smoothpar,6) + "\n");
    optionsp->out("  (Smoothing parameter = scale / variance)\n");
    optionsp->out("  Degrees of freedom: "
                + ST::doubletostring(thetareml(thetapos,3),6) + "\n");
    }
  else
    {
    smoothpar = 1/thetareml(thetapos,0);
    optionsp->out("  Smoothing parameter: "
                + ST::doubletostring(smoothpar,6) + "\n");
    optionsp->out("  (Smoothing parameter = 1 / variance)\n");
    optionsp->out("  Degrees of freedom: "
                + ST::doubletostring(thetareml(thetapos,3),6) + "\n");
    }
  if(thetareml(thetapos,1)==1)
    {
    optionsp->out("  NOTE: Estimation of the variance was stopped after iteration "
                  + ST::doubletostring(thetareml(thetapos,2),0) + "\n");
    optionsp->out("        because the corresponding penalized part was small relative to the linear predictor.");
    }
  if(thetareml(thetapos,1)==-1)
    {
    optionsp->out("  NOTE: Estimation of the variance was stopped after iteration "
                  + ST::doubletostring(thetareml(thetapos,2),0) + "\n");
    optionsp->out("        to avoid numerical problems due to large variances.");
    }
  ST::string varpath=pathcurrent.substr(0,pathcurrent.length()-4) + "_var.res";
  if(ismultinomial)
    {
    varpath=varpath.insert_after_string(ST::doubletostring(category,6)+"_","_f_");
    }
  optionsp->out("\n");
  optionsp->out("  Variance and smoothing parameter are stored in file\n");
  optionsp->out("  " + varpath + "\n");

  ofstream outvarres(varpath.strtochar());
  outvarres << "variance  ";
  outvarres << "smoothpar  ";
  outvarres << "df  ";
  outvarres << "stopped  " <<endl;

  outvarres << thetareml(thetapos,0) <<"  ";
  outvarres << smoothpar <<"  ";
  outvarres << thetareml(thetapos,3) <<"  ";
  outvarres << (thetareml(thetapos,1)==1);
  outvarres << endl;
  outvarres.close();

  ST::string outest=pathcurrent;
  if(ismultinomial)
    {
    outest = pathcurrent.insert_after_string(ST::doubletostring(category,6)+"_","_f_");
    }
  ofstream outres(outest.strtochar());

  optionsp->out("\n");
  optionsp->out("  Results are stored in file\n");
  optionsp->out("  " + outest + "\n");
  optionsp->out("\n");

  if(mapexisting)
    {
#if defined(JAVA_OUTPUT_WINDOW)
    optionsp->out("  Postscript files are stored in files\n");
    ST::string psfile;
    psfile = outest.substr(0,outest.length()-4) + ".ps";
    optionsp->out("  " + psfile + "\n");
    psfile = outest.substr(0,outest.length()-4) + "_pcatbig" + ".ps";
    optionsp->out("  " + psfile + "\n");
    psfile = outest.substr(0,outest.length()-4) + "_pcatsmall" + ".ps";
    optionsp->out("  " + psfile + "\n");
    optionsp->out("\n");
    optionsp->out("  Results may be visualized using method 'drawmap'\n");
    optionsp->out("  Type for example:\n");
    optionsp->out("  objectname.drawmap " + ST::inttostring(plotpos) + "\n");
#else
    optionsp->out("  Results may be visualized using the R function 'drawmap' \n");
#endif
    optionsp->out("\n");
    }
  else
    {
    optionsp->out("  Results may be visualized using the R function 'plotsurf'\n");
    ST::string doublebackslash = "/";
    ST::string spluspath = outest.insert_string_char('\\',doublebackslash);
    optionsp->out("  Type for example:\n");
    optionsp->out("  plotsurf(\"" + spluspath + "\")");
    optionsp->out("\n");
    optionsp->out("\n");
    }

  assert(!outres.fail());
  outres << "intnr" << "   ";
  if(mapexisting)
    {
    outres << datanames[0] << "   ";
    outres << "x   " << "y   ";
    }
  else
    {
    outres << datanames[1] << "   ";
    outres << datanames[0] << "   ";
    }
  outres << "pmode   ";
  outres << "ci"  << level1  << "lower   ";
  outres << "ci"  << level2  << "lower   ";
  outres << "std   ";
  outres << "ci"  << level2  << "upper   ";
  outres << "ci"  << level1  << "upper   ";
  outres << "pcat" << level1 << "   ";
  outres << "pcat" << level2 << "   ";
  outres << endl;

  double * workmean = betamean.getV();
  double * workstd = betastd.getV();
  double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
  double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
  double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
  double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
  vector<double>::iterator workxvalues = effectvaluesx.begin();
  vector<double>::iterator workyvalues = effectvaluesy.begin();

  for(i=0;i<effectvaluesx.size();i++,workmean++,workstd++,
                     workbetaqu_l1_lower_p++,workbetaqu_l2_lower_p++,
                     workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                     workxvalues++,workyvalues++)
    {
    outres << (i+1) << "   ";
    if(mapexisting)
      outres << regionnames[i] << "   ";
    outres << *workxvalues << "   ";
    outres << *workyvalues << "   ";
    outres << *workmean << "   ";
//    outres << *workbetapval << "   ";
    outres << *workbetaqu_l1_lower_p << "   ";
    outres << *workbetaqu_l2_lower_p << "   ";
    outres << *workstd << "   ";
    outres << *workbetaqu_l2_upper_p << "   ";
    outres << *workbetaqu_l1_upper_p << "   ";
    if (*workbetaqu_l1_lower_p > 0)
      outres << "1   ";
    else if (*workbetaqu_l1_upper_p < 0)
      outres << "-1   ";
    else
      outres << "0   ";

    if (*workbetaqu_l2_lower_p > 0)
      outres << "1   ";
    else if (*workbetaqu_l2_upper_p < 0)
      outres << "-1   ";
    else
      outres << "0   ";

    outres << endl;
    }

  if(gridsize>0)
    {
    unsigned nr = gridsize;
    betamean=datamatrix(nr,1,0);
    datamatrix betastd=datamatrix(nr,1,0);
    betaqu_l1_lower=datamatrix(nr,1,0);
    betaqu_l1_upper=datamatrix(nr,1,0);
    betaqu_l2_lower=datamatrix(nr,1,0);
    betaqu_l2_upper=datamatrix(nr,1,0);

    if(type==mrflinear)
      {
      betamean = Z_grid*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1);
      for(i=0; i<gridsize; i++)
        {
        betamean(i,0)=betamean(i,0)-mean;
        betastd(i,0) = sqrt((Z_grid.getRow(i)*
                 betacov.getBlock(betaZpos,betaZpos,betaZpos+nrpar-1,betaZpos+nrpar-1)*
                 Z_grid.getRow(i).transposed())(0,0));
        }
      }
    else if(type==mrfquadratic8 || type==mrfquadratic12)
      {
      betamean = X_grid*betareml.getBlock(betaXpos,0,betaXpos+dimX,1) + Z_grid*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1);
      for(i=0; i<gridsize; i++)
        {
        betamean(i,0)=betamean(i,0)-mean;
        betastd(i,0) = sqrt(
                           ((
                            X_grid.getRow(i)*betacov.getBlock(betaXpos,betaXpos,betaXpos+dimX,betaXpos+dimX)
                            +
                            (Z_grid.getRow(i)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+dimX))
                           )*X_grid.getRow(i).transposed())(0,0)
                           +
                           (
                            (
                             X_grid.getRow(i)*betacov.getBlock(betaXpos,betaZpos,betaXpos+dimX,betaZpos+dimZ)
                             +
                             Z_grid.getRow(i)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                            )*(Z_grid.getRow(i).transposed())
                           )(0,0)
                          );
        }
      }

     for(j=0; j<nr; j++)
       {
       betaqu_l1_lower(j,0) = betamean(j,0)+randnumbers::invPhi2(lower1/100)*betastd(j,0);
       betaqu_l1_upper(j,0) = betamean(j,0)+randnumbers::invPhi2(upper2/100)*betastd(j,0);
       betaqu_l2_lower(j,0) = betamean(j,0)+randnumbers::invPhi2(lower2/100)*betastd(j,0);
       betaqu_l2_upper(j,0) = betamean(j,0)+randnumbers::invPhi2(upper1/100)*betastd(j,0);
       }

    outest = outest.substr(0,outest.length()-4) + "_grid.res";
    ofstream outgrid(outest.strtochar());

    optionsp->out("  Results on a grid are stored in file\n");
    optionsp->out("  " + outest + "\n");
    optionsp->out("\n");

    optionsp->out("  Results may be visualized using the R function 'plotsurf'\n");
    ST::string doublebackslash = "/";
    ST::string spluspath = outest.insert_string_char('\\',doublebackslash);
    optionsp->out("  Type for example:\n");
    optionsp->out("  plotsurf(\"" + spluspath + "\")");
    optionsp->out("\n");
    optionsp->out("\n");

    assert(!outgrid.fail());
    outgrid << "intnr" << "   ";
    if(mapexisting)
     {
     outgrid << "x   " << "y   ";
     }
    else
      {
      outgrid << datanames[1] << "   ";
      outgrid << datanames[0] << "   ";
      }
    outgrid << "pmode   ";
    outgrid << "ci"  << level1  << "lower   ";
    outgrid << "ci"  << level2  << "lower   ";
    outgrid << "std   ";
    outgrid << "ci"  << level2  << "upper   ";
    outgrid << "ci"  << level1  << "upper   ";
    outgrid << "pcat" << level1 << "   ";
    outgrid << "pcat" << level2 << "   ";
    outgrid << endl;

    workmean = betamean.getV();
    workstd = betastd.getV();
    workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
    workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
    workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
    workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
    workxvalues = effectvaluesxgrid.begin();
    workyvalues = effectvaluesygrid.begin();

    for(i=0;i<gridsize;i++,workmean++,workstd++,
                       workbetaqu_l1_lower_p++,workbetaqu_l2_lower_p++,
                       workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                       workxvalues++,workyvalues++)
      {
      outgrid << (i+1) << "   ";
      outgrid << *workxvalues << "   ";
      outgrid << *workyvalues << "   ";
      outgrid << *workmean << "   ";
      outgrid << *workbetaqu_l1_lower_p << "   ";
      outgrid << *workbetaqu_l2_lower_p << "   ";
      outgrid << *workstd << "   ";
      outgrid << *workbetaqu_l2_upper_p << "   ";
      outgrid << *workbetaqu_l1_upper_p << "   ";
      if (*workbetaqu_l1_lower_p > 0)
        outgrid << "1   ";
      else if (*workbetaqu_l1_upper_p < 0)
        outgrid << "-1   ";
      else
        outgrid << "0   ";

      if (*workbetaqu_l2_lower_p > 0)
        outgrid << "1   ";
      else if (*workbetaqu_l2_upper_p < 0)
        outgrid << "-1   ";
      else
        outgrid << "0   ";

      outgrid << endl;
      }
    }

/*
//BEGIN: Susanne
  datamatrix betaout = Z_VCM*betareml.getRowBlock(X.cols()+Zpos,X.cols()+Zpos+dimZ);
  ofstream outbetares((pathcurrent.substr(0,pathcurrent.length()-4) + "_beta.res").strtochar());
  betaout.prettyPrint(outbetares);
  outbetares.close();
//END: Susanne*/

  return mean;
  }

void spline_basis_surf::outresultsgrid()
  {
  if(gridsize>0)
    {
    ST::string outest = pathcurrent;
    outest = outest.substr(0,outest.length()-4) + "_randomdesign.raw";

    ofstream outresZ(outest.strtochar());
    Z_grid.prettyPrint(outresZ);
    outresZ.close();

    outest = pathcurrent;
    outest = outest.substr(0,outest.length()-4) + "_fixeddesign.raw";

    ofstream outresX(outest.strtochar());
    X_grid.prettyPrint(outresX);
    outresX.close();
    }
  }

void spline_basis_surf::outoptionsreml()
  {
  optionsp->out("OPTIONS FOR P-SPLINE TERM:: " + title + "\n",true);
  optionsp->out("\n");

  ST::string typestr;
  if(type==mrflinear)
    {
    optionsp->out("  Prior: 2 dimensional first order random walk\n");
    optionsp->out("         (Kronecker sum of two first order random walks)\n");
    }
  else if(type==mrfquadratic8)
    {
    optionsp->out("  Prior: 2 dimensional second order random walk\n");
    optionsp->out("         (Kronecker sum of two second order random walks)\n");
    }
  else if(type==mrfquadratic12)
    {
    optionsp->out("  Prior: 2 dimensional second order random walk\n");
    optionsp->out("         (Approximation to the biharmonic differential operator)\n");
    }

  optionsp->out("  Number of knots: " + ST::inttostring(nrknots) + "\n" );
  optionsp->out("  Degree of Splines: " + ST::inttostring(degree) + "\n" );
  optionsp->out("  Starting value for lambda: " + ST::doubletostring(startlambda,6) + "\n" );
  optionsp->out("\n");
  }

} // end: namespace MCMC



//---------------------------------------------------------------------------
#if !defined(__BUILDING_GNU)
#pragma package(smart_init)
#endif




