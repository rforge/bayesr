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



#include "kriging.h"

namespace MCMC
{

FULLCOND_kriging::FULLCOND_kriging(MCMCoptions * o, const datamatrix & v1,
               const datamatrix & v2, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu,const fieldtype & ft, const ST::string & ti,
               const ST::string & fp,const ST::string & pres, const double & l,
               const double & sl, const bool & catsp, const unsigned & gsx,
               const unsigned & gsy)
  : FULLCOND_nonp_basis(o,ti)
  {
  catspecific = catsp;

  mapexisting=false;
  plotstyle=noplot;

  varcoeff=false;
  onedim=false;

  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;
  samplepath=fp;

  nu=n;
  maxdist=maxd;

  p=pval;
  q=qval;
  maxsteps=maxst;

  type = ft;

  lambda = l;
  startlambda = sl;

  xorig=v1;
  yorig=v2;

  make_index(v1,v2);
  make_xy_values(v1,v2);

  full=fu;
  if(full)
    {
    nrknots=nrdiffobs;
    }
  else
    {
    nrknots = nrk;
    }

  xknots.clear();
  yknots.clear();
  if(knotdata.cols()>1)
    {
    spacefill=false;
    unsigned i;
    nrknots=knotdata.rows();
    for(i=0; i<nrknots; i++)
      {
      xknots.push_back(knotdata(i,0));
      yknots.push_back(knotdata(i,1));
      }
    }
  else
    {
    spacefill=true;
    compute_knots(xvalues,yvalues);
    }

  nrpar = nrknots;
  dimX = 0;
  dimZ = nrknots;

  // berechne rho
  unsigned i,j;
  rho=0;
  double norm2;
  for(i=0; i<xvalues.size(); i++)
    {
    for(j=0; j<xvalues.size(); j++)
      {
      norm2=(xvalues[i]-xvalues[j])*(xvalues[i]-xvalues[j])+(yvalues[i]-yvalues[j])*(yvalues[i]-yvalues[j]);
      if(norm2>rho)
        {
        rho=norm2;
        }
      }
    }
  rho=sqrt(rho)/maxdist;

  gridsizex = gsx;
  gridsizey = gsy;
  gridsize = gsx*gsy;
  if(gridsize>0)
    {
    X_grid = datamatrix(gridsize,dimX,1.0);
    Z_grid = datamatrix(gridsize,dimZ,0.0);
    }
  make_xy_values_grid(v1,v2);
  }

FULLCOND_kriging::FULLCOND_kriging(MCMCoptions * o,
               const datamatrix & intact, const datamatrix & v1,
               const datamatrix & v2, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu,const fieldtype & ft, const ST::string & ti,
               const ST::string & fp,const ST::string & pres, const double & l,
               const double & sl, const bool & catsp)
  : FULLCOND_nonp_basis(o,ti)
  {
  catspecific = catsp;

  mapexisting=false;
  plotstyle=noplot;

  varcoeff=true;
  onedim=false;

  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;
  samplepath=fp;

  nu=n;
  maxdist=maxd;

  p=pval;
  q=qval;
  maxsteps=maxst;

  type = ft;

  lambda = l;
  startlambda = sl;

  xorig=v1;
  yorig=v2;

  make_index(v1,v2);
  make_xy_values(v1,v2);

  full=fu;
  if(full)
    {
    nrknots=nrdiffobs;
    }
  else
    {
    nrknots = nrk;
    }


  xknots.clear();
  yknots.clear();
  if(knotdata.cols()>1)
    {
    spacefill=false;
    unsigned i;
    nrknots=knotdata.rows();
    for(i=0; i<nrknots; i++)
      {
      xknots.push_back(knotdata(i,0));
      yknots.push_back(knotdata(i,1));
      }
    }
  else
    {
    spacefill=true;
    compute_knots(xvalues,yvalues);
    }

  nrpar = nrknots;
  dimX = 0;
  dimZ = nrknots;

  Z_VCM = datamatrix(intact.rows(),dimZ,0.0);
  data_forfixed = intact;

  // berechne rho
  unsigned i,j;
  rho=0;
  double norm2;
  for(i=0; i<xvalues.size(); i++)
    {
    for(j=0; j<xvalues.size(); j++)
      {
      norm2=(xvalues[i]-xvalues[j])*(xvalues[i]-xvalues[j])+(yvalues[i]-yvalues[j])*(yvalues[i]-yvalues[j]);
      if(norm2>rho)
        {
        rho=norm2;
        }
      }
    }
  rho=sqrt(rho)/maxdist;
  }

FULLCOND_kriging::FULLCOND_kriging(MCMCoptions * o, const datamatrix & region,
               const MAP::map & mp, const ST::string & mn, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp, const unsigned & gsx,
               const unsigned & gsy)
  : FULLCOND_nonp_basis(o,ti)
  {
  catspecific = catsp;

  m = mp;
  mapexisting = true;
  mapname = mn;
  if(mp.polygones_existing() == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

  varcoeff=false;
  onedim=false;

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
  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;
  samplepath=fp;

  nu=n;
  maxdist=maxd;

  p=pval;
  q=qval;
  maxsteps=maxst;

  type = ft;

  lambda = l;
  startlambda = sl;

  xorig=v1;
  yorig=v2;

  make_index(v1,v2);
  make_xy_values(v1,v2);

  full=fu;
  if(full)
    {
    nrknots=nrdiffobs;
    }
  else
    {
    nrknots = nrk;
    }

  xknots.clear();
  yknots.clear();
  if(knotdata.cols()>1)
    {
    spacefill=false;
    unsigned i;
    nrknots=knotdata.rows();
    for(i=0; i<nrknots; i++)
      {
      xknots.push_back(knotdata(i,0));
      yknots.push_back(knotdata(i,1));
      }
    }
  else
    {
    spacefill=true;
    compute_knots(xvalues,yvalues);
    }

  nrpar = nrknots;
  dimX = 0;
  dimZ = nrknots;

  // berechne rho
  unsigned i,j;
  rho=0;
  double norm2;
  for(i=0; i<xvalues.size(); i++)
    {
    for(j=0; j<xvalues.size(); j++)
      {
      norm2=(xvalues[i]-xvalues[j])*(xvalues[i]-xvalues[j])+(yvalues[i]-yvalues[j])*(yvalues[i]-yvalues[j]);
      if(norm2>rho)
        {
        rho=norm2;
        }
      }
    }
  rho=sqrt(rho)/maxdist;

  gridsizex = gsx;
  gridsizey = gsy;
  gridsize = gsx*gsy;
  if(gridsize>0)
    {
    X_grid = datamatrix(gridsize,dimX,1.0);
    Z_grid = datamatrix(gridsize,dimZ,0.0);
    }
  make_xy_values_grid(v1,v2);
  }

  // Constructor 3: geokriging (varcoeff)

FULLCOND_kriging::FULLCOND_kriging(MCMCoptions * o, const datamatrix & intact,
               const datamatrix & region, const MAP::map & mp,
               const ST::string & mn, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp)
  : FULLCOND_nonp_basis(o,ti)
  {
  catspecific = catsp;
  gridsize=0;

  m = mp;
  mapexisting = true;
  mapname = mn;
  if(mp.polygones_existing() == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

  varcoeff=true;
  data_forfixed = intact;
  onedim=false;

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
  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;
  samplepath=fp;

  nu=n;
  maxdist=maxd;

  p=pval;
  q=qval;
  maxsteps=maxst;

  type = ft;

  lambda = l;
  startlambda = sl;

  xorig=v1;
  yorig=v2;

  make_index(v1,v2);
  make_xy_values(v1,v2);

  full=fu;
  if(full)
    {
    nrknots=nrdiffobs;
    }
  else
    {
    nrknots = nrk;
    }

  xknots.clear();
  yknots.clear();
  if(knotdata.cols()>1)
    {
    spacefill=false;
    unsigned i;
    nrknots=knotdata.rows();
    for(i=0; i<nrknots; i++)
      {
      xknots.push_back(knotdata(i,0));
      yknots.push_back(knotdata(i,1));
      }
    }
  else
    {
    spacefill=true;
    compute_knots(xvalues,yvalues);
    }

  nrpar = nrknots;
  dimX = 0;
  dimZ = nrknots;
  Z_VCM = datamatrix(region.rows(),dimZ,0.0);

  // berechne rho
  unsigned i,j;
  rho=0;
  double norm2;
  for(i=0; i<xvalues.size(); i++)
    {
    for(j=0; j<xvalues.size(); j++)
      {
      norm2=(xvalues[i]-xvalues[j])*(xvalues[i]-xvalues[j])+(yvalues[i]-yvalues[j])*(yvalues[i]-yvalues[j]);
      if(norm2>rho)
        {
        rho=norm2;
        }
      }
    }
  rho=sqrt(rho)/maxdist;
  }

  // Constructor 4: 1-dimensional kriging

  FULLCOND_kriging::FULLCOND_kriging(MCMCoptions * o, const datamatrix & v,
               const double & n, const double & maxd,
               const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp)
  : FULLCOND_nonp_basis(o,ti)
  {
  catspecific = catsp;

  mapexisting=false;
  plotstyle=plotnonp;

  varcoeff=false;
  onedim=true;

  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;
  samplepath=fp;

  nu=n;
  maxdist=maxd;

  type = ft;

  lambda = l;
  startlambda = sl;

  xorig=v;

  make_index(v);

  xknots = xvalues;

  // compute incidence vector (contains for each observation the position in the knots-vector)

  unsigned i,j;
  incidence = vector<unsigned>(v.rows(),0);
  for(i=0; i<v.rows(); i++)
    {
    j=0;
    while(xknots[j]!=v(i,0))
      {
      j++;
      }
    incidence[i]=j;
    }

  nrknots = nrdiffobs;
  nrpar = nrknots;
  dimX = 0;
  dimZ = nrknots;

  // berechne rho
  rho=0;
  double absv;
  for(i=0; i<xvalues.size(); i++)
    {
    for(j=0; j<xvalues.size(); j++)
      {
      absv=fabs(xvalues[i]-xvalues[j]);
      if(absv>rho)
        {
        rho=absv;
        }
      }
    }
  rho=rho/maxdist;
  }


FULLCOND_kriging::FULLCOND_kriging(const FULLCOND_kriging & kr)
  : FULLCOND_nonp_basis(FULLCOND_nonp_basis(kr))
  {
  nrknots=kr.nrknots;
  nu=kr.nu;
  full=kr.full;
  xknots=kr.xknots;
  yknots=kr.yknots;
  xvalues=kr.xvalues;
  yvalues=kr.yvalues;
  xorig=kr.xorig;
  yorig=kr.yorig;
  index2=kr.index2;
  freq=kr.freq;
  freqoutput=kr.freqoutput;
  nrdiffobs=kr.nrdiffobs;
  rho=kr.rho;
  p=kr.p;
  q=kr.q;
  maxsteps=kr.maxsteps;
  spacefill=kr.spacefill;
  m = kr.m;
  mapexisting = kr.mapexisting;
  mapname = kr.mapname;
  regionnames = kr.regionnames;
  Z_VCM = kr.Z_VCM;
  X_VCM = kr.X_VCM;
  data_forfixed = kr.data_forfixed;
  onedim = kr.onedim;
  incidence = kr.incidence;
  gridsize = kr.gridsize;
  gridsizex = kr.gridsizex;
  gridsizey = kr.gridsizey;
  xvaluesgrid=kr.xvaluesgrid;
  yvaluesgrid=kr.yvaluesgrid;
  effectvaluesxgrid=kr.effectvaluesxgrid;
  effectvaluesygrid=kr.effectvaluesygrid;
  Z_grid = kr.Z_grid;
  X_grid = kr.X_grid;
  }

const FULLCOND_kriging & FULLCOND_kriging::operator=(const FULLCOND_kriging & kr)
  {
  if(&kr==this)
    return *this;
  FULLCOND_nonp_basis::operator=(FULLCOND_nonp_basis(kr));

  nrknots=kr.nrknots;
  nu=kr.nu;
  full=kr.full;
  xknots=kr.xknots;
  yknots=kr.yknots;
  xvalues=kr.xvalues;
  yvalues=kr.yvalues;
  xorig=kr.xorig;
  yorig=kr.yorig;
  index2=kr.index2;
  freq=kr.freq;
  freqoutput=kr.freqoutput;
  nrdiffobs=kr.nrdiffobs;
  rho=kr.rho;
  p=kr.p;
  q=kr.q;
  maxsteps=kr.maxsteps;
  spacefill=kr.spacefill;
  m = kr.m;
  mapexisting = kr.mapexisting;
  mapname = kr.mapname;
  regionnames = kr.regionnames;
  Z_VCM = kr.Z_VCM;
  X_VCM = kr.X_VCM;
  data_forfixed = kr.data_forfixed;
  onedim = kr.onedim;
  incidence = kr.incidence;
  gridsize = kr.gridsize;
  gridsizex = kr.gridsizex;
  gridsizey = kr.gridsizey;
  xvaluesgrid=kr.xvaluesgrid;
  yvaluesgrid=kr.yvaluesgrid;
  effectvaluesxgrid=kr.effectvaluesxgrid;
  effectvaluesygrid=kr.effectvaluesygrid;
  Z_grid = kr.Z_grid;
  X_grid = kr.X_grid;
  return *this;
  }


void FULLCOND_kriging::createreml(datamatrix & X,datamatrix & Z,
                  const unsigned & Xpos,const unsigned & Zpos)
  {
  unsigned i,j;
  double r;

  if(!onedim)
    {
    // Korrelationen zwischen Daten und Knoten berechnen
    for(i=0; i<Z.rows(); i++)
      {
      for(j=0; j<nrknots; j++)
        {
        r=sqrt((xorig(i,0)-xknots[j])*(xorig(i,0)-xknots[j]) +
                 (yorig(i,0)-yknots[j])*(yorig(i,0)-yknots[j]))/rho;
        if(nu==0.5)
          {
          Z(i,Zpos+j)=exp(-r);
          }
        else if(nu==1.5)
          {
          Z(i,Zpos+j)=exp(-r)*(1+r);
          }
        else if(nu==2.5)
          {
          Z(i,Zpos+j)=exp(-r)*(1+r+r*r/3);
          }
        else if(nu==3.5)
          {
          Z(i,Zpos+j)=exp(-r)*(1+r+2*r*r/5+r*r*r/15);
          }
        }
      }
    if(gridsize>0)
      {
      for(i=0; i<gridsize; i++)
        {
        for(j=0; j<nrknots; j++)
          {
          r=sqrt((effectvaluesxgrid[i]-xknots[j])*(effectvaluesxgrid[i]-xknots[j]) +
                   (effectvaluesygrid[i]-yknots[j])*(effectvaluesygrid[i]-yknots[j]))/rho;
          if(nu==0.5)
            {
            Z_grid(i,j)=exp(-r);
            }
          else if(nu==1.5)
            {
            Z_grid(i,j)=exp(-r)*(1+r);
            }
          else if(nu==2.5)
            {
            Z_grid(i,j)=exp(-r)*(1+r+r*r/3);
            }
          else if(nu==3.5)
            {
            Z_grid(i,j)=exp(-r)*(1+r+2*r*r/5+r*r*r/15);
            }
          }
        }
      }

    // Korrelationen zwischen Knoten berechnen
    datamatrix cov(xknots.size(),xknots.size(),0);
    for(i=0; i<cov.rows(); i++)
      {
      for(j=0; j<cov.cols(); j++)
        {
        r=sqrt((xknots[i]-xknots[j])*(xknots[i]-xknots[j]) + (yknots[i]-yknots[j])*(yknots[i]-yknots[j]))/rho;
        if(nu==0.5)
          {
          cov(i,j)=exp(-r);
          }
        else if(nu==1.5)
          {
          cov(i,j)=exp(-r)*(1+r);
          }
        else if(nu==2.5)
          {
          cov(i,j)=exp(-r)*(1+r+r*r/3);
          }
        else if(nu==3.5)
          {
          cov(i,j)=exp(-r)*(1+r+2*r*r/5+r*r*r/15);
          }
        }
      }

    // Z modifizieren

    // Cholesky-Wurzel der Inversen Kovarianzmatrix

    cov=cov.inverse();
    cov=cov.root();

/*  // Eigenwertwurzel
  datamatrix vals(cov.rows(),1,0);
  bool test=eigen2(cov,vals);
  eigensort(vals,cov);
  for(i=0; i<vals.rows(); i++)
    {
    vals(i,0)=sqrt(vals(i,0));
    }
  cov.multdiagback(vals);*/

    if(!varcoeff)
      {
      Z.putColBlock(Zpos,Zpos+nrknots,Z.getColBlock(Zpos,Zpos+nrknots)*cov);
      if(gridsize>0)
        {
        Z_grid = Z_grid*cov;
        }
      }
    else
      {
      Z_VCM = Z.getColBlock(Zpos,Zpos+nrknots)*cov;
      Z.putColBlock(Zpos,Zpos+nrknots,multdiagfront(Z_VCM,data_forfixed));
      }
    }
  else
    {
    // Korrelationen zwischen Knoten berechnen
    datamatrix cov(xknots.size(),xknots.size(),0);
    for(i=0; i<cov.rows(); i++)
      {
      for(j=0; j<cov.cols(); j++)
        {
        r=fabs(xknots[i]-xknots[j])/rho;
        if(nu==0.5)
          {
          cov(i,j)=exp(-r);
          }
        else if(nu==1.5)
          {
          cov(i,j)=exp(-r)*(1+r);
          }
        else if(nu==2.5)
          {
          cov(i,j)=exp(-r)*(1+r+r*r/3);
          }
        else if(nu==3.5)
          {
          cov(i,j)=exp(-r)*(1+r+2*r*r/5+r*r*r/15);
          }
        }
      }
    cov=cov.root();
    for(i=0;i<Z.rows();i++)
      {
      for(j=0; j<xknots.size(); j++)
        {
        Z(i,Zpos+j)=cov(incidence[i],j);
        }
      }
    }
  }

void FULLCOND_kriging::make_index(const datamatrix & var1,const datamatrix & var2)
  {
  unsigned i,j,k,beg,end;
  unsigned nrobs = var1.rows();

  index = statmatrix<int>(nrobs,1);
  index.indexinit();
  var1.indexsort(index,0,nrobs-1,0,0);

  i = 0;
  j = 1;
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

  index2.push_back(index(0,0));
  for(i=1;i<index.rows();i++)
    index2.push_back(index(i,0)-index(i-1,0));

  }

void FULLCOND_kriging::make_xy_values(const datamatrix & var1,const datamatrix & var2)
  {
  unsigned i;
  vector<int>::iterator indexit = index2.begin();
  unsigned k = *indexit;
  vector<int>::iterator freqwork = freqoutput.begin();

  for(i=0;i<var1.rows();i++,indexit++,freqwork++,k+=*indexit)
    {
    if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
      {
      xvalues.push_back(var1(k,0));
      yvalues.push_back(var2(k,0));
      }
    }
  }

void FULLCOND_kriging::compute_knots(const vector<double> & xvals,
                                     const vector<double> & yvals)
  {
  if(nrknots>xvals.size())
    {
    errors.push_back("ERROR: More knots requested than different locations observed");
    }
  else if(nrknots==xvals.size())
    {
    xknots=xvals;
    yknots=yvals;
    }
  else
    {
    optionsp->out("\n");
    optionsp->out("\n");
    optionsp->out("Computing knots (this may take some time)\n",true);
    optionsp->out("\n");
    //Space-Filling-Algorithmus
    unsigned i,j,k;

    unsigned nrcand = nrdiffobs-nrknots; // Anzahl der Kandidatenpunkte

    vector<unsigned> knotindex(nrknots,0);
    vector<unsigned> candindex(nrcand,0);
    datamatrix distmat(nrcand, nrknots,0);     // Distanzmatrix
    datamatrix rowsums(nrcand,1,0);            // Zeilensummen der Distanzmatrix

    datamatrix rowsumswithouti(nrcand,1,0);    // Zeilensummen ohne den Beitrag von Knoten i
    double partialnewrow;   // Zeilensumme der neuen Zeile der Distanzmatrix, die sich ergibt, wenn Knoten i
                            // zum Kandidaten wird (ohne Beitrag des neuen Knotens)

    int swapindex;                        // Index des zu tauschenden Kandidaten
                                          // (Kandidat mit dem minimalen Kriterium beim Loop über den aktuellen Knoten)

    double covcritold=0; // Coverage-Kriterium (ändert sich nur nach Durchlauf über Knoten)
    double covcritnew=-1;   // Neues Coverage-Kriterium (ändert sich nur nach Durchlauf über Knoten)
    double covcritoldi;  // Coverage-Kriterium vorm Swappen von Knoten i (ändert sich nach jedem Swap)
    double covcritnewi;  // Neues Coverage-Kriterium beim Swappen von Knoten i (ändert sich nach jedem Swap)

    // Startdesign und Kandidaten zufällig bestimmen

    datamatrix u(nrdiffobs,1,0);
    statmatrix<int> ind(nrdiffobs,1,0);
    ind.indexinit();
    for(i=0; i<nrdiffobs; i++)
      {
      u(i,0) = uniform();
      }
    u.indexsort(ind,0,nrdiffobs-1,0,0);
    for(i=0; i<nrknots; i++)
      {
      knotindex[i] = ind(i,0);
      }
    for(i=nrknots, j=0; i<nrdiffobs; i++, j++)
      {
      candindex[j] = ind(i,0);
      }

    // Zentriere xvals und yvals
    double xmean=0;
    double ymean=0;
    for(i=0; i<nrdiffobs; i++)
      {
      xmean += xvals[i];
      ymean += yvals[i];
      }
    vector<double> xvalscentered(nrdiffobs);
    vector<double> yvalscentered(nrdiffobs);
    for(i=0; i<nrdiffobs; i++)
      {
      xvalscentered[i] = xvals[i]-xmean;
      yvalscentered[i] = yvals[i]-ymean;
      }


/*    datamatrix kandidaten(nrdiffobs,2,0);
    statmatrix<int> knotenindex(nrknots,1,0);
    for(i=0; i<nrdiffobs; i++)
      {
      kandidaten(i,0)=xvalscentered[i];
      kandidaten(i,1)=yvalscentered[i];
      }
    for(i=0; i<nrknots; i++)
      {
      knotenindex(i,0)=knotindex[i];
      }
    ofstream out1("c:\\temp\\kandidaten.raw");
    kandidaten.prettyPrint(out1);
    out1.close();
    ofstream out2("c:\\temp\\knotenindex.raw");
    knotenindex.prettyPrint(out2);
    out2.close();*/

    // Distanzmatrix mit Startdesign und Zeilensummen ausrechnen

    for(i=0; i<nrcand; i++)
      {
      for(j=0; j<nrknots; j++)
        {
        distmat(i,j) = pow(
               (xvalscentered[knotindex[j]]-xvalscentered[candindex[i]])*(xvalscentered[knotindex[j]]-xvalscentered[candindex[i]])
               + (yvalscentered[knotindex[j]]-yvalscentered[candindex[i]])*(yvalscentered[knotindex[j]]-yvalscentered[candindex[i]]),p/2);
        rowsums(i,0) += distmat(i,j);
        }
      }
    for(i=0; i<nrcand; i++)
      {
      covcritold += pow(rowsums(i,0),q/p);
      }
    covcritold = pow(covcritold,1/q);
    covcritoldi = covcritold;

    // Swappen
    unsigned steps=1;
    while(covcritnew<covcritold && steps<=maxsteps)
      {
      if(steps>1)
        {
        covcritold=covcritnew;
        }
      for(i=0; i<nrknots; i++)
        {
        swapindex = -1;
        for(k=0; k<nrcand; k++)
          {
          rowsumswithouti(k,0)=rowsums(k,0)-distmat(k,i);
          }
        partialnewrow=0;
        for(k=0; k<i; k++)
          {
          partialnewrow += pow(
               (xvalscentered[knotindex[k]]-xvalscentered[knotindex[i]])*(xvalscentered[knotindex[k]]-xvalscentered[knotindex[i]])
               + (yvalscentered[knotindex[k]]-yvalscentered[knotindex[i]])*(yvalscentered[knotindex[k]]-yvalscentered[knotindex[i]]),p/2);
          }
        for(k=i+1; k<nrknots; k++)
          {
          partialnewrow += pow(
               (xvalscentered[knotindex[k]]-xvalscentered[knotindex[i]])*(xvalscentered[knotindex[k]]-xvalscentered[knotindex[i]])
               + (yvalscentered[knotindex[k]]-yvalscentered[knotindex[i]])*(yvalscentered[knotindex[k]]-yvalscentered[knotindex[i]]),p/2);
          }

        for(j=0; j<nrcand; j++)
          {
          // Loop über die Kandidaten. Jeweil das Coverage-Kriterium bei Swap i vs. j ausrechnen
          // und gegebenenfalls swapindex umsetzen.
          covcritnewi = 0;
          // Loopen über Zeilensummen ohne j
          for(k=0; k<j; k++)
            {
            covcritnewi += pow(rowsumswithouti(k,0) +
               pow((xvalscentered[candindex[k]]-xvalscentered[candindex[j]])*(xvalscentered[candindex[k]]-xvalscentered[candindex[j]])
               + (yvalscentered[candindex[k]]-yvalscentered[candindex[j]])*(yvalscentered[candindex[k]]-yvalscentered[candindex[j]]),p/2),q/p);
            }
          for(k=j+1; k<nrcand; k++)
            {
            covcritnewi += pow(rowsumswithouti(k,0) +
               pow((xvalscentered[candindex[k]]-xvalscentered[candindex[j]])*(xvalscentered[candindex[k]]-xvalscentered[candindex[j]])
               + (yvalscentered[candindex[k]]-yvalscentered[candindex[j]])*(yvalscentered[candindex[k]]-yvalscentered[candindex[j]]),p/2),q/p);
            }
          // Neue Zeile
          covcritnewi += pow(partialnewrow + distmat(j,i),q/p);
          covcritnewi = pow(covcritnewi,1/q);

          if(covcritnewi<covcritoldi)
            {
            swapindex=j;
            covcritnew=covcritnewi;
            covcritoldi=covcritnew;
            }
          }
        if(swapindex!=-1)
          {
          // Swappen
          unsigned h = knotindex[i];
          knotindex[i] = candindex[swapindex];
          candindex[swapindex] = h;

//          unsigned test1=h;
//          unsigned test2=knotindex[i];

          // Distanzmatrix und Zeilensummen updaten
          // Neue Spalte ausrechnen (Index i)
          for(j=0; j<nrcand; j++)
            {
            distmat(j,i) = pow(
               (xvalscentered[knotindex[i]]-xvalscentered[candindex[j]])*(xvalscentered[knotindex[i]]-xvalscentered[candindex[j]])
               + (yvalscentered[knotindex[i]]-yvalscentered[candindex[j]])*(yvalscentered[knotindex[i]]-yvalscentered[candindex[j]]),p/2);
            }

          // Neue Zeile ausrechnen (Index swapindex)
          for(j=0; j<nrknots; j++)
            {
            distmat(swapindex,j) = pow(
               (xvalscentered[knotindex[j]]-xvalscentered[candindex[swapindex]])*(xvalscentered[knotindex[j]]-xvalscentered[candindex[swapindex]])
               + (yvalscentered[knotindex[j]]-yvalscentered[candindex[swapindex]])*(yvalscentered[knotindex[j]]-yvalscentered[candindex[swapindex]]),p/2);
            }

        // Zeilensummen
          for(k=0; k<nrcand; k++)
            {
            rowsums(k,0)=0;
            for(j=0; j<nrknots; j++)
              {
              rowsums(k,0) += distmat(k,j);
              }
            }
          }
        }
      steps++;
      }


    // Knoten speichern
    for(i=0; i<nrknots; i++)
      {
      xknots.push_back(xvals[knotindex[i]]);
      yknots.push_back(yvals[knotindex[i]]);
      }
    }
  }

void FULLCOND_kriging::outoptionsreml()
  {
  optionsp->out("OPTIONS FOR KRIGING TERM:: " + title + "\n",true);
  optionsp->out("\n");
  optionsp->out("  Correlation function: Matern\n");
  optionsp->out("    Parameter nu: " + ST::doubletostring(nu,2) + "\n");
  optionsp->out("    Parameter rho: " + ST::doubletostring(rho,4) + "\n");
  if(!full)
    {
    optionsp->out("  Number of knots: " + ST::inttostring(nrknots) + "\n");
    optionsp->out("    Corresponds to " + ST::doubletostring(100*(double)nrknots/(double)nrdiffobs,3) + "% of the " + ST::inttostring(nrdiffobs) + " different observation points\n");
    }
  if(nrknots<nrdiffobs && spacefill)
    {
    optionsp->out("  Options for the space-filling algorithm:\n");
    optionsp->out("    Maximum number of iteration steps: " + ST::inttostring(maxsteps) + "\n");
    optionsp->out("    p: " + ST::inttostring(p) + "\n");
    optionsp->out("    q: " + ST::inttostring(q) + "\n");
    }
  optionsp->out("  Starting value for lambda: " + ST::doubletostring(startlambda,6) + "\n" );
  }

double FULLCOND_kriging::outresultsreml(datamatrix & X,datamatrix & Z,
                                        datamatrix & betareml,
                                        datamatrix & betacov,
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

  unsigned nr = xvalues.size();
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

  vector<ST::string> regionnameshelp=regionnames;

  if(!varcoeff)
    {
    for(i=0,j=0;i<Z.rows();i++,indexit++,freqwork++,k+=*indexit)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        betamean(j,0) = (Z.getBlock(k,Zpos,k+1,Zpos+nrpar)*betareml.getBlock(betaZpos,0,betaZpos+nrpar,1))(0,0);
        betastd(j,0) = sqrt((Z.getBlock(k,Zpos,k+1,Zpos+nrpar)*
                 betacov.getBlock(betaZpos,betaZpos,betaZpos+nrpar,betaZpos+nrpar)*
                 Z.getBlock(k,Zpos,k+1,Zpos+nrpar).transposed())(0,0));
        if(mapexisting)
         {
          regionnameshelp[j]=regionnames[k];
          }
        j++;
        }
      }
    }
  else
    {
    for(i=0,j=0;i<Z.rows();i++,indexit++,freqwork++,k+=*indexit)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        betamean(j,0) = (Z_VCM.getRow(k)*betareml.getBlock(betaZpos,0,betaZpos+nrpar,1))(0,0);
        betastd(j,0) = sqrt((Z_VCM.getRow(k)*
                 betacov.getBlock(betaZpos,betaZpos,betaZpos+nrpar,betaZpos+nrpar)*
                 Z_VCM.getRow(k).transposed())(0,0));
        if(mapexisting)
         {
          regionnameshelp[j]=regionnames[k];
          }
        j++;
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

  FULLCOND::outresults();

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

  if(!onedim)
    {
    optionsp->out("\n");
    optionsp->out("  Knots are stored in file\n");
    optionsp->out("  " + pathcurrent.substr(0,pathcurrent.length()-4)+"_knots.raw" + "\n");

    ofstream outknots((pathcurrent.substr(0,pathcurrent.length()-4)+"_knots.raw").strtochar());
    if(mapexisting)
      {
      outknots << "x" << "  " << "y" << endl;
      }
    else
      {
      outknots << datanames[1] << "  " << datanames[0] << endl;
      }
    for(i=0; i<nrknots; i++)
      {
      outknots << xknots[i] << "   " << yknots[i] << endl;
      }
    outknots.close();
    }

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
  else if(onedim)
    {
    #if defined(JAVA_OUTPUT_WINDOW)
    optionsp->out("  Postscript file is stored in file\n");
    ST::string psfile = outest.substr(0,outest.length()-4) + ".ps";
    optionsp->out("  " + psfile + "\n");
    optionsp->out("\n");
    optionsp->out("  Results may be visualized using method 'plotnonp'\n");
    optionsp->out("  Type for example: objectname.plotnonp " + ST::inttostring(plotpos) + "\n");
    #else
    optionsp->out("  Results may be visualized using the R function 'plotnonp'\n");
    ST::string doublebackslash = "/";
    ST::string spluspath = outest.insert_string_char('\\',doublebackslash);
    optionsp->out("  Type for example:\n");
    optionsp->out("  plotnonp(\"" + spluspath + "\")");
    optionsp->out("\n");
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
  else if(onedim)
    {
    outres << datanames[0] << "   ";
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
  vector<double>::iterator workxvalues = xvalues.begin();
  vector<double>::iterator workyvalues = yvalues.begin();

  for(i=0;i<xvalues.size();i++,workmean++,workstd++,
                     workbetaqu_l1_lower_p++,workbetaqu_l2_lower_p++,
                     workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                     workxvalues++,workyvalues++)
    {
    outres << (i+1) << "   ";
    if(mapexisting)
      outres << regionnameshelp[i] << "   ";
    outres << *workxvalues << "   ";
    if(!onedim)
      outres << *workyvalues << "   ";
    outres << *workmean << "   ";
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
    nr = gridsize;

    betamean=datamatrix(nr,1,0);
    datamatrix betastd=datamatrix(nr,1,0);
    betaqu_l1_lower=datamatrix(nr,1,0);
    betaqu_l1_upper=datamatrix(nr,1,0);
    betaqu_l2_lower=datamatrix(nr,1,0);
    betaqu_l2_upper=datamatrix(nr,1,0);

    betamean = Z_grid*betareml.getBlock(betaZpos,0,betaZpos+nrpar,1);
    for(i=0; i<gridsize; i++)
      {
      betamean(i,0)=betamean(i,0)-mean;
      betastd(i,0) = sqrt((Z_grid.getRow(i)*
               betacov.getBlock(betaZpos,betaZpos,betaZpos+nrpar,betaZpos+nrpar)*
               Z_grid.getRow(i).transposed())(0,0));
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

  return mean;
  }

void FULLCOND_kriging::init_names(const vector<ST::string> & na)
  {
  FULLCOND::init_names(na);
  ST::string underscore = "\\_";
  if(!varcoeff)
    {
    if(mapexisting || onedim)
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
      ST::string helpname1 = na[1].insert_string_char('_',underscore);
      ST::string helpname2 = na[2].insert_string_char('_',underscore);
      term_symbolic = "f_{" +  helpname0 + "," + helpname1 + "}(" + helpname0 + "," + helpname1 + ")" + " \\cdot " + helpname2;
      }
    }


  priorassumptions.push_back("$" + term_symbolic + "$");
  priorassumptions.push_back("Stationary Gaussian Random Field");

  priorassumptions.push_back("Correlation function: Matern\n");
  priorassumptions.push_back("Parameter nu: " + ST::doubletostring(nu,2));
  priorassumptions.push_back("Parameter rho: " + ST::doubletostring(rho,4));
  if(!full)
    {
    priorassumptions.push_back("Number of knots: " + ST::inttostring(nrknots));
    priorassumptions.push_back("Corresponds to " + ST::doubletostring(100*(double)nrknots/(double)nrdiffobs,3) + "\\% of the " + ST::inttostring(nrdiffobs) + " different observation points");
    }
  if(nrknots<nrdiffobs && spacefill)
    {
    priorassumptions.push_back("Options for the space-filling algorithm:");
    priorassumptions.push_back("Maximum number of iteration steps: " + ST::inttostring(maxsteps));
    priorassumptions.push_back("p: " + ST::inttostring(p));
    priorassumptions.push_back("q: " + ST::inttostring(q));
    }
  }

ST::string FULLCOND_kriging::getinfo(void)
  {
  if(mapexisting)
    return mapname;
  else
    return title;
  }

void FULLCOND_kriging::make_index(const datamatrix & moddata)
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
  for(unsigned i=1;i<moddata.rows();i++)
    index2.push_back(index(i,0)-index(i-1,0));

  vector<int>::iterator freqwork = freqoutput.begin();
  workindex = index.getV();

  xvalues = vector<double>(nrdiffobs,0);
  for(j=0;j<moddata.rows();j++,freqwork++,workindex++)
    if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
      xvalues[*freqwork] = moddata(*workindex,0);
  }

void FULLCOND_kriging::make_xy_values_grid(const datamatrix & var1,const datamatrix & var2)
  {
  int i,j;

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

  effectvaluesxgrid = vector<double>(gridsize);
  effectvaluesygrid = vector<double>(gridsize);

  for(i=0;i<gridsizex;i++)
    {
    for(j=0;j<gridsizey;j++)
      {
      effectvaluesxgrid[i*gridsizey + j] = xvaluesgrid(i,0);
      effectvaluesygrid[i*gridsizey + j] = yvaluesgrid(j,0);
      }
    }
  }

} // end: namespace MCMC

//---------------------------------------------------------------------------
#if !defined(__BUILDING_GNU)
#pragma package(smart_init)
#endif


