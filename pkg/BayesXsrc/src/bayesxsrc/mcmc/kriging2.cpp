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



#include "kriging2.h"

namespace MCMC
{

FULLCOND_kriging2::FULLCOND_kriging2(MCMCoptions * o, DISTRIBUTION * dp, FULLCOND_const * fcc,
               const datamatrix & v1, const datamatrix & v2, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu,const fieldtype & ft, const ST::string & ti,
               const ST::string & fp,const ST::string & pres, const double & l,
               const double & sl, const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,1,MCMC::equidistant,-1,fp,pres,false,0.0,0.0,0.0,0.0,c)
  {

  mapexisting=false;
  plotstyle=noplot;

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

  X = datamatrix(likep->get_nrobs(),xknots.size(),0);
  create();

//

  }

FULLCOND_kriging2::FULLCOND_kriging2(MCMCoptions * o, DISTRIBUTION * dp, FULLCOND_const * fcc,
               const datamatrix & region,
               const MAP::map & mp, const ST::string & mn, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,1,MCMC::equidistant,-1,fp,pres,false,0.0,0.0,0.0,0.0,c)
  {

  utype = gaussian;

  m = mp;
  mapexisting = true;
  mapname = mn;
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

  setbeta(nrpar,1,0);

  prec_env = envmatdouble(0.0,xknots.size(),xknots.size()-1);
  mu = datamatrix(likep->get_nrobs(),1,0);
  muy = datamatrix(nrpar,1,0);
  betahelp = muy;
  standnormal = datamatrix(nrpar,1,0);

// fchelp initialisieren
  ST::string pnt = samplepath.substr(0,samplepath.length()-4)+"_fchelp.raw";
  if(gridsize < 0)
    {
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"help",nrdiffobs,1,pnt);
    splinehelp = datamatrix(likep->get_nrobs(),1,0);
    }
  else
    {
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"help",gridsize,1,pnt);
    splinehelp = datamatrix(gridsize,1,0);
    }
  fchelp.setflags(MCMC::norelchange | MCMC::nooutput);
  fchelp.set_transform(transform);

  X = datamatrix(likep->get_nrobs(),xknots.size(),0);
  create();

  identifiable = true;

  betaold = datamatrix(nrpar,1,0);

  }


FULLCOND_kriging2::FULLCOND_kriging2(MCMCoptions * o, DISTRIBUTION * dp, FULLCOND_const * fcc,
               const datamatrix & region,
               const MAP::map & mp, const ST::string & mn, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & mode, const unsigned & upW,
               const bool & updatetau, const double & fstart, const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,1,MCMC::equidistant,-1,fp,pres,false,0.0,0.0,0.0,0.0,c)
  {

  updateW = upW;
  f = fstart;

  if(mode)
    {
    if(updatetau)
      utype = hyperblockmode;
    else
      utype = iwlsmode;
    }
  else
    {
    if(updatetau)
      utype = hyperblock;
    else
      utype = iwls;
    }

  m = mp;
  mapexisting = true;
  mapname = mn;
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

  setbeta(nrpar,1,0.0);

  prec_env = envmatdouble(0.0,xknots.size(),xknots.size()-1);
  mu = datamatrix(likep->get_nrobs(),1,0);
  muy = datamatrix(nrpar,1,0);
  betahelp = muy;
  standnormal = datamatrix(nrpar,1,0);
  W = datamatrix(likep->get_nrobs(),1,0);

// fchelp initialisieren
  ST::string pnt = samplepath.substr(0,samplepath.length()-4)+"_fchelp.raw";
  if(gridsize < 0)
    {
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"help",nrdiffobs,1,pnt);
    splinehelp = datamatrix(likep->get_nrobs(),1,0);
    }
  else
    {
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"help",gridsize,1,pnt);
    splinehelp = datamatrix(gridsize,1,0);
    }
  fchelp.setflags(MCMC::norelchange | MCMC::nooutput);
  fchelp.set_transform(transform);

  X = datamatrix(likep->get_nrobs(),xknots.size(),0);
  create();

  identifiable = true;

  betaold = datamatrix(nrpar,1,0);
  betaprop = datamatrix(nrpar,1,0);

  }


FULLCOND_kriging2::FULLCOND_kriging2(const FULLCOND_kriging2 & kr)
  : spline_basis(spline_basis(kr))
  {
  X=kr.X;
  updateW=kr.updateW;
  utype=kr.utype;
  nu=kr.nu;
  full=kr.full;
  xknots=kr.xknots;
  yknots=kr.yknots;
  xvalues=kr.xvalues;
  yvalues=kr.yvalues;
  xorig=kr.xorig;
  yorig=kr.yorig;
  rho=kr.rho;
  p=kr.p;
  q=kr.q;
  maxsteps=kr.maxsteps;
  spacefill=kr.spacefill;
  m = kr.m;
  mapexisting = kr.mapexisting;
  mapname = kr.mapname;
  regionnames = kr.regionnames;
  }

const FULLCOND_kriging2 & FULLCOND_kriging2::operator=(const FULLCOND_kriging2 & kr)
  {
  if(&kr==this)
    return *this;
  FULLCOND_nonp_basis::operator=(FULLCOND_nonp_basis(kr));

  X=kr.X;
  updateW=kr.updateW;
  utype=kr.utype;
  nu=kr.nu;
  full=kr.full;
  xknots=kr.xknots;
  yknots=kr.yknots;
  xvalues=kr.xvalues;
  yvalues=kr.yvalues;
  xorig=kr.xorig;
  yorig=kr.yorig;
  rho=kr.rho;
  p=kr.p;
  q=kr.q;
  maxsteps=kr.maxsteps;
  spacefill=kr.spacefill;
  m = kr.m;
  mapexisting = kr.mapexisting;
  mapname = kr.mapname;
  regionnames = kr.regionnames;
  return *this;
  }


void FULLCOND_kriging2::create()
  {

  unsigned Zpos = 0;

  datamatrix K(xknots.size(),xknots.size(),0);

  unsigned i,j;

  // Korrelationen zwischen Daten und Knoten berechnen
  double r;
  for(i=0; i<X.rows(); i++)
    {
    for(j=0; j<nrknots; j++)
      {
      r=sqrt((xorig(i,0)-xknots[j])*(xorig(i,0)-xknots[j]) +
             (yorig(i,0)-yknots[j])*(yorig(i,0)-yknots[j]))/rho;
      if(nu==0.5)
        {
        X(i,Zpos+j)=exp(-r);
        }
      else if(nu==1.5)
        {
        X(i,Zpos+j)=exp(-r)*(1+r);
        }
      else if(nu==2.5)
        {
        X(i,Zpos+j)=exp(-r)*(1+r+r*r/3);
        }
      else if(nu==3.5)
        {
        X(i,Zpos+j)=exp(-r)*(1+r+2*r*r/5+r*r*r/15);
        }
      }
    }

  // Korrelationen zwischen Knoten berechnen
  for(i=0; i<K.rows(); i++)
    {
    for(j=0; j<K.cols(); j++)
      {
      r=sqrt((xknots[i]-xknots[j])*(xknots[i]-xknots[j]) + (yknots[i]-yknots[j])*(yknots[i]-yknots[j]))/rho;
      if(nu==0.5)
        {
        K(i,j)=exp(-r);
        }
      else if(nu==1.5)
        {
        K(i,j)=exp(-r)*(1+r);
        }
      else if(nu==2.5)
        {
        K(i,j)=exp(-r)*(1+r+r*r/3);
        }
      else if(nu==3.5)
        {
        K(i,j)=exp(-r)*(1+r+2*r*r/5+r*r*r/15);
        }
      }
    }
/*
  K=K.inverse();
  K=K.root();
  X=X*K;
  Kenv = envmatdouble(1.0,nrpar);
*/
  Kenv = envmatdouble(K);
  rankK = K.rows();

  betaweight = datamatrix(nrpar,1,0.0);
  for(i=0;i<X.cols();i++)
    {
    double mean = X.mean(i);
    betaweight(i,0) = mean;
    for(j=0;j<X.rows();j++)
      X(j,i) -= mean;
    }

  }

void FULLCOND_kriging2::make_index(const datamatrix & var1,const datamatrix & var2)
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

void FULLCOND_kriging2::make_xy_values(const datamatrix & var1,const datamatrix & var2)
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

void FULLCOND_kriging2::compute_knots(const vector<double> & xvals,
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


void FULLCOND_kriging2::update()
  {

  unsigned i;

  if(lambdaconst == true)
    sigma2 = likep->get_scale(column)/lambda;

  if(optionsp->get_nriter()==1)       // posterior mode Schätzung übernehmen
    betaold.assign(beta);

  if(utype == gaussian)
    update_gaussian();
  else if(utype == iwls)
    update_iwls();
  else if(utype == iwlsmode)
    update_iwlsmode();

  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1) % optionsp->get_step() == 0) )
    {
    double * fchelpbetap = fchelp.getbetapointer();

    vector<int>::iterator freqwork = freqoutput.begin();
    int * workindex = index.getV();
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        *fchelpbetap = spline(*workindex,0);
        fchelpbetap++;
        }
      }

    }

  fchelp.update();
  FULLCOND::update();

  }


void FULLCOND_kriging2::update_gaussian()
  {

  if (optionsp->get_nriter()==1)
    XX_env = envmatdouble(X.sscp());    // XWX

  if(lambdaconst == true)
    sigma2 = likep->get_scale(column)/lambda;

  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);

  unsigned i;

  likep->substr_linearpred(spline);

  if(changingweight)  // für t-link
    {
    datamatrix help = X.transposed();
    help.multdiagback(likep->get_weight());
    help = help*X;
    XX_env = envmatdouble(help);
    }

  double scaleinv = 1.0/likep->get_scale(column);

  prec_env.addto(XX_env,Kenv,scaleinv,1.0/sigma2);

  double * work = standnormal.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  likep->compute_respminuslinpred(mu,column);   // nicht ändern wegen multgaussian

//  compute_XWtildey(likep->get_weight(),scaleinv);
  double * weightp = likep->get_weightp();
  for(i=0;i<mu.rows();i++,weightp++)
    {
    mu(i,0) = scaleinv * *weightp * mu(i,0);
    }
  muy.mult(X.transposed(),mu);

  beta.assign(standnormal);
  prec_env.solve(muy,betahelp);
  prec_env.solveU(beta,betahelp);

  spline.mult(X,beta);
  likep->add_linearpred(spline);

  acceptance++;

  }


void FULLCOND_kriging2::update_iwls()
  {

  unsigned i,j;

  datamatrix help;
  datamatrix help2;
  double * weightp;

  double logold = - 0.5*Kenv.compute_quadform(betaold,0)/sigma2;
  logold += likep->loglikelihood();

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {

    likep->compute_weight(W,column,true);
    likep->tilde_y(mu,spline,column,true,W);

    help = X;
    help2 = X.transposed();
    weightp = W.getV();
    for(i=0;i<likep->get_nrobs();i++,weightp++)
      {
      for(j=0;j<nrpar;j++)
        {
        help(i,j) *= sqrt(*weightp);
        help2(j,i) *= *weightp;
        }
      }

    help = help.sscp();
    XX_env = envmatdouble(help);

    muy.mult(help2,mu);
    }
  else
    {
    likep->tilde_y(mu,spline,column,true,W);
//    compute_XWtildey(W,1.0);
    help2 = X.transposed();
    help2.multdiagback(W);
    muy.mult(help2,mu);
    }

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

  prec_env.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.solveU(beta,betahelp);

//  add_linearpred_multBS(beta,betaold,true);
  spline.mult(X,betaold);
  likep->substr_linearpred_m(spline,column);
  spline.mult(X,beta);
  likep->add_linearpred_m(spline,column);

  betahelp.minus(beta,betahelp);

  double qold = - 0.5*prec_env.compute_quadform(betahelp,0);

  double lognew = - 0.5*Kenv.compute_quadform(beta,0)/sigma2;
  lognew += likep->loglikelihood();

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    qold += 0.5*prec_env.getLogDet();

    likep->compute_weight(W,column,true);
    likep->tilde_y(mu,spline,column,true,W);

    help = X;
    help2 = X.transposed();
    weightp = W.getV();
    for(i=0;i<likep->get_nrobs();i++,weightp++)
      {
      for(j=0;j<nrpar;j++)
        {
        help(i,j) *= sqrt(*weightp);
        help2(j,i) *= *weightp;
        }
      }

    help = help.sscp();
    XX_env = envmatdouble(help);

    muy.mult(help2,mu);

    prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);
    }
  else
    {
    likep->tilde_y(mu,spline,column,true,W);
//    compute_XWtildey(W,1.0);
    help2 = X.transposed();
    help2.multdiagback(W);
    muy.mult(help2,mu);
    }

  prec_env.solve(muy,betahelp);

  betahelp.minus(betaold,betahelp);
  double qnew = - 0.5*prec_env.compute_quadform(betahelp,0);

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    qnew += 0.5*prec_env.getLogDet();
    }

  double alpha = lognew + qnew - logold - qold;
  double u = log(uniform());

  if(u<=alpha)
    {
    acceptance++;
    betaold.assign(beta);
    }
  else
    {
//    add_linearpred_multBS(betaold,beta,true);
    spline.mult(X,beta);
    likep->substr_linearpred_m(spline,column);
    spline.mult(X,betaold);
    likep->add_linearpred_m(spline,column);
    beta.assign(betaold);
    }

  }


void FULLCOND_kriging2::update_iwlsmode()
  {

  unsigned i,j;

  datamatrix help;
  datamatrix help2;
  double * weightp;

  double logold = likep->loglikelihood(true)
                - 0.5*Kenv.compute_quadform(betaold,0)/sigma2;

//  add_linearpred_multBS(beta_mode,betaold,true);
  spline.mult(X,betaold);
  likep->substr_linearpred_m(spline,column);
  spline.mult(X,beta_mode);
  likep->add_linearpred_m(spline,column);

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    likep->compute_IWLS_weight_tildey(W,mu,column,true);
    mu.plus(spline,mu);

//    compute_XWXenv(W);
    help = X;
    help2 = X.transposed();
    weightp = W.getV();
    for(i=0;i<likep->get_nrobs();i++,weightp++)
      {
      for(j=0;j<nrpar;j++)
        {
        help(i,j) *= sqrt(*weightp);
        help2(j,i) *= *weightp;
        }
      }

    help = help.sscp();
    XX_env = envmatdouble(help);

    muy.mult(help2,mu);
    }
  else
    {
    likep->tilde_y(mu,spline,column,true,W);
//    compute_XWtildey(W,1.0);
    help2 = X.transposed();
    help2.multdiagback(W);
    muy.mult(help2,mu);
    }

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

  prec_env.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.solveU(beta,betahelp);

//  add_linearpred_multBS(beta,beta_mode,true);
  spline.mult(X,beta_mode);
  likep->substr_linearpred_m(spline,column);
  spline.mult(X,beta);
  likep->add_linearpred_m(spline,column);

  beta_mode.assign(betahelp);

  betahelp.minus(beta,beta_mode);
  double qold = - 0.5*prec_env.compute_quadform(betahelp,0);

  betahelp.minus(betaold,beta_mode);
  double lognew = likep->loglikelihood(true)
                - 0.5*Kenv.compute_quadform(beta,0)/sigma2;
  double qnew = - 0.5*prec_env.compute_quadform(betahelp,0);

  double alpha = lognew + qnew - logold - qold;
  double u = log(uniform());

  if(u<=alpha)
    {
    acceptance++;
    betaold.assign(beta);
    }
  else
    {
//    add_linearpred_multBS(betaold,beta,true);
    spline.mult(X,beta);
    likep->substr_linearpred_m(spline,column);
    spline.mult(X,betaold);
    likep->add_linearpred_m(spline,column);

    beta.assign(betaold);
    }

  }


bool FULLCOND_kriging2::posteriormode_converged(const unsigned & itnr)
  {
  return likep->posteriormode_converged_fc(beta,beta_mode,itnr);
  }


bool FULLCOND_kriging2::posteriormode(void)
  {

  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);

  unsigned i,j;

  datamatrix help = X;
  datamatrix help2 = X.transposed();
  double * weightp = likep->get_weightiwls().getV()+column;
  unsigned columns = likep->get_weightiwls().cols();
  for(i=0;i<likep->get_nrobs();i++,weightp+=columns)
    {
    for(j=0;j<nrpar;j++)
      {
      help(i,j) *= sqrt(*weightp);
      help2(j,i) *= *weightp;
      }
    }

  help = help.sscp();
  XX_env = envmatdouble(help);

  prec_env.addto(XX_env,Kenv,1.0,lambda);

  likep->substr_linearpred_m(spline,column);

  lambda_prec = lambda;

  likep->compute_workingresiduals(column);

  muy.mult(help2,likep->get_workingresiduals());

  prec_env.solve(muy,beta);

  spline.mult(X,beta);
  likep->add_linearpred_m(spline,column);

  double * fchelpbetap = fchelp.getbetapointer();
  vector<int>::iterator freqwork = freqoutput.begin();
  int * workindex = index.getV();
  for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
    {
    if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
      {
      *fchelpbetap = spline(*workindex,0);
      fchelpbetap++;
      }
    }

  fchelp.posteriormode();
  return FULLCOND_nonp_basis::posteriormode();

  }


double FULLCOND_kriging2::compute_quadform()
  {
  return Kenv.compute_quadform(beta,0);
  }

void FULLCOND_kriging2::outoptions()
  {
  optionsp->out("OPTIONS FOR KRIGING TERM: " + title + "\n",true);
  optionsp->out("\n");
  optionsp->out("  Correlation function: Matern\n");
  optionsp->out("  Parameter nu: " + ST::doubletostring(nu,2) + "\n");
  optionsp->out("  Parameter rho: " + ST::doubletostring(rho,4) + "\n");
  optionsp->out("\n");
  if(!full)
    {
    optionsp->out("  Number of knots: " + ST::inttostring(nrknots) + "\n");
    optionsp->out("  Corresponds to " + ST::doubletostring(100*(double)nrknots/(double)nrdiffobs,3) + "% of the " + ST::inttostring(nrdiffobs) + " different observation points\n");
    optionsp->out("\n");
    }
  if(nrknots<nrdiffobs && spacefill)
    {
    optionsp->out("  Options for the space-filling algorithm:\n");
    optionsp->out("  Maximum number of iteration steps: " + ST::inttostring(maxsteps) + "\n");
    optionsp->out("  p: " + ST::doubletostring(p,5) + "\n");
    optionsp->out("  q: " + ST::doubletostring(q,5) + "\n");
    optionsp->out("\n");
    }
  }

void FULLCOND_kriging2::outresults()
  {
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

  for(i=0,j=0;i<X.rows();i++,indexit++,freqwork++,k+=*indexit)
    {
    if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
      {

      if(mapexisting)
        {
        regionnameshelp[j]=regionnames[k];
        }
      j++;
      }
    }

  fchelp.outresults();
  FULLCOND::outresults();

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

  ST::string outest=pathcurrent;
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
    optionsp->out("  objectname.drawmap " + ST::inttostring(fcnumber) + "\n");
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
  outres << "pmean   ";
  outres << "pqu"  << l1 << "   ";
  outres << "pqu"  << l2 << "   ";
  outres << "pmed   ";
  outres << "pqu"  << u1 << "   ";
  outres << "pqu"  << u2 << "   ";
  outres << "pcat" << level1 << "   ";
  outres << "pcat" << level2 << "   ";
  outres << endl;

  double * workmean = fchelp.get_betameanp();
  double * workbetaqu_l1_lower_p = fchelp.get_beta_lower1_p();
  double * workbetaqu_l2_lower_p = fchelp.get_beta_lower2_p();
  double * workbetaqu50 = fchelp.get_betaqu50p();
  double * workbetaqu_l1_upper_p = fchelp.get_beta_upper1_p();
  double * workbetaqu_l2_upper_p = fchelp.get_beta_upper2_p();
  vector<double>::iterator workxvalues = xvalues.begin();
  vector<double>::iterator workyvalues = yvalues.begin();

  for(i=0;i<xvalues.size();i++,workmean++,workbetaqu50++,
                     workbetaqu_l1_lower_p++,workbetaqu_l2_lower_p++,
                     workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                     workxvalues++,workyvalues++)
    {
    outres << (i+1) << "   ";
    if(mapexisting)
      outres << regionnameshelp[i] << "   ";
    outres << *workxvalues << "   ";
    outres << *workyvalues << "   ";
    outres << *workmean << "   ";
    outres << *workbetaqu_l1_lower_p << "   ";
    outres << *workbetaqu_l2_lower_p << "   ";
    outres << *workbetaqu50 << "   ";
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

  }

void FULLCOND_kriging2::init_names(const vector<ST::string> & na)
  {
  FULLCOND::init_names(na);
  ST::string underscore = "\\_";
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
    priorassumptions.push_back("p: " + ST::doubletostring(p,5));
    priorassumptions.push_back("q: " + ST::doubletostring(q,5));
    }
  }

ST::string FULLCOND_kriging2::getinfo(void)
  {
  if(mapexisting)
    return mapname;
  else
    return title;
  }

} // end: namespace MCMC

//---------------------------------------------------------------------------
#if !defined(__BUILDING_GNU)
#pragma package(smart_init)
#endif


