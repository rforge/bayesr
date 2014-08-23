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



#include "spline_basis.h"

namespace MCMC
{
                            

void spline_basis::change_K(void)
  {

  unsigned i;

  if(type==RW1)
    {
    if(predictright)
      {
      for(i=0;i<nrparpredictright;i++)
        {
        K.set(nrpar-2-i,nrpar-2-i,1.0);
        K.set(nrpar-1-i,nrpar-1-i-1,0.0);

        Kenv.setDiag(nrpar-2-i,1.0);
        Kenv.set(nrpar-1-i,nrpar-1-i-1,0.0);
        }
      }
    if(predictleft)
      {
      for(i=0;i<nrparpredictleft;i++)
        {
        K.set(1+i,1+i,1.0);
        K.set(i,i+1,0.0);

        Kenv.setDiag(1+i,1.0);
        Kenv.set(i,i+1,0.0);
        }
      K.set(nrparpredictleft+1,nrparpredictleft+1+1,-1.0);

      Kenv.set(nrparpredictleft+1,nrparpredictleft+1+1,-1.0);
      }
    }
  else if(type==RW2)
    {
    if(predictright)
      {
      for(i=0;i<nrparpredictright;i++)
        {
        K.set(nrpar-2-i,nrpar-2-i,1.0);
        K.set(nrpar-1-i,nrpar-1-i-1,0.0);
        K.set(nrpar-1-i,nrpar-1-i-2,0.0);

        Kenv.setDiag(nrpar-2-i,1.0);
        Kenv.set(nrpar-1-i,nrpar-1-i-1,0.0);
        Kenv.set(nrpar-1-i,nrpar-1-i-2,0.0);
        }
      K.set(nrpar-1-nrparpredictright-1,nrpar-1-nrparpredictright-1,5.0);
      K.set(nrpar-1-nrparpredictright,nrpar-1-nrparpredictright-1,-2.0);

      Kenv.setDiag(nrpar-1-nrparpredictright-1,5.0);
      Kenv.set(nrpar-1-nrparpredictright,nrpar-1-nrparpredictright-1,-2.0);
      }
    if(predictleft)
      {
      for(i=0;i<nrparpredictleft;i++)
        {
        K.set(1+i,1+i,1.0);
        K.set(0+i,1+i,0.0);
        K.set(0+i,2+i,0.0);

        Kenv.setDiag(1+i,1.0);
        Kenv.set(0+i,1+i,0.0);
        Kenv.set(0+i,2+i,0.0);
        }
      K.set(nrparpredictleft+1,nrparpredictleft+1,5.0);
      K.set(nrparpredictleft,nrparpredictleft+1,-2.0);

      Kenv.setDiag(nrparpredictleft+1,5.0);
      Kenv.set(nrparpredictleft,nrparpredictleft+1,-2.0);
      }
    }
/*
  ofstream out("c:\\cprog\\K.raw");
  K.print2(out);
  out.close();

  ofstream out2("c:\\cprog\\Kenv.raw");
  Kenv.print2(out2);
  out2.close();
*/
  }


void spline_basis::make_index2(void)
  {
  unsigned i;

  index2.push_back(index(0,0));
  for(i=1;i<likep->get_nrobs();i++)
    index2.push_back(index(i,0)-index(i-1,0));
  }


void spline_basis::init_data_varcoeff(const datamatrix & intvar,
   double add)
  {

  // HIER GEHTS WEITER

  unsigned i;
  datamatrix ivar(intvar.rows(),1);
  for (i=0;i<intvar.rows();i++)
    {
    ivar(i,0) = intvar(i,0)+add;
    }


  make_index(data,ivar);
  make_Bspline(data);
  make_BS(ivar);


/*

  data = datamatrix(intvar.rows(),1);
  data2 = datamatrix(intvar.rows(),1);

  unsigned i;
  double * workdata = data.getV();
  double * workdata2 = data2.getV();
  int * workindex = index.getV();
  for (i=0;i<data.rows();i++,workdata++,workdata2++,workindex++)
    {
    *workdata = intvar(*workindex,0)+add;
    *workdata2 = (*workdata) * (*workdata);
    }

  */
  }


unsigned spline_basis::get_nreffects(effecttype t)
  {
  if (varcoeff)
    {
    if (t==MCMC::fvar_current || t==MCMC::fvar_mean  || t==MCMC::fvar_median)
      return 1;
    else
      return 3;
    }
  else
    return 2;
  }


void spline_basis::get_effectmatrix(datamatrix & e, vector<ST::string> & enames,
                                    unsigned be, unsigned en, effecttype t)
  {

  int * workindex = index.getV();

  double * worksplinehelp;
  double * worke = e.getV()+be;
  int cole = e.cols();
  if (t==MCMC::current || t==MCMC::fvar_current)
    worksplinehelp = splinehelp.getV();
//  else if (t==MCMC::mean || t==MCMC::fvar_mean)
//    workbeta = betamean.getV();
//  else
//    workbeta = betaqu50.getV();

  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();

  int j;
  unsigned i;

  if (varcoeff)
    {

    if (t==MCMC::fvar_current || t==MCMC::fvar_mean  || t==MCMC::fvar_median)
      {

      for (i=0;i<spline.rows();i++,workindex++,worke+=cole)
        {
        e(i,be) = spline(i,0)-intercept;
        }

      }
    else
      {

      /*
      double * workdata=data.getV();

      vector<ST::string>::iterator effit = effectvalues.begin();
      int t;
      enames.push_back("f_"+datanames[0]+"_"+datanames[1]);
      enames.push_back(datanames[0]);
      enames.push_back(datanames[1]);

      for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend,++effit)
        {
        if (*itbeg != -1)
          {
          for(j=*itbeg;j<=*itend;j++,workindex++,workdata++)
            {
            e(*workindex,be) = *workbeta*(*workdata);
            t = (*effit).strtodouble(e(*workindex,be+1));
            e(*workindex,be+2) = *workdata;
            }
          }
        }

      */
      }
    }
  else
    {
    /*
    vector<ST::string>::iterator effit = effectvalues.begin();
    int t;

    enames.push_back("f_"+datanames[0]);
    enames.push_back(datanames[0]);


    for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend,++effit)
      {
      if (*itbeg != -1)
        {
        for (k=(*itbeg);k<=(*itend);k++,workindex++)
          {
          e(*workindex,be) = *workbeta;
          t = (*effit).strtodouble(e(*workindex,be+1));
          }
        }
      }

    */
    }

  }


  // CONSTRUCTOR

spline_basis::spline_basis(MCMCoptions * o, DISTRIBUTION * dp,
                FULLCOND_const * fcc, const fieldtype & ft,
                const ST::string & ti, const unsigned & nrk, const unsigned & degr,
                const MCMC::knotpos & kp, const int & gs, const ST::string & fp,
                const ST::string & pres, const bool & deriv, const double & lk,
                const double & uk, const double & lg, const double & ug,
                const unsigned & c)
  : FULLCOND_nonp_basis(o,dp,ft,ti,fp,pres,c)
  {

  fcconst = fcc;

  lambdaconst = false;
  outbsplines = false;
  lambda_prec = -1.0;

  lowerknot = lk;
  upperknot = uk;

  lowergrid = lg;
  uppergrid = ug;

  fctype = nonparametric;

  predictright = false;
  predictleft = false;
  nrparpredictright = 0;
  nrparpredictleft = 0;

  derivative = deriv;

  increasing = false;
  decreasing = false;

  plotstyle = plotnonp;
  pathresult = pres;

  nrknots = nrk;
  degree = degr;
  knpos = kp;
  gridsize = gs;

  setbeta(nrknots-1+degree,1,0);
//  betaold = datamatrix(nrpar,1,0);
//  betahelp = datamatrix(nrpar,1,0);

  intercept = 0.0;
  spline = datamatrix(likep->get_nrobs(),1,0);

  pseudocontourprob = false;

  XX_env = envmatdouble();
  prec_env = envmatdouble();

  }

  // CONSTRUCTOR für REML

spline_basis::spline_basis(MCMCoptions * o,
                      const datamatrix & d, const unsigned & nrk, const unsigned & degr,
                      const knotpos & kp, const fieldtype & ft, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres, const double & l,
                      const double & sl, const bool & catsp, const double & lg,
                      const double & ug, const double & lk, const double & uk,
                      const int & gs, const double & rv)
  : FULLCOND_nonp_basis(o,ti)
  {
  catspecific = catsp;

  pseudocontourprob = false;

  gridsize = gs;

  lowergrid = lg;
  uppergrid = ug;
  lowerknot = lk;
  upperknot = uk;

  reference = rv;
  if(reference == -9999)
    refcheck=false;
  else
    refcheck = true;

  fctype = nonparametric;

  lambdaconst = false;
  outbsplines = false;
  lambda_prec = -1.0;

  derivative = false;

  increasing = false;
  decreasing = false;

  plotstyle = plotnonp;

  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;

  nrknots = nrk;
  degree = degr;
  knpos = kp;

  varcoeff = false;

  transformnonlinear = false;
  transformed =  false;

  type = ft;
  nrknots = nrk;
  degree = degr;
  knpos = kp;

  nrpar = nrknots-1+degree;

//  samplepath = pres;
  samplepath=fp;

  likep = NULL;

  spline = datamatrix(d.rows(),1,0);

  lambda = l;
  startlambda = sl;

  make_index(d);
  make_Bspline(d);

  if(type == RW1)
    dimX = 0;
  else
    dimX = 1;

  dimZ = nrpar-dimX-1;

  if(gridsize>0)
    {
    X_grid = datamatrix(gridsize,dimX,1.0);
    Z_grid = datamatrix(gridsize,dimZ,0.0);
    }

  index2.push_back(index(0,0));
  for(unsigned i=1;i<d.rows();i++)
    index2.push_back(index(i,0)-index(i-1,0));

  unsigned j;
  vector<int>::iterator freqwork = freqoutput.begin();
  int * workindex = index.getV();

  xvalues = datamatrix(nrdiffobs,1,0);
  for(j=0;j<d.rows();j++,freqwork++,workindex++)
    if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
      xvalues(*freqwork,0) = d(*workindex,0);

  if(gridsize>0)
    {
    if(lowergrid<uppergrid)
      {
      }
    else
      {
      lowergrid = d.min(0);
      uppergrid = d.max(0);
      }
    xvaluesgrid = datamatrix(gridsize,1,lowergrid);
    double step = (uppergrid-lowergrid)/(gridsize-1);
    for(j=1; j<gridsize; j++)
      {
      xvaluesgrid(j,0) = xvaluesgrid(j-1,0) + step;
      }
    make_DG_REML();
    }
  }

  // CONSTRUCTOR für REML VCM

spline_basis::spline_basis(MCMCoptions * o,const datamatrix & d1,
                      const datamatrix & d2, const unsigned & nrk, const unsigned & degr,
                      const knotpos & kp, const fieldtype & ft, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres, const double & l,
                      const double & sl, const bool & catsp, const bool & ctr, const double & rv)
  : FULLCOND_nonp_basis(o,ti)
  {
  catspecific = catsp;
  centervcm = ctr;

  lowerknot=0.0;
  upperknot=0.0;

  pseudocontourprob = false;

//------------------------------------------------------------------------------

  reference = rv;
  if(reference == -9999)
    refcheck=false;
  else
    refcheck = true;

  fctype = nonparametric;

  data_forfixed = d2;

  lambdaconst = false;
  outbsplines = false;
  lambda_prec = -1.0;

  derivative = false;

  increasing = false;
  decreasing = false;

  plotstyle = plotnonp;

  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;

  nrknots = nrk;
  degree = degr;
  knpos = kp;

  gridsize = -1;

  varcoeff = true;

  transformnonlinear = false;
  transformed =  false;

  type = ft;
  nrknots = nrk;
  degree = degr;
  knpos = kp;

  nrpar = nrknots-1+degree;

//------------------------------------------------------------------------------

//  samplepath = pres;
  samplepath=fp;

  likep = NULL;

  if(type == RW1)
    dimX = 1;
  else
    dimX = 2;

  dimZ = nrpar-dimX;

  if(centervcm)
    {
    dimX = dimX-1;
    }

  X_VCM = datamatrix(d1.rows(),dimX,1.0);
  Z_VCM = datamatrix(d1.rows(),dimZ,0.0);

//------------------------------------------------------------------------------

  spline = datamatrix(d1.rows(),1,0);

  lambda = l;
  startlambda = sl;

  make_index(d1,d2);
  make_Bspline(d1);
  make_BS(d2);

  index2.push_back(index(0,0));
  for(unsigned i=1;i<d1.rows();i++)
    index2.push_back(index(i,0)-index(i-1,0));

  unsigned j;
  vector<int>::iterator freqwork = freqoutput.begin();
  int * workindex = index.getV();

  xvalues = datamatrix(nrdiffobs,1,0);
  for(j=0;j<d1.rows();j++,freqwork++,workindex++)
    if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
      xvalues(*freqwork,0) = d1(*workindex,0);
  }

  // COPY CONSTRUCTOR

spline_basis::spline_basis(const spline_basis & sp)
  : FULLCOND_nonp_basis(FULLCOND_nonp_basis(sp))
  {

  fcconst = sp.fcconst;

  pseudocontourprob = sp.pseudocontourprob;
  approx = sp.approx;
  lengthstart = sp.lengthstart;

  lambdaconst = sp.lambdaconst;
  outbsplines = sp.outbsplines;
  lambda_prec = sp.lambda_prec;

  predictright = sp.predictright;
  nrparpredictright = sp.nrparpredictright;
  predictleft = sp.predictleft;
  nrparpredictleft = sp.nrparpredictleft;

  derivative = sp.derivative;
  index2 = sp.index2;

  increasing = sp.increasing;
  decreasing = sp.decreasing;

  interactvar = sp.interactvar;

  W = sp.W;
  betaold = sp.betaold;
  betaprop = sp.betaprop;

  mu = sp.mu;
  muy = sp.muy;
  betahelp = sp.betahelp;
  standnormal = sp.standnormal;
  XX = sp.XX;
  prec = sp.prec;
  prec_env = sp.prec_env;
  XX_env = sp.XX_env;

  fchelp = sp.fchelp;
  fcderivative = sp.fcderivative;

  Bderivative = sp.Bderivative;
  splinederivative = sp.splinederivative;

  Bcolmean = sp.Bcolmean;

  nrknots = sp.nrknots;
  degree = sp.degree;
  nrdiffobs = sp.nrdiffobs;

  gridsize = sp.gridsize;
  intercept = sp.intercept;
  knpos = sp.knpos;

  freq = sp.freq;
  freqoutput = sp.freqoutput;
  firstnonzero = sp.firstnonzero;
  lastnonzero = sp.lastnonzero;
  knot = sp.knot;

  xvalues = sp.xvalues;
  spline = sp.spline;
  splinehelp = sp.splinehelp;
  betaweight = sp.betaweight;

  B = sp.B;
  BS = sp.BS;
  G = sp.G;

  X_VCM = sp.X_VCM;
  Z_VCM = sp.Z_VCM;

  begcol = sp.begcol;
  DG = sp.DG;
  DGfirst = sp.DGfirst;

  lowergrid = sp.lowergrid;
  uppergrid = sp.uppergrid;
  lowerknot = sp.lowerknot;
  upperknot = sp.upperknot;

  X_grid = sp.X_grid;
  Z_grid = sp.Z_grid;

  xvaluesgrid = sp.xvaluesgrid;

  reference = sp.reference;
  refcheck = sp.refcheck;
  X_ref = sp.X_ref;
  Z_ref = sp.Z_ref;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const spline_basis & spline_basis::operator=(const spline_basis & sp)
  {
  if(&sp==this)
    return *this;
  FULLCOND_nonp_basis::operator=(FULLCOND_nonp_basis(sp));

  fcconst = sp.fcconst;

  pseudocontourprob = sp.pseudocontourprob;
  approx = sp.approx;
  lengthstart = sp.lengthstart;

  lambdaconst = sp.lambdaconst;
  outbsplines = sp.outbsplines;
  lambda_prec = sp.lambda_prec;

  predictright = sp.predictright;
  nrparpredictright = sp.nrparpredictright;
  predictleft = sp.predictleft;
  nrparpredictleft = sp.nrparpredictleft;

  derivative = sp.derivative;
  index2 = sp.index2;

  interactvar = sp.interactvar;
  W = sp.W;
  betaold = sp.betaold;
  betaprop = sp.betaprop;

  mu = sp.mu;
  muy = sp.muy;
  betahelp = sp.betahelp;
  standnormal = sp.standnormal;
  XX = sp.XX;
  prec = sp.prec;
  prec_env = sp.prec_env;
  XX_env = sp.XX_env;

  fchelp = sp.fchelp;
  fcderivative = sp.fcderivative;

  Bderivative = sp.Bderivative;
  splinederivative = sp.splinederivative;

  Bcolmean = sp.Bcolmean;

  nrknots = sp.nrknots;
  degree = sp.degree;
  nrdiffobs = sp.nrdiffobs;

  gridsize = sp.gridsize;
  intercept = sp.intercept;
  knpos = sp.knpos;

  freq = sp.freq;
  freqoutput = sp.freqoutput;
  firstnonzero = sp.firstnonzero;
  lastnonzero = sp.lastnonzero;
  knot = sp.knot;

  xvalues = sp.xvalues;
  spline = sp.spline;
  splinehelp = sp.splinehelp;
  betaweight = sp.betaweight;

  B = sp.B;
  BS = sp.BS;
  G = sp.G;

  X_VCM = sp.X_VCM;
  Z_VCM = sp.Z_VCM;

  begcol = sp.begcol;
  DG = sp.DG;
  DGfirst = sp.DGfirst;

  lowergrid = sp.lowergrid;
  uppergrid = sp.uppergrid;
  lowerknot = sp.lowerknot;
  upperknot = sp.upperknot;

  X_grid = sp.X_grid;
  Z_grid = sp.Z_grid;

  xvaluesgrid = sp.xvaluesgrid;

  reference = sp.reference;
  refcheck = sp.refcheck;
  X_ref = sp.X_ref;
  Z_ref = sp.Z_ref;


  return *this;
  }


void spline_basis::make_index(const datamatrix & moddata)
  {
// index berechnen
  index = statmatrix<int>(moddata.rows(),1);
  index.indexinit();
  moddata.indexsort(index,0,moddata.rows()-1,0,0);
// freq berechnen
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
// freqoutput, nrdiffobs
  freqoutput = freq;
  nrdiffobs = i+1;

  }

void spline_basis::make_index(const datamatrix & em,const datamatrix & ia)
  {
// index berechnen
  index = statmatrix<int>(em.rows(),1);
  index.indexinit();
  em.indexsort(index,0,em.rows()-1,0,0);
// freq berechnen
  unsigned i,j,k,beg,end;

  int *workindex = index.getV();
  freq.reserve(em.rows());

  workindex++;
  freq.push_back(0);
  i = 0;
  j = 1;
  while(j<em.rows())
    {
    while(j<em.rows() && em(*workindex,0)!=em(*(workindex-1),0))
      {
      i++;
      freq.push_back(i);
      j++;
      workindex++;
      }
    beg = j-1;
    while(j<em.rows() && em(*workindex,0)==em(*(workindex-1),0))
      {
      workindex++;
      j++;
      }
    end = j-1;
    if(end!=beg)
      {
      ia.indexsort(index,beg,end,0,0);
      for(k=beg+1;k<=end;k++)
        {
        if(ia(index(k,0),0) != ia(index(k-1,0),0))
          {
          i++;
          }
        freq.push_back(i);
        }
      }
    }

// initialisiere 'freqoutput'
  workindex = index.getV();
  freqoutput.reserve(em.rows());
  workindex++;
  freqoutput.push_back(0);
  i = 0;
  for(j=1;j<em.rows();j++,workindex++)
    {
    if ( em(*workindex,0) != em(*(workindex-1),0))
      i++;
    freqoutput.push_back(i);
    }
// nrdiffobs
  nrdiffobs = i+1;

  }

/*
void spline_basis::make_Bspline(const datamatrix & md)
  {

  unsigned i,j,l,posmd;

  double min = md(index(0,0),0);
  double max = md(index(md.rows()-1,0),0);
  double dist = max - min;
  SparseMatrix M(*(freq.end()-1)+1,nrpar,degree+1);
  vector<int>::iterator freqwork;

  min -= 0.01*dist;
  max += 0.01*dist;

  if(knpos == equidistant)
    {

    dist = (max - min)/(nrknots-1);
    knot.push_back(min - degree*dist);
    for(i=1;i<nrknots+2*degree;i++)
      knot.push_back(knot[i-1] + dist);

    for(j=0;j<nrpar;j++)
      {
      posmd = 0;
      freqwork = freq.begin();
      while(posmd<md.rows() && md(index(posmd,0),0)<knot[j])
        {
        freqwork++;
        posmd++;
        }
      firstnonzero.push_back(posmd);
      while(posmd<md.rows() && knot[j]<=md(index(posmd,0),0) && md(index(posmd,0),0)<knot[j+1])
        {
        if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
          M.put(*freqwork,j,(md(index(posmd,0),0)-knot[j])/dist);
        freqwork++;
        posmd++;
        }
      while(posmd<md.rows() && knot[j+1]<=md(index(posmd,0),0) && md(index(posmd,0),0)<knot[j+2])
        {
        if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
          M.put(*freqwork,j,(knot[j+2]-md(index(posmd,0),0))/dist);
        freqwork++;
        posmd++;
        }
      lastnonzero.push_back(posmd-1);
      } // end: for(j=0;j<nrpar;j++)

    } // end: if(knpos == equidistant)

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
    for(j=0;j<nrpar;j++)
      {
      posmd = 0;
      freqwork = freq.begin();
      while(posmd<md.rows() && md(index(posmd,0),0)<knot[j])
        {
        freqwork++;
        posmd++;
        }
      firstnonzero.push_back(posmd);
      while(posmd<md.rows() && knot[j]<=md(index(posmd,0),0) && md(index(posmd,0),0)<knot[j+1])
        {
        if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
          M.put(*freqwork,j,(md(index(posmd,0),0)-knot[j])/(knot[j+1]-knot[j]));
        freqwork++;
        posmd++;
        }
      while(posmd<md.rows() && knot[j+1]<=md(index(posmd,0),0) && md(index(posmd,0),0)<knot[j+2])
        {
        if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
          M.put(*freqwork,j,(knot[j+2]-md(index(posmd,0),0))/(knot[j+2]-knot[j+1]));
        freqwork++;
        posmd++;
        }
      lastnonzero.push_back(posmd-1);
      } // end: for(j=0;j<nrpar;j++)
    } // end: else if(knpos == quantiles)

  if(degree > 1)
    {
    datamatrix help2(*(freq.end()-1)+1,nrpar,0);
    SparseMatrix help(*(freq.end()-1)+1,nrpar,degree+1);
    for(l=2;l<degree+1;l++)
      {
      lastnonzero.pop_front();
      lastnonzero.push_back(*(lastnonzero.end()-1));
      for(j=0;j<nrpar;j++)
        {
        posmd = 0;
        freqwork = freq.begin();
        for(posmd=0;posmd<md.rows();posmd++,freqwork++)
          {
          if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
            help.put(*freqwork,j,(md(index(posmd,0),0)-knot[j])*M(*freqwork,j)/(knot[j+l]-knot[j])
                   +(knot[j+l+1]-md(index(posmd,0),0))*M(*freqwork,j+1)/(knot[j+l+1]-knot[j+1]));
          }
        }
        M = help;
        help = SparseMatrix(help2);
      } // end: for(l=2;l<degree+1;l++)
    } // end: if(degree > 1)

// compute Bcolmean
  Bcolmean = datamatrix(nrpar,1,0);
  Bcolmean.mult(datamatrix(M).transposed(),datamatrix(*(freq.end()-1)+1,1,1/double(*(freq.end()-1)+1)));
//  Bcolmean = datamatrix(1,*(freq.end()-1)+1,1/double(*(freq.end()-1)+1)) * datamatrix(M);

  double *work;
  if(varcoeff)
    {
    B = datamatrix(*(freq.end()-1)+1,degree+1,0);
    work = B.getV();
    }
  else
    {
    BS = datamatrix(*(freq.end()-1)+1,degree+1,0);
    work = BS.getV();
    }
  posmd = 0;
  unsigned k = 0;
  freqwork = freq.begin();
  while(k<nrpar)
    {
    while(posmd<lastnonzero[k]+1)
      {
      if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
        {
        for(j=0;j<degree+1;j++,work++)
          *work = M(*freqwork,j+k);
        begcol.push_back(k);
        }
      freqwork++;
      posmd++;
      }
    k++;
    }

  }
*/

void spline_basis::make_Bspline(const datamatrix & md, const bool & minnull)
  {

  unsigned i=0,j=0,k=0;
  double value=0.0;
  double * work = NULL;

  vector<int>::iterator freqwork;
  datamatrix help;

  unsigned mindomain = 0;
  int maxdomain = md.rows()-1;
// Prädiktion ja/nein
  if(likep != NULL)
    {
    while(likep->get_weight(index(mindomain,0),0)==0 && mindomain<likep->get_nrobs())
      {
      predictleft = true;
      mindomain++;
      }

    while(likep->get_weight(index(maxdomain,0),0)==0 && maxdomain>=0)
      {
      predictright = true;
      maxdomain--;
      }
    }
// berechne x_min, x_max
  double min;
  if(predictleft)
    min = md(index(mindomain,0),0);
  else
    min = md(index(0,0),0);
  double max;
  if(predictright)
    max = md(index(maxdomain,0),0);
  else
    max = md(index(md.rows()-1,0),0);
  double dist = max-min;

  min -= 0.01*dist;
  max += 0.01*dist;

  if(minnull)
    min = 0.0;

// SIMEX

  if(lowerknot<upperknot)
    {
    min=lowerknot;
    max=upperknot;
    }

// Knoten berechnen
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
/*
    knot.push_back(min);
    for(i=1;i<nrknots-1;i++)
      knot.push_back(md.quantile((i*100)/double(nrknots-1),0));
    knot.push_back(max);
*/

    datamatrix mdhelp = md;//.getRowBlock(0,500);

    knot.push_back(min);
    for(i=1;i<nrknots-1;i++)
      {
      double q1 = mdhelp.quantile((i*100)/double(nrknots-1),0);
      double q2 = mdhelp.quantile(((i-1)*100)/double(nrknots-1),0);
      if(q1 != q2)
        knot.push_back(q1);
      }
    knot.push_back(max);

    nrknots = knot.size();

    distfirst = knot[1] - knot[0];
    distlast = knot[nrknots-1] - knot[nrknots-2];

    for(i=0;i<degree;i++)
      {
      knot.push_front(min - (i+1)*distfirst);
      knot.push_back(max + (i+1)*distlast);
      }
    }
// Knoten, falls predictright==true || predictleft==true
  if(predictright)
    {
    while(md(index(md.rows()-1,0),0) > knot[knot.size()-1-degree])
      {
      if(knpos==equidistant)
        knot.push_back(knot[knot.size()-1] + dist);
      else if(knpos==quantiles)
        knot.push_back(knot[knot.size()-1] + (knot[nrknots-1] - knot[nrknots-2]));
      nrknots++;
      nrparpredictright++;
      }
    }
  if(predictleft)
    {
    while(md(index(0,0),0) < knot[degree])
      {
      if(knpos==equidistant)
        knot.push_front(knot[0] - dist);
      else if(knpos==quantiles)
        knot.push_front(knot[0] - (knot[1] - knot[0]));
      nrknots++;
      nrparpredictleft++;
      }
    }

  if(knpos==quantiles)
    {
    if(nrpar > nrknots-1+degree)
      {
      optionsp->outerror("\n");
      optionsp->outerror("WARNING: Reducing the number of basis functions for term " + title + "\n");
      optionsp->outerror("         due to equal quantiles for the knot positions.\n");
      }
    }

  setbeta(nrknots-1+degree,1,0);
// Designmatrix BS bzw. B (bei VCM), lastnonzero, firstnonzero und Bcolmean berechnen
  help = datamatrix(nrpar,1,0.0);
  Bcolmean = datamatrix(nrpar,1,0.0);

  for(i=0;i<nrpar;i++)
    {
    lastnonzero.push_back(-1);
    firstnonzero.push_back(0);
    }

  if(varcoeff)
    {
    B = datamatrix(*(freq.end()-1)+1,degree+1,0.0);
    work = B.getV();
    }
  else
    {
    BS = datamatrix(nrdiffobs,degree+1,0.0);
    work = BS.getV();
    }

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


void spline_basis::make_BS(const datamatrix & ia)
  {

  unsigned i,j;
  vector<int>::iterator freqwork = freq.begin();

  BS = B;
  for(i=0;i<ia.rows();i++,freqwork++)
    {
    if(freqwork == freq.begin() || *freqwork != *(freqwork-1))
      {
      for(j=0;j<degree+1;j++)
        {
        BS(*freqwork,j) *= ia(index(i,0),0);
        }
      }
    }

  }


void spline_basis::make_DG(void)
  {
// DG und DGfirst berechnen
  int i;
  unsigned j,k;
  datamatrix betahelp(nrpar,1,0);
  DG = datamatrix(gridsize,degree+1,0);
  DGfirst = vector<int>(gridsize);

  for(i=0;i<gridsize;i++)
    {
    betahelp.assign( bspline(xvalues(i,0)) );
    j=degree+1;
    while(knot[j] <= xvalues(i,0) && j<nrknots+degree)
      j++;
    for(k=0;k<degree+1;k++)
      DG(i,k) = betahelp(k+j-(degree+1),0);
    DGfirst[i] = j-(degree+1);
    }
  }

void spline_basis::make_DG_REML(void)
  {
// DG und DGfirst berechnen
  int i;
  unsigned j,k;
  datamatrix betahelp(nrpar,1,0);
  DG = datamatrix(gridsize,degree+1,0);
  DGfirst = vector<int>(gridsize);

  for(i=0;i<gridsize;i++)
    {
    betahelp.assign( bspline(xvaluesgrid(i,0)) );
    j=degree+1;
    while(knot[j] <= xvaluesgrid(i,0) && j<nrknots+degree)
      j++;
    for(k=0;k<degree+1;k++)
      DG(i,k) = betahelp(k+j-(degree+1),0);
    DGfirst[i] = j-(degree+1);
    }
  }

void spline_basis::init_fchelp(const datamatrix & d)
  {

  int i;
  unsigned j;
// fchelp und xvalues initialisieren
  ST::string path = samplepath.substr(0,samplepath.length()-4)+"_fchelp.raw";
  vector<int>::iterator freqwork = freqoutput.begin();
  int * workindex = index.getV();
  if(gridsize < 0)
    {
    xvalues = datamatrix(nrdiffobs,1,0);
    for(j=0;j<d.rows();j++,freqwork++,workindex++)
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        xvalues(*freqwork,0) = d(*workindex,0);
    fchelp = FULLCOND(optionsp,xvalues,title+"fchelp",nrdiffobs,1,path);
    splinehelp = datamatrix(d.rows(),1,0);
    }
  else
    {
    double xmin, xmax;
    if(lowergrid<uppergrid)
      {
      xmin = lowergrid;
      xmax = uppergrid;
      }
    else
      {
      xmin = d.min(0);
      xmax = d.max(0);
      }

    xvalues = datamatrix(gridsize,1);
    for(i=0;i<gridsize;i++)
      xvalues(i,0) = xmin + i*(xmax-xmin)/double(xvalues.rows()-1);
    fchelp = FULLCOND(optionsp,xvalues,title+"fchelp",gridsize,1,path);
    splinehelp = datamatrix(gridsize,1,0);
    make_DG();
    }
  fchelp.setflags(MCMC::norelchange | MCMC::nooutput);
  fchelp.set_transform(transform);
// fcderivative initialisieren
  if(derivative)
    {
    ST::string pnt = path.substr(0,path.length()-11)+"_fcderivate.raw";
    if(gridsize < 0)
      {
      fcderivative = FULLCOND(optionsp,datamatrix(1,1,0),title+"fcderivative",nrdiffobs,1,pnt);
      splinederivative = datamatrix(nrdiffobs,1,0);
      }
    else
      {
      fcderivative = FULLCOND(optionsp,datamatrix(1,1,0),title+"fcderivative",gridsize,1,pnt);
      splinederivative = datamatrix(gridsize,1,0);
      }
    Bderivative = bsplinemat(derivative,xvalues,nrknots,degree,knpos);
    fcderivative.setflags(MCMC::norelchange | MCMC::nooutput);
    fcderivative.set_transform(transform);
    }


  }


double spline_basis::bspline_rek(unsigned l, unsigned nu, const datamatrix & X)
  {
// Rekursionsformel aus Hämmerlin/Hoffmann
  if(l==0)
    {
    if(knot[nu] <= X(0,0) && X(0,0) < knot[nu+1])
      return 1.0;
    else
      return 0.0;
    }
  else
    {
    return (X(0,0)-knot[nu])*bspline_rek(l-1,nu,X)/(knot[nu+l]-knot[nu])
      +(knot[nu+l+1]-X(0,0))*bspline_rek(l-1,nu+1,X)/(knot[nu+l+1]-knot[nu+1]);
    }
  }


void spline_basis::compute_betaweight(void)
  {

  unsigned i;

  if(knpos == equidistant && degree==1)         // Integrieren (betaweight*beta = \int_a^b Xbeta dx = intercept)
    {
    betaweight = datamatrix(nrpar,1,1);

    betaweight(0,0) = 0.5;
    betaweight(nrpar-1,0) = 0.5;

    for(i=0;i<betaweight.rows();i++)
      betaweight(i,0) /= (nrknots-1);
    }
  else if(knpos == equidistant && degree==2)    // Integrieren
    {
    betaweight = datamatrix(nrpar,1,1);

    betaweight(0,0) = 1/6.0;
    betaweight(nrpar-1,0) = 1/6.0;
    betaweight(1,0) = 5/6.0;
    betaweight(nrpar-2,0) = 5/6.0;

    for(i=0;i<betaweight.rows();i++)
      betaweight(i,0) /= (nrknots-1);
    }
  else if(knpos == equidistant && degree==3)    // Integrieren
    {
    betaweight = datamatrix(nrpar,1,1);

    betaweight(0,0) = 1/24.0;
    betaweight(nrpar-1,0) = 1/24.0;
    betaweight(1,0) = 12/24.0;
    betaweight(nrpar-2,0) = 12/24.0;
    betaweight(2,0) = 23/24.0;
    betaweight(nrpar-3,0) = 23/24.0;

    for(i=0;i<betaweight.rows();i++)
      betaweight(i,0) /= (nrknots-1);
    }
  else
    {
    betaweight = datamatrix(nrpar,1,1.0/double(nrpar));  // intercept = 1/nrpar * \sum_i beta_i
    }

  }


void spline_basis::compute_intercept(const datamatrix & beta)
  {

  double *workbeta = beta.getV();
  double *workbetaweight = betaweight.getV();
  unsigned i;

  intercept = 0.0;
  for(i=0;i<nrpar;i++,workbetaweight++,workbeta++)
    intercept += *workbetaweight * *workbeta;
  }


void spline_basis::compute_intercept(void)
  {

  double *workbeta = beta.getV();
  double *workbetaweight = betaweight.getV();
  unsigned i;

  intercept = 0.0;
  for(i=0;i<nrpar;i++,workbetaweight++,workbeta++)
    intercept += *workbetaweight * *workbeta;

  }


void spline_basis::subtr_spline(void)
  {

  unsigned i;
  double *workspline = spline.getV();
/*
  datamatrix *workl = &(likep->get_linearpred(true));
  for(i=0;i<freq.size();i++,workspline++)
    (*workl)(i,column) -= (*workspline - intercept);
*/
  unsigned col = (likep->get_linearpred(true)).cols();
  double *lp = (likep->get_linearpred(true)).getV() + column;
  if ((varcoeff==true) && (center==true))
    {
    double * workint = interactvar.getV();
    for(i=0;i<likep->get_nrobs();i++,workspline++,lp+=col,workint++)
      *lp -= (*workspline - intercept*(*workint));

    }
  else
    {
    for(i=0;i<likep->get_nrobs();i++,workspline++,lp+=col)
      *lp -= (*workspline - intercept);
    }

  }


void spline_basis::add_linearpred_multBS(const bool & current)
  {

  double *workBS;
  double *workbeta;
  double *lp;

  unsigned j,k;
  unsigned col = degree+1;
  int lpcols = likep->get_linearpred(current).cols();
  int i,stop;

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator workindex2 = index2.begin();
// spline = 0 setzen
  double * workspline = spline.getV();
  for(j=0;j<spline.rows();j++,workspline++)
    *workspline = 0.0;

  lp = likep->get_linearpred(current).getV() + column + *workindex2*lpcols;
  workspline = spline.getV() + *workindex2;
// Spaltenweise multiplizieren
  i = 0;
  k = 0;
  workBS = BS.getV();
  while (k<nrpar)
    {
    stop = lastnonzero[k];
//    while (i<lastnonzero[k]+1)
    while (i <= stop)
      {
      workbeta = beta.getV() + k;
      for(j=0;j<col;j++,workBS++,workbeta++)
        {
        *lp += *workBS * *workbeta;
        *workspline += *workBS * *workbeta;
        }
      if((freqwork+1)!=freq.end() && *freqwork==*(freqwork+1))
        {
        workBS -= col;
        workbeta -= col;
        }
      i++;
      freqwork++;
      workindex2++;
      workspline += *workindex2;
      lp += *workindex2*lpcols;
      }
    k++;
    }

  }


void spline_basis::add_linearpred_multBS(const datamatrix & beta,const bool & current)
  {

  double *workBS;
  double *workbeta;
  double *lp;

  unsigned j,k;
  unsigned col = degree+1;
  int lpcols = likep->get_linearpred(current).cols();
  int i,stop;
// spline = 0 setzen
  double * workspline = spline.getV();
  for(j=0;j<spline.rows();j++,workspline++)
    *workspline = 0.0;

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator workindex2 = index2.begin();

  lp = likep->get_linearpred(current).getV() + column + *workindex2*lpcols;
  workspline = spline.getV() + *workindex2;
// Spaltenweise multiplizieren
  i = 0;
  k = 0;
  workBS = BS.getV();
  while (k<nrpar)
    {
    stop = lastnonzero[k];
//    while (i<lastnonzero[k]+1)
    while (i <= stop)
      {
      workbeta = beta.getV() + k;
      for(j=0;j<col;j++,workBS++,workbeta++)
        {
        *lp += *workBS * *workbeta;
        *workspline += *workBS * *workbeta;
        }
      if((freqwork+1)!=freq.end() && *freqwork==*(freqwork+1))
        {
        workBS -= col;
        workbeta -= col;
        }
      i++;
      freqwork++;
      workindex2++;
      workspline += *workindex2;
      lp += *workindex2*lpcols;
      }
    k++;
    }

  }


void spline_basis::add_linearpred_multBS(const datamatrix & beta1,const datamatrix & beta2,
                                          const bool & current)
  {

  double *workBS;
  double *workbeta1;
  double *workbeta2;
  double *lp;

  unsigned j,k;
  unsigned col = degree+1;
  int lpcols = likep->get_linearpred(current).cols();
  int i,stop;
// spline = 0 setzen
  double * workspline = spline.getV();
  for(j=0;j<spline.rows();j++,workspline++)
    *workspline = 0.0;

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator workindex2 = index2.begin();

  lp = likep->get_linearpred(current).getV() + column + *workindex2*lpcols;
  workspline = spline.getV() + *workindex2;
// Spaltenweise multiplizieren
  i = 0;
  k = 0;
  workBS = BS.getV();
  while (k<nrpar)
    {
    stop = lastnonzero[k];
//    while (i<lastnonzero[k]+1)
    while (i <= stop)
      {
      workbeta1 = beta1.getV() + k;
      workbeta2 = beta2.getV() + k;
      for(j=0;j<col;j++,workBS++,workbeta1++,workbeta2++)
        {
        *lp += *workBS * (*workbeta1 - *workbeta2);
        *workspline += *workBS * *workbeta1;
        }
      if((freqwork+1)!=freq.end() && *freqwork==*(freqwork+1))
        {
        workBS -= col;
        workbeta1 -= col;
        workbeta2 -= col;
        }
      i++;
      freqwork++;
      workindex2++;
      workspline += *workindex2;
      lp += *workindex2*lpcols;
      }
    k++;
    }

  }


void spline_basis::add_linearpred_multBS_Block(const unsigned a,const unsigned e,const datamatrix & b)
  {

  assert(e < nrpar);

  datamatrix *workl = &(likep->get_linearpred(false));
  double *workBS;
  double *workbeta;
  double *workfcrand;
  int *workindex;

  unsigned l;
  unsigned col = degree+1;
  unsigned k = a;
  int first,last,ende;

  vector<int>::iterator freqwork = freq.begin();
// zeilenweise multiplizieren !
  while(k<e+1 && k<firstnonzero.size())
    {
    first = firstnonzero[k];
    ende = lastnonzero[k]+1;
    workindex = index.getV() + first;
    freqwork = freq.begin();
    if(k<col)
      {
      last = lastnonzero[0]+1;
      workBS = BS.getV() + k;
      }
    else
      {
      freqwork += first;
      last = lastnonzero[k-col+1]+1;
      workBS = BS.getV() + col*(*freqwork) + col-1;
      }
    l=0;
    while(l<col && first<ende)
      {
      workbeta = beta.getV() + a;
      workfcrand = b.getV();
      while(first<last && first<ende && first<(lastnonzero[e]+1))
        {
        (*workl)(*workindex,column) += *(workBS - l) * (*(workfcrand + k-a) - *(workbeta + k-a));
        if((freqwork+1)!=freq.end() && *freqwork!=*(freqwork+1))
          workBS += col;
        first++;
        workindex++;
        freqwork++;
        }
      l++;
      if(k<col)
        last = lastnonzero[l]+1;
      else
        last = lastnonzero[k-col+1+l]+1;
      }
    k++;
    }

  }

/*
void spline_basis::add_linearpred_multBS_Block2(const unsigned a,const unsigned e,const datamatrix & b)
  {

  unsigned i,j;
  double lambda_i;
  double lambda_iprop;

  datamatrix *workl = &(likep->get_linearpred(false));

  unsigned col = degree + 1;
  unsigned beg = firstnonzero[a];
  unsigned end = lastnonzero[e]+1;

  vector<int>::iterator freqwork = freq.begin()+beg;

  betaprop = beta;
  for(i=0;i<b.rows();i++)
    betaprop(a+i,0) = b(i,0);

  for(i=beg;i<end;i++,++freqwork)
    {
    lambda_i = 0.0;
    lambda_iprop = 0.0;
    for(j=0;j<col;j++)
      {
      lambda_i += BS(*freqwork,j) * beta(begcol[*freqwork]+j,0);
      lambda_iprop += BS(*freqwork,j) * betaprop(begcol[*freqwork]+j,0);
      }
    (*workl)(index(i,0),column) += log(lambda_iprop/lambda_i);
    }

  }
*/

datamatrix spline_basis::bspline(const double & x)
  {
// nach Hämmerlin/Hoffmann
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


datamatrix spline_basis::bspline(const double & x, const unsigned & d)
  {
// nach Hämmerlin/Hoffmann
  datamatrix b(nrpar+degree-d,1,0.0);
  datamatrix help(nrpar+d,1,0.0);
  unsigned j;
  double * bwork;
  double * helpwork;

// Grad 0

  for(j=0;j<nrpar;j++)
    if( knot[j]<=x && x<knot[j+1])
      b(j,0) = 1.0;

  for(unsigned l=1;l<=d;l++)
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


double spline_basis::deriv_f(const double & x)
  {
  double res = 0.0;
  double h = knot[degree+1]-knot[degree];
  datamatrix X = bspline(x,degree-1);

  for(unsigned i=0;i<nrpar;i++)
    res += beta(i,0)*(X(i,0)-X(i+1,0));

  return res/h;
  }


double spline_basis::sample_monotonic(const unsigned i, const double m, const double s)
  {
  bool standard = true;
//  bool standard = false;

  double sample;

  if(standard)
    {
    if(i==0)
      {
      if(increasing)
        sample = trunc_normal2(-20,beta(1,0),m,s);
      else
        sample = trunc_normal2(beta(1,0),20,m,s);
      }
    else if(i==nrpar-1)
      {
      if(increasing)
        sample = trunc_normal2(beta(nrpar-2,0),20,m,s);
      else
        sample = trunc_normal2(-20,beta(nrpar-2,0),m,s);
      }
    else
      {
      if(increasing)
        sample = trunc_normal2(beta(i-1,0),beta(i+1,0),m,s);
      else
        sample = trunc_normal2(beta(i+1,0),beta(i-1,0),m,s);
      }
    }
  else
    {
    double betaold,oben,unten;
    bool monotonic = false;

    if(increasing)
      {
      if(i < nrpar-2)
        oben = beta(i+2,0);
      else
        oben = 20;
      if(i >= 2)
        unten = beta(i-2,0);
      else
        unten = -20;
      }
    else
      {
      if(i < nrpar-2)
        unten = beta(i+2,0);
      else
        unten = -20;
      if(i >= 2)
        oben = beta(i-2,0);
      else
        oben = 20;
      }

    betaold = beta(i,0);

    while(monotonic == false)
      {

      sample = randnumbers::trunc_normal3(unten,oben,m,s);
      beta(i,0) = sample;

      if(i==0)
        {
        if(is_monotonic(0))
          monotonic = true;
        }
      else if(i==1)
        {
        if(is_monotonic(0) && is_monotonic(1))
          monotonic = true;
        }
      else if(i==2)
        {
        if(is_monotonic(0) && is_monotonic(1) && is_monotonic(2))
          monotonic = true;
        }
      else if(i==nrpar-1)
        {
        if(is_monotonic(nrpar-4))
          monotonic = true;
        }
     else if(i==nrpar-2)
        {
        if(is_monotonic(nrpar-4) && is_monotonic(nrpar-5))
          monotonic = true;
        }
      else if(i==nrpar-3)
        {
        if(is_monotonic(nrpar-4) && is_monotonic(nrpar-5) && is_monotonic(nrpar-6))
          monotonic = true;
        }
      else
        {
        if(is_monotonic(i-3) && is_monotonic(i-2) && is_monotonic(i-1) && is_monotonic(i))
          monotonic = true;
        }

      if(sample < betaold)
        unten = sample;
      else
        oben = sample;

      }

    }

  return sample;
  }


bool spline_basis::is_monotonic(const unsigned i)
  {
  double a,c,x,deriv;
  double h = knot[degree+1]-knot[degree];

  a = - beta(i,0) + 3*beta(i+1,0) - 3*beta(i+2,0) + beta(i+3,0);
  c = beta(i,0)*knot[degree+i+1] - beta(i+1,0)*(3*knot[degree+i+1] - h)
      + beta(i+2,0)*(h + 3*knot[degree+i]) - beta(i+3,0)*knot[degree+i];

  x = -c/a;

  if(a!=0 && knot[degree+i] < x && x < knot[degree+i+1])
    {
    deriv = deriv_f(x);
    if(increasing)
      {
      if(deriv >= -0.0000000001)
        return true;
      else
        return false;
      }
    else
      {
      if(deriv <= 0.0000000001)
        return true;
      else
        return false;
      }
    }
  else
    {
    return true;
    }

  }


bool spline_basis::is_monotonic(void)
  {
  bool monotonic = true;

  datamatrix help = spline;
  multBS(help,beta);

  unsigned i=1;
  if(increasing)
    {
    while(i<likep->get_nrobs() && monotonic==true)
      {
      if(help(i,0)+0.001 < help(i-1,0))
        monotonic = false;
      i++;
      }
    }
  else
    {
    while(i<likep->get_nrobs() && monotonic==true)
      {
      if(help(i,0) > help(i-1,0)+0.001)
        monotonic = false;
      i++;
      }
    }

  return monotonic;

  }


void spline_basis::multBS(datamatrix & res, const datamatrix & beta)
  {

  double *workres;
  double *workbeta;
  double *workBS;

  vector<int>::iterator freqwork = freq.begin();

  if(varcoeff)
    workBS = B.getV();
  else
    workBS = BS.getV();

  unsigned col = degree+1;
  unsigned j,k;
  int i,stop;
// res = 0 setzen
  workres = res.getV();
  for(j=0;j<res.rows()*res.cols();j++,workres++)
    *workres = 0.0;
// spaltenweise multiplizieren
  i = 0;
  k = 0;
  workres = res.getV();
  while (k<nrpar)
    {
    stop = lastnonzero[k];
//    while (i<lastnonzero[k]+1)
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


void spline_basis::multBS_index(datamatrix & res, const datamatrix & beta)
  {

  double *workres;
  double *workbeta;
  double *workBS;
  int *workindex;

  vector<int>::iterator freqwork = freq.begin();

  if(varcoeff)
    workBS = B.getV();
  else
    workBS = BS.getV();

  unsigned col = degree+1;
  unsigned j,k;
  int i,stop;
// res = 0 setzen
  workres = res.getV();
  for(j=0;j<res.rows()*res.cols();j++,workres++)
    *workres = 0.0;
// spaltenweise multiplizieren
  i = 0;
  k = 0;
  workres = res.getV();
  workindex = index.getV();
  while (k<nrpar)
    {
    stop = lastnonzero[k];
//    while (i<lastnonzero[k]+1)
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



void spline_basis::multDG(datamatrix & res, const datamatrix & b)
  {

  int i;
  unsigned j;
  double *workres;
  double *workDG;

  workres = res.getV();
  for(j=0;j<res.rows()*res.cols();j++,workres++)
    *workres = 0.0;

  workres = res.getV();
  workDG = DG.getV();
  for(i=0;i<gridsize;i++,workres++)
    for(j=0;j<degree+1;j++,workDG++)
      *workres += *workDG * b(DGfirst[i] + j,0);

  }


void spline_basis::update_prediction(void)
  {

  int j;

  if(predictright)
    {
    if(type==RW1)
      for(j=nrpar-nrparpredictright;j<int(nrpar);j++)
        beta(j,0) += beta(j-1,0);
    else if(type==RW2)
      for(j=nrpar-nrparpredictright;j<int(nrpar);j++)
        beta(j,0) += 2*beta(j-1,0) - beta(j-2,0);
    }
  if(predictleft)
    {
    if(type==RW1)
      for(j=nrparpredictleft-1;j>=0;j--)
        beta(j,0) += beta(j+1,0);
    else if(type==RW2)
      for(j=nrparpredictleft-1;j>=0;j--)
        beta(j,0) += 2*beta(j+1,0) - beta(j+2,0);
    }

  }


void spline_basis::outoptions(void)
  {
  ST::string typestr;
  ST::string knotstr;

  if (type == RW1)
    typestr = "first order random walk";
  else if (type == RW2)
    typestr = "second order random walk";
  else if (type == seasonal)
    typestr = "seasonal component";
  else if (type==mrf)
    typestr = "spatial Markov random field";
  else if (type==mrfkronecker)
    typestr= "Kronecker product interaction";
  else if (type==mrflinear)
    typestr = "2 dimensional first order random walk";
  else if (type==mrfkr1)
    typestr = "Kronecker product interaction (RW1*RW1)";
  else if (type==mrfkr2)
    typestr = "Kronecker product interaction (RW2*RW2)";
  else if (type == smoothspline)
    typestr = "Smoothing Splines";
  else if (type == npspline)
    typestr = "Natural P-Splines";

  if (knpos == equidistant)
    knotstr = "equidistant";
  else if (knpos == quantiles)
    knotstr = "quantiles";
  if (type == smoothspline)
    knotstr = "data points";

  optionsp->out("\n");
  optionsp->out("  Prior: " + typestr + "\n");
  optionsp->out("  Number of knots: " + ST::inttostring(nrknots) + "\n" );
  optionsp->out("  Knot choice: " + knotstr + "\n");
  optionsp->out("  Degree of Splines: " + ST::inttostring(degree) + "\n" );
  optionsp->out("\n");

  }


void spline_basis::outresults(void)
  {

  ST::string pathderiv = pathcurrent.substr(0,pathcurrent.length()-4)+"_derivative.res";


  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');


  if (transformnonlinear)
    {
    if (transformtype=="elasticity")
      {
      if (derivative)
        {
        ST::string suffix = "";
        fchelp.set_transform(suffix,transformtype);
        fcderivative.set_transform(suffix,transformtype);
        vector<FULLCOND*> fcvec(2);
        fcvec[0] = &fchelp;
        fcvec[1] = &fcderivative;
        likep->transform_nonlinear(fcvec,transformtype);
        pathderiv = pathcurrent;
        }
      else
        {
        optionsp->out("  Results for elasticities could not be computed because of\n");
        optionsp->out("  missing derivatives\n");
        optionsp->out("\n");
        }
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


  if ( !((transformnonlinear) && (transformtype=="elasticity") ) )
    {

    optionsp->out("  Results are stored in file\n");
    optionsp->out("  " + pathcurrent + "\n");
    optionsp->out("\n");
    #if defined(JAVA_OUTPUT_WINDOW)
    optionsp->out("  Postscript file is stored in file\n");
    ST::string psfile = pathcurrent.substr(0,pathcurrent.length()-4) + ".ps";
    optionsp->out("  " + psfile + "\n");
    optionsp->out("\n");
    optionsp->out("  Results may be visualized using method 'plotnonp'\n");
    optionsp->out("  Type for example: objectname.plotnonp " + ST::inttostring(fcnumber) + "\n");
    #else
    optionsp->out("  Results may be visualized using the R function 'plotnonp'\n");
    ST::string doublebackslash = "/";
    ST::string spluspath = pathcurrent.insert_string_char('\\',doublebackslash);
    optionsp->out("  Type for example:\n");
    optionsp->out("  plotnonp(\"" + spluspath + "\")");
    optionsp->out("\n");
    #endif
    optionsp->out("\n");


    fchelp.outresults();

    unsigned i;

    ofstream outres(pathcurrent.strtochar());

    outres << "intnr" << "   ";
    outres << datanames[0] << "   ";
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
    double * workbetaqu50 = fchelp.get_betaqu50p();
    double * workbetaqu_l1_upper_p = fchelp.get_beta_upper1_p();
    double * workbetaqu_l2_upper_p = fchelp.get_beta_upper2_p();
    double * workxvalues = xvalues.getV();

    for(i=0;i<xvalues.rows();i++,workmean++,
                     workbetaqu_l1_lower_p++,workbetaqu_l2_lower_p++,
                     workbetaqu50++,
                     workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                     workxvalues++)
      {
      outres << (i+1) << "   ";
      outres << *workxvalues << "   ";
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

    }

//---------------------------- Test-Ergebnisse ausgeben ------------------------

  bool test = false;
  if(test)          // Globale Option abfragen!!!
    {
    vector<ST::string> out;
//  ST::string pathref = make_pathref();    // Stefan!
    ST::string pathref = "c:\\cprog\\test\\results\\gaussian_nonp_1fkt_f_x_pspline_ref.res";
    compare_nonp(pathref,pathcurrent,0.02,out);
    for(unsigned i=0;i<out.size();i++)
      optionsp->out(out[i]);
    }

//------------------------------ Ableitung ausgeben ----------------------------

  if(derivative && ( !((transformnonlinear) && (transformtype!="elasticity"))  ) )
    {

    if ( (transformnonlinear) && (transformtype=="elasticity") )
      {
      optionsp->out("  Results for elasticities are stored in file\n");
      }
    else
      {
      optionsp->out("  Results for derivatives are stored in file\n");
      }

    optionsp->out("  " + pathderiv + "\n");
    optionsp->out("\n");
    #if defined(JAVA_OUTPUT_WINDOW)
    optionsp->out("  Postscript file is stored in file\n");
    ST::string psfile = pathderiv.substr(0,pathcurrent.length()-4) + ".ps";
    optionsp->out("  " + psfile + "\n");
    optionsp->out("\n");
    optionsp->out("  Results may be visualized using method 'plotnonp'\n");
    optionsp->out("  Type for example: objectname.plotnonp " + ST::inttostring(fcnumber) + "\n");
    #else
    optionsp->out("  Results may be visualized using the R function 'plotnonp'\n");
    ST::string doublebackslash = "/";
    ST::string spluspath = pathderiv.insert_string_char('\\',doublebackslash);
    optionsp->out("  Type for example:\n");
    optionsp->out("  plotnonp(\"" + spluspath + "\")");
    optionsp->out("\n");
    #endif
    optionsp->out("\n");

    fcderivative.outresults();

    ofstream outres2(pathderiv.strtochar());

    l1 = ST::doubletostring(lower1,4);
    l2 = ST::doubletostring(lower2,4);
    u1 = ST::doubletostring(upper1,4);
    u2 = ST::doubletostring(upper2,4);
    l1 = l1.replaceallsigns('.','p');
    l2 = l2.replaceallsigns('.','p');
    u1 = u1.replaceallsigns('.','p');
    u2 = u2.replaceallsigns('.','p');

    outres2 << "intnr" << "   ";
    outres2 << datanames[0] << "   ";
    outres2 << "pmean   ";
    outres2 << "pqu"  << l1  << "   ";
    outres2 << "pqu"  << l2  << "   ";
    outres2 << "pmed   ";
    outres2 << "pqu"  << u1  << "   ";
    outres2 << "pqu"  << u2  << "   ";
    outres2 << "pcat" << level1 << "   ";
    outres2 << "pcat" << level2 << "   ";

    outres2 << endl;

    double * workmean = fcderivative.get_betameanp();
    double * workbetaqu_l1_lower_p = fcderivative.get_beta_lower1_p();
    double * workbetaqu_l2_lower_p = fcderivative.get_beta_lower2_p();
    double * workbetaqu50 = fcderivative.get_betaqu50p();
    double * workbetaqu_l1_upper_p = fcderivative.get_beta_upper1_p();
    double * workbetaqu_l2_upper_p = fcderivative.get_beta_upper2_p();
    double * workxvalues = xvalues.getV();

    unsigned i;
    for(i=0;i<xvalues.rows();i++,workmean++,
                   workbetaqu_l1_lower_p++,workbetaqu_l2_lower_p++,
                   workbetaqu50++,
                   workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                   workxvalues++)
      {
      outres2 << (i+1) << "   ";
      outres2 << *workxvalues << "   ";
      outres2 << *workmean << "   ";
      outres2 << *workbetaqu_l1_lower_p << "   ";
      outres2 << *workbetaqu_l2_lower_p << "   ";
      outres2 << *workbetaqu50 << "   ";
      outres2 << *workbetaqu_l2_upper_p << "   ";
      outres2 << *workbetaqu_l1_upper_p << "   ";

      if (*workbetaqu_l1_lower_p > 0)
        outres2 << 1 << "   ";
      else if (*workbetaqu_l1_upper_p < 0)
        outres2 << -1 << "   ";
      else
        outres2 << 0 << "   ";

      if (*workbetaqu_l2_lower_p > 0)
        outres2 << 1 << "   ";
      else if (*workbetaqu_l2_upper_p < 0)
        outres2 << -1 << "   ";
      else
        outres2 << 0 << "   ";

      outres2 << endl;
      }

    }   // END: if(derivative)

// ------------------------ B-Spline Basis Funktionen ausgeben -----------------

//  outbsplines = true;
  if(outbsplines)
    {

    ST::string pathbsplines = pathcurrent.substr(0,pathcurrent.length()-4)
    +"_bsplines.res";
    datamatrix bsplines;
    write_bsplinefunctions(betamean,bsplines);

    ofstream outb(pathbsplines.strtochar());
    unsigned i,j;

    outb << datanames[0] << "   ";
    for (i=0;i<bsplines.cols();i++)
      outb << "B" << (i+1) << "   ";

    outb << endl;

    for(i=0;i<bsplines.rows();i++)
      {
      outb << xvalues(i,0) << "   ";
      for(j=0;j<bsplines.cols();j++)
        outb << bsplines(i,j) << "   ";

      outb << endl;
      }

    }

  }


void spline_basis::compute_XWX(const datamatrix & weight)
  {

/*
  unsigned i,j,k,l;
  vector<int>::iterator freqwork = freq.begin();
  int *workindex = index.getV();
  datamatrix xh1(nrpar,1,0);
  datamatrix xh2(nrpar,degree,0);

  for(i=0;i<nrknots-1;i++)
    {
    for(j=0;j<degree+1;j++)
      {
      for(k=0;k<degree+1;k++)
        {
        if(j+i<=k+i)
          {
          l = firstnonzero[i+degree];
          freqwork = freq.begin() + l;
          workindex = index.getV() + l;
          while(l<lastnonzero[i]+1)
            {
            if(j+i == k+i)
              xh1(j+i,0) += BS(*freqwork,j)*BS(*freqwork,k)*weight(*workindex,0);
            else
              xh2(j+i,k-j-1) += BS(*freqwork,j)*BS(*freqwork,k)*weight(*workindex,0);
            l++;
            freqwork++;
            workindex++;
            }
          }
        }
      }
    }

  XX.assign(xh1,xh2);
*/

  unsigned i,j,k,l;
  unsigned stop;
  unsigned BScols = degree+1;

  double *upper = XX.getupperpointer();
  double *diag = XX.getdiagpointer();
  double *workBS;
  double *workweight;
  double *workupper;
  double *workdiag;

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator workindex2 = index2.begin();

  deque<int>::iterator firstit = firstnonzero.begin();
  deque<int>::iterator lastit = lastnonzero.begin();

  for(i=0;i<nrpar;i++,diag++)
    {
    *diag = 0.0;
    for(j=0;j<degree;j++,upper++)
      *upper = 0.0;
    }

  upper = XX.getupperpointer();
  diag = XX.getdiagpointer();

  firstit += degree;
  for(i=0;i<nrknots-1;i++,diag++,upper+=degree,++firstit,++lastit)
    {
    stop = *lastit;
    for(j=0;j<BScols;j++)
      {
      workdiag = diag + j;
      for(k=j;k<BScols;k++)
        {

        l = *firstit;
        freqwork = freq.begin() + l;
        workindex2 = index2.begin() + l;
        workweight = weight.getV() + index(l,0);
        workBS = BS.getV() + BScols**freqwork;
        workupper = upper+j*degree+k-j-1;
        while(l <= stop)
          {
          if(j == k)
            *workdiag += *(workBS+j) * *workweight * *(workBS+k);
          else
            *workupper += *(workBS+j) * *workweight * *(workBS+k);
          l++;
          freqwork++;
          workBS += BScols*(*freqwork-*(freqwork-1));
          workindex2++;
          workweight += *workindex2;
          }

        }
      }
    }

  XX.set_decomposed();

  }


void spline_basis::compute_XWXenv(const datamatrix & weight, const unsigned & c)
  {

  unsigned i,j,k;
  int stop, l;
  unsigned BScols = degree+1;
  int weightcols = weight.cols();

  vector<double>::iterator diag = XX_env.getDiagIterator();
  vector<double>::iterator env = XX_env.getEnvIterator();
  vector<unsigned>::iterator xenv = XX_env.getXenvIterator();
  vector<double>::iterator workupper;
  vector<double>::iterator workdiag;

  double *workBS;
  double *workweight;

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator workindex2 = index2.begin();

  deque<int>::iterator firstit = firstnonzero.begin();
  deque<int>::iterator lastit = lastnonzero.begin();

  for(i=0;i<nrpar;++i,++diag)
    *diag = 0.0;

  unsigned envsize = XX_env.getEnv().size();
  for(j=0;j<envsize;++j,++env)
    *env = 0.0;

  env = XX_env.getEnvIterator();
  diag = XX_env.getDiagIterator();

  firstit += degree;
  for(i=0;i<nrknots-1;i++,diag++,++firstit,++lastit,++xenv)
    {
    stop = *lastit;
    for(j=0;j<BScols;j++)
      {
      workdiag = diag + j;
      for(k=j;k<BScols;k++)
        {

        l = *firstit;
        if(l <= stop)
          {
          freqwork = freq.begin() + l;
          workindex2 = index2.begin() + l;
          workweight = weight.getV() + c + weightcols*index(l,0);
          workBS = BS.getV() + BScols**freqwork;
          workupper = XX_env.getEnvIterator() + ( *(xenv+k+1) - (k-j) );
//          workupper = XX_env.getEnvIterator() + ( *(xenv+i+k+1) - (k-j) );
//          workupper = XX_env.getEnvIterator() + ( XX_env.getXenv(i+k+1) - (k-j) );
            while(l <= stop)
            {
            if(j == k)
              *workdiag += *(workBS+j) * *workweight * *(workBS+k);
            else
              *workupper += *(workBS+j) * *workweight * *(workBS+k);
            l++;
            if ( l >= stop ) break;
            freqwork++;
            workBS += BScols*(*freqwork-*(freqwork-1));
            workindex2++;
            workweight += *workindex2*weightcols;
            }

          }

        }
      }
    }

  XX_env.setDecomposed(false);

  }


void spline_basis::compute_XWXenv_XWtildey(const datamatrix & weight, const double & scale, const unsigned & c)
  {

  unsigned i,j,k;
  int stop, l;
  unsigned ind;
  unsigned BScols = degree+1;
  int weightcols = weight.cols();

  vector<double>::iterator diag = XX_env.getDiagIterator();
  vector<double>::iterator env = XX_env.getEnvIterator();
  vector<unsigned>::iterator xenv = XX_env.getXenvIterator();
  vector<double>::iterator workupper;
  vector<double>::iterator workdiag;

  double *workmu=NULL;
  double *workmuy=NULL;
  double *workmuy2=NULL;

  double *workBS=NULL;
  double *workweight=NULL;

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator workindex2 = index2.begin();

  deque<int>::iterator firstit = firstnonzero.begin();
  deque<int>::iterator lastit = lastnonzero.begin();

  workmuy = muy.getV();
  for(i=0;i<nrpar;i++,++workmuy)
    *workmuy = 0.0;

  for(i=0;i<nrpar;++i,++diag)
    *diag = 0.0;

  unsigned envsize = XX_env.getEnv().size();
  for(j=0;j<envsize;++j,++env)
    *env = 0.0;

  env = XX_env.getEnvIterator();
  diag = XX_env.getDiagIterator();

  firstit += degree;
  workmuy = muy.getV();
  for(i=0;i<nrknots-1;i++,++workmuy,diag++,++firstit,++lastit,++xenv)
    {
    stop = *lastit;
    for(j=0;j<BScols;j++)
      {
      workdiag = diag + j;
      for(k=j;k<BScols;k++)
        {
        l = *firstit;
        if(l <= stop)
          {
          ind = index(l,0);
          freqwork = freq.begin() + l;
          workindex2 = index2.begin() + l;
          workweight = weight.getV() + c + weightcols*ind;
          if(j == k)
            {
            workmu = mu.getV() + ind;
            workmuy2 = workmuy + j;
            }
          workBS = BS.getV() + BScols**freqwork;
          workupper = XX_env.getEnvIterator() + ( *(xenv+k+1) - (k-j) );
          while(l <= stop)
            {
            if(j == k)
              {
              *workmuy2 += *(workBS+j) * *workweight * *workmu;
              *workdiag += *(workBS+j) * *workweight * *(workBS+k);
              }
            else
              {
              *workupper += *(workBS+j) * *workweight * *(workBS+k);
              }
            l++;
            freqwork++;
            workBS += BScols*(*freqwork-*(freqwork-1));
            workindex2++;
            workweight += *workindex2*weightcols;
            workmu += *workindex2;
            }
          }
        }
      }
    }

  XX_env.setDecomposed(false);

  workmuy = muy.getV();
  for(i=0;i<nrpar;i++,++workmuy)
    *workmuy *= scale;

  }


void spline_basis::compute_XWtildey(const datamatrix & weight, const double & scale)
  {

/*
  unsigned i,j,l;
  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();

  for(i=0;i<nrpar;i++)
    muy(i,0) = 0.0;

  for(i=0;i<nrknots-1;i++)
    {
    for(j=0;j<degree+1;j++)
      {
      l = firstnonzero[i+degree];
      freqwork = freq.begin() + l;
      workindex = index.getV() + l;
      while(l<lastnonzero[i]+1)
        {
        muy(j+i,0) += weight(*workindex,0)*BS(*freqwork,j)*mu(*workindex,0);
        l++;
        freqwork++;
        workindex++;
        }
      }
    }

  for(i=0;i<nrpar;i++)
    muy(i,0) *= scale;
*/

  unsigned i,j;
  int stop, l;
  unsigned ind;
  unsigned BScols = degree+1;

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator workindex2 = index2.begin();

  deque<int>::iterator firstit = firstnonzero.begin();
  deque<int>::iterator lastit = lastnonzero.begin();

  double * workmu;
  double * workmuy;
  double * workmuy2;
  double * workweight;
  double * workBS;

  workmuy = muy.getV();
  for(i=0;i<nrpar;i++,++workmuy)
    *workmuy = 0.0;

  firstit += degree;
  workmuy = muy.getV();
  for(i=0;i<nrknots-1;i++,++workmuy,++firstit,++lastit)
    {
    stop = *lastit;
    for(j=0;j<BScols;j++)
      {
      l = *firstit;
      if(l<=stop)
        {
        ind = index(l,0);
        freqwork = freq.begin() + l;
        workindex2 = index2.begin() + l;
        workweight = weight.getV() + ind;
        workmu = mu.getV() + ind;
        workBS = BS.getV() + BScols**freqwork + j;
        workmuy2 = workmuy + j;
        while(l <= stop)
          {
          *workmuy2 += *workBS * *workweight * *workmu;
          l++;
          freqwork++;
          workindex2++;
          workweight += *workindex2;
          workmu += *workindex2;
          workBS += BScols*(*freqwork - *(freqwork-1));
          }
        }
      }
    }

  workmuy = muy.getV();
  for(i=0;i<nrpar;i++,++workmuy)
    *workmuy *= scale;

  }


void spline_basis::compute_XWtildey(const datamatrix & weight, const datamatrix & tildey, const double & scale, const unsigned & c)
  {

  unsigned i,j;
  int stop, l;
  unsigned ind;
  unsigned BScols = degree+1;
  int weightcols = weight.cols();

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator workindex2 = index2.begin();

  deque<int>::iterator firstit = firstnonzero.begin();
  deque<int>::iterator lastit = lastnonzero.begin();

  double * workmu;
  double * workmuy;
  double * workmuy2;
  double * workweight;
  double * workBS;

  workmuy = muy.getV();
  for(i=0;i<nrpar;i++,++workmuy)
    *workmuy = 0.0;

  firstit += degree;
  workmuy = muy.getV();
  for(i=0;i<nrknots-1;i++,++workmuy,++firstit,++lastit)
    {
    stop = *lastit;
    for(j=0;j<BScols;j++)
      {
      l = *firstit;
      if(l<=stop)
        {
        ind = index(l,0);
        freqwork = freq.begin() + l;
        workindex2 = index2.begin() + l;
        workweight = weight.getV() + c + weightcols*ind;
        workmu = tildey.getV() + ind;
        workBS = BS.getV() + BScols**freqwork + j;
        workmuy2 = workmuy + j;
        while(l <= stop)
          {
          *workmuy2 += *workBS * *workweight * *workmu;
          l++;
          freqwork++;
          workindex2++;
          workweight += *workindex2*weightcols;
          workmu += *workindex2;
          workBS += BScols*(*freqwork - *(freqwork-1));
          }
        }
      }
    }

  workmuy = muy.getV();
  for(i=0;i<nrpar;i++,++workmuy)
    *workmuy *= scale;

  }


void spline_basis::change(const datamatrix & main,const double & inter)
  {
  unsigned i;

  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();

// spline ändern
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    spline(*workindex,0) += main(*freqwork,0);

// Intercept ändern
  intercept += inter;

// fchelp ändern
  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
    {

    double * fchelpbetap = fchelp.getbetapointer();

    freqwork = freq.begin();
    workindex = index.getV();
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      {
      if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
        {
        *fchelpbetap = spline(*workindex,0) - intercept;
        fchelpbetap++;
        }
      }

    write_derivative();
    }
// update
  if(derivative)
    fcderivative.update();

  fchelp.update();
  FULLCOND::update();

  }


void spline_basis::change(const datamatrix & main)
  {

  unsigned i;
// beta, spline ändern
  beta.plus(beta,main);
  multBS_index(spline,beta);
  compute_intercept();
  for(i=0;i<nrpar;i++)
    beta(i,0) -= intercept;
  for(i=0;i<likep->get_nrobs();i++)
    spline(i,0) -= intercept;
  betaold.assign(beta);
  intercept = 0.0;

// fchelp ändern
  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
    {
    write_spline();
    write_derivative();
    }

  if(derivative)
    fcderivative.update();

  fchelp.update();
  FULLCOND::update();

  }


bool spline_basis::changeposterior(const datamatrix & main,const double & inter)
  {
  unsigned i;

  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();

// spline ändern
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    spline(*workindex,0) += main(*freqwork,0);

// Intercept ändern
  intercept += inter;

// fchelp ändern
  double * fchelpbetap = fchelp.getbetapointer();

  freqwork = freq.begin();
  workindex = index.getV();
  for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
    {
    if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
      {
      *fchelpbetap = spline(*workindex,0) - intercept;
      fchelpbetap++;
      }
    }

  write_derivative();
// posteriormode
  if(derivative)
    fcderivative.posteriormode();

  fchelp.posteriormode();
  return FULLCOND_nonp_basis::posteriormode();

  }


bool spline_basis::changeposterior(const datamatrix & main)
  {

  unsigned i;
// beta, spline ändern
  beta.plus(beta,main);
  multBS_index(spline,beta);
  compute_intercept();
  for(i=0;i<nrpar;i++)
    beta(i,0) -= intercept;
  for(i=0;i<likep->get_nrobs();i++)
    spline(i,0) -= intercept;
  betaold.assign(beta);
  intercept = 0.0;

// fchelp ändern
  write_spline();
  write_derivative();

  if(derivative)
    fcderivative.posteriormode();

  fchelp.posteriormode();
  return FULLCOND_nonp_basis::posteriormode();

  }


void spline_basis::sample_centered(datamatrix & beta)
    {

    unsigned i;
    double help;
    double * v;
    double * work;

    prec.solve(betaweight,betahelp,0,0);

    v = betahelp.getV();
    work = betaweight.getV();

    help = 0.0;
    for(i=0;i<nrpar;++i,++work,++v)
      help += *work * *v;

    compute_intercept();
    help = intercept/help;

    v = betahelp.getV();
    work = beta.getV();

    for(i=0;i<nrpar;++i,++work,++v)
      *work -= *v * help;

    intercept = 0.0;

    }


void spline_basis::sample_centered_env(datamatrix & beta)
    {

    unsigned i;
    double help;
    double * v;
    double * work;

    prec_env.solve(betaweight,betahelp);

    v = betahelp.getV();
    work = betaweight.getV();

    help = 0.0;
    for(i=0;i<nrpar;++i,++work,++v)
      help += *work * *v;

    compute_intercept();
    help = intercept/help;

    v = betahelp.getV();
    work = beta.getV();

    for(i=0;i<nrpar;++i,++work,++v)
      *work -= *v * help;

    intercept = 0.0;

    }


void spline_basis::getX(datamatrix & X)
  {

  double *workres;
  double *workBS;

  vector<int>::iterator freqwork = freq.begin();

  if(varcoeff)
    workBS = B.getV();
  else
    workBS = BS.getV();

  int i;
  unsigned j,k;
  unsigned col = degree+1;

  workres = X.getV();
  for(j=0;j<X.rows()*X.cols();j++,workres++)
    *workres = 0.0;

  i = 0;
  k = 0;
  while (k<nrpar)
    {
    while (i<lastnonzero[k]+1)
      {
      for(j=0;j<col;j++,workBS++)
        X(i,k+j) = *workBS;
      if((freqwork+1)!=freq.end() && *freqwork==*(freqwork+1))
        workBS -= col;
      i++;
      freqwork++;
      }
    k++;
    }

  }


void spline_basis::write_spline(const datamatrix & b)
  {

  double * splinep = splinehelp.getV();
  double * fchelpbetap = fchelp.getbetapointer();

  if(gridsize < 0)
    {
    multBS(splinehelp,b);
    vector<int>::iterator freqwork = freqoutput.begin();
    for(unsigned i=0;i<likep->get_nrobs();i++,freqwork++,splinep++)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        *fchelpbetap = *splinep;
        fchelpbetap++;
        }
      }
    }
  else
    {
    multDG(splinehelp,b);
    for(int i=0;i<gridsize;i++,fchelpbetap++,splinep++)
      *fchelpbetap = *splinep;
    }

  }


void spline_basis::write_bsplinefunctions(
const datamatrix & beta,datamatrix & bsplines)
  {

  unsigned k;

  datamatrix b(beta.rows(),1,0);
  if (gridsize < 0)
    bsplines = datamatrix(likep->get_nrobs(),beta.rows());
  else
    bsplines = datamatrix(gridsize,beta.rows());


  for (k=0;k<beta.rows();k++)
    {

    b(k,0) = beta(k,0);
    if (k > 0)
      b(k-1,0) = 0;

    double * splinep = splinehelp.getV();


    if(gridsize < 0)
      {
      multBS(splinehelp,b);
      vector<int>::iterator freqwork = freqoutput.begin();
      for(unsigned i=0;i<likep->get_nrobs();i++,freqwork++,splinep++)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          bsplines(i,k) = *splinep;
          }
        }
      }
    else
      {
      multDG(splinehelp,b);
      for(int i=0;i<gridsize;i++,splinep++)
        bsplines(i,k) = *splinep;
      }

    }

  }


void spline_basis::write_spline(void)
  {
  write_spline(beta);
  }


double spline_basis::compute_df(void)
  {
  if(prec_env.getDim()==0)
    return -1.0;
  if(lambda != lambda_prec || likep->iwlsweights_constant() == false)
    prec_env.addto(XX_env,Kenv,1.0,lambda);
  invprec = envmatdouble(0,nrpar,prec_env.getBandwidth());
  prec_env.inverse_envelope(invprec);
  if(identifiable)
    return invprec.traceOfProduct(XX_env);
  else
    return invprec.traceOfProduct(XX_env)-1;
  }


double spline_basis::compute_df_eigen(void)
  {
  if(prec_env.getDim()==0)
    return -1.0;

  unsigned i,j;
  datamatrix L(nrpar,nrpar,0.0);
  datamatrix help(nrpar,nrpar,0.0);
  datamatrix ev(nrpar,1,0.0);

  XX_env.decomp();
  for(i=0;i<nrpar;i++)
    for(j=0;j<=i;j++)
      L(i,j) = XX_env.getL(i,j);

  L = L.inverse();
  K.mult(L.transposed(),help);
  help = L*help;

  eigen2(help,ev);

  double df=0.0;
  for(i=0;i<nrpar;i++)
    df += 1.0/(1.0+lambda*ev(i,0));

  if(identifiable)
    return df;
  else
    return df-1;
  }


void spline_basis::write_derivative(void)
  {

  if(derivative)
    {
    unsigned i;
    fcderivative.set_transform(transform);
    Bderivative.mult(splinederivative,beta);
    double * fcderivbetap = fcderivative.getbetapointer();
    double * splinederivp = splinederivative.getV();
    for(i=0;i<splinederivative.rows();i++,fcderivbetap++,splinederivp++)
      *fcderivbetap = *splinederivp;
    } // END: if(derivative)

  }


void spline_basis::compute_Kweights(void)
  {

/*  if(type == RW1)
    weight = vector<double>(nrpar,1.0);
  else if(type == RW2)
    weight = vector<double>(nrpar,0.5);*/

  unsigned i,j;

// compute weights
  weight.reserve(nrpar);
  if(knpos == quantiles)
    {
    for(i=0;i<nrpar;i++)
      weight.push_back(knot[i+degree+1]-knot[i]);
    }
  else if(knpos == equidistant)
    {
    weight = vector<double>(nrpar,1.0);
    }

  double sum = 0.0;
  for(j=1;j<weight.size();j++)
    sum += weight[j];
  sum = double(weight.size()-1)/sum;

  if (type == RW2)
    sum *= 0.5;

  for(j=1;j<weight.size();j++)
    weight[j] *= sum;

  }

void spline_basis::init_name(const ST::string & na)
  {

  FULLCOND::init_name(na);

  ST::string underscore = "\\_";
  ST::string helpname = na.insert_string_char('_',underscore);
  term_symbolic = "f_{" + helpname + "}(" + helpname + ")";
  priorassumptions.push_back("$" + term_symbolic + "$:");

  if(column > 0)
    //term_symbolic = term_symbolic + " (" + ST::inttostring(column+1) + ". response category)";
    priorassumptions.push_back("$" + term_symbolic
           + " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)$:");

  if(type==MCMC::RW1)
     priorassumptions.push_back("P-spline with first order random walk penalty");
  else if(type==MCMC::RW2)
     priorassumptions.push_back("P-spline with second order random walk penalty");

  ST::string knotstr;
  if (knpos == equidistant)
    knotstr = "equidistant";
  else if (knpos == quantiles)
    knotstr = "quantiles";
  priorassumptions.push_back("Number of knots: " + ST::inttostring(nrknots));
  priorassumptions.push_back("Knot choice: " + knotstr);
  priorassumptions.push_back("Degree of Splines: " + ST::inttostring(degree));
  }


void spline_basis::init_names(const vector<ST::string> & na)
  {

  FULLCOND::init_names(na);

  ST::string underscore = "\\_";
  ST::string helpname0 = na[0].insert_string_char('_',underscore);
  ST::string helpname1 = na[1].insert_string_char('_',underscore);
  term_symbolic = "f_{" + helpname0 + "}(" + helpname0 + ")" + " \\cdot " + helpname1;
  priorassumptions.push_back("$" + term_symbolic + "$");

  if(column > 0)
    //term_symbolic = term_symbolic + " (" + ST::inttostring(column+1) + ". response category)";
    priorassumptions.push_back("$" + term_symbolic
           + " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)$:");

  if(type==MCMC::RW1)
     priorassumptions.push_back("P-spline with first order random walk penalty");
  else if(type==MCMC::RW2)
     priorassumptions.push_back("P-spline with second order random walk penalty");

  ST::string knotstr;
  if (knpos == equidistant)
    knotstr = "equidistant";
  else if (knpos == quantiles)
    knotstr = "quantiles";
  priorassumptions.push_back("Number of knots: " + ST::inttostring(nrknots));
  priorassumptions.push_back("Knot choice: " + knotstr);
  priorassumptions.push_back("Degree of Splines: " + ST::inttostring(degree));
  }


void spline_basis::set_lambdaconst(double la)
  {
  lambda=la;
  lambdaconst = true;
  }


void spline_basis::set_contour(int cp, bool pseudocp, bool app, int ls, const datamatrix & b)
  {
  contourprob = cp;
  lengthstart = ls;
  approx = app;
  pseudocontourprob = pseudocp;

  datamatrix beta_0;
  if(b.rows() < nrpar)
    beta_0 = datamatrix(nrpar,1,0);
  else
    beta_0 = b;

// Initialisierung von fc_contour
  if(contourprob >= 0)
    {
    ST::string path = samplepath.substr(0,samplepath.length()-4)+"_contour.raw";
//    fc_contour = FULLCOND(optionsp,beta_0,title+"_contour",nrpar+2,1,path);
    fc_contour = FULLCOND(optionsp,beta_0,title+"_contour",2*nrpar+6,1,path);
    fc_contour.setflags(MCMC::norelchange | MCMC::nooutput);
    }

  }


void spline_basis::createreml(datamatrix & X,datamatrix & Z,
                                const unsigned & Xpos, const unsigned & Zpos)
  {
  unsigned i,j,k;

  double * workdata;
  double * workZ;

  datamatrix refhelp;
  if(refcheck)
    {
    refhelp = bspline(reference);
    if(!varcoeff)
      X_ref = datamatrix(1,1);
    else
      X_ref = datamatrix(1,2);
    }

// X berechen

  datamatrix spline2;
  if(gridsize>0)
    {
    spline2 = datamatrix(gridsize,1,0);
    }

  if(type == RW1)
    {
    if(varcoeff && !centervcm)
      {
      if(X.rows()<spline.rows())
      // category specific interacting var
        {
        unsigned nrcat2 = spline.rows()/X.rows()-1;
        unsigned nrobs=X.rows();
        for(i=0; i<X.rows(); i++)
          {
          for(j=0; j<nrcat2; j++)
            {
            X(i,Xpos+j) = data_forfixed(j*nrobs+i,0) - data_forfixed(nrcat2*nrobs+i,0);
            }
          }
        }
      else
      // no category interacting var
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
      }
    }
  else if(type == RW2)
    {
    double * workX;
    unsigned Xcols = X.cols();

    datamatrix knoten = datamatrix(nrpar,1,0.0);
    for(i=0;i<nrpar;i++)
      knoten(i,0) = knot[i];

    multBS_index(spline,knoten);

    double splinemean=spline.mean(0);
    if(!varcoeff)
      {
      for(i=0; i<spline.rows(); i++)
        {
        spline(i,0)=spline(i,0)-splinemean;
        }
      if(gridsize>0)
        {
        multDG(spline2,knoten);
        for(i=0; i<spline2.rows(); i++)
          {
          X_grid(i,0) = spline2(i,0)-splinemean;
          }
        }
      }

      if(refcheck)
        {
        if(!varcoeff)
          {
          for(i=0; i<knoten.rows(); i++)
            X_ref(0,0) += knoten(i,0)*refhelp(i,0);
          X_ref(0,0) -= splinemean;
          }
        else
          {
          X_ref(0,0) = 1.0;
          for(i=0; i<knoten.rows(); i++)
            X_ref(0,1) += knoten(i,0)*refhelp(i,0);
          }
        }

    if(X.rows()<spline.rows())
    // category specific covariates
      {
      unsigned nrcat2 = spline.rows()/X.rows()-1;
      unsigned nrobs=X.rows();
      if(varcoeff)
        {
        X_VCM.putCol(1,spline);

        for(i=0; i<X.rows(); i++)
          {
          for(j=0; j<nrcat2; j++)
            {
            X(i,Xpos+2*j) = data_forfixed(j*nrobs+i,0) - data_forfixed(nrcat2*nrobs+i,0);
            X(i,Xpos+2*j+1) = spline(i,0)*(data_forfixed(j*nrobs+i,0) - data_forfixed(nrcat2*nrobs+i,0));
            }
          }
        }
      else
        {
        X_VCM = spline;

        for(i=0; i<X.rows(); i++)
          {
          for(j=0; j<nrcat2; j++)
            {
            X(i,Xpos+j) = spline(j*nrobs+i,0) - spline(nrcat2*nrobs+i,0);
            }
          }
        }
      }
    else
    // no category specific covariates
      {
      workdata = spline.getV();
      workX = X.getV()+Xpos;

      if(varcoeff)
        {
        double * workintact = data_forfixed.getV();
        if(!centervcm)
          {
          double * workX_VCM = X_VCM.getV()+1;
          for (i=0;i<spline.rows();i++,workdata++,workintact++,workX+=Xcols,workX_VCM+=2)
            {
            *workX = *workintact;
            *(workX+1) = *workdata**workintact;
            *workX_VCM = *workdata;
            }
          }
        else
          {
          double * workX_VCM = X_VCM.getV();
          for (i=0;i<spline.rows();i++,workdata++,workintact++,workX+=Xcols,workX_VCM++)
            {
            *workX = *workdata**workintact;
            *workX_VCM = *workdata;
            }
          }
        }
      else
        {
        for (i=0;i<spline.rows();i++,workdata++,workX+=Xcols)
          {
          *workX = *workdata;
          }
        }
      }
    }

// Z berechnen

  compute_Kweights();
//  datamatrix diffmatrix = diffmat(dimX+1,nrpar);
  datamatrix diffmatrix = weighteddiffmat(type==RW1?1:2,weight);
  diffmatrix = diffmatrix.transposed()*diffmatrix.transposed().sscp().inverse();

  unsigned Zcols = Z.cols();

  if(Z.rows()<spline.rows())
    // category specific covariates
    {
    Z_VCM = datamatrix(spline.rows(), dimZ, 0);
    }

  if(refcheck)
    Z_ref = datamatrix(1,dimZ);

  for(j=0;j<dimZ;j++)
    {
    multBS_index(spline,diffmatrix.getCol(j));

    if(Z.rows()<spline.rows())
      // category specific covariates
      {
      unsigned nrcat2 = spline.rows()/Z.rows()-1;
      unsigned nrobs=Z.rows();
      if(varcoeff)
        {
        for(i=0; i<nrobs; i++)
          {
          for(k=0; k<nrcat2; k++)
            {
            Z(i, Zpos + k*dimZ + j) = spline(i,0)*(data_forfixed(k*nrobs+i,0) - data_forfixed(nrcat2*nrobs+i,0));
            }
          }
        for(i=0; i<spline.rows(); i++)
          {
          Z_VCM(i,j) = spline(i,0);
          }
        }
      else
        {
        for(i=0; i<nrobs; i++)
          {
          for(k=0; k<nrcat2; k++)
            {
            Z(i, Zpos + k*dimZ + j) = spline(k*nrobs+i,0) - spline(nrcat2*nrobs+i,0);
            }
          }
        for(i=0; i<spline.rows(); i++)
          {
          Z_VCM(i,j) = spline(i,0);
          }
        }
      }
    else
      // no category specific covariates
      {
      if(refcheck)
        {
        Z_ref(0,j) = (refhelp.transposed()*diffmatrix.getCol(j))(0,0);
        }

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
    }

  if(gridsize>0)
    {
    for(j=0;j<dimZ;j++)
      {
      multDG(spline2,diffmatrix.getCol(j));
      workdata = spline2.getV();
      workZ = Z_grid.getV()+j;
      for (i=0;i<spline2.rows();i++,workdata++,workZ+=dimZ)
        {
        *workZ = *workdata;
        }
      }
    }
  }


double spline_basis::outresultsreml(datamatrix & X,datamatrix & Z,
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

  unsigned nr = xvalues.rows();
  unsigned i,j;

  betamean=datamatrix(nr,1,0);
  datamatrix betastd=datamatrix(nr,1,0);
  betaqu_l1_lower=datamatrix(nr,1,0);
  betaqu_l1_upper=datamatrix(nr,1,0);
  betaqu_l2_lower=datamatrix(nr,1,0);
  betaqu_l2_upper=datamatrix(nr,1,0);

  datamatrix betarefmean;
  datamatrix betarefstd;
  datamatrix betarefqu_l1_lower;
  datamatrix betarefqu_l1_upper;
  datamatrix betarefqu_l2_lower;
  datamatrix betarefqu_l2_upper;

  if(refcheck)
    {
    betarefmean=datamatrix(nr,1,0);
    betarefstd=datamatrix(nr,1,0);
    betarefqu_l1_lower=datamatrix(nr,1,0);
    betarefqu_l1_upper=datamatrix(nr,1,0);
    betarefqu_l2_lower=datamatrix(nr,1,0);
    betarefqu_l2_upper=datamatrix(nr,1,0);
    }

  vector<int>::iterator indexit = index2.begin();
  unsigned k = *indexit;
  vector<int>::iterator freqwork = freqoutput.begin();

  if(varcoeff || X.rows()<X_VCM.rows())
    {
    for(i=0,j=0;i<Z_VCM.rows();i++,indexit++,freqwork++,k+=*indexit)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        if(type == RW1)
          {
          if(!centervcm)
            {
            betamean(j,0) = betareml(betaXpos,0)*X_VCM(k,0) + (Z_VCM.getRow(k)*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1))(0,0);
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
            }
          else
            {
            betamean(j,0) = (Z_VCM.getRow(k)*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1))(0,0);
            betastd(j,0) = sqrt(
                                (Z_VCM.getRow(k)*
                                 betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)*
                                 (Z_VCM.getRow(k).transposed()))(0,0));
            }
          }
        else
          {
          betamean(j,0) = (X_VCM.getRow(k)*betareml.getBlock(betaXpos,0,betaXpos+dimX,1))(0,0) + (Z_VCM.getRow(k)*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1))(0,0);
          betastd(j,0) = sqrt(
                              ((
                              X_VCM.getRow(k)*betacov.getBlock(betaXpos,betaXpos,betaXpos+dimX,betaXpos+dimX)
                              +
                              (Z_VCM.getRow(k)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+dimX))
                              )*X_VCM.getRow(k).transposed())(0,0)
                              +
                              (
                              (
                               X_VCM.getRow(k)*betacov.getBlock(betaXpos,betaZpos,betaXpos+dimX,betaZpos+dimZ)
                               +
                               Z_VCM.getRow(k)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                              )*(Z_VCM.getRow(k).transposed())
                              )(0,0)
                             );
          }
        j++;
        }
      }
    }
  else
    {
    for(i=0,j=0;i<X.rows();i++,indexit++,freqwork++,k+=*indexit)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        if(type == RW1)
          {
          betamean(j,0) = (Z.getBlock(k,Zpos,k+1,Zpos+nrpar-1)*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1))(0,0);
          betastd(j,0) = sqrt((Z.getBlock(k,Zpos,k+1,Zpos+nrpar-1)*
                   betacov.getBlock(betaZpos,betaZpos,betaZpos+nrpar-1,betaZpos+nrpar-1)*
                   Z.getBlock(k,Zpos,k+1,Zpos+nrpar-1).transposed())(0,0));
          }
        else
          {
          betamean(j,0) = betareml(betaXpos,0)*X(k,Xpos) + (Z.getBlock(k,Zpos,k+1,Zpos+dimZ)*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1))(0,0);
          betastd(j,0) = sqrt(
                              (
                               X(k,Xpos)*betacov(betaXpos,betaXpos)
                               +
                               (Z.getBlock(k,Zpos,k+1,Zpos+dimZ)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+1))(0,0)
                              )*X(k,Xpos)
                              +
                              (
                               (
                                X(k,Xpos)*betacov.getBlock(betaXpos,betaZpos,betaXpos+1,betaZpos+dimZ)
                                +
                                Z.getBlock(k,Zpos,k+1,Zpos+dimZ)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                               )*(Z.getBlock(k,Zpos,k+1,Zpos+dimZ).transposed())
                              )(0,0)
                             );
          }
        j++;
        }
      }
    }

  if(!varcoeff)
    {
    mean = betamean.mean(0);
//    betareml(0,0)=betareml(0,0)+mean;
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

  indexit = index2.begin();
  k = *indexit;
  freqwork = freqoutput.begin();

  if(refcheck)
    {
    if(varcoeff || X.rows()<X_VCM.rows())
      {
      for(i=0,j=0;i<Z_VCM.rows();i++,indexit++,freqwork++,k+=*indexit)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          if(type == RW1)
            {
            if(!centervcm)
              {
              betarefmean(j,0) = betareml(betaXpos,0)*(X_VCM(k,0)-X_ref(0,0)) + ((Z_VCM.getRow(k)-Z_ref)*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1))(0,0);
              betarefstd(j,0) = sqrt(
                                  (
                                   (X_VCM(k,0)-X_ref(0,0))*betacov(betaXpos,betaXpos)
                                   +
                                   ((Z_VCM.getRow(k)-Z_ref)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+1))(0,0)
                                  )*(X_VCM(k,0)-X_ref(0,0))
                                  +
                                  (
                                   (
                                    (X_VCM(k,0)-X_ref(0,0))*betacov.getBlock(betaXpos,betaZpos,betaXpos+1,betaZpos+dimZ)
                                    +
                                    (Z_VCM.getRow(k)-Z_ref)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                                   )*((Z_VCM.getRow(k)-Z_ref).transposed())
                                  )(0,0)
                                 );
              }
            else
              {
              betarefmean(j,0) = ((Z_VCM.getRow(k)-Z_ref)*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1))(0,0);
              betarefstd(j,0) = sqrt(
                                  ((Z_VCM.getRow(k)-Z_ref)*
                                   betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)*
                                   ((Z_VCM.getRow(k)-Z_ref).transposed()))(0,0));
              }
            }
          else
            {
            betarefmean(j,0) = ((X_VCM.getRow(k)-X_ref)*betareml.getBlock(betaXpos,0,betaXpos+dimX,1))(0,0) + ((Z_VCM.getRow(k)-Z_ref)*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1))(0,0);
            betarefstd(j,0) = sqrt(
                                ((
                                (X_VCM.getRow(k)-X_ref)*betacov.getBlock(betaXpos,betaXpos,betaXpos+dimX,betaXpos+dimX)
                                +
                                ((Z_VCM.getRow(k)-Z_ref)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+dimX))
                                )*(X_VCM.getRow(k)-X_ref).transposed())(0,0)
                                +
                                (
                                (
                                 (X_VCM.getRow(k)-X_ref)*betacov.getBlock(betaXpos,betaZpos,betaXpos+dimX,betaZpos+dimZ)
                                 +
                                 (Z_VCM.getRow(k)-Z_ref)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                                )*((Z_VCM.getRow(k)-Z_ref).transposed())
                                )(0,0)
                               );
            }
          j++;
          }
        }
      }
    else
      {
      for(i=0,j=0;i<X.rows();i++,indexit++,freqwork++,k+=*indexit)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          if(type == RW1)
            {
            betarefmean(j,0) = ((Z.getBlock(k,Zpos,k+1,Zpos+nrpar-1)-Z_ref)*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1))(0,0);
            betarefstd(j,0) = sqrt(((Z.getBlock(k,Zpos,k+1,Zpos+nrpar-1)-Z_ref)*
                     betacov.getBlock(betaZpos,betaZpos,betaZpos+nrpar-1,betaZpos+nrpar-1)*
                     (Z.getBlock(k,Zpos,k+1,Zpos+nrpar-1)-Z_ref).transposed())(0,0));
            }
          else
            {
            betarefmean(j,0) = betareml(betaXpos,0)*(X(k,Xpos)-X_ref(0,0)) + ((Z.getBlock(k,Zpos,k+1,Zpos+dimZ)-Z_ref)*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1))(0,0);
            betarefstd(j,0) = sqrt(
                                (
                                 (X(k,Xpos)-X_ref(0,0))*betacov(betaXpos,betaXpos)
                                 +
                                 ((Z.getBlock(k,Zpos,k+1,Zpos+dimZ)-Z_ref)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+1))(0,0)
                                )*(X(k,Xpos)-X_ref(0,0))
                                +
                                (
                                 (
                                  (X(k,Xpos)-X_ref(0,0))*betacov.getBlock(betaXpos,betaZpos,betaXpos+1,betaZpos+dimZ)
                                  +
                                  (Z.getBlock(k,Zpos,k+1,Zpos+dimZ)-Z_ref)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                                 )*((Z.getBlock(k,Zpos,k+1,Zpos+dimZ)-Z_ref).transposed())
                                )(0,0)
                               );
            }
          j++;
          }
        }
      }  
    for(j=0; j<nr; j++)
      {
//      betarefmean(j,0)=betarefmean(j,0)-mean;
      betarefqu_l1_lower(j,0) = betarefmean(j,0)+randnumbers::invPhi2(lower1/100)*betarefstd(j,0);
      betarefqu_l1_upper(j,0) = betarefmean(j,0)+randnumbers::invPhi2(upper2/100)*betarefstd(j,0);
      betarefqu_l2_lower(j,0) = betarefmean(j,0)+randnumbers::invPhi2(lower2/100)*betarefstd(j,0);
      betarefqu_l2_upper(j,0) = betarefmean(j,0)+randnumbers::invPhi2(upper1/100)*betarefstd(j,0);
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
  assert(!outres.fail());

  optionsp->out("\n");
  optionsp->out("  Results are stored in file\n");
  optionsp->out("  " + outest + "\n");
  optionsp->out("\n");
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

  outres << "intnr" << "   ";
  outres << datanames[0] << "   ";
  outres << "pmode   ";
  outres << "ci"  << level1  << "lower   ";
  outres << "ci"  << level2  << "lower   ";
  outres << "std   ";
  outres << "ci"  << level2  << "upper   ";
  outres << "ci"  << level1  << "upper   ";
  outres << "pcat" << level1 << "   ";
  outres << "pcat" << level2 << "   ";

  if(refcheck)
    {
    outres << "pmoderef   ";
    outres << "ci"  << level1  << "lowerref   ";
    outres << "ci"  << level2  << "lowerref   ";
    outres << "stdref   ";
    outres << "ci"  << level2  << "upperref   ";
    outres << "ci"  << level1  << "upperref   ";
    }

  outres << endl;

  double * workmean = betamean.getV();
  double * workstd = betastd.getV();
  double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
  double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
  double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
  double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
  double * workxvalues = xvalues.getV();

  double * workrefmean = betarefmean.getV();
  double * workrefstd = betarefstd.getV();
  double * workrefbetaqu_l1_lower_p = betarefqu_l1_lower.getV();
  double * workrefbetaqu_l2_lower_p = betarefqu_l2_lower.getV();
  double * workrefbetaqu_l1_upper_p = betarefqu_l1_upper.getV();
  double * workrefbetaqu_l2_upper_p = betarefqu_l2_upper.getV();

  for(i=0;i<xvalues.rows();i++,workmean++,workstd++,
                     workbetaqu_l1_lower_p++,workbetaqu_l2_lower_p++,
                     workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                     workxvalues++)
    {
    outres << (i+1) << "   ";
    outres << *workxvalues << "   ";
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

    if(refcheck)
      {
      outres << *workrefmean << "   ";
      outres << *workrefbetaqu_l1_lower_p << "   ";
      outres << *workrefbetaqu_l2_lower_p << "   ";
      outres << *workrefstd << "   ";
      outres << *workrefbetaqu_l2_upper_p << "   ";
      outres << *workrefbetaqu_l1_upper_p << "   ";

      workrefmean++;
      workrefstd++;
      workrefbetaqu_l1_lower_p++;
      workrefbetaqu_l2_lower_p++;
      workrefbetaqu_l1_upper_p++;
      workrefbetaqu_l2_upper_p++;
      }

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

    if(type == RW1)
      {
      betamean = Z_grid*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1);
      for(i=0; i<gridsize; i++)
        {
        betamean(i,0)=betamean(i,0)-mean;
        betastd(i,0) = sqrt((Z_grid.getRow(i)*
                 betacov.getBlock(betaZpos,betaZpos,betaZpos+nrpar-1,betaZpos+nrpar-1)*
                 Z_grid.getRow(i).transposed())(0,0));
        }
      }
    else
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

    assert(!outgrid.fail());

    outgrid << "intnr" << "   ";
    outgrid << datanames[0] << "   ";
    outgrid << "pmode   ";
    outgrid << "ci"  << level1  << "lower   ";
    outgrid << "ci"  << level2  << "lower   ";
    outgrid << "std   ";
    outgrid << "ci"  << level2  << "upper   ";
    outgrid << "ci"  << level1  << "upper   ";
    outgrid << "pcat" << level1 << "   ";
    outgrid << "pcat" << level2 << "   ";
    outgrid << endl;

    double * workmean = betamean.getV();
    double * workstd = betastd.getV();
    double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
    double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
    double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
    double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
    double * workxvalues = xvaluesgrid.getV();

    for(i=0;i<gridsize;i++,workmean++,workstd++,
                       workbetaqu_l1_lower_p++,workbetaqu_l2_lower_p++,
                       workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                       workxvalues++)
      {
      outgrid << (i+1) << "   ";
      outgrid << *workxvalues << "   ";
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


void spline_basis::outresultsgrid()
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


void spline_basis::outoptionsreml()
  {
  optionsp->out("OPTIONS FOR P-SPLINE TERM: " + title + "\n",true);
  optionsp->out("\n");

  ST::string typestr;
  ST::string knotstr;

  if (type == RW1)
    typestr = "first order random walk";
  else if (type == RW2)
    typestr = "second order random walk";

  if (knpos == equidistant)
    knotstr = "equidistant";
  else if (knpos == quantiles)
    knotstr = "quantiles";

  optionsp->out("  Prior: " + typestr + "\n");
  optionsp->out("  Number of knots: " + ST::inttostring(nrknots) + "\n" );
  optionsp->out("  Knot choice: " + knotstr + "\n");
  optionsp->out("  Degree of Splines: " + ST::inttostring(degree) + "\n" );
  optionsp->out("  Starting value for lambda: " + ST::doubletostring(startlambda,6) + "\n" );
  optionsp->out("\n");
  }

// ------------------------- FOR MERROR ----------------------------------------

void spline_basis::update_merror_varcoef(datamatrix & effmod, datamatrix & newintact)
  {
//  unsigned i,j,k;
//  double value;

  freq = vector<int>();
  freqoutput = vector<int>();
  begcol = vector<int>();
  index2 = vector<int>();

  firstnonzero = deque<int>();
  lastnonzero = deque<int>();
  knot = deque<double>();

  nrdiffobs=0;

  make_index(effmod,newintact);
  make_index2();
  make_Bspline(effmod);
  make_BS(newintact);
  }

void spline_basis::update_merror(datamatrix & newdata)
  {
  unsigned i,j,k;
  double value;

  freq = vector<int>();
  begcol = vector<int>();
  index2 = vector<int>();

  firstnonzero = deque<int>();
  lastnonzero = deque<int>();
  for(i=0;i<nrpar;i++)
    {
    lastnonzero.push_back(-1);
    firstnonzero.push_back(0);
    }

  make_index(newdata);
  make_index2();

  datamatrix help(nrpar,1,0.0);
  Bcolmean = datamatrix(nrpar,1,0.0);

  BS = datamatrix(nrdiffobs,degree+1,0.0);
  double * work = BS.getV();

  vector<int>::iterator freqwork = freq.begin();
  for(i=0;i<newdata.rows();i++,++freqwork)
    {
    value = newdata(index(i,0),0);
    if(freqwork == freq.begin() || *freqwork != *(freqwork-1))
      {
      j=0;
      while(knot[degree+j+1] <= value)
        j++;
      begcol.push_back(j);
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

void spline_basis::update_merror_discrete(datamatrix & newdata)
  {
  unsigned i;
  freq = vector<int>();
  index2 = vector<int>();

  make_index(newdata);
  make_index2();

  for(i=0; i<spline.rows(); i++)
    spline(i,0) -= intercept;
  }

datamatrix spline_basis::get_spline_merror(void)
  {
  datamatrix help(spline.rows(),1,0);
  multBS_index(help, beta);
  return help;
  }


// -------------------------END: FOR MERROR ------------------------------------


} // end: namespace MCMC

//---------------------------------------------------------------------------
#if !defined(__BUILDING_GNU)
#pragma package(smart_init)
#endif




