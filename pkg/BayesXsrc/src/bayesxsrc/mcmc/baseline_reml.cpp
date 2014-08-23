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



#include "baseline_reml.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//---------- CLASS: baseline_reml (implementation of member functions) ---------
//------------------------------------------------------------------------------

baseline_reml::baseline_reml(MCMCoptions * o,
              const datamatrix & d, const datamatrix & leftint,
              const datamatrix & lefttrunc, const unsigned & nrk,
              const unsigned & degr, const unsigned & tgr, const unsigned & nrq,
              const unsigned & nrb, const knotpos & kp, const fieldtype & ft,
              const ST::string & ti, const ST::string & fp,
              const ST::string & pres, const double & l, const double & sl,
              const knotpos & gp, const int & gs, const bool & catsp, const double & rv)
  : spline_basis(o,d,nrk,degr,kp,ft,ti,fp,pres,l,sl,catsp,0.0,0.0,0.0,0.0,gs,rv)
  {
  unsigned i,j,k;

  baseline=true;
  varcoeff=false;

//  refcheck=false;

  gridpos = gp;
  double tmax=d.max(0);

  if(gridpos == MCMC::equidistant)
    {
    tgrid = tgr;
    tvalues = datamatrix(tgrid+1,1);
    tstep = tmax / tgrid;
    for(i=0;i<tvalues.rows();i++)
      tvalues(i,0) = i*tstep;
    }
  else if(gridpos == MCMC::quantiles)
    {
    tgrid = nrq*nrb;
    tvalues = datamatrix(tgrid+1,1);
    nrquant = nrq;
    nrbetween = nrb;
    datamatrix tquantiles = datamatrix(nrquant+1,1,0);
    for(i=1; i<nrquant; i++)
      {
      tquantiles(i,0) = d.quantile(((double)i/nrquant)*100,0);
      }
    tquantiles(nrquant,0) = tmax;

    double intmax, intmin, intstep;
    for(i=0; i<nrquant; i++)
      {
      intmin=tquantiles(i,0);
      intmax=tquantiles(i+1,0);
      intstep=(intmax-intmin)/nrbetween;
      for(j=0; j<nrbetween; j++)
        {
        tvalues(i*nrbetween+j,0) = intmin + j*intstep;
        }
      }
    tvalues(tgrid,0) = tmax;
    }
  else
    {
    make_index(d);
    vector<int>::iterator freqwork = freqoutput.begin();
    int * workindex = index.getV();
    tvalues = datamatrix(nrdiffobs,1,0);
    for(j=0;j<d.rows();j++,freqwork++,workindex++)
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        tvalues(*freqwork,0) = d(*workindex,0);
    }

  tsteps = datamatrix(tvalues.rows()-1,1,0);
  for(i=0; i<tsteps.rows(); i++)
    tsteps(i,0) = tvalues(i+1,0)-tvalues(i,0);

  interact_var = datamatrix(d.rows(),1,1);

  datamatrix betahelp(nrpar,1,0);
  DG = datamatrix(tvalues.rows(),degree+1,0);
  DGfirst = vector<int>(tvalues.rows());

  for(i=0;i<tvalues.rows();i++)
    {
    betahelp.assign( bspline(tvalues(i,0)) );
    j=degree+1;
    while(knot[j] <= tvalues(i,0) && j<nrknots+degree)
      j++;
    for(k=0;k<degree+1;k++)
      DG(i,k) = betahelp(k+j-(degree+1),0);
    DGfirst[i] = j-(degree+1);
    }

  tleft = vector<unsigned>(d.rows(),0);
  tright = vector<unsigned>(d.rows(),1);
  ttrunc = vector<unsigned>(d.rows(),0);

  // indices for truncation times

  if(lefttrunc.rows()>1)
    {
    for(i=0; i<lefttrunc.rows(); i++)
      {
      j=1;
      while(j<tvalues.rows() && tvalues(j,0) < lefttrunc(i,0))
        {
        ttrunc[i]++;
        j++;
        }
      }
    }

  for(i=0; i<d.rows(); i++)
    {
    j=0;
    while(j<tvalues.rows()-2 && tvalues(j,0)<d(i,0))
      {
      tright[i]++;
      j++;
      }
    }

  if(leftint.rows()>1)
    {
    for(i=0; i<d.rows(); i++)
      {
      if( leftint(i,0) < d(i,0))
        {
        j=0;
        while(j<tvalues.rows()-1 && tvalues(j,0)<leftint(i,0))
          {
          tleft[i]++;
          j++;
          }
        }
      else if(leftint(i,0)> d(i,0))
        {
        tleft[i] = tright[i]+1;
        }
      else
        {
        tleft[i] = tright[i];
        }
      }
    }
  else
    {
    for(i=0; i<d.rows(); i++)
      {
      tleft[i] = tright[i];
      }
    }

  }

baseline_reml::baseline_reml(MCMCoptions * o,const datamatrix & d1,
                      const datamatrix & d2, const unsigned & nrk,
                      const unsigned & degr, const unsigned & tgr,
                      const knotpos & kp, const fieldtype & ft, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres, const double & l,
                      const double & sl, const int & gs, const bool & catsp, const double & rv)
  : spline_basis(o,d1,d2,nrk,degr,kp,ft,ti,fp,pres,l,sl,catsp,false,rv)
  {
  baseline=true;
  interact_var = d2;
  }

baseline_reml::baseline_reml(const baseline_reml & fc)
  :spline_basis(spline_basis(fc))
  {
  tstep=fc.tstep;
  tvalues=fc.tvalues;
  tleft=fc.tleft;
  tright=fc.tright;
  ttrunc=fc.ttrunc;
  t_X=fc.t_X;
  t_Z=fc.t_Z;
  interact_var=fc.interact_var;
  tgrid=fc.tgrid;
  tsteps=fc.tsteps;
  gridpos=fc.gridpos;
  nrquant=fc.nrquant;
  nrbetween=fc.nrbetween;
  }

const baseline_reml & baseline_reml::operator=(const baseline_reml & fc)
  {
  if (this == &fc)
    return *this;
  spline_basis::operator=(spline_basis(fc));

  tstep=fc.tstep;
  tvalues=fc.tvalues;
  tleft=fc.tleft;
  tright=fc.tright;
  ttrunc=fc.ttrunc;
  t_X=fc.t_X;
  t_Z=fc.t_Z;
  interact_var=fc.interact_var;
  tgrid=fc.tgrid;
  tsteps=fc.tsteps;
  gridpos=fc.gridpos;
  nrquant=fc.nrquant;
  nrbetween=fc.nrbetween;

  return *this;
  }

void baseline_reml::createreml(datamatrix & X,datamatrix & Z,
                                const unsigned & Xpos, const unsigned & Zpos)
  {
  unsigned i,j;

  double * workdata;
  double * workZ;
  double * workX;
  unsigned Xcols = X.cols();

  datamatrix refhelp;
  if(refcheck)
    {
    refhelp = bspline(reference);
    if(!varcoeff)
      X_ref = datamatrix(1,1);
    else
      X_ref = datamatrix(1,2);
    }

  // X für Daten berechen

  datamatrix knoten = datamatrix(nrpar,1,0.0);
  for(i=0;i<nrpar;i++)
    knoten(i,0) = knot[i];

  multBS_index(spline,knoten);

  workdata = spline.getV();
  workX = X.getV()+Xpos;
  if(varcoeff)
    {
    double * workintact = data_forfixed.getV();
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
    for (i=0;i<spline.rows();i++,workdata++,workX+=Xcols)
      {
      *workX = *workdata;
      }
    }

  if(refcheck)
    {
    if(!varcoeff)
      {
      for(i=0; i<knoten.rows(); i++)
        X_ref(0,0) += knoten(i,0)*refhelp(i,0);
//      X_ref(0,0) -= splinemean;
      }
    else
      {
      X_ref(0,0) = 1.0;
      for(i=0; i<knoten.rows(); i++)
        X_ref(0,1) += knoten(i,0)*refhelp(i,0);
      }
    }

// Z für Daten berechnen

  compute_Kweights();
  datamatrix diffmatrix = weighteddiffmat(2,weight);
  diffmatrix = diffmatrix.transposed()*diffmatrix.transposed().sscp().inverse();

  if(refcheck)
    Z_ref = datamatrix(1,dimZ);

  unsigned Zcols = Z.cols();
  for(j=0;j<dimZ;j++)
    {
    multBS_index(spline,diffmatrix.getCol(j));

    workdata = spline.getV();
    workZ = Z.getV()+Zpos+j;

    if(refcheck)
      {
      Z_ref(0,j) = (refhelp.transposed()*diffmatrix.getCol(j))(0,0);
      }

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

  if(!varcoeff)
    {
// X für tvalues berechen

    t_X = datamatrix(tvalues.rows(),2,1);

    spline = datamatrix(tvalues.rows(),1,0);
    multDG(spline,knoten);

    workdata = spline.getV();
    for (i=0;i<spline.rows();i++,workdata++)
      {
      t_X(i,1) = *workdata;
      }

// Z für tvalues berechnen

    t_Z = datamatrix(tvalues.rows(),dimZ,0);
    for(j=0;j<dimZ;j++)
      {
      multDG(spline,diffmatrix.getCol(j));
      workdata = spline.getV();

      for (i=0;i<spline.rows();i++,workdata++)
        {
        t_Z(i,j) = *workdata;
        }
      }
    }

  if(gridsize>0)
    {
    X_grid=t_X.getCol(1);
    Z_grid=t_Z;
    }


  }

void baseline_reml::multDG(datamatrix & res, const datamatrix & b)
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
  for(i=0;i<res.rows();i++,workres++)
    for(j=0;j<degree+1;j++,workDG++)
      *workres += *workDG * b(DGfirst[i] + j,0);

  }

void baseline_reml::initialize_baseline(unsigned j, datamatrix & tx,
               datamatrix & tz, vector<unsigned> & ts, vector<unsigned> & te,
               vector<unsigned> & tt,  datamatrix & iv,
               statmatrix<double> & steps, statmatrix<int> & ind)
  {
  if(j==0)
    {
    tx = t_X;
    tz = t_Z;
    ts = tleft;
    te = tright;
    tt = ttrunc;
    steps = tsteps;
    ind = index;
    }
  iv.putCol(j,interact_var);

  }

void baseline_reml::outoptionsreml()
  {
  if(!varcoeff)
    optionsp->out("OPTIONS FOR BASELINE TERM:: " + title + " (log(baseline))\n",true);
  else
    optionsp->out("OPTIONS FOR PSPLINE TERM:: " + title + "\n",true);
  optionsp->out("\n");

  optionsp->out("  Prior: second order random walk\n");
  if(!varcoeff)
    {
    optionsp->out("  Number of knots: " + ST::inttostring(nrknots) + "\n" );
    optionsp->out("  Degree of Splines: " + ST::inttostring(degree) + "\n" );
    optionsp->out("  Starting value for lambda: " + ST::doubletostring(startlambda,6) + "\n" );
    if(gridpos==MCMC::equidistant)
      {
      optionsp->out("  Grid choice for numerical integration: equidistant");
      optionsp->out("  Number of grid points: " + ST::inttostring(tgrid) +"\n");
      }
    else if(gridpos==MCMC::quantiles)
      {
      optionsp->out("  Grid choice for numerical integration: quantiles");
      optionsp->out("  Number of quantiles: " + ST::inttostring(nrquant) +"\n");
      optionsp->out("  Number of points between quantiles: " + ST::inttostring(nrbetween) +"\n");
      }
    else
      {
      optionsp->out("  Grid choice for numerical integration: all");
      }
    }
  optionsp->out("\n");
  }

void baseline_reml::init_name(const ST::string & na)
  {

  FULLCOND::init_name(na);

  ST::string underscore = "\\_";
  ST::string helpname = na.insert_string_char('_',underscore);
  term_symbolic = "f_{" + helpname + "}(" + helpname + ")";
  priorassumptions.push_back("$" + term_symbolic + "$:");
  priorassumptions.push_back("P-spline with second order random walk penalty");

  ST::string knotstr;
  if (knpos == equidistant)
    knotstr = "equidistant";
  else if (knpos == quantiles)
    knotstr = "quantiles";
  priorassumptions.push_back("Number of knots: " + ST::inttostring(nrknots));
  priorassumptions.push_back("Knot choice: " + knotstr);
  priorassumptions.push_back("Degree of Splines: " + ST::inttostring(degree));

  if(gridpos==MCMC::equidistant)
    {
    priorassumptions.push_back("Grid choice for numerical integration: equidistant");
    priorassumptions.push_back("Number of grid points: " + ST::inttostring(tgrid) +"\n");
    }
  else if(gridpos==MCMC::quantiles)
    {
    priorassumptions.push_back("Grid choice for numerical integration: quantiles");
    priorassumptions.push_back("Number of quantiles: " + ST::inttostring(nrquant) +"\n");
    priorassumptions.push_back("Number of points between quantiles: " + ST::inttostring(nrbetween) +"\n");
    }
  else
    {
    priorassumptions.push_back("Grid choice for numerical integration: all");
    priorassumptions.push_back("Number of grid points: " + ST::inttostring(nrdiffobs) +"\n");
    }
  }

} // end: namespace MCMC
