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





#include "baseline.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//---------- CLASS: pspline_baseline (implementation of member functions) ------
//------------------------------------------------------------------------------





//------------------------------------------------------------------------------
//----------------Constructor für zeitlich variierende Effekte------------------
//------------------------------------------------------------------------------
pspline_baseline::pspline_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & time, const datamatrix & z,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres,
                    const int & gs,const unsigned & c,const datamatrix & anfang)
 : FULLCOND_pspline(o,dp,fcc,ft,ti,nrk,degr,kp,fp,pres,false,gs,c)

  {

  unsigned i;
  gauss_n=9;
  vc_dummy1 = false;
  baselinep = vector<pspline_baseline*>(0);
  Weibull = false;
// NEW FOR PARTIALLIKELIHOOD
  PartialLikelihood = false;

  lambda = l;
  sigma2 = 1.0/l;

  zi = time;
  z_vc = z;

//------------Linkstrunkierung oder zeitlich variierende Kovariablen?-----------
  if(anfang.rows()==1)
    {
    begin0 = true;
    beg_i = datamatrix(zi.rows(),1,0);
    }
  else
    {
    begin0 = false;
    beg_i = anfang;
    }

//--bei linkstrunkierten Daten oder zeitlich var. Kovariablen Beginn der Beob.zeit berücksichtigen,
//--d.h. statt zi -> zi_ges=(z_i,beg_i)
//--statt index -> ges_index
//--und statt "multBS" -> testmat verwenden------------------------------------------------------------

/*  zi_ges = datamatrix(2*zi.rows(),1,0);

  vector<datamatrix> gaussy(gauss_n);
  for(i=0;i<gauss_n;i++)
    {
    gaussy[i] = datamatrix(zi.rows()+1,1,0);
    }

  for(i=0;i<zi.rows();i++)
    {
    zi_ges(i,0)=zi(i,0);
    zi_ges(zi.rows()+i,0)=beg_i(i,0);
    }

  ges_index = statmatrix<int>(zi_ges.rows(),1);
  ges_index.indexinit();
  zi_ges.indexsort(ges_index,0,zi_ges.rows()-1,0,0);

//  testmat = MCMC::bsplinemat(zi_ges,nrk,degr,kp,true,knot);
  testmat = MCMC::bsplinemat(zi_ges,nrk,degr,kp,true);

  coeff = datamatrix(gauss_n,1,0);
  coeff(0,0) = 0.330239355001260;
  coeff(1,0) = 0.312347077040003;
  coeff(2,0) = 0.312347077040003;
  coeff(3,0) = 0.260610696402935;
  coeff(4,0) = 0.260610696402935;
  coeff(5,0) = 0.180648160694857;
  coeff(6,0) = 0.180648160694857;
  coeff(7,0) = 0.081274388361574;
  coeff(8,0) = 0.081274388361574;

  double help1,help2;
  for(i=0;i<zi.rows();i++)
    {
    help1 = (zi(i,0)-beg_i(i,0))*0.5;
    help2 = (zi(i,0)+beg_i(i,0))*0.5;
//    if(gauss_n==3)
//      {
//      gaussy[0](i,0)= help2;
//      gaussy[1](i,0)= help2 + 0.774596669241483*help1;
//      gaussy[2](i,0)= help2 - 0.774596669241483*help1;
//      }
//    if(gauss_n==5)
//      {
//      gaussy[0](i,0) = help2;
//      gaussy[1](i,0) = help2 + 0.538469310105683*help1;
//      gaussy[2](i,0) = help2 - 0.538469310105683*help1;
//      gaussy[3](i,0) = help2 + 0.906179845938664*help1;
//      gaussy[4](i,0) = help2 - 0.906179845938664*help1;
//      }
    gaussy[0](i,0) = help2;
    gaussy[1](i,0) = help2 + 0.324253423403809*help1;
    gaussy[2](i,0) = help2 - 0.324253423403809*help1;
    gaussy[3](i,0) = help2 + 0.613371432700590*help1;
    gaussy[4](i,0) = help2 - 0.613371432700590*help1;
    gaussy[5](i,0) = help2 + 0.836031107326636*help1;
    gaussy[6](i,0) = help2 - 0.836031107326636*help1;
    gaussy[7](i,0) = help2 + 0.968160239507626*help1;
    gaussy[8](i,0) = help2 - 0.968160239507626*help1;
    }

  double maxzi=0.0;
  for(i=0;i<zi.rows();i++)
    if (zi(i,0)>maxzi) maxzi=zi(i,0);

  gaussmat = vector<MCMC::bsplinemat>(gauss_n);
  for(i=0;i<gauss_n;i++)
    {
    gaussy[i](zi.rows(),0)=maxzi;
//    gaussmat[i] = MCMC::bsplinemat(gaussy[i],nrk,degr,kp,true,knot);
    gaussmat[i] = MCMC::bsplinemat(gaussy[i],nrk,degr,kp,true);
    }*/

//------------------------------------------------------------------------------------------------------*/

  oldacceptance = 0;
  oldnrtrials = 0;

  min = minb;
  max = maxb;
  mintoobig = false;
  maxtoobig = false;

  varcoeff = true;
  setbeta(nrknots-1+degree,1,0);
  betaold = datamatrix(nrpar,1,0);

  make_index(time,z);
  make_index2();
  make_Bspline(time,true);
  make_BS(z);


// xvalues und fchelp initialisieren
  ST::string pnt = fp.substr(0,fp.length()-4)+"_fchelp.raw";
  vector<int>::iterator freqwork = freqoutput.begin();
  int * workindex = index.getV();
  if(gridsize < 0)
    {
    xvalues = datamatrix(nrdiffobs,1,0);
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        xvalues(*freqwork,0) = time(*workindex,0);
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"fchelp",nrdiffobs,1,pnt);
    splinehelp = datamatrix(likep->get_nrobs(),1,0);
    }
  else
    {
    double xmin = time.min(0);
    double xmax = time.max(0);
//    xmin=0.0;
//    xmax=2.0;
    xvalues = datamatrix(gridsize,1);
    for(i=0;i<gridsize;i++)
      xvalues(i,0) = xmin + i*(xmax-xmin)/double(xvalues.rows()-1);
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"fchelp",gridsize,1,pnt);
    splinehelp = datamatrix(gridsize,1,0);
    make_DG();
    }
  fchelp.setflags(MCMC::norelchange | MCMC::nooutput);

  compute_Kweights();

  if (type == RW1)
    {
    K = Krw1(weight);
    rankK = K.get_rows()-1;
    }
  else if (type == RW2)
    {
    K = Krw2(weight);
    rankK = K.get_rows()-2;
    }

  if(minb == 0 && maxb == 0)
    {
    automatic = true;

    min = 1;
    max = rankK;
    minauto = nrpar/5;
    maxauto = nrpar/3;
    if(minauto < 1)
      minauto = 1;
    }
  else
    {
    automatic = false;

    if(max > rankK || max == 0)
      {
      maxtoobig = true;
      max = rankK;
      }
    if(min > max || min == 0)
      {
      mintoobig = true;
      min = 1;
      }
    }

  for(i=0;i<max;i++)
    {
    fc_random.push_back(datamatrix(i+1,1,0));
    randnorm.push_back(datamatrix(i+1,1,0));
    }

  make_Kab_list();

  identifiable = true;

  compute_betaweight();


//------------------Designmatrix int_D für P-Spline an Knoten-------------------
  double knot_min = 0.0;
  double knot_max = zi.max(0);
  int_knots=datamatrix (50,1,0);
  unsigned j;
  for(j=0;j<int_knots.rows();j++)
    int_knots(j,0) = knot_min + j*(knot_max-knot_min)/double(int_knots.rows()-1);
  int_D = datamatrix(int_knots.rows(),nrpar,0.0);
  datamatrix bsp;
  for(i=0;i<int_knots.rows();i++)
    {
    bsp = bspline(int_knots(i,0));
    for(j=0;j<nrpar;j++)
      {
      int_D(i,j) = bsp(j,0);
      }
    }
//-------------------------------------------------------------------------------------------
  spline_ges = datamatrix(2*likep->get_nrobs(),1,0);
  spline_ges2 = datamatrix(2*likep->get_nrobs(),1,0);
  gaussspline = datamatrix(zi.rows()+1,gauss_n,0);
  int_ti_help = datamatrix(2*likep->get_nrobs(),1,0);
  spline_zi = datamatrix(likep->get_nrobs(),1,0);
//-----------------------
  }


//------------------------------------------------------------------------------
//-----------------------Baseline-----------------------------------------------
//------------------------------------------------------------------------------
pspline_baseline::pspline_baseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d, const double & a, const double & b,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres,
                    const int & gs, const unsigned & c,const datamatrix & anfang, const bool & wb,
                    const bool & partlik)    // NEW FOR PARTIALLIKELIHOOD
  : FULLCOND_pspline(o,dp,fcc,ft,ti,nrk,degr,kp,fp,pres,false,gs,c)
  {
  unsigned i,j;
  gauss_n = 9;
  vc_dummy1 = false;
  baseline = true;
  baselinep = vector<pspline_baseline*>(0);

  lambda = l;
  sigma2 = 1.0/l;

  zi = d;

  Weibull = wb;

  if(anfang.rows()==1)
    {
    begin0 = true;
    beg_i = datamatrix(zi.rows(),1,0);
    }
  else
    {
    begin0 = false;
    beg_i = anfang;
    }

//-----------------linkstrunkiert oder zeitl. variierende Kovariablen?--------
/*  zi_ges = datamatrix(2*zi.rows(),1,0);

  vector<datamatrix> gaussy(gauss_n);
  for(i=0;i<gauss_n;i++)
    {
    gaussy[i] = datamatrix(zi.rows()+1,1,0);
    }

//  make_index(d);
//  make_index2();
//  make_Bspline(d,true);


  for(i=0;i<zi.rows();i++)
    {
    zi_ges(i,0) = zi(i,0);
    zi_ges(zi.rows()+i,0) = beg_i(i,0);
    }
  ges_index = statmatrix<int>(zi_ges.rows(),1);
  ges_index.indexinit();
  zi_ges.indexsort(ges_index,0,zi_ges.rows()-1,0,0);

//  testmat = MCMC::bsplinemat(zi_ges,nrk,degr,kp,true,knot);
  testmat = MCMC::bsplinemat(zi_ges,nrk,degr,kp,true);

  coeff = datamatrix(gauss_n,1,0);
  coeff(0,0) = 0.330239355001260;
  coeff(1,0) = 0.312347077040003;
  coeff(2,0) = 0.312347077040003;
  coeff(3,0) = 0.260610696402935;
  coeff(4,0) = 0.260610696402935;
  coeff(5,0) = 0.180648160694857;
  coeff(6,0) = 0.180648160694857;
  coeff(7,0) = 0.081274388361574;
  coeff(8,0) = 0.081274388361574;

  double help1,help2;
  for(i=0;i<zi.rows();i++)
    {
    help1 = (zi(i,0)-beg_i(i,0))*0.5;
    help2 = (zi(i,0)+beg_i(i,0))*0.5;
//    if(gauss_n==3)
//      {
//      gaussy[0](i,0)= help2;
//      gaussy[1](i,0)= help2 + 0.774596669241483*help1;
//      gaussy[2](i,0)= help2 - 0.774596669241483*help1;
//      }
//    if(gauss_n==5)
//      {
//      gaussy[0](i,0) = help2;
//      gaussy[1](i,0) = help2 + 0.538469310105683*help1;
//      gaussy[2](i,0) = help2 - 0.538469310105683*help1;
//      gaussy[3](i,0) = help2 + 0.906179845938664*help1;
//      gaussy[4](i,0) = help2 - 0.906179845938664*help1;
//      }
    gaussy[0](i,0) = help2;
    gaussy[1](i,0) = help2 + 0.324253423403809*help1;
    gaussy[2](i,0) = help2 - 0.324253423403809*help1;
    gaussy[3](i,0) = help2 + 0.613371432700590*help1;
    gaussy[4](i,0) = help2 - 0.613371432700590*help1;
    gaussy[5](i,0) = help2 + 0.836031107326636*help1;
    gaussy[6](i,0) = help2 - 0.836031107326636*help1;
    gaussy[7](i,0) = help2 + 0.968160239507626*help1;
    gaussy[8](i,0) = help2 - 0.968160239507626*help1;
    }

  double maxzi=0.0;
  for(i=0;i<zi.rows();i++)
    if (zi(i,0)>maxzi) maxzi=zi(i,0);

  gaussmat = vector<MCMC::bsplinemat>(gauss_n);
  for(i=0;i<gauss_n;i++)
    {
    gaussy[i](zi.rows(),0)=maxzi;
//    gaussmat[i] = MCMC::bsplinemat(gaussy[i],nrk,degr,kp,true,knot);
    gaussmat[i] = MCMC::bsplinemat(gaussy[i],nrk,degr,kp,true);
    }*/
//-----------------------------------------------------------------------------

  oldacceptance = 0;
  oldnrtrials = 0;

  min = minb;
  max = maxb;
  mintoobig = false;
  maxtoobig = false;

  varcoeff = false;
  setbeta(nrknots-1+degree,1,0);
  betaold = datamatrix(nrpar,1,0);

  make_index(d);
  make_index2();
  make_Bspline(d,true);


// xvalues und fchelp initialisieren
  ST::string pnt = fp.substr(0,fp.length()-4)+"_fchelp.raw";
  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();
  if(gridsize < 0)
    {
    xvalues = datamatrix(nrdiffobs,1,0);
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
        xvalues(*freqwork,0) = d(*workindex,0);
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"fchelp",nrdiffobs,1,pnt);
    splinehelp = datamatrix(likep->get_nrobs(),1,0);
    }
  else
    {
    double xmin = d.min(0);
    double xmax = d.max(0);
//    xmin = 0.0;
//    xmax = 2.0;
    xvalues = datamatrix(gridsize,1);
    for(i=0;i<gridsize;i++)
      xvalues(i,0) = xmin + i*(xmax-xmin)/double(xvalues.rows()-1);
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"fchelp",gridsize,1,pnt);
    splinehelp = datamatrix(gridsize,1,0);
    make_DG();
    }
  fchelp.setflags(MCMC::norelchange | MCMC::nooutput);

  compute_Kweights();

  if (type == RW1)
    {
    K = Krw1(weight);
    rankK = K.get_rows()-1;
    }
  else if (type == RW2)
    {
    K = Krw2(weight);
    rankK = K.get_rows()-2;
    }

  if(minb == 0 && maxb == 0)
    {
    automatic = true;

    min = 1;
    max = rankK;
    minauto = nrpar/5;
    maxauto = nrpar/3;
    if(minauto < 1)
      minauto = 1;
    }
  else
    {
    automatic = false;

    if(max > rankK || max == 0)
      {
      maxtoobig = true;
      max = rankK;
      }
    if(min > max || min == 0)
      {
      mintoobig = true;
      min = 1;
      }
    }

  for(i=0;i<max;i++)
    {
    fc_random.push_back(datamatrix(i+1,1,0));
    randnorm.push_back(datamatrix(i+1,1,0));
    }

  make_Kab_list();

  identifiable = false;

  compute_betaweight();

//------------------Designmatrix int_D für P-Spline an Knoten-------------------

  double knot_min = 0.0;
  double knot_max = zi.max(0);
  int_knots=datamatrix (50,1,0);
  for(j=0;j<int_knots.rows();j++)
    int_knots(j,0) = knot_min + j*(knot_max-knot_min)/double(int_knots.rows()-1);

  int_D = datamatrix(int_knots.rows(),nrpar,0.0);
  datamatrix bsp;
  for(i=0;i<int_knots.rows();i++)
    {
    bsp = bspline(int_knots(i,0));
    for(j=0;j<nrpar;j++)
      {
      int_D(i,j) = bsp(j,0);
      }
    }
//------------------------------------------------------------------------------

  spline_ges = datamatrix(2*likep->get_nrobs(),1,0);
  spline_ges2 = datamatrix(2*likep->get_nrobs(),1,0);
  int_ti_help = datamatrix(2*likep->get_nrobs(),1,0);
  gaussspline = datamatrix(zi.rows()+1,gauss_n,0);
  spline_zi = datamatrix(likep->get_nrobs(),1,0);
//------------------------------------------------------------------------

  double sum_logti = 0.0;
  for(i=0;i<zi.rows();i++)
    {
    sum_logti = sum_logti + log(zi(i,0))*likep->get_response(i,0);
    }

// For the weibull-model
  if(a>0) weibullprior_alpha = a;
  else weibullprior_alpha = 0.001;
  weibullproposal_a1 = weibullprior_alpha + (likep->get_response()).sum(0);
  weibullproposal_a2 = 1.0/weibullprior_alpha - sum_logti;
  create_lgamma();
  b_prop = b;
  acceptance_between = 0.0;


// NEW FOR PARTIALLIKELIHOOD
//---------------------------
  PartialLikelihood = partlik;

  if(partlik)
    {
    // compute the Riskset
    PartialLikelihood_Riskset = datamatrix(likep->get_nrobs(), likep->get_nrobs(), 0);
    for(unsigned int i=0;i<likep->get_nrobs();i++)
      {
      for(unsigned int j=0;j<likep->get_nrobs();j++)
        {
         if(zi(i,0)<=zi(j,0))
           {
           PartialLikelihood_Riskset(i,j) = 1;     //upper triangle matrix if the timepoints zi are orderd
           }
        }
      }

    // matrix for the baselinecomponents
    firstevent = 0;
    lastevent = likep->get_nrobs();
    breslowdeltatime = datamatrix(likep->get_nrobs(),1,0);
    breslowbaseline = datamatrix(likep->get_nrobs(),1,0);
    breslowcumbaseline = datamatrix(likep->get_nrobs(),1,0);

    ST::string path = fp.substr(0,fp.length()-4)+"_breslowcumbaseline.raw";
    fc_breslowcumbaseline = FULLCOND(o,datamatrix(likep->get_nrobs(),1),title+"_breslow",likep->get_nrobs(),1,path);
    fc_breslowcumbaseline.setflags(MCMC::norelchange | MCMC::nooutput);
    }
  }


pspline_baseline::pspline_baseline(const pspline_baseline & fc)
  :FULLCOND_pspline(FULLCOND_pspline(fc))
  {

  begin0 = fc.begin0;
  int_knots = fc.int_knots;
  int_D = fc.int_D;
  testmat = fc.testmat;
  gaussmat = fc.gaussmat;
  coeff = fc.coeff;
  gauss_n = fc.gauss_n;
  zi = fc.zi;
  vc_dummy1 = fc.vc_dummy1;
  beg_i = fc.beg_i;
  zi_ges = fc.zi_ges;
  z_vc = fc.z_vc;
  spline_ges = fc.spline_ges;
  spline_ges2 = fc.spline_ges2;
  spline_zi = fc.spline_zi;
  gaussspline=fc.gaussspline;
  ges_index = fc.ges_index;
  int_ti_help = fc.int_ti_help;
  baselinep = fc.baselinep;
  weibullprior_alpha = fc.weibullprior_alpha;
  weibullproposal_a1 = fc.weibullproposal_a1;
  weibullproposal_a2 = fc.weibullproposal_a2;
  lgamma = fc.lgamma;
  Weibull = fc.Weibull;
  b_prop = fc.b_prop;
  acceptance_between = fc.acceptance_between;
  
  // NEW FOR PARTIALLIKELIHOOD
  PartialLikelihood = fc.PartialLikelihood;
  PartialLikelihood_Riskset = fc.PartialLikelihood_Riskset;
  firstevent = fc.firstevent;
  lastevent =fc.lastevent;
  breslowdeltatime = fc.breslowdeltatime;
  breslowbaseline = fc.breslowbaseline; 
  breslowcumbaseline = fc.breslowcumbaseline;
  fc_breslowcumbaseline = fc.fc_breslowcumbaseline;
  }


const pspline_baseline & pspline_baseline::operator=(const pspline_baseline & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_pspline::operator=(FULLCOND_pspline(fc));

  begin0 = fc.begin0;
  int_knots = fc.int_knots;
  int_D = fc.int_D;
  testmat = fc.testmat;
  gaussmat = fc.gaussmat;
  coeff = fc.coeff;
  gauss_n = fc.gauss_n;
  zi=fc.zi;
  beg_i = fc.beg_i;
  zi_ges = fc.zi_ges;
  z_vc = fc.z_vc;
  vc_dummy1 = fc.vc_dummy1;
  spline_ges = fc.spline_ges;
  spline_ges2 = fc.spline_ges2;
  spline_zi = fc.spline_zi;
  gaussspline = fc.gaussspline;
  ges_index = fc.ges_index;
  int_ti_help = fc.int_ti_help;
  baselinep = fc.baselinep;
  weibullprior_alpha = fc.weibullprior_alpha;
  weibullproposal_a1 = fc.weibullproposal_a1;
  weibullproposal_a2 = fc.weibullproposal_a2;
  lgamma = fc.lgamma;
  Weibull = fc.Weibull;
  b_prop = fc.b_prop;
  acceptance_between = fc.acceptance_between;

  // NEW FOR PARTIALLIKELIHOOD
  PartialLikelihood = fc.PartialLikelihood;
  PartialLikelihood_Riskset = fc.PartialLikelihood_Riskset; 
  firstevent = fc.firstevent;
  lastevent = fc.lastevent;
  breslowdeltatime =  fc.breslowdeltatime;
  breslowbaseline = fc.breslowbaseline;   
  breslowcumbaseline = fc.breslowcumbaseline;
  fc_breslowcumbaseline = fc.fc_breslowcumbaseline;

  return *this;
  }


void pspline_baseline::outoptions(void)
  {

  if(Weibull)
    {
    optionsp->out("  OPTIONS FOR Weibull-BASELINE: " + title + " (log(baseline))\n",true);
    optionsp->out("\n");
    }
  
  // NEW FOR PARTIALLIKELIHOOD
  if(PartialLikelihood)
    {
    optionsp->out("  Partial Likelihood is used for estimation\n",true);
    optionsp->out("\n");
    } 

  if(!Weibull && !PartialLikelihood)
    {
    if(varcoeff)
      optionsp->out("  OPTIONS FOR P-SPLINE TERM: " + title + "\n",true);
    else
      optionsp->out("  OPTIONS FOR P-SPLINE TERM: " + title + " (log(baseline))\n",true);

    if(maxtoobig || mintoobig)
      optionsp->out("\n");

    if(maxtoobig)
      optionsp->out("NOTE:  Maximum blocksize is missing or too big, "
                    + ST::inttostring(max) + " has been used\n");
    if(mintoobig)
      optionsp->out("NOTE:  Minimum blocksize is missing or too big, "
                      + ST::inttostring(min) + " has been used\n");

    spline_basis::outoptions();

    if(automatic)
      {
      optionsp->out("  Initial minimum blocksize for automatic tuning: " + ST::inttostring(minauto) + "\n");
      optionsp->out("  Initial maximum blocksize for automatic tuning: " + ST::inttostring(maxauto) + "\n");
      }
    else
      {
      optionsp->out("  Minimum blocksize: " + ST::inttostring(min) + "\n");
      optionsp->out("  Maximum blocksize: " + ST::inttostring(max) + "\n");
      }

    optionsp->out("\n");

    }
  }


/*void pspline_baseline::outresults(void)
  {

  spline_basis::outresults();

  unsigned i(),j;

  ST::string pnt = pathresult.substr(0,pathresult.length()-15)+"baseline.res";
  ofstream outres(pnt.strtochar());

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  outres << "intnr" << "   ";
  if(varcoeff)
    outres << datanames[1] << "   ";
  else
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

  datamatrix sample(optionsp->get_samplesize(),1);

  double* workxvalues = xvalues.getV();
  double* workmean = fchelp.get_betameanp();
  double* wqu1l = fchelp.get_beta_lower1_p();
  double* wqu2l = fchelp.get_beta_lower2_p();
  double* wqu50 = fchelp.get_betaqu50p();
  double* wqu1u = fchelp.get_beta_upper1_p();
  double* wqu2u = fchelp.get_beta_upper2_p();

  for(i=0;i<xvalues.rows();i++,workxvalues++,workmean++,wqu1l++,wqu2l++,wqu50++,wqu1u++,wqu2u++)
    {
    fchelp.readsample(sample,i);

    for(j=0;j<sample.rows();j++)
      sample(j,0) = exp(sample(j,0));

    outres << (i+1) << "   ";
    outres << *workxvalues << "   ";
    outres << exp(*workmean) << "   ";

    outres << sample.quantile(lower1,0) << "   ";
    outres << sample.quantile(lower2,0) << "   ";
    outres << sample.quantile(50,0) << "   ";
    outres << sample.quantile(upper1,0) << "   ";
    outres << sample.quantile(upper2,0) << "   ";
// wegen plotnonp!
    outres << 1 << "   ";
    outres << 1 << "   ";

    outres << endl;

    }

  }     */





void pspline_baseline::update(void)
  {

//  ofstream oflinpred("d:\\temp\\linpred.txt");
//  for(unsigned i=0;i<4007;i++)
//    oflinpred<<likep->get_linearpred(i,0)<<endl;
//  oflinpred.close();

if(Weibull)
{
// "...f_time_logbaseline_sample.raw" enthält in der 1. Spalte das Sample der
// alpha-Werte aus der baseline=alpha*t^(alpha-1)
// als priori für alpha wird eine Gamma(a,a)-Verteilung gewählt
// als proposal für alpha wird eine Gamma(alpha_current*b,b)-Verteilung gewählt
  unsigned i;
  if(optionsp->get_nriter()==1)
    {
    beta(0,0) = 1.0;   // startwert alpha=1
    for(i=0;i<zi.rows();i++)
      {
      spline(i,0)=log(beta(0,0))+(beta(0,0)-1.0)*log(zi(i,0));
      }
    likep->add_linearpred_m(spline,column,true);
    }

// ---------------- Weibull baseline Effekt --------------------------------//
//  unsigned i;
  double rho;
  double loglold;
  double priorold;
  double logold;
  double loglprop;
  double priorprop;
  double logprop;

  betaold(0,0) = beta(0,0);

//---------Integral berechnen----------

  rho = betaold(0,0);

  for(i=0;i<zi.rows();i++)
  {
  spline(i,0)=log(rho)+(rho-1.0)*log(zi(i,0));
  }

  compute_int_ti_weibull(rho);

//-------------------------------------

  loglold = likep->loglikelihood(true);
//  priorold = (weibullprior_alpha-1.0)*log(rho) - rho/weibullprior_alpha; //gamma(alpha,1/alpha)-Priori
  priorold = (weibullprior_alpha-1.0)*log(rho) - rho*weibullprior_alpha; //G(alpha,alpha)-priori
  logold = loglold + priorold;



  likep->substr_linearpred_m(spline,column,true);          //alten Spline abziehen


//  beta(0,0) = rand_gamma(weibullproposal_a1,weibullproposal_a2);
  beta(0,0) = rand_gamma(betaold(0,0)*b_prop,b_prop);         //=1.0 for linear model

//  beta(0,0) = betaold(0,0) + sqrt(0.01)*rand_normal();      //random-walk proposal
//  while(beta(0,0)<0.0)
//   beta(0,0) = betaold(0,0) + sqrt(0.01)*rand_normal();


//  logold = logold + (weibullproposal_a1-1.0)*log(beta(0,0)) - beta(0,0)*weibullproposal_a2; //gamma(a1,a2)-Proposal
  logold = logold + ((betaold(0,0)*b_prop)-1.0)*log(beta(0,0)) - beta(0,0)*b_prop + betaold(0,0)*b_prop*log(b_prop) - lgammafunc(betaold(0,0)*b_prop); //gamma(beta+b,b)-Proposal

 //---------Integral berechnen----------
//  double test = betaold(0,0);
  rho = beta(0,0);

  for(i=0;i<zi.rows();i++)
  {
  spline(i,0)=log(rho)+(rho-1.0)*log(zi(i,0));
  }
  compute_int_ti_weibull(rho);

  likep->add_linearpred_m(spline,column,true);    //neuen Spline addieren

 //-------------------------------------

  loglprop = likep->loglikelihood(true);
//  priorprop = (weibullprior_alpha-1.0)*log(rho) - rho/weibullprior_alpha; //gamma(alpha,1/alpha)-Priori
  priorprop = (weibullprior_alpha-1.0)*log(rho) - rho*weibullprior_alpha; //gamma(alpha,alpha)-Priori
  logprop = loglprop + priorprop;
//  logprop = logprop + (weibullproposal_a1-1.0)*log(betaold(0,0)) - betaold(0,0)*weibullproposal_a2; //gamma(a1,a2)-Proposal
  logprop = logprop + ((beta(0,0)*b_prop)-1.0)*log(betaold(0,0)) - betaold(0,0)*b_prop + beta(0,0)*b_prop*log(b_prop) - lgammafunc(beta(0,0)*b_prop);//gamma(beta*b,b)-prior

  double u = log(uniform());

  if (u <= (logprop - logold))
    {
    acceptance++;
    }
  else
    {
    beta(0,0) = betaold(0,0);
    for(i=0;i<zi.rows();i++)
      {
      spline(i,0)=log(betaold(0,0)/rho)+log(zi(i,0))*(betaold(0,0)-rho);
      }

    likep->add_linearpred_m(spline,column,true);

 //---------Integral berechnen----------

      rho = beta(0,0);

      for(i=0;i<zi.rows();i++)
        {
        spline(i,0)=log(rho)+(rho-1.0)*log(zi(i,0));
        }
      compute_int_ti_weibull(rho);

 //-------------------------------------

    }

  if( (optionsp->get_nriter() < optionsp->get_burnin()) &&
      (optionsp->get_nriter() % (200) == 0) )
    {
//    if(optionsp->get_nriter()== 200)
//      acceptance_between = 0.0;
//    else
//      acceptance_between = acceptance - acceptance_between;
    if((acceptance-acceptance_between)/2.0 < 20.0)
      b_prop = b_prop*1.5;
    if((acceptance-acceptance_between)/2.0 >= 20.0 && (acceptance-acceptance_between)/2.0 < 30.0)
      b_prop = b_prop*1.25;
    if((acceptance-acceptance_between)/2.0 >= 40.0 && (acceptance-acceptance_between)/2.0 < 50.0)
      b_prop = b_prop/1.25;
    if((acceptance-acceptance_between)/2.0 >= 50.0)
      b_prop = b_prop/1.5;
    acceptance_between = acceptance;
    }

  if((optionsp->get_nriter() == optionsp->get_burnin()))
    {
    optionsp->out("\n");
    optionsp->out("b-parameter is set to " + ST::inttostring(b_prop));
    optionsp->out("\n");
    }

// Übergabe der Baselinehazardfunktion alpha*t^(alpha-1) an das fullcondobjekt fchelp.
// Die Ausgabe erfolgt in die Datei "...f_time_logbaseline.res" zu den
// entsprechenden Zeitpunkten
  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      (optionsp->get_nriter() % (optionsp->get_step()) == 0) )
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
    
}//-------------------------Ende Weibull-baseline-----------------------------------/

//##############################################################################
// NEW FOR PARTIALLIKELIHOOD
//##############################################################################
if(PartialLikelihood)
{
  // Hilfsvariablen
  unsigned i, j;
  int helpindex = -1;
  double workintercept;
  double worklinpred;
  double workresponse;
  double riskset_linpred;
  double sumbaseline = 0.0;
  double helpcumbaseline;
  double * workbcb = fc_breslowcumbaseline.getbetapointer();
  
    
  // Startwerte fuer 1. Iteration + Initialisieren der baseline-Matrizen
  if(optionsp->get_nriter()==1)
    {
    beta(0,0) = 1.0;
    for(i=0;i<zi.rows();i++,workbcb++)
      {
      spline(i,0)=0.0;      
      breslowbaseline(i,0) = 1.0;
      breslowcumbaseline(i,0) = 1.0;
      *workbcb = 0.0;
      }
    likep->add_linearpred_m(spline,column,true);
    fc_breslowcumbaseline.update();
    workbcb = fc_breslowcumbaseline.getbetapointer();
    
    // Get index of first event
    i = 0;
    while(helpindex<0)
      {
       breslowdeltatime(i,0) = zi(i,0);
       workresponse = likep->get_response(i,0);
       if(workresponse==1)
         {
         helpindex = i;
         }
       i = i+1;  
      }
 
    firstevent = helpindex;    

    // Get index of last event
    for(i =(helpindex+1); i<zi.rows();i++)
      {
      breslowdeltatime(i,0) = zi(i,0) - zi(helpindex,0);
      workresponse = likep->get_response(i,0);  
      if(workresponse==1)
        {
        helpindex = i;
        lastevent = helpindex;
        }
      }  

//TEMP:BEGIN--------------------------------------------------------------------
//ofstream output_dt("c:/bayesx/test/test_deltatime.txt", ios::out|ios::app);
//output_dt << firstevent <<"\n";
//output_dt << lastevent <<"\n";
//output_dt << helpindex <<"\n";
//for(i =0; i<zi.rows();i++)
//{
//output_dt << i << " " << likep->get_response(i,0) << " " << zi(i,0) << " " << breslowdeltatime(i,0) << "\n";
//}
//TEMP:END----------------------------------------------------------------------


    }

  workintercept = fcconst->getbeta(0,0);

//TEMP:BEGIN--------------------------------------------------------------------
//ofstream output_1("c:/bayesx/test/test_cbh.txt", ios::out|ios::app);

//TEMP:END----------------------------------------------------------------------

  // compute current values of the log-baselinehazard
  for(i=0;i<zi.rows();i++)
  {
  spline(i,0) = log(breslowbaseline(i,0));
  }

  // subtraction of the current log-baelinnehazard from the predictor
  likep->substr_linearpred_m(spline,column,true);

  // compute new value of the Breslows log-baselinehazard
  for(i=0;i<zi.rows();i++,workbcb++)
    {    
    workresponse = likep->get_response(i,0);                   // value of the response delta_i {0,1}
    riskset_linpred = 0.0;
  
    for(j=0;j<zi.rows();j++)
      {
      worklinpred = likep->get_linearpred(j,0) - workintercept;
      riskset_linpred = riskset_linpred + PartialLikelihood_Riskset(i,j) * exp(worklinpred);
      }
    
    helpcumbaseline = workresponse/riskset_linpred;  
    sumbaseline = sumbaseline + helpcumbaseline;
    breslowcumbaseline(i,0) = sumbaseline;  
    *workbcb = sumbaseline;
  
    if(workresponse==1.0)
      {
      breslowbaseline(i,0) = workresponse/(breslowdeltatime(i,0)*riskset_linpred);
      }

    //breslowbaseline(i,0) = 1;//fuer weibull mit alpha=1
    //breslowcumbaseline(i,0) = zi(i,0);//fuer weibull mit alpha=1  
    }
    
    if(firstevent>0)
    {
      for(i=0;i<firstevent;i++)
       {
       breslowbaseline(i,0) = breslowbaseline(firstevent,0);
       breslowcumbaseline(i,0) = 0.0;//breslowcumbaseline(firstevent,0);
       }
    }
    
    for(i=(lastevent-1);i>firstevent;i--)
    {
      workresponse = likep->get_response(i,0);
      if(workresponse==0.0)
        {
        breslowbaseline(i,0) = breslowbaseline(i+1,0);
//        output_1 << i << "\n";
        }
    }
    if(lastevent<(zi.rows()-1))
    {
      for(i=(lastevent+1);i<zi.rows();i++)
       {
       breslowbaseline(i,0) = breslowbaseline(lastevent,0);
       }
    }

    for(i=0;i<zi.rows();i++)
    {
    spline(i,0) = log(breslowbaseline(i,0));
    }
    
    fc_breslowcumbaseline.update();
    
  // compute the new ratio log-cumbaseline(zi)/logbaseline(zi) to enter the likelihood
  compute_int_ti_partiallikelihood(breslowcumbaseline, breslowbaseline);
//  compute_int_ti_weibull(beta(0,0));

  // ddition of the new log-baelinnehazard from the predictor
  likep->add_linearpred_m(spline,column,true);    //neuen Spline addieren

 //-------------------------------------

  acceptance++;

// übergabe der baselinehazardfunktion an das fullcondobjekt fchelp
// die Ausgabe erfolgt in die Datei "...f_time_logbaseline.res" zu den
// entsprechenden zeitpunkten
  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      (optionsp->get_nriter() % (optionsp->get_step()) == 0) )
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

  if(optionsp->get_nriter()==optionsp->get_iterations())
    {
      fc_breslowcumbaseline.outresults();
    
      ST::string l1 = ST::doubletostring(lower1,4);
      ST::string l2 = ST::doubletostring(lower2,4);
      ST::string u1 = ST::doubletostring(upper1,4);
      ST::string u2 = ST::doubletostring(upper2,4);
    
      ST::string nl1 = ST::doubletostring(lower1,4);
      ST::string nl2 = ST::doubletostring(lower2,4);
      ST::string nu1 = ST::doubletostring(upper1,4);
      ST::string nu2 = ST::doubletostring(upper2,4);
      nl1 = nl1.replaceallsigns('.','p');
      nl2 = nl2.replaceallsigns('.','p');
      nu1 = nu1.replaceallsigns('.','p');
      nu2 = nu2.replaceallsigns('.','p');
    
    
      // Datei-Ausgabe Ergebnisse
      ST::string breslowcumbaseline_pathresults = pathresult.substr(0,pathresult.length()-15) + "breslowcumbaseline.res";
      ofstream ou(breslowcumbaseline_pathresults.strtochar());

    //  ou << "varname  pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
    //  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
      ou << "time  delta  pmean   pqu"  << nl1 << "   pqu" << nu2 << endl;  
      for(unsigned i=0; i<zi.rows(); i++)
        {
        ou << zi(i,0) << "  " << likep->get_response(i,0) << "  " << fc_breslowcumbaseline.get_betamean(i,0) << "  " 
           << fc_breslowcumbaseline.get_beta_lower1(i,0) << "  "
           << fc_breslowcumbaseline.get_beta_upper1(i,0) << "  "
           << endl;
        }
      
      optionsp->out("  Results for the breslowcumbaseline parameter are also stored in file\n");
      optionsp->out("  " + breslowcumbaseline_pathresults + "\n");
    
      optionsp->out("\n");
    }
    
}


/*void pspline_baseline::outresults_breslowcumbaseline(void)
  {

  fc_breslowcumbaseline.outresults();

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);

  ST::string nl1 = ST::doubletostring(lower1,4);
  ST::string nl2 = ST::doubletostring(lower2,4);
  ST::string nu1 = ST::doubletostring(upper1,4);
  ST::string nu2 = ST::doubletostring(upper2,4);
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');


  // Datei-Ausgabe Ergebnisse
  ST::string breslowcumbaseline_pathresults = fp.substr(0,fp.length()-7) + "breslowcumbaseline.res";
  ofstream ou(breslowcumbaseline_pathresults.strtochar());

//  ou << "varname  pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
//  nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
  ou << "varname" << endl;  
  for(unsigned i=0; i<beta.rows(); i++)
    {
    ou << fc_breslowcumbaseline.get_betamean(i,0) << "  "<< endl;
    }
  
  optionsp->out("  Results for the breslowcumbaseline parameter are also stored in file\n");
  optionsp->out("  " + breslowcumbaseline_pathresults + "\n");

  optionsp->out("\n");
  } 
*/
  
// END FOR PARTIALLIKELIHOOD
//##############################################################################


/*/ ---------------- linearer baseline Effekt --------------------------------//
  unsigned i;
  double m;
  double var;
  double logold;
  double logprop;

  betaold(0,0) = beta(0,0);

//---------Integral berechnen----------

  m = betaold(0,0);

  for(i=0;i<zi.rows();i++)
  {
  spline(i,0)=m*zi(i,0);
  }

  compute_int_ti_linear(m);

//-------------------------------------

  logold = likep->loglikelihood(true);

  likep->substr_linearpred_m(spline,column,true);

  var = 0.0005;
  beta(0,0) = betaold(0,0) + sqrt(var)*rand_normal();

 //---------Integral berechnen----------

  m = beta(0,0);

  for(i=0;i<zi.rows();i++)
  {
  spline(i,0)=m*zi(i,0);
  }
  compute_int_ti_linear(m);

  likep->add_linearpred_m(spline,column,true);

 //-------------------------------------

  logprop = likep->loglikelihood(true);

  double u = log(uniform());

  if (u <= (logprop-logold))
    {
    acceptance++;
    }
  else
    {
    beta(0,0) = betaold(0,0);
    for(i=0;i<zi.rows();i++)
    {
    spline(i,0)=(betaold(0,0)-m)*zi(i,0);
    }

    likep->add_linearpred_m(spline,column,true);

 //---------Integral berechnen----------

      m = beta(0,0);

      for(i=0;i<zi.rows();i++)
      {
      spline(i,0)=m*zi(i,0);
      }
      compute_int_ti_linear(m);

 //-------------------------------------

    }

  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      (optionsp->get_nriter() % (optionsp->get_step()) == 0) )
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

// ---------------- ENDE: linearer baseline Effekt ---------------------------*/


//-handelt es sich bei der Kov. mit zeitl. var. Effekt um eine dummykodierte
//-und liegt keine Linkstrunkierung bzw. zeitl. var. Kov. vor?-//
//else
if(!Weibull && !PartialLikelihood)
{

  if(optionsp->get_nriter()==1 && baselinep.size()>1 && begin0==true)
    {
    unsigned i=0;

    statmatrix<double*> z_vc_help;
    z_vc_help = statmatrix<double*>(baselinep.size()-1,1);

    double sum_z_vc = 0.0;

    for(unsigned i_vc=0;i_vc<baselinep.size()-1;i_vc++)
      z_vc_help(i_vc,0) = baselinep[i_vc+1]->get_z_vc();

    while( i<zi.rows() && (sum_z_vc==0.0||sum_z_vc==1.0))
      {
      sum_z_vc = 0.0;
      for(unsigned i_vc=0;i_vc<baselinep.size()-1;i_vc++)
        {
        if(*z_vc_help(i_vc,0)==0.0||*z_vc_help(i_vc,0)==1.0)
          sum_z_vc = sum_z_vc + *z_vc_help(i_vc,0);
        else sum_z_vc = 2.0;
          z_vc_help(i_vc,0)++;
        }
      i++;
      }
    if(i==zi.rows()) vc_dummy1=true;
    }
//------------------------------------------------------------------------------


  unsigned blocksize;

  if(automatic)
    {
     if(optionsp->get_nriter()%100==0 && optionsp->get_nriter()<optionsp->get_burnin())
      adjust_blocksize(30,70);
//      blocksize = minauto + random(maxauto-minauto+1);
      blocksize = minauto + rand() % (maxauto-minauto+1);
    }
  else
    {
//    blocksize = min + random(max-min+1);
    blocksize = min + rand() % (max-min+1);
    }

  double u;
  unsigned an = 1;
  unsigned en = blocksize;

  unsigned beg;

  double logold;
  double logprop;
  double * workbeta;
  unsigned i,j,k;



  for(j=0;j<matquant[blocksize-min];j++)
    {
    nrtrials++;

    compute_fc(beta,blocksize,an,en,sqrt(sigma2));

    logold = 0;
    logprop = 0;

    beg = firstnonzero[an-1];

    logold += likep->loglikelihood(beg,likep->get_nrobs()-1,index);

    likep->assign(false);

    add_linearpred_multBS_Block(an-1,en-1,fc_random[en-an]);

    double * workbetaold = betaold.getV()+an-1;
    workbeta = beta.getV()+an-1;          // change beta
    for(k=an-1;k<en;k++,workbeta++,workbetaold++)
      {
      *workbetaold = *workbeta;
      *workbeta = fc_random[en-an](k-an+1,0);
      }
//----------------------------------------------------------
    update_baseline();
//----------------------------------------------------------

    logprop += likep->loglikelihood(beg,likep->get_nrobs()-1,index,false);

    u = log(uniform());

    if (u <= (logprop-logold))
      {
	  acceptance++;
      likep->swap_linearpred();
      }
    else
      {
      workbetaold = betaold.getV()+an-1;
      workbeta = beta.getV()+an-1;
      for(k=an-1;k<en;k++,workbeta++,workbetaold++)
        *workbeta = *workbetaold;

//---------Integral berechnen für vorgeschlagenes beta--------
      update_baseline();
//------------------------------------/

      }

    an+=blocksize;
    if (j == matquant[blocksize-min]-2)
      en = nrpar;
    else
      en+=blocksize;

    } // end: for(j=0;j<matquant[blocksize-min];j++)  */

  if (center)
    {
    compute_intercept();

    for(i=0;i<nrpar;i++)
      beta(i,0) -= intercept;

//-------------------Spline zentrieren------------------------

    if(baselinep.size()>1)
      {
      if(vc_dummy1==true)
        {
        for(i=0;i<likep->get_nrobs();i++)
          {
          spline(i,0) -= intercept;
          spline_zi(i,0) -= intercept;
          }
        }
//      else
//        {
        for(i=0;i<2.0*likep->get_nrobs();i++)
          {
          spline_ges(i,0) -= intercept;
          spline_ges2(i,0) -= intercept;
          }
/*       for(i=0;i<likep->get_nrobs();i++)
       {
       spline1(i,0) -= intercept;
       spline2(i,0) -= intercept;
       spline3(i,0) -= intercept;
       spline4(i,0) -= intercept;
       spline5(i,0) -= intercept;
       spline6(i,0) -= intercept;
       spline7(i,0) -= intercept;
       spline8(i,0) -= intercept;
       spline9(i,0) -= intercept;
       }*/
//        }
      }
    else
      {
      if(begin0==false)
        {
        for(i=0;i<2.0*likep->get_nrobs();i++)
          {
          spline_ges(i,0) -= intercept;
          spline_ges2(i,0) -= intercept;
/*        if(i<likep->get_nrobs())
            {
            spline1(i,0) -= intercept;
            spline2(i,0) -= intercept;
            spline3(i,0) -= intercept;
            spline4(i,0) -= intercept;
            spline5(i,0) -= intercept;
            spline6(i,0) -= intercept;
            spline7(i,0) -= intercept;
            spline8(i,0) -= intercept;
            spline9(i,0) -= intercept;
            } */
          }
        }
      else
        {
        for(i=0;i<likep->get_nrobs();i++)
          spline(i,0) -= intercept;
        }
      }
    //-------------------------------------/

//    likep->add_linearpred_m(-intercept,column);
    fcconst->update_intercept(intercept);
    }

  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      (optionsp->get_nriter() % (optionsp->get_step()) == 0) )
    {

    double * fchelpbetap = fchelp.getbetapointer();

    if(gridsize < 0)
      {
      multBS(splinehelp,beta);
      vector<int>::iterator freqwork = freqoutput.begin();
      int * workindex = index.getV();
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          *fchelpbetap = splinehelp(i,0);
          fchelpbetap++;
          }
        }
      }
    else
      {
      multDG(splinehelp,beta);
      for(i=0;i<gridsize;i++,fchelpbetap++)
        *fchelpbetap = splinehelp(i,0);
      }
    }
}// Ende not Weibull
  fchelp.update();
  FULLCOND::update();

  } // end function update

//--------Für's DIC-------------------------------------------------------------
void pspline_baseline::compute_int_ti_mean(void)
  {
  unsigned i;
  if(baselinep.size()>1)
    {
    if(vc_dummy1==true)
      {
      vector <double *> splinevec;
      vector <double *> betavec;
      for(i=0;i<(baselinep.size());i++)
        splinevec.push_back(baselinep[i]->get_spline_zi_mean());
      for(i=0;i<(baselinep.size());i++)
        betavec.push_back(baselinep[i]->get_betamean());
      compute_int_ti_vc_di0(splinevec,betavec);
      for(i=1;i<baselinep.size();i++)
        {
        compute_int_ti_vc_di(i,splinevec,betavec);
        }
      }
    else
      {
      compute_int_gauss_DIC();
      }
    }
  else
    {
    if(begin0==false)
      {
      testmat.mult(spline_ges,betamean);
      testmat.mult_index(spline_ges2,betamean);
      compute_int_ti(betamean);
      }
    else
      {
      if(Weibull)
        compute_int_ti_weibull(betamean(0,0));  //für Weibull baseline
      // NEW FOR PARTIALLIKELIHOOD  
      if(PartialLikelihood)
        compute_int_ti_weibull(betamean(0,0));  //für Breslowbaseline
      //else
      if(!Weibull && !PartialLikelihood)
        {
        multBS(spline,betamean);
        compute_int_ti(betamean);
        }
      }
    }
  }


void pspline_baseline::compute_int_ti_linear(const double & b)
  {
  double * int_ti_p=likep->get_integral_ti();
  for(unsigned i=0;i<zi.rows();i++,int_ti_p++)
    {
    if(b==0)
      *int_ti_p = zi(i,0)/(exp(b*zi(i,0)));
    else
      *int_ti_p =(1/b*(exp(b*zi(i,0))-1.0))/(exp(b*zi(i,0)));
    }
  }


void pspline_baseline::compute_int_ti_weibull(const double & r)

  {
  double * int_ti_p=likep->get_integral_ti();
  for(unsigned i=0;i<zi.rows();i++,int_ti_p++)
    {
    if(r==0)
      *int_ti_p = 0.0;
    else
      {
      if(begin0==true)
        *int_ti_p = zi(i,0)/r;
      else
        *int_ti_p = (pow(zi(i,0),r)-pow(beg_i(i,0),r))/(r*pow(zi(i,0),(r-1.0)));
      }
    }
  }

// NEW FOR PARTIALLIKELIHOOD
void pspline_baseline::compute_int_ti_partiallikelihood(const datamatrix & cumbaseline, const datamatrix & baseline)
  {
  double * int_ti_p=likep->get_integral_ti();
  for(unsigned i=0;i<zi.rows();i++,int_ti_p++)
    {
    if(begin0==true)
      {
      if(cumbaseline(i,0)==0.0)
        {
        *int_ti_p = 1/baseline(i,0);
        }
      if(cumbaseline(i,0)!=0.0)
        {
        *int_ti_p = cumbaseline(i,0)/baseline(i,0);
        }
      }
      
    }
  }

//--------berechnet int_0^{t_i}/exp(logbaseline(t_i))---------------------------
//--------mit Hilfe von geordneten t_i und Knoten-------------------------------
void pspline_baseline::compute_int_ti(const datamatrix & b)
{

//------------------------left truncation----------------------------
if(begin0==false)
  {
  double * int_D_help;
  double * betap;
  double dist_knots = int_knots(1,0)-int_knots(0,0);
  unsigned i,j,k;
  k=1;
  double erg,spline_u,spline_o;
  erg = 0.0;
  double * int_ti_p = likep->get_integral_ti();
  double * int_ti_help_p = int_ti_help.getV();
  double * int_ti_help_p2 = int_ti_help.getV();

  spline_o=0.0;
  spline_u=0.0;
  int_D_help = int_D.getV();
  betap=b.getV();

  for(j=0;j<nrpar;j++,int_D_help++,betap++)
    spline_o += *betap* *int_D_help;
  spline_u=spline_o;

//--------------------------------erster Integralwert---------------------------

  while(k<int_knots.rows() && int_knots(k,0)<=zi_ges(ges_index(0,0),0) )
    {
    spline_u=spline_o;
    spline_o=0.0;
    betap = b.getV();
    for(j=0;j<nrpar;j++,int_D_help++,betap++)
      spline_o += *betap* *int_D_help;
    erg=erg+(exp(spline_u)+exp(spline_o));

    k=k+1;
    }

  erg=erg*dist_knots;
  erg=erg+(exp(spline_ges(0,0))+exp(spline_o))*(zi_ges(ges_index(0,0),0)-int_knots(k-1,0));

  int_ti_p=likep->get_integral_ti()+ges_index(0,0);
  *int_ti_p =erg*0.5/(exp(spline_ges(0,0)));

  int_ti_help_p=int_ti_help.getV()+ges_index(0,0);
     *int_ti_help_p = erg*0.5;

//------------------------------------------------------------

  for(i=1;i<zi_ges.rows();i++)
    {
    if(k==int_knots.rows())
      k=int_knots.rows()-1;
    if(k<int_knots.rows() && zi_ges(ges_index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi_ges(ges_index(i,0),0)-zi_ges(ges_index(i-1,0),0))*(exp(spline_ges(i-1,0))+exp(spline_ges(i,0)));
    else
      {
      spline_u=spline_o;
      spline_o=0.0;
      betap = b.getV();
      for(j=0;j<nrpar;j++,int_D_help++,betap++)
        spline_o += *betap* *int_D_help;
      erg=erg+(int_knots(k,0)-zi_ges(ges_index(i-1,0),0))*(exp(spline_ges(i-1,0))+exp(spline_o)) ;

      k++;

      while(k<int_knots.rows() && int_knots(k,0)<=zi_ges(ges_index(i,0),0) )
        {
        spline_u=spline_o;
        spline_o=0.0;
        betap = b.getV();
        for(j=0;j<nrpar;j++,int_D_help++,betap++)
          spline_o += *betap* *int_D_help;
        erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));

        k++;
        }
      erg=erg+(exp(spline_ges(i,0))+exp(spline_o))*(zi_ges(ges_index(i,0),0)-int_knots(k-1,0));
      }

    int_ti_p=likep->get_integral_ti()+ges_index(i,0);
    *int_ti_p =erg*0.5/(exp(spline_ges(i,0)));

    int_ti_help_p=int_ti_help.getV()+ges_index(i,0);
    *int_ti_help_p =erg*0.5;
    }
//------------------------------------------------------------------------------

  i=0;
  for(i=likep->get_nrobs();i<2*likep->get_nrobs();i++)
    {
    if(zi_ges(i,0)!=0)
      {
      int_ti_p=likep->get_integral_ti()+i-likep->get_nrobs();
      int_ti_help_p=int_ti_help.getV()+i-likep->get_nrobs();
      int_ti_help_p2=int_ti_help.getV()+i;

//      double help = (*int_ti_help_p-*int_ti_help_p2)/exp(spline_ges(ges_index(i-likep->get_nrobs(),0),0));
//      *int_ti_p = (*int_ti_help_p-*int_ti_help_p2)/exp(spline_ges(ges_index(i-likep->get_nrobs(),0),0));
      *int_ti_p = (*int_ti_help_p-*int_ti_help_p2)/exp(spline_ges2(i-likep->get_nrobs(),0));
      assert(*int_ti_p>=0.0);
      }
    }
  }//left_trunc


//--------------------Beginn=0 ---------------------
else
  {
  double * int_D_help;
  double * betap;
  double dist_knots=int_knots(1,0)-int_knots(0,0);
  unsigned i,j,k;
  k=1;


  double erg,spline_u,spline_o;
  erg = 0.0;
  double * int_ti_p=likep->get_integral_ti();
  double * int_ti_help_p=int_ti_help.getV();

  spline_o=0.0;
  spline_u=0.0;

  int_D_help =int_D.getV();
  betap=b.getV();

  for(j=0;j<nrpar;j++,int_D_help++,betap++)
    spline_o += *betap* *int_D_help;
  spline_u=spline_o;

//------------------erster Integralwert------------------------

  while(k<int_knots.rows() && int_knots(k,0)<=zi(index(0,0),0) )
    {
    spline_u=spline_o;
    spline_o=0.0;
    betap = b.getV();
    for(j=0;j<nrpar;j++,int_D_help++,betap++)
      spline_o += *betap* *int_D_help;
    erg=erg+(exp(spline_u)+exp(spline_o));
    k=k+1;
    }

  erg=erg*dist_knots;
  erg=erg+(exp(spline(0,0))+exp(spline_o))*(zi(index(0,0),0)-int_knots(k-1,0));

  int_ti_p=likep->get_integral_ti()+index(0,0);
  *int_ti_p =erg*0.5/(exp(spline(0,0)));

  int_ti_help_p=int_ti_help.getV()+index(0,0);
  *int_ti_help_p =erg*0.5;

//------------------------------------------------------------


  for(i=1;i<zi.rows();i++)
    {
    if(k==int_knots.rows())
      k=int_knots.rows()-1;
    if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi(index(i,0),0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline(i,0))) ;
    else
      {
      spline_u=spline_o;
      spline_o=0.0;
      betap = b.getV();
      for(j=0;j<nrpar;j++,int_D_help++,betap++)
        spline_o += *betap* *int_D_help;
      erg=erg+(int_knots(k,0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline_o)) ;

      k++;

      while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
        {
        spline_u=spline_o;
        spline_o=0.0;
        betap = b.getV();
        for(j=0;j<nrpar;j++,int_D_help++,betap++)
          spline_o += *betap* *int_D_help;
        erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));

        k++;
        }
      erg=erg+(exp(spline(i,0))+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));
      }

    int_ti_p=likep->get_integral_ti()+index(i,0);
    *int_ti_p =erg*0.5/(exp(spline(i,0)));

    int_ti_help_p=int_ti_help.getV()+index(i,0);
    *int_ti_help_p =erg*0.5;
    }
  }//else, i.e. not begin0==false

} //compute_ti




void pspline_baseline::compute_int_ti(unsigned beg)
{
double * int_D_help;
double * betap;
double dist_knots=int_knots(1,0)-int_knots(0,0);
unsigned i,j,k;
k=1;
double erg,spline_u,spline_o;
erg = 0.0;
double * int_ti_p=likep->get_integral_ti();
double * int_ti_help_p=int_ti_help.getV();
spline_o=0.0;
spline_u=0.0;
int_D_help =int_D.getV();
betap=beta.getV();
if(beg==0)
  {
  for(j=0;j<nrpar;j++,int_D_help++,betap++)
    spline_o += *betap* *int_D_help;
  spline_u=spline_o;

//------------------erster Integralwert------------------------

  while(k<int_knots.rows() && int_knots(k,0)<=zi(index(0,0),0) )
    {
    spline_u=spline_o;
    spline_o=0.0;

    betap = beta.getV();
    for(j=0;j<nrpar;j++,int_D_help++,betap++)
      spline_o += *betap* *int_D_help;
    erg=erg+(exp(spline_u)+exp(spline_o));

    k++;
    }

  erg=erg*dist_knots;
  erg=erg+(exp(spline(0,0))+exp(spline_o))*(zi(index(0,0),0)-int_knots(k-1,0));

  int_ti_p=likep->get_integral_ti()+index(0,0);
  *int_ti_p =erg*0.5/(exp(spline(0,0)));

  int_ti_help_p=int_ti_help.getV()+index(0,0);
  *int_ti_help_p =erg*0.5;

//------------------------------------------------------------
  for(i=1;i<zi.rows();i++)
    {
    if(k==int_knots.rows())
      k=int_knots.rows()-1;
    if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi(index(i,0),0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline(i,0))) ;
    else
      {
      spline_u=spline_o;
      spline_o=0.0;

      betap = beta.getV();
      for(j=0;j<nrpar;j++,int_D_help++,betap++)
        spline_o += *betap* *int_D_help;
      erg=erg+(int_knots(k,0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline_o)) ;

      k++;

      while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
        {
        spline_u=spline_o;
        spline_o=0.0;

        betap = beta.getV();
        for(j=0;j<nrpar;j++,int_D_help++,betap++)
          spline_o += *betap* *int_D_help;
        erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));

        k++;
        }

      erg=erg+(exp(spline(i,0))+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));
      }

    int_ti_p=likep->get_integral_ti()+index(i,0);
    *int_ti_p =erg*0.5/(exp(spline(i,0)));

    int_ti_help_p=int_ti_help.getV()+index(i,0);
    *int_ti_help_p =erg*0.5;
    }
  }//----------------if beg==0-------------------------

    //-----------------if beg!=0-----------------------------
else
  {
//---------------------erg=Integral bei "beg-1"------------------
  int_ti_help_p=int_ti_help.getV()+index(beg-1,0);
  erg=*int_ti_help_p*2.0;
//ersten Knoten finden, der größer ist als "beg-1"
  while(k<int_knots.rows() && int_knots(k,0)<=zi(index(beg-1,0),0) )
    {
    for(j=0;j<nrpar;j++)
      int_D_help++;
    k++;
    }
//--------Wert des Splines an diesem Knoten ausrechnen----------
  for(j=0;j<nrpar;j++,int_D_help++,betap++)
    spline_o += *betap* *int_D_help;
  spline_u=spline_o;
//-----------------------------------------------------------------------

  for(i=beg;i<zi.rows();i++)
    {
    if(k==int_knots.rows())
      k=int_knots.rows()-1;
    if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi(index(i,0),0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline(i,0))) ;
    else
      {
      spline_u=spline_o;
      spline_o=0.0;

      betap = beta.getV();
      for(j=0;j<nrpar;j++,int_D_help++,betap++)
        spline_o += *betap* *int_D_help;
      erg=erg+(int_knots(k,0)-zi(index(i-1,0),0))*(exp(spline(i-1,0))+exp(spline_o)) ;

      k++;

      while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
        {
        spline_u=spline_o;
        spline_o=0.0;

        betap = beta.getV();
        for(j=0;j<nrpar;j++,int_D_help++,betap++)
          spline_o += *betap* *int_D_help;
        erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));

        k++;
        }
      erg=erg+(exp(spline(i,0))+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));
      }

    int_ti_p=likep->get_integral_ti()+index(i,0);
    *int_ti_p =erg*0.5/(exp(spline(i,0)));

    int_ti_help_p=int_ti_help.getV()+index(i,0);
    *int_ti_help_p =erg*0.5;
    }
  }// ----------else, d.h. beg!=0
} //compute_ti




void pspline_baseline::compute_int_ti_vc_di0(const vector<double *> splinevector,const vector<double *> betavector)
{
//ofstream oftest_spline ("f:\\baseline\\spline_di0.txt");
double * int_D_help;
double * betap;
double dist_knots=int_knots(1,0)-int_knots(0,0);
unsigned i,j,k,i_help,i_vc;
i=0;
i_help=0;
double erg,spline_u,spline_o,spline_ti_help,spline_ti;
double * int_ti_p=likep->get_integral_ti();
double * int_ti_help_p=int_ti_help.getV();
double * spline_help;
statmatrix<double*>z_vc_help;
z_vc_help = statmatrix<double*>(baselinep.size()-1,1);
for(i_vc=1;i_vc<baselinep.size();i_vc++)
  z_vc_help(i_vc-1,0) = baselinep[i_vc]->get_z_vc();
k=1;
erg=0.0;
spline_o=0.0;
spline_u=0.0;
spline_help = splinevector[0];

int_D_help =baselinep[0]->get_int_D();
betap = betavector[0];

for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
  spline_o += *betap* *int_D_help;
spline_u=spline_o;
//------i0 berechnen---------
double z_vc_sum =0.0;
for(i_vc=1;i_vc<baselinep.size();i_vc++)
  z_vc_sum = z_vc_sum + *(z_vc_help(i_vc-1,0) + index(i,0));
while(z_vc_sum!=0.0 && i < index.rows()-1)
  {
  i++;
  spline_help++;
  z_vc_sum =0.0;
  for(i_vc=1;i_vc<baselinep.size();i_vc++)
    {
//    double test1= *(z_vc_help(i_vc-1,0) + index(i,0));
//    double test2 = index(i,0);
    z_vc_sum = z_vc_sum + *(z_vc_help(i_vc-1,0) + index(i,0));
    }
  }

i_help=i;
//------------------erster Integralwert------------------------

while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
  {
  spline_u=spline_o;
  spline_o=0.0;
  betap = betavector[0];
  for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
    spline_o += *betap* *int_D_help;
  erg=erg+(exp(spline_u)+exp(spline_o));

  k++;
  }
erg=erg*dist_knots;
spline_ti=*spline_help;
spline_ti_help=spline_ti;
erg=erg+(exp(spline_ti)+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));

int_ti_p=likep->get_integral_ti()+index(i,0);
*int_ti_p =erg*0.5/(exp(spline_ti));

int_ti_help_p=int_ti_help.getV()+index(i,0);
*int_ti_help_p =erg*0.5;
//------------------------------------------------------------

i++;
spline_help++;
while(i<zi.rows())
  {
  z_vc_sum =0.0;
  for(i_vc=1;i_vc<baselinep.size();i_vc++)
    z_vc_sum = z_vc_sum + *(z_vc_help(i_vc-1,0) + index(i,0));
  if(z_vc_sum!=0.0)
    {
    i++;
    spline_help++;
    }
  else
    {
    if(k == int_knots.rows())
      k=int_knots.rows()-1;
    spline_ti=*spline_help;
    if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi(index(i,0),0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_ti));
    else
      {
      spline_u=spline_o;
      spline_o=0.0;
      betap = betavector[0];
      for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
        spline_o += *betap* *int_D_help;
      erg=erg+(int_knots(k,0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_o)) ;

      k++;

      while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
        {
        spline_u=spline_o;
        spline_o=0.0;
        betap = betavector[0];
        for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
          spline_o += *betap* *int_D_help;
        erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));

        k++;
        }
      erg=erg+(exp(spline_ti)+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));
      }

    int_ti_p=likep->get_integral_ti()+index(i,0);
    *int_ti_p =erg*0.5/(exp(spline_ti));

    int_ti_help_p=int_ti_help.getV()+index(i,0);
    *int_ti_help_p =erg*0.5;

    i_help=i;
    i++;
    spline_help++;
    spline_ti_help=spline_ti;
    }//--------------- else: d.h. z_vc(index(i,0),0)==0 -----------------------
  } //while
//oftest_spline.close();
}



void pspline_baseline::compute_int_ti_vc_di(const int dummy, const vector<double *> splinevector,const vector<double *> betavector)
{
double * betap;
double dist_knots=int_knots(1,0)-int_knots(0,0);
unsigned i,j,k,i_help;
i=0;
i_help=0;
double erg,spline_u,spline_o,spline_ti_help,spline_ti;
double * int_ti_p=likep->get_integral_ti();
double * int_ti_help_p=int_ti_help.getV();
double * z_vc_help;
spline_ti=0.0;
spline_ti_help=0.0;
statmatrix<double*> int_D_help_1;
int_D_help_1=statmatrix<double*>(2,1);
statmatrix<double*> spline_zi_help;
spline_zi_help=statmatrix<double*>(2,1);
double help;
int i_vc;
k=1;
erg=0.0;
spline_o=0.0;
spline_u=0.0;
z_vc_help = baselinep[dummy]->get_z_vc();
spline_zi_help(0,0) = splinevector[0];
int_D_help_1(0,0)= baselinep[0]->get_int_D();
spline_zi_help(1,0) = splinevector[dummy];
int_D_help_1(1,0)= baselinep[dummy]->get_int_D();

betap = betavector[0];
help=0.0;
for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help_1(0,0)++,betap++)
  help+=*betap* *int_D_help_1(0,0);
spline_o +=help;

betap = betavector[dummy];
help=0.0;
for(j=0;j<baselinep[dummy]->get_nrpar();j++,int_D_help_1(1,0)++,betap++)
  help+=*betap* *int_D_help_1(1,0);
spline_o +=help;

spline_u=spline_o;

//------kleinstes ti mit z_vc==1 und splines an der Stelle suchen---------
while(*(z_vc_help+index(i,0))==0.0)
  {
  i++;
  for(i_vc=0;i_vc<2;i_vc++)
    {
    spline_zi_help(i_vc,0)++;
    }
  }
i_help=i;

//------------------erster Integralwert------------------------
while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
  {
  spline_u=spline_o;
  spline_o=0.0;

  betap = betavector[0];
  help=0.0;
  for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help_1(0,0)++,betap++)
    help+=*betap* *int_D_help_1(0,0);
  spline_o +=help;

  betap = betavector[dummy];
  help=0.0;
  for(j=0;j<baselinep[dummy]->get_nrpar();j++,int_D_help_1(1,0)++,betap++)
    help+=*betap* *int_D_help_1(1,0);
  spline_o +=help;

  erg=erg+(exp(spline_u)+exp(spline_o));

  k++;
  }
erg=erg*dist_knots;

//---------------Spline an der Stelle ti ---------------
for(i_vc=0;i_vc<2;i_vc++)
  {
  help= *spline_zi_help(i_vc,0);
  spline_ti +=help;
  }
spline_ti_help=spline_ti;
erg=erg+(exp(spline_ti)+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));

int_ti_p=likep->get_integral_ti()+index(i,0);
*int_ti_p =erg*0.5/(exp(spline_ti));

int_ti_help_p=int_ti_help.getV()+index(i,0);
*int_ti_help_p =erg*0.5;

//------------------------------------------------------------

i++;
for(i_vc=0;i_vc<2;i_vc++)
  spline_zi_help(i_vc,0)++;
while(i<zi.rows())
  {
//-----------falls z_vc==0: zu nächster Beobachtung---------------
  if(*(z_vc_help+index(i,0))==0)
    {
    i++;
    for(i_vc=0;i_vc<2;i_vc++)
      spline_zi_help(i_vc,0)++;
    }
//-------------falls z_vc==1: Integral berechenen---------------------
  else
    {
    if(k == int_knots.rows())
      k=int_knots.rows()-1;

//Spline an der Stelle ti
    spline_ti=0.0;
    for(i_vc=0;i_vc<2;i_vc++)
      {
      help= *spline_zi_help(i_vc,0);
      spline_ti +=help;
      }
    if(k<int_knots.rows() && zi(index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi(index(i,0),0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_ti)) ;
    else
      {
      spline_u=spline_o;
      spline_o = 0.0;

      betap = betavector[0];
      help = 0.0;
      for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help_1(0,0)++,betap++)
        help+=*betap* *int_D_help_1(0,0);
      spline_o +=help;

      betap = betavector[dummy];
      help=0.0;
      for(j=0;j<baselinep[dummy]->get_nrpar();j++,int_D_help_1(1,0)++,betap++)
        help+=*betap* *int_D_help_1(1,0);
      spline_o +=help;

      erg=erg+(int_knots(k,0)-zi(index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_o)) ;

      k++;

      while(k<int_knots.rows() && int_knots(k,0)<=zi(index(i,0),0) )
        {
        spline_u=spline_o;
        spline_o=0.0;

        betap = betavector[0];
        help=0.0;
        for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help_1(0,0)++,betap++)
          help+=*betap* *int_D_help_1(0,0);
        spline_o +=help;

        betap = betavector[dummy];
        help=0.0;
        for(j=0;j<baselinep[dummy]->get_nrpar();j++,int_D_help_1(1,0)++,betap++)
          help+=*betap* *int_D_help_1(1,0);
        spline_o +=help;

        erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));

        k++;
        }

      erg=erg+(exp(spline_ti)+exp(spline_o))*(zi(index(i,0),0)-int_knots(k-1,0));
      }

    int_ti_p=likep->get_integral_ti()+index(i,0);
    *int_ti_p =erg*0.5/(exp(spline_ti));

    int_ti_help_p=int_ti_help.getV()+index(i,0);
    *int_ti_help_p =erg*0.5;

    i_help=i;
    spline_ti_help=spline_ti;
    i++;
    for(i_vc=0;i_vc<2;i_vc++)
      spline_zi_help(i_vc,0)++;
    }//--------------- else: d.h. z_vc(index(i,0),0)==1 -----------------------
  } //while
}//--------------compute_int_ti_vc_di---------------------



void pspline_baseline::compute_int_gauss(void)
{

vector<double *> betavec;
unsigned i,j,k;
vector<double *> splinevec;
vector<double *> gausssplinevec;
vector<datamatrix > z_vc_help;

for(i=0;i<(baselinep.size());i++)
  splinevec.push_back(baselinep[i]->get_spline_ges());
for(i=0;i<(baselinep.size());i++)
  gausssplinevec.push_back(baselinep[i]->get_gaussspline());
for(i=0;i<(baselinep.size()-1);i++)
  z_vc_help.push_back(baselinep[i+1]->get_z_vc_np());

double * int_ti_p=likep->get_integral_ti();

double help,splinehelp1,splinehelp2;

for(i=0;i<zi.rows();i++,int_ti_p++)
  {
  help = 0.0;
  splinehelp2 = 0.0;
  for(j=0;j<gauss_n;j++)
    {
    splinehelp2=0.0;
    for(k=0;k<baselinep.size();k++)
      {
      splinehelp1 = *(gausssplinevec[k]);
      if(k>0) splinehelp1 = splinehelp1* (z_vc_help[k-1](i,0));
      splinehelp2 = splinehelp2 + splinehelp1;
      gausssplinevec[k]++;
      }
    help = help + coeff(j,0)*exp(splinehelp2);
    }
  splinehelp2 = 0.0;
  for(k=0;k<baselinep.size();splinevec[k]++,k++)
    {
    splinehelp1 = *(splinevec[k]);
    if(k>0)
      {
      splinehelp1 = splinehelp1* (z_vc_help[k-1](i,0));
      }
    splinehelp2 = splinehelp2 + splinehelp1;
    }
  *int_ti_p = ((zi(i,0)-beg_i(i,0))*0.5)*help/(exp(splinehelp2));
  }
}


void pspline_baseline::compute_int_gauss_DIC(void)
{

unsigned i,j,k;
vector<double *> splinevec;
vector<double *> gausssplinevec;
vector<double *> z_vc_help;

for(i=0;i<(baselinep.size());i++)
  splinevec.push_back(baselinep[i]->get_spline_ges_mean());
for(i=0;i<(baselinep.size());i++)
  gausssplinevec.push_back(baselinep[i]->get_gaussspline_mean());
for(i=0;i<(baselinep.size()-1);i++)
  z_vc_help.push_back(baselinep[i+1]->get_z_vc());

double * int_ti_p=likep->get_integral_ti();

double help,splinehelp1,splinehelp2;

for(i=0;i<zi.rows();i++,int_ti_p++)
  {
  help = 0.0;
  splinehelp2 = 0.0;
  for(j=0;j<gauss_n;j++)
    {
    splinehelp2 = 0.0;
    for(k=0;k<baselinep.size();gausssplinevec[k]++,k++)
      {
      splinehelp1 = *(gausssplinevec[k]);
      if(k>0) splinehelp1 = splinehelp1* *(z_vc_help[k-1]);
      splinehelp2 = splinehelp2 + splinehelp1;
      }

    help = help + coeff(j,0)*exp(splinehelp2);
    }
  splinehelp2 = 0.0;
  for(k=0;k<baselinep.size();splinevec[k]++,k++)
    {
    splinehelp1 = *(splinevec[k]);
    if(k>0)
      {
      splinehelp1 = splinehelp1* *(z_vc_help[k-1]);
      z_vc_help[k]++;
      }
    splinehelp2 = splinehelp2 + splinehelp1;
    }
  *int_ti_p = ((zi(i,0)-beg_i(i,0))*0.5)*help/(exp(splinehelp2));
  }
}


void pspline_baseline::update_baseline()
{
//---------Integral berechnen---------------------------------
unsigned i;
if(baselinep.size()>=1)
  {
  if(vc_dummy1==true)   //keine Linkstrunkierung, zeitl. var. Effekt für dummykod. Variable
    {
    vector <double *> splinevec;
    vector <double *> betavec;
    for(i=0;i<(baselinep.size());i++)
      splinevec.push_back(baselinep[i]->get_spline_zi());
    for(i=0;i<(baselinep.size());i++)
      betavec.push_back(baselinep[i]->getbetapointer());
    compute_int_ti_vc_di0(splinevec,betavec);
    for(i=1;i<baselinep.size();i++)
      {
      compute_int_ti_vc_di(i,splinevec,betavec);
      }
    }
  else    //zeitl. var. Effekt für beliebige Kovariablen, Linkstrunkierung
    {
    compute_int_gauss();
    }
  }
else    //kein zeitl. var. Effekt
  {
  if(begin0==false)     //Linkstrunkierung
    {
    testmat.mult(spline_ges,beta);
    testmat.mult_index(spline_ges2,beta);
    compute_int_ti(beta);
    }
  else    //keine Linkstrunkierung
    {
    multBS(spline,beta);
    compute_int_ti(beta);
    }
  }
}
//------------------------------------------------------------------------------


double * pspline_baseline::get_gaussspline()
  {
  datamatrix egon (gaussspline.rows(),1,0);
  for(unsigned i=0;i<gauss_n;i++)
    {
    egon = gaussspline.getCol(i);
    gaussmat[i].mult_index(egon,beta);
    gaussspline.putCol(i,egon);
    }
  return gaussspline.getV();
  }

double * pspline_baseline::get_gaussspline_mean()
  {
  datamatrix egon (gaussspline.rows(),1,0);
  for(unsigned i=0;i<gauss_n;i++)
    {
    egon = gaussspline.getCol(i);
    gaussmat[i].mult_index(egon,betamean);
    gaussspline.putCol(i,egon);
    }
  return gaussspline.getV();
  }



void pspline_baseline::create_lgamma(void)
  {
  lgamma = datamatrix(101,1);

  lgamma(0,0) = 0.000000000;
  lgamma(1,0) = -0.005690308;
  lgamma(2,0) = -0.011218489;
  lgamma(3,0) = -0.016586854;
  lgamma(4,0) = -0.021797651;
  lgamma(5,0) = -0.026853073;
  lgamma(6,0) = -0.031755254;
  lgamma(7,0) = -0.036506276;
  lgamma(8,0) = -0.041108170;
  lgamma(9,0) = -0.045562915;
  lgamma(10,0) = -0.049872441;
  lgamma(11,0) = -0.054038634;
  lgamma(12,0) = -0.058063333;
  lgamma(13,0) = -0.061948332;
  lgamma(14,0) = -0.065695387;
  lgamma(15,0) = -0.069306209;
  lgamma(16,0) = -0.072782472;
  lgamma(17,0) = -0.076125811;
  lgamma(18,0) = -0.079337824;
  lgamma(19,0) = -0.082420074;
  lgamma(20,0) = -0.085374090;
  lgamma(21,0) = -0.088201365;
  lgamma(22,0) = -0.090903362;
  lgamma(23,0) = -0.093481511;
  lgamma(24,0) = -0.095937212;
  lgamma(25,0) = -0.098271836;
  lgamma(26,0) = -0.100486725;
  lgamma(27,0) = -0.102583193;
  lgamma(28,0) = -0.104562527;
  lgamma(29,0) = -0.106425987;
  lgamma(30,0) = -0.108174810;
  lgamma(31,0) = -0.109810204;
  lgamma(32,0) = -0.111333359;
  lgamma(33,0) = -0.112745436;
  lgamma(34,0) = -0.114047576;
  lgamma(35,0) = -0.115240897;
  lgamma(36,0) = -0.116326498;
  lgamma(37,0) = -0.117305454;
  lgamma(38,0) = -0.118178821;
  lgamma(39,0) = -0.118947635;
  lgamma(40,0) = -0.119612914;
  lgamma(41,0) = -0.120175656;
  lgamma(42,0) = -0.120636841;
  lgamma(43,0) = -0.120997431;
  lgamma(44,0) = -0.121258371;
  lgamma(45,0) = -0.121420591;
  lgamma(46,0) = -0.121485001;
  lgamma(47,0) = -0.121452498;
  lgamma(48,0) = -0.121323962;
  lgamma(49,0) = -0.121100259;
  lgamma(50,0) = -0.120782238;
  lgamma(51,0) = -0.120370735;
  lgamma(52,0) = -0.119866573;
  lgamma(53,0) = -0.119270560;
  lgamma(54,0) = -0.118583490;
  lgamma(55,0) = -0.117806145;
  lgamma(56,0) = -0.116939293;
  lgamma(57,0) = -0.115983691;
  lgamma(58,0) = -0.114940083;
  lgamma(59,0) = -0.113809201;
  lgamma(60,0) = -0.112591766;
  lgamma(61,0) = -0.111288486;
  lgamma(62,0) = -0.109900061;
  lgamma(63,0) = -0.108427177;
  lgamma(64,0) = -0.106870510;
  lgamma(65,0) = -0.105230728;
  lgamma(66,0) = -0.103508486;
  lgamma(67,0) = -0.101704430;
  lgamma(68,0) = -0.099819197;
  lgamma(69,0) = -0.097853413;
  lgamma(70,0) = -0.095807697;
  lgamma(71,0) = -0.093682657;
  lgamma(72,0) = -0.091478893;
  lgamma(73,0) = -0.089196995;
  lgamma(74,0) = -0.086837546;
  lgamma(75,0) = -0.084401121;
  lgamma(76,0) = -0.081888285;
  lgamma(77,0) = -0.079299595;
  lgamma(78,0) = -0.076635603;
  lgamma(79,0) = -0.073896851;
  lgamma(80,0) = -0.071083873;
  lgamma(81,0) = -0.068197197;
  lgamma(82,0) = -0.065237343;
  lgamma(83,0) = -0.062204825;
  lgamma(84,0) = -0.059100148;
  lgamma(85,0) = -0.055923813;
  lgamma(86,0) = -0.052676312;
  lgamma(87,0) = -0.049358131;
  lgamma(88,0) = -0.045969750;
  lgamma(89,0) = -0.042511642;
  lgamma(90,0) = -0.038984276;
  lgamma(91,0) = -0.035388112;
  lgamma(92,0) = -0.031723605;
  lgamma(93,0) = -0.027991206;
  lgamma(94,0) = -0.024191358;
  lgamma(95,0) = -0.020324499;
  lgamma(96,0) = -0.016391062;
  lgamma(97,0) = -0.012391474;
  lgamma(98,0) = -0.008326158;
  lgamma(99,0) = -0.004195529;
  lgamma(100,0) = 0.000000000;

  }


double pspline_baseline::lgammafunc(const double & nu) const
    {

    if (fmod(nu,1)==0)
      return lfac(nu-1);
    else if (nu<1)
      return lgammafunc(nu+1) - log(nu);
    else if (nu>2)
      return log(nu-1) + lgammafunc(nu-1);
    else
      return lgamma(int(nu*100)-100,0);
    }


double pspline_baseline::lfac(const double & nu) const
    {
    if (nu==0 || nu==1) return 0;
    else return log(nu) + lfac(nu-1);
    }


} // end: namespace MCMC



