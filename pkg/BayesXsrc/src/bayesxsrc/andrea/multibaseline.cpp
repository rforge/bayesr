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





#include "multibaseline.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//-----------------------Multi-Baseline-----------------------------------------
//------------------------------------------------------------------------------
pspline_multibaseline::pspline_multibaseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d, const double & a, const double & b,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres,
                    const int & gs, const unsigned & c,const datamatrix & zustand, const datamatrix & anfang, const bool & gl)
  : FULLCOND_pspline(o,dp,fcc,ft,ti,nrk,degr,kp,fp,pres,false,gs,c)
  {
  unsigned i,j,k;

  baseline = true;
  baselinep = vector<pspline_multibaseline*>(0);

  lambda = l;
  sigma2 = 1.0/l;

  zi = d;
  col = c;

  global = gl;

  state_i = zustand;

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

  zi_ges = datamatrix(2*zi.rows(),1,0);

  for(i=0;i<zi.rows();i++)
    {
    zi_ges(i,0) = zi(i,0);
    zi_ges(zi.rows()+i,0) = beg_i(i,0);
    }
  ges_index = statmatrix<int>(zi_ges.rows(),1);
  ges_index.indexinit();
  zi_ges.indexsort(ges_index,0,zi_ges.rows()-1,0,0);

  if(global==false)
    {
    zi_teil = vector<datamatrix>(likep->get_responsedim());
    j=0;
    for(i=0;i<zi.rows();i++)
      if(likep->get_transition(i,col)==1) j=j+1;
    zi_teil[col] = datamatrix(2*j,1,0);

    teil_index = vector< statmatrix<int> >(likep->get_responsedim());
    teil_index[col] = statmatrix<int>(zi_teil[col].rows(),1);

    k=0;
    for(i=0;i<zi.rows();i++)
      {
      if(likep->get_transition(i,col)==1)
        {
        zi_teil[col](k,0) = zi(i,0);
        zi_teil[col](j+k,0) = beg_i(i,0);
        k=k+1;
        }
      }

    k=0;
    for(i=0;i<2*zi.rows();i++)
      {
//      unsigned tin = ges_index(i,0);
      if(likep->get_transition(ges_index(i,0),col)==1)
        {
        teil_index[col](k,0)=ges_index(i,0);
        k = k+1;
        }
      }
    }

  testmat = MCMC::bsplinemat(zi_ges,nrk,degr,kp,true);

  if(global==false)
    {
    testmat_l = vector<MCMC::bsplinemat>(likep->get_responsedim());
    testmat_l[col] = MCMC::bsplinemat(zi_teil[col],nrk,degr,kp,true);
    }

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

  if(global==false)
    {
    knot_max = zi_teil[col].max(0);
    }

  if(global == false)
    {
    int_knots_l = vector<datamatrix>(likep->get_responsedim());
    int_knots_l[col]=datamatrix (50,1,0);
    for(j=0;j<int_knots_l[col].rows();j++)
      int_knots_l[col](j,0) = knot_min + j*(knot_max-knot_min)/double(int_knots_l[col].rows()-1);
    int_D_l = vector<datamatrix>(likep->get_responsedim());
    int_D_l[col] = datamatrix(int_knots_l[col].rows(),nrpar,0.0);
    datamatrix bsp;
    for(i=0;i<int_knots_l[col].rows();i++)
      {
      bsp = bspline(int_knots_l[col](i,0));
      for(j=0;j<nrpar;j++)
        {
        int_D_l[col](i,j) = bsp(j,0);
        }
      }
    }
  else
    {
    int_knots = datamatrix (50,1,0);
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
    }
//------------------------------------------------------------------------------
     double * int_ti_out_p=likep->get_integral_ti()+col;
     for(i=0;i<zi.rows();i++)
       {
         *int_ti_out_p=0.1;
          int_ti_out_p= int_ti_out_p+likep->get_responsedim();
       }

//---------------------------------------------------------------------------/

//------------------------------------------------------------------------------
  int_ti_help = datamatrix(2*likep->get_nrobs(),1,0);

  if(global==true)
    {
    spline_ges = datamatrix(2*likep->get_nrobs(),1,0);
    spline_ges2 = datamatrix(2*likep->get_nrobs(),1,0);
    spline_zi = datamatrix(likep->get_nrobs(),1,0);
    }
  else
    {
    spline_teil = vector<datamatrix>(likep->get_responsedim());
    spline_teil2 = vector<datamatrix>(likep->get_responsedim());
    spline_teil[col] = datamatrix(zi_teil[col].rows(),1,0);
    spline_teil2[col] = datamatrix(zi_teil[col].rows(),1,0);
    }
//------------------------------------------------------------------------

  }


//------------------------------------------------------------------------------
//-----------------------Multi-Baseline-varcoeff--------------------------------
//------------------------------------------------------------------------------
pspline_multibaseline::pspline_multibaseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & time, const datamatrix & z,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres,
                    const int & gs, const unsigned & c,const datamatrix & zustand, const datamatrix & anfang, const bool & gl)
  : FULLCOND_pspline(o,dp,fcc,ft,ti,nrk,degr,kp,fp,pres,false,gs,c)
  {
  unsigned i,j;

  baselinep = vector<pspline_multibaseline*>(0);

  lambda = l;
  sigma2 = 1.0/l;

  zi = time;
  z_vc = datamatrix(2*z.rows(),1,0);
  for(i=0;i<z.rows();i++)
    {
    z_vc(i,0) = z(i,0);
    z_vc(z.rows()+i,0) = z(i,0);
    }

  col = c;

  global = true;

  state_i = zustand;

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

  zi_ges = datamatrix(2*zi.rows(),1,0);

  for(i=0;i<zi.rows();i++)
    {
    zi_ges(i,0) = zi(i,0);
    zi_ges(zi.rows()+i,0) = beg_i(i,0);
    }
  ges_index = statmatrix<int>(zi_ges.rows(),1);
  ges_index.indexinit();
  zi_ges.indexsort(ges_index,0,zi_ges.rows()-1,0,0);

  testmat = MCMC::bsplinemat(zi_ges,nrk,degr,kp,true);

//-----------------------------------------------------------------------------

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


  int_knots = datamatrix (50,1,0);
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
     double * int_ti_out_p=likep->get_integral_ti()+col;
     for(i=0;i<zi.rows();i++)
       {
         *int_ti_out_p=0.1;
          int_ti_out_p= int_ti_out_p+likep->get_responsedim();
       }

//---------------------------------------------------------------------------/

//------------------------------------------------------------------------------
  int_ti_help = datamatrix(2*likep->get_nrobs(),1,0);

  spline_ges = datamatrix(2*likep->get_nrobs(),1,0);
  spline_ges2 = datamatrix(2*likep->get_nrobs(),1,0);
  spline_zi = datamatrix(likep->get_nrobs(),1,0);


  }





pspline_multibaseline::pspline_multibaseline(const pspline_multibaseline & fc)
  :FULLCOND_pspline(FULLCOND_pspline(fc))
  {

  begin0 = fc.begin0;
  int_knots = fc.int_knots;
  int_knots_l = fc.int_knots_l;
  int_D = fc.int_D;
  int_D_l = fc.int_D_l;
  testmat = fc.testmat;
  testmat_l = fc.testmat_l;
  zi = fc.zi;
  vc_dummy1 = fc.vc_dummy1;
  beg_i = fc.beg_i;
  state_i = fc.state_i;
  zi_ges = fc.zi_ges;
  zi_teil = fc.zi_teil;
  z_vc = fc.z_vc;
  spline_ges = fc.spline_ges;
  spline_teil = fc.spline_teil;
  spline_teil2 = fc.spline_teil2;
  spline_ges2 = fc.spline_ges2;
  spline_zi = fc.spline_zi;
  ges_index = fc.ges_index;
  teil_index = fc.teil_index;
  int_ti_help = fc.int_ti_help;
  baselinep = fc.baselinep;
  col = fc.col;
  global = fc.global;
  }


const pspline_multibaseline & pspline_multibaseline::operator=(const pspline_multibaseline & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_pspline::operator=(FULLCOND_pspline(fc));

  begin0 = fc.begin0;
  int_knots = fc.int_knots;
  int_knots_l = fc.int_knots_l;
  int_D = fc.int_D;
  int_D_l = fc.int_D_l;
  testmat = fc.testmat;
  testmat_l = fc.testmat_l;
  zi=fc.zi;
  beg_i = fc.beg_i;
  state_i = fc.state_i;
  zi_ges = fc.zi_ges;
  zi_teil = fc.zi_teil;
  z_vc = fc.z_vc;
  vc_dummy1 = fc.vc_dummy1;
  spline_ges = fc.spline_ges;
  spline_teil = fc.spline_teil;
  spline_teil2 = fc.spline_teil2;
  spline_ges2 = fc.spline_ges2;
  spline_zi = fc.spline_zi;
  ges_index = fc.ges_index;
  teil_index = fc.teil_index;
  int_ti_help = fc.int_ti_help;
  baselinep = fc.baselinep;
  col = fc.col;
  global = fc.global;
  return *this;
  }


void pspline_multibaseline::outoptions(void)
  {

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







void pspline_multibaseline::update(void)
  {



  unsigned blocksize;

  if(automatic)
    {
     if(optionsp->get_nriter()%100==0 && optionsp->get_nriter()<optionsp->get_burnin())
       adjust_blocksize(30,70);
//     blocksize = minauto + random(maxauto-minauto+1);
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
    update_multibaseline();
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
      update_multibaseline();
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

//-------------------Spline zentrieren------------------------------------------


        if(global==true)
          {
          for(i=0;i<2.0*likep->get_nrobs();i++)
            {
            spline_ges(i,0) -= intercept;
            spline_ges2(i,0) -= intercept;
            }
          }
        else
          {
          for(i=0;i<zi_teil[col].rows();i++)
            {
            spline_teil[col](i,0) -= intercept;
            spline_teil2[col](i,0) -=intercept;
            }
          }


//------------------------------------------------------------------------------

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

  fchelp.update();
  FULLCOND::update();

  } // end function update

//--------Für's DIC-------------------------------------------------------------
void pspline_multibaseline::compute_int_ti_mean(void)
  {
  if(baselinep.size()>1)
    {
      unsigned i;
      vector <double *> splinevec;
      vector <double *> splinevec2;
      vector <double *> betavec;
      for(i=0;i<(baselinep.size());i++)
        {
        splinevec.push_back(baselinep[i]->get_spline_ges_mean());
        splinevec2.push_back(baselinep[i]->get_spline_ges2_mean());
        }
      for(i=0;i<(baselinep.size());i++)
        betavec.push_back(baselinep[i]->get_betamean());

      compute_int_ti_vc_di0(splinevec,splinevec2,betavec);
      for(i=1;i<baselinep.size();i++)
        {
        compute_int_ti_vc_di(i,splinevec,splinevec2,betavec);
        }
    }
  else
    {
    if(global==true)
      {
      testmat.mult(spline_ges,betamean);
      testmat.mult_index(spline_ges2,betamean);
      compute_int_ti_global(betamean);
      }
    else
      {
      testmat_l[col].mult(spline_teil[col],betamean);
      testmat_l[col].mult_index(spline_teil2[col],betamean);
      compute_int_ti_nonglobal(betamean);
      }
    }
  }





//--------berechnet int_0^{t_i}/exp(logbaseline(t_i))---------------------------
//--------mit Hilfe von geordneten t_i und Knoten-------------------------------
void pspline_multibaseline::compute_int_ti_nonglobal(const datamatrix & b)
{

//------------------------left truncation----------------------------
if(begin0==false)
  {
  double * int_D_help;
  double * betap;
  double dist_knots = int_knots_l[col](1,0)-int_knots_l[col](0,0);
  unsigned i,j,k;
  k=1;
  double erg,spline_u,spline_o;
  erg = 0.0;
  double * int_ti_p = likep->get_integral_ti()+col;
  double * int_ti_help_p = int_ti_help.getV();
  double * int_ti_help_p2 = int_ti_help.getV();

  spline_o=0.0;
  spline_u=0.0;
  int_D_help = int_D_l[col].getV();
  betap=b.getV();

  for(j=0;j<nrpar;j++,int_D_help++,betap++)
    spline_o += *betap* *int_D_help;
  spline_u=spline_o;

//--------------------------------erster Integralwert---------------------------

  while(k<int_knots_l[col].rows() && int_knots_l[col](k,0)<=zi_ges(teil_index[col](0,0),0))
    {
    spline_u=spline_o;
    spline_o=0.0;
    betap = b.getV();
    for(j=0;j<nrpar;j++,int_D_help++,betap++)
      spline_o += *betap* *int_D_help;
    erg=erg+(exp(spline_u)+exp(spline_o));

    k=k+1;
    }
//  double indextest=teil_index[col](0,0);
//  double zigestest=zi_ges(teil_index[col](0,0),0);
  erg=erg*dist_knots;
  erg=erg+(exp(spline_teil[col](0,0))+exp(spline_o))*(zi_ges(teil_index[col](0,0),0)-int_knots_l[col](k-1,0));

  int_ti_p=likep->get_integral_ti()+(teil_index[col](0,0)*likep->get_responsedim()+col);
  *int_ti_p =erg*0.5/(exp(spline_teil[col](0,0)));

  int_ti_help_p=int_ti_help.getV()+teil_index[col](0,0);
     *int_ti_help_p = erg*0.5;

//------------------------------------------------------------

  for(i=1;i<(zi_teil[col].rows());i++)
    {
    if(k==int_knots_l[col].rows())
      k=int_knots_l[col].rows()-1;
    if(k<int_knots_l[col].rows() && zi_ges(teil_index[col](i,0),0)<=int_knots_l[col](k,0))
      erg=erg+(zi_ges(teil_index[col](i,0),0)-zi_ges(teil_index[col](i-1,0),0))*(exp(spline_teil[col](i-1,0))+exp(spline_teil[col](i,0)));
    else
      {
      spline_u=spline_o;
      spline_o=0.0;
      betap = b.getV();
      for(j=0;j<nrpar;j++,int_D_help++,betap++)
        spline_o += *betap* *int_D_help;
      erg=erg+(int_knots_l[col](k,0)-zi_ges(teil_index[col](i-1,0),0))*(exp(spline_teil[col](i-1,0))+exp(spline_o)) ;

      k++;

      while(k<int_knots_l[col].rows() && int_knots_l[col](k,0)<=zi_ges(teil_index[col](i,0),0) )
        {
        spline_u=spline_o;
        spline_o=0.0;
        betap = b.getV();
        for(j=0;j<nrpar;j++,int_D_help++,betap++)
          spline_o += *betap* *int_D_help;
        erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));

        k++;
        }
//      double test= spline_o;
      erg=erg+(exp(spline_teil[col](i,0))+exp(spline_o))*(zi_ges(teil_index[col](i,0),0)-int_knots_l[col](k-1,0));
      }

    int_ti_p=likep->get_integral_ti()+teil_index[col](i,0)*likep->get_responsedim()+col;
    *int_ti_p =erg*0.5/(exp(spline_teil[col](i,0)));

    int_ti_help_p=int_ti_help.getV()+teil_index[col](i,0);
    *int_ti_help_p =erg*0.5;
    }
//------------------------------------------------------------------------------



  i=0;
  k=0;
  for(i=likep->get_nrobs();i<2*likep->get_nrobs();i++)
    {
    if(likep->get_transition(i-likep->get_nrobs(),col)==1)
      {
      if(zi_ges(i,0)!=0)
        {
        int_ti_p=likep->get_integral_ti()+(i-likep->get_nrobs())*likep->get_responsedim()+col;
        int_ti_help_p=int_ti_help.getV()+i-likep->get_nrobs();
        int_ti_help_p2=int_ti_help.getV()+i;
        *int_ti_p = (*int_ti_help_p-*int_ti_help_p2)/exp(spline_teil2[col](k,0));
        assert(*int_ti_p>=0.0);
        }
      k = k+1;
      }
    }
  }//left_trunc



} //compute_ti



void pspline_multibaseline::compute_int_ti_global(const datamatrix & b)
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
  double * int_ti_p = likep->get_integral_ti()+col;
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

  int_ti_p=likep->get_integral_ti()+(ges_index(0,0)*likep->get_responsedim()+col);
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
//      double test= spline_o;
      erg=erg+(exp(spline_ges(i,0))+exp(spline_o))*(zi_ges(ges_index(i,0),0)-int_knots(k-1,0));
      }

    int_ti_p=likep->get_integral_ti()+ges_index(i,0)*likep->get_responsedim()+col;
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
      int_ti_p=likep->get_integral_ti()+(i-likep->get_nrobs())*likep->get_responsedim()+col;
      int_ti_help_p=int_ti_help.getV()+i-likep->get_nrobs();
      int_ti_help_p2=int_ti_help.getV()+i;
      *int_ti_p = (*int_ti_help_p-*int_ti_help_p2)/exp(spline_ges2(i-likep->get_nrobs(),0));
      assert(*int_ti_p>=0.0);
      }
    }
  }//left_trunc



} //compute_ti


void pspline_multibaseline::compute_int_ti_vc_di0(const vector<double *> splinevector,const vector<double *> splinevector2,const vector<double *> betavector)
{
//ofstream oftest_spline ("f:\\baseline\\spline_di0.txt");
double * int_D_help;
double * betap;
double dist_knots=int_knots(1,0)-int_knots(0,0);
unsigned i,j,k,i_help,i_vc;
i=0;
i_help=0;
double erg,spline_u,spline_o,spline_ti_help,spline_ti;
double * int_ti_p=likep->get_integral_ti()+col;
double * int_ti_help_p=int_ti_help.getV();
double * int_ti_help_p2=int_ti_help.getV();
double * spline_help;
double * spline2;
statmatrix<double*>z_vc_help;
z_vc_help = statmatrix<double*>(baselinep.size()-1,1);
for(i_vc=1;i_vc<baselinep.size();i_vc++)
  z_vc_help(i_vc-1,0) = baselinep[i_vc]->get_z_vc();
k=1;
erg=0.0;
spline_o=0.0;
spline_u=0.0;
spline_help = splinevector[0];
spline2 = splinevector2[0];

int_D_help =baselinep[0]->get_int_D();
betap = betavector[0];

for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
  spline_o += *betap* *int_D_help;
spline_u=spline_o;
//------i0 berechnen---------
double z_vc_sum =0.0;
for(i_vc=1;i_vc<baselinep.size();i_vc++)
  z_vc_sum = z_vc_sum + *(z_vc_help(i_vc-1,0) + ges_index(i,0));
while(z_vc_sum!=0.0 && i < ges_index.rows()-1)
  {
  i++;
  spline_help++;
  z_vc_sum =0.0;
  for(i_vc=1;i_vc<baselinep.size();i_vc++)
    {
//    double test1= *(z_vc_help(i_vc-1,0) + index(i,0));
//    double test2 = index(i,0);
    z_vc_sum = z_vc_sum + *(z_vc_help(i_vc-1,0) + ges_index(i,0));
    }
  }

i_help=i;
//------------------erster Integralwert------------------------

while(k<int_knots.rows() && int_knots(k,0)<=zi_ges(ges_index(i,0),0) )
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
erg=erg+(exp(spline_ti)+exp(spline_o))*(zi_ges(ges_index(i,0),0)-int_knots(k-1,0));

int_ti_p=likep->get_integral_ti()+ges_index(i,0)*likep->get_responsedim()+col;
*int_ti_p =erg*0.5/(exp(spline_ti));

int_ti_help_p=int_ti_help.getV()+ges_index(i,0);
*int_ti_help_p =erg*0.5;
//------------------------------------------------------------

i++;
spline_help++;
while(i<zi_ges.rows())
  {
  z_vc_sum =0.0;
  for(i_vc=1;i_vc<baselinep.size();i_vc++)
    z_vc_sum = z_vc_sum + *(z_vc_help(i_vc-1,0) + ges_index(i,0));
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
    if(k<int_knots.rows() && zi_ges(ges_index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi_ges(ges_index(i,0),0)-zi_ges(ges_index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_ti));
    else
      {
      spline_u=spline_o;
      spline_o=0.0;
      betap = betavector[0];
      for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
        spline_o += *betap* *int_D_help;
      erg=erg+(int_knots(k,0)-zi_ges(ges_index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_o)) ;

      k++;

      while(k<int_knots.rows() && int_knots(k,0)<=zi_ges(ges_index(i,0),0) )
        {
        spline_u=spline_o;
        spline_o=0.0;
        betap = betavector[0];
        for(j=0;j<baselinep[0]->get_nrpar();j++,int_D_help++,betap++)
          spline_o += *betap* *int_D_help;
        erg=erg+dist_knots*(exp(spline_u)+exp(spline_o));

        k++;
        }
      erg=erg+(exp(spline_ti)+exp(spline_o))*(zi_ges(ges_index(i,0),0)-int_knots(k-1,0));
      }

    int_ti_p=likep->get_integral_ti()+ges_index(i,0)*likep->get_responsedim()+col;
    *int_ti_p =erg*0.5/(exp(spline_ti));

    int_ti_help_p=int_ti_help.getV()+ges_index(i,0);
    *int_ti_help_p =erg*0.5;

    i_help=i;
    i++;
    spline_help++;
    spline_ti_help=spline_ti;
    }//--------------- else: d.h. z_vc(index(i,0),0)==0 -----------------------
  } //while
//oftest_spline.close();
//------------------------------------------------------------------------------

  i=0;
  for(i=likep->get_nrobs();i<2*likep->get_nrobs();i++)
    {
    z_vc_sum=0.0;
    for(i_vc=1;i_vc<baselinep.size();i_vc++)
      z_vc_sum = z_vc_sum + *(z_vc_help(i_vc-1,0)+i);
    if(zi_ges(i,0)!=0 && z_vc_sum==0)
      {
      int_ti_p=likep->get_integral_ti()+(i-likep->get_nrobs())*likep->get_responsedim()+col;
      int_ti_help_p=int_ti_help.getV()+i-likep->get_nrobs();
      int_ti_help_p2=int_ti_help.getV()+i;
      *int_ti_p = (*int_ti_help_p-*int_ti_help_p2)/exp(*(spline2+i-likep->get_nrobs()));
      assert(*int_ti_p>=0.0);
      }
    }
}



void pspline_multibaseline::compute_int_ti_vc_di(const int dummy, const vector<double *> splinevector,const vector<double *> splinevector2,const vector<double *> betavector)
{
double * betap;
double dist_knots=int_knots(1,0)-int_knots(0,0);
unsigned i,j,k,i_help;
i=0;
i_help=0;
double erg,spline_u,spline_o,spline_ti_help,spline_ti;
double * int_ti_p=likep->get_integral_ti()+col;
double * int_ti_help_p=int_ti_help.getV();
double * int_ti_help_p2=int_ti_help.getV();
double * z_vc_help;
spline_ti=0.0;
spline_ti_help=0.0;
statmatrix<double*> int_D_help_1;
int_D_help_1=statmatrix<double*>(2,1);
statmatrix<double*> spline_zi_help;
spline_zi_help=statmatrix<double*>(2,1);
statmatrix<double*> spline_zi2_help;
spline_zi2_help=statmatrix<double*>(2,1);
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
spline_zi2_help(0,0) = splinevector2[0];
spline_zi2_help(1,0) = splinevector2[dummy];

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
while(*(z_vc_help+ges_index(i,0))==0.0)
  {
  i++;
  for(i_vc=0;i_vc<2;i_vc++)
    {
    spline_zi_help(i_vc,0)++;
    }
  }
i_help=i;

//------------------erster Integralwert------------------------
while(k<int_knots.rows() && int_knots(k,0)<=zi_ges(ges_index(i,0),0) )
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
erg=erg+(exp(spline_ti)+exp(spline_o))*(zi_ges(ges_index(i,0),0)-int_knots(k-1,0));

int_ti_p=likep->get_integral_ti()+ges_index(i,0)*likep->get_responsedim()+col;
*int_ti_p =erg*0.5/(exp(spline_ti));

int_ti_help_p=int_ti_help.getV()+ges_index(i,0);
*int_ti_help_p =erg*0.5;

//------------------------------------------------------------

i++;
for(i_vc=0;i_vc<2;i_vc++)
  spline_zi_help(i_vc,0)++;
while(i<zi_ges.rows())
  {
//-----------falls z_vc==0: zu nächster Beobachtung---------------
  if(*(z_vc_help+ges_index(i,0))==0)
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
    if(k<int_knots.rows() && zi_ges(ges_index(i,0),0)<=int_knots(k,0))
      erg=erg+(zi_ges(ges_index(i,0),0)-zi_ges(ges_index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_ti)) ;
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

      erg=erg+(int_knots(k,0)-zi_ges(ges_index(i_help,0),0))*(exp(spline_ti_help)+exp(spline_o)) ;

      k++;

      while(k<int_knots.rows() && int_knots(k,0)<=zi_ges(ges_index(i,0),0) )
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

      erg=erg+(exp(spline_ti)+exp(spline_o))*(zi_ges(ges_index(i,0),0)-int_knots(k-1,0));
      }

    int_ti_p=likep->get_integral_ti()+ges_index(i,0)*likep->get_responsedim()+col;
    *int_ti_p =erg*0.5/(exp(spline_ti));

    int_ti_help_p=int_ti_help.getV()+ges_index(i,0);
    *int_ti_help_p =erg*0.5;

    i_help=i;
    spline_ti_help=spline_ti;
    i++;
    for(i_vc=0;i_vc<2;i_vc++)
      spline_zi_help(i_vc,0)++;
    }//--------------- else: d.h. z_vc(index(i,0),0)==1 -----------------------
  } //while
  //------------------------------------------------------------------------------

  i=0;
  for(i=likep->get_nrobs();i<2*likep->get_nrobs();i++)
    {
    if(zi_ges(i,0)!=0 && *(z_vc_help+i-likep->get_nrobs())==1)
      {
      spline_ti=0.0;
      for(i_vc=0;i_vc<2;i_vc++)
        {
        help= *(spline_zi2_help(i_vc,0)+i-likep->get_nrobs());
        spline_ti +=help;
        }
      int_ti_p=likep->get_integral_ti()+(i-likep->get_nrobs())*likep->get_responsedim()+col;
      int_ti_help_p=int_ti_help.getV()+i-likep->get_nrobs();
      int_ti_help_p2=int_ti_help.getV()+i;
      *int_ti_p = (*int_ti_help_p-*int_ti_help_p2)/exp(spline_ti);
      assert(*int_ti_p>=0.0);
      }
    }
}//--------------compute_int_ti_vc_di---------------------









void pspline_multibaseline::update_multibaseline()
{
//---------Integral berechnen---------------------------------
    unsigned i;

    if(baselinep.size()>1)
      {
      vector <double *> splinevec;
      vector <double *> splinevec2;
      vector <double *> betavec;
      for(i=0;i<(baselinep.size());i++)
        {
        splinevec.push_back(baselinep[i]->get_spline_ges());
        splinevec2.push_back(baselinep[i]->get_spline_ges2());
        }
      for(i=0;i<(baselinep.size());i++)
        betavec.push_back(baselinep[i]->getbetapointer());

      compute_int_ti_vc_di0(splinevec,splinevec2,betavec);
      for(i=1;i<baselinep.size();i++)
        {
        compute_int_ti_vc_di(i,splinevec,splinevec2,betavec);
        }
      }
    else
      {
      if(global == true)
        {
        testmat.mult(spline_ges,beta);
        testmat.mult_index(spline_ges2,beta);
        compute_int_ti_global(beta);
        }
      else
        {
        testmat_l[col].mult(spline_teil[col],beta);
        testmat_l[col].mult_index(spline_teil2[col],beta);
        compute_int_ti_nonglobal(beta);
        }
      }
}
//------------------------------------------------------------------------------









} // end: namespace MCMC






