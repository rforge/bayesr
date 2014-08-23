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



#include "mcmc_pspline.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//---------- CLASS: FULLCOND_pspline (implementation of member functions) ------
//------------------------------------------------------------------------------


void FULLCOND_pspline::make_Kab_list(void)
  {

  unsigned i,j;
  unsigned nrupdate;
  unsigned a,b;
  datamatrix kab;
  datamatrix help;
  datamatrix null(1,1,0);

  unsigned sum=0;
  for(i=min;i<=max;i++)
    {
    nrupdate = K.get_rows()/i;
    if ((nrupdate*i) < K.get_rows())
      nrupdate++;
    sum+= nrupdate;
    }

  KAB.reserve(sum);
  KABroot.reserve(sum);
  KABr_sp.reserve(sum);
  KABl_sp.reserve(sum);

  for (i=min;i<=max;i++)
    {

    nrupdate = K.get_rows()/i;
    if ((nrupdate*i) < K.get_rows())
      nrupdate++;

    begin.push_back(KAB.size());
    matquant.push_back(nrupdate);

    for(j=1;j <= nrupdate;j++)
      {

      a = 1+(j-1)*i;

      if (j == nrupdate)
        b = K.get_rows();
      else
        b = j*i;

      kab = (K.getBlock(a-1,a-1,b,b)).inverse();
      KAB.push_back(kab);
      KABroot.push_back(kab.root());

      if (b != K.get_cols())
        {
        help = kab*K.getBlock(a-1,b,b,K.get_cols());
        KABr_sp.push_back(SparseMatrix(help));
        }
      else
        KABr_sp.push_back(SparseMatrix());

      if (a != 1)
        {
        help = kab*K.getBlock(a-1,0,b,a-1);
        KABl_sp.push_back(SparseMatrix(help));
        }
      else
        KABl_sp.push_back(SparseMatrix());

      } // end: for(j=1;j <= nrupdate;j++)

    } // end: for (i=min;i<=max;i++)

  } // end: function make_Kab_list


void FULLCOND_pspline::compute_mu(const datamatrix & beta,
      const unsigned & bs,const unsigned & a,const unsigned & b,
      const unsigned & v)
  {
  unsigned matnr = begin[bs-min]+ ((a-1)/bs);

  if (a==1)
    KABr_sp[matnr].substr_mult(beta,b,v,fc_random[b-a],0);
  else if (b == nrpar)
    KABl_sp[matnr].substr_mult(beta,0,v,fc_random[b-a],0);
  else
    {
    KABr_sp[matnr].substr_mult(beta,b,v,fc_random[b-a],0);
    KABl_sp[matnr].substr_mult(beta,0,v,fc_random[b-a],0);
    }
  }


void FULLCOND_pspline::compute_fc(const datamatrix & beta, const unsigned & bs,
                               const unsigned & a,const unsigned & b,
                               const double & Q,const unsigned & v)
  {
  unsigned matnr = begin[bs-min]+ ((a-1)/bs );
  unsigned l = b-a+1;
  register unsigned i,j;

  double * fc_randwork = fc_random[b-a].getV();
  double * KABrootwork = KABroot[matnr].getV();

  double * randnormwork = randnorm[b-a].getV();

  for(i=0;i<l;i++,randnormwork++)
    *randnormwork = rand_normal();

  for (i=0;i<l;i++,fc_randwork++)
    {
    *fc_randwork = 0;
    randnormwork = randnorm[b-a].getV();
    for(j=0;j<l;j++,KABrootwork++,randnormwork++)
      *fc_randwork += *KABrootwork * *randnormwork;

    *fc_randwork *= Q;
    }

  compute_mu(beta,bs,a,b,v);

  }


FULLCOND_pspline::FULLCOND_pspline(MCMCoptions * o,DISTRIBUTION * dp,
                    FULLCOND_const * fcc,const datamatrix & d,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const double & lk, const double & uk,
                    const double & lg, const double & ug,
                    const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,degr,kp,gs,fp,pres,deriv,lk,uk,lg,ug,c)
  {

  oldacceptance = 0;
  oldnrtrials = 0;

  lambda = l;
  sigma2 = 1.0/l;

  min = minb;
  max = maxb;
  mintoobig = false;
  maxtoobig = false;

  unsigned i;

  varcoeff = false;

  compute_betaweight();

  make_index(d);
  make_Bspline(d);

  compute_Kweights();

  if(predictright || predictleft)
    {
    datamatrix help = betaweight;
    betaweight = datamatrix(nrpar,1,0);
    for(i=0;i<nrparpredictleft;i++)
      betaweight(i,0) = 0.0;
    for(;i<nrpar-nrparpredictright;i++)
      betaweight(i,0) = help(i-nrparpredictleft,0);
    for(;i<nrpar;i++)
      betaweight(i,0) = 0.0;
    }

  init_fchelp(d);

  if (type == RW1)
    {
    K = Krw1(weight);
    Kenv = Krw1env(weight);
    rankK = nrpar-nrparpredictleft-nrparpredictright-1;
    }
  else if (type == RW2)
    {
    K = Krw2(weight);
    Kenv = Krw2env(weight);
    rankK = nrpar-nrparpredictleft-nrparpredictright-2;
    }

  if(predictleft || predictright)
    change_K();

  betaold = datamatrix(nrpar,1,0); // wegen "change" bei Interaktion

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

//  compute_betaweight();

  }


FULLCOND_pspline::FULLCOND_pspline(MCMCoptions * o,DISTRIBUTION * dp,
                    FULLCOND_const * fcc,
                    const datamatrix & effmod,const datamatrix & intact,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,degr,kp,gs,fp,pres,deriv,0.0,0.0,0.0,0.0,c)
  {

  oldacceptance = 0;
  oldnrtrials = 0;

  lambda = l;
  sigma2 = 1.0/l;

  min = minb;
  max = maxb;
  mintoobig = false;
  maxtoobig = false;

  unsigned i;

  varcoeff = true;

  compute_betaweight();

  make_index(effmod,intact);
  make_Bspline(effmod);
  make_BS(intact);

  compute_Kweights();

  if(predictright || predictleft)
    {
    datamatrix help = betaweight;
    betaweight = datamatrix(nrpar,1,0);
    for(i=0;i<nrparpredictleft;i++)
      betaweight(i,0) = 0.0;
    for(;i<nrpar-nrparpredictright;i++)
      betaweight(i,0) = help(i-nrparpredictleft,0);
    for(;i<nrpar;i++)
      betaweight(i,0) = 0.0;
    }

  init_fchelp(effmod);

  if (type == RW1)
    {
    K = Krw1(weight);
    rankK = nrpar-nrparpredictleft-nrparpredictright-1;
    }
  else if (type == RW2)
    {
    K = Krw2(weight);
    rankK = nrpar-nrparpredictleft-nrparpredictright-2;
    }

  if(predictleft || predictright)
    change_K();

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

//  compute_betaweight();

  }


FULLCOND_pspline::FULLCOND_pspline(MCMCoptions * o,DISTRIBUTION * dp,
                    FULLCOND_const * fcc,const fieldtype & ft,
                    const ST::string & ti,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,degr,kp,gs,fp,pres,deriv,0.0,0.0,0.0,0.0,c)
  {
  }


FULLCOND_pspline::FULLCOND_pspline(const FULLCOND_pspline & fc)
  :spline_basis(spline_basis(fc))
  {

  min = fc.min;
  max = fc.max;
  minauto = fc.minauto;
  maxauto = fc.maxauto;
  automatic = fc.automatic;
  mintoobig = fc.mintoobig;
  maxtoobig = fc.maxtoobig;
  oldacceptance = fc.oldacceptance;
  oldnrtrials = fc.oldnrtrials;

  K = fc.K;
  fc_random = fc.fc_random;
  randnorm = fc.randnorm;
  KAB = fc.KAB;
  KABl_sp = fc.KABl_sp;
  KABr_sp = fc.KABr_sp;
  KABroot = fc.KABroot;
  begin = fc.begin;
  matquant = fc.matquant;

  }


const FULLCOND_pspline & FULLCOND_pspline::operator=(const FULLCOND_pspline & fc)
  {
  if (this == &fc)
    return *this;
  spline_basis::operator=(spline_basis(fc));

  min = fc.min;
  max = fc.max;
  minauto = fc.minauto;
  maxauto = fc.maxauto;
  automatic = fc.automatic;
  mintoobig = fc.mintoobig;
  maxtoobig = fc.maxtoobig;
  oldacceptance = fc.oldacceptance;
  oldnrtrials = fc.oldnrtrials;

  K = fc.K;
  fc_random = fc.fc_random;
  randnorm = fc.randnorm;
  KAB = fc.KAB;
  KABl_sp = fc.KABl_sp;
  KABr_sp = fc.KABr_sp;
  KABroot = fc.KABroot;
  begin = fc.begin;
  matquant = fc.matquant;

  return *this;
  }


void FULLCOND_pspline::outoptions(void)
  {

  optionsp->out("  OPTIONS FOR P-SPLINE TERM: " + title + "\n",true);

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




void FULLCOND_pspline::update(void)
  {

  unsigned blocksize;

  if(automatic)
    {
    if(optionsp->get_nriter()%100==0 && optionsp->get_nriter()<optionsp->get_burnin())
      adjust_blocksize(30,70);
    if(optionsp->get_nriter() == optionsp->get_burnin())
      {
      optionsp->out("\n");
      optionsp->out("NOTE: Minimum blocksize for " + title + " set to " + ST::inttostring(minauto) + "\n");
      optionsp->out("NOTE: Maximum blocksize for " + title + " set to " + ST::inttostring(maxauto) + "\n");
      optionsp->out("\n");
      }
#ifndef __BUILDING_GNU
    blocksize = minauto + random(maxauto-minauto+1);
#else
    blocksize = minauto + int((maxauto-minauto+1)*rand()/(RAND_MAX + 1.0));
#endif
    }
  else
    {
#ifndef __BUILDING_GNU
    blocksize = min + random(max-min+1);
#else
    blocksize = min + int((max-min+1)*rand()/(RAND_MAX + 1.0));
#endif
    }

  unsigned i;

  double u;
  unsigned an = 1;
  unsigned en = blocksize;

  unsigned beg;
  unsigned end;

  double logold;
  double logprop;
  double * workbeta;
  double * workfcrand;
  unsigned j,k;

  if(lambdaconst == true)
    sigma2 = likep->get_scale(column)/lambda;

  for(j=0;j<matquant[blocksize-min];j++)
    {

    nrtrials++;

    compute_fc(beta,blocksize,an,en,sqrt(sigma2));

    logold = 0;
    logprop = 0;

    beg = firstnonzero[an-1];
    end = lastnonzero[en-1];

    logold += likep->loglikelihood(beg,end,index);

    likep->assign(false);

    add_linearpred_multBS_Block(an-1,en-1,fc_random[en-an]);

    logprop += likep->loglikelihood(beg,end,index,false);

    u = log(uniform());

    if (u <= (logprop-logold))
      {
      workbeta = beta.getV()+an-1;          // change beta
      workfcrand = fc_random[en-an].getV();
      for(k=an-1;k<en;k++,workbeta++,workfcrand++)
        *workbeta = *workfcrand;
	  acceptance++;
      likep->swap_linearpred();
      }

    an+=blocksize;
    if (j == matquant[blocksize-min]-2)
      en = nrpar;
    else
      en+=blocksize;

    } // end: for(j=0;j<matquant[blocksize-min];j++)

  if(predictright || predictleft)
    update_prediction();

  if (center)
    {
    compute_intercept();

    workbeta = beta.getV();
    for(i=0;i<nrpar;i++,workbeta++)
      *workbeta -= intercept;

//    likep->add_linearpred_m(-intercept,column);
    fcconst->update_intercept(intercept);
    }

  if(interaction == false)
    {

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

    }  // end: if(interaction == false)

  } // end: function update


void FULLCOND_pspline::outresults(void)
  {
  FULLCOND::outresults();
  spline_basis::outresults();
  }


void FULLCOND_pspline::predict(const datamatrix & newX, datamatrix & linpred)
  {

  unsigned i,j;
  datamatrix betac(beta.rows(),beta.cols());
  datamatrix bspline(1,nrpar,0);
  double * worklin = linpred.getV();

  for(i=0;i<nrpar;i++)
    bspline(0,i) = bspline_rek(degree,i,newX);

  if(varcoeff)
    {
    double help;
    for(i=0;i<optionsp->get_samplesize();i++,worklin++)
      {
      help = 0.0;
      readsample2(betac,i);
      for(j=0;j<nrpar;j++)
        {
        help += betac(j,0) * bspline(0,j);
        }
      *worklin += help * newX(0,1);
      }
    }
  else
    {
    for(i=0;i<optionsp->get_samplesize();i++,worklin++)
      {
      readsample2(betac,i);
      for(j=0;j<nrpar;j++)
        {
        *worklin += betac(j,0) * bspline(0,j);
        }
      }
    }

  }


double FULLCOND_pspline::compute_quadform(void)
  {
  if(predictright || predictleft)
    {
    double quadform = K.compute_quadform(beta,0);
    if(predictright)
      for(unsigned i=0;i<nrparpredictright;i++)
        quadform -= beta(nrpar-1-i,0)*beta(nrpar-1-i,0);
    if(predictleft)
      for(unsigned i=0;i<nrparpredictleft;i++)
        quadform -= beta(i,0)*beta(i,0);
    return quadform;
    }
  else
    return K.compute_quadform(beta,0);
  }


void FULLCOND_pspline::confidencebands(const ST::string & filename)
  {

   double min = data(index(0,0),0);
   double max = data(index(likep->get_nrobs()-1,0),0);
   double dist = (max - min)/99.0;
   datamatrix value(1,1);

   ofstream confidout(filename.strtochar());

   confidout << "intnr" << "   ";
   confidout << title << "   ";
   confidout << "pmean   ";
   confidout << "pqu10   ";
   confidout << "pmed   ";
   confidout << "pqu90   ";
   confidout << endl;

   for(unsigned i=0;i<100;i++,min+=dist)
     {
     value(0,0) = data(index(i,0),0);
     datamatrix linpred(optionsp->get_samplesize(),1,0);

     predict(value,linpred);

     confidout << (i+1) << "   ";
     confidout << value(0,0) << "   ";
     confidout << linpred.sum(0)/double(optionsp->get_samplesize()) << "   ";
     confidout << linpred.quantile(10,0) << "   ";
     confidout << linpred.quantile(50,0) << "   ";
     confidout << linpred.quantile(90,0) << "   ";
     confidout << endl;
     }

  }

double FULLCOND_pspline::compute_mse(const datamatrix & m)
  {

  double mse = 0.0;
  datamatrix help(m.rows(),1);
  datamatrix mean = help;

  multBS(mean,betamean);
  help.minus(m,mean);

  mse = norm(help,0)*norm(help,0);
  return mse/double(m.rows());

  }


void FULLCOND_pspline::adjust_blocksize(const unsigned & alphamin,const unsigned & alphamax)
  {

  int min = minauto;
  int max = maxauto;

  double rate;
  if (nrtrials == 0)
//    rate = (double(acceptance)/double(optionsp->get_nriter()) )*100;
    rate = (double(acceptance-oldacceptance)/double(100))*100;
  else
//    rate = (double(acceptance)/double(nrtrials))*100;
    rate = (double(acceptance-oldacceptance)/double(nrtrials-oldnrtrials))*100;

  oldacceptance = acceptance;
  oldnrtrials = nrtrials;

  int limit = 1;
  int span = (nrpar/10)>1?(nrpar/10):2;

  if(rate<alphamin)
    {
    if(max-min<span)
      {
      if(rate<(alphamin-15))
        min -= span;
      else
        min--;
      if(min<limit)
        min = limit;
//      if(min<minauto)
//        optionsp->out("NOTE: Minimum blocksize for " + title + " changed: " + ST::inttostring(min) + "\n");
      }
    else
      {
      if(rate<(alphamin-15))
        max -= span;
      else
        max--;
      if(max<min)
        max = min;
//      if(max<maxauto)
//        optionsp->out("NOTE: Maximum blocksize for " + title + " changed: " + ST::inttostring(max) + "\n");
      }
    }

  if(rate>alphamax)
    {
    if(max-min<span)
      {
      if(rate>(alphamax+15))
        max += span;
      else
        max++;
      if(max>rankK)
        max = rankK;
//      if(max>maxauto)
//        optionsp->out("NOTE: Maximum blocksize for " + title + " changed: " + ST::inttostring(max) + "\n");
      }
    else
      {
      if(rate>(alphamax+15))
        min += span;
      else
        min++;
      if(min>max)
        min = max;
//      if(min>minauto)
//        optionsp->out("NOTE: Minimum blocksize for " + title + " changed: " + ST::inttostring(min) + "\n");
      }
    }

  minauto = min;
  maxauto = max;

  }


} // end: namespace MCMC










