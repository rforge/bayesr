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



#include "FC_hrandom_variance.h"


//------------------------------------------------------------------------------
//--------- CLASS: FC_hrandom_variance implementation of member functions ------
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_hrandom_variance::read_options(vector<ST::string> & op,
vector<ST::string> & vn)
  {

  int f;

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  7       center
  8       map
  9       lambda_re
  10      a_re
  11      b_re
  12      internal_mult
  */

  FC_nonp_variance::read_options(op,vn);

  if (op[12] == "true")
    {
    mult =true;
    f = op[9].strtodouble(lambdastart);
    f = op[10].strtodouble(a_invgamma);
    f = op[11].strtodouble(b_invgamma);
    }
  else
    {
    mult = false;
    }

  if (op[38] == "true")
    {
    lambdaconst = true;
    nosamples =true;
    }

  }


FC_hrandom_variance::FC_hrandom_variance(void)
  {

  }


FC_hrandom_variance::FC_hrandom_variance(MASTER_OBJ * mp,unsigned & enr,
                 GENERAL_OPTIONS * o,DISTR * lp,
                  DISTR * lpRE,
                 const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                 vector<ST::string> & vn)
     : FC_nonp_variance(mp,enr,o,lp,t,fp,Dp,FCn,op,vn)
  {
  read_options(op,vn);
  likepRE = lpRE;
  }


FC_hrandom_variance::FC_hrandom_variance(const FC_hrandom_variance & m)
  : FC_nonp_variance(FC_nonp_variance(m))
  {
  likepRE = m.likepRE;
  mult = m.mult;
  }


const FC_hrandom_variance & FC_hrandom_variance::operator=(
const FC_hrandom_variance & m)
  {

  if (this==&m)
	 return *this;
  FC_nonp_variance::operator=(FC_nonp_variance(m));
  likepRE = m.likepRE;
  mult = m.mult;
  return *this;
  }


double FC_hrandom_variance::compute_quadform(void)
  {

  unsigned n;
  double sum = 0;
  double * workbeta = FCnonpp->beta.getV();
  register unsigned i;

  n = FCnonpp->beta.rows();


  double * linpredREp;
  if (likepRE->linpred_current==1)
    linpredREp = likepRE->linearpred1.getV();
  else
    linpredREp = likepRE->linearpred2.getV();

  for(i=0;i<n;i++,workbeta++,linpredREp++)
    {
    sum += pow(*workbeta-(*linpredREp),2);
    }

  return sum;

  }


void FC_hrandom_variance::update(void)
  {

  b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

  if (lambdaconst==false)
    {
    beta(0,0) = rand_invgamma(a_invgamma+0.5*designp->rankK,
                                    b_invgamma+0.5*compute_quadform());

    if (beta(0,0) > 5)
      beta(0,0) = 5;


    beta(0,1) = likep->get_scale()/beta(0,0);

    FCnonpp->tau2 = beta(0,0);
    likepRE->sigma2=beta(0,0);

    acceptance++;
    FC::update();
    }

  }


bool FC_hrandom_variance::posteriormode(void)
  {
  return  FC_nonp_variance::posteriormode();
  }


void FC_hrandom_variance::outresults(ofstream & out_stata, ofstream & out_R,
                         const ST::string & pathresults)
  {
  if (lambdaconst==false)
    FC_nonp_variance::outresults(out_stata,out_R,pathresults);
  }


//------------------------------------------------------------------------------
//---------------------- CLASS: FC_hrandom_variance_ssvs -----------------------
//------------------------------------------------------------------------------

void FC_hrandom_variance_ssvs::read_options(vector<ST::string> & op,
vector<ST::string> & vn)
  {

  int f;

  /*

  39 abeta
  40 bbeta
  41 r
  42 v
  43 aQ
  44 bQ

  */

  FC_nonp_variance::read_options(op,vn);
  FC_hrandom_variance::read_options(op,vn);

  f = op[39].strtodouble(abeta);
  f = op[40].strtodouble(bbeta);
  f = op[41].strtodouble(r);
  f = op[45].strtolong(regiterates);
  }


FC_hrandom_variance_ssvs::FC_hrandom_variance_ssvs(void)
  {

  }



FC_hrandom_variance_ssvs::FC_hrandom_variance_ssvs(MASTER_OBJ * mp,unsigned & enr,
                  GENERAL_OPTIONS * o,DISTR * lp,
                  DISTR * lpRE,
                  const ST::string & t,const ST::string & fp,
                  DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                  vector<ST::string> & vn)
     : FC_hrandom_variance(mp,enr, o,lp,lpRE,t,fp,Dp,FCn,op,vn)
  {
  read_options(op,vn);

  FC_delta = FC(o,"",FCn->beta.rows(),1,"");
  FC_delta.setbeta(FCn->beta.rows(),1,1);

  FC_omega = FC(o,"",1,1,"");
  FC_omega.setbeta(1,1,0.5);

  pen = datamatrix(FCn->beta.rows(),1,1);
  }


FC_hrandom_variance_ssvs::FC_hrandom_variance_ssvs(
   const FC_hrandom_variance_ssvs & m)
  : FC_hrandom_variance(FC_hrandom_variance(m))
  {
  FC_delta = m.FC_delta;
  FC_omega = m.FC_omega;
  abeta = m.abeta;
  bbeta = m.bbeta;
  r = m.r;
  regiterates = m.regiterates;
  pen = m.pen;
  }


const FC_hrandom_variance_ssvs & FC_hrandom_variance_ssvs::operator=(
const FC_hrandom_variance_ssvs & m)
  {

  if (this==&m)
	 return *this;
  FC_hrandom_variance::operator=(FC_hrandom_variance(m));
  FC_delta = m.FC_delta;
  FC_omega = m.FC_omega;
  abeta = m.abeta;
  bbeta = m.bbeta;
  r = m.r;
  regiterates = m.regiterates;
  pen = m.pen;
  return *this;
  }


void FC_hrandom_variance_ssvs::update(void)
  {

  register unsigned i;

  double * workbetafcn = FCnonpp->beta.getV();
  double * workdelta = FC_delta.beta.getV();

  double * linpredREp;
  if (likepRE->linpred_current==1)
    linpredREp = likepRE->linearpred1.getV();
  else
    linpredREp = likepRE->linearpred2.getV();

  double * ww = likepRE->workingweight.getV();

  double * penp = pen.getV();

  double L;
  double sumdelta=0;
  double diff2;
  double quadform = 0;
  double rsqrt = sqrt(r);
  double rinv = 1/r;
  double pr;

  double omega = FC_omega.beta(0,0);

  int nrcluster = FCnonpp->beta.rows();

  // update deltas
  for (i=0;i<nrcluster;i++,penp++,workbetafcn++,ww++,workdelta++,linpredREp++)
    {

    diff2 = pow((*workbetafcn) - (*linpredREp),2);

    L = 1/rsqrt*exp(1/(2*(beta(0,0)))*diff2*(1-1/r)) ;
    pr = 1/(1+ (1-omega)/omega*L);

    *workdelta = randnumbers::bernoulli(pr);

    sumdelta+=(*workdelta);

    if (*workdelta == 0)
      {
      *ww = rinv;
      *penp = r;
      quadform += rinv*diff2;
      }
    else
      {
      *ww = 1;
      *penp = 1;
       quadform += diff2;
      }
    }

  b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

  beta(0,0) = rand_invgamma(a_invgamma+0.5*nrcluster,b_invgamma+0.5*quadform);

  beta(0,1) = likep->get_scale()/beta(0,0);

  // update omega
  FC_omega.beta(0,0) = randnumbers::rand_beta(abeta+sumdelta,bbeta+nrcluster-sumdelta);

  FCnonpp->tau2 = beta(0,0);

  designp->compute_penalty2(pen);

  likepRE->sigma2 = beta(0,0);

  acceptance++;

  FC_delta.update();
  FC_omega.update();
  FC::update();

  }


void FC_hrandom_variance_ssvs::outoptions(void)
  {

  optionsp->out("  Options for spike and slap prior\n");
  optionsp->out("\n");

  FC_nonp_variance::outoptions();

  optionsp->out("  Hyperparameter a for beta distribution: " +
                ST::doubletostring(abeta) + "\n" );
  optionsp->out("  Hyperparameter b for beta distribution: " +
                ST::doubletostring(bbeta) + "\n" );

  optionsp->out("  Spike and slap parameter r: " + ST::doubletostring(r) + "\n");

  optionsp->out("\n");

  }


void FC_hrandom_variance_ssvs::outresults(ofstream & out_stata,
                                              ofstream & out_R,
                                              const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    ST::string pathresults_delta = pathresults.substr(0,pathresults.length()-4) + "_delta.res";
    ST::string pathresults_omega = pathresults.substr(0,pathresults.length()-4) + "_omega.res";

    FC_hrandom_variance::outresults(out_stata,out_R,pathresults);

    FC_delta.outresults(out_stata,out_R,"");
    FC_omega.outresults(out_stata,out_R,pathresults_omega);

    optionsp->out("    Results for the inclusion probabilities are stored in file\n");
    optionsp->out("    " +  pathresults_delta + "\n");
    optionsp->out("\n");
    optionsp->out("\n");

    optionsp->out("    Inclusion probability parameter omega:\n");
    optionsp->out("\n");
    FC_omega.outresults_singleparam(out_stata,out_R,"");
    optionsp->out("    Results for the inclusion probability parameter omega are also stored in file\n");
    optionsp->out("    " +  pathresults_omega + "\n");
    optionsp->out("\n");
    optionsp->out("\n");

    // deltas
    ofstream ou(pathresults_delta.strtochar());

    ou << "intnr   " << designp->datanames[0]  << "   pmean" << endl;

    unsigned i;
    double * deltameanp = FC_delta.betamean.getV();

    for (i=0;i<FC_delta.beta.rows();i++,deltameanp++)
      ou << (i+1) << "  " << designp->effectvalues[i]  << "  " << (*deltameanp) << endl;
    }

  }


void FC_hrandom_variance_ssvs::get_samples(
   const ST::string & filename,ofstream & outg) const
  {
  FC_hrandom_variance::get_samples(filename,outg);

  ST::string filename_delta = filename.substr(0,filename.length()-4) + "_delta.raw";
  FC_delta.get_samples(filename_delta,outg);

  ST::string filename_omega = filename.substr(0,filename.length()-4) + "_omega.raw";
  FC_omega.get_samples(filename_omega,outg);
  }


void FC_hrandom_variance_ssvs::compute_autocorr_all(const ST::string & path,
                                           unsigned lag, ofstream & outg) const
  {
  FC_hrandom_variance::compute_autocorr_all(path,lag,outg);

//  ST::string path_delta = path.substr(0,path.length()-4) + "_delta.raw";
//  FC_delta.compute_autocorr_all(path_delta,lag,outg);

  ST::string path_omega = path.substr(0,path.length()-4) + "_omega.raw";
  FC_omega.compute_autocorr_all(path_omega,lag,outg);
  }



} // end: namespace MCMC



