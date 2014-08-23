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



#include "FC_hrandom_variance_vec_nmig.h"


//------------------------------------------------------------------------------
//--- CLASS: FC_hrandom_variance_vec_nmig implementation of member functions ---
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_hrandom_variance_vec_nmig::read_options(vector<ST::string> & op,
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

  f = op[39].strtodouble(abeta);
  f = op[40].strtodouble(bbeta);
  f = op[41].strtodouble(r);
  f = op[42].strtodouble(v);
  f = op[43].strtodouble(aQ);
  f = op[44].strtodouble(bQ);
  f = op[45].strtolong(regiterates);
  }


FC_hrandom_variance_vec_nmig::FC_hrandom_variance_vec_nmig(void)
  {

  }


FC_hrandom_variance_vec_nmig::FC_hrandom_variance_vec_nmig(MASTER_OBJ * mp,unsigned & enr,
                  GENERAL_OPTIONS * o,DISTR * lp,
                  DISTR * lpRE,
                  const ST::string & t,const ST::string & fp,
                  DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                  vector<ST::string> & vn)
     : FC_hrandom_variance_vec(mp,enr,o,lp,lpRE,t,fp,Dp,FCn,op,vn)
  {
  read_options(op,vn);

  FC_delta = FC(o,"",beta.rows(),1,"");
  FC_delta.setbeta(beta.rows(),1,1);

  FC_omega = FC(o,"",1,1,"");
  FC_omega.setbeta(1,1,0.5);

  FC_Q = FC(o,"",1,1,"");
  FC_Q.setbeta(1,1,0.1);

  }


FC_hrandom_variance_vec_nmig::FC_hrandom_variance_vec_nmig(
   const FC_hrandom_variance_vec_nmig & m)
  : FC_hrandom_variance_vec(FC_hrandom_variance_vec(m))
  {
  FC_delta = m.FC_delta;
  FC_omega = m.FC_omega;
  FC_Q = m.FC_Q;
  abeta = m.abeta;
  bbeta = m.bbeta;
  r = m.r;
  v = m.v;
  aQ=m.aQ;
  bQ=m.aQ;
  regiterates = m.regiterates;
  }


const FC_hrandom_variance_vec_nmig & FC_hrandom_variance_vec_nmig::operator=(
const FC_hrandom_variance_vec_nmig & m)
  {

  if (this==&m)
	 return *this;
  FC_hrandom_variance_vec::operator=(FC_hrandom_variance_vec(m));
  FC_delta = m.FC_delta;
  FC_omega = m.FC_omega;
  FC_Q = m.FC_Q;
  abeta = m.abeta;
  bbeta = m.bbeta;
  r = m.r;
  v = m.v;
  aQ=m.aQ;
  bQ=m.aQ;
  regiterates = m.regiterates;
  return *this;
  }



void FC_hrandom_variance_vec_nmig::update(void)
  {

  register unsigned i;
  double * workbeta = beta.getV();               // psi
  double * workbetafcn = FCnonpp->beta.getV();   // beta, i.e. random effects
  double * workdelta = FC_delta.beta.getV();       // deltas

  double * linpredREp;
  if (likepRE->linpred_current==1)
    linpredREp = likepRE->linearpred1.getV();
  else
    linpredREp = likepRE->linearpred2.getV();

  double * ww = likepRE->workingweight.getV();

  double nup05 = v+0.5;
  double rdelta;
  double L;
  double sumdelta=0;
  double suminvpsi = 0;
  double diff2;
  double rsqrt = sqrt(r);
  double pr;
  double omega;

  for (i=0;i<beta.rows();i++,workbeta++,workbetafcn++,ww++,workdelta++,
                         linpredREp++)
    {

    if (*workdelta==1)
      rdelta = 1;
    else
      rdelta = r;

    // update psi's

    diff2 = pow((*workbetafcn) - (*linpredREp),2);

    *workbeta = rand_invgamma(nup05,FC_Q.beta(0,0)
                + 0.5*diff2/rdelta);
    suminvpsi += 1/(*workbeta);

    // update delta's

    L = 1/rsqrt*exp(1/(2*(*workbeta))*diff2*(1-1/r)) ;
    omega = FC_omega.beta(0,0);
    pr = 1/(1+ (1-omega)/omega*L);

    *workdelta = randnumbers::bernoulli(pr);

    sumdelta+=(*workdelta);


    if (*workdelta == 0)
      *workbeta *= r;

    *ww = 1/(*workbeta);
    }


  // update omega
  FC_omega.beta(0,0) = randnumbers::rand_beta(abeta+sumdelta,bbeta+beta.rows()-sumdelta);

  // update Q
  FC_Q.beta(0,0) = rand_gamma(aQ+beta.rows()*v,bQ+suminvpsi);
//  FC_Q.beta(0,0) = randnumbers::GIG(v*beta.rows()-aQ,2*suminvpsi,2*bQ);


  FCnonpp->tau2 = 1;
  designp->compute_penalty2(beta);

  likepRE->sigma2 = 1;

  acceptance++;

  FC_delta.update();
  FC_omega.update();
  FC_Q.update();
  FC::update();
  }


void FC_hrandom_variance_vec_nmig::outresults(ofstream & out_stata,
                                              ofstream & out_R,
                                              const ST::string & pathresults)
  {


  if (pathresults.isvalidfile() != 1)
    {

    ST::string pathresults_delta = pathresults.substr(0,pathresults.length()-4) + "_delta.res";
    ST::string pathresults_omega = pathresults.substr(0,pathresults.length()-4) + "_omega.res";
    ST::string pathresults_Q = pathresults.substr(0,pathresults.length()-4) + "_Q.res";

    FC_nonp_variance_vec::outresults(out_stata,out_R,pathresults);

    FC_delta.outresults(out_stata,out_R,"");
    FC_omega.outresults(out_stata,out_R,pathresults_omega);
    FC_Q.outresults(out_stata,out_R,pathresults_Q);

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

    optionsp->out("    Variance parameter Q:\n");
    optionsp->out("\n");
    FC_Q.outresults_singleparam(out_stata,out_R,"");
    optionsp->out("    Results for the variance parameter Q are also stored in file\n");
    optionsp->out("    " +  pathresults_Q + "\n");
    optionsp->out("\n");

    // deltas
    ofstream ou(pathresults_delta.strtochar());

    ou << "intnr   " << designp->datanames[0]  << "   pmean" << endl;

    unsigned i;
    double * deltameanp = FC_delta.betamean.getV();

    for (i=0;i<beta.rows();i++,deltameanp++)
      ou << (i+1) << "  " << designp->effectvalues[i]  << "  " << (*deltameanp) << endl;
    }

  }


void FC_hrandom_variance_vec_nmig::get_samples(
   const ST::string & filename,ofstream & outg) const
  {
  FC_hrandom_variance_vec::get_samples(filename,outg);

  ST::string filename_delta = filename.substr(0,filename.length()-4) + "_delta.raw";
  FC_delta.get_samples(filename_delta,outg);

  ST::string filename_omega = filename.substr(0,filename.length()-4) + "_omega.raw";
  FC_omega.get_samples(filename_omega,outg);

  ST::string filename_Q = filename.substr(0,filename.length()-4) + "_Q.raw";
  FC_Q.get_samples(filename_Q,outg);
  }


void FC_hrandom_variance_vec_nmig::compute_autocorr_all(const ST::string & path,
                                           unsigned lag, ofstream & outg) const
  {
  FC_hrandom_variance_vec::compute_autocorr_all(path,lag,outg);

//  ST::string path_delta = path.substr(0,path.length()-4) + "_delta.raw";
//  FC_delta.compute_autocorr_all(path_delta,lag,outg);

  ST::string path_omega = path.substr(0,path.length()-4) + "_omega.raw";
  FC_omega.compute_autocorr_all(path_omega,lag,outg);

  ST::string path_Q = path.substr(0,path.length()-4) + "_Q.raw";
  FC_Q.compute_autocorr_all(path_Q,lag,outg);
  }


} // end: namespace MCMC




