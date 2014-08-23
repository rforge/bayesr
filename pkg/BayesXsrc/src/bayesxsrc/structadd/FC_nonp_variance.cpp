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



#include "FC_nonp_variance.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//---------- CLASS: FC_nonp_variance implementation of member functions ---------
//------------------------------------------------------------------------------

void FC_nonp_variance::read_options(vector<ST::string> & op,
                                   vector<ST::string> & vn)
  {

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  */

  int f;

  f = op[4].strtodouble(lambdastart);
  f = op[5].strtodouble(a_invgamma);
  f = op[6].strtodouble(b_invgamma_orig);

  if (op[38] == "true")
    {
    lambdaconst = true;
    nosamples =true;
    }


  f = op[47].strtodouble(tildea);
  f = op[48].strtodouble(tildeb);

  if (op[49] == "true")
    {
    cauchy = true;
    }
  else
    cauchy = false;

  }


FC_nonp_variance::FC_nonp_variance(void)
  {

  }



FC_nonp_variance::FC_nonp_variance(MASTER_OBJ * mp, unsigned & enr, GENERAL_OPTIONS * o,
                 DISTR * lp, const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                 vector<ST::string> & vn)
     : FC(o,t,1,2,fp)
  {
  FCnonpp = FCn;
  likep = lp;
  designp = Dp;
  masterp = mp;
  equationnr = enr,
  lambdaconst = false;

  read_options(op,vn);

  datamatrix betanew(1,2);
  betanew(0,0) = likep->get_scale()/lambdastart;
  betanew(0,1) = lambdastart;
  setbeta(betanew);

  FCnonpp->tau2 = beta(0,0);
  FCnonpp->lambda = beta(0,1);

  }


FC_nonp_variance::FC_nonp_variance(const FC_nonp_variance & m)
    : FC(FC(m))
  {
  FCnonpp = m.FCnonpp;
  likep = m.likep;
  designp = m.designp;
  masterp = m.masterp;
  equationnr = m.equationnr;
  a_invgamma = m.a_invgamma;
  b_invgamma_orig = m.b_invgamma_orig;
  b_invgamma = m.b_invgamma;
  lambdastart = m.lambdastart;
  lambdaconst = m.lambdaconst;
  tildea = m.tildea;
  tildeb = m.tildeb;
  cauchy = m.cauchy;
  }


const FC_nonp_variance & FC_nonp_variance::operator=(const FC_nonp_variance & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  FCnonpp = m.FCnonpp;
  likep = m.likep;
  designp = m.designp;
  masterp = m.masterp;
  equationnr = m.equationnr;
  a_invgamma = m.a_invgamma;
  b_invgamma_orig = m.b_invgamma_orig;
  b_invgamma = m.b_invgamma;
  lambdastart = m.lambdastart;
  lambdaconst = m.lambdaconst;
  tildea = m.tildea;
  tildeb = m.tildeb;
  cauchy = m.cauchy;
  return *this;
  }


void FC_nonp_variance::update(void)
  {

  // TEST
  //  ofstream out("c:\\bayesx\\test\\results\\param.res");
  //  (FCnonpp->param).prettyPrint(out);
  // END: TEST

  b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

  if (lambdaconst == false)
    {

    if (cauchy == true)
      {
      double quadf = designp->penalty_compute_quadform(FCnonpp->param);
      double gamma = rand_invgamma(designp->rankK/2 +tildea ,0.5*quadf+tildeb);

      double u = log(uniform());

      double fcold = -(0.5*designp->rankK+0.5)*log(beta(0,0))-1/(2*beta(0,0))*quadf-log(1+beta(0,0));
      double fcnew = -(0.5*designp->rankK+0.5)*log(gamma)-1/(2*gamma)*quadf-log(1+gamma);

      if (u <= (fcnew - fcold ))
        {

        beta(0,0) = gamma;
        acceptance++;
        }

      }
    else
      {
      beta(0,0) = rand_invgamma(a_invgamma+0.5*designp->rankK,
                  b_invgamma+0.5*designp->penalty_compute_quadform(FCnonpp->param));
      acceptance++;
      }

    beta(0,1) = likep->get_scale()/beta(0,0);

    FCnonpp->tau2 = beta(0,0);


    FC::update();

    }

  }


bool FC_nonp_variance::posteriormode(void)
  {

  b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

  beta(0,0) = likep->get_scale()/beta(0,1);

  FCnonpp->tau2 = beta(0,0);

  posteriormode_betamean();

  return true;
  }



void FC_nonp_variance::outresults(ofstream & out_stata,ofstream & out_R,
                                  const ST::string & pathresults)
  {

  FC::outresults(out_stata,out_R,"");

  ST::string l1 = ST::doubletostring(optionsp->lower1,4);
  ST::string l2 = ST::doubletostring(optionsp->lower2,4);
  ST::string u1 = ST::doubletostring(optionsp->upper1,4);
  ST::string u2 = ST::doubletostring(optionsp->upper2,4);

  ST::string nl1 = ST::doubletostring(optionsp->lower1,4);
  ST::string nl2 = ST::doubletostring(optionsp->lower2,4);
  ST::string nu1 = ST::doubletostring(optionsp->upper1,4);
  ST::string nu2 = ST::doubletostring(optionsp->upper2,4);
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');

  ST::string vstr;


  if (optionsp->samplesize > 1)
    {

    FC::outresults_acceptance();

    vstr = "    Mean:         ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betamean(0,0),6) + "\n");

    vstr = "    Std. dev.:    ";

    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(sqrt(betavar(0,0)),6) + "\n");

    vstr = "    " + l1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_lower(0,0),6) + "\n");

    vstr = "    " + l2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_lower(0,0),6) + "\n");

    vstr = "    50% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu50(0,0),6) + "\n");

    vstr = "    " + u1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_upper(0,0),6) + "\n");

    vstr = "    " + u2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_upper(0,0),6) + "\n");

    optionsp->out("\n");

    }
  else
    {
    optionsp->out("    Smoothing parameter: " +
    ST::doubletostring(betamean(0,1),6) + "\n");

    optionsp->out("\n");
    }

//  out_R << "term=" << title <<  ";" << endl;

  if (pathresults.isvalidfile() != 1)
    {

    optionsp->out("    Results for the variance component are also stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    ST::string paths = pathresults.substr(0,pathresults.length()-4) +
                                 "_sample.raw";

    out_R << "pathvarsample=" << paths << endl;
//    out_R << "filetype=param; path=" << pathresults << ";" <<  endl;

    ofstream ou(pathresults.strtochar());

    if (optionsp->samplesize > 1)
      {
      ou << "pmean  pstd  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
      nu1 << "   pqu" << nu2 << endl;
      }
    else
      {
      ou << "pmean" << endl;
      }

    ou << betamean(0,0) << "  ";
    if (optionsp->samplesize > 1)
      {
      if (betavar(0,0) < 0.0000000000001)
        ou << 0 << "  ";
      else
        ou << sqrt(betavar(0,0)) << "  ";
      ou << betaqu_l1_lower(0,0) << "  ";
      ou << betaqu_l2_lower(0,0) << "  ";
      ou << betaqu50(0,0) << "  ";
      ou << betaqu_l2_upper(0,0) << "  ";
      ou << betaqu_l1_upper(0,0) << "  " << endl;
      }

    optionsp->out("\n");
    }

  }


void FC_nonp_variance::outoptions(void)
  {
  if (cauchy)
    {
    optionsp->out("  Cauchy prior\n");

    optionsp->out("  Hyperparameter tildea for proposal density: " +
                ST::doubletostring(tildea) + "\n" );

    optionsp->out("  Hyperparameter tildeb for proposal density: " +
                ST::doubletostring(tildeb) + "\n" );

    }
  else
    {

    optionsp->out("  Inverse gamma prior\n");

    b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

    optionsp->out("  Hyperprior a for variance parameter: " +
                ST::doubletostring(a_invgamma) + "\n" );
    optionsp->out("  Hyperprior b for variance parameter: " +
                ST::doubletostring(b_invgamma) + "\n" );
    optionsp->out("\n");
    }
  }


void FC_nonp_variance::reset(void)
  {

  datamatrix betanew(1,2);
  betanew(0,0) = likep->get_scale()/lambdastart;
  betanew(0,1) = lambdastart;
  setbeta(betanew);

  FCnonpp->tau2 = beta(0,0);
  FCnonpp->lambda = beta(0,1);
//  transform(0,0) = 1;
//  transform(1,0) = 1;

  }



//------------------------------------------------------------------------------
//--- CLASS: FC_nonp_variance_varselection implementation of member functions --
//------------------------------------------------------------------------------

void FC_nonp_variance_varselection::read_options(vector<ST::string> & op,
                                   vector<ST::string> & vn)
  {
  FC_nonp_variance::read_options(op,vn);

  int f;

  f = op[41].strtodouble(r);
  }


FC_nonp_variance_varselection::FC_nonp_variance_varselection(void)
  {

  }



FC_nonp_variance_varselection::FC_nonp_variance_varselection(MASTER_OBJ * mp,
                 unsigned & enr, GENERAL_OPTIONS * o,
                 DISTR * lp, const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                 vector<ST::string> & vn)
     : FC_nonp_variance(mp,enr,o,lp,t,fp,Dp,FCn,op,vn)
  {

  read_options(op,vn);

  FC_delta = FC(o,"",1,1,"");
  FC_delta.setbeta(1,1,1);

  FC_psi2 = FC(o,"",1,1,"");
  FC_psi2.setbeta(1,1,0.5);


  FC_omega = FC(o,"",1,1,"");

  FC_omega.setbeta(1,1,0.5);


  a_omega = 1;

  b_omega = 1;


  v = 5;

  Q = 25;


  }


FC_nonp_variance_varselection::FC_nonp_variance_varselection(const FC_nonp_variance_varselection & m)
    : FC_nonp_variance(FC_nonp_variance(m))
  {
  FC_delta = m.FC_delta;
  FC_psi2 = m.FC_psi2;
  FC_omega = m.FC_omega;
  a_omega = m.a_omega;
  b_omega = m.b_omega;
  v = m.v;
  Q = m.Q;
  r = m.r;
  X = m.X;
  }


const FC_nonp_variance_varselection & FC_nonp_variance_varselection::operator=(const FC_nonp_variance_varselection & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  FC_delta = m.FC_delta;
  FC_psi2 = m.FC_psi2;
  FC_omega = m.FC_omega;
  a_omega = m.a_omega;
  b_omega = m.b_omega;
  v = m.v;
  Q = m.Q;
  r = m.r;
  X = m.X;
  return *this;
  }


void FC_nonp_variance_varselection::update(void)
  {

  unsigned i;

  // updating psi2

  double r_delta;
  if (FC_delta.beta(0,0) == 0)
    r_delta = r;
  else
    r_delta = 1;

  FC_psi2.beta(0,0) = rand_invgamma(v+0.5,Q+0.5*beta(0,0)*r_delta);

  FC_psi2.update();

  // end: updating psi2

  // updating delta

  double u = uniform();
  double L = 1/sqrt(r)*exp(- beta(0,0)/(2*FC_psi2.beta(0,0))*(1/r-1));
  double pr1 = 1/(1+ ((1-FC_omega.beta(0,0))/FC_omega.beta(0,0))*L);
  if (u <=pr1)
    {
    FC_delta.beta(0,0) = 1;
    r_delta = 1;
    }
  else
    {
    FC_delta.beta(0,0) = 0;
    r_delta = r;
    }

  FC_delta.update();

  // end: updating delta


  // updating w

  FC_omega.beta(0,0) = randnumbers::rand_beta(a_omega+FC_delta.beta(0,0),
                                          b_omega+1-FC_delta.beta(0,0));

  FC_omega.update();

  // end: updating w


  // updating tau2

  FCnonpp->designp->compute_effect(X,FCnonpp->beta);

  double * worklin;
  if (likep->linpred_current==1)
    worklin = likep->linearpred1.getV();
  else
    worklin = likep->linearpred2.getV();

  double Sigmatau;
  double mutau = 0;
  double * Xp = X.getV();
  double * responsep = likep->workingresponse.getV();
  double varinv = 1/(likep->get_scale()*beta(0,0));
  double xtx=0;
  for (i=0;i<X.rows();i++,Xp++,responsep++,worklin++)
    {
    xtx += pow(*Xp,2);
    mutau += (*Xp) * ((*responsep) - (*worklin)+(*Xp));
    }

  Sigmatau = 1/(varinv*xtx + 1/(r_delta*FC_psi2.beta(0,0)));

  mutau *= Sigmatau/(likep->get_scale()*sqrt(beta(0,0)));

  double tau = mutau + sqrt(Sigmatau) * rand_normal();

  double tau2 = tau*tau;
  if (tau2 < 0.000000001)
    tau2 = 0.000000001;

  beta(0,0) = tau2;

  beta(0,1) = likep->get_scale()/beta(0,0);

  FCnonpp->tau2 = beta(0,0);

  // end: updating tau2

  acceptance++;
  FC::update();

  }


bool FC_nonp_variance_varselection::posteriormode(void)
  {
  bool t = FC_nonp_variance::posteriormode();

  FC_psi2.beta(0,0) = beta(0,0);

  return true;
  }



void FC_nonp_variance_varselection::outresults(ofstream & out_stata,ofstream & out_R,
                                  const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    ST::string pathresults_delta = pathresults.substr(0,pathresults.length()-4) + "_delta.res";
    ST::string pathresults_omega = pathresults.substr(0,pathresults.length()-4) + "_omega.res";

    FC_nonp_variance::outresults(out_stata,out_R,pathresults);

    FC_delta.outresults(out_stata,out_R,"");
    FC_omega.outresults(out_stata,out_R,pathresults_omega);


    optionsp->out("    Inclusion probability: " + ST::doubletostring(FC_delta.betamean(0,0),6)  + "\n");
    optionsp->out("\n");
    optionsp->out("    Results for the inclusion probabilities are also stored in file\n");
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

    ou << "pmean" << endl;
    ou << FC_delta.betamean(0,0) << endl;
    }


//  FC_nonp_variance::outresults(out_stata,out_R,pathresults);

  }


void FC_nonp_variance_varselection::get_samples(
   const ST::string & filename,ofstream & outg) const
  {
  FC_nonp_variance::get_samples(filename,outg);

  ST::string filename_delta = filename.substr(0,filename.length()-4) + "_delta.raw";
  FC_delta.get_samples(filename_delta,outg);

/*
  ST::string filename_omega = filename.substr(0,filename.length()-4) + "_omega.raw";
  FC_omega.get_samples(filename_omega,outg);

  ST::string filename_Q = filename.substr(0,filename.length()-4) + "_Q.raw";
  FC_Q.get_samples(filename_Q,outg);
*/
  }


void FC_nonp_variance_varselection::outoptions(void)
  {
  FC_nonp_variance::outoptions();
  }


void FC_nonp_variance_varselection::reset(void)
  {

  FC_nonp_variance::reset();

  }


} // end: namespace MCMC



