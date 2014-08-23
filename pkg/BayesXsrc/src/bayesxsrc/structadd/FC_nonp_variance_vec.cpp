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



#include "FC_nonp_variance_vec.h"


//------------------------------------------------------------------------------
//------ CLASS: FC_non_variance_vec implementation of member functions ---------
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_nonp_variance_vec::read_options(vector<ST::string> & op,
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

  FC_nonp_variance::read_options(op,vn);

  }


FC_nonp_variance_vec::FC_nonp_variance_vec(void)
  {

  }



FC_nonp_variance_vec::FC_nonp_variance_vec(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,
                 DISTR * lp, const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,FC_nonp * FCn,vector<ST::string> & op,
                 vector<ST::string> & vn)
     : FC_nonp_variance(mp,enr, o,lp,t,fp,Dp,FCn,op,vn)
  {


  read_options(op,vn);

  datamatrix betanew(FCn->beta.rows(),1,likep->get_scale()/lambdastart);
  setbeta(betanew);

  FCnonpp->tau2 = 1;
  FCnonpp->lambda = likep->get_scale();
  Dp->compute_penalty2(beta);


  }


FC_nonp_variance_vec::FC_nonp_variance_vec(const FC_nonp_variance_vec & m)
    : FC_nonp_variance(FC_nonp_variance(m))
  {
  }


const FC_nonp_variance_vec & FC_nonp_variance_vec::operator=(
        const FC_nonp_variance_vec & m)
  {

  if (this==&m)
	 return *this;
  FC_nonp_variance::operator=(FC_nonp_variance(m));
  return *this;
  }


void FC_nonp_variance_vec::update(void)
  {

  }


bool FC_nonp_variance_vec::posteriormode(void)
  {


  b_invgamma = masterp->level1_likep[equationnr]->trmult*b_invgamma_orig;

  unsigned i;
  for(i=0;i<beta.rows();i++)
    beta(i,0) = likep->get_scale()/lambdastart;

  FCnonpp->lambda = likep->get_scale();
  FCnonpp->tau2 = 1;
  designp->compute_penalty2(beta);

  posteriormode_betamean();


  return true;
  }



void FC_nonp_variance_vec::outresults(ofstream & out_stata,ofstream & out_R,
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

  if (pathresults.isvalidfile() != 1)
    {

    optionsp->out("    Results for the variance component are also stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    ofstream ou(pathresults.strtochar());

    if (optionsp->samplesize > 1)
      {
      ou << "intnr  pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
      nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
      }
    else
      {
      ou << "intnr  pmean" << endl;
      }


    unsigned i;
    for (i=0;i<beta.rows();i++)
      {
      ou << (i+1) << "  "  << betamean(i,0) << "  ";
      if (optionsp->samplesize > 1)
        {
        ou << (betavar(i,0)<0.0?0.0:sqrt(betavar(i,0))) << "  ";
        ou << betaqu_l1_lower(i,0) << "  ";
        ou << betaqu_l2_lower(i,0) << "  ";
        ou << betaqu50(i,0) << "  ";
        ou << betaqu_l2_upper(i,0) << "  ";
        ou << betaqu_l1_upper(i,0) << "  ";
        ou << betamin(i,0) << "  ";
        ou << betamax(i,0) << "  " << endl;
        }
      }

    optionsp->out("\n");
    }


  }


void FC_nonp_variance_vec::outoptions(void)
  {

  FC_nonp_variance::outoptions();

  }


void FC_nonp_variance_vec::reset(void)
  {
  }


} // end: namespace MCMC



