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



#include "FC_predictive_check.h"
#include "clstring.h"

//------------------------------------------------------------------------------
//---------------- CLASS: FC_predict implementation of member functions --------
//------------------------------------------------------------------------------


namespace MCMC
{



void FC_predictive_check::read_options(vector<ST::string> & op,
vector<ST::string> & vn)
  {

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
  13      samplemult
  14      constraints
  */

  }


FC_predictive_check::FC_predictive_check(void)
  {
  }


FC_predictive_check::FC_predictive_check(GENERAL_OPTIONS * o,DISTR * lp,
                 const ST::string & t, const ST::string & fp,datamatrix & dm,
                 vector<ST::string> & dn)
  : FC(o,t,1,1,fp)
  {

  likep = lp;

  sampled_responses = datamatrix(likep->nrobs,o->compute_samplesize(),0);

  designmatrix= dm;
  varnames = dn;

  nosamples=true;

  }


FC_predictive_check::FC_predictive_check(const FC_predictive_check & m)
  : FC(FC(m))
  {
  sampled_responses = m.sampled_responses;
  likep = m.likep;
  designmatrix = m.designmatrix;
  varnames = m.varnames;
  }


const FC_predictive_check & FC_predictive_check::operator=(
const FC_predictive_check & m)
  {
  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  sampled_responses = m.sampled_responses;
  likep = m.likep;
  designmatrix = m.designmatrix;
  varnames = m.varnames;
  return *this;
  }


void  FC_predictive_check::update(void)
  {

  if(
     (optionsp->nriter > optionsp->burnin)
     &&
     ((optionsp->nriter-optionsp->burnin-1) % (optionsp->step) == 0)
    )
    {

    unsigned samplesize = optionsp->samplesize;

    likep->sample_responses(samplesize-1,sampled_responses); 

    }

  }



bool FC_predictive_check::posteriormode(void)
  {

  return true;
  }


void FC_predictive_check::outoptions(void)
  {

  }


void FC_predictive_check::outresults(ofstream & out_stata, ofstream & out_R,
                            const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    optionsp->out("  POSTERIOR PREDICTIVE CHECKS: \n",true);
    optionsp->out("\n");

    optionsp->out("    Samples for posterior predictive checks are stored in\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    likep->outresults_predictive_check(designmatrix,sampled_responses);

    ofstream outres(pathresults.strtochar());

    unsigned nrobs = designmatrix.rows();
    unsigned i,j;

    for(j=0;j<varnames.size();j++)
      outres << varnames[j] << "  ";
    for(j=0;j<sampled_responses.cols();j++)
      outres << "s" << (j+1) << "  ";
    outres << endl;

    for (i=0;i<nrobs;i++)
      {

      for(j=0;j<designmatrix.cols();j++)
        outres << designmatrix(i,j) << "  ";

      for(j=0;j<sampled_responses.cols();j++)
        outres << sampled_responses(i,j) << "  ";

      outres << endl;
      }

    }   // end if (pathresults.isvalidfile() != 1)

  }


void FC_predictive_check::reset(void)
  {

  }



} // end: namespace MCMC



