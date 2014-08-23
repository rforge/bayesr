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



#include "FC_predict_predictor.h"
#include "clstring.h"

//------------------------------------------------------------------------------
//---------------- CLASS: FC_predict implementation of member functions --------
//------------------------------------------------------------------------------


namespace MCMC
{

void FC_predict_predictor::compute_autocorr_all(const ST::string & path, unsigned lag,
                                    ofstream & outg) const
  {

  }


void FC_predict_predictor::read_options(vector<ST::string> & op,vector<ST::string> & vn)
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

FC_predict_predictor::FC_predict_predictor(void)
  {
  }


FC_predict_predictor::FC_predict_predictor(GENERAL_OPTIONS * o,DISTR * lp,
     const ST::string & t,
     const ST::string & fp, const ST::string & fpd, datamatrix & dm,
      vector<ST::string> & dn)
  : FC(o,t,1,1,fp)
  {
  nosamples = true;
  likep = lp;
  designmatrix= dm;
  varnames = dn;
  setbeta(lp->nrobs,1,0);

  }


FC_predict_predictor::FC_predict_predictor(const FC_predict_predictor & m)
  : FC(FC(m))
  {
  likep = m.likep;
  designmatrix = m.designmatrix;
  varnames = m.varnames;
  }


const FC_predict_predictor &
FC_predict_predictor::operator=(const FC_predict_predictor & m)
  {
  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  likep = m.likep;
  designmatrix = m.designmatrix;
  varnames = m.varnames;
  return *this;
  }


void  FC_predict_predictor::update(void)
  {

  likep->FCpredict_betamean = &betamean;

  get_predictor();

  acceptance++;

  FC::update();

  }


void FC_predict_predictor::get_predictor(void)
  {

  unsigned i;
  double * betap = beta.getV();

  double * worklinp;
  if (likep->linpred_current==1)
    worklinp = likep->linearpred1.getV();
  else
    worklinp = likep->linearpred2.getV();

  for(i=0;i<likep->nrobs;i++,worklinp++,betap++)
    {
    *betap = *worklinp;
    }
  }


bool FC_predict_predictor::posteriormode(void)
  {

  get_predictor();

  posteriormode_betamean();

  return true;
  }


void FC_predict_predictor::outresults(ofstream & out_stata, ofstream & out_R,
                            const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    FC::outresults(out_stata,out_R,"");

    optionsp->out("  PREDICTED VALUES: \n",true);
    optionsp->out("\n");

    optionsp->out("    Results for the predictor are stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");
    out_R << "predict=" << pathresults << ";" <<  endl;

    ofstream outres(pathresults.strtochar());

    optionsp->out("\n");

    unsigned i,j;

    ST::string l1 = ST::doubletostring(optionsp->lower1,4);
    ST::string l2 = ST::doubletostring(optionsp->lower2,4);
    ST::string u1 = ST::doubletostring(optionsp->upper1,4);
    ST::string u2 = ST::doubletostring(optionsp->upper2,4);
    l1 = l1.replaceallsigns('.','p');
    l2 = l2.replaceallsigns('.','p');
    u1 = u1.replaceallsigns('.','p');
    u2 = u2.replaceallsigns('.','p');

    outres << "intnr" << "   ";

    for (i=1;i<varnames.size();i++)
      outres << varnames[i] << "   ";

    outres << "pmean_pred   ";

    if (optionsp->samplesize > 1)
      {
      outres << "pqu"  << l1  << "_pred   ";
      outres << "pqu"  << l2  << "_pred   ";
      outres << "pmed_pred   ";
      outres << "pqu"  << u1  << "_pred   ";
      outres << "pqu"  << u2  << "_pred   ";
      }

    outres << endl;

    double * workmean = betamean.getV();
    double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
    double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
    double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
    double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
    double * workbetaqu50 = betaqu50.getV();

   for(i=0;i<designmatrix.rows();i++,
          workmean++,
          workbetaqu_l1_lower_p++,
          workbetaqu_l2_lower_p++,
          workbetaqu50++,
          workbetaqu_l1_upper_p++,
          workbetaqu_l2_upper_p++)
      {

      outres << (i+1) << "   ";

      for (j=1;j<designmatrix.cols();j++)
        outres << designmatrix(i,j) << "   ";

      outres << *workmean << "   ";

      if (optionsp->samplesize > 1)
        {
        outres << *workbetaqu_l1_lower_p << "   ";
        outres << *workbetaqu_l2_lower_p << "   ";
        outres << *workbetaqu50 << "   ";
        outres << *workbetaqu_l2_upper_p << "   ";
        outres << *workbetaqu_l1_upper_p << "   ";
        }

      outres << endl;

     }

    }   // end if (pathresults.isvalidfile() != 1)

  }


void FC_predict_predictor::reset(void)
  {

  }



} // end: namespace MCMC



