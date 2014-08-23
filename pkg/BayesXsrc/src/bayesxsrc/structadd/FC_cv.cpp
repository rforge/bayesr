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



#include "FC_cv.h"
#include "clstring.h"

//------------------------------------------------------------------------------
//----------------- CLASS: FC_cv implementation of member functions ------------
//------------------------------------------------------------------------------


namespace MCMC
{



void FC_cv::read_options(vector<ST::string> & op,
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



void FC_cv::get_ind(void)
  {

  vector<FC_hrandom>::iterator it = hrandoms->begin();
  ind = (*it).designp->ind;
  nrcat = (*it).beta.rows();
  effectvalues = (*it).designp->effectvalues;

  }



FC_cv::FC_cv(void)
  {
  }


FC_cv::FC_cv(GENERAL_OPTIONS * o,DISTR * lp,
                 const ST::string & t, const ST::string & fp,
                 vector<FC_hrandom> * FChs)
  : FC(o,t,lp->nrobs,1,fp)
  {

  likep = lp;

  hrandoms = FChs;
  sampled_etas = datamatrix(likep->nrobs,o->compute_samplesize(),0);
  sampled_responses = datamatrix(likep->nrobs,o->compute_samplesize(),0);
  sampled_likelihood = datamatrix(likep->nrobs,o->compute_samplesize(),0);

  FC_sampled_l = FC(o,"",likep->nrobs,1,fp+"_like");

  effect = datamatrix(likep->nrobs,1,0);
  linpred = datamatrix(likep->nrobs,1,0);
  size =   hrandoms->size();
  get_ind();
  }


FC_cv::FC_cv(const FC_cv & m)
  : FC(FC(m))
  {
  ind = m.ind;
  nrcat = m.nrcat;
  sampled_etas = m.sampled_etas;
  sampled_responses = m.sampled_responses;
  sampled_likelihood = m.sampled_likelihood;
  FC_sampled_l = m.FC_sampled_l;
  likep = m.likep;
  hrandoms = m.hrandoms;
  effect = m.effect;
  linpred = m.linpred;
  size = m.size;
  e_score = m.e_score;
  log_score = m.log_score;
  effectvalues = m.effectvalues;
  }


const FC_cv & FC_cv::operator=(
const FC_cv & m)
  {
  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  ind = m.ind;
  nrcat = m.nrcat;
  sampled_etas = m.sampled_etas;
  sampled_responses = m.sampled_responses;
  sampled_likelihood = m.sampled_likelihood;
  FC_sampled_l = m.FC_sampled_l;
  likep = m.likep;
  hrandoms = m.hrandoms;
  effect = m.effect;
  linpred = m.linpred;
  size = m.size;
  e_score = m.e_score;
  log_score = m.log_score;
  effectvalues = m.effectvalues;
  return *this;
  }


void  FC_cv::update(void)
  {

  if(
     (optionsp->nriter > optionsp->burnin)
     &&
     ((optionsp->nriter-optionsp->burnin-1) % (optionsp->step) == 0)
    )
    {

    unsigned samplesize = optionsp->samplesize;

    unsigned i;

    vector<FC_hrandom>::iterator it = hrandoms->begin();

    if (likep->linpred_current==1)
      linpred.assign(likep->linearpred1);
    else
      linpred.assign(likep->linearpred2);


    for(i=0;i<size;i++,++it)
      {
      (*it).compute_effect_cv(effect);

      linpred.minus(linpred,effect);

      (*it).sample_for_cv(effect);

      linpred.plus(linpred,effect);
      }

    double * workse = sampled_etas.getV()+samplesize-1;
    double * workli = sampled_likelihood.getV()+samplesize-1;
    double * worklin = linpred.getV();
    double * wr = likep->response.getV();
    double * ww = likep->weight.getV();

    double muhelp;
//    double scalehelp=likep->get_scale(true);
    double scalehelp=likep->get_scale();
    for (i=0;i<sampled_etas.rows();i++,workse+=sampled_etas.cols(),
    workli+=sampled_likelihood.cols(),worklin++,wr++,ww++)
      {
      // *worklin *= likep->trmult;
      *workse = *worklin;

//      likep->compute_mu(worklin,&muhelp,true);
      likep->compute_mu(worklin,&muhelp);

      likep->compute_deviance(wr,ww,&muhelp, workli,&scalehelp);

      *workli *= -0.5;
      }

    likep->sample_responses_cv(samplesize-1,linpred,sampled_responses);

    double * workbeta = beta.getV();
    double * sr = sampled_responses.getV()+samplesize-1;

    for(i=0;i<beta.rows();i++,workbeta++,sr+=sampled_responses.cols())
      *workbeta = *sr;


    workbeta = FC_sampled_l.beta.getV();
    double * sl = sampled_likelihood.getV()+samplesize-1;

    for(i=0;i<FC_sampled_l.beta.rows();i++,workbeta++,sl+=sampled_likelihood.cols())
      *workbeta = *sl;

    }


  FC::update();
  FC_sampled_l.update();

  }



bool FC_cv::posteriormode(void)
  {

  return true;
  }


void FC_cv::outoptions(void)
  {

  }


void FC_cv::outresults(ofstream & out_stata, ofstream & out_R,
                            const ST::string & pathresults)
  {


  if (pathresults.isvalidfile() != 1)
    {

    FC::outresults(out_stata,out_R,pathresults);

    optionsp->out("  Marshall-Spiegelhalter Cross Validation: \n",true);
    optionsp->out("\n");

    optionsp->out("    Estimated individual observation samples are stored in\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");


    ST::string pathresults_like = pathresults.substr(0,pathresults.length()-4)+
                                  "_like.res";

    FC_sampled_l.outresults(out_stata,out_R,pathresults_like);

    optionsp->out("    Estimated individual observation likelihoods are stored in\n");
    optionsp->out("    " +  pathresults_like + "\n");
    optionsp->out("\n");



//    unsigned nrobs = sampled_etas.rows();
    unsigned i;

    /*
    ofstream outres(pathresults.strtochar());

    for(j=0;j<sampled_etas.cols();j++)
      outres << "s_eta_" << (j+1) << "  ";

    for(j=0;j<sampled_etas.cols();j++)
      outres << "s_resp_" << (j+1) << "  ";

    outres << endl;

    for (i=0;i<nrobs;i++)
      {

      for(j=0;j<sampled_etas.cols();j++)
        outres << sampled_etas(i,j) << "  ";

      for(j=0;j<sampled_responses.cols();j++)
        outres << sampled_responses(i,j) << "  ";

      outres << endl;

      }
     */

    // Energy score

    double es = compute_energyscore();


    ST::string pathresults_e = pathresults.substr(0,pathresults.length()-4)+
                               "_energy.res";


    ofstream out2(pathresults_e.strtochar());
    out2 << "id   score" << endl;

    for (i=0;i<e_score.rows();i++)
      out2 << effectvalues[i] << "  " << e_score(i,0) << endl;

    // Log-score


    double ls = compute_logscore();

    ST::string pathresults_l = pathresults.substr(0,pathresults.length()-4)+
                               "_logscore.res";

    ofstream out3(pathresults_l.strtochar());
    out3 << "id   score" << endl;

    for (i=0;i<e_score.rows();i++)
      out3 << effectvalues[i] << "  " << log_score(i,0) << endl;

    optionsp->out("    Estimated energy scores are stored in\n");
    optionsp->out("    " +  pathresults_e + "\n");
    optionsp->out("\n");

    optionsp->out("    Estimated log-scores are stored in\n");
    optionsp->out("    " +  pathresults_l + "\n");
    optionsp->out("\n");

    optionsp->out("    Mean energy score: " + ST::doubletostring(es,8) + "\n");
    optionsp->out("    Mean log score: " + ST::doubletostring(ls,8) + "\n");


    }   // end if (pathresults.isvalidfile() != 1)


  }


double FC_cv::compute_energyscore(void)
  {


  unsigned s,i,j;

  unsigned S = sampled_responses.cols();
  unsigned I = sampled_responses.rows();

  //double * srp = sampled_responses.getV();
  // double * esp1;
  // double * esp2;


  datamatrix es1 = datamatrix(nrcat,S,0);
  datamatrix es2 = datamatrix(nrcat,S,0);

/*
  for (i=0;i<I;i++)
    {
    esp1 = es1.getV()+ind(i,0)*S;
    esp2 = es2.getV()+ind(i,0)*S;

    for (s=0;s<S-1;s++,srp++,esp1++,esp2++)
      {
      *esp1 += pow(sampled_responses(i,s)-likep->response_untransformed(i,0),2);
      *esp2 += pow(sampled_responses(i,s+1)-sampled_responses(i,s),2);
      }
    esp1++;
    *esp1 += pow(sampled_responses(i,s)-likep->response_untransformed(i,0),2);

    }
*/


  unsigned in;

  for (i=0;i<I;i++)
    {

    in = ind(i,0);

    for (s=0;s<S-1;s++)
      {
      es1(in,s) += pow(sampled_responses(i,s)-likep->response(i,0),2);
      es2(in,s) += pow(sampled_responses(i,s+1)-sampled_responses(i,s),2);
      }

    es1(in,S-1) += pow(sampled_responses(i,s)-likep->response(i,0),2);

    }

//  ofstream out("c:\\bayesx\\testh\\results\\es1.res");
//  es1.prettyPrint(out);



  if (e_score.rows() != nrcat)
    e_score = datamatrix(nrcat,1,0);

  double h;
  for (j=0;j<nrcat;j++)
    {
    for (s=0;s<S;s++)
      e_score(j,0) += sqrt(es1(j,s));

    e_score(j,0) /= S;

    h=0;
    for(s=0;s<S-1;s++)
      h += sqrt(es2(j,s));

    e_score(j,0) -= h/(2*(S-1));
    }


  return e_score.mean(0);
  }



double FC_cv::compute_logscore(void)
  {

  unsigned s,i;

//  ofstream out("c:\\bayesx\\testh\\results\\sampled_likelihood.res");
//  sampled_likelihood.prettyPrint(out);

  unsigned S = sampled_likelihood.cols();       // S = number of samples
  unsigned I = sampled_likelihood.rows();       // I = number of observations
                                                // i-th row, s-th column contains
                                                // log-likelihood of i-th
                                                // observation in sample s

  if (log_score.rows() != nrcat)
    log_score = datamatrix(nrcat,1,0);

  datamatrix lhelp(nrcat,S,0);

  unsigned in;

  for (i=0;i<I;i++)
    {

    in = ind(i,0);

    for (s=0;s<S;s++)
      {
      lhelp(in,s) += sampled_likelihood(i,s);
      }



//    log_score(in,0) = log(S) - log(log_score(in,0));
    }


  for (i=0;i<nrcat;i++)
    {
    log_score(i,0) = 0;
    for (s=0;s<S;s++)
      log_score(i,0) +=  exp(lhelp(i,s));

    log_score(i,0) = log(static_cast<double>(S)) - log(log_score(i,0));
    }

  return log_score.mean(0);

  }


void FC_cv::reset(void)
  {

  }



} // end: namespace MCMC



