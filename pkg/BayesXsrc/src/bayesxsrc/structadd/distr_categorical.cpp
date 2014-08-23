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


#include "distr_categorical.h"
//#include "gsl/gsl_randist.h"
//#include "gsl/gsl_cdf.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//------------------ CLASS: DISTRIBUTION_logit_fruehwirth ----------------------
//------------------------------------------------------------------------------

DISTR_logit_fruehwirth::DISTR_logit_fruehwirth(const int h, GENERAL_OPTIONS * o,
                                             const datamatrix r,
                                             const datamatrix & w)
  : DISTR_binomial(o, r, w), H(h)
  {

  predictor_name = "pi";
  outexpectation = true;

  family = "Binomial_l1";
  updateIWLS = false;


  SQ = datamatrix(6,5,0);
  SQ(0,1) = 1/1.2131;
  SQ(1,1) = 1/2.9955;
  SQ(2,1) = 1/7.5458;

  SQ(0,4) = 1/0.68159;
  SQ(1,4) = 1/1.2419;
  SQ(2,4) = 1/2.2388;
  SQ(3,4) = 1/4.0724;
  SQ(4,4) = 1/7.4371;
  SQ(5,4) = 1/13.772;


  weights_mixed = datamatrix(6,5,0);
  weights_mixed(0,1) = 0.2522;
  weights_mixed(1,1) = 0.58523;
  weights_mixed(2,1) = 0.16257;


  weights_mixed(0,4) = 0.018446;
  weights_mixed(1,4) = 0.17268;
  weights_mixed(2,4) = 0.37393;
  weights_mixed(3,4) = 0.31697;
  weights_mixed(4,4) = 0.1089;
  weights_mixed(5,4) = 0.0090745;

  linpredminlimit=-10;
  linpredmaxlimit= 10;

  check_errors();

  }



void DISTR_logit_fruehwirth::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight == 1)
        {
        if (*workresp != 0 && *workresp!=1)
          {
          errors=true;
          errormessages.push_back("ERROR: response must be either 0 or 1\n");
          }

        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: weighted regression not allowed\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


DISTR_logit_fruehwirth::DISTR_logit_fruehwirth(const DISTR_logit_fruehwirth & nd)
	: DISTR_binomial(DISTR_binomial(nd)) , H(nd.H), SQ(nd.SQ), weights_mixed(nd.weights_mixed)
	{
	}

const DISTR_logit_fruehwirth & DISTR_logit_fruehwirth::operator=(const DISTR_logit_fruehwirth & nd)
	{
	if (this==&nd)
  	return *this;
  DISTR_binomial::operator=(DISTR_binomial(nd));
  H = nd.H;
  SQ = nd.SQ;
  weights_mixed = nd.weights_mixed;
  return *this;
	}

/*
double DISTR_logit_fruehwirth::compute_MSE()
  {

  }

void DISTR_logit_fruehwirth::compute_mu()
  {
  // datamatrix test(3,3,0);
  // ofstream out("c:\\bayesx\\temp\\r.res");
  // test.prettyPrint(out);
  }
*/




void DISTR_logit_fruehwirth::outoptions()
  {
  DISTR::outoptions();
  optionsp->out("  Response function: logistic distribution function\n");
  optionsp->out("  Number of mixture components: " + ST::inttostring(H) + "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_logit_fruehwirth::update(void)
  {

  double * workresp;
  double * workwresp;
  double * weightwork;
  double * wweightwork;

  workresp = response.getV();
  workwresp = workingresponse.getV();
  weightwork = weight.getV();
  wweightwork = workingweight.getV();

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double lambda;
  double U;
  datamatrix weights_aux(H,1);


  for(unsigned i=0;i<nrobs;i++,worklin++,workresp++,weightwork++,workwresp++,wweightwork++)
    {
    lambda= exp(*worklin);
    U = uniform();
    *workwresp = log(lambda*U+*workresp)-log(1-U+lambda*(1-*workresp));

    //weights_mixed
    for(int j=0; j < H; j++)
    	{
      weights_aux(j,0) = weights_mixed(j,H-2)*sqrt(SQ(j,H-2)) * exp(-1/2*  pow((*workresp - *worklin), 2)*SQ(j,H-2) );
      }

    //distribution function
    for(int j=1; j <H; j++)
    	{
       weights_aux(j,0) = weights_aux(j-1,0) + weights_aux(j,0);
      }

    U = uniform();
    U = U*weights_aux(H-1,0);	//scale to [0, max]


    int iaux = 0;
    while (U > weights_aux(iaux,0))
      {
      iaux++;
      }

    *wweightwork =  SQ(iaux,H-2);

    }
  }


bool DISTR_logit_fruehwirth::posteriormode(void)
  {
  return DISTR_binomial::posteriormode();
  }

//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_binomial --------------------------
//------------------------------------------------------------------------------


DISTR_binomial::DISTR_binomial(GENERAL_OPTIONS * o, const datamatrix & r,
                               const datamatrix & w)
  : DISTR(o,r,w)

  {

  predictor_name = "pi";
  outexpectation = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  family = "Binomial";
  updateIWLS = true;

  linpredminlimit=-10;
  linpredmaxlimit= 10;

  check_errors();
  }


void DISTR_binomial::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {
        if (*workresp != int(*workresp))
          {
          errors=true;
          errormessages.push_back("ERROR: response cannot be binomial; values must be integer numbers\n");
          }

        if (*workresp < 0)
          {
          errors=true;
          errormessages.push_back("ERROR: response cannot be binomial; some values are negative\n");
          }

        if (*workresp > *workweight)
          {
          errors = true;
          errormessages.push_back("ERROR: response cannot be binomial;\n");
          errormessages.push_back("       number of successes larger than number of trials for some values\n");
          }

        *workresp = *workresp/(*workweight);
        }
      else if (*workweight < 0)
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


const DISTR_binomial & DISTR_binomial::operator=(
                                      const DISTR_binomial & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  return *this;
  }


DISTR_binomial::DISTR_binomial(const DISTR_binomial & nd)
   : DISTR(DISTR(nd))
  {
  }


void DISTR_binomial::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: logistic distribution function\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_binomial::get_intercept_start(void)
  {
  double m = response.mean(0);
  return log(m/(1-m));
  }


double DISTR_binomial::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {

  if (*linpred >= 10)
    return *weight *(*response * *linpred - *linpred);
  else
    return *weight *(*response * *linpred - log(1+exp(*linpred)));
  }


double DISTR_binomial::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {

  if (*linpred >= 10)
    return *response * (*linpred) - *linpred;
  else
    return *response * (*linpred) - log(1+exp(*linpred));
  }


void DISTR_binomial::compute_mu(const double * linpred,double * mu)
  {
  double el = exp(*linpred);
  *mu = el/(1+el);
  }


void DISTR_binomial::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  double el = exp(*linpred[predstart_mumult]);
  *mu = el/(1+el);
  }



void DISTR_binomial::compute_deviance(const double * response,
                   const double * weight,const double * mu,double * deviance,
                   double * scale) const
  {

  if (*response==0)
    *deviance = -2* *weight * log(1-*mu);
  else if (*response == 1)
    *deviance = -2* *weight*log(*mu);
  else
    *deviance = -2* *weight*( *response*log(*mu)+(1-*response)*log(1-*mu) );

  }


double DISTR_binomial::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  double el = exp(*linpred);
  double mu = el/(1+el);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;
  double v = mu*(1-mu);

  *workingweight = *weight * v;

  *workingresponse = *linpred + (*response - mu)/v;

  if (like)
    {
    if (*linpred >= 10)
      return *weight *(*response * *linpred - *linpred);
    else
      return *weight *(*response * *linpred - log(1+el));
    }
  else
    {
    return 0;
    }

  }


void DISTR_binomial::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {


  double el = exp(*linpred);
  double mu = el/(1+el);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;
  double v = mu*(1-mu);

  *workingweight =  v;

  *workingresponse = *linpred + (*response - mu)/v;

  if (compute_like)
    {
    if (*linpred >= 10)
      like += *response * (*linpred) - *linpred;
    else
      like += *response * (*linpred) - log(1+el);
    }

  }


void DISTR_binomial::compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  }

void DISTR_binomial::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {


  }


void DISTR_binomial::sample_responses(unsigned i,datamatrix & sr)
  {
  double * linpredp;

  if (linpred_current==1)
    linpredp = linearpred1.getV();
  else
    linpredp = linearpred2.getV();

  double * rp = sr.getV()+i;
  double mu;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
    {
    compute_mu(linpredp,&mu);

    *rp = randnumbers::rand_binom(1,mu);

    }

  }


void DISTR_binomial::sample_responses_cv(unsigned i,datamatrix & linpred,
                                         datamatrix & sr)
  {

  double * linpredp;

  linpredp = linpred.getV();

  double * rp = sr.getV()+i;
  double mu;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
    {
    compute_mu(linpredp,&mu);

    *rp = randnumbers::rand_binom(1,mu);
    }

  }



//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_binomialprobit --------------------
//------------------------------------------------------------------------------


DISTR_binomialprobit::DISTR_binomialprobit(GENERAL_OPTIONS * o,
                                           const datamatrix & r,const bool ut,
                                           const ST::string & ps,
                                           const datamatrix & w)
  : DISTR(o,r,w)

  {

  outexpectation = true;
  predictor_name = "pi";

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  family = "Binomial";
  updateIWLS = false;

  utilities = ut;
  if (utilities)
    FC_latentutilities = FC(o,"",nrobs,1,ps);

  linpredminlimit=-10;
  linpredmaxlimit= 10;

  }


void DISTR_binomialprobit::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight == 1)
        {
        if (*workresp != 0 && *workresp!=1)
          {
          errors=true;
          errormessages.push_back("ERROR: response must be either 0 or 1\n");
          }

        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: weighted regression not allowed\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }


const DISTR_binomialprobit & DISTR_binomialprobit::operator=(
                                      const DISTR_binomialprobit & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  utilities = nd.utilities;
  FC_latentutilities = nd.FC_latentutilities;
  return *this;
  }


DISTR_binomialprobit::DISTR_binomialprobit(const DISTR_binomialprobit & nd)
   : DISTR(DISTR(nd))
  {
  utilities = nd.utilities;
  FC_latentutilities = nd.FC_latentutilities;
  }


void DISTR_binomialprobit::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: standard normal (probit link)\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_binomialprobit::get_intercept_start(void)
  {
//  double m = response.mean(0);
//  return invPhi2(m);
  return 0;
  }


double DISTR_binomialprobit::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {

  if (*weight!=0)
    {
    double mu = randnumbers::Phi2(*linpred);
    if (*response > 0)
      return log(mu);
    else
      return log(1-mu);
    }
  else
    return 0;

  }



double DISTR_binomialprobit::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {

  double mu = randnumbers::Phi2(*linpred);
  if (*response > 0)
    return log(mu);
  else
    return log(1-mu);

  }


void DISTR_binomialprobit::compute_mu(const double * linpred,double * mu)
  {
  *mu = randnumbers::Phi2(*linpred);
  }


void DISTR_binomialprobit::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {
  *mu = randnumbers::Phi2(*linpred[predstart_mumult]);
  }


void DISTR_binomialprobit::compute_deviance(const double * response,
                   const double * weight,const double * mu,double * deviance,
                   double * scale) const
  {

  if (*weight !=  0)
    {
    if (*response<=0)
      *deviance = -2*log(1-*mu);

    else if (*response > 0)
      *deviance = -2*log(*mu);
    }
  else
    *deviance = 0;

  }


double DISTR_binomialprobit::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  double  mu = randnumbers::Phi2(*linpred);

  double h = 0.39894228*exp(-0.5 * *linpred * *linpred);
  double g = 1/pow(h,2);

  *workingweight = *weight / (mu*(1-mu) * g);


  *workingresponse = *linpred + (*response - mu)/h;

  if (like)
    {

    if (*response > 0)
      return log(mu);
    else
      return log(1-mu);
    }
  else
    {
    return 0;
    }


  }



void DISTR_binomialprobit::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  double  mu = randnumbers::Phi2(*linpred);
  double h = 0.39894228*exp(-0.5 * *linpred * *linpred);
  double g = 1/pow(h,2);

  *workingweight = 1.0 / (mu*(1-mu) * g);

  *workingresponse = *linpred + (*response - mu)/h;

  if (compute_like)
    {

    if (*response > 0)
      like+= log(mu);
    else
      like+= log(1-mu);
    }

  }


void DISTR_binomialprobit::compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  }

void DISTR_binomialprobit::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  }


void DISTR_binomialprobit::update(void)
  {

  double * worklin;
  double * workresp;
  double * workwresp;
  double * weightwork;
  double * workingweightwork;

  register unsigned i;


  if (optionsp->nriter==1)
    {

    weightwork = weight.getV();
    workingweightwork = workingweight.getV();

    for (i=0;i<nrobs;i++,weightwork++,workingweightwork++)
      {
      *workingweightwork = *weightwork;
      }

    }
  else
    {
    if (check_weightsone())
      wtype = wweightsnochange_one;
    else
      wtype = wweightsnochange_constant;
    }



  workresp = response.getV();
  workwresp = workingresponse.getV();
  weightwork = weight.getV();

  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();


  for(i=0;i<nrobs;i++,worklin++,workresp++,weightwork++,workwresp++)
    {

    if (*weightwork != 0)
      {
      if (*workresp > 0)
        *workwresp = trunc_normal2(0,20,*worklin,1);
      else
        *workwresp = trunc_normal2(-20,0,*worklin,1);
      }

    }


  if (utilities)
    {
    workwresp = workingresponse.getV();
    double * wb = FC_latentutilities.beta.getV();

    for (i=0;i<nrobs;i++,workwresp++,wb++)
      {
      *wb = *workwresp;
      }
    FC_latentutilities.update();
    }

  DISTR::update();

  }


void DISTR_binomialprobit::outresults(ofstream & out_stata, ofstream & out_R,
                                      ST::string pathresults)
  {
  if (utilities && pathresults.isvalidfile() != 1)
    {
    unsigned i;

    FC_latentutilities.outresults(out_stata,out_R,"");

    ofstream out(pathresults.strtochar());

    double * wb = FC_latentutilities.betamean.getV();

    for (i=0;i<nrobs;i++,wb++)
      {
      out << *wb << endl;
      }

    }

  }


void DISTR_binomialprobit::get_samples(const ST::string & filename,ofstream & outg) const
  {
  if ((utilities) && (filename.isvalidfile() != 1))
    {
    FC_latentutilities.get_samples(filename,outg);
    }

  }

//------------------------------------------------------------------------------
//--------------------- CLASS DISTRIBUTION_binomialsvm -------------------------
//------------------------------------------------------------------------------


DISTR_binomialsvm::DISTR_binomialsvm(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR(o,r,w)

  {

  predictor_name = "pi";
  outexpectation = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  family = "Binomial_SVM";
  updateIWLS = false;
  }


const DISTR_binomialsvm & DISTR_binomialsvm::operator=(
                                      const DISTR_binomialsvm & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  return *this;
  }


DISTR_binomialsvm::DISTR_binomialsvm(const DISTR_binomialsvm & nd)
   : DISTR(DISTR(nd))
  {
  }


void DISTR_binomialsvm::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Bayesian support vector machine\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_binomialsvm::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {
/*
  if (*weight!=0)
    {
    double mu = randnumbers::Phi2(*linpred);
    if (*response > 0)
      return log(mu);
    else
      return log(1-mu);
    }
  else
    return 0;
*/
  return 0;
  }



double DISTR_binomialsvm::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {
  /*
  double mu = randnumbers::Phi2(*linpred);
  if (*response > 0)
    return log(mu);
  else
    return log(1-mu);
  */
  return 0;
  }


void DISTR_binomialsvm::compute_mu(const double * linpred,double * mu)
  {
//  *mu = randnumbers::Phi2(*linpred);
  }


void DISTR_binomialsvm::compute_deviance(const double * response,
                   const double * weight,const double * mu,double * deviance,
                   double * scale) const
  {
  /*
  if (*weight !=  0)
    {

    if (*response<=0)
      {
      *deviance = -2*log(1-*mu);
      *deviancesat = *deviance;
      }
    else if (*response > 0)
      {
      *deviance = -2*log(*mu);
      *deviancesat = *deviance;
      }

    }
  else
    {
    *deviance = 0;
    *deviancesat = 0;
    }
  */

  }


double DISTR_binomialsvm::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  /*
  double  mu = randnumbers::Phi2(*linpred);

  double h = 0.39894228*exp(-0.5 * *linpred * *linpred);
  double g = 1/pow(h,2);

  *workingweight = *weight / (mu*(1-mu) * g);


  *workingresponse = *linpred + (*response - mu)/h;

  if (like)
    {

    if (*response > 0)
      return log(mu);
    else
      return log(1-mu);
    }
  else
    {
    return 0;
    }
  */
  return 0;
  }



void DISTR_binomialsvm::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  /*
  double  mu = randnumbers::Phi2(*linpred);
  double h = 0.39894228*exp(-0.5 * *linpred * *linpred);
  double g = 1/pow(h,2);

  *workingweight = 1.0 / (mu*(1-mu) * g);

  *workingresponse = *linpred + (*response - mu)/h;

  if (compute_like)
    {

    if (*response > 0)
      like+= log(mu);
    else
      like+= log(1-mu);
    }
  */
  }



void DISTR_binomialsvm::update(void)
  {

  double * worklin;
  double * workwresp;
  double * weightwork;

  double * wweightwork;

  register unsigned i;

  workwresp = workingresponse.getV();
  weightwork = weight.getV();
  wweightwork = workingweight.getV();

  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double a;
  double b=1;
  double lambda;

  for(i=0;i<nrobs;i++,worklin++,weightwork++,workwresp++,
                      wweightwork++)
    {

    if (*weightwork != 0)
      {
      a= 1/fabs(1-(*worklin));
      lambda= rand_inv_gaussian(a,b);

      *wweightwork = 1/lambda;

      *workwresp = 1+lambda;
      }

    }

  DISTR::update();

  }


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_poisson ---------------------------
//------------------------------------------------------------------------------

void DISTR_poisson::check_errors(void)
  {

  if (errors==false)
    {
    unsigned i=0;
    double * workresp = response.getV();
    double * workweight = weight.getV();
    while ( (i<nrobs) && (errors==false) )
      {

      if (*workweight > 0)
        {
        if (*workresp != int(*workresp))
          {
          errors=true;
          errormessages.push_back("ERROR: response must be integer values\n");
          }

        if (*workresp < 0)
          {
          errors=true;
          errormessages.push_back("ERROR: negative response values encountered\n");
          }


        }
      else if (*workweight == 0)
        {
        }
      else
        {
        errors=true;
        errormessages.push_back("ERROR: negative weights encountered\n");
        }

      i++;
      workresp++;
      workweight++;

      }

    }

  }



DISTR_poisson::DISTR_poisson(GENERAL_OPTIONS * o, const datamatrix & r,
                               const datamatrix & w)
  : DISTR(o,r,w)

  {

  predictor_name = "lambda";
  outexpectation = true;

  if (check_weightsone() == true)
    wtype = wweightschange_weightsone;
  else
    wtype = wweightschange_weightsneqone;

  family = "Poisson";
  updateIWLS = true;

  linpredminlimit=-10;
  linpredmaxlimit= 15;

  check_errors();
  }


const DISTR_poisson & DISTR_poisson::operator=(
                                      const DISTR_poisson & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  return *this;
  }


DISTR_poisson::DISTR_poisson(const DISTR_poisson & nd)
   : DISTR(DISTR(nd))
  {
  }


void DISTR_poisson::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_poisson::get_intercept_start(void)
  {
//  return log(response.mean(0));
  return 0;
  }

double DISTR_poisson::cdf(double * res,double * param,double * weight,double * scale)
    {
  /* double cl = gsl_cdf_poisson_P(((*res)-1), *param);
    if ((*res) == 0)
    {
        cl = 0;
    }
    double cr = gsl_cdf_poisson_P((*res), *param);
    double u = randnumbers::uniform_ab(cl, cr);

    return u;
    */
    return 0;
    }

double DISTR_poisson::pdf(double * res,double * param,double * weight,double * scale)
    {
 //   double p = gsl_ran_poisson_pdf((*res), (*param));
//    std::ofstream out;
////  // helpmat1.prettyPrint(out);
//    out.open ("C:\\tmp\\p.raw", std::ofstream::out | std::ofstream::app);
//
//    out << p ;
//    out << " " ;
//    out << *res ;
//    out << " " ;
//    out << *param << endl;
 //   return p;
      return 0;
    }


double DISTR_poisson::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {

  double lambda;
  if (*linpred <= linpredminlimit)
    lambda = exp(linpredminlimit);
  else if (*linpred >= linpredmaxlimit)
    lambda = exp(linpredmaxlimit);
  else
    lambda = exp(*linpred);

  if (*response==0)
    return -(*weight)* lambda;
  else
    return *weight * ((*response) * (*linpred) - lambda);
  }


double DISTR_poisson::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {

  double lambda;
  if (*linpred <= linpredminlimit)
    lambda = exp(linpredminlimit);
  else if (*linpred >= linpredmaxlimit)
    lambda = exp(linpredmaxlimit);
  else
    lambda = exp(*linpred);

  if (*response==0)
    return  -  lambda;
  else
    return  (*response) * (*linpred) - lambda;

  }


void DISTR_poisson::compute_mu(const double * linpred,double * mu)
  {

  if (*linpred <= linpredminlimit)
    *mu = exp(linpredminlimit);
  else if (*linpred >= linpredmaxlimit)
    *mu = exp(linpredmaxlimit);
  else
    *mu = exp(*linpred);

  }



void DISTR_poisson::compute_deviance(const double * response,
                   const double * weight,const double * mu,double * deviance,
                   double * scale) const
  {

  if (*weight==0)
    *deviance = 0;
  else
    {
    if (*response==0)
      {
      *deviance = 2* *weight * *mu;
      }
    else
      {
      double rplusone = *response+1;
      *deviance = -2* *weight*(*response*log(*mu)-*mu-
                      randnumbers::lngamma_exact(rplusone));
      }
    }

  }


double DISTR_poisson::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  double mu;
  if (*linpred <= linpredminlimit)
    mu = exp(linpredminlimit);
  else if (*linpred >= linpredmaxlimit)
    mu = exp(linpredmaxlimit);
  else
    mu = exp(*linpred);

  *workingweight = *weight * mu;

  if (*response==0)
    {
    *workingresponse = *linpred -1;

    if (like)
       return -(*weight) *mu;
     else
       return 0;
    }
  else
    {
    *workingresponse = *linpred + (*response - mu)/mu;

    if (like)
       return *weight * (*response * (*linpred) - mu);
     else
       return 0;
    }

  }


void DISTR_poisson::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {


  if (*linpred <= linpredminlimit)
    *workingweight = exp(linpredminlimit);
  else if (*linpred >= linpredmaxlimit)
    *workingweight = exp(linpredmaxlimit);
  else
    *workingweight = exp(*linpred);

  if (*response==0)
    {
    *workingresponse = *linpred - 1;

     if (compute_like)
       like -=   (*workingweight);
    }
  else
    {
    *workingresponse = *linpred + (*response - (*workingweight))/(*workingweight);

     if (compute_like)
       like+=  *response * (*linpred) - (*workingweight);
    }

  }




void DISTR_poisson::sample_responses(unsigned i,datamatrix & sr)
  {
  double * linpredp;

  if (linpred_current==1)
    linpredp = linearpred1.getV();
  else
    linpredp = linearpred2.getV();

  double * rp = sr.getV()+i;
  double mu;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
    {
    compute_mu(linpredp,&mu);

    *rp = randnumbers::rand_pois(mu);

    }

  }


void DISTR_poisson::sample_responses_cv(unsigned i,datamatrix & linpred,
                                         datamatrix & sr)
  {

  double * linpredp;

  linpredp = linpred.getV();

  double * rp = sr.getV()+i;
  double mu;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
    {
    compute_mu(linpredp,&mu);

    *rp = randnumbers::rand_pois(mu);
    }

  }



//------------------------------------------------------------------------------
//--------------------- CLASS DISTRIBUTION_poisson_ext -------------------------
//------------------------------------------------------------------------------


DISTR_poisson_ext::DISTR_poisson_ext(GENERAL_OPTIONS * o, const datamatrix & r,
                                 double ap, double bp, bool ada,
                               const datamatrix & w)
  : DISTR_poisson(o,r,w)

  {
  a=ap;
  b=bp;
  adapt = ada;
  if (adapt)
    {
    double rmax = response.max(0);
    }
  }


const DISTR_poisson_ext & DISTR_poisson_ext::operator=(
                                      const DISTR_poisson_ext & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_poisson::operator=(DISTR_poisson(nd));
  a = nd.a;
  b = nd.b;
  adapt = nd.adapt;
  return *this;
  }


DISTR_poisson_ext::DISTR_poisson_ext(const DISTR_poisson_ext & nd)
   : DISTR_poisson(DISTR_poisson(nd))
  {
  a = nd.a;
  b = nd.b;
  adapt = nd.adapt;
  }


void DISTR_poisson_ext::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: extended exponential exp(a+b eta)\n");
  optionsp->out("  a = " + ST::doubletostring(a) + "\n");
  optionsp->out("  b = " + ST::doubletostring(b) + "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_poisson_ext::get_intercept_start(void)
  {
//  return log(response.mean(0));
  return 0;
  }


double DISTR_poisson_ext::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {
  return (*weight)*loglikelihood_weightsone(response,linpred);
  }


double DISTR_poisson_ext::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {
  if (*response==0)
    return  - exp(a+b*(*linpred));
  else
    return  b*(*response) * (*linpred) - exp(a+b*(*linpred));
  }


void DISTR_poisson_ext::compute_mu(const double * linpred,double * mu)
  {
  *mu = exp(a+b*(*linpred));
  }


double DISTR_poisson_ext::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  double lambda = exp(a+b*(*linpred));

  *workingweight = *weight * pow(b,2)* lambda;

  if (*response==0)
    {
    *workingresponse = *linpred -1/b;

    if (like)
       return -(*weight) *lambda;
     else
       return 0;
    }
  else
    {
    *workingresponse = *linpred + (*response - lambda)/(b*lambda);

    if (like)
       return (*weight)*(b*(*response) * (*linpred) - lambda);
     else
       return 0;
    }

  }


void DISTR_poisson_ext::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  double lambda = exp(a+b*(*linpred));
  *workingweight = b*b*lambda;

  if (*response==0)
    {
    *workingresponse = *linpred - 1/b;

     if (compute_like)
       like -=   lambda;
    }
  else
    {
    *workingresponse = *linpred + (*response - lambda)/(b*lambda);

     if (compute_like)
       like+=  b*(*response) * (*linpred) - lambda;
    }

  }



//------------------------------------------------------------------------------
//--------------------- CLASS DISTRIBUTION_poisson_extlin ----------------------
//------------------------------------------------------------------------------


DISTR_poisson_extlin::DISTR_poisson_extlin(GENERAL_OPTIONS * o,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR_poisson(o,r,w)

  {
  }


const DISTR_poisson_extlin & DISTR_poisson_extlin::operator=(
                                      const DISTR_poisson_extlin & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_poisson::operator=(DISTR_poisson(nd));
  return *this;
  }


DISTR_poisson_extlin::DISTR_poisson_extlin(const DISTR_poisson_extlin & nd)
   : DISTR_poisson(DISTR_poisson(nd))
  {
  }


void DISTR_poisson_extlin::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: exp(eta) for eta < 0 and eta+1 for eta >= 0\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_poisson_extlin::get_intercept_start(void)
  {
//  return log(response.mean(0));
  return 0;
  }


double DISTR_poisson_extlin::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {
  return (*weight)*loglikelihood_weightsone(response,linpred);
  }


double DISTR_poisson_extlin::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {
  if (*linpred < 0)
    {
    if (*response==0)
      return  - exp(*linpred);
    else
      return  (*response) * (*linpred) - exp(*linpred);
    }
  else
    {
    if (*response==0)
      return  - (*linpred+1);
    else
      return  (*response) * log(*linpred + 1) - (*linpred+1);
    }

  }


void DISTR_poisson_extlin::compute_mu(const double * linpred,double * mu)
  {
  if (*linpred < 0)
    *mu = exp(*linpred);
  else
    *mu = *linpred + 1;

  }


double DISTR_poisson_extlin::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  if (*linpred < 0)
    {

    *workingweight = (*weight)*exp(*linpred);

    if (*response==0)
      {
      *workingresponse = *linpred - 1;

       if (like)
         return -(*workingweight);
      }
    else
      {
      *workingresponse = *linpred + (*response - (*workingweight))/(*workingweight);

       if (like)
         return  (*response) * (*linpred) - (*workingweight);
      }

    }
  else
    {
    *workingweight = (*weight)/(*linpred+1);

    *workingresponse = *response-1;

    if (like)
      {

      if (*response==0)
        return  - (*weight)*(*linpred+1);
      else
        return  (*weight)*((*response) * log(*linpred + 1) - (*linpred+1));
      }

    }

  return 0;
  }


void DISTR_poisson_extlin::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  if (*linpred < 0)
    {

    *workingweight = exp(*linpred);

    if (*response==0)
      {
      *workingresponse = *linpred - 1;

       if (compute_like)
         like -=   (*workingweight);
      }
    else
      {
      *workingresponse = *linpred + (*response - (*workingweight))/(*workingweight);

       if (compute_like)
         like+=  (*response) * (*linpred) - (*workingweight);
      }

    }
  else
    {
    *workingweight = 1/(*linpred+1);

    *workingresponse = *response-1;

    if (compute_like)
      {

      if (*response==0)
        like+=  - (*linpred+1);
      else
        like+=  (*response) * log(*linpred + 1) - (*linpred+1);
      }

    }

  }



} // end: namespace MCMC



