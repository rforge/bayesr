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




#include "distr_mixture.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//------------------------- CLASS DISTR_gaussianmixture ------------------------
//------------------------------------------------------------------------------

void DISTR_gaussianmixture::define_knots(double left, double dist)
  {

  unsigned nr = alpha.rows();

  unsigned i;
  for(i=0;i<nr;i++)
    m(i,0) = left+dist*i;

  // ofstream out("c:\\bayesx\\testh\\results\\m.res");
  // m.prettyPrint(out);

  }


DISTR_gaussianmixture::DISTR_gaussianmixture(const double & a,const double & b,
                               GENERAL_OPTIONS * o, const datamatrix & r,
                               const ST::string & ps,
                               const datamatrix & w)
  : DISTR_gaussian(a,b,o,r,ps,w)
  {

  family = "Mixture of Gaussians";

  nrknots = 5;

  alpha = datamatrix(nrknots,1,1/double(nrknots));
  alphasample = FC(o,"",nrknots,1,ps + "_alpha.raw");
  alpha_prob = datamatrix(nrknots,1,1);

  rho = statmatrix<int>(nrobs,1,0);

  s2 = 1;
  sigma2 = 1;

  m = datamatrix(alpha.rows(),1,0);
  define_knots(-4.5,9/(double(nrknots)-1));

  wtype = wweightschange_weightsneqone;
  }


DISTR_gaussianmixture::DISTR_gaussianmixture(const DISTR_gaussianmixture & nd)
  : DISTR_gaussian(DISTR_gaussian(nd))
  {
  nrknots = nd.nrknots;
  alpha = nd.alpha;
  alphasample = nd.alphasample;
  alpha_prob = nd.alpha_prob;
  m = nd.m;
  s2 = nd.s2;
  rho = nd.rho;
  }


const DISTR_gaussianmixture & DISTR_gaussianmixture::operator=(
                                      const DISTR_gaussianmixture & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gaussian::operator=(DISTR_gaussian(nd));
  nrknots = nd.nrknots;
  alpha = nd.alpha;
  alphasample = nd.alphasample;
  alpha_prob = nd.alpha_prob;
  m = nd.m;
  s2 = nd.s2;
  rho = nd.rho;
  return *this;
  }


/*
void DISTR_gaussianmixture::compute_mu(const double * linpred,double * mu,
                                bool notransform)
  {

  }
*/


/*
void DISTR_gaussianmixture::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           double * scale) const
{
}
*/

/*
double DISTR_gaussianmixture::loglikelihood(double * res, double * lin,
                                     double * w) const
  {

  }
*/

/*
double DISTR_gaussianmixture::loglikelihood_weightsone(double * res,
 double * lin) const
  {

  }


double DISTR_gaussianmixture::compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like)
  {

  }


void DISTR_gaussianmixture::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  }


void DISTR_gaussianmixture::compute_iwls_wweightsnochange_constant(
                                              double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  }


void DISTR_gaussianmixture::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  }
*/


void DISTR_gaussianmixture::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: identity\n");
  optionsp->out("  Hyperparameter a: " + ST::doubletostring(a_invgamma,6) + "\n");
  optionsp->out("  Hyperparameter b: " + ST::doubletostring(b_invgamma,6) + "\n");
  optionsp->out("\n");
  optionsp->out("\n");

  }


void DISTR_gaussianmixture::update(void)
  {

  // Zuorndung samplen

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * workresp;
  workresp = workingresponse.getV();

  double * workorigresp;
  workorigresp = response.getV();



  unsigned i;
  int j;
  double sigma = sqrt(sigma2);
  double help=0;
  double u;
  for (i=0;i<nrobs;i++,workresp++,worklin++,workorigresp++)
    {
    help=0;

    for (j=0;j<nrknots;j++)
      {
      help += exp(alpha(j,0) - 1.0/(2*s2*sigma2) *
                        pow(*workresp-*worklin,2)       );

      alpha_prob(j,0) = help;
      }



    u = uniform()*alpha_prob(nrknots-1,0);

    j = 0;
    while (u >= alpha_prob(j,0))
      {
      j++;
      }

    rho(i,0) = j;


    *workresp = *workorigresp - sigma*m(rho(i,0),0);

    }

  }


bool DISTR_gaussianmixture::posteriormode(void)
  {
  return true;
  }


/*
bool DISTR_gaussianmixture::posteriormode(void)
  {

  }


void DISTR_gaussianmixture::outresults(ST::string pathresults="")
  {

  }

double DISTR_gaussianmixture::get_scalemean(void)
  {

  }


void DISTR_gaussianmixture::sample_responses(unsigned i,datamatrix & sr)
  {

  }


void DISTR_gaussianmixture::outresults_predictive_check(datamatrix & D,datamatrix & sr)
  {

  }
*/





} // end: namespace MCMC



