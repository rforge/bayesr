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




#if !defined (GAUSSIANHETEROSKEDASTIC_INCLUDED)

#define GAUSSIANHETEROSKEDASTIC_INCLUDED

#include"../export_type.h"
#include"distribution.h"

#if !defined(M_PI)
#define M_PI		3.14159265358979323846
#endif

namespace MCMC
{

class __EXPORT_TYPE DISTRIBUTION_gaussianh : public DISTRIBUTION
  {

  protected:

  unsigned nrcat;               // number of categories of the response

  // FUNCTION: standardise
  // TASK: standardises the response and the offset
  //       sets scalesave.transform = trmult*trmult (!!!)

  void standardise(void);



  public:

  // DEFAULT CONSTRUCTOR

  DISTRIBUTION_gaussianh(void) : DISTRIBUTION() {}

   // CONSTRUCTOR

   DISTRIBUTION_gaussianh(const double & a,const datamatrix & b,
                   MCMCoptions * o, const datamatrix & r,
                   const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTRIBUTION_gaussianh(const DISTRIBUTION_gaussianh & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_gaussianh &
   operator=(const DISTRIBUTION_gaussianh & nd);

   // DESTRUCTOR

   ~DISTRIBUTION_gaussianh() {}


   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * response,double * linpred,
                       double * weight,const int & i) const;

  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred' and stores
  //       the result in 'mu'

  void compute_mu(const double * linpred,double * mu) const;

  void compute_mu_notransform(const double * linpred,double * mu) const;


  // FUNCTION: compute_deviance

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * deviancesat,
                        const datamatrix & scale,const int & i) const;


  double compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col=0) const;


 double compute_IWLS(double * response,double * linpred, double * weight,
                      const int & i,double * weightiwls,double * tildey,
                      bool weightyes, const unsigned & col=0);

 void compute_IWLS_weight_tildey(double * response,double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col=0);

double compute_gmu(double * linpred,const unsigned & col=0) const;


  // FUNCTION: outoptions
  // TASK: writing options of the distribution

  void outoptions(void);


  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  void update_predict(void);

  void outresults(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode for the scale parameter and the
  //       intercept

  bool posteriormode(void);

/*
  bool posteriormode_converged_fc(const datamatrix & beta,
                                  const datamatrix & beta_mode,
                                  const unsigned & itnr);
*/

  void compute_iwls(void);


  };


 }




 // end: namespace MCMC

#endif

