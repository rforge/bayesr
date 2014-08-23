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



#if !defined (DISTRIBUTION_nbinomial_INCLUDED)

#define DISTRIBUTION_nbinomial_INCLUDED

#include"../export_type.h"
#include "distribution.h"
#include "Random.h"




namespace MCMC
{

////////////////////////////////////////////////////////////////////////////////
// Distribution family for count data with overdispersion.
// Implemented are three models:
//
//      * NB:   negative binomial model
//              y ~ NB(mu, s)
//              p(y) = G(y+s)/(G(s)*G(y+1))(s/(s+mu))^s (mu/(s+mu))^y
//              family = nbinomial
//              distopt = nb
//
//      * POGA: Poisson-Gamma model
//              y ~ PO(nu*mu)
//              nu ~ G(s, s)
//              family = nbinomial
//              distopt = poga
//
//      * POIG: Poisson-Inverse Gauss model
//              y ~ PO(nu*mu)
//              nu ~ IGauss(1, s)
//              family = nbinomial
//              distopt = poig
//
//  Hierarchical versions of the POGA and POIG model are also implemented.
//  Just add "hierarchical" to the input line of the programm.
//
//  You have two options for the proposal of s and both work very similar.
//  propopt = unif or propopt = gamma
//
//  Summary: NB and POGA works well, POIG does not work properly (why?).
//  Only NB is described in the BayesX manual.
//
//  Literature:
//
//      * Cameron and Trivedi (1998), "Regression analysis for count data"
//        Cambridge University Press, New York.
//
//      * Osuna, Leyre (2004), "Semiparametric Bayesian count data models"
//        Dissertation an der LMU, München.
//
////////////////////////////////////////////////////////////////////////////////



enum vertopt {nb,poga,poig};
enum propscale {unif,gam};


class __EXPORT_TYPE DISTRIBUTION_nbinomial : public DISTRIBUTION
  {

   protected:

   // include here additional private variables needed for the new
   // distribution

   bool oversize;           // Boolean variable. If True, then number of
                            // observations is larger than 500.
                            // Concerns only the number of multiplicative
                            // effects, for which samplings will be saved.

   datamatrix accept;       // Vector of accepted update-steps. It has length
                            // nrobs+2. The first entry is reserved for the
                            // scale parameter in the POGA model and the
                            // following nrobs-entries are reserved for the
                            // multiplicative random effects of the POIG model.
                            // The last entry is reserved for the hierarchical
                            // intercept.

   datamatrix nu;           // Vector of lenght nrobs to store the multiplicative
                            // random effects of the POGA and POIG models.

   FULLCOND nusave;         // Fullconditional object to store the samples of the
                            // random effects of the POGA and POIG models.
                            // Only when oversize=False!!!

   FULLCOND nusavekfz;      // Fullconditional object to store the samples of the
                            // random effects of the POGA and POIG models.
                            // Only when oversize=True!!!

   datamatrix hierint;      // To store the intercept in a Hierarchcial model.
                            // See Section 2.3 from diss for the Theory.

   FULLCOND hierintsave;    // Fullconditional object to store the samples of the
                            // intercept in a Hierarchcial model.

   datamatrix pvar;         // Vector of proposal windows/variances. It has legnth
                            // nrobs+2. The first entry is reserved for the scale
                            // parameter in the POGA model and the following
                            // nrobs-entries are reserved for the multiplicative
                            // random effects of the POIG model. The last entry
                            // is reserved for the hierarchical intercept.


   double a_pri;            // Stores the hyperparameter a for the Gamma Priori
                            // of the scale parameter. It is no sampled in the model.

   datamatrix b_pri;        // To store the hyperparameter b for the Gamma Priori
                            // of the scale parameter. It will be sampled in the
                            // model.

   FULLCOND b_pri_save;     // Fullconditional object to store the samples of the
                            // hyperparameter b for the Gamma Priori of the scale
                            // parameter.

   double prop_var;         // Option. Start value for the proposal
                            // window/variance of the scale parameter. Not
                            // used...

   vertopt ver;             // Option. It stores the distributional assumption
                            // for the response variable. Possible values:
                            //          nb = negative binomial
                            //          poga = poisson-gamma
                            //          poig = poisson-inverse gauss

   propscale pscale;        // Option. It stores the proposal assumption for
                            // the scale parameter. Possible values:
                            //          gam = gamma proposal
                            //          unif = uniform proposal

   bool hierarchical;       // Option. Boolean variable. Possible values:
                            //          True = hierarchical model
                            //          False = "normal" model


   datamatrix sum_nu;       // Helpvariable for the calculations.

   datamatrix sum2_nu;      // Helpvariable for the calculations.




   public:

   //  the following public functions must be implemented


   // DEFAULT CONSTRUCTOR

    DISTRIBUTION_nbinomial(void);


   // CONSTRUCTOR WITHOUT OFFSET

   DISTRIBUTION_nbinomial(const double & a,const double & b,
                          const double & pv,const vertopt & vo,
                          const propscale & psc,bool hie,
                          MCMCoptions * o,
                          const datamatrix & r,
                          const ST::string & p,const ST::string & ps,
                          const datamatrix & w=datamatrix());


   // CONSTRUCTOR WITH OFFSET
   DISTRIBUTION_nbinomial(const double & a,const double & b,
                          const double & pv,const vertopt & vo,
                          const propscale & psc,bool hie,
                          const datamatrix & offset, MCMCoptions * o,
                          const datamatrix & r,
                          const ST::string & p,const ST::string & ps,
                          const datamatrix & w=datamatrix());


void create(MCMCoptions * o, const double & a,
                                    const double & b, const double & pv,
                                    const vertopt & vo, const propscale & psc,
                                    bool hie, const ST::string & ps);




   // COPY CONSTRUCTOR

   DISTRIBUTION_nbinomial(const DISTRIBUTION_nbinomial & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_nbinomial & operator=(const DISTRIBUTION_nbinomial & nd);

   // DESTRUCTOR

   ~DISTRIBUTION_nbinomial() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a !!single observation!! depending
   // on the model we have chosen and omiting all the terms that do not depend
   // on the predictor. From diss:
   // negative binomial = Formula (2.4)
   // POGA and POIG = Formula (2.11)

  double loglikelihood(double * response,double * linpred,
                       double * weight,const int & i) const;

  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred' and stores
  //       the result in 'mu'
  // The link function is allways the loglink!

  void compute_mu(const double * linpred,double * mu) const;

  // FUNCTION: compute_deviance
  // TASK: computes the individual deviance for the Negative Binomial
  // model, independently of the chosen model!

  void compute_deviance(const double * response,const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale, const int & i) const;


   // FUNCTION: compute_weight
   // TASK: computes the weights for iteratively weighted least squares
   //       w_i = weight_i/scale*[b''(theta_i) g'^2(mu_i)]^{-1}
   // It differentiate between posteriormode and MCMC iterations.
   // Posteriormode calculations are bases only on the NB model.
   // Therefore only NB-weights are used.
   // MCMC iterations used the corresponding weights to the chosen models.
   // See for the derivation of the used values: diss!!!
   // NB model: Formula (B.5)
   // POGA and POIG models: Formula (B.7)

  double compute_weight(double * linpred,double * weight,
                        const int & i,const unsigned & col=0) const;

  // FUNCTION: compute_gmu
  // TASK: compute g'(eta_i) = 1/h'(eta_i)

  double compute_gmu(double * linpred,const unsigned & col=0) const;

  void compute_mu_notransform(const double * linpred,double * mu) const;

  // FUNCTION: compute_IWLS
  // There are three blocks in this function.
  // First block: calculate IWLS weights. The same Formulae as in the
  // compute_weight function.
  // Second block: calculate tildey = term/weight. The "terms" are
  // the first derivatives of the loglikelihood without the x_ij.
  // See diss, Appendix B!!!
  // Third block: it returns the loglikelihood of the models.

  double compute_IWLS(double * response,double * linpred,double * weight,
                      const int & i,double * weightiwls,double * tildey,
                      bool weightyes, const unsigned & col=0);

  // FUNCTION: compute_IWLS_weight_tildey
  // the two first blocks of the function compute_IWLS!!

  void compute_IWLS_weight_tildey(double * response,double * linpred,
                              double * weight,const int & i,double * weightiwls,
                              double * tildey,const unsigned & col=0);


  // FUNCTION: outoptions
  // TASK: writing options of the distribution

  void outoptions(void);

  void outresults(void);

  // FUNCTION: update
  // TASK: calls the update functions for the scale parameter and its hyperparameter
  // b, and for the multiplicative random effects and hierarchical intercept,
  // when present in the model.

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode for the scale parameter but
  // only when oversize = True!!! otherwise there were numerical problems
  // with the calculations.
  // For oversize=True: Cameron & Trivedi (1998)

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  // Updates the multiplicative random effects, depending on the model chosen.

  double update_nu(void);

  // Updates the hierarchical intercept, when hierarchical=True. Aktually
  // it only stores the current value in the current place. The true atuallisation
  // happens in the file "mcmc_const.cpp", in the function
  // "FULLCOND_const_nbinomial::update_hierint(void)".

  double update_hierint(void);

  // Updates the scale parameter, depending on the model chosen.

  double update_scale(void) const;

  // Updates the hyperparameter b of the prior for the scale parameter.
  // The update step is independent of the model chosen!

  double update_b_pri(void);

  // This function samples a new value for the scale parameter, depending on the
  // proposal, that we have chosen, and calculates the proposal_ratio for this
  // new value.

  double proposal_scale(void) const;

  // Tuning function for the acceptance rates of the M-H algorithms.
  // It will be called only in the burnin phase and each 100 iterations.
  // It tries to achive acdeptance rates between 0.3 and 0.6.

  double pwork_tunin(unsigned i) const;

  // Logarithm of the density f(x)~G(a, b). Not used...

  double log_gamma(double & x, double & a, double & b) const;

  // Vector of multiplicative random effects: nu
  // nu ~ product from 1 to nrobs of G(s, s) (s=scale)
  // Function calculates log(g(nu|s_neu))-log(g(nu|s))

  double log_gamma_likelihood(double &s, double &s_neu) const;

  // Vector of multiplicative random effects: nu
  // nu ~ product from 1 to nrobs of G(s, s/hierint) (s=scale)
  // Function calculates log(g(nu|s_neu, hierint))-log(g(nu|s, hierint))

  double log_gamma_likelihood_hier(double &s, double &s_neu) const;

  // Logarithm of the Quotient: G_(a,b)(x_prop)/G_(a, b)(x)
  // Nor used...

  double log_gamma_quot(double & x, double & x_prop, const double & a,
                        const double & b) const;

  // not implemented and not used...

  double log_gamma_prop(double & x, double & x_prop) const;

    // loggamma-function

  double lgamma(const double & xx) const;

  //Loglikelihood-difference of a NB distribution with:
  // log(NB(lambda, s_neu))-log(NB(lambda, s))

  double log_nbin(const double & s_neu, const double & s) const;

  // log[Gamma(g_prop+y)/Gamma(g_prop)]-log[Gamma(g+y)/Gamma(g)] with y integer!

  double log_gamma_function_quot(const double & g, const double & g_prop,
                                 const unsigned & y) const;

  // The following functions are needed for the hierarchical versions of the models,
  // in order to "conect" this subclass with the file "mcmc_const.cpp".

  const double & get_sum_nu(void) const
    {
    return sum_nu(0,0);
    }

  const double & get_sum2_nu(void) const
    {
    return sum2_nu(0,0);
    }

  const double & get_hierint(void) const
    {
    return hierint(0,0);
    }

  const int get_distopt(void) const
    {
    if(ver==poga) return 1;
    else return 2;
    }

  const double & get_pvar(void) const
    {
    return pvar(nrobs+1, 0);
    }

  void initialize_hierint(double & inter)
    {
    hierint(0,0) = inter;
    }

  void exchange_hierint(double & inter)
    {
    hierint(0,0) += inter;
    }

  void exchange_accept(void)
    {
        double *acceptwork = accept.getV();
        acceptwork += nrobs +1;
        *acceptwork += 1;

    }

  void add_nu(double m) const;


  };   // end: DISTRIBUTION_nbinomial


}  // end: namespace MCMC

#endif
