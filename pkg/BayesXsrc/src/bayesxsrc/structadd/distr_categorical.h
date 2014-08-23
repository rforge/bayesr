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



#if !defined (DISTRcategorical_INCLUDED)
#define DISTRcategorical_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include"distr.h"


namespace MCMC
{


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_binomial ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_binomial : public DISTR
  {

  protected:

  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

   DISTR_binomial(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_binomial(GENERAL_OPTIONS * o, const datamatrix & r,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_binomial(const DISTR_binomial & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_binomial & operator=(const DISTR_binomial & nd);

   // DESTRUCTOR

   ~DISTR_binomial() {}

  void compute_mu(const double * linpred,double * mu);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * scale) const;

  double get_intercept_start(void);

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void outoptions(void);

  void sample_responses(unsigned i,datamatrix & sr);

  void sample_responses_cv(unsigned i,datamatrix & linpred,
                                   datamatrix & sr);

  };


//------------------------------------------------------------------------------
//------------------ CLASS: DISTRIBUTION_logit_fruehwirth-----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_logit_fruehwirth : public DISTR_binomial
{
 protected:

	int H;
  datamatrix SQ;
  datamatrix weights_mixed;


 public:

 void check_errors(void);

 	// DEFAULT CONSTRUCTOR

 	DISTR_logit_fruehwirth(void) : DISTR_binomial()
 		{
 		}

 	// CONSTRUCTOR1
 	DISTR_logit_fruehwirth(const int h, GENERAL_OPTIONS * o,
  											const datamatrix r,
                        const datamatrix & w=datamatrix());


 	// COPY CONSTRUCTOR
 	DISTR_logit_fruehwirth(const DISTR_logit_fruehwirth & nd);


 	// OVERLOADED ASSIGNMENT OPERATOR
 	const DISTR_logit_fruehwirth & operator=(const DISTR_logit_fruehwirth & nd);


 	// DESTRUCTOR
 	~DISTR_logit_fruehwirth()
 	{
 	}
////////////////////

/*
 	double compute_MSE();
*/

//  basis class implementation
// 	void compute_mu(const double * linpred,double * mu);

// basis class implementation
//  void compute_deviance(const double * response, const double * weight,
//                        const double * mu,double * deviance,
//                        double * deviancesat, double * scale) const;

// basis class implementation
//  double loglikelihood(double * response, double * linpred,
//                       double * weight) const;


// basis class implementation
//  double loglikelihood_weightsone(double * response, double * linpred) const;


// basis class implementation
/*
  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

*/


	void outoptions();

	// FUNCTION: update
	// TASK: uptdates the scale parameter

	void update(void);

	bool posteriormode(void);

  // no results
  // void outresults();

  // not required
	// double get_scalemean(void);

//  basis class implementation
// 	void sample_responses(unsigned i,datamatrix & sr);

//  basis class implementation
//	void sample_responses_cv();

//  basis class implementation
//	void outresults_predictive_check();

};



//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_binomialprobit ----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_binomialprobit : public DISTR
  {

  protected:

  bool utilities;
  FC FC_latentutilities;

  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

   DISTR_binomialprobit(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_binomialprobit(GENERAL_OPTIONS * o, const datamatrix & r,const bool ut,
                        const ST::string & ps,const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_binomialprobit(const DISTR_binomialprobit & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_binomialprobit & operator=(const DISTR_binomialprobit & nd);

   // DESTRUCTOR

   ~DISTR_binomialprobit() {}

  void compute_mu(const double * linpred,double * mu);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * scale) const;

  double get_intercept_start(void);

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void outoptions(void);

  void update(void);

  void outresults(ofstream & out_stata, ofstream & out_R,
                  ST::string pathresults);

  void get_samples(const ST::string & filename,ofstream & outg) const;

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_binomialsvm -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_binomialsvm : public DISTR
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

   DISTR_binomialsvm(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_binomialsvm(GENERAL_OPTIONS * o, const datamatrix & r,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_binomialsvm(const DISTR_binomialsvm & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_binomialsvm & operator=(const DISTR_binomialsvm & nd);

   // DESTRUCTOR

   ~DISTR_binomialsvm() {}

  void compute_mu(const double * linpred,double * mu);

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance, double * scale)
                        const;

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);



  void outoptions(void);

  void update(void);

  };



//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_poisson -----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_poisson : public DISTR
  {

  protected:

  public:

  void check_errors(void);

   // DEFAULT CONSTRUCTOR

   DISTR_poisson(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_poisson(GENERAL_OPTIONS * o, const datamatrix & r,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_poisson(const DISTR_poisson & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_poisson & operator=(const DISTR_poisson & nd);

   // DESTRUCTOR

   ~DISTR_poisson() {}

  void compute_mu(const double * linpred,double * mu);

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance,
                        double * scale) const;

  double get_intercept_start(void);

  double cdf(double * res,double * param,double * weight,double * scale);

  double pdf(double * res,double * param,double * weight,double * scale);

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void sample_responses(unsigned i,datamatrix & sr);

  void sample_responses_cv(unsigned i,datamatrix & linpred,
                                   datamatrix & sr);

  void outoptions(void);

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_poisson_ext -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_poisson_ext : public DISTR_poisson
  {

  protected:

  double a;
  double b;
  bool adapt;

  public:

   // DEFAULT CONSTRUCTOR

   DISTR_poisson_ext(void) : DISTR_poisson()
     {
     }

   // CONSTRUCTOR

   DISTR_poisson_ext(GENERAL_OPTIONS * o, const datamatrix & r,
                     double ap, double bp, bool ada,
                     const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_poisson_ext(const DISTR_poisson_ext & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_poisson_ext & operator=(const DISTR_poisson_ext & nd);

   // DESTRUCTOR

   ~DISTR_poisson_ext() {}

  void compute_mu(const double * linpred,double * mu);

  double get_intercept_start(void);

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void outoptions(void);

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_poisson_extlin ----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_poisson_extlin : public DISTR_poisson
  {

  protected:

  public:

   // DEFAULT CONSTRUCTOR

   DISTR_poisson_extlin(void) : DISTR_poisson()
     {
     }

   // CONSTRUCTOR

   DISTR_poisson_extlin(GENERAL_OPTIONS * o, const datamatrix & r,
                        const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_poisson_extlin(const DISTR_poisson_extlin & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_poisson_extlin & operator=(const DISTR_poisson_extlin & nd);

   // DESTRUCTOR

   ~DISTR_poisson_extlin() {}

  void compute_mu(const double * linpred,double * mu);

  double get_intercept_start(void);

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void outoptions(void);

  };


} // end: namespace MCMC


#endif
