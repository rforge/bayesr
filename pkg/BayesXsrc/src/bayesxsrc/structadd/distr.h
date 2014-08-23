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



#if !defined (DISTR_INCLUDED)
#define DISTR_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include <fstream> 	//used for class logit_fruehwirt

#if !defined(M_PI)
#define M_PI		3.14159265358979323846
#endif

namespace MCMC
{

using randnumbers::rand_invgamma;
using randnumbers::rand_normal;
using randnumbers::uniform;
using randnumbers::trunc_normal;
using randnumbers::trunc_normal2;
using randnumbers::truncnormal;
using randnumbers::kssample;
using randnumbers::rand_gamma;
using randnumbers::rand_inv_gaussian;

using randnumbers::invPhi2;

/*
1. workingweights ändern sich, weights ungleich 1

2. workingweights ändern sich, weights gleich eins

3. workingweights ändern sich nicht und sind konstant

4. workingweights ändern sich nicht und sind eins
*/

enum weighttype{wweightschange_weightsneqone,wweightschange_weightsone,
wweightsnochange_constant,wweightsnochange_one};

enum msetype{noMSE,quadraticMSE,checkMSE};

enum auxiliarytype{auxcurrent,auxpostmean};

class __EXPORT_TYPE DISTR
  {

  protected:

  // FUNCTION: check_workingweights_one
  // TASK: checks if all workingweights are one (returns true if this is the
  //       case)

  bool check_weightsone(void);

  // FUNCTION: compute_nrzeroweights
  // TASK: determines the number of zero weights and returns the result

  unsigned compute_nrzeroweights(void);


  GENERAL_OPTIONS * optionsp;         // pointer to general MCMC options object


  public:

  bool maindistribution;
  bool predict_mult;

  datamatrix * FCpredict_betamean;

  bool optionbool1;
  ST::string option1;


  double sigma2;

  bool updateIWLS;

  ST::string family;              // name of the distribution
  ST::string familyshort;
  unsigned hlevel;
  ST::string equationtype;

  unsigned nrobs;                 // Number of observations


  datamatrix response;                // Response

  datamatrix workingresponse;         // Working response, tilde y
  ST::string responsename;            // Name of the response

  ST::string offsetname;          // name of offset variable


  datamatrix weight;              // Weightvariable for weighted regression

  ST::string weightname;          // Name of the weightvariable

  datamatrix workingweight;       // Working weight (workingweight = weight
                                  // in the constructor)

  weighttype wtype;               // weight type: default is
                                  // wweightschange_weightsneqone, i.e.
                                  // workingweights change and weights are
                                  // not equal to one
  bool weightsone;                // true if weights are one for all
                                  // observations

  unsigned nrzeroweights;         // number of zero weights

  datamatrix linearpred1;          // Linear predictor
  datamatrix linearpred2;          // Proposed linear predictor
  int linpred_current;

  bool outpredictor;
  bool outexpectation;
  ST::string predictor_name;
  unsigned predstart_mumult;

  double meaneffect;

  //----------------------------------------------------------------------------
  //---------------------- linpredlimits for save estimation -------------------
  //----------------------------------------------------------------------------

  double linpredminlimit;
  double linpredmaxlimit;

  // FUNCTION: check_linpred
  // TASK: checks whether current predictor vector is within linpredlimits

  bool check_linpred(bool current = true);

  void changelimits(double min,double max);

  //----------------------------------------------------------------------------
  //---------------------- linpredlimits for save estimation -------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //------------------------- auxiliary variables ------------------------------
  //----------------------------------------------------------------------------

  datamatrix helpmat1;              // Stores auxiliary quantities
  datamatrix helpmat2;              // Stores auxiliary quantities
  datamatrix helpmat3;              // Stores auxiliary quantities

  double helpquantity1;             // Stores auxiliary quantities
  double helpquantity2;             // Stores auxiliary quantities
  double helpquantity3;             // Stores auxiliary quantities

  //----------------------------------------------------------------------------
  //---------------------- end: auxiliary variables ----------------------------
  //----------------------------------------------------------------------------


  void swap_linearpred(void);


  double trmult;                 // multiplicative constant for hyperparameters


//------------------------------------------------------------------------------
//------------------------------- ERRORS ---------------------------------------

  bool errors;

  vector<ST::string> errormessages;

  virtual void check_errors(void);



//------------------------------------------------------------------------------
//--------------------------- CONSTRUCTORS -------------------------------------
//------------------------------------------------------------------------------

  // DEFAULT CONSTRUCTOR

  DISTR(void)
    {
    }

  // CONSTRUCTOR1
  // TASK: initializes data
  //       response = r
  //       weight = w
  //       nrobs = r.rows()

  DISTR(GENERAL_OPTIONS * o,const datamatrix & r,
               const datamatrix & w=datamatrix());


  // COPY CONSTRUCTOR

  DISTR(const DISTR & d);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR & operator=(const DISTR & d);

  // DESTRUCTOR

  ~DISTR() {}

  //----------------------------------------------------------------------------
  //------------------------------ WRITING OPTIONS -----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: outoptions
  // TASK: writing options

  virtual void outoptions(void);

  //----------------------------------------------------------------------------
  //-------------- OBTAINING SAMPLES OF DISTRIBUTION PARAMETERS  ---------------
  //----------------------------------------------------------------------------

  virtual void get_samples(const ST::string & filename,ofstream & outg) const;

  //----------------------------------------------------------------------------
  //-------------- OBTAINING SAMPLES OF DISTRIBUTION PARAMETERS  ---------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //----------------------- INITIALIZE AN INTERCEPT ----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: get_intercept_start
  // TASK: returns starting value for the intercept (if specified)

  virtual double get_intercept_start(void);


  //----------------------------------------------------------------------------
  //---------------------------- COMPUTING THE CDF -----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: cdf
  // TASK: computes the cdf for a single observation

  virtual double cdf(double * res,double * param,double * weight,double * scale)
    {
    return 0;
    }

  virtual double cdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    return 0;
    }


  virtual double pdf(double * res,double * param,double * weight,double * scale)
    {
    return 0;
    }

  virtual double pdf_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux)


    {
    return 0;
    }

  virtual double compute_quantile_residual(double * res,double * param,double * weight,
                                    double * scale);

  virtual double compute_quantile_residual_mult(vector<double *> response,
                                         vector<double *> param,
                                         vector<double *> weight,
                                          vector<datamatrix *> aux);


   double compute_quadr(void);
   double compute_quadr_mult(void);
   double compute_log(double * res,double * param,double * weight,
                                        double * scale);
   double compute_log_mult(vector<double *> response,
                          vector<double *> param,
                          vector<double *> weight,
                          vector<datamatrix *> aux);
   double compute_spherical(void);
   double compute_spherical_mult(void);
   double compute_CRPS(void);
   double compute_CRPS_mult(void);

  //----------------------------------------------------------------------------
  //----------------------- COMPUTING THE LOGLIKELIHOOD ------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: loglikelihood
  // TASK: computes the loglikelihood for a single observation

  virtual double loglikelihood(double * res,double * lin,double * weight)
    {
    return 0;
    }

  virtual double loglikelihood_weightsone(double * res,double * lin)
    {
    return 0;
    }

  // FUNCTION: loglikelihood
  // TASK: computes the complete loglikelihood for all observations

  double loglikelihood(const bool & current=true);

  // FUNCTION: loglikelihood
  // TASK: computes the loglikelihood for observations between begin and end
  //       response, weights, predicor stored in responsep,workingweightp,
  //       linpredp

  double loglikelihood(int & begin,
                       int & end, statmatrix<double *> & responsep,
                       statmatrix<double *> & workingweightp,
                       statmatrix<double *> & linpredp);


  //----------------------------------------------------------------------------
  //------------------------------- COMPUTE mu ---------------------------------
  //----------------------------------------------------------------------------

  virtual void compute_mu(const double * linpred,double * mu);

  virtual void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);


  //----------------------------------------------------------------------------
  //------------------------------- COMPUTE param ------------------------------
  //----------------------------------------------------------------------------

  virtual void compute_param(const double * linpred,double * param);

  virtual void compute_param_mult(vector<double *>  linpred,double * param);

  //----------------------------------------------------------------------------
  //------------------------------- COMPUTE deviance ---------------------------
  //----------------------------------------------------------------------------

  virtual void compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * scale) const;


  virtual void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix *> aux);


  //----------------------------------------------------------------------------
  //------------------------------- COMPUTE MSE --------------------------------
  //----------------------------------------------------------------------------

  virtual double compute_MSE(const double * response, const double * weight,
                             const double * linpred, msetype t, double v);

  virtual void compute_MSE_all(datamatrix & meanpred, double & MSE,
                               double & MSEzeroweight, unsigned & nrzeroweights,
                               msetype & t, double & v);

  //----------------------------------------------------------------------------
  //----------------------------- IWLS Algorithm -------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: compute_IWLS (for one observation)
  // TASK: computes the iwls weights (will be stored in workingweight),
  //       tildey=predicor+(y-mu)g'(mu) (stored in workingresponse) and
  //       the loglikelihood (will be returned)

  //       type: wweightschange_weightsneqone

  virtual double compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse,const bool & like)
    {
    return 0;
    }

  // FUNCTION: compute_IWLS (for one observation)
  // TASK: computes the iwls weights (will be stored in workingweight),
  //       tildey=predicor+(y-mu)g'(mu) (stored in workingresponse) and
  //       the loglikelihood stored in like (only if compute_like = true)
  //       assumes that weighs=1 (for all observations)

  virtual void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
    {
    }


  // FUNCTION: compute_IWLS (for one observation)
  // TASK: computes tildey=predicor+(y-mu)g'(mu) (stored in workingresponse) and
  //       the loglikelihood stored in like (only if compute_like = true)
  //       assumes that workingweighs=constant (for all observations), i.e.
  //       they are not recomputed in the function
  //       wweightsnochange_constant

  virtual void compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
    {
    }


  // FUNCTION: compute_IWLS (for one observation)
  // TASK: computes tildey=predicor+(y-mu)g'(mu) (stored in workingresponse) and
  //       the loglikelihood stored in like (only if compute_like = true)
  //       assumes that workingweighs=1 (for all observations), must be set
  //       to one in advance

  virtual void compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
    {
    }


  // FUNCTION: compute_IWLS (for the whole dataset
  // TASK:

  double compute_iwls(const bool & current,const bool & like);

  void compute_iwls(const bool & current,datamatrix & likelihood,
                    statmatrix<unsigned> & ind);

  // FUNCTION: compute_IWLS (
  // TASK: computes the iwls weights (will be stored in workingweight),
  //       tildey=(y-mu)g'(mu) (stored in workingresponse) and
  //       the loglikelihood (will be returned) for the begin - end observation
  //       in the pointer vectors

  double compute_iwls_loglikelihood(int & begin,
                                 int & end, statmatrix<double *> & responsep,
                                 statmatrix<double *> & workingresponsep,
                                 statmatrix<double *> & weightp,
                                 statmatrix<double *> & workingweightp,
                                 statmatrix<double *> & linpredp);


  double compute_iwls_loglikelihood_sumworkingweight(
         int & begin,int & end, statmatrix<double *> & responsep,
         statmatrix<double *> & workingresponsep,statmatrix<double *> & weightp,
         statmatrix<double *> & workingweightp, statmatrix<double *> & linpredp,
         datamatrix & intvar2,double & sumworkingweight);


  //----------------------------------------------------------------------------
  //----------------------- ACCESS TO SCALE PARAMETER --------------------------
  //----------------------------------------------------------------------------

  virtual double get_scale(void);
  virtual double get_scalemean(void);
  virtual void update_scale_hyperparameters(datamatrix & h);
  virtual datamatrix * get_auxiliary_parameter(auxiliarytype t = auxcurrent);

  //----------------------------------------------------------------------------
  //----------------------- POSTERIORMODE FUNCTIONS ----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  virtual bool posteriormode(void);

  virtual void posteriormode_end(void);


  //----------------------------------------------------------------------------
  //--------------------------- UPDATE FUNCTIONS -------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: update
  // TASK: base function for inherited classes,
  //       should update the scale parameter
  //       the base function updates the estimated mean and variance
  //       of the scale parameter only

  virtual void update(void);


  // FUNCTION: update
  // TASK: base function for inherited classes,
  //       may be used to update quantities that have been changed while updating
  //       FC's and that are required for other equations (e.g. in ZIP models)

  virtual void update_end(void);


  //----------------------------------------------------------------------------
  //---------------------------- SAMPLE RESPONSES ------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: sample responses
  // TASK: samples vector of responses based on current predictor and scale
  //       parameter,stores results in the i-th col of sr

  virtual void sample_responses(unsigned i,datamatrix & sr);

  virtual void sample_responses_cv(unsigned i,datamatrix & linpred,
                                   datamatrix & sr);

  virtual void outresults_predictive_check(datamatrix & D,datamatrix & sr);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: outresults
  // TASK: writes estimation results for the scale parameter
  //       estimated mean and variance

  virtual void outresults(ofstream & out_stata, ofstream & out_R,ST::string pathresults="");

  // FUNCTION: reset
  // TASK: resets linpred (all values to 0)

  void reset(void);


  }; // end: class DISTR


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussian : public DISTR
  {

  protected:

  double a_invgamma;                    // hyperparameter a (for the inverse
                                        // gamma distribution of the scale
                                        // parameter, i.e. sigma^2
  double b_invgamma;                    // hyperparameter b

  FC FCsigma2;

  double nrlasso;
  double nrridge;
  double lassosum;
  double ridgesum;

  public:


   // DEFAULT CONSTRUCTOR

   DISTR_gaussian(void) : DISTR()
     {
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b
   // N(linpred,sigma2/weight)

   DISTR_gaussian(const double & a,const double & b,GENERAL_OPTIONS * o,
                  const datamatrix & r,const ST::string & ps,
                         const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_gaussian(const DISTR_gaussian & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_gaussian & operator=(const DISTR_gaussian & nd);

   // DESTRUCTOR

   ~DISTR_gaussian() {}


   void get_samples(const ST::string & filename,ofstream & outg) const;

   double compute_MSE(const double * response,
                          const double * weight,
                          const double * linpred, msetype t, double v);

  void compute_mu(const double * linpred,double * mu);

  void compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * scale) const;

  double get_intercept_start(void);

  double loglikelihood(double * res,
                       double * lin,
                       double * w);

  double loglikelihood_weightsone(double * res,double * lin);

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

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  bool posteriormode(void);

  void outresults(ofstream & out_stata, ofstream & out_R,ST::string pathresults="");

  double get_scalemean(void);

  void sample_responses(unsigned i,datamatrix & sr);

  void sample_responses_cv(unsigned i,datamatrix & linpred, datamatrix & sr);

  void outresults_predictive_check(datamatrix & D,datamatrix & sr);

  // FUNCTION: update_scale_hyperparameters
  // TASK: updates parameters for lasso, ridge etc.
  //       h(0,0) = type, 1 =ridge, 2=lasso
  //       h(1,0) = nrridge/nrlasso
  //       h(2,0) = lassosum/ridgesum

  void update_scale_hyperparameters(datamatrix & h);

  };

//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_quantreg ----------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE DISTR_quantreg : public DISTR_gaussian
  {

  protected:

  double quantile;
  double xi,xi2;
  double num;
  double sigma02;

  public:

   // DEFAULT CONSTRUCTOR

   DISTR_quantreg(void) : DISTR_gaussian()
     {
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b

   DISTR_quantreg(const double & a,const double & b,GENERAL_OPTIONS * o,
                  const datamatrix & r,const ST::string & ps,double & quant,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_quantreg(const DISTR_quantreg & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_quantreg & operator=(const DISTR_quantreg & nd);

   // DESTRUCTOR

   ~DISTR_quantreg() {}


   double compute_MSE(const double * response,const double * weight,
                                      const double * linpred, msetype t,double v);


/*
  void compute_mu(const double * linpred,double * mu, bool notransform);

  void compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           double * scale) const;

  double loglikelihood(double * res,
                       double * lin,
                       double * w) const;

  double loglikelihood_weightsone(double * res,double * lin) const;

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

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  /*
  bool posteriormode(void);

  void outresults(ST::string pathresults="");

  double get_scalemean(void);

  void sample_responses(unsigned i,datamatrix & sr);

  void outresults_predictive_check(datamatrix & D,datamatrix & sr);
  */
  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_hetgaussian -------------------------
//------------------------------------------------------------------------------

// heteroscedastic gaussian for simultaneous estimation of mean and variance
// variance estimation via family vargaussian below

class __EXPORT_TYPE DISTR_hetgaussian : public DISTR_gaussian
  {

  protected:

  bool sigma2const;

  public:

  datamatrix * FCpredict_betamean_vargaussian;

  datamatrix weightoriginal;

  // DEFAULT CONSTRUCTOR

  DISTR_hetgaussian(void) : DISTR_gaussian()
    {
    }

  // CONSTRUCTOR1

  DISTR_hetgaussian(double a, double b, GENERAL_OPTIONS * o,
                     const datamatrix & r,
                     const ST::string & ps, const bool sc,
                     const datamatrix & w=datamatrix());

  // COPY CONSTRUCTOR

  DISTR_hetgaussian(const DISTR_hetgaussian & nd);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_hetgaussian & operator=(const DISTR_hetgaussian & nd);

  // DESTRUCTOR

  ~DISTR_hetgaussian() {}

  double compute_MSE(const double * response, const double * weight,
                         const double * linpred, msetype t, double v);

  void compute_MSE_all(datamatrix & meanpred, double & MSE,
                               double & MSEzeroweight, unsigned & nrzeroweights,
                               msetype & t, double & v);

  void update(void);

  bool posteriormode(void);

  void outresults(ofstream & out_stata, ofstream & out_R,
                  ST::string pathresults="");

  };


//------------------------------------------------------------------------------
//--------------------------- DISTR_vargaussian --------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_vargaussian  : public DISTR
  {

  protected:

  double sigma2old;

  public:

  DISTR_hetgaussian * dgaussian;

//------------------------------------------------------------------------------
//------------------------------- ERRORS ---------------------------------------

  // void check_errors(void);

//------------------------------------------------------------------------------
//--------------------------- CONSTRUCTORS -------------------------------------
//------------------------------------------------------------------------------

  // DEFAULT CONSTRUCTOR

  DISTR_vargaussian(void)
    {
    }

  // CONSTRUCTOR1
  // TASK: initializes data
  //       response = r
  //       weight = w
  //       nrobs = r.rows()

  DISTR_vargaussian(GENERAL_OPTIONS * o,const datamatrix & r);


  // COPY CONSTRUCTOR

  DISTR_vargaussian(const DISTR_vargaussian & d);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_vargaussian & operator=(const DISTR_vargaussian & d);

  // DESTRUCTOR

  ~DISTR_vargaussian() {}

  //----------------------------------------------------------------------------
  //----------------------- COMPUTING THE LOGLIKELIHOOD ------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: loglikelihood
  // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * res,double * lin,double * weight);

  //----------------------------------------------------------------------------
  //------------------------------- COMPUTE mu ---------------------------------
  //----------------------------------------------------------------------------

  void compute_mu(const double * linpred,double * mu);

  //----------------------------------------------------------------------------
  //------------------------------- COMPUTE MSE --------------------------------
  //----------------------------------------------------------------------------

//  double compute_MSE(const double * response, const double * weight,
//                             const double * linpred, msetype t, double v);



  //----------------------------------------------------------------------------
  //----------------------------- IWLS Algorithm -------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: compute_IWLS (for one observation)
  // TASK: computes the iwls weights (will be stored in workingweight),
  //       tildey=predicor+(y-mu)g'(mu) (stored in workingresponse) and
  //       the loglikelihood (will be returned)

  //       type: wweightschange_weightsneqone

  double compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse,const bool & like);


  void outoptions(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  // FUNCTION: update
  // TASK: base function for inherited classes,
  //       should update the scale parameter
  //       the base function updates the estimated mean and variance
  //       of the scale parameter only

  void update(void);

  }; // end: class DISTR_vargaussian


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_loggaussian -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_loggaussian : public DISTR_gaussian
  {

  protected:


  public:

   // DEFAULT CONSTRUCTOR

   DISTR_loggaussian(void) : DISTR_gaussian()
     {
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b

   DISTR_loggaussian(const double & a,const double & b,GENERAL_OPTIONS * o,
                  const datamatrix & r,const ST::string & ps,
                         const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_loggaussian(const DISTR_loggaussian & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_loggaussian & operator=(const DISTR_loggaussian & nd);

   // DESTRUCTOR

   ~DISTR_loggaussian() {}

  void compute_mu(const double * linpred,double * mu);

  double compute_MSE(const double * response, const double * weight,
                     const double * linpred, msetype t,double v);


  void compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * scale) const;

  void sample_responses(unsigned i,datamatrix & sr);

  void sample_responses_cv(unsigned i,datamatrix & linpred, datamatrix & sr);

  void outresults_predictive_check(datamatrix & D,datamatrix & sr);

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian_exp ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussian_exp : public DISTR_gaussian
  {

  protected:

  // FUNCTION: standardise
  // TASK: standardises the response and the offset
  //       sets scalesave.transform = trmult*trmult (!!!)

  // void standardise(void);


  public:

   // DEFAULT CONSTRUCTOR

   DISTR_gaussian_exp(void) : DISTR_gaussian()
     {
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b

   DISTR_gaussian_exp(const double & a,const double & b,GENERAL_OPTIONS * o,
                  const datamatrix & r,const ST::string & ps,
                         const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_gaussian_exp(const DISTR_gaussian_exp & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_gaussian_exp & operator=(const DISTR_gaussian_exp & nd);

   // DESTRUCTOR

   ~DISTR_gaussian_exp() {}


  void compute_mu(const double * linpred,double * mu);

  void compute_param(const double * linpred,double * param);

  double loglikelihood(double * res,
                       double * lin,
                       double * w);

  double compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like);

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  bool posteriormode(void);

  void sample_responses(unsigned i,datamatrix & sr);

  void outresults_predictive_check(datamatrix & D,datamatrix & sr);

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian_mult -----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussian_mult : public DISTR_gaussian_exp
  {

  protected:

  // FUNCTION: standardise
  // TASK: standardises the response and the offset
  //       sets scalesave.transform = trmult*trmult (!!!)

//  void standardise(void);


  public:

  void set_mult(bool & m);

   // DEFAULT CONSTRUCTOR

   DISTR_gaussian_mult(void) : DISTR_gaussian_exp()
     {
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b

   DISTR_gaussian_mult(const double & a,const double & b,GENERAL_OPTIONS * o,
                  const datamatrix & r,const ST::string & ps,
                         const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_gaussian_mult(const DISTR_gaussian_mult & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_gaussian_mult & operator=(const DISTR_gaussian_mult & nd);

   // DESTRUCTOR

   ~DISTR_gaussian_mult() {}

//  void compute_mu(const double * linpred,double * mu, bool notransform);
  void compute_mu(const double * linpred,double * mu);


  double loglikelihood(double * res,
                       double * lin,
                       double * w);

  double compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like);

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  bool posteriormode(void);

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian_re -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_gaussian_re : public DISTR_gaussian
  {

  protected:

  public:

   // DEFAULT CONSTRUCTOR

   DISTR_gaussian_re(void) : DISTR_gaussian()
     {
     }

   // CONSTRUCTOR1

   DISTR_gaussian_re(GENERAL_OPTIONS * o,
                  const datamatrix & r,
                  const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_gaussian_re(const DISTR_gaussian_re & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_gaussian_re & operator=(const DISTR_gaussian_re & nd);

   // DESTRUCTOR

   ~DISTR_gaussian_re() {}

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  bool posteriormode(void);

  void outresults(ofstream & out_stata, ofstream & out_R,ST::string pathresults="");

  void outoptions(void);

  void get_samples(const ST::string & filename,ofstream & outg) const;

  void check_errors(void);

  };



} // end: namespace MCMC


#endif
