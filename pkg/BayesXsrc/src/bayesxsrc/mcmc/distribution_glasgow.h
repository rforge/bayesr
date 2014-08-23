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


#if !defined (DISTRIBUTION_INCLUDED)

#define DISTRIBUTION_INCLUDED

#include<statmat.h>
#include<random.h>
#include<mcmc.h>
#include<fullcond.h>


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



class DISTRIBUTION
  {

  protected:

  bool nosamples;                 // true if samples should not be stored
                                  //

  MCMCoptions * optionsp;         // pointer to general MCMC options object

  ST::string family;              // name of the distribution

  //----------------------------------------------------------------------------
  //-------------------- Variables for the scale parameter ---------------------
  //----------------------------------------------------------------------------

  bool scaleexisting;             // true, if scale parameter is present
  datamatrix scale;               // current value of the scale parameter
                                  // set to 0.1


  double  acceptancescale;        // number of accepted iterations


  datamatrix scale_mode;

  FULLCOND Scalesave;            // handling and saving sampled scaleparameters

  ST::string pathresultsscale;

  //----------------------------------------------------------------------------


  unsigned nrobs;                 // Number of observations
  unsigned nrobsmweightzero;      // Number of observations-observations with
                                  // weight zero

  datamatrix response;            // Response
  double * responsep;             // Pointer to a particular element in
                                  // response
  ST::string responsename;        // Name of the response
  datamatrix trmult;              // multiplicative constant with which
                                  // the response has been transformed

  double addinterceptsample;      // constant to be added to  the samples of the
                                  // intercept

  datamatrix weight;              // Weightvariable for weighted regression
  double sumweight;               // sum of weights
  bool changingweight;            // is true, if weights change in every
                                  // iteration
  ST::string weightname;          // Name of the weightvariable
  double * weightp;               // Pointer to a particular element in
                                  // weight. access via functions set_weightp
                                  // and get_weight

  ST::string offsetname;

  datamatrix linearpred;          // Linear predictor
  datamatrix linearpredprop;      // Proposed linear predictor
  datamatrix * linpred_current;   // Pointer that contains adress of current
                                  // predictor
  datamatrix * linpred_proposed;  // Pointer that contains adress of proposed
                                  // predictor
  double * linpredp_current;      // Pointer to a particular element in
                                  // linpred_current
  double * linpredp_proposed;     // Pointer to a particular element in
                                  // linpred_proposed


//------------------------------------------------------------------------------
//------------------ VARIABLES FOR COMPUTING THE POSTERIOR MODE ----------------
//------------------------------------------------------------------------------

  datamatrix tildey;              // working response
  datamatrix weightiwls;          // weights for iwls
                                  // w_i = weight_i/scale*[b''(theta_i)
                                  // g'^2(mu_i)]^{-1}
  datamatrix workingres;          // Working residual

  double * workingresp;           // Pointer to a particular element in
                                  // workingres. Access via functions
                                  // set_weightp and get_weight

//------------------------------------------------------------------------------
// ------------------------- For prediction ------------------------------------
//------------------------------------------------------------------------------

  bool predict;                   // if predict = true the sampling mean of the
                                  // predictor will be updated in every
                                  // iteration
  bool predictfull;               // if predictfull = true, for the mean mu
                                  // credible intervals etc. will be
                                  // additionally computed
  unsigned firstobs;              // computes credible intervals for mu only
                                  // for the first 'firstobs' observations.
  ST::string predictpath;         // Path for storing predicted mu, residuals
                                  // etc.
  ST::string predictfullpath;     // Path for storing predicted mu's incl.
                                  // credible intervals etc.
  ST::string deviancepath;        // Path for storing the deviance and DIC
  datamatrix linpredmean;         // sample mean of the linear predictor
  datamatrix mumean;              // sample mean of mu
  datamatrix deviancemean;        // sample mean of individual deviances
  datamatrix deviancemean_sat;    // sample mean of individual deviances
  datamatrix deviance;            // sample of the deviance

  FULLCOND musave;                // handles samples for mu

  datamatrix * Dp;
  vector<ST::string> Dnames;

  bool predictresponse;
  datamatrix predictindicator;

  vector<ST::string> results_latex;

//------------------------------------------------------------------------------
//----------------------------- FOR OUTRESULTS ---------------------------------
//------------------------------------------------------------------------------

  datamatrix interceptsample;     // stores sampled intercept values
                                  // used for function outresults (transformed
                                  // parameters) and cumulative probit models

  double interceptold;

//------------------------------------------------------------------------------
//----------------------------- FOR MISSING VALUES -----------------------------
//------------------------------------------------------------------------------

  vector<FULLCOND *> fcmissing;
  statmatrix<unsigned> missingpos;
  FULLCOND MissingSave;
  datamatrix missingind;
  ST::string pathmissing;


//------------------------------------------------------------------------------
//--------------------------- ERRORMESSAGES ------------------------------------
//------------------------------------------------------------------------------

  vector<ST::string> errors;      // contains errormesages


  void create(MCMCoptions * o, const datamatrix & r, const datamatrix & w);

  public:

//------------------------------------------------------------------------------
//--------------------------- CONSTRUCTORS -------------------------------------
//------------------------------------------------------------------------------

  // DEFAULT CONSTRUCTOR

  DISTRIBUTION(void)
    {
    nrobs = 0;

    family = "unknown";

    scaleexisting = true;
    scale = datamatrix(1,1,0.1);
    }

  // CONSTRUCTOR1
  // TASK: initializes data
  //       scale = 1 x 1 matrix
  //       scale(0,0) = 0.1 (may be reset with set_scale)
  //       response = r
  //       weight = w
  //       nrobs = r.rows()
  //       set names of the variables to their default
  //       Y for response, W for weightvariable
  //       default names may be changes with init_names

  DISTRIBUTION(MCMCoptions * o,const datamatrix & r,
               const datamatrix & w=datamatrix(),const ST::string & pr="",
               const ST::string & ps="");

  // CONSTRUCTOR2
  // TASK: initializes also an offset

  DISTRIBUTION(const datamatrix & offset,MCMCoptions * o,const datamatrix & r,
               const datamatrix & w=datamatrix(),const ST::string & pr="",
               const ST::string & ps="");

  // COPY CONSTRUCTOR

  DISTRIBUTION(const DISTRIBUTION & d);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DISTRIBUTION & operator=(const DISTRIBUTION & d);

  // DESTRUCTOR

  ~DISTRIBUTION() {}


  //----------------------------------------------------------------------------
  //----------------------- ACCESS TO ERROR MESSAGES ---------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: geterrors
  // TASK: returns errormessages

  const vector<ST::string> & geterrors(void)
    {
    return errors;
    }

  // FUNCTION: outerrors
  // TASK: writes errors

  void outerrors(void);


  //----------------------------------------------------------------------------
  //------------------------------ WRITING OPTIONS -----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: outoptions
  // TASK: writing options

  virtual void outoptions(void);

  //----------------------------------------------------------------------------
  //--------------- ACCSESS TO NAME OF THE FAMILY, e.g. poisson ----------------
  //----------------------------------------------------------------------------

  // FUNCTION: get_family
  // TASK: returns the name of the distribution

  const ST::string & get_family(void) const
    {
    return family;
    }

  // FUNCTION: get_responsename
  // TASK: returns the name of the response variable

  const ST::string & get_responsename(void) const
    {
    return responsename;
    }

  // FUNCTION: get_nrpar
  // TASK: returns the number of parameters

  unsigned get_nrpar(void);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: set_predict
  // TASK: indicates that the mean of the linearpredictor should be updated
  //       in every iteration, mean linearpred, deviance, etc. are written
  //       to 'path'

  void set_predict(const ST::string & path,const ST::string & pathdev,
                   datamatrix * p,vector<ST::string> & Dn);

  // FUNCTION: set_predictfull
  // TASK: indicates that predicted mu's for the first 'fo' observations
  //       should be computed (incl. credible intervals)
  //       results are written to 'path'

  void set_predictfull(const ST::string & pathsample,
                       const ST::string & path,const unsigned & fo);

  virtual void update_predict(void);

  void set_predictresponse(const datamatrix & pr);

  // FUNCTION: init_names
  // TASK: initializes the names of the response variable
  //      (stored in datamatrix response),
  //      the names of the covariates (stored in datamatrix data)
  //      and the name of the weight variable
  //      stored in datamarix weight
  //      and the name of the offset variable

  void init_names(const ST::string & rn, const ST::string & wn = "",
                  const ST::string & on="");


  void init_offset(const datamatrix & o);

  //----------------------------------------------------------------------------
  //-------------------- ACCESS TO THE SCALE PARAMETER -------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: set_scale
  // TASK: reinitializes the scale parameter
  //       if scale is a scalar

  void set_scale(const double & newscale)
    {
    scale(0,0) = newscale;
    }


  // FUNCTION: get_scale_rows
  // TASK: returns the number of rows of scale

  unsigned get_scale_rows(void)
    {
    return scale.rows();
    }


  // FUNCTION: get_scale_cols
  // TASK: returns the number of columns of scale

  unsigned get_scale_cols(void)
    {
    return scale.cols();
    }


  // FUNCTION: get_scale
  // TASK: returns the current value of the scale parameter

  virtual const double & get_scale(const unsigned & r=0,
                                   const unsigned & c=0) const
    {
    return scale(0,0);
    }

  // FUNCTION: get_scale_sample
  // TASK: stores the sampled scale parameters in file 'file'

  virtual ST::string get_scale_sample(void) const;

  // FUNCTION: get_scaleexisting
  // TASK: returns the value of scaleexisting, if scaleexisting == true
  //       the distribution has a nonconstant scale parameter

  const bool & get_scaleexisting(void)
    {
    return scaleexisting;
    }

  // FUNCTION: compute_autocor_scale
  // TASK: computes autocorrelation function for lags 1 - 'lag'.

  datamatrix compute_autocor_scale(const unsigned & lag,const unsigned & row,
                                   const unsigned & col) const;

  //----------------------------------------------------------------------------
  //----------------------- Interceptsample ------------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: set_interceptsample
  // TASK: sets column th column of the interceptsample
  //       is done from a fullcond_const object

  void set_interceptsample(datamatrix & s,unsigned & column);


  //----------------------------------------------------------------------------
  //----------------------- COMPUTING THE LOGLIKELIHOOD ------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: loglikelihood
  // TASK: computes the loglikelihood for a single observation
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  virtual double loglikelihood(double * res,double * lin,double * weight,
                               const int & i) const
    {
    return 0;
    }

  // FUNCTION: loglikelihood
  // TASK: computes the complete loglikelihood for all observations
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  double loglikelihood(const bool & current=true) const;

  // FUNCTION: loglikelihood
  // TASK: computes the loglikelihood for observations between 'beg' and 'end'
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  double loglikelihood(const unsigned & beg,const unsigned & end,
                       const statmatrix<int> & index,
                       const bool & current=true) ;

  double loglikelihood2(const unsigned & beg,const unsigned & end,
                        const statmatrix<int> & index,
                        const statmatrix<int> & index2,
                        const bool & current=true) ;

  //----------------------------------------------------------------------------
  //------------------------------ COMUPUTING mu -------------------------------
  //----------------------------------------------------------------------------

   // FUNCTION: compute_mu
   // TASK: computes mu for a new linear predictor 'linpred' and stores
   //       the result in 'mu'

   virtual void compute_mu(const double * linpred,double * mu) const
     {
     }

   // FUNCTION: compute_mu
   // TASK: computes mu for a new linear predictor 'linpred' and stores
   //       the result in 'mu'

  void compute_mu(const datamatrix & linpred, datamatrix & mu) const;

  //----------------------------------------------------------------------------
  //-------------------------- COMPUTING the deviance  -------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: compute_deviance
  // TASK: computes the individual deviance and saturated deviance;
  //       deviance = -2log(Liklihood(y|mu))
  //       deviancesat = -2log(likelihood(y|mu)+2log(Likelihood(y|y))
  //       if the response is transformed before or during MCMC, i.e.
  //       response = trmult*response_current+addintercept
  //       scale = trmult^2*scale_current
  //       the residual will be automatically and correctly retransformed
  // IMPORTANT: scale will be retransformed during the computation of the
  //            deviance, i.e. scale must be passed in the form scale_current

  virtual void compute_deviance(const double * response,const double * weight,
                                   const double * mu, double * deviance,
                                   double * deviancesat,
                                   const datamatrix & scale,const int & i) const
    {
    }

  void compute_deviance(double & deviance,double & deviancesat);


  virtual double compute_gcv(double & df)
    {
    return 0;
    }

  //----------------------------------------------------------------------------
  //--------------- functions for maximizing the loglikelihood -----------------
  //----------------------------------------------------------------------------

   // FUNCTION: compute_weight
   // TASK: computes the weights for iteratively weighted least squares
   //       proposals
   //       w_i = weight_i/scale*[b''(theta_i) g'^2(mu_i)]^{-1}

  virtual double compute_weight(double * linpred,double * weight,
                                const int & i, const unsigned & col=0) const
    {
    return 0;
    }

  void compute_weight(datamatrix & w,const unsigned & col,const bool
                      & current=true) const;

  void compute_weight(datamatrix & w,const unsigned & beg, const unsigned & end,
                       const statmatrix<int> & index, const unsigned & col = 0);

  double compute_sumweight(const unsigned & col,const bool & current=true) const;

  double compute_sumweight(const unsigned & beg, const unsigned & end,
                       const statmatrix<int> & index, const unsigned & col,
                       const bool & current=true);

  double compute_sumweight2(const unsigned & beg, const unsigned & end,
                            const statmatrix<int> & index,
                            const statmatrix<int> & index2,
                            const unsigned & col,
                            const bool & current=true);


  // FUNCTION: compute_gmu
  // TASK: compute g'(eta_i) = 1/h'(eta_i)

  virtual double compute_gmu(double * linpred,const unsigned & col=0) const
    {
    return 0;
    }

  // FUNCTION: fisher
  // TASK: computes  data' W data  and stores the result in 'XWX' , where
  //       'data' is the designmatrix,
  //       'w' = diag [ weight_i/scale*[b''(theta_i) g'^2(mu_i)]^{-1} ]
  //       (will be computed in the function)
  // IMPORTANT: 'XWX' and 'w' must have proper dimensions, i.e. XWX must have
  //             dimension data.cols() x data.cols() and w must have dimension
  //             nrconst x response.cols()

  void fisher(datamatrix & XWX,datamatrix & w,datamatrix & data,
              const unsigned & col,const bool & current=true) const;

  void fisher(datamatrix & XWX,datamatrix & w,vector<unsigned> & posbeg,
              vector<unsigned> & posend, statmatrix<int> & index,
              unsigned & refind,const unsigned & c=0) const;

  double fisher2(const unsigned & beg, const unsigned & end,
                            const statmatrix<int> & index,
                            datamatrix & data, const unsigned & col,
                            const bool & current=true) const;

  // FUNCTION: tilde_y
  // TASK: computes tildey = eta + (y-mu)g'(mu)

  void tilde_y(datamatrix & tildey,const bool & current = true);

  // FUNCTION: tilde_y
  // TASK: computes tildey = m + (y-mu)g'(mu)      ( für IWLS )

  virtual void tilde_y(datamatrix & tildey,datamatrix & m,
              const unsigned & col,const bool & current,const datamatrix & w);

  //----------------------------------------------------------------------------
  //--------------- ACCESS TO DATA (RESPONSE,WEIGHT,COVARIATES) ----------------
  //----------------------------------------------------------------------------

  // FUNCTION: get_nrobs
  // TASK: returns the number of observation

  const unsigned & get_nrobs(void) const
    {
    return nrobs;
    }

  // FUNCTION: get_trmult
  // TASK: returns trmult

  const double & get_trmult(const unsigned & c=0)
    {
    return trmult(c,0);
    }

  // FUNCTION: get_addinterceptsample
  // TASK: returns addinterceptsample

  const double & get_addinterceptsample(void)
    {
    return addinterceptsample;
    }

  void update_intercept(double m,unsigned c)
    {
    }
 

  void set_responsep(const unsigned & r, const unsigned & c)
    {
    responsep = response.getV()+r*response.cols()+c;
    }


  // FUNCTION: get_response
  // TASK: returns the i,j th element of the response matrix

  const double & get_response(const unsigned & i,const unsigned & j) const
    {
    return response(i,j);
    }

  // FUNCTION: get_response
  // TASK: returns the response matrix

  const datamatrix & get_response(void) const
    {
    return response;
    }

  // FUNCTION: get_responsedim
  // TASK: returns the dimension of the response variable

  unsigned get_responsedim(void) const
    {
    return response.cols();
    }

  // FUNCTION: set_weightp
  // TASK: sets the pointer weightp to the first element in weight

  void set_weightp(const unsigned & c=0)
    {
    weightp = weight.getV()+c;
    }


  double * get_weightp(void)
    {
    return weight.getV();
    }


  // FUNCTION: get_weight
  // TASK: returns the i,j th element of the weight matrix

  const double & get_weight(const unsigned & i, const unsigned & j) const
    {
    return weight(i,j);
    }

  inline const double & get_weight(const int & i)
    {
    weightp+=i;
    return *weightp;
    }

  // FUNCTION: get_weight
  // TASK: returns the weight matrix

  const datamatrix & get_weight(void) const
    {
    return weight;
    }

  // FUNCTION: get_changingweight
  // TASK: returns true, if weights change in every iteration of the sampling
  //       scheme

  bool get_changingweight(void)
    {
    return changingweight;
    }


  // FUNCTION: set_changingweight
  // TASK: returns true, if weights change in every iteration of the sampling
  //       scheme

  void set_changingweight(void)
    {
    changingweight=true;
    }


  virtual double * get_integral_ti(void)
    {
    double * temp;
    return temp;
    }


  //----------------------------------------------------------------------------
  //---------------------- ACCESS TO THE LINEAR PREDICTOR ----------------------
  //----------------------------------------------------------------------------

  void set_linpredp_current(const unsigned & c)
    {
    linpredp_current = (*linpred_current).getV()+c;
    }

  //

  void set_linpredp_current(const unsigned & r, const unsigned & c)
    {
    linpredp_current = (*linpred_current).getV()+r*linearpred.cols()+c;
    }

  // FUNCTION: get_linearpred
  // TASK: returns the i,j th. element of linearpred
  //       if current = true the value(s) of the current linear predictor
  //       will be returned, else the value(s) of the proposed
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  const double & get_linearpred(const unsigned & i,const unsigned & j,
                                const bool & current = true) const
    {
    if (current)
      return (*linpred_current)(i,j);
    else
      return (*linpred_proposed)(i,j);
    }

  // FUNCTION: get_linearpred
  // TASK:  returns the linear predictor
  //       if current = true the value(s) of the current linear predictor
  //       will be returned, else the value(s) of the proposed
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  datamatrix & get_linearpred(const bool & current = true) const
    {
    if (current)
      return *linpred_current;
    else
      return *linpred_proposed;
    }

  // FUNCTION: assign
  // TASK: if current = true, then the proposed linear predictor will be
  //       assigned to the current
  //       if current = false, then the current linear predictor will be
  //       assigned to the proposed
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  void assign(const bool & current = true);

  void assigncol(const unsigned & col,const bool & current = true);

  // FUNCTION: assign
  // TASK: if current = true, then the proposed linear predictor will be
  //       assigned to the current
  //       if current = false, then the current linear predictor will be
  //       assigned to the proposed
  // IMPORTANT: for univariate response only

  void assign(const unsigned & beg,const unsigned & end,
              const statmatrix<int> & index,const bool & current=true);

  // FUNCTION: substr_linearpred
  // TASK: substracts the matrix m from the linear predictor
  //       if current = true the value(s) of the current linear predictor
  //       will be returned, else the value(s) of the proposed
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  void substr_linearpred(const datamatrix & m,const bool & current = true);

  // FUNCTION: substr_linearpred
  // TASK: substracts the matrix m from the col th-column of the linear
  //       predictor
  //       if current = true the value(s) of the current linear predictor
  //       will be returned, else the value(s) of the proposed
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  void substr_linearpred_m(const datamatrix & m,const unsigned & col,
                           const bool & current = true);

  // FUNCTION: add_linearpred
  // TASK: adds the matrix m to the linear predictor
  //       if current = true the value(s) of the current linear predictor
  //       will be used, else the value(s) of the proposed
  //       col is the column of the linearpredictor to change
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  void add_linearpred(const datamatrix & m,const bool & current = true);


  // FUNCTION: add_linearpred
  // TASK: adds the value m to the linear predictor
  //       if current = true the value(s) of the current linear predictor
  //       will be used, else the value(s) of the proposed
  //       row and col are the row and the colum of the linearpredictor to change
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  void add_linearpred(const double & m, unsigned & row, unsigned & col,
  const bool & current = true);

  // FUNCTION: add_linearpred
  // TASK: adds the matrix m to the col th column of the linear predictor
  //       if current = true the value(s) of the current linear predictor
  //       will be used, else the value(s) of the proposed
  //       col is the column of the linearpredictor to change
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  void add_linearpred_m(const double & m,const unsigned & col,
                        const bool & current = true);

  void add_linearpred_m(const datamatrix & m,const unsigned & col,
                        const bool & current = true);


  void add_linearpred(const double & m,const bool & current = true);

  // FUNCTION: add_linearpred
  // TASK: adds the matrix m to the linear predictor between beg and end
  //       if current = true the value(s) of the current linear predictor
  //       will be used, else the value(s) of the proposed
  //       col is the column of the linear predictor to change
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  void add_linearpred(const datamatrix & m,const unsigned & beg,
                        const unsigned & end,const statmatrix<int> & index,
                        const unsigned & col,const bool & current = true);

  // FUNCTION: addtocurrent
  // TASK: adds to the current linear predictor (pointer linpred_current)
  //       the values of 'm' and stores the result in linpred_proposed
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  void addtocurrent(const datamatrix & m);

  void addtocurrent(const double & m);

  // FUNCTION: addtocurrentcol
  // TASK: adds to the current linear predictor (pointer linpred_current)
  //       between 'beg' and 'end' the value 'm' and stores the result in
  //       linpred_proposed
  //       col is the column of the linear predictor to change
  // FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  void addtocurrentcol(const double & m,const unsigned & beg,
                       const unsigned & end, const statmatrix<int> & index,
                       const unsigned & col = 0);

  void addtocurrentcol(const datamatrix & m,const unsigned & col=0);

  void addtocurrentcol(const double & m,const unsigned & col=0);

  void addtocurrentcol_single(const double & m,const unsigned & r,
                       const unsigned & col=0);

  // FUNCTION: add_linearpred
  // TASK: adds the value m to the element of the i-th row and j-th column of
  //       the linear predictor
  //       if current = true the value(s) of the current linear predictor
  //       will be returned, else the value(s) of the proposed
  // IMPORTANT: FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  void add_linearpred(const double & m,const unsigned & i,const unsigned & j,
                      const bool & current = true)
    {
    if (current)
      (*linpred_current)(i,j)+= m;
    else
      (*linpred_proposed)(i,j)+= m;
    }

  void add_linearpred2(const double & m,const int index,
                      const bool & current = true)
    {
    if (current)
      {
      linpredp_current += index*linearpred.cols();
      *linpredp_current += m;
      }
    else
      {
      linpredp_proposed += index*linearpred.cols();
      *linpredp_proposed += m;
      }

    }


  // FUNCTION: add_linearpred
  // TASK: adds the value 'm' between 'beg' and 'end' to the linear predictor
  //       if current = true the value(s) of the current linear predictor
  //       will be returned, else the value(s) of the proposed
  // IMPORTANT: FOR UNIVARIATE AND MULTIVARIATE RESPONSE

  void add_linearpred(const double & m,const unsigned & beg,
                      const unsigned & end,const statmatrix<int> & index,
                      const unsigned & col,
                      const bool & current = true);

  void add_linearpred2(const double & m,const unsigned & beg,
                       const unsigned & end,const statmatrix<int> & index2,
                       const unsigned & col,
                       const bool & current = true);

  // FUNCTION: substr_linearpred
  // TASK: substracts the value m to the element of the i-th row and 0-th column
  //       of the linear predictor
  //       if current = true the value(s) of the current linear predictor
  //       will be returned, else the value(s) of the proposed
  // IMPORTANT: for univariate response only


  void substr_linearpred(const double & m,const unsigned & i,
                         const bool & current = true)
    {
    if (current)
      (*linpred_current)(i,0)-= m;
    else
      (*linpred_proposed)(i,0)-= m;
    }

  // FUNCTION: swap_linearpred
  // TASK: swaps the pointers linpred_current and linpred_proposed, that is
  //       the proposed new predictor will be the current

  void swap_linearpred(void);

  virtual void compute_respminuslinpred(datamatrix & res,
                                        const unsigned & co);

  //----------------------------------------------------------------------------
  //----------------- FUNCTIONS FOR COMPUTING THE POSTERIOR MODE ---------------
  //----------------------------------------------------------------------------

  // FUNCTION: compute_iwls
  // TASK: computes the working responses (stored in tildey) and the working
  //       weights for iwls (stored in weightiwls) based on the current
  //       linear predictor

  virtual void compute_iwls(void);

  // FUNCTION: fisher for iwls
  // TASK: computes  data' W data  and stores the result in 'XWX' , where
  //       'data' is the designmatrix,
  //       'w' = diag [ weight_i/scale*[b''(theta_i) g'^2(mu_i)]^{-1} ]
  //       (the function uses weightiwls as weights)
  // IMPORTANT: 'XWX', i.e. XWX must have dimension data.cols() x data.cols()

  void fisher(datamatrix & XWX,datamatrix & data,const unsigned & col=0) const;


  // FUNCTION: set_weightp
  // TASK: sets the pointer weightp to the first element in weight

  void set_workingresp(const unsigned & c=0)
    {
    workingresp = workingres.getV()+c;
    }


  inline const double & get_workingres(const int & i)
    {
    workingresp+=i;
    return *workingresp;
    }

  void compute_workingresiduals(const unsigned & col=0);

  void compute_weightiwls_workingresiduals(const unsigned & col=0);

  const datamatrix & get_workingresiduals(void)
    {
    return workingres;
    }

  const datamatrix & get_weightiwls(void)
    {
    return weightiwls;
    }

  const datamatrix & get_tildey(void)
    {
    return tildey;
    }

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  virtual bool posteriormode(void)
    {
    return true;
    }


  virtual bool posteriormode_converged(const unsigned & itnr)
    {
    return true;
    }

  virtual bool posteriormode_converged_fc(const datamatrix & beta,
                                          const datamatrix & beta_mode,
                                          const unsigned & itnr);

  void posteriormode_set_beta_mode(void);

  //----------------------------------------------------------------------------
  //--------------------------- UPDATE FUNCTIONS -------------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: update
  // TASK: base function for inherited classes,
  //       should update the scale parameter
  //       the base function updates the estimated mean and variance
  //       of the scale parameter only

  virtual void update(void);

  // FUNCTION: outresults
  // TASK: writes estimation results for the scale parameter
  //       estimated mean and variance

  virtual void outresults(void);

  // FUNCTION: reset
  // TASK: resets linpred (all values to 0) and scale (0.1)

  void reset(void);


  vector<ST::string> & get_results_latex(void)
    {
    return results_latex;
    }

  //----------------------------------------------------------------------------
  //------------------ Transformations of parameters  --------------------------
  //----------------------------------------------------------------------------

  void transform_nonlinear(vector<FULLCOND *> & fc,ST::string & trtype);

  virtual void tr_nonlinear(vector<double *> b,vector<double *> br,
                            vector<FULLCOND*> & fcp,unsigned & nr,
                            unsigned & it,ST::string & trtype)
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      *(br[i]) = exp(*(b[i]));
    }

  //----------------------------------------------------------------------------
  //----------------------------- HELP VARIABLES -------------------------------
  //----------------------------------------------------------------------------

  datamatrix helpmat1;            // help matrix that may be used for
                                  // calculations
                                  // size: nrobs times 1

  //----------------------------------------------------------------------------
  //-------------------------- FOR MISSING VALUES ------------------------------
  //----------------------------------------------------------------------------

  virtual void DISTRIBUTION::update_missings(void);

  void set_missings(vector<FULLCOND *> & fcm,unsigned & begin,
                    unsigned & end,datamatrix & mi,
                    ST::string & pt,ST::string & pr);

  //----------------------------------------------------------------------------
  //-------------------------- FOR LEYRE--------- ------------------------------
  //----------------------------------------------------------------------------

  void exchange_intercept(double & interconst)
    {
    interconst=interceptold;
    }

  void exchange_intercept2(double & inter)
    {
    interceptold = inter;
    }


  }; // end: class DISTRIBUTION



//------------------------------------------------------------------------------
//--------------- A CLASS FRAGMENT FOR A NEW DISTRIBUTION OBJECT----------------
//------------------------------------------------------------------------------

/*

// full conditional classes using conditional prior proposals work properly
// provided that the following virtual functions are implemented:

//  double loglikelihood(double * response,double * linpred,
//                       double * weight,const int & i) const


// full conditional class for fixed effects (non Gaussian case) work properly
// provided that the following virtual functions are implemented:

//  double compute_weight(double * linpred,double * weight,
//                        const unsigned & col=0) const

//  double compute_gmu(double * linpred) const


// predicted linearpred, deviance, mu etc. will be computed correctly if the
// following functions are implemented:

//  void compute_mu(const double * linpred,double * mu) const

//  void compute_devresidual(const double * response,const double * weight,
//                           const double * mu,
//                           double * residual, const datamatrix & scale,const int & i) const


class newdistribution : public DISTRIBUTION
  {

   protected:

   // include here additional private variables needed for the new
   // distribution

   public:

   //  the following public functions must be implemented


   // DEFAULT CONSTRUCTOR

   newdistribution(void) : DISTRIBUTION()
     {

     // include program code

     }

   // CONSTRUCTOR

   newdistribution(// include here additional variables
                   MCMCoptions * o, const datamatrix & r,
                   const datamatrix & w=datamatrix())

   : DISTRIBUTION(o,r,w)
     {
     // include program code
     // important: string variable 'family' must be declared
     // e.g. family = gaussian;
     // set scaleexisting = false, if the scale parameter is constant
     // e.g. binomial, poisson
     }

   // COPY CONSTRUCTOR

   newdistribution(const newdistribution & nd) : DISTRIBUTION(DISTRIBUTION(nd))
     {

     // include additional code
     // additional variables specified in the protected part must be
     // assigned, e.g. X = dn.X

     }

   // OVERLOADED ASSIGNMENT OPERATOR

   const newdistribution & operator=(const newdistribution & nd)
     {

     if (this==&nd)
	   return *this;
     DISTRIBUTION::operator=(DISTRIBUTION(nd));

     // include additional code
     // additional variables specified in the protected part must be
     // assigned, e.g. X = dn.X

     return *this;
     }

   // DESTRUCTOR

   ~newdistribution() {}


   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * response,double * linpred,
                       double * weight,const int & i) const
    {

    // include program code to compute the loglikelihood for the distribution

    }

  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred' and stores
  //       the result in 'mu'

  void compute_mu(const double * linpred,double * mu) const
     {

    // include program code

     }

  // FUNCTION: compute_devresidual
  // TASK: computes the deviance residual

  void compute_devresidual(const double * response,const double * weight,
                           const double * mu, double * residual, const datamatrix & scale,const int & i) const
    {

    // include program code

    }


   // FUNCTION: compute_weight
   // TASK: computes the weights for iteratively weighted least squares
   //       w_i = weight_i/scale*[b''(theta_i) g'^2(mu_i)]^{-1}

  double compute_weight(double * linpred,double * weight,
                        const int & i, const unsigned & col=0) const
    {

    // include program code

    }


  // FUNCTION: compute_gmu
  // TASK: compute g'(eta_i) = 1/h'(eta_i)

  double compute_gmu(double * linpred) const
    {

    // include program code

    }

  // FUNCTION: compute_g
  // TASK: compute g(mu_i)

  double compute_g(double * mu) const
    {

    // include program code

    }


  // FUNCTION: outoptions
  // TASK: writing options of the distribution

  void outoptions(void)
    {

    // include program code here

    DISTRIBUTION::outoptions();
    }


  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void)
    {

    // update program code to update the intercept

    // include program code to update the scale parameter (e.g. sigma2 for
    // gaussian response. Leave this free, if the special distribution has
    // a constant scale parameter (e.g. poisson)

    DISTRIBUTION::update(); // this command should be inserted, because
                            // the update method of the base class
                            // updates the estimaed mean and the variance of
                            // the scale paramter
                            // updates also intercept samples
    }

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode for the scale parameter

  bool posteriormode(void)
    {

    // include program code

    }



  };
*/


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gamma -------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// MODELLIERUNG DES SKALENPARAMETERS SCALE=PHI
//
// Implementiert sind drei Möglichkeiten:
//
//
// 1. FIXER SKALENPARAMETER (Konstruktor 0)
//
// Phi wird im Konstruktor 0 einem bekannten, über alle Iterationen konstant
// bleibenden Wert gleichgesetzt.
//
//
// 2. KONSISTENTE SCHÄTZUNG DES SKALENPARAMETERS (Konstruktor 1)
//
// Phi wird in jeder Iteration der Simulation mit phi_hut
// (vgl. Fahrmeir/Tutz S.44) konsistent geschätzt.
//
// Zur Verbesserung und Beschleunigung von Konvergenz und Mixing des
// Gesamtmodells kann phi zu Beginn für eine bestimmte Anzahl von Iterationen
// (günstig ca. 500) konstant gehalten werden, siehe Konstruktor 1.
//
//
// 3. UPDATE DES SKALENPARAMETERS MIT MH-ALGORITHMUS (Konstruktor 2)
//
// Der FORMPARAMETER nu = 1/phi wird in die MCMC-Simulation mitaufgenommen und
// in jeder Iteration mit dem MH-Algorithmus upgedatet.
// Der Skalenparameter phi ergibt sich dann jeweils zu phi = 1/nu.
//
// PRIORI: Gamma(a,b)
// PROPOSAL: Gamma(a_nu,b_nu) mit fester Varianz var = (a_nu/b_nu*b_nu)
// Die feste Varianz muß als Tuning Parameter im Konstruktor 2 geeignet
// gewählt werden. Dabei gilt: Für große Varianzen ist die Akzeptanz des
// Formparameters gering u.u.
// FULL CONDITIONAL UND AKZEPTANZWAHRSCHEINLICHKEIT:
// siehe Diplomarbeit Petra Kragler S.93 und S.101ff
//
// Zur Verbesserung und Beschleunigung von Konvergenz und Mixing des
// Gesamtmodells beginnt das Update des Formparameters erst beim Burnin
// der Gesamtsimulation. Bis dahin wird der Skalenparameter wie in 2. in jeder
// Iteration konsistent geschätzt. (Dazu ist ebenfalls wieder ein Konstanthalten
// in den allerersten Iterationen möglich.)
//------------------------------------------------------------------------------


class DISTRIBUTION_gamma : public DISTRIBUTION
  {

  protected:

  DISTRIBUTION * dgaussian;

  datamatrix lgamma;


  double a_gamma;                   // hyperparameter a (for the (priori)
                                     // gamma distribution of the shape
                                     // parameter, i.e. nu ,
                                     // nu = 1/phi = 1/scale)
  double b_gamma;                   // hyperparameter b

  double scaleold;

  bool scalefixed;

  unsigned nriterations;
  unsigned acceptance;

  bool mh;                          // estimation of scale parameter with
                                     // Metropolis-Hastings or not

  unsigned const_it;                // number of "burnin" iterations with
                                     // constant scale parameter

//   with MH-algorithm for scale:

  double var_nu;                     // fixed variance of shape parameter nu
                                     //(used for proposal after burnin)


  bool variance_gaussian;

  void check(void);

  void standardize(void);

  public:


  // DEFAULT CONSTRUCTOR

   DISTRIBUTION_gamma(void) : DISTRIBUTION()
     {
     family = "gamma";
     a_gamma = 1;
     b_gamma = 0.005;
     scale(0,0) = 0.1;
     const_it = 500;
     }


   // CONSTRUCTOR 0 (with constant scale)
   // scale = scale_initial

   DISTRIBUTION_gamma(const double & scale_initial,
                         MCMCoptions * o,
                         const datamatrix & r,
                         const ST::string & p,
                         const ST::string & ps,
                         const datamatrix & w=datamatrix());

   // CONSTRUCTOR 1 (with consistent estimation for scale)
   // a_gamma = a
   // b_gamma = b
   // const_it = cit;

   DISTRIBUTION_gamma(const double & a,
                         const double & b,
                         const unsigned & cit,
                         MCMCoptions * o,
                         const datamatrix & r,
                         const ST::string & p,
                         const ST::string & ps,
                         const datamatrix & w=datamatrix());


   // CONSTRUCTOR 2  (with MH-algorithm for scale)
   // a_gamma = a
   // b_gamma = b
   // var_nu = var;
   // const_it = cit;

   DISTRIBUTION_gamma(const double & a,
                         const double & b,
                         const double & var,
                         const unsigned & cit,
                         MCMCoptions * o,
                         const datamatrix & r,
                         const ST::string & p,
                         const ST::string & ps,
                         const datamatrix & w=datamatrix());


   // COPY CONSTRUCTOR

   DISTRIBUTION_gamma(const DISTRIBUTION_gamma & ga);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_gamma & operator=(const DISTRIBUTION_gamma & ga);

   // DESTRUCTOR

   ~DISTRIBUTION_gamma() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * res, double * lin, double * w,
                       const int & i) const;

  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred'

  void compute_mu(const double * linpred,double * mu) const;

  // FUNCTION: compute_devresidual
  // TASK: computes the deviance residual
  // weight NICHT berücksichtigt

  void compute_deviance(const double * response,const double * weight,
                        const double * mu,double * deviance,
                        double * deviancesat,const datamatrix & scale,
                        const int & i) const;

  // FUNCTION: update
  // TASK: updates the scale parameter and the interecept

  void update(void);

  void update_predict(void);

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void outoptions(void);

  void outresults(void);

   // FUNCTION: compute_weight
   // TASK: computes the weights for iteratively weighted least squares

  double compute_weight(double * worklin,double * weight,
                        const int & i,const unsigned & col=0) const;

  double compute_gmu(double * linpred,const unsigned & col=0) const
    {
    return 1.0/exp(*linpred);
    }


  void create_lgamma(void);

  // FUNCTION lgammafunc
  // TASK: computes the log gammafunction for a special value nu

  double lgammafunc(const double & nu) const;

  // FUNCTION lfac
  // TASK: computes the log faculty for a special value nu

  double lfac(const double & nu) const;

  // FUNCTION log_prop
  // TASK: computes the log proposal for a special value nu and hyperparameters a, b

  double log_prop(const double & nu, double & a, double & b) const;

  // FUNCTION phi_hat
  // TASK: computes phi_hat as a consistent estimation for the scale parameter phi

  double phi_hat() const;

  // FUNCTION log_fullcond
  // TASK: computes the log full conditional for a particular value nu

  double log_fullcond(const double & nu, double & bnew,double & sw) const;

  void update_variance(datamatrix & we);
//  void update_variance(const double & in,datamatrix & res,
//                       datamatrix & linp, datamatrix & we);

  void set_gaussian(DISTRIBUTION * g)
    {
    dgaussian = g;
    variance_gaussian=true;
    }


  void tr_nonlinear(vector<double *> b,vector<double *> br,
                    vector<FULLCOND*> & fcp,unsigned & nr,
                    unsigned & it,ST::string & trtype);

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_gaussian ----------------------------
//------------------------------------------------------------------------------

class DISTRIBUTION_gaussian : public DISTRIBUTION
  {

  protected:

  double a_invgamma;                   // hyperparameter a (for the inverse
                                        // gamma distribution of the scale
                                        // parameter, i.e. sigma^2
  double b_invgamma;                   // hyperparameter b

  DISTRIBUTION_gamma * dgamma;
  bool varianceest;

  bool constscale;

  // FUNCTION: standardise
  // TASK: standardises the response and the offset
  //       sets scalesave.transform = trmult*trmult (!!!)

  void standardise(void);


  public:


   // DEFAULT CONSTRUCTOR

   DISTRIBUTION_gaussian(void) : DISTRIBUTION()
     {
     family = "gaussian";
     a_invgamma = 1;
     b_invgamma = 0.005;
     varianceest = false;
     constscale = false;
     }

   // CONSTRUCTOR1
   // a_invgamma = a
   // b_invgamma = b

   DISTRIBUTION_gaussian(const double & a,
                         const double & b,
                         MCMCoptions * o,
                         const datamatrix & r,const ST::string & p,
                         const ST::string & ps,
                         const datamatrix & w=datamatrix());

   // CONSTRUCTOR2
   // a_invgamma = a
   // b_invgamma = b

   DISTRIBUTION_gaussian(const datamatrix & offset, const double & a,
                         const double & b,
                         MCMCoptions * o,
                         const datamatrix & r,
                         const ST::string & p,const ST::string & ps,
                         const datamatrix & w=datamatrix());



   // COPY CONSTRUCTOR

   DISTRIBUTION_gaussian(const DISTRIBUTION_gaussian & nd)
   : DISTRIBUTION(DISTRIBUTION(nd))
     {
     a_invgamma = nd.a_invgamma;
     b_invgamma = nd.b_invgamma;
     varianceest = nd.varianceest;
     constscale = nd.constscale;
     }

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_gaussian & operator=(const DISTRIBUTION_gaussian & nd);

   // DESTRUCTOR

   ~DISTRIBUTION_gaussian() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * res,
                       double * lin,
                       double * w,
                       const int & i) const;


  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred'

  void compute_mu(const double * linpred,double * mu) const;

  // FUNCTION: compute_deviance
  // TASK: computes the retransformed individual deviance
  //       scale and response is assumed to be NOT RETRANSFORMED
  //       but will be retransformed when computing the residual
  //       mu is assumed to be already restransformed

  void compute_deviance(const double * response,const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const;

  double compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col=0) const;

  void compute_iwls(void)
    {
    tildey.assign(response);
    DISTRIBUTION::compute_weight(weightiwls,0);
    }

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  void update_predict(void)
    {
    DISTRIBUTION::update_predict();
    }

  bool posteriormode(void);

  bool posteriormode_converged_fc(const datamatrix & beta,
                                  const datamatrix & beta_mode,
                                  const unsigned & itnr)
    {
    return true;
    }

  void outresults(void)
    {
    DISTRIBUTION::outresults();
    }


  void set_variance(DISTRIBUTION_gamma * dg);

  void get_residuals(datamatrix & r);

  void tr_nonlinear(vector<double *> b,vector<double *> br,
                    vector<FULLCOND*> & fcp,unsigned & nr,
                    unsigned & it,ST::string & trtype);

  void set_constscale(double s);

  void update_missings(void);

  double compute_gcv(double & df);

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_binomial ----------------------------
//------------------------------------------------------------------------------


class DISTRIBUTION_binomial : public DISTRIBUTION
  {

   protected:

   void create(void);

   public:


   // DEFAULT CONSTRUCTOR

   DISTRIBUTION_binomial(void) : DISTRIBUTION()
     {
     }

   // CONSTRUCTOR1

   DISTRIBUTION_binomial(MCMCoptions * o,
                         const datamatrix & r,
                         const datamatrix & w=datamatrix());

   // CONSTRUCTOR2

   DISTRIBUTION_binomial(const datamatrix & offset,MCMCoptions * o,
                         const datamatrix & r,
                         const datamatrix & w=datamatrix());


   // COPY CONSTRUCTOR

   DISTRIBUTION_binomial(const DISTRIBUTION_binomial & nd)
   : DISTRIBUTION(DISTRIBUTION(nd))
     {
     }

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_binomial & operator=(const DISTRIBUTION_binomial & nd)
     {
     if (this==&nd)
	   return *this;
     DISTRIBUTION::operator=(DISTRIBUTION(nd));

     return *this;
     }

   // DESTRUCTOR

   ~DISTRIBUTION_binomial() {}


   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * response,
                       double * linpred,
                       double * weight,
                       const int & i) const;


  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred'

  void compute_mu(const double * linpred,double * mu) const;

  // FUNCTION: compute_deviance
  // TASK: computes the deviance

  void compute_deviance(const double * response,const double * weight,
                           const double * mu,double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const;

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  void update_predict(void)
    {
    DISTRIBUTION::update_predict();
    }

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void outresults(void)
    {
    DISTRIBUTION::outresults();
    }

  double compute_weight(double * linpred,double * weight,
                        const int & i,const unsigned & col) const;

  double compute_gmu(double * linpred,const unsigned & col=0) const;


  void tr_nonlinear(vector<double *> b,vector<double *> br,
                    vector<FULLCOND*> & fcp,unsigned & nr,
                    unsigned & it,ST::string & trtype);


  };



//------------------------------------------------------------------------------
//------------------- CLASS: DISTRIBUTION_binomial_latent ----------------------
//------------------------------------------------------------------------------


class DISTRIBUTION_binomial_latent : public DISTRIBUTION
  {

   protected:

   datamatrix res;

   double nu;
   bool tlink;                             // response function is t-cdf with
                                           // nu degrees of freedom

   void create(const bool & tl,const unsigned & n);

   public:


   // DEFAULT CONSTRUCTOR

   DISTRIBUTION_binomial_latent(void) : DISTRIBUTION()
     {
     }

   // CONSTRUCTOR1

   DISTRIBUTION_binomial_latent(MCMCoptions * o,const datamatrix & r,
                                const datamatrix & w,
                                const bool & tl,
                                const unsigned & n = 8);

   // CONSTRUCTOR2

   DISTRIBUTION_binomial_latent(const datamatrix & offset,MCMCoptions * o,
                                const datamatrix & r,
                                const datamatrix & w,
                                const bool & tl, const unsigned & n = 8);


   // COPY CONSTRUCTOR

   DISTRIBUTION_binomial_latent(const DISTRIBUTION_binomial_latent & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_binomial_latent &
   operator=(const DISTRIBUTION_binomial_latent & nd);

   // DESTRUCTOR

   ~DISTRIBUTION_binomial_latent() {}


   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * response,
                       double * linpred,
                       double * weight,
                       const int & i) const;


  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred'

  void compute_mu(const double * linpred,double * mu) const;

  // FUNCTION: compute_deviance
  // TASK: computes the indiviudal deviance

  void compute_deviance(const double * response,const double * weight,
                           const double * mu,double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const;

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the latent variables

  void update(void);

  void update_predict(void)
    {
    DISTRIBUTION::update_predict();
    }

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void outresults(void);


  double compute_weight(double * linpred,double * weight,
                        const int & i,const unsigned & col) const;

  double compute_gmu(double * linpred,const unsigned & col=0) const;


  void tr_nonlinear(vector<double *> b,vector<double *> br,
                    vector<FULLCOND*> & fcp,unsigned & nr,
                    unsigned & it,ST::string & trtype);

  };


//------------------------------------------------------------------------------
//--------------- CLASS: DISTRIBUTION_binomial_logit_latent --------------------
//------------------------------------------------------------------------------

class DISTRIBUTION_binomial_logit_latent : public DISTRIBUTION
  {

   protected:

   datamatrix res;
   datamatrix as;
   unsigned acceptancescale;
   unsigned nrtrials;


   void create(const bool & tl);

   public:


   // DEFAULT CONSTRUCTOR

   DISTRIBUTION_binomial_logit_latent(void) : DISTRIBUTION()
     {
     }

   // CONSTRUCTOR1

   DISTRIBUTION_binomial_logit_latent(MCMCoptions * o,const datamatrix & r,
                                const datamatrix & w,
                                const bool & tl);

   // CONSTRUCTOR2

   DISTRIBUTION_binomial_logit_latent(const datamatrix & offset,MCMCoptions * o,
                                const datamatrix & r,
                                const datamatrix & w,
                                const bool & tl);


   // COPY CONSTRUCTOR

   DISTRIBUTION_binomial_logit_latent(
   const DISTRIBUTION_binomial_logit_latent & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_binomial_logit_latent &
   operator=(const DISTRIBUTION_binomial_logit_latent & nd);

   // DESTRUCTOR

   ~DISTRIBUTION_binomial_logit_latent() {}


   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * response,
                       double * linpred,
                       double * weight,
                       const int & i) const;


  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred'

  void compute_mu(const double * linpred,double * mu) const;

  // FUNCTION: compute_deviance
  // TASK: computes the indiviudal deviance

  void compute_deviance(const double * response,const double * weight,
                           const double * mu,double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const;

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the latent variables

  void update(void);

  void update_predict(void)
    {
    DISTRIBUTION::update_predict();
    }

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void outresults(void);


  double compute_weight(double * linpred,double * weight,
                        const int & i,const unsigned & col) const;

  double compute_gmu(double * linpred,const unsigned & col=0) const;


  void tr_nonlinear(vector<double *> b,vector<double *> br,
                    vector<FULLCOND*> & fcp,unsigned & nr,
                    unsigned & it,ST::string & trtype);

  };


//------------------------------------------------------------------------------
//------------------------ CLASS: DISTRIBUTION_poisson -------------------------
//------------------------------------------------------------------------------

class DISTRIBUTION_poisson : public DISTRIBUTION
  {

   protected:

   public:

   // DEFAULT CONSTRUCTOR

   DISTRIBUTION_poisson(void) : DISTRIBUTION()
     {
     }

   // CONSTRUCTOR1 (without offset)

   DISTRIBUTION_poisson(MCMCoptions * o, const datamatrix & r,
                        const datamatrix & w=datamatrix());

   // CONSTRUCTOR2

   DISTRIBUTION_poisson(const datamatrix & offset, MCMCoptions * o,
                        const datamatrix & r,
                        const datamatrix & w=datamatrix());


   // COPY CONSTRUCTOR

   DISTRIBUTION_poisson(const DISTRIBUTION_poisson & nd)
   : DISTRIBUTION(DISTRIBUTION(nd))
     {
     }

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_poisson & operator=(const DISTRIBUTION_poisson & nd)
     {

     if (this==&nd)
	   return *this;

     DISTRIBUTION::operator=(DISTRIBUTION(nd));

     return *this;
     }

   // DESTRUCTOR

   ~DISTRIBUTION_poisson() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * response,double * linpred,
                       double * weight,const int & i) const
    {
    return *weight * (*response * *linpred - exp(*linpred));
    }

  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred'

  void compute_mu(const double * linpred,double * mu) const
    {
    *mu = exp(*linpred);
    }


  // FUNCTION: compute_deviance
  // TASK: computes the individual deviance

  void compute_deviance(const double * response,const double * weight,
                           const double * mu,double * deviance,double *
                           deviancesat,
                           const datamatrix & scale,const int & i) const;

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  void update_predict(void)
    {
    DISTRIBUTION::update_predict();
    }

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void outresults(void)
    {
    DISTRIBUTION::outresults();
    }

  double compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col) const;

  double compute_gmu(double * linpred,const unsigned & col=0) const
    {
    return 1.0/exp(*linpred);
    }


  void tr_nonlinear(vector<double *> b,vector<double *> br,
                    vector<FULLCOND*> & fcp,unsigned & nr,
                    unsigned & it,ST::string & trtype);


  };  // end: class DISTRIBUTION_poisson


//------------------------------------------------------------------------------
//-------------------------- CLASS: multinomial --------------------------------
//------------------------------------------------------------------------------

class DISTRIBUTION_multinom : public DISTRIBUTION
  {

   protected:

   ST::string reference;

   datamatrix muhelp;

   public:

   // DEFAULT CONSTRUCTOR

   DISTRIBUTION_multinom(void) : DISTRIBUTION()
     {
     }

   // CONSTRUCTOR

   DISTRIBUTION_multinom(MCMCoptions * o, const datamatrix & r,
                         const double & refvalue,
                         const datamatrix & w=datamatrix());


   // COPY CONSTRUCTOR

   DISTRIBUTION_multinom(const DISTRIBUTION_multinom & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_multinom & operator=(const DISTRIBUTION_multinom & nd);

   // DESTRUCTOR

   ~DISTRIBUTION_multinom() {}

  void outoptions(void);

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * resp,double * linpred,
                       double * weight,const int & i) const;

  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred' and returns the
  //        result

  void compute_mu(const double * linpred,double * mu) const;

  void compute_deviance(const double * response,const double * weight,
                           const double * mu,double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const;

  // FUNCTION: compute_weight
  // TASK: computes the weights for iteratively weighted least squares and
  //       stores the result in w, w must be of size nrobs X 1
  //       see e.g. Fahrmeir, Tutz (1997) page 39 ff.

  double compute_weight(double * linpred,double * weight,
                        const int & i,const unsigned & col=0) const;

  double compute_gmu(double * linpred,const unsigned & col=0) const;

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

  void compute_iwls(void);

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void update_predict(void)
    {
    DISTRIBUTION::update_predict();
    }


  };



//------------------------------------------------------------------------------
//------------------ CLASS: DISTRIBUTION_multinomial_latent --------------------
//------------------------------------------------------------------------------



class DISTRIBUTION_multinomial_latent : public DISTRIBUTION
  {

  protected:

  datamatrix responsecat;

  vector<unsigned> posbeg;
  vector<unsigned> posend;

  double refvalue;
  unsigned refcat;

  unsigned nrcat;

  double maxutility(double * r,const unsigned & cat);

  public:


   // DEFAULT CONSTRUCTOR

   DISTRIBUTION_multinomial_latent(void) : DISTRIBUTION()
     {
     }

   // CONSTRUCTOR

   DISTRIBUTION_multinomial_latent(MCMCoptions * o, const datamatrix & r,
                                   const double & rv);

   // COPY CONSTRUCTOR

   DISTRIBUTION_multinomial_latent(const DISTRIBUTION_multinomial_latent & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_multinomial_latent &
   operator=(const DISTRIBUTION_multinomial_latent & nd);

   // DESTRUCTOR

   ~DISTRIBUTION_multinomial_latent() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * response,
                       double * linpred,
                       double * weight,
                       const int & i) const;

  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred'

  void compute_mu(const double * linpred,double * mu) const;

  // FUNCTION: compute_deviance
  // TASK: computes the individual deviance

  void compute_deviance(const double * response,const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const;

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the latent variables

  void update(void);

  void update_predict(void)
    {
    DISTRIBUTION::update_predict();
    }

  bool posteriormode(void)
    {
    return true;
    }

  void outresults(void);

  double compute_weight(double * linpred,double * weight,
                        const int & i,const unsigned & col) const;

  double compute_gmu(double * linpred,const unsigned & col=0) const;


  unsigned get_nrcat(void)
    {
    return nrcat;
    }


  };


//------------------------------------------------------------------------------
//------------------ CLASS: DISTRIBUTION_cumulative_latent3 --------------------
//------------------------------------------------------------------------------

// Important: values of the response are assumed to be sorted !!!

class DISTRIBUTION_cumulative_latent3 : public DISTRIBUTION
  {

  protected:

  vector<unsigned> posbeg;             // vector whose i-th element indicates
                                       // the first position of the
                                       // i-th category in 'response'
  vector<unsigned> posend;             // vector whose i-th element indicates
                                       // the last position of the
                                       // i-th category in 'response'

  double refvalue;                     // reference category
                                       // (always the largest value of the
                                       //  response)

  double a_invgamma;
  double b_invgamma;

  double sumweight;

  void update_utilities(void);

  public:

   // DEFAULT CONSTRUCTOR

   DISTRIBUTION_cumulative_latent3(void) : DISTRIBUTION()
     {
     }

   // CONSTRUCTOR

   DISTRIBUTION_cumulative_latent3(MCMCoptions * o, const datamatrix & r,
                                  const datamatrix & w,
                                  const double & a,const double & b,
                                  const ST::string & p,const ST::string & ps);

   // COPY CONSTRUCTOR

   DISTRIBUTION_cumulative_latent3(const DISTRIBUTION_cumulative_latent3 & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_cumulative_latent3 &
   operator=(const DISTRIBUTION_cumulative_latent3 & nd);

  // FUNCTION: set_predict
  // TASK: indicates that the mean of the linearpredictor should be updated
  //       in every iteration, mean linearpred, deviance, etc. are written
  //       to 'path'

  void set_predict_cum(const ST::string & path,
                  const ST::string & pathdev,datamatrix * p,
                   vector<ST::string> & Dn);

   // DESTRUCTOR

   ~DISTRIBUTION_cumulative_latent3() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation

  double loglikelihood(double * response,
                       double * linpred,
                       double * weight,
                       const int & i) const;

  // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred'

  void compute_mu(const double * linpred,double * mu) const;

  // FUNCTION: compute_devresidual
  // TASK: computes the deviance residual

  void compute_deviance(const double * response,const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const;

  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the latent variables
  //        and the intercept

  void update(void);

  void update_predict(void);

  bool posteriormode(void)
    {
    return true;
    }

  void outresults(void);

  };


} // end: namespace MCMC


#endif
