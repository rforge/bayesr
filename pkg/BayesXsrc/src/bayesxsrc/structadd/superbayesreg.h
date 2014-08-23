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



#if !defined (superBAYESREG_INCLUDED)

#define superBAYESREG_INCLUDED

#include"../export_type.h"
#include"statobj.h"
#include"dataobj.h"
#include"MASTER_obj.h"
#include"GENERAL_OPTIONS.h"

#include"model_parameters.h"

#include"distr.h"
#include"distr_categorical.h"
#include"distr_categorical_mult.h"
#include"distr_mixture.h"
#include"distr_gamlss.h"
#include"distr_zeroadjusted.h"
#include"distr_gamlss_nadja.h"

#include"design.h"
#include"design_pspline.h"
#include"design_hrandom.h"
#include"design_mrf.h"
#include"design_kriging.h"

#include"FC.h"
#include"FC_predict.h"
#include"FC_predict_mult.h"
#include"FC_predict_predictor.h"
#include"FC_predictive_check.h"
#include"FC_nonp.h"
#include"FC_linear.h"
#include"FC_hrandom.h"
#include"FC_mult.h"
#include"FC_nonp_variance.h"
#include"FC_nonp_variance_vec.h"
#include"FC_variance_pen_vector.h"
#include"FC_hrandom_variance.h"
#include"FC_hrandom_variance_vec.h"
#include"FC_hrandom_variance_vec_nmig.h"
#include"FC_cv.h"

#include"mcmcsim.h"

using MCMC::MASTER_OBJ;

using randnumbers::uniform;
using randnumbers::rand_normvek;
using randnumbers::rand_normal;
using randnumbers::rand_invgamma;
using MCMC::GENERAL_OPTIONS;

using MCMC::DISTR;
using MCMC::DISTR_gaussian;
using MCMC::DISTR_vargaussian;
using MCMC::DISTR_hetgaussian;
using MCMC::DISTR_gaussianmixture;
using MCMC::DISTR_quantreg;
using MCMC::DISTR_loggaussian;
using MCMC::DISTR_gaussian_re;
using MCMC::DISTR_gaussian_exp;
using MCMC::DISTR_gaussian_mult;
using MCMC::DISTR_binomial;
using MCMC::DISTR_poisson;
using MCMC::DISTR_poisson_ext;
using MCMC::DISTR_poisson_extlin;
using MCMC::DISTR_binomialprobit;
using MCMC::DISTR_cloglog;
using MCMC::DISTR_binomialsvm;
using MCMC::DISTR_logit_fruehwirth;
using MCMC::DISTR_multinomprobit;
using MCMC::DISTR_multgaussian;
using MCMC::DISTR_multinomlogit;
using MCMC::DISTR_ziplambda;
using MCMC::DISTR_zippi;
using MCMC::DISTR_hurdle_lambda;
using MCMC::DISTR_hurdle_pi;
using MCMC::DISTR_hurdle_mu;
using MCMC::DISTR_hurdle_delta;
using MCMC::DISTR_negbinzip_mu;
using MCMC::DISTR_negbinzip_pi;
using MCMC::DISTR_negbinzip_delta;
using MCMC::DISTR_zip_cloglog_mu;
using MCMC::DISTR_zip_cloglog_pi;
using MCMC::DISTR_negbin_mu;
using MCMC::DISTR_negbin_delta;
using MCMC::DISTR_beta_mu;
using MCMC::DISTR_beta_sigma2;
using MCMC::DISTR_lognormal_mu;
using MCMC::DISTR_lognormal_sigma2;
using MCMC::DISTR_lognormal2_mu;
using MCMC::DISTR_lognormal2_sigma;
using MCMC::DISTR_normal_mu;
using MCMC::DISTR_normal_sigma2;
using MCMC::DISTR_normal2_mu;
using MCMC::DISTR_normal2_sigma;
using MCMC::DISTR_truncnormal2_mu;
using MCMC::DISTR_truncnormal2_sigma;
using MCMC::DISTR_gamma_mu;
using MCMC::DISTR_gamma_sigma;
using MCMC::DISTR_pareto_b;
using MCMC::DISTR_pareto_p;
using MCMC::DISTR_invgaussian_mu;
using MCMC::DISTR_invgaussian_sigma2;
using MCMC::DISTR_gengamma_mu;
using MCMC::DISTR_gengamma_sigma;
using MCMC::DISTR_gengamma_tau;
using MCMC::DISTR_t_mu;
using MCMC::DISTR_t_sigma2;
using MCMC::DISTR_t_df;
using MCMC::DISTR_zeroadjusted;
using MCMC::DISTR_zeroadjusted_mult;
using MCMC::DISTR_weibull_lambda;
using MCMC::DISTR_weibull_alpha;
using MCMC::DISTR_dagum_a;
using MCMC::DISTR_dagum_b;
using MCMC::DISTR_dagum_p;
using MCMC::DISTR_betainf_mu;
using MCMC::DISTR_betainf_sigma2;
using MCMC::DISTR_betainf_nu;
using MCMC::DISTR_betainf0_nu;
using MCMC::DISTR_betainf1_tau;
using MCMC::DISTR_betainf_tau;
using MCMC::DISTR_bivt_mu;
using MCMC::DISTR_bivt_sigma;
using MCMC::DISTR_bivt_df;
using MCMC::DISTR_bivt_rho;
using MCMC::DISTR_bivnormal_mu;
using MCMC::DISTR_bivnormal_sigma;
using MCMC::DISTR_bivnormal_rho;
using MCMC::DISTR_bivnormal_mufz;
using MCMC::DISTR_bivnormal_rhofz;
using MCMC::DISTR_bivprobit_mu;
using MCMC::DISTR_bivprobit_rho;
using MCMC::DISTR_bivlogit_mu;
using MCMC::DISTR_bivlogit_or;
using MCMC::DISTR_dirichlet;
using MCMC::DISTR_BCCG_mu;
using MCMC::DISTR_BCCG_sigma;
using MCMC::DISTR_BCCG_nu;
using MCMC::DISTR_copula;
using MCMC::DISTR_gumbelcopula_rho;
using MCMC::DISTR_gumbelcopula2_rho;
using MCMC::DISTR_gumbelcopula2_normal_mu;
using MCMC::DISTR_gumbelcopula2_normal_sigma2;
using MCMC::DISTR_claytoncopula_rho;
using MCMC::DISTR_claytoncopula2_rho;
using MCMC::DISTR_claytoncopula2_normal_mu;
using MCMC::DISTR_claytoncopula2_normal_sigma2;
using MCMC::DISTR_sfa0_mu_y;
using MCMC::DISTR_sfa0_sigma_u;
using MCMC::DISTR_sfa0_sigma_v;
using MCMC::DISTR_sfa_mu_y_id;
using MCMC::DISTR_sfa_mu_u_id;
using MCMC::DISTR_sfa_mu_y;
using MCMC::DISTR_sfa_mu_u;
using MCMC::DISTR_sfa_sigma_u;
using MCMC::DISTR_sfa_sigma_v;
using MCMC::DISTR_sfa2_mu_y_id;
using MCMC::DISTR_sfa2_mu_u_id;
using MCMC::DISTR_sfa2_mu_y;
using MCMC::DISTR_sfa2_mu_u;
using MCMC::DISTR_sfa2_sigma_u;
using MCMC::DISTR_sfa2_sigma_v;
using MCMC::DISTR_sfa_alpha;
using MCMC::DISTR_copula;
using MCMC::DISTR_gaussiancopula_rho;
using MCMC::DISTR_gaussiancopula_rhofz;
using MCMC::DISTR_frankcopula_rho;
using MCMC::DISTR_frankcopula2_rho;
using MCMC::DISTR_frankcopula2_normal_mu;
using MCMC::DISTR_frankcopula2_normal_sigma2;
using MCMC::DISTR_tcopula_df;
using MCMC::DISTR_tcopula_rho;

using MCMC::DESIGN_pspline;
using MCMC::DESIGN_hrandom;
using MCMC::DESIGN_mrf;
using MCMC::DESIGN_kriging;
using MCMC::equation;

using MCMC::FC;
using MCMC::FC_predict;
using MCMC::FC_predict_mult;
using MCMC::FC_predict_predictor;
using MCMC::FC_predictive_check;
using MCMC::FC_cv;
using MCMC::FC_nonp;
using MCMC::FC_linear;
using MCMC::FC_linear_pen;
using MCMC::FC_mult;
using MCMC::FC_hrandom;
using MCMC::FC_nonp_variance;
using MCMC::FC_nonp_variance_varselection;
using MCMC::FC_hrandom_variance;
using MCMC::FC_variance_pen_vector;
using MCMC::FC_nonp_variance_vec;
using MCMC::FC_hrandom_variance_vec;
using MCMC::FC_hrandom_variance_vec_nmig;
using MCMC::FC_hrandom_variance_ssvs;


using MCMC::MCMCsim;


class __EXPORT_TYPE superbayesreg : public statobject
  {

  private :

  ST::string pathres;
  ST::string title;
  ST::string pathnonp;

  void make_header(unsigned & modnr);

  bool find_binomial(DISTR* & b);
  bool find_continuous_singleparam(DISTR* & m);
  bool find_continuous_multparam(vector<DISTR*> & m);

  void make_paths(ST::string & pathnonp, ST::string & pathres,
                  ST::string & title, vector<ST::string> vn,
                  ST::string  endingraw,
                  ST::string  endingres, ST::string  endingtitle);

  void extract_data(unsigned i, datamatrix & d,datamatrix & iv,unsigned dim_dm);

  bool findREdistr(ST::string & na,equation & maine,unsigned & fnr);

  bool check_errors(void);

  void clear(void);
//  void initpointers(void);

  //------------------------- PRIVATE VARIABLES --------------------------------

  // vector of pointers to current statobjects

  vector<statobject*> * statobj;

  // pointer to functions

  typedef void (* runpointer )(superbayesreg & b);

  runpointer functions[10];

  datamatrix D;

  vector<ST::string> modelvarnamesv;

  // global options

  fileoption outfile;

  optionlist globaloptions;

  // ---------------------  for method 'hregress'  ------------------------------

  vector<equation> equations;          // Vector of equations
  unsigned nrlevel1;
  MCMCsim simobj;                      // Simulation object;

  GENERAL_OPTIONS generaloptions;
  bool generaloptions_yes;
  bool create_generaloptions(void);

  bool errors;

  // OPTIONS for method hregress

  simpleoption modeonly;               // Computes the posterior mode only
  intoption setseed;
  intoption modemaxit;                 // maxmimum number of iterations for
                                      // posterior mode

  // general MCMC options
  intoption iterations;                // Number of iterations
  intoption burnin;                    // Number of burnin iterations
  intoption step;                      // Thinning parameter
  doubleoption level1;                 // Nominal level 1 of credible intervals
  doubleoption level2;                 // Nominal level 2 of credible intervals

  // general MCMC options

  vector<ST::string> families;          // Response families
  stroption family;                     // specifies the response distribution
  doubleoption aresp;                   // Hyperparameter a of overal variance
                                        // (Gaussian response)
  doubleoption bresp;                   // Hyperparameter b of overal variance
                                        // (Gaussian response)

  simpleoption standardize;             // standardize response before
                                        // estimation?

  intoption H;                          // Number of mixture components in
                                        // binomial_logit_l1

  intoption equationnr;                 // Equationnumber for multivariate
                                        // responses
  intoption hlevel;                     // hierarchy
  vector<ST::string> equationtypes;
  stroption equationtype;               // type of equation, e.g. mean, variance


  simpleoption includescale;            // multiply error variance to penalty
                                        // matrices

  simpleoption scaleconst;
  simpleoption utilities;

  simpleoption slow;

  // prediction
  vector<ST::string> predictop;
  stroption predict;

  simpleoption pred_check;

  simpleoption cv;

  simpleoption imeasures;

  vector<ST::string> MSEop;
  stroption mse;
  doubleoption mseparam;

  doubleoption quantile;

  // linear effects

  simpleoption centerlinear;

  optionlist regressoptions;

  // Zero Inflated Negative Binomial - delta parameter

  doubleoption stopsum;
  intoption stoprmax;
  doubleoption fraclimit;
  intoption nrbetween;

  // extended Poisson

  doubleoption aexp;
  doubleoption bexp;
  simpleoption adaptexp;

  // options for linpredlimits

  simpleoption changelinpredlimits;
  doubleoption linpredminlimit;
  doubleoption linpredmaxlimit;
  simpleoption saveestimation;

  //dirichlet regression

  intoption nrcat;

  // end: OPTIONS for method regress

 // ------------------------------- MASTER_OBJ ---------------------------------

 MASTER_OBJ master;


//---------------------------------- DISTR  ------------------------------------

  vector<DISTR_gaussian> distr_gaussians;
  vector<DISTR_hetgaussian> distr_hetgaussians;
  vector<DISTR_vargaussian> distr_vargaussians;
  vector<DISTR_quantreg> distr_quantregs;
  vector<DISTR_gaussianmixture> distr_gaussianmixtures;
  vector<DISTR_loggaussian> distr_loggaussians;
  vector<DISTR_gaussian_re> distr_gaussian_res;
  vector<DISTR_gaussian_exp> distr_gaussian_exps;
  vector<DISTR_gaussian_mult> distr_gaussian_mults;
  vector<DISTR_binomial> distr_binomials;
  vector<DISTR_poisson> distr_poissons;
  vector<DISTR_poisson_ext> distr_poisson_exts;
  vector<DISTR_poisson_extlin> distr_poisson_extlins;
  vector<DISTR_binomialprobit> distr_binomialprobits;
  vector<DISTR_cloglog> distr_cloglogs;
  vector<DISTR_binomialsvm> distr_binomialsvms;
  vector<DISTR_logit_fruehwirth> distr_logit_fruehwirths;
  vector<DISTR_multinomprobit> distr_multinomprobits;
  vector<DISTR_multgaussian> distr_multgaussians;
  vector<DISTR_multinomlogit> distr_multinomlogits;
  vector<DISTR_ziplambda> distr_ziplambdas;
  vector<DISTR_zippi> distr_zippis;
  vector<DISTR_hurdle_lambda> distr_hurdle_lambdas;
  vector<DISTR_hurdle_pi> distr_hurdle_pis;
  vector<DISTR_hurdle_mu> distr_hurdle_mus;
  vector<DISTR_hurdle_delta> distr_hurdle_deltas;
  vector<DISTR_negbinzip_mu> distr_negbinzip_mus;
  vector<DISTR_negbinzip_pi> distr_negbinzip_pis;
  vector<DISTR_negbinzip_delta> distr_negbinzip_deltas;
  vector<DISTR_zip_cloglog_pi> distr_zip_cloglog_pis;
  vector<DISTR_zip_cloglog_mu> distr_zip_cloglog_mus;
  vector<DISTR_negbin_mu> distr_negbin_mus;
  vector<DISTR_negbin_delta> distr_negbin_deltas;
  vector<DISTR_beta_mu> distr_beta_mus;
  vector<DISTR_beta_sigma2> distr_beta_sigma2s;
  vector<DISTR_lognormal_mu> distr_lognormal_mus;
  vector<DISTR_lognormal_sigma2> distr_lognormal_sigma2s;
  vector<DISTR_lognormal2_mu> distr_lognormal2_mus;
  vector<DISTR_lognormal2_sigma> distr_lognormal2_sigmas;
  vector<DISTR_normal_mu> distr_normal_mus;
  vector<DISTR_normal_sigma2> distr_normal_sigma2s;
  vector<DISTR_normal2_mu> distr_normal2_mus;
  vector<DISTR_normal2_sigma> distr_normal2_sigmas;
  vector<DISTR_truncnormal2_mu> distr_truncnormal2_mus;
  vector<DISTR_truncnormal2_sigma> distr_truncnormal2_sigmas;
  vector<DISTR_gamma_mu> distr_gamma_mus;
  vector<DISTR_gamma_sigma> distr_gamma_sigmas;
  vector<DISTR_pareto_b> distr_pareto_bs;
  vector<DISTR_pareto_p> distr_pareto_ps;
  vector<DISTR_invgaussian_mu> distr_invgaussian_mus;
  vector<DISTR_invgaussian_sigma2> distr_invgaussian_sigma2s;
  vector<DISTR_gengamma_mu> distr_gengamma_mus;
  vector<DISTR_gengamma_sigma> distr_gengamma_sigmas;
  vector<DISTR_gengamma_tau> distr_gengamma_taus;
  vector<DISTR_t_mu> distr_t_mus;
  vector<DISTR_t_sigma2> distr_t_sigma2s;
  vector<DISTR_t_df> distr_t_dfs;
  vector<DISTR_weibull_lambda> distr_weibull_lambdas;
  vector<DISTR_weibull_alpha> distr_weibull_alphas;
  vector<DISTR_dagum_a> distr_dagum_as;
  vector<DISTR_dagum_b> distr_dagum_bs;
  vector<DISTR_dagum_p> distr_dagum_ps;
  vector<DISTR_zeroadjusted> distr_zeroadjusteds;
  vector<DISTR_zeroadjusted_mult> distr_zeroadjusted_mults;
  vector<DISTR_betainf_mu> distr_betainf_mus;
  vector<DISTR_betainf_sigma2> distr_betainf_sigma2s;
  vector<DISTR_betainf_nu> distr_betainf_nus;
  vector<DISTR_betainf_tau> distr_betainf_taus;
  vector<DISTR_betainf0_nu> distr_betainf0_nus;
  vector<DISTR_betainf1_tau> distr_betainf1_taus;
  vector<DISTR_bivt_mu> distr_bivt_mus;
  vector<DISTR_bivt_sigma> distr_bivt_sigmas;
  vector<DISTR_bivt_rho> distr_bivt_rhos;
  vector<DISTR_bivt_df> distr_bivt_dfs;
  vector<DISTR_bivnormal_mu> distr_bivnormal_mus;
  vector<DISTR_bivnormal_sigma> distr_bivnormal_sigmas;
  vector<DISTR_bivnormal_rho> distr_bivnormal_rhos;
  vector<DISTR_bivnormal_mufz> distr_bivnormal_mufzs;
  vector<DISTR_bivnormal_rhofz> distr_bivnormal_rhofzs;
  vector<DISTR_bivprobit_mu> distr_bivprobit_mus;
  vector<DISTR_bivprobit_rho> distr_bivprobit_rhos;
  vector<DISTR_bivlogit_mu> distr_bivlogit_mus;
  vector<DISTR_bivlogit_or> distr_bivlogit_ors;
  vector<DISTR_dirichlet> distr_dirichlets;
  vector<DISTR_BCCG_mu> distr_BCCG_mus;
  vector<DISTR_BCCG_sigma> distr_BCCG_sigmas;
  vector<DISTR_BCCG_nu> distr_BCCG_nus;
  vector<DISTR_copula> distr_copulas;
  vector<DISTR_tcopula_df> distr_tcopula_dfs;
  vector<DISTR_tcopula_rho> distr_tcopula_rhos;
  vector<DISTR_gumbelcopula_rho> distr_gumbelcopula_rhos;
  vector<DISTR_gumbelcopula2_rho> distr_gumbelcopula2_rhos;
  vector<DISTR_gumbelcopula2_normal_mu> distr_gumbelcopula2_normal_mus;
  vector<DISTR_gumbelcopula2_normal_sigma2> distr_gumbelcopula2_normal_sigma2s;
  vector<DISTR_claytoncopula_rho> distr_claytoncopula_rhos;
  vector<DISTR_claytoncopula2_rho> distr_claytoncopula2_rhos;
  vector<DISTR_claytoncopula2_normal_mu> distr_claytoncopula2_normal_mus;
  vector<DISTR_claytoncopula2_normal_sigma2> distr_claytoncopula2_normal_sigma2s;
  vector<DISTR_gaussiancopula_rho> distr_gaussiancopula_rhos;
  vector<DISTR_gaussiancopula_rhofz> distr_gaussiancopula_rhofzs;
  vector<DISTR_frankcopula_rho> distr_frankcopula_rhos;
  vector<DISTR_frankcopula2_rho> distr_frankcopula2_rhos;
  vector<DISTR_frankcopula2_normal_mu> distr_frankcopula2_normal_mus;
  vector<DISTR_frankcopula2_normal_sigma2> distr_frankcopula2_normal_sigma2s;
  vector<DISTR_sfa0_mu_y> distr_sfa0_mu_ys;
  vector<DISTR_sfa0_sigma_u> distr_sfa0_sigma_us;
  vector<DISTR_sfa0_sigma_v> distr_sfa0_sigma_vs;
  vector<DISTR_sfa_mu_y> distr_sfa_mu_ys;
  vector<DISTR_sfa_mu_u> distr_sfa_mu_us;
  vector<DISTR_sfa_mu_y_id> distr_sfa_mu_y_ids;
  vector<DISTR_sfa_mu_u_id> distr_sfa_mu_u_ids;
  vector<DISTR_sfa_sigma_u> distr_sfa_sigma_us;
  vector<DISTR_sfa_sigma_v> distr_sfa_sigma_vs;
  vector<DISTR_sfa_alpha> distr_sfa_alphas;
  vector<DISTR_sfa2_mu_y_id> distr_sfa2_mu_y_ids;
  vector<DISTR_sfa2_mu_u_id> distr_sfa2_mu_u_ids;
  vector<DISTR_sfa2_mu_y> distr_sfa2_mu_ys;
  vector<DISTR_sfa2_mu_u> distr_sfa2_mu_us;
  vector<DISTR_sfa2_sigma_u> distr_sfa2_sigma_us;
  vector<DISTR_sfa2_sigma_v> distr_sfa2_sigma_vs;

  bool create_distribution(void);

  bool resultsyesno;
  bool posteriormode;
  bool computemodeforstartingvalues;

  //----------------------------------------------------------------------------

  // --------------------------- offset handling -------------------------------

  void create_offset(unsigned i);

  // ------------------------- end: offset handling ----------------------------

  //----------------------------------------------------------------------------

  use udata;

  modelterm modreg;
  vector<basic_termtype*> termtypes;

  vector<term> terms;

  term_nonp tnonp;

  //---------------------------- for predict -----------------------------------

  vector<FC_predict> FC_predicts;

  vector<FC_predict_mult> FC_predicts_mult;
  vector<DISTR *> predict_mult_distrs;

  vector<FC_predict_predictor> FC_predict_predictors;

  bool create_predict(void);

  //----------------------------------------------------------------------------


  //-------------------------- for predictive_checks ---------------------------

  vector<FC_predictive_check> FC_predictive_checks;

  void create_predictive_check(void);

  FC_cv FCcv;
  void create_cv(void);

  //----------------------------------------------------------------------------


  //---------------------------- for linear terms ------------------------------

  basic_termtype lineareffects;

  vector<FC_linear> FC_linears;

  bool create_linear(void);

  //----------------------------------------------------------------------------

  //----------------------- for penalized linear terms -------------------------

  int ridge;
  int lasso;
  int ridge_linear;
  int lasso_linear;

  vector<FC_linear_pen> FC_linear_pens;
  vector<FC_variance_pen_vector> FC_variance_pen_vectors;

  bool create_ridge_lasso(unsigned i);


  //----------------------- end for penalized linear terms ---------------------

  //----------------------------------------------------------------------------

  //----------------------- for nonparametric terms ----------------------------

  vector<DESIGN_pspline> design_psplines;
  vector<DESIGN_mrf> design_mrfs;
  vector<DESIGN_kriging> design_krigings;
  vector<FC_nonp> FC_nonps;
  vector<FC_nonp_variance> FC_nonp_variances;
  vector<FC_nonp_variance_varselection> FC_nonp_variance_varselections;

  bool create_nonp(void);
  void create_pspline(unsigned i);
  bool create_mrf(unsigned i);
  bool create_kriging(unsigned i);
  bool create_geokriging(unsigned i);

  bool find_map(unsigned i,MAP::map & m);

//------------------------ end for nonparametric terms -------------------------

//----------------------- hierarchical random effects --------------------------

  vector<DESIGN_hrandom>  design_hrandoms;
  vector<FC_hrandom> FC_hrandoms;
  vector<FC_hrandom_variance> FC_hrandom_variances;
  vector<FC_hrandom_variance_vec> FC_hrandom_variance_vecs;
  vector<FC_hrandom_variance_vec_nmig> FC_hrandom_variance_vec_nmigs;
  vector<FC_hrandom_variance_ssvs> FC_hrandom_variance_ssvss;

  bool create_hrandom(unsigned i);

//------------------- end for hierarchical random effects ----------------------

//---------------------- multiplicative random effects -------------------------

  vector<FC_mult> FC_mults;

  bool create_random_pspline(unsigned i);

//-------------------- end: multiplicative random effects ----------------------

  friend void __EXPORT_TYPE hregressrun(superbayesreg & b);

  bool run_yes;

//------------------------ end: for method regress -----------------------------

//--------------------------- for method autocorr ------------------------------

  intoption maxlag;

  optionlist autocorroptions;

  modelStandard ma;

  useDataset ad;

  friend void __EXPORT_TYPE autocorrrun(superbayesreg & b);

//------------------------ end: for method autocorr ----------------------------

// ------------------------- for method getsample ------------------------------

  optionlist getsampleoptions;

  modelStandard mgetsample;

  useDataset usegetsample;

  friend void __EXPORT_TYPE getsamplerun(superbayesreg & b);

//----------------------- end: for method getsample ----------------------------

  void create(void);
  void create_hregress(void);
  void create_autocorr(void);
  void create_getsample(void);

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_pointer * adminp_p;
  #endif

  public:


  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR

  superbayesreg (void)  : statobject()
    {
    type = "mcmcreg";
    resultsyesno = false;
    }

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - name = n

  #if defined(JAVA_OUTPUT_WINDOW)
  superbayesreg (administrator_basic * adb, administrator_pointer * adp,
                 const ST::string & n,ofstream * lo,istream * i,
                                ST::string p,vector<statobject*> * st);
  #else
  superbayesreg (const ST::string & n,ofstream * lo,istream * i,
                 ST::string p,vector<statobject*> * st);
  #endif

  // COPY CONSTRUCTOR

  superbayesreg (const superbayesreg & b);

  // DESTRUCTOR

  ~superbayesreg()
         {
         }

  // OVERLOADED ASSIGNMENT OPERATOR

  const superbayesreg & operator=(const superbayesreg & b);


  int parse(const ST::string & c);

  void describe(const optionlist & globaloptions = optionlist());

  };

#endif

