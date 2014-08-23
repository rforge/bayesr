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



#if !defined (MCMCconst_INCLUDED)

#define MCMCconst_INCLUDED

#include"../export_type.h"
#include"mcmc.h"
#include"fullcond.h"
#include"distribution.h"
#include"nbinomial.h"
#include"zip.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FULLCOND_const ----------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_const : public FULLCOND
  {

  protected:

  double lambda;

  bool negbin;

  bool interceptyes;
  int interceptpos;
  double interceptadd;
  datamatrix effectsadd;

  // Shrinkage
  bool shrinkage;
  datamatrix variances;
  bool use_effectstart;
  datamatrix effectstart;
  // Shrinkage

  unsigned nrconst;                // number of fixed effects paramters

  datamatrix linold;
  datamatrix linnew;
  datamatrix * linnewp;
  datamatrix * linoldp;

  DISTRIBUTION * likep;

  double sumold;

  vector<ST::string> table_results;

  void transfer_interceptsample(void);

  // REML

  vector<bool> catspecific_fixed;
  unsigned nrvars;   // no. of variables involved (equal to no.of fixed effects
                     // if no category specific covariates are present)
  datamatrix cats;   // names of the categories
  int catspecific_effects; // no. of category specific covariates
  bool ismultinomialcatsp; // indicator for multinomial models with cat specific covariates

  // End: REML

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_const(void) : FULLCOND()
    {
    }


  // CONSTRUCTOR_1 (diffuse prior)


  FULLCOND_const(MCMCoptions * o,DISTRIBUTION * dp,const datamatrix & d,
                 const ST::string & t, const int & constant,
                 const ST::string & fs,const ST::string & fr,
                 const unsigned & c=0);

  // CONSTRUCTOR_2 (multivariate gaussian prior for beta with mean 'prmean' and
  //                covariance matrix 'prcov')

  FULLCOND_const(MCMCoptions * o,DISTRIBUTION * dp,const datamatrix & d,
                 const datamatrix & prmean, const datamatrix & prcov,
                 const ST::string & t, const int & constant, const ST::string & fs);

  // CONSTRUCTOR_3 (REML)

  FULLCOND_const(MCMCoptions * op, const datamatrix & d, const ST::string & t,
                 const int & constant, const ST::string & fs,
                 const ST::string & fr, const vector<bool> & catsp,
                 const unsigned & np, const unsigned & nrpar,
                 const datamatrix & c, const bool & ismcatsp);


  // COPY CONSTRUCTOR

  FULLCOND_const(const FULLCOND_const & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_const & operator=(const FULLCOND_const & m);

  // shrinkage
  void update_variances(datamatrix & v);
  datamatrix get_variances(void);

   double * getvariancespointer(void)
     {
     return variances.getV();
     }
  // shrinkage

  void update(void);

  virtual void update_intercept(double & m)
    {
    }

  virtual void update_interceptold(double & m)
    {
    }

  bool posteriormode(void);

  virtual void posteriormode_intercept(double & m)
    {
    }

  void init_name(const ST::string & na);

  void init_names(const vector<ST::string> & na);

  void get_effectmatrix(datamatrix & e,vector<ST::string> & enames,
                        unsigned be, unsigned en, effecttype t);

  unsigned get_nreffects(effecttype t);

  void outresults(void);

  void outoptions(void);

  bool is_missing(const ST::string & na);

  bool set_negbin(void)
    {
    negbin=true;
    return negbin;
    }

  virtual void update_fix_effect(const unsigned & pos, double & value, datamatrix fix)
    {
    }

  virtual void update_fix_varcoeff(double & value,ST::string & name);

  virtual void posteriormode_fix_varcoeff(double & value,ST::string & name);


  virtual void posteriormode_const_varcoeff(datamatrix newx)
    {
    }

  // ------------------------------- for REML ----------------------------------

  void createreml(datamatrix & X,datamatrix & Z,const unsigned & Xpos,
                  const unsigned & Zpos);

  double outresultsreml(datamatrix & X,datamatrix & Z,datamatrix & betareml,
                      datamatrix & betacov, datamatrix & thetareml,
                      const unsigned & Xpos, const unsigned & Zpos,
                      const unsigned & thetapos, const bool & dispers,
                      const unsigned & betaXpos,
                      const unsigned & betaZpos,
                      const double & category,
                      const bool & ismultinomial,
                      const unsigned plotpos);

  void outresultsreml_ordinal(datamatrix & X,datamatrix & Z,datamatrix & betareml,
                      datamatrix & betacov, unsigned nrcat2);

  void outoptionsreml()
    {
    }

  vector<bool> get_catspecific_fixed();

  // end: for reml

  double compute_df(void);

  };


//------------------------------------------------------------------------------
//---------------------- CLASS: FULLCOND_const_gaussian ------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_const_gaussian : public FULLCOND_const
  {

  protected:

  bool changingweight;

  datamatrix X1;                   // (X'WX)^-0.5
  datamatrix X2;                   // (X'WX)^-1X'W

  datamatrix help;

  datamatrix mu1;

  // FUNCTION: compute_matrices
  // TASK: computes X1 = (X'WX)^-0.5
  //       computes X2 = (X'WX)^-1X'W

  void compute_matrices(void);

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_const_gaussian(void) : FULLCOND_const()
    {
    }

  //CONSTRUCTOR1

  FULLCOND_const_gaussian(MCMCoptions * o,DISTRIBUTION * dp,const datamatrix & d,
                          const ST::string & t,const int & constant,
                          const ST::string & fs,const ST::string & fr,
                          const bool & r, const datamatrix vars,
                          const bool & useeff, const datamatrix eff,            //NEW
                          const unsigned & c=0);

  //CONSTRUCTOR2 (for factor)

  FULLCOND_const_gaussian(MCMCoptions * o,
                 DISTRIBUTION * dp,const datamatrix & d,
                 const ST::string & code, int & ref,
                 const ST::string & t,const ST::string & fs,
                 const ST::string & fr,const unsigned & c=0);


  // COPY CONSTRUCTOR

  FULLCOND_const_gaussian(const FULLCOND_const_gaussian & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_const_gaussian & operator=(const FULLCOND_const_gaussian & m);

  void update(void);

  void update_intercept(double & m);

  void posteriormode_intercept(double & m);

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void outresults(void);

  void outoptions(void)
    {
    FULLCOND_const::outoptions();
    }


  void update_missings(datamatrix & x,const datamatrix & linpredx,
                       const statmatrix<unsigned> & w,
                       const ST::string & na,double & scalex);

  };



//------------------------------------------------------------------------------
//---------------------- CLASS: FULLCOND_const_gaussian_re ---------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_const_gaussian_re : public FULLCOND_const_gaussian
  {

  protected:

  FULLCOND_const * fc_intercept;


  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_const_gaussian_re(void) : FULLCOND_const_gaussian()
    {
    }

  //CONSTRUCTOR1

  FULLCOND_const_gaussian_re(MCMCoptions * o,DISTRIBUTION * dp,
                             const datamatrix & d,
                          const ST::string & t,const int & constant,
                          const ST::string & fs,const ST::string & fr,
                          const bool & r, const datamatrix vars,
                          const bool & useeff, const datamatrix eff,            //NEW
                          const unsigned & c=0);


  // COPY CONSTRUCTOR

  FULLCOND_const_gaussian_re(const FULLCOND_const_gaussian_re & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_const_gaussian_re & operator=(const FULLCOND_const_gaussian_re & m);

  void update_intercept(double & m);

  void posteriormode_intercept(double & m);


  bool posteriormode_converged(const unsigned & itnr);

  void set_fcintercept(FULLCOND_const * fci)
    {
    fc_intercept= fci;
    }

  void update(void);

  bool posteriormode(void);

  void outresults(void);

  void outoptions(void);

  };



//------------------------------------------------------------------------------
//------------------- CLASS: FULLCOND_const_nongaussian ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FULLCOND_const_nongaussian : public FULLCOND_const
  {

  protected:

  datamatrix proposal;
  datamatrix weightiwls;
  datamatrix XWX;
  datamatrix diff;
  datamatrix tildey;

  datamatrix XWXold;
  datamatrix muy;
  datamatrix mode;
  datamatrix linmode;
  datamatrix help;

  unsigned step;

  ofstream out;

  void compute_XWX(datamatrix & XWXw);

  void compute_XWtildey(datamatrix * b);

  void update_iwls(void);

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_const_nongaussian(void) : FULLCOND_const() {}


  // CONSTRUCTOR

  FULLCOND_const_nongaussian(MCMCoptions* o,DISTRIBUTION * dp,const datamatrix & d,
                             const ST::string & t, const int & constant,
                             const ST::string & fs, const ST::string & fr,
                             const bool & r, const datamatrix vars,
                             const bool & useeff, const datamatrix eff,            //NEW
                             const unsigned & c=0);

  //CONSTRUCTOR2 (for factor)

  FULLCOND_const_nongaussian(MCMCoptions * o,
                 DISTRIBUTION * dp,const datamatrix & d,
                 const ST::string & code, int & ref,
                 const ST::string & t,const ST::string & fs,
                 const ST::string & fr,const unsigned & c=0);


  // COPY CONSTRUCTOR

  FULLCOND_const_nongaussian(const FULLCOND_const_nongaussian & m);

  //OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_const_nongaussian & operator=(
                           const FULLCOND_const_nongaussian & m);

  void update(void);

  void update_intercept(double & m);

  bool posteriormode(void);

  void posteriormode_intercept(double & m);

  bool posteriormode_converged(const unsigned & itnr);

  void outresults(void)
    {
    FULLCOND_const::outresults();
    }

  void outoptions(void)
    {
    FULLCOND_const::outoptions();
    }


  ~FULLCOND_const_nongaussian() {}

  };


class FULLCOND_const_gamma : public FULLCOND_const
  {
  };


//------------------------------------------------------------------------------
//--------------------- CLASS: FULLCOND_const_nbinomial ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FULLCOND_const_nbinomial : public FULLCOND_const_nongaussian
  {
  protected:

    DISTRIBUTION_nbinomial * nblikep;


  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_const_nbinomial(void) : FULLCOND_const_nongaussian(){}


  // CONSTRUCTOR

  FULLCOND_const_nbinomial(MCMCoptions* o,DISTRIBUTION * dp,DISTRIBUTION_nbinomial * nb,
                             const datamatrix & d, const ST::string & t,
                             const int & constant, const ST::string & fs,
                             const ST::string & fr, const bool & r, const datamatrix & vars,
                             const bool & useeff, const datamatrix eff,            //NEW
                             const unsigned & c=0);

  // COPY CONSTRUCTOR

  FULLCOND_const_nbinomial(const FULLCOND_const_nbinomial & m);

  //OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_const_nbinomial & operator=(
                           const FULLCOND_const_nbinomial & m);

  bool posteriormode_converged(const unsigned & itnr);

  void update_intercept(double & m);

  void update(void);

  double update_hierint(void) const;

  };


} // end: namespace MCMC

#endif
