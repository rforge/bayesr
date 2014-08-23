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



#if !defined (MCMCrandom_INCLUDED)

#define MCMCrandom_INCLUDED

#include"../export_type.h"
#include <iomanip>
using std::setprecision;

#include"mcmc.h"
#include"fullcond.h"
#include"distribution.h"
#include"mcmc_nonpbasis.h"
#include"mcmc_nonp.h"


namespace MCMC
{


class __EXPORT_TYPE FULLCOND_random : public FULLCOND
  {

  protected:

  datamatrix muy;

  FULLCOND_const * fcconst;

  DISTRIBUTION * likep;

  statmatrix<int> index;
  statmatrix<int> index2;

  vector<unsigned>     posbeg;
  vector<unsigned>     posend;

  datamatrix XX;
  datamatrix effvalues;
  double sigma2;                            // prior variance parameter
  double lambda;
  double lambdaold1;
  double lambdaold2;
  double df_lambdaold1;
  double df_lambdaold2;
  bool lambdaconst;

  bool randomslope;
  bool includefixed;

  bool notransform;                         // Results will not be transformed
                                            // for multiplicative random effects

  datamatrix data2;

  bool spatialtotal;
  statmatrix<int> indextotal;
  ST::string pathsample_total;


  FULLCOND ftotal;

  // BEGIN: DSB //
  #if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
  bool mscheck;
  void update_linpred_mscheck(datamatrix & priorSamples, datamatrix & posteriorSamples);
  #endif
  // END: DSB //

  void init_spatialtotal(vector<ST::string> & ev, const ST::string & pnt,
                         const ST::string & prt);

  bool changingweight;

  double centerbeta(void);

  void update_linpred(const bool & add);

  void update_linpred_diff(datamatrix & b1,datamatrix & b2);

  public:

// Begin: DSB
  #if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
  void set_mscheck(const bool & val);
  #endif
// End: DSB

  // DEFAULT CONSTRUCTOR:

  FULLCOND_random(void) : FULLCOND()
    {
    }

  // CONSTRUCTOR1
  // random intercept

  FULLCOND_random(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & d, const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const double & la, const unsigned & c=0);

  // CONSTRUCTOR2
  // random slope

  FULLCOND_random(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & intvar,const datamatrix & effmod,
                  const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const ST::string & prf,
                  const double & la, const bool & inclfix,
                  const unsigned & c=0);

  // CONSTRUCTOR3
  // random intercept (FOR REML)

  FULLCOND_random(MCMCoptions * op,const datamatrix & d, const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const double & la, const double & las, const bool & catsp);

  // CONSTRUCTOR4
  // random slope (FOR REML)

  FULLCOND_random(MCMCoptions * o,
                  const datamatrix & intvar,const datamatrix & effmod,
                  const ST::string & t,const ST::string & fp,
                  const ST::string & pr,const double & la, const double & las,
                  const bool & catsp);

  void init_name(const ST::string & na);

  void init_names(const vector<ST::string> & na);

  // COPY CONSTRUCTOR

  FULLCOND_random(const FULLCOND_random & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_random & operator=(
                        const FULLCOND_random & fc);

  // DESTRUCTOR

  ~FULLCOND_random() {}

  void compute_XWX(const datamatrix & weightmat,const unsigned & col);

  // FUNCTION update
  // TASK: updates parameters (i.e. matrix beta)

  void update(void);

  bool posteriormode(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file (after estimation)

  void outresults(void);

  // FUNCTION: outoptions

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters for a new simulation

  void reset(void)
    {
    FULLCOND::reset();
    sigma2 = 10;
    }

  double compute_quadform(void);

  void update_sigma2(const double & s)
    {
    sigma2 = s;
    }

  unsigned get_rankK(void);

  unsigned get_rankK2(void);

  double getlambda(void)
    {
    return lambda;
    }

  double get_sigma2(void)
    {
    return sigma2;
    }

  void set_lambdaconst(double la);

  void set_notransform(void)
    {
    notransform=true;
    }

  void set_changingweight(void)
    {
    changingweight=true;
    }

  void init_data_varcoeff(const datamatrix & intvar,double add=0);

  void get_effectmatrix(datamatrix & e,vector<ST::string> & enames,
                                unsigned be,unsigned en, effecttype t);



  // FOR REML

  void createreml(datamatrix & X,datamatrix & Z,
                  const unsigned & Xpos,const unsigned & Zpos);

  double outresultsreml(datamatrix & X,datamatrix & Z,datamatrix & betareml,
                      datamatrix & betacov, datamatrix & thetareml,
                      const unsigned & Xpos,const unsigned & Zpos,
                      const unsigned & thetapos, const bool & dispers,
                      const unsigned & betaXpos,
                      const unsigned & betaZpos,
                      const double & category,
                      const bool & ismultinomial,
                      const unsigned plotpos);

  void outoptionsreml();


  };     // end: class FULLCOND_random


//------------------------------------------------------------------------------
//--------------------- class: FULLCOND_random_gaussian ------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_random_gaussian : public FULLCOND_random
  {

  protected:

  datamatrix mu;

  FULLCOND_nonp_basis * fbasisp;

  public:

  // DEFAULT CONSTRUCTOR:

  FULLCOND_random_gaussian(void) : FULLCOND_random()
    {
    }

  // CONSTRUCTOR1
  // random intercept

  FULLCOND_random_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                           FULLCOND_const * fcc,
                           const datamatrix & d, const ST::string & t,
                           const ST::string & fp,const ST::string & pr,
                           const double & la,
                           const unsigned & c = 0)
                           : FULLCOND_random(o,dp,fcc,d,t,fp,pr,la,c)
    {
    mu = datamatrix(index.rows(),1);
    }

  // CONSTRUCTOR2
  // random slope

  FULLCOND_random_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                           FULLCOND_const * fcc,
                           const datamatrix & intvar,
                           const datamatrix & effmod,
                           const ST::string & t,
                           const ST::string & fp,const ST::string & pr,
                           const ST::string & prf,
                           const double & la,
                           const bool & inclfixed,
                           const unsigned & c = 0)
                           : FULLCOND_random(o,dp,fcc,intvar,effmod,t,
                                             fp,pr,prf,la,inclfixed,c)
    {
    mu = datamatrix(index.rows(),1);
    }


  // COPY CONSTRUCTOR

  FULLCOND_random_gaussian(const FULLCOND_random_gaussian & fc)
  : FULLCOND_random(FULLCOND_random(fc))
    {
    mu = fc.mu;
    muy = fc.muy;
    fbasisp = fc.fbasisp;
    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_random_gaussian & operator=(
                        const FULLCOND_random_gaussian & fc)
    {
    if (this==&fc)
      return *this;
    FULLCOND_random::operator=(FULLCOND_random(fc));
    mu = fc.mu;
    muy = fc.muy;
    fbasisp = fc.fbasisp;
    return *this;
    }

  // DESTRUCTOR

  ~FULLCOND_random_gaussian() {}

  // FUNCTION: update
  // TASK: updates parameters (i.e. matrix beta)

  void update(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file (after estimation)

  void outresults(void)
    {
    FULLCOND_random::outresults();
    }

  // FUNCTION: outoptions

  void outoptions(void)
    {
    FULLCOND_random::outoptions();
    }

  void init_spatialtotal(FULLCOND_nonp_basis * sp,const ST::string & pnt,
                         const ST::string & prt);

  // FUNCTION: reset
  // TASK: resets all parameters for a new simulation

  void reset(void)
    {
    FULLCOND_random::reset();
    }

  };     // end: class FULLCOND_random_gaussian


//------------------------------------------------------------------------------
//--------------------- class: FULLCOND_random_nongaussian ---------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FULLCOND_random_nongaussian : public FULLCOND_random
  {

  protected:

  bool iwlsmode;

  datamatrix mode;
  datamatrix modeold;
  datamatrix proposal;
  datamatrix w;
  FULLCOND_nonp * fnonp;
  FULLCOND_nonp_basis * fbasisp;
  bool nongaussian;
  datamatrix weightiwls;
  datamatrix tildey;


  unsigned oldacceptance;
  unsigned oldnrtrials;
  double lambdaprop;
  double a_invgamma;
  double b_invgamma;
  double f;


  void update_spatialtotal(void);


  void update_random_intercept(void);

  void update_random_intercept_singleblock(void);

  void update_random_intercept_iwls_singleblock(void);


  void update_random_slope(void);

  void update_random_slope_singleblock(void);

  void update_random_slope_iwls_singleblock(void);


  void update_random_slope_includefixed(void);

  void update_random_slope_includefixed_singleblock(void);

  void update_random_slope_includefixed_iwls(void);

  void update_random_slope_includefixed_iwls_singleblock(void);


  double scale_proposal();

  void tune_updatetau(const rate & r);


  public:

  // DEFAULT CONSTRUCTOR:

  FULLCOND_random_nongaussian(void) : FULLCOND_random()
    {
    }

  // CONSTRUCTOR1
  // random intercept

  FULLCOND_random_nongaussian(MCMCoptions * o,DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & d, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,
                              const double & la,const bool & im,
                              const unsigned & c=0);

  // CONSTRUCTOR2
  // random slope

  FULLCOND_random_nongaussian(MCMCoptions * o,DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & intvar,
                              const datamatrix & effmod, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,
                              const ST::string & prf,
                              const double & la,
                              const bool & im,
                              const bool & inclfixed,
                              const unsigned & c=0);


  // COPY CONSTRUCTOR

  FULLCOND_random_nongaussian(const FULLCOND_random_nongaussian & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_random_nongaussian & operator=(
                        const FULLCOND_random_nongaussian & fc);

  // DESTRUCTOR

  ~FULLCOND_random_nongaussian() {}

  // FUNCTION update
  // TASK: updates parameters (i.e. matrix beta)

  void update(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file (after estimation)

  void outresults(void);

  // FUNCTION: outoptions

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters for a new simulation

  void reset(void)
    {
    FULLCOND_random::reset();
    }

  void init_spatialtotal(FULLCOND_nonp * sp, const ST::string & pnt,
                                        const ST::string & prt);

  void init_spatialtotal(FULLCOND_nonp_basis * sp,const ST::string & pnt,
                         const ST::string & prt);


  };     // end: class FULLCOND_random_nongaussian


}   // end: namespace MCMC


#endif

