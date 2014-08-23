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



#if !defined (MCMChrandom_INCLUDED)

#define MCMChrandom_INCLUDED


#include"../export_type.h"
#include"mcmc.h"
#include"fullcond.h"
#include"distribution.h"
#include"mcmc_nonpbasis.h"
#include"mcmc_nonp.h"


namespace MCMC
{


class __EXPORT_TYPE FULLCOND_hrandom : public FULLCOND
  {

  protected:

  datamatrix muy;
  datamatrix mu;

  DISTRIBUTION * likep;
  DISTRIBUTION * likep_RE;

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

  datamatrix data2;

  bool changingweight;

  double centerbeta(void);

  void update_linpred(const bool & add);

//  void update_linpred_diff(datamatrix & b1,datamatrix & b2);

  public:

  // DEFAULT CONSTRUCTOR:

  FULLCOND_hrandom(void) : FULLCOND()
    {
    }

  // CONSTRUCTOR1
  // random intercept

  FULLCOND_hrandom(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & d, const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const double & la, const unsigned & c=0);

  void init_name(const ST::string & na);

  void init_names(const vector<ST::string> & na);

  // COPY CONSTRUCTOR

  FULLCOND_hrandom(const FULLCOND_hrandom & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_hrandom & operator=(
                        const FULLCOND_hrandom & fc);

  // DESTRUCTOR

  ~FULLCOND_hrandom() {}

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


  void set_changingweight(void)
    {
    changingweight=true;
    }

  void get_effectmatrix(datamatrix & e,vector<ST::string> & enames,
                                unsigned be,unsigned en, effecttype t);


  void set_REdistr(DISTRIBUTION * lik_RE)
    {
    likep_RE = lik_RE;
    }

  };     // end: class FULLCOND_random


//------------------------------------------------------------------------------
//--------------------- class: FULLCOND_random_gaussian ------------------------
//------------------------------------------------------------------------------

/*
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
*/



}   // end: namespace MCMC


#endif

