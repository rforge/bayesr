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



#if !defined (VARIANCENONP_INCLUDED)
#define VARIANCENONP_INCLUDED

#include"../export_type.h"
#include"mcmc_nonpbasis.h"
#include"randomeffect.h"
#include"hrandom.h"
#include"mcmc_nonp.h"

namespace MCMC
{


class __EXPORT_TYPE FULLCOND_variance_nonp : public FULLCOND
  {

  protected:

  bool constlambda;
  bool uniformprior;
  bool Laplace;
  bool stationary;
  bool alphafix;

  FULLCOND_nonp_basis * Kp;

  FULLCOND_random * REp;
  bool randomeffect;


  FULLCOND_hrandom * hREp;
  bool hrandom;


  FULLCOND_nonp * Fnp;
  bool fullcondnonp;

  bool update_sigma2;

  DISTRIBUTION * distrp;

  double a_invgamma;
  double b_invgamma;
  unsigned rankK;

  ST::string pathresults;

  bool average;

  unsigned column;

  double scale;
  bool discrete;
  unsigned df;
  vector<double> tau;
  vector<double> lambda;

  FULLCOND fc_lambda;
  void outresults_lambda(void);

  public:


  // DEFAULT CONSTRUCTOR

  FULLCOND_variance_nonp(void) : FULLCOND()
    {
    }


  // CONSTRUCTOR1

  FULLCOND_variance_nonp(MCMCoptions * o,FULLCOND_nonp_basis * p,
                         DISTRIBUTION * d,
                         const double & a,const double & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,const bool & av,
                         const unsigned & c);

  // CONSTRUCTOR2

  FULLCOND_variance_nonp(MCMCoptions * o,FULLCOND_random * p,
                         DISTRIBUTION * d,
                         const double & a,const double & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,const bool & av,
                         const unsigned & c);

  // CONSTRUCTOR3

  FULLCOND_variance_nonp(MCMCoptions * o,FULLCOND_nonp * p,
                         DISTRIBUTION * d,
                         const double & a,const double & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,const bool & av,
                         const unsigned & c);

  // CONSTRUCTOR4

  FULLCOND_variance_nonp(MCMCoptions * o,FULLCOND_hrandom * p,
                         DISTRIBUTION * d,
                         const double & a,const double & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,const bool & av,
                         const unsigned & c);


  // COPY CONSTRUCTOR

  FULLCOND_variance_nonp(const FULLCOND_variance_nonp & t);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_variance_nonp & operator=(const FULLCOND_variance_nonp & t);

  void update(void);

  void update_stationary(void);

  void set_stationary(double alphastart,bool afix=false);

  bool posteriormode(void);

  void outresults(void);

  void outoptions(void);

  void set_update_sigma2(void)
    {
    update_sigma2 = false;
    average = false;
    }

  void set_Laplace(void)
    {
    Laplace = true;
    }

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND::reset();
    setbeta(1,1,0.1);
    }

  // FUNCTION: predict (virtual)
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred)
    {
    }

  // FUNCTION: set_constlambda
  // TASK: variance/smoothingparameter is held constant

  void set_constlambda(void)
    {
    constlambda=true;
    }

  void set_uniformprior(void)
    {
    uniformprior=true;
    }

  void set_discrete(const unsigned & dof)
    {
    discrete=true;
    df=dof;
    setbeta(2,1,distrp->get_scale(column,column)/Kp->getlambda());
    }

  // DESTRUCTOR

  ~FULLCOND_variance_nonp() {}

  }; // end: class FULLCOND_variance_nonp



} // end: namespace MCMC

#endif


