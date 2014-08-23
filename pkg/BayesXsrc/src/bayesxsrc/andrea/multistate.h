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



#ifndef multistateH
#define multistateH

#include"../export_type.h"
#include "distribution.h"
#include "mcmc_pspline.h"
#include "multibaseline.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//------------------------ CLASS: DISTRIBUTION_multistatemodel ------------------------
//------------------------------------------------------------------------------



class __EXPORT_TYPE DISTRIBUTION_multistatemodel : public DISTRIBUTION
  {

  protected:

  unsigned nrcat;
  unsigned nrtransition;

  datamatrix int_ti;
  datamatrix mean_int_ti;

  datamatrix ti;
  datamatrix state_i;
  datamatrix beg;
  datamatrix transition;
  datamatrix transition_help;
//  datamatrix relrisk;
//  bool offsetexisting;


  public:


  double*  get_integral_ti (void )
  {
  return int_ti.getV();
  }

  double  get_transition (unsigned i, unsigned j)
  {
  if(i<transition.rows())
    return transition(i,j);
  else
  {
//    double testrow= transition.rows();
    return transition(i-transition.rows(),j);
  }
  }


   // DEFAULT CONSTRUCTOR


   DISTRIBUTION_multistatemodel(void) : DISTRIBUTION()
     {
     }

   // CONSTRUCTOR

   DISTRIBUTION_multistatemodel(MCMCoptions * o, const datamatrix & r,
   const datamatrix & t,  const datamatrix & dbeg,
   const datamatrix & state , const datamatrix & w=datamatrix());


   // COPY CONSTRUCTOR

   DISTRIBUTION_multistatemodel(const DISTRIBUTION_multistatemodel & nd)
   : DISTRIBUTION(DISTRIBUTION(nd))
     {

     nrcat = nd.nrcat;
     nrtransition = nd.nrtransition;
     ti=nd.ti;
     state_i = nd.state_i;
     beg = nd.beg;
     transition_help = nd.transition_help;
     transition = nd.transition;

     int_ti = nd.int_ti;
     mean_int_ti = nd.mean_int_ti;
//     offsetexisting = nd.offsetexisting;
//     relrisk = nd.relrisk;

     }

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTRIBUTION_multistatemodel & operator=(const DISTRIBUTION_multistatemodel & nd)
     {

     if (this==&nd)
	   return *this;

     DISTRIBUTION::operator=(DISTRIBUTION(nd));

     nrcat = nd.nrcat;
     nrtransition = nd.nrtransition;
     ti=nd.ti;
     state_i = nd.state_i;
     beg = nd.beg;
     transition_help = nd.transition_help;
     transition = nd.transition;

     int_ti= nd.int_ti;
     mean_int_ti = nd.mean_int_ti;
//     offsetexisting = nd.offsetexisting;
//     relrisk = nd.relrisk;
     return *this;
     }

   // DESTRUCTOR

   ~DISTRIBUTION_multistatemodel() {}

   // FUNCTION: loglikelihood
   // TASK: computes the loglikelihood for a single observation



  double loglikelihood(double * response,double * linpred,
                       double * weight,const int & i) const;




 // FUNCTION: compute_mu
  // TASK: computes mu for a new linear predictor 'linpred'

  void compute_mu(const double * linpred,double * mu) const;



  void compute_mu_notransform(const double * linpred,double * mu) const
    {
    for(unsigned i=0;i<linearpred.cols();i++,linpred++,mu++)
      *mu = exp(*linpred);
    }


  // FUNCTION: compute_devresidual
  // TASK: computes the deviance residual

  void compute_deviance(const double * response,const double * weight,
                           const double * mu,double * deviance,double * deviancesat,
                           const datamatrix & scale,const int & i) const;


  void outoptions(void);

  // FUNCTION: update
  // TASK: updates the scale parameter

  void update(void);

//  void update_predict(void);

  bool posteriormode(void);

  void outresults(void);

  double compute_weight(double * linpred, double * weight, const int & i,
                                                   const unsigned & col) const;

  void tilde_y(datamatrix & tildey,datamatrix & m, const unsigned & col,
                const bool & current, const datamatrix & w);

  void compute_iwls(void);

  void compute_IWLS_weight_tildey(double * response,double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col=0);

  double compute_IWLS(double * response,double * linpred,double * weight,const int & i,
                              double * weightiwls, double * tildey,
                              bool weightyes,const unsigned & col=0);

  void assign_int_ti(const datamatrix & m)
    {
    int_ti = m;
    }

  void assign_mean_int_ti(void)
    {
    mean_int_ti = int_ti;
    }

  datamatrix & get_int_ti(void)
    {
    return int_ti;
    }



  unsigned get_nrcat(void)
    {
    return nrcat;
    }

  unsigned get_nrtransition(void)
    {
    return nrtransition;
    }


  };  // end: class DISTRIBUTION_multistatemodel

} // END: namespace MCMC


//---------------------------------------------------------------------------
#endif
