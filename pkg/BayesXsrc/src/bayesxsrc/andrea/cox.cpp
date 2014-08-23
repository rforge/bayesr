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





#include "cox.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_coxmodel ---------------------------
//------------------------------------------------------------------------------

DISTRIBUTION_coxmodel::DISTRIBUTION_coxmodel(MCMCoptions * o,
                                          const datamatrix & r,
                                          const datamatrix & t,
                                          const datamatrix & dbeg,
                                          const datamatrix & w)

   : DISTRIBUTION(o,r,w)
     {

     unsigned i;

     nrcat = 1; // ändern bei competing risk

     ti = t;

     int_ti=datamatrix(2*t.rows(),1,0.0);
     for(i=0;i<t.rows();i++)
     {
     int_ti(i,0) = t(i,0)- dbeg(i,0);
     int_ti(t.rows()+i,0)=0.0;
     }

     family = "cox";
     scale(0,0) = 1;
     scaleexisting = false;
     offsetexisting = false;

     }

DISTRIBUTION_coxmodel::DISTRIBUTION_coxmodel(const datamatrix & offset,
                   MCMCoptions * o, const datamatrix & r, const datamatrix & t,
                   const datamatrix & dbeg,
                        const datamatrix & w)
  : DISTRIBUTION(datamatrix(offset.rows(),1,0),o,r,w)                        
{
     unsigned i;

     nrcat = 1; // ändern bei competing risk

     ti = t;
     relrisk = offset;
     
     int_ti=datamatrix(2*t.rows(),1,0.0);
     for(i=0;i<t.rows();i++)
     {
     int_ti(i,0) = t(i,0)-dbeg(i,0);
     int_ti(t.rows()+i,0)=0.0;
     }

     family = "cox";
     scale(0,0) = 1;
     scaleexisting = false;
     offsetexisting = true;

}


double DISTRIBUTION_coxmodel::loglikelihood(double * response,double * linpred,
                       double * weight,const int & i) const
    {
     if(offsetexisting==false)
       return *weight * (*response * (*linpred) - exp(*linpred)* *(int_ti.getV()+i)); //int_ti = 1/(exp(beta0_ti)) *integral(0 bis ti) exp(beta_0u)du

     else
      {
//      double test = relrisk(i,0);
      return *weight * (*response * log(relrisk(i,0) + exp(*linpred)) - exp(*linpred)* *(int_ti.getV()+i));
      }
    }





double DISTRIBUTION_coxmodel::compute_weight(double * linpred, double * weight,const int & i,
                                             const unsigned & col) const
  {
  double weighthelp = exp(*linpred)* *(int_ti.getV()+i);
  if(offsetexisting == false)
    return  *weight * weighthelp;
  else
    {
    double help = *weight * (weighthelp - relrisk(i,0)*  *(response.getV()+i) *exp(*linpred)/((relrisk(i,0)+exp(*linpred)) * (relrisk(i,0)+exp(*linpred))));
    if (help<0.0) help = 0.000001;
    return help;
    }
  }



void  DISTRIBUTION_coxmodel::tilde_y(datamatrix & tildey,datamatrix & m,const unsigned & col,
                                     const bool & current,const datamatrix & w)
 {
  unsigned i;

  double * workspline = m.getV();
  double * workresponse = response.getV();
  double * workweight = w.getV();
  double * ywork = tildey.getV();

  if(offsetexisting == false)
    {
    for (i=0;i<nrobs;i++,workspline++,ywork++,workresponse++,workweight++)
      {
      if(*workweight == 0.0)
        *ywork = 0.0;
      else
        *ywork = *workspline + *workresponse/(*workweight)-1.0;
      }
    }
  else
    {
    double * relrisk_work = relrisk.getV();
    double * linpred = linearpred.getV();
    double * int_ti_work = int_ti.getV();
    for (i=0;i<nrobs;i++,workspline++,ywork++,workresponse++,workweight++,relrisk_work++,linpred++,int_ti_work++)
      {
      if(*workweight == 0.0)
        *ywork = 0.0;
      else
        {
        double deltastar = *workresponse* exp(*linpred)/(*relrisk_work+ exp(*linpred));
        double weighthelp = exp(*linpred)* *int_ti_work;
        *ywork = *workspline + (deltastar - weighthelp)/(*workweight);
        }
      }
    }
 }


void DISTRIBUTION_coxmodel::compute_iwls(void)
  {
  unsigned i;

  double * linpred = linearpred.getV();
  double * workresponse = response.getV();
  double * workweightiwls = weightiwls.getV();
  double * ywork = tildey.getV();
  double * int_ti_work = int_ti.getV();
  double * workweight = weight.getV();

  if(offsetexisting == false)
    {
    for (i=0;i<nrobs;i++,linpred++,ywork++,workresponse++,workweightiwls++,int_ti_work++,workweight++)
      {
      *workweightiwls = *workweight * exp(*linpred)* *int_ti_work;
      if (*workweightiwls==0.0)
        *ywork = 0.0;
      else
        *ywork = *linpred + *workresponse/(*workweightiwls)-1.0;
      }
    }
  else
    {
    double * relrisk_work = relrisk.getV();
    for (i=0;i<nrobs;i++,linpred++,ywork++,workresponse++,workweightiwls++,int_ti_work++,relrisk_work++,workweight++)
      {
      double weighthelp = exp(*linpred)* *int_ti_work;
      double deltastar = *workresponse*exp(*linpred)/(*relrisk_work + exp(*linpred));
      *workweightiwls = *workweight * (weighthelp - (*relrisk_work* *workresponse*exp(*linpred))/
                                        ((*relrisk_work+exp(*linpred))*(*relrisk_work+exp(*linpred))));
      if (*workweightiwls<0.0) *workweightiwls = 0.000001;
      if (*workweightiwls==0.0)
          *ywork = 0.0;
      else
          *ywork = *linpred + (deltastar - weighthelp)/(*workweightiwls);
      }
    }
  }


void DISTRIBUTION_coxmodel::compute_IWLS_weight_tildey(double * response,double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col)
  {
  double weighthelp = exp(*linpred)* *(int_ti.getV() + i);
  if(offsetexisting==false)
    {
    *weightiwls = *weight * weighthelp;
    if (*weightiwls==0.0)
        *tildey = 0.0;
    else
        *tildey = *response/(*weightiwls)-1.0;
    }
  else
    {
    *weightiwls = *weight * (weighthelp - relrisk(i,0)* *response *exp(*linpred)/
                            ((relrisk(i,0)+exp(*linpred))*(relrisk(i,0)+exp(*linpred))));
    if(*weightiwls<0.0) *weightiwls=0.000001;
    double deltastern = *response * exp(*linpred)/(relrisk(i,0)+exp(*linpred));
    if (*weightiwls==0.0)
        *tildey = 0.0;
    else
        *tildey = (deltastern - weighthelp)/(*weightiwls);
    }
  }


double DISTRIBUTION_coxmodel::compute_IWLS(double * response,double * linpred,double * weight,const int & i,
                              double * weightiwls, double * tildey,
                              bool weightyes,const unsigned & col)
  {
  double weighthelp = exp(*linpred)* *(int_ti.getV() + i);

  if(offsetexisting == false)
    {
    if(weightyes)
      *weightiwls = *weight * weighthelp;
    if (*weightiwls==0.0)
       *tildey = 0.0;
    else
       *tildey = *response/(*weightiwls)-1.0;
    return *weight * (*response * (*linpred) - *weightiwls);
    }
  else
    {
    if(weightyes)
      {
      *weightiwls = *weight * (weighthelp- relrisk(i,0)* *response*exp(*linpred)/
                                ((relrisk(i,0)+exp(*linpred))*(relrisk(i,0)+exp(*linpred))));
      if (*weightiwls<0.0) *weightiwls = 0.000001;
      }
    double deltastern = *response * exp(*linpred)/(relrisk(i,0)+exp(*linpred));
    if (*weightiwls==0.0)
        *tildey = 0.0;
    else
        *tildey = (deltastern - weighthelp)/(*weightiwls);
    return *weight * (*response * log(relrisk(i,0) + exp(*linpred)) - weighthelp);
    }
  }


void DISTRIBUTION_coxmodel::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTRIBUTION_coxmodel::update(void)
  {
  DISTRIBUTION::update();
  }


bool DISTRIBUTION_coxmodel::posteriormode(void)
  {
  return true;
  }

/*
void DISTRIBUTION_coxmodel::update_predict(void)
  {

  }
*/

void DISTRIBUTION_coxmodel::outresults(void)
  {
  DISTRIBUTION::outresults();
  }

void DISTRIBUTION_coxmodel::compute_deviance(const double * response, const double * weight,
                           const double * mu, double * deviance, double * deviancesat,
                           const datamatrix & scale, const int & i) const
    {
    double help = *mu * *(int_ti.getV()+i);

    if(offsetexisting == false)
      {
      *deviance = -2.0 * *weight *( *response * log(*mu) - help );
      *deviancesat = -2.0* *weight * (*response - help + *response * log(help));
      }
    else
      {
      *deviance = -2.0 * *weight *( *response * log(relrisk(i,0)+ *mu) - help); 
      }
    }


} // END: namespace MCMC

//---------------------------------------------------------------------------
#if !defined(__BUILDING_GNU)
#pragma package(smart_init)
#endif









