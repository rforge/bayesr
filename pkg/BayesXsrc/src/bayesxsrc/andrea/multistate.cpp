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





#include "multistate.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_multistatemodel -------------------
//------------------------------------------------------------------------------

DISTRIBUTION_multistatemodel::DISTRIBUTION_multistatemodel(MCMCoptions * o,
                                          const datamatrix & r,
                                          const datamatrix & t,
                                          const datamatrix & dbeg,
                                          const datamatrix & state,
                                          const datamatrix & w)

   : DISTRIBUTION(o,r,w)
     {

     init_names("Transition indicators");

     unsigned i,j,k;

     nrtransition = response.cols();
     ti = t;
     state_i = state;
     beg = dbeg;

     nrcat = state_i.max(0);

     transition_help = datamatrix(nrcat,nrtransition,0.0);
     for(i=0;i<ti.rows();i++)
       {
       for(j=1;j<=nrcat;j++)
         {
         for(k=0;k<nrtransition;k++)
           {
           if(state_i(i,0)==j && response(i,k)==1.0) transition_help(j-1,k)+=1.0;
           }
         }
       }

     optionsp->out("\n");
     optionsp->out("Matrix of possible transitions:\n");
     optionsp->out("\n");

     ST::string str = "\tTransition\t";
     for(j=0;j<transition_help.cols();j++)
       {
       str = str + ST::inttostring(j+1) + "\t";
       }
     optionsp->out(str);
     optionsp->out("State");

     for(i=0;i<transition_help.rows();i++)
       {
       ST::string str = ST::inttostring(i+1) + "\t\t\t";
       for(j=0;j<transition_help.cols();j++)
         {
         str = str + ST::doubletostring(transition_help(i,j)) + "\t";
         }
       optionsp->out(str);
       }
     optionsp->out("\n");

     transition = datamatrix(ti.rows(),nrtransition,0.0);
     for(i=0;i<ti.rows();i++)
       {
       for(j=1;j<=nrcat;j++)
         {
         if(state_i(i,0)==j)
           {
           for(k=0;k<nrtransition;k++)
             {
             if(transition_help(j-1,k)>0.0) transition(i,k)=1.0;
             }
           }
         }
       }

     int_ti=datamatrix(2*t.rows(),nrtransition,0.0);
/*     for(i=0;i<t.rows();i++)
       {
       for(j=0;j<nrtransition;j++)
         {
         if(transition(i,j)==1.0)
           {
           int_ti(i,j) = ti(i,0)-beg(i,0);
           double testti= ti(i,0);
           double testbeg=beg(i,0);
           double testint=int_ti(i,j);
           }
         else
           int_ti(i,j) = 0.0;
         int_ti(t.rows()+i,j)=0.0;
         }
       }      */

/*     ofstream intout("d:\\temp\\int_ti_start.txt");
     for(i=0;i<t.rows();i++)
       {
       for(j=0;j<nrtransition;j++)
         {
         intout<<int_ti(i,j)<<" ";
         }
       intout<<endl;
       }
     intout.close();     */

     family = "multistate";
     scale(0,0) = 1;
     scaleexisting = false;
//     offsetexisting = false;

     }

/*DISTRIBUTION_multistatemodel::DISTRIBUTION_multistatemodel(const datamatrix & offset,
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

     family = "multistate";
     scale(0,0) = 1;
     scaleexisting = false;
     offsetexisting = true;

} */


double DISTRIBUTION_multistatemodel::loglikelihood(double * response,double * linpred,
                       double * weight,const int & i) const
    {
    double erg=0.0;
    for(unsigned j=0; j<nrtransition;j++,response++,linpred++)
      {
      erg = erg + *response * *(linpred);
      if( transition(i,j)==1.0)
        erg = erg - exp(*linpred)* *(int_ti.getV()+nrtransition*i+j);
      }
    return *weight * erg;
    }


void DISTRIBUTION_multistatemodel::compute_mu(const double * linpred,double * mu)
                                           const
  {
  for(unsigned i=0;i<linearpred.cols();i++,linpred++,mu++)
    *mu = exp( *linpred);
  }



double DISTRIBUTION_multistatemodel::compute_weight(double * linpred, double * weight,const int & i,
                                             const unsigned & col) const
  {

  return  *weight *(transition(i,col)* exp(*(linpred+col))* *(int_ti.getV()+i*nrtransition+col));

  }



void  DISTRIBUTION_multistatemodel::tilde_y(datamatrix & tildey,datamatrix & m,const unsigned & col,
                                     const bool & current,const datamatrix & w)
 {
  unsigned i;

  double * workspline = m.getV();
  double * workresponse = response.getV()+col;
  double * workweight = w.getV();
  double * ywork = tildey.getV();

  for (i=0;i<nrobs;i++,workspline++,ywork++,workresponse+nrtransition,workweight++)
    {
    if(*workweight == 0.0)
      *ywork = 0.0;
    else
      *ywork = *workspline + *workresponse/(*workweight)-1.0;
    }
 }


void DISTRIBUTION_multistatemodel::compute_iwls(void)
  {
  unsigned i,j;

  double * linpred = linearpred.getV();
  double * workresponse = response.getV();
  double * workweightiwls = weightiwls.getV();
  double * ywork = tildey.getV();
  double * int_ti_work = int_ti.getV();


  for (i=0;i<nrobs;i++,linpred++,ywork++,workresponse++,workweightiwls++,int_ti_work++)
    {
    for(j=0;j<nrtransition;j++,linpred++,ywork++,workresponse++,workweightiwls++,int_ti_work++)
      {
      *workweightiwls = transition(i,j)*exp(*linpred)* *int_ti_work;
      if(*workweightiwls==0.0)
        *ywork = 0.0;
      else
        *ywork = *linpred + *workresponse/(*workweightiwls)-1.0;
      }
    }
  }


void DISTRIBUTION_multistatemodel::compute_IWLS_weight_tildey(double * response,double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col)
  {

  double weighthelp = transition(i,col)* exp(*(linpred+col))* *(int_ti.getV()+i*nrtransition+col);

  *weightiwls = weighthelp;
  if(*weightiwls==0.0)
    *tildey = 0.0;
  else
    *tildey = *(response+col)/(*weightiwls)-1.0;

  }


double DISTRIBUTION_multistatemodel::compute_IWLS(double * response,double * linpred,double * weight,const int & i,
                              double * weightiwls, double * tildey,
                              bool weightyes,const unsigned & col)
  {
  double weighthelp = transition(i,col)*exp(*(linpred+col))* *(int_ti.getV()+nrtransition*i+col);


  if(weightyes)
    *weightiwls = weighthelp;
  if(weighthelp==0.0)
    *tildey = 0.0;
  else
    *tildey = *(response+col)/(*weightiwls)-1.0;

  double erg=0.0;
  for(unsigned j=0; j<nrtransition;j++,response++,linpred++)
    {
    erg = erg + *response * *(linpred);
    if( transition(i,j)==1.0)
      erg = erg - exp(*linpred)* *(int_ti.getV()+nrtransition*i+j);
    }
  return *weight * erg;
  }


void DISTRIBUTION_multistatemodel::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTRIBUTION_multistatemodel::update(void)
  {
  DISTRIBUTION::update();
  }


bool DISTRIBUTION_multistatemodel::posteriormode(void)
  {
  return true;
  }

/*
void DISTRIBUTION_multistatemodel::update_predict(void)
  {

  }
*/

void DISTRIBUTION_multistatemodel::outresults(void)
  {
  /*DISTRIBUTION::outresults();*/
  }

void DISTRIBUTION_multistatemodel::compute_deviance(const double * response, const double * weight,
                           const double * mu, double * deviance, double * deviancesat,
                           const datamatrix & scale, const int & i) const
    {
    double erg = 0.0;
    for(unsigned j=0; j<nrtransition;j++,response++,mu++)
      {
      erg = erg + *mu * *(int_ti.getV()+i*nrtransition+j);
      if(transition(i,j)==1.0){
        erg = erg - *mu * *(int_ti.getV()+i*nrtransition+j);
        }
      }
    *deviance = -2.0 * *weight *erg;
    }


} // END: namespace MCMC

//---------------------------------------------------------------------------
#if !defined(__BUILDING_GNU)
#pragma package(smart_init)
#endif

