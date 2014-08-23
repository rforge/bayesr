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

#include "distr_zeroadjusted.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS DISTR_zeroadjusted -------------------------
//------------------------------------------------------------------------------


DISTR_zeroadjusted::DISTR_zeroadjusted(GENERAL_OPTIONS * o,DISTR* dpi,
                                       DISTR* dmu)


  {
  family="zero adjusted";

  predict_mult = true;
  outpredictor = false;
  outexpectation = true;
  predictor_name = "overall";


  optionsp = o;
  distrp_pi = dpi;
  distrp_mu = dmu;

  response = distrp_pi->response;
  workingresponse = response;

  maindistribution=true;

  nrobs = response.rows();

  weight = datamatrix(response.rows(),1,1);

  workingweight = weight;

  linearpred1 = datamatrix(nrobs,1,0);
  linearpred2 = datamatrix(nrobs,1,0);

  linpred_current = 1;

  helpmat1 = datamatrix(1,1,1);

  }


const DISTR_zeroadjusted & DISTR_zeroadjusted::operator=(
const DISTR_zeroadjusted & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  distrp_pi = nd.distrp_pi;
  distrp_mu = nd.distrp_mu;
  return *this;
  }


DISTR_zeroadjusted::DISTR_zeroadjusted(const DISTR_zeroadjusted & nd)
   : DISTR(DISTR(nd))
  {
  distrp_pi = nd.distrp_pi;
  distrp_mu = nd.distrp_mu;
  }


void DISTR_zeroadjusted::outoptions(void)
  {

  }


void DISTR_zeroadjusted::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {

  double pi;
  double E;
  distrp_pi->compute_mu(linpred[0],&pi);
  distrp_mu->compute_mu(linpred[1],&E);

  *mu = pi*E;

  }


datamatrix * DISTR_zeroadjusted::get_auxiliary_parameter(
                                           auxiliarytype t)
  {
  if (t == auxcurrent)
    helpmat1(0,0) = distrp_mu->get_scale();
  else
    helpmat1(0,0) = distrp_mu->get_scalemean();

  return &helpmat1;
  }



void DISTR_zeroadjusted::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

  double s=1;
  double pi;
  distrp_pi->compute_mu(linpred[0],&pi);

  double d_pi;
  distrp_pi->compute_deviance(response[0],weight[0],&pi,&d_pi,&s);

  if (*response[0] == 0)
    {

    *deviance = d_pi;

    }
  else
    {
    double mu;
    distrp_mu->compute_mu(linpred[1],&mu);

    double s_mu= (*aux[2])(0,0);

    double d_mu;
    distrp_mu->compute_deviance(response[1],weight[1],&mu,&d_mu,&s_mu);

    *deviance = d_pi + d_mu;

    }
  }


//------------------------------------------------------------------------------
//------------------------ CLASS DISTR_zeroadjusted_mult -----------------------
//------------------------------------------------------------------------------


DISTR_zeroadjusted_mult::DISTR_zeroadjusted_mult(GENERAL_OPTIONS * o,DISTR* dpi,
                                       vector<DISTR*> dmu)


  {
  family="zero adjusted";

  predict_mult = true;
  outpredictor = false;
  outexpectation = true;
  predictor_name = "overall";

  optionsp = o;
  distrp_pi = dpi;
  distrp_mu = dmu;

  distrp_mu[1]->predstart_mumult=1;

  response = distrp_pi->response;
  workingresponse = response;

  maindistribution=true;

  nrobs = response.rows();

  weight = datamatrix(response.rows(),1,1);

  workingweight = weight;

  linearpred1 = datamatrix(nrobs,1,0);
  linearpred2 = datamatrix(nrobs,1,0);

  linpred_current = 1;

  helpmat1 = datamatrix(1,1,1);
  linpredvec = vector<double*>(distrp_mu.size());
  responsevec = vector<double*>(distrp_mu.size());
  weightvec = vector<double*>(distrp_mu.size());

  }


const DISTR_zeroadjusted_mult & DISTR_zeroadjusted_mult::operator=(
const DISTR_zeroadjusted_mult & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  distrp_pi = nd.distrp_pi;
  distrp_mu = nd.distrp_mu;
  linpredvec = nd.linpredvec;
  weightvec = nd.weightvec;
  responsevec = nd.responsevec;
  return *this;
  }


DISTR_zeroadjusted_mult::DISTR_zeroadjusted_mult(
                         const DISTR_zeroadjusted_mult & nd)
   : DISTR(DISTR(nd))
  {
  distrp_pi = nd.distrp_pi;
  distrp_mu = nd.distrp_mu;
  linpredvec = nd.linpredvec;
  weightvec = nd.weightvec;
  responsevec = nd.responsevec;
  }


void DISTR_zeroadjusted_mult::outoptions(void)
  {

  }


void DISTR_zeroadjusted_mult::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {

  double pi;
  double E;
  distrp_pi->compute_mu(linpred[0],&pi);

/*
  unsigned i;
  for (i=0;i<linpredvec.size();i++)
    linpredvec[i] = linpred[i+1];
*/

  distrp_mu[distrp_mu.size()-1]->compute_mu_mult(linpred,response,&E);

  *mu = pi*E;

  }



void DISTR_zeroadjusted_mult::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux)
  {

  double s=1;
  double pi;
  distrp_pi->compute_mu(linpred[0],&pi);

  double d_pi;
  distrp_pi->compute_deviance(response[0],weight[0],&pi,&d_pi,&s);

  if (*response[0] == 0)
    {
    *deviance = d_pi;
    }
  else
    {

    double d_mu;
    unsigned i;
    for (i=0;i<linpredvec.size();i++)
      {
      linpredvec[i] = linpred[i+1];
      weightvec[i] = weight[i+1];
      responsevec[i] = response[i+1];
      }

    distrp_mu[distrp_mu.size()-1]->compute_deviance_mult(responsevec,weightvec,
                                                         linpredvec,&d_mu,aux);
    *deviance = d_pi + d_mu;

    }

  }


} // end: namespace MCMC



