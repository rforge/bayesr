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



#if !defined (DISTRzeroadjusted_INCLUDED)
#define DISTRzeroadjusted_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include"distr.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//-------------------------- CLASS: DISTR_zeroadjusted -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_zeroadjusted : public DISTR
  {

  protected:


  public:

  DISTR* distrp_pi;
  DISTR* distrp_mu;

   // DEFAULT CONSTRUCTOR

  DISTR_zeroadjusted(void) : DISTR()
    {
    }

   // CONSTRUCTOR

  DISTR_zeroadjusted(GENERAL_OPTIONS * o,DISTR* dpi,DISTR* dmu);

   // COPY CONSTRUCTOR

  DISTR_zeroadjusted(const DISTR_zeroadjusted & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_zeroadjusted & operator=(const DISTR_zeroadjusted & nd);

   // DESTRUCTOR

  ~DISTR_zeroadjusted() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  datamatrix * get_auxiliary_parameter(auxiliarytype t = auxcurrent);

  void outoptions(void);

  };


//------------------------------------------------------------------------------
//---------------------- CLASS: DISTR_zeroadjusted_mult ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_zeroadjusted_mult : public DISTR
  {

  protected:

  DISTR* distrp_pi;
  vector<DISTR*> distrp_mu;
  vector<double*> linpredvec;
  vector<double*> responsevec;
  vector<double*> weightvec;


  public:


   // DEFAULT CONSTRUCTOR

  DISTR_zeroadjusted_mult(void) : DISTR()
    {
    }

   // CONSTRUCTOR

  DISTR_zeroadjusted_mult(GENERAL_OPTIONS * o,DISTR* dpi,vector<DISTR*> dmu);

   // COPY CONSTRUCTOR

  DISTR_zeroadjusted_mult(const DISTR_zeroadjusted_mult & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_zeroadjusted_mult & operator=(const DISTR_zeroadjusted_mult & nd);

   // DESTRUCTOR

  ~DISTR_zeroadjusted_mult() {}

  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix*> aux);

  void compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu);

  void outoptions(void);

  };



} // end: namespace MCMC


#endif
