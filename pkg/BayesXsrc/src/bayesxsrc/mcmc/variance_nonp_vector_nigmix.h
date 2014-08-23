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



#if !defined (VARIANCENONP_VECTOR_NIG_INCLUDED)
#define VARIANCENONP_VECTOR_NIG_INCLUDED

#include"../export_type.h"
#include"mcmc_const.h"

namespace MCMC
{


class __EXPORT_TYPE FULLCOND_variance_nonp_vector_nigmix : public FULLCOND
  {

  protected:

  vector<double> tau;               // Varianceparameters
  vector<double> lambda;            // Inverse Varianceparameter: lambda=1/tau^2

  bool update_sigma2;

  unsigned column;                  //  Category for for fullcond if multivariate response

  ST::string pathresults;           //  File path for storing sampled parameters

  vector<FULLCOND_const *> Cp;

  DISTRIBUTION * distrp;

  FULLCOND fc_shrinkage;
  FULLCOND fc_indicator;
  FULLCOND fc_t2;

  vector<double> v_0;               //  Hyperparameter for Varianceparameterkomponent: Indicator
  vector<double> v_1;               //  Hyperparameter for Varianceparameterkomponent: Indicator
  vector<double> a_t2;              //  Hyperparameter for Varianceparameterkomponent: t2
  vector<double> b_t2;              //  Hyperparameter for Varianceparameterkomponent: t2
  vector<double> a_omega;           //  Hyperparameter for Prameter w
  vector<double> b_omega;           //  Hyperparameter for Prameter w
  bool omega_fix;                   //  Mixingparameter fix
  bool omega_adaptive;              //  Mixingparameter fix

  datamatrix indicator;             // Matrix for 1. Varianceparameterkomponent: Indicators
  datamatrix t2;                    // Matrix for 1. Varianceparameterkomponent: t2
  datamatrix startdata;             //  Matrix with Starting values and hyperparameters

  double nigmixsum;                 //  sum(beta^2/tau^2) for update of scaleparameter

  vector<unsigned> cut;             //  Blocks of regression coefficients

  void outresults_shrinkage(void);  //  Function to write results of omega to output window and files
  void outresults_indicator(void);  //  Function to write results of indicator to output window and files
  void outresults_t2(void);         //  Function to write results of t2 to output window and files

  public:

  //____________________________________________________________________________
  //
  // DEFAULT CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector_nigmix(void) : FULLCOND()
    {
    }


  //____________________________________________________________________________
  //
  // CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector_nigmix(MCMCoptions * o, vector<FULLCOND_const*> & p,
                         DISTRIBUTION * d,const ST::string & ti,
                         const ST::string & fp, const ST::string & fr,
                         const vector<unsigned long> & indicators,
                         const vector<double> & v0, const vector<double> & v1,
                         const vector<double> & t2s,
                         const vector<double> & at2, const vector<double> & bt2,
                         const vector<double> & omegas,
                         const vector<double> & aomega, const vector<double> & bomega,
                         const vector<bool> & omegaf,
//                         const datamatrix start_data,
                         const vector<bool> & omegaad,
                         const vector<unsigned> & ct,
                         const unsigned & c);

  //____________________________________________________________________________
  //
  // COPY CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector_nigmix(const FULLCOND_variance_nonp_vector_nigmix & t);


  //____________________________________________________________________________
  //
  // OVERLOADED ASSIGNMENT OPERATOR
  //____________________________________________________________________________

  const FULLCOND_variance_nonp_vector_nigmix & operator=(const FULLCOND_variance_nonp_vector_nigmix & t);


  // Pointer auf das shrinkage-Parameter Fullcond-Objekt
  FULLCOND * get_shrinkagepointer();


  void get_samples(const ST::string & filename, const unsigned & step = 1) const;

  void get_startvalues(void);

  //____________________________________________________________________________
  //
  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //____________________________________________________________________________

  void update(void);


  //____________________________________________________________________________
  //
  // FUNCTION: outresults
  // TASK: - write results to output window and files
  //____________________________________________________________________________

  void outresults(void);

  //____________________________________________________________________________
  //
  // FUNCTION: outoptions
  // TASK: - write options to output window
  //____________________________________________________________________________

  void outoptions(void);


  //____________________________________________________________________________
  //
  // FUNCTION: reset
  // TASK: resets all parameters
  //____________________________________________________________________________

  void reset(void)
    {
    FULLCOND::reset();
    setbeta(beta.rows(),1,0.1);
    }


  //____________________________________________________________________________
  //
  // DESTRUCTOR
  //____________________________________________________________________________

  ~FULLCOND_variance_nonp_vector_nigmix() {}

  }; // end: class FULLCOND_variance_nonp_vector_nigmix



} // end: namespace MCMC

#endif

