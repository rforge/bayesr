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



#if !defined (VARIANCENONP_VECTOR_INCLUDED)
#define VARIANCENONP_VECTOR_INCLUDED

#include"../export_type.h"
#include"mcmc_const.h"

namespace MCMC
{


class __EXPORT_TYPE FULLCOND_variance_nonp_vector : public FULLCOND
  {

  protected:

  vector<double> tau;               // Varianceparameters
  vector<double> lambda;            // Inverse Varianceparameter: lambda=1/tau^2

  bool update_sigma2;

  unsigned column;                  //  Category for for fullcond if multivariate response

  ST::string pathresults;           //  File path for storing sampled parameters

  vector<FULLCOND_const *> Cp;

  DISTRIBUTION * distrp;

  FULLCOND fc_shrinkage;            //  Fullcondobjekt shrinkageparameter
  bool shrinkage_fix;               //  Shrinkageparameter fix
  bool shrinkage_adaptive;          //  variancespecific shrinkageparameter
  vector<double> a_shrinkage;       //  Hyperparameter for Shrinkageparameter
  vector<double> b_shrinkage;       //  Hyperparameter for Shrinkageparameter
  vector<double> weight;            //  Weights for Shrinkage
  datamatrix startdata;             //  Matrix with Starting values and hyperparameters

  double lassosum;                  //  sum(beta^2/tau^2)
  double ridgesum;                  //  sum(beta^2/tau^2)

  vector<unsigned> cut;             //  Blocks of regression coefficients
  bool is_ridge;                    //  The Components indicates if "true" the L2-penalty
                                    //  and if "false" the L1-penalty is used

//  datamatrix variances;              // current values of the variances

  void outresults_shrinkage(void);  //  Function to write results to output window and files

  public:

  //____________________________________________________________________________
  //
  // DEFAULT CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector(void) : FULLCOND()
    {
    }


  //____________________________________________________________________________
  //
  // CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector(MCMCoptions * o, vector<FULLCOND_const*> & p,
                         DISTRIBUTION * d,const ST::string & ti,
                         const ST::string & fp, const ST::string & fr,
                         const vector<double> & shrinkagestart, const vector<double> & ashrinkage,
                         const vector<double> & bshrinkage, const vector<bool> & shrinkagefix,
                         const vector<double> & shrinkageweight,
//                         const datamatrix start_data,
                         const vector<bool> & shrinkageadaptive,
                         const bool & isridge, const vector<unsigned> & ct,
                         const unsigned & c);


  //____________________________________________________________________________
  //
  // COPY CONSTRUCTOR
  //____________________________________________________________________________

  FULLCOND_variance_nonp_vector(const FULLCOND_variance_nonp_vector & t);


  //____________________________________________________________________________
  //
  // OVERLOADED ASSIGNMENT OPERATOR
  //____________________________________________________________________________

  const FULLCOND_variance_nonp_vector & operator=(const FULLCOND_variance_nonp_vector & t);


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

  ~FULLCOND_variance_nonp_vector() {}

  }; // end: class FULLCOND_variance_nonp_vector



} // end: namespace MCMC

#endif

