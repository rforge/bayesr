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



#if !defined (FULLCOND_MERROR_INCLUDED)

#define FULLCOND_MERROR_INCLUDED

#include"../export_type.h"
#include "fullcond.h"
#include "mcmc_nonpbasis.h"
#include "fullcond_nonp_gaussian.h"
#include "spline_basis.h"
#include "mcmc_nonp.h"
#include "mcmc.h"

namespace MCMC
{
class __EXPORT_TYPE fullcond_merror : public FULLCOND
  {

  protected:

  FULLCOND_nonp_basis * designp;  // Pointer wird für IWLS-proposal benötigt
//  FULLCOND_nonp * designp;          // Pointer wird für Conditional prior proposal benötigt
  DISTRIBUTION * likep;

// BEGIN: merror
  spline_basis * splinep;
  bool varcoeff;

  double maxx;
  double minx;

  datamatrix meandata;      // mean of the observed covariate values
  datamatrix old;           // sampled values from the previous iteration

  datamatrix currentspline; // stores current f(x) (i.e. f(x^{old})
  datamatrix diffspline;    // stores f(x^{prop})-f(x^{old})
  datamatrix logfcold;      // full conditional evaluated for old values
  datamatrix logfcnew;      // full conditional evaluated for propsed values

  int merror;               // counts number of replicated measurements for each covariate value

  statmatrix<int> index;

  FULLCOND fc_merrorvar;    // full conditional for the variance of the measurement error
  FULLCOND fc_ximu;         // full conditional for the expectation of the true covariate values
  FULLCOND fc_xivar;        // full conditional for the variance of the true covariate values

  unsigned generrcount;      // counts the number of generated values out of range
  unsigned generrtrial;

  bool discretize;           // generation of rounded covariate values
  unsigned digits;           // no. of digits for rounding
  unsigned nbeta;            // no of observations with nonzero weight

  ST::string pathresults;
// END: merror

// BEGIN: Susi

  // SUSI: Add help fullcond object

//  FULLCOND whatsoever;
  FULLCOND fc_bias;

  unsigned drows;
  unsigned dcols;

  double biasmean;
  double biasvar;
  double sigma12;

  datamatrix effmod;

  datamatrix P;                     // Präzisionsmatrix der wahren Werte
  //envmatdouble precenv;           // envmatrix zum Speichern der Präzisionsmatrixx der wahren Werte
                                    // wird nur für IWLS-proposal benötigt

  datamatrix rhs;                   // datamatrix zum Speichern von Z'Sigma^{-1}X + Omega^{-1}mu_xi
 //datamatrix betahelp1;            // wird nur für IWLS-proposal benötigt

  datamatrix mmu;                   // P^{-1}*rhs = Erwartungswert der wahren Werte;

  datamatrix linold;
  datamatrix linnew;
  datamatrix diff;

  datamatrix proposalold;
  datamatrix proposal;

  datamatrix PABn;                  // sqrt(P_ab^-1)

  datamatrix PABl;                  // P_left matrices
  datamatrix PABr;                  // P_right matrices

  unsigned minmerror;               // Minimum Blocksize
  unsigned maxmerror;               // Maximum Blocksize

  vector<unsigned> matquant;         // matquant[size-min] gives the number of
                                     // blocks for blocksize 'size'

  datamatrix merror_random;
  datamatrix randnorm;

  datamatrix randnormal;

  datamatrix xi;

// END: Susi

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  fullcond_merror(void) : FULLCOND()
    {
    }

  // CONSTRUCTOR : Susi (Measurement error in the interaction variable of a VCM)
  // o    : pointer to MCMCoptions object
  // t    : title of the full conditional (for example "fixed effects")
  //        (i.e. number of categories of the response variable)
  // fp   : file path for storing sampled parameters

  fullcond_merror(MCMCoptions * o, /*FULLCOND_nonp * p,*/ FULLCOND_nonp_basis * p, DISTRIBUTION * dp,
           const datamatrix & d, const datamatrix & em, const ST::string & t, const ST::string & fp,
           const double & mvar1, const double & mvar2, const double & arvar,
           const double & arpar1, const double & arpar2, const double & bmean,
           const double & bvar);

  // BEGIN: merror
  // CONSTRUCTOR : Thomas (Measurement error in a nonparametric effect)
  fullcond_merror(MCMCoptions * o, spline_basis * p, DISTRIBUTION * dp,
           const datamatrix & d, const ST::string & t, const ST::string & fp,
           const ST::string & pres, const double & lk, const double & uk,
           const double & mvar, const bool & disc, const int & dig,
           const unsigned & nb);
// END: merror

  // COPY CONSTRUCTOR

  fullcond_merror(const fullcond_merror & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const fullcond_merror & operator=(const fullcond_merror & m);

  // DESTRUCTOR

  ~fullcond_merror()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  void compute_mu(datamatrix & muexi,
      const unsigned & blocks,const unsigned & a,const unsigned & b);

  void compute_proposal(const datamatrix & xi, const unsigned & blocks,
                               const unsigned & a,const unsigned b);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void posteriormode_set_beta_mode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

//  vector<ST::string> & get_results_latex(void);

  };

} // end: namespace MCMC

#endif

