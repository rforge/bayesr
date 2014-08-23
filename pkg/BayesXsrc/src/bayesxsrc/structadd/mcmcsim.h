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



#if !defined (MCMCsim_INCLUDED)

#define MCMCsim_INCLUDED

#include"../export_type.h"
#include"GENERAL_OPTIONS.h"
#include"distr.h"
#include"FC.h"

namespace MCMC
{


class __EXPORT_TYPE equation
  {

  
  protected:

  public:

  int hlevel;
  int equationnr;
  ST::string equationtype;

  unsigned nrfc;

  ST::string header;
  ST::string paths;

  DISTR * distrp;
  ST::string pathd;

  vector<FC*> FCpointer;
  vector<ST::string> FCpaths;

  // DEFAULT CONSTRUCTOR

  equation(void);

  // CONSTRUCTOR1

  equation(int enr, int hl,ST::string t);

  // CONSTRUCTOR2

  equation(unsigned enr, const ST::string & h, DISTR * dp,
           const vector<FC*> fcp,
           const ST::string & pd, const vector<ST::string> & ps);

  // COPY CONSTRUCTOR

  equation(const equation & s);

  // OVERLOADED ASSIGNMENT CONSTRUCTOR

  const equation & operator=(const equation & s);

  void add_FC(FC * FCp,const ST::string & p);

  // DESTRUCTOR

  ~equation() {}

  };


class __EXPORT_TYPE MCMCsim
  {

  protected:

  GENERAL_OPTIONS * genoptions;

  vector<equation> equations;

  unsigned maxiterations;             // for posteriormode, maximum number of
                                      // iterations, default = 1000

  public:

  // DEFAULT CONSTRUCTOR

  MCMCsim(void)
    {
    }

  // CONSTRUCTOR
  // TASK: initializes the MCMC simulation object with general MCMC options 'go'
  //       a vector of equations 'equ'

  MCMCsim(GENERAL_OPTIONS * go,vector<equation> & equ,unsigned maxit=1000);

  // COPY CONSTRUCTOR

  MCMCsim(const MCMCsim & s);

  // OVERLOADED ASSIGNMENT CONSTRUCTOR

  const MCMCsim & operator=(const MCMCsim & s);

  // FUNCTION: simulate
  // TASK: runs a MCMC simulation
  //       returns true, if simulation error or user break occured

  bool simulate(ST::string & pathgraphs, const int & seed,
                const bool & computemode=true);

  bool posteriormode(ST::string & pathgraphs,const bool & presim=false);

  void out_effects(const vector<ST::string> & paths);


  // FUNCTION: autocorr
  // TASK: computes autocorrelations for all samples parameters
  //      (i.e. for all beta matrices) and stores the result in datamatrix
  //      'cmat'

  void autocorr(const unsigned & lag,datamatrix & cmat);

  // FUNCTION: autocorr
  // TASK: computes autocorrelations for all samples parameters
  //      (i.e. for all beta matrices) and stores the result in file 'path'

  void autocorr(const unsigned & lag,ST::string & pathgraphs);

  // FUNCTION: compute_nrpar
  // TASK: computes the total number of parameters

  unsigned compute_nrpar(void);

  // FUNCTION: get_samples
  // TASK: stores sampled parameters of all full conditionals in ASCII format
  //       for each full conditional one file will be created with filename
  //       'path' + title of the full conditional + "_sample.raw"

  void get_samples(ST::string & pathgraphs
  #if defined(JAVA_OUTPUT_WINDOW)
  , vector<ST::string> & newc
  #endif
  );

  // DESTRUCTOR

  ~MCMCsim() {}

  };



} // end: namespace MCMC

#endif
