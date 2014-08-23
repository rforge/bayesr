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



#include"../export_type.h"
#if !defined (MCMCsimulate_INCLUDED)

#define MCMCsimulate_INCLUDED

#include"../export_type.h"
#include"mcmc.h"
#include"distribution.h"
#include"fullcond.h"
#include"baseline.h"
#include"multibaseline.h"

namespace MCMC
{

class __EXPORT_TYPE MCMCsimulate
  {

  protected:

  vector<MCMCoptions *> genoptions_mult;  // Pointer to general MCMC options
  vector<DISTRIBUTION *> likep_mult;      // Pointer to distribution object

  vector<FULLCOND*> fullcondp;        // Vector of pointers to full conditionals

  bool likepexisting;


  vector<unsigned> begin;
  vector<unsigned> end;

  void set_center(DISTRIBUTION * lp,vector<FULLCOND *> fp,
                 unsigned beg, unsigned en);

  unsigned compute_nrpar(void);

  bool checkerrors(DISTRIBUTION * lp,vector<FULLCOND *> fp,
                   unsigned beg, unsigned en);


  public:


  // DEFAULT CONSTRUCTOR

  MCMCsimulate(void)
    {
    }

  // CONSTRUCTOR1
  // TASK: initializes the MCMC simulation object with general MCMC options 'go'
  //       distribuiton object dp and a vector of full conditionals 'fc'

  MCMCsimulate(MCMCoptions * go,DISTRIBUTION * dp,vector<FULLCOND*> & fc);

  // CONSTRUCTOR2
  // TASK: initializes the MCMC simulation object with general MCMC options 'go'
  //       distribution dp and full conditional 'fc'

  MCMCsimulate(MCMCoptions * go,DISTRIBUTION * dp,FULLCOND* fc);

  // CONSTRUCTOR3
  // TASK: initializes the MCMC simulation object with general MCMC options 'go'
  //       and a vector of full conditionals 'fc', NO distribution object

  MCMCsimulate(MCMCoptions * go,vector<FULLCOND*> & fc);

  // CONSTRUCTOR4
  // TASK: initializes the MCMC simulation object with general MCMC
  //       options 'gom' distribuiton objects dp and a vector of vectors of
  //       full conditionals 'fc'

  MCMCsimulate(vector<MCMCoptions *> go,vector<DISTRIBUTION *> dp,
               vector<FULLCOND*>  & fc,vector<unsigned> & be,
                                       vector<unsigned> & en);

  // COPY CONSTRUCTOR

  MCMCsimulate(const MCMCsimulate & s);

  // OVERLOADED ASSIGNMENT CONSTRUCTOR

  const MCMCsimulate & operator=(const MCMCsimulate & s);

  // FUNCTION: simulate
  // TASK: runs a MCMC simulation
  //       returns true, if simulation error or user break occured

  bool simulate(const vector<ST::string> & header, const int & seed,
                const bool & computemode=true);

  bool posteriormode(const vector<ST::string> & header,
                     const bool & presim=false);

  void out_effects(const vector<ST::string> & paths);

  void make_graphics(const vector<ST::string> & title,
                     const vector<ST::string> & path_batch,
                     const vector<ST::string> & path_tex,
                     const vector<ST::string> & path_splus,
                     const vector<ST::string> & path_stata);

  void make_options(ofstream & o,const unsigned & nr);

  void make_predictor(ofstream & outtex,const unsigned & nr);

  void make_model(ofstream & o,const unsigned & nr);

  void make_prior(ofstream & o,const unsigned & nr);

  void make_fixed_table(ofstream & o,const unsigned & nr);

  void make_plots(ofstream & outtex,const unsigned nr,
                  const ST::string & path_batch,const ST::string & path_splus,
                  const ST::string & path_stata);

  // FUNCTION: autocorr
  // TASK: computes autocorrelations for all samples parameters
  //      (i.e. for all beta matrices) and stores the result in datamatrix
  //      'cmat'

  void autocorr(const unsigned & lag,datamatrix & cmat);

  // FUNCTION: autocorr
  // TASK: computes autocorrelations for all samples parameters
  //      (i.e. for all beta matrices) and stores the result in file 'path'

  void autocorr(const unsigned & lag,const ST::string & path);

  // FUNCTION: get_samples
  // TASK: stores sampled parameters of all full conditionals in ASCII format
  //       for each full conditional one file will be created with filename
  //       'path' + title of the full conditional + "_sample.raw"

  void get_samples(
  #if defined(JAVA_OUTPUT_WINDOW)
  vector<ST::string> & newc,
  #endif
  const ST::string & path,const unsigned & step=1);

  void setflags(const bitset<flagnr> & newflags);

  // DESTRUCTOR

  ~MCMCsimulate() {}

  };


//------------------------------------------------------------------------------
//------------------ FUNCTIONS FOR COMPARING DIFFERENT MODELS ------------------
//------------------------------------------------------------------------------

// FUNCTION: compare
// TASK: compares estimation results stored in 'files'
//       assumes that each file in 'files' has the same structure
//       assumes that the first line contains variable names
//       computes for each column the relative differences compared to the
//       first file in 'files' and writes the result to 'out'

// RELATIVE DIFFERENCE BETWEEN TWO VECTORS:
// a = a_1,...a_n
// b = b_1,...,b_n
// RELDIFF = sqrt{ a_1-b_1)^2 + ... + (a_n-b_n)^2 } / sqrt{ a_1^2+...+a_n^2 }
//         = ||a-b|| / ||a||
// if a and b are two statmatrix objects (see statmat.h) the norm may be
// computed using member function norm

void compare(const vector<ST::string> & files,ostream & out);


} // end: namespace MCMC

#endif
