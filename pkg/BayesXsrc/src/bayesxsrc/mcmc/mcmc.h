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



#if !defined (MCMC_INCLUDED)

#define MCMC_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"statmat_penalty.h"
#include"sparsemat.h"
#include"Random.h"
#include<fstream>
#include<vector>

#if defined(JAVA_OUTPUT_WINDOW)
#include"adminparse_basic.h"
#endif

using std::cout;

namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: MCMCoptions -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE MCMCoptions
  {

  protected:

  unsigned iterations;            // Number of iterations
  unsigned burnin;                // Number of burnin iterations
  unsigned step;                  // Thinning parameter
                                  // (step = 100: every 100th sampled parameter
                                  //  will be stored and used for computing
                                  //  charkteristics of the posterior)

  double level1;                  // level1 for credible intervals
  double level2;                  // level2 for credible intervals

  unsigned nrbetween;
  unsigned nrout;

  unsigned nriter;                // current iteration number
  unsigned samplesize;            // current number of parameters used to
                                  // compute characteristics of the
                                  // posterior

  ostream * logout;               // Pointer to filestream for writing
                                  // output

  public:

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adminb_p;
  #endif


  // DEFAULT CONSTRUCTOR
  // Defines:
  // iterations = 52000
  // burnin = 2000
  // step = 50
  // nrbetween = 1000    (may be changes via set_nrbetween)
  // nrout = 100         (may be changes via set_nrout)
  // nriter = 0
  // samplesize = 0
  // logout = &cout

  MCMCoptions(void);

  // CONSTRUCTOR
  // Defines
  // iterations = it
  // burnin = bu
  // step = st
  // nrbetween = (iterations-burnin)/10
  // if (nrbetween < 1000) nrbetween = 1000
  // nrout = 100
  // nriter = 0
  // samplesize = 0
  // logout = lo

  #if defined(JAVA_OUTPUT_WINDOW)
  MCMCoptions(administrator_basic * abp,
              const unsigned & it,const unsigned & bu,const unsigned & st,
              ostream * lo=&cout,const double & l1=95,const double & l2=80);
  #else
  MCMCoptions(const unsigned & it,const unsigned & bu,const unsigned & st,
              ostream * lo=&cout,const double & l1=95,const double & l2=80);
  #endif

  // COPY CONSTRUCTOR

  MCMCoptions(const MCMCoptions & o);

  // OVERLOADED ASSIGNMENT OPERATOR

  const MCMCoptions & operator=(const MCMCoptions & o);

  // FUNCTION: out
  // TASK: writes results to outputstream or
  //       in Output window if BORLAND_OUTPUT_WINDOW is defined

  void out(const ST::string & s,bool thick=false,bool italic = false,
           unsigned size = 12,int r=0,int g=0, int b=0);

  void outerror(const ST::string & s);


  void set_nrbetween(const unsigned & n)
    {
    nrbetween = n;
    }


  void set_nrout(const unsigned & n)
    {
    nrout = n;
    }

  const unsigned & get_iterations(void)
    {
    return iterations;
    }

  const unsigned & get_burnin(void)
    {
    return burnin;
    }

  const unsigned & get_step(void)
    {
    return step;
    }

  const double & get_level1(void)
    {
    return level1;
    }

  const double & get_level2(void)
    {
    return level2;
    }

  const unsigned & get_nriter(void)
    {
    return nriter;
    }

  const unsigned & get_nrbetween(void)
    {
    return nrbetween;
    }

  const unsigned & get_samplesize(void)
    {
    return samplesize;
    }


  unsigned compute_samplesize(void);


  // FUNCTION: reset

  void reset(void)
    {
    nriter = 0;
    samplesize = 0;
    }

  void update(void);

  void update_bootstrap(void);

  // FUNCTION: outoptions
  // TASK: writes general options (iterations,burnin,step,family) to
  //       output stream

  virtual void outoptions(void);

  // DESTRUCTOR

  ~MCMCoptions() {}


  };


//------------------------------------------------------------------------------
//------------------------ End: CLASS MCMCoptions ------------------------------
//------------------------------------------------------------------------------



} // end: namespace MCMC

#endif
