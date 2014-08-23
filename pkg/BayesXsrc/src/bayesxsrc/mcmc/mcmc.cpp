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



#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"

#endif


#include "mcmc.h"
#include "clstring.h"

#include <iostream>
using std::flush;

namespace MCMC
{

//------------------------------------------------------------------------------
//-------------------------- CLASS: MCMCoptions --------------------------------
//------------------------------------------------------------------------------


MCMCoptions::MCMCoptions(void)
  {
  iterations = 52000;
  burnin = 2000;
  step = 50;
  nrbetween = 1000;
  nrout = 100;
  nriter = 0;
  samplesize = 0;
  logout = &cout;
  }


MCMCoptions::MCMCoptions(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * abp,
#endif
const unsigned & it,const unsigned & bu,
                         const unsigned & st, ostream * lo,
                         const double & l1,const double & l2)
  {
  iterations = it;
  burnin = bu;
  step = st;
  level1 = l1;
  level2 = l2;
  assert(iterations > 0);
  assert(burnin < iterations);
  assert(step < iterations);

  nrbetween = (iterations-burnin)/3;
  if (nrbetween < 1000)
    nrbetween = 1000;
  nrout = 100;
  nriter = 0;
  samplesize = 0;
  logout = lo;

  (*logout) << flush;
#if defined(BORLAND_OUTPUT_WINDOW)

#elif defined(JAVA_OUTPUT_WINDOW)
adminb_p = abp;
#else
//  if (logout->fail())
//    logout = &cout;
#endif
  }


MCMCoptions::MCMCoptions(const MCMCoptions & o)
  {
  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = o.adminb_p;
  #endif
  iterations = o.iterations;
  burnin = o.burnin;
  step = o.step;
  level1 = o.level1;
  level2 = o.level2;
  nrbetween = o.nrbetween;
  nrout = o.nrout;
  nriter = o.nriter;
  samplesize = o.samplesize;
  logout = o.logout;
  }


const MCMCoptions & MCMCoptions::operator=(const MCMCoptions & o)
  {
  if (this == &o)
    return *this;
  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = o.adminb_p;
  #endif
  iterations = o.iterations;
  burnin = o.burnin;
  step = o.step;
  level1 = o.level1;
  level2 = o.level2;
  nrbetween = o.nrbetween;
  nrout = o.nrout;
  nriter = o.nriter;
  samplesize = o.samplesize;
  logout = o.logout;
  return *this;
  }


void MCMCoptions::out(const ST::string & s,bool thick,bool italic,
                      unsigned size,int r,int g, int b)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  ST::string sh = s;
  sh = sh.replaceallsigns('\n',' ');
  if (!Frame->suppoutput)
    Results->ResultsRichEdit->Lines->Append(sh.strtochar());
 if (!(logout->fail()))
    (*logout) << s << flush;
#elif defined(JAVA_OUTPUT_WINDOW)

  ST::string sh = s;
  sh = sh.replaceallsigns('\n',' ');
  sh = sh+"\n";

  if (!adminb_p->get_suppressoutput())
    adminb_p->Java->CallVoidMethod(adminb_p->BayesX_obj, adminb_p->javaoutput,
    adminb_p->Java->NewStringUTF(sh.strtochar()),
    thick, italic, size,r,g,b);

  if (!(logout->fail()))
    (*logout) << s << flush;
#else
  cout << s << flush;
  if (!(logout->fail()))
   (*logout) << s << flush;
#endif
  }


void MCMCoptions::outerror(const ST::string & s)
  {
  out(s,true,true,12,255,0,0);
  }



void MCMCoptions::outoptions(void)
  {

  out("GENERAL OPTIONS:\n",true);
  out("\n");
  out("  Number of iterations:  " + ST::inttostring(iterations) + "\n");
  out("  Burn-in period:        " + ST::inttostring(burnin)+ "\n");
  out("  Thinning parameter:    " + ST::inttostring(step)+ "\n");
  out("\n");
  }


void MCMCoptions::update(void)
  {

  nriter++;

  if (nriter % nrout == 0 || nriter == 1)
    {
    out("  ITERATION: " + ST::inttostring(nriter) + "\n");
    }


  if( (nriter > burnin) && ((nriter-burnin-1) % step == 0) )
      samplesize++;
  }


void MCMCoptions::update_bootstrap(void)
  {
  if(nriter==0)
    {
    burnin=0;
    step=1;
    }
  nriter++;
  //if(nriter > burnin+1 || nriter==1)
    samplesize++;
  }


unsigned MCMCoptions::compute_samplesize(void)
  {
  return 1+(iterations-burnin-1)/step;
  }

} // end: namespace MCMC



