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


#include "GENERAL_OPTIONS.h"

using std::flush;

namespace MCMC
{

//------------------------------------------------------------------------------
//-------------------------- CLASS: GENERAL_OPTIONS ----------------------------
//------------------------------------------------------------------------------


GENERAL_OPTIONS::GENERAL_OPTIONS(void)
  {
  iterations = 22000;
  burnin = 2000;
  step = 20;
  nrbetween = 1000;
  nrout = 1000;
  nriter = 0;
  samplesize = 0;
  logout = &cout;
  set_level1(95);
  set_level2(80);
  saveestimation = false;
  }


GENERAL_OPTIONS::GENERAL_OPTIONS(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * abp,
#endif
const unsigned & it,const unsigned & bu,
                         const unsigned & st, const bool & sa, ostream * lo,
                         const double & l1,const double & l2)
  {
  iterations = it;
  burnin = bu;
  step = st;
  set_level1(l1);
  set_level2(l2);

  nrbetween = (iterations-burnin)/3;
  if (nrbetween < 1000)
    nrbetween = 1000;
  nrout = 1000;
  nriter = 0;
  samplesize = 0;
  logout = lo;
  saveestimation = sa;

  (*logout) << flush;
#if defined(BORLAND_OUTPUT_WINDOW)

#elif defined(JAVA_OUTPUT_WINDOW)
adminb_p = abp;
#else
  if (logout->fail())
    logout = &cout;
#endif
  }


GENERAL_OPTIONS::GENERAL_OPTIONS(const GENERAL_OPTIONS & o)
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
  lower1 = o.lower1;
  lower2 = o.lower2;
  upper1 = o.upper1;
  upper2 = o.upper2;
  saveestimation = o.saveestimation;
  }


const GENERAL_OPTIONS & GENERAL_OPTIONS::operator=(const GENERAL_OPTIONS & o)
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
  lower1 = o.lower1;
  lower2 = o.lower2;
  upper1 = o.upper1;
  upper2 = o.upper2;
  saveestimation = o.saveestimation;
  return *this;
  }


void GENERAL_OPTIONS::out(const ST::string & s,bool thick,bool italic,
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


void GENERAL_OPTIONS::outoptions(void)
  {

  out("GENERAL OPTIONS:\n",true);
  out("\n");
  out("  Number of iterations:  " + ST::inttostring(iterations) + "\n");
  out("  Burn-in period:        " + ST::inttostring(burnin)+ "\n");
  out("  Thinning parameter:    " + ST::inttostring(step)+ "\n");
  if (saveestimation)
    out("  Saveestimation:        enabled\n");
  else
    out("  Saveestimation:        disabled\n");
  out("\n");
  }


void GENERAL_OPTIONS::update(void)
  {

  nriter++;

  if (nriter % nrout == 0 || nriter == 1)
    {
    out("  ITERATION: " + ST::inttostring(nriter) + "\n");
    }

  if( (nriter > burnin) && ((nriter-burnin-1) % step == 0) )
      samplesize++;
  }


unsigned GENERAL_OPTIONS::compute_samplesize(void)
  {
  return 1+(iterations-burnin-1)/step;
  }


void GENERAL_OPTIONS::reset(void)
  {
  nriter = 0;
  samplesize = 0;
  }


void GENERAL_OPTIONS::set_level1(double l1)
  {
  level1 = l1;

  lower1 = (100.0-level1)/2;
  upper2 = 100.0 - lower1;
  }

void GENERAL_OPTIONS::set_level2(double l2)
  {
  level2 = l2;
  lower2 = (100.0-level2)/2;
  upper1 = 100.0 - lower2;
  }


} // end: namespace MCMC



