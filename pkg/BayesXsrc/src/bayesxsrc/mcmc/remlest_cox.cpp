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

#include"remlest.h"
#include"remlest_cox.h"


//------------------------------------------------------------------------------
//----------------- AFT models with smooth error distribution ------------------
//------------------------------------------------------------------------------

  bool remlest::estimate_aft(datamatrix resp, const datamatrix & offset,
                    const datamatrix & weight)
  {
  unsigned i;

  outoptions();
  out("\n");

  for(i=0;i<fullcond.size();i++)
    fullcond[i]->outoptionsreml();

  out("\n");
  out("REML ESTIMATION STARTED\n",true);
  out("\n");

  bool stop = check_pause();
  if (stop)
    return true;


  // implement AFTs here



  out("ESTIMATION RESULTS:\n",true);
  out("\n");

  for(i=1;i<fullcond.size();i++)
    {
//    beta(0,0) += fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],i-1,false,xcut[i],X.cols()+zcut[i-1],0,false,i);
    }
//  beta(0,0) += fullcond[0]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[0],0,0,false,xcut[0],0,0,false,0);

  return false;
  }

//------------------------------------------------------------------------------
//----------------- AFT models with smooth error distribution ------------------
//---------------------------- fixed effects only ------------------------------
//------------------------------------------------------------------------------


  bool remlest::estimate_aft_glm(datamatrix resp, const datamatrix & offset,
                    const datamatrix & weight)
  {
  unsigned i;

  outoptions();
  out("\n");

  for(i=0;i<fullcond.size();i++)
    fullcond[i]->outoptionsreml();

  out("\n");
  out("REML ESTIMATION STARTED\n",true);
  out("\n");

  bool stop = check_pause();
  if (stop)
    return true;


  // implement AFTs here



  out("ESTIMATION RESULTS:\n",true);
  out("\n");

//  beta(0,0) += fullcond[0]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[0],0,0,false,xcut[0],0,0,false,0);

  return false;
  }

