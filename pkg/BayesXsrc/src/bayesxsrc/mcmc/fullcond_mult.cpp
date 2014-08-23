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



#include "fullcond_mult.h"

namespace MCMC
{



FULLCOND_mult::FULLCOND_mult(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_random * rp,
                         FULLCOND_nonp_basis * ba,
                         bool fi,
                         const ST::string & ti,
                         const ST::string & fp, const ST::string & pres,
                         const unsigned & c)
  : FULLCOND(o,datamatrix(dp->get_nrobs(),1),ti,1,1,fp)
  {

  reffectp = rp;
  basis1p = ba;
  first = fi;
  ttype = re_rw;
  }


FULLCOND_mult::FULLCOND_mult(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_nonp_basis * ba1,
                         FULLCOND_nonp_basis * ba2,
                         bool fi,
                         const ST::string & ti,
                         const ST::string & fp, const ST::string & pres,
                         const unsigned & c)
  : FULLCOND(o,datamatrix(dp->get_nrobs(),1),ti,1,1,fp)
  {

  basis1p = ba1;
  basis2p = ba2;
  first = fi;
  ttype = mrf_rw;
  }

  // COPY CONSTRUCTOR

FULLCOND_mult::FULLCOND_mult(const FULLCOND_mult & fc)
  : FULLCOND(FULLCOND(fc))
  {
  basis1p = fc.basis1p;
  basis2p = fc.basis2p;
  reffectp = fc.reffectp;
  first = fc.first;
  ttype = fc.ttype;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const FULLCOND_mult & FULLCOND_mult::operator=(const FULLCOND_mult & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND::operator=(FULLCOND(fc));
  basis1p = fc.basis1p;
  basis2p = fc.basis2p;
  reffectp = fc.reffectp;
  first = fc.first;
  ttype = fc.ttype;
  return *this;
  }

void FULLCOND_mult::update(void)
  {

  vector<ST::string> enames;

  if (ttype==re_rw)
    {
    if (first)
      {
      basis1p->get_effectmatrix(data,enames,0,0,MCMC::fvar_current);
      reffectp->init_data_varcoeff(data,0);
      }
    else
      {
      reffectp->get_effectmatrix(data,enames,0,0,MCMC::fvar_current);
      basis1p->init_data_varcoeff(data,1);
      }
    }
  else if (ttype == mrf_rw)
    {
    if (first)
      {
      basis2p->get_effectmatrix(data,enames,0,0,MCMC::fvar_current);
      basis1p->init_data_varcoeff(data,0);
      }
    else
      {
      basis1p->get_effectmatrix(data,enames,0,0,MCMC::fvar_current);
      basis2p->init_data_varcoeff(data,1);
      }
    }

  }

//  void update_linpred(const bool & add)
//    {
//    }



bool FULLCOND_mult::posteriormode(void)
  {
  vector<ST::string> enames;

  if (ttype==re_rw)
    {
    if (first)
      {
      basis1p->get_effectmatrix(data,enames,0,0,MCMC::fvar_current);
      reffectp->init_data_varcoeff(data,0);
      }
    else
      {
      reffectp->get_effectmatrix(data,enames,0,0,MCMC::fvar_current);
      basis1p->init_data_varcoeff(data,1);
      }
    }
  else if (ttype==mrf_rw)
    {
    if (first)
      {
      basis2p->get_effectmatrix(data,enames,0,0,MCMC::fvar_current);
      basis1p->init_data_varcoeff(data,0);
      }
    else
      {
      basis1p->get_effectmatrix(data,enames,0,0,MCMC::fvar_current);
      basis2p->init_data_varcoeff(data,1);
      }
    }

  return true;
  }

/*
bool FULLCOND_mult::posteriormode_converged(const unsigned & itnr)
  {


  }
*/


void FULLCOND_mult::outresults(void)
  {

  }

void FULLCOND_mult::get_effectmatrix(datamatrix & e,vector<ST::string> & enames,
                        unsigned be, unsigned en,effecttype t)
  {

  }


// unsigned FULLCOND_mult::get_nreffects(effecttype t)
//  {


//  }


void FULLCOND_mult::outoptions(void)
  {


  }

ST::string FULLCOND_mult::getinfo(void)
  {
    return "";
  }


void FULLCOND_mult::init_name(const ST::string & na)
  {

  }

void FULLCOND_mult::init_names(const vector<ST::string> & na)
  {


  }

void FULLCOND_mult::init_priorassumptions(const ST::string & na)
  {


  }

  // FUNCTION: reset
  // TASK: resets all parameters

void FULLCOND_mult::reset(void)
  {

  }


  } // end: namespace MCMC
