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



#include "FC_nonp_nongaussian.h"


//------------------------------------------------------------------------------
//-------- CLASS: FC_nonp_nongaussian implementation of member functions -------
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_nonp_nongaussian::read_options(vector<ST::string> & op,
                                      vector<ST::string> & vn)
  {

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  7       center
  8       map
  9       lambda_re
  10      a_re
  11      b_re
  12      internal_mult
  13      samplemult
  14      constraints
  */

  FC_nonp::read_options(op,vn);

  }


FC_nonp_nongaussian::FC_nonp_nongaussian(void)
  {
  }


FC_nonp_nongaussian::FC_nonp_nongaussian(GENERAL_OPTIONS * o,DISTR * lp,
                 const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,vector<ST::string> & op, vector<ST::string> & vn)
     : FC_nonp(o,lp,t,fp,Dp,op,vn)
  {

  }


FC_nonp_nongaussian::FC_nonp_nongaussian(const FC_nonp_nongaussian & m)
  : FC_nonp(FC_nonp(m))
  {

  }


const FC_nonp_nongaussian & FC_nonp_nongaussian::operator=(
               const FC_nonp_nongaussian & m)
  {

  if (this==&m)
	 return *this;
  FC_nonp::operator=(FC_nonp(m));

  return *this;
  }


void FC_nonp_nongaussian::update(void)
  {
  if ((stype == increasing) || (stype==decreasing))
    {

    }
  else
    {

    update_IWLS();
    }
  }


void FC_nonp_nongaussian::update_IWLS(void)
  {


  transform_beta();

  FC::update();

  }




bool FC_nonp_nongaussian::posteriormode(void)
  {


  }

void FC_nonp_nongaussian::outoptions(void)
  {
  FC_nonp::outoptions();
  }



void FC_nonp_nongaussian::reset(void)
  {

  }


} // end: namespace MCMC



