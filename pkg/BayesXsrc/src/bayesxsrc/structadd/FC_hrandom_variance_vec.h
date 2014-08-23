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



#if !defined (FChrandomVARIANCEVECINCLUDED)

#define FChrandomVARIANCEVECINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC_nonp_variance_vec.h"
#include"design.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- CLASS: FC_hrandom_variance_vec -----------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_hrandom_variance_vec  : public FC_nonp_variance_vec
  {

  protected:

  DISTR * likepRE;

  bool mult;

  double hyperLambda;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_hrandom_variance_vec(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_hrandom_variance_vec(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
                          DISTR * lpRE, const ST::string & t,
                          const ST::string & fp,DESIGN * dp,
                          FC_nonp * FCn,vector<ST::string> & op,
                          vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_hrandom_variance_vec(const FC_hrandom_variance_vec & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_hrandom_variance_vec & operator=(const FC_hrandom_variance_vec & m);

  // DESTRUCTOR

  ~FC_hrandom_variance_vec()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  // void transform_beta(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  };


} // end: namespace MCMC

#endif


