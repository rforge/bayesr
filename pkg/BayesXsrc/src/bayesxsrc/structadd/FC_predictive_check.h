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



#if !defined (FCpredictivecheckINCLUDED)

#define FCpredictivecheckINCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"sparsemat.h"

#include"Random.h"
#include"../values.h"
#include<fstream>
#include<vector>
#include<bitset>
#include"GENERAL_OPTIONS.h"
#include"distr.h"
#include"clstring.h"
#include<cmath>

namespace MCMC
{

using std::vector;
using std::bitset;

//------------------------------------------------------------------------------
//------------------------- CLASS: FC_predictive_check -------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FC_predictive_check   : public FC
  {

  protected:


  datamatrix sampled_responses;

  DISTR * likep;

  datamatrix designmatrix;
  vector<ST::string> varnames;

  void get_predictor(void);


  public:

  // DEFAULT CONSTRUCTOR

  FC_predictive_check(void);

  // CONSTRUCTOR

  FC_predictive_check(GENERAL_OPTIONS * o,DISTR * lp,const ST::string & t,
     const ST::string & fp,datamatrix & dm, vector<ST::string> & dn);

  // COPY CONSTRUCTOR

  FC_predictive_check(const FC_predictive_check & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_predictive_check & operator=(const FC_predictive_check & m);

  // DESTRUCTOR

  ~FC_predictive_check()
    {
    }


  void update(void);

  bool posteriormode(void);

  void outoptions(void);

  void outresults(ofstream & out_stata, ofstream & out_R,
                  const ST::string & pathresults);

  void reset(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  };


} // end: namespace MCMC

#endif


