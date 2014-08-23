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



#if !defined (FCcvINCLUDED)

#define FCcvINCLUDED

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
#include"FC_hrandom.h"
#include"clstring.h"
#include<cmath>

namespace MCMC
{

using std::vector;
using std::bitset;

//------------------------------------------------------------------------------
//------------------------------ CLASS: FC_cv ----------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FC_cv   : public FC
  {

  protected:

  statmatrix<unsigned> ind;

  vector<ST::string> effectvalues;


  unsigned nrcat;

  void get_ind(void);

  datamatrix sampled_etas;
  datamatrix sampled_responses;
  datamatrix sampled_likelihood;

  FC FC_sampled_l;

  DISTR * likep;

  vector<FC_hrandom> * hrandoms;
  unsigned size;

  datamatrix effect;
  datamatrix linpred;


  datamatrix e_score;                                        // energy score
  datamatrix log_score;                                       // log score

  // FUNCTION: compute_energyscore
  // TASK: computes the energy score and stores the result in e_score

  double compute_energyscore(void);

  // FUNCTION: compute_logscore
  // TASK: computes the log score and stores the result in log_score

  double compute_logscore(void);


  public:

  // DEFAULT CONSTRUCTOR

  FC_cv(void);

  // CONSTRUCTOR

  FC_cv(GENERAL_OPTIONS * o,DISTR * lp,const ST::string & t,
     const ST::string & fp,vector<FC_hrandom> * FChs);

  // COPY CONSTRUCTOR

  FC_cv(const FC_cv & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_cv & operator=(const FC_cv & m);

  // DESTRUCTOR

  ~FC_cv()
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


