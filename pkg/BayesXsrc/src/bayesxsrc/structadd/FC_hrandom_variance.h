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



#if !defined (FChrandomVARIANCEINCLUDED)

#define FChrandomVARIANCEINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC_nonp_variance.h"
#include"design.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_hrandom_variance -----------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_hrandom_variance  : public FC_nonp_variance
  {

  protected:

  DISTR * likepRE;

  bool mult;

  double compute_quadform(void);

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_hrandom_variance(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_hrandom_variance(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp, DISTR * lpRE,
                      const ST::string & t, const ST::string & fp,DESIGN * dp,
                      FC_nonp * FCn,vector<ST::string> & op,
                      vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_hrandom_variance(const FC_hrandom_variance & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_hrandom_variance & operator=(const FC_hrandom_variance & m);

  // DESTRUCTOR

  ~FC_hrandom_variance()
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


  void outresults(ofstream & out_stata, ofstream & out_R,
                         const ST::string & pathresults);


  };


//------------------------------------------------------------------------------
//---------------------- CLASS: FC_hrandom_variance_ssvs -----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FC_hrandom_variance_ssvs
      : public FC_hrandom_variance
  {

  protected:

  FC FC_delta;
  FC FC_omega;

  double abeta;
  double bbeta;
  double r;
  long regiterates;

  datamatrix pen;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_hrandom_variance_ssvs(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_hrandom_variance_ssvs(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
                          DISTR * lpRE, const ST::string & t,
                          const ST::string & fp,DESIGN * dp,
                          FC_nonp * FCn,vector<ST::string> & op,
                          vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_hrandom_variance_ssvs(const FC_hrandom_variance_ssvs & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_hrandom_variance_ssvs & operator=(
   const FC_hrandom_variance_ssvs & m);

  // DESTRUCTOR

  ~FC_hrandom_variance_ssvs()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata,ofstream & out_R,
                  const ST::string & pathresults);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  void get_samples(const ST::string & filename,ofstream & outg) const;

  void compute_autocorr_all(const ST::string & path,
                                      unsigned lag, ofstream & outg) const;


  };


} // end: namespace MCMC

#endif


