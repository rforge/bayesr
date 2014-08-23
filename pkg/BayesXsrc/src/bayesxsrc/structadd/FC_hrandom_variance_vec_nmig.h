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


#if !defined (FChrandomVARIANCEVECNMIGINCLUDED)

#define FChrandomVARIANCEVECNMIGINCLUDED



#include"FC_hrandom_variance_vec.h"


namespace MCMC
{

//------------------------------------------------------------------------------
//---------------------- CLASS: FC_hrandom_variance_vec_nmig -------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_hrandom_variance_vec_nmig
      : public FC_hrandom_variance_vec
  {

  protected:

  FC FC_delta;
  FC FC_omega;
  FC FC_Q;

  double abeta;
  double bbeta;
  double r;
  double v;
  double aQ;
  double bQ;
  long regiterates;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_hrandom_variance_vec_nmig(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_hrandom_variance_vec_nmig(MASTER_OBJ * mp,unsigned & enr,GENERAL_OPTIONS * o,DISTR * lp,
                          DISTR * lpRE, const ST::string & t,
                          const ST::string & fp,DESIGN * dp,
                          FC_nonp * FCn,vector<ST::string> & op,
                          vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_hrandom_variance_vec_nmig(const FC_hrandom_variance_vec_nmig & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_hrandom_variance_vec_nmig & operator=(
   const FC_hrandom_variance_vec_nmig & m);

  // DESTRUCTOR

  ~FC_hrandom_variance_vec_nmig()
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


  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  void get_samples(const ST::string & filename,ofstream & outg) const;

  void compute_autocorr_all(const ST::string & path,
                                      unsigned lag, ofstream & outg) const;


  };


} // end: namespace MCMC

#endif


