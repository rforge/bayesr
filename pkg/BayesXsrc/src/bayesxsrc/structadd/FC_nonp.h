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



#if !defined (FCNONPINCLUDED)

#define FCNONPINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"FC.h"
#include"design.h"
#include"MASTER_obj.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_nonp -----------------------------------
//------------------------------------------------------------------------------

enum sampletype {unconstrained,increasing,decreasing};

class __EXPORT_TYPE FC_nonp  : public FC
  {

  protected:

  MASTER_OBJ * masterp;
  unsigned equationnr;

  FC fsample;
  bool samplef;

  FC paramsample;

  bool derivative;
  bool samplederivative;
  FC derivativesample;


  bool IWLS;


  bool orthogonal;
  datamatrix acuteparam;

  sampletype stype;

  DISTR * likep;                             // Pointer to DISTR obejct

  datamatrix betaold;
  datamatrix betadiff;

  double s2;
  void perform_centering(void);
  void centerparam(void);
  void centerparam_weight(void);
  void centerparam_sum2(double & s2);
  void centerparam_sample(void);
  void centerparam_sample_var(void);

  bool computemeaneffect;
  double meaneffectconstant;
  FC meaneffect_sample;


  public:

  DESIGN * designp;                          // Pointer to design object


  datamatrix param;                          // Parameters, beta stores hatf
  datamatrix paramlin;
  datamatrix paramold;
  datamatrix paramhelp;
  datamatrix parammode;
  double paramKparam;

  datamatrix partres;                        // sum of partial residuals

  double lambda;
  double tau2;

  //---------------------------- centering -------------------------------------

  datamatrix Vcenter;
  datamatrix Vcentert;
  datamatrix Wcenter;
  datamatrix Ucenter;
  datamatrix Utc;
  datamatrix ccenter;
  datamatrix helpcenter;

  //--------------------------- importance measures ---------------------------

  bool imeasures;


  void get_linparam(void);

  void initialize_center(void);

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_nonp(void);

  // CONSTRUCTOR
  // o    : pointer to GENERAL_OPTIONS object
  // t    : title of the full conditional (for example "fixed effects")
  // fp   : file path for storing sampled parameters

  FC_nonp(MASTER_OBJ * mp,unsigned & enr ,GENERAL_OPTIONS * o,DISTR * lp,
          const ST::string & t,
           const ST::string & fp,DESIGN * dp,vector<ST::string> & op,
             vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  FC_nonp(const FC_nonp & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_nonp & operator=(const FC_nonp & m);

  // DESTRUCTOR


  ~FC_nonp()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  void update_gaussian_transform(void);
  void update_gaussian(void);
  void update_IWLS(void);
  void update_isotonic(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);

  bool posteriormode_transform(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void outoptions(void);

  // FUNCTION: outgraphs
  // TASK: writes batch files for STATA and R for visualizing results

  void outgraphs(ofstream & out_stata, ofstream & out_R,const ST::string & path);


  double kernel_density(const double & x, const double & h);

  double compute_importancemeasure(bool absolute);
  double compute_importancemeasure_discrete(bool absolute);

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void outresults(ofstream & out_stata, ofstream & out_R,
                  const ST::string & pathresults);

  // FUNCTION: outresults_derivative
  // TASK: writes estimated first derivatives

  void outresults_derivative(ofstream & out_stata, ofstream & out_R,
                        const ST::string & pathresults);

  void outbasis_R(const ST::string & pathbasis);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void);

  void compute_autocorr_all(const ST::string & path,
                                      unsigned lag, ofstream & outg) const;

  void get_samples(const ST::string & filename,ofstream & outg) const;

  void check_errors(void);

  void get_effect(datamatrix & effect);


  };


} // end: namespace MCMC

#endif


