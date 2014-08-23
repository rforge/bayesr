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



#if !defined (FClinearINCLUDED)

#define FClinearINCLUDED

#include"../export_type.h"
#include"../values.h"
#include<fstream>
#include"GENERAL_OPTIONS.h"
#include"distr.h"
#include"clstring.h"
#include"FC.h"
#include"MASTER_obj.h"
#include<cmath>

namespace MCMC
{

//------------------------------------------------------------------------------
//--------------------------- CLASS: FC_linear ---------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_linear  : public FC
  {

  protected:

  int constposition;
  bool constwarning;
  bool center;

  MASTER_OBJ * masterp;
  unsigned equationnr;

  void compute_meaneffect_design(void);

  bool initialize;
  bool IWLS;

  DISTR * likep;                             // Pointer to DISTR obejct
  datamatrix design;                         // Designmatrix
  datamatrix mean_designcols;
  vector<datamatrix> designhelp;             // help vector for constructing the
                                             // designmatrix
  datamatrix meaneffectdesign;


  datamatrix Xt;                             // transposed designmatrix
  datamatrix XWX;
  bool rankXWX_ok;
  datamatrix XWXold;
  datamatrix XWXroot;
  datamatrix residual;
  datamatrix Xtresidual;

  datamatrix betam;
  datamatrix help;
  datamatrix betaold;
  datamatrix betadiff;
  datamatrix mode;
  datamatrix proposal;

  datamatrix linold;
  datamatrix linnew;
  datamatrix linmode;
  datamatrix diff;
  datamatrix * linoldp;
  datamatrix * linnewp;

  void find_const(datamatrix & design);

  void create_matrices(void);

  // FUNCTION: compute_XWX
  // TASK: computes XWX on the basis of the current working weight and stores
  //       the result in r

  virtual void compute_XWX(datamatrix & r);
  virtual void compute_XWXroot(datamatrix & r);
  void compute_Wpartres(datamatrix & linpred);
  double compute_XtWpartres(double & mo);

  void add_linpred(datamatrix & l);

  public:

  vector<ST::string> datanames;              // names of covariates

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_linear(void);

  // CONSTRUCTOR

  FC_linear(MASTER_OBJ * mp, unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
            datamatrix & d, vector<ST::string> & vn, const ST::string & t,
           const ST::string & fp,bool cent);

  // COPY CONSTRUCTOR

  FC_linear(const FC_linear & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_linear & operator=(const FC_linear & m);

  // DESTRUCTOR

  ~FC_linear()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  void update_gaussian(void);
  void update_IWLS(void);

  // FUNCTION: posteriormode

  bool posteriormode(void);

  // FUNCTION: outoptions

  void outoptions(void);

  // FUNCTION: outresults

  void outresults(ofstream & out_stata,ofstream & out_R,
                  const ST::string & pathresults);

  void compute_autocorr_all(const ST::string & path,
                              unsigned lag, ofstream & outg) const;


  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  // FUNCTION: reset

  void reset(void);

  // FUNCTION: add_variable

  int add_variable(const datamatrix & d,ST::string & name);

  };


//------------------------------------------------------------------------------
//------------------------ CLASS: FC_linear_pen --------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE FC_linear_pen  : public FC_linear
  {

  protected:

  public:

  datamatrix  tau2;
  datamatrix  tau2oldinv;

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FC_linear_pen(void);

  // CONSTRUCTOR

  FC_linear_pen(MASTER_OBJ * mp,unsigned & enr,
                GENERAL_OPTIONS * o,DISTR * lp, datamatrix & d,
                vector<ST::string> & vn, const ST::string & t,
                const ST::string & fp,bool cent);

  // COPY CONSTRUCTOR

  FC_linear_pen(const FC_linear_pen & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FC_linear_pen & operator=(const FC_linear_pen & m);

  // DESTRUCTOR

  ~FC_linear_pen()
    {
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void update(void);

  // FUNCTION: posteriormode

  bool posteriormode(void);

  void compute_XWX(datamatrix & r);
  void compute_XWXroot(datamatrix & r);

  // FUNCTION: outoptions

  void outoptions(void);

  // FUNCTION: outresults

  void outresults(ofstream & out_stata,ofstream & out_R,
                  const ST::string & pathresults);

  // FUNCTION: reset

  void reset(void);

  };



} // end: namespace MCMC

#endif


