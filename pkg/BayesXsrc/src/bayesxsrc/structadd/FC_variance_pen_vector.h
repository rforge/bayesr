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



#if !defined (VARIANCEPEN_VECTOR_INCLUDED)
#define VARIANCEPEN_VECTOR_INCLUDED

#include"../export_type.h"
#include"FC_linear.h"

namespace MCMC
{


class __EXPORT_TYPE FC_variance_pen_vector : public FC
  {

  protected:

  vector<double> tau2;                 //  tau^2

  bool update_sigma2;                 //  update sigma2

  FC_linear_pen * Cp;

  DISTR * distrp;
  datamatrix shelp;

  FC FC_shrinkage;
  vector<bool> shrinkagefix;          //  Shrinkageparameter fix, at value shrinkagestart
  vector<bool> adaptiveshrinkage;     //  each variance get own Shrinkageparameter
  vector<double> a_shrinkagegamma;    //  Hyperparameter for Shrinkageparameter
  vector<double> b_shrinkagegamma;    //  Hyperparameter for Shrinkageparameter
  vector<double> shrinkagestart;      //  Startvalues for shrinkageparameters
  vector<double> shrinkageweight;     //  Weights for shrinkage


  double pensum;                    //  sum(beta^2/tau^2)
  int nrpen;

  bool is_ridge;          //  indicates if "true" the L2-penalty
                          //  and if "false" the L1-penalty is used
  bool is_fix;            //  indicates if "true" that the Shrinkageparameter is fixed
  bool is_adaptive;       //  indicates if "true" that the Shrinkage is adaptive


//  void outresults_shrinkage(void);  //  Function to write results to output window and files

  void outresults_shrinkage(const ST::string & pathresults);

  public:

  //____________________________________________________________________________
  //
  // DEFAULT CONSTRUCTOR
  //____________________________________________________________________________

  FC_variance_pen_vector(void) : FC()
    {
    }


  //____________________________________________________________________________
  //
  // CONSTRUCTOR
  //____________________________________________________________________________

  FC_variance_pen_vector(MASTER_OBJ * mp,GENERAL_OPTIONS * o, FC_linear_pen * p,
                         DISTR * d,const ST::string & ti,
                         const ST::string & fp, bool isr);

  //____________________________________________________________________________
  //
  // COPY CONSTRUCTOR
  //____________________________________________________________________________

  FC_variance_pen_vector(const FC_variance_pen_vector & t);


  //____________________________________________________________________________
  //
  // OVERLOADED ASSIGNMENT OPERATOR
  //____________________________________________________________________________

  const FC_variance_pen_vector & operator=(const FC_variance_pen_vector & t);


  void add_variable(datamatrix & x,vector<ST::string> & op,
                         vector<ST::string> & vn);

  //____________________________________________________________________________
  //
  // OVERLOADED ASSIGNMENT OPERATOR
  //____________________________________________________________________________


  // Pointer auf das shrinkage-Parameter Fullcond-Objekt
  FC * get_shrinkagepointer();

  void get_samples(const ST::string & filename,ofstream & outg) const;
  //  void get_samples(const ST::string & filename, const unsigned & step = 1) const;

  //____________________________________________________________________________
  //
  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //____________________________________________________________________________

  void update(void);

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool posteriormode(void);


  //____________________________________________________________________________
  //
  // FUNCTION: outresults
  // TASK: - write results to output window and files
  //____________________________________________________________________________

  void outresults(ofstream & out_stata, ofstream & out_R,
                  const ST::string & pathresults);

  //____________________________________________________________________________
  //
  // FUNCTION: outoptions
  // TASK: - write options to output window
  //____________________________________________________________________________

  void outoptions(void);


  //____________________________________________________________________________
  //
  // FUNCTION: reset
  // TASK: resets all parameters
  //____________________________________________________________________________

  void reset(void)
    {
    /*
    FC::reset();
    setbeta(beta.rows(),1,0.1);
    */
    }


  //____________________________________________________________________________
  //
  // DESTRUCTOR
  //____________________________________________________________________________

  ~FC_variance_pen_vector() {}

  }; // end: class FC_variance_pen_vector


} // end: namespace MCMC

#endif

