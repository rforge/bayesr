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



#if !defined (MODELREMLREG_INCLUDED)
#define MODELREMLREG_INCLUDED

#include"../export_type.h"
#include"model.h"

//------------------------------------------------------------------------------
//------------------------- class term_fixed_catspecific -----------------------
//------------------------------------------------------------------------------

// Category-specific fixed effects in cumulative or sequential models

class __EXPORT_TYPE term_fixed_catspecific : public basic_termtype
  {
  protected:

  void setdefault(void)
    {
    }

  public:

  // DEFAULT CONSTRUCTOR

  term_fixed_catspecific(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a category-specific fixed effect

  bool checkvector(const vector<term>  & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_fixed_catspecific() {}
  };

//------------------------------------------------------------------------------
//-------------------------- class term_autoreg_remlreg ------------------------
//------------------------------------------------------------------------------

// First or second order random walk

class __EXPORT_TYPE term_autoreg_remlreg : public basic_termtype
  {
  protected:

  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_autoreg_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a first or second order random walk

  bool checkvector(const vector<term>  & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_autoreg_remlreg() {}

  };

//------------------------------------------------------------------------------
//-------------------- class term_autoreg_varcoef_remlreg ----------------------
//------------------------------------------------------------------------------

// First or second order random walk as effect modifier in a VCM

class __EXPORT_TYPE term_autoreg_varcoef_remlreg : public basic_termtype
  {
  protected:

  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;
  simpleoption center;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_autoreg_varcoef_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a VCM with first or second order random
  //       walk as effect modifier

  bool checkvector(const vector<term>  & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_autoreg_varcoef_remlreg() {}

  };

//------------------------------------------------------------------------------
//-------------------------- class term_season_remlreg -------------------------
//------------------------------------------------------------------------------

// Seasonal Effect

class __EXPORT_TYPE term_season_remlreg : public basic_termtype
  {
  protected:

  intoption period;
  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_season_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a seasonal component

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_season_remlreg() {}

  };

//------------------------------------------------------------------------------
//------------------- class term_season_varcoef_remlreg ------------------------
//------------------------------------------------------------------------------

// Seasonal effect as effect modifier in a VCM

class __EXPORT_TYPE term_season_varcoef_remlreg : public basic_termtype
  {
  protected:

  intoption period;
  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_season_varcoef_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a VCM with a seasonal component as
  //       effect modifier

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_season_varcoef_remlreg() {}

  };

//------------------------------------------------------------------------------
//------------------------- class term_pspline_remlreg -------------------------
//------------------------------------------------------------------------------

// P-spline with first or second order random walk penalty

class __EXPORT_TYPE term_pspline_remlreg : public basic_termtype
  {
  protected:

  intoption degree;                // degree of the b-spline basis
  intoption numberknots;           // number of knots
  doubleoption lambda;
  intoption gridsize;
  simpleoption diagtransform;
  simpleoption derivative;
  doubleoption lambdastart;
  simpleoption catspecific;
  doubleoption lowergrid;
  doubleoption uppergrid;
  doubleoption lowerknot;
  doubleoption upperknot;
  doubleoption reference;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_pspline_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a psline with rw1 oder rw2 penalty

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_pspline_remlreg() {}

  };

//------------------------------------------------------------------------------
//------------------- class term_varcoeff_pspline_remlreg ----------------------
//------------------------------------------------------------------------------

// P-spline with first or second order random walk prior as effect modifier in
// a VCM

class __EXPORT_TYPE term_varcoeff_pspline_remlreg : public basic_termtype
  {
  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;
  simpleoption center;
  doubleoption reference;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_varcoeff_pspline_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a VCM with a P-spline as effect modifier

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_varcoeff_pspline_remlreg() {}

  };

//------------------------------------------------------------------------------
//--------------------- class term_baseline_remlreg ----------------------------
//------------------------------------------------------------------------------

// p-spline for the baseline in a cox model

class __EXPORT_TYPE term_baseline_remlreg : public basic_termtype
  {
  protected:

  intoption degree;           // degree of the b-spline basis
  intoption numberknots;      // number of knots for the b-spline
  intoption tgrid;            // number of equidistant gridpoints
  stroption gridchoice;       // equidistant or quantile grid for numerical
                              // integration
  intoption numberquantiles;  // numberof quantiles
  intoption numberbetween;    // number of gridpoints between quantiles
  doubleoption lambda;
  doubleoption lambdastart;
  stroption lower;            // lower boundary of the interval for interval
                              // censored data
  simpleoption catspecific;
  intoption gridsize;
  doubleoption reference;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_baseline_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a p-spline for the baseline in a cox model

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_baseline_remlreg() {}

  };

//------------------------------------------------------------------------------
//------------------ class term_baseline_varcoef_remlreg -----------------------
//------------------------------------------------------------------------------

// p-spline for time-varying effects in a cox model

class __EXPORT_TYPE term_baseline_varcoeff_remlreg : public basic_termtype
  {
  protected:

  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;
  intoption gridsize;
  doubleoption reference;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_baseline_varcoeff_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a p-spline for time-varying effects in
  //       a cox model

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_baseline_varcoeff_remlreg() {}

  };

//------------------------------------------------------------------------------
//-------------------- class term_interactpspline_remlreg ----------------------
//------------------------------------------------------------------------------

// P-spline interaction surface

class __EXPORT_TYPE term_interactpspline_remlreg : public basic_termtype
  {
  protected:

  intoption degree;               // degree of the b-spline basis
  intoption numberknots;          // number of knots per dimension
  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;
  intoption gridsizex;
  intoption gridsizey;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_interactpspline_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a p-spline interaction surface

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_interactpspline_remlreg() {}

  };

//------------------------------------------------------------------------------
//-------------- class term_interactpspline_varcoeff_remlreg -------------------
//------------------------------------------------------------------------------

// VCM with p-spline interaction surface as effect modifier

class __EXPORT_TYPE term_interactpspline_varcoeff_remlreg : public basic_termtype
  {
  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;
  simpleoption center;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_interactpspline_varcoeff_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a VCM with p-spline interaction surface
  //       as effect modifier

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_interactpspline_varcoeff_remlreg() {}

  };

//------------------------------------------------------------------------------
//--------------------------- class term_geospline_remlreg ---------------------
//------------------------------------------------------------------------------

// p-spline interaction surface basedon centroids of regions

class __EXPORT_TYPE term_geospline_remlreg : public basic_termtype
  {
  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  stroption map;               // name of the map-object containing the
                               // geographical information
  doubleoption lambdastart;
  simpleoption catspecific;
  intoption gridsizex;
  intoption gridsizey;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_geospline_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a p-spline interaction surface based on
  //       centroids of regions

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_geospline_remlreg() {}

  };

//------------------------------------------------------------------------------
//----------------- class term_geopspline_varcoeff_remlreg ---------------------
//------------------------------------------------------------------------------

// VCM with p-spline interaction surface based on centroids of regions as
// effect modifier

class __EXPORT_TYPE term_geospline_varcoeff_remlreg : public basic_termtype
  {
  protected:

  intoption degree;
  intoption numberknots;
  stroption map;
  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;
  simpleoption center;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_geospline_varcoeff_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a VCM with p-spline interaction surface
  //       as effect modifier

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_geospline_varcoeff_remlreg() {}

  };

//------------------------------------------------------------------------------
//-------------------------- class term_spatial_remlreg ------------------------
//------------------------------------------------------------------------------

// Spatial effect with Markov random field prior

class __EXPORT_TYPE term_spatial_remlreg : public basic_termtype
  {
  protected:

  stroption map;                    // Name of the map-object containing the
                                    // geographical information
  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_spatial_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a spatial component with MRF prior

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_spatial_remlreg() {}

  };

//------------------------------------------------------------------------------
//--------------------- class term_spatial_varcoef_remlreg ---------------------
//------------------------------------------------------------------------------

// Spatial effect with MRF prior as effect modifier in a VCM

class __EXPORT_TYPE term_spatial_varcoef_remlreg : public basic_termtype
  {
  protected:

  stroption map;
  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;
  simpleoption center;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_spatial_varcoef_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a VCM with spatial component as effect
  //       modifier

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_spatial_varcoef_remlreg() {}

  };

//------------------------------------------------------------------------------
//---------------------- class term_random_remlreg -----------------------------
//------------------------------------------------------------------------------

// iid random effect

class __EXPORT_TYPE term_random_remlreg : public basic_termtype
  {
  protected:

  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_random_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a random intercept

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_random_remlreg() {}

  };

//------------------------------------------------------------------------------
//----------------------- class term_randomslope_remlreg -----------------------
//------------------------------------------------------------------------------

// Random slope

class __EXPORT_TYPE term_randomslope_remlreg : public basic_termtype
  {
  protected:

  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_randomslope_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a random slope

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_randomslope_remlreg() {}

  };

//------------------------------------------------------------------------------
//--------------------- class term_kriging_1dim_remlreg ------------------------
//------------------------------------------------------------------------------

// one-dimensional stationary gaussian random field

class __EXPORT_TYPE term_kriging_1dim_remlreg : public basic_termtype
  {
  protected:

  doubleoption nu;              // Parameter of the matern correlation function
  doubleoption maxdist;         // Distance involved in the determination of rho
  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_kriging_1dim_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a one-dimensional stationary gaussian
  //       random field

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_kriging_1dim_remlreg() {}

  };

//------------------------------------------------------------------------------
//------------------------ class term_kriging_remlreg --------------------------
//------------------------------------------------------------------------------

// two-dimensional stationary gaussian random field

class __EXPORT_TYPE term_kriging_remlreg : public basic_termtype
  {
  protected:

  intoption numberknots;    // number of knots to be computed by the space
                            // filling algorithm
  doubleoption nu;          // parameter of the matern correlation function
  doubleoption maxdist;     // distance involved in the determination of rho
  simpleoption full;        // use all different points as knots
  stroption knotdata;       // dataset-object providing the knots
  doubleoption p;           // Options for the space-filling algorithm
  doubleoption q;
  intoption maxsteps;       // maximum number of steps for space filling
  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;
  intoption gridsizex;
  intoption gridsizey;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_kriging_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a two dimensional stationary gaussian
  //       random field

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_kriging_remlreg() {}

  };

//------------------------------------------------------------------------------
//------------------ class term_kriging_varcoeff_remlreg -----------------------
//------------------------------------------------------------------------------

// VCM with two-dimensional stationary gaussian random field as effect modifier

class __EXPORT_TYPE term_kriging_varcoeff_remlreg : public basic_termtype
  {
  protected:

  intoption numberknots;
  doubleoption nu;
  doubleoption maxdist;
  simpleoption full;
  stroption knotdata;
  doubleoption p;
  doubleoption q;
  intoption maxsteps;
  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption catspecific;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_kriging_varcoeff_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a VCM with two-dimensional stationary
  //       gaussian random field as effect modifier

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_kriging_varcoeff_remlreg() {}

  };

//------------------------------------------------------------------------------
//----------------------- class term_geokriging_remlreg ------------------------
//------------------------------------------------------------------------------

// Two-dimensional stationary gaussian random field based on the centroids of
// regions

class __EXPORT_TYPE term_geokriging_remlreg : public basic_termtype
  {
  protected:

  intoption numberknots;
  doubleoption nu;
  doubleoption maxdist;
  simpleoption full;
  stroption knotdata;
  doubleoption p;
  doubleoption q;
  intoption maxsteps;
  doubleoption lambda;
  doubleoption lambdastart;
  stroption map;            // name of the map-object containing the
                            // geographical information
  simpleoption catspecific;
  intoption gridsizex;
  intoption gridsizey;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_geokriging_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a two-dimensional stationary gaussian
  //       random field based on centroids of regions

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_geokriging_remlreg() {}

  };

//------------------------------------------------------------------------------
//---------------- class term_geokriging_varcoeff_remlreg ----------------------
//------------------------------------------------------------------------------

// VCM with two-dimensional stationary gaussian random field based on centroids
// as effect modifier

class __EXPORT_TYPE term_geokriging_varcoeff_remlreg : public basic_termtype
  {
  protected:

  intoption numberknots;
  doubleoption nu;
  doubleoption maxdist;
  simpleoption full;
  stroption knotdata;
  doubleoption p;
  doubleoption q;
  intoption maxsteps;
  doubleoption lambda;
  doubleoption lambdastart;
  stroption map;
  simpleoption catspecific;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_geokriging_varcoeff_remlreg(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a VCM with two-dimensional stationary
  //       gaussian random field based on centroids as effect modifier

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_geokriging_varcoeff_remlreg() {}

  };

#endif
