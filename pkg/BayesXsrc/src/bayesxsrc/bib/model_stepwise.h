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



#if !defined (MODELSTEPWISE_INCLUDED)
#define MODELSTEPWISE_INCLUDED

#include"../export_type.h"
#include"model.h"

class __EXPORT_TYPE term_nonlinearf_stepwise : public basic_termtype
  {

  protected:

  doubleoption lambda;
  doubleoption lambdastart;
  simpleoption forced_into;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_nonlinearf_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a nonlinearf component

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( terms[i].type == "nonlinearf" )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_nonlinearf_stepwise() {}

  };


class __EXPORT_TYPE term_autoreg_stepwise : public basic_termtype
  {

  protected:

  doubleoption lambda;
  doubleoption spmin;
  doubleoption spmax;
  doubleoption spstart;
  simpleoption forced_into;
  doubleoption dfmax;
  doubleoption dfmin;
  doubleoption dfstart;
  //impleoption sp;
  stroption sp;
  intoption number;
  simpleoption logscale;
  doubleoption df_accuracy;
  simpleoption center;
  simpleoption nofixed;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_autoreg_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a first or second order random walk

  bool checkvector(const vector<term>  & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( (terms[i].type == "rw1") || (terms[i].type == "rw2") ||
         (terms[i].type == "varcoeffrw1") || (terms[i].type == "varcoeffrw2")
       )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_autoreg_stepwise() {}

  };


class __EXPORT_TYPE term_season_stepwise : public basic_termtype
  {

  protected:

  intoption period;
  doubleoption lambda;
  doubleoption spmin;
  doubleoption spmax;
  doubleoption spstart;
  simpleoption forced_into;
  doubleoption dfmin;
  doubleoption dfmax;
  doubleoption dfstart;
  //simpleoption sp;
  stroption sp;
  intoption number;
  simpleoption logscale;
  doubleoption df_accuracy;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_season_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a seasonal component

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( ( terms[i].type == "season" ) ||
         ( terms[i].type == "varcoeffseason")
       )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_season_stepwise() {}

  };


class __EXPORT_TYPE term_pspline_stepwise : public basic_termtype
  {

  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  intoption gridsize;
  simpleoption diagtransform;
  simpleoption derivative;
  doubleoption spmin;
  doubleoption spmax;
  doubleoption spstart;
  simpleoption forced_into;
  doubleoption dfmin;
  doubleoption dfmax;
  doubleoption dfstart;
  //simpleoption sp;
  stroption sp;
  intoption number;
  simpleoption logscale;
  doubleoption df_accuracy;
  stroption monotone;
  simpleoption center;
  stroption knots;
  simpleoption nofixed;
  doubleoption spmonotone;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_pspline_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a psline with rw1 oder rw2 Penalty

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( (terms[i].type == "psplinerw1") || (terms[i].type == "psplinerw2") ||
         (terms[i].type == "varpsplinerw1") || (terms[i].type == "varpsplinerw2") ||
         (terms[i].type == "psplinerw1rw2") || (terms[i].type == "varpsplinerw1rw2"))
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_pspline_stepwise() {}

  };


class __EXPORT_TYPE term_spatial_stepwise : public basic_termtype
  {

  protected:

  stroption map;
  doubleoption lambda;
  doubleoption spmin;
  doubleoption spmax;
  doubleoption spstart;
  simpleoption forced_into;
  doubleoption dfmin;
  doubleoption dfmax;
  doubleoption dfstart;
  //simpleoption sp;
  stroption sp;
  intoption number;
  simpleoption logscale;
  doubleoption df_accuracy;
  simpleoption center;
  simpleoption nofixed;
  stroption map2;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_spatial_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a spatial component

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( (terms[i].type == "spatial") ||  (terms[i].type == "varcoeffspatial")
         || (terms[i].type == "spatialrandom") || (terms[i].type == "twospatialrandom") )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_spatial_stepwise() {}

  };


class __EXPORT_TYPE term_randomslope_stepwise : public basic_termtype
  {

  protected:

  simpleoption center;
  doubleoption lambda;
  doubleoption spmin;
  doubleoption spmax;
  doubleoption spstart;
  simpleoption forced_into;
  doubleoption dfmin;
  doubleoption dfmax;
  doubleoption dfstart;
  //impleoption sp;
  stroption sp;
  intoption number;
  simpleoption logscale;
  doubleoption df_accuracy;
  simpleoption nofixed;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_randomslope_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a random effect

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( terms[i].type == "randomslope" )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_randomslope_stepwise() {}

  };


class __EXPORT_TYPE term_random_stepwise : public basic_termtype
  {

  protected:

  doubleoption lambda;
  doubleoption spmin;
  doubleoption spmax;
  doubleoption spstart;
  simpleoption forced_into;
  doubleoption dfmin;
  doubleoption dfmax;
  doubleoption dfstart;
  //simpleoption sp;
  stroption sp;
  intoption number;
  simpleoption logscale;
  doubleoption df_accuracy;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_random_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a random effect

  bool checkvector(const vector<term> & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if ( terms[i].type == "random" )
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_random_stepwise() {}

  };


class __EXPORT_TYPE term_factor_stepwise : public basic_termtype
  {

  protected:

  stroption coding;
  doubleoption reference;
  intoption spstart;
  simpleoption forced_into;
  intoption dfstart;
  simpleoption nofixed;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_factor_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a first or second order random walk

  bool checkvector(const vector<term>  & terms,const unsigned & i)
    {

    assert(i< terms.size());

    if (terms[i].type == "factor")
      return true;

    return false;
    }

  // DESTRUCTOR

  ~term_factor_stepwise() {}

  };


//------------------------------------------------------------------------------
//---------------------- class term_interactpspline  ---------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_interactpspline_stepwise : public basic_termtype
  {

  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  intoption gridsize;
  doubleoption spmin;
  doubleoption spmax;
  doubleoption spstart;
  simpleoption forced_into;
  doubleoption dfmin;
  doubleoption dfmax;
  doubleoption dfstart;
  //simpleoption sp;
  stroption sp;
  intoption number;
  simpleoption logscale;
  doubleoption df_accuracy;
  simpleoption center;
  simpleoption nofixed;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_interactpspline_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_interactpspline_stepwise() {}

  };


//------------------------------------------------------------------------------
//--------------------------- class term_geospline  ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE term_geospline_stepwise : public basic_termtype
  {

  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  stroption map;
  doubleoption spmin;
  doubleoption spmax;
  doubleoption spstart;
  simpleoption forced_into;
  doubleoption dfmin;
  doubleoption dfmax;
  doubleoption dfstart;
  //simpleoption sp;
  stroption sp;
  intoption number;
  simpleoption logscale;
  doubleoption df_accuracy;
  simpleoption nofixed;
  simpleoption center;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_geospline_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is an interaction term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_geospline_stepwise() {}

  };


class __EXPORT_TYPE term_projection_stepwise : public basic_termtype
  {

  protected:

  intoption degree;
  intoption numberknots;
  doubleoption lambda;
  intoption gridsize;
  simpleoption diagtransform;
  simpleoption derivative;
  doubleoption spmin;
  doubleoption spmax;
  doubleoption spstart;
  simpleoption forced_into;
  doubleoption dfmin;
  doubleoption dfmax;
  doubleoption dfstart;
  //simpleoption sp;
  stroption sp;
  intoption number;
  simpleoption logscale;
  doubleoption df_accuracy;
  stroption monotone;
  simpleoption center;
  intoption nterms;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_projection_stepwise(void);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a psline with rw1 oder rw2 Penalty

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_projection_stepwise() {}

  };


#endif
