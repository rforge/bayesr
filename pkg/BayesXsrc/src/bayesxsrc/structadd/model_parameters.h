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



#if !defined (MODELparameters_INCLUDED)
#define MODELparameters_INCLUDED

#include"../export_type.h"
#include"clstring.h"
#include<vector>
#include"data.h"
#include"option.h"
#include"model.h"



//------------------------------------------------------------------------------
//--------------------------- class term_nonp ----------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE term_nonp : public basic_termtype
  {
  protected:

  intoption degree;
  intoption numberknots;
  intoption difforder;
  doubleoption lambda;
  doubleoption a;
  doubleoption b;
  simpleoption nocenter;
  stroption map;
  doubleoption lambda_re;
  doubleoption a_re;
  doubleoption b_re;
  simpleoption internal_mult;
  simpleoption samplemult;
  stroption constraints;
  doubleoption round;
  stroption centermethod;
  simpleoption internal_multexp;
  simpleoption pvalue;
  simpleoption meaneffect;
  doubleoption binning;
  stroption update;
  stroption nu;
  doubleoption maxdist;
  simpleoption ccovariate;
  doubleoption sum2;
  simpleoption derivative;
  simpleoption samplederivative;
  simpleoption samplef;
  doubleoption shrinkage;
  simpleoption shrinkagefix;
  doubleoption shrinkageweight;
  simpleoption adaptiveshrinkage;
  doubleoption tau2;
  doubleoption meaneffectconst;
  stroption prior;
  fileoption knotpath;
  stroption datasetref;
  simpleoption lambdaconst;

  doubleoption abeta;
  doubleoption bbeta;
  doubleoption r;
  doubleoption v;
  doubleoption aQ;
  doubleoption bQ;
  intoption regiterates;
  simpleoption center;

  doubleoption tildea;
  doubleoption tildeb;
  simpleoption cauchy;

  vector<ST::string> termnames;

  void setdefault(void);

  public:

  // DEFAULT CONSTRUCTOR

  term_nonp(void) {}

  // CONSTRUCTOR

  term_nonp(vector<ST::string> & na);

  // FUNCTION: check

  bool check(term & t);

  // FUNCTION: checkvector
  // TASK: returns true if term 'i' is a nonparametric term

  bool checkvector(const vector<term> & terms,const unsigned & i);

  // DESTRUCTOR

  ~term_nonp() {}

  };

#endif
