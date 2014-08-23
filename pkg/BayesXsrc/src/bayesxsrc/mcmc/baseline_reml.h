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



#ifndef baseline_remlH
#define baseline_remlH

#include"../export_type.h"
#include"mcmc_pspline.h"
#include"spline_basis.h"
#include<vector>

namespace MCMC
{

//------------------------------------------------------------------------------
//---------------------------- class: baseline_reml ----------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE baseline_reml : public spline_basis
  {

  protected:

  double tstep;
  unsigned tgrid;
  datamatrix tsteps;

  knotpos gridpos;
  unsigned nrquant;
  unsigned nrbetween;

  // interval censoring and left truncation
  vector<unsigned>tright;
  vector<unsigned>tleft;
  vector<unsigned>ttrunc;

  datamatrix t_X;
  datamatrix t_Z;

  datamatrix interact_var;

  datamatrix tvalues;

  public:

  // DEFAULT CONSTRUCTOR

  baseline_reml(void) : spline_basis()
    {
    }

  // CONSTRUCTOR 1

  baseline_reml(MCMCoptions * o, const datamatrix & d,
               const datamatrix & leftint, const datamatrix & lefttrunc,
               const unsigned & nrk, const unsigned & degr, const unsigned & tgr,
               const unsigned & nrq, const unsigned & nrb, const knotpos & kp,
               const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const knotpos & gp, const int & gs,
               const bool & catsp, const double & rv);

  // CONSTRUCTOR 2 (VCM)

  baseline_reml(MCMCoptions * o,const datamatrix & d1,
                      const datamatrix & d2, const unsigned & nrk,
                      const unsigned & degr, const unsigned & tgr,
                      const knotpos & kp, const fieldtype & ft,
                      const ST::string & ti, const ST::string & fp,
                      const ST::string & pres, const double & l,
                      const double & sl, const int & gs, const bool & catsp,
                      const double & rv);

  // COPY CONSTRUCTOR

  baseline_reml(const baseline_reml & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const baseline_reml & operator=(const baseline_reml & fc);

  // DESTRUCTOR

  ~baseline_reml() {}

  void createreml(datamatrix & X,datamatrix & Z,const unsigned & Xpos,
                  const unsigned & Zpos);

  void multDG(datamatrix & res, const datamatrix & b);

  void initialize_baseline(unsigned j, datamatrix & tx, datamatrix & tz,
               vector<unsigned> & ts, vector<unsigned> & te,
               vector<unsigned> & tt, datamatrix & iv,
               statmatrix<double> & steps, statmatrix<int> & ind);

  void outoptionsreml();

  void init_name(const ST::string & na);

  unsigned & get_tgrid(void)
    {
    return tgrid;
    }

  };

}   // end: namespace MCMC

//---------------------------------------------------------------------------
#endif
