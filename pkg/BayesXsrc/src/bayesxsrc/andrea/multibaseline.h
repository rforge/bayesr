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



#ifndef multibaselineH
#define multibaselineH

#include"../export_type.h"
#include"multistate.h"
#include"mcmc_pspline.h"
#include"spline_basis.h"
#include<vector>

namespace MCMC
{


//------------------------------------------------------------------------------
//---------------------------- class: pspline_multibaseline --------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE pspline_multibaseline : public FULLCOND_pspline
  {

  protected:

   bool begin0;
   unsigned col;
//   datamatrix int_knots;
//   datamatrix int_D;
   vector<datamatrix> int_D_l;
   datamatrix int_D;
   vector<datamatrix> int_knots_l;
   datamatrix int_knots;
   MCMC::bsplinemat testmat;
   vector<MCMC::bsplinemat> testmat_l;
   vector<pspline_multibaseline*> baselinep;
   datamatrix zi;
   datamatrix z_vc;
   datamatrix zi_ges;
   datamatrix beg_i;
   datamatrix state_i;
   statmatrix<int> zi_index;
   statmatrix<int> ges_index;
   vector< statmatrix<int> > teil_index;
   vector<datamatrix> zi_teil;
   vector<datamatrix> spline_teil;
   vector<datamatrix> spline_teil2;
   datamatrix spline_ges;
   datamatrix spline_ges2;
   datamatrix spline_zi;
   datamatrix int_ti_help;
   bool vc_dummy1;
   bool global;


  public:


  // DEFAULT CONSTRUCTOR

  pspline_multibaseline(void) : FULLCOND_pspline()
    {
    }

  // CONSTRUCTOR 1

  pspline_multibaseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d, const double & a, const double & b,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres,
                    const int & gs,const unsigned & c,const datamatrix & zustand,
                    const datamatrix & anfang,const bool & wb);


// CONSTRUCTOR 2 (für zeitabhängige Effekte)

pspline_multibaseline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & time, const datamatrix & z,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const unsigned & minb,const unsigned & maxb,
                    const fieldtype & ft,const ST::string & ti,
                    const ST::string & fp, const ST::string & pres,
                    const int & gs,const unsigned & c,const datamatrix & zustand,
                    const datamatrix & anfang,const bool & wb);


  // COPY CONSTRUCTOR

  pspline_multibaseline(const pspline_multibaseline & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const pspline_multibaseline & operator=(const pspline_multibaseline & fc);


  void update(void);


  void outoptions(void);
//  void outresults(void);
  void compute_int_ti_nonglobal(const datamatrix & b);
  void compute_int_ti_global(const datamatrix & b);
//void compute_int_ti_linear(const double & b);
//void compute_int_ti_weibull(const double & r);
  void compute_int_ti(unsigned beg);
  void compute_int_ti_vc_di0(const vector<double *>,const vector<double *>,const vector<double *>);
  void compute_int_ti_vc_di(const int,const vector<double *>,const vector<double *>,const vector<double *>);
//void compute_int_gauss(void);
//void compute_int_gauss_DIC(void);
  void update_multibaseline(void);
  void compute_int_ti_mean(void);
  void set_baselinep(vector<pspline_multibaseline*> bp)
    {
    baselinep=bp;
    }

  double * get_int_D(void)
    {
    return int_D.getV();
    }

  double * get_z_vc(void)
    {
    return z_vc.getV();
    }

  datamatrix get_z_vc_np(void)
    {
    return z_vc;
    }

  double * get_spline_zi(void)
    {
    multBS(spline_zi,beta);
    return spline_zi.getV();
    }

  double * get_spline_ges(void)
    {
    testmat.mult(spline_ges,beta);
    return spline_ges.getV();
    }

  double * get_spline_ges2(void)
    {
    testmat.mult_index(spline_ges2,beta);
    return spline_ges2.getV();
    }

  double * get_spline_ges_mean(void)
    {
    testmat.mult(spline_ges,betamean);
    return spline_ges.getV();
    }

  double * get_spline_ges2_mean(void)
    {
    testmat.mult_index(spline_ges2,betamean);
    return spline_ges2.getV();
    }

//double * get_gaussspline(void);


//double * get_gaussspline_mean(void);


  double * get_betamean(void)
    {
    return betamean.getV();
    }

  double * get_beta(void)
    {
    return beta.getV();
    }

  double * get_spline_zi_mean(void)
    {
    multBS(spline_zi,betamean);
    return spline_zi.getV();
    }

  void set_fcconst(FULLCOND_const * fcc)
    {
    fcconst = fcc;
    }


/*datamatrix lgamma;
  void create_lgamma(void);
  double lgammafunc(const double & nu) const;
  double lfac(const double & nu) const;*/



  // DESTRUCTOR

  ~pspline_multibaseline() {}

  };



}   // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
