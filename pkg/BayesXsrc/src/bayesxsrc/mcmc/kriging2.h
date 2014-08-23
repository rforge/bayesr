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



#if !defined (kriging2_INCLUDED)
#define kriging2_INCLUDED

#include"../export_type.h"
#include"fullcond.h"
#include"mcmc_nonpbasis.h"
#include"spline_basis.h"
#include"fullcond_nonp_gaussian.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//------------------------------- class: kriging2 ------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FULLCOND_kriging2 : public spline_basis
  {

  protected:

  datamatrix X;
  updatetype utype;
  unsigned updateW;

  double nu;
  double rho;
  double maxdist;
  bool full;
  bool spacefill;

  MAP::map m;                         // Variablen für geokriging
  bool mapexisting;
  ST::string mapname;
  vector<ST::string> regionnames;

  double p;                           // p und q für Space-Filling-Algorithmus
  double q;
  unsigned maxsteps;

  vector<double> xknots;              // x-und y-Koordinaten der Knoten
  vector<double> yknots;

  vector<double> xvalues;             // unterschiedliche Werte der Kovariablen
  vector<double> yvalues;

  datamatrix xorig;                   // Original x- und y-Variable.
  datamatrix yorig;

  void make_index(const datamatrix & var1,const datamatrix & var2);
  void make_xy_values(const datamatrix & var1,const datamatrix & var2);
  void compute_knots(const vector<double> & xvals,const vector<double> & yvals);

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_kriging2(void) : spline_basis()
    {
    }

  // COPY CONSTRUCTOR

  FULLCOND_kriging2(const FULLCOND_kriging2 & kr);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_kriging2 & operator=(const FULLCOND_kriging2 & kr);

  // Constructor1

  FULLCOND_kriging2(MCMCoptions * o, DISTRIBUTION * dp, FULLCOND_const * fcc,
               const datamatrix & v1, const datamatrix & v2, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const unsigned & c=0);

  // Constructor2: geokriging

  FULLCOND_kriging2(MCMCoptions * o, DISTRIBUTION * dp, FULLCOND_const * fcc, const datamatrix & region,
               const MAP::map & mp, const ST::string & mn, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const unsigned & c=0);

  // Constructor3: geokriging IWLS

  FULLCOND_kriging2(MCMCoptions * o, DISTRIBUTION * dp, FULLCOND_const * fcc, const datamatrix & region,
               const MAP::map & mp, const ST::string & mn, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & mode, const unsigned & upW,
               const bool & updatetau, const double & fstart, const unsigned & c=0);

  // DESTRUCTOR

  ~FULLCOND_kriging2(){}

  void create(void);

  void update(void);

  void update_gaussian(void);

  void update_iwls(void);

  void update_iwlsmode(void);

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  double compute_quadform(void);

  void outresults(void);

  void outoptions();

  void init_names(const vector<ST::string> & na);

  ST::string getinfo(void);

  };


} // end: namespace MCMC

//---------------------------------------------------------------------------
#endif
