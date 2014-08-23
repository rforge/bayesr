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

#include "../export_type.h"

#if !defined (DESIGNkrigingINCLUDED)

#define DESIGNkrigingINCLUDED

#include"statmat.h"
#include"design.h"
#include"Random.h"
#include"envmatrix_penalty.h"
#include"../values.h"
#include<fstream>
#include<vector>
#include<bitset>
#include"clstring.h"
#include<cmath>


namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_kriging ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_kriging : public DESIGN
  {

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  double compute_matern(double & nu,double & r);

  void compute_tildeZ(void);  

  protected:

  double rho;
  double maxdist;
  double nu;

  ST::string knotdatapath;

  vector<double> xknots;              // x-und y-Koordinaten der Knoten
  vector<double> yknots;

  vector<double> xvalues;             // unterschiedliche Werte der Kovariablen
  vector<double> yvalues;

  void read_knots_from_data(void);

  void compute_knots(const vector<double> & xvals,
                     const vector<double> & yvals,
                     unsigned nrknots,double p,double q,
                     vector<double> & xknots,
                     vector<double> & yknots);

  long nrknots;

  // ---------------------- Variablen für geokriging ---------------------------

  MAP::map ma;
  bool mapexisting;
  datamatrix regions;


  // ---------------------- Variablen für geokriging ---------------------------

  datamatrix tildeZ_t;

  datamatrix XWXfull;
  datamatrix WsumtildeZ;

  datamatrix Kfull;

  public:

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_kriging(void);

  void help_construct(const datamatrix & dmy,const datamatrix & iv,
                 vector<ST::string> & op, vector<ST::string> & vn,
                 datamatrix & kd);

  // CONSTRUCTOR 1
  // x,y covariates

  DESIGN_kriging(const datamatrix & dm, const datamatrix & iv,
                 GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl,
                 vector<ST::string> & op,
                 vector<ST::string> & vn, datamatrix & kd);

  // CONSTRUCTOR 2
  // spatial covariate

  DESIGN_kriging(const datamatrix & dm,const datamatrix & iv,
                 const MAP::map & mp,
                 GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl,
                 vector<ST::string> & op, vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  DESIGN_kriging(const DESIGN_kriging & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_kriging & operator=(const DESIGN_kriging & m);

  // VIRTUAL FUNCTIONS

  void init_data(const datamatrix & dm,const datamatrix & iv);

  void compute_penalty(void);

  void compute_XtransposedWres(datamatrix & partres, double l);

  void compute_XtransposedWX(void);

  void compute_precision(double l);

  void outoptions(GENERAL_OPTIONS * op);

  void compute_orthogonaldecomp(void);

  double penalty_compute_quadform(datamatrix & beta);

  // DESTRUCTOR

  ~DESIGN_kriging() {}

  };


} // end: namespace MCMC

#endif


