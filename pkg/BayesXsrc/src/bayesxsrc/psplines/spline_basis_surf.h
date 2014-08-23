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



#ifndef spline_basis_surfH
#define spline_basis_surfH

#include"../export_type.h"
#include "fullcond.h"
#include "mcmc_nonpbasis.h"
#include "spline_basis.h"
#include<deque>
#include "statmat_penalty.h"

namespace MCMC
{

//---------------------------------------------------------------------------
//----------------------- class: spline_basis_surf --------------------------
//---------------------------------------------------------------------------


class __EXPORT_TYPE spline_basis_surf : public FULLCOND_nonp_basis
  {

  protected:

  FULLCOND_const * fcconst;

  bool lambdaconst;

  vector<int> index2;

  MAP::map m;
  bool mapexisting;
  ST::string mapname;
  vector<ST::string> regionnames;

  FULLCOND fchelp;
  FULLCOND fctotal;

  spline_basis * mainp1;
  spline_basis * mainp2;

  unsigned nrdiffobs;                // Anzahl verschiedender Beobachtungspaare
  unsigned nrpar1dim;                // Anzahl Parameter in einer Dimension (d.h. nrpar=nrpar1dim^2)

  ST::string fctotalrespath;
  ST::string outfile;

  bool centertotal;

  int gridsize;
  int gridsizex;
  int gridsizey;

  vector<double> xv;                 // geordnete, beobachtete x-Werte
  vector<double> yv;                 // geordnete, beobachtete y-Werte
  datamatrix xvalues;                // geordnete, äquidistante Werte zwischen Min(x/y) und Max(x/y)
  datamatrix yvalues;                // falls gridsize > 0
  vector<double> effectvaluesx;      // bildet zusammen mit effectvaluesy die Daten-Paare
  vector<double> effectvaluesy;      // für die Ausgabe

  datamatrix interactvar;

  datamatrix beta1;
  datamatrix beta2;

  datamatrix beta_uncentered;

  datamatrix he1;
  datamatrix he2;

  datamatrix betaweight;
  datamatrix betaweight_main;

  datamatrix betaweightx;
  datamatrix betaweighty;

  double intercept;                  // aktueller Intercept
  datamatrix spline;                 // geordnet wie im linearen Prädiktor
  datamatrix splinehelp;

  unsigned nrknots;                  // number of 'visible' knots
  unsigned degree;                   // degree of B-splines
  knotpos knpos;                     // knot choice (equidistant or quantiles)

  vector<int>     freq;              // Vektor der Länge N (=Anzahl Beobachtungen), dessen Werte sich
                                     // - beginnend bei 0 - immer, wenn die aktuelle Beobachtung der GEORDNETEN Daten
                                     // von der vorangehenden verschieden ist um 1 erhöhen und ansonsten gleich bleiben
                                     // (d.h. das letzte Element enthält die Anzahl VERSCHIEDENER Beobachtungen -1)
  vector<int>     freqoutput;
  deque<double>   knot1;             // Position of knots in x1 direction
  deque<double>   knot2;             // Position of knots in x2 direction

  datamatrix B;                      // Design-Matrix (enthält nur Werte ungleich 0)
  datamatrix Bout;                   // für Output bei VCM
  vector<int> first;                 // bezeichnet für jede Zeile der 'vollständigen' Design-Matrix
                                     // die Position der ersten Elements, das ungleich 0 ist

  datamatrix X_VCM;                  // für REML VCM
  datamatrix Z_VCM;                  // für REML VCM
  datamatrix X_grid;                  // für Ausgabe auf einem Gitter
  datamatrix Z_grid;                  // für Ausgabe auf einem Gitter
  vector<double> effectvaluesxgrid;      // bildet zusammen mit effectvaluesy die Daten-Paare
  vector<double> effectvaluesygrid;      // für die Ausgabe auf einem Gitter
  datamatrix xvaluesgrid;                // geordnete, äquidistante Werte zwischen Min(x/y) und Max(x/y)
  datamatrix yvaluesgrid;                // falls gridsize > 0

  datamatrix DG;                     // B-Spline-Matrix für regelmäßiges 50*50-Gitter
  vector<int> DGfirst;


  // FUNCTION: add_linearpred_multBS
  // TASK: multipliziert die Spalten 'a' bis 'e' von 'X' (Desing-Matrix) mit 'b'
  //       und addiert das Ergebnis zum linearen Prädiktor ('spline' wird nicht berechnet!!!)

  void add_linearpred_multBS_Block(const datamatrix & b,const unsigned a,const unsigned e,
                                   const unsigned beg,const unsigned end);

  // FUNCTION: multBS
  // TASK: multipliziert 'X' (Desing-Matrix) mit 'b' und weist das Ergebnis 'res'
  // NOTE: b und res müssen die richtigen Dimensionen haben (N x 1 bzw. nrpar x 1)

  void multBS(datamatrix & res, const datamatrix & b);
  void multBout(datamatrix & res, const datamatrix & b);

  void multBS_index(datamatrix & res, const datamatrix & b);

  void multBS_index_Block(datamatrix & res,const datamatrix & b,const unsigned a,const unsigned e,
                                           const unsigned beg,const unsigned end);

  // FUNCTION: multDG
  // TASK: multipliziert 'DG' (Desing-Matrix) mit 'b' und weist das Ergebnis 'res'
  // NOTE: b und res müssen die richtigen Dimensionen haben (N x 1 bzw. nrpar x 1)

  void multDG(datamatrix & res, const datamatrix & b);

  // FUNCTION: make_index
  // TASK: initializes 'freq', 'index'

  void make_index(const datamatrix & var1,const datamatrix & var2);

  // FUNCTION: compute_knots
  // TASK: initializes 'knot1' and 'knot2'

  void compute_knots(const datamatrix & x,const datamatrix & y);

  void make_xy_v(datamatrix var1,datamatrix var2);

  void make_xy_values(const datamatrix & var1,const datamatrix & var2);

  void make_xy_values_REML(const datamatrix & var1,const datamatrix & var2);

  // FUNCTION: make_B
  // TASK: creates matrix of B-Splines basis functions for (x,y)-values

  void make_B(const datamatrix & x,const datamatrix & y);

  void make_BVC(const datamatrix & intact);

  void make_DG(void);

  void make_DG_REML(void);

  datamatrix bspline(const double x, const double y);

  datamatrix bspline(const double x, const deque<double> knot);

  // FUNCTION: bspline_rek
  // TASK: needed in function 'predict' to compute B-Splines for a single observation

  double bspline_rek(unsigned l, unsigned knot, const double value, const bool xvar);

  bool breakpause(void);

  void make_index2(void);

  void compute_beta(void);

  void compute_main(void);


  public:


  // DEFAULT CONSTRUCTOR

  spline_basis_surf(void) : FULLCOND_nonp_basis()
    {
    }

  // CONSTRUCTOR

  spline_basis_surf(MCMCoptions * o, DISTRIBUTION * dp,
                 FULLCOND_const * fcc, const fieldtype & ft,
                const ST::string & ti, const unsigned & nrk, const unsigned & degr,
                const MCMC::knotpos & kp, const double & l, const int & gsx, const int & gsy,
                const ST::string & fp, const ST::string & pres, const unsigned & c);

  // CONSTRUCTOR for REML

 spline_basis_surf(MCMCoptions * o, const datamatrix & v1, const datamatrix & v2,
               const unsigned & nrk, const unsigned & degr,
               const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp, const unsigned & grx,
               const unsigned & gry);

  // CONSTRUCTOR for REML (varcoeff)

 spline_basis_surf(MCMCoptions * o, const datamatrix & intact,
               const datamatrix & v1, const datamatrix & v2,
               const unsigned & nrk, const unsigned & degr,
               const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp, const bool & ctr);

  // CONSTRUCTOR for REML (geospline)

  spline_basis_surf(MCMCoptions * o, const datamatrix & region, const MAP::map & mp,
               const ST::string & mn,
               const unsigned & nrk, const unsigned & degr,
               const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp, const unsigned & grx,
               const unsigned & gry);

  // CONSTRUCTOR for REML (geospline_varcoeff)

  spline_basis_surf(MCMCoptions * o, const datamatrix & intact,
               const datamatrix & region, const MAP::map & mp,
               const ST::string & mn,
               const unsigned & nrk, const unsigned & degr,
               const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp, const bool & ctr);

  // COPY CONSTRUCTOR

  spline_basis_surf(const spline_basis_surf & sp);

  // OVERLOADED ASSIGNMENT OPERATOR

  const spline_basis_surf & operator=(const spline_basis_surf & sp);

  void outresults(void);

  // FUNCTION: predict
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred);

  void compute_intercept(void);

  void compute_betaweight(void);

  void compute_betaweightxy(void);

  // falls gridsize > 0
  void compute_betaweightx2(void);
  void compute_betaweighty2(void);

  ST::string getinfo(void)
    {
    if(mapexisting)
      return mapname;
    else
      return title;
    }

  void init_names(const vector<ST::string> & na);

  void set_lambdaconst(double la);


  // REML

  void createreml(datamatrix & X,datamatrix & Z,const unsigned & Xpos,
                  const unsigned & Zpos);

  double outresultsreml(datamatrix & X,datamatrix & Z,
                                     datamatrix & betareml,
                                     datamatrix & betacov,
                                     datamatrix & thetareml,
                                     const unsigned & Xpos,
                                     const unsigned & Zpos,
                                     const unsigned & thetapos,
                                     const bool & dispers,
                                     const unsigned & betaXpos,
                                     const unsigned & betaZpos,
                                     const double & category,
                                     const bool & ismultinomial,
                                     const unsigned plotpos);

  void outresultsgrid();

  void outoptionsreml();


  // DESTRUCTOR

  ~spline_basis_surf(){}

  };

} // end: namespace MCMC


//---------------------------------------------------------------------------
#endif
