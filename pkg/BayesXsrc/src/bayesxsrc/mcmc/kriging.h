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



#if !defined (kriging_INCLUDED)
#define kriging_INCLUDED

#include"../export_type.h"
#include"fullcond.h"
#include"mcmc_nonpbasis.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//------------------------------- class: kriging -------------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE FULLCOND_kriging : public FULLCOND_nonp_basis
  {
  protected:

  unsigned nrknots;
  double nu;
  double rho;
  double maxdist;
  bool full;
  bool spacefill;

  bool onedim;
  vector<unsigned> incidence;         // contains for each observation the position in the knots-vector

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

  unsigned gridsize;                 // Variablen zur Auswertung des Effekts auf einem vorgegebenen Gitter
  unsigned gridsizex;
  unsigned gridsizey;

  datamatrix X_grid;
  datamatrix Z_grid;

  vector<double> effectvaluesxgrid;      // bildet zusammen mit effectvaluesy die Daten-Paare
  vector<double> effectvaluesygrid;      // für die Ausgabe auf einem Gitter
  datamatrix xvaluesgrid;                // geordnete, äquidistante Werte zwischen Min(x/y) und Max(x/y)
  datamatrix yvaluesgrid;                // falls gridsize > 0

  vector<int> index2;

  vector<int>freq;
  vector<int>freqoutput;
  unsigned nrdiffobs;

  datamatrix X_VCM;                  // für REML VCM
  datamatrix Z_VCM;                  // für REML VCM

  void make_index(const datamatrix & var1,const datamatrix & var2);
  void make_xy_values(const datamatrix & var1,const datamatrix & var2);
  void compute_knots(const vector<double> & xvals,const vector<double> & yvals);

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_kriging(void) : FULLCOND_nonp_basis()
    {
    }

  // COPY CONSTRUCTOR

  FULLCOND_kriging(const FULLCOND_kriging & kr);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_kriging & operator=(const FULLCOND_kriging & kr);

  // Constructor 1

  FULLCOND_kriging(MCMCoptions * o, const datamatrix & v1,
               const datamatrix & v2, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp, const unsigned &gsx,
               const unsigned & gsy);


  // Constructor 2: varcoef

  FULLCOND_kriging(MCMCoptions * o, const datamatrix & intact,
               const datamatrix & v1,
               const datamatrix & v2, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp);

  // Constructor 3: geokriging

  FULLCOND_kriging(MCMCoptions * o, const datamatrix & region,
               const MAP::map & mp, const ST::string & mn, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp, const unsigned & gsx,
               const unsigned & gsy);

  // Constructor 4: geokriging (varcoeff)

  FULLCOND_kriging(MCMCoptions * o, const datamatrix & intact,
               const datamatrix & region, const MAP::map & mp,
               const ST::string & mn, const datamatrix & knotdata,
               const unsigned & nrk, const double & n, const double & maxd,
               const double & pval, const double & qval, const unsigned & maxst,
               const bool & fu, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp);

  // Constructor 4

  FULLCOND_kriging(MCMCoptions * o, const datamatrix & v,
               const double & n, const double & maxd,
               const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp);

  // DESTRUCTOR

  ~FULLCOND_kriging(){}

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

  void outoptionsreml();

  void init_names(const vector<ST::string> & na);

  ST::string getinfo(void);

  void make_index(const datamatrix & moddata);

  void make_xy_values_grid(const datamatrix & var1,const datamatrix & var2);

  };


} // end: namespace MCMC

//---------------------------------------------------------------------------
#endif
