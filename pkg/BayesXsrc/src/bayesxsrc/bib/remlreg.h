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



#if !defined (remlreg_INCLUDED)

#define remlreg_INCLUDED

#include"../export_type.h"
#include"statobj.h"
#include"dataobj.h"
#include"map.h"
#include"mapobject.h"

#include"remlest.h"
#include"remlest_multi.h"
#include"remlest_multi2.h"
#include"remlest_multi3.h"
//#include<remlest_multistate.h>
#include"model_remlreg.h"

#include"mcmc.h"
#include"mcmc_const.h"
#include"fullcond_nonp_gaussian.h"
#include"spline_basis.h"
#include"spline_basis_surf.h"
#include"randomeffect.h"
#include"kriging.h"
#include"baseline_reml.h"

using MCMC::MCMCoptions;

using MCMC::FULLCOND;
using MCMC::FULLCOND_const;
using MCMC::FULLCOND_nonp_gaussian;
using MCMC::spline_basis;
using MCMC::spline_basis_surf;
using MCMC::FULLCOND_random;
using MCMC::FULLCOND_kriging;
using MCMC::baseline_reml;

class __EXPORT_TYPE remlreg : public statobject
  {
  private :

  void make_paths(unsigned collinpred,ST::string & pathnonp,
                          ST::string & pathres,ST::string & title,
                          ST::string varname1,ST::string varname2,
                          ST::string endingraw,
                          ST::string endingres,ST::string endingtitle) ;

  void initpointers(void);
  void create(void);
  void clear(void);

//------------------------------------------------------------------------------
//--------------------------- PRIVATE VARIABLES --------------------------------
//------------------------------------------------------------------------------

  // vector of pointers to current statobjects

  vector<statobject*> * statobj;

  // pointer to functions

  typedef void (* runpointer )(remlreg & b);

  runpointer functions[10];

  datamatrix D;
  vector<ST::string> modelvarnamesv;

  // global options

  fileoption outfile;
  optionlist globaloptions;

  // general variables for plotting results

  unsigned nrterms;                         // Anzahl Modellterme
  vector<unsigned> fullcondnr;              // Welches fullcond-Objekt gehört zu dem Term
  datamatrix catnr;                         // Zu welcher Kategorie gehört der Term
  vector<bool> needscat;                    // ist der Effekt kategorienspezifisch

//------------------------------------------------------------------------------
// ------------------------- for method drawmap --------------------------------
//------------------------------------------------------------------------------

  modelStandard mdrawmap;

  simpleoption replace;
  simpleoption swapcolors;
  simpleoption nolegend;
  simpleoption color;
  stroption title2;
  stroption outfile4;
  doubleoption upperlimit;
  doubleoption lowerlimit;
  intoption nrcolors;
  stroption plotvar;
  simpleoption pcat;
  simpleoption drawnames;
  simpleoption hclcolors;

  optionlist drawmapoptions;

  use udrawmap;

  friend void drawmaprun(remlreg & b);

//------------------------------------------------------------------------------
// ----------------------------- for method plotnonp ---------------------------
//------------------------------------------------------------------------------

  vector<ST::string> resultfiles;

  modelStandard mplotnonp;

  stroption xlab;
  stroption ylab;
  stroption connect;
  intoption height;
  intoption width;
  doubleoption ylimtop;
  doubleoption ylimbottom;
  doubleoption ystep;
  doubleoption ystart;
  stroption levels;
  simpleoption median;
  stroption outfile2;
  stroption title;
  simpleoption replace2;
  doubleoption xlimtop;
  doubleoption xlimbottom;
  doubleoption xstep;
  doubleoption xstart;
  intoption linewidth;
  intoption fontsize;
  intoption pointsize;
  stroption linecolor;
  doubleoption titlescale;

  optionlist plotnonpoptions;

  use uplotnonp;

  friend void plotnonprun(remlreg & b);

  // for method texsummary

  modelStandard mtexsummary;

  optionlist texsummaryoptions;

  use utexsummary;

  friend void texsummaryrun(remlreg & b);


//------------------------------------------------------------------------------
// -----------------------  for method 'remlrun'  ------------------------------
//------------------------------------------------------------------------------

  remlest RE;
  remlest_multinomial RE_M;
  remlest_multinomial_catsp RE_M_catsp;
  remlest_ordinal RE_O;
  remlest_multistate RE_MSM;

  friend void remlrun(remlreg & b);

  vector<FULLCOND*> fullcond;       // Vector of pointers to full conditionals

  MCMCoptions generaloptions;

  bool create_data(datamatrix & weight);

  bool create_response(datamatrix & response, datamatrix & weight);

  // OPTIONS for method regress

  ST::string add_name;

  doubleoption level1;                 // Nominal level 1 of credible intervals
  doubleoption level2;                 // Nominal level 2 of credible intervals

  doubleoption reference;              // Options and variables for
  datamatrix cats;                     // multicategorical models
  vector<int> allcats;                 // categories including reference category (at last position)
  bool ismultinomial;

  intoption maxint;
  vector<ST::string> families;          // Response families
  stroption family;                     // specifies the response distribution
  stroption knots;                      // equidistant knots or non equidistant
  vector<ST::string> knotsdef;          // knots (P-splines)

  intoption maxit;                      // Properties of the estimation
  doubleoption lowerlim;
  doubleoption eps;
  doubleoption maxchange;               // maximal value allowed for relative changes
  doubleoption maxvar;                  // maximal value allowed for a variance

  simpleoption aiccontrol;              // control estimation of survival models by AIC

  simpleoption noconst;                 // exclude intercept
  simpleoption fisher;                  // request storage of the fisher information matrix

  simpleoption constlambda;             // fix values of the smoothing parameters at the starting values
  simpleoption constscale;              // fix values of the scale parameter at the starting value

  stroption leftint;                    // Cox Model: left interval boundary
  stroption lefttrunc;                  //            left truncation time
  int leftintpos;                       // Positionen der Variablen in der
  int lefttruncpos;                     // Datenmatrix

  stroption state;                      // multistate: aktueller Zustand
  int statepos;                         // Position in der Datenmatrix

  stroption binomweight;                // Gewichtung von Binomialverteilung
  int binomweightpos;                   // Position in der Datenmatrix

  stroption naindicator;
  datamatrix naind;

  stroption globalfrailty;              // multistate: Globale Frailty-Variable
  doubleoption gflambdastart;           // starting value for variance of the global frailty
  int gfrailtypos;                      // Position in der Datenmatrix

  optionlist regressoptions;

  // end: OPTIONS for method regress

  vector<basic_termtype*> termtypes;
  modelterm modreg;
  vector<term> terms;

  use udata;

  bool resultsyesno;

//------------------------------------------------------------------------------
// -----------------------  for method 'mremlrun'  -----------------------------
//------------------------------------------------------------------------------

  friend void mremlrun(remlreg & b);

  modeltermmult modregmult;
  vector < vector <term> > termsmult;

  unsigned nrtransitions;             // no. of possible transitions
  vector<unsigned> nrfullconds;       // no. of fullconds per transition

//------------------------------------------------------------------------------
//----------------------------- Model terms ------------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//--------------------------- for the offset -----------------------------------
//------------------------------------------------------------------------------

  term_offset offset;
  bool create_offset(datamatrix & o);

//------------------------------------------------------------------------------
//-------------------------  for fixed effects ---------------------------------
//------------------------------------------------------------------------------

  FULLCOND_const * interceptpointer;
  vector<FULLCOND_const> fcconst;

  basic_termtype fixedeffects;
  term_fixed_catspecific fixed_catsp;

  bool create_const(const unsigned & colllinpred=0);

//------------------------------------------------------------------------------
//------------------------ for nonparametric terms -----------------------------
//------------------------------------------------------------------------------

  vector<FULLCOND_nonp_gaussian> fcnonpgaussian;

  term_autoreg_remlreg nonprw1rw2;
  bool create_nonprw1rw2(const unsigned & collinpred=0);

  term_autoreg_varcoef_remlreg nonprw1rw2_varcoef;
  bool create_nonprw1rw2_varcoef(const unsigned & collinpred=0);

  term_season_remlreg nonpseason;
  bool create_nonpseason(const unsigned & collinpred=0);

  term_season_varcoef_remlreg nonpseason_varcoef;
  bool create_nonpseason_varcoef(const unsigned & collinpred=0);

  term_spatial_remlreg nonpspatial;
  bool create_spatial(const unsigned & collinpred=0);

  term_spatial_varcoef_remlreg nonpspatial_varcoef;
  bool create_spatial_varcoef(const unsigned & collinpred=0);

//------------------------------------------------------------------------------
//-------------------------- for p-spline terms --------------------------------
//------------------------------------------------------------------------------

  vector<spline_basis> fcpspline;
  vector<spline_basis_surf> fcpsplinesurf;

  term_pspline_remlreg nonppspline;
  bool create_pspline(const unsigned & collinpred=0);

  term_varcoeff_pspline_remlreg nonpvarcoeffpspline;
  bool create_varcoeffpspline(const unsigned & collinpred=0);

  term_interactpspline_remlreg nonpinteractpspline;
  bool create_interactionspspline(const unsigned & collinpred=0);

  term_interactpspline_varcoeff_remlreg nonpvarcoeffinteractpspline;
  bool create_varcoeffinteractionspspline(const unsigned & collinpred=0);

  term_geospline_remlreg nonpgeospline;
  bool create_geospline(const unsigned & collinpred=0);

  term_geospline_varcoeff_remlreg nonpvarcoeffgeospline;
  bool create_geospline_varcoeff(const unsigned & collinpred=0);

//------------------------------------------------------------------------------
//-------------------------- for kriging terms ---------------------------------
//------------------------------------------------------------------------------

  vector<FULLCOND_kriging> fckriging;

  term_kriging_remlreg nonpspatial_kriging;
  bool create_kriging(const unsigned & collinpred=0);

  term_kriging_1dim_remlreg nonp_kriging;
  bool create_kriging_1dim(const unsigned & collinpred=0);

  term_kriging_varcoeff_remlreg nonpspatial_kriging_varcoeff;
  bool create_kriging_varcoeff(const unsigned & collinpred=0);

  term_geokriging_remlreg nonpspatial_geokriging;
  bool create_geokriging(const unsigned & collinpred=0);

  term_geokriging_varcoeff_remlreg nonpspatial_geokriging_varcoeff;
  bool create_geokriging_varcoeff(const unsigned & collinpred=0);

//------------------------------------------------------------------------------
//-------------------------- for baseline terms --------------------------------
//------------------------------------------------------------------------------

  vector<baseline_reml> fcbaseline;
  vector<baseline_reml> fcbaseline_varcoeff;

  term_baseline_remlreg nonp_baseline;
  bool create_baseline(const unsigned & collinpred=0);

  term_baseline_varcoeff_remlreg nonp_baseline_varcoeff;
  bool create_baseline_varcoeff(const unsigned & collinpred=0);

//------------------------------------------------------------------------------
//------------------------ for random effect terms -----------------------------
//------------------------------------------------------------------------------

  vector<FULLCOND_random> fcrandom;

  term_random_remlreg randomeff;
  bool create_random(const unsigned & collinpred=0);

  term_randomslope_remlreg randomeffslope;
  bool create_randomslope(const unsigned & collinpred=0);


  public:


//------------------------------------------------------------------------------
//--------------------------- PUBLIC FUNCTIONS ---------------------------------
//------------------------------------------------------------------------------

  // DEFAULT CONSTRUCTOR

  remlreg (void) : statobject()
         {
         type = "remlreg";
         resultsyesno = false;
         }

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - name = n

  remlreg (
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adb,
  #endif
  const ST::string & n,ofstream * lo,istream * i,
                                ST::string p,vector<statobject*> * st);

  // COPY CONSTRUCTOR

  remlreg (const remlreg & b);

  // DESTRUCTOR

  ~remlreg()
         {
         }


  // OVERLOADED ASSIGNMENT OPERATOR

  const remlreg & operator=(const remlreg & b);

  int parse(const ST::string & c);

  void describe(const optionlist & globaloptions = optionlist());

  };

#if defined (__BUILDING_GNU)
// ------------------ forward friend decls --------------------------

void drawmaprun(remlreg & b);
void plotnonprun(remlreg & b);
void mremlrun(remlreg & b);
void remlrun(remlreg & b);
void texsummaryrun(remlreg & b);
#endif

#endif





