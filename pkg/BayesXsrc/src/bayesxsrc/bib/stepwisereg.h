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




#if !defined (stepwisereg_INCLUDED)

#define stepwisereg_INCLUDED

#include"../export_type.h"
#include"statobj.h"
#include"dataobj.h"
#include"map.h"
#include"mapobject.h"

#include"mcmc.h"

#include"distribution.h"
//#include"nbinomial.h"

#include"mcmc_const_stepwise.h"

#include"fullcond_nonp_gaussian_stepwise.h"
#include"variance_nonp.h"
#include"fullcond_surf_gaussian.h"

#include"fullcond_pspline_stepwise.h"
#include"mcmc_pspline_surf.h"
#include"fullcond_pspline_surf_stepwise.h"
#include"mcmc_pspline.h"
//#include"fullcond_projection.h"

#include"randomeffect_stepwise.h"

#include"mcmcsimul2.h"
#include"mcmcsimul2_multi.h"

#include"model_stepwise.h"


using MCMC::MCMCoptions;
using MCMC::DISTRIBUTION;
using MCMC::DISTRIBUTION_gaussian;
using MCMC::DISTRIBUTION_binomial;
using MCMC::DISTRIBUTION_binomial_latent;
using MCMC::DISTRIBUTION_poisson;
using MCMC::DISTRIBUTION_gamma2;
using MCMC::DISTRIBUTION_vargaussian;
//using MCMC::DISTRIBUTION_nbinomial;
using MCMC::DISTRIBUTION_multinom2;
//using MCMC::DISTRIBUTION_multinomial_latent;
using MCMC::DISTRIBUTION_cumulative_latent3;
using MCMC::FULLCOND;
using MCMC::FULLCOND_const;
using MCMC::FULLCOND_const_stepwise;
using MCMC::FULLCOND_const_gaussian_special;
using MCMC::FULLCOND_nonp_gaussian_stepwise;
using MCMC::FULLCOND_variance_nonp;
using MCMC::FULLCOND_pspline;
using MCMC::FULLCOND_pspline_stepwise;
using MCMC::FULLCOND_pspline_surf;
using MCMC::FULLCOND_pspline_surf_stepwise;
using MCMC::FULLCOND_random_stepwise;
//using MCMC::FULLCOND_projection;
using MCMC::STEPWISErun;
using MCMC::STEPMULTIrun;


class __EXPORT_TYPE stepwisereg : public statobject
  {


  private :

  ST::string pathres;
  ST::string title;
  ST::string pathnonp;

  void make_paths(unsigned collinpred,ST::string & pathnonp,
                          ST::string & pathres,ST::string & title,
                          ST::string  varname1,ST::string  varname2,
                          ST::string  endingraw,
                          ST::string  endingres,ST::string  endingtitle) ;

  bool check_gaussian(void);

  bool check_nongaussian(void);

  void clear(void);
  void initpointers(void);

  //------------------------- PRIVATE VARIABLES --------------------------------

  // vector of pointers to current statobjects

  vector<statobject*> * statobj;

  // pointer to functions

  typedef void (* runpointer )(stepwisereg & b);

  runpointer functions[10];

  datamatrix D;
  vector<ST::string> modelvarnamesv;

  // global options

  fileoption outfile;
  optionlist globaloptions;

  // for stepwise

  stroption algorithm;

  stroption procedure;

  stroption minimum;

  stroption criterion;
  doubleoption gcvfactor;

  doubleoption proportion;

  intoption steps;

  stroption trace;

  intoption number;

  stroption startmodel;

  intoption increment;

  intoption bootstrap;
  simpleoption unconditional;
  intoption setseed;

  stroption CI;
  intoption iterations;                // Number of iterations
  intoption burnin;                    // Number of burnin iterations
  intoption step;                      // Thinning parameter
  doubleoption level1;
  doubleoption level2;

  simpleoption hier;

  intoption maxint;

  vector<ST::string> outfiles;

  vector<FULLCOND*> fullcond;       // Vector of pointers to full conditionals
  STEPWISErun runobj;
  STEPMULTIrun runobjm;

  vector<MCMCoptions> generaloptions;
  bool create_generaloptions(void);

  // OPTIONS for method regress

  ST::string add_name;

  bool varianceest;

  simpleoption constscale;

  // options gamma distributed response
  stroption scalegamma;
  doubleoption scalevalue;
  doubleoption gamvar;
  intoption cit;
  // options gamma distributed response

  // options negative binomial distributed response
  doubleoption propvar;
  stroption distopt;
  stroption propopt;
  simpleoption hierarchical;
  // options negative binomial distributed response

  // options for cumprobit
  simpleoption nosort;
  // end: options for cumprobit

  simpleoption predict;                 // indicates that predicted values,
                                        // deviances, etc. should be computed
  simpleoption predictmu;
  intoption predictuntil;

  vector<ST::string> families;          // Response families
  stroption family;                     // specifies the response distribution

  stroption knots;                      // equidistant knots or non equidistant
                                        // knots (P-splines)
  vector<ST::string> knotsdef;


  optionlist regressoptions;


  vector<unsigned> begin_fc;
  vector<unsigned> end_fc;


  friend void regressrun(stepwisereg & b);

  friend void mregressrun(stepwisereg & b);

  // end for method stepwise

  // for method drawmap

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

  optionlist drawmapoptions;

  use udrawmap;

  friend void drawmaprun(stepwisereg & b);

  // for method plotnonp

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
  stroption title0;
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

  friend void plotnonprun(stepwisereg & b);

  // for method texsummary

  modelStandard mtexsummary;

  optionlist texsummaryoptions;

  use utexsummary;

  friend void texsummaryrun(stepwisereg & b);


  // ---------------------  for method 'getsample' -----------------------------

  optionlist getsampleoptions;

  modelStandard mgetsample;

  useDataset usegetsample;

  friend void __EXPORT_TYPE getsamplerun(stepwisereg & b);


//------------------------------ DISTRIBUTION ----------------------------------

  vector<ST::string> distrstring;
  vector<unsigned> distrposition;

  unsigned nrcategories;           // number of categories of the response

  vector<DISTRIBUTION_gaussian> distr_gaussian;
  DISTRIBUTION_binomial distr_binomial;
  DISTRIBUTION_binomial_latent distr_binomlat;
  DISTRIBUTION_poisson distr_poisson;
  DISTRIBUTION_gamma2 distr_gamma;
  DISTRIBUTION_vargaussian distr_vargaussian;
//  DISTRIBUTION_nbinomial distr_nbinomial;
  //DISTRIBUTION_multgaussian distr_multgaussian;
  DISTRIBUTION_multinom2 distr_multinom;
//  DISTRIBUTION_multinomial_latent distr_multinom_latent;
  DISTRIBUTION_cumulative_latent3 distr_cumlat3;

  doubleoption reference;

  vector<DISTRIBUTION *> distr;              // Pointer to distribution objects


  bool create_distribution(ST::string method);

  vector<basic_termtype*> termtypes;
  modelterm modreg;
  vector<term> terms;

  // for multivariate regression

  modeltermmult modregmult;
  vector < vector <term> > termsmult;

  use udata;

  bool resultsyesno;
  bool bootyesno;
  bool hierarchical_model_yesno;

//--------------------------- for the offset -----------------------------------

  term_offset offset;
  bool create_offset(datamatrix & o);

//-------------------------  for fixed effects ---------------------------------

  vector<FULLCOND_const_stepwise> factor;
  vector<FULLCOND_const_gaussian_special> normalconst_special;
  vector<FULLCOND_const_stepwise> normalconst;
//  vector<FULLCOND_const_gamma> gammaconst;
  FULLCOND_const * fcconst_intercept;

  basic_termtype fixedeffects;

  term_factor_stepwise termfactor;
  term_nonlinearf_stepwise termnonlinearf;

  bool create_const(const unsigned & collinpred=0);
  bool create_factor(const unsigned & collinpred=0);
  bool create_nonlinearf(const unsigned & collinpred=0);

// ----------------------- end: for fixed effects ------------------------------

//------------------------ for nonparametric terms -----------------------------

  vector<FULLCOND_nonp_gaussian_stepwise> fcnonpgaussian;
  term_autoreg_stepwise nonprw1rw2;
  term_season_stepwise nonpseason;

  bool create_nonprw1rw2(const unsigned & collinpred=0);
  bool create_nonpseason(const unsigned & collinpred=0);

  term_spatial_stepwise nonpspatial;

  bool create_spatial(const unsigned & collinpred=0);

  vector<FULLCOND_pspline> fcpspline;
  //vector<FULLCOND_pspline_gaussian> fcpsplinegaussian;
  vector<FULLCOND_pspline_stepwise> fcpsplinestep;
  //vector<IWLS_pspline> fciwlspspline;
  vector<FULLCOND_pspline_surf> fcpsplinesurf;
  vector<FULLCOND_pspline_surf_stepwise> fcpsplinesurfstep;
  term_pspline_stepwise nonppspline;
  bool create_pspline(const unsigned & collinpred=0);

  term_interactpspline_stepwise nonpinteractpspline;                //neu
  term_geospline_stepwise nonpgeospline;                            //neu

  bool create_interactionspspline(const unsigned & collinpred=0);
  bool create_geospline(const unsigned & collinpred=0);

//  vector<FULLCOND_projection> fcprojection;
//  term_projection_stepwise termprojection;
//  bool create_projection(const unsigned & collinpred=0);

//------------------------ for nonparametric terms -----------------------------

//------------------------- for random effects ---------------------------------

  term_random_stepwise randomeff;
  term_randomslope_stepwise randomeffslope;

  vector<FULLCOND_random_stepwise> fcrandomgaussian;

  bool create_random(const unsigned & collinpred=0);
  bool create_randomslope(const unsigned & collinpred=0);


  //-------------------- end: for random effects -------------------------------

  void create(void);

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_pointer * adminp_p;
  #endif

  public:

  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR

  stepwisereg (void) : statobject()
         {
         type = "stepwisereg";
         resultsyesno = false;
         }

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - name = n


  #if defined(JAVA_OUTPUT_WINDOW)
  stepwisereg (administrator_basic * adb, administrator_pointer * adp,
               const ST::string & n,ofstream * lo,istream * i,
                                ST::string p,vector<statobject*> * st);
  #else
  stepwisereg (const ST::string & n,ofstream * lo,istream * i,
                                ST::string p,vector<statobject*> * st);
  #endif

  // COPY CONSTRUCTOR

  stepwisereg (const stepwisereg & b);

  // DESTRUCTOR

  ~stepwisereg()
         {
         }


  // OVERLOADED ASSIGNMENT OPERATOR

  const stepwisereg & operator=(const stepwisereg & b);

  int parse(const ST::string & c);

  void describe(const optionlist & globaloptions = optionlist());


  };

#if defined (__BUILDING_GNU)
// ----------------------- forward friend decls ---------------------

void __EXPORT_TYPE getsamplerun(stepwisereg & b);
void plotnonprun(stepwisereg & b);
void texsummaryrun(stepwisereg & b);
void drawmaprun(stepwisereg & b);
void regressrun(stepwisereg & b);
void mregressrun(stepwisereg & b);
#endif

#endif





