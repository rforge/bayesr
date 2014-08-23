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



#if !defined (MCMCnonpgaussian_INCLUDED)

#define MCMCnonpgaussian_INCLUDED

#include"../export_type.h"
#include "mcmc_nonpbasis.h"
#include "statmat_penalty.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//----------------------- class: FULLCOND_nonp_gaussian ------------------------
//------------------------------------------------------------------------------

enum updatetype {gaussian,iwls,iwlsmode,hyperblock,hyperblockmode,gaussianlaplace};

class __EXPORT_TYPE FULLCOND_nonp_gaussian : public FULLCOND_nonp_basis
  {

  protected:

  updatetype utype;

  double lambdaprop;
  double lambda_prec;
  double a_invgamma;
  double b_invgamma;
  bool lambdaconst;
  bool Laplace;

  bool notransform;                         // Results will not be transformed
                                            // for multiplicative random effects


  datamatrix delta;
  vector< vector<unsigned> > neighbors;

// For REML
  datamatrix X_VCM;                     // für REML VCM
  datamatrix Z_VCM;                     // für REML VCM

  datamatrix remlspatialdesign;

//End: For REML

  double betaKbeta;

  unsigned updateW;

  datamatrix tildey;
  datamatrix weightiwls;

  datamatrix data2;

  FULLCOND_const * fcconst;


  datamatrix betaold;


  bool updatelinpred;

  envmatdouble XXenv;
  envmatdouble precenv;

  datamatrix mu;
  datamatrix muy;

  datamatrix betahelp;
  datamatrix diff;

  datamatrix betamode;
  datamatrix betamodeold;

  ST::string mapname;

  datamatrix xyvalues;
  bool lattice;
  ST::string pathmap;


  // FUNCTION: make_categories
  // TASK: devides the data in moddata into categories, maximal number of
  //       categories is 'maxnrint'
  //       Initialices 'index', 'posbeg', 'posend', 'weight' and
  //                   'effectvalues'
  // ADDITIONAL NOTE: Implementation is independent of the type of MRF

  void make_categories(const datamatrix & moddata,const unsigned & maxintnr);

  // FUNCTION: compute_XWX_env
  // TASK: Computes X'WX where W are the weights specified in weightmat
  //       stores the result in XX_env

  void compute_XWX_env(const datamatrix & weightmat,const unsigned & col=0);

  // FUNCTION: compute_XWX_varcoeff_env
  // TASK: Computes X'WX where W are the weights specified in weightmat
  //       stores the result in XX_env

  void compute_XWX_varcoeff_env(const datamatrix & weightmat,
                                const unsigned & col=0);

  // FUNCTION: compute_muy
  // TASK: computes W_i*(tildey_i+f_i)         (additive models)
  //                W_i*(tildey_i+f_i)*data_i  (varying coefficient models)
  //        stores the result in muy

  void compute_muy(double * workbeta);

  // FUNCTION: compute_XWX_XWtildey_env
  // TASK: computes X'WX (stored in XXenv) and X'W(tildey+Xbeta) (stored in muy)

  void compute_XWX_XWtildey_env(
  const datamatrix & weightmat,const datamatrix & tildey,double * workbeta,
  const unsigned & col);

  // FUNCTION: compute_XWX_XWtildey_varcoeff_env
  // TASK: computes X'WX (stored in XXenv) and X'W(tildey+Xbeta) (stored in muy)

  void compute_XWX_XWtildey_varcoeff_env(
  const datamatrix & weightmat,const datamatrix & tildey,double * workbeta,
  const unsigned & col);

  void update_linpred_diff(datamatrix & b1,datamatrix & b2);

  void update_linpred_diff(const unsigned & beg,const unsigned & end,const double & beta);

  double scale_proposal(void);

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_nonp_gaussian(void) : FULLCOND_nonp_basis()
    {
    }

  // CONSTRUCTOR 1  (for additive models)
  // o    : pointer to MCMCoptions object
  // dp   : pointer to distribution object
  // minb : Minimum blocksize (minblock)
  // maxb : Maximum blocksize (maxblock)
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored

  FULLCOND_nonp_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         const datamatrix & d,
                         FULLCOND_const * fcc,
                         const unsigned & maxint, const fieldtype & ft,
                         const ST::string & ti,
                         const ST::string & fp, const ST::string & pres,
                         const unsigned & c,const double & l,
                         const unsigned & per=12);

  // CONSTRUCTOR 2  (for spatial covariates)
  // o    : pointer to MCMCoptions object
  // dp   : pointer to distribution object
  // minb : Minimum blocksize (minblock)
  // maxb : Maximum blocksize (maxblock)
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored

  FULLCOND_nonp_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         const datamatrix & d, FULLCOND_const * fcc,
                         const MAP::map & m, const ST::string & mn,
                         const ST::string & ti,
                         const ST::string & fp, const ST::string & pres,
                         const unsigned & c,const double & l);

  // CONSTRUCTOR 3  (for varying coefficients models)
  // o    : pointer to MCMCoptions object
  // dp   : pointer to distribution object
  // minb : Minimum blocksize (minblock)
  // maxb : Maximum blocksize (maxblock)
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored

  FULLCOND_nonp_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         const datamatrix & d,const datamatrix & intvar,
                         FULLCOND_const * fcc,const unsigned & maxint,
                         const fieldtype & ft, const ST::string & ti,
                         const ST::string & fp, const ST::string & pres,
                         const unsigned & c,const double & l,bool  ce=false,
                         const unsigned & per=12);

  // CONSTRUCTOR 4  (spatial effect with x and y)
  // o    : pointer to MCMCoptions object
  // dp   : pointer to distribution object
  // minb : Minimum blocksize (minblock)
  // maxb : Maximum blocksize (maxblock)
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored

  FULLCOND_nonp_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         const datamatrix & dx,const datamatrix & dy,
                         FULLCOND_const * fcc,
                         const double & lambda, const double & md,
                         const ST::string & mp, const ST::string & ti,
                         const ST::string & fp,
                         const ST::string & pres,const ST::string & pmap,
                         const unsigned & c);

  // CONSTRUCTOR 5  (for spatial covariates, varying coefficients)
  // o    : pointer to MCMCoptions object
  // dp   : pointer to distribution object
  // minb : Minimum blocksize (minblock)
  // maxb : Maximum blocksize (maxblock)
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored

  FULLCOND_nonp_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,
                         const MAP::map & m, const ST::string & mn,
                         const datamatrix & d1, // effect modifier (map)
                         const datamatrix & d2, // interacting variable
                         const ST::string & ti,
                         const ST::string & fp, const ST::string & pres,
                         const unsigned & c,const double & l,bool ce);


  //------------------------ REML CONSTRUCTOR(s) -------------------------------

  // CONSTRUCTOR 6: REML: RW1, RW2, seasonal

  // d      : datamatrix
  // maxint : Maximum number of intervals
  // ft     : Field type (e.g RW1,Rw2, etc.)
  // ti     : Title of fullcond
  // pres   : path for results
  // l      : Smoothing parameter
  // sl     : Starting value for smoothing parameter
  // per    : Period of seasonal effect

  // Nach der Initalisierung stehen folgende Variablen zur Verfügung:

  // lambda: Glättungsparameter
  // period: Periode eine möglichen Saisoneffekts
  // title:  Titel der Full conditional
  // type: Type (z.B. RW1, RW2 etc.), vom Typ fieldtype, in mcmc_nonpbasis.h
  //       definiert.
  // plotstyle: Visualisierbar ja/nein. Werte sind plotnonp oder drawmap
  // index: Indexmatrix (Anzahl Beobachtungen x 1)
  // posbeg: Vektor, der für die Kategorien
  // posend:
  // weight: Differenzen benachbarter Kovariablenwerte
  // effectvalues: Verschiedenen Kovariablenwerte (als vector von string)
  // effectvaluesdouble: Verschiedenen Kovariablenwerte (als vector von doubles)

  // nach init_names:
  // datanames: vector von strings mit Variablennamen
  //            erst mal nur datanames[0]

  FULLCOND_nonp_gaussian(MCMCoptions * o,const datamatrix & d,
                         const unsigned & maxint, const fieldtype & ft,
                         const ST::string & ti, const ST::string & pres,
                         const double & l, const double & sl,
                         const bool & catsp, const unsigned & per=12);

  // CONSTRUCTOR 7: REML: RW1, RW2, seasonal (VCM)

  FULLCOND_nonp_gaussian(MCMCoptions * o,
                        const datamatrix & d1, const datamatrix & d2,
                        const unsigned & maxint,
                        const fieldtype & ft,const ST::string & ti,
                        const ST::string & pres,const double & l,
                        const double & sl,
                        const bool & catsp, const bool & ctr,
                        const unsigned & per=12
                        );

  // Constructor 8: REML spatial

  FULLCOND_nonp_gaussian(MCMCoptions * o,
                         const datamatrix & d,const MAP::map & m,
                         const ST::string & mn,const ST::string & ti,
                         const ST::string & fp, const ST::string & pres,
                         const double & l, const double & sl,
                         const bool & catsp);

  // Constructor 9: REML spatial (VCM)

  FULLCOND_nonp_gaussian(MCMCoptions * o, const datamatrix & d1,
                         const datamatrix & d2,const MAP::map & m,
                         const ST::string & mn,const ST::string & ti,
                         const ST::string & fp, const ST::string & pres,
                         const double & l, const double & sl,
                         const bool & catsp, const bool & ctr);

  void set_IWLS(const unsigned & uw,bool mode=false);

  void set_IWLS_hyperblock(const unsigned & uw,const double & ai,
                           const double & bi,bool mode=false);

  void set_lambdaconst(double la);

  void set_Laplace(void)
    {
    Laplace = true;
    utype = gaussianlaplace;
    betaold = beta;
    tildey=datamatrix(likep->get_nrobs(),1);
    weightiwls=datamatrix(likep->get_nrobs(),1,0);
    }

  void set_deltadim(const unsigned & d)
    {
    delta = datamatrix(d,1,1.0);
    }

  void set_delta(const datamatrix & d)
    {
    delta.assign(d);
    }

  // COPY CONSTRUCTOR

  FULLCOND_nonp_gaussian(const FULLCOND_nonp_gaussian & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_nonp_gaussian & operator=(const FULLCOND_nonp_gaussian & fc);

  const envmatdouble & get_K(void)
    {
    return Kenv;
    }

  void setK(unsigned i, unsigned j, double t)
    {
    Kenv.set(i,j,t);
    }

  double get(unsigned i, unsigned j)
    {
    return Kenv(i,j);
    }

  double compute_squareddiff(unsigned i,unsigned j)
    {
    double diff = beta(i,0)-beta(j,0);

    return (diff*diff)/sigma2;
    }

  double compute_squareddiff2(unsigned i,unsigned j)
    {
    double diff = beta(i,0)-beta(j,0);

    return (diff*diff);
    }

  double compute_quadform(void)
    {
    return Kenv.compute_quadform(beta,0);
    }

  double compute_fabsdiff(unsigned i,unsigned j)
    {
    double diff = beta(i,0)-beta(j,0);

    return fabs(diff)/sigma2;
    }

  double compute_sumfabsdiff(void)
    {
    return Kenv.compute_sumfabsdiff(beta,0);
    }

  // FUNCTION: updateK
  // TASK: updates the penalty matrix K

  void updateK(const datamatrix & q)
    {
    FULLCOND_nonp_basis::updateK(q);
    }

  // FUNCTION: compute_u
  // TASK: computes sum { u_t^2 }

  void compute_u(datamatrix & u)
    {
    FULLCOND_nonp_basis::compute_u(u);
    }

  // FUNCTION: update_linearpred
  // TASK: adds (add=true) or substracts the current function
  //       (additive model term) or the function*interacting
  //        variable (varying coefficients)

  void update_linpred(const bool & add);

  void update_linpred_current(const bool & add);

  // FUNCTION: init_data_varcoeff
  // TASK: initializes data and data2 (data^2) for varying coefficient model

  void init_data_varcoeff(const datamatrix & intvar, double add=0);


  void update(void);

  void update_IWLS(void);

  void update_IWLS_mode(void);

  void update_IWLS_hyperblock(void);

  void update_IWLS_hyperblock_mode(void);

  void update_gaussian_laplace(void);

  void update_gaussian_gemanreynolds(void);

  void update_lambdaconst(void);

  bool posteriormode(void);

  bool posteriormode_converged(const unsigned & itnr);

  void outresults(void);

  void get_effectmatrix(datamatrix & e,vector<ST::string> & enames,
                        unsigned be, unsigned en,effecttype t);


  unsigned get_nreffects(effecttype t);


  void outoptions(void);

  ST::string getinfo(void)
    {
    if (type == MCMC::mrf)
      return mapname;
    else
      return title;
    }

  void init_name(const ST::string & na);

  void init_names(const vector<ST::string> & na);

  void init_priorassumptions(const ST::string & na);

  void set_notransform(void)
    {
    notransform=true;
    }


  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND::reset();
    sigma2 = 10;
    }

  // DESTRUCTOR

  ~FULLCOND_nonp_gaussian() {}


  //------------------------------- FOR REML -----------------------------------

  void createreml(datamatrix & X,datamatrix & Z,
                                     const unsigned & Xpos,
                                     const unsigned & Zpos);


  double outresultsreml(datamatrix & X,datamatrix & Z,
                                     datamatrix & betareml,datamatrix & betacov,
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

  };


} // end: namespace MCMC

#endif

