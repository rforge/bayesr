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



#if !defined (FULLCOND_INCLUDED)

#define FULLCOND_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"Random.h"
#if defined(MICROSOFT_VISUAL)
#include<limits>
#else
#include"../values.h"
#endif
#include<fstream>
#include<vector>
#include<bitset>
#include"mcmc.h"
#include"clstring.h"
#include<cmath>

namespace MCMC
{

using randnumbers::rand_normvek;
using randnumbers::rand_invgamma;
using randnumbers::uniform;
using randnumbers::rand_normal;
using std::vector;
using std::bitset;
using std::ofstream;

//------------------------------------------------------------------------------
//--------------------------- CLASS: FULLCOND ----------------------------------
//------------------------------------------------------------------------------


const int flagnr = 3;

const std::bitset<flagnr> nosamples(1ul);
//const std::bitset<flagnr> nosamples(001);   // samples should not be stored in
                                              // a file
const std::bitset<flagnr> norelchange(2ul);
//const std::bitset<flagnr> norelchange(010); // relative changes of parameters
                                              // should not be printed
const std::bitset<flagnr> nooutput(4ul);
//const std::bitset<flagnr> nooutput(100);    // acceptance rates and title
                                              // should not be printed

enum covstyle {covariance,precision,correlation};

enum plotstyles {noplot,plotnonp,drawmap,drawmapgraph};

enum fullcondtype {fixed,variance,nonparametric,spatial,randomeffects,
                   season,randomslopes,factor,nonlinearf};

enum effecttype {current,mean,median,fvar_current,fvar_mean,fvar_median};

class __EXPORT_TYPE FULLCOND
  {

  protected:

  fullcondtype fctype;

  MCMCoptions * optionsp;        // Pointer to general MCMC options

  ST::string title;              // Title/name of the full conditional

  ST::string samplepath;         // filename for storing sampled parameters
  ofstream samplestream;         // stream object for storing sampled parameters

  ST::string pathresult;         // filename for storing results;
  ST::string pathresult2;
  ST::string pathresult3;

  ST::string pathcurrent;
  ST::string pathcurrent2;
  ST::string pathcurrent3;

  plotstyles plotstyle;

  ST::string term_symbolic;

  vector<ST::string> priorassumptions;

  vector<ST::string> results_latex;
  ST::string results_type;

  unsigned fcnumber;

  std::bitset<flagnr> flags;          // some flags see bitsets above


  datamatrix data;                   // Matrix of Covariates
  vector<ST::string> datanames;      // Names of Covariates

  vector<ST::string> errors;


  unsigned nrpar;                // total number of parameters
                                 // (beta.rows()*beta.cols())
  datamatrix beta;               // Matrix of current parameters
  datamatrix beta_mode;          // for posteriormode, the parameters
                                 // of the last iteration of iwls
                                 // used to check if convergence has been
                                 // already achieved
  datamatrix betamean;           // Sampling mean of parameters
  datamatrix betas2;             // Sampling sum of squares of parameters
  datamatrix betavar;            // Sampling variance of parameters
  datamatrix betamin;            // Sampling minimum
  datamatrix betamax;            // Sampling maximum

  double level1;                 // level1 of credible intervals
                                 // DEFAULT: 95 %
  double level2;                 // level2 of credible intervals
                                 // DEFAULT: 80 %
  double lower1;
  double lower2;
  double upper1;
  double upper2;

  datamatrix betaqu_l1_lower;    // (100-level1)/2 percent quantile
  datamatrix betaqu_l2_lower;    // (100-level2)/2 percent quantile
  datamatrix betaqu50;           //  50 percent quantile
  datamatrix betaqu_l1_upper;    // corresponding upper quantile (level1)
  datamatrix betaqu_l2_upper;    // corresponding upper quantile (level2)

  datamatrix betameanold;        // Will be initialized after the burnin period
  datamatrix betavarold;
  datamatrix betaminold;
  datamatrix betamaxold;

  double transform;              // The factor with which all beta's will be
                                 // multiplied before storing them and computing
                                 // means, std, etc.
                                 // DEFAULT: transform = 1
  double addon;                  // An additive constant that will be added
                                 // on each component of beta before storing
                                 // DEFAULT: addon = 0;

  datamatrix transformmult;

  bool transformnonlinear;
  bool transformed;
  ST::string transformtype;

  //unsigned long acceptance;      // number of accepted iterations
  unsigned long nrtrials;        // number of trials

  bool identifiable;             // true, if term is identifiable
  bool center;                   // true, if beta should be centered
  bool baseline;

  unsigned column;               // the response category the fc belongs to

  vector<double> weight;

  //----------------------------------------------------------------------------
  //------------------------- Variables for stepwise ---------------------------
  //----------------------------------------------------------------------------

  double lambdastart;
  double lambdamin;
  double lambdamax;
  datamatrix data_forfixed;
  bool forced_into;
  double df_for_lambdamax;
  double df_for_lambdamin;
  double dfstart;
  //bool spfromdf;
  ST::string spfromdf;
  int number;
  bool df_equidist;
  double df_accuracy;
  bool inthemodel;    //gibt an, ob Fullc-Obj. im aktuellen Modell enthalten ist
  bool fixornot;      // gibt an, ob statt Fullc-Obj. fixer Effekt da ist
  vector<FULLCOND*> interactions_pointer;
  datamatrix betaright;
  datamatrix beta_average;
  bool calculate_xwx;
  bool calculate_xwx_vc;
  bool nofixed;
  bool kombimatrix;
  unsigned numberofmatrices;
  unsigned matrixnumber;

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  //-------------------------- Variables for REML ------------------------------
  //----------------------------------------------------------------------------

  unsigned dimX;
  unsigned dimZ;
  double startlambda;
  bool isnonparametric;
  bool catspecific;
  bool centervcm;

  //----------------------------------------------------------------------------
  //------------------------ end: Variables for REML ---------------------------
  //----------------------------------------------------------------------------

  public:

  unsigned long acceptance;      // number of accepted iterations

  // FUNCTION: get_fctype
  // TASK: returns the type of the full conditional

  fullcondtype get_fctype(void)
    {
    return fctype;
    }

  // FUNCTION: get_fctype
  // TASK: returns the type of the full conditional

  void set_fctype(const fullcondtype & fc)
    {
    fctype = fc;
    }

  void set_optionsp(MCMCoptions * o)
    {
    optionsp = o;
    }

  // FUNCTION: readsample
  // TASK: reads sample of parameter 'nr' and stores the sample in datamatrix
  //       sample (sample must be optionsp->samplesize x 1 matrix)
  //       nr >= 0 and nr < nrpar required

  void readsample(datamatrix & sample,const unsigned & nr,
                  const unsigned & col=0) const;


  // FUNCTION readsample
  // TASK: reads the nr-th stored parameter matrix of the sample
  //       'b' must have proper size, i.e. b.rows() == beta.rows() and
  //       b.cols() == b.cols()

  void readsample2(datamatrix & b,const unsigned & nr) const;


  // FUNCTION: readsample
  // TASK: reads ALL stored paramters of the sample
  //       'b' must have proper size, i.e. b.rows() = samplesize
  //       b.cols() = nrpar

  void readsample3(datamatrix & b) const;

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  FULLCOND(void);

  // CONSTRUCTOR
  // o    : pointer to MCMCoptions object
  // t    : title of the full conditional (for example "fixed effects")
  // rows : number of rows of the beta matrix (i.e. number of parameters)
  // cols : number of columns of the beta matrix
  //        (i.e. number of categories of the response variable)
  // fp   : file path for storing sampled parameters

  FULLCOND(MCMCoptions * o,const datamatrix & d,const ST::string & t,
           const unsigned & rows, const unsigned & cols,
           const ST::string & fp);

  // CONSTRUCTOR (REML)

  FULLCOND(MCMCoptions * o,const ST::string & t);

  // COPY CONSTRUCTOR

  FULLCOND(const FULLCOND & m);


  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND & operator=(const FULLCOND & m);


  // DESTRUCTOR

  ~FULLCOND()
    {
    if (flags[0] == 0)
      remove(samplepath.strtochar());
    }


  //----------------------------------------------------------------------------
  //------------------ FUNCTIONS FOR CENTERING BETA ----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: centerbeta
  // TASK: centers the current parametermatrix about the mean

  double centerbeta(void);

  double centerbeta2(datamatrix & sumx1,datamatrix & sumx2);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //------------------------ACCESS TO ERROR MESSAGES ---------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: geterrors
  // TASK: returns errormessages

  const vector<ST::string> & geterrors(void)
    {
    return errors;
    }

  // FUNCTION: outerrors
  // TASK: writes errors

  void outerrors(void);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //----------------- ACCSESS TO NAMES OF THE DATA VARIABLES -------------------
  //----------------------------------------------------------------------------

  // FUNCTION: init_names
  // TASK: initializes names of the covariates

  virtual void init_names(const vector<ST::string> & na)
    {
//    assert (na.size() == data.cols());
    datanames=na;
    }

  virtual void init_name(const ST::string & n)
    {
//    assert (data.cols() == 1);
    datanames = vector<ST::string>(1,n);
    }

  // FUNCTION: get_datanames
  // TASK: returns datanames

  const vector<ST::string> & get_datanames(void) const
    {
    return datanames;
    }


  virtual ST::string getinfo(void)
    {
    return title;
    }


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  // FUNCTION: get_col
  // TASK: returns the column the full conditional belongs to

  const unsigned & get_col(void) const
    {
    return column;
    }


  //--------------------------- ACCESS TO EFFECTS ------------------------------

  // FUNCTION: get_effectmatrix
  // TASK: returns true if full conditional contains predictor effects
  //       returns false if full conditional contains no predictor effects
  //       (e.g. variance parameters)
  //       if true, the effect will be stored in datamatrix e

  virtual void get_effectmatrix(datamatrix & e,vector<ST::string> & enames,
                                unsigned be,unsigned en, effecttype t)
    {
    }

  virtual unsigned get_nreffects(effecttype t)
    {
    return 0;
    }


  //------------- ACCESS TO PARAMETERS AND OTHER CHARACTERISTICS ---------------

   // FUNCTION: setbeta
   // TASK: initializes beta matrices (i.e. beta, betaold, etc)

   void setbeta(const unsigned & rows,const unsigned & cols,const double & v);

   void setbeta(const datamatrix & betanew);

   void setbetavalue(const unsigned & row,const unsigned & col,const double & v);

   // FUNCTION: getbetapointer
   // TASK: returns the pointer of the first element in beta

   double * getbetapointer(void)
     {
     return beta.getV();
     }

  // FUNCTION: getbeta
  // TASK: returns the current state 'beta' of the markov chain

  const datamatrix & getbeta(void) const
    {
    return beta;
    }


  // FUNCTION: getbeta
  // TASK: returns the current state 'beta' of the markov chain

  const double & getbeta(const unsigned & row,const unsigned & col) const
    {
    return beta(row,col);
    }

  // FUNCTION: get_betamean
  // TASK: returns 'betamean'

  const double & get_betamean(const unsigned & row,const unsigned & col) const
    {
    return betamean(row,col);
    }

  const datamatrix & get_betamean(void) const
    {
    return betamean;
    }

  double * get_betameanp(void)
    {
    return betamean.getV();
    }

  // FUNCTION: get_betavar
  // TASK: returns 'betavar'

  const double & get_betavar(const unsigned & row,const unsigned & col) const
    {
    return betavar(row,col);
    }


  double * get_betaavep(void)
    {
    return beta_average.getV();
    }


  double * get_betavarp(void)
    {
    return betavar.getV();
    }


  // FUNCTION: get_beta_lower1
  // TASK: returns 'beta_lower1'

  const double & get_beta_lower1(const unsigned & row,
                                 const unsigned & col) const
    {
    return betaqu_l1_lower(row,col);
    }


  double * get_beta_lower1_p(void)
    {
    return betaqu_l1_lower.getV();
    }


  // FUNCTION: get_beta_lower2
  // TASK: returns 'beta_lower2

  const double & get_beta_lower2(const unsigned & row,
                                 const unsigned & col) const
    {
    return betaqu_l2_lower(row,col);
    }


  double * get_beta_lower2_p(void)
    {
    return betaqu_l2_lower.getV();
    }


  // FUNCTION: get_betaqu50
  // TASK: returns 'betaq50'

  const double & get_betaqu50(const unsigned & row,const unsigned & col) const
    {
    return betaqu50(row,col);
    }

  double * get_betaqu50p(void)
    {
    return betaqu50.getV();
    }


  // FUNCTION: get_beta_upper1
  // TASK: returns 'beta_upper1'

  const double & get_beta_upper1(const unsigned & row,
                                 const unsigned & col) const
    {
    return betaqu_l1_upper(row,col);
    }


  double * get_beta_upper1_p(void)
    {
    return betaqu_l1_upper.getV();
    }


  // FUNCTION: get_beta_upper2
  // TASK: returns 'beta_upper2

  const double & get_beta_upper2(const unsigned & row,
                                 const unsigned & col) const
    {
    return betaqu_l2_upper(row,col);
    }


  double * get_beta_upper2_p(void)
    {
    return betaqu_l2_upper.getV();
    }

  const double & get_betamin(const unsigned & row, const unsigned & col) const
    {
    return betamin(row,col);
    }

  const double & get_betamax(const unsigned & row, const unsigned & col) const
    {
    return betamax(row,col);
    }

  double get_lower1(void)
    {
    return lower1;
    }

  double get_lower2(void)
    {
    return lower2;
    }

  double get_upper1(void)
    {
    return upper1;
    }

  double get_upper2(void)
    {
    return upper2;
    }

  double get_level1(void)
    {
    return level1;
    }

  double get_level2(void)
    {
    return level2;
    }

  // FUNCTION: get_nrpar
  // TASK: returns the total number of parameters (i.e. the number of elements
  //       in beta = beta.rows()*beta.cols()

  const unsigned & get_nrpar(void) const
    {
    return nrpar;
    }

  // FUNCTION: set_transform
  // TASK: sets transform

  void set_transform(const double & t)
    {
    transform = t;
    }


  // FUNCTION: set_addon

  void set_addon(const double & ao)
    {
    addon = ao;
    }


  // FUNCTION: set_transformmult
  // TASK: sets transformmult

  void set_transformmult(const datamatrix & tr)
    {
    transformmult=tr;
    }

  // FUNCTION: get_title
  // TASK: returns the name/title of the full conditional

  const ST::string & get_title(void)
    {
    return title;
    }


  void setflags(const bitset<flagnr> & newflags);

  void resetflags(void)
    {
    flags.reset();
    }

  bool is_identifiable(void)
    {
    return identifiable;
    }

  void set_center(const bool & c)
    {
    center = c;
    }

  bool is_baseline(void)
    {
    return baseline;
    }

  // FUNCTION: sample_stored
  // TASK: returns true, if sampled parameters are/will be stored in a file

  bool sample_stored(void)
    {
    if (flags[0] == 0)
      return true;
    else
      return false;
    }


  ST::string get_pathresult(void) const
    {
    return pathresult;
    }

  ST::string get_pathcurrent(void) const
    {
    return pathcurrent;
    }

  unsigned get_fcnumber(void) const
    {
    return fcnumber;
    }

  void set_fcnumber(const unsigned & n)
    {
    fcnumber = n;
    }


  // FUNCTION: get_errors
  // TASK: returns a vector of strings with errormessages
  //       if the vector is empty no errors during initialization occured

  const vector<ST::string> & get_errors(void)
    {
    return errors;
    }


  const plotstyles & get_plotstyle(void)
    {
    return plotstyle;
    }

  const ST::string & get_term_symbolic(void)
    {
    return term_symbolic;
    }

  const vector<ST::string> & get_priorassumptions(void)
    {
    return priorassumptions;
    }

  // FUNCTION: compute_autocorr
  // TASK: computes autocorrelation function for lags 1 - 'lag' for parameter
  //       beta(row,col). returns the result as  a column vector of
  //       autocorrelations

  datamatrix compute_autocorr(const unsigned & lag,const unsigned & row,
                              const unsigned & col) const;

  // FUNCTION: get_samples
  // TASK: stores the sampled parameters in ASCII format

  virtual void get_samples(const ST::string & filename, const unsigned & step = 1) const;

  // FUNCTION: get_covmatrix
  // TASK: computes the covariance matrix of the sampled parameters
  //       when MCMC simulation is completed

  void get_covmatrix(datamatrix & r);

  void get_covmatrix(const ST::string & file,const covstyle & st=covariance);

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  virtual void update(void);

  void updatemult(void);

  virtual void update_linpred(const bool & add)
    {
    }

  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  virtual bool posteriormode(void);

  virtual bool posteriormode_converged(const unsigned & itnr);

  void posteriormode_set_beta_mode(void);

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  virtual void outoptions(void)
    {
    }

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  virtual void outresults(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  virtual void reset(void);

   //---------------------------------------------------------------------------
   //---------------------- Handling of missing values -------------------------
   //---------------------------------------------------------------------------

  virtual bool is_missing(const ST::string & na)
    {
      return false;
    }

  virtual void update_missings(datamatrix & x,const datamatrix & linpred,
                               const statmatrix<unsigned> & w,
                               const ST::string & na,double & scalex)
    {

    }


  double get_data(const unsigned & row,const unsigned & col) const
    {
    return data(row,col);
    }

  void set_transform(ST::string & suffix,ST::string & trtype);

  bool transform_yes(void)
    {
    return transformnonlinear;
    }


  vector<ST::string> & get_results_latex(void)
    {
    return results_latex;
    }


  ST::string & get_results_type(void)
    {
    return results_type;
    }


  ST::string & get_samplepath(void)
    {
    return samplepath;
    }


  // ---------------------------------------------------------------------------
  // ------------------------ FOR STEPWISE SELECTION ---------------------------
  // ---------------------------------------------------------------------------

  void set_calculate_xwx(void)
    {
    calculate_xwx = true;
    calculate_xwx_vc = true;
    }

  virtual void set_nofixed(bool fix)
    {
    nofixed = fix;
    }

  virtual void create_weight(datamatrix & w)
    {
    }

  virtual double get_lambda(void)
    {
      return 0.1;
    }

  virtual void set_lambdaconst(double la)
    {
    }

  void set_inthemodel(double modell);
  void get_inthemodel(bool & drin, bool & fix);

  virtual void get_interactionspointer(vector<FULLCOND*> & inter)
    {
    }

  virtual void split_data(const vector<ST::string> & names)
    {
    }

  virtual void merge_data(void)
    {
    }

  virtual void hierarchical(ST::string & possible)  // neu! f\FCr hierarchisches Modell
    {
    }

  virtual void const_varcoeff(void)
    {
    }

  virtual void set_pointer_to_interaction(FULLCOND * inter)
    {
    }

  virtual void update_linold(void)
    {
    }

  virtual void update_linold_vc(void)
    {
    }

  virtual unsigned get_rankK2(void)
     {
     return 0;
     }

    // FUNCTION: compute_df
    // TASK: returns the approximate degrees of freedom of a smoother

  virtual double compute_df(void)
    {
    return 0.;
    }

  virtual void set_dfunstruct(const double & df_unstr)  // f\FCr spatialtotal
    {
    }

    // FUNCTION: set_stepwise_options
    // TASK: \FCbergibt die Optionen der einzelnen Funktionen

  void set_stepwise_options(double lstart, double lmax, double lmin, bool forced,
                            double df_lmax, double df_lmin, ST::string spdf,
                            double numb, bool df_equi)
     {
     lambdamin=lmin;
     lambdamax=lmax;
     lambdastart=lstart;
     forced_into = forced;
     df_for_lambdamax = df_lmax;
     df_for_lambdamin = df_lmin;
     spfromdf = spdf;
     number = numb;
     df_equidist = df_equi;
     }

  void set_dfstart(double df_start)
    {
    dfstart = df_start;
    }

  void set_stepwise_accuracy(double df_accu)
     {
     df_accuracy = df_accu;
     }

  void set_number(double numb)
     {
     number = numb;
     }

  void set_df_lambdamin(double & df_lmin)
    {
    df_for_lambdamin = df_lmin;
    }

  void set_df_lambdamax(double & df_lmax)
    {
    df_for_lambdamax = df_lmax;
    }

  double get_lambdastart(void)
     {
     return lambdastart;
     }

  double get_lambdamin(void)
     {
     return lambdamin;
     }

  double get_lambdamax(void)
     {
     return lambdamax;
     }

  bool get_forced(void)
     {
     return forced_into;
     }

  double get_df_lambdamax(void)
     {
     return df_for_lambdamax;
     }

  double get_df_lambdamin(void)
     {
     return df_for_lambdamin;
     }

  ST::string get_spfromdf(void)
     {
     return spfromdf;
     }

  double get_dfstart(void)
    {
    return dfstart;
    }

  double get_number(void)
     {
     return number;
     }

  bool get_df_equidist(void)
     {
     return df_equidist;
     }

  double get_accuracy(void)
     {
     return df_accuracy;
     }

  bool get_kombimatrix(void)
    {
    return kombimatrix;
    }

  unsigned get_numberofmatrices(void)
    {
    return numberofmatrices;
    }

  unsigned get_matrixnumber(void)
    {
    return matrixnumber;
    }

  void set_matrixnumber(unsigned mno)    // gibt an, f\FCr welche Strafmatrix das fullcond-Objekt zust\E4ndig ist
    {
    matrixnumber = mno;
    }

  // FUNCTION: update_stepwise
  // TASK: returns (usually) the current smoothing parameter

  virtual void update_stepwise(double la)
    {
    }

  // FUNCTION: get_effect
  // TASK: returns a string of the estimated effect

  virtual ST::string get_effect(void)
    {
    return "";
    }

  // FUNCTION: reset_effect
  // TASK: resets the effect, subtracts the current effect from linearpred

  virtual void reset_effect(const unsigned & pos)
    {
    }

  virtual void set_effect_zero(void)
    {
    }


  // FUNCTION: include_effect
  // TASK: includes an effect at position 'pos' with name 'name'

  virtual void include_effect(const vector<ST::string> & name, const datamatrix & xnew)
    {
    }

  // FUNCTION: posteriormode_single
  // TASK: estimates only the coefficient of the fix effect 'name'

  virtual void posteriormode_single(const vector<ST::string> & name, datamatrix xnew,
                                    const bool include)
    {
    }

  virtual void safe_splines(bool & interact)
    {
    }

  virtual void set_splines_old(void)
    {
    }

  virtual void safe_const(void)
    {
    }

  virtual void set_const_old(void)
    {
    }

  virtual void posteriormode_const(void)
    {
    }

  virtual const datamatrix & get_data_forfixedeffects(void)
    {
    return data_forfixed;
    }


  virtual void compute_lambdavec(vector<double> & lvec, int & nr);
  void compute_lambdavec_equi(vector<double> & lvec, int & number);

  double lambda_from_df(double & df_wunsch, double & lambda_vorg);

  virtual void update_bootstrap(const bool & uncond=false);

  virtual void update_beta_average(unsigned & samplesize);

  virtual void update_bootstrap_betamean(void);

  virtual void save_betamean(void);

  virtual void update_bootstrap_df(void);

  virtual void outresults_df(unsigned & size)
    {
    }

  void readsample_df(datamatrix & sample,const unsigned & nr,
                          const unsigned & col=0) const;

  virtual void change_Korder(double lam)
    {
    }

  virtual void undo_Korder(void)
    {
    }

  // ---------------------------------------------------------------------------
  // ------------------------------- FOR REML ----------------------------------
  // ---------------------------------------------------------------------------

  double get_startlambda(void)
    {
    return startlambda;
    }

  unsigned get_dimX(void)
    {
    return dimX;
    }

  unsigned get_dimZ(void)
    {
    return dimZ;
    }

  bool get_isnonparametric(void)
    {
    return isnonparametric;
    }

  bool get_catspecific(void)
    {
    return catspecific;
    }

  virtual void createreml(datamatrix & X,datamatrix & Z,
                          const unsigned & Xpos, const unsigned & Zpos)
    {
    }

  virtual void initialize_baseline(unsigned j, datamatrix & tx, datamatrix & tz,
               vector<unsigned> & ts, vector<unsigned> & te,
               vector<unsigned> & tt, datamatrix & iv,
               statmatrix<double> & steps, statmatrix<int> & index)
    {
    }

  virtual double outresultsreml(datamatrix & X,datamatrix & Z,
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
                              const unsigned plotpos)
    {
      return 0.0;
    }

  virtual void outresultsgrid()
    {
    }

    virtual void outoptionsreml()
    {
    }

  // ---------------------------------------------------------------------------
  // --------------------------- end: FOR REML ---------------------------------
  // ---------------------------------------------------------------------------


  };


//------------------------------------------------------------------------------
//----------------- A CLASS FRAGMENT FOR A NEW FULLCOND OBJEKT -----------------
//------------------------------------------------------------------------------


/*
The following variables are already defined In the base class FULLCOND:


Important variables:


MCMCoptions * optionsp;        Pointer to general MCMC options


datamatrix beta;               Matrix of current parameters
datamatrix weight;             Weight matrix for centering beta

unsigned nrpar;                total number of parameters
                               (beta.rows()*beta.cols())

datamatrix data;                Covariates
vector<ST::string> datanames;   Names of the covariates

unsigned long acceptance;       number of accepted iterations
unsigned long nrtrials;         number of trials

bool identifiable;              true, if term is identifiable

unsigned column;                the response category the full conditional
                                belongs to (important for multivariate response)

Not so important:


The following variables are more or less administered by the base class.
There should be no reason to use them in the program code of the inherited
class.

ST::string title;               Title/name of the full conditional

ST::string samplepath;          filename for storing sampled parameters
ofstream samplestream;          stream object for storing sampled parameters

bitset<flagnr> flags;           some flags

datamatrix betamean;            Sampling mean of parameters
datamatrix betas2;              Sampling sum of squares of parameters
datamatrix betavar;             Sampling variance of parameters
datamatrix betamin;             Sampling minimum
datamatrix betamax;             Sampling maximum
datamatrix betaqu10;            Sampling 10 percent quantile
datamatrix betaqu50;            Sampling 50 percent quantile
datamatrix betaqu90;            Sampling 90 percent quantile

datamatrix betameanold;
datamatrix betavarold;
datamatrix betaminold;
datamatrix betamaxold;

bool center;                    true, if beta should be centered



class newfullcond : public FULLCOND
  {

  protected:

  // add here additional variables needed for the new class
  // add here additional functions/methods needed for the inherited class
  // (e.g. a method for centering beta)

  public:

  // DEFAULT CONSTRUCTOR:

  newfullcond(void) : FULLCOND()
    {

    // add program code here

    }

  // CONSTRUCTOR

  newfullcond( // include here variables additionaly needed to initialize the
               // new class
               MCMCoptions * o,const datamatrix & d,
               const ST::string & t, const unsigned & r, const unsigned & c,
               const ST::string & fp)
               : FULLCOND(o,d,t,r,c,fp)
    {

    // initialize here additional variables of the new class
    // initialize the variable weight (if centering of beta will be necessary)
    // initialize variable identifiable (set true if term is ALWAYS identifiable
    //                                   set false, if term may be unidentifiable

    }

  // COPY CONSTRUCTOR

  newfullcond(const newfullcond & fc) : FULLCOND(FULLCOND(fc))
    {
    // assign here ONLY additional variables of the inherited class
    // e.g. y = fc.y
    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const newfullcond & operator=(const newfullcond & fc)
    {
    if (this==&fc)
      return *this;

    // assign here ONLY additional variables of the inherited class
    // e.g. y = fc.y

    return *this;
    }

  // DESTRUCTOR

  ~newfullcond() {}

  // FUNCTION update
  // TASK: updates parameters (i.e. matrix beta)

  void update(void)
    {

    // it is very likely that you will need the following variables from the base
    // class:
    // beta               current parameter vector
    // nrpar              total number of parameters (beta.cols()*beta.rows())

    // data               datamatrix of covariates

    // likep              pointer to DISTRIBUTION object, especially functions
    //                    to modify the linear predictor and to compute the
    //                    likelihood

    // weight             weight matrix for centering (must be specified in the
    //                    constructor or here

    // update here the parameter matrix 'beta'
    //    step 1: compute the likelihood with the old beta
    //            (using likep->loglikelihood)
    //    step 2: compute a new proposed beta
    //    step 3: change the linear predictor according to the new proposed beta
    //            (using likep->add_linearpred), i.e. substract the old part and
    //            add the new part
    //    step 4: compute the likelihood with the new beta
    //            (using likep->loglikelihood)
    //    step 4: if the proposed beta is accepted, set beta =  newbeta
    //            but use the faster command beta.assign(newbeta)
    //            if the proposed beta is recected don't forget to modify
    //            again the linear predictor

    // update the variables 'acceptance' and 'nrtrials', i.e. raise acceptance
    // if your proposed beta has been accepted, raise also nrtrials

    // The term your full conditional belongs to may be unidentifiable
    // the MCMCsimulate object that runs the simulation may have
    // set the boolean variable 'center'. In that case you must center
    // your beta matrix appropriately.
    // DON'T FORGET: centering of beta CHANGES the linear predictor, thus
    // in your centering method the linear predictor must be modified
    // appropriately

    FULLCOND::update();    // this command should be included (at the end of the
                           // function) to update automatically the curent
                           // mean, variance,minimum,maximum etc.
    }

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file (after estimation)

  void outresults(void)
    {
    FULLCOND::outresults(); // this command should be included at the beginning
                            // of the function, because the outresults function
                            // of the base class automatically computes
                            // estimated 10, 50, 90 percent quantiles of the
                            // parameters. In addition the acceptance rate will
                            // be computed and written to (*optionsp->logout)
                            // as well as the Title/Name of the full conditional

    // It is very likely that the following variables from the base class will
    // be needed:

    // betamean             estimated posteriori mean
    // betavar              posteriori variance
    // betaqu10             posteriori 10% quantile
    // betaqu50             posteriori 50% quantile
    // betaqu90             posteriori 90% quantile
    // betamin              posteriori Minimum
    // betamax              posteriori Maximum


    // Write the posterior characteristics of beta either to a
    // file or to (*optionsp->logout)
    // e.g. (*optionsp->logout) << "The mean is: " << betamean(0,0) << endl;

    }


  // FUNCTION: outoptions
  // TASK: writes estimation options to outputstream

  void outoptions(void)
    {

    }


  // FUNCTION: reset
  // TASK: resets all parameters for a new simulation

  void reset(void)
    {
    FULLCOND::reset();   // include this command, because the reset function
                         // of the base class automatically resets all beta
                         // matrices to their starting values (beta = 0)
                         // Sets in addition acceptance =  0;
                         //                  nrtrials = 0;

    // reset here additional variables of the inherited class
    // (e.g. additional scale parameters, auxiliary variables etc.)

    }

  // FUNCTION: predict
  // TASK: used to predict mu for a new observation Xnew
  //       computes the part of the linear predictor belonging
  //       to the full conditional and adds this part to 'linpred'

  void predict(const datamatrix & newX, datamatrix & linpred)
    {

    assert(Xnew.cols() == data.cols());

    // add the part of the linear predictor that belongs to this full
    // conditional

    }


  };

*/


} // end: namespace MCMC

#endif


