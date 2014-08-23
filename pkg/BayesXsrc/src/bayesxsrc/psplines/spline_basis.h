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



#ifndef spline_basisH
#define spline_basisH

#include"../export_type.h"
#include<deque>
#include "fullcond.h"
#include "mcmc_nonpbasis.h"
#include "bsplinemat.h"

namespace MCMC
{

//---------------------------------------------------------------------------
//----------------------- class: spline_basis -------------------------------
//---------------------------------------------------------------------------


class __EXPORT_TYPE spline_basis : public FULLCOND_nonp_basis
  {

  protected:

  FULLCOND_const * fcconst;             // Zeiger auf Fullcond-Objekt für fixe Effekte + Intercept

  bool pseudocontourprob;               // Sollen Pseudo-Contour-Wahrscheinlichkeiten nach Besag
                                        // berechnet werden?
  bool approx;                          // Sollen die Quantile beim Berechnen der Contour-
                                        // Wahrscheinlichkeiten mit dem Verfahren von
                                        // Tierney approximiert werden?
  int lengthstart;                      // Anzahl der Iteration die verwendet werden sollen
                                        // um Startwerte für die Approximation nach Tierney
                                        // zu berechnen.

  bool lambdaconst;                     // Soll der Glättungsparameter geschätzt werden oder nicht?
  bool outbsplines;                     // Sollen die B-Spline Basisfunktionen ausgegeben werden?
  double lambda_prec;                   // lambda mit dem die aktuelle Präzisionsmatrix berechnet wurde

  bool predictright;                    // Sollen auf der rechten Seite Parameter prognostiziert werden?
  unsigned nrparpredictright;           // Anzahl der Parameter die rechts prognostiziert werden sollen
  bool predictleft;                     // Sollen auf der linken Seite Parameter prognostiziert werden?
  unsigned nrparpredictleft;            // Anzahl der Parameter die links prognostiziert werden sollen

  bool derivative;                      // Soll die erste Ableitung berechnet werden?
  vector<int> index2;                   // Indexvariable zum schnelleren durchlaufen
                                        // von spline, linpred, etc.

  bool increasing;                      // für monotone Regression: monoton steigender Effekt
  bool decreasing;                      // für monotone Regression: monoton fallender Effekt

  datamatrix interactvar;               // speichert bei variierden Koeffizienten
                                        // die Werte der Interaktionsvariable

  datamatrix W;                         // Gewichtsmatrix
  datamatrix betaold;                   // altes beta (aus der letzten Iteration)
  datamatrix betaprop;                  // Vorgeschlagenes beta

  bandmatdouble XX;                     // X'X
  bandmatdouble prec;                   // Präzisionsmatrix (X'WX + 1/sigma2*K)

  envmatdouble prec_env;                // Präzisionsmatrix als Envelope-Matrix (X'WX + 1/sigma2*K)
  envmatdouble XX_env;                  // X'X als Envelope-Matrix

  datamatrix mu;                        // tildey
  datamatrix muy;                       // X'W*tildey
  datamatrix standnormal;               // N(0,1)-verteilte ZV
  datamatrix betahelp;                  // Hilfsvariable


  FULLCOND fchelp;                      // Fullcond-Objekt zum Abspeichern der Samples
                                        // des Splines in jeder Iteration
  FULLCOND fcderivative;                // Fullcond-Objekt zum Abspeichern der Samples
                                        // für die erste Ableitung in jeder Iteration

  bsplinemat Bderivative;               // Designmatrix für die erste Ableitung
                                        // (f'(x) = Bderivative*beta)
  datamatrix splinederivative;          // Erste Ableitung des Splines

  datamatrix Bcolmean;                  // Spaltenmittelwert der B-Spline-Designmatrix

  unsigned nrknots;                     // Anzahl der (sichtbaren) Knoten
  unsigned degree;                      // Grad des Splines
  unsigned nrdiffobs;                   // Anzahl verschiedener Kovariablen-Werte

  int gridsize;                         // -1 falls Schätzungen für alle verschiedenen
                                        // Kovariablen-Werte ausgegeben werden sollen
                                        // sonst: Anzahl von äquidistanten Punkten, an
                                        // denen Schätzungen ausgegeben werden sollen
  double intercept;                     // Zentrierungskonstante
  knotpos knpos;                        // 'equidistant' oder 'quantiles'

  vector<int> freq;                     // Wert wird um 1 erhöht, wenn ein Kovariablen-Wert
                                        // der sortierten Beobachtungen verschieden vom
                                        // vorherigen ist. Wichtig für die Berechnung von X*beta,
                                        // X'WX etc., da in X nur Einträge für VERSCHIEDENE
                                        // Kovariablen-Werte gespeichert werden.
                                        // Beobachtungen ausgegeben werden
  vector<int> freqoutput;               // Für die Ausgabe. Wichtig damit nur verschiedene
                                        // Beobachtungen ausgegeben werden. Identisch mit
                                        // freq, falls varcoeff==false

  deque<int> firstnonzero;              // Gibt pro Spalte die Position des ersten Elements
                                        // der B-Spline-Designmatrix das ungleich Null ist
  deque<int> lastnonzero;               // Gibt pro Spalte die Position des letzten Elements
                                        // der B-Spline-Designmatrix das ungleich Null ist
  deque<double> knot;                   // Vektor der Knoten (sichtbare und unsichtbare)

  datamatrix xvalues;                   // geordente, verschiedene Kovariablen-
                                        // Werte (für die Ausgabe)
  datamatrix spline;                    // aktueller Spline
  datamatrix splinehelp;                // Hilfmatrix
  datamatrix betaweight;                // Gewichte für compute_intercept

  datamatrix B;                         // 0 für nonp, X für varcoeff==true
                                        // (enthält nur Einträge ungleich 0 und nur
                                        // für verschiedene Kovariablen-Werte.)
  datamatrix BS;                        // X für nonp, XZ für varcoeff==true
                                        // (enthält nur Einträge ungleich 0 und nur
                                        // für verschiedene Kovariablen-Werte.)
  datamatrix G;                         // Transformationsmatrix für diagtransform

  datamatrix X_VCM;                     // für REML VCM
  datamatrix Z_VCM;                     // für REML VCM

  datamatrix X_grid;                    // Evaulation on a predefined grid
  datamatrix Z_grid;                    //

  datamatrix xvaluesgrid;

  double lowergrid;                     // Variable for SIMEX
  double uppergrid;
  double lowerknot;
  double upperknot;

  double reference;                     // Relative effect with respect to a reference value
  bool refcheck;
  datamatrix X_ref;
  datamatrix Z_ref;

  vector<int> begcol;                   // Gibt pro Zeile die Position des ersten Elements
                                        // der B-Spline-Designmatrix an das ungleich Null ist

  datamatrix DG;                        // Hilfs-Designmatrix, falls gridsize>0
  vector<int> DGfirst;                  // Gibt pro Zeile die Position des ersten Elements
                                        // von DG an das ungleich Null ist


  // FUNCTION: make_index
  // TASK: berechnet index, freq und freqoutput

  void make_index(const datamatrix & moddata);

  // FUNCTION: make_index
  // TASK: berechnet index, freq und freqoutput bei VCM-Modellen

  void make_index(const datamatrix & em,const datamatrix & ia);

  // FUNCTION: make_index2
  // TASK: computes index2

  void make_index2(void);

  // FUNCTION: make_Bspline
  // TASK: berechnet B, BS, knot, nrparpredictleft, nrparpredictright, Bcolmean
  // md:        Datenmatrix
  // minnull:   false - erster sichtbarer Knoten ist bei x_min
  //            true  - erster sichtbarer Knoten ist bei 0 (für baseline bei Cox)

  void make_Bspline(const datamatrix & md,const bool & minnull = false);

  // FUNCTION: make_BS
  // TASK: berechnet die B-Spline-Designmatrix für VCM-Modelle (B=BS, BS = diag(ia)*BS)
  // NOTE: make_Bspline muss bereits ausgeführt sein
  // ia:   Interaktionsvariable

  void make_BS(const datamatrix & ia);

  // FUNCTION: make_DG
  // TASK: berechnet die B-Spline-Designmatrix falls gridsize>0

  void make_DG(void);

  // FUNCTION: make_DG_REML
  // TASK: berechnet die B-Spline-Designmatrix falls gridsize>0

  void make_DG_REML(void);

  // FUNCTION: init_fchelp
  // TASK: initialisiert xvalues, fchelp und fcderivative

  void init_fchelp(const datamatrix & d);

  // FUNCTION: compute_betaweight
  // TASK: berechnet die nötigen Gewichte zur Berechnung der Zentrierungskonstanten

  void compute_betaweight(void);

  // FUNCTION: compute_intercept
  // TASK: berechnet die Zentrierungskonstante: intercept = betaweight*beta

  void compute_intercept(void);

  // FUNCTION: compute_intercept
  // TASK: berechnet die Zentrierungskonstante mit der übergebenen Matrix

  void compute_intercept(const datamatrix & beta);

  // FUNCTION: subtr_spline
  // TASK: berechnet eta = eta - spline + intercept

  void subtr_spline(void);

  // FUNCTION add_linearpred_multBS
  // TASK: multiplies the X-matrix with beta and adds the result to
  //       the current (proposed) linear predictor
  //       The result is also assigned to 'spline'

  void add_linearpred_multBS(const bool & current = true);

  // FUNCTION add_linearpred_multBS
  // TASK: multiplies the X-matrix with its argument and adds the result to
  //       the current (proposed) linear predictor
  //       The result is also assigned to 'spline'

  void add_linearpred_multBS(const datamatrix & beta,const bool & current = true);

  // FUNCTION add_linearpred_multBS
  // TASK: multiplies the X-matrix with (beta1-beta2) and adds the result to
  //       the current (proposed) linear predictor and computes spline = Xbeta1

  void add_linearpred_multBS(const datamatrix & beta1,const datamatrix & beta2, const bool & current = true);

  // FUNCTION add_linearpred_multBS
  // TASK: multiplies the X-matrix with beta[a,e] and adds the result to
  // NOTE: for conditional prior proposals

  void add_linearpred_multBS_Block(const unsigned a, const unsigned e, const datamatrix & b);

  // für Cox-Modell, wenn baseline direkt modelliert wird anstatt log(baseline)  (condprior)

//  void add_linearpred_multBS_Block2(const unsigned a, const unsigned e, const datamatrix & b);

  // FUNCTION compute_XWX
  // TASK: computes X'diag(weight)X and assigns the result to XX

  void compute_XWX(const datamatrix & weight);     // nur für diagtransform

  // FUNCTION compute_XWXenv
  // TASK: computes X'diag(weight)X and assigns the result to XX_env

  void compute_XWXenv(const datamatrix & weight, const unsigned & c=0);

  // FUNCTION compute_XWtildey
  // TASK: computes scale * X'diag(weight)muy and assigns the result to muy

  void compute_XWtildey(const datamatrix & weight, const double & scale);

  // FUNCTION compute_XWXenv_XWtildey
  // TASK: computes X'diag(weight)X and assigns the result to XX_env
  //       computes scale * X'diag(weight)muy and assigns the result to muy

  void compute_XWXenv_XWtildey(const datamatrix & weight, const double & scale, const unsigned & c=0);

  // FUNCTION compute_XWtildey
  // TASK: computes scale * X'diag(weight)tildey and assigns the result to muy
  // NOTE: für posteriormode

  void compute_XWtildey(const datamatrix & weight, const datamatrix & tildey, const double & scale, const unsigned & c=0);

  // FUNCTION: sample_centered
  // TASK: Sample under condition x|Ax=0
  //       x - Q^-1AT(AQ^-1AT)^-1(Ax), V=Q^-1AT
  // NOTE: assumes that the precision matrix is stored in prec

  void sample_centered(datamatrix & beta);

  // FUNCTION: sample_centered_env
  // TASK: Sample under condition x|Ax=0
  //       x - Q^-1AT(AQ^-1AT)
  // NOTE: assumes that the precision matrix is stored in prec_env

  void sample_centered_env(datamatrix & beta);

  // FUNCTION: compute_Kweights
  // TASK: computes the weights 'weight' to be used when computing the penalty matrix K
  //       K = D'diag(weight)D

  void compute_Kweights(void);

  // FUNCTION: bspline_rek
  // TASK: Berechnet B-Spline Basisfunktion and der Stelle X
  // NOTE: needed in function 'predict' to compute B-Splines for a single observation

  double bspline_rek(unsigned l, unsigned knot, const datamatrix & X);

  // Functions to write samples to full conditional objects

  void write_spline(void);

  void write_spline(const datamatrix & beta);

  void write_derivative(void);

  // Function to write B-Spline basis functions to a file

  void write_bsplinefunctions(const datamatrix & beta,datamatrix & bsplines);

  // Functions for prediction

  // FUNCTION: change_K
  // TASK: changes the entries of K and Kenv
  // NOTE: only necessary, when predictleft==true || predictright==true

  void change_K(void);

  // FUNCTION: update_prediction
  // TASK: updates the parameters to be predicted

  void update_prediction(void);

  // FUNCTION: sample_monotonic
  // TASK: samples beta under monotonicity constraints

  double sample_monotonic(const unsigned i, const double m, const double s);

  bool is_monotonic(const unsigned i);
  bool is_monotonic(void);


  public:


  // DEFAULT CONSTRUCTOR

  spline_basis(void) : FULLCOND_nonp_basis()
    {
    }

  // CONSTRUCTOR

  spline_basis(MCMCoptions * o, DISTRIBUTION * dp,
               FULLCOND_const * fcc, const fieldtype & ft,
                const ST::string & ti, const unsigned & nrk, const unsigned & degr,
                const MCMC::knotpos & kp, const int & gs, const ST::string & fp,
                const ST::string & pres, const bool & deriv, const double & lk,
                const double & uk, const double & lg, const double & ug,
                const unsigned & c);

  // CONSTRUCTOR für REML

  spline_basis(MCMCoptions * o, const datamatrix & d, const unsigned & nrk, const unsigned & degr,
               const knotpos & kp, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp, const double & lg,
               const double & ug, const double & lk, const double & uk,
               const int & gs, const double & rv);

  // CONSTRUCTOR für REML VCM

  spline_basis(MCMCoptions * o, const datamatrix & d1, const datamatrix & d2,
               const unsigned & nrk, const unsigned & degr,
               const knotpos & kp, const fieldtype & ft, const ST::string & ti,
               const ST::string & fp, const ST::string & pres, const double & l,
               const double & sl, const bool & catsp, const bool & ctr, const double & rv);

  // COPY CONSTRUCTOR

  spline_basis(const spline_basis & sp);

  // OVERLOADED ASSIGNMENT OPERATOR

  const spline_basis & operator=(const spline_basis & sp);

  // Funktionen zum ändern der Haupteffekte bei Interaktionen

  // main: Matrix mit den Haupteffekten
  //       Im Gauss-Fall:       'rausgerechneter' Haupteffekt
  //       Im Nicht-Gauss-Fall: 'rausgerechnetes' beta
  // inter: Zentrierungskonstante

  // für die MCMC-Simulation

  void change(const datamatrix & main, const double & inter);            // Gauss-Fall
  void change(const datamatrix & main);                                  // Nicht-Gauss-Fall

  // für die Posteriori-Modus-Schätzung

  bool changeposterior(const datamatrix & main, const double & inter);   // Gauss-Fall
  bool changeposterior(const datamatrix & main);                         // Nicht-Gauss-Fall

  // FUNCTION: bspline
  // TASK: berechnet Werte der B-Spline Basisfunktionen an der Stelle x

  datamatrix bspline(const double & x);

  // FUNCTION: bspline
  // TASK: berechnet Werte der B-Spline Basisfunktionen vom Grad d an der Stelle x

  datamatrix bspline(const double & x, const unsigned & d);

  // FUNCTION: deriv_f
  // TASK: berechnet den Spline an der Stelle x mit dem aktuellen beta

  double deriv_f(const double & x);

  // FUNCTION: multBS
  // TASK: computes BS*beta (B*beta respectively) and stores the result in res
  // NOTE: the result is sorted according to the covariate

  void multBS(datamatrix & res, const datamatrix & beta);

  // FUNCTION: multBS_index
  // TASK: computes BS*beta (B*beta respectively) and stores the result in res
  // NOTE: the result is sorted according to the linear predictor

  void multBS_index(datamatrix & res, const datamatrix & beta);

  // FUNCTION: multDG
  // TASK: computes DG*b and stores the result in res

  void multDG(datamatrix & res, const datamatrix & b);

  void outoptions(void);

  void outresults(void);

  unsigned & get_nrknots(void)
    {
    return nrknots;
    }

  deque<double> & get_knots(void)
    {
    return knot;
    }

  datamatrix & get_spline(void)
    {
    return spline;
    }

  datamatrix & get_splinehelp(void)
    {
    return splinehelp;
    }

  int * get_indexp(void)
    {
    return index.getV();
    }

  vector<int>::iterator get_freqit(void)
    {
    return freq.begin();
    }

  double get_intercept(void)
    {
    return intercept;
    }

  double * get_fchelpbetap(void)
    {
    return fchelp.getbetapointer();
    }

  void fchelpupdate(void)
    {
    fchelp.update();
    }

  bandmatdouble & get_XX(void)
    {
    return XX;
    }

  int & get_gridsize(void)
    {
    return gridsize;
    }

  unsigned & get_degree(void)
    {
    return degree;
    }

  knotpos & get_knotpos(void)
    {
    return knpos;
    }

  // FUNCTION: getX
  // TASK: schreibt BS (B bei VCM) in datamatrix X
  // NOTE: X muss die richtige Dimension haben

  void getX(datamatrix & X);

  void init_name(const ST::string & na);

  void init_names(const vector<ST::string> & na);

  void set_lambdaconst(double la);

  void set_contour(int cp, bool pseudocp, bool app, int ls,
                    const datamatrix & b = datamatrix(1,1,0.0));

  void set_outbsplines(void)
    {
    outbsplines = true;
    }

  double compute_df(void);

  double compute_df_eigen(void);

  void get_effectmatrix(datamatrix & e,vector<ST::string> & enames,
                        unsigned be, unsigned en,effecttype t);


  unsigned get_nreffects(effecttype t);

  // FUNCTION: init_data_varcoeff
  // TASK: initializes data and data2 (data^2) for varying coefficient model

  void init_data_varcoeff(const datamatrix & intvar, double add=0);

  // ------------------------- FOR MERROR --------------------------------------

  void update_merror_varcoef(datamatrix & effmod, datamatrix & newintact);

  void update_merror(datamatrix & newdata);
  void update_merror_discrete(datamatrix & newdata);
  datamatrix get_spline_merror(void);

//  void make_index_discrete(const datamatrix & moddata, const datamatrix & grid);

//  datamatrix discretise(datamatrix & moddata);

//  void init_fchelp(const datamatrix & d, datamatrix & grid);

  // -------------------------END: FOR MERROR ----------------------------------

  // ------------------------- FOR REML ----------------------------------------

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

  ~spline_basis(){}

  };


} // end: namespace MCMC

//---------------------------------------------------------------------------
#endif
