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



#if !defined (STEPWISErun_INCLUDED)

#define STEPWISErun_INCLUDED

#include"../export_type.h"
#include"mcmcsimul.h"


namespace MCMC
{

class __EXPORT_TYPE STEPWISErun : public MCMCsimulate
  {

  protected:

// -----------------------------------------------------------------------------
// ------------ Variablen, die in der ganzen Datei bekannt sind: ---------------
// -----------------------------------------------------------------------------

  datamatrix D;
  vector<ST::string> modelv;
  vector<FULLCOND*> fullcond_alle;     // Fullcond-Vektor, wie er zu Beginn übergeben wird
  ST::string algorithm;                // Minimierungsalgorithmus (stepwise, Koordinatenabstieg)
  ST::string minim;                    // Art der univariaten Minimierung bei Koordinatenabstieg (z.B. adaptiv, exact)
  ST::string minim2;
  ST::string criterion;
  int increment;                       // Anzahl der zu testenden Nachbar-Alternativen bei stepwise
  int steps;
  ST::string startmodel;
  ST::string trace;                    // Ausführlichkeit der Ausgabe am Bildschirm
  double kriterium_tex;
  ofstream outmodels;
  ofstream outcriterium;
  ofstream outtex;
  bool hierarchical;
  int bootstrap;
  bool isboot;                         // ist der Algorithmus gerade beim Bootstrap-Teil?
  bool unconditional;

  vector<vector<double> > lambdavec;  // enthält für jedes fullcond-Objekt (außer fixen Effekten) alle Modellierungs-Alternativen
  vector<ST::string> names_fixed;     // Namen der fixen Effekte (einschließlich "const")
  vector<vector<ST::string> > names_nonp; // enthält für jedes fullcond-Objekt (außer fixen Effekten) alle Modellierungs-Alternativen
  vector<double> modell_neu;        // enthält für jede Variable/Funktion die Modellierungsalternative
  vector<double> modell_alt;
  double kriterium_alt;
  double kriterium_neu;
  ST::string text_alt;             // für Bildschirm-Ausgabe
  vector<vector<vector<double> > > modellematrix;  // speichert bereits ausprobierte Modelle ab
  bool fertig;
  int steps_aktuell;               // Laufindex für Iteration
  vector<ST::string> posttitle;

  void schaetzen(int z, double & kriterium, bool neu, ST::string variante);

// -----------------------------------------------------------------------------
// -------------- Funktionen, für Stepwise / Stepmin ---------------------------
// -----------------------------------------------------------------------------

  bool stepfunctions(void);

  unsigned stepwise_fixfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration);

  void stepwise_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z);

  void stepmin_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z);

  void minexact_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z);

// -----------------------------------------------------------------------------
// ------------------ Funktionen für Stepmin -----------------------------------
// -----------------------------------------------------------------------------

  void step_minfix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration);

  void stepmin_fix_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration, unsigned & i);

  void stepmin_leer_fix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration, unsigned & i);

  unsigned step_minfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration);

  void stepmin_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z);

  void stepmin_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z);

  void stepmin_nonp_nonp(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  void stepmin_nonp_fix(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  void stepmin_nonp_leer(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  void minexact_nonp_nonp(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  void minexact_nonp_fix(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  void minexact_nonp_leer(unsigned & z, vector<double> & krit_fkt, double & kriterium);

  double criterion_min(const double & df);

  double criterion_min(const double & df, const ST::string & auswahl);

// -----------------------------------------------------------------------------
// ------------------ Funktionen für Koordinatenmethode ------------------------
// -----------------------------------------------------------------------------

  bool koordabstieg(void);

  // Funktion für Modellierungsänderung bei fixen Effekten
  void koord_minfix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell);

  // fixe Effekte; alter Zustand: Variable drin, auszuprobieren: Variable raus
  void koord_fix_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & i);

  // fixe Effekte; alter Zustand: Variable raus, auszuprobieren: Variable rein
  void koord_leer_fix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & i);

  // Funktion für Modellierungsänderung bei Faktor-Variablen
  unsigned koord_minfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell);

  // Faktor-Variable; alter Zustand: Faktoren drin, auszuprobieren: Faktoren raus
  void koord_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z);

  // Faktor-Variable; alter Zustand: Faktoren raus, auszuprobieren: Faktoren rein
  void koord_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z);

  void koord_minnonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z, double & kriterium_aktuell);

  unsigned koordexact_fixfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell);

  void koordexact_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z, double & kriterium_aktuell);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des Startmodels -----------------------
// -----------------------------------------------------------------------------

  // Fehler-Abfrage vor Start des Algorithmus
  bool vcm_doppelt(void);

  // Hier werden die Modellierungsalternativen bestimmt
  void initialise_lambdas(vector<vector<ST::string> > & namen_nonp,
       vector<ST::string> & namen_fix, vector<vector<double> > & lambdavector,
       const int & number, const bool & gewichte);

  // Bestimmung der 0/1 Gewichte bei MSEP
  void initialise_weights(double prop);

  // sucht den Index von "lambdavec" bei gegebenem Lambda
  unsigned search_lambdaindex(const double & m, const vector<double> lam,
                                            bool & b) const;

  // sucht ein dem vorgegebenen Startwert ähnliches Lambda raus
  unsigned search_lambdastartindex(const double & start,
                           const vector<double> & lambdas) const;

  // bestimmt Lambdas für's Startmodell
  void startwerte(const ST::string & startmodel,
       vector<vector<unsigned> > & startindex, vector<vector<double> > & startfix);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Berechnung neuer Modelle bei Stepwise ------------
// -----------------------------------------------------------------------------

  void newmodel(vector<double> & krit, vector<vector<double> > & mi,
    vector<ST::string> & textit);

  void newmodel_fix(const double & mo, vector<double> & krit,
      vector<vector<double> > & mi, vector<ST::string> & textit,
      const ST::string & name);

  void newmodel_factor(const double & mo,  const unsigned & index,
      vector<double> & krit, vector<vector<double> > & mi,
      vector<ST::string> & textit, const vector<ST::string> & name);

  void newmodel_nonp(const unsigned & index, vector<double> & krit,
     vector<vector<double> > & mi, vector<ST::string> & textit);

  bool modelcomparison(const vector<double> & m,
       const vector<vector<vector<double> > > & mmatrix);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des fullcondp-Vektors -----------------
// -----------------------------------------------------------------------------

  // ändert Fullcond-Vektor für einen nichtlinearen Effekt "i" von modell2 --> modell1
  void fullcond_einzeln(const vector<double> & modell1,
         const vector<double> & modell2, const unsigned & index);

  // stellt den Fullcond-Vektor neu auf, passend zu Modell-Vektor "m" (außer für fixe Effekte)
  void fullcond_komplett(const vector<double> & m);

  // stellt die fixen Effekte richtig zusammen, passend zu Modell-Vektor "m";
  // die Reihenfolge der fixen Effekte wird hier verändert
  void fix_komplett(const vector<double> & modell);

  // stellt die fixen Effekte richtig zusammen, passend zu Modell-Vektor "m";
  // die Reihenfolge der fixen Effekte ist hier immer wie im Regressions-Befehl
  void fix_ganz_komplett(const vector<double> &  modell);

  // entfernt einen bestimmten fixen Effekt aus fullcond-Objekt der fixen Effekte
  void reset_fix(const ST::string & name);

  // fügt einen bestimmten fixen Effekt dem fullcond-Objekt der fixen Effekte hinzu
  void include_fix(const ST::string & name);

  // sucht die Spalte der Datenmatrix für einen bestimmten fixen Effekt heraus
  int column_for_fix(const ST::string & name);

  // passt den Intercept nach Entfernen / hinzufügen von fixen Effekten an
  // so an, dass Mittelwert(Prädiktor) = Mittelwert((Arbeits-)Beobachtungen)
  void korrektur(void);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Ausgabe im Output-Fenster ------------------------
// -----------------------------------------------------------------------------

  bool make_pause(void);

  void maketext(const ST::string & h, const vector<double> & m,
                const double & a, ST::string & text, const bool & neutext,
                const ST::string & tr,const bool & datei);

  void options_text(const int & number, const vector<vector<double> > & startfix,
      const vector<vector<unsigned> > & startindex, const ST::string & name);

// -----------------------------------------------------------------------------
// ------- Funktionen für die Ausgabe im Tex-File ------------------------------
// -----------------------------------------------------------------------------

  void make_graphics(const ST::string & name, vector<vector<unsigned> > & startindex);

  void make_tex_end(ST::string & path, vector<double> & modell,const ST::string & CI);

  void make_options(void);

  void make_predictor(void);

  void make_model(void);

  void make_prior(vector<vector<unsigned> > & startindex);

  void make_fixed_table(const ST::string & CI);

  void make_plots(ST::string & path_batch,ST::string & path_splus);
                 // ST::string & path_stata);


// -----------------------------------------------------------------------------
// ------------- Model Averaging -----------------------------------------------
// -----------------------------------------------------------------------------

  void update_bootstrap(unsigned & zaehler);

  bool confidence_intervals(const ST::string & CI,
          const vector<double> & modell_final,const double & kriterium_final,
          vector<FULLCOND*> & fullcond_z);

  bool confidence_bootstrap(const vector<double> & modell_final,const double & kriterium_final,
                                          vector<FULLCOND*> & fullcond_z);

  bool confidence_MCMCbootstrap(const vector<double> & modell_final,const double & kriterium_final,
                                          vector<FULLCOND*> & fullcond_z);

  bool confidence_MCMCselect(const vector<double> & modell_final,const double & kriterium_final,
                                          vector<FULLCOND*> & fullcond_z);

  // für Ziehen von Zufallszahlen (ähnlich wie Funktion in "mcmcsimul", aber es muß unterschieden werden,
  // ob Bootstrap-Algorithmus oder nicht)
  bool simulate(const vector<ST::string> & header, const int & seed,
                           const unsigned & startit, const unsigned & endit);



  public:

  // DEFAULT CONSTRUCTOR

  STEPWISErun(void)
    {
    }

  // CONSTRUCTOR1
  // TASK: initializes the MCMC simulation object with general MCMC options 'go'
  //       distribuiton object dp and a vector of full conditionals 'fc'

  STEPWISErun(MCMCoptions * go,DISTRIBUTION * dp,vector<FULLCOND*> & fc);


  // COPY CONSTRUCTOR

  STEPWISErun(const STEPWISErun & s);

  // OVERLOADED ASSIGNMENT CONSTRUCTOR

  const STEPWISErun & operator=(const STEPWISErun & s);

  // ähnlich wie in "mcmcsimul.cpp", aber mit Anpassung an Bootstrap-Algorithmus
  bool posteriormode(const vector<ST::string> & header, const bool & presim);

  bool single_stepwise(const vector<unsigned> & start,
                         const vector<double> & startfix, const bool & tex);

  bool stepwise(const ST::string & procedure, const ST::string & minimum,
         const ST::string & crit, const int & stp, const ST::string & trac,
         const int & number, const ST::string & stam, const int & inc,
         const int & boot, const bool & uncond,
         const datamatrix & D,const vector<ST::string> & modelv,
         const ST::string & name, vector<FULLCOND*> & fullcond_z, ST::string & path,
         const ST::string & CI, bool & hier, const double & prop);

  double compute_criterion(void);

  // DESTRUCTOR

  ~STEPWISErun() {}

  };


} // end: namespace MCMC

#endif










