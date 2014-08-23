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



#include"../export_type.h"

#if !defined (DAGOBJECT_INCLUDED)

#define DAGOBJECT_INCLUDED

#include"statobj.h"
#include"dataobj.h"
#include"adjacency.h"
#include"fullcond_dag.h"
#include"fullcond_dag_d.h"
#include"fullcond_dag_ia.h"
#include"fullcond_dag_ia_mixed.h"
#include"fullcond_rj.h"
#include"fullcond_rj_int.h"
#include"fullcond_rj_mix.h"
#include"ia.h"
#include"mcmcsimul.h"

using MCMC::MCMCoptions;
using MCMC::DISTRIBUTION;
using MCMC::DISTRIBUTION_gaussian;
using MCMC::FULLCOND;
using MCMC::FULLCOND_dag;
using MCMC::FULLCOND_dag_d;
using MCMC::FULLCOND_dag_ia;
using MCMC::FULLCOND_dag_ia_mixed;
using MCMC::FULLCOND_rj;
//using MCMC::FULLCOND_rj_c;
using MCMC::FULLCOND_rj_int;
using MCMC::FULLCOND_rj_mix;
using MCMC::IA;
using MCMC::MCMCsimulate;

// HINZUFÜGEN EINER NEUEN METHODE

// 1. private - Teil ergänzen
//   - optionlist für die neue Methode definieren
//   - optionen definieren
//   - Modell Klasse definieren
//   - use Klasse definieren
//   - run Funktion definieren
// 2. create - Teil ergänzen (in bayesreg.cpp)
// 3. run funktion schreiben


class __EXPORT_TYPE dagobject : public statobject
  {

  private:

  vector<statobject *> * statobj;

  typedef void (*runpointer) (dagobject & d);

  runpointer functions[10];

  // for method estimate

  optionlist estimateoptions;







  // ************* possible options *******************

  intoption iterations;
  intoption burnin;
  intoption step;
  intoption iterationsprint;
  intoption typ;
  intoption number;

  stroption prior_sig;		// dag: variance of prior
  stroption res_file;		// rj: path for aggregated results
  stroption fix_file;		// rj: path with the conditions on d
  stroption family;			// rj: kind of variable and wether there are interactions or not
  stroption print_models;	// rj: criterion for number of models in the output
  stroption switch_typ;		// rj: kind of switch type


  simpleoption detail_ia;	// dag: detailed output abeout the regression models
  simpleoption print_dags;	// rj: not really necessary....

  doubleoption value_a;		// dag: a of IG(a,b)
  doubleoption value_b;		// dag: b of IG(a,b)
  doubleoption alpha;		// rj: criterion for number of models in the output










  modelStandard mod;

  use udata;

  friend void estimaterun(dagobject & d);

  void create(void);

  public:

  dagobject(void)
    {
    type="dag";
    }

  dagobject(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adb,
  #endif
  const ST::string & n,ofstream * lo,istream * i,ST::string p,
            vector<statobject*> * st);

  dagobject(const dagobject & d);

  const dagobject & operator=(const dagobject & d);

  ~dagobject() {}

  int parse(const ST::string & c);


  };


#endif
