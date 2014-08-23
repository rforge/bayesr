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





#include"model.h"
#include<algorithm>

//------------------------------------------------------------------------------
//-------------- CLASS model: implementation of member functions ---------------
//------------------------------------------------------------------------------


model::model(const model & m)
  {                                                  
  modelexisting = m.modelexisting;
  modeltext = m.modeltext;
  errormessages = m.errormessages;
  modelVarnames = m.modelVarnames;
  }


const model & model::operator=(const model & m)
  {
  if (this == &m)
	 return *this;
  modelexisting = m.modelexisting;
  modeltext = m.modeltext;
  errormessages = m.errormessages;
  modelVarnames = m.modelVarnames;
  return *this;
  }


vector<ST::string> model::getModelVarnamesAsVector()
  {
  vector<ST::string> help;
  if (! modelVarnames.empty())
    {
    list<ST::string>::iterator it;
    for (it = modelVarnames.begin();it != modelVarnames.end();++it)
      help.push_back(*it);
    }
  return help;
  }

void model::makeModelMatrix_j(dataset & ds,datamatrix & d,const unsigned & j)
  {
  if (modelexisting == true)
    {
//    int j2 = j;
//    list<ST::string>::iterator it = modelVarnames.begin()+j2;
//    ST::string var_j = *it;
    vector<ST::string> varn = getModelVarnamesAsVector();
    ST::string var_j = varn[j];
    ds.makematrix(var_j,d);
    errormessages = ds.geterrormessages();
    if (! errormessages.empty())
      clear();
    }
  }



//------------------------------------------------------------------------------
//---------- CLASS modelStandard: implementation of member functions -----------
//------------------------------------------------------------------------------


const modelStandard & modelStandard::operator=(const modelStandard & m)
  {
  if (this == &m)
	 return *this;
  model::operator=(model(m));
  return *this;
  }


void modelStandard::parse (const ST::string & m)

  {
  model::parse(m);
  modelVarnames = m.strtokenlist(" ");
  modelexisting =  true;
  modeltext = m;
  }


//------------------------------------------------------------------------------
//------------ CLASS expression: implementation of member functions ------------
//------------------------------------------------------------------------------


expression::expression(const expression & e) : model(model(e))
  {
  varname = e.varname;
  expr = e.expr;
  }


const expression & expression::operator=(const expression & e)
  {
  if (this == &e)
	 return *this;
  model::operator=(model(e));
  varname = e.varname;
  expr = e.expr;
  return *this;
  }


void expression::parse(const ST::string & e)
  {

  model::parse(e);

  int equalsignpos = e.checksign('=');
  if (equalsignpos == -1)
	 errormessages.push_back("ERROR: \"=\" expected\n");
  else if (e.length() <= equalsignpos+1)
	 errormessages.push_back("ERROR: expression expected\n");
  else
	 {
     if (equalsignpos > 0)
       {
	   varname = e.substr(0,equalsignpos);
	   varname = varname.eatwhitespace();
	   if (varname.isvarname() == 1)
		 errormessages.push_back("ERROR: " + varname + " invalid varname\n");
       }
     else
       {
       errormessages.push_back("ERROR: new varname expected\n");
       }

     if (e.length()-equalsignpos-1>0)
       {
	   expr = e.substr(equalsignpos+1,e.length()-equalsignpos-1);
	   expr = expr.eatwhitespace();
       }
     else
       errormessages.push_back("ERROR: expression expected\n");  

	 }

  if (errormessages.empty())
	 {
	 modelexisting = true;
	 modeltext = e;
	 }
  else
	 clear();

  }


//------------------------------------------------------------------------------
//--------------- CLASS term: implementation of member functions ---------------
//------------------------------------------------------------------------------


void term::parse(const ST::string & t)
  {

  clear();

  ST::string te;
  te = t.eatallwhitespace();
  te = te.eatallcarriagereturns();
  if (te.length() == 0)
    {
    errormessages.push_back("ERROR: invalid term specification");
    }
  else
    {
    ST::string functionname;
    ST::string argument;
    int isfunc = te.isfunction(functionname,argument);
    if (isfunc == 0)
      {
      vector<ST::string> token = te.strtoken(" *",false);
      unsigned i;
      for(i=0;i<token.size();i++)
        {
        if (token[i].isvarname() == 0)
          varnames.push_back(token[i]);
        else
          errormessages.push_back("ERROR: " + token[i] +
                                  " is not a valid varname\n");
        }
      }
    else if (isfunc==-1)
      errormessages.push_back("ERROR: missing bracket(s) in " + te + "\n");
    else
      {
      vector<ST::string> token = functionname.strtoken(" *",false);
      unsigned i;
      for(i=0;i<token.size();i++)
        {
        if (token[i].isvarname() == 0)
          varnames.push_back(token[i]);
        else
          errormessages.push_back("ERROR: " + token[i] +
                                  " is not a valid varname\n");
        }

      if (argument.length() > 0)
        options = argument.strtoken(",",false);

      }
    }
  }


//------------------------------------------------------------------------------
//---------- CLASS basic_termtype: implementation of member functions ----------
//------------------------------------------------------------------------------


vector<ST::string> basic_termtype::get_constvariables(vector<term> & terms)
  {
  vector<ST::string> res;
  unsigned i;
  for(i=0;i<terms.size();i++)
    {
    if (terms[i].type == "basic_termtype")
      res.push_back(terms[i].varnames[0]);

    }
  return res;
  }

//------------------------------------------------------------------------------
//---------- class term_shrinkage: implementation of member functions ----------
//------------------------------------------------------------------------------

// DEFAULT CONSTRUCTOR
//---------------------
term_shrinkage::term_shrinkage(void)
  {
  type = "term_shrinkage";

  // Einlesen der Starwerte aus externer Datei
  startdata = stroption("startdata");

  // Startwerte für die shhrinkgeeffekte
  effectstart = doubleoption("effect",1E8,-1E7,1E7) ;
  
  // Startwert fuer inverse varianzparameter tau^2
  tau2 = doubleoption("tau2",0.1,0,10000000);

  // Startwert für den Shrinkageparameter
  shrinkagestart = doubleoption("shrinkage",1,0,10000000);

  // Hyperparameter der Priori fuer Shrinkageparameter
  a_shrinkage = doubleoption("a",0.001,0,500);
  b_shrinkage = doubleoption("b",0.001,0,500);

  // Gewichte für die Shrinkageeffekte
  weight = doubleoption("weight",1,0,10000000);

  // Feste Werte für den Shrinkageparameter
  shrinkagefix = simpleoption("shrinkagefix",false);

  // Varianzspezifische Werte für den Shrinkageparameter
  shrinkageadaptive = simpleoption("adaptive",false);

  }


// FUNCTION: setdefault
//---------------------
void term_shrinkage::setdefault(void)
  {
  // call setdefault-methods of the options
  startdata.setdefault();
  effectstart.setdefault();
  tau2.setdefault();
  shrinkagestart.setdefault();
  weight.setdefault();
  a_shrinkage.setdefault();
  b_shrinkage.setdefault();
  shrinkagefix.setdefault();
  shrinkageadaptive.setdefault();

  }


// FUNCTION: check
//-----------------
bool term_shrinkage::check(term & t)
  {
  if ( (t.varnames.size() == 1) && (t.options.size()<=10) && (t.options.size()>=1))   // SET: Anzahl Optionen
    {
    // extract options
    // prior coorespponds to L2-Norm df coefficients
    if (t.options[0] == "ridge")
      t.type = "ridge";

    // prior coorespponds to L1-Norm df coefficients
    else if (t.options[0] == "lasso")
      t.type = "lasso";
    else
      {
      setdefault();
      return false;
      }

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&startdata);
    optlist.push_back(&effectstart);
    optlist.push_back(&tau2);
    optlist.push_back(&shrinkagestart);
    optlist.push_back(&weight);
    optlist.push_back(&a_shrinkage);
    optlist.push_back(&b_shrinkage);
    optlist.push_back(&shrinkagefix);
    optlist.push_back(&shrinkageadaptive);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(10);                        // SET: Anzahl Optionen
    t.options[0] = t.type;
    t.options[1] = startdata.getvalue();
    t.options[2] = ST::doubletostring(effectstart.getvalue());
    t.options[3] = ST::doubletostring(tau2.getvalue());
    t.options[4] = ST::doubletostring(shrinkagestart.getvalue());
    t.options[5] = ST::doubletostring(weight.getvalue());
    t.options[6] = ST::doubletostring(a_shrinkage.getvalue());
    t.options[7] = ST::doubletostring(b_shrinkage.getvalue());
    if (shrinkagefix.getvalue()==false)
       t.options[8] = "false";
    else
       t.options[8] = "true";
    if (shrinkageadaptive.getvalue()==false)
       t.options[9] = "false";
    else
       t.options[9] = "true";
       
    setdefault();
    return true;
    }

  else
    {
    setdefault();
    return false;
    }
  }


// FUNCTION: checkvector
//----------------------
bool term_shrinkage::checkvector(const vector<term> & terms, const unsigned & i)
  {
  assert(i< terms.size());
  if (terms[i].type == "ridge" || terms[i].type == "lasso")
    return true;
  return false;
  }


//------------------------------------------------------------------------------
//---------- class term_nigmix: implementation of member functions -------------
//------------------------------------------------------------------------------

// DEFAULT CONSTRUCTOR
//---------------------
term_nigmix::term_nigmix(void)
  {
  type = "term_nigmix";

  // Einlesen der Starwerte aus externer Datei
  startdata = stroption("startdata");

  // Startwerte für die shhrinkgeeffekte
  effectstart = doubleoption("effect",1E8,-1E7,1E7) ;

  // Startwert für Indicator (1. Komponente des Varianzparameters)
  indicatorstart = intoption("I",1,0,1);

  // Startwert für t^2 (2. Komponente des Varianzparameters)
  t2start = doubleoption("t2",11,0,10000000);

  // Startwert fuer die Mischungskomponente 
  omegastart = doubleoption("w",0.5,0,1);  
  
  // Lage der Punktmassen des Indikators
  v0 = doubleoption("v0",0.005,0,10000000);
  v1 = doubleoption("v1",1,0,10000000);
  
  // Hyperparameter der Priori fuer Shrinkageparameter
  a_t2 = doubleoption("a",5,0,500);
  b_t2 = doubleoption("b",50,0,500);
  
  // Hyperparameter der Priori fuer Mischungskomponente
  a_omega = doubleoption("aw",1,0,500);
  b_omega = doubleoption("bw",1,0,500);

  // Feste Werte für die Komponenten
  omegafix = simpleoption("wfix",false);

  // Varinazspezifische Werte für die Komponenten
  omegaadaptive = simpleoption("adaptive",false);
 
  }


// FUNCTION: setdefault
//---------------------
void term_nigmix::setdefault(void)
  {
  // call setdefault-methods of the options
  startdata.setdefault();
  effectstart.setdefault();
  indicatorstart.setdefault();
  t2start.setdefault();
  omegastart.setdefault();
  v0.setdefault();
  v1.setdefault();
  a_t2.setdefault();
  b_t2.setdefault();
  a_omega.setdefault();
  b_omega.setdefault();
  omegafix.setdefault();
  omegaadaptive.setdefault();
  }


// FUNCTION: check
//-----------------
bool term_nigmix::check(term & t)
  {
  if ( (t.varnames.size() == 1) && (t.options.size()<=14) && (t.options.size()>=1))   // SET: Anzahl Optionen
    {
    // extract options
    // prior coorespponds to Normal-Indicator-Mixing
    if (t.options[0] == "nigmix")
      t.type = "nigmix";
    
    else
      {
      setdefault();
      return false;
      }

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&startdata);
    optlist.push_back(&effectstart);
    optlist.push_back(&indicatorstart);
    optlist.push_back(&t2start);
    optlist.push_back(&omegastart);
    optlist.push_back(&v0);
    optlist.push_back(&v1);
    optlist.push_back(&a_t2);
    optlist.push_back(&b_t2);
    optlist.push_back(&a_omega);
    optlist.push_back(&b_omega);
    optlist.push_back(&omegafix);
    optlist.push_back(&omegaadaptive);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(14);                        // SET: Anzahl Optionen
    t.options[0] = t.type;
    t.options[1] = startdata.getvalue();
    t.options[2] = ST::doubletostring(effectstart.getvalue());
    t.options[3] = ST::inttostring(indicatorstart.getvalue());
    t.options[4] = ST::doubletostring(t2start.getvalue());
    t.options[5] = ST::doubletostring(omegastart.getvalue());    
    t.options[6] = ST::doubletostring(v0.getvalue());
    t.options[7] = ST::doubletostring(v1.getvalue());
    t.options[8] = ST::doubletostring(a_t2.getvalue());
    t.options[9] = ST::doubletostring(b_t2.getvalue());
    t.options[10] = ST::doubletostring(a_omega.getvalue());
    t.options[11] = ST::doubletostring(b_omega.getvalue());

    if (omegafix.getvalue()==false)
       t.options[12] = "false";
     else
       t.options[12] = "true";
    if (omegaadaptive.getvalue()==false)
       t.options[13] = "false";
     else
       t.options[13] = "true";

    setdefault();
    return true;
    }

  else
    {
    setdefault();
    return false;
    }
  }


// FUNCTION: checkvector
//----------------------
bool term_nigmix::checkvector(const vector<term> & terms, const unsigned & i)
  {
  assert(i< terms.size());
  if (terms[i].type == "nigmix")
    return true;
  return false;
  }

//------------------------------------------------------------------------------
//---------- class term_autoreg: implementation of member functions ------------
//------------------------------------------------------------------------------


term_autoreg::term_autoreg(void)
  {

  type = "term_autoreg";
  min=intoption("min",1,1,500);
  max=intoption("max",1,1,500);
  minvar=intoption("minvar",1,1,500);
  maxvar=intoption("maxvar",1,1,500);
  startv = doubleoption("startv",0.05,0.00001,1000);
  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("cp");
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,10000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,10000000);
  lambdastart = doubleoption("lambdastart",-1,-1,10000000);
  stationary = simpleoption("stationary",false);
  alpha = doubleoption("alpha",0.9,-1.0,1.0);
  alphafix = simpleoption("alphafix",false);
  center=simpleoption("center",false);
  vector<ST::string> adm_centerm;
  adm_centerm.push_back("mean");
  adm_centerm.push_back("samplecentered");
  adm_centerm.push_back("meanintercept");
  centermethod = stroption("centermethod",adm_centerm,"mean");
  }

void term_autoreg::setdefault(void)
  {
  min.setdefault();
  max.setdefault();
  minvar.setdefault();
  maxvar.setdefault();
  startv.setdefault();
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  alpha.setdefault();
  stationary.setdefault();
  alphafix.setdefault();
  center.setdefault();
  centermethod.setdefault();
  }


bool term_autoreg::check(term & t)
  {

  if ( (t.varnames.size() <= 2) && (t.varnames.size() >= 1) &&
       (t.options.size()<=21) && (t.options.size() >= 1) )
    {

    if (t.options[0] == "rw1" && t.varnames.size() == 1)
      t.type = "rw1";
    else if (t.options[0] == "rw2" && t.varnames.size() == 1)
      t.type = "rw2";
    else if (t.options[0] == "trw1" && t.varnames.size() == 1)
      t.type = "trw1";
    else if (t.options[0] == "trw2" && t.varnames.size() == 1)
      t.type = "trw2";
    else if (t.options[0] == "rw1vrw1" && t.varnames.size() == 1)
      t.type = "rw1vrw1";
    else if (t.options[0] == "rw2vrw1" && t.varnames.size() == 1)
      t.type = "rw2vrw1";
    else if (t.options[0] == "rw1vrw2" && t.varnames.size() == 1)
      t.type = "rw1vrw2";
    else if (t.options[0] == "rw2vrw2" && t.varnames.size() == 1)
      t.type = "rw2vrw2";
    else if (t.options[0] == "rw1" && t.varnames.size() == 2)
      t.type = "varcoeffrw1";
    else if (t.options[0] == "rw2" && t.varnames.size() == 2)
      t.type = "varcoeffrw2";
    else
      {
      setdefault();
      return false;
      }

    long minim,maxim;
    double minl, maxl, startl;

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&minvar);
    optlist.push_back(&maxvar);
    optlist.push_back(&startv);
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&stationary);
    optlist.push_back(&alpha);
    optlist.push_back(&alphafix);
    optlist.push_back(&center);
    optlist.push_back(&centermethod);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(21);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(min.getvalue());
    t.options[2] = ST::inttostring(max.getvalue());
    t.options[3] = ST::inttostring(minvar.getvalue());
    t.options[4] = ST::inttostring(maxvar.getvalue());
    t.options[5] = ST::doubletostring(startv.getvalue());
    t.options[6] = ST::doubletostring(lambda.getvalue());
    t.options[7] = ST::doubletostring(a.getvalue());
    t.options[8] = ST::doubletostring(b.getvalue());
    t.options[9] = proposal.getvalue();
    t.options[10] = ST::inttostring(updateW.getvalue());
    if (updatetau.getvalue()==false)
      t.options[11] = "false";
    else
      t.options[11] = "true";
    t.options[12] = ST::doubletostring(f.getvalue());
    t.options[13] = ST::doubletostring(lambdamin.getvalue());
    t.options[14] = ST::doubletostring(lambdamax.getvalue());
    t.options[15] = ST::doubletostring(lambdastart.getvalue());
    if(stationary.getvalue() == false)
      t.options[16] = "false";
    else
      t.options[16] = "true";
    t.options[17] = ST::doubletostring(alpha.getvalue());
    if(alphafix.getvalue() == false)
      t.options[18] = "false";
    else
      t.options[18] = "true";

    if(center.getvalue() == false)
      t.options[19] = "false";
    else
      t.options[19] = "true";

    t.options[20] = centermethod.getvalue();

    if (t.options[1].strtolong(minim) == 1)
      {
      setdefault();
      return false;
      }


    if (minim < 1)
      {
      setdefault();
      return false;
      }

    if (t.options[2].strtolong(maxim) == 1)
      {
      setdefault();
      return false;
      }

    if (maxim < minim)
      {
      setdefault();
      return false;
      }

    if (t.options[3].strtolong(minim) == 1)
      {
      setdefault();
      return false;
      }

    if (minim < 1)
      {
      setdefault();
      return false;
      }

    if (t.options[4].strtolong(maxim) == 1)
      {
      setdefault();
      return false;
      }

    if (maxim < minim)
      {
      setdefault();
      return false;
      }

    // stepwise

    int b = t.options[13].strtodouble(minl);
    b = t.options[14].strtodouble(maxl);
    b = t.options[15].strtodouble(startl);

    if (b==1)
      {
      setdefault();
      return false;
      }

    if (minl >= maxl)
      {
      setdefault();
      return false;
      }

    if (maxl < startl)
      {
      setdefault();
      return false;
      }

    // END: stepwise

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//---------- class term_season: implementation of member functions ------------
//------------------------------------------------------------------------------

term_season::term_season(void)
  {

  type = "term_season";
  period = intoption("period",12,2,72);
  min=intoption("min",1,1,500);
  max=intoption("max",1,1,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("cp");
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,10000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,10000000);
  lambdastart = doubleoption("lambdastart",10000,0,10000000);
  uniformprior = simpleoption("uniformprior",false);

  }


void term_season::setdefault(void)
  {
  period.setdefault();
  min.setdefault();
  max.setdefault();
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  uniformprior.setdefault();
  }


bool term_season::check(term & t)
  {
  if ( (t.varnames.size()<=2)  && (t.options.size()<=16) &&
       (t.options.size() >= 1) )
    {

    if ( (t.options[0] == "season") && (t.varnames.size() == 1) )
      t.type = "season";
    else if( (t.options[0] == "season") && (t.varnames.size() == 2) )
      t.type = "varcoeffseason";
    else
      {
      setdefault();
      return false;
      }

    long minim,maxim,per;
    double l,minl,maxl;

    vector<ST::string> opt;
    optionlist optlist;

    optlist.push_back(&period);
    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&uniformprior);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(16);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(period.getvalue());
    t.options[2] = ST::inttostring(min.getvalue());
    t.options[3] = ST::inttostring(max.getvalue());
    t.options[4] = ST::doubletostring(lambda.getvalue());
    t.options[5] = ST::doubletostring(a.getvalue());
    t.options[6] = ST::doubletostring(b.getvalue());
    t.options[7] = proposal.getvalue();
    t.options[8] = ST::inttostring(updateW.getvalue());
    if (updatetau.getvalue()==false)
      t.options[9] = "false";
    else
      t.options[9] = "true";
    t.options[10] = ST::doubletostring(f.getvalue());
    t.options[11] = ST::doubletostring(lambdamin.getvalue());
    t.options[12] = ST::doubletostring(lambdamax.getvalue());
    t.options[13] = ST::doubletostring(lambdastart.getvalue());
    if (uniformprior.getvalue() == false)
      t.options[14] = "false";
    else
      t.options[14] = "true";


    if (t.options[1].strtolong(per) == 1)
      {
      setdefault();
      return false;
      }

    if (t.options[2].strtolong(minim) == 1)
      {
      setdefault();
      return false;
      }

    if (minim < 1)
      {
      setdefault();
      return false;
      }

    if (t.options[3].strtolong(maxim) == 1)
      {
      setdefault();
      return false;
      }

    if (maxim < minim)
      {
      setdefault();
      return false;
      }

    if (t.options[4].strtodouble(l) == 1)
      {
      setdefault();
      return false;
      }

    // stepwise

    int b = t.options[11].strtodouble(minl);
    b = t.options[12].strtodouble(maxl);

    if (b==1)
      {
      setdefault();
      return false;
      }


    if (minl >= maxl)
      {
      setdefault();
      return false;
      }

    // END: stepwise

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//----------- class term_pspline: implementation of member functions -----------
//------------------------------------------------------------------------------

term_pspline::term_pspline(void)
  {
  type = "term_psline";
  min=intoption("min",0,1,100);
  max=intoption("max",0,1,100);
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  uniformb = simpleoption("uniformb",false);
  gridsize = intoption("gridsize",-1,10,2000);
  minvar=intoption("minvar",1,1,500);
  maxvar=intoption("maxvar",1,1,500);
  startv = doubleoption("startv",0.05,0.00001,1000);
  vector<ST::string> adm_prop;
  adm_prop.push_back("cp");
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  vector<ST::string> adm;
  adm.push_back("unrestricted");
  adm.push_back("increasing");
  adm.push_back("decreasing");
  monotone = stroption("monotone",adm,"unrestricted");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  diagtransform = simpleoption("diagtransform",false);
  derivative = simpleoption("derivative",false);
  bsplinebasis = simpleoption("bsplinebasis",false);
  contourprob = intoption("contourprob",-1,0,6);
  uniformprior = simpleoption("uniformprior",false);
  beta_0 = stroption("beta_0");
  discrete = simpleoption("discrete",false);
  df = intoption("df",20,3,50);
  stationary = simpleoption("stationary",false);
  alpha = doubleoption("alpha",0.9,-1.0,1.0);
  alphafix = simpleoption("alphafix",false);
  vector<ST::string> knotsdef;
  knotsdef.push_back("equidistant");
  knotsdef.push_back("quantiles");
  knots = stroption("knots",knotsdef,"equidistant");
//  lambdamin = doubleoption("lambdamin",0.0001,0.000001,10000000);
//  lambdamax = doubleoption("lambdamax",10000,0.000001,10000000);
//  lambdastart = doubleoption("lambdastart",-1,-1,10000000);
  lowerknot = doubleoption("lowerknot",0,-10000000,10000000);
  upperknot = doubleoption("upperknot",0,-10000000,10000000);
  merrorvar = doubleoption("merrorvar",0,0,10000000);
  lowergrid = doubleoption("lowergrid",0,-10000000,10000000);
  uppergrid = doubleoption("uppergrid",0,-10000000,10000000);
  discretize = simpleoption("discretize", false);
  digits = intoption("digits",2,0,5);
  nobs = intoption("nobs",0,0,10000000);
  vector<ST::string> adm_centerm;
  adm_centerm.push_back("mean");
  adm_centerm.push_back("samplecentered");
  adm_centerm.push_back("meanintercept");
  centermethod = stroption("centermethod",adm_centerm,"mean");
  }

void term_pspline::setdefault(void)
  {
  min.setdefault();
  max.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  uniformb.setdefault();
  gridsize.setdefault();
  minvar.setdefault();
  maxvar.setdefault();
  startv.setdefault();
  proposal.setdefault();
  monotone.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  diagtransform.setdefault();
  derivative.setdefault();
  bsplinebasis.setdefault();
  contourprob.setdefault();
  uniformprior.setdefault();
  beta_0.setdefault();
  discrete.setdefault();
  df.setdefault();
  alpha.setdefault();
  stationary.setdefault();
  alphafix.setdefault();
  knots.setdefault();
  lowerknot.setdefault();
  upperknot.setdefault();
  merrorvar.setdefault();
  lowergrid.setdefault();
  uppergrid.setdefault();
  discretize.setdefault();
  digits.setdefault();
  nobs.setdefault();
  centermethod.setdefault();
  }

bool term_pspline::check(term & t)
  {

  if ( (t.varnames.size()==1)  && (t.options.size() >= 1)
        && (t.options.size() <= 39) )
    {

    if (t.options[0] == "psplinerw1")
      t.type = "psplinerw1";
    else if (t.options[0] == "psplinerw2")
      t.type = "psplinerw2";
    else if (t.options[0] == "tpsplinerw1")
      t.type = "tpsplinerw1";
    else if (t.options[0] == "tpsplinerw2")
      t.type = "tpsplinerw2";
    else if (t.options[0] == "psplinerw1vrw1")
      t.type = "psplinerw1vrw1";
    else if (t.options[0] == "psplinerw1vrw2")
      t.type = "psplinerw1vrw2";
    else if (t.options[0] == "psplinerw2vrw1")
      t.type = "psplinerw2vrw1";
    else if (t.options[0] == "psplinerw2vrw2")
      t.type = "psplinerw2vrw2";
    else
      {
      setdefault();
      return false;
      }

    long minim,maxim;
//    double maxl,minl,startl;

    optionlist optlist;
    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&uniformb);
    optlist.push_back(&gridsize);
    optlist.push_back(&minvar);
    optlist.push_back(&maxvar);
    optlist.push_back(&startv);
    optlist.push_back(&proposal);
    optlist.push_back(&monotone);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&diagtransform);
    optlist.push_back(&derivative);
    optlist.push_back(&bsplinebasis);
    optlist.push_back(&contourprob);
    optlist.push_back(&uniformprior);
    optlist.push_back(&beta_0);
    optlist.push_back(&discrete);
    optlist.push_back(&df);
    optlist.push_back(&stationary);
    optlist.push_back(&alpha);
    optlist.push_back(&alphafix);
    optlist.push_back(&knots);
    optlist.push_back(&lowerknot);
    optlist.push_back(&upperknot);
    optlist.push_back(&merrorvar);
    optlist.push_back(&lowergrid);
    optlist.push_back(&uppergrid);
    optlist.push_back(&discretize);
    optlist.push_back(&digits);
    optlist.push_back(&nobs);
    optlist.push_back(&centermethod);    


    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

   t.options.erase(t.options.begin(),t.options.end());
   t.options = vector<ST::string>(39);
   t.options[0] = t.type;
   t.options[1] = ST::inttostring(min.getvalue());
   t.options[2] = ST::inttostring(max.getvalue());
   t.options[3] = ST::inttostring(degree.getvalue());
   t.options[4] = ST::inttostring(numberknots.getvalue());
   t.options[5] = ST::doubletostring(lambda.getvalue());
   t.options[6] = ST::doubletostring(a.getvalue());
   t.options[7] = ST::doubletostring(b.getvalue());
   if (uniformb.getvalue() == false)
     t.options[8] = "false";
   else
     t.options[8] = "true";
   t.options[9] = ST::inttostring(gridsize.getvalue());
   t.options[10] = ST::inttostring(minvar.getvalue());
   t.options[11] = ST::inttostring(maxvar.getvalue());
   t.options[12] = ST::doubletostring(startv.getvalue());
   t.options[13] = proposal.getvalue();
   t.options[14] = monotone.getvalue();
   t.options[15] = ST::inttostring(updateW.getvalue());
   if (updatetau.getvalue() == false)
     t.options[16] = "false";
   else
     t.options[16] = "true";
   t.options[17] = ST::doubletostring(f.getvalue());
   if (diagtransform.getvalue() == false)
     t.options[18] = "false";
   else
     t.options[18] = "true";
   if (derivative.getvalue() == false)
     t.options[19] = "false";
   else
     t.options[19] = "true";
   if (bsplinebasis.getvalue() == false)
     t.options[20] = "false";
   else
     t.options[20] = "true";
   t.options[21] = ST::inttostring(contourprob.getvalue());
   if(uniformprior.getvalue() == false)
     t.options[22] = "false";
   else
     t.options[22] = "true";
   t.options[23] = beta_0.getvalue();
   if(discrete.getvalue() == false)
     t.options[24] = "false";
   else
     t.options[24] = "true";
   t.options[25] = ST::inttostring(df.getvalue());
   if(stationary.getvalue() == false)
     t.options[26] = "false";
   else
     t.options[26] = "true";
   t.options[27] = ST::doubletostring(alpha.getvalue());
   if(alphafix.getvalue() == false)
     t.options[28] = "false";
   else
     t.options[28] = "true";
   t.options[29] = knots.getvalue();
    t.options[30] = ST::doubletostring(lowerknot.getvalue());
    t.options[31] = ST::doubletostring(upperknot.getvalue());
    t.options[32] = ST::doubletostring(merrorvar.getvalue());
    t.options[33] = ST::doubletostring(lowergrid.getvalue());
    t.options[34] = ST::doubletostring(uppergrid.getvalue());
   if (discretize.getvalue() == false)
     t.options[35] = "false";
   else
     t.options[35] = "true";
   t.options[36] = ST::inttostring(digits.getvalue());
   t.options[37] = ST::inttostring(nobs.getvalue());

   t.options[38] = centermethod.getvalue();

   if (t.options[1].strtolong(minim) == 1)
     {
     setdefault();
     return false;
     }

   if (t.options[2].strtolong(maxim) == 1)
     {
     setdefault();
     return false;
     }

   if (maxim < minim)
     {
     setdefault();
     return false;
     }

   if (lambda.getvalue() < 0)
     {
     setdefault();
     return false;
     }

   if ( contourprob.getvalue()-1 > degree.getvalue())
     {
     setdefault();
     return false;
     }


    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//----------- class term_spatial: implementation of member functions -----------
//------------------------------------------------------------------------------

term_spatial::term_spatial(void)
  {
  type = "term_spatial";
  map=stroption("map");
  min=intoption("min",1,1,500);
  max=intoption("max",1,1,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("cp");
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"cp");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,10000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,10000000);
  lambdastart = doubleoption("lambdastart",10000,0,10000000);
  uniformprior = simpleoption("uniformprior",false);
  nrrows = intoption("nrrows",2,0,100);
  Laplace = simpleoption("Laplace",false);
  stationary = simpleoption("stationary",false);
  alpha = doubleoption("alpha",0.9,-1.0,1.0);
  alphafix = simpleoption("alphafix",false);
  center = simpleoption("center",false);
  vector<ST::string> adm_centerm;
  adm_centerm.push_back("mean");
  adm_centerm.push_back("samplecentered");
  adm_centerm.push_back("meanintercept");
  centermethod = stroption("centermethod",adm_centerm,"mean");
  }

void term_spatial::setdefault(void)
  {
  map.setdefault();
  min.setdefault();
  max.setdefault();
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  uniformprior.setdefault();
  nrrows.setdefault();
  Laplace.setdefault();
  alpha.setdefault();
  stationary.setdefault();
  alphafix.setdefault();
  center.setdefault();
  centermethod.setdefault();
  }


bool term_spatial::check(term & t)
  {

  if ( (t.varnames.size()<=2)  && (t.varnames.size()>=1) &&
       (t.options.size()<=22) && (t.options.size() >= 1) )
    {

    if (t.options[0] == "spatial" && t.varnames.size()==1)
      t.type = "spatial";
    else  if (t.options[0] == "tspatial" && t.varnames.size()==1)
      t.type = "tspatial";
    else  if (t.options[0] == "spatial" && t.varnames.size()==2)
      t.type = "varcoeffspatial";
    else  if (t.options[0] == "tspatial" && t.varnames.size()==2)
      t.type = "tvarcoeffspatial";
    else
      {
      setdefault();
      return false;
      }

    long minim,maxim;
    double minl,maxl;

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&map);
    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&uniformprior);
    optlist.push_back(&nrrows);
    optlist.push_back(&Laplace);
    optlist.push_back(&stationary);
    optlist.push_back(&alpha);
    optlist.push_back(&alphafix);
    optlist.push_back(&center);
    optlist.push_back(&centermethod);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(22);
    t.options[0] = t.type;
    t.options[1] = map.getvalue();
    t.options[2] = ST::inttostring(min.getvalue());
    t.options[3] = ST::inttostring(max.getvalue());
    t.options[4] = ST::doubletostring(lambda.getvalue());
    t.options[5] = ST::doubletostring(a.getvalue());
    t.options[6] = ST::doubletostring(b.getvalue());
    t.options[7] = proposal.getvalue();
    t.options[8] = ST::inttostring(updateW.getvalue());
    if (updatetau.getvalue()==false)
      t.options[9] = "false";
    else
      t.options[9] = "true";
    t.options[10] = ST::doubletostring(f.getvalue());
    t.options[11] = ST::doubletostring(lambdamin.getvalue());
    t.options[12] = ST::doubletostring(lambdamax.getvalue());
    t.options[13] = ST::doubletostring(lambdastart.getvalue());
    if (uniformprior.getvalue() == false)
      t.options[14] = "false";
    else
      t.options[14] = "true";
    t.options[15] = ST::inttostring(nrrows.getvalue());
    if (Laplace.getvalue()==false)
      t.options[16] = "false";
    else
      t.options[16] = "true";
    if(stationary.getvalue() == false)
      t.options[17] = "false";
    else
      t.options[17] = "true";
    t.options[18] = ST::doubletostring(alpha.getvalue());
    if(alphafix.getvalue() == false)
      t.options[19] = "false";
    else
      t.options[19] = "true";
    if(center.getvalue() == false)
      t.options[20] = "false";
    else
      t.options[20] = "true";

    t.options[21] = centermethod.getvalue();

    if (t.options[2].strtolong(minim) == 1)
      {
      setdefault();
      return false;
      }

    if (t.options[3].strtolong(maxim) == 1)      {
      setdefault();
      return false;
      }

    if (maxim < minim)
      {
      setdefault();
      return false;
      }

    // stepwise

    int b = t.options[11].strtodouble(minl);
    b = t.options[12].strtodouble(maxl);

    if (b==1)
      {
      setdefault();
      return false;
      }

    if (minl >= maxl)
      {
      setdefault();
      return false;
      }

    // END: stepwise

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }



//------------------------------------------------------------------------------
//----------- class term_spatial: implementation of member functions -----------
//------------------------------------------------------------------------------

term_spatialxy::term_spatialxy(void)
  {
  type = "term_spatialxy";
  min=intoption("min",1,1,500);
  max=intoption("max",1,1,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  maxdist = doubleoption("maxdist",1,0,100000000);
  vector<ST::string> adm_centerm;
  adm_centerm.push_back("mean");
  adm_centerm.push_back("samplecentered");
  adm_centerm.push_back("meanintercept");
  centermethod = stroption("centermethod",adm_centerm,"mean");
  }

void term_spatialxy::setdefault(void)
  {
  min.setdefault();
  max.setdefault();
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  maxdist.setdefault();
  centermethod.setdefault();
  }

bool term_spatialxy::check(term & t)
  {

  if ( (t.varnames.size()==2)  && (t.options.size()<=8) &&
       (t.options.size() >= 1) )
    {

    if (t.options[0] == "spatialxy")
      t.type = "spatialxy";
    else
      {
      setdefault();
      return false;
      }

    long minim,maxim;

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&maxdist);
    optlist.push_back(&centermethod);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(8);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(min.getvalue());
    t.options[2] = ST::inttostring(max.getvalue());
    t.options[3] = ST::doubletostring(lambda.getvalue());
    t.options[4] = ST::doubletostring(a.getvalue());
    t.options[5] = ST::doubletostring(b.getvalue());
    t.options[6] = ST::doubletostring(maxdist.getvalue());
    t.options[7] = centermethod.getvalue();

    maxim=max.getvalue();
    minim=min.getvalue();

    if (maxim < minim)
      {
      setdefault();
      return false;
      }

    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//----- class term_geokriging: implementation of member functions ------
//------------------------------------------------------------------------------


term_geokriging::term_geokriging(void)
  {
  type = "term_geokriging";
  numberknots=intoption("nrknots",50,5,500);
  nu=doubleoption("nu",1.5,0.5,3.5);
  maxdist=doubleoption("maxdist",-1,0.00001,10000);
  full=simpleoption("full",false);
  knotdata=stroption("knotdata");
  p=doubleoption("p",-20,-1000,-0.0001);
  q=doubleoption("q",20,0.0001,1000);
  maxsteps=intoption("maxsteps",100,1,10000);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",0.1,0,10000000);
  map=stroption("map");

  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  uniformprior = simpleoption("uniformprior",false);
  }


void term_geokriging::setdefault(void)
  {
  numberknots.setdefault();
  nu.setdefault();
  maxdist.setdefault();
  full.setdefault();
  knotdata.setdefault();
  p.setdefault();
  q.setdefault();
  maxsteps.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  map.setdefault();

  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  uniformprior.setdefault();
  }


bool term_geokriging::check(term & t)
  {

  if ( (t.varnames.size()==1)  && (t.options.size() >=1)
        && (t.options.size() <= 19) )
    {

    if (t.options[0] == "geokriging")
      t.type = "geokriging";
    else
      {
      setdefault();
      return false;
      }

    optionlist optlist;
    optlist.push_back(&numberknots);
    optlist.push_back(&nu);
    optlist.push_back(&maxdist);
    optlist.push_back(&full);
    optlist.push_back(&knotdata);
    optlist.push_back(&p);
    optlist.push_back(&q);
    optlist.push_back(&maxsteps);
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&map);

    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&uniformprior);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(19);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(numberknots.getvalue());
    t.options[2] = ST::doubletostring(nu.getvalue());
    t.options[3] = ST::doubletostring(maxdist.getvalue());
    if(full.getvalue()==true)
      {
      t.options[4] = "true";
      }
    else
      {
      t.options[4] = "false";
      }
    t.options[5] = knotdata.getvalue();
    t.options[6] = ST::doubletostring(p.getvalue());
    t.options[7] = ST::doubletostring(q.getvalue());
    t.options[8] = ST::inttostring(maxsteps.getvalue());
    t.options[9] = ST::doubletostring(lambda.getvalue());
    t.options[10] = ST::doubletostring(lambdastart.getvalue());
    t.options[11] = map.getvalue();

    t.options[12] = ST::doubletostring(a.getvalue());
    t.options[13] = ST::doubletostring(b.getvalue());
    t.options[14] = proposal.getvalue();
    t.options[15] = ST::inttostring(updateW.getvalue());
    if(updatetau.getvalue()==true)
      t.options[16] = "true";
    else
      t.options[16] = "false";
    t.options[17] = ST::doubletostring(f.getvalue());
    if (uniformprior.getvalue() == false)
      t.options[18] = "false";
    else
      t.options[18] = "true";

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//------ class term_interactpspline: implementation of member functions --------
//------------------------------------------------------------------------------

term_interactpspline::term_interactpspline(void)
  {
  type = "term_interactpspline";
  min=intoption("min",0,1,100);
  max=intoption("max",0,1,100);
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  reduced = simpleoption("reduced",false);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  singleblock = simpleoption("singleblock",false);
  gridsize = intoption("gridsize",-1,10,35);
  vector<ST::string> adm_prop;
  adm_prop.push_back("cp");
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  uniformprior = simpleoption("uniformprior",false);
  blocksize=intoption("blocksize",6,2,100);
  center=simpleoption("center",false);
  vector<ST::string> adm_centerm;
  adm_centerm.push_back("mean");
  adm_centerm.push_back("samplecentered");
  adm_centerm.push_back("meanintercept");
  centermethod = stroption("centermethod",adm_centerm,"mean");

  }


void term_interactpspline::setdefault(void)
  {
  min.setdefault();
  max.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  reduced.setdefault();
  a.setdefault();
  b.setdefault();
  singleblock.setdefault();
  gridsize.setdefault();
  proposal.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  uniformprior.setdefault();
  blocksize.setdefault();
  center.setdefault();
  centermethod.setdefault();
  }

bool term_interactpspline::check(term & t)
  {

  optionlist optlist;
  optlist.push_back(&min);
  optlist.push_back(&max);
  optlist.push_back(&degree);
  optlist.push_back(&numberknots);
  optlist.push_back(&lambda);
  optlist.push_back(&reduced);
  optlist.push_back(&a);
  optlist.push_back(&b);
  optlist.push_back(&singleblock);
  optlist.push_back(&gridsize);
  optlist.push_back(&proposal);
  optlist.push_back(&updateW);
  optlist.push_back(&updatetau);
  optlist.push_back(&f);
  optlist.push_back(&uniformprior);
  optlist.push_back(&blocksize);
  optlist.push_back(&center);
  optlist.push_back(&centermethod);


  if ( (t.varnames.size()<=3)  && (t.options.size() >= 1)
        && (t.options.size() <= 19) )
    {

    if (t.options[0] == "pspline2dimrw1" && t.varnames.size()==2)
      t.type = "pspline2dimrw1";
    else if (t.options[0] == "pspline2dimrw2" && t.varnames.size()==2)
      t.type = "pspline2dimrw2";
    else if (t.options[0] == "tpspline2dimrw1" && t.varnames.size()==2)
      t.type = "tpspline2dimrw1";
    else if (t.options[0] == "pspline2dimband" && t.varnames.size()==2)
      t.type = "pspline2dimband";
    else if (t.options[0] == "tpspline2dimband" && t.varnames.size()==2)
      t.type = "tpspline2dimband";
    else if (t.options[0] == "psplinekrrw1" && t.varnames.size()==2)
      t.type = "psplinekrrw1";
    else if (t.options[0] == "psplinekrrw2" && t.varnames.size()==2)
      t.type = "psplinekrrw2";
    else if (t.options[0] == "pspline2dimrw1" && t.varnames.size()==3)
      t.type = "varpspline2dimrw1";
    else if (t.options[0] == "pspline2dimrw2" && t.varnames.size()==3)
      t.type = "varpspline2dimrw2";
    else if (t.options[0] == "pspline2dimband" && t.varnames.size()==3)
      t.type = "varpspline2dimband";
    else if (t.options[0] == "psplinekrrw1" && t.varnames.size()==3)
      t.type = "varpsplinekrrw1";
    else if (t.options[0] == "psplinekrrw2" && t.varnames.size()==3)
      t.type = "varpsplinekrrw2";
    else
      {
      setdefault();
      return false;
      }

    long minim,maxim;

    unsigned i;

    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

   t.options.erase(t.options.begin(),t.options.end());
   t.options = vector<ST::string>(19);
   t.options[0] = t.type;
   t.options[1] = ST::inttostring(min.getvalue());
   t.options[2] = ST::inttostring(max.getvalue());
   t.options[3] = ST::inttostring(degree.getvalue());
   t.options[4] = ST::inttostring(numberknots.getvalue());
   t.options[5] = ST::doubletostring(lambda.getvalue());
   if (reduced.getvalue() == false)
     t.options[6] = "false";
   else
     t.options[6] = "true";
   t.options[7] = ST::doubletostring(a.getvalue());
   t.options[8] = ST::doubletostring(b.getvalue());
   if (singleblock.getvalue() == false)
     t.options[9] = "false";
   else
     t.options[9] = "true";
   t.options[10] = ST::inttostring(gridsize.getvalue());
   t.options[11] = proposal.getvalue();
   t.options[12] = ST::inttostring(updateW.getvalue());
   if (updatetau.getvalue() == false)
     t.options[13] = "false";
   else
     t.options[13] = "true";
   t.options[14] = ST::doubletostring(f.getvalue());
   if (uniformprior.getvalue() == false)
     t.options[15] = "false";
   else
     t.options[15] = "true";
   t.options[16] = ST::inttostring(blocksize.getvalue());
   if (center.getvalue() == false)
     t.options[17] = "false";
   else
     t.options[17] = "true";

   t.options[18] = centermethod.getvalue();

   if (t.options[1].strtolong(minim) == 1)
     {
     setdefault();
     return false;
     }

   if (t.options[2].strtolong(maxim) == 1)
     {
     setdefault();
     return false;
     }

   if (maxim < minim)
     {
     setdefault();
     return false;
     }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


bool term_interactpspline::checkvector(const vector<term> & terms,
                                       const unsigned & i)
  {

  assert(i< terms.size());

  if ((terms[i].type == "pspline2dimrw1") || (terms[i].type == "pspline2dimrw2")
     || (terms[i].type == "psplinekrrw1") || (terms[i].type == "psplinekrrw2")
     || (terms[i].type == "varpspline2dimrw1") || (terms[i].type == "varpspline2dimrw2")
     || (terms[i].type == "varpsplinekrrw1") || (terms[i].type == "varpsplinekrrw2")
     || (terms[i].type == "tpspline2dimrw1")
     || (terms[i].type == "pspline2dimband") || (terms[i].type == "tpspline2dimband")
     ||  (terms[i].type == "varpspline2dimband")
     )
    return true;

  return false;
  }


//------------------------------------------------------------------------------
//---------- class term_geospline: implementation of member functions ----------
//------------------------------------------------------------------------------

term_geospline::term_geospline(void)
  {
  type = "term_geospline";
  map=stroption("map");
  min=intoption("min",0,1,100);
  max=intoption("max",0,1,100);
  degree=intoption("degree",3,1,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  reduced = simpleoption("reduced",false);
  singleblock = simpleoption("singleblock",false);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("cp");
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  uniformprior = simpleoption("uniformprior",false);
  vector<ST::string> adm_centerm;
  adm_centerm.push_back("mean");
  adm_centerm.push_back("samplecentered");
  adm_centerm.push_back("meanintercept");
  centermethod = stroption("centermethod",adm_centerm,"mean");
  }


void term_geospline::setdefault(void)
  {
  map.setdefault();
  min.setdefault();
  max.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  reduced.setdefault();
  singleblock.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  uniformprior.setdefault();
  centermethod.setdefault();
  }

bool term_geospline::check(term & t)
  {

  if ( (t.varnames.size()==1)  && (t.options.size() >= 1)
        && (t.options.size() <= 17) )
    {

    if (t.options[0] == "geospline")
      t.type = "geospline";
    else if (t.options[0] == "geosplinerw1")
      t.type = "geospline";
    else if (t.options[0] == "geosplinerw2")
      t.type = "geosplinerw2";
    else
      {
      setdefault();
      return false;
      }

    long minim,maxim;

    optionlist optlist;
    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&reduced);
    optlist.push_back(&map);
    optlist.push_back(&singleblock);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&uniformprior);
    optlist.push_back(&centermethod);

    unsigned i;

    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

   t.options.erase(t.options.begin(),t.options.end());
   t.options = vector<ST::string>(17);
   t.options[0] = t.type;
   t.options[1] = ST::inttostring(min.getvalue());
   t.options[2] = ST::inttostring(max.getvalue());
   t.options[3] = ST::inttostring(degree.getvalue());
   t.options[4] = ST::inttostring(numberknots.getvalue());
   t.options[5] = ST::doubletostring(lambda.getvalue());
   if (reduced.getvalue() == false)
     t.options[6] = "false";
   else
     t.options[6] = "true";
   t.options[7] = map.getvalue();
   if (singleblock.getvalue() == false)
     t.options[8] = "false";
   else
     t.options[8] = "true";
   t.options[9] = ST::doubletostring(a.getvalue());
   t.options[10] = ST::doubletostring(b.getvalue());
   t.options[11] = proposal.getvalue();
   t.options[12] = ST::inttostring(updateW.getvalue());
   if (updatetau.getvalue() == false)
     t.options[13] = "false";
   else
     t.options[13] = "true";
   t.options[14] = ST::doubletostring(f.getvalue());
   if (uniformprior.getvalue() == false)
     t.options[15] = "false";
   else
     t.options[15] = "true";

   t.options[16] = centermethod.getvalue();

   t.options[1].strtolong(minim);
   t.options[2].strtolong(maxim);

   if (maxim < minim)
     {
     setdefault();
     return false;
     }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


bool term_geospline::checkvector(const vector<term> & terms, const unsigned & i)
  {

  assert(i< terms.size());

  if (terms[i].type == "geospline" || terms[i].type == "geosplinerw1" || terms[i].type == "geosplinerw2")
    return true;

  return false;
  }


//------------------------------------------------------------------------------
//------ class term_varcoeff_interactpspline: implementation of member functions
//------------------------------------------------------------------------------

/*
term_varcoeff_interactpspline::term_varcoeff_interactpspline(void)
  {
  type = "term_varcoeff_interactpspline";
  min=intoption("min",0,1,100);
  max=intoption("max",0,1,100);
  degree=intoption("degree",3,1,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  reduced = simpleoption("reduced",false);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  singleblock = simpleoption("singleblock",false);
  gridsize = intoption("gridsize",-1,10,35);
  vector<ST::string> adm_prop;
  adm_prop.push_back("cp");
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  uniformprior = simpleoption("uniformprior",false);
  }


void term_varcoeff_interactpspline::setdefault(void)
  {
  min.setdefault();
  max.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  reduced.setdefault();
  a.setdefault();
  b.setdefault();
  singleblock.setdefault();
  gridsize.setdefault();
  proposal.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  uniformprior.setdefault();
  }

bool term_varcoeff_interactpspline::check(term & t)
  {

  optionlist optlist;
  optlist.push_back(&min);
  optlist.push_back(&max);
  optlist.push_back(&degree);
  optlist.push_back(&numberknots);
  optlist.push_back(&lambda);
  optlist.push_back(&reduced);
  optlist.push_back(&a);
  optlist.push_back(&b);
  optlist.push_back(&singleblock);
  optlist.push_back(&gridsize);
  optlist.push_back(&proposal);
  optlist.push_back(&updateW);
  optlist.push_back(&updatetau);
  optlist.push_back(&f);
  optlist.push_back(&uniformprior);

  if ( (t.varnames.size()==2)  && (t.options.size() >= 1)
        && (t.options.size() <= 16) )
    {

    if (t.options[0] == "pspline2dimrw1")
      t.type = "varpspline2dimrw1";
    else if (t.options[0] == "tpspline2dimrw1")
      t.type = "vartpspline2dimrw1";
    else if (t.options[0] == "pspline2dimband")
      t.type = "varpspline2dimband";
    else if (t.options[0] == "tpspline2dimband")
      t.type = "vartpspline2dimband";
    else if (t.options[0] == "psplinekrrw1")
      t.type = "varpsplinekrrw1";
    else if (t.options[0] == "psplinekrrw2")
      t.type = "varpsplinekrrw2";
    else
      {
      setdefault();
      return false;
      }

    long minim,maxim;

    unsigned i;

    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

   t.options.erase(t.options.begin(),t.options.end());
   t.options = vector<ST::string>(16);
   t.options[0] = t.type;
   t.options[1] = ST::inttostring(min.getvalue());
   t.options[2] = ST::inttostring(max.getvalue());
   t.options[3] = ST::inttostring(degree.getvalue());
   t.options[4] = ST::inttostring(numberknots.getvalue());
   t.options[5] = ST::doubletostring(lambda.getvalue());
   if (reduced.getvalue() == false)
     t.options[6] = "false";
   else
     t.options[6] = "true";
   t.options[7] = ST::doubletostring(a.getvalue());
   t.options[8] = ST::doubletostring(b.getvalue());
   if (singleblock.getvalue() == false)
     t.options[9] = "false";
   else
     t.options[9] = "true";
   t.options[10] = ST::inttostring(gridsize.getvalue());
   t.options[11] = proposal.getvalue();
   t.options[12] = ST::inttostring(updateW.getvalue());
   if (updatetau.getvalue() == false)
     t.options[13] = "false";
   else
     t.options[13] = "true";
   t.options[14] = ST::doubletostring(f.getvalue());
   if (uniformprior.getvalue() == false)
     t.options[15] = "false";
   else
     t.options[15] = "true";

   if (t.options[1].strtolong(minim) == 1)
     {
     setdefault();
     return false;
     }

   if (t.options[2].strtolong(maxim) == 1)
     {
     setdefault();
     return false;
     }

   if (maxim < minim)
     {
     setdefault();
     return false;
     }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


bool term_varcoeff_interactpspline::checkvector(const vector<term> & terms,
                                       const unsigned & i)
  {

  assert(i< terms.size());

  if ((terms[i].type == "varpspline2dimrw1") || (terms[i].type == "varpsplinekrrw1")
     || (terms[i].type == "varpsplinekrrw2") || (terms[i].type == "vartpspline2dimrw1")
     || (terms[i].type == "varpspline2dimband") || (terms[i].type == "vartpspline2dimband")
     )
    return true;

  return false;
  }
*/

//------------------------------------------------------------------------------
//---------- class term_varcoeff_geospline: implementation of member functions -
//------------------------------------------------------------------------------

term_varcoeff_geospline::term_varcoeff_geospline(void)
  {
  type = "term_varcoeff_geospline";
  map=stroption("map");
  min=intoption("min",0,1,100);
  max=intoption("max",0,1,100);
  degree=intoption("degree",3,1,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  reduced = simpleoption("reduced",false);
  singleblock = simpleoption("singleblock",false);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("cp");
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  uniformprior = simpleoption("uniformprior",false);
  center = simpleoption("center",false);
  vector<ST::string> adm_centerm;
  adm_centerm.push_back("mean");
  adm_centerm.push_back("samplecentered");
  adm_centerm.push_back("meanintercept");
  centermethod = stroption("centermethod",adm_centerm,"mean");
  }


void term_varcoeff_geospline::setdefault(void)
  {
  map.setdefault();
  min.setdefault();
  max.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  reduced.setdefault();
  singleblock.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  uniformprior.setdefault();
  center.setdefault();
  centermethod.setdefault();
  }

bool term_varcoeff_geospline::check(term & t)
  {

  if ( (t.varnames.size()==2)  && (t.options.size() >= 1)
        && (t.options.size() <= 18) )
    {

    if (t.options[0] == "geospline")
      t.type = "vargeospline";
    else
      {
      setdefault();
      return false;
      }

    long minim,maxim;

    optionlist optlist;
    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&reduced);
    optlist.push_back(&map);
    optlist.push_back(&singleblock);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&uniformprior);
    optlist.push_back(&center);
    optlist.push_back(&centermethod);

    unsigned i;

    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

   t.options.erase(t.options.begin(),t.options.end());
   t.options = vector<ST::string>(18);
   t.options[0] = t.type;
   t.options[1] = ST::inttostring(min.getvalue());
   t.options[2] = ST::inttostring(max.getvalue());
   t.options[3] = ST::inttostring(degree.getvalue());
   t.options[4] = ST::inttostring(numberknots.getvalue());
   t.options[5] = ST::doubletostring(lambda.getvalue());
   if (reduced.getvalue() == false)
     t.options[6] = "false";
   else
     t.options[6] = "true";
   t.options[7] = map.getvalue();
   if (singleblock.getvalue() == false)
     t.options[8] = "false";
   else
     t.options[8] = "true";
   t.options[9] = ST::doubletostring(a.getvalue());
   t.options[10] = ST::doubletostring(b.getvalue());
   t.options[11] = proposal.getvalue();
   t.options[12] = ST::inttostring(updateW.getvalue());
   if (updatetau.getvalue() == false)
     t.options[13] = "false";
   else
     t.options[13] = "true";
   t.options[14] = ST::doubletostring(f.getvalue());
   if (uniformprior.getvalue() == false)
     t.options[15] = "false";
   else
     t.options[15] = "true";

   if (center.getvalue() == false)
     t.options[16] = "false";
   else
     t.options[16] = "true";

    t.options[17] = centermethod.getvalue();     

   t.options[1].strtolong(minim);
   t.options[2].strtolong(maxim);

   if (maxim < minim)
     {
     setdefault();
     return false;
     }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


bool term_varcoeff_geospline::checkvector(const vector<term> & terms,
                                          const unsigned & i)
  {

  assert(i< terms.size());

  if (terms[i].type == "vargeospline")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//------ class term_varcoeff_pspline: implementation of member functions -------
//------------------------------------------------------------------------------


term_varcoeff_pspline::term_varcoeff_pspline(void)
  {
  type = "term_varcoeff";
  min=intoption("min",0,1,100);
  max=intoption("max",0,1,100);
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  gridsize = intoption("gridsize",-1,10,500);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("cp");
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  vector<ST::string> adm;
  adm.push_back("unrestricted");
  adm.push_back("increasing");
  adm.push_back("decreasing");
  monotone = stroption("monotone",adm,"unrestricted");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  diagtransform = simpleoption("diagtransform",false);
  derivative = simpleoption("derivative",false);
  contourprob = intoption("contourprob",-1,0,6);
  uniformprior = simpleoption("uniformprior",false);
  beta_0 = stroption("beta_0");
  vector<ST::string> knotsdef;
  knotsdef.push_back("equidistant");
  knotsdef.push_back("quantiles");
  knots = stroption("knots",knotsdef,"equidistant");
  center = simpleoption("center",false);
  vector<ST::string> adm_centerm;
  adm_centerm.push_back("mean");
  adm_centerm.push_back("samplecentered");
  adm_centerm.push_back("meanintercept");
  centermethod = stroption("centermethod",adm_centerm,"mean");

  }


void term_varcoeff_pspline::setdefault(void)
  {
  min.setdefault();
  max.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  gridsize.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  monotone.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  diagtransform.setdefault();
  derivative.setdefault();
  contourprob.setdefault();
  uniformprior.setdefault();
  beta_0.setdefault();
  knots.setdefault();
  center.setdefault();
  centermethod.setdefault();
  }


bool term_varcoeff_pspline::check(term & t)
  {

  if ( (t.varnames.size()==2)  && (t.options.size() >=1)
        && (t.options.size() <= 22) )
    {

    if (t.options[0] == "psplinerw1")
      t.type = "varpsplinerw1";
    else if (t.options[0] == "psplinerw2")
      t.type = "varpsplinerw2";
    else
      {
      setdefault();
      return false;
      }

    optionlist optlist;
    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&gridsize);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&monotone);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&diagtransform);
    optlist.push_back(&derivative);
    optlist.push_back(&contourprob);
    optlist.push_back(&uniformprior);
    optlist.push_back(&beta_0);
    optlist.push_back(&knots);
    optlist.push_back(&center);
    optlist.push_back(&centermethod);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(22);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(min.getvalue());
    t.options[2] = ST::inttostring(max.getvalue());
    t.options[3] = ST::inttostring(degree.getvalue());
    t.options[4] = ST::inttostring(numberknots.getvalue());
    t.options[5] = ST::doubletostring(lambda.getvalue());
    t.options[6] = ST::inttostring(gridsize.getvalue());
    t.options[7] = ST::doubletostring(a.getvalue());
    t.options[8] = ST::doubletostring(b.getvalue());
    t.options[9] = proposal.getvalue();
    t.options[10] = monotone.getvalue();
    t.options[11] = ST::inttostring(updateW.getvalue());
    if (updatetau.getvalue() == false)
      t.options[12] = "false";
    else
      t.options[12] = "true";
    t.options[13] = ST::doubletostring(f.getvalue());
    if (diagtransform.getvalue() == false)
      t.options[14] = "false";
    else
      t.options[14] = "true";
    if (derivative.getvalue() == false)
      t.options[15] = "false";
    else
      t.options[15] = "true";
    t.options[16] = ST::inttostring(contourprob.getvalue());
    if (uniformprior.getvalue() == false)
      t.options[17] = "false";
    else
      t.options[17] = "true";
    t.options[18] = beta_0.getvalue();
    t.options[19] = knots.getvalue();
    if (center.getvalue() == false)
      t.options[20] = "false";
    else
      t.options[20] = "true";

    t.options[21] = centermethod.getvalue();      

    if ( contourprob.getvalue()-1 > degree.getvalue())
      {
      setdefault();
      return false;
      }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//------ class term_varcoeff_merror: implementation of member functions --------
//------------------------------------------------------------------------------

term_varcoeff_merror::term_varcoeff_merror(void)
  {
  type = "term_varcoeff_merror";
  min=intoption("min",0,1,100);
  max=intoption("max",0,1,100);
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  gridsize = intoption("gridsize",-1,10,500);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("cp");
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  vector<ST::string> adm;
  adm.push_back("unrestricted");
  adm.push_back("increasing");
  adm.push_back("decreasing");
  monotone = stroption("monotone",adm,"unrestricted");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  diagtransform = simpleoption("diagtransform",false);
  derivative = simpleoption("derivative",false);
  contourprob = intoption("contourprob",-1,0,6);
  uniformprior = simpleoption("uniformprior",false);
  beta_0 = stroption("beta_0");
  vector<ST::string> knotsdef;
  knotsdef.push_back("equidistant");
  knotsdef.push_back("quantiles");
  knots = stroption("knots",knotsdef,"equidistant");
  center = simpleoption("center",false);
  vector<ST::string> adm_centerm;
  adm_centerm.push_back("mean");
  adm_centerm.push_back("samplecentered");
  adm_centerm.push_back("meanintercept");
  centermethod = stroption("centermethod",adm_centerm,"mean");
  

  // SUSI: initialize new option
  // Syntax : doubleoption("name", default, lower limit, upper limit)
  merrorvar1 = doubleoption("merrorvar1",20,0.0,10000000);
  merrorvar2 = doubleoption("merrorvar2",20,0.0,10000000);
  arvar = doubleoption("arvar",20,0.0,10000000);
  arpar1 = doubleoption("arpar1",0.5,0.0,1);
  arpar2 = doubleoption("arpar2",0.25,0.0,1);
  biasmean = doubleoption("biasmean",0,-10000000,10000000);
  biasvar = doubleoption("biasvar",1000,0.0,10000000);
  }


void term_varcoeff_merror::setdefault(void)
  {
  min.setdefault();
  max.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  gridsize.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  monotone.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  diagtransform.setdefault();
  derivative.setdefault();
  contourprob.setdefault();
  uniformprior.setdefault();
  beta_0.setdefault();
  knots.setdefault();
  center.setdefault();
  centermethod.setdefault();

  // SUSI: call setdefault() for new option
  merrorvar1.setdefault();
  merrorvar2.setdefault();
  arvar.setdefault();
  arpar1.setdefault();
  arpar2.setdefault();
  biasmean.setdefault();
  biasvar.setdefault();
  }


bool term_varcoeff_merror::check(term & t)
  {

  if ( (t.varnames.size()==2)  && (t.options.size() >=1)
// SUSI: Adjust options.size()
        && (t.options.size() <= 29) )
    {

    if (t.options[0] == "merrorrw1")
      t.type = "varcoeffmerrorrw1";
    else if (t.options[0] == "merrorrw2")
      t.type = "varcoeffmerrorrw2";
    else
      {
      setdefault();
      return false;
      }

    optionlist optlist;
    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&gridsize);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&monotone);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&diagtransform);
    optlist.push_back(&derivative);
    optlist.push_back(&contourprob);
    optlist.push_back(&uniformprior);
    optlist.push_back(&beta_0);
    optlist.push_back(&knots);
    optlist.push_back(&center);


    // SUSI: add option to options list
    optlist.push_back(&merrorvar1);
    optlist.push_back(&merrorvar2);
    optlist.push_back(&arvar);
    optlist.push_back(&arpar1);
    optlist.push_back(&arpar2);
    optlist.push_back(&biasmean);
    optlist.push_back(&biasvar);

    optlist.push_back(&centermethod);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    // SUSI: Adjust length of options
    t.options = vector<ST::string>(29);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(min.getvalue());
    t.options[2] = ST::inttostring(max.getvalue());
    t.options[3] = ST::inttostring(degree.getvalue());
    t.options[4] = ST::inttostring(numberknots.getvalue());
    t.options[5] = ST::doubletostring(lambda.getvalue());
    t.options[6] = ST::inttostring(gridsize.getvalue());
    t.options[7] = ST::doubletostring(a.getvalue());
    t.options[8] = ST::doubletostring(b.getvalue());
    t.options[9] = proposal.getvalue();
    t.options[10] = monotone.getvalue();
    t.options[11] = ST::inttostring(updateW.getvalue());
    if (updatetau.getvalue() == false)
      t.options[12] = "false";
    else
      t.options[12] = "true";
    t.options[13] = ST::doubletostring(f.getvalue());
    if (diagtransform.getvalue() == false)
      t.options[14] = "false";
    else
      t.options[14] = "true";
    if (derivative.getvalue() == false)
      t.options[15] = "false";
    else
      t.options[15] = "true";
    t.options[16] = ST::inttostring(contourprob.getvalue());
    if (uniformprior.getvalue() == false)
      t.options[17] = "false";
    else
      t.options[17] = "true";
    t.options[18] = beta_0.getvalue();
    t.options[19] = knots.getvalue();
    if (center.getvalue() == false)
      t.options[20] = "false";
    else
      t.options[20] = "true";

    // SUSI: Add new option
    t.options[21] = ST::doubletostring(merrorvar1.getvalue());
    t.options[22] = ST::doubletostring(merrorvar2.getvalue());
    t.options[23] = ST::doubletostring(arvar.getvalue());
    t.options[24] = ST::doubletostring(arpar1.getvalue());
    t.options[25] = ST::doubletostring(arpar2.getvalue());
    t.options[26] = ST::doubletostring(biasmean.getvalue());
    t.options[27] = ST::doubletostring(biasvar.getvalue());

    t.options[28] = centermethod.getvalue();

    if ( contourprob.getvalue()-1 > degree.getvalue())
      {
      setdefault();
      return false;
      }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }



//------------------------------------------------------------------------------
//--------- class term_randomslope: implementation of member functions ---------
//------------------------------------------------------------------------------

term_randomslope::term_randomslope(void)
  {
  type = "term_randomslope";
  nofixed = simpleoption("nofixed",false);
  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  updatetau = simpleoption("updatetau",false);
  uniformprior = simpleoption("uniformprior",false);
  constlambda = simpleoption("constlambda",false);
  }

void term_randomslope::setdefault(void)
  {
  nofixed.setdefault();
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  updatetau.setdefault();
  uniformprior.setdefault();
  constlambda.setdefault();
  }


bool term_randomslope::check(term & t)
  {

  if ( (t.varnames.size()==2)  && (t.options.size()<=9) )
    {

    if (t.options[0] == "random")
      t.type = "randomslope";
    else
      {
      setdefault();
      return false;
      }

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&nofixed);
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&updatetau);
    optlist.push_back(&uniformprior);
    optlist.push_back(&constlambda);


    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(9);
    t.options[0] = t.type;
    if (nofixed.getvalue() == true)
      t.options[1] = "true";
    else
      t.options[1] = "false";
    t.options[2] = ST::doubletostring(lambda.getvalue());
    t.options[3] = ST::doubletostring(a.getvalue());
    t.options[4] = ST::doubletostring(b.getvalue());
    t.options[5] = proposal.getvalue();
    if (updatetau.getvalue() == false)
      t.options[6] = "false";
    else
      t.options[6] = "true";
    if (uniformprior.getvalue() == false)
      t.options[7] = "false";
    else
      t.options[7] = "true";
    if (constlambda.getvalue() == false)
      t.options[8] = "false";
    else
      t.options[8] = "true";

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }



//------------------------------------------------------------------------------
//----------- class term_random: implementation of member functions ------------
//------------------------------------------------------------------------------

term_random::term_random(void)
  {
  type = "term_random";
  lambda = doubleoption("lambda",100,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  updatetau = simpleoption("updatetau",false);
  uniformprior = simpleoption("uniformprior",false);
  constlambda = simpleoption("constlambda",false);
  }


void term_random::setdefault(void)
  {
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  updatetau.setdefault();
  uniformprior.setdefault();
  constlambda.setdefault();
  }


bool term_random::check(term & t)
  {

  if ( (t.varnames.size()==1)  && (t.options.size()<=8) )
    {

    if (t.options[0] == "random")
      t.type = "random";

    else
      {
      setdefault();
      return false;
      }


    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&updatetau);
    optlist.push_back(&uniformprior);
    optlist.push_back(&constlambda);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(8);
    t.options[0] = t.type;
    t.options[1] = ST::doubletostring(lambda.getvalue());
    t.options[2] = ST::doubletostring(a.getvalue());
    t.options[3] = ST::doubletostring(b.getvalue());
    t.options[4] = proposal.getvalue();
    if (updatetau.getvalue() == false)
      t.options[5] = "false";
    else
      t.options[5] = "true";
    if (uniformprior.getvalue() == false)
      t.options[6] = "false";
    else
      t.options[6] = "true";
    if (constlambda.getvalue() == false)
      t.options[7] = "false";
    else
      t.options[7] = "true";


    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }



//------------------------------------------------------------------------------
//----------- class term_hrandom: implementation of member functions ------------
//------------------------------------------------------------------------------

term_hrandom::term_hrandom(void)
  {
  type = "term_hrandom";
  lambda = doubleoption("lambda",100,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"iwls");
  updatetau = simpleoption("updatetau",false);
  uniformprior = simpleoption("uniformprior",false);
  constlambda = simpleoption("constlambda",false);
  }


void term_hrandom::setdefault(void)
  {
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  updatetau.setdefault();
  uniformprior.setdefault();
  constlambda.setdefault();
  }


bool term_hrandom::check(term & t)
  {

  if ( (t.varnames.size()==1)  && (t.options.size()<=8) )
    {

    if (t.options[0] == "hrandom")
      t.type = "hrandom";

    else
      {
      setdefault();
      return false;
      }


    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&updatetau);
    optlist.push_back(&uniformprior);
    optlist.push_back(&constlambda);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(8);
    t.options[0] = t.type;
    t.options[1] = ST::doubletostring(lambda.getvalue());
    t.options[2] = ST::doubletostring(a.getvalue());
    t.options[3] = ST::doubletostring(b.getvalue());
    t.options[4] = proposal.getvalue();
    if (updatetau.getvalue() == false)
      t.options[5] = "false";
    else
      t.options[5] = "true";
    if (uniformprior.getvalue() == false)
      t.options[6] = "false";
    else
      t.options[6] = "true";
    if (constlambda.getvalue() == false)
      t.options[7] = "false";
    else
      t.options[7] = "true";


    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }



//------------------------------------------------------------------------------
//----------- class term_mixture: implementation of member functions ------------
//------------------------------------------------------------------------------

term_mixture::term_mixture(void)
  {
  type = "term_mixture";
  nrcomp = intoption("nrcomp",1,1,50); //name,vorbelegung,minerlaubt,maxerlaubt
  wprior = doubleoption("wprior",1.0,0.0,100.0);
  mpriorm = doubleoption("mpriorm",0.0,-100.0,100.0);
  mpriorv = doubleoption("mpriorv",100,0.000001,1000);
  vpriora = doubleoption("vpriora",2.0,0.000001,100.0);
  vpriorb = doubleoption("vpriorb",1.0,0.000001,100.0);
  nosamples = simpleoption("nosamples",false);
  aclag = intoption("aclag",0,0,500);
  vector<ST::string> adm_order;
  adm_order.push_back("n");
  adm_order.push_back("w");
  order = stroption("order",adm_order,"n");
  vpriorbunif = simpleoption("vpriorbunif",false);
  vpriorbgamma = simpleoption("vpriorbgamma",false);
  }


void term_mixture::setdefault(void)
  {
  nrcomp.setdefault();
  wprior.setdefault();
  mpriorm.setdefault();
  mpriorv.setdefault();
  vpriora.setdefault();
  vpriorb.setdefault();
  nosamples.setdefault();
  aclag.setdefault();
  order.setdefault();
  vpriorbunif.setdefault();
  vpriorbgamma.setdefault();
  }


bool term_mixture::check(term & t)
  {

  if ( (t.varnames.size()==1)  && (t.options.size()<=12) ) // 12, da 12 optionen ("mixture",nrcomp,...)
    {

    if (t.options[0] == "mixture")
      t.type = "mixture";

    else
      {
      setdefault();
      return false;
      }

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&nrcomp);
    optlist.push_back(&wprior);
    optlist.push_back(&mpriorm);
    optlist.push_back(&mpriorv);
    optlist.push_back(&vpriora);
    optlist.push_back(&vpriorb);
    optlist.push_back(&nosamples);
    optlist.push_back(&aclag);
    optlist.push_back(&order);
    optlist.push_back(&vpriorbunif);
    optlist.push_back(&vpriorbgamma);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(12); // 12, s.o.
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(nrcomp.getvalue());
    t.options[2] = ST::doubletostring(wprior.getvalue());
    t.options[3] = ST::doubletostring(mpriorm.getvalue());
    t.options[4] = ST::doubletostring(mpriorv.getvalue());
    t.options[5] = ST::doubletostring(vpriora.getvalue());
    t.options[6] = ST::doubletostring(vpriorb.getvalue());
    if (nosamples.getvalue() == false)
      t.options[7] = "false";
    else
      t.options[7] = "true";
    t.options[8] = ST::inttostring(aclag.getvalue());
    t.options[9] = order.getvalue();
    if (vpriorbunif.getvalue() == false)
      t.options[10] = "false";
    else
      t.options[10] = "true";
    if (vpriorbgamma.getvalue() == false)
      t.options[11] = "false";
    else
      t.options[11] = "true";

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }



//------------------------------------------------------------------------------
//----------- class term_offset: implementation of member functions ------------
//------------------------------------------------------------------------------

bool term_offset::check(term & t)
  {

  if ( (t.varnames.size()==1)  && (t.options.size()==1) )
    {

    if (t.options[0] == "offset")
      t.type = "offset";
    else
      return false;

    return true;

    }
  else
    return false;

  }



//------------------------------------------------------------------------------
//----------- class term_baseline: implementation of member functions -----------
//------------------------------------------------------------------------------


term_baseline::term_baseline(void)
  {
  type = "term_baseline";
  min=intoption("min",0,1,100);
  max=intoption("max",0,1,100);
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.01,0,500);
  uniformb = simpleoption("uniformb",false);
  gridsize = intoption("gridsize",-1,10,500);
  uniformprior = simpleoption("uniformprior",false);
  vector<ST::string> adm_prop;
  adm_prop.push_back("cp");
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal = stroption("proposal",adm_prop,"cp");
  weibull = simpleoption("weibull",false);
  begin = stroption("begin");
  vector<ST::string> knotsdef;
  knotsdef.push_back("equidistant");
  knotsdef.push_back("quantiles");
  knots = stroption("knots",knotsdef,"equidistant");
  // NEW FOR PARTIALLIKELIHOOD
  partialLikelihood = simpleoption("partialLikelihood",false);
  }

void term_baseline::setdefault(void)
  {
  min.setdefault();
  max.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  uniformb.setdefault();
  gridsize.setdefault();
  uniformprior.setdefault();
  proposal.setdefault();
  weibull.setdefault();
  begin.setdefault();
  knots.setdefault();
  // NEW FOR PARTIALLIKELIHOOD
  partialLikelihood.setdefault();
  }

bool term_baseline::check(term & t)
  {

  if ( (t.varnames.size()==1)  && (t.options.size() >= 1)
        && (t.options.size() <= 16) )
    {

    if (t.options[0] == "baseline")
      t.type = "baseline";
    else
      {
      setdefault();
      return false;
      }

    long minim,maxim;

    optionlist optlist;
    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&uniformb);
    optlist.push_back(&gridsize);
    optlist.push_back(&uniformprior);
    optlist.push_back(&proposal);
    optlist.push_back(&weibull);
    optlist.push_back(&begin);
    optlist.push_back(&knots);
    // NEW FOR PARTIALLIKELIHOOD
    optlist.push_back(&partialLikelihood);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

   t.options.erase(t.options.begin(),t.options.end());
   t.options = vector<ST::string>(16);
   t.options[0] = t.type;
   t.options[1] = ST::inttostring(min.getvalue());
   t.options[2] = ST::inttostring(max.getvalue());
   t.options[3] = ST::inttostring(degree.getvalue());
   t.options[4] = ST::inttostring(numberknots.getvalue());
   t.options[5] = ST::doubletostring(lambda.getvalue());
   t.options[6] = ST::doubletostring(a.getvalue());
   t.options[7] = ST::doubletostring(b.getvalue());
   if (uniformb.getvalue() == false)
     t.options[8] = "false";
   else
     t.options[8] = "true";
   t.options[9] = ST::inttostring(gridsize.getvalue());
   if (uniformprior.getvalue() == false)
     t.options[10] = "false";
   else
     t.options[10] = "true";
   t.options[11] = proposal.getvalue();
   if (weibull.getvalue() == false)
     t.options[12] = "false";
   else
     t.options[12] = "true";

   t.options[13] = begin.getvalue();
   t.options[14] = knots.getvalue();
   // NEW FOR PARTIALLIKELIHOOD
   if (partialLikelihood.getvalue() == false)
     t.options[15] = "false";
   else
     t.options[15] = "true";

   if (t.options[1].strtolong(minim) == 1)
     {
     setdefault();
     return false;
     }

   if (t.options[2].strtolong(maxim) == 1)
     {
     setdefault();
     return false;
     }

   if (maxim < minim)
     {
     setdefault();
     return false;
     }

   if (lambda.getvalue() < 0)
     {
     setdefault();
     return false;
     }

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//------ class term_varcoeff_baseline: implementation of member functions -------
//------------------------------------------------------------------------------


term_varcoeff_baseline::term_varcoeff_baseline(void)
  {
  type = "term_varcoeff";
  min=intoption("min",0,1,100);
  max=intoption("max",0,1,100);
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  gridsize = intoption("gridsize",-1,10,500);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  uniformprior = simpleoption("uniformprior",false);
  vector<ST::string> knotsdef;
  knotsdef.push_back("equidistant");
  knotsdef.push_back("quantiles");
  knots = stroption("knots",knotsdef,"equidistant");
  }


void term_varcoeff_baseline::setdefault(void)
  {
  min.setdefault();
  max.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  gridsize.setdefault();
  a.setdefault();
  b.setdefault();
  uniformprior.setdefault();
  knots.setdefault();
  }


bool term_varcoeff_baseline::check(term & t)
  {

  if ( (t.varnames.size()==2)  && (t.options.size() >=1)
        && (t.options.size() <= 11) )
    {

    if (t.options[0] == "baseline")
      t.type = "varbaseline";
    else
      {
      setdefault();
      return false;
      }

    optionlist optlist;
    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&gridsize);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&uniformprior);
    optlist.push_back(&knots);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(11);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(min.getvalue());
    t.options[2] = ST::inttostring(max.getvalue());
    t.options[3] = ST::inttostring(degree.getvalue());
    t.options[4] = ST::inttostring(numberknots.getvalue());
    t.options[5] = ST::doubletostring(lambda.getvalue());
    t.options[6] = ST::inttostring(gridsize.getvalue());
    t.options[7] = ST::doubletostring(a.getvalue());
    t.options[8] = ST::doubletostring(b.getvalue());
    if(uniformprior.getvalue() == false)
      t.options[9] = "false";
    else
      t.options[9] = "true";
    t.options[10] = knots.getvalue();

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }




//------------------------------------------------------------------------------
//------------- CLASS modelterm: implementation of member functions ------------
//------------------------------------------------------------------------------


const modelterm & modelterm::operator=(const modelterm & m)
  {
  if (this==&m)
    return *this;
  model::operator=(model(m));
  terms = m.terms;
  termtypes=m.termtypes;
  responsevar = m.responsevar;
  return *this;
  }


void modelterm::parse(const ST::string & m)
  {

  terms.erase(terms.begin(),terms.end());

  model::parse(m);
  ST::string mod;
  mod = m.eatallwhitespace();
  mod = mod.eatallcarriagereturns();

  bool bracketmiss=false;
  vector<ST::string> token = mod.strtoken2("=",bracketmiss);
  if (bracketmiss==true)
    errormessages.push_back("ERROR: missing brackets\n");
  else if (token.size() > 2 || token.size() < 1)
    errormessages.push_back("ERROR: invalid model specification\n");
  else
    {
    if (token[0].isvarname()== 1)
      errormessages.push_back("ERROR: " + token[0] + " invalid variable name\n");
    else
      {
      responsevar = token[0];
      modelVarnames.push_back(token[0]);
      }
    }


  if (errormessages.empty() && token.size()==2)
    {
    token = token[1].strtoken("+",false);
    terms = vector<term>(token.size());
    unsigned i,k;
    bool found;
    for(i=0;i<token.size();i++)
      {

      terms[i].parse(token[i]);
      if (terms[i].geterrormessages().size() > 0)
        errormessages.insert(errormessages.end(),
        terms[i].geterrormessages().begin(),
        terms[i].geterrormessages().end() );
      else
        {
        k=0;
        found = false;
        while ( (k<termtypes->size()) && (!found) )
          {
          found = (*termtypes)[k]->check(terms[i]);
          if (found)
            {
            unsigned j;
            for (j=0;j<terms[i].varnames.size();j++)
              modelVarnames.push_back(terms[i].varnames[j]);
            }
          k++;
          }
        if (!found)
          errormessages.push_back("ERROR: syntax error in term " + token[i] + "\n");
        }

      } // end: for(i=0;i<token.size();i++)

    }

  if (errormessages.empty())
    {

    modeltext=m;
    modelexisting = true;
    }
  else
    {
    model::clear();
    terms.erase(terms.begin(),terms.end());
    }

  }




//------------------------------------------------------------------------------
//------------- CLASS modeltermmult: implementation of member functions --------
//------------------------------------------------------------------------------

void modeltermmult::clear(void)
  {
  model::clear();
  responsevar.erase(responsevar.begin(),responsevar.end());
  responsecol.erase(responsecol.begin(),responsecol.end()); 
  }


const modeltermmult & modeltermmult::operator=(const modeltermmult & m)
  {
  if (this==&m)
    return *this;
  model::operator=(model(m));
  terms = m.terms;
  termtypes=m.termtypes;
  responsevar = m.responsevar;
  responsecol = m.responsecol;
  return *this;
  }


void modeltermmult::parse(const ST::string & m)
  {

  model::parse(m);
  ST::string mod;
  mod = m.eatallwhitespace();
  mod = mod.eatallcarriagereturns();

  vector<ST::string> equations = mod.strtoken(":",false);

  terms = vector< vector<term> >(equations.size());

  unsigned e;
  for(e=0;e<equations.size();e++)
    {

    vector<ST::string> token;

    if (errormessages.empty())
      {
      bool bracketmiss=false;
      token = equations[e].strtoken2("=",bracketmiss);
      if (bracketmiss==true)
        errormessages.push_back("ERROR: missing brackets\n");
      else if (token.size() != 2)
        errormessages.push_back("ERROR: invalid model specification\n");
      else
        {
        if (token[0].isvarname()== 1)
          errormessages.push_back("ERROR: " + token[0] + " invalid variable name\n");
        else
          {
          responsevar.push_back(token[0]);
          modelVarnames.push_back(token[0]);
          responsecol.push_back(modelVarnames.size()-1);
          }
        }
      }


    if (errormessages.empty())
      {
      token = token[1].strtoken("+",false);
      terms[e] = vector<term>(token.size());
      unsigned i,k;
      bool found;
      for(i=0;i<token.size();i++)
        {

        terms[e][i].parse(token[i]);
        if (terms[e][i].geterrormessages().size() > 0)
          errormessages.insert(errormessages.end(),
          terms[e][i].geterrormessages().begin(),
          terms[e][i].geterrormessages().end() );
        else
          {
          k=0;
          found = false;
          while ( (k<termtypes->size()) && (!found) )
            {
            found = (*termtypes)[k]->check(terms[e][i]);
            if (found)
              {
              unsigned j;
              for (j=0;j<terms[e][i].varnames.size();j++)
                modelVarnames.push_back(terms[e][i].varnames[j]);
              }
          k++;
          }
          if (!found)
            errormessages.push_back("ERROR: syntax error in term " + token[i] + "\n");
          }

        } // end: for(i=0;i<token.size();i++)

      }

    } // end: for(e=0;e<equations.size();e++)


  if (errormessages.empty())
    {
    modeltext=m;
    modelexisting = true;
    }
  else
    {
    model::clear();
    terms.erase(terms.begin(),terms.end());
    }


  }




//------------------------------------------------------------------------------
//--------- class term_random_autoreg: implementation of member functions ------
//------------------------------------------------------------------------------


term_random_autoreg::term_random_autoreg(void)
  {
  type = "term_random_autoreg";

  lambda_r = doubleoption("lambda_r",100000,0,10000000);
  a_r = doubleoption("a_r",0.001,-1.0,500);
  b_r = doubleoption("b_r",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal_r = stroption("proposal_r",adm_prop,"iwls");
  updatetau_r = simpleoption("updatetau_r",false);
  uniformprior_r = simpleoption("uniformprior_r",false);
  constlambda_r = simpleoption("constlambda_r",false);

  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  proposal = stroption("proposal",adm_prop,"iwls");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,10000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,10000000);
  lambdastart = doubleoption("lambdastart",-1,-1,10000000);
  stationary = simpleoption("stationary",false);
  alpha = doubleoption("alpha",0.9,-1.0,1.0);
  alphafix = simpleoption("alphafix",false);

  }


void term_random_autoreg::setdefault(void)
  {
  lambda_r.setdefault();
  a_r.setdefault();
  b_r.setdefault();
  proposal_r.setdefault();
  updatetau_r.setdefault();
  uniformprior_r.setdefault();
  constlambda_r.setdefault();  // 7

  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  stationary.setdefault();
  alpha.setdefault();
  alphafix.setdefault();     // 13 (total 13+7=20)
  }


bool term_random_autoreg::check(term & t)
  {

  if ( (t.varnames.size()==2)  && (t.options.size()<=21) )
    {

    if (t.options[0] == "random_rw1")
      t.type = "random_rw1";
    else if   (t.options[0] == "random_rw2")
      t.type = "random_rw2";
    else
      {
      setdefault();
      return false;
      }

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&lambda_r);
    optlist.push_back(&a_r);
    optlist.push_back(&b_r);
    optlist.push_back(&proposal_r);
    optlist.push_back(&updatetau_r);
    optlist.push_back(&uniformprior_r);
    optlist.push_back(&constlambda_r);

    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&stationary);
    optlist.push_back(&alpha);
    optlist.push_back(&alphafix);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(21);
    t.options[0] = t.type;
    t.options[1] = ST::doubletostring(lambda_r.getvalue());
    t.options[2] = ST::doubletostring(a_r.getvalue());
    t.options[3] = ST::doubletostring(b_r.getvalue());
    t.options[4] = proposal_r.getvalue();
    if (updatetau_r.getvalue() == false)
      t.options[5] = "false";
    else
      t.options[5] = "true";
    if (uniformprior_r.getvalue() == false)
      t.options[6] = "false";
    else
      t.options[6] = "true";
    if (constlambda_r.getvalue() == false)
      t.options[7] = "false";
    else
      t.options[7] = "true";

    t.options[8] = ST::doubletostring(lambda.getvalue());
    t.options[9] = ST::doubletostring(a.getvalue());
    t.options[10] = ST::doubletostring(b.getvalue());
    t.options[11] = proposal.getvalue();
    t.options[12] = ST::inttostring(updateW.getvalue());
    if (updatetau.getvalue()==false)
      t.options[13] = "false";
    else
      t.options[13] = "true";
    t.options[14] = ST::doubletostring(f.getvalue());
    t.options[15] = ST::doubletostring(lambdamin.getvalue());
    t.options[16] = ST::doubletostring(lambdamax.getvalue());
    t.options[17] = ST::doubletostring(lambdastart.getvalue());
    if(stationary.getvalue() == false)
      t.options[18] = "false";
    else
      t.options[18] = "true";
    t.options[19] = ST::doubletostring(alpha.getvalue());
    if(alphafix.getvalue() == false)
      t.options[20] = "false";
    else
      t.options[20] = "true";


    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//--------- class term_spatial_autoreg: implementation of member functions -----
//------------------------------------------------------------------------------

term_spatial_autoreg::term_spatial_autoreg(void)
  {
  type = "term_spatial_autoreg";

  map_s=stroption("map_s");
  lambda_s = doubleoption("lambda_s",0.1,0,10000000);
  a_s = doubleoption("a_s",0.001,-1.0,500);
  b_s = doubleoption("b_s",0.001,0,500);
  vector<ST::string> adm_prop_s;
  adm_prop_s.push_back("cp");
  adm_prop_s.push_back("iwls");
  adm_prop_s.push_back("iwlsmode");
  proposal_s = stroption("proposal_s",adm_prop_s,"iwls");
  updateW_s = intoption("updateW_s",1,0,100);
  updatetau_s = simpleoption("updatetau_s",false);
  f_s = doubleoption("f_s",2,0,10000000);
  uniformprior_s = simpleoption("uniformprior_s",false);
  nrrows_s = intoption("nrrows_s",2,0,100);
  Laplace_s = simpleoption("Laplace_s",false);
  stationary_s = simpleoption("stationary_s",false);
  alpha_s = doubleoption("alpha_s",0.9,-1.0,1.0);
  alphafix_s = simpleoption("alphafix_s",false);
  center_s = simpleoption("center_s",false);

  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  proposal = stroption("proposal",adm_prop_s,"iwls");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  lambdamin = doubleoption("lambdamin",0.0001,0.000001,10000000);
  lambdamax = doubleoption("lambdamax",10000,0.000001,10000000);
  lambdastart = doubleoption("lambdastart",-1,-1,10000000);
  stationary = simpleoption("stationary",false);
  alpha = doubleoption("alpha",0.9,-1.0,1.0);
  alphafix = simpleoption("alphafix",false);

  }


void term_spatial_autoreg::setdefault(void)
  {

  map_s.setdefault();
  lambda_s.setdefault();
  a_s.setdefault();
  b_s.setdefault();
  proposal_s.setdefault();
  updateW_s.setdefault();
  updatetau_s.setdefault();
  f_s.setdefault();
  uniformprior_s.setdefault();
  nrrows_s.setdefault();
  Laplace_s.setdefault();
  stationary_s.setdefault();
  alpha_s.setdefault();
  alphafix_s.setdefault();
  center_s.setdefault();   // 15

  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  proposal.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  lambdamin.setdefault();
  lambdamax.setdefault();
  lambdastart.setdefault();
  stationary.setdefault();
  alpha.setdefault();
  alphafix.setdefault();     // 13 (total 15+13=28)
  }


bool term_spatial_autoreg::check(term & t)
  {

  if ( (t.varnames.size()==2)  && (t.options.size()<=29) )
    {

    if (t.options[0] == "spatial_rw1")
      t.type = "spatial_rw1";
    else if   (t.options[0] == "spatial_rw2")
      t.type = "spatial_rw2";
    else
      {
      setdefault();
      return false;
      }

    vector<ST::string> opt;
    optionlist optlist;

    optlist.push_back(&map_s);
    optlist.push_back(&lambda_s);
    optlist.push_back(&a_s);
    optlist.push_back(&b_s);
    optlist.push_back(&proposal_s);
    optlist.push_back(&updateW_s);
    optlist.push_back(&updatetau_s);
    optlist.push_back(&f_s);
    optlist.push_back(&uniformprior_s);
    optlist.push_back(&nrrows_s);
    optlist.push_back(&Laplace_s);
    optlist.push_back(&stationary_s);
    optlist.push_back(&alpha_s);
    optlist.push_back(&alphafix_s);
    optlist.push_back(&center_s);

    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&proposal);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&lambdamin);
    optlist.push_back(&lambdamax);
    optlist.push_back(&lambdastart);
    optlist.push_back(&stationary);
    optlist.push_back(&alpha);
    optlist.push_back(&alphafix);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(29);
    t.options[0] = t.type;

    t.options[1] = map_s.getvalue();
    t.options[2] = ST::doubletostring(lambda_s.getvalue());
    t.options[3] = ST::doubletostring(a_s.getvalue());
    t.options[4] = ST::doubletostring(b_s.getvalue());
    t.options[5] = proposal_s.getvalue();
    t.options[6] = ST::inttostring(updateW_s.getvalue());
    if (updatetau_s.getvalue()==false)
      t.options[7] = "false";
    else
      t.options[7] = "true";
    t.options[8] = ST::doubletostring(f_s.getvalue());
    if (uniformprior_s.getvalue() == false)
      t.options[9] = "false";
    else
      t.options[9] = "true";
    t.options[10] = ST::inttostring(nrrows_s.getvalue());
    if (Laplace_s.getvalue()==false)
      t.options[11] = "false";
    else
      t.options[11] = "true";
    if(stationary_s.getvalue() == false)
      t.options[12] = "false";
    else
      t.options[12] = "true";
    t.options[13] = ST::doubletostring(alpha_s.getvalue());
    if(alphafix_s.getvalue() == false)
      t.options[14] = "false";
    else
      t.options[14] = "true";
    if(center_s.getvalue() == false)
      t.options[15] = "false";
    else
      t.options[15] = "true";

    t.options[16] = ST::doubletostring(lambda.getvalue());
    t.options[17] = ST::doubletostring(a.getvalue());
    t.options[18] = ST::doubletostring(b.getvalue());
    t.options[19] = proposal.getvalue();
    t.options[20] = ST::inttostring(updateW.getvalue());
    if (updatetau.getvalue()==false)
      t.options[21] = "false";
    else
      t.options[21] = "true";
    t.options[22] = ST::doubletostring(f.getvalue());
    t.options[23] = ST::doubletostring(lambdamin.getvalue());
    t.options[24] = ST::doubletostring(lambdamax.getvalue());
    t.options[25] = ST::doubletostring(lambdastart.getvalue());
    if(stationary.getvalue() == false)
      t.options[26] = "false";
    else
      t.options[26] = "true";
    t.options[27] = ST::doubletostring(alpha.getvalue());
    if(alphafix.getvalue() == false)
      t.options[28] = "false";
    else
      t.options[28] = "true";

    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }


//------------------------------------------------------------------------------
//-------- class term_random_pspline: implementation of member functions -------
//------------------------------------------------------------------------------

term_random_pspline::term_random_pspline(void)
  {
  type = "term_random_pspline";

  lambda_r = doubleoption("lambda_r",100000,0,10000000);
  a_r = doubleoption("a_r",0.001,-1.0,500);
  b_r = doubleoption("b_r",0.001,0,500);
  vector<ST::string> adm_prop;
  adm_prop.push_back("iwls");
  adm_prop.push_back("iwlsmode");
  proposal_r = stroption("proposal_r",adm_prop,"iwls");
  updatetau_r = simpleoption("updatetau_r",false);
  uniformprior_r = simpleoption("uniformprior_r",false);
  constlambda_r = simpleoption("constlambda_r",false);


  min=intoption("min",0,1,100);
  max=intoption("max",0,1,100);
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  a = doubleoption("a",0.001,-1.0,500);
  b = doubleoption("b",0.001,0,500);
  uniformb = simpleoption("uniformb",false);
  gridsize = intoption("gridsize",-1,10,500);
  minvar=intoption("minvar",1,1,500);
  maxvar=intoption("maxvar",1,1,500);
  startv = doubleoption("startv",0.05,0.00001,1000);
  proposal = stroption("proposal",adm_prop,"iwls");
  vector<ST::string> adm;
  adm.push_back("unrestricted");
  adm.push_back("increasing");
  adm.push_back("decreasing");
  monotone = stroption("monotone",adm,"unrestricted");
  updateW = intoption("updateW",1,0,100);
  updatetau = simpleoption("updatetau",false);
  f = doubleoption("f",2,0,10000000);
  diagtransform = simpleoption("diagtransform",false);
  derivative = simpleoption("derivative",false);
  bsplinebasis = simpleoption("bsplinebasis",false);
  contourprob = intoption("contourprob",-1,0,6);
  uniformprior = simpleoption("uniformprior",false);
  beta_0 = stroption("beta_0");
  discrete = simpleoption("discrete",false);
  df = intoption("df",20,3,50);
  stationary = simpleoption("stationary",false);
  alpha = doubleoption("alpha",0.9,-1.0,1.0);
  alphafix = simpleoption("alphafix",false);
  vector<ST::string> knotsdef;
  knotsdef.push_back("equidistant");
  knotsdef.push_back("quantiles");
  knots = stroption("knots",knotsdef,"equidistant");
  lowerknot = doubleoption("lowerknot",0,-10000000,10000000);
  upperknot = doubleoption("upperknot",0,-10000000,10000000);
  }

void term_random_pspline::setdefault(void)
  {

  lambda_r.setdefault();
  a_r.setdefault();
  b_r.setdefault();
  proposal_r.setdefault();
  updatetau_r.setdefault();
  uniformprior_r.setdefault();
  constlambda_r.setdefault();  // 7

  min.setdefault();
  max.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  a.setdefault();
  b.setdefault();
  uniformb.setdefault();
  gridsize.setdefault();
  minvar.setdefault();
  maxvar.setdefault();
  startv.setdefault();
  proposal.setdefault();
  monotone.setdefault();
  updateW.setdefault();
  updatetau.setdefault();
  f.setdefault();
  diagtransform.setdefault();
  derivative.setdefault();
  bsplinebasis.setdefault();
  contourprob.setdefault();
  uniformprior.setdefault();
  beta_0.setdefault();
  discrete.setdefault();
  df.setdefault();
  alpha.setdefault();
  stationary.setdefault();
  alphafix.setdefault();
  knots.setdefault();
  lowerknot.setdefault();
  upperknot.setdefault();   // 31 (31+7=38)
  }

bool term_random_pspline::check(term & t)
  {

  if ( (t.varnames.size()==2)  && (t.options.size() >= 1)
        && (t.options.size() <= 39) )
    {

    if (t.options[0] == "random_psplinerw1")
      t.type = "random_psplinerw1";
    else if (t.options[0] == "random_psplinerw2")
      t.type = "random_psplinerw2";
    else if (t.options[0] == "random_tpsplinerw1")
      t.type = "random_tpsplinerw1";
    else if (t.options[0] == "random_tpsplinerw2")
      t.type = "random_tpsplinerw2";
    else if (t.options[0] == "random_psplinerw1vrw1")
      t.type = "random_psplinerw1vrw1";
    else if (t.options[0] == "random_psplinerw1vrw2")
      t.type = "random_psplinerw1vrw2";
    else if (t.options[0] == "random_psplinerw2vrw1")
      t.type = "random_psplinerw2vrw1";
    else if (t.options[0] == "random_psplinerw2vrw2")
      t.type = "random_psplinerw2vrw2";
    else
      {
      setdefault();
      return false;
      }

    long minim,maxim;

    optionlist optlist;

    optlist.push_back(&lambda_r);
    optlist.push_back(&a_r);
    optlist.push_back(&b_r);
    optlist.push_back(&proposal_r);
    optlist.push_back(&updatetau_r);
    optlist.push_back(&uniformprior_r);
    optlist.push_back(&constlambda_r);

    optlist.push_back(&min);
    optlist.push_back(&max);
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&a);
    optlist.push_back(&b);
    optlist.push_back(&uniformb);
    optlist.push_back(&gridsize);
    optlist.push_back(&minvar);
    optlist.push_back(&maxvar);
    optlist.push_back(&startv);
    optlist.push_back(&proposal);
    optlist.push_back(&monotone);
    optlist.push_back(&updateW);
    optlist.push_back(&updatetau);
    optlist.push_back(&f);
    optlist.push_back(&diagtransform);
    optlist.push_back(&derivative);
    optlist.push_back(&bsplinebasis);
    optlist.push_back(&contourprob);
    optlist.push_back(&uniformprior);
    optlist.push_back(&beta_0);
    optlist.push_back(&discrete);
    optlist.push_back(&df);
    optlist.push_back(&stationary);
    optlist.push_back(&alpha);
    optlist.push_back(&alphafix);
    optlist.push_back(&knots);
    optlist.push_back(&lowerknot);
    optlist.push_back(&upperknot);

    unsigned i;
    bool rec = true;
    for (i=1;i<t.options.size();i++)
      {

      if (optlist.parse(t.options[i],true) == 0)
        rec = false;

      if (optlist.geterrormessages().size() > 0)
        {
        setdefault();
        return false;
        }

      }

    if (rec == false)
      {
      setdefault();
      return false;
      }

   t.options.erase(t.options.begin(),t.options.end());
   t.options = vector<ST::string>(39);
   t.options[0] = t.type;

    t.options[1] = ST::doubletostring(lambda_r.getvalue());
    t.options[2] = ST::doubletostring(a_r.getvalue());
    t.options[3] = ST::doubletostring(b_r.getvalue());
    t.options[4] = proposal_r.getvalue();
    if (updatetau_r.getvalue() == false)
      t.options[5] = "false";
    else
      t.options[5] = "true";
    if (uniformprior_r.getvalue() == false)
      t.options[6] = "false";
    else
      t.options[6] = "true";
    if (constlambda_r.getvalue() == false)
      t.options[7] = "false";
    else
      t.options[7] = "true";

   t.options[8] = ST::inttostring(min.getvalue());
   t.options[9] = ST::inttostring(max.getvalue());
   t.options[10] = ST::inttostring(degree.getvalue());
   t.options[11] = ST::inttostring(numberknots.getvalue());
   t.options[12] = ST::doubletostring(lambda.getvalue());
   t.options[13] = ST::doubletostring(a.getvalue());
   t.options[14] = ST::doubletostring(b.getvalue());
   if (uniformb.getvalue() == false)
     t.options[15] = "false";
   else
     t.options[15] = "true";
   t.options[16] = ST::inttostring(gridsize.getvalue());
   t.options[17] = ST::inttostring(minvar.getvalue());
   t.options[18] = ST::inttostring(maxvar.getvalue());
   t.options[19] = ST::doubletostring(startv.getvalue());
   t.options[20] = proposal.getvalue();
   t.options[21] = monotone.getvalue();
   t.options[22] = ST::inttostring(updateW.getvalue());
   if (updatetau.getvalue() == false)
     t.options[23] = "false";
   else
     t.options[23] = "true";
   t.options[24] = ST::doubletostring(f.getvalue());
   if (diagtransform.getvalue() == false)
     t.options[25] = "false";
   else
     t.options[25] = "true";
   if (derivative.getvalue() == false)
     t.options[26] = "false";
   else
     t.options[26] = "true";
   if (bsplinebasis.getvalue() == false)
     t.options[27] = "false";
   else
     t.options[27] = "true";
   t.options[28] = ST::inttostring(contourprob.getvalue());
   if(uniformprior.getvalue() == false)
     t.options[29] = "false";
   else
     t.options[29] = "true";
   t.options[30] = beta_0.getvalue();
   if(discrete.getvalue() == false)
     t.options[31] = "false";
   else
     t.options[31] = "true";
   t.options[32] = ST::inttostring(df.getvalue());
   if(stationary.getvalue() == false)
     t.options[33] = "false";
   else
     t.options[33] = "true";
   t.options[34] = ST::doubletostring(alpha.getvalue());
   if(alphafix.getvalue() == false)
     t.options[35] = "false";
   else
     t.options[35] = "true";
   t.options[36] = knots.getvalue();
   t.options[37] = ST::doubletostring(lowerknot.getvalue());
   t.options[38] = ST::doubletostring(upperknot.getvalue());

   if (t.options[8].strtolong(minim) == 1)
     {
     setdefault();
     return false;
     }

   if (t.options[9].strtolong(maxim) == 1)
     {
     setdefault();
     return false;
     }

   if (maxim < minim)
     {
     setdefault();
     return false;
     }

   if (lambda.getvalue() < 0)
     {
     setdefault();
     return false;
     }

   if ( contourprob.getvalue()-1 > degree.getvalue())
     {
     setdefault();
     return false;
     }


    setdefault();
    return true;

    }
  else
    {
    setdefault();
    return false;
    }

  }



