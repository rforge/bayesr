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





#include "model_remlreg.h"

//------------------------------------------------------------------------------
//----- class term_fixed_catspecific: implementation of member functions -------
//------------------------------------------------------------------------------

term_fixed_catspecific::term_fixed_catspecific(void)
  {
  type = "term_fixed_cat";
  }

bool term_fixed_catspecific::check(term & t)
  {
  if ( (t.varnames.size() == 1) && t.options.size()==1)
    {
    if (t.options[0] == "catspecific")
      t.type = "catspecific";
    else
      {
      setdefault();
      return false;
      }

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(1);
    t.options[0] = t.type;

    setdefault();
    return true;
    }
  else
    {
    setdefault();
    return false;
    }
  }

bool term_fixed_catspecific::checkvector(const vector<term> & terms,
                                       const unsigned & i)
  {
  assert(i< terms.size());
  if (terms[i].type == "catspecific")
    return true;
  return false;
  }

//------------------------------------------------------------------------------
//----- class term_autoreg_remlreg: implementation of member functions ---------
//------------------------------------------------------------------------------

term_autoreg_remlreg::term_autoreg_remlreg(void)
  {
  type = "term_autoreg";
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  }

void term_autoreg_remlreg::setdefault(void)
  {
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  }

bool term_autoreg_remlreg::check(term & t)
  {
  if ( (t.varnames.size() == 1) &&
       (t.options.size()<=4) && (t.options.size() >= 1) )
    {
    if (t.options[0] == "rw1")
      t.type = "rw1";
    else if (t.options[0] == "rw2")
      t.type = "rw2";
    else
      {
      setdefault();
      return false;
      }

    double startl;

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);

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

    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(4);
    t.options[0] = t.type;
    t.options[1] = ST::doubletostring(lambda.getvalue());
    t.options[2] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue()==true)
      {
      t.options[3] = "true";
      }
    else
      {
      t.options[3] = "false";
      }


    int b = t.options[2].strtodouble(startl);
    if (b==1)
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

bool term_autoreg_remlreg::checkvector(const vector<term> & terms,
                                       const unsigned & i)
  {
  assert(i< terms.size());
  if ( (terms[i].type == "rw1") || (terms[i].type == "rw2"))
    return true;
  return false;
  }

//------------------------------------------------------------------------------
//-- class term_autoreg_varcoef_remlreg: implementation of member functions ----
//------------------------------------------------------------------------------

term_autoreg_varcoef_remlreg::term_autoreg_varcoef_remlreg(void)
  {
  type = "term_autoreg_varcoef";
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  center = simpleoption("center",false);
  }

void term_autoreg_varcoef_remlreg::setdefault(void)
  {
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  center.setdefault();
  }

bool term_autoreg_varcoef_remlreg::check(term & t)
  {
  if ( (t.varnames.size() == 2) &&
       (t.options.size()<=5) && (t.options.size() >= 1) )
    {

    if (t.options[0] == "rw1")
      t.type = "varcoeffrw1";
    else if (t.options[0] == "rw2")
      t.type = "varcoeffrw2";
    else
      {
      setdefault();
      return false;
      }

    double startl;

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);
    optlist.push_back(&center);

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

    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);

    t.options.erase(t.options.begin(),t.options.end());
    t.options = vector<ST::string>(5);
    t.options[0] = t.type;
    t.options[1] = ST::doubletostring(lambda.getvalue());
    t.options[2] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      t.options[3] = "true";
    else
      t.options[3] = "false";
    if(center.getvalue())
      t.options[4] = "true";
    else
      t.options[4] = "false";

    int b = t.options[2].strtodouble(startl);
    if (b==1)
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

bool term_autoreg_varcoef_remlreg::checkvector(const vector<term> & terms,
                                               const unsigned & i)
  {
  assert(i< terms.size());

  if ((terms[i].type == "varcoeffrw1") || (terms[i].type == "varcoeffrw2"))
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//------- class term_season_remlreg: implementation of member functions --------
//------------------------------------------------------------------------------

term_season_remlreg::term_season_remlreg(void)
  {
  type = "term_season";
  period = intoption("period",12,2,72);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  }

void term_season_remlreg::setdefault(void)
  {
  period.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  }

bool term_season_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==1)  && (t.options.size()<=5) &&
       (t.options.size() >= 1) )
    {

    if (t.options[0] == "season")
      t.type = "season";
    else
      {
      setdefault();
      return false;
      }

    long per;

    vector<ST::string> opt;
    optionlist optlist;

    optlist.push_back(&period);
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);

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
    t.options = vector<ST::string>(5);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(period.getvalue());
    t.options[2] = ST::doubletostring(lambda.getvalue());
    t.options[3] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options[4] = "true";
      }
    else
      {
      t.options[4] = "false";
      }

    if (t.options[1].strtolong(per) == 1)
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

bool term_season_remlreg::checkvector(const vector<term> & terms,
                                      const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "season")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//-- class term_season_varcoef_remlreg: implementation of member functions -----
//------------------------------------------------------------------------------

term_season_varcoef_remlreg::term_season_varcoef_remlreg(void)
  {
  type = "term_season_varcoef";
  period = intoption("period",12,2,72);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  }

void term_season_varcoef_remlreg::setdefault(void)
  {
  period.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  }

bool term_season_varcoef_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==2)  && (t.options.size()<=5) &&
       (t.options.size() >= 1) )
    {

    if( (t.options[0] == "season") && (t.varnames.size() == 2) )
      t.type = "varcoeffseason";
    else
      {
      setdefault();
      return false;
      }

    long per;

    vector<ST::string> opt;
    optionlist optlist;

    optlist.push_back(&period);
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);

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
    t.options = vector<ST::string>(5);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(period.getvalue());
    t.options[2] = ST::doubletostring(lambda.getvalue());
    t.options[3] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options[4] = "true";
      }
    else
      {
      t.options[4] = "false";
      }

    if (t.options[1].strtolong(per) == 1)
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

bool term_season_varcoef_remlreg::checkvector(const vector<term> & terms,
                                              const unsigned & i)
  {
  assert(i< terms.size());

  if ( terms[i].type == "varcoeffseason")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//-------- class term_pspline_remlreg: implementation of member functions ------
//------------------------------------------------------------------------------

term_pspline_remlreg::term_pspline_remlreg(void)
  {
  type = "term_pspline";
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  gridsize = intoption("gridsize",-1,10,500);
  diagtransform = simpleoption("diagtransform",false);
  derivative = simpleoption("derivative",false);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  lowergrid = doubleoption("lowergrid",0,-10000000,10000000);
  uppergrid = doubleoption("uppergrid",0,-10000000,10000000);
  lowerknot = doubleoption("lowerknot",0,-10000000,10000000);
  upperknot = doubleoption("upperknot",0,-10000000,10000000);
  reference = doubleoption("reference",-9999,-10000000,10000000);
  }

void term_pspline_remlreg::setdefault(void)
  {
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  gridsize.setdefault();
  diagtransform.setdefault();
  derivative.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  lowergrid.setdefault();
  uppergrid.setdefault();
  lowerknot.setdefault();
  upperknot.setdefault();
  reference.setdefault();
  }

bool term_pspline_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==1)  && (t.options.size() >= 1)
        && (t.options.size() <= 14) )
    {
    if (t.options[0] == "psplinerw1")
      t.type = "psplinerw1";
    else if (t.options[0] == "psplinerw2")
      t.type = "psplinerw2";
    else
      {
      setdefault();
      return false;
      }

    double startl;

    optionlist optlist;
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&gridsize);
    optlist.push_back(&diagtransform);
    optlist.push_back(&derivative);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);
    optlist.push_back(&lowergrid);
    optlist.push_back(&uppergrid);
    optlist.push_back(&lowerknot);
    optlist.push_back(&upperknot);
    optlist.push_back(&reference);

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
   t.options = vector<ST::string>(14);
   t.options[0] = t.type;
   t.options[1] = ST::inttostring(degree.getvalue());
   t.options[2] = ST::inttostring(numberknots.getvalue());
   t.options[3] = ST::doubletostring(lambda.getvalue());
   t.options[4] = ST::inttostring(gridsize.getvalue());
   if (diagtransform.getvalue() == false)
     t.options[5] = "false";
   else
     t.options[5] = "true";
   if (derivative.getvalue() == false)
     t.options[6] = "false";
   else
     t.options[6] = "true";
   t.options[7] = ST::doubletostring(lambdastart.getvalue());

   if (lambda.getvalue() < 0)
     {
     setdefault();
     return false;
     }

    int b = t.options[7].strtodouble(startl);

    if (b==1)
      {
      setdefault();
      return false;
      }
    if(catspecific.getvalue())
      {
      t.options[8] = "true";
      }
    else
      {
      t.options[8] = "false";
      }
    t.options[9] = ST::doubletostring(lowergrid.getvalue());
    t.options[10] = ST::doubletostring(uppergrid.getvalue());
    t.options[11] = ST::doubletostring(lowerknot.getvalue());
    t.options[12] = ST::doubletostring(upperknot.getvalue());
    t.options[13] = ST::doubletostring(reference.getvalue());

    setdefault();
    return true;
    }
  else
    {
    setdefault();
    return false;
    }
  }

bool term_pspline_remlreg::checkvector(const vector<term> & terms,
                                       const unsigned & i)
  {
  assert(i< terms.size());

  if ( (terms[i].type == "psplinerw1") || (terms[i].type == "psplinerw2"))
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//-- class term_varcoeff_pspline_remlreg: implementation of member functions ---
//------------------------------------------------------------------------------

term_varcoeff_pspline_remlreg::term_varcoeff_pspline_remlreg(void)
  {
  type = "term_varcoeff";
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",0.1,0,10000000);
  catspecific = simpleoption("catspecific",false);
  center = simpleoption("center",false);
  reference = doubleoption("reference",-9999,-10000000,10000000);
  }

void term_varcoeff_pspline_remlreg::setdefault(void)
  {
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  center.setdefault();
  reference.setdefault();
  }

bool term_varcoeff_pspline_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==2)  && (t.options.size() >=1)
        && (t.options.size() <= 8) )
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
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);
    optlist.push_back(&center);
    optlist.push_back(&reference);

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
    t.options[1] = ST::inttostring(degree.getvalue());
    t.options[2] = ST::inttostring(numberknots.getvalue());
    t.options[3] = ST::doubletostring(lambda.getvalue());
    t.options[4] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options[5] = "true";
      }
    else
      {
      t.options[5] = "false";
      }
    if(center.getvalue())
      {
      t.options[6] = "true";
      }
    else
      {
      t.options[6] = "false";
      }
    t.options[7] = ST::doubletostring(reference.getvalue());

    setdefault();
    return true;
    }
  else
    {
    setdefault();
    return false;
    }
  }

bool term_varcoeff_pspline_remlreg::checkvector(const vector<term> & terms,
                                                const unsigned & i)
  {
  assert(i< terms.size());

  if ((terms[i].type == "varpsplinerw1") || (terms[i].type == "varpsplinerw2"))
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//------- class term_baseline_remlreg: implementation of member functions ------
//------------------------------------------------------------------------------

term_baseline_remlreg::term_baseline_remlreg(void)
  {
  type = "term_baseline";
  degree=intoption("degree",3,0,5);
  numberknots=intoption("nrknots",20,5,500);
  tgrid = intoption("tgrid",300,5,10000);
  vector<ST::string> gridchoices;
  gridchoices.push_back("equidistant");
  gridchoices.push_back("quantiles");
  gridchoices.push_back("all");
  gridchoice = stroption("gridchoice",gridchoices,"quantiles");
  numberquantiles = intoption("nrquantiles",60,5,1000);
  numberbetween = intoption("nrbetween",5,1,100);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",1000,0,10000000);
  lower = stroption("lower");
  catspecific = simpleoption("catspecific",false);
  gridsize = intoption("gridsize",-1,10,500);
  reference = doubleoption("reference",-9999,-10000000,10000000);
  }

void term_baseline_remlreg::setdefault(void)
  {
  degree.setdefault();
  numberknots.setdefault();
  tgrid.setdefault();
  gridchoice.setdefault();
  numberquantiles.setdefault();
  numberbetween.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  gridsize.setdefault();
  reference.setdefault();
  }

bool term_baseline_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==1)  && (t.options.size() >= 1)
        && (t.options.size() <= 13) )
    {
    if (t.options[0] == "baseline")
      t.type = "baseline";
    else
      {
      setdefault();
      return false;
      }

    double startl;

    optionlist optlist;
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&tgrid);
    optlist.push_back(&gridchoice);
    optlist.push_back(&numberquantiles);
    optlist.push_back(&numberbetween);
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&lower);
    optlist.push_back(&catspecific);
    optlist.push_back(&gridsize);
    optlist.push_back(&reference);

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
   t.options = vector<ST::string>(13);
   t.options[0] = t.type;
   t.options[1] = ST::inttostring(degree.getvalue());
   t.options[2] = ST::inttostring(numberknots.getvalue());
   t.options[3] = ST::inttostring(tgrid.getvalue());
   t.options[4] = gridchoice.getvalue();
   t.options[5] = ST::inttostring(numberquantiles.getvalue());
   t.options[6] = ST::inttostring(numberbetween.getvalue());
   t.options[7] = ST::doubletostring(lambda.getvalue());
   t.options[8] = ST::doubletostring(lambdastart.getvalue());
   t.options[9] = lower.getvalue();
   if(catspecific.getvalue())
     {
     t.options[10] = "true";
     }
   else
     {
     t.options[10] = "false";
     }
   t.options[11] = ST::inttostring(gridsize.getvalue());
   t.options[12] = ST::inttostring(reference.getvalue());

   if (lambda.getvalue() < 0)
     {
     setdefault();
     return false;
     }

    int b = t.options[8].strtodouble(startl);
    if (b==1)
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

bool term_baseline_remlreg::checkvector(const vector<term> & terms,
                                        const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "baseline")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//--- class term_baseline_varcoef_remlreg: implementation of member functions --
//------------------------------------------------------------------------------

term_baseline_varcoeff_remlreg::term_baseline_varcoeff_remlreg(void)
  {
  type = "term_baseline_varcoeff";
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",1000,0,10000000);
  catspecific = simpleoption("catspecific",false);
  gridsize = intoption("gridsize",-1,10,500);
  reference = doubleoption("reference",-9999,-10000000,10000000);
  }

void term_baseline_varcoeff_remlreg::setdefault(void)
  {
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  gridsize = intoption("gridsize",-1,10,500);
  reference.setdefault();
  }

bool term_baseline_varcoeff_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==2)  && (t.options.size() >= 1)
        && (t.options.size() <= 6) )
    {
    if (t.options[0] == "baseline")
      t.type = "baseline_varcoeff";
    else
      {
      setdefault();
      return false;
      }

    double startl;

    optionlist optlist;
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);
    optlist.push_back(&gridsize);
    optlist.push_back(&reference);

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
   t.options = vector<ST::string>(6);
   t.options[0] = t.type;
   t.options[1] = ST::doubletostring(lambda.getvalue());
   t.options[2] = ST::doubletostring(lambdastart.getvalue());
   if(catspecific.getvalue())
     {
     t.options[3] = "true";
     }
   else
     {
     t.options[3] = "false";
     }
   t.options[4] = ST::inttostring(gridsize.getvalue());
   t.options[5] = ST::doubletostring(reference.getvalue());

   if (lambda.getvalue() < 0)
     {
     setdefault();
     return false;
     }

    int b = t.options[2].strtodouble(startl);
    if (b==1)
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

bool term_baseline_varcoeff_remlreg::checkvector(const vector<term> & terms,
                                                 const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "baseline_varcoeff")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//-- class term_interactpspline_remlreg: implementation of member functions ----
//------------------------------------------------------------------------------

term_interactpspline_remlreg::term_interactpspline_remlreg(void)
  {
  type = "term_interactpspline";
  degree=intoption("degree",3,1,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  gridsizex=intoption("gridsizex",0,0,500);
  gridsizey=intoption("gridsizey",0,0,500);
  }

void term_interactpspline_remlreg::setdefault(void)
  {
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  gridsizex.setdefault();
  gridsizey.setdefault();
  }

bool term_interactpspline_remlreg::check(term & t)
  {
  optionlist optlist;
  optlist.push_back(&degree);
  optlist.push_back(&numberknots);
  optlist.push_back(&lambda);
  optlist.push_back(&lambdastart);
  optlist.push_back(&catspecific);
  optlist.push_back(&gridsizex);
  optlist.push_back(&gridsizey);

  if ( (t.varnames.size()==2)  && (t.options.size() >= 1)
        && (t.options.size() <= 8) )
    {

    if (t.options[0] == "pspline2dimrw1")
      t.type = "pspline2dimrw1";
    else if (t.options[0] == "pspline2dimrw2")
      t.type = "pspline2dimrw2";
    else if (t.options[0] == "pspline2dimbiharmonic")
      t.type = "pspline2dimbiharmonic";
    else
      {
      setdefault();
      return false;
      }

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
    t.options[1] = ST::inttostring(degree.getvalue());
    t.options[2] = ST::inttostring(numberknots.getvalue());
    t.options[3] = ST::doubletostring(lambda.getvalue());
    t.options[4] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options[5] = "true";
      }
    else
      {
      t.options[5] = "false";
      }
    t.options[6] = ST::inttostring(gridsizex.getvalue());
    t.options[7] = ST::inttostring(gridsizey.getvalue());

    setdefault();
    return true;
    }
  else
    {
    setdefault();
    return false;
    }
  }

bool term_interactpspline_remlreg::checkvector(const vector<term> & terms,
                                               const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "pspline2dimrw1" || terms[i].type == "pspline2dimrw2" ||
      terms[i].type == "pspline2dimbiharmonic")
    return true;

  return false;
  }

//-----------------------------------------------------------------------------------
//- class term_interactpspline_varcoeff_remlreg: implementation of member functions -
//-----------------------------------------------------------------------------------

term_interactpspline_varcoeff_remlreg::term_interactpspline_varcoeff_remlreg(void)
  {
  type = "term_interactpspline_varcoeff";
  degree=intoption("degree",3,1,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  center = simpleoption("center",false);
  }


void term_interactpspline_varcoeff_remlreg::setdefault(void)
  {
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  center.setdefault();
  }

bool term_interactpspline_varcoeff_remlreg::check(term & t)
  {
  optionlist optlist;
  optlist.push_back(&degree);
  optlist.push_back(&numberknots);
  optlist.push_back(&lambda);
  optlist.push_back(&lambdastart);
  optlist.push_back(&catspecific);
  optlist.push_back(&center);

  if ( (t.varnames.size()==3)  && (t.options.size() >= 1)
        && (t.options.size() <= 7) )
    {
    if (t.options[0] == "pspline2dimrw1")
      t.type = "varpspline2dimrw1";
    else
      {
      setdefault();
      return false;
      }

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
    t.options = vector<ST::string>(7);
    t.options[0] = t.type;
    t.options[1] = ST::inttostring(degree.getvalue());
    t.options[2] = ST::inttostring(numberknots.getvalue());
    t.options[3] = ST::doubletostring(lambda.getvalue());
    t.options[4] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options[5] = "true";
      }
    else
      {
      t.options[5] = "false";
      }
    if(center.getvalue())
      {
      t.options[6] = "true";
      }
    else
      {
      t.options[6] = "false";
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


bool term_interactpspline_varcoeff_remlreg::checkvector(const vector<term> & terms,
                                       const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "varpspline2dimrw1")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//------ class term_geospline_remlreg: implementation of member functions ------
//------------------------------------------------------------------------------

term_geospline_remlreg::term_geospline_remlreg(void)
  {
  type = "term_geospline";
  map=stroption("map");
  degree=intoption("degree",3,1,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  gridsizex = intoption("gridsizex",0,0,500);
  gridsizey = intoption("gridsizey",0,0,500);
  }


void term_geospline_remlreg::setdefault(void)
  {
  map.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  gridsizex.setdefault();
  gridsizey.setdefault();
  }

bool term_geospline_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==1)  && (t.options.size() >= 1)
        && (t.options.size() <= 9) )
    {
    if (t.options[0] == "geosplinerw1" || t.options[0] == "geospline")
      t.type = "geosplinerw1";
    else if (t.options[0] == "geosplinerw2")
      t.type = "geosplinerw2";
    else if (t.options[0] == "geosplinebiharmonic")
      t.type = "geosplinebiharmonic";
    else
      {
      setdefault();
      return false;
      }

    optionlist optlist;
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&map);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);
    optlist.push_back(&gridsizex);
    optlist.push_back(&gridsizey);

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
    t.options[1] = ST::inttostring(degree.getvalue());
    t.options[2] = ST::inttostring(numberknots.getvalue());
    t.options[3] = ST::doubletostring(lambda.getvalue());
    t.options[4] = map.getvalue();
    t.options[5] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options[6] = "true";
      }
    else
      {
      t.options[6] = "false";
      }
    t.options[7] = ST::inttostring(gridsizex.getvalue());
    t.options[8] = ST::inttostring(gridsizey.getvalue());

    setdefault();
    return true;
    }
  else
    {
    setdefault();
    return false;
    }
  }


bool term_geospline_remlreg::checkvector(const vector<term> & terms,
                                         const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "geosplinerw1" || terms[i].type == "geosplinerw2" ||
      terms[i].type == "geosplinebiharmonic")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//-- class term_geospline_varcoeff_remlreg: implementation of member functions -
//------------------------------------------------------------------------------

term_geospline_varcoeff_remlreg::term_geospline_varcoeff_remlreg(void)
  {
  type = "term_geospline_varcoeff";
  map=stroption("map");
  degree=intoption("degree",3,1,5);
  numberknots=intoption("nrknots",20,5,500);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  center = simpleoption("center",false);
  }


void term_geospline_varcoeff_remlreg::setdefault(void)
  {
  map.setdefault();
  degree.setdefault();
  numberknots.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  center.setdefault();
  }

bool term_geospline_varcoeff_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==2)  && (t.options.size() >= 1)
        && (t.options.size() <= 8) )
    {
    if (t.options[0] == "geospline")
      t.type = "vargeospline";
    else
      {
      setdefault();
      return false;
      }

    optionlist optlist;
    optlist.push_back(&degree);
    optlist.push_back(&numberknots);
    optlist.push_back(&lambda);
    optlist.push_back(&map);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);
    optlist.push_back(&center);

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
    t.options[1] = ST::inttostring(degree.getvalue());
    t.options[2] = ST::inttostring(numberknots.getvalue());
    t.options[3] = ST::doubletostring(lambda.getvalue());
    t.options[4] = map.getvalue();
    t.options[5] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options [6] = "true";
      }
    else
      {
      t.options [6] = "false";
      }
    if(center.getvalue())
      {
      t.options [7] = "true";
      }
    else
      {
      t.options [7] = "false";
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

bool term_geospline_varcoeff_remlreg::checkvector(const vector<term> & terms,
                                                  const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "vargeospline")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//------ class term_spatial_remlreg: implementation of member functions --------
//------------------------------------------------------------------------------

term_spatial_remlreg::term_spatial_remlreg(void)
  {
  type = "term_spatial";
  map=stroption("map");
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  }

void term_spatial_remlreg::setdefault(void)
  {
  map.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  }


bool term_spatial_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==1) &&
       (t.options.size()<=5) && (t.options.size() >= 1) )
    {
    if (t.options[0] == "spatial")
      t.type = "spatial";
    else
      {
      setdefault();
      return false;
      }

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&map);
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);

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
    t.options = vector<ST::string>(5);
    t.options[0] = t.type;
    t.options[1] = map.getvalue();
    t.options[2] = ST::doubletostring(lambda.getvalue());
    t.options[3] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options[4] = "true";
      }
    else
      {
      t.options[4] = "false";
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

bool term_spatial_remlreg::checkvector(const vector<term> & terms,
                                       const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "spatial")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//--- class term_spatialvarcoef_remlreg: implementation of member functions ----
//------------------------------------------------------------------------------

term_spatial_varcoef_remlreg::term_spatial_varcoef_remlreg(void)
  {
  type = "term_spatial_varcoef";
  map=stroption("map");
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  center = simpleoption("center",false);
  }

void term_spatial_varcoef_remlreg::setdefault(void)
  {
  map.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  center.setdefault();
  }


bool term_spatial_varcoef_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==2) &&
       (t.options.size()<=6) && (t.options.size() >= 1) )
    {
    if (t.options[0] == "spatial")
      t.type = "varcoeffspatial";
    else
      {
      setdefault();
      return false;
      }

    vector<ST::string> opt;
    optionlist optlist;
    optlist.push_back(&map);
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);
    optlist.push_back(&center);

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
    t.options = vector<ST::string>(6);
    t.options[0] = t.type;
    t.options[1] = map.getvalue();
    t.options[2] = ST::doubletostring(lambda.getvalue());
    t.options[3] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options[4] = "true";
      }
    else
      {
      t.options[4] = "false";
      }
    if(center.getvalue())
      {
      t.options[5] = "true";
      }
    else
      {
      t.options[5] = "false";
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

bool term_spatial_varcoef_remlreg::checkvector(const vector<term> & terms,
                                               const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "varcoeffspatial")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//------ class term_random_remlreg: implementation of member functions ---------
//------------------------------------------------------------------------------

term_random_remlreg::term_random_remlreg(void)
  {
  type = "term_random";
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  }


void term_random_remlreg::setdefault(void)
  {
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  }

bool term_random_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==1)  && (t.options.size()<=4) )
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
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);

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
    t.options = vector<ST::string>(4);
    t.options[0] = t.type;
    t.options[1] = ST::doubletostring(lambda.getvalue());
    t.options[2] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options[3] = "true";
      }
    else
      {
      t.options[3] = "false";
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

bool term_random_remlreg::checkvector(const vector<term> & terms,
                                      const unsigned & i)
  {
  assert(i< terms.size());

  if ( terms[i].type == "random" )
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//----- class term_randomslope_remlreg: implementation of member functions -----
//------------------------------------------------------------------------------

term_randomslope_remlreg::term_randomslope_remlreg(void)
  {
  type = "term_randomslope";
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",10,0,10000000);
  catspecific = simpleoption("catspecific",false);
  }

void term_randomslope_remlreg::setdefault(void)
  {
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  }


bool term_randomslope_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==2)  && (t.options.size()<=4) )
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
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);

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
    t.options = vector<ST::string>(4);
    t.options[0] = t.type;
    t.options[1] = ST::doubletostring(lambda.getvalue());
    t.options[2] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options[3] = "true";
      }
    else
      {
      t.options[3] = "false";
      }

    double startl;

    int b = t.options[2].strtodouble(startl);
    if (b==1)
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

bool term_randomslope_remlreg::checkvector(const vector<term> & terms,
                                           const unsigned & i)
  {
  assert(i< terms.size());

  if ( terms[i].type == "randomslope" )
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//---- class term_kriging_1dim_remlreg: implementation of member functions -----
//------------------------------------------------------------------------------

term_kriging_1dim_remlreg::term_kriging_1dim_remlreg(void)
  {
  type = "term_kriging";
  nu=doubleoption("nu",1.5,0.5,3.5);
  maxdist=doubleoption("maxdist",-1,0.00001,10000);
  lambda = doubleoption("lambda",0.1,0,10000000);
  lambdastart = doubleoption("lambdastart",0.1,0,10000000);
  catspecific = simpleoption("catspecific",false);
  }

void term_kriging_1dim_remlreg::setdefault(void)
  {
  nu.setdefault();
  maxdist.setdefault();
  lambda.setdefault();
  lambdastart.setdefault();
  catspecific.setdefault();
  }

bool term_kriging_1dim_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==1)  && (t.options.size() >=1)
        && (t.options.size() <= 6) )
    {
    if (t.options[0] == "kriging")
      t.type = "1dimkriging";
    else
      {
      setdefault();
      return false;
      }

    optionlist optlist;
    optlist.push_back(&nu);
    optlist.push_back(&maxdist);
    optlist.push_back(&lambda);
    optlist.push_back(&lambdastart);
    optlist.push_back(&catspecific);

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
    t.options = vector<ST::string>(6);
    t.options[0] = t.type;
    t.options[1] = ST::doubletostring(nu.getvalue());
    t.options[2] = ST::doubletostring(maxdist.getvalue());
    t.options[3] = ST::doubletostring(lambda.getvalue());
    t.options[4] = ST::doubletostring(lambdastart.getvalue());
    if(catspecific.getvalue())
      {
      t.options[5] = "true";
      }
    else
      {
      t.options[5] = "false";
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

bool term_kriging_1dim_remlreg::checkvector(const vector<term> & terms,
                                            const unsigned & i)
    {
    assert(i< terms.size());

    if (terms[i].type == "1dimkriging")
      return true;

    return false;
    }

//------------------------------------------------------------------------------
//------- class term_kriging_remlreg: implementation of member functions -------
//------------------------------------------------------------------------------

term_kriging_remlreg::term_kriging_remlreg(void)
  {
  type = "term_kriging";
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
  catspecific = simpleoption("catspecific",false);
  gridsizex=intoption("gridsizex",0,0,500);
  gridsizey=intoption("gridsizey",0,0,500);
  }

void term_kriging_remlreg::setdefault(void)
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
  catspecific.setdefault();
  gridsizex.setdefault();
  gridsizey.setdefault();
  }

bool term_kriging_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==2)  && (t.options.size() >=1)
        && (t.options.size() <= 14) )
    {
    if (t.options[0] == "kriging")
      t.type = "kriging";
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
    optlist.push_back(&catspecific);
    optlist.push_back(&gridsizex);
    optlist.push_back(&gridsizey);

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
    t.options = vector<ST::string>(14);
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
    if(catspecific.getvalue()==true)
      {
      t.options[11] = "true";
      }
    else
      {
      t.options[11] = "false";
      }
    t.options[12] = ST::inttostring(gridsizex.getvalue());
    t.options[13] = ST::inttostring(gridsizey.getvalue());

    setdefault();
    return true;
    }
  else
    {
    setdefault();
    return false;
    }
  }

bool term_kriging_remlreg::checkvector(const vector<term> & terms,
                                       const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "kriging")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//-- class term_kriging_varcoeff_remlreg: implementation of member functions ---
//------------------------------------------------------------------------------

term_kriging_varcoeff_remlreg::term_kriging_varcoeff_remlreg(void)
  {
  type = "term_kriging_varcoeff";
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
  catspecific = simpleoption("catspecific",false);
  }

void term_kriging_varcoeff_remlreg::setdefault(void)
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
  catspecific.setdefault();
  }

bool term_kriging_varcoeff_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==3)  && (t.options.size() >=1)
        && (t.options.size() <= 12) )
    {
    if (t.options[0] == "kriging")
      t.type = "varkriging";
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
    optlist.push_back(&catspecific);

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
    t.options = vector<ST::string>(12);
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
    if(catspecific.getvalue()==true)
      {
      t.options[11] = "true";
      }
    else
      {
      t.options[11] = "false";
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

bool term_kriging_varcoeff_remlreg::checkvector(const vector<term> & terms,
                                                const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "varkriging")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//----- class term_geokriging_remlreg: implementation of member functions ------
//------------------------------------------------------------------------------

term_geokriging_remlreg::term_geokriging_remlreg(void)
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
  catspecific = simpleoption("catspecific",false);
  gridsizex=intoption("gridsizex",0,0,500);
  gridsizey=intoption("gridsizey",0,0,500);
  }

void term_geokriging_remlreg::setdefault(void)
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
  gridsizex.setdefault();
  gridsizey.setdefault();
  }

bool term_geokriging_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==1)  && (t.options.size() >=1)
        && (t.options.size() <= 15) )
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
    optlist.push_back(&catspecific);
    optlist.push_back(&gridsizex);
    optlist.push_back(&gridsizey);

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
    t.options = vector<ST::string>(15);
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
    if(catspecific.getvalue()==true)
      {
      t.options[12] = "true";
      }
    else
      {
      t.options[12] = "false";
      }
    t.options[13] = ST::inttostring(gridsizex.getvalue());
    t.options[14] = ST::inttostring(gridsizey.getvalue());

    setdefault();
    return true;
    }
  else
    {
    setdefault();
    return false;
    }
  }

bool term_geokriging_remlreg::checkvector(const vector<term> & terms,
                                          const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "geokriging")
    return true;

  return false;
  }

//------------------------------------------------------------------------------
//- class term_geokriging_varcoeff_remlreg: implementation of member functions -
//------------------------------------------------------------------------------

term_geokriging_varcoeff_remlreg::term_geokriging_varcoeff_remlreg(void)
  {
  type = "term_geokriging_varcoeff";
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
  catspecific = simpleoption("catspecific",false);
  }

void term_geokriging_varcoeff_remlreg::setdefault(void)
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
  catspecific.setdefault();
  }

bool term_geokriging_varcoeff_remlreg::check(term & t)
  {
  if ( (t.varnames.size()==2)  && (t.options.size() >=1)
        && (t.options.size() <= 13) )
    {
    if (t.options[0] == "geokriging")
      t.type = "vargeokriging";
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
    optlist.push_back(&catspecific);

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
    t.options = vector<ST::string>(13);
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
    if(catspecific.getvalue()==true)
      {
      t.options[12] = "true";
      }
    else
      {
      t.options[12] = "false";
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

bool term_geokriging_varcoeff_remlreg::checkvector(const vector<term> & terms,
                                                   const unsigned & i)
  {
  assert(i< terms.size());

  if (terms[i].type == "vargeokriging")
    return true;

  return false;
  }


