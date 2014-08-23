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





#include "option.h"

//------------------------------------------------------------------------------
//------------ CLASS option: implementation of member functions ----------------
//------------------------------------------------------------------------------


option:: option(void)
  {
  optionname = "";
  valuechanged = false;
  maxtoken = 3;
  }


option::option(const ST::string & n)
  {
  optionname = n;
  valuechanged = false;
  maxtoken=3;
  }


option::option(const option & b)
  {
  optionname = b.optionname;
  errormessages = b.errormessages;
  valuechanged = b.valuechanged;
  maxtoken = b.maxtoken;
  }


const option & option::operator=(const option & b)
  {
  if (this == &b)
	 return *this;
  optionname = b.optionname;
  errormessages = b.errormessages;
  valuechanged = b.valuechanged;
  maxtoken = b.maxtoken;
  return *this;
  }


//------------------------------------------------------------------------------
//---------- CLASS simpleoption: implementation of member functions ------------
//------------------------------------------------------------------------------


simpleoption::simpleoption(const ST::string & n,bool v) : option(n)
  {
  defaultvalue = v;
  value = v;
  maxtoken=1;
  }


simpleoption::simpleoption(const simpleoption & o) : option(option(o))
  {
  value = o.value;
  defaultvalue = o.defaultvalue;
  }


const simpleoption & simpleoption::operator=(const simpleoption & o)
  {
  if (this == &o)
	 return *this;
  option::operator=(option(o));
  value = o.value;
  defaultvalue = o.defaultvalue;
  return *this;
  }


int simpleoption::parse(const ST::string & c)
  {

  errormessages.clear();

  if (c==optionname)
	 {
	 value = true;
     valuechanged = true;
	 return 1;
	 }
  else
    return 0;

  }


//------------------------------------------------------------------------------
//----------- CLASS fileoption: implementation of member functions -------------
//------------------------------------------------------------------------------


fileoption::fileoption (const ST::string & n,const ST::string & v,bool t) : option(n)
  {
  defaultvalue = v;
  value = v;
  filetype = t;
  }


fileoption::fileoption(const fileoption & o) : option(option(o))
  {
  value = o.value;
  defaultvalue = o.defaultvalue;
  filetype = o.filetype;
  }


const fileoption & fileoption::operator=(const fileoption & o)
  {
  if (this==&o)
	 return *this;
  optionname = o.optionname;
  errormessages = o.errormessages;
  value = o.value;
  defaultvalue = o.defaultvalue;
  filetype = o.filetype;
  return *this;
  }


int fileoption::parse(const ST::string & c)
  {

  errormessages.clear();

  vector<ST::string> token = c.strtoken(" =");

  if ( (! token.empty()) && (token[0] == optionname ) )
	 {

	 if (token.size() < 2 || token[1] != "=")
		errormessages.push_back("ERROR in option " + optionname + ": \"=\" expected\n");
	 if (token.size() < 3)
		errormessages.push_back ("ERROR in option " + optionname + ": filename specification expected\n");
	 if (token.size() > 3)
		errormessages.push_back("ERROR in option " + optionname + ": invalid option specification");

	 if (errormessages.empty())
		{

        ST::string t = token[2];

		int isex = token[2].isvalidfile();

		if (isex == 1)
		  errormessages.push_back("ERROR in option " + optionname + ": " + token[2] +
										  " is not a valid filename\n");
		else
          {
		  value = token[2];
          valuechanged = true;
          }

		}

	 return 1;

	 }
  else
    return 0;

  }


//------------------------------------------------------------------------------
//----------- CLASS fileoption2: implementation of member functions ------------
//------------------------------------------------------------------------------


fileoption2::fileoption2 (const ST::string & n,
const ST::string & v,bool t) : fileoption(n,v,t)
  {
  }


fileoption2::fileoption2(const fileoption2 & o) : fileoption(fileoption(o))
  {
  }


const fileoption2 & fileoption2::operator=(const fileoption2 & o)
  {
  if (this==&o)
	 return *this;
  fileoption::operator=(fileoption(o));
  return *this;
  }


int fileoption2::parse(const ST::string & c)
  {

  errormessages.clear();

  vector<ST::string> token = c.strtoken(" =");

  if ( (! token.empty()) && (token[0] == optionname ) )
	 {

	 if (token.size() < 2 || token[1] != "=")
		errormessages.push_back("ERROR in option " + optionname + ": \"=\" expected\n");
	 if (token.size() < 3)
		errormessages.push_back ("ERROR in option " + optionname + ": filename specification expected\n");
	 if (token.size() > 3)
		errormessages.push_back("ERROR in option " + optionname + ": invalid option specification");

	 if (errormessages.empty())
		{

		int isex = token[2].isexistingfile();

		if (isex == 1)
		  errormessages.push_back("ERROR in option " + optionname + ": " + token[2] +
										  " is not an existing filename\n");
		else
          {
		  value = token[2];
          valuechanged = true;
          }

		}

	 return 1;

	 }
  else
    return 0;

  }



//------------------------------------------------------------------------------
//----------- CLASS stroption: implementation of member functions --------------
//------------------------------------------------------------------------------


stroption::stroption(void) : option()
  {
  allallowed = false;
  admissible = vector<ST::string>();
  value = "";
  defaultvalue = "";
  }


stroption::stroption(const ST::string & n,const vector<ST::string> & adm,
							const ST::string & v) : option(n)
  {
  admissible = adm;
  allallowed = false;
  defaultvalue = v;
  value = v;
  }

stroption::stroption(const ST::string & n,const ST::string & v) : option(n)
  {
  allallowed = true;
  defaultvalue = v;
  value = v;
  }

stroption::stroption(const ST::string & n) : option(n)
  {
  allallowed  = true;
  value = "";
  defaultvalue="";
  }


stroption::stroption(const stroption & b)
  {
  allallowed = b.allallowed;
  optionname = b.optionname;
  errormessages = b.errormessages;
  admissible = b.admissible;
  defaultvalue = b.defaultvalue;
  value = b.value;
  }


const stroption & stroption::operator=(const stroption & b)
  {
  if (this == &b)
	 return *this;
  allallowed = b.allallowed;
  optionname = b.optionname;
  errormessages = b.errormessages;
  admissible = b.admissible;
  defaultvalue = b.defaultvalue;
  value = b.value;
  return *this;
  }


int stroption::parse(const ST::string & c)
  {

  errormessages.clear();

  vector<ST::string> token = c.strtoken("=");

  if(!token.empty())
    {
    token[0] = token[0].eatallwhitespace();
    token[0] = token[0].eatallcarriagereturns();
    }

  if (token.size() >= 3)
    {
    token[1] = token[1].eatallwhitespace();
    token[1] = token[1].eatallcarriagereturns();
    token[2] = token[2].eatwhitespace();
    token[2] = token[2].eatallcarriagereturns();
    }

  if ( (! token.empty()) && (token[0] == optionname) )
	 {

     if (token.size() < 2 || token[1] != "=")
		errormessages.push_back("ERROR in option " + optionname + ": \"=\" expected\n");

	 if (token.size() < 3)
		errormessages.push_back("ERROR in option " + optionname + ": new value expected\n");

	 if (token.size() > 3)
		errormessages.push_back("ERROR in option " + optionname + ": invalid option specification\n");

	 if (errormessages.empty())
		{
		if (allallowed == false)
		  {
		  int i = 0;
		  int h = -1;
		  while ((i < admissible.size()) && (h == -1))
			 {
			 if (token[2] == admissible[i])
				h=i;
			 i++;
			 }
		  if (h == -1)
			 errormessages.push_back("ERROR in option " + optionname + ": " + token[2] + " unknown value\n");
		  } // end: if {allalllowed == false)

		} // end: if (errormessages.empty())

	 if (errormessages.empty())
       {
       value = token[2];
       valuechanged = true;
       }


	 return 1;

	 }
  else
	 return 0;

  }


//------------------------------------------------------------------------------
//----------- CLASS intoption: implementation of member functions --------------
//------------------------------------------------------------------------------


intoption::intoption(void) : option()
  {
  lowerbound = INT_MIN;
  upperbound = INT_MAX;
  defaultvalue = 0;
  value = 0;
  }


intoption::intoption(const ST::string & n,const int v,const int lb,const int ub)
							: option(n)
  {
  lowerbound = lb;
  upperbound = ub;
  defaultvalue = v;
  value = v;
  }


intoption::intoption(const intoption & b)
  {
  optionname = b.optionname;
  errormessages = b.errormessages;
  lowerbound = b.lowerbound;
  upperbound = b.upperbound;
  defaultvalue = b.defaultvalue;
  value = b.value;
  }


const intoption & intoption::operator=(const intoption & b)
  {
  if (this == & b)
	 return *this;
  optionname = b.optionname;
  errormessages = b.errormessages;
  lowerbound = b.lowerbound;
  upperbound = b.upperbound;
  defaultvalue = b.defaultvalue;
  value = b.value;
  return *this;
  }


int intoption::parse(const ST::string & c)
  {

  errormessages.clear();

  vector<ST::string> token = c.strtoken(" =");

  if ( (! token.empty()) && (token[0] == optionname) )
	 {

	 long help;

	 if (token.size() < 2 || token[1] != "=")
		errormessages.push_back("ERROR in option " + optionname + ": \"=\" expected\n");

	 if (token.size() < 3)
		errormessages.push_back("ERROR in option " + optionname + ": new value expected\n");

	 if (token.size() > 3)
		errormessages.push_back("ERROR in option " + optionname + ": invalid option specification\n");

	 if (errormessages.empty())
		{
		if (token[2].strtolong(help) == 1)
		  errormessages.push_back("ERROR in option " + optionname + ": integer value expected\n");
		else if ((help < lowerbound) || (help > upperbound))
		  errormessages.push_back("ERROR in option " + optionname + ": value between " +
		  ST::inttostring(lowerbound) + " and " + ST::inttostring(upperbound) +
		  " expected\n");
		}

	 if (errormessages.empty())
       {
	   value = help;
       valuechanged = true;
       }


	 return 1;

	 }
  else
	 return 0;

  }


//------------------------------------------------------------------------------
//------------ CLASS doubleoption: implementation of member functions ----------
//------------------------------------------------------------------------------

void doubleoption::setdefault(void)
  {
  value = defaultvalue;
  valuechanged = false;
  }


doubleoption::doubleoption(void) : option()
  {
  lowerbound = MINDOUBLE;
  upperbound = MAXDOUBLE;
  defaultvalue = 0;
  value =0;
  }


doubleoption::doubleoption(const ST::string & n,const double v,const double lb,
									const double ub) : option(n)
  {
  lowerbound = lb;
  upperbound = ub;
  defaultvalue = v;
  value = v;
  }


doubleoption::doubleoption(const doubleoption & b)
  {
  optionname = b.optionname;
  errormessages = b.errormessages;
  lowerbound = b.lowerbound;
  upperbound = b.upperbound;
  defaultvalue = b.defaultvalue;
  value = b.value;
  }


const doubleoption & doubleoption::operator=(const doubleoption & b)
  {
  if (this == & b)
	 return *this;
  optionname = b.optionname;
  errormessages = b.errormessages;
  lowerbound = b.lowerbound;
  upperbound = b.upperbound;
  value = b.value;
  defaultvalue = b.defaultvalue;
  return *this;
  }


int doubleoption::parse(const ST::string & c)
  {

  errormessages.clear();

  vector<ST::string> token = c.strtoken(" =");

  if ( (! token.empty()) && (token[0] == optionname) )
	 {

	 double help;

	 if (token.size() < 2 || token[1] != "=")
		errormessages.push_back("ERROR in option " + optionname + ": \"=\" expected");

	 if (token.size() < 3)
		errormessages.push_back("ERROR in option " + optionname + ": new value expected\n");

	 if (token.size() > 3)
		errormessages.push_back("ERROR in option " + optionname + ": invalid option specification\n");

	 if (errormessages.empty())
		{
		if (token[2].strtodouble(help) == 1)
		  errormessages.push_back("ERROR in option " + optionname + ": real value expected\n");
		else if ((help < lowerbound) || (help > upperbound))
          {
		  errormessages.push_back("ERROR in option " + optionname + ": value between " +
		  ST::doubletostring(lowerbound) + " and " + ST::doubletostring(upperbound) +
		  " expected\n");
          }
		}

	 if (errormessages.empty())
       {
       value = help;
       valuechanged = true;
       }

	 return 1;

	 }
  else
	 return 0;

  }


//------------------------------------------------------------------------------
//------------ CLASS optionlist: implementation of member functions ------------
//------------------------------------------------------------------------------


optionlist::iterator optionlist::find(const ST::string & optionname)
  {
  optionlist::iterator f = end();
  optionlist::iterator i = begin();
  while ((i != end()) && (f == end()))
	 {
	 if ((*i)->getname() == optionname)
		f = i;
	 i++;
	 }
  return f;
  }


const optionlist & optionlist::operator=(const optionlist & ol)
  {
  if (this == &ol)
	 return *this;
  vector<option*>::operator=(vector<option*>(ol));
  return *this;
  }


int optionlist::parse(const ST::string & c, bool clearerrors)
  {
  if (clearerrors == true)
	 errormessages.clear();

  int recognized = 0;
  optionlist::iterator i = begin();
  while ((i != end()) && (recognized == 0))
	 {
	 if ((*i)->parse(c) == 1)
		{
		recognized = 1;
//		errormessages.insert_back((*i)->geterrormessages());
        if (! ((*i)->geterrormessages()).empty())
          errormessages.insert(errormessages.end(),
          ((*i)->geterrormessages()).begin(),((*i)->geterrormessages()).end());

		}
	 i++;
	 }
  return recognized;
  }


void optionlist::parsemultiple(const ST::string & c)
  {

  errormessages.clear();
  setdefault();

  if (c.length() > 0)
	 {
	 vector<ST::string> token;
     int ok = c.strtoken_quot(token," =");
     if (ok == 1)
       {
	   vector<ST::string> options;
       unsigned maxtoken=1;
       unsigned t=1;
//	   options.push_back(token[0]);
//	   int j=1;
	   int j=0;
	   int k=-1;
       optionlist::iterator it;
	   while (j < token.size())
         {
         ST::string to = token[j];
         it = find(token[j]);
         if (it == end())
           {
           if (t < maxtoken)
             {
             options[k] = options[k] + " " + token[j];
             t++;
             }
           else
             {
   		     options.push_back(token[j]);
  	         errormessages.push_back("ERROR: " + token[j] + " unknown option\n");
             maxtoken = 1;
             t=1;
		     k++;
             }
           }
         else
           {
		   options.push_back(token[j]);
           maxtoken = (*it)->getmaxtoken();
           t=1;
		   k++;
           }
		 j++;
		 }

	   j=0;
	   while ( (j<options.size()) && (errormessages.empty()) )
		  {
		  if (parse(options[j],false) == 0)
		    errormessages.push_back("ERROR: unknown option\n");
		  j++;
		  }
       }
     else
       errormessages.push_back("ERROR: \" required\n");
     }

//   unsigned r = errormessages.size();
     
  }


void optionlist::setdefault(void)
  {
  if (! empty())
	 {
	 optionlist::iterator i;
	 for(i=begin();i!=end();++i)
		(*i)->setdefault();
	 }
  }


