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





#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#include<statwin_haupt.h>
#endif

#include"data.h"
#include<time.h>

using std::ios;

//------------------------------------------------------------------------------
//------------- CLASS data: implementation of member functions -----------------
//------------------------------------------------------------------------------


void data::indexcreate(void)
  {
  index = vector<int>(obs());
  unsigned i;
  for (i=0;i<index.size();i++)
	 index[i] = i;
  }


int data::makeitlist(list<ST::string> & names,
					  vector<list<realvar>::iterator> & itl)
  {
  int correct = 0;
  itl = vector<list<realvar>::iterator>();
  list<ST::string>::iterator i,k;
  list<realvar>::iterator j;
  k= names.begin();
  while (k != names.end() && (correct == 0))
	 {
	 if ( (*k) == "const")
		k++;
	 else
		{
		if (findvar(*k,i,j) == 0)
		  itl.push_back(j);
		else
		  correct = 1;
		k++;
		}
	 }
  return correct;
  }


data::data(list <ST::string> & vn,list<realvar> & va)
  {
  varnames = vn;
  variables = va;
  if (!(*variables.begin()).empty())
	 empty = false;
  else
	 empty = true;
  indexcreate();
  }


int data::findvar(const ST::string & name,list<ST::string>::iterator & stit,
				  list<realvar>::iterator & varit)
  {
  if (empty)
   return 1;
  else
	 {
	 stit = varnames.begin();
	 varit = variables.begin();
	 while (stit != varnames.end())
		{
		if ((*stit) == name)
		  break;
		stit++;
		varit++;
		}
	 if (stit == varnames.end())
		return 1;
	 else
		return 0;
	 }
  }


data::data(const data & d)
  {
  varnames = d.varnames;
  variables = d.variables;
  empty = d.empty;
  index = d.index;
  }


const data & data::operator=(const data & d)
  {
  if (this == &d)
	 return *this;
  varnames = d.varnames;
  variables = d.variables;
  empty = d.empty;
  index = d.index;
  return *this;
  }


ostream & operator<<(ostream & out,data & d)
  {
  if (! d.empty)
	 {
	 out.setf(ios::left);
	 list<ST::string>::iterator p;
	 for(p=d.varnames.begin();p != d.varnames.end();p++)
		{
		out.width(10);
		out << *p << " ";
		}
	 out << endl;

	 unsigned i;
	 list<realvar>::iterator j;
	 for(i=0;i<d.obs();i++)
		{
		for(j=d.variables.begin();j!=d.variables.end();j++)
		  {
		  out.width(10);
		  out << (*j)[d.index[i]] << " ";
		  }
		out << endl;
		}
	 }
  else
	 out << "dataset is empty" << endl;
  return out;
  }


unsigned data::obs() const
  {
  if (empty)
	 return 0;
  else
	 return (*variables.begin()).size();
  }


void data::clear()
  {
  empty = true;
  list<realvar>::iterator i;
  for(i=variables.begin();i!=variables.end();i++)
	 (*i).clear();
  if (variables.size() > 0)
	 variables.erase(variables.begin(),variables.end());
  if (varnames.size() > 0)
	 varnames.erase(varnames.begin(),varnames.end());
  index = vector<int>();
  }


//------------------------------------------------------------------------------
//------------- CLASS filter: implementation of member functions  --------------
//------------------------------------------------------------------------------


filter::filter(vector<bool> & v) : vector<bool>(v)
  {
  sum = 0;
  vector<bool>::iterator i;
  for(i=v.begin();i!=v.end();++i)
	 if (*i == false)
		sum++;
  }


const filter & filter::operator=(vector<bool> & v)
  {
  vector<bool>::operator=(v);
  vector<bool>::iterator i;
  sum = 0;
  for(i=v.begin();i!=v.end();++i)
	 if (*i == false)
		sum++;
  return *this;
  }


void filter::filterNA(::data & d, list<ST::string> & names)
  {

  if (size() != d.obs())
	 *this = filter(d.obs());

  vector<list<realvar>::iterator> l;
  if (d.makeitlist(names,l) == 0)
	 {

	 unsigned i,j;
	 bool missing;
	 sum = 0;
	 for(i=0;i<d.obs();i++)
		{

		missing = false;
		j = 0;
		while ((j < l.size()) && (missing == false))
		  {
		  if ((*l[j])[i] == NA)
			 missing = true;

		  j++;
		  }
		if (missing == true)
		  {
		  sum++;
		  filter::operator[](i) = 0;
		  }
		else
		  filter::operator[](i) = 1;
		}
	 }
  }


filter filter::operator+(filter & f)
  {
  filter help(size());
  if (! empty())
	 {
	 filter::iterator pos1,pos2,pos3;
	 for(pos1 = begin(),pos2 = f.begin(),pos3=help.begin();pos1 != end();
		  ++pos1,++pos2,++pos3)
		if ( (*pos1 == false) || (*pos2 == false) )
		  {
		  help.sum++;
		  *pos3 = false;
		  }
		else
		 *pos3 = true;
	 }
  return help;
  }


filter filter::operator+(realvar & v)
  {
  filter help(size());
  if (! empty())
	 {
	 filter::iterator pos1,pos3;
	 realvar::iterator pos2;
	 for(pos1 = begin(),pos2 = v.begin(),pos3=help.begin();pos1 != end();
		  ++pos1,++pos2,++pos3)
		if ( (*pos1 == false) || ((*pos2) <= 0.0) )
		  {
		  help.sum++;
		  *pos3 = false;
		  }
		else
		 *pos3 = true;
	 }
  return help;
  }


//------------------------------------------------------------------------------
//-------------- CLASS dataset: implementation of member functions -------------
//------------------------------------------------------------------------------


dataset::dataset(const dataset & d)
  {
  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = d.adminb_p;
  #endif
  name = d.name;
  nrobs = d.nrobs;
  datarep = d.datarep;
  f = d.f;
  errormessages = d.errormessages;
  }


const dataset & dataset::operator=(const dataset & d)
  {
  if (this == &d)
	 return *this;
  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = d.adminb_p;
  #endif
  name = d.name;
  nrobs = d.nrobs;
  datarep = d.datarep;
  f = d.f;
  errormessages = d.errormessages;
  return *this;
  }


void dataset::filldata(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adminb_p,
#endif
istream & in,ST::string & m,const unsigned & maxobs)
  {
  datarep.variables = list<realvar>(datarep.varnames.size());

  list<realvar>::iterator i;

  bool stop;
  i = datarep.variables.begin();
  while (i != datarep.variables.end())
    {
    (*i).reserve(maxobs);
    i++;
    }


  ST::string obs;
  ST::string missing;
  if (m.length() == 0)
	 missing = ".";
  else
	 missing = m;
  double v;

  while ((!in.eof()) && (errormessages.empty()))
	 {
	 i = datarep.variables.begin();
	 while ((!in.eof()) && (i != datarep.variables.end())
				  && (errormessages.empty()))
		{
		in >> obs;
		if (obs.length() > 0 && !in.eof())
		  {
		  if ((obs == ".") || (obs == "NA") || (obs == missing) )
			 (*i).push_back(NA);
		  else
			 {
			 if (obs.strtodouble(v) == 1)
				errormessages.push_back("ERROR: "  + obs +
				" cannot be read as a number\n");
			 else
				(*i).push_back(v);
			 }
		  }
		i++;
		}  // end: while ....

     #if defined(BORLAND_OUTPUT_WINDOW)
     stop = hauptformular->breakcommand();
     #elif defined(JAVA_OUTPUT_WINDOW)
     stop = adminb_p->breakcommand();
     #else
     stop = false;
     #endif
     if (stop)
       errormessages.push_back("ERROR: No observations read\n");
     if (stop)
       break;


	 } // end: while ((! in.eof()) && (errormessages.empty()))
  i = datarep.variables.begin();
  int o = (*i).size();
  i++;


  while ( (i != datarep.variables.end()) && (errormessages.empty()) )
    {
    if ((*i).size() != o)
      errormessages.push_back(
      "ERROR: missing observations for one or more variable\n");
      i++;
    }

  }


void dataset::checkvarnames(void)
  {
  list<ST::string>::iterator i = datarep.varnames.begin();
  while ((i != datarep.varnames.end()) && (errormessages.empty()))
	 {
	 if ((*i).isvarname() == 1)
        {
		errormessages.push_back("ERROR: " + (*i) + " invalid variable name\n");
        }
	 i++;
	 }
  }


void dataset::read(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adminb_p,
#endif
ifstream & in,ST::string & missing,
                   const unsigned & maxobs, const list<ST::string> & names)
  {
  datarep.clear();
  errormessages.clear();

  if (names.empty())
	 {
	 ST::string h;
	 ST::getline(in,1000000,h,'\n');
     h = h.eatallcarriagereturns();
	 datarep.varnames = h.strtokenlist(" \t",false);
	 }
  else
	 datarep.varnames = names;

  checkvarnames();


  if (errormessages.empty())
	 filldata(
     #if defined(JAVA_OUTPUT_WINDOW)
     adminb_p,
     #endif
     in,missing,maxobs);


  if (! errormessages.empty())
	 datarep.clear();
  else
	 {
	 datarep.empty = false;
	 f = filter(datarep.obs());
	 nrobs = datarep.obs();
	 datarep.indexcreate();
	 }

  }


unsigned dataset::write (
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adminb_p,
#endif
ostream & out, list<ST::string> & names,const bool header,const realvar & v0)
  {
  errormessages.clear();

  unsigned nrwritten=0;
  realvar v = v0;

  if (v.empty())
    v = realvar(nrobs,1);

  if (names.empty())
	 names = datarep.varnames;

  list<ST::string>::iterator i;
  for (i=names.begin();i!=names.end();++i)
	 if (datarep.findvar(*i) == 1)
		 errormessages.push_back(
		 "ERROR: variable " + (*i) + " can not be found\n");

	if (errormessages.empty())
	  {
	  if (header == true)
		 {
		 for (i=names.begin();i!=names.end();++i)
			out << (*i) << "   ";
		 out << endl;
		 }


	  vector<list<realvar>::iterator> itl;
	  datarep.makeitlist(names,itl);
	  unsigned n = 0;
	  unsigned j;
      bool stop;
	  for (n=0;n<nrobs;n++)
		 {
		 if (v[datarep.index[n]] == 1)
			{
            nrwritten++;
			for (j=0;j<itl.size();j++)
			  out << (*itl[j])[datarep.index[n]] << "   ";
			out << endl;
			}
#if defined(BORLAND_OUTPUT_WINDOW)
         stop = hauptformular->breakcommand();
#elif defined(JAVA_OUTPUT_WINDOW)
         stop = adminb_p->breakcommand();
#else
         stop = false;
#endif
         if (stop)
           {
           errormessages.push_back("ERROR: no observations written to external file\n");
           nrwritten = 0;
           }

         if (stop)
           break;


		 }

	  }

  return nrwritten;
  }


realvar dataset::eval_exp(ST::string  expression, bool clearerrors)
  {

//  randomize();

  if (clearerrors == true)
	 errormessages.clear();

  expression = expression.eatallwhitespace();

  list<ST::string>::iterator stit;
  list<realvar>::iterator varit;

  ST::string functionname;
  ST::string argument;
  int isfunc = expression.isfunction(functionname,argument);

  ST::string varname;
  ST::string argument2;
  int issubscr = expression.issubscribing(varname,argument2);

  ST::string sign;
  int precedencepos = expression.lowestprecedencepos(sign);

  double value;
  realvar valuevek;
  realvar valuevek2;

  vector<ST::string> arglist;

  if (datarep.findvar(expression,stit,varit) == 0)  // exp = valid variable name
	 return (*varit);
  else if ( (expression[0] == '(') &&
				(expression.closingbracketpos(0) == expression.length()-1) )
	 return eval_exp(expression.substr(1,expression.length()-2),false);

  else if (expression == "_n")
	 {
	 realvar v(nrobs);
	 for (unsigned i =0;i<nrobs;i++)
		v[datarep.index[i]] = i+1;
	 return v;
	 }
  else if (expression == "_pi")
	 return realvar(nrobs,3.14159265);
  else if ( (expression == "NA") || (expression == ".") )
	 return realvar(nrobs,NA);
  else if (expression.strtodouble(value) == 0)
	 return realvar(nrobs,value);
  else if (expression == "_N")
	 return realvar(nrobs,double(nrobs));
  else if (precedencepos >=0)
	 {
	 if (precedencepos == 0)
		{
		if (expression.length() > 1)
		  {
		  if (sign == "+")
			 return eval_exp(expression.substr(1,expression.length()-1),false);
		  else if (sign == "-")
			 return -eval_exp(expression.substr(1,expression.length()-1),false);
		  else
			 {
			 errormessages.push_back(
			 "ERROR: expression syntax error in \"" + expression + "\"\n");
			 return realvar(nrobs);
			 }
		  }
		else
		  {
		  errormessages.push_back(
		  "ERROR: expression syntax error in \"" + expression + "\"\n");
		  return realvar(nrobs);
		  }
		}
	 else if (precedencepos == expression.length()-1)
       {
       errormessages.push_back(
       "ERROR: expression syntax error in \"" + expression + "\"\n");
       return realvar(nrobs);
       }
	 else
		{

		ST::string leftexp = expression.substr(0,precedencepos);
		ST::string rightexp;
		if (sign.length() == 1)
		  rightexp = expression.substr(precedencepos+1,
						 expression.length()-precedencepos-1);
		else
		  rightexp = expression.substr(precedencepos+2,
		  expression.length()-precedencepos-2);

		if (sign == "+")
		  return eval_exp(leftexp,false) + eval_exp(rightexp,false);
		else if (sign == "-")
		  return eval_exp(leftexp,false) - eval_exp(rightexp,false);
		else if (sign == "*")
		  return eval_exp(leftexp,false) * eval_exp(rightexp,false);
		else if (sign == "/")
		  return eval_exp(leftexp,false) / eval_exp(rightexp,false);
		else if (sign == "^")
		  return realob::power(eval_exp(leftexp,false),eval_exp(rightexp,false));
		else if (sign == "<")
		  return (eval_exp(leftexp,false) < eval_exp(rightexp,false));
		else if (sign == ">")
		  return (eval_exp(leftexp,false) > eval_exp(rightexp,false));
		else if (sign == "<=")
		  return (eval_exp(leftexp,false) <= eval_exp(rightexp,false));
		else if (sign == ">=")
		  return (eval_exp(leftexp,false) >= eval_exp(rightexp,false) );
		else if (sign == "=")
//		  return (eval_exp(leftexp,false).isequal(eval_exp(rightexp,false)));
          {
		  realvar help1,help2;
		  help1 = eval_exp(leftexp,false);
		  help2 =eval_exp(rightexp,false);
          return (help1.isequal(help2));
          }
		else if (sign == "!=")
//		  return (eval_exp(leftexp,false).isnotequal(eval_exp(rightexp,false)));
          {
          realvar help1,help2;
          help1= eval_exp(leftexp,false);
          help2 =eval_exp(rightexp,false);
          return (help1.isnotequal(help2));
          }
		else if (sign == "&")
		  return (eval_exp(leftexp,false) && eval_exp(rightexp,false));
		else if (sign == "|")
		  return (eval_exp(leftexp,false) || eval_exp(rightexp,false));
        else
          {
          errormessages.push_back("ERROR: " + sign + " invalid sign\n");
          return realvar(nrobs);
          }

		}
	 } // end: else if (precedencepos >=0)

  else if (issubscr == 1)
	 {
	 if (datarep.findvar(varname,stit,varit) == 0)
		{
		realvar v(nrobs);
        realvar::iterator vit = v.begin();
		valuevek = eval_exp(argument2,false);
        realvar::iterator vvekit=valuevek.begin();
        double ind;
		for (unsigned i=0;i<nrobs;i++,++vit,++vvekit)
		  {
		  ind = (*vvekit).getvalue();
		  if ( (ind < 1) || (ind > nrobs) )
			 *vit = NA;
		  else
			 *vit = (*varit)[long(ind)-1];
		  }
		return v;
		}
	 else
       {
       errormessages.push_back("ERROR: variable " + varname + " not found\n");
       return realvar(nrobs);
       }
	 }

  else if (isfunc == 1)
	 {


	 if (functionname == "sqrt")
		return realob::sqrt(eval_exp(argument,false));
	 else if (functionname == "abs")
		return realob::abs(eval_exp(argument,false));
	 else if (functionname == "exp")
		return realob::exp(eval_exp(argument,false));
	 else if (functionname == "sin")
		return realob::sin(eval_exp(argument,false));
	 else if (functionname == "cos")
		return realob::cos(eval_exp(argument,false));
	 else if (functionname == "floor")
		return realob::floor(eval_exp(argument,false));
	 else if (functionname == "log")
		return realob::log(eval_exp(argument,false));
     else if (functionname == "cumulnorm")
//       return realob::cumulnorm(eval_exp(argument,false));
        {
        realvar help;
        help = eval_exp(argument,false);
        return realob::cumulnorm(help);
        }
	 else if (functionname == "log10")
		return realob::log10(eval_exp(argument,false));
	 else if (functionname == "lag")
//		return lagrealvar(eval_exp(argument,false),datarep.index);
        {
		realvar temp1;
		temp1 = eval_exp(argument,false);
		return lagrealvar(temp1,datarep.index);
        }
	 else if (functionname == "cumul")
//		return cumul(eval_exp(argument,false),datarep.index);
		{
		realvar temp2;
		temp2 = eval_exp(argument,false);
		return cumul(temp2,datarep.index);
		}
	 else if (functionname == "uniform")
		{
		if (argument.length() == 0)
          {
#if defined(JAVA_OUTPUT_WINDOW)
          realvar h(nrobs);
          register unsigned i;
          realvar::iterator it = h.begin();

          jmethodID javauniform = adminb_p->Java->GetMethodID(
          adminb_p->BayesX_cls, "juniform", "()D");

          for (i=0;i<nrobs;i++,++it)
        	 *it = adminb_p->Java->CallDoubleMethod(adminb_p->BayesX_obj, javauniform);
          return h;
#else
		  return realob::uniform(nrobs);
#endif
          }
		else
		  {
		  errormessages.push_back(
		  "ERROR: argument not allowed in function uniform\n");
    	  return realvar(nrobs);
		  }
		}  // end: uniform
	 else if (functionname == "normal")
		{
		if (argument.length() == 0)
          {
#if defined(JAVA_OUTPUT_WINDOW)
          realvar h(nrobs);
          register unsigned i;
          realvar::iterator it = h.begin();

          jmethodID javanormal = adminb_p->Java->GetMethodID(
          adminb_p->BayesX_cls, "jnormal", "()D");

          for (i=0;i<nrobs;i++,++it)
        	 *it = adminb_p->Java->CallDoubleMethod(adminb_p->BayesX_obj, javanormal);
          return h;
#else
		  return realob::normal(nrobs);
#endif
          }
		else
		  {
		  errormessages.push_back(
		  "ERROR: argument not allowed in function normal\n");
		  return realvar(nrobs);
		  }
		}   // end: normal
	 else if (functionname == "exponential")
		{
        bool bracketmiss;
		arglist = argument.strtoken2(",",bracketmiss);

        if (bracketmiss)
          {
          errormessages.push_back("ERROR: missing bracket(s)\n");
		  return realvar(nrobs);
          }

		if (arglist.size() != 1)
		  {
		  errormessages.push_back("ERROR: invalid number of arguments for function 'exponential'\n");
		  return realvar(nrobs);
		  }
		else
          {
   		  valuevek = eval_exp(arglist[0],false);
#if defined(JAVA_OUTPUT_WINDOW)
          realvar h(nrobs);
          register unsigned i;
          realvar::iterator it = h.begin();

          jmethodID javaexponential = adminb_p->Java->GetMethodID(
          adminb_p->BayesX_cls, "jexponential", "(D)D");

          for (i=0;i<nrobs;i++,++it)
            {
            if(valuevek[i]<=0 || valuevek[i]==NA)
              *it =  NA;
            else
              *it = adminb_p->Java->CallDoubleMethod(adminb_p->BayesX_obj, javaexponential, valuevek[i].getvalue());
            }
          return h;
#else
		  return exponential(valuevek);
#endif
          }
		}   // end: exponential
	 else if (functionname == "poisson")
		{
        bool bracketmiss;
		arglist = argument.strtoken2(",",bracketmiss);

        if (bracketmiss)
          {
          errormessages.push_back("ERROR: missing bracket(s)\n");
		  return realvar(nrobs);
          }

		if (arglist.size() != 1)
		  {
		  errormessages.push_back("ERROR: invalid number of arguments for function 'poisson'\n");
		  return realvar(nrobs);
		  }
		else
          {
   		  valuevek = eval_exp(arglist[0],false);
#if defined(JAVA_OUTPUT_WINDOW)
          realvar h(nrobs);
          register unsigned i;
          realvar::iterator it = h.begin();

          jmethodID javapoisson = adminb_p->Java->GetMethodID(
          adminb_p->BayesX_cls, "jpoisson", "(D)D");

          for (i=0;i<nrobs;i++,++it)
            {
            if(valuevek[i]<=0 || valuevek[i]==NA)
              *it = NA;
            else
        	  *it = adminb_p->Java->CallDoubleMethod(adminb_p->BayesX_obj, javapoisson, valuevek[i].getvalue());
            }
          return h;
#else
		  errormessages.push_back("ERROR: Function 'poisson' not available in this version.\n");
		  return realvar(nrobs);
#endif
          }
  		}   // end: poisson
	 else if (functionname == "weibull")
		{
        bool bracketmiss;
		arglist = argument.strtoken2(",",bracketmiss);

        if (bracketmiss)
          {
          errormessages.push_back("ERROR: missing bracket(s)\n");
		  return realvar(nrobs);
          }

		if (arglist.size() != 2)
		  {
		  errormessages.push_back("ERROR: invalid number of arguments for function 'weibull'\n");
		  return realvar(nrobs);
		  }
		else
          {
		  valuevek = eval_exp(arglist[0],false);
		  valuevek2 = eval_exp(arglist[1],false);
#if defined(JAVA_OUTPUT_WINDOW)
          realvar h(nrobs);
          register unsigned i;
          realvar::iterator it = h.begin();

          jmethodID javaweibull = adminb_p->Java->GetMethodID(
          adminb_p->BayesX_cls, "jweibull", "(DD)D");

          for (i=0;i<nrobs;i++,++it)
            {
            if(valuevek[i]<=0 || valuevek[i]==NA || valuevek2[i]<=0 || valuevek2[i]==NA)
              *it =  NA;
            else
              *it = adminb_p->Java->CallDoubleMethod(adminb_p->BayesX_obj, javaweibull, valuevek[i].getvalue(), valuevek2[i].getvalue());
            }
          return h;
#else
		  errormessages.push_back("ERROR: Function 'weibull' not available in this version.\n");
		  return realvar(nrobs);
#endif
          }
  		}   // end: weibull
	 else if (functionname == "bernoulli")
		{
        bool bracketmiss;
		arglist = argument.strtoken2(",",bracketmiss);

        if (bracketmiss)
          {
          errormessages.push_back("ERROR: missing bracket(s)\n");
		  return realvar(nrobs);
          }

		if (arglist.size() != 1)
		  {
		  errormessages.push_back("ERROR: invalid number of arguments for function 'bernoulli'\n");
		  return realvar(nrobs);
		  }
		else
          {
          valuevek = eval_exp(arglist[0],false);
#if defined(JAVA_OUTPUT_WINDOW)
          realvar h(nrobs);
          register unsigned i;
          realvar::iterator it = h.begin();

          jmethodID javabernoulli = adminb_p->Java->GetMethodID(
          adminb_p->BayesX_cls, "jbernoulli", "(D)D");

          for (i=0;i<nrobs;i++,++it)
            {
            if(valuevek[i]<0 || valuevek[i]>1 || valuevek[i]==NA)
              *it =  NA;
            else
              *it = adminb_p->Java->CallDoubleMethod(adminb_p->BayesX_obj, javabernoulli, valuevek[i].getvalue());
            }
          return h;
#else
          return bernoulli(valuevek);
#endif
          }
		}   // end: bernoulli
	 else if (functionname == "binomial")
		{
        bool bracketmiss;
		arglist = argument.strtoken2(",",bracketmiss);

        if (bracketmiss)
          {
          errormessages.push_back("ERROR: missing bracket(s)\n");
		  return realvar(nrobs);
          }

		if (arglist.size() != 2)
		  {
		  errormessages.push_back("ERROR: invalid number of arguments for function 'binomial'\n");
		  return realvar(nrobs);
		  }
		else
		  {
		  valuevek = eval_exp(arglist[0],false);
		  valuevek2 = eval_exp(arglist[1],false);
#if defined(JAVA_OUTPUT_WINDOW)
          realvar h(nrobs);
          register unsigned i;
          realvar::iterator it = h.begin();

          jmethodID javabinomial = adminb_p->Java->GetMethodID(
          adminb_p->BayesX_cls, "jbinomial", "(DD)D");

          for (i=0;i<nrobs;i++,++it)
            {
            if(valuevek2[i]<0 || valuevek2[i]>1 || valuevek2[i]==NA
                              || valuevek[i]<1 || valuevek[i]==NA)
              *it = NA;
            else
      	      *it = adminb_p->Java->CallDoubleMethod(adminb_p->BayesX_obj, javabinomial, valuevek[i].getvalue(), valuevek2[i].getvalue());
            }
          return h;
#else
		  return binomial(valuevek,valuevek2);
#endif
		  }
		}   // end: binomial
      else if (functionname == "gamma")
        {

        bool bracketmiss;
		arglist = argument.strtoken2(",",bracketmiss);

        if (bracketmiss)
          {
          errormessages.push_back("ERROR: missing bracket(s)\n");
		  return realvar(nrobs);
          }

		if (arglist.size() != 2)
		  {
		  errormessages.push_back("ERROR: invalid number of arguments for function 'gamma'\n");
		  return realvar(nrobs);
		  }
		else
		  {
		  valuevek = eval_exp(arglist[0],false);
		  valuevek2 = eval_exp(arglist[1],false);
#if defined(JAVA_OUTPUT_WINDOW)
          realvar h(nrobs);
          register unsigned i;
          realvar::iterator it = h.begin();

          jmethodID javagamma = adminb_p->Java->GetMethodID(
          adminb_p->BayesX_cls, "jgamma", "(DD)D");

          for (i=0;i<nrobs;i++,++it)
            {
            if(valuevek[i]<=0 || valuevek[i]==NA || valuevek2[i]<=0 || valuevek2[i]==NA)
              *it = NA;
            else
      	      *it = adminb_p->Java->CallDoubleMethod(adminb_p->BayesX_obj, javagamma, valuevek[i].getvalue(), valuevek2[i].getvalue());
            }
          return h;
#else
		  return gamma(valuevek,valuevek2);
#endif
		  }

        }
	 else
		{
		errormessages.push_back("ERROR: " + functionname + " unknown function\n");
		return realvar(nrobs);
		}
	 }  // end: functions
  else
	 {
	 errormessages.push_back(
	 "ERROR: expression syntax error in \"" + expression + "\"\n");
	 return realvar(nrobs);
	 }

  }


void dataset::sort(list<realvar>::iterator varit,int left,int right)
  {
  int hilfe;
  int i = left;
  int j = right;
  realobs x = (*varit)[datarep.index[(left+right)/2]];
  do
	 {
	 while ((*varit)[datarep.index[i]] < x)
		i++;
	 while (x < (*varit)[datarep.index[j]])
	 j--;
	 if (i <= j)
		{
		hilfe = datarep.index[i];
		datarep.index[i] = datarep.index[j];
		datarep.index[j] = hilfe;
		i++;
		j--;
		}
	 }
  while (i <= j);

  if (left < j)
	 sort(varit,left,j);
  if (i < right)
	 sort(varit,i,right);
  }


void dataset::sort (list<ST::string> names,int left,int right)
  {
  errormessages.clear();

  list<ST::string>::iterator i;

  for (i=names.begin();i!=names.end();++i)
	 if (datarep.findvar(*i) == 1)
		 errormessages.push_back(
		 "ERROR: variable " + (*i) + " can not be found\n");


  if (errormessages.empty())
	 {
	 list<realvar>::iterator varit;
	 datarep.findvar(*(names.begin()),i,varit);

	 sort(varit,left,right);

	 if (names.size() > 1)
		{
		names.pop_front();
		int beg = left;
		int j = left;
		realobs pred = (*varit)[datarep.index[left]];
		while (j <= right)
		  {
		  if ((*varit)[datarep.index[j]] != pred)
			 {
			 sort(names,beg,j-1);
			 beg = j;
			 }
		  else if (j == right)
			 sort(names,beg,j);
		  pred = (*varit)[datarep.index[j]];
		  j++;
		  }
		}
	 }

  }


void dataset::rename(const ST::string & oldname, const ST::string & newname)
  {
  errormessages.clear();

  list<ST::string>::iterator stit;
  list<realvar>::iterator varit;
  list<ST::string>::iterator stit2;
  list<realvar>::iterator varit2;


  if (datarep.findvar(oldname,stit,varit) == 0)
	 {
	 if (datarep.findvar(newname,stit2,varit2) != 0)
		{
		if (newname.isvarname() == 0)
		  (*stit) = newname;
		else
		  errormessages.push_back("ERROR: " + newname + " invalid varname\n");
		}
	 else
		errormessages.push_back(
		"ERROR: " + newname + " is already existing\n");


	 }
  else
	 errormessages.push_back(
	 "ERROR: variable " + oldname + " can not be found\n");
  }


void dataset::addvariable(const ST::string & name, const realvar & v)
  {

  errormessages.clear();

  if (datarep.findvar(name) == 0)
	 errormessages.push_back("ERROR: variable " + name + " is already existing\n");
  else if (name.isvarname() != 0)
	 errormessages.push_back("ERROR: invalid variable name specification\n");
  else
	 {
	 datarep.empty = false;
	 datarep.varnames.push_back(name);
	 datarep.variables.push_back(v);
	 }

  }


void dataset::addvariable(const ST::string & name, const ST::string & e)
  {

  realvar v = eval_exp(e);

  if (errormessages.empty())
	 addvariable(name,v);

  }


unsigned dataset::replace(const ST::string & name, const realvar & v)
  {

  errormessages.clear();

  list<ST::string>::iterator stit;
  list<realvar>::iterator varit;

  if (datarep.findvar(name,stit,varit) == 0)
    {
    *varit = v;
    return nrobs;
    }
  else
    {
    errormessages.push_back("ERROR: variable " + name + "not found\n");
    return 0;
    }


  }


unsigned dataset::replace(const ST::string & name,const ST::string & e,const ST::string & boole)
  {

  unsigned changed = 0;
  errormessages.clear();

  list<ST::string>::iterator stit;
  list<realvar>::iterator varit;

  realvar vneu;
  realvar vbool;
  realvar vboolneg;


  if (datarep.findvar(name,stit,varit) != 0)
	 errormessages.push_back("ERROR: variable " + name + " not found\n");
  else
	 vneu = eval_exp(e);

  if (errormessages.empty() && (boole.length() != 0))
	 vbool = eval_exp(boole);

  if (errormessages.empty())
	 if (boole.length() != 0)
		{
        unsigned i;
        realvar::iterator boolit = vbool.begin();
        for (i=0;i<nrobs;i++,++boolit)
          {
          if (*boolit == 1)
            {
            (*varit)[i] = vneu[i];
            changed++;
            }
          }

		}
	 else
       {
       *varit = vneu;
       changed = nrobs;
       }

  return changed;
  }


void dataset::dropvariable(const ST::string & name)
  {
  errormessages.clear();

  list<ST::string>::iterator stit;
  list<realvar>::iterator varit;

  if (datarep.findvar(name,stit,varit) == 0)
	 {
	 datarep.varnames.erase(stit);
	 (*varit).clear();
	 datarep.variables.erase(varit);
	 if (datarep.variables.empty())
		datarep.clear();
	 }
  else
	 errormessages.push_back("ERROR: variable " + name + " can not be found\n");
  }


void dataset::dropvariables(list<ST::string> & names)
  {
  errormessages.clear();

  if (names.empty())
	 names = datarep.varnames;

  list<ST::string>::iterator i;

  for (i=names.begin();i != names.end();++i)
	 if (datarep.findvar(*i) == 1)
	  errormessages.push_back("ERROR: variable " + (*i) + " can not be found\n");

  i=names.begin();
  while ( (i != names.end()) && (errormessages.empty()) )
	 {
	 dropvariable(*i);
	 ++i;
	 }
  }


unsigned dataset::dropobservations(const ST::string & e)
  {
  errormessages.clear();
  unsigned nrelim = 0;

  if (e.length() != 0)
	 {

	 realvar vbool;
	 vbool = eval_exp(e);

	 list<realvar>::iterator varit;

//	 long i,k,anf,end;

//	 anf = 0;
//	 end = nrobs-1;

     unsigned i;
     unsigned pos=0;

     statmatrix<int> rank(nrobs,1);
     vector<int>::iterator indit = datarep.index.begin();
     for(i=0;i<nrobs;i++,++indit)
       rank(*indit,0) = i;


     vector<realobs>::iterator vboolit=vbool.begin();
     int* workrank_i = rank.getV();
     int* workrank_pos = rank.getV();
	 for(i=0;i<nrobs;i++,++vboolit,workrank_i++)
		{

        if (*vboolit == 0)
          {
          if (pos < i)
            {
            for (varit = datarep.variables.begin();varit!=datarep.variables.end();varit++)
              {
              (*varit)[pos] = (*varit)[i];
			  }
            *workrank_pos = *workrank_i;
            }
          pos++;
          workrank_pos++;
          }
        else
          nrelim++;

		} // end: for(i=0;i<nrobs;i++)

	 for (varit = datarep.variables.begin();varit!=datarep.variables.end();varit++)
       (*varit).erase((*varit).begin()+nrobs-nrelim,(*varit).end());

     datarep.index.erase(datarep.index.begin()+nrobs-nrelim,datarep.index.end());

	 nrobs-=nrelim;

     // repairing the index:

     if (nrobs > 0)
       {
       statmatrix<int> ind(nrobs,1);
       ind.indexinit();

       rank.indexsort(ind,0,nrobs-1,0,0);

       int * workind = ind.getV();
       indit = datarep.index.begin();
       for(i=0;i<nrobs;i++,workind++,++indit)
         *indit = *workind;
       }

     return nrelim;

	 } // end:  if (e.length() != 0)

  return 0;
  }


void dataset::setobs(unsigned newobs)
  {
  errormessages.clear();
  if (newobs < nrobs)
     {
	 errormessages.push_back("ERROR: new number of observations must be greater\n");
     errormessages.push_back("       than current number of observations\n");
     }
  else
	 {

	 unsigned i;
	 list<realvar>::iterator varit;
	 for (varit = datarep.variables.begin();varit != datarep.variables.end();
			++varit)
        {

        (*varit).resize(newobs,NA);
//        (*varit).reserve(newobs);

//		for (i=nrobs;i<newobs;i++)
//		  (*varit).push_back(NA);

        }

     datarep.index.reserve(newobs);
	 for (i=nrobs;i<newobs;i++)
		datarep.index.push_back(i);

	 nrobs = newobs;
	 }
  }


void dataset::makematrix(list<ST::string> & names,datamatrix & d,ST::string boole)
  {

  errormessages.clear();

  if (names.empty())
	 names = datarep.varnames;

  if (names.size() == 0)
    errormessages.push_back("ERROR: dataset contains no variables\n");

  list<ST::string>::iterator i;
  for (i=names.begin();i!=names.end();++i)
//	 if (*i != "const")
		if (datarep.findvar(*i) == 1)
		  errormessages.push_back(
		  "ERROR: variable " + (*i) + " can not be found\n");

  if (errormessages.empty())
	 {

	 f.filterNA(datarep,names);

	 realvar vbool;
	 if (boole.length() > 0)
		{
		vbool = eval_exp(boole);
		f = f + vbool;
		}

     if (f.nrunfiltered() == 0)
       {
       errormessages.push_back("ERROR: no valid (nonmissing) observations\n");
       }

     if (errormessages.empty())
       {
	   d = datamatrix(f.nrunfiltered(),names.size());
	   unsigned col = 0;
	   unsigned row;
	   unsigned j;
       double * dp;


	   for (i=names.begin();i != names.end();++i)
		 {
/*
		 if ( (*i) == "const")
           {
           dp = d.getV()+col;
           for(k=0;k<d.rows();k++,dp+=d.cols())
             *dp = 1;
           }
		 else
		   {
*/
		   list<ST::string>::iterator stit;
		   list<realvar>::iterator vit;
		   datarep.findvar(*i,stit,vit);
		   row  = 0;
           vector<int>::iterator indexp = datarep.index.begin();
           dp = d.getV()+col;
		   for (j = 0;j<datarep.obs();j++,++indexp)
			 if (f[*indexp] == 1)
               {
               *dp = (*vit)[*indexp].getvalue();
               row++;
               dp+=d.cols();
               }
//           }
         col++;
         } // end: for (i=names.begin();i != names.end();++i)
       }
	 } // end: if (errormessages.empty())

  }


void dataset::makematrix(vector<ST::string> & names, datamatrix & d,
								 ST::string boole)
  {
  list<ST::string> namesl(names.begin(),names.end());
  makematrix(namesl,d,boole);
  }


void dataset::makematrix(ST::string & name, datamatrix & d, ST::string boole)
  {
  list<ST::string> namesl;
  namesl.push_back(name);
  makematrix(namesl,d,boole);
  }


void dataset::reverseorder(void)
  {
  register unsigned i;
  vector<int> help = datarep.index;
  vector<int>::iterator hit = help.end()-1;
  vector<int>::iterator it = datarep.index.begin();

  for (i=0;i<nrobs;i++,--hit,++it)
    *it = *hit;

  }


void dataset::marketing(vector<ST::string> & names, ST::string & defs, int & lak, double & alph)
  {

  list<ST::string>::iterator stit;
  list<realvar>::iterator varitwo;
  list<realvar>::iterator varitmp;
  list<realvar>::iterator varitp;
  list<realvar>::iterator varitm;
  list<realvar>::iterator regit;

//  ofstream out("c:\\cprog\\test\\osaft\\marke.raw");  //File nur zum Überprüfen

  int nrout=0;
  int nrwoche=0;

  unsigned i;
  for(i=0;i<names.size()-2;i++)    //durchläuft die Variablen outlet, wochenin
    {
    list<ST::string>namen;
    namen.push_back(names[i]);
    int found = datarep.findvar(names[i],stit,varitwo);  //damit man varit bekommt und später Werte auslesen kann
    sort(namen,0,nrobs-1);   //sortiert nach Variable "i"
    unsigned k = 0;      //Anzahl der Beobachtungen der einzelnen Variablen feststellen
    int zaehler = 0;
    while(k<datarep.index.size())
      {
      int anz = 0;
      realobs q = getvalue(k,varitwo);
      for(unsigned j=0;j<nrobs;j++)
         {
         realobs p = getvalue(j,varitwo);
         if (p == q)
           anz = anz+1;
         }
         k = k + anz;
         zaehler = zaehler + 1;
      }
    if(i==0)
      nrout = zaehler;
    else
      nrwoche = zaehler;
    }

  list<ST::string> namen;  //zum Sortieren
  namen.push_back(names[2]);
  namen.push_back(names[0]);
  namen.push_back(names[1]);
  namen.push_back(names[3]);
  sort(namen,0,nrobs-1);   // endgültiges Sortieren; Reihenfolge: markenin, outlet, wochenin, preis
  realvar wert(nrobs,NA);      //erzeugt neue Variable, die später den minimalen Preis pro Marke enthält

  int found = datarep.findvar(names[3],stit,varitp);    //Zeiger auf "preis"
  found = datarep.findvar(names[1],stit,varitwo);       //Zeiger auf "wochenin"
  found = datarep.findvar(names[2],stit,varitm);        //Zeiger auf "markenin"

  int w = 0;
  realobs wo = getvalue(0,varitwo);  //liefert die einzelnen Werte als "realobs"
  realobs min = getvalue(0,varitp);
  realobs ma = getvalue(0,varitm);
  vector<int> anzpreis;  //um zu bekommen, wie oft jede Marke pro Woche u. Geschäft vorkommt

  for(i=0;i<nrobs;i++)    //Schleife, die Werte zum "Weiterhüpfen" liefert
    {
    if(getvalue(i,varitm) != ma)
       {
       anzpreis.push_back(w);
       ma = getvalue(i,varitm);
       }
    if(getvalue(i,varitwo) != wo)
       {
       wo = getvalue(i,varitwo);
       min = getvalue(i,varitp);
       w = 1;
       }
    else
       {
         {
         w = w + 1;
         }
       }
    wert[datarep.index[i]] = min;
    }
  ST::string minpreis = "minpreis";
  addvariable(minpreis, wert); //fügt neue Variable hinzu
  found = datarep.findvar(minpreis,stit,varitmp);  //Zeiger auf "minpreis"

  double sum = 0;
  for(i=0;i<anzpreis.size();i++)    //berechnet Anzahl der letzten Marke pro Geschäft und Woche
    {
    sum = sum + nrwoche*nrout*anzpreis[i];
    }
  w = (nrobs-sum) / (nrwoche*nrout);
  anzpreis.push_back(w);            //Vektor "anzpreis" ist fertig

  unsigned zma = 0;
  int zout = 1;
  double oben = nrwoche*anzpreis[zma];
  double unten = 0;
  vector<realobs> armit;
  if(defs=="priceindex")             //Durchschnittspreis einer Marke pro Geschäft
    {                                //Werte sind in Vektor "armit"
    while(oben<=nrobs && zma<anzpreis.size())
       {
       realobs sum = 0;
       double nenner = oben - unten;
       for(unsigned y=unten;y<oben;y++)  //Schleife durchläuft die Werte einer Marke
         {
         if(getvalue(y,varitmp)==NA)
           nenner = nenner - 1;
         }
       for(unsigned y=unten;y<oben;y++)  //Schleife durchläuft die Werte einer Marke
         {
         if(getvalue(y,varitmp)!=NA)
           sum = sum + getvalue(y,varitmp) / nenner;
         }
       armit.push_back(sum);
       unten = oben;
       zout = zout + 1;
       if(zout>nrout)
         {
         zma = zma + 1;
         zout = 1;
         }
       oben = oben + nrwoche*anzpreis[zma];
       }
    }

  if(defs=="regular")
    {
    realvar wert(nrobs,NA);      //erzeugt neue Variable, die später den regulären Preis einer Marke enthält
    zma = 0;
    realobs reg = getvalue(0,varitmp);
    realobs akt1 = getvalue(0,varitmp);
    realobs akt2 = getvalue(0,varitmp);
    realobs regm;
    realobs aktm;
    realvar eintrag(nrobs,NA);
    int j = 0;
    int jalt = 0;
    int wo = 0;
    vector<int> erstreg;
    vector<realobs> erp;

    while(j<nrobs)          // Initialisierung
      {
      while(akt1==NA && j < jalt + anzpreis[zma]*nrwoche)
        {
        j = j + 1;
        akt1 = getvalue(j,varitmp);
        }
      if(j == jalt + anzpreis[zma]*nrwoche)
        {
        erstreg.push_back(jalt);
        erp.push_back(NA);
        }
      else
        {
        unsigned jhilf = j;
        while(jhilf < jalt + anzpreis[zma]*nrwoche - anzpreis[zma])
          {
          akt2 = getvalue(j+anzpreis[zma],varitmp);
          if(akt2<(1-alph)*akt1 && akt2!=NA)
            {
            erp.push_back(akt1);
            erstreg.push_back(j+anzpreis[zma]);
            jhilf = jalt + anzpreis[zma]*nrwoche;
            }
          else if(akt2>(1-alph)*akt1 && akt2!=NA)
            {
            erp.push_back(akt2);
            erstreg.push_back(j+anzpreis[zma]);
            jhilf = jalt + anzpreis[zma]*nrwoche;
            }
          else if(akt2==NA)
            {
            jhilf = jhilf + anzpreis[zma];
            j = j + anzpreis[zma];
            }
          else
            {
            jhilf = jhilf + anzpreis[zma];
            j = j + anzpreis[zma];
            akt1 = akt2;
            }
          }
        if(j == jalt + anzpreis[zma]*nrwoche - anzpreis[zma])
          {
          erstreg.push_back(jalt);
          erp.push_back(getvalue(0,varitmp));
          }
        }
      j = jalt + anzpreis[zma]*nrwoche;
      jalt = j;
      if(j < nrobs)
        akt1 = getvalue(j,varitmp);
      wo = wo + anzpreis[zma]*nrwoche;
      if(wo>=nrwoche*nrout*anzpreis[zma])
        {
        zma = zma + 1;
        wo = 0;
        }
      }

    wo = 0;
    zma = 0;
    j = 0;
    int y = 0;
    jalt = 0;
    for(i=0;i<erstreg.size();i++)          // Eintragen der Werte
      {
      reg = erp[i];
      akt1 = getvalue(erstreg[i],varitmp);
      regm = reg;
      aktm = akt1;
      for(int k=erstreg[i];k<(erstreg[i]+anzpreis[zma]);k++)
        {
        wert[datarep.index[k]] = reg;
        }
      j = erstreg[i] + anzpreis[zma];
      y = erstreg[i] - 1;
      while(j < jalt + nrwoche*anzpreis[zma])   // Eintragen der Werte nach "oben"!!!
        {
        akt2 = getvalue(j,varitmp);
        realobs reg2;
        if(akt1==reg && akt1!=NA)      //NA hat immer den größten Wert, deshalb Sonderfall!!
          {
          if(akt2!=NA)
            {
            aktm = akt2;
            if(akt2>=(1-alph)*akt1)
              {
              reg = akt2;
              regm = reg;
              }
            }
          else
            reg = NA;
          }
        else if(akt1==NA)       //wenn alter Wert fehlend ist
          {
          if(akt2!=NA)
            {
            if(aktm==regm)
              {
              if(akt2>=(1-alph)*aktm)
                reg = akt2;
              else
                reg = regm;
              }
            else
              {
              if(akt2>(1+alph)*aktm)
                {
                reg2 = akt2;
                if(reg2>=(1-alph)*regm)
                  reg = reg2;
                else
                  reg = regm;
                }
              else
                reg = regm;
              }
            aktm = akt2;
            regm = reg;
            }
          else
            reg = NA;
          }
        else                   //Bedingung 2 (AP_t < RP_t)
          {
          if(akt2!=NA)
            {
            aktm = akt2;
            if(akt2>(1+alph)*akt1)
              {
              //reg = akt2;            // wenn innere Bedingung weggelassen wird
              reg2 = akt2;
              //if(reg2>=(1-alph)*reg)      // Original-Bed.
              //if(reg2<=(1-alph)*reg)      // umgedrehtes Zeichen
              if(reg2<=(1-alph)*reg || reg2>=reg)
                reg = reg2;
              }
            regm = reg;
            }
          else
            reg = NA;
          }
        for(int k=j;k<(j+anzpreis[zma]);k++)
          {
          wert[datarep.index[k]] = reg;
          }
        j = j + anzpreis[zma];
        akt1 = akt2;
        }
      reg = erp[i];
      akt1 = getvalue(erstreg[i],varitmp);
      regm = reg;
      aktm = akt1;
      while(y >= jalt)              // Eintragen der Werte nach "unten"!!!
        {
        akt2 = getvalue(y,varitmp);
        realobs reg2;
        if(akt1==reg && akt1!=NA)      //NA hat immer den größten Wert, deshalb Sonderfall!!
          {
          if(akt2!=NA)
            {
            aktm = akt2;
            if(akt2>=(1-alph)*akt1)
              {
              reg = akt2;
              regm = reg;
              }
            }
          else
            reg = NA;
          }
        else if(akt1==NA)       //wenn alter Wert fehlend ist
          {
          if(akt2!=NA)
            {
            if(aktm==regm)
              {
              if(akt2>=(1-alph)*aktm)
                reg = akt2;
              else
                reg = regm;
              }
            else
              {
              if(akt2>(1+alph)*aktm)
                {
                reg2 = akt2;
                if(reg2>=(1-alph)*regm)
                  reg = reg2;
                else
                  reg = regm;
                }
              else
                reg = regm;
              }
            aktm = akt2;
            regm = reg;
            }
          else
            reg = NA;
          }
        else                   //Bedingung 2 (AP_t < RP_t)
          {
          if(akt2!=NA)
            {
            aktm = akt2;
            if(akt2>(1+alph)*akt1)
              {
              //reg = akt2;
              reg2 = akt2;
              //if(reg2>=(1-alph)*reg)
              //if(reg2<=(1-alph)*reg)
              if(reg2<=(1-alph)*reg || reg2>=reg)
                reg = reg2;
              }
            regm = reg;
            }
          else
            reg = NA;
          }
          for(int k=y;k>(y-anzpreis[zma]);k--)
            {
            wert[datarep.index[k]] = reg;
            }
        y = y - anzpreis[zma];
        akt1 = akt2;
        }
      jalt = jalt + nrwoche*anzpreis[zma];
      wo = wo + anzpreis[zma]*nrwoche;
      if(wo>=nrwoche*nrout*anzpreis[zma])
        {
        zma = zma + 1;
        wo = 0;
        }
      }
    ST::string regpreis = "regpreis";
    addvariable(regpreis, wert); //fügt neue Variable hinzu
    found = datarep.findvar(regpreis,stit,regit);  //Zeiger auf "regpreis"
    for(i=0;i<nrobs;i++)         //berechnet den Quotienten "aktueller Preis / regulärer Preis"
      {
      if(getvalue(i,regit)!=NA)
        eintrag[datarep.index[i]] = getvalue(i,varitmp) / getvalue(i,regit);
      else
       eintrag[datarep.index[i]] = NA;
      }
    ST::string regquot = "regquot";
    addvariable(regquot, eintrag); //fügt neue Variable hinzu
    found = datarep.findvar(regquot,stit,regit);  //Zeiger auf "regquot"
    }

  unsigned z1 = 0;
  zma = 0;
  unsigned zko = 0;
  oben = nrwoche*nrout*anzpreis[zma]; //festlegen der obersten Beob., bis zu der Werte eingtrag. werden sollen
  unten = 0;    //unterste Beob., ab der Werte eingetragen werden sollen
  int l = 0;
  i = 0;
  while(i<nrobs)         //Schleife, die die Variablen erzeugt
     {
     realobs ma = getvalue(i,varitm);
     double marke = ma.getvalue(); //wandelt die "realobs"-Werte in "double" um
     realvar vorwo(nrobs,NA);   //für die lag-Variable
     unsigned j = 0;
     unsigned z2 = 0;
     int mitind = 0;
     while(j<nrobs)
        {
        realobs ko = getvalue(j,varitm);   //liefert die einzelnen Werte als "realobs"
        double konk = ko.getvalue();  //wandelt die "realobs"-Werte in "double" um
        realvar neuvar(nrobs,NA);     //für die Preisindizes
        j = j + nrwoche*nrout*anzpreis[z2];
        z2 = z2 + 1;
        int xy = 1;
        int k = 0;
        int mit = 0;
        int zaehlermarke = 0;
        for(unsigned y=unten;y<oben;y++)  //Schleife durchläuft die Werte der Marke und weist die Beobachtungen zu
           {
           if(anzpreis[zko]==anzpreis[zma])      //Anzahl pro Woche und Geschäft bei Marke und Konkurrenz gleich
             {
             if(defs=="notransform")
               neuvar[datarep.index[y]] = getvalue(y+l,varitmp);
             else if(defs=="priceindex")
               {
               if(getvalue(y+l,varitmp)==NA)
                 neuvar[datarep.index[y]] = NA;
               else
                 neuvar[datarep.index[y]] = getvalue(y+l,varitmp) / armit[mitind];
               mit = mit + 1;
               if(mit==nrwoche*anzpreis[zma])
                 {
                 mitind = mitind + 1;
                 mit = 0;
                 }
               }
             else
               neuvar[datarep.index[y]] = getvalue(y+l,regit);
             }
           else if(anzpreis[zma]<anzpreis[zko])  //Anzahl pro Woche und Geschäft bei Marke kleiner als bei Konkurrenz
             {
             zaehlermarke = zaehlermarke + 1;
             if(defs=="notransform")
               neuvar[datarep.index[y]] = getvalue(y+l+k,varitmp);
             else if(defs=="priceindex")
               {
               if(getvalue(y+l+k,varitmp)==NA)
                 neuvar[datarep.index[y]] = NA;
               else
                 neuvar[datarep.index[y]] = getvalue(y+l+k,varitmp) / armit[mitind];
               mit = mit + 1;
               if(mit==nrwoche*anzpreis[zma])
                 {
                 mitind = mitind + 1;
                 mit = 0;
                 }
               }
             else
               neuvar[datarep.index[y]] = getvalue(y+l+k,regit);
             if(zaehlermarke>=anzpreis[zma])
             {
             zaehlermarke = 0;
             k = k + anzpreis[zko] - anzpreis[zma];  //k sorgt dafür, daß die überflüssigen Beobachtungen übersprungen werden
             }
             }
           else                  //Anzahl pro Woche und Geschäft bei Marke größer als bei Konkurrenz
             {
             if(defs=="notransform")
               neuvar[datarep.index[y]] = getvalue(y+l+k,varitmp);
             else if(defs=="priceindex")
               {
               if(getvalue(y+l+k,varitmp)==NA)
                 neuvar[datarep.index[y]] = NA;
               else
                 neuvar[datarep.index[y]] = getvalue(y+l+k,varitmp) / armit[mitind];
               mit = mit + 1;
               if(mit==nrwoche*anzpreis[zma])
                 {
                 mitind = mitind + 1;
                 mit = 0;
                 }
               }
             else
               neuvar[datarep.index[y]] = getvalue(y+l+k,regit);
             if(xy<anzpreis[zma])
               {
               k = k - 1;    //k sorgt dafür, daß die Beobachtungen immer mehrmals übernommen werden
               xy = xy + 1;  //xy zählt mit, wie oft ein Wert genommen wird
               }
             else
               {
               xy = 1;
               k = k + anzpreis[zko] - 1;
               }
             }
           if(getvalue(y,varitwo)-getvalue(0,varitwo) >= lak && zma==zko)  //für Spalte mit Werten einer Vorwoche
             {
             int stelle = y-(anzpreis[zma]*lak);
                vorwo[datarep.index[y]] = neuvar[datarep.index[stelle]];
             }
           }
        if(zko<anzpreis.size()-1)
          {
          l = l + nrwoche*nrout*anzpreis[zko];  //berechnet das neue "l" zum Weiterspringen
          zko = zko + 1;
          }
        else               //die Werte für die nächste Marke werden angegeben
          {
          unten = oben;
          zma = zma + 1;
          oben = oben + nrwoche*nrout*anzpreis[zma];
          zko = 0;
          l = 0 - unten;
          }
        ST::string neu;
        if(defs=="notransform")
           neu = "price_" + ST::doubletostring(marke,0)
             + "_" + ST::doubletostring(konk,0);
        else if(defs=="priceindex")
           neu = "price_ind_" + ST::doubletostring(marke,0)
             + "_" + ST::doubletostring(konk,0);
        else
           neu = "price_regquot_" + ST::doubletostring(marke,0)
             + "_" + ST::doubletostring(konk,0);
        addvariable(neu,neuvar);
        }                                      //Ende innere while-Schleife
     i = i + nrwoche*nrout*anzpreis[z1];
     z1 = z1 + 1;
     ST::string vor;
     if(defs=="notransform")
       vor = "price" + ST::doubletostring(marke,0) + "_weekbef" + ST::inttostring(lak);
     else if(defs=="priceindex")
       vor = "price_ind" + ST::doubletostring(marke,0) + "_weekbef" + ST::inttostring(lak);
     else
       vor = "price_regquot" + ST::doubletostring(marke,0) + "_weekbef" + ST::inttostring(lak);
     addvariable(vor,vorwo);
     }                                        //Ende äußere while-Schleife
  }




#if defined(BORLAND_OUTPUT_WINDOW)
//---------------------------------------------------------------------------
#pragma package(smart_init)
#endif




