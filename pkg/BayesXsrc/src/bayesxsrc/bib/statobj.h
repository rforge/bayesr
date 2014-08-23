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



#if !defined (STATOBJECT_INCLUDED)

#define STATOBJECT_INCLUDED

#include"../export_type.h"
#include<vector>
#include"clstring.h"
#include"option.h"
#include"command.h"

#if defined(JAVA_OUTPUT_WINDOW)
#include "adminparse_basic.h"
#endif


class __EXPORT_TYPE statobject
  {

  protected:

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adminb_p;
  #endif

  //------------------------- PROTECTED VARIABLES ------------------------------

  // name of the statobject

  ST::string name;

  // type of the statobject

  ST::string type;

  // contains methods

  vector<command> methods;

  // contains current errormessages
  // errormessages will be deleted before (re-)parsing

  vector<ST::string> errormessages;

  // contains new commands
  // to be executed after parse

  vector<ST::string> newcommands;


  // contains text for method describe

  vector<ST::string> describetext;

  // stream object to log-file

  ofstream * logout;

  // stream object to input file

  istream * input;

  // defaultpath for output

  ST::string defaultpath;

  // FUNCTION: parsecommands
  // TASK: parses command 'c'
  //       returns: -2 if an error occured
  //                -1 if in 'c' a 'globaloption'
  //                   method 'describe' is specified
  //                if 'c' specifies a method, the position of method in
  //                'methods'

  int parsecom(const ST::string & c, vector<command> & methods,
					optionlist & globaloptions);// = optionlist());

  // FUNCTIONS: out
  // TASK: writes the string 'c' or strings 'm' to cout and logout


  void out(const ST::string & c,bool thick=false,bool italic = false,
           unsigned size = 12,int r=0, int g=0, int b=0, bool descr=false);

  void out(const vector<ST::string> & m,bool thick=false,bool italic = false,
           unsigned size = 12,int r=0, int g=0, int b=0,bool descr=false);

  // FUNCTIONS: out
  // TASK: writes string 'c' or 'strings 'c' to the output windows
  //       uses function out with thick=true and italic = true

  void outerror(const ST::string & c);

  void outerror(const vector<ST::string> & c);

  public:


  //--------------------------- PUBLIC FUNCTIONS -------------------------------

  // DEFAULT CONSTRUCTOR

  statobject(void) {}

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - n = name
  // - t = type
  // - logout = lo
  // - input = in
  // - defaultpath = p

  #if defined(JAVA_OUTPUT_WINDOW)
  statobject(administrator_basic * adb,const ST::string & n,
             const ST::string t,ofstream * lo,istream * in,ST::string p = "");
  #else
  statobject(const ST::string & n,const ST::string t,ofstream * lo,istream * in,
				 ST::string p = "");
  #endif

  // COPY CONSTRUCTOR

  statobject(const statobject & so);

  // OVERLOADED ASSIGNMENT OPERATOR

  const statobject & operator=(const statobject & o);

  // DESTRUCTOR

  ~statobject() {}

  // FUNCTION: geterrormessages
  // TASK: returns current errormessages

  const vector<ST::string> & geterrormessages(void) const
	 {
	 return errormessages;
	 }

  // FUNCTION: getname
  // TASK: returns the name of the statobject

  const ST::string & getname(void) const
	 {
	 return name;
	 }

  // FUNCTION: getnameAsBasicString
  // TASK: returns the name of the staobject as an AnsiString

  std::string getnameAsBasisString(void) const
    {
    return name.to_bstr();
    }

  // FUNCTION: gettype
  // TASK: returns the type of the statobject

  const ST::string & gettype(void) const
	 {
	 return type;
	 }

  // FUNCTION: getdescription
  // TASK: returns current describe text

  const vector <ST::string> & getdescription (void) const
    {
    return describetext;
    }

  // VIRTUAL FUNCTION: parse
  // TASK: base function for inherited classes

  virtual int parse(const ST::string & c);

  // FUNCTION: describe
  // TASK: gives a description of the current state of the statobject

  virtual void describe(const optionlist & globaloptions = optionlist());

  // FUNCTION: get_newcommands
  // TASK: returns additional commands to be executed

  const vector<ST::string> & get_newcommands(void) const
    {
    return newcommands;
    }

  friend int operator<(statobject & s1, statobject & s2)
	 {
	 return 0;
	 }

  friend int operator==(statobject & s1, statobject & s2)
	 {
	 return 0;
	 }


  };


int __EXPORT_TYPE findstatobject(
const vector<statobject*> & stats,const ST::string & name,
                   const ST::string & type);



#endif
