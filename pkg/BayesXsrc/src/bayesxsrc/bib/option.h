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



#if !defined (OPTION_INCLUDED)

#define OPTION_INCLUDED

#include"../export_type.h"
#include<fstream>
#include<vector>
#include"clstring.h"
//#include<errorm.h>
#include"../values.h"
#include <climits>

// POSSIBLE OPTIONS:

// simple option (class simpleoption):
// SYNTAX: optionname
// ADDITIONAL INFORMATION:
// - if the optionname occurs in an optionlist, the option is activated

// stroption (class stroption):
// SYNTAX: optionname = newvalue
// ADDITIONAL INFORMATION: type of newvalue should be double

// intoption (class intoption):
// SYNTAX: optionname = newvalue
// ADDITIONAL INFORMATION: type of newvalue should be integer

// doubleoption (class doubleoption):
// SYNTAX: optionname = newvalue
// ADDITIONAL INFORMATION: type of newvalue should be double


//------------------------------------------------------------------------------
//---------------------------- CLASS : option ----------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE option

  {

  protected:


  // PROTECTED VARIABLES:

  // name of the option

  ST::string optionname;

  // true if current value has been changed

  bool valuechanged;

  // contains current errormessages

  // errorm::messages errormessages;
  vector<ST::string> errormessages;

  unsigned maxtoken;


  public:


  // PUBLIC FUNCTIONS:

  // DEFAULT CONSTRUCTOR

  option(void);

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - optionname = n

  option(const ST::string & n);

  // COPY CONSTRUCTOR

  option(const option & b);

  // DESTRUCTOR

  ~option(void) {}

  // OVERLOADED ASSIGNMENT OPERATOR

  const option & operator=(const option & b);

  // FUNCTION: geterrormessages
  // TASK: returns current errormessages

//  const errorm::messages & geterrormessages(void)
  const vector<ST::string> & geterrormessages(void)
	 {
	 return errormessages;
	 }

  // FUNCTION: getname
  // TASK: returns the name of the option

  const ST::string & getname(void)
	 {
	 return optionname;
	 }

   const unsigned getmaxtoken(void)
     {
     return maxtoken;
     }

  // VIRTUAL FUNCTION: parse
  // TASK: base function for inherited classes
  //       parses string c

  virtual int parse(const ST::string & c)
	 {
	 return 0;
	 }

  const bool & changed(void)
    {
    return valuechanged;
    }

  // VIRTUAL FUNCTION: setdefault
  // TASK: base function for inherited classes
  //       sets current value of the option to its default value

  virtual void setdefault(void) {}

  // VIRTUAL FUNCTION: getValueAsString
  // TASK: basefunction for inherited classes
  //       returns the current value of the option as a string

  virtual ST::string getValueAsString(void)
	 {
	 return "";
	 }


  };


//------------------------------------------------------------------------------
//-------------------------- CLASS simpleoption --------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE simpleoption : public option
  {

  protected:


  // ------------------------ PROTECTED VARIABLES ------------------------------

  // default value of the option

  bool defaultvalue;

  // value of the option, true = option activated, false = option not activated

  bool value;


  public:


  // ------------------------- PUBLIC FUNCTIONS --------------------------------

  // DEFAULT CONSTRUCTOR

  simpleoption(void) : option() {}

  // CONSTRUCTOR

  simpleoption(const ST::string & n,bool v);

  // COPY CONSTRUCTOR

  simpleoption(const simpleoption & o);

  // OVERLOADED ASSIGNMENT OPERATOR

  const simpleoption & operator=(const simpleoption & o);

  // FUNCTION: parse
  // TASK: parses the string c

  int parse(const ST::string & c);

  // FUNCTION: setdefault
  // TASK: sets current value 'value' to default value 'defaultvalue'

  void setdefault(void)
	 {
     valuechanged = false;
	 value = defaultvalue;
	 }

  // FUNCTION: getvalue
  // TASK: returns current value of the option

  bool getvalue(void)
	 {
	 return value;
	 }

  virtual ST::string getValueAsString(void)
	 {
	 if (value==true)
		return "true";
	 else
		return "false";
	 }


  };


//------------------------------------------------------------------------------
//------------------------- CLASS: fileoption ----------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE fileoption : public option
  {

  protected:

  // ------------------------ PROTECTED VARIABLES ------------------------------

  // default value of the option (= default filename)

  ST::string defaultvalue;

  // value of the option(= filename)

  ST::string value;

  // filtetype = true = specification of file type allowed
  // filtetype = false = specification of file type not allowed

  bool filetype;


  public:

  // -------------------------- PUBLIC FUNCTIONS -------------------------------


  // DEFAULT CONSTRUCTOR

  fileoption(void) : option() {}

  // CONSTRUCTOR

  fileoption(const ST::string & n,const ST::string & v,bool t=true);

  // COPY CONSTRUCTOR

  fileoption(const fileoption & o);

  // OVERLOADED ASSIGNMENT OPERATOR

  const fileoption & operator=(const fileoption & o);

  // FUNCTION: parse
  // TASK: parse the string c

  int parse (const ST::string & c);

  // FUNCTION: setdefault
  // TASK: sets current value 'value' to default value 'defaultvalue'

  void setdefault(void)
	 {
     valuechanged = false;
	 value = defaultvalue;
	 }

  // FUNCTION: getvalue
  // TASK: returns current value of the option

  ST::string getvalue(void)
	 {
	 return value;
	 }

  // VIRTUAL FUNCTION: getValueAsString
  // TASK: returns the current value of the option as a string

  virtual ST::string getValueAsString(void)
	 {
	 return value;
	 }


  };


//------------------------------------------------------------------------------
//------------------------ END class fileoption --------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//------------------------- CLASS: fileoption2 ----------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE fileoption2 : public fileoption
  {

  protected:


  public:

  // -------------------------- PUBLIC FUNCTIONS -------------------------------


  // DEFAULT CONSTRUCTOR

  fileoption2(void) : fileoption() {}

  // CONSTRUCTOR

  fileoption2(const ST::string & n,const ST::string & v,bool t=true);

  // COPY CONSTRUCTOR

  fileoption2(const fileoption2 & o);

  // OVERLOADED ASSIGNMENT OPERATOR

  const fileoption2 & operator=(const fileoption2 & o);

  // FUNCTION: parse
  // TASK: parse the string c

  int parse (const ST::string & c);

  };


//------------------------------------------------------------------------------
//------------------------ END class fileoption2 --------------------------------
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
//------------------------- CLASS: stroption -----------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE stroption : public option

  {

  protected:


  // PROTECTED VARIABLES:


  // allallowed is true, if all strings are allowed

  bool allallowed;

  // contains admissible values

  vector<ST::string> admissible;

  // contains default value

  ST::string defaultvalue;

  // contains current value

  ST::string value;


  public:


  // PUBLIC FUNCTIONS:

  // DEFAULT CONSTRUCTOR

  stroption(void);

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - optionname = n,
  // - admissible = adm,
  // - allallowed = false,
  // - defaultvalue = value = v

  stroption(const ST::string & n,const vector<ST::string> & adm,const ST::string & v);

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - optionname = n,
  // - allallowed = true,
  // - defaultvalue = value = v

  stroption(const ST::string & n,const ST::string & v);

  // CONSTRUCTOR
  // ADDITION INFORMATION
  // - optionname = n
  // - allallowed = true

  stroption(const ST::string & n);

  // COPY CONSTRUCTOR

  stroption(const stroption & b);

  // OVERLOADED ASSIGNMENT OPERATOR

  const stroption & operator=(const stroption & b);

  // FUNCTION: parse
  // TASK: parses the option
  //       returns 1, if option is recognized
  //       (i.e. first token of c equals optionname
  //       returns 0, if option is not recognized
  // POSSIBLE ERRORS:
  // - second token is not '='
  // - command contains more than three token (invalid option specification)
  // - specified new value is unknown (i.e. is not in admissible)
  // - new value is missing (i.e. less than three token)
  // ADDITIONAL INFORMATION:
  // - errormessages will be deleted before (re-)parsing,
  //   if 'clearerrors = true'

  int parse(const ST::string & c);

  // FUNCTION: setdefault
  // TASK: sets current value 'value' to default value 'defaultvalue'

  void setdefault(void)
	 {
     valuechanged = false;
	 value = defaultvalue;
	 }

  // FUNCTION: getvalue
  // TASK: returns current value of the option

  const ST::string & getvalue()
	 {
	 return value;
	 }

  // VIRTUAL FUNCTION: getValueAsString
  // TASK: returns the current value of the option as a string

  virtual ST::string getValueAsString(void)
	 {
	 return value;
	 }


  };


//------------------------------------------------------------------------------
//------------------------- END class stroption --------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//------------------------- CLASS: intoption -----------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE intoption : public option

  {

  protected:


  // PROTECTED VARIABLES:

  int lowerbound;

  int upperbound;

  // default value of the option

  int defaultvalue;

  // current value of the option, value must be in [lowerbound,upperbound]

  int value;


  public:


  // PUBLIC FUNCTIONS:

  // DEFAULT CONSTRUCTOR

  intoption(void);

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - optionname = n,
  // - defaultvalue = value = v
  // - lowerbound = lb
  // - upperbound = ub

  intoption(const ST::string & n,const int v,
				const int lb = INT_MIN,const int ub = INT_MAX);

  // COPY CONSTRUCTOR

  intoption(const intoption & b);

  // OVERLOADED ASSIGNMENT OPERATOR

  const intoption & operator=(const intoption & b);

  // FUNCTION: parse
  // TASK: parses the option
  //       returns 1, if option is recognized
  //       (i.e. first token of c equals optionname
  //       returns 0, if option is not recognized
  // POSSIBLE ERRORS:
  // - second token is not '='
  // - command contains more than three token (invalid option specification)
  // - specified new value is not in [lowerbound, upperbound]
  // - new value is missing (i.e. less than three token)
  // ADDITIONAL INFORMATION:
  // - errormessages will be deleted before (re-)parsing
  //   if clearerrors = true

  int parse(const ST::string & c);

  // FUNCTION: setdefault
  // TASK: changes value to defaultvalue

  void setdefault(void)
	 {
     valuechanged = false;
	 value = defaultvalue;
	 }

  void setvalue(const int & v)
     {
     assert(v >= lowerbound);
     assert(v <= upperbound);
     value = v;
     }

  // FUNCTION: getvalue
  // TASK: returns current value of the option

  int getvalue()
	 {
	 return value;
	 }

  // VIRTUAL FUNCTION: getValueAsString
  // TASK: returns the current value of the option as a string

  virtual ST::string getValueAsString(void)
	 {
	 return ST::inttostring(value);
	 }


  };


//------------------------------------------------------------------------------
//------------------------- END class intoption --------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//------------------------- CLASS: doubleoption --------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE doubleoption : public option

  {

  protected:


  // PROTECTED VARIABLES:

  double lowerbound;

  double upperbound;

  // contains default value of the option

  double defaultvalue;

  // value of the option, value must be in [lowerbound,upperbound]

  double value;


  public:


  // PUBLIC FUNCTIONS:

  // DEFAULT CONSTRUCTOR

  doubleoption(void);

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - optionname = n,
  // - defaultvalue = value = v
  // - lowerbound = lb
  // - upperbound = ub

  doubleoption(const ST::string & n,const double v,const double lb = MINDOUBLE,
					const double ub = MAXDOUBLE);

  // COPY CONSTRUCTOR

  doubleoption(const doubleoption & b);

  // OVERLOADED ASSIGNMENT OPERATOR

  const doubleoption & operator=(const doubleoption & b);

  // FUNCTION: parse
  // TASK: parses the option
  //       returns 1, if option is recognized
  //       (i.e. first token of c equals optionname
  //       returns 0, if option is not recognized
  // POSSIBLE ERRORS:
  // - second token is not '='
  // - command contains more than three token (invalid option specification)
  // - specified new value is not in [lowerbound, upperbound]
  // - new value is missing (i.e. less than three token)
  // ADDITIONAL INFORMATION:
  // - errormessages will be deleted before (re-)parsing
  //   if clearerrors = true

  int parse(const ST::string & c);

  // FUNCTION: setdefault
  // TASK: changes value to defaultvalue

  void setdefault(void);

  void setvalue(const double & v)
     {
     assert(v >= lowerbound);
     assert(v <= upperbound);
     value = v;
     }

  // FUNCTION: getvalue
  // TASK: returns current value of the option

  double getvalue()
	 {
	 return value;
	 }

  // VIRTUAL FUNCTION: getValueAsString
  // TASK: returns the current value of the option as a string

  virtual ST::string getValueAsString(void)
	 {
	 return ST::doubletostring(value);
	 }


  };


//------------------------------------------------------------------------------
//----------------------- END class doubleoption -------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//--------------------------- CLASS: optionlist --------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE optionlist : public vector<option *>
  {

  protected:


  // PROTECTED VARIABLES:

  // contains current errormessages

//  errorm::messages errormessages;
  vector<ST::string> errormessages;


  // PROTECTED FUNCTIONS:

  // FUNCTION: find
  // TASK: returns the position of the option called optionname
  //       end() is returned if option is not existing

  optionlist::iterator find (const ST::string & optionname);


  public:


  // PUBLIC FUNCTIONS:

  // DEFAULT CONSTRUCTOR

  optionlist(void) : vector<option *>()
	 {
	 }

  // COPY CONSTRUCTOR

  optionlist(const optionlist & ol) : vector<option *>(vector<option*>(ol)) {}

  // OVERLOADED ASSIGNMENT OPERATOR

  const optionlist & operator=(const optionlist & ol);

  // FUNCTION: geterrormessages
  // TASK: returns current errormessages

  const vector<ST::string> & geterrormessages(void)
	 {
	 return errormessages;
	 }

  // FUNCTION: parse
  // TASK: parses the string c
  //       if an option is recognized the current value of the option will
  //       be changed
  //       returns 1, if an option is recognized
  //               0, otherwise
  // ADDITIONAL INFORMATION:
  // - errormessages will be deleted before (re-)parsing,
  //   if clearerrors = true

  int parse(const ST::string & c,bool clearerrors);

  // FUNCTION: parsemultiple
  // TASK: parses all options stored in c
  // POSSIBLE ERRORS:
  // - a specified option is unknown
  // - see also the different versions of function parse above

  void parsemultiple(const ST::string & c);

  // FUNCTION: setdafault
  // TASK: sets current value of all options in the optionlist to their
  //       default values

  void setdefault(void);


  };


#endif
