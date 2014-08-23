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



#if !defined (DATA_INCLUDED)

#define DATA_INCLUDED

#include"../export_type.h"
#include<fstream>
#include"clstring.h"
#include"realobs.h"
#include"realvar.h"
#include"statmat.h"
#include"vectorn.h"
#include<list>

#if defined(JAVA_OUTPUT_WINDOW)
#include"adminparse_basic.h"
#endif


//typedef statmatrix<double> __EXPORT_TYPE datamatrix;

using realob::NA;
using realob::realobs;
using realob::realvar;
using realob::binomial;
using realob::bernoulli;
using realob::normal;
using realob::exponential;
using realob::cumul;
using realob::lagrealvar;
using realob::gamma;

using std::ifstream;
using std::ofstream;
using std::ostream;

//------------------------------------------------------------------------------
//---------------------------- CLASS: data -------------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE data

  {

  friend class filter;
  friend class dataset;


  protected:


  //----------------------- PROTECTED VARIABLES --------------------------------

  // true = no variables (and no observations) present

  bool empty;

  // contains variablenames

  list<ST::string> varnames;

  // contains data vectors

  list< realvar > variables;

  // contains the current order of the observations

  vector<int> index;


  //--------------------- PROTECTED FUNCTIONS ----------------------------------

  // FUNCTION: indexcreate

  void indexcreate(void);

  // FUNCTION: makeitlist
  // TASK: creates a list of iterators (stored in 'itl'), that point to
  //       variables specified in 'names'
  //       returns: 1 if one or more variables in 'name' could not be found
  //                0 if no error occured
  // ADDITIONAL INFORMATION: the special variablename "const" will produce no
  //                         error

  int makeitlist(list<ST::string> & names, vector<list<realvar>::iterator> & itl);

  // FUNCTION: findvar
  // TASK: find variable with variablename 'name'
  //       returns: 1 if variable name is not found
  //                0 if variable name is found (no error occured)
  //                in that case stit and varit contain the position of the
  //                variable in varnames and variables

  int findvar(const ST::string & name,list<ST::string>::iterator & stit,
			 list<realvar>::iterator & varit);

  // FUNCTION: findvar
  // TASK: find varaible with varname 'name'
  //       returns: 1 if variable name is not found
  //                0 if variable name is found (no error occured)

  int findvar(const ST::string & name)// const
	 {
	 list<ST::string>::iterator stit;
	 list<realvar>::iterator varit;
	 return findvar(name,stit,varit);
	 }

  // DEFAULT CONSTRUCTOR

  data(void)
	 {
	 empty = true;
	 }

  // CONSTRUCTOR

  data(list <ST::string> & vn,list<realvar> & va);

  // COPY CONSTRUCTOR

  data(const data & d);

 // OVERLOADED ASSIGNMENT OPERATOR

  const data & operator=(const data & d);

  // DESTRUCTOR

  ~data()
  {}

  // OVERLOADED << OPERATOR

  friend ostream & operator<<(ostream & out,data & d);

  // FUNCTION: obs
  // TASK: returns number of observations in the dataset

  unsigned obs() const;

  // FUNCTION: clear
  // drops all variables and observations, empty is set to true

  void clear();


  };

//------------------------------------------------------------------------------
//--------------------------- END CLASS: data ----------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//---------------------------- CLASS: filter -----------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE filter : public vector<bool>
  {


  // number of filtered observations (number of false's in filter)

  int sum;

  public:

  // DEFAULT CONSTRUCTOR

  filter(void) : vector<bool>()
	 {
	 sum =0;
	 }

  // CONSTRUCTOR
  // TASK: creates a filter with dimension nr,
  //       filter values will be set to true

  filter(int nr) : vector<bool>(nr,true)
	 {
	 sum = 0;
	 }

  // COPY CONSTRUCTOR

  filter(const filter & f) : vector<bool>(vector<bool>(f))
	 {
	 sum = f.sum;
	 }

  // COPY CONSTRUCTOR

  filter(vector<bool> & v);

  // OVERLOADED ASSIGNMENT OPERATOR

  const filter & operator=(const filter & f)
	 {
	 vector<bool>::operator=(vector<bool>(f));
	 sum = f.sum;
	 return *this;
	 }

  // OVERLOADED ASSIGNMENT OPERATOR

  const filter & operator=(vector<bool> & v);

  // FUNCTION: nrunfiltered
  // TASK: returns number of unfiltered observations (number of false's)

  int nrunfiltered()
	 {
	 return (size()-sum);
	 }

  // FUNCTION: filter NA
  // TASK: returns a filter
  //       a value is true, if one or more variables specified in names
  //       have missing values
  // ADDITIONAL INFORMATION: filter will not be changed if messages occur
  //                         i.e. 'names' contain unknown variablenames

  void filterNA(::data & d, list<ST::string> & names);

  // OVERLOADED + OPERATOR
  // TASK: a value of the resulting filter is true if
  //       either the corresponding element of the first filter or of the
  //       second filter is true, false otherwise

  filter operator+(filter & f);


  filter operator+(realvar & v);


  };


//------------------------------------------------------------------------------
//-------------------------- END: class filter ---------------------------------
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
//----------------------------- CLASS dataset ----------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE dataset
  {

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adminb_p;
  #endif

  protected:


  //-------------------------- PROTECTED VARIABLES -----------------------------

  // name of the dataset

  ST::string name;

  // Datarepresentation (see class data above)

  data datarep;

  list<realvar>::iterator varit;

  // number of observations

  unsigned nrobs;

  filter f;

  // contains current errormessages

  vector<ST::string> errormessages;
//  errorm::messages errormessages;


  //------------------------- PROTECTED FUNCTIONS ------------------------------

  void filldata(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adminb_p,
  #endif
  istream & in, ST::string & m,const unsigned & maxobs);

  // FUNCTION: checkvarnames
  // TASK: checks, if the variable names stored in 'varnames' are valid

  void checkvarnames(void);


  public:


  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR
  // TASK: creates an empty dataset

  dataset(void)
	 {
	 datarep = data();
	 name = "noname";
	 nrobs = 0;
	 }

  // CONSTRUCTOR

#if defined(JAVA_OUTPUT_WINDOW)
  dataset(const ST::string & n,administrator_basic * adb)
#else
  dataset(const ST::string & n)
#endif
	 {
     #if defined(JAVA_OUTPUT_WINDOW)
     adminb_p = adb;
     #endif
	 datarep = data();
	 name = n;
	 nrobs = 0;
	 }

  // COPYCONSTRUCTOR

  dataset(const dataset & d);

  // OVERLOADED ASSIGNMENT OPERATOR

  const dataset & operator=(const dataset & d);

  // DESTRUCTOR

  ~dataset() {}

  // FUNCTION: read
  // TASK: fill the dataset with data
  //       if there are no variable names specified, the function expects, that
  //       the first line of the input file contains the variable names
  // POSSIBLE ERRORS:
  // - Header contains invalid variable specification
  //   (some variable names are not valid)
  // - some observations are no double numbers
  // - some variables have more observations than others
  // ADDITIONAL INFORMATION:
  // - errormessages will be cleared before reading new data
  // - if an error occurs, when trying to read new data,
  //   the dataset will be cleared (dataset is empty)
  // - maximal linesize: 8192 byte

  void read(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adminb_p,
  #endif
  ifstream & in,ST::string & missing,const unsigned & maxobs,
            const list<ST::string> & names = list<ST::string>());

  // FUNCTION: write
  // TASK: writes variables 'names' in ASCII-format to an external file
  //       if 'v' = 1
  // ADDITIONAL INFORMATION:
  // - if 'v' is missing all observations will be written to the external file
  // - if 'header' = true (default = false) the first line of the external file
  //   consists of a header, containing the variable names
  // - if names is empty, all variables in the dataset will be wirtten to the
  //   file

  unsigned write (
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adminb_p,
  #endif
  ostream & out,list<ST::string> & names, // = list<ST::string>(),
				  const bool header = false, const realvar & v0 = realvar());

  // FUNCTION: eval_exp
  // TASK: evaluates an expression and returns the resulting variable
  // VALID EXPRESSION:
  //  - Mathematical functions: sin,cos,log,log10,exp,sqrt,abs,floor
  //  - Statistical functions:
  //     - uniform() : in (0,1) uniformly distributed random numbers
  //     - normal()  : standard normal distributed random numbers
  //     - cumul(v)  : empirical distribution function of v
  //  - unary operators: + -
  //  - binary operators:
  //     - arithmetic: * / + - ^
  //     - comparison: = < > <= >=
  //     - logical:    & |
  //  - Special Symbols:
  //      - _n  : returns the vector (1,2,3,4,5,...,nr of obs)
  //      - _pi : returns a vector containing pi
  // ADDITIONAL INFORMATION:
  // - errmormessages will be deleted, if clearerrors = true

  realvar eval_exp(ST::string e, bool clearerrors = true);

  // FUNCTION: sort
  // TASK: sorts the dataset with respect to the variable at iterator position
  //       'varit' between observations 'left' and 'right'

  void sort(list<realvar>::iterator varit,int left,int right);

  // FUNCTION: sort
  // TASK: sorts the dataset with respect to variables 'names' between
  //       observations 'left' and 'right'
  // POSSIBLE ERRORS:
  // - at least one variable specified in 'names' is not existing
  // ADDITIONAL INFORMATION:
  // - errormessages will be deleted

  void sort (list<ST::string> names,int left,int right);

  // OVERLOADED << OPERATOR
  // TASK: Print the dataset
  // ADDITIONAL INFORMATION:
  // - in case of an empty dataset the,
  //   message "dataset is empty" will be displayed

  friend ostream & operator<<(ostream & out,dataset & d)
	 {
	 return out << d.datarep;
	 }

  // FUNCTION: geterrormessages
  // TASK: returns possible errormessages produced by the latest called
  //       dataset method

  const vector<ST::string> & geterrormessages() const
	 {
	 return errormessages;
	 }

  // FUNCTION: getVarnames
  // TASK: returns list of variable names

  const list<ST::string> & getVarnames() const
	 {
	 return datarep.varnames;
	 }

  // FUNCTION: replaceName
  // TASK: changes the name of the dataset to 'newname'

  void replaceName(ST::string & newname)
	 {
	 name = newname;
	 }

  // FUNCTION: rename
  // TASK: changes the name of variable 'oldname' to 'newname'

  void rename(const ST::string & oldname, const ST::string & newname);

  // FUNCTION: addvariable
  // TASK: adds variable 'v' with variable name 'name' to the dataset
  // POSSIBLE ERRORS:
  // - 'name' is not a valid variable name
  // ADDITIONAL INFORMATION:
  // - errormessages will be cleared before adding a new variable

  void addvariable(const ST::string & name, const realvar & v);

  // FUNCTION: addvariable
  // TASK: adds variable 'name' = e to the dataset (e is an expression)
  // POSSIBLE ERRORS:
  // - 'name' is not a valid variable name
  // ADDITIONAL INFORMATION:
  // - errormessages will be cleared before adding a new variable

  void addvariable(const ST::string & name, const ST::string & e);

  // FUNCTION: replace
  // TASK: replaces the observations of variable 'name' with the observations
  //       in 'v', returns the number of changed observations
  // POSSIBLE ERRORS:
  // - 'name' is an unknown variable
  // ADDITIONAL INFOEMATION:
  // - errormessages will be cleared before replacing the variable

  unsigned replace(const ST::string & name, const realvar & v);

  // FUNCTION: replace
  // TASK: changes variable 'name' to 'name' = e (where e is an expression)
  //       if a boolean expression 'boole' is specified, the variable
  //       will be changed only for observations with boole > 0, i.e.
  //       boole = true.
  //       returns the number of changed observations
  // POSSIBLE ERRORS:
  // - 'name' is an unknown variable
  // - expression syntax error
  // ADDITIONAL INFORMATION:
  // - errormessages will be cleared before replacing the variable

  unsigned replace(const ST::string & name, const ST::string & e,const ST::string & boole = "");

  // FUNCTION: dropvariable
  // TASK: drops variable with variablename 'name' from the dataset
  // POSSIBLE ERRORS:
  // - variable 'name' can not be found in the dataset
  // ADDITIONAL INFORMATION:
  // - errormessages will be cleared before dropping a variable

  void dropvariable(const ST::string & name);

  // FUNCTION: dropvariables
  // TASK: drops all variables specified in names
  // POSSIBLE ERRORS:
  // - one of the variables specified in 'names' can not be found in the
  //   dataset
  // ADDITIONAL INFORMATION:
  // - errormessages will be cleared before dropping variables
  // - if 'names' is empty all variables in the dataset will be deleted
  // - if the dataset is empty and 'names' is empty, nothing will happen

  void dropvariables(list<ST::string> & names);


  // FUNCTION: dropobservations
  // TASK: drops all observations with e = 1 (i.e. e = true)
  //       if e is empty, the procedure will do nothing
  //       'nrelim' is the number of eliminated observations
  // POSSIBLE ERRORS:
  // - expression syntax error
  // ADDITIONAL INFORMATION:
  // - errormessages will be cleared before dropping observations.

  unsigned dropobservations(const ST::string & e);


  // FUNCTION: makematrix
  // TASK: create a designmatrix from the dataset containing the variables with
  //       variablenames 'names'
  // POSSIBLE ERRORS:
  // - one or more variables specified in names are not existing
  // ADDITIONAL INFORMATION:
  // - if names is empty, the datamatrix consists of all variables in the
  // - dataset
  // - errormessages will be cleared before creating the datamatrix
  // - if the special name "const" occurs, a column with 1's will be created

  void makematrix(list<ST::string> & names, datamatrix & d,ST::string boole = "");

  void makematrix(vector<ST::string> & names, datamatrix & d, ST::string boole = "");

  void makematrix(ST::string & name, datamatrix & d, ST::string boole = "");

  void marketing(vector<ST::string> & names, ST::string & defs, int & lak, double & alph);

  const ST::string & getname(void)
	 {
	 return name;
	 }


  // ------------------------- ACCESS TO THE DATA ------------------------------

  // FUNCTION: obs
  // TASK: returns number of observations in the dataset

  const unsigned & obs(void) const
	 {
	 return nrobs;
	 }

  // FUNCTION: setobs
  // TASK: changes number of observations in the dataset
  //       newobs >= nrobs (current number of observations) is recommended
  // POSSIBLE ERRORS:
  //
  // ADDITIONAL INFORMATION:
  // - errormessages will be deleted

  void setobs(unsigned newobs);


  // FUNCTION: findvar
  // TASK: find varaible with varname 'name'
  //       returns: 1 if variable name is not found
  //                0 if variable name is found (no error occured)

  int findvar(ST::string & name)
	 {
	 return datarep.findvar(name);
	 }


  // FUNCTION: varnr
  // TASK: returns number of variables in the dataset

  unsigned varnr(void)
	 {
	 return datarep.varnames.size();
	 }

  // FUNCTION: set_iterator
  // TASK: sets an iterator to the nr. variable in the dataset
  //       after setting the iterator you can get access to the observations
  //       of the nr. variable with funtion getvalue (see below)

  void set_iterator(const unsigned & nr)
    {
    varit = datarep.variables.begin();
    unsigned i;
    for(i=1;i<nr;i++)
      ++varit;
    }

  // FUNCTION: getvalue
  // TASK: gets the nr th observation
  // IMPORTANT: set first an iterator to the position of the variable you
  //            want observations

  const double & getvalue(const unsigned & nr)
    {

    return ((*varit)[datarep.index[nr]]).getvalue();

    }


  const double & getvalue(const unsigned & nr,list<realvar>::iterator & vit)
    {
    return ((*vit)[datarep.index[nr]]).getvalue();
    }



  // FUNCTION: reverseorder
  // TASK: reverses the order of the data

  void reverseorder(void);

  };


#endif

