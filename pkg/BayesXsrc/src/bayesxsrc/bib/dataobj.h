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



#if !defined (DATAOBJECT_INCLUDED)

#define DATAOBJECT_INCLUDED

#include"../export_type.h"
#include"statobj.h"
#if defined(JAVA_OUTPUT_WINDOW)
#include "adminparse_pointers.h"
#endif


// #include "describe_dataset.h"


// HINZUFÜGEN EINER METHODE ZUR KLASSE DATAOBJECT:
// 1. runpointer um eins erhöhen
// 2. run (friend)function zur neuen Methode schreiben
// 3. Funktion create ändern


// dataobjects methods:

// METHOD: infile
// SYNTAX: infile [var1 var2 ... varn] using filename (path)
// TASK: fills the dataobject with data stored in file filename
// ADDITIONAL INFORMATION:
// - if there are no variable names specified, a fileheader
//  (first line of the file) containing the variable names is expected

// METHOD: outfile
// SYNTAX: outfile [var1 var2 ... varn] , options using filename (path)
// OPTIONS:
// - header = true|false (default = false) : if header = true, the first
//   line of the file will contain the names of the variables
// TASK: writes data in ASCII file 'filename'
// ADDITIONAL INFORMATION:
// - if no variables are specified, all variables will be writen to filename

// METHOD: generate

// METHOD: replace

// METHOD: drop

// METHOD: rename

// METHOD: set obs

// METHOD: sort


class __EXPORT_TYPE dataobject : public statobject
  {


  private:

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_pointer * adminp_p;
  #endif

  //----------------------- PRIVATE VARIABLES ----------------------------------

  // contains data

  dataset d;

  typedef void (* runpointer)(dataobject &);

  runpointer functions[15];

  modelStandard m;

  expression e;

  usePathRead uread;

  usePathWrite uwrite;

  optionlist emptyoptions;

  // for method 'infile'

  stroption missing;
  intoption maxobs;
  simpleoption nonote;

  optionlist infileoptions;

  friend void infilerun(dataobject & o);

  // for method 'drop'

  friend void droprun(dataobject & o);

  // for method 'rename'

  friend void renamerun(dataobject & o);

  // for method 'generate'

  friend void generaterun(dataobject & o);

  // for method 'replace'

  friend void replacerun(dataobject & o);

  // for method 'set obs'

  friend void setrun(dataobject & o);

  // for method 'outfile'

  optionlist outfileoptions;
  simpleoption header;
  simpleoption replace;
  friend void outfilerun(dataobject & o);

  // for method 'sort'

  optionlist sortoptions;
  simpleoption descending;
  friend void sortrun(dataobject & o);

  // for method 'descriptive'

  optionlist descriptiveoptions;

  friend void descriptiverun(dataobject & o);

  // for method 'tabulate'

  optionlist tabulateoptions;

  friend void tabulaterun(dataobject & o);

  // for method 'pctile'

  optionlist pctileoptions;

  friend void pctilerun(dataobject & o);

  // for method 'marketing'

  optionlist marketingoptions;
  stroption pricedef;
  intoption lag;
  doubleoption alpha;
  friend void marketingrun(dataobject & o);


  //------------------------ PRIVATE FUNCTIONS ---------------------------------

  void create(void);

  void changedescription(void);


  public:


  //------------------------- PUBLIC FUNCTIONS ---------------------------------


  // DEFAULT CONSTRUCTOR

  dataobject(void) : statobject()
	 {
	 type = "dataset";
//     describewindow = datasetwindow;
	 }

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - name = n

  #if defined(JAVA_OUTPUT_WINDOW)
  dataobject(administrator_basic * adb, administrator_pointer * adp,
              const ST::string & n,ofstream * lo,istream * in);
  #else
  dataobject(const ST::string & n,ofstream * lo,istream * in);
  #endif

  // COPY CONSTRUCTOR

  dataobject(const dataobject & o);

  // OVERLOADED ASSIGNMENT OPERATOR

  const dataobject & operator=(const dataobject & o);

  // FUNCTION: allexisting
  // TASK: returns true if all variables specified in varnames are existing
  //       variables that are not existing are stored in vector notex

  bool allexisting(vector<ST::string> & varnames,vector<ST::string> & notex);

  // FUNCTIONS: makematrix
  // TASK: create a designmatrix from the dataset containing the variables with
  //       variablenames 'names'
  // POSSIBLE ERRORS:
  // - one or more variables specified in names are not existing
  // ADDITIONAL INFORMATION:
  // - if names is empty, the datamatrix consists of all variables in the
  // - dataset
  // - if the special name "const" occurs, a column with 1's will be created

  void makematrix(list<ST::string> & names, datamatrix & da,ST::string boole = "")
	 {
	 d.makematrix(names,da,boole);
	 errormessages = d.geterrormessages();
	 }

  void makematrix(vector<ST::string> & names, datamatrix & da, ST::string boole = "")
	 {
	 d.makematrix(names,da,boole);
	 errormessages = d.geterrormessages();
	 }

  void makematrix(ST::string & name, datamatrix & da, ST::string boole = "")
	  {
	  d.makematrix(name,da,boole);
     errormessages = d.geterrormessages();
	  }

  // FUNCTION: parse
  // TASK: parses command c

  int parse(const ST::string & c);

  void describe(const optionlist & globaloptions = optionlist());


  // FUNCTION: getVarnames
  // TASK: returns list of all variable names

  const list<ST::string> & getVarnames() const
	 {
	 return d.getVarnames();
	 }

  // FUNCTION: obs
  // TASK: returns the number ob observations

   const unsigned & obs(void) const
    {
    return d.obs();
    }

  // FUNCTION: set_iterator
  // TASK: sets an iterator to the nr. variable in the dataset
  //       after setting the iterator you can get access to the observations
  //       of the nr. variable with funtion getvalue (see below)

  void set_iterator(const unsigned & nr)
    {
    d.set_iterator(nr);
    }

  // FUNCTION: getvalue
  // TASK: gets the nr th observation
  // IMPORTANT: set first an iterator to the position of the variable you
  //            want observations

  const double & getvalue(const unsigned & nr)
    {

    return d.getvalue(nr);

    }

  dataset getdata(void) const
    {
    return d;
    }


  };

// ------------------------ forward friends decls --------------------------

#if defined (__BUILDING_GNU)
void infilerun(dataobject & o);
void droprun(dataobject & o);
void renamerun(dataobject & o);
void generaterun(dataobject & o);
void replacerun(dataobject & o);
void setrun(dataobject & o);
void outfilerun(dataobject & o);
void sortrun(dataobject & o);
void descriptiverun(dataobject & o);
void tabulaterun(dataobject & o);
void pctilerun(dataobject & o);
void marketingrun(dataobject & o);
#endif

#endif
