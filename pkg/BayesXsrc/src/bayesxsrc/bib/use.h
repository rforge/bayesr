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



#if !defined (USE_INCLUDED)

#define USE_INCLUDED

#include"../export_type.h"

#include<vector>
//#include<errorm.h>
#include"clstring.h"
#include"data.h"


//------------------------------------------------------------------------------
//------------------------------ CLASS use -------------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE use
  {

  protected:


  //------------------------ PROTECTED VARIABLES -------------------------------

  // contains current errormessages

//  errorm::messages errormessages;
  vector<ST::string> errormessages;


  // 'notext' = true, if usetext = ""

  bool notext;

  ST::string usingtext;

  public:


  // ------------------------ PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR

   use(void) {notext = true;}

  // COPY CONSTRUCTOR

   use(const use & u)
	 {
	 errormessages = u.errormessages;
	 notext = u.notext;
	 usingtext = u.usingtext;
	 }

  // OVERLOADED ASSIGNMENT OPERATOR

  const use &  operator=(const use & u);

  // DESTRUCTOR

  ~use(void) {}

  // FUNCTION: geterrormessages
  // TASK: returns current errormessages

  const vector<ST::string> &  geterrormessages(void)
	 {
	 return errormessages;
	 }

  // FUNCTION: getnotext
  // TASK: returns current value of variable 'notext'

  bool  getnotext(void)
	 {
	 return notext;
	 }

  // FUNCTION getusetext
  // TASK: returns current usetext;

  ST::string  getusingtext(void)
	 {
	 return usingtext;
	 }

  // VIRTUAL FUNCTION: parse
  // TASK: base function for inherited classes

  virtual void  parse(const ST::string & usetext)
	 {
	 usingtext=usetext;
	 if (usetext.length() != 0)
		notext = false;
	 }


  };


//------------------------------------------------------------------------------
//--------------------------- CLASS usePathRead --------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE usePathRead : public use
  {


  protected:


  // ------------------------ PROTECTED VARIABLES ------------------------------

  // contains the current path to a file, that can be opened  for reading

  ST::string path;


  public:


  // ------------------------- PUBLIC FUNCTIONS --------------------------------

  // DEFAULT CONSTRUCTOR

   usePathRead(void) {notext = true;}

  // COPY CONSTRUCTOR

   usePathRead(const usePathRead & u)
	 {
	 errormessages = u.errormessages;
	 notext = u.notext;
	 path = u.path;
	 }

  // OVERLOADED ASSIGNMENT OPERATOR

  const usePathRead & __EXPORT_TYPE operator=(const usePathRead & u);

  // FUNCTION: getPath
  // TASK: returns the current path

  const ST::string &  getPath(void)
	 {
	 return path;
	 }

  // FUNCTION: parse
  // TASK: parses usetext, checks, if usetext contains a valid path to a file,
  //       that can be opened for reading

  void  parse(const ST::string & usetext);


  };


//------------------------------------------------------------------------------
//--------------------------- CLASS usePathWrite -------------------------------
//------------------------------------------------------------------------------


class  __EXPORT_TYPE usePathWrite : public use
  {

  protected:


  // ------------------------ PROTECTED VARIABLES ------------------------------

  ST::string path;

  bool alreadyexisting;

  public:


  // -------------------------- PUBLIC FUNCTIONS -------------------------------

  // DEFAULT CONSTRUCTOR

   usePathWrite(void)
	 {
	 alreadyexisting = false;
	 notext = true;
	 }

  // COPY CONSTRUCTOR

   usePathWrite(const usePathWrite & u)
	 {
	 errormessages = u.errormessages;
	 notext = u.notext;
	 path = u.path;
	 alreadyexisting = u.alreadyexisting;
	 }

  // OVERLOADED ASSIGNMENT OPERATOR

  const usePathWrite &  operator=(const usePathWrite & u);

  // FUNCTION: getPath
  // TASK: returns the current path

  const ST::string &  getPath(void)
	 {
	 return path;
	 }

  // FUNCTION: isexisting
  // TASK: returns true if the current path (file) is already existing,
  //       else false

  bool  isexisting(void)
	 {
	 return alreadyexisting;
	 }

  // FUNCTION: parse
  // TASK: parses usetext, checks, if usetext contains a valid path to a file,
  //       that can be opened for reading

  void  parse(const ST::string & usetext);


  };


//------------------------------------------------------------------------------
//-------------------------- CLASS useDataset ----------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE useDataset : public use
  {


  protected:

  // ---------------------- PROTECTED VARIABLES --------------------------------

  vector<dataset*> * datasets;

  dataset * datasetpointer;

  public:

  // ----------------------- PUBLIC FUNCTIONS ----------------------------------

  // DEFAULT CONSTRUCTOR

   useDataset(void) {notext = false;}

  // CONSTRUCTOR

   useDataset(vector<dataset*> * d)
	 {
	 datasets = d;
	 }

  // COPY CONSTRUCTOR

   useDataset(const useDataset & d);

  // OVERLOADED ASSIGNMENT OPERATOR

  const useDataset &  operator=(const useDataset & d);


  void  parse(const ST::string & usetext);


  dataset *  getDatasetpointer(void)
	 {
	 return datasetpointer;
	 }


  };


#endif





