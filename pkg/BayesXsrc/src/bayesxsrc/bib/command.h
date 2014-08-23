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



#if !defined (COMMAND_INCLUDED)

#define COMMAND_INCLUDED

#include"../export_type.h"
#include<vector>
#include"clstring.h"
#include"model.h"
#include"option.h"
#include"use.h"


enum specification_allowed {required,optional,notallowed};


// SYNTAX OF COMMANDS:
// name [model] [weight varname] [by varname] [if expression]
//      [, options] [using usingtext]


class __EXPORT_TYPE command
  {

  protected:


  //------------------------- PROTECTED VARIABLES ------------------------------

  ST::string weighttext;

  ST::string bytext;

  ST::string expression;

  ST::string optionstext;

  vector<ST::string> parsingtoken;

  // name of the command

  ST::string name;

  // stores if a model is required, optional or not allowed

  specification_allowed spec_model;

  // stores, if specification of a weight variable is required, optional or not
  // allowed

  specification_allowed spec_weight;

  // stores if a by statement is required, optional or not allowed

  specification_allowed spec_by;

  // stores if a if statement is required, optional or not allowed

  specification_allowed spec_exp;

  // stores if options specification is required, optional or not allowed

  specification_allowed spec_opt;

  // stores if a using statement is required, optional or not allowed

  specification_allowed spec_us;

  // reference to a specific model type (see file model.h for different types)

  model* modelref;

  // reference to a list of options (options specific to this command)

  optionlist* lokaloptions;

  // reference to a using type (see file use.h for different types)

  use * useref;

  // stores current errormessages
  // errormessages will be deleted, if command is (re-)parsed

  vector<ST::string> errormessages;
//  errorm::messages errormessages;


  public:


  //----------------------- PUBLIC FUNCTIONS -----------------------------------

  // DEFAULT CONSTRUCTOR

  command(void);

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - defaults:
  //   using required
  //   model specification required
  //   weight variable optional
  //   boolean expression optional
  //   options optional

  command(const ST::string & n,model * mr,optionlist * ol,use * ur,
			 specification_allowed spec_m  = required,
			 specification_allowed spec_w = optional,
			 specification_allowed spec_b = optional,
			 specification_allowed spec_e = optional,
			 specification_allowed spec_o = optional,
			 specification_allowed spec_u = required);

  // COPY CONSTRUCTOR

  command(const command & c);

  // OVERLOADED ASSIGNMENT OPERATOR

  const command & operator=(const command & c);

  // FUNCTION: parse
  // TASK: parses the command c
  //       returns 1, if command (with commandname 'name') is recognized
  //               0, if not recognized
  // ADDITIONAL INFORMATION:
  // - current errormessages will be deleted
  // - option values are set to their default values before (re-)parsing

  int parse(const ST::string & c1);

  // FUNCTION: get_weight_variable
  // TASK: returns the weight variable name (if specified)
  //       if no by variable has been specified, an empty string will be
  //       returned

  ST::string get_weight_variable(void)
	 {
	 return weighttext;
	 }

  // FUNCTION: get_by_variable
  // TASK: returns the by variable name (if specified)
  //       if no by variable has been specified, an empty string will be
  //       returned

  ST::string get_by_variable(void)
	 {
	 return bytext;
	 }

	// FUNCTION: getexpression
   // TASK: returns expression text

  ST::string getexpression(void)
	 {
	 return expression;
	 }

   // FUNCTION: getoptionstext
   // TASK: returns the optionstext

   ST::string getoptionstext(void)
     {
     return optionstext;
     }

  // FUNCTION: geterrormessages
  // TASK: returns errormessages

  const vector<ST::string> & geterrormessages(void)
	 {
	 return errormessages;
	 }

  ~command() {}


  };


#endif
