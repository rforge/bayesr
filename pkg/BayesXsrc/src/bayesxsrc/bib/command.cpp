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





#include"command.h"


//------------------------------------------------------------------------------
//------------ CLASS command: implementation of member functions ---------------
//------------------------------------------------------------------------------


command::command(void)
  {
  spec_us    = required;
  spec_model = required;
  spec_exp   = optional;
  spec_opt   = optional;
  spec_by = optional;
  }


command::command(const ST::string & n,model * mr,optionlist * ol,use * ur,
					  specification_allowed spec_m,
					  specification_allowed spec_w,
					  specification_allowed spec_b,
					  specification_allowed spec_e,
					  specification_allowed spec_o,
					  specification_allowed spec_u)
  {
  spec_us = spec_u;
  spec_model = spec_m;
  spec_weight = spec_w;
  spec_by = spec_b;
  spec_exp = spec_e;
  spec_opt = spec_o;
  name = n;
  modelref = mr;
  useref = ur;
  lokaloptions = ol;
  parsingtoken.push_back(name);
  parsingtoken.push_back("weight");
  parsingtoken.push_back("by");
  parsingtoken.push_back("if");
  parsingtoken.push_back(",");
  parsingtoken.push_back("using");
  }


command::command(const command & c)
  {
  name = c.name;
  weighttext = c.weighttext;
  bytext = c.bytext;
  optionstext = c.optionstext;  
  spec_model = c.spec_model;
  spec_weight = c.spec_weight;
  spec_by = c.spec_by;
  spec_exp = c.spec_exp;
  spec_opt = c.spec_opt;
  spec_us = c.spec_us;
  modelref = c.modelref;
  useref = c.useref;
  lokaloptions = c.lokaloptions;
  errormessages = c.errormessages;
  parsingtoken = c.parsingtoken;
  }


const command & command::operator=(const command & c)
  {
  if (this == &c)
	 return *this;
  name = c.name;
  weighttext = c.weighttext;
  bytext = c.bytext;
  optionstext = c.optionstext;
  spec_model = c.spec_model;
  spec_weight = c.spec_weight;
  spec_by = c.spec_by;
  spec_exp = c.spec_exp;
  spec_opt = c.spec_opt;
  spec_us = c.spec_us;
  modelref = c.modelref;
  useref = c.useref;
  lokaloptions = c.lokaloptions;
  errormessages = c.errormessages;
  parsingtoken = c.parsingtoken;
  return *this;
  }


int command::parse(const ST::string & c1)
  {

  ST::string c=c1;

  errormessages.clear();

  lokaloptions->setdefault();

  weighttext = "";
  bytext = "";
  expression = "";
  optionstext = "";

  ST::string modeltext;
//  ST::string optionstext;
  ST::string usetext;

  int modelpos = -1;
  int weightpos = -1;
  int bypos = -1;
  int expressionpos = -1;
  int optionspos = -1;
  int usepos = -1;

  int i=0;
  ST::string cnew;
  while(i<c.length())
	 {
	 if (c[i]==',')
		cnew = cnew + " , ";
	 else
		cnew = cnew+c.substr(i,1);
	 i++;
	 }

  c = cnew;

  bool bracketmiss;
  bool quotmiss;
  vector<ST::string> token = c.strtoken2_quot(parsingtoken,bracketmiss,quotmiss);
  if (bracketmiss)
    {
    errormessages.push_back("ERROR: missing bracket(s)\n");
    return 1;
    }

  if (quotmiss)
    {
    errormessages.push_back("ERROR: missing quotation marks\n");
    return 1;
    }


  if (token[0] == name)
	 {

	 unsigned i = 1;
	 while ( (i < token.size()) && (errormessages.empty()) )
		{
		if (token[i] == "weight")
		  {
		  if (spec_weight == notallowed)
			 errormessages.push_back("ERROR: weight statement not allowed\n");
		  else
			 {
			 if ( (expressionpos != -1) || (optionspos != -1) || (usepos != -1)
					|| (bypos != -1) )
				errormessages.push_back("ERROR: invalid syntax\n");
			 else
				{
				if (weightpos == -1)
				  {
				  if (i+1 < token.size())
					 {
					 if (token[i+1].isinlist(parsingtoken) == -1)
						{
						weightpos = i;
						weighttext = token[i+1];
						}
					 else
						errormessages.push_back
						("ERROR: weight statement incomplete\n");
					 }
				  }
				else
				  errormessages.push_back
				  ("ERROR: too many weight variables specified\n");
				}
			 }
          i = i + 2;
		  }
		else if (token[i] == "by")
		  {
		  if (spec_by == notallowed)
			 errormessages.push_back("ERROR: by statement not allowed\n");
		  else
			 {
			 if ( (expressionpos != -1) || (optionspos != -1) || (usepos != -1) )
				errormessages.push_back("ERROR: invalid syntax\n");
			 else
				{
				if (bypos == -1)
				  {
				  if (i+1 < token.size())
					 {
					 if (token[i+1].isinlist(parsingtoken) == -1)
						{
						bypos = i;
						bytext = token[i+1];
						}
					 else
						errormessages.push_back
						("ERROR: by statement incomplete\n");
					 }
				  }
				else
				  errormessages.push_back
				  ("ERROR: too many by statements specified\n");
				}
			 }
        i = i + 2;
		  }
		else if (token[i] == "if")
		  {
		  if (spec_exp == notallowed)
			 errormessages.push_back("ERROR: boolean expression not allowed\n");
		  else
			 {
			 if ((optionspos != -1) || (usepos != -1))
				errormessages.push_back("ERROR: invalid syntax\n");
			 else
				{
				if (expressionpos == -1)
				  {
				  if (i+1 < token.size())
					 {
					 if (token[i+1].isinlist(parsingtoken) == -1)
						{
						expressionpos = i;
						expression = token[i+1];
						}
					 else
						errormessages.push_back
						("ERROR: boolean expression expected\n");
					 }
				  else
					 errormessages.push_back("ERROR: boolean expression expected\n");
				  }
				else
				  errormessages.push_back
				  ("ERROR: more than one boolean expression specified\n");
				}
			 }
		  i = i + 2;
		  }
		else if (token[i] == ",")
		  {
		  if (spec_opt == notallowed)
			 errormessages.push_back("ERROR: options not allowed\n");
		  else
			 {
			 if (usepos != -1)
				errormessages.push_back("ERROR: invalid syntax\n");
			 else
				{
				if (optionspos == -1)
				  {
				  if (i+1 < token.size())
					 {
					 if (token[i+1].isinlist(parsingtoken) == -1)
						{
						optionspos = i;
						optionstext = token[i+1];
						}
					 else
						errormessages.push_back("ERROR: options expected\n");
					 }
				  else
					 errormessages.push_back("ERROR: options expected\n");
				  }
				else
				  errormessages.push_back
				  ("ERROR: options more than once specified\n");
				}
			 }
		  i = i + 2;
		  }
		else if (token[i] == "using")
		  {
		  if (spec_us == notallowed)
			 errormessages.push_back("ERROR: using not allowed\n");
		  else
			 {
			 if (usepos == -1)
				{
				if (i+1 < token.size())
				  {
				  if (token[i+1].isinlist(parsingtoken) == -1)
					 {
					 usepos = i;
					 usetext = token[i+1];
					 }
				  else
					errormessages.push_back("ERROR: invalid using specification\n");
				  }
				else
				  errormessages.push_back("ERROR: invalid using specification\n");
				}
			 else
				errormessages.push_back("ERROR: using more than once specified\n");
			 }
		  i = i + 2;
		  }
		else
		  {
		  if (spec_model == notallowed)
			 errormessages.push_back("ERROR: model specification not allowed\n");
		  else
			 {
			 if (modelpos == -1)
				{
				if (i == 1)
				  {
				  if (token[i].isinlist(parsingtoken) == -1)
					 {
					 modelpos = 1;
					 modeltext = token[i];
					 }
				  else
					 errormessages.push_back("ERROR: invalid syntax\n");
				  }
				else
				  errormessages.push_back("ERROR: invalid syntax\n");
				}
			 else
				errormessages.push_back("ERROR: more than one model specified\n");
			 }
		  i++;
		  }
		}                              // end while

	 if (errormessages.empty())
		{
		if ( (spec_model == required) && (modeltext.length() == 0) )
		  errormessages.push_back("ERROR: model specification missing\n");
		if ( (spec_weight == required) && (weighttext.length() == 0) )
		  errormessages.push_back("ERROR: weight statement missing\n");
		if ( (spec_by == required) && (bytext.length() == 0) )
		  errormessages.push_back("ERROR: by statement missing\n");
		if ( (spec_exp == required) && (expression.length() == 0) )
		  errormessages.push_back("ERROR: boolean expression required\n");
		if ( (spec_opt == required) && (optionstext.length() == 0) )
		  errormessages.push_back("ERROR: options required\n");
		if ( (spec_us == required) && (usetext.length() == 0) )
		  errormessages.push_back("ERROR: using required\n");
		}

	 if (errormessages.empty())
		{
		if (spec_model != notallowed)       // parse model
		  {
		  modelref->parse(modeltext);
		  errormessages = modelref->geterrormessages();
		  }

		if (spec_opt != notallowed)         // parse options
		  {
		  lokaloptions->parsemultiple(optionstext);
//		  errormessages.insert_back(lokaloptions->geterrormessages());
	      if (! (lokaloptions->geterrormessages()).empty())
            errormessages.insert(errormessages.end(),
            (lokaloptions->geterrormessages()).begin(),
            (lokaloptions->geterrormessages()).end());
		  }

		if (spec_us != notallowed)          // parse using
		  {
		  useref->parse(usetext);
//		  errormessages.insert_back(useref->geterrormessages());

	      if (! (useref->geterrormessages()).empty())
            errormessages.insert(errormessages.end(),
            (useref->geterrormessages()).begin(),
            (useref->geterrormessages()).end());

		  }
		}

	 return 1;
	 }
  else
	 return 0;
  }

