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





#include"use.h"


//------------------------------------------------------------------------------
//-------------- CLASS use: implementation of member functions -----------------
//------------------------------------------------------------------------------


const use & use::operator=(const use & u)
  {
  if (this == &u)
	 return *this;
  errormessages = u.errormessages;
  notext = u.notext;
  usingtext = u.usingtext;
  return *this;
  }


//------------------------------------------------------------------------------
//----------- CLASS usePathRead: implementation of member functions ------------
//------------------------------------------------------------------------------


const usePathRead & usePathRead::operator=(const usePathRead & u)
  {
  if (this == &u)
	 return *this;
  errormessages = u.errormessages;
  path = u.path;
  notext = u.notext;
  return *this;
  }


void usePathRead::parse(const ST::string & usetext)
  {

  path = "";
  errormessages.clear();
  notext = true;

  if (usetext.length() > 0)
	 {
	 notext = false;
	 if (usetext.isexistingfile() == 1)
		errormessages.push_back(
		"ERROR: file " + usetext + " could not be opened for reading\n");
	 if (errormessages.empty())
		path = usetext;
	 }

  }


//------------------------------------------------------------------------------
//----------- CLASS usePathWrite: implementation of member functions -----------
//------------------------------------------------------------------------------


const usePathWrite & usePathWrite::operator=(const usePathWrite & u)
  {
  if (this == &u)
	 return *this;
  errormessages = u.errormessages;
  notext = u.notext;
  path = u.path;
  alreadyexisting = u.alreadyexisting;
  return *this;
  }


void usePathWrite::parse(const ST::string & usetext)
  {

  path = "";
  errormessages.clear();
  notext = true;

  if (usetext.length() > 0)
	 {
	 notext = false;
	 int k = usetext.isvalidfile();

	 if (k == 1)
		{
		errormessages.push_back("ERROR: file " + usetext +
		" could not be opened for writing\n");
		alreadyexisting = false;
		}
	 else if (k == 0)
		alreadyexisting = false;
	 else
		alreadyexisting = true;

	 if (errormessages.empty())
		path = usetext;
	 }

  }


//------------------------------------------------------------------------------
//------------ CLASS useDataset: implementation of member functions ------------
//------------------------------------------------------------------------------


useDataset::useDataset(const useDataset & u)
  {
  errormessages = u.errormessages;
  notext = u.notext;
  datasets = u.datasets;
  }


const useDataset & useDataset::operator=(const useDataset & u)
  {
  if (this == &u)
	 return *this;
  errormessages = u.errormessages;
  notext = u.notext;
  datasets = u.datasets;
  return *this;
  }


void useDataset::parse(const ST::string & usetext)
  {

  errormessages.clear();
  notext = true;

  if (usetext.length() > 0)
	 {
	 notext = false;
	 if (! datasets->empty())
		{
		int i = 0;
		bool existing = false;
		while ( (i < datasets->size()) && (existing == false) )
		  {
		  if (usetext == (*((*datasets)[i])).getname())
			 {
			 existing = true;
			 datasetpointer = (*datasets)[i];
			 }
		  i++;
		  }
		if (existing == false)
		  errormessages.push_back(
		  "ERROR: dataset " + usetext + " is not existing\n");

		}
	 else
		errormessages.push_back(
		"ERROR: dataset " + usetext + " is not existing\n");
	 }

  }

