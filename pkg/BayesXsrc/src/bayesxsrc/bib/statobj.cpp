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

#pragma package(smart_init)

#include "StatResults.h"
#include "statwinframe.h"

#endif

#include"statobj.h"

using std::flush;
using std::cout;

//------------------------------------------------------------------------------
//----------- CLASS statobject: implementation of member functions -------------
//------------------------------------------------------------------------------


int statobject::parsecom(const ST::string & c, vector<command> & methods,
								 optionlist & globaloptions)
  {

  errormessages.clear();

  if (c=="describe")
	 {
	 describe(globaloptions);
	 return -1;
	 }

  if (globaloptions.parse(c,true) == 1)
	 {
	 errormessages = globaloptions.geterrormessages();
	 if (errormessages.empty())
		return -1;
	 else
		return -2;
	 }


  unsigned i=0;

  while (i < methods.size())
	 {

	 if (methods[i].parse(c) == 1)
		{
		errormessages = methods[i].geterrormessages();
		if (errormessages.empty())
		  return i;
		else
		  return -2;
		}
	 i++;
	 }

  errormessages.push_back("ERROR: unknown command\n");

  return -2;

  }

#if defined(BORLAND_OUTPUT_WINDOW)

void statobject::out(const ST::string & c,bool thick,bool italic,unsigned size,
                     int r, int g, int b,
                     bool descr)
  {
  ST::string sh = c;
  sh = sh.replaceallsigns('\n',' ');
  if (!Frame->suppoutput)
    Results->ResultsRichEdit->Lines->Append(sh.strtochar());
  if (logout->is_open())
    (*logout) << c << flush;

  if (descr==true)
	 describetext.push_back(c);
  }


void statobject::out(const vector<ST::string> & m,bool thick,bool italic,
                     unsigned size,int r, int g, int b,bool descr)
  {
  for (unsigned i=0;i<m.size();i++)
	 out(m[i],descr);
  }

#elif defined (JAVA_OUTPUT_WINDOW)

void statobject::out(const ST::string & c,bool thick,bool italic,
                     unsigned size,int r, int g, int b,bool descr)
  {
  ST::string sh = c;
  sh = sh.replaceallsigns('\n',' ');
  sh = sh+"\n";
  if (!adminb_p->get_suppressoutput())
    adminb_p->Java->CallVoidMethod(adminb_p->BayesX_obj, adminb_p->javaoutput,
    adminb_p->Java->NewStringUTF(sh.strtochar()),thick,italic,size,r,g,b);
  if (logout->is_open())
    (*logout) << c << flush;

  if (descr==true)
	 describetext.push_back(c);
  }


void statobject::out(const vector<ST::string> & m,bool thick,
                     bool italic, unsigned size,int r, int g, int b,bool descr)
  {
  for (unsigned i=0;i<m.size();i++)
	 out(m[i],thick,italic,size,r,g,b,descr);
  }

#else

void statobject::out(const ST::string & c,bool thick,bool italic,unsigned size,
                     int r, int g, int b,bool descr)
  {
  cout << c << flush;
  if (logout->is_open())
    (*logout) << c << flush;

  if (descr==true)
	 describetext.push_back(c);
  }


void statobject::out(const vector<ST::string> & m,bool thick,bool italic,
                     unsigned size,int r, int g, int b,bool descr)
  {
  for (unsigned i=0;i<m.size();i++)
	 out(m[i],descr);
  }

#endif

void statobject::outerror(const ST::string & c)
  {
  out(c,true,true,12,255,0,0);
  }

void statobject::outerror(const vector<ST::string> & c)
  {
  out(c,true,true,12,255,0,0);
  }

#if defined (JAVA_OUTPUT_WINDOW)
statobject::statobject(administrator_basic * adb,
                       const ST::string & n,const ST::string t,ofstream * lo,
                       istream * in,ST::string p)
  {
  adminb_p = adb;
  name = n;
  type = t;
  logout = lo;
  input = in;
  defaultpath = p;
  }
#else
statobject::statobject(const ST::string & n,const ST::string t,ofstream * lo,
                       istream * in,ST::string p)
  {
  name = n;
  type = t;
  logout = lo;
  input = in;
  defaultpath = p;
  }
#endif

statobject::statobject(const statobject & so)
  {
  #if defined (JAVA_OUTPUT_WINDOW)
  adminb_p = so.adminb_p;
  #endif
  name = so.name;
  type = so.type;
  describetext = so.describetext;
  errormessages = so.errormessages;
  logout = so.logout;
  input = so.input;
  defaultpath = so.defaultpath;
  newcommands = so.newcommands;

  }


const statobject & statobject::operator=(const statobject & so)
  {
  if (this == &so)
	 return *this;
  #if defined (JAVA_OUTPUT_WINDOW)
  adminb_p = so.adminb_p;
  #endif
  name = so.name;
  type = so.type;
  describetext = so.describetext;
  errormessages = so.errormessages;
  logout = so.logout;
  input = so.input;
  defaultpath = so.defaultpath;
  newcommands = so.newcommands;
  return *this;
  }


int statobject::parse(const ST::string & c)
  {
  newcommands.erase(newcommands.begin(),newcommands.end());
  return 0;
  }


int findstatobject(const vector<statobject*> & stats,const ST::string & name,
                   const ST::string & type)
  {
  int i=0;
  bool found=false;
  bool t = false;



  while ( (i<stats.size()) && (found==false) )
	 {
	 if (stats[i]->getname() == name)
		{
		found=true;
		if (stats[i]->gettype() == type)
		  t = true;
		}
	 i++;
	 }

  if (found == false)
	 return -1;
  else
	 {
	 if (t == false)
		return -2;
	 else
		return i-1;
	 }

  }


void statobject::describe(const optionlist & globaloptions)
  {

  out("\n");
  out(type + " object: " + name + "\n",true,false,16);
  out("\n");
  if (globaloptions.empty())
    {
    out("GLOBAL OPTIONS: none\n",true);
    out("\n");
    }
  else
	 {
	 out("GLOBAL OPTIONS:\n",true);
     out("\n");
	 unsigned i;
	 ST::string text;
	 for(i=0;i < globaloptions.size();i++)
       {
       text = globaloptions[i]->getname() + "= " +
			  globaloptions[i]->getValueAsString() + "\n";
       out(text);
       }
	 out("\n");
	 }

  out(describetext);
  out ("\n");

  }


#if defined(BORLAND_OUTPUT_WINDOW)
//---------------------------------------------------------------------------
#pragma package(smart_init)
#endif

