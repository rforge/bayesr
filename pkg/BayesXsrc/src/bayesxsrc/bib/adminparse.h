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



#include"../export_type.h"

#ifndef ADMINPARSE
#define ADMINPARSE

#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include "StatReview.h"
#include<StatwinFrame.h>
#include<statwin_haupt.h>

#elif defined(JAVA_OUTPUT_WINDOW)
#include<jni.h>
#endif

//#include <vcl.h>
#if !defined(__BUILDING_GNU)
#include <FileCtrl.hpp>
#include<dir.h>
#endif

#include<fstream>
#include"adminparse_basic.h"
#include"adminparse_pointers.h"
#include"data.h"
#include"statobj.h"
#include"dataobj.h"
#include"bayesreg.h"
#include"remlreg.h"
#include"mapobject.h"
#include"dagobject.h"
#include"graphobj.h"
#include"stepwisereg.h"
//------------------------------------------------------------------------------

using namespace std;

class __EXPORT_TYPE administrator
  {

  private:

  //------------------------ PRIVATE VARIABLES ---------------------------------


  ST::string cmd;          // command

  char delim;              // delimiter

  ofstream logout;

  istream * input;

  ST::string defaultpath;

  bool logfileopen;

  ST::string logfilepath;

  // contains (valid) objecttyps
  // valid types:
  // - dataset
  // - bayesreg
  // - remlreg
  // - stepwisereg
  // - map
  // - dag
  // - graph

  vector<ST::string> objecttyps;

  // contains pointers to current (stat-)objects

  vector<statobject*> objects;

  // contains current errmormessages

//  errorm::messages errormessages;
  vector<ST::string> errormessages;

  // 'dataobjects' contains current dataobjects

  vector<dataobject> dataobjects;

  // 'bayesregobjects' contains current bayesreg objects

  vector<bayesreg> bayesregobjects;

  // 'stepwiseregobjects' contains current stepwisereg objects

  vector<stepwisereg> stepwiseregobjects;

  // 'remlregobjects' contains current remlreg objects

  vector<remlreg> remlregobjects;


  // 'mapobjects' contains current map objects

  vector<mapobject> mapobjects;

  // 'dagobjects' contains current dag objects

  vector<dagobject> dagobjects;

  // 'graphobj' contains current graph objects

  vector<graphobj> graphobjects;


  //------------------------ PRIVATE FUNCTIONS ---------------------------------

  void out(const ST::string & c,bool thick=false,bool italic = false,
           unsigned size = 12,int r=0,int g=0, int b=0);

  void out(const vector<ST::string> & m,bool thick=false,bool italic = false,
           unsigned size = 12,int r=0,int g=0, int b=0);

  void outerror(const ST::string & c);

  void outerror(const vector<ST::string> & m);


  // FUNCTION: alreadyexisting
  // TASK: returns 'true', if object with objectname 'name' is already existing

  bool alreadyexisting(const ST::string & name);

  // FUNCTION: create

  ST::string create(const ST::string & in);

  void adjustobjects(void);

  // FUNCTION: drop

  void dropobjects(ST::string name,ST::string type);


  // FUNCTION: parseexisting
  // TASK: parses command 'com' for object with name 'objectname'
  //       objectname should be an object, that is still existing
  // POSSIBLE ERRORS:
  // - object with name 'objectname' is not existing
  // - command 'com' is invalid (i.e. contains errors) (depending on the
  //   special structure of object with name 'objectname')

  void parseexisting(const ST::string & objectname,const ST::string & com);

  void parsespecial(const ST::string & com);


  public:

  administrator_basic adminb;
  administrator_pointer adminp;

  // CONSTRUCTOR

  administrator(void);

  // DESTRUCTOR

  ~administrator() {}

  bool parse(ST::string & in);

  char get_delim(void)
    {
    return delim;
    }

  vector<ST::string> & get_objecttype(void)
    {
    return objecttyps;
    }

  vector<statobject*> & get_objects(void)
    {
    return objects;
    }

  vector<mapobject> & get_mapobjects(void)
    {
    return mapobjects;
    }

  vector<dataobject> & get_dataobjects(void)
    {
    return dataobjects;
    }



  };


#endif
