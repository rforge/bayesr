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


// DATE: 16.01.98


#if !defined (ADMINISTRATOR_INCLUDED)

#define ADMINISTRATOR_INCLUDED

#include<fstream>
#include<errorm.h>
#include<data.h>
#include<statobj.h>
#include<dataobj.h>
#include<bayesreg.h>
#include<map.h>


// HINZUFÜGEN EINES NEUEN OBJEKTTYPS
// 1. Neue header Datei mit include hinzufügen
// 2. Im default constructor den neuen Objekttyp initialisieren
// 3. Vektor definieren, der die neuen Objekttypen speichern kann (im private
//    Teil von administrator)
// 4. In der Funktion 'create' den neuen Objekttyp initialisieren
// 5. function 'dropobjects' abändern
// 6. function 'adjustobjects' abändern


class administrator
  {

  private:


  //------------------------ PRIVATE VARIABLES ---------------------------------

  char delim;              // sign that indicates the end of a command

  ofstream logout;

  istream * input;

  ST::string defaultpath;

  bool logfileopen;

  ST::string logfilepath;

  // contains (valid) objecttyps
  // valid types:
  // - dataset
  // - bayesreg

  vector<ST::string> objecttyps;

  // contains pointers to current (stat-)objects

  vector<statobject*> objects;

  // contains current errmormessages

  errorm::messages errormessages;

  // 'dataobjects' contains current dataobjects

  vector<dataobject> dataobjects;

  // 'bayesregobjects' contains current bayesreg objects

  vector<bayesreg> bayesregobjects;

  // 'mapobjects' contains current map objects

  vector<map> mapobjects;


  //------------------------ PRIVATE FUNCTIONS ---------------------------------

  void out(const ST::string & c);

  void out(const vector<ST::string> & m);

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

  bool parse(ST::string & in);

  public:

  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR

  administrator(void)
	 {
	 ST::string line;
	 ifstream fin("statprog.ini");
	 ST::getline(fin,line);
	 defaultpath = line;
	 fin.close();
	 logfileopen = false;
	 input = &cin;
	 objecttyps.push_back("dataset");
	 objecttyps.push_back("bayesreg");
	 objecttyps.push_back("map");
    delim = '\n';
	 }

  // DESTRUCTOR

  ~administrator() {}

  void run(void);


  };


#endif

