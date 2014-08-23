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


#include<admin.h>

using std::flush;
using std::ios;

void administrator::out(const ST::string & c)
  {
  cout << c << flush;
  if (logout.is_open())
    logout << c << flush;
  }


void administrator::out(const vector<ST::string> & m)
  {
  for (int i=0;i<m.size();i++)
	 out(m[i]);
  }


void administrator::dropobjects(ST::string name, ST::string type)
  {

  int recognized = 0;
  int i=0;

  if (type == "dataset")
	 {
	 while ( (i < dataobjects.size()) && (recognized == 0) )
		{

		if ( name == dataobjects[i].getname())
		  {

		  dataobjects.erase(dataobjects.begin()+i,dataobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 } // end: type dataset
  else if (type == "bayesreg")
	 {
	 while ( (i < bayesregobjects.size()) && (recognized == 0) )
		{

		if ( name == bayesregobjects[i].getname())
		  {

		  bayesregobjects.erase(bayesregobjects.begin()+i,bayesregobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 }  // end: type bayesreg
  else if (type == "map")
	 {
	 while ( (i < mapobjects.size()) && (recognized == 0) )
		{

		if ( name == mapobjects[i].getname())
		  {

		  mapobjects.erase(mapobjects.begin()+i,mapobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 } // end: type map

  adjustobjects();

  }


bool administrator::alreadyexisting(const ST::string & name)
  {

  int i = 0;
  bool existing = false;
  while ( (i < objects.size()) && (existing == false) )
	 {
	 if (name == objects[i]->getname())
		existing = true;
	 i++;
	 }
  return existing;
  }


ST::string administrator::create(const ST::string & in)
  {

//  int test;



  ST::string name;
  vector<ST::string> token = in.strtoken(" ");

  if (token.size() != 2)
	 {
	 errormessages.push_back("ERROR: invalid creation of new object\n");
	 return "";
	 }
  else // token.size() == 2
	 {

	 if (token[1].isvarname() == 1)
		{
		errormessages.push_back("ERROR: invalid object name\n");
		return "";
		}
	 else
		{
		if (alreadyexisting(token[1]) == true)
		  errormessages.push_back(
		  "ERROR: object " + token[1] + " is already existing\n");
		else
		  {

		  if (token[0] == "dataset")
			 {
			 dataobject newobject(token[1],&logout,input);
			 dataobjects.push_back(newobject);
			 }
		  else if (token[0] == "bayesreg")
			 {
			 bayesreg newobject(token[1],&logout,input,defaultpath,&objects);
			 bayesregobjects.push_back(newobject);
			 }
		  else if (token[0] == "map")
			 {
			 MAP::map newobject(token[1],&logout,input,defaultpath);
			 mapobjects.push_back(newobject);
			 }


		  adjustobjects();


		  }
		return token[1];
		}
	 } // end: if token.size() == 2
  }


void administrator::parseexisting(const ST::string & objectname,const ST::string & com)
  {

  int recognized = 0;
  int i=0;


  while ( (i < objects.size()) && (recognized == 0) )
	 {

	 if ( objectname == objects[i]->getname())
		{
		objects[i]->parse(com);
		errormessages = objects[i]->geterrormessages();
		recognized = 1;
		}
	 i++;
	 }

  if (recognized == 0)
	 errormessages.push_back(
		"ERROR: object " + objectname + " is not existing\n");

  }


void administrator::adjustobjects(void)
  {
  objects.erase(objects.begin(),objects.end());
  int i;

  for (i=0;i<dataobjects.size();i++)
	 objects.push_back(&dataobjects[i]);

  for (i=0;i<bayesregobjects.size();i++)
	 objects.push_back(&bayesregobjects[i]);

  for (i=0;i<mapobjects.size();i++)
	 objects.push_back(&mapobjects[i]);


  }


bool administrator::parse(ST::string & in)
  {

  errormessages.clear();

  ST::string objectname;
  ST::string firsttoken = in.getFirstToken(" .");
  int pointpos = in.checksign('.');

  if (firsttoken.length() > 0)
	 {
	 if ( (firsttoken == "quit") || (firsttoken == "exit") )
		return true;
	 else if (firsttoken == "delimeter")
		{
		vector<ST::string> token = in.strtoken(" ");
		if (token.size() != 3)
		  errormessages.push_back("ERROR: invalid syntax\n");
		else if (token[1] != "=")
		  errormessages.push_back("ERROR: \"=\" expected\n");
		else
		  {
		  if (token[2] == "newline")
			 delim = '\n';
		  else if (token[2].length() > 1)
			 errormessages.push_back("ERROR: invalid delimeter symbol\n");
		  else
			 delim = token[2][0];
		  }

		return false;
		} // end: delimeter
	 else if (firsttoken == "usefile")
		{

		vector<ST::string> token = in.strtoken(" ");
		if (token.size() < 2)
		  errormessages.push_back("ERROR: filename expected\n");
		else if (token.size() > 2)
		  errormessages.push_back("ERROR: invalid syntax\n");

		if (errormessages.empty())
		  {
		  ST::string path = token[1];
		  if (path.isexistingfile() == 1)
			 errormessages.push_back("ERROR: file " + path +
											 " could not be opened\n");
		  else
			 {
			 ST::string in;
			 ifstream infile;
			 input = &infile;
			 ST::open(infile,path);
			 while (! infile.eof())
				{
				ST::getline(infile,10000,in,delim);
				if (delim != '\n')
				 in = in.replaceallsigns('\n',' ');
	         in = in.eatwhitespace();


				out("> " + in + "\n");
				parse(in);
				}

			 }
		  }

		input = &cin;
		out(errormessages);
		return false;

		}
	 else if (firsttoken == "logopen")
		{
		model m;
		simpleoption replace("replace",false);
		optionlist logoptions;
		logoptions.push_back(&replace);
		usePathWrite uw;

		command logopen("logopen",&m,&logoptions,&uw,notallowed,notallowed,
							 notallowed,notallowed,optional,required);

		logopen.parse(in);
		errormessages = logopen.geterrormessages();
		if (logfileopen == true)
		  errormessages.push_back("ERROR: logfile is already open\n");
		if (errormessages.empty())
		  {
		  logfileopen = true;
		  logfilepath = uw.getPath();
		  if ((replace.getvalue() == false) && (uw.isexisting() == true))
            {
            ST::open(logout,logfilepath,ios::app);
            }
		  else
            {
            ST::open(logout,logfilepath);
            }
		  }
		else
		  out(errormessages);
		return false;
		}  // end: logopen
	 else if (firsttoken == "logclose")
		{
		if (logfileopen == false)
		  {
		  errormessages.push_back("ERROR: currently no logfile open\n");
		  out(errormessages);
		  }
		else
		  {
		  logfileopen = false;
		  logout.close();
		  out("NOTE: logfile " + logfilepath + " closed\n");
		  }
		return false;
		}  // end: logclose
	 else if (firsttoken == "drop")
		{

		modelStandard m;
		optionlist dropoptions;
		usePathWrite uw;

		command drop("drop",&m,&dropoptions,&uw,required,notallowed,
							 notallowed,notallowed,notallowed,notallowed);

		drop.parse(in);
		errormessages = drop.geterrormessages();
		vector<ST::string> objectnames =  m.getModelVarnamesAsVector();
		if (objectnames.size() == 0)
		  errormessages.push_back("ERROR: objectlist required\n");

		if (errormessages.empty())
		  {

		  int j;
		  for (j=0;j<objectnames.size();j++)
			 {

			 int recognized = 0;
			 int i=0;
			 while ( (i < objects.size()) && (recognized == 0) )
				{

				if ( objectnames[j] == objects[i]->getname())
				  {
				  ST::string type = objects[i]->gettype();
				  dropobjects(objectnames[j],type);
				  recognized = 1;
				  }
				i++;
				}

			  if (recognized == 0)
				 errormessages.push_back(
				 "ERROR: object " + objectnames[j] + " is not existing\n");

			 } // end: for (j=0;j<objectnames.size();j++)

		  }  // end: if (errormessages.empty())

		out(errormessages);
		return false;
		} // end: drop
	 else if (firsttoken.isinlist(objecttyps) >= 0)      // create a new object
		{


		if (pointpos == -1)
		  objectname = create(in);
		else
		  objectname = create(in.substr(0,pointpos));


		if ( (errormessages.empty()) && (pointpos > 0) )
		  {
		  if (in.length()-1-pointpos <= 0)
			 errormessages.push_back("ERROR: invalid syntax\n");
		  else
			parseexisting(objectname,in.substr(pointpos+1,in.length()-1-pointpos));
		  }

		out(errormessages);
		return false;

		}               // end: create a new object
	 else     // existing object
		{
		if (pointpos != firsttoken.length())
		  errormessages.push_back("ERROR: invalid syntax\n");
		else
		  if (in.length() > pointpos+1)
			parseexisting(firsttoken,in.substr(pointpos+1,in.length()-pointpos-1));
		  else
			 errormessages.push_back("ERROR: invalid syntax\n");

		out(errormessages);
		return false;

		}

	 }  // end: if (firsttoken.length() > 0)
  else                                                  // empty command
	 return false;

  }


void administrator::run(void)
  {

  ST::string inp;
  ST::string h;
  bool stop = false;

  while (stop == false)
	 {
	 cout << "> ";
	 ST::getline (cin,10000,inp,delim);
	 if (delim != '\n')
		inp = inp.replaceallsigns('\n',' ');
     inp = inp.eatwhitespace();

     if (logout.is_open())
	   logout << "> " << inp << endl;

	 stop = parse(inp);
	 }

  bayesregobjects.erase(bayesregobjects.begin(),bayesregobjects.end());

  }


