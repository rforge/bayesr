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



#include "adminparse_gnu.h"

using std::cin;
using std::flush;
using std::ios;

admin_gnu::admin_gnu(void)
  {

//  adminb = administrator_basic();
//  adminp = administrator_pointer();

  char path[100] = "";
  getcwd(path, 100);

  defaultpath = ST::string(path);

//  bool error = false;

  logfileopen = false;
  input = &cin;

  objecttyps.reserve(10);

  objecttyps.push_back("dataset");
  objecttyps.push_back("bayesreg");
  objecttyps.push_back("mcmcreg");
  objecttyps.push_back("stepwisereg");
  objecttyps.push_back("remlreg");
  objecttyps.push_back("map");
  objecttyps.push_back("dag");
  objecttyps.push_back("graph");

  delim = '\r';

  }
//------------------------------------------------------------------------------

void admin_gnu::out(const ST::string & c,
                        bool thick,bool italic,
                        unsigned size,int r,int g, int b)
  {
  ST::string sh = c;
  sh = sh.replaceallsigns('\n',' ');
  sh = sh+"\n";

  std::cout << c << flush;
  if (logfileopen)
    logout << c << flush;
  }


void admin_gnu::out(const vector<ST::string> & m,bool thick,bool italic,
                        unsigned size,int r,int g, int b)
  {
  unsigned i;
  for (i=0;i<m.size();i++)
	 out(m[i],thick,italic,size,r,g,b);
  }


void admin_gnu::outerror(const ST::string & c)
  {
  out(c,true,true,12,255,0,0);
  }

void admin_gnu::outerror(const vector<ST::string> & m)
  {
  out(m,true,true,12,255,0,0);
  }



void admin_gnu::dropobjects(ST::string name, ST::string type)
  {

  int recognized = 0;
  unsigned i=0;

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
  else if (type == "mcmcreg")
	 {
	 while ( (i < mcmcregobjects.size()) && (recognized == 0) )
		{

		if ( name == mcmcregobjects[i].getname())
		  {

		  mcmcregobjects.erase(mcmcregobjects.begin()+i,mcmcregobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 }  // end: type mcmcreg
  else if (type == "stepwisereg")
	 {
	 while ( (i < stepwiseregobjects.size()) && (recognized == 0) )
		{

		if ( name == stepwiseregobjects[i].getname())
		  {

		  stepwiseregobjects.erase(stepwiseregobjects.begin()+i,
          stepwiseregobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 }  // end: type stepwisereg
  else if (type == "remlreg")
	 {
	 while ( (i < remlregobjects.size()) && (recognized == 0) )
		{

		if ( name == remlregobjects[i].getname())
		  {

		  remlregobjects.erase(remlregobjects.begin()+i,remlregobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}
	 }  // end: type remlreg
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

  else if (type == "dag")
     {
	 while ( (i < dagobjects.size()) && (recognized == 0) )
		{

		if ( name == dagobjects[i].getname())
		  {

		  dagobjects.erase(dagobjects.begin()+i,dagobjects.begin()+i+1);
		  recognized = 1;
		  }
		i++;
		}

     }

   else if (type == "graph")
     {
     }

  adjustobjects();

  }


bool admin_gnu::alreadyexisting(const ST::string & name)
  {

  unsigned i = 0;
  bool existing = false;
  while ( (i < objects.size()) && (existing == false) )
	 {
	 if (name == objects[i]->getname())
		existing = true;
	 i++;
	 }
  return existing;
  }


ST::string admin_gnu::create(const ST::string & in)
  {

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

          else if (token[0] == "mcmcreg")
			 {
			 superbayesreg newobject(token[1],&logout,input,defaultpath,&objects);
			 mcmcregobjects.push_back(newobject);
			 }

		  else if (token[0] == "stepwisereg")
			 {
			 stepwisereg newobject(token[1],&logout,input,defaultpath,&objects);
			 stepwiseregobjects.push_back(newobject);
			 }

		  else if (token[0] == "remlreg")
			 {
			 remlreg newobject(token[1],&logout,input,defaultpath,&objects);
			 remlregobjects.push_back(newobject);
			 }

		  else if (token[0] == "map")
			 {
			 mapobject newobject(token[1],&logout,input,defaultpath,&objects);
			 mapobjects.push_back(newobject);
			 }

           else if (token[0] == "dag")
             {
             dagobject newobject(token[1],&logout,input,defaultpath,&objects);
             dagobjects.push_back(newobject);
             }

           else if (token[0] == "graph")
             {
             errormessages.push_back("ERROR: graph objects are not included in this version of BayesX\n");

//             graphobj newobject(token[1],&logout,input,&objects);
//             graphobjects.push_back(newobject);
             }

		  adjustobjects();

		  }
		return token[1];
		}
	 } // end: if token.size() == 2
  }


void admin_gnu::parseexisting(const ST::string & objectname,const ST::string & com)
  {

  int recognized = 0;
  unsigned i=0;

  while ( (i < objects.size()) && (recognized == 0) )
	 {

	 if ( objectname == objects[i]->getname())
		{
		objects[i]->parse(com);
		errormessages = objects[i]->geterrormessages();
		recognized = 1;
        if (errormessages.size() == 0)
          {
          vector<ST::string> newc =objects[i]->get_newcommands();
          unsigned i;
          for (i=0;i<newc.size();i++)
            parse(newc[i]);

          }
		}
	 i++;
	 }

  if (recognized == 0)
	 errormessages.push_back(
		"ERROR: object " + objectname + " is not existing\n");

  }


void admin_gnu::adjustobjects(void)
  {
  objects.erase(objects.begin(),objects.end());
  unsigned i;

  for (i=0;i<dataobjects.size();i++)
	 objects.push_back(&dataobjects[i]);

  for (i=0;i<bayesregobjects.size();i++)
	 objects.push_back(&bayesregobjects[i]);

  for (i=0;i<mcmcregobjects.size();i++)
	 objects.push_back(&mcmcregobjects[i]);

  for (i=0;i<stepwiseregobjects.size();i++)
	 objects.push_back(&stepwiseregobjects[i]);

  for (i=0;i<remlregobjects.size();i++)
	 objects.push_back(&remlregobjects[i]);

  for (i=0;i<mapobjects.size();i++)
	 objects.push_back(&mapobjects[i]);

  for (i=0;i<dagobjects.size();i++)
	 objects.push_back(&dagobjects[i]);

/*  for (i=0;i<graphobjects.size();i++)
    objects.push_back(&graphobjects[i]);*/

  }


bool admin_gnu::parse(ST::string & in)
  {

//  bool stop;

  errormessages.clear();

  ST::string objectname;
  ST::string firsttoken = in.getFirstToken(" ,.=");
  int pointpos = in.checksign('.');

  if (firsttoken.length() > 0 && firsttoken[0] != '%')
	 {
	 if ( (firsttoken == "quit") || (firsttoken == "exit") )
		return true;
	 else if (firsttoken == "delimiter")
		{
		vector<ST::string> token = in.strtoken(" =");
		if (token.size() != 3)
		  errormessages.push_back("ERROR: invalid syntax\n");
		else if (token[1] != "=")
		  errormessages.push_back("ERROR: \"=\" expected\n");
		else
		  {
		  if (token[2] == "return")
            {
			delim = '\r';
            }
		  else if (token[2] == ";")
            {
            delim = token[2][0];
            }
		  else
            {
            errormessages.push_back("ERROR: invalid delimiter symbol\n");
            }
		  }
		return false;
		} // end: delimiter
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

/*                stop = adminb.breakcommand();

                if(stop)
                  break;*/

				if (delim != '\r')
                  {
                  bool end = false;
                  ST::string help="";
                  in = "";
                  while ( (! infile.eof()) && (end == false) )
                    {
                    ST::getline(infile,10000,help,'\n');
                    help = help.eatwhitespace();
                    help = help.eatallcarriagereturns();
                    if (help.length() > 0 && help[0] != '%')
                      {
                      in = in + " " + help;
                      if (help[help.length()-1] == ';')
                        end = true;
                      }
                    }

                  if (in.length() > 0)
                    {
                    if (in[in.length()-1] == delim)
                      in = in.deletesign(in.length()-1);
                    in = in.replaceallsigns('\n',' ');
                    }

                  }
                else
                  ST::getline(infile,100000,in,'\n');

	            in = in.eatwhitespace();
	            in = in.eatallcarriagereturns();

                if (in.length() > 0 && in[0] != '%')
                  {
				  out("> " + in + "\n");
                  }
				parse(in);
				} // end: while (! infile.eof())

			 }

		  }

		input = &cin;

        outerror(errormessages);

        errormessages.clear();
		return false;

		}
/*     else if (firsttoken == "saveoutput")
        {

		model m;
		simpleoption replace("replace",false);
        vector<ST::string> types;
        types.push_back("rtf");
        types.push_back("txt");
        stroption type("type",types,"rtf");
		optionlist saveoutoptions;
		saveoutoptions.push_back(&replace);
        saveoutoptions.push_back(&type);
		usePathWrite uw;

		command saveoutput("saveoutput",&m,&saveoutoptions,&uw,notallowed,notallowed,
        notallowed,notallowed,optional,required);

		saveoutput.parse(in);
		errormessages = saveoutput.geterrormessages();
		if (errormessages.empty())
		  {

		  if ((replace.getvalue() == false) && (uw.isexisting() == true))
            {
            outerror("ERROR: file " + uw.getPath() +
                                    " is already existing\n");
            }
		  else
            {

#if defined(BORLAND_OUTPUT_WINDOW)
            if (type.getvalue() == "rtf")
              Results->ResultsRichEdit->PlainText = false;
            else
              Results->ResultsRichEdit->PlainText = true;

            Results->ResultsRichEdit->Lines->SaveToFile(uw.getPath().strtochar());
            Results->ResultsRichEdit->Modified = false;
            Results->ResultsRichEdit->Tag = true;
            Results->Caption = uw.getPath().strtochar();
#elif defined (JAVA_OUTPUT_WINDOW)
            jmethodID javasaveoutput = adminb.Java->GetMethodID(adminb.BayesX_cls, "JavaSaveOutput", "()V");
            adminb.Java->CallVoidMethod(adminb.BayesX_obj, javasaveoutput);
#endif

            }
          }
        else
          {
          outerror(errormessages);
          errormessages.clear();
          }
        return false;
        }*/
/*     else if (firsttoken == "clearoutput")
       {
#if defined(BORLAND_OUTPUT_WINDOW)
       Results->ResultsRichEdit->Clear();
#elif defined (JAVA_OUTPUT_WINDOW)
       jmethodID clearoutput = adminb.Java->GetMethodID(adminb.BayesX_cls, "ClearOutput", "()V");
       adminb.Java->CallVoidMethod(adminb.BayesX_obj, clearoutput);
#endif
       }*/
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
		  errormessages.push_back("ERROR: log-file is already open\n");
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
          {
		  outerror(errormessages);
          errormessages.clear();
          }
		return false;
		}  // end: logopen
	 else if (firsttoken == "logclose")
		{
		if (logfileopen == false)
		  {
		  errormessages.push_back("ERROR: currently no log-file open\n");
		  outerror(errormessages);
          errormessages.clear();
		  }
		else
		  {
		  logfileopen = false;
		  logout.close();
		  out("NOTE: log-file " + logfilepath + " closed\n");
		  }
		return false;
		}  // end: logclose
     else if (firsttoken == "defaultpath")
        {
         vector<ST::string> token = in.strtoken(" =");
         if (token.size() != 3)
           errormessages.push_back("ERROR: invalid syntax\n");
         else if (token[1] != "=")
           errormessages.push_back("ERROR: \"=\" expected\n");
         else
           {
           int test = access(token[2].strtochar(), 06);
           if(test==-1)
             {
             if(errno==ENOENT)
               errormessages.push_back("ERROR: " + token[2] + " does not exist\n");
             else if(errno==EACCES)
               errormessages.push_back("ERROR: no write access to " + token[2] + "\n");
             }
           else
             {
             defaultpath = token[2];
             bool error = false;

             #if defined(__BUILDING_LINUX)
               ST::string tempstring = defaultpath + "/temp";
             #else
               ST::string tempstring = defaultpath + "\\temp";
             #endif
               char* pathtemp = tempstring.strtochar();
               int testtemp = access(pathtemp, 06);
               if(testtemp==-1)
                 {
                 if(errno==ENOENT)
                   {
             #if defined(__BUILDING_LINUX)
                   mkdir(pathtemp, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
             #else
                   mkdir(pathtemp);
             #endif
                   errormessages.push_back("NOTE: created directory " + tempstring + "\n");
                   error=true;
                   }
                 else if(errno==EACCES)
                   {
                   errormessages.push_back("ERROR: no write access to " + tempstring + "\n");
                   error=true;
                   }
                 }

             #if defined(__BUILDING_LINUX)
               ST::string outputstring = defaultpath + "/output";
             #else
               ST::string outputstring = defaultpath + "\\output";
             #endif
               char* pathoutput = outputstring.strtochar();
               int testoutput = access(pathoutput, 00);
               if(testoutput==-1)
                 {
                 if(errno==ENOENT)
                   {
             #if defined(__BUILDING_LINUX)
                   mkdir(pathoutput, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
             #else
                   mkdir(pathoutput);
             #endif
                   errormessages.push_back("NOTE: created directory " + outputstring + "\n");
                   error=true;
                   }
                 else if(errno==EACCES)
                   {
                   errormessages.push_back("ERROR: no write access to " + outputstring + "\n");
                   error=true;
                   }
                 }
            }
          }
	    outerror(errormessages);
        errormessages.clear();
        return false;
        }
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

		  unsigned j;
		  for (j=0;j<objectnames.size();j++)
			 {

			 int recognized = 0;
			 unsigned i=0;
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

		outerror(errormessages);
        errormessages.clear();
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

		outerror(errormessages);
        errormessages.clear();
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

		outerror(errormessages);
        errormessages.clear();
		return false;

		}

	 }  // end: if (firsttoken.length() > 0)
  else                                                  // empty command
    {
    if (in.length() > 0 && in[0] != '%')
      outerror("ERROR: invalid syntax\n");
    return false;
    }

  }








