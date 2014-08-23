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
#include <vcl.h>
#pragma hdrstop

#include <describe_map.h>
#endif


#include"mapobject.h"


void mapobject::create(void)
  {

  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  // method infile

  mod = modelStandard();
  uread = usePathRead();

  weightdefs = vector<ST::string>(3);
  weightdefs[0] = "adjacency";
  weightdefs[1] = "centroid";
  weightdefs[2] = "combnd";
  weightdef = stroption("weightdef",weightdefs,"adjacency");

  neighbors = simpleoption("neighbors",false);

  graf = simpleoption("graph",false);
  centroids = simpleoption("centroids",false);

  infileoptions.push_back(&weightdef);
  infileoptions.push_back(&neighbors);
  infileoptions.push_back(&graf);
  infileoptions.push_back(&centroids);

  methods.push_back(command("infile",&mod,&infileoptions,&uread,notallowed,
						  notallowed,notallowed,notallowed,optional,required));

  functions[0] = infilerun;

  // method computeneighbors

  modn = modelStandard();
  uwrite = usePathWrite();
  replace = simpleoption("replace",false);
  computeneighborsoptions.push_back(&replace);

  methods.push_back(command("computeneighbors",&modn,&computeneighborsoptions,
									 &uwrite,notallowed,notallowed,notallowed,notallowed,
									 optional,optional));

  functions[1] = computeneighborsrun;


  // method createmap

  modcreate = modelStandard();
  udata = use();

  weightdefscr = vector<ST::string>(2);
  weightdefscr[0] = "adjacency";
  weightdefscr[1] = "centroid";
  weightdefcr = stroption("weightdef",weightdefscr,"adjacency");

  maxdif  = doubleoption("maxdif",1,0.0000001,100000000);

  createmapoptions.push_back(&weightdefcr);
  createmapoptions.push_back(&maxdif);

  methods.push_back(command("createmap",&modcreate,&createmapoptions,&udata,
                            required,notallowed,notallowed,notallowed,optional,
                            required));

  functions[2] = createmaprun;

  // method outfile

  outfileoptions.reserve(5);

  modo = modelStandard();
  uwriteo = usePathWrite();
  grafo = simpleoption("graph",false);
  centroidso = simpleoption("centroids",false);
  replaceo = simpleoption("replace",false);
  includeweights = simpleoption("includeweights",false);

  outfileoptions.push_back(&grafo);
  outfileoptions.push_back(&centroidso);
  outfileoptions.push_back(&replaceo);
  outfileoptions.push_back(&includeweights);

  methods.push_back(command("outfile",&modo,&outfileoptions,&uwriteo,notallowed,
						  notallowed,notallowed,notallowed,optional,required));

  functions[3] = outfilerun;

  // method reorder

  modre = modelStandard();
  udatare = use();

  methods.push_back(command("reorder",&modre,&reorderoptions,&udatare,notallowed,
						  notallowed,notallowed,notallowed,notallowed,notallowed));

  functions[4] = reorderrun;


  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  }

#if defined(JAVA_OUTPUT_WINDOW)
mapobject::mapobject(administrator_basic * adb, administrator_pointer * adp,
const ST::string & n,ofstream * lo,istream * in,ST::string p,
                     vector<statobject*> * st)
	                 : statobject(adb,n,"map",lo,in,p)
  {
  adminp_p = adp;
  statobj = st;
  create();
  mapexisting = false;
  describetext.push_back("Number of polygones currently in memory: none \n");
  }
#else
mapobject::mapobject(const ST::string & n,ofstream * lo,istream * in,ST::string p,
                     vector<statobject*> * st)
	                 : statobject(n,"map",lo,in,p)
  {
  statobj = st;
  create();
  mapinfo = MAP::map();
  mapexisting = false;
  describetext.push_back("Number of polygones currently in memory: none \n");
  }
#endif


mapobject::mapobject(const mapobject & m) : statobject(statobject(m))
  {
  statobj = m.statobj;
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = m.adminp_p;
  #endif
  mapinfo = m.mapinfo;
  mapexisting = m.mapexisting;
  }


const mapobject & mapobject::operator=(const mapobject & m)
  {
  if (this == &m)
	 return *this;
  statobject::operator=(statobject(m));
  statobj = m.statobj;
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = m.adminp_p;
  #endif
  mapexisting = m.mapexisting;
  mapinfo = m.mapinfo;
  return *this;
  }


int mapobject::parse(const ST::string & c)
  {
  optionlist globaloptions = optionlist();
  int pos = statobject::parsecom(c,methods,globaloptions);

  if (pos >= 0)
	 (*functions[pos])(*this);

  return(pos);
  }


void infilerun(mapobject & m)
  {

  vector<ST::string> errormessages;

  m.mapexisting = false;
  m.describetext.erase(m.describetext.begin(),m.describetext.end());

  ST::string path = m.uread.getPath();

  MAP::metric weightmode;
  if (m.weightdef.getvalue() == "adjacency")
    weightmode = MAP::adjacent;
  else if (m.weightdef.getvalue() == "centroid")
    weightmode = MAP::centroid;
  else
    weightmode = MAP::combnd;

  if (m.graf.getvalue() == true)
    {
    graph g(path);
    if (g.geterror() == true)
      errormessages.push_back(g.geterrormessage());
    else
      {
      #if defined(JAVA_OUTPUT_WINDOW)
      m.mapinfo = MAP::map(m.adminb_p,g);
      #else
      m.mapinfo = MAP::map(g);
      #endif

      if (g.get_nrgraphs() > 1)
        m.out("NOTE: The graph is disconnected (" +
        ST::inttostring(g.get_nrgraphs()) + " parts).\n");
      }
    }
  else if (m.centroids.getvalue() == true)
    {
    #if defined(JAVA_OUTPUT_WINDOW)
    m.mapinfo = MAP::map(m.adminb_p,path);
    #else
    m.mapinfo = MAP::map(path);
    #endif
    }
  else
    {
    #if defined(JAVA_OUTPUT_WINDOW)
    m.mapinfo = MAP::map(m.adminb_p,path,weightmode);
    #else
    m.mapinfo = MAP::map(path,weightmode);
    #endif
    }

  if (errormessages.size() == 0)
    {
    errormessages = m.mapinfo.get_errormessages();

    if (errormessages.size() == 0)
	   {
	   m.out("NOTE: " + ST::inttostring(m.mapinfo.get_nrregions())
              + " regions read from file " + path + "\n");
	   m.mapexisting = true;
	   m.describetext.push_back("Number of regions currently in memory: " +
 						  ST::inttostring(m.mapinfo.get_nrregions()) + "\n");
       }
    else
	   {
       m.outerror(errormessages);
	   m.describetext.erase(m.describetext.begin(),m.describetext.end());
	   m.describetext.push_back("Number of polygones currently in memory: none \n");
	   }

    }
  else
    {
    m.outerror(errormessages);
    m.describetext.erase(m.describetext.begin(),m.describetext.end());
    m.describetext.push_back("Number of polygones currently in memory: none \n");
    }


  }


void createmaprun(mapobject & m)
  {

  unsigned i;

  bool failure=false;

  m.mapexisting = false;
  m.describetext.erase(m.describetext.begin(),m.describetext.end());

  MAP::metric weightmode;
  if (m.weightdefcr.getvalue() == "adjacency")
    weightmode = MAP::adjacent;
  else
    weightmode = MAP::centroid;

  datamatrix X;

  vector<ST::string> modelvarnamesv;

  //--------------------- reading dataset information --------------------------

  dataobject * datap;               // pointer to dataset

  int objpos = findstatobject(*m.statobj,m.udata.getusingtext(),"dataset");
  statobject * s;
  if (objpos >= 0)
    {
    s = m.statobj->at(objpos);
    datap = dynamic_cast<dataobject*>(s);
    }
  else
    {
    if (objpos == -1)
      m.outerror("ERROR: " + m.udata.getusingtext() + " is not existing\n");
    else
      m.outerror("ERROR: " + m.udata.getusingtext() + " is not a dataset object\n");
    failure = true;
    }

  //------------------ end: reading dataset information ------------------------


  //---------------- reading data, creating designmatrices ---------------------

  if (failure == false)
  {

  modelvarnamesv = m.modcreate.getModelVarnamesAsVector();

  if (modelvarnamesv.size() != 2)
    {
    m.outerror("ERROR: number of variables must be two\n");
    failure = true;
    }

  // testing, wether all variables specified are already existing
  vector<ST::string> notex;
  if ((datap->allexisting(modelvarnamesv,notex)) == false)
    {
    for (i=0;i<notex.size();i++)
      m.outerror("ERROR: variable " + notex[i] + " is not existing\n");
    failure = true;
    } // end: if ((datap->allexisting(modelvarnamesv,notex)) == false)

  }

  if (failure == false)
    {
    datap->makematrix(modelvarnamesv,X);

    m.errormessages = datap->geterrormessages();
    if (!m.errormessages.empty())
      {
      failure = true;
      }
    }

  if (failure == false)
    {
    #if defined(JAVA_OUTPUT_WINDOW)
    m.mapinfo = MAP::map(m.adminb_p,X,m.maxdif.getvalue(),weightmode);
    #else
    m.mapinfo = MAP::map(X,m.maxdif.getvalue(),weightmode);
    #endif
//    m.errormessages = m.mapinfo.get_errormessages();
//    if (!m.errormessages.empty())
//      {
//      failure = true;
//      }

    }

  if (failure == false)
    {
    m.mapexisting = true;
    m.out("\n");
    m.out("NOTE: map created\n");
    m.out("Total number of sites: " + ST::inttostring(m.mapinfo.get_nrregions())
           + "\n");
    m.out("Minimum number of neighbors: " +
          ST::inttostring(m.mapinfo.get_minn()) + "\n");
    m.out("Maximum number of neighbors: " +
          ST::inttostring(m.mapinfo.get_maxn()) + "\n");
    m.out("Minimum distance between two neighbors: " +
          ST::doubletostring(m.mapinfo.get_mindistance(),6) + "\n");
    m.out("Maximum distance between two neighbors: " +
          ST::doubletostring(m.mapinfo.get_maxdistance(),6) + "\n");
    m.out("\n");

    }

  }


void computeneighborsrun(mapobject & m)
  {

/*
  ST::string path = m.uwrite.getPath();

  if ( (path.length() > 0) &&  (m.uwrite.isexisting() == true) && (m.replace.getvalue() == false) )
	 m.out("ERROR: file " + path + " is already existing\n");
  else
	 {
	 if (m.mapexisting == true)
       {

       m.mapinfo.computeneighbors(path);

		m.out("NOTE: minimum number of neighbors is " + ST::inttostring(m.mapinfo.getminn()) +
			 ",\n      maximum number of neighbors is " + ST::inttostring(m.mapinfo.getmaxn()) + "\n" );

       } // end: if m.mapexisting == true;
	 else  // end: if (mapexisting == true)
		m.out("ERROR: no map currently in memormy\n");
	 }

  */

  }


void outfilerun(mapobject & m)
  {
  if (m.mapexisting)
    {
    if (m.grafo.getvalue()==true)
      {
      ST::string path = m.uwriteo.getPath();
      if ( (m.uwriteo.isexisting() == true) && (m.replaceo.getvalue() == false) )
	    m.errormessages.push_back(
        "ERROR: file " + path + " is already existing\n");
      else
        {
        if (m.includeweights.getvalue() == true)
          m.mapinfo.outgraph(path,true);
        else
          m.mapinfo.outgraph(path);
        }
      }
    else if (m.centroidso.getvalue()==true)
      {
      ST::string path = m.uwriteo.getPath();
      if ( (m.uwriteo.isexisting() == true) && (m.replaceo.getvalue() == false) )
	    m.errormessages.push_back(
        "ERROR: file " + path + " is already existing\n");
      else
        {
        m.mapinfo.outcentroids(path);
        }
      }
    else
      {
      ST::string path = m.uwriteo.getPath();
      if ( (m.uwriteo.isexisting() == true) && (m.replaceo.getvalue() == false) )
	    m.errormessages.push_back(
        "ERROR: file " + path + " is already existing\n");
      else
        {
        if (m.mapinfo.polygones_existing())
          m.mapinfo.outmap(path);
        else
          {
          if (m.includeweights.getvalue() == true)
            m.mapinfo.outgraph(path,true);
          else
            m.mapinfo.outgraph(path);
          m.out("NOTE: polygones of map did not exist\n");
          m.out("graph file was written to file " + path + " instead\n");
          m.out("\n");
          }

        }
      }

    } // end: if (m.mapexisting)
  else
    m.out("NOTE: map object does not contain any data\n");

  }


void reorderrun(mapobject & m)
  {

  if (m.mapexisting)
    {
    unsigned bandsizeold=m.mapinfo.get_bandsize();
    m.mapinfo.reorderopt();
    if (m.mapinfo.get_errormessages().size() == 0)
      {
      unsigned bandsizenew=m.mapinfo.get_bandsize();
      m.out("NOTE: bandsize changed from " + ST::inttostring(bandsizeold) +
            " to " + ST::inttostring(bandsizenew) + "\n");
      }
    else
      m.out(m.mapinfo.get_errormessages());
    }
  else
    m.out("NOTE: map object does not contain any data\n");

  }



void mapobject::describe(const optionlist & globaloptions)
  {
  if(mapexisting)
    {
#if defined (BORLAND_OUTPUT_WINDOW)
    mapform->mapinfo = &mapinfo;
    mapform->mapname = getname();
#elif defined (JAVA_OUTPUT_WINDOW)
    adminp_p->set_mapinfo(&mapinfo);
#endif
    out("\n");
    out("\n");
    out("MAP " + name + "\n",true,false,16);
    out("Number of regions: " + ST::inttostring(mapinfo.get_nrregions()) + "\n");
    out("Minimum number of neighbors: " + ST::inttostring(mapinfo.get_minn()) + "\n");
    out("Maximum number of neighbors: " + ST::inttostring(mapinfo.get_maxn()) + "\n");
    out("Bandsize of corresponding adjacency matrix: " + ST::inttostring(mapinfo.get_bandsize()) + "\n");
    out("\n");
    if (mapinfo.polygones_existing())
      {
#if defined (BORLAND_OUTPUT_WINDOW)
      mapform->ShowModal();
#elif defined (JAVA_OUTPUT_WINDOW)

      int fontsize = 0;
      ST::string psname="";
      jmethodID javashowmap = adminb_p->Java->GetMethodID(adminb_p->BayesX_cls, "JavaDescribeMap", "(Z)V");
      adminb_p->Java->CallVoidMethod(adminb_p->BayesX_obj, javashowmap, false);

//      jmethodID javashowmap = adminb_p->Java->GetMethodID(
//      adminb_p->BayesX_cls, "JavaShowMap", "(Z)V");
//      adminb_p->Java->CallVoidMethod(adminb_p->BayesX_obj, javashowmap, false);
#endif
      }
    else
      out("NOTE: map object does not contain any data\n");
    }
  else
    out("NOTE: map object does not contain any data\n");
  }


#if defined(BORLAND_OUTPUT_WINDOW)
//---------------------------------------------------------------------------
#pragma package(smart_init)
#endif





