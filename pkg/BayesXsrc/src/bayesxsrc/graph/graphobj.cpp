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
#endif


#include"graphobj.h"


void graphobj::create (void)
  {

  m = modelStandard();
  e = expression();
  u = use();

// SYNTAX OF COMMANDS:
// name [model] [weight varname] [by varname] [if expression]
//      [, options] [using usingtext]

  // method drawmap

  mapname = stroption("map","m");
  psname = stroption("outfile");
  title = stroption("title");
//  upperlimit = doubleoption("upperlimit",1.0E7,-1.0E7,1.0E7);
//  lowerlimit = doubleoption("lowerlimit",-1.0E7,-1.0E7,1.0E7);
  upperlimit = doubleoption("upperlimit",1,-MAXDOUBLE,MAXDOUBLE);
  lowerlimit = doubleoption("lowerlimit",0,-MAXDOUBLE,MAXDOUBLE);
  nrcolors = intoption("nrcolors",256,1,256);
  color = simpleoption("color",false);
  nolegend = simpleoption("nolegend",false);
  swapcolors = simpleoption("swapcolors",false);
  drawnames = simpleoption("drawnames",false);
  hclcolors = simpleoption("hcl",false);
  replace = simpleoption("replace",false);
  pcat = simpleoption("pcat",false);
  fontsize = intoption("fontsize",12,0,100);
  titlescale = doubleoption("titlesize",1.5,0.0,MAXDOUBLE);

  drawmapoptions.push_back(&mapname);
  drawmapoptions.push_back(&psname);
  drawmapoptions.push_back(&title);
  drawmapoptions.push_back(&lowerlimit);
  drawmapoptions.push_back(&upperlimit);
  drawmapoptions.push_back(&nrcolors);
  drawmapoptions.push_back(&color);
  drawmapoptions.push_back(&nolegend);
  drawmapoptions.push_back(&swapcolors);
  drawmapoptions.push_back(&drawnames);
  drawmapoptions.push_back(&hclcolors);
  drawmapoptions.push_back(&replace);
  drawmapoptions.push_back(&pcat);
  drawmapoptions.push_back(&fontsize);
  drawmapoptions.push_back(&titlescale);

  methods.push_back(command("drawmap",&m,&drawmapoptions,&u,
//              required,notallowed,notallowed,optional,optional,required));
			 optional,notallowed,notallowed,notallowed,required,optional));

  functions[0] = drawmaprun;

  // method plotnonp

//  vector<ST::string> admiss;
//  admiss.push_back("line");
//  admiss.push_back("points");

  xlab = stroption("xlab");
  ylab = stroption("ylab");
//  connect = stroption("connect",admiss,"line");
  connect = stroption("connect");
  height = intoption("height",210,0,500);
  width = intoption("width",356,0,500);
  xlimtop = doubleoption("xlimtop",-MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xlimbottom = doubleoption("xlimbottom",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  ylimtop = doubleoption("ylimtop",-MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  ylimbottom = doubleoption("ylimbottom",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xstep = doubleoption("xstep",0.0,-MAXDOUBLE,MAXDOUBLE);
  xstart = doubleoption("xstart",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  ystep = doubleoption("ystep",0.0,-MAXDOUBLE,MAXDOUBLE);
  ystart = doubleoption("ystart",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  year = intoption("year",0,0,3000);
  month = intoption("month",0,1,12);
  linewidth = intoption("linewidth",5,0,100);
//  fontsize = intoption("fontsize",12,0,100);
  pointsize = intoption("pointsize",20,0,100);
  linecolor = stroption("linecolor");

  plotnonpoptions.push_back(&psname);
  plotnonpoptions.push_back(&title);
  plotnonpoptions.push_back(&xlab);
  plotnonpoptions.push_back(&ylab);
  plotnonpoptions.push_back(&connect);
  plotnonpoptions.push_back(&replace);
  plotnonpoptions.push_back(&height);
  plotnonpoptions.push_back(&width);
  plotnonpoptions.push_back(&xlimtop);
  plotnonpoptions.push_back(&xlimbottom);
  plotnonpoptions.push_back(&ylimtop);
  plotnonpoptions.push_back(&ylimbottom);
  plotnonpoptions.push_back(&xstep);
  plotnonpoptions.push_back(&xstart);
  plotnonpoptions.push_back(&ystep);
  plotnonpoptions.push_back(&ystart);
  plotnonpoptions.push_back(&year);
  plotnonpoptions.push_back(&month);
  plotnonpoptions.push_back(&linewidth);
  plotnonpoptions.push_back(&fontsize);
  plotnonpoptions.push_back(&pointsize);
  plotnonpoptions.push_back(&linecolor);
  plotnonpoptions.push_back(&titlescale);

  methods.push_back(command("plot",&m,&plotnonpoptions,&u,
			 required,notallowed,notallowed,optional,optional,required));

  functions[1] = plotnonprun;

  plotsampleoptions.push_back(&psname);
  plotsampleoptions.push_back(&connect);
  plotsampleoptions.push_back(&replace);

  methods.push_back(command("plotsample",&m,&plotsampleoptions,&u,
                    notallowed,notallowed,notallowed,optional,optional,required));

  functions[2] = plotsamplerun;

  mean = simpleoption("mean",false);

  plotautocoroptions.push_back(&psname);
  plotautocoroptions.push_back(&connect);
  plotautocoroptions.push_back(&replace);
  plotautocoroptions.push_back(&mean);

  methods.push_back(command("plotautocor",&m,&plotautocoroptions,&u,
                    notallowed,notallowed,notallowed,optional,optional,required));

  functions[3] = plotautocorrun;

  zlab = stroption("zlab");
  zlimtop = doubleoption("zlimtop",-MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  zlimbottom = doubleoption("zlimbottom",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  zstep = doubleoption("zstep",0.0,-MAXDOUBLE,MAXDOUBLE);
  zstart = doubleoption("zstart",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xrot = doubleoption("xrot",0.0,-360,360);
  yrot = doubleoption("yrot",0.0,-360,360);
  zrot = doubleoption("zrot",0.0,-360,360);
  gridsize = intoption("gridsize",30,0,500);

  plotsurfoptions.push_back(&psname);
  plotsurfoptions.push_back(&replace);
  plotsurfoptions.push_back(&title);
  plotsurfoptions.push_back(&xlab);
  plotsurfoptions.push_back(&ylab);
  plotsurfoptions.push_back(&zlab);
  plotsurfoptions.push_back(&xrot);
  plotsurfoptions.push_back(&yrot);
  plotsurfoptions.push_back(&zrot);
  plotsurfoptions.push_back(&height);
  plotsurfoptions.push_back(&width);
  plotsurfoptions.push_back(&xlimtop);
  plotsurfoptions.push_back(&ylimtop);
  plotsurfoptions.push_back(&zlimtop);
  plotsurfoptions.push_back(&xlimbottom);
  plotsurfoptions.push_back(&ylimbottom);
  plotsurfoptions.push_back(&zlimbottom);
  plotsurfoptions.push_back(&xstep);
  plotsurfoptions.push_back(&ystep);
  plotsurfoptions.push_back(&zstep);
  plotsurfoptions.push_back(&xstart);
  plotsurfoptions.push_back(&ystart);
  plotsurfoptions.push_back(&zstart);
  plotsurfoptions.push_back(&gridsize);
  plotsurfoptions.push_back(&linewidth);
  plotsurfoptions.push_back(&fontsize);
  plotsurfoptions.push_back(&pointsize);
  plotsurfoptions.push_back(&linecolor);
  plotsurfoptions.push_back(&titlescale);

  methods.push_back(command("plotsurf",&m,&plotsurfoptions,&u,
                    required,notallowed,notallowed,optional,optional,required));

  functions[4] = plotsurfrun;

  }


#if defined(JAVA_OUTPUT_WINDOW)
graphobj::graphobj(administrator_basic * adb, administrator_pointer * adp,
                   const ST::string & n,ofstream * lo,istream * in,
                   vector<statobject*> * st)
			 : statobject(adb,n,"graph",lo,in)
  {
  adminp_p=adp;
  statobj = st;
  create();
  }
#else
graphobj::graphobj(const ST::string & n,
                   ofstream * lo,istream * in,
                   vector<statobject*> * st)
			 : statobject(n,"graph",lo,in)
  {
  statobj = st;
  create();
  }
#endif


graphobj::graphobj(const graphobj & o) : statobject(statobject(o))
  {
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = o.adminp_p;
  #endif
  statobj = o.statobj;
  mapinfo = o.mapinfo;
  D = o.D;
  varnames = o.varnames;
  }



const graphobj & graphobj::operator=(const graphobj & o)
  {
  if (this == &o)
    return *this;
  statobject::operator=(statobject(o));
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = o.adminp_p;
  #endif
  statobj = o.statobj;
  mapinfo = o.mapinfo;
  D = o.D;
  varnames = o.varnames;
  return *this;
  }


int graphobj::parse(const ST::string & c)
  {
  optionlist globaloptions = optionlist();
  int pos = statobject::parsecom(c,methods,globaloptions);

  if (pos >= 0)
	 (*functions[pos])(*this);

  }


void graphobj::describe(const optionlist & globaloptions)
  {

  vector<ST::string> varnames = m.getModelVarnamesAsVector();

  if(varnames.size()>0 || mapname.getvalue()!="")
    {
    out("\n");
    out("\n");
    out("GRAPH " + name + "\n",true,false,16);

    if(varnames.size()>0)
      {
      ST::string str = "";
      for(int i=0;i<varnames.size();i++)
        str = str + "'" + varnames[i] + "' ";
      out("Variable(s) " + str + "of dataset '" + u.getusingtext() + "' used for last graph\n");
      out("\n");
      }
    if(mapname.getvalue() != "")
      {
      out("Map used for last graph: '" + mapname.getvalue() + "'\n");
      out("\n");
      }

    }
  else
    out("NOTE: graph object does not contain any data\n");

  }


void drawmaprun(graphobj & o)
  {

  bool error = false;
  bool drawonlyboundaries = false;

  // ----------------------------- reading map info ----------------------------

  mapobject * mapp;                           // pointer to mapobject

  if (error == false)
    {

    int objpos = findstatobject(*(o.statobj),o.mapname.getvalue(),"map");

    if (objpos >= 0)
      {
      statobject * s = o.statobj->at(objpos);
      mapp = dynamic_cast<mapobject*>(s);
      o.mapinfo = mapp->getmap();
      }
    else
      {
      if (objpos == -1)
        o.outerror("ERROR: map object " + o.mapname.getvalue() + " is not existing\n");
      else
        o.outerror("ERROR: " + o.mapname.getvalue() + " is not a map object\n");
      error = true;
      }

    }

  if (error == false)
    {
    if (o.mapinfo.polygones_existing() == false)
      {
      error=true;
      o.outerror("ERROR: boundaries of the regions are not available\n");
      }
    }

  // ------------------------- end: reading map info ---------------------------

  vector<ST::string> varnames = o.m.getModelVarnamesAsVector();

  if (error==false)
    {

    if(varnames.size() == 0)
      {
      if(o.u.getnotext())
        {
        drawonlyboundaries = true;
        }
      else
        {
        error = true;
        o.outerror("ERROR: 2 variable names required\n");
        }
      }
    else if (varnames.size() == 1)
      {
      error = true;
      o.outerror("ERROR: 2 variable names required\n");
      }
    else if (varnames.size() == 2)
      {
      if(o.u.getnotext())
        {
        error = true;
        o.outerror("ERROR: using required\n");
        }
      }
    else if (varnames.size() > 2)
      {
      error = true;
      o.outerror("ERROR: only 2 variables allowed\n");
      }

    }

if(drawonlyboundaries)
  {
  if(o.color.changed())
    {
    o.outerror("ERROR: option " + o.color.getname() + " not meaningful for drawing boundaries only\n");
    error = true;
    }
  if(o.swapcolors.changed())
    {
    o.outerror("ERROR: option " + o.swapcolors.getname() + " not meaningful for drawing boundaries only\n");
    error = true;
    }
  if(o.upperlimit.changed())
    {
    o.outerror("ERROR: option " + o.upperlimit.getname() + " not meaningful for drawing boundaries only\n");
    error = true;
    }
  if(o.lowerlimit.changed())
    {
    o.outerror("ERROR: option " + o.lowerlimit.getname() + " not meaningful for drawing boundaries only\n");
    error = true;
    }
  if(o.nolegend.changed())
    {
    o.outerror("ERROR: option " + o.nolegend.getname() + " not meaningful for drawing boundaries only\n");
    error = true;
    }
  if(o.nrcolors.changed())
    {
    o.outerror("ERROR: option " + o.nrcolors.getname() + " not meaningful for drawing boundaries only\n");
    error = true;
    }


  if(error==false)
    {

    ST::string path = o.psname.getvalue();
    if ( path.isvalidfile() == 1 && path != "")
        {
        o.outerror("ERROR: file " + path + " is not a valid file name\n");
        error = true;
        }
    if ( (path.isexistingfile() == 0) && (o.replace.getvalue() == false) )
        {
        o.outerror("ERROR: file " + path + " is already existing\n");
        error = true;
        }

#if defined(JAVA_OUTPUT_WINDOW)
    o.adminp_p->set_Dp(&o.D);
    o.adminp_p->set_mapinfo(&o.mapinfo);

    if(drawonlyboundaries)
      {
      jmethodID javashowmap = o.adminb_p->Java->GetMethodID(o.adminb_p->BayesX_cls, "JavaShowMap", "(ZILjava/lang/String;)V");
      o.adminb_p->Java->CallVoidMethod(o.adminb_p->BayesX_obj, javashowmap, o.drawnames.getvalue(),
                                       o.fontsize.getvalue(),
                                       o.adminb_p->Java->NewStringUTF(o.psname.getvalue().strtochar()));
      }
    else
      {
      jmethodID javadrawmap = o.adminb_p->Java->GetMethodID(o.adminb_p->BayesX_cls, "Javadrawmap",
                                                    "(ZZZZZDDSZILjava/lang/String;Ljava/lang/String;)V");
      o.adminb_p->Java->CallVoidMethod(o.adminb_p->BayesX_obj, javadrawmap, false, false, false, o.drawnames.getvalue(),
                               false, 0.0, 1.0, 0, false, o.fontsize.getvalue(),
                               o.adminb_p->Java->NewStringUTF(o.psname.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.title.getvalue().strtochar()));
      }
    bool stop=o.adminb_p->breakcommand();
#endif

/*
    admin.set_mapinfo(&o.mapinfo);

    jmethodID javashowmap = o.adminb_p->Java->GetMethodID(o.adminb_p->BayesX_cls, "JavaShowMap", "(Z)V");
    o.adminb_p->Java->CallVoidMethod(o.adminb_p->BayesX_obj, javashowmap, o.drawnames.getvalue());
*/
    }

  }
else
  {
  //--------------------- reading dataset information --------------------------

  dataobject * datap;               // pointer to dataset

  if(error == false)
    {

    int objpos = findstatobject(*(o.statobj),o.u.getusingtext(),"dataset");

    statobject * s;
    if (objpos >= 0)
      {
      s = o.statobj->at(objpos);
      datap = dynamic_cast<dataobject*>(s);
      }
    else
      {
      if (objpos == -1)
        {
        o.outerror("ERROR: dataset object " + o.u.getusingtext() + " is not existing\n");
        error = true;
        }
      else
        {
        o.outerror("ERROR: " + o.u.getusingtext() + " is not a dataset object\n");
        error = true;
        }
      }
    }

  //------------------ end: reading dataset information ------------------------

  if (error == false)
    {

    vector<ST::string> notex;

    if (datap->allexisting(varnames,notex) == false)
      {
      error = true;
      unsigned i;
      for (i=0;i<notex.size();i++)
        o.outerror("ERROR: variable " + notex[i] + " is not existing\n");
      }

    }

  if (error == false)
    {

    ST::string ifexpression = o.methods[0].getexpression();

    datap->makematrix(varnames,o.D,ifexpression);

    vector<ST::string> em = datap->geterrormessages();

    if (em.size() > 0)
      {
      o.outerror(em);
      error = true;
      }

    }

  double upperlim;
  double lowerlim;

  if (error == false)
    {
/*
    if(o.lowerlimit.getvalue()!=-1.0E7)
        lowerlim = o.lowerlimit.getvalue();
    else
        lowerlim = o.D.min(0);
    if(o.upperlimit.getvalue()!=1.0E7)
        upperlim = o.upperlimit.getvalue();
    else
        upperlim = o.D.max(0);
*/
    if(o.lowerlimit.changed())
        lowerlim = o.lowerlimit.getvalue();
    else
        lowerlim = o.D.min(0);

    if(o.upperlimit.changed())
        upperlim = o.upperlimit.getvalue();
    else
        upperlim = o.D.max(0);

    if(upperlim < lowerlim)
        {
        o.outerror("ERROR: 'lowerlimit' is above 'upperlimit'\n");
        error = true;
        }

    if(o.D.min(0) == o.D.max(0))
        {
        o.out("NOTE: variable " + varnames[0] + " does not vary\n");
        o.nrcolors.setvalue(1);
        }

    ST::string path = o.psname.getvalue();
    if ( path.isvalidfile() == 1 && path != "")
        {
        o.outerror("ERROR: file " + path + " is not a valid file name\n");
        error = true;
        }
    if ( (path.isexistingfile() == 0) && (o.replace.getvalue() == false) )
        {
        o.outerror("ERROR: file " + path + " is already existing\n");
        error = true;
        }
    }


  if (error == false)
    {

    // Landkarte zeichnen einfügen
    // Name der plotvariable:   plotvar
    // Name der Regionvariable: regionvar
    // Daten: erste Spalte von D enthält plotvar, zweite Spalte von D enthält
    //        regionvar
    // Map: m
    // Farbe ja/nein: o.color.getvalue() == true bedeutet Farbe, ansonsten
    //                schwarzweiss

//    if(o.lowerlimit.getvalue() > o.D.min(0))
    if(lowerlim > o.D.min(0))
        {
        o.out("NOTE: 'lowerlimit' is above minimum value (minimum value is " + ST::doubletostring(o.D.min(0)) + ")\n");
        }
//    if(o.upperlimit.getvalue() < o.D.max(0))
    if(upperlim < o.D.max(0))
        {
        o.out("NOTE: 'upperlimit' is below maximum value (maximum value is " + ST::doubletostring(o.D.max(0)) + ")\n");
        }

#if defined(JAVA_OUTPUT_WINDOW)
    o.adminp_p->set_Dp(&o.D);
    o.adminp_p->set_mapinfo(&o.mapinfo);

    bool legend;
    if(o.nolegend.getvalue())
      legend = false;
    else
      legend = true;

    jmethodID javadrawmap = o.adminb_p->Java->GetMethodID(o.adminb_p->BayesX_cls, "Javadrawmap",
                                                    "(ZZZZZDDSZILjava/lang/String;Ljava/lang/String;D)V");
    o.adminb_p->Java->CallVoidMethod(o.adminb_p->BayesX_obj, javadrawmap, o.color.getvalue(), legend,
                               o.swapcolors.getvalue(), o.drawnames.getvalue(), o.hclcolors.getvalue(),
                               lowerlim, upperlim, o.nrcolors.getvalue(), o.pcat.getvalue(),
                               o.fontsize.getvalue(),
                               o.adminb_p->Java->NewStringUTF(o.psname.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.title.getvalue().strtochar()),
                               o.titlescale.getvalue());

    bool stop=o.adminb_p->breakcommand();
#endif
    }

  } // END: if(drawonlyboundaries)

  o.u = use();

  }


void plotnonprun(graphobj & o)
  {

  bool error = false;

  //--------------------- reading dataset information --------------------------

  dataobject * datap;               // pointer to dataset

  int objpos = findstatobject(*(o.statobj),o.u.getusingtext(),"dataset");

  statobject * s;
  if (objpos >= 0)
    {
    s = o.statobj->at(objpos);
    datap = dynamic_cast<dataobject*>(s);
    }
  else
    {
    if (objpos == -1)
      {
      o.outerror("ERROR: " + o.u.getusingtext() + " is not existing\n");
      error = true;
      }
    else
      {
      o.outerror("ERROR: " + o.u.getusingtext() + " is not a dataset object\n");
      error = true;
      }

    }

  //------------------ end: reading dataset information ------------------------

  vector<ST::string> varnames;

  if (error==false)
    {

    varnames = o.m.getModelVarnamesAsVector();

    if (varnames.size() < 2)
      {
      error = true;
      o.outerror("ERROR: at least two variable names expected\n");
      }

    }


  if (error == false)
    {

    vector<ST::string> notex;

    if (datap->allexisting(varnames,notex) == false)
      {
      error = true;
      unsigned i;
      for (i=0;i<notex.size();i++)
        o.outerror("ERROR: variable " + notex[i] + " is not existing\n");
      }

    }

  if (error == false)
    {

    ST::string ifexpression = o.methods[0].getexpression();

    datap->makematrix(varnames,o.D,ifexpression);

    vector<ST::string> em = datap->geterrormessages();

    if (em.size() > 0)
      {
      o.outerror(em);
      error = true;
      }

    }

  if (error == false)
    {
    ST::string path = o.psname.getvalue();
    if ( path.isvalidfile() == 1 && path != "" )
        {
        o.outerror("ERROR: file " + path + " is not a valid file name\n");
        error = true;
        }
    if ( (path.isexistingfile() == 0) && (o.replace.getvalue() == false) )
        {
        o.outerror("ERROR: file " + path + " is already existing\n");
        error = true;
        }
    }

  if(error == false)
    {

    if( o.xlimtop.getvalue() > -MAXDOUBLE && o.xlimbottom.getvalue() < MAXDOUBLE
        && o.xlimtop.getvalue() <= o.xlimbottom.getvalue() )
      {
      o.outerror("ERROR: 'xlimtop' must be bigger than 'xlimbottom'\n");
      error = true;
      }

    if(o.ylimtop.getvalue() > -MAXDOUBLE && o.ylimbottom.getvalue() < MAXDOUBLE
        && o.ylimtop.getvalue() <= o.ylimbottom.getvalue() )
      {
      o.outerror("ERROR: 'ylimtop' must be bigger than 'ylimbottom'\n");
      error = true;
      }

    if(!o.xlimtop.changed())
      o.xlimtop.setvalue(-MAXDOUBLE);
    if(!o.xlimbottom.changed())
      o.xlimbottom.setvalue(MAXDOUBLE);

    if(!o.ylimtop.changed())
      o.ylimtop.setvalue(-MAXDOUBLE);
    if(!o.ylimbottom.changed())
      o.ylimbottom.setvalue(MAXDOUBLE);

    if(o.year.getvalue() != 0 || o.month.getvalue() != 0)
      {
      if(o.year.getvalue() == 0)
        {
        o.outerror("ERROR: option 'year' expected\n");
        error = true;
        }
      if(o.month.getvalue() == 0)
        {
        o.month.setvalue(1);
        }
      }

    }

  if (error == false)
    {

#if defined(JAVA_OUTPUT_WINDOW)
    o.adminp_p->set_Dp(&o.D);

    jmethodID javaplotnonp = o.adminb_p->Java->GetMethodID(o.adminb_p->BayesX_cls, "Javaplotnonp",
                "(Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;IIDDDDDDDDIIIIID)V");
    o.adminb_p->Java->CallVoidMethod(o.adminb_p->BayesX_obj, javaplotnonp,
                               o.adminb_p->Java->NewStringUTF(o.psname.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.title.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.xlab.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.ylab.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.connect.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.linecolor.getvalue().strtochar()),
                               o.height.getvalue(),o.width.getvalue(),
                               o.xlimtop.getvalue(),o.xlimbottom.getvalue(),
                               o.ylimtop.getvalue(),o.ylimbottom.getvalue(),
                               o.xstep.getvalue(),o.xstart.getvalue(),
                               o.ystep.getvalue(),o.ystart.getvalue(),
                               o.year.getvalue(),o.month.getvalue(),
                               o.linewidth.getvalue(),o.pointsize.getvalue(),
                               o.fontsize.getvalue(),o.titlescale.getvalue()
                               );

    bool stop=o.adminb_p->breakcommand();
#endif
    }

  }


void plotsamplerun(graphobj & o)
  {

  bool error = false;

  //--------------------- reading dataset information --------------------------

  dataobject * datap;               // pointer to dataset

  int objpos = findstatobject(*(o.statobj),o.u.getusingtext(),"dataset");

  statobject * s;
  if (objpos >= 0)
    {
    s = o.statobj->at(objpos);
    datap = dynamic_cast<dataobject*>(s);
    }
  else
    {
    if (objpos == -1)
      {
      o.outerror("ERROR: " + o.u.getusingtext() + " is not existing\n");
      error = true;
      }
    else
      {
      o.outerror("ERROR: " + o.u.getusingtext() + " is not a dataset object\n");
      error = true;
      }

    }

  //------------------ end: reading dataset information ------------------------

  list<ST::string> varnames;

  if (error == false)
    {

    varnames = datap->getVarnames();

    ST::string ifexpression = o.methods[0].getexpression();

    datap->makematrix(varnames,o.D,ifexpression);

    vector<ST::string> em = datap->geterrormessages();

    if (em.size() > 0)
      {
      o.outerror(em);
      error = true;
      }

    }

  if (error == false)
    {
    ST::string path = o.psname.getvalue();
    if ( path.isvalidfile() == 1 && path != "" )
        {
        o.outerror("ERROR: file " + path + " is not a valid file name\n");
        error = true;
        }
    if ( (path.isexistingfile() == 0) && (o.replace.getvalue() == false) )
        {
        o.outerror("ERROR: file " + path + " is already existing\n");
        error = true;
        }
    }


  if (error == false)
    {

#if defined(JAVA_OUTPUT_WINDOW)
    o.adminp_p->set_Dp(&o.D);

    jmethodID javaplotsample = o.adminb_p->Java->GetMethodID(o.adminb_p->BayesX_cls, "Javaplotsample",
                "(Ljava/lang/String;Ljava/lang/String;)V");
    o.adminb_p->Java->CallVoidMethod(o.adminb_p->BayesX_obj, javaplotsample,
                               o.adminb_p->Java->NewStringUTF(o.psname.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.connect.getvalue().strtochar())
                               );

    bool stop=o.adminb_p->breakcommand();
#endif
    }

  }


void plotautocorrun(graphobj & o)
  {

  bool error = false;

  //--------------------- reading dataset information --------------------------

  dataobject * datap;               // pointer to dataset

  int objpos = findstatobject(*(o.statobj),o.u.getusingtext(),"dataset");

  statobject * s;
  if (objpos >= 0)
    {
    s = o.statobj->at(objpos);
    datap = dynamic_cast<dataobject*>(s);
    }
  else
    {
    if (objpos == -1)
      {
      o.outerror("ERROR: " + o.u.getusingtext() + " is not existing\n");
      error = true;
      }
    else
      {
      o.outerror("ERROR: " + o.u.getusingtext() + " is not a dataset object\n");
      error = true;
      }

    }

  //------------------ end: reading dataset information ------------------------

//  list<ST::string> varnames;

  if (error == false)
    {

    list<ST::string> helpvarnames = datap->getVarnames();
    list<ST::string>::iterator it = helpvarnames.begin();

    o.varnames = vector<ST::string>();

    if(o.mean.getvalue() && helpvarnames.size()>0)
      {
      o.varnames.push_back(*it);
//      varnames.push_back(*it);
      it++;
      while(it != helpvarnames.end())
        {
        if((*it).endswith("mean") || (*it).endswith("min") || (*it).endswith("max"))
          {
          o.varnames.push_back(*it);
//          varnames.push_back(*it);
          }
        it++;
        }
      }
    else
      {
      while(it != helpvarnames.end())
        {
        if(!(*it).endswith("mean") && !(*it).endswith("min") && !(*it).endswith("max"))
          {
          o.varnames.push_back(*it);
//          varnames.push_back(*it);
          }
        it++;
        }
      }

    ST::string ifexpression = o.methods[0].getexpression();

    datap->makematrix(o.varnames,o.D,ifexpression);

    vector<ST::string> em = datap->geterrormessages();

    if (em.size() > 0)
      {
      o.outerror(em);
      error = true;
      }

    }

  if (error == false)
    {
    ST::string path = o.psname.getvalue();
    if ( path.isvalidfile() == 1 && path != "" )
        {
        o.outerror("ERROR: file " + path + " is not a valid file name\n");
        error = true;
        }
    if ( (path.isexistingfile() == 0) && (o.replace.getvalue() == false) )
        {
        o.outerror("ERROR: file " + path + " is already existing\n");
        error = true;
        }
    }


  if (error == false)
    {

#if defined(JAVA_OUTPUT_WINDOW)
    o.adminp_p->set_Dp(&o.D);
    o.adminp_p->set_varnamesp(&o.varnames);

    jmethodID javaplotautocor = o.adminb_p->Java->GetMethodID(o.adminb_p->BayesX_cls, "Javaplotautocor",
                "(Ljava/lang/String;Ljava/lang/String;Z)V");
    o.adminb_p->Java->CallVoidMethod(o.adminb_p->BayesX_obj, javaplotautocor,
                               o.adminb_p->Java->NewStringUTF(o.psname.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.connect.getvalue().strtochar()),
                               o.mean.getvalue()
                               );

    bool stop=o.adminb_p->breakcommand();
#endif
    }

  }


void plotsurfrun(graphobj & o)
  {

  bool error = false;

  //--------------------- reading dataset information --------------------------

  dataobject * datap;               // pointer to dataset

  int objpos = findstatobject(*(o.statobj),o.u.getusingtext(),"dataset");

  statobject * s;
  if (objpos >= 0)
    {
    s = o.statobj->at(objpos);
    datap = dynamic_cast<dataobject*>(s);
    }
  else
    {
    if (objpos == -1)
      {
      o.outerror("ERROR: " + o.u.getusingtext() + " is not existing\n");
      error = true;
      }
    else
      {
      o.outerror("ERROR: " + o.u.getusingtext() + " is not a dataset object\n");
      error = true;
      }

    }

  //------------------ end: reading dataset information ------------------------

  vector<ST::string> varnames;

  if (error==false)
    {

    varnames = o.m.getModelVarnamesAsVector();

    if (varnames.size() != 3)
      {
      error = true;
      o.outerror("ERROR: 3 variable names expected\n");
      }

    }


  if (error == false)
    {

    vector<ST::string> notex;

    if (datap->allexisting(varnames,notex) == false)
      {
      error = true;
      unsigned i;
      for (i=0;i<notex.size();i++)
        o.outerror("ERROR: variable " + notex[i] + " is not existing\n");
      }

    }

  if (error == false)
    {

    ST::string ifexpression = o.methods[0].getexpression();

    datap->makematrix(varnames,o.D,ifexpression);

    vector<ST::string> em = datap->geterrormessages();

    if (em.size() > 0)
      {
      o.outerror(em);
      error = true;
      }

    }

  if (error == false)
    {
    ST::string path = o.psname.getvalue();
    if ( path.isvalidfile() == 1 && path != "" )
        {
        o.outerror("ERROR: file " + path + " is not a valid file name\n");
        error = true;
        }
    if ( (path.isexistingfile() == 0) && (o.replace.getvalue() == false) )
        {
        o.outerror("ERROR: file " + path + " is already existing\n");
        error = true;
        }
    }

  if (error == false)
    {

#if defined(JAVA_OUTPUT_WINDOW)
    o.adminp_p->set_Dp(&o.D);

    jmethodID javaplotsurf = o.adminb_p->Java->GetMethodID(o.adminb_p->BayesX_cls, "Javaplotsurf",
                "(Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;DDDLjava/lang/String;IIDDDDDDDDDDDDIIIID)V");
    o.adminb_p->Java->CallVoidMethod(o.adminb_p->BayesX_obj, javaplotsurf,
                               o.adminb_p->Java->NewStringUTF(o.psname.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.title.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.xlab.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.ylab.getvalue().strtochar()),
                               o.adminb_p->Java->NewStringUTF(o.zlab.getvalue().strtochar()),
                               o.xrot.getvalue(),o.yrot.getvalue(),o.zrot.getvalue(),
                               o.adminb_p->Java->NewStringUTF(o.linecolor.getvalue().strtochar()),
                               o.height.getvalue(),o.width.getvalue(),
                               o.xlimtop.getvalue(),o.xlimbottom.getvalue(),
                               o.ylimtop.getvalue(),o.ylimbottom.getvalue(),
                               o.zlimtop.getvalue(),o.zlimbottom.getvalue(),
                               o.xstep.getvalue(),o.xstart.getvalue(),
                               o.ystep.getvalue(),o.ystart.getvalue(),
                               o.zstep.getvalue(),o.zstart.getvalue(),
                               o.gridsize.getvalue(),
                               o.linewidth.getvalue(),o.pointsize.getvalue(),
                               o.fontsize.getvalue(),o.titlescale.getvalue()
                               );

    bool stop=o.adminb_p->breakcommand();
#endif
    }

  }



void graphobj::changedescription(void)
  {



  }





