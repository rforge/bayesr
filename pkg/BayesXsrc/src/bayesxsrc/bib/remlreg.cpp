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

#include<StatwinFrame.h>
#include<statwin_haupt.h>

#endif

#include"remlreg.h"
#include"fullcond.h"
//#include<typeinfo.h>
#include<stddef.h>


//------------------------------------------------------------------------------
//------------- CLASS remlreg: implementation of member functions -------------
//------------------------------------------------------------------------------

void remlreg::make_paths(unsigned  collinpred,ST::string & pathnonp,
                          ST::string & pathres,ST::string & title,
                          ST::string varname1,ST::string varname2,
                          ST::string endingraw,ST::string endingres,
                          ST::string endingtitle)
  {
  if (collinpred == 0)
    {
    if (varname2=="")
      {
#if defined(__BUILDING_LINUX)
      pathnonp = defaultpath + "/temp/" + name +  add_name + "_f_" +
                 varname1 + endingraw;
#else
      pathnonp = defaultpath + "\\temp\\" + name +  add_name + "_f_" +
                 varname1 + endingraw;
#endif

      pathres = outfile.getvalue() + add_name + "_f_" + varname1 + endingres;

      title = "f_"+ varname1 + endingtitle +  add_name;
      }
    else
      {
#if defined(__BUILDING_LINUX)
      pathnonp = defaultpath + "/temp/" + name +  add_name + "_" +
                 varname2 +  "_f_" + varname1 + endingraw;
#else
      pathnonp = defaultpath + "\\temp\\" + name +  add_name + "_" +
                 varname2 +  "_f_" + varname1 + endingraw;
#endif
      pathres = outfile.getvalue() + add_name + "_" + varname2 +
                "_f_" + varname1 + endingres;

      title = varname2 + "_f_"+ varname1 + endingtitle +  add_name;
      }
    }
  else
    {
    if (varname2=="")
      {
#if defined(__BUILDING_LINUX)
      pathnonp = defaultpath + "/temp/" + name +  add_name + "_f_" +
                 ST::inttostring(collinpred+1) + "_" +
                     varname1 + endingraw;
#else
      pathnonp = defaultpath + "\\temp\\" + name +  add_name + "_f_" +
                 ST::inttostring(collinpred+1) + "_" +
                     varname1 + endingraw;
#endif
      pathres = outfile.getvalue() + add_name + "_f_" +
                ST::inttostring(collinpred+1) + "_" +
                varname1 + endingres;

      title = "f_"+ ST::inttostring(collinpred+1) + "_" +
               varname1 + endingtitle + add_name;
      }
    else
      {
#if defined(__BUILDING_LINUX)
      pathnonp = defaultpath + "/temp/" + name +  add_name + "_" + varname2
                 + "_f_" + ST::inttostring(collinpred+1) + "_" +
                     varname1 + endingraw;
#else
      pathnonp = defaultpath + "\\temp\\" + name +  add_name + "_" + varname2
                 + "_f_" + ST::inttostring(collinpred+1) + "_" +
                     varname1 + endingraw;
#endif
      pathres = outfile.getvalue() + add_name + "_" + varname2 + "_f_" +
                ST::inttostring(collinpred+1) + "_" +
                varname1 + endingres;

      title = varname2 + "_f_"+ ST::inttostring(collinpred+1) + "_" +
               varname1 + endingtitle + add_name;
      }
    }
  }

void remlreg::initpointers(void)
  {
  unsigned i;

  for(i=0;i<fcpsplinesurf.size();i++)
    fullcond.push_back(&fcpsplinesurf[i]);
  for(i=0;i<fcconst.size();i++)
    fullcond.push_back(&fcconst[i]);
  for(i=0;i<fcnonpgaussian.size();i++)
    fullcond.push_back(&fcnonpgaussian[i]);
  for(i=0;i<fcpspline.size();i++)
    fullcond.push_back(&fcpspline[i]);
  for(i=0;i<fcrandom.size();i++)
    fullcond.push_back(&fcrandom[i]);
  for(i=0; i<fckriging.size();i++)
    fullcond.push_back(&fckriging[i]);
  for(i=0; i<fcbaseline.size(); i++)
    fullcond.push_back(&fcbaseline[i]);
  for(i=0; i<fcbaseline_varcoeff.size(); i++)
    fullcond.push_back(&fcbaseline_varcoeff[i]);
  }

void remlreg::clear(void)
  {
  fullcond.erase(fullcond.begin(),fullcond.end());
  fullcond.reserve(50);
  fcconst.erase(fcconst.begin(),fcconst.end());
  fcconst.reserve(20);
  fcnonpgaussian.erase(fcnonpgaussian.begin(),fcnonpgaussian.end());
  fcnonpgaussian.reserve(20);
  fcpspline.erase(fcpspline.begin(),fcpspline.end());
  fcpspline.reserve(20);
  fcpsplinesurf.erase(fcpsplinesurf.begin(),fcpsplinesurf.end());
  fcpsplinesurf.reserve(20);
  fcrandom.erase(fcrandom.begin(),fcrandom.end());
  fcrandom.reserve(20);
  fckriging.erase(fckriging.begin(),fckriging.end());
  fckriging.reserve(20);
  fcbaseline.erase(fcbaseline.begin(),fcbaseline.end());
  fcbaseline.reserve(20);
  fcbaseline_varcoeff.erase(fcbaseline_varcoeff.begin(),fcbaseline_varcoeff.end());
  fcbaseline_varcoeff.reserve(20);
  }

void remlreg::create(void)
  {
  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  add_name="";
#if defined(__BUILDING_LINUX)
  ST::string h = defaultpath+"/output/"+name;
#else
  ST::string h = defaultpath+"\\output\\"+name;
#endif

  outfile = fileoption("outfile",h,false);

  globaloptions.push_back(&outfile);

//------------------------------------------------------------------------------
// ----------------------------- method regress --------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------- for the offset -----------------------------------
//------------------------------------------------------------------------------

  offset = term_offset();
  termtypes.push_back(&offset);

//------------------------------------------------------------------------------
//-------------------------  for fixed effects ---------------------------------
//------------------------------------------------------------------------------

  fixedeffects = basic_termtype();
  termtypes.push_back(&fixedeffects);

  fixed_catsp = term_fixed_catspecific();
  termtypes.push_back(&fixed_catsp);

//------------------------------------------------------------------------------
//------------------------ for nonparametric terms -----------------------------
//------------------------------------------------------------------------------

  nonprw1rw2 = term_autoreg_remlreg();
  termtypes.push_back(&nonprw1rw2);

  nonprw1rw2_varcoef = term_autoreg_varcoef_remlreg();
  termtypes.push_back(&nonprw1rw2_varcoef);

  nonpseason = term_season_remlreg();
  termtypes.push_back(&nonpseason);

  nonpseason_varcoef = term_season_varcoef_remlreg();
  termtypes.push_back(&nonpseason_varcoef);

  nonpspatial = term_spatial_remlreg();
  termtypes.push_back(&nonpspatial);

  nonpspatial_varcoef = term_spatial_varcoef_remlreg();
  termtypes.push_back(&nonpspatial_varcoef);

//------------------------------------------------------------------------------
//-------------------------- for p-spline terms --------------------------------
//------------------------------------------------------------------------------

  nonppspline = term_pspline_remlreg();
  termtypes.push_back(&nonppspline);

  nonpvarcoeffpspline = term_varcoeff_pspline_remlreg();
  termtypes.push_back(&nonpvarcoeffpspline);

  nonpinteractpspline = term_interactpspline_remlreg();
  termtypes.push_back(&nonpinteractpspline);

  nonpvarcoeffinteractpspline = term_interactpspline_varcoeff_remlreg();
  termtypes.push_back(&nonpvarcoeffinteractpspline);

  nonpgeospline = term_geospline_remlreg();
  termtypes.push_back(&nonpgeospline);

  nonpvarcoeffgeospline = term_geospline_varcoeff_remlreg();
  termtypes.push_back(&nonpvarcoeffgeospline);

//------------------------------------------------------------------------------
//-------------------------- for kriging terms ---------------------------------
//------------------------------------------------------------------------------

  nonpspatial_kriging = term_kriging_remlreg();
  termtypes.push_back(&nonpspatial_kriging);

  nonp_kriging = term_kriging_1dim_remlreg();
  termtypes.push_back(&nonp_kriging);

  nonpspatial_kriging_varcoeff = term_kriging_varcoeff_remlreg();
  termtypes.push_back(&nonpspatial_kriging_varcoeff);

  nonpspatial_geokriging = term_geokriging_remlreg();
  termtypes.push_back(&nonpspatial_geokriging);

  nonpspatial_geokriging_varcoeff = term_geokriging_varcoeff_remlreg();
  termtypes.push_back(&nonpspatial_geokriging_varcoeff);

//------------------------------------------------------------------------------
//-------------------------- for baseline terms --------------------------------
//------------------------------------------------------------------------------

  nonp_baseline = term_baseline_remlreg();
  termtypes.push_back(&nonp_baseline);

  nonp_baseline_varcoeff = term_baseline_varcoeff_remlreg();
  termtypes.push_back(&nonp_baseline_varcoeff);

//------------------------------------------------------------------------------
//------------------------ for random effect terms -----------------------------
//------------------------------------------------------------------------------

  randomeff = term_random_remlreg();
  termtypes.push_back(&randomeff);

  randomeffslope = term_randomslope_remlreg();
  termtypes.push_back(&randomeffslope);

//------------------------------------------------------------------------------
//--------------------------- General options ----------------------------------
//------------------------------------------------------------------------------

  modreg = modelterm(&termtypes);

  udata = use();

  knotsdef.push_back("equidistant");
  knotsdef.push_back("quantiles");
  knots = stroption("knots",knotsdef,"equidistant");

  level1 = doubleoption("level1",95,40,99);
  level2 = doubleoption("level2",80,40,99);
  maxint = intoption("maxint",150,0,20000);

  families.reserve(25);
  families.push_back("gaussian");
  families.push_back("binomial");
  families.push_back("binomialprobit");
  families.push_back("binomialcomploglog");
  families.push_back("poisson");
  families.push_back("gamma");
  families.push_back("poissondispers");
  families.push_back("binomialdispers");
  families.push_back("binomialprobitdispers");
  families.push_back("multinomial");
  families.push_back("multinomialcatsp");
  families.push_back("cumlogit");
  families.push_back("cumprobit");
  families.push_back("seqlogit");
  families.push_back("seqprobit");
  families.push_back("cox");
  families.push_back("coxinterval");
  families.push_back("coxold");
  families.push_back("aft");
  families.push_back("multistate");
  family = stroption("family",families,"binomial");

  maxit = intoption("maxit",400,1,100000);
  lowerlim = doubleoption("lowerlim",0.001,0,1);
  eps = doubleoption("eps",0.00001,0,1);
  maxchange = doubleoption("maxchange",1000000,0,100000000);
  maxvar = doubleoption("maxvar",100000,0,100000000);

  reference = doubleoption("reference",0,-10000,10000);

  noconst = simpleoption("noconst",false);

  aiccontrol = simpleoption("aiccontrol",false);
  fisher = simpleoption("fisher",false);

  constlambda = simpleoption("constlambda",false);
  constscale = simpleoption("constscale",false);

  leftint = stroption("leftint");
  lefttrunc = stroption("lefttrunc");
  state = stroption("state");
  binomweight = stroption("binomweight");
  naindicator = stroption("naindicator");
  globalfrailty = stroption("globalfrailty");
  gflambdastart = doubleoption("gflambdastart",1000,0,10000000);

  regressoptions.reserve(100);

  regressoptions.push_back(&level1);
  regressoptions.push_back(&level2);
  regressoptions.push_back(&maxint);
  regressoptions.push_back(&family);

  regressoptions.push_back(&knots);
  regressoptions.push_back(&maxit);
  regressoptions.push_back(&lowerlim);
  regressoptions.push_back(&eps);
  regressoptions.push_back(&maxchange);
  regressoptions.push_back(&maxvar);

  regressoptions.push_back(&reference);

  regressoptions.push_back(&noconst);
  regressoptions.push_back(&fisher);

  regressoptions.push_back(&aiccontrol);

  regressoptions.push_back(&constlambda);
  regressoptions.push_back(&constscale);

  regressoptions.push_back(&leftint);
  regressoptions.push_back(&lefttrunc);
  regressoptions.push_back(&state);
  regressoptions.push_back(&binomweight);
  regressoptions.push_back(&naindicator);
  regressoptions.push_back(&globalfrailty);
  regressoptions.push_back(&gflambdastart);

//------------------------------------------------------------------------------
//-------------------------- methods[0]: remlrun -------------------------------
//------------------------------------------------------------------------------

  methods.push_back(command("regress",&modreg,&regressoptions,&udata,required,
			 optional,optional,optional,optional,required));

  functions[0] = remlrun;

//------------------------------------------------------------------------------
// --------------------------- methods[1] plotnonp -----------------------------
//------------------------------------------------------------------------------

  uplotnonp = use();
  mplotnonp = modelStandard();

  xlab = stroption("xlab");
  ylab = stroption("ylab");
  connect = stroption("connect");
  height = intoption("height",210,0,500);
  width = intoption("width",356,0,500);
  ylimtop = doubleoption("ylimtop",-MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  ylimbottom = doubleoption("ylimbottom",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xlimtop = doubleoption("xlimtop",-MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xlimbottom = doubleoption("xlimbottom",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  xstep = doubleoption("xstep",0.0,-MAXDOUBLE,MAXDOUBLE);
  ystep = doubleoption("ystep",0.0,-MAXDOUBLE,MAXDOUBLE);
  xstart = doubleoption("xstart",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  ystart = doubleoption("ystart",MAXDOUBLE,-MAXDOUBLE,MAXDOUBLE);
  linewidth = intoption("linewidth",5,0,100);
  fontsize = intoption("fontsize",12,0,100);
  pointsize = intoption("pointsize",20,0,100);
  linecolor = stroption("linecolor");
  titlescale = doubleoption("titlesize",1.5,0.0,MAXDOUBLE);

  vector<ST::string> levelchoice;
  levelchoice.reserve(4);
  levelchoice.push_back("all");
  levelchoice.push_back("1");
  levelchoice.push_back("2");
  levelchoice.push_back("none");
  outfile2=stroption("outfile");
  title = stroption("title");
  replace2 = simpleoption("replace",false);

  levels = stroption("levels",levelchoice,"all");
  median = simpleoption("median",false);

  plotnonpoptions.push_back(&xlab);
  plotnonpoptions.push_back(&ylab);
  plotnonpoptions.push_back(&connect);
  plotnonpoptions.push_back(&height);
  plotnonpoptions.push_back(&width);
  plotnonpoptions.push_back(&ylimtop);
  plotnonpoptions.push_back(&ylimbottom);
  plotnonpoptions.push_back(&ystep);
  plotnonpoptions.push_back(&ystart);
  plotnonpoptions.push_back(&xlimtop);
  plotnonpoptions.push_back(&xlimbottom);
  plotnonpoptions.push_back(&xstep);
  plotnonpoptions.push_back(&xstart);
  plotnonpoptions.push_back(&levels);
  plotnonpoptions.push_back(&median);
  plotnonpoptions.push_back(&title);
  plotnonpoptions.push_back(&outfile2);
  plotnonpoptions.push_back(&replace2);
  plotnonpoptions.push_back(&linewidth);
  plotnonpoptions.push_back(&fontsize);
  plotnonpoptions.push_back(&pointsize);
  plotnonpoptions.push_back(&linecolor);
  plotnonpoptions.push_back(&titlescale);

  // methods[1]:

  methods.push_back(command("plotnonp",&mplotnonp,&plotnonpoptions,&uplotnonp,
                   required,notallowed,notallowed,notallowed,optional,
                   notallowed));

  functions[1] = plotnonprun;

//------------------------------------------------------------------------------
// -------------------------- methods[2] drawmaprun ----------------------------
//------------------------------------------------------------------------------

  udrawmap = use();

  mdrawmap = modelStandard();

  outfile4 = stroption("outfile");
  title2 = stroption("title");
  upperlimit = doubleoption("upperlimit",1,-MAXDOUBLE,MAXDOUBLE);
  lowerlimit = doubleoption("lowerlimit",0,-MAXDOUBLE,MAXDOUBLE);
  nrcolors = intoption("nrcolors",256,1,256);
  color = simpleoption("color",false);
  nolegend = simpleoption("nolegend",false);
  swapcolors = simpleoption("swapcolors",false);
  replace = simpleoption("replace",false);
  plotvar = stroption("plotvar","pmode");
  pcat = simpleoption("pcat",false);
  drawnames = simpleoption("drawnames",false);
  hclcolors = simpleoption("hcl",false);

  drawmapoptions.push_back(&outfile4);
  drawmapoptions.push_back(&title2);
  drawmapoptions.push_back(&upperlimit);
  drawmapoptions.push_back(&lowerlimit);
  drawmapoptions.push_back(&nrcolors);
  drawmapoptions.push_back(&color);
  drawmapoptions.push_back(&nolegend);
  drawmapoptions.push_back(&swapcolors);
  drawmapoptions.push_back(&replace);
  drawmapoptions.push_back(&plotvar);
  drawmapoptions.push_back(&pcat);
  drawmapoptions.push_back(&drawnames);
  drawmapoptions.push_back(&hclcolors);
  drawmapoptions.push_back(&fontsize);
  drawmapoptions.push_back(&titlescale);

  // methods[2]:
  methods.push_back(command("drawmap",&mdrawmap,&drawmapoptions,&udrawmap,
                   required,notallowed,notallowed,notallowed,optional,
                   notallowed));

  functions[2] = drawmaprun;
//------------------------------------------------------------------------------
// ------------------------ methods[3] texsummary ------------------------------
//------------------------------------------------------------------------------

  utexsummary = use();

  mtexsummary = modelStandard();

  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  // methods 2
  methods.push_back(command("texsummary",&mtexsummary,&texsummaryoptions,&utexsummary,
                   notallowed,notallowed,notallowed,notallowed,notallowed,
                   notallowed));

  functions[3] = texsummaryrun;

//------------------------------------------------------------------------------
// -------------------------- methods[4] mremlrun ------------------------------
//------------------------------------------------------------------------------

  modregmult = modeltermmult(&termtypes);
  methods.push_back(command("mregress",&modregmult,&regressoptions,&udata,
                   required,optional,optional,optional,optional,required));
//  methods.push_back(command("mregress",&modreg,&regressoptions,&udata,required,
//			 optional,optional,optional,optional,required));

  functions[4] = mremlrun;

  }

//------------------------------------------------------------------------------
// ------------------------------- Constructor ---------------------------------
//------------------------------------------------------------------------------

  // CONSTRUCTOR 1:
  // ADDITIONAL INFORMATION: name = n

remlreg::remlreg(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adb,
  #endif
  const ST::string & n,ofstream * lo,istream * in,
						 ST::string p,vector<statobject*> * st)
						 : statobject(
                         #if defined(JAVA_OUTPUT_WINDOW)
                         adb,
                         #endif
                         n,"remlreg",lo,in,p)
  {
  statobj = st;
  create();
  resultsyesno = false;
  describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
  }

  // COPY CONSTRUCTOR

remlreg::remlreg(const remlreg & b) : statobject(statobject(b))
  {
  create();
  statobj = b.statobj;
  D = b.D;
  modelvarnamesv = b.modelvarnamesv;
  terms = b.terms;
  fcpspline = b.fcpspline;
  resultsyesno = b.resultsyesno;
  leftintpos = b.leftintpos;
  lefttruncpos = b.lefttruncpos;
  statepos = b.statepos;
  binomweightpos = b.binomweightpos;
  initpointers();
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const remlreg & remlreg::operator=(const remlreg & b)
  {
  if (this == & b)
	 return *this;
  statobject::operator=(statobject(b));
  create();
  statobj = b.statobj;
  D = b.D;
  modelvarnamesv = b.modelvarnamesv;
  terms = b.terms;
  fcpspline = b.fcpspline;
  resultsyesno = b.resultsyesno;
  statepos = b.statepos;
  binomweightpos = b.binomweightpos;
  initpointers();
  return *this;
  }

int remlreg::parse(const ST::string & c)
  {
  int u = statobject::parse(c);
  int pos = statobject::parsecom(c,methods,globaloptions);

  if (pos >= 0)
	 (*functions[pos])(*this);

  return(pos);
  }

void remlreg::describe(const optionlist & globaloptions)
  {
  statobject::describe(globaloptions);
  }

// -----------------------------------------------------------------------------
// ---------------------------- Create data ------------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_data(datamatrix & weight)
  {
  unsigned i,j;
  bool failure=false;

  // find the dataset object

  dataobject * datap;               // pointer to dataset
  int objpos = findstatobject(*statobj,udata.getusingtext(),"dataset");
  statobject * s;
  if (objpos >= 0)
    {
    s = statobj->at(objpos);
    datap = dynamic_cast<dataobject*>(s);
    }
  else
    {
    if (objpos == -1)
      outerror("ERROR: " + udata.getusingtext() + " is not existing\n");
    else
      outerror("ERROR: " + udata.getusingtext() + " is not a dataset object\n");
    return true;
    }

  // create data matrix and weights

  ST::string rname;
  ST::string wn;
  ST::string ifexpression;
  ST::string predictindicator;
  ST::string missingindicator;
  int wpos=-1;

  if(family.getvalue()=="multistate")
    {
    modelvarnamesv = modregmult.getModelVarnamesAsVector();
    wn = methods[4].get_weight_variable().to_bstr();
    if (wn.length() != 0)
      {
      modelvarnamesv.push_back(wn);
      wpos = modelvarnamesv.size()-1;
      }
    ifexpression = methods[4].getexpression();

    // variable 'state' anfügen
    statepos=-1;
    if(state.getvalue() != "")
      {
      modelvarnamesv.push_back(state.getvalue());
      statepos = modelvarnamesv.size()-1;
      }
    else
      {
      outerror("ERROR: Variable state has to be specified as a global option!\n");
      return(true);
      }
    // 'lefttrunc' anfügen
    lefttruncpos=-1;
    if(lefttrunc.getvalue() != "")
      {
      modelvarnamesv.push_back(lefttrunc.getvalue());
      lefttruncpos = modelvarnamesv.size()-1;
      }
    // variable 'globalfrailty' anfügen
    gfrailtypos=-1;
    if(globalfrailty.getvalue() != "")
      {
      modelvarnamesv.push_back(globalfrailty.getvalue());
      gfrailtypos = modelvarnamesv.size()-1;
      }
    }
  else
    {
    modelvarnamesv = modreg.getModelVarnamesAsVector();
    rname = modelvarnamesv[0].to_bstr();

    ifexpression = methods[0].getexpression();

    // Für Cox model: variable 'leftint' an 2. Stelle setzen
    leftintpos=-1;
    if(leftint.getvalue() != "")
      {
      modelvarnamesv.push_back(leftint.getvalue());
      leftintpos = modelvarnamesv.size()-1;
      }

    // Für Cox model: variable 'lefttrunc' an 3. Stelle setzen
    lefttruncpos=-1;
    if(lefttrunc.getvalue() != "")
      {
      modelvarnamesv.push_back(lefttrunc.getvalue());
      lefttruncpos = modelvarnamesv.size()-1;
      }

    // Für Binomial-Verteilung: variable 'binomweightpos' anfügen
    binomweightpos=-1;
    if(binomweight.getvalue() != "")
      {
      modelvarnamesv.push_back(binomweight.getvalue());
      binomweightpos = modelvarnamesv.size()-1;
      }

    // Für multinomiale Modelle: NA-Indikator anfügen
    if(naindicator.getvalue() != "")
      {
      ST::string test;
      if(naindicator.getvalue().length()>1500)
        {
        test = naindicator.getvalue().substr(naindicator.getvalue().length()-1500,1500);
        }
      else
        {
        outerror("ERROR: naindicator has to be category-specific\n");
        return true;
        }
      if(test=="_catspecific")
        {
        modelvarnamesv.push_back(naindicator.getvalue());
        }
      else
        {
        outerror("ERROR: naindicator has to be category-specific\n");
        return true;
        }
      }

    // Für multinomiale Modelle: Kategorienspezifische Kovariablen extrahieren

    if(family.getvalue()=="multinomialcatsp")
      {
      vector<ST::string> helpstring1;
      helpstring1.push_back(modelvarnamesv[0]);
      vector<ST::string> helpstring2;

      if ((datap->allexisting(helpstring1,helpstring2)) == false)
        {
        outerror("ERROR: variable " + helpstring2[0] + " is not existing\n");
        return(true);
        }

      datamatrix resphelp;
      datap->makematrix(modelvarnamesv[0],resphelp,ifexpression);

      // extract categories
      resphelp.sort(0,resphelp.rows()-1,0);
      allcats.erase(allcats.begin(),allcats.end());
      allcats.push_back(resphelp(0,0));
//      unsigned refpos;
      for(i=1; i<resphelp.rows(); i++)
        {
        if(resphelp(i,0)!=resphelp(i-1,0))
          {
          allcats.push_back(resphelp(i,0));
          }
        }

      ST::string test;
      vector<ST::string> modelvarnamesvhelp;
      for(i=0; i<modelvarnamesv.size(); i++)
        {
        if(modelvarnamesv[i].length()>1500)
          {
          test = modelvarnamesv[i].substr(modelvarnamesv[i].length()-1500,1500);
          if(test=="_catspecific")
            {
            test = modelvarnamesv[i].substr(0,modelvarnamesv[i].length()-1500);
            for(j=0; j<allcats.size(); j++)
              {
              modelvarnamesvhelp.push_back(test + ST::inttostring(allcats[j]));
              }
            }
          else
            {
            modelvarnamesvhelp.push_back(modelvarnamesv[i]);
            }
          }
        else
          {
          modelvarnamesvhelp.push_back(modelvarnamesv[i]);
          }
        }
      modelvarnamesv = modelvarnamesvhelp;
      }

    wn = methods[0].get_weight_variable().to_bstr();

    if (wn.length() != 0)
      {
      modelvarnamesv.push_back(wn);
      wpos = modelvarnamesv.size()-1;
      }

    }

  // testing, wether all variables specified are already existing
  vector<ST::string> notex;
  if ((datap->allexisting(modelvarnamesv,notex)) == false)
    {
    for (i=0;i<notex.size();i++)
      if (notex[i] != "const")
        {
        outerror("ERROR: variable " + notex[i] + " is not existing\n");
        failure = true;
        }
    if (failure)
      return true;
    }

  datap->makematrix(modelvarnamesv,D,ifexpression);

  // Check for error messages

  errormessages = datap->geterrormessages();

  if (!errormessages.empty())
    return true;

  // extract weights

  if (wpos==-1)
    {
    // No weights specified
    weight = datamatrix(D.rows(),1,1);
    }
  else
    {
    // check for correct distributions
/*    if(family.getvalue()=="cox" || family.getvalue()=="coxold")
      {
      outerror("ERROR: weight not allowed for family=cox\n");
      return true;
      }*/
    if(family.getvalue()=="multistate")
      {
      outerror("ERROR: weight not allowed for family=multistate\n");
      return true;
      }

    weight = D.getCol(wpos);

    // check for negative weights
    if(weight.min(0)<0)
      {
      outerror("ERROR: negative weights encountered\n");
      return true;
      }
    // only weights 1 or 0 allowed for multinomial responses and survival times
    if(family.getvalue()=="multinomial" ||
       family.getvalue()=="multinomialcatsp" ||
       family.getvalue()=="cumlogit" ||
       family.getvalue()=="cumprobit" ||
       family.getvalue()=="seqlogit" ||
       family.getvalue()=="seqprobit")
      {
      for(int i=0; i<weight.rows(); i++)
        {
        if(weight(i,0)!=0 && weight(i,0)!=1)
          {
          outerror("ERROR: weights have to equal 0 or 1 for multicategorical response\n");
          return true;
          }
        }
      }
    if(family.getvalue()=="cox" ||
       family.getvalue()=="coxold")
      {
      for(int i=0; i<weight.rows(); i++)
        {
        if(weight(i,0)!=0 && weight(i,0)!=1)
          {
          outerror("ERROR: weights have to equal 0 or 1 for survival times\n");
          return true;
          }
        }
      }

    }
  return false;
  }

// -----------------------------------------------------------------------------
// -------------------------- Create Response ----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_response(datamatrix & response, datamatrix & weight)
  {

  // extract response variable(s)
  if(family.getvalue()=="multistate")
    {
    vector<unsigned> rescol  = modregmult.getresponsecol();
    response = datamatrix(D.rows(),rescol.size());
    unsigned j,k,c;
    for(j=0;j<response.cols();j++)
      {
      c = rescol[j];
      for(k=0; k<response.rows(); k++)
        {
        response(k,j) = D(k,c);
        }
      }
    nrtransitions = rescol.size();
    }
  else
    {
    response = D.getCol(0);
    }

  // family=binomial
  if(family.getvalue()=="binomial" ||
     family.getvalue()=="binomialprobit" ||
     family.getvalue()=="binomialcomploglog" ||
     family.getvalue()=="binomialdispers" ||
     family.getvalue()=="binomialprobitdispers")
    {
    unsigned i;
    for(i=0; i<response.rows(); i++)
      {
      if(weight(i,0)>0)
        {
        if (response(i,0) != int(response(i,0)))
          {
          outerror("ERROR: response cannot be binomial; values must be integer numbers\n");
          return true;
          }

        if (response(i,0) < 0)
          {
          outerror("ERROR: response cannot be binomial; some values are negative\n");
          return true;
          }

        if (response(i,0) > weight(i,0))
          {
          outerror("ERROR: response cannot be binomial;\n");
          outerror("       number of successes larger than number of trials for some values\n");
          return true;
          }
        // Transform response for binomial families  (to match usual GLM definition)
        response(i,0) = response(i,0)/weight(i,0);
        if(binomweightpos>-1)
          {
          weight(i,0) *= D(i,binomweightpos);
          }
        }
      }
    }

  // family=gamma
  if(family.getvalue()=="gamma")
    {
    if(response.min(0)<0)
      {
      outerror("ERROR: response cannot be gamma distributed; some values are negative\n");
      return true;
      }
    }

  // family=poisson
  if(family.getvalue()=="poisson")
    {
    if(response.min(0)<0)
      {
      outerror("ERROR: response cannot be poisson distributed; some values are negative\n");
      return true;
      }
    }

  // family=cox
  if(family.getvalue()=="cox" || family.getvalue()=="coxold" ||
     family.getvalue()=="coxinterval")
    {
    unsigned i;
    for(i=0; i<response.rows(); i++)
      {
      if (response(i,0)!=0 && response(i,0)!=1)
        {
        outerror("ERROR: response must be either zero or one\n");
        return true;
        }
      }

    // check whether there is a baseline effect
    bool baselineexisting = false;
    for(i=0;i<terms.size();i++)
      {
      if(nonp_baseline.checkvector(terms,i) == true )
        {
        baselineexisting = true;
        }
      }
    if(baselineexisting == false)
      {
      outerror("ERROR: no baseline specified\n");
      return true;
      }
    }

  // check whether a baseline effect is specified if family != cox
  if(family.getvalue()!="cox" && family.getvalue()!="coxinterval" &&
     family.getvalue()!="coxold" && family.getvalue()!="multistate")
    {
    unsigned i;
    bool baselineexisting = false;
    for(i=0;i<terms.size();i++)
      {
      if(nonp_baseline.checkvector(terms,i) == true )
        {
        baselineexisting = true;
        }
      }
    if(baselineexisting == true)
      {
      outerror("ERROR: term baseline can only be used with family=cox\n");
      return true;
      }
    }

  // check whether leftint or lefttrunc are specified if family != cox
  if(family.getvalue()!="cox" && family.getvalue()!="coxinterval" &&
     family.getvalue()!="coxold" && family.getvalue()!="multistate")
    {
    if(leftint.getvalue() != "")
      {
      outerror("ERROR: variable leftint can only be used with family=cox\n");
      return true;
      }
    if(lefttrunc.getvalue() != "")
      {
      outerror("ERROR: variable lefttrunc can only be used with family=cox\n");
      return true;
      }
    }

  // family=multistate

  if(family.getvalue()!="multistate")
    {
    if(state.getvalue() != "")
      {
      outerror("ERROR: variable state can only be used with family=multistate\n");
      return true;
      }
    }
  // family=multinomial / family=cumlogit / family=cumprobit
  if (family.getvalue()=="multinomial" || family.getvalue()=="multinomialcatsp")
    {
    ismultinomial=true;
    }
  else
    {
    ismultinomial=false;
    }

  if (family.getvalue()=="multinomial" || family.getvalue()=="multinomialcatsp" ||
      family.getvalue()=="cumlogit" || family.getvalue()=="cumprobit" ||
      family.getvalue()=="seqlogit" || family.getvalue()=="seqprobit")
    {
    // extract categories
    datamatrix resphelp=response;
    resphelp.sort(0,resphelp.rows()-1,0);
    vector<double> categories;
    categories.push_back(resphelp(0,0));
    unsigned i;
    unsigned refpos;
    for(i=1; i<resphelp.rows(); i++)
      {
      if(resphelp(i,0)!=resphelp(i-1,0))
        {
        categories.push_back(resphelp(i,0));
        }
      }

    // check for proper categories
    if (categories.size() == 1)
      {
      outerror("ERROR: response variable does not vary\n");
      return true;
      }
    if (categories.size() > 1500)
      {
      outerror("ERROR: too many values for the response variable\n");
      return true;
      }

    // Define / extract the reference category
    if(family.getvalue()=="multinomial" || family.getvalue()=="multinomialcatsp")
      {
      refpos=0; //categories.size()-1;
      if(reference.changed())
        {
        bool existing = false;
        double ref=reference.getvalue();
        i=0;
        while(!existing && i<categories.size())
          {
          if(categories[i] == ref)
            {
            existing=true;
            refpos=i;
            }
          i++;
          }
        if(!existing)
          {
          outerror("ERROR: reference category is not existing\n");
          return true;
          }
        }
      }
    else
      {
      if(reference.changed())
        {
        outerror("ERROR: Option reference is not allowed in cumulative models\n");
        return true;
        }
      else
        {
        refpos=categories.size()-1;
        }
      }

    if(family.getvalue()=="multinomialcatsp")
      {
      // put reference category at last position in allcats
      int helpint = categories[refpos];//reference.getvalue();
      for(i=refpos; i<allcats.size(); i++)
        {
        allcats[i] = allcats[i+1];
        }
      allcats[allcats.size()-1] = helpint;

    // extract NA-Indicator

      naind = datamatrix(D.rows(),allcats.size(),0);
      if(naindicator.getvalue()!="" && family.getvalue()=="multinomialcatsp")
        {
        unsigned k,j;
        ST::string help = naindicator.getvalue().substr(0,naindicator.getvalue().length()-1500);
        for(k=0; k<allcats.size(); k++)
          {
          j = (help+ST::inttostring(allcats[k])).isinlist(modelvarnamesv);
          naind.putCol(k, D.getCol(j));
          }
        }
      }

    // define cats (categories without reference category)
    cats=datamatrix(categories.size()-1,1,0);
    unsigned j=0;
    for(i=0; i<categories.size(); i++)
      {
      if(i!=refpos)
        {
        cats(j,0)=categories[i];
        j++;
        }
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// --------------------------- Create_offset -----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_offset(datamatrix & o)
  {
  unsigned i;
  unsigned nroffs=0;
  unsigned index = 0;
  for(i=0;i<terms.size();i++)
    {
    if ( offset.checkvector(terms,i) == true)
      {
      nroffs ++;
      index = i;
      }
    }

  if(nroffs > 1)
    {
    outerror("ERROR: multiple offsets are not allowed\n");
    return true;
    }
  else if(nroffs==0)
    {
    // Default offset (=0)
    if(family.getvalue()=="multinomial" ||
       family.getvalue()=="multinomialcatsp" ||
       family.getvalue()=="cumlogit" ||
       family.getvalue()=="cumprobit" ||
       family.getvalue()=="seqlogit" ||
       family.getvalue()=="seqprobit")
      {
      o = datamatrix(cats.rows()*D.rows(),1,0);
      }
    else
      {
      o = datamatrix(D.rows(),1,0);
      }
    }
  else if(nroffs==1)
    {
    // check for right distributions
    if(family.getvalue()=="multinomial" ||
       family.getvalue()=="cumlogit" ||
       family.getvalue()=="cumprobit" ||
       family.getvalue()=="seqlogit" ||
       family.getvalue()=="seqprobit")
      {
      outerror("ERROR: offset not allowed for multinomial response\n");
      return true;
      }
    else if(family.getvalue()=="cox" || family.getvalue()=="coxold" ||
       family.getvalue()=="coxinterval")
      {
      outerror("ERROR: offset not allowed for family=cox\n");
      return true;
      }
    else if(family.getvalue()=="multinomialcatsp")
      {
      unsigned j, k;
      o = datamatrix(cats.rows()*D.rows(),1,0);
      datamatrix ohelp = datamatrix(allcats.size()*D.rows(),1,0);
      ST::string test ="test";
      if(terms[index].varnames[0].length()>1500)
        {
        test = terms[index].varnames[0].substr(terms[index].varnames[0].length()-1500,1500);
        }
      else
        {
        outerror("ERROR: offset has to be category-specific if family=multinomial\n");
        return true;
        }
      if(test=="_catspecific")
        // category specific covariates
        {
        test = terms[index].varnames[0].substr(0,terms[index].varnames[0].length()-1500);
        for(k=0; k<allcats.size(); k++)
          {
          j = (test+ST::inttostring(allcats[k])).isinlist(modelvarnamesv);
          ohelp.putRowBlock(k*D.rows(), (k+1)*D.rows(), D.getCol(j));
          }

        for(i=0; i<D.rows(); i++)
          {
          for(k=0; k<cats.rows(); k++)
            {
            o(i*cats.rows()+k,0) = ohelp(k*D.rows()+i,0)-ohelp(cats.rows()*D.rows()+i,0);
            }
          }

        terms[index].varnames[0] = test;
        test="_catspecific";
        }
      else
        {
        outerror("ERROR: offset has to be category-specific if family=multinomial\n");
        return true;
        }
      }
    else
      {
      // Extract offset if specified
      unsigned j = terms[index].varnames[0].isinlist(modelvarnamesv);
      if (o.rows() < D.rows())
        o = datamatrix(D.rows(),1,0);

      register unsigned k;
      double * worko = o.getV();
      double * workD = D.getV()+j;
      unsigned size = D.cols();
      for(k=0;k<D.rows();k++,worko++,workD+=size)
        *worko += *workD;
      }
    }

  return false;
  }

// -----------------------------------------------------------------------------
// --------------------------- Create_const -----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_const(const unsigned & collinpred)
  {
  unsigned i;
  int j, k, l, m;

  vector<ST::string> varnames;
  vector<ST::string> varnamesh =  fixedeffects.get_constvariables(terms);

  if(noconst.getvalue()==false)
    {
    varnames.push_back("const");
    }

  // add catspecific terms for ordinal and sequential models
  for(i=0;i<terms.size();i++)
    {
    if(fixed_catsp.checkvector(terms,i) == true)
      {
      if(family.getvalue()!="cumlogit" && family.getvalue()!="cumprobit" &&
         family.getvalue()!="seqlogit" && family.getvalue()!="seqprobit")
        {
        outerror("ERROR: category specific effects not allowed for family=" + family.getvalue() + "\n");
        return true;
        }
      varnames.push_back(terms[i].varnames[0]);
      }
    }

  for(i=0;i<varnamesh.size();i++)
    varnames.push_back(varnamesh[i]);

  unsigned nr = varnames.size();

  vector<bool>catsp(nr,false);
  catsp[0] = true;
  j=1;
  for(i=0;i<terms.size();i++)
    {
    if(fixed_catsp.checkvector(terms,i) == true)
      {
      catsp[j]=true;
      j++;
      }
    }

  vector<ST::string> varnames2 = varnames;
  // Category specific covariates for the multinomiallogit model
  if(family.getvalue()=="multinomialcatsp")
    {
    // set values in catsp (catsp[i]==true if corresponding covariate is category-specific)
    catsp[0]=false;
    ST::string test;
    for(i=0; i<varnames.size(); i++)
      {
      if(varnames[i].length()>1500)
        {
        test = varnames[i].substr(varnames[i].length()-1500,1500);
        if(test=="_catspecific")
          {
          catsp[i] = true;
          varnames2[i] = varnames[i].substr(0,varnames[i].length()-1500);
          }
        }
      }

    // modify variable names
    vector<ST::string> varnameshelp;
    for(i=0; i<varnames.size(); i++)
      {
      if(varnames[i].length()>1500)
        {
        test = varnames[i].substr(varnames[i].length()-1500,1500);
        if(test=="_catspecific")
          {
          test = varnames[i].substr(0,varnames[i].length()-1500);
          for(j=0; j<(int)allcats.size(); j++)
            {
            varnameshelp.push_back(test + ST::inttostring(allcats[j]));
            }
          }
        else
          {
          varnameshelp.push_back(varnames[i]);
          }
        }
      else
        {
        varnameshelp.push_back(varnames[i]);
        }
      }

    // redefine varnames
    varnames=varnameshelp;
    }

  ST::string title;
  ST::string pathconst;
  ST::string pathconstres;

  if (collinpred == 0)
    {
    title = "FixedEffects" + add_name;
#if defined(__BUILDING_LINUX)
    pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                           + add_name + "_FixedEffects" + ".raw";
#else
    pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                           + add_name + "_FixedEffects" + ".raw";
#endif
    pathconstres = outfile.getvalue() + add_name +
                     "_FixedEffects"  + ".res";
    }
  else
    {
    title = "FixedEffects"  "_" +
                            ST::inttostring(collinpred+1) + add_name;
#if defined(__BUILDING_LINUX)
    pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                           + add_name + "_FixedEffects"  +
                           "_" + ST::inttostring(collinpred+1) + ".raw";
#else
    pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                           + add_name + "_FixedEffects"  +
                           "_" + ST::inttostring(collinpred+1) + ".raw";
#endif
    pathconstres = outfile.getvalue() + add_name + "_FixedEffects" + "_" +
                   ST::inttostring(collinpred+1) + ".res";
    }

  if (pathconst.isvalidfile() == 1)
    {
    errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
    return true;
    }

  datamatrix X(D.rows(),varnames.size(),1);

  for(i=0;i<varnames.size();i++)
    {
    if (varnames[i] != "const")
      {
      j = varnames[i].isinlist(modelvarnamesv);

      if (j != -1)
        {
        unsigned l;
        double * workX=X.getV()+i;
        double * workD=D.getV()+j;
        for (l=0;l<X.rows();l++,workX+=X.cols(),workD+=D.cols())
          *workX = *workD;
        }
      }
    }

  // modify X matrix if category-specific covariates are present

  unsigned nrcatsp=0;
  if(family.getvalue()=="multinomialcatsp")
    {
    for(i=0; i<catsp.size(); i++)
      {
      if(catsp[i])
        nrcatsp++;
      }

    datamatrix Xhelp(D.rows(),nr + nrcatsp*(cats.rows()-1),0);
    k=0;
    m=0;
    for(i=0; i<varnames.size();)
      {
      if(catsp[m])
        {
        for(l=0; l<(int)cats.rows(); l++)
          {
          for(j=0; j<(int)D.rows(); j++)
            {
            Xhelp(j,k) = X(j,i+l)-X(j,i+cats.rows());
            }
          k++;
          }
        i += allcats.size();
        }
      else
        {
        for(j=0; j<(int)D.rows(); j++)
          {
          Xhelp(j,k) = X(j,i);
          }
        k++;
        i++;
        }
      m++;
      }
    X = Xhelp;
    }

  // define indicator for multinomial models with category specific covariates
  bool ismultinomialcatsp=false;
  if(family.getvalue()=="multinomialcatsp")
    {
    ismultinomialcatsp=true;
    }

  fcconst.push_back(FULLCOND_const(&generaloptions,X,title,0,pathconst,
                                   pathconstres,catsp,nrcatsp,nr,cats,
                                   ismultinomialcatsp));

  fcconst[fcconst.size()-1].init_names(varnames2);
  fcconst[fcconst.size()-1].set_fcnumber(fullcond.size());
  fullcond.push_back(&fcconst[fcconst.size()-1]);

  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_nonprw1rw2 ---------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_nonprw1rw2(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;

  double hd;
  double lambda, startlambda;
  bool catsp;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonprw1rw2.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "rw1")
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[2]).strtodouble(hd);
      startlambda = hd;
      if(terms[i].options[3] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }


      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_rw.raw","_rw.res","_rw");

      fcnonpgaussian.push_back(FULLCOND_nonp_gaussian(&generaloptions,D.getCol(j),
                         unsigned(maxint.getvalue()),type,
                         title,pathres,lambda,startlambda,catsp));

      fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// --------------------- create_nonprw1rw2_varcoef -----------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_nonprw1rw2_varcoef(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;

  double hd;
  double lambda, startlambda;
  bool catsp, ctr;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonprw1rw2_varcoef.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "varcoeffrw1")
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[2]).strtodouble(hd);
      startlambda = hd;
      if(terms[i].options[3] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      if(terms[i].options[4] == "true")
        ctr=true;
      else
        ctr=false;

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],"_rw.raw","_rw.res","_rw");

      fcnonpgaussian.push_back(FULLCOND_nonp_gaussian(&generaloptions,
                         D.getCol(j2),D.getCol(j1),
                         unsigned(maxint.getvalue()),type,
                         title,pathres,lambda,startlambda,catsp,ctr));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_nonpseason ---------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_nonpseason(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;

  long h;
  double hd;
  double lambda, startlambda;
  unsigned per;
  bool catsp;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpseason.checkvector(terms,i) == true )
      {
      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      per = unsigned(h);
      f = (terms[i].options[2]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[3]).strtodouble(hd);
      startlambda = hd;
      if(terms[i].options[4] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }

      if (f==1)
        return true;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_season.raw","_season.res","_season");

      fcnonpgaussian.push_back(
      FULLCOND_nonp_gaussian(&generaloptions,D.getCol(j),
                                         unsigned(maxint.getvalue()),
                                         MCMC::seasonal,title,pathres,
                                         lambda,startlambda,catsp,per));

      fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_nonpseason_varcoef ---------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_nonpseason_varcoef(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double hd;
  double lambda, startlambda;
  unsigned per;
  bool catsp;
  int f;

  unsigned i;
  int j1, j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpseason_varcoef.checkvector(terms,i) == true )
      {
      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtolong(h);
      per = unsigned(h);
      f = (terms[i].options[2]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[3]).strtodouble(hd);
      startlambda = hd;
      if(terms[i].options[4] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }

      if (f==1)
        return true;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],"_season.raw","_season.res",
                 "_season");

      fcnonpgaussian.push_back(
      FULLCOND_nonp_gaussian(&generaloptions,D.getCol(j2),D.getCol(j1),
                                         unsigned(maxint.getvalue()),
                                         MCMC::seasonal,title,pathres,
                                         lambda,startlambda,catsp,per));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// --------------------------- create_spatial ---------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_spatial(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  double hd;
  int f;
  double lambda, startlambda;
  bool catsp;
  unsigned i;
  int j;

  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial.checkvector(terms,i) == true )
      {
      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[1],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          {
          if ((terms[i].options[1] == "") || (terms[i].options[1] == " "))
            outerror("ERROR: map object must be specified to estimate a spatial effect\n");
          else
            outerror("ERROR: map object " + terms[i].options[1] + " is not existing\n");
          }
        else
          outerror("ERROR: " + terms[i].options[1] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();

      if (m.isconnected()==false)
        {
        outerror("ERROR: map is disconnected, spatial effect cannot be estimated\n");
        return true;
        }

      f = (terms[i].options[2]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[3]).strtodouble(hd);
      startlambda = hd;
      if(terms[i].options[4] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }

      if (f==1)
        return true;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_spatial.raw","_spatial.res","_spatial");

      fcnonpgaussian.push_back(
      FULLCOND_nonp_gaussian(&generaloptions,D.getCol(j),m,terms[i].options[1],
                             title,pathnonp,pathres,lambda,startlambda,catsp)
                           );

      if (fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size() > 0)
        {
        unsigned i;
        for(i=0;i<fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size();i++)
          errormessages.push_back(fcnonpgaussian[fcnonpgaussian.size()-1].get_errors()[i]);
        return true;
        }

      fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_spatial_varcoef ----------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_spatial_varcoef(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  double hd;
  int f;
  double lambda, startlambda;
  bool catsp, ctr;
  unsigned i;
  int j1, j2;

  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_varcoef.checkvector(terms,i) == true )
      {
      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[1],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          {
          if ((terms[i].options[1] == "") || (terms[i].options[1] == " "))
            outerror("ERROR: map object must be specified to estimate a spatial effect\n");
          else
            outerror("ERROR: map object " + terms[i].options[1] + " is not existing\n");
          }
        else
          outerror("ERROR: " + terms[i].options[1] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();

      f = (terms[i].options[2]).strtodouble(hd);
      lambda = hd;
      f = (terms[i].options[3]).strtodouble(hd);
      startlambda = hd;
      if(terms[i].options[4] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      if(terms[i].options[5] == "true")
        ctr=true;
      else
        ctr=false;

      if (f==1)
        return true;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],"_spatial.raw","_spatial.res","_spatial");

      fcnonpgaussian.push_back(
      FULLCOND_nonp_gaussian(&generaloptions,D.getCol(j2),D.getCol(j1),m,
                             terms[i].options[1],title,pathnonp,pathres,lambda,
                             startlambda,catsp,ctr)
                           );

      if (fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size() > 0)
        {
        unsigned i;
        for(i=0;i<fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size();i++)
          errormessages.push_back(fcnonpgaussian[fcnonpgaussian.size()-1].get_errors()[i]);
        return true;
        }

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);
      fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// -------------------------- create_pspline -----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_pspline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;
  long h;
  unsigned degree,nrknots;
  double lambda, startlambda, lowergrid, uppergrid, lowerknot, upperknot, refval;
  bool catsp;
  int f, gridsize;

  unsigned i,k;
  int j;

  datamatrix data;

  for(i=0;i<terms.size();i++)
    {
    if ( nonppspline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "psplinerw1")
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      ST::string test ="test";
      if(terms[i].varnames[0].length()>1500)
        {
        test = terms[i].varnames[0].substr(terms[i].varnames[0].length()-1500,1500);
        }
      if(test=="_catspecific")
        // category specific covariates
        {
        test = terms[i].varnames[0].substr(0,terms[i].varnames[0].length()-1500);
        data = datamatrix(allcats.size()*D.rows(),1,0);
        for(k=0; k<allcats.size(); k++)
          {
          j = (test+ST::inttostring(allcats[k])).isinlist(modelvarnamesv);
          data.putRowBlock(k*D.rows(), (k+1)*D.rows(), D.getCol(j));
          }
        terms[i].varnames[0] = test;
        test="_catspecific";
        }
      else
        // no category specific covariates
        {
        j = terms[i].varnames[0].isinlist(modelvarnamesv);
        data = D.getCol(j);
        }

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[4]).strtolong(h);
      gridsize = unsigned(h);
      f = (terms[i].options[7]).strtodouble(startlambda);
      if(terms[i].options[8] == "true" || test=="_catspecific")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      f = (terms[i].options[9]).strtodouble(lowergrid);
      f = (terms[i].options[10]).strtodouble(uppergrid);
      f = (terms[i].options[11]).strtodouble(lowerknot);
      f = (terms[i].options[12]).strtodouble(upperknot);
      f = (terms[i].options[13]).strtodouble(refval);

      if (f==1)
        return true;

      MCMC::knotpos po;
      if (knots.getvalue() == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline.raw","_pspline.res","_pspline");

      fcpspline.push_back( spline_basis(&generaloptions,
                                              data,
                                              nrknots,
                                              degree,
                                              po,
                                              type,
                                              title,
                                              pathnonp,
                                              pathres,
                                              lambda,
                                              startlambda,
                                              catsp,
                                              lowergrid,
                                              uppergrid,
                                              lowerknot,
                                              upperknot,
                                              gridsize,
                                              refval
                                             )
                           );
      fcpspline[fcpspline.size()-1].init_name(terms[i].varnames[0]);
      fcpspline[fcpspline.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpspline[fcpspline.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_varcoeffpspline ----------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_varcoeffpspline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string monotone;
  long h;
  unsigned degree,nrknots;
  double lambda,startlambda,refval;
  bool catsp, ctr;
  int f;

//  unsigned i;
  int j1,j2;

  datamatrix data1;
  datamatrix data2;
  for(unsigned i=0;i<terms.size();i++)
    {
    if ( nonpvarcoeffpspline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "varpsplinerw1")
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      ST::string test ="test";
      if(terms[i].varnames[0].length()>1500)
        {
        test = terms[i].varnames[0].substr(terms[i].varnames[0].length()-1500,1500);
        }
      if(test=="_catspecific")
        // category specific interacting variable
        {
        test = terms[i].varnames[0].substr(0,terms[i].varnames[0].length()-1500);
        data1 = datamatrix(allcats.size()*D.rows(),1,0);
        data2 = datamatrix(allcats.size()*D.rows(),1,0);
        for(unsigned k=0; k<allcats.size(); k++)
          {
          j1 = (test+ST::inttostring(allcats[k])).isinlist(modelvarnamesv);
          data1.putRowBlock(k*D.rows(), (k+1)*D.rows(), D.getCol(j1));
          data2.putRowBlock(k*D.rows(), (k+1)*D.rows(), D.getCol(j2));
          }
        terms[i].varnames[0] = test;
        test="_catspecific";
        }
      else
        // no category specific interacting variable
        {
        j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
        data1 = D.getCol(j1);
        data2 = D.getCol(j2);
        }

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[4]).strtodouble(startlambda);
      if(terms[i].options[5] == "true" || test=="_catspecific")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      if(terms[i].options[6] == "true")
        ctr=true;
      else
        ctr=false;

      f = (terms[i].options[7]).strtodouble(refval);

      if (f==1)
        return true;

      MCMC::knotpos po;
      if (knots.getvalue() == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_pspline.raw","_pspline.res","_pspline");

      fcpspline.push_back( spline_basis(&generaloptions,
                                              data2,
                                              data1,
                                              nrknots,
                                              degree,
                                              po,
                                              type,
                                              title,
                                              pathnonp,
                                              pathres,
                                              lambda,
                                              startlambda,
                                              catsp,
                                              ctr,
                                              refval
                                             )
                           );
      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcpspline[fcpspline.size()-1].init_names(na);
      fcpspline[fcpspline.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpspline[fcpspline.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------- create_interactionpspline -------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_interactionspspline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double lambda,startlambda;
  unsigned nrknots,degree,gridsizex,gridsizey;
  bool catsp;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpinteractpspline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "pspline2dimrw1")
        type = MCMC::mrflinear;
      else if (terms[i].options[0] == "pspline2dimrw2")
        type = MCMC::mrfquadratic8;
      else if( terms[i].options[0] =="pspline2dimbiharmonic")
        type = MCMC::mrfquadratic12;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[4]).strtodouble(startlambda);
      if(terms[i].options[5] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      f = (terms[i].options[6]).strtolong(h);
      gridsizex = unsigned(h);
      f = (terms[i].options[7]).strtolong(h);
      gridsizey = unsigned(h);

      if (f==1)
        return true;

      ST::string title;

      ST::string help  = terms[i].varnames[0] + "_" + terms[i].varnames[1];

      make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline.raw","_pspline.res","_pspline");

      fcpsplinesurf.push_back(
      spline_basis_surf(&generaloptions,
                                      D.getCol(j1),
                                      D.getCol(j2),
                                      nrknots,degree,type,
                                      title,
                                      pathnonp,
                                      pathres,
                                      lambda,
                                      startlambda,
                                      catsp,
                                      gridsizex,gridsizey
                                      ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
      fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ---------------------- create_varcoeffinteractionpspline --------------------
// -----------------------------------------------------------------------------

bool remlreg::create_varcoeffinteractionspspline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double lambda,startlambda;
  unsigned nrknots,degree;
  bool catsp, ctr;
  int f;

  unsigned i;
  int j1,j2,j3;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpvarcoeffinteractpspline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "varpspline2dimrw1")
        type = MCMC::mrflinear;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);  //interaction variable
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);
      j3 = terms[i].varnames[2].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[4]).strtodouble(startlambda);
      if(terms[i].options[5] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      if(terms[i].options[6] == "true")
        ctr=true;
      else
        ctr=false;

      if (f==1)
        return true;

      ST::string title;

      ST::string help  = terms[i].varnames[1] + "_" + terms[i].varnames[2];

      make_paths(collinpred,pathnonp,pathres,title,help,terms[i].varnames[0],
                 "_pspline.raw","_pspline.res","_pspline");

      fcpsplinesurf.push_back(
      spline_basis_surf(&generaloptions,
                                      D.getCol(j1),
                                      D.getCol(j2),
                                      D.getCol(j3),
                                      nrknots,degree,type,
                                      title,
                                      pathnonp,
                                      pathres,
                                      lambda,
                                      startlambda,
                                      catsp,
                                      ctr
                                      ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[2]);
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
      fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_geospline ----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_geospline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;

  long h;
  double lambda,startlambda;
  unsigned nrknots,degree,gridsizex,gridsizey;
  bool catsp;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpgeospline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "geosplinerw1")
        type = MCMC::mrflinear;
      else if (terms[i].options[0] == "geosplinerw2")
        type = MCMC::mrfquadratic8;
      else if (terms[i].options[0] == "geosplinebiharmonic")
        type = MCMC::mrfquadratic12;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      mapobject * mapp;                           // pointer to mapobject
      statobject * s;
      int objpos = findstatobject(*statobj,terms[i].options[4],"map");
      if (objpos >= 0)
        {
        s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[4] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[4] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();

      if(!m.centroids_existing())
        {
        outerror("ERROR: map object does not contain centroids\n");
        return true;
        }

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[5]).strtodouble(startlambda);
      if(terms[i].options[6] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      f = (terms[i].options[7]).strtolong(h);
      gridsizex = unsigned(h);
      f = (terms[i].options[8]).strtolong(h);
      gridsizey = unsigned(h);

      if (f==1)
        return true;

//      ST::string test = terms[i].varnames[0];
//      ST::string test = "region";

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline.raw","_geospline.res","_geospline");

      fcpsplinesurf.push_back(
      spline_basis_surf(&generaloptions,
                                      D.getCol(j),m,terms[i].options[4],
                                      nrknots,degree,
                                      type,
                                      title,
                                      pathnonp,
                                      pathres,
                                      lambda,
                                      startlambda,
                                      catsp,
                                      gridsizex,gridsizey
                                      ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[0]);
      fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
      fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------- create_geospline_varcoef --------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_geospline_varcoeff(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double lambda,startlambda;
  unsigned nrknots,degree;
  bool catsp, ctr;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpvarcoeffgeospline.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type = MCMC::mrflinear;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[5]).strtodouble(startlambda);
      if(terms[i].options[6] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      if(terms[i].options[7] == "true")
        ctr=true;
      else
        ctr=false;

      if (f==1)
        return true;

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[4],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[4] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[4] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();

      if(!m.centroids_existing())
        {
        outerror("ERROR: map object does not contain centroids\n");
        return true;
        }

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_geospline.raw","_geospline.res","_geospline");

      fcpsplinesurf.push_back(
      spline_basis_surf(&generaloptions,
                                      D.getCol(j1),D.getCol(j2),
                                      m,terms[i].options[4],
                                      nrknots,degree,
                                      type,
                                      title,
                                      pathnonp,
                                      pathres,
                                      lambda,
                                      startlambda,
                                      catsp,
                                      ctr
                                      ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
      fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_kriging --------------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_kriging(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double nu,maxdist,p,q,lambda,startlambda;
  unsigned nrknots, maxsteps,gridsizex,gridsizey;
  bool full, catsp;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_kriging.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "kriging")
        type = MCMC::kriging;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[2]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[3]).strtodouble(maxdist);
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mit maximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;
          }
        }
      if(terms[i].options[4] == "true")
        {
        full=true;
        }
      else
        {
        full=false;
        }
      f = (terms[i].options[6]).strtodouble(p);
      f = (terms[i].options[7]).strtodouble(q);
      f = (terms[i].options[8]).strtolong(h);
      maxsteps = unsigned(h);
      f = (terms[i].options[9]).strtodouble(lambda);
      f = (terms[i].options[10]).strtodouble(startlambda);
      if(terms[i].options[11] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      f = (terms[i].options[12]).strtolong(h);
      gridsizex = unsigned(h);
      f = (terms[i].options[13]).strtolong(h);
      gridsizey = unsigned(h);

      if (f==1)
        return true;

      // read knots from a specified dataset object
      datamatrix knotdata;
      if(terms[i].options[5]!="" && !full)
        {
        dataobject * datap;                           // pointer to datsetobject
        int objpos = findstatobject(*statobj,terms[i].options[5],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>2)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " contains more than two variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[5] + " is not existing\n");
          return true;
          }
        list<ST::string> knotnames = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(knotnames,knotdata,expr);
        }
      else
        {
        knotdata = datamatrix(1,1,0);
        }

      ST::string title;

      ST::string help  = terms[i].varnames[0] + "_" + terms[i].varnames[1];

      make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_kriging.raw","_kriging.res","_kriging");

      fckriging.push_back(
      FULLCOND_kriging(&generaloptions,
              D.getCol(j1),
              D.getCol(j2),
              knotdata,
              nrknots,nu,maxdist,p,q,maxsteps,full,type,
              title,
              pathnonp,
              pathres,
              lambda,
              startlambda,
              catsp,
              gridsizex,gridsizey
              ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_kriging_1dim ---------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_kriging_1dim(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  double nu,maxdist,lambda,startlambda;
  bool catsp;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonp_kriging.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "1dimkriging")
        type = MCMC::kriging;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[2]).strtodouble(maxdist);
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;
          }
        }
      f = (terms[i].options[3]).strtodouble(lambda);
      f = (terms[i].options[4]).strtodouble(startlambda);
      if(terms[i].options[5] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }

      if (f==1)
        return true;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_kriging.raw","_kriging.res","_kriging");

      fckriging.push_back(
      FULLCOND_kriging(&generaloptions,
              D.getCol(j),
              nu,maxdist,type,
              title,
              pathnonp,
              pathres,
              lambda,
              startlambda,
              catsp
              ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[0]);
      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_kriging_varcoeff -----------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_kriging_varcoeff(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double nu,maxdist,p,q,lambda,startlambda;
  unsigned nrknots, maxsteps;
  bool full, catsp;
  int f;

  unsigned i;
  int j1,j2,j3;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_kriging_varcoeff.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "varkriging")
        type = MCMC::kriging;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); //interaction variable
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);
      j3 = terms[i].varnames[2].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[2]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[3]).strtodouble(maxdist);
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;
          }
        }
      if(terms[i].options[4] == "true")
        {
        full=true;
        }
      else
        {
        full=false;
        }

      f = (terms[i].options[6]).strtodouble(p);
      f = (terms[i].options[7]).strtodouble(q);
      f = (terms[i].options[8]).strtolong(h);
      maxsteps = unsigned(h);

      f = (terms[i].options[9]).strtodouble(lambda);
      f = (terms[i].options[10]).strtodouble(startlambda);
      if(terms[i].options[11] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }

      if (f==1)
        return true;

      datamatrix knotdata;
      if(terms[i].options[5]!="" && !full)
        {
        dataobject * datap;                           // pointer to datsetobject
        int objpos = findstatobject(*statobj,terms[i].options[5],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>2)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " contains more than two variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[5] + " is not existing\n");
          return true;
          }
        list<ST::string> knotnames = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(knotnames,knotdata,expr);
        }
      else
        {
        knotdata = datamatrix(1,1,0);
        }

      ST::string title;

      ST::string help  = terms[i].varnames[1] + "_" + terms[i].varnames[2];

      make_paths(collinpred,pathnonp,pathres,title,help,terms[i].varnames[0],
                 "_kriging.raw","_kriging.res","_kriging");

      fckriging.push_back(
      FULLCOND_kriging(&generaloptions,
              D.getCol(j1),
              D.getCol(j2),
              D.getCol(j3),
              knotdata,
              nrknots,nu,maxdist,p,q,maxsteps,full,type,
              title,
              pathnonp,
              pathres,
              lambda,
              startlambda,
              catsp
              ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[2]);
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_geokriging -----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_geokriging(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double nu,maxdist,p,q,lambda,startlambda;
  unsigned nrknots, maxsteps, gridsizex, gridsizey;
  bool full, catsp;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_geokriging.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "geokriging")
        type = MCMC::kriging;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[2]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[3]).strtodouble(maxdist);
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;
          }
        }
      if(terms[i].options[4] == "true")
        {
        full=true;
        }
      else
        {
        full=false;
        }
      f = (terms[i].options[6]).strtodouble(p);
      f = (terms[i].options[7]).strtodouble(q);
      f = (terms[i].options[8]).strtolong(h);
      maxsteps = unsigned(h);
      f = (terms[i].options[9]).strtodouble(lambda);
      f = (terms[i].options[10]).strtodouble(startlambda);
      if(terms[i].options[12] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      f = (terms[i].options[13]).strtolong(h);
      gridsizex = unsigned(h);
      f = (terms[i].options[14]).strtolong(h);
      gridsizey = unsigned(h);

      if (f==1)
        return true;

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[11],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[11] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[11] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();
      if(!m.centroids_existing())
        {
        outerror("ERROR: map object does not contain centroids\n");
        return true;
        }

      datamatrix knotdata;
      if(terms[i].options[5]!="" && !full)
        {
        dataobject * datap;                           // pointer to datasetobject
        int objpos = findstatobject(*statobj,terms[i].options[5],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>2)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " contains more than two variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[5] + " is not existing\n");
          return true;
          }
        list<ST::string> knotnames = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(knotnames,knotdata,expr);
        }
      else
        {
        knotdata = datamatrix(1,1,0);
        }

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geokriging.raw","_geokriging.res","_geokriging");

      fckriging.push_back(
      FULLCOND_kriging(&generaloptions,
              D.getCol(j),m,terms[i].options[11],
              knotdata,
              nrknots,nu,maxdist,p,q,maxsteps,full,type,
              title,
              pathnonp,
              pathres,
              lambda,
              startlambda,
              catsp,
              gridsizex,gridsizey
              ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[0]);
      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_geokriging_varcoeff --------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_geokriging_varcoeff(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  long h;
  double nu,maxdist,p,q,lambda,startlambda;
  unsigned nrknots, maxsteps;
  bool full, catsp;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_geokriging_varcoeff.checkvector(terms,i) == true )
      {
      MCMC::fieldtype type;
      if (terms[i].options[0] == "vargeokriging")
        type = MCMC::kriging;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[2]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[3]).strtodouble(maxdist);
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;
          }
        }
      if(terms[i].options[4] == "true")
        {
        full=true;
        }
      else
        {
        full=false;
        }
      f = (terms[i].options[6]).strtodouble(p);
      f = (terms[i].options[7]).strtodouble(q);
      f = (terms[i].options[8]).strtolong(h);
      maxsteps = unsigned(h);
      f = (terms[i].options[9]).strtodouble(lambda);
      f = (terms[i].options[10]).strtodouble(startlambda);
      if(terms[i].options[12] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }

      if (f==1)
        return true;

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[11],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[11] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[11] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();
      if(!m.centroids_existing())
        {
        outerror("ERROR: map object does not contain centroids\n");
        return true;
        }

      datamatrix knotdata;
      if(terms[i].options[5]!="" && !full)
        {
        dataobject * datap;                           // pointer to datasetobject
        int objpos = findstatobject(*statobj,terms[i].options[5],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>2)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " contains more than two variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[5] + " is not existing\n");
          return true;
          }
        list<ST::string> knotnames = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(knotnames,knotdata,expr);
        }
      else
        {
        knotdata = datamatrix(1,1,0);
        }

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_geokriging.raw","_geokriging.res","_geokriging");

      fckriging.push_back(
      FULLCOND_kriging(&generaloptions,
              D.getCol(j1),D.getCol(j2),m,terms[i].options[11],
              knotdata,
              nrknots,nu,maxdist,p,q,maxsteps,full,type,
              title,
              pathnonp,
              pathres,
              lambda,
              startlambda,
              catsp
              ));

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_baseline -------------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_baseline(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;
  long h;
  unsigned degree,nrknots,tgrid,nrquant,nrbetween;
  double lambda, startlambda, refval;
  bool catsp;
  int f, gridsize;
  MCMC::knotpos gridpo;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonp_baseline.checkvector(terms,i) == true )
      {
      if(family.getvalue()!="multistate" && fcbaseline.size()>0)
        {
        outerror("ERROR: More than one baseline term specified!\n");
        return true;
        }

      MCMC::fieldtype type;
      if (terms[i].options[0] == "baseline")
        type = MCMC::RW2;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      // read lower interval boundary
      datamatrix lowerint;
      if(leftint.getvalue()!="")
        {
        lowerint = D.getCol(leftintpos);
        }
      else
        {
        lowerint = datamatrix(1,1,0);
        }

      // read left truncation time
      datamatrix lowertrunc;
      if(lefttrunc.getvalue()!="")
        {
        lowertrunc = D.getCol(lefttruncpos);
        }
      else
        {
        lowertrunc = datamatrix(1,1,0);
        }

/*      datamatrix lower;
      if(terms[i].options[9]!="")
        {
        dataobject * datap;                           // pointer to datsetobject
        int objpos = findstatobject(*statobj,terms[i].options[9],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[9] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>1)
            {
            outerror("ERROR: dataset object " + terms[i].options[9] + " contains more than one variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[9] + " is not existing\n");
          return true;
          }
        list<ST::string> lowname = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(lowname,lower,expr);
        }
      else
        {
        lower = datamatrix(1,1,0);
        }*/

      f = (terms[i].options[1]).strtolong(h);
      degree = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      nrknots = unsigned(h);
      f = (terms[i].options[3]).strtolong(h);
      tgrid = unsigned(h);
      if(terms[i].options[4] == "equidistant")
        {
        gridpo = MCMC::equidistant;
        }
      else if(terms[i].options[4] == "quantiles")
        {
        gridpo = MCMC::quantiles;
        }
      else
        {
        gridpo = MCMC::all;
        }
      f = (terms[i].options[5]).strtolong(h);
      nrquant = unsigned(h);
      f = (terms[i].options[6]).strtolong(h);
      nrbetween = unsigned(h);
      f = (terms[i].options[7]).strtodouble(lambda);
      f = (terms[i].options[8]).strtodouble(startlambda);
      if(terms[i].options[10] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      f = (terms[i].options[11]).strtolong(h);
      gridsize = int(h);
      f = (terms[i].options[12]).strtodouble(refval);

      if (f==1)
        return true;

      MCMC::knotpos po = MCMC::equidistant;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_logbaseline.raw","_logbaseline.res","_logbaseline");

      fcbaseline.push_back( baseline_reml(&generaloptions,
                                              D.getCol(j),
                                              lowerint,
                                              lowertrunc,
                                              nrknots,
                                              degree,
                                              tgrid,
                                              nrquant,
                                              nrbetween,
                                              po,
                                              type,
                                              title,
                                              pathnonp,
                                              pathres,
                                              lambda,
                                              startlambda,
                                              gridpo,
                                              gridsize,
                                              catsp,
                                              refval
                                             )
                           );

      fcbaseline[fcbaseline.size()-1].set_fcnumber(fullcond.size());
      fcbaseline[fcbaseline.size()-1].init_name(terms[i].varnames[0]);
      fullcond.push_back(&fcbaseline[fcbaseline.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_baseline_varcoeff-----------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_baseline_varcoeff(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;
  double lambda, startlambda,refval;
  unsigned degree,nrknots,tgrid;
  bool catsp;
  int f, gridsize;
  long h;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonp_baseline_varcoeff.checkvector(terms,i) == true )
      {
      if(fcbaseline.size()<1)
        {
        outerror("ERROR: Time-varying effects without baseline effect!\n");
        return true;
        }

      MCMC::fieldtype type = MCMC::RW2;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtodouble(lambda);
      f = (terms[i].options[2]).strtodouble(startlambda);
      if(terms[i].options[3] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }
      f = (terms[i].options[4]).strtolong(h);
      gridsize = int(h);
      f = (terms[i].options[5]).strtodouble(refval);

      if (f==1)
        return true;

      tgrid = fcbaseline[0].get_tgrid();
      degree = fcbaseline[0].get_degree();
      nrknots = fcbaseline[0].get_nrknots();

      MCMC::knotpos po = MCMC::equidistant;

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_logbaseline.raw","_logbaseline.res","_logbaseline");

      fcbaseline_varcoeff.push_back( baseline_reml(&generaloptions,
                                              D.getCol(j2),
                                              D.getCol(j1),
                                              nrknots,
                                              degree,
                                              tgrid,
                                              po,
                                              type,
                                              title,
                                              pathnonp,
                                              pathres,
                                              lambda,
                                              startlambda,
                                              gridsize,
                                              catsp,
                                              refval
                                             )
                           );

      vector<ST::string> na;
      na.push_back(terms[i].varnames[1]);
      na.push_back(terms[i].varnames[0]);
      fcbaseline_varcoeff[fcbaseline_varcoeff.size()-1].init_names(na);
      fcbaseline_varcoeff[fcbaseline_varcoeff.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcbaseline_varcoeff[fcbaseline_varcoeff.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ------------------------ create_random --------------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_random(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;
  double lambda, startlambda;
  bool catsp;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( randomeff.checkvector(terms,i) == true )
      {
      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtodouble(lambda);
      f = (terms[i].options[2]).strtodouble(startlambda);
      if(terms[i].options[3] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_random.raw","_random.res","_random");

      fcrandom.push_back(FULLCOND_random(&generaloptions,
                         D.getCol(j),title,pathnonp,
                         pathres,lambda,startlambda,catsp));

      fcrandom[fcrandom.size()-1].init_name(terms[i].varnames[0]);
      fcrandom[fcrandom.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcrandom[fcrandom.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------- create_randomslope ----------------------------------
// -----------------------------------------------------------------------------

bool remlreg::create_randomslope(const unsigned & collinpred)
  {
  ST::string pathnonp;
  ST::string pathres;
  ST::string title;
  unsigned i,j;
  int j1,j2;
  double lambda, startlambda;
  bool catsp;
  int f;
  for(i=0;i<terms.size();i++)
    {
    if ( randomeffslope.checkvector(terms,i) == true )
      {
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

      datamatrix intvar;
      ST::string test ="test";
      if(terms[i].varnames[0].length()>1500)
        {
        test = terms[i].varnames[0].substr(terms[i].varnames[0].length()-1500,1500);
        }
      if(test=="_catspecific")
        {
        test = terms[i].varnames[0].substr(0,terms[i].varnames[0].length()-1500);
        intvar = datamatrix(D.rows(),allcats.size(),0);
        for(j=0; j<allcats.size(); j++)
          {
          j1 = (test+ST::inttostring(allcats[j])).isinlist(modelvarnamesv);
          intvar.putCol(j,D.getCol(j1));
          }
        terms[i].varnames[0] = terms[i].varnames[0].substr(0,terms[i].varnames[0].length()-1500);
        }
      else
        {
        j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
        intvar = D.getCol(j1);
        }

      f = (terms[i].options[1]).strtodouble(lambda);
      f = (terms[i].options[2]).strtodouble(startlambda);
      if(terms[i].options[3] == "true")
        {
        catsp=true;
        }
      else
        {
        catsp=false;
        }

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_random.raw","_random.res","_random");

      fcrandom.push_back(FULLCOND_random(&generaloptions,intvar,
                          D.getCol(j2),title,pathnonp,pathres,lambda,
                          startlambda,catsp));

      vector<ST::string> varnameshelp;
      varnameshelp.push_back(terms[i].varnames[1]);
      varnameshelp.push_back(terms[i].varnames[0]);
      fcrandom[fcrandom.size()-1].init_names(varnameshelp);
      fcrandom[fcrandom.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcrandom[fcrandom.size()-1]);
      }
    }
  return false;
  }

// -----------------------------------------------------------------------------
// ----------------------------- remlrun ---------------------------------------
// -----------------------------------------------------------------------------

void remlrun(remlreg & b)
  {
  b.resultsyesno = false;
  b.terms = b.modreg.getterms();

  b.describetext.erase(b.describetext.begin(),b.describetext.end());
  b.describetext.push_back("LAST ESTIMATED MODEL: \n");
  b.describetext.push_back("\n");
  b.describetext.push_back(b.modreg.getModelText());
  b.describetext.push_back("\n");

  b.clear();

  #if defined(JAVA_OUTPUT_WINDOW)
    bool failure = b.adminb_p->breakcommand();
  #else
    bool failure = false;
  #endif

  b.generaloptions = MCMCoptions(
  #if defined(JAVA_OUTPUT_WINDOW)
  b.adminb_p,
  #endif
  12000,2000,100,b.logout,b.level1.getvalue(),
                               b.level2.getvalue());

  ST::string header;
  bool dispers;

// Read design matrix, compute weights, etc.
  datamatrix weight;
  if (!failure)
    failure = b.create_data(weight);

// Define and check response
  datamatrix response;
  if(!failure)
    failure = b.create_response(response,weight);

// Compute offset
  datamatrix offset;
  if(!failure)
    failure = b.create_offset(offset);

// Compute different model terms
  if (!failure)
    failure = b.create_const(0);
  if( !failure)
    failure = b.create_baseline(0);
  if( !failure)
    failure = b.create_baseline_varcoeff(0);
  if (!failure)
    failure = b.create_nonprw1rw2(0);
  if (!failure)
    failure = b.create_nonprw1rw2_varcoef(0);
  if (!failure)
    failure = b.create_pspline(0);
  if (!failure)
    failure = b.create_nonpseason(0);
  if (!failure)
    failure = b.create_nonpseason_varcoef(0);
  if (!failure)
    failure = b.create_spatial(0);
  if (!failure)
    failure = b.create_spatial_varcoef(0);
  if (!failure)
    failure = b.create_geospline(0);
  if (!failure)
    failure = b.create_geospline_varcoeff(0);
  if (!failure)
    failure = b.create_varcoeffpspline(0);
  if (!failure)
    failure = b.create_random(0);
  if (!failure)
    failure = b.create_randomslope(0);
  if (!failure)
    failure = b.create_interactionspspline(0);
  if (!failure)
    failure = b.create_varcoeffinteractionspspline(0);
  if (!failure)
    failure = b.create_kriging(0);
  if (!failure)
    failure = b.create_kriging_1dim(0);
  if (!failure)
    failure = b.create_kriging_varcoeff(0);
  if (!failure)
    failure = b.create_geokriging(0);
  if (!failure)
    failure = b.create_geokriging_varcoeff(0);

// Setup for plotting results

  if(!failure)
    {
    unsigned i,j,k;
    if(b.family.getvalue()=="multinomial")
      {
      b.nrterms = b.cats.rows() * b.fullcond.size();
      b.needscat = vector<bool>(b.nrterms,true);
      b.fullcondnr = vector<unsigned>(b.nrterms,0);
      b.catnr = statmatrix<double>(b.nrterms,1,0);
      for(j=0; j<b.cats.rows(); j++)
        {
        for(i=0; i<b.fullcond.size(); i++)
          {
          b.fullcondnr[j*b.fullcond.size()+i] = i;
          b.catnr(j*b.fullcond.size()+i,0) = b.cats(j,0);
          }
        }
      }
    else if(b.family.getvalue()=="multinomialcatsp")
      {
      b.nrterms=1;
      for(i=1; i<b.fullcond.size(); i++)
        {
        if(b.fullcond[i]->get_catspecific())
          {
          b.nrterms++;
          }
        else
          {
          b.nrterms += b.cats.rows();
          }
        }
      b.needscat = vector<bool>(b.nrterms,false);
      b.fullcondnr = vector<unsigned>(b.nrterms,0);
      b.catnr = statmatrix<double>(b.nrterms,1,0);
      k=1;
      for(i=1; i<b.fullcond.size(); i++)
        {
        if(b.fullcond[i]->get_catspecific())
          {
          b.fullcondnr[k] = i;
          b.needscat[k] = false;
          k++;
          }
        else
          {
          for(j=0; j<b.cats.rows(); j++)
            {
            b.fullcondnr[k] = i;
            b.needscat[k] = true;
            b.catnr(k,0) = b.cats(j,0);
            k++;
            }
          }
        }
      }
    else if(b.family.getvalue()=="cumlogit" || b.family.getvalue()=="cumprobit" ||
            b.family.getvalue()=="seqlogit" || b.family.getvalue()=="seqprobit")
      {
      b.nrterms=1;
      for(i=1; i<b.fullcond.size(); i++)
        {
        if(b.fullcond[i]->get_catspecific())
          {
          b.nrterms += b.cats.rows();
          }
        else
          {
          b.nrterms++;
          }
        }
      b.needscat = vector<bool>(b.nrterms,false);
      b.fullcondnr = vector<unsigned>(b.nrterms,0);
      b.catnr = statmatrix<double>(b.nrterms,1,0);
      k=1;
      for(i=1; i<b.fullcond.size(); i++)
        {
        if(b.fullcond[i]->get_catspecific())
          {
          for(j=0; j<b.cats.rows(); j++)
            {
            b.fullcondnr[k] = i;
            b.needscat[k] = true;
            b.catnr(k,0) = b.cats(j,0);
            k++;
            }
          }
        else
          {
          b.fullcondnr[k] = i;
          b.needscat[k] = false;
          k++;
          }
        }
      }
    else
      {
      b.nrterms = b.fullcond.size();
      b.needscat = vector<bool>(b.nrterms,false);
      b.fullcondnr = vector<unsigned>(b.nrterms,0);
      for(i=0; i<b.nrterms; i++)
        {
        b.fullcondnr[i] = i;
        }
      }
    }

  if (!failure)
    {
    header= "remlreg object " + b.name.to_bstr() + ": reml procedure" ;

// Nominale Modelle
    if (b.family.getvalue()=="multinomial")
      {
      b.RE_M = remlest_multinomial(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),
      b.maxchange.getvalue(), b.maxvar.getvalue(),b.cats,weight,
      b.fisher.getvalue(),b.logout);
      if (b.fullcond.size() == 1)    // fixed effects only
        failure = b.RE_M.estimate_glm(response,offset,weight);
      else
        failure = b.RE_M.estimate(response,offset,weight);
      }
// Nominale Modelle mit kategorienspezifischen Kovariablen
    else if (b.family.getvalue()=="multinomialcatsp")
      {
      b.RE_M_catsp = remlest_multinomial_catsp(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond, response, b.family.getvalue(), b.outfile.getvalue(),
      b.maxit.getvalue(), b.lowerlim.getvalue(), b.eps.getvalue(),
      b.maxchange.getvalue(), b.maxvar.getvalue(), b.cats, weight,
      b.fisher.getvalue(), b.logout);
      if (b.fullcond.size() == 1)    // fixed effects only
        failure = b.RE_M_catsp.estimate_glm(response,offset,weight,b.naind);
      else
        failure = b.RE_M_catsp.estimate(response,offset,weight,b.naind);
      }
// Ordinale Modelle
    else if (b.family.getvalue()=="cumlogit" ||
        b.family.getvalue()=="cumprobit" ||
        b.family.getvalue()=="seqlogit" ||
        b.family.getvalue()=="seqprobit")
      {
      b.RE_O = remlest_ordinal(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),
      b.maxchange.getvalue(), b.maxvar.getvalue(),b.cats,weight,
      b.fisher.getvalue(),b.logout);
      if(b.RE_O.get_catspec())
        {
        if (b.fullcond.size() == 1)    // fixed effects only
          failure = b.RE_O.estimate_glm2(response,offset,weight);
        else
          failure = b.RE_O.estimate2(response,offset,weight);
        }
      else
        {
        if (b.fullcond.size() == 1)    // fixed effects only
          failure = b.RE_O.estimate_glm(response,offset,weight);
        else
          failure = b.RE_O.estimate(response,offset,weight);
        }
      }
// Univariate Modelle ohne Dispersionsparameter
    else if (b.family.getvalue()=="binomial" ||
        b.family.getvalue()=="binomialprobit" ||
        b.family.getvalue()=="binomialcomploglog" ||
        b.family.getvalue()=="poisson")
      {
      dispers=false;
      b.RE = remlest(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,dispers,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),
      b.maxchange.getvalue(), b.maxvar.getvalue(), b.fisher.getvalue(),
      b.constlambda.getvalue(),b.constscale.getvalue(),b.logout);
      if (b.fullcond.size() == 1)    // fixed effects only
        failure = b.RE.estimate_glm(response,offset,weight);
      else
        failure = b.RE.estimate(response,offset,weight);
      }
// Cox-Modell (alte Implementierung: Rechtszensierung + zeitvariierende Effekte
    else if (b.family.getvalue()=="coxold")
      {
      dispers=false;
      b.RE = remlest(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,dispers,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),
      b.maxchange.getvalue(), b.maxvar.getvalue(), b.fisher.getvalue(),
      b.constlambda.getvalue(),b.constscale.getvalue(),b.logout);
      failure = b.RE.estimate_survival(response,offset,weight);
      }
// Cox-Modell mit Intervallzensierung (ohne zeitvariierende Effekte)
    else if (b.family.getvalue()=="coxinterval")
      {
      dispers=false;
      b.RE = remlest(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,dispers,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),
      b.maxchange.getvalue(), b.maxvar.getvalue(), b.fisher.getvalue(),
      b.constlambda.getvalue(),b.constscale.getvalue(),b.logout);
      failure = b.RE.estimate_survival_interval(response,offset,weight);
      }
// Cox-Modell mit Intervallzensierung & Linkstrunkierung  & zeitvariierenden Effekten
    else if (b.family.getvalue()=="cox")
      {
      dispers=false;
      b.RE = remlest(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond, response, dispers, b.family.getvalue(), b.outfile.getvalue(),
      b.maxit.getvalue(), b.lowerlim.getvalue(), b.eps.getvalue(),
      b.maxchange.getvalue(), b.maxvar.getvalue(), b.fisher.getvalue(),
      b.constlambda.getvalue(),b.constscale.getvalue(),b.logout);
      failure = b.RE.estimate_survival_interval2(response,offset,weight,
                                                 b.aiccontrol.getvalue());
      }
// AFT-models with smoothed error distribution
/*    else if (b.family.getvalue()=="aft")
      {
      dispers=false;
      b.RE = remlest(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,dispers,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),
      b.maxchange.getvalue(),b.logout);
      if (b.fullcond.size() == 1)    // fixed effects only
        failure = b.RE.estimate_aft_glm(response,offset,weight);
      else
        failure = b.RE.estimate_aft(response,offset,weight);
      }*/
// Univariate Modelle mit Dispersionsparameter
    else
      {
      dispers=true;
      b.RE = remlest(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond,response,dispers,b.family.getvalue(),b.outfile.getvalue(),
      b.maxit.getvalue(),b.lowerlim.getvalue(),b.eps.getvalue(),
      b.maxchange.getvalue(), b.maxvar.getvalue(), b.fisher.getvalue(),
      b.constlambda.getvalue(),b.constscale.getvalue(),
      b.logout);
      if (b.fullcond.size() == 1)    // fixed effects only
        failure = b.RE.estimate_glm_dispers(response,offset,weight);
      else
        failure = b.RE.estimate_dispers(response,offset,weight);
      }
    }

// Compute graphics
#if defined(JAVA_OUTPUT_WINDOW)
  if(!failure)
    {
    unsigned j;
    for(j=0; j<b.nrterms; j++)
      {
      MCMC::plotstyles plst = b.fullcond[b.fullcondnr[j]]->get_plotstyle();
      if(plst != MCMC::noplot)
        {
        vector<ST::string> varnames = b.fullcond[b.fullcondnr[j]]->get_datanames();
        ST::string xvar = varnames[0];
        ST::string effect = xvar;
        if(varnames.size()>1)
          {
          effect = varnames[1] + "*" + effect;
          }
        ST::string pathresult = b.fullcond[b.fullcondnr[j]]->get_pathresult();
        if(b.needscat[j])
          {
          pathresult = pathresult.insert_after_string(ST::doubletostring(b.catnr(j,0),6)+"_","_f_");
          }
        ST::string pathps = pathresult.substr(0, pathresult.length()-4);
        if(plst == MCMC::plotnonp)
          {
          b.newcommands.push_back(b.name + ".plotnonp " + ST::inttostring(j)
          + ", title = \"Effect of " + effect +"\" xlab = " + xvar
          + " ylab = \" \" outfile = " + pathps + ".ps replace");
          }
        else if(plst==MCMC::drawmap)
          {
          double u = b.fullcond[b.fullcondnr[j]]->get_level1();
          double o = b.fullcond[b.fullcondnr[j]]->get_level2();
          ST::string u_str = ST::doubletostring(u,0);
          ST::string o_str = ST::doubletostring(o,0);
          b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", color outfile = " + pathps + "_pmode.ps replace");
          b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", plotvar = pcat" + u_str + " nolegend  pcat outfile = " + pathps
            + "_pcat" + u_str + ".ps replace");
          b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", plotvar = pcat" + o_str + " nolegend  pcat outfile = " + pathps
            + "_pcat" + o_str + ".ps replace");
          }
        }
      }

/*    for(unsigned j=0;j<b.fullcond.size();j++)
      {
      MCMC::plotstyles plst = b.fullcond[j]->get_plotstyle();
      if(plst != MCMC::noplot)
        {
        vector<ST::string> varnames = b.fullcond[j]->get_datanames();
        ST::string xvar = varnames[0];
        ST::string effect = xvar;
        if(varnames.size()>1)
          {
          effect = varnames[1] + "*" + effect;
          }
        if(b.family.getvalue()=="multinomial")
          {
          for(unsigned i=0; i<b.cats.rows(); i++)
            {
            ST::string pathresult = b.fullcond[j]->get_pathresult();
            pathresult = pathresult.insert_after_string(ST::doubletostring(b.cats(i,0),6)+"_","_f_");
            ST::string pathps = pathresult.substr(0, pathresult.length()-4);
            if(plst == MCMC::plotnonp)
              {
              b.newcommands.push_back(b.name + ".plotnonp " + ST::inttostring(i*b.fullcond.size()+j)
              + ", title = \"Effect of " + effect +"\" xlab = " + xvar
              + " ylab = \" \" outfile = " + pathps + ".ps replace");
              }
            else if(plst==MCMC::drawmap)
              {
              double u = b.fullcond[j]->get_level1();
              double o = b.fullcond[j]->get_level2();
              ST::string u_str = ST::doubletostring(u,0);
              ST::string o_str = ST::doubletostring(o,0);
              b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(i*b.fullcond.size()+j)
              + ", color outfile = " + pathps + "_pmode.ps replace");
              b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(i*b.fullcond.size()+j)
              + ", plotvar = pcat" + u_str + " nolegend  pcat outfile = " + pathps
              + "_pcat" + u_str + ".ps replace");
              b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(i*b.fullcond.size()+j)
              + ", plotvar = pcat" + o_str + " nolegend  pcat outfile = " + pathps
              + "_pcat" + o_str + ".ps replace");
              }
            }
          }
        else
          {
          ST::string pathresult = b.fullcond[j]->get_pathresult();
          ST::string pathps = pathresult.substr(0, pathresult.length()-4);
          if(plst == MCMC::plotnonp)
            {
            b.newcommands.push_back(b.name + ".plotnonp " + ST::inttostring(j)
            + ", title = \"Effect of " + effect +"\" xlab = " + xvar
            + " ylab = \" \" outfile = " + pathps + ".ps replace");
            }
          else if(plst==MCMC::drawmap)
            {
            double u = b.fullcond[j]->get_level1();
            double o = b.fullcond[j]->get_level2();
            ST::string u_str = ST::doubletostring(u,0);
            ST::string o_str = ST::doubletostring(o,0);
            b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", color outfile = " + pathps + "_pmode.ps replace");
            b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", plotvar = pcat" + u_str + " nolegend  pcat outfile = " + pathps
            + "_pcat" + u_str + ".ps replace");
            b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", plotvar = pcat" + o_str + " nolegend  pcat outfile = " + pathps
            + "_pcat" + o_str + ".ps replace");
            }
          }
        }
      }*/

    b.newcommands.push_back(b.name + ".texsummary");
    }
#endif

// Produce batch-file for graphics and model summary in tex
  if(!failure)
    {
    ST::string path = b.outfile.getvalue() + "_graphics.prg";
    ST::string path2 = b.outfile.getvalue() + "_model_summary.tex";
    ST::string path3 = b.outfile.getvalue() +  "_r.R";

    if(b.family.getvalue()=="multinomial")
      {
      b.RE_M.make_graphics(header,path,path2,path3,
                         b.modreg.getModelVarnamesAsVector()[0].to_bstr());
      }
    else if(b.family.getvalue()=="multinomialcatsp")
      {
      b.RE_M_catsp.make_graphics(header,path,path2,path3,
                         b.modreg.getModelVarnamesAsVector()[0].to_bstr());
      }
    else if(b.family.getvalue()=="cumlogit" || b.family.getvalue()=="cumprobit" ||
            b.family.getvalue()=="seqlogit" || b.family.getvalue()=="seqprobit")
      {
      b.RE_O.make_graphics(header,path,path2,path3,
                         b.modreg.getModelVarnamesAsVector()[0].to_bstr());
      }
    else
      {
      b.RE.make_graphics(header,path,path2,path3,
                         b.modreg.getModelVarnamesAsVector()[0].to_bstr(),
                         dispers);
      }
    }

  if (!failure)
    {
    b.resultsyesno = true;
    }
  else
    {
    b.describetext.erase(b.describetext.begin(),b.describetext.end());
    b.describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
    b.resultsyesno = false;
    }
  }

// -----------------------------------------------------------------------------
// ----------------------------- drawmaprun ------------------------------------
// -----------------------------------------------------------------------------

void drawmaprun(remlreg & b)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  b.outerror("ERROR: method drawmap is not available in this version\n");
#elif defined(__BUILDING_GNU)
  b.outerror("ERROR: method drawmap is not available in this version\n");
#elif defined(JAVA_OUTPUT_WINDOW)
  bool error = false;

  vector<ST::string> varnames = b.mdrawmap.getModelVarnamesAsVector();
  if (varnames.size() != 1)
    {
    b.outerror("ERROR: syntax error for method drawmap\n");
    error = true;
    }

  long nr;
  if (varnames[0].strtolong(nr) != 0)
    {
    b.outerror("ERROR: syntax error for method drawmap\n");
    error = true;
    }

  if(nr < 0 || nr>=b.nrterms)
    {
    b.outerror("ERROR: syntax error for method drawmap\n");
    error = true;
    }
/*  unsigned catnum;
  if(b.ismultinomial)
    {
    catnum = nr / b.fullcond.size();
    nr = nr % b.fullcond.size();
    if(catnum >= b.cats.rows())
      {
      b.outerror("ERROR: syntax error for method plotnonp\n");
      error = true;
      }
    }

  if (nr < 0 || nr >= b.fullcond.size())
    {
    b.outerror("ERROR: syntax error for method drawmap\n");
    error = true;
    }*/

  if (error == false)
    {
//    if (b.fullcond[nr]->get_plotstyle() != MCMC::drawmap)
    if (b.fullcond[b.fullcondnr[nr]]->get_plotstyle() != MCMC::drawmap
                           && b.fullcond[b.fullcondnr[nr]]->get_plotstyle() != MCMC::drawmapgraph)
      {
      error = true;
      b.outerror("ERROR: results cannot be visualized with method drawmap\n");
      }
    else if (b.fullcond[b.fullcondnr[nr]]->get_plotstyle() == MCMC::drawmapgraph)
      {
      error = true;
      b.outerror("ERROR: boundaries of the regions are not available from the graph-file \n");
      }
    }

  if (error==false)
    {
/*    ST::string path = b.fullcond[nr]->get_pathresult();
    if(b.ismultinomial)
      {
      path = path.insert_after_string(ST::doubletostring(b.cats(catnum,0),6)+"_","_f_");
      }*/
    ST::string path = b.fullcond[b.fullcondnr[nr]]->get_pathresult();
    if(b.needscat[nr])
      {
      path = path.insert_after_string(ST::doubletostring(b.catnr(nr,0),6)+"_","_f_");
      }

    vector<ST::string> vnames;
    ifstream in(path.strtochar());
    ST::string h;
    ST::getline(in,10000,h);
    h = h.eatallcarriagereturns();
    vnames = 	h.strtoken(" ");

    ST::string graphname = "_" + b.name + "_graph";
    b.newcommands.push_back("graph " + graphname);

    ST::string datasetname = "_" + b.name + "_r0";
    b.newcommands.push_back("dataset " + datasetname);
    b.newcommands.push_back(datasetname + ".infile , nonote using " + path);

    ST::string plotvar;
      plotvar = b.plotvar.getvalue() + " " + vnames[1] + " ";

//    ST::string ot="map=" + b.fullcond[nr]->getinfo() + " ";
    ST::string ot="map=" + b.fullcond[b.fullcondnr[nr]]->getinfo() + " ";

    ot = ot + "nrcolors="+b.nrcolors.getValueAsString()+" ";
    ot = ot + "title=\""+b.title2.getvalue() + "\" ";
    if (b.outfile4.getvalue().length() > 0)
      ot = ot + "outfile=\""+b.outfile4.getvalue() + "\" ";
    if (b.nolegend.getvalue() == true)
      ot = ot + "nolegend ";
    if (b.color.getvalue() == true)
      ot = ot + "color ";
    if (b.swapcolors.getvalue() == true)
      ot = ot + "swapcolors ";
    if (b.replace.getvalue() == true)
      ot = ot + "replace ";
    if (b.lowerlimit.changed() == true)
      ot = ot + "lowerlimit="  + b.lowerlimit.getValueAsString() + " ";
    if (b.upperlimit.changed() == true)
      ot = ot + "upperlimit=" + b.upperlimit.getValueAsString() + " ";
    if (b.pcat.getvalue() == true)
      ot = ot + "pcat ";
    if (b.drawnames.getvalue() == true)
      ot = ot + "drawnames ";
    if (b.hclcolors.getvalue() == true)
      ot = ot + "hcl ";
    if (b.fontsize.changed() == true)
      ot = ot + "fontsize=" + b.fontsize.getValueAsString() + " ";
    if (b.titlescale.changed() == true)
      ot = ot + "titlesize=" + b.titlescale.getValueAsString() + " ";

    if (ot.length() == 0)
      b.newcommands.push_back(graphname + ".drawmap " + plotvar + " using " + datasetname);
    else
      b.newcommands.push_back(graphname + ".drawmap " + plotvar + "," + ot + " using "
      + datasetname);

    b.newcommands.push_back("drop " + graphname + " " + datasetname);
    }
#endif
  }

// -----------------------------------------------------------------------------
// -------------------------- plotnonprun --------------------------------------
// -----------------------------------------------------------------------------

void plotnonprun(remlreg & b)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  b.outerror("ERROR: method plotnonp is not available in this version\n");
#elif defined(__BUILDING_GNU)
  b.outerror("ERROR: method plotnonp is not available in this version\n");
#elif defined(JAVA_OUTPUT_WINDOW)
  bool error = false;

  vector<ST::string> varnames = b.mplotnonp.getModelVarnamesAsVector();
  if (varnames.size() != 1)
    {
    b.outerror("ERROR: syntax error for method plotnonp\n");
    error = true;
    }

  long nr;
  if (varnames[0].strtolong(nr) != 0)
    {
    b.outerror("ERROR: syntax error for method plotnonp\n");
    error = true;
    }

/*  unsigned catnum;
  if(b.ismultinomial)
    {
    catnum = nr / b.fullcond.size();
    nr = nr % b.fullcond.size();
    if(catnum >= b.cats.rows())
      {
      b.outerror("ERROR: syntax error for method plotnonp\n");
      error = true;
      }
    }

  if (nr < 0 || nr >= b.fullcond.size())
    {
    b.outerror("ERROR: syntax error for method plotnonp\n");
    error = true;
    }*/

  if (nr < 0 || nr >= b.nrterms)
    {
    b.outerror("ERROR: syntax error for method plotnonp\n");
    error = true;
    }

  if (error == false)
    {
//    if (b.fullcond[nr]->get_plotstyle() != MCMC::plotnonp)
    if (b.fullcond[b.fullcondnr[nr]]->get_plotstyle() != MCMC::plotnonp)
      {
      error = true;
      b.outerror("ERROR: results cannot be visualized with method plotnonp\n");
      }
    }

  if (error==false)
    {
/*    ST::string path = b.fullcond[nr]->get_pathresult();
    if(b.ismultinomial)
      {
      path = path.insert_after_string(ST::doubletostring(b.cats(catnum,0),6)+"_","_f_");
      }*/
    ST::string path = b.fullcond[b.fullcondnr[nr]]->get_pathresult();
    if(b.needscat[nr])
      {
      path = path.insert_after_string(ST::doubletostring(b.catnr(nr,0),6)+"_","_f_");
      }

    vector<ST::string> vnames;
    ifstream in(path.strtochar());
    ST::string h;
    ST::getline(in,10000,h);
    h = h.eatallcarriagereturns();
    vnames = 	h.strtoken(" ");

    ST::string graphname = "_" + b.name + "_graph";
    b.newcommands.push_back("graph " + graphname);

    ST::string datasetname = "_" + b.name + "_r0";
    b.newcommands.push_back("dataset " + datasetname);
    b.newcommands.push_back(datasetname + ".infile , nonote using " + path);

    ST::string plotvar;
    if (b.levels.getvalue()=="all")
      {
      plotvar = vnames[1] + " ";
      if (b.median.getvalue() == true)
        plotvar = plotvar + vnames[5] + " ";
      else
        plotvar = plotvar + vnames[2] + " ";
      plotvar = plotvar + vnames[3] + " " +
                         vnames[4] + " " +
                         vnames[6] + " " +
                         vnames[7] + " ";
      }
    else if (b.levels.getvalue()=="none")
      {
      if (b.median.getvalue() == true)
        plotvar = vnames[1] + " " + vnames[5];
      else
        plotvar = vnames[1] + " " + vnames[2];

      }
    else if (b.levels.getvalue()=="1")
      {
      if (b.median.getvalue() == true)
        plotvar = vnames[1] + " " + vnames[5] + " " + vnames[3] + " " +
                  vnames[7] + " ";
      else
        plotvar = vnames[1] + " " + vnames[2] + " " + vnames[3] + " " +
                  vnames[7] + " ";
      }
    else
      {
      if (b.median.getvalue() == true)
        plotvar = vnames[1] + " " + vnames[5] + " " + vnames[4] + " " +
                  vnames[6] + " ";
      else
        plotvar = vnames[1] + " " + vnames[2] + " " + vnames[4] + " " +
                  vnames[6] + " ";
      }


    ST::string ot;
    ot = "xlab=\""+b.xlab.getvalue() + "\" ";
    ot = ot + "ylab=\""+b.ylab.getvalue() + "\" ";
    ot = ot + "title=\""+b.title.getvalue() + "\" ";
    if (b.outfile2.getvalue().length() > 0)
      ot = ot + "outfile=\""+b.outfile2.getvalue() + "\" ";
    ot = ot + "height="+b.height.getValueAsString() + " ";
    ot = ot + "width="+b.width.getValueAsString() + " ";
    if (b.replace2.getvalue() == true)
      ot = ot + " replace ";
    if (b.connect.changed() == true)
      ot = ot + "connect="+b.connect.getvalue() + " ";
    if (b.ylimbottom.changed() == true)
      ot = ot + "ylimbottom="+b.ylimbottom.getValueAsString() + " ";
    if (b.ylimtop.changed() == true)
      ot = ot + "ylimtop="+b.ylimtop.getValueAsString() + " ";
    if (b.ystep.changed() == true)
      ot = ot + "ystep="+b.ystep.getValueAsString() + " ";
    if (b.ystart.changed() == true)
      ot = ot + "ystart="+b.ystart.getValueAsString() + " ";
    if (b.xlimbottom.changed() == true)
      ot = ot + "xlimbottom="+b.xlimbottom.getValueAsString() + " ";
    if (b.xlimtop.changed() == true)
      ot = ot + "xlimtop="+b.xlimtop.getValueAsString() + " ";
    if (b.xstep.changed() == true)
      ot = ot + "xstep="+b.xstep.getValueAsString() + " ";
    if (b.xstart.changed() == true)
      ot = ot + "xstart="+b.xstart.getValueAsString() + " ";
    if (b.linewidth.changed() == true)
      ot = ot + "linewidth="+b.linewidth.getValueAsString() + " ";
    if (b.fontsize.changed() == true)
      ot = ot + "fontsize="+b.fontsize.getValueAsString() + " ";
    if (b.pointsize.changed() == true)
      ot = ot + "pointsize="+b.pointsize.getValueAsString() + " ";
    if (b.linecolor.changed() == true)
      ot = ot + "linecolor="+b.linecolor.getValueAsString() + " ";
    if (b.titlescale.changed() == true)
      ot = ot + "titlesize="+b.titlescale.getValueAsString() + " ";

    if (ot.length() == 0)
      b.newcommands.push_back(graphname + ".plot " + plotvar + " using " + datasetname);
    else
      b.newcommands.push_back(graphname + ".plot " + plotvar + "," + ot + " using "
      + datasetname);

    b.newcommands.push_back("drop " + graphname + " " + datasetname);
    }
#endif
  b.plotnonpoptions.setdefault();
  }

// -----------------------------------------------------------------------------
// -------------------------- texsummaryrun ------------------------------------
// -----------------------------------------------------------------------------

void texsummaryrun(remlreg & b)
  {

#if defined(BORLAND_OUTPUT_WINDOW)

  b.outerror("ERROR: method texsummary is not available in this version\n");

#elif defined(JAVA_OUTPUT_WINDOW)

  bool error = false;

  if (error==false)
    {

    ST::string path = b.outfile.getvalue();
    ST::string path2 = path;

    int i = path2.length()-1;
    bool gefunden = false;
    while(i>=0 && gefunden == false)
      {
      if(path2[i] == '\\' || path2[i]=='/')
        gefunden = true;
      path2 = path2.substr(0,i);
      i--;
      }

    ST::string helpbat = path2 + "_latexcommands.bat";
    ofstream outbat(helpbat.strtochar());
    outbat << "cd " << path2 << endl;
    outbat << path.substr(0,1) << ":" << endl;
    outbat << "latex " << path << "_model_summary.tex" << endl;
    outbat << "dvips " << path << "_model_summary.dvi" << endl;
    outbat.close();
    system(helpbat.strtochar());
    remove(helpbat.strtochar());
    }

#endif
  }

// -----------------------------------------------------------------------------
// ----------------------------- mremlrun --------------------------------------
// -----------------------------------------------------------------------------

void mremlrun(remlreg & b)
  {
  b.termsmult = b.modregmult.getterms();

  b.describetext.erase(b.describetext.begin(),b.describetext.end());
  b.describetext.push_back("LAST ESTIMATED MODEL: \n");
  b.describetext.push_back("\n");
  b.describetext.push_back(b.modregmult.getModelText());
  b.describetext.push_back("\n");

  b.clear();

  #if defined(JAVA_OUTPUT_WINDOW)
    bool failure = b.adminb_p->breakcommand();
  #else
    bool failure = false;
  #endif

  if ((b.family.getvalue() != "multistate"))
    {
    failure = true;
    b.out("ERROR: family " + b.family.getvalue() + " is not allowed for method mregress\n");
    }

  b.generaloptions = MCMCoptions(
  #if defined(JAVA_OUTPUT_WINDOW)
  b.adminb_p,
  #endif
  12000,2000,100,b.logout,b.level1.getvalue(),
                               b.level2.getvalue());

  ST::string header;

  unsigned i;

// Read design matrix, compute weights, etc.
  datamatrix weight;
  if (!failure)
    failure = b.create_data(weight);

// Define and check response
  datamatrix response;
  if(!failure)
    failure = b.create_response(response,weight);

// Compute offset
  datamatrix offset;
  if(!failure)
    failure = b.create_offset(offset);

// Extract state variable
  datamatrix state;
  if (!failure)
    {
    if(b.statepos==-1)
      {
      b.errormessages.push_back("  Variable state has to be specified as global option!");
      failure=true;
      }
    else
      {
      state = b.D.getCol(b.statepos);
      }
    }

// Compute model terms
  b.nrfullconds = vector<unsigned>(b.nrtransitions,0);
  bool glfrailty = false;
  if (!failure)
    {
    for (i=0;i<b.nrtransitions;i++)
      {
      b.terms = b.termsmult[i];
      if (!failure)
        failure = b.create_const(i);
      if( !failure)
        failure = b.create_baseline(i);
      if( !failure)
        failure = b.create_baseline_varcoeff(i);
      if (!failure)
        failure = b.create_nonprw1rw2(i);
      if (!failure)
        failure = b.create_nonprw1rw2_varcoef(i);
      if (!failure)
        failure = b.create_pspline(i);
      if (!failure)
       failure = b.create_nonpseason(i);
      if (!failure)
        failure = b.create_nonpseason_varcoef(i);
      if (!failure)
        failure = b.create_spatial(i);
      if (!failure)
        failure = b.create_spatial_varcoef(i);
      if (!failure)
         failure = b.create_geospline(i);
      if (!failure)
        failure = b.create_geospline_varcoeff(i);
      if (!failure)
        failure = b.create_varcoeffpspline(i);
      if (!failure)
        failure = b.create_random(i);
      if (!failure)
        failure = b.create_randomslope(i);
      if (!failure)
        failure = b.create_interactionspspline(i);
      if (!failure)
        failure = b.create_varcoeffinteractionspspline(i);
      if (!failure)
        failure = b.create_kriging(i);
      if (!failure)
        failure = b.create_kriging_1dim(i);
      if (!failure)
        failure = b.create_kriging_varcoeff(i);
      if (!failure)
        failure = b.create_geokriging(i);
      if (!failure)
        failure = b.create_geokriging_varcoeff(i);

      b.nrfullconds[i] = b.fullcond.size();
      }

    for(i=b.nrtransitions-1; i>0; i--)
      {
      b.nrfullconds[i] = b.nrfullconds.at(i) - b.nrfullconds.at(i-1);
      }

    // Add global frailty term (if specified)
    if(!failure)
      {
      if(b.gfrailtypos>-1)
        {
        glfrailty=true;
        ST::string pathnonp, pathres, title;
        double lambda = 1000;
        double startlambda = b.gflambdastart.getvalue();
        b.make_paths(0,pathnonp,pathres,title,b.globalfrailty.getvalue(),"",
                 "_random.raw","_random.res","_random");
        b.fcrandom.push_back(FULLCOND_random(&b.generaloptions,
                         b.D.getCol(b.gfrailtypos),title,pathnonp,
                         pathres,lambda,startlambda,false));
        b.fcrandom[b.fcrandom.size()-1].init_name(b.globalfrailty.getvalue());
        b.fcrandom[b.fcrandom.size()-1].set_fcnumber(b.fullcond.size());
        b.fullcond.push_back(&b.fcrandom[b.fcrandom.size()-1]);
        }
      }
    }

// Setup for plotting results

  if(!failure)
    {
    unsigned i;
    b.nrterms = b.fullcond.size();
    b.needscat = vector<bool>(b.nrterms,false);
    b.fullcondnr = vector<unsigned>(b.nrterms,0);
    for(i=0; i<b.nrterms; i++)
      {
      b.fullcondnr[i] = i;
      }
    }


  if (!failure)
    {
    header= "remlreg object " + b.name.to_bstr() + ": reml procedure" ;
    if (b.family.getvalue()=="multistate")
      {
      b.RE_MSM = remlest_multistate(
      #if defined(JAVA_OUTPUT_WINDOW)
      b.adminb_p,
      #endif
      b.fullcond, response, b.family.getvalue(), b.outfile.getvalue(),
      b.maxit.getvalue(), b.lowerlim.getvalue(), b.eps.getvalue(),
      b.maxchange.getvalue(), b.maxvar.getvalue(), b.constlambda.getvalue(),
      glfrailty, b.nrfullconds, weight, b.logout);

      failure = b.RE_MSM.estimate(response,offset,weight,state);
      }
    }

// Compute graphics
#if defined(JAVA_OUTPUT_WINDOW)
  if(!failure)
    {
    unsigned j;
    for(j=0; j<b.nrterms; j++)
      {
      MCMC::plotstyles plst = b.fullcond[b.fullcondnr[j]]->get_plotstyle();
      if(plst != MCMC::noplot)
        {
        vector<ST::string> varnames = b.fullcond[b.fullcondnr[j]]->get_datanames();
        ST::string xvar = varnames[0];
        ST::string effect = xvar;
        if(varnames.size()>1)
          {
          effect = varnames[1] + "*" + effect;
          }
        ST::string pathresult = b.fullcond[b.fullcondnr[j]]->get_pathresult();
        if(b.needscat[j])
          {
          pathresult = pathresult.insert_after_string(ST::doubletostring(b.catnr(j,0),6)+"_","_f_");
          }
        ST::string pathps = pathresult.substr(0, pathresult.length()-4);
        if(plst == MCMC::plotnonp)
          {
          b.newcommands.push_back(b.name + ".plotnonp " + ST::inttostring(j)
          + ", title = \"Effect of " + effect +"\" xlab = " + xvar
          + " ylab = \" \" outfile = " + pathps + ".ps replace");
          }
        else if(plst==MCMC::drawmap)
          {
          double u = b.fullcond[b.fullcondnr[j]]->get_level1();
          double o = b.fullcond[b.fullcondnr[j]]->get_level2();
          ST::string u_str = ST::doubletostring(u,0);
          ST::string o_str = ST::doubletostring(o,0);
          b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", color outfile = " + pathps + "_pmode.ps replace");
          b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", plotvar = pcat" + u_str + " nolegend  pcat outfile = " + pathps
            + "_pcat" + u_str + ".ps replace");
          b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
            + ", plotvar = pcat" + o_str + " nolegend  pcat outfile = " + pathps
            + "_pcat" + o_str + ".ps replace");
          }
        }
      }
    b.newcommands.push_back(b.name + ".texsummary");
    }
#endif

// Produce batch-file for graphics and model summary in tex
  if(!failure)
    {
    ST::string path = b.outfile.getvalue() + "_graphics.prg";
    ST::string path2 = b.outfile.getvalue() + "_model_summary.tex";
    ST::string path3 = b.outfile.getvalue() +  "_r.R";
    vector<ST::string> rnames;
    vector<unsigned> rescol  = b.modregmult.getresponsecol();
    for(i=0; i<b.nrtransitions; i++)
      {
      rnames.push_back(b.modregmult.getModelVarnamesAsVector()[rescol[i]].to_bstr());
      }

    b.RE_MSM.make_graphics(header,path,path2,path3,rnames);
    }

  if (!failure)
    {
    b.resultsyesno = true;
    }
  else
    {
    b.describetext.erase(b.describetext.begin(),b.describetext.end());
    b.describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
    b.resultsyesno = false;
    }
  }


#if defined(BORLAND_OUTPUT_WINDOW)
#pragma package(smart_init)
#endif
























