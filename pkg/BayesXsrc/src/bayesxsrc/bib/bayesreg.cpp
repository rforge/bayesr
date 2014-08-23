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

#include"bayesreg.h"
#include"dataobj.h"

// Vorschlag:
//#include<typeinfo.h>

#include<stddef.h>


//------------------------------------------------------------------------------
//------------- CLASS bayesreg: implementation of member functions -------------
//------------------------------------------------------------------------------

// Vorschlag:
//bool bayesreg::check_iwls(bool & iwls)
bool bayesreg::check_iwls(bool iwls,const unsigned & collinpred)
  {
  if ( ( (family.getvalue() == "binomial") && (iwls==true) ) ||
       ( (family.getvalue() == "poisson") && (iwls==true) ) ||
       ( (family.getvalue() == "gamma") && (iwls==true) ) ||
       ( (family.getvalue() == "vargaussian") && (iwls==true) ) ||
       ( (family.getvalue() == "nbinomial") && (iwls==true) ) ||
       ( (family.getvalue() == "zip") && (iwls==true) ) ||
       ((family.getvalue() == "gaussianh") && (collinpred==1) && (iwls==true)) ||
       ( (family.getvalue() == "multinomial") && (iwls==true) ) ||
       ( (family.getvalue() == "cox") && (iwls==true) ) ||
       ( (family.getvalue() == "multistate") && (iwls==true) )
      )
    return true;
  else
    return false;
  }

bool bayesreg::check_gaussian(const unsigned & collinpred)
  {

  if ( (family.getvalue() == "gaussian") ||
       (family.getvalue() == "gaussian_re") ||
     (family.getvalue() == "multgaussian") ||
     (family.getvalue() == "lognormal") ||
     (family.getvalue() == "binomialprobit") ||
     (family.getvalue() == "bernoullilogit") ||
     (family.getvalue() == "binomialtlink") ||
     (family.getvalue() == "multinomialprobit") ||
     ((family.getvalue() == "gaussianh") && (collinpred==0)) ||
     (family.getvalue() == "cumprobit") ||
     (family.getvalue() == "aft")
#if !defined (__BUILDING_THE_DLL)
     ||
     (family.getvalue() == "quantreg")
#endif
     )
     return true;
  else
    return false;
  }


bool  bayesreg::check_nongaussian(const unsigned & collinpred)
  {
  if ( (family.getvalue() == "binomial") || (family.getvalue() == "poisson") ||
       (family.getvalue() == "gamma") ||
       (family.getvalue() == "vargaussian") ||
       (family.getvalue() == "nbinomial") ||
       (family.getvalue() == "zip") ||
       ((family.getvalue() == "gaussianh") && (collinpred==1)) ||
       (family.getvalue() == "multinomial") || (family.getvalue() == "cox") ||
       (family.getvalue() == "multistate") )
     return true;
  else
    return false;
  }


void bayesreg::make_paths(unsigned  collinpred,ST::string & pathnonp,
                          ST::string & pathres,
                          ST::string & title,ST::string  varname1,
                          ST::string varname2,
                          ST::string  endingraw,ST::string endingres,
                          ST::string  endingtitle)
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



void bayesreg::create(void)
  {

  varianceest=false;
  RE_est=false;
  missingest=false;
  add_name="";

#if defined(__BUILDING_LINUX)
  ST::string h = defaultpath+"/output/"+name;
#else
  ST::string h = defaultpath+"\\output\\"+name;
#endif

  outfile = fileoption("outfile",h,false);

  iterationsprint = intoption("iterationsprint",1000,1,10000000);

  globaloptions.push_back(&outfile);
  globaloptions.push_back(&iterationsprint);

  distrstring = vector<ST::string>();

  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  // --------------------------- method regress --------------------------------

  offset = term_offset();
  fixedeffects = basic_termtype();
  nonprw1rw2 = term_autoreg();
  nonpseason = term_season();
  nonppspline = term_pspline();
  nonpspatial = term_spatial();
  nonpspatialxy = term_spatialxy();
  randomeff = term_random();
  hrandomeff = term_hrandom();
  randomeffslope = term_randomslope();
  mixtureeff = term_mixture();
  nonpvarcoeffpspline = term_varcoeff_pspline();
  nonpinteractpspline = term_interactpspline();
  nonpgeospline = term_geospline();
  nonpvarcoeffgeospline = term_varcoeff_geospline();
  nonpspatial_geokriging = term_geokriging();
  baseline = term_baseline();
  varcoeffbaseline = term_varcoeff_baseline();
  nonpvarcoeffmerror = term_varcoeff_merror();
  shrinkage = term_shrinkage();
  randomrw = term_random_autoreg();
  spatialrw = term_spatial_autoreg();
  randompspline = term_random_pspline();
  nigmix = term_nigmix();

  termtypes.push_back(&offset);
  termtypes.push_back(&fixedeffects);
  termtypes.push_back(&nonprw1rw2);
  termtypes.push_back(&nonpseason);
  termtypes.push_back(&nonppspline);
  termtypes.push_back(&nonpspatial);
  termtypes.push_back(&randomeff);
  termtypes.push_back(&hrandomeff);
  termtypes.push_back(&randomeffslope);
  termtypes.push_back(&mixtureeff);
  termtypes.push_back(&nonpvarcoeffpspline);
  termtypes.push_back(&nonpinteractpspline);
  termtypes.push_back(&nonpspatialxy);
  termtypes.push_back(&nonpgeospline);
  termtypes.push_back(&nonpvarcoeffgeospline);
  termtypes.push_back(&nonpspatial_geokriging);
  termtypes.push_back(&baseline);
  termtypes.push_back(&varcoeffbaseline);
  termtypes.push_back(&nonpvarcoeffmerror);
  termtypes.push_back(&shrinkage);
  termtypes.push_back(&nigmix);
  termtypes.push_back(&randomrw);
  termtypes.push_back(&spatialrw);
  termtypes.push_back(&randompspline);

  modreg = modelterm(&termtypes);

  udata = use();

  vector<ST::string> scalega;
  scalega.push_back("fixed");
  scalega.push_back("phi");
  scalega.push_back("random");
  scalegamma = stroption("scalegamma",scalega,"phi");
  gamvar = doubleoption("gammavar",0.001,0,1000);
  cit = intoption("cit",500,0,10000000);
  scalevalue = doubleoption("scale",1,0,1000000);
  constscale = simpleoption("constscale",false);

  knotsdef.push_back("equidistant");
  knotsdef.push_back("quantiles");
  knots = stroption("knots",knotsdef,"equidistant");

  modeonly = simpleoption("modeonly",false);
  noposteriormode = simpleoption("noposteriormode",false);

  blocksize = intoption("blocksize",20,5,500);

  setseed = intoption("setseed",-1,0,MAXINT);
  nographs = simpleoption("nographs",false);

  pseudocontourprob = simpleoption("pseudocontourprob",false);
  uniformprior = simpleoption("uniformprior",false);
  approx = simpleoption("approx",false);
  lengthstart = intoption("lengthstart",200,0,1000);

  iterations = intoption("iterations",52000,1,10000000);
  burnin = intoption("burnin",2000,0,500000);
  step = intoption("step",50,1,1000);
  level1 = doubleoption("level1",95,40,99);
  level2 = doubleoption("level2",80,40,99);
  maxint = intoption("maxint",150,0,20000);
// BEGIN: merror
  merror = intoption("merror",0,0,10);
// END: merror
  families.reserve(30);
  families.push_back("gaussian");
  families.push_back("gaussian_re");
  families.push_back("multgaussian");
  families.push_back("lognormal");
  families.push_back("binomial");
  families.push_back("binomialprobit");
  families.push_back("binomialtlink");
  families.push_back("bernoullilogit");
  families.push_back("poisson");
  families.push_back("gamma");
  families.push_back("vargaussian");
  families.push_back("multinomial");
  families.push_back("multinomialprobit");
  families.push_back("cumprobit");
  families.push_back("nbinomial");
  families.push_back("zip");
  families.push_back("cox");
  families.push_back("multistate");
  families.push_back("gaussianh");
  families.push_back("aft");
#if !defined (__BUILDING_THE_DLL)
  families.push_back("quantreg");
#endif
  family = stroption("family",families,"binomial");

#if !defined (__BUILDING_THE_DLL)
  quantile = doubleoption("quantile",0.5,0.001,0.999);
  mscheck = simpleoption("mscheck",false);
#endif

  aresp = doubleoption("aresp",0.001,-1.0,500);
  bresp = doubleoption("bresp",0.001,0.0,500);
  reference = doubleoption("reference",0,-10000,10000);


  vector<ST::string> dop;
  dop.push_back("nb");
  dop.push_back("poga");
  dop.push_back("poig");
  distopt = stroption("distopt",dop,"nb");

  vector<ST::string> pro;
  pro.push_back("uniform");
  pro.push_back("gamma");
  propopt = stroption("propopt",pro,"uniform");

  vector<ST::string> zipdop;
  zipdop.push_back("zip");
  zipdop.push_back("zinb");
  zipdop.push_back("zipga");
  zipdop.push_back("zipig");
  zipdistopt = stroption("zipdistopt",zipdop,"zip");

//vector<ST::string> zippro;
//zippro.push_back("uniform");
//zippro.push_back("gamma");
//zippropopt = stroption("propopt",zippro,"uniform");

  propvar = doubleoption("propvar",0.1,0,500);

  noconst = simpleoption("noconst",false);

  predict = simpleoption("predict",false);
  predictmu = simpleoption("predictmu",false);

// Vorschlag:
//  predictuntil=intoption("predictuntil",0,1,100000000000);
  predictuntil=intoption("predictuntil",0,1,1000000000);

  nutlink = intoption("nutlink",8,1,40);

  nu = intoption("nu",1,1,100);

  nosamples = simpleoption("nosamples",false);

  begin = stroption("begin");

  censoring = stroption("censoring");

  state = stroption("state");

  predictind = stroption("predictind");

  constlambda = simpleoption("constlambda",false);

  missingreg = simpleoption("missingreg",false);

  missingind = stroption("missingind");

  nosort = simpleoption("nosort",false);

  hierarchical = simpleoption("hierarchical",false);

  regressoptions.reserve(100);

  regressoptions.push_back(&modeonly);
  regressoptions.push_back(&noposteriormode);
  regressoptions.push_back(&setseed);
  regressoptions.push_back(&nographs);
  regressoptions.push_back(&pseudocontourprob);
  regressoptions.push_back(&uniformprior);
  regressoptions.push_back(&approx);
  regressoptions.push_back(&lengthstart);
  regressoptions.push_back(&iterations);
  regressoptions.push_back(&burnin);
  regressoptions.push_back(&level1);
  regressoptions.push_back(&level2);
  regressoptions.push_back(&step);
  regressoptions.push_back(&maxint);
  regressoptions.push_back(&family);
  regressoptions.push_back(&aresp);
  regressoptions.push_back(&bresp);
#if !defined (__BUILDING_THE_DLL)
  regressoptions.push_back(&quantile);
  regressoptions.push_back(&mscheck);
#endif
  regressoptions.push_back(&gamvar);
  regressoptions.push_back(&cit);
  regressoptions.push_back(&scalevalue);
  regressoptions.push_back(&scalegamma);
  regressoptions.push_back(&constscale);

  regressoptions.push_back(&merror);

  regressoptions.push_back(&nutlink);

  regressoptions.push_back(&reference);

  regressoptions.push_back(&knots);

  regressoptions.push_back(&predict);
  regressoptions.push_back(&predictmu);
  regressoptions.push_back(&predictuntil);

  regressoptions.push_back(&nu);

  regressoptions.push_back(&noconst);

  regressoptions.push_back(&propvar);
  regressoptions.push_back(&propopt);
  regressoptions.push_back(&distopt);
  regressoptions.push_back(&zipdistopt);
  regressoptions.push_back(&nosamples);

  regressoptions.push_back(&censoring);
  regressoptions.push_back(&begin);
  regressoptions.push_back(&state);

  regressoptions.push_back(&predictind);

  regressoptions.push_back(&constlambda);

  regressoptions.push_back(&missingreg);
  regressoptions.push_back(&missingind);

  regressoptions.push_back(&nosort);

  regressoptions.push_back(&hierarchical);

  regressoptions.push_back(&blocksize);

  // methods 0
  methods.push_back(command("regress",&modreg,&regressoptions,&udata,required,
			 optional,optional,optional,optional,required));

  functions[0] = regressrun;

  // --------------------------- method autocor --------------------------------

  maxlag = intoption("maxlag",250,1,500);

  autocorroptions.push_back(&maxlag);

  // methods 1
  methods.push_back(command("autocor",&ma,&autocorroptions,&ad,notallowed,
						  notallowed,notallowed,notallowed,optional,notallowed));

  functions[1] = autocorrrun;

  // -------------------------- method getsample -------------------------------

  // methods 2
  methods.push_back(command("getsample",&mgetsample,&getsampleoptions,
						  &usegetsample,notallowed,notallowed,notallowed,
                          notallowed,notallowed,notallowed));

  functions[2] = getsamplerun;

  //---------------------------- method mregress -------------------------------

  // methods 3
  modregmult = modeltermmult(&termtypes);
  methods.push_back(command("mregress",&modregmult,&regressoptions,&udata,
                   required,optional,optional,optional,optional,required));

  functions[3] = mregressrun;

  // --------------------------- method plotnonp -------------------------------

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
  title0 = stroption("title");
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
  plotnonpoptions.push_back(&title0);
  plotnonpoptions.push_back(&outfile2);
  plotnonpoptions.push_back(&replace2);
  plotnonpoptions.push_back(&linewidth);
  plotnonpoptions.push_back(&fontsize);
  plotnonpoptions.push_back(&pointsize);
  plotnonpoptions.push_back(&linecolor);
  plotnonpoptions.push_back(&titlescale);

  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  // methods 4
  methods.push_back(command("plotnonp",&mplotnonp,&plotnonpoptions,&uplotnonp,
                   required,notallowed,notallowed,notallowed,optional,
                   notallowed));

  functions[4] = plotnonprun;

  // --------------------------- method plotautocor ----------------------------

  uplotautocor = use();
  mplotautocor = modelStandard();

  maxlag2 = intoption("maxlag",250,1,500);
  outfile3 = stroption("outfile");
  meanonly = simpleoption("mean",false);
  replaceautocor = simpleoption("replace",false);

  plotautocoroptions.push_back(&maxlag2);
  plotautocoroptions.push_back(&outfile3);
  plotautocoroptions.push_back(&meanonly);
  plotautocoroptions.push_back(&replaceautocor);

  // methods 5
  methods.push_back(command("plotautocor",&mplotautocor,&plotautocoroptions,
                        &uplotautocor,notallowed,notallowed,notallowed,
                        notallowed,optional,notallowed));

  functions[5] = plotautocorrun;


  // -------------------------- method drawmaprun ------------------------------

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
  plotvar = stroption("plotvar","pmean");
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

  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  // methods 6
  methods.push_back(command("drawmap",&mdrawmap,&drawmapoptions,&udrawmap,
                   required,notallowed,notallowed,notallowed,optional,
                   notallowed));

  functions[6] = drawmaprun;

  // method outresults

  useoutresults = use();
  moutresults = modelStandard();

  vector<ST::string> trtypes;
  trtypes.reserve(10);
  trtypes.push_back("exp");
  trtypes.push_back("logit");
  trtypes.push_back("probit");
  trtypes.push_back("lognormal");
  trtypes.push_back("elasticity");
  trtypes.push_back("marginal");
  trtypes.push_back("oddsratio");

  transformtype = stroption("transformtype",trtypes,"exp");

  outresultsoptions.push_back(&transformtype);
  outresultsoptions.push_back(&nographs);

  // methods 7
  methods.push_back(command("outresults",&moutresults,
                    &outresultsoptions,&useoutresults,notallowed,
                    notallowed,notallowed,notallowed,optional,notallowed));

  functions[7] = outresultsrun;

  // -------------------------- method texsummary ------------------------------

  utexsummary = use();

  mtexsummary = modelStandard();

  // SYNTAX OF COMMANDS:
  // name [model] [weight varname] [by varname] [if expression]
  //      [, options] [using usingtext]

  // methods 8
  methods.push_back(command("texsummary",&mtexsummary,&texsummaryoptions,&utexsummary,
                   notallowed,notallowed,notallowed,notallowed,notallowed,
                   notallowed));

  functions[8] = texsummaryrun;

  // --------------------------- method hregress -------------------------------

  // methods 9
  methods.push_back(command("hregress",&modreg,&regressoptions,&udata,required,
			 optional,optional,optional,optional,required));

  functions[9] = hregressrun;

  }


void bayesreg::initpointers(void)
  {

  unsigned i;

//  unsigned r = distrstring.size();

  for(i=0;i<distrstring.size();i++)
    {
    if (distrstring[i] == "gaussian")
      distr.push_back(&distr_gaussian[distrposition[i]]);
    else if (distrstring[i] == "gaussian_RE")
      distr.push_back(&distr_gaussian_re[distrposition[i]]);
    else if (distrstring[i] == "gaussianh")
      distr.push_back(&distr_gaussianh);
    else if (distrstring[i] == "multgaussian")
      distr.push_back(&distr_multgaussian);
    else if (distrstring[i] == "lognormal")
      distr.push_back(&distr_lognormal);
    else if (distrstring[i] == "binomial")
      distr.push_back(&distr_binomial);
    else if (distrstring[i] == "poisson")
      distr.push_back(&distr_poisson);
    else if (distrstring[i] == "gamma")
      distr.push_back(&distr_gamma);
    else if (distrstring[i] == "vargaussian")
      distr.push_back(&distr_vargaussian);
    else if (distrstring[i] == "multinom")
      distr.push_back(&distr_multinom);
    else if (distrstring[i] == "multinom_latent")
      distr.push_back(&distr_multinom_latent);
    else if (distrstring[i] == "binomlat")
      distr.push_back(&distr_binomlat);
    else if (distrstring[i] == "cumlat3")
      distr.push_back(&distr_cumlat3);
    else if (distrstring[i] == "nbinomial")
      distr.push_back(&distr_nbinomial);
    else if (distrstring[i] == "zip")
      distr.push_back(&distr_zip);
    else if (distrstring[i] == "cox")
      distr.push_back(&distr_cox);
    else if (distrstring[i] == "multistate")
      distr.push_back(&distr_multistatemodel);
    else if (distrstring[i] == "aft")
      distr.push_back(&distr_aft);
#if !defined (__BUILDING_THE_DLL)
    else if (distrstring[i] == "quantreg")
      distr.push_back(&distr_quantreg);
#endif
    }


  for(i=0;i<fcbaseline.size();i++)
    fullcond.push_back(&fcbaseline[i]);

//  for(i=0;i<fcbaselineiwls.size();i++)
//    fullcond.push_back(&fcbaselineiwls[i]);

  for(i=0;i<fcmultibaseline.size();i++)
    fullcond.push_back(&fcmultibaseline[i]);

  for(i=0;i<fcpsplinesurfgaussian.size();i++)
    fullcond.push_back(&fcpsplinesurfgaussian[i]);

  for(i=0;i<normalconst.size();i++)
    fullcond.push_back(&normalconst[i]);

  for(i=0;i<normalconst_re.size();i++)
    fullcond.push_back(&normalconst_re[i]);

  for(i=0;i<nongaussianconst.size();i++)
    fullcond.push_back(&nongaussianconst[i]);

  for(i=0;i<nbinomialconst.size();i++)
    fullcond.push_back(&nbinomialconst[i]);

  for(i=0;i<fcnonp.size();i++)
    fullcond.push_back(&fcnonp[i]);

  for(i=0;i<fcnonpgaussian.size();i++)
    fullcond.push_back(&fcnonpgaussian[i]);

  for(i=0;i<fcpspline.size();i++)
    fullcond.push_back(&fcpspline[i]);

  for(i=0;i<fciwlspspline.size();i++)
    fullcond.push_back(&fciwlspspline[i]);

  for(i=0;i<fcpsplinegaussian.size();i++)
    fullcond.push_back(&fcpsplinegaussian[i]);

  for(i=0;i<fcpsplinesurf.size();i++)
    fullcond.push_back(&fcpsplinesurf[i]);

  for(i=0;i<fcrandom.size();i++)
    fullcond.push_back(&fcrandom[i]);

  for(i=0;i<fchrandom.size();i++)
    fullcond.push_back(&fchrandom[i]);

  for(i=0;i<fcrandomgaussian.size();i++)
    fullcond.push_back(&fcrandomgaussian[i]);

  for(i=0;i<fcmixture.size();i++)
    fullcond.push_back(&fcmixture[i]);

  for(i=0;i<fcvarnonp.size();i++)
    fullcond.push_back(&fcvarnonp[i]);

  for(i=0;i<fcvarnonpvec.size();i++)
    fullcond.push_back(&fcvarnonpvec[i]);

  for(i=0;i<fcvarnonpvecnigmix.size();i++)
    fullcond.push_back(&fcvarnonpvecnigmix[i]);

  for(i=0;i<fctvariance.size();i++)
    fullcond.push_back(&fctvariance[i]);

  for(i=0;i<fcadaptiv.size();i++)
    fullcond.push_back(&fcadaptiv[i]);

  for(i=0;i<fctvariance2dim.size();i++)
    fullcond.push_back(&fctvariance2dim[i]);

  for(i=0;i<fckriging.size();i++)
    fullcond.push_back(&fckriging[i]);

  for(i=0;i<fcmerror.size();i++)
    fullcond.push_back(&fcmerror[i]);

  for(i=0;i<normalshrinkage.size();i++)
    fullcond.push_back(&normalshrinkage[i]);

  for(i=0;i<nongaussianshrinkage.size();i++)
    fullcond.push_back(&nongaussianshrinkage[i]);

  }


void bayesreg::clear(void)
  {

  mind.erase(mind.begin(),mind.end());
  mind.reserve(5);

  outfiles.erase(outfiles.begin(),outfiles.end());
  outfiles.reserve(10);

  generaloptions.erase(generaloptions.begin(),generaloptions.end());
  generaloptions.reserve(10);

  distrstring.erase(distrstring.begin(),distrstring.end());
  distrstring.reserve(10);

  distrposition.erase(distrposition.begin(),distrposition.end());
  distrposition.reserve(10);

  distr.erase(distr.begin(),distr.end());
  distr.reserve(10);

  distr_gaussian.erase(distr_gaussian.begin(),distr_gaussian.end());
  distr_gaussian.reserve(5);

  distr_gaussian_re.erase(distr_gaussian_re.begin(),distr_gaussian_re.end());
  distr_gaussian_re.reserve(5);


  fullcond.erase(fullcond.begin(),fullcond.end());
  fullcond.reserve(200);

  normalconst.erase(normalconst.begin(),normalconst.end());
  normalconst.reserve(20);

  normalconst_re.erase(normalconst_re.begin(),normalconst_re.end());
  normalconst_re.reserve(20);

  nongaussianconst.erase(nongaussianconst.begin(),nongaussianconst.end());
  nongaussianconst.reserve(20);

  nbinomialconst.erase(nbinomialconst.begin(),nbinomialconst.end());
  nbinomialconst.reserve(20);

  Pmatrices.erase(Pmatrices.begin(),Pmatrices.end());
  Pmatrices.reserve(20);

  fcvarnonp.erase(fcvarnonp.begin(),fcvarnonp.end());
  fcvarnonp.reserve(40);

  fcvarnonpvec.erase(fcvarnonpvec.begin(),fcvarnonpvec.end());
  fcvarnonpvec.reserve(5);

  fcvarnonpvecnigmix.erase(fcvarnonpvecnigmix.begin(),fcvarnonpvecnigmix.end());
  fcvarnonpvecnigmix.reserve(5);

  fcnonp.erase(fcnonp.begin(),fcnonp.end());
  fcnonp.reserve(20);

  fcnonpgaussian.erase(fcnonpgaussian.begin(),fcnonpgaussian.end());
  fcnonpgaussian.reserve(20);

  fcpspline.erase(fcpspline.begin(),fcpspline.end());
  fcpspline.reserve(20);

  fciwlspspline.erase(fciwlspspline.begin(),fciwlspspline.end());
  fciwlspspline.reserve(20);

  fcpsplinegaussian.erase(fcpsplinegaussian.begin(),fcpsplinegaussian.end());
  fcpsplinegaussian.reserve(20);

  fcpsplinesurf.erase(fcpsplinesurf.begin(),fcpsplinesurf.end());
  fcpsplinesurf.reserve(20);

  fcpsplinesurfgaussian.erase(fcpsplinesurfgaussian.begin(),
                              fcpsplinesurfgaussian.end());
  fcpsplinesurfgaussian.reserve(20);


  fcrandom.erase(fcrandom.begin(),fcrandom.end());
  fcrandom.reserve(20);

  fchrandom.erase(fchrandom.begin(),fchrandom.end());
  fchrandom.reserve(20);

  fcrandomgaussian.erase(fcrandomgaussian.begin(),fcrandomgaussian.end());
  fcrandomgaussian.reserve(20);

  fcmixture.erase(fcmixture.begin(),fcmixture.end());
  fcmixture.reserve(20);

  fctvariance.erase(fctvariance.begin(),fctvariance.end());
  fctvariance.reserve(20);

  fctvariance2dim.erase(fctvariance2dim.begin(),fctvariance2dim.end());
  fctvariance2dim.reserve(20);

  fcadaptiv.erase(fcadaptiv.begin(),fcadaptiv.end());
  fcadaptiv.reserve(20);

  fcbaseline.erase(fcbaseline.begin(),fcbaseline.end());
  fcbaseline.reserve(20);

//  fcbaselineiwls.erase(fcbaselineiwls.begin(),fcbaselineiwls.end());
//  fcbaselineiwls.reserve(20);

  fcmultibaseline.erase(fcmultibaseline.begin(),fcmultibaseline.end());
  fcmultibaseline.reserve(20);

  fckriging.erase(fckriging.begin(),fckriging.end());
  fckriging.reserve(20);

  fcmerror.erase(fcmerror.begin(),fcmerror.end());
  fcmerror.reserve(2);

  normalshrinkage.erase(normalshrinkage.begin(),normalshrinkage.end());
  normalshrinkage.reserve(20);

  nongaussianshrinkage.erase(nongaussianshrinkage.begin(),nongaussianshrinkage.end());
  nongaussianshrinkage.reserve(20);

  fcmult.erase(fcmult.begin(),fcmult.end());
  fcmult.reserve(20);


  }


#if defined(JAVA_OUTPUT_WINDOW)
bayesreg::bayesreg(
administrator_basic * adb, administrator_pointer * adp,
const ST::string & n,ofstream * lo,istream * in,
						 ST::string p,vector<statobject*> * st)
						 : statobject(adb,n,"bayesreg",lo,in,p)
  {
  adminp_p = adp;
  statobj = st;
  create();
  resultsyesno = false;
  posteriormode = false;
  describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
  }
#else
bayesreg::bayesreg(const ST::string & n,ofstream * lo,istream * in,
						 ST::string p,vector<statobject*> * st)
						 : statobject(n,"bayesreg",lo,in,p)
  {
  statobj = st;
  create();
  resultsyesno = false;
  posteriormode = false;
  describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
  }
#endif


bayesreg::bayesreg(const bayesreg & b) : statobject(statobject(b))
  {
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = b.adminp_p;
  #endif
  statobj = b.statobj;
  D = b.D;
  distrstring = b.distrstring;
  modelvarnamesv = b.modelvarnamesv;
  simobj = b.simobj;
  distr_gaussian = b.distr_gaussian;
  distr_gaussian_re = b.distr_gaussian_re;
  distr_binomial = b.distr_binomial;
  distr_poisson = b.distr_poisson;
  distr_gamma = b.distr_gamma;
  distr_vargaussian = b.distr_vargaussian;
  distr_nbinomial = b.distr_nbinomial;
  distr_zip = b.distr_zip;
  distr_cox = b.distr_cox;
  distr_aft = b.distr_aft;
#if !defined (__BUILDING_THE_DLL)
  distr_quantreg = b.distr_quantreg;
#endif
  distr_gaussianh = b.distr_gaussianh;
  terms = b.terms;
  normalconst = b.normalconst;
  normalconst_re = b.normalconst_re;
  nongaussianconst = b.nongaussianconst;
  nbinomialconst = b.nbinomialconst;
  Pmatrices = b.Pmatrices;
  fcnonp = b.fcnonp;
  fcpspline = b.fcpspline;
  fckriging = b.fckriging;
  fcbaseline = b.fcbaseline;
//  fcbaselineiwls = b.fcbaselineiwls;
  fcmultibaseline = b.fcmultibaseline;
  normalshrinkage = b.normalshrinkage;
  nongaussianshrinkage = b.nongaussianshrinkage;
  fcmult = b.fcmult;
  resultsyesno = b.resultsyesno;
  posteriormode = b.posteriormode;
//  initpointers();
  }


const bayesreg & bayesreg::operator=(const bayesreg & b)
  {
  if (this == & b)
	 return *this;
  statobject::operator=(statobject(b));
  create();
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = b.adminp_p;
  #endif
  statobj = b.statobj;
  D = b.D;
  distrstring = b.distrstring;
  modelvarnamesv = b.modelvarnamesv;
  simobj = b.simobj;
  distr_gaussian = b.distr_gaussian;
  distr_gaussian_re = b.distr_gaussian_re;
  distr_binomial = b.distr_binomial;
  distr_poisson = b.distr_poisson;
  distr_gamma = b.distr_gamma;
  distr_vargaussian = b.distr_vargaussian;
  distr_nbinomial = b.distr_nbinomial;
  distr_zip = b.distr_zip;
  distr_cox = b.distr_cox;
  distr_aft = b.distr_aft;
#if !defined (__BUILDING_THE_DLL)
  distr_quantreg = b.distr_quantreg;
#endif
  distr_gaussianh = b.distr_gaussianh;
  terms = b.terms;
  normalconst = b.normalconst;
  normalconst_re = b.normalconst_re;
  nongaussianconst = b.nongaussianconst;
  nbinomialconst = b.nbinomialconst;
  Pmatrices = b.Pmatrices;
  fcnonp = b.fcnonp;
  fcpspline = b.fcpspline;
  fckriging = b.fckriging;
  fcbaseline = b.fcbaseline;
//  fcbaselineiwls = b.fcbaselineiwls;
  fcmultibaseline = b.fcmultibaseline;
  normalshrinkage = b.normalshrinkage;
  nongaussianshrinkage = b.nongaussianshrinkage;
  fcmult = b.fcmult;
  resultsyesno = b.resultsyesno;
  posteriormode = b.posteriormode;
//  initpointers();
  return *this;
  }


int bayesreg::parse(const ST::string & c)
  {

  int u = statobject::parse(c);

  int pos = statobject::parsecom(c,methods,globaloptions);

  if (pos >= 0)
	 (*functions[pos])(*this);

  return(pos);
  }


bool bayesreg::create_generaloptions(void)
  {

  if (iterations.getvalue()- burnin.getvalue() < 100)
    {
    outerror("ERROR: number of iterations must exceed number of burnin iterations about 100\n");
    return true;
    }

  if (step.getvalue() >= iterations.getvalue() - burnin.getvalue())
    {
    outerror("ERROR: thinning parameter too large\n");
    return true;
    }


  generaloptions.push_back(MCMCoptions(
  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p,
  #endif
  iterations.getvalue(),burnin.getvalue(),
                               step.getvalue(),logout,
                               level1.getvalue(),level2.getvalue()));

  describetext.push_back("ESTIMATION OPTIONS:\n");
  describetext.push_back("\n");
  describetext.push_back("Number of Iterations: "
                            + ST::inttostring(iterations.getvalue()) + "\n");
  describetext.push_back("Burnin: " + ST::inttostring(burnin.getvalue()) + "\n");
  describetext.push_back("Thinning parameter: " +
                            ST::inttostring(step.getvalue()) + "\n");

  generaloptions[generaloptions.size()-1].set_nrout(iterationsprint.getvalue());

  return false;

  }



bool bayesreg::create_distribution(void)
  {

  unsigned i;
  bool failure=false;


  //--------------------- reading dataset information --------------------------

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

  //------------------ end: reading dataset information ------------------------


  //---------------- reading data, creating designmatrices ---------------------

  ST::string rname;
  ST::string wn;
  ST::string ifexpression;
  ST::string predictindicator;
  ST::string missingindicator;
  unsigned weightpos=0;
  unsigned predictindpos=0;
  unsigned missingindpos=0;

  if (family.getvalue() =="multgaussian" ||
      family.getvalue() == "multistate"  ||
      family.getvalue() == "gaussianh"
     )
    {
    modelvarnamesv = modregmult.getModelVarnamesAsVector();
    wn = methods[5].get_weight_variable().to_bstr();
    if (predictind.changed())
      predictindicator = predictind.getvalue();
    if (wn.length() != 0)
      {
      modelvarnamesv.push_back(wn);
      weightpos = modelvarnamesv.size()-1;
      }
    if (predictindicator.length() != 0)
      {
      modelvarnamesv.push_back(predictindicator);
      predictindpos = modelvarnamesv.size()-1;
      }

    if (family.getvalue() == "multistate")
      {

      unsigned k;
      for(i=0;i<termsmult.size();i++)
        {
        for(k=0;k<termsmult[i].size();k++)
          {
          if ( baseline.checkvector(termsmult[i],k) == true )
            {
            ST::string beg = termsmult[i][k].options[13];
            if(beg.length() != 0)
              modelvarnamesv.push_back(beg);
            }
          }
        }

      } // end: if (family.getvalue() == "multistate")

    ifexpression = methods[5].getexpression();
    }
  else
    {
    modelvarnamesv = modreg.getModelVarnamesAsVector();

  // add variables for measurement error terms

  unsigned j;
  ST::string test;
  vector<ST::string> modelvarnamesvhelp;
  for(i=0; i<modelvarnamesv.size(); i++)
    {
    if(modelvarnamesv[i].length()>7)
      {
      test = modelvarnamesv[i].substr(modelvarnamesv[i].length()-7,7);
      if(test=="_merror")
        {
        test = modelvarnamesv[i].substr(0,modelvarnamesv[i].length()-7);
        for(j=0; j<merror.getvalue(); j++)
          {
          modelvarnamesvhelp.push_back(test + ST::inttostring(j+1));
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

    rname = modelvarnamesv[0].to_bstr();
    wn = methods[0].get_weight_variable().to_bstr();
    if (predictind.changed())
      predictindicator = predictind.getvalue();
    if (missingind.changed())
      missingindicator = missingind.getvalue();
    if (wn.length() != 0)
      {
      modelvarnamesv.push_back(wn);
      weightpos = modelvarnamesv.size()-1;
      }
    if (predictindicator.length() != 0)
      {
      modelvarnamesv.push_back(predictindicator);
      predictindpos = modelvarnamesv.size()-1;
      }
    if (missingindicator.length() != 0)
      {
      modelvarnamesv.push_back(missingind.getvalue());
      missingindpos = modelvarnamesv.size()-1;
      }

    ifexpression = methods[0].getexpression();
    }

  // Für Cox/Multistate model
  if(begin.getvalue() != "")
    {
    modelvarnamesv.push_back(begin.getvalue());
    begpos = modelvarnamesv.size()-1;
    }

  if(censoring.getvalue() != "")
    {
    modelvarnamesv.push_back(censoring.getvalue());
    censpos = modelvarnamesv.size()-1;
    }

  // Für Multistate model
  if(state.getvalue() != "")
    {
    modelvarnamesv.push_back(state.getvalue());
    statepos = modelvarnamesv.size()-1;
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
    } // end: if ((datap->allexisting(modelvarnamesv,notex)) == false)


  datap->makematrix(modelvarnamesv,D,ifexpression);

  errormessages = datap->geterrormessages();

  if (!errormessages.empty())
    return true;



  datamatrix offs(1,1,0);
  failure = create_offset(offs);


  datamatrix w;
  datamatrix pind;

  if (wn.length() > 0)
    {
    w = D.getCol(weightpos);
    }
  else
    w = datamatrix(1,1);

  if (predictindicator.length()>0)
    {
    pind = D.getCol(predictindpos);
    }
  else
    {
    pind = datamatrix(1,1);
    }

  if (missingindicator.length() != 0)
    {
    mind.push_back(D.getCol(missingindpos));
    }
  else
    {
    mind.push_back(datamatrix(1,1));
    }

  datamatrix beg;
  datamatrix statemat;

  if (begin.getvalue() == "")
    beg = datamatrix(D.rows(),1,0);
  else
    beg = D.getCol(begpos);

  if (state.getvalue() == "")
    statemat = datamatrix(1,1);
  else
    statemat = D.getCol(statepos);

  datamatrix cens;
  if(censoring.getvalue() == "")
    cens = datamatrix(D.rows(),1,1);
  else
    {
    cens = D.getCol(censpos);
    unsigned k=0;
    bool censerr = false;
    while(k<cens.rows() && !censerr)
      {
      if(cens(k,0)!=1 && cens(k,0)!=0)
        censerr=true;
      k++;
      }
    if(censerr)
      {
      outerror("ERROR: censoring indicator must be either zero or one\n");
      return true;
      }
    }


  describetext.push_back("Response distribution: "
                           + family.getvalue() + "\n");

  ST::string path = outfile.getvalue() + add_name + "_predictmean.raw";
  ST::string pathfull = outfile.getvalue() + add_name + "_predictmu.raw";
#if defined(__BUILDING_LINUX)
  ST::string pathfullsample = defaultpath + "/temp/" + name + add_name + "_predictmu.raw";
#else
  ST::string pathfullsample = defaultpath + "\\temp\\" + name + add_name + "_predictmu.raw";
#endif

  ST::string pathdev = outfile.getvalue() + add_name + "_deviance_sample.raw";


//---------------------------- Gaussian response -------------------------------
  if (family.getvalue() == "gaussian")
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
#if defined(__BUILDING_LINUX)
    ST::string path3 = defaultpath + "/temp/" + name + add_name + "_scale.raw";
#else
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";
#endif

    if (offs.rows() == 1)
      distr_gaussian.push_back(DISTRIBUTION_gaussian(aresp.getvalue(),bresp.getvalue(),
                                             &generaloptions[generaloptions.size()-1],
                                             D.getCol(0),path2,
                                             path3,w));
    else
      distr_gaussian.push_back(DISTRIBUTION_gaussian(offs,aresp.getvalue(),
                                             bresp.getvalue(),&generaloptions[generaloptions.size()-1],
                                             D.getCol(0),path2,path3,
                                             w));

    if (uniformprior.getvalue()==true)
      {
      distr_gaussian[distr_gaussian.size()-1].set_uniformprior();
      }

    if (constscale.getvalue()==true)
      {
      distr_gaussian[distr_gaussian.size()-1].set_constscale(scalevalue.getvalue());
      }


    distr_gaussian[distr_gaussian.size()-1].init_names(rname,wn);


    // prediction stuff

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_gaussian[distr_gaussian.size()-1].set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_gaussian[distr_gaussian.size()-1].set_predictfull(pathfullsample,pathfull,n);
      }

    if (pind.rows() > 1)
      distr_gaussian[distr_gaussian.size()-1].set_predictresponse(pind);

    // end: prediction stuff

    if (varianceest==true)
      {

      distr_gaussian[distr_gaussian.size()-1].set_variance(&distr_vargaussian);
      distr_vargaussian.set_gaussian(&distr_gaussian[distr_gaussian.size()-1]);

      }

    distr.push_back(&distr_gaussian[distr_gaussian.size()-1]);
    distrstring.push_back("gaussian");
    distrposition.push_back(distr_gaussian.size()-1);
    nrcategories = 1;
    }
//-------------------------- END: Gaussian response ----------------------------

//--------------------------------------- AFT  ---------------------------------
  else if (family.getvalue() == "aft")
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
#if defined(__BUILDING_LINUX)
    ST::string path3 = defaultpath + "/temp/" + name + add_name + "_scale.raw";
#else
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";
#endif

    if (offs.rows() == 1)
      distr_aft = DISTRIBUTION_AFT(aresp.getvalue(),bresp.getvalue(),
                                             &generaloptions[generaloptions.size()-1],
                                             D.getCol(0),cens,path2,
                                             path3,w);
    else
      distr_aft = DISTRIBUTION_AFT(offs,aresp.getvalue(),
                                             bresp.getvalue(),&generaloptions[generaloptions.size()-1],
                                             D.getCol(0),cens,path2,path3,
                                             w);

    if (uniformprior.getvalue()==true)
      {
      distr_aft.set_uniformprior();
      }

    if (constscale.getvalue()==true)
      {
      distr_aft.set_constscale(scalevalue.getvalue());
      }

    distr_aft.init_names(rname,wn);

    // prediction stuff

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_aft.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_aft.set_predictfull(pathfullsample,pathfull,n);
      }

    if (pind.rows() > 1)
      distr_aft.set_predictresponse(pind);

    // end: prediction stuff

    distr.push_back(&distr_aft);
    distrstring.push_back("AFT");
    distrposition.push_back(0);
    nrcategories = 1;
    }
//----------------------------------- END: AFT ---------------------------------

//----------------------------------- quantreg ---------------------------------
#if !defined (__BUILDING_THE_DLL)
  else if (family.getvalue() == "quantreg")
    {
    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
#if defined(__BUILDING_LINUX)
    ST::string path3 = defaultpath + "/temp/" + name + add_name + "_scale.raw";
#else
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";
#endif

    double quant = quantile.getvalue();

    if (offs.rows() == 1)
      distr_quantreg = DISTRIBUTION_QUANTREG(aresp.getvalue(),bresp.getvalue(),
                                             &generaloptions[generaloptions.size()-1],
                                             D.getCol(0),path2,
                                             path3,quant,w);
    else
      distr_quantreg = DISTRIBUTION_QUANTREG(offs,aresp.getvalue(),
                                             bresp.getvalue(),&generaloptions[generaloptions.size()-1],
                                             D.getCol(0),path2,path3,
                                             quant,w);

    distr_quantreg.init_names(rname,wn);
    distr_quantreg.set_constscale(scalevalue.getvalue());

    // prediction stuff

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_quantreg.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_quantreg.set_predictfull(pathfullsample,pathfull,n);
      }

    if (pind.rows() > 1)
      distr_quantreg.set_predictresponse(pind);

    // end: prediction stuff

    distr.push_back(&distr_quantreg);
    distrstring.push_back("QuantileReg");
    distrposition.push_back(0);
    nrcategories = 1;
    }
//-------------------------------- END: quantreg -------------------------------
#endif
//---------------------------- Gaussian response RE  ---------------------------
  else if (family.getvalue() == "gaussian_re")
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
#if defined(__BUILDING_LINUX)
    ST::string path3 = defaultpath + "/temp/" + name + add_name + "_scale.raw";
#else
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";
#endif

    if (offs.rows() == 1)
      distr_gaussian_re.push_back(DISTRIBUTION_gaussian_re(aresp.getvalue(),bresp.getvalue(),
                                             &generaloptions[generaloptions.size()-1],
                                             D.getCol(0),path2,
                                             path3,w));
    else
      distr_gaussian_re.push_back(DISTRIBUTION_gaussian_re(offs,aresp.getvalue(),
                                             bresp.getvalue(),&generaloptions[generaloptions.size()-1],
                                             D.getCol(0),path2,path3,
                                             w));

    if (uniformprior.getvalue()==true)
      {
      distr_gaussian_re[distr_gaussian_re.size()-1].set_uniformprior();
      }

    if (constscale.getvalue()==true)
      {
      distr_gaussian[distr_gaussian.size()-1].set_constscale(scalevalue.getvalue());
      }


    distr_gaussian_re[distr_gaussian_re.size()-1].init_names(rname,wn);

    // prediction stuff

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_gaussian_re[distr_gaussian_re.size()-1].set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_gaussian_re[distr_gaussian_re.size()-1].set_predictfull(pathfullsample,pathfull,n);
      }

    if (pind.rows() > 1)
      distr_gaussian_re[distr_gaussian_re.size()-1].set_predictresponse(pind);

    // end: prediction stuff

    distr.push_back(&distr_gaussian_re[distr_gaussian_re.size()-1]);
    distrstring.push_back("gaussian_re");
    distrposition.push_back(distr_gaussian_re.size()-1);
    nrcategories = 1;
    }
//-------------------------- END: Gaussian_RE response -------------------------

//---------------------- Gaussian with heteroscedastic variances ---------------
  else if (family.getvalue() == "gaussianh")
    {



    datamatrix res(D.rows(),2);
    unsigned k;
    for(k=0;k<D.rows();k++)
      {
      res(k,0) = D(k,0);
      res(k,1) = D(k,0);
      }



    distr_gaussianh = DISTRIBUTION_gaussianh(aresp.getvalue(),bresp.getvalue(),
                                       &generaloptions[generaloptions.size()-1],
                                             res,w);

    distr_gaussianh.init_names(rname,wn);


    // prediction stuff

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_gaussianh.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_gaussianh.set_predictfull(pathfullsample,pathfull,n);
      }

    if (pind.rows() > 1)
      distr_gaussianh.set_predictresponse(pind);

    // end: prediction stuff


    distr.push_back(&distr_gaussianh);
    distrstring.push_back("gaussianh");
    distrposition.push_back(0);
    nrcategories = 2;
    }
//-------------- END: Gaussian with heteroscedastic variances ------------------
  else if (family.getvalue() == "lognormal")
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
#if defined(__BUILDING_LINUX)
    ST::string path3 = defaultpath + "/temp/" + name + add_name + "_scale.raw";
#else
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";
#endif

    if (offs.rows() == 1)
      distr_lognormal = DISTRIBUTION_lognormal(aresp.getvalue(),bresp.getvalue(),
                                             &generaloptions[generaloptions.size()-1],
                                             D.getCol(0),path2,
                                             path3,w);
    else
      distr_lognormal = DISTRIBUTION_lognormal(offs,aresp.getvalue(),
                                             bresp.getvalue(),&generaloptions[generaloptions.size()-1],
                                             D.getCol(0),path2,path3,
                                             w);

    if (uniformprior.getvalue()==true)
      {
      distr_lognormal.set_uniformprior();
      }

    if (constscale.getvalue()==true)
      {
      distr_lognormal.set_constscale(scalevalue.getvalue());
      }


    distr_lognormal.init_names(rname,wn);


    // prediction stuff

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_lognormal.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_lognormal.set_predictfull(pathfullsample,pathfull,n);
      }

    if (pind.rows() > 1)
      distr_lognormal.set_predictresponse(pind);

    // end: prediction stuff

    if (varianceest==true)
      {

      distr_lognormal.set_variance(&distr_vargaussian);
      distr_vargaussian.set_gaussian(&distr_lognormal);

      }

    distr.push_back(&distr_lognormal);
    distrstring.push_back("lognormal");
    distrposition.push_back(0);
    nrcategories = 1;
    }
//-------------------------- END: lognormal ------------------------------------
  else if (family.getvalue() == "multgaussian")
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
#if defined(__BUILDING_LINUX)
    ST::string path3 = defaultpath + "/temp/" + name + add_name + "_scale.raw";
#else
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";
#endif

    vector<unsigned> rescol  = modregmult.getresponsecol();

    // prior for B

    datamatrix B(rescol.size(),rescol.size(),0);
    unsigned j;
    for(j=0;j<B.rows();j++)
      B(j,j) = bresp.getvalue();

    // Responsematrix

    unsigned c;

    datamatrix res(D.rows(),rescol.size());
    unsigned k;
    for(j=0;j<res.cols();j++)
      {
      c = rescol[j];
      for(k=0;k<res.rows();k++)
        res(k,j) = D(k,c);
      }

    if (offs.rows() == 1)
      distr_multgaussian = DISTRIBUTION_multgaussian(aresp.getvalue(),B,
                           &generaloptions[generaloptions.size()-1],res,
                           path2,path3);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_multgaussian.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_multgaussian.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_multgaussian);
    distrstring.push_back("multgaussian");
    distrposition.push_back(0);
    nrcategories = res.cols();

    }
//---------------------- binomial response, probit link ------------------------
  else if ( (family.getvalue() == "binomialprobit") ||
             (family.getvalue() == "binomialtlink") )
    {


    if (offs.rows() == 1)
      {
      if (family.getvalue() == "binomialtlink")
        distr_binomlat = DISTRIBUTION_binomial_latent(&generaloptions[generaloptions.size()-1],
                            D.getCol(0),w,true,nutlink.getvalue());
      else
        distr_binomlat = DISTRIBUTION_binomial_latent(&generaloptions[generaloptions.size()-1],
                            D.getCol(0),w,false);
      }
    else
      {

      if (family.getvalue() == "binomialtlink")
        distr_binomlat = DISTRIBUTION_binomial_latent(offs,&generaloptions[generaloptions.size()-1],
                            D.getCol(0),w,true,nutlink.getvalue());
      else
        distr_binomlat = DISTRIBUTION_binomial_latent(offs,&generaloptions[generaloptions.size()-1],
                            D.getCol(0),w,false);

      }

    if (distr_binomlat.geterrors().size() > 0)
      {
      outerror(distr_binomlat.geterrors());
      return true;
      }

    distr_binomlat.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_binomlat.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_binomlat.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_binomlat);
    distrstring.push_back("binomlat");
    distrposition.push_back(0);
    nrcategories = 1;
    }
//--------------------- end: binomial response, probit link --------------------

//----------------------- binomial response, logit link ------------------------
  else if (family.getvalue() == "binomial")
    {

    if (offs.rows() == 1)
      distr_binomial = DISTRIBUTION_binomial(
      &generaloptions[generaloptions.size()-1],D.getCol(0),w);
    else
      distr_binomial = DISTRIBUTION_binomial(offs,
      &generaloptions[generaloptions.size()-1],D.getCol(0),w);

    if (distr_binomial.geterrors().size() > 0)
      {
      outerror(distr_binomial.geterrors());
      return true;
      }

    distr_binomial.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_binomial.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_binomial.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_binomial);
    distrstring.push_back("binomial");
    distrposition.push_back(0);
    nrcategories = 1;
    }
//--------------------- END: binomial response, logit link ---------------------

//------------ bernoulli response, logit link, scale mixture of  normals -------
  else if (family.getvalue() == "bernoullilogit")
    {

    if (offs.rows() == 1)
      {
      distr_binomlogitlat = DISTRIBUTION_binomial_logit_latent(&generaloptions[generaloptions.size()-1],
                                  D.getCol(0),w,false);
      }
    else
      {
      distr_binomlogitlat = DISTRIBUTION_binomial_logit_latent(offs,
                            &generaloptions[generaloptions.size()-1],
                            D.getCol(0),w,false);
      }

    if (distr_binomlogitlat.geterrors().size() > 0)
      {
      outerror(distr_binomlogitlat.geterrors());
      return true;
      }

    distr_binomlogitlat.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_binomlogitlat.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
        n = D.rows();
      distr_binomlogitlat.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_binomlogitlat);
    distrstring.push_back("binomlogitlat");
    distrposition.push_back(0);
    nrcategories = 1;
    }
//--------------------- END: binomial response, logit link ---------------------

//------------------------------ gamma response --------------------------------
  else if (family.getvalue() == "gamma")
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
#if defined(__BUILDING_LINUX)
    ST::string path3 = defaultpath + "/temp/" + name + add_name + "_scale.raw";
#else
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";
#endif

    int st = cit.getvalue();
    double v1 = gamvar.getvalue();

    double sv = scalevalue.getvalue();

    ST::string type = scalegamma.getvalue();

    if (type=="fixed")
      {
      distr_gamma = DISTRIBUTION_gamma(sv,&generaloptions[generaloptions.size()-1],D.getCol(0),path2,
                                       path3,w);
      }
    else if (type=="phi")
      {
      distr_gamma = DISTRIBUTION_gamma(aresp.getvalue(),bresp.getvalue(),st,
               &generaloptions[generaloptions.size()-1],D.getCol(0),path2,path3,
                                     w);
      }
    else
      {
      distr_gamma = DISTRIBUTION_gamma(aresp.getvalue(),bresp.getvalue(),
                                       v1,st,&generaloptions[generaloptions.size()-1],D.getCol(0),
                                       path2,path3,w);
      }

    if (offs.rows() > 1)
      {
      distr_gamma.init_offset(offs);
      }


    distr_gamma.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_gamma.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
        n = D.rows();
      distr_gamma.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_gamma);
    distrstring.push_back("gamma");
    distrposition.push_back(0);
    nrcategories = 1;
    }
//----------------------------- END: gamma response ----------------------------

//----------------------------- vargaussian response ---------------------------
  else if (family.getvalue() == "vargaussian")
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
#if defined(__BUILDING_LINUX)
    ST::string path3 = defaultpath + "/temp/" + name + add_name + "_scale.raw";
#else
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";
#endif

    int st = cit.getvalue();
    double v1 = gamvar.getvalue();

    double sv = scalevalue.getvalue();

    ST::string type = scalegamma.getvalue();

    if (type=="fixed")
      {
      distr_vargaussian = DISTRIBUTION_vargaussian(
      sv,&generaloptions[generaloptions.size()-1],D.getCol(0),path2,
                                       path3,w);
      }
    else if (type=="phi")
      {
      distr_vargaussian = DISTRIBUTION_vargaussian(
               aresp.getvalue(),bresp.getvalue(),st,
               &generaloptions[generaloptions.size()-1],D.getCol(0),path2,path3,
                                     w);
      }
    else
      {
      distr_vargaussian = DISTRIBUTION_vargaussian(
                                       aresp.getvalue(),bresp.getvalue(),
                                       v1,st,
                                       &generaloptions[generaloptions.size()-1],
                                       D.getCol(0),path2,path3,w);
      }

    if (offs.rows() > 1)
      {
      distr_vargaussian.init_offset(offs);
      }


    distr_vargaussian.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_vargaussian.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
        n = D.rows();
      distr_vargaussian.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_vargaussian);
    distrstring.push_back("vargaussian");
    distrposition.push_back(0);
    nrcategories = 1;
    }
//------------------------- END: vargaussian response --------------------------

//-------------------- multinomial response, probit link -----------------------
  else if (family.getvalue() == "multinomialprobit")
    {

    if (wn.length() != 0)
      {
      outerror("ERROR: weight not allowed for multivariate response\n");
      return true;
      }

    if (offs.rows() > 1)
      {
      outerror("ERROR: offset not allowed for family=multinomialprobit\n");
      return true;
      }

    D.sort(0,D.rows()-1,0);

    if (reference.changed() == true)
      {
      bool existing=false;
      unsigned i;
      double * workD = D.getV();
      unsigned c = D.cols();
      i=0;
      double ref = reference.getvalue();
      while ( (i<D.rows()) && (existing == false) )
        {
        if (*workD == ref)
          existing = true;
        i++;
        if (i<D.rows())
          workD+=c;
        }

      if (existing == false)
        {
        outerror("ERROR: reference category is not existing\n");
        return true;
        }

      }


    distr_multinom_latent =
    DISTRIBUTION_multinomial_latent(
    &generaloptions[generaloptions.size()-1],D.getCol(0),
                                    reference.getvalue());
    distr_multinom_latent.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_multinom_latent.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_multinom_latent.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_multinom_latent);
    distrstring.push_back("multinom_latent");
    distrposition.push_back(0);
    nrcategories = distr_multinom_latent.get_nrcat();

    }
//------------------- END: multinomial response, probit link -------------------

//--------------------- multinomial response, logit link -----------------------
  else if (family.getvalue() == "multinomial")
    {

    if (wn.length() != 0)
      {
      outerror("ERROR: weight not allowed for multivariate response\n");
      return true;
      }

    if (offs.rows() > 1)
      {
      outerror("ERROR: offset not allowed for family=multinomial\n");
      return true;
      }

    statmatrix<int> index(D.rows(),1);
    index.indexinit();
    D.indexsort(index,0,D.rows()-1,0,0);

    unsigned nrcat = 1;
    unsigned refcat;
    double refvalue;

    unsigned i,j;
    vector<unsigned> beg;
    beg.push_back(0);
    bool existing = false;
    if (reference.getvalue() == D(index(0,0),0))
      {
      refvalue = D(index(0,0),0);
      refcat = 0;
      existing = true;
      }

    for (i=1;i<D.rows();i++)
      {
      if ( D(index(i,0),0) != D(index(i-1,0),0) )
        {
        beg.push_back(i);
        if (reference.getvalue() == D(index(i,0),0))
          {
          refcat = nrcat;
          refvalue = D(index(i,0),0);
          existing = true;
          }
        nrcat++;
        }

      }

    if (existing == false)
      {
      if (reference.changed() == true)
        {
        outerror("ERROR: reference category is not existing\n");
        return true;
        }
      else
        {
        refvalue = D(index(0,0),0);
        refcat = 0;
        }

      }

    if (nrcat == 1)
      {
      outerror("ERROR: response variable does not vary\n");
      return true;
      }

    if (nrcat > 10)
      {
      outerror("ERROR: too many values for the response variable\n");
      return true;
      }


    datamatrix response(D.rows(),nrcat-1,0);

    unsigned end;
    unsigned c=0;
    for(i=0;i<nrcat;i++)
      {
      if (i == nrcat-1)
        end = D.rows()-1;
      else
        end = beg[i+1]-1;

      if (i != refcat)
        {
        for (j=beg[i];j<=end;j++)
          response(index(j,0),c) = 1;
        c++;
        }

      }

    distr_multinom = DISTRIBUTION_multinom(
    &generaloptions[generaloptions.size()-1],response,refvalue,w);

    distr_multinom.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_multinom.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_multinom.set_predictfull(pathfullsample,pathfull,n);
      }


    distr.push_back(&distr_multinom);
    distrstring.push_back("multinom");
    distrposition.push_back(0);
    nrcategories = nrcat-1;
    }
//------------------- END: multinomial response, logit link -------------------

//------------------- cumulative threshold model, probit link ------------------
  else if (family.getvalue() == "cumprobit")
    {

    if (nosort.getvalue() == false)
      D.sort(0,D.rows()-1,0);

    if (wn.length() == 0)
      w = datamatrix(1,1);
    else
      w = D.getCol(D.cols()-1);

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
#if defined(__BUILDING_LINUX)
    ST::string path3 = defaultpath + "/temp/" + name + add_name + "_scale.raw";
#else
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";
#endif

    distr_cumlat3 = DISTRIBUTION_cumulative_latent3(
    &generaloptions[generaloptions.size()-1],D.getCol(0),
                    w,aresp.getvalue(),bresp.getvalue(),path2,path3);

    distr_cumlat3.init_names(rname);

    if (offs.rows() > 1)
      {
      distr_cumlat3.init_offset(offs);
      }


    if (predict.getvalue() == true || (predictmu.getvalue() == true) )
      distr_cumlat3.set_predict_cum(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_cumlat3.set_predictfull(pathfullsample,pathfull,n);
      }



    distr.push_back(&distr_cumlat3);
    distrstring.push_back("cumulat3");
    distrposition.push_back(0);
    nrcategories = 1;

    }
//----------------- END: cumulative threshold model, probit link ---------------

//--------------------------- poisson response ---------------------------------
  else if (family.getvalue() == "poisson")
    {

    if (offs.rows() == 1)
      distr_poisson = DISTRIBUTION_poisson(
      &generaloptions[generaloptions.size()-1],D.getCol(0),w);
    else
      distr_poisson = DISTRIBUTION_poisson(offs,
      &generaloptions[generaloptions.size()-1],D.getCol(0),w);
    distr_poisson.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_poisson.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_poisson.set_predictfull(pathfullsample,pathfull,n);
      }


    distr.push_back(&distr_poisson);
    distrstring.push_back("poisson");
    distrposition.push_back(0);
    nrcategories = 1;

    }
//------------------------ END: poisson response -------------------------------

//---------------------------- Cox model ---------------------------------------
  else if (family.getvalue() == "cox")
    {

    bool baselineexisting = false;

    unsigned j=0;
    for(i=0;i<terms.size();i++)
      {
      if ( baseline.checkvector(terms,i) == true )
        {
        j = terms[i].varnames[0].isinlist(modelvarnamesv);
        baselineexisting = true;
        }
      }

    if(baselineexisting == false)
      {
      outerror("ERROR: no baseline specified\n");
      return true;
      }

    if (offs.rows() == 1)
      distr_cox = DISTRIBUTION_coxmodel(
      &generaloptions[generaloptions.size()-1],
      D.getCol(0),D.getCol(j),beg,w);
    else
      distr_cox = DISTRIBUTION_coxmodel(
      offs,&generaloptions[generaloptions.size()-1],D.getCol(0),
      D.getCol(j),beg,w);

    distr_cox.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_cox.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_cox.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_cox);
    distrstring.push_back("cox");
    distrposition.push_back(0);
    nrcategories = distr_cox.get_nrcat();

    }
//--------------------- END: Cox model -----------------------------------------
//---------------------Multistate model-----------------------------------------
  else if (family.getvalue() == "multistate")
    {

    bool baselineexisting;

    unsigned j=0,k=0;
    for(i=0;i<termsmult.size();i++)
      {
      baselineexisting=false;
      for(k=0;k<termsmult[i].size();k++)
        {
        if ( baseline.checkvector(termsmult[i],k) == true )
          {
          j = termsmult[i][k].varnames[0].isinlist(modelvarnamesv);
          baselineexisting = true;
          }
        }

      if(baselineexisting == false)
        {
        outerror("ERROR: no baseline specified\n");
        return true;
        }
      }

    vector<unsigned> rescol  = modregmult.getresponsecol();

    // Responsematrix

    unsigned c;

    datamatrix res(D.rows(),rescol.size());
//    unsigned k;
    for(i=0;i<res.cols();i++)
      {
      c = rescol[i];
      for(k=0;k<res.rows();k++)
        res(k,i) = D(k,c);
      }
    if (offs.rows() == 1)
      distr_multistatemodel = DISTRIBUTION_multistatemodel(
                              &generaloptions[generaloptions.size()-1],res,
                              D.getCol(j),beg,statemat,w);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_multistatemodel.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_multistatemodel.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_multistatemodel);
    distrstring.push_back("multistate");
    distrposition.push_back(0);
    nrcategories = res.cols();
    }

//------------------------- END: multistate model ------------------------------

//----------------------- negative binomial response ---------------------------
  else if (family.getvalue() == "nbinomial")
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
#if defined(__BUILDING_LINUX)
    ST::string path3 = defaultpath + "/temp/" + name + add_name + "_scale.raw";
#else
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";
#endif

    MCMC::vertopt vo;

    if (distopt.getvalue() == "nb")
      vo = MCMC::nb;
    else if (distopt.getvalue() == "poga")
      vo = MCMC::poga;
    else
      vo = MCMC::poig;

    MCMC::propscale po;

    if (propopt.getvalue() == "uniform")
      po = MCMC::unif;
    else
      po = MCMC::gam;

    if (offs.rows() == 1)       // without offset
      distr_nbinomial = DISTRIBUTION_nbinomial(aresp.getvalue(),bresp.getvalue(),
      propvar.getvalue(),vo,po,hierarchical.getvalue(),
      &generaloptions[generaloptions.size()-1],D.getCol(0),path2,path3, w);
    else
      distr_nbinomial = DISTRIBUTION_nbinomial(aresp.getvalue(),bresp.getvalue(),
                                     propvar.getvalue(),vo,po,hierarchical.getvalue(),
                                   offs,&generaloptions[generaloptions.size()-1],D.getCol(0),path2,path3,
                                   w);
    distr_nbinomial.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_nbinomial.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_nbinomial.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_nbinomial);
    distrstring.push_back("nbinomial");
    distrposition.push_back(0);
    nrcategories = 1;

    }
//--------------------- END: negative binomial response ------------------------

//----------------------- zip response ---------------------------
  else
    {

    ST::string path2 = outfile.getvalue() + add_name + "_scale.res";
#if defined(__BUILDING_LINUX)
    ST::string path3 = defaultpath + "/temp/" + name + add_name + "_scale.raw";
#else
    ST::string path3 = defaultpath + "\\temp\\" + name + add_name + "_scale.raw";
#endif

    MCMC::zipvertopt vo;

    if (zipdistopt.getvalue() == "zip")
      vo = MCMC::zip;
    else if (zipdistopt.getvalue() == "zinb")
      vo = MCMC::zinb;
    else if (zipdistopt.getvalue() == "zipga")
      vo = MCMC::zipga;
    else
      vo = MCMC::zipig;

    MCMC::zippropscale po;

    if (propopt.getvalue() == "uniform")
      po = MCMC::unifzip;
    else
      po = MCMC::gamzip;

    if (offs.rows() == 1)       // without offset
      distr_zip = DISTRIBUTION_zip(aresp.getvalue(),bresp.getvalue(),
      propvar.getvalue(),vo,po,hierarchical.getvalue(),
      &generaloptions[generaloptions.size()-1],D.getCol(0),path2,path3, w);
    else
      distr_zip = DISTRIBUTION_zip(aresp.getvalue(),bresp.getvalue(),
                                     propvar.getvalue(),vo,po,hierarchical.getvalue(),
                                   offs,&generaloptions[generaloptions.size()-1],D.getCol(0),path2,path3,
                                   w);
    distr_zip.init_names(rname,wn);

    if ((predict.getvalue() == true) || (predictmu.getvalue() == true) )
      distr_zip.set_predict(path,pathdev,&D,modelvarnamesv);

    if (predictmu.getvalue() == true)
      {
      unsigned n;
      if (predictuntil.changed())
        {
        n = predictuntil.getvalue();
        if (n > D.rows())
          n = D.rows();
        }
      else
         n = D.rows();
      distr_zip.set_predictfull(pathfullsample,pathfull,n);
      }

    distr.push_back(&distr_zip);
    distrstring.push_back("zip");
    distrposition.push_back(0);
    nrcategories = 1;

    }
//--------------------- END: zip response --------------------------------------



//----------------- end: reading data, creating designmatrices -----------------

  return false;

  }


bool bayesreg::create_offset(datamatrix & o)
  {
  unsigned i;
  for(i=0;i<terms.size();i++)
    {

    if ( offset.checkvector(terms,i) == true)
      {
      unsigned j = terms[i].varnames[0].isinlist(modelvarnamesv);
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


bool bayesreg::create_const(const unsigned & collinpred)
  {
  unsigned i;
  int j;

  vector<ST::string> varnames;
  vector<ST::string> varnamesh =  fixedeffects.get_constvariables(terms);

  varnames.push_back("const");

  for(i=0;i<varnamesh.size();i++)
    varnames.push_back(varnamesh[i]);

  unsigned nr = varnames.size();
  unsigned blocksize = 10;

  bool shrinkage = false;
  vector<double> variances;
  bool use_effectstart = false;
  vector<double> effectstart;

  if (nr > 0)
    {
    vector< vector<ST::string> > varnamesvec;

    unsigned c=0;
    varnamesvec.push_back( vector<ST::string> () );
    while (c<varnames.size())
      {

      if (c>0 && c%blocksize == 0)
        {
        varnamesvec.push_back( vector<ST::string>() );
        }

      varnamesvec[varnamesvec.size()-1].push_back(varnames[c]);

      c++;
      }

    for(unsigned k=0;k<varnamesvec.size();k++)
      {

      bool constincl = false;
      int constpos = -1;

      nr = varnamesvec[k].size();

      ST::string title;
      ST::string pathconst;
      ST::string pathconstres;

      if (collinpred == 0)
        {
        title = "FixedEffects" + ST::inttostring(k+1) + add_name;
#if defined(__BUILDING_LINUX)
        pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                           + add_name + "_FixedEffects" + ST::inttostring(k+1) + ".raw";
#else
        pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                           + add_name + "_FixedEffects" + ST::inttostring(k+1) + ".raw";
#endif
        pathconstres = outfile.getvalue() + add_name + "_FixedEffects" + ST::inttostring(k+1)
                       + ".res";

        }
      else
        {
        title = "FixedEffects" + ST::inttostring(k+1) + "_" +
                            ST::inttostring(collinpred+1) + add_name;
#if defined(__BUILDING_LINUX)
        pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                           + add_name + "_FixedEffects" + ST::inttostring(k+1) +
                           "_" + ST::inttostring(collinpred+1) + ".raw";
#else
        pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                           + add_name + "_FixedEffects" + ST::inttostring(k+1) +
                           "_" + ST::inttostring(collinpred+1) + ".raw";
#endif
        pathconstres = outfile.getvalue() + add_name + "_FixedEffects" + ST::inttostring(k+1)
                        + "_" + ST::inttostring(collinpred+1) + ".res";

        }

      if (pathconst.isvalidfile() == 1)
        {
        errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
        return true;
        }

      datamatrix X(D.rows(),nr,1);

      if (varnamesvec[k].size() != 0)
        {

        for(i=0;i<varnamesvec[k].size();i++)
          {

          if (varnamesvec[k][i] != "const")
            {

            j = varnamesvec[k][i].isinlist(modelvarnamesv);

            if (j != -1)
              {
              unsigned l;
              double * workX=X.getV()+i;
              double * workD=D.getV()+j;
              for (l=0;l<X.rows();l++,workX+=X.cols(),workD+=D.cols())
                *workX = *workD;
              }

            }
          else
            {
            constincl = true;
            constpos = 0;
            }

          }

        if ( check_gaussian(collinpred))
          {

          if (family.getvalue() == "gaussian_re")
            {

            normalconst_re.push_back(FULLCOND_const_gaussian_re(
            &generaloptions[generaloptions.size()-1],distr[distr.size()-1],X,
                                  title,constpos,pathconst,pathconstres,
                                  shrinkage, variances, use_effectstart, effectstart, collinpred));
            normalconst_re[normalconst_re.size()-1].init_names(varnamesvec[k]);

            normalconst_re[normalconst_re.size()-1].set_fcnumber(fullcond.size());


            if (normalconst_re[normalconst_re.size()-1].get_errors().size() > 0)
              {

              unsigned i;
              for(i=0;i<normalconst_re[normalconst_re.size()-1].get_errors().size();i++)
              errormessages.push_back(
              normalconst_re[normalconst_re.size()-1].get_errors()[i]);
              return true;
              }

            if (constincl == true)
              fcconst_intercept = &normalconst_re[normalconst_re.size()-1];
            fullcond.push_back(&normalconst_re[normalconst_re.size()-1]);

            }
          else
            {
            normalconst.push_back(FULLCOND_const_gaussian(
            &generaloptions[generaloptions.size()-1],distr[distr.size()-1],X,
                                  title,constpos,pathconst,pathconstres,
                                  shrinkage, variances, use_effectstart, effectstart, collinpred));
            normalconst[normalconst.size()-1].init_names(varnamesvec[k]);

            normalconst[normalconst.size()-1].set_fcnumber(fullcond.size());


            if (normalconst[normalconst.size()-1].get_errors().size() > 0)
              {

              unsigned i;
              for(i=0;i<normalconst[normalconst.size()-1].get_errors().size();i++)
              errormessages.push_back(
              normalconst[normalconst.size()-1].get_errors()[i]);
              return true;
              }


            if (constincl == true)
              fcconst_intercept = &normalconst[normalconst.size()-1];
            fullcond.push_back(&normalconst[normalconst.size()-1]);
            }

          } // end: if (family.getvalue() == "gaussian")
        else
          {

          if(family.getvalue() == "nbinomial" && hierarchical.getvalue() == true)
          {
            nbinomialconst.push_back(FULLCOND_const_nbinomial(&generaloptions[generaloptions.size()-1],
                                distr[distr.size()-1],&distr_nbinomial,X,title,constpos,pathconst,
                                pathconstres, shrinkage, variances, use_effectstart, effectstart, collinpred));

            nbinomialconst[nbinomialconst.size()-1].init_names(varnamesvec[k]);

            nbinomialconst[nbinomialconst.size()-1].set_fcnumber(fullcond.size());

          if (constincl == true)
            fcconst_intercept = &nbinomialconst[nbinomialconst.size()-1];

          fullcond.push_back(&nbinomialconst[nbinomialconst.size()-1]);

          } //end: if (family.getvalue() == "nbinomial" && hierarchical.getvalue() == true)

          else
          {
            nongaussianconst.push_back(FULLCOND_const_nongaussian(&generaloptions[generaloptions.size()-1],
                                   distr[distr.size()-1],X,title,constpos,pathconst,pathconstres,
                                   shrinkage, variances, use_effectstart, effectstart, collinpred));
            nongaussianconst[nongaussianconst.size()-1].init_names(varnamesvec[k]);

            nongaussianconst[nongaussianconst.size()-1].set_fcnumber(fullcond.size());


          if (constincl == true)
            fcconst_intercept = &nongaussianconst[nongaussianconst.size()-1];

          if (family.getvalue() == "cox")
            {
            for(unsigned ii=0;ii<fcbaseline.size();ii++)
              fcbaseline[ii].set_fcconst(fcconst_intercept);
//            for(unsigned ii=0;ii<fcbaselineiwls.size();ii++)
//              fcbaselineiwls[ii].set_fcconst(fcconst_intercept);
            }

          fullcond.push_back(&nongaussianconst[nongaussianconst.size()-1]);
          }//end: else "(family.getvalue() == "nbinomial" && hierarchical.getvalue() == true)"
          }// end: else "(family.getvalue() == "gaussian")"

        } // end: if (varnamesvec[k].size() != 0)

      }  // end:   for(k=0;k<varnamesvec.size();k++)

    } // end: if (nr > 0)

  return false;
  }

bool bayesreg::create_nonprw1rw2(const unsigned & collinpred)
  {

  long h;
  double hd;
  unsigned min,max;
  double lambda,a1,b1,alpha;
  bool updatetau;
  bool iwls;
  double ftune;
  unsigned updateW;
  ST::string proposal;
  int f;
  bool varcoeff;
  bool center;

  unsigned i;
  int j1=0,j2=0;

  for(i=0;i<terms.size();i++)
    {
    if ( nonprw1rw2.checkvector(terms,i) == true)
      {

      // -------------- reading options, term information ----------------------

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);

      MCMC::fieldtype type;

      if (terms[i].varnames.size()==1)
        {
        varcoeff=false;
        if (terms[i].options[0] == "rw1" || terms[i].options[0] == "trw1" ||
           terms[i].options[0] == "rw1vrw1" || terms[i].options[0] == "rw1vrw2")
           type = MCMC::RW1;
        else
          type = MCMC::RW2;
        }
      else
        {
        varcoeff=true;
        if (terms[i].options[0] == "varcoeffrw1")
          type = MCMC::RW1;
        else
          type = MCMC::RW2;
        j2 = terms[i].varnames[1].isinlist(modelvarnamesv);
        }


      f = (terms[i].options[1]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[2]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[6]).strtodouble(hd);
      lambda = hd;

      f = (terms[i].options[7]).strtodouble(a1);

      f = (terms[i].options[8]).strtodouble(b1);

      proposal = terms[i].options[9];

      f = (terms[i].options[10]).strtolong(h);
      updateW = unsigned(h);

      if (terms[i].options[11] == "true")
        updatetau=true;
      else
        updatetau=false;

      f = (terms[i].options[12]).strtodouble(ftune);

      f = (terms[i].options[17]).strtodouble(alpha);

      if (terms[i].options[19] == "true")
        center=true;
      else
        center=false;

      if (f==1)
        return true;


      // -------------- reading options, term information ----------------------

      // -------- creating paths for samples and results, titles ---------------

      if (varcoeff==false)
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_rw.raw","_rw.res","_rw");
      else
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                terms[i].varnames[0],
                 "_rw.raw","_rw.res","_rw");

      // -------- end: creating paths for samples and results, titles ----------


      if (proposal != "cp")
        iwls=true;
      else
        iwls=false;


        //-------------------- gaussian response, etc. -------------------------
        if ( (check_gaussian(collinpred)) || (check_iwls(iwls,collinpred)) )
          {

          if (varcoeff==false)
            {
            fcnonpgaussian.push_back(
            FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
            distr[distr.size()-1],D.getCol(j1),fcconst_intercept,
            unsigned(maxint.getvalue()),type,title,pathnonp,pathres,collinpred,
            lambda));

            fcnonpgaussian[fcnonpgaussian.size()-1].init_name(
            terms[i].varnames[0]);
            }
          else
            {
            fcnonpgaussian.push_back(
            FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
            distr[distr.size()-1],D.getCol(j2),D.getCol(j1),fcconst_intercept,
            unsigned(maxint.getvalue()),type,title,pathnonp,pathres,collinpred,
            lambda,center));

            vector<ST::string> na;
            na.push_back(terms[i].varnames[1]);
            na.push_back(terms[i].varnames[0]);
            fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);

            }


          if (constlambda.getvalue() == true)
            {
            if (check_nongaussian(collinpred))
              fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(updateW,true);
            fcnonpgaussian[fcnonpgaussian.size()-1].set_lambdaconst(lambda);
            }
          else
            {

            if ( (check_nongaussian(collinpred)) && (proposal == "iwls")
                && (updatetau==false) )
              fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(updateW);

            if ( (check_nongaussian(collinpred)) && (proposal == "iwlsmode")
               && (updatetau==false) )
               fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(updateW,true);

            if ( (check_nongaussian(collinpred)) && (proposal == "iwls")
                && (updatetau==true) )
              fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS_hyperblock(
                                                                 updateW,a1,b1);

            if ( (check_nongaussian(collinpred)) && (proposal == "iwlsmode")
                && (updatetau==true) )
              fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS_hyperblock(
                                                      updateW,a1,b1,true);
            }

          if (terms[i].options[16] == "true")
            fcnonpgaussian[fcnonpgaussian.size()-1].set_stationary(alpha);

          if (terms[i].options[0] == "rw1vrw1" || terms[i].options[0] == "rw2vrw1"
          || terms[i].options[0] == "rw1vrw2" || terms[i].options[0] == "rw2vrw2"
          || terms[i].options[0] == "trw1" || terms[i].options[0] == "trw2")
            fcnonpgaussian[fcnonpgaussian.size()-1].set_adaptiv();

          fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);

          if (
               (terms[i].options[0] != "rw1vrw1") &&
               (terms[i].options[0] != "rw2vrw1") &&
               (terms[i].options[0] != "rw1vrw2") &&
               (terms[i].options[0] != "rw2vrw2") &&
               (constlambda.getvalue() == false)
             )
            {

            if (varcoeff==false)
              make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_rw_var.raw","_rw_var.res","_rw_variance");
            else
              make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                  terms[i].varnames[0],"_rw_var.raw","_rw_var.res","_rw_variance");

            fcvarnonp.push_back(
            FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
            &fcnonpgaussian[fcnonpgaussian.size()-1],distr[distr.size()-1],
            a1,b1,title,pathnonp,pathres,false,collinpred));

            if ( (check_nongaussian(collinpred)) && (proposal == "iwls")
              && (updatetau==true) )
              fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();

            if ( (check_nongaussian(collinpred)) && (proposal == "iwlsmode")
              && (updatetau==true) )
              fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();

            bool alphafix = false;
            if (terms[i].options[18] == "true")
              alphafix = true;
            if (terms[i].options[16] == "true")
              fcvarnonp[fcvarnonp.size()-1].set_stationary(alpha,alphafix);

            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

            }

          if (terms[i].options[0] == "trw1" || terms[i].options[0] == "trw2")
            {

            if (varcoeff==false)
              make_paths(collinpred,pathnonp,pathres,title,
                         terms[i].varnames[0],"",
                         "_rw_tvar.raw","_rw_tvar.res","_rw_tvariance");
            else
              make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                  terms[i].varnames[0],"_rw_tvar.raw","_rw_tvar.res",
                  "_rw_tvariance");

            unsigned v = nu.getvalue();

            fctvariance.push_back(FULLCOND_tvariance(
            &generaloptions[generaloptions.size()-1],
            &fcnonpgaussian[fcnonpgaussian.size()-1],v,title,pathnonp,pathres)
                                );
            fctvariance[fctvariance.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fctvariance[fctvariance.size()-1]);

            }

          if (terms[i].options[0] == "rw1vrw1" || terms[i].options[0] == "rw2vrw1"
          || terms[i].options[0] == "rw1vrw2" || terms[i].options[0] == "rw2vrw2")
            {

            if (varcoeff==false)
              {
              make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_rw_vvar.raw","_rw_vvar.res","_rw_vvariance");
              }
            else
              make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                  terms[i].varnames[0],"_rw_vvar.raw","_rw_vvar.res",
                  "_rw_vvariance");


            MCMC::fieldtype tv;
            if (terms[i].options[0] == "rw1vrw1" || terms[i].options[0] == "rw2vrw1")
              tv = MCMC::RW1;
            else
              tv = MCMC::RW2;

            f= (terms[i].options[3]).strtolong(h);
            unsigned minvar = unsigned(h);

            f= (terms[i].options[4]).strtolong(h);
            unsigned maxvar = unsigned(h);

            f= (terms[i].options[5]).strtodouble(hd);
            double startv = hd;

           if (f==1)
             return true;

            fcadaptiv.push_back(
            FULLCOND_adaptiv(&generaloptions[generaloptions.size()-1],
            &fcnonpgaussian[fcnonpgaussian.size()-1],tv,title,a1,b1,false,startv,
            minvar,maxvar,pathnonp,pathres));

            fcadaptiv[fcadaptiv.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcadaptiv[fcadaptiv.size()-1]);

            }

          //------------------- end: gaussian response, etc. -------------------
          }
        else
          {
          //------------------- non-gaussian response, etc. --------------------

          if (terms[i].options[0] == "rw1vrw1" || terms[i].options[0] == "rw2vrw1"
          || terms[i].options[0] == "rw1vrw2" || terms[i].options[0] == "rw2vrw2"
          || terms[i].options[0] == "trw1" || terms[i].options[0] == "trw2")
            {
            outerror("ERROR: '" + terms[i].options[0] + "' not available\n");
            return true;
            }

          if (varcoeff==false)
            Pmatrices.push_back(PenaltyMatrix(D.getCol(j1),title,
            unsigned(maxint.getvalue()),min,max,type));
          else
            Pmatrices.push_back(PenaltyMatrix(D.getCol(j2),terms[i].varnames[1],
            unsigned(maxint.getvalue()),min,max,type));

          fcnonp.push_back(
          FULLCOND_nonp(&generaloptions[generaloptions.size()-1],
          distr[distr.size()-1],&Pmatrices[Pmatrices.size()-1],
          fcconst_intercept,lambda,pathnonp,pathres,title," ",collinpred));

          fcnonp[fcnonp.size()-1].init_name(terms[i].varnames[0]);

          fcnonp[fcnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcnonp[fcnonp.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_rw_var.raw","_rw_var.res","_rw_variance");

          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
          &fcnonp[fcnonp.size()-1],distr[distr.size()-1],a1,b1,title,pathnonp,
          pathres,false,collinpred));

          if (constlambda.getvalue() == true)
            fcvarnonp[fcvarnonp.size()-1].set_constlambda();


          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

          //------------------- end: non-gaussian response, etc. ---------------
          }

      } // end: if ( nonprw1rw2.checkvector(terms,i) == true )

    }

  return false;

  }


bool bayesreg::create_pspline(const unsigned & collinpred)
  {

  ST::string monotone;
  ST::string proposal;

  long h;
  unsigned min,max,degree,nrknots;
  double lambda,a1,b1,alpha,merrorvar;
  bool ub,diagtransform,derivative,bsplinebasis,discretize;
  int gridsize=0,contourprob=0,digits=0,nobs=0;
  int f=0;
  ST::string test ="test";
  double lowerknot=0;
  double upperknot=0;
  double lowergrid=0;
  double uppergrid=0;
  bool testmerror=false;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonppspline.checkvector(terms,i) == true )
      {

      // --------------- reading options, term information ---------------------
      MCMC::fieldtype type;
      if ( (terms[i].options[0] == "psplinerw1") ||
           (terms[i].options[0] == "tpsplinerw1") ||
           (terms[i].options[0] == "psplinerw1vrw1") ||
           (terms[i].options[0] == "psplinerw1vrw2") )
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

// BEGIN: merror
      datamatrix meandata;
      datamatrix medata;
      if(terms[i].varnames[0].length()>7)
        {
        test = terms[i].varnames[0].substr(terms[i].varnames[0].length()-7,7);
        }
      if(test=="_merror")
        // covariates measured with measurement error
        {
        testmerror=true;
        if(merror.getvalue()==0)
          {
          outerror("ERROR: Option merror has not been specified\n");
          return true;
          }
        unsigned me=merror.getvalue();
        test = terms[i].varnames[0].substr(0,terms[i].varnames[0].length()-7);
        medata = datamatrix(D.rows(),me,0);
        meandata = datamatrix(D.rows(),1,0);
        unsigned k;
        for(k=0; k<me; k++)
          {
          j = (test+ST::inttostring(k+1)).isinlist(modelvarnamesv);
          medata.putCol(k, D.getCol(j));
          }
        terms[i].varnames[0] = test;
        for(k=0; k<D.rows(); k++)
          {
          for(j=0; j<me; j++)
            {
            meandata(k,0) += medata(k,j);
            }
          meandata(k,0) /= me;
          }
        if (terms[i].options[35] == "false")
          discretize = false;
        else
          discretize = true;

        f = terms[i].options[36].strtolong(h);
        digits = unsigned(h);
        f = terms[i].options[37].strtolong(h);
        nobs = unsigned(h);
        if(discretize && nobs==0)
          {
          outerror("ERROR: Option nobs has to be specified if discretize is true\n");
          return true;
          }
        if(nobs==0)
          nobs=D.rows();
        meandata.round(digits,0,1,0,nobs);
        }
      else
        {
        j = terms[i].varnames[0].isinlist(modelvarnamesv);
        meandata = D.getCol(j);
        }
// END: merror

      f = (terms[i].options[1]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[2]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[3]).strtolong(h);
      degree = unsigned(h);

      f = (terms[i].options[4]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[5]).strtodouble(lambda);

      f = (terms[i].options[6]).strtodouble(a1);

      f = (terms[i].options[7]).strtodouble(b1);

      if (terms[i].options[8] == "false")
        ub = false;
      else
        ub = true;

      f = (terms[i].options[9]).strtolong(h);
      gridsize = unsigned(h);

      if (f==1)
        return true;

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant" && terms[i].options[29] == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      proposal = terms[i].options[13];
      monotone = terms[i].options[14];

      if (terms[i].options[18] == "false")
        diagtransform = false;
      else
        diagtransform = true;

      if (terms[i].options[19] == "false")
        derivative = false;
      else
        derivative = true;

      if (terms[i].options[20] == "false")
        bsplinebasis = false;
      else
        bsplinebasis = true;

      f = (terms[i].options[21]).strtolong(h);
      contourprob = unsigned(h);

      f = (terms[i].options[27]).strtodouble(alpha);

// BEGIN: merror
     f = (terms[i].options[30]).strtodouble(lowerknot);
     f = (terms[i].options[31]).strtodouble(upperknot);

     f = (terms[i].options[33]).strtodouble(lowergrid);
     f = (terms[i].options[34]).strtodouble(uppergrid);

      if(testmerror)
        // covariates measured with measurement error
        // Overwrite some of the options
        {
//        derivative=true;               // for IWLS-proposal
        if(gridsize<5)
          gridsize = 500;              // evaluate the function on a grid

       f = (terms[i].options[32]).strtodouble(merrorvar);

       if(lowerknot==upperknot)
         {
         double xmin = meandata.min(0);
         double xmax = meandata.max(0);
         lowerknot = xmin - 3*sqrt(merrorvar);
         upperknot = xmax + 3*sqrt(merrorvar);
         }
       }

// END: merror

      // -------------end: reading options, term information -------------------

      //--------- creating path for samples and and results, creating title ----

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline.raw","_pspline.res","_pspline");

      //----- end: creating path for samples and and results, creating title ---

      // ---------------------- gaussian response, etc. ------------------------
      if (check_gaussian(collinpred))
        {

        fcpsplinegaussian.push_back(
        FULLCOND_pspline_gaussian(&generaloptions[generaloptions.size()-1],
        distr[distr.size()-1],fcconst_intercept,meandata,nrknots,degree,po,
        type,monotone,title,pathnonp,pathres,derivative,lambda,gridsize,
        diagtransform,lowerknot,upperknot,lowergrid,uppergrid,collinpred));

        //
        datamatrix beta_0;
        if(terms[i].options[23]!="")
          {
          dataobject * datap;                           // pointer to datasetobject
          int objpos = findstatobject(*statobj,terms[i].options[23],"dataset");
          if (objpos >= 0)
            {
            statobject * s = statobj->at(objpos);
            datap = dynamic_cast<dataobject*>(s);
            if (datap->obs()==0 || datap->getVarnames().size()==0)
              {
              outerror("ERROR: dataset object " + terms[i].options[23] + " does not contain any data\n");
              return true;
              }
            else if (datap->getVarnames().size()>1)
              {
              outerror("ERROR: dataset object " + terms[i].options[23] + " contains more than one variable\n");
              return true;
              }
            }
          else
            {
            outerror("ERROR: dataset object " + terms[i].options[23] + " is not existing\n");
            return true;
            }
          list<ST::string> names = datap->getVarnames();
          ST::string expr = "";
          datap->makematrix(names,beta_0,expr);
          }
        else
          {
          beta_0 = datamatrix(1,1,0);
          }
        //

        fcpsplinegaussian[fcpsplinegaussian.size()-1].set_contour(contourprob,
            pseudocontourprob.getvalue(),approx.getvalue(),lengthstart.getvalue(),beta_0);

        if (constlambda.getvalue() == true)
          fcpsplinegaussian[fcpsplinegaussian.size()-1].set_lambdaconst(lambda);

        if (bsplinebasis == true)
          fcpsplinegaussian[fcpsplinegaussian.size()-1].set_outbsplines();

        if (terms[i].options[26] == "true")
          fcpsplinegaussian[fcpsplinegaussian.size()-1].set_stationary(alpha);

        fcpsplinegaussian[fcpsplinegaussian.size()-1].init_name(terms[i].varnames[0]);
        fcpsplinegaussian[fcpsplinegaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcpsplinegaussian[fcpsplinegaussian.size()-1]);

        if(terms[i].options[0] != "psplinerw1vrw1" && terms[i].options[0] != "psplinerw1vrw2" &&
           terms[i].options[0] != "psplinerw2vrw1" && terms[i].options[0] != "psplinerw2vrw2")
          {

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                     "_pspline_var.raw","_pspline_var.res","_pspline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(
          &generaloptions[generaloptions.size()-1],
          &fcpsplinegaussian[fcpsplinegaussian.size()-1],distr[distr.size()-1],
          a1,b1,title,pathnonp,pathres,ub,collinpred));

          if (constlambda.getvalue() == false)
            {

            if(terms[i].options[22]=="true")
               fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            if(terms[i].options[24]=="true")
              {
              f = terms[i].options[25].strtolong(h);
              fcvarnonp[fcvarnonp.size()-1].set_discrete(unsigned(h));
              }
            bool alphafix = false;
            if (terms[i].options[28] == "true")
              alphafix = true;
            if (terms[i].options[26] == "true")
              fcvarnonp[fcvarnonp.size()-1].set_stationary(alpha,alphafix);

            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          }

        if (terms[i].options[0] == "tpsplinerw1" || terms[i].options[0] == "tpsplinerw2")
          {
          fcpsplinegaussian[fcpsplinegaussian.size()-1].set_adaptiv();

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_tvar.raw","_pspline_tvar.res","_pspline_tvariance");

          unsigned v = nu.getvalue();

          fctvariance.push_back(FULLCOND_tvariance(&generaloptions[generaloptions.size()-1],
                                &fcpsplinegaussian[fcpsplinegaussian.size()-1],
                                v,title,pathnonp,pathres));

          fctvariance[fctvariance.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fctvariance[fctvariance.size()-1]);

          } // end: if (terms[i].options[0] == "tpsplinerw1" ...

        if ( terms[i].options[0] == "psplinerw1vrw1" || terms[i].options[0] == "psplinerw1vrw2" ||
             terms[i].options[0] == "psplinerw2vrw1" || terms[i].options[0] == "psplinerw2vrw2" )
          {
          fcpsplinegaussian[fcpsplinegaussian.size()-1].set_adaptiv();

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_vvar.raw","_pspline_vvar.res","_pspline_vvariance");


          f = (terms[i].options[10]).strtolong(h);
          unsigned minvar = unsigned(h);
          f = (terms[i].options[11]).strtolong(h);
          unsigned maxvar = unsigned(h);
          double startv;
          f = (terms[i].options[12]).strtodouble(startv);
          MCMC::fieldtype typevar;
          if (terms[i].options[0] == "psplinerw1vrw1" ||
          terms[i].options[0] == "psplinerw2vrw1")
            typevar = MCMC::RW1;
          else
            typevar = MCMC::RW2;

          fcadaptiv.push_back(FULLCOND_adaptiv(&generaloptions[generaloptions.size()-1],
                              &fcpsplinegaussian[fcpsplinegaussian.size()-1],
                              typevar,title,a1,b1,ub,startv,minvar,maxvar,pathnonp,pathres));

          fcadaptiv[fcadaptiv.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcadaptiv[fcadaptiv.size()-1]);

          } // end: if (terms[i].options[0] == "psplinerw1vrw1" ...

//BEGIN: merror
        if(testmerror)
        // Fullcond-Objekt zur Generierung der wahren Kovariablenwerte
          {
          fcpsplinegaussian[fcpsplinegaussian.size()-1].set_changingweight();

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                       "_merror.raw","_merror.res","_merror");

          fcmerror.push_back(fullcond_merror(&generaloptions[generaloptions.size()-1],
                         &fcpsplinegaussian[fcpsplinegaussian.size()-1],
                         distr[distr.size()-1],
                                   medata,
                                   title,
                                   pathnonp,
                                   pathres,
                                   lowerknot, upperknot,
                                   merrorvar,
                                   discretize, digits, nobs)
                         );
          fcmerror[fcmerror.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcmerror[fcmerror.size()-1]);
          }
//END: merror

        // ---------------------- end: gaussian response, etc. -----------------
        }
      else
        {
        // -------------------- non-gaussian response, etc. --------------------

        if (proposal == "cp")
          {

          if (terms[i].options[26] == "true")
            {
            outerror("ERROR: option 'alpha' not available\n");
            return true;
            }

          if (terms[i].options[0] == "tpsplinerw1" || terms[i].options[0] == "tpsplinerw2" ||
              terms[i].options[0] == "psplinerw1vrw1" || terms[i].options[0] == "psplinerw1vrw2" ||
              terms[i].options[0] == "psplinerw2vrw1" || terms[i].options[0] == "psplinerw2vrw2" )
            {
            outerror("ERROR: '" + terms[i].options[0] + "' not available\n");
            return true;
            }

          fcpspline.push_back( FULLCOND_pspline(&generaloptions[generaloptions.size()-1],
          distr[distr.size()-1],fcconst_intercept,meandata,nrknots,degree,po,
          lambda,min,max,type,title,pathnonp,pathres,derivative,lowerknot,
          upperknot,lowergrid,uppergrid,gridsize,
          collinpred));

          if (constlambda.getvalue() == true)
            fcpspline[fcpspline.size()-1].set_lambdaconst(lambda);

          fcpspline[fcpspline.size()-1].init_name(terms[i].varnames[0]);
          fcpspline[fcpspline.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpspline[fcpspline.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_var.raw","_pspline_var.res","_pspline_variance");


          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
          &fcpspline[fcpspline.size()-1],distr[distr.size()-1],a1,b1,title,
          pathnonp,pathres,ub,collinpred));

          if (constlambda.getvalue() == false)
            {

            if(terms[i].options[22]=="true")
                fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

//BEGIN: merror
          if(testmerror)
          // Fullcond-Objekt zur Generierung der wahren Kovariablenwerte
            {
            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                       "_merror.raw","_merror.res","_merror");

            fcmerror.push_back(fullcond_merror(&generaloptions[generaloptions.size()-1],
                           &fcpspline[fcpspline.size()-1],
                           distr[distr.size()-1],
                                   medata,
                                   title,
                                   pathnonp,
                                   pathres,
                                   lowerknot, upperknot,
                                   merrorvar,
                                   discretize, digits, nobs)
                           );
            fcmerror[fcmerror.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcmerror[fcmerror.size()-1]);
            }
//END: merror
          }
        else // iwls
          {

          bool iwlsmode;
          if(proposal == "iwlsmode")
            iwlsmode = true;
          else
            iwlsmode = false;
          f = (terms[i].options[15]).strtolong(h);
          unsigned updateW;
          updateW = unsigned(h);
          bool updatetau;
          if(terms[i].options[16] == "false" || constlambda.getvalue() == true)
            updatetau = false;
          else
            updatetau = true;

          double fstart;
            f = (terms[i].options[17]).strtodouble(fstart);

          fciwlspspline.push_back(
          IWLS_pspline(&generaloptions[generaloptions.size()-1],
          distr[distr.size()-1],fcconst_intercept,meandata,iwlsmode,nrknots,degree,po,
          lambda,type,monotone,updateW,updatetau,fstart,a1,b1,title,pathnonp,
          pathres,derivative,gridsize,diagtransform,lowerknot,upperknot,
          lowergrid,uppergrid,
          collinpred));

          if (constlambda.getvalue() == true)
            fciwlspspline[fciwlspspline.size()-1].set_lambdaconst(lambda);

          if (terms[i].options[26] == "true")
            fciwlspspline[fciwlspspline.size()-1].set_stationary(alpha);
          if (bsplinebasis == true)
            fciwlspspline[fciwlspspline.size()-1].set_outbsplines();

          fciwlspspline[fciwlspspline.size()-1].init_name(terms[i].varnames[0]);
          fciwlspspline[fciwlspspline.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fciwlspspline[fciwlspspline.size()-1]);

          if(terms[i].options[0] != "psplinerw1vrw1" && terms[i].options[0] != "psplinerw1vrw2" &&
             terms[i].options[0] != "psplinerw2vrw1" && terms[i].options[0] != "psplinerw2vrw2")
            {

            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_var.raw","_pspline_var.res","_pspline_variance");

            fcvarnonp.push_back(
            FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
            &fciwlspspline[fciwlspspline.size()-1],distr[distr.size()-1],a1,b1,
            title,pathnonp,pathres,ub,collinpred));

            bool alphafix = false;
            if (terms[i].options[28] == "true")
              alphafix = true;
            if (terms[i].options[26] == "true")
              fcvarnonp[fcvarnonp.size()-1].set_stationary(alpha,alphafix);

            if (constlambda.getvalue() == false)
              {

              if(terms[i].options[22]=="true")
                fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
              if(updatetau)
                fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();
              fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
              fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
              }

            }

          if (terms[i].options[0] == "tpsplinerw1" || terms[i].options[0] == "tpsplinerw2")
            {
            fciwlspspline[fciwlspspline.size()-1].set_adaptiv();

            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_tvar.raw","_pspline_tvar.res","_pspline_tvariance");

            unsigned v = nu.getvalue();

            fctvariance.push_back(FULLCOND_tvariance(&generaloptions[generaloptions.size()-1],
            &fciwlspspline[fciwlspspline.size()-1],v,title,pathnonp,pathres)
                                );

            fctvariance[fctvariance.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fctvariance[fctvariance.size()-1]);

            } // end: if (terms[i].options[0] == "tpsplinerw1" ...

          if (terms[i].options[0] == "psplinerw1vrw1" || terms[i].options[0] == "psplinerw1vrw2" ||
              terms[i].options[0] == "psplinerw2vrw1" || terms[i].options[0] == "psplinerw2vrw2")
            {
            fciwlspspline[fciwlspspline.size()-1].set_adaptiv();

            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_vvar.raw","_pspline_vvar.res","_pspline_vvariance");

            f = (terms[i].options[10]).strtolong(h);
            unsigned minvar = unsigned(h);
            f = (terms[i].options[11]).strtolong(h);
            unsigned maxvar = unsigned(h);
            double startv;
            f = (terms[i].options[12]).strtodouble(startv);
            MCMC::fieldtype typevar;
            if (terms[i].options[0] == "psplinerw1vrw1" || terms[i].options[0] == "psplinerw2vrw1")
              typevar = MCMC::RW1;
            else
              typevar = MCMC::RW2;

            fcadaptiv.push_back(FULLCOND_adaptiv(&generaloptions[generaloptions.size()-1],
                                &fciwlspspline[fciwlspspline.size()-1],
                                typevar,title,a1,b1,ub,startv,minvar,maxvar,pathnonp,pathres));

            fcadaptiv[fcadaptiv.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcadaptiv[fcadaptiv.size()-1]);

            } // end: if (terms[i].options[0] == "psplinerw1vrw1" ...

          }

//BEGIN: merror
          if(testmerror)
          // Fullcond-Objekt zur Generierung der wahren Kovariablenwerte
            {
            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                       "_merror.raw","_merror.res","_merror");

            fcmerror.push_back(fullcond_merror(&generaloptions[generaloptions.size()-1],
                           &fciwlspspline[fciwlspspline.size()-1],
                           distr[distr.size()-1],
                                   medata,
                                   title,
                                   pathnonp,
                                   pathres,
                                   lowerknot, upperknot,
                                   merrorvar,
                                   discretize, digits, nobs)
                           );
            fcmerror[fcmerror.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcmerror[fcmerror.size()-1]);
            }
//END: merror
        // ----------------- end:  non-gaussian response, etc. -----------------
        }

      }

    }

  return false;
  }

#if defined(BORLAND_OUTPUT_WINDOW)
//------------------------------------------------------------------------------
#pragma package(smart_init)
#endif











