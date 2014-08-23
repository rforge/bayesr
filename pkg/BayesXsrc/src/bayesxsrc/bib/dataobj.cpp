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

#include<describe_dataset.h>
#endif

#include "dataobj.h"


void dataobject::create(void)
  {

  srand((unsigned)time(0));
//  randomize();

  m = modelStandard();
  e = expression();

  emptyoptions = optionlist();

  uread = usePathRead();

  uwrite = usePathWrite();

// SYNTAX OF COMMANDS:
// name [model] [weight varname] [by varname] [if expression]
//      [, options] [using usingtext]

  // method infile

  missing = stroption("missing");
  maxobs = intoption("maxobs",10000,10,10000000);
  nonote = simpleoption("nonote",false);

  infileoptions.push_back(&missing);
  infileoptions.push_back(&maxobs);
  infileoptions.push_back(&nonote);

  methods.push_back(command("infile",&m,&infileoptions,&uread,
			 optional,notallowed,notallowed,notallowed,optional,required));
  functions[0] = infilerun;

  // method drop

  methods.push_back(command("drop",&m,&emptyoptions,&uread,optional,notallowed,
									 notallowed,optional,notallowed,notallowed));

  functions[1] = droprun;

  // method rename

  methods.push_back(command("rename",&m,&emptyoptions,&uread,required,
						  notallowed,notallowed,notallowed,notallowed,notallowed));
  functions[2] = renamerun;

  // method generate

  methods.push_back(command("generate",&e,&emptyoptions,&uread,required,
									 notallowed,notallowed,notallowed,notallowed,
									 notallowed));

  functions[3] = generaterun;

  // method replace

  methods.push_back(command("replace",&e,&emptyoptions,&uread,required,
									 notallowed,notallowed,optional,notallowed,
									 notallowed));

  functions[4] = replacerun;

  // method set obs

  methods.push_back(command("set",&e,&emptyoptions,&uread,required,notallowed,
									 notallowed,notallowed,notallowed,notallowed));

  functions[5] = setrun;

  // method outfile

  header = simpleoption("header",false);
  replace = simpleoption("replace",false);
  outfileoptions.push_back(&header);
  outfileoptions.push_back(&replace);

  methods.push_back(command("outfile",&m,&outfileoptions,&uwrite,optional,
									 notallowed,notallowed,optional,optional,required));

  functions[6] = outfilerun;

// SYNTAX OF COMMANDS:
// name [model] [weight varname] [by varname] [if expression]
//      [, options] [using usingtext]

  // method sort

  descending = simpleoption("descending",false);
  sortoptions.push_back(&descending);

  methods.push_back(command("sort",&m,&sortoptions,&uread,required,notallowed,
									 notallowed,notallowed,optional,notallowed));

  functions[7] = sortrun;

  // method descriptive

  methods.push_back(command("descriptive",&m,&descriptiveoptions,&uread,required,
                    notallowed,notallowed,optional,notallowed,notallowed));

  functions[8] = descriptiverun;

  //method tabulate

  methods.push_back(command("tabulate",&m,&tabulateoptions,&uread,required,
                    notallowed,notallowed,optional,notallowed,notallowed));

  functions[9] = tabulaterun;

  //method pctile

  methods.push_back(command("pctile",&m,&pctileoptions,&uread,required,
                    notallowed,notallowed,optional,notallowed,notallowed));

  functions[10] = pctilerun;

  //method marketing

  lag = intoption("lag",1,1,100);
  alpha = doubleoption("alpha",0.05,0,1);

  vector<ST::string> prdefs;
  prdefs.push_back("notransform");
  prdefs.push_back("regular");
  prdefs.push_back("priceindex");

  pricedef = stroption("pricedef",prdefs,"regular");

  marketingoptions.push_back(&lag);
  marketingoptions.push_back(&alpha);
  marketingoptions.push_back(&pricedef);


  methods.push_back(command("marketing",&m,&marketingoptions,&uread,required,
  notallowed,notallowed,notallowed,optional,notallowed));

  functions[11] = marketingrun;

  }


void dataobject::changedescription(void)
  {
  unsigned long s = d.varnr()*d.obs()*(sizeof(double));
  describetext.erase(describetext.begin(),describetext.end());
  describetext.push_back("Number of variables:    " + ST::inttostring(d.varnr()));
  describetext.push_back("Number of observations: " + ST::inttostring(d.obs()));
  if (s >= 1024)
    {
    double sd = double(s)/1024.0;
    describetext.push_back("Size of dataset:        "
                           + ST::doubletostring(sd,6) + " Kb");
    }
  else
    describetext.push_back("Size of dataset:        " + ST::inttostring(s));
  }


#if defined(JAVA_OUTPUT_WINDOW)
dataobject::dataobject(administrator_basic * adb, administrator_pointer * adp,
                       const ST::string & n,ofstream * lo,istream * in)
			 : statobject(adb,n,"dataset",lo,in)
	 {
     adminp_p = adp;
	 d = dataset(n,adb);
	 create();
	 }
#else
dataobject::dataobject(const ST::string & n,ofstream * lo,istream * in)
			 : statobject(n,"dataset",lo,in)
	 {
	 d = dataset(n);
	 create();
	 }
#endif


dataobject::dataobject(const dataobject & o) : statobject(statobject(o))
  {
  create();
  d = o.d;
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = o.adminp_p;
  #endif
  }


const dataobject & dataobject::operator=(const dataobject & o)
  {
  if (this == &o)
	 return *this;
  statobject::operator=(statobject(o));
  create();
  d = o.d;
  #if defined(JAVA_OUTPUT_WINDOW)
  adminp_p = o.adminp_p;
  #endif
  return *this;
  }


bool dataobject::allexisting(vector<ST::string> & varnames,
                             vector<ST::string> & notex)
  {
  int i=0;
  bool notfound = false;
  while (i<varnames.size())
	 {
	 if (d.findvar(varnames[i]) == 1)
		{
		notfound=true;
		notex.push_back(varnames[i]);
		}

	 i++;
	 }

  if (notfound == false)
	 return true;
  else
	 return false;
  }


int dataobject::parse(const ST::string & c)
  {
  optionlist globaloptions = optionlist();
  int pos = statobject::parsecom(c,methods,globaloptions);

  if (pos >= 0)
	 (*functions[pos])(*this);

  return(pos);
  }


void infilerun(dataobject & o)
  {

  ST::string path = o.uread.getPath();
  ST::string missingvalue = o.missing.getvalue();
  list<ST::string> names = o.m.getModelVarnames();
  ifstream fin;
  ST::open(fin,path);

  #if defined(JAVA_OUTPUT_WINDOW)
  o.d.read(o.adminb_p,fin,missingvalue,o.maxobs.getvalue(),names);
  #else
  o.d.read(fin,missingvalue,o.maxobs.getvalue(),names);
  #endif
  fin.close();
  o.errormessages = o.d.geterrormessages();
  if ((o.errormessages.empty()) && (o.nonote.getvalue()==false))
     {
	 o.out(
	 "NOTE: " + ST::inttostring(o.d.varnr()) + " variables with " +
	 ST::inttostring(o.d.obs()) + " observations read from file\n");
     o.out(path + "\n");
     o.out("\n");
     o.changedescription();
     }

  }


void droprun(dataobject & o)
  {
  list<ST::string> names = o.m.getModelVarnames();
  ST::string boolexp = o.methods[1].getexpression();

  if ( (names.size() > 0) && (boolexp.length() == 0) )
	 {
	 o.d.dropvariables(names);
	 o.errormessages = o.d.geterrormessages();
	if (o.errormessages.empty())
	  o.out("NOTE: " + ST::inttostring(names.size()) + " variables dropped\n");
	 }
  else if ( (names.size() == 0) && (boolexp.length() > 0) )
	 {
     unsigned nrelim;
	 nrelim = o.d.dropobservations(boolexp);
	 o.errormessages = o.d.geterrormessages();
	 if (o.errormessages.empty())
		o.out("NOTE: " + ST::inttostring(nrelim) + " observations dropped\n");
	 }
  else if ( (names.size() == 0) && (boolexp.length() == 0) )
	 {
	 o.outerror("ERROR: varlist or boolean expression expected\n");
	 }
  else
    o.outerror("ERROR: dropping variables and observations in one step not allowed\n");

  o.changedescription();
  }


void renamerun(dataobject & o)
  {
  list<ST::string> names = o.m.getModelVarnames();
  if (names.size() != 2)
	 o.errormessages.push_back("ERROR: invalid syntax\n");
  else
	 {
	 list<ST::string>::iterator i1,i2;
	 i1 = names.begin();
	 i2 = names.begin();
	 i2++;
	 o.d.rename(*i1,*i2);
	 o.errormessages = o.d.geterrormessages();
	 }
  }


void generaterun(dataobject & o)
  {
  o.errormessages = o.e.geterrormessages();

  if (o.errormessages.empty())
	 {
	 o.e.addvariable(o.d);
	 o.errormessages = o.e.geterrormessages();
     }

  if ( (o.errormessages.empty()) && (o.d.obs() == 0) )
    o.out("WARNING: number of observations in dataset is zero\n",true,true);

  o.changedescription();
  }


void replacerun(dataobject & o)
  {
  ST::string name = o.e.getvarname();
  ST::string expression = o.e.getexpression();
  ST::string ifexpression = o.methods[4].getexpression();
  unsigned changed = o.d.replace(name,expression,ifexpression);
  o.errormessages = o.d.geterrormessages();
  if (o.errormessages.empty())
    o.out("NOTE: " + ST::inttostring(changed) + " observations changed\n");
  o.changedescription();
  }


void setrun(dataobject & o)
  {
  int old = o.d.obs();
  if (o.e.getvarname() != "obs")
	 o.errormessages.push_back("ERROR: invalid syntax\n");
  else
	 {
	 long value;
	 if ((o.e.getexpression()).strtolong(value) == 1)
		o.errormessages.push_back("ERROR: integer value expected\n");
	 else
		{
		o.d.setobs(int(value));
		o.errormessages = o.d.geterrormessages();
		}
	 }
  if (o.errormessages.empty())
	 o.out(
	 "NOTE: number of observations raised from " + ST::inttostring(old) +
	 " to " + ST::inttostring(o.d.obs()) + "\n" );
  o.changedescription();
  }


void outfilerun(dataobject & o)
  {

  unsigned nrwritten;
  ST::string path  = o.uwrite.getPath();
  list<ST::string> names = o.m.getModelVarnames();
  ST::string expression = o.methods[6].getexpression();

  if ( (o.uwrite.isexisting() == true) && (o.replace.getvalue() == false) )
	 o.errormessages.push_back(
	 "ERROR: file " + path + " is already existing\n");
  else
	 {
	 ofstream fout;
	 ST::open(fout,path);
	 if (expression.length() > 0)
		{
		realvar v = o.d.eval_exp(expression);
        #if defined(JAVA_OUTPUT_WINDOW)
		nrwritten = o.d.write(o.adminb_p,fout,names,o.header.getvalue(),v);
        #else
        nrwritten = o.d.write(fout,names,o.header.getvalue(),v);
        #endif
		}
	 else
       {
        #if defined(JAVA_OUTPUT_WINDOW)
	   nrwritten = o.d.write(o.adminb_p,fout,names,o.header.getvalue());
       #else
	   nrwritten = o.d.write(fout,names,o.header.getvalue());
       #endif
       }

     o.errormessages = o.d.geterrormessages();
     if (o.errormessages.empty())
       {
       o.out("NOTE: " + ST::inttostring(names.size()) + " variable(s) with " +
             ST::inttostring(nrwritten) +  " observations written to file\n");
       o.out("      " + path + "\n");
       }
     else
       {
       fout.close();
       remove(path.strtochar());
       }


	 }
  }


void sortrun(dataobject & o)
  {
  list<ST::string> names = o.m.getModelVarnames();
  if (o.d.obs() > 0)
    {
    o.d.sort(names,0,o.d.obs()-1);
    o.errormessages = o.d.geterrormessages();
    if ( (o.errormessages.empty()) && (o.descending.getvalue() == true))
      o.d.reverseorder();
    }
  else
    o.out("NOTE: dataset contains no data\n");
  }


void descriptiverun(dataobject & o)
  {

  // datamatrix X;

  unsigned i;

  vector<ST::string> vnames = o.m.getModelVarnamesAsVector();

  // testing, wether all variables specified are already existing
  vector<ST::string> notex;
  bool failure=false;
  if ((o.allexisting(vnames,notex)) == false)
    {
    for (i=0;i<notex.size();i++)
      o.outerror("ERROR: variable " + notex[i] + " is not existing\n");

    failure = true;

    } // end: if ((datap->allexisting(modelvarnamesv,notex)) == false)


  ST::string ifexpression = o.methods[8].getexpression();

  if (failure == false)
    {

    // o.makematrix(vnames,X,ifexpression);

    o.out("\n");
    ST::string varst = "Variable";
    ST::string obsst = "Obs";
    ST::string meanst = "Mean";
    ST::string medst = "Median";
    ST::string stdst = "Std";
    ST::string minst = "Min";
    ST::string maxst ="Max";
    o.out(varst.helpfill(11) + obsst.helpfill(10) + meanst.helpfill(15) + medst.helpfill(15)
           + stdst.helpfill(15) + minst.helpfill(15) + maxst.helpfill(15) + "\n",true);

    for(i=0;i<vnames.size();i++)
      {
//      o.m.makeModelMatrix_j(o.d,X,i);
      datamatrix Xi;
      o.makematrix(vnames[i],Xi,ifexpression);      // für jede Variable eigene Matrix erzeugen?

      unsigned int obs = Xi.rows();
      double mean = Xi.mean(0);
      double var = 0;
      if(obs > 1)
        var = Xi.var(0) * obs / (obs-1);

      ST::string str_std;
      if (var>0)
        {
        double std = sqrt(var);
        str_std = ST::doubletostring(std,8);
        }
      else
        {
        int std = 0;
        str_std = ST::inttostring(std);
        }
      double min = Xi.min(0);
      double max = Xi.max(0);
      double med = Xi.quantile(50,0);

      int l_vn = vnames[i].length();
      if (l_vn > 9)
         {
         vnames[i] = vnames[i].substr(0, 9) + "~";
         }
      ST::string str_obs = ST::inttostring(obs);
      ST::string str_mean = ST::doubletostring(mean,8);
      ST::string str_med = ST::doubletostring(med,8);
      ST::string str_min = ST::doubletostring(min,8);
      ST::string str_max = ST::doubletostring(max,8);

      o.out(vnames[i].helpfill(11) + str_obs.helpfill(10) + str_mean.helpfill(15) +
            str_med.helpfill(15) + str_std.helpfill(15) + str_min.helpfill(15) +
            str_max.helpfill(15) + "\n");

      }
    o.out("\n");
    }

  }


void tabulaterun(dataobject & o)
  {

  // datamatrix X;

  unsigned i;

  vector<ST::string> vnames = o.m.getModelVarnamesAsVector();

  // testing, wether all variables specified are already existing
  vector<ST::string> notex;
  bool failure=false;
  if ((o.allexisting(vnames,notex)) == false)
    {
    for (i=0;i<notex.size();i++)
      o.outerror("ERROR: variable " + notex[i] + " is not existing\n");

    failure = true;

    } // end: if ((datap->allexisting(modelvarnamesv,notex)) == false)

  ST::string ifexpression = o.methods[9].getexpression();


  if (failure==false)
    {

    // o.makematrix(vnames,X,ifexpression);

    // statmatrix<int> index(X.rows(),1);

    for(i=0;i<vnames.size();i++)
        {
//        o.m.makeModelMatrix_j(o.d,X,i);
        datamatrix Xi;
        o.makematrix(vnames[i],Xi,ifexpression);
        statmatrix<int> index(Xi.rows(),1);

        index.indexinit();

        Xi.indexsort(index,0,Xi.rows()-1,0,0);

        o.out("\n");
        o.out("Variable: " + vnames[i] + "\n",true, false, 14);

        o.out("\n");

        unsigned k = 0;
        double cum = 0;
        unsigned zaehler = 0;

        while(k<index.rows() && zaehler<101)
            {
            double anz=0;

            int* p = index.getV() + k;
            int* q = index.getV() + k;
            for(unsigned j=k;j<index.rows();j++,p++)
                 {
                 if (Xi.get(*p,0) == Xi.get(*q,0))
                    anz = anz+1;
                 }

                 k = k + anz;
                 zaehler = zaehler + 1;

            }

        if(zaehler>=101)
           o.outerror("ERROR: too many values");
        else
           {
           ST::string value = "Value";
           ST::string obs = "Obs";
           ST::string freqst = "Freq";
           ST::string cumst = "Cum";
           o.out(value.helpfill(15) + obs.helpfill(10) +
                 freqst.helpfill(15) + cumst.helpfill(15) + "\n",true);
           k = 0;

           while(k<index.rows())
                 {
                 double anz=0;

                 int* p = index.getV() + k;
                 int* q = index.getV() + k;
                 for(unsigned j=k;j<index.rows();j++,p++)
                    {
                    if (Xi.get(*p,0) == Xi.get(*q,0))
                       anz = anz+1;
                    }

                 double wert;
                 wert = Xi.get(*q,0);
                 double freq = anz / Xi.rows();
                 cum = cum + freq;

                 ST::string str_wert = ST::doubletostring(wert,8);
                 ST::string str_anz = ST::inttostring(anz);
                 ST::string str_freq = ST::doubletostring(freq,4);
                 ST::string str_cum = ST::doubletostring(cum,4);

                 o.out(str_wert.helpfill(15) + str_anz.helpfill(10) +
                       str_freq.helpfill(15) + str_cum.helpfill(15) + "\n");

                 k = k + anz;
                 }
           }

        }
    o.out("\n");
    }
  }

void pctilerun(dataobject & o)
  {

  datamatrix X;

  unsigned i;

  vector<ST::string> vnames = o.m.getModelVarnamesAsVector();

  ST::string ifexpression = o.methods[10].getexpression();

  // testing, wether all variables specified are already existing
  vector<ST::string> notex;
  bool failure=false;
  if ((o.allexisting(vnames,notex)) == false)
    {
    for (i=0;i<notex.size();i++)
      o.outerror("ERROR: variable " + notex[i] + " is not existing\n");

    failure = true;

    } // end: if ((datap->allexisting(modelvarnamesv,notex)) == false)

  if (failure==false)
    {

    o.makematrix(vnames,X,ifexpression);

    o.out("\n");

    //Ausgabe:
    //ST::string varst = "Variable";
    //ST::string

    for(i=0;i<vnames.size();i++)
       {
 //       o.m.makeModelMatrix_j(o.d,X,i);
        o.out("Variable: " + vnames[i] + "\n",true, false, 14);
        o.out("\n");

        for(unsigned j=1; j<6; j++)
           {
           double quant = X.quantile(j,i);
           ST::string quantst = ST::doubletostring(quant,8);
           o.out(" " + ST::inttostring(j) + "%" + quantst.helpfill(15) + "\n");
           }
        unsigned j = 25;
        while(j<=75)
           {
           double quant = X.quantile(j,i);
           ST::string quantst = ST::doubletostring(quant,8);
           o.out(ST::inttostring(j) + "%" + quantst.helpfill(15) + "\n");
           j = j + 25;
           }
        for(unsigned j=95; j<100; j++)
           {
           double quant = X.quantile(j,i);
           ST::string quantst = ST::doubletostring(quant,8);
           o.out(ST::inttostring(j) + "%" + quantst.helpfill(15) + "\n");
           }
           o.out("\n");
        }
   }

 }


void marketingrun(dataobject & o)    //Reihenfolge zum Einlesen: outlet, wochenin, markenin, preis!!!
  {

  ST::string defs = o.pricedef.getvalue(); //Herauslesen, welcher Preisindex berechnet wird

  int lak = o.lag.getvalue();  //Herauslesen, aus welcher Vorwoche die Preise sein sollen
  double alph = o.alpha.getvalue();

      //Herauslesen der Variablennamen
  vector<ST::string> vnames = o.m.getModelVarnamesAsVector();
          //Überprüfen der Anzahl der übergebenen Variablen
  if(vnames.size()>4)
    o.errormessages.push_back("ERROR: too many variables");
  else if(vnames.size()<4)
    o.errormessages.push_back("ERROR: not enough variables");
  else     //Überprüfen, ob die Variablen existieren
    {
    bool richtig = true;
    for(unsigned i=0; i<vnames.size();i++)
      {
      int test = o.d.findvar(vnames[i]);
      if(test == 1)
        {
        o.errormessages.push_back("ERROR: variable '" + vnames[i] + "' is not existing \n");
        richtig = false;
        }
      }
    if(richtig==true)
      {
      o.d.marketing(vnames,defs,lak,alph);  //Aufruf der Funktion "marketing" aus Datei "data"
      }
    }
  }

void dataobject::describe(const optionlist & globaloptions)
  {
  if(d.getVarnames().size()>0 && d.obs() > 0)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
    datasetform->datap = &d;
    datasetform->dataname = getname();
    datasetform->ShowModal();
#elif defined(JAVA_OUTPUT_WINDOW)
    adminp_p->set_datap(&d);
    jmethodID javashowdata = adminb_p->Java->GetMethodID(
    adminb_p->BayesX_cls, "JavaShowData", "()V");
    adminb_p->Java->CallVoidMethod(adminb_p->BayesX_obj, javashowdata);
#elif(__BUILDING_GNU)
    out("ERROR: method describe is not available in this version\n");
#endif
  }
  else
    out("NOTE: dataset does not contain any data\n");
  }


#if defined(BORLAND_OUTPUT_WINDOW)
//---------------------------------------------------------------------------
#pragma package(smart_init)
#endif




