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



#include "randomeffect_stepwise.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------- class FULLCOND_random_stepwise -----------------------------
//------------------------------------------------------------------------------

FULLCOND_random_stepwise::FULLCOND_random_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & d, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,
                              const double & la, const unsigned & c)
                            : FULLCOND_random(o,dp,fcc,d,t,fp,pr,la,c)
  {
  utype = "gaussian";

  identifiable = false;
  intercept = 0.0;
  }

// randomslope
FULLCOND_random_stepwise::FULLCOND_random_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & intvar,const datamatrix & effmod,
                  const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const ST::string & prf,
                  const double & la,
                  const bool & inclfixed,const unsigned & c)
                  : FULLCOND_random(o,dp,fcc,intvar,effmod,t,fp,pr,prf,la,false,c)
  {
  utype = "gaussian";

  if(inclfixed == false)
    {
    identifiable = false;
    }

  includefixed = false;
  intercept = 0.0;
  get_data_forfixedeffects();
  }


FULLCOND_random_stepwise::FULLCOND_random_stepwise(const FULLCOND_random_stepwise & fc)
                            : FULLCOND_random(FULLCOND_random(fc))
  {
  intercept = fc.intercept;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  fbasisp = fc.fbasisp;
  df_unstruct = fc.df_unstruct;
  fc_df = fc.fc_df;
  utype = fc.utype;
  //gleichwertig = fc.gleichwertig;
  }


const FULLCOND_random_stepwise & FULLCOND_random_stepwise::
         operator=(const FULLCOND_random_stepwise & fc)
  {
  if (this==&fc)
    return *this;

  FULLCOND_random::operator=(FULLCOND_random(fc));

  intercept = fc.intercept;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  fbasisp = fc.fbasisp;
  df_unstruct = fc.df_unstruct;
  fc_df = fc.fc_df;
  utype = fc.utype;
  //gleichwertig = fc.gleichwertig;

  return *this;
  }


void FULLCOND_random_stepwise::set_nofixed(bool fix)
  {
  FULLCOND::set_nofixed(fix);

  if(!nofixed && identifiable)
    {
    includefixed = true;
    setbeta(beta.rows()+1,1,0);
    }

  nrpar = beta.rows();
  }


bool FULLCOND_random_stepwise::posteriormode(void)
  {
  unsigned n = nrpar;
  if (includefixed)
    n = nrpar-1;

  update_linpred(false);

  if(calculate_xwx == true)
    {
    calculate_xwx = false;
    compute_XWX(likep->get_weightiwls(),column);
    }

  likep->compute_weightiwls_workingresiduals(column);

  unsigned i,j;
  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  int * workindex2 = index2.getV();
  itbeg = posbeg.begin();
  itend = posend.begin();
  double * workmuy = muy.getV();
  likep->set_workingresp();

  if (!randomslope)
    {
    for(i=0;i<nrpar;i++,workmuy++,++itbeg,++itend)
      {
      *workmuy = 0;
      for(j=*itbeg;j<=*itend;j++,workindex2++)
        {
        *workmuy+= likep->get_workingres(*workindex2);
        }
      }
    }
  else
    {
    double * datap = data.getV();
    if(includefixed)
      {
      double ms = beta(nrpar-1,0);
      likep->set_linpredp_current(column);
      for (i=0;i<n;i++,workmuy++,++itbeg,++itend)
        {
        *workmuy = 0;
        for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
          *workmuy += likep->get_workingres(*workindex2)* (*datap);

        *workmuy+= lambda* ms;
        }
      }
    else
      {
      for(i=0;i<n;i++,workmuy++,++itbeg,++itend)
        {
        *workmuy = 0;
        for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
          {
          *workmuy+= likep->get_workingres(*workindex2)* (*datap);
          }
        }
      }
    }

  itbeg = posbeg.begin();
  itend = posend.begin();
  workmuy = muy.getV();
  double * workbeta = beta.getV();
  double * workXX = XX.getV();

  for(i=0;i<n;i++,workmuy++,++itbeg,++itend,workbeta++,workXX++)
    {
    *workbeta = (*workmuy)/(*workXX+lambda);
    }

  if (randomslope && (center || includefixed))
    {
    double * workbeta = beta.getV();
    double sum=0;
    for (i=0;i<n;i++,workbeta++)
      {
      sum += *workbeta;
      }
    intercept = sum/double(n);

    if(includefixed)
      beta(nrpar-1,0) = intercept;
    else
      {
      update_linpred(true);
      update_fix_effect(intercept);
      }

    workbeta = beta.getV();
    for (i=0;i<n;i++,workbeta++)
      *workbeta -= intercept;
    intercept = 0.0;

    if(includefixed)
      update_linpred(true);
    }
  else
    update_linpred(true);

  transform = likep->get_trmult(column);

  return FULLCOND::posteriormode();

  }


// BEGIN: For Varying Coefficients ---------------------------------------------

void FULLCOND_random_stepwise::update_fix_effect(double & intercept)
  {
  bool raus = false;
  unsigned j = 1;
  ST::string name_richtig = datanames[1];  // (VC alt) 0
  while(j<fcconst->get_datanames().size() && raus==false)
     {
     if(fcconst->get_datanames()[j] == datanames[1])  // (VC alt) 0
        {
        raus = true;
        }
     if(fcconst->get_datanames()[j] == (datanames[1]+"_1"))   // (VC alt) 0
        {
        raus = true;
        name_richtig = datanames[1] + "_1";     // (VC alt) 0
        }
     j = j + 1;
     }
  if(raus == true)
    {
    fcconst->update_fix_effect(j-1,intercept,data_forfixed);
    }
  else
    {
    vector<ST::string> names;
    names.push_back(name_richtig);
    fcconst->include_effect(names,data_forfixed);
    interactions_pointer[0]->set_inthemodel(-1);
    fcconst->update_fix_effect(j,intercept,data_forfixed);
    }
  }


void FULLCOND_random_stepwise::set_pointer_to_interaction(FULLCOND * inter)
  {
  interactions_pointer.push_back(inter);
  }

void FULLCOND_random_stepwise::get_interactionspointer(vector<FULLCOND*> & inter)
  {
  inter = interactions_pointer;
  }


void FULLCOND_random_stepwise::hierarchical(ST::string & possible)
  {
  if(!randomslope)
    {
    possible = "alles";
    }
  else
    {
    possible = "valles";
    }
  }


void FULLCOND_random_stepwise::const_varcoeff(void)
  {
  if(randomslope)
    fcconst->posteriormode_const_varcoeff(data_forfixed);
  }

// END: For Varying Coefficients -----------------------------------------------


void FULLCOND_random_stepwise::create_weight(datamatrix & w)
  {
  if(!spatialtotal)
    {
    unsigned n = nrpar;
    if (includefixed)
      n = nrpar-1;

    vector<unsigned>::iterator itbeg = posbeg.begin();
    vector<unsigned>::iterator itend = posend.begin();
    int * workindex = index.getV();
    unsigned i;
    unsigned j;
    for(i=0;i<n;i++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        w(*workindex,0) = 1;
        for(j=*itbeg;j<=*itend;j++)
          workindex++;
        }
      }
    }
  }


void FULLCOND_random_stepwise::compute_lambdavec(vector<double> & lvec, int & number)
  {
//ST::string sp = "automatic";
//ST::string sp = "nonautomatic";
  if(spfromdf=="automatic")
    {
    df_equidist = true;
    double maxi = floor(nrpar/4.0*3.0);

    if(maxi <= 30)
      {
      df_for_lambdamin = maxi;
      if(!varcoeff || !identifiable)
        {
        df_for_lambdamax = 1;
        number = maxi;
        }
      else
        {
        df_for_lambdamax = 2;
        number = maxi - 1;
        }
      }
    else if(maxi > 30 && maxi<=60)
      {
      number = floor(maxi/2);
      df_for_lambdamax = 2;
      df_for_lambdamin = number*2;
      }
    else if(maxi > 60 && maxi<=100)
      {
      df_for_lambdamax = 3;
      number = floor(maxi/3);
      df_for_lambdamin = number*3;
      }
    else if(maxi > 100 && maxi<=180)
      {
      df_for_lambdamax = 5;
      number = floor(maxi/5);
      df_for_lambdamin = number*5;
      }
    else if(maxi > 180)
      {
      df_for_lambdamax = 10;
      number = floor(maxi/10);
      df_for_lambdamin = number*10;
      }
    }

  if (df_equidist==true && spfromdf!="direct" && number>1)
     FULLCOND::compute_lambdavec_equi(lvec,number);
  else
     FULLCOND::compute_lambdavec(lvec,number);
  if (!nofixed && randomslope && identifiable)
    hierarchie_fix(lvec,1);
  if(forced_into==false)
     lvec.push_back(0);

  // Startwert für lambda aus df:
  if(spfromdf!="direct")
    {
    double lambdavorg = 1000;
    if(!randomslope)
      {
      if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    else
      {
      if(!nofixed && dfstart==1 && identifiable)
        lambdastart = -1;
      else if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    if(lambdastart==-9 || lambdastart==1000000000)    // falls dfstart nicht erreicht werden kann
      lambdastart = 0;
    }

  }


void FULLCOND_random_stepwise::hierarchie_fix(vector<double> & untervector, int dfo)
  {

  unsigned number = untervector.size()-1;

  update_stepwise(untervector[0]);
  double df_max = compute_df();

  update_stepwise(untervector[number]);
  double df_min = compute_df();

  if(df_max > dfo && df_min < dfo)
     {
     bool geordnet = false;
     unsigned stelle_oben = number;
     unsigned stelle_unten = 0;
     while(geordnet==false)
        {
        unsigned stelle = stelle_oben + stelle_unten;
        update_stepwise(untervector[stelle/2]);
        double df_mitteunten = compute_df();
        update_stepwise(untervector[stelle/2 + 1]);
        double df_mitteoben = compute_df();

        if(df_mitteunten > dfo && df_mitteoben > dfo)
          stelle_unten = stelle/2;
        else if(df_mitteunten < dfo && df_mitteoben < dfo)
          stelle_oben = stelle/2 + 1;
        else
          {
          geordnet = true;
          vector<double> hilf;
          unsigned i;
          stelle_unten = stelle/2;
          stelle_oben = stelle/2 + 1;
          for(i=0;i<=stelle_unten;i++)
             hilf.push_back(untervector[i]);
          hilf.push_back(-1);
          for(i=stelle_oben;i<untervector.size();i++)
            hilf.push_back(untervector[i]);
          untervector = hilf;
          }
        }
     }
  else if(df_min >= dfo)
     untervector.push_back(-1);
  else
     {
     vector<double> hilf;
     hilf.push_back(-1);
     unsigned i;
     for(i=0;i<untervector.size();i++)
        hilf.push_back(untervector[i]);
     untervector = hilf;
     }
  }


const datamatrix & FULLCOND_random_stepwise::get_data_forfixedeffects(void)
  {
  // useful for randomslopes only
  if ( (data_forfixed.rows() < data.rows()) && (randomslope==true) )
    {
    data_forfixed=datamatrix(data.rows(),1);
    unsigned i;
    int * workindex = index.getV();
    double * workdata = data.getV();
    for (i=0;i<data.rows();i++,workindex++,workdata++)
      {
      data_forfixed(*workindex,0) = *workdata;
      }
    }

  return data_forfixed;
  }


double FULLCOND_random_stepwise::compute_df(void)
  {
  double df = 0;
  if(inthemodel == true)
    {
    unsigned i;
    bool struct_included = false;
    if(spatialtotal)
      {
      bool fix;
      fbasisp->get_inthemodel(struct_included,fix);
      if(struct_included == true)
        {
        df = df_unstruct;
        }
      }

    if(!spatialtotal || !struct_included)
      {
      if (lambdaold1==lambda && likep->get_iwlsweights_notchanged() == true && !spatialtotal)
        {
        df = df_lambdaold1;
        }
      else
        {
        if(calculate_xwx == true)
          {
          calculate_xwx = false;
          compute_XWX(likep->get_weightiwls(),column);
          }
        double * workXX=XX.getV();
        unsigned n = nrpar;
        if(includefixed)
          n = nrpar - 1;

        if(!identifiable)
          {
          double c = 0;
          double w = 0;
          for(i=0;i<n;i++,workXX++)
            {
            c += (*workXX * *workXX)/(*workXX+lambda);
            w += *workXX;
            }
         c = 1/(w - c);

         workXX = XX.getV();
         for(i=0;i<n;i++,workXX++)
           {
           df += (*workXX * (lambda + *workXX * (-c * (*workXX + 2*lambda) + 1)))/((*workXX+lambda) * (*workXX+lambda));
           }
         df += w*c - 1;
         }
        else
          {
          for(i=0;i<n;i++,workXX++)
            df += (*workXX)/(*workXX+lambda);
          if(includefixed)
            df += 1;
          }

        df_lambdaold1 = df;
        lambdaold1 = lambda;
        }
      }
    }

  return df;
  }


void FULLCOND_random_stepwise::set_dfunstruct(const double & df_unstr)
  {
  df_unstruct = df_unstr;
  }


void FULLCOND_random_stepwise::update_stepwise(double la)
  {
  lambda = la;
  }

double FULLCOND_random_stepwise::get_lambda(void)
  {
  return lambda;
  }


ST::string FULLCOND_random_stepwise::get_effect(void)
  {

  ST::string h;

  if(randomslope)
    h = datanames[1] + "*" + datanames[0];   // (VC alt) 0 - 1
  else
    h = datanames[0];

  h = h + "(random,df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;

  }


void FULLCOND_random_stepwise::init_names(const vector<ST::string> & na)
    {
    char charh = '_';
    ST::string stringh = "\\_";

    FULLCOND::init_names(na);
    if (na.size()==1)
      {
      ST::string helpname = na[0].insert_string_char(charh,stringh);
      term_symbolic = "f_{" +  helpname + "}("+helpname+")";
      }
    else
      {
      ST::string helpname1 = na[0].insert_string_char(charh,stringh);     // (VC alt) 1
      ST::string helpname2 = na[1].insert_string_char(charh,stringh);     // (VC alt) 0
      term_symbolic = "f_{" +  helpname1 + "}("+helpname1+") \\cdot "
                        + helpname2;
      }

    if (column > 0)
      priorassumptions.push_back("$" + term_symbolic + "$" +
       " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");

    priorassumptions.push_back("i.i.d. Gaussian random effects");
    }


void FULLCOND_random_stepwise::reset_effect(const unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  intercept = 0.0;
  }


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_random_stepwise::update_bootstrap(const bool & uncond)
  {
  update_bootstrap_df();

  datamatrix betaold = beta;

  if(fixornot==true)
    {
    bool raus = false;
    unsigned j = 1;
    ST::string name_richtig;
    if(!randomslope)
      name_richtig = datanames[0];
    else
      name_richtig = datanames[1]; // (VC alt) 0
    while(j<fcconst->get_datanames().size() && raus==false)
      {
      if(fcconst->get_datanames()[j] == name_richtig)
        raus = true;
      j = j + 1;
      }
    unsigned index_fix = j-1;
    double fix = fcconst->getbeta(index_fix,0);
    unsigned i;
    double * workbeta = beta.getV();
    vector<unsigned>::iterator itbeg = posbeg.begin();
    vector<unsigned>::iterator itend = posend.begin();
    //int * workindex = index2.getV();
    //int k;
    for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        *workbeta = fix;
        //for(k=*itbeg;k<=*itend;k++)
        //  workindex++;
        }
      }
    workbeta = beta.getV();

    FULLCOND::update_bootstrap();
    }
  else if(inthemodel==false && fixornot==false)
    {
    beta = datamatrix(nrpar,1,0);
    FULLCOND::update_bootstrap();
    }
  else
    {
    FULLCOND::update_bootstrap();
    }
  beta = betaold;
  }


void FULLCOND_random_stepwise::update_bootstrap_df(void)
  {
  if(optionsp->get_samplesize()<=1)
    {
    ST::string path = samplepath.substr(0,samplepath.length()-4)+"_df.raw";
    fc_df = FULLCOND(optionsp,datamatrix(1,1),"title?",1,1,path);
    fc_df.setflags(MCMC::norelchange | MCMC::nooutput);
    }
  if(fixornot==true)
    {
    fc_df.setbetavalue(0,0,-1);
    fc_df.update_bootstrap_df();
    }
  else if(inthemodel==false && fixornot==false)
    {
    fc_df.setbetavalue(0,0,0);
    fc_df.update_bootstrap_df();
    }
  else
    {
    fc_df.setbetavalue(0,0,lambda);
    fc_df.update_bootstrap_df();
    }
  }


void FULLCOND_random_stepwise::save_betamean(void)
  {
  datamatrix betaold = beta;

  if(fixornot==true)
    {
    bool raus = false;
    unsigned j = 1;
    ST::string name_richtig;
    if(!randomslope)
      name_richtig = datanames[0];
    else
      name_richtig = datanames[1];  // (VC alt) 0
    while(j<fcconst->get_datanames().size() && raus==false)
      {
      if(fcconst->get_datanames()[j] == name_richtig)
        raus = true;
      j = j + 1;
      }
    unsigned index_fix = j-1;
    double fix = fcconst->getbeta(index_fix,0);
    unsigned i;
    double * workbeta = beta.getV();
    vector<unsigned>::iterator itbeg = posbeg.begin();
    vector<unsigned>::iterator itend = posend.begin();
    //int * workindex = index2.getV();
    //int k;
    for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        *workbeta = fix;
        //for(k=*itbeg;k<=*itend;k++)
        //  workindex++;
        }
      }
    workbeta = beta.getV();

    FULLCOND::save_betamean();
    }
  else if(inthemodel==false && fixornot==false)
    {
    beta = datamatrix(nrpar,1,0);
    FULLCOND::save_betamean();
    }
  else
    {
    FULLCOND::save_betamean();
    }
  beta = betaold;
  }

void FULLCOND_random_stepwise::update_bootstrap_betamean(void)
  {
  FULLCOND::update_bootstrap_betamean();
  FULLCOND::setflags(MCMC::norelchange);
  }

void FULLCOND_random_stepwise::outresults_df(unsigned & size)
  {
  fc_df.update_bootstrap_betamean();
  //fc_df.outresults();
  double * workmean = fc_df.get_betameanp();

  ST::string pathdf = pathcurrent.substr(0,pathcurrent.length()-4)+"_df.res";
  ofstream outres(pathdf.strtochar());

  outres << "value   ";
  outres << "frequency  ";
  outres << "selected  " << endl;

// Häufigkeitstabelle:

  //samplestream.close();
  datamatrix sample(size,1);
  fc_df.readsample_df(sample,0);
  unsigned i;

  vector<unsigned> number;
  vector<unsigned> number1;
  vector<unsigned> number2;
  vector<unsigned> cumnumber1;
  vector<unsigned> cumnumber;

  statmatrix<int> index(sample.rows(),1);
  index.indexinit();
  sample.indexsort(index,0,sample.rows()-1,0,0);

  i = 0;
  unsigned j,anz;
  while(i<index.rows())
     {
     anz=0;
     int* p = index.getV() + i;
     int* q = index.getV() + i;
     j=i;
     while(j<index.rows() && (sample.get(*p,0) == sample.get(*q,0)))
        {
        anz = anz+1;
        j++;
        p++;
        }
     if(sample.get(*q,0) <= 0)
       number1.push_back(anz);
     else if(sample.get(*q,0) > 0)
       number2.push_back(anz);
     if(cumnumber1.size()>0)
       cumnumber1.push_back(cumnumber1[cumnumber1.size()-1]+anz);
     else
       cumnumber1.push_back(anz);
     i = i + anz;
     }

  int k;
  for(k=number1.size()-1;k>=0;k--)
    {
    cumnumber.push_back(cumnumber1[k]);
    number.push_back(number1[k]);
    }
  for(k=number2.size()-1;k>=0;k--)
    {
    cumnumber.push_back(cumnumber1[k+number1.size()]);
    number.push_back(number2[k]);
    }

  for(i=0;i<number.size();i++)
    {
    double help = sample.get(index(cumnumber[i]-1,0),0);
    double dfs = -1*help;
    if(help>0)
      {
      update_stepwise(help);
      set_inthemodel(help);
      dfs = compute_df();
      }
    outres << ST::doubletostring(dfs,6) << "   " << ST::inttostring(number[i]) << "   ";
    if(*workmean == help)
      outres << "+"; // ST::doubletostring(*workmean,6);
    else
      outres << "-";
    outres << endl;
    }
  }


void FULLCOND_random_stepwise::update(void)
  {
  if(lambda==0)
    {
    beta = datamatrix(beta.rows(),beta.cols(),0);
    FULLCOND_random::update();
    }
  else if(lambda==1000000000 && randomslope && identifiable && !includefixed)
    {
    update_linpred(false);
    unsigned i;
    double slope = 0;
    double xwx = 0;
    double xwy = 0;
    double * wresp = likep->get_tildey().getV();
    double * wlin = likep->get_linearpred().getV();
    double * wdat = data_forfixed.getV();
    double * workw = likep->get_weightiwls().getV();
    for(i=0;i<data_forfixed.rows();i++,wresp++,wlin++,wdat++,workw++)
      {
      xwx += *wdat * *wdat * *workw;
      xwy += *wdat * *workw * (*wresp - *wlin);
      }
    xwx = 1/xwx;
    slope = sqrt(likep->get_scale(column)) * sqrt(xwx) * rand_normal();
    slope += xwx * xwy;
    double * workbeta = beta.getV();
    vector<unsigned>::iterator itbeg = posbeg.begin();
    vector<unsigned>::iterator itend = posend.begin();
    for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        *workbeta = slope;
        }
      }
    update_linpred(true);
    FULLCOND_random::update();
    }
  else
    {
    if(utype == "gaussian")
      update_gauss();
    else
      update_nongauss();
    }
  }


void FULLCOND_random_stepwise::update_gauss(void)
  {
  double var;
  double m;
  unsigned i,j;
  unsigned n = nrpar;
  if(randomslope && includefixed)
    n = nrpar-1;

  if (optionsp->get_nriter()==1 || changingweight)
    compute_XWX(likep->get_weight(),0);

  sigma2 = likep->get_scale(column)/lambda;  // nur kontantes Lambda
  double sqrtscale = sqrt(likep->get_scale(column));
  update_linpred(false);

  datamatrix mu = datamatrix(index.rows(),1,0);
  // nicht verändern wegen SUR-Modellen
  likep->compute_respminuslinpred(mu,column);

  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  double * workbeta = beta.getV();

  int * workindex2 = index2.getV();
  double * workmuy = muy.getV();
  double * mup = mu.getV();
  likep->set_weightp();

  if (!randomslope)
    {
    for(i=0;i<nrpar;i++,workmuy++,++itbeg,++itend)
      {
      *workmuy = 0;
      for(j=*itbeg;j<=*itend;j++,workindex2++)
        {
        mup += *workindex2;
        *workmuy+= likep->get_weight(*workindex2)* *mup;
        }
      }
    }
  else
    {
    double * datap = data.getV();
    for(i=0;i<n;i++,workmuy++,++itbeg,++itend)
      {
      *workmuy = 0;
      for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
        {
        mup += *workindex2;
        *workmuy+= likep->get_weight(*workindex2)* (*mup) * (*datap);
        }
      if (includefixed)
        *workmuy += beta(n,0)*lambda;
      }
    }

  workbeta = beta.getV();
  workmuy = muy.getV();
  double * workXX = XX.getV();
  for(i=0;i<n;i++,workbeta++,workmuy++,workXX++)
    {
    var = 1.0/(*workXX  + lambda);
    m = var * *workmuy;
    *workbeta = m + sqrtscale*sqrt(var)*rand_normal();
    }

  if(randomslope && (includefixed || center))
    {
    workbeta = beta.getV();
    double s=0;
    for (i=0;i<n;i++,workbeta++)
      s += *workbeta;
    s /= double(n);

    if(includefixed)
      {
      double v = sigma2/double(nrpar-1);
      beta(nrpar-1,0) = s+sqrt(v)*rand_normal();
      s = beta(nrpar-1,0);
      }
    else
      {
      update_linpred(true);
      fcconst->update_fix_varcoeff(s,datanames[1]);
      }

    workbeta = beta.getV();
    for (i=0;i<n;i++,workbeta++)
      *workbeta -= s;

    if(includefixed)
      update_linpred(true);
    }
  else
    update_linpred(true);

  acceptance++;
  transform = likep->get_trmult(column);

  FULLCOND_random::update();
  update_spatialtotal();
  }


void FULLCOND_random_stepwise::update_nongauss(void)  // entspricht "update_random_intercept_iwls_singleblock"
  {                                                   //            "update_random_slope_includefixed_iwls_singleblock"
  unsigned n = nrpar;                                 //            "update_random_slope_iwls_singleblock"
  if(randomslope && includefixed)
    n = nrpar-1;
  double ms = 0;
  if(randomslope && includefixed)
    ms = beta(nrpar-1,0);

  unsigned i;
  double sumw;
  double sumy;
  double var;
  double proposal;
  double u;

  double logold;
  double lognew;
  double qnew;
  double qold;

  double diff;
  double mode;

  double * workbeta = beta.getV();
  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();

  if (lambdaconst == false)
    lambda=1.0/sigma2;
  else
    sigma2 = likep->get_scale(column)/lambda;  // ist nicht "sigma2 = phi/lambda" richtig???

  for (i=0;i<n;i++,workbeta++,++itbeg,++itend)
    {
    nrtrials++;

    // Compute loglikelihood based on current beta
    // compute sum weightiwls  (corresponds to X'WX) and sum tildey
    double help = *workbeta;
    if(randomslope && includefixed)
      help += ms;
    if(!randomslope)
      logold = likep->compute_loglikelihood_sumweight_sumy(help,sumw,sumy,
                                             *itbeg,*itend,index,index2,column);
    else
      logold = likep->compute_loglikelihood_sumweight_sumy(help,sumw,sumy,
                                        *itbeg,*itend,data,index,index2,column);

    sumw = sumw/likep->get_scale(column);
    sumy = sumy/likep->get_scale(column);
    // log-prior densities
    logold -= 0.5*(*workbeta)*(*workbeta)/sigma2;

    if(randomslope && includefixed)
      sumy += ms/sigma2;

    // compute proposal
    var = 1.0/(sumw + 1.0/sigma2);
    mode = var*sumy;
    proposal = mode + sqrt(var)*rand_normal();

    // log q(beta_c,proposal)
    diff = proposal-mode;
    diff = diff*diff;
    qold = -(1.0/(2*var))*diff-0.5*log(var);

    // add proposed beta and substract current beta from linear predictor
    diff = proposal-*workbeta;
    if(randomslope && includefixed)
      diff -= ms;
    if(!randomslope)
      likep->add_linearpred2(diff,*itbeg,*itend,index,index2,column);
    else
      likep->add_linearpred2(diff,*itbeg,*itend,data,index,index2,column);

    // Compute log-likelihood based on proposed beta
    // compute sum weightiwls  (corresponds to X'WX) and sum
    // tildey based on proposed beta
    if(!randomslope)
      lognew = likep->compute_loglikelihood_sumweight_sumy(proposal,sumw,sumy,
                                             *itbeg,*itend,index,index2,column);
    else
      lognew = likep->compute_loglikelihood_sumweight_sumy(proposal,sumw,sumy,
                                        *itbeg,*itend,data,index,index2,column);

    sumw = sumw/likep->get_scale(column);
    sumy = sumy/likep->get_scale(column);
    // log-prior densities
    if(randomslope && includefixed)
      {
      lognew -= 0.5*(proposal-ms)*(proposal-ms)/sigma2;
      sumy += ms/sigma2;
      }
    else
      lognew -= 0.5*(proposal)*(proposal)/sigma2;

    // log q(proposal,beta_c)
    var = 1.0/(sumw + 1.0/sigma2);
    mode = var*sumy;
    diff = *workbeta-mode;
    if(randomslope && includefixed)
      diff += ms;
    diff = diff*diff;
    qnew = -(1.0/(2*var))*diff-0.5*log(var);

    // accept/reject proposal
    u = log(uniform());
    if (u <= (lognew + qnew - logold - qold) )
      {
      acceptance++;
      *workbeta = proposal;
      }
    else
      {
      if(randomslope && includefixed)
        *workbeta += ms;

      if(!randomslope)
        likep->add_linearpred2(*workbeta-proposal,*itbeg,*itend,index,index2,column);
      else
        likep->add_linearpred2(*workbeta-proposal,*itbeg,*itend,data,index,index2,column);
      }
    } // end:   for (i=0;i<nrpar;i++,workbeta++)

  /*if (center)
    {
    double m = centerbeta();
    fcconst->update_intercept(m);
    }*/
  if(randomslope && (includefixed || center))
    {
    workbeta = beta.getV();
    double s=0;
    for (i=0;i<n;i++,workbeta++)
      s += *workbeta;
    s /= double(n);

    if(includefixed)
      {
      double v = sigma2/double(nrpar-1);
      beta(nrpar-1,0) = s+sqrt(v)*rand_normal();
      s = beta(nrpar-1,0);
      }
    else
      {
      fcconst->update_fix_varcoeff(s,datanames[1]);
      }

    workbeta = beta.getV();
    for (i=0;i<n;i++,workbeta++)
      *workbeta -= s;
    }

  FULLCOND_random::update();
  update_spatialtotal();
  }


void FULLCOND_random_stepwise::update_spatialtotal(void)
  {
  if (spatialtotal)
    {
    double * ftotal_bp = ftotal.getbetapointer();
    double * workbeta=beta.getV();
    double * workbetaspat;
    //if (nongaussian)
    //  workbetaspat = fnonp->getbetapointer();
    //else
      workbetaspat = fbasisp->getbetapointer();

    int * indexp = indextotal.getV();
    unsigned i;
    for (i=0;i<nrpar;i++,workbeta++,ftotal_bp++,indexp++)
      {
      workbetaspat+= *indexp;
      *ftotal_bp = *workbeta + *workbetaspat;
      }

    ftotal.set_transform(likep->get_trmult(column));
    ftotal.update();
    }
  }


void FULLCOND_random_stepwise::change_Korder(double lamb)
  {
  set_lambdaconst(1000000000);
  }

void FULLCOND_random_stepwise::undo_Korder(void)
  {
  }


void FULLCOND_random_stepwise::init_spatialtotal(FULLCOND_nonp_basis * sp,
                                        const ST::string & pnt,
                                        const ST::string & prt)
  {
  df_unstruct = 0;
  fbasisp = sp;
  vector<ST::string> ev = sp->get_effectvalues();

  FULLCOND_random::init_spatialtotal(ev,pnt,prt);

  }


} // end: namespace MCMC




