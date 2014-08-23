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



#include "randomeffect.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------- class FULLCOND_random --------------------------------------
//------------------------------------------------------------------------------

// FOR STEPWISE

void FULLCOND_random::compute_lambdavec(vector<double> & lvec,unsigned & number)
  {
  FULLCOND::compute_lambdavec(lvec,number);
  if (randomslope)
    lvec.push_back(-1);

  lvec.push_back(0);
  }


const datamatrix & FULLCOND_random::get_data_forfixedeffects(void)
  {
  // useful for randomslopes only
  if (data_forfixed.rows() < data.rows())
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


// END: FOR STEPWISE

void FULLCOND_random::init_name(const ST::string & na)
    {
    char charh = '_';
    ST::string stringh = "\\_";

    FULLCOND::init_name(na);

    ST::string helpname = na.insert_string_char(charh,stringh);
    term_symbolic = "f_{" +  helpname + "}("+helpname+")";

    if (column > 0)
      priorassumptions.push_back("$" + term_symbolic + "$" +
       " (" + ST::inttostring(column+1) + ". response category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");

    priorassumptions.push_back("i.i.d. Gaussian random effects");
    }


void FULLCOND_random::init_names(const vector<ST::string> & na)
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
      ST::string helpname1 = na[0].insert_string_char(charh,stringh);
      ST::string helpname2 = na[1].insert_string_char(charh,stringh);
      term_symbolic = "f_{" +  helpname1 + "}("+helpname1+") \\cdot "
                        + helpname2;
      }

    if (column > 0)
      priorassumptions.push_back("$" + term_symbolic + "$" +
       " (" + ST::inttostring(column+1) + ". response category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");
    priorassumptions.push_back(term_symbolic);

    priorassumptions.push_back("i.i.d. Gaussian random effects");
    }


void FULLCOND_random::init_spatialtotal(vector<ST::string> & ev,
                                        const ST::string & pnt,
                                        const ST::string & prt)
{

spatialtotal = true;
pathresult3 = prt;
pathcurrent3 = prt;
pathsample_total = pnt;

ftotal =
FULLCOND(optionsp,datamatrix(1,1),"spatialtotal",effvalues.rows(),1,pnt);
ftotal.setflags(MCMC::norelchange | MCMC::nooutput);

unsigned i;
datamatrix evdouble(ev.size(),1);
statmatrix<int> index_ev(ev.size(),1);
index_ev.indexinit();
double v;
int h;
double * evdoublep = evdouble.getV();
vector<ST::string>::iterator it = ev.begin();
for(i=0;i<evdouble.rows();i++,evdoublep++,++it)
  {
  h = (*it).strtodouble(v);
  *evdoublep = v;
  }

evdouble.indexsort(index_ev,0,evdouble.rows()-1,0,0);

indextotal = statmatrix<int>(effvalues.rows(),1);

double * effvp = effvalues.getV();
int * indexp = indextotal.getV();
int indexalt = 0;
int * indexp2 = index_ev.getV();
for(i=0;i<ev.size();i++,indexp2++)
  {

  if (evdouble(*indexp2,0)== *effvp)
    {
    *indexp = *indexp2-indexalt;
    indexalt = *indexp2;
    effvp++;
    indexp++;
    }

  }

}


double FULLCOND_random::centerbeta(void)
  {
  unsigned i;

  double sum=0;
  double * workbeta = beta.getV();

  for (i=0;i<nrpar;i++,workbeta++)
    sum+= *workbeta;

  sum /= nrpar;


  double v = sigma2/double(nrpar);

  sum = sum+sqrt(v)*rand_normal();


  workbeta = beta.getV();

  for (i=0;i<nrpar;i++,workbeta++)
    *workbeta-= sum;

  return sum;

  }


void FULLCOND_random::compute_XWX(const datamatrix & weightmat,
                                  const unsigned & col)
  {

  register unsigned j,i;

  double * workXX = XX.getV();
  int *  workindex = index.getV();
  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  unsigned n = posbeg.size();

  if (!randomslope)
    {
    for(j=0;j<n;j++,workXX++,++itbeg,++itend)
      {
      *workXX = 0;
      for (i=*itbeg;i<=*itend;i++,workindex++)
        *workXX += weightmat(*workindex,col);
      }
    }
  else
    {

    double * datap = data.getV();
    for(j=0;j<n;j++,workXX++,++itbeg,++itend)
      {
      *workXX = 0;
      for (i=*itbeg;i<=*itend;i++,workindex++,datap++)
        {
        *workXX += weightmat(*workindex,col) * (*datap) * (*datap);
        }
      }
    }

  }


double FULLCOND_random::compute_quadform(void)
  {

  unsigned n;
  double sum = 0;
  double * workbeta = beta.getV();
  register unsigned i;

  if (randomslope && includefixed)
    n = nrpar-1;
  else
    n = nrpar;

  for(i=0;i<n;i++,workbeta++)
    sum += *workbeta * *workbeta;

  return sum;

  }


unsigned FULLCOND_random::get_rankK(void)
  {
  if (randomslope && includefixed)
    return nrpar-1;
  else
    return nrpar;
  }


void FULLCOND_random::set_lambdaconst(double la)
  {
  lambda=la;
  lambdaconst = true;
  }


// REML: random slope
FULLCOND_random::FULLCOND_random(MCMCoptions * o,
                  const datamatrix & intvar,const datamatrix & effmod,
                  const ST::string & t,const ST::string & fp,
                  const ST::string & pr,const double & la)
                  : FULLCOND(o,t)
  {

  fctype = randomslopes;

  spatialtotal=false;
  randomslope = true;

  pathresult = pr;
  pathcurrent = pr;

  lambda = la;
  lambdaconst=false;

  index = statmatrix<int>(effmod.rows(),1);
  index2 = statmatrix<int>(effmod.rows(),1);
  index.indexinit();
  effmod.indexsort(index,0,effmod.rows()-1,0,0);

  data = intvar;

  posbeg = vector<unsigned>();
  posend = vector<unsigned>();

  posbeg.push_back(0);
  int * workindex=index.getV()+1;
  int help = index(0,0);
  unsigned j;
  for(j=1;j<effmod.rows();j++,workindex++)
    {
    if ( effmod(*workindex,0) != effmod(help,0) )
      {
      posbeg.push_back(j);
      posend.push_back(j-1);
      }

    help = *workindex;
    }

  posend.push_back(effmod.rows()-1);

  effvalues = datamatrix(posbeg.size(),1);
  double * workeffvalues = effvalues.getV();
  for(j=0;j<posbeg.size();j++,workeffvalues++)
    *workeffvalues = effmod(index(posbeg[j],0),0);

  dimX = 0;
  dimZ = posbeg.size();

  nrpar = posbeg.size();

  }


// REML: random intercept
FULLCOND_random::FULLCOND_random(MCMCoptions * op,const datamatrix & d,
                                 const ST::string & t,const ST::string & fp,
                                 const ST::string & pr,const double & la)
                            : FULLCOND(op,t)
  {

  fctype = randomeffects;

  spatialtotal=false;
  randomslope = false;
  includefixed = false;

  pathresult = pr;
  pathcurrent = pr;

  lambda = la;
  lambdaconst=false;

  index = statmatrix<int>(d.rows(),1);
  index2 = statmatrix<int>(d.rows(),1);
  index.indexinit();
  d.indexsort(index,0,d.rows()-1,0,0);

  unsigned j;
  int * workindex = index.getV();
  int * workindex2 = index2.getV();
  *workindex2 = *workindex;
  int help = *workindex;
  workindex++;
  workindex2++;

  for(j=1;j<d.rows();j++,workindex++,workindex2++)
    {
    *workindex2 = *workindex-help;
    help = *workindex;
    }

  posbeg = vector<unsigned>();
  posend = vector<unsigned>();

  posbeg.push_back(0);
  workindex=index.getV()+1;
  help = index(0,0);
  for(j=1;j<d.rows();j++,workindex++)
    {
    if ( d(*workindex,0) != d(help,0) )
      {
      posbeg.push_back(j);
      posend.push_back(j-1);
      }

    help = *workindex;

    }

  posend.push_back(d.rows()-1);

  effvalues = datamatrix(posbeg.size(),1);
  double * workeffvalues = effvalues.getV();
  for(j=0;j<posbeg.size();j++,workeffvalues++)
    *workeffvalues = d(index(posbeg[j],0),0);

  dimX = 0;
  dimZ = posbeg.size();

  nrpar = posbeg.size();

  }


void FULLCOND_random::createreml(datamatrix & X,datamatrix & Z,
                                 const unsigned & Xpos, const unsigned & Zpos)
  {
  unsigned i,j;

  if (randomslope)
    {

    for(i=0;i<posbeg.size();i++)
      {
      for (j=posbeg[i];j<=posend[i];j++)
        Z(index(j,0),Zpos+i) = data(index(j,0),0);
      }

    }
  else
    {
    for(i=0;i<posbeg.size();i++)
      {
      for (j=posbeg[i];j<=posend[i];j++)
        Z(index(j,0),Zpos+i) = 1;
      }
    }

  }


void FULLCOND_random::outresultsreml(datamatrix & X,datamatrix & Z,
                                     datamatrix & betareml,datamatrix & betacov,
                                     datamatrix & thetareml,
                                     const unsigned & Xpos,
                                     const unsigned & Zpos,
                                     const unsigned & fcpos,
                                     const bool & dispers)
  {

  betamean=datamatrix(nrpar,1,0);
  betavar=datamatrix(nrpar,1,0);
  datamatrix betastd=datamatrix(nrpar,1,0);
  betaqu_l1_lower=datamatrix(nrpar,1,0);
  betaqu_l1_upper=datamatrix(nrpar,1,0);
  betaqu_l2_lower=datamatrix(nrpar,1,0);
  betaqu_l2_upper=datamatrix(nrpar,1,0);
  betaqu50=datamatrix(nrpar,1,0);
//  datamatrix betapval=datamatrix(nrpar,1,0);

  unsigned i;
  for(i=0;i<nrpar;i++)
    {
    betamean(i,0) = betareml(X.cols()+Zpos+i,0);
    betastd(i,0) = sqrt(betacov(X.cols()+Zpos+i,X.cols()+Zpos+i));
    betaqu_l1_lower(i,0) = betamean(i,0)+randnumbers::invPhi2(lower1/100)*betastd(i,0);
    betaqu_l1_upper(i,0) = betamean(i,0)+randnumbers::invPhi2(upper2/100)*betastd(i,0);
    betaqu_l2_lower(i,0) = betamean(i,0)+randnumbers::invPhi2(lower2/100)*betastd(i,0);
    betaqu_l2_upper(i,0) = betamean(i,0)+randnumbers::invPhi2(upper1/100)*betastd(i,0);
/*    if (betamean(i,0)>0)
      {
      betapval(i,0) = 1-randnumbers::Phi2(betamean(i,0)/betastd(i,0));
      }
    else
      {
      betapval(i,0) = randnumbers::Phi2(betamean(i,0)/betastd(i,0));
      }*/
    }

  FULLCOND::outresults();

  ST::string vstr;

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);

  ST::string nl1 = l1;
  ST::string nl2 = l2;
  ST::string nu1 = u1;
  ST::string nu2 = u2;
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');

  optionsp->out("  Estimated variance: "
                + ST::doubletostring(thetareml(fcpos-1,0),6) + "\n");
  double smoothpar;
  if(dispers==true)
    {
    smoothpar = thetareml(thetareml.rows()-1,0)/thetareml(fcpos-1,0);
    optionsp->out("  Estimated smoothing parameter: "
                + ST::doubletostring(smoothpar,6) + "\n");
    }
  else
    {
    smoothpar = 1/thetareml(fcpos-1,0);
    optionsp->out("  Estimated smoothing parameter: "
                + ST::doubletostring(smoothpar,6) + "\n");
    }
  if(thetareml(fcpos-1,1)==1)
    {
    optionsp->out("  NOTE: Estimation of the variance was stopped after iteration "
                  + ST::doubletostring(thetareml(fcpos-1,2),0) + "\n");
    optionsp->out("        because the corresponding penalized part was small relative to the linear predictor.");
    }
  optionsp->out("\n");
  optionsp->out("  Variance and smoothing parameter are stored in file\n");
  optionsp->out("  " + pathcurrent.substr(0,pathcurrent.length()-4) + "_var.res\n");

  ofstream outvarres((pathcurrent.substr(0,pathcurrent.length()-4) + "_var.res").strtochar());
  outvarres << "variance  ";
  outvarres << "smoothpar  ";
  outvarres << "stopped  " <<endl;

  outvarres << thetareml(fcpos-1,0) <<"  ";
  outvarres << smoothpar <<"  ";
  outvarres << (thetareml(fcpos-1,1)==1);
  outvarres << endl;
  outvarres.close();

  optionsp->out("\n");
  if (randomslope)
    {
    optionsp->out("  Results for random slopes are stored in file\n");
    optionsp->out("  " + pathcurrent + "\n");
    }
  else
    {
    optionsp->out("  Results for random effects are stored in file\n");
    optionsp->out("  " + pathcurrent + "\n");
    }

  optionsp->out("\n");

  ofstream outres(pathcurrent.strtochar());
  assert(!outres.fail());

  ST::string name = datanames[0];

  outres << "intnr" << "   ";
  outres << name << "   ";
  outres << "pmode   ";
//  outres << "pvalue    ";
  outres << "ci"  << level1  << "lower   ";
  outres << "ci"  << level2  << "lower   ";
  outres << "std   ";
  outres << "ci"  << level2  << "upper   ";
  outres << "ci"  << level1  << "upper   ";
  outres << "pcat" << level1 << "   ";
  outres << "pcat" << level2 << "   ";
  outres << endl;

  double * workmean = betamean.getV();
  double * workstd = betastd.getV();
  double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
  double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
  double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
  double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
//  double * workbetapval = betapval.getV();

  for(i=0;i<nrpar;i++,workmean++,workstd++,workbetaqu_l1_lower_p++,
      workbetaqu_l2_lower_p++,workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
//      workbetapval++)
    {
    outres << (i+1) << "   ";
    outres << effvalues(i,0) << "   ";
    outres << *workmean << "   ";
//    outres << *workbetapval << "   ";
    outres << *workbetaqu_l1_lower_p << "   ";
    outres << *workbetaqu_l2_lower_p << "   ";
    outres << *workstd << "   ";
    outres << *workbetaqu_l2_upper_p << "   ";
    outres << *workbetaqu_l1_upper_p << "   ";

    if (*workbetaqu_l1_lower_p > 0)
      outres << "1   ";
    else if (*workbetaqu_l1_upper_p < 0)
      outres << "-1   ";
    else
      outres << "0   ";

    if (*workbetaqu_l2_lower_p > 0)
      outres << "1   ";
    else if (*workbetaqu_l2_upper_p < 0)
      outres << "-1   ";
    else
      outres << "0   ";

    outres << endl;
    }
  }

// end reml


FULLCOND_random::FULLCOND_random (MCMCoptions * o,DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & d, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,
                              const double & la, const unsigned & c)
                            : FULLCOND(o,datamatrix(1,1),t,1,1,fp)
  {

  fcconst = fcc;

  fctype = randomeffects;

  spatialtotal=false;
  randomslope = false;
  includefixed = false;

  changingweight = dp->get_changingweight();

  column = c;

  likep = dp;
  pathresult = pr;
  pathcurrent = pr;

  lambda = la;
  lambdaold1 = -1;
  lambdaold2 = -1;
  lambdaconst=false;

  index = statmatrix<int>(d.rows(),1);
  index2 = statmatrix<int>(d.rows(),1);
  index.indexinit();
  d.indexsort(index,0,d.rows()-1,0,0);

  unsigned j;
  int * workindex = index.getV();
  int * workindex2 = index2.getV();
  *workindex2 = *workindex;
  int help = *workindex;
  workindex++;
  workindex2++;

  for(j=1;j<d.rows();j++,workindex++,workindex2++)
    {
    *workindex2 = *workindex-help;
    help = *workindex;
    }

  posbeg = vector<unsigned>();
  posend = vector<unsigned>();

  posbeg.push_back(0);
  workindex=index.getV()+1;
  help = index(0,0);
  for(j=1;j<d.rows();j++,workindex++)
    {
    if ( d(*workindex,0) != d(help,0) )
      {
      posbeg.push_back(j);
      posend.push_back(j-1);
      }

    help = *workindex;

    }

  posend.push_back(d.rows()-1);

  effvalues = datamatrix(posbeg.size(),1);
  double * workeffvalues = effvalues.getV();
  for(j=0;j<posbeg.size();j++,workeffvalues++)
    *workeffvalues = d(index(posbeg[j],0),0);

  XX = datamatrix(posbeg.size());
  compute_XWX(likep->get_weight(),0);


  setbeta(posbeg.size(),1,0);

  identifiable = false;

  muy = datamatrix(nrpar,1);

//  identifiable =true;

  }


// randomslope
FULLCOND_random::FULLCOND_random(MCMCoptions * o,DISTRIBUTION * dp,
                  FULLCOND_const * fcc,
                  const datamatrix & intvar,const datamatrix & effmod,
                  const ST::string & t,
                  const ST::string & fp,const ST::string & pr,
                  const ST::string & prf,
                  const double & la,
                  const bool & inclfixed,const unsigned & c)
                  : FULLCOND(o,datamatrix(1,1),t,1,1,fp)
  {


  fcconst = fcc;

  fctype = randomslopes;

  spatialtotal=false;
  randomslope = true;
  includefixed = inclfixed;

  changingweight = dp->get_changingweight();

  column = c;

  likep = dp;
  pathresult = pr;
  pathcurrent = pr;
  pathresult2 = prf;
  pathcurrent2 = prf;


  lambda = la;
  lambdaold1 = -1;
  lambdaold2 = -1;
  lambdaconst=false;

  index = statmatrix<int>(effmod.rows(),1);
  index2 = statmatrix<int>(effmod.rows(),1);
  index.indexinit();
  effmod.indexsort(index,0,effmod.rows()-1,0,0);

  unsigned j;
  int * workindex = index.getV();
  int * workindex2 = index2.getV();
  *workindex2 = *workindex;
  int help = *workindex;
  workindex++;
  workindex2++;

  for(j=1;j<effmod.rows();j++,workindex++,workindex2++)
    {
    *workindex2 = *workindex-help;
    help = *workindex;
    }


  data = datamatrix(effmod.rows(),1);
  data2 = datamatrix(effmod.rows(),1);
  double * datap = data.getV();
  double * datap2 = data2.getV();
  workindex = index.getV();
  for(j=0;j<effmod.rows();j++,datap++,datap2++,workindex++)
    {
    *datap = intvar(*workindex,0);
    *datap2 = (*datap) * (*datap);
    }

  posbeg = vector<unsigned>();
  posend = vector<unsigned>();

  posbeg.push_back(0);
  workindex=index.getV()+1;
  help = index(0,0);
  for(j=1;j<effmod.rows();j++,workindex++)
    {
    if ( effmod(*workindex,0) != effmod(help,0) )
      {
      posbeg.push_back(j);
      posend.push_back(j-1);
      }

    help = *workindex;
    }

  posend.push_back(effmod.rows()-1);

  effvalues = datamatrix(posbeg.size(),1);
  double * workeffvalues = effvalues.getV();
  for(j=0;j<posbeg.size();j++,workeffvalues++)
    *workeffvalues = effmod(index(posbeg[j],0),0);

  XX = datamatrix(posbeg.size());
  compute_XWX(likep->get_weight(),0);

  if (includefixed)
    setbeta(posbeg.size()+1,1,0);
  else
    setbeta(posbeg.size(),1,0);

  muy = datamatrix(nrpar,1);

  identifiable = true;

  }


FULLCOND_random::FULLCOND_random(const FULLCOND_random & fc)
                            : FULLCOND(FULLCOND(fc))
  {
  muy = fc.muy;
  fcconst = fc.fcconst;
  randomslope = fc.randomslope;
  includefixed = fc.includefixed;
  XX = fc.XX;
  likep = fc.likep;
  index = fc.index;
  index2 = fc.index2;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effvalues = fc.effvalues;
  sigma2 = fc.sigma2;
  pathsample_total = fc.pathsample_total;
  ftotal = fc.ftotal;
  spatialtotal = fc.spatialtotal;
  lambda = fc.lambda;
  lambdaold1 = fc.lambdaold1;
  lambdaold2 = fc.lambdaold2;
  df_lambdaold1 = fc.df_lambdaold1;
  df_lambdaold2 = fc.df_lambdaold2;
  lambdaconst=fc.lambdaconst;
  data2 = fc.data2;
  }


const FULLCOND_random & FULLCOND_random::
         operator=(const FULLCOND_random & fc)
  {
  if (this==&fc)
    return *this;

  FULLCOND::operator=(FULLCOND(fc));

  muy = fc.muy;
  fcconst = fc.fcconst;
  randomslope = fc.randomslope;
  includefixed = fc.includefixed;
  XX = fc.XX;
  likep = fc.likep;
  index = fc.index;
  index2 = fc.index2;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effvalues = fc.effvalues;
  sigma2 = fc.sigma2;
  pathsample_total = fc.pathsample_total;
  ftotal = fc.ftotal;
  spatialtotal = fc.spatialtotal;
  lambda = fc.lambda;
  lambdaold1 = fc.lambdaold1;
  lambdaold2 = fc.lambdaold2;
  df_lambdaold1 = fc.df_lambdaold1;
  df_lambdaold2 = fc.df_lambdaold2;
  lambdaconst=fc.lambdaconst;
  data2 = fc.data2;

  return *this;
  }


void FULLCOND_random::update(void)
  {

  transform = likep->get_trmult(column);

  FULLCOND::update();

  }


void FULLCOND_random::outresults(void)
  {
  FULLCOND::outresults();

  ST::string vstr;

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);

  ST::string nl1 = l1;
  ST::string nl2 = l2;
  ST::string nu1 = u1;
  ST::string nu2 = u2;
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');

  if (randomslope && includefixed)
    {
    optionsp->out("  Fixed effect:\n");
    optionsp->out("\n");

    ST::string help =  ST::doubletostring(lower1,4) + "% quant.";
    ST::string levell = help + ST::string(' ',15-help.length());
    help = ST::doubletostring(upper2,4) + "% quant.";
    ST::string levelu = help + ST::string(' ',15-help.length());
    help = ST::string(' ',+2);

    optionsp->out(help +
                    "mean           " +
                    "Std. Dev.      " +
                    levell +
                    "median         " +
                    levelu + "\n");

    optionsp->out(ST::outresults(0,"",betamean(nrpar-1,0),
                  sqrt(betavar(nrpar-1,0)),betaqu_l1_lower(nrpar-1,0),
                  betaqu50(nrpar-1,0),betaqu_l1_upper(nrpar-1,0)));

    optionsp->out("\n");

    optionsp->out("  Results for the fixed effect are also stored in file \n");
    optionsp->out("  " + pathcurrent2 + "\n");

    optionsp->out("\n");

    ofstream outfixed(pathcurrent2.strtochar());

    outfixed << "pmean   ";
    outfixed << "pqu"  << nl1  << "   ";
    outfixed << "pqu"  << nl2  << "   ";
    outfixed << "pmed   ";
    outfixed << "pqu"  << nu1  << "   ";
    outfixed << "pqu"  << nu2  << "   ";
    outfixed << "pcat" << level1 << "   ";
    outfixed << "pcat" << level2 << "   ";
    outfixed << endl;

    outfixed << betamean(nrpar-1,0) << "   ";
    outfixed << betaqu_l1_lower(nrpar-1,0) << "   ";
    outfixed << betaqu_l2_lower(nrpar-1,0) << "   ";
    outfixed << betaqu50(nrpar-1,0) << "   ";
    outfixed << betaqu_l2_upper(nrpar-1,0) << "   ";
    outfixed << betaqu_l1_upper(nrpar-1,0) << "   ";
    if (betaqu_l1_lower(nrpar-1,0) > 0)
      outfixed << "1   ";
    else if (betaqu_l1_upper(nrpar-1,0) < 0)
      outfixed << "-1   ";
    else
      outfixed << "0   ";
    if (betaqu_l2_lower(nrpar-1,0) > 0)
      outfixed << "1   ";
    else if (betaqu_l2_upper(nrpar-1,0) < 0)
      outfixed << "-1   ";
    else
      outfixed << "0   ";

    outfixed << endl;
    }

  if (randomslope)
    {
    optionsp->out("  Results for random slopes are stored in file\n");
    optionsp->out("  " + pathcurrent + "\n");
    }
  else
    {
    optionsp->out("  Results for random effects are stored in file\n");
    optionsp->out("  " + pathcurrent + "\n");
    }

  if (lambdaconst==true)
    {
    optionsp->out("\n");
    optionsp->out("  Constant smoothing parameter: " +
    ST::doubletostring(lambda,6) + "\n");
    optionsp->out("\n");
    }

  if (optionsp->get_samplesize() == 0)
    {
    optionsp->out("\n");
    double df = compute_df();
    optionsp->out("  Approximate degrees of freedom: "
                    + ST::doubletostring(df,6) + "\n");
    }

  optionsp->out("\n");

  unsigned i;

  ofstream outres(pathcurrent.strtochar());
  assert(!outres.fail());

  ST::string name = datanames[0];

  outres << "intnr" << "   ";
  outres << name << "   ";
  outres << "pmean   ";
  outres << "pqu"  << nl1  << "   ";
  outres << "pqu"  << nl2  << "   ";
  outres << "pmed   ";
  outres << "pqu"  << nu1  << "   ";
  outres << "pqu"  << nu2  << "   ";
  outres << "pcat" << level1 << "   ";
  outres << "pcat" << level2 << "   ";
  outres << endl;

  double * workmean = betamean.getV();
  double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
  double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
  double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
  double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
  double * workbetaqu50 = betaqu50.getV();

  for(i=0;i<nrpar;i++,workmean++,workbetaqu_l1_lower_p++,
      workbetaqu_l2_lower_p++,workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
      workbetaqu50++)
    {
    if (randomslope && includefixed && i == nrpar-1)
      {
      }
    else
      {
      outres << (i+1) << "   ";
      outres << effvalues(i,0) << "   ";
      outres << *workmean << "   ";
      outres << *workbetaqu_l1_lower_p << "   ";
      outres << *workbetaqu_l2_lower_p << "   ";
      outres << *workbetaqu50 << "   ";
      outres << *workbetaqu_l2_upper_p << "   ";
      outres << *workbetaqu_l1_upper_p << "   ";

      if (*workbetaqu_l1_lower_p > 0)
        outres << "1   ";
      else if (*workbetaqu_l1_upper_p < 0)
        outres << "-1   ";
      else
        outres << "0   ";

      if (*workbetaqu_l2_lower_p > 0)
        outres << "1   ";
      else if (*workbetaqu_l2_upper_p < 0)
        outres << "-1   ";
      else
        outres << "0   ";

      outres << endl;
      }
    }

  if (spatialtotal)
    {

    if (transformnonlinear)
      {
      ST::string suffix="";
      vector<FULLCOND *> fcp(1);
      fcp[0] = &ftotal;
      likep->transform_nonlinear(fcp,transformtype);
      ftotal.set_transform(suffix,transformtype);
      }

    ftotal.outresults();

    optionsp->out("  Results for the sum of the structured and unstructured \n");
    optionsp->out("  spatial effects are stored in file \n");
    optionsp->out("  " + pathcurrent3 + "\n");
    optionsp->out("\n");

    ofstream outres2(pathcurrent3.strtochar());
    assert(!outres2.fail());

    outres2 << "intnr" << "   ";
    outres2 << name << "   ";
    outres2 << "pmean   ";
    outres2 << "pqu"  << nl1  << "   ";
    outres2 << "pqu"  << nl2  << "   ";
    outres2 << "pmed   ";
    outres2 << "pqu"  << nu1  << "   ";
    outres2 << "pqu"  << nu2  << "   ";
    outres2 << "pcat" << level1 << "   ";
    outres2 << "pcat" << level2 << "   ";
    outres2 << endl;

    workmean = ftotal.get_betameanp();
    workbetaqu_l1_lower_p = ftotal.get_beta_lower1_p();
    workbetaqu_l2_lower_p = ftotal.get_beta_lower2_p();
    workbetaqu_l1_upper_p = ftotal.get_beta_upper1_p();
    workbetaqu_l2_upper_p = ftotal.get_beta_upper2_p();
    workbetaqu50 = ftotal.get_betaqu50p();

    for(i=0;i<nrpar;i++,workmean++,workbetaqu_l1_lower_p++,
        workbetaqu_l2_lower_p++,workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
        workbetaqu50++)
      {
      outres2 << (i+1) << "   ";
      outres2 << effvalues(i,0) << "   ";
      outres2 << *workmean << "   ";
      outres2 << *workbetaqu_l1_lower_p << "   ";
      outres2 << *workbetaqu_l2_lower_p << "   ";
      outres2 << *workbetaqu50 << "   ";
      outres2 << *workbetaqu_l2_upper_p << "   ";
      outres2 << *workbetaqu_l1_upper_p << "   ";

      if (*workbetaqu_l1_lower_p > 0)
        outres2 << "1   ";
      else if (*workbetaqu_l1_upper_p < 0)
        outres2 << "-1   ";
      else
        outres2 << "0   ";

      if (*workbetaqu_l2_lower_p > 0)
        outres2 << "1   ";
      else if (*workbetaqu_l2_upper_p < 0)
        outres2 << "-1   ";
      else
        outres2 << "0   ";

      outres2 << endl;
      }

    }

  }


void FULLCOND_random::outoptions(void)
  {

  optionsp->out("  OPTIONS FOR RANDOM EFFECT: " + title + "\n",true);
  optionsp->out("\n");
  if (lambdaconst==true)
    {
    optionsp->out("  Constant smoothing parameter: " +
    ST::doubletostring(lambda,6) + "\n");
    optionsp->out("\n");
    }

  }


void FULLCOND_random::update_linpred_diff(datamatrix & b1,datamatrix & b2)
  {

  unsigned n=nrpar;
  if (includefixed)
    n = nrpar-1;

  int * workindex;
  double * workb1 = b1.getV();
  double * workb2 = b2.getV();
  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();

  unsigned i,j;

  if (randomslope)
    {

    workindex = index.getV();
    double * workdata=data.getV();

    if (includefixed)
      {

      n = nrpar-1;
      double ms1 = b1(nrpar-1,0);
      double ms2 = b2(nrpar-1,0);
      double h;

      for (i=0;i<n;i++,workb1++,workb2++,++itbeg,++itend)
        {
        if (*itbeg != -1)
          {
          h =*workb1+ms1-*workb2-ms2;
          for(j=*itbeg;j<=*itend;j++,workindex++,workdata++)
            {
            likep->add_linearpred(h * (*workdata),
                                  unsigned(*workindex),column);
            }
          }
        }

      }
    else
      {
      for (i=0;i<n;i++,workb1++,workb2++,++itbeg,++itend)
        {
        if (*itbeg != -1)
          {
          for(j=*itbeg;j<=*itend;j++,workindex++,workdata++)
            {
            likep->add_linearpred((*workb1-*workb2) * (*workdata),
                                  unsigned(*workindex),column);
            }

          }

        }

      }

    }
  else
    {
    for (i=0;i<n;i++,workb1++,workb2++,++itbeg,++itend)
      {
      if (*itbeg != -1)
        likep->add_linearpred(*workb1-*workb2,*itbeg,*itend,index,column);
      }
    }

  }



void FULLCOND_random::update_linpred(const bool & add)
  {
  unsigned i,j;

  unsigned n = nrpar;

  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  double * workbeta = beta.getV();

  int * workindex2;

  if (add==false)
    {
    if (!randomslope)
      {
      likep->set_linpredp_current(column);
      for (i=0;i<nrpar;i++,++itbeg,++itend,workbeta++)
        if (*itbeg != -1)
          likep->add_linearpred2(-(*workbeta),*itbeg,*itend,index2,column);
      }
    else
      {
      workindex2 = index2.getV();
      double * datap = data.getV();
      if (includefixed)
        {
        n = nrpar-1;
        double ms = beta(nrpar-1,0);
        double h;
        likep->set_linpredp_current(column);
        for (i=0;i<n;i++,++itbeg,++itend,workbeta++)
          {
          if (*itbeg != -1)
            {
            h = *workbeta+ms;
            for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
              likep->add_linearpred2(-h*(*datap),*workindex2);
            }
          }
        }
      else
        {
        n = nrpar;
        likep->set_linpredp_current(column);
        for (i=0;i<n;i++,++itbeg,++itend,workbeta++)
          {
          if (*itbeg != -1)
            {
            for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
              {
              likep->add_linearpred2(-*workbeta*(*datap),*workindex2);
              }
            }

          }
        }

      }


    } // end: if (add == false)
  else
    {

    if (!randomslope)
      {
      likep->set_linpredp_current(column);
      for (i=0;i<nrpar;i++,++itbeg,++itend,workbeta++)
        if (*itbeg != -1)
          likep->add_linearpred2(*workbeta,*itbeg,*itend,index2,column);
      }
    else
      {
      workindex2 = index2.getV();
      double * datap = data.getV();
      if (includefixed)
        {
        n = nrpar-1;
        double ms = beta(nrpar-1,0);
        double h;
        likep->set_linpredp_current(column);
        for (i=0;i<n;i++,++itbeg,++itend,workbeta++)
          {
          if (*itbeg != -1)
            {
            h = *workbeta+ms;
            for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
              likep->add_linearpred2(h*(*datap),*workindex2);
            }
          }
        }
      else
        {
        n = nrpar;
        likep->set_linpredp_current(column);
        for (i=0;i<n;i++,++itbeg,++itend,workbeta++)
          {
          if (*itbeg != -1)
            {
            for(j=*itbeg;j<=*itend;j++,workindex2++,datap++)
              {
              likep->add_linearpred2(*workbeta*(*datap),*workindex2);
              }
            }
          }
        }

      }

    }

  }


bool FULLCOND_random::posteriormode(void)
  {


  unsigned n = nrpar;
  if (includefixed)
    n = nrpar-1;

  update_linpred(false);

  compute_XWX(likep->get_weightiwls(),column);

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
    if (includefixed)
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


  if (randomslope && includefixed)
    {
    double * workbeta = beta.getV();
    double sum=0;
    for (i=0;i<n;i++,workbeta++)
      {
      sum += *workbeta;
      }

    beta(nrpar-1,0) = sum/double(n);

    workbeta = beta.getV();
    double ms = beta(nrpar-1,0);
    for (i=0;i<n;i++,workbeta++)
      *workbeta -= ms;

    }

  update_linpred(true);

  transform = likep->get_trmult(column);

  return FULLCOND::posteriormode();

  }


//------------------------------------------------------------------------------
//------------------- class FULLCOND_random_gaussian ---------------------------
//------------------------------------------------------------------------------

void FULLCOND_random_gaussian::init_spatialtotal(FULLCOND_nonp_basis * sp,
                                        const ST::string & pnt,
                                        const ST::string & prt)
{

fbasisp = sp;
vector<ST::string> ev = sp->get_effectvalues();

FULLCOND_random::init_spatialtotal(ev,pnt,prt);

}


void FULLCOND_random_gaussian::update(void)
  {

  double var;
  double m;
  unsigned i,j;
  unsigned n = nrpar;
  if (randomslope && includefixed)
  n = nrpar-1;


  if (optionsp->get_nriter()==1 || changingweight)
    compute_XWX(likep->get_weight(),0);


  if (lambdaconst == false)
    lambda = likep->get_scale(column)/sigma2;
  else
    sigma2 = likep->get_scale(column)/lambda;

  double sqrtscale = sqrt(likep->get_scale(column));


  update_linpred(false);


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
  for (i=0;i<n;i++,workbeta++,workmuy++,workXX++)
    {

    var = 1.0/(*workXX  + lambda);

    m = var * *workmuy;

    *workbeta = m + sqrtscale*sqrt(var)*rand_normal();

    }


  if (randomslope && includefixed)
    {

    workbeta = beta.getV();
    double s=0;
    for (i=0;i<nrpar-1;i++,workbeta++)
      s += *workbeta;
    s /= double(nrpar-1);

    double v = sigma2/double(nrpar-1);

    beta(nrpar-1,0) = s+sqrt(v)*rand_normal();

    workbeta = beta.getV();
    double ms = beta(nrpar-1,0);
    for (i=0;i<nrpar-1;i++,workbeta++)
      *workbeta -= ms;

    }


  update_linpred(true);


  if (center)
    {
    double m = centerbeta();
    fcconst->update_intercept(m);
    }


  acceptance++;

  transform = likep->get_trmult(column);


  FULLCOND_random::update();

  if (spatialtotal)
    {
    double * ftotal_bp = ftotal.getbetapointer();
    workbeta=beta.getV();
    double * workbetaspat = fbasisp->getbetapointer();
    int * indexp = indextotal.getV();
    for (i=0;i<nrpar;i++,workbeta++,ftotal_bp++,indexp++)
      {
      workbetaspat+= *indexp;
      *ftotal_bp = *workbeta + *workbetaspat;
      }

    ftotal.set_transform(likep->get_trmult(column));

    ftotal.update();
    }

  }


double FULLCOND_random_gaussian::compute_df(void)
  {
  unsigned i;
  double df=0;
  double * workXX=XX.getV();
  unsigned n;

  if (randomslope && includefixed)
    {
    n = nrpar-1;
    df=1;
    }
  else
    n = nrpar;


  if ((lambdaold1==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold1;
    }
  else if ((lambdaold2==lambda) && (likep->iwlsweights_constant() == true) )
    {
    df = df_lambdaold2;
    }
  else
    {
    for(i=0;i<n;i++,workXX++)
      {
      df += (*workXX)/(*workXX+lambda);
      }

    df_lambdaold2 = df_lambdaold1;
    lambdaold2 = lambdaold1;
    df_lambdaold1 = df;
    lambdaold1 = lambda;

    }

  return df;
  }


void FULLCOND_random_gaussian::update_stepwise(double la)
  {
  lambda = la;
  }


ST::string FULLCOND_random_gaussian::get_effect(void)
  {

  ST::string h;

  if(randomslope)
    h = datanames[0] + "*" + datanames[1];
  else
    h = datanames[0];

  h = h + "(random,lambda=" + ST::doubletostring(lambda,6) + ")";

  return h;

  }


void FULLCOND_random_gaussian::reset_effect(unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  }

//------------------------------------------------------------------------------
//----------------- class FULLCOND_random_nongaussian --------------------------
//------------------------------------------------------------------------------

void FULLCOND_random_nongaussian::init_spatialtotal(FULLCOND_nonp * sp,
                                        const ST::string & pnt,
                                        const ST::string & prt)
{

fnonp = sp;
nongaussian = true;
vector<ST::string> ev = sp->get_effectvalues();

FULLCOND_random::init_spatialtotal(ev,pnt,prt);

}

void FULLCOND_random_nongaussian::init_spatialtotal(FULLCOND_nonp_basis * sp,
                                        const ST::string & pnt,
                                        const ST::string & prt)
{

fbasisp = sp;
nongaussian=false;
vector<ST::string> ev = sp->get_effectvalues();

FULLCOND_random::init_spatialtotal(ev,pnt,prt);

}


FULLCOND_random_nongaussian::FULLCOND_random_nongaussian
                              (MCMCoptions * o,DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & d, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,
                              const double & la,
                              const unsigned & c)
                            : FULLCOND_random(o,dp,fcc,d,t,fp,pr,la,c)
  {
  }


FULLCOND_random_nongaussian::FULLCOND_random_nongaussian(MCMCoptions * o,
                              DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & intvar,
                              const datamatrix & effmod, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,
                              const ST::string & prf,
                              const double & la,
                              const bool & inclfixed, const unsigned & c)
                            : FULLCOND_random(o,dp,fcc,intvar,effmod,t,
                                              fp,pr,prf,la,inclfixed,c)
  {
  w = datamatrix(effmod.rows(),1,0);
  }


FULLCOND_random_nongaussian::FULLCOND_random_nongaussian
                           (const FULLCOND_random_nongaussian & fc)
                            : FULLCOND_random(FULLCOND_random(fc))
  {
  w = fc.w;
  mode = fc.mode;
  modeold = fc.modeold;
  proposal = fc.proposal;
  weightiwls = fc.weightiwls;
  tildey = fc.tildey;
  fnonp = fc.fnonp;
  fbasisp = fc.fbasisp;
  nongaussian = fc.nongaussian;
  }


const FULLCOND_random_nongaussian & FULLCOND_random_nongaussian::
         operator=(const FULLCOND_random_nongaussian & fc)
  {
  if (this==&fc)
    return *this;
  FULLCOND_random::operator=(FULLCOND_random(fc));
  w = fc.w;
  mode = fc.mode;
  modeold = fc.modeold;
  proposal = fc.proposal;
  weightiwls = fc.weightiwls;
  tildey = fc.tildey;
  fnonp = fc.fnonp;
  fbasisp = fc.fbasisp;
  nongaussian = fc.nongaussian;
  return *this;
  }


void FULLCOND_random_nongaussian::update_spatialtotal(void)
  {

  if (spatialtotal)
    {
    double * ftotal_bp = ftotal.getbetapointer();
    double * workbeta=beta.getV();
    double * workbetaspat;
    if (nongaussian)
      workbetaspat = fnonp->getbetapointer();
    else
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


void FULLCOND_random_nongaussian::update_random_intercept(void)
  {
  unsigned i,j;
  double sumw;
  double sumy;
  double var;
  double u;

  double logold;
  double lognew=0;

  double qbetanew=0;
  double qbetaold=0;

  double diff;
  double help;

  nrtrials++;

  // Initializing variables

  if (optionsp->get_nriter()==1)
    {
    mode = beta;
    modeold = beta;
    weightiwls=datamatrix(likep->get_nrobs(),1);
    tildey=datamatrix(likep->get_nrobs(),1);
    proposal = datamatrix(beta.rows(),1);
    }

  if (lambdaconst == false)
    lambda=1.0/sigma2;
  else
    sigma2 = 1.0/lambda;

  // Compute log-likelihood for old beta

  logold = likep->loglikelihood();

  // Use mode in the predictor

  update_linpred_diff(mode,beta);

  // Compute weightiwls ans tildey

  likep->compute_IWLS_weight_tildey(weightiwls,tildey,column);

  int * workindex2 = index2.getV();
  double * workmode = mode.getV();
  double * workbeta = beta.getV();
  double * workproposal = proposal.getV();

  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();

  double * workweightiwls = weightiwls.getV()+(*workindex2);
  double * worktildey = tildey.getV()+(*workindex2);

  modeold.assign(mode);

  // Compute proposal

  for (i=0;i<nrpar;i++,workbeta++,workmode++,workproposal++,++itbeg,++itend)
    {

    sumw=0;
    sumy=0;
    for (j=*itbeg;j<=*itend;j++,workindex2++,
         workweightiwls+=*workindex2,worktildey+=*workindex2)
      {
      help = *workweightiwls;
      sumw += help;
      sumy += help*( *worktildey+(*workmode));
      }

    var = 1.0/(sumw + lambda);

    *workmode = var*sumy;

    *workproposal = *workmode + sqrt(var)*rand_normal();

    logold -= 0.5*(*workbeta)*(*workbeta)/sigma2;
    lognew -= 0.5*(*workproposal)*(*workproposal)/sigma2;

    diff = *workbeta-*workmode;
    qbetaold -= 0.5*diff*diff/var;

    diff = *workproposal-*workmode;
    qbetanew -= 0.5*diff*diff/var;

    }

  update_linpred_diff(proposal,modeold);

  lognew += likep->loglikelihood();

  u = log(uniform());

  if (u <= (lognew - logold  + qbetaold - qbetanew) )
    {
    acceptance++;
    beta.assign(proposal);

    if (center)
      {
      double m = centerbeta();
     fcconst->update_intercept(m);
      }

    }
  else
    {
    update_linpred_diff(beta,proposal);
    }


  FULLCOND_random::update();


  update_spatialtotal();  


  }


void FULLCOND_random_nongaussian::update_random_intercept_singleblock(void)
  {

  unsigned i;
  double sumw=0;
  double sumy;
  double var;
  double u;

  double logold;
  double lognew=0;

  double qbetanew=0;
  double qbetaold=0;

  double diff;

  // Initializing variables

  if (optionsp->get_nriter()==1)
    {
    mode = beta;
    modeold = beta;
    }

  if (lambdaconst == false)
    lambda=1.0/sigma2;
  else
    sigma2 = 1.0/lambda;


  double * workmode = mode.getV();
  double * workbeta = beta.getV();

  double proposal;

  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();


  modeold.assign(mode);
  double * workmodeold = modeold.getV();


  for (i=0;i<nrpar;i++,workbeta++,workmode++,workmodeold++,++itbeg,++itend)
    {

    nrtrials++;

    // Compute loglikelihood based on current beta
    logold = likep->loglikelihood2(*itbeg,*itend,index,index2);


    // adds the current mode to linpred and substracts the current beta
    if (*itbeg==0)
      likep->set_linpredp_current(column);
    else
      likep->set_linpredp_current(index(*itbeg-1,0),column);
    likep->add_linearpred2(*workmode-*workbeta,*itbeg,*itend,index2,column);



    // compute sum weightiwls  (corresponds to X'WX) and sum tildey

    sumy = likep->compute_sumweight_sumy(
    *workmode,sumw,*itbeg,*itend,index,index2,column);

    // compute proposal
    var = 1.0/(sumw + lambda);

    *workmode = var*sumy;

    proposal = *workmode + sqrt(var)*rand_normal();

    logold -= 0.5*(*workbeta)*(*workbeta)/sigma2;
    lognew = 0.5*(proposal)*(proposal)/sigma2;

    diff = *workbeta-*workmode;
    qbetaold = -0.5*diff*diff/var;

    diff = proposal-*workmode;
    qbetanew = -0.5*diff*diff/var;


    // adds the proposal to linpred and substracts the old mode
    if (*itbeg==0)
      likep->set_linpredp_current(column);
    else
      likep->set_linpredp_current(index(*itbeg-1,0),column);
    likep->add_linearpred2(proposal-*workmodeold,*itbeg,*itend,index2,column);


    // Compute loglikelihood based on proposal
    lognew += likep->loglikelihood2(*itbeg,*itend,index,index2);

    // accept/reject proposal
    u = log(uniform());

    if (u <= (lognew - logold  + qbetaold - qbetanew) )
      {
      acceptance++;
      *workbeta = proposal;

      }
    else
      {

      // adds the old beta  to linpred and substracts the proposal

      if (*itbeg==0)
        likep->set_linpredp_current(column);
      else
        likep->set_linpredp_current(index(*itbeg-1,0),column);

      likep->add_linearpred2(*workbeta - proposal,*itbeg,*itend,index2,column);

      }


    }


  // centers beta around zero
 if (center)
   {
   double m = centerbeta();
   fcconst->update_intercept(m);
   }


  FULLCOND_random::update();

  update_spatialtotal();


  }



void FULLCOND_random_nongaussian::update_random_slope_includefixed(void)
  {

    nrtrials++;

  // Initializing variables

  if (optionsp->get_nriter()==1)
    {
    mode = beta;
    modeold = beta;
    weightiwls=datamatrix(likep->get_nrobs(),1);
    tildey=datamatrix(likep->get_nrobs(),1);
    proposal = datamatrix(beta.rows(),1,0);
    }


  if (lambdaconst == false)
    lambda=1.0/sigma2;
  else
    sigma2 = 1.0/lambda;


  unsigned i,j;
  double sumw;
  double sumy;
  double var;

  double u;

  double logold=0;
  double lognew=0;

  double qbetanew=0;
  double qbetaold=0;

  double diff;

  unsigned n=nrpar-1;

  // Compute log-likelihood for old beta

  logold = likep->loglikelihood();

  // Use mode in the predictor

  update_linpred_diff(mode,beta);

  // Compute weightiwls and tildey

  likep->compute_IWLS_weight_tildey(weightiwls,tildey,column);


  int * workindex2 = index2.getV();

  double * workmode = mode.getV();
  double * workbeta = beta.getV();
  double * workproposal = proposal.getV();

  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();

  double * workweightiwls = weightiwls.getV()+(*workindex2);
  double * worktildey = tildey.getV()+(*workindex2);

  double * workdata2 = data2.getV();
  double * workdata = data.getV();

  modeold.assign(mode);

  for (i=0;i<n;i++,workbeta++,workmode++,workproposal++,++itbeg,++itend)
    {

    sumw=0;
    sumy=0;
    for (j=*itbeg;j<=*itend;j++,workindex2++,
         workweightiwls+=*workindex2,worktildey+=*workindex2,
         workdata2++,workdata++)
      {
      sumw += (*workweightiwls) * (*workdata2);
      sumy += (*workdata)*(*workweightiwls)*
                 ((*worktildey)+(*workdata)*(*workmode +mode(n,0)) );
      }


    sumy += mode(n,0)*lambda;

    var = 1.0/(sumw + lambda);

    *workmode = var*sumy;

    *workproposal = *workmode + sqrt(var)*rand_normal();


    logold -= 0.5*(*workbeta)*(*workbeta)*lambda;
    lognew -= 0.5*(*workproposal-beta(n,0))*(*workproposal)-beta(n,0)*lambda;


    diff = *workbeta+beta(n,0)-*workmode;
    qbetaold -= 0.5*diff*diff/var;

    diff = *workproposal-*workmode;
    qbetanew -= 0.5*diff*diff/var;

    } // end:   for (i=0;i<nrpar;i++,workbeta++)


  update_linpred_diff(proposal,modeold);

  lognew += likep->loglikelihood();

  u = log(uniform());

  if (u <= (lognew - logold  + qbetaold - qbetanew) )
    {
    acceptance++;
    beta.assign(proposal);
    }
  else
    {
    update_linpred_diff(beta,proposal);
    }


  workbeta = beta.getV();
  workmode = mode.getV();
  double s=0;
  double smode=0;
  for (i=0;i<nrpar-1;i++,workbeta++,workmode++)
    {
    s += *workbeta;
    smode += *workmode;
    }

  s /= double(nrpar-1);
  smode /= double(nrpar-1);

  double v = sigma2/double(nrpar-1);

  beta(n,0) = s+sqrt(v)*rand_normal();

  mode(n,0) = smode;

  workbeta = beta.getV();
  workmode = mode.getV();
  double ms = beta(n,0);

  for (i=0;i<nrpar-1;i++,workbeta++,workmode++)
    {
    *workbeta -= ms;
    *workmode -= smode;
    }


  FULLCOND_random::update();

  }


void FULLCOND_random_nongaussian::update_random_slope(void)
  {
  if (includefixed)
    update_random_slope_includefixed();
  else
    {
    nrtrials++;

    // Initializing variables

    if (optionsp->get_nriter()==1)
      {
      mode = beta;
      modeold = beta;
      weightiwls=datamatrix(likep->get_nrobs(),1);
      tildey=datamatrix(likep->get_nrobs(),1);
      proposal = datamatrix(beta.rows(),1,0);
      }

    if (lambdaconst == false)
      lambda=1.0/sigma2;
    else
      sigma2 = 1.0/lambda;

    unsigned i,j;
    double sumw;
    double sumy;
    double var;

    double u;

    double logold=0;
    double lognew=0;

    double qbetanew=0;
    double qbetaold=0;

    double diff;

    // Compute log-likelihood for old beta

    logold = likep->loglikelihood();

    // Use mode in the predictor

    update_linpred_diff(mode,beta);

    // Compute weightiwls and tildey

    likep->compute_IWLS_weight_tildey(weightiwls,tildey,column);

    int * workindex2 = index2.getV();

    double * workmode = mode.getV();
    double * workbeta = beta.getV();
    double * workproposal = proposal.getV();

    vector<unsigned>::iterator itbeg = posbeg.begin();
    vector<unsigned>::iterator itend = posend.begin();

    double * workweightiwls = weightiwls.getV()+(*workindex2);
    double * worktildey = tildey.getV()+(*workindex2);

    double * workdata2 = data2.getV();
    double * workdata = data.getV();

    modeold.assign(mode);

    for (i=0;i<nrpar;i++,workbeta++,workmode++,workproposal++,++itbeg,++itend)
      {

      sumw=0;
      sumy=0;
      for (j=*itbeg;j<=*itend;j++,workindex2++,
           workweightiwls+=*workindex2,worktildey+=*workindex2,
           workdata2++,workdata++)
        {
        sumw += (*workweightiwls) * (*workdata2);
        sumy += (*workdata)*(*workweightiwls)*
                 ((*worktildey)+(*workdata)*(*workmode));
        }

      var = 1.0/(sumw + lambda);

      *workmode = var*sumy;

      *workproposal = *workmode + sqrt(var)*rand_normal();

      logold -= 0.5*(*workbeta)*(*workbeta)*lambda;
      lognew -= 0.5*(*workproposal)*(*workproposal)*lambda;


      diff = *workbeta-*workmode;
      qbetaold -= 0.5*diff*diff/var;

      diff = *workproposal-*workmode;
      qbetanew -= 0.5*diff*diff/var;

      } // end:   for (i=0;i<nrpar;i++,workbeta++)

    update_linpred_diff(proposal,modeold);

    lognew += likep->loglikelihood();

    u = log(uniform());

    if (u <= (lognew - logold  + qbetaold - qbetanew) )
      {
      acceptance++;
      beta.assign(proposal);
      }
    else
      {
      update_linpred_diff(beta,proposal);
      }


    FULLCOND_random::update();


    }

  }


void FULLCOND_random_nongaussian::update(void)
  {

  if (!randomslope)
    {
//    if (nrpar >20)
      update_random_intercept_singleblock();
//    else
//      update_random_intercept();
    }
  else
    update_random_slope();
/*
  {

  unsigned i,j;
  double sumw;
  double std;
  double proposal;
  double u;

  double logold;
  double lognew;

  double qnew;
  double qold;

  double diff;

  double * workbeta = beta.getV();

  double ms=0;
  unsigned n;
  if (randomslope && includefixed)
    {
    n = nrpar-1;
    ms = beta(nrpar-1,0);
    }
  else
    n = nrpar;

  double h;
  datamatrix xwx(1,1);

  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  int * workindex;

  for (i=0;i<n;i++,workbeta++,++itbeg,++itend)
    {

    nrtrials++;

    if (randomslope && includefixed)
      h = *workbeta+ms;
    else
      h = *workbeta;

    if (!randomslope)
      {
      if (*itbeg==0)
        likep->set_linpredp_current(column);
      else
        likep->set_linpredp_current(index(*itbeg-1,0),column);
      sumw = likep->compute_sumweight2(*itbeg,*itend,index,index2,column);
      }
    else
      {
      sumw = likep->fisher2(*itbeg,*itend,index,data,column,true);
      }

    std = sqrt(1.0/(sumw + 1.0/sigma2));

    proposal = h + std*rand_normal();

    if (*itbeg==0)
      {
      likep->set_linpredp_current(column);
      }
    else
      {
      likep->set_linpredp_current(index(*itbeg-1,0),column);
      }

    logold = likep->loglikelihood2(*itbeg,*itend,index,index2);

    if (!randomslope)
      {
      if (*itbeg==0)
        likep->set_linpredp_current(column);
      else
        likep->set_linpredp_current(index(*itbeg-1,0),column);
      likep->add_linearpred2(proposal-h,*itbeg,*itend,index2,column);
      }
//      likep->add_linearpred(proposal - h,*itbeg,*itend,index,column);
    else
      {
      double * datap = data.getV()+(*itbeg);
      workindex = index.getV()+(*itbeg);
      for(j=*itbeg;j<=*itend;j++,datap++,workindex++)
        likep->add_linearpred((proposal - h)*(*datap),
          unsigned(*workindex),column);
      }


    if (*itbeg==0)
      {
      likep->set_linpredp_current(column);
      }
    else
      {
      likep->set_linpredp_current(index(*itbeg-1,0),column);
      }
    lognew = likep->loglikelihood2(*itbeg,*itend,index,index2);

    logold -= 0.5*(*workbeta-ms)*(*workbeta-ms)/sigma2;
    lognew -= 0.5*(proposal-ms)*(proposal-ms)/sigma2;


    diff = h-proposal;
    diff = diff*diff;
    qold = -(1.0/(2*std*std))*diff-log(std);

    if (!randomslope)
      {
      if (*itbeg==0)
        likep->set_linpredp_current(column);
      else
        likep->set_linpredp_current(index(*itbeg-1,0),column);
      sumw = likep->compute_sumweight2(*itbeg,*itend,index,index2,column);
      }
    else
      {
      sumw = likep->fisher2(*itbeg,*itend,index,data,column,true);
      }

    std = sqrt(1.0/(sumw + 1.0/sigma2));

    qnew = -(1.0/(2*std*std))*diff-log(std);

    u = log(uniform());

    if (u <= (lognew + qnew - logold - qold) )
      {
      acceptance++;
      *workbeta = proposal;
      }
    else
      {
      *workbeta = h;
      if (!randomslope)
        {
        if (*itbeg==0)
          likep->set_linpredp_current(column);
        else
          likep->set_linpredp_current(index(*itbeg-1,0),column);
        likep->add_linearpred2(h-proposal,*itbeg,*itend,index2,column);
        }
//        likep->add_linearpred(h-proposal,*itbeg,*itend,index,column);
      else
        {
        double * datap = data.getV()+(*itbeg);
        workindex = index.getV()+(*itbeg);

        for(j=*itbeg;j<=*itend;j++,workindex++,datap++)
          likep->add_linearpred((h-proposal)*(*datap),unsigned(*workindex),
                                column);
        }

      }

    } // end:   for (i=0;i<nrpar;i++,workbeta++)

  if (randomslope && includefixed)
    {

    workbeta = beta.getV();
    double s=0;
    for (i=0;i<nrpar-1;i++,workbeta++)
      s += *workbeta;
    s /= double(nrpar-1);

    double v = sigma2/double(nrpar-1);

    beta(nrpar-1,0) = s+sqrt(v)*rand_normal();

    workbeta = beta.getV();
    ms = beta(nrpar-1,0);
    for (i=0;i<nrpar-1;i++,workbeta++)
      *workbeta -= ms;

    }

  if (center)
    {
    double m = centerbeta();
    fcconst->update_intercept(m);
    }

  FULLCOND_random::update();

  update_spatialtotal();

  }

*/
  }


void FULLCOND_random_nongaussian::outresults(void)
  {
  FULLCOND_random::outresults();
  }


void FULLCOND_random_nongaussian::outoptions(void)
  {
  FULLCOND_random::outoptions();

  }

} // end: namespace MCMC



