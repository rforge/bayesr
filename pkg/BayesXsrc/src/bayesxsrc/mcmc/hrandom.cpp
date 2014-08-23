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



#include "hrandom.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------- class FULLCOND_hrandom -------------------------------------
//------------------------------------------------------------------------------

void FULLCOND_hrandom::init_name(const ST::string & na)
    {
    char charh = '_';
    ST::string stringh = "\\_";

    FULLCOND::init_name(na);

    ST::string helpname = na.insert_string_char(charh,stringh);
    term_symbolic = "f_{" +  helpname + "}("+helpname+")";

    if (column > 0)
      priorassumptions.push_back("$" + term_symbolic + "$" +
       " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");

    priorassumptions.push_back("structured Gaussian random effect");
    }


void FULLCOND_hrandom::init_names(const vector<ST::string> & na)
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
       " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");

    priorassumptions.push_back("structured Gaussian random effects");
    }




double FULLCOND_hrandom::centerbeta(void)
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


void FULLCOND_hrandom::compute_XWX(const datamatrix & weightmat,
                                  const unsigned & col)
  {

  register unsigned j,i;

  double * workXX = XX.getV();
  int *  workindex = index.getV();
  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  unsigned n = posbeg.size();

  for(j=0;j<n;j++,workXX++,++itbeg,++itend)
    {
    *workXX = 0;
    for (i=*itbeg;i<=*itend;i++,workindex++)
      *workXX += weightmat(*workindex,col);
    }

  }


double FULLCOND_hrandom::compute_quadform(void)
  {

  unsigned n;
  double sum = 0;
  double * workbeta = beta.getV();
  register unsigned i;

  n = nrpar;

  for(i=0;i<n;i++,workbeta++)
    sum += *workbeta * *workbeta;

  return sum;

  }


unsigned FULLCOND_hrandom::get_rankK(void)
  {
  return nrpar;
  }

unsigned FULLCOND_hrandom::get_rankK2(void)
  {
  return get_rankK();
  }


void FULLCOND_hrandom::set_lambdaconst(double la)
  {
  lambda=la;
  lambdaconst = true;
  }


FULLCOND_hrandom::FULLCOND_hrandom (MCMCoptions * o,DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & d, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,
                              const double & la, const unsigned & c)
                            : FULLCOND(o,datamatrix(1,1),t,1,1,fp)
  {

  fctype = randomeffects;

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

  muy = datamatrix(nrpar,1);

  mu = datamatrix(index.rows(),1);

  identifiable =true;

  }


FULLCOND_hrandom::FULLCOND_hrandom(const FULLCOND_hrandom & fc)
                            : FULLCOND(FULLCOND(fc))
  {
  mu = fc.mu;
  muy = fc.muy;
  XX = fc.XX;
  likep = fc.likep;
  likep_RE = fc.likep_RE;
  index = fc.index;
  index2 = fc.index2;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effvalues = fc.effvalues;
  sigma2 = fc.sigma2;
  lambda = fc.lambda;
  lambdaold1 = fc.lambdaold1;
  lambdaold2 = fc.lambdaold2;
  df_lambdaold1 = fc.df_lambdaold1;
  df_lambdaold2 = fc.df_lambdaold2;
  lambdaconst=fc.lambdaconst;
  data2 = fc.data2;
  }


const FULLCOND_hrandom & FULLCOND_hrandom::
         operator=(const FULLCOND_hrandom & fc)
  {
  if (this==&fc)
    return *this;

  FULLCOND::operator=(FULLCOND(fc));

  mu = fc.mu;  
  muy = fc.muy;
  XX = fc.XX;
  likep = fc.likep;
  likep_RE = fc.likep_RE;
  index = fc.index;
  index2 = fc.index2;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effvalues = fc.effvalues;
  sigma2 = fc.sigma2;
  lambda = fc.lambda;
  lambdaold1 = fc.lambdaold1;
  lambdaold2 = fc.lambdaold2;
  df_lambdaold1 = fc.df_lambdaold1;
  df_lambdaold2 = fc.df_lambdaold2;
  lambdaconst=fc.lambdaconst;
  data2 = fc.data2;

  return *this;
  }


void FULLCOND_hrandom::update(void)
  {

  double var;
  double m;
  unsigned i,j;
  unsigned n = nrpar;

  if (optionsp->get_nriter()==1 || changingweight)
    compute_XWX(likep->get_weight(),0);

  double scale = likep->get_scale(column);

  if (lambdaconst == false)
    lambda = scale/sigma2;
  else
    sigma2 = scale/lambda;


//  double sqrtscale = sqrt(scale);


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


  // speichert in muy X'W (y-tilde(eta))

  for(i=0;i<nrpar;i++,workmuy++,++itbeg,++itend)
    {
    *workmuy = 0;
    for(j=*itbeg;j<=*itend;j++,workindex2++)
      {
      mup += *workindex2;
      *workmuy+= likep->get_weight(*workindex2)* *mup;
      }

    }


  workbeta = beta.getV();
  workmuy = muy.getV();
  double * workXX = XX.getV();
  for (i=0;i<n;i++,workbeta++,workmuy++,workXX++)
    {

    var = 1.0/(*workXX/scale  + 1/sigma2);

    m = var * (*workmuy/scale+ likep_RE->get_linearpred(i,column)/sigma2);

    *workbeta = m + sqrt(var)*rand_normal();

    }

  update_linpred(true);


  acceptance++;


  transform = likep->get_trmult(column);

  likep_RE->set_response(beta);

  FULLCOND::update();

  }


void FULLCOND_hrandom::outresults(void)
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

  optionsp->out("  Results for random effects are stored in file\n");
  optionsp->out("  " + pathcurrent + "\n");

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


void FULLCOND_hrandom::outoptions(void)
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


void FULLCOND_hrandom::update_linpred(const bool & add)
  {
  unsigned i;

  vector<unsigned>::iterator itbeg = posbeg.begin();
  vector<unsigned>::iterator itend = posend.begin();
  double * workbeta = beta.getV();

  if (add==false)
    {
//  likep->set_linpredp_current(column);
    for (i=0;i<nrpar;i++,++itbeg,++itend,workbeta++)
      if (*itbeg != -1)
        likep->add_linearpred2(-(*workbeta),*itbeg,*itend,index,index2,column);
    } // end: if (add == false)
  else
    {

//      likep->set_linpredp_current(column);
    for (i=0;i<nrpar;i++,++itbeg,++itend,workbeta++)
      if (*itbeg != -1)
        likep->add_linearpred2(*workbeta,*itbeg,*itend,index,index2,column);

    }

  }


void FULLCOND_hrandom::get_effectmatrix(datamatrix & e,
                                vector<ST::string> & enames,
                                unsigned be,unsigned en, effecttype t)
  {


  }



bool FULLCOND_hrandom::posteriormode(void)
  {

  sigma2=0.1;
  unsigned n = nrpar;

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


  for(i=0;i<nrpar;i++,workmuy++,++itbeg,++itend)
    {
    *workmuy = 0;
    for(j=*itbeg;j<=*itend;j++,workindex2++)
      {
      *workmuy+= likep->get_workingres(*workindex2);
      }
    }


  itbeg = posbeg.begin();
  itend = posend.begin();
  workmuy = muy.getV();
  double * workbeta = beta.getV();
  double * workXX = XX.getV();

  double var;

  double scale = likep->get_scale(column);

  for(i=0;i<n;i++,workmuy++,++itbeg,++itend,workbeta++,workXX++)
    {
    var = 1.0/(*workXX/scale  + 1/sigma2);
    *workbeta =  var*( (*workmuy)/scale  + likep_RE->get_linearpred(i,column)/sigma2);
    }


  update_linpred(true);

  transform = likep->get_trmult(column);

  likep_RE->set_response(beta);

  return FULLCOND::posteriormode();

  }





} // end: namespace MCMC




