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



#include<mcmc_nonpbasis.h>

namespace MCMC
{

double FULLCOND_nonp_basis::compute_quadform(void)
  {
  return K.compute_quadform(beta,0);
  }

double FULLCOND_nonp_basis::compute_sumfabsdiff(void)
  {
  return 0.0;
  }

void FULLCOND_nonp_basis::updateK(const datamatrix & q)
  {


  updateKenv(q);

/*
  register unsigned i;

  if (type==RW1)
    {
    double * workdiag=K.getdiagpointer();
    double * workupper = K.getupperpointer();
    double * workg = g.getV()+1;
    double * workq = q.getV()+1;
    double wold=1.0/(*workq * *workg);
    double wnew;

    *workdiag = wold;
    *workupper = -wold;

    workdiag++;
    workupper++;
    workq++;
    workg++;
    for (i=1;i<K.getdim()-1;i++,workdiag++,workupper++,workq++,workg++)
      {
      wnew = 1.0/(*workq * *workg);
      *workdiag = wold+wnew;
      *workupper = -wnew;
      wold = wnew;
      }

    *workdiag = wold;


//  NEU

//    ofstream out("c:\\temp\\K.raw");
//      K.print2(out);

// ENDE NEU


    }
  else if (type==RW2)
    {

    double * workdiag=K.getdiagpointer();
    double * workupper = K.getupperpointer();
//    double * workq = q.getV()+2;


    *workdiag = (F2(2,0)*F2(2,0))/(g(2,0)*q(2,0));
    workdiag++;
    *workdiag = (F1(2,0)*F1(2,0))/(g(2,0)*q(2,0)) +
                (F2(3,0)*F2(3,0))/(g(3,0)*q(3,0));
    workdiag++;

    *workupper = (F1(2,0)*F2(2,0))/(g(2,0)*q(2,0));                               // (1,2)
    workupper++;
    *workupper = F2(2,0)/(g(2,0)*q(2,0));                                         // (1,3)
    workupper++;

    *workupper = F1(2,0)/(g(2,0)*q(2,0)) + (F1(3,0)*F2(3,0))/(g(3,0)*q(3,0));     // (2,3)
    workupper++;

    *workupper = F2(3,0)/(g(3,0)*q(3,0));   // (2,4)
    workupper++;

    for(i=2;i<nrpar-2;i++,workdiag++,workupper++)
      {
      *workdiag = (F1(i+1,0)*F1(i+1,0))/(g(i+1,0)*q(i+1,0)) +
                (F2(i+2,0)*F2(i+2,0))/(g(i+2,0)*q(i+2,0)) + 1/(g(i,0)*q(i,0));

      *workupper = F1(i+1,0)/(g(i+1,0)*q(i+1,0)) +
                   (F1(i+2,0)*F2(i+2,0))/(g(i+2,0)*q(i+2,0));  // (i,i+1)
      workupper++;

      *workupper = F2(i+2,0)/(g(i+2,0)*q(i+2,0));          // (i,i+2)

      }


    *workdiag = (F1(nrpar-1,0)*F1(nrpar-1,0))/(g(nrpar-1,0)*q(nrpar-1,0)) +
                1/(g(nrpar-2,0)*q(nrpar-2,0));
    workdiag++;
    *workdiag = 1/(g(nrpar-1,0)*q(nrpar-1,0));

    *workupper = F1(nrpar-1,0)/(g(nrpar-1,0)*q(nrpar-1,0));

    }

  ofstream out1("c:\\tmp\\Kenv");
  Kenv.print2(out1);
  out1.close();


  ofstream out2("c:\\tmp\\K");
  K.print2(out2);
  out2.close();
*/

  }

void FULLCOND_nonp_basis::updateKenv(const datamatrix & q)
  {

//ofstream out1("c:\\cprog\\test\\results\\Kenv_alt.txt");
//Kenv.print2(out1);
//out1.close();

  register unsigned i;

  if (type==RW1)
    {
    vector<double>::iterator workdiag = Kenv.getDiagIterator();
    vector<double>::iterator workenv = Kenv.getEnvIterator();
    double * workg = g.getV()+1;
    double * workq = q.getV()+1;
    double wold=1.0/(*workq * *workg);
    double wnew;

    *workdiag = wold;
    *workenv = -wold;

    workdiag++;
    workenv++;
    workq++;
    workg++;
    for (i=1;i<nrpar-1;i++,workdiag++,workenv++,workq++,workg++)
      {
      wnew = 1.0/(*workq * *workg);
      *workdiag = wold+wnew;
      *workenv = -wnew;
      wold = wnew;
      }

    *workdiag = wold;
    }
  else if (type==RW2)
    {

    vector<double>::iterator workdiag=Kenv.getDiagIterator();
    vector<double>::iterator workenv = Kenv.getEnvIterator();
//    double * workq = q.getV()+2;


    *workdiag = (F2(2,0)*F2(2,0))/(g(2,0)*q(2,0));
    workdiag++;
    *workdiag = (F1(2,0)*F1(2,0))/(g(2,0)*q(2,0)) +
                (F2(3,0)*F2(3,0))/(g(3,0)*q(3,0));
    workdiag++;

    *workenv = (F1(2,0)*F2(2,0))/(g(2,0)*q(2,0));                 //(2,1)
    workenv++;
    *workenv = F2(2,0)/(g(2,0)*q(2,0));                           //(3,1)
    workenv++;

    *workenv = F1(2,0)/(g(2,0)*q(2,0)) + (F1(3,0)*F2(3,0))/(g(3,0)*q(3,0));     //(3,2)
    workenv++;

    *workenv = F2(3,0)/(g(3,0)*q(3,0));                           //(4,2)
    workenv++;

    for(i=2;i<nrpar-2;i++,workdiag++,workenv++)
      {
      *workdiag = (F1(i+1,0)*F1(i+1,0))/(g(i+1,0)*q(i+1,0)) +
                (F2(i+2,0)*F2(i+2,0))/(g(i+2,0)*q(i+2,0)) + 1/(g(i,0)*q(i,0));

      *workenv = F1(i+1,0)/(g(i+1,0)*q(i+1,0)) +
                   (F1(i+2,0)*F2(i+2,0))/(g(i+2,0)*q(i+2,0));      //(i+1,i)
      workenv++;

      *workenv = F2(i+2,0)/(g(i+2,0)*q(i+2,0));                    //(i+2,i)

      }


    *workdiag = (F1(nrpar-1,0)*F1(nrpar-1,0))/(g(nrpar-1,0)*q(nrpar-1,0)) +
                1/(g(nrpar-2,0)*q(nrpar-2,0));
    workdiag++;
    *workdiag = 1/(g(nrpar-1,0)*q(nrpar-1,0));

    *workenv = F1(nrpar-1,0)/(g(nrpar-1,0)*q(nrpar-1,0));

    }

//ofstream out2("c:\\cprog\\test\\results\\Kenv_neu.txt");
//Kenv.print2(out2);
//out2.close();

  }

void FULLCOND_nonp_basis::updateKenv_alpha(const double alpha1, const double alpha2)
  {
  unsigned i;

  if (type==RW1)
    {
    double alpha1_2 = alpha1*alpha1;

    vector<double>::iterator workdiag=Kenv.getDiagIterator();
    vector<double>::iterator workenv =Kenv.getEnvIterator();

    *workenv  = -alpha1;

    workdiag++;
    workenv++;

    for (i=1;i<nrpar-1;i++,workdiag++,workenv++)
      {
      *workdiag = 1.0 + alpha1_2;
      *workenv  = -alpha1;
      }
    }
  else if (type==RW2)
    {
    double alpha1_2 = alpha1*alpha1;
    double alpha2_2 = alpha2*alpha2;

    vector<double>::iterator workdiag=Kenv.getDiagIterator();
    vector<double>::iterator workenv = Kenv.getEnvIterator();

    workdiag++;

    *workdiag = 1.0 + alpha1_2;
    workdiag++;

    *workenv = alpha1;                 //(2,1)
    workenv++;
    *workenv = alpha2;                           //(3,1)
    workenv++;

    *workenv = alpha1*(1.0+alpha2);     //(3,2)
    workenv++;

    *workenv = alpha2;                           //(4,2)
    workenv++;

    for(i=2;i<nrpar-2;i++,workdiag++,workenv++)
      {
      *workdiag = 1.0 + alpha1_2 + alpha2_2;

      *workenv = alpha1*(1.0+alpha2);      //(i+1,i)
      workenv++;

      *workenv = alpha2;                    //(i+2,i)
      }

    *workdiag = 1.0 + alpha1_2;
    *workenv = alpha1;

    }
  else if (type==mrf)
    {
    vector<double>::iterator workenv =Kenv.getEnvIterator();
    for (i=0;i<Kenv.getXenv(Kenv.getDim());i++,workenv++)
      {
      if(*workenv != 0.0)
        *workenv  = -alpha1;
      }
    }

  }


void FULLCOND_nonp_basis::set_stationary(double alphastart)
  {
  if(type==RW2)
    updateKenv_alpha(-(alphastart+alphastart),alphastart*alphastart);
  else
    updateKenv_alpha(alphastart);
  rankK = nrpar;
  }


double FULLCOND_nonp_basis::getLogDet(void)
  {
  return Kenv.getLogDet();
  }

void FULLCOND_nonp_basis::compute_u(datamatrix & u)
  {
  register unsigned i;
  if (type== RW1)
    {
    double * worku = u.getV()+1;
    double * workg = g.getV()+1;
    double * workbeta1 = beta.getV();
    double * workbeta2 = beta.getV()+1;
    double v;
    for(i=1;i<nrpar;i++,worku++,workbeta1++,workbeta2++,workg++)
      {
      v = *workbeta2-*workbeta1;
      *worku = (v*v)/(sigma2* * workg);
      }
    }
  else if (type==RW2)
    {

    double * worku = u.getV()+2;
    double * workbeta1 = beta.getV();
    double * workbeta2 = beta.getV()+1;
    double * workbeta3 = beta.getV()+2;
    double * workg = g.getV()+2;
    double * workF1 = F1.getV()+2;
    double * workF2 = F2.getV()+2;

    double v;
    for(i=2;i<nrpar;i++,worku++,workbeta1++,workbeta2++,workbeta3++,workg++,
                    workF1++,workF2++)
      {
      v = *workbeta3 + *workF1 * *workbeta2 + *workF2 * *workbeta1;
      *worku = (v*v)/(sigma2* *workg);
      }

    }

  }


double FULLCOND_nonp_basis::compute_ui(unsigned i)
  {
  double v=0.0;

  if (type ==RW1)
    v = beta(i,0)-beta(i-1,0);
  else if (type==RW2)
    v = beta(i,0) + F1(i,0)*beta(i-1,0)+F2(i,0)*beta(i-2,0);

  return v*v;
  }


void FULLCOND_nonp_basis::set_adaptiv(void)
  {
  adaptiv = true;
  if (type==RW1)
    {
    g = datamatrix(nrpar,1,0);
    unsigned s;
    for (s=1;s<nrpar;s++)
      {
      g(s,0) = 1.0;
      }
    }
  else if (type==RW2)
    {

    F1 = datamatrix(nrpar,1,0);
    F2 = datamatrix(nrpar,1,0);
    g  = datamatrix(nrpar,1,0);
    unsigned s;
    for (s=2;s<nrpar;s++)
      {
      F1(s,0) = -(1+weight[s]/weight[s-1]);
      F2(s,0) = weight[s]/weight[s-1];
      g(s,0) = 1.0;
      }

    } // end: type = RW2


  }



FULLCOND_nonp_basis::FULLCOND_nonp_basis(MCMCoptions * o,DISTRIBUTION * dp,
                            const fieldtype & ft,const ST::string & ti,
                            const ST::string & fp,const ST::string & pres,
                            const unsigned & c,const unsigned & per)
  : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
  {
  interaction=false;
  changingweight = dp->get_changingweight();
  tildey = datamatrix(1,1);
  column = c;
  type = ft;
  period = per;
  pathresults = pres;
  pathresult = pres;
  pathcurrent = pres;
  sigma2 = 10;
  likep = dp;
  adaptiv = false;
  f = 2.0;
  oldacceptance = 0;
  oldnrtrials = 0;
  contourprob = -1;
  fc_contour = FULLCOND();
  }




  // COPY CONSTRUCTOR

FULLCOND_nonp_basis::FULLCOND_nonp_basis(const FULLCOND_nonp_basis & fc)
    : FULLCOND(FULLCOND(fc))
  {
  invprec = fc.invprec;
  adaptiv = fc.adaptiv;
  F1 = fc.F1;
  F2 = fc.F2;
  g = fc.g;
  type = fc.type;
  period = fc.period;
  distance = fc.distance;
  index = fc.index;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effectvalues = fc.effectvalues;
  effectvdouble = fc.effectvdouble;
  K = fc.K;
  Kenv = fc.Kenv;
  Ksp = fc.Ksp;
  rankK = fc.rankK;
  likep = fc.likep;
  sigma2 = fc.sigma2;
  pathresults = fc.pathresults;
  varcoeff = fc.varcoeff;
  changingweight = fc.changingweight;
  tildey = fc.tildey;
  interaction=fc.interaction;
  lambda = fc.lambda;
  polex = fc.polex;
  f = fc.f;
  oldacceptance = fc.oldacceptance;
  oldnrtrials = fc.oldnrtrials;
  fc_contour = fc.fc_contour;
  contourprob = fc.contourprob;
  }


const FULLCOND_nonp_basis & FULLCOND_nonp_basis::operator=(
const FULLCOND_nonp_basis & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND::operator=(FULLCOND(fc));
  invprec = fc.invprec;
  adaptiv = fc.adaptiv;
  F1 = fc.F1;
  F2 = fc.F2;
  g = fc.g;
  type = fc.type;
  period = fc.period;
  distance = fc.distance;
  index = fc.index;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effectvalues = fc.effectvalues;
  effectvdouble = fc.effectvdouble;
  K = fc.K;
  Kenv = fc.Kenv;
  Ksp = fc.Ksp;
  rankK = fc.rankK;
  likep = fc.likep;
  sigma2 = fc.sigma2;
  pathresults = fc.pathresults;
  varcoeff = fc.varcoeff;
  changingweight = fc.changingweight;
  tildey = fc.tildey;
  interaction=fc.interaction;
  lambda = fc.lambda;
  polex = fc.polex;
  f = fc.f;
  oldacceptance = fc.oldacceptance;
  oldnrtrials = fc.oldnrtrials;
  fc_contour = fc.fc_contour;
  contourprob = fc.contourprob;
  return *this;
  }


void FULLCOND_nonp_basis::update(void)
  {
  FULLCOND::update();
  }

bool FULLCOND_nonp_basis::posteriormode(void)
  {
  return FULLCOND::posteriormode();
  }

void FULLCOND_nonp_basis::outresults(void)
  {

  FULLCOND::outresults();

  optionsp->out("  Results are stored in file\n");
  optionsp->out("  " +   pathcurrent + "\n");
  optionsp->out("\n");

  if (type == MCMC::mrf)
    {
    #if defined(JAVA_OUTPUT_WINDOW)

    if (polex == true)
      {
      optionsp->out("  Postscript file is stored in file\n");
      ST::string psfile = pathcurrent.substr(0,pathcurrent.length()-4) + ".ps";
      optionsp->out("  " + psfile + "\n");
      optionsp->out("\n");
      }

    optionsp->out("  Results may be visualized in BayesX using method 'drawmap'\n");
    optionsp->out("  Type for example: objectname.drawmap " +
    ST::inttostring(fcnumber) + "\n");
    optionsp->out("\n");
    #elif defined(BORLAND_OUTPUT_WINDOW)
    optionsp->out("  Results may be visualized using the R function");
    optionsp->out(" 'drawmap'\n");
    optionsp->out("\n");
    #endif
    }
  else
    {
    #if defined(JAVA_OUTPUT_WINDOW)
    optionsp->out("  Postscript file is stored in file\n");
    ST::string psfile = pathcurrent.substr(0,pathcurrent.length()-4) + ".ps";
    optionsp->out("  " + psfile + "\n");
    optionsp->out("\n");
    optionsp->out("  Results may be visualized in BayesX using method 'plotnonp'\n");
    optionsp->out("  Type for example: objectname.plotnonp " +
    ST::inttostring(fcnumber) + "\n");
    optionsp->out("\n");
    #else
    char hchar = '\\';
    ST::string hstring = "/";
    ST::string pathresultsplus = pathcurrent.insert_string_char(hchar,hstring);
    ST::string psfile = pathresultsplus.substr(0,pathresultsplus.length()-4)
    + ".ps";
    optionsp->out("  Results may be visualized using the R function 'plotnonp'");
    optionsp->out("\n");
    optionsp->out("  Type for example:\n");
    optionsp->out("\n");
    optionsp->out("  plotnonp(\""+ pathresultsplus + "\")\n");
    optionsp->out("\n");
    #endif
    }
  optionsp->out("\n");


  if (optionsp->get_samplesize() == 0)
    {
    double df = compute_df();
    optionsp->out("  Approximate degrees of freedom: "
                  + ST::doubletostring(df,6) + "\n");
    optionsp->out("\n");
    }

  unsigned i;

  ofstream outres(pathcurrent.strtochar());
//  ST::string name = title;
  ST::string name = datanames[0];

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  outres << "intnr" << "   ";
  outres << name << "   ";
  outres << "pmean   ";
  outres << "pqu"  << l1  << "   ";
  outres << "pqu"  << l2  << "   ";
  outres << "pmed   ";
  outres << "pqu"  << u1  << "   ";
  outres << "pqu"  << u2  << "   ";
  outres << "pcat" << level1 << "   ";
  outres << "pcat" << level2 << "   ";

  outres << endl;

    double * workmean = betamean.getV();
    double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
    double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
    double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
    double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
    double * workbetaqu50 = betaqu50.getV();

    vector<ST::string>::iterator effit = effectvalues.begin();

    for(i=0;i<nrpar;i++,++effit,workmean++,workbetaqu_l1_lower_p++,
                            workbetaqu_l2_lower_p++,workbetaqu50++,
                            workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
      {
      outres << (i+1) << "   ";
      outres << *effit << "   ";
      outres << *workmean << "   ";
      outres << *workbetaqu_l1_lower_p << "   ";
      outres << *workbetaqu_l2_lower_p << "   ";
      outres << *workbetaqu50 << "   ";
      outres << *workbetaqu_l2_upper_p << "   ";
      outres << *workbetaqu_l1_upper_p << "   ";

      if (*workbetaqu_l1_lower_p > 0)
        outres << 1 << "   ";
      else if (*workbetaqu_l1_upper_p < 0)
        outres << -1 << "   ";
      else
        outres << 0 << "   ";

      if (*workbetaqu_l2_lower_p > 0)
        outres << 1 << "   ";
      else if (*workbetaqu_l2_upper_p < 0)
        outres << -1 << "   ";
      else
        outres << 0 << "   ";

      outres << endl;
      }
  }


void FULLCOND_nonp_basis::outoptions(void)
  {
  optionsp->out("  OPTIONS FOR NONPARAMETRIC TERM: " + title + "\n",true);
  optionsp->out("\n");
  ST::string typestr;
  if (type == RW1)
    typestr = "first order random walk";
  else if (type == RW2)
    typestr = "second order random walk";
  else if (type == seasonal)
    typestr = "seasonal component";
  else if (type==mrf)
    typestr = "spatial Markov random field";
  else if (type==mrfkronecker)
    typestr= "Kronecker product interaction";
  else if (type==mrflinear)
    typestr = "2 dimensional first order random walk";
  else if (type==mrfkr1)
    typestr = "Kronecker product interaction (RW1*RW1)";
  else if (type==mrfkr2)
    typestr = "Kronecker product interaction (RW2*RW2)";

  optionsp->out("  Prior: " + typestr + "\n");
  if (type==seasonal)
    optionsp->out("  Period: " + ST::inttostring(period) + "\n");

  }


void FULLCOND_nonp_basis::write_contour(const datamatrix & m, const double & scaleinv,
                         const double & sigma2inv,
                         const double & bXXb, const double & bKb, const double & mPm, const double & logDetP,
                         const envmatdouble * prec_envp)
  {

  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
    {
    unsigned i,j;

    double * contourp = fc_contour.getbetapointer();
    for(i=0;i<nrpar;i++,contourp++)
      *contourp = m(i,0);

    *contourp = scaleinv;                           // contour(,nrpar)
    contourp++;
    *contourp = sigma2inv;                          // contour(,nrpar+1)

    contourp++;
    *contourp = bXXb;                               // contour(,nrpar+2)
    contourp++;
    *contourp = bKb;                                // contour(,nrpar+3)
    contourp++;
    *contourp = mPm;                                // contour(,nrpar+4)
    contourp++;
    *contourp = logDetP;                            // contour(,nrpar+5)

    datamatrix mPhelp = datamatrix(nrpar,1,0);
    for(i=0;i<nrpar;i++)
      for(j=0;j<nrpar;j++)
        mPhelp(i,0) += m(j,0)*(*prec_envp)(j,i);

    contourp++;
    for(i=0;i<nrpar;i++,contourp++)
      *contourp = mPhelp(i,0);

    }

  }


void FULLCOND_nonp_basis::reset(void)
  {
  FULLCOND::reset();
  sigma2 = 10;
  }


void FULLCOND_nonp_basis::predict(const datamatrix & newX, datamatrix & linpred)
  {
  }


void FULLCOND_nonp_basis::tune_updatetau(const rate & r)
  {

  double acceptrate;
  if (nrtrials == 0)
    acceptrate = (double(acceptance-oldacceptance)/double(100))*100;
  else
    acceptrate = (double(acceptance-oldacceptance)/double(nrtrials-oldnrtrials))*100;

  oldacceptance = acceptance;
  oldnrtrials = nrtrials;

#define DOITU(rate,scale) if (acceptrate < rate) { f = pow(f,1.0/scale);}
#define DOITD(rate,scale) if (acceptrate > rate) { f = pow(f,scale);}

  switch(r)
    {
    case 10: DOITU(2, 5.0);
             DOITU(4, 2.0);
        	 DOITU(6, 1.5);
        	 DOITU(8, 1.2);
        	 DOITU(10, 1.05);

             DOITD(80, 5.0);
             DOITD(60, 3.0);
        	 DOITD(40, 2.0);
         	 DOITD(30, 1.7);
             DOITD(20, 1.5);
        	 DOITD(15, 1.3);
          	 DOITD(14, 1.1);
         	 DOITD(13, 1.05);
             break;
    case 30: DOITU(10, 5.0);
         	 DOITU(15, 2.0);
         	 DOITU(20, 1.2);
        	 DOITU(25, 1.05);

        	 DOITD(80, 3.0);
        	 DOITD(70, 2.0);
        	 DOITD(60, 1.5);
        	 DOITD(50, 1.3);
        	 DOITD(45, 1.2);
        	 DOITD(40, 1.1);
        	 DOITD(35, 1.05);
        	 DOITD(30, 1.025);
    break;
    case 50: DOITD(90, 5.0);
        	 DOITD(80, 3.0);
        	 DOITD(70, 2.0);
        	 DOITD(60, 1.2);
        	 DOITD(55, 1.05);

        	 DOITU(10, 2.0);
        	 DOITU(20, 1.5);
        	 DOITU(30, 1.3);
        	 DOITU(40, 1.2);
        	 DOITU(45, 1.1);
        	 DOITU(50, 1.05);
    break;
    case 60: DOITD(90, 5.0);
        	 DOITD(80, 2.0);
        	 DOITD(70, 1.2);
        	 DOITD(65, 1.1);
        	 DOITD(60, 1.05);

        	 DOITU(20, 3.0);
        	 DOITU(30, 2.0);
        	 DOITU(40, 1.5);
        	 DOITU(50, 1.2);
        	 DOITU(55, 1.1);
        	 DOITU(60, 1.05);
    break;
    case 70: DOITD(90, 5.0);
        	 DOITD(85, 2.0);
        	 DOITD(80, 1.2);
        	 DOITD(75, 1.05);

        	 DOITU(20, 3.0);
        	 DOITU(30, 2.0);
        	 DOITU(40, 1.5);
        	 DOITU(50, 1.3);
        	 DOITU(55, 1.2);
        	 DOITU(60, 1.1);
        	 DOITU(65, 1.05);
        	 DOITU(70, 1.025);
    break;
    case 80: DOITD(95, 2.0);
        	 DOITD(90, 1.5);
        	 DOITD(85, 1.1);

        	 DOITU(30, 3.0);
        	 DOITU(40, 2.0);
        	 DOITU(50, 1.5);
        	 DOITU(60, 1.2);
        	 DOITU(70, 1.1);
        	 DOITU(75, 1.05);
        	 DOITU(80, 1.025);
    break;
    };

#undef DOITU
#undef DOITD

  if(f < 1.1)
    f = 1.1;

  }


} // end: namespace MCMC



