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



#include"tvariance.h"

namespace MCMC
{


FULLCOND_tvariance::FULLCOND_tvariance(MCMCoptions * o,FULLCOND_nonp_basis * p,
                    unsigned & v,
                    const ST::string & ti, const ST::string & fp,
                    const ST::string & pres)
                     : FULLCOND(o,datamatrix(1,1),ti,p->get_nrpar(),1,fp)
  {
  Kp = p;
  pathresult = pres;
  setbeta(nrpar,1,1);
  u = datamatrix(nrpar,1);
  nu = v;
  if (p->get_type()==RW1)
    start=1;
  else if (p->get_type()==RW2)
    start = 2;
  plotstyle=plotnonp;
  }


void FULLCOND_tvariance::update(void)
  {


  register unsigned i;


  Kp->compute_u(u);

  double * workbeta=beta.getV()+start;
  double * worku = u.getV()+start;
  double v = double(nu)/2.0;
  for(i=start;i<nrpar;i++,workbeta++,worku++)
    {
    *workbeta = rand_invgamma(v+0.5,v+0.5* *worku);
    }


  acceptance++;
  Kp->updateK(beta);
  FULLCOND::update();
  }



void FULLCOND_tvariance::outresults(void)
  {

  FULLCOND::outresults();

  optionsp->out("  Results are stored in file " + pathresult + "\n");
  #if defined(JAVA_OUTPUT_WINDOW)
  optionsp->out("  Results may be visualized using method 'plotnonp'\n");
  optionsp->out("  Type for example: objectname.plotnonp " +
                   ST::inttostring(fcnumber) + " , median levels=none\n");
  #else
  optionsp->out("  Results may be visualized using the R function");
  optionsp->out(" 'plotnonp'\n");
  #endif
  optionsp->out("\n");

  unsigned i;

  ofstream outres(pathresult.strtochar());

  ST::string name = title;

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');


  outres << "intnr" << "   ";
  outres << "paramnr" << "   ";
  outres << "pmean   ";
  outres << "qu" << l1 << "   ";
  outres << "qu" << l2 << "   ";
  outres << "pmed   ";
  outres << "qu" << u1 << "   ";
  outres << "qu" << u2 << "   ";
  outres << endl;

  double * workmean = betamean.getV()+start;
  double * workbetaqu_l1 = betaqu_l1_lower.getV()+start;
  double * workbetaqu_l2 = betaqu_l2_lower.getV()+start;
  double * workbetaqu50 = betaqu50.getV()+start;
  double * workbetaqu_u1 = betaqu_l2_upper.getV()+start;
  double * workbetaqu_u2 = betaqu_l1_upper.getV()+start;

  for(i=start;i<nrpar;i++,workmean++,workbetaqu_l1++,workbetaqu_l2++,
                         workbetaqu50++,workbetaqu_u1++,workbetaqu_u2++)
    {
    outres << (i+1) << "   ";
    outres << (i+1) << "   ";
    outres << *workmean << "   ";
    outres << *workbetaqu_l1 << "   ";
    outres << *workbetaqu_l2 << "   ";
    outres << *workbetaqu50 << "   ";
    outres << *workbetaqu_u1 << "   ";
    outres << *workbetaqu_u2 << "   ";
    outres << endl;
    }

  }


void FULLCOND_tvariance::outoptions(void)
  {
  optionsp->out("  OPTIONS FOR NONPARAMETRIC TERM: " + title +
                " (variance parameters)\n",true);
  optionsp->out("\n");

  optionsp->out("  Hyperprior nu for variance parameter: " +
                ST::inttostring(nu) + "\n" );
  optionsp->out("\n");

  }




} // end: namespace MCMC





