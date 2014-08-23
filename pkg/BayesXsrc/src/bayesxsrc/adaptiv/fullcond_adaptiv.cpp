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





#include "fullcond_adaptiv.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//---------- CLASS: FULLCOND_adaptiv (implementation of member fucntions) ------
//------------------------------------------------------------------------------

FULLCOND_adaptiv::FULLCOND_adaptiv(MCMCoptions * o, FULLCOND_nonp_basis * p,
                             const fieldtype & ft,
                             const ST::string & ti,const double & a,
                             const double & b,
                             const bool & uniformb, const double & startvalue,
                             const unsigned & minb,const unsigned & maxb,
                             const ST::string & fp, const ST::string & pres)
  : FULLCOND(o,datamatrix(1,1),ti,p->get_nrpar(),1,fp)
  {
  unifb = uniformb;

  Gp = p;
  type=ft;

  if (Gp->get_type() == RW1)
    start=1;
  else
    start=2;

  Gp->update_sigma2(1);
//  sigma2beta=1;

  identifiable = true;
  sigma2 = startvalue;
  sigma2sum = 0;
  a_invgamma = a;
  b_invgamma = b;

  pathresults = pres;

  setbeta(p->get_nrpar(),2,startvalue);

  h = datamatrix(p->get_nrpar()-start,1,log(startvalue));

  datamatrix md(p->get_nrpar()-start,1);
  unsigned i;
  for(i=0;i<md.rows();i++)
    md(i,0) = i+1;

  minblocksize = minb;
  maxblocksize = maxb;
  Pm = PenaltyMatrix(md,"varhelp",p->get_nrpar(),minb,maxb,ft);
  }



FULLCOND_adaptiv::FULLCOND_adaptiv(const FULLCOND_adaptiv & fc)
  : FULLCOND(FULLCOND(fc))
  {
  unifb = fc.unifb;
  sigma2 = fc.sigma2;
  sigma2sum = fc.sigma2sum;
  a_invgamma = fc.a_invgamma;
  b_invgamma = fc.b_invgamma;
  type = fc.type;
  pathresults = fc.pathresults;
  Gp = fc.Gp;
//  sigma2beta = fc.sigma2beta;
  start = fc.start;
  h = fc.h;
  Pm = fc.Pm;
  minblocksize = fc.minblocksize;
  maxblocksize = fc.maxblocksize;
  }


const FULLCOND_adaptiv & FULLCOND_adaptiv::operator=(const FULLCOND_adaptiv & fc)
  {
  if (this == &fc)
    return *this;
  unifb = fc.unifb;
  sigma2 = fc.sigma2;
  sigma2sum = fc.sigma2sum;
  a_invgamma = fc.a_invgamma;
  b_invgamma = fc.b_invgamma;
  type = fc.type;
  pathresults = fc.pathresults;
  Gp = fc.Gp;
//  sigma2beta = fc.sigma2beta;
  start = fc.start;
  h = fc.h;
  Pm = fc.Pm;
  minblocksize = fc.minblocksize;
  maxblocksize = fc.maxblocksize;
  return *this;
  }


void FULLCOND_adaptiv::outoptions(void)
  {
  ST::string typestr;
  if (type == RW1)
    typestr = "first order random walk";
  else if (type == RW2)
    typestr = "second order random walk";

  optionsp->out("\n");
  optionsp->out("  OPTIONS FOR NONPARAMETRIC TERM: " + title + "\n",true);
  optionsp->out("\n");
  optionsp->out("  Prior: " + typestr + "\n");
  optionsp->out("  Hyperprior a for variance parameter: " +
                ST::doubletostring(a_invgamma) + "\n" );
  if (unifb==true)
    optionsp->out("  Hyperprior b for variance parameter: uniform\n" );
  else
    optionsp->out("  Hyperprior b for variance parameter: " +
                     ST::doubletostring(b_invgamma) + "\n" );

  optionsp->out("  Minimum blocksize for blockmove updates: " +
                   ST::inttostring(minblocksize) + "\n");
  optionsp->out("  Maximum blocksize for blockmove updates: " +
                   ST::inttostring(maxblocksize) + "\n");
  optionsp->out("  Variance: " + ST::doubletostring(sigma2,6)
                   + "\n");
  }


/*
double FULLCOND_adaptiv::compute_hprop(unsigned & i)
  {
  double mu;
  double sigma = sqrt(sigma2/2.0);
  if (i==start)
    {
    mu = 0.5*h(1,0);
    return rand_normal()*sigma + mu;
    }
  else if (i==beta.rows()-1)
    {
    mu = h(i-1-start,0);
    return rand_normal()*sigma + mu;
    }
  else
    {
    mu = 0.5*h(i+1-start,0) + 0.5*h(i-1-start,0);
    return rand_normal()*sigma + mu;
    }
  }
*/


double FULLCOND_adaptiv::compute_denquot(unsigned i,double hp)
  {
  double hdiff,qt,qtp;
  hdiff = h(i-start,0) - hp;
  qt = beta(i,0);
  qtp = exp(hp);
  return 0.5*hdiff + 0.5*Gp->compute_ui(i)*(1/qt-1/qtp);
  }


void FULLCOND_adaptiv::update(void)
  {
  unsigned j,k;
  double denquot,u;

  unsigned blocksize = int(Pm.get_minsize() +
                       uniform()*(Pm.get_maxsize()-Pm.get_minsize()+1));
  unsigned an = 1;
  unsigned en = blocksize;

  for(j=0;j<Pm.get_nrblocks(blocksize);j++)
    {
    nrtrials++;

// Vorschlag:
# if defined(__BUILDING_GNU)
    Pm.compute_fc(h,blocksize,an,en,sqrt(sigma2));
# else
    Pm.compute_fc(h,blocksize,an,en,sqrtl(sigma2));
#endif

    denquot = 0;
    for (k=an;k<=en;k++)
      denquot += compute_denquot(k+start-1,Pm.get_fc_random()[en-an](k-an,0));

    u = log(uniform());
    if (u <= denquot)    // accept
      {
      for(k=an;k<=en;k++)
        {
        h(k-1,0) = Pm.get_fc_random()[en-an](k-an,0);
        if (h(k-1,0) < -20)
          h(k-1,0) = -20;
        beta(k+start-1,0) = exp(h(k-1,0));
        beta(k+start-1,1) = h(k-1,0);
        }
      acceptance++;
      }

    an+=blocksize;
    if (j ==  Pm.get_nrblocks(blocksize)-2)
      en = Pm.get_sizeK();
    else
      en+=blocksize;
    }

  if (unifb==true)
    {
    double help=0;
    while (help < 0.0000005 || help > 0.05)
      help = randnumbers::rand_gamma(a_invgamma+1,1.0/beta(0,0));
    b_invgamma = help;
//    beta(1,0) = help;
    }

  if (type==RW1)
    sigma2 = rand_invgamma(a_invgamma+0.5*(h.rows()-1),
                           b_invgamma+0.5*Pm.compute_quadform(h,0));
  else
    sigma2 = rand_invgamma(a_invgamma+0.5*(h.rows()-2),
                           b_invgamma+0.5*Pm.compute_quadform(h,0));

  if(   (optionsp->get_nriter() > optionsp->get_burnin())
    &&
    (optionsp->get_nriter() % (optionsp->get_step()) == 0)
    )
    {
    sigma2sum += sigma2;
    }

  Gp->updateK(beta.getCol(0));

  FULLCOND::update();
  } // end: function update


void FULLCOND_adaptiv::outresults(void)
  {
  FULLCOND::outresults();

  char hchar = '\\';
  ST::string hstring = "/";
  ST::string pathresultsplus = pathresults.insert_string_char(hchar,hstring);
  ST::string psfile = pathresultsplus.substr(0,pathresultsplus.length()-4) + ".ps";

  optionsp->out("  Results are stored in file\n");
  optionsp->out("  " +   pathresults + "\n");
  optionsp->out("\n");
  optionsp->out("  Results may be visualized using the R function 'plotnonp'\n");
  optionsp->out("\n");

  optionsp->out("  Type for example:\n");
  optionsp->out("\n");
  optionsp->out("  plotnonp(\""+ pathresultsplus + "\")\n");
  optionsp->out("\n");
  optionsp->out("  Estimated variance parameter for variance random walk (sample mean):\n"
                + ST::doubletostring(sigma2sum/optionsp->get_samplesize(),6) + "\n");
  optionsp->out("\n");
  optionsp->out("\n");

  ST::string pathhyper = pathresults.substr(0,pathresults.length()-4)+"_var.res";
  ofstream outreshyper(pathhyper.strtochar());
  outreshyper << "pmean" << endl;
  outreshyper << ST::doubletostring(sigma2sum/optionsp->get_samplesize(),6);

  unsigned i;

  ofstream outrespar(pathresults.strtochar());
  ST::string covname = Gp->get_title();
//  vector<ST::string> effvalues = Gp->get_effectvalues();

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  outrespar << "intnr" << "   ";
  outrespar << covname << "   ";
  outrespar << title << "mean   ";
  outrespar << title << "qu" << l1 << "   ";
  outrespar << title << "qu" << l2 << "   ";
  outrespar << title << "med   ";
  outrespar << title << "qu" << u1 << "   ";
  outrespar << title << "qu" << u2 << "   ";
  outrespar << endl;

     for(i=start;i<beta.rows();i++)
       {
       outrespar << (i+1) << "   ";
//       outrespar << effvalues[i] << "   ";
       outrespar << (i+1) << "   ";
       outrespar << betamean(i,0) << "   ";
       outrespar << betaqu_l1_lower(i,0) << "   ";
       outrespar << betaqu_l2_lower(i,0) << "   ";
       outrespar << betaqu50(i,0) << "   ";
       outrespar << betaqu_l2_upper(i,0) << "   ";
       outrespar << betaqu_l1_upper(i,0) << "   ";
       outrespar << endl;
       }
  }


} // end: namespace MCMC





