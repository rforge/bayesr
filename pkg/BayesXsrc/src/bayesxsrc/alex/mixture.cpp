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





#include "mixture.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------- class FULLCOND_mixture --------------------------------------
//------------------------------------------------------------------------------

void FULLCOND_mixture::init_name(const ST::string & na)
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

    priorassumptions.push_back("Gaussian mixture random effect");
    }


void FULLCOND_mixture::init_names(const vector<ST::string> & na)
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

    priorassumptions.push_back("Gaussian mixture random effect");
    }



FULLCOND_mixture::FULLCOND_mixture(MCMCoptions * o,DISTRIBUTION * dp,
                              FULLCOND_const * fcc,
                              const datamatrix & d, const ST::string & t,
                              const ST::string & fp,const ST::string & pr,
                              const unsigned & nrc,const double & pw,
                              const double & pmm,const double & pmv,
                              const double & pva,const double & pvb,
                              const bool & s,const unsigned & acl,
                              const ST::string & ot,
                              const bool & pvbu,const bool & pvbg,
                              const unsigned & c)
                            : FULLCOND(o,datamatrix(1,1),t,1,1,fp)
  {
  nrcomp = nrc;
  compweight = datamatrix(nrcomp,1,1.0/nrcomp);
  compmean = datamatrix(nrcomp,1,0.0);
  compvar = datamatrix(nrcomp,1,1.0);

  cwprior = datamatrix(nrcomp,1,pw);
  cmpriorm = pmm;
  cmpriorv = pmv;
  cvpriora = pva;
  cvpriorb = pvb;
  nosamples = s;
  cvpriorbunif = pvbu;
  cvpriorbgamma = pvbg;
  if(cvpriorb==1.0 && cvpriorbunif==true && cvpriorbgamma==false)
    cvpriorb=10.0;
  aclag = acl;
  ordertype = ot;

  csize = statmatrix<unsigned>(nrcomp,1,0);
  temp = datamatrix(nrcomp,1,0.0);
  checkorder = false;

  fcconst = fcc;
  fctype = randomeffects;

  column = c;

  likep = dp;
  pathresult = pr;
  pathcurrent = pr;

  index = statmatrix<int>(d.rows(),1);
  index2 = statmatrix<int>(d.rows(),1);
  index.indexinit();
  d.indexsort(index,0,d.rows()-1,0,0);

  unsigned j,k;
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

  setbeta(posbeg.size(),1,0.0);

  compind = datamatrix(beta.rows(),nrcomp+1,0.0);
  double u,inc,incsum;
  for(j=0;j<compind.rows();j++)
    {
    u=uniform();
    incsum=0.0;
    inc=1.0/nrcomp;
    for(k=0;k<nrcomp;k++)
      {
      if( (incsum<u) && (u<=incsum+inc) )
        compind(j,nrcomp)=k+1;
      incsum+=inc;
      }
    }

  ST::string path = samplepath.substr(0,samplepath.length()-4)+"_compparsample.raw";
  cpar_fc = FULLCOND(o,datamatrix(1,1),t+"_cpar_fc",nrcomp,3,path);
  cpar_fc.setflags(MCMC::norelchange | MCMC::nooutput);

  ST::string path2 = samplepath.substr(0,samplepath.length()-4)+"_compindsample.raw";
  cind_fc = FULLCOND(o,datamatrix(1,1),t+"_cind_fc",beta.rows(),nrcomp+1,path2);
  cind_fc.setflags(MCMC::norelchange | MCMC::nooutput);

  transform = likep->get_trmult(column);

  setflags(MCMC::norelchange);
  identifiable = false;
  }


FULLCOND_mixture::FULLCOND_mixture(const FULLCOND_mixture & fc)
                            : FULLCOND(FULLCOND(fc))
  {
  nrcomp = fc.nrcomp;
  compweight=fc.compweight;
  cwprior=fc.cwprior;
  csize=fc.csize;
  compind=fc.compind;
  compmean=fc.compmean;
  compvar=fc.compvar;
  cmpriorm=fc.cmpriorm;
  cmpriorv=fc.cmpriorv;
  cvpriora=fc.cvpriora;
  cvpriorb=fc.cvpriorb;
  nosamples=fc.nosamples;
  aclag=fc.aclag;
  ordertype=fc.ordertype;
  cvpriorbunif=fc.cvpriorbunif;
  cvpriorbgamma=fc.cvpriorbgamma;

  temp=fc.temp;
  checkorder = fc.checkorder;

  fcconst = fc.fcconst;
  cpar_fc = fc.cpar_fc;
  cind_fc = fc.cind_fc;

  likep = fc.likep;
  index = fc.index;
  index2 = fc.index2;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effvalues = fc.effvalues;
  }


const FULLCOND_mixture & FULLCOND_mixture::
         operator=(const FULLCOND_mixture & fc)
  {
  if (this==&fc)
    return *this;

  FULLCOND::operator=(FULLCOND(fc));

  nrcomp = fc.nrcomp;
  compweight=fc.compweight;
  compind=fc.compind;
  compmean=fc.compmean;
  compvar=fc.compvar;
  cmpriorm=fc.cmpriorm;
  cmpriorv=fc.cmpriorv;
  cvpriora=fc.cvpriora;
  cvpriorb=fc.cvpriorb;
  nosamples=fc.nosamples;
  aclag=fc.aclag;
  ordertype=fc.ordertype;
  cvpriorbunif=fc.cvpriorbunif;
  cvpriorbgamma=fc.cvpriorbgamma;

  temp=fc.temp;
  checkorder = fc.checkorder;

  fcconst = fc.fcconst;
  cpar_fc = fc.cpar_fc;
  cind_fc = fc.cind_fc;

  likep = fc.likep;
  index = fc.index;
  index2 = fc.index2;
  posbeg = fc.posbeg;
  posend = fc.posend;
  effvalues = fc.effvalues;

  return *this;
  }



void FULLCOND_mixture::outoptions(void)
  {
  optionsp->out("  OPTIONS FOR MIXTURE RANDOM EFFECT: " + title + "\n",true);
  optionsp->out("\n");

  optionsp->out("  Number of components: " + ST::inttostring(nrcomp) + "\n",true);
  optionsp->out("  Type of Mixture: Gaussian\n",true);
  if(ordertype=="w")
    optionsp->out("  Labelling restriction: ordered weights\n",true);
  else
    optionsp->out("  Labelling restriction: none\n",true);
  optionsp->out("  Prior parameter for component weights: " + ST::doubletostring(cwprior(0,0)) + "\n",true);
  optionsp->out("  Prior parameters for component means:\n",true);
  optionsp->out("    Hyperprior for means: " + ST::doubletostring(cmpriorm,4) + "\n",true);
  optionsp->out("    Hyperprior for variances: " + ST::doubletostring(cmpriorv,4) + "\n",true);
  optionsp->out("  Prior parameters for component variances:\n",true);
  optionsp->out("    Hyperprior for a: " + ST::doubletostring(cvpriora,4) + "\n",true);
  if((cvpriorbunif==false && cvpriorbgamma==false) || (cvpriorbunif==true && cvpriorbgamma==true))
    optionsp->out("    Hyperprior for b: " + ST::doubletostring(cvpriorb,4) + "\n",true);
  if(cvpriorbunif==true && cvpriorbgamma==false)
    optionsp->out("    Hyperprior for b: U(0,"+ ST::doubletostring(cvpriorb,4) +")\n",true);
  if(cvpriorbunif==false && cvpriorbgamma==true)
    optionsp->out("    Hyperprior for b: G(" + ST::doubletostring(cvpriorb,4)+ "," +
                                               ST::doubletostring((100.0*cvpriorb)/(cvpriora*cmpriorv),4) +")\n",true);
  if(nosamples==false)
    optionsp->out("  Samples: yes\n",true);
  else
    optionsp->out("  Samples: no\n",true);
  if(aclag>0)
    optionsp->out("  Autocorrelations: yes\n",true);
  else
    optionsp->out("  Autocorrelations: no\n",true);
  optionsp->out("\n");
  }


bool FULLCOND_mixture::posteriormode(void)
  {
  return true;
  }


double FULLCOND_mixture::centerbeta(void)
  {
  unsigned i;
  double sum = beta.sum(0);
  sum/=beta.rows();
  for (i=0;i<beta.rows();i++)
    beta(i,0) -= sum;
  return sum;
 }


void FULLCOND_mixture::update_weights(void)
  {
  double cwtempsum;
  unsigned k,nzc=0,nzcindex=0;
  for(k=0;k<nrcomp;k++)
    {
    temp(k,0)=rand_gamma(cwprior(k,0)+csize(k,0),1.0);
    if(csize(k,0)>0)
      nzc++;
    }
  cwtempsum=temp.sum(0);
  temp=(1.0/cwtempsum)*temp;
  compweight.assign(temp);

  datamatrix nzcmat=datamatrix(nzc,1,0.0);
  for(k=0;k<nrcomp;k++)
    {
    if(csize(k,0)>0)
      {
      nzcmat(nzcindex,0)=compweight(k,0);
      nzcindex++;
      }
    }

  checkorder=true;
  for(k=0;k<nzcmat.rows()-1;k++)
    {
    if( (nzcmat(k,0)<nzcmat(k+1,0)) )
      checkorder=false;
    }
  }



void FULLCOND_mixture::update(void)
  {
  unsigned i,k;

// Update component indicators, class probabilities
  datamatrix cprob(nrcomp,1,1.0/nrcomp); // probabilities psi_{ik}
  double cptempsum=0.0;
  for(i=0;i<compind.rows();i++)
    {
    // calculate psi_{ik}
    for(k=0;k<nrcomp;k++)
      temp(k,0)=compweight(k,0) * (1.0/(sqrt(compvar(k,0))))
                                * exp(-0.5*(beta(i,0)-compmean(k,0))*(1.0/compvar(k,0))*(beta(i,0)-compmean(k,0)));
    cptempsum=temp.sum(0);
    temp=(1.0/cptempsum)*temp;
    cprob.assign(temp);
    for(k=0;k<nrcomp;k++)
      compind(i,k) = cprob(k,0);

    // sample component indicator
    double cprobsum=0.0, u=uniform();
    for(k=0;k<nrcomp;k++)
      {
      if ( (cprobsum<u) && (u<=cprobsum+cprob(k,0)))
        compind(i,nrcomp) = k+1;
      cprobsum+=cprob(k,0);
      }
    }


// Update component sizes
  csize=statmatrix<unsigned>(nrcomp,1,0);
  for(k=0;k<nrcomp;k++)
    {
      for(i=0;i<beta.rows();i++)
        {
        if(compind(i,nrcomp)==k+1) csize(k,0)+=1;
        }
    }


// Update random effects
  double scaletemp=likep->get_scale(column); // scale parameter sigma_i^2
  double remtemp,revtemp,indobs,sumworkres;
  for(i=0;i<beta.rows();i++)
    {
    likep->add_linearpred(-1.0*beta(i,0),posbeg[i],posend[i],index,0,true);

    indobs=posend[i]-posbeg[i]+1; // number of observations for individual i=X_i'X_i
    sumworkres=0;
    for(k=posbeg[i];k<=posend[i];k++)
      {
      sumworkres += likep->get_response(k,0)-likep->get_linearpred(k,0);  // sum of X_i'(y_i-eta_i)
      }
    revtemp = 1.0/( indobs*(1.0/scaletemp)+(1.0/compvar(compind(i,nrcomp)-1,0)) );
    remtemp = revtemp*( sumworkres*(1.0/scaletemp)+(1.0/compvar(compind(i,nrcomp)-1,0))*compmean(compind(i,nrcomp)-1,0) );
    beta(i,0)=remtemp+sqrt(revtemp)*rand_normal();

    likep->add_linearpred(beta(i,0),posbeg[i],posend[i],index,0,true);
    }



// Update component means, variances
  double mtemp,vtemp,remean,resum,cvb=0.0,cvsum=0.0;
  for(k=0;k<nrcomp;k++)
    {
    resum=0.0;
    for(i=0;i<beta.rows();i++)
      {
      if(compind(i,nrcomp)==k+1)
        resum+=beta(i,0);
      }
    if(csize(k,0)==0)
      remean=0.0;
    else
      remean=resum/csize(k,0);

    vtemp = 1.0 / ( csize(k,0)*(1.0/compvar(k,0))+(1.0/cmpriorv) ) ;
    mtemp = vtemp*( csize(k,0)*(1.0/compvar(k,0))*remean+(1.0/cmpriorv)*cmpriorm );
    compmean(k,0)=mtemp+sqrt(vtemp)*rand_normal();
    }

  if(cvpriorbunif==false && cvpriorbgamma==true)
    {
    for(k=0;k<nrcomp;k++)
      cvsum+=(1.0/compvar(k,0));
    cvb=0.0;
    while (cvb < 1e-3)
       cvb=rand_gamma(cvpriorb+nrcomp*cvpriora,(100.0*cvpriorb)/(cvpriora*cmpriorv)+cvsum);
    }
  if(cvpriorbunif==true && cvpriorbgamma==false)
    {
    for(k=0;k<nrcomp;k++)
      cvsum+=(1.0/compvar(k,0));
    cvb=cvpriorb+1.0;
    while (cvb < 1e-6 || cvb > cvpriorb)
      cvb = rand_gamma(1+nrcomp*cvpriora,cvsum);
    }
  for(k=0;k<nrcomp;k++)
    {
    resum=0.0;
    for(i=0;i<beta.rows();i++)
      {
      if(compind(i,nrcomp)==k+1)
        resum+=(beta(i,0)-compmean(k,0))*(beta(i,0)-compmean(k,0));
      }
    if((cvpriorbunif==true && cvpriorbgamma==false) || (cvpriorbunif==false && cvpriorbgamma==true))
      {
      compvar(k,0)=rand_invgamma(cvpriora+0.5*csize(k,0),cvb+0.5*resum);
      }
    else
      compvar(k,0)=rand_invgamma(cvpriora+0.5*csize(k,0),cvpriorb+0.5*resum);
    }



// Update component weights
  if(ordertype=="w")
    {
    checkorder=false;
    unsigned loopcount=0;
    while(checkorder==false && loopcount<5000)
      {
      ++loopcount;
      update_weights();
      }
    }
  else
    {
    double cwtempsum=0.0;
    for(unsigned k=0;k<nrcomp;k++)
      temp(k,0)=rand_gamma(cwprior(k,0)+csize(k,0),1.0);
    cwtempsum=temp.sum(0);
    temp=(1.0/cwtempsum)*temp;
    compweight.assign(temp);
    }



  double m = centerbeta();
  fcconst->update_intercept(m);


  double * cp_p = cpar_fc.getbetapointer();
  for(k=0;k<nrcomp;k++)
    {
    *cp_p = compweight(k,0);
    cp_p++;
    *cp_p = compmean(k,0)*transform;
    cp_p++;
    *cp_p = compvar(k,0)*transform*transform;
    cp_p++;
    }
  cpar_fc.update();

  double * ci_p = cind_fc.getbetapointer();
  for(i=0;i<beta.rows();i++)
    {
    for(k=0;k<nrcomp+1;k++)
      {
      *ci_p = compind(i,k);
      ci_p++;
      }
    }
  cind_fc.update();

  acceptance++;

  FULLCOND::update();

/*
  ST::string pt2 = pathcurrent.substr(0,pathcurrent.length()-4)+"_test.txt";
  ofstream ot2(pt2.strtochar(),ios::app);
  for(k=0;k<nrcomp;k++)
    ot2 << compvar(k,0)*transform*transform << "   ";
  ot2 << cvb << "   " << cvsum << endl;
//  compvar.prettyPrint(ot2);
*/
 }


void FULLCOND_mixture::outresults(void)
  {
  FULLCOND::outresults();
  cpar_fc.outresults();
  cind_fc.outresults();

  unsigned i,k;

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

  ST::string help =  ST::doubletostring(lower1,4) + "% quant.";
  ST::string levell = help + ST::string(' ',15-help.length());
  help = ST::doubletostring(upper2,4) + "% quant.";
  ST::string levelu = help + ST::string(' ',15-help.length());

  // computing autocorrelations/sampling inefficiency factors of mixture component parameters
  unsigned aclagloc,acsig;
  if(aclag==0 && optionsp->get_samplesize()<=1000)      // to ensure calculation of cdmat for output
    aclagloc=floor(optionsp->get_samplesize()*0.1);
  else if(aclag==0 && optionsp->get_samplesize()>1000)
    aclagloc=100;
  else
    aclagloc=aclag;
  datamatrix acmat = datamatrix(aclagloc,3*nrcomp,0.0);
  for(k=0;k<acmat.cols();k++)
    {
    acmat.putCol(k,cpar_fc.compute_autocorr(aclagloc,0,k));
    }
  datamatrix cdmat(2,acmat.cols(),1.0);

  double siftemp,cv=2.0/sqrt(static_cast<double>(optionsp->get_samplesize()));
  for(i=0;i<acmat.cols();i++)
    {
    acsig=0;
    siftemp=1.0;
    for(k=0;k<acmat.rows();k++)
      {
      acsig++;
      if(k<acmat.rows()-2 && fabs(acmat(k,i))<=cv && fabs(acmat(k+1,i))<=cv && fabs(acmat(k+2,i))<=cv)
        {
        acsig=k;
        break;
        }
      else if(k==acmat.rows()-3 && fabs(acmat(k,i))>cv && fabs(acmat(k+1,i))<=cv && fabs(acmat(k+2,i))<=cv)
        {
        acsig=k+1;
        break;
        }
      else if(k==acmat.rows()-2 && fabs(acmat(k,i))>cv && fabs(acmat(k+1,i))<=cv)
        {
        acsig=k+1;
        break;
        }
      }
    cdmat(0,i)=acsig;

    if(acsig>0)
      {
      for(k=1;k<=acsig;k++)
        siftemp+=2.0*(1.0-k/(acsig+1.0))*acmat(k-1,i);
      cdmat(1,i)=siftemp;
      }
    }

  // calculating results for component indicators
  datamatrix cindres(beta.rows(),nrcomp+1,0.0);
  unsigned cmax=0;
  for(i=0;i<cindres.rows();i++)
    {
    for(k=0;k<nrcomp;k++)
      {
      cindres(i,k)=cind_fc.get_betamean(i,k);
      if(cind_fc.get_betamean(i,k)>cind_fc.get_betamean(i,cmax))
        cmax=k;
      }
    cindres(i,nrcomp)=cmax+1;
    }
  datamatrix cindout(nrcomp,2,0);
  for(i=0;i<cindres.rows();i++)
    {
    for(k=0;k<nrcomp;k++)
      {
      if(cindres(i,nrcomp)==k+1)
        {
        cindout(k,0)+=1;
        cindout(k,1)+=cindres(i,k);
        }
      }
    }

  // Bildschirm-Ausgabe
  optionsp->out("  Sampling diagnostics for mixture component parameters:\n");
  optionsp->out("\n");
  for(k=0;k<nrcomp;k++)
    {
    optionsp->out("   Component " + ST::inttostring(k+1) + "\n");
    optionsp->out("    Parameter   SIF-lag   SIF\n");
    optionsp->out("    Weight         " + ST::doubletostring(cdmat(0,3*k),3) +
                  "      " + ST::doubletostring(cdmat(1,3*k),3) + "\n" );
    optionsp->out("    Mean           " + ST::doubletostring(cdmat(0,3*k+1),3) +
                  "      " + ST::doubletostring(cdmat(1,3*k+1),3) + "\n");
    optionsp->out("    Variance       " + ST::doubletostring(cdmat(0,3*k+2),3) +
                  "      " + ST::doubletostring(cdmat(1,3*k+2),3) + "\n");
    optionsp->out("\n");
    }
  optionsp->out("\n");
  optionsp->out("\n");

  optionsp->out("  Results for mixture component parameters:\n");
  optionsp->out("\n");
  for(k=0;k<nrcomp;k++)
    {
    optionsp->out("   Component " + ST::inttostring(k+1) + "\n");
    optionsp->out("    Parameter    Post. Mean     Std. Dev.      " + levell + "50% quant.     " + levelu + "\n");
    if(nrcomp==1)
      {
      optionsp->out("  " + ST::outresults(5,"Weight  ",cpar_fc.get_betamean(k,0),
                     0,cpar_fc.get_beta_lower1(k,0),cpar_fc.get_betaqu50(k,0),cpar_fc.get_beta_upper1(k,0)) + "\n");
      }
    else
      {
      optionsp->out("  " + ST::outresults(5,"Weight  ",cpar_fc.get_betamean(k,0),
                     sqrt(cpar_fc.get_betavar(k,0)),cpar_fc.get_beta_lower1(k,0),
                     cpar_fc.get_betaqu50(k,0),cpar_fc.get_beta_upper1(k,0)) + "\n");
      }
    optionsp->out("  " + ST::outresults(5,"Mean    ",cpar_fc.get_betamean(k,1),
                     sqrt(cpar_fc.get_betavar(k,1)),cpar_fc.get_beta_lower1(k,1),
                     cpar_fc.get_betaqu50(k,1),cpar_fc.get_beta_upper1(k,1)) + "\n");
    optionsp->out("  " + ST::outresults(5,"Variance",cpar_fc.get_betamean(k,2),
                     sqrt(cpar_fc.get_betavar(k,2)),cpar_fc.get_beta_lower1(k,2),
                     cpar_fc.get_betaqu50(k,2),cpar_fc.get_beta_upper1(k,2)) + "\n");
    optionsp->out("\n");
    }
  optionsp->out("\n");

  optionsp->out("  Results for mixture component indicators:\n");
  optionsp->out("\n");
  optionsp->out("   Component   Assigned subjects   Mean of comp. prob.\n");
  for(k=0;k<nrcomp;k++)
    {
    if(cindout(k,0)>=10)
      optionsp->out("       " + ST::inttostring(k+1) +
                    "               " + ST::doubletostring(cindout(k,0)) +
                    "                  " + ST::doubletostring(cindout(k,1)/cindout(k,0),4) + "\n");
    else if(cindout(k,0)==0)
      optionsp->out("       " + ST::inttostring(k+1) +
                    "                " + ST::doubletostring(cindout(k,0)) +
                    "                  " + ST::doubletostring(0) + "\n");
    else
      optionsp->out("       " + ST::inttostring(k+1) +
                    "                " + ST::doubletostring(cindout(k,0)) +
                    "                  " + ST::doubletostring(cindout(k,1)/cindout(k,0),4) + "\n");
    }
  optionsp->out("\n");
  optionsp->out("\n");


  // Datei-Ausgabe Ergebnisse random effects
  ST::string name = datanames[0];
  ST::string pathre = pathcurrent.substr(0,pathcurrent.length()-4)+".res";
  optionsp->out("  Results for estimated random effects are stored in file\n");
  optionsp->out("  " + pathre + "\n");
  optionsp->out("\n");
  ofstream outres(pathre.strtochar());
  assert(!outres.fail());

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

  for(i=0;i<beta.rows();i++)
    {
    outres << (i+1) << "   ";
    outres << effvalues(i,0) << "   ";
    outres << betamean(i,0) << "   ";
    outres << betaqu_l1_lower(i,0) << "   ";
    outres << betaqu_l2_lower(i,0) << "   ";
    outres << betaqu50(i,0) << "   ";
    outres << betaqu_l2_upper(i,0) << "   ";
    outres << betaqu_l1_upper(i,0) << "   ";
    if (betaqu_l1_lower(i,0) > 0)
      outres << "1   ";
    else if (betaqu_l1_upper(i,0) < 0)
      outres << "-1   ";
    else
      outres << "0   ";
    if (betaqu_l2_lower(i,0) > 0)
      outres << "1   ";
    else if (betaqu_l2_upper(i,0) < 0)
      outres << "-1   ";
    else
      outres << "0   ";
    outres << endl;
    }


  // Datei-Ausgabe Ergebnisse mixture component parameters
  ST::string pathcpar = pathcurrent.substr(0,pathcurrent.length()-4)+"_comppar.res";
  optionsp->out("  Results for mixture component parameters are stored in file\n");
  optionsp->out("  " + pathcpar + "\n");
  optionsp->out("\n");
  ofstream outres2(pathcpar.strtochar());
  assert(!outres2.fail());

  outres2 << "comp   ";
  outres2 << "par   ";
  outres2 << "pmean    ";
  outres2 << "pstddev    ";
  outres2 << "pqu"  << nl1  << "   ";
  outres2 << "pqu"  << nl2  << "   ";
  outres2 << "pmed   ";
  outres2 << "pqu"  << nu1  << "   ";
  outres2 << "pqu"  << nu2  << "   ";
  outres2 << "siflag    ";
  outres2 << "sif";
  outres2 << endl;

  for(k=0;k<nrcomp;k++)
    {
    outres2 << (k+1) << "   ";
    outres2 << "weight" << (k+1) << "   ";
    outres2 << cpar_fc.get_betamean(k,0) << "   ";
    if(nrcomp==1)
      outres2 << 0 << "   ";
    else
      outres2 << sqrt(cpar_fc.get_betavar(k,0)) << "   ";
    outres2 << cpar_fc.get_beta_lower1(k,0) << "   ";
    outres2 << cpar_fc.get_beta_lower2(k,0) << "   ";
    outres2 << cpar_fc.get_betaqu50(k,0) << "   ";
    outres2 << cpar_fc.get_beta_upper2(k,0) << "   ";
    outres2 << cpar_fc.get_beta_upper1(k,0) << "   ";
    outres2 << ST::doubletostring(cdmat(0,3*k),3) << "   ";
    outres2 << ST::doubletostring(cdmat(1,3*k),3) << "   ";
    outres2 << endl;

    outres2 << (k+1) << "   ";
    outres2 << "mean" << (k+1) << "   ";
    outres2 << cpar_fc.get_betamean(k,1) << "   ";
    if(nrcomp==1)
      outres2 << 0 << "   ";
    else
      outres2 << sqrt(cpar_fc.get_betavar(k,1)) << "   ";
    outres2 << cpar_fc.get_beta_lower1(k,1) << "   ";
    outres2 << cpar_fc.get_beta_lower2(k,1) << "   ";
    outres2 << cpar_fc.get_betaqu50(k,1) << "   ";
    outres2 << cpar_fc.get_beta_upper2(k,1) << "   ";
    outres2 << cpar_fc.get_beta_upper1(k,1) << "   ";
    outres2 << ST::doubletostring(cdmat(0,3*k+1),3) << "   ";
    outres2 << ST::doubletostring(cdmat(1,3*k+1),3) << "   ";
    outres2 << endl;

    outres2 << (k+1) << "   ";
    outres2 << "var" << (k+1) << "   ";
    outres2 << cpar_fc.get_betamean(k,2) << "   ";
    if(nrcomp==1)
      outres2 << 0 << "   ";
    else
      outres2 << sqrt(cpar_fc.get_betavar(k,2)) << "   ";
    outres2 << cpar_fc.get_beta_lower1(k,2) << "   ";
    outres2 << cpar_fc.get_beta_lower2(k,2) << "   ";
    outres2 << cpar_fc.get_betaqu50(k,2) << "   ";
    outres2 << cpar_fc.get_beta_upper2(k,2) << "   ";
    outres2 << cpar_fc.get_beta_upper1(k,2) << "   ";
    outres2 << ST::doubletostring(cdmat(0,3*k+2),3) << "   ";
    outres2 << ST::doubletostring(cdmat(1,3*k+2),3) << "   ";
    outres2 << endl;
    }


  // Datei-Ausgabe Ergebnisse mixture component indicators
  ST::string pathcind = pathcurrent.substr(0,pathcurrent.length()-4)+"_compind.res";
  optionsp->out("  Results for mixture component indicators are stored in file\n");
  optionsp->out("  " + pathcind + "\n");
  optionsp->out("\n");
  ofstream outres3(pathcind.strtochar());
  assert(!outres3.fail());

  outres3 << "intnr" << "   ";
  outres3 << name << "   ";
  for(k=0;k<nrcomp;k++)
    outres3 << "cprob" << (k+1) << "   ";
  outres3 << "compind";
  outres3 << endl;

  for(i=0;i<cindres.rows();i++)
    {
    outres3 << (i+1) << "   ";
    outres3 << effvalues(i,0) << "   ";
    for(k=0;k<cindres.cols();k++)
      outres3 << cindres(i,k) << "   ";
    outres3 << endl;
    }


  // Datei-Ausgabe sampling mixture random effects
  if(nosamples==false)
    {
    ST::string file3 = pathcurrent.substr(0,pathcurrent.length()-4) + "_sample.raw";
    get_samples(file3);
    optionsp->out("  Sampling paths for random effects are stored in file\n");
    optionsp->out("  " + file3 + "\n");
    optionsp->out("\n");
    }


  // Datei-Ausgabe sampling/autocorrelations mixture component parameters
  if(nosamples==false)
    {
    ST::string file = pathcurrent.substr(0,pathcurrent.length()-4) + "_comppar_sample.raw";
    cpar_fc.get_samples(file);
    optionsp->out("  Sampling paths for mixture component parameters are stored in file\n");
    optionsp->out("  " + file + "\n");
    optionsp->out("\n");

    ST::string file2 = pathcurrent.substr(0,pathcurrent.length()-4) + "_compind_sample.raw";
    cind_fc.get_samples(file2);
    optionsp->out("  Sampling paths for mixture component indicators are stored in file\n");
    optionsp->out("  " + file2 + "\n");
    optionsp->out("\n");
    }
  if(aclag>0)
    {
    ST::string acfile = pathcurrent.substr(0,pathcurrent.length()-4)+"_comppar_autocor.raw";
    ofstream outac(acfile.strtochar());

    outac << "lag   ";
    for(k=1;k<=nrcomp;k++)
      for(i=1;i<=3;i++)
        outac << "b_" << k << "_" << i << " ";
    outac << endl;

    for(i=0;i<acmat.rows();i++)
      {
      outac << (i+1) << " ";
      for(k=0;k<acmat.cols();k++)
        outac << acmat(i,k) << " ";
      outac << endl;
      }
    optionsp->out("  Autocorrelation functions for mixture component parameters are stored in file\n");
    optionsp->out("  " + acfile + "\n");
    optionsp->out("\n");
    }


  if(nosamples==false && aclag>0)
    {
    optionsp->out("  Sampled parameters may be visualized using the R function 'plotsample'.");
    optionsp->out("  Autocorrelation functions may be visualized using the R function 'plotautocor'.");
    optionsp->out("\n");
    optionsp->out("\n");
    }
  if(nosamples==false && aclag==0)
    {
    optionsp->out("  Sampled parameters may be visualized using the R function 'plotsample'.");
    optionsp->out("\n");
    optionsp->out("\n");
    }
  if(nosamples==true && aclag>0)
    {
    optionsp->out("  Autocorrelation functions may be visualized using the R function 'plotautocor'.");
    optionsp->out("\n");
    optionsp->out("\n");
    }
  }

} // end: namespace MCMC




