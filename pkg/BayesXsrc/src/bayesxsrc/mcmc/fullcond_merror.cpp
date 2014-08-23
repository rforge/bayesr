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



#include "fullcond_merror.h"

//------------------------------------------------------------------------------
//------------ CLASS: FULLCOND_merror implementation of member functions -------
//------------------------------------------------------------------------------


namespace MCMC
{

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // CONSTRUCTOR (Susi)
  // o    : pointer to MCMCoptions object
  // t    : title of the full conditional (for example "fixed effects")
  //        (i.e. number of categories of the response variable)
  // fp   : file path for storing sampled parameters

  fullcond_merror::fullcond_merror(MCMCoptions * o,
//               FULLCOND_nonp_basis * p, const datamatrix & d,        //wird bei IWLS-propsal benötigt
               /*FULLCOND_nonp * p,*/ FULLCOND_nonp_basis * p, DISTRIBUTION * dp, const datamatrix & d,               //wird bei Conditional prior propsal benötigt
               const datamatrix & em, const ST::string & t, const ST::string & fp,
               const double & mvar1, const double & mvar2, const double & arvar,
               const double & arpar1, const double & arpar2, const double & bmean,
               const double & bvar)
                         : FULLCOND(o,d,t,d.rows(),1,fp)
    {

    designp = p;

    effmod = em;

    varcoeff = true;

    likep = dp;

    drows = d.rows();
    dcols = d.cols();

//    double a = 0.75;                // AR-parameter
//    double a2 = a*a;                // AR-parameter^2
//    double sige = 20;               // sigma_e^2 = Varianz AR-Prozess
//    double sigma12 = 20;            // sigma_1^2 = Messfehlervarianz 1
//    double sigma22 = 20;            // sigma_2^2 = Messfehlervarianz 2
    double a = arpar1;                // AR-parameter
    double a2 = arpar2;                // AR-parameter^2
    double sige = arvar;               // sigma_e^2 = Varianz AR-Prozess
    sigma12 = mvar1;                   // sigma_1^2 = Messfehlervarianz 1
    double sigma22 = mvar2;            // sigma_2^2 = Messfehlervarianz 2
    minmerror = 1;                  // Minimum Blocksize
    maxmerror = 1;                  // Maximum Blocksize

    biasmean = bmean;
    biasvar = bvar;

  ST::string path = samplepath.substr(0,samplepath.length()-4);

    // SUSI: Initialize help fullcond
//    whatsoever = FULLCOND(o,datamatrix(1,1),p->get_title()+"_whatsoever",1,1,path);
//    whatsoever.setflags(MCMC::norelchange | MCMC::nooutput);
    double startbias = d.mean(0)-d.mean(1);
    fc_bias = FULLCOND(o,datamatrix(1,1,startbias),p->get_title()+"_bias",1,1,path);
    fc_bias.setflags(MCMC::norelchange | MCMC::nooutput);

    merror_random = datamatrix(drows,1,0);
    randnorm = datamatrix(drows,1,0);

    unsigned i,j;
    datamatrix meandata = datamatrix(drows,1,0);
    for(i=0; i<drows; i++)
      {
       meandata(i,0) = d(i,0)+d(i,1);
       meandata(i,0) /= dcols;
      }

    setbeta(meandata);
    nrpar=drows;

    // Die Matrix P, also die Präzisionsmatrix für die
    // wahren Werte xi muss initialisiert werden.
    // Die Matrix setzt sich aus zwei Teilen zusammen:
    // der Kovarianzmatrix des Messfehlermodells (beobachtete
    // Werte) und der Kovarianzmatrix aus dem Exposure-Modell
    // (AR1-Struktur)

    // First the inverse AR1-Covariance-Matrix:

    datamatrix AR1inv (drows,drows,0);

    // Jetzt muss die Matrix belegt werden

    unsigned k;
    for (k=1;k<drows-1;k++)
        {
        AR1inv(k,k) = 1+a2;
        AR1inv(k,k-1) = -a;
        AR1inv(k,k+1) = -a;
        }
    AR1inv(0,0) = 1;
    AR1inv(drows-1,drows-1) = 1;
    AR1inv(0,1) = -a;
    AR1inv(drows-1,drows-2) = -a;

    // Now we compute the inverse Covariance matrix
    // for the exposure model of xi: Omega^{-1} = (sigma_e^2*AR)^{-1}

    datamatrix Omegainv (drows,drows,0);
    Omegainv = AR1inv /= sige;

    // Als nächstes brauchen wir die Matrix
    // der Messfehlervarianzen -> dabei handelt es sich
    // um eine Diagonalmatrix

    unsigned s,v;
    statmatrix<double> Sigeps1(2*drows,2*drows,0);

    for(s=0; s<drows; s++)
       {
        Sigeps1(s,s) = sigma12;
        }
    for(v=drows; v<2*drows; v++)
        {
        Sigeps1(v,v) = sigma22;
        }

    // Um die richtige Dimension für Sigma^{-1} zu bekommen,
    // brauchen wir die Hilfsmatrix Z = [I_I]'

    datamatrix Z = datamatrix(2*drows,drows,0);

    for(s=0; s<drows; s++)
       {
        Z(s,s) = 1.0;
        }
    for(v=drows; v<2*drows; v++)
        {
        Z(v,v-drows) = 1.0;
        }

    // Now we compute Z' Sigma^{-1} Z = Covariance matrix of measurement model

    datamatrix ZSZ_help;
    datamatrix ZSZ;
    ZSZ_help = datamatrix(drows,2*drows,0);
    ZSZ = datamatrix(drows,drows,0);
    ZSZ_help.mult(Z.transposed(),Sigeps1.inverse());
    ZSZ.mult(ZSZ_help,Z);

    // Now obtain Matrix P = Covariance matrix of true values = Z' Sigma^{-1} Z + Omega^{-1}

    datamatrix P;
    P = datamatrix(drows,drows,0);
    P.plus(ZSZ,Omegainv);

    // Und nun noch die Matrix als envelope-Matrix   // wird nur für IWLS-proposal benötigt

    //precenv = envmatdouble(P);   // wird nur für IWLS-proposal benötigt

    datamatrix Pinv;
    Pinv = P.inverse();

    PABn = datamatrix(drows,1,0);         // wird für Conditional prior proposal benötigt
    for(i=0; i<drows; i++)
       {
        PABn(i,0) = sqrt(Pinv(i,i));
        }

    PABl = datamatrix(drows-1,1,0);       // wird für Conditional prior proposal benötigt
    for(i=0; i<drows-1; i++)
       {
        PABl(i,0) = P(i+1,i);
       }
    PABr = datamatrix(drows-1,1,0);       // wird für Conditional prior proposal benötigt
    for(i=0; i<drows-1; i++)
       {
        PABr(i,0) = P(i,i+1);
        }

    unsigned nrupdate;

    for (i=minmerror;i<=maxmerror;i++)
    {
    nrupdate = drows/i;
    if ((nrupdate*i) < drows)
      nrupdate++;

    matquant.push_back(nrupdate);
    }

    // Berechne rhs

    // Dazu brauchen wir zuerst die Datenmatrix X = (X1_X2), also
    // die Daten untereinander

    datamatrix X = datamatrix(2*drows,1,0);
    for(i=0;i<drows;i++)
        {
        X(i,0) = d(i,0);
        }
    for(j=drows;j<2*drows;j++)
        {
        X(j,0) = d(j-drows,1);
        }

    // Nun berechnen wir Z' Sigeps^{-1} X

    datamatrix mu1_help;
    datamatrix mu1;
    mu1_help = datamatrix(drows,2*drows,0);
    mu1 = datamatrix(drows,1,0);
    mu1_help.addmult(Z.transposed(),Sigeps1.inverse());
    mu1.addmult(mu1_help,X);

    // E(mu_xi)

    datamatrix muxi = datamatrix(drows,1,0);
    unsigned l;
    for (l=0;l<drows;l++)
        {
        muxi(l,0) = 20 + (-0.01*l);
        }

    // Omega^{-1}*muxi

    datamatrix mu2;
    mu2 = datamatrix(drows,1,0);
    mu2.addmult(Omegainv,muxi);

    // Nun rhs = Z'Sigma^{-1}X + Omega^{-1}mu_xi

    rhs = datamatrix(drows,1,0);
    rhs.plus(mu2,mu1);

    mmu = datamatrix(drows,1,0);
    mmu.mult(P.inverse(),rhs);    // P^{-1}*rhs = Erwartungswert der wahren Werte;

    proposalold = meandata;
    proposal = meandata;

    xi = datamatrix(drows,1,0);

    diff = datamatrix(drows,1,0);
     }


// BEGIN: merror
  // CONSTRUCTOR (Thomas)
  fullcond_merror::fullcond_merror(MCMCoptions * o, spline_basis * p,
           DISTRIBUTION * dp, const datamatrix & d, const ST::string & t,
           const ST::string & fp, const ST::string & pres, const double & lk,
           const double & uk, const double & mvar, const bool & disc,
           const int & dig, const unsigned & nb)
           : FULLCOND(o,d,t,d.rows(),1,fp)
  {
  splinep = p;
  likep= dp;
  varcoeff=false;
  unsigned i, j;
  data = d;
  meandata = datamatrix(d.rows(),1,0);
  merror = d.cols();
  for(i=0; i<d.rows(); i++)
    {
    for(j=0; j<merror; j++)
      {
      meandata(i,0) += d(i,j);
      }
    meandata(i,0) /= merror;
    }
  discretize = disc;
  digits = dig;
  nbeta = nb;

  minx = lk+1/pow(10.0,static_cast<double>(digits));
  maxx = uk-1/pow(10.0,static_cast<double>(digits));

  old = meandata;
//  if(discretize)
    old.round(digits,0,1,0,nbeta);
  setbeta(old);

/*  ofstream out0("c:\\temp\\beta.raw");
  beta.prettyPrint(out0);
  out0.close();*/

  logfcold = datamatrix(d.rows(),1,0);
  logfcnew = datamatrix(d.rows(),1,0);

  generrcount=0;
  generrtrial=0;

  pathresults = pres;
  ST::string path = samplepath.substr(0,samplepath.length()-4);

  fc_merrorvar = FULLCOND(o,datamatrix(1,1,mvar),title+"_merror_var",1,1,path+"_merror_var.raw");
  fc_merrorvar.setflags(MCMC::norelchange | MCMC::nooutput);
  double * merrorvarp = fc_merrorvar.getbetapointer();
  *merrorvarp = mvar;

  fc_ximu = FULLCOND(o,datamatrix(1,1,0.0),title+"_truecov_expectation",1,1,path+"_truecov_expectation.raw");
  fc_ximu.setflags(MCMC::norelchange | MCMC::nooutput);
//  double * ximup = fc_ximu.getbetapointer();
//  *ximup = 0.0;

  fc_xivar = FULLCOND(o,datamatrix(1,1,1.0),title+"_truecov_var",1,1,path+"truecov_var.raw");
  fc_xivar.setflags(MCMC::norelchange | MCMC::nooutput);
//  double * xivarp = fc_xivar.getbetapointer();
//  *xivarp = 1;

  index = statmatrix<int>(beta.rows(),1,0);
  for(i=0; i<beta.rows(); i++)
    {
    index(i,0)=i;
    }
  }
// END: merror

  // COPY CONSTRUCTOR

  fullcond_merror::fullcond_merror(const fullcond_merror & m)
    : FULLCOND(FULLCOND(m))
    {
    designp = m.designp;
    likep = m.likep;
    varcoeff = m.varcoeff;
// BEGIN: merror
    splinep = m.splinep;
    minx = m.minx;
    maxx = m.maxx;
    meandata = m.meandata;
    old = m.old;
    index = m.index;
    merror = m.merror;
    currentspline = m.currentspline;
    diffspline = m.diffspline;
    logfcold = m.logfcold;
    logfcnew = m.logfcnew;
    fc_merrorvar = m.fc_merrorvar;
    fc_ximu = m.fc_ximu;
    fc_xivar = m.fc_xivar;
    generrcount = m.generrcount;
    pathresults = m.pathresults;
    discretize = m.discretize;
    digits = m.digits;
    nbeta = m.nbeta;
// END: merror

// BEGIN: Susi
    // SUSI: add help fullcond to copy constructor
    //whatsoever = m.whatsoever;
    fc_bias = m.fc_bias;
    sigma12 = m.sigma12;
    biasmean = m.biasmean;
    biasvar = m.biasvar;
    drows = m.drows;
    minmerror = m.minmerror;
    maxmerror = m.maxmerror;
    matquant = m.matquant;
    merror_random = m.merror_random;
    randnorm = m.randnorm;
    P = m.P;
    PABn = m.PABn;
    PABl = m.PABl;
    PABr = m.PABr;
    rhs = m.rhs;
    proposal = m.proposal;
    proposalold = m.proposalold;
    linold = m.linold;
    linnew = m.linnew;
    diff = m.diff;
    mmu = m.mmu;
    xi = m.xi;
    effmod = m.effmod;
// END: Susi

    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const fullcond_merror & fullcond_merror::operator=(const fullcond_merror & m)
    {
    if (this == &m)
      return *this;
    FULLCOND::operator=(FULLCOND(m));

    designp = m.designp;
    likep = m.likep;
    varcoeff = m.varcoeff;
// BEGIN: merror
    splinep = m.splinep;
    minx = m.minx;
    maxx = m.maxx;
    meandata = m.meandata;
    old = m.old;
    index = m.index;
    merror = m.merror;
    fc_merrorvar = m.fc_merrorvar;
    fc_ximu = m.fc_ximu;
    fc_xivar = m.fc_xivar;
    generrcount = m.generrcount;
    pathresults = m.pathresults;
    discretize = m.discretize;
    digits = m.digits;
    nbeta = m.nbeta;
// END: merror

// BEGIN: Susi
    // SUSI: add help fullcond to copy constructor
    //whatsoever = m.whatsoever;
    fc_bias = m.fc_bias;
    sigma12 = m.sigma12;
    biasmean = m.biasmean;
    biasvar = m.biasvar;
    drows = m.drows;
    minmerror = m.minmerror;
    maxmerror = m.maxmerror;
    matquant = m.matquant;
    merror_random = m.merror_random;
    randnorm = m.randnorm;
    P = m.P;
    PABn = m.PABn;
    PABl = m.PABl;
    PABr = m.PABr;
    rhs = m.rhs;
    proposal = m.proposal;
    proposalold = m.proposalold;
    linold = m.linold;
    linnew = m.linnew;
    diff = m.diff;
    mmu = m.mmu;
    xi = m.xi;
    effmod = m.effmod;
 // END: Susi

    return *this;
    }

  // FUNCTION: update
  // TASK: - stores sampled parameters in file 'samplepath'
  //         storing order: first row, second row, ...

  void fullcond_merror::update(void)
    {
    if(varcoeff)
      {
      datamatrix e = datamatrix(drows,1,0);
//      vector<ST::string> enames;
//      designp->get_effectmatrix(e,enames,0,1,MCMC::fvar_current);
      e = (dynamic_cast <spline_basis*> (designp) )->get_spline();

      // Conditional prior proposal

      unsigned blocks = minmerror + uniform() * (maxmerror-minmerror+1);

      unsigned an = 1;
      unsigned en = blocks;

      double logold;
      double logprop;
      double propold;
      double proposal;
      double proposal1;
      double * workxi;
      double * workxiold;
      unsigned j,k;
      double help;
      double u;
      double acchelp=0;

      statmatrix<int> inde(drows,1);
      for(j=0; j<drows; j++)
        inde(j,0) = j;

      datamatrix proposalmat(drows,1,0);


      for(j=0;j<matquant[blocks-minmerror];j++)
       {

       compute_mu(xi,blocks,an,en);
       compute_proposal(xi,blocks,an,en);

       logold = 0;
       logprop = 0;

       // Compute predictor based on proposed values

       for(k=an-1;k<en;k++)
        {

        logold += likep->loglikelihood(k,k,inde,true);

 	      help = e(k,0);
        help = help/drows;

        propold = proposalold(k,0);
        propold = propold*help;

        proposal = merror_random(k,0);
        proposal = proposal*help;

        proposal1 = proposal - propold;

        likep->addtocurrentcol_single(proposal1,k,column);

        logprop += likep->loglikelihood(k,k,inde,false);
        }

      u = log(uniform());

      if (u <= (logprop-logold))    // accept
        {
        workxiold = proposalold.getV()+an-1;
        workxi = xi.getV()+an-1;
        for(k=an-1;k<en;k++,workxi++,workxiold++)
          {
          *workxiold = *workxi;
          *workxi = merror_random(k,0);
          proposalmat(k,0) = proposal1;
          acchelp++;
          }
        }
      else
        {
        workxiold = proposalold.getV()+an-1;
        workxi = xi.getV()+an-1;
        for(k=an-1;k<en;k++,workxi++,workxiold++)
            {
            *workxi = *workxiold;
            proposalmat(k,0) = 0.0;
            }
        }

      acceptance += acchelp/xi.rows();

      an+=blocks;
      if (j == matquant[blocks-minmerror]-2)
        {
        en = drows;
        }
      else
        en+=blocks;
        }    // end: for(j=0;j<drows;j++)

      likep->addtocurrent(proposalmat);
      likep->swap_linearpred();
//      (dynamic_cast <FULLCOND_nonp*> (designp) )->set_effectmod(xi);
      (dynamic_cast <spline_basis*> (designp) )->update_merror_varcoef(effmod, xi);
      beta.assign(xi);

      FULLCOND::update();

      // SUSI: update help-fullcond object

      // pointer to the current value of the fullcond object
//      double * whatsoeverp = whatsoever.getbetapointer();
      // assign a new value
//      *whatsoeverp = sqrt(rand_gamma(0.5, 0.5));
      // and update
//      whatsoever.update();

      // pointer to the current value of the fullcond object
      double * biasp = fc_bias.getbetapointer();
      // assign a new value
      double empbias=0;
      unsigned i;
      for(i=0; i<drows; i++)
        {
        empbias += data(i,0)-xi(i,0);
        }
      double propvar = 1/biasvar + drows/sigma12;
      double propmean = (biasmean/biasvar + empbias/sigma12)/propvar;

      *biasp = propmean + sqrt(propvar)*rand_normal();
      
      // and update
      fc_bias.update();

      }
    else
      {
      // some miscellaneous variables
      unsigned i, j;
      double u;

      // store results from previous iteration in old
      old = beta;

      // sampling step:
      // standard deviation of measurement error
      double * merrorvarp = fc_merrorvar.getbetapointer();
      double anew, bnew;
/*    currently disabled: measurement error variance is treated as fixed
      anew = 0.001 + 0.5*nbeta*data.cols();
      bnew = 0.0;
      for(i=0; i<nbeta; i++)
        for(j=0; j<data.cols(); j++)
          bnew += (data(i,j)-beta(i,0))*(data(i,j)-beta(i,0));
      bnew = 0.001 + 0.5*bnew;
      *merrorvarp = rand_invgamma(anew,bnew);*/
      double mesd = sqrt(*merrorvarp);
      fc_merrorvar.update();

      // sampling step:
      // expectation of true covariate values
      // sampling step:
      // standard deviation of true covariate values
      double * xivarp = fc_xivar.getbetapointer();
      double * ximup = fc_ximu.getbetapointer();

      anew = 0.001 + 0.5*nbeta;
      bnew = 0.0;
      for(i=0; i<nbeta; i++)
        bnew += (beta(i,0)-*ximup)*(beta(i,0)-*ximup);
      bnew = 0.001 + 0.5*bnew;
      *xivarp = rand_invgamma(anew,bnew);
      double priorsd = sqrt(*xivarp);
      fc_xivar.update();

      double hypersd = 1000;
      double betasum = 0;
      for(i=0; i<nbeta; i++)
        betasum += beta(i,0);
      double muhelp = betasum*hypersd*hypersd/(nbeta*hypersd*hypersd + *xivarp);
      double sdhelp = sqrt(*xivarp *hypersd*hypersd / (nbeta*hypersd*hypersd + *xivarp));
      *ximup = muhelp + sdhelp*rand_normal();
      double priormean = *ximup;
      fc_ximu.update();

      // generate proposed values (random walk proposal according to Berry et al.)
      double * work = beta.getV();
      double * workold = old.getV();
      for(i=0;i<nbeta;i++,work++,workold++)
        {
        *work = *workold + 2*mesd*rand_normal()/(double)merror;
        generrtrial++;
        while(*work<minx || *work>maxx)
          {
          generrcount++;
//          optionsp->out("  WARNING in "+title+":");
//          optionsp->out("          Generated true covariate value out of range!");
          *work = *workold + 2*mesd*rand_normal()/(double)merror;
          generrtrial++;
          }
        }

//      if(discretize)
        beta.round(digits,0,1,0,nbeta);

      // extract current f(x) from spline_basis
      currentspline = splinep->get_spline_merror();

      // call update_merror and compute new values for f(x).
      if(discretize)
        splinep->update_merror_discrete(beta);
      else
        splinep->update_merror(beta);

      diffspline = splinep->get_spline_merror()-currentspline;

      // full conditional for old values
      for(i=0; i<nbeta; i++)
        {
        logfcold(i,0) = likep->loglikelihood(i,i,index,true) - 0.5*((old(i,0)-priormean)/priorsd)*((old(i,0)-priormean)/priorsd);
        for(j=0; j<merror; j++)
          {
          logfcold(i,0) -= 0.5*((data(i,j)-old(i,0))/mesd)*((data(i,j)-old(i,0))/mesd);
          }
        }

      // full conditional for proposed values
      likep->addtocurrent(diffspline);
      for(i=0; i<nbeta; i++)
        {
        logfcnew(i,0) = likep->loglikelihood(i,i,index,false) - 0.5*((beta(i,0)-priormean)/priorsd)*((beta(i,0)-priormean)/priorsd);
        for(j=0; j<merror; j++)
          {
          logfcnew(i,0) -= 0.5*((data(i,j)-beta(i,0))/mesd)*((data(i,j)-beta(i,0))/mesd);
          }
        }

/*      ofstream out1("c:\\temp\\beta.raw");
      beta.prettyPrint(out1);
      out1.close();
      ofstream out2("c:\\temp\\old.raw");
      old.prettyPrint(out2);
      out2.close();
      ofstream out3("c:\\temp\\logfcold.raw");
      logfcold.prettyPrint(out3);
      out3.close();
      ofstream out4("c:\\temp\\logfcnew.raw");
      logfcnew.prettyPrint(out4);
      out4.close();
      ofstream out5("c:\\temp\\diffspline.raw");
      diffspline.prettyPrint(out5);
      out5.close();
      ofstream out6("c:\\temp\\currentspline.raw");
      currentspline.prettyPrint(out6);
      out6.close(); */

      datamatrix accmat(nbeta,1,0);

      // compute acceptance probabilities; overwrite non-accepted values with old values
      for(i=0; i<nbeta; i++)
        {
        u = log(uniform());
        nrtrials++;
//        double test = (logfcnew(i,0)-logfcold(i,0));
        if(u <= (logfcnew(i,0)-logfcold(i,0)))
          {
          acceptance += 1.0;
          accmat(i,0) = 1;
          }
        else
          {
          beta(i,0) = old(i,0);
          }
        }

      // update spline_basis with the final information
      if(discretize)
        splinep->update_merror_discrete(beta);
      else
        splinep->update_merror(beta);

      // Update the linear predictor
      diffspline = splinep->get_spline_merror()-currentspline;
      likep->addtocurrent(diffspline);
      likep->swap_linearpred();

/*      ofstream out1("c:\\temp\\beta.raw");
      beta.prettyPrint(out1);
      out1.close();
      ofstream out2("c:\\temp\\old.raw");
      old.prettyPrint(out2);
      out2.close();
      ofstream out7("c:\\temp\\diffspline2.raw");
      diffspline.prettyPrint(out7);
      out7.close();
      ofstream out8("c:\\temp\\accmat.raw");
      accmat.prettyPrint(out8);
      out8.close();*/

      FULLCOND::update();
      }
    }

  void fullcond_merror::compute_mu(datamatrix & xi,
            const unsigned & blocks,const unsigned & a,const unsigned & b)
  {
  unsigned t;
  datamatrix xihelp1 = datamatrix(blocks,1,0);
  datamatrix xihelp2 = datamatrix(blocks,1,0);
  datamatrix xihelp3 = datamatrix(blocks,1,0);
  datamatrix xihelp4 = datamatrix(blocks,1,0);
  datamatrix xihelp5 = datamatrix(blocks,1,0);
  datamatrix xihelp6 = datamatrix(blocks,1,0);

  if (a==1)
    {
    for (t=a-1; t<b; t++)
      {
      xihelp1(0,0) = PABn(t,0);
      xihelp2(0,0) = PABr(t,0);
      xihelp1.elemmult(xihelp2);
      xihelp2(0,0) = proposalold(t+1,0) - mmu(t+1,0);
      xihelp1.elemmult(xihelp2);
      xi(t,0) = mmu(t,0)-xihelp1;
      }
    }
  else if (b == drows)
    {
    for (t=a-1; t<b; t++)
      {
      xihelp1(0,0) = PABn(t,0);
      xihelp2(0,0) = PABl(t-1,0);
      xihelp1.elemmult(xihelp2);
      xihelp3(0,0) = proposalold(t-1,0) - mmu(t-1,0);
      xihelp1.elemmult(xihelp3);
      xi(t,0) = mmu(t,0)-xihelp2;
      }
    }
  else
    {
    for (t=a-1; t<b; t++)
      {
      xihelp1(0,0) = PABn(t,0);
      xihelp2(0,0) = PABl(t,0);
      xihelp3(0,0) = PABr(t,0);
      xihelp4(0,0) = proposalold(t-1,0) - mmu(t-1,0);
      xihelp5(0,0) = proposalold(t+1,0) - mmu(t+1,0);
      xihelp2.elemmult(xihelp4);
      xihelp3.elemmult(xihelp5);
      xihelp6.plus(xihelp2,xihelp3);
      xihelp1.elemmult(xihelp6);
      xi(t,0) = mmu(t,0) - xihelp1;
      }
    }
  }

  void fullcond_merror::compute_proposal(const datamatrix & xi, const unsigned & blocks,
                               const unsigned & a,const unsigned b)
  {

  unsigned l = b-a+1;
  register unsigned i;

  double * merror_randwork = merror_random.getV()+a-1;
  double * randnormwork = randnorm.getV()+a-1;
  double * PABrootwork = PABn.getV()+a-1;
  double * workxi = xi.getV()+a-1;

  for (i=0;i<l;i++,merror_randwork++,PABrootwork++,randnormwork++)
    {
    *randnormwork = rand_normal();
    *merror_randwork = 0;
    *merror_randwork += *PABrootwork * *randnormwork;
    *merror_randwork += *workxi;
    }

  }


  // FUNCTION: posteriormode
  // TASK: computes the posterior mode

  bool fullcond_merror::posteriormode(void)
    {
      return true;
    }

  bool fullcond_merror::posteriormode_converged(const unsigned & itnr)
    {
      return true;
    }

  void fullcond_merror::posteriormode_set_beta_mode(void)
    {
    }

  // FUNCTION: outoptions
  // TASK: writes estimation options (hyperparameters, etc.) to outputstream

  void fullcond_merror::outoptions(void)
    {
    }

  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file

  void fullcond_merror::outresults(void)
    {
    FULLCOND::outresults();

    if(!varcoeff)
      {
    ST::string l1 = ST::doubletostring(lower1,4);
    ST::string l2 = ST::doubletostring(lower2,4);
    ST::string u1 = ST::doubletostring(upper1,4);
    ST::string u2 = ST::doubletostring(upper2,4);

    ST::string nl1 = ST::doubletostring(lower1,4);
    ST::string nl2 = ST::doubletostring(lower2,4);
    ST::string nu1 = ST::doubletostring(upper1,4);
    ST::string nu2 = ST::doubletostring(upper2,4);
    nl1 = nl1.replaceallsigns('.','p');
    nl2 = nl2.replaceallsigns('.','p');
    nu1 = nu1.replaceallsigns('.','p');
    nu2 = nu2.replaceallsigns('.','p');

    ST::string vstr;

    ofstream ou(pathresults.strtochar());

    unsigned i;
    ou << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
    nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
    for(i=0; i<nbeta; i++)
      {
      ou << betamean(i,0) << "  ";
      ou << (betavar(i,0)<0.0?0.0:sqrt(betavar(i,0))) << "  ";
      ou << betaqu_l1_lower(i,0) << "  ";
      ou << betaqu_l2_lower(i,0) << "  ";
      ou << betaqu50(i,0) << "  ";
      ou << betaqu_l2_upper(i,0) << "  ";
      ou << betaqu_l1_upper(i,0) << "  ";
      ou << betamin(i,0) << "  ";
      ou << betamax(i,0) << "  " << endl;
      }

    optionsp->out("  Results for the covariate values are stored in file\n");
    optionsp->out("  " + pathresults + "\n");

    optionsp->out("\n");

    if(generrcount>0)
      {
      optionsp->out("  WARNING: "+ST::doubletostring((double)generrcount/(double)generrtrial,5)+"% of the generated covariate values were out of the specified range!");
      optionsp->out("           Consider making the grid more wide.\n\n");
      }

// Variance of the measurement error

    fc_merrorvar.outresults();
//    ST::string pathhelp = pathresults.substr(0,pathresults.length()-7)+"merrorvar_sample.raw";
//    fc_merrorvar.get_samples(pathhelp);

    optionsp->out("\n");
    optionsp->out("  "+fc_merrorvar.get_title()+"\n");
    optionsp->out("\n");
    optionsp->out("\n");

    vstr = "  Mean:         ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_betamean(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  Std. dev.:    ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring((fc_merrorvar.get_betavar(0,0)<0.0?0.0:sqrt(fc_merrorvar.get_betavar(0,0))),6) + "\n");

    vstr = "  " + l1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_beta_lower1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + l2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_beta_lower2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  50% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_betaqu50(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_beta_upper2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_merrorvar.get_beta_upper1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    optionsp->out("\n");

    ST::string merrorvar_pathresults = pathresults.substr(0,pathresults.length()-4) + "_merror_var.res";

    ofstream ou1(merrorvar_pathresults.strtochar());

    ou1 << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
    nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
    ou1 << fc_merrorvar.get_betamean(0,0) << "  ";
    ou1 << (fc_merrorvar.get_betavar(0,0)<0.0?0.0:sqrt(fc_merrorvar.get_betavar(0,0))) << "  ";
    ou1 << fc_merrorvar.get_beta_lower1(0,0) << "  ";
    ou1 << fc_merrorvar.get_beta_lower2(0,0) << "  ";
    ou1 << fc_merrorvar.get_betaqu50(0,0) << "  ";
    ou1 << fc_merrorvar.get_beta_upper2(0,0) << "  ";
    ou1 << fc_merrorvar.get_beta_upper1(0,0) << "  ";
    ou1 << fc_merrorvar.get_betamin(0,0) << "  ";
    ou1 << fc_merrorvar.get_betamax(0,0) << "  " << endl;

    optionsp->out("  Results for the variance of the measurement error are also stored in file\n");
    optionsp->out("  " + merrorvar_pathresults + "\n");

    optionsp->out("\n");

// mean of the true covariate values

    fc_ximu.outresults();
//    pathhelp = pathresults.substr(0,pathresults.length()-7)+"ximu_sample.raw";
//    fc_ximu.get_samples(pathhelp);

    optionsp->out("\n");
    optionsp->out("  "+fc_ximu.get_title()+"\n");
    optionsp->out("\n");
    optionsp->out("\n");

    vstr = "  Mean:         ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_betamean(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  Std. dev.:    ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring((fc_ximu.get_betavar(0,0)<0.0?0.0:sqrt(fc_ximu.get_betavar(0,0))),6) + "\n");

    vstr = "  " + l1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_beta_lower1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + l2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_beta_lower2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  50% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_betaqu50(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_beta_upper2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_ximu.get_beta_upper1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    optionsp->out("\n");

    ST::string ximu_pathresults = pathresults.substr(0,pathresults.length()-4) + "_truecov_expectation.res";

    ofstream ou2(ximu_pathresults.strtochar());

    ou2 << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
    nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
    ou2 << fc_ximu.get_betamean(0,0) << "  ";
    ou2 << (fc_ximu.get_betavar(0,0)<0.0?0.0:sqrt(fc_ximu.get_betavar(0,0))) << "  ";
    ou2 << fc_ximu.get_beta_lower1(0,0) << "  ";
    ou2 << fc_ximu.get_beta_lower2(0,0) << "  ";
    ou2 << fc_ximu.get_betaqu50(0,0) << "  ";
    ou2 << fc_ximu.get_beta_upper2(0,0) << "  ";
    ou2 << fc_ximu.get_beta_upper1(0,0) << "  ";
    ou2 << fc_ximu.get_betamin(0,0) << "  ";
    ou2 << fc_ximu.get_betamax(0,0) << "  " << endl;

    optionsp->out("  Results for the expectation of the true covariate values are also stored in file\n");
    optionsp->out("  " + ximu_pathresults + "\n");

    optionsp->out("\n");

// variance of the true covariate values

    fc_xivar.outresults();
//    pathhelp = pathresults.substr(0,pathresults.length()-7)+"xivar_sample.raw";
//    fc_xivar.get_samples(pathhelp);

    optionsp->out("\n");
    optionsp->out("  "+fc_xivar.get_title()+"\n");
    optionsp->out("\n");
    optionsp->out("\n");

    vstr = "  Mean:         ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_betamean(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  Std. dev.:    ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring((fc_xivar.get_betavar(0,0)<0.0?0.0:sqrt(fc_xivar.get_betavar(0,0))),6) + "\n");

    vstr = "  " + l1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_beta_lower1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + l2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_beta_lower2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  50% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_betaqu50(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_beta_upper2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    vstr = "  " + u2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_xivar.get_beta_upper1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length()) + "\n");

    optionsp->out("\n");

    ST::string xivar_pathresults = pathresults.substr(0,pathresults.length()-4) + "_truecov_var.res";

    ofstream ou3(xivar_pathresults.strtochar());

    ou3 << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
    nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
    ou3 << fc_xivar.get_betamean(0,0) << "  ";
    ou3 << (fc_xivar.get_betavar(0,0)<0.0?0.0:sqrt(fc_xivar.get_betavar(0,0))) << "  ";
    ou3 << fc_xivar.get_beta_lower1(0,0) << "  ";
    ou3 << fc_xivar.get_beta_lower2(0,0) << "  ";
    ou3 << fc_xivar.get_betaqu50(0,0) << "  ";
    ou3 << fc_xivar.get_beta_upper2(0,0) << "  ";
    ou3 << fc_xivar.get_beta_upper1(0,0) << "  ";
    ou3 << fc_xivar.get_betamin(0,0) << "  ";
    ou3 << fc_xivar.get_betamax(0,0) << "  " << endl;

    optionsp->out("  Results for the variance of the true covariate values are also stored in file\n");
    optionsp->out("  " + xivar_pathresults + "\n");

    optionsp->out("\n");
      }
    }

  // FUNCTION: reset
  // TASK: resets all parameters

  void fullcond_merror::reset(void)
    {
    }

//  vector<ST::string> & fullcond_merror::get_results_latex(void)
//    {
//      return vector<ST::string>();
//    }

} // end: namespace MCMC



