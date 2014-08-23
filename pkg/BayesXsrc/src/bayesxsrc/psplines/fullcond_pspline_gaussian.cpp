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



#include "fullcond_pspline_gaussian.h"

using std::ifstream;
using std::ios;
using std::min;

namespace MCMC
{

  // CONSTRUCTOR 1  (for additive models)

FULLCOND_pspline_gaussian::FULLCOND_pspline_gaussian(MCMCoptions * o,
                      DISTRIBUTION * dp, FULLCOND_const * fcc, const datamatrix & d,
                      const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                      const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres, const bool & deriv,
                      const double & l, const int & gs, const bool & diag,
                      const double & lk, const double & uk, const double & lg,
                      const double & ug, const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,degr,kp,gs,fp,pres,deriv,lk,uk,lg,ug,c)
  {

  diagtransform = diag;

  if(monotone == "increasing")
    increasing = true;
  else if(monotone == "decreasing")
    decreasing = true;

  hierarchical = false;

  samplecentered = false;

  varcoeff = false;

  transform = likep->get_trmult(c);

  compute_betaweight();

  unsigned i;

  lambda = l;

  make_index(d);
  make_Bspline(d);

  compute_Kweights();

  if(predictright || predictleft)
    {
    datamatrix help = betaweight;
    betaweight = datamatrix(nrpar,1,0);
    for(i=0;i<nrparpredictleft;i++)
      betaweight(i,0) = 0.0;
    for(;i<nrpar-nrparpredictright;i++)
      betaweight(i,0) = help(i-nrparpredictleft,0);
    for(;i<nrpar;i++)
      betaweight(i,0) = 0.0;
    }

// index2 initialisieren

  index2.push_back(index(0,0));
  for(i=1;i<likep->get_nrobs();i++)
    index2.push_back(index(i,0)-index(i-1,0));

  init_fchelp(d);

// Varianz für die priori des linearen Anteils bei hierachical centering

  double priorvar = 100;

// Penalty Matrix erstellen

  if (type == RW1)
    {
    K = Krw1band(weight);
    Kenv = Krw1env(weight);
    rankK = nrpar-nrparpredictleft-nrparpredictright-1;
    if(hierarchical)
      {
      *K.getdiagpointer() += 1.0/priorvar;
      *Kenv.getDiagIterator() += 1.0/priorvar;
      }
    }
  else if (type == RW2)
    {
    K = Krw2band(weight);
    Kenv = Krw2env(weight);
    rankK = nrpar-nrparpredictleft-nrparpredictright-2;
    if(hierarchical)
      {
      vector<double>::iterator it;
      double * p;
      p = K.getdiagpointer();
      it = Kenv.getDiagIterator();
      *p += 1.0/priorvar;
      *it += 1.0/priorvar;
      p++;it++;
      *p += 1.0/priorvar;
      *it += 1.0/priorvar;
      }
    }

  if(predictleft || predictright)
    change_K();

  standnormal = datamatrix(nrpar,1,0);

  betaprop = datamatrix(nrpar,1,0);

  XX_env = envmatdouble(bandmatdouble(nrpar,degree,0));
  compute_XWXenv(likep->get_weight());

  if (type==RW1)
    prec_env = envmatdouble(0.0,nrpar,degree>1?degree:1);
  else if (type == RW2)
    prec_env = envmatdouble(0.0,nrpar,degree>2?degree:2);

  mu = datamatrix(likep->get_nrobs(),1,0);
  muy = datamatrix(nrpar,1,0);
  betahelp = muy;

// gamma für hierarchical centering initialisieren

  if(hierarchical)
    {
    double h = knot[1]-knot[0];
    double gamma0 = -0.5*(nrpar-1)*h;
    gamma = datamatrix(nrpar,1,0);
    for(i=0;i<nrpar;i++)
      gamma(i,0) = gamma0 + i*h;
    }

  identifiable = false;

  }


// CONSTRUCTOR 2  (for varying coefficients term)

FULLCOND_pspline_gaussian::FULLCOND_pspline_gaussian(MCMCoptions * o, DISTRIBUTION * dp,
                      FULLCOND_const * fcc,const datamatrix & effmod, const datamatrix & intact,
                      const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                      const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres, const bool & deriv,
                      const double & l, const int & gs,bool  ce, const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,degr,kp,gs,fp,pres,false,0.0,0.0,0.0,0.0,c)
  {

  assert(effmod.rows() == intact.rows());

  data_forfixed = intact;

  interactvar=intact;

  if(monotone == "increasing")
    increasing = true;
  else if(monotone == "decreasing")
    decreasing = true;

  hierarchical = false;

  samplecentered = false;

  varcoeff = true;
  diagtransform = false;

  transform = likep->get_trmult(c);

  compute_betaweight();

  unsigned i;

  lambda = l;

  make_index(effmod,intact);
  make_Bspline(effmod);
  make_BS(intact);

  data = effmod;

  compute_Kweights();

  if(predictright || predictleft)
    {
    datamatrix help = betaweight;
    betaweight = datamatrix(nrpar,1,0);
    for(i=0;i<nrparpredictleft;i++)
      betaweight(i,0) = 0.0;
    for(;i<nrpar-nrparpredictright;i++)
      betaweight(i,0) = help(i-nrparpredictleft,0);
    for(;i<nrpar;i++)
      betaweight(i,0) = 0.0;
    }

// index2 initialisieren

  index2.push_back(index(0,0));
  for(i=1;i<likep->get_nrobs();i++)
    index2.push_back(index(i,0)-index(i-1,0));

  init_fchelp(effmod);

// Penalty Matrix erstellen

  if (type == RW1)
    {
    K = Krw1band(weight);
    Kenv = Krw1env(weight);
    rankK = nrpar-nrparpredictleft-nrparpredictright-1;
    }
  else if (type == RW2)
    {
    K = Krw2band(weight);
    Kenv = Krw2env(weight);
    rankK = nrpar-nrparpredictleft-nrparpredictright-2;
    }

  if(predictleft || predictright)
    change_K();

  standnormal = datamatrix(nrpar,1,0);

  XX_env = envmatdouble(bandmatdouble(nrpar,degree,0));
  compute_XWXenv(likep->get_weight());

  if (type==RW1)
    prec_env = envmatdouble(0.0,nrpar,degree>1?degree:1);
  else if (type == RW2)
    prec_env = envmatdouble(0.0,nrpar,degree>2?degree:2);

  mu = datamatrix(likep->get_nrobs(),1,0);
  muy = datamatrix(nrpar,1,0);
  betahelp = muy;

// gamma für hierarchical centering initialisieren

  if(hierarchical)
    {
    double h = knot[1]-knot[0];
    double gamma0 = -0.5*(nrpar-1)*h;
    gamma = datamatrix(nrpar,1,0);
    for(i=0;i<nrpar;i++)
      gamma(i,0) = gamma0 + i*h;
    }

  if (ce==false)
    identifiable = true;
  else
    identifiable = false;

  }

  // COPY CONSTRUCTOR

FULLCOND_pspline_gaussian::FULLCOND_pspline_gaussian(const FULLCOND_pspline_gaussian & fc)
  : spline_basis(spline_basis(fc))
  {
  diagtransform = fc.diagtransform;
  samplecentered = fc.samplecentered;
  hierarchical = fc.hierarchical;
  lineff = fc.lineff;
  lineffsum = fc.lineffsum;
  lineffsamples = fc.lineffsamples;
  gamma = fc.gamma;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const FULLCOND_pspline_gaussian & FULLCOND_pspline_gaussian::operator=(
                                            const FULLCOND_pspline_gaussian & fc)
  {
  if (this == &fc)
    return *this;
  spline_basis::operator=(spline_basis(fc));

  diagtransform = fc.diagtransform;
  samplecentered = fc.samplecentered;
  hierarchical = fc.hierarchical;
  lineff = fc.lineff;
  lineffsum = fc.lineffsum;
  lineffsamples = fc.lineffsamples;
  gamma = fc.gamma;

  return *this;
  }


void FULLCOND_pspline_gaussian::update(void)
  {
// XWX initialisieren
  if (optionsp->get_nriter()==1)
    compute_XWXenv(likep->get_weight());
// sigma2 so setzen, dass scale/sigma2 = const
  if(lambdaconst == true)
    sigma2 = likep->get_scale(column)/lambda;

  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);

  unsigned i;

  if(increasing || decreasing)
    {
    update_isotonic();
    }
  else if(diagtransform)
    {
    update_diagtransform();
    }
  else
    {

    if(samplecentered)
      likep->substr_linearpred(spline);               // eta = eta - spline
    else
      subtr_spline();                                 // eta = eta - spline + intercept

    if(changingweight)                                // für t-link
      compute_XWXenv(likep->get_weight());

    double scaleinv = 1.0/likep->get_scale(column);   // scaleinv = 1/scale

    prec_env.addto(XX_env,Kenv,scaleinv,1.0/sigma2);  // prec_env = (scaleinv*XX_env + 1.0/sigma*Kenv)

    double * work = standnormal.getV();               // standnormal ~ N(0,I)
    for(i=0;i<nrpar;i++,work++)
      *work = rand_normal();

    likep->compute_respminuslinpred(mu,column);       // nicht ändern wegen multgaussian
    compute_XWtildey(likep->get_weight(),scaleinv);   // muy = scaleinv * X'W*mu

    if(hierarchical)                                  // muy = muy + lineff/sigma2 * K'gamma
      {
      datamatrix muy2 = datamatrix(nrpar,1,0);
      K.mult(gamma,muy2);
      muy2 = (lineff/sigma2)*muy2;
      muy.plus(muy,muy2);
      }

    beta.assign(standnormal);
    prec_env.solve(muy,betahelp);
    prec_env.solveU(beta,betahelp);                   // betahelp = P^(-1) * muy
                                                      // beta ~ N(betahelp,P^(-1))
    if(predictright || predictleft)
      update_prediction();

    if(hierarchical)                                  // linearen Anteil updaten
      {
      double N01 = rand_normal();
      double lineffsigma2 = Kenv.compute_quadform(gamma,0)/sigma2;
      double lineffmu = 0.0;

      datamatrix res = datamatrix(nrpar,1,0);
      K.mult(beta,res);
      for(i=0;i<nrpar;i++)
        lineffmu += gamma(i,0)*res(i,0);
      lineffmu /= sigma2;

      lineffmu /= lineffsigma2;

      lineff = lineffmu + sqrt(lineffsigma2)*N01;

      if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
          ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
        {
        lineffsum += transform*lineff;
        lineffsamples.push_back(transform*lineff);
        }
      }

    add_linearpred_multBS();

    }    // ENDE: update

  if(center)                                          // zentrieren
    {
    if(samplecentered)
      {
      sample_centered_env(beta);
      }
    else
      {
      if(diagtransform)
        compute_intercept(betaprop);
      else
        compute_intercept();

      if (varcoeff)
        fcconst->update_fix_varcoeff(intercept,datanames[1]);
      else
        fcconst->update_intercept(intercept);

      }
    }


  if(contourprob >= 0)                                // für contour probabilities
    {
    for(i=0;i<nrpar;i++)
      beta(i,0) -= intercept;

//    write_contour();            // für S-PLUS
//    FULLCOND_nonp_basis::write_contour(betahelp,1.0/likep->get_scale(column),1.0/sigma2);
    FULLCOND_nonp_basis::write_contour(betahelp,1.0/likep->get_scale(column),1.0/sigma2,
                        XX_env.compute_quadform(beta,0),Kenv.compute_quadform(beta,0),
                        prec_env.compute_quadform(betahelp,0),prec_env.getLogDet(),&prec_env);
    fc_contour.update();
    }

  acceptance++;

  if(interaction == false)
    {
// wird bei interaction==true in der full conditional des Interaktionseffekts gemacht
// spline in fchelp schreiben
    if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
        ((optionsp->get_nriter()-optionsp->get_burnin()-1) % optionsp->get_step() == 0) )
      {
      if(samplecentered)
        {
        write_spline();
        write_derivative();
        }
      else
        {
        double * splinep;
        double * fchelpbetap = fchelp.getbetapointer();

        if(gridsize < 0)                              // alle verschiedene Beobachtungen
          {
          if(varcoeff)
            multBS(splinehelp,beta);

          vector<int>::iterator freqwork = freqoutput.begin();
          int * workindex = index.getV();
          for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
            {
            if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
              {
              if(varcoeff)
                {
                *fchelpbetap = splinehelp(i,0) - intercept;
                }
              else
                *fchelpbetap = spline(*workindex,0) - intercept;
              fchelpbetap++;
              }
            }
          }
        else                                          // Gitterpunkte
          {
          multDG(splinehelp,beta);
          splinep = splinehelp.getV();
          for(i=0;i<gridsize;i++,fchelpbetap++,splinep++)
            *fchelpbetap = *splinep - intercept;
          }

        write_derivative();                           // 1. Ableitung rausschreiben
        }
      }

    if(derivative)
      fcderivative.update();

    fchelp.update();
    FULLCOND::update();

    }  // end: if(interaction == false)

  }


void FULLCOND_pspline_gaussian::update_isotonic(void)
  {

  unsigned i,j;
  double help, m;

  subtr_spline();

  if(changingweight)                                  // für t-link
    compute_XWXenv(likep->get_weight());

  double scaleinv = 1.0/likep->get_scale(column);
  prec_env.addto(XX_env,Kenv,scaleinv,1.0/sigma2);

  likep->compute_respminuslinpred(mu,column);         // nicht ändern wegen multgaussian
  compute_XWtildey(likep->get_weight(),scaleinv);

  int count = 0;
//  int maxit = 100;
  int maxit = 20;
// interner Gibbs-Sampler mit maxit Iterationen (siehe Geweke, Robert)
  while(count < maxit)
    {

    help = 0.0;
    for(j=1;j<nrpar;j++)
      help += prec_env(0,j)*beta(j,0);
    m = (muy(0,0) - help)/prec_env(0,0);

    beta(0,0) = sample_monotonic(0,m,sqrt(1.0/prec_env(0,0)));

    for(i=1;i<nrpar-1;i++)
      {

      help = 0.0;
      for(j=0;j<i;j++)
        help += prec_env(i,j)*beta(j,0);
      for(j=i+1;j<nrpar;j++)
        help += prec_env(i,j)*beta(j,0);
      m = (muy(i,0) - help)/prec_env(i,i);

      beta(i,0) = sample_monotonic(i,m,sqrt(1.0/prec_env(i,i)));
      }

    help = 0.0;
    for(j=0;j<nrpar-1;j++)
      help += prec_env(nrpar-1,j)*beta(j,0);
    m = (muy(nrpar-1,0) - help)/prec_env(nrpar-1,nrpar-1);

    beta(nrpar-1,0) = sample_monotonic(nrpar-1,m,sqrt(1.0/prec_env(nrpar-1,nrpar-1)));

    count++;
    }  // END: while

  add_linearpred_multBS();

  }


void FULLCOND_pspline_gaussian::update_diagtransform(void)
  {

  unsigned i;

  if(optionsp->get_nriter() == 1)
    {
/*
// L rausschreiben (für S-Plus)
    XX = bandmatdouble(nrpar,degree,0);
    compute_XWX(likep->get_weight());

    bandmatdouble LL = XX;
    LL.decomp();
    datamatrix L = datamatrix(nrpar,nrpar,0);
    LL.getL(L);

    ofstream outL("c:\\cprog\\testmcmc\\L.raw");
    L.prettyPrint(outL);
    outL.close();
*/
// B, Kenv einlesen (aus S-Plus)
    B = datamatrix(nrpar,nrpar);
    ifstream in("c:\\cprog\\testmcmc\\G.raw");
    B.prettyScan(in);
    in.close();

    betahelp = datamatrix(1,nrpar);
    ifstream in2("c:\\cprog\\testmcmc\\D.raw");
    betahelp.prettyScan(in2);
    in2.close();

    Kenv = envmatdouble(bandmatdouble(betahelp.transposed()));
    rankK = nrpar;

/*
datamatrix XB = datamatrix(likep->get_nrobs(),nrpar);
getX(XB);
XB = XB*B;
ofstream out("c:\\cprog\\testmcmc\\XB.raw");
XB.prettyPrint(out);
out.close();
*/
    }

  subtr_spline();

  double scaleinv = 1.0/likep->get_scale(column);

  likep->compute_respminuslinpred(mu,column);      // nicht ändern wegen multgaussian
  compute_XWtildey(likep->get_weight(),scaleinv);
  muy = B.transposed()*muy;                        // B ist volle Matrix (nrpar x nrpar)

  double prec_ii;                                  // prec_env ist diagonal!
  for(i=0;i<beta.rows();i++)
    {
    prec_ii = scaleinv + betahelp(0,i)/sigma2;
    beta(i,0) = rand_normal()/sqrt(prec_ii) + muy(i,0)/prec_ii;
    }

  betaprop.mult(B,beta);                           // beta rücktransformieren
  add_linearpred_multBS(betaprop);

  }


bool FULLCOND_pspline_gaussian::posteriormode_converged(const unsigned & itnr)
  {
  return likep->posteriormode_converged_fc(beta,beta_mode,itnr);
  }


bool FULLCOND_pspline_gaussian::posteriormode(void)
  {

  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);

  unsigned i;

  if(samplecentered)
    likep->substr_linearpred(spline);
  else
    subtr_spline();

  compute_XWXenv(likep->get_weightiwls(),column);
  prec_env.addto(XX_env,Kenv,1.0,lambda);
  lambda_prec = lambda;

  likep->compute_workingresiduals(column);
  compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);

  prec_env.solve(muy,beta);
// monotone Regression!
  if(decreasing)
    {
    bool ok = false;
    while(!ok)
      {
      bool ok2 = true;
      for(unsigned i=1;i<nrpar;i++)
        {
        double diff = beta(i,0)-beta(i-1,0);
        if(diff > 0.0001)
          {
          ok2 = false;
          double mean = 0.5*( beta(i-1,0)+beta(i,0) );
          beta(i-1,0) = mean;
          beta(i,0) = mean;
          }
        if(diff > 0)
          {
          double help = beta(i,0);
          beta(i,0) = beta(i-1,0);
          beta(i-1,0) = help;
          }
        }
      ok = ok2;
      }
    beta.sortcol(0,nrpar-1,0);
    datamatrix bsort = beta;
    for(unsigned j=0;j<nrpar;j++)
      beta(j,0) = bsort(nrpar-1-j,0);
    }

  if(increasing)
    {
    bool ok = false;
    while(!ok)
      {
      bool ok2 = true;
      for(unsigned i=1;i<nrpar;i++)
        {
        double diff = beta(i-1,0)-beta(i,0);
        if(diff > 0.0001)
          {
          ok2 = false;
          double mean = 0.5*(beta(i-1,0)+beta(i,0));
          beta(i-1,0) = mean;
          beta(i,0) = mean;
          }
        if(diff > 0)
          {
          double help = beta(i,0);
          beta(i,0) = beta(i-1,0);
          beta(i-1,0) = help;
          }
        }
      ok = ok2;
      }
    beta.sortcol(0,nrpar-1,0);
    }
// ENDE: monoton
  add_linearpred_multBS();

  if(center)
    {
    if(samplecentered)
      {
      sample_centered_env(beta);
      }
    else
      {
      compute_intercept();

      if (varcoeff)
        fcconst->posteriormode_fix_varcoeff(intercept,datanames[1]);
      else
        fcconst->posteriormode_intercept(intercept);
      }
    }

  if(interaction == false)
    {

    if(samplecentered)
      {
      write_spline();
      write_derivative();
      }
    else
      {
      double * fchelpbetap = fchelp.getbetapointer();

      if(gridsize < 0)
        {
        if(varcoeff)
          multBS(splinehelp,beta);

        vector<int>::iterator freqwork = freqoutput.begin();
        int * workindex = index.getV();
        for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
          {
          if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
            {
            if(varcoeff)
              *fchelpbetap = splinehelp(i,0) - intercept;
            else
              *fchelpbetap = spline(*workindex,0) - intercept;
            fchelpbetap++;
            }
          }
        }
      else
        {
        multDG(splinehelp,beta);
        for(i=0;i<gridsize;i++,fchelpbetap++)
          *fchelpbetap = splinehelp(i,0) - intercept;
        }

      write_derivative();
      }

    if(derivative)
      fcderivative.posteriormode();

    fchelp.posteriormode();
    return FULLCOND_nonp_basis::posteriormode();
    }  // end: if(interaction == false)
  else
    {
    return true;
    }

  }


void FULLCOND_pspline_gaussian::outresults(void)
  {

  FULLCOND::outresults();
  spline_basis::outresults();

  if(hierarchical)
    {
    ST::string fclineffpathsample = pathresults.substr(0,pathresults.length()-12)+"_lineff_sample.raw";
    ST::string fclineffpath = pathresults.substr(0,pathresults.length()-12)+"_lineff.res";

    ofstream out(fclineffpathsample.strtochar());
    out << "intnr   lineff" << endl;
    for(unsigned i=0;i<lineffsamples.size();i++)
      out << (i+1) << "   " << lineffsamples[i] << endl;
    out.close();

    ofstream out2(fclineffpath.strtochar());
    out2 << "lineff" << endl;
    out2 << lineffsum/lineffsamples.size() << endl;
    out2.close();
    }

// compute contour probabilities

  if(pseudocontourprob)
    {
    for(int i=0;i<=contourprob;i++)
      compute_pseudocontourprob(i);

    optionsp->out("\n");
    optionsp->out("\n");
    }

  if(contourprob >= 0)
    {
    fc_contour.outresults();

    compute_contourprob();
    for(int i=1;i<=contourprob;i++)
      compute_contourprob(i);

    optionsp->out("\n");
    optionsp->out("\n");
    }

  }


void FULLCOND_pspline_gaussian::predict(const datamatrix & newX, datamatrix & linpred)
  {

  unsigned i,j;
  datamatrix betac(beta.rows(),beta.cols());
  datamatrix bspline(1,nrpar,0);
  double * worklin = linpred.getV();

  for(i=0;i<nrpar;i++)
    bspline(0,i) = bspline_rek(degree,i,newX);

  if(varcoeff)
    {
    double help;
    for(i=0;i<optionsp->get_samplesize();i++,worklin++)
      {
      help = 0.0;
      readsample2(betac,i);
      for(j=0;j<nrpar;j++)
        {
        help += betac(j,0) * bspline(0,j);
        }
      *worklin += help * newX(0,1);
      }
    }   // end: if
  else
    {
    for(i=0;i<optionsp->get_samplesize();i++,worklin++)
      {
      readsample2(betac,i);
      for(j=0;j<nrpar;j++)
        {
        *worklin += betac(j,0) * bspline(0,j);
        }
      }
    }   // end: else

  }


double FULLCOND_pspline_gaussian::compute_quadform(void)
  {
  if(hierarchical)
    {
    datamatrix b = beta;
    for(unsigned i=0;i<nrpar;i++)
      b(i,0) = beta(i,0) - lineff*gamma(i,0);
    return Kenv.compute_quadform(b,0);
    }
  else if(predictright || predictleft)
    {
    return Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1);
    }
  else
    return Kenv.compute_quadform(beta,0);
  }


void FULLCOND_pspline_gaussian::outoptions(void)
  {
  optionsp->out("  OPTIONS FOR P-SPLINE TERM: " + title + "\n",true);
  spline_basis::outoptions();
  }


void FULLCOND_pspline_gaussian::write_contour(void)
  {

  unsigned i,j;

  if( optionsp->get_nriter() == 1)
    {
    ofstream outm("c:\\cprog\\m.raw");
    ofstream out5("c:\\cprog\\mP.raw");
    ofstream out2("c:\\cprog\\quadform.raw");
    ofstream outsample("c:\\cprog\\sample.raw");
    outm.close();
    out5.close();
    out2.close();
    outsample.close();

    datamatrix beta0(nrpar,1,0);
    ofstream out3("c:\\cprog\\beta0.raw");
    out3 << Kenv.compute_quadform(beta0,0) << "   ";
    out3 << XX_env.compute_quadform(beta0,0) << "   ";
    out3.close();
    }

  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
    {

    datamatrix mP = datamatrix(1,nrpar,0);
    for(i=0;i<nrpar;i++)
      for(j=0;j<nrpar;j++)
        mP(0,i) += betahelp(j,0)*prec_env(j,i);

    ofstream out5("c:\\cprog\\mP.raw",ios::app);
    mP.prettyPrint(out5);
    out5.close();

    ofstream outm("c:\\cprog\\m.raw",ios::app);
    betahelp.transposed().prettyPrint(outm);
    outm.close();

    ofstream outsample("c:\\cprog\\sample.raw",ios::app);
    beta.transposed().prettyPrint(outsample);
    outsample.close();

    ofstream out2("c:\\cprog\\quadform.raw",ios::app);
    out2 << 1.0/sigma2 << "   " << Kenv.compute_quadform(beta,0) << "   ";
    out2 << 1.0/likep->get_scale(column) << "   " << XX_env.compute_quadform(beta,0) << "   ";
    out2 << prec_env.compute_quadform(betahelp,0) << "   ";
    out2 << prec_env.getLogDet() << endl;
    out2.close();

    }

  }


void FULLCOND_pspline_gaussian::compute_contourprob(void)
  {

//  double LIMIT = log(MAXDOUBLE);
  double LIMIT = 700;

  unsigned i,j,k;

  datamatrix mPhelp;

  datamatrix contour(optionsp->get_samplesize(),fc_contour.get_nrpar(),0);
  fc_contour.readsample3(contour);

  mPhelp = datamatrix(nrpar,1,0);

  datamatrix beta_0(nrpar,1,0);

  for(i=0;i<nrpar;i++)
    beta_0(i,0) = fc_contour.get_data(i,0);

  double exponent,mPbeta;

  datamatrix pbeta_0(1,3,0);
  datamatrix pbeta_j(optionsp->get_samplesize(),3,0);

  unsigned step = optionsp->get_samplesize()/1000;
  datamatrix RB(optionsp->get_samplesize()/step,1,0);

  datamatrix beta(optionsp->get_samplesize(),nrpar,0);
  readsample3(beta);
  beta = beta.transposed();

  datamatrix m(optionsp->get_samplesize(),nrpar,0);
  m = contour.getColBlock(0,nrpar);
  m = m.transposed();

  for(j=0;j<optionsp->get_samplesize();j++)
    {
// compute p(beta_j)
    for(i=0;i<optionsp->get_samplesize();i+=step)
      {
      mPbeta = 0.0;
      for(k=0;k<nrpar;k++)
        mPbeta += contour(i,nrpar+6+k)*beta(k,j);

      exponent = contour(i,nrpar)   * contour(j,nrpar+2)
               + contour(i,nrpar+1) * contour(j,nrpar+3)
               + contour(i,nrpar+4) - 2*mPbeta/transform;
      RB(i/step,0) = 0.5*contour(i,nrpar+5) - 0.5*(exponent);
      }

    pbeta_j(j,0) = RB.quantile(50,0);
    pbeta_j(j,1) = RB.mean(0);

    double help = 0.0;
    for(i=0;i<RB.rows();i++)
      help += RB(i,0)<LIMIT?exp( RB(i,0) ):exp(LIMIT);
    pbeta_j(j,2) = help/RB.rows();

    }

// compute p(beta_0)
  for(i=0;i<optionsp->get_samplesize();i+=step)
    {
    prec_env.addto(XX_env,Kenv,contour(i,nrpar),contour(i,nrpar+1));

    mPbeta = 0.0;
    for(k=0;k<nrpar;k++)
      mPbeta += contour(i,nrpar+6+k)*beta_0(k,0);

    exponent = prec_env.compute_quadform(beta_0,0)/(transform*transform)
             + contour(i,nrpar+4) - 2*mPbeta/transform;
    RB(i/step,0) = 0.5*contour(i,nrpar+5) - 0.5*(exponent);
    }

  pbeta_0(0,0) = RB.quantile(50,0);
  pbeta_0(0,1) = RB.mean(0);

  double help = 0.0;
  for(i=0;i<RB.rows();i++)
    help += RB(i,0)<LIMIT?exp( RB(i,0) ):exp(LIMIT);
  pbeta_0(0,2) = help/RB.rows();

//------------------------------------------------------------------------------

  datamatrix contourprob(1,3,0.0);
  for(i=0;i<contourprob.cols();i++)
    {
    for(j=0;j<optionsp->get_samplesize();j++)
      if(pbeta_j(j,i) < pbeta_0(0,i))
        contourprob(0,i)++;
    contourprob(0,i) = contourprob(0,i)/optionsp->get_samplesize();
    }

  optionsp->out("  Contour probability                                : " + ST::doubletostring(contourprob(0,0),2) + "\n");

  ST::string path = pathresult.substr(0,pathresult.length()-4)+"_contour.res";

  ofstream out(path.strtochar());
  out << "difforder   contourprob   mean(log)  mean" << endl;
  out << (ST::inttostring(0) + "   " + ST::doubletostring(contourprob(0,0)) + "   ");
  out << (ST::doubletostring(contourprob(0,1)) + "   " + ST::doubletostring(contourprob(0,2))) << endl;
  out.close();

  }


void FULLCOND_pspline_gaussian::compute_contourprob(const int & diff)
  {

  unsigned i,j,k;

  datamatrix help,mPhelp;

  datamatrix contour(optionsp->get_samplesize(),fc_contour.get_nrpar(),0);
  fc_contour.readsample3(contour);

  datamatrix A = diffmat_k(diff,nrpar);

  unsigned rankA = A.rows();

  mPhelp = datamatrix(rankA,1,0);

  datamatrix beta_0(rankA,1,0.0);
  for(i=0;i<rankA;i++)
    beta_0(i,0) = fc_contour.get_data(i,0);

  envmatdouble tildeP_env;
  datamatrix tildeP;

  double exponent,mPbeta;
  double value;

  unsigned step;                             // jedes wievielte Sample soll verwendet werden?
  if(approx)
    step = optionsp->get_samplesize()/1000;
  else
    step = optionsp->get_samplesize()/1000;

// beta samples einlesen und Differenzen bilden
  datamatrix beta(optionsp->get_samplesize(),nrpar,0);
  readsample3(beta);
  beta = beta.transposed();
  beta = A*beta;

  datamatrix m(optionsp->get_samplesize(),nrpar,0);
  m = contour.getColBlock(0,nrpar);
  m = m.transposed();
  m = A*m;

  double I,Z;
  int n;
  ST::string RBpath;
  ofstream RBfile;
  datamatrix start,dn,d0,fxi,var;       // für approx: pbeta_0, pbeta_j =^ xi
  datamatrix pbeta_0,pbeta_j,RB;

  pbeta_0 = datamatrix(1,3,0);
  pbeta_j = datamatrix(optionsp->get_samplesize(),3,0);

  if(approx)
    {
    start = datamatrix(lengthstart,optionsp->get_samplesize()+1,0);
    d0 = datamatrix(optionsp->get_samplesize()+1,1,0.0);
    dn = datamatrix(optionsp->get_samplesize()+1,1,0.0);
    fxi = datamatrix(optionsp->get_samplesize()+1,1,0.0);
    var = datamatrix(optionsp->get_samplesize()+1,1,0.0);
    }
  else
    {
    RB = datamatrix(optionsp->get_samplesize()/step,1,0);
    RBpath = pathresult.substr(0,pathresult.length()-4) + "_RBfile.raw";
    RBfile.open(RBpath.strtochar(),ios::binary);
    }

  for(i=0;i<optionsp->get_samplesize();i+=step)
    {
// Grenze für 1/sigma2 damit Cholesky-Zerlegung von prec_env funktioniert
    double lim;
    if( contour(i,nrpar+1) > 100000)
      lim = 100000;
    else
      lim = contour(i,nrpar+1);
    prec_env.addto(XX_env,Kenv,contour(i,nrpar),lim);

    if(diff>0)
      {
      tildeP = datamatrix(nrpar,rankA,0);
      for(j=0;j<rankA;j++)                           //   P^-1A^T
        {
        help = A.transposed().getCol(j);
        prec_env.solve(help);
        for(k=0;k<nrpar;k++)
          tildeP(k,j) = help(k,0);
        }
      tildeP = A*tildeP;                              //  AP^-1A^T
      tildeP = tildeP.inverse();                      // (AP^-1A^T)^-1
//      tildeP_env = envmatdouble(tildeP,0.000001);
      tildeP_env = envmatdouble(tildeP,0.01);
      }
    else
      {
      tildeP_env = prec_env;
      }

    mPhelp = datamatrix(rankA,1,0);
    for(j=0;j<rankA;j++)
      for(k=0;k<rankA;k++)
        mPhelp(j,0) += m(k,i)*tildeP_env(k,j);

// compute p(beta_j)
    for(j=0;j<optionsp->get_samplesize();j++)
      {
      mPbeta = 0.0;
      for(k=0;k<rankA;k++)
        mPbeta += mPhelp(k,0)*beta(k,j);

      exponent = tildeP_env.compute_quadform(beta,j)/(transform*transform)
               + tildeP_env.compute_quadform(m,i) - 2* mPbeta/transform;
      value = 0.5*tildeP_env.getLogDet() - 0.5*(exponent);

      if(i==0)
        {
        pbeta_j(j,1) = value;
        pbeta_j(j,2) = exp(value);
        }
      else
        {
        pbeta_j(j,1) = step*( (i-1)/double(step)*pbeta_j(j,1) + value)/double(i);
        pbeta_j(j,2) = step*( (i-1)/double(step)*pbeta_j(j,2) + exp(value))/double(i);
        }

      if(!approx)
        {
        RBfile.write((char *) &value,sizeof(double));
        }
      else  // Tierney
        {
        if(i<lengthstart*step)
          {
          start(i/step,j) = value;

          if(i==(lengthstart-1)*step)
            {
            pbeta_j(j,0) = start.quantile(50,j);
            d0(j,0) = 1.0/(start.quantile(75,j)-start.quantile(25,j));
            fxi(j,0) = 0.0;
            var(j,0) = start.var(j);
            }

          }
        else
          {
          n = i/step-lengthstart;

          double help = fabs(value-pbeta_j(j,0));
          help<=var(j,0)*pow(n+1,-0.5)?I=1:I=0;
          if(n==0)
            dn(j,0) = d0(j,0);
          else
            {
            double help2 = 1.0/(n+1)*(n*fxi(j,0) + 0.5*I/(var(j,0)*pow(n+1,-0.5)));
            fxi(j,0) = help2;
            if(fxi(j,0)==0.0)
              dn(j,0) = d0(j,0)*pow(n,0.5);
            else
              dn(j,0) = min(1.0/fxi(j,0),d0(j,0)*pow(n,0.5));
            }
          value<=pbeta_j(j,0)?Z=1:Z=0;
          pbeta_j(j,0) = pbeta_j(j,0) - dn(j,0)/(n+1)*(Z-0.5);
          }
        }

      } // END:     for(j=0;j<optionsp->get_samplesize();j++)

// compute p(beta_0)
    mPbeta = 0.0;
    for(k=0;k<rankA;k++)
      mPbeta += mPhelp(k,0)*beta_0(k,0);

    exponent = tildeP_env.compute_quadform(beta_0,0)/(transform*transform)
             + tildeP_env.compute_quadform(m,i) - 2* mPbeta/transform;
    value = 0.5*tildeP_env.getLogDet() - 0.5*(exponent);

    if(i==0)
      {
      pbeta_0(0,1) = value;
      pbeta_0(0,2) = exp(value);
      }
    else
      {
      pbeta_0(0,1) = step*( (i-1)/double(step)*pbeta_0(0,1) + value)/double(i);
      pbeta_0(0,2) = step*( (i-1)/double(step)*pbeta_0(0,2) + exp(value))/double(i);
      }

    if(!approx)
      {
      RBfile.write((char *) &value,sizeof(double));
      }
    else    // Tierney
      {
      if(i<lengthstart*step)
        {
        start(i/step,j) = value;

        if(i==(lengthstart-1)*step)
          {
          pbeta_0(0,0) = start.quantile(50,j);
          d0(j,0) = 1.0/(start.quantile(75,j)-start.quantile(25,j));
          fxi(j,0) = 0.0;
          var(j,0) = start.var(j);
          }

        }
      else
        {
        n = i/step-lengthstart;

        double help = fabs(value-pbeta_0(0,0));
        help<=var(j,0)*pow(n+1,-0.5)?I=1:I=0;
        if(n==0)
          dn(j,0) = d0(j,0);
        else
          {
          double help2 = 1.0/(n+1)*(n*fxi(j,0) + 0.5*I/(var(j,0)*pow(n+1,-0.5)));
          fxi(j,0) = help2;
          if(fxi(j,0)==0.0)
            dn(j,0) = d0(j,0)*pow(n,0.5);
          else
            dn(j,0) = min(1.0/fxi(j,0),d0(j,0)*pow(n,0.5));
          }
        value<=pbeta_0(0,0)?Z=1:Z=0;
        pbeta_0(0,0) = pbeta_0(0,0) - dn(j,0)/(n+1)*(Z-0.5);
        }
      }

    }   // END:   for(i=0;i<optionsp->get_samplesize();i+=step)

  if(!approx)
    RBfile.close();

//------------------------------------------------------------------------------

  if(!approx)
    {
    unsigned size = sizeof(double);
    ifstream in;

    double* work;

    for(j=0;j<optionsp->get_samplesize();j++)
      {
      in.open(RBpath.strtochar(),ios::binary);
      work = RB.getV();
      in.seekg(size*j);
      for (i=0;i<optionsp->get_samplesize()/step;i++,work++)
        {
        in.read((char*) work,size);
        in.seekg(size*(optionsp->get_samplesize()),ios::cur);
        }
      in.close();
      pbeta_j(j,0) = RB.quantile(50,0);
      }

    in.open(RBpath.strtochar(),ios::binary);
    work = RB.getV();
    in.seekg(size*j);
    for (i=0;i<optionsp->get_samplesize()/step;i++,work++)
      {
      in.read((char*) work,size);
      in.seekg(size*(optionsp->get_samplesize()),ios::cur);
      }
    pbeta_0(0,0) = RB.quantile(50,0);

    in.close();
    remove(RBpath.strtochar());
    }

// p-Werte berechnen

  datamatrix contourprob(1,3,0.0);

  for(i=0;i<contourprob.cols();i++)
    {
    for(j=0;j<optionsp->get_samplesize();j++)
      if(pbeta_j(j,i) < pbeta_0(0,i))
        contourprob(0,i)++;
    contourprob(0,i) = contourprob(0,i)/optionsp->get_samplesize();
    }

// Ausgabe

  if(diff==0)
    optionsp->out("  Contour probability                                : " + ST::doubletostring(contourprob(0,0),2) + "\n");
  else if(diff==1)
    optionsp->out("  Contour probability for first differences (const)  : " + ST::doubletostring(contourprob(0,0),2) + "\n");
  else if(diff==2)
    optionsp->out("  Contour probability for second differences (linear): " + ST::doubletostring(contourprob(0,0),2) + "\n");
  else if(diff==3)
    optionsp->out("  Contour probability for " + ST::inttostring(diff) + ". differences (quadratic) : " + ST::doubletostring(contourprob(0,0),2) + "\n");
  else if(diff==4)
    optionsp->out("  Contour probability for " + ST::inttostring(diff) + ". differences (cubic)     : " + ST::doubletostring(contourprob(0,0),2) + "\n");
  else
    optionsp->out("  Contour probability for " + ST::inttostring(diff) + ". differences             : " + ST::doubletostring(contourprob(0,0),2) + "\n");

  ST::string path = pathresult.substr(0,pathresult.length()-4)+"_contour.res";

  if(diff==0)
    {
    ofstream out(path.strtochar());
    out << "difforder   contourprob   mean(log)   mean" << endl;
    out.close();
    }

  ofstream out(path.strtochar(),ios::app);
  out << (ST::inttostring(diff) + "   " + ST::doubletostring(contourprob(0,0)) + "   ");
  out << (ST::doubletostring(contourprob(0,1)) + "   " + ST::doubletostring(contourprob(0,2))) << endl;
  out.close();

  }


void FULLCOND_pspline_gaussian::compute_pseudocontourprob(const int & diff)
  {
  unsigned i,j;

// Differenzenmatrix erzeugen
  datamatrix A = diffmat_k(diff,nrpar);

  unsigned rankA = A.rows();
// beta_0 setzen
  datamatrix beta_0(rankA,1,0);
  for(i=0;i<rankA;i++)
    beta_0(i,0) = fc_contour.get_data(i,0);
// beta samples einlesen und Differenzen bilden
  datamatrix beta(optionsp->get_samplesize(),nrpar,0);
  readsample3(beta);
  beta = beta*(A.transposed());

  datamatrix betasort = beta;
// sortieren
  for(i=0;i<rankA;i++)
    betasort.sortcol(0,optionsp->get_samplesize()-1,i);

  double p;
  bool inside;
  unsigned count,t,tstar;

// Intervall berechnen in dem beta_0 gerade noch drin liegt

  tstar = 0;
  while( tstar<optionsp->get_samplesize()-1 &&
         ( betasort(tstar,0) <= beta_0(0,0) && beta_0(0,0) <= betasort(optionsp->get_samplesize()-tstar-1,0) )
       )
    tstar++;

  for(j=1;j<rankA;j++)
    {
    while( 0<tstar &&
           !( betasort(tstar,j) <= beta_0(j,0) && beta_0(j,0) <= betasort(optionsp->get_samplesize()-tstar-1,j) )
         )
      tstar--;
    }

// Wieviele liegen drin?

  count = 0;
  for(t=0;t<optionsp->get_samplesize();t++)
    {
    inside = true;
    for(j=0;j<rankA;j++)
      {
      if( !( betasort(tstar,j) <= beta(t,j) && beta(t,j) <= betasort(optionsp->get_samplesize()-tstar-1,j) ) )
        inside = false;
      }
    if(inside)           // liegt beta^t drin?
      count++;
    }

// p-Wert ausrechnen

  p = 1.0 - (double)count/optionsp->get_samplesize();

// Ausgabe

  if(diff==0)
    optionsp->out("  Pseudo contour probability                                : " + ST::doubletostring(p,2) + "\n");
  else if(diff==1)
    optionsp->out("  Pseudo contour probability for first differences (const)  : " + ST::doubletostring(p,2) + "\n");
  else if(diff==2)
    optionsp->out("  Pseudo contour probability for second differences (linear): " + ST::doubletostring(p,2) + "\n");
  else if(diff==3)
    optionsp->out("  Pseudo contour probability for " + ST::inttostring(diff) + ". differences (quadratic) : " + ST::doubletostring(p,2) + "\n");
  else if(diff==4)
    optionsp->out("  Pseudo contour probability for " + ST::inttostring(diff) + ". differences (cubic)     : " + ST::doubletostring(p,2) + "\n");
  else
    optionsp->out("  Pseudo contour probability for " + ST::inttostring(diff) + ". differences             : " + ST::doubletostring(p,2) + "\n");

  ST::string path = pathresult.substr(0,pathresult.length()-4)+"_pseudocontour.res";

  if(diff==0)
    {
    ofstream out(path.strtochar());
    out << "difforder   pseudocontourprob" << endl;
    out << (ST::inttostring(diff) + "   " + ST::doubletostring(p)) << endl;
    out.close();
    }
  else
    {
    ofstream out(path.strtochar(),ios::app);
    out << (ST::inttostring(diff) + "   " + ST::doubletostring(p)) << endl;
    out.close();
    }

  }



} // end: namespace MCMC



