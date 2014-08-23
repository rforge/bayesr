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
#endif
#pragma hdrstop

#include "IWLS_pspline.h"

using std::ifstream;

namespace MCMC
{


//------------------------------------------------------------------------------
//---------- CLASS: IWLS_pspline (implementation of member functions) ------
//------------------------------------------------------------------------------


void IWLS_pspline::create_iwls(void)
  {

  unsigned i;

  compute_Kweights();

  if(predictright || predictleft)                          // betaweight erweitern, damit in die Zentrierungskonstante
    {                                                      // nur die zu schätzenden Parameter einfliessen, nicht die
    datamatrix help = betaweight;                          // zu prognostizierenden
    betaweight = datamatrix(nrpar,1,0);
    for(i=0;i<nrparpredictleft;i++)
      betaweight(i,0) = 0.0;
    for(;i<nrpar-nrparpredictright;i++)
      betaweight(i,0) = help(i-nrparpredictleft,0);
    for(;i<nrpar;i++)
      betaweight(i,0) = 0.0;
    }

  W = datamatrix(likep->get_nrobs(),1,0);
  mu = datamatrix(likep->get_nrobs(),1,0);

  XX_env = envmatdouble(0.0,nrpar,degree);

  if (type==RW1)
    prec_env = envmatdouble(0.0,nrpar,degree>1?degree:1);
  else if (type == RW2)
    prec_env = envmatdouble(0.0,nrpar,degree>2?degree:2);

  if (type == RW1)                                         // Penalty Matrix erstellen
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

  if(predictleft || predictright)                          // K-matrix anpassen für Prädiktion
    change_K();

  betaold = datamatrix(nrpar,1,0);
  betaprop = datamatrix(nrpar,1,0);
  betahelp = datamatrix(nrpar,1,0);
  standnormal = datamatrix(nrpar,1,0);
  muy = datamatrix(nrpar,1,0);

  }

  // CONSTRUCTOR 1

IWLS_pspline::IWLS_pspline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const fieldtype & ft,const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const double & lk,
                    const double & uk, const double & lg, const double & ug,
                    const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,degr,kp,gs,fp,pres,deriv,lk,uk,lg,ug,c)
  {

  a_invgamma = a;
  b_invgamma = b;

  varcoeff = false;
  identifiable = false;

  diagtransform = diag;

  if(monotone == "increasing")
    increasing = true;
  else if(monotone == "decreasing")
    decreasing = true;

  if(mode)
    {
    if(updatetau)
      utype = hyperblockmode;
    else
      utype = iwlsmode;
    }
  else
    {
    if(updatetau)
      utype = hyperblock;
    else
      utype = iwls;
    }

  updateW = upW;

  f = fstart;

  lambda = l;
  sigma2 = 1.0/l;
  kappa = l;
  kappamode = l;

  compute_betaweight();

  make_index(d);
  make_index2();
  make_Bspline(d);

  create_iwls();
  init_fchelp(d);

  }

  // CONSTRUCTOR 2 (für variierende Koeffizienten)

IWLS_pspline::IWLS_pspline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & effmod,const datamatrix & intact,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l, const fieldtype & ft, const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const double & a, const double & b, const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,degr,kp,gs,fp,pres,deriv,0.0,0.0,0.0,0.0,c)
  {

  diagtransform = false;

  a_invgamma = a;
  b_invgamma = b;

  varcoeff = true;
  identifiable = true;

  diagtransform = diag;

  if(monotone == "increasing")
    increasing = true;
  else if(monotone == "decreasing")
    decreasing = true;

  if(mode)
    {
    if(updatetau)
      utype = hyperblockmode;
    else
      utype = iwlsmode;
    }
  else
    {
    if(updatetau)
      utype = hyperblock;
    else
      utype = iwls;
    }

  updateW = upW;

  f = fstart;

  lambda = l;
  sigma2 = 1.0/l;
  kappa = l;
  kappamode = l;

  compute_betaweight();

  make_index(effmod,intact);
  make_index2();
  make_Bspline(effmod);
  make_BS(intact);

  create_iwls();
  init_fchelp(effmod);

  }

  // CONSTRUCTOR 3 (für Cox)

IWLS_pspline::IWLS_pspline(MCMCoptions * o,DISTRIBUTION * dp,FULLCOND_const * fcc,
                    const datamatrix & d,const bool & mode,
                    const unsigned & nrk,const unsigned & degr,const knotpos & kp,
                    const double & l,const fieldtype & ft,const ST::string & monotone,
                    const unsigned & upW, const bool & updatetau, const double & fstart,
                    const ST::string & ti,
                    const ST::string & fp, const ST::string & pres, const bool & deriv,
                    const int & gs, const bool & diag, const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,degr,kp,gs,fp,pres,deriv,0.0,0.0,0.0,0.0,c)
  {

  varcoeff = false;
  identifiable = false;

  diagtransform = diag;

  if(monotone == "increasing")
    increasing = true;
  else if(monotone == "decreasing")
    decreasing = true;

  if(mode)
    {
    if(updatetau)
      utype = hyperblockmode;
    else
      utype = iwlsmode;
    }
  else
    {
    if(updatetau)
      utype = hyperblock;
    else
      utype = iwls;
    }

  updateW = upW;

  f = fstart;

  lambda = l;
  sigma2 = 1.0/l;
  kappa = l;
  kappamode = l;

  compute_betaweight();

  make_index(d);
  make_index2();
  make_Bspline(d,true);                               // Unterschied zu Construktor 1 und 2:
                                                      // Intervall beginnt bei 0 ( [0,x_max] )
  create_iwls();
  init_fchelp(d);

  }

  // COPY CONSTRUCTOR

IWLS_pspline::IWLS_pspline(const IWLS_pspline & fc)
  :spline_basis(spline_basis(fc))
  {

  utype = fc.utype;

  a_invgamma = fc.a_invgamma;
  b_invgamma = fc.b_invgamma;

  kappa = fc.kappa;
  kappamode = fc.kappamode;
  kappaprop = fc.kappaprop;
  kappamean = fc.kappamean;

  diagtransform = fc.diagtransform;
  updateW = fc.updateW;

  }

  // OVERLOADED ASSIGNMENT OPERATOR

const IWLS_pspline & IWLS_pspline::operator=(const IWLS_pspline & fc)
  {
  if (this == &fc)
    return *this;
  spline_basis::operator=(spline_basis(fc));

  utype = fc.utype;

  a_invgamma = fc.a_invgamma;
  b_invgamma = fc.b_invgamma;

  kappa = fc.kappa;
  kappamode = fc.kappamode;
  kappaprop = fc.kappaprop;
  kappamean = fc.kappamean;

  diagtransform = fc.diagtransform;
  updateW = fc.updateW;

  return *this;
  }


void IWLS_pspline::outoptions(void)
  {

  optionsp->out("  OPTIONS FOR P-SPLINE TERM: " + title + "\n",true);

  spline_basis::outoptions();

  if(utype == iwlsmode || utype == hyperblockmode)
    optionsp->out("  Proposal: IWLS based on posterior mode estimation\n");
  else
    optionsp->out("  Proposal: IWLS\n");

  if(updateW == 0)
    optionsp->out("  Weight matrix W is fixed for the whole simulation\n");
  else if(updateW == 1)
    optionsp->out("  Weight matrix W is updated in every iteration\n");
  else if(updateW == 2)
    optionsp->out("  Weight matrix W is updated in every 2nd iteration\n");
  else if(updateW == 3)
    optionsp->out("  Weight matrix W is updated in every 3rd iteration\n");
  else
    optionsp->out("  Weight matrix W is updated in every " + ST::inttostring(updateW) + "th iteration\n");

  optionsp->out("\n");

  if(utype == hyperblock || utype == hyperblockmode)
    {
    optionsp->out("  Updating scheme: single block (including variance parameter)\n");
    optionsp->out("  Starting value for tuning parameter f: " + ST::doubletostring(f) + "\n");
    }
  else
    {
    optionsp->out("  Updating scheme: single block\n");
    }

  optionsp->out("\n");

  }


bool IWLS_pspline::posteriormode_converged(const unsigned & itnr)
  {
//  if(increasing || decreasing || diagtransform)
  if(diagtransform)
    return true;
  else
    return likep->posteriormode_converged_fc(beta,beta_mode,itnr);
  }


bool IWLS_pspline::posteriormode(void)
  {

//  if(increasing || decreasing || diagtransform)
  if(diagtransform)
    {
    return true;
    }
  else
    {

    compute_XWXenv(likep->get_weightiwls(),column);

    prec_env.addto(XX_env,Kenv,1.0,lambda);
    lambda_prec = lambda;

    likep->substr_linearpred_m(spline,column,true);
    likep->compute_workingresiduals(column);

    compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);

    prec_env.solve(muy,beta);

    if(decreasing)                                          // monotone Regression
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
      }                                                     // ENDE: monotone Regression

    add_linearpred_multBS(beta,true);

    if(center)                                              // zentrieren
      {
      compute_intercept();
      for(unsigned i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
      fcconst->posteriormode_intercept(intercept);
      for(unsigned i=0;i<likep->get_nrobs();i++)
        spline(i,0) -= intercept;
      intercept = 0.0;
      }

    if(interaction == false)                                // wird bei interaction==true in der full conditional
      {                                                     // des Interaktionseffekts gemacht
      write_spline();
      write_derivative();

      if(derivative)
        fcderivative.posteriormode();

      fchelp.posteriormode();
      return FULLCOND::posteriormode();
      }  // end: if(interaction == false)
    else
      {
      return true;
      }
    }

  }


void IWLS_pspline::update(void)
  {

  if(lambdaconst == true)
    sigma2 = likep->get_scale(column)/lambda;

  if(optionsp->get_nriter()==1)                             // posterior mode Schätzung übernehmen
    betaold.assign(beta);

  if(increasing || decreasing)                              // update bei Monotonie-Restriktion
    {
    update_isotonic();
    }
  else if(diagtransform)                                    // update bei Diagonal-Transformation
    {
    update_diagtransform();
    }
  else if(utype == iwls)                                    // update nach IWLS
    {
    update_IWLS();
    }
  else if(utype == iwlsmode)                                // update nach IWLS basierend auf posterior mode
    {
    update_IWLS_mode();
    }
  else if(utype == hyperblock)                              // kompletten Block inklusive Hyperparameter updaten
    {
    update_IWLS_hyperblock();
    }
  else if(utype == hyperblockmode)                          // kompletten Block inklusive Hyperparameter updaten (basierend auf posterior mode)
    {
    update_IWLS_hyperblock_mode();
    }


  if(predictright || predictleft)                           // Prädiktion
    update_prediction();

  if(interaction == false)                                  // wird bei interaction==true in der full conditional des Interaktionseffekts gemacht
    {                                                       // spline in fchelp schreiben, etc.
    if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
        ((optionsp->get_nriter()-optionsp->get_burnin()-1) % optionsp->get_step() == 0) )
      {
      if(diagtransform)
        write_spline(G*beta);
      else
        write_spline();

      write_derivative();
      }

    if(derivative)
      fcderivative.update();

    fchelp.update();
    FULLCOND::update();

    }  // end: if(interaction == false)

  } // end: function update


void IWLS_pspline::update_IWLS_hyperblock_mode(void)
  {

  unsigned i;

  double aprop,bprop;

// Akzeptanzwahrscheinlichkeit in der burnin Phase auf 60% tunen
  if(optionsp->get_nriter()%100==0 && optionsp->get_nriter()<optionsp->get_burnin()/2)
    tune_updatetau(alpha_60);

  if(optionsp->get_nriter() == optionsp->get_burnin())
    optionsp->out("NOTE: Tuning constant 'f' for term " + title + " set to " + ST::doubletostring(f) + "\n");

// startwert für kappamode berechnen

  aprop = a_invgamma + 0.5*rankK;
  bprop = b_invgamma + 0.5*Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1);

  if(optionsp->get_nriter() == optionsp->get_burnin()/2)
    {
    kappamean = 0.5*Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1);
    }
  else if(optionsp->get_nriter() > optionsp->get_burnin()/2 && optionsp->get_nriter() < optionsp->get_burnin() )
    {
    kappamean += 0.5*Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1);
    }
  else if(optionsp->get_nriter() == optionsp->get_burnin())
    {
    kappamean = 2*kappamean/optionsp->get_burnin();
    }
  else if(optionsp->get_nriter() > optionsp->get_burnin())
    {
    bprop = b_invgamma + kappamean;
    }

  kappamode = (aprop-1)/bprop;                               // kappamode setzen
  sigma2 = 1.0/kappamode;

  kappaprop = kappamode*randnumbers::rand_variance(f);       // neues kappa vorschlagen

  double logold = likep->loglikelihood(true);
  logold += - 0.5*Kenv.compute_quadformblock(betaold,0,nrparpredictleft,nrpar-nrparpredictright-1)*kappa;
  logold += 0.5*rankK*log(kappa);                            // Normalisierungs Konstante
  logold += (a_invgamma-1)*log(kappa) - b_invgamma*kappa;    // gamma prior

  add_linearpred_multBS(beta_mode,betaold,true);
//  multBS_index(spline,beta_mode);

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
//    likep->compute_weight(W,column,true);
//    likep->tilde_y(mu,spline,column,true,W);
    likep->compute_IWLS_weight_tildey(W,mu,column,true);
    mu.plus(spline,mu);

//    compute_XWXenv(W);
//    compute_XWtildey(W,1.0);
    compute_XWXenv_XWtildey(W,1.0);
    }
  else
    {
    likep->tilde_y(mu,spline,column,true,W);
    compute_XWtildey(W,1.0);
    }

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);                //  fixes sigma2 für beta_mode (m)

  prec_env.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.addto(XX_env,Kenv,1.0,kappaprop);                 //  P(kappaprop)

  prec_env.solveU(beta,betahelp);

  add_linearpred_multBS(beta,beta_mode,true);
  beta_mode.assign(betahelp);
  betahelp.minus(beta,beta_mode);

  double qold = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);
  qold += 0.5*prec_env.getLogDet();
//  qold += (aprop-1)*log(kappaprop) - bprop*kappaprop;        // Gamma proposal
  qold += log(1.0/kappamode + 1.0/kappaprop);                // kappa ~ 1+1/kappa

  double lognew = likep->loglikelihood(true);
  lognew += - 0.5*Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1)*kappaprop;
  lognew += 0.5*rankK*log(kappaprop);
  lognew += (a_invgamma-1)*log(kappaprop) - b_invgamma*kappaprop;

  prec_env.addto(XX_env,Kenv,1.0,kappa);                     //  P(kappa)

  betahelp.minus(betaold,beta_mode);
  double qnew = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);
  qnew += 0.5*prec_env.getLogDet();
//  qnew += (aprop-1)*log(kappa) - bprop*kappa;                // Gamma proposal
  qnew += log(1.0/kappamode + 1.0/kappa);                    // kappa ~ 1+1/kappa

  double alpha = lognew + qnew - logold - qold;
  double u = log(uniform());

  if(center)
    {
    compute_intercept();
    for(i=0;i<nrpar;i++)
      beta_mode(i,0) -= intercept;
    intercept = 0.0;
    }

  if(u <= alpha)
    {
    acceptance++;
    kappa = kappaprop;
    sigma2 = 1.0/kappa;                                      // für variance_nonp
    if(center)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
      fcconst->update_intercept(intercept);
      for(i=0;i<likep->get_nrobs();i++)
        spline(i,0) -= intercept;
      intercept = 0.0;
      }
    betaold.assign(beta);
    }
  else
    {
    add_linearpred_multBS(betaold,beta,true);
    beta.assign(betaold);
    }

  sigma2 = 1.0/kappa;

  }


void IWLS_pspline::update_IWLS_hyperblock(void)
  {

  unsigned i;

  if(optionsp->get_nriter()%100==0 && optionsp->get_nriter()<optionsp->get_burnin())
    tune_updatetau(alpha_60);

  if(optionsp->get_nriter() == optionsp->get_burnin())
    optionsp->out("NOTE: Tuning constant 'f' for term " + title + " set to " + ST::doubletostring(f) + "\n");

  kappaprop = kappa*randnumbers::rand_variance(f);                  // kappa ~ (1+1/z)

  double logold = - 0.5*Kenv.compute_quadformblock(betaold,0,nrparpredictleft,nrpar-nrparpredictright-1)*kappa;
  logold += 0.5*rankK*log(kappa);             // Normalisierungs Konstante
  logold += (a_invgamma-1)*log(kappa) - b_invgamma*kappa;           // gamma prior

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
//    likep->compute_weight(W,column,true);
//    likep->tilde_y(mu,spline,column,true,W);
    logold += likep->compute_IWLS(W,mu,true,column,true);
    mu.plus(spline,mu);

//    compute_XWXenv(W);
//    compute_XWtildey(W,1.0);
    compute_XWXenv_XWtildey(W,1.0);
    }
  else
    {
//    likep->tilde_y(mu,spline,column,true,W);
    logold += likep->compute_IWLS(W,mu,false,column,true);
    mu.plus(spline,mu);

    compute_XWtildey(W,1.0);
    }

  prec_env.addto(XX_env,Kenv,1.0,kappaprop);

  prec_env.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.solveU(beta,betahelp);

  add_linearpred_multBS(beta,betaold,true);

  betahelp.minus(beta,betahelp);
  double qold = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);
  qold += 0.5*prec_env.getLogDet();

  double lognew = - 0.5*Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1)*kappaprop;
  lognew += 0.5*rankK*log(kappaprop);
  lognew += (a_invgamma-1)*log(kappaprop) - b_invgamma*kappaprop;

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
//    likep->compute_weight(W,column,true);
//    likep->tilde_y(mu,spline,column,true,W);
    lognew += likep->compute_IWLS(W,mu,true,column,true);
    mu.plus(spline,mu);
//    compute_XWXenv(W);
//    compute_XWtildey(W,1.0);
    compute_XWXenv_XWtildey(W,1.0);
    }
  else
    {
//    likep->tilde_y(mu,spline,column,true,W);
    lognew += likep->compute_IWLS(W,mu,false,column,true);
    mu.plus(spline,mu);

    compute_XWtildey(W,1.0);
    }

  prec_env.addto(XX_env,Kenv,1.0,kappa);

  prec_env.solve(muy,betahelp);

  betahelp.minus(betaold,betahelp);
  double qnew = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);
  qnew += 0.5*prec_env.getLogDet();

  double alpha = lognew + qnew - logold - qold;
  double u = log(uniform());

  if(u<=alpha)
    {
    kappa = kappaprop;
    sigma2 = 1.0/kappa;                                             // für variance_nonp
    acceptance++;
    if(center)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
      fcconst->update_intercept(intercept);
      for(i=0;i<likep->get_nrobs();i++)
        spline(i,0) -= intercept;
      intercept = 0.0;
      }
    betaold.assign(beta);
    }
  else
    {
    add_linearpred_multBS(betaold,beta,true);
    beta.assign(betaold);
    }

  }


void IWLS_pspline::update_IWLS_mode(void)
  {

  unsigned i;

  double logold = likep->loglikelihood(true)
                - 0.5*Kenv.compute_quadformblock(betaold,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;

  add_linearpred_multBS(beta_mode,betaold,true);
//  multBS_index(spline,beta_mode);

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
//    likep->compute_weight(W,column,true);
//    likep->tilde_y(mu,spline,column,true,W);
    likep->compute_IWLS_weight_tildey(W,mu,column,true);
    mu.plus(spline,mu);
//    compute_XWXenv(W);
//    compute_XWtildey(W,1.0);
    compute_XWXenv_XWtildey(W,1.0);
    }
  else
    {
    likep->tilde_y(mu,spline,column,true,W);
    compute_XWtildey(W,1.0);
    }

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

  prec_env.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.solveU(beta,betahelp);

  add_linearpred_multBS(beta,beta_mode,true);
  beta_mode.assign(betahelp);

  betahelp.minus(beta,beta_mode);
  double qold = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  betahelp.minus(betaold,beta_mode);
  double lognew = likep->loglikelihood(true)
                - 0.5*Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;
  double qnew = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  double alpha = lognew + qnew - logold - qold;
  double u = log(uniform());

  if(center)
    {
    compute_intercept();
    for(i=0;i<nrpar;i++)
      beta_mode(i,0) -= intercept;
    intercept = 0.0;
    }

  if(u<=alpha)
    {
    acceptance++;
    if(center)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
      fcconst->update_intercept(intercept);
      for(i=0;i<likep->get_nrobs();i++)
        spline(i,0) -= intercept;
      intercept = 0.0;
      }
    betaold.assign(beta);
    }
  else
    {
    add_linearpred_multBS(betaold,beta,true);
    beta.assign(betaold);
    }

  }


void IWLS_pspline::update_IWLS(void)
  {

  unsigned i;

  double logold = - 0.5*Kenv.compute_quadformblock(betaold,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
//    likep->compute_weight(W,column,true);
//    likep->tilde_y(mu,spline,column,true,W);
    logold += likep->compute_IWLS(W,mu,true,column,true);
    mu.plus(spline,mu);

//    compute_XWXenv(W);
//    compute_XWtildey(W,1.0);
    compute_XWXenv_XWtildey(W,1.0);
    }
  else
    {
//    likep->tilde_y(mu,spline,column,true,W);
    logold += likep->compute_IWLS(W,mu,false,column,true);
    mu.plus(mu,spline);

    compute_XWtildey(W,1.0);
    }

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

  prec_env.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.solveU(beta,betahelp);

  add_linearpred_multBS(beta,betaold,true);
  betahelp.minus(beta,betahelp);

  double qold = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  double lognew = - 0.5*Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    qold += 0.5*prec_env.getLogDet();

//    likep->compute_weight(W,column,true);
//    likep->tilde_y(mu,spline,column,true,W);
    lognew += likep->compute_IWLS(W,mu,true,column,true);
    mu.plus(spline,mu);

//    compute_XWXenv(W);
//    compute_XWtildey(W,1.0);
    compute_XWXenv_XWtildey(W,1.0);
    prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);
    }
  else
    {
//    likep->tilde_y(mu,spline,column,true,W);
    lognew += likep->compute_IWLS(W,mu,false,column,true);
    mu.plus(mu,spline);

    compute_XWtildey(W,1.0);
    }

  prec_env.solve(muy,betahelp);

  betahelp.minus(betaold,betahelp);
  double qnew = - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    qnew += 0.5*prec_env.getLogDet();
    }

  double alpha = lognew + qnew - logold - qold;
  double u = log(uniform());

  if(u<=alpha)
    {
    acceptance++;
    if(center)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
      fcconst->update_intercept(intercept);
      for(i=0;i<likep->get_nrobs();i++)
        spline(i,0) -= intercept;
      intercept = 0.0;
      }
    betaold.assign(beta);
    }
  else
    {
    add_linearpred_multBS(betaold,beta,true);
    beta.assign(betaold);
    }

  }


void IWLS_pspline::update_isotonic(void)
  {

  unsigned i,j;
  double m,help;
  double unten = -MAXDOUBLE;
  double oben = MAXDOUBLE;

//*------------------------ single move --------------------------------------//
  double logold,lognew,qold,qnew,alpha,u;

  for(i=0;i<nrpar;i++)
    {

    nrtrials++;
/*
  if(optionsp->get_nriter() < optionsp->get_burnin())
    {
    likep->compute_weight(W,column,true);
    compute_XWXenv(W);

    prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

    likep->tilde_y(mu,spline,column,true,W);
    compute_XWtildey(W,1.0);
    }
*/
    beta.assign(betaold);

    help = 0.0;
    for(j=0;j<i;j++)
      help += prec_env(i,j)*beta(j,0);
    for(j=i+1;j<nrpar;j++)
      help += prec_env(i,j)*beta(j,0);
    m = (muy(i,0) - help)/prec_env(i,i);

    beta(i,0) = sample_monotonic(i,m,sqrt(1.0/prec_env(i,i)));

    logold = likep->loglikelihood(true) - 0.5*Kenv.compute_quadform(betaold,0)/sigma2;
    qold = 0.5*log(prec_env(i,i)) - 0.5*(beta(i,0)-m)*prec_env(i,i)*(beta(i,0)-m);

    if(i>0)
      unten = beta(i-1,0);
    if(i<nrpar-1)
      oben = beta(i+1,0);

    double helpold;
    if(i==0)
      helpold = randnumbers::Phi2( (oben-m) * sqrt(prec_env(i,i)) );
    else if (i==nrpar-1)
      helpold = 1.0 - randnumbers::Phi2( (unten-m) * sqrt(prec_env(i,i)) );
    else
      helpold = randnumbers::Phi2( (oben-m) * sqrt(prec_env(i,i)) ) - randnumbers::Phi2( (unten-m) * sqrt(prec_env(i,i)) );

    add_linearpred_multBS(beta,betaold,true);
/*
  if(optionsp->get_nriter() < optionsp->get_burnin())
    {
    likep->compute_weight(W,column,true);
    compute_XWXenv(W);

    prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

    likep->tilde_y(mu,spline,column,true,W);
    compute_XWtildey(W,1.0);
    }
*/
    help = 0.0;
    for(j=0;j<i;j++)
      help += prec_env(i,j)*beta(j,0);
    for(j=i+1;j<nrpar;j++)
      help += prec_env(i,j)*beta(j,0);
    m = (muy(i,0) - help)/prec_env(i,i);

    lognew = likep->loglikelihood(true) - 0.5*Kenv.compute_quadform(beta,0)/sigma2;
    qnew = 0.5*log(prec_env(i,i)) - 0.5*(betaold(i,0)-m)*prec_env(i,i)*(betaold(i,0)-m);

    if(i>0)
      unten = betaold(i-1,0);
    if(i<nrpar-1)
      oben = betaold(i+1,0);

    double helpnew;
    if(i==0)
      helpnew = randnumbers::Phi2( (oben-m) * sqrt(prec_env(i,i)) );
    else if(i==nrpar-1)
      helpnew = 1.0 - randnumbers::Phi2( (unten-m) * sqrt(prec_env(i,i)) );
    else
      helpnew = randnumbers::Phi2( (oben-m) * sqrt(prec_env(i,i)) ) - randnumbers::Phi2( (unten-m) * sqrt(prec_env(i,i)) );

    alpha = lognew + qnew - logold - qold;

    if(helpnew < 1e-100)
      helpnew = 1e-100;
    if(helpold < 1e-100)
      helpold = 1e-100;

    alpha += log( helpold/helpnew );

    u = log(uniform());

    if(u<=alpha)
      {
      acceptance++;
      if(center)
        {
        compute_intercept();
        for(j=0;j<nrpar;j++)
          beta(j,0) -= intercept;
        fcconst->update_intercept(intercept);
        for(j=0;j<likep->get_nrobs();j++)
          spline(j,0) -= intercept;
        intercept = 0.0;
        }
      betaold.assign(beta);
      }
    else
      {
      add_linearpred_multBS(betaold,beta,true);
      beta.assign(betaold);
      }

    }  // END:   for(i=0;i<nrpar;i++)

//------------------------ END: single move ----------------------------------*/

/*---------------------------- Geweke ---------------------------------------//     FALSCH!!! wegen Normierungskonstante!!!!
  likep->compute_weight(W,column,true);
  compute_XWXenv(W);

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

  likep->tilde_y(mu,spline,column,true,W);
  compute_XWtildey(W,1.0);
  prec_env.solve(muy,betahelp);

  int count = 0;
  int maxit = 250;

  beta.assign(betaold);
  count = 0;

  while(count < maxit)
    {

    help = 0.0;
    for(j=1;j<nrpar;j++)
      help += prec_env(0,j)*(beta(j,0)-betahelp(j,0));
    m = betahelp(0,0) - help/prec_env(0,0);
//    if(increasing)
//      beta(0,0) = trunc_normal2(-20,beta(1,0),m,sqrt(1.0/prec_env(0,0)));
//    else
//      beta(0,0) = trunc_normal2(beta(1,0),20,m,sqrt(1.0/prec_env(0,0)));

    beta(0,0) = sample_monotonic(0,m,sqrt(1.0/prec_env(0,0)));

    for(i=1;i<nrpar-1;i++)
      {

      help = 0.0;
      for(j=0;j<i;j++)
        help += prec_env(i,j)*(beta(j,0)-betahelp(j,0));
      for(j=i+1;j<nrpar;j++)
        help += prec_env(i,j)*(beta(j,0)-betahelp(j,0));
      m = betahelp(i,0) - help/prec_env(i,i);
//      if(increasing)
//        beta(i,0) = trunc_normal2(beta(i-1,0),beta(i+1,0),m,sqrt(1.0/prec_env(i,i)));
//      else
//        beta(i,0) = trunc_normal2(beta(i+1,0),beta(i-1,0),m,sqrt(1.0/prec_env(i,i)));

      beta(i,0) = sample_monotonic(i,m,sqrt(1.0/prec_env(i,i)));
      }

    help = 0.0;
    for(j=0;j<nrpar-1;j++)
      help += prec_env(nrpar-1,j)*(beta(j,0)-betahelp(j,0));
    m = betahelp(nrpar-1,0) - help/prec_env(nrpar-1,nrpar-1);
//    if(increasing)
//      beta(nrpar-1,0) = trunc_normal2(beta(nrpar-2,0),20,m,sqrt(1.0/prec_env(nrpar-1,nrpar-1)));
//    else
//      beta(nrpar-1,0) = trunc_normal2(-20,beta(nrpar-2,0),m,sqrt(1.0/prec_env(nrpar-1,nrpar-1)));

    beta(nrpar-1,0) = sample_monotonic(nrpar-1,m,sqrt(1.0/prec_env(nrpar-1,nrpar-1)));

    count++;

    }  // END: while

  betahelp.minus(beta,betahelp);
  double logold = likep->loglikelihood(true)
                - 0.5*Kenv.compute_quadformblock(betaold,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;
  double qold = 0.5*prec_env.getLogDet()
                - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  add_linearpred_multBS(beta,betaold,true);

  likep->compute_weight(W,column,true);
  compute_XWXenv(W);

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

  likep->tilde_y(mu,spline,column,true,W);
  compute_XWtildey(W,1.0);

  prec_env.solve(muy,betahelp);

  betahelp.minus(betaold,betahelp);
  double lognew = likep->loglikelihood(true)
                - 0.5*Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1)/sigma2;
  double qnew = 0.5*prec_env.getLogDet()
                - 0.5*prec_env.compute_quadformblock(betahelp,0,nrparpredictleft,nrpar-nrparpredictright-1);

  double alpha = lognew + qnew - logold - qold;
  double u = log(uniform());

  if(u<=alpha)
    {
    acceptance++;
    if(center)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
      fcconst->update_intercept(intercept);
      for(i=0;i<likep->get_nrobs();i++)
        spline(i,0) -= intercept;
      intercept = 0.0;
      }
    betaold.assign(beta);
    }
  else
    {
    add_linearpred_multBS(betaold,beta,true);
    beta.assign(betaold);
    }
//----------------------- END: Geweke ----------------------------------------*/

/*----------------------- END: Lui & Sabatti ---------------------------------//
//----------------------- END: Lui & Sabatti ---------------------------------*/

  }


void IWLS_pspline::update_diagtransform(void)
  {

  if(optionsp->get_nriter() == 1)
    {

    ifstream inw("c:\\cprog\\testmcmc\\W0.25.raw");
    W.prettyScan(inw);
    inw.close();

/*
    XX = bandmatdouble(nrpar,degree,0);
    compute_XWX(W,column);

    bandmatdouble LL = XX;
    LL.decomp();
    datamatrix L = datamatrix(nrpar,nrpar,0);
    LL.getL(L);

    ofstream outL("c:\\cprog\\testmcmc\\L.raw");
    L.prettyPrint(outL);
    outL.close();
*/

    G = datamatrix(nrpar,nrpar);
    ifstream in("c:\\cprog\\testmcmc\\G.raw");
    G.prettyScan(in);
    in.close();

    betahelp = datamatrix(1,nrpar);
    ifstream in2("c:\\cprog\\testmcmc\\D.raw");
    betahelp.prettyScan(in2);
    in2.close();

    Kenv = envmatdouble(bandmatdouble(betahelp.transposed()));
    rankK = nrpar;
    }

  unsigned i;
  double m;

  double logold,lognew,qold,qnew,alpha,u;
  double prec_ii;

  for(i=0;i<nrpar;i++)
    {

    nrtrials++;

    likep->tilde_y(mu,spline,column,true,W);
    compute_XWtildey(W,1.0);
    muy = G.transposed()*muy;

    beta.assign(betaold);

    prec_ii = 1.0 + betahelp(0,i)/sigma2;

    beta(i,0) = rand_normal()/sqrt(prec_ii);
    m = muy(i,0)/prec_ii;
    beta(i,0) += m;

    logold = likep->loglikelihood(true) - 0.5*Kenv.compute_quadform(betaold,0)/sigma2;
    qold = 0.5*log(prec_ii) - 0.5*(beta(i,0)-m)*prec_ii*(beta(i,0)-m);

    add_linearpred_multBS(G*beta,G*betaold,true);

    likep->tilde_y(mu,spline,column,true,W);
    compute_XWtildey(W,1.0);
    muy = G.transposed()*muy;

    m = muy(i,0)/prec_ii;

    lognew = likep->loglikelihood(true) - 0.5*Kenv.compute_quadform(beta,0)/sigma2;
    qnew = 0.5*log(prec_ii) - 0.5*(betaold(i,0)-m)*prec_ii*(betaold(i,0)-m);

    alpha = lognew + qnew - logold - qold;
    u = log(uniform());

    if(u<=alpha)
      {
      acceptance++;
      betaold.assign(beta);
      }
    else
      {
      add_linearpred_multBS(G*betaold,G*beta,true);
      beta.assign(betaold);
      }

    }

  if(center)
    {

    compute_intercept(G*beta);
    datamatrix help = G.inverse();
    for(i=0;i<nrpar;i++)
      {
      double intercepthelp = 0.0;
      for(unsigned j=0;j<nrpar;j++)
        intercepthelp += help(i,j)*intercept;
      beta(i,0) -= intercepthelp;
      }
    fcconst->update_intercept(intercept);
    for(i=0;i<likep->get_nrobs();i++)
      spline(i,0) -= intercept;
    intercept = 0.0;
    betaold.assign(beta);
    }

  }


void IWLS_pspline::outresults(void)
  {
  FULLCOND::outresults();
  spline_basis::outresults();
  }


void IWLS_pspline::predict(const datamatrix & newX, datamatrix & linpred)
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
    }
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
    }

  }


double IWLS_pspline::compute_quadform(void)
  {
  if(predictright || predictleft)
    {
/*
    double quadform = Kenv.compute_quadform(beta,0);
    if(predictright)
      for(unsigned i=0;i<nrparpredictright;i++)
        quadform -= beta(nrpar-1-i,0)*beta(nrpar-1-i,0);
    if(predictleft)
      for(unsigned i=0;i<nrparpredictleft;i++)
        quadform -= beta(i,0)*beta(i,0);
    return quadform;
*/
    return Kenv.compute_quadformblock(beta,0,nrparpredictleft,nrpar-nrparpredictright-1);
    }
  else
    return Kenv.compute_quadform(beta,0);
  }


} // end: namespace MCMC


//---------------------------------------------------------------------------
#if !defined(__BUILDING_GNU)
#pragma package(smart_init)
#endif


