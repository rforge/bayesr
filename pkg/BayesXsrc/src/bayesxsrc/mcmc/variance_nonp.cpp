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



#include"variance_nonp.h"

namespace MCMC
{

FULLCOND_variance_nonp::FULLCOND_variance_nonp(MCMCoptions * o,
                         FULLCOND_nonp_basis * p,DISTRIBUTION * d,
                         const double & a, const double & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,const bool & av,
                         const unsigned & c)
                         : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
    {

    fctype = MCMC::variance;
    constlambda=false;
    uniformprior=false;
    Laplace=false;
    stationary=false;
    alphafix=false;
    discrete=false;
    update_sigma2 = true;
    randomeffect=false;
    hrandom=false;
    fullcondnonp = false;
    average = av;
    column = c;
    pathresults = fr;
    Kp = p;
    distrp = d;
    rankK= Kp->get_rankK();
    a_invgamma = a;
    b_invgamma = b;
    priorassumptions.push_back("Inverse gamma prior for variance component with hyperparameters a="
    +ST::doubletostring(a,6)+ " and b=" + ST::doubletostring(b,6) );
    priorassumptions.push_back("\\\\");

    if (average==true)
      setbeta(2,1,distrp->get_scale(column,column)/Kp->getlambda());
    else
      setbeta(1,1,distrp->get_scale(column,column)/Kp->getlambda());

    scale = 1.0;
    df = 0;

    ST::string path = samplepath.substr(0,samplepath.length()-4)+"_lambda.raw";
    fc_lambda = FULLCOND(o,datamatrix(1,1),Kp->get_title()+"_lambda",1,1,path);
    fc_lambda.setflags(MCMC::norelchange | MCMC::nooutput);
    }


FULLCOND_variance_nonp::FULLCOND_variance_nonp(MCMCoptions * o,
                         FULLCOND_nonp * p,DISTRIBUTION * d,
                         const double & a, const double & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,const bool & av,
                         const unsigned & c)
                         : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
    {
    fctype = MCMC::variance;
    constlambda=false;
    uniformprior=false;
    Laplace=false;
    stationary=false;
    alphafix=false;
    discrete=false;
    update_sigma2 = true;
    randomeffect=false;
    hrandom=false;
    fullcondnonp = true;
    average = av;
    column = c;
    pathresults = fr;
    Fnp = p;
    distrp = d;
    rankK= Fnp->get_rankK();
    a_invgamma = a;
    b_invgamma = b;
    priorassumptions.push_back("Inverse gamma prior for variance component with hyperparameters a="
    +ST::doubletostring(a,6)+ " and b=" + ST::doubletostring(b,6) );
    priorassumptions.push_back("\\\\");

    if (average==true)
      setbeta(2,1,distrp->get_scale(column,column)/Fnp->getlambda());
    else
      setbeta(1,1,distrp->get_scale(column,column)/Fnp->getlambda());
    }


FULLCOND_variance_nonp::FULLCOND_variance_nonp(MCMCoptions * o,
                         FULLCOND_random * p,DISTRIBUTION * d,
                         const double & a, const double & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,const bool & av,
                         const unsigned & c)
                         : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
    {
    fctype = MCMC::variance;
    constlambda=false;
    uniformprior=false;
    Laplace=false;
    stationary=false;
    alphafix=false;
    discrete=false;
    update_sigma2 = true;
    randomeffect = true;
    hrandom=false;
    fullcondnonp = false;
    average = av;
    column = c;
    pathresults = fr;
    REp = p;
    distrp = d;
    rankK= REp->get_rankK();
    a_invgamma = a;
    b_invgamma = b;
    priorassumptions.push_back("Inverse gamma prior for variance component with hyperparameters a="
    +ST::doubletostring(a,6)+ " and b=" + ST::doubletostring(b,6) );
    priorassumptions.push_back("\\\\");

    if (average==true)
      setbeta(2,1,distrp->get_scale(column,column)/REp->getlambda());
    else
      setbeta(1,1,distrp->get_scale(column,column)/REp->getlambda());
    REp->update_sigma2(distrp->get_scale(column,column)/REp->getlambda());
    }


FULLCOND_variance_nonp::FULLCOND_variance_nonp(MCMCoptions * o,
                         FULLCOND_hrandom * p,DISTRIBUTION * d,
                         const double & a, const double & b,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr,const bool & av,
                         const unsigned & c)
                         : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
    {
    fctype = MCMC::variance;
    constlambda=false;
    uniformprior=false;
    Laplace=false;
    stationary=false;
    alphafix=false;
    discrete=false;
    update_sigma2 = true;
    randomeffect = false;
    hrandom=true;
    fullcondnonp = false;
    average = av;
    column = c;
    pathresults = fr;
    hREp = p;
    distrp = d;
    rankK= hREp->get_rankK();
    a_invgamma = a;
    b_invgamma = b;
    priorassumptions.push_back("Inverse gamma prior for variance component with hyperparameters a="
    +ST::doubletostring(a,6)+ " and b=" + ST::doubletostring(b,6) );
    priorassumptions.push_back("\\\\");

    if (average==true)
      setbeta(2,1,distrp->get_scale(column,column)/hREp->getlambda());
    else
      setbeta(1,1,distrp->get_scale(column,column)/hREp->getlambda());
    hREp->update_sigma2(distrp->get_scale(column,column)/hREp->getlambda());
    }



FULLCOND_variance_nonp::FULLCOND_variance_nonp(const FULLCOND_variance_nonp & t)
    : FULLCOND(FULLCOND(t))
  {
  constlambda=t.constlambda;
  uniformprior=t.uniformprior;
  Laplace=t.Laplace;
  stationary=t.stationary;
  alphafix=t.alphafix;
  randomeffect = t.randomeffect;
  fullcondnonp = t.fullcondnonp;
  average = t.average;
  update_sigma2 = t.update_sigma2;
  column = t.column;
  pathresults = t.pathresults;
  Kp = t.Kp;
  REp = t.REp;
  hREp = t.hREp;
  hrandom = t.hrandom;
  Fnp = t.Fnp;
  distrp = t.distrp;
  a_invgamma = t.a_invgamma;
  b_invgamma = t.b_invgamma;
  rankK = t.rankK;
  scale = t.scale;
  discrete = t.discrete;
  df = t.df;
  tau = t.tau;
  lambda = t.lambda;
  fc_lambda = t.fc_lambda;
  }


const FULLCOND_variance_nonp & FULLCOND_variance_nonp::operator=(
                                            const FULLCOND_variance_nonp & t)
  {
  if (this == &t)
    return *this;
  FULLCOND::operator=(FULLCOND(t));
  constlambda=t.constlambda;
  uniformprior=t.uniformprior;
  Laplace=t.Laplace;
  stationary=t.stationary;
  alphafix=t.alphafix;
  randomeffect = t.randomeffect;
  fullcondnonp = t.fullcondnonp;
  average = t.average;
  update_sigma2 = t.update_sigma2;
  column = t.column;
  pathresults = t.pathresults;
  Kp = t.Kp;
  REp = t.REp;
  hREp = t.hREp;
  hrandom = t.hrandom;
  Fnp = t.Fnp;
  distrp = t.distrp;
  a_invgamma = t.a_invgamma;
  b_invgamma = t.b_invgamma;
  rankK = t.rankK;
  scale = t.scale;
  discrete = t.discrete;
  df = t.df;
  tau = t.tau;
  lambda = t.lambda;
  fc_lambda = t.fc_lambda;
  return *this;
  }


void FULLCOND_variance_nonp::update(void)
  {

  transform = pow(distrp->get_trmult(column),2);

  acceptance++;

  if (fullcondnonp)
    {
    if (constlambda==false)
      {
      if(uniformprior)
        beta(0,0) = rand_invgamma(-0.5+0.5*rankK,0.5*Fnp->compute_quadform());
      else
        beta(0,0) = rand_invgamma(a_invgamma+0.5*rankK,
                                  b_invgamma+0.5*Fnp->compute_quadform());
      }
    else
      {
      beta(0,0) = distrp->get_scale(column,column)/Fnp->getlambda();
      }

    Fnp->update_sigma2(beta(0,0));
    }
  else if (randomeffect)
    {

    if (update_sigma2)
      {

      if (constlambda==false)
        {
        if(uniformprior==false)
          beta(0,0) = rand_invgamma(a_invgamma+0.5*rankK,
                                    b_invgamma+0.5*REp->compute_quadform());
        else
          beta(0,0) = rand_invgamma(-0.5+0.5*rankK,0.5*REp->compute_quadform());
        }
      else
        {
        beta(0,0) = distrp->get_scale(column,column)/REp->getlambda();
        }

      REp->update_sigma2(beta(0,0));

      }
    else
      {
      beta(0,0) = REp->get_sigma2();
      }

    }
  else if (hrandom)
    {

    if (update_sigma2)
      {

      if (constlambda==false)
        {
        if(uniformprior==false)
          beta(0,0) = rand_invgamma(a_invgamma+0.5*rankK,
                                    b_invgamma+0.5*hREp->compute_quadform());
        else
          beta(0,0) = rand_invgamma(-0.5+0.5*rankK,0.5*hREp->compute_quadform());
        }
      else
        {
        beta(0,0) = distrp->get_scale(column,column)/hREp->getlambda();
        }

      hREp->update_sigma2(beta(0,0));

      }
    else
      {
      beta(0,0) = hREp->get_sigma2();
      }

    }
  else
    {
    if (update_sigma2)
      {
      if (constlambda==false)
        {
        if(discrete)
          {
          unsigned i;
          double help,sum;
          vector<double> cumtau;

          int start = 1;

          if(Kp->get_type() == MCMC::RW1)
            start = 0;
          else if(Kp->get_type() == MCMC::RW2)
            start = 1;

          if(optionsp->get_nriter() == 1)
            {
            double init;
            double help1=(double)start;
            double help2=100.0;
            lambda.push_back(Kp->lambda_from_df(help1,help2));
            init = lambda[lambda.size()-1]/2;
            for(i=start+1;i<=df;i++)
              {
              if(init < 0.0000001)
                break;
              double help1=i-0.5;
              lambda.push_back(Kp->lambda_from_df(help1,init));
              init = lambda[lambda.size()-1]/2;
              if(init < 0.0000001)
                break;
              help1=(double)i;
              lambda.push_back(Kp->lambda_from_df(help1,init));
              init = lambda[lambda.size()-1]/2;
              }
            }

          scale = distrp->get_scale(column,column);

          tau = vector<double>(0);
          for(i=0;i<lambda.size();i++)
            {
            help = scale/lambda[i];
            tau.push_back(help);
            }

          i = 0;
          sum = 0.0;
          help = -0.5 * Kp->compute_quadform();
          while(i<tau.size())
            {
            if(uniformprior)
              sum += 1.0/pow(tau[i],0.5*rankK) * exp(help/tau[i]);
            else
              sum += 1.0/pow(tau[i],0.5*rankK) * exp(help/tau[i]) / double(start + 0.5 + i/2.0);
            cumtau.push_back(sum);
            i++;
            }

          help = uniform()*sum;

          i = 0;
          while(cumtau[i] < help)
            i++;

          beta(0,0) = tau[i];
          beta(1,0) = (start + 0.5 + i/2.0)/transform;            // --> df
          }
        else if(uniformprior)
          {
          beta(0,0) = rand_invgamma(-0.5+0.5*rankK,0.5*Kp->compute_quadform());
          }
        else if(Laplace)
          {
          beta(0,0) = rand_invgamma(a_invgamma+rankK,
                                    b_invgamma+Kp->compute_sumfabsdiff());
          }
        else
          {
          if(stationary)
            update_stationary();
          else
            beta(0,0) = rand_invgamma(a_invgamma+0.5*rankK,
                                      b_invgamma+0.5*Kp->compute_quadform());
          }

/*
        beta(0,0) = Kp->get_sigma2();

        double kappa = 1.0/beta(0,0);

        double lold = (a_invgamma-1)*log(kappa) - b_invgamma*kappa;           // gamma prior
        lold += -0.5*Kp->compute_quadform()*kappa;                            // priori für betas
        lold += 0.5*rankK*log(kappa);                                         // Normalisierungs Konstante

        double kappaprop = kappa*randnumbers::rand_variance(2.0);             // kappa ~ (1+1/z)

        double lnew = (a_invgamma-1)*log(kappaprop) - b_invgamma*kappaprop;   // gamma prior
        lnew += -0.5*Kp->compute_quadform()*kappaprop;                        // priori für betas
        lnew += 0.5*rankK*log(kappaprop);                                     // Normalisierungs Konstante

        double u = log(uniform());
        double alpha = lnew - lold;

        if(u <= alpha)
          beta(0,0) = 1.0/kappaprop;
        else
          acceptance--;
*/
        }
      else
        {
        beta(0,0) = distrp->get_scale(column,column)/Kp->getlambda();
        }

      Kp->update_sigma2(beta(0,0));
      }
    else
      {
      beta(0,0) = Kp->get_sigma2();
      }


    }

  if (average==true)
    {
    double help=0;
    while (help < 0.0000005 || help > 0.05)
      help = randnumbers::rand_gamma(a_invgamma+1,1.0/beta(0,0));

    b_invgamma = help;
    beta(1,0) = help;
    }

  FULLCOND::update();

  if(!fullcondnonp && !randomeffect && !hrandom)
    {
    double * lambdap = fc_lambda.getbetapointer();
    *lambdap = distrp->get_scale(column)/beta(0,0);
    fc_lambda.update();
    }

  }


void FULLCOND_variance_nonp::update_stationary(void)
  {

  double betaold = beta(0,0);
  beta(0,0) = rand_invgamma(a_invgamma+0.5*rankK,
                            b_invgamma+0.5*Kp->compute_quadform());

  if(Kp->get_type()==RW1)
    {
    double alphawork = beta(1,0)*transform;
//    double alphaprop = uniform()*2.0-1.0;
    double alphaprop;
    double help = uniform();
    if(help < (alphawork+1.0)/2.0)
      alphaprop = sqrt(2.0*help*(alphawork+1.0)) - 1.0;
    else
      alphaprop = 1.0 - sqrt(2.0*(help-1.0)*(alphawork-1.0));

    if(alphafix)
      alphaprop = alphawork;

    double u = log(uniform());
    double alpha = 0.5*log(1.0-alphaprop*alphaprop) - 0.5*log(1.0-alphawork*alphawork);

    alpha -= -0.5*Kp->compute_quadform()/beta(0,0);
    Kp->updateKenv_alpha(alphaprop);
    alpha += -0.5*Kp->compute_quadform()/beta(0,0);

    alpha += log( (pow(alphawork,2)-1.0) / (pow(alphaprop,2)-1.0) );

    if(u<=alpha)
      {
      beta(1,0) = alphaprop/transform;
      }
    else
      {
      Kp->updateKenv_alpha(alphawork);
//      beta(0,0) = betaold;
      }
    }
  else if(Kp->get_type()==RW2)
    {
    double alpha1work = beta(1,0)*transform;
    double alpha2work = beta(2,0)*transform;

//    double alpha1help = uniform()*2.0-1.0;
    double alphawork = sqrt(alpha2work);
    double alpha1help;
    double help = uniform();
    if(help < (alphawork+1.0)/2.0)
      alpha1help = sqrt(2.0*help*(alphawork+1.0)) - 1.0;
    else
      alpha1help = 1.0 - sqrt(2.0*(help-1.0)*(alphawork-1.0));
//    double alpha2help = uniform()*2.0-1.0;       // falls alpha1 und alpha2 unabhängig
    double alpha2help = alpha1help;

    double alpha1prop = -(alpha1help+alpha2help);
    double alpha2prop = alpha1help*alpha2help;

    if(alphafix)
      {
      alpha1prop = alpha1work;
      alpha2prop = alpha2work;
      }

    double u = log(uniform());

    double helpprop = 1.0-alpha2prop*alpha2prop;
    double helpwork = 1.0-alpha2work*alpha2work;

    double detKwork = helpwork * (helpwork - alpha1work*alpha1work * (1.0-alpha2work)/(1.0+alpha2work) );
    double detKprop = helpprop * (helpprop - alpha1prop*alpha1prop * (1.0-alpha2prop)/(1.0+alpha2prop) );

    double alpha = 0.5*log(detKprop)-0.5*log(detKwork);

    alpha -= -0.5*Kp->compute_quadform()/beta(0,0);
    Kp->updateKenv_alpha(alpha1prop,alpha2prop);
    alpha += -0.5*Kp->compute_quadform()/beta(0,0);

    alpha += log( (pow(alphawork,2)-1.0) / (pow(alpha1help,2)-1.0) );

    if(u<=alpha)
      {
      beta(1,0) = alpha1prop/transform;
      beta(2,0) = alpha2prop/transform;
      }
    else
      {
      Kp->updateKenv_alpha(alpha1work,alpha2work);
//      beta(0,0) = betaold;
      }
    }
  else if (Kp->get_type()==mrf)
    {
    double alphawork = beta(1,0)*transform;
//    double alphaprop = uniform()*2.0-1.0;
    double alphaprop;
    double help = uniform();
    if(help < (alphawork+1.0)/2.0)
      alphaprop = sqrt(2.0*help*(alphawork+1.0)) - 1.0;
    else
      alphaprop = 1.0 - sqrt(2.0*(help-1.0)*(alphawork-1.0));

    if(alphafix)
      alphaprop = alphawork;

    double u = log(uniform());
    double alpha = 0.0;

    alpha -= -0.5*Kp->compute_quadform()/beta(0,0);
    alpha -= 0.5*Kp->getLogDet();
    Kp->updateKenv_alpha(alphaprop);
    alpha += -0.5*Kp->compute_quadform()/beta(0,0);
    alpha += 0.5*Kp->getLogDet();

    alpha += log( (pow(alphawork,2)-1.0) / (pow(alphaprop,2)-1.0) );

    if(u<=alpha)
      {
      beta(1,0) = alphaprop/transform;
      }
    else
      {
      Kp->updateKenv_alpha(alphawork);
//      beta(0,0) = betaold;
      }
    }

  }


void FULLCOND_variance_nonp::set_stationary(double alphastart,bool afix)
  {
  alphafix = afix;
  stationary = true;
  if(Kp->get_type()==RW1)
    {
    setbeta(2,1,0.0);
    setbetavalue(0,0,distrp->get_scale(column,column)/Kp->getlambda());
    setbetavalue(1,0,alphastart/(distrp->get_trmult(column)*distrp->get_trmult(column)));
    }
  else if(Kp->get_type()==RW2)
    {
    setbeta(3,1,0.0);
    setbetavalue(0,0,distrp->get_scale(column,column)/Kp->getlambda());
    setbetavalue(1,0,-(alphastart+alphastart)/(distrp->get_trmult(column)*distrp->get_trmult(column)));
    setbetavalue(2,0,alphastart*alphastart/(distrp->get_trmult(column)*distrp->get_trmult(column)));
    }
  else if(Kp->get_type()==mrf)
    {
    setbeta(2,1,0.0);
    setbetavalue(0,0,distrp->get_scale(column,column)/Kp->getlambda());
    setbetavalue(1,0,alphastart/(distrp->get_trmult(column)*distrp->get_trmult(column)));
    }
  }


bool FULLCOND_variance_nonp::posteriormode(void)
  {

  if (fullcondnonp)
    {
    beta(0,0) = distrp->get_scale(column,column)/Fnp->getlambda();
    Fnp->update_sigma2(beta(0,0));
    }
  else if (randomeffect)
    {
    beta(0,0) = distrp->get_scale(column,column)/REp->getlambda();
    REp->update_sigma2(beta(0,0));
    }
  else if (hrandom)
    {
    beta(0,0) = distrp->get_scale(column,column)/hREp->getlambda();
    hREp->update_sigma2(beta(0,0));
    }
  else
    {
    beta(0,0) = distrp->get_scale(column,column)/Kp->getlambda();
    Kp->update_sigma2(beta(0,0));
    }

  return true;

  }


void FULLCOND_variance_nonp::outresults(void)
  {
  FULLCOND::outresults();

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

  if (optionsp->get_samplesize()>0)
    {
    optionsp->out("  Estimation results for the variance component:\n");
    optionsp->out("\n");

    vstr = "  Mean:         ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betamean(0,0),6) + "\n");

    vstr = "  Std. dev.:    ";
    if (constlambda==false)
      {
      optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
      ST::doubletostring(sqrt(betavar(0,0)),6) + "\n");
      }
    else
      {
      optionsp->out(vstr + ST::string(' ',20-vstr.length()) + "0 \n");
      }

    vstr = "  " + l1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_lower(0,0),6) + "\n");

    vstr = "  " + l2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_lower(0,0),6) + "\n");

    vstr = "  50% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu50(0,0),6) + "\n");

    vstr = "  " + u1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_upper(0,0),6) + "\n");

    vstr = "  " + u2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_upper(0,0),6) + "\n");

    optionsp->out("\n");

    ofstream ou(pathresults.strtochar());

    ou << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
    nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
    ou << betamean(0,0) << "  ";
    ou << (betavar(0,0)<0.0?0.0:sqrt(betavar(0,0))) << "  ";
    ou << betaqu_l1_lower(0,0) << "  ";
    ou << betaqu_l2_lower(0,0) << "  ";
    ou << betaqu50(0,0) << "  ";
    ou << betaqu_l2_upper(0,0) << "  ";
    ou << betaqu_l1_upper(0,0) << "  ";
    ou << betamin(0,0) << "  ";
    ou << betamax(0,0) << "  " << endl;

    if(beta.rows()>1)
      {
      ou << betamean(1,0) << "  ";
      ou << (betavar(1,0)<0.0?0.0:sqrt(betavar(1,0))) << "  ";
      ou << betaqu_l1_lower(1,0) << "  ";
      ou << betaqu_l2_lower(1,0) << "  ";
      ou << betaqu50(1,0) << "  ";
      ou << betaqu_l2_upper(1,0) << "  ";
      ou << betaqu_l1_upper(1,0) << "  ";
      ou << betamin(1,0) << "  ";
      ou << betamax(1,0) << "  " << endl;
      }
    if(beta.rows()>2)
      {
      ou << betamean(2,0) << "  ";
      ou << (betavar(2,0)<0.0?0.0:sqrt(betavar(2,0))) << "  ";
      ou << betaqu_l1_lower(2,0) << "  ";
      ou << betaqu_l2_lower(2,0) << "  ";
      ou << betaqu50(2,0) << "  ";
      ou << betaqu_l2_upper(2,0) << "  ";
      ou << betaqu_l1_upper(2,0) << "  ";
      ou << betamin(2,0) << "  ";
      ou << betamax(2,0) << "  " << endl;
      }

    optionsp->out("  Results for the variance component are also stored in file\n");
    optionsp->out("  " + pathresults + "\n");

    optionsp->out("\n");


    }
  else
    {

    transform = pow(distrp->get_trmult(column),2);

    optionsp->out("  Variance parameter: "
                  + ST::doubletostring(transform*beta(0,0),6) + "\n");
    optionsp->out("\n");
    }

  if(!fullcondnonp && !randomeffect && !hrandom)
    outresults_lambda();

  }


void FULLCOND_variance_nonp::outresults_lambda(void)
  {

  fc_lambda.outresults();

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

  if (optionsp->get_samplesize()>0)
    {
    optionsp->out("  Estimation results for the smoothing parameter:\n");
    optionsp->out("\n");

    Kp->update_stepwise(fc_lambda.get_betamean(0,0));

    vstr = "  Mean:         ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_lambda.get_betamean(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length())
                       + "(df: " + ST::doubletostring(Kp->compute_df(),6) + ")\n");

    vstr = "  Std. dev.:    ";
    if (constlambda==false)
      {
      optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
      ST::doubletostring(sqrt(fc_lambda.get_betavar(0,0)),6) + "\n");
      }
    else
      {
      optionsp->out(vstr + ST::string(' ',20-vstr.length()) + "0 \n");
      }

    Kp->update_stepwise(fc_lambda.get_beta_lower1(0,0));

    vstr = "  " + l1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_lambda.get_beta_lower1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length())
                       + "(df: " + ST::doubletostring(Kp->compute_df(),6) + ")\n");

    Kp->update_stepwise(fc_lambda.get_beta_lower2(0,0));

    vstr = "  " + l2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_lambda.get_beta_lower2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length())
                       + "(df: " + ST::doubletostring(Kp->compute_df(),6) + ")\n");

    Kp->update_stepwise(fc_lambda.get_betaqu50(0,0));

    vstr = "  50% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_lambda.get_betaqu50(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length())
                       + "(df: " + ST::doubletostring(Kp->compute_df(),6) + ")\n");

    Kp->update_stepwise(fc_lambda.get_beta_upper2(0,0));

    vstr = "  " + u1 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_lambda.get_beta_upper2(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length())
                       + "(df: " + ST::doubletostring(Kp->compute_df(),6) + ")\n");

    Kp->update_stepwise(fc_lambda.get_beta_upper1(0,0));

    vstr = "  " + u2 + "% Quantile: ";
    vstr = vstr + ST::string(' ',20-vstr.length())
                + ST::doubletostring(fc_lambda.get_beta_upper1(0,0),6);
    optionsp->out(vstr + ST::string(' ',40-vstr.length())
                       + "(df: " + ST::doubletostring(Kp->compute_df(),6) + ")\n");

    optionsp->out("\n");

    ST::string lambda_pathresults = pathresults.substr(0,pathresults.length()-7) + "lambda.res";

    ofstream ou(lambda_pathresults.strtochar());

    ou << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
    nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
    ou << fc_lambda.get_betamean(0,0) << "  ";
    if (constlambda==false)
      {
      ou << sqrt(fc_lambda.get_betavar(0,0)) << "  ";
      }
    else
      {
      ou << 0 << "  ";
      }
    ou << fc_lambda.get_beta_lower1(0,0) << "  ";
    ou << fc_lambda.get_beta_lower2(0,0) << "  ";
    ou << fc_lambda.get_betaqu50(0,0) << "  ";
    ou << fc_lambda.get_beta_upper2(0,0) << "  ";
    ou << fc_lambda.get_beta_upper1(0,0) << "  ";
    ou << fc_lambda.get_betamin(0,0) << "  ";
    ou << fc_lambda.get_betamax(0,0) << "  " << endl;

    optionsp->out("  Results for the smoothing parameter are also stored in file\n");
    optionsp->out("  " + lambda_pathresults + "\n");

    optionsp->out("\n");


    }
  else
    {

    transform = pow(distrp->get_trmult(column),2);

    optionsp->out("  Smoothing parameter: "
                  + ST::doubletostring(transform*distrp->get_scale(0,0)/beta(0,0),6) + "\n");
    optionsp->out("\n");
    }

/*
  double lambda;
  for(int i=1;i<=3;i++)
    {
    lambda = Kp->lambda_from_df(i,exp(5.0/i));
    optionsp->out("  " + ST::inttostring(i) + " degree(s) of freedom approximately correspond(s) to lambda = "
                    + ST::doubletostring(lambda,6) + "\n");
    Kp->update_stepwise(lambda);
    optionsp->out("  (" + ST::doubletostring(Kp->compute_df(),6) + ")\n");
    }
  optionsp->out("\n");
*/  
/*
  ST::string file = pathresults.substr(0,pathresults.length()-7) + "lambda_sample.raw";
  fc_lambda.get_samples(file);
*/
  }


void FULLCOND_variance_nonp::outoptions(void)
  {

  if(discrete)
    {
    if(uniformprior)
      optionsp->out("  Discrete uniform prior on approximate degrees of freedom\n");
    else
      optionsp->out("  Discrete geometric prior on approximate degrees of freedom\n");
    }
  else if(uniformprior)
    {
    optionsp->out("  Uniform prior on tau\n");
    }
  else
    {
    optionsp->out("  Hyperprior a for variance parameter: " +
                ST::doubletostring(a_invgamma) + "\n" );
    if (!average)
      {
      optionsp->out("  Hyperprior b for variance parameter: " +
                  ST::doubletostring(b_invgamma) + "\n" );
      }
    else
      optionsp->out("  Hyperprior b for variance parameter: uniform prior\n");
    }

  optionsp->out("\n");
  }


} // end: namespace MCMC




