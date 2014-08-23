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



#include "fullcond_nonp_gaussian_stepwise.h"

namespace MCMC
{

// additive Effekte, RW1 RW2 und season

FULLCOND_nonp_gaussian_stepwise::FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                      DISTRIBUTION * dp,
                      const datamatrix & d,
                      FULLCOND_const * fcc,
                      const unsigned & maxint,const fieldtype & ft,
                      const ST::string & ti,
                      const ST::string & fp, const ST::string & pres,
                      const unsigned & c,const double & l,
                      const unsigned & per)
  : FULLCOND_nonp_gaussian(o,dp,d,fcc,maxint,ft,ti,fp,pres,c,l,per)
  {

  isbootstrap = false;
  intercept = 0.0;

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  spatialtotal = false;
  }

// varying coefficients , RW1 RW2 und season

FULLCOND_nonp_gaussian_stepwise::FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,DISTRIBUTION * dp,
                       const datamatrix & d,
                       const datamatrix & intvar,
                       FULLCOND_const * fcc,
                       const unsigned & maxint,
                       const fieldtype & ft,const ST::string & ti,
                       const ST::string & fp, const ST::string & pres,
                       const unsigned & c,const double & l, const bool & vccent,
                       const unsigned & per)
  : FULLCOND_nonp_gaussian(o,dp,d,intvar,fcc,maxint,ft,ti,fp,pres,c,l,vccent,per)

  {

  isbootstrap = false;
  intercept = 0.0;
  //VCM_neu
  if(vccent == true)
    identifiable = false;

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  get_data_forfixedeffects();
  effmodi = d;
  unsigned nrobs = index.rows();

  if (data_varcoeff_fix.rows() < nrobs)
    {
    data_varcoeff_fix = datamatrix(nrobs,2,1);
    for(unsigned i=0;i<nrobs;i++)
      {
      data_varcoeff_fix(i,0) = intvar(i,0);
      data_varcoeff_fix(i,1) = d(i,0)*intvar(i,0);
      }
    }

  XVX = datamatrix(2,2,0);

  spatialtotal = false;
  }

// spatial covariates

FULLCOND_nonp_gaussian_stepwise::FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                        DISTRIBUTION * dp,const datamatrix & d,
                        FULLCOND_const * fcc,
                        const MAP::map & m, const ST::string & mn,
                        const ST::string & ti,
                        const ST::string & fp, const ST::string & pres,
                        const unsigned & c,const double & l,
                        const fieldtype & ft, const MAP::map & m2)
  : FULLCOND_nonp_gaussian(o,dp,d,fcc,m,mn,ti,fp,pres,c,l)

  {

  isbootstrap = false;
  intercept = 0.0;

  kombimatrix = false;
  matrixnumber = 1;

  type = ft;
  if(type == mrfI)
    {
    nofixed = true;
    forced_into = true;
    kombimatrix = true;
    numberofmatrices = 2;
    kappa = vector<double>(1,1);
    kappaold = vector<double>(1,-2);
    kappa_prec = vector<double>(1,-1);
    Kenv2 = Krw0env(nrpar);
    }

  if(type == twomrfI)
    {
    nofixed = true;
    forced_into = true;
    kombimatrix = true;
    numberofmatrices = 3;
    kappa = vector<double>(1,1);
    kappaold = vector<double>(1,-2);
    kappa_prec = vector<double>(1,-1);
    Kenv2 = Kmrfenv(m2);
    Kenv3 = Krw0env(nrpar);
    }

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  spatialtotal = false;
  }

// varying coefficients , spatial covariates as effect modifier

FULLCOND_nonp_gaussian_stepwise::FULLCOND_nonp_gaussian_stepwise(MCMCoptions * o,
                        DISTRIBUTION * dp,
                        FULLCOND_const * fcc,
                        const MAP::map & m,
                        const ST::string & mn,
                        const datamatrix & d,
                        const datamatrix & d2,
                        const ST::string & ti,
                        const ST::string & fp, const ST::string & pres,
                        const unsigned & c, const double & l, const bool & vccent,
                        const fieldtype & ft)
  : FULLCOND_nonp_gaussian(o,dp,fcc,m,mn,d,d2,ti,fp,pres,c,l,vccent)

  {

  isbootstrap = false;
  intercept = 0.0;
  //VCM_neu
  if(vccent == true)
    identifiable = false;

  get_data_forfixedeffects();

  type = ft;
  if(type == mrfI)
    {
    nofixed = true;
    forced_into = true;
    kombimatrix = true;
    numberofmatrices = 2;
    kappa = vector<double>(1,1);
    kappaold = vector<double>(1,-2);
    kappa_prec = vector<double>(1,-1);
    Kenv2 = Krw0env(nrpar);
    }

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  spatialtotal = false;
  }


 void FULLCOND_nonp_gaussian_stepwise::init_spatialtotal(FULLCOND * unstructp)
  {
  spatialtotal = true;
  fcunstruct = unstructp;
  }

  // COPY CONSTRUCTOR

FULLCOND_nonp_gaussian_stepwise::FULLCOND_nonp_gaussian_stepwise(const FULLCOND_nonp_gaussian_stepwise & fc)
  : FULLCOND_nonp_gaussian(FULLCOND_nonp_gaussian(fc))
  {
  intercept = fc.intercept;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  XVX = fc.XVX;
  fcunstruct = fc.fcunstruct;
  spatialtotal = fc.spatialtotal;
  df_lambdaold = fc.df_lambdaold;
  lambdaold = fc.lambdaold;
  df_lambdaold_unstr = fc.df_lambdaold_unstr;
  lambdaold_unstr = fc.lambdaold_unstr;
  lambdavec = fc.lambdavec;
  all_precenv = fc.all_precenv;
  fc_df = fc.fc_df;
  isbootstrap = fc.isbootstrap;
  Kenv2 = fc.Kenv2;
  Kenv3 = fc.Kenv3;
  kappa = fc.kappa;
  kappaold = fc.kappaold;
  kappa_prec = fc.kappa_prec;
  otherfullcond = fc.otherfullcond;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const FULLCOND_nonp_gaussian_stepwise & FULLCOND_nonp_gaussian_stepwise::operator=(
                                            const FULLCOND_nonp_gaussian_stepwise & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_nonp_gaussian::operator=(FULLCOND_nonp_gaussian(fc));

  intercept = fc.intercept;
  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  XVX = fc.XVX;
  fcunstruct = fc.fcunstruct;
  spatialtotal = fc.spatialtotal;
  df_lambdaold = fc.df_lambdaold;
  lambdaold = fc.lambdaold;
  df_lambdaold_unstr = fc.df_lambdaold_unstr;
  lambdaold_unstr = fc.lambdaold_unstr;
  lambdavec = fc.lambdavec;
  all_precenv = fc.all_precenv;
  fc_df = fc.fc_df;
  isbootstrap = fc.isbootstrap;
  Kenv2 = fc.Kenv2;
  Kenv3 = fc.Kenv3;
  kappa = fc.kappa;
  kappaold = fc.kappaold;
  kappa_prec = fc.kappa_prec;
  otherfullcond = fc.otherfullcond;

  return *this;
  }


bool FULLCOND_nonp_gaussian_stepwise::posteriormode(void)
  {
if(kombimatrix == true)
  {
  return posteriormode_kombi();
  }
else
  {
  int j;
  unsigned i;

  int * workindex;

  update_linpred(false);

  // NEU!!!
  if(varcoeff && lambda == -2)
    {
    //datamatrix X = datamatrix(2,2,0);
    datamatrix betas = datamatrix(2,1,0);

    if(calculate_xwx_vc == true || (XVX(1,1)==0 && XVX(0,0)==0))
      {
      calculate_xwx_vc = false;
      likep->fisher(XVX,data_varcoeff_fix,column);            // recomputes X1 = (newx' W newx)^{-1}
      XVX.assign((XVX.cinverse()));               // continued
      }
    likep->compute_weightiwls_workingresiduals(column); // computes W(y-linpred)
    betas = XVX*data_varcoeff_fix.transposed()*likep->get_workingresiduals();

    double * workbeta = beta.getV();
    vector<int>::iterator itbeg = posbeg.begin();
    vector<int>::iterator itend = posend.begin();
    int * workindex = index.getV();
    for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        *workbeta = betas(0,0) + effmodi(*workindex,0)*betas(1,0);
        for(j=*itbeg;j<=*itend;j++)
          workindex++;
        }
      }

    update_linpred(true);  // addiert Fkt. zum Gesamtprädiktor

    if(center)
      {
      intercept = centerbeta();
      update_fix_effect(intercept);
      intercept = 0.0;
      }
    }
  else
    {
    if (lambda_prec != lambda || calculate_xwx == true)
      {
      if (calculate_xwx == true)
        {
        calculate_xwx = false;
        if (varcoeff)
          compute_XWX_varcoeff_env(likep->get_weightiwls(),column);
        else
          compute_XWX_env(likep->get_weightiwls(),column);
        }
      precenv.addtodiag(XXenv,Kenv,1.0,lambda);
      lambda_prec = lambda;
      }

    likep->compute_weightiwls_workingresiduals(column);

    workindex = index.getV();
    double * workmuy = beta.getV();

    if (varcoeff)
      {
      double * workdata=data.getV();
      for(i=0;i<nrpar;i++,workmuy++)
        {
        *workmuy = 0;
        if (posbeg[i] != -1)
          for(j=posbeg[i];j<=posend[i];j++,workindex++,workdata++)
            *workmuy+=
            likep->get_workingresiduals()(*workindex,0)*(*workdata);
        }
      }
    else  // else additive
      {
      for(i=0;i<nrpar;i++,workmuy++)
        {
        *workmuy = 0;
        if (posbeg[i] != -1)
          for(j=posbeg[i];j<=posend[i];j++,workindex++)
            *workmuy+= likep->get_workingresiduals()(*workindex,0);
        }
      }

    precenv.solve(beta);
    update_linpred(true);

    if (center)
      {
      intercept = centerbeta();
      if(varcoeff == false)
        {
        fcconst->posteriormode_intercept(intercept);
        }
      else if(varcoeff == true)
        {
        update_fix_effect(intercept);
        intercept = 0.0;
        }
      }
    } // END: else if(!varcoeff || lambda != -2)

  transform = likep->get_trmult(column);

  return FULLCOND_nonp_basis::posteriormode();
  }
  }


bool FULLCOND_nonp_gaussian_stepwise::posteriormode_converged(const unsigned & itnr)
  {
  if(kombimatrix == false || matrixnumber==1)
    return FULLCOND_nonp_gaussian::posteriormode_converged(itnr);
  else
    return true;
  }


const datamatrix & FULLCOND_nonp_gaussian_stepwise::get_data_forfixedeffects(void)
  {

  if ( (data_forfixed.rows() < index.rows()) &&
       (!varcoeff) &&
       ( (type==RW1) || (type==RW2) )
     )
    {
    data_forfixed=datamatrix(index.rows(),1);
    unsigned i;
    int j;
    int * workindex = index.getV();
    double h;
    for(i=0;i<posbeg.size();i++)
      {
      h = effectvdouble[i];
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++)
          {
          data_forfixed(*workindex,0) = h;
          }
      }

    }
  else if ( (data_forfixed.rows() < data.rows()) && (varcoeff==true) )
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


void FULLCOND_nonp_gaussian_stepwise::update_stepwise(double la)
  {
  lambda=la;

  if(likep->iwlsweights_constant() == true && kombimatrix==false)
    {
    bool gefunden = false;
    unsigned i = 0;
    while(i<lambdavec.size() && gefunden == false)
      {
      if(lambda == lambdavec[i])
        gefunden = true;
      i++;
      }
    if(gefunden == true)
      {
      precenv = all_precenv[i-1];
      lambda_prec = lambda;
      }
    }
  }


double FULLCOND_nonp_gaussian_stepwise::compute_df(void)
  {
if(kombimatrix == true)
  {
  return compute_df_kombi();
  }
else
  {
  double df = 0;
  if(inthemodel == true)
    {
    if(varcoeff && lambda == -2)
      {
      if(identifiable)
        df = 2;
      else
        df = df + 1;
      }
    else
      {
      bool unstr_included = false;
      if(spatialtotal)
        {
        bool fix;
        fcunstruct->get_inthemodel(unstr_included,fix);
        if(unstr_included == true)
          {
          double lambda_unstr = fcunstruct->get_lambda();
          double df_str = 0;
          double df_unstr = 0;

          if (lambdaold == lambda && lambdaold_unstr == lambda_unstr && likep->get_iwlsweights_notchanged() == true)
            {
            df = df_lambdaold;
            fcunstruct->set_dfunstruct(df_lambdaold_unstr);
            }
          else
            {
            if (calculate_xwx == true)
              {
              calculate_xwx = false;
              compute_XWX_env(likep->get_weightiwls(),column);
              }

            envmatdouble Diag_neu = envmatdouble(0,nrpar);
            vector<double>::iterator d = XXenv.getDiagIterator();
            vector<double>::iterator d2 = Diag_neu.getDiagIterator();
            unsigned i;
            for(i=0;i<nrpar;i++,d++,d2++)
              *d2 = *d * lambda_unstr / (*d + lambda_unstr);
            envmatdouble precenv_neu = envmatdouble(Kenv.getXenv(),0,nrpar);
            precenv_neu.addtodiag(Diag_neu,Kenv,1.0,lambda);
            invprec = envmatdouble(precenv_neu.getXenv(),0,precenv_neu.getDim());
            precenv_neu.inverse_envelope(invprec);

            df_str += invprec.traceOfProduct(XXenv);

            d = XXenv.getDiagIterator();
            d2 = Diag_neu.getDiagIterator();
            for(i=0;i<nrpar;i++,d++,d2++)
              *d2 = *d * *d / (*d + lambda_unstr);

            df_unstr -= invprec.traceOfProduct(Diag_neu);
            df_str += df_unstr;

            d = XXenv.getDiagIterator();
            d2 = Diag_neu.getDiagIterator();
            for(i=0;i<nrpar;i++,d++,d2++)
              *d2 = *d2 * *d / (*d + lambda_unstr);

            df_unstr += invprec.traceOfProduct(Diag_neu);

            d = XXenv.getDiagIterator();
            for(i=0;i<nrpar;i++,d++,d2++)
              df_unstr += *d / (*d + lambda_unstr);

            fcunstruct->set_dfunstruct(df_unstr);
            df = df_str-1;
            df_lambdaold = df;
            df_lambdaold_unstr = df_unstr;
            lambdaold = lambda;
            lambdaold_unstr = lambda_unstr;
          //  }
            }
          }
        }

      if(!spatialtotal || !unstr_included)
        {

        if (lambdaold == lambda && likep->get_iwlsweights_notchanged() == true && !spatialtotal)
          {
          df = df_lambdaold;
          }
        else if (spatialtotal || lambdaold != lambda || likep->get_iwlsweights_notchanged() == false)
          {
          if (calculate_xwx == true)
            {
            if (varcoeff)
              compute_XWX_varcoeff_env(likep->get_weightiwls(),column);
            else
              compute_XWX_env(likep->get_weightiwls(),column);
            }
          if(lambda != lambda_prec || calculate_xwx == true)
            {
            calculate_xwx = false;
            precenv.addtodiag(XXenv,Kenv,1.0,lambda);
            lambda_prec = lambda;
            }

          if(type!=MCMC::seasonal)
            {
            if (type==MCMC::mrf || type==MCMC::mrfI || type==MCMC::twomrfI)
              invprec = envmatdouble(precenv.getXenv(),0,precenv.getDim());
            else
              invprec = envmatdouble(0,nrpar,Kenv.getBandwidth());

            precenv.inverse_envelope(invprec);

            if(identifiable)
              df = invprec.traceOfProduct(XXenv);
            else
              df = df + invprec.traceOfProduct(XXenv)-1;
            }
          else    // Kombination: MRF + I
            {
            datamatrix Z = datamatrix(nrpar,nrpar,0);
            vector<double>::iterator d = XXenv.getDiagIterator();
            vector<double>::iterator d2 = XXenv.getDiagIterator();
            vector<double>::iterator p = Kenv.getDiagIterator();
            double sumn = 0; // enthält Summe der Elemente von XXenv

            unsigned i,j;
            for(i=0;i<nrpar;i++,d++)
              {
              sumn += *d;
              }

            d = XXenv.getDiagIterator();
            double * workz = Z.getV();

            for(i=0;i<nrpar;i++,d++,p++,workz++)
              {
              *workz = *d - 1/sumn * *d * *d + lambda * *p;
              d2 = XXenv.getDiagIterator() + i+1;
              workz++;
              for(j=i+1;j<nrpar;j++,workz++,d2++)
                {
                *workz = -1/sumn * *d * *d2;
                }
             workz += i;
             }

            p = Kenv.getEnvIterator();
            double bandw = Kenv.getBandwidth();

            for(i=1;i<nrpar;i++)
              {
              unsigned oben = unsigned(bandw);
              if(i <= unsigned(bandw))
                {
                oben = i;
                workz = Z.getV() + i;
                }
              else
                workz = Z.getV() + i + (i - unsigned(bandw))*nrpar;
              for(j=0;j<oben;j++,workz+=nrpar,p++)
                {
                *workz = *workz + lambda * *p;
                }
              }

            for(i=0;i<nrpar;i++)    // schreibt untere Hälfte von Z voll
              {
              for(j=i+1;j<nrpar;j++)
                {
                Z(j,i) = Z(i,j);
                }
              }

            Z.assign(Z.inverse());    // berechnet Inverse von (n1-n1^2/n, -n1*n2/n,..., -n1*np/n) + lambda*P

            workz = Z.getV();
            d = XXenv.getDiagIterator();

            df = 1;
            for(i=0;i<nrpar;i++,d++,workz++)
              {
              df += *d * *workz - 1/sumn * *d * *d * *workz;
              p = XXenv.getDiagIterator() + i+1;
              workz++;
              for(j=i+1;j<nrpar;j++,p++,workz++)
                {
                df -= 2/sumn * *d * *p * *workz;
                }
              workz += i;
              }

            df -= 1;
            }

          df_lambdaold = df;
          lambdaold = lambda;
          }
        }
      }
    }

  return df;
  }
  }


ST::string FULLCOND_nonp_gaussian_stepwise::get_effect(void)
  {
  ST::string h = "";

if(matrixnumber == 1)
  {
  ST::string t;
  if (type==MCMC::RW1)
    t = "rw1";
  else if (type==MCMC::RW2)
    t = "rw2";
  else if (type==MCMC::seasonal)
    t = "seasonal";
  else if (type==MCMC::mrf)
    t = "spatial";
  else if (type==MCMC::mrfI)
    t = "spatialrandom";
  else if (type==MCMC::twomrfI)
    t = "twospatialrandom";

  if(varcoeff)
    h = datanames[1] + "*" + datanames[0];   // (VC alt) 0 - 1
  else
    h = datanames[0];

  if(type != mrfI && type != twomrfI && type != twomrfI)
    h = h + "(" + t + ",df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";
  else if(type == mrfI)
    h = h + "(" + t + ",df=" + ST::doubletostring(compute_df(),6) +
        ",(lambda1=" + ST::doubletostring(lambda,6) + "),(lambda2=" + ST::doubletostring(otherfullcond[0]->get_lambda(),6) + "))";
  else
    h = h + "(" + t + ",df=" + ST::doubletostring(compute_df(),6) +
        ",(lambda1=" + ST::doubletostring(lambda,6) + "),(lambda2=" +
        ST::doubletostring(otherfullcond[0]->get_lambda(),6) + "),(lambda3=" +
        ST::doubletostring(otherfullcond[1]->get_lambda(),6) + "))";
  }

  return h;
  }


void FULLCOND_nonp_gaussian_stepwise::init_names(const vector<ST::string> & na)
    {
    FULLCOND::init_names(na);

    char hchar = '_';
    ST::string hstring =  "\\_";

    if (na.size()==1)
      {
      ST::string helpname = na[0].insert_string_char(hchar,hstring);
      if (type==MCMC::seasonal)
        term_symbolic = "f^{Season}_{" +  helpname + "}("+helpname+")";
      else
        term_symbolic = "f_{" +  helpname + "}("+helpname+")";
      }
    else
      {
      ST::string helpname1 = na[0].insert_string_char(hchar,hstring);   // (VC alt) 1
      ST::string helpname2 = na[1].insert_string_char(hchar,hstring);   // (VC alt) 0
      if (type==MCMC::seasonal)
        term_symbolic = "f^{Season}_{" +  helpname1 + "}("+helpname1+") \\cdot "
                        + helpname2;
      else
        term_symbolic = "f_{" +  helpname1 + "}("+helpname1+") \\cdot "
                        + helpname2;
      }

    if (column > 0)
      priorassumptions.push_back("$" + term_symbolic + "$" +
       " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");
//    priorassumptions.push_back(term_symbolic);
    init_priorassumptions(na[0]);
    }


void FULLCOND_nonp_gaussian_stepwise::init_priorassumptions(const ST::string & na)
    {
    if(type==MCMC::RW1)
       priorassumptions.push_back("first order difference penalty");
    else if(type==MCMC::RW2)
       priorassumptions.push_back("second order difference penalty");
    else if(type==MCMC::mrf)
       priorassumptions.push_back("spatial pairwise difference penalty");
    else if(type==MCMC::RE)
       priorassumptions.push_back("random effect");
    else if(type==MCMC::seasonal)
       priorassumptions.push_back("time varying seasonal component");
    else if(type==MCMC::smoothspline)
       priorassumptions.push_back("smoothing spline");
    else if(type==MCMC::mrfkronecker)
       priorassumptions.push_back("Kronecker product interaction");
    else if(type==MCMC::mrflinear)
       priorassumptions.push_back("2 dimensional first order difference penalty");
    else if(type==MCMC::mrfkr1)
       priorassumptions.push_back("Kronecker product interaction (RW1*RW1)");
    else if(type==MCMC::mrfkr2)
       priorassumptions.push_back("Kronecker product interaction (RW2*RW2)");
    }




void FULLCOND_nonp_gaussian_stepwise::reset_effect(const unsigned & pos)
  {

  update_linpred(false);

  unsigned i;
  double * work;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  intercept = 0.0;
  }


void FULLCOND_nonp_gaussian_stepwise::hierarchie_rw1(vector<double> & untervector, int dfo)
  {

  unsigned number = untervector.size()-1;

  update_stepwise(untervector[0]);
  double df_max = compute_df();

  update_stepwise(untervector[number]);
  double df_min = compute_df();

  if(df_max > dfo && df_min < dfo)
     {
     bool geordnet = false;
     unsigned stelle_oben = number;
     unsigned stelle_unten = 0;
     while(geordnet==false)
        {
        unsigned stelle = stelle_oben + stelle_unten;
        update_stepwise(untervector[stelle/2]);
        double df_mitteunten = compute_df();
        update_stepwise(untervector[stelle/2 + 1]);
        double df_mitteoben = compute_df();

        if(df_mitteunten > dfo && df_mitteoben > dfo)
          stelle_unten = stelle/2;
        else if(df_mitteunten < dfo && df_mitteoben < dfo)
          stelle_oben = stelle/2 + 1;
        else
          {
          geordnet = true;
          vector<double> hilf;
          unsigned i;
          stelle_unten = stelle/2;
          stelle_oben = stelle/2 + 1;
          for(i=0;i<=stelle_unten;i++)
             hilf.push_back(untervector[i]);
          if(!varcoeff)
            hilf.push_back(-1);
          else
            hilf.push_back(-2);
          for(i=stelle_oben;i<untervector.size();i++)
            hilf.push_back(untervector[i]);
          untervector = hilf;
          }
        }
     }
  else if(df_min >= dfo)
     {
     if(!varcoeff)
       untervector.push_back(-1);
     else
       untervector.push_back(-2);
     }
  else
     {
     vector<double> hilf;
     if(!varcoeff)
       hilf.push_back(-1);
     else
       hilf.push_back(-2);
     unsigned i;
     for(i=0;i<untervector.size();i++)
        hilf.push_back(untervector[i]);
     untervector = hilf;
     }
  }


void FULLCOND_nonp_gaussian_stepwise::compute_lambdavec(
vector<double> & lvec, int & number)
  {
  lambdaold = -1;
  unsigned i;
  bool kchanged = false;
  if(kombimatrix==true)
    {
    nofixed = true;
    forced_into = true;
    for(i=0;i<otherfullcond.size();i++)
      otherfullcond[i]->update_stepwise(0.0000001);
    kchanged = true;
    kombimatrix = false;    // damit das Auswählen der Lambdas schneller geht
    }

  if(spfromdf=="automatic")
    {
    df_equidist = true;
    double maxi = floor(nrpar/4.0*3.0);

    if(type != seasonal)
      {
      if(maxi <= 30)
        {
        df_for_lambdamin = maxi;
        if( (!varcoeff || !identifiable) && type==mrf)
          {
          df_for_lambdamax = 1;
          number = maxi;
          }
        else if(varcoeff && identifiable && type!=mrf)
          {
          df_for_lambdamax = 3;
          number = maxi - 2;
          }
        else
          {
          df_for_lambdamax = 2;
          number = maxi - 1;
          }
        }
      else if(maxi > 30 && maxi<=60)
        {
        number = floor(maxi/2);
        if( (!varcoeff || !identifiable) && type==mrf)
          {
          df_for_lambdamax = 2;
          df_for_lambdamin = number*2;
          }
        else
          {
          df_for_lambdamax = 3;
          df_for_lambdamin = number*2+1;
          }
        }
      else if(maxi > 60 && maxi<=100)
        {
        df_for_lambdamax = 3;
        number = floor(maxi/3);
        df_for_lambdamin = number*3;
        }
      else if(maxi > 100 && maxi<=180)
        {
        df_for_lambdamax = 5;
        number = floor(maxi/5);
        df_for_lambdamin = number*5;
        }
      else if(maxi > 180)
        {
        df_for_lambdamax = 10;
        number = floor(maxi/10);
        df_for_lambdamin = number*10;
        }
      }
    else
      {
      double diff = maxi-period+1;
      if(diff <= 30)
        {
        df_for_lambdamin = maxi;
        df_for_lambdamax = period-1+0.001;
        number = diff;
        }
      else if(diff > 30 && diff<=60)
        {
        number = floor(diff/2);
        df_for_lambdamin = period-1 + number*2;
        df_for_lambdamax = period-1+0.001;
        }
      else if(diff > 60 && diff<=100)
        {
        number = floor(diff/3);
        df_for_lambdamin = period-1 + number*3;
        df_for_lambdamax = period-1+0.001;
        }
      else if(diff > 100 && diff<=180)
        {
        number = floor(diff/5);
        df_for_lambdamin = period-1 + number*5;
        df_for_lambdamax = period-1+0.001;
        }
      else if(diff > 180)
        {
        number = floor(diff/10);
        df_for_lambdamin = period-1 + number*10;
        df_for_lambdamax = period-1+0.001;
        }
      }
    }

  if (df_equidist==true && spfromdf!="direct" && number>1)
     FULLCOND::compute_lambdavec_equi(lvec,number);
  else
     FULLCOND::compute_lambdavec(lvec,number);

  if(likep->iwlsweights_constant() == true && kombimatrix == false)
    {
    lambdavec = lvec;
    if (varcoeff)
      compute_XWX_varcoeff_env(likep->get_weightiwls(),column);
    else
      compute_XWX_env(likep->get_weightiwls(),column);
    for(unsigned i=0;i<lambdavec.size();i++)
      {
      precenv.addtodiag(XXenv,Kenv,1.0,lambdavec[i]);
      precenv.decomp();
      all_precenv.push_back(precenv);
      }
    }

  if (!nofixed && (type==RW1) && (!varcoeff) )
    {
    hierarchie_rw1(lvec,1);
    }
  else if (!nofixed && (type==RW1) && (varcoeff) )
    {
    if(identifiable)
      {
      hierarchie_rw1(lvec,2);
      lvec.push_back(-1);
      }
    else
      hierarchie_rw1(lvec,1);
    }
  else if (!nofixed && (type==RW2) && (!varcoeff) )
    {
    lvec.push_back(-1);
    }
  else if (!nofixed && (type==RW2) && (varcoeff) )
    {
    lvec.push_back(-2);
    if(identifiable)    //VCM_neu
      lvec.push_back(-1);
    }
  else if ( (type==mrf) && (varcoeff) )
    {
    if(!nofixed && identifiable)
      lvec.push_back(-1);
    }

  if(forced_into==false)
     lvec.push_back(0);

  // Startwert für lambda aus df:
  if(spfromdf!="direct")
    {
    double lambdavorg = 1000;
    if(!varcoeff)
      {
      if(!nofixed && dfstart==1 && (type==RW1 || type==RW2))
        lambdastart = -1;
      else if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    else
      {
      if(!nofixed && dfstart==1 && identifiable)
        lambdastart = -1;
      else if(!nofixed && ((dfstart==2 && identifiable) || (dfstart==1 && !identifiable)) && (type==RW1 || type==RW2))
        lambdastart = -2;
      else if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    if(lambdastart==-9 || lambdastart==1000000000)    // falls dfstart nicht erreicht werden kann
      lambdastart = 0;
    }

  if(kchanged == true)
    kombimatrix = true;
  }


// BEGIN: For Varying Coefficients ---------------------------------------------

void FULLCOND_nonp_gaussian_stepwise::update_fix_effect(double & intercept)
  {
  bool raus = false;
  unsigned j = 1;
  ST::string name_richtig = datanames[1];  // (VC alt) 0
  while(j<fcconst->get_datanames().size() && raus==false)
     {
     if(fcconst->get_datanames()[j] == datanames[1])   // (VC alt) 0
        {
        raus = true;
        }
     if(fcconst->get_datanames()[j] == (datanames[1]+"_1"))  // (VC alt) 0
        {
        raus = true;
        name_richtig = datanames[1] + "_1";            // (VC alt) 0
        }
     j = j + 1;
     }
  if(raus == true)
    {
    fcconst->update_fix_effect(j-1,intercept,data_forfixed);
    }
  else
    {
    vector<ST::string> names;
    names.push_back(name_richtig);
    fcconst->include_effect(names,data_forfixed);
    interactions_pointer[0]->set_inthemodel(-1);
    fcconst->update_fix_effect(j,intercept,data_forfixed);
    }
  }


void FULLCOND_nonp_gaussian_stepwise::set_pointer_to_interaction(FULLCOND * inter)
  {
  interactions_pointer.push_back(inter);
  }

void FULLCOND_nonp_gaussian_stepwise::get_interactionspointer(vector<FULLCOND*> & inter)
  {
  inter = interactions_pointer;
  }


void FULLCOND_nonp_gaussian_stepwise::hierarchical(ST::string & possible)
  {
  if(!varcoeff)
    {
    possible = "alles";
    }
  else
    {
    possible = "valles";
    }
  }


void FULLCOND_nonp_gaussian_stepwise::const_varcoeff(void)
  {
  if(varcoeff)
    fcconst->posteriormode_const_varcoeff(data_forfixed);
  }

// END: For Varying Coefficients -----------------------------------------------


void FULLCOND_nonp_gaussian_stepwise::create_weight(datamatrix & w)
  {
  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
  int * workindex = index.getV();
  unsigned i;

  if(type != mrf)
    {
    w(*workindex,0) = 1;

    itbeg = posbeg.begin() + posbeg.size()-1;
    workindex = index.getV() + *itbeg;
    w(*workindex,0) = 1;
    }
  else
    {
    int j;
    for(i=0;i<nrpar;i++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        w(*workindex,0) = 1;
        for(j=*itbeg;j<=*itend;j++)
          workindex++;
        }
      }
    }
  }


// BEGIN: MODEL-AVERAGING ------------------------------------------------------

void FULLCOND_nonp_gaussian_stepwise::update_bootstrap(const bool & uncond)
  {
if(kombimatrix == false || matrixnumber==1)
  {
  update_bootstrap_df();

  datamatrix betaold = beta;

  if(fixornot==true)
    {
    bool raus = false;
    unsigned j = 1;
    ST::string name_richtig;
    if(!varcoeff)
      name_richtig = datanames[0];
    else
      name_richtig = datanames[1];       // (VC alt) 0
    while(j<fcconst->get_datanames().size() && raus==false)
      {
      if(fcconst->get_datanames()[j] == name_richtig)
        raus = true;
      j = j + 1;
      }
    unsigned index_fix = j-1;
    double fix = fcconst->getbeta(index_fix,0);
    unsigned i;
    double * workbeta = beta.getV();
    vector<int>::iterator itbeg = posbeg.begin();
    vector<int>::iterator itend = posend.begin();
    int * workindex = index.getV();
    double sum = 0;
    int k;
    for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        if(!varcoeff)
          {
          *workbeta = data_forfixed(*workindex,0)*fix;
          sum += *workbeta;
          }
        else
          *workbeta = fix;
        for(k=*itbeg;k<=*itend;k++)
          workindex++;
        }
      }
    workbeta = beta.getV();
    sum /= double(nrpar);
    if(!center)
      sum = 0;
    for(i=0;i<nrpar;i++,workbeta++)
      *workbeta -= sum;
    double chelp = sum*double(nrpar);
    fcconst->update_intercept(chelp);

    FULLCOND::update_bootstrap();
    //fc_df.setbetavalue(0,0,1);
    }
  else if(inthemodel==false && fixornot==false)
    {
    beta = datamatrix(nrpar,1,0);
    FULLCOND::update_bootstrap();
    }
  else
    {
    FULLCOND::update_bootstrap();
    }
  beta = betaold;
  }
  }


void FULLCOND_nonp_gaussian_stepwise::update_bootstrap_df(void)
  {
if(kombimatrix == false || matrixnumber==1)
  {

  if(optionsp->get_nriter()<=1)
    {
    ST::string path = samplepath.substr(0,samplepath.length()-4)+"_df.raw";
    fc_df = FULLCOND(optionsp,datamatrix(1,1),"title?",1,1,path);
    fc_df.setflags(MCMC::norelchange | MCMC::nooutput);
    isbootstrap = true;
    }

  if(fixornot==true)
    {
    fc_df.setbetavalue(0,0,-1.0);
    fc_df.update_bootstrap_df();
    }
  else if(inthemodel==false && fixornot==false)
    {
    fc_df.setbetavalue(0,0,0.0);
    fc_df.update_bootstrap_df();
    }
  else // if(inthemodel==true && fixornot==false)
    {
    fc_df.setbetavalue(0,0,lambda);
    fc_df.update_bootstrap_df();
    }

  betaold = datamatrix(1,1,0); // löscht "beta" vom letzten Bootstrap-Sample
  }
  }


void FULLCOND_nonp_gaussian_stepwise::save_betamean(void)
  {
if(kombimatrix == false || matrixnumber==1)
  {
  datamatrix betaold = beta;

  if(fixornot==true)
    {
    bool raus = false;
    unsigned j = 1;
    ST::string name_richtig;
    if(!varcoeff)
      name_richtig = datanames[0];
    else
      name_richtig = datanames[1];       // (VC alt) 0
    while(j<fcconst->get_datanames().size() && raus==false)
      {
      if(fcconst->get_datanames()[j] == name_richtig)
        raus = true;
      j = j + 1;
      }
    unsigned index_fix = j-1;
    double fix = fcconst->getbeta(index_fix,0);
    unsigned i;
    double * workbeta = beta.getV();
    vector<int>::iterator itbeg = posbeg.begin();
    vector<int>::iterator itend = posend.begin();
    int * workindex = index.getV();
    double sum = 0;
    int k;
    for(i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
      {
      if(*itbeg != -1)
        {
        if(!varcoeff)
          {
          *workbeta = data_forfixed(*workindex,0)*fix;
          sum += *workbeta;
          }
        else
          *workbeta = fix;
        for(k=*itbeg;k<=*itend;k++)
          workindex++;
        }
      }
    workbeta = beta.getV();
    sum /= double(nrpar);
    if(!center)
      sum = 0;
    for(i=0;i<nrpar;i++,workbeta++)
      *workbeta -= sum;
    double chelp = sum*double(nrpar);
    fcconst->update_intercept(chelp);

    FULLCOND::save_betamean();
    //fc_df.setbetavalue(0,0,1);
    }
  else if(inthemodel==false && fixornot==false)
    {
    beta = datamatrix(nrpar,1,0);
    FULLCOND::save_betamean();
    }
  else
    {
    FULLCOND::save_betamean();
    }
  beta = betaold;
  }
  }


void FULLCOND_nonp_gaussian_stepwise::update_bootstrap_betamean(void)
  {
if(kombimatrix == false || matrixnumber==1)
  {
  FULLCOND::update_bootstrap_betamean();
  FULLCOND::setflags(MCMC::norelchange);
  }
  }


void FULLCOND_nonp_gaussian_stepwise::outresults_df(unsigned & size)
  {
if(kombimatrix == false || matrixnumber==1)
  {

  fc_df.update_bootstrap_betamean();
  //fc_df.outresults();
  double * workmean = fc_df.get_betameanp();

  ST::string pathdf = pathcurrent.substr(0,pathcurrent.length()-4)+"_df.res";
  ofstream outres(pathdf.strtochar());

  outres << "df_value   ";
  outres << "sp_value  ";
  outres << "frequency  ";
  outres << "selected  " << endl;

// Häufigkeitstabelle:

  //samplestream.close();
  datamatrix sample(size,1);
  fc_df.readsample_df(sample,0);
  unsigned i;

  vector<unsigned> number;
  vector<unsigned> number1;
  vector<unsigned> number2;
  vector<unsigned> cumnumber1;
  vector<unsigned> cumnumber;

  statmatrix<int> index(sample.rows(),1);
  index.indexinit();
  sample.indexsort(index,0,sample.rows()-1,0,0);

  i = 0;
  unsigned j,anz;
  while(i<index.rows())
     {
     anz=0;
     int* p = index.getV() + i;
     int* q = index.getV() + i;
     j=i;
     while(j<index.rows() && (sample.get(*p,0) == sample.get(*q,0)))
        {
        anz = anz+1;
        j++;
        p++;
        }
     if(sample.get(*q,0) <= 0)
       number1.push_back(anz);
     else if(sample.get(*q,0) > 0)
       number2.push_back(anz);
     if(cumnumber1.size()>0)
       cumnumber1.push_back(cumnumber1[cumnumber1.size()-1]+anz);
     else
       cumnumber1.push_back(anz);
     i = i + anz;
     }

  int k;
  for(k=number1.size()-1;k>=0;k--)
    {
    cumnumber.push_back(cumnumber1[k]);
    number.push_back(number1[k]);
    }
  for(k=number2.size()-1;k>=0;k--)
    {
    cumnumber.push_back(cumnumber1[k+number1.size()]);
    number.push_back(number2[k]);
    }

  double merken = 0;
  for(i=0;i<number.size();i++)
    {
    double help = sample.get(index(cumnumber[i]-1,0),0);
    double dfs = -1*help;
    if(help!=0 && help!=-1)
      {
      if(merken==0 && help!=-2)
        merken = i;
      update_stepwise(help);
      set_inthemodel(help);
      dfs = compute_df();
      }
    outres << ST::doubletostring(dfs,6) << "   " << ST::doubletostring(help,6) << "   " << ST::inttostring(number[i]) << "   ";
    if(*workmean == help)
      outres << "+"; // ST::doubletostring(*workmean,6);
    else
      outres << "-";
    outres << endl;
    }
  }
  }

void FULLCOND_nonp_gaussian_stepwise::update(void)
  {
  if(kombimatrix == false || matrixnumber==1)
  {
  if (betaold.rows() == 1)
    {
    betaold = beta;
    betaKbeta=Kenv.compute_quadform(beta,0);
    }

  if(lambda==0)
    {
    beta = datamatrix(beta.rows(),beta.cols(),0);
    FULLCOND::update();
    }
  else
    {
    if(utype == gaussian)
      update_gauss();
    else
      update_IWLS();
    }
  }
  }


void FULLCOND_nonp_gaussian_stepwise::update_gauss(void)
  {
  int j;
  unsigned i;

  int * workindex;
  update_linpred(false);

  if(optionsp->get_nriter()==1 || changingweight)
    {
    if(varcoeff)
      compute_XWX_varcoeff_env(likep->get_weight());
    else
      compute_XWX_env(likep->get_weight());
    }

  precenv.addtodiag(XXenv,Kenv,1.0,lambda);
  if(kombimatrix==true)
    precenv.addto(precenv,Kenv2,1.0,otherfullcond[0]->get_lambda());
    if(numberofmatrices==3)
      precenv.addto(precenv,Kenv3,1.0,otherfullcond[1]->get_lambda());

  double sigmaresp = sqrt(likep->get_scale(column));

  double * work = betahelp.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = sigmaresp*rand_normal();

  precenv.solveU(betahelp);

  likep->compute_respminuslinpred(mu,column);

  workindex = index.getV();
  double * workmuy = muy.getV();

  if (varcoeff)
    {
    double * workdata=data.getV();
    for(i=0;i<nrpar;i++,workmuy++)
      {
      *workmuy = 0;
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++,workdata++)
          *workmuy+= likep->get_weight(*workindex,0)*mu(*workindex,0)*
          (*workdata);
      }
    }
  else  // else additive
    {
    for(i=0;i<nrpar;i++,workmuy++)
      {
      *workmuy = 0;
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++)
          *workmuy+= likep->get_weight(*workindex,0)*mu(*workindex,0);
      }
    }

  precenv.solve(muy,betahelp,beta);

  update_linpred(true);

  if (center)
    {
    double m = centerbeta();
    if (varcoeff)
      fcconst->update_fix_varcoeff(m,datanames[1]);
    else
      fcconst->update_intercept(m);
    }

  acceptance++;

  transform = likep->get_trmult(column);

  FULLCOND::update();
  }


void FULLCOND_nonp_gaussian_stepwise::update_IWLS(void)
  {
  unsigned i;

  if(betaold.rows()!=beta.rows())
    betaold = datamatrix(beta.rows(),1,0);

  betaold.assign(beta);

  double sigma2 = likep->get_scale(column)/lambda;;
  double invscale = 1.0/likep->get_scale(column);

  envmatdouble Ksum;
  if(kombimatrix==true)
    {
    if(numberofmatrices==2)
      {
      Ksum = envmatdouble(0.0,nrpar,Kenv.getBandwidth());
      Ksum.addto(Kenv,Kenv2,1.0/sigma2,otherfullcond[0]->get_lambda()*invscale);
      betaKbeta = Ksum.compute_quadform(beta,0);
      }
    else if(numberofmatrices==3)
      {
      Ksum = envmatdouble(0.0,nrpar,Kenv2.getBandwidth());
      Ksum.addto(Kenv3,Kenv,otherfullcond[1]->get_lambda()*invscale,1.0/sigma2);
      Ksum.addto(Ksum,Kenv2,1.0,otherfullcond[0]->get_lambda()*invscale);
      betaKbeta = Ksum.compute_quadform(beta,0);
      }
    }
  else
    betaKbeta = Kenv.compute_quadform(beta,0);

  double * workbeta;

  // Compute log-likelihood with old beta

  double logold;
  if (  (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) )
     )
    {
    logold = likep->compute_IWLS(weightiwls,tildey,true,column);
    }
  else
    {
    logold = likep->compute_IWLS(weightiwls,tildey,false,column);
    }

  if(kombimatrix==true)
    logold -= 0.5*betaKbeta;
  else
    logold -= 0.5*betaKbeta/sigma2;

  workbeta = betaold.getV();
  if (  (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) )
      )
    {
    if (varcoeff)
      compute_XWX_XWtildey_varcoeff_env(weightiwls,tildey,workbeta,0);
    else
      compute_XWX_XWtildey_env(weightiwls,tildey,workbeta,0);

    precenv.addtodiag(XXenv,Kenv,invscale,1.0/sigma2);
    if(kombimatrix==true)
      {
      precenv.addto(precenv,Kenv2,1.0,otherfullcond[0]->get_lambda()*invscale);
      }
    }
  else
    {
    compute_muy(workbeta);
    }

  double * workmuy = muy.getV();
  if(invscale != 1.0)
    {
    for(i=0;i<likep->get_nrobs();i++,workmuy++)
      *workmuy = *workmuy*invscale;
    }

  precenv.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  precenv.solveU(beta,betahelp);

  betahelp.minus(beta,betahelp);

  double qold = 0.5*precenv.getLogDet()- 0.5*precenv.compute_quadform(betahelp,0);

  update_linpred_diff(beta,betaold);

  // Proposal computed and stored in beta, linear predictor with new beta
  // Compute new log-likelihood

  double lognew;
  if (  (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) )
     )
    {
    lognew = likep->compute_IWLS(weightiwls,tildey,true,column);
    }
  else
    {
    lognew = likep->compute_IWLS(weightiwls,tildey,false,column);
    }

  if(kombimatrix==true)
    lognew  -= 0.5*Ksum.compute_quadform(beta,0);
  else
    lognew  -= 0.5*Kenv.compute_quadform(beta,0)/sigma2;

  workbeta = beta.getV();
  if (  (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) )
      )
    {
    if (varcoeff)
      compute_XWX_XWtildey_varcoeff_env(weightiwls,tildey,workbeta,0);
    else
      compute_XWX_XWtildey_env(weightiwls,tildey,workbeta,0);

    precenv.addtodiag(XXenv,Kenv,invscale,1.0/sigma2);
    if(kombimatrix==true)
      {
      precenv.addto(precenv,Kenv2,1.0,otherfullcond[0]->get_lambda()*invscale);
      }
    }
  else
    {
    compute_muy(workbeta);
    }

  workmuy = muy.getV();
  if(invscale != 1.0)
    {
    for(i=0;i<likep->get_nrobs();i++,workmuy++)
      *workmuy = *workmuy*invscale;
    }

  precenv.solve(muy,betahelp);

  betahelp.minus(betaold,betahelp);
  double qnew = 0.5*precenv.getLogDet() - 0.5*precenv.compute_quadform(betahelp,0);

  double u = log(uniform());
  if (u <= (lognew - logold  + qnew - qold) )
    {
    acceptance++;
    if (center)
      {
      double m = centerbeta();
      if (varcoeff)
        fcconst->update_fix_varcoeff(m,datanames[1]);
      else
        fcconst->update_intercept(m);
      }
    betaold.assign(beta);

    if(!adaptiv)
      betaKbeta=Kenv.compute_quadform(beta,0);
    }
  else
    {
    update_linpred_diff(betaold,beta);
    beta.assign(betaold);
    }

  transform = likep->get_trmult(column);
  FULLCOND::update();
  }


void FULLCOND_nonp_gaussian_stepwise::change_Korder(double lamb)
  {
  set_lambdaconst(1000000000);
  if(lamb==-1)
    {
    if(!varcoeff)
      {
      if(type==RW1)
        {
        Kenv = Krw2env(weight);
        rankK = nrpar-2;
        }
      }
    else
      {
      if(type==RW2)
        {
        Kenv = Krw1env(weight);
        rankK = nrpar-1;
        }
      }
    }
  else if(lamb==-2)
    {
    if(!varcoeff)
      {
      // MONOTON: fehlt noch!!!
      }
    else
      {
      if(type==RW1)
        {
        Kenv = Krw2env(weight);
        rankK = nrpar-2;
        }
      }
    }
  }

void FULLCOND_nonp_gaussian_stepwise::undo_Korder(void)
  {
  if(type==RW1 && rankK==nrpar-2)
    {
    Kenv = Krw1env(weight);
    rankK = nrpar-1;
    }
  else if(type==RW2 && rankK==nrpar-1)
    {
    Kenv = Krw2env(weight);
    rankK = nrpar-2;
    }
  }


void FULLCOND_nonp_gaussian_stepwise::outresults(void)
  {
if(matrixnumber == 1 || kombimatrix == false)
  {

  FULLCOND_nonp_gaussian::outresults();

  /*
  if(!isbootstrap)
    {
    FULLCOND_nonp_gaussian::outresults();
    }
  else
    {
    if (lattice==false)
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
        #else
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
      outres << "paverage   ";
      outres << "pqu"  << l1  << "   ";
      outres << "pqu"  << l2  << "   ";
      outres << "pmed   ";
      outres << "pqu"  << u1  << "   ";
      outres << "pqu"  << u2  << "   ";
      outres << "pcat" << level1 << "   ";
      outres << "pcat" << level2 << "   ";
      outres << endl;

      double * workmean = betamean.getV();
      double * workave = beta_average.getV();
      double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
      double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
      double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
      double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
      double * workbetaqu50 = betaqu50.getV();

      vector<ST::string>::iterator effit = effectvalues.begin();

      for(i=0;i<nrpar;i++,++effit,workmean++,workave++,workbetaqu_l1_lower_p++,
                            workbetaqu_l2_lower_p++,workbetaqu50++,
                            workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
        {
        outres << (i+1) << "   ";
        outres << *effit << "   ";
        outres << *workmean << "   ";
        outres << *workave << "   ";
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
    else
      {
      FULLCOND::outresults();

      optionsp->out("  Results are stored in file " + pathresults + "\n");
      optionsp->out("  Corresponding boundary-file is stored in " + pathmap + "\n");
      ST::string pathgraph = pathmap.substr(0, pathmap.length()-4);
      pathgraph = pathgraph + "gra";
      optionsp->out("  Corresponding graph-file is stored in " + pathgraph + "\n");
      #if defined(JAVA_OUTPUT_WINDOW)
      optionsp->out("  Results may be visualized using method 'drawmap'\n");
      optionsp->out("  Type for example: objectname.drawmap " +
      ST::inttostring(fcnumber) + "\n");
      #else
      optionsp->out("  Results may be visualized using the R function 'drawmap'\n");
      #endif
      optionsp->out("\n");

      unsigned i;

      ofstream outres(pathresults.strtochar());

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
      outres << "xcoord" << "   ";
      outres << "ycoord" << "   ";

      outres << name << "   ";
      outres << "pmean   ";
      outres << "paverage   ";
      outres << "pqu"  << l1  << "   ";
      outres << "pqu"  << l2  << "   ";
      outres << "pmed   ";
      outres << "pqu"  << u1  << "   ";
      outres << "pqu"  << u2  << "   ";
      outres << "pcat" << level1 << "   ";
      outres << "pcat" << level2 << "   ";

      outres << endl;

      double * workmean = betamean.getV();
      double * workave = beta_average.getV();
      double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
      double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
      double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
      double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
      double * workbetaqu50 = betaqu50.getV();
      double * workxyvalues = xyvalues.getV();

      vector<ST::string>::iterator effit = effectvalues.begin();

      for(i=0;i<nrpar;i++,++effit,workmean++,workave++,workbetaqu_l1_lower_p++,
                             workbetaqu_l2_lower_p++,workbetaqu50++,
                             workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                             workxyvalues++)
        {
        outres << (i+1) << "   ";
        outres << *workxyvalues << "   ";
        workxyvalues++;
        outres << *workxyvalues << "   ";
        outres << *effit << "   ";
        outres << *workmean << "   ";
        outres << *workave << "   ";
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
    }
    */
  }
  }



/*     in envmatrix_penalty:
envmatrix<double> Krw0env(const unsigned & nrpar)
  {
  vector<double> diag(nrpar,1);
  vector<double> env;
  vector<unsigned> xenv(nrpar+1,0);

  return envmatrix<double>(env, diag, xenv, 0);
  }*/


bool FULLCOND_nonp_gaussian_stepwise::posteriormode_kombi(void)
  {
if(matrixnumber>1)
  return otherfullcond[0]->posteriormode();
else
  {
  unsigned i;
  int j;
  int * workindex;

  update_linpred(false);

  kappa.erase(kappa.begin(),kappa.end());
  kappa.push_back(otherfullcond[0]->get_lambda());
  if(numberofmatrices==3)
      kappa.push_back(otherfullcond[1]->get_lambda());
  if ( (lambda_prec != lambda) || (kappa_prec != kappa) || (calculate_xwx == true))
    {
    if(calculate_xwx == true)
      {
      calculate_xwx = false;
      if (varcoeff)
        compute_XWX_varcoeff_env(likep->get_weightiwls(),column);
      else
        compute_XWX_env(likep->get_weightiwls(),column);
      }
    precenv.addto(XXenv,Kenv,1.0,lambda);       // hier gehört lambda zu "MRF" und kappa zu "I"
    precenv.addto(precenv,Kenv2,1.0,kappa[0]);
    if(type == twomrfI)
      precenv.addto(precenv,Kenv3,1.0,kappa[1]);
//precenv.addto(XXenv,Kenv,1.0,kappa);
//precenv.addto(precenv,Kenv2,1.0,lambda);
    lambda_prec = lambda;
    kappa_prec = kappa;
    }

  likep->compute_weightiwls_workingresiduals(column);

  workindex = index.getV();
  double * workmuy = beta.getV();

  if (varcoeff)
    {
    double * workdata=data.getV();
    for(i=0;i<nrpar;i++,workmuy++)
      {
      *workmuy = 0;
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++,workdata++)
          *workmuy+=
          likep->get_workingresiduals()(*workindex,0)*(*workdata);
      }
    }
  else  // else additive
    {
    for(i=0;i<nrpar;i++,workmuy++)
      {
      *workmuy = 0;
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++)
          *workmuy+= likep->get_workingresiduals()(*workindex,0);
      }
    }

  precenv.solve(beta);
  update_linpred(true);

  if (center)
    {
    intercept = centerbeta();
    if(varcoeff == false)
      {
      fcconst->posteriormode_intercept(intercept);
      }
    else if(varcoeff == true)
      {
      update_fix_effect(intercept);
      intercept = 0.0;
      }
    }

  transform = likep->get_trmult(column);
  return FULLCOND_nonp_basis::posteriormode();
  }
  }

double FULLCOND_nonp_gaussian_stepwise::compute_df_kombi(void)
  {
  double df = 0;
if(matrixnumber==1)
  {
  kappa.erase(kappa.begin(),kappa.end());
  kappa.push_back(otherfullcond[0]->get_lambda());
  if(type == twomrfI)
    kappa.push_back(otherfullcond[1]->get_lambda());
  if(inthemodel == true)
    {
    if (lambdaold == lambda && kappaold == kappa && likep->get_iwlsweights_notchanged() == true)
      {
      df = df_lambdaold;
      }
    else if (lambdaold != lambda || kappa != kappaold || likep->get_iwlsweights_notchanged() == false)
      {
      if (calculate_xwx == true)
        {
        if (varcoeff)
          compute_XWX_varcoeff_env(likep->get_weightiwls(),column);
        else
          compute_XWX_env(likep->get_weightiwls(),column);
        calculate_xwx = false;
        }
      datamatrix Z = datamatrix(nrpar,nrpar,0);
      vector<double>::iterator d = XXenv.getDiagIterator();
      vector<double>::iterator d2 = XXenv.getDiagIterator();
      vector<double>::iterator p = Kenv.getDiagIterator();
      vector<double>::iterator p2 = Kenv2.getDiagIterator();
      vector<double>::iterator p3;
      if(type == twomrfI)
        p3 = Kenv3.getDiagIterator();
      double sumn = 0; // enthält Summe der Elemente von XXenv
      unsigned i,j;

      for(i=0;i<nrpar;i++,d++)
        {
        sumn += *d;
        }

      d = XXenv.getDiagIterator();
      double * workz = Z.getV();
      if(type != twomrfI)
        {
        for(i=0;i<nrpar;i++,d++,p++,p2++,workz++)
          {
          *workz = *d - 1/sumn * *d * *d + lambda * *p + kappa[0] * *p2;
          d2 = XXenv.getDiagIterator() + i+1;
          workz++;
          for(j=i+1;j<nrpar;j++,workz++,d2++)
            {
            *workz = -1/sumn * *d * *d2;
            }
          workz += i;
          }
        }
      else
        {
        for(i=0;i<nrpar;i++,d++,p++,p2++,p3++,workz++)
          {
          *workz = *d - 1/sumn * *d * *d + lambda * *p + kappa[0] * *p2 + kappa[1] * *p3;
          d2 = XXenv.getDiagIterator() + i+1;
          workz++;
          for(j=i+1;j<nrpar;j++,workz++,d2++)
            {
            *workz = -1/sumn * *d * *d2;
            }
          workz += i;
          }

        }


      int kl,ku,k;
      vector<unsigned>::iterator a;
      if(type != twomrfI)
        {
        for(i=0;i<nrpar;i++)
          {
          a = Kenv.getXenvIterator() + i;
          kl = *a;
//          kl=xenv[i];
          a = Kenv.getXenvIterator() + i + 1;
          ku = *a;
//          ku=xenv[i+1];
          for(k=kl; k<ku; k++)
            {
            d = Kenv.getEnvIterator() + k;
            Z(i-ku+k,i) = Z(i-ku+k,i) + lambda * *d;
          //Z(i,i-ku+k) = Z(i-ku+k,i);
            }
          }
        }
      else
        {
        vector<double>::iterator b;
        for(i=0;i<nrpar;i++)
          {
          a = Kenv.getXenvIterator() + i;
          kl = *a;
//          kl=xenv[i];
          a = Kenv.getXenvIterator() + i + 1;
          ku = *a;
//          ku=xenv[i+1];
          for(k=kl; k<ku; k++)
            {
            d = Kenv.getEnvIterator() + k;
            b = Kenv2.getEnvIterator() + k;
            Z(i-ku+k,i) = Z(i-ku+k,i) + lambda * *d + kappa[0] * *b;
          //Z(i,i-ku+k) = Z(i-ku+k,i);
            }
          }
        }

      for(i=0;i<nrpar;i++)    // schreibt untere Hälfte von Z voll
        {
        for(j=i+1;j<nrpar;j++)
          {
          Z(j,i) = Z(i,j);
          }
        }

//ofstream out("c:\\cprog\\test\\results\\Z.txt");
//Z.prettyPrint(out);

      Z.assign(Z.inverse());    // berechnet Inverse von (n1-n1^2/n, -n1*n2/n,..., -n1*np/n) + lambda*P

      workz = Z.getV();
      d = XXenv.getDiagIterator();
      df = 1;
      for(i=0;i<nrpar;i++,d++,workz++)
        {
        df += *d * *workz - 1/sumn * *d * *d * *workz;
        p = XXenv.getDiagIterator() + i+1;
        workz++;
        for(j=i+1;j<nrpar;j++,p++,workz++)
          {
          df -= 2/sumn * *d * *p * *workz;
          }
        workz += i;
        }
      df -= 1;

      df_lambdaold = df;
      lambdaold = lambda;
      kappaold = kappa;
      }
    }

  }
  return df;
  }


} // end: namespace MCMC




