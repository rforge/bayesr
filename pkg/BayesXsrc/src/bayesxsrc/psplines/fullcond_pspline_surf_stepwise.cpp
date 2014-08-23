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



#include "fullcond_pspline_surf_stepwise.h"

namespace MCMC
{


void FULLCOND_pspline_surf_stepwise::init_maineffects(
                        FULLCOND_pspline_stepwise * mp1,FULLCOND_pspline_stepwise * mp2,
                        const ST::string & pnt,const ST::string & prt)
  {

  mainpoi1 = mp1;
  mainpoi2 = mp2;

  interactions_pointer.push_back(mp1);
  interactions_pointer.push_back(mp2);

  interaction = true;

  centertotal = false;
  maineffectsexisting = 11;

  fctotalrespath = prt;

  datamatrix h(1,1,0);
  if(gridsize < 0)
    fctotal = FULLCOND(optionsp,h,title+"total",nrdiffobs,1,pnt);
  else
    fctotal = FULLCOND(optionsp,h,title+"total",gridsize,1,pnt);
  fctotal.setflags(MCMC::norelchange | MCMC::nooutput);
  fctotal.set_transform(transform);

  beta1 = datamatrix(nrpar1dim,1,0);
  beta2 = datamatrix(nrpar1dim,1,0);
  he1 = datamatrix(xv.size(),1,0);
  he2 = datamatrix(yv.size(),1,0);

  // 2) Berechnen der Strafmatrix für 1. Haupteffekt
  if(type == mrfkr1)
    {
    datamatrix K1 = datamatrix(nrpar1dim,nrpar1dim,0.0);
    unsigned i,j,bands;
    datamatrix de(nrpar,1);
    datamatrix ud;

    for(i=0;i<nrpar1dim;i++)
      {
      if(mainpoi2->get_type() == MCMC::RW1)
        K1(i,i) = 1.0;
      else
        K1(i,i) = 2.0;
      }
    if(mainpoi2->get_type() == MCMC::RW1)
      {
      bands=1;
      Ksp = Krw1(vector<double>(nrpar1dim,1.0)).kronecker(SparseMatrix(K1,true));
      ud = datamatrix(nrpar,(nrpar1dim)*bands,0);
      }
    else
      {
      bands=2;
      Ksp = Krw2(vector<double>(nrpar1dim,1.0)).kronecker(SparseMatrix(K1,true));
      ud = datamatrix(nrpar,(nrpar1dim)*bands,0);
      }

    for(i=0;i<nrpar;i++)
      {
      de(i,0) = Ksp(i,i);
      for(j=0;j<ud.cols();j++)
        {
        if (i+j+1 < nrpar)
          ud(i,j) = Ksp(i,i+j+1);
        }
      } // end: for(i=0;i<sizeK;i++)

    Ky = bandmatdouble(de,ud);
    Kyenv = envmatdouble(Ky);

  // 3) Berechnen der Strafmatrix für 2. Haupteffekt

    for(i=0;i<nrpar1dim;i++)
      {
      if(mainpoi1->get_type() == MCMC::RW1)
        K1(i,i) = 1.0;
      else
        K1(i,i) = 2.0;
      }
    if(mainpoi1->get_type() == MCMC::RW1)
      {
      bands=1;
      Ksp = SparseMatrix(K1,true).kronecker(Krw1(vector<double>(nrpar1dim,1.0)));
      ud = datamatrix(nrpar,(nrpar1dim)*bands,0);
      }
    else
      {
      bands=2;
      Ksp = SparseMatrix(K1,true).kronecker(Krw2(vector<double>(nrpar1dim,1.0)));
      ud = datamatrix(nrpar,(nrpar1dim)*bands,0);
      }

    for(i=0;i<nrpar;i++)
      {
      de(i,0) = Ksp(i,i);
      for(j=0;j<ud.cols();j++)
        {
        if (i+j+1 < nrpar)
          ud(i,j) = Ksp(i,i+j+1);
        }
      } // end: for(i=0;i<sizeK;i++)

    Kx = bandmatdouble(de,ud);
    Kxenv = envmatdouble(Kx);

    double bandw = Kenv.getBandwidth();
    KHenv = envmatdouble(0.0,nrpar,bandw);
    KHenv = Kenv;

    if(Kxenv.getBandwidth() > bandw)
      bandw = Kxenv.getBandwidth();
    if(Kyenv.getBandwidth() > bandw)
      bandw = Kyenv.getBandwidth();
    Kenv = envmatdouble(0.0,nrpar,bandw);
    }

/*ofstream out("c:\\cprog\\test\\results\\Kenv.txt");
Kenv.print4(out);
ofstream out2("c:\\cprog\\test\\results\\Kxenv.txt");
Kxenv.print4(out2);
ofstream out3("c:\\cprog\\test\\results\\Kyenv.txt");
Kyenv.print4(out3); */

  }


void FULLCOND_pspline_surf_stepwise::create(const datamatrix & v1, const datamatrix & v2, const datamatrix & intact)
  {

  unsigned i;
  data_forfixed = datamatrix(v1.rows(),1,0);

  // neu:
  for(i=0;i<v2.rows();i++)
    {
    data_forfixed(i,0) = v1(i,0) * v2(i,0);
    }

  }

  // CONSTRUCTOR

FULLCOND_pspline_surf_stepwise::FULLCOND_pspline_surf_stepwise(MCMCoptions * o,
                      DISTRIBUTION * dp, FULLCOND_const * fcc,
                      const datamatrix & v1, const datamatrix & v2,
                      const ST::string & ti,
                      const unsigned & nrk, const unsigned & degr,
                      const knotpos & kp, const double & l, const int & gs,
                      const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const ST::string & of, const bool & gauss, const unsigned & c)
  : FULLCOND_pspline_surf_gaussian(o,dp,fcc,v1,v2,ti,nrk,degr,kp,l,gs,ft,fp,pres,of,true,c)
  {

  fctype = nonparametric;

  if(gauss==false)
    utype = iwls;

  create(v1,v2);

  maineffectsexisting = 0;

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  lambdaxold = -1;
  lambdayold = -1;
  lambdax_prec = -1;
  lambday_prec = -1;
  lambdaold = -1;
  }


  // CONSTRUCTOR 2: varying coefficients

FULLCOND_pspline_surf_stepwise::FULLCOND_pspline_surf_stepwise(MCMCoptions * o,
                      DISTRIBUTION * dp, FULLCOND_const * fcc,
                      const datamatrix & intact,
                      const datamatrix & v1, const datamatrix & v2,
                      const ST::string & ti, const unsigned & nrk, const unsigned & degr,
                      const knotpos & kp, const double & l, const int & gs,
                      const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const ST::string & of, const bool & sb, const bool & vccent, const unsigned & c)
  : FULLCOND_pspline_surf_gaussian(o,dp,fcc,intact,v1,v2,ti,nrk,degr,kp,l,gs,ft,fp,pres,of,sb,c)
  {

  fctype = nonparametric;

  centertotal = true;
  //VCM_neu
  if(vccent == true)
    identifiable = false;

  data_forfixed = intact;

  unsigned nrobs = intact.rows();

  if (data_varcoeff_fix.rows() < nrobs)
    {
    data_varcoeff_fix = datamatrix(nrobs,2,1);
    effmodi = datamatrix(nrobs,1,0);
    for(unsigned i=0;i<nrobs;i++)
      {
      effmodi(i,0) = v1(i,0)*v2(i,0);
      data_varcoeff_fix(i,0) = intact(i,0);
      data_varcoeff_fix(i,1) = effmodi(i,0)*intact(i,0);
      }
    }

  centervalue = (v1.max(0) + v1.min(0)) * (v2.max(0) + v2.min(0));
  XVX = datamatrix(2,2,0);

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  lambdaold = -1;
  }

  // CONSTRUCTOR 3: geosplines

FULLCOND_pspline_surf_stepwise::FULLCOND_pspline_surf_stepwise(
                      MCMCoptions * o, DISTRIBUTION * dp,
                      FULLCOND_const * fcc,
                      const datamatrix & region, const MAP::map & mp, const ST::string & mn,
                      const ST::string & ti, const unsigned & nrk,
                      const unsigned & degr, const knotpos & kp, const double & l,
                      const int & gs, const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const bool & gauss, const unsigned & c)
  : FULLCOND_pspline_surf_gaussian(o,dp,fcc,region,mp,mn,ti,nrk,degr,kp,l,gs,ft,fp,pres,true,c)
  {

  fctype = nonparametric;

  if(gauss==false)
    utype = iwls;

  datamatrix v1 = datamatrix(likep->get_nrobs(),1,0.0);
  datamatrix v2 = datamatrix(likep->get_nrobs(),1,0.0);

  ST::string regname;
  for(unsigned i=0;i<likep->get_nrobs();i++)
    {
    regname = ST::doubletostring(region(i,0));
    regionnames.push_back(regname);
    v1(i,0) = m.get_region(m.getnr(regname)).get_xcenter();
    v2(i,0) = m.get_region(m.getnr(regname)).get_ycenter();
    }

  double mean1 = v1.mean(0);
  double mean2 = v2.mean(0);
  double * v1z = v1.getV();
  double * v2z = v2.getV();
  for(unsigned j=0;j<v1.rows();j++,v1z++,v2z++)
    {
    *v1z -= mean1;
    *v2z -= mean2;
    }

  create(v1,v2);

  all_precenv.erase(all_precenv.begin(),all_precenv.end());
  lambdavec.erase(lambdavec.begin(),lambdavec.end());

  lambdaold = -1;
  }


  // CONSTRUCTOR 7: geosplines varying coefficients

FULLCOND_pspline_surf_stepwise::FULLCOND_pspline_surf_stepwise(
                      MCMCoptions * o, DISTRIBUTION * dp, FULLCOND_const * fcc,
                      const datamatrix & intact,
                      const datamatrix & region, const MAP::map & mp, const ST::string & mn,
                      const ST::string & ti, const unsigned & nrk,
                      const unsigned & degr, const knotpos & kp, const double & l,
                      const int & gs, const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const bool & gauss, const bool & vccent, const unsigned & c)
  : FULLCOND_pspline_surf_gaussian(o,dp,fcc,intact,region,mp,mn,ti,nrk,degr,kp,l,gs,ft,fp,pres,true,c)
  {

  fctype = nonparametric;

  if(gauss==false)
    utype = iwls;

  datamatrix v1 = datamatrix(likep->get_nrobs(),1,0.0);
  datamatrix v2 = datamatrix(likep->get_nrobs(),1,0.0);

  ST::string regname;
  for(unsigned i=0;i<likep->get_nrobs();i++)
    {
    regname = ST::doubletostring(region(i,0));
    regionnames.push_back(regname);
    v1(i,0) = m.get_region(m.getnr(regname)).get_xcenter();
    v2(i,0) = m.get_region(m.getnr(regname)).get_ycenter();
    }

  double mean1 = v1.mean(0);
  double mean2 = v2.mean(0);
  double * v1z = v1.getV();
  double * v2z = v2.getV();
  for(unsigned j=0;j<v1.rows();j++,v1z++,v2z++)
    {
    *v1z -= mean1;
    *v2z -= mean2;
    }

  centertotal = true;
  //VCM_neu
  if(vccent == true)
    identifiable = false;

  data_forfixed = intact;

  unsigned nrobs = intact.rows();

  if (data_varcoeff_fix.rows() < nrobs)
    {
    data_varcoeff_fix = datamatrix(nrobs,2,1);
    effmodi = datamatrix(nrobs,1,0);
    for(unsigned i=0;i<nrobs;i++)
      {
      effmodi(i,0) = v1(i,0)*v2(i,0);
      data_varcoeff_fix(i,0) = intact(i,0);
      data_varcoeff_fix(i,1) = effmodi(i,0)*intact(i,0);
      }
    }

  centervalue = (v1.max(0) + v1.min(0)) * (v2.max(0) + v2.min(0));
  XVX = datamatrix(2,2,0);

  lambdaold = -1;
  }


FULLCOND_pspline_surf_stepwise::FULLCOND_pspline_surf_stepwise(const FULLCOND_pspline_surf_stepwise & fc)
  : FULLCOND_pspline_surf_gaussian(FULLCOND_pspline_surf_gaussian(fc))
  {
  mainpoi1 = fc.mainpoi1;
  mainpoi2 = fc.mainpoi2;
  maineffectsexisting = fc.maineffectsexisting;
  df_lambdaold = fc.df_lambdaold;
  lambdaold = fc.lambdaold;
  all_precenv = fc.all_precenv;
  lambdavec = fc.lambdavec;
  lambdaxold = fc.lambdaxold;
  lambdayold = fc.lambdayold;
  lambdax_prec = fc.lambdax_prec;
  lambday_prec = fc.lambday_prec;

  KHenv = fc.KHenv;
  Kxenv = fc.Kxenv;
  Kyenv = fc.Kyenv;
  KH = fc.KH;
  Kx = fc.Kx;
  Ky = fc.Ky;
  Kgrenz = fc.Kgrenz;
  Kalt = fc.Kalt;

  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  XVX = fc.XVX;
  centervalue = fc.centervalue;
  fc_df = fc.fc_df;
  betaold = fc.betaold;

  splineo1 = fc.splineo1;
  splineo2 = fc.splineo2;
  }


const FULLCOND_pspline_surf_stepwise & FULLCOND_pspline_surf_stepwise::operator=(
                                            const FULLCOND_pspline_surf_stepwise & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_pspline_surf_gaussian::operator=(FULLCOND_pspline_surf_gaussian(fc));

  mainpoi1 = fc.mainpoi1;
  mainpoi2 = fc.mainpoi2;
  maineffectsexisting = fc.maineffectsexisting;
  df_lambdaold = fc.df_lambdaold;
  lambdaold = fc.lambdaold;
  all_precenv = fc.all_precenv;
  lambdavec = fc.lambdavec;
  lambdaxold = fc.lambdaxold;
  lambdayold = fc.lambdayold;
  lambdax_prec = fc.lambdax_prec;
  lambday_prec = fc.lambday_prec;

  KHenv = fc.KHenv;
  Kxenv = fc.Kxenv;
  Kyenv = fc.Kyenv;
  KH = fc.KH;
  Kx = fc.Kx;
  Ky = fc.Ky;
  Kgrenz = fc.Kgrenz;
  Kalt = fc.Kalt;

  data_varcoeff_fix = fc.data_varcoeff_fix;
  effmodi = fc.effmodi;
  XVX = fc.XVX;
  centervalue = fc.centervalue;
  fc_df = fc.fc_df;
  betaold = fc.betaold;

  splineo1 = fc.splineo1;
  splineo2 = fc.splineo2;

  return *this;
  }


bool FULLCOND_pspline_surf_stepwise::posteriormode(void)
  {
  bool converged = false;
  bool converged1 = true;
  bool converged2 = true;

  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);
  fctotal.set_transform(transform);

  unsigned i;

  likep->substr_linearpred_m(spline,column,true);

  if(varcoeff && lambda == -2)      // Gerade beim VC
    {
    datamatrix betas = datamatrix(2,1,0);
    if(calculate_xwx_vc == true || (XVX(1,1)==0 && XVX(0,0)==0))
      {
      calculate_xwx_vc = false;
      likep->fisher(XVX,data_varcoeff_fix,column);
      XVX.assign((XVX.cinverse()));
      }
    likep->compute_weightiwls_workingresiduals(column);
    betas = XVX*data_varcoeff_fix.transposed()*likep->get_workingresiduals();
    spline.mult(data_varcoeff_fix,betas);
    likep->add_linearpred_m(spline,column);     // addiert durchschnittl. Fkt. zum Gesamtprädiktor

    if(center)
      {
      intercept = betas(0,0) + 0.25*betas(1,0)*centervalue;
      for(i=0;i<spline.rows();i++)
        spline(i,0) -= intercept*data_forfixed(i,0);
      betas(0,0) -= intercept;
      update_fix_effect();
      intercept = 0.0;
      }

    double * fchelpbetap = fchelp.getbetapointer();
    datamatrix help = datamatrix(beta.rows(),1,0);
    unsigned j = 0;
    if(gridsize<0)
      {
      vector<int>::iterator freqwork = freqoutput.begin();
      int * workindex = index.getV();
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
          {
          *fchelpbetap = betas(0,0) + effmodi(*workindex,0)*betas(1,0);
          while(j<beta.rows())
            {
            help(j,0) = *fchelpbetap;
            j += 1;
            }
          fchelpbetap++;
          }
        }
      }
    else
      {
      vector<double>::iterator effitx = effectvaluesx.begin();
      vector<double>::iterator effity = effectvaluesy.begin();
      for(i=0;i<unsigned(gridsize);i++,fchelpbetap++,effitx++,effity++)
        {
        *fchelpbetap = betas(0,0) + *effitx * *effity * betas(1,0);   // keine Fallunterscheidung VC / nicht VC nötig!
        while(j<beta.rows())
          {
          help(j,0) = *fchelpbetap;
          j += 1;
          }
        }
      }
    beta.assign(help);
    return fchelp.posteriormode();
    }

  else  // if(lambda!=-2)
    {
    if(type == mrfkr1 && rankK==(nrpar1dim-1)*(nrpar1dim-1))   // hier passen "type" und "rankK" zusammen, d.h. kein Spezialfall für Bootstrap
      {                                                        // außerdem nicht für "rw2" oder "rw1"
      mainpoi1->reset_effect(0);
      mainpoi2->reset_effect(0);

      double lambdax = mainpoi1->get_lambda() / nrpar1dim;
      double lambday = mainpoi2->get_lambda() / nrpar1dim;

      if (lambda_prec != lambda || lambdax_prec != lambdax || lambday_prec != lambday
            || calculate_xwx == true)
        {
        if (calculate_xwx == true)
          {
          calculate_xwx = false;
          compute_XWXenv(likep->get_weightiwls(),column);
          }
        if(Kyenv.getBandwidth() > Kxenv.getBandwidth())
          {
          Kenv.addto(Kyenv,KHenv,lambday,lambda);
          Kenv.addto(Kenv,Kxenv,1.0,lambdax);
          }
        else
          {
          Kenv.addto(Kxenv,KHenv,lambdax,lambda);
          Kenv.addto(Kenv,Kyenv,1.0,lambday);
          }
        //Kenv.addto(Kxenv,Kyenv,lambdax,lambday);
        //Kenv.addto(Kenv,KHenv,1.0,lambda);
        prec_env.addto(Kenv,XX_env,1.0,1.0);
        lambda_prec = lambda;
        lambdax_prec = lambdax;
        lambday_prec = lambday;
        }
      }
    else
      {
      if (lambda_prec != lambda || calculate_xwx == true)
        {
        if (calculate_xwx == true)
          {
          calculate_xwx = false;
          compute_XWXenv(likep->get_weightiwls(),column);
          }
        prec_env.addto(XX_env,Kenv,1.0,lambda);
        lambda_prec = lambda;
        }
      }

    likep->compute_workingresiduals(column);

    compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);

    prec_env.solve(muy,beta);

    if((lambda-1000000000)==0)   // MCMCbootstrap: entspricht fixem Effekt  (lambda -> unendlich und rw2)
      {
      compute_intercept();
      compute_beta();
      for(i=0;i<nrpar;i++)
        beta(i,0) += intercept;
      intercept = 0.0;
      }

    add_linearpred_multBS(beta);

    double intercept_save = 0.0;
    if(type == mrfkr1 && rankK==(nrpar1dim-1)*(nrpar1dim-1))
      {
      compute_intercept();
      if(!varcoeff)
        {
        compute_main();
        if(utype != gaussian)
          {
          beta_uncentered = beta;
          //betaold = beta;
          //for(i=0;i<nrpar;i++)
          //  betaold(i,0) -= intercept;
          }
        compute_beta();
        if(center)
          fcconst->posteriormode_intercept(intercept);
        mainpoi1->changeposterior3(beta1,he1,intercept);
        mainpoi2->changeposterior3(beta2,he2,intercept);
        if(!center)
          {
          double * work = spline.getV();
          for(i=0;i<spline.rows();i++,work++)
            *work += intercept;
          work = beta.getV();
          for(i=0;i<beta.rows();i++,work++)
            *work += intercept;
          }
        intercept_save = -intercept;
        }
      else
        {
        multBS_index(splinehelp,beta);
        compute_main_varcoeff();
        if(utype != gaussian)
          {
          beta_uncentered = beta;
          }
        compute_beta();
        if(center)
          update_fix_effect();
        mainpoi1->changeposterior_varcoeff(beta1,he1,intercept);
        mainpoi2->changeposterior_varcoeff(beta2,he2,intercept);
        if(!center)
          {
          double * work = beta.getV();
          for(i=0;i<beta.rows();i++,work++)
            *work += intercept;
          }
        }
      intercept = 0.0;
      }

    if(center)
      {
      if(centertotal || (lambda-1000000000)==0)  // Grenzfall bei Bootstrap oder "rw1" / "rw2"
        {
        compute_intercept();
        if(!varcoeff)
          {
          for(i=0;i<nrpar;i++)
            beta(i,0) -= intercept;
          for(i=0;i<likep->get_nrobs();i++)
            spline(i,0) -= intercept;
          }
        else
          {
          for(i=0;i<spline.rows();i++)
            spline(i,0) -= intercept*data_forfixed(i,0);
          }

        if(!varcoeff)
          fcconst->posteriormode_intercept(intercept);
        else
          update_fix_effect();

        intercept_save = intercept;
        intercept = 0.0;
        converged = FULLCOND_nonp_basis::posteriormode();
        }
      if(!centertotal)
        {
        // Gesamteffekt in fctotal schreiben
        converged = FULLCOND_nonp_basis::posteriormode();
        if(converged && converged1 && converged2)
          {
          double * fctotalbetap = fctotal.getbetapointer();

          if(gridsize < 0)
            {
            vector<int>::iterator freqwork = freq.begin();
            int * workindex = index.getV();
            for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
              {
              if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
                {
                if(!varcoeff)
                  *fctotalbetap = spline(*workindex,0)
                              + mainpoi1->get_spline()(*workindex,0)
                              + mainpoi2->get_spline()(*workindex,0);
                else
                  *fctotalbetap = splinehelp(*workindex,0)
                              + mainpoi1->get_splinehelp()(*workindex,0)
                              + mainpoi2->get_splinehelp()(*workindex,0);
                fctotalbetap++;
                }
              }
            }
          else
            {
            multDG(splinehelp,beta);
            unsigned k,l;
            for(k=0;k<unsigned(gridsizex);k++)
              for(l=0;l<unsigned(gridsizey);l++,fctotalbetap++)
                {
                *fctotalbetap = splinehelp(k*gridsizey + l,0)
                            + mainpoi1->get_splinehelp()(k,0)
                            + mainpoi2->get_splinehelp()(l,0);
                }
            }

          fctotal.posteriormode();
          }
        } // END: interaction
      } // END: if(center)
    else
      converged = FULLCOND_nonp_basis::posteriormode();

  // Interaktionseffekt in fchelp schreiben
      if(converged && converged1 && converged2)
        {
        double * fchelpbetap = fchelp.getbetapointer();
        if(gridsize < 0)
          {
          if(varcoeff && (type!=mrfkr1 || rankK!=(nrpar1dim-1)*(nrpar1dim-1)))
            {
            multBout(splinehelp,beta);
              if(center)
                {
                int * workindex = index.getV();
                for(i=0;i<splinehelp.rows();i++,workindex++)
                  splinehelp(i,0) -= intercept_save;
                }
             }

          vector<int>::iterator freqwork = freq.begin();
          int * workindex = index.getV();
          for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
            {
            if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
              {
              if(varcoeff)
                *fchelpbetap = splinehelp(i,0);
              else
                *fchelpbetap = spline(*workindex,0);
              fchelpbetap++;
              }
            }
          }
        else
          {
          multDG(splinehelp,beta);  // NICHT zeilen- und spaltenweise zentriert!!!
                                // dazu: if(center && !centertotal){multDG(splinehelp,beta);}
                                // und: 'splinehelp' ändern 'in compute_maineffects'
          for(i=0;i<unsigned(gridsize);i++,fchelpbetap++)
            *fchelpbetap = splinehelp(i,0) - intercept_save;
          }
        fchelp.posteriormode();
        }

    if(converged && converged1 && converged2)
      return true;
    else
      return false;
    }
  }


void FULLCOND_pspline_surf_stepwise::reset_effect(const unsigned & pos)
  {
  likep->substr_linearpred_m(spline,column,true);

  unsigned i;
  double * work;
  work = spline.getV();
  for(i=0;i<spline.rows();i++,work++)
    *work = 0.0;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;
  intercept = 0.0;
  }


void FULLCOND_pspline_surf_stepwise::reset(void)
  {
  spline = datamatrix(spline.rows(),spline.cols(),0);
  FULLCOND::reset();
  }

void FULLCOND_pspline_surf_stepwise::remove_centering_fix(void)
  {
  }


void FULLCOND_pspline_surf_stepwise::hierarchical(ST::string & possible)
  {
  bool raus = false;
  bool fix = false;
  bool spline1, fix1, fix2;

  if(!varcoeff)
    {
    if(maineffectsexisting == 11)
      {
      mainpoi1->get_inthemodel(spline1,fix1);
      if(spline1 == false && fix1 == false)
        raus = true;
      mainpoi2->get_inthemodel(spline1,fix2);
      if(spline1 == false && fix2 == false)
        raus = true;
      if(raus == false && (fix1 == true || fix2 == true))
        fix = true;
      }

    if(raus == true)
      possible = "raus";
    else if(fix == false && raus == false)
      possible = "alles";
    else
      possible = "rfix";
    }
  else if(varcoeff)    // nur für Unterscheidung VC / !VC
    {
    if(maineffectsexisting == 11)
      {
      mainpoi1->get_inthemodel(spline1,fix1);
      if(spline1 == false && fix1 == false)
        raus = true;
      mainpoi2->get_inthemodel(spline1,fix2);
      if(spline1 == false && fix2 == false)
        raus = true;
      if(raus == false && (mainpoi1->get_lambda()==-2 || mainpoi2->get_lambda()==-2))
        fix = true;

      if(raus == true)
        possible = "vraus";
      else if(fix == false && raus == false)
        possible = "valles";
      else
        possible = "vrfix";
      }
    else
      possible = "valles";
    }
  }


void FULLCOND_pspline_surf_stepwise::set_pointer_to_interaction(FULLCOND * inter)
  {
  interactions_pointer.push_back(inter);
  }


void FULLCOND_pspline_surf_stepwise::get_interactionspointer(vector<FULLCOND*> & inter)
  {
  inter = interactions_pointer;
  }


void FULLCOND_pspline_surf_stepwise::hierarchie_rw1(vector<double> & untervector, int dfo)
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
          hilf.push_back(-1);
          for(i=stelle_oben;i<untervector.size();i++)
            hilf.push_back(untervector[i]);
          untervector = hilf;
          }
        }
     }
  else if(df_min >= dfo)
     {
     if(varcoeff)
       untervector.push_back(-2);
     else
       untervector.push_back(-1);
     }
  else
     {
     vector<double> hilf;
     if(varcoeff)
       hilf.push_back(-2);
     else
       hilf.push_back(-1);
     unsigned i;
     for(i=0;i<untervector.size();i++)
        hilf.push_back(untervector[i]);
     untervector = hilf;
     }
  }


void FULLCOND_pspline_surf_stepwise::compute_lambdavec(
vector<double> & lvec, int & number)
  {
  if(spfromdf=="automatic")
    {
    df_equidist = true;
    double maxi = floor(nrpar/4.0*3.0);

    if(maxi<=60)
      {
      number = floor(maxi/2);
      df_for_lambdamax = 4;
      df_for_lambdamin = number*2+2;
      }
    else if(maxi > 60 && maxi<=100)
      {
      df_for_lambdamax = 4;
      number = floor(maxi/3);
      df_for_lambdamin = number*3+1;
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

  if(number>0)
    {
    if (df_equidist==true && spfromdf!="direct" && number>1)
       FULLCOND::compute_lambdavec_equi(lvec,number);
    else
       FULLCOND::compute_lambdavec(lvec,number);
    }

  if(centertotal && likep->iwlsweights_constant() == true)
    {
    lambdavec = lvec;
    compute_XWXenv(likep->get_weightiwls(),column);
    for(unsigned i=0;i<lambdavec.size();i++)
      {
      prec_env.addto(XX_env,Kenv,1.0,lambdavec[i]);
      prec_env.decomp();
      all_precenv.push_back(prec_env);
      }
    }

  if(!nofixed && !varcoeff)
    {
    if((type==RW1 || type==mrflinear) && number>0)
      hierarchie_rw1(lvec,1);
    else  // if(varcoeff || type==RW2 || (type==RW1 && number==-1))
      lvec.push_back(-1);
    }
  else
    {
    if(!nofixed && (type==RW1 || type==mrflinear) && number>0)
      {
      if(identifiable)
        {
        hierarchie_rw1(lvec,2);
        lvec.push_back(-1);
        }
      else
        hierarchie_rw1(lvec,1);
      }
    else if(!nofixed)
      {
      lvec.push_back(-2);
      if(identifiable)
        lvec.push_back(-1);
      }
    }

  if(forced_into==false)
     lvec.push_back(0);

  // Startwert für lambda aus df:
  if(spfromdf!="direct")
    {
    double lambdavorg = 1000;
    if(!nofixed && !varcoeff)
      {
      if(dfstart==1)
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
      else if(!nofixed && (dfstart==2 && identifiable) || (dfstart==1 && !identifiable))
        lambdastart = -2;
      else if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    if(lambdastart==-9 || lambdastart==1000000000)    // falls dfstart nicht erreicht werden kann
      lambdastart = 0;
    }

  }


void FULLCOND_pspline_surf_stepwise::update_stepwise(double la)
  {
  lambda=la;

  if(centertotal)
    {
    if(likep->iwlsweights_constant() == true)
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
        prec_env = all_precenv[i-1];
        lambda_prec = lambda;
        }
      }
    }
  }


double FULLCOND_pspline_surf_stepwise::compute_df(void)
  {
  double df = 0;
  //if(inthemodel == false && fixornot == true)
  //  df = 1;
  if(inthemodel == true)
    {
    if(varcoeff && lambda == -2)
      {
      if(identifiable)
        df = 2;
      else
        df = df + 1;
      }
    else if(type == mrfkr1 && rankK==(nrpar1dim-1)*(nrpar1dim-1))   // d.h. Typ und Rang passen zusammen --> kein Spezialfall für Bootstrap
      {
      double lambdax = 0;
      double lambday = 0;
      bool fix,drin;    // für das Initialisieren der Startwerte!
      mainpoi1->get_inthemodel(drin,fix);
      if(drin==true)
        {
        lambdax = mainpoi1->get_lambda() / nrpar1dim;
        lambday = mainpoi2->get_lambda() / nrpar1dim;
        }
      if(lambda != lambdaold || lambdax != lambdaxold || lambday != lambdayold
           || likep->get_iwlsweights_notchanged() == false)
        {
        if (calculate_xwx == true)
          compute_XWXenv(likep->get_weightiwls(),column);
        if(lambda != lambda_prec || lambdax != lambdax_prec || lambday != lambday_prec
             || calculate_xwx == true)
          {
          calculate_xwx = false;
          if(Kyenv.getBandwidth() > Kxenv.getBandwidth())
            {
            Kenv.addto(Kyenv,KHenv,lambday,lambda);
            Kenv.addto(Kenv,Kxenv,1.0,lambdax);
            }
          else
            {
            Kenv.addto(Kxenv,KHenv,lambdax,lambda);
            Kenv.addto(Kenv,Kyenv,1.0,lambday);
            }
          //Kenv.addto(Kxenv,Kyenv,lambdax,lambday);
          //Kenv.addto(Kenv,KHenv,1.0,lambda);
          prec_env.addto(Kenv,XX_env,1.0,1.0);
          lambda_prec = lambda;
          lambdax_prec = lambdax;
          lambday_prec = lambday;
          }
        invprec = envmatdouble(0,nrpar,prec_env.getBandwidth());
        prec_env.inverse_envelope(invprec);
        df = df + invprec.traceOfProduct(XX_env);
        if(!identifiable)
          df -= 1;
        //df = FULLCOND_pspline_surf_gaussian::compute_df();
        df_lambdaold = df;
        lambdaold = lambda;
        lambdaxold = lambdax;
        lambdayold = lambday;
        }
      else
        df = df_lambdaold;

      if(maineffectsexisting != 0)
        {
        double gesamt = df;
        double df_h1 = 0;
        double df_h2 = 0;
        if(drin==true)
          {
          df_h1 += mainpoi1->compute_df();
          df_h2 += mainpoi2->compute_df();
          }
        df = gesamt - df_h1 - df_h2;
        }
      }
    else
      {
      if(lambda != lambdaold || likep->get_iwlsweights_notchanged() == false)
        {
        if (calculate_xwx == true)
          compute_XWXenv(likep->get_weightiwls(),column);
        if(lambda != lambda_prec || calculate_xwx == true)
          {
          calculate_xwx = false;
          prec_env.addto(XX_env,Kenv,1.0,lambda);
          lambda_prec = lambda;
          }
        invprec = envmatdouble(0,nrpar,prec_env.getBandwidth());
        prec_env.inverse_envelope(invprec);
        df = df + invprec.traceOfProduct(XX_env);
        if(!identifiable)
          df -= 1;
        df_lambdaold = df;
        lambdaold = lambda;
        }
      else
        df = df_lambdaold;
      }
    }
  return df;
  }


const datamatrix & FULLCOND_pspline_surf_stepwise::get_data_forfixedeffects(void)
  {
  return data_forfixed;
  }


ST::string FULLCOND_pspline_surf_stepwise::get_effect(void)
  {
  ST::string h;

  if(varcoeff)
    {
    if(datanames.size()==5)
      {
      if(type==mrflinear)
        h = datanames[1] + "*" + datanames[3] + "(pspline2dimrw1";
      else if(type==mrfquadratic8)
        h = datanames[1] + "*" + datanames[3] + "(pspline2dimrw2";
      else if(type==mrfkr1)
        h = datanames[1] + "*" + datanames[3] + "(psplineinteract";
      }
    else
      {
      if(type==mrflinear)
        h = datanames[1] + "*" + datanames[0] + "(geosplinerw1";
      else if(type==mrfquadratic8)
        h = datanames[1] + "*" + datanames[0] + "(geosplinerw2";
      }
    }
  else
    {
    if(type==mrflinear)
      {
      if(datanames.size()>1)
        h = datanames[2] + "(pspline2dimrw1";
      else
        h = datanames[0] + "(geosplinerw1";
      }
    else if(type==mrfquadratic8)
      {
      if(datanames.size()>1)
        h = datanames[2] + "(pspline2dimrw2";
      else
        h = datanames[0] + "(geosplinerw2";
      }
    else
      h = datanames[2] + "(psplineinteract";
    }

  h = h + ",df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;
  }


// BEGIN: For Varying Coefficients ---------------------------------------------

void FULLCOND_pspline_surf_stepwise::update_fix_effect(void)
  {
  bool raus = false;
  unsigned j = 1;
  ST::string name_richtig = datanames[1];
  while(j<fcconst->get_datanames().size() && raus==false)
     {
     if(fcconst->get_datanames()[j] == datanames[1])
        {
        raus = true;
        }
     if(fcconst->get_datanames()[j] == (datanames[1]+"_1"))
        {
        raus = true;
        name_richtig = datanames[1] + "_1";
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

void FULLCOND_pspline_surf_stepwise::const_varcoeff(void)
  {
  if(varcoeff)
    fcconst->posteriormode_const_varcoeff(data_forfixed);
  }

// END: For Varying Coefficients -----------------------------------------------

// BEGIN: For Bootstrap --------------------------------------------------------

void FULLCOND_pspline_surf_stepwise::update_bootstrap(const bool & uncond)
  {
  update_bootstrap_df();

  if(fixornot==true)
    {
    bool raus = false;
    unsigned j = 1;
    ST::string name_richtig;
    if(!varcoeff)
      name_richtig = datanames[datanames.size()-1];
    else
      name_richtig = datanames[1];
    while(j<fcconst->get_datanames().size() && raus==false)
      {
      if(fcconst->get_datanames()[j] == name_richtig)
        raus = true;
      j = j + 1;
      }
    unsigned index_fix = j-1;
    double fix = fcconst->getbeta(index_fix,0);
    unsigned i;
    double korrektur = 0;
    if(center)
      korrektur = -0.25*fix*centervalue;
    double * fchelpbetap = fchelp.getbetapointer();
    if(gridsize < 0)                              // alle verschiedene Beobachtungen
      {
      vector<int>::iterator freqwork = freq.begin();
      int * workindex = index.getV();
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
          {
          if(!varcoeff)
            *fchelpbetap = fix * data_forfixed(*workindex,0) + korrektur;
          else
            *fchelpbetap = fix;
          fchelpbetap++;
          }
        }
      }
    else //if(gridsize>0) // Gitterpunkte
      {
      vector<double>::iterator effitx = effectvaluesx.begin();
      vector<double>::iterator effity = effectvaluesy.begin();

      for(i=0;i<unsigned(gridsize);i++,fchelpbetap++,effitx++,effity++)
        {
        if(!varcoeff)
          *fchelpbetap = fix * *effitx * *effity + korrektur;
        else
          *fchelpbetap = fix;
        }
      }
    double help = -korrektur;
    fcconst->update_intercept(help);
    fchelp.update_bootstrap();
    }
  else if(inthemodel==false && fixornot==false)
    {
    double * fchelpbetap = fchelp.getbetapointer();
    for(unsigned i=0;i<fchelp.getbeta().rows();i++,fchelpbetap++)
        *fchelpbetap = 0;
    fchelp.update_bootstrap();
    }
  else
    {
    fchelp.update_bootstrap();
    }
  }


void FULLCOND_pspline_surf_stepwise::update_bootstrap_df(void)
  {
  if(optionsp->get_nriter()<=1)
    {
    ST::string path = samplepath.substr(0,samplepath.length()-4)+"_df.raw";
if(centertotal == false)
  fc_df = FULLCOND(optionsp,datamatrix(1,1),"title?",2,1,path);
else
    fc_df = FULLCOND(optionsp,datamatrix(1,1),"title?",1,1,path);
    fc_df.setflags(MCMC::norelchange | MCMC::nooutput);

    //if(centertotal == false)
    //  {
    //  path = samplepath.substr(0,samplepath.length()-4)+"_total_df.raw";
    //  fctotal_df = FULLCOND(optionsp,datamatrix(1,1),"title?",3,1,path);
    //  fctotal_df.setflags(MCMC::norelchange | MCMC::nooutput);
    //  }
    }

if(centertotal == false)
  {
  double df = compute_df() + mainpoi1->compute_df() + mainpoi2->compute_df();
  fc_df.setbetavalue(1,0,df);
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
  else
    {
    fc_df.setbetavalue(0,0,lambda);
    fc_df.update_bootstrap_df();
    }

  //if(centertotal == false)
  //  {
  //  double df = compute_df() + mainpoi1->compute_df() + mainpoi2->compute_df();
  //  fctotal_df.setbetavalue(0,0,df);
  //  fctotal_df.update_bootstrap_df();
  //  }
  }


void FULLCOND_pspline_surf_stepwise::save_betamean(void)
  {
  if(fixornot==true)
    {
    bool raus = false;
    unsigned j = 1;
    ST::string name_richtig;
    if(!varcoeff)
      name_richtig = datanames[datanames.size()-1];
    else
      name_richtig = datanames[1];
    while(j<fcconst->get_datanames().size() && raus==false)
      {
      if(fcconst->get_datanames()[j] == name_richtig)
        raus = true;
      j = j + 1;
      }
    unsigned index_fix = j-1;
    double fix = fcconst->getbeta(index_fix,0);
    unsigned i;
    double korrektur = 0;
    if(center)
      korrektur = -0.25*fix*centervalue;
    double * fchelpbetap = fchelp.getbetapointer();
    if(gridsize < 0)                              // alle verschiedene Beobachtungen
      {
      vector<int>::iterator freqwork = freq.begin();
      int * workindex = index.getV();
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
          {
          if(!varcoeff)
            *fchelpbetap = fix * data_forfixed(*workindex,0) + korrektur;
          else
            *fchelpbetap = fix;
          fchelpbetap++;
          }
        }
      }
    else //if(gridsize>0) // Gitterpunkte
      {
      vector<double>::iterator effitx = effectvaluesx.begin();
      vector<double>::iterator effity = effectvaluesy.begin();

      for(i=0;i<unsigned(gridsize);i++,fchelpbetap++,effitx++,effity++)
        {
        if(!varcoeff)
          *fchelpbetap = fix * *effitx * *effity + korrektur;
        else
          *fchelpbetap = fix;
        }
      }
    double help = -korrektur;
    fcconst->update_intercept(help);
    fchelp.save_betamean();
    }
  else if(inthemodel==false && fixornot==false)
    {
    double * fchelpbetap = fchelp.getbetapointer();
    for(unsigned i=0;i<fchelp.getbeta().rows();i++,fchelpbetap++)
        *fchelpbetap = 0;
    fchelp.save_betamean();
    }
  else
    {
    fchelp.save_betamean();
    }
  }


void FULLCOND_pspline_surf_stepwise::update_bootstrap_betamean(void)
  {
  fchelp.update_bootstrap_betamean();
  FULLCOND::setflags(MCMC::norelchange);
  }

void FULLCOND_pspline_surf_stepwise::outresults_df(unsigned & size)
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

//pathdf = pathcurrent.substr(0,pathcurrent.length()-4)+"_df_sample.raw";
//ofstream out(pathdf.strtochar());
//sample.prettyPrint(out);

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

  for(i=0;i<number.size();i++)
    {
    double help = sample.get(index(cumnumber[i]-1,0),0);
    double dfs = -1*help;
    if(help>0)
      {
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

if(centertotal == false)
  {
  fc_df.readsample_df(sample,1);
  pathdf = pathcurrent.substr(0,pathcurrent.length()-4)+"_total_df.raw";
  ofstream out(pathdf.strtochar());
  sample.prettyPrint(out);
  }

  }


void FULLCOND_pspline_surf_stepwise::update(void)
  {
  if(type==mrfkr1 && optionsp->get_nriter()==1)
    {
    mainp1 = mainpoi1;
    mainp2 = mainpoi2;
    }
  if(utype != gaussian && optionsp->get_nriter()==1)
    {
    updateW = 1;
    proposal = datamatrix(nrpar,1,0);
    }

  if(lambda==0)     // Nullfunktion rausschreiben
    {
    beta = datamatrix(beta.rows(),beta.cols(),0);
    FULLCOND::update();
    double * fchelpbetap = fchelp.getbetapointer();
    for(unsigned i=0;i<fchelp.getbeta().rows();i++,fchelpbetap++)
        *fchelpbetap = 0;
    fchelp.update();

    if(center && !centertotal)     // Gesamteffekt in fctotal schreiben
      {
      if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
          ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
        {
        unsigned i;
        double * fctotalbetap = fctotal.getbetapointer();
        if(gridsize < 0)
          {
          vector<int>::iterator freqwork = freqoutput.begin();
          int * workindex = index.getV();
          for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
            {
            if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
              {
              if(!varcoeff)
                *fctotalbetap = mainpoi1->get_spline()(*workindex,0)
                              + mainpoi2->get_spline()(*workindex,0);
              else
                *fctotalbetap = mainpoi1->get_splinehelp()(*workindex,0)
                              + mainpoi2->get_splinehelp()(*workindex,0);
              fctotalbetap++;
              }
            }
          }
        else
          {
          multDG(splinehelp,beta);
          unsigned k,l;
          for(k=0;k<unsigned(gridsizex);k++)
            for(l=0;l<unsigned(gridsizey);l++,fctotalbetap++)
              *fctotalbetap = mainp1->get_splinehelp()(k,0)
                            + mainp2->get_splinehelp()(l,0);
          }
        }   // ENDE: Gesamteffekt in fctotal schreiben
      fctotal.update();
      }
    }
  else if((lambda-1000000000)==0)     // lineare Interaktion wird erreicht durch Grenzfall (lambda = 1000000000 und "rw2")
    {
    update_linear_function();
    }
  else
    {
    if(type==mrfkr1 && rankK==(nrpar1dim-1)*(nrpar1dim-1))
      {
      if(utype == gaussian)
        {
        // Bandmatrizen nötig für Update-Funktion!
        // "K" ist gesamte Matrix, deshalb KH = rw1#rw1
        if(optionsp->get_nriter() == 1 || KH.bandsize()!=(nrpar1dim+1))
          {
          KH = bandmatdouble(nrpar,K.bandsize(),0);
          KH = K;
          K = bandmatdouble(nrpar,Kenv.getBandwidth(),0);
          }
        double lambdax = mainpoi1->get_lambda() / nrpar1dim;
        double lambday = mainpoi2->get_lambda() / nrpar1dim;
        if(Ky.bandsize() > Kx.bandsize())
          {
          if(Ky.bandsize() < KH.bandsize())
            K.addto2(Ky,KH,lambday,lambda);
          else
            K.addto2(KH,Ky,lambda,lambday);
          K.addto2(Kx,K,lambdax,1.0);
          }
        else
          {
          if(Kx.bandsize() < KH.bandsize())
            K.addto2(Kx,KH,lambdax,lambda);
          else
            K.addto2(KH,Kx,lambda,lambdax);
          K.addto2(Ky,K,lambday,1.0);
          }
        }
      else
        {
        beta.assign(beta_uncentered);
        compute_intercept();
        for(unsigned i=0;i<nrpar;i++)
          beta(i,0) -= intercept;
        intercept = 0.0;
        }

      set_lambdaconst(1);
      mainpoi1->reset_effect(0);  // Haupteffekte werden hier mitgeschätzt
      mainpoi2->reset_effect(0);
      }
    /*else
      {
      if(type==mrfkr1 && rankK!=(nrpar1dim-1)*(nrpar1dim-1))
        centertotal = true;
      }   */

    if(!varcoeff || centertotal)
      {
      if(utype == gaussian)
        FULLCOND_pspline_surf_gaussian::update();
      else
        update_IWLS();
      }
    else
      update_vc_anova();

    /*if(type==mrfkr1 && rankK!=(nrpar1dim-1)*(nrpar1dim-1))
      centertotal = false; */
    }
  }


void FULLCOND_pspline_surf_stepwise::update_IWLS(void)
  {
  double invscale = 1.0/likep->get_scale(column);
  sigma2 = likep->get_scale(column)/lambda;

  if(betahelp.rows() != beta.rows() || W.rows() != likep->get_nrobs())
    {
    betahelp = datamatrix(nrpar,1,0);
    W = datamatrix(likep->get_nrobs(),1,0);
    }

  double logold = - 0.5*Kenv.compute_quadform(beta,0)/sigma2;

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    logold += likep->compute_IWLS(W,mu,true,column,true);
    mu.plus(spline,mu);
    compute_XWXenv(W,0);  // W hat nur eine Spalte!!!
    }
  else
    {
    logold += likep->compute_IWLS(W,mu,false,column,true);
    mu.plus(mu,spline);
    }

  compute_XWtildey(W,invscale);

  prec_env.addto(XX_env,Kenv,invscale,1.0/sigma2);

  double * work = proposal.getV();
  for(unsigned i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.solve(muy,betahelp);
  prec_env.solveU(proposal,betahelp);

  add_linearpred_multBS2(proposal);

  betahelp.minus(proposal,betahelp);
  double qold = - 0.5*prec_env.compute_quadform(betahelp,0);

  double lognew = - 0.5*Kenv.compute_quadform(proposal,0)/sigma2;

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    qold += 0.5*prec_env.getLogDet();

    lognew += likep->compute_IWLS(W,mu,true,column,true);
    mu.plus(spline,mu);
    compute_XWXenv(W,0);  // W hat nur eine Spalte!!!
    prec_env.addto(XX_env,Kenv,invscale,1.0/sigma2);
    }
  else
    {
    lognew += likep->compute_IWLS(W,mu,false,column,true);
    mu.plus(mu,spline);
    }

  compute_XWtildey(W,invscale);

  prec_env.solve(muy,betahelp);

  betahelp.minus(beta,betahelp);
  double qnew = - 0.5*prec_env.compute_quadform(betahelp,0);

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
    beta.assign(proposal);
    }
  else
    {
    add_linearpred_multBS2(beta);
    }

  if(center)
    {
    if(!centertotal)
      {
// Gesamteffekt in fctotal schreiben
      if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
          ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
        {
        unsigned i;
        double * fctotalbetap = fctotal.getbetapointer();

        if(gridsize < 0)
          {
          vector<int>::iterator freqwork = freqoutput.begin();
          int * workindex = index.getV();
          for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
            {
            if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
              {
              *fctotalbetap = spline(*workindex,0)
                            + mainp1->get_spline()(*workindex,0)
                            + mainp2->get_spline()(*workindex,0);
              fctotalbetap++;
              }
            }
          }
        else
          {
          multDG(splinehelp,beta);
          unsigned k,l;
          for(k=0;k<gridsizex;k++)
            for(l=0;l<gridsizey;l++,fctotalbetap++)
              *fctotalbetap = splinehelp(k*gridsizey + l,0)
                            + mainp1->get_splinehelp()(k,0)
                            + mainp2->get_splinehelp()(l,0);
          }
        }   // ENDE: Gesamteffekt in fctotal schreiben
      fctotal.update();
      } // END: interaction
    } // END: if(center)

// Interaktionseffekt in fchelp schreiben
  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
    {
    double * fchelpbetap = fchelp.getbetapointer();
    unsigned i;
    if(gridsize < 0)
      {
      if(varcoeff)
        multBout(splinehelp,beta);

      vector<int>::iterator freqwork = freqoutput.begin();
      int * workindex = index.getV();
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          if(varcoeff)
            *fchelpbetap = splinehelp(i,0);
          else
            *fchelpbetap = spline(*workindex,0);
          fchelpbetap++;
          }
        }
      }
    else
      {
      multDG(splinehelp,beta);
      for(i=0;i<gridsize;i++,fchelpbetap++)
        *fchelpbetap = splinehelp(i,0);
      }
    }

  fchelp.update();
  FULLCOND::update();
  }


void FULLCOND_pspline_surf_stepwise::update_vc_anova(void)
  {
  if(lambdaconst == true)
    sigma2 = likep->get_scale(column)/lambda;

  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);
  fctotal.set_transform(transform);

  unsigned i;

  if(utype == gaussian)
    {
    double scaleinv = 1.0/likep->get_scale(column);

    if(changingweight || optionsp->get_nriter()==1)
      compute_XWXenv(likep->get_weight());

    likep->substr_linearpred_m(spline,column,true);

    if(XX.bandsize()<=K.bandsize())
      prec.addto2(XX,K,scaleinv,1.0/sigma2);
    else
      prec.addto2(K,XX,1.0/sigma2,scaleinv);

    double * work = standnormal.getV();
    for(i=0;i<nrpar;i++,work++)
      *work = rand_normal();

    prec.solveL(standnormal,beta);
    likep->compute_respminuslinpred(mu,column);   // nicht ändern wegen multgaussian
    compute_XWtildey(likep->get_weight(),scaleinv);
    prec.solve(muy,betahelp,0,0);
    beta.plus(beta,betahelp);

    if(samplecentered)
      if(center)
        sample_centered(beta);

    add_linearpred_multBS(beta);

    acceptance++;
    }  // end: gauss
  else  // IWLS
    {
    update_vc_anova_nongauss();
    }

  if(center)
    {
    if(!centertotal)
      {
      beta_uncentered.assign(beta);

      compute_intercept();
      multBS_index(splinehelp,beta);
      compute_main_varcoeff();
      compute_beta();
      fcconst->update_fix_varcoeff(intercept,datanames[1]);
      mainpoi1->change_varcoeff(beta1,he1,intercept);
      mainpoi2->change_varcoeff(beta2,he2,intercept);
      intercept = 0.0;

// Gesamteffekt in fctotal schreiben
      if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
          ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
        {
        double * fctotalbetap = fctotal.getbetapointer();
        if(gridsize < 0)
          {
          vector<int>::iterator freqwork = freqoutput.begin();
          int * workindex = index.getV();
          for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
            {
            if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
              {
              *fctotalbetap = splinehelp(*workindex,0)
                            + mainp1->get_splinehelp()(*workindex,0)
                            + mainp2->get_splinehelp()(*workindex,0);
              fctotalbetap++;
              }
            }
          }
        else
          {
          multDG(splinehelp,beta);
          unsigned k,l;
          for(k=0;k<unsigned(gridsizex);k++)
            for(l=0;l<unsigned(gridsizey);l++,fctotalbetap++)
              *fctotalbetap = splinehelp(k*gridsizey + l,0)
                            + mainp1->get_splinehelp()(k,0)
                            + mainp2->get_splinehelp()(l,0);
          }
        }   // ENDE: Gesamteffekt in fctotal schreiben

      fctotal.update();
      } // END: interaction
    } // END: if(center)

// Interaktionseffekt in fchelp schreiben
  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
    {
    double * fchelpbetap = fchelp.getbetapointer();

    if(gridsize < 0)
      {
      if(varcoeff)
        multBout(splinehelp,beta);

      vector<int>::iterator freqwork = freqoutput.begin();
      int * workindex = index.getV();
      for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          *fchelpbetap = splinehelp(i,0);
          fchelpbetap++;
          }
        }
      }
    else
      {
      multDG(splinehelp,beta);
      for(i=0;i<unsigned(gridsize);i++,fchelpbetap++)
        *fchelpbetap = splinehelp(i,0);
      }
    }

  fchelp.update();
  FULLCOND::update();
  }


void FULLCOND_pspline_surf_stepwise::update_vc_anova_nongauss(void)
  {
  sigma2 = likep->get_scale(column)/lambda;
  double invscale = 1.0/likep->get_scale(column);

  double logold = - 0.5*Kenv.compute_quadform(beta,0)/sigma2;

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
        ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    logold += likep->compute_IWLS(W,mu,true,column,true);
    mu.plus(spline,mu);
    compute_XWXenv(W,0);    // W hat nur eine Spalte!!!
    }
  else
    {
    logold += likep->compute_IWLS(W,mu,false,column,true);
    mu.plus(mu,spline);
    }

  compute_XWtildey(W,invscale);
  //prec_env.addto(XX_env,Kenv,iwlsscale,iwlsscale/sigma2);
  prec_env.addto(XX_env,Kenv,invscale,1.0/sigma2);

  double * work = proposal.getV();
  for(unsigned i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.solve(muy,betahelp);
  prec_env.solveU(proposal,betahelp);

  add_linearpred_multBS2(proposal);

  betahelp.minus(proposal,betahelp);
  double qold = - 0.5*prec_env.compute_quadform(betahelp,0);

  double lognew = - 0.5*Kenv.compute_quadform(proposal,0)/sigma2;

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    qold += 0.5*prec_env.getLogDet();

    lognew += likep->compute_IWLS(W,mu,true,column,true);
    mu.plus(spline,mu);
    compute_XWXenv(W,0);     // W hat nur eine Spalte!!!
    prec_env.addto(XX_env,Kenv,invscale,1.0/sigma2);
    }
  else
    {
    lognew += likep->compute_IWLS(W,mu,false,column,true);
    mu.plus(mu,spline);
    }

  compute_XWtildey(W,invscale);
  prec_env.solve(muy,betahelp);

  betahelp.minus(beta,betahelp);
  double qnew = - 0.5*prec_env.compute_quadform(betahelp,0);

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
    beta.assign(proposal);
    }
  else
    {
    add_linearpred_multBS2(beta);
    }
  }


void FULLCOND_pspline_surf_stepwise::update_linear_function(void)
  {
  if(lambdaconst == true)
    sigma2 = likep->get_scale(column)/lambda;

  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);
  fctotal.set_transform(transform);

//  unsigned i;

  if(utype == iwls)
    {
    double invscale = 1.0/likep->get_scale(column);

    double logold = - 0.5*Kenv.compute_quadform(beta,0)/sigma2;

    if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
        ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
      {
      logold += likep->compute_IWLS(W,mu,true,column,true);
      mu.plus(spline,mu);
      compute_XWXenv(W,0);    // W hat nur eine Spalte!!!
      }
    else
      {
      logold += likep->compute_IWLS(W,mu,false,column,true);
      mu.plus(mu,spline);
      }

    compute_XWtildey(W,invscale);

    prec_env.addto(XX_env,Kenv,invscale,1.0/sigma2);

    double * work = proposal.getV();
    for(unsigned int i=0;i<nrpar;i++,work++){
      *work = rand_normal();
      }

    prec_env.solve(muy,betahelp);
    prec_env.solveU(proposal,betahelp);

    // hier bezüglich Haupteffekten zentrieren und diese "wegschmeißen"! (damit nur Anteil "x1*x2" bleibt)
    datamatrix betasave = beta;
    beta.assign(proposal);
    compute_intercept();
    compute_beta();
    for(unsigned int i=0;i<nrpar;i++){
      beta(i,0) += intercept;
      }
    intercept = 0.0;
    proposal.assign(beta);
    beta.assign(betahelp);
    compute_intercept();
    compute_beta();
    for(unsigned int i=0;i<nrpar;i++){
      beta(i,0) += intercept;
      }
    intercept = 0.0;
    betahelp.assign(beta);
    beta.assign(betasave);

    add_linearpred_multBS2(proposal);

    betahelp.minus(proposal,betahelp);
    double qold = - 0.5*prec_env.compute_quadform(betahelp,0);

    double lognew = - 0.5*Kenv.compute_quadform(proposal,0)/sigma2;

    if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
        ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
      {
      qold += 0.5*prec_env.getLogDet();

      lognew += likep->compute_IWLS(W,mu,true,column,true);
      mu.plus(spline,mu);
      compute_XWXenv(W,0);    // W hat nur eine Spalte!!!
      prec_env.addto(XX_env,Kenv,invscale,1.0/sigma2);
      }
    else
      {
      lognew += likep->compute_IWLS(W,mu,false,column,true);
      mu.plus(mu,spline);
      }

    compute_XWtildey(W,invscale);

    prec_env.solve(muy,betahelp);

    // hier bezüglich Haupteffekten zentrieren und diese "wegschmeißen"!
    betasave.assign(beta);
    beta.assign(betahelp);
    compute_intercept();
    compute_beta();
    for(unsigned int i=0;i<nrpar;i++){
      beta(i,0) += intercept;
      }
    intercept = 0.0;
    betahelp.assign(beta);
    beta.assign(betasave);

    betahelp.minus(beta,betahelp);
    double qnew = - 0.5*prec_env.compute_quadform(betahelp,0);

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
      beta.assign(proposal);
      compute_intercept();
      fcconst->update_intercept(intercept);
      for(unsigned int i=0;i<nrpar;i++){
        beta(i,0) -= intercept;
        }
      for(unsigned int i=0;i<likep->get_nrobs();i++){
        spline(i,0) -= intercept;
        }
      intercept = 0.0;
      }
    else
      {
      add_linearpred_multBS2(beta);
      }
    }

  else if(utype == gaussian)
    {
    double scaleinv = 1.0/likep->get_scale(column);

    if(changingweight || optionsp->get_nriter()==1)
      compute_XWXenv(likep->get_weight());

    likep->substr_linearpred_m(spline,column,true);

    if(XX.bandsize()<=K.bandsize())
      prec.addto2(XX,K,scaleinv,1.0/sigma2);
    else
      prec.addto2(K,XX,1.0/sigma2,scaleinv);

    double * work = standnormal.getV();
    for(unsigned int i=0;i<nrpar;i++,work++){
      *work = rand_normal();
      }

    prec.solveL(standnormal,beta);
    likep->compute_respminuslinpred(mu,column);   // nicht ändern wegen multgaussian
    compute_XWtildey(likep->get_weight(),scaleinv);
    prec.solve(muy,betahelp,0,0);
    beta.plus(beta,betahelp);

    // hier bezüglich Haupteffekten zentrieren und diese "wegschmeißen"!
    compute_intercept();
    compute_beta();
    fcconst->update_intercept(intercept);
    for(unsigned int i=0;i<nrpar;i++){
      beta(i,0) = beta(i,0) + intercept;
      }

    add_linearpred_multBS(beta);

    for(unsigned int i=0;i<nrpar;i++){
      beta(i,0) -= intercept;
      }
    for(unsigned int i=0;i<spline.rows();i++){
      spline(i,0) -= intercept;
      }
    intercept = 0.0;

    acceptance++;
    }  // end: gauss

  if(center && !centertotal)   // Gesamteffekt in fctotal schreiben
      {
      if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
          ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
        {

        double * fctotalbetap = fctotal.getbetapointer();

        if(gridsize < 0)
          {
          vector<int>::iterator freqwork = freqoutput.begin();
          int * workindex = index.getV();
          for(unsigned int i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
            {
            if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
              {
              *fctotalbetap = spline(*workindex,0)
                            + mainp1->get_spline()(*workindex,0)
                            + mainp2->get_spline()(*workindex,0);
              fctotalbetap++;
              }
            }
          }
        else
          {
          multDG(splinehelp,beta);
          unsigned k,l;
          for(k=0;k<unsigned(gridsizex);k++)
            for(l=0;l<unsigned(gridsizey);l++,fctotalbetap++)
              *fctotalbetap = splinehelp(k*gridsizey + l,0)
                            + mainp1->get_splinehelp()(k,0)
                            + mainp2->get_splinehelp()(l,0);
          }
        }   // ENDE: Gesamteffekt in fctotal schreiben

      fctotal.update();
      } // END: interaction

// Interaktionseffekt in fchelp schreiben
  if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
    {
    double * fchelpbetap = fchelp.getbetapointer();

    if(gridsize < 0)
      {
      if(varcoeff)
        multBout(splinehelp,beta);

      vector<int>::iterator freqwork = freqoutput.begin();
      int * workindex = index.getV();
      for(unsigned int i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
        {
        if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
          {
          if(varcoeff)
            *fchelpbetap = splinehelp(i,0);
          else
            *fchelpbetap = spline(*workindex,0);
          fchelpbetap++;
          }
        }
      }
    else
      {
      multDG(splinehelp,beta);
      for(unsigned int i=0;i<unsigned(gridsize);i++,fchelpbetap++)
        *fchelpbetap = splinehelp(i,0);
      }
    }

  fchelp.update();
  FULLCOND::update();
  }


void FULLCOND_pspline_surf_stepwise::change_Korder(double lamb)
  {
  set_lambdaconst(1000000000);
  if(lamb==-1)
    {
    if(!varcoeff)
      {
      if(type==mrflinear || type==mrfkr1)
        {
        if(Kalt.bandsize()!=(nrpar1dim+1)*2)
          {
          Kalt = K;
          datamatrix Kstat = STATMAT_PENALTY::K2dim_pspline_rw2(nrpar1dim,2,2);
          datamatrix de(nrpar,1);
          datamatrix ud = datamatrix(nrpar,(nrpar1dim+1)*2);
          for(unsigned int i=0;i<nrpar;i++)
           {
           de(i,0) = Kstat(i,i);
           for(unsigned j=0;j<ud.cols();j++)
             {
             if (i+j+1 < nrpar)
               ud(i,j) = Kstat(i,i+j+1);
             }
           } // end: for(i=0;i<sizeK;i++)
          K = bandmatdouble(de,ud);
          if(utype!=gaussian)
            Kenv = envmatdouble(K);
          }
        else
          {
          bandmatdouble Khelp = K;
          K = Kalt;
          Kalt = Khelp;
          if(utype!=gaussian)
            Kenv = envmatdouble(K);
          }
        rankK = nrpar-2;
        }
      }
    else
      {
      if(type==mrfquadratic8 || type==mrfkr1)
        {
        if(Kalt.bandsize()!=nrpar1dim)
          {
          Kalt = K;
          Ksp = Kmrflinear(nrpar1dim,nrpar1dim);
          datamatrix de(nrpar,1);
          datamatrix ud = datamatrix(nrpar,nrpar1dim);
          for(unsigned i=0;i<nrpar;i++)
            {
            de(i,0) = Ksp(i,i);
            for(unsigned j=0;j<ud.cols();j++)
              {
              if (i+j+1 < nrpar)
                ud(i,j) = Ksp(i,i+j+1);
              }
            } // end: for(i=0;i<sizeK;i++)
          K = bandmatdouble(de,ud);
          if(utype!=gaussian)
            Kenv = envmatdouble(K);
          }
        else
          {
          bandmatdouble Khelp = K;
          K = Kalt;
          Kalt = Khelp;
          if(utype!=gaussian)
            Kenv = envmatdouble(K);
          }
        rankK = nrpar-1;
        }
      }
    }
  else if(lamb==-2)
    {
    if(type==mrflinear || type==mrfkr1)
      {
      if(Kalt.bandsize()!=(nrpar1dim+1)*2)
        {
        Kalt = K;
        datamatrix Kstat = STATMAT_PENALTY::K2dim_pspline_rw2(nrpar1dim,2,2);
        datamatrix de(nrpar,1);
        datamatrix ud = datamatrix(nrpar,(nrpar1dim+1)*2);
        for(unsigned i=0;i<nrpar;i++)
         {
         de(i,0) = Kstat(i,i);
         for(unsigned j=0;j<ud.cols();j++)
           {
           if (i+j+1 < nrpar)
             ud(i,j) = Kstat(i,i+j+1);
           }
         } // end: for(i=0;i<sizeK;i++)
        K = bandmatdouble(de,ud);
        if(utype!=gaussian)
          Kenv = envmatdouble(K);
        }
      else
        {
        bandmatdouble Khelp = K;
        K = Kalt;
        Kalt = Khelp;
        if(utype!=gaussian)
          Kenv = envmatdouble(K);
        }
      rankK = nrpar-2;
      }
    }
  }

void FULLCOND_pspline_surf_stepwise::undo_Korder(void)
  {
  if((type==mrflinear || type==mrfkr1) && rankK==nrpar-2)
    {
    bandmatdouble Khelp = K;
    K = Kalt;
    Kalt = Khelp;
    if(utype!=gaussian && type==mrflinear)
      Kenv = envmatdouble(K);
    else if(utype!=gaussian && type==mrfkr1)
      {
      double bandw = KHenv.getBandwidth();
      if(Kxenv.getBandwidth() > bandw)
        bandw = Kxenv.getBandwidth();
      if(Kyenv.getBandwidth() > bandw)
        bandw = Kyenv.getBandwidth();
      Kenv = envmatdouble(0.0,nrpar,bandw);
      }

    if(type==mrflinear)
      rankK = nrpar-1;
    else
      rankK = (nrpar1dim-1)*(nrpar1dim-1);
    }
  else if((type==mrfquadratic8 || type==mrfkr1) && rankK==nrpar-1)
    {
    bandmatdouble Khelp = K;
    K = Kalt;
    Kalt = Khelp;
    if(utype!=gaussian && type==mrfquadratic8)
      Kenv = envmatdouble(K);
    else if(utype!=gaussian && type==mrfkr1)
      {
      double bandw = KHenv.getBandwidth();
      if(Kxenv.getBandwidth() > bandw)
        bandw = Kxenv.getBandwidth();
      if(Kyenv.getBandwidth() > bandw)
        bandw = Kyenv.getBandwidth();
      Kenv = envmatdouble(0.0,nrpar,bandw);
      }

    if(type==mrfquadratic8)
      rankK = nrpar-2;
    else
      rankK = (nrpar1dim-1)*(nrpar1dim-1);
    }
  }


void FULLCOND_pspline_surf_stepwise::get_samples(const ST::string & filename,const unsigned & step) const
  {
  fchelp.get_samples(filename,step);
  }


void FULLCOND_pspline_surf_stepwise::compute_main(void)
  {
  unsigned i;
  double *workspline;
  int *workindex = index.getV();
  vector<int>::iterator freqwork;

// Haupteffekte berechnen
  he1.mult(betaweightx,beta);
  he2.mult(betaweighty,beta);

// 'spline' ändern
  freqwork = mainpoi1->get_freqit();
  workindex = mainpoi1->get_indexp();
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    spline(*workindex,0) -= he1(*freqwork,0);
  freqwork = mainpoi2->get_freqit();
  workindex = mainpoi2->get_indexp();
  for(i=0;i<spline.rows();i++,freqwork++,workindex++)
    spline(*workindex,0) -= he2(*freqwork,0);

  workspline = spline.getV();
  for(i=0;i<spline.rows();i++,workspline++)
    *workspline += intercept;
  }


void FULLCOND_pspline_surf_stepwise::safe_splines(bool & interact)
  {
  if(type==mrfkr1 && rankK==(nrpar1dim-1)*(nrpar1dim-1))
    {
    if(splineo1.rows()<spline.rows())
      {
      splineo1 = datamatrix(spline.rows(),1,0);
      splineo2 = datamatrix(spline.rows(),1,0);
      }
    splineo1.assign(mainpoi1->get_spline());
    splineo2.assign(mainpoi2->get_spline());
    interact = true;
    }
  else
    interact = false;
  }


void FULLCOND_pspline_surf_stepwise::set_splines_old(void)
  {
  likep->substr_linearpred_m(mainpoi1->get_spline(),column,true);
  likep->substr_linearpred_m(mainpoi2->get_spline(),column,true);
  mainpoi1->set_spline(splineo1);
  mainpoi2->set_spline(splineo2);
  likep->add_linearpred_m(mainpoi1->get_spline(),column,true);
  likep->add_linearpred_m(mainpoi2->get_spline(),column,true);
  }

void FULLCOND_pspline_surf_stepwise::compute_main_varcoeff(void)
  {
  unsigned i;
  double *workspline;
  int *workindex = index.getV();
  vector<int>::iterator freqwork;

// Haupteffekte berechnen
  he1.mult(betaweightx,beta);
  he2.mult(betaweighty,beta);

// 'splinehelp' ändern
  freqwork = mainpoi1->get_freqoutputit();
  workindex = mainpoi1->get_indexp();
  for(i=0;i<splinehelp.rows();i++,freqwork++,workindex++)
    splinehelp(*workindex,0) -= he1(*freqwork,0);
  freqwork = mainpoi2->get_freqoutputit();
  workindex = mainpoi2->get_indexp();
  for(i=0;i<splinehelp.rows();i++,freqwork++,workindex++)
    splinehelp(*workindex,0) -= he2(*freqwork,0);

  workspline = splinehelp.getV();
  if(center)
    {
    for(i=0;i<splinehelp.rows();i++,workspline++)
      *workspline += intercept;
    }
  else
    {
    for(i=0;i<splinehelp.rows();i++,workspline++)
      *workspline += 2*intercept;
    }

  workspline = spline.getV();
  double * worksph = splinehelp.getV();
  double * workint = data_forfixed.getV();
  for(i=0;i<spline.rows();i++,workspline++,worksph++,workint++)
      *workspline = *worksph * *workint;
  }


void FULLCOND_pspline_surf_stepwise::multBS_index(datamatrix & res, const datamatrix & b)
  {
  int i;
  unsigned j,k;
  double val=0.0;

  double *workbeta;
  double *workB;
  if(varcoeff)
    workB = Bout.getV();
  else
    workB = B.getV();
  int *workindex = index.getV();

  int maxfirst = first[res.rows()-1];

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator firstit = first.begin();

  i = 0;
  while(i<=maxfirst)
    {
    workbeta = b.getV() + i;
    while(*firstit==i)
      {
      if( freqwork==freq.begin() || *freqwork!=*(freqwork-1) )
        {
        val = 0.0;
        for(j=0;j<degree+1;j++)
          for(k=0;k<degree+1;k++,workB++)
            val += *workB * *(workbeta+k+j*nrpar1dim);
        }
      res(*workindex,0) = val;
      workindex++;
      freqwork++;
      firstit++;
      }
    i++;
    }
  }


} // end: namespace MCMC











