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

#include "fullcond_pspline_surf_gaussian.h"

namespace MCMC
{


double FULLCOND_pspline_surf_gaussian::compute_squareddiff(unsigned i,
                                         unsigned j,unsigned k, unsigned l,
                                         unsigned nr)
  {
  unsigned row1 = i*nr + j;
  unsigned row2 = k*nr + l;

  double diff = beta(row1,0)-beta(row2,0);

  return (diff*diff)/sigma2;

  }


void FULLCOND_pspline_surf_gaussian::compute_squareddiff(datamatrix & u)
  {

  double * worku = u.getV()+1;
  double * workbeta = beta.getV();
  double * workbeta2 = beta.getV()+1;
  double v;
  unsigned i;
  for (i=1;i<beta.rows();i++,worku++,workbeta++,workbeta2++)
    {
    v = *workbeta2-*workbeta;
    *worku = (v*v)/sigma2;
    }

/*
   ofstream out("c:\\tmp\\u.raw");
   u.prettyPrint(out);
   out.close();


   ofstream out2("c:\\tmp\\beta.raw");
   beta.prettyPrint(out2);
   out2.close();
*/

  }





void FULLCOND_pspline_surf_gaussian::init_maineffects(
                        spline_basis * mp1,spline_basis * mp2,
                        const ST::string & pnt,const ST::string & prt)
  {

  mainp1 = mp1;
  mainp2 = mp2;

  interaction = true;

  centertotal = false;

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

  }


void FULLCOND_pspline_surf_gaussian::create(const datamatrix & v1, const datamatrix & v2, const datamatrix & intact)
  {

  samplecentered = false;

  if(!singleblock || utype != gaussian)
    samplecentered = false;

  unsigned i=0,j=0,bands=0;

  lambda_prec = -1;

  make_index(v1,v2);

  datamatrix Kstat;

  if(type == mrflinear)
    {
    Ksp = Kmrflinear(nrpar1dim,nrpar1dim);
    rankK = nrpar-1;
    bands = 1;
    }
  else if(type == mrfquadratic8)          //  RW2 * RW2
    {
    Kstat = STATMAT_PENALTY::K2dim_pspline_rw2(nrpar1dim,2,2);
    rankK = nrpar-2;
    bands = 2;
    }
  else if(type == mrfkr1)          //  RW1 * RW1
    {
    Ksp = Krw1(vector<double>(nrpar1dim,1.0)).kronecker(Krw1(vector<double>(nrpar1dim,1.0)));
/*	SparseMatrix mathur;              // Rajan ???
	mathur = Krw1(vector<double>(nrpar1dim,1.0));
    Ksp = mathur.kronecker(mathur); */
    rankK = (nrpar1dim-1)*(nrpar1dim-1);
    bands = 1;
    }
  else if(type == mrfkr2)          //  RW2 * RW2
    {
    Ksp = Krw2(vector<double>(nrpar1dim,1.0)).kronecker(Krw2(vector<double>(nrpar1dim,1.0)));
    rankK = (nrpar1dim-2)*(nrpar1dim-2);
    bands = 2;
    }
  // NEU STEFAN
  else if(type == mrflinearband)
    {
    K = Kmrflinearband(nrpar1dim,nrpar1dim);
    rankK = nrpar1dim*nrpar1dim-nrpar1dim;
    bands=1;
    }
  // ENDE: NEU STEFAN

  if(type == mrflinear || type == mrfkr1 || type == mrfkr2)
    {
    datamatrix de(nrpar,1);
    datamatrix ud;
    if (type==mrflinear)
      ud = datamatrix(nrpar,(nrpar1dim)*bands);
    else
      ud = datamatrix(nrpar,(nrpar1dim+1)*bands);
    for(i=0;i<nrpar;i++)
      {
      de(i,0) = Ksp(i,i);
      for(j=0;j<ud.cols();j++)
        {
        if (i+j+1 < nrpar)
          ud(i,j) = Ksp(i,i+j+1);
        }
      } // end: for(i=0;i<sizeK;i++)

    K = bandmatdouble(de,ud);
    }
  else if(type == mrfquadratic8)
    {
    datamatrix de(nrpar,1);
    datamatrix ud;
    if (type==mrflinear)
      ud = datamatrix(nrpar,(nrpar1dim)*bands);
    else
      ud = datamatrix(nrpar,(nrpar1dim+1)*bands);
    for(i=0;i<nrpar;i++)
      {
      de(i,0) = Kstat(i,i);
      for(j=0;j<ud.cols();j++)
        {
        if (i+j+1 < nrpar)
          ud(i,j) = Kstat(i,i+j+1);
        }
      } // end: for(i=0;i<sizeK;i++)

    K = bandmatdouble(de,ud);
    }

  Kenv = envmatdouble(K);

  compute_knots(v1,v2);
  make_B(v1,v2);
  if(varcoeff)
    make_BVC(intact);

// fchelp initialisieren
  ST::string pnt = samplepath.substr(0,samplepath.length()-4)+"_fchelp.raw";
  if(gridsize < 0)
    {
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"help",nrdiffobs,1,pnt);
    splinehelp = datamatrix(likep->get_nrobs(),1,0);
    }
  else
    {
    fchelp = FULLCOND(optionsp,datamatrix(1,1,0),title+"help",gridsize,1,pnt);
    splinehelp = datamatrix(gridsize,1,0);
    }
  fchelp.setflags(MCMC::norelchange | MCMC::nooutput);
  fchelp.set_transform(transform);

  make_xy_v(v1,v2);

  if(gridsize > 0)
    {
    make_xy_values(v1,v2);
    make_DG();
    }

  W = datamatrix(likep->get_nrobs(),1,1.0);
  XX = bandmatdouble(nrpar,(nrpar1dim+1)*degree,0);
  XX_env = envmatdouble(0.0,nrpar,(nrpar1dim+1)*degree);
  compute_XWX(likep->get_weight());
  compute_XWXenv(likep->get_weight());

  if(degree > bands)
    {
    prec = bandmatdouble(nrpar,(nrpar1dim+1)*degree,0);
    prec_env = envmatdouble(0.0,nrpar,(nrpar1dim+1)*degree);
    }
  else
    {
    prec = bandmatdouble(nrpar,(nrpar1dim+1)*bands,0);
    prec_env = envmatdouble(0.0,nrpar,(nrpar1dim+1)*degree);
    }

  if(singleblock)
    {
    betahelp = datamatrix(nrpar,1,0);
    standnormal = datamatrix(nrpar,1,0);
    }
  else
    {
    betahelp = datamatrix(nrpar1dim,1,0);
    standnormal = datamatrix(nrpar1dim,1,0);
    betahelp2 = datamatrix(nrpar1dim,1,0);
    muyhelp = datamatrix(nrpar1dim,1,0);
    beta_ab = datamatrix(nrpar1dim,1,0);
    prop_ab = datamatrix(nrpar1dim,1,0);
    beta_mode_ab = datamatrix(nrpar1dim,1,0);    
    }  // end if(!singleblock)

  mu = datamatrix(likep->get_nrobs(),1,0);
  muy = datamatrix(nrpar,1,0);

  weight = vector<double>(nrpar,1.0/double(nrpar));

  compute_betaweight();
  compute_betaweightxy();

  if(varcoeff)
    {
    identifiable = true;
    interactvar=intact;
    }
  else
    identifiable = false;

  }

  // CONSTRUCTOR

FULLCOND_pspline_surf_gaussian::FULLCOND_pspline_surf_gaussian(MCMCoptions * o,
                      DISTRIBUTION * dp, FULLCOND_const * fcc,
                      const datamatrix & v1, const datamatrix & v2,
                      const ST::string & ti,
                      const unsigned & nrk, const unsigned & degr,
                      const knotpos & kp, const double & l, const int & gs,
                      const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const ST::string & of, const bool & sb, const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {

  utype = gaussian;

  varcoeff = false;

  singleblock = true;
  transform = likep->get_trmult(c);
  outfile = of;

  create(v1,v2);

  }

  // CONSTRUCTOR 2: IWLS

FULLCOND_pspline_surf_gaussian::FULLCOND_pspline_surf_gaussian(MCMCoptions * o,
                      DISTRIBUTION * dp, FULLCOND_const * fcc,
                      const datamatrix & v1, const datamatrix & v2,
                      const bool & mode, const ST::string & ti,
                      const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                      const double & l, const unsigned & upW, const bool & updatetau,
                      const double & fstart, const double & a, const double & b,
                      const int & gs, const fieldtype & ft,
                      const ST::string & fp, const ST::string & pres, const ST::string & of,
                      const bool & iw, const bool & sb,const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {


  if(mode)
    {
    if(updatetau)
      {
      utype = hyperblockmode;
      f2 = fstart;
      kappaburnin = datamatrix(o->get_burnin()/2.0,1,0.0);
      }
    else
      {
      utype = iwlsmode;
      }
    }
  else
    {
    if(updatetau)
      utype = hyperblock;
    else
      utype = iwls;
    }

  varcoeff = false;

  kappa = 1.0/sigma2;
  kappacurrent = 1.0/sigma2;

  singleblock = true;

  f = 2.0;
  updateW = upW;

  a_invgamma = a;
  b_invgamma = b;

  proposal = beta;
  outfile = of;

  create(v1,v2);

/*
// Xblock und Kblock initialisieren (für zeilenweise updaten)
  envmatdouble help;
  vector<double>::iterator di;

  help = Krw1env(vector<double>(nrpar1dim,1.0));

  di = help.getDiagIterator();
  for(unsigned i=0;i<nrpar1dim;i++,di++)
    (*di)++;
  Kblock.push_back(help);

  di = help.getDiagIterator();
  for(unsigned i=0;i<nrpar1dim;i++,di++)
    (*di)++;
  Kblock.push_back(help);

//  for(unsigned i=0;i<nrpar1dim;i++,di++)
//    Xblock.push_back();
*/
  }

  // CONSTRUCTOR 3: geosplines

FULLCOND_pspline_surf_gaussian::FULLCOND_pspline_surf_gaussian(
                      MCMCoptions * o, DISTRIBUTION * dp,
                      FULLCOND_const * fcc,
                      const datamatrix & region, const MAP::map & mp, const ST::string & mn,
                      const ST::string & ti, const unsigned & nrk,
                      const unsigned & degr, const knotpos & kp, const double & l,
                      const int & gs, const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const bool & sb, const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {

  utype = gaussian;

  varcoeff = false;

  singleblock = true;

  transform = likep->get_trmult(c);

  m = mp;
  mapexisting = true;
  mapname = mn;
  if(mp.polygones_existing() == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

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

  create(v1,v2);

  }

  // CONSTRUCTOR 4: IWLS geosplines

FULLCOND_pspline_surf_gaussian::FULLCOND_pspline_surf_gaussian(
                      MCMCoptions * o, DISTRIBUTION * dp, FULLCOND_const * fcc,
                      const datamatrix & region, const MAP::map & mp, const ST::string & mn,
                      const bool & mode, const ST::string & ti, const unsigned & nrk,
                      const unsigned & degr, const knotpos & kp, const double & l, const unsigned & upW,
                      const bool & updatetau, const double & fstart, const double & a, const double & b,
                      const int & gs, const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const bool & iw, const bool & sb, const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {

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

  varcoeff = false;

  singleblock = true;

  f = fstart;
  updateW = upW;

  a_invgamma = a;
  b_invgamma = b;

  proposal = beta;

  m = mp;
  mapexisting = true;
  mapname = mn;
  if(mp.polygones_existing() == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

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

  create(v1,v2);

  }

  // CONSTRUCTOR 5: varying coefficients

FULLCOND_pspline_surf_gaussian::FULLCOND_pspline_surf_gaussian(MCMCoptions * o,
                      DISTRIBUTION * dp, FULLCOND_const * fcc,
                      const datamatrix & intact,
                      const datamatrix & v1, const datamatrix & v2,
                      const ST::string & ti, const unsigned & nrk, const unsigned & degr,
                      const knotpos & kp, const double & l, const int & gs,
                      const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const ST::string & of,
                      const bool & sb,const bool & ce, const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {

  utype = gaussian;

  varcoeff = true;

  singleblock = true;

  transform = likep->get_trmult(c);
  outfile = of;

  create(v1,v2,intact);

  if (ce==false)
    identifiable = true;
  else
    identifiable = false;
  }

  // CONSTRUCTOR 6: IWLS varying coefficients

FULLCOND_pspline_surf_gaussian::FULLCOND_pspline_surf_gaussian(MCMCoptions * o,
                      DISTRIBUTION * dp, FULLCOND_const * fcc,const datamatrix & intact,
                      const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                      const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                      const double & l, const unsigned & upW, const bool & updatetau,
                      const double & fstart, const double & a, const double & b,
                      const int & gs, const fieldtype & ft,
                      const ST::string & fp, const ST::string & pres,
                      const ST::string & of,
                      const bool & iw, const bool & sb,
                      const bool & ce, const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {

/*  if(mode)
    {
    if(updatetau)
      utype = hyperblockmode;
    else
      utype = iwlsmode;
    }
  else
*/    {
    if(updatetau)
      utype = hyperblock;
    else
      utype = iwls;
    }

  varcoeff = true;

  kappa = 1.0/sigma2;
  kappacurrent = 1.0/sigma2;

  singleblock = true;

  f = fstart;
  updateW = upW;

  a_invgamma = a;
  b_invgamma = b;

  proposal = beta;
  outfile = of;

  create(v1,v2,intact);

  if (ce==false)
    identifiable = true;
  else
    identifiable = false;


  }

  // CONSTRUCTOR 7: geosplines varying coefficients

FULLCOND_pspline_surf_gaussian::FULLCOND_pspline_surf_gaussian(
                      MCMCoptions * o, DISTRIBUTION * dp, FULLCOND_const * fcc,
                      const datamatrix & intact,
                      const datamatrix & region, const MAP::map & mp, const ST::string & mn,
                      const ST::string & ti, const unsigned & nrk,
                      const unsigned & degr, const knotpos & kp, const double & l,
                      const int & gs, const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const bool & sb,
                      const bool & ce, const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {


  utype = gaussian;

  varcoeff = true;

  singleblock = true;

  transform = likep->get_trmult(c);

  m = mp;
  mapexisting = true;
  mapname = mn;
  if(mp.polygones_existing() == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

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

  create(v1,v2,intact);

  if (ce==false)
    identifiable = true;
  else
    identifiable = false;

  }

  // CONSTRUCTOR 8: IWLS geosplines varying coefficients

FULLCOND_pspline_surf_gaussian::FULLCOND_pspline_surf_gaussian(
                      MCMCoptions * o, DISTRIBUTION * dp,FULLCOND_const * fcc,
                      const datamatrix & intact,
                      const datamatrix & region, const MAP::map & mp, const ST::string & mn,
                      const ST::string & ti, const unsigned & nrk,
                      const unsigned & degr, const knotpos & kp, const double & l, const unsigned & upW,
                      const bool & updatetau, const double & fstart, const double & a, const double & b,
                      const int & gs, const fieldtype & ft, const ST::string & fp,
                      const ST::string & pres, const bool & iw, const bool & sb,
                      const bool & ce, const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {


/*  if(mode)
    {
    if(updatetau)
      utype = hyperblockmode;
    else
      utype = iwlsmode;
    }
  else
*/    {
    if(updatetau)
      utype = hyperblock;
    else
      utype = iwls;
    }

  varcoeff = true;

  singleblock = true;

  f = fstart;
  updateW = upW;

  a_invgamma = a;
  b_invgamma = b;

  proposal = beta;

  m = mp;
  mapexisting = true;
  mapname = mn;
  if(mp.polygones_existing() == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

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

  create(v1,v2,intact);

  if (ce==false)
    identifiable = true;
  else
    identifiable = false;


  }


FULLCOND_pspline_surf_gaussian::FULLCOND_pspline_surf_gaussian(const FULLCOND_pspline_surf_gaussian & fc)
  : spline_basis_surf(spline_basis_surf(fc))
  {

  lambda_prec = fc.lambda_prec;

  f2 = fc.f2;
  kappaburnin = fc.kappaburnin;

  utype = fc.utype;

  kappa = fc.kappa;
  kappaprop = fc.kappaprop;
  kappacurrent = fc.kappacurrent;
  kappamean = fc.kappamean;
  kappamode = fc.kappamode;
  kappavar = fc.kappavar;

  samplecentered = fc.samplecentered;

  updateW = fc.updateW;

  a_invgamma = fc.a_invgamma;
  b_invgamma = fc.b_invgamma;

  beta_ab = fc.beta_ab;
  prop_ab = fc.prop_ab;
  beta_mode_ab = fc.beta_mode_ab;

  W = fc.W;
  proposal = fc.proposal;
  prec2 = fc.prec2;

  XX = fc.XX;
  XX_env = fc.XX_env;
  prec_env = fc.prec_env;

  Xblock = fc.Xblock;
  Kblock = fc.Kblock;

  singleblock = fc.singleblock;

  prec = fc.prec;
  mu = fc.mu;
  muy = fc.muy;
  muyhelp = fc.muyhelp;
  betahelp = fc.betahelp;
  betahelp2 = fc.betahelp2;
  standnormal = fc.standnormal;

  }


const FULLCOND_pspline_surf_gaussian & FULLCOND_pspline_surf_gaussian::operator=(
                                            const FULLCOND_pspline_surf_gaussian & fc)
  {
  if (this == &fc)
    return *this;
  spline_basis_surf::operator=(spline_basis_surf(fc));

  lambda_prec = fc.lambda_prec;

  f2 = fc.f2;
  kappaburnin = fc.kappaburnin;

  utype = fc.utype;

  kappa = fc.kappa;
  kappaprop = fc.kappaprop;
  kappacurrent = fc.kappacurrent;
  kappamean = fc.kappamean;
  kappamode = fc.kappamode;
  kappavar = fc.kappavar;

  samplecentered = fc.samplecentered;

  updateW = fc.updateW;

  a_invgamma = fc.a_invgamma;
  b_invgamma = fc.b_invgamma;

  beta_ab = fc.beta_ab;
  prop_ab = fc.prop_ab;
  beta_mode_ab = fc.beta_mode_ab;

  W = fc.W;
  proposal = fc.proposal;
  prec2 = fc.prec2;

  XX = fc.XX;
  XX_env = fc.XX_env;
  prec_env = fc.prec_env;

  Xblock = fc.Xblock;
  Kblock = fc.Kblock;

  singleblock = fc.singleblock;

  prec = fc.prec;
  mu = fc.mu;
  muy = fc.muy;
  muyhelp = fc.muyhelp;
  betahelp = fc.betahelp;
  betahelp2 = fc.betahelp2;
  standnormal = fc.standnormal;

  return *this;
  }


void FULLCOND_pspline_surf_gaussian::update_IWLS(void)
  {

  double iwlsscale = 1.0;

  if(singleblock)
    {

    double logold = - 0.5*Kenv.compute_quadform(beta,0)/sigma2;

    if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
        ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
      {
      logold += likep->compute_IWLS(W,mu,true,column,true);
      mu.plus(spline,mu);
      compute_XWXenv(W,column);
      }
    else
      {
      logold += likep->compute_IWLS(W,mu,false,column,true);
      mu.plus(mu,spline);
      }

    compute_XWtildey(W,1.0);

    prec_env.addto(XX_env,Kenv,iwlsscale,iwlsscale/sigma2);

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
      compute_XWXenv(W,column);
      prec_env.addto(XX_env,Kenv,iwlsscale,iwlsscale/sigma2);
      }
    else
      {
      lognew += likep->compute_IWLS(W,mu,false,column,true);
      mu.plus(mu,spline);
      }

    compute_XWtildey(W,1.0);

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
  else  // rowwise
    {
/*
    unsigned beg,end;
    vector<int>::iterator firstit;

    double * work;
    double * workbeta;
    unsigned i,j;

    unsigned an = 0;
    unsigned en = nrpar1dim;

    if(iwlsmode)
      {
      add_linearpred_multBS2(beta_mode);
      likep->compute_weight(W,column,true);
      compute_XWX(W);
      likep->tilde_y_minus_eta(mu,column,true);
      compute_XWtildey(W,1.0);
      add_linearpred_multBS2(beta);
      }

    for(i=0;i<nrpar1dim;i++,an+=nrpar1dim,en+=nrpar1dim)
      {

      beg = 0;
      firstit = first.begin();
      while(firstit != first.end() && *firstit+nrpar1dim*degree+(degree+1) < an+1)
        {
        beg++;
        firstit++;
        }
      end = beg;
      while(firstit != first.end() && end < likep->get_nrobs() && *firstit < en)
        {
        end++;
        firstit++;
        }
      if(end > 0)
        end--;

      nrtrials++;

      if(iwlsmode)
      {

      prec.addtoblock2(XX,K,1.0,1.0/sigma2,an,en-1);
      prec2 = prec;

      work = standnormal.getV();
      for(j=0;j<nrpar1dim;j++,work++)
        *work = rand_normal();

      XX.multBlock(beta_mode,muyhelp,an,an,en-1,en-1,an);

      work = muy.getV() + an;
      for(j=0;j<nrpar1dim;j++,++work)
        muyhelp(j,0) += *work;

      if(an>0)
        {
        work = beta_mode.getV()+an-nrpar1dim;
        for(j=0;j<nrpar1dim;j++,++work)
          muyhelp(j,0) += (*work/sigma2);
        }
      if(en<nrpar)
        {
        work = beta_mode.getV()+en;
        for(j=0;j<nrpar1dim;j++,++work)
          muyhelp(j,0) += (*work/sigma2);
        }

      prec.solve(muyhelp,beta_mode_ab,0,0);

      for(j=0;j<nrpar1dim;j++)
        beta_mode(j+an,0) = beta_mode_ab(j,0);

      prec.solveL(standnormal,prop_ab);

      prop_ab.plus(prop_ab,beta_mode_ab);

      proposal.assign(beta);
      for(j=0;j<nrpar1dim;j++)
        {
        proposal(an+j,0) = prop_ab(j,0);
        beta_ab(j,0) = beta(an+j,0);
        }

      betahelp.minus(prop_ab,beta_mode_ab);
      double qold = - 0.5*prec2.compute_quadform(betahelp,0);
      double logold = likep->loglikelihood(beg,end,index) - 0.5*K.compute_quadform(beta,0)/sigma2;

      betahelp.minus(prop_ab,beta_ab);
      add_linearpred_multBS_Block(betahelp,an,en-1,beg,end);

      betahelp.minus(beta_ab,beta_mode_ab);
      double qnew = - 0.5*prec2.compute_quadform(betahelp,0);
      double lognew = likep->loglikelihood(beg,end,index,true) - 0.5*K.compute_quadform(proposal,0)/sigma2;

      double alpha = lognew + qnew - logold - qold;
      double u = log(uniform());

      if (u <= alpha)
        {
        workbeta = beta.getV()+an;          // change beta
        for(j=0;j<nrpar1dim;j++,workbeta++)
          *workbeta = prop_ab(j,0);
        acceptance++;
        }
      else
        {
        betahelp.minus(beta_ab,prop_ab);
        add_linearpred_multBS_Block(betahelp,an,en-1,beg,end);
        }

      }
      else
      {

      likep->compute_weight(W,beg,end,index,column);
      compute_XWX_Block(W,an,en-1,beg,end);

      compute_q(beta,an,en,beg,end,sigma2,true);

      double * work = standnormal.getV();
      for(j=0;j<nrpar1dim;j++,work++)
        *work = rand_normal();

      prec.solveL(standnormal,prop_ab);

      prop_ab.plus(prop_ab,betahelp);

      proposal.assign(beta);
      for(j=0;j<nrpar1dim;j++)
        {
        proposal(an+j,0) = prop_ab(j,0);
        beta_ab(j,0) = beta(an+j,0);
        }

      betahelp.minus(prop_ab,betahelp);
      double qold = 0.5*prec.get_det() - 0.5*prec2.compute_quadform(betahelp,0);
      double logold = likep->loglikelihood(beg,end,index) - 0.5*K.compute_quadform(beta,0)/sigma2;

      betahelp.minus(prop_ab,beta_ab);
      add_linearpred_multBS_Block(betahelp,an,en-1,beg,end);

      likep->compute_weight(W,beg,end,index,column);
      compute_XWX_Block(W,an,en-1,beg,end);

      compute_q(proposal,an,en,beg,end,sigma2,true);

      betahelp.minus(beta_ab,betahelp);
      double qnew = 0.5*prec.get_det() - 0.5*prec2.compute_quadform(betahelp,0);
      double lognew = likep->loglikelihood(beg,end,index,true) - 0.5*K.compute_quadform(proposal,0)/sigma2;

      double alpha = lognew + qnew - logold - qold;
      double u = log(uniform());

      if (u <= alpha)
        {
        workbeta = beta.getV()+an;          // change beta
        for(j=0;j<nrpar1dim;j++,workbeta++)
          *workbeta = prop_ab(j,0);
        acceptance++;
        }
      else
        {
        betahelp.minus(beta_ab,prop_ab);
        add_linearpred_multBS_Block(betahelp,an,en-1,beg,end);
        }

      }

      }   // end:    for(j=0;j<nrpar1dim;j++,an+=nrpar1dim,en+=nrpar1dim)

    if(center)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta_mode(i,0) -= intercept;
      }
*/
    }  // end: rowwise

  }


void FULLCOND_pspline_surf_gaussian::update_IWLS_mode(void)
  {

  double iwlsscale = 1.0;

  if(singleblock)
    {

    double logold = likep->loglikelihood(true) - 0.5*Kenv.compute_quadform(beta,0)/sigma2;

    add_linearpred_multBS2(beta_mode);

    if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
        ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
      {
      likep->compute_IWLS_weight_tildey(W,mu,column,true);
      mu.plus(spline,mu);
      compute_XWXenv(W,column);
      }
    else
      {
      likep->tilde_y(mu,spline,column,true,W);
      }

    compute_XWtildey(W,1.0);

    prec_env.addto(XX_env,Kenv,iwlsscale,iwlsscale/sigma2);

    double * work = proposal.getV();
    for(unsigned i=0;i<nrpar;i++,work++)
      *work = rand_normal();

    prec_env.solve(muy,betahelp);
    prec_env.solveU(proposal,betahelp);

    add_linearpred_multBS2(proposal);
    beta_mode.assign(betahelp);

    betahelp.minus(proposal,beta_mode);
    double qold = - 0.5*prec_env.compute_quadform(betahelp,0);

    double lognew = likep->loglikelihood(true) - 0.5*Kenv.compute_quadform(proposal,0)/sigma2;

    betahelp.minus(beta,beta_mode);
    double qnew = - 0.5*prec_env.compute_quadform(betahelp,0);

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
  else  // rowwise
    {
/*
    unsigned beg,end;
    vector<int>::iterator firstit;

    double * workbeta;
    unsigned i,j;

    unsigned an = 0;
    unsigned en = nrpar1dim;

    for(i=0;i<nrpar1dim;i++,an+=nrpar1dim,en+=nrpar1dim)
      {

      beg = 0;
      firstit = first.begin();
      while(firstit != first.end() && *firstit+nrpar1dim*degree+(degree+1) < an+1)
        {
        beg++;
        firstit++;
        }
      end = beg;
      while(firstit != first.end() && end < likep->get_nrobs() && *firstit < en)
        {
        end++;
        firstit++;
        }
      if(end > 0)
        end--;

      nrtrials++;

      compute_q(beta,an,en,beg,end,sigma2,true);

      double * work = standnormal.getV();
      for(j=0;j<nrpar1dim;j++,work++)
        *work = rand_normal();

      prec.solveL(standnormal,prop_ab);

      prop_ab.plus(prop_ab,betahelp);

      proposal.assign(beta);
      for(j=0;j<nrpar1dim;j++)
        {
        proposal(an+j,0) = prop_ab(j,0);
        beta_ab(j,0) = beta(an+j,0);
        }

      betahelp.minus(prop_ab,betahelp);
      double qold = - 0.5*prec2.compute_quadform(betahelp,0);
      double logold = likep->loglikelihood(beg,end,index) - 0.5*K.compute_quadform(beta,0)/sigma2;

      betahelp.minus(prop_ab,beta_ab);
      add_linearpred_multBS_Block(betahelp,an,en-1,beg,end);

      compute_q(proposal,an,en,beg,end,sigma2,true);

      betahelp.minus(beta_ab,betahelp);
      double qnew = - 0.5*prec2.compute_quadform(betahelp,0);
      double lognew = likep->loglikelihood(beg,end,index,true) - 0.5*K.compute_quadform(proposal,0)/sigma2;

      double alpha = lognew + qnew - logold - qold;
      double u = log(uniform());

      if (u <= alpha)
        {
        workbeta = beta.getV()+an;          // change beta
        for(j=0;j<nrpar1dim;j++,workbeta++)
          *workbeta = prop_ab(j,0);
        acceptance++;
        }
      else
        {
        betahelp.minus(beta_ab,prop_ab);
        add_linearpred_multBS_Block(betahelp,an,en-1,beg,end);
        }

      }   // end:    for(j=0;j<nrpar1dim;j++,an+=nrpar1dim,en+=nrpar1dim)
*/
    }  // end: rowwise

  if(center)
    {
    compute_intercept();
    for(unsigned i=0;i<nrpar;i++)
      beta_mode(i,0) -= intercept;
    intercept = 0.0;
    }

  }


void FULLCOND_pspline_surf_gaussian::update_IWLS_hyperblock()
  {

  if(optionsp->get_nriter()%100==0 && optionsp->get_nriter()<optionsp->get_burnin())
    tune_updatetau(alpha_30);

  if(f<1.1)
    f = 1.1;

  if(optionsp->get_nriter() == optionsp->get_burnin())
    optionsp->out("NOTE: Tuning constant 'f' for term " + title + " set to " + ST::doubletostring(f) + "\n");

  kappaprop = kappa*randnumbers::rand_variance(f);                             // kappa ~ (1+1/z)

  if(singleblock)
    {

    double logold = - 0.5*Kenv.compute_quadform(beta,0)*kappa;
    logold += 0.5*rankK*log(kappa);             // Normalisierungs Konstante
    logold += (a_invgamma-1)*log(kappa) - b_invgamma*kappa;           // gamma prior

    if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
        ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
      {
      logold += likep->compute_IWLS(W,mu,true,column,true);
      mu.plus(spline,mu);
      compute_XWXenv(W,column);
      }
    else
      {
      logold += likep->compute_IWLS(W,mu,false,column,true);
      mu.plus(spline,mu);
      }

    compute_XWtildey(W,1.0);

    prec_env.addto(XX_env,Kenv,1.0,kappaprop);

    double * work = proposal.getV();
    for(unsigned i=0;i<nrpar;i++,work++)
      *work = rand_normal();

    prec_env.solve(muy,betahelp);
    prec_env.solveU(proposal,betahelp);

    add_linearpred_multBS2(proposal);

    betahelp.minus(proposal,betahelp);
    double qold = 0.5*prec_env.getLogDet() - 0.5*prec_env.compute_quadform(betahelp,0);

    double lognew = - 0.5*Kenv.compute_quadform(proposal,0)*kappaprop;
    lognew += 0.5*rankK*log(kappaprop);
    lognew += (a_invgamma-1)*log(kappaprop) - b_invgamma*kappaprop;

    if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
        ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
      {
      lognew += likep->compute_IWLS(W,mu,true,column,true);
      mu.plus(spline,mu);
      compute_XWXenv(W,column);
      }
    else
      {
      lognew += likep->compute_IWLS(W,mu,false,column,true);
      mu.plus(spline,mu);
      }

    compute_XWtildey(W,1.0);

    prec_env.addto(XX_env,Kenv,1.0,kappa);

    prec_env.solve(muy,betahelp);

    betahelp.minus(beta,betahelp);
    double qnew = 0.5*prec_env.getLogDet() - 0.5*prec_env.compute_quadform(betahelp,0);

    double alpha = lognew + qnew - logold - qold;

    double u = log(uniform());

    if(u<=alpha)
      {
      kappa = kappaprop;
      sigma2 = 1.0/kappa;
      acceptance++;
      beta.assign(proposal);
      }
    else
      {
      add_linearpred_multBS2(beta);
      }

    }
  else   // rowwise
    {
/*
    double sigma2mean = 0.0;

    unsigned beg,end;
    vector<int>::iterator firstit;

    double * workbeta;
    unsigned i,j;

    unsigned an = 0;
    unsigned en = nrpar1dim;

    for(i=0;i<nrpar1dim;i++,an+=nrpar1dim,en+=nrpar1dim)
      {

      beg = 0;
      firstit = first.begin();
      while(firstit != first.end() && *firstit+nrpar1dim*degree+(degree+1) < an+1)
        {
        beg++;
        firstit++;
        }
      end = beg;
      while(firstit != first.end() && end < likep->get_nrobs() && *firstit < en)
        {
        end++;
        firstit++;
        }
      if(end > 0)
        end--;

      nrtrials++;

      kappaprop = kappacurrent*randnumbers::rand_variance(f);

      likep->compute_weight(W,column,true);
      compute_XWX_Block(W,an,en-1,beg,end);

      compute_q(beta,an,en,beg,end,1.0/kappaprop,true);

      double * work = standnormal.getV();
      for(j=0;j<nrpar1dim;j++,work++)
        *work = rand_normal();

      prec.solveL(standnormal,prop_ab);

      prop_ab.plus(prop_ab,betahelp);

      proposal.assign(beta);
      for(j=0;j<nrpar1dim;j++)
        {
        proposal(an+j,0) = prop_ab(j,0);
        beta_ab(j,0) = beta(an+j,0);
        }

      betahelp.minus(prop_ab,betahelp);
      double qold = 0.5*prec.get_det() - 0.5*prec2.compute_quadform(betahelp,0);
      double logold = likep->loglikelihood(beg,end,index) - 0.5*K.compute_quadform(beta,0)*kappacurrent;
      logold += 0.5*rankK*log(kappacurrent);             // Normalisierungs Konstante
      logold += (a_invgamma-1)*log(kappacurrent) - b_invgamma*kappacurrent;           // gamma prior

      betahelp.minus(prop_ab,beta_ab);
      add_linearpred_multBS_Block(betahelp,an,en-1,beg,end);

      likep->compute_weight(W,column,true);
      compute_XWX_Block(W,an,en-1,beg,end);

      compute_q(proposal,an,en,beg,end,1.0/kappacurrent,true);

      betahelp.minus(beta_ab,betahelp);
      double qnew = 0.5*prec.get_det() - 0.5*prec2.compute_quadform(betahelp,0);
      double lognew = likep->loglikelihood(beg,end,index,true) - 0.5*K.compute_quadform(proposal,0)*kappaprop;
      lognew += 0.5*rankK*log(kappaprop);             // Normalisierungs Konstante
      lognew += (a_invgamma-1)*log(kappaprop) - b_invgamma*kappaprop;           // gamma prior

      double alpha = lognew + qnew - logold - qold;
      double u = log(uniform());

      if (u <= alpha)
        {
        kappacurrent = kappaprop;
        workbeta = beta.getV()+an;          // change beta
        for(j=0;j<nrpar1dim;j++,workbeta++)
          *workbeta = prop_ab(j,0);
        acceptance++;
        }
      else
        {
        betahelp.minus(beta_ab,prop_ab);
        add_linearpred_multBS_Block(betahelp,an,en-1,beg,end);
        }

      sigma2mean += (1.0/kappacurrent)/nrpar1dim;

      }

    sigma2 = sigma2mean;
*/
    }   // END: rowwise

  if(utype == hyperblockmode && optionsp->get_nriter() < optionsp->get_burnin()/2 )
    kappaburnin(optionsp->get_nriter()-1,0) = kappa;

  }


void FULLCOND_pspline_surf_gaussian::update_IWLS_hyperblock_mode(void)
  {

  unsigned i;

  double aprop,bprop;

// Tuning von f bei proposal nach Rue/Held
//  if(optionsp->get_nriter()%100==0 && optionsp->get_nriter()<optionsp->get_burnin()/2)
//    tune_updatetau(alpha_30);

//  if(optionsp->get_nriter() == optionsp->get_burnin())
//    optionsp->out("NOTE: Tuning constant 'f' for term " + title + " set to " + ST::doubletostring(f) + "\n");

// für aprop,bprop (Gamma proposal)
  if(optionsp->get_nriter() == optionsp->get_burnin()/2)
    {
    kappamean = 2.0*kappaburnin.mean(0);
    kappavar = kappaburnin.var(0);
//    optionsp->out("  mu for the Gamma proposal:  " + ST::doubletostring(kappamean) + "\n");
//    optionsp->out("  var for the Gamma proposal: " + ST::doubletostring(kappavar)  + "\n");
    }

// für Gamma proposal
  aprop = kappamean*kappamean/(f2*kappavar);
  bprop = kappamean/(f2*kappavar);
  kappamode = (aprop-1)/bprop;                                        // Mode Gamma proposal

  sigma2 = 1.0/kappamode;                                             // für prec_env!

//  kappaprop = kappamode*randnumbers::rand_variance(f);              // Rue/Held
  kappaprop = randnumbers::rand_gamma(aprop,bprop);                   // Gamma proposal

  double logold = likep->loglikelihood(true) - 0.5*Kenv.compute_quadform(beta,0)*kappa;
  logold += 0.5*rankK*log(kappa);                                     // Normalisierungskonstante der Priori von beta ( exp(-1/2 kappa * beta'K beta) )
  logold += (a_invgamma-1)*log(kappa) - b_invgamma*kappa;             // Gamma prior

  add_linearpred_multBS2(beta_mode);

  if( (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) ) )
    {
    likep->compute_IWLS_weight_tildey(W,mu,column,true);
    mu.plus(spline,mu);
    compute_XWXenv(W,column);
    }
  else
    {
    likep->tilde_y(mu,spline,column,true,W);
    }

  compute_XWtildey(W,1.0);

  prec_env.addto(XX_env,Kenv,1.0,1.0/sigma2);

  prec_env.solve(muy,betahelp);

  double * work = proposal.getV();
  for(unsigned i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  prec_env.addto(XX_env,Kenv,1.0,kappaprop);

  prec_env.solveU(proposal,betahelp);

  add_linearpred_multBS2(proposal);
  beta_mode.assign(betahelp);

  betahelp.minus(proposal,beta_mode);
  double qold = 0.5*prec_env.getLogDet() - 0.5*prec_env.compute_quadform(betahelp,0);
//  qold += log(1.0/kappamode + 1.0/kappaprop);                         // Rue/Held
  qold += (aprop-1)*log(kappaprop) - bprop*kappaprop;                 // Gamma proposal

  double lognew = likep->loglikelihood(true) - 0.5*Kenv.compute_quadform(proposal,0)*kappaprop;
  lognew += 0.5*rankK*log(kappaprop);                                 // Normalisierungskonstante der Priori von beta ( exp(-1/2 kappa * beta'K beta) )
  lognew += (a_invgamma-1)*log(kappaprop) - b_invgamma*kappaprop;     // Gamma prior

  prec_env.addto(XX_env,Kenv,1.0,kappa);

  betahelp.minus(beta,beta_mode);
  double qnew = 0.5*prec_env.getLogDet() - 0.5*prec_env.compute_quadform(betahelp,0);
//      qnew += log(1.0/kappamode + 1.0/kappa);                         // Rue/Held
  qnew += (aprop-1)*log(kappa) - bprop*kappa;                         // Gamma proposal

  double alpha = lognew + qnew - logold - qold;

  double u = log(uniform());

  if(u<=alpha)
    {
    kappa = kappaprop;
    sigma2 = 1.0/kappa;
    acceptance++;
    beta.assign(proposal);
    }
  else
    {
    add_linearpred_multBS2(beta);
    }

  if(center)
    {
    compute_intercept();
    for(i=0;i<nrpar;i++)
      beta_mode(i,0) -= intercept;
    intercept = 0.0;
    }

  sigma2 = 1.0/kappa;

  }


void FULLCOND_pspline_surf_gaussian::update(void)
  {

  if(lambdaconst == true)
    sigma2 = likep->get_scale(column)/lambda;

  if( !singleblock && optionsp->get_nriter() == 1)
    {
    unsigned bands;
    if(type == mrfkr2)
      bands = 2;
    else
      bands = 1;

    if(degree > bands)
      {
      prec = bandmatdouble(nrpar1dim,degree,0);
      prec_env = envmatdouble(0.0,nrpar1dim,degree);
      }
    else
      {
      prec = bandmatdouble(nrpar1dim,bands,0);
      prec_env = envmatdouble(0.0,nrpar1dim,degree);
      }
    }

  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);
  fctotal.set_transform(transform);

  unsigned i;

  if(utype == iwls)
    {
    update_IWLS();
    }
  else if(utype == iwlsmode)
    {
    update_IWLS_mode();
    }
  else if(utype == hyperblock)
    {
    update_IWLS_hyperblock();
    }
  else if(utype == hyperblockmode)
    {
    if(optionsp->get_nriter()<optionsp->get_burnin()/2)
      update_IWLS_hyperblock();                          // Startwert für kappa ermitteln: getrennt updaten
    else
      update_IWLS_hyperblock_mode();                     // Simulation: gemeinsam updaten

//    update_IWLS_hyperblock_mode();
    }
  else if(utype == gaussian)
    {

    double scaleinv = 1.0/likep->get_scale(column);

    if(changingweight || optionsp->get_nriter()==1)
      compute_XWX(likep->get_weight());

    if(singleblock)
      {

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

      }
    else    // rowwise
      {
/*
      unsigned beg,end;
      vector<int>::iterator firstit;

      double * work;
      double * workbeta;
      unsigned j;

      unsigned an = 0;
      unsigned en = nrpar1dim;

      for(i=0;i<nrpar1dim;i++,an+=nrpar1dim,en+=nrpar1dim)
        {

        beg = 0;
        firstit = first.begin();
        while(firstit != first.end() && *firstit+nrpar1dim*degree+(degree+1) < an+1)
          {
          beg++;
          firstit++;
          }
        end = beg;
        while(firstit != first.end() && end < likep->get_nrobs() && *firstit < en)
          {
          end++;
          firstit++;
          }
        if(end > 0)
          end--;


        likep->compute_respminuslinpred(mu,column);   // nicht ändern wegen multgaussian
        compute_XWtildey_Block(likep->get_weight(),scaleinv,beg,end);        // --> muy

        prec.addtoblock2(XX,K,1.0/likep->get_scale(column),1.0/sigma2,an,en-1);

        work = standnormal.getV();
        for(j=0;j<nrpar1dim;j++,work++)
          *work = rand_normal();

        prec.solveL(standnormal,betahelp2);

        XX.multBlock(beta,muyhelp,an,an,en-1,en-1,an);

        work = muyhelp.getV();
        workbeta = muy.getV()+an;
        for(j=0;j<nrpar1dim;j++,++work,++workbeta)
          *work = scaleinv * *work + *workbeta;

        if(i>0)
          {
          work = muyhelp.getV();
          workbeta = beta.getV()+an-nrpar1dim;
          for(j=0;j<nrpar1dim;j++,++work,++workbeta)
            *work += (*workbeta/sigma2);
          }
        if(i<nrpar1dim-1)
          {
          work = muyhelp.getV();
          workbeta = beta.getV()+en;
          for(j=0;j<nrpar1dim;j++,++work,++workbeta)
            *work += (*workbeta/sigma2);
          }

        prec.solve(muyhelp,betahelp,0,0);
        betahelp.plus(betahelp,betahelp2);

        work = betahelp.getV();
        workbeta = beta.getV()+an;
        for(j=0;j<nrpar1dim;j++,++work,++workbeta)
          {
          beta_ab(j,0) = *work - *workbeta;
          *workbeta = *work;
          }

        add_linearpred_multBS_Block(beta_ab,an,en-1,beg,end);

        } // end: for(i=0;i<nrpar1dim;i++)
*/
      }  // end: blocked

    acceptance++;

    }  // end: gauss

  if(!singleblock)
    multBS_index(spline,beta);

  if(center)
    {
    if(centertotal)
      {
      if(!samplecentered)
        {
        compute_intercept();
        if (varcoeff)
          {
          for(i=0;i<nrpar;i++)
            beta(i,0) -= intercept;

          if (center)
            {
            for(i=0;i<likep->get_nrobs();i++)
              spline(i,0) -= intercept*interactvar(i,0);
            }
          else
            {
            for(i=0;i<likep->get_nrobs();i++)
              spline(i,0) -= intercept;
            }

          fcconst->update_fix_varcoeff(intercept,datanames[1]);

          intercept = 0.0;
          }
        else
          {
          for(i=0;i<nrpar;i++)
            beta(i,0) -= intercept;
          for(i=0;i<likep->get_nrobs();i++)
            spline(i,0) -= intercept;
          fcconst->update_intercept(intercept);
          intercept = 0.0;
          }
        }
      }
    else
      {

      beta_uncentered.assign(beta);

      if(utype != gaussian)
        {
        compute_intercept();
        compute_main();
        compute_beta();
        fcconst->update_intercept(intercept);
        mainp1->change(beta1);
        mainp2->change(beta2);
        intercept = 0.0;
        }
      else
        {
        compute_intercept();
        compute_main();
        compute_beta();
        fcconst->update_intercept(intercept);
        mainp1->change(he1,intercept);
        mainp2->change(he2,intercept);
        intercept = 0.0;
        }

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


void FULLCOND_pspline_surf_gaussian::compute_q(const datamatrix & b, const unsigned & an, const unsigned & en,
                                               const unsigned & beg, const unsigned & end, const double & sigma2,
                                               const bool & current)
  {

  unsigned j;

  double * work;
  double * workmuyhelp;

  prec.addtoblock2(XX,K,1.0,1.0/sigma2,an,en-1);
  prec2 = prec;                       // wegen decomposedonly !!!

  likep->tilde_y_minus_eta(mu,column,current);
  compute_XWtildey_Block(W,1.0,beg,end);

  XX.multBlock(b,muyhelp,an,an,en-1,en-1,an);

  workmuyhelp = muyhelp.getV();
  work = muy.getV()+an;
  for(j=0;j<nrpar1dim;j++,++workmuyhelp,++work)
    *workmuyhelp += *work;

  if(an>0)
    {
    workmuyhelp = muyhelp.getV();
    work = b.getV()+an-nrpar1dim;
    for(j=0;j<nrpar1dim;j++,++workmuyhelp,++work)
      *workmuyhelp += (*work/sigma2);
    }
  if(en<nrpar)
    {
    workmuyhelp = muyhelp.getV();
    work = b.getV()+en;
    for(j=0;j<nrpar1dim;j++,++workmuyhelp,++work)
      *workmuyhelp += (*work/sigma2);
    }

  prec.solve(muyhelp,betahelp,0,0);

  }


bool FULLCOND_pspline_surf_gaussian::posteriormode_converged(const unsigned & itnr)
  {
  return likep->posteriormode_converged_fc(beta,beta_mode,itnr);
  }


bool FULLCOND_pspline_surf_gaussian::posteriormode(void)
  {
  bool converged = false;
  bool converged1 = false;
  bool converged2 = false;

  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);
  fctotal.set_transform(transform);

  unsigned i;

  if(utype == gaussian)
    likep->substr_linearpred_m(spline,column,true);

  compute_XWXenv(likep->get_weightiwls(),column);

  prec_env.addto(XX_env,Kenv,1.0,lambda);

  if(utype != gaussian)
    likep->substr_linearpred_m(spline,column,true);

  likep->compute_workingresiduals(column);
  compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);

  prec_env.solve(muy,beta);

  add_linearpred_multBS(beta);

  if(center)
    {
    if(centertotal)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;


      if (varcoeff)
        {
        if (center)
          {
          for(i=0;i<likep->get_nrobs();i++)
            spline(i,0) -= intercept*interactvar(i,0);
          }
        else
          {
          for(i=0;i<likep->get_nrobs();i++)
            spline(i,0) -= intercept;
          }

        fcconst->update_fix_varcoeff(intercept,datanames[1]);
        }
      else
        {
        for(i=0;i<likep->get_nrobs();i++)
          spline(i,0) -= intercept;
        fcconst->posteriormode_intercept(intercept);
        }

      intercept = 0.0;
      }
    else
      {

      beta_uncentered.assign(beta);

      if(utype != gaussian)
        {
        compute_intercept();
        compute_main();
        compute_beta();
        fcconst->posteriormode_intercept(intercept);
        converged1 = mainp1->changeposterior(beta1);
        converged2 = mainp2->changeposterior(beta2);
        intercept = 0.0;
        }
      else
        {
        compute_intercept();
        compute_main();
        compute_beta();
        fcconst->posteriormode_intercept(intercept);
        converged1 = mainp1->changeposterior(he1,intercept);
        converged2 = mainp2->changeposterior(he2,intercept);
        intercept = 0.0;
        }

// Gesamteffekt in fctotal schreiben

        double * fctotalbetap = fctotal.getbetapointer();

        if(gridsize < 0)
          {
          vector<int>::iterator freqwork = freq.begin();
          int * workindex = index.getV();
          for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
            {
            if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
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

      fctotal.posteriormode();

      } // END: interaction
    } // END: if(center)

// Interaktionseffekt in fchelp schreiben
    double * fchelpbetap = fchelp.getbetapointer();

    if(gridsize < 0)
      {
      if(varcoeff)
        multBout(splinehelp,beta);

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
      for(i=0;i<gridsize;i++,fchelpbetap++)
        *fchelpbetap = splinehelp(i,0);
      }

  fchelp.posteriormode();
  converged = FULLCOND_nonp_basis::posteriormode();

  if(converged && converged1 && converged2)
    return true;
  else
    return false;

  }


void FULLCOND_pspline_surf_gaussian::outresults(void)
  {
  FULLCOND::outresults();
  spline_basis_surf::outresults();
  }


void FULLCOND_pspline_surf_gaussian::outoptions(void)
  {

  ST::string typestr;
  ST::string knotstr;

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
  else if (type==mrfquadratic8)
    typestr = "2 dimensional second order random walk";
  else if (type==mrfkr1)
    typestr = "Kronecker product interaction (RW1*RW1)";
  else if (type==mrfkr2)
    typestr = "Kronecker product interaction (RW2*RW2)";
  else if (type == smoothspline)
    typestr = "Smoothing Splines";
  else if (type == npspline)
    typestr = "Natural P-Splines";

  if (knpos == equidistant)
    knotstr = "equidistant";
  else if (knpos == quantiles)
    knotstr = "quantiles";

  optionsp->out("  OPTIONS FOR P-SPLINE TERM: " + title + "\n",true);
  optionsp->out("\n");
  optionsp->out("  Prior: " + typestr + "\n");
  optionsp->out("  Number of knots: " + ST::inttostring(nrknots) + "\n" );
  optionsp->out("  Knot choice: " + knotstr + "\n");
  optionsp->out("  Degree of Splines: " + ST::inttostring(degree) + "\n" );
  optionsp->out("\n");

  if(utype != gaussian)
    {

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
    }

  if(utype == hyperblock || utype == hyperblockmode)
    {
    if(singleblock)
      optionsp->out("  Updating scheme: single block (including variance parameter)\n");
    else
      optionsp->out("  Updating scheme: rowwise (including variance parameter)\n");
    optionsp->out("  Starting value for tuning parameter f: " + ST::doubletostring(f) + "\n");
    }
  else
    {
    if(singleblock)
      optionsp->out("  Updating scheme: single block\n");
    else
      optionsp->out("  Updating scheme: rowwise\n");
    }

  optionsp->out("\n");

  }


void FULLCOND_pspline_surf_gaussian::add_linearpred_multBS(const datamatrix & b)
  {
/*
  unsigned i,j,k;

  int * workindex = index.getV();
  vector<int>::iterator freqwork = freq.begin();
  datamatrix *workl = &(likep->get_linearpred(true));
  double * workspline = spline.getV();

  for(i=0;i<spline.rows();i++,workspline++)
    *workspline = 0.0;

  for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
    {
    for(j=0;j<degree+1;j++)
      {
      for(k=0;k<degree+1;k++)
        {
        (*workl)(*workindex,0) += B(*freqwork,j+k*(degree+1))
                                      * b(first[i] + j + k*nrpar1dim,0);
        spline(*workindex,0) += B(*freqwork,j+k*(degree+1))
                                      * b(first[i] + j + k*nrpar1dim,0);
        }
      }
    }
*/

  unsigned i,j,k,l;
  double val=0.0;

  double *workbeta;
  double *workB;
  int *workindex;

  unsigned maxfirst = first[likep->get_nrobs()-1];
  unsigned nrobs = likep->get_nrobs();
  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator firstit = first.begin();

  i=0;
  j=0;
  workB = B.getV();
  workindex = index.getV();
  workbeta = b.getV();
  while(i<=maxfirst)
    {

    while(*firstit==i && j<nrobs)
      {
      if( freqwork==freq.begin() || *freqwork!=*(freqwork-1) )
        {
        val = 0.0;
        for(k=0;k<degree+1;k++)
          for(l=0;l<degree+1;l++,workB++)
            val += *workB * *(workbeta+l+k*nrpar1dim);
        }
      spline(*workindex,0) = val;
      workindex++;
      freqwork++;
      firstit++;
      j++;
      }
    i++;
    workbeta++;

    }

  likep->add_linearpred_m(spline,column,true);

  }


void FULLCOND_pspline_surf_gaussian::add_linearpred_multBS2(const datamatrix & b)
  {

  unsigned i,j,k,l;
  unsigned col = degree+1;
  double val=0.0;

  double *workbeta;
  double *workB;
  int *workindex;

  unsigned maxfirst = first[likep->get_nrobs()-1];
  unsigned nrobs = likep->get_nrobs();
  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator firstit = first.begin();

  likep->substr_linearpred_m(spline,column,true);

  i=0;
  j=0;
  workB = B.getV();
  workindex = index.getV();
  workbeta = b.getV();
  while(i<=maxfirst)
    {

    while(*firstit==i && j<nrobs)
      {
      if( freqwork==freq.begin() || *freqwork!=*(freqwork-1) )
        {
        val = 0.0;
        for(k=0;k<col;k++)
          for(l=0;l<col;l++,workB++)
            val += *workB * *(workbeta+l+k*nrpar1dim);
        }
      spline(*workindex,0) = val;
      workindex++;
      freqwork++;
      firstit++;
      j++;
      }
    i++;
    workbeta++;

    }

  likep->add_linearpred_m(spline,column,true);

  }


void FULLCOND_pspline_surf_gaussian::compute_XWX(const datamatrix & W,const unsigned & col)
  {
/*
  unsigned i,k,l,m,n;

  datamatrix XH1(nrpar,1,0);
  datamatrix XH2(nrpar,(nrpar1dim+1)*degree,0);

  vector<int>::iterator firstwork;
  vector<int>::iterator freqwork;

  int *workindex = index.getV();

  firstwork =  first.begin();
  freqwork = freq.begin();

  for(i=0;i<likep->get_nrobs();i++,++firstwork,++freqwork,++workindex)
    {
    for(k=0;k<degree+1;k++)
      {
      for(l=0;l<degree+1;l++)
        {
        for(m=0;m<degree+1;m++)
          {
          for(n=0;n<degree+1;n++)
            {
            if((*firstwork + k + l*nrpar1dim) == (*firstwork + m + n*nrpar1dim))
              XH1(*firstwork + k + l*nrpar1dim,0) +=
                B(*freqwork,k+l*(degree+1)) * B(*freqwork,m+n*(degree+1))
                    * W(*workindex,0);
            else if(k+l*nrpar1dim < m+n*nrpar1dim)
              XH2(*firstwork + k + l*nrpar1dim,m+n*nrpar1dim - (k+l*nrpar1dim)-1) +=
                B(*freqwork,k+l*(degree+1)) * B(*freqwork,m+n*(degree+1))
                    * W(*workindex,0);
            }
          }
        }
      }
    }

  XX.assign(XH1,XH2);
*/

  unsigned i,k,l,m,n;
  unsigned cols = degree+1;
  unsigned Bcols = B.cols();
  unsigned uppercols = (nrpar1dim+1)*degree;
  int dim = W.cols();

  vector<int>::iterator firstwork;
  vector<int>::iterator freqwork;
  vector<int>::iterator workindex2;

  int help;

  double * diag;
  double * upper;
  double * workweight;
  double * work1;
  double * work2;

  diag = XX.getdiagpointer();
  upper = XX.getupperpointer();

  for(m=0;m<nrpar;m++,diag++)
    {
    *diag = 0.0;
    for(n=0;n<uppercols;n++,upper++)
      *upper = 0.0;
    }

  firstwork = first.begin();
  freqwork = freq.begin();
  workindex2 = index2.begin();
  workweight = W.getV() + col + *workindex2*dim;

  for(i=0;i<likep->get_nrobs();i++,++firstwork,++freqwork,++workindex2,workweight+=(*workindex2*dim))
    {
    diag = XX.getdiagpointer() + *firstwork;
    upper = XX.getupperpointer() + *firstwork*uppercols;
    work1 = B.getV() + *freqwork*Bcols;
    for(k=0;k<cols;k++)
      {
      for(l=0;l<cols;l++,work1++)
        {
        work2 = B.getV() + *freqwork*Bcols;
        for(m=0;m<cols;m++)
          {
          for(n=0;n<cols;n++,work2++)
            {
            help = (n+m*nrpar1dim) - (l+k*nrpar1dim);
            if(help == 0)
              *(diag + l+k*nrpar1dim) += *work1 * *workweight * *work2;
            else if(help > 0)
              *(upper + (l+k*nrpar1dim)*uppercols + (help-1)) += *work1 * *workweight * *work2;
            }
          }
        }
      }
    }

  XX.set_decomposed();

  }


void FULLCOND_pspline_surf_gaussian::compute_XWXenv(const datamatrix & W,const unsigned & col)
  {

  compute_XWX(W,col);
  XX_env = envmatdouble(XX);

  }


void FULLCOND_pspline_surf_gaussian::compute_XWX_Block(const datamatrix & W,const unsigned a,const unsigned e,
                                                       const unsigned beg,const unsigned end,const unsigned & col)
  {

  unsigned i,k,l,m,n;
  datamatrix XH1(nrpar,1,0);
  datamatrix XH2(nrpar,(nrpar1dim+1)*degree,0);

  vector<int>::iterator firstwork = first.begin() + beg;
  vector<int>::iterator freqwork = freq.begin() + beg;
  int *workindex = index.getV() + beg;

  for(i=beg;i<=end;i++,++firstwork,++freqwork,++workindex)
    {
    for(k=0;k<degree+1;k++)
      {
      for(l=0;l<degree+1;l++)
        {
        for(m=0;m<degree+1;m++)
          {
          for(n=0;n<degree+1;n++)
            {

            if( a<=(*firstwork + k + l*nrpar1dim) && (*firstwork + k + l*nrpar1dim)<=e &&
                a<=(*firstwork + m + n*nrpar1dim) && (*firstwork + m + n*nrpar1dim)<=e )
              {

              if((*firstwork + k + l*nrpar1dim) == (*firstwork + m + n*nrpar1dim))
                XH1(*firstwork + k + l*nrpar1dim,0) +=
                  B(*freqwork,k+l*(degree+1)) * B(*freqwork,m+n*(degree+1))
                      * W(*workindex,col);
              else if(k+l*nrpar1dim < m+n*nrpar1dim)
                XH2(*firstwork + k + l*nrpar1dim,m+n*nrpar1dim - (k+l*nrpar1dim)-1) +=
                  B(*freqwork,k+l*(degree+1)) * B(*freqwork,m+n*(degree+1))
                      * W(*workindex,col);

              }

            }
          }
        }
      }
    }

  XX.assign(XH1,XH2);

  }


void FULLCOND_pspline_surf_gaussian::compute_XWXenv_Block(const datamatrix & W,const unsigned a,const unsigned e,
                                                       const unsigned beg,const unsigned end,const unsigned & col)
  {

  compute_XWX_Block(W,a,e,beg,end,col);
  XX_env = envmatdouble(XX);

  }


void FULLCOND_pspline_surf_gaussian::compute_XWtildey(const datamatrix & W,const double & scale)
  {
/*
  unsigned i,j,k;

  vector<int>::iterator freqwork;
  vector<int>::iterator firstwork;

  int * workindex = index.getV();

  double * muyp;

  muyp = muy.getV();
  for(i=0;i<nrpar;i++,muyp++)
    *muyp = 0.0;

  freqwork = freq.begin();
  firstwork = first.begin();

  for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
    {
    for(j=0;j<degree+1;j++)
      {
      for(k=0;k<degree+1;k++)
        {
        muy(first[i]+j+k*nrpar1dim,0) += B(*freqwork,j + k*(degree+1)) * W(*workindex,0)
                                         * mu(*workindex,0);
        }
      }
    }
*/

  unsigned i,j,k;
  unsigned cols = degree+1;

  vector<int>::iterator freqwork;
  vector<int>::iterator firstwork;
  vector<int>::iterator workindex2;

  double * muyp;
  double * Bwork;
  double * workmu;
  double * workweight;

  muyp = muy.getV();
  for(i=0;i<nrpar;i++,muyp++)
    *muyp = 0.0;

  freqwork = freq.begin();
  firstwork = first.begin();
  workindex2 = index2.begin();
  workweight = W.getV() + *workindex2;
  workmu = mu.getV() + *workindex2;

  for(i=0;i<likep->get_nrobs();i++,freqwork++,firstwork++)
    {
    Bwork = B.getV() + *freqwork*B.cols();
    for(j=0;j<cols;j++)
      {
      for(k=0;k<cols;k++,Bwork++)
        {
        muy(*firstwork+k+j*nrpar1dim,0) += *Bwork * *workweight * *workmu;
        }
      }
    workindex2++;
    workmu += *workindex2;
    workweight += *workindex2;
    }

  muyp = muy.getV();
  for(i=0;i<nrpar;i++,muyp++)
    *muyp *= scale;

  }


void FULLCOND_pspline_surf_gaussian::compute_XWtildey_Block(const datamatrix & W,const double & scale,
                                                            const unsigned & beg,const unsigned & end)
  {
/*
  unsigned i,j,k;
  unsigned cols = degree+1;

  vector<int>::iterator freqwork;
  int * workindex;
  double * muyp;

  muyp = muy.getV();
  for(i=0;i<nrpar;i++,muyp++)
    *muyp = 0.0;

  workindex = index.getV() + beg;
  freqwork = freq.begin() + beg;

  for(i=beg;i<=end;i++,freqwork++,workindex++)
    {
    for(j=0;j<cols;j++)
      {
      for(k=0;k<cols;k++)
        {
        muy(first[i]+j+k*nrpar1dim,0) += B(*freqwork,j + k*(degree+1)) * W(*workindex,0)
                                         * mu(*workindex,0);
        }
      }
    }
*/

  unsigned i,j,k;
  unsigned cols = degree+1;

  vector<int>::iterator freqwork;
  vector<int>::iterator firstwork;
  vector<int>::iterator workindex2;

  double * muyp;
  double * Bwork;
  double * workmu;
  double * workweight;

  muyp = muy.getV();
  for(i=0;i<nrpar;i++,muyp++)
    *muyp = 0.0;

  freqwork = freq.begin() + beg;
  firstwork = first.begin() + beg;
  workindex2 = index2.begin() + beg;
  workweight = W.getV() + index(beg,0);
  workmu = mu.getV() + index(beg,0);

  for(i=beg;i<=end;i++,freqwork++,firstwork++)
    {
    Bwork = B.getV() + *freqwork*B.cols();
    for(j=0;j<cols;j++)
      {
      for(k=0;k<cols;k++,Bwork++)
        {
        muy(*firstwork+k+j*nrpar1dim,0) += *Bwork * *workweight * *workmu;
        }
      }
    workindex2++;
    workmu += *workindex2;
    workweight += *workindex2;
    }

  muyp = muy.getV();
  for(i=0;i<nrpar;i++,muyp++)
    *muyp *= scale;


  }


void FULLCOND_pspline_surf_gaussian::compute_XWtildey(const datamatrix & W,const datamatrix & tildey,
                                                        const double & scale,const unsigned & col)
  {

  unsigned i,j,k;
  vector<int>::iterator freqwork = freq.begin();
  int * workindex = index.getV();
  double * muyp;

  muyp = muy.getV();
  for(i=0;i<nrpar;i++,muyp++)
    *muyp = 0.0;

  for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
    {
    for(j=0;j<degree+1;j++)
      {
      for(k=0;k<degree+1;k++)
        {
        muy(first[i]+j+k*nrpar1dim,0) += B(*freqwork,j + k*(degree+1)) * W(*workindex,col)
                                         * tildey(*workindex,0);
        }
      }
    }

  muyp = muy.getV();
  for(i=0;i<nrpar;i++,muyp++)
    *muyp *= scale;

  }


void FULLCOND_pspline_surf_gaussian::sample_centered(datamatrix & beta)
  {

  if(centertotal)
    {

    unsigned i;
    double help;
    double * v;
    double * work;

    datamatrix V = datamatrix(nrpar,1,0.0);
    prec.solve(betaweight,V,0,0);

    v = V.getV();
    work = betaweight.getV();

    help = 0.0;
    for(i=0;i<nrpar;++i,++work,++v)
      help += *work * *v;

    compute_intercept();
    help = intercept/help;

    v = V.getV();
    work = beta.getV();

    for(i=0;i<nrpar;++i,++work,++v)
      *work -= *v * help;

    intercept = 0.0;

    }
  else
    {

    unsigned i,j;
    datamatrix help = datamatrix(nrpar1dim,2,0.0);

// Intercept
    intercept = 0.0;
    for(i=0;i<nrpar;i++)
      intercept += beta(i,0)/nrpar;

// betax1mean, betax2mean
    for(i=0;i<nrpar1dim;i++)
      for(j=0;j<nrpar1dim;j++)
        help(i,0) += beta(i+j*nrpar1dim,0)/nrpar1dim;

    for(i=0;i<nrpar1dim;i++)
      for(j=0;j<nrpar1dim;j++)
        help(i,1) += beta(j+i*nrpar1dim,0)/nrpar1dim;

// betas ändern
    for(i=0;i<nrpar1dim;i++)
      for(j=0;j<nrpar1dim;j++)
        beta(i+j*nrpar1dim,0) -= help(i,0);

    for(i=0;i<nrpar1dim;i++)
      for(j=0;j<nrpar1dim;j++)
        beta(j+i*nrpar1dim,0) -= help(i,1);

    for(i=0;i<nrpar;i++)
      beta(i,0) += intercept;

    intercept = 0.0;

    }

  }


double FULLCOND_pspline_surf_gaussian::compute_df(void)
  {
  if(utype == gaussian)
    {
    if(prec.getdim()==0)
      return -1.0;
    if(lambda != lambda_prec || likep->iwlsweights_constant() == false)
      {
      if(XX.bandsize()<=K.bandsize())
        prec.addto2(XX,K,1.0,lambda);
      else
        prec.addto2(K,XX,lambda,1.0);
      }
    prec_env = envmatdouble(prec);
    XX_env = envmatdouble(XX);
    }
  else
    {
    if(prec_env.getDim()==0)
      return -1.0;
    if(lambda != lambda_prec || likep->iwlsweights_constant() == false)
      prec_env.addto(XX_env,Kenv,1.0,lambda);
    }

  invprec = envmatdouble(0,nrpar,prec_env.getBandwidth());
  prec_env.inverse_envelope(invprec);
  if(identifiable)
    return invprec.traceOfProduct(XX_env);
  else
    return invprec.traceOfProduct(XX_env)-1;
  }


} // end: namespace MCMC

//---------------------------------------------------------------------------
#if !defined(__BUILDING_GNU)
#pragma package(smart_init)
#endif

