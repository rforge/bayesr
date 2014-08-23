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



#include "mcmc_pspline_surf.h"


namespace MCMC
{


void FULLCOND_pspline_surf::init_maineffects(spline_basis * mp1,spline_basis * mp2,
                            const ST::string & pnt,const ST::string & prt)
  {

  mainp1 = mp1;
  mainp2 = mp2;

  assert(mainp1->get_nrknots() == nrknots);
  assert(mainp2->get_nrknots() == nrknots);

  centertotal = false;

  fctotalrespath = prt;

  datamatrix h(1,1,0);
  if(gridsize < 0)
    fctotal = FULLCOND(optionsp,h,title+"total",nrdiffobs,1,pnt);
  else
    fctotal = FULLCOND(optionsp,h,title+"total",gridsize,1,pnt);
  fctotal.setflags(MCMC::norelchange | MCMC::nooutput);

  beta1 = datamatrix(nrpar1dim,1,0);
  beta2 = datamatrix(nrpar1dim,1,0);

  }


FULLCOND_pspline_surf::FULLCOND_pspline_surf(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const unsigned & minb, const unsigned & maxb,
                         const double & l, const int & gs, const fieldtype & ft,
                         const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {

  assert(v1.rows() == v2.rows());

  varcoeff = false;
  mapexisting = false;

  unsigned i;

  outfile = of;

  centertotal = true;

  likep = dp;
  type = ft;
  pathresult = pres;

  min = minb;
  max = maxb;
  mintoobig = false;
  maxtoobig = false;
  oldnrtrials = 0;
  oldacceptance = 0;

  lambda = l;
  sigma2 = 1.0/l;

  intercept = 0.0;
  spline = datamatrix(likep->get_nrobs(),1,0);

  setbeta((nrk-1+degr)*(nrk-1+degr),1,0);

  nrknots = nrk;
  degree = degr;
  knpos = kp;

  nrpar1dim = nrknots-1+degree;

  make_index(v1,v2);
  compute_knots(v1,v2);
  make_B(v1,v2);

// fchelp initialisieren
  ST::string pnt = fp.substr(0,fp.length()-4)+"_fchelp.raw";
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
  make_xy_v(v1,v2);

  if(gridsize > 0)
    {
    make_xy_values(v1,v2);
    make_DG();
    }

  if (type == mrflinear)
    {
    Ksp = Kmrflinear(nrpar1dim,nrpar1dim);
    rankK = nrpar-1;

    if(minb == 0)
      {
      automatic = true;

      min = 1;
      if(minb == 0)
        max = 50;
      else
        max = maxb;

      minauto = int(sqrt(static_cast<double>(nrpar)/5.0));
      maxauto = int(sqrt(static_cast<double>(nrpar)/3.0));
      if(minauto < 1)
        minauto = 1;
      }
    else
      {
      automatic = false;

      if(max > rankK || max == 0)
        {
        maxtoobig = true;
        max = rankK;
        }
      if(min > max || min == 0)
        {
        mintoobig = true;
        min = 1;
        }
      }

    }
  else if(type == mrfkr1)          //  RW1 * RW1
    {
    Ksp = Krw1(vector<double>(nrpar1dim,1.0)).kronecker(Krw1(vector<double>(nrpar1dim,1.0)));
    rankK = (nrpar1dim-1)*(nrpar1dim-1);
    if(max > nrpar1dim-1)
      {
      maxtoobig = true;
      max = nrpar1dim-1;
      }
    }
  else if(type == mrfkr2)          //  RW2 * RW2
    {
    Ksp = Krw2(vector<double>(nrpar1dim,1.0)).kronecker(Krw2(vector<double>(nrpar1dim,1.0)));
    rankK = (nrpar1dim-2)*(nrpar1dim-2);
    if(max > nrpar1dim-2)
      {
      maxtoobig = true;
      max = nrpar1dim-2;
      }
    }

  weight = vector<double>(nrpar,1.0/double(nrpar));

  for(i=0;i<max;i++)
    {
    fc_random.push_back(datamatrix(i+1,1,0));
    randnorm.push_back(datamatrix(i+1,1,0));
    }

  make_Kab_list();

  identifiable = false;

  compute_betaweight();
  compute_betaweightxy();

  }


// geosplines
FULLCOND_pspline_surf::FULLCOND_pspline_surf(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc,
                         const datamatrix & region, const MAP::map & mp,
                         const ST::string & mn, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const unsigned & minb, const unsigned & maxb,
                         const double & l, const int & gs, const fieldtype & ft,
                         const ST::string & fp, const ST::string & pres, const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {

  unsigned i;

  varcoeff = false;

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

  centertotal = true;

  likep = dp;
  type = ft;
  pathresult = pres;

  min = minb;
  max = maxb;
  mintoobig = false;
  maxtoobig = false;
  oldnrtrials = 0;
  oldacceptance = 0;

  lambda = l;
  sigma2 = 1.0/l;

  intercept = 0.0;
  spline = datamatrix(likep->get_nrobs(),1,0);

  setbeta((nrk-1+degr)*(nrk-1+degr),1,0);

  nrknots = nrk;
  degree = degr;
  knpos = kp;

  nrpar1dim = nrknots-1+degree;

  make_index(v1,v2);
  compute_knots(v1,v2);
  make_B(v1,v2);

// fchelp initialisieren
  ST::string pnt = fp.substr(0,fp.length()-4)+"_fchelp.raw";
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
  make_xy_v(v1,v2);

  if(gridsize > 0)
    {
    make_xy_values(v1,v2);
    make_DG();
    }

  if (type == mrflinear)
    {
    Ksp = Kmrflinear(nrpar1dim,nrpar1dim);
    rankK = nrpar-1;

    if(minb == 0)
      {
      automatic = true;

      min = 1;
      if(minb == 0)
        max = 50;
      else
        max = maxb;

      minauto = int(sqrt(static_cast<double>(nrpar)/5.0));
      maxauto = int(sqrt(static_cast<double>(nrpar)/3.0));
      if(minauto < 1)
        minauto = 1;
      }
    else
      {
      automatic = false;

      if(max > rankK || max == 0)
        {
        maxtoobig = true;
        max = rankK;
        }
      if(min > max || min == 0)
        {
        mintoobig = true;
        min = 1;
        }
      }

    }
  else if(type == mrfkr1)          //  RW1 * RW1
    {
    Ksp = Krw1(vector<double>(nrpar1dim,1.0)).kronecker(Krw1(vector<double>(nrpar1dim,1.0)));
    rankK = (nrpar1dim-1)*(nrpar1dim-1);
    if(max > nrpar1dim-1)
      {
      maxtoobig = true;
      max = nrpar1dim-1;
      }
    }
  else if(type == mrfkr2)          //  RW2 * RW2
    {
    Ksp = Krw2(vector<double>(nrpar1dim,1.0)).kronecker(Krw2(vector<double>(nrpar1dim,1.0)));
    rankK = (nrpar1dim-2)*(nrpar1dim-2);
    if(max > nrpar1dim-2)
      {
      maxtoobig = true;
      max = nrpar1dim-1;
      }
    }

  if(min > max)
    {
    mintoobig = true;
    min = max;
    }

  weight = vector<double>(nrpar,1.0/double(nrpar));

  for(i=0;i<max;i++)
    {
    fc_random.push_back(datamatrix(i+1,1,0));
    randnorm.push_back(datamatrix(i+1,1,0));
    }

  make_Kab_list();

  identifiable = false;

  compute_betaweight();
  compute_betaweightxy();

  }


// varying coefficients
FULLCOND_pspline_surf::FULLCOND_pspline_surf(MCMCoptions * o,DISTRIBUTION * dp,
                         FULLCOND_const * fcc, const datamatrix & intact,
                         const datamatrix & v1, const datamatrix & v2, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const unsigned & minb, const unsigned & maxb,
                         const double & l, const int & gs, const fieldtype & ft,
                         const ST::string & fp, const ST::string & pres,
                         const ST::string & of, const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {

  assert(v1.rows() == v2.rows());

  varcoeff = true;
  mapexisting = false;

  unsigned i;

  outfile = of;

  centertotal = true;

  likep = dp;
  type = ft;
  pathresult = pres;

  min = minb;
  max = maxb;
  mintoobig = false;
  maxtoobig = false;
  oldnrtrials = 0;
  oldacceptance = 0;

  lambda = l;
  sigma2 = 1.0/l;

  intercept = 0.0;
  spline = datamatrix(likep->get_nrobs(),1,0);

  setbeta((nrk-1+degr)*(nrk-1+degr),1,0);

  nrknots = nrk;
  degree = degr;
  knpos = kp;

  nrpar1dim = nrknots-1+degree;

  make_index(v1,v2);
  compute_knots(v1,v2);
  make_B(v1,v2);
  if(varcoeff)
    make_BVC(intact);

// fchelp initialisieren
  ST::string pnt = fp.substr(0,fp.length()-4)+"_fchelp.raw";
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
  make_xy_v(v1,v2);

  if(gridsize > 0)
    {
    make_xy_values(v1,v2);
    make_DG();
    }

  if (type == mrflinear)
    {
    Ksp = Kmrflinear(nrpar1dim,nrpar1dim);
    rankK = nrpar-1;

    if(minb == 0)
      {
      automatic = true;

      min = 1;
      if(minb == 0)
        max = 50;
      else
        max = maxb;

      minauto = int(sqrt(static_cast<double>(nrpar)/5.0));
      maxauto = int(sqrt(static_cast<double>(nrpar)/3.0));
      if(minauto < 1)
        minauto = 1;
      }
    else
      {
      automatic = false;

      if(max > rankK || max == 0)
        {
        maxtoobig = true;
        max = rankK;
        }
      if(min > max || min == 0)
        {
        mintoobig = true;
        min = 1;
        }
      }

    }
  else if(type == mrfkr1)          //  RW1 * RW1
    {
    Ksp = Krw1(vector<double>(nrpar1dim,1.0)).kronecker(Krw1(vector<double>(nrpar1dim,1.0)));
    rankK = (nrpar1dim-1)*(nrpar1dim-1);
    if(max > nrpar1dim-1)
      {
      maxtoobig = true;
      max = nrpar1dim-1;
      }
    }
  else if(type == mrfkr2)          //  RW2 * RW2
    {
    Ksp = Krw2(vector<double>(nrpar1dim,1.0)).kronecker(Krw2(vector<double>(nrpar1dim,1.0)));
    rankK = (nrpar1dim-2)*(nrpar1dim-2);
    if(max > nrpar1dim-2)
      {
      maxtoobig = true;
      max = nrpar1dim-2;
      }
    }

  weight = vector<double>(nrpar,1.0/double(nrpar));

  for(i=0;i<max;i++)
    {
    fc_random.push_back(datamatrix(i+1,1,0));
    randnorm.push_back(datamatrix(i+1,1,0));
    }

  make_Kab_list();

  identifiable = true;

  compute_betaweight();
  compute_betaweightxy();

  }


// geosplines varying coefficients
FULLCOND_pspline_surf::FULLCOND_pspline_surf(MCMCoptions * o,DISTRIBUTION * dp,
                          FULLCOND_const * fcc,const datamatrix & intact,
                         const datamatrix & region, const MAP::map & mp, const ST::string & mn, const ST::string & ti,
                         const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                         const unsigned & minb, const unsigned & maxb,
                         const double & l, const int & gs, const fieldtype & ft,
                         const ST::string & fp, const ST::string & pres, const unsigned & c)
  : spline_basis_surf(o,dp,fcc,ft,ti,nrk,degr,kp,l,gs,gs,fp,pres,c)
  {

  unsigned i;

  varcoeff = true;

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

  centertotal = true;

  likep = dp;
  type = ft;
  pathresult = pres;

  min = minb;
  max = maxb;
  mintoobig = false;
  maxtoobig = false;
  oldnrtrials = 0;
  oldacceptance = 0;

  lambda = l;
  sigma2 = 1.0/l;

  intercept = 0.0;
  spline = datamatrix(likep->get_nrobs(),1,0);

  setbeta((nrk-1+degr)*(nrk-1+degr),1,0);

  nrknots = nrk;
  degree = degr;
  knpos = kp;

  nrpar1dim = nrknots-1+degree;

  make_index(v1,v2);
  compute_knots(v1,v2);
  make_B(v1,v2);
  if(varcoeff)
    make_BVC(intact);

// fchelp initialisieren
  ST::string pnt = fp.substr(0,fp.length()-4)+"_fchelp.raw";
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
  make_xy_v(v1,v2);

  if(gridsize > 0)
    {
    make_xy_values(v1,v2);
    make_DG();
    }

  if (type == mrflinear)
    {
    Ksp = Kmrflinear(nrpar1dim,nrpar1dim);
    rankK = nrpar-1;

    if(minb == 0)
      {
      automatic = true;

      min = 1;
      if(minb == 0)
        max = 50;
      else
        max = maxb;

      minauto = int(sqrt(static_cast<double>(nrpar)/5.0));
      maxauto = int(sqrt(static_cast<double>(nrpar)/3.0));
      if(minauto < 1)
        minauto = 1;
      }
    else
      {
      automatic = false;

      if(max > rankK || max == 0)
        {
        maxtoobig = true;
        max = rankK;
        }
      if(min > max || min == 0)
        {
        mintoobig = true;
        min = 1;
        }
      }

    }
  else if(type == mrfkr1)          //  RW1 * RW1
    {
    Ksp = Krw1(vector<double>(nrpar1dim,1.0)).kronecker(Krw1(vector<double>(nrpar1dim,1.0)));
    rankK = (nrpar1dim-1)*(nrpar1dim-1);
    if(max > nrpar1dim-1)
      {
      maxtoobig = true;
      max = nrpar1dim-1;
      }
    }
  else if(type == mrfkr2)          //  RW2 * RW2
    {
    Ksp = Krw2(vector<double>(nrpar1dim,1.0)).kronecker(Krw2(vector<double>(nrpar1dim,1.0)));
    rankK = (nrpar1dim-2)*(nrpar1dim-2);
    if(max > nrpar1dim-2)
      {
      maxtoobig = true;
      max = nrpar1dim-1;
      }
    }

  if(min > max)
    {
    mintoobig = true;
    min = max;
    }

  weight = vector<double>(nrpar,1.0/double(nrpar));

  for(i=0;i<max;i++)
    {
    fc_random.push_back(datamatrix(i+1,1,0));
    randnorm.push_back(datamatrix(i+1,1,0));
    }

  make_Kab_list();

  identifiable = true;

  compute_betaweight();
  compute_betaweightxy();

  }


FULLCOND_pspline_surf::FULLCOND_pspline_surf(const FULLCOND_pspline_surf & fc)
  : spline_basis_surf(spline_basis_surf(fc))
  {

  min = fc.min;
  max = fc.max;
  mintoobig = fc.mintoobig;
  maxtoobig = fc.maxtoobig;

  minauto = fc.minauto;
  maxauto = fc.maxauto;
  automatic = fc.automatic;
  oldacceptance = fc.oldacceptance;
  oldnrtrials = fc.oldnrtrials;

  fc_random = fc.fc_random;
  randnorm = fc.randnorm;

  KAB = fc.KAB;
  KABl_sp = fc.KABl_sp;
  KABr_sp = fc.KABr_sp;
  KABroot = fc.KABroot;
  begin = fc.begin;
  matquant = fc.matquant;

  }


const FULLCOND_pspline_surf & FULLCOND_pspline_surf::operator=(const FULLCOND_pspline_surf & fc)
  {
  if (this == &fc)
    return *this;
  spline_basis_surf::operator=(spline_basis_surf(fc));

  min = fc.min;
  max = fc.max;
  mintoobig = fc.mintoobig;
  maxtoobig = fc.maxtoobig;

  minauto = fc.minauto;
  maxauto = fc.maxauto;
  automatic = fc.automatic;
  oldacceptance = fc.oldacceptance;
  oldnrtrials = fc.oldnrtrials;

  fc_random = fc.fc_random;
  randnorm = fc.randnorm;

  KAB = fc.KAB;
  KABl_sp = fc.KABl_sp;
  KABr_sp = fc.KABr_sp;
  KABroot = fc.KABroot;
  begin = fc.begin;
  matquant = fc.matquant;

  return *this;
  }


void FULLCOND_pspline_surf::make_Kab_list(void)
  {

  unsigned i,j;
  unsigned nrupdate;
  unsigned a,b;
  datamatrix kab;
  datamatrix help;
  datamatrix null(1,1,0);

  unsigned sum=0;
  for(i=min;i<=max;i++)
    {
    nrupdate = Ksp.get_rows()/i;
    if ((nrupdate*i) < Ksp.get_rows())
      nrupdate++;
    sum+= nrupdate;
    }

  KAB.reserve(sum);
  KABroot.reserve(sum);
  KABr_sp.reserve(sum);
  KABl_sp.reserve(sum);

  for (i=min;i<=max;i++)
    {

    nrupdate = Ksp.get_rows()/i;
    if ((nrupdate*i) < Ksp.get_rows())
      nrupdate++;

    begin.push_back(KAB.size());
    matquant.push_back(nrupdate);

    for(j=1;j <= nrupdate;j++)
      {

      a = 1+(j-1)*i;

      if (j == nrupdate)
        b = Ksp.get_rows();
      else
        b = j*i;

      kab = (Ksp.getBlock(a-1,a-1,b,b)).inverse();
      KAB.push_back(kab);
      KABroot.push_back(kab.root());

      if (b != Ksp.get_cols())
        {
        help = kab*Ksp.getBlock(a-1,b,b,Ksp.get_cols());
        KABr_sp.push_back(SparseMatrix(help));
        }
      else
        KABr_sp.push_back(SparseMatrix());

      if (a != 1)
        {
        help = kab*Ksp.getBlock(a-1,0,b,a-1);
        KABl_sp.push_back(SparseMatrix(help));
        }
      else
        KABl_sp.push_back(SparseMatrix());

      } // end: for(j=1;j <= nrupdate;j++)

    } // end: for (i=min;i<=max;i++)

  } // end: function make_Kab_list


void FULLCOND_pspline_surf::compute_mu(const datamatrix & beta,
      const unsigned & bs,const unsigned & a,const unsigned & b,
      const unsigned & v)
  {
  unsigned matnr = begin[bs-min]+ ((a-1)/bs);

  if (a==1)
    KABr_sp[matnr].substr_mult(beta,b,v,fc_random[b-a],0);
  else if (b == nrpar)
    KABl_sp[matnr].substr_mult(beta,0,v,fc_random[b-a],0);
  else
    {
    KABr_sp[matnr].substr_mult(beta,b,v,fc_random[b-a],0);
    KABl_sp[matnr].substr_mult(beta,0,v,fc_random[b-a],0);
    }
  }


void FULLCOND_pspline_surf::compute_fc(const datamatrix & beta, const unsigned & bs,
                               const unsigned & a,const unsigned & b,
                               const double & Q,const unsigned & v)
  {
  unsigned matnr = begin[bs-min]+ ((a-1)/bs );
  unsigned l = b-a+1;
  register unsigned i,j;

  double * fc_randwork = fc_random[b-a].getV();
  double * KABrootwork = KABroot[matnr].getV();

  double * randnormwork = randnorm[b-a].getV();

  for(i=0;i<l;i++,randnormwork++)
    *randnormwork = rand_normal();

  for (i=0;i<l;i++,fc_randwork++)
    {
    *fc_randwork = 0;
    randnormwork = randnorm[b-a].getV();
    for(j=0;j<l;j++,KABrootwork++,randnormwork++)
      *fc_randwork += *KABrootwork * *randnormwork;

    *fc_randwork *= Q;
    }

  compute_mu(beta,bs,a,b,v);

  }


void FULLCOND_pspline_surf::outoptions(void)
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

  if(mintoobig || maxtoobig)
    optionsp->out("\n");

  if(maxtoobig)
    optionsp->out("NOTE:  Maximum blocksize is too big, "
                    + ST::inttostring(max) + " has been used instead\n");
  if(mintoobig)
    optionsp->out("NOTE:  Minimum blocksize is too big, "
                      + ST::inttostring(min) + " has been used instead\n");

  optionsp->out("\n");
  optionsp->out("  Prior: " + typestr + "\n");
  optionsp->out("  Number of knots: " + ST::inttostring(nrknots) + "\n" );
  optionsp->out("  Knot choice: " + knotstr + "\n");
  optionsp->out("  Degree of Splines: " + ST::inttostring(degree) + "\n" );
  optionsp->out("\n");

  if(automatic)
    {
    optionsp->out("  Initial minimum blocksize for automatic tuning: " + ST::inttostring(minauto) + "\n");
    optionsp->out("  Initial maximum blocksize for automatic tuning: " + ST::inttostring(maxauto) + "\n");
    }
  else
    {
    optionsp->out("  Minimum blocksize: " + ST::inttostring(min) + "\n");
    optionsp->out("  Maximum blocksize: " + ST::inttostring(max) + "\n");
    }

  optionsp->out("\n");

  }


void FULLCOND_pspline_surf::update(void)
  {

  if(lambdaconst == true)
    sigma2 = likep->get_scale(column)/lambda;

  unsigned blocksize;

  if(automatic)
    {
    if(optionsp->get_nriter()%100==0 && optionsp->get_nriter()<optionsp->get_burnin())
      adjust_blocksize(30,70);
    if(optionsp->get_nriter() == optionsp->get_burnin())
      {
      optionsp->out("\n");
      optionsp->out("NOTE: Minimum blocksize for " + title + " set to " + ST::inttostring(minauto) + "\n");
      optionsp->out("NOTE: Maximum blocksize for " + title + " set to " + ST::inttostring(maxauto) + "\n");
      optionsp->out("\n");
      }
#ifndef __BUILDING_GNU
    blocksize = minauto + random(maxauto-minauto+1);
#else
    blocksize = minauto + int((maxauto-minauto+1)*rand()/(RAND_MAX + 1.0));
#endif
    }
  else
    {
#ifndef __BUILDING_GNU
    blocksize = min + random(max-min+1);
#else
    blocksize = min + int((max-min+1)*rand()/(RAND_MAX + 1.0));
#endif
    }

  double u;
  unsigned an = 1;
  unsigned en = blocksize;

  unsigned beg;
  unsigned end;
  vector<int>::iterator firstit;

  double logold;
  double logprop;
  double * workbeta;
  unsigned i,j,k;

  for(j=0;j<matquant[blocksize-min];j++)
    {

    beg = 0;
    firstit = first.begin();
    while(firstit != first.end() && *firstit+nrpar1dim*degree+(degree+1) < an)
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

    compute_fc(beta,blocksize,an,en,sqrt(sigma2));

    logold = 0;
    logprop = 0;

    likep->assign(false);

    logold += likep->loglikelihood(beg,end,index);

    add_linearpred_multBS_Block(an-1,en-1,beg,end);

    logprop += likep->loglikelihood(beg,end,index,false);

    u = log(uniform());

    if (u <= (logprop-logold))
      {
      workbeta = beta.getV()+an-1;          // change beta
      for(k=an-1;k<en;k++,workbeta++)
        *workbeta = fc_random[en-an](k-an+1,0);
	  acceptance++;
      likep->swap_linearpred();
      }

    an+=blocksize;
    if (j == matquant[blocksize-min]-2)
      en = nrpar;
    else
      en+=blocksize;

    } // end: for(j=0;j<matquant[blocksize-min];j++)

  if(center)
    {
    if(centertotal)
      {
      compute_intercept();
      for(i=0;i<nrpar;i++)
        beta(i,0) -= intercept;
//      likep->add_linearpred_m(-intercept,column);
      fcconst->update_intercept(intercept);
      }
    else
      {

      beta_uncentered.assign(beta);

      compute_intercept();
      compute_beta();
//      likep->add_linearpred_m(-intercept,column);
      fcconst->update_intercept(intercept);
      mainp1->change(beta1);
      mainp2->change(beta2);

// Gesamteffekt in fctotal schreiben
      if( (optionsp->get_nriter() > optionsp->get_burnin()) &&
          ((optionsp->get_nriter()-optionsp->get_burnin()-1) % (optionsp->get_step()) == 0) )
        {

        double * fctotalbetap = fctotal.getbetapointer();

        if(gridsize < 0)
          {
          multBS(splinehelp,beta);
          vector<int>::iterator freqwork = freq.begin();
          int * workindex = index.getV();
          for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
            {
            if(freqwork==freq.begin() || *freqwork!=*(freqwork-1))
              {
              *fctotalbetap = splinehelp(i,0)
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
      else
        multBS(splinehelp,beta);

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
      for(i=0;i<gridsize;i++,fchelpbetap++)
        *fchelpbetap = splinehelp(i,0);
      }
    }

  fchelp.update();

  FULLCOND::update();

  } // end: function update


void FULLCOND_pspline_surf::outresults(void)
  {
  FULLCOND::outresults();
  spline_basis_surf::outresults();
  }


void FULLCOND_pspline_surf::add_linearpred_multBS_Block(const unsigned a,const unsigned e,
                                                        const unsigned beg,const unsigned end)
  {

  datamatrix *workl = &(likep->get_linearpred(false));
  int *workindex;

  unsigned i,j,k;

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator firstit = first.begin();

  workindex = index.getV() + beg;
  freqwork = freq.begin() + beg;
  firstit = first.begin() + beg;

  for(i=beg;i<=end;i++,freqwork++,workindex++,firstit++)
    {
    for(j=0;j<degree+1;j++)
      {
      for(k=0;k<degree+1;k++)
        {
        if(a <= *firstit+j+k*nrpar1dim && *firstit+j+k*nrpar1dim <= e && *firstit+j+k*nrpar1dim < nrpar)
          {
          (*workl)(*workindex,column) += B(*freqwork,j+k*(degree+1))
                                 * (fc_random[e-a]((*firstit+j+k*nrpar1dim)-a,0) - beta(*firstit + j + k*nrpar1dim,0));
          }
        }
      }
    }

  }


void FULLCOND_pspline_surf::adjust_blocksize(const unsigned & alphamin,const unsigned & alphamax)
  {

  int min = minauto;
  int max = maxauto;

  double rate;
  if (nrtrials == 0)
    rate = (double(acceptance-oldacceptance)/double(100))*100;
  else
    rate = (double(acceptance-oldacceptance)/double(nrtrials-oldnrtrials))*100;

  oldacceptance = acceptance;
  oldnrtrials = nrtrials;

  int limit = 1;
  int span = sqrt((static_cast<double>(nrpar)/10.0)>1?(static_cast<double>(nrpar)/10.0):2.0);

  if(rate<alphamin)
    {
    if(max-min<span)
      {
      if(rate<(alphamin-15))
        min -= span;
      else
        min -= 2;
      if(min<limit)
        min = limit;
      }
    else
      {
      if(rate<(alphamin-15))
        max -= span;
      else
        max -= 2;
      if(max<min)
        max = min;
      }
    }

  if(rate>alphamax)
    {
    if(max-min<span)
      {
      if(rate>(alphamax+15))
        max += span;
      else
        max += 2;
      if(max>FULLCOND_pspline_surf::max)
        max = FULLCOND_pspline_surf::max;
      }
    else
      {
      if(rate>(alphamax+15))
        min += span;
      else
        min += 2;
      if(min>max)
        min = max;
      }
    }

  minauto = min;
  maxauto = max;

  }


} // end: namespace MCMC

















