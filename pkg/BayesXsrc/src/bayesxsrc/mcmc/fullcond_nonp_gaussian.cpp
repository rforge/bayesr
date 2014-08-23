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



#include "fullcond_nonp_gaussian.h"

namespace MCMC
{

//---------------------------FOR REML ------------------------------------------

// REML: spatial covariates

FULLCOND_nonp_gaussian::FULLCOND_nonp_gaussian(MCMCoptions * o,
                        const datamatrix & d,const MAP::map & m,
                        const ST::string & mn,const ST::string & ti,
                        const ST::string & fp, const ST::string & pres,
                        const double & l, const double & sl, const bool & catsp)
  : FULLCOND_nonp_basis(o,ti)

  {

  notransform=false;
  catspecific = catsp;

  pathcurrent = pres;
  pathresult = pres;

  fctype = spatial;
  type = mrf;
  varcoeff=false;

  utype = gaussian;

  mapname = mn;

  lambda = l;
  startlambda = sl;

  unsigned i;

  if (m.polygones_existing())
    polex = true;
  else
    polex = false;

  if(polex == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

#ifdef __BUILDING_GNU
  MAP::map type15;
  type15 = m;
  type15.compute_reg(d,posbeg,posend,effectvalues,index);
#else
  m.compute_reg(d,posbeg,posend,effectvalues,index);
#endif

  if (m.get_errormessages().size() > 0)
    errors = m.get_errormessages();
  else
    {

    datamatrix Kstat=STATMAT_PENALTY::Kmrf(m);
    datamatrix vals(Kstat.rows(),1,0);

    nrpar=Kstat.rows();

    bool eigentest=eigen2(Kstat,vals);
    if(eigentest==false)
      {
      errors.push_back("ERROR: Unable to compute eigen decomposition for structured spatial effect.\n");
      }
    else
      {
      eigensort(vals,Kstat);

      for(i=0; i<vals.rows()-1; i++)
        {
        vals(i,0)=1/sqrt(vals(i,0));
        }
      vals(vals.rows()-1,0)=0;
      remlspatialdesign = multdiagback(Kstat,vals).getColBlock(0,Kstat.cols()-1);

      for (i=0;i<posbeg.size();i++)
        {
        if (posbeg[i] == -1)
          optionsp->out("NOTE: no observations for region "
          + effectvalues[i] + "\n");
        }
      }

    } // end: if (error==false)

  dimX=0;
  dimZ=posbeg.size()-1;

  }

// REML: spatial covariates (VCM)

FULLCOND_nonp_gaussian::FULLCOND_nonp_gaussian(MCMCoptions * o,
                        const datamatrix & d1, const datamatrix & d2,
                        const MAP::map & m, const ST::string & mn,
                        const ST::string & ti, const ST::string & fp,
                        const ST::string & pres, const double & l,
                        const double & sl, const bool & catsp,
                        const bool & ctr)
  : FULLCOND_nonp_basis(o,ti)

  {
  notransform=false;
  catspecific = catsp;
  centervcm=ctr;

  data_forfixed = d2;

  pathcurrent = pres;
  pathresult = pres;

  fctype = spatial;
  type = mrf;
  varcoeff=true;

  utype = gaussian;

  mapname = mn;

  lambda = l;
  startlambda = sl;

  unsigned i;

  if (m.polygones_existing())
    polex = true;
  else
    polex = false;

  if(polex == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

#ifdef __BUILDING_GNU
  MAP::map temp16;
  temp16 = m;
  temp16.compute_reg(d1,posbeg,posend,effectvalues,index);
#else
  m.compute_reg(d1,posbeg,posend,effectvalues,index);
#endif

  if (m.get_errormessages().size() > 0)
    errors = m.get_errormessages();
  else
    {

    datamatrix Kstat=STATMAT_PENALTY::Kmrf(m);
    datamatrix vals(Kstat.rows(),1,0);

    nrpar=Kstat.rows();

    bool eigentest=eigen2(Kstat,vals);
    if(eigentest==false)
      {
      errors.push_back("ERROR: Unable to compute eigen decomposition for structured spatial effect.\n");
      }
    else
      {
      eigensort(vals,Kstat);

      for(i=0; i<vals.rows()-1; i++)
        {
        vals(i,0)=1/sqrt(vals(i,0));
        }
      vals(vals.rows()-1,0)=0;
      remlspatialdesign = multdiagback(Kstat,vals).getColBlock(0,Kstat.cols()-1);

      for (i=0;i<posbeg.size();i++)
        {
        if (posbeg[i] == -1)
          optionsp->out("NOTE: no observations for region "
          + effectvalues[i] + "\n");
        }
      }

    } // end: if (error==false)

  dimX=1;
  dimZ=posbeg.size()-1;

  if(centervcm)
    {
    dimX = dimX-1;
    }

  X_VCM = datamatrix(d1.rows(),dimX,1.0);
  }

// REML: additive, RW1, RW2, seasonal

FULLCOND_nonp_gaussian::FULLCOND_nonp_gaussian(MCMCoptions * o,
                        const datamatrix & d,const unsigned & maxint,
                        const fieldtype & ft,const ST::string & ti,
                        const ST::string & pres,const double & l,
                        const double & sl, const bool & catsp,
                        const unsigned & per)
  : FULLCOND_nonp_basis(o,ti)
  {
  notransform=false;
  catspecific = catsp;

  fctype = nonparametric;

  type = ft;

  isnonparametric = true;
  varcoeff=false;

  title = ti;

  period=per;

  plotstyle = plotnonp;
  pathresult = pres;
  pathcurrent = pres;

  lambda = l;
  startlambda = sl;

  make_categories(d,maxint);

  if (errors.size() == 0)
    {

    setbeta(posbeg.size(),1,0);

    if (type==MCMC::RW1)
      {
      dimX = 0;
      dimZ = nrpar-1;
      }
    else if (type==MCMC::RW2)
      {
      dimX = 1;
      dimZ = nrpar-2;
      }
    else if (type==MCMC::seasonal)
      {
      dimX = per-1;
      dimZ = nrpar-per+1;
      }

    } // end: errors.size() == 0


  }

// REML: additive, RW1, RW2, seasonal (VCM)

FULLCOND_nonp_gaussian::FULLCOND_nonp_gaussian(MCMCoptions * o,
                        const datamatrix & d1, const datamatrix & d2,
                        const unsigned & maxint,
                        const fieldtype & ft,const ST::string & ti,
                        const ST::string & pres,const double & l,
                        const double & sl,
                        const bool & catsp, const bool & ctr,
                        const unsigned & per)
  : FULLCOND_nonp_basis(o,ti)
  {
  notransform=false;
  catspecific = catsp;
  centervcm = ctr;

  fctype = nonparametric;

  data_forfixed = d2;

  type = ft;

  isnonparametric = true;
  varcoeff=true;

  title = ti;

  period=per;

  plotstyle = plotnonp;
  pathresult = pres;
  pathcurrent = pres;

  lambda = l;
  startlambda = sl;

  make_categories(d1,maxint);

  if (errors.size() == 0)
    {

    setbeta(posbeg.size(),1,0);

    if (type==MCMC::RW1)
      {
      dimX = 1;
      dimZ = nrpar-1;
      }
    else if (type==MCMC::RW2)
      {
      dimX = 2;
      dimZ = nrpar-2;
      }
    else if (type==MCMC::seasonal)
      {
      dimX = per-1;
      dimZ = nrpar-per+1;
      }

    if(centervcm)
      {
      dimX = dimX-1;
      }

    X_VCM = datamatrix(d1.rows(),dimX,1.0);
    Z_VCM = datamatrix(d1.rows(),dimZ,0.0);

    } // end: errors.size() == 0
  }

void FULLCOND_nonp_gaussian::createreml(datamatrix & X,datamatrix & Z,
                                  const unsigned & Xpos, const unsigned & Zpos)
  {
  unsigned i,k;
  int j;
  if(!varcoeff)
    {
    if(fctype==MCMC::spatial)
      {
      for(i=0;i<posbeg.size();i++)
        {
        for (j=posbeg[i];j<=posend[i];j++)
          {
          if(j!=-1)
            {
            for(k=0; k<remlspatialdesign.cols(); k++)
              {
              Z(index(j,0),Zpos+k)=remlspatialdesign(i,k);
              }
            }
          }
        }
      }
    else if(type==MCMC::seasonal)
      {
      datamatrix factor = seasonalfactor(period,weight.size());
      factor = factor*factor.sscp().inverse();
      datamatrix smallx = seasonalX(period,weight.size());
      for(i=0;i<posbeg.size();i++)
        {
        for (j=posbeg[i];j<=posend[i];j++)
          {
          for(k=0; k<dimZ; k++)
            {
            Z(index(j,0),Zpos+k)=factor(i,k);
            }
          for(k=0; k<dimX; k++)
            {
            X(index(j,0),Xpos+k)=smallx(i,k);
            }
          }
        }
      }
    else if(type==MCMC::RW1)
      {
      datamatrix diffmatrix = weighteddiffmat(1,weight);
      diffmatrix = diffmatrix.transposed()*diffmatrix.transposed().sscp().inverse();
      for(i=0;i<posbeg.size();i++)
        {
        for (j=posbeg[i];j<=posend[i];j++)
          {
          for(k=0; k<nrpar-1; k++)
            {
            Z(index(j,0),Zpos+k)=diffmatrix(i,k);
            }
          }
        }
      }
    else if(type==MCMC::RW2)
      {
      datamatrix diffmatrix = weighteddiffmat(2,weight);
      diffmatrix = diffmatrix.transposed()*diffmatrix.transposed().sscp().inverse();

      datamatrix xhelp(weight.size(),1,0);
      xhelp(0,0)=0;
      for(i=1; i<xhelp.rows(); i++)
        {
        xhelp(i,0)=xhelp(i-1,0)+weight[i];
        }

      for(i=0;i<posbeg.size();i++)
        {
        for (j=posbeg[i];j<=posend[i];j++)
          {
          X(index(j,0),Xpos)=xhelp(i,0);
          for(k=0; k<nrpar-2; k++)
            {
            Z(index(j,0),Zpos+k)=diffmatrix(i,k);
            }
          }
        }
      double xmean=X.mean(Xpos);
      for(i=0; i<X.rows(); i++)
        {
       X(i,Xpos)=X(i,Xpos)-xmean;
        }
      }
    }
  else
    {
    if(fctype==MCMC::spatial)
      {
      if(!centervcm)
        {
        for(i=0; i<X.rows(); i++)
          {
          X(i,Xpos)=data_forfixed(i,0);
          }
        }
      for(i=0;i<posbeg.size();i++)
        {
        for (j=posbeg[i];j<=posend[i];j++)
          {
          if(j!=-1)
            {
            for(k=0; k<remlspatialdesign.cols(); k++)
              {
              Z(index(j,0),Zpos+k)=remlspatialdesign(i,k)*data_forfixed(index(j,0),0);
              }
            }
          }
        }
      }
    else if(type==MCMC::seasonal)
      {
      datamatrix factor = seasonalfactor(period,weight.size());
      factor = factor*factor.sscp().inverse();
      datamatrix smallx = seasonalX(period,weight.size());
      for(i=0;i<posbeg.size();i++)
        {
        for (j=posbeg[i];j<=posend[i];j++)
          {
          for(k=0; k<dimZ; k++)
            {
            Z(index(j,0),Zpos+k)=factor(i,k)*data_forfixed(index(j,0),0);
            Z_VCM(index(j,0),k)=factor(i,k);
            }
          for(k=0; k<dimX; k++)
            {
            X(index(j,0),Xpos+k)=smallx(i,k)*data_forfixed(index(j,0),0);
            X_VCM(index(j,0),k)=smallx(i,k);
            }
          }
        }
      }
    else if(type==MCMC::RW1)
      {
      if(!centervcm)
        {
        for(i=0; i<X.rows(); i++)
          {
          X(i,Xpos)=data_forfixed(i,0);
          }
        }
      datamatrix diffmatrix = weighteddiffmat(1,weight);
      diffmatrix = diffmatrix.transposed()*diffmatrix.transposed().sscp().inverse();
      for(i=0;i<posbeg.size();i++)
        {
        for (j=posbeg[i];j<=posend[i];j++)
          {
          for(k=0; k<nrpar-1; k++)
            {
            Z(index(j,0),Zpos+k)=diffmatrix(i,k)*data_forfixed(index(j,0),0);
            Z_VCM(index(j,0),k)=diffmatrix(i,k);
            }
          }
        }
      }
    else if(type==MCMC::RW2)
      {
      datamatrix diffmatrix = weighteddiffmat(2,weight);
      diffmatrix = diffmatrix.transposed()*diffmatrix.transposed().sscp().inverse();

      datamatrix xhelp(weight.size(),1,0);
      xhelp(0,0)=0;
      for(i=1; i<xhelp.rows(); i++)
        {
        xhelp(i,0)=xhelp(i-1,0)+weight[i];
        }

      for(i=0;i<posbeg.size();i++)
        {
        for (j=posbeg[i];j<=posend[i];j++)
          {
          if(!centervcm)
            {
            X(index(j,0),Xpos)=data_forfixed(index(j,0),0);
            X_VCM(index(j,0),1)=xhelp(i,0);
            X(index(j,0),Xpos+1)=xhelp(i,0)*data_forfixed(index(j,0),0);
            }
          else
            {
            X_VCM(index(j,0),0)=xhelp(i,0);
            X(index(j,0),Xpos)=xhelp(i,0)*data_forfixed(index(j,0),0);
            }
          for(k=0; k<nrpar-2; k++)
            {
            Z(index(j,0),Zpos+k)=diffmatrix(i,k)*data_forfixed(index(j,0),0);
            Z_VCM(index(j,0),k)=diffmatrix(i,k);
            }
          }
        }
      }
    }
  }

double FULLCOND_nonp_gaussian::outresultsreml(datamatrix & X,datamatrix & Z,
                                     datamatrix & betareml,datamatrix & betacov,
                                     datamatrix & thetareml,
                                     const unsigned & Xpos,
                                     const unsigned & Zpos,
                                     const unsigned & thetapos,
                                     const bool & dispers,
                                     const unsigned & betaXpos,
                                     const unsigned & betaZpos,
                                     const double & category,
                                     const bool & ismultinomial,
                                     const unsigned plotpos)
  {
  unsigned i,j;
  double mean=0;

  betamean=datamatrix(nrpar,1,0);
  betavar=datamatrix(nrpar,1,0);
  datamatrix betastd=datamatrix(nrpar,1,0);
  betaqu_l1_lower=datamatrix(nrpar,1,0);
  betaqu_l1_upper=datamatrix(nrpar,1,0);
  betaqu_l2_lower=datamatrix(nrpar,1,0);
  betaqu_l2_upper=datamatrix(nrpar,1,0);

  if(!varcoeff)
    {
    if(fctype==MCMC::spatial)
      {
//      betamean=remlspatialdesign*betareml.getBlock(X.cols()+Zpos,0,X.cols()+Zpos+nrpar-1,1);
      betamean=remlspatialdesign*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1);
      for(i=0; i<nrpar; i++)
        {
//        betastd(i,0)=sqrt((remlspatialdesign.getRow(i)*
//                     betacov.getBlock(X.cols()+Zpos,X.cols()+Zpos,X.cols()+Zpos+nrpar-1,X.cols()+Zpos+nrpar-1)*
//                     remlspatialdesign.getRow(i).transposed())(0,0));
        betastd(i,0)=sqrt((remlspatialdesign.getRow(i)*
                     betacov.getBlock(betaZpos,betaZpos,betaZpos+nrpar-1,betaZpos+nrpar-1)*
                     remlspatialdesign.getRow(i).transposed())(0,0));
        }
      }
    else if(type==MCMC::RW1)
      {
      for(i=0;i<nrpar;i++)
        {
        j=index(posbeg[i],0);
//        betamean(i,0)=(Z.getBlock(j,Zpos,j+1,Zpos+nrpar-1)*betareml.getBlock(X.cols()+Zpos,0,X.cols()+Zpos+nrpar-1,1))(0,0);
//        betastd(i,0)=sqrt((Z.getBlock(j,Zpos,j+1,Zpos+nrpar-1)*
//                     betacov.getBlock(X.cols()+Zpos,X.cols()+Zpos,X.cols()+Zpos+nrpar-1,X.cols()+Zpos+nrpar-1)*
//                     Z.getBlock(j,Zpos,j+1,Zpos+nrpar-1).transposed())(0,0));
        betamean(i,0)=(Z.getBlock(j,Zpos,j+1,Zpos+nrpar-1)*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1))(0,0);
        betastd(i,0)=sqrt((Z.getBlock(j,Zpos,j+1,Zpos+nrpar-1)*
                     betacov.getBlock(betaZpos,betaZpos,betaZpos+nrpar-1,betaZpos+nrpar-1)*
                     Z.getBlock(j,Zpos,j+1,Zpos+nrpar-1).transposed())(0,0));
        }
      }
    else if(type==MCMC::RW2)
      {
      for(i=0;i<nrpar;i++)
        {
        j=index(posbeg[i],0);
/*        betamean(i,0) = betareml(Xpos,0)*X(j,Xpos) + (Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betareml.getBlock(X.cols()+Zpos,0,X.cols()+Zpos+dimZ,1))(0,0);
        betastd(i,0) = sqrt(
                           (
                            X(j,Xpos)*betacov(Xpos,Xpos)
                            +
                            (Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betacov.getBlock(X.cols()+Zpos,Xpos,X.cols()+Zpos+dimZ,Xpos+1))(0,0)
                           )*X(j,Xpos)
                           +
                           (
                            (
                             X(j,Xpos)*betacov.getBlock(Xpos,X.cols()+Zpos,Xpos+1,X.cols()+Zpos+dimZ)
                             +
                             Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betacov.getBlock(X.cols()+Zpos,X.cols()+Zpos,X.cols()+Zpos+dimZ,X.cols()+Zpos+dimZ)
                           )*(Z.getBlock(j,Zpos,j+1,Zpos+dimZ).transposed())
                           )(0,0)
                          );*/
        betamean(i,0) = betareml(betaXpos,0)*X(j,Xpos) + (Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1))(0,0);
        betastd(i,0) = sqrt(
                           (
                            X(j,Xpos)*betacov(betaXpos,betaXpos)
                            +
                            (Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+1))(0,0)
                           )*X(j,Xpos)
                           +
                           (
                            (
                             X(j,Xpos)*betacov.getBlock(betaXpos,betaZpos,betaXpos+1,betaZpos+dimZ)
                             +
                             Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                           )*(Z.getBlock(j,Zpos,j+1,Zpos+dimZ).transposed())
                           )(0,0)
                          );
        }
      }
    else if(type==MCMC::seasonal)
      {
      for(i=0;i<nrpar;i++)
        {
        j=index(posbeg[i],0);
/*        betamean(i,0) = (X.getBlock(j,Xpos,j+1,Xpos+dimX)*betareml.getBlock(Xpos,0,Xpos+dimX,1))(0,0) + (Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betareml.getBlock(X.cols()+Zpos,0,X.cols()+Zpos+dimZ,1))(0,0);
        betastd(i,0) = sqrt(
                           ((
                            X.getBlock(j,Xpos,j+1,Xpos+dimX)*betacov.getBlock(Xpos,Xpos,Xpos+dimX,Xpos+dimX)
                            +
                            (Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betacov.getBlock(X.cols()+Zpos,Xpos,X.cols()+Zpos+dimZ,Xpos+dimX))
                           )*X.getBlock(j,Xpos,j+1,Xpos+dimX).transposed())(0,0)
                           +
                           (
                             (
                             X.getBlock(j,Xpos,j+1,Xpos+dimX)*betacov.getBlock(Xpos,X.cols()+Zpos,Xpos+dimX,X.cols()+Zpos+dimZ)
                             +
                             Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betacov.getBlock(X.cols()+Zpos,X.cols()+Zpos,X.cols()+Zpos+dimZ,X.cols()+Zpos+dimZ)
                           )*(Z.getBlock(j,Zpos,j+1,Zpos+dimZ).transposed())
                           )(0,0)
                          );*/
        betamean(i,0) = (X.getBlock(j,Xpos,j+1,Xpos+dimX)*betareml.getBlock(betaXpos,0,betaXpos+dimX,1))(0,0) + (Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1))(0,0);
        betastd(i,0) = sqrt(
                           ((
                            X.getBlock(j,Xpos,j+1,Xpos+dimX)*betacov.getBlock(betaXpos,betaXpos,betaXpos+dimX,betaXpos+dimX)
                            +
                            (Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+dimX))
                           )*X.getBlock(j,Xpos,j+1,Xpos+dimX).transposed())(0,0)
                           +
                           (
                             (
                             X.getBlock(j,Xpos,j+1,Xpos+dimX)*betacov.getBlock(betaXpos,betaZpos,betaXpos+dimX,betaZpos+dimZ)
                             +
                             Z.getBlock(j,Zpos,j+1,Zpos+dimZ)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                           )*(Z.getBlock(j,Zpos,j+1,Zpos+dimZ).transposed())
                           )(0,0)
                          );
        }
      }

    mean = betamean.mean(0);
//    betareml(0,0)=betareml(0,0)+mean;
    for(j=0; j<nrpar; j++)
      {
      betamean(j,0)=betamean(j,0)-mean;
      betaqu_l1_lower(j,0) = betamean(j,0)+randnumbers::invPhi2(lower1/100)*betastd(j,0);
      betaqu_l1_upper(j,0) = betamean(j,0)+randnumbers::invPhi2(upper2/100)*betastd(j,0);
      betaqu_l2_lower(j,0) = betamean(j,0)+randnumbers::invPhi2(lower2/100)*betastd(j,0);
      betaqu_l2_upper(j,0) = betamean(j,0)+randnumbers::invPhi2(upper1/100)*betastd(j,0);
      }
    }
  else
    {
    if(fctype==MCMC::spatial)
      {
      if(!centervcm)
        {
        betamean=datamatrix(nrpar,1,betareml(betaXpos,0))+remlspatialdesign*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1);
        for(i=0; i<nrpar; i++)
          {
          betastd(i,0) = sqrt(
                              (
                              X_VCM(0,0)*betacov(betaXpos,betaXpos)
                               +
                               (remlspatialdesign.getRow(i)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+1))(0,0)
                              )*X_VCM(0,0)
                              +
                              (
                               (
                                X_VCM(0,0)*betacov.getBlock(betaXpos,betaZpos,betaXpos+1,betaZpos+dimZ)
                                +
                                remlspatialdesign.getRow(i)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                               )*(remlspatialdesign.getRow(i).transposed())
                              )(0,0)
                             );
          }
        }
      else
        {
        betamean=remlspatialdesign*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1);
        for(i=0; i<nrpar; i++)
          {
          betastd(i,0) = sqrt((remlspatialdesign.getRow(i)*
                              betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)*
                              (remlspatialdesign.getRow(i).transposed()))(0,0));
          }
        }
      }
    else if(type==MCMC::RW1)
      {
      for(i=0;i<nrpar;i++)
        {
        j=index(posbeg[i],0);
        if(!centervcm)
          {
          betamean(i,0) = betareml(betaXpos,0)*X_VCM(j,0) + (Z_VCM.getRow(j)*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1))(0,0);
          betastd(i,0) = sqrt(
                              (
                               X_VCM(j,0)*betacov(betaXpos,betaXpos)
                               +
                               (Z_VCM.getRow(j)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+1))(0,0)
                              )*X_VCM(j,0)
                              +
                              (
                               (
                                X_VCM(j,0)*betacov.getBlock(betaXpos,betaZpos,betaXpos+1,betaZpos+dimZ)
                                +
                                Z_VCM.getRow(j)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                               )*(Z_VCM.getRow(j).transposed())
                              )(0,0)
                             );
          }
        else
          {
          betamean(i,0) = (Z_VCM.getRow(j)*betareml.getBlock(betaZpos,0,betaZpos+nrpar-1,1))(0,0);
          betastd(i,0) = sqrt((Z_VCM.getRow(j)*
                              betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)*
                              (Z_VCM.getRow(j).transposed()))(0,0));
          }
        }
      }
    else if(type==MCMC::RW2)
      {
      for(i=0;i<nrpar;i++)
        {
        j=index(posbeg[i],0);
/*        betamean(i,0) = (X_VCM.getRow(j)*betareml.getBlock(Xpos,0,Xpos+dimX,1))(0,0) + (Z_VCM.getRow(j)*betareml.getBlock(X.cols()+Zpos,0,X.cols()+Zpos+dimZ,1))(0,0);
        betastd(i,0) = sqrt(
                            ((
                            X_VCM.getRow(j)*betacov.getBlock(Xpos,Xpos,Xpos+dimX,Xpos+dimX)
                            +
                            (Z_VCM.getRow(j)*betacov.getBlock(X.cols()+Zpos,Xpos,X.cols()+Zpos+dimZ,Xpos+dimX))
                            )*X_VCM.getRow(j).transposed())(0,0)
                            +
                            (
                            (
                             X_VCM.getRow(j)*betacov.getBlock(Xpos,X.cols()+Zpos,Xpos+dimX,X.cols()+Zpos+dimZ)
                             +
                             Z_VCM.getRow(j)*betacov.getBlock(X.cols()+Zpos,X.cols()+Zpos,X.cols()+Zpos+dimZ,X.cols()+Zpos+dimZ)
                            )*(Z_VCM.getRow(j).transposed())
                            )(0,0)
                           );*/
        betamean(i,0) = (X_VCM.getRow(j)*betareml.getBlock(betaXpos,0,betaXpos+dimX,1))(0,0) + (Z_VCM.getRow(j)*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1))(0,0);
        betastd(i,0) = sqrt(
                            ((
                            X_VCM.getRow(j)*betacov.getBlock(betaXpos,betaXpos,betaXpos+dimX,betaXpos+dimX)
                            +
                            (Z_VCM.getRow(j)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+dimX))
                            )*X_VCM.getRow(j).transposed())(0,0)
                            +
                            (
                            (
                             X_VCM.getRow(j)*betacov.getBlock(betaXpos,betaZpos,betaXpos+dimX,betaZpos+dimZ)
                             +
                             Z_VCM.getRow(j)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                            )*(Z_VCM.getRow(j).transposed())
                            )(0,0)
                           );
        }
      }
    else if(type==MCMC::seasonal)
      {
      for(i=0;i<nrpar;i++)
        {
        j=index(posbeg[i],0);
/*        betamean(i,0) = (X_VCM.getRow(j)*betareml.getBlock(Xpos,0,Xpos+dimX,1))(0,0) + (Z_VCM.getRow(j)*betareml.getBlock(X.cols()+Zpos,0,X.cols()+Zpos+dimZ,1))(0,0);
        betastd(i,0) = sqrt(
                           ((
                            X_VCM.getRow(j)*betacov.getBlock(Xpos,Xpos,Xpos+dimX,Xpos+dimX)
                            +
                            (Z_VCM.getRow(j)*betacov.getBlock(X.cols()+Zpos,Xpos,X.cols()+Zpos+dimZ,Xpos+dimX))
                           )*X_VCM.getRow(j).transposed())(0,0)
                           +
                           (
                             (
                             X_VCM.getRow(j)*betacov.getBlock(Xpos,X.cols()+Zpos,Xpos+dimX,X.cols()+Zpos+dimZ)
                             +
                             Z_VCM.getRow(j)*betacov.getBlock(X.cols()+Zpos,X.cols()+Zpos,X.cols()+Zpos+dimZ,X.cols()+Zpos+dimZ)
                           )*(Z_VCM.getRow(j).transposed())
                           )(0,0)
                          );*/
        betamean(i,0) = (X_VCM.getRow(j)*betareml.getBlock(betaXpos,0,betaXpos+dimX,1))(0,0) + (Z_VCM.getRow(j)*betareml.getBlock(betaZpos,0,betaZpos+dimZ,1))(0,0);
        betastd(i,0) = sqrt(
                           ((
                            X_VCM.getRow(j)*betacov.getBlock(betaXpos,betaXpos,betaXpos+dimX,betaXpos+dimX)
                            +
                            (Z_VCM.getRow(j)*betacov.getBlock(betaZpos,betaXpos,betaZpos+dimZ,betaXpos+dimX))
                           )*X_VCM.getRow(j).transposed())(0,0)
                           +
                           (
                             (
                             X_VCM.getRow(j)*betacov.getBlock(betaXpos,betaZpos,betaXpos+dimX,betaZpos+dimZ)
                             +
                             Z_VCM.getRow(j)*betacov.getBlock(betaZpos,betaZpos,betaZpos+dimZ,betaZpos+dimZ)
                           )*(Z_VCM.getRow(j).transposed())
                           )(0,0)
                          );
        }
      }
    for(j=0; j<nrpar; j++)
      {
      betaqu_l1_lower(j,0) = betamean(j,0)+randnumbers::invPhi2(lower1/100)*betastd(j,0);
      betaqu_l1_upper(j,0) = betamean(j,0)+randnumbers::invPhi2(upper2/100)*betastd(j,0);
      betaqu_l2_lower(j,0) = betamean(j,0)+randnumbers::invPhi2(lower2/100)*betastd(j,0);
      betaqu_l2_upper(j,0) = betamean(j,0)+randnumbers::invPhi2(upper1/100)*betastd(j,0);
      }
    }

  if(ismultinomial)
    {
    optionsp->out("  " + title + " (cat."+ST::doubletostring(category,6)+")\n",true);
    }
  else
    {
    optionsp->out("  " + title + "\n",true);
    }
  optionsp->out("\n");
  optionsp->out("\n");

  optionsp->out("  Estimated variance:  "
                + ST::doubletostring(thetareml(thetapos,0),6) + "\n");
  double smoothpar;
  if(dispers==true)
    {
    smoothpar = thetareml(thetareml.rows()-1,0)/thetareml(thetapos,0);
    optionsp->out("  Inverse variance:    "
                  + ST::doubletostring(1/thetareml(thetapos,0),6) + "\n");
    optionsp->out("  Smoothing parameter: "
                + ST::doubletostring(smoothpar,6) + "\n");
    optionsp->out("  (Smoothing parameter = scale / variance)\n");
    optionsp->out("  Degrees of freedom: "
                + ST::doubletostring(thetareml(thetapos,3),6) + "\n");
    }
  else
    {
    smoothpar = 1/thetareml(thetapos,0);
    optionsp->out("  Smoothing parameter: "
                + ST::doubletostring(smoothpar,6) + "\n");
    optionsp->out("  (Smoothing parameter = 1 / variance)\n");
    optionsp->out("  Degrees of freedom: "
                + ST::doubletostring(thetareml(thetapos,3),6) + "\n");
    }
  if(thetareml(thetapos,1)==1)
    {
    optionsp->out("  NOTE: Estimation of the variance was stopped after iteration "
                  + ST::doubletostring(thetareml(thetapos,2),0) + "\n");
    optionsp->out("        because the corresponding penalized part was small relative to the linear predictor.");
    }
  if(thetareml(thetapos,1)==-1)
    {
    optionsp->out("  NOTE: Estimation of the variance was stopped after iteration "
                  + ST::doubletostring(thetareml(thetapos,2),0) + "\n");
    optionsp->out("        to avoid numerical problems due to large variances.");
    }
  ST::string varpath=pathcurrent.substr(0,pathcurrent.length()-4) + "_var.res";
  if(ismultinomial)
    {
    varpath=varpath.insert_after_string(ST::doubletostring(category,6)+"_","_f_");
    }
  optionsp->out("\n");
  optionsp->out("  Variance and smoothing parameter are stored in file\n");
  optionsp->out("  " + varpath + "\n");

  ofstream outvarres(varpath.strtochar());
  outvarres << "variance  ";
  outvarres << "smoothpar  ";
  outvarres << "df  ";
  outvarres << "stopped  " <<endl;

  outvarres << thetareml(thetapos,0) <<"  ";
  outvarres << smoothpar <<"  ";
  outvarres << thetareml(thetapos,3) <<"  ";
  outvarres << (thetareml(thetapos,1)==1);
  outvarres << endl;
  outvarres.close();

  ST::string outest=pathcurrent;
  if(ismultinomial)
    {
    outest = pathcurrent.insert_after_string(ST::doubletostring(category,6)+"_","_f_");
    }
  ofstream outres(outest.strtochar());

  optionsp->out("\n");
  optionsp->out("  Results are stored in file\n");
  optionsp->out("  " +   outest + "\n");
  optionsp->out("\n");

  if (fctype == MCMC::spatial)
    {
    #if defined(JAVA_OUTPUT_WINDOW)

    if (polex == true)
      {
      optionsp->out("  Postscript file is stored in file\n");
      ST::string psfile = outest.substr(0,outest.length()-4) + ".ps";
      optionsp->out("  " + psfile + "\n");
      optionsp->out("\n");
      }

    optionsp->out("  Results may be visualized in BayesX using method 'drawmap'\n");
    optionsp->out("  Type for example: objectname.drawmap " +
    ST::inttostring(plotpos) + "\n");
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
    ST::string psfile = outest.substr(0,outest.length()-4) + ".ps";
    optionsp->out("  " + psfile + "\n");
    optionsp->out("\n");
    optionsp->out("  Results may be visualized in BayesX using method 'plotnonp'\n");
    optionsp->out("  Type for example: objectname.plotnonp " +
    ST::inttostring(plotpos) + "\n");
    optionsp->out("\n");
    #else
    char hchar = '\\';
    ST::string hstring = "/";
    ST::string pathresultsplus = outest.insert_string_char(hchar,hstring);
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
  outres << "pmode   ";
  outres << "ci"  << level1  << "lower   ";
  outres << "ci"  << level2  << "lower   ";
  outres << "std   ";
  outres << "ci"  << level2  << "upper   ";
  outres << "ci"  << level1  << "upper   ";
  outres << "pcat" << level1 << "   ";
  outres << "pcat" << level2 << "   ";
  outres << endl;

  double * workmean = betamean.getV();
  double * workstd = betastd.getV();
  double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
  double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
  double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
  double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();

  vector<ST::string>::iterator effit = effectvalues.begin();

  for(i=0;i<nrpar;i++,++effit,workmean++,workbetaqu_l1_lower_p++,
                           workbetaqu_l2_lower_p++,workstd++,
                           workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
    {
    outres << (i+1) << "   ";
    outres << *effit << "   ";
    outres << *workmean << "   ";
    outres << *workbetaqu_l1_lower_p << "   ";
    outres << *workbetaqu_l2_lower_p << "   ";
    outres << *workstd << "   ";
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
  return mean;
  }

void FULLCOND_nonp_gaussian::outoptionsreml()
  {
  optionsp->out("OPTIONS FOR NONPARAMETRIC TERM: " + title + "\n",true);
  optionsp->out("\n");
  ST::string typestr;
  if (type == RW1)
    typestr = "first order random walk";
  else if (type == RW2)
    typestr = "second order random walk";
  else if (type == seasonal)
    typestr = "seasonal component";
  else if (type==mrf)
    typestr = "spatial Markov random field";

  optionsp->out("  Prior: " + typestr + "\n");
  if (type==seasonal)
    optionsp->out("  Period: " + ST::inttostring(period) + "\n");
  optionsp->out("  Starting value for lambda: " + ST::doubletostring(startlambda,6) + "\n" );
  optionsp->out("\n");
  }

// -------------------------------- END: FOR REML ------------------------------

void FULLCOND_nonp_gaussian::set_lambdaconst(double la)
  {
  lambda=la;
  lambdaconst = true;
  }


void FULLCOND_nonp_gaussian::init_priorassumptions(const ST::string & na)
    {
    if(type==MCMC::RW1)
       priorassumptions.push_back("first order random walk");
    else if(type==MCMC::RW2)
       priorassumptions.push_back("second order random walk");
    else if(type==MCMC::mrf)
       priorassumptions.push_back("Markov random field");
    else if(type==MCMC::RE)
       priorassumptions.push_back("random effect");
    else if(type==MCMC::seasonal)
       priorassumptions.push_back("time varying seasonal component");
    else if(type==MCMC::smoothspline)
       priorassumptions.push_back("smoothing spline");
    else if(type==MCMC::mrfkronecker)
       priorassumptions.push_back("Kronecker product interaction");
    else if(type==MCMC::mrflinear)
       priorassumptions.push_back("2 dimensional first order random walk");
    else if(type==MCMC::mrfkr1)
       priorassumptions.push_back("Kronecker product interaction (RW1*RW1)");
    else if(type==MCMC::mrfkr2)
       priorassumptions.push_back("Kronecker product interaction (RW2*RW2)");
    }

void FULLCOND_nonp_gaussian::init_name(const ST::string & na)
    {
    FULLCOND::init_name(na);

    char hchar = '_';
    ST::string hstring =  "\\_";

    ST::string helpname = na.insert_string_char(hchar,hstring);
    if (type==MCMC::seasonal)
      term_symbolic = "f^{Season}_{" +  helpname + "}("+helpname+")";
    else
      term_symbolic = "f_{" +  helpname + "}("+helpname+")";

    if (column > 0)
      priorassumptions.push_back("$" + term_symbolic + "$" +
       " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");
    init_priorassumptions(na);
    }

void FULLCOND_nonp_gaussian::init_names(const vector<ST::string> & na)
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
      ST::string helpname1 = na[0].insert_string_char(hchar,hstring);
      ST::string helpname2 = na[1].insert_string_char(hchar,hstring);
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





// additive Effekte, RW1 RW2 und season

FULLCOND_nonp_gaussian::FULLCOND_nonp_gaussian(MCMCoptions * o,
                      DISTRIBUTION * dp,
                      const datamatrix & d,
                      FULLCOND_const * fcc,
                      const unsigned & maxint,const fieldtype & ft,
                      const ST::string & ti,
                      const ST::string & fp, const ST::string & pres,
                      const unsigned & c,const double & l,
                      const unsigned & per)
  : FULLCOND_nonp_basis(o,dp,ft,ti,fp,pres,c,per)
  {

  notransform=false;

  lambdaconst=false;
  Laplace=false;

  fcconst = fcc;

  if (ft == seasonal)
    fctype = season;
  else
    fctype = nonparametric;

  utype = gaussian;

  lattice=false;

  plotstyle = plotnonp;
  pathresult = pres;

  lambda = l;

  transform = likep->get_trmult(c);

  updatelinpred = true;

  make_categories(d,maxint);

  if (errors.size() == 0)
    {

    if (type == RW1)
      {
      Kenv = Krw1env(weight);
      rankK = Kenv.getDim()-1;
      identifiable = false;
      }
    else if (type == RW2)
      {
      Kenv = Krw2env(weight);
      rankK = Kenv.getDim()-2;
      identifiable = false;
      }
    else if (type == seasonal)
      {
      Kenv = Kseasonenv(period,weight.size());
      rankK = Kenv.getDim()-period+1;
      identifiable = true;
      }

    setbeta(Kenv.getDim(),1,0);


    XXenv = envmatdouble(0,nrpar);
    compute_XWX_env(likep->get_weight());

    precenv = envmatdouble(0,nrpar,Kenv.getBandwidth());

    mu = datamatrix(likep->get_nrobs(),1,0);
    muy = datamatrix(nrpar,1);
    betahelp = muy;

    precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    lambda_prec=lambda;

    varcoeff = false;

    } // end: errors.size() == 0

  }


void FULLCOND_nonp_gaussian::init_data_varcoeff(const datamatrix & intvar,
   double add)
  {

  data = datamatrix(intvar.rows(),1);
  data2 = datamatrix(intvar.rows(),1);

  unsigned i;
  double * workdata = data.getV();
  double * workdata2 = data2.getV();
  int * workindex = index.getV();
  for (i=0;i<data.rows();i++,workdata++,workdata2++,workindex++)
    {
    *workdata = intvar(*workindex,0)+add;
    *workdata2 = (*workdata) * (*workdata);
    }

  }


// varying coefficients , RW1 RW2 und season

FULLCOND_nonp_gaussian::FULLCOND_nonp_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                                     const datamatrix & d,
                                     const datamatrix & intvar,
                                     FULLCOND_const * fcc,
                                     const unsigned & maxint,
                                     const fieldtype & ft,const ST::string & ti,
                               const ST::string & fp, const ST::string & pres,
                               const unsigned & c,const double & l,bool ce,
                               const unsigned & per)
  : FULLCOND_nonp_basis(o,dp,ft,ti,fp,pres,c,per)

  {

  notransform=false;

  lambdaconst=false;
  Laplace=false;

  fcconst = fcc;

  if (ft == seasonal)
    fctype = season;
  else
    fctype = nonparametric;

  lattice=false;

  utype = gaussian;

  plotstyle = plotnonp;
  lambda = l;

  transform = likep->get_trmult(c);

  updatelinpred = true;

  make_categories(d,maxint);

  init_data_varcoeff(intvar);

  if (errors.size() == 0)
    {

    if (ce==false)
      identifiable = true;
    else
      identifiable = false;

    if (type == RW1)
      {
      Kenv = Krw1env(weight);
      rankK = Kenv.getDim()-1;
      }
    else if (type == RW2)
      {
      Kenv = Krw2env(weight);
      rankK = Kenv.getDim()-2;
      }
    else if (type == seasonal)
      {
      Kenv = Kseasonenv(period,weight.size());
      rankK = Kenv.getDim()-period+1;
      }

    setbeta(Kenv.getDim(),1,0);

    XXenv = envmatdouble(0,nrpar);
    compute_XWX_varcoeff_env(likep->get_weight());

    precenv = envmatdouble(0,nrpar,Kenv.getBandwidth());

    mu = datamatrix(likep->get_nrobs(),1,0);
    muy = datamatrix(nrpar,1);
    betahelp = muy;

    precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    lambda_prec=lambda;

    varcoeff = true;

    }
  }


// spatial covariates

FULLCOND_nonp_gaussian::FULLCOND_nonp_gaussian(MCMCoptions * o,
                        DISTRIBUTION * dp,const datamatrix & d,
                        FULLCOND_const * fcc,
                        const MAP::map & m, const ST::string & mn,
                        const ST::string & ti,
                        const ST::string & fp, const ST::string & pres,
                        const unsigned & c,const double & l)
  : FULLCOND_nonp_basis(o,dp,mrf,ti,fp,pres,c)

  {

  notransform=false;

  lambdaconst=false;
  Laplace=false;

  MAP::map ma=m;

  fcconst = fcc;

  fctype = spatial;

  utype = gaussian;

  lattice=false;

  mapname = mn;

  lambda = l;

  transform = likep->get_trmult(c);

  updatelinpred = true;

  unsigned i;

  if (ma.polygones_existing())
    polex = true;
  else
    polex = false;

  if(polex == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

  if (ma.get_bandsize() > 40)
    ma.reorderopt();

  ma.compute_reg(d,posbeg,posend,effectvalues,index);

  if (ma.get_errormessages().size() > 0)
    errors = ma.get_errormessages();
  else
    {

    Kenv = Kmrfenv(ma);

    rankK = Kenv.getDim()-1;

    setbeta(Kenv.getDim(),1,0);

    identifiable = false;
    varcoeff = false;

    XXenv = envmatdouble(0,nrpar);
    compute_XWX_env(likep->get_weight());

    precenv = envmatdouble(Kenv.getXenv(),0,nrpar);

    precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    lambda_prec=lambda;

    mu = datamatrix(likep->get_nrobs(),1,0);
    muy = datamatrix(nrpar,1);
    betahelp = muy;

    for (i=0;i<posbeg.size();i++)
      {
      if (posbeg[i] == -1)
        optionsp->out("NOTE: no observations for region " + effectvalues[i] + "\n");
      }


    // NIKI

/*    datamatrix Kstat=STATMAT_PENALTY::Kmrf(ma);
    datamatrix vals(Kstat.rows(),1,0);

    bool eigentest=eigen2(Kstat,vals);

    ofstream out("d:\\temp\\Kstat.raw");
    Kstat.prettyPrint(out);

    ofstream out2("d:\\temp\\vals.raw");
    vals.prettyPrint(out2);
*/
    // NIKI


    } // end: if (error==false)

  unsigned j;
  neighbors=vector< vector<unsigned> >(nrpar);
  for(i=0;i<nrpar;i++)
    {
    for(j=0;j<nrpar;j++)
      {
      if(Kenv(i,j)!=0.0 && i!=j)
        {
        neighbors[i].push_back(j);
        }
      }
    }

  }


// varying coefficients , spatial covariates as effect modifier

FULLCOND_nonp_gaussian::FULLCOND_nonp_gaussian(MCMCoptions * o,
                        DISTRIBUTION * dp,
                        FULLCOND_const * fcc,
                        const MAP::map & m,
                        const ST::string & mn,
                        const datamatrix & d,
                        const datamatrix & d2,
                        const ST::string & ti,
                        const ST::string & fp, const ST::string & pres,
                        const unsigned & c,const double & l,bool ce)
  : FULLCOND_nonp_basis(o,dp,mrf,ti,fp,pres,c)

  {



  notransform=false;
  lambdaconst=false;
  Laplace=false;

  MAP::map ma = m;

  fcconst = fcc;

  fctype = spatial;

  utype = gaussian;

  lattice=false;

  lambda = l;
  mapname = mn;

  transform = likep->get_trmult(c);

  updatelinpred = true;

  unsigned i;

  if (ma.polygones_existing())
    polex = true;
  else
    polex = false;

  if(polex == true)
    plotstyle = drawmap;
  else
    plotstyle = drawmapgraph;

  if (ma.get_bandsize() > 40)
    ma.reorderopt();

  ma.compute_reg(d,posbeg,posend,effectvalues,index);

  init_data_varcoeff(d2);

  if (ma.get_errormessages().size() > 0)
    errors = ma.get_errormessages();
  else
    {

    Kenv = Kmrfenv(ma);

    rankK = Kenv.getDim()-1;

    setbeta(Kenv.getDim(),1,0);

    if (ce==false)
      identifiable = true;
    else
      identifiable = false;  

    varcoeff = true;

    XXenv = envmatdouble(0,nrpar);
    compute_XWX_varcoeff_env(likep->get_weight());

    precenv = envmatdouble(Kenv.getXenv(),0,nrpar);

    precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    lambda_prec=lambda;

    mu = datamatrix(likep->get_nrobs(),1,0);
    muy = datamatrix(nrpar,1);
    betahelp=muy;

    for (i=0;i<posbeg.size();i++)
      {
      if (posbeg[i] == -1)
        optionsp->out("NOTE: no observations for region " + effectvalues[i] + "\n");
      }

    } // end: if (error == false)

  unsigned j;
  neighbors=vector< vector<unsigned> >(nrpar);
  for(i=0;i<nrpar;i++)
    {
    for(j=0;j<nrpar;j++)
      {
      if(Kenv(i,j)!=0.0 && i!=j)
        {
        neighbors[i].push_back(j);
        }
      }
    }

  }



// spatial effect with x and y
FULLCOND_nonp_gaussian::FULLCOND_nonp_gaussian(MCMCoptions * o,
                         DISTRIBUTION * dp,const datamatrix & dx,
                         const datamatrix & dy,
                         FULLCOND_const * fcc,
                         const double & l,
                         const double & md,const ST::string & mp,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & pres, const ST::string & pmap,
                         const unsigned & c)
  : FULLCOND_nonp_basis(o,dp,mrf,ti,fp,pres,c)
  {

  notransform=false;

  lambdaconst=false;
  Laplace=false;

  fcconst = fcc;

  fctype = spatial;

  utype = gaussian;

  lattice = true;

  lambda = l;

  unsigned i;

  datamatrix X(dx.rows(),2);
  double * workx = dx.getV();
  double * worky = dy.getV();
  double * workX = X.getV();
  for (i=0;i<dx.rows();i++,workx++,worky++)
    {
    *workX = *workx;
    workX++;
    *workX = *worky;
    *workX++;
    }

  MAP::map m(
  #if defined(JAVA_OUTPUT_WINDOW)
  optionsp->adminb_p,
  #endif
  X,md,MAP::adjacent);

  if (m.isconnected()==false)
    errors.push_back("ERROR: maxdist=" + ST::doubletostring(md) +
                     " leads to a disconnected graph\n");
  if (m.get_nrregions() <=3)
    errors.push_back("ERROR: not enough regions to estimate spatial effect\n");

  if (errors.size() == 0)
    {

    m.reorderopt();

    m.outmap(pmap);
    pathmap = pmap;

    ST::string pathgraph = pathmap.substr(0, pathmap.length()-4);
    pathgraph = pathgraph + "gra";
    m.outgraph(pathgraph);


    mapname = mp;

    updatelinpred = true;

    transform = likep->get_trmult(c);

    // Initialization of the index, Indexsort of moddata
    index = statmatrix<int>(dx.rows(),1);
    index.indexinit();
    dx.indexsort(index,0,dx.rows()-1,0,0);

    unsigned beg = 0;
    unsigned end;
    for (i=1;i<dx.rows();i++)
      {
      if ( dx(index(i,0),0) != dx(index(i-1,0),0) )
        {
        end = i-1;
        if (end-beg > 0)
          dy.indexsort(index,beg,end,0,0);

        beg = i;
        }
      else if (i==dx.rows()-1)
        {
        end = i;
        if (end-beg > 0)
          dy.indexsort(index,beg,end,0,0);
        }

      }

    // end: indexsort

    datamatrix d(dx.rows(),1);


    double n = 1;
    d(index(0,0),0) = n;
    for (i=1;i<d.rows();i++)
      {
      if ( (dx(index(i,0),0) != dx(index(i-1,0),0)) ||
           (dy(index(i,0),0) != dy(index(i-1,0),0))
         )
        n++;
      d(index(i,0),0) = n;
      }

    m.compute_reg(d,posbeg,posend,effectvalues,index);

    nrpar = posbeg.size();

    xyvalues = datamatrix(nrpar,2);
    for (i=0;i<nrpar;i++)
      {
      xyvalues(i,0) = dx(index(posbeg[i],0),0);
      xyvalues(i,1) = dy(index(posbeg[i],0),0);
      }

    Kenv = Kmrfenv(m);

    rankK = Kenv.getDim()-1;

    setbeta(Kenv.getDim(),1,0);

    identifiable = false;
    varcoeff = false;

    XXenv = envmatdouble(0,nrpar);
    compute_XWX_env(likep->get_weight());

    precenv = envmatdouble(Kenv.getXenv(),0,nrpar);

    precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    lambda_prec=lambda;

    mu = datamatrix(likep->get_nrobs(),1,0);
    muy = datamatrix(nrpar,1);
    betahelp = muy;
    }

  }


FULLCOND_nonp_gaussian::FULLCOND_nonp_gaussian(const FULLCOND_nonp_gaussian & fc)
  : FULLCOND_nonp_basis(FULLCOND_nonp_basis(fc))
  {
  notransform = fc.notransform;
  lambdaconst = fc.lambdaconst;
  Laplace=fc.Laplace;
  delta=fc.delta;
  neighbors=fc.neighbors;
  diff = fc.diff;
  betaKbeta = fc.betaKbeta;
  data2 = fc.data2;
  betaold = fc.betaold;
  tildey = fc.tildey;
  weightiwls = fc.weightiwls;
  a_invgamma=fc.a_invgamma;
  b_invgamma=fc.b_invgamma;
  oldacceptance = fc.oldacceptance;
  oldnrtrials=fc.oldnrtrials;
  lambdaprop = fc.lambdaprop;
  lambda_prec = fc.lambda_prec;
  f=fc.f;
  utype = fc.utype;
  updateW=fc.updateW;
  fcconst = fc.fcconst;
  updatelinpred = fc.updatelinpred;
  mu = fc.mu;
  muy = fc.muy;
  betahelp = fc.betahelp;
  betamode = fc.betamode;
  betamodeold = fc.betamodeold;
  XXenv = fc.XXenv;
  precenv = fc.precenv;
  mapname = fc.mapname;
  lattice = fc.lattice;
  xyvalues = fc.xyvalues;
  pathmap = fc.pathmap;
  X_VCM=fc.X_VCM;
  Z_VCM=fc.Z_VCM;
  remlspatialdesign=fc.remlspatialdesign;
  }


const FULLCOND_nonp_gaussian & FULLCOND_nonp_gaussian::operator=(
                                            const FULLCOND_nonp_gaussian & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_nonp_basis::operator=(FULLCOND_nonp_basis(fc));
  notransform = fc.notransform;
  lambdaconst = fc.lambdaconst;
  Laplace=fc.Laplace;
  delta=fc.delta;
  neighbors=fc.neighbors;
  diff = fc.diff;
  betaKbeta = fc.betaKbeta;
  data2 = fc.data2;
  betaold = fc.betaold;
  tildey = fc.tildey;
  weightiwls = fc.weightiwls;
  a_invgamma=fc.a_invgamma;
  b_invgamma=fc.b_invgamma;
  oldacceptance = fc.oldacceptance;
  oldnrtrials=fc.oldnrtrials;
  lambdaprop = fc.lambdaprop;
  lambda_prec = fc.lambda_prec;
  f=fc.f;
  utype = fc.utype;
  updateW=fc.updateW;
  fcconst = fc.fcconst;
  updatelinpred = fc.updatelinpred;
  mu = fc.mu;
  muy = fc.muy;
  betahelp = fc.betahelp;
  betamode = fc.betamode;
  betamodeold = fc.betamodeold;
  XXenv = fc.XXenv;
  precenv = fc.precenv;
  mapname = fc.mapname;
  lattice = fc.lattice;
  xyvalues = fc.xyvalues;
  pathmap = fc.pathmap;
  X_VCM=fc.X_VCM;
  Z_VCM=fc.Z_VCM;
  remlspatialdesign=fc.remlspatialdesign;
  return *this;
  }


void FULLCOND_nonp_gaussian::make_categories(const datamatrix & moddata,
                                             const unsigned & maxnrint)
  {

  unsigned j;                              // Loopingvariable

  // Initialization of the index, Indexsort of moddata

  index = statmatrix<int>(moddata.rows(),1);
  index.indexinit();
  moddata.indexsort(index,0,moddata.rows()-1,0,0);

  double diff;                             // difference of two succeeding
                                           // observations
  double maxdiff = 1.0/double(maxnrint);
  double beg;                              // first observation in the
                                           // current intervall
  double last;                             // (j-1). observation
  double intbeg;

  double distance =
  moddata(index(moddata.rows()-1,0),0) - moddata(index(0,0),0);

  if (distance == 0)
    errors.push_back("ERROR: not enough different covariate values (" + title +
                     ")\n");
  else
    {

    beg = moddata(index(0,0),0);
    last = moddata(index(0,0),0);
    intbeg = moddata(index(0,0),0);

    posbeg.push_back(0);

    for(j=1;j<moddata.rows();j++)
      {


      diff = (moddata(index(j,0),0)-beg)/distance;

      if (diff > maxdiff)
        {
        weight.push_back(last - intbeg);
        effectvalues.push_back(ST::doubletostring(last));
        effectvdouble.push_back(last);
        posbeg.push_back(j);
        posend.push_back(j-1);

        beg = moddata(index(j,0),0);
        intbeg = last;
        }

      if (j == moddata.rows()-1)
        {
        weight.push_back(moddata(index(j,0),0)-intbeg);
        effectvalues.push_back((ST::doubletostring(moddata(index(j,0),0))));
        effectvdouble.push_back(moddata(index(j,0),0));
        }

      last = moddata(index(j,0),0);


      }  // end: for(j=1;j<moddata.rows();j++)

    posend.push_back(moddata.rows()-1);

    double sum=0;
    for(j=1;j<weight.size();j++)
      sum+= weight[j];

    if (type==RW1)
      sum = double(weight.size()-1)/sum;
    else  // RW2
      sum = 0.5*double(weight.size()-1)/sum;

    for(j=1;j<weight.size();j++)
      weight[j] *= sum;

    if (posbeg.size() < 6)
      errors.push_back("ERROR: not enough different covariate values (" +
                       title + ")\n");
    }


  }  // end: function make_categories


void FULLCOND_nonp_gaussian::compute_muy(double * workbeta)
  {
  int * workindex = index.getV();
  double * workmuy = muy.getV();
  unsigned i;
  int j;

  if (varcoeff)
    {
    double * workdata=data.getV();
    for(i=0;i<nrpar;i++,workmuy++,workbeta++)
      {
      *workmuy = 0;
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++,workdata++)
          *workmuy+=   weightiwls(*workindex,0) *
          (tildey(*workindex,0) + *workbeta) * (*workdata);
      }
    }
  else  // else additive
    {

    for(i=0;i<nrpar;i++,workmuy++,workbeta++)
      {
      *workmuy = 0;
      if (posbeg[i] != -1)
        for(j=posbeg[i];j<=posend[i];j++,workindex++)
          *workmuy+=  weightiwls(*workindex,0)
                      * (tildey(*workindex,0) + *workbeta);
      }

    }

  }


void FULLCOND_nonp_gaussian::compute_XWX_XWtildey_env(
const datamatrix & weightmat,const datamatrix & tildey,double * workbeta,
const unsigned & col)
  {
  unsigned i;
  int j;

  int *  workindex = index.getV();
  double * workmuy = muy.getV();
  vector<double>::iterator d = XXenv.getDiagIterator();
  double cw;

  for(i=0;i<nrpar;i++,++d,workmuy++,workbeta++)
    {
    *d=0;
    *workmuy=0;
    if (posbeg[i] != -1)
      {
      for (j=posbeg[i];j<=posend[i];j++,workindex++)
        {
        cw =weightmat(*workindex,col);
        *d += cw;
        *workmuy+=  cw * (tildey(*workindex,0) + *workbeta);
        }
      }

    }

  }


void FULLCOND_nonp_gaussian::compute_XWX_XWtildey_varcoeff_env(
  const datamatrix & weightmat,const datamatrix & tildey,double * workbeta,
  const unsigned & col)
  {

  unsigned i;
  int j;

  int *  workindex = index.getV();
  double * workmuy = muy.getV();
  double * workdata=data.getV();
  double * workdata2=data2.getV();
  vector<double>::iterator d = XXenv.getDiagIterator();
  double cw;

  for(i=0;i<nrpar;i++,++d,workmuy++,workbeta++)
    {
    *d=0;
    *workmuy=0;
    if (posbeg[i] != -1)
      {
      for (j=posbeg[i];j<=posend[i];j++,workindex++,workdata++,workdata2++)
        {
        cw =weightmat(*workindex,col);
        *d += cw* (*workdata2);
        *workmuy+=  cw * (tildey(*workindex,0) + *workbeta) * (*workdata);
        }
      }

    }

  }



void FULLCOND_nonp_gaussian::compute_XWX_env(const datamatrix & weightmat,
const unsigned & col)
  {
  unsigned i;
  int j;
  int *  workindex = index.getV();

  vector<double>::iterator d = XXenv.getDiagIterator();

  for(i=0;i<posbeg.size();i++,++d)
    {
    *d=0;
    if (posbeg[i] != -1)
      {

      for (j=posbeg[i];j<=posend[i];j++,workindex++)
        *d += weightmat(*workindex,col);
      }

    }

  }



void FULLCOND_nonp_gaussian::compute_XWX_varcoeff_env
                             (const datamatrix & weightmat,const unsigned & col)
  {
  int j;
  unsigned i;
  int *  workindex = index.getV();

  vector<double>::iterator d = XXenv.getDiagIterator();
  double * workdata=data.getV();
  double * workdata2=data2.getV();

  for(i=0;i<posbeg.size();i++,++d)
    {
    *d=0;
    if (posbeg[i] != -1)
      {

      for (j=posbeg[i];j<=posend[i];j++,workindex++,workdata++,workdata2++)
        {

        *d += weightmat(*workindex,col)*(*workdata2);
        }
      }

    }

  }


void FULLCOND_nonp_gaussian::update_linpred_diff(datamatrix & b1,
                                                 datamatrix & b2)
  {

  int * workindex;
  double * workb1 = b1.getV();
  double * workb2 = b2.getV();
  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();

  int j;
  unsigned i;

  if (varcoeff)
    {

    workindex = index.getV();
    double * workdata=data.getV();

    for (i=0;i<nrpar;i++,workb1++,workb2++,++itbeg,++itend)
      {
      if (*itbeg != -1)
        {
        for(j=*itbeg;j<=*itend;j++,workindex++,workdata++)
          {
          likep->add_linearpred((*workb1-*workb2) * (*workdata),
                                unsigned(*workindex),column);
          }
        }
      }

    }
  else
    {
    for (i=0;i<nrpar;i++,workb1++,workb2++,++itbeg,++itend)
      {
      if (*itbeg != -1)
        likep->add_linearpred(*workb1-*workb2,*itbeg,*itend,index,column);
      }
    }

  }


void FULLCOND_nonp_gaussian::update_linpred_diff(const unsigned & beg,const unsigned & end,
                                                 const double & beta)
  {

  int * workindex;
  int j;

  if (varcoeff)
    {

    workindex = index.getV()+beg;
    double * workdata=data.getV()+beg;

    for(j=beg;j<=end;j++,workindex++,workdata++)
      {
      likep->add_linearpred((beta) * (*workdata),
                                unsigned(*workindex),column);
      }

    }
  else
    {
    likep->add_linearpred(beta,beg,end,index,column);
    }

  }


void FULLCOND_nonp_gaussian::update_linpred(const bool & add)
  {

  int * workindex;
  double * workbeta = beta.getV();
  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();

  int j;
  unsigned i;

  if (add==false)
    {

    if (varcoeff)
      {
      workindex = index.getV();
      double * workdata=data.getV();
      for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
        {
        if (*itbeg != -1)
          {
          for(j=*itbeg;j<=*itend;j++,workindex++,workdata++)
            {
            likep->add_linearpred(-*workbeta*(*workdata),unsigned(*workindex),
            column);
            }
          }
        }

      }
    else
      {
      for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
        {
        if (*itbeg != -1)
          likep->add_linearpred(-*workbeta,*itbeg,*itend,index,column);
        }
      }

    } // end: if (add==false)
  else
    {   // add==true

    if (varcoeff)
      {

      workindex = index.getV();
      double * workdata=data.getV();

      for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
        {
        if (*itbeg != -1)
          {
          for(j=*itbeg;j<=*itend;j++,workindex++,workdata++)
            {
            likep->add_linearpred(*workbeta*(*workdata),unsigned(*workindex),
            column);
            }
          }
        }

      }
    else
      {
      for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
        {
        if (*itbeg != -1)
          likep->add_linearpred(*workbeta,*itbeg,*itend,index,column);
        }
      }

    }   // END: add==true


  }


void FULLCOND_nonp_gaussian::update_linpred_current(const bool & add)
  {

  int * workindex;
  double * workbeta = beta.getV();
  double * workbetaold = betaold.getV();
  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
  int j;
  unsigned i;

  if (add==false)                   // substract
    {

    if (varcoeff)
      {
      workindex = index.getV();
      double * workdata=data.getV();
      for (i=0;i<nrpar;i++,workbeta++,workbetaold++,++itbeg,++itend)
        {
        if (*itbeg != -1)
          {
          for(j=*itbeg;j<=*itend;j++,workindex++,workdata++)
            {
            likep->addtocurrentcol_single((*workbetaold-*workbeta)*(*workdata),*workindex,column);
            }
          }
        }

      }
    else
      {
      for (i=0;i<nrpar;i++,workbeta++,workbetaold++,++itbeg,++itend)
        {
        if (*itbeg != -1)
          likep->addtocurrentcol(*workbetaold-*workbeta,*itbeg,*itend,index,column);
        }
      }

    } // end: if (add==false)
  else
    {   // add==true

    if (varcoeff)
      {

      workindex = index.getV();
      double * workdata=data.getV();

      for (i=0;i<nrpar;i++,workbeta++,workbetaold++,++itbeg,++itend)
        {
        if (*itbeg != -1)
          {
          for(j=*itbeg;j<=*itend;j++,workindex++,workdata++)
            {
            likep->addtocurrentcol_single((*workbeta-*workbetaold)*
            (*workdata),*workindex,column);
            }
          }
        }

      }
    else
      {
      for (i=0;i<nrpar;i++,workbeta++,workbetaold++,++itbeg,++itend)
        {
        if (*itbeg != -1)
          likep->addtocurrentcol(*workbeta-*workbetaold,*itbeg,*itend,index,column);
        }
      }

    }   // END: add==true


  }


void FULLCOND_nonp_gaussian::set_IWLS(const unsigned & uw,bool mode)
  {
  if (mode)
    {
    utype=iwlsmode;
    diff = beta;
    betamodeold = beta;
    betamode=beta;
    }
  else
    utype=iwls;

  betaold = beta;
  updateW = uw;
  tildey=datamatrix(likep->get_nrobs(),1);
  weightiwls=datamatrix(likep->get_nrobs(),1,0);
  }


void FULLCOND_nonp_gaussian::set_IWLS_hyperblock(const unsigned & uw,
                                      const double & ai,const double & bi,
                                      bool mode)
  {
  f = 2;
  if (mode)
    {
    utype=hyperblockmode;
    diff = beta;
    betamodeold = beta;
    betamode=beta;
    }
  else
    utype=hyperblock;

  betaold = beta;
  updateW = uw;
  a_invgamma = ai;
  b_invgamma = bi;
  lambda = 1/sigma2;
  oldacceptance = 0;
  oldnrtrials = 0;
  tildey=datamatrix(likep->get_nrobs(),1);
  weightiwls=datamatrix(likep->get_nrobs(),1,0);
  }


double FULLCOND_nonp_gaussian::scale_proposal()
  {
  double len = f - 1/f;
  if (f == 1.0)
    return 1.0;
  if (uniform() < len/(len+2*log(f)))
    return (1/f + len*uniform());
  else
    return pow(f, 2.0*uniform()-1.0);
  }


void FULLCOND_nonp_gaussian::update_IWLS_hyperblock(void)
  {

  unsigned i;
  double * workbeta;

  betaold.assign(beta);

  // Compute log-likelihood with old beta, compute weights, tildey

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
  logold -= 0.5*Kenv.compute_quadform(beta,0)*lambda;
  logold += 0.5*rankK*log(lambda);
  logold += (a_invgamma-1)*log(lambda) - b_invgamma*lambda;

  // updating the tuning constant f

  if(optionsp->get_nriter()<optionsp->get_burnin() &&
     optionsp->get_nriter()%100==0)
    tune_updatetau(alpha_50);

  if(optionsp->get_nriter() == optionsp->get_burnin())
    optionsp->out("  NOTE: Tuning constant 'f' for term "
    + title + " set to " + ST::doubletostring(f) + "\n");  // lambda ~ (1+1/z)

  // computing new proposal for lambda

  lambdaprop = lambda*scale_proposal();

  // Compute XWX and XW(tildey+f)

  workbeta = betaold.getV();

  if (  (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) )
      )
    {
    if (varcoeff)
      compute_XWX_XWtildey_varcoeff_env(weightiwls,tildey,workbeta,0);
    else
      compute_XWX_XWtildey_env(weightiwls,tildey,workbeta,0);
    }
  else
    {
    compute_muy(workbeta);
    }

  // Compute iwls proposal

  precenv.addtodiag(XXenv,Kenv,1.0,lambdaprop);
  precenv.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  precenv.solveU(beta,betahelp);

  /// compute -1/2 (beta^p-m(beta^c))' P(beta^c,lambdaprop) (beta^p-m(beta^c))
  // corresponds to log(q(lambda^c,beta^c,lambdaprop,beta^p))

  betahelp.minus(beta,betahelp);
  double qold = 0.5*(precenv.getLogDet()-precenv.compute_quadform(betahelp,0));

  // Compute predictor based on proposed beta

  update_linpred_diff(beta,betaold);

  // Compute new log-likelihood, weights, tildey based on lambdaprop and beta^p

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

  lognew -= 0.5*Kenv.compute_quadform(beta,0)*lambdaprop;
  lognew += 0.5*(rankK)*log(lambdaprop);
  lognew += (a_invgamma-1)*log(lambdaprop) - b_invgamma*lambdaprop;

  // Compute XWX and XW(tildey+f) based on new beta^p

  workbeta = beta.getV();

  if (  (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) )
      )
    {
    if (varcoeff)
      compute_XWX_XWtildey_varcoeff_env(weightiwls,tildey,workbeta,0);
    else
      compute_XWX_XWtildey_env(weightiwls,tildey,workbeta,0);
    }
  else
    {
    compute_muy(workbeta);
    }

  precenv.addtodiag(XXenv,Kenv,1.0,lambda);
  precenv.solve(muy,betahelp);

  /// compute -1/2 (beta^c-m(beta^p))' P(beta^p,lambda) (beta^c-m(beta^p))
  // corresponds to log(q(lambdaprop,beta^p,lambda^c,beta^c))

  betahelp.minus(betaold,betahelp);
  double qnew = 0.5*precenv.getLogDet() -
         0.5*precenv.compute_quadform(betahelp,0);

  double u = log(uniform());
  if (u <= (lognew - logold  + qnew - qold) )
    {
    lambda = lambdaprop;
    sigma2 = 1.0/lambda;
    acceptance++;

    if (center)
      {
      double m = centerbeta();
      if (varcoeff)
        fcconst->update_fix_varcoeff(m,datanames[1]);
      else
        fcconst->update_intercept(m);
      }

    }
  else
    {

    update_linpred_diff(betaold,beta);

    beta.assign(betaold);
    }


  if (notransform==false)
    transform = likep->get_trmult(column);
  else
    transform = 1;


  FULLCOND::update();


  }


void FULLCOND_nonp_gaussian::update_IWLS_hyperblock_mode(void)
  {

  unsigned i;

  double * workbeta;


  if (optionsp->get_nriter()==1)
    {
    betamode.assign(beta);
    betaold.assign(beta);
    betaKbeta =Kenv.compute_quadform(beta,0);
    }

  // Tuning of constant f

  if (optionsp->get_nriter()<optionsp->get_burnin() &&
      optionsp->get_nriter() % 100==0)
    tune_updatetau(alpha_50);

  if(optionsp->get_nriter() == optionsp->get_burnin())
    optionsp->out("  NOTE: Tuning constant 'f' for term "
    + title + " set to " + ST::doubletostring(f) + "\n");  // lambda ~ (1+1/z)


  // Compute proposal for lambda=1/sigma2

  lambdaprop = lambda*scale_proposal();

  // Compute log-likelihood with old beta and old lambda:
  // log(l)-0.5*lambda*beta^c'K beta^c

  double logold = likep->loglikelihood();
  logold -= 0.5*betaKbeta*lambda;
  logold += 0.5*(nrpar-1)*log(lambda);       // normalising constant
  logold += (a_invgamma-1)*log(lambda)       // gamma prior
             - b_invgamma*lambda;

  // Compute new weights and tildey's based on old betamode

  betamodeold.assign(betamode);
  update_linpred_diff(betamodeold,beta);

  if (  (optionsp->get_nriter() < optionsp->get_burnin())  ||
        ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) )
     )
    {
    likep->compute_IWLS_weight_tildey(weightiwls,tildey,column);
    workbeta = betamodeold.getV();
    if (varcoeff)
      compute_XWX_XWtildey_varcoeff_env(weightiwls,tildey,workbeta,0);
    else
      compute_XWX_XWtildey_env(weightiwls,tildey,workbeta,0);
    }
  else
    {
    likep->tilde_y_minus_eta(tildey,column);
    workbeta = betamodeold.getV();
    compute_muy(workbeta);
    }

  // Compute new betamode

  precenv.addtodiag(XXenv,Kenv,1.0,lambdaprop);
  precenv.solve(muy,betamode);

  // compute new proposal for beta conditional to the proposed lambdaprop

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  precenv.solveU(beta,betamode);

  // compute -1/2 (beta^p-betamode)'P(lambdaprop) (beta^p-betamode)
  // corresponds to log(q(beta^c,beta^p))

  diff.minus(beta,betamode);
  double qnewbeta = 0.5*(precenv.getLogDet()-precenv.compute_quadform(diff,0));

  // compute -1/2 (beta^c-betamode)'P(lambda) (beta^c-betamode)
  // corresponds to log(q(beta^c,beta^p))

  precenv.addtodiag(XXenv,Kenv,1.0,lambda);

  diff.minus(betaold,betamode);
  double qoldbeta = 0.5*(precenv.getLogDet()-precenv.compute_quadform(diff,0));

  // compute new predictor with proposed beta

  update_linpred_diff(beta,betamodeold);

  // Compute log-likelihood with proposed beta^p:
  // log(l)-1/(2tau^2)*beta^p'K beta^p

  double lognew = likep->loglikelihood();
  lognew  -= 0.5*Kenv.compute_quadform(beta,0)*lambdaprop;
  lognew += 0.5*(nrpar-1)*log(lambdaprop);      // normalising constant
  lognew += (a_invgamma-1)*log(lambdaprop)
             - b_invgamma*lambdaprop;           // gamma prior

  // Accept/reject proposed beta and lambda

  double u = log(uniform());
  if (u <= (lognew - logold  + qoldbeta - qnewbeta) )
    {
    acceptance++;

    lambda = lambdaprop;
    sigma2 = 1.0/lambda;


  if (center)
    {
    double m = centerbeta();
    if (varcoeff)
      fcconst->update_fix_varcoeff(m,datanames[1]);
    else
      fcconst->update_intercept(m);
    }


    betaKbeta=Kenv.compute_quadform(beta,0);
    betaold.assign(beta);
    }
  else
    {
    update_linpred_diff(betaold,beta);
    beta.assign(betaold);
    }

  if (notransform==false)
    transform = likep->get_trmult(column);
  else
    transform = 1;

  FULLCOND::update();

  }


void FULLCOND_nonp_gaussian::update_IWLS_mode(void)
  {

  unsigned i;

  double * workbeta;

  if (optionsp->get_nriter()==1)
    {
    betamode.assign(beta);
    betaold.assign(beta);
    betaKbeta =Kenv.compute_quadform(beta,0);
    }

  if (!lambdaconst)
    lambda = 1.0/sigma2;

  // Compute log-likelihood with old beta: log(l)-1/(2tau^2)*beta^c'K beta^c

  double logold = likep->loglikelihood();

  if(adaptiv)
    betaKbeta =Kenv.compute_quadform(beta,0);

  logold -= 0.5*betaKbeta*lambda;

  // Compute new weights and tildey's based on old betamode

  betamodeold.assign(betamode);
  update_linpred_diff(betamodeold,beta);

  if (  (optionsp->get_nriter() < optionsp->get_burnin())  ||
        ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) )
     )
    {
    likep->compute_IWLS_weight_tildey(weightiwls,tildey,column);
    workbeta = betamodeold.getV();
    if (varcoeff)
      compute_XWX_XWtildey_varcoeff_env(weightiwls,tildey,workbeta,0);
    else
      compute_XWX_XWtildey_env(weightiwls,tildey,workbeta,0);

    precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    }
  else
    {
    likep->tilde_y_minus_eta(tildey,column);
    workbeta = betamodeold.getV();
    compute_muy(workbeta);

    if (!lambdaconst)
      precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    }

  // Compute new betamode

  precenv.solve(muy,betamode);

  // compute -1/2 (beta^c-betamode)'P (beta^c-betamode)
  // corresponds to log(q(beta^p,beta^c))

  diff.minus(beta,betamode);
  double qoldbeta = -0.5*precenv.compute_quadform(diff,0);

  // compute new proposal

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  precenv.solveU(beta,betamode);

  // compute -1/2 (beta^p-betamode)'P (beta^p-betamode)
  // corresponds to log(q(beta^c,beta^p))

  diff.minus(beta,betamode);
  double qnewbeta = -0.5*precenv.compute_quadform(diff,0);

  // compute new predictor with proposed beta

  update_linpred_diff(beta,betamodeold);

  // Compute log-likelihood with proposed beta^p:
  // log(l)-1/(2tau^2)*beta^p'K beta^p

  double lognew = likep->loglikelihood();
  lognew  -= 0.5*Kenv.compute_quadform(beta,0)*lambda;

  // Accept/reject proposed beta

  double u = log(uniform());
  if (u <= (lognew - logold  + qoldbeta - qnewbeta) )
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


    if(!adaptiv)
      betaKbeta=Kenv.compute_quadform(beta,0);

    betaold.assign(beta);
    }
  else
    {
    update_linpred_diff(betaold,beta);
    beta.assign(betaold);
    }

  if (notransform==false)
    transform = likep->get_trmult(column);
  else
    transform = 1;  


  FULLCOND::update();

  }


void FULLCOND_nonp_gaussian::update_IWLS(void)
  {

//  int j;
  unsigned i;

  if (optionsp->get_nriter() == 1)
    {
    betaold.assign(beta);
    betaKbeta=Kenv.compute_quadform(beta,0);
    }

  double * workbeta;

  if (!lambdaconst)
    lambda = 1.0/sigma2;

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

  if(adaptiv)
    betaKbeta=Kenv.compute_quadform(beta,0);

  logold -= 0.5*betaKbeta*lambda;

  workbeta = betaold.getV();
  if (  (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) )
      )
    {
    if (varcoeff)
      compute_XWX_XWtildey_varcoeff_env(weightiwls,tildey,workbeta,0);
    else
      compute_XWX_XWtildey_env(weightiwls,tildey,workbeta,0);

    precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    }
  else
    {
    compute_muy(workbeta);

    if (!lambdaconst)
      precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    }


  precenv.solve(muy,betahelp);

  double * work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = rand_normal();

  precenv.solveU(beta,betahelp);

  betahelp.minus(beta,betahelp);

  double qold = 0.5*precenv.getLogDet()-
                0.5*precenv.compute_quadform(betahelp,0);

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
  lognew  -= 0.5*Kenv.compute_quadform(beta,0)*lambda;

  workbeta = beta.getV();
  if (  (optionsp->get_nriter() < optionsp->get_burnin()) ||
      ( (updateW != 0) && ((optionsp->get_nriter()-1) % updateW == 0) )
      )
    {
    if (varcoeff)
      compute_XWX_XWtildey_varcoeff_env(weightiwls,tildey,workbeta,0);
    else
      compute_XWX_XWtildey_env(weightiwls,tildey,workbeta,0);

    precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    }
  else
    {
    compute_muy(workbeta);

    if (!lambdaconst)
      precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    }


  precenv.solve(muy,betahelp);

  betahelp.minus(betaold,betahelp);
  double qnew = 0.5*precenv.getLogDet() -
                0.5*precenv.compute_quadform(betahelp,0);


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

  if (notransform==false)
    transform = likep->get_trmult(column);
  else
    transform = 1;

  FULLCOND::update();

  }


void FULLCOND_nonp_gaussian::update_lambdaconst(void)
  {

  int j;
  unsigned i;

  int * workindex;

  update_linpred(false);

  if (optionsp->get_nriter()==1)
    {
    if (varcoeff)
      compute_XWX_varcoeff_env(likep->get_weight());
    else
      compute_XWX_env(likep->get_weight());

    precenv.addtodiag(XXenv,Kenv,1.0,lambda);
    }


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

  if (notransform==false)
    transform = likep->get_trmult(column);
  else
    transform = 1;


  FULLCOND::update();

  }


void FULLCOND_nonp_gaussian::update(void)
  {

  if (utype==iwls)
    {
    update_IWLS();
    }
  else if (utype==iwlsmode)
    {
    update_IWLS_mode();
    }
  else if (utype==hyperblock)
    {
    if (lambdaconst == true && changingweight == false)
      update_IWLS();
    else
      update_IWLS_hyperblock();
    }
  else if (utype==hyperblockmode)
    {
    if (lambdaconst == true && changingweight == false)
      update_IWLS_mode();
    else
      update_IWLS_hyperblock_mode();
    }
  else if (utype==gaussianlaplace)
    {
//    update_gaussian_laplace();
    update_gaussian_gemanreynolds();
    }
  else
  {
  if (lambdaconst == true && changingweight == false)
    {
    update_lambdaconst();
    }
  else
    {
    int j;
    unsigned i;

    int * workindex;

    update_linpred(false);

    lambda = likep->get_scale(column)/sigma2;

    if (optionsp->get_nriter()==1 || changingweight)
      {
      if (varcoeff)
        compute_XWX_varcoeff_env(likep->get_weight());
      else
        compute_XWX_env(likep->get_weight());
      }

    precenv.addtodiag(XXenv,Kenv,1.0,lambda);

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

// Multiplikative Effekte: Zentrieren des Effekts
    if (notransform)
      double m = centerbeta();

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

    if (notransform==false)
      transform = likep->get_trmult(column);
    else
      transform = 1;


    FULLCOND::update();

    }

  }

  }


void FULLCOND_nonp_gaussian::update_gaussian_laplace(void)
  {
  unsigned i,j;

  double sqrtscale = 0.3;
  double help;
  double sum;
  double w;

  unsigned beg,end;

  for(i=0;i<nrpar;i++)
    {

    beg=posbeg[i];
    end=posend[i];

    betaold.assign(beta);

    double logold;
    logold = likep->loglikelihood(beg,end,index);

    sum = 0.0;
    if(delta.rows()>1)
      for(j=0;j<neighbors[i].size();j++)
        {
        w = -Kenv(i,neighbors[i][j]);
        sum += w * fabs(betaold(i,0)-betaold(neighbors[i][j],0));
        }
    else
      for(j=0;j<neighbors[i].size();j++)
        sum += fabs(betaold(i,0)-betaold(neighbors[i][j],0));

//      logold -= Kenv.compute_sumfabsdiff(betaold,0)/sigma2;
    logold -= sum/sigma2;

    help = sqrtscale*rand_normal();
    beta(i,0) = betaold(i,0) + help;

//    update_linpred_diff(beta,betaold);
    update_linpred_diff(beg,end,help);

    double lognew;
    lognew = likep->loglikelihood(beg,end,index);

    sum = 0.0;
    if(delta.rows()>1)
      for(j=0;j<neighbors[i].size();j++)
        {
        w = -Kenv(i,neighbors[i][j]);
        sum += w * fabs(beta(i,0)-beta(neighbors[i][j],0));
        }
    else
      for(j=0;j<neighbors[i].size();j++)
        sum += fabs(beta(i,0)-beta(neighbors[i][j],0));

//    lognew -= Kenv.compute_sumfabsdiff(beta,0)/sigma2;
    lognew -= sum/sigma2;

    double alpha = lognew - logold;
    double u = log(uniform());

    nrtrials++;

    if ( u <= alpha )
      {
      acceptance++;
      }
    else
      {
//      update_linpred_diff(betaold,beta);
      update_linpred_diff(beg,end,-help);
      beta.assign(betaold);
      }

    }


  if (center)
    {
    double m = centerbeta();
    if (varcoeff)
      fcconst->update_fix_varcoeff(m,datanames[1]);
    else
      fcconst->update_intercept(m);
    }


  if (notransform==false)
    transform = likep->get_trmult(column);
  else
    transform = 1;


  FULLCOND::update();


/*
  int j;
  unsigned i;

  int * workindex;

  if (optionsp->get_nriter() == 1)
    betaold.assign(beta);

  double logold;
  logold = likep->loglikelihood();
  logold -= Kenv.compute_sumfabsdiff(betaold,0)/sigma2;

  update_linpred(false);

  lambda = likep->get_scale(column)/sigma2;
  double scale = likep->get_scale(column);

  if (optionsp->get_nriter()==1 || changingweight)
    {
    if (varcoeff)
      compute_XWX_varcoeff_env(likep->get_weight());
    else
      compute_XWX_env(likep->get_weight());
    }

  precenv.addtodiag(XXenv,Kenv,1.0,lambda);

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

  datamatrix mhelp = betahelp;
  precenv.solve(muy,mhelp);

  betahelp.minus(beta,mhelp);
  double qold = - 0.5*precenv.compute_quadform(betahelp,0)/scale;

  double lognew;
  lognew = likep->loglikelihood();
  lognew -= Kenv.compute_sumfabsdiff(beta,0)/sigma2;

  betahelp.minus(betaold,mhelp);
  double qnew = - 0.5*precenv.compute_quadform(betahelp,0)/scale;

  double u = log(uniform());
  double alpha = (lognew - logold  + qnew - qold);
  if ( u <= alpha )
    {
    acceptance++;

    if (center)
      {
      double m = centerbeta();
      fcconst->update_intercept(m);
      }

    betaold.assign(beta);
    }
  else
    {
    update_linpred_diff(betaold,beta);
    beta.assign(betaold);
    }


  if (notransform==false)
    transform = likep->get_trmult(column);
  else
    transform = 1;


  FULLCOND::update();
*/
  }

void FULLCOND_nonp_gaussian::update_gaussian_gemanreynolds(void)
  {
  unsigned i,j;

  double sqrtscale = 0.3;
  double help;
  double sum;
  double w;
  double u;

  double p = 2.0;
  double lambda = 3.0;


  unsigned beg,end;

  for(i=0;i<nrpar;i++)
    {

    beg=posbeg[i];
    end=posend[i];

    betaold.assign(beta);

    double logold;
    logold = likep->loglikelihood(beg,end,index);

    sum = 0.0;
    if(delta.rows()>1)
      for(j=0;j<neighbors[i].size();j++) // adaptive Gewichte!!!
        {
        w = -Kenv(i,neighbors[i][j]);
        sum += w * fabs(betaold(i,0)-betaold(neighbors[i][j],0));
        }
    else
      for(j=0;j<neighbors[i].size();j++)
        {
//        sum += fabs(betaold(i,0)-betaold(neighbors[i][j],0));
        u = (betaold(i,0)-betaold(neighbors[i][j],0))/sigma2;
        sum += lambda/(1.0 + pow(fabs(u),p));
        }

//    logold -= sum/sigma2;
    logold -= sum;

    help = sqrtscale*rand_normal();
    beta(i,0) = betaold(i,0) + help;

    update_linpred_diff(beg,end,help);

    double lognew;
    lognew = likep->loglikelihood(beg,end,index);

    sum = 0.0;
    if(delta.rows()>1)
      for(j=0;j<neighbors[i].size();j++) // adaptive Gewichte!!!
        {
        w = -Kenv(i,neighbors[i][j]);
        sum += w * fabs(beta(i,0)-beta(neighbors[i][j],0));
        }
    else
      for(j=0;j<neighbors[i].size();j++)
        {
//        sum += fabs(beta(i,0)-beta(neighbors[i][j],0));
        u = (beta(i,0)-beta(neighbors[i][j],0))/sigma2;
        sum += lambda/(1.0 + pow(fabs(u),p));
        }

//    lognew -= sum/sigma2;
    lognew -= sum/sigma2;

    double alpha = lognew - logold;
    double u = log(uniform());

    nrtrials++;

    if ( u <= alpha )
      {
      acceptance++;
      }
    else
      {
      update_linpred_diff(beg,end,-help);
      beta.assign(betaold);
      }

    }


  if (center)
    {
    double m = centerbeta();
    if (varcoeff)
      fcconst->update_fix_varcoeff(m,datanames[1]);
    else
      fcconst->update_intercept(m);
    }


  if (notransform==false)
    transform = likep->get_trmult(column);
  else
    transform = 1;


  FULLCOND::update();

  }



bool FULLCOND_nonp_gaussian::posteriormode(void)
  {
  int j;
  unsigned i;

  int * workindex;

  update_linpred(false);

  if ( (lambda_prec != lambda) ||
     (likep->iwlsweights_constant() == false) || (changingweight==true) )
    {

    if ( (likep->iwlsweights_constant() == false ) || (changingweight==true) )
      {
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

// Multiplikative Effekte: Zentrieren des Effekts
  if (notransform)
    double m = centerbeta();

  update_linpred(true);

  if (center)
    {
    double m = centerbeta();
    if (varcoeff)
      fcconst->posteriormode_fix_varcoeff(m,datanames[1]);
    else
    fcconst->posteriormode_intercept(m);
    }

  if (notransform==false)
    transform = likep->get_trmult(column);
  else
    transform = 1;


  return FULLCOND_nonp_basis::posteriormode();

  }


bool FULLCOND_nonp_gaussian::posteriormode_converged(const unsigned & itnr)
  {
  return likep->posteriormode_converged_fc(beta,beta_mode,itnr);
  }


void FULLCOND_nonp_gaussian::get_effectmatrix(datamatrix & e,
                                        vector<ST::string> & enames,unsigned be,
                                        unsigned en, effecttype t)
  {

  int * workindex = index.getV();

  double * workbeta;
  if (t==MCMC::current || t==MCMC::fvar_current)
    workbeta = beta.getV();
  else if (t==MCMC::mean || t==MCMC::fvar_mean)
    workbeta = betamean.getV();
  else
    workbeta = betaqu50.getV();

  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();

  int j;
  unsigned i,k;

  if (varcoeff)
    {

    if (t==MCMC::fvar_current || t==MCMC::fvar_mean  || t==MCMC::fvar_median)
      {

      for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
        {
        if (*itbeg != -1)
          {
          for(j=*itbeg;j<=*itend;j++,workindex++)
            e(*workindex,be) = *workbeta;
          }
        }

      }
    else
      {

      double * workdata=data.getV();

      vector<ST::string>::iterator effit = effectvalues.begin();
//      int t;
      enames.push_back("f_"+datanames[0]+"_"+datanames[1]);
      enames.push_back(datanames[0]);
      enames.push_back(datanames[1]);

      for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend,++effit)
        {
        if (*itbeg != -1)
          {
          for(j=*itbeg;j<=*itend;j++,workindex++,workdata++)
            {
            e(*workindex,be) = *workbeta*(*workdata);
            t = (MCMC::effecttype)((*effit).strtodouble(e(*workindex,be+1)));
            e(*workindex,be+2) = *workdata;
            }
          }
        }

      }
    }
  else
    {
    vector<ST::string>::iterator effit = effectvalues.begin();
//    int t;

    enames.push_back("f_"+datanames[0]);
    enames.push_back(datanames[0]);


    for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend,++effit)
      {
      if (*itbeg != -1)
        {
        for (k=(*itbeg);k<=(*itend);k++,workindex++)
          {
          e(*workindex,be) = *workbeta;
          t = (MCMC::effecttype)((*effit).strtodouble(e(*workindex,be+1)));
          }
        }
      }
    }

  }


unsigned FULLCOND_nonp_gaussian::get_nreffects(effecttype t)
  {
  if (varcoeff)
    {
    if (t==MCMC::fvar_current || t==MCMC::fvar_mean  || t==MCMC::fvar_median)
      return 1;
    else
      return 3;
    }
  else
    return 2;
  }


void FULLCOND_nonp_gaussian::outresults(void)
  {
  if (lattice==false)
    FULLCOND_nonp_basis::outresults();
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
    outres << "pqu"  << l1  << "   ";
    outres << "pqu"  << l2  << "   ";
    outres << "pmed   ";
    outres << "pqu"  << u1  << "   ";
    outres << "pqu"  << u2  << "   ";
    outres << "pcat" << level1 << "   ";
    outres << "pcat" << level2 << "   ";

    outres << endl;

      double * workmean = betamean.getV();
      double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
      double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
      double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
      double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
      double * workbetaqu50 = betaqu50.getV();
      double * workxyvalues = xyvalues.getV();

      vector<ST::string>::iterator effit = effectvalues.begin();

      for(i=0;i<nrpar;i++,++effit,workmean++,workbetaqu_l1_lower_p++,
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


void FULLCOND_nonp_gaussian::outoptions(void)
  {
  FULLCOND_nonp_basis::outoptions();

  if (utype!=gaussian)
    {
    if (utype == iwls)
      optionsp->out("  Proposal: IWLS based on current regression coefficients\n");
    else if (utype == iwlsmode)
      optionsp->out("  Proposal: IWLS based on current mode\n");
    else if (utype == hyperblock)
      {
      optionsp->out("  Proposal: IWLS based on current regression coefficients\n");
      optionsp->out("            Variance parameter is updated jointly with regression coefficients\n");
      }
    else if (utype == hyperblockmode)
      {
      optionsp->out("  Proposal: IWLS based on current mode\n");
      optionsp->out("            Variance parameter is updated jointly with regression coefficients\n");
      }
    }

  }





} // end: namespace MCMC


