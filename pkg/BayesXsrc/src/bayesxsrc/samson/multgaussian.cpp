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



#include "multgaussian.h"

namespace MCMC
{


void DISTRIBUTION_multgaussian::standardise(void)
  {
  unsigned i,j;
  double s;
  trmult = datamatrix(nrcat,1);

  for (j=0;j<nrcat;j++)
    {
    s = sqrt(response.var(j,weight));
    trmult(j,0) = s;
    }

  datamatrix tr(nrcat,nrcat);
  for (j=0;j<nrcat;j++)
    for(i=0;i<nrcat;i++)
      tr(j,i) = trmult(j,0)*trmult(i,0);

  Scalesave.set_transformmult(tr);

  double * workresp = response.getV();
  double * worklin = (*linpred_current).getV();
  for (i=0;i<nrobs;i++)
   {
   for (j=0;j<nrcat;j++,workresp++,worklin++)
     {
     *workresp = *workresp/trmult(j,0);
     *worklin = *worklin/trmult(j,0);
     }
   }

  }


DISTRIBUTION_multgaussian::DISTRIBUTION_multgaussian(const double & a,
                   const datamatrix & b, MCMCoptions * o, const datamatrix & r,
                   const ST::string & fp,const ST::string & fs,
                   const datamatrix & w)
  : DISTRIBUTION(o,r,w,fp,fs)
  {

  nrcat = response.cols();

  Scalesave = FULLCOND(o,datamatrix(1,1),"Sigmasave",nrcat,nrcat,fs);
  Scalesave.setflags(MCMC::norelchange | MCMC::nooutput);

  standardise();

  sumB = datamatrix(r.cols(),r.cols(),0);
  diff = datamatrix(r.rows(),r.cols(),0);

  A = a;
  B = b;
  family = "Multivariate Gaussian";
  scaleexisting = true;

  scale = datamatrix(nrcat,nrcat,0);
  scale_mode = datamatrix(nrcat,nrcat,0);
  sigma_rmr = datamatrix(nrcat,1);
  unsigned i;
  for (i=0;i<nrcat;i++)
    {
    scale(i,i) = response.var(i);
    sigma_rmr(i,0) = scale(i,i);
    }

  offset = datamatrix(r.rows(),r.cols(),0);

  SIGMA_rmr = datamatrix(1,nrcat-1);
  SIGMA_mr = datamatrix(nrcat-1,nrcat-1);

  }


DISTRIBUTION_multgaussian::DISTRIBUTION_multgaussian(
const DISTRIBUTION_multgaussian & nd)
: DISTRIBUTION(DISTRIBUTION(nd))
  {
  sumB = nd.sumB;
  diff = nd.diff;
  offset = nd.offset;
  sigma_rmr = nd.sigma_rmr;
  SIGMA_mr = nd.SIGMA_mr;
  SIGMA_rmr = nd.SIGMA_rmr;
  A = nd.A;
  B = nd.B;
  nrcat = nd.nrcat;
  }


const DISTRIBUTION_multgaussian & DISTRIBUTION_multgaussian::operator=(
const DISTRIBUTION_multgaussian & nd)
  {

  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));
  diff = nd.diff;
  sumB = nd.sumB;
  offset = nd.offset;
  sigma_rmr = nd.sigma_rmr;
  SIGMA_mr = nd.SIGMA_mr;
  SIGMA_rmr = nd.SIGMA_rmr;
  A = nd.A;
  B = nd.B;
  nrcat = nd.nrcat;
  return *this;
  }


void DISTRIBUTION_multgaussian::compute_mu(const double * linpred,double * mu)
                                           const
  {
  unsigned i;
  double * worktrmult= trmult.getV();
  for(i=0;i<linearpred.cols();i++,linpred++,mu++,worktrmult++)
    *mu = *worktrmult * *linpred;
  }

void DISTRIBUTION_multgaussian::compute_deviance(const double * response,
                             const double * weight,const double * mu,
                             double * deviance,double * deviancesat,
                             const datamatrix & scale,const int & i) const
  {


  if (*weight != 0)
    {
    unsigned i,j;
    datamatrix scaleinv(nrcat,nrcat);
    double * workscale = scale.getV();
    double * workscaleinv = scaleinv.getV();
    datamatrix responsetr(nrcat,1);
    for (i=0;i<nrcat;i++,response++,mu++)
      {
      responsetr(i,0)= *response * trmult(i,0) - *mu;
      for(j=0;j<nrcat;j++,workscale++,workscaleinv++)
        *workscaleinv = trmult(i,0)*trmult(j,0) * *workscale;
      }

    double de=scaleinv.det();

    scaleinv=scaleinv.inverse();


    *deviancesat = (responsetr.transposed()*scaleinv*responsetr)(0,0);
    *deviance =  *deviancesat + log(de);

    }
  else
    {
    *deviance = 0;
    *deviancesat = 0;
    }



  }


void DISTRIBUTION_multgaussian::compute_SIGMA_rmr(unsigned r)
  {

  unsigned k,l;

  l=0;
  for (k=0;k<nrcat;k++)
    {
    if (k != r)
      {
      SIGMA_rmr(0,l) = scale(r,k);
      l++;
      } //end: if (k != i)
    } //end: for (k=0;k<nrcat;k++)

  }


void DISTRIBUTION_multgaussian::compute_SIGMA_mr(unsigned r)
  {

  unsigned n,l,k,m;

  n = 0;
  for (k=0;k<nrcat;k++)
    {
    if (k != r)
      {
      l = 0;
      for (m=0;m<nrcat;m++)
        {
        if (m != r)
          {
          SIGMA_mr(n,l) = scale(k,m);
          l++;
          }
        } // end: for (m=0;m<nrcat;m++)
        n++;
      } // end: if (k != i)

    } // end: for (k=0;k<nrcat;k++)

  SIGMA_mr = SIGMA_mr.inverse();

  }


void DISTRIBUTION_multgaussian::compute_offset(void)
  {

  if (nrcat==2)
    {

    unsigned i;
    double f1 = scale(0,1)/scale(1,1);
    double f2 = scale(0,1)/scale(0,0);
    for (i=0;i<nrobs;i++)
      {
      offset(i,0) = f1*(response(i,1)-(*linpred_current)(i,1));
      offset(i,1) = f2*(response(i,0)-(*linpred_current)(i,0));
      }

    }
  else
    {

    unsigned i,k,l,r;

    datamatrix help(1,nrcat-1);

    for (r=0;r<nrcat;r++)
      {

      compute_SIGMA_rmr(r);

      compute_SIGMA_mr(r);

      help.mult(SIGMA_rmr,SIGMA_mr);

      double di;

      for (i=0;i<nrobs;i++)
        {

        offset(i,r) = 0;

        l = 0;
        for (k=0;k<nrcat;k++)
          {
          if (k != r)
            {
            di = response(i,k)-(*linpred_current)(i,k);
            offset(i,r) += help(0,l)*di;
            l++;
            } // end: if (k != r)
          } //end: for (k=0;k<nrcat;k++)
        } // end: for (i=0;i<nrobs;i++)
      } // end: for (r=0;r<offset.cols();r++)
    }

  } // end: compute_offset


void DISTRIBUTION_multgaussian::compute_IWproduct(void)
  {

  diff.minus(response,*linpred_current);

  unsigned i,j,k;

  double * sumBwork = sumB.getV();
  double * diff1;
  double * diff2;

  for(i=0;i<sumB.rows();i++)
    {
    for(j=0;j<sumB.cols();j++,sumBwork++)
      {
      *sumBwork = 0;
      diff1 = diff.getV()+i;
      diff2 = diff.getV()+j;

      for(k=0;k<nrobs;k++,diff1+=sumB.cols(),diff2+=sumB.cols())
        *sumBwork +=  *diff1 * *diff2;

      *sumBwork *= 0.5;
      }
    }

  }


void DISTRIBUTION_multgaussian::compute_respminuslinpred(datamatrix & res,
const unsigned & co)
  {

  compute_offset();

  unsigned i;
  double * reswork = res.getV();
  double * responsework = response.getV() + co;
  double * linwork = (*linpred_current).getV()+ co;
  double * offsetwork = offset.getV()+co;
  unsigned c = response.cols();
  for (i=0;i<nrobs;i++,reswork++,responsework+=c,linwork+=c,offsetwork+=c)
    *reswork = *responsework - *linwork - *offsetwork;
  }


void DISTRIBUTION_multgaussian::compute_sigmarmr(void)
  {

  if (nrcat==2)
    {
    sigma_rmr(0,0) = scale(0,0)-(scale(0,1)*scale(0,1))/scale(1,1);
    sigma_rmr(1,0) = scale(1,1)-(scale(0,1)*scale(0,1))/scale(0,0);
    }
  else
    {

    unsigned r;
    datamatrix help;

    for (r=0;r<nrcat;r++)
      {
      compute_SIGMA_mr(r);
      compute_SIGMA_rmr(r);
      help = SIGMA_rmr*SIGMA_mr*SIGMA_rmr.transposed();
      sigma_rmr(r,0) = scale(r,r) - help(0,0);
      } // end: for (r=0;r<nrcat;r++)

    }

  }


bool DISTRIBUTION_multgaussian::posteriormode(void)
  {

  register unsigned i,j;

  double help;

  double * worklin = (*linpred_current).getV();
  double * workresp = response.getV();
  double * workweight = weight.getV();

  datamatrix sum(nrcat,1,0);

  for (i=0;i<nrobs;i++,workweight++)
    {
    for(j=0;j<nrcat;j++,worklin++,workresp++)
       {
       help = *workresp - *worklin;
       sum(j,0) += *workweight*help*help;
       }
    }


  for(j=0;j<nrcat;j++)
    {
    scale(j,j) = (1.0/nrobs)*sum(j,0);
    sigma_rmr(j,0) = scale(j,j);
    }

  return true;

  }


void DISTRIBUTION_multgaussian::update(void)
  {

  compute_IWproduct();

  sumB.plus(B,sumB);

  sumB = 0.5*sumB.inverse();

  randnumbers::rand_wishart(sumB,2.0*A+nrobs,scale);

  scale = scale.inverse();

  compute_sigmarmr();

  DISTRIBUTION::update();

  }


} // end: namespace MCMC

