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




#include "distr_categorical_mult.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_multgaussian ----------------------
//------------------------------------------------------------------------------


DISTR_multgaussian::DISTR_multgaussian(const double & a,const double & b,
                                       unsigned & cnr,
                                       GENERAL_OPTIONS * o,
                                       bool mast,
                                       const datamatrix & r,
                                       const ST::string & ps,
                                       const datamatrix & w)
  : DISTR_multinomprobit(o,cnr,mast,r,w)

  {
  family = "Multivariate Gaussian";
  pathscale = ps;
  predict_mult=true;
  A = a;
  double h = sqrt(response.var(0,weight));
  helpquantity1  = h*b;
  trmult = h;

  wtype = wweightsnochange_constant;

  }




const DISTR_multgaussian & DISTR_multgaussian::operator=(
                                      const DISTR_multgaussian & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_multinomprobit::operator=(DISTR_multinomprobit(nd));
  devianceconst = nd.devianceconst;
  FC_scale = nd.FC_scale;
  pathscale = nd.pathscale;
  FC_corr = nd.FC_corr;
  sumB = nd.sumB;
  diff = nd.diff;
  offset = nd.offset;
  SIGMA_mr = nd.SIGMA_mr;
  SIGMA_rmr = nd.SIGMA_rmr;
  A = nd.A;
  B = nd.B;
  return *this;
  }


DISTR_multgaussian::DISTR_multgaussian(const DISTR_multgaussian & nd)
   : DISTR_multinomprobit(DISTR_multinomprobit(nd))
  {
  devianceconst = nd.devianceconst;
  FC_scale = nd.FC_scale;
  pathscale = nd.pathscale;
  FC_corr = nd.FC_corr;
  sumB = nd.sumB;
  diff = nd.diff;
  offset = nd.offset;
  SIGMA_mr = nd.SIGMA_mr;
  SIGMA_rmr = nd.SIGMA_rmr;
  A = nd.A;
  B = nd.B;
  }


void DISTR_multgaussian::outoptions(void)
  {
  DISTR::outoptions();

  optionsp->out("  Hyperparameter a: " + ST::doubletostring(A,6) + "\n");
  optionsp->out("  Hyperparameter b: " + ST::doubletostring(helpquantity1,6) + "\n");

  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_multgaussian::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {
  return 0;
  }



double DISTR_multgaussian::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {
  return 0;
  }


void DISTR_multgaussian::compute_mu(const double * linpred,double * mu)
  {
  *mu = *linpred;
  }

void DISTR_multgaussian::compute_mu_mult(vector<double *> linpred,double * mu)
  {
  *mu = *linpred[catnr];
  }

datamatrix * DISTR_multgaussian::get_auxiliary_parameter(auxiliarytype t)
  {

   return &helpmat2;

  /*
  if (master)
    {

    if (t == current)
      return &FC_scale.beta;
    else
      return &FC_scale.betamean;
    }
  else
    return &helpmat1;
  */  }


void DISTR_multgaussian::compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix *> aux)
  {

  if (*weight[nrcat-1] != 0)
    {

    unsigned j,k;
    double qf=0;
    double diffj;
    for (j=0;j<nrcat;j++)
      {
      diffj = (*response[j])-(*linpred[j]);
      qf += (*aux[nrcat-1])(j,j)* pow(diffj,2);
      if (j<nrcat-1)
        {
        for(k=j+1;k<nrcat;k++)
          {
          qf += 2*(*aux[nrcat-1])(j,k)*diffj*((*response[k])-(*linpred[k]));
          }
        }
      }

    double l = -devianceconst+0.5*log((*aux[nrcat-1]).det())-0.5*qf;

    *deviance = -2*l;

    }
  else
    *deviance = 0;


  }




double DISTR_multgaussian::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  *workingweight=*weight;
  *workingresponse = *response;
  if (like && (*weight != 0))
    return  - *weight * (pow(*response-(*linpred),2))/(2* sigma2);
  else
    return 0;

  }



void DISTR_multgaussian::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  *workingweight=1;
  *workingresponse = *response;
  if (compute_like)
    like -=  (pow(*response-(*linpred),2))/(2* sigma2);


  }


void DISTR_multgaussian::compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  *workingresponse = *response;
  if (compute_like && *workingweight!=0)
    like  -= *workingweight * (pow(*response-(*linpred),2))/(2* sigma2);
  }

void DISTR_multgaussian::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  *workingresponse = *response;
  if (compute_like)
    like -=  (pow(*response-(*linpred),2))/(2* sigma2);
  }


void DISTR_multgaussian::compute_sigmarmr(void)
  {

  unsigned r;
  datamatrix help;

//  ofstream out("d:\\_sicher\\anett\\scale.raw");
//  FC_scale.beta.prettyPrint(out);

  for (r=0;r<nrcat;r++)
    {
    compute_SIGMA_mr(r);
    compute_SIGMA_rmr(r);
    othercat[r]->helpmat1 = SIGMA_rmr * SIGMA_mr;
    help = othercat[r]->helpmat1 * SIGMA_rmr.transposed();
    othercat[r]->sigma2 = FC_scale.beta(r,r) - help(0,0);
    } // end: for (r=0;r<nrcat;r++)

  }


void DISTR_multgaussian::compute_SIGMA_rmr(unsigned r)
  {

  unsigned k,l;

  l=0;
  for (k=0;k<nrcat;k++)
    {
    if (k != r)
      {
      SIGMA_rmr(0,l) = FC_scale.beta(r,k);
      l++;
      } //end: if (k != i)
    } //end: for (k=0;k<nrcat;k++)

  }


void DISTR_multgaussian::compute_SIGMA_mr(unsigned r)
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
          SIGMA_mr(n,l) = FC_scale.beta(k,m);
          l++;
          }
        } // end: for (m=0;m<nrcat;m++)
        n++;
      } // end: if (k != i)

    } // end: for (k=0;k<nrcat;k++)

  SIGMA_mr = SIGMA_mr.inverse();

  }


void DISTR_multgaussian::compute_offset(void)
  {
  unsigned i,j,k;

  vector<double *> workresp;
  vector<double *> worklin;

  double * workresp_c;
  double * workrespcat = response.getV();

  for (j=0;j<nrcat;j++)
    {
    if (j != catnr)
      {

      if (othercat[j]->linpred_current==1)
        worklin.push_back(othercat[j]->linearpred1.getV());
      else
        worklin.push_back(othercat[j]->linearpred2.getV());

      workresp.push_back(othercat[j]->response.getV());

      }
    else
      {
      workresp_c = workingresponse.getV();
      }

    }

  double o;

  for (i=0;i<nrobs;i++,workresp_c++,workrespcat++)
    {

    o = 0;
    for (k=0;k<nrcat-1;k++)
      {
      o += helpmat1(0,k)*((*workresp[k]) - (*worklin[k]));
      if (i<nrobs-1)
        {
        workresp[k]++;
        worklin[k]++;
        }
      }

    *workresp_c = *workrespcat-o;

    } // end: for (i=0;i<nrobs;i++)

  } // end: compute_offset


void DISTR_multgaussian::initpointer(unsigned j,double* & worklin, double* & workresp)
  {

  if (othercat[j]->linpred_current==1)
    worklin = othercat[j]->linearpred1.getV();
  else
    worklin = othercat[j]->linearpred2.getV();

  workresp = othercat[j]->response.getV();

  }


void DISTR_multgaussian::compute_IWproduct(void)
  {

  unsigned i,j,k;

  double * sumBwork = sumB.getV();

  double * worklin1;
  double * worklin2;

  double * workresp1;
  double * workresp2;


  for(i=0;i<sumB.rows();i++)
    {
    for(j=0;j<sumB.cols();j++,sumBwork++)
      {
      *sumBwork = 0;

      initpointer(i,worklin1,workresp1);
      initpointer(j,worklin2,workresp2);

      for(k=0;k<nrobs;k++,worklin1++,worklin2++,workresp1++,workresp2++)
        *sumBwork +=  (*workresp1-(*worklin1)) * (*workresp2-(*worklin2));

      *sumBwork *= 0.5;
      }
    }

  }


void DISTR_multgaussian::update(void)
  {
  if (master)
    {
    compute_IWproduct();

    sumB.plus(B,sumB);

    sumB = 0.5*sumB.inverse();

    randnumbers::rand_wishart(sumB,2.0*A+nrobs,helpmat2);

    FC_scale.beta = helpmat2.inverse();

    compute_sigmarmr();

    unsigned i,j;
    for(i=0;i<nrcat;i++)
      {
      for(j=0;j<nrcat;j++)
        FC_corr.beta(i,j) = FC_scale.beta(i,j)/sqrt(FC_scale.beta(i,i)*FC_scale.beta(j,j));
      }

    FC_scale.update();

    FC_corr.update();

    DISTR::update();
    }

  compute_offset();

  }


bool DISTR_multgaussian::posteriormode(void)
  {

  if (offset.rows() ==  1)
    {

    offset = datamatrix(nrobs,1,0);

    helpmat1 = datamatrix(1,nrcat-1,0);

    devianceconst = double(nrcat)/2*log(6.2831853);

    if (master)
      {

      FC_scale = FC(optionsp,"",nrcat,nrcat,pathscale);
      FC_corr = FC(optionsp,"",nrcat,nrcat,"");

      sumB = datamatrix(nrcat,nrcat,0);
      diff = datamatrix(nrcat,nrcat,0);

      SIGMA_rmr = datamatrix(1,nrcat-1);
      SIGMA_mr = datamatrix(nrcat-1,nrcat-1);


      B = datamatrix(nrcat,nrcat,0);
      unsigned j;
      for (j=0;j<B.rows();j++)
        B(j,j) = othercat[j]->helpquantity1;

      }

    }


  if (master)
    {
    register unsigned i,j;

    double * worklin;
    double * workresp;
    double * workweight;

    double sum;

    for(j=0;j<nrcat;j++)
      {
      sum = 0;

      workweight = weight.getV();

      initpointer(j,worklin,workresp);

      for (i=0;i<nrobs;i++,workweight++,workresp++,worklin++)
        sum += *workweight*pow(*workresp - *worklin,2);

      FC_scale.beta(j,j) = (1.0/nrobs)*sum;

      othercat[j]->sigma2 = FC_scale.beta(j,j);

      }

    helpmat2 = FC_scale.beta.inverse();

    FC_scale.posteriormode_betamean();
    }

  return true;

  }


void DISTR_multgaussian::outresults_help(ST::string t,datamatrix & r)
  {

  ST::string help;

  optionsp->out("  " + t + ":\n");
  optionsp->out("\n");

  unsigned i,j;
  for (i=0;i<nrcat;i++)
    {
    help="  ";
    for(j=0;j<nrcat;j++)
      {
      help = help + ST::doubletostring(r(i,j),6) + "   ";
      }
    optionsp->out(help + "\n");
    }

  optionsp->out("\n");
  }


void DISTR_multgaussian::outresults(ofstream & out_stata,ofstream & out_R,
                                   ST::string pathresults)
  {

  if (master)
    {

    FC_scale.outresults(out_stata,out_R,"");

    ST::string l1 = ST::doubletostring(optionsp->lower1,4);
    ST::string l2 = ST::doubletostring(optionsp->lower2,4);
    ST::string u1 = ST::doubletostring(optionsp->upper1,4);
    ST::string u2 = ST::doubletostring(optionsp->upper2,4);
    l1 = l1.replaceallsigns('.','p');
    l2 = l2.replaceallsigns('.','p');
    u1 = u1.replaceallsigns('.','p');
    u2 = u2.replaceallsigns('.','p');

    optionsp->out("  Estimation results for the covariance matrix\n",true);
    optionsp->out("\n");

    outresults_help("Posterior mean",FC_scale.betamean);
//    outresults_help("Posterior variance",FC_scale.betavar);
//    outresults_help("Posterior " + l1 + " percent quantile",
//                  FC_scale.betaqu_l1_lower);
//    outresults_help("Posterior median",FC_scale.betaqu50);
//    outresults_help("Posterior " + u2 + " percent quantile",
//                    FC_scale.betaqu_l1_upper);

    optionsp->out("\n");
    optionsp->out("\n");

    optionsp->out("  Estimation results for correlation matrix\n",true);
    optionsp->out("\n");

    outresults_help("Posterior mean",FC_corr.betamean);
  //  outresults_help("Posterior variance",FC_corr.betavar);

    }

  }


void DISTR_multgaussian::get_samples(const ST::string & filename,ofstream & outg) const
  {
  if (master && (filename.isvalidfile() != 1))
    FC_scale.get_samples(filename,outg);
  }


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_multinomprobit --------------------
//------------------------------------------------------------------------------


void DISTR_multinomprobit::assign_distributions(vector<DISTR * > dp)
  {
  othercat = dp;
  nrcat = dp.size();
  }

void DISTR_multinomprobit::assign_othercat(DISTR* o)
  {
  othercat.push_back(o);
  nrcat++;
  nrothercat = othercat.size();
  }



void DISTR_multinomprobit::create_responsecat(void)
  {

  responsecat = datamatrix(nrobs,1,nrothercat);

  unsigned i,j;

  bool found;
  for (i=0;i<nrobs;i++)
    {
    found = false;
    j=0;
    while ((found==false) && (j < nrothercat) )
      {
      if (othercat[j]->response(i,0) == 1)
        {
        responsecat(i,0) = j;
        found=true;
        }

      j++;
      }

    if ((found==false) && response(i,0)==1)    // master
      {
//      responsecat(i,0) = nrothercat;
      found=true;
      }

    if (found==false)   // reference
      responsecat(i,0) = -1;

    }

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\responsecat.raw");
  // responsecat.prettyPrint(out);
  // ENDE: TEST

  }


DISTR_multinomprobit::DISTR_multinomprobit(GENERAL_OPTIONS * o,
                                           unsigned & cnr,
                                           bool mast,
                                           const datamatrix & r,
                                           const datamatrix & w)
  : DISTR(o,r,w)

  {

  master = mast;
  catnr = cnr;

  if (master==true)
    {

    nrcat=2;
    nrothercat = 0;

    if (check_weightsone() == true)
      wtype = wweightschange_weightsone;
    else
      wtype = wweightschange_weightsneqone;

    }
  else
    wtype = wweightschange_weightsone;

  family = "Multinomial probit";

  updateIWLS = false;
  }




const DISTR_multinomprobit & DISTR_multinomprobit::operator=(
                                      const DISTR_multinomprobit & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  responsecat = nd.responsecat;
  master=nd.master;
  othercat = nd.othercat;
  nrcat = nd.nrcat;
  nrothercat = nd.nrothercat;
  catnr = nd.catnr;
  return *this;
  }


DISTR_multinomprobit::DISTR_multinomprobit(const DISTR_multinomprobit & nd)
   : DISTR(DISTR(nd))
  {
  responsecat = nd.responsecat;
  master=nd.master;
  othercat = nd.othercat;
  nrcat = nd.nrcat;
  nrothercat = nd.nrothercat;
  catnr = nd.catnr;
  }


void DISTR_multinomprobit::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: standard normal (probit link)\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTR_multinomprobit::loglikelihood(double * response, double * linpred,
                                     double * weight)
  {
/*
  if (*weight!=0)
    {
    double mu = randnumbers::Phi2(*linpred);
    if (*response > 0)
      return log(mu);
    else
      return log(1-mu);
    }
  else
*/    return 0;

  }



double DISTR_multinomprobit::loglikelihood_weightsone(
                                  double * response, double * linpred)
  {
/*
  double mu = randnumbers::Phi2(*linpred);
  if (*response > 0)
    return log(mu);
  else
    return log(1-mu);
  */
  return 0;
  }


void DISTR_multinomprobit::compute_mu(const double * linpred,double * mu)
  {
  // *mu = randnumbers::Phi2(*linpred);
  }


void DISTR_multinomprobit::compute_deviance(const double * response,
                   const double * weight,const double * mu,double * deviance,
                   double * scale) const
  {
/*
  if (*weight !=  0)
    {

    if (*response<=0)
      {
      *deviance = -2*log(1-*mu);
      *deviancesat = *deviance;
      }
    else if (*response > 0)
      {
      *deviance = -2*log(*mu);
      *deviancesat = *deviance;
      }

    }
  else
    {
    *deviance = 0;
    *deviancesat = 0;
    }
  */
  }


double DISTR_multinomprobit::compute_iwls(double * response, double * linpred,
                           double * weight, double * workingweight,
                           double * workingresponse, const bool & like)
  {

  double  mu = randnumbers::Phi2(*linpred);

  double h = 0.39894228*exp(-0.5 * *linpred * *linpred);
  double g = 1/pow(h,2);

  *workingweight = *weight / (mu*(1-mu) * g);


  *workingresponse = *linpred + (*response - mu)/h;

  if (like)
    {

    if (*response > 0)
      return log(mu);
    else
      return log(1-mu);
    }
  else
    {
    return 0;
    }

  }



void DISTR_multinomprobit::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  double  mu = randnumbers::Phi2(*linpred);
  double h = 0.39894228*exp(-0.5 * *linpred * *linpred);
  double g = 1/pow(h,2);

  *workingweight = 1.0 / (mu*(1-mu) * g);

  *workingresponse = *linpred + (*response - mu)/h;

  if (compute_like)
    {

    if (*response > 0)
      like+= log(mu);
    else
      like+= log(1-mu);
    }

  }


void DISTR_multinomprobit::compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  }

void DISTR_multinomprobit::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {

  }





double DISTR_multinomprobit::maxutility(vector<datamatrix*> responsep,
const unsigned & i, const unsigned & cat)
  {
  unsigned j;
  double max = 0;
  double help;


  for (j=0;j<=nrothercat;j++)
    {
    help = (*responsep[j])(i,0);
    if ( (j != cat) && (help > max) )
      max = help;
    }

  return max;
  }



void DISTR_multinomprobit::update(void)
  {

  if (optionsp->nriter==1)
    {
    workingweight = weight;
    }

  if (optionsp->nriter==2)
    {

    if (check_weightsone() == true)
      wtype = wweightsnochange_one;
    else
      wtype = wweightsnochange_constant;

    }

  if (master==true)
    {
    unsigned i,j;


    vector<datamatrix *> worklin;
    for (j=0;j<nrothercat;j++)
      {
      if (othercat[j]->linpred_current==1)
        worklin.push_back(&othercat[j]->linearpred1);
      else
        worklin.push_back(&othercat[j]->linearpred2);
      }

    if (linpred_current==1)
      worklin.push_back(&linearpred1);
    else
      worklin.push_back(&linearpred2);


    vector<datamatrix *> responsep;
    for (j=0;j<nrothercat;j++)
      {
      responsep.push_back(&othercat[j]->workingresponse);
      }
    responsep.push_back(&workingresponse);

    double * wresponsecat = responsecat.getV();

    double lin;

    for (i=0;i<nrobs;i++,wresponsecat++)
      {
      if (*wresponsecat == -1)   // reference category
        {

        for (j=0;j<=nrothercat;j++)
          {
          lin = (*worklin[j])(i,0);
          (*responsep[j])(i,0) = lin+truncnormal(-20-lin,-lin);
          }

        }
      else
        {
        lin = (*worklin[*wresponsecat])(i,0);
        (*responsep[*wresponsecat])(i,0) = lin + truncnormal(maxutility(responsep,i,*wresponsecat) - lin,20-lin);

        for (j=0;j<=nrothercat;j++)
          {
          if (j != (*wresponsecat))
            {
            lin = (*worklin[j])(i,0);
            (*responsep[j])(i,0) = lin + truncnormal(-20-lin,(*responsep[*wresponsecat])(i,0) - lin);
            }
          }

        }

      }

    // TEST
    /*
    ofstream out("c:\\bayesx\\testh\\results\\utility.raw");
    for (i=0;i<nrobs;i++)
      {
      for (j=0;j<nrcat-1;j++)
        out << (*responsep[j])(i,0) << "   ";
      out << endl;
      }
    */
    // TEST

    } // end: if (master==true)


  DISTR::update();

  }



//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_multinomlogit- --------------------
//------------------------------------------------------------------------------

DISTR_multinomlogit::DISTR_multinomlogit(GENERAL_OPTIONS * o,unsigned & cnr,
                                         bool mast,
                                         const datamatrix & r,
                                         const datamatrix & w)
 	: DISTR_multinomprobit(o, cnr, mast, r, w)
	{

  family = "Multinomial logit";
  H = 6; // vorläufig fix vorgegeben
// Tabellenwerte aus Frühwirth-Schnatter
  SQ = datamatrix(6,5,0); // entspricht 1/s_r^2

// H = 2:
  SQ(0,0) = 1/1.6927;
  SQ(1,0) = 1/5.2785;

// H = 3:
  SQ(0,1) = 1/1.2131;
  SQ(1,1) = 1/2.9955;
  SQ(2,1) = 1/7.5458;

// H = 4:
  SQ(0,2) = 1/0.95529;
  SQ(1,2) = 1/2.048;
  SQ(2,2) = 1/4.4298;
  SQ(3,2) = 1/9.701;

// H = 5:
  SQ(0,3) = 1/0.79334;
  SQ(1,3) = 1/1.5474;
  SQ(2,3) = 1/3.012;
  SQ(3,3) = 1/5.9224;
  SQ(4,3) = 1/11.77;

// H = 6:
  SQ(0,4) = 1/0.68159;
  SQ(1,4) = 1/1.2419;
  SQ(2,4) = 1/2.2388;
  SQ(3,4) = 1/4.0724;
  SQ(4,4) = 1/7.4371;
  SQ(5,4) = 1/13.772;


  weights_mixed = datamatrix(6,5,0);  // enstpricht w_r

// H = 2:
  weights_mixed(0,0) = 0.56442;
  weights_mixed(1,0) = 0.43558;

// H = 3:
  weights_mixed(0,1) = 0.2522;
  weights_mixed(1,1) = 0.58523;
  weights_mixed(2,1) = 0.16257;

// H = 4:
  weights_mixed(0,2) = 0.1065;
  weights_mixed(1,2) = 0.45836;
  weights_mixed(2,2) = 0.37419;
  weights_mixed(3,2) = 0.060951;

// H = 5:
  weights_mixed(0,3) = 0.044333;
  weights_mixed(1,3) = 0.29497;
  weights_mixed(2,3) = 0.42981;
  weights_mixed(3,3) = 0.20759;
  weights_mixed(4,3) = 0.023291;

// H = 6:
  weights_mixed(0,4) = 0.018446;
  weights_mixed(1,4) = 0.17268;
  weights_mixed(2,4) = 0.37393;
  weights_mixed(3,4) = 0.31697;
  weights_mixed(4,4) = 0.1089;
  weights_mixed(5,4) = 0.0090745;

	}


DISTR_multinomlogit::DISTR_multinomlogit(const DISTR_multinomlogit & nd)
	: DISTR_multinomprobit(DISTR_multinomprobit(nd))
	{
    H = nd.H;
	SQ = nd.SQ;
	weights_mixed = nd.weights_mixed;
	}

const DISTR_multinomlogit & DISTR_multinomlogit::operator=(const DISTR_multinomlogit & nd)
	{
	if (this==&nd)
  	return *this;
  DISTR_multinomprobit::operator=(DISTR_multinomprobit(nd));
  H = nd.H;
  SQ = nd.SQ;
  weights_mixed = nd.weights_mixed;
  return *this;
	}

void DISTR_multinomlogit::outoptions()   //output for logfile
  {
  DISTR::outoptions();
  optionsp->out("  Response function: logistic distribution function\n");
  optionsp->out("  Number of mixture components: " + ST::inttostring(H) + "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_multinomlogit::update(void)
  {


// zum testen was drinsteht Test Anfang
  if (master)
    {

 /* ofstream out1("c:\\testhn\\lin_master.raw");
  if (linpred_current==1)
  linearpred1.prettyPrint(out1);
  else
  linearpred2.prettyPrint(out1);
  }
  else
    {
  ofstream out1("c:\\testhn\\lin.raw");
  if (linpred_current==1)
  linearpred1.prettyPrint(out1);
  else
  linearpred2.prettyPrint(out1);
  }
*/

  ofstream out1("c:\\testhn\\res_master.raw");
  workingresponse.prettyPrint(out1);
  }
  else
    {
  ofstream out1("c:\\testhn\\res.raw");
  workingresponse.prettyPrint(out1);
  }

  //Test Ende


  if (master == true)
  {
  unsigned j;
// linear predictor
  vector<datamatrix *> worklin;
  for (j=0;j<nrothercat;j++)
    {
    if (othercat[j]->linpred_current==1)
      worklin.push_back(&othercat[j]->linearpred1);
    else
      worklin.push_back(&othercat[j]->linearpred2);
    }

  if (linpred_current==1)
    worklin.push_back(&linearpred1);
  else
    worklin.push_back(&linearpred2);

// response
  vector<datamatrix *> responsep;
  for (j=0;j<nrothercat;j++)
    {
    responsep.push_back(&othercat[j]->workingresponse);
    }
  responsep.push_back(&workingresponse);

// working weights
  vector<datamatrix *> workingweightp;
  for (j=0;j<nrothercat;j++)
    {
    workingweightp.push_back(&othercat[j]->workingweight);
    }
  workingweightp.push_back(&workingweight);

// original response
  vector<datamatrix *> origrespp;
  for (j=0;j<nrothercat;j++)
    {
    origrespp.push_back(&othercat[j]->response);
    }
  origrespp.push_back(&response);  // response 0,1 je nach Kategorie

  unsigned i;  // nrobs
  int k;
  int l;
  double sumpred; // sum exp(x_i* beta_k)
  double uni;
  double uni2;
  double ratio;
  double elin;
  double hresp;
  double hprop;
  datamatrix rvektor(H,1,0);


  for (i=0;i<nrobs;i++)
    {
    sumpred = 0;
    for (j=0;j<nrcat-1;j++)
      {
      sumpred += exp((*worklin[j])(i,0));   // sum Prediktoren
      }
    for (j=0;j<nrcat-1;j++)
      {
      uni = uniform();
      elin = exp((*worklin[j])(i,0)); // xi*betak
      ratio = elin / (sumpred - elin);
      hresp = (*origrespp[j])(i,0);  // 0 or 1
      (*responsep[j])(i,0) = log(ratio*uni + hresp) - log(1 - uni + ratio*(1-hresp));
      hprop = -0.5 *  (pow((*responsep[j])(i,0) - (*worklin[j])(i,0),2));

/*
     ofstream out("c:\\testhn\\rvektor.raw");  // Test
     cout << "sumpred " << sumpred << "\n"; // test
     cout << "elin " <<elin<< "\n"; // test
     cout << "ratio "<<ratio<< "\n"; // test
     cout << "uni " << uni << "\n"; // test
     cout << "responsep " << (*responsep[j])(i,0) << "\n"; // test
     cout << "origrespp[j])(i,0)" << (*origrespp[j])(i,0) << "\n"; // test
     cout << "hresp "<<hresp<<"\n"; // test
     cout << "hprop "<<hprop<<"\n"; // test
*/

      for (k=0;k<H;k++)
        {
        rvektor(k,0) = weights_mixed(k,H-2) * sqrt(SQ(k,H-2)) * exp(hprop*SQ(k,H-2));
        }


//     rvektor.prettyPrint(out);  // test
//      out <<   endl;  // Test

      for(int k=1;k<H;k++)
        {
        rvektor(k,0) = rvektor(k-1,0) + rvektor(k,0);
        }

//      rvektor.prettyPrint(out); // Test
//      out <<   endl;  // Test

      for(int k=0;k<H;k++)
        {
        rvektor(k,0) = rvektor(k,0)/rvektor(H-1,0);
        }


 //     rvektor.prettyPrint(out); // Test


      uni2 = uniform();
      l = 0;
      while (uni2> rvektor(l,0))
        {
        l++;
        }
      (*workingweightp[j])(i,0) = SQ(l,H-2);
 //     cout << "ww " << (*workingweightp[j])(i,0) << "\n"; // test
      }
    }
  }
  }

/*

   double lin;

    for (i=0;i<nrobs;i++)
      {
      if (responsecat(i,0) == -1)   // reference category
        {

        for (j=0;j<=nrothercat;j++)
          {
          lin = (*worklin[j])(i,0);
          (*responsep[j])(i,0) = lin+truncnormal(-20-lin,-lin);
          }

        }
      else
        {
        lin = (*worklin[responsecat(i,0)])(i,0);
        (*responsep[responsecat(i,0)])(i,0) = lin + truncnormal(maxutility(responsep,i,responsecat(i,0)) - lin,20-lin);

        for (j=0;j<=nrothercat;j++)
          {
          if (j != responsecat(i,0))
            {
            lin = (*worklin[j])(i,0);
            (*responsep[j])(i,0) = lin + truncnormal(-20-lin,(*responsep[responsecat(i,0)])(i,0) - lin);
            }
          }

        }
      }


    }


*/

} // end: namespace MCMC



