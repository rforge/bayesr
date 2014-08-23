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




#include "distr.h"

namespace MCMC
{

void DISTR::get_samples(const ST::string & filename,ofstream & outg) const
  {


  }


void DISTR::check_errors(void)
  {
  if (response.var(0) < 0.0000000001)
    {
    errormessages.push_back("ERROR: response is not varying \n");
    errors=true;
    }
  else
    errors=false;

  }

bool DISTR::check_weightsone(void)
  {
  unsigned i=0;
  double * work_weight = weight.getV();
  bool one = true;
  while (i<nrobs && one == true)
    {
    if (*work_weight != 1)
      one = false;
    work_weight++;
    i++;
    }
  return one;
  }


unsigned DISTR::compute_nrzeroweights(void)
  {
  unsigned i=0;
  double * work_weight = weight.getV();
  unsigned nr=0;
  for (i=0;i<weight.rows();i++,work_weight++)
    {
    if (*work_weight == 0)
      nr++;
    }
  return nr;
  }


DISTR::DISTR(GENERAL_OPTIONS * o, const datamatrix & r,
             const datamatrix & w)
  {

  maindistribution=true;
  predict_mult=false;

  option1 = "";
  offsetname = "";

  sigma2=1;

  optionsp = o;

  family = "unknown";
  familyshort = "unknown";
  updateIWLS = false;

  response = r;
  workingresponse = r;
  responsename = "Y";

  nrobs = response.rows();

  if (w.rows() == 1)
    {
    weight = datamatrix(r.rows(),1,1);
    }
  else
    {
    weight = w;
    }

  workingweight = weight;

  weightsone = check_weightsone();

  nrzeroweights = compute_nrzeroweights();

  wtype = wweightschange_weightsneqone;

  weightname = "W";

  linearpred1 = datamatrix(nrobs,1,0);
  linearpred2 = datamatrix(nrobs,1,0);

  linpred_current = 1;

  trmult=1;

  meaneffect = 0;

  linpredminlimit = -1000000000;
  linpredmaxlimit =  1000000000;

  check_errors();

  outpredictor = true;
  outexpectation = false;
  predictor_name = "";
  predstart_mumult = 0;

  }


DISTR::DISTR(const DISTR & d)
  {

  FCpredict_betamean = d.FCpredict_betamean;

  maindistribution = d.maindistribution;
  predict_mult = d.predict_mult;

  option1 = d.option1;
  optionbool1 = d.optionbool1;

  sigma2 = d.sigma2;

  optionsp = d.optionsp;
  nrobs = d.nrobs;

  response = d.response;
  workingresponse = d.workingresponse;
  responsename = d.responsename;
  offsetname = d.offsetname;

  weight = d.weight;
  weightname = d.weightname;
  nrzeroweights = d.nrzeroweights;

  workingweight = d.workingweight;
  weightsone = d.weightsone;
  wtype = d.wtype;

  linearpred1 = d.linearpred1;
  linearpred2 = d.linearpred2;
  linpred_current = d.linpred_current;

  helpmat1 = d.helpmat1;
  helpmat2 = d.helpmat2;
  helpmat3 = d.helpmat3;

  helpquantity1 = d.helpquantity1;
  helpquantity2 = d.helpquantity2;
  helpquantity3 = d.helpquantity3;

  updateIWLS = d.updateIWLS;
  family = d.family;
  familyshort = d.familyshort;
  equationtype = d.equationtype;
  hlevel = d.hlevel;

  trmult=d.trmult;

  meaneffect = d.meaneffect;

  errors = d.errors;
  errormessages = d.errormessages;

  linpredminlimit = d.linpredminlimit;
  linpredmaxlimit = d.linpredmaxlimit;

  outpredictor = d.outpredictor;
  outexpectation = d.outexpectation;
  predictor_name = d.predictor_name;
  predstart_mumult = d.predstart_mumult;
  }


const DISTR & DISTR::operator=(const DISTR & d)
  {
  if (this == &d)
    return *this;

  FCpredict_betamean = d.FCpredict_betamean;

  maindistribution = d.maindistribution;
  predict_mult = d.predict_mult;

  option1 = d.option1;
  optionbool1 = d.optionbool1;

  sigma2 = d.sigma2;

  optionsp = d.optionsp;
  nrobs = d.nrobs;

  response = d.response;
  workingresponse = d.workingresponse;
  responsename = d.responsename;
  offsetname = d.offsetname;

  weight = d.weight;
  weightname = d.weightname;
  nrzeroweights = d.nrzeroweights;

  workingweight = d.workingweight;
  weightsone = d.weightsone;
  wtype = d.wtype;

  linearpred1 = d.linearpred1;
  linearpred2 = d.linearpred2;
  linpred_current = d.linpred_current;

  helpmat1 = d.helpmat1;
  helpmat2 = d.helpmat2;
  helpmat3 = d.helpmat3;

  helpquantity1 = d.helpquantity1;
  helpquantity2 = d.helpquantity2;
  helpquantity3 = d.helpquantity3;

  updateIWLS = d.updateIWLS;
  family = d.family;
  equationtype = d.equationtype;
  hlevel = d.hlevel;

  trmult=d.trmult;

  meaneffect = d.meaneffect;

  errors = d.errors;
  errormessages = d.errormessages;

  linpredminlimit = d.linpredminlimit;
  linpredmaxlimit = d.linpredmaxlimit;

  outpredictor = d.outpredictor;
  outexpectation = d.outexpectation;
  predictor_name = d.predictor_name;
  predstart_mumult = d.predstart_mumult;

  return *this;
  }


bool DISTR::check_linpred(bool current)
  {

  bool ok = true;

  double* worklin;
  if (current)
    {
    if (linpred_current==1)
      worklin = linearpred1.getV();
    else
      worklin = linearpred2.getV();
    }
  else
    {
    if (linpred_current==1)
      worklin = linearpred2.getV();
    else
      worklin = linearpred1.getV();
    }

  unsigned i=0;
  while (ok && (i<nrobs))
    {
    if (*worklin > linpredmaxlimit)
      ok = false;
    if (*worklin < linpredminlimit)
      ok = false;

    worklin++;
    i++;
    }

  return ok;
  }


void DISTR::changelimits(double min,double max)
  {
  linpredminlimit=min;
  linpredmaxlimit=max;
  }


void DISTR::outoptions(void)
  {
  optionsp->out("RESPONSE DISTRIBUTION:\n",true);
  optionsp->out("\n");
  optionsp->out("  Family: " + family + "\n");
  optionsp->out("  Number of observations: " + ST::inttostring(nrobs) + "\n");
  if (offsetname.length() > 0)
    optionsp->out("  Offset: " + offsetname + "\n");


    optionsp->out("  Number of observations with positive weights: " +
      ST::inttostring(nrobs-nrzeroweights) + "\n");

      optionsp->out("\n");


  if (optionsp->saveestimation)
    {
    optionsp->out("  Limits for predictor (save estimation mode):\n");
    optionsp->out("    Minimum: " + ST::doubletostring(linpredminlimit,6) + "\n");
    optionsp->out("    Maximum: " + ST::doubletostring(linpredmaxlimit,6) + "\n");
    }
  }


double DISTR::get_intercept_start(void)
  {
  return 0;
  }


double DISTR::compute_quantile_residual(double * res,double * param,double * weight,
                                        double * scale)
  {
  if ((*weight) == 0)
  {
    return  0;
  }
  else
  {
  double u_est = cdf(res,param,weight,scale);
  double res_est = randnumbers::invPhi2(u_est);
  return res_est;
  }

  }


double DISTR::compute_quantile_residual_mult(vector<double *> response,
                                             vector<double *> param,
                                             vector<double *> weight,
                                             vector<datamatrix *> aux)
  {

    double u_est = cdf_mult(response,param,weight,aux);
    double res_est = randnumbers::invPhi2(u_est);
    return res_est;

  }


double DISTR::compute_quadr(void)
  {
  return 0;
  }


double DISTR::compute_quadr_mult(void)
  {
  return 0;
  }

double DISTR::compute_log(double * res,double * param,double * weight,
                                        double * scale)
  {

        double result = log(pdf(res,param,weight,scale));
        return result;


  }


double DISTR::compute_log_mult(vector<double *> response,
                               vector<double *> param,
                               vector<double *> weight,
                               vector<datamatrix *> aux)
  {

    double result = log(pdf_mult(response,param,weight,aux));
    return result;


  }
double DISTR::compute_spherical(void)
  {
  return 0;
  }


double DISTR::compute_spherical_mult(void)
  {
  return 0;
  }
double DISTR::compute_CRPS(void)
  {
  return 0;
  }


double DISTR::compute_CRPS_mult(void)
  {
  return 0;
  }



double DISTR::loglikelihood(const bool & current)
  {

  register unsigned  i;
  double* workweight = weight.getV();
  double* workres = response.getV();
  double help = 0;

  double* worklin;
  if (current)
    {
    if (linpred_current==1)
      worklin = linearpred1.getV();
    else
      worklin = linearpred2.getV();
    }
  else
    {
    if (linpred_current==1)
      worklin = linearpred2.getV();
    else
      worklin = linearpred1.getV();
    }

  if (weightsone==true)
    {
    for (i=0;i<nrobs;i++,worklin++,workres++)
      help += loglikelihood_weightsone(workres,worklin);
    }
  else
    {
    for (i=0;i<nrobs;i++,workweight++,worklin++,workres++)
      help += loglikelihood(workres,worklin,workweight);
    }

  return help;

  }


double DISTR::loglikelihood(int & begin,
int & end, statmatrix<double *> & responsep,
statmatrix<double *> & workingweightp, statmatrix<double *> & linpredp)
  {
  double help=0;
  int i;

  double * * workresponsep = responsep.getV()+begin;
  double * * work_workingweightp = workingweightp.getV()+begin;
  double * * work_linpredp  = linpredp.getV()+begin;

  for (i=begin;i<=end;i++,work_workingweightp++,work_linpredp++,workresponsep++)
    help += loglikelihood(*workresponsep,*work_linpredp,*work_workingweightp);

  return help;


  }


double DISTR::compute_iwls_loglikelihood(int & begin,
int & end, statmatrix<double *> & responsep,
statmatrix<double *> & workingresponsep,
statmatrix<double *> & weightp,
statmatrix<double *> & workingweightp, statmatrix<double *> & linpredp
)
  {
  double help=0;
  int i;

  double * * workresponsep = responsep.getV()+begin;
  double * * work_workingresponsep = workingresponsep.getV()+begin;
  double * * work_weightp = weightp.getV()+begin;
  double * * work_workingweightp = workingweightp.getV()+begin;
  double * * work_linpredp  = linpredp.getV()+begin;

  for (i=begin;i<=end;i++,work_workingweightp++,work_linpredp++,
       workresponsep++,work_weightp++,work_workingresponsep++)
    {

    help += compute_iwls(*workresponsep,*work_linpredp,*work_weightp,
                         *work_workingweightp,*work_workingresponsep,true);
    }

  return help;


  }



double DISTR::compute_iwls_loglikelihood_sumworkingweight(int & begin,
int & end, statmatrix<double *> & responsep,
statmatrix<double *> & workingresponsep,
statmatrix<double *> & weightp,
statmatrix<double *> & workingweightp, statmatrix<double *> & linpredp,
datamatrix & intvar2,
double & sumworkingweight)
  {
  double help=0;
  int i;

  double * * workresponsep = responsep.getV()+begin;
  double * * work_workingresponsep = workingresponsep.getV()+begin;
  double * * work_weightp = weightp.getV()+begin;
  double * * work_workingweightp = workingweightp.getV()+begin;
  double * * work_linpredp  = linpredp.getV()+begin;

  sumworkingweight = 0;


  if (wtype==wweightschange_weightsneqone)
    {

    for (i=begin;i<=end;i++,work_workingweightp++,work_linpredp++,
         workresponsep++,work_weightp++,work_workingresponsep++)
      {

      help += compute_iwls(*workresponsep,*work_linpredp,*work_weightp,
                         *work_workingweightp,*work_workingresponsep,true);

      sumworkingweight+= *(*work_workingweightp);
      }

    }
  else if (wtype==wweightschange_weightsone)
    {

    for (i=begin;i<=end;i++,work_workingweightp++,work_linpredp++,
         workresponsep++,work_workingresponsep++)
      {

      compute_iwls_wweightschange_weightsone(*workresponsep,
                         *work_linpredp, *work_workingweightp,
                         *work_workingresponsep,help,true);

      sumworkingweight+= *(*work_workingweightp);
      }

    }
  else if (wtype==wweightsnochange_constant)
    {

    for (i=begin;i<=end;i++,work_workingweightp++,work_linpredp++,
         workresponsep++,work_workingresponsep++)
      {

      compute_iwls_wweightsnochange_constant(*workresponsep,
                         *work_linpredp, *work_workingweightp,
                         *work_workingresponsep,help,true);

      sumworkingweight+= *(*work_workingweightp);
      }

    }
  else if (wtype==wweightsnochange_one)
    {


    for (i=begin;i<=end;i++,work_linpredp++,
         workresponsep++,work_workingresponsep++)
      {

      compute_iwls_wweightsnochange_one(*workresponsep,
                         *work_linpredp, *work_workingresponsep,help,true);

      }

    sumworkingweight = begin-end+1;

    }

  return help;

  }


void DISTR::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * scale) const
  {

  }


void DISTR::compute_deviance_mult(vector<double *> response,
                                  vector<double *> weight,
                                  vector<double *> linpred,
                                  double * deviance,
                                  vector<datamatrix *> aux)
  {

  }



void DISTR::compute_mu(const double * linpred,double * mu)
  {

  }


void DISTR::compute_mu_mult(vector<double *> linpred,vector<double *> response,double * mu)
  {

  }


void DISTR::compute_param(const double * linpred,double * param)
  {
  compute_mu(linpred,param);
  }


void DISTR::compute_param_mult(vector<double *> linpred,double * param)
  {

  }

double DISTR::compute_MSE(const double * response, const double * weight,
                          const double * linpred, msetype t, double v)
  {
  if (t == quadraticMSE)
    return pow(*response-*linpred,2);
  else
    {
    double u = *response-*linpred;
    if (u >= 0)
      return u*v;
    else
      return u*(v-1);
    }
  }


void DISTR::compute_MSE_all(datamatrix & meanpred, double & MSE,
                            double & MSEzeroweight, unsigned & nrzeroweights,
                            msetype & t, double & v)
  {
  unsigned i;
  nrzeroweights = 0;
  MSE = 0;
  MSEzeroweight=0;
  double * responsep = response.getV();
  double * weightp = weight.getV();
  double * linpredp = meanpred.getV();
  for(i=0;i<nrobs;i++,responsep++,weightp++,linpredp+=2)
    if (*weightp==0)
      {
      MSEzeroweight += compute_MSE(responsep,weightp,linpredp,t,v);
      nrzeroweights++;
      }
    else
      MSE += compute_MSE(responsep,weightp,linpredp,t,v);
  }


void DISTR::swap_linearpred(void)
  {

  if (linpred_current==1)
    linpred_current=2;
  else
    linpred_current=1;
  }


void DISTR::update(void)
  {
  } // end: update


void DISTR::update_end(void)
  {
  } // end: update_end


bool DISTR::posteriormode(void)
  {
  double h = compute_iwls(true,false);
  return true;
  }

void DISTR::posteriormode_end(void)
  {
  }

void DISTR::update_scale_hyperparameters(datamatrix & h)
  {

  }


double DISTR::compute_iwls(const bool & current, const bool & like)
  {

  register unsigned  i;

  double * workweight = weight.getV();
  double * workresponse = response.getV();

  double * worklin;
  if (current)
    {
    if (linpred_current == 1)
      worklin = linearpred1.getV();
    else
      worklin = linearpred2.getV();
    }
  else  // use porposed
    {
    if (linpred_current == 1)
      worklin = linearpred2.getV();
    else
      worklin = linearpred1.getV();
    }

  double * work_workingresponse=workingresponse.getV();
  double * work_workingweight = workingweight.getV();

  double likelihood = 0;

  if (wtype==wweightschange_weightsneqone)
    {

    for (i=0;i<nrobs;i++,workweight++,work_workingweight++,workresponse++,
          work_workingresponse++,worklin++)
      {

      likelihood += compute_iwls(workresponse,worklin,
                                 workweight,work_workingweight,
                                 work_workingresponse,like);
      }

    }
  else if (wtype==wweightschange_weightsone)
    {

    for (i=0;i<nrobs;i++,work_workingweight++,workresponse++,
          work_workingresponse++,worklin++)
      {

      compute_iwls_wweightschange_weightsone(workresponse,worklin,
                                 work_workingweight,work_workingresponse,
                                 likelihood,like);
      }
/*
    ofstream out("d:\\_sicher\\papzip\\wweight.raw");
    workingweight.prettyPrint(out);
    out.close();

    ofstream out2("d:\\_sicher\\papzip\\wresponse.raw");
    workingresponse.prettyPrint(out2);
    out2.close();
*/
    }
  else if (wtype==wweightsnochange_constant)
    {

    for (i=0;i<nrobs;i++,work_workingweight++,workresponse++,
          work_workingresponse++,worklin++)
      {

      compute_iwls_wweightsnochange_constant(workresponse,worklin,
                                 work_workingweight,work_workingresponse,
                                 likelihood,like);
      }

    }
  else if (wtype==wweightsnochange_one)
    {

    for (i=0;i<nrobs;i++,workresponse++,
          work_workingresponse++,worklin++)
      {

      compute_iwls_wweightsnochange_one(workresponse,worklin,
                                        work_workingresponse,
                                        likelihood,like);
      }

    }

  // TEST

/*
  ofstream out("c:\\bayesx\\testh\\results\\workresponse.res");
  workingresponse.prettyPrint(out);

  ofstream out2("c:\\bayesx\\testh\\results\\workweight.res");
  workingweight.prettyPrint(out2);

  ofstream out3("c:\\bayesx\\testh\\results\\linpred.res");
  linearpred1.prettyPrint(out3);
  */

  // TEST

  return likelihood;
  }




void DISTR::compute_iwls(const bool & current,datamatrix & likelihood,
                    statmatrix<unsigned> & ind)
  {

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\ind.res");
  // ind.prettyPrint(out);
  // TEST

  register unsigned  i;

  double * workweight = weight.getV();
  double * workresponse = response.getV();

  double * worklin;
  if (current)
    {
    if (linpred_current == 1)
      worklin = linearpred1.getV();
    else
      worklin = linearpred2.getV();
    }
  else  // use porposed
    {
    if (linpred_current == 1)
      worklin = linearpred2.getV();
    else
      worklin = linearpred1.getV();
    }

  double * work_workingresponse=workingresponse.getV();
  double * work_workingweight = workingweight.getV();

  unsigned * workind = ind.getV();

  double * worklikelihood=likelihood.getV();
  for (i=0;i<likelihood.rows();i++,worklikelihood++)
    *worklikelihood=0;

  if (wtype==wweightschange_weightsneqone)
    {

    for (i=0;i<nrobs;i++,workweight++,work_workingweight++,workresponse++,
          work_workingresponse++,worklin++,workind++)
      {

      likelihood(*workind,0) += compute_iwls(workresponse,worklin,
                                 workweight,work_workingweight,
                                 work_workingresponse,true);
      }

    }
  else if (wtype==wweightschange_weightsone)
    {

    for (i=0;i<nrobs;i++,work_workingweight++,workresponse++,
          work_workingresponse++,worklin++,workind++)
      {

      compute_iwls_wweightschange_weightsone(workresponse,worklin,
                                 work_workingweight,work_workingresponse,
                                 likelihood(*workind,0),true);
      }


    }
  else if (wtype==wweightsnochange_constant)
    {

    for (i=0;i<nrobs;i++,work_workingweight++,workresponse++,
          work_workingresponse++,worklin++,workind++)
      {

      compute_iwls_wweightsnochange_constant(workresponse,worklin,
                                 work_workingweight,work_workingresponse,
                                 likelihood(*workind,0),true);
      }

    }
  else if (wtype==wweightsnochange_one)
    {

    for (i=0;i<nrobs;i++,workresponse++,
          work_workingresponse++,worklin++,workind++)
      {

      compute_iwls_wweightsnochange_one(workresponse,worklin,
                                        work_workingresponse,
                                        likelihood(*workind,0),true);
      }

    }

  // TEST
  /*
  ofstream out("c:\\bayesx\\test\\results\\workresponse.res");
  workingresponse.prettyPrint(out);

  ofstream out2("c:\\bayesx\\test\\results\\workweight.res");
  workingweight.prettyPrint(out2);
  */
  // TEST

  }



void DISTR::outresults(ofstream & out_stata, ofstream & out_R,ST::string pathresults)
  {
  optionsp->out("\n");
  }


void DISTR::reset(void)
  {
  linearpred1 = datamatrix(nrobs,1,0);
  linearpred2 = datamatrix(nrobs,1,0);
  linpred_current = 1;
  }


double DISTR::get_scale(void)
  {
  return sigma2;
  }

datamatrix * DISTR::get_auxiliary_parameter(auxiliarytype t)
  {
  return &helpmat1;
  }

double DISTR::get_scalemean(void)
  {
  return sigma2;
  }


void DISTR::sample_responses(unsigned i,datamatrix & sr)
  {
  }


void DISTR::sample_responses_cv(unsigned i,datamatrix & linpred, datamatrix & sr)
  {
  }


void DISTR::outresults_predictive_check(datamatrix & D,datamatrix & sr)
  {
  unsigned j;
  datamatrix h(sr.cols(),1,0);

  optionsp->out("\n");

  for (j=1;j<h.rows();j++)
    h(j,0) = sr.mean(j);

  optionsp->out("    Mean           " + ST::doubletostring(D.mean(0)) +  "  " +
                ST::doubletostring(h.quantile(5,0)) + " - "
                + ST::doubletostring(h.quantile(95,0)) + "\n" ) ;


  for (j=1;j<h.rows();j++)
    h(j,0) = sqrt(sr.var(j));

  optionsp->out("    Std.Dev        " + ST::doubletostring(sqrt(D.var(0))) +  "  " +
                ST::doubletostring(h.quantile(5,0)) + " - "
                + ST::doubletostring(h.quantile(95,0)) + "\n" ) ;


  for (j=1;j<h.rows();j++)
    h(j,0) = sr.min(j);
  optionsp->out("    Minimum        " + ST::doubletostring(D.min(0)) +  "  " +
                ST::doubletostring(h.quantile(5,0)) + " - "
                + ST::doubletostring(h.quantile(95,0)) + "\n" ) ;


  for (j=1;j<h.rows();j++)
    h(j,0) = sr.quantile(25,j);
  optionsp->out("    25\% Quantile  " + ST::doubletostring(D.quantile(25,0)) +  "  " +
                ST::doubletostring(h.quantile(5,0)) + " - "
                + ST::doubletostring(h.quantile(95,0)) + "\n" ) ;

  for (j=1;j<h.rows();j++)
    h(j,0) = sr.quantile(50,j);
  optionsp->out("    50\% Quantile  " + ST::doubletostring(D.quantile(50,0)) +  "  " +
                ST::doubletostring(h.quantile(5,0)) + " - "
                + ST::doubletostring(h.quantile(95,0)) + "\n" ) ;

  for (j=1;j<h.rows();j++)
    h(j,0) = sr.quantile(75,j);
  optionsp->out("    75\% Quantile  " + ST::doubletostring(D.quantile(75,0)) +  "  " +
                ST::doubletostring(h.quantile(5,0)) + " - "
                + ST::doubletostring(h.quantile(95,0)) + "\n" ) ;

  for (j=1;j<h.rows();j++)
    h(j,0) = sr.max(j);
  optionsp->out("    Maximum        " + ST::doubletostring(D.max(0)) +  "  " +
                ST::doubletostring(h.quantile(5,0)) + " - "
                + ST::doubletostring(h.quantile(95,0)) + "\n" ) ;


  optionsp->out("\n");
  }


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_gaussian --------------------------
//------------------------------------------------------------------------------


DISTR_gaussian::DISTR_gaussian(const double & a,
                                             const double & b,
                                             GENERAL_OPTIONS * o,
                                             const datamatrix & r,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTR(o,r,w)

  {

  predictor_name = "mu";
  outexpectation = true;

  lassosum = 0;
  ridgesum = 0;
  nrlasso=0;
  nrridge=0;

  if (check_weightsone())
    wtype = wweightsnochange_one;
  else
    wtype = wweightsnochange_constant;

  a_invgamma = a;
  double h = sqrt(response.var(0,weight));
  b_invgamma = b*h;
  trmult = h;
  family = "Normal distribution with homoscedastic variance";

  FCsigma2 = FC(o,"",1,1,ps);

  }


const DISTR_gaussian & DISTR_gaussian::operator=(
                                      const DISTR_gaussian & nd)
  {
  if (this==&nd)
    return *this;
  DISTR::operator=(DISTR(nd));
  a_invgamma = nd.a_invgamma;
  b_invgamma = nd.b_invgamma;
  lassosum = nd.lassosum;
  ridgesum = nd.ridgesum;
  nrlasso = nd.nrlasso;
  nrridge = nd.nrridge;
  FCsigma2 = nd.FCsigma2;
  return *this;
  }



DISTR_gaussian::DISTR_gaussian(const DISTR_gaussian & nd)
   : DISTR(DISTR(nd))
  {
  a_invgamma = nd.a_invgamma;
  b_invgamma = nd.b_invgamma;
  lassosum = nd.lassosum;
  ridgesum = nd.ridgesum;
  nrlasso = nd.nrlasso;
  nrridge = nd.nrridge;
  FCsigma2 = nd.FCsigma2;
  }



void DISTR_gaussian::get_samples(const ST::string & filename,ofstream & outg) const
  {
  if (filename.isvalidfile() != 1)
    {
    FCsigma2.get_samples(filename,outg);
    }

  }



void DISTR_gaussian::sample_responses(unsigned i,datamatrix & sr)
  {

  double * linpredp;

  if (linpred_current==1)
    linpredp = linearpred1.getV();
  else
    linpredp = linearpred2.getV();

  double * rp = sr.getV()+i;

  double * wweightp  = workingweight.getV();

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,wweightp++,rp+=sr.cols())
    *rp = *linpredp+sqrt(sigma2)/sqrt(*wweightp)*rand_normal();
  }

// linpred already transformed
void DISTR_gaussian::sample_responses_cv(unsigned i,datamatrix & linpred,
                                         datamatrix & sr)
  {

  double * linpredp;

  linpredp = linpred.getV();

  double * rp = sr.getV()+i;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
    *rp = *linpredp + sqrt(sigma2)*rand_normal();


  }


  // FUNCTION: update_scale_hyperparameters
  // TASK: updates parameters for lasso, ridge etc.
  //       h(0,0) = type, 1 =ridge, 2=lasso
  //       h(1,0) = nrridge/nrlasso
  //       h(2,0) = lassosum/ridgesum

void DISTR_gaussian::update_scale_hyperparameters(datamatrix & h)
  {
   if (h(0,0) == 1)        //  ridge
     {
     nrridge = h(1,0);
     ridgesum = h(2,0);
     }
   else if (h(0,0) == 2)   // lasso
     {
     nrlasso = h(1,0);
     lassosum = h(2,0);
     }
  }


void DISTR_gaussian::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: identity\n");

  optionsp->out("  Hyperparameter a: " + ST::doubletostring(a_invgamma,6) + "\n");
  optionsp->out("  Hyperparameter b: " + ST::doubletostring(b_invgamma,6) + "\n");

  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussian::outresults_predictive_check(datamatrix & D,datamatrix & sr)
  {

  DISTR::outresults_predictive_check(D,sr);
  }


void DISTR_gaussian::update(void)
  {

  register unsigned i;

  double help;

  double * worklin;
  double * workresp;
  double * workweight;

  // scaleparameter

  double sum = 0;

  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  workresp = workingresponse.getV();
  workweight = weight.getV();

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    if (*workweight !=0)
      {
      help = *workresp - *worklin;
      sum += *workweight*pow(help,2);
      }
    }

  sigma2  = rand_invgamma(a_invgamma+0.5*((nrobs-nrzeroweights)+nrlasso+nrridge),
                          b_invgamma+0.5*(sum+lassosum+ridgesum));


  FCsigma2.beta(0,0) = sigma2;
  FCsigma2.acceptance++;
  FCsigma2.update();

  DISTR::update();

  }


void DISTR_gaussian::compute_mu(const double * linpred,double * mu)
  {
  *mu = (*linpred);
  }


void DISTR_gaussian::compute_deviance(const double * response,
                                 const double * weight, const double * mu,
                                 double * deviance,
                                 double * scale) const
  {
  if (*weight == 0)
    {
    *deviance = 0;
    }
  else
    {
    double r = *response-*mu;
    *deviance =  (*weight/(*scale))*r*r+log(2*M_PI*(*scale)/(*weight));
    }
  }


double DISTR_gaussian::compute_MSE(const double * response,
                          const double * weight,
                          const double * linpred, msetype t, double v)
  {
  if (t == quadraticMSE)
    return pow(*response-*linpred,2);
  else
    {
    double u;

    if (*weight == 0)
      {
      u = *response - ( *linpred +
      sqrt(FCsigma2.betamean(0,0))* randnumbers::invPhi2(v) );
      }
    else
      {
      u = *response - ( *linpred +
      sqrt(FCsigma2.betamean(0,0)/(*weight))* randnumbers::invPhi2(v) );
      }

    if (u >= 0)
      return u*v;
    else
      return u*(v-1);
    }
  }


double DISTR_gaussian::get_intercept_start(void)
  {
  return response.mean(0);
  }


double DISTR_gaussian::loglikelihood(double * res, double * lin,
                                     double * w)
  {
  if (*w==0)
    return 0;
  else
    {
    double help = *res-*lin;
    return  - *w * (pow(help,2))/(2* sigma2);
    }
  }


double DISTR_gaussian::loglikelihood_weightsone(double * res, double * lin)
  {
  double help = *res-*lin;
  return  - (pow(help,2))/(2* sigma2);
  }


double DISTR_gaussian::compute_iwls(double * response, double * linpred,
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


void DISTR_gaussian::compute_iwls_wweightschange_weightsone(
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


void DISTR_gaussian::compute_iwls_wweightsnochange_constant(double * response,
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


void DISTR_gaussian::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  *workingresponse = *response;
  if (compute_like)
    like -=  (pow(*response-(*linpred),2))/(2* sigma2);
  }


bool DISTR_gaussian::posteriormode(void)
  {

  unsigned i;

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * workresp = workingresponse.getV();
  double * workweight = weight.getV();

  double sum = 0;
  double sumweight=0;
  double help;

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - *worklin;
    sum += *workweight*pow(help,2);
    sumweight+=*workweight;
    }

  sigma2 = (1.0/sumweight)*sum;

  FCsigma2.beta(0,0) = sigma2;

  FCsigma2.posteriormode_betamean();

  return true;

  }


void DISTR_gaussian::outresults(ofstream & out_stata, ofstream & out_R,
                                ST::string pathresults)
  {
  DISTR::outresults(out_stata,out_R);

  FCsigma2.outresults(out_stata,out_R,"");

  ST::string l1 = ST::doubletostring(optionsp->lower1,4);
  ST::string l2 = ST::doubletostring(optionsp->lower2,4);
  ST::string u1 = ST::doubletostring(optionsp->upper1,4);
  ST::string u2 = ST::doubletostring(optionsp->upper2,4);

  ST::string nl1 = ST::doubletostring(optionsp->lower1,4);
  ST::string nl2 = ST::doubletostring(optionsp->lower2,4);
  ST::string nu1 = ST::doubletostring(optionsp->upper1,4);
  ST::string nu2 = ST::doubletostring(optionsp->upper2,4);
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');

  double help;

  optionsp->out("  SCALE PARAMETER:\n",true);
  optionsp->out("\n");


  ST::string vstr;

  vstr = "    Mean:         ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(FCsigma2.betamean(0,0),6) + "\n");

  if (optionsp->samplesize > 1)
    {
    vstr = "    Std. dev.:    ";
    if (FCsigma2.betavar(0,0) < 0)
      help = 0;
    else
      help = sqrt(FCsigma2.betavar(0,0));
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(help,6) + "\n");

    vstr = "    " + l1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(FCsigma2.betaqu_l1_lower(0,0),6) + "\n");

    vstr = "    " + l2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(FCsigma2.betaqu_l2_lower(0,0),6) + "\n");

    vstr = "    50% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(FCsigma2.betaqu50(0,0),6) + "\n");

    vstr = "    " + u1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(FCsigma2.betaqu_l2_upper(0,0),6) + "\n");

    vstr = "    " + u2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(FCsigma2.betaqu_l1_upper(0,0),6) + "\n");
    }


  if (pathresults.isvalidfile() != 1)
    {

    optionsp->out("\n");

    optionsp->out("    Results for variance parameter are also stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");
//    out_R << "scale=" << pathresults << ";" <<  endl;

    ofstream outscale(pathresults.strtochar());

    if (optionsp->samplesize > 1)
      {
      outscale << "pmean   pstddev   pqu" << nl1 << "   pqu" << nl2 <<
                "   pqu50   pqu" << nu1 << "   pqu" << nu2 << endl;
      }
    else
      {
      outscale << "pmean" << endl;
      }

    outscale << FCsigma2.betamean(0,0) << "  ";

    if (optionsp->samplesize > 1)
      {
      if (FCsigma2.betavar(0,0) < 0)
        help = 0;
      else
        help = sqrt(FCsigma2.betavar(0,0));
      outscale << help << "  ";

      outscale << FCsigma2.betaqu_l1_lower(0,0) << "  ";

      outscale << FCsigma2.betaqu_l2_lower(0,0) << "  ";

      outscale << FCsigma2.betaqu50(0,0) << "  ";

      outscale << FCsigma2.betaqu_l2_upper(0,0) << "  ";

      outscale << FCsigma2.betaqu_l1_upper(0,0) << "  ";
      }

    outscale << endl;
    }

  optionsp->out("\n");

  }


double DISTR_gaussian::get_scalemean(void)
  {
  return FCsigma2.betamean(0,0);
  }


//------------------------------------------------------------------------------
//----------------------------- DISTR_vargaussian ------------------------------
//------------------------------------------------------------------------------

//  void check_errors(void);

  // CONSTRUCTOR1
  // TASK: initializes data
  //       response = r
  //       weight = w
  //       nrobs = r.rows()

  DISTR_vargaussian::DISTR_vargaussian(GENERAL_OPTIONS * o,const datamatrix & r)
   : DISTR(o,r)
    {

    predictor_name = "sigma2";
    maindistribution=false;
    family="heteroscedastic Gaussian, variance component";
    wtype = wweightschange_weightsneqone;
    weightsone=false;
    updateIWLS = true;
    sigma2old=0;

    linpredminlimit=-10;
    linpredmaxlimit= 15;

    }

  // COPY CONSTRUCTOR

  DISTR_vargaussian::DISTR_vargaussian(const DISTR_vargaussian & d)
     : DISTR(DISTR(d))
    {
    dgaussian = d.dgaussian;
    sigma2old=d.sigma2old;
    }

  // OVERLOADED ASSIGNMENT OPERATOR

  const DISTR_vargaussian & DISTR_vargaussian::operator=(
                                               const DISTR_vargaussian & d)
    {
    if (this==&d)
      return *this;
    DISTR::operator=(DISTR(d));
    dgaussian = d.dgaussian;
    sigma2old=d.sigma2old;
    return *this;
    }


  double DISTR_vargaussian::loglikelihood(double * res,double * lin,
                                          double * weight)
    {
    if (*weight !=0)
      {
      double m = exp(*lin);
      return  -  (*res)/(2*m) - 0.5* (*lin) ;
      }
    else
      return 0;
    }


void DISTR_vargaussian::compute_mu(const double * linpred,double * mu)
  {
  *mu = exp(*linpred);
  }


// double compute_MSE(const double * response, const double * weight,
//                             const double * linpred, msetype t, double v);


  void DISTR_vargaussian::outoptions(void)
    {
    DISTR::outoptions();

    optionsp->out("\n");
    optionsp->out("\n");

    }


  double DISTR_vargaussian::compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse,const bool & like)
    {
    double m = exp(*linpred);

    *workingweight=(*weight)*0.5;

    *workingresponse = *linpred + (*response - m)/m;

    if (like && (*weight != 0))
      return  -  (*response)/(2*m) - 0.5* (*linpred);
    else
      return 0;
    }


bool DISTR_vargaussian::posteriormode(void)
  {

  dgaussian->FCpredict_betamean_vargaussian = FCpredict_betamean;
  weight = dgaussian->weightoriginal;

  double s = dgaussian->sigma2;

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\linpredvar.res");
  // linearpred1.prettyPrint(out);
  // END: TEST

  double * worklin_mean;
  if (dgaussian->linpred_current==1)
    worklin_mean = dgaussian->linearpred1.getV();
  else
    worklin_mean = dgaussian->linearpred2.getV();

  double * wweight_mean = dgaussian->weight.getV();

  double * wweight_orig = dgaussian->weightoriginal.getV();

  double * workresponse = response.getV();
  double * workresponse_mean = dgaussian->response.getV();

  unsigned i;
  double m;
  if (sigma2old==0)
    {
    for (i=0;i<nrobs;i++, wweight_mean++, worklin++,wweight_orig++,
    worklin_mean++,workresponse++,workresponse_mean++)
      {
      *workresponse = pow(*workresponse_mean-(*worklin_mean),2);
      if (*wweight_orig==0)
        *worklin = log(s);
      else
        *worklin = log(s/(*wweight_orig));
      }
    }
    else
      {
      for (i=0;i<nrobs;i++, wweight_mean++, worklin++,wweight_orig++,
           worklin_mean++,workresponse++,workresponse_mean++)
        {
        *workresponse = pow(*workresponse_mean-(*worklin_mean),2);
        if (*wweight_orig==0)
          {
          *worklin -= log(sigma2old);
          *wweight_mean = 0;
          *worklin += log(s);
          }
        else
          {
          *worklin -= log(sigma2old/(*wweight_orig));
          m = exp(*worklin);
          *wweight_mean = 1/m;
          *worklin += log(s/(*wweight_orig));
          }
        }
      }

    sigma2old = s;

    return true;
    }


void DISTR_vargaussian::update(void)
  {

  dgaussian->FCpredict_betamean_vargaussian = FCpredict_betamean;

  unsigned i;
  double * workresponse = response.getV();
  double * workresponse_mean = dgaussian->response.getV();

  double * worklin_mean;
  if (dgaussian->linpred_current==1)
    worklin_mean = dgaussian->linearpred1.getV();
  else
    worklin_mean = dgaussian->linearpred2.getV();


  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * wweight_mean = dgaussian->weight.getV();

  double * wweight_orig = dgaussian->weightoriginal.getV();

  double s = dgaussian->sigma2;
  double m;

  for (i=0;i<nrobs;i++,workresponse++,workresponse_mean++,worklin_mean++,
        wweight_mean++,wweight_orig++,worklin++)
    {
    *workresponse = pow(*workresponse_mean-(*worklin_mean),2);
    if (*wweight_orig==0)
      {
      *worklin-=log(sigma2old);
      *wweight_mean = 0;
      *worklin+=log(s);
      }
    else
      {
      *worklin-=log(sigma2old/(*wweight_orig));
      m = exp(*worklin);
      *wweight_mean = 1/m;
      *worklin+=log(s/(*wweight_orig));
      }


    }

    sigma2old=s;

    DISTR::update();

    }


//------------------------------------------------------------------------------
//-------------------- CLASS DISTRIBUTION_hetgaussian --------------------------
//------------------------------------------------------------------------------

DISTR_hetgaussian::DISTR_hetgaussian(double a,double b, GENERAL_OPTIONS * o,
                                     const datamatrix & r,
                                     const ST::string & ps, const bool sc,
                                     const datamatrix & w)
  : DISTR_gaussian(a,b,o,r,ps,w)

  {

  predictor_name = "mu";
  outexpectation = true;

  wtype = wweightschange_weightsneqone;

  family = "Heteroscedastic Gaussian";

  weightoriginal = weight;

  sigma2const=sc;
  sigma2=1;

  }


const DISTR_hetgaussian & DISTR_hetgaussian::operator=(
                                      const DISTR_hetgaussian & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gaussian::operator=(DISTR_gaussian(nd));
  weightoriginal = nd.weightoriginal;
  FCpredict_betamean_vargaussian = nd.FCpredict_betamean_vargaussian;
  sigma2const = nd.sigma2const;
  return *this;
  }


DISTR_hetgaussian::DISTR_hetgaussian(const DISTR_hetgaussian & nd)
   : DISTR_gaussian(DISTR_gaussian(nd))
  {
  weightoriginal = nd.weightoriginal;
  FCpredict_betamean_vargaussian = nd.FCpredict_betamean_vargaussian;
  sigma2const = nd.sigma2const;
  }


void DISTR_hetgaussian::compute_MSE_all(datamatrix & meanpred, double & MSE,
                               double & MSEzeroweight, unsigned & nrzeroweights,
                               msetype & t, double & v)
  {
  unsigned i;
  nrzeroweights = 0;
  MSE = 0;
  MSEzeroweight=0;
  double * responsep = response.getV();
  double * weightp = FCpredict_betamean_vargaussian->getV();
  weightp++;
  double * weightorigp = weightoriginal.getV();

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\betamean.raw");
  // FCpredict_betamean_vargaussian->prettyPrint(out);
  // TEST

  double w;

  double * linpredp = meanpred.getV();
  for(i=0;i<nrobs;i++,responsep++,linpredp+=2,weightp+=2,
      weightorigp++)
    {
    w=1/(*weightp);
    if (*weightorigp==0)
      {
      MSEzeroweight += compute_MSE(responsep,&w,linpredp,t,v);
      nrzeroweights++;
      }
    else
      MSE += compute_MSE(responsep,&w,linpredp,t,v);
    }
  }


double DISTR_hetgaussian::compute_MSE(const double * response,
                          const double * weight,
                          const double * linpred, msetype t, double v)
  {
  if (t == quadraticMSE)
    return pow(*response-*linpred,2);
  else
    {
    double u;

    u = *response - ( *linpred +
    sqrt(1/(*weight)) * randnumbers::invPhi2(v) );

    if (u >= 0)
      return u*v;
    else
      return u*(v-1);
    }
  }


void DISTR_hetgaussian::update(void)
  {
  if (!sigma2const)
    DISTR_gaussian::update();
  }


bool DISTR_hetgaussian::posteriormode(void)
  {
  if (!sigma2const)
    return DISTR_gaussian::posteriormode();
  else
    return true;
  }


void DISTR_hetgaussian::outresults(ofstream & out_stata, ofstream & out_R,
                                   ST::string pathresults)
  {
  if (!sigma2const)
    DISTR_gaussian::outresults(out_stata,out_R,pathresults);
  }

//------------------------------------------------------------------------------
//------------------------- CLASS DISTR_quantreg -------------------------------
//------------------------------------------------------------------------------

DISTR_quantreg::DISTR_quantreg(const double & a,const double & b,
                               GENERAL_OPTIONS * o, const datamatrix & r,
                               const ST::string & ps,double & quant,
                               const datamatrix & w)
  : DISTR_gaussian(a,b,o,r,ps,w)
  {
  family="Quantile regression based on asymmetric Laplace distribution";
  predictor_name = "quantile";
  outexpectation = true;

  quantile = quant;
  xi = (1-2*quantile)/(quantile * (1-quantile));
  xi2 = xi*xi;
  sigma02 = 2 / (quantile*(1-quantile));
  num = sqrt(xi2 + 2*sigma02);
  wtype = wweightschange_weightsneqone;
  }


DISTR_quantreg::DISTR_quantreg(const DISTR_quantreg & nd)
  : DISTR_gaussian(DISTR_gaussian(nd))
  {
  quantile = nd.quantile;
  xi = nd.xi;
  xi2 = nd.xi2;
  num = nd.num;
  sigma02 = nd.sigma02;
  }


const DISTR_quantreg & DISTR_quantreg::operator=(
                                      const DISTR_quantreg & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gaussian::operator=(DISTR_gaussian(nd));
  quantile = nd.quantile;
  xi = nd.xi;
  xi2 = nd.xi2;
  num = nd.num;
  sigma02 = nd.sigma02;
  return *this;
  }


/*
void DISTR_quantreg::compute_mu(const double * linpred,double * mu,
                                bool notransform)
  {

  }
*/


/*
void DISTR_quantreg::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           double * scale) const
{
}
*/

/*
double DISTR_quantreg::loglikelihood(double * res, double * lin,
                                     double * w) const
  {

  }
*/

/*
double DISTR_quantreg::loglikelihood_weightsone(double * res,double * lin) const
  {

  }


double DISTR_quantreg::compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like)
  {

  }


void DISTR_quantreg::compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like)
  {

  }


void DISTR_quantreg::compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  }


void DISTR_quantreg::compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like)
  {
  }
*/


double DISTR_quantreg::compute_MSE(const double * response,
                                   const double * weight,
                                   const double * linpred, msetype t,double v)
  {
  if (t==quadraticMSE)
    return pow(*response-*linpred,2);
  else
    {
    double u = *response-*linpred;
    if (u >= 0)
      return u*quantile;
    else
      return u*(quantile-1);
    }
  }



void DISTR_quantreg::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: identity\n");
  optionsp->out("  Hyperparameter a: " + ST::doubletostring(a_invgamma,6) + "\n");
  optionsp->out("  Hyperparameter b: " + ST::doubletostring(b_invgamma,6) + "\n");
  optionsp->out("  Quantile: " + ST::doubletostring(quantile,6) + "\n");
  optionsp->out("\n");
  optionsp->out("\n");

  }


void DISTR_quantreg::update(void)
  {

  // sigma2 corresponds to 1/delta_0
  // sigma02 corresponds to sigma^2 in the paper

  unsigned i;

  double * worklin;
  double * workresp;
  double * workorigresp;
  double * workweight;
  double * workw;

  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  workresp = workingresponse.getV();
  workorigresp = response.getV();
  workweight = workingweight.getV();
  workw = weight.getV();

  double lambda = (xi2 + 2*sigma02) / sigma2;

  double mu;

  double sumres=0;
  double sumw=0;


  //  sigma02 = delta2

  for(i=0; i<nrobs; i++, workw++, workweight++, workresp++, workorigresp++,
      worklin++)
    {

    mu = num/(fabs(*workorigresp-*worklin));

    if (*workw == 0)
      {
      *workweight=0;
      *workresp = 0;
      }
    else
      {
      *workweight = rand_inv_gaussian(mu, lambda);         // 1/z_i

      *workresp = *workorigresp  - xi/ (*workweight);      // y_i - offset = y_i - xi * z_i
      sumw += 1 / (*workweight);                           // sum z_i
      sumres += *workweight * pow(*workresp-*worklin,2);   // 1/z_i * (y_i-eta_i - xi * z_i)
      }
    }

  sumres /= sigma02;                  // 1/(delta2 *z_i) * (y_i-eta_i - xi * z_i)

  sigma2 = rand_invgamma(a_invgamma + 1.5*(nrobs-nrzeroweights),
                         b_invgamma + 0.5*sumres + sumw);

//  DISTR_gaussian::update();  // ?????

  FCsigma2.beta(0,0) = sigma2;
  FCsigma2.acceptance++;
  FCsigma2.update();

  sigma2 *= sigma02;

  }




/*
bool DISTR_quantreg::posteriormode(void)
  {

  }


void DISTR_quantreg::outresults(ST::string pathresults="")
  {

  }

double DISTR_quantreg::get_scalemean(void)
  {

  }


void DISTR_quantreg::sample_responses(unsigned i,datamatrix & sr)
  {

  }


void DISTR_quantregoutresults_predictive_check(datamatrix & D,datamatrix & sr)
  {

  }
*/

//------------------------------------------------------------------------------
//---------------------- CLASS DISTRIBUTION_loggaussian ------------------------
//------------------------------------------------------------------------------


DISTR_loggaussian::DISTR_loggaussian(const double & a,
                                             const double & b,
                                             GENERAL_OPTIONS * o,
                                             const datamatrix & r,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTR_gaussian(a,b,o,r,ps,w)

  {
  family = "log-Gaussian";
  outexpectation = true;
  }


const DISTR_loggaussian & DISTR_loggaussian::operator=(
                                      const DISTR_loggaussian & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gaussian::operator=(DISTR_gaussian(nd));
  return *this;
  }


DISTR_loggaussian::DISTR_loggaussian(const DISTR_loggaussian & nd)
   : DISTR_gaussian(DISTR_gaussian(nd))
  {
  }


void DISTR_loggaussian::compute_mu(const double * linpred,double * mu)
  {
  *mu = exp((*linpred) + 0.5*sigma2);
  }


double DISTR_loggaussian::compute_MSE(const double * response,
                                      const double * weight,
                                      const double * linpred,msetype t, double v)
  {
  double diff = exp(*response) - exp(*linpred +  FCsigma2.betamean(0,0)/2);
  if (t == quadraticMSE)
    return pow(diff,2);
  else
    {
    if (diff >= 0)
      return diff*v;
    else
      return diff*(v-1);
    }
  }




void DISTR_loggaussian::compute_deviance(const double * response,
                                 const double * weight, const double * mu,
                                 double * deviance,
                                 double * scale) const
  {


  double pred = log(*mu)-(*scale)/(2*(*weight));
  double deviancehelp;

  if (*weight != 0)
    {
    DISTR_gaussian::compute_deviance(response,weight,&pred,&deviancehelp,scale);
    *deviance =  2*(*response)  + deviancehelp;
    }

  else
    {
    *deviance = 0;
    }

  }


void DISTR_loggaussian::sample_responses(unsigned i,datamatrix & sr)
  {

  double * linpredp;

  if (linpred_current==1)
    linpredp = linearpred1.getV();
  else
    linpredp = linearpred2.getV();

  double * rp = sr.getV()+i;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
    *rp = exp(*linpredp+sqrt(sigma2)*rand_normal());
  }


void DISTR_loggaussian::sample_responses_cv(unsigned i,datamatrix & linpred,
                                            datamatrix & sr)
  {

  double * linpredp;

  linpredp = linpred.getV();

  double * rp = sr.getV()+i;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
    *rp = exp(*linpredp+sqrt(sigma2)*rand_normal());
  }



void DISTR_loggaussian::outresults_predictive_check(datamatrix & D,datamatrix & sr)
  {

  datamatrix Dh = D;
  unsigned i;
  for (i=0;i<D.rows();i++)
    Dh(i,0) = exp(D(i,0));

  DISTR::outresults_predictive_check(Dh,sr);
  }

//------------------------------------------------------------------------------
//--------------------- CLASS DISTRIBUTION_gaussian_exp ------------------------
//------------------------------------------------------------------------------


DISTR_gaussian_exp::DISTR_gaussian_exp(const double & a,
                                             const double & b,
                                             GENERAL_OPTIONS * o,
                                             const datamatrix & r,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTR_gaussian(a,b,o,r,ps,w)

  {
  updateIWLS = true;
  outexpectation = true;
  }



const DISTR_gaussian_exp & DISTR_gaussian_exp::operator=(
                                      const DISTR_gaussian_exp & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gaussian::operator=(DISTR_gaussian(nd));
  return *this;
  }


DISTR_gaussian_exp::DISTR_gaussian_exp(const DISTR_gaussian_exp & nd)
   : DISTR_gaussian(DISTR_gaussian(nd))
  {
  }


void DISTR_gaussian_exp::outoptions(void)
  {
  DISTR::outoptions();
  optionsp->out("  Response function: exponential\n");

  optionsp->out("  Hyperparameter a: " + ST::doubletostring(a_invgamma,6) + "\n");
  optionsp->out("  Hyperparameter b: " + ST::doubletostring(b_invgamma,6) + "\n");

  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTR_gaussian_exp::update(void)
  {

  register unsigned i;

  double help;

  double * worklin;
  double * workresp;
  double * workweight;


  // scaleparameter

  double sum = 0;

  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  workresp = response.getV();
  workweight = weight.getV();

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - exp(*worklin);
    sum += *workweight*pow(help,2);
    }

  sigma2  = rand_invgamma(a_invgamma+0.5*nrobs,
                          b_invgamma+0.5*sum);

  FCsigma2.beta(0,0) = sigma2;
  FCsigma2.acceptance++;
  FCsigma2.update();

  DISTR::update();

  }


void DISTR_gaussian_exp::compute_mu(const double * linpred,double * mu)
  {
  *mu = exp(*linpred);
  }


void DISTR_gaussian_exp::compute_param(const double * linpred,double * param)
  {
  *param = (*linpred);
  }


double DISTR_gaussian_exp::loglikelihood(double * res, double * lin,
                                         double * w)
  {
  double help = *res-exp(*lin);
  return  - *w * (pow(help,2))/(2* sigma2);
  }


double DISTR_gaussian_exp::compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like)
  {
  double mu = exp(*linpred);
  double mu2 = pow(mu,2);
  *workingweight=mu2 * (*weight)/sigma2;

  *workingresponse = (*response-mu)/mu + (*linpred);

  if (like==true)
    {
    double h = *response-mu;
    return  - (*weight) * pow(h,2)/(2* sigma2);
    }
  else
    return 0;
  }


bool DISTR_gaussian_exp::posteriormode(void)
  {

  unsigned i;

  double * worklin;
  if (linpred_current==1)
    worklin = linearpred1.getV();
  else
    worklin = linearpred2.getV();

  double * workresp = response.getV();
  double * workweight = weight.getV();

  double sum = 0;
  double sumweight=0;
  double help;

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - exp(*worklin);
    sum += *workweight*pow(help,2);
    sumweight+=*workweight;
    }

  sigma2 = (1.0/sumweight)*sum;

  FCsigma2.beta(0,0) = sigma2;

  FCsigma2.posteriormode_betamean();

  return true;

  }


void DISTR_gaussian_exp::sample_responses(unsigned i,datamatrix & sr)
  {

  double * linpredp;

  if (linpred_current==1)
    linpredp = linearpred1.getV();
  else
    linpredp = linearpred2.getV();

  double * rp = sr.getV()+i;

  unsigned j;
  for (j=0;j<nrobs;j++,linpredp++,rp+=sr.cols())
//    *rp = trmult*( exp(*linpredp)+sqrt(sigma2)*rand_normal() );
    *rp =  exp(*linpredp)+sqrt(sigma2)*rand_normal();
  }


void DISTR_gaussian_exp::outresults_predictive_check(datamatrix & D,
                                                     datamatrix & sr)
  {
  DISTR::outresults_predictive_check(D,sr);
  }


//------------------------------------------------------------------------------
//--------------------- CLASS DISTRIBUTION_gaussian_mult -----------------------
//------------------------------------------------------------------------------


DISTR_gaussian_mult::DISTR_gaussian_mult(const double & a,
                                             const double & b,
                                             GENERAL_OPTIONS * o,
                                             const datamatrix & r,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTR_gaussian_exp(a,b,o,r,ps,w)

  {
  //standardise();
  optionbool1 = false;
//  changingworkingweights = false;
  updateIWLS = false;
  }


void DISTR_gaussian_mult::set_mult(bool & m)
  {

  if (m==true)
    {
    optionbool1 = true;
//    changingworkingweights = true;
    updateIWLS = true;

    }
  else
    {
    optionbool1 = false;
//    changingworkingweights = false;
    updateIWLS = false;
    }

  }


const DISTR_gaussian_mult & DISTR_gaussian_mult::operator=(
                                      const DISTR_gaussian_mult & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gaussian_exp::operator=(DISTR_gaussian_exp(nd));
  return *this;
  }


DISTR_gaussian_mult::DISTR_gaussian_mult(const DISTR_gaussian_mult & nd)
   : DISTR_gaussian_exp(DISTR_gaussian_exp(nd))
  {
  }


/*
void DISTR_gaussian_mult::standardise(void)
  {

  trmult = 1;

  unsigned i;
  double * workresp = workingresponse.getV();
  double * worklin = linearpred1.getV();
  for (i=0;i<nrobs;i++,workresp++,worklin++)
    {
    *workresp = response(i,0);
    *worklin = 0;
    }

  FCsigma2.transform(0,0) = pow(trmult,2);

  }
*/


void DISTR_gaussian_mult::outoptions(void)
  {
  DISTR_gaussian::outoptions();
  }


void DISTR_gaussian_mult::update(void)
  {
  if (optionbool1==true)
    DISTR_gaussian_exp::update();
  else
    DISTR_gaussian::update();
  }

/*
void DISTR_gaussian_mult::compute_mu(const double * linpred,double * mu,
                                    bool notransform)
  {
  if (!optionbool1)
    {
    DISTR_gaussian::compute_mu(linpred,mu,notransform);
    }
 else
   {
   DISTR_gaussian_exp::compute_mu(linpred,mu,notransform);
   }

  }
*/

void DISTR_gaussian_mult::compute_mu(const double * linpred,double * mu)
  {
  if (!optionbool1)
    {
    DISTR_gaussian::compute_mu(linpred,mu);
    }
 else
   {
   DISTR_gaussian_exp::compute_mu(linpred,mu);
   }

  }


double DISTR_gaussian_mult::loglikelihood(double * res, double * lin,
                                         double * w)
  {

  if (!optionbool1)
    {
    return DISTR_gaussian::loglikelihood(res,lin,w);
    }
 else
   {
   return DISTR_gaussian_exp::loglikelihood(res,lin,w);
   }

  }


double DISTR_gaussian_mult::compute_iwls(double * response, double * linpred,
                              double * weight, double * workingweight,
                              double * workingresponse, const bool & like)
  {


  if (!optionbool1)
    {
    return DISTR_gaussian::compute_iwls(response,linpred,weight,workingweight,
                                         workingresponse,like);
    }
 else
   {
    return DISTR_gaussian_exp::compute_iwls(response,linpred,weight,workingweight,
                                         workingresponse,like);
   }

  }


bool DISTR_gaussian_mult::posteriormode(void)
  {

  if (!optionbool1)
    {
    return DISTR_gaussian::posteriormode();

    }
 else
   {
    return DISTR_gaussian_exp::posteriormode();
   }


  }


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_gaussian_re -----------------------
//------------------------------------------------------------------------------


DISTR_gaussian_re::DISTR_gaussian_re(GENERAL_OPTIONS * o,const datamatrix & r,
                                     const datamatrix & w)
  : DISTR_gaussian(1,1,o,r,"",w)

  {
  maindistribution=false;

  family = "Gaussian_random_effect";

  check_errors();
  }



void DISTR_gaussian_re::check_errors(void)
  {
  DISTR::check_errors();
  bool helperror;


  unsigned col = 0;
  helperror = response.check_ascending(col);
  if (helperror == false)
    {
    errors = true;
    errormessages.push_back("ERROR: group indicator values must be sorted in ascending order for distribution gaussian_re\n");
    }


  if (!weightsone)
    {
    errors = true;
    errormessages.push_back("ERROR: weights not allowed for distribution gaussian_re\n");
    }

  }


const DISTR_gaussian_re & DISTR_gaussian_re::operator=(
                                      const DISTR_gaussian_re & nd)
  {
  if (this==&nd)
    return *this;
  DISTR_gaussian::operator=(DISTR_gaussian(nd));
  return *this;
  }



DISTR_gaussian_re::DISTR_gaussian_re(const DISTR_gaussian_re & nd)
   : DISTR_gaussian(DISTR_gaussian(nd))
  {
  }



void DISTR_gaussian_re::update(void)
  {

  // TEST
  //  ofstream out("c:\\bayesx\\test\\results\\response_RE.res");
  //  response.prettyPrint(out);
  // ENDE TEST

  }



bool DISTR_gaussian_re::posteriormode(void)
  {

  // TEST
  /*
    ofstream out("c:\\bayesx\\test\\results\\response_RE.res");
    response.prettyPrint(out);
  */
  // ENDE TEST

  return true;
  }

void DISTR_gaussian_re::outoptions(void)
  {
  optionsp->out("RANDOM EFFECTS DISTRIBUTION:\n",true);
  optionsp->out("\n");
  optionsp->out("  Family: " + family + "\n");
  optionsp->out("  Number of clusters: " + ST::inttostring(nrobs) + "\n");
  optionsp->out("\n");
  }


void DISTR_gaussian_re::outresults(ofstream & out_stata, ofstream & out_R,
                                   ST::string pathresults)
  {
  }


void DISTR_gaussian_re::get_samples(const ST::string &
                                    filename,ofstream & outg) const
  {
  }


} // end: namespace MCMC
