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



#include "FC_nonp.h"


//------------------------------------------------------------------------------
//----------------- CLASS: FC implementation of member functions ---------------
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_nonp::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  7       center
  8       map
  9       lambda_re
  10      a_re
  11      b_re
  12      internal_mult
  13      samplemult
  14      constraints
  15      round
  16      centermethod
  17      internal_mult
  18      pvalue
  19      meaneffect
  20      binning
  21      update
  */

  if (op[14] == "increasing")
    stype = increasing;
  else if (op[14] == "decreasing")
    stype = decreasing;
  else
    stype = unconstrained;

  if (op[19] == "true")
    computemeaneffect = true;
  else
    computemeaneffect = false;

  if (op[21] == "direct")
    orthogonal = false;
  else
    orthogonal = true;

  int f;

  f = op[25].strtodouble(s2);


  if (op[26] == "false")
     derivative = false;
  else
    derivative = true;

  if (op[27] == "false")
     samplederivative = false;
  else
    samplederivative = true;

  if (op[28] == "false")
     samplef = false;
  else
    samplef = true;

   f = op[34].strtodouble(meaneffectconstant);

  }


void FC_nonp::check_errors(void)
  {
  FC::check_errors();
  }


FC_nonp::FC_nonp(void)
  {
  }


FC_nonp::FC_nonp(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,
                 const ST::string & t,const ST::string & fp,
                 DESIGN * Dp,vector<ST::string> & op,
                 vector<ST::string> & vn)
     : FC(o,t,Dp->Zout.rows(),1,fp)
  {

  read_options(op,vn);

  imeasures=false;

  masterp = mp;
  equationnr = enr;
  likep = lp;
  designp = Dp;

  if (designp->errors==false)
    {
    param = datamatrix(designp->nrpar,1,0);
    paramold = param;
    parammode = param;
    paramhelp = param;
    betaold = beta;
    betadiff = beta;
    partres = datamatrix(designp->posbeg.size(),1,0);
    lambda=1;
    tau2 = likep->get_scale()/lambda;
    IWLS = likep->updateIWLS;

    if (Dp->position_lin != -1)
      {
      fsample = FC(o,"",beta.rows(),beta.cols(),fp + ".lin");
      paramlin = datamatrix(Dp->designlinear.cols(),1,0);
      }

    paramsample = FC(o,"",param.rows(),1,fp + ".param");

    helpcenter = datamatrix(designp->nrpar,1);

    if (computemeaneffect==true)
      {
      meaneffect_sample = FC(o,"",beta.rows(),1,fp+".meaneffect");
      }

    if (derivative == true)
      {
      if (designp->type==Mrf)
        derivative = false;
      else
        {
        derivativesample = FC(o,"",beta.rows(),1,fp+".derivative");
        }
      }

    } // end: if designp->errors() == false
  else
    {
    errors=true;
    }

//  check_errors();

  }


FC_nonp::FC_nonp(const FC_nonp & m)
  : FC(FC(m))
  {

  imeasures=m.imeasures;

  masterp = m.masterp;
  equationnr = m.equationnr;

  fsample = m.fsample;

  paramsample = m.paramsample;

  computemeaneffect = m.computemeaneffect;
  meaneffectconstant = m.meaneffectconstant;
  meaneffect_sample = m.meaneffect_sample;

  stype = m.stype;
  likep = m.likep;
  designp = m.designp;
  param = m.param;
  paramlin = m.paramlin;
  parammode = m.parammode;
  paramold = m.paramold;
  paramhelp = m.paramhelp;
  paramKparam = paramKparam;
  betaold = m.betaold;
  betadiff = m.betadiff;
  partres = m.partres;
  lambda=m.lambda;
  tau2 = m.tau2;
  IWLS = m.IWLS;
  orthogonal = m.orthogonal;

  Vcenter = m.Vcenter;
  Vcentert = m.Vcentert;
  Wcenter = m.Wcenter;
  Ucenter = m.Ucenter;
  Utc = m.Utc;
  ccenter = m.ccenter;
  helpcenter = m.helpcenter;

  acuteparam = m.acuteparam;
  s2 = m.s2;

  derivative = m.derivative;
  derivativesample = m.derivativesample;
  samplederivative = m.samplederivative;
  samplef = m.samplef;
  }


const FC_nonp & FC_nonp::operator=(const FC_nonp & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));

  imeasures=m.imeasures;

  masterp = m.masterp;
  equationnr = m.equationnr;

  fsample = m.fsample;

  paramsample = m.paramsample;

  computemeaneffect = m.computemeaneffect;
  meaneffectconstant = m.meaneffectconstant;
  meaneffect_sample = m.meaneffect_sample;

  stype = m.stype;
  likep = m.likep;
  designp = m.designp;
  param = m.param;
  paramlin = m.paramlin;
  paramKparam = paramKparam;
  parammode = m.parammode;
  paramold = m.paramold;
  paramhelp = m.paramhelp;
  betaold = m.betaold;
  betadiff = m.betadiff;
  partres = m.partres;
  lambda=m.lambda;
  tau2 = m.tau2;
  IWLS = m.IWLS;
  orthogonal = m.orthogonal;

  Vcenter = m.Vcenter;
  Vcentert = m.Vcentert;
  Wcenter = m.Wcenter;
  Ucenter = m.Ucenter;
  Utc = m.Utc;
  ccenter = m.ccenter;
  helpcenter = m.helpcenter;

  acuteparam = m.acuteparam;
  s2 = m.s2;

  derivative = m.derivative;
  derivativesample = m.derivativesample;
  samplederivative = m.samplederivative;
  samplef = m.samplef;

  return *this;
  }


void FC_nonp::get_linparam(void)
  {
  int pos = designp->position_lin;
  int i;
  for (i=0;i<paramlin.rows();i++)
    paramlin(i,0) = designp->FClinearp->beta(pos+i,0);
  }


void FC_nonp::perform_centering(void)
  {

  if(designp->center)
    {
    if (designp->centermethod==meansimple)
      centerparam();
    else if (designp->centermethod==integralsimple)
      centerparam_weight();
    else if (designp->centermethod==meansum2)
      centerparam_sum2(s2);
    else if (designp->centermethod==meansimplevar)
      centerparam_sample_var();
    else
      centerparam_sample();
    }

  if (designp->position_lin!=-1)
    {
    get_linparam();
    }

  }


void FC_nonp::update_IWLS(void)
  {

  // TEST
//  ofstream out("c:\\bayesx\\testh\\results\\beta.res");
//  beta.prettyPrint(out);
  // TEST


  unsigned i;
  double * workparam;

//  lambda = likep->get_scale()/tau2;
  lambda = 1/tau2;

  if (optionsp->nriter == 1)
    {
    paramold.assign(param);
    betaold.assign(beta);
    paramKparam=designp->penalty_compute_quadform(param);
    }


  // Compute log-likelihood with old param, computes workingweight and
  // workingresponse
  double logold = likep->compute_iwls(true,true);
  logold -= 0.5*paramKparam*lambda;

  designp->compute_partres(partres,beta);
  designp->compute_XtransposedWX();
  designp->compute_XtransposedWres(partres,lambda);

  designp->compute_precision(lambda);

//  bool error = designp->precision.decomp_save();

//  if (error == false)
//    {
    designp->precision.solve(*(designp->XWres_p),paramhelp);

    // TEST
    // ofstream out("c:\\bayesx\\testh\\results\\paramhelp_v.res");
    // paramhelp.prettyPrint(out);
    // TEST

    workparam = param.getV();
    unsigned nrpar = param.rows();
    for(i=0;i<nrpar;i++,workparam++)
      *workparam = rand_normal();

    designp->precision.solveU(param,paramhelp); // param contains now the proposed
                                                // new parametervector

    perform_centering();

    paramhelp.minus(param,paramhelp);

    double qold = 0.5*designp->precision.getLogDet()-
                0.5*designp->precision.compute_quadform(paramhelp,0);

    designp->compute_f(param,paramlin,beta,fsample.beta);

    betadiff.minus(beta,betaold);

    bool ok;
    if (optionsp->saveestimation)
      {
      ok = designp->update_linpred_save(betadiff);
      if (!ok)
        outsidelinpredlimits++;
      }
    else
      {
      designp->update_linpred(betadiff);
      ok = true;
      }

    // Compute new log-likelihood

    double lognew=0;
    double qnew=0;
    if (ok)
      {
      lognew   = likep->compute_iwls(true,true);
      lognew  -= 0.5*designp->penalty_compute_quadform(param)*lambda;

      designp->compute_partres(partres,beta);
      designp->compute_XtransposedWX();
      designp->compute_XtransposedWres(partres,lambda);

      designp->compute_precision(lambda);

      designp->precision.solve(*(designp->XWres_p),paramhelp);

      // TEST
      // ofstream out2("c:\\bayesx\\testh\\results\\paramhelp_n.res");
      // paramhelp.prettyPrint(out2);
      // TEST


      paramhelp.minus(paramold,paramhelp);
      qnew = 0.5*designp->precision.getLogDet() -
             0.5*designp->precision.compute_quadform(paramhelp,0);
      }

    double u = log(uniform());

    if (ok && (u <= (lognew - logold  + qnew - qold)) )
      {
      acceptance++;

      paramKparam=designp->penalty_compute_quadform(param);

      betaold.assign(beta);
      paramold.assign(param);
      }
    else
      {

      betadiff.minus(betaold,beta);
      designp->update_linpred(betadiff);


      param.assign(paramold);
      beta.assign(betaold);
      }

//    } // end if (error==false)

  if (derivative)
    {
    designp->compute_f_derivative(param,paramlin,derivativesample.beta,
                                    derivativesample.beta);
    }


  // TEST

  // ofstream out("c:\\bayesx\\test\\results\\param.res");
  // param.prettyPrint(out);

  // ofstream out2("c:\\bayesx\\test\\results\\beta.res");
  // beta.prettyPrint(out2);

  // TEST

  // transform_beta();

  if (designp->position_lin!=-1)
    {
    fsample.update();
    }


  paramsample.beta.assign(param);
  paramsample.update();

  if (derivative)
    derivativesample.update();

  FC::update();

  }


void FC_nonp::update(void)
  {
  if (IWLS)
    {
    update_IWLS();
    }
  else
    {
    if ((stype == increasing) || (stype==decreasing))
      {
      update_isotonic();
      }
    else
      {
      update_gaussian();
      }
    }


  designp->compute_meaneffect(masterp->level1_likep[equationnr],meaneffect,beta,
                             meaneffect_sample.beta,computemeaneffect,
                             meaneffectconstant);

  if (computemeaneffect == true)
    {
    meaneffect_sample.update();
    }

  }


void FC_nonp::update_gaussian_transform(void)
  {

  betaold.assign(beta);
  if (optionsp->saveestimation)
    paramold.assign(param);

  double sigma2resp = likep->get_scale();
  lambda = sigma2resp/tau2;

  designp->compute_partres(partres,beta);

  designp->compute_XtransposedWres(partres, lambda);

  designp->u.mult(designp->QtRinv,*(designp->XWres_p));

  unsigned j;
  double * up = designp->u.getV();
  double * sp = designp->s.getV();
  double * acuteparamp = acuteparam.getV();
  double mu,var;
  double h;
  for (j=0;j<param.rows();j++,up++,sp++,acuteparamp++)
    {
    h = 1/(1+lambda* (*sp));
    mu = h * (*up);
    var = sigma2resp * h;
    *acuteparamp = mu + sqrt(var)*rand_normal();
    }

  param.mult(designp->RtinvQ,acuteparam);

  perform_centering();

  designp->compute_f(param,paramlin,beta,fsample.beta);

  betadiff.minus(beta,betaold);

  bool ok;
  if (optionsp->saveestimation)
    {
    ok = designp->update_linpred_save(betadiff);

    if (!ok)
      {
      outsidelinpredlimits++;
      betadiff.minus(betaold,beta);
      designp->update_linpred(betadiff);
      beta.assign(betaold);
      param.assign(paramold);
      }
    else
      acceptance++;
    }
  else
    {
    designp->update_linpred(betadiff);
    ok = true;
    acceptance++;
    }

/*
  if (designp->position_lin!=-1)
    {
    fsample.update();
    }
*/

  paramsample.beta.assign(param);
  paramsample.update();

  FC::update();

  }


void FC_nonp::update_gaussian(void)
  {

  // TEST
  // ofstream out0("c:\\bayesx\\testh\\results\\beta_re.res");
  // beta.prettyPrint(out0);

  // ofstream out00("c:\\bayesx\\testh\\results\\intvar_re.res");
  // designp->intvar.prettyPrint(out00);


  // TEST

  if (orthogonal)
    update_gaussian_transform();
  else
    {
    bool lambdaconst = false;

    betaold.assign(beta);
    if (optionsp->saveestimation)
      paramold.assign(param);

    double sigmaresp = sqrt(likep->get_scale());
    lambda = likep->get_scale()/tau2;

    // TEST
    //  ofstream out("c:\\bayesx\\testh\\results\\responseRE.res");
    //  likep->response.prettyPrint(out);
    // TEST

    designp->compute_partres(partres,beta);

    if ( (likep->wtype==wweightschange_weightsneqone)  ||
         (likep->wtype==wweightschange_weightsone) ||
         (designp->changingdesign)
       )
       designp->compute_XtransposedWX();

    designp->compute_XtransposedWres(partres, lambda);

    if ((likep->wtype==wweightschange_weightsneqone) ||
        (likep->wtype==wweightschange_weightsone) ||
        (designp->changingdesign) ||
        (!lambdaconst)
        )
      {
      designp->compute_precision(lambda);
      }

    double * work = paramhelp.getV();
    unsigned i;
    unsigned nrpar = param.rows();
    for(i=0;i<nrpar;i++,work++)
      *work = sigmaresp*rand_normal();

    designp->precision.solveU(paramhelp);

    designp->precision.solve(*(designp->XWres_p),paramhelp,param);

    perform_centering();

    designp->compute_f(param,paramlin,beta,fsample.beta);

    if (derivative)
      designp->compute_f_derivative(param,paramlin,derivativesample.beta,
                                    derivativesample.beta);


    betadiff.minus(beta,betaold);

    bool ok;
    if (optionsp->saveestimation)
      {
      ok = designp->update_linpred_save(betadiff);

      if (!ok)
        {
        outsidelinpredlimits++;
        betadiff.minus(betaold,beta);
        designp->update_linpred(betadiff);
        beta.assign(betaold);
        param.assign(paramold);
        }
      else
        acceptance++;
      }
    else
      {
      designp->update_linpred(betadiff);
      ok = true;
      acceptance++;
      }


    // TEST
    // ofstream out("c:\\bayesx\\testh\\results\\beta_re.res");
    // beta.prettyPrint(out);
    // TEST

    if (designp->position_lin!=-1)
      {
      fsample.update();
      }

    paramsample.beta.assign(param);

    paramsample.update();

    if (derivative)
      derivativesample.update();

    FC::update();
    }

  }


void FC_nonp::update_isotonic(void)
  {

  // TEST

   // ofstream out0("c:\\bayesx\\testh\\results\\beta_f.res");
   // beta.prettyPrint(out0);
   // out0.close();

   // ofstream out00("c:\\bayesx\\testh\\results\\intvar_f.res");
   // designp->intvar.prettyPrint(out00);
   // out00.close();

  // TEST

  unsigned i,j;

  bool lambdaconst = false;

  double sigma2resp = likep->get_scale();
  lambda = likep->get_scale()/tau2;

  betaold.assign(beta);
  if (optionsp->saveestimation)
    paramold.assign(param);


  designp->compute_partres(partres,beta);

  if ( (likep->wtype==wweightschange_weightsneqone)  ||
       (likep->wtype==wweightschange_weightsone) ||
       (designp->changingdesign)
     )
    designp->compute_XtransposedWX();

  designp->compute_XtransposedWres(partres, lambda);

  if ((likep->wtype==wweightschange_weightsneqone) ||
      (likep->wtype==wweightschange_weightsone) ||
      (designp->changingdesign) ||
      (!lambdaconst)
      )
    {
    designp->compute_precision(lambda);
    }

  int count = 0;
  int maxit = 20;
  double mu;
  double s;
  double * paramp;
  double * parampi;
  double * XWresp;


  while(count < maxit)
    {

    XWresp = (*(designp->XWres_p)).getV();
    parampi = param.getV();
    for (i=0;i<param.rows();i++,XWresp++,parampi++)
      {

      mu = 0;
      paramp = param.getV();
      for (j=0;j<i;j++,paramp++)    // links
        {
        // mu+= param(j,0)*designp->precision(i,j);
        mu+= (*paramp) * designp->precision(i,j);
        }



      paramp = param.getV()+i+1;
      for (j=i+1;j<param.rows();j++,paramp++)  // rechts
        {
        // mu+= param(j,0)*designp->precision(i,j);
        mu+= (*paramp)*designp->precision(i,j);
        }

      // mu = ((*(designp->XWres_p))(i,0) -mu)/designp->precision(i,i);
      mu = (*XWresp -mu)/designp->precision(i,i);

      s = sqrt(sigma2resp/designp->precision(i,i));

      if(i==0)
        {
        if(stype==increasing)
          *parampi = trunc_normal2(-20,param(1,0),mu,s);
        else
          *parampi = trunc_normal2(param(1,0),20,mu,s);
        }
      else if(i==param.rows()-1)
        {
        if(stype==increasing)
          *parampi = trunc_normal2(param(param.rows()-2,0),20,mu,s);
        else
          *parampi = trunc_normal2(-20,param(param.rows()-2,0),mu,s);
        }
      else
        {
        if(stype==increasing)
          {
          *parampi = trunc_normal2(param(i-1,0),param(i+1,0),mu,s);
          }
        else
          *parampi = trunc_normal2(param(i+1,0),param(i-1,0),mu,s);
        }

      }

    count++;
    }

  /*
  TEST
  ofstream out("c:\\bayesx\\test\\results\\paramhelp.res");
  paramhelp.prettyPrint(out);

  ofstream out2("c:\\bayesx\\test\\results\\param.res");
  param.prettyPrint(out2);
  TEST
  */

  perform_centering();

  designp->compute_f(param,paramlin,beta,fsample.beta);

  betadiff.minus(beta,betaold);

  bool ok;
  if (optionsp->saveestimation)
    {
    ok = designp->update_linpred_save(betadiff);

    if (!ok)
      {
      outsidelinpredlimits++;
      betadiff.minus(betaold,beta);
      designp->update_linpred(betadiff);
      beta.assign(betaold);
      param.assign(paramold);
      }
    else
      acceptance++;
    }
  else
    {
    designp->update_linpred(betadiff);
    ok = true;
    acceptance++;
    }

  if (designp->position_lin!=-1)
    {
    fsample.update();
    }

  paramsample.beta.assign(param);
  paramsample.update();

  if (derivative)
    {
    designp->compute_f_derivative(param,paramlin,derivativesample.beta,
                                    derivativesample.beta);

    derivativesample.update();
    }

  FC::update();

  }


bool FC_nonp::posteriormode_transform(void)
  {

  double h = likep->compute_iwls(true,false);

  betaold.assign(beta);

  designp->compute_partres(partres,beta,true);

  if (designp->QtRinv.rows() <= 1)
    {
    designp->compute_orthogonaldecomp();
    acuteparam = datamatrix(param.rows(),1,0);
    }

  designp->compute_XtransposedWres(partres, lambda);

  designp->u.mult(designp->QtRinv,*(designp->XWres_p));

  unsigned j;
  double * up = designp->u.getV();
  double * sp = designp->s.getV();
  double * acuteparamp = acuteparam.getV();
  for (j=0;j<param.rows();j++,up++,sp++,acuteparamp++)
    {
    h = 1/(1+lambda* (*sp));
    *acuteparamp = h * (*up);
    }

  param.mult(designp->RtinvQ,acuteparam);

  if(designp->center)
    centerparam();

  if (designp->position_lin!=-1)
    {
    get_linparam();
    }


  designp->compute_f(param,paramlin,beta,fsample.beta);

  betadiff.minus(beta,betaold);

  designp->update_linpred(betadiff);

  if (designp->position_lin!=-1)
    {
    fsample.posteriormode_betamean();
    }


  ST::string n = designp->datanames[0];

  designp->compute_meaneffect(masterp->level1_likep[equationnr],meaneffect,beta,
                              meaneffect_sample.beta,computemeaneffect,
                              meaneffectconstant);

  return FC::posteriormode();

  }



bool FC_nonp::posteriormode(void)
  {

  if (orthogonal)
    return posteriormode_transform();
  else
    {

    double h = likep->compute_iwls(true,false);

    betaold.assign(beta);

    designp->compute_partres(partres,beta,true);

    designp->compute_XtransposedWX();
    designp->compute_XtransposedWres(partres, lambda);

    designp->compute_precision(lambda);

    designp->precision.solve(*(designp->XWres_p),param);

    if(designp->center)
      {
      if ((designp->centermethod==integralsimple) ||
          (designp->centermethod==meanf) ||
          (designp->centermethod==meanfd)
         )
        centerparam_weight();
      else
        centerparam();
      }

    if (designp->position_lin!=-1)
      {
      get_linparam();
      }

    designp->compute_f(param,paramlin,beta,fsample.beta);

    betadiff.minus(beta,betaold);

    bool ok;

    if (optionsp->saveestimation)
      {
      ok = designp->update_linpred_save(betadiff);
      }
    else
      {
      ok = true;
      designp->update_linpred(betadiff);
      }


    if (!ok)
      {
      betadiff.minus(betaold,beta);

      designp->update_linpred(betadiff);

      beta.assign(betaold);
      }


    if (designp->position_lin!=-1)
      {
      fsample.posteriormode_betamean();
      }


    designp->compute_meaneffect(masterp->level1_likep[equationnr],meaneffect,
                                beta,meaneffect_sample.beta,computemeaneffect,
                                meaneffectconstant);

    return FC::posteriormode();
    }

  }

void FC_nonp::outoptions(void)
  {
  optionsp->out("  " + title + "\n",true);
  optionsp->out("\n");
  designp->outoptions(optionsp);
  }


void FC_nonp::outgraphs(ofstream & out_stata, ofstream & out_R,const ST::string & path)
  {

  ST::string pathps = path.substr(0,path.length()-4) + "_statagraph";

  double u = optionsp->level1;
  double o = optionsp->level2;
  double u1 = optionsp->lower1;
  double u2 = optionsp->upper2;
  double o1 = optionsp->lower2;
  double o2 = optionsp->upper1;
  ST::string u_str = ST::doubletostring(u,0);
  ST::string o_str = ST::doubletostring(o,0);
  ST::string u1_str = ST::doubletostring(u1,5);
  ST::string u2_str = ST::doubletostring(u2,5);
  ST::string o1_str = ST::doubletostring(o1,5);
  ST::string o2_str = ST::doubletostring(o2,5);

  ST::string pu1_str = u1_str.replaceallsigns('.','p');
  ST::string pu2_str = u2_str.replaceallsigns('.','p');
  ST::string po1_str = o1_str.replaceallsigns('.','p');
  ST::string po2_str = o2_str.replaceallsigns('.','p');
  ST::string pu_str = u_str.replaceallsigns('.','p');
  ST::string po_str = o_str.replaceallsigns('.','p');

  ST::string xvar = designp->datanames[designp->datanames.size()-1];

  if (derivative)
    {

    ST::string pathderivative = path.substr(0,path.length()-4) +
                                 "_derivative.res";

    ST::string pathderivativeps = path.substr(0,path.length()-4) +
                                  "_derivative_statagraph";



    out_stata << "clear" << endl
              << "infile intnr " << xvar

              << " pmean pstd pqu" << pu1_str
              << " pqu" << po1_str << " pmed pqu" << po2_str << " pqu" << pu2_str
              << " pcat" << pu_str << " pcat" << po_str
              << " using "
              << pathderivative << endl
              << "drop in 1" << endl;

    out_stata << "graph twoway rarea pqu" << pu1_str << " pqu" << pu2_str
              << " " << xvar << ", bcolor(gs13) || rarea pqu" << po1_str
              << " pqu" << po2_str << " " << xvar << " , bcolor(gs10) || /*"
              << endl << " */ scatter pmean "
              << xvar << ", c(l) m(i) clpattern(l) clcolor(gs0) /* "
              << endl << " */ ytitle(\"Derivative of effect of "
              << xvar << "\") xtitle(\"" << xvar
              << "\") xlab(,grid) ylab(,grid) legend(off)"
              << endl << "graph export " << pathderivativeps << ".eps, replace"
                      << endl << endl;

   out_stata << "sleep 1000" << endl << endl;


    }

  out_stata << "clear" << endl
            << "infile intnr " << xvar
            << " pmean pstd pqu" << pu1_str
            << " pqu" << po1_str << " pmed pqu" << po2_str << " pqu" << pu2_str
            << " pcat" << pu_str << " pcat" << po_str
            << " sim_pqu" << pu1_str
            << " sim_pqu" << po1_str << " sim_pqu" << po2_str << " sim_pqu" << pu2_str
            << " sim_pcat" << pu_str << " sim_pcat" << po_str;



  if (designp->position_lin!=-1)
    {
    out_stata << " pmean_d pstd_d pqu"
              << pu1_str << "_d"
              << " pqu" << po1_str << "_d"
              << " pmed_d pqu" << po2_str << "_d"
              << " pqu" << pu2_str << "_d"
              << " pcat" << pu_str << "_d"
              << " pcat" << po_str << "_d";
    }


  if (computemeaneffect==true)
    {
    out_stata << " pmean_mu pstd_mu pqu"
              << pu1_str << "_mu"
              << " pqu" << po1_str << "_mu"
              << " pmed_d pqu" << po2_str << "_mu"
              << " pqu" << pu2_str << "_mu";
    }


  out_stata << " using "
            << path << endl
            << "drop in 1" << endl;
            if (designp->type == Mrf)
              {
              out_stata << "kdensity pmean" << endl
              << "graph export " << pathps << ".eps, replace"
                      << endl << endl;
              }
            else
              {
              out_stata << "graph twoway rarea pqu" << pu1_str << " pqu" << pu2_str
              << " " << xvar << ", bcolor(gs13) || rarea pqu" << po1_str
              << " pqu" << po2_str << " " << xvar << " , bcolor(gs10) || /*"
              << endl << " */ scatter pmean "
              << xvar << ", c(l) m(i) clpattern(l) clcolor(gs0) /* "
              << endl << " */ ytitle(\"Effect of "
              << xvar << "\") xtitle(\"" << xvar
              << "\") xlab(,grid) ylab(,grid) legend(off)"
              << endl << "graph export " << pathps << ".eps, replace"
                      << endl << endl;
              }

   out_stata << "sleep 1000" << endl << endl;


  if (designp->position_lin!=-1)
    {
    if (designp->type == Mrf)
      {
      out_stata << "kdensity pmean_d" << endl
                << "graph export " << pathps << "_linexluded.eps, replace"
                << endl << endl;

      }
    else
      {

      out_stata << "graph twoway rarea pqu" << pu1_str << "_d"
                << " pqu" << pu2_str << "_d"
                << " " << xvar << ", bcolor(gs13) || rarea pqu" << po1_str << "_d"
                << " pqu" << po2_str << "_d"
                << " " << xvar << " , bcolor(gs10) || /*"
                << endl << " */ scatter pmean_d "
                << xvar << ", c(l) m(i) clpattern(l) clcolor(gs0) /* "
                << endl << " */ ytitle(\"Effect of "
                << xvar << "\") xtitle(\"" << xvar
                << "\") xlab(,grid) ylab(,grid) legend(off)"
                << endl << "graph export " << pathps << "_linexluded.eps, replace"
                      << endl << endl;
      }

    out_stata << "sleep 1000" << endl << endl;
    }


  if (computemeaneffect==true)
    {

    if (designp->type == Mrf)
      {
      out_stata << "kdensity pmean_mu" << endl
                << "graph export " << pathps << "_mu.eps, replace"
                << endl << endl;
      }
    else
      {
      out_stata << "graph twoway rarea pqu" << pu1_str << "_mu"
                << " pqu" << pu2_str << "_mu"
                << " " << xvar << ", bcolor(gs13) || rarea pqu" << po1_str << "_mu"
                << " pqu" << po2_str << "_mu"
                << " " << xvar << " , bcolor(gs10) || /*"
                << endl << " */ scatter pmean_mu "
                << xvar << ", c(l) m(i) clpattern(l) clcolor(gs0) /* "
                << endl << " */ ytitle(\"Effect of "
                << xvar << "\") xtitle(\"" << xvar
                << "\") xlab(,grid) ylab(,grid) legend(off)"
                << endl << "graph export " << pathps << "_mu.eps, replace"
                      << endl << endl;
      }

    out_stata << "sleep 1000" << endl << endl;
    }

  }



void FC_nonp::outresults_derivative(ofstream & out_stata, ofstream & out_R,
                        const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    derivativesample.outresults(out_stata,out_R,"");

    optionsp->out("    Estimated first derivatives are stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    ofstream outres(pathresults.strtochar());

    optionsp->out("\n");

    unsigned i;

    ST::string l1 = ST::doubletostring(optionsp->lower1,4);
    ST::string l2 = ST::doubletostring(optionsp->lower2,4);
    ST::string u1 = ST::doubletostring(optionsp->upper1,4);
    ST::string u2 = ST::doubletostring(optionsp->upper2,4);
    l1 = l1.replaceallsigns('.','p');
    l2 = l2.replaceallsigns('.','p');
    u1 = u1.replaceallsigns('.','p');
    u2 = u2.replaceallsigns('.','p');

    outres << "intnr" << "   ";

    unsigned start=0;
    if (designp->intvar.rows()==designp->data.rows())
      start = 1;

    for (i=start;i<designp->datanames.size();i++)
      outres << designp->datanames[i] << "   ";
    outres << "pmean   ";

    outres << "pqu"  << l1  << "   ";
    outres << "pqu"  << l2  << "   ";
    outres << "pmed   ";
    outres << "pqu"  << u1  << "   ";
    outres << "pqu"  << u2  << "   ";
    outres << "pcat" << optionsp->level1 << "   ";
    outres << "pcat" << optionsp->level2 << "   ";


    outres << endl;


    double * workmean;
    double * workbetaqu_l1_lower_p;
    double * workbetaqu_l2_lower_p;
    double * workbetaqu_l1_upper_p;
    double * workbetaqu_l2_upper_p;
    double * workbetaqu50;

    workmean = derivativesample.betamean.getV();
    workbetaqu_l1_lower_p = derivativesample.betaqu_l1_lower.getV();
    workbetaqu_l2_lower_p = derivativesample.betaqu_l2_lower.getV();
    workbetaqu_l1_upper_p = derivativesample.betaqu_l1_upper.getV();
    workbetaqu_l2_upper_p = derivativesample.betaqu_l2_upper.getV();
    workbetaqu50 = derivativesample.betaqu50.getV();


    unsigned nrpar = beta.rows();
    for(i=0;i<nrpar;i++,
        workmean++,workbetaqu_l1_lower_p++,
                              workbetaqu_l2_lower_p++,workbetaqu50++,
                              workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++
                              )
      {
      outres << (i+1) << "   ";
      outres << designp->effectvalues[i] << "   ";
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


void FC_nonp::outbasis_R(const ST::string & pathbasis)
  {
  ofstream out(pathbasis.strtochar());
  designp->outbasis_R(out);
  }


double FC_nonp::kernel_density(const double & x, const double & h)
  {
  double sum=0;
  double n = designp->data.rows();
  unsigned j;
  double n_j;
  double u;
  double x_j;
  double Kd;

  for (j=0;j<designp->posbeg.size();j++)
    {
    x_j = designp->data(designp->posbeg[j],0);
    u = (x-x_j)/h;
    if ( (u <= 1) && (u >= -1) )
      {
      n_j = designp->posend[j]-designp->posbeg[j]+1;
      Kd = 0.75*(1.0-u*u);
      sum += n_j * Kd;

     }
    }

  return 1/(n*h)*sum;
  }


double FC_nonp::compute_importancemeasure(bool absolute)
  {
  double n = designp->data.rows();

  double si=sqrt(designp->data.var(0));
  double R = (designp->data.quantile(75,0) - designp->data.quantile(25,0))/1.34;

  double h;
  if (si < R)
    h = 0.9*si*pow(n,-0.2);
  else
    h = 0.9*R*pow(n,-0.2);

  unsigned j;
  double x_jm1 = designp->data(designp->posbeg[0],0);
  double x_j;
  double diff;
  double kd_xjm1;
  double kd_xj;
  double sum = 0;

  for (j=1;j<designp->posbeg.size();j++)
    {
    x_j = designp->data(designp->posbeg[j],0);
    diff = x_j - x_jm1;
    kd_xjm1 = kernel_density(x_jm1,h);
    kd_xj = kernel_density(x_j,h);
    if (absolute)
      sum += diff* (fabs(betamean(j-1,0))* kd_xjm1 + fabs(betamean(j,0))*kd_xj)/2;
    else
      sum += diff* (pow(betamean(j-1,0),2)* kd_xjm1 + pow(betamean(j,0),2)*kd_xj)/2;

    x_jm1 = x_j;
    }
  return sum;
  }


double FC_nonp::compute_importancemeasure_discrete(bool absolute)
  {
  double n = designp->data.rows();

  unsigned j;
  double sum = 0;
  double pr;

  for (j=0;j<designp->posbeg.size();j++)
    {
    pr = (designp->posend[j]-designp->posbeg[j]+1)/n;
    if (absolute)
      sum +=  fabs(betamean(j,0))*pr;
    else
      sum +=  pow(betamean(j,0),2)*pr;

    }
  return sum;
  }


void FC_nonp::outresults(ofstream & out_stata, ofstream & out_R,
                        const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    outgraphs(out_stata,out_R,pathresults);

    FC::outresults(out_stata,out_R,"");
    if (designp->position_lin != -1)
      fsample.outresults(out_stata,out_R,"");

    if (computemeaneffect==true)
       meaneffect_sample.outresults(out_stata,out_R,"");

    ST::string pathresultsbeta = pathresults.substr(0,pathresults.length()-4) +
                                 "_param.res";


    ST::string pathbasis = pathresults.substr(0,pathresults.length()-4) +
                                 "_basisR.res";

    outbasis_R(pathbasis);

    outresults_acceptance();

    paramsample.outresults(out_stata,out_R,pathresultsbeta);

    ST::string paths = pathresults.substr(0,pathresults.length()-4) +
                                 "_sample.raw";

    out_R << "family=" << likep->familyshort.strtochar() << ",";
    out_R << "hlevel=" << likep->hlevel << ",";
    out_R << "equationtype=" << likep->equationtype.strtochar() << ",";
    if (designp->intvar.rows()==designp->data.rows())
      out_R << "term=sx("  << designp->datanames[1].strtochar() << ",by=" << designp->datanames[0].strtochar() << "),";
    else
      out_R << "term=sx("  << designp->datanames[0].strtochar()  << "),";
    out_R << "filetype=nonlinear,";
    out_R << "pathsamples=" << paths.strtochar() << ",";
    out_R << "pathbasis=" << pathbasis.strtochar() << ",";

    optionsp->out("    Estimated parameters are stored in file\n");
    optionsp->out("    " +  pathresultsbeta + "\n");
    optionsp->out("\n");

    optionsp->out("    Function estimates are stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

    // ---------------------------- compute importancemeasures ------------------------

    if (imeasures)
      {
      double im_absolute;
      double im_var;
      if (designp->discrete)
        {
        im_var = compute_importancemeasure_discrete(false);
        im_absolute = compute_importancemeasure_discrete(true);
        }
      else
        {
        im_var = compute_importancemeasure(false);
        im_absolute = compute_importancemeasure(true);
        }

      if (designp->intvar.rows()==designp->data.rows())
        {
        double kintvar = designp->compute_kernel_intvar(true);
        im_absolute *= kintvar;
        kintvar = designp->compute_kernel_intvar(false);
        im_var *= kintvar;
        }

      optionsp->out("    Importance measures\n");
      optionsp->out("      based on absolute function values: " + ST::doubletostring(im_absolute,6) + "\n");
      optionsp->out("      based on squared function values:  " + ST::doubletostring(im_var,6) + "\n");
      optionsp->out("\n");
      }

    // ---------------------------- end: compute importancemeasures -----------------------------------------

    if (designp->intvar.rows()==designp->data.rows())
      {
      optionsp->out("    Mean effects evaluated at " +
                  designp->datanames[0] + "=" + ST::doubletostring(designp->meaneffectintvar,6)
                   + " and " + designp->datanames[1] + "=" +
                   designp->effectvalues[designp->meaneffectnr]
                  + " with value " +
                  ST::doubletostring(designp->meaneffectintvar * betamean(designp->meaneffectnr,0),6) + "\n");
      }
    else
      {
      optionsp->out("    Mean effects evaluated at " +
                  designp->datanames[designp->datanames.size()-1] + "=" +
                  designp->effectvalues[designp->meaneffectnr]
                  + " with value " +
                  ST::doubletostring(betamean(designp->meaneffectnr,0)) + "\n");
      }

    double s_level1=0;
    double s_level2=0;
    if (optionsp->samplesize > 0)
      {
      optionsp->out("\n");
      s_level1 = simconfBand(true);
      s_level2 = simconfBand(false);

      optionsp->out("    Scaling factor to blow up pointwise " +
                   ST::inttostring(optionsp->level1) + " percent credible intervals\n");
      optionsp->out("    to obtain simultaneous credible intervals: " +
           ST::doubletostring(s_level2,6) + "\n");

      optionsp->out("\n");

      optionsp->out("    Scaling factor to blow up pointwise " +
                   ST::inttostring(optionsp->level2) + " percent credible intervals\n");
      optionsp->out("    to obtain simultaneous credible intervals: " +
           ST::doubletostring(s_level1,6) + "\n");

      optionsp->out("\n");
      }


    if ((optionsp->samplesize > 1) && (derivative==true))
      {
      ST::string pathresultsderivative =
                         pathresults.substr(0,pathresults.length()-4) +
                                 "_derivative.res";

      outresults_derivative(out_stata,out_R,pathresultsderivative);
      }


    ofstream outres(pathresults.strtochar());

    optionsp->out("\n");

    unsigned i;

    ST::string l1 = ST::doubletostring(optionsp->lower1,4);
    ST::string l2 = ST::doubletostring(optionsp->lower2,4);
    ST::string u1 = ST::doubletostring(optionsp->upper1,4);
    ST::string u2 = ST::doubletostring(optionsp->upper2,4);
    l1 = l1.replaceallsigns('.','p');
    l2 = l2.replaceallsigns('.','p');
    u1 = u1.replaceallsigns('.','p');
    u2 = u2.replaceallsigns('.','p');

    outres << "intnr" << "   ";

    unsigned start=0;
    if (designp->intvar.rows()==designp->data.rows())
      start = 1;

    for (i=start;i<designp->datanames.size();i++)
      outres << designp->datanames[i] << "   ";
    outres << "pmean   ";

    outres << "pstd   ";

      outres << "pqu"  << l1  << "   ";
      outres << "pqu"  << l2  << "   ";
      outres << "pmed   ";
      outres << "pqu"  << u1  << "   ";
      outres << "pqu"  << u2  << "   ";
      outres << "pcat" << optionsp->level1 << "   ";
      outres << "pcat" << optionsp->level2 << "   ";

      outres << "pqu"  << l1  << "_sim   ";
      outres << "pqu"  << l2  << "_sim   ";
      outres << "pqu"  << u1  << "_sim   ";
      outres << "pqu"  << u2  << "_sim   ";
      outres << "pcat" << optionsp->level1 << "_sim   ";
      outres << "pcat" << optionsp->level2 << "_sim   ";



    if (designp->position_lin!=-1)
      {

      outres << "pmean_d   ";
      outres << "pstd_d   ";

      if (optionsp->samplesize > 1)
        {
        outres << "pqu"  << l1  << "_d   ";
        outres << "pqu"  << l2  << "_d   ";
        outres << "pmed_d   ";
        outres << "pqu"  << u1  << "_d   ";
        outres << "pqu"  << u2  << "_d   ";
        outres << "pcat" << optionsp->level1 << "_d   ";
        outres << "pcat" << optionsp->level2 << "_d   ";
        }
      }


    if (computemeaneffect==true)
      {

      outres << "pmean_mu   ";
      outres << "pstd_mu   ";

      if (optionsp->samplesize > 1)
        {
        outres << "pqu"  << l1  << "_mu   ";
        outres << "pqu"  << l2  << "_mu   ";
        outres << "pmed_mu   ";
        outres << "pqu"  << u1  << "_mu   ";
        outres << "pqu"  << u2  << "_mu   ";
        }

      }


    outres << endl;


    double * workmean=NULL;
    double * workstd=NULL;
    double * workbetaqu_l1_lower_p=NULL;
    double * workbetaqu_l2_lower_p=NULL;
    double * workbetaqu_l1_upper_p=NULL;
    double * workbetaqu_l2_upper_p=NULL;
    double * workbetaqu50=NULL;

    double * dworkmean=NULL;
    double * dworkstd=NULL;
    double * dworkbetaqu_l1_lower_p=NULL;
    double * dworkbetaqu_l2_lower_p=NULL;
    double * dworkbetaqu_l1_upper_p=NULL;
    double * dworkbetaqu_l2_upper_p=NULL;
    double * dworkbetaqu50=NULL;


    double * mu_workmean=NULL;
    double * mu_workstd=NULL;
    double * mu_workbetaqu_l1_lower_p=NULL;
    double * mu_workbetaqu_l2_lower_p=NULL;
    double * mu_workbetaqu_l1_upper_p=NULL;
    double * mu_workbetaqu_l2_upper_p=NULL;
    double * mu_workbetaqu50=NULL;


    if (designp->position_lin!=-1)
      {
      workmean = fsample.betamean.getV();
      workstd = fsample.betavar.getV();
      workbetaqu_l1_lower_p = fsample.betaqu_l1_lower.getV();
      workbetaqu_l2_lower_p = fsample.betaqu_l2_lower.getV();
      workbetaqu_l1_upper_p = fsample.betaqu_l1_upper.getV();
      workbetaqu_l2_upper_p = fsample.betaqu_l2_upper.getV();
      workbetaqu50 = fsample.betaqu50.getV();

      dworkmean = betamean.getV();
      dworkstd = betavar.getV();
      dworkbetaqu_l1_lower_p = betaqu_l1_lower.getV();
      dworkbetaqu_l2_lower_p = betaqu_l2_lower.getV();
      dworkbetaqu_l1_upper_p = betaqu_l1_upper.getV();
      dworkbetaqu_l2_upper_p = betaqu_l2_upper.getV();
      dworkbetaqu50 = betaqu50.getV();

      }
    else
      {
      workmean = betamean.getV();
      workstd = betavar.getV();
      workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
      workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
      workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
      workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
      workbetaqu50 = betaqu50.getV();
      }


    if (computemeaneffect==true)
      {
      mu_workmean = meaneffect_sample.betamean.getV();
      mu_workstd = meaneffect_sample.betavar.getV();
      mu_workbetaqu_l1_lower_p = meaneffect_sample.betaqu_l1_lower.getV();
      mu_workbetaqu_l2_lower_p = meaneffect_sample.betaqu_l2_lower.getV();
      mu_workbetaqu_l1_upper_p = meaneffect_sample.betaqu_l1_upper.getV();
      mu_workbetaqu_l2_upper_p = meaneffect_sample.betaqu_l2_upper.getV();
      mu_workbetaqu50 = meaneffect_sample.betaqu50.getV();
      }

    double l1_sim,l2_sim,u1_sim,u2_sim;

    unsigned nrpar = beta.rows();
    for(i=0;i<nrpar;i++,workmean++,workstd++,workbetaqu_l1_lower_p++,
                              workbetaqu_l2_lower_p++,workbetaqu50++,
                              workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
      {
      outres << (i+1) << "   ";
      outres << designp->effectvalues[i] << "   ";
      outres << *workmean << "   ";

      if (optionsp->samplesize > 1)
        {
        if (*workstd < 0.0000000000001)
          outres << 0 << "   ";
        else
          outres << sqrt(*workstd) << "   ";
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

        l1_sim = *workmean - s_level1*(*workmean- *workbetaqu_l1_lower_p);
        l2_sim = *workmean - s_level2*(*workmean- *workbetaqu_l2_lower_p);
        u1_sim = *workmean + s_level1*(*workbetaqu_l1_upper_p - *workmean);
        u2_sim = *workmean + s_level2*(*workbetaqu_l2_upper_p - *workmean);

        outres << l1_sim << "   ";
        outres << l2_sim << "   ";
        outres << u2_sim << "   ";
        outres << u1_sim << "   ";

        if (l1_sim > 0)
          outres << 1 << "   ";
        else if (u1_sim < 0)
          outres << -1 << "   ";
        else
          outres << 0 << "   ";

        if (l2_sim > 0)
          outres << 1 << "   ";
        else if (u2_sim < 0)
          outres << -1 << "   ";
        else
          outres << 0 << "   ";

        }
      else
        outres << "0  0  0  0  0  0  0  0  0  0  0  0  0  0  ";


      if (designp->position_lin!=-1)
        {

        outres << *dworkmean << "   ";

        if (optionsp->samplesize > 1)
          {

          if ((*dworkstd) < 0.0000000000001)
            outres << 0 << "   ";
          else
            outres << sqrt(*dworkstd) << "   ";

          outres << *dworkbetaqu_l1_lower_p << "   ";
          outres << *dworkbetaqu_l2_lower_p << "   ";
          outres << *dworkbetaqu50 << "   ";
          outres << *dworkbetaqu_l2_upper_p << "   ";
          outres << *dworkbetaqu_l1_upper_p << "   ";

          if (*dworkbetaqu_l1_lower_p > 0)
            outres << 1 << "   ";
          else if (*dworkbetaqu_l1_upper_p < 0)
            outres << -1 << "   ";
          else
            outres << 0 << "   ";

          if (*dworkbetaqu_l2_lower_p > 0)
            outres << 1 << "   ";
          else if (*dworkbetaqu_l2_upper_p < 0)
            outres << -1 << "   ";
          else
            outres << 0 << "   ";

          }

        if (i <nrpar-1)
          {
          dworkmean++;
          dworkstd++;
          dworkbetaqu_l1_lower_p++;
          dworkbetaqu_l2_lower_p++;
          dworkbetaqu50++;
          dworkbetaqu_l1_upper_p++;
          dworkbetaqu_l2_upper_p++;
          }

        }


      if (computemeaneffect==true)
        {

        outres << *mu_workmean << "   ";

        if (optionsp->samplesize > 1)
          {
          if (*mu_workstd < 0.0000000000001)
            outres << 0 << "   ";
          else
            outres << sqrt(*mu_workstd) << "   ";

          outres << *mu_workbetaqu_l1_lower_p << "   ";
          outres << *mu_workbetaqu_l2_lower_p << "   ";
          outres << *mu_workbetaqu50 << "   ";
          outres << *mu_workbetaqu_l2_upper_p << "   ";
          outres << *mu_workbetaqu_l1_upper_p << "   ";
          }

        if (i <nrpar-1)
          {
          mu_workmean++;
          mu_workstd++;
          mu_workbetaqu_l1_lower_p++;
          mu_workbetaqu_l2_lower_p++;
          mu_workbetaqu50++;
          mu_workbetaqu_l1_upper_p++;
          mu_workbetaqu_l2_upper_p++;
          }

        }

      outres << endl;
      }

    }

  }


void FC_nonp::initialize_center(void)
  {
  int nrrest = designp->basisNull.rows();
  int nrpar = param.rows();
  Vcenter = datamatrix(nrpar,nrrest);
  Vcentert = datamatrix(nrrest,nrpar);
  Wcenter = datamatrix(nrrest,nrrest);
  Ucenter = datamatrix(nrrest,nrpar);
  ccenter = datamatrix(nrrest,1);
  Utc = datamatrix(nrpar,1);
  }


void FC_nonp::centerparam_sample_var(void)
  {
  unsigned j;
  double * paramp = param.getV();
  double sumparam = 0;
  double sumvar =  0;
  for (j=0;j<param.rows();j++,paramp++)
    {
    sumparam += (*paramp);
    sumvar += 1/designp->precision.getDiag(j);
    }

  double c = sumparam/sumvar;

  paramp = param.getV();
  for (j=0;j<param.rows();j++,paramp++)
    {
    *paramp -= c/designp->precision.getDiag(j);
    }

  }

void FC_nonp::centerparam_sample(void)
  {

  int nrrest = designp->basisNull.rows();
  int nrpar = param.rows();

  if ((int(Vcenter.rows()) != nrpar) ||
      (Vcenter.cols() != designp->basisNull.rows()))
    initialize_center();

  int i,j;
  double * helpcenterp = helpcenter.getV();

  double * Vcentertp = Vcentert.getV();
  double * Vcenterp;
  for (i=0;i<nrrest;i++)
    {
    Vcenterp = Vcenter.getV()+i;

    designp->precision.solve(designp->basisNullt[i],helpcenter);

    helpcenterp = helpcenter.getV();

    for (j=0;j<nrpar;j++,helpcenterp++,Vcentertp++,Vcenterp+=nrrest)
      {
      *Vcenterp = *helpcenterp;
      *Vcentertp = *helpcenterp;
      }
    }

  // TEST
   /*
   ofstream out0("c:\\bayesx\\testh\\results\\basisnull.res");
   (designp->basisNullt[0]).prettyPrint(out0);

   ofstream out("c:\\bayesx\\testh\\results\\Vcenter.res");
   Vcenter.prettyPrint(out);

   ofstream out2("c:\\bayesx\\testh\\results\\praecision.res");
   designp->precision.print4(out2);
   */
  // TEST

  Wcenter.mult(designp->basisNull,Vcenter);
  Ucenter = Wcenter.inverse()*Vcentert;
  ccenter.mult(designp->basisNull,param);
  Utc = Ucenter.transposed()*ccenter;

  //  TEST

  // ofstream out4("c:\\bayesx\\testh\\results\\param.res");
  // param.prettyPrint(out4);

  //  TEST

  param.minus(param,Utc);

  //  TEST
  // ofstream out5("c:\\bayesx\\testh\\results\\paramneu.res");
  // param.prettyPrint(out5);

  // ofstream out6("c:\\bayesx\\testh\\results\\Utc.res");
  // Utc.prettyPrint(out6);
  //  TEST


  //  TEST
  //  ofstream out2("c:\\bayesx\\test\\results\\param.res");
  //  param.prettyPrint(out2);
  //  TEST
  }



void FC_nonp::get_effect(datamatrix & effect)
  {

  designp->compute_effect(effect,beta,Varcoefftotal);

  }


void FC_nonp::centerparam(void)
  {

  unsigned i;
  double sum=0;
  double * workparam = param.getV();
  unsigned nrparam = param.rows();

  for (i=0;i<nrparam;i++,workparam++)
    {
    sum+= *workparam;
    }

  workparam = param.getV();

  sum /= double(nrparam);

  for (i=0;i<nrparam;i++,workparam++)
    *workparam-= sum;

  }


void FC_nonp::centerparam_weight(void)
  {

  unsigned i;
  double sum=0;
  double * workparam = param.getV();
  unsigned nrparam = param.rows();
  double * wp = designp->basisNull.getV();

  double sumw=0;

  for (i=0;i<nrparam;i++,workparam++,wp++)
    {
    sum+= (*workparam)*(*wp);
    sumw+= (*wp);
    }

  sum /= sumw;

  workparam = param.getV();

  for (i=0;i<nrparam;i++,workparam++)
    *workparam-= sum;

  }



void FC_nonp::centerparam_sum2(double & s2)
  {

  unsigned i;
  double sum=0;
  double sum2=0;
  double * workparam = param.getV();
  unsigned nrparam = param.rows();

  for (i=0;i<nrparam;i++,workparam++)
    {
    sum+= *workparam;
    }

  workparam = param.getV();

  sum /= double(nrparam);


  workparam = param.getV();

  for (i=0;i<nrparam;i++,workparam++)
    {
    sum2+= pow(*workparam-sum,2);
    }


  workparam = param.getV();
  for (i=0;i<nrparam;i++,workparam++)
    *workparam = sqrt(s2/sum2) * (*workparam-sum);

//   ofstream out("c:\\bayesx\\testh\\results\\param.res");
//   param.prettyPrint(out);

  }



void FC_nonp::compute_autocorr_all(const ST::string & path,
                                      unsigned lag, ofstream & outg) const
  {
  paramsample.compute_autocorr_all(path,lag,outg);
  }


void FC_nonp::get_samples(const ST::string & filename,ofstream & outg) const
  {
  paramsample.get_samples(filename,outg);

  if (derivative && samplederivative)
    {
    ST::string filed = filename.substr(0,filename.length()-11) +
                                 "_derivative_sample.raw";

    derivativesample.get_samples(filed,outg);
    }

  if (samplef)
    {
    ST::string filef = filename.substr(0,filename.length()-11) +
                                 "_function_sample.raw";

    FC::get_samples(filef,outg);
    }

  }



void FC_nonp::reset(void)
  {

  }


} // end: namespace MCMC





