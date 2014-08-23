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



#include "FC_hrandom.h"


//------------------------------------------------------------------------------
//-------------- CLASS: FC_hrandom implementation of member functions ----------
//------------------------------------------------------------------------------


namespace MCMC
{


FC_hrandom::FC_hrandom(void)
  {
  }


void FC_hrandom::read_options(vector<ST::string> & op,vector<ST::string> & vn)
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

  17     internal_multexp
  */

/*  if (op[14] == "increasing")
    stype = increasing;
  else if (op[14] == "decreasing")
    stype = decreasing;
  else
*/
    stype = unconstrained;

  rtype = additive;
  if (op[12] == "true")
    rtype = mult;

  if (op[17] == "true")
    rtype = multexp;

  }

FC_hrandom::FC_hrandom(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,DISTR * lp,DISTR * lp_RE,
                 const ST::string & t,const ST::string & fp,
                 const ST::string & fp2, DESIGN * Dp,
                 vector<ST::string> & op, vector<ST::string> & vn)
     : FC_nonp(mp,enr,o,lp,t,fp,Dp,op,vn)
  {
  read_options(op,vn);
  likep_RE = lp_RE;
  FCrcoeff = FC(o,"",beta.rows(),beta.cols(),fp2);
  derivative=false;
  }


FC_hrandom::FC_hrandom(const FC_hrandom & m)
  : FC_nonp(FC_nonp(m))
  {
  rtype = m.rtype;
  likep_RE = m.likep_RE;
  FCrcoeff = m.FCrcoeff;
  response_o = m.response_o;
  linpred_o = m.linpred_o;
  likelihoodc = m.likelihoodc;
  likelihoodn = m.likelihoodn;
  beta_prior = m.beta_prior;
  }


const FC_hrandom & FC_hrandom::operator=(const FC_hrandom & m)
  {

  if (this==&m)
	 return *this;
  FC_nonp::operator=(FC_nonp(m));
  rtype = m.rtype;
  likep_RE = m.likep_RE;
  FCrcoeff = m.FCrcoeff;
  response_o = m.response_o;
  linpred_o = m.linpred_o;
  likelihoodc = m.likelihoodc;
  likelihoodn = m.likelihoodn;
  beta_prior = m.beta_prior;
  return *this;
  }


void FC_hrandom::set_rcoeff(void)
  {
  unsigned i;
  double * betap = beta.getV();
  double * betarcoeffp = FCrcoeff.beta.getV();


  double * linpredREp;
  if (likep_RE->linpred_current==1)
    linpredREp = likep_RE->linearpred1.getV();
  else
    linpredREp = likep_RE->linearpred2.getV();

  for (i=0;i<beta.rows();i++,betap++,betarcoeffp++,linpredREp++)
    *betarcoeffp = *betap - *linpredREp;

  }


void FC_hrandom::update_IWLS(void)
  {
  unsigned i;

  lambda = likep->get_scale()/tau2;

  if (optionsp->nriter == 1)
    {
    betaold.assign(beta);
    }

  double * betap = beta.getV();
  double * betaoldp = betaold.getV();

  double * linpredREp;
  if (likep_RE->linpred_current==1)
    linpredREp = likep_RE->linearpred1.getV();
  else
    linpredREp = likep_RE->linearpred2.getV();

  if (likelihoodc.rows() <=1)
    {
    likelihoodc = datamatrix(beta.rows(),1,0);
    likelihoodn = datamatrix(beta.rows(),1,0);
    }

  double postmode;
  double diff;
  double var;
  double u;
  double xwres;


  likep->compute_iwls(true,likelihoodc,designp->ind);

  designp->compute_partres(partres,beta);

  double * workpartres = partres.getV();
  double * worklikelihoodc = likelihoodc.getV();
  double * workWsum = designp->Wsum.getV();

  for (i=0;i<beta.rows();i++,betap++,linpredREp++,
       workpartres++,worklikelihoodc++,workWsum++)

    {

    *worklikelihoodc  -= 0.5*pow((*betap)-(*linpredREp),2)/tau2;

    xwres =  lambda*(*linpredREp)+ (*workpartres);

    var = 1/(*workWsum+lambda);
    postmode =  var * xwres;
    *betap = postmode + sqrt(var)*rand_normal();
    diff = *betap - postmode;
    *worklikelihoodc += -1.0/(2*var)* pow(diff,2)-0.5*log(var);
    }


  betadiff.minus(beta,betaold);

  bool ok;
  if (optionsp->saveestimation)
    {
    ok = designp->update_linpred_save(betadiff);
    if (!ok)
      {
      int i = optionsp->nriter;
      outsidelinpredlimits++;
      }
    }
  else
    {
    designp->update_linpred(betadiff);
    ok = true;
    }


  if (ok)
    {
    likep->compute_iwls(true,likelihoodn,designp->ind);
    designp->compute_partres(partres,beta);

    workpartres = partres.getV();
    double * worklikelihoodn = likelihoodn.getV();
    worklikelihoodc = likelihoodc.getV();
    workWsum = designp->Wsum.getV();

    betap = beta.getV();
    betaoldp = betaold.getV();

    if (likep_RE->linpred_current==1)
      linpredREp = likep_RE->linearpred1.getV();
    else
      linpredREp = likep_RE->linearpred2.getV();

    double * betadiffp = betadiff.getV();

    for (i=0;i<beta.rows();i++,betap++,linpredREp++,betadiffp++,
         betaoldp++,workpartres++,worklikelihoodn++,workWsum++,worklikelihoodc++)

      {

      *worklikelihoodn -= 0.5*pow((*betap)-(*linpredREp),2)/tau2;

      xwres =  lambda*(*linpredREp)+ (*workpartres);


      var = 1/(*workWsum+lambda);
      diff = *betaoldp - var * xwres;

      *worklikelihoodn += -1.0/(2*var)* pow(diff,2)-0.5*log(var);


      nrtrials++;
      u = log(uniform());
      if (u <= (*worklikelihoodn) - (*worklikelihoodc))
        {
        acceptance++;
        *betaoldp = *betap;
        *betadiffp = 0;
        }
      else
        {
        *betadiffp = *betaoldp - *betap;
        *betap = *betaoldp;
        }

      }

    } // if ok
  else
    betadiff.minus(betaold,beta);

  designp->update_linpred(betadiff);

  FC::update();

  }



void FC_hrandom::update(void)
  {

  if (IWLS)
    {
    update_IWLS();
    }
  else
    FC_nonp::update();

  set_rcoeff();

  FCrcoeff.acceptance++;
  FCrcoeff.update();

  likep_RE->workingresponse.assign(beta);
  likep_RE->response.assign(beta);

  // TEST
  //  ofstream out("c:\\bayesx\\testh\\results\\response_h.res");
  //  beta.prettyPrint(out);
  // END TEST

  }


void FC_hrandom::sample_for_cv(datamatrix & pred)
  {
  if (beta_prior.rows() == 1)
    beta_prior = beta;

  unsigned i;
  double * workbeta = beta_prior.getV();
  double tau = sqrt(tau2);
  for (i=0;i<beta_prior.rows();i++,workbeta++)
    *workbeta =  tau*rand_normal();

  designp->compute_effect(pred,beta_prior,Varcoefftotal);

  }


void FC_hrandom::compute_effect_cv(datamatrix & effect)
  {

  int i;

  if (effect.rows() != designp->data.rows())
    effect = datamatrix(designp->data.rows(),1,0);

  double * workintvar = designp->intvar.getV();
  unsigned * workind = designp->ind.getV();
  double * workeffect = effect.getV();

  int size = designp->ind.rows();

  if (likep_RE->linpred_current==1)
    {
    if (designp->intvar.rows() != designp->data.rows())
      {
      for (i=0;i<size;i++,workind++,workeffect++)
        *workeffect = beta(*workind,0)- likep_RE->linearpred1(*workind,0);
      }
    else
      {
      for (i=0;i<size;i++,workind++,workeffect++,workintvar++)
        *workeffect = *workintvar * (beta(*workind,0) -
        likep_RE->linearpred1(*workind,0));
      }
    }
  else
    {

    if (designp->intvar.rows() != designp->data.rows())
      {
      for (i=0;i<size;i++,workind++,workeffect++)
        *workeffect = beta(*workind,0)- likep_RE->linearpred2(*workind,0);
      }
    else
      {
      for (i=0;i<size;i++,workind++,workeffect++,workintvar++)
        *workeffect = *workintvar * (beta(*workind,0) -
        likep_RE->linearpred2(*workind,0));
      }

    }

  }



void FC_hrandom::update_response_multexp(void)
  {
/*
  unsigned i,j;

  int size = designp->posbeg.size();

  vector<int>::iterator itbeg = designp->posbeg.begin();
  vector<int>::iterator itend = designp->posend.begin();

  double * * linpredp;

  if (likep->linpred_current==1)
    {
    linpredp = designp->linpredp1.getV();
    linpred_o.assign(likep->linearpred1);
    }
  else
    {
    linpredp = designp->linpredp2.getV();
    linpred_o.assign(likep->linearpred2);
    }


  double ** responsepp = designp->responsep.getV();

  double * workintvar = designp->intvar.getV();

  double * betap = beta.getV();

  for (j=0;j<size;j++,++itbeg,++itend,betap++)
    {
    for (i=*itbeg;i<=*itend;i++,linpredp++,responsepp++,workintvar++)
      {
      *(*responsepp) = *(*responsepp) - *(*linpredp) + exp((*betap)  * (*workintvar));
      *(*linpredp) = *betap * (*workintvar);
      }
    }

  */
  }


void FC_hrandom::update_linpred_multexp(void)
  {

  /*
  unsigned i,j;

  int size = designp->posbeg.size();

  vector<int>::iterator itbeg = designp->posbeg.begin();
  vector<int>::iterator itend = designp->posend.begin();

  double * * linpredp;

  if (likep->linpred_current==1)
    {
    linpredp = designp->linpredp1.getV();
    }
  else
    {
    linpredp = designp->linpredp2.getV();
    }

  double * workintvar = designp->intvar.getV();

  double * betap = beta.getV();
  double * betaoldp = betaold.getV();

  for (j=0;j<size;j++,++itbeg,++itend,betap++,betaoldp++)
    {
    for (i=*itbeg;i<=*itend;i++,linpredp++,workintvar++)
      {
      *(*linpredp) =  linpred_o(designp->index_data(i,0),0) + exp((*betap)  * (*workintvar))
                      - exp(*(betaoldp)  * (*workintvar));

      }
    }
   */
  }


bool FC_hrandom::posteriormode_multexp(void)
  {
/*
  if (response_o.rows()==1)
    {
    response_o = likep->response;
    linpred_o = datamatrix(response_o.rows(),1);
    }
  // intvar = log(f)
  // linpred = etarest + exp(random_effect)* intvar

  likep->optionbool1 = true;
//  likep->changingworkingweights = true;
  likep->updateIWLS = true;

  update_response_multexp();     // linpred = random_effect*intvar

  bool h = posteriormode_additive();

  update_linpred_multexp();

  likep->optionbool1 = false;
//  likep->changingworkingweights = false;
  likep->updateIWLS = false;

  likep->response.assign(response_o);
 */

  return true;

  }


bool FC_hrandom::posteriormode_additive(void)
  {

  bool conv;
  conv= FC_nonp::posteriormode();

  set_rcoeff();

  bool conv2 = FCrcoeff.posteriormode();

  likep_RE->workingresponse.assign(beta);
  likep_RE->response.assign(beta);


  // TEST
  /*
  ofstream out5("c:\\bayesx\\test\\results\\fhrandom.res");
  beta.prettyPrint(out5);
  */
  // TEST

  return conv;
  }


bool FC_hrandom::posteriormode(void)
  {

  if (rtype==multexp)
    {
    return posteriormode_multexp();
    }
  else
    {
    return posteriormode_additive();
    }
  }


void FC_hrandom::compute_autocorr_all(const ST::string & path,
                                      unsigned lag, ofstream & outg) const
  {
  FC::compute_autocorr_all(path,lag,outg);
  ST::string path2 = path.substr(0,path.length()-4) + "2.raw";

  FCrcoeff.compute_autocorr_all(path2,lag,outg);
  }


void FC_hrandom::get_samples(const ST::string & filename,ofstream & outg) const
  {
  FC::get_samples(filename,outg);
  ST::string path2 = filename.substr(0,filename.length()-4) + "2.raw";
  FCrcoeff.get_samples(path2,outg);
  }


void FC_hrandom::outgraphs(ofstream & out_stata, ofstream & out_R,
                          const ST::string & path)
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


  out_stata << "clear" << endl
            << "infile intnr " << xvar
            << " pmean_tot pqu" << pu1_str << "_tot "
            << " pqu" << po1_str << "_tot" << " pmed_tot pqu" << po2_str
            << "_tot" << " pqu" << pu2_str << "_tot"
            << " pcat" << pu_str << "_tot" << " pcat" << po_str << "_tot pqu"
            << pu1_str << "sim_tot "
            << " pqu" << po1_str << "sim_tot"
            << " pqu" << po2_str << "sim_tot"
            << " pqu" << pu2_str << "sim_tot"
            << " pcat" << pu_str << "sim_tot" << " pcat" << po_str << "sim_tot"
            << " pmean pqu" << pu1_str
            << " pqu" << po1_str << " pmed pqu" << po2_str << " pqu" <<
            pu2_str << " pcat" << pu_str  << " pcat" << po_str
            << " pqu" << pu1_str
            << "_sim pqu" << po1_str << "_sim pqu" << po2_str << "_sim pqu" <<
            pu2_str << "_sim pcat" << pu_str  << "_sim pcat" << po_str << "_sim";



  if (computemeaneffect==true)
    {
    out_stata << " pmean_mu pqu"
              << pu1_str << "_mu"
              << " pqu" << po1_str << "_mu"
              << " pmed_d pqu" << po2_str << "_mu"
              << " pqu" << pu2_str << "_mu";
    }


  out_stata << " using "
            << path << endl
            << "drop in 1" << endl;

  out_stata << "kdensity pmean_tot" << endl
            << "graph export " << pathps << "_tot.eps, replace"
            << endl << endl;

  out_stata << "kdensity pmean" << endl
            << "graph export " << pathps << ".eps, replace"
            << endl << endl;

  if (computemeaneffect==true)
    {
    out_stata << "kdensity pmean_mu" << endl
              << "graph export " << pathps << "_mu.eps, replace"
              << endl << endl;
    }

  }


void FC_hrandom::outresults(ofstream & out_stata,ofstream & out_R,
                            const ST::string & pathresults)
  {

  if (pathresults.isvalidfile() != 1)
    {

    ST::string pathbasis = pathresults.substr(0,pathresults.length()-4) +
                                 "_basisR.res";

    outbasis_R(pathbasis);

    ST::string paths = pathresults.substr(0,pathresults.length()-4) +
                                 "_sample.raw";

    out_R << "family=" << likep->family.strtochar() << ",";
    out_R << "hlevel=" << likep->hlevel << ",";
    out_R << "equationtype=" << likep->equationtype.strtochar() << ",";
    out_R << "term=sx("  << designp->datanames[0].strtochar()  << "),";
    out_R << "filetype=nonlinear,";
    out_R << "pathsamples=" << paths.strtochar() << ",";
    out_R << "pathbasis=" << pathbasis.strtochar() << ",";


    outgraphs(out_stata,out_R,pathresults);

    FC::outresults(out_stata,out_R,"");
    FCrcoeff.outresults(out_stata,out_R,"");

    outresults_acceptance();

   if (computemeaneffect==true)
      meaneffect_sample.outresults(out_stata,out_R,"");

    optionsp->out("    Results are stored in file\n");
    optionsp->out("    " +  pathresults + "\n");
    optionsp->out("\n");

//   ofstream out("c:\\bayesx\\trunk\\testh\\results\\datare.raw");
//   designp->data.prettyPrint(out);

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

    optionsp->out("    Mean effects evaluated at " +
                  designp->datanames[designp->datanames.size()-1] + "=" +
                  designp->effectvalues[designp->meaneffectnr]);


    double s_level1 = 0;
    double s_level2 = 0;

    double s_level1_rcoeff = 0;
    double s_level2_rcoeff = 0;

    if (optionsp->samplesize > 0)
      {
      s_level1 = simconfBand(true);
      s_level2 = simconfBand(false);

      s_level1_rcoeff = FCrcoeff.simconfBand(true);
      s_level2_rcoeff = FCrcoeff.simconfBand(false);


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
    outres << designp->datanames[designp->datanames.size()-1] << "   ";
    outres << "pmean_tot   ";

    outres << "pqu"  << l1  << "_tot   ";
    outres << "pqu"  << l2  << "_tot   ";
    outres << "pmed_tot   ";
    outres << "pqu"  << u1  << "_tot   ";
    outres << "pqu"  << u2  << "_tot   ";
    outres << "pcat" << optionsp->level1 << "_tot   ";
    outres << "pcat" << optionsp->level2 << "_tot   ";

    outres << "pqu"  << l1  << "tot_sim   ";
    outres << "pqu"  << l2  << "tot_sim   ";
    outres << "pqu"  << u1  << "tot_sim   ";
    outres << "pqu"  << u2  << "tot_sim   ";
    outres << "pcat" << optionsp->level1 << "tot_sim   ";
    outres << "pcat" << optionsp->level2 << "tot_sim   ";

    outres << "pmean   ";

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

    if (computemeaneffect==true)
      {

      outres << "pmean_mu   ";

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

    double * workmean = betamean.getV();
    double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
    double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
    double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
    double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
    double * workbetaqu50 = betaqu50.getV();

    double * workmean_rcoeff = FCrcoeff.betamean.getV();
    double * workbetaqu_l1_lower_p_rcoeff = FCrcoeff.betaqu_l1_lower.getV();
    double * workbetaqu_l2_lower_p_rcoeff = FCrcoeff.betaqu_l2_lower.getV();
    double * workbetaqu_l1_upper_p_rcoeff = FCrcoeff.betaqu_l1_upper.getV();
    double * workbetaqu_l2_upper_p_rcoeff = FCrcoeff.betaqu_l2_upper.getV();
    double * workbetaqu50_rcoeff = FCrcoeff.betaqu50.getV();

    double * mu_workmean=NULL;
    double * mu_workbetaqu_l1_lower_p=NULL;
    double * mu_workbetaqu_l2_lower_p=NULL;
    double * mu_workbetaqu_l1_upper_p=NULL;
    double * mu_workbetaqu_l2_upper_p=NULL;
    double * mu_workbetaqu50=NULL;

    if (computemeaneffect==true)
      {
      mu_workmean = meaneffect_sample.betamean.getV();
      mu_workbetaqu_l1_lower_p = meaneffect_sample.betaqu_l1_lower.getV();
      mu_workbetaqu_l2_lower_p = meaneffect_sample.betaqu_l2_lower.getV();
      mu_workbetaqu_l1_upper_p = meaneffect_sample.betaqu_l1_upper.getV();
      mu_workbetaqu_l2_upper_p = meaneffect_sample.betaqu_l2_upper.getV();
      mu_workbetaqu50 = meaneffect_sample.betaqu50.getV();
      }

    double l1_sim,l2_sim,u1_sim,u2_sim;

    unsigned nrpar = beta.rows();
    for(i=0;i<nrpar;i++,workmean++,workbetaqu_l1_lower_p++,
                              workbetaqu_l2_lower_p++,workbetaqu50++,
                              workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++,
                              workmean_rcoeff++,workbetaqu_l1_lower_p_rcoeff++,
                              workbetaqu_l2_lower_p_rcoeff++,
                              workbetaqu_l1_upper_p_rcoeff++,
                              workbetaqu_l2_upper_p_rcoeff++,
                              workbetaqu50_rcoeff++)
      {
      outres << (i+1) << "   ";
      outres << designp->effectvalues[i] << "   ";
      outres << *workmean << "   ";

      if (optionsp->samplesize > 1)
        {
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
        outres << "0   0   0   0   0   0   0   0   0   0   0   0   0   ";

      outres << *workmean_rcoeff << "   ";

      if (optionsp->samplesize > 1)
        {
        outres << *workbetaqu_l1_lower_p_rcoeff << "   ";
        outres << *workbetaqu_l2_lower_p_rcoeff << "   ";
        outres << *workbetaqu50_rcoeff << "   ";
        outres << *workbetaqu_l2_upper_p_rcoeff << "   ";
        outres << *workbetaqu_l1_upper_p_rcoeff << "   ";

        if (*workbetaqu_l1_lower_p_rcoeff > 0)
          outres << 1 << "   ";
        else if (*workbetaqu_l1_upper_p_rcoeff < 0)
          outres << -1 << "   ";
        else
          outres << 0 << "   ";

        if (*workbetaqu_l2_lower_p_rcoeff > 0)
          outres << 1 << "   ";
        else if (*workbetaqu_l2_upper_p_rcoeff < 0)
          outres << -1 << "   ";
        else
          outres << 0 << "   ";

        l1_sim = *workmean_rcoeff - s_level1_rcoeff*(*workmean_rcoeff- *workbetaqu_l1_lower_p_rcoeff);
        l2_sim = *workmean_rcoeff - s_level2_rcoeff*(*workmean_rcoeff- *workbetaqu_l2_lower_p_rcoeff);
        u1_sim = *workmean_rcoeff + s_level1_rcoeff*(*workbetaqu_l1_upper_p_rcoeff - *workmean_rcoeff);
        u2_sim = *workmean_rcoeff + s_level2_rcoeff*(*workbetaqu_l2_upper_p_rcoeff - *workmean_rcoeff);

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
        outres << "0   0   0   0   0   0   0   0   0   0   0   0   0   ";


      if (computemeaneffect==true)
        {

        outres << *mu_workmean << "   ";

        if (optionsp->samplesize > 1)
          {
          outres << *mu_workbetaqu_l1_lower_p << "   ";
          outres << *mu_workbetaqu_l2_lower_p << "   ";
          outres << *mu_workbetaqu50 << "   ";
          outres << *mu_workbetaqu_l2_upper_p << "   ";
          outres << *mu_workbetaqu_l1_upper_p << "   ";
          }
        else
          outres << "0   0   0   0   0   ";

        if (i <nrpar-1)
          {
          mu_workmean++;
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

} // end: namespace MCMC




