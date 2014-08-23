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



#include "FC_linear.h"


//------------------------------------------------------------------------------
//----------------- CLASS: FC implementation of member functions ---------------
//------------------------------------------------------------------------------


namespace MCMC
{


void FC_linear::read_options(vector<ST::string> & op,vector<ST::string> & vn)
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
  */

  }


FC_linear::FC_linear(void)
  {
  }


int FC_linear::add_variable(const datamatrix & d,ST::string & name)
  {
  datanames.push_back(name);
  designhelp.push_back(d);
  return designhelp.size()-1;
  }


FC_linear::FC_linear(MASTER_OBJ * mp,unsigned & enr,GENERAL_OPTIONS * o,DISTR * lp,
                    datamatrix & d,
                 vector<ST::string> & vn, const ST::string & t,
                 const ST::string & fp,bool cent)
     : FC(o,t,1,1,fp)
  {

  masterp = mp;
  equationnr = enr;
  likep = lp;
  unsigned i;
  datanames = vn;
  if (datanames.size() > 0)
    {
    for (i=0;i<d.cols();i++)
      designhelp.push_back(d.getCol(i));
    }
  initialize = false;
  IWLS = likep->updateIWLS;

  center = cent;
  rankXWX_ok = true;
  constwarning=false;
  }


FC_linear::FC_linear(const FC_linear & m)
  : FC(FC(m))
  {
  constwarning=m.constwarning;
  constposition = m.constposition;
  masterp = m.masterp;
  equationnr = m.equationnr;
  IWLS = m.IWLS;
  likep = m.likep;
  design = m.design;
  designhelp = m.designhelp;
  meaneffectdesign = m.meaneffectdesign;
  XWX = m.XWX;
  rankXWX_ok = m.rankXWX_ok;
  XWXold = m.XWXold;
  XWXroot = m.XWXroot;
  Xt = m.Xt;
  initialize = m.initialize;
  residual = m.residual;
  Xtresidual = m.Xtresidual;
  betaold = m.betaold;
  betadiff = m.betadiff;
  betam = m.betam;
  mode = m.mode;
  help = m.help;
  linold = m.linold;
  linnew = m.linnew;
  linmode = m.linmode;
  proposal=m.proposal;
  diff = m.diff;
  linoldp = m.linoldp;
  linnewp = m.linnewp;
  datanames = m.datanames;
  mean_designcols = m.mean_designcols;
  center = m.center;
  }


const FC_linear & FC_linear::operator=(const FC_linear & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  constwarning=m.constwarning;
  constposition = m.constposition;
  masterp = m.masterp;
  equationnr = m.equationnr;
  IWLS = m.IWLS;
  likep = m.likep;
  design = m.design;
  designhelp = m.designhelp;
  meaneffectdesign = m.meaneffectdesign;
  XWX = m.XWX;
  rankXWX_ok = m.rankXWX_ok;
  XWXold = m.XWXold;
  XWXroot = m.XWXroot;
  Xt = m.Xt;
  initialize = m.initialize;
  residual = m.residual;
  Xtresidual = m.Xtresidual;
  betaold = m.betaold;
  betadiff = m.betadiff;
  betam = m.betam;
  mode = m.mode;
  proposal=m.proposal;
  help = m.help;
  linold = m.linold;
  linnew = m.linnew;
  linmode = m.linmode;
  diff = m.diff;
  linoldp = m.linoldp;
  linnewp = m.linnewp;
  mean_designcols = m.mean_designcols;
  datanames = m.datanames;
  center = m.center;
  return *this;
  }


void FC_linear::add_linpred(datamatrix & l)
  {

  if (likep->linpred_current==1)
    likep->linearpred1.plus(l);
  else
    likep->linearpred2.plus(l);
  }



void FC_linear::update_IWLS(void)
  {


/*
      ofstream out3("c:\\bayesx\\testh\\results\\linpred.res");
      if (likep->linpred_current==1)
        likep->linearpred1.prettyPrint(out3);
      else
        likep->linearpred2.prettyPrint(out3);
*/

  double qoldbeta;
  double qnewbeta;

  if (!initialize)
    create_matrices();

  if (optionsp->nriter == 1)
    {
    linold.mult(design,beta);
    mode.assign(beta);
    }

    double logold = likep->loglikelihood(true);

    linmode.mult(design,mode);
    diff.minus(linmode,*linoldp);
    add_linpred(diff);

    double h = likep->compute_iwls(true,false);

    compute_XWXroot(XWXold);

    compute_Wpartres(linmode);
    Xtresidual.mult(Xt,residual);

    XWXroot.solveroot(Xtresidual,help,mode);

    help.minus(beta,mode);
    qoldbeta = -0.5*XWXold.compute_quadform(help);


    unsigned i;
    double * workh = help.getV();
    for(i=0;i<help.rows();i++,workh++)
      *workh = rand_normal();

    XWXroot.solveroot_t(help,proposal);

    proposal.plus(mode);

    help.minus(proposal,mode);

    qnewbeta = -0.5*XWXold.compute_quadform(help);

    linnewp->mult(design,proposal);

    diff.minus(*linnewp,linmode);

    add_linpred(diff);                           // (mit proposed)

    bool ok;
    if (optionsp->saveestimation)
      {
      ok = likep->check_linpred();
      if (!ok)
        outsidelinpredlimits++;
      }
    else
      ok = true;

    double logprop=0;
    if (ok)
      logprop = likep->loglikelihood();     // mit proposed

    double u = log(uniform());

    if (ok && (u <= (logprop + qoldbeta - logold - qnewbeta)) )
      {
      datamatrix * mp = linoldp;
      linoldp = linnewp;
      linnewp = mp;

      beta.assign(proposal);

      acceptance++;
      }
    else
      {
      diff.minus(*linoldp,*linnewp);
      add_linpred(diff);
      }


/*
      ofstream out4("c:\\bayesx\\testh\\results\\linprednew.res");
      if (likep->linpred_current==1)
        likep->linearpred1.prettyPrint(out4);
      else
        likep->linearpred2.prettyPrint(out4);
*/

    FC::update();

//    }


  }

void FC_linear::update(void)
  {
  if ((datanames.size() > 0) && (rankXWX_ok==true))
    {
    if (IWLS)
      update_IWLS();
    else
      update_gaussian();

    masterp->level1_likep[equationnr]->meaneffect -= meaneffect;
    meaneffect = (meaneffectdesign*beta)(0,0);
    masterp->level1_likep[equationnr]->meaneffect += meaneffect;

    }
  else
    nosamples = true;
  }


void FC_linear::update_gaussian(void)
  {
  if ((datanames.size() > 0) && (rankXWX_ok==true))
    {

    if (!initialize)
      create_matrices();

    compute_XWXroot(XWX);

    linold.mult(design,beta);
    compute_Wpartres(linold);
    Xtresidual.mult(Xt,residual);

    XWXroot.solveroot(Xtresidual,help,betam);

    double sigmaresp = sqrt(likep->get_scale());
    unsigned i;
    double * workh = help.getV();
    for(i=0;i<help.rows();i++,workh++)
      *workh = sigmaresp*rand_normal();

    XWXroot.solveroot_t(help,beta);
    beta.plus(betam);

    betadiff.minus(beta,betaold);

    if (likep->linpred_current==1)
      likep->linearpred1.addmult(design,betadiff);
    else
      likep->linearpred2.addmult(design,betadiff);

    bool ok;
    if (optionsp->saveestimation)
      {
      ok = likep->check_linpred();
      if (!ok)
        outsidelinpredlimits++;
      }
    else
      ok = true;

    if (ok)
      {
      betaold.assign(beta);

      acceptance++;
      }
    else
      {
      betadiff.minus(betaold,beta);

      if (likep->linpred_current==1)
        likep->linearpred1.addmult(design,betadiff);
      else
        likep->linearpred2.addmult(design,betadiff);

      beta.assign(betaold);
      }

    FC::update();
    }
  }


void FC_linear::compute_XWXroot(datamatrix & r)
  {

  compute_XWX(r);

  if ((likep->wtype==wweightschange_weightsneqone) ||
      (likep->wtype==wweightschange_weightsone) ||
      (optionsp->nriter<=1))
    {
    XWXroot = r.root();
    }

  }


void FC_linear::compute_XWX(datamatrix & r)
  {

  if ((likep->wtype==wweightschange_weightsneqone) ||
      (likep->wtype==wweightschange_weightsone) ||
      (optionsp->nriter<=1))
    {

    unsigned i,j,k;
    unsigned nrconst = beta.rows();
    unsigned nrobs = Xt.cols();
    double * Xt_ip;
    double * Xt_jp;
    double * workingweightp;
    double help;

    if (likep->wtype==wweightsnochange_one)
      {
      for (i=0;i<nrconst;i++)
        for (j=i;j<nrconst;j++)
          {
          help = 0;
          Xt_ip = Xt.getV()+i*nrobs;
          Xt_jp = Xt.getV()+j*nrobs;

          for (k=0;k<nrobs;k++,Xt_ip++,Xt_jp++)
            help += (*Xt_ip)*(*Xt_jp);

          r(i,j) = help;
          if (i!=j)
            r(j,i) = help;
          }
      }
    else
      {
      for (i=0;i<nrconst;i++)
        for (j=i;j<nrconst;j++)
          {
          help = 0;
          Xt_ip = Xt.getV()+i*nrobs;
          Xt_jp = Xt.getV()+j*nrobs;
          workingweightp = likep->workingweight.getV();

          for (k=0;k<nrobs;k++,Xt_ip++,Xt_jp++,workingweightp++)
            help += (*workingweightp) * (*Xt_ip)*(*Xt_jp);

          r(i,j) = help;
          if (i!=j)
            r(j,i) = help;
          }
      }

    }

  }


void FC_linear::compute_meaneffect_design(void)
  {
  unsigned i,j;

  meaneffectdesign = datamatrix(1,design.cols(),0);

  double  mhelp;

  double bestdiff;
  double currentdiff;

  for (j=0;j<design.cols();j++)
    {
    mhelp = design.mean(j);
    bestdiff = fabs(design(0,j) - mhelp);
    meaneffectdesign(0,j) = design(0,j);
    for (i=1;i<design.rows();i++)
      {
      currentdiff = fabs(design(i,j) - mhelp);
      if (currentdiff < bestdiff)
        {
        bestdiff = currentdiff;
        meaneffectdesign(0,j) = design(i,j);
        }

      }

    }

  }


void FC_linear::find_const(datamatrix & design)
  {
  constposition = -1;
  bool constfound = false;
  unsigned i=0;
  unsigned j;
  while (constfound==false && i < design.cols())
    {
    j = 0;
    bool allone = true;
    while (allone == true && j < design.rows())
      {
      if (design(j,i) != 1)
        allone = false;
      j++;
      }
    if (allone == true)
      {
      constfound = true;
      constposition = i;
      }
    i++;
    }

  if (constposition==-1)
    {
    optionsp->out("\n");
    optionsp->out("WARNING: AT LEAST ONE EQUATION CONTAINS NO INTERCEPT\n");
    optionsp->out("         Intercept may be specified using const in linear effects term\n");
    optionsp->out("\n");
    }

  }


void FC_linear::create_matrices(void)
  {

  unsigned i,j;
  design = datamatrix(designhelp[0].rows(),designhelp.size());
  for(i=0;i<designhelp.size();i++)
    design.putCol(i,designhelp[i]);

  find_const(design);

  if (center == true)
    {
    double m;
    int i;
    mean_designcols = datamatrix(1,design.cols(),1);
    for (i=0;i<design.cols();i++)
      {
      if (i!= constposition)
        {
        m = design.mean(i);
        for (j=0;j<design.rows();j++)
          design(j,i) -= m;
        mean_designcols(0,i) = m;
        }
      }

    }

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\design.res");
  // design.prettyPrint(out);
  // TEST

  compute_meaneffect_design();

  // TEST
  // ofstream out("c:\\bayesx\\testh\\results\\meandesign.res");
  // meaneffectdesign.prettyPrint(out);
  // TEST

  Xt = design.transposed();
  XWX = datamatrix(design.cols(),design.cols(),0);


  residual = datamatrix(design.rows(),1,0);
  Xtresidual = datamatrix(design.cols(),1,0);

  setbeta(design.cols(),1,0);
  betaold=datamatrix(beta.rows(),1,0);

  /*
  if (constposition != -1)
    {
    double m = likep->get_intercept_start();
    beta(constposition,0) = m;
    betaold(constposition,0) = m;
    double * linpred;
    if (likep->linpred_current==1)
      linpred = likep->linearpred1.getV();
    else
      linpred = likep->linearpred2.getV();
    unsigned i;
    for (i=0;i<likep->nrobs;i++,linpred++)
      *linpred += m;
    }
   */

  betadiff = betaold;
  betam = beta;
  help = beta;
  linold = datamatrix(design.rows(),1,0);
  initialize=true;


  linnew = datamatrix(design.rows(),1,0);
  linmode = datamatrix(design.rows(),1,0);
  diff = datamatrix(design.rows(),1,0);
  linnewp = &linnew;
  linoldp = &linold;
  mode = beta;
  proposal = beta;
  XWXold = datamatrix(design.cols(),design.cols(),0);

  }


void FC_linear::compute_Wpartres(datamatrix & linpred)
  {
  unsigned i;
  double * workingweightp = likep->workingweight.getV();
  double * workingresponsep = likep->workingresponse.getV();
  double * residualp = residual.getV();
  double * linpredp = linpred.getV();

  double * predictorp;
  if (likep->linpred_current==1)
    predictorp = likep->linearpred1.getV();
  else
    predictorp = likep->linearpred2.getV();

  if (likep->wtype == wweightsnochange_one)
    {
    for (i=0;i<likep->nrobs;i++,workingresponsep++,
                                residualp++,linpredp++,predictorp++)
      *residualp = ((*workingresponsep)  - (*predictorp)
                   + (*linpredp));
    }
  else
    {
    for (i=0;i<likep->nrobs;i++,workingweightp++,workingresponsep++,
                                residualp++,linpredp++,predictorp++)
      if (*workingweightp==0)
        *residualp==0;
      else
        *residualp = *workingweightp * ((*workingresponsep)  - (*predictorp)
                     + (*linpredp));
    }
  }


double FC_linear::compute_XtWpartres(double & mo)
  {

  unsigned i;
  double * workingweightp = likep->workingweight.getV();
  double * workingresponsep = likep->workingresponse.getV();
  double res=0;


  double * predictorp;
  if (likep->linpred_current==1)
    predictorp = likep->linearpred1.getV();
  else
    predictorp = likep->linearpred2.getV();

  double * workdesign = design.getV();

  for (i=0;i<likep->nrobs;i++,workingweightp++,workingresponsep++,
                              predictorp++,workdesign++)
    res += *workdesign * (*workingweightp) *
           ((*workingresponsep)  - (*predictorp) + *workdesign*mo);

  return res;
  }



bool FC_linear::posteriormode(void)
  {

//  ofstream out("d:\\_sicher\\papzip\\resultsconst\\linpredvorher.raw");
//  likep->linearpred1.prettyPrint(out);

  if (datanames.size() > 0)
    {
    if (rankXWX_ok == true)
      {
      if (!initialize)
        create_matrices();

      double h = likep->compute_iwls(true,false);

      compute_XWX(XWX);
      datamatrix test = XWX.cinverse();
      if (test.rows() < XWX.rows())
        {
        rankXWX_ok = false;
        optionsp->out("    WARNING: Cross product matrix for linear effects is rank deficient\n");
        optionsp->out("             linear effects are not estimated\n");
        optionsp->out("\n");
        }
      }

    if (rankXWX_ok == true)
      {
      linold.mult(design,beta);
      compute_Wpartres(linold);

      Xtresidual.mult(Xt,residual);

      beta = XWX.solve(Xtresidual);

      betadiff.minus(beta,betaold);

      if (likep->linpred_current==1)
        likep->linearpred1.addmult(design,betadiff);
      else
        likep->linearpred2.addmult(design,betadiff);


      bool ok;
      if (optionsp->saveestimation)
        ok = likep->check_linpred();
      else
        ok = true;

      if (ok)
        {
        betaold.assign(beta);

        masterp->level1_likep[equationnr]->meaneffect -= meaneffect;
        meaneffect = (meaneffectdesign*beta)(0,0);
        masterp->level1_likep[equationnr]->meaneffect += meaneffect;
        }
      else
        {
        betadiff.minus(betaold,beta);

        if (likep->linpred_current==1)
          likep->linearpred1.addmult(design,betadiff);
        else
          likep->linearpred2.addmult(design,betadiff);

        beta.assign(betaold);
        }

      return FC::posteriormode();
      }
    }
  else
    {
    if (constwarning==false)
      {
      constwarning=true;
      optionsp->out("\n");
      optionsp->out("WARNING: AT LEAST ONE EQUATION CONTAINS NO INTERCEPT\n");
      optionsp->out("         Intercept may be specified using const in linear effects term\n");
      optionsp->out("\n");
      }
    }

  return true;
  }


void FC_linear::compute_autocorr_all(const ST::string & path,
                              unsigned lag, ofstream & outg) const
  {
  if ((datanames.size() > 0) && (rankXWX_ok==true))
    {
    FC::compute_autocorr_all(path,lag,outg);
    }
  }



void FC_linear::outoptions(void)
  {
//  optionsp->out("  OPTIONS FOR TERM: " + title + "\n",true);
//  optionsp->out("\n");
  }



void FC_linear::outresults(ofstream & out_stata,ofstream & out_R,
                           const ST::string & pathresults)
  {
  if ((datanames.size() > 0) && (rankXWX_ok==true))
    {

    FC::outresults(out_stata,out_R,pathresults);
    FC::outresults_help(out_stata,out_R,pathresults,datanames);
    FC::outresults_acceptance();

    optionsp->out("    Results for fixed effects are also stored in file\n");
    optionsp->out("    " + pathresults + "\n");

    ST::string paths = pathresults.substr(0,pathresults.length()-4) +
                                 "_sample.raw";

    out_R << "family=" << likep->familyshort.strtochar() << ",";
    out_R << "hlevel=" << likep->hlevel << ",";
    out_R << "equationtype=" << likep->equationtype.strtochar() << ",";
    out_R << "term=";
    unsigned k;
    for (k=0;k<datanames.size();k++)
      out_R << datanames[k].strtochar() << " ";
    out_R  << ",";
    out_R << "filetype=linear,";
    out_R << "pathsamples=" << paths.strtochar() << ",";
    out_R << "pathbasis=" << endl;

    if (center==true)
      {
      optionsp->out("\n");
      optionsp->out("    Note: Covariates with linear effects are centered around zero before estimation\n");
      optionsp->out("          Centering of covariates may improve the mixing of the MCMC sampler while\n");
      optionsp->out("          the regression coefficents are unchanged\n");
      optionsp->out("          However the intercept is changed due to the centering of covariates.\n");
      optionsp->out("          The means of the covariates are:\n");
      unsigned k;
      for (k=0;k<mean_designcols.cols();k++)
        {
        if (k != constposition)
          {
          optionsp->out("          " + datanames[k] + ": " + ST::doubletostring(mean_designcols(0,k),6) + "\n");
          }

        }

      } // end: if (center==true)

    optionsp->out("\n");

    }

  }


void FC_linear::reset(void)
  {

  }


//------------------------------------------------------------------------------
//------------------------------- FC_linear_pen --------------------------------
//------------------------------------------------------------------------------


FC_linear_pen::FC_linear_pen(void)
  {
  }




FC_linear_pen::FC_linear_pen(MASTER_OBJ * mp,unsigned & enr,
                            GENERAL_OPTIONS * o,DISTR * lp, datamatrix & d,
                 vector<ST::string> & vn, const ST::string & t,
                 const ST::string & fp,bool cent)
     : FC_linear(mp,enr,o,lp,d,vn,t,fp,cent)
  {


  }


FC_linear_pen::FC_linear_pen(const FC_linear_pen & m)
  : FC_linear(FC_linear(m))
  {
  tau2 = m.tau2;
  tau2oldinv = m.tau2oldinv;
  }


const FC_linear_pen & FC_linear_pen::operator=(const FC_linear_pen & m)
  {

  if (this==&m)
	 return *this;
  FC_linear::operator=(FC_linear(m));
  tau2 = m.tau2;
  tau2oldinv = m.tau2oldinv;
  return *this;
  }



void FC_linear_pen::update(void)
  {

  FC_linear::update();
  }



bool FC_linear_pen::posteriormode(void)
  {
  return FC_linear::posteriormode();
  }




void FC_linear_pen::outoptions(void)
  {
//  optionsp->out("  OPTIONS FOR TERM: " + title + "\n",true);
//  optionsp->out("\n");
  }

void FC_linear_pen::outresults(ofstream & out_stata,ofstream & out_R,
                           const ST::string & pathresults)
  {
  FC_linear::outresults(out_stata,out_R,pathresults);
  }


void FC_linear_pen::compute_XWXroot(datamatrix & r)
  {
  compute_XWX(r);
  XWXroot = r.root();
  }


void FC_linear_pen::compute_XWX(datamatrix & r)
  {

  unsigned i;
  double * tau2p = tau2.getV();
  double * tau2oldinvp = tau2oldinv.getV();
  unsigned nrpar = beta.rows();

  if ((likep->wtype==wweightschange_weightsneqone) ||
      (likep->wtype==wweightschange_weightsone) ||
      (optionsp->nriter<=1))
    {
    FC_linear::compute_XWX(r);
    for (i=0;i<nrpar;i++,tau2p++,tau2oldinvp++)
      {
      XWX(i,i) += (1/(*tau2p));
      *tau2oldinvp =  1/(*tau2p);
      }

    }
  else
    {

    for (i=0;i<nrpar;i++,tau2p++,tau2oldinvp++)
      {
      XWX(i,i) += (1/(*tau2p) - *tau2oldinvp);
      *tau2oldinvp =  1/(*tau2p);
      }

    }

  }


void FC_linear_pen::reset(void)
  {

  }




} // end: namespace MCMC



