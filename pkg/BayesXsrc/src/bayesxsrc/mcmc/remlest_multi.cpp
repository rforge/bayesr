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



#include "remlest_multi.h"

#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"

#endif

using std::ofstream;
using std::flush;

//------------------------------------------------------------------------------
//------------------------ CLASS: remlest_multinomial --------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------------- Constructor ------------------------------------
//------------------------------------------------------------------------------

remlest_multinomial::remlest_multinomial(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adb,
#endif
vector<MCMC::FULLCOND*> & fc,datamatrix & re,
                const ST::string & family, const ST::string & ofile,
                const int & maxiter, const double & lowerlimit,
                const double & epsi, const double & maxch, const double & maxv,
                const datamatrix & categories,
                const datamatrix & weight, const bool & fi, ostream * lo)
  {

  nrcat2=categories.rows();
  nrcat=nrcat2+1;
  cats=categories;

  nrobs=re.rows();
  nrobspos=nrobs;
  for(int i=0; i<nrobs; i++)
    {
    if(weight(i,0)==0)
      {
      nrobspos--;
      }
    }

  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = adb;
  #endif

  fisher=fi;

  logout = lo;
  respfamily=family;
  outfile=ofile;

  maxit=maxiter;
  lowerlim=lowerlimit;
  eps=epsi;
  maxchange=maxch;
  maxvar = maxv;

  fullcond = fc;
  unsigned i,j;

  xcut.push_back(0);
  xcutbeta.push_back(0);
  zcut.push_back(0);
  zcutbeta.push_back(0);

  for(i=0;i<fullcond.size();i++)
    {
    xcut.push_back(xcut[i]+fullcond[i]->get_dimX());
    if (i>0)
      {
      zcut.push_back(zcut[i-1]+fullcond[i]->get_dimZ());
      }
    }

  for(j=0; j<nrcat2; j++)
    {
    for(i=0;i<fullcond.size();i++)
      {
      xcutbeta.push_back(xcutbeta[xcutbeta.size()-1]+fullcond[i]->get_dimX());
      if (i>0)
        {
        zcutbeta.push_back(zcutbeta[zcutbeta.size()-1]+fullcond[i]->get_dimZ());
        }
      }
    }


  X = datamatrix(re.rows(),xcut[xcut.size()-1],0);
  Z = datamatrix(re.rows(),zcut[zcut.size()-1],0);

  fullcond[0]->createreml(X,Z,xcut[0],0);

  for(i=1;i<fullcond.size();i++)
    {
    fullcond[i]->createreml(X,Z,xcut[i],zcut[i-1]);
    }

  partialnrpar=X.cols()+Z.cols();
  partialnrfixed=X.cols();
  totalnrfixed=X.cols()*nrcat2;
  partialnrrandom=Z.cols();
  totalnrpar=partialnrpar*nrcat2;
  partialvar=zcut.size()-1;

  beta=statmatrix<double>(partialnrpar*(nrcat2),1,0);

  theta=statmatrix<double>(zcutbeta.size()-1,1,0);

  for(j=0; j<nrcat2; j++)
    {
    for(i=1; i<fullcond.size(); i++)
      {
      theta(j*(fullcond.size()-1)+i-1,0) = fullcond[i]->get_startlambda();
      }
    }
  }

//------------------------------------------------------------------------------
//----------------------------- Estimation -------------------------------------
//------------------------------------------------------------------------------

bool remlest_multinomial::estimate(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight)
  {
  unsigned i, j, k, l;

  outoptions();
  out("\n");

  for(i=0;i<fullcond.size();i++)
    fullcond[i]->outoptionsreml();

  out("\n");
  out("REML ESTIMATION STARTED\n",true);
  out("\n");

  bool stop = check_pause();
  if (stop)
    return true;

  double help;
  datamatrix helpmat(nrcat2,1,0);
  datamatrix helpmat2(resp.rows(),1,0);

  // Matrix to store old versions of beta and theta
  statmatrix<double>betaold(beta.rows(),1,0);
  statmatrix<double>thetaold(theta.rows(),1,0);

  // Score-function and expected Fisher information for theta
  statmatrix<double>score(theta.rows(),1,0);
  statmatrix<double>Fisher(theta.rows(),theta.rows(),0);

  // Number of iterations
  unsigned it=1;

  // Criteria to detemine convergence
  double crit1=1;                //relative changes in regression parameters
  double crit2=1;                //relative changes in variance parameters
  bool test=true;

  vector<double>stopcrit(theta.rows(),10);
  vector<int>its(theta.rows(),0);
  vector<int>thetastop(theta.rows(),0);
  vector<int>signs(theta.rows(),1);

  // Linear predictor and indicator response
  statmatrix<double>respind((nrcat2)*resp.rows(),1,0);
  statmatrix<double>eta(respind.rows(),1,0);
  compute_respind(resp,respind);

  // Working observations and weights
  statmatrix<double>worky(respind.rows(),1,0);
  statmatrix<double>workweight(respind.rows(),nrcat2,0);
  statmatrix<double>mu(respind.rows(),1,0);

  // Matrix containing the inverse covariance matrix of the random effects
  statmatrix<double>Qinv(Z.cols()*nrcat2,1,0);

  // Matrices for Fisher scoring (regression parameters)
  statmatrix<double>H(beta.rows(),beta.rows(),0);
  statmatrix<double>H1(beta.rows(),1,0);

  // Matrices for Fisher scoring (variance parameters)
  statmatrix<double>Hinv(beta.rows(),beta.rows(),0);
  statmatrix<double>wresid(nrcat2,resp.rows(),0);                  //matrix containing row vectors !!

  // Transform smoothing paramater starting values to variances
  for(i=0; i<theta.rows(); i++)
    {
    theta(i,0)=1/theta(i,0);
    }

  while(test==true)
    {

    // store current values in betaold and thetaold and compute Qinv
    betaold=beta;
    thetaold=theta;
    for(i=0, l=0; i<theta.rows(); i++)
      {
      for(k=zcutbeta[i]; k<zcutbeta[i+1]; k++, l++)
        {
        Qinv(l,0)=1/theta(i,0);
        }
      }

    compute_eta2(eta);

    compute_weights(mu,workweight,worky,eta,respind,weight);
    compute_sscp2(H,workweight);

    stop = check_pause();
    if (stop)
      return true;

    compute_sscp_resp2(H1,workweight,worky);
    H.addtodiag(Qinv,totalnrfixed,totalnrpar);

    // Fisher-Scoring for beta
    beta=H.solve(H1);

    // update linear predictor and compute residuals
    compute_eta2(eta);
    worky = worky - eta;
    for(i=0; i<resp.rows(); i++)
      {
      for(j=0; j<nrcat2; j++)
        {
        wresid(j,i)=(workweight.getRow(i*nrcat2+j)*worky.getRowBlock(i*nrcat2,(i+1)*nrcat2))(0,0);
        }
      }

    // transform theta
    for(i=0; i<theta.rows(); i++)
      {
      thetaold(i,0)=signs[i]*sqrt(thetaold(i,0));
      }

    Hinv=H.inverse();
    H.subfromdiag(Qinv,totalnrfixed,totalnrpar);

    stop = check_pause();
    if (stop)
      return true;

    // compute score-function and expected fisher information

    for(i=0; i<theta.rows(); i++)
      {
      score(i,0)=-1*(H.getBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1],totalnrfixed+zcutbeta[i+1])*thetaold(i,0)).trace()+
                 ((H.getRowBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1]))*Hinv*(H.getColBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1]))*thetaold(i,0)).trace();
      for(k=0; k<theta.rows(); k++)
        {
        Fisher(i,k)=2*(H.getBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[i+1],totalnrfixed+zcutbeta[k+1])*H.getBlock(totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[k+1],totalnrfixed+zcutbeta[i+1])*thetaold(i,0)*thetaold(k,0)).trace()-
                    4*(H.getRowBlock(totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[k+1])*Hinv*H.getColBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1])*H.getBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[i+1],totalnrfixed+zcutbeta[k+1])*thetaold(i,0)*thetaold(k,0)).trace()+
                    2*(H.getRowBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1])*Hinv*H.getColBlock(totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[k+1])*H.getRowBlock(totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[k+1])*Hinv*H.getColBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1])*thetaold(i,0)*thetaold(k,0)).trace();
        Fisher(k,i)=Fisher(i,k);
        }
      }

    for(j=0; j<nrcat2; j++)
      {
      for(i=0; i<partialvar; i++)
        {
        for(l=zcut[i]; l<zcut[i+1]; l++)
          {
          help = (wresid.getRow(j)*(Z.getCol(l)))(0,0);
          score(j*partialvar+i,0) += help*help*thetaold(j*partialvar+i,0);
          }
        }
      }

    // fisher scoring for theta
    theta = thetaold + Fisher.solve(score);

    // transform theta back to original parameterisation

    for(i=0; i<theta.rows(); i++)
      {
      signs[i] = -1*(theta(i,0)<0)+1*(theta(i,0)>=0);
      theta(i,0) *= theta(i,0);
      thetaold(i,0) *= thetaold(i,0);
      }

    // test whether to stop estimation of theta[i]

    //compute norm of eta for the different catetgories
    helpmat=datamatrix(nrcat2,1,0);
    for(i=0; i<resp.rows(); i++)
      {
      for(j=0; j<nrcat2; j++)
        {
        helpmat(j,0) += eta(i*nrcat2+j,0)*eta(i*nrcat2+j,0);
        }
      }
    for(j=0; j<nrcat2; j++)
      {
      helpmat(j,0) = sqrt(helpmat(j,0));
      }
    // compute norm of the random parts
    for(j=0; j<nrcat2; j++)
      {
      for(i=0; i<partialvar; i++)
        {
        helpmat2 = Z.getColBlock(zcut[i],zcut[i+1])*beta.getRowBlock(totalnrfixed+zcutbeta[j*partialvar+i],totalnrfixed+zcutbeta[j*partialvar+i+1]);
        stopcrit[j*partialvar+i]=helpmat2.norm(0)/helpmat(j,0);
        if(stopcrit[j*partialvar+i]<lowerlim || theta(j*partialvar+i,0)>maxvar)
          {
          theta(j*partialvar+i,0)=thetaold(j*partialvar+i,0);
          if(stopcrit[j*partialvar+i]<lowerlim)
            {
            thetastop[j*partialvar+i]=1;
            }
          else
            {
            thetastop[j*partialvar+i]=-1;
            }
          }
        else
          {
          its[j*partialvar+i]=it;
          }
        }
      }

    // compute convergence criteria
    help=betaold.norm(0);
    if(help==0)
      {
      help=0.00001;
      }
    betaold.minus(betaold,beta);
    crit1 = betaold.norm(0)/help;

    help=thetaold.norm(0);
    if(help==0)
      {
      help=0.00001;
      }
    thetaold.minus(thetaold,theta);
    crit2 = thetaold.norm(0)/help;

    stop = check_pause();
    if (stop)
      return true;

    out("  iteration "+ST::inttostring(it)+"\n");
    out("  relative changes in the regression coefficients: "+
         ST::doubletostring(crit1,6)+"\n");
    out("  relative changes in the variance parameters:     "+
         ST::doubletostring(crit2,6)+"\n");
    out("\n");

    // test criterion
    test=((crit1>eps) || (crit2>eps)) && (it<(unsigned)maxit);
    if(it>2)
      {
      test = test && (crit1<maxchange && crit2<maxchange);
      }

    // count iteration
    it=it+1;
    }

  if(crit1>=maxchange || crit2>=maxchange)
    {
    out("\n");
    outerror("ERROR: numerical problems due to large relative changes\n");
    outerror("       REML ESTIMATION DID NOT CONVERGE\n");
    out("\n");
    }
  else if(it>=(unsigned)maxit)
    {
    out("\n");
    outerror("WARNING: Number of iterations reached " + ST::inttostring(maxit) + "\n");
    outerror("         REML ESTIMATION DID NOT CONVERGE\n");
    out("\n");
    }
  else
    {
    out("\n");
    out("REML ESTIMATION CONVERGED\n",true);
    out("\n");
    }
  out("ESTIMATION RESULTS:\n",true);
  out("\n");

/*  ofstream outit((outfile+"_it.raw").strtochar());
  outit << it-1;
  outit.close();*/

  datamatrix thetareml(theta.rows(),4,0);
  thetareml.putCol(0,theta);
  datamatrix Hhelp = (H*Hinv);
  for(i=0; i<theta.rows(); i++)
    {
    thetareml(i,1)=thetastop[i];
    thetareml(i,2)=its[i];
    thetareml(i,3)=xcutbeta[i+2]-xcutbeta[i+1]+(Hhelp.getBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1],totalnrfixed+zcutbeta[i+1])).trace();
    }

// store inverse Fisher-Info and design matrices
  if(fisher)
    {
    ofstream outbeta((outfile+"_coef.raw").strtochar());
    beta.prettyPrint(outbeta);
    outbeta.close();
    ofstream outfisher((outfile+"_inversefisher.raw").strtochar());
    Hinv.prettyPrint(outfisher);
    outfisher.close();
    ofstream outx((outfile+"_fixeddesign.raw").strtochar());
    X.prettyPrint(outx);
    outx.close();
    ofstream outz((outfile+"_randomdesign.raw").strtochar());
    Z.prettyPrint(outz);
    outz.close();

    for(i=1;i<fullcond.size();i++)
      {
      fullcond[i]->outresultsgrid();
      }
    }

  for(j=0; j<nrcat2; j++)
    {
//    out("\n");
//    out("RESULTS FOR CATEGORY "+ST::doubletostring(cats(j,0),6)+":\n",true);
    out("\n");
    for(i=1; i<fullcond.size(); i++)
      {
      beta(j*partialnrfixed,0) += fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],j*partialvar+i-1,false,xcutbeta[j*fullcond.size()+i],totalnrfixed+zcutbeta[j*partialvar+i-1],cats(j,0),true,j*fullcond.size()+i);
      }
    beta(j*partialnrfixed,0) += fullcond[0]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[0],0,0,false,xcutbeta[j*fullcond.size()],0,cats(j,0),true,0);
    }

  loglike=0;
  aic=0;
  bic=0;
  gcv=0;
  df=(H*Hinv).trace();
  double refprob;

  for(i=0; i<resp.rows(); i++)
    {
    if(weight(i,0)>0)
      {
      k=0;
      refprob=0;
      for(j=0; j<nrcat2; j++)
        {
        if(respind(i*nrcat2+j,0)==1)
          {
          loglike += log(mu(i*nrcat2+j,0));
          k=1;
          }
        else
          {
          refprob += mu(i*nrcat2+j,0);
          }
        }
      if(k==0)
        {
        loglike += log(1-refprob);
        }
      }
    }
  loglike *= -2;
  gcv = loglike/(double)nrobspos*(1-(double)df/(double)nrobspos)*(1-(double)df/(double)nrobspos);
  aic = loglike + 2*df;
  bic = loglike + log(static_cast<double>(nrobspos))*df;

  out("\n");
  out("  Model Fit\n",true);
  out("\n");
  out("\n");
  out("  -2*log-likelihood:                 " + ST::doubletostring(loglike,6) + "\n");
  out("  Degrees of freedom:                " + ST::doubletostring(df,6) + "\n");
  out("  (conditional) AIC:                 " + ST::doubletostring(aic,6) + "\n");
  out("  (conditional) BIC:                 " + ST::doubletostring(bic,6) + "\n");
  out("  GCV (based on deviance residuals): " + ST::doubletostring(gcv,6) + "\n");
  out("\n");
  out("  Results on the model fit are stored in file\n");
  out("  "+outfile+"_modelfit.raw");
  out("\n");

  ofstream outfit((outfile+"_modelfit.raw").strtochar());
  outfit << "loglike df aic bic gcv" << endl;
  outfit << loglike << " " << df << " " << aic << " " << bic << " " << gcv << endl;
  outfit.close();

  out("\n");
  out("  Additive predictors and expectations\n",true);
  out("\n");
  out("\n");
  out("  Additive predictors and expectations for each observation\n");
  out("  and category are stored in file\n");
  out("  "+outfile+"_predict.raw\n");
  out("\n");
  out("\n");

  ofstream outpredict((outfile+"_predict.raw").strtochar());
  for(j=0; j< nrcat2; j++)
    {
    outpredict << "eta" << cats(j,0) << " ";
    outpredict << "mu" << cats(j,0) << " ";
    }
  outpredict << endl;
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      outpredict << eta(i*nrcat2+j,0) << " " << mu(i*nrcat2+j,0) << " ";
      }
    outpredict << endl;
    }
  outpredict.close();

  return false;
  }

bool remlest_multinomial::estimate_glm(const datamatrix resp,
                  const datamatrix & offset, const datamatrix & weight)
  {
  unsigned i,j;

  outoptions();
  out("\n");

  for(i=0;i<fullcond.size();i++)
    fullcond[i]->outoptionsreml();

  out("\n");
  out("REML ESTIMATION STARTED\n",true);
  out("\n");

  bool stop = check_pause();
  if (stop)
    return true;

  double help;

  // Matrix to store old version of beta
  statmatrix<double>betaold(beta.rows(),1,0);

  // Number of iterations
  unsigned it=1;

  // Criteria to detemine convergence
  double crit1=1;                //relative changes in regression parameters
  bool test=true;

  // Linear predictor and indicator response
  statmatrix<double>respind((nrcat2)*resp.rows(),1,0);
  statmatrix<double>eta(respind.rows(),1,0);
  compute_respind(resp,respind);

  // Working observations and weights
  statmatrix<double>worky(respind.rows(),1,0);
  statmatrix<double>workweight(respind.rows(),nrcat2,0);
  statmatrix<double>mu(respind.rows(),1,0);

  // Matrices for Fisher scoring (regression parameters)
  statmatrix<double>H(beta.rows(),beta.rows(),0);
  statmatrix<double>H1(beta.rows(),1,0);

  // Estimation loop
  while(test==true)
    {
    // store current values in betaold
    betaold=beta;

    compute_eta(eta);

    compute_weights(mu,workweight,worky,eta,respind,weight);
    compute_sscp(H,workweight);
    compute_sscp_resp(H1,workweight,worky);

    stop = check_pause();
    if (stop)
      return true;

    // Fisher-Scoring for beta
    beta=H.solve(H1);

    // compute convergence criteria
    help=betaold.norm(0);
    if(help==0)
      {
      help=0.00001;
      }
    betaold.minus(betaold,beta);
    crit1 = betaold.norm(0)/help;

    stop = check_pause();
    if (stop)
      return true;

    out("  iteration "+ST::inttostring(it)+"\n");
    out("  relative changes in the regression coefficients: "+
         ST::doubletostring(crit1,6)+"\n");
    out("\n");

    // test criterion
    test=(crit1>eps) && (it<(unsigned)maxit);
    if(it>2)
      {
      test = test && crit1<maxchange;
      }

    // count iteration
    it=it+1;

    }

  if(crit1>=maxchange)
    {
    out("\n");
    outerror("ERROR: numerical problems due to large relative changes\n");
    outerror("       REML ESTIMATION DID NOT CONVERGE\n");
    out("\n");
    }
  else if(it>=(unsigned)maxit)
    {
    out("\n");
    outerror("WARNING: Number of iterations reached " + ST::inttostring(maxit) + "\n");
    outerror("         REML ESTIMATION DID NOT CONVERGE\n");
    out("\n");
    }
  else
    {
    out("\n");
    out("REML ESTIMATION CONVERGED\n",true);
    out("\n");
    }
  out("ESTIMATION RESULTS:\n",true);
  out("\n");

  H=H.inverse();

  datamatrix helpmat = datamatrix(1,1,0);

// store inverse Fisher-Info and design matrices
  if(fisher)
    {
    ofstream outbeta((outfile+"_coef.raw").strtochar());
    beta.prettyPrint(outbeta);
    outbeta.close();
    ofstream outfisher((outfile+"_inversefisher.raw").strtochar());
    H.prettyPrint(outfisher);
    outfisher.close();
    ofstream outx((outfile+"_fixeddesign.raw").strtochar());
    X.prettyPrint(outx);
    outx.close();
    ofstream outz((outfile+"_randomdesign.raw").strtochar());
    Z.prettyPrint(outz);
    outz.close();

    for(i=1;i<fullcond.size();i++)
      {
      fullcond[i]->outresultsgrid();
      }
    }

  for(j=0; j<nrcat2; j++)
    {
//    out("\n");
//    out("RESULTS FOR CATEGORY "+ST::doubletostring(cats(j,0),6)+":\n",true);
    out("\n");
    for(i=1; i<fullcond.size(); i++)
      {
      beta(j*partialnrfixed,0) += fullcond[i]->outresultsreml(X,Z,beta,H,helpmat,xcut[i],zcut[i-1],j*partialvar+i-1,false,xcutbeta[j*fullcond.size()+i],totalnrfixed+zcutbeta[j*partialvar+i-1],cats(j,0),true,j*fullcond.size()+i);
      }
    beta(j*partialnrfixed,0) += fullcond[0]->outresultsreml(X,Z,beta,H,helpmat,xcut[0],0,0,false,xcutbeta[j*fullcond.size()],0,cats(j,0),true,0);
    }

  loglike=0;
  aic=0;
  bic=0;
  gcv=0;
  df=beta.rows();
  double refprob;
  unsigned k;

  for(i=0; i<resp.rows(); i++)
    {
    if(weight(i,0)>0)
      {
      k=0;
      refprob=0;
      for(j=0; j<nrcat2; j++)
        {
        if(respind(i*nrcat2+j,0)==1)
          {
          loglike += log(mu(i*nrcat2+j,0));
          k=1;
          }
        else
          {
          refprob += mu(i*nrcat2+j,0);
          }
        }
      if(k==0)
        {
        loglike += log(1-refprob);
        }
      }
    }
  loglike *= -2;
  gcv = loglike/(double)nrobspos*(1-(double)df/(double)nrobspos)*(1-(double)df/(double)nrobspos);
  aic = loglike + 2*df;
  bic = loglike + log(static_cast<double>(nrobspos))*df;

  out("\n");
  out("  Model Fit\n",true);
  out("\n");
  out("\n");
  out("  -2*log-likelihood:                 " + ST::doubletostring(loglike,6) + "\n");
  out("  Degrees of freedom:                " + ST::doubletostring(df,6) + "\n");
  out("  (conditional) AIC:                 " + ST::doubletostring(aic,6) + "\n");
  out("  (conditional) BIC:                 " + ST::doubletostring(bic,6) + "\n");
  out("  GCV (based on deviance residuals): " + ST::doubletostring(gcv,6) + "\n");
  out("\n");
  out("  Results on the model fit are stored in file\n");
  out("  "+outfile+"_modelfit.raw");
  out("\n");

  ofstream outfit((outfile+"_modelfit.raw").strtochar());
  outfit << "loglike df aic bic gcv" << endl;
  outfit << loglike << " " << df << " " << aic << " " << bic << " " << gcv << endl;
  outfit.close();


  out("\n");
  out("  Linear predictors and expectations\n",true);
  out("\n");
  out("\n");
  out("  Linear predictors and expectations for each observation\n");
  out("  and category are stored in file\n");
  out("  "+outfile+"_predict.raw\n");
  out("\n");
  out("\n");

  ofstream outpredict((outfile+"_predict.raw").strtochar());
  for(j=0; j< nrcat2; j++)
    {
    outpredict << "eta" << cats(j,0) << " ";
    outpredict << "mu" << cats(j,0) << " ";
    }
  outpredict << endl;
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      outpredict << eta(i*nrcat2+j,0) << " " << mu(i*nrcat2+j,0) << " ";
      }
    outpredict << endl;
    }
  outpredict.close();

  return false;
  }

//------------------------------------------------------------------------------
//------------- Weights, expectation, linear predictor, etc --------------------
//------------------------------------------------------------------------------

void remlest_multinomial::compute_respind(const datamatrix & re, datamatrix & respind)
  {
  unsigned i,j;
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      if(re(i,0)==cats(j,0))
        {
        respind(i*nrcat2+j,0)=1;
        }
      }
    }
  }

void remlest_multinomial::compute_eta(datamatrix & eta)
  {
  unsigned i,j;
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      eta(i*nrcat2+j,0)=((X.getRow(i))*(beta.getRowBlock(j*partialnrpar,(j+1)*partialnrpar)))(0,0);
      }
    }
  }

void remlest_multinomial::compute_eta2(datamatrix & eta)
  {
  unsigned i,j;
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      eta(i*nrcat2+j,0)=((X.getRow(i))*(beta.getRowBlock(j*partialnrfixed,(j+1)*partialnrfixed)))(0,0)+
                        ((Z.getRow(i))*(beta.getRowBlock(totalnrfixed+j*partialnrrandom,totalnrfixed+(j+1)*partialnrrandom)))(0,0);
      }
    }
  }

void remlest_multinomial::compute_weights(datamatrix & mu, datamatrix & workweights,
                  datamatrix & worky, datamatrix & eta, datamatrix & respind,
                  const datamatrix & weight)
  {
  unsigned i,j,k1,k2,l;

// Compute mu
  datamatrix expos(nrcat2,1,0);
  double exposum;
  for(i=0; i<nrobs; i++)
    {
    exposum=0;
    for(j=0; j<nrcat2; j++)
      {
      expos(j,0)=exp(eta(i*nrcat2+j,0));
      exposum+=expos(j,0);
      }
    exposum+=1;
    for(j=0; j<nrcat2; j++)
      {
      mu(i*nrcat2+j,0)=expos(j,0)/exposum;
      }
    }

// Compute weights

  for(i=0; i<nrobs; i++)
    {
    if(weight(i,0)>0)
      {
      for(j=i*nrcat2,l=0; l<nrcat2; j++, l++)
        {
        workweights(j,l)=mu(j,0)*(1-mu(j,0));
        for(k1=j+1,k2=l+1; k2<nrcat2; k1++, k2++)
          {
          workweights(j,k2)=-mu(j,0)*mu(k1,0);
          workweights(k1,l)=workweights(j,k2);
          }
        }
      }
    }

// Compute worky;

  for(i=0; i<nrobs; i++)
    {
    if(weight(i,0)>0)
      {
      worky.putRowBlock(i*nrcat2,(i+1)*nrcat2,eta.getRowBlock(i*nrcat2,(i+1)*nrcat2)+
                        workweights.getRowBlock(i*nrcat2,(i+1)*nrcat2).inverse()*
                        (respind.getRowBlock(i*nrcat2,(i+1)*nrcat2)-mu.getRowBlock(i*nrcat2,(i+1)*nrcat2)));
      }
    }

  }

void remlest_multinomial::compute_sscp(datamatrix & H, datamatrix & workweight)
  {
  unsigned i;
  H=datamatrix(H.rows(),H.cols(),0);
  datamatrix Htemp=datamatrix(H.rows(),H.cols(),0);
  for(i=0; i<nrobs; i++)
    {
    Htemp=kronecker(workweight.getRowBlock(i*nrcat2,(i+1)*nrcat2),(X.getRow(i).transposed())*(X.getRow(i)));
    H.plus(Htemp);
    }
  }

void remlest_multinomial::compute_sscp2(datamatrix & H, datamatrix & workweight)
  {
  unsigned i;
  H=datamatrix(H.rows(),H.cols(),0);
  datamatrix Htemp=datamatrix(H.rows(),H.cols(),0);
  datamatrix weighttemp=datamatrix(nrcat2,nrcat2,0);
  for(i=0; i<nrobs; i++)
    {
    weighttemp=workweight.getRowBlock(i*nrcat2,(i+1)*nrcat2);
    Htemp.putBlock(kronecker(weighttemp,(X.getRow(i).transposed())*(X.getRow(i))),0,0,totalnrfixed,totalnrfixed);
    Htemp.putBlock(kronecker(weighttemp,(Z.getRow(i).transposed())*(Z.getRow(i))),totalnrfixed,totalnrfixed,totalnrpar,totalnrpar);
    Htemp.putBlock(kronecker(weighttemp,(X.getRow(i).transposed())*(Z.getRow(i))),0,totalnrfixed,totalnrfixed,totalnrpar);
    Htemp.putBlock(Htemp.getBlock(0,totalnrfixed,totalnrfixed,totalnrpar).transposed(),totalnrfixed,0,totalnrpar,totalnrfixed);
    H.plus(Htemp);
    }
  }

void remlest_multinomial::compute_sscp_resp(datamatrix & H1, datamatrix & workweight, datamatrix & worky)
  {
  unsigned i;
  H1=datamatrix(H1.rows(),1,0);
  datamatrix H1temp = datamatrix(nrcat2,1,0);
  datamatrix xtemp;
  for(i=0; i<nrobs;i++)
    {
    H1temp = workweight.getRowBlock(i*nrcat2,(i+1)*nrcat2)*worky.getRowBlock(i*nrcat2,(i+1)*nrcat2);
    xtemp = X.getRow(i).transposed();
    H1.plus(kronecker(H1temp,xtemp));
    }
  }

void remlest_multinomial::compute_sscp_resp2(datamatrix & H1, datamatrix & workweight, datamatrix & worky)
  {
  unsigned i;
  H1=datamatrix(H1.rows(),1,0);
  datamatrix H1temp = datamatrix(H1.rows(),1,0);
  datamatrix weightytemp=datamatrix(nrcat2,1,0);
  for(i=0; i<nrobs;i++)
    {
    weightytemp=workweight.getRowBlock(i*nrcat2,(i+1)*nrcat2)*worky.getRowBlock(i*nrcat2,(i+1)*nrcat2);
    H1temp.putRowBlock(0,totalnrfixed,kronecker(weightytemp,X.getRow(i).transposed()));
    H1temp.putRowBlock(totalnrfixed,totalnrpar,kronecker(weightytemp,Z.getRow(i).transposed()));
    H1.plus(H1temp);
    }
  }

//------------------------------------------------------------------------------
//----------------------------- Object description -----------------------------
//------------------------------------------------------------------------------

  void remlest_multinomial::outoptions()
    {
    out("\n");
    out("GENERAL OPTIONS:\n",true);
    out("\n");
    out("  Maxmimum number of iterations:          "+ST::inttostring(maxit)+"\n");
    out("  Termination criterion:                  "+ST::doubletostring(eps,7)+"\n");
    out("  Stopping criterion for small variances: "+ST::doubletostring(lowerlim,6)+"\n");
    out("\n");
    out("RESPONSE DISTRIBUTION:\n",true);
    out("\n");
    ST::string familyname = "multinomial logit";
    out("  Family:                 "+familyname+"\n");
    out("  Number of observations: "+ST::inttostring(X.rows())+"\n");
    out("  Number of observations with positive weight: "+ST::inttostring(nrobspos)+"\n");
    }

//------------------------------------------------------------------------------
//----------------------------- Writing results --------------------------------
//------------------------------------------------------------------------------

void remlest_multinomial::make_plots(ofstream & outtex,ST::string path_batch,
                         ST::string path_splus)
  {

  char hcharu = '_';
  ST::string hstringu = "\\_";

  unsigned i,j;
  ST::string pathresult;
  bool stil = false;

// Schleife überprüft, ob es ein fullcond-Object
// gibt, bei dem Effekt gezeichnet werden kann
  MCMC::plotstyles plst;
  for(j=0;j<fullcond.size();j++)
    {
    plst = fullcond[j]->get_plotstyle();
    if(plst != MCMC::noplot)
      stil = true;
    }


  if(stil == true)
    {
//erzeugt File, das Plot-Befehle für Java-Version enthält
    ofstream outbatch(path_batch.strtochar());

//erzeugt File, das SPlus-Befehle zum Plotten enthält
    ofstream outsplus(path_splus.strtochar());

    outtex << "\n\\newpage" << "\n\\noindent {\\bf \\large Plots:}" << endl;

    outsplus << "library(\"BayesX\")\n\n";

/*    outsplus << "# NOTE: 'directory' has to be substituted by the directory"
             << " where the functions are stored \n"
             << endl
             << "# In S-PLUS the file extension in the source command has to be changed"
             << " to '.s' \n"
             << endl
    // einlesen der Source-Files für S-Plus
             << "source(\"'directory'\\\\sfunctions\\\\plotsample.r\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\plotnonp.r\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\plotsurf.r\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\drawmap.r\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\readbndfile.r\")\n" << endl;*/

#if defined(JAVA_OUTPUT_WINDOW)
    out("  --------------------------------------------------------------------------- \n");
    out("\n");
    out("  Batch file for visualizing effects of nonlinear functions is stored in file \n");
    out("  " + path_batch + "\n");
    out("\n");
#endif

    bool stil2 = true;
    for(j=0;j<fullcond.size();j++)  //Schleife überprüft, ob es map-Objekt gibt
      {
      plst = fullcond[j]->get_plotstyle();
      if(plst == MCMC::drawmap || plst == MCMC::drawmapgraph)
        stil2 = false;
      }

    if(stil2 == true)
      {
      out("  --------------------------------------------------------------------------- \n");
      out("\n");
      out("  Batch file for visualizing effects of nonlinear functions ");
      out("  in R is stored in file \n");
      out("  " + path_splus + "\n");
      out("\n");
      }

    if(stil2 == false)
      {
      out("  NOTE: 'input filename' must be substituted by the filename of the boundary-file \n");
      out("\n");
      out("  --------------------------------------------------------------------------- \n");
      out("\n");
      out("  Batch file for visualizing effects of nonlinear functions ");
      out("  in R is stored in file \n");
      out("  " + path_splus + "\n");
      out("\n");
      out("  NOTE: 'input filename' must be substituted by the filename of the boundary-file \n");
      out("\n");
      }


    outbatch << "% usefile " << path_batch << endl;

    // falls andere Quantile gewünscht werden
    double u = fullcond[0]->get_level1();
    double o = fullcond[0]->get_level2();
    ST::string u_str = ST::doubletostring(u,0);
    ST::string o_str = ST::doubletostring(o,0);

    for(i=0; i<nrcat2; i++)
      {
    // durchlaufen der Fullconditionals
    for(j=0;j<fullcond.size();j++)
      {

      // Pfad der Regr.-Ergebnisse
      pathresult = fullcond[j]->get_pathresult();
      pathresult = pathresult.insert_after_string(ST::doubletostring(cats(i,0),6)+"_","_f_");

      // Plotstyle: noplot, plotnonp, drawmap, drawmapgraph
      plst = fullcond[j]->get_plotstyle();

      if (plst != MCMC::noplot)
        {

        // Pfade für ps-, tex-, SPlus-files
        ST::string pathps = pathresult.substr(0, pathresult.length()-4);
        ST::string pathgr = pathps.replaceallsigns('\\', '/');

        char hchar = '\\';
        ST::string hstring = "/";

        ST::string pathps_spl = pathps.insert_string_char(hchar,hstring);
        ST::string pathres_spl = pathresult.insert_string_char(hchar,hstring);

        if (plst == MCMC::plotnonp)
          {
          outbatch << "\n";                // Befehle f. d. batch-file
          outbatch << "dataset _dat" << endl;
          outbatch << "_dat.infile using " << pathresult << endl;
          outbatch << "graph _g" << endl;
          vector<ST::string> varnames = fullcond[j]->get_datanames();
          ST::string xvar = varnames[0];
          outbatch << "_g.plot " << xvar
                   << " pmode ci" << u_str << "lower ci"
                   << o_str.replaceallsigns('.','p') << "lower ci"
                   << o_str.replaceallsigns('.','p') << "upper ci"
                   << u_str.replaceallsigns('.','p') << "upper, "
                   << "title = \"Effect of " << xvar << "\" xlab = " << xvar
                   << " ylab = \" \" " << "outfile = " << pathps
                   << ".ps replace using _dat" << endl;
          outbatch << "drop _dat" << endl;
          outbatch << "drop _g" << endl;
          // Plot-Befehle f. d. SPlus-file
          outsplus << "plotnonp(\"" << pathres_spl << "\")" << endl;
          // Plot-Befehle f. d. tex-file
          ST::string effect = xvar;
          if(varnames.size()>1)
            {
            effect = varnames[1] + "*" + effect;
            }
          outtex << "\n\\begin{figure}[h!]" << endl
                  << "\\centering" << endl
                  << "\\includegraphics[scale=0.6]{" << pathgr << "}" << endl
                  << "\\caption{Non--linear Effect of '" <<
                  effect.insert_string_char(hcharu,hstringu) << "'";
          outtex << " (Category " << cats(i,0) << ")." << endl << "Shown are the posterior modes together with "
                 << u_str << "\\% and " << o_str
                 << "\\% pointwise credible intervals.}" << endl
                 << "\\end{figure}" << endl;
          }
        // für map-Funktionen
        else if (plst == MCMC::drawmap || plst == MCMC::drawmapgraph)
          {
          outbatch << "\n";                 // Befehle f. d. batch-file
          outbatch << "dataset _dat" << endl;
          outbatch << "_dat.infile using " << pathresult << endl;
          outbatch << "map _map" << endl;
          outbatch << "_map.infile using input_filename" << endl;
          outbatch << "graph _g" << endl;
          vector<ST::string> varnames = fullcond[j]->get_datanames();
          ST::string regionvar = varnames[0];
          outbatch << "_g.drawmap " << "pmode" << " " << regionvar
                   << ", map = _map color outfile = " << pathps
                   << "_pmode.ps replace using _dat" << endl;
          outbatch << "_g.drawmap " << "pcat" << u_str << " " << regionvar
                     << ", map = _map nolegend pcat outfile = " << pathps
                     << "_pcat" << u_str << ".ps replace using _dat" << endl;
          outbatch << "_g.drawmap " << "pcat" << o_str << " " << regionvar
                     << ", map = _map nolegend pcat outfile = " << pathps
                     << "_pcat" << o_str << ".ps replace using _dat" << endl;
          outbatch << "drop _dat" << endl;
          outbatch << "drop _g" << endl;
          outbatch << "drop _map" << endl;
          // Plot-Befehle f. d. SPlus-file
          outsplus << "# NOTE: 'input_filename' must be substituted by the "
                   << "filename of the boundary-file \n"
//                   << "# NOTE: choose a 'name' for the map \n" << endl
                     << "m <- read.bnd(\"'input_filename'\")" << endl
                   << "drawmap(data = \"" << pathres_spl << "\", map = m, plotvar = \"pmode\", regionvar = \""
                   << regionvar << "\")" << endl;
          outsplus << "drawmap(data = \"" << pathres_spl << "\", map = m, plotvar = \"pcat" << u_str << "\", regionvar = \""
                  << regionvar << "\", legend = F, pcat = T)" << endl;
          outsplus << "drawmap(data = \"" << pathres_spl << "\", map = m, plotvar = \"pcat" << o_str << "\", regionvar = \""
                    << regionvar << "\", legend = F, pcat = T)" << endl;
            // Plot-Befehle f. d. tex-file
          ST::string effect = regionvar;
          if(varnames.size()>1)
            {
            effect = varnames[1] + "*" + effect;
            }

          if(plst == MCMC::drawmap)
            {
            outtex << "\n\\begin{figure}[h!]" << endl
                   << "\\centering" << endl
                   << "\\includegraphics[scale=0.6]{" << pathgr << "_pmode}"
                   << endl
                   << "\\caption{Non--linear Effect of '" <<
                   effect.insert_string_char(hcharu,hstringu) << "'";
            outtex << " (Category " << cats(i,0) << "). Shown are the posterior modes.}" << endl
                   << "\\end{figure}" << endl;
            outtex << "\n\\begin{figure}[htb]" << endl
                   << "\\centering" << endl
                   << "\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                   << u_str << "}" << endl
                   << "\\caption{Non--linear Effect of '" << effect << "'";
            outtex << " (Category " << cats(i,0) << "). Posterior probabilities for a nominal level of "
                   << u_str << "\\%." << endl
                   << "Black denotes regions with strictly negative credible intervals,"
                   << endl
                   << "white denotes regions with strictly positive credible intervals.}"
                   << endl << "\\end{figure}" << endl;
            outtex << "\n\\begin{figure}[htb]" << endl
                   << "\\centering" << endl
                   << "\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                   << o_str << "}" << endl
                   << "\\caption{Non--linear Effect of '" << effect << "'";
            outtex << " (Category " << cats(i,0) << "). Posterior probabilities for a nominal level of "
                   << o_str << "\\%." << endl
                   << "Black denotes regions with strictly negative credible intervals,"
                   << endl
                   << "white denotes regions with strictly positive credible intervals.}"
                   << endl << "\\end{figure}" << endl;
            }
          else if(plst == MCMC::drawmapgraph)
            {
            outtex << "\n%\\begin{figure}[h!]" << endl
                   << "%\\centering" << endl
                   << "%\\includegraphics[scale=0.6]{" << pathgr << "_pmode}"
                   << endl
                   << "%\\caption{Non--linear Effect of '" <<
                   effect.insert_string_char(hcharu,hstringu) << "'";
            outtex << " (Category " << cats(i,0) << "). Shown are the posterior modes.}" << endl
                   << "%\\end{figure}" << endl;
            outtex << "\n%\\begin{figure}[htb]" << endl
                   << "%\\centering" << endl
                   << "%\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                   << u_str << "}" << endl
                   << "%\\caption{Non--linear Effect of '" << effect << "'";
            outtex << " (Category " << cats(i,0) << "). Posterior probabilities for a nominal level of "
                   << u_str << "\\%." << endl
                   << "%Black denotes regions with strictly negative credible intervals,"
                   << endl
                   << "%white denotes regions with strictly positive credible intervals.}"
                   << endl << "%\\end{figure}" << endl;
            outtex << "\n%\\begin{figure}[htb]" << endl
                   << "%\\centering" << endl
                   << "%\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                   << o_str << "}" << endl
                   << "%\\caption{Non--linear Effect of '" << effect << "'";
            outtex << " (Category " << cats(i,0) << "). Posterior probabilities for a nominal level of "
                   << o_str << "\\%." << endl
                   << "%Black denotes regions with strictly negative credible intervals,"
                   << endl
                   << "%white denotes regions with strictly positive credible intervals.}"
                   << endl << "%\\end{figure}" << endl;
            }

          } // end: else if
        } // end: if
      } // end: for
      }
    }
  }

void remlest_multinomial::make_model(ofstream & outtex, const ST::string & rname)
  {
  ST::string familyname;
  if(respfamily=="multinomial")
    {
    familyname="multinomial logit";
    }

  //Anz. Beob. wird übergeben
  unsigned obs = X.rows();

  char charh = '_';
  ST::string stringh = "\\_";
  ST::string helprname = rname.insert_string_char(charh,stringh);

  //schreibt das Modell und die Priori-annahmen ins Tex-File
  outtex << "\n\\noindent {\\bf \\large Response:}" << endl
         << "\\begin{tabbing}\n"
         << "Number of observations with positive weight: \\= \\kill" << endl
         << "Number of observations: \\> " << obs << "\\\\" << endl
         << "Number of observations with positive weight: \\> " << nrobspos << "\\\\" << endl
         << "Response Variable: \\> " << helprname << "\\\\" << endl
         << "Family: \\> " << familyname << "\\\\" << endl
         << "\\end{tabbing}" << endl
         << "\n\\noindent {\\bf \\large Predictor:}\\\\" << endl;
  }

void remlest_multinomial::make_predictor(ofstream & outtex)
  {

  unsigned j;

  ST::string term2 = fullcond[0]->get_term_symbolic();
  ST::string term = "$\\eta$ & $=$ & $" + term2;    //linearer Prädiktor wird erweitert
  for(j=1;j<fullcond.size();j++)
    {
    out(fullcond[j]->get_results_type());
    term2 = fullcond[j]->get_term_symbolic();
    term = term + " + " + term2;    //linearer Prädiktor wird erweitert
    }

  outtex << endl << "\n\\begin{tabular}{ccp{12cm}}\n";
  for(j=0; j<nrcat2; j++)
    {
    term2 = term.insert_after_string("^{("+ST::doubletostring(cats(j,0),6)+")}","\\eta");
    term2 = term2.insert_after_all_string("^{("+ST::doubletostring(cats(j,0),6)+")}","\\gamma");
    term2 = term2.insert_after_all_string("^{("+ST::doubletostring(cats(j,0),6)+")}","+ f");

    outtex << term2 << "$\\\\\n";
    }
  outtex << "\\end{tabular}\n\\\\ \n\\\\" << endl;
  }

void remlest_multinomial::make_prior(ofstream & outtex)
  {
  unsigned i,j;
  outtex << "\n\\noindent {\\bf \\large Priors:}\\\\" << endl << "\\\\" << endl;
  for(j=0;j<fullcond.size();j++)
    {
    vector<ST::string> prior = fullcond[j]->get_priorassumptions();
    if(prior.size() != 0)// nur wenn Priors da sind (d.h. Vektor hat Elemente)
      {
      if(fullcond[j]->get_results_type()!="fixed")
        {
        prior[0] = prior[0].insert_after_string("^{(j)}","f");
        }
      for(i=0;i<prior.size();i++)
        {
        if( j!=0 || i<prior.size()-1)
          {
          outtex << prior[i] << "\\\\" << endl;
          }
        }
      outtex << "\\\\" <<endl;
      }
    }
  }

void remlest_multinomial::make_options(ofstream & outtex)
  {
  double l1 = fullcond[0]->get_level1();
  double l2 = fullcond[0]->get_level2();

  //schreibt REML options ins Tex-File
  outtex << "\n\\noindent {\\bf \\large General Options:}" << endl
         << "\\begin{tabbing}" << endl
         << "Levels for credible intervals: \\hspace{2cm}\\= \\\\" << endl
         << "Level 1: \\> " << l1 << "\\\\" << endl
         << "Level 2: \\> " << l2 << "\\\\" << endl
         << "Maxmimum number of iterations: \\> " << maxit << "\\\\" << endl
         << "Termination criterion: \\> " << eps << "\\\\" << endl
         << "Stopping criterion for small variances: \\> " << lowerlim << endl
         << "\\end{tabbing}\n"  << "\\vspace{0.5cm}" <<  endl;
  }

void remlest_multinomial::make_fixed_table(ofstream & outtex)
  {

  // falls andere Quantile gewünscht werden
  double u = fullcond[0]->get_level1();
  ST::string u_str = ST::doubletostring(u,0);

  vector<ST::string> h;
  h = fullcond[0]->get_results_latex();

  unsigned i,j, r;
  unsigned partialhsize=h.size()/nrcat2;

  for(i=0; i<nrcat2; i++)
    {
    r=2;
    // Tabelle im Tex-File mit fixen Effekten
    outtex << "\n\\newpage \n" << endl << "\n\\noindent {\\bf \\large Fixed Effects (Category "
         << ST::doubletostring(cats(i,0),6) << "):}\\\\"
         << endl << "\\\\" << endl;

    outtex << "\\begin{tabular}{|r|rrrrr|}" << endl << "\\hline" << endl
           << "Variable & Post. Mode & Std. Dev. & p-value & \\multicolumn{2}{r|}{" << u << "\\% confidence interval}\\\\"
           << endl << "\\hline" << endl;

    for (j=i*partialhsize;j<(i+1)*partialhsize;j++)
      {
      r++;
      if (r < 39)
        {
        outtex << h[j] << endl;
        }
      else
        {
        r=1;
        outtex << "\\hline \n\\end{tabular}" << endl;

        outtex << "\n\\newpage \n" << endl
               << "\n\\noindent {\\bf \\large Fixed Effects (Category "
               << ST::doubletostring(cats(i,0),6) << "continued):}\\\\"
               << endl << "\\\\" << endl;

        outtex << "\\begin{tabular}{|r|rrrrr|}" << endl << "\\hline" << endl
               << "Variable & Post. Mode & Std. Dev. & p-value & \\multicolumn{2}{r|}{" << u << "\\% confidence interval}\\\\"
               << endl << "\\hline" << endl;

        outtex << h[j] << endl;
        }
      }
    outtex << "\\hline \n\\end{tabular}" << endl;
    }
  }

void remlest_multinomial::make_graphics(const ST::string & title,
                     const ST::string & path_batch,
                     const ST::string & path_tex,
                     const ST::string & path_splus,
                     const ST::string & rname)
  {
  ST::string pathresult;                 //Pfad des Ergebnis-Files

  vector<ST::string> distr_results;

 // erzeugt Tex-File
  ofstream outtex(path_tex.strtochar());

  char charh = '_';
  ST::string stringh = "\\_";
  ST::string helptitle = title.insert_string_char(charh,stringh);

  //erzeugt den Kopf des Tex-Files
  outtex << "\\documentclass[a4paper, 12pt]{article}" << endl
         << "\n" << "\\usepackage{graphicx}" << endl
         << "\\parindent0em" << endl
         << "\n\\begin{document}" << endl
         << "\\begin{center}" << endl
         << "\\LARGE{\\bf " << helptitle << "}"
         << endl << "\\end{center} \n\\vspace{1cm}" << endl;

  make_model(outtex,rname);

  make_predictor(outtex);

  make_prior(outtex);

  make_options(outtex);

  outtex << "\n\\noindent {\\bf \\large Model Fit:}" << endl
         << "\\begin{tabbing}\n";
  outtex << "GCV (based on deviance residuals): \\= \\kill" << endl;
  outtex << "-2*log-likelihood: \\> " << loglike << "\\\\" << endl;
  outtex << "Degrees of freedom: \\> " << df << "\\\\" << endl;
  outtex << "(conditional) AIC: \\> " << aic << "\\\\" << endl;
  outtex << "(conditional) BIC: \\> " << bic << "\\\\" << endl;
  outtex << "GCV (based on deviance residuals): \\> " << gcv << "\\\\" << endl;
  outtex << "\\end{tabbing}" << endl;

  make_fixed_table(outtex);

  // Pfade der Files
  //werden im BayesX-Output angegeben
  out("  Files of model summary: \n" , true);
  out("\n");

  make_plots(outtex,path_batch,path_splus);

  out("  --------------------------------------------------------------------------- \n");
  out("\n");
  out("  Latex file of model summaries is stored in file \n");
  out("  " + path_tex + "\n");
  out("\n");
  out("  --------------------------------------------------------------------------- \n");
  out("\n");


  outtex << "\\end{document}" << endl;

  }

bool remlest_multinomial::check_pause()
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  Application->ProcessMessages();
  if (Frame->stop)
    {
    return true;
    }

  if (Frame->pause)
    {
    out("\n");
    out("ESTIMATION PAUSED\n");
    out("Click CONTINUE to proceed\n");
    out("\n");

    while (Frame->pause)
      {
      Application->ProcessMessages();
      }

    out("ESTIMATION CONTINUED\n");
    out("\n");
    }
  return false;
#elif defined(JAVA_OUTPUT_WINDOW)
  return adminb_p->breakcommand();
#else
  return false;
#endif
  }

void remlest_multinomial::out(const ST::string & s,bool thick,bool italic,
                      unsigned size,int r,int g, int b)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  ST::string sh = s;
  sh = sh.replaceallsigns('\n',' ');
  if (!Frame->suppoutput)
    Results->ResultsRichEdit->Lines->Append(sh.strtochar());
 if (!(logout->fail()))
    (*logout) << s << flush;
#elif defined(JAVA_OUTPUT_WINDOW)
  ST::string sh = s;
  sh = sh.replaceallsigns('\n',' ');
  sh = sh+"\n";
  if (!adminb_p->get_suppressoutput())
    adminb_p->Java->CallVoidMethod(adminb_p->BayesX_obj, adminb_p->javaoutput,
    adminb_p->Java->NewStringUTF(sh.strtochar()),thick,italic,size,r,g,b);
  if (!(logout->fail()))
    (*logout) << s << flush;
#else
  cout << s << flush;
  if (!(logout->fail()))
    (*logout) << s << flush;
#endif
  }


void remlest_multinomial::outerror(const ST::string & s)
  {
  out(s,true,true,12,255,0,0);
  }


//------------------------------------------------------------------------------
//----------------------------- Constructor ------------------------------------
//------------------------------------------------------------------------------

remlest_multistate::remlest_multistate(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adb,
#endif
vector<MCMC::FULLCOND*> & fc,datamatrix & re,
                const ST::string & family, const ST::string & ofile,
                const int & maxiter, const double & lowerlimit,
                const double & epsi, const double & maxch, const double & maxv,
                const bool & cl, const bool & glfrailty,
                const vector<unsigned> & nrfullc,
                const datamatrix & weight, ostream * lo)
  {

    #if defined(JAVA_OUTPUT_WINDOW)
    adminb_p = adb;
    #endif

    logout = lo;
    respfamily=family;
    outfile=ofile;

    maxit=maxiter;
    lowerlim=lowerlimit;
    eps=epsi;
    maxchange=maxch;
    maxvar=maxv;

    constlambda = cl;

    globalfrailty=glfrailty;

    nrtransitions = re.cols();
    nrfullconds = nrfullc;

    fullcond = fc;
    unsigned i, j, k, l;

    xcut.push_back(0);
    zcut.push_back(0);
    xcuttrans.push_back(0);
    zcuttrans.push_back(0);

    k=0;
    for(i=0; i<nrfullconds.size(); i++)
      {
      for(j=0; j<nrfullconds[i]; j++)
        {
        xcut.push_back(xcut[k]+fullcond[k]->get_dimX());
        k++;
        }
      xcuttrans.push_back(xcut[k]);
      }
    if(globalfrailty)
      {
      xcut.push_back(xcut[k]+fullcond[k]->get_dimX());
      }
    k=l=0;
    for(i=0; i<nrfullconds.size(); i++)
      {
      k++;
      for(j=1; j<nrfullconds[i]; j++)
        {
        zcut.push_back(zcut[l]+fullcond[k]->get_dimZ());
        k++;
        l++;
        }
      zcuttrans.push_back(zcut[l]);
      }
    if(globalfrailty)
      {
      zcut.push_back(zcut[l]+fullcond[k]->get_dimZ());
      }

    X = datamatrix(re.rows(),xcut[xcut.size()-1],0);
    Z = datamatrix(re.rows(),zcut[zcut.size()-1],0);

    k=l=0;
    for(i=0; i<nrfullconds.size(); i++)
      {
      fullcond[k]->createreml(X,Z,xcut[k],0);
      k++;
      for(j=1; j<nrfullconds[i]; j++)
        {
        fullcond[k]->createreml(X,Z,xcut[k],zcut[l]);
        k++;
        l++;
        }
      }
    if(globalfrailty)
      {
      fullcond[k]->createreml(X,Z,xcut[k],zcut[l]);
      }

    beta=statmatrix<double>(X.cols()+Z.cols(),1,0);
    theta=statmatrix<double>(zcut.size()-1,1,0);

    k=l=0;
    for(i=0; i<nrfullconds.size(); i++)
      {
      k++;
      for(j=1; j<nrfullconds[i]; j++)
        {
        theta(l,0) = fullcond[k]->get_startlambda();
        k++;
        l++;
        }
      }
    if(globalfrailty)
      {
      theta(l,0) = fullcond[k]->get_startlambda();
      }
    }

//------------------------------------------------------------------------------
//----------------------------- REML estimation --------------------------------
//------------------------------------------------------------------------------

  // Function: estimate
  // Task: Perform REML-estimation for multi state models
  //       returns true if an error or user break occured

bool remlest_multistate::estimate(const datamatrix resp,
               const datamatrix & offset, const datamatrix & weight,
               const datamatrix & state)
  {
  unsigned i, j, k, l, m, n;
  double help, former;

  outoptions();
  out("\n");

  for(i=0;i<fullcond.size();i++)
    fullcond[i]->outoptionsreml();

  out("\n");
  out("REML ESTIMATION STARTED\n",true);
  out("\n");

  bool stop = check_pause();
  if (stop)
    return true;

  // Matrix to store old versions of beta and theta
  statmatrix<double>betaold(beta.rows(),1,0);
  statmatrix<double>thetaold(theta.rows(),1,0);

  // Score-function and expected Fisher information for theta
  statmatrix<double>score(theta.rows(),1,0);
  statmatrix<double>Fisher(theta.rows(),theta.rows(),0);

  // Matrices for Fisher scoring (regression parameters)
  statmatrix<double>H(beta.rows(),beta.rows(),0);
  statmatrix<double>Hinv(beta.rows(),beta.rows(),0);
  statmatrix<double>H1(beta.rows(),1,0);

  // Number of iterations, nr of observations
  unsigned it=1;
  unsigned nrobs = Z.rows();
  unsigned xcols = X.cols();
  unsigned zcols = Z.cols();

  // Criteria to detemine convergence
  double crit1=1;                //relative changes in regression parameters
  double crit2=1;                //relative changes in variance parameters
  bool test=true;

  vector<double>stopcrit(theta.rows(),10);
  vector<int>its(theta.rows(),0);
  vector<int>thetastop(theta.rows(),0);
  vector<int>signs(theta.rows(),1);

  // Matrix containing the inverse covariance matrix of the random effects
  statmatrix<double>Qinv(zcols,1,0);

  // Inzidenzvektor, die für jeden Wert in fullcond bzw. beta angibt, ob er zur Baseline-HR beiträgt
  vector<int>isbaseline(fullcond.size(),0);
  vector<bool> timevarying(nrtransitions,false);
  unsigned nrbaseline=0;
  vector<int> basecut;
  basecut.push_back(0);

  k=0;
  for(i=0; i<nrtransitions; i++)
    {
    l=0;
    for(j=0; j<nrfullconds[i]; j++)
      {
      if(fullcond[k]->is_baseline()==true)
        {
        isbaseline[k]=1;
        nrbaseline++;
        l++;
        }
      k++;
      }
    basecut.push_back(basecut[i]+l);
    if(l>1)
      {
      timevarying[i] = true;
      }
    }

  vector<int>cumtimevar(nrtransitions,0);
  if(nrbaseline>nrtransitions)
    {
    k=0;
    for(j=1; j<nrtransitions; j++)
      {
      cumtimevar[j] = k;
      if(timevarying[j])
        {
        k++;
        }
      }
    }

  vector<int>isbaselinebeta(beta.rows(),0);
  vector<int>dm_pos(beta.rows(),0);

  for(i=0; i<fullcond.size(); i++)
    {
    if(isbaseline[i]==1)
      {
      k=1;
      for(j=xcut[i]; j<xcut[i+1]; j++, k++)
        {
        isbaselinebeta[j]=1;
        dm_pos[j] = k;
        }
      }
    }

  k=0;
  for(i=0; i<nrtransitions; i++)
    {
    l=0;
    for(j=0; j<nrfullconds[i]; j++)
      {
      if(isbaseline[k])
        {
        if(l==0)
          {
          // baseline
          dm_pos[xcut[k]] = 1;
          l++;
          }
        else
          {
          // time-varying effect
          dm_pos[xcut[k]] = 0;
          dm_pos[xcut[k]+1] = 1;
          }
        }
      k++;
      }
    }

  k=l=0;
  for(i=0; i<nrtransitions; i++)
    {
    k++;
    for(j=1; j<nrfullconds[i]; j++)
      {
      if(isbaseline[k]==1)
        {
        m=0;
        for(n=zcut[l]; n<zcut[l+1]; n++)
          {
          isbaselinebeta[xcols+n]=1;
          dm_pos[xcols+n]=m;
          m++;
          }
        }
      k++;
      l++;
      }
    }

// Inzidenzvector defining the assignment between beta and position in the interactvar-matrix

  vector<int>fc_pos;
  if(nrbaseline>nrtransitions)
    {
    j=0;
    fc_pos = vector<int>(beta.rows(),0);
    for(i=0; i<fullcond.size(); i++)
      {
      if(isbaseline[i])
        {
        for(k=xcut[i]; k<xcut[i+1]; k++)
          {
          fc_pos[k] = j;
          }
        j++;
        }
      }
    l=m=n=0;;
    for(i=0; i<nrtransitions; i++)
      {
      l++;
      for(j=1; j<nrfullconds[i]; j++)
        {
        if(isbaseline[l])
          {
          for(k=zcut[m]; k<zcut[m+1]; k++)
            {
            fc_pos[xcols+k] = n;
            }
          n++;
          }
        m++;
        l++;
        }
      }
    }

/*datamatrix help1(fc_pos.size(),1,0);
for(i=0; i<fc_pos.size(); i++)
  {
  help1(i,0) = fc_pos[i];
  }
ofstream out1("c:\\temp\\fc_pos.raw");
help1.prettyPrint(out1);
out1.close();*/

// Inzidenzvector defining the assignment between beta and transition

  vector<int>tr_pos(beta.rows(),1);
  k=0;
  for(i=0; i<nrtransitions; i++)
    {
    for(j=xcuttrans[i]; j<xcuttrans[i+1]; j++)
      {
      tr_pos[k] = i;
      k++;
      }
    }
  for(i=0; i<nrtransitions; i++)
    {
    for(j=zcuttrans[i]; j<zcuttrans[i+1]; j++)
      {
      tr_pos[k] = i;
      k++;
      }
    }
  // -1 denotes global effects
  if(globalfrailty)
    {
    for(j=zcut[zcut.size()-2]; j<zcut[zcut.size()-1]; j++)
      {
      tr_pos[k] = -1;
      k++;
      }
    }

// Matrices and variables for baseline effects
  datamatrix tsteps;
  datamatrix t_X;
  datamatrix t_Z;
  vector<unsigned> tstart;
  vector<unsigned> tend;
  vector<unsigned> ttrunc;
  datamatrix interactvar(nrobs,nrbaseline,0);
  statmatrix<int> index(nrobs,1,0);
  j=0;
  for(i=0; i<fullcond.size(); i++)
    {
    if(isbaseline[i])
      {
      fullcond[i]->initialize_baseline(j,t_X,t_Z,tstart,tend,ttrunc,interactvar,tsteps,index);
      j++;
      }
    }

  // Baseline-HR, linear predictor
  datamatrix basef(t_X.rows(),nrbaseline,0);
  statmatrix<double>baseline(t_X.rows(),nrtransitions,0);
  statmatrix<double>cumbaseline(t_X.rows(),nrtransitions,0);
  statmatrix<double>baseline_tv;
  statmatrix<double>cumbaseline_tv;

  if(nrbaseline>nrtransitions)
    {
    baseline_tv=datamatrix(nrobs,t_X.rows()*(cumtimevar[nrtransitions-1]+1),0);
    cumbaseline_tv=datamatrix(nrobs,nrtransitions,0);
    }
  statmatrix<double>eta(nrobs,nrtransitions,0);
  statmatrix<double>baseline_eta(nrobs,nrtransitions,0);
  statmatrix<double>mult_eta(nrobs,nrtransitions,0);
  statmatrix<double>mult_hazard(nrobs,nrtransitions,0);
  statmatrix<double>mult_global(nrobs,1,1);
  datamatrix helpmat;

  // incidence matrix of possible transitions
  helpmat=state;
  helpmat.sort(0,nrobs-1,0);
  vector<int> states;
  states.push_back(helpmat(0,0));
  for(i=1; i<nrobs; i++)
    {
    if(helpmat(i-1,0)!=helpmat(i,0))
      {
      states.push_back(helpmat(i,0));
      }
    }
  statmatrix<bool>transition(states.size(), nrtransitions, false);
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrtransitions; j++)
      {
      for(k=0; k<states.size(); k++)
        {
        if(state(i,0)==states[k] && resp(i,j)==1)
          {
          transition(k,j)=true;
          }
        }
      }
    }

  // matrix defining the transitions an observation is under risk for
  datamatrix riskset(nrobs, nrtransitions, 0);
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrtransitions; j++)
      {
      for(k=0; k<states.size(); k++)
        {
        if(state(i,0)==states[k] && transition(k,j))
          {
          riskset(i,j)=true;
          }
        }
      }
    }

// check for incorrect times

  for(i=0; i<nrobs; i++)
    {
    if(resp(i,0)==1 && tstart[i] < tend[i])
      {
      outerror("ERROR: observation "+ST::inttostring(i+1)+" has censoring indicator 1\n");
      outerror("       and left interval time smaller than right interval time\n");
      return true;
      }
    if(ttrunc[i] > tend[i])
      {
      outerror("ERROR: left-truncation time larger then observed survival time\n");
      outerror("       for observation "+ST::inttostring(i+1)+"\n");
      return true;
      }
    if(ttrunc[i] > tstart[i])
      {
      outerror("ERROR: left-truncation time larger then left interval time\n");
      outerror("       for observation "+ST::inttostring(i+1)+"\n");
      return true;
      }
   if(tstart[i]>tend[i])
       {
      outerror("ERROR: left interval time larger then right interval time");
      outerror("       for observation "+ST::inttostring(i+1)+"\n");
      return true;
      }
    }

  // Transform smoothing paramater starting values to variances
  for(i=0; i<theta.rows(); i++)
    {
    theta(i,0)=1/theta(i,0);
    }

  // Startwerte für beta0
  j=k=0;
  for(i=0; i<fullcond.size(); i++)
    {
    if(fullcond[i]->get_dimZ()==0)
      {
      beta(j,0) = log(100/t_X(t_X.rows()-1,1));
      k++;
      }
    j += fullcond[i]->get_dimX();
    }

  while(test==true)
    {

    // store current values in betaold and thetaold
    betaold=beta;
    thetaold=theta;

    // compute Qinv
    for(i=0, j=0; i<theta.rows(); i++)
      {
      for(k=zcut[i]; k<zcut[i+1]; k++, j++)
        {
        Qinv(j,0)=1/theta(i,0);
        }
      }

/*ofstream out8("c:\\temp\\Qinv.raw");
Qinv.prettyPrint(out8);
out8.close();*/

    // compute basef
    j=k=m=0;
    for(i=0; i<nrtransitions; i++)
      {
      k++;
      for(l=1; l<nrfullconds[i]; l++)
        {
        if(isbaseline[k]==1)
          {
          if(xcut[k+1]==xcut[k]+1)
            {
            basef.putCol(j,beta(xcut[k],0)*t_X.getCol(1)+t_Z*beta.getRowBlock(xcols+zcut[m],xcols+zcut[m+1]));
            }
          else
            {
            basef.putCol(j,t_X*beta.getRowBlock(xcut[k],xcut[k+1])+t_Z*beta.getRowBlock(xcols+zcut[m],xcols+zcut[m+1]));
            }
          j++;
          }
        k++;
        m++;
        }
      }

    // compute baseline & cumulative baseline

    for(j=0; j<nrtransitions; j++)
      {
      if(!timevarying[j])
        {
        for(i=0; i<baseline.rows(); i++)
          {
          baseline(i,j) = exp(basef(i,basecut[j]));
          }
        }
      else
        {
        helpmat = interactvar.getColBlock(basecut[j],basecut[j+1])*(basef.getColBlock(basecut[j],basecut[j+1])).transposed();
        unsigned off=cumtimevar[j]*t_X.rows();
        for(i=0; i<nrobs; i++)
          {
          for(k=ttrunc[i]; k<tend[i]; k++)
            {
            baseline_tv(i,off+k) = exp(helpmat(i,k));
            }
          }
        }
      }
    cumbaseline = datamatrix(t_X.rows(),nrtransitions,0);
    if(nrbaseline>nrtransitions)
      {
      cumbaseline_tv = datamatrix(nrobs,nrtransitions,0);
      }
    for(j=0; j<nrtransitions; j++)
      {
      if(!timevarying[j])
        {
        former=0;
        for(i=0; i<t_X.rows()-1; i++)
          {
          cumbaseline(i,j) = former + 0.5*tsteps(i,0)*(baseline(i,j)+baseline(i+1,j));
          former = cumbaseline(i,j);
          }
        cumbaseline(t_X.rows()-1,j) = former;
        }
      else
        {
        unsigned off=cumtimevar[j]*t_X.rows();
        for(i=0; i<nrobs; i++)
          {
          for(k=ttrunc[i]; k<tend[i]-1; k++)
            {
            cumbaseline_tv(i,j) += 0.5*tsteps(k,0)*(baseline_tv(i,off+k)+baseline_tv(i,off+k+1));
            }
          }
        }
      }

    // compute eta, mult. hazard, ...

    k=m=0;
    baseline_eta = datamatrix(nrobs,nrtransitions,0);
    for(i=0; i<nrtransitions; i++)
      {
      k++;
      for(l=1; l<nrfullconds[i]; l++)
        {
        if(isbaseline[k]==1)
          {
          baseline_eta.putCol(i, baseline_eta.getCol(i) + X.getColBlock(xcut[k],xcut[k+1])*beta.getRowBlock(xcut[k],xcut[k+1])+Z.getColBlock(zcut[m],zcut[m+1])*beta.getRowBlock(xcols+zcut[m],xcols+zcut[m+1]));
//          j++;
          }
        k++;
        m++;
        }
      }
    for(i=0; i<nrtransitions; i++)
      {
      eta.putCol(i,X.getColBlock(xcuttrans[i],xcuttrans[i+1])*beta.getRowBlock(xcuttrans[i],xcuttrans[i+1]) + Z.getColBlock(zcuttrans[i],zcuttrans[i+1])*beta.getRowBlock(xcols+zcuttrans[i],xcols+zcuttrans[i+1]));
      }
    mult_eta = eta - baseline_eta;
    if(globalfrailty)
      {
      mult_global = Z.getColBlock(zcut[zcut.size()-2],zcut[zcut.size()-1])*beta.getRowBlock(xcols+zcut[zcut.size()-2],xcols+zcut[zcut.size()-1]);
      for(i=0; i<nrobs; i++)
        {
        mult_global(i,0) = exp(mult_global(i,0));
        }
      }

    for(j=0; j<nrtransitions; j++)
      {
      for(i=0; i<nrobs; i++)
        {
        mult_hazard(i,j)=exp(mult_eta(i,j))*mult_global(i,0);
        }
      }

/*ofstream out5("c:\\temp\\baseline.raw");
baseline.prettyPrint(out5);
out5.close();
ofstream out6("c:\\temp\\mult_hazard.raw");
mult_hazard.prettyPrint(out6);
out6.close();
ofstream out7("c:\\temp\\eta.raw");
eta.prettyPrint(out7);
out7.close();*/

  // Score-function for beta

    //clear H1-matrix
    H1 = datamatrix(H1.rows(),1,0);

    // X

    for(j=0; j<X.cols(); j++)
      {
      if(isbaselinebeta[j]==1)
        {
        if(!timevarying[tr_pos[j]])
          {
          // not time-varying
          helpmat = datamatrix(t_X.rows(),1,0);
          former=0;
          for(i=0; i<t_X.rows()-1; i++)
            {
            helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_X(i,dm_pos[j])*baseline(i,tr_pos[j])+t_X(i+1,dm_pos[j])*baseline(i+1,tr_pos[j]));
            former=helpmat(i,0);
            }
          helpmat(t_X.rows()-1,0) = former;

          for(i=0; i<nrobs; i++)
            {
            H1(j,0) += resp(i,tr_pos[j])*X(i,j);
            }
          for(i=0; i<nrobs; i++)
            {
            if(riskset(i,tr_pos[j])>0)
              {
              if(ttrunc[i]>0)
                {
                H1(j,0) -= (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * mult_hazard(i,tr_pos[j]);
                }
              else
                {
                H1(j,0) -= helpmat(tend[i]-1,0) * mult_hazard(i,tr_pos[j]);
                }
              }
            }
          }
        else
          {
          // time-varying
          helpmat = datamatrix(nrobs,1,0);
          unsigned off = cumtimevar[tr_pos[j]]*t_X.rows();
          for(i=0; i<nrobs; i++)
            {
            if(riskset(i,tr_pos[j])>0)
              {
              for(k=ttrunc[i]; k<tend[i]-1; k++)
                {
                helpmat(i,0) += 0.5*tsteps(k,0) * interactvar(i,fc_pos[j]) * (t_X(k,dm_pos[j])*baseline_tv(i,off+k)+t_X(k+1,dm_pos[j])*baseline_tv(i,off+k+1));
                }
              }
            }
          for(i=0; i<nrobs; i++)
            {
            if(riskset(i,tr_pos[j])>0)
              {
              H1(j,0) += resp(i,tr_pos[j])*X(i,j) - helpmat(i,0) * mult_hazard(i,tr_pos[j]);
              }
            else
              {
              H1(j,0) += resp(i,tr_pos[j])*X(i,j);
              }
            }
          }
        }
      else
        {
        if(!timevarying[tr_pos[j]])
          {
          // not time-varying
          for(i=0; i<nrobs; i++)
            {
            H1(j,0) += resp(i,tr_pos[j])*X(i,j);
            }
          for(i=0; i<nrobs; i++)
            {
            if(riskset(i,tr_pos[j])>0)
              {
              if(ttrunc[i]>0)
                {
                H1(j,0) -= (cumbaseline(tend[i]-1,tr_pos[j])-cumbaseline(ttrunc[i]-1,tr_pos[j])) * X(i,j) * mult_hazard(i,tr_pos[j]);
                }
              else
                {
                H1(j,0) -= cumbaseline(tend[i]-1,tr_pos[j]) * X(i,j) * mult_hazard(i,tr_pos[j]);
                }
              }
            }
          }
        else
          {
          // time-varying
          for(i=0; i<nrobs; i++)
            {
            if(riskset(i,tr_pos[j])>0)
              {
              H1(j,0) += resp(i,tr_pos[j])*X(i,j) - cumbaseline_tv(i,tr_pos[j]) * X(i,j) * mult_hazard(i,tr_pos[j]);
              }
            else
              {
              H1(j,0) += resp(i,tr_pos[j])*X(i,j);
              }
            }
          }
        }
      }

/*ofstream out1("c:\\temp\\H1.raw");
H1.prettyPrint(out1);
out1.close();*/

    // Z

    for(j=0; j<Z.cols(); j++)
      {
      if(isbaselinebeta[xcols+j]==1)
        {
        if(!timevarying[tr_pos[xcols+j]])
          {
          // not time-varying
          helpmat = datamatrix(t_Z.rows(),1,0);
          former=0;
          for(i=0; i<t_Z.rows()-1; i++)
            {
            helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+j])*baseline(i,tr_pos[xcols+j])+t_Z(i+1,dm_pos[xcols+j])*baseline(i+1,tr_pos[xcols+j]));
            former=helpmat(i,0);
            }
          helpmat(t_Z.rows()-1,0) = former;

          for(i=0; i<nrobs; i++)
            {
            H1(xcols+j,0) += resp(i,tr_pos[xcols+j])*Z(i,j);
            }
          for(i=0; i<nrobs; i++)
            {
            if(riskset(i,tr_pos[xcols+j])>0)
              {
              if(ttrunc[i]>0)
                {
                H1(xcols+j,0) -= (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * mult_hazard(i,tr_pos[xcols+j]);
                }
              else
                {
                H1(xcols+j,0) -= helpmat(tend[i]-1,0) * mult_hazard(i,tr_pos[xcols+j]);
                }
              }
            }
          }
        else
          {
          // time-varying
          helpmat = datamatrix(nrobs,1,0);
          unsigned off = cumtimevar[tr_pos[xcols+j]]*t_X.rows();
          for(i=0; i<nrobs; i++)
            {
            if(riskset(i,tr_pos[xcols+j])>0)
              {
              for(k=ttrunc[i]; k<tend[i]-1; k++)
                {
                helpmat(i,0) += 0.5*tsteps(k,0) * interactvar(i,fc_pos[xcols+j]) * (t_Z(k,dm_pos[xcols+j])*baseline_tv(i,off+k)+t_Z(k+1,dm_pos[xcols+j])*baseline_tv(i,off+k+1));
                }
              }
            }
          for(i=0; i<nrobs; i++)
            {
            if(riskset(i,tr_pos[xcols+j])>0)
              {
              H1(xcols+j,0) += resp(i,tr_pos[xcols+j])*Z(i,j) - helpmat(i,0) * mult_hazard(i,tr_pos[xcols+j]);
              }
            else
              {
              H1(xcols+j,0) += resp(i,tr_pos[xcols+j])*Z(i,j);
              }
            }
          }
        }
      else
        {
        if(tr_pos[xcols+j]==-1)
          {
          // global frailty
          for(i=0; i<nrobs; i++)
            {
            for(k=0; k<nrtransitions; k++)
              {
              H1(xcols+j,0) += resp(i,k)*Z(i,j);
              }
            }

          for(i=0; i<nrobs; i++)
            {
            for(k=0; k<nrtransitions; k++)
              {
              if(riskset(i,k)>0)
                {
                if(!timevarying[k])
                  {
                  if(ttrunc[i]>0)
                    {
                    H1(xcols+j,0) -= (cumbaseline(tend[i]-1,k)-cumbaseline(ttrunc[i]-1,k)) * Z(i,j) * mult_hazard(i,k);
                    }
                  else
                    {
                    H1(xcols+j,0) -= cumbaseline(tend[i]-1,k) * Z(i,j) * mult_hazard(i,k);
                    }
                  }
                else
                  {
                  H1(xcols+j,0) -= cumbaseline_tv(i,k) * Z(i,j) * mult_hazard(i,k);
                  }
                }
              }
            }
          }
        else if(!timevarying[tr_pos[xcols+j]])
          {
          // not time-varying
          for(i=0; i<nrobs; i++)
            {
            H1(xcols+j,0) += resp(i,tr_pos[xcols+j])*Z(i,j);
            }
          for(i=0; i<nrobs; i++)
            {
            if(riskset(i,tr_pos[xcols+j])>0)
              {
              if(ttrunc[i]>0)
                {
                H1(xcols+j,0) -= (cumbaseline(tend[i]-1,tr_pos[xcols+j])-cumbaseline(ttrunc[i]-1,tr_pos[xcols+j])) * Z(i,j) * mult_hazard(i,tr_pos[xcols+j]);
                }
              else
                {
                H1(xcols+j,0) -= cumbaseline(tend[i]-1,tr_pos[xcols+j]) * Z(i,j) * mult_hazard(i,tr_pos[xcols+j]);
                }
              }
            }
          }
        else
          {
          // time-varying
          for(i=0; i<nrobs; i++)
            {
            if(riskset(i,tr_pos[xcols+j])>0)
              {
              H1(xcols+j,0) += resp(i,tr_pos[xcols+j])*Z(i,j) - cumbaseline_tv(i,tr_pos[xcols+j]) * Z(i,j) * mult_hazard(i,tr_pos[xcols+j]);
              }
            else
              {
              H1(xcols+j,0) += resp(i,tr_pos[xcols+j])*Z(i,j);
              }
            }
          }
        }
      }
    for(j=0; j<zcols;j++)
      {
      H1(xcols+j,0) -= Qinv(j,0)*beta(xcols+j,0);
      }

  // Fisher-Info for beta

    //clear H-matrix
    H = datamatrix(H.rows(),H.cols(),0);

  // X & X

    for(j=0; j<xcols; j++)
      {
      for(k=j; k<xcols; k++)
        {
        if(tr_pos[j]==tr_pos[k])
          {
          if(!timevarying[tr_pos[j]])
            {
            // not time-varying
            if(isbaselinebeta[j]==1 && isbaselinebeta[k]==1)
              {
              helpmat = datamatrix(t_X.rows(),1,0);
              former=0;
              for(i=0; i<t_X.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_X(i,dm_pos[j])*t_X(i,dm_pos[k])*baseline(i,tr_pos[j])+t_X(i+1,dm_pos[j])*t_X(i+1,dm_pos[k])*baseline(i+1,tr_pos[j]));
                former=helpmat(i,0);
                }
              helpmat(t_X.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(j,k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * mult_hazard(i,tr_pos[j]);
                    }
                  else
                    {
                    H(j,k) += helpmat(tend[i]-1,0) * mult_hazard(i,tr_pos[j]);
                    }
                  }
                }
              }
            else if(isbaselinebeta[j]==1 && isbaselinebeta[k]==0)
              {
              helpmat = datamatrix(t_X.rows(),1,0);
              former=0;
              for(i=0; i<t_X.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_X(i,dm_pos[j])*baseline(i,tr_pos[j])+t_X(i+1,dm_pos[j])*baseline(i+1,tr_pos[j]));
                former=helpmat(i,0);
                }
              helpmat(t_X.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(j,k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * X(i,k) * mult_hazard(i,tr_pos[j]);
                    }
                  else
                    {
                    H(j,k) += helpmat(tend[i]-1,0) * X(i,k) * mult_hazard(i,tr_pos[j]);
                    }
                  }
                }
              }
            else if(isbaselinebeta[j]==0 && isbaselinebeta[k]==1)
              {
              helpmat = datamatrix(t_X.rows(),1,0);
              former=0;
              for(i=0; i<t_X.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_X(i,dm_pos[k])*baseline(i,tr_pos[j])+t_X(i+1,dm_pos[k])*baseline(i+1,tr_pos[j]));
                former=helpmat(i,0);
                }
              helpmat(t_X.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(j,k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * X(i,j) * mult_hazard(i,tr_pos[j]);
                    }
                  else
                    {
                    H(j,k) += helpmat(tend[i]-1,0) * X(i,j) * mult_hazard(i,tr_pos[j]);
                    }
                  }
                }
              }
            else if(isbaselinebeta[j]==0 && isbaselinebeta[k]==0)
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(j,k) += (cumbaseline(tend[i]-1,tr_pos[j])-cumbaseline(ttrunc[i]-1,tr_pos[j])) * X(i,j)*X(i,k) * mult_hazard(i,tr_pos[j]);
                    }
                  else
                    {
                    H(j,k) += cumbaseline(tend[i]-1,tr_pos[j]) * X(i,j)*X(i,k) * mult_hazard(i,tr_pos[j]);
                    }
                  }
                }
              }
            }
          else
            {
            //time-varying
            if(isbaselinebeta[j]==1 && isbaselinebeta[k]==1)
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[j]]*t_X.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[j]) * interactvar(i,fc_pos[k])
                     * (t_X(l,dm_pos[j])*t_X(l,dm_pos[k])*baseline_tv(i,off+l)+t_X(l+1,dm_pos[j])*t_X(l+1,dm_pos[k])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  H(j,k) += helpmat(i,0) * mult_hazard(i,tr_pos[j]);
                  }
                }
              }
            else if(isbaselinebeta[j]==1 && isbaselinebeta[k]==0)
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[j]]*t_X.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[j])
                     * (t_X(l,dm_pos[j])*baseline_tv(i,off+l)+t_X(l+1,dm_pos[j])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  H(j,k) += X(i,k) * helpmat(i,0) * mult_hazard(i,tr_pos[j]);
                  }
                }
              }
            else if(isbaselinebeta[j]==0 && isbaselinebeta[k]==1)
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[j]]*t_X.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[k])
                     * (t_X(l,dm_pos[k])*baseline_tv(i,off+l)+t_X(l+1,dm_pos[k])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  H(j,k) += X(i,j) * helpmat(i,0) * mult_hazard(i,tr_pos[j]);
                  }
                }
              }
            else if(isbaselinebeta[j]==0 && isbaselinebeta[k]==0)
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  H(j,k) += cumbaseline_tv(i,tr_pos[j]) * X(i,j)*X(i,k) * mult_hazard(i,tr_pos[j]);
                  }
                }
              }
            }
          H(k,j) = H(j,k);
          }
        }
      }

  // Z & Z
    for(j=0; j<zcols; j++)
      {
      for(k=j; k<zcols; k++)
        {
        if(tr_pos[xcols+j]==-1 && tr_pos[xcols+k]==-1)
          // both terms belong to global frailty
          {
          for(i=0; i<nrobs; i++)
            {
            for(l=0; l<nrtransitions; l++)
              {
              if(riskset(i,l)>0)
                {
                if(!timevarying[l])
                  {
                  if(ttrunc[i]>0)
                    {
                    H(xcols+j,xcols+k) += (cumbaseline(tend[i]-1,l)-cumbaseline(ttrunc[i]-1,l)) * Z(i,j)*Z(i,k) * mult_hazard(i,l);
                    }
                  else
                    {
                    H(xcols+j,xcols+k) += cumbaseline(tend[i]-1,l) * Z(i,j)*Z(i,k) * mult_hazard(i,l);
                    }
                  }
                else
                  {
                  H(xcols+j,xcols+k) += cumbaseline_tv(i,l) * Z(i,j)*Z(i,k) * mult_hazard(i,l);
                  }
                }
              }
            }
          }
        else if(tr_pos[xcols+j]==-1 && tr_pos[xcols+k]>-1)
          // only term j belongs to global frailty
          {
          if(isbaselinebeta[xcols+k])
            {
            if(!timevarying[tr_pos[xcols+k]])
              {
              helpmat = datamatrix(t_Z.rows(),1,0);
              former=0;
              for(i=0; i<t_Z.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+k])*baseline(i,tr_pos[xcols+k])+t_Z(i+1,dm_pos[xcols+k])*baseline(i+1,tr_pos[xcols+k]));
                former=helpmat(i,0);
                }
              helpmat(t_Z.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+k])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(xcols+j,xcols+k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * Z(i,j) * mult_hazard(i,tr_pos[xcols+k]);
                   }
                  else
                    {
                    H(xcols+j,xcols+k) += helpmat(tend[i]-1,0) * Z(i,j) * mult_hazard(i,tr_pos[xcols+k]);
                    }
                  }
                }
              }
            else
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[xcols+k]]*t_Z.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+k])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[xcols+k])
                     * (t_Z(l,dm_pos[xcols+k])*baseline_tv(i,off+l)+t_Z(l+1,dm_pos[xcols+k])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+k])>0)
                  {
                  H(xcols+j,xcols+k) += Z(i,j) * helpmat(i,0) * mult_hazard(i,tr_pos[xcols+k]);
                  }
                }
              }
            }
          else
            {
            if(!timevarying[tr_pos[xcols+k]])
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+k])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(xcols+j,xcols+k) += (cumbaseline(tend[i]-1,tr_pos[xcols+k])-cumbaseline(ttrunc[i]-1,tr_pos[xcols+k])) * Z(i,j)*Z(i,k) * mult_hazard(i,tr_pos[xcols+k]);
                    }
                  else
                    {
                    H(xcols+j,xcols+k) += cumbaseline(tend[i]-1,tr_pos[xcols+k]) * Z(i,j)*Z(i,k) * mult_hazard(i,tr_pos[xcols+k]);
                    }
                  }
                }
              }
            else
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+k])>0)
                  {
                  H(xcols+j,xcols+k) += cumbaseline_tv(i,tr_pos[xcols+k]) * Z(i,j)*Z(i,k) * mult_hazard(i,tr_pos[xcols+k]);
                  }
                }
              }
            }
          }
        else if(tr_pos[xcols+j]>-1 && tr_pos[xcols+k]==-1)
          // only term k belongs to global frailty
          {
          if(isbaselinebeta[xcols+j])
            {
            if(!timevarying[tr_pos[xcols+j]])
              {
              helpmat = datamatrix(t_Z.rows(),1,0);
              former=0;
              for(i=0; i<t_Z.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+j])*baseline(i,tr_pos[xcols+j])+t_Z(i+1,dm_pos[xcols+j])*baseline(i+1,tr_pos[xcols+j]));
                former=helpmat(i,0);
                }
              helpmat(t_Z.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(xcols+j,xcols+k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * Z(i,k) * mult_hazard(i,tr_pos[xcols+j]);
                   }
                  else
                    {
                    H(xcols+j,xcols+k) += helpmat(tend[i]-1,0) * Z(i,k) * mult_hazard(i,tr_pos[xcols+j]);
                    }
                  }
                }
              }
            else
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[xcols+j]]*t_Z.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[xcols+j])
                     * (t_Z(l,dm_pos[xcols+j])*baseline_tv(i,off+l)+t_Z(l+1,dm_pos[xcols+j])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  H(xcols+j,xcols+k) += Z(i,k) * helpmat(i,0) * mult_hazard(i,tr_pos[xcols+j]);
                  }
                }
              }
            }
          else
            {
            if(!timevarying[tr_pos[xcols+j]])
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(xcols+j,xcols+k) += (cumbaseline(tend[i]-1,tr_pos[xcols+j])-cumbaseline(ttrunc[i]-1,tr_pos[xcols+j])) * Z(i,j)*Z(i,k) * mult_hazard(i,tr_pos[xcols+j]);
                    }
                  else
                    {
                    H(xcols+j,xcols+k) += cumbaseline(tend[i]-1,tr_pos[xcols+j]) * Z(i,j)*Z(i,k) * mult_hazard(i,tr_pos[xcols+j]);
                    }
                  }
                }
              }
            else
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  H(xcols+j,xcols+k) += cumbaseline_tv(i,tr_pos[xcols+j]) * Z(i,j)*Z(i,k) * mult_hazard(i,tr_pos[xcols+j]);
                  }
                }
              }
            }
          }
        else if(tr_pos[xcols+j]==tr_pos[xcols+k])
          // no global frailty term involved
          {
          if(!timevarying[tr_pos[xcols+j]])
            {
            // not time-varying
            if(isbaselinebeta[xcols+j]==1 && isbaselinebeta[xcols+k]==1)
              {
              helpmat = datamatrix(t_Z.rows(),1,0);
              former=0;
              for(i=0; i<t_Z.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+j])*t_Z(i,dm_pos[xcols+k])*baseline(i,tr_pos[xcols+j])+t_Z(i+1,dm_pos[xcols+j])*t_Z(i+1,dm_pos[xcols+k])*baseline(i+1,tr_pos[xcols+j]));
                former=helpmat(i,0);
                }
              helpmat(t_Z.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(xcols+j,xcols+k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * mult_hazard(i,tr_pos[xcols+j]);
                   }
                  else
                    {
                    H(xcols+j,xcols+k) += helpmat(tend[i]-1,0) * mult_hazard(i,tr_pos[xcols+j]);
                    }
                  }
                }
              }
            else if(isbaselinebeta[xcols+j]==0 && isbaselinebeta[xcols+k]==1)
              {
              helpmat = datamatrix(t_Z.rows(),1,0);
              former=0;
              for(i=0; i<t_Z.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+k])*baseline(i,tr_pos[xcols+k])+t_Z(i+1,dm_pos[xcols+k])*baseline(i+1,tr_pos[xcols+k]));
                former=helpmat(i,0);
                }
              helpmat(t_Z.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+k])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(xcols+j,xcols+k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * Z(i,j) * mult_hazard(i,tr_pos[xcols+k]);
                    }
                  else
                    {
                    H(xcols+j,xcols+k) += helpmat(tend[i]-1,0) * Z(i,j) * mult_hazard(i,tr_pos[xcols+k]);
                    }
                  }
                }
              }
            else if(isbaselinebeta[xcols+j]==1 && isbaselinebeta[xcols+k]==0)
              {
              helpmat = datamatrix(t_Z.rows(),1,0);
              former=0;
              for(i=0; i<t_Z.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+j])*baseline(i,tr_pos[xcols+j])+t_Z(i+1,dm_pos[xcols+j])*baseline(i+1,tr_pos[xcols+j]));
                former=helpmat(i,0);
                }
              helpmat(t_Z.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(xcols+j,xcols+k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * Z(i,k) * mult_hazard(i,tr_pos[xcols+j]);
                    }
                  else
                    {
                    H(xcols+j,xcols+k) += helpmat(tend[i]-1,0) * Z(i,k) * mult_hazard(i,tr_pos[xcols+j]);
                    }
                  }
                }
              }
            else if(isbaselinebeta[xcols+j]==0 && isbaselinebeta[xcols+k]==0)
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(xcols+j,xcols+k) += (cumbaseline(tend[i]-1,tr_pos[xcols+j])-cumbaseline(ttrunc[i]-1,tr_pos[xcols+j])) * Z(i,j)*Z(i,k) * mult_hazard(i,tr_pos[xcols+j]);
                    }
                  else
                    {
                    H(xcols+j,xcols+k) += cumbaseline(tend[i]-1,tr_pos[xcols+j]) * Z(i,j)*Z(i,k) * mult_hazard(i,tr_pos[xcols+j]);
                    }
                  }
                }
              }
            }
          else
            // time-varying
            {
            if(isbaselinebeta[xcols+j]==1 && isbaselinebeta[xcols+k]==1)
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[xcols+j]]*t_X.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[xcols+j]) * interactvar(i,fc_pos[xcols+k])
                     * (t_Z(l,dm_pos[xcols+j])*t_Z(l,dm_pos[xcols+k])*baseline_tv(i,off+l)+t_Z(l+1,dm_pos[xcols+j])*t_Z(l+1,dm_pos[xcols+k])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  H(xcols+j,xcols+k) += helpmat(i,0) * mult_hazard(i,tr_pos[xcols+j]);
                  }
                }
              }
            else if(isbaselinebeta[xcols+j]==1 && isbaselinebeta[xcols+k]==0)
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[xcols+j]]*t_X.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[xcols+j])
                     * (t_Z(l,dm_pos[xcols+j])*baseline_tv(i,off+l)+t_Z(l+1,dm_pos[xcols+j])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  H(xcols+j,xcols+k) += Z(i,k) * helpmat(i,0) * mult_hazard(i,tr_pos[xcols+j]);
                  }
                }
              }
            else if(isbaselinebeta[xcols+j]==0 && isbaselinebeta[xcols+k]==1)
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[xcols+j]]*t_X.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[xcols+k])
                     * (t_Z(l,dm_pos[xcols+k])*baseline_tv(i,off+l)+t_Z(l+1,dm_pos[xcols+k])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  H(xcols+j,xcols+k) += Z(i,j) * helpmat(i,0) * mult_hazard(i,tr_pos[xcols+j]);
                  }
                }
              }
            else if(isbaselinebeta[xcols+j]==0 && isbaselinebeta[xcols+k]==0)
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[xcols+j])>0)
                  {
                  H(xcols+j,xcols+k) += cumbaseline_tv(i,tr_pos[xcols+j]) * Z(i,j)*Z(i,k) * mult_hazard(i,tr_pos[xcols+j]);
                  }
                }
              }
            }
          }
        H(xcols+k,xcols+j)=H(xcols+j,xcols+k);
        }
      }

  // X & Z
    for(j=0; j<xcols; j++)
      {
      for(k=0; k<zcols; k++)
        {
        if(tr_pos[xcols+k]==-1)
          // only term k belongs to global frailty
          // cases 'both terms belong to frailty' and 'only term j belongs to
          // frailty are not possible, since tr_pos[j]>-1 for all 0<=j<=xcols
          {
          if(isbaselinebeta[j])
            {
            if(!timevarying[tr_pos[j]])
              {
              helpmat = datamatrix(t_X.rows(),1,0);
              former=0;
              for(i=0; i<t_X.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_X(i,dm_pos[j])*baseline(i,tr_pos[j])+t_X(i+1,dm_pos[j])*baseline(i+1,tr_pos[j]));
                former=helpmat(i,0);
                }
              helpmat(t_X.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(j,xcols+k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * Z(i,k) * mult_hazard(i,tr_pos[j]);
                   }
                  else
                    {
                    H(j,xcols+k) += helpmat(tend[i]-1,0) * Z(i,k) * mult_hazard(i,tr_pos[j]);
                    }
                  }
                }
              }
            else
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[j]]*t_X.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[j])
                     * (t_X(l,dm_pos[j])*baseline_tv(i,off+l)+t_X(l+1,dm_pos[j])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  H(j,xcols+k) += Z(i,k) * helpmat(i,0) * mult_hazard(i,tr_pos[j]);
                  }
                }
              }
            }
          else
            {
            if(!timevarying[tr_pos[j]])
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(j,xcols+k) += (cumbaseline(tend[i]-1,tr_pos[j])-cumbaseline(ttrunc[i]-1,tr_pos[j])) * X(i,j)*Z(i,k) * mult_hazard(i,tr_pos[j]);
                    }
                  else
                    {
                    H(j,xcols+k) += cumbaseline(tend[i]-1,tr_pos[j]) * X(i,j)*Z(i,k) * mult_hazard(i,tr_pos[j]);
                    }
                  }
                }
              }
            else
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  H(j,xcols+k) += cumbaseline_tv(i,tr_pos[j]) * X(i,j)*Z(i,k) * mult_hazard(i,tr_pos[j]);
                  }
                }
              }
            }
          }
        else if(tr_pos[j]==tr_pos[xcols+k])
          // no global frailty term involved
          {
          if(!timevarying[tr_pos[j]])
            {
            // not time-varying
            if(isbaselinebeta[j]==1 && isbaselinebeta[xcols+k]==1)
              {
              helpmat = datamatrix(t_X.rows(),1,0);
              former=0;
              for(i=0; i<t_X.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_X(i,dm_pos[j])*t_Z(i,dm_pos[xcols+k])*baseline(i,tr_pos[j])+t_X(i+1,dm_pos[j])*t_Z(i+1,dm_pos[xcols+k])*baseline(i+1,tr_pos[j]));
                former=helpmat(i,0);
                }
              helpmat(t_X.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(j,xcols+k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * mult_hazard(i,tr_pos[j]);
                    }
                  else
                    {
                    H(j,xcols+k) += helpmat(tend[i]-1,0) * mult_hazard(i,tr_pos[j]);
                    }
                  }
                }
              }
            else if(isbaselinebeta[j]==0 && isbaselinebeta[xcols+k]==1)
              {
              helpmat = datamatrix(t_X.rows(),1,0);
              former=0;
              for(i=0; i<t_X.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_Z(i,dm_pos[xcols+k])*baseline(i,tr_pos[j])+t_Z(i+1,dm_pos[xcols+k])*baseline(i+1,tr_pos[j]));
                former=helpmat(i,0);
                }
              helpmat(t_X.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(j,xcols+k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * X(i,j) * mult_hazard(i,tr_pos[j]);
                    }
                  else
                    {
                    H(j,xcols+k) += helpmat(tend[i]-1,0) * X(i,j) * mult_hazard(i,tr_pos[j]);
                    }
                  }
                }
              }
            else if(isbaselinebeta[j]==1 && isbaselinebeta[xcols+k]==0)
              {
              helpmat = datamatrix(t_X.rows(),1,0);
              former=0;
              for(i=0; i<t_X.rows()-1; i++)
                {
                helpmat(i,0) = former + 0.5*tsteps(i,0)*(t_X(i,dm_pos[j])*baseline(i,tr_pos[j])+t_X(i+1,dm_pos[j])*baseline(i+1,tr_pos[j]));
                former=helpmat(i,0);
                }
              helpmat(t_X.rows()-1,0) = former;

              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(j,xcols+k) += (helpmat(tend[i]-1,0)-helpmat(ttrunc[i]-1,0)) * Z(i,k) * mult_hazard(i,tr_pos[j]);
                    }
                  else
                    {
                    H(j,xcols+k) += helpmat(tend[i]-1,0) * Z(i,k) * mult_hazard(i,tr_pos[j]);
                    }
                  }
                }
              }
            else if(isbaselinebeta[j]==0 && isbaselinebeta[xcols+k]==0)
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  if(ttrunc[i]>0)
                    {
                    H(j,xcols+k) += (cumbaseline(tend[i]-1,tr_pos[j])-cumbaseline(ttrunc[i]-1,tr_pos[j])) * X(i,j)*Z(i,k) * mult_hazard(i,tr_pos[j]);
                    }
                  else
                    {
                    H(j,xcols+k) += cumbaseline(tend[i]-1,tr_pos[j]) * X(i,j)*Z(i,k) * mult_hazard(i,tr_pos[j]);
                    }
                  }
                }
              }
            }
          else
            {
            //time-varying
            if(isbaselinebeta[j]==1 && isbaselinebeta[xcols+k]==1)
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[j]]*t_X.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[j]) * interactvar(i,fc_pos[xcols+k])
                     * (t_X(l,dm_pos[j])*t_Z(l,dm_pos[xcols+k])*baseline_tv(i,off+l)+t_X(l+1,dm_pos[j])*t_Z(l+1,dm_pos[xcols+k])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  H(j,xcols+k) += helpmat(i,0) * mult_hazard(i,tr_pos[j]);
                  }
                }
              }
            else if(isbaselinebeta[j]==1 && isbaselinebeta[xcols+k]==0)
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[j]]*t_X.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[j])
                     * (t_X(l,dm_pos[j])*baseline_tv(i,off+l)+t_X(l+1,dm_pos[j])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  H(j,xcols+k) += Z(i,k) * helpmat(i,0) * mult_hazard(i,tr_pos[j]);
                  }
                }
              }
            else if(isbaselinebeta[j]==0 && isbaselinebeta[xcols+k]==1)
              {
              helpmat = datamatrix(nrobs,1,0);
              unsigned off = cumtimevar[tr_pos[j]]*t_X.rows();
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  for(l=ttrunc[i]; l<tend[i]-1; l++)
                    {
                    helpmat(i,0) += 0.5*tsteps(l,0) * interactvar(i,fc_pos[xcols+k])
                     * (t_Z(l,dm_pos[xcols+k])*baseline_tv(i,off+l)+t_Z(l+1,dm_pos[xcols+k])*baseline_tv(i,off+l+1));
                    }
                  }
                }
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  H(j,xcols+k) += X(i,j) * helpmat(i,0) * mult_hazard(i,tr_pos[j]);
                  }
                }
              }
            else if(isbaselinebeta[j]==0 && isbaselinebeta[xcols+k]==0)
              {
              for(i=0; i<nrobs; i++)
                {
                if(riskset(i,tr_pos[j])>0)
                  {
                  H(j,xcols+k) += cumbaseline_tv(i,tr_pos[j]) * X(i,j)*Z(i,k) * mult_hazard(i,tr_pos[j]);
                  }
                }
              }
            }
          }
        H(xcols+k,j)=H(j,xcols+k);
        }
      }
    H.addtodiag(Qinv,xcols,beta.rows());

    // Fisher-scoring für beta
    beta = betaold + H.solve(H1);


    stop = check_pause();
    if (stop)
      return true;

    Hinv=H.inverse();

    if(!constlambda)
      {
      // transform theta
      for(i=0; i<theta.rows(); i++)
        {
        thetaold(i,0)=log(sqrt(thetaold(i,0)));
        theta(i,0)=log(sqrt(theta(i,0)));
        }

      // Score-Funktion für theta

     for(j=0; j<theta.rows(); j++)
        {
/*        score(j,0)=-1*exp(theta(j,0))*((zcut[j+1]-zcut[j])/exp(theta(j,0))-
                         (Hinv.getBlock(xcols+zcut[j],X.cols()+zcut[j],xcols+zcut[j+1],xcols+zcut[j+1])).trace()/(exp(theta(j,0))*exp(theta(j,0))*exp(theta(j,0)))-
                         (beta.getRowBlock(xcols+zcut[j],xcols+zcut[j+1]).transposed()*beta.getRowBlock(xcols+zcut[j],xcols+zcut[j+1]))(0,0)/(exp(theta(j,0))*exp(theta(j,0))*exp(theta(j,0))));*/
        score(j,0)=-1*((zcut[j+1]-zcut[j])-
                         (Hinv.getBlock(xcols+zcut[j],X.cols()+zcut[j],xcols+zcut[j+1],xcols+zcut[j+1])).trace()/(exp(theta(j,0))*exp(theta(j,0)))-
                         (beta.getRowBlock(xcols+zcut[j],xcols+zcut[j+1]).transposed()*beta.getRowBlock(xcols+zcut[j],xcols+zcut[j+1]))(0,0)/(exp(theta(j,0))*exp(theta(j,0))));
        }

      // Fisher-Info für theta

      for(j=0; j<theta.rows(); j++)
        {
        for(k=j; k< theta.rows(); k++)
          {
//          Fisher(j,k) = 2*exp(theta(j,0))*exp(theta(k,0))*((Hinv.getBlock(xcols+zcut[j],xcols+zcut[k],xcols+zcut[j+1],xcols+zcut[k+1])*Hinv.getBlock(xcols+zcut[k],xcols+zcut[j],xcols+zcut[k+1],xcols+zcut[j+1])).trace())/(exp(theta(j,0))*exp(theta(j,0))*exp(theta(j,0))*exp(theta(k,0))*exp(theta(k,0))*exp(theta(k,0)));
          Fisher(j,k) = 2*((Hinv.getBlock(xcols+zcut[j],xcols+zcut[k],xcols+zcut[j+1],xcols+zcut[k+1])*Hinv.getBlock(xcols+zcut[k],xcols+zcut[j],xcols+zcut[k+1],xcols+zcut[j+1])).trace())/(exp(theta(j,0))*exp(theta(j,0))*exp(theta(k,0))*exp(theta(k,0)));
          Fisher(k,j) = Fisher(j,k);
          }
        }

      //Fisher-scoring für theta

      theta = thetaold + Fisher.solve(score);

      // transform theta back to original parameterisation

      for(i=0; i<theta.rows(); i++)
        {
        theta(i,0) = exp(theta(i,0))*exp(theta(i,0));
        thetaold(i,0) = exp(thetaold(i,0))*exp(thetaold(i,0));
        }

      // stop estimation of small variances

      datamatrix helpmat2=datamatrix(nrtransitions,1,0);
      for(i=0; i<nrtransitions; i++)
        {
        helpmat2(i,0) = eta.norm(i);
        }
      k=0;
      for(i=0; i<nrtransitions; i++)
        {
        for(j=1; j<nrfullconds[i]; j++)
          {
          helpmat = Z.getColBlock(zcut[k],zcut[k+1])*beta.getRowBlock(xcols+zcut[k],xcols+zcut[k+1]);
          stopcrit[k] = helpmat.norm(0)/helpmat2(i,0);
          k++;
          }
        }
      if(globalfrailty)
        {
        helpmat = Z.getColBlock(zcut[k],zcut[k+1])*beta.getRowBlock(xcols+zcut[k],xcols+zcut[k+1]);
        stopcrit[k] = helpmat.norm(0)/helpmat2.mean(0);
        }

      for(i=0; i<theta.rows(); i++)
        {
        if(stopcrit[i]<lowerlim || theta(i,0)>maxvar)
          {
          theta(i,0)=thetaold(i,0);
          if(stopcrit[i]<lowerlim)
            {
            thetastop[i]=1;
            }
          else
            {
            thetastop[i]=-1;
            }
          }
        else
          {
          its[i]=it;
          }
        }
      }

    // compute convergence criteria
    help=betaold.norm(0);
    if(help==0)
      {
      help=0.00001;
      }
    betaold.minus(betaold,beta);
    crit1 = betaold.norm(0)/help;

    help=thetaold.norm(0);
    if(help==0)
      {
      help=0.00001;
      }
    thetaold.minus(thetaold,theta);
    crit2 = thetaold.norm(0)/help;

    // test criterion
    test=((crit1>eps) || (crit2>eps)) && (it<(unsigned)maxit);
    if(it>2)
      {
      test = test && (crit1<maxchange && crit2<maxchange);
      }

    out("  iteration "+ST::inttostring(it)+"\n");
    out("  relative changes in the regression coefficients: "+
         ST::doubletostring(crit1,6)+"\n");
    out("  relative changes in the variance parameters:     "+
         ST::doubletostring(crit2,6)+"\n");
    out("\n");

    // count iteration
    it=it+1;
    }

  if(crit1>=maxchange || crit2>=maxchange)
    {
    out("\n");
    outerror("ERROR: numerical problems due to large relative changes\n");
    outerror("       REML ESTIMATION DID NOT CONVERGE\n");
    out("\n");
    }
  else if(it>=(unsigned)maxit)
    {
    out("\n");
    outerror("WARNING: Number of iterations reached " + ST::inttostring(maxit) + "\n");
    outerror("         REML ESTIMATION DID NOT CONVERGE\n");
    out("\n");
    }
  else
    {
    out("\n");
    out("REML ESTIMATION CONVERGED\n",true);
    out("\n");
    }
  out("ESTIMATION RESULTS:\n",true);
  out("\n");

  datamatrix thetareml(theta.rows(),4,0);
  thetareml.putCol(0,theta);
  H.addtodiag(-Qinv,xcols,beta.rows());
  datamatrix Hhelp = (H*Hinv);
  for(i=0; i<theta.rows(); i++)
    {
    thetareml(i,1)=thetastop[i];
    thetareml(i,2)=its[i];
    thetareml(i,3)=xcut[i+2]-xcut[i+1]+(Hhelp.getBlock(xcols+zcut[i],xcols+zcut[i],xcols+zcut[i+1],xcols+zcut[i+1])).trace();
    }

  k=l=0;
  unsigned m1,m2,l1,l2;
  m1=m2=l1=l2=0;
  for(i=0; i<nrfullconds.size(); i++)
    {
    m2 += fullcond[k]->get_dimX();
    l2++;
    k++;
    for(j=1; j<nrfullconds[i]; j++)
      {
      beta(m1,0) += fullcond[k]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[k],zcut[l],l,false,xcut[k],xcols+zcut[l],0,false,k);
      m2 += fullcond[k]->get_dimX();
      k++;
      l++;
      l2++;
      }
    beta(m1,0) += fullcond[l1]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[l1],0,0,false,xcut[l1],0,0,false,0);
    m1=m2;
    l1=l2;
    }
  if(globalfrailty)
    {
    former = fullcond[k]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[k],zcut[l],l,false,xcut[k],xcols+zcut[l],0,false,k);
    }

  return false;
  }


//------------------------------------------------------------------------------
//----------------------------- Object description -----------------------------
//------------------------------------------------------------------------------

void remlest_multistate::outoptions()
  {
  out("\n");
  out("GENERAL OPTIONS:\n",true);
  out("\n");
  out("  Maxmimum number of iterations:          "+ST::inttostring(maxit)+"\n");
  out("  Termination criterion:                  "+ST::doubletostring(eps,7)+"\n");
  out("  Stopping criterion for small variances: "+ST::doubletostring(lowerlim,6)+"\n");
  out("\n");
  out("RESPONSE DISTRIBUTION:\n",true);
  out("\n");
  out("  Family:                 multistate\n");
  out("  Number of observations: "+ST::inttostring(X.rows())+"\n");
  }

//------------------------------------------------------------------------------
//----------------------------- Writing results --------------------------------
//------------------------------------------------------------------------------

void remlest_multistate::make_plots(ofstream & outtex,ST::string path_batch,
                         ST::string path_splus)
  {
  char hcharu = '_';
  ST::string hstringu = "\\_";

  unsigned i,j,k;
  ST::string pathresult;
  bool stil = false;

// Schleife überprüft, ob es ein fullcond-Object
// gibt, bei dem Effekt gezeichnet werden kann
  MCMC::plotstyles plst;
  for(j=0;j<fullcond.size();j++)
    {
    plst = fullcond[j]->get_plotstyle();
    if(plst != MCMC::noplot)
      stil = true;
    }


  if(stil == true)
    {
//erzeugt File, das Plot-Befehle für Java-Version enthält
    ofstream outbatch(path_batch.strtochar());

//erzeugt File, das SPlus-Befehle zum Plotten enthält
    ofstream outsplus(path_splus.strtochar());

    outsplus << "library(\"BayesX\")\n\n";
/*    outsplus << "# NOTE: 'directory' has to be substituted by the directory"
             << " where the functions are stored \n"
             << endl
             << "# In S-PLUS the file extension in the source command has to be changed"
             << " to '.s' \n"
             << endl
    // einlesen der Source-Files für S-Plus
             << "source(\"'directory'\\\\sfunctions\\\\plotsample.r\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\plotnonp.r\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\plotsurf.r\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\drawmap.r\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\readbndfile.r\")\n" << endl;*/

#if defined(JAVA_OUTPUT_WINDOW)
    out("  --------------------------------------------------------------------------- \n");
    out("\n");
    out("  Batch file for visualizing effects of nonlinear functions is stored in file \n");
    out("  " + path_batch + "\n");
    out("\n");
#endif

    bool stil2 = true;
    for(j=0;j<fullcond.size();j++)  //Schleife überprüft, ob es map-Objekt gibt
      {
      plst = fullcond[j]->get_plotstyle();
      if(plst == MCMC::drawmap || plst == MCMC::drawmapgraph)
        stil2 = false;
      }

    if(stil2 == true)
      {
      out("  --------------------------------------------------------------------------- \n");
      out("\n");
      out("  Batch file for visualizing effects of nonlinear functions ");
      out("  in R is stored in file \n");
      out("  " + path_splus + "\n");
      out("\n");
      }

    if(stil2 == false)
      {
      out("  NOTE: 'input filename' must be substituted by the filename of the boundary-file \n");
      out("\n");
      out("  --------------------------------------------------------------------------- \n");
      out("\n");
      out("  Batch file for visualizing effects of nonlinear functions ");
      out("  in R is stored in file \n");
      out("  " + path_splus + "\n");
      out("\n");
      out("  NOTE: 'input filename' must be substituted by the filename of the boundary-file \n");
      out("\n");
      }


    outbatch << "% usefile " << path_batch << endl;

    // falls andere Quantile gewünscht werden
    double u = fullcond[0]->get_level1();
    double o = fullcond[0]->get_level2();
    ST::string u_str = ST::doubletostring(u,0);
    ST::string o_str = ST::doubletostring(o,0);

    // durchlaufen der Fullconditionals
    j=0;
    for(i=0; i<nrtransitions; i++)
      {
      outtex << "\n\\newpage" << "\n\\noindent {\\bf \\large Plots (Transition " << ST::inttostring(i+1) << "):}" << endl;
      for(k=0; k<nrfullconds[i]; k++, j++)
        {

        // Pfad der Regr.-Ergebnisse
        pathresult = fullcond[j]->get_pathresult();

        // Plotstyle: noplot, plotnonp, drawmap, drawmapgraph
        plst = fullcond[j]->get_plotstyle();

        if (plst != MCMC::noplot)
          {

          // Pfade für ps-, tex-, SPlus-files
          ST::string pathps = pathresult.substr(0, pathresult.length()-4);
          ST::string pathgr = pathps.replaceallsigns('\\', '/');

          char hchar = '\\';
          ST::string hstring = "/";

          ST::string pathps_spl = pathps.insert_string_char(hchar,hstring);
          ST::string pathres_spl = pathresult.insert_string_char(hchar,hstring);

          if (plst == MCMC::plotnonp)
            {
            outbatch << "\n";                // Befehle f. d. batch-file
            outbatch << "dataset _dat" << endl;
            outbatch << "_dat.infile using " << pathresult << endl;
            outbatch << "graph _g" << endl;
            vector<ST::string> varnames = fullcond[j]->get_datanames();
            ST::string xvar = varnames[0];
            outbatch << "_g.plot " << xvar
                     << " pmode ci" << u_str << "lower ci"
                     << o_str.replaceallsigns('.','p') << "lower ci"
                     << o_str.replaceallsigns('.','p') << "upper ci"
                     << u_str.replaceallsigns('.','p') << "upper, "
                     << "title = \"Effect of " << xvar << "\" xlab = " << xvar
                     << " ylab = \" \" " << "outfile = " << pathps
                     << ".ps replace using _dat" << endl;
            outbatch << "drop _dat" << endl;
            outbatch << "drop _g" << endl;
            // Plot-Befehle f. d. SPlus-file
            outsplus << "plotnonp(\"" << pathres_spl << "\")" << endl;
            // Plot-Befehle f. d. tex-file
            ST::string effect = xvar;
            if(varnames.size()>1)
              {
              effect = varnames[1] + "*" + effect;
              }
            outtex << "\n\\begin{figure}[h!]" << endl
                    << "\\centering" << endl
                    << "\\includegraphics[scale=0.6]{" << pathgr << "}" << endl
                    << "\\caption{Non--linear Effect of '" <<
                    effect.insert_string_char(hcharu,hstringu) << "'";
            outtex << "." << endl << "Shown are the posterior modes together with "
                   << u_str << "\\% and " << o_str
                   << "\\% pointwise credible intervals.}" << endl
                   << "\\end{figure}" << endl;
            }
          // für map-Funktionen
          else if (plst == MCMC::drawmap || plst == MCMC::drawmapgraph)
            {
            outbatch << "\n";                 // Befehle f. d. batch-file
            outbatch << "dataset _dat" << endl;
            outbatch << "_dat.infile using " << pathresult << endl;
            outbatch << "map _map" << endl;
            outbatch << "_map.infile using input_filename" << endl;
            outbatch << "graph _g" << endl;
            vector<ST::string> varnames = fullcond[j]->get_datanames();
            ST::string regionvar = varnames[0];
            outbatch << "_g.drawmap " << "pmode" << " " << regionvar
                     << ", map = _map color outfile = " << pathps
                     << "_pmode.ps replace using _dat" << endl;
            outbatch << "_g.drawmap " << "pcat" << u_str << " " << regionvar
                       << ", map = _map nolegend pcat outfile = " << pathps
                       << "_pcat" << u_str << ".ps replace using _dat" << endl;
            outbatch << "_g.drawmap " << "pcat" << o_str << " " << regionvar
                       << ", map = _map nolegend pcat outfile = " << pathps
                       << "_pcat" << o_str << ".ps replace using _dat" << endl;
            outbatch << "drop _dat" << endl;
            outbatch << "drop _g" << endl;
            outbatch << "drop _map" << endl;
            // Plot-Befehle f. d. SPlus-file
            outsplus << "# NOTE: 'input_filename' must be substituted by the "
                     << "filename of the boundary-file \n"
//                     << "# NOTE: choose a 'name' for the map \n" << endl
                     << "m <- read.bnd(\"'input_filename'\")" << endl
                     << "drawmap(data = \"" << pathres_spl << "\", map = m, plotvar = \"pmode\", regionvar = \""
                     << regionvar << "\")" << endl;
            outsplus << "drawmap(data = \"" << pathres_spl << "\", map = m, plotvar = \"pcat" << u_str << "\", regionvar = \""
                    << regionvar << "\", legend = F, pcat = T)" << endl;
            outsplus << "drawmap(data = \"" << pathres_spl << "\", map = m, plotvar = \"pcat" << o_str << "\", regionvar = \""
                      << regionvar << "\", legend = F, pcat = T)" << endl;
              // Plot-Befehle f. d. tex-file
            ST::string effect = regionvar;
            if(varnames.size()>1)
              {
              effect = varnames[1] + "*" + effect;
              }

            if(plst == MCMC::drawmap)
              {
              outtex << "\n\\begin{figure}[h!]" << endl
                     << "\\centering" << endl
                     << "\\includegraphics[scale=0.6]{" << pathgr << "_pmode}"
                     << endl
                     << "\\caption{Non--linear Effect of '" <<
                     effect.insert_string_char(hcharu,hstringu) << "'";
              outtex << ". Shown are the posterior modes.}" << endl
                     << "\\end{figure}" << endl;
              outtex << "\n\\begin{figure}[htb]" << endl
                     << "\\centering" << endl
                     << "\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                     << u_str << "}" << endl
                     << "\\caption{Non--linear Effect of '" << effect << "'";
              outtex << ". Posterior probabilities for a nominal level of "
                     << u_str << "\\%." << endl
                     << "Black denotes regions with strictly negative credible intervals,"
                     << endl
                     << "white denotes regions with strictly positive credible intervals.}"
                     << endl << "\\end{figure}" << endl;
              outtex << "\n\\begin{figure}[htb]" << endl
                     << "\\centering" << endl
                     << "\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                     << o_str << "}" << endl
                     << "\\caption{Non--linear Effect of '" << effect << "'";
              outtex << ". Posterior probabilities for a nominal level of "
                   << o_str << "\\%." << endl
                     << "Black denotes regions with strictly negative credible intervals,"
                     << endl
                     << "white denotes regions with strictly positive credible intervals.}"
                     << endl << "\\end{figure}" << endl;
              }
            else if(plst == MCMC::drawmapgraph)
              {
              outtex << "\n%\\begin{figure}[h!]" << endl
                     << "%\\centering" << endl
                     << "%\\includegraphics[scale=0.6]{" << pathgr << "_pmode}"
                     << endl
                     << "%\\caption{Non--linear Effect of '" <<
                     effect.insert_string_char(hcharu,hstringu) << "'";
              outtex << ". Shown are the posterior modes.}" << endl
                     << "%\\end{figure}" << endl;
              outtex << "\n%\\begin{figure}[htb]" << endl
                     << "%\\centering" << endl
                     << "%\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                     << u_str << "}" << endl
                     << "%\\caption{Non--linear Effect of '" << effect << "'";
              outtex << ". Posterior probabilities for a nominal level of "
                     << u_str << "\\%." << endl
                     << "%Black denotes regions with strictly negative credible intervals,"
                     << endl
                     << "%white denotes regions with strictly positive credible intervals.}"
                     << endl << "%\\end{figure}" << endl;
              outtex << "\n%\\begin{figure}[htb]" << endl
                     << "%\\centering" << endl
                     << "%\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                     << o_str << "}" << endl
                     << "%\\caption{Non--linear Effect of '" << effect << "'";
              outtex << ". Posterior probabilities for a nominal level of "
                     << o_str << "\\%." << endl
                     << "%Black denotes regions with strictly negative credible intervals,"
                     << endl
                     << "%white denotes regions with strictly positive credible intervals.}"
                     << endl << "%\\end{figure}" << endl;
              }
            } // end: else if
          } // end: if
        } // end: for
      }
    }
  }

void remlest_multistate::make_model(ofstream & outtex,
                                    const vector<ST::string> & rnames)
  {
  unsigned i;
  ST::string familyname ="multistate";

  //Anz. Beob. wird übergeben
  unsigned obs = X.rows();

  char charh = '_';
  ST::string stringh = "\\_";
  vector<ST::string> helprname;
  for(i=0; i<nrtransitions; i++)
    {
    helprname.push_back(rnames[i].insert_string_char(charh,stringh));
    }

  //schreibt das Modell und die Priori-annahmen ins Tex-File
  outtex << "\n\\noindent {\\bf \\large Response:}" << endl
         << "\\begin{tabbing}\n";
  outtex << "Number of observations: \\= \\kill" << endl;
  outtex << "Number of observations: \\> " << obs << "\\\\" << endl;
  outtex << "Response Variables: \\> ";
  for(i=0; i < (nrtransitions-1); i++)
    {
    outtex << helprname[i] << ", ";
    }
  outtex << helprname[i] << "\\\\" << endl;

  outtex << "Family: \\> " << familyname << "\\\\" << endl;
  outtex << "\\end{tabbing}" << endl;

  outtex << "\n\\noindent {\\bf \\large Predictor:}\\\\" << endl;
  }

void remlest_multistate::make_predictor(ofstream & outtex)
  {
  unsigned i,j,k;
  ST::string term, term2;
  outtex << endl << "\n\\begin{tabular}{ccp{12cm}}\n";
  j=0;
  for(i=0; i<nrtransitions; i++)
    {
    term2 = fullcond[j]->get_term_symbolic();
    term = "$\\eta_" + ST::inttostring(i+1) + "$ & $=$ & $" + term2;    //linearer Prädiktor wird erweitert
    j++;
    for(k=1; k<nrfullconds[i]; k++, j++)
      {
      term2 = fullcond[j]->get_term_symbolic();
      term = term + " + " + term2;    //linearer Prädiktor wird erweitert
      }
    term = term.insert_after_all_string("^{("+ST::inttostring(i+1)+")}","\\gamma");
    term = term.insert_after_all_string("^{("+ST::inttostring(i+1)+")}","+ f");
    outtex << term << "$\\\\\n";
    }
  outtex << "\n\\end{tabular}\n\\\\ \n\\\\" << endl;

  }

void remlest_multistate::make_prior(ofstream & outtex)
  {
  unsigned i,j,k,l;
  j=0;
  outtex << "\n\\noindent {\\bf \\large Priors:}\\\\" << endl;
  for(i=0; i<nrtransitions; i++)
    {
    outtex << "\n\\noindent {\\bf Transition " << i+1 << ":}\\\\" << endl << "\\\\" << endl;
    for(k=0; k<nrfullconds[i]; k++, j++)
      {
      vector<ST::string> prior = fullcond[j]->get_priorassumptions();
      if(prior.size() != 0)// nur wenn Priors da sind (d.h. Vektor hat Elemente)
        {
        for(l=0;l<prior.size();l++)
          {
          if( k!=0 || l<prior.size()-1)
            {
            outtex << prior[l] << "\\\\" << endl;
            }
          }
        outtex << "\\\\" <<endl;
        }
      }
    }
  }

void remlest_multistate::make_options(ofstream & outtex)
  {
  double l1 = fullcond[0]->get_level1();
  double l2 = fullcond[0]->get_level2();

  //schreibt REML options ins Tex-File
  outtex << "\n\\noindent {\\bf \\large General Options:}" << endl
         << "\\begin{tabbing}" << endl
         << "Levels for credible intervals: \\hspace{2cm}\\= \\\\" << endl
         << "Level 1: \\> " << l1 << "\\\\" << endl
         << "Level 2: \\> " << l2 << "\\\\" << endl
         << "Maxmimum number of iterations: \\> " << maxit << "\\\\" << endl
         << "Termination criterion: \\> " << eps << "\\\\" << endl
         << "Stopping criterion for small variances: \\> " << lowerlim << endl
         << "\\end{tabbing}\n"  << "\\vspace{0.5cm}" <<  endl;
  }

void remlest_multistate::make_fixed_table(ofstream & outtex)
  {
  // falls andere Quantile gewünscht werden
  double u = fullcond[0]->get_level1();
  ST::string u_str = ST::doubletostring(u,0);

  unsigned i,j,k,l;
  j=0;
  for(i=0; i<nrtransitions; i++)
    {
    vector<ST::string> h;
    unsigned r=2;

    // Tabelle im Tex-File mit fixen Effekten
    outtex << "\n\\newpage \n" << endl << "\n\\noindent {\\bf \\large Fixed Effects (Transition " << ST::inttostring(i+1) << "):}\\\\"
         << endl << "\\\\" << endl;

    outtex << "\\begin{tabular}{|r|rrrrr|}" << endl << "\\hline" << endl
           << "Variable & Post. Mode & Std. Dev. & p-value & \\multicolumn{2}{r|}{" << u << "\\% confidence interval}\\\\"
           << endl << "\\hline" << endl;

    h = fullcond[j]->get_results_latex();
    for (l=0;l<h.size();l++)
      {
      r++;
      if (r < 39)
        {
        outtex << h[l] << endl;
        }
      else
        {
        r=1;
        outtex << "\\hline \n\\end{tabular}" << endl;

        outtex << "\n\\newpage \n" << endl
               << "\n\\noindent {\\bf \\large Fixed Effects (continued):}\\\\"
               << endl << "\\\\" << endl;

        outtex << "\\begin{tabular}{|r|rrrrr|}" << endl << "\\hline" << endl
               << "Variable & Post. Mode & Std. Dev. & p-value & \\multicolumn{2}{r|}{" << u << "\\% confidence interval}\\\\"
               << endl << "\\hline" << endl;

        outtex << h[l] << endl;
        }
      }
    outtex << "\\hline \n\\end{tabular}" << endl;
    for(k=0; k<nrfullconds[i]; k++, j++)
      {
      }
    }
  }

void remlest_multistate::make_graphics(const ST::string & title,
                     const ST::string & path_batch,
                     const ST::string & path_tex,
                     const ST::string & path_splus,
                     const vector<ST::string> & rnames)
  {
  ST::string pathresult;                 //Pfad des Ergebnis-Files

  vector<ST::string> distr_results;

 // erzeugt Tex-File
  ofstream outtex(path_tex.strtochar());

  char charh = '_';
  ST::string stringh = "\\_";
  ST::string helptitle = title.insert_string_char(charh,stringh);

  //erzeugt den Kopf des Tex-Files
  outtex << "\\documentclass[a4paper, 12pt]{article}" << endl
         << "\n" << "\\usepackage{graphicx}" << endl
         << "\\parindent0em" << endl
         << "\n\\begin{document}" << endl
         << "\\begin{center}" << endl
         << "\\LARGE{\\bf " << helptitle << "}"
         << endl << "\\end{center} \n\\vspace{1cm}" << endl;

  make_model(outtex,rnames);

  make_predictor(outtex);

  make_prior(outtex);

  make_options(outtex);

/*  outtex << "\n\\noindent {\\bf \\large Model Fit:}" << endl
         << "\\begin{tabbing}\n";
  if(respfamily=="cox" || respfamily=="coxold" || respfamily=="coxinterval" || respfamily=="gaussian")
    {
    outtex << "Degrees of freedom: \\= \\kill" << endl;
    }
  else
    {
    outtex << "GCV (based on deviance residuals): \\= \\kill" << endl;
    }
  outtex << "-2*log-likelihood: \\> " << loglike << "\\\\" << endl;
  outtex << "Degrees of freedom: \\> " << df << "\\\\" << endl;
  outtex << "(conditional) AIC: \\> " << aic << "\\\\" << endl;
  outtex << "(conditional) BIC: \\> " << bic << "\\\\" << endl;
  if(respfamily!="cox" && respfamily!="coxold" && respfamily!="coxinterval" && respfamily!="gaussian")
    {
    outtex << "GCV (based on deviance residuals): \\> " << gcv << "\\\\" << endl;
    }
  else if (respfamily=="gaussian")
    {
    outtex << "GCV: \\> " << gcv << "\\\\" << endl;
    }
  outtex << "\\end{tabbing}" << endl;*/

  make_fixed_table(outtex);

  // Pfade der Files
  //werden im BayesX-Output angegeben
  out("  Files of model summary: \n" , true);
  out("\n");

  make_plots(outtex,path_batch,path_splus);

  out("  --------------------------------------------------------------------------- \n");
  out("\n");
  out("  Latex file of model summaries is stored in file \n");
  out("  " + path_tex + "\n");
  out("\n");
  out("  --------------------------------------------------------------------------- \n");
  out("\n");


  outtex << "\\end{document}" << endl;
  }

bool remlest_multistate::check_pause()
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  Application->ProcessMessages();
  if (Frame->stop)
    {
    return true;
    }

  if (Frame->pause)
    {
    out("\n");
    out("ESTIMATION PAUSED\n");
    out("Click CONTINUE to proceed\n");
    out("\n");

    while (Frame->pause)
      {
      Application->ProcessMessages();
      }

    out("ESTIMATION CONTINUED\n");
    out("\n");
    }
  return false;
#elif defined(JAVA_OUTPUT_WINDOW)
  return adminb_p->breakcommand();
#else
  return false;
#endif
  }

void remlest_multistate::out(const ST::string & s,bool thick,bool italic,
                      unsigned size,int r,int g, int b)
  {
#if defined(BORLAND_OUTPUT_WINDOW)
  ST::string sh = s;
  sh = sh.replaceallsigns('\n',' ');
  if (!Frame->suppoutput)
    Results->ResultsRichEdit->Lines->Append(sh.strtochar());
 if (!(logout->fail()))
    (*logout) << s << flush;
#elif defined(JAVA_OUTPUT_WINDOW)
  ST::string sh = s;
  sh = sh.replaceallsigns('\n',' ');
  sh = sh+"\n";
  if (!adminb_p->get_suppressoutput())
    adminb_p->Java->CallVoidMethod(adminb_p->BayesX_obj, adminb_p->javaoutput,
    adminb_p->Java->NewStringUTF(sh.strtochar()),thick,italic,size,r,g,b);
  if (!(logout->fail()))
    (*logout) << s << flush;
#else
  cout << s << flush;
  if (!(logout->fail()))
    (*logout) << s << flush;
#endif
  }


void remlest_multistate::outerror(const ST::string & s)
  {
  out(s,true,true,12,255,0,0);
  }






















