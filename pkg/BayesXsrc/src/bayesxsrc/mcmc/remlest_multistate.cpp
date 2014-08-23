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



#include "remlest_multistate.h"

#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"
#endif

//------------------------------------------------------------------------------
//----------------------------- Constructor ------------------------------------
//------------------------------------------------------------------------------

/*remlest_multistate::remlest_multistate(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adb,
#endif
vector<MCMC::FULLCOND*> & fc,datamatrix & re,
                const ST::string & family, const ST::string & ofile,
                const int & maxiter, const double & lowerlimit,
                const double & epsi, const double & maxch,
                const vector<unsigned> & nrfullconds,
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

    nrtransitions = re.cols();

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
      zcuttrans.push_back(zcut[k]);
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

    }

//------------------------------------------------------------------------------
//----------------------------- REML estimation --------------------------------
//------------------------------------------------------------------------------

  // Function: estimate
  // Task: Perform REML-estimation for multi state models
  //       returns true if an error or user break occured

bool remlest_multistate::estimate(const datamatrix resp,
               const datamatrix & offset, const datamatrix & weight)
  {
  unsigned i, j, k, l;
  double help;

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
  unsigned nrobs=Z.rows();
  unsigned xcols = X.cols();
  unsigned zcols = Z.cols();

  // Criteria to detemine convergence
  double crit1=1;                //relative changes in regression parameters
  double crit2=1;                //relative changes in variance parameters
  bool test=true;

  vector<double>stopcrit(theta.rows(),10);
  vector<int>its(theta.rows(),0);
  vector<int>signs(theta.rows(),1);

  // Matrix containing the inverse covariance matrix of the random effects
  statmatrix<double>Qinv(zcols,1,0);

  // Inzidenzmatrix, die für jeden Wert in fullcond bzw. beta angibt, ob er zur Baseline-HR beiträgt
  vector<int>isbaseline(fullcond.size(),0);
  int nrbaseline=0;
  for(i=0; i<fullcond.size(); i++)
    {
    if(fullcond[i]->is_baseline()==true)
      {
      isbaseline[i]=1;
      nrbaseline++;
      }
    }

  vector<int>isbaselinebeta(beta.rows(),0);
  vector<int>fc_pos(beta.rows(),0);
  vector<int>dm_pos(beta.rows(),0);
  l=k=0;
  for(i=0; i<fullcond.size(); i++)
    {
    if(isbaseline[i]==1)
      {
      for(j=zcut[i-1]; j<zcut[i]; j++, k++)
        {
        isbaselinebeta[xcols+j]=1;
        fc_pos[xcols+j]=l;
        dm_pos[xcols+j]=k;
        }
      l++;
      }
    }
  l=k=0;
  for(i=0; i<fullcond.size(); i++)
    {
    if(isbaseline[i]==1)
      {
      for(j=xcut[i]; j<xcut[i+1]; j++, k++)
        {
        isbaselinebeta[j]=1;
        fc_pos[j]=l;
        dm_pos[j]=k;
        if(xcut[i+1]==xcut[i]+1)
          {
          dm_pos[j]=1;
          }
        }
      l++;
      }
    }

// Matrices and variables for baseline effects
  datamatrix tsteps;
  datamatrix t_X;
  datamatrix t_Z;
  statmatrix<unsigned> tstart(nrobs,nrtransitions,0);
  statmatrix<unsigned> tend(nrobs,nrtransitions,0);
  statmatrix<unsigned> ttrunc(nrobs,nrtransitions,0);
  datamatrix interactvar(nrobs,1,0);
  statmatrix<int> index(nrobs,nrtransitions,0);
  datamatrix tstepshelp;
  datamatrix t_Xhelp;
  datamatrix t_Zhelp;
  vector<unsigned> tstarthelp;
  vector<unsigned> tendhelp;
  vector<unsigned> ttrunchelp;
  statmatrix<int> indexhelp;
  j=0;
  for(i=0; i<fullcond.size(); i++)
    {
    if(isbaseline[i]==1)
      {
      fullcond[i]->initialize_baseline(0,t_Xhelp,t_Zhelp,tstarthelp,tendhelp,ttrunchelp,interactvar,tstepshelp,indexhelp);
      if(j==0)
        {
        tsteps = datamatrix(tstepshelp.rows(),nrtransitions,0);
        t_X = datamatrix(t_X.rows(),t_X.cols()*nrtransitions,0);
        t_Z = datamatrix(t_Z.rows(),t_Z.cols()*nrtransitions,0);
        index = statmatrix<int>(indexhelp.rows(),nrtransitions,0);
        }
      tsteps.putCol(j,tstepshelp);
      t_X.putColBlock(j*t_Xhelp.cols(),(j+1)*t_Xhelp.cols(),t_Xhelp);
      t_Z.putColBlock(j*t_Zhelp.cols(),(j+1)*t_Zhelp.cols(),t_Zhelp);
      index.putCol(j,indexhelp);
      for(k=0; k<nrobs; k++)
        {
        tstart(k,j) = tstarthelp[k];
        tend(k,j) = tendhelp[k];
        ttrunc(k,j) = ttrunchelp[k];
        }
      }
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
  }

void remlest_multistate::make_model(ofstream & outtex, const ST::string & rname)
  {
  }

void remlest_multistate::make_predictor(ofstream & outtex)
  {
  }

void remlest_multistate::make_prior(ofstream & outtex)
  {
  }

void remlest_multistate::make_options(ofstream & outtex)
  {
  }

void remlest_multistate::make_fixed_table(ofstream & outtex)
  {
  }

void remlest_multistate::make_graphics(const ST::string & title,
                     const ST::string & path_batch,
                     const ST::string & path_tex,
                     const ST::string & path_splus,
                     const ST::string & rname,
                     const bool & dispers)
  {
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
  (*logout) << s << flush;
#endif
  }


void remlest_multistate::outerror(const ST::string & s)
  {
  out(s,true,true,12,255,0,0);
  }
*/












