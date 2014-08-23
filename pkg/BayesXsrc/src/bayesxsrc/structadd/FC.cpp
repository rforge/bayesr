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



#include "FC.h"
#include "clstring.h"

using std::ofstream;
using std::ifstream;
using std::ios;

//------------------------------------------------------------------------------
//----------------- CLASS: FC implementation of member functions ---------------
//------------------------------------------------------------------------------


namespace MCMC
{

void FC::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {

  }

FC::FC(void)
  {
  }


void FC::check_errors(void)
  {
  errors = false;
  }


FC::FC(GENERAL_OPTIONS * o,const ST::string & t,const unsigned & rows,
                   const unsigned & cols,const ST::string & fp)
  {

  nosamples = false;
  nosamplessave = false;
  optionsp = o;

  title = t;

  samplepath = fp;

  setbeta(rows,cols,0);

  acceptance = 0;
  nrtrials = 0;
  outsidelinpredlimits = 0;

  column = 0;

  addon = 0;

  meaneffect = 0;

  check_errors();
  }


FC::FC(const FC & m)
  {

  nosamples = m.nosamples;
  nosamplessave = m.nosamplessave;

  optionsp = m.optionsp;

  title = m.title;

  samplepath = m.samplepath;

  priorassumptions = m.priorassumptions;

  beta = m.beta;
  beta_mode = m.beta_mode;
  betamean = m.betamean;
  betas2 = m.betas2;
  betavar = m.betavar;
  betamin = m.betamin;
  betamax = m.betamax;

  betaqu_l1_lower = m.betaqu_l1_lower;
  betaqu_l2_lower = m.betaqu_l2_lower;
  betaqu_l1_upper = m.betaqu_l1_upper;
  betaqu_l2_upper = m.betaqu_l2_upper;
  betaqu50 = m.betaqu50;

  betameanold = m.betameanold;
  betavarold = m.betavarold;
  betaminold = m.betaminold;
  betamaxold = m.betamaxold;

  sampled_beta = m.sampled_beta;

  addon = m.addon;

  acceptance = m.acceptance;
  nrtrials = m.nrtrials;
  outsidelinpredlimits = m.outsidelinpredlimits;

  column = m.column;

  meaneffect = m.meaneffect;

  errors = m.errors;
  errormessages = m.errormessages;

  }



const FC & FC::operator=(const FC & m)
  {

  if (this==&m)
	 return *this;

  nosamples = m.nosamples;
  nosamplessave = m.nosamplessave;

  optionsp = m.optionsp;

  title = m.title;

  samplepath = m.samplepath;

  priorassumptions = m.priorassumptions;

  beta = m.beta;
  beta_mode = m.beta_mode;
  betamean = m.betamean;
  betas2 = m.betas2;
  betavar = m.betavar;
  betamin = m.betamin;
  betamax = m.betamax;

  betaqu_l1_lower = m.betaqu_l1_lower;
  betaqu_l2_lower = m.betaqu_l2_lower;
  betaqu_l1_upper = m.betaqu_l1_upper;
  betaqu_l2_upper = m.betaqu_l2_upper;
  betaqu50 = m.betaqu50;

  betameanold = m.betameanold;
  betavarold = m.betavarold;
  betaminold = m.betaminold;
  betamaxold = m.betamaxold;

  sampled_beta = m.sampled_beta;

  addon = m.addon;


  acceptance = m.acceptance;
  nrtrials = m.nrtrials;
  outsidelinpredlimits = m.outsidelinpredlimits;

  column = m.column;

  meaneffect = m.meaneffect;

  errors = m.errors;
  errormessages = m.errormessages;

  return *this;

  }



void FC::setbeta(const unsigned & rows,const unsigned & cols,
                       const double & v)
  {

  beta = datamatrix(rows,cols,v);
  betamean = datamatrix(rows,cols,0);
  beta_mode = datamatrix(rows,cols,0);
  betas2 = datamatrix(rows,cols,0);
  betavar = datamatrix(rows,cols,0);
  betamin = datamatrix(rows,cols,0);
  betamax = datamatrix(rows,cols,0);
  betaqu_l1_lower = datamatrix(rows,cols,0);
  betaqu_l2_lower = datamatrix(rows,cols,0);
  betaqu_l1_upper = datamatrix(rows,cols,0);
  betaqu_l2_upper = datamatrix(rows,cols,0);
  betaqu50 = datamatrix(rows,cols,0);
  betameanold = datamatrix(rows,cols,0);
  betavarold = datamatrix(rows,cols,0);
  betaminold = datamatrix(rows,cols,0);
  betamaxold = datamatrix(rows,cols,0);

  }


void FC::setbetavalue(const unsigned & row,const unsigned & col,
                            const double & v)
  {
  beta(row,col) = v;
  }


void FC::setbeta(const datamatrix & betanew)
  {
  beta = betanew;
  betamean = datamatrix(beta.rows(),beta.cols(),0);
  beta_mode = datamatrix(beta.rows(),beta.cols(),0);
  betas2 = datamatrix(beta.rows(),beta.cols(),0);
  betavar = datamatrix(beta.rows(),beta.cols(),0);
  betamin = datamatrix(beta.rows(),beta.cols(),0);
  betamax = datamatrix(beta.rows(),beta.cols(),0);
  betaqu_l1_lower = datamatrix(beta.rows(),beta.cols(),0);
  betaqu_l2_lower = datamatrix(beta.rows(),beta.cols(),0);
  betaqu_l1_upper = datamatrix(beta.rows(),beta.cols(),0);
  betaqu_l2_upper = datamatrix(beta.rows(),beta.cols(),0);
  betaqu50 = datamatrix(beta.rows(),beta.cols(),0);
  betameanold = datamatrix(beta.rows(),beta.cols(),0);
  betavarold = datamatrix(beta.rows(),beta.cols(),0);
  betaminold = datamatrix(beta.rows(),beta.cols(),0);
  betamaxold = datamatrix(beta.rows(),beta.cols(),0);

  }


void FC::readsample(datamatrix & sample,const unsigned & nr,
                          const unsigned & col) const
  {
  unsigned nrpar = beta.cols()*beta.rows();

  unsigned i;

  unsigned s = sample.cols();

  double* work=sample.getV()+col;

  double * sampled_betap = sampled_beta.getV()+nr;
  for (i=0;i<optionsp->samplesize;i++,work+=s,sampled_betap+=nrpar)
    {
    *work = *sampled_betap;
    }

  }


void FC::readsample2(datamatrix & b,const unsigned & nr) const
  {

  unsigned nrpar = beta.cols()*beta.rows();
  unsigned size = sizeof(double);
  ifstream in;
  in.open(samplepath.strtochar(),ios::binary);
  double * work = b.getV();
  unsigned i;
  for (i=0;i<nrpar;i++,work++)
    {
    in.seekg(size*nrpar*nr + i*size);
    in.read((char*) work,size);
    }

  }


void FC::readsample3(datamatrix & b) const
  {

  if ((nosamples == false) && (nosamplessave == false))
    {
    b.assign(sampled_beta);
    } // end: if (nosamples == false)

  }


datamatrix FC::compute_autocorr(const unsigned & lag,const unsigned & row,
                                      const unsigned & col) const
  {
  if ((nosamples == false) && (nosamplessave == false))
    {

    unsigned nr = row*(beta.cols())+col;

    return sampled_beta.autocorr(1,lag,nr);

    }
  else
    return datamatrix(1,1);
  }


double FC::compute_autocorr_single(const unsigned & lag,const unsigned & row,
                                      const unsigned & col) const
  {
  if ((nosamples == false) && (nosamplessave == false))
    {

    unsigned nr = row*(beta.cols())+col;

    return sampled_beta.autocorr(lag,nr);

    }
  else
    return 0;
  }


void FC::compute_autocorr(const ST::string & path, unsigned lag) const
  {
  if (nosamples == false)
    {
    unsigned k,i,c,l,r;
    unsigned nrpar;
    bool misstot;
    double min,max,mean,autoc;

    ofstream out(path.strtochar());

    optionsp->out(path);
    optionsp->out("\n");

    out << "lag  ";

    for (k=0;k<beta.cols();k++)
      for(i=0;i<beta.rows();i++)
        {
        if (beta.cols() == 1)
          out << "b_" << (i+1) << " ";
        else
          out << "b_" << (i+1) << "_" << (k+1) << " ";
        }

    out  << "b_min " << "b_mean " << "b_max " << endl;

    misstot = false;

    for(l=1;l<=lag;l++)
      {

      nrpar = 0;

      out << l << "  ";

      min = 1;
      max = -1;
      mean = 0;

      for(c=0;c<beta.cols();c++)
        for (r=0;r<beta.rows();r++)
          {
          autoc = compute_autocorr_single(l,r,c);
          if ( (autoc <= 1) && (autoc >= -1) )
            {
            nrpar++;
            if (autoc < min)
              min = autoc;
            if (autoc > max)
              max = autoc;
            mean += autoc;
            out << autoc << "  ";
            }
          else
            {
            out << "NA  " << endl;
            misstot = true;
            }
          }

        out << min << "  ";
        out << max << "  ";
        out << mean/nrpar << "  ";
        out << endl;
      }  // end: for(l=0;l<lag;l++)

    if (misstot==true)
      {
      optionsp->out("WARNING: There were undefined autocorrelations\n",true,true);
      optionsp->out("\n");
      }
    }
  }


void FC::compute_autocorr_all(const ST::string & path,
                              unsigned lag, ofstream & outg) const
  {
  compute_autocorr(path,lag);
  outg << "_d.infile using " << path << endl;
  ST::string pathps = path.substr(0,path.length()-4) + ".ps";
  outg << "_g.plotautocor , outfile=" <<  pathps.strtochar() <<  " using _d" << endl;
  outg << endl;
  }


void FC::get_samples(const ST::string & filename,ofstream & outg) const
  {
  if ((nosamples == false) && (nosamplessave == false))
    {

    unsigned i,j,k;

    unsigned nrpar = beta.rows()*beta.cols();

    ofstream out(filename.strtochar());
    assert(!out.fail());

    out << "intnr " << " ";
    if (beta.cols() > 1)
      {
      for (j=0;j<beta.rows();j++)
        for(k=0;k<beta.cols();k++)
          out << "b_" << (j+1) << "_" << (k+1) << " ";
      }
    else
      {
      for (j=0;j<nrpar;j++)
        out << "b_" << (j+1) << " ";
      }

    out << endl;

    double * sampled_betap = sampled_beta.getV();
    for(i=0;i<optionsp->samplesize;i++)
      {
      out << (i+1) << " ";
      for (j=0;j<nrpar;j++,sampled_betap++)
        out << *sampled_betap << " ";

      out << endl;
      }

    out.close();

    optionsp->out(filename + "\n");

    outg << "_d.infile using " << filename << endl;
    ST::string pathps = filename.substr(0,filename.length()-4) + ".ps";
    outg << "_g.plotsample , outfile=" <<  pathps.strtochar() <<  " using _d" << endl;
    outg << endl;

    } // end: if (nosamples == false)
  }



void FC::update(void)
  {
  double diffmean;
  double diffvar;
  double diffmin;
  double diffmax;
  double normold;
  double rate;

  if(
     (optionsp->nriter > optionsp->burnin)
     &&
     ((optionsp->nriter-optionsp->burnin-1) % (optionsp->step) == 0)
    )
    {

    register unsigned i,j;
    double* workbeta = beta.getV();
    double* workbetamean = betamean.getV();
    double* workbetas2 = betas2.getV();
    double* workbetavar = betavar.getV();
    double* workbetamin = betamin.getV();
    double* workbetamax = betamax.getV();


    unsigned samplesize = optionsp->samplesize;

    if ((samplesize==1) && (nosamplessave == false))
     {
     unsigned ssize = optionsp->compute_samplesize();
     unsigned npar = beta.rows()*beta.cols();
     sampled_beta = datamatrix(ssize,npar,0);
     }


    double betatransform;

    if (nosamplessave==false)
      {
      double * sbetap = sampled_beta.getV() +  (samplesize-1)*sampled_beta.cols();

      for(i=0;i<beta.rows();i++)
        {
        for (j=0;j<beta.cols();j++,workbeta++,workbetamean++,workbetas2++,
             workbetavar++,workbetamin++,workbetamax++,sbetap++)

          {

          betatransform = (*workbeta)+addon;

          // storing sampled parameters in binary mode

          *sbetap = betatransform;

          // updating betamean
          if (samplesize==1)
            *workbetamean = betatransform;
          else
            *workbetamean = (1.0/(samplesize))*
                         ((samplesize-1)*(*workbetamean) + betatransform);

          // updating betavar
          if (samplesize==1)
            *workbetas2 = betatransform*betatransform;
          else
            *workbetas2 += betatransform*betatransform;
          *workbetavar = (1.0/samplesize)*
                         (*workbetas2)-(*workbetamean)*(*workbetamean);

          // updating betamin, betamax
          if (samplesize==1)
            {
            *workbetamin = betatransform;
            *workbetamax = betatransform;
            }
          else
            {
            if (betatransform < *workbetamin)
              *workbetamin = betatransform;
            if (betatransform > *workbetamax)
              *workbetamax = betatransform;
            }

          }

        }  // end: for i=0; ...

      }
    else
      {

      for(i=0;i<beta.rows();i++)
        {
        for (j=0;j<beta.cols();j++,workbeta++,workbetamean++,workbetas2++,
             workbetavar++,workbetamin++,workbetamax++)

          {

          betatransform = (*workbeta)+addon;

          // updating betamean
          if (samplesize==1)
            *workbetamean = betatransform;
          else
            *workbetamean = (1.0/(samplesize))*
                         ((samplesize-1)*(*workbetamean) + betatransform);

          // updating betavar
          if (samplesize==1)
            *workbetas2 = betatransform*betatransform;
          else
            *workbetas2 += betatransform*betatransform;
          *workbetavar = (1.0/samplesize)*
                         (*workbetas2)-(*workbetamean)*(*workbetamean);

          // updating betamin, betamax
          if (samplesize==1)
            {
            *workbetamin = betatransform;
            *workbetamax = betatransform;
            }
          else
            {
            if (betatransform < *workbetamin)
              *workbetamin = betatransform;
            if (betatransform > *workbetamax)
              *workbetamax = betatransform;
            }

          }

        }  // end: for i=0; ...

      }

    } // end: if ( (nriter > optionsp->burnin) && ((nriter-burnin-1) % optionsp->step == 0) )

  if(
        (optionsp->nriter > optionsp->burnin)
     &&
        ((optionsp->nriter-optionsp->burnin) %
          optionsp->nrbetween == 0)
     && (title!="")
    )
    {

    optionsp->out("\n");
    optionsp->out("  " + title + "\n");
    optionsp->out("\n");

    if (nrtrials == 0)
      rate = (double(acceptance)/double(optionsp->nriter) )*100;
    else
      rate = (double(acceptance)/double(nrtrials) )*100;
    optionsp->out("  Acceptance rate:    "  + ST::doubletostring(rate,4)
                  + " %\n");
    optionsp->out("\n");


    normold = norm(betameanold);
    if (normold==0)
	  diffmean = MAXDOUBLE;
    else
      diffmean = norm(betamean-betameanold)/normold;

    normold = norm(betavarold);
    if (normold==0)
	  diffvar = MAXDOUBLE;
	else
      diffvar = norm(betavar-betavarold)/normold;

    normold = norm(betaminold);
    if (normold==0)
	  diffmin = MAXDOUBLE;
    else
      diffmin = norm(betamin-betaminold)/normold;

    normold = norm(betamaxold);
    if (normold==0)
	  diffmax = MAXDOUBLE;
    else
      diffmax = norm(betamax-betamaxold)/normold;

    optionsp->out("  Relative Changes in  \n");
    optionsp->out("\n");

    optionsp->out("  Mean:               " + ST::doubletostring(diffmean,6) + "\n");
    optionsp->out("  Variance:           " + ST::doubletostring(diffvar,6) + "\n");
    optionsp->out("  Minimum:            " + ST::doubletostring(diffmin,6) + "\n");
    optionsp->out("  Maximum:            " + ST::doubletostring(diffmax,6) + "\n");

    optionsp->out("\n");
    optionsp->out("\n");

    betameanold.assign(betamean);
    betavarold.assign(betavar);
    betaminold.assign(betamin);
    betamaxold.assign(betamax);


    }  // end: if ((it > burnin) && ((it-burnin) % nrbetween == 0))

  }





bool FC::posteriormode(void)
  {

  double diffmean;
  double normold;

  normold = norm(betameanold);

  if (normold==0)
      diffmean = MAXDOUBLE;
  else
    diffmean = norm(beta-betameanold)/normold;

  betameanold.assign(beta);

  posteriormode_betamean();

  if (diffmean <= 0.00001)
    {
    return true;
    }
  else
    return false;

  }


void FC::posteriormode_betamean(void)
  {

  unsigned i,j;
  double* workbetamean = betamean.getV();
  double* workbeta = beta.getV();

  for(i=0;i<beta.rows();i++)
    for (j=0;j<beta.cols();j++,workbetamean++,workbeta++)
      *workbetamean =  (*workbeta)+addon;

  }


void FC::outresults_acceptance(void)
  {
  if (optionsp->samplesize > 0)
    {
    double rate;
    double rateoutside;
    rateoutside = (double(outsidelinpredlimits)/double(optionsp->nriter))*100;
    if (nrtrials == 0)
      {
      rate = (double(acceptance)/double(optionsp->nriter))*100;
      }
    else
      {
      rate = (double(acceptance)/double(nrtrials))*100;
      }
    optionsp->out("    Acceptance rate:    "  + ST::doubletostring(rate,4) + " %\n");
    if (optionsp->saveestimation)
      {
      if (outsidelinpredlimits > 0)
        optionsp->out("    WARNING: In " + ST::doubletostring(rateoutside,4) +
                      " percent of iterations the predictor was outside limits\n");

      }
    optionsp->out("\n");
    }
  }



double FC::simconfBand(bool l1)
  {
  unsigned i,j;

  datamatrix maxscaling(sampled_beta.rows(),1,1);

  double shelp;
  double * lop;
  double * upp;
  double * meanp;

  for (i=0;i<sampled_beta.rows();i++)
    {
    if (l1==true)
      {
      lop = betaqu_l1_lower.getV();
      upp = betaqu_l1_upper.getV();
      }
    else
      {
      lop = betaqu_l2_lower.getV();
      upp = betaqu_l2_upper.getV();
      }

    meanp = betamean.getV();

    for (j=0;j<sampled_beta.cols();j++,lop++,upp++,meanp++)
      {

      if ( sampled_beta(i,j) < *lop)  // funktion drunter
        {
        shelp = (*meanp -  sampled_beta(i,j)) / (*meanp - *lop);
        if (shelp > maxscaling(i,0))
          maxscaling(i,0) = shelp;
        }
      else if ( sampled_beta(i,j) > *upp ) // funktion drüber
        {
        shelp = (sampled_beta(i,j) -  *meanp) / (*upp - *meanp);
        if (shelp > maxscaling(i,0))
          maxscaling(i,0) = shelp;
        }

      }

    }


  // TEST

  /*
  double s;
  if (l1==true)
    s= maxscaling.quantile(optionsp->level1,0);
  else
    s= maxscaling.quantile(optionsp->level2,0);

  int rate = 0;

  bool o = true;
  for (i=0;i<sampled_beta.rows();i++)
    {
    if (l1==true)
      {
      lop = betaqu_l1_lower.getV();
      upp = betaqu_l1_upper.getV();
      }
    else
      {
      lop = betaqu_l2_lower.getV();
      upp = betaqu_l2_upper.getV();
      }

    meanp = betamean.getV();
    o = true;
    for (j=0;j<sampled_beta.cols();j++,lop++,upp++,meanp++)
      {
      if (sampled_beta(i,j) < *meanp - (*meanp-*lop)*s)
        o = false;

      if (sampled_beta(i,j) > *meanp + (*upp-*meanp)*s)
        o = false;

      }
    if (o==false)
      rate++;

    }

  ofstream out("c:\\bayesx\\testh\\results\\f.raw");
  sampled_beta.prettyPrint(out);
  ENDE: TEST
  */

  if (l1==true)
    return  maxscaling.quantile(optionsp->level1,0);
  else
    return  maxscaling.quantile(optionsp->level2,0);

  }


void FC::outresults_help(ofstream & out_stata, ofstream & out_R,
                    const ST::string & pathresults,
                    const vector<ST::string> & datanames)
  {

  unsigned i;
  unsigned nr = beta.rows();

  vector<ST::string> vnames;
  if (datanames.size() == beta.rows())
    vnames = datanames;
  else
    {
    nr = datanames.size();                            //Kommentar Susanne: das hab ich hier eigefügt
    unsigned i;
    for (i=0;i<nr;i++)
//      vnames.push_back("_x_" + ST::inttostring(i)); //Kommentar Susanne: das hab ich auskommentiert
        vnames.push_back(datanames[i]);               //Kommentar Susanne: das hab ich ersetzt
    }

//  outresults(out_stata,out_R,"");

  ST::string l1 = ST::doubletostring(optionsp->lower1,4);
  ST::string l2 = ST::doubletostring(optionsp->lower2,4);
  ST::string u1 = ST::doubletostring(optionsp->upper1,4);
  ST::string u2 = ST::doubletostring(optionsp->upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  ofstream outp(pathresults.strtochar());

  if (pathresults.isvalidfile() != 1)
    outp << "paramnr varname pmean pstd pqu" << l1 << " pqu" << l2 <<
              " pmed pqu" << u1 << " pqu" << u2 << " pcat" << optionsp->level1
               << " pcat" << optionsp->level2 << endl;

  optionsp->out("\n");

  ST::string l;
  int maxvarnamelength = 0;
  int len;

  for(i=0;i<nr;i++)
    {
    len = vnames[i].length();

    if (len > maxvarnamelength)
      maxvarnamelength = len;
    }

  if (maxvarnamelength>10)
    l = ST::string(' ',maxvarnamelength-4);
  else
    l = "  ";

  ST::string help =  ST::doubletostring(optionsp->lower1,4) + "% quant.";
  ST::string levell = help + ST::string(' ',15-help.length());
  help = ST::doubletostring(optionsp->upper2,4) + "% quant.";
  ST::string levelu = help + ST::string(' ',15-help.length());

  optionsp->out("    Variable" + l +
                    "mean           " +
                    "Std. Dev.      " +
                    levell +
                    "median         " +
                    levelu + "\n");

  ST::string mean;
  ST::string std;
  ST::string qu10;
  ST::string qu50;
  ST::string qu90;

  double m,stddouble;

  unsigned nsp;

  for (i=0;i<nr;i++)
    {

    if (maxvarnamelength  > 10)
      nsp = 4+maxvarnamelength-vnames[i].length();
    else
      nsp = 10-vnames[i].length();

    m= betamean(i,0);

    if (betavar(i,0) < 0.0000000000001)
      stddouble = 0;
    else
      stddouble = sqrt(betavar(i,0));

    if (pathresults.isvalidfile() != 1)
      {
      outp << (i+1) << "   ";
      outp << vnames[i] << "   ";
      outp << m << "   ";
      outp << stddouble << "   ";

      outp << betaqu_l1_lower(i,0) << "   ";
      outp << betaqu_l2_lower(i,0) << "   ";
      outp << betaqu50(i,0) << "   ";
      outp << betaqu_l2_upper(i,0) << "   ";
      outp << betaqu_l1_upper(i,0) << "   ";
      if (betaqu_l1_lower(i,0) > 0)
        outp << "1   ";
      else if (betaqu_l1_upper(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";

      if (betaqu_l2_lower(i,0) > 0)
        outp << "1   ";
      else if (betaqu_l2_upper(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";

      outp << endl;


      optionsp->out(ST::outresults(nsp,vnames[i],m,
                        stddouble,betaqu_l1_lower(i,0),
                        betaqu50(i,0),betaqu_l1_upper(i,0)) + "\n");

      }

    }

  optionsp->out("\n");

  }



void FC::outresults(ofstream & out_stata, ofstream & out_R,
                    const ST::string & pathresults)
  {

  unsigned nrpar=beta.rows()*beta.cols();

  if (title != "")
    {
    optionsp->out("\n");
    optionsp->out("  " + title + "\n",true);
    optionsp->out("\n");
    }

  if (optionsp->samplesize > 0)
    {

    unsigned i;

    statmatrix<int> index(sampled_beta.rows(),1);

    double* wqu1l = betaqu_l1_lower.getV();
    double* wqu2l = betaqu_l2_lower.getV();
    double* wqu50 = betaqu50.getV();
    double* wqu1u = betaqu_l2_upper.getV();
    double* wqu2u = betaqu_l1_upper.getV();

    if (nosamplessave==false)
      {

      for(i=0;i<nrpar;i++,wqu1l++,wqu2l++,wqu50++,wqu1u++,wqu2u++)
        {

        index.indexinit();

        sampled_beta.indexsort(index,0,sampled_beta.rows()-1,i,0);

        *wqu1l = sampled_beta.quantile(optionsp->lower1,i,index);
        *wqu2l = sampled_beta.quantile(optionsp->lower2,i,index);
        *wqu50 = sampled_beta.quantile(50,i,index);
        *wqu1u = sampled_beta.quantile(optionsp->upper1,i,index);
        *wqu2u = sampled_beta.quantile(optionsp->upper2,i,index);
        }
      }

    if (pathresults.isvalidfile() != 1)
      {

      ofstream outres(pathresults.strtochar());

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
      outres << "pmean   ";

      if (nosamplessave==false)
        {

        if (optionsp->samplesize > 1)
          {
          outres << "pstd   ";
          outres << "pqu"  << l1  << "   ";
          outres << "pqu"  << l2  << "   ";
          outres << "pmed   ";
          outres << "pqu"  << u1  << "   ";
          outres << "pqu"  << u2  << "   ";
          }

        outres << endl;

        double * workmean;
        double * workstd;
        double * workbetaqu_l1_lower_p;
        double * workbetaqu_l2_lower_p;
        double * workbetaqu_l1_upper_p;
        double * workbetaqu_l2_upper_p;
        double * workbetaqu50;

        workmean = betamean.getV();
        workstd =  betavar.getV();
        workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
        workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
        workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
        workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
        workbetaqu50 = betaqu50.getV();

        unsigned nrpar = beta.rows();
        for(i=0;i<nrpar;i++,workmean++,workstd++,workbetaqu_l1_lower_p++,
                                   workbetaqu_l2_lower_p++,workbetaqu50++,
                                   workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
          {
          outres << (i+1) << "   ";
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

            }

          outres << endl;

          }

        }
      else
        {
        if (optionsp->samplesize > 1)
          outres << "pstd   ";

        outres << endl;

        double * workmean;
        double * workstd;

        workmean = betamean.getV();
        workstd =  betavar.getV();

        unsigned nrpar = beta.rows();
        for(i=0;i<nrpar;i++,workmean++,workstd++)
          {
          outres << (i+1) << "   ";
          outres << *workmean << "   ";

          if (optionsp->samplesize > 1)
            {
            if (*workstd < 0.0000000000001)
              outres << 0 << "   ";
            else
              outres << sqrt(*workstd) << "   ";
            }

          outres << endl;

          }

        }

      }

    }  // if (optionsp->get_samplesize() > 0)
  }


void FC::outresults_singleparam(ofstream & out_stata,ofstream & out_R,
                                  const ST::string & pathresults)
  {

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

  ST::string vstr;

  vstr = "    Mean:         ";
  optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
  ST::doubletostring(betamean(0,0),6) + "\n");

  if (optionsp->samplesize > 1)
    {

    vstr = "    Std. dev.:    ";

    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(sqrt(betavar(0,0)),6) + "\n");

    vstr = "    " + l1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_lower(0,0),6) + "\n");

    vstr = "    " + l2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_lower(0,0),6) + "\n");

    vstr = "    50% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu50(0,0),6) + "\n");

    vstr = "    " + u1 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l2_upper(0,0),6) + "\n");

    vstr = "    " + u2 + "% Quantile: ";
    optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
    ST::doubletostring(betaqu_l1_upper(0,0),6) + "\n");

    optionsp->out("\n");
    }


  if (pathresults.isvalidfile() != 1)
    {

    ofstream ou(pathresults.strtochar());

    if (optionsp->samplesize > 1)
      {
      ou << "pmean  pstddev  pqu"  << nl1 << "   pqu" << nl2 << "  pmed pqu" <<
      nu1 << "   pqu" << nu2 << "  pmin  pmax" << endl;
      }
    else
      {
      ou << "pmean" << endl;
      }

    ou << betamean(0,0) << "  ";
    if (optionsp->samplesize > 1)
      {
      ou << (betavar(0,0)<0.0?0.0:sqrt(betavar(0,0))) << "  ";
      ou << betaqu_l1_lower(0,0) << "  ";
      ou << betaqu_l2_lower(0,0) << "  ";
      ou << betaqu50(0,0) << "  ";
      ou << betaqu_l2_upper(0,0) << "  ";
      ou << betaqu_l1_upper(0,0) << "  ";
      ou << betamin(0,0) << "  ";
      ou << betamax(0,0) << "  " << endl;
      }

    optionsp->out("\n");
    }

  }



void FC::reset(void)
  {

  setbeta(beta.rows(),beta.cols(),0);
  acceptance = 0;
  nrtrials = 0;
  }


} // end: namespace MCMC



