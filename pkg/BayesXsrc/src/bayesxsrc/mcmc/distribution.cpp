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



#include "distribution.h"

namespace MCMC
{

using randnumbers::rand_inv_gaussian;

// BEGIN: DSB //
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
// initialize mschecking with a given path for the output
void
DISTRIBUTION::initialise_mscheck(const ST::string & path,
      const statmatrix<int> & index, const unsigned & nrpar,
      const vector<unsigned> & posbeg, const vector<unsigned> & posend)
{
    // only initialize if we have not yet a full predchange container
    if (predchange.rows() < nrobs)
    {
        // activate mscheck
        mscheck = true;

        msnrind = nrpar;
        msindex = index;
        msposbeg = posbeg;
        msposend = posend;

        // initialize predictor change vector into which the
        // random effects objects put their individual changes
        predchange = datamatrix(nrobs, 1, 0);

        // initialize the full conditional object for the mscheck predictor
        ST::string path1 = path + "_mscheck_pred.raw";
        msc_pred = FULLCOND(optionsp, datamatrix(1, 1), "mscheck_pred",
                            nrobs, 1, path1);
        msc_pred.setflags(MCMC::norelchange | MCMC::nooutput);

        // initialize the full conditional object for the mscheck likelihood contributions
        ST::string path2 = path + "_mscheck_like.raw";
        msc_like = FULLCOND(optionsp, datamatrix(1, 1), "mscheck_like",
                            msnrind, 1, path2);
        msc_like.setflags(MCMC::norelchange | MCMC::nooutput);

        // initialize the full conditional object for the mscheck predictive observation samples
        ST::string path3 = path + "_mscheck_obssamples.raw";
        msc_obssamples = FULLCOND(optionsp, datamatrix(1, 1), "mscheck_obssamples",
                            nrobs, 1, path3);
        msc_obssamples.setflags(MCMC::norelchange | MCMC::nooutput);
    }
}

// add the double "m" to the position "row" of "predchange"
void
DISTRIBUTION::add_linearpred_mscheck(const double m, const unsigned int row)
{
    // only if mscheck is activated
    if (mscheck)
    {
        // then go to the position "row" of "predchange" and add the double "m"
        double * workl = predchange.getV() + row;
        * workl += m;
    }
}

// add the double "m" to the range "index"["beg":"end"] of "predchange".
// This is an analogue to the function
// DISTRIBUTION::add_linearpred(const double & m,const unsigned & beg,
//                              const unsigned & end,const statmatrix<int> & index,
//                              const unsigned & col, const bool & current)
void
DISTRIBUTION::add_linearpred_mscheck(const double m,
                                     const unsigned int beg,
                                     const unsigned int end,
                                     const statmatrix<int> & index)
{
    // only if mscheck is activated
    if (mscheck)
    {
        // then get the start position given by the "beg" element of "index"
        int * workindex = index.getV() + beg;

        // and add the "m" to "predchange" until the "end" of the "index" range
        for (unsigned register i = beg; i <= end; i++, workindex++)
        {
            predchange(* workindex, 0) += m;
        }
    }
}

void DISTRIBUTION::get_mssamples()
  {
  ST::string pathhelp =
        pathresultsscale.substr(0, pathresultsscale.length() - 10)
                        + "_mscheck_like_sample.raw";
  optionsp->out(pathhelp + "\n\n");
  msc_like.get_samples(pathhelp);
  pathhelp = pathresultsscale.substr(0, pathresultsscale.length()
                - 10) + "_mscheck_pred_sample.raw";
  optionsp->out(pathhelp + "\n\n");
  msc_pred.get_samples(pathhelp);
  pathhelp = pathresultsscale.substr(0, pathresultsscale.length() - 10)
                       + "_mscheck_obssamples_sample.raw";
  optionsp->out(pathhelp + "\n\n");
  msc_obssamples.get_samples(pathhelp);
  }
#endif
// END: DSB


void DISTRIBUTION::compute_overall_deviance(double & deviance,double & deviancesat)
  {

  unsigned i;
  double * workresp=response.getV();
  double * worklin = (*linpred_current).getV();
  double * workweight = weight.getV();
  double dev=0;
  double devsat=0;
  double mu;            // Fehler bei multinomial!

  for(i=0;i<nrobs;i++,workresp++,worklin++,workweight++)
    {
    if (*workweight != 0)
      {
      compute_mu(worklin,&mu);
      compute_deviance(workresp,workweight,&mu,&dev,&devsat,scale,0);
      deviance+=dev;
      deviancesat+=devsat;
      }
    }
  }

double DISTRIBUTION::compute_msep(void)
  {
  unsigned i;
  unsigned nlcols = linearpred.cols();

  double * worklin = (*linpred_current).getV();
  double * workresp = response.getV();
  double * workweight = weight.getV();
  double dev=0;
  double devsat=0;
  double deviancesat = 0;
  double deviance = 0;

  datamatrix mu = datamatrix(nlcols,1,0);
  double * mum = mu.getV();
  // double mu;   // siehe "compute_deviance"

  double nr = 0;
  double * work2 = weight2.getV();

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++,work2++)
    {
    if (*workweight == 0)
      {
      compute_mu(worklin,mum);
      compute_deviance(workresp,work2,mum,&dev,&devsat,scale,0);
      deviance+=dev;
      deviancesat+=devsat;
      worklin += nlcols - 1;
      workresp += nlcols - 1;

      nr += 1;
      }
    else       // überprüfen!!!!!
      {
      worklin += nlcols - 1;
      workresp += nlcols - 1;
      }
    }

  return  deviancesat;   //      /nr;
  }


double DISTRIBUTION::compute_rss(void)
  {
  unsigned i;

  double * worklin = (*linpred_current).getV();
  double * workresp = tildey.getV();
  double * workweight = weightiwls.getV();

  double sum = 0;
  double help;

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    if (*workweight !=0)
      {
      help = *workresp - *worklin;
      sum += *workweight*help*help;
      }
    }
  sum = trmult(0,0)*trmult(0,0)*sum;
  return  sum;
  }


double DISTRIBUTION::compute_gcv(const double & df)
  {
  double dev=0;
  double devsat=0;
  compute_overall_deviance(dev,devsat);

  double help2 = 1- gcvfactor * df/get_nrobs_wpw();
  double gcv = devsat / (get_nrobs_wpw()*help2*help2);
  return gcv;
  }

double DISTRIBUTION::compute_gcv2(const double & df)
  {
  //double dev=0;
  //double devsat=0;
  //compute_overall_deviance(dev,devsat);

  double help2 = 1 - gcvfactor * df/get_nrobs_wpw();
  //double gcv = dev  - 2*get_nrobs_wpw() * log(help2);
  double gcv = compute_rss() / (get_nrobs_wpw()*help2*help2);
  return gcv;
  }

double DISTRIBUTION::compute_aic(const double & df)
  {
  double dev=0;
  double devsat=0;
  compute_overall_deviance(dev,devsat);

  double aic = dev + 2*df;
  return aic;
  }

double DISTRIBUTION::compute_improvedaic(const double & df)
  {
  double dev=0;
  double devsat=0;
  compute_overall_deviance(dev,devsat);

  double impaic = dev + 2*df + 2*df*(df+1)/(get_nrobs_wpw()-df-1);
  return impaic;
  }

double DISTRIBUTION::compute_bic(const double & df)
  {
  double dev=0;
  double devsat=0;
  compute_overall_deviance(dev,devsat);

  double bic = dev + log(static_cast<double>(get_nrobs_wpw()))*df;
  return bic;
  }


void DISTRIBUTION::set_nrobs_wpw(void)
  {
  unsigned i;
  double * workweight = weight.getV();
  double s = 0;

  for (i=0;i<nrobs;i++,workweight++)
    {
    if (*workweight == 0)
      s++;
    }

  nrobs_wpw = nrobs-s;
  }

unsigned DISTRIBUTION::get_nrpar(void)
  {
  unsigned nrpar=0;
  if (scaleexisting)
    nrpar += scale.rows()*scale.cols();

  return nrpar;
  }


void DISTRIBUTION::create_bootstrap_weights(void)   // wie bei MSEP und CV???
  {

  if(isbootstrap==false)
    isbootstrap = true;
  if(response_ori.rows()==1)
    {
    srand(seed);
    response_ori = response;  // speichert Response vom Original-Datensatz
    linearpred_ori = linearpred;     // speichert Prädiktor vom Original-Datensatz
    }

  double * workresp = response.getV();
  double * worklin = linearpred_ori.getV();
  double * workw = weight.getV();

  unsigned i;

  for(i=0;i<nrobs;i++,workresp+=linearpred.cols(),worklin+=linearpred.cols(),workw++)
    {
    compute_bootstrap_data(worklin,workw,workresp);
    }
  }


void DISTRIBUTION::set_original_response(void)
  {
  response = response_ori;
  }


void DISTRIBUTION::update_bootstrap_betamean(void)
  {
  linearpred = linearpred_ori;
  }

void DISTRIBUTION::save_betamean(void)
  {
  srand(seed);
  response_ori = response;  // speichert Response vom Original-Datensatz
  linearpred_ori = linearpred;     // speichert Residuen vom Original-Datensatz
  }


void DISTRIBUTION::create_weight(datamatrix & w, const double & p1, const bool & fertig, const bool & CV)
  {
  weight2 = weight;          // speichert die alten Gewichte ab

  srand(seed);

  if(fertig == true)        // erstellt neue Gewichte für MSEP / AUC
    {
    constant_iwlsweights = false;
    iwlsweights_notchanged_df = false;
    double a = 0;
    unsigned i;
    double * workweight = w.getV();
    for(i=0;i<nrobs;i++,workweight++)
      a += *workweight;
    double p = (p1 * nrobs - a) / (nrobs - a);
    double u;
    workweight = w.getV();
    for(i=0;i<nrobs;i++,workweight++)
      {
      u = uniform();
      if(u<p && *workweight==0)
        *workweight = 1;
      }

    double * workw = w.getV();
    workweight = weight.getV();
    for(i=0;i<nrobs;i++,workweight++,workw++)
      {
      if(*workw == 0)
        *workweight = 0;
      }
    }

  if(CV == true)
    {
    unsigned i;
    double anzahl = floor(double(nrobs)/p1);    // p1 hier =5 oder =10, gibt Aufteilung an
    double rest = fmod(double(nrobs),p1);      // "anzahl" = Anzahl der Beobachtungen pro Teil,
                                              // "rest" = Anzahl Teile mit einer Beobachtung mehr
    constant_iwlsweights = false;
    iwlsweights_notchanged_df = false;
    weightcv = datamatrix(nrobs,p1,1);

    if(weight.max(0) != 1 || weight.max(0) != 1)  // falls Gewichte im Befehl angegeben wurden, werden diese berücksuchtigt
      {
      double * workw = weight.getV();
      double * workcv = weightcv.getV();
      for(i=0;i<nrobs;i++,workw++)
        {
        for(unsigned j=0;j<p1;j++,workcv++)
          *workcv = *workw;
        }
      }

    datamatrix u = datamatrix(nrobs,1,0);
    for(i=0;i<nrobs;i++)
      u(i,0) = uniform();

    statmatrix<int> index(nrobs,1);
    index.indexinit();
    u.indexsort(index,0,nrobs-1,0,0);

    unsigned pos = 0;
    unsigned anz=0;
    int* p = index.getV();
    for(i=0;i<index.rows();i++,p++)
      {
      double plus = 0;
      if(rest > 0)
        plus += 1;
      if(anz < anzahl + plus)
        {
        weightcv(*p,pos) = 0;
        anz = anz+1;
        }
      else
        {
        pos += 1;
        rest -= 1;
        weightcv(*p,pos) = 0;
        anz = 1;
        }
      }
    }
  }


void DISTRIBUTION::weight_for_all(void)
  {
  constant_iwlsweights = false;
  iwlsweights_notchanged_df = false;
  datamatrix save = datamatrix(weight.rows(),1,0);
  save.assign(weight);
  weight = weight2;
  weight2 = save;
  if(weightiwls2.rows()>1)
    {
    save.assign(weightiwls);
    weightiwls = weightiwls2;
    weightiwls2 = save;
    }
  }


void DISTRIBUTION::compute_cvweights(int pos)
  {
  if(pos < 0)
    {
    weight = weight2;
    weightiwls = weightiwls2;
    }
  if(pos >= 0)
    {
    double * w1 = weight.getV();
    double * w2 = weightiwls.getV();
    double * mw = weightcv.getV()+pos;
    unsigned cv = weightcv.cols();
    unsigned i;
    for(i=0;i<nrobs;i++,w1++,w2++,mw+=cv)
      {
      *w1 *= *mw;
      *w2 *= *mw;
      }
    }
  }


void DISTRIBUTION::save_weightiwls(void)
  {
  weightiwls2 = weightiwls;
  }


void DISTRIBUTION::create(MCMCoptions * o, const datamatrix & r,
                           const datamatrix & w)
  {

  optionsp = o;

  scaleexisting = true;
  scale=datamatrix(1,1,0.1);
  scale_mode=scale;

  family = "unknown";

  response = r;
  nrobs = response.rows();
  trmult = datamatrix(r.cols(),1,1.0);

  helpmat1 = datamatrix(nrobs,1);

  if (w.rows() == 1)
    {
    weight = datamatrix(r.rows(),1,1);
    nrobsmweightzero=nrobs;
    }
  else
    {
    assert(w.rows()==r.rows());
    weight = w;
    unsigned i;
    unsigned nz=0;
    double * workweight=weight.getV();
    bool error = false;
    for(i=0;i<nrobs;i++,workweight++)
      {
      if (*workweight == 0)
        nz++;
      else if (*workweight < 0)
        error = true;
      }
    nrobsmweightzero = nrobs-nz;
    if (error)
      errors.push_back("ERROR: negative weights encountered\n");
    }

  linearpred = datamatrix(nrobs,r.cols(),0);
  linearpredprop = linearpred;
  linpred_current = &linearpred;
  linpred_proposed = &linearpredprop;

  sumweight = weight.sum(0);

  tildey = datamatrix(nrobs,r.cols());
  weightiwls = datamatrix(nrobs,r.cols(),1);
  workingres = datamatrix(nrobs,1);

  predict = false;
  predictfull = false;
  linpredmean = datamatrix(1,1);
  mumean = datamatrix(1,1);
  deviancemean = datamatrix(1,1);

  predictindicator = datamatrix(1,1);
  predictresponse = false;

  responsename = "Y";

  weightname = "W";
  changingweight = false;

  interceptsample = datamatrix(optionsp->compute_samplesize(),r.cols(),0);

  addinterceptsample=0;

  constant_iwlsweights=false;
  iwlsweights_notchanged_df = false;
  nrobs_wpw = -1;
  isbootstrap = false;
  }


DISTRIBUTION::DISTRIBUTION(MCMCoptions * o, const datamatrix & r,
                           const datamatrix & w,const ST::string & pr,
                           const ST::string & ps)
  {
  nosamples = false;
// Begin: DSB //
  #if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
  mscheck=false;
  #endif
// End: DSB //
  pathresultsscale = pr;
  Scalesave = FULLCOND(o,datamatrix(1,1),"Scaleparameter",1,1,ps);
  Scalesave.setflags(MCMC::norelchange | MCMC::nooutput);

  shrinkage=false;
  nrridge=0;      //new
  ridgesum=0.0;   //new
  nrlasso=0;      //new
  lassosum=0.0;   //new
  nrnigmix=0;     //new
  nigmixsum=0.0;  //new
  create(o,r,w);
  }


DISTRIBUTION::DISTRIBUTION(const datamatrix & offset,MCMCoptions * o,
                           const datamatrix & r,const datamatrix & w,
                           const ST::string & pr,const ST::string & ps)
  {
  nosamples = false;
// Begin: DSB
  #if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
  mscheck=false;
  #endif
// End : DSB
  pathresultsscale = pr;
  Scalesave = FULLCOND(o,datamatrix(1,1),"Scaleparameter",1,1,ps);
  Scalesave.setflags(MCMC::norelchange | MCMC::nooutput);

  shrinkage=false;
  nrridge=0;      //new
  ridgesum=0.0;   //new
  nrlasso=0;      //new
  lassosum=0.0;   //new
  nrnigmix=0;     //new
  nigmixsum=0.0;  //new

  create(o,r,w);
  add_linearpred(offset);

  }


void DISTRIBUTION::set_predict(const ST::string & path,
                               const ST::string & pathdev,
                               datamatrix * p, vector<ST::string> & dn)
  {
  Dnames = dn;
  Dp = p;
  predict = true;
  predictpath = path;
  deviancepath = pathdev;
  linpredmean = datamatrix(nrobs,response.cols(),0);
  mumean = datamatrix(nrobs,response.cols(),0);
  deviancemean = datamatrix(nrobs,1,0);
  deviancemean_sat = datamatrix(nrobs,1,0);
  deviance = datamatrix(optionsp->compute_samplesize(),2,0);
  }

void DISTRIBUTION::set_predictfull(const ST::string & pathsample,
                                   const ST::string & path,const unsigned & fo)
  {
  predictfullpath = path;
  predictfull = true;
  firstobs=fo;
  musave = FULLCOND(optionsp,datamatrix(1,1),"Predictmu",
           fo,mumean.cols(),pathsample);
  musave.setflags(MCMC::norelchange | MCMC::nooutput);

  }

void DISTRIBUTION::set_predictresponse(const datamatrix & pr)
  {
  predictindicator = pr;
  predictresponse = true;
  }




void DISTRIBUTION::init_names(const ST::string & rn, const ST::string & wn,const ST::string & on)
  {
  responsename = rn;
  weightname = wn;
  offsetname = on;
  }


void DISTRIBUTION::init_offset(const datamatrix & o)
  {
  add_linearpred(o);
  }


ST::string DISTRIBUTION::get_scale_sample(void) const
  {

  ST::string file;

  if (scaleexisting)
    {

    file = pathresultsscale.substr(0,pathresultsscale.length()-4) + "_sample.raw";
    Scalesave.get_samples(file);

    }

/*
    file = pathresultsscale.substr(0,pathresultsscale.length()-9) + "mean_sample.raw";
    musave.get_samples(file);

    file = pathresultsscale.substr(0,pathresultsscale.length()-9) + "_predictor_sample.raw";
    responsesave.get_samples(file);*/

  return file;

  }

ST::string DISTRIBUTION::get_mean_sample(void) const
  {

  ST::string file;

  if(predictfull)
    {

    file = predictfullpath.substr(0,predictfullpath.length()-4) + "_mean_sample.raw";
    musave.get_samples(file);

/*  file = pathresultsscale.substr(0,pathresultsscale.length()-9) + "_predictor_sample.raw";
    responsesave.get_samples(file);*/

    }

  return file;

  }


datamatrix DISTRIBUTION::compute_autocor_scale(
                                      const unsigned & lag,const unsigned & row,
                                      const unsigned & col) const
  {
  assert(optionsp->get_samplesize() > 0);
  datamatrix sample;
  sample = Scalesave.compute_autocorr(lag,row,col);
  return sample;
  }


DISTRIBUTION::DISTRIBUTION(const DISTRIBUTION & d)
  {

  constant_iwlsweights=d.constant_iwlsweights;
  iwlsweights_notchanged_df = d.iwlsweights_notchanged_df;

  nosamples = d.nosamples;

  optionsp = d.optionsp;

  scaleexisting = d.scaleexisting;

  scale = d.scale;
  scale_mode = d.scale_mode;
  pathresultsscale = d.pathresultsscale;
  Scalesave = d.Scalesave;
  acceptancescale = d.acceptancescale;

  nrobs = d.nrobs;
  nrobsmweightzero = d.nrobsmweightzero;

  response = d.response;
  tildey = d.tildey;
  weightiwls = d.weightiwls;
  workingres = d.workingres;
  responsename = d.responsename;
  trmult = d.trmult;

  weight = d.weight;
  weightname = d.weightname;
  weight2 = d.weight2;
  weightiwls2 = d.weightiwls2;
  nrobs_wpw = d.nrobs_wpw;

  isbootstrap = d.isbootstrap;
  linearpred_ori = d.linearpred_ori;
  response_ori = d.response_ori;
  gcvfactor = d.gcvfactor;
  seed = d.seed;

  offsetname = d.offsetname;

  linearpred = d.linearpred;
  linearpredprop = linearpred;
  linpred_current = &linearpred;
  linpred_proposed = &linearpredprop;

  predict = d.predict;
  predictfull = d.predictfull;
  predictpath = d.predictpath;
  predictfullpath = d.predictfullpath;
  deviancepath = d.deviancepath;
  deviance = d.deviance;
  linpredmean = d.linpredmean;

  results_latex = d.results_latex;

  predictindicator=d.predictindicator;
  predictresponse = d.predictresponse;

  mumean = d.mumean;
  musave = d.musave;
  firstobs = d.firstobs;
  deviancemean = d.deviancemean;
  deviancemean_sat = d.deviancemean_sat;
  Dp = d.Dp;
  Dnames = d.Dnames;

  sumweight = d.sumweight;

  changingweight = d.changingweight;

  errors = d.errors;

  family = d.family;

  helpmat1 = d.helpmat1;

  fcmissing = d.fcmissing;
  missingpos = d.missingpos;
  MissingSave = d.MissingSave;
  missingind = d.missingind;
  pathmissing = d.pathmissing;

  interceptsample = d.interceptsample;
  interceptold = d.interceptold;
  addinterceptsample = d.addinterceptsample;

  shrinkage = d.shrinkage;
  nrridge = d.nrridge;
  ridgesum = d.ridgesum;
  nrlasso = d.nrlasso;
  lassosum = d.lassosum;
  nrnigmix = d.nrnigmix;
  nigmixsum = d.nigmixsum;

// Begin: DSB
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
  mscheck = d.mscheck;
  msc_pred = d.msc_pred;
  msc_like = d.msc_like;
  msc_obssamples = d.msc_obssamples;
  predchange = d.predchange;
  msindex = d.msindex;
  msnrind = d.msnrind;
  msposbeg = d.msposbeg;
  msposend = d.msposend;
#endif
// End: DSB

  }


const DISTRIBUTION & DISTRIBUTION::operator=(const DISTRIBUTION & d)
  {
  if (this == &d)
    return *this;

  constant_iwlsweights=d.constant_iwlsweights;
  iwlsweights_notchanged_df = d.iwlsweights_notchanged_df;

  nosamples = d.nosamples;

  optionsp = d.optionsp;

  scaleexisting = d.scaleexisting;
  scale = d.scale;
  scale_mode = d.scale_mode;
  Scalesave = d.Scalesave;
  pathresultsscale = d.pathresultsscale;
  acceptancescale = d.acceptancescale;

  nrobs = d.nrobs;
  nrobsmweightzero = d.nrobsmweightzero;

  response = d.response;
  tildey = d.tildey;
  weightiwls = d.weightiwls;
  workingres = d.workingres;
  responsename = d.responsename;
  trmult = d.trmult;

  weight = d.weight;
  weightname = d.weightname;
  weight2 = d.weight2;
  weightiwls2 = d.weightiwls2;

  nrobs_wpw = d.nrobs_wpw;
  isbootstrap = d.isbootstrap;
  linearpred_ori = d.linearpred_ori;
  response_ori = d.response_ori;
  gcvfactor = d.gcvfactor;
  seed = d.seed;

  offsetname = d.offsetname;

  linearpred = d.linearpred;
  linearpredprop = linearpred;
  linpred_current = &linearpred;
  linpred_proposed = &linearpredprop;

  results_latex = d.results_latex;

  linpredmean = d.linpredmean;
  predict = d.predict;
  predictfull = d.predictfull;
  predictpath = d.predictpath;
  predictfullpath = d.predictfullpath;
  deviancepath = d.deviancepath;
  deviance = d.deviance;
  mumean = d.mumean;
  deviancemean = d.deviancemean;
  deviancemean_sat = d.deviancemean_sat;
  musave = d.musave;
  firstobs = d.firstobs;
  Dp = d.Dp;
  Dnames = d.Dnames;

  predictindicator=d.predictindicator;
  predictresponse = d.predictresponse;

  sumweight = d.sumweight;
  changingweight = d.changingweight;

  errors = d.errors;

  family = d.family;

  helpmat1 = d.helpmat1;

  fcmissing = d.fcmissing;
  missingpos = d.missingpos;
  MissingSave = d.MissingSave;
  missingind = d.missingind;
  pathmissing = d.pathmissing;

  interceptsample = d.interceptsample;
  interceptold = d.interceptold;
  addinterceptsample = d.addinterceptsample;

  shrinkage = d.shrinkage;
  nrridge = d.nrridge;
  ridgesum = d.ridgesum;
  nrlasso = d.nrlasso;
  lassosum = d.lassosum;
  nrnigmix = d.nrnigmix;
  nigmixsum = d.nigmixsum;

// Begin: DSB
  #if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
  mscheck = d.mscheck;
  msc_pred = d.msc_pred;
  msc_like = d.msc_like;
  msc_obssamples = d.msc_obssamples;
  predchange = d.predchange;
  msindex = d.msindex;
  msnrind = d.msnrind;
  msposbeg = d.msposbeg;
  msposend = d.msposend;
  #endif
// End: DSB
  return *this;
  }


void DISTRIBUTION::outerrors(void)
  {
  unsigned i;
  for (i=0;i<errors.size();i++)
    optionsp->outerror(errors[i]);
  }


void DISTRIBUTION::outoptions(void)
  {
  optionsp->out("RESPONSE DISTRIBUTION:\n",true);
  optionsp->out("\n");
  optionsp->out("  Family: " + family + "\n");
  optionsp->out("  Number of observations: " + ST::inttostring(nrobs) + "\n");
  optionsp->out("  Number of observations with positive weights: "
                  + ST::inttostring(nrobsmweightzero) + "\n");
  }


void DISTRIBUTION::set_interceptsample(datamatrix & s,unsigned & column)
  {
  unsigned i;
  double * work = interceptsample.getV()+column;
  double * works = s.getV();
  unsigned size = interceptsample.cols();
  for (i=0;i<interceptsample.rows();i++,work+=size,works++)
    {
    *work = *works;
    }
  }


double DISTRIBUTION::compute_IWLS(datamatrix & weightiwls,datamatrix & tildey,
                                  bool weightyes,
                                  const unsigned & col,const bool & current)
  {

  register unsigned  i;
  double* workweight = weight.getV();
  double* workres = response.getV();
  double help = 0;

  double* worklin;
  if (current)
    worklin = (*linpred_current).getV();
  else
    worklin = (*linpred_proposed).getV();

  double * worktildey=tildey.getV();

  double * worklinhelp;
  double * workresphelp;

  double * workweightiwls = weightiwls.getV();

  if (weightyes)
    {

    for (i=0;i<nrobs;i++,workweight++,worktildey++,workweightiwls++,
         worklin+=linearpred.cols(),workres+=linearpred.cols())
      {
      worklinhelp = worklin;
      workresphelp = workres;
      help += compute_IWLS(workresphelp,worklinhelp,workweight,i,workweightiwls,
                           worktildey,true,col);
      }
    }
  else
    {
    for (i=0;i<nrobs;i++,workweight++,worktildey++,
         worklin+=linearpred.cols(),workres+=linearpred.cols())
      {
      worklinhelp = worklin;
      workresphelp = workres;
      help += compute_IWLS(workresphelp,worklinhelp,workweight,i,workweightiwls,
                           worktildey,false,col);
      }
    }

  return help;

  }


void DISTRIBUTION::set_response(datamatrix & newresponse)
  {
  response.assign(newresponse);
  }


void DISTRIBUTION::compute_IWLS_weight_tildey(datamatrix & weightiwls,
                                      datamatrix & tildey,
                                      const unsigned & col,
                                      const bool & current)
  {

  register unsigned  i;
  double* workweight = weight.getV();
  double* workres = response.getV();

  double* worklin;
  if (current)
    worklin = (*linpred_current).getV();
  else
    worklin = (*linpred_proposed).getV();

  double * worktildey=tildey.getV();
  double * workweightiwls = weightiwls.getV();

  double * worklinhelp;
  double * workresphelp;

  for (i=0;i<nrobs;i++,workweight++,worktildey++,workweightiwls++,
       worklin+=linearpred.cols(),workres+=linearpred.cols())
    {
    worklinhelp = worklin;
    workresphelp = workres;
    compute_IWLS_weight_tildey(workresphelp,worklinhelp,workweight,i,
                            workweightiwls,worktildey,col);
    }

  }



double DISTRIBUTION::loglikelihood(const bool & current) const
  {

  register unsigned  i;
  double* workweight = weight.getV();
  double* workres = response.getV();
  double help = 0;

  double* worklin;
  if (current)
    worklin = (*linpred_current).getV();
  else
    worklin = (*linpred_proposed).getV();

  double * worklinhelp;
  double * workresphelp;
  for (i=0;i<nrobs;i++,workweight++,worklin+=linearpred.cols(),
       workres+=linearpred.cols())
    {
    worklinhelp = worklin;
    workresphelp = workres;
    help += loglikelihood(workresphelp,worklinhelp,workweight,i);
    }

  return help;

  }


double DISTRIBUTION::loglikelihood(const unsigned & beg,const unsigned & end,
                                   const statmatrix<int> & index,
                                   const bool & current)
  {

  unsigned i;
  double help = 0;

  int* workind = index.getV()+beg;

  if (current)
    {
    for (i=beg;i<=end;i++,workind++)
      help+=loglikelihood(&response(*workind,0),&((*linpred_current)(*workind,0)),
                          &weight(*workind,0),*workind);

    }
  else
    {
    for (i=beg;i<=end;i++,workind++)
      help+=loglikelihood(&response(*workind,0),&((*linpred_proposed)(*workind,0)),
                          &weight(*workind,0),*workind);
    }

  return help;

  }


double DISTRIBUTION::loglikelihood2(const unsigned & beg,const unsigned & end,
                                    const statmatrix<int> & index,
                                    const statmatrix<int> & index2,
                                    const bool & current)
  {

  unsigned i;
  double help = 0;

  int* workind2 = index2.getV()+beg;

  double * workresponse = &response(index(beg,0),0);
  double * workresponsehelp;
  double * workweight = &weight(index(beg,0),0);
  int obsnr = index(beg,0);

  double * worklinpred;

  if (current)
    worklinpred = &((*linpred_current)(index(beg,0),0));
  else
    worklinpred = &((*linpred_proposed)(index(beg,0),0));

  for (i=beg;i<=end;i++,workind2++,
       workresponse += *workind2*response.cols(),workweight += *workind2,
       obsnr += *workind2,worklinpred+=*workind2*linearpred.cols())
    {
    workresponsehelp = workresponse;
    if (*workweight != 0)
      help+=loglikelihood(workresponsehelp,worklinpred,workweight,obsnr);

    }

  return help;

  }

// BEGIN: DSB //
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
// FUNCTION: loglikelihood_from_deviance
// TASK: computes the loglikelihood for a single observation for univariate
// response, by using the compute_deviance function.

    double
    DISTRIBUTION::loglikelihood_from_deviance(const double res, // response
                                              const double mu, // mean sample
                                              const double weight // weight
                                              ) const
    {
        // compute the deviance
        double deviance = 0;
        double unused = 0;

        compute_deviance(&res, &weight, &mu, &deviance, &unused, scale, 0);

        // and return the corresponding loglikelihood
        return - 0.5 * deviance;
    }

#endif
// END: DSB //

void DISTRIBUTION::compute_weight(datamatrix & w,const unsigned & col,
                                  const bool & current) const
  {

  register unsigned i;

  double* workw = w.getV();

  double * worklin;
  if (current)
    worklin = (*linpred_current).getV();
  else
    worklin = (*linpred_proposed).getV();

  double * workweight = weight.getV();

  double * worklinhelp;

  unsigned size = linearpred.cols();
  for (i=0;i<nrobs;i++,workw++,worklin+=size,workweight++)
    {
    worklinhelp = worklin;
    *workw = compute_weight(worklinhelp,workweight,i,col);
    }

  }

double DISTRIBUTION::compute_sumweight(const unsigned & col,
                                       const bool & current) const
  {
  double sumw = 0;

  register unsigned i;

  double * worklin;
  if (current)
    worklin = (*linpred_current).getV();
  else
    worklin = (*linpred_proposed).getV();

  double * workweight = weight.getV();

  double * worklinhelp;

  unsigned size = linearpred.cols();
  for (i=0;i<nrobs;i++,worklin+=size,workweight++)
    {
    if (*workweight != 0)
      {
      worklinhelp = worklin;
      sumw +=compute_weight(worklinhelp,workweight,i,col);
      }
    }

  return sumw;
  }


void DISTRIBUTION::compute_weight(datamatrix & w,const unsigned & beg,
                                  const unsigned & end,
                       const statmatrix<int> & index, const unsigned & col)
  {

  register unsigned i;
  int * workind = index.getV() + beg;

  for (i=beg;i<=end;i++,workind++)
    {
    w(*workind,0) = compute_weight( &(*linpred_current)(*workind,0),
                                    &weight(*workind,0),*workind,col );

    }

  }


double DISTRIBUTION::compute_sumweight(const unsigned & beg,
                                       const unsigned & end,
                                       const statmatrix<int> & index,
                                       const unsigned & col,
                                       const bool & current)
  {

  register unsigned i;

  int * workind = index.getV() + beg;

  double sum=0;

  if (current==true)
    {
    for (i=beg;i<=end;i++,workind++)
      {
      sum += compute_weight( &(*linpred_current)(*workind,0),
                             &weight(*workind,0),*workind,col );

      }
    }
  else
    {
    for (i=beg;i<=end;i++,workind++)
      {
      sum += compute_weight( &(*linpred_proposed)(*workind,0),
                             &weight(*workind,0),*workind,col );

      }

    }

  return sum;

  }


double DISTRIBUTION::compute_sumweight2(const unsigned & beg,
                                        const unsigned & end,
                                        const statmatrix<int> & index,
                                        const statmatrix<int> & index2,
                                        const unsigned & col,
                                        const bool & current)
  {

  register unsigned i;

  int * workind = index.getV() + beg;
  int * workind2 = index2.getV() + beg;

  double sum=0;

  if (current==true)
    {
    for (i=beg;i<=end;i++,workind2++,workind++)
      {
      linpredp_current += *workind2*linearpred.cols();
      if (weight(*workind,0) != 0)
        sum += compute_weight(linpredp_current,&weight(*workind,0),*workind,col);
      }
    }
  else
    {
    for (i=beg;i<=end;i++,workind2++,workind++)
      {
      linpredp_proposed += *workind2*linearpred.cols();
      if (weight(*workind,0) != 0)
        sum += compute_weight(linpredp_proposed, &weight(*workind,0),
               *workind,col);
      }

    }

  return sum;

  }


double DISTRIBUTION::compute_sumweight_sumy(double beta,
                                            double & sumweight,
                                            const unsigned & beg,
                                            const unsigned & end,
                                            const statmatrix<int> & index,
                                            const statmatrix<int> & index2,
                                            const unsigned & col,
                                            const bool & current)
  {

  register unsigned i;

  int * workind = index.getV() + beg;
  int * workind2 = index2.getV() + beg;

  double * workresponse = &response(index(beg,0),0);
  double * workresponsehelp;
  double * workweight = &weight(index(beg,0),0);
  int obsnr = index(beg,0);
  double * worklinpred;

  sumweight=0;
  double sumy = 0;
  double weightiwls_help;
  double tildey_help;

  if (current==true)
    worklinpred = &((*linpred_current)(index(beg,0),0));
  else
    worklinpred = &((*linpred_proposed)(index(beg,0),0));

  for (i=beg;i<=end;i++,workind2++,workind++,
       workresponse += *workind2*response.cols(),
       workweight += *workind2,obsnr += *workind2,
       worklinpred+=*workind2*linearpred.cols())
    {
    workresponsehelp = workresponse;
    compute_IWLS_weight_tildey(
    workresponsehelp,worklinpred,workweight,obsnr,
    &weightiwls_help,&tildey_help,col);

    sumweight += weightiwls_help;
    sumy += weightiwls_help*(tildey_help+beta);

    }

  return sumy;

  }


double DISTRIBUTION::compute_sumweight_sumy(
double beta,double & sumweight,const unsigned & beg, const unsigned & end,
                                const datamatrix & data,
                                const statmatrix<int> & index,
                                const statmatrix<int> & index2,
                                const unsigned & col,
                                const bool & current)

  {

  register unsigned i;

  int * workind = index.getV() + beg;
  int * workind2 = index2.getV() + beg;

  double * workdata= data.getV()+beg;

  double * workresponse = &response(*workind,0);
  double * workresponsehelp;
  double * workweight = &weight(*workind,0);
  int obsnr = *workind;
  double * worklinpred;

  sumweight=0;
  double sumy = 0;
  double weightiwls_help;
  double tildey_help;

  if (current==true)
    worklinpred = &((*linpred_current)(*workind,0));
  else
    worklinpred = &((*linpred_proposed)(*workind,0));

  for (i=beg;i<=end;i++,workind2++,
       workresponse += *workind2*response.cols(),
       workweight += *workind2,obsnr += *workind2,
       worklinpred+=*workind2*linearpred.cols(),workdata++)
    {
    workresponsehelp = workresponse;
    compute_IWLS_weight_tildey(
    workresponsehelp,worklinpred,workweight,obsnr,
    &weightiwls_help,&tildey_help,col);

    sumweight += weightiwls_help* (*workdata) * (*workdata);

    sumy += (*workdata)*weightiwls_help*(tildey_help + *workdata * beta);

    }

  return sumy;

  }


double DISTRIBUTION::compute_loglikelihood_sumweight_sumy(
                                double beta,double & sumweight,
                                double & sumy,
                                const unsigned & beg, const unsigned & end,
                                const statmatrix<int> & index,
                                const statmatrix<int> & index2,
                                const unsigned & col,
                                const bool & current)
  {

  sumweight=0;
  sumy = 0;
  double loglikel = 0;


  register unsigned i;

  int * workind = index.getV() + beg;
  int * workind2 = index2.getV() + beg;

  double * workresponse = &response(index(beg,0),0);
  double * workresponsehelp;
  double * workweight = &weight(index(beg,0),0);
  int obsnr = index(beg,0);
  double * worklinpred;

  double weightiwls_help;
  double tildey_help;

  if (current==true)
    worklinpred = &((*linpred_current)(index(beg,0),0));
  else
    worklinpred = &((*linpred_proposed)(index(beg,0),0));

  for (i=beg;i<=end;i++,workind2++,workind++,
       workresponse += *workind2*response.cols(),
       workweight += *workind2,obsnr += *workind2,
       worklinpred+=*workind2*linearpred.cols())
    {
    workresponsehelp = workresponse;

    loglikel += compute_IWLS(workresponsehelp,worklinpred,workweight,obsnr,
                             &weightiwls_help,&tildey_help,true,col);

    sumweight += weightiwls_help;
    sumy += weightiwls_help*(tildey_help+beta);

    }

  return loglikel;

  }

double DISTRIBUTION::compute_loglikelihood_sumweight_sumy(
                                double beta,double & sumweight,
                                double & sumy, const unsigned & beg,
                                const unsigned & end,
                                const datamatrix & data,
                                const statmatrix<int> & index,
                                const statmatrix<int> & index2,
                                const unsigned & col,
                                const bool & current)

  {

  sumweight=0;
  sumy = 0;
  double loglikel=0;

  register unsigned i;

  int * workind = index.getV() + beg;
  int * workind2 = index2.getV() + beg;

  double * workdata= data.getV()+beg;

  double * workresponse = &response(*workind,0);
  double * workresponsehelp;
  double * workweight = &weight(*workind,0);
  int obsnr = *workind;
  double * worklinpred;

  double weightiwls_help;
  double tildey_help;

  if (current==true)
    worklinpred = &((*linpred_current)(*workind,0));
  else
    worklinpred = &((*linpred_proposed)(*workind,0));

  for (i=beg;i<=end;i++,workind2++,
       workresponse += *workind2*response.cols(),
       workweight += *workind2,obsnr += *workind2,
       worklinpred+=*workind2*linearpred.cols(),workdata++)
    {
    workresponsehelp = workresponse;

    loglikel += compute_IWLS(workresponsehelp,worklinpred,workweight,obsnr,
                             &weightiwls_help,&tildey_help,true,col);

    sumweight += weightiwls_help* (*workdata) * (*workdata);

    sumy += (*workdata)*weightiwls_help*(tildey_help + *workdata * beta);

    }

  return loglikel;

  }



void DISTRIBUTION::fisher(datamatrix & XWX,datamatrix & data,
                          const unsigned & col) const
  {
  register unsigned i;
  unsigned p,k;

  unsigned nrconst = data.cols();


  double* workw;
  double*workXp;
  double* workXk;

  for (p=0;p<nrconst;p++)
    for (k=p;k<nrconst;k++)
      {
      XWX(p,k)=0;
      workw = weightiwls.getV()+col;
      workXp = data.getV()+p;
      workXk = data.getV()+k;
      for(i=0;i<nrobs;i++,workw+=weightiwls.cols(),
        workXp+=nrconst,workXk+=nrconst)
        XWX(p,k)+= *workw  *  *workXp * *workXk;
      XWX(k,p) = XWX(p,k);
      }

  }


void DISTRIBUTION::fisher(datamatrix & XWX,datamatrix & w,datamatrix & data,
                          const unsigned & col,const bool & current) const
  {

  register unsigned i;
  unsigned p,k;

  unsigned nrconst = data.cols();

  compute_weight(w,col,current);

  double* workw = w.getV();
  double*workXp;
  double* workXk;

  for (p=0;p<nrconst;p++)
    for (k=p;k<nrconst;k++)
      {
      XWX(p,k)=0;
      workw = w.getV();
      workXp = data.getV()+p;
      workXk = data.getV()+k;
      for(i=0;i<nrobs;i++,workw++,workXp+=nrconst,workXk+=nrconst)
        XWX(p,k)+= *workw  *  *workXp * *workXk;
      XWX(k,p) = XWX(p,k);
      }

  }


double DISTRIBUTION::fisher2(const unsigned & beg, const unsigned & end,
                            const statmatrix<int> & index,
                            datamatrix & data, const unsigned & col,
                            const bool & current) const
  {

  register unsigned i;

  double * workdata = data.getV()+beg;

  double * worklin;
  if (current)
    worklin = (*linpred_current).getV()+col;
  else
    worklin = (*linpred_proposed).getV()+col;

  int * workind = index.getV()+beg;

  double XWX=0;
  double w;
  for(i=beg;i<=end;i++,workdata++,workind++)
    {
    w = weight(*workind,0);
    if (w > 0)
      XWX+= compute_weight(worklin+(*workind)*linearpred.cols(),
            &w,*workind,col) * (*workdata)*(*workdata);
    }
  return XWX;
  }

void DISTRIBUTION::compute_iwls(void)
  {
  iwlsweights_notchanged_df = false;
  tilde_y(tildey);
  compute_weight(weightiwls,0);
  }

void DISTRIBUTION::fisher(datamatrix & XWX,datamatrix & w,
                          vector<unsigned> & posbeg, vector<unsigned> & posend,
                          statmatrix<int> & index, unsigned & refind,
                          const unsigned & c) const
  {


  compute_weight(w,c);


  register unsigned i;
  unsigned p,k,l;

  unsigned nrconst = XWX.cols();

  for (p=0;p<nrconst;p++)
    for (k=p;k<nrconst;k++)
      {
      XWX(p,k)=0;

      unsigned ind=0;
      for(i=0;i<posbeg.size();i++)
        {
        if (i!=refind)
          {
          if ( (p==k) && (ind==k) )
            for(l=posbeg[i];l<=posend[i];l++)
              XWX(p,k) += w(index(l,0),0);
          ind++;
          }
        else
          {
          for(l=posbeg[i];l<=posend[i];l++)
            XWX(p,k) += w(index(l,0),0);
          }

        }  // end: for(i=0;i<posbeg.size();i++)

      XWX(k,p) = XWX(p,k);
      }

  }


void DISTRIBUTION::tilde_y(datamatrix & tildey,const bool & current)
  {
  register unsigned i;

  double * worklin;
  if (current == true)
    worklin = (*linpred_current).getV();
  else
    worklin = (*linpred_proposed).getV();

  double * workres = response.getV();
  double * ywork = tildey.getV();
  double mu;
  for (i=0;i<nrobs;i++,worklin++,ywork++,workres++)
    {
    compute_mu(worklin,&mu);
    *ywork = *worklin + (*workres - mu)*compute_gmu(worklin);
    }
  }


void DISTRIBUTION::tilde_y(datamatrix & tildey,datamatrix & m,
const unsigned & col,const bool & current,const datamatrix & w)
  {
  register unsigned i;
  unsigned dim = linearpred.cols();

  double * worklin;
  if (current == true)
    worklin = (*linpred_current).getV();
  else
    worklin = (*linpred_proposed).getV();

  double * workres = response.getV()+col;
  double * ywork = tildey.getV();
  double * mwork = m.getV();
  datamatrix muhelp = datamatrix(dim,1,0.0);
  double * muhelpp = muhelp.getV();
  for (i=0;i<nrobs;i++,worklin+=dim,ywork++,workres+=dim,mwork++)
    {
    compute_mu(worklin,muhelpp);
    *ywork = *mwork + (*workres - muhelp(col,0))*compute_gmu(worklin,col);
    }
  }


void DISTRIBUTION::tilde_y_minus_eta(datamatrix & tildey,
const unsigned & col,const bool & current)
  {
  register unsigned i;
  unsigned dim = linearpred.cols();

  double * worklin;
  if (current == true)
    worklin = (*linpred_current).getV();
  else
    worklin = (*linpred_proposed).getV();

  double * workres = response.getV()+col;
  double * ywork = tildey.getV();
  datamatrix muhelp = datamatrix(dim,1,0.0);
  double * muhelpp = muhelp.getV();
  for (i=0;i<nrobs;i++,worklin+=dim,ywork++,workres+=dim)
    {
    compute_mu(worklin,muhelpp);
    *ywork = (*workres - muhelp(col,0))*compute_gmu(worklin,col);
    }
  }


void DISTRIBUTION::compute_mu(const datamatrix & linpred, datamatrix & mu) const
  {

  assert (linpred.rows() == mu.rows());

  unsigned i;
  double * worklin = linpred.getV();
  double * workmu  = mu.getV();
  unsigned size = response.cols();
  double * muh;
  double * linh;
  for (i=0;i<linpred.rows();i++,workmu+=size,worklin+=size)
    {
    linh = worklin;
    muh = workmu;
    compute_mu(linh,muh);
    }

  }


void DISTRIBUTION::substr_linearpred(const datamatrix & m,const bool & current)
  {
  unsigned register i;
  double* workl;
  double* workm = m.getV();
  unsigned size=linearpred.rows()*linearpred.cols();

  if (current)
    workl = (*linpred_current).getV();
  else
    workl = (*linpred_proposed).getV();

    for(i=0;i<size;i++,workl++,workm++)
      *workl-=*workm;

  }


void DISTRIBUTION::substr_linearpred_m(const datamatrix & m,const unsigned & col,
                           const bool & current)
  {

  unsigned register i;
  double* workl;
  double* workm = m.getV();
  unsigned size=linearpred.cols();

  if (current)
    workl = (*linpred_current).getV()+col;
  else
    workl = (*linpred_proposed).getV()+col;

    for(i=0;i<nrobs;i++,workl+=size,workm++)
      *workl-=*workm;

  }


void DISTRIBUTION::add_linearpred(const datamatrix & m,const bool & current)
  {

  assert(m.rows() == linearpred.rows());
  assert(m.cols() == linearpred.cols());

  unsigned register i;
  double* workl;
  double* workm = m.getV();
  unsigned size=linearpred.rows()*linearpred.cols();

  if (current)
    workl = (*linpred_current).getV();
  else
    workl = (*linpred_proposed).getV();

  for(i=0;i<size;i++,workl++,workm++)
    *workl+=*workm;
  }


void DISTRIBUTION::add_linearpred(const double & m, unsigned & row, unsigned & col,
  const bool & current)
  {
  double * workl;
  if (current)
    workl = (*linpred_current).getV()+linearpred.cols()*row+col;
  else
    workl = (*linpred_proposed).getV()+linearpred.cols()*row+col;

  *workl += m;

  }


void DISTRIBUTION::add_linearpred_m(const double & m,const unsigned & col,
                        const bool & current)
  {
  unsigned register i;
  double* workl;
  unsigned size=linearpred.cols();

  if (current)
    workl = (*linpred_current).getV()+col;
  else
    workl = (*linpred_proposed).getV()+col;

  for(i=0;i<nrobs;i++,workl+=size)
    *workl+=m;

  }


void DISTRIBUTION::add_linearpred_m(const datamatrix & m,const unsigned & col,
                        const bool & current)
  {
  unsigned register i;
  double* workl;
  double* workm = m.getV();
  unsigned size=linearpred.cols();

  if (current)
    workl = (*linpred_current).getV()+col;
  else
    workl = (*linpred_proposed).getV()+col;

  for(i=0;i<nrobs;i++,workl+=size,workm++)
    *workl+=*workm;

  }


void DISTRIBUTION::add_linearpred(const double & m,const bool & current)
  {

  unsigned register i;
  double* workl;

  unsigned size=linearpred.rows()*linearpred.cols();

  if (current)
    workl = (*linpred_current).getV();
  else
    workl = (*linpred_proposed).getV();

  for(i=0;i<size;i++,workl++)
    *workl+=m;

  }




void DISTRIBUTION::add_linearpred(const datamatrix & m,const unsigned & beg,
                        const unsigned & end,const statmatrix<int> & index,
                        const unsigned & col, const bool & current)
  {
  unsigned register i;
  int * workindex = index.getV()+beg;
  double * workm = m.getV();
  if (current)
    for (i=beg;i<=end;i++,workindex++,workm++)
      (*linpred_current)(*workindex,col)+=*workm;
  else
    for (i=beg;i<=end;i++,workindex++,workm++)
      (*linpred_proposed)(*workindex,col)+=*workm;
  }


void DISTRIBUTION::add_linearpred(const double & m,const unsigned & beg,
                         const unsigned & end,const statmatrix<int> & index,
                         const unsigned & col, const bool & current)
  {
  unsigned register i;
  int * workindex = index.getV()+beg;
  if (current)
    for (i=beg;i<=end;i++,workindex++)
      (*linpred_current)(*workindex,col)+=m;
  else
    for (i=beg;i<=end;i++,workindex++)
      (*linpred_proposed)(*workindex,col)+=m;

  }


void DISTRIBUTION::add_linearpred2(const double & m,const unsigned & beg,
                         const unsigned & end,
                         const statmatrix<int> & index,
                         const statmatrix<int> & index2,
                         const unsigned & col,const bool & current)
  {

  unsigned register i;

  double * worklin;
  if (current)
    {
    worklin = &((*linpred_current)(index(beg,0),col));
    }
  else
    worklin = &((*linpred_proposed)(index(beg,0),col));

  int * workindex2 = index2.getV()+beg;
  int linearpred_cols = linearpred.cols();

  for (i=beg;i<=end;i++,workindex2++,worklin+=*workindex2*linearpred_cols)
    *worklin += m;

  }


void DISTRIBUTION::add_linearpred2(
                     const double & m,const unsigned & beg,
                     const unsigned & end,
                     const datamatrix & data,
                     const statmatrix<int> & index,
                     const statmatrix<int> & index2,
                     const unsigned & col,const bool & current)
  {


  double * worklin;
  if (current)
    {
    worklin = &((*linpred_current)(index(beg,0),col));
    }
  else
    worklin = &((*linpred_proposed)(index(beg,0),col));

  unsigned register i;

  int * workindex2 = index2.getV()+beg;
  double * workdata = data.getV()+beg;

  int * workindex = index.getV()+beg;

  for (i=beg;i<=end;i++,workindex2++,workdata++,workindex++,
       worklin+=*workindex2*linearpred.cols())
    *worklin += *workdata*m;

  }



void DISTRIBUTION::addtocurrent(const datamatrix & m)
  {

  unsigned register i;
  double* workc = (*linpred_current).getV();
  double * workp = (*linpred_proposed).getV();

  double* workm = m.getV();
  unsigned size=linearpred.rows()*linearpred.cols();

  for(i=0;i<size;i++,workc++,workm++,workp++)
    *workp = *workc + *workm;

  }


void DISTRIBUTION::addtocurrent(const double & m)
  {

  unsigned register i;
  double* workc = (*linpred_current).getV();
  double * workp = (*linpred_proposed).getV();

  unsigned size=linearpred.rows()*linearpred.cols();

  for(i=0;i<size;i++,workc++,workp++)
    *workp = *workc + m;

  }


void DISTRIBUTION::addtocurrentcol(const datamatrix & m,const unsigned & col)
  {

  unsigned register i,j;
  double* workc = (*linpred_current).getV();
  double * workp = (*linpred_proposed).getV();

  double* workm = m.getV();
  unsigned size=linearpred.cols();

  for(i=0;i<nrobs;i++,workm++)
    {
    for(j=0;j<size;j++,workc++,workp++)
      {
      if (j==col)
        *workp = *workc + *workm;
      else
        *workp = *workc;
      }

    }

  }


void DISTRIBUTION::addtocurrentcol(const double & m,const unsigned & col)
  {

  unsigned register i,j;
  double* workc = (*linpred_current).getV();
  double * workp = (*linpred_proposed).getV();


  unsigned size=linearpred.cols();

  for(i=0;i<nrobs;i++)
    {
    for(j=0;j<size;j++,workc++,workp++)
      {
      if (j==col)
        *workp = *workc + m;
      else
        *workp = *workc;
      }

    }

  }



void DISTRIBUTION::addtocurrentcol_single(const double & m,const unsigned & r,
                                   const unsigned & col)
  {
  unsigned j;
  unsigned size=linearpred.cols();

  for(j=0;j<size;j++)
    {
    if (j==col)
      (*linpred_proposed)(r,col) = (*linpred_current)(r,col)+m;
    else
      (*linpred_proposed)(r,col) = (*linpred_current)(r,col);
    }

  }



void DISTRIBUTION::addtocurrentcol(const double & m,const unsigned & beg,
                                   const unsigned & end,
                                   const statmatrix<int> & index,
                                   const unsigned & col)
  {

  unsigned register i,j;
  int * workindex = index.getV()+beg;
  unsigned size=linearpred.cols();

  for (i=beg;i<=end;i++,workindex++)
    {

    for(j=0;j<size;j++)
      {
      if (j==col)
        (*linpred_proposed)(*workindex,j) =
        (*linpred_current)(*workindex,j) + m;
      else
        (*linpred_proposed)(*workindex,j) =
        (*linpred_current)(*workindex,j);
      }

    }

  }


void DISTRIBUTION::assign(const bool & current)
  {
  unsigned register i;

  double * workc;
  double * workp;

  if (current)
    {
    workc = linpred_current->getV();
    workp = linpred_proposed->getV();
    }
  else
    {
    workc = linpred_proposed->getV();
    workp = linpred_current->getV();
    }

  unsigned size = linearpred.cols()*linearpred.rows();

  for (i=0;i<size;i++,workc++,workp++)
    *workc = *workp;

  }


void DISTRIBUTION::assigncol(const unsigned & col,const bool & current)
  {

  unsigned register i;

  double * workc;
  double * workp;

  if (current)
    {
    workc = linpred_current->getV()+col;
    workp = linpred_proposed->getV()+col;
    }
  else
    {
    workc = linpred_proposed->getV()+col;
    workp = linpred_current->getV()+col;
    }

  unsigned size = linearpred.rows();

  for (i=0;i<size;i++,workc+=linearpred.cols(),workp+=linearpred.cols())
    *workc = *workp;

  }


void DISTRIBUTION::assign(const unsigned & beg,const unsigned & end,
              const statmatrix<int> & index,const bool & current)
  {

  unsigned register i;
  int * workindex = index.getV()+beg;
  if (current)
    {
    for (i=beg;i<=end;i++,workindex++)
      (*linpred_current)(*workindex,0) = (*linpred_proposed)(*workindex,0);
    }
  else
    {
    for (i=beg;i<=end;i++,workindex++)
      (*linpred_proposed)(*workindex,0) = (*linpred_current)(*workindex,0);
    }

  }


void DISTRIBUTION::compute_respminuslinpred(datamatrix & res,
                                            const unsigned & co)
  {
  res.minus(response,(*linpred_current),co,co);
  }


void DISTRIBUTION::compute_workingresiduals(const unsigned & co)
  {
  workingres.minus(tildey,(*linpred_current),co,co);
  }




void DISTRIBUTION::compute_weightiwls_workingresiduals(const unsigned & co)
  {
  unsigned i;
  double * worktildey = tildey.getV()+co;
  double * workweightiwls = weightiwls.getV()+co;
  double * worklin = (*linpred_current).getV()+co;
  double * workres = workingres.getV();
  unsigned c = weightiwls.cols();

  for(i=0;i<nrobs;i++,worktildey+=c,workweightiwls+=c,worklin+=c,workres++)
    *workres = *workweightiwls * (*worktildey - *worklin);

  }


void DISTRIBUTION::compute_weightiwls_workingresiduals(datamatrix & tildey,
                datamatrix & weightiwls, const unsigned & co)
  {
  unsigned i;
  double * worktildey = tildey.getV();
  double * workweightiwls = weightiwls.getV();
  double * worklin = (*linpred_current).getV()+co;
  double * workres = workingres.getV();
  unsigned c = (*linpred_current).cols();

  for(i=0;i<nrobs;i++,worktildey++,workweightiwls++,worklin+=c,workres++)
    *workres = *workweightiwls * (*worktildey - *worklin);

  }


void DISTRIBUTION::swap_linearpred(void)
  {
  datamatrix * help = linpred_current;
  linpred_current = linpred_proposed;
  linpred_proposed = help;
  }


    void
    DISTRIBUTION::update(void)
    {

        if (scaleexisting)
        {

            double * bp = Scalesave.getbetapointer();
            double * sc = scale.getV();
            unsigned i, j;
            for (i = 0; i < scale.rows(); i++)
                for (j = 0; j < scale.cols(); j++, bp++, sc++)
                    * bp = * sc;

            Scalesave.updatemult();

        } // end: if (scaleexisting)

        // BEGIN: DSB //
        #if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
        // do we need to do mschecking?
        if (mscheck)
        {

            // add the current linear predictor values to the predictor changes values
            {
                double * prc = predchange.getV();
                double * prcu = linpred_current->getV();
                for (unsigned int i = 0; i < nrobs; i++, prc++, prcu++)
                {
                    * prc += * prcu;
                }
            }
            // so now predchange contains the approximate leave-one-out predictor samples!

            // write the approximate LOO predictor values into the output matrix,
            // and compute the corresponding log-likelihood contributions,
            // for all (single) observations.
            datamatrix singleLogLikelihoods = datamatrix(nrobs, 1);
            {
                double * msc_predp = msc_pred.getbetapointer();
                double * msc_obssamplesp = msc_obssamples.getbetapointer();
                double * resp = response.getV();
                double * wei = weight.getV();
                double * prc = predchange.getV();
                double * singleLogLikelihoodsp = singleLogLikelihoods.getV();

                for (unsigned int i = 0; i < nrobs;
                     i++, singleLogLikelihoodsp++, msc_predp++, msc_obssamplesp++, resp++, wei++, prc++)
                {
                    // first compute the corresponding mu (with transformation onto the original scale,
                    // as the functions below expect that.)
                    double mu = 0;
                    compute_mu(prc, &mu);

                    // compute the loglikelihood of the observation *resp
                    * singleLogLikelihoodsp = loglikelihood_from_deviance(*resp, mu, *wei);

                    // sample once from the corresponding predictive distribution
                    * msc_obssamplesp = sample_from_likelihood(*wei, mu);

                    // and finally compute the linear predictor sample
                    double factor = trmult(0,0);
                    double shift = addinterceptsample;
                    * msc_predp = factor * (*prc) + shift;
                }
            }


            double * msc_likep = msc_like.getbetapointer();

            for(unsigned int i=0; i<msnrind; i++, msc_likep++)
                {

                // first sum all loglikelihood contributions for individual i
                int * workindex = msindex.getV() + msposbeg[i];
                *msc_likep = 0.0;

                for (unsigned register j = msposbeg[i]; j <= msposend[i]; j++, workindex++)
                   {
                   *msc_likep += singleLogLikelihoods(*workindex,0);
                   }

                // and then exponentiate to get the likelihood of individual i
                *msc_likep = exp(*msc_likep);
                }

            // write into files
            msc_like.update();
            msc_pred.update();
            msc_obssamples.update();

            // refresh predchange vector with zeros
            predchange = datamatrix(nrobs, 1, 0);
        }
        #endif
        // END: DSB //

    } // end: update


void DISTRIBUTION::update_predict(void)
  {

  if (predict)
    {

    unsigned samplesize = optionsp->get_samplesize();

    if(
      (optionsp->get_nriter() > optionsp->get_burnin())
      &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1)
       % (optionsp->get_step()) == 0)
      )
      {

      register unsigned i,j;
      double * worklin = (*linpred_current).getV();
      double * workmean = linpredmean.getV();
      double * workmumean = mumean.getV();
      datamatrix muhelp(mumean.cols(),1,0);
      double * workdevmean = deviancemean.getV();
      double * workdevmean_sat = deviancemean_sat.getV();
      double * workresponse = response.getV();
      double * workw = weight.getV();

      unsigned size = linearpred.cols();
      unsigned size2 = mumean.cols();
      unsigned respsize = response.cols();

      double * mumeanhelp;
      double reshelp;

      double devhelp;

      double * musavep = musave.getbetapointer();



      if (samplesize==1)
        {
        for (i=0;i<nrobs;i++,workdevmean++,workdevmean_sat++,workmumean+=size2,
                                      workresponse+=respsize,workw++)
          {
          mumeanhelp = workmumean;
          compute_mu(worklin,mumeanhelp);
          compute_deviance(workresponse,workw,workmumean,workdevmean,
                           &devhelp,scale,i);
          deviance(0,0) += *workdevmean;
          deviance(0,1) += devhelp;

          *workdevmean_sat = devhelp;
          for(j=0;j<size;j++,worklin++,workmean++)
            *workmean = trmult(j,0) * *worklin;

          if ((predictfull) && (i<firstobs))
            {
            mumeanhelp = workmumean;
            for(j=0;j<size2;j++,musavep++,mumeanhelp++)
              {
              *musavep = *mumeanhelp;
              }
            }

          }

        }  // end: (samplesize==1)
      else
        {
        for (i=0;i<nrobs;i++,workresponse+=respsize,workw++,workdevmean++,
                             workdevmean_sat++)
          {
          mumeanhelp = muhelp.getV();
          compute_mu(worklin,mumeanhelp);
          mumeanhelp = muhelp.getV();

          compute_deviance(workresponse,workw,mumeanhelp,&reshelp,&devhelp,scale,i);
          deviance(samplesize-1,0) += reshelp;
          deviance(samplesize-1,1) += devhelp;

          mumeanhelp = muhelp.getV();

          for(j=0;j<size;j++,workmean++,worklin++)
            {
            *workmean = (1.0)/samplesize * ( (samplesize-1)* *workmean +
                         *worklin*trmult(j,0));
            }

          for(j=0;j<size2;j++,workmumean++,mumeanhelp++)
            {
            *workmumean = (1.0)/samplesize * ( (samplesize-1)* *workmumean +
                           *mumeanhelp);
            }


          *workdevmean = (1.0)/samplesize * ( (samplesize-1)* *workdevmean +
                          reshelp);

          *workdevmean_sat = (1.0)/samplesize * ( (samplesize-1)*
                             *workdevmean_sat + devhelp);

          if ((predictfull) && (i<firstobs))
            {
            mumeanhelp = muhelp.getV();
            for(j=0;j<size2;j++,musavep++,mumeanhelp++)
              {
              *musavep = *mumeanhelp;
              }
            }  // end: if (predictfull)


          }

        }

      }

     if (predictfull)
       musave.update();

    } // end: if predict

  }


void DISTRIBUTION::update_bootstrap(void)
  {
  if (scaleexisting)
    {
    DISTRIBUTION::update();
    } // end: if (scaleexisting)
  }


void DISTRIBUTION::update_predict_bootstrap(int & bootstrapsamples)
  {
  if (predict)
    {
    unsigned samplesize = optionsp->get_samplesize();
    if(samplesize == 1)
       deviance = datamatrix((bootstrapsamples + 1),2,0);

    DISTRIBUTION::update_predict();
    }
  }


void DISTRIBUTION::transform_nonlinear(vector<FULLCOND *> & fc,ST::string & trtype)
  {

  unsigned nrfc = fc.size();

  vector<datamatrix> sample(nrfc);

  unsigned i,j,k;

  vector<double *>  wmean(nrfc);
  vector<double *>  wqu1l(nrfc);
  vector<double *>  wqu2l(nrfc);
  vector<double *>  wqu50(nrfc);
  vector<double *>  wqu1u(nrfc);
  vector<double *>  wqu2u(nrfc);
  vector<double *>  wvar(nrfc);

  for(i=0;i<fc.size();i++)
    {
    sample[i] = datamatrix(optionsp->get_samplesize(),1);
    wmean[i] = fc[i]->get_betameanp();
    wqu1l[i] = fc[i]->get_beta_lower1_p();
    wqu2l[i] = fc[i]->get_beta_lower2_p();
    wqu50[i] = fc[i]->get_betaqu50p();
    wqu1u[i] = fc[i]->get_beta_upper2_p();
    wqu2u[i] = fc[i]->get_beta_upper1_p();
    wvar[i]  = fc[i]->get_betavarp();
    }

  vector<double *> wsample(nrfc);

  unsigned nrpar = fc[0]->get_nrpar();
  double lower1 = fc[0]->get_lower1();
  double lower2 = fc[0]->get_lower2();
  double upper1 = fc[0]->get_upper1();
  double upper2 = fc[0]->get_upper2();

  for(i=0;i<nrpar;i++)
    {
    for(k=0;k<nrfc;k++)
      {
      fc[k]->readsample(sample[k],i);
      wsample[k] = sample[k].getV();
      }

    for (j=0;j<sample[0].rows();j++)
      {
      tr_nonlinear(wsample,wsample,fc,i,j,trtype);

      if (j <sample[0].rows()-1)
        {
        for (k=0;k<nrfc;k++)
          {
          wsample[k]++;
          }

        }

      }

    for(k=0;k<nrfc;k++)
      {

      *(wmean[k]) = sample[k].mean(0);
      *(wqu1l[k]) = sample[k].quantile(lower1,0);
      *(wqu2l[k]) = sample[k].quantile(lower2,0);
      *(wqu50[k]) = sample[k].quantile(50,0);
      *(wqu1u[k]) = sample[k].quantile(upper1,0);
      *(wqu2u[k]) = sample[k].quantile(upper2,0);
      *(wvar[k]) = sample[k].var(0);

      }

    if (i < nrpar-1)
      {
      for (k=0;k<nrfc;k++)
        {
        wmean[k]++;
        wqu1l[k]++;
        wqu2l[k]++;
        wqu50[k]++;
        wqu1u[k]++;
        wqu2u[k]++;
        wvar[k]++;
        }

      }

    } // end: for(i=0;i<nrpar;i++)


  }


void DISTRIBUTION::outresults(void)
  {

  char hchar = '%';
  ST::string helpstring = "\\%";

  double lower1;
  double lower2;
  double upper1;
  double upper2;

  if (scaleexisting)
    {
    Scalesave.outresults();
    }

  lower1 = Scalesave.get_lower1();
  lower2 = Scalesave.get_lower2();
  upper1 = Scalesave.get_upper1();
  upper2 = Scalesave.get_upper2();

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);

  ST::string nl1 = ST::doubletostring(lower1,4);
  ST::string nl2 = ST::doubletostring(lower2,4);
  ST::string nu1 = ST::doubletostring(upper1,4);
  ST::string nu2 = ST::doubletostring(upper2,4);
  nl1 = nl1.replaceallsigns('.','p');
  nl2 = nl2.replaceallsigns('.','p');
  nu1 = nu1.replaceallsigns('.','p');
  nu2 = nu2.replaceallsigns('.','p');

  if (fcmissing.size() > 0)
    {
    MissingSave.outresults();

    if (optionsp->get_samplesize() > 0)
      {

      ofstream outmissing(pathmissing.strtochar());

      optionsp->out("  Missing values:\n",true);
      optionsp->out("\n");
      optionsp->out("  Estimates for missing values are stored in file\n");
      optionsp->out("  " + pathmissing + "\n");
      optionsp->out("\n");

      outmissing << "obsnr   pmean   pstddev   pqu" << nl1 << "   pqu" << nl2 <<
                    "   pqu50   pqu" << nu1 << "   pqu" << nu2 << endl;

      double * workmean = MissingSave.get_betameanp();
      double * workvar = MissingSave.get_betavarp();
      double * workbetaqu_l1_lower_p = MissingSave.get_beta_lower1_p();
      double * workbetaqu_l2_lower_p = MissingSave.get_beta_lower2_p();
      double * workbetaqu_l1_upper_p = MissingSave.get_beta_upper1_p();
      double * workbetaqu_l2_upper_p = MissingSave.get_beta_upper2_p();
      double * workbetaqu50 = MissingSave.get_betaqu50p();
      unsigned nr=1;
      unsigned i;
      for(i=0;i<missingpos.rows();i++,workmean++,workvar++,
                                  workbetaqu_l1_lower_p++,
                                  workbetaqu_l2_lower_p++,
                                  workbetaqu_l1_upper_p++,
                                  workbetaqu_l2_upper_p++,
                                  workbetaqu50++)
        {
        nr += missingpos(i,0);
        outmissing << nr << "   ";
        outmissing << *workmean << "   ";
        outmissing << *workvar << "   ";
        outmissing << *workbetaqu_l1_lower_p << "   ";
        outmissing << *workbetaqu_l2_lower_p << "   ";
        outmissing << *workbetaqu50 << "   ";
        outmissing << *workbetaqu_l2_upper_p << "   ";
        outmissing << *workbetaqu_l1_upper_p << "   ";

        outmissing << endl;
        }

      }

    }

    // BEGIN: DSB //
    #if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
    // if we do mschecking
  if (mscheck)
    {
        msc_like.outresults();
        // and then the predictor samples
        msc_pred.outresults();
        // then the predictive observation samples
        msc_obssamples.outresults();

    double avescore = 0;

    optionsp->out("  Marshall-Spiegelhalter Cross Validation:\n",true);
    optionsp->out("\n");
    optionsp->out("  Estimated individual likelihoods, predictors and\n");
    optionsp->out("  observation samples are stored in files\n");
    optionsp->out("\n");

    // Likelihood results
    ST::string pathhelp =
        pathresultsscale.substr(0, pathresultsscale.length() - 10)
                        + "_mscheck_like.res";
    optionsp->out("  " + pathhelp + "\n");
    optionsp->out("\n");

    ofstream outres(pathhelp.strtochar());

    outres << "intnr" << "   ";
    outres << "pmean   ";
    outres << "pqu"  << l1  << "   ";
    outres << "pqu"  << l2  << "   ";
    outres << "pmed   ";
    outres << "pqu"  << u1  << "   ";
    outres << "pqu"  << u2  << "   ";
    outres << endl;

    // here is the pointer to the means of the individuals likelihood samples:
    double * workmean = msc_like.get_betameanp();

    double * workbetaqu_l1_lower_p = msc_like.get_beta_lower1_p();
    double * workbetaqu_l2_lower_p = msc_like.get_beta_lower2_p();
    double * workbetaqu50 = msc_like.get_betaqu50p();
    double * workbetaqu_l1_upper_p = msc_like.get_beta_upper1_p();
    double * workbetaqu_l2_upper_p = msc_like.get_beta_upper2_p();

    for(unsigned int i=0; i<msnrind; i++, workmean++, workbetaqu_l1_lower_p++,
    workbetaqu_l2_lower_p++, workbetaqu50++, workbetaqu_l2_upper_p++, workbetaqu_l1_upper_p++)
      {
        // for the average log-score, we first have to logarithmize the individuals mean likelihoods,
        // and then average this over all individuals
      avescore += log(*workmean);

      outres << (i+1) << "   ";
      outres << *workmean << "   ";
      outres << *workbetaqu_l1_lower_p << "   ";
      outres << *workbetaqu_l2_lower_p << "   ";
      outres << *workbetaqu50 << "   ";
      outres << *workbetaqu_l2_upper_p << "   ";
      outres << *workbetaqu_l1_upper_p << "   ";
      outres << endl;
      }

    // divide by the number of individuals to get the average log-score,
    // and also add a minus to get the same orientation as for the DIC (lower scores are then better)
    avescore /= - (double)msnrind;

    // Predictor results
    pathhelp = pathresultsscale.substr(0, pathresultsscale.length()
                - 10) + "_mscheck_pred.res";
    optionsp->out("  " + pathhelp + "\n");
    optionsp->out("\n");

    ofstream outres2(pathhelp.strtochar());

    outres2 << "intnr" << "   ";
    outres2 << "pmean   ";
    outres2 << "pqu"  << l1  << "   ";
    outres2 << "pqu"  << l2  << "   ";
    outres2 << "pmed   ";
    outres2 << "pqu"  << u1  << "   ";
    outres2 << "pqu"  << u2  << "   ";
    outres2 << endl;

    workmean = msc_pred.get_betameanp();
    workbetaqu_l1_lower_p = msc_pred.get_beta_lower1_p();
    workbetaqu_l2_lower_p = msc_pred.get_beta_lower2_p();
    workbetaqu50 = msc_pred.get_betaqu50p();
    workbetaqu_l1_upper_p = msc_pred.get_beta_upper1_p();
    workbetaqu_l2_upper_p = msc_pred.get_beta_upper2_p();

    for(unsigned int i=0; i<nrobs; i++, workmean++, workbetaqu_l1_lower_p++,
    workbetaqu_l2_lower_p++, workbetaqu50++, workbetaqu_l2_upper_p++, workbetaqu_l1_upper_p++)
      {
      outres2 << (i+1) << "   ";
      outres2 << *workmean << "   ";
      outres2 << *workbetaqu_l1_lower_p << "   ";
      outres2 << *workbetaqu_l2_lower_p << "   ";
      outres2 << *workbetaqu50 << "   ";
      outres2 << *workbetaqu_l2_upper_p << "   ";
      outres2 << *workbetaqu_l1_upper_p << "   ";
      outres2 << endl;
      }

    // Observation samples

    pathhelp = pathresultsscale.substr(0, pathresultsscale.length() - 10)
                       + "_mscheck_obssamples.res";
    optionsp->out("  " + pathhelp + "\n\n\n");

    ofstream outres3(pathhelp.strtochar());

    outres3 << "intnr" << "   ";
    outres3 << "pmean   ";
    outres3 << "pqu"  << l1  << "   ";
    outres3 << "pqu"  << l2  << "   ";
    outres3 << "pmed   ";
    outres3 << "pqu"  << u1  << "   ";
    outres3 << "pqu"  << u2  << "   ";
    outres3 << endl;

    workmean = msc_obssamples.get_betameanp();
    workbetaqu_l1_lower_p = msc_obssamples.get_beta_lower1_p();
    workbetaqu_l2_lower_p = msc_obssamples.get_beta_lower2_p();
    workbetaqu50 = msc_obssamples.get_betaqu50p();
    workbetaqu_l1_upper_p = msc_obssamples.get_beta_upper1_p();
    workbetaqu_l2_upper_p = msc_obssamples.get_beta_upper2_p();

    for(unsigned int i=0; i<nrobs; i++, workmean++, workbetaqu_l1_lower_p++,
    workbetaqu_l2_lower_p++, workbetaqu50++, workbetaqu_l2_upper_p++, workbetaqu_l1_upper_p++)
      {
      outres3 << (i+1) << "   ";
      outres3 << *workmean << "   ";
      outres3 << *workbetaqu_l1_lower_p << "   ";
      outres3 << *workbetaqu_l2_lower_p << "   ";
      outres3 << *workbetaqu50 << "   ";
      outres3 << *workbetaqu_l2_upper_p << "   ";
      outres3 << *workbetaqu_l1_upper_p << "   ";
      outres3 << endl;
      }

    // Results for the log-score
    optionsp->out("\n");
    optionsp->out("  Resulting average leave one out log-score: " + ST::doubletostring(avescore, 6) + "\n");
    optionsp->out("\n");
    optionsp->out("\n");

    }
    #endif
    // END: DSB //


  if ( (predict) && (optionsp->get_samplesize() > 0) && (!isbootstrap) )
    {



    register unsigned i,j;

    optionsp->out("  Predicted values:\n",true);
    optionsp->out("\n");
    optionsp->out("  Estimated mean of predictors, expectation of response and\n");
    optionsp->out("  individual deviances are stored in file\n");
    optionsp->out("  " + predictpath + "\n");
    optionsp->out("\n");

    if (predictfull)
      {

      optionsp->out("  Estimation results for the expectation of responses are stored in file\n");
      optionsp->out("  " + predictfullpath + "\n");
      optionsp->out("\n");

      musave.outresults();

/*
// speichert die samples der mu's

      datamatrix s(optionsp->get_samplesize(),nrobs*mumean.cols());

      musave.readsample3(s);

      ST::string predictfullpathsample = predictfullpath.substr(0,predictfullpath.length()-4)+"_sample.raw";

      ofstream outmusample(predictfullpathsample.strtochar());

      for (i=0;i< nrobs;i++)
        {
        for (j=0;j<mumean.cols();j++)
          {

          outmusample << "mu" << (j+1) << "_" << (i+1) << "  ";

          }

        }

      outmusample << endl;

      for (i=0;i<s.rows();i++)
        {
        for (j=0;j<s.cols();j++)
          {

          outmusample << s(i,j) << "  ";

          }
         outmusample << endl;
        }

*/

      ofstream outmu(predictfullpath.strtochar());

      outmu << "intnr   ";

      for (j=0;j<Dnames.size();j++)
        outmu << Dnames[j] << "   ";

      if (mumean.cols() > 1)
        {
        for(j=1;j<=mumean.cols();j++)
          {
          outmu << "pmean_" << j << "   ";
          outmu << "stddev_" << j << "   ";
          outmu << "pqu" << nl1  << "_" << j << "   ";
          outmu << "pqu" << nl2  << "_" << j << "   ";
          outmu << "pqu50_"  << j << "   ";
          outmu << "pqu" << nu1  << "_" << j << "   ";
          outmu << "pqu" << nu2  << "_" << j << "   ";
          }
        }
      else
        {
        for(j=1;j<=mumean.cols();j++)
          {
          outmu << "pmean" << "   ";
          outmu << "stddev" << "   ";
          outmu << "pqu" << nl1  << "   ";
          outmu << "pqu" << nl2  << "   ";
          outmu << "pqu50" <<  "   ";
          outmu << "pqu" << nu1  <<  "   ";
          outmu << "pqu" << nu2  << "   ";
          }
        }

      outmu << endl;

      double * workmean = musave.get_betameanp();
      double * workbetaqu_l1_lower_p = musave.get_beta_lower1_p();
      double * workbetaqu_l2_lower_p = musave.get_beta_lower2_p();
      double * workbetaqu_l1_upper_p = musave.get_beta_upper1_p();
      double * workbetaqu_l2_upper_p = musave.get_beta_upper2_p();
      double * workbetaqu50 = musave.get_betaqu50p();
      double * workbetastd = musave.get_betavarp();

      double * datap = Dp->getV();
      unsigned sD = Dp->cols();

      for(i=0;i<firstobs;i++)
        {
        outmu << (i+1) << "   ";

        for(j=0;j<sD;j++,datap++)
          outmu << (*datap) << "   ";

        for (j=0;j<mumean.cols();j++,workmean++,workbetaqu_l1_lower_p++,
             workbetaqu_l2_lower_p++,workbetaqu50++, workbetaqu_l1_upper_p++,
             workbetaqu_l2_upper_p++,workbetastd++)
          {

          outmu << *workmean << "   ";
          outmu << (sqrt(*workbetastd)) << "   ";
          outmu << *workbetaqu_l1_lower_p << "   ";
          outmu << *workbetaqu_l2_lower_p << "   ";
          outmu << *workbetaqu50 << "   ";
          outmu << *workbetaqu_l2_upper_p << "   ";
          outmu << *workbetaqu_l1_upper_p << "   ";

          }

        outmu << endl;
        }

      } // end: if (predictfull)


    ofstream out(predictpath.strtochar());
    ofstream out2(deviancepath.strtochar());

    for(i=0;i<Dnames.size();i++)
      out << Dnames[i] << "   ";

    unsigned size1 = linearpred.cols();
    unsigned size2 = mumean.cols();

    if (size1 > 1)
      {
      for (i=0;i<size1;i++)
        {
        out << "linpred_" << (i+1) << "   ";
        }
      }
    else
      {
      out << "linpred" << "   ";
      }

    if (size2 > 1)
      {
      for (i=0;i<size2;i++)
        {
        out << "mu_" << (i+1) << "   ";
        }
      }
    else
      {
      out << "mu" << "   ";
      }

    if(family!="cox")
      out << "saturated_deviance" << "   ";
    else
      out << "unstandardized_deviance" << "   ";

    out << "leverage" << "   ";

    out << endl;

    double * workmean = linpredmean.getV();
    double * workmumean = mumean.getV();
    double * workdevmean = deviancemean.getV();

    double deviance2=0;
    double deviance2_sat=0;

    double * datap = Dp->getV();
    unsigned sD = Dp->cols();

    // used for computing DIC
    datamatrix scalehelp;
    if (scaleexisting)
      {
      scalehelp = Scalesave.get_betamean();

      for(i=0;i<scalehelp.rows();i++)
        for (j=0;j<scalehelp.cols();j++)
          scalehelp(i,j) = (1.0/(trmult(i,0)*trmult(j,0)))*scalehelp(i,j);
      }

    double reshelp;
    double devhelp;
    double * workmean_firstcol;
    datamatrix mu_meanlinpred(1,size2);


    for (i=0;i<nrobs;i++,workdevmean++)
      {

      for(j=0;j<sD;j++,datap++)
        out << *datap << "   ";

      workmean_firstcol = workmean;
      for (j=0;j<size1;j++,workmean++)
        {
        out << *workmean << "   ";
        }

      for (j=0;j<size2;j++,workmumean++)
        {
        out << *workmumean << "   ";
        }

      compute_mu_notransform(workmean_firstcol,&mu_meanlinpred(0,0));

      compute_deviance(&response(i,0),&weight(i,0),&mu_meanlinpred(0,0),
      &reshelp,&devhelp,scalehelp,i);

//      compute_deviance(&response(i,0),&weight(i,0),&mumean(i,0),
//      &reshelp,&devhelp,scalehelp,i);

      deviance2 += reshelp;
      deviance2_sat += devhelp;

      if (weight(i,0) != 0 && family!="cox")
        {
        out << devhelp  << "   ";

        out << (*workdevmean-reshelp) << "   ";
        }
      else if (weight(i,0) != 0 && family=="cox")
        {
        out << reshelp << "   ";

        out << ".   ";
        }
      else
        {
        out << ".   .   ";
        }

      out << endl;
      }

    double devhelpm;

    ST::string meanstr = "  Mean:          ";
    unsigned l_meanstr = meanstr.length();

    ST::string stdstr =  "  Std. Dev:      ";
    unsigned l_stdstr = stdstr.length();

    ST::string l1str = "  " + l1 + "% Quantile: ";
    unsigned l_l1str = l1str.length();

    ST::string l2str = "  " + l2 + "% Quantile: ";
    unsigned l_l2str = l2str.length();

    ST::string medianstr = "  50% Quantile: ";
    unsigned l_medianstr = medianstr.length();

    ST::string u1str = "  " + u1 + "% Quantile: ";
    unsigned l_u1str = u1str.length();

    ST::string u2str = "  " + u2 + "% Quantile: ";
    unsigned l_u2str = u2str.length();


    optionsp->out("  Estimation results for the deviance: \n",true);
    optionsp->out("\n");

    results_latex.push_back("\n {\\bf \\large Estimation results for the deviance: }\\\\ \n");
    results_latex.push_back("{\\bf Unstandardized deviance } \n");
    results_latex.push_back("\\vspace{-0.4cm}");
    results_latex.push_back("\\begin{tabbing}");
    results_latex.push_back("\\hspace{3cm} \\= \\\\");


    optionsp->out("  Unstandardized Deviance (-2*Loglikelihood(y|mu))\n");
    optionsp->out("\n");

    devhelpm = deviance.mean(0);

    unsigned d;
    if (devhelpm > 1000000000)
      d = 14;
    else if (devhelpm > 1000000)
      d = 11;
    else
      d = 8;
    optionsp->out(meanstr + ST::string(' ',20-l_meanstr) +
    ST::doubletostring(devhelpm,d) + "\n");

    results_latex.push_back(meanstr + " \\> " + ST::doubletostring(devhelpm,d)
                           + " \\\\");

    devhelp = sqrt(deviance.var(0));
    optionsp->out(stdstr + ST::string(' ',20-l_stdstr) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(stdstr + " \\> " + ST::doubletostring(devhelp,d)
                           + " \\\\");


    devhelp = deviance.quantile(lower1,0);
    optionsp->out(l1str +  ST::string(' ',20-l_l1str) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(l1str.insert_string_char(hchar,helpstring) + " \\> "
                           + ST::doubletostring(devhelp,d)
                           + " \\\\");


    devhelp = deviance.quantile(lower2,0);
    optionsp->out(l2str +  ST::string(' ',20-l_l2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(l2str.insert_string_char(hchar,helpstring)
                           + " \\> " + ST::doubletostring(devhelp,d)
                           + " \\\\");

    devhelp = deviance.quantile(50,0);
    optionsp->out(medianstr +  ST::string(' ',20-l_medianstr) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(medianstr.insert_string_char(hchar,helpstring)
                           + " \\> " + ST::doubletostring(devhelp,d)
                           + " \\\\");

    devhelp = deviance.quantile(upper1,0);
    optionsp->out(u1str +  ST::string(' ',20-l_u1str) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(u1str.insert_string_char(hchar,helpstring)
                           + " \\> " + ST::doubletostring(devhelp,d)
                           + " \\\\");

    devhelp = deviance.quantile(upper2,0);
    optionsp->out(u2str +  ST::string(' ',20-l_u2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(u2str.insert_string_char(hchar,helpstring)
                           + " \\> " + ST::doubletostring(devhelp,d)
                           + " \\\\");

    optionsp->out("\n");

    results_latex.push_back("\\end{tabbing}\n");

    if(family!="cox")
    {
    results_latex.push_back("{\\bf Saturated deviance } \n");
    results_latex.push_back("\\vspace{-0.4cm}");
    results_latex.push_back("\\begin{tabbing}");
    results_latex.push_back("\\hspace{3cm} \\= \\\\");

    optionsp->out("  Saturated Deviance (-2*Loglikelihood(y|mu) + 2*Loglikelihood(y|mu=y))\n");
    optionsp->out("\n");

    devhelpm = deviance.mean(1);

    optionsp->out(meanstr + ST::string(' ',20-l_meanstr) +
    ST::doubletostring(devhelpm,d) + "\n");

    results_latex.push_back(meanstr + " \\> "
                           + ST::doubletostring(devhelpm,d)
                           + " \\\\");

    devhelp = sqrt(deviance.var(1));
    optionsp->out(stdstr + ST::string(' ',20-l_stdstr) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(stdstr + " \\> "
                           + ST::doubletostring(devhelp,d)
                           + " \\\\");


    devhelp = deviance.quantile(lower1,1);
    optionsp->out(l1str +  ST::string(' ',20-l_l1str) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(l1str.insert_string_char(hchar,helpstring) + " \\> "
                           + ST::doubletostring(devhelp,d)
                           + " \\\\");


    devhelp = deviance.quantile(lower2,1);
    optionsp->out(l2str +  ST::string(' ',20-l_l2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(l2str.insert_string_char(hchar,helpstring) + " \\> "
                           + ST::doubletostring(devhelp,d)
                           + " \\\\");

    devhelp = deviance.quantile(50,1);
    optionsp->out(medianstr +  ST::string(' ',20-l_medianstr) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(medianstr.insert_string_char(hchar,helpstring) + " \\> "
                           + ST::doubletostring(devhelp,d)
                           + " \\\\");


    devhelp = deviance.quantile(upper1,1);
    optionsp->out(u1str +  ST::string(' ',20-l_u1str) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(u1str.insert_string_char(hchar,helpstring) + " \\> "
                           + ST::doubletostring(devhelp,d)
                           + " \\\\");


    devhelp = deviance.quantile(upper2,1);
    optionsp->out(u2str +  ST::string(' ',20-l_u2str) +
    ST::doubletostring(devhelp,d) +  "\n");

    results_latex.push_back(u2str.insert_string_char(hchar,helpstring) + " \\> "
                           + ST::doubletostring(devhelp,d)
                           + " \\\\");

    results_latex.push_back("\\end{tabbing}\n");

    optionsp->out("\n");
    }

    optionsp->out("  Samples of the deviance are stored in file\n");
    optionsp->out("  " + deviancepath + "\n");
    optionsp->out("\n");


    optionsp->out("  Estimation results for the DIC: \n",true);
    optionsp->out("\n");

    results_latex.push_back("\n {\\bf \\large Estimation results for the DIC: }\\\\ \n");
    results_latex.push_back("{\\bf DIC based on the unstandardized deviance } \n");
    results_latex.push_back("\\vspace{-0.4cm}");
    results_latex.push_back("\\begin{tabbing}");
    results_latex.push_back("\\hspace{3cm} \\= \\\\");


    optionsp->out("  DIC based on the unstandardized deviance\n");
    optionsp->out("\n");
    devhelpm = deviance.mean(0);
    optionsp->out("  Deviance(bar_mu):           " +
    ST::doubletostring(deviance2,d) + "\n");
    results_latex.push_back( "deviance($\\bar{\\mu}$) \\> "
                           + ST::doubletostring(deviance2,d)
                           + " \\\\");

    optionsp->out("  pD:                         " +
    ST::doubletostring(devhelpm-deviance2,d) + "\n");
    results_latex.push_back("pD  \\> " +
                            ST::doubletostring(devhelpm-deviance2,d) + " \\\\");


    optionsp->out("  DIC:                        " +
    ST::doubletostring(2*devhelpm-deviance2,d) + "\n");
    optionsp->out("\n");
    results_latex.push_back("DIC  \\> " +
                            ST::doubletostring(2*devhelpm-deviance2,d)
                            + " \\\\");

    results_latex.push_back("\\end{tabbing}\n");
    if(family!="cox")
    {
    optionsp->out("  DIC based on the saturated deviance\n");
    optionsp->out("\n");

    results_latex.push_back("{\\bf DIC based on the saturated deviance } \n");
    results_latex.push_back("\\vspace{-0.4cm}");
    results_latex.push_back("\\begin{tabbing}");
    results_latex.push_back("\\hspace{3cm} \\= \\\\");
    }

    double devhelpm_sat = deviance.mean(1);
    if(family!="cox")
    {
    optionsp->out("  Deviance(bar_mu):           " +
    ST::doubletostring(deviance2_sat,d) + "\n");
    results_latex.push_back( "deviance($\\bar{\\mu}$) \\> "
                           + ST::doubletostring(deviance2_sat,d)
                           + " \\\\");

    optionsp->out("  pD:                         " +
    ST::doubletostring(devhelpm_sat-deviance2_sat,d) + "\n");
    results_latex.push_back("pD \\> " +
                            ST::doubletostring(devhelpm_sat-deviance2_sat,d) + " \\\\");


    optionsp->out("  DIC:                        " +
    ST::doubletostring(2*devhelpm_sat-deviance2_sat,d) + "\n");
    results_latex.push_back("DIC \\> " +
                            ST::doubletostring(2*devhelpm_sat-deviance2_sat,d)
                            + " \\\\");

    optionsp->out("\n");

    results_latex.push_back("\\end{tabbing}\n");
    }
    if(family!="cox")
    out2 << "intnr  unstandardized_deviance  saturated_deviance" << endl;
    else out2 << "intnr  unstandardized_deviance" << endl;
    double * workdev = deviance.getV();
    for (i=0;i<deviance.rows();i++,workdev++)
      {
      out2 << (i+1) << "   ";
      out2 << *workdev << "   ";
      workdev++;
      if (family!="cox") out2 << *workdev;
      out2 << endl;
      }
    if(family!="cox")    out2 << "p_D   " << (devhelpm-deviance2) << "   " <<
                          (devhelpm_sat-deviance2_sat) << endl;
    else                 out2 << "p_D   " << (devhelpm-deviance2) << "   " <<endl;
    if(family!="cox")    out2 << "DIC   " << (2*devhelpm-deviance2) << "   " <<
                          (2*devhelpm_sat-deviance2_sat) << endl;
    else                 out2 << "DIC   " << (2*devhelpm-deviance2) << "   " << endl;


    } // end: if ( (predict) && (optionsp->get_samplesize() > 0) )

  else if ( (predict) && (optionsp->get_samplesize() > 0) && (isbootstrap) )
    {

    register unsigned i,j;

    optionsp->out("  Predicted values:\n",true);
    optionsp->out("\n");
    optionsp->out("  Estimated mean of predictors, expectation of response and\n");
    optionsp->out("  individual deviances are stored in file\n");
    optionsp->out("  " + predictpath + "\n");
    optionsp->out("\n");

    ofstream out(predictpath.strtochar());
    ofstream out2(deviancepath.strtochar());

    for(i=0;i<Dnames.size();i++)
      out << Dnames[i] << "   ";

    unsigned size1 = linearpred.cols();
    unsigned size2 = mumean.cols();

    if (size1 > 1)
      {
      for (i=0;i<size1;i++)
        {
        out << "linpred_" << (i+1) << "   ";
        out << "average_linpred_" << (i+1) << "   ";
        }
      }
    else
      {
      out << "linpred" << "   ";
      out << "average_linpred" << "   ";
      }

    if (size2 > 1)
      {
      for (i=0;i<size2;i++)
        {
        out << "mu_" << (i+1) << "   ";
        out << "average_mu_" << (i+1) << "   ";
        }
      }
    else
      {
      out << "mu" << "   ";
      out << "average_mu" << "   ";
      }

    out << "saturated_deviance" << "   ";

    out << endl;

    double * workmeanave = linpredmean.getV();
    double * workmean = linearpred.getV();
    double * workmuave = mumean.getV();
    double * workdevmean = deviancemean.getV();

    double deviance2=0;
    double deviance2_sat=0;

    double * datap = Dp->getV();
    unsigned sD = Dp->cols();

    // used for computing DIC
    datamatrix scalehelp;
    if (scaleexisting)
      {
      scalehelp = Scalesave.get_betamean();

      for(i=0;i<scalehelp.rows();i++)
        for (j=0;j<scalehelp.cols();j++)
          scalehelp(i,j) = (1.0/(trmult(i,0)*trmult(j,0)))*scalehelp(i,j);
      }

    double reshelp;
    double devhelp;
    double * workmean_firstcol;
    datamatrix mu_meanlinpred(1,size2);

    for (i=0;i<nrobs;i++,workdevmean++)
      {

      for(j=0;j<sD;j++,datap++)
        out << *datap << "   ";

      workmean_firstcol = workmean;
      for (j=0;j<size1;j++,workmean++,workmeanave++)
        {
        *workmean = trmult(j,0)* (*workmean) + addinterceptsample;
        out << *workmean << "   ";
        out << *workmeanave << "   ";
        }

      compute_mu_notransform(workmean_firstcol,&mu_meanlinpred(0,0));

      for (j=0;j<size2;j++,workmuave++)
        {
        out << mu_meanlinpred(0,j) << "   ";
        out << *workmuave << "   ";
        }

      compute_deviance(&response(i,0),&weight(i,0),&mu_meanlinpred(0,0),
      &reshelp,&devhelp,scalehelp,i);

//      compute_deviance(&response(i,0),&weight(i,0),&mumean(i,0),
//      &reshelp,&devhelp,scalehelp,i);

      deviance2 += reshelp;
      deviance2_sat += devhelp;

      if (weight(i,0) != 0)
        {
        out << devhelp  << "   ";
        }
      else
        {
        out << ".   ";
        }

      out << endl;
      }
    } // end: if ( (predict) && stepwise)

  else if ( (predict) && (optionsp->get_samplesize() == 0) )
    {

    register unsigned i,j;

    optionsp->out("  Predicted values:\n",true);
    optionsp->out("\n");
    optionsp->out("  Estimated mean of predictors, expectation of response and\n");
    optionsp->out("  individual deviances are stored in file\n");
    optionsp->out("  " + predictpath + "\n");
    optionsp->out("\n");

    ofstream out(predictpath.strtochar());


    for(i=0;i<Dnames.size();i++)
      out << Dnames[i] << "   ";

    unsigned size1 = linearpred.cols();
    unsigned size2 = mumean.cols();

    if (size1 > 1)
      {
      for (i=0;i<size1;i++)
        {
        out << "linpred_" << (i+1) << "   ";
        }
      }
    else
      {
      out << "linpred" << "   ";
      }

    if (size2 > 1)
      {
      for (i=0;i<size2;i++)
        {
        out << "mu_" << (i+1) << "   ";
        }
      }
    else
      {
      out << "mu" << "   ";
      }

    out << "saturated_deviance" << "   ";

    out << endl;


    double * worklinpred = (*linpred_current).getV();


    double * datap = Dp->getV();
    unsigned sD = Dp->cols();

    double reshelp;
    double devhelp;
    double * workmean_firstcol;
    datamatrix mu_meanlinpred(1,size2);

    double deviance_sat=0;

    for (i=0;i<nrobs;i++)
      {

      for(j=0;j<sD;j++,datap++)
        out << *datap << "   ";

      workmean_firstcol = worklinpred;
      for (j=0;j<size1;j++,worklinpred++)
        {
        *worklinpred = trmult(j,0)* (*worklinpred) + addinterceptsample;
        out << *worklinpred << "   ";
        }

      compute_mu_notransform(workmean_firstcol,&mu_meanlinpred(0,0));

      for (j=0;j<size2;j++)
        {
        out << mu_meanlinpred(0,j) << "   ";
        }

      compute_deviance(&response(i,0),&weight(i,0),&mu_meanlinpred(0,0),
      &reshelp,&devhelp,scale,i);

      deviance_sat += devhelp;

      if (weight(i,0) != 0)
        {
        out << devhelp  << "   ";

        }
      else
        out << ".  ";

      out << endl;
      }

    optionsp->out("\n");
    optionsp->out("  Saturated deviance: " + ST::doubletostring(deviance_sat,6) +
    "\n");

    optionsp->out("\n");

    }


  if (scaleexisting)
    {

    if ( (scale.rows() == 1) && (scale.cols() == 1) )
      {

      double help;

      optionsp->out("  Estimation results for the scale parameter:\n",true);
      optionsp->out("\n");

      results_latex.push_back(
      "\n {\\bf \\large Estimation results for the scale parameter: }\\\\ \n");
      results_latex.push_back("\\vspace{-0.4cm}");
      results_latex.push_back("\\begin{tabbing}");
      results_latex.push_back("\\hspace{3cm} \\= \\\\");


      if (optionsp->get_samplesize() > 0)
        {

        optionsp->out("  Acceptance rate:   " +
        ST::doubletostring(acceptancescale,4) + " %\n");
        optionsp->out("\n");

        ofstream outscale(pathresultsscale.strtochar());

        outscale << "pmean   pstddev   pqu" << nl1 << "   pqu" << nl2 <<
                    "   pqu50   pqu" << nu1 << "   pqu" << nu2 << endl;

        ST::string vstr;

        vstr = "  Mean:         ";
        help = Scalesave.get_betamean(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back("Mean  \\> " +
                            ST::doubletostring(help,6) + " \\\\");


        vstr = "  Std. dev.:    ";
        if (Scalesave.get_betavar(0,0) < 0)
          help = 0;
        else
          help = sqrt(Scalesave.get_betavar(0,0));
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back("Std. dev.:  \\> " +
                            ST::doubletostring(help,6) + " \\\\");

        vstr = "  " + l1 + "% Quantile: ";
        help = Scalesave.get_beta_lower1(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back(vstr.insert_string_char(hchar,helpstring) + " \\>" +
                            ST::doubletostring(help,6) + " \\\\");


        vstr = "  " + l2 + "% Quantile: ";
        help = Scalesave.get_beta_lower2(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back(vstr.insert_string_char(hchar,helpstring) + " \\>" +
                            ST::doubletostring(help,6) + " \\\\");

        vstr = "  50% Quantile: ";
        help = Scalesave.get_betaqu50(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back(vstr.insert_string_char(hchar,helpstring) + " \\>" +
                            ST::doubletostring(help,6) + " \\\\");

        vstr = "  " + u1 + "% Quantile: ";
        help = Scalesave.get_beta_upper2(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back(vstr.insert_string_char(hchar,helpstring) + " \\>" +
                            ST::doubletostring(help,6) + " \\\\");

        vstr = "  " + u2 + "% Quantile: ";
        help = Scalesave.get_beta_upper1(0,0);
        optionsp->out(vstr + ST::string(' ',20-vstr.length()) +
        ST::doubletostring(help,6) + "\n");
        outscale << help << "  ";
        results_latex.push_back(vstr.insert_string_char(hchar,helpstring) + " \\>" +
                            ST::doubletostring(help,6) + " \\\\");

        outscale << endl;

        optionsp->out("\n");

        }
      else  // posteriormode
        {

        ofstream outscale(pathresultsscale.strtochar());

        outscale << "sigma2  "<< endl;


        double m = pow(trmult(0,0),2)*scale(0,0);
        optionsp->out("  sigma2:         " + ST::doubletostring(m,6) + "\n");
        optionsp->out("\n");
        outscale << m << endl;

        results_latex.push_back("$\\sigma^2$:  \\> " +
                            ST::doubletostring(m,6) + " \\\\");

        }

      results_latex.push_back("\\end{tabbing}\n");

      }  // end: if ( (scale.rows() == 1) && (scale.cols() == 1) )
    else  // scale is a matrix
      {
      if (optionsp->get_samplesize() > 0)
        {

        ofstream outscale(pathresultsscale.strtochar());

        unsigned i,j;

        outscale << "row   col   pmean   pstddev   pqu" << nl1 << "   pqu" <<
                     nl2 << "   pqu50   pqu" << nu1 << "   pqu" << nu2 <<
                     "   pmean_corr" <<  endl;

        datamatrix corr(optionsp->get_samplesize(),scale.cols()*scale.cols(),1);

        datamatrix corr_i(scale.cols(),scale.cols(),1);

        unsigned l,r;
        for (i=0;i<optionsp->get_samplesize();i++)
          {
          Scalesave.readsample2(corr_i,i);
          j = 0;
          for (l=0;l<corr_i.rows();l++)
            for(r=0;r<corr_i.cols();r++)
              {
              corr(i,j) = corr_i(l,r)/sqrt(corr_i(l,l)*corr_i(r,r));
              j++;
              }

          }

        l=0;
        for (i=0;i<scale.rows();i++)
          {

          for(j=0;j<scale.cols();j++)
            {
            outscale << (i+1) << "   " << (j+1) << "   ";
            outscale << Scalesave.get_betamean(i,j) << "   ";
            outscale << sqrt(Scalesave.get_betavar(i,j)) << "   ";

            outscale << Scalesave.get_beta_lower1(i,j) <<  "   ";
            outscale << Scalesave.get_beta_lower2(i,j) <<  "   ";
            outscale << Scalesave.get_betaqu50(i,j) << "   ";
            outscale << Scalesave.get_beta_upper1(i,j) <<  "   ";
            outscale << Scalesave.get_beta_upper2(i,j) <<  "   ";
            outscale << corr.mean(l) << "   ";
            l++;
            outscale << endl;
            }

          }

        optionsp->out("  Estimation results for the scale parameter\n",true);
        optionsp->out("\n");

        ST::string help;

        optionsp->out("  Posterior mean:\n");
        optionsp->out("\n");

        for (i=0;i<scale.rows();i++)
          {
          help="  ";
          for(j=0;j<scale.cols();j++)
            {
            help = help + ST::doubletostring(Scalesave.get_betamean(i,j),6)
                   + "   ";
            }
          optionsp->out(help + "\n");
          }

        optionsp->out("\n");
        optionsp->out("  Posterior mean (correlations):\n");
        optionsp->out("\n");

        l=0;
        for (i=0;i<scale.rows();i++)
          {
          help="  ";
          for(j=0;j<scale.cols();j++)
            {
            help = help + ST::doubletostring(corr.mean(l),6)
                   + "   ";
            l++;
            }
          optionsp->out(help + "\n");
          }


        optionsp->out("\n");
        optionsp->out("  Posterior standard deviation:\n");
        optionsp->out("\n");

        for (i=0;i<scale.rows();i++)
          {
          help="  ";
          for(j=0;j<scale.cols();j++)
            {
            help = help + ST::doubletostring(sqrt(Scalesave.get_betavar(i,j)),6)
                   + "   ";
            }
          optionsp->out(help + "\n");
          }

        optionsp->out("\n");
        optionsp->out("  Posterior " + l1 + " percent quantile:\n");
        optionsp->out("\n");

        for (i=0;i<scale.rows();i++)
          {
          help="  ";
          for(j=0;j<scale.cols();j++)
            {
            help = help + ST::doubletostring(Scalesave.get_beta_lower1(i,j),6)
                   + "   ";
            }
          optionsp->out(help + "\n");
          }

        optionsp->out("\n");
        optionsp->out("  Posterior 50\% quantile:\n");
        optionsp->out("\n");

        for (i=0;i<scale.rows();i++)
          {
          help="  ";
          for(j=0;j<scale.cols();j++)
            {
            help = help + ST::doubletostring(Scalesave.get_betaqu50(i,j),6)
                   + "   ";
            }
          optionsp->out(help + "\n");
          }

        optionsp->out("\n");
        optionsp->out("  Posterior " + u1 + " percent quantile:\n");
        optionsp->out("\n");

        for (i=0;i<scale.rows();i++)
          {
          help="  ";
          for(j=0;j<scale.cols();j++)
            {
            help = help + ST::doubletostring(Scalesave.get_beta_upper1(i,j),6)
                   + "   ";
            }
          optionsp->out(help + "\n");
          }

        optionsp->out("\n");
        optionsp->out("\n");

        } // end: if samplesize > 0
      else // posterior mode
        {
        ofstream outscale(pathresultsscale.strtochar());

        unsigned i,j;

        outscale << "row   col   pmean" << endl;

        for (i=0;i<scale.rows();i++)
          {
          for(j=0;j<scale.cols();j++)
            {
            outscale << (i+1) << "   " << (j+1) << "   ";
            outscale << trmult(i,0)*trmult(j,0)*scale(i,j) << "   ";
            outscale << endl;
            }
          }

        optionsp->out("  Estimation results for the scale parameter\n",true);
        optionsp->out("\n");

        ST::string help;

        optionsp->out("  Posterior mean:\n");
        optionsp->out("\n");

        for (i=0;i<scale.rows();i++)
          {
          help="  ";
          for(j=0;j<scale.cols();j++)
            {
            help = help + ST::doubletostring(trmult(i,0)*trmult(j,0)*scale(i,j),6)
                   + "   ";
            }
          optionsp->out(help + "\n");
          }

        optionsp->out("\n");
        optionsp->out("\n");

        }

      }

    } // end: if (scaleexisting)

  optionsp->out("\n");

  }


void DISTRIBUTION::reset(void)
  {
  linearpred = datamatrix(nrobs,linearpred.cols(),0);
  linearpredprop = linearpred;
  linpred_current = &linearpred;
  linpred_proposed = &linearpredprop;
  if (scaleexisting)
    {
    scale = datamatrix(1,1,0.1);
    }
  Scalesave.reset();
  if (predictfull)
    musave.reset();
  }

void DISTRIBUTION::posteriormode_set_beta_mode(void)
  {
  scale_mode.assign(scale);
  }


bool DISTRIBUTION::posteriormode_converged_fc(const datamatrix & beta,
                                              const datamatrix & beta_mode,
                                              const unsigned & itnr)
  {
  if (itnr > 1)
    {
    double diff;
    double normold;

    normold = norm(beta_mode);

    diff = norm(beta-beta_mode)/normold;

    if (diff <= 0.00001)
      return true;
    else
      return false;
    }
  else
    return false;
  }


void DISTRIBUTION::update_missings(void)
  {
  }


void DISTRIBUTION::set_missings(vector<FULLCOND *> & fcm,
                                unsigned & begin,unsigned & end,
                                datamatrix & mi,ST::string & pt,
                                ST::string & pr)
  {
  unsigned i;
  bool found=false;
  for(i=begin;i<=end;i++)
    {
    if (fcm[i]->is_missing(responsename))
      {
      found=true;
      fcmissing.push_back(fcm[i]);
      }
    }

  if (found==true)
    {

    pathmissing = pr;

    missingind = mi;

    double * workm = missingind.getV();
    unsigned nrmissing=0;

    for (i=0;i<nrobs;i++,workm++)
      {
      if (*workm == 0)
        nrmissing ++;
      }

    missingpos = statmatrix<unsigned>(nrmissing,1);

    workm = missingind.getV();
    unsigned start = 0;
    unsigned im = 0;

    for (i=0;i<nrobs;i++,workm++)
      {
      if (*workm==0)
        {
        missingpos(im,0) = i-start;
        im++;
        start = i;
        }

      }


    MissingSave = FULLCOND(optionsp,datamatrix(1,1),
                  "Missingvalues",nrmissing,1,pt);

    MissingSave.set_transform(trmult(0,0));

    MissingSave.setflags(MCMC::norelchange | MCMC::nooutput);



    }  // end: if (found==true)

  }

//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_gamma -----------------------------
//------------------------------------------------------------------------------

double DISTRIBUTION_gamma::compute_weight(double * worklin,double * weight,
                                          const int & i,const unsigned & col)
                                          const
  {
  return 1/scale(0,0) * *weight;
  }


void DISTRIBUTION_gamma::create_lgamma(void)
  {
  lgamma = datamatrix(101,1);

  lgamma(0,0) = 0.000000000;
  lgamma(1,0) = -0.005690308;
  lgamma(2,0) = -0.011218489;
  lgamma(3,0) = -0.016586854;
  lgamma(4,0) = -0.021797651;
  lgamma(5,0) = -0.026853073;
  lgamma(6,0) = -0.031755254;
  lgamma(7,0) = -0.036506276;
  lgamma(8,0) = -0.041108170;
  lgamma(9,0) = -0.045562915;
  lgamma(10,0) = -0.049872441;
  lgamma(11,0) = -0.054038634;
  lgamma(12,0) = -0.058063333;
  lgamma(13,0) = -0.061948332;
  lgamma(14,0) = -0.065695387;
  lgamma(15,0) = -0.069306209;
  lgamma(16,0) = -0.072782472;
  lgamma(17,0) = -0.076125811;
  lgamma(18,0) = -0.079337824;
  lgamma(19,0) = -0.082420074;
  lgamma(20,0) = -0.085374090;
  lgamma(21,0) = -0.088201365;
  lgamma(22,0) = -0.090903362;
  lgamma(23,0) = -0.093481511;
  lgamma(24,0) = -0.095937212;
  lgamma(25,0) = -0.098271836;
  lgamma(26,0) = -0.100486725;
  lgamma(27,0) = -0.102583193;
  lgamma(28,0) = -0.104562527;
  lgamma(29,0) = -0.106425987;
  lgamma(30,0) = -0.108174810;
  lgamma(31,0) = -0.109810204;
  lgamma(32,0) = -0.111333359;
  lgamma(33,0) = -0.112745436;
  lgamma(34,0) = -0.114047576;
  lgamma(35,0) = -0.115240897;
  lgamma(36,0) = -0.116326498;
  lgamma(37,0) = -0.117305454;
  lgamma(38,0) = -0.118178821;
  lgamma(39,0) = -0.118947635;
  lgamma(40,0) = -0.119612914;
  lgamma(41,0) = -0.120175656;
  lgamma(42,0) = -0.120636841;
  lgamma(43,0) = -0.120997431;
  lgamma(44,0) = -0.121258371;
  lgamma(45,0) = -0.121420591;
  lgamma(46,0) = -0.121485001;
  lgamma(47,0) = -0.121452498;
  lgamma(48,0) = -0.121323962;
  lgamma(49,0) = -0.121100259;
  lgamma(50,0) = -0.120782238;
  lgamma(51,0) = -0.120370735;
  lgamma(52,0) = -0.119866573;
  lgamma(53,0) = -0.119270560;
  lgamma(54,0) = -0.118583490;
  lgamma(55,0) = -0.117806145;
  lgamma(56,0) = -0.116939293;
  lgamma(57,0) = -0.115983691;
  lgamma(58,0) = -0.114940083;
  lgamma(59,0) = -0.113809201;
  lgamma(60,0) = -0.112591766;
  lgamma(61,0) = -0.111288486;
  lgamma(62,0) = -0.109900061;
  lgamma(63,0) = -0.108427177;
  lgamma(64,0) = -0.106870510;
  lgamma(65,0) = -0.105230728;
  lgamma(66,0) = -0.103508486;
  lgamma(67,0) = -0.101704430;
  lgamma(68,0) = -0.099819197;
  lgamma(69,0) = -0.097853413;
  lgamma(70,0) = -0.095807697;
  lgamma(71,0) = -0.093682657;
  lgamma(72,0) = -0.091478893;
  lgamma(73,0) = -0.089196995;
  lgamma(74,0) = -0.086837546;
  lgamma(75,0) = -0.084401121;
  lgamma(76,0) = -0.081888285;
  lgamma(77,0) = -0.079299595;
  lgamma(78,0) = -0.076635603;
  lgamma(79,0) = -0.073896851;
  lgamma(80,0) = -0.071083873;
  lgamma(81,0) = -0.068197197;
  lgamma(82,0) = -0.065237343;
  lgamma(83,0) = -0.062204825;
  lgamma(84,0) = -0.059100148;
  lgamma(85,0) = -0.055923813;
  lgamma(86,0) = -0.052676312;
  lgamma(87,0) = -0.049358131;
  lgamma(88,0) = -0.045969750;
  lgamma(89,0) = -0.042511642;
  lgamma(90,0) = -0.038984276;
  lgamma(91,0) = -0.035388112;
  lgamma(92,0) = -0.031723605;
  lgamma(93,0) = -0.027991206;
  lgamma(94,0) = -0.024191358;
  lgamma(95,0) = -0.020324499;
  lgamma(96,0) = -0.016391062;
  lgamma(97,0) = -0.012391474;
  lgamma(98,0) = -0.008326158;
  lgamma(99,0) = -0.004195529;
  lgamma(100,0) = 0.000000000;

  }


double DISTRIBUTION_gamma::lgammafunc(const double & nu) const
    {

    if (fmod(nu,1)==0)
      return lfac(nu-1);
    else if (nu<1)
      return lgammafunc(nu+1) - log(nu);
    else if (nu>2)
      return log(nu-1) + lgammafunc(nu-1);
    else
      return lgamma(int(nu*100)-100,0);
    }


double DISTRIBUTION_gamma::lfac(const double & nu) const
    {
    if (nu==0 || nu==1) return 0;
    else return log(nu) + lfac(nu-1);
    }


bool DISTRIBUTION_gamma::posteriormode(void)
  {

  if (!scalefixed)
    scale(0,0) = phi_hat();

  return true;

  }


bool DISTRIBUTION_gamma::posteriormode_converged(const unsigned & itnr)
  {

  return true;

  }


void DISTRIBUTION_gamma::update(void)
  {


  if (scalefixed)                     // with constant scale
    {

    }
  else if (!mh)                     // with consistent estimation for scale
    {
    nriterations++;

    if (nriterations>const_it)
      {
      scale(0,0) = phi_hat();
      }

    }

  else if (mh)                             // with MH-algorithm for scale
    {

    nriterations++;

    if (nriterations <= const_it)
      {
      }
    else if ( (nriterations > const_it) &&
              (nriterations <= optionsp->get_burnin()) )
      {
      scale(0,0) = phi_hat();
      }

    else if (nriterations > optionsp->get_burnin())
      {


      // random number, to be drawn from the proposal of nu: gamma(a_nu,b_nu)
      double nu_prop = 0;
      // log acceptance
      double logalpha = 0;


      double var = var_nu;



      register unsigned i;
      double * worklin = (*linpred_current).getV();;
      double * workresp = response.getV();
      double * workweight = weight.getV();
      double sumresp = 0;

      // new hyperparameter b for the full conditional of
      // the shape parameter nu
      double bnew_gamma = b_gamma;

      double sumweight = 0;

      for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
        {
        if (*workweight != 0)
          {
          bnew_gamma += *worklin - log(*workresp) + (*workresp / exp(*worklin)) ;
          sumresp -= log(*workresp);
          sumweight+= *workweight;
          }
        }

      double nu = 1/scale(0,0);
      double a_nu,b_nu;
      a_nu = nu*nu / var;
      b_nu = a_nu / nu;

      nu_prop = rand_gamma(a_nu,b_nu);

      logalpha = log_fullcond(nu_prop,bnew_gamma,sumweight) -
                 log_fullcond(nu,bnew_gamma,sumweight)      +
                 (a_nu - 1)*(log(nu) - log(nu_prop)) - b_nu*(nu-nu_prop);

      double u = log(uniform());

      if (nriterations==optionsp->get_burnin()+1)
        {
        acceptance = 0;
        }

      if (u<=logalpha)
        {
        nu = nu_prop;
        scale(0,0) = 1/nu;
        acceptance++;
        }
      } // else if



    }  //else if (MH)

  DISTRIBUTION::update();

  }


void DISTRIBUTION_gamma::check(void)
  {

  unsigned i=0;
  bool error=false;
  double * workr=response.getV();
  while ( (i<nrobs) && (error==false) )
    {
    if (*workr <= 0)
      {
      error=true;
      errors.push_back(
    "ERROR: response cannot be gamma distributed; some values are negative\n");
      }
    i++;
    workr++;
    }

  }



void DISTRIBUTION_gamma::standardize(void)
  {

  double s = sqrt(response.var(0,weight));
  unsigned i;
  double * workresp = response.getV();
  for (i=0;i<nrobs;i++,workresp++)
   {
   *workresp = *workresp/s;
   }

  addinterceptsample = log(s);

  }


// CONSTRUCTOR 0

DISTRIBUTION_gamma::DISTRIBUTION_gamma(const double & scale_initial,
                                       MCMCoptions * o,
                                       const datamatrix & r,
                                       const ST::string & p,
                                       const ST::string & ps,
                                       const datamatrix & w)
  : DISTRIBUTION(o,r,w,p,ps)
  {

  create_lgamma();

  nriterations = 0;
  acceptance = 0;
  optionsp = o;

  a_gamma = 0;
  b_gamma = 0;
  const_it = 0;

  mh = false;

  scaleold = 0;

  scalefixed=true;
  scaleexisting=true;
  scale(0,0) = scale_initial;

  family = "Gamma";

  check();

  standardize();

  }



// CONSTRUCTOR 1

DISTRIBUTION_gamma::DISTRIBUTION_gamma(const double & a, const double & b,
                                       const unsigned & cit, MCMCoptions * o,
                                       const datamatrix & r,
                                       const ST::string & p,
                                       const ST::string & ps,
                                       const datamatrix & w)
  : DISTRIBUTION(o,r,w,p,ps)
  {

  assert (a > 0);
  assert (b > 0);

  create_lgamma();

  mh = false;

  nriterations = 0;
  acceptance = 0;
  optionsp = o;

  a_gamma = a;
  b_gamma = b;
  const_it = cit;
  if (o->get_burnin() < const_it)
    const_it = o->get_burnin();

  var_nu = 0;
  scaleold = 0;
  family = "Gamma";

  check();

  scalefixed=false;

  standardize();

  }


// CONSTRUCTOR 2

DISTRIBUTION_gamma::DISTRIBUTION_gamma(const double & a, const double & b,
                                       const double & var, const unsigned & cit,
                                       MCMCoptions * o, const datamatrix & r,
                                       const ST::string & p,
                                       const ST::string & ps,
                                       const datamatrix & w)
  : DISTRIBUTION(o,r,w,p,ps)
  {

  assert (a > 0);
  assert (b > 0);

  mh = true;

  create_lgamma();

  nriterations = 0;
  acceptance = 0;
  optionsp = o;

  a_gamma = a;
  b_gamma = b;
  var_nu = var;
  const_it = cit;
  if (o->get_burnin() < const_it)
    const_it = o->get_burnin();

  scaleold = 0;

  family = "Gamma";

  check();

  scalefixed=false;

  standardize();

  }


DISTRIBUTION_gamma::DISTRIBUTION_gamma(const DISTRIBUTION_gamma & ga)
   : DISTRIBUTION(DISTRIBUTION(ga))
   {
   a_gamma = ga.a_gamma;
   b_gamma = ga.b_gamma;
   scale = ga.scale;
   scaleold = ga.scaleold;
   mh = ga.mh;
   const_it = ga.const_it;
   var_nu = ga.var_nu;
   lgamma=ga.lgamma;
   nriterations=ga.nriterations;
   acceptance=ga.acceptance;
   scalefixed = ga.scalefixed;
   }


const DISTRIBUTION_gamma & DISTRIBUTION_gamma::operator=(
                                      const DISTRIBUTION_gamma & ga)
  {

  if (this==&ga)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(ga));
  a_gamma = ga.a_gamma;
  b_gamma = ga.b_gamma;
  scale = ga.scale;
  scaleold = ga.scaleold;
  mh = ga.mh;
  const_it = ga.const_it;
  var_nu = ga.var_nu;
  lgamma=ga.lgamma;
  nriterations=ga.nriterations;
  acceptance=ga.acceptance;
  scalefixed = ga.scalefixed;
  return *this;
  }


double DISTRIBUTION_gamma::loglikelihood(double * res, double * lin,
                                         double * w,const int & i) const
  {

  if (*w != 0)
    return  - *w * (*res / exp(*lin) + *lin) / (scale(0,0));
  else
    return 0;

  }


void DISTRIBUTION_gamma::compute_mu(const double * linpred,double * mu) const
  {
  *mu = exp(*linpred);
  }

void DISTRIBUTION_gamma::compute_mu_notransform(
const double * linpred,double * mu) const
  {
  *mu = exp(*linpred);
  }


double DISTRIBUTION_gamma::compute_gmu(double * linpred,
                                       const unsigned & col) const
  {
  return 1.0/exp(*linpred);
  }


double DISTRIBUTION_gamma::compute_IWLS(double * response,
                       double * linpred, double * weight,
                       const int & i,double * weightiwls,double * tildey,
                       bool weightyes, const unsigned & col)
  {

  double m = exp(*linpred);

  if (weightyes)
    *weightiwls = 1/scale(0,0) * *weight;

  *tildey = (*response - m)/m;

  if (*weight != 0)
    return  - *weight * (*response / m + *linpred) / (scale(0,0));
  else
    return 0;

  }


void DISTRIBUTION_gamma::compute_IWLS_weight_tildey(
                               double * response,double * linpred,
                               double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col)
  {

  double m = exp(*linpred);

  *weightiwls = 1/scale(0,0) * *weight;

  *tildey = (*response - m)/m;

  }


// BEGIN: DSB //
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
double
    DISTRIBUTION_gamma::sample_from_likelihood(const double weight,
                                               const double mu) const
    {
        // get the current "scale" parameter of the Gamma likelihood
        double nu = (weight != 0.0) ? weight / scale(0, 0) : 1.0 / scale(0, 0);
        // todo: is this correct?? Even if the observation has weight 0,
        // we should be able to produce samples from the predictive distribution
        // which conditions on the observations with positive weights!

        // transform into the Gamma(a, b) parametrization
        double a = nu;
        double b = nu / mu;

        // and then use the sample function from Random.cpp
        return randnumbers::rand_gamma(a, b);
    }
#endif
// END: DSB //


void DISTRIBUTION_gamma::compute_deviance(const double * response,
                                          const double * weight,
                                          const double * mu,double * deviance,
                        double * deviancesat,const datamatrix & scale,const int & i) const
  {
  if (*weight!=0)
    {

    double rtr = *response*exp(addinterceptsample);
    double nu   = (*weight)/scale(0,0);
    double help = nu/(*mu);
    double help2 = (rtr)/(*mu);
    if (scalefixed)
      {
      *deviance   = 2*(nu*log(*mu)
                  -(nu-1)*log(rtr)+ help*(rtr));
      *deviancesat = 2*nu*(-log(help2) + help2 - 1);
      }
    else
      {
      *deviance   = 2*(lgammafunc(nu)-nu*log(help)
                  -(nu-1)*log(rtr)+ help*(rtr));
      *deviancesat = 2*nu*(-log(help2) + help2 - 1);
      }
    }
  else
    {
    *deviance=0;
    *deviancesat=0;
    }

  }


void DISTRIBUTION_gamma::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("  Response function: exp\n");

  if (mh)
    {
    optionsp->out("  Update of scale parameter by MH-algorithm \n");
    optionsp->out("  Fixed variance: " + ST::doubletostring(var_nu,6) + "\n");
    optionsp->out("  Hyperparameter a: " + ST::doubletostring(a_gamma,6) + "\n");
    optionsp->out("  Hyperparameter b: " + ST::doubletostring(b_gamma,6) + "\n");
    }

  else if (!mh && !scalefixed)
    {
    optionsp->out("  Update of scale parameter by consistent estimation \n");
    }

  else
    {
    optionsp->out("  Fixed scale parameter: " + ST::doubletostring(scale(0,0),6) + "\n");
    }

  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTRIBUTION_gamma::update_predict(void)
  {
    if(
      (optionsp->get_nriter() > optionsp->get_burnin())
      &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1)
       % (optionsp->get_step()) == 0)
      )
      {
      add_linearpred(addinterceptsample);
      }

  DISTRIBUTION::update_predict();

  if(
    (optionsp->get_nriter() > optionsp->get_burnin())
    &&
    ((optionsp->get_nriter()-optionsp->get_burnin()-1)
    % (optionsp->get_step()) == 0)
    )
    {
    add_linearpred(-addinterceptsample);
    }
  }


void DISTRIBUTION_gamma::outresults(void)
    {

    if (mh)
      {
      acceptancescale = double(acceptance)/
                        double(nriterations-optionsp->get_burnin())*100;

      }

    else
      {
      acceptancescale=100;
      }


    DISTRIBUTION::outresults();

    }


double DISTRIBUTION_gamma::phi_hat() const
  {

  register unsigned i;
  double* workweight = weight.getV();
  double* worklin = (*linpred_current).getV();
  double* workres = response.getV();
  double phi = 0;

  double explin;
  double diff;
  double sumweight = 0;
  for (i=0;i<nrobs;i++,workweight++,worklin++,workres++)
     {
     if (*workweight != 0)
       {
       sumweight += *workweight;
       explin = exp(*worklin);
       diff = (*workres - explin);
       phi += diff*diff / (explin*explin / *workweight);
       }
     }

  return phi / sumweight;    // ohne p!

  }


double DISTRIBUTION_gamma::log_fullcond(const double & nu,
                                        double & bnew, double & sw) const
  {
  return - (sw * lgammafunc(nu)) +
           log(nu)*(sw*nu + a_gamma -1) - nu*bnew;
  }


double DISTRIBUTION_gamma::log_prop(const double & nu,
                                    double & a, double & b) const
  {
  return a*log(b) - lgammafunc(a) + (a-1)*log(nu) - b*nu;
  }


void DISTRIBUTION_gamma::tr_nonlinear(vector<double *> b,vector<double *> br,
                    vector<FULLCOND*> & fcp,unsigned & nr,
                    unsigned & it,ST::string & trtype)
  {

  if (trtype == "exp")
    DISTRIBUTION::tr_nonlinear(b,br,fcp,nr,it,trtype);
  else if (trtype == "elasticity")
    {
    if (b.size() == 2)
      {
      *br[1] = *b[1] * fcp[0]->get_data(nr,0);
      }
    }
  else if (trtype=="marginal")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0)+ *b[i]);
      }
    }
  else if (trtype=="marginalintercept")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0));
      }
    }

  }


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_gamma (stepwise) ------------------
//------------------------------------------------------------------------------

void DISTRIBUTION_gamma2::set_constscale(double s)
  {
  scale(0,0) = s/(trmult(0,0)*trmult(0,0));
  scalefixed=true;
  }

void DISTRIBUTION_gamma2::undo_constscale(void)
  {
  scalefixed = false;
  }

double DISTRIBUTION_gamma2::compute_weight(double * worklin,double * weight,
                                          const int & i,const unsigned & col)
                                          const
  {
  //return 1/scale(0,0) * *weight;
  return *weight;
  }


void DISTRIBUTION_gamma2::create_lgamma(void)
  {
  lgamma = datamatrix(101,1);

  lgamma(0,0) = 0.000000000;
  lgamma(1,0) = -0.005690308;
  lgamma(2,0) = -0.011218489;
  lgamma(3,0) = -0.016586854;
  lgamma(4,0) = -0.021797651;
  lgamma(5,0) = -0.026853073;
  lgamma(6,0) = -0.031755254;
  lgamma(7,0) = -0.036506276;
  lgamma(8,0) = -0.041108170;
  lgamma(9,0) = -0.045562915;
  lgamma(10,0) = -0.049872441;
  lgamma(11,0) = -0.054038634;
  lgamma(12,0) = -0.058063333;
  lgamma(13,0) = -0.061948332;
  lgamma(14,0) = -0.065695387;
  lgamma(15,0) = -0.069306209;
  lgamma(16,0) = -0.072782472;
  lgamma(17,0) = -0.076125811;
  lgamma(18,0) = -0.079337824;
  lgamma(19,0) = -0.082420074;
  lgamma(20,0) = -0.085374090;
  lgamma(21,0) = -0.088201365;
  lgamma(22,0) = -0.090903362;
  lgamma(23,0) = -0.093481511;
  lgamma(24,0) = -0.095937212;
  lgamma(25,0) = -0.098271836;
  lgamma(26,0) = -0.100486725;
  lgamma(27,0) = -0.102583193;
  lgamma(28,0) = -0.104562527;
  lgamma(29,0) = -0.106425987;
  lgamma(30,0) = -0.108174810;
  lgamma(31,0) = -0.109810204;
  lgamma(32,0) = -0.111333359;
  lgamma(33,0) = -0.112745436;
  lgamma(34,0) = -0.114047576;
  lgamma(35,0) = -0.115240897;
  lgamma(36,0) = -0.116326498;
  lgamma(37,0) = -0.117305454;
  lgamma(38,0) = -0.118178821;
  lgamma(39,0) = -0.118947635;
  lgamma(40,0) = -0.119612914;
  lgamma(41,0) = -0.120175656;
  lgamma(42,0) = -0.120636841;
  lgamma(43,0) = -0.120997431;
  lgamma(44,0) = -0.121258371;
  lgamma(45,0) = -0.121420591;
  lgamma(46,0) = -0.121485001;
  lgamma(47,0) = -0.121452498;
  lgamma(48,0) = -0.121323962;
  lgamma(49,0) = -0.121100259;
  lgamma(50,0) = -0.120782238;
  lgamma(51,0) = -0.120370735;
  lgamma(52,0) = -0.119866573;
  lgamma(53,0) = -0.119270560;
  lgamma(54,0) = -0.118583490;
  lgamma(55,0) = -0.117806145;
  lgamma(56,0) = -0.116939293;
  lgamma(57,0) = -0.115983691;
  lgamma(58,0) = -0.114940083;
  lgamma(59,0) = -0.113809201;
  lgamma(60,0) = -0.112591766;
  lgamma(61,0) = -0.111288486;
  lgamma(62,0) = -0.109900061;
  lgamma(63,0) = -0.108427177;
  lgamma(64,0) = -0.106870510;
  lgamma(65,0) = -0.105230728;
  lgamma(66,0) = -0.103508486;
  lgamma(67,0) = -0.101704430;
  lgamma(68,0) = -0.099819197;
  lgamma(69,0) = -0.097853413;
  lgamma(70,0) = -0.095807697;
  lgamma(71,0) = -0.093682657;
  lgamma(72,0) = -0.091478893;
  lgamma(73,0) = -0.089196995;
  lgamma(74,0) = -0.086837546;
  lgamma(75,0) = -0.084401121;
  lgamma(76,0) = -0.081888285;
  lgamma(77,0) = -0.079299595;
  lgamma(78,0) = -0.076635603;
  lgamma(79,0) = -0.073896851;
  lgamma(80,0) = -0.071083873;
  lgamma(81,0) = -0.068197197;
  lgamma(82,0) = -0.065237343;
  lgamma(83,0) = -0.062204825;
  lgamma(84,0) = -0.059100148;
  lgamma(85,0) = -0.055923813;
  lgamma(86,0) = -0.052676312;
  lgamma(87,0) = -0.049358131;
  lgamma(88,0) = -0.045969750;
  lgamma(89,0) = -0.042511642;
  lgamma(90,0) = -0.038984276;
  lgamma(91,0) = -0.035388112;
  lgamma(92,0) = -0.031723605;
  lgamma(93,0) = -0.027991206;
  lgamma(94,0) = -0.024191358;
  lgamma(95,0) = -0.020324499;
  lgamma(96,0) = -0.016391062;
  lgamma(97,0) = -0.012391474;
  lgamma(98,0) = -0.008326158;
  lgamma(99,0) = -0.004195529;
  lgamma(100,0) = 0.000000000;

  }


double DISTRIBUTION_gamma2::lgammafunc(const double & nu) const
    {

    if (fmod(nu,1)==0)
      return lfac(nu-1);
    else if (nu<1)
      return lgammafunc(nu+1) - log(nu);
    else if (nu>2)
      return log(nu-1) + lgammafunc(nu-1);
    else
      return lgamma(int(nu*100)-100,0);
    }

double DISTRIBUTION_gamma2::lfac(const double & nu) const
    {
    if (nu==0 || nu==1) return 0;
    else return log(nu) + lfac(nu-1);
    }


bool DISTRIBUTION_gamma2::posteriormode(void)
  {

  if (!scalefixed)
    scale(0,0) = phi_hat();

  return true;

  }


bool DISTRIBUTION_gamma2::posteriormode_converged(const unsigned & itnr)
  {

  return true;

  }


void DISTRIBUTION_gamma2::check(void)
  {

  unsigned i=0;
  bool error=false;
  double * workr=response.getV();
  while ( (i<nrobs) && (error==false) )
    {
    if (*workr <= 0)
      {
      error=true;
      errors.push_back(
    "ERROR: response cannot be gamma distributed; some values are negative\n");
      }
    i++;
    workr++;
    }

  }


void DISTRIBUTION_gamma2::standardize(void)
  {

  double s = sqrt(response.var(0,weight));
  unsigned i;
  double * workresp = response.getV();
  for (i=0;i<nrobs;i++,workresp++)
   {
   *workresp = *workresp/s;
   }

  addinterceptsample = log(s);

  }


// CONSTRUCTOR 0

DISTRIBUTION_gamma2::DISTRIBUTION_gamma2(const double & scale_initial,
                                       MCMCoptions * o,
                                       const datamatrix & r,
                                       const ST::string & p,
                                       const ST::string & ps,
                                       const datamatrix & w)
  : DISTRIBUTION(o,r,w,p,ps)
  {

  create_lgamma();

  optionsp = o;

  scaleold = 0;

  scalefixed=true;
  scaleexisting=true;
  scale(0,0) = scale_initial;

  family = "Gamma";

  check();

  standardize();

  }



// CONSTRUCTOR 1

DISTRIBUTION_gamma2::DISTRIBUTION_gamma2(const double & a, const double & b,
                                       const unsigned & cit, MCMCoptions * o,
                                       const datamatrix & r,
                                       const ST::string & p,
                                       const ST::string & ps,
                                       const datamatrix & w)
  : DISTRIBUTION(o,r,w,p,ps)
  {

  assert (a > 0);
  assert (b > 0);

  create_lgamma();

  optionsp = o;

  scaleold = 0;
  family = "Gamma";

  check();

  scalefixed=false;

  standardize();

  }


DISTRIBUTION_gamma2::DISTRIBUTION_gamma2(const DISTRIBUTION_gamma2 & ga)
   : DISTRIBUTION(DISTRIBUTION(ga))
   {
   scale = ga.scale;
   scaleold = ga.scaleold;
   lgamma=ga.lgamma;
   scalefixed = ga.scalefixed;
   }


const DISTRIBUTION_gamma2 & DISTRIBUTION_gamma2::operator=(
                                      const DISTRIBUTION_gamma2 & ga)
  {

  if (this==&ga)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(ga));
  scale = ga.scale;
  scaleold = ga.scaleold;
  lgamma=ga.lgamma;
  scalefixed = ga.scalefixed;
  return *this;
  }


void DISTRIBUTION_gamma2::compute_mu(const double * linpred,double * mu) const
  {
  *mu = exp(*linpred);
  }

void DISTRIBUTION_gamma2::compute_mu_notransform(
const double * linpred,double * mu) const
  {
  *mu = exp(*linpred);
  }


double DISTRIBUTION_gamma2::compute_gmu(double * linpred,
                                       const unsigned & col) const
  {
  return 1.0/exp(*linpred);
  }


// BEGIN: DSB //
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
// this is just a copy of DISTRIBUTION_gamma::sample_from_likelihood
double
    DISTRIBUTION_gamma2::sample_from_likelihood(const double weight,
                                                const double mu) const
    {
        // get the current "scale" parameter of the Gamma likelihood
        double nu = (weight != 0.0) ? weight / scale(0, 0) : 1.0 / scale(0, 0);
        // todo: is this correct?

        // transform into the Gamma(a, b) parametrization
        double a = nu;
        double b = nu / mu;

        // and then use the sample function from Random.cpp
        return randnumbers::rand_gamma(a, b);
    }
#endif
// END: DSB //

void DISTRIBUTION_gamma2::compute_deviance(const double * response,
                                          const double * weight,
                                          const double * mu,double * deviance,
                        double * deviancesat,const datamatrix & scale,const int & i) const
  {
  if (*weight!=0)
    {

    double rtr = *response;
    double nu   = (*weight)/scale(0,0);
    double help = nu/(*mu);
    double help2 = (rtr)/(*mu);
    if (scalefixed)
      {
      *deviance   = 2*(nu*log(*mu)
                  -(nu-1)*log(rtr)+ help*(rtr));
      *deviancesat = 2* *weight *(-log(help2) + help2 - 1);
      }
    else
      {
      *deviance   = 2*(lgammafunc(nu)-nu*log(help)
                  -(nu-1)*log(rtr)+ help*(rtr));
      *deviancesat = 2* *weight *(-log(help2) + help2 - 1);
      }
    *deviance -= 2*addinterceptsample;
    }
  else
    {
    *deviance=0;
    *deviancesat=0;
    }

  }


void DISTRIBUTION_gamma2::outresults(void)
    {
    DISTRIBUTION::outresults();
    }


double DISTRIBUTION_gamma2::phi_hat() const
  {

  register unsigned i;
  double* workweight = weight.getV();
  double* worklin = (*linpred_current).getV();
  double* workres = response.getV();
  double phi = 0;

  double explin;
  double diff;
  double sumweight = 0;
  for (i=0;i<nrobs;i++,workweight++,worklin++,workres++)
     {
     if (*workweight != 0)
       {
       sumweight += 1; // *workweight;
       explin = exp(*worklin);
       diff = (*workres - explin);
       phi += diff*diff / (explin*explin / *workweight);
       }
     }

  return phi / sumweight;    // ohne p!
  }


double DISTRIBUTION_gamma2::compute_IWLS(double * response,
                       double * linpred, double * weight,
                       const int & i,double * weightiwls,double * tildey,
                       bool weightyes, const unsigned & col)
  {
  double m = exp(*linpred);
  if (weightyes)
    *weightiwls = *weight;
    //*weightiwls = 1/scale(0,0) * *weight;

  *tildey = (*response - m)/m;

  if (*weight != 0)
    return  - *weight * (*response / m + *linpred) / (scale(0,0));
  else
    return 0;
  }


double DISTRIBUTION_gamma2::loglikelihood(double * res, double * lin,
                                         double * w,const int & i) const
  {
  if (*w != 0)
    return  - *w * (*res / exp(*lin) + *lin) / (scale(0,0));
  else
    return 0;
  }


void DISTRIBUTION_gamma2::compute_bootstrap_data(const double * linpred,const double * weight,double * wresp)
  {
  double mu = exp(*linpred);
  double pres = 0;
  if(*weight>0)
    {
    mu *= *weight;
    double nu = *weight/scale(0,0);
    double a = nu;
    double b = nu/mu;
    pres = rand_gamma(a,b);
    pres /= *weight;
    }
  *wresp = pres;
  }


void DISTRIBUTION_gamma2::update_predict(void)
  {
    if(
      (optionsp->get_nriter() > optionsp->get_burnin())
      &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1)
       % (optionsp->get_step()) == 0)
      )
      {
      add_linearpred(addinterceptsample);
      }

  DISTRIBUTION::update_predict();

  if(
    (optionsp->get_nriter() > optionsp->get_burnin())
    &&
    ((optionsp->get_nriter()-optionsp->get_burnin()-1)
    % (optionsp->get_step()) == 0)
    )
    {
    add_linearpred(-addinterceptsample);
    }
  }


void DISTRIBUTION_gamma2::compute_IWLS_weight_tildey(
                               double * response,double * linpred,
                               double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col)
  {
  double m = exp(*linpred);
  //*weightiwls = 1/scale(0,0) * *weight;
  *weightiwls = *weight;
  *tildey = (*response - m)/m;
  }

//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_vargaussian -----------------------
//------------------------------------------------------------------------------


bool DISTRIBUTION_vargaussian::posteriormode(void)
  {

  if (usegamma==true)
    {
    dgaussian->compute_respminuslinpred(response,0);
    unsigned i;
    double * workresponse = response.getV();
    for (i=0;i<nrobs;i++,workresponse++)
      *workresponse = *workresponse * *workresponse;

    return DISTRIBUTION_gamma::posteriormode();
    }
  else
    return true;

  }


void DISTRIBUTION_vargaussian::update(void)
  {


  dgaussian->compute_respminuslinpred(response,0);
  unsigned i;
  double * workresponse = response.getV();
  for (i=0;i<nrobs;i++,workresponse++)
  *workresponse = *workresponse * *workresponse;

  if (usegamma==true)
    {
    DISTRIBUTION_gamma::update();
    }
  else
    DISTRIBUTION::update();

  }


void DISTRIBUTION_vargaussian::update_variance(datamatrix & we)
  {

  unsigned i;
  double h = dgaussian->get_trmult(0);
  addinterceptsample = 2*log(h);

  double * workwe =  we.getV();

  double * work_lin_gamma = (*linpred_current).getV();

  for(i=0;i<nrobs;i++,work_lin_gamma++,workwe++)
    {
    *workwe = 1.0/exp(*work_lin_gamma);
    }

  }




// CONSTRUCTOR 0

DISTRIBUTION_vargaussian::DISTRIBUTION_vargaussian(const double & scale_initial,
                                                   MCMCoptions * o,
                                                   const datamatrix & r,
                                                   const ST::string & p,
                                                   const ST::string & ps,
                                                   const datamatrix & w)
  : DISTRIBUTION_gamma(scale_initial,o,r,p,ps,w)
  {

  family = "VarianceGaussian";
  usegamma=false;

  if (usegamma==false)
    scaleexisting=false;
  }



// CONSTRUCTOR 1

DISTRIBUTION_vargaussian::DISTRIBUTION_vargaussian(
                                       const double & a, const double & b,
                                       const unsigned & cit, MCMCoptions * o,
                                       const datamatrix & r,
                                       const ST::string & p,
                                       const ST::string & ps,
                                       const datamatrix & w)
  : DISTRIBUTION_gamma(a,b,cit,o,r,p,ps,w)
  {

  family = "VarianceGaussian";
  usegamma=false;

  if (usegamma==false)
    scaleexisting=false;
  }



// CONSTRUCTOR 2

DISTRIBUTION_vargaussian::DISTRIBUTION_vargaussian(
                                       const double & a, const double & b,
                                       const double & var, const unsigned & cit,
                                       MCMCoptions * o, const datamatrix & r,
                                       const ST::string & p,
                                       const ST::string & ps,
                                       const datamatrix & w)
  : DISTRIBUTION_gamma(a,b,var,cit,o,r,p,ps,w)
  {

  family = "VarianceGaussian";
  usegamma=false;

  if (usegamma==false)
    scaleexisting=false;
  }


DISTRIBUTION_vargaussian::DISTRIBUTION_vargaussian(
                          const DISTRIBUTION_vargaussian & ga)
   : DISTRIBUTION_gamma(DISTRIBUTION_gamma(ga))
  {
  dgaussian = ga.dgaussian;
  usegamma=ga.usegamma;
  }


const DISTRIBUTION_vargaussian & DISTRIBUTION_vargaussian::operator=(
                                      const DISTRIBUTION_vargaussian & ga)
  {

  if (this==&ga)
    return *this;
  DISTRIBUTION_gamma::operator=(DISTRIBUTION_gamma(ga));
  dgaussian = ga.dgaussian;
  usegamma=ga.usegamma;
  return *this;
  }


double DISTRIBUTION_vargaussian::loglikelihood(double * res, double * lin,
                                         double * w,const int & i) const
  {

  double m = exp(*lin);
  return  -  (*res)/(2*m) - 0.5* (*lin) ;

  }


void DISTRIBUTION_vargaussian::set_gaussian(DISTRIBUTION * g)
  {
  dgaussian = g;
  }


double DISTRIBUTION_vargaussian::compute_IWLS(double * response,
                                              double * linpred, double * weight,
                                              const int & i,double * weightiwls,
                                              double * tildey,bool weightyes,
                                              const unsigned & col)
  {

  double m = exp(*linpred);

  if (weightyes)
    {
    if (usegamma==true)
      *weightiwls = 1/scale(0,0) * *weight;
    else
      *weightiwls = 0.5;
    }

  *tildey = (*response - m)/m;

  if (*weight != 0)
    return  -  (*response)/(2*m) - 0.5* (*linpred);
  else
    return 0;

  }


double DISTRIBUTION_vargaussian::compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col) const
  {
  if (usegamma==true)
    return DISTRIBUTION_gamma::compute_weight(linpred,weight,i,col);
  else
    {
    return 0.5;
    }

  }


void DISTRIBUTION_vargaussian::compute_IWLS_weight_tildey(double * response,
                              double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col)
  {
  if (usegamma==true)
    {
    DISTRIBUTION_gamma::compute_IWLS_weight_tildey(
    response,linpred,weight,i,weightiwls,tildey,col);
    }
  else
    {
    *weightiwls = 0.5;

    double m = exp(*linpred);

    *tildey = (*response - m)/m;

    }

  }


void DISTRIBUTION_vargaussian::compute_iwls(void)
  {

  DISTRIBUTION_gamma::compute_iwls();
  }


double DISTRIBUTION_vargaussian::compute_gmu(double * linpred,
const unsigned & col) const
  {
  return DISTRIBUTION_gamma::compute_gmu(linpred,col);
  }



void DISTRIBUTION_vargaussian::compute_deviance(const double * response,
                                          const double * weight,
                                          const double * mu,double * deviance,
                                          double * deviancesat,
                                          const datamatrix & scale,
                                          const int & i) const
  {
  *deviance=0;
  *deviancesat=0;
  }


void DISTRIBUTION_vargaussian::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("  Response function: exp\n");

  if (usegamma==true)
    {
    if (mh)
      {
      optionsp->out("  Update of scale parameter by MH-algorithm \n");
      optionsp->out("  Fixed variance: " + ST::doubletostring(var_nu,6) + "\n");
      optionsp->out("  Hyperparameter a: " + ST::doubletostring(a_gamma,6) + "\n");
      optionsp->out("  Hyperparameter b: " + ST::doubletostring(b_gamma,6) + "\n");
      }

    else if (!mh && !scalefixed)
      {
      optionsp->out("  Update of scale parameter by consistent estimation \n");
      }

    else
      {
      optionsp->out("  Fixed scale parameter: " + ST::doubletostring(scale(0,0),6) + "\n");
      }

    optionsp->out("\n");
    optionsp->out("\n");
    }

  }


void DISTRIBUTION_vargaussian::outresults(void)
    {
    if (usegamma==true)
      {
      if (mh)
        {
        acceptancescale = double(acceptance)/
                          double(nriterations-optionsp->get_burnin())*100;

        }

      else
        {
        acceptancescale=100;
        }
      }

    DISTRIBUTION::outresults();

    }



void DISTRIBUTION_vargaussian::update_predict(void)
  {
    DISTRIBUTION::update_predict();
  }


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_gaussian --------------------------
//------------------------------------------------------------------------------

void DISTRIBUTION_gaussian::set_constscale(double s)
  {
  scale(0,0) = s/(trmult(0,0)*trmult(0,0));
  constscale=true;
  }

void DISTRIBUTION_gaussian::undo_constscale(void)
  {
  constscale = false;
  }

void DISTRIBUTION_gaussian::set_uniformprior(void)
  {
  uniformprior=true;
  }

void DISTRIBUTION_gaussian::set_variance(DISTRIBUTION_vargaussian * dg)
  {
  dgamma = dg;
  varianceest = true;
  scale(0,0)=1;
  changingweight=true;
  scaleexisting=false;
  }

void DISTRIBUTION_gaussian::standardise(void)
  {

  double s = sqrt(response.var(0,weight));

  trmult = datamatrix(1,1,s);

  unsigned i;
  double * workresp = response.getV();
  double * worklin = (*linpred_current).getV();
  for (i=0;i<nrobs;i++,workresp++,worklin++)
   {
   *workresp = *workresp/trmult(0,0);
   *worklin = *worklin/trmult(0,0);
   }


  datamatrix tr(1,1,trmult(0,0)*trmult(0,0));
  Scalesave.set_transformmult(tr);

  }

void DISTRIBUTION_gaussian::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("  Response function: identity\n");
  if(uniformprior)
    {
    optionsp->out("  Uniform prior on sigma\n");
    }
  else
    {
    optionsp->out("  Hyperparameter a: " + ST::doubletostring(a_invgamma,6) + "\n");
    optionsp->out("  Hyperparameter b: " + ST::doubletostring(b_invgamma,6) + "\n");
    }
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTRIBUTION_gaussian::update(void)
  {

  if (varianceest==true)
    {
    scale(0,0) = 1;
    dgamma->update_variance(weight);
    }


  register unsigned i;

  double help;

  double * worklin;
  double * workresp;
  double * workweight;


  // scaleparameter

  if ( (varianceest==true) || (constscale==true) )
    {
/* ------------------------- für sigma2_i bei fMRI !!! ------------------------
    unsigned j;
    double sum;

    worklin = (*linpred_current).getV();
    workresp = response.getV();
    workweight = weight.getV();

    i=0;
    while(i<nrobs)
      {
      j=0;
      sum=0.0;
      while(j<70)
        {
        help = *workresp - *worklin;
        sum += help*help;
        worklin++;workresp++;
        j++;
        i++;
        }

      help = rand_invgamma(a_invgamma+0.5*70,b_invgamma+0.5*sum);

      j=0;
      while(j<70)
        {
        *workweight = 1.0/help;
        workweight++;
        j++;
        }

      }
//----------------------------------------------------------------------------*/
    }

  else
    {
    double sum = 0;

    worklin = (*linpred_current).getV();
    workresp = response.getV();
    workweight = weight.getV();

    for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
      {

      help = *workresp - *worklin;
      sum += *workweight*help*help;
      }

    if(uniformprior==true)
      {
      double help = 1000000;
      while (help > 200000)
        help = rand_invgamma(-0.5+0.5*nrobsmweightzero,0.5*sum);
      scale(0,0) = help;
      }
    else
      {
      if(shrinkage)
        {
        scale(0,0) = rand_invgamma(a_invgamma+0.5*nrobsmweightzero + 0.5*nrridge + 0.5*nrlasso + 0.5*nrnigmix,
                     b_invgamma+0.5*sum + 0.5*ridgesum + 0.5*lassosum + 0.5*nigmixsum);
//TEMP:BEGIN--------------------------------------------------------------------
// Pruefen der uebergebenen Optionen
/*
int iteration = optionsp->get_nriter();
ofstream output("c:/bayesx/test/scale.txt", ios::out|ios::app);
if (iteration ==1)
{output <<"iter " << "a_ig " << "b_ig " << "obswzero " << "sum " << "scale "<< "nrridge "<< "nrlasso " << "ridgesum " << "lassosum " << "\n" ;
}
output << iteration << " " << a_invgamma << " " << b_invgamma << " " << nrobsmweightzero << " " << sum << " " << scale(0,0) << " " << nrridge << " " << nrlasso << " " << ridgesum << " " << lassosum << "\n" ;
*/
//TEMP:END----------------------------------------------------------------------
        }
      else
        {
        scale(0,0) = rand_invgamma(a_invgamma+0.5*nrobsmweightzero,
                     b_invgamma+0.5*sum);
        }
      }
    }


  // Prediction

  if (predictresponse==true)
    {
    worklin = (*linpred_current).getV();
    workresp = response.getV();
    workweight = weight.getV();
    double * workpredictind = predictindicator.getV();
    double sscale = sqrt(scale(0,0));

    for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++,workpredictind++)
      {
      if (*workpredictind == 0)
        {
        *workresp = *worklin + (sscale/(*workweight))*rand_normal();
        }

      }
    }


  DISTRIBUTION::update();


  if (fcmissing.size() > 0)
    update_missings();

  }


void DISTRIBUTION_gaussian::update_missings(void)
  {
  unsigned i;
  for (i=0;i<fcmissing.size();i++)
    {
    fcmissing[i]->update_missings(response,
    (*linpred_current),missingpos,responsename,scale(0,0));
    }

  unsigned * workmissingpos = missingpos.getV();

  double * workresponse = response.getV()+(*workmissingpos);

  double * betap = MissingSave.getbetapointer();

  for(i=0;i<missingpos.rows();i++,betap++,workmissingpos++,
      workresponse+=(*workmissingpos))
    {
    *betap = *workresponse;
    }

  MissingSave.update();

  }


double DISTRIBUTION_gaussian::loglikelihood(double * res,
                       double * lin,
                       double * w,
                       const int & i) const
  {
  double help = *res-*lin;
  return  - *w * ( help * help )/(2* scale(0,0));
  }


// ------------------- For Stepwise --------------------------------------------

double DISTRIBUTION_gaussian::compute_rss(void)
    {
  unsigned i;

  double * worklin = (*linpred_current).getV();
  double * workresp = response.getV();
  double * workweight = weight.getV();

  double sum = 0;
  double help;

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    if (*workweight !=0)
      {
      help = *workresp - *worklin;
      sum += *workweight*help*help;
      }
    }
  sum = trmult(0,0)*trmult(0,0)*sum;

  return  sum;
  }



//------------------------------------------------------------------------------
//------------------------- CLASS: DISTRIBUTION_AFT ----------------------------
//------------------------------------------------------------------------------


DISTRIBUTION_AFT::DISTRIBUTION_AFT(void) : DISTRIBUTION_gaussian()
  {
  }



DISTRIBUTION_AFT::DISTRIBUTION_AFT(const double & a,
                                             const double & b,
                                             MCMCoptions * o,
                                             const datamatrix & r,
                                             const datamatrix & cens,
                                             const ST::string & p,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTRIBUTION_gaussian(a,b,o,r,p,ps,w)

  {
  censoring = cens;
  responseorig = response;
  }


// constructor with offset
DISTRIBUTION_AFT::DISTRIBUTION_AFT(const datamatrix & offset,
                      const double & a, const double & b, MCMCoptions * o,
                      const datamatrix & r,const datamatrix & cens,
                      const ST::string & p,
                      const ST::string & ps,const datamatrix & w)
  : DISTRIBUTION_gaussian(offset,a,b,o,r,p,ps,w)

  {
  censoring = cens;
  responseorig = response;
  }


// COPY CONSTRUCTOR

DISTRIBUTION_AFT::DISTRIBUTION_AFT(const DISTRIBUTION_AFT & nd)
   : DISTRIBUTION_gaussian(DISTRIBUTION_gaussian(nd))
  {
  censoring = nd.censoring;
  responseorig = nd.responseorig;
  }


const DISTRIBUTION_AFT & DISTRIBUTION_AFT::operator=(
                                      const DISTRIBUTION_AFT & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION_gaussian::operator=(DISTRIBUTION_gaussian(nd));
  censoring = nd.censoring;
  responseorig = nd.responseorig;
  return *this;
  }


void DISTRIBUTION_AFT::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const
    {

    }


void DISTRIBUTION_AFT::outoptions(void)
  {
  DISTRIBUTION_gaussian::outoptions();
  }


void DISTRIBUTION_AFT::update(void)
  {
  unsigned i;

//  double test = trmult(0,0);

  double * rp = response.getV();
  double * rpo = responseorig.getV();
  double * cp = censoring.getV();
  double * worklin = (*linpred_current).getV();
  double sqrtscale = sqrt(scale(0,0));
  for(i=0; i< responseorig.rows(); i++, rp++, rpo++, cp++,worklin++)
    {
    if(*cp==0)
      *rp = trunc_normal4(*rpo, *worklin, sqrtscale);
    }
  DISTRIBUTION_gaussian::update();
  }


bool DISTRIBUTION_AFT::posteriormode(void)
  {

  return true;
  }

bool DISTRIBUTION_AFT::posteriormode_converged_fc(const datamatrix & beta,
                                  const datamatrix & beta_mode,
                                  const unsigned & itnr)
  {
  return true;
  }


void DISTRIBUTION_AFT::outresults(void)
  {
  DISTRIBUTION_gaussian::outresults();
  }


  /*
  void set_constscale(double s);

  void undo_constscale(void);

  void set_uniformprior(void);

  void update_missings(void);
  */


//------------------------------------------------------------------------------
//---------------------- CLASS: DISTRIBUTION_QUANTREG --------------------------
//------------------------------------------------------------------------------
#if !defined (__BUILDING_THE_DLL)

DISTRIBUTION_QUANTREG::DISTRIBUTION_QUANTREG(void) : DISTRIBUTION_gaussian()
  {
  }



DISTRIBUTION_QUANTREG::DISTRIBUTION_QUANTREG(const double & a,
                                             const double & b,
                                             MCMCoptions * o,
                                             const datamatrix & r,
                                             const ST::string & p,
                                             const ST::string & ps,
                                             const double & quant,
                                             const datamatrix & w)
  : DISTRIBUTION_gaussian(a,b,o,r,p,ps,w)

  {
  responseorig = response;
  weightorig = w;
  changingweight = true;
  quantile = quant;
  xi = (1-2*quant)/(quant * (1-quant));
  sigma02 = 2 / (quant*(1-quant));
  }


// constructor with offset
DISTRIBUTION_QUANTREG::DISTRIBUTION_QUANTREG(const datamatrix & offset,
                      const double & a, const double & b, MCMCoptions * o,
                      const datamatrix & r,
                      const ST::string & p,
                      const ST::string & ps,
                      const double & quant,
                      const datamatrix & w)
  : DISTRIBUTION_gaussian(offset,a,b,o,r,p,ps,w)

  {
  responseorig = response;
  weightorig = w;
  changingweight = true;
  quantile = quant;
  xi = (1-2*quant)/(quant * (1-quant));
  sigma02 = 2 / (quant*(1-quant));
  }


// COPY CONSTRUCTOR

DISTRIBUTION_QUANTREG::DISTRIBUTION_QUANTREG(const DISTRIBUTION_QUANTREG & nd)
   : DISTRIBUTION_gaussian(DISTRIBUTION_gaussian(nd))
  {
  responseorig = nd.responseorig;
  weightorig = nd.weightorig;
  quantile = nd.quantile;
  xi = nd.xi;
  sigma02 = nd.sigma02;
  }


const DISTRIBUTION_QUANTREG & DISTRIBUTION_QUANTREG::operator=(
                                      const DISTRIBUTION_QUANTREG & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION_gaussian::operator=(DISTRIBUTION_gaussian(nd));
  weightorig = nd.weightorig;
  responseorig = nd.responseorig;
  quantile = nd.quantile;
  xi = nd.xi;
  sigma02 = nd.sigma02;
  return *this;
  }


void DISTRIBUTION_QUANTREG::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const
    {

    }


void DISTRIBUTION_QUANTREG::outoptions(void)
  {
  DISTRIBUTION_gaussian::outoptions();
  }


void DISTRIBUTION_QUANTREG::update(void)
  {
  unsigned i;

  double * rp = response.getV();
  double * rpo = responseorig.getV();
  double * w = weight.getV();
  double * wo = weightorig.getV();
  double * worklin = (*linpred_current).getV();

  double xi2 = xi*xi;
  double lambda = (xi2 + 2*sigma02)/ (scale(0,0));
  double num = sqrt(xi2 + 2*sigma02);
  double mu;

  double sumres=0;
  double sumw=0;

  for(i=0; i<weight.rows(); i++, w++, rp++, worklin++, rpo++, wo++)
    {
    if(*wo != 0)
      {
      mu = num/(fabs(*rpo-*worklin));
      *w = rand_inv_gaussian(mu, lambda);
      *rp = *rpo - xi / *w;

      sumw += 1 / *w;
      sumres += *w * (*rp-*worklin)*(*rp-*worklin);
      }
    }

  sumres /= sigma02;

  if(shrinkage)
    {
    scale(0,0) = rand_invgamma(a_invgamma + 1.5*nrobsmweightzero + 0.5*nrridge + 0.5*nrlasso + 0.5*nrnigmix,
                 b_invgamma + 0.5*sumres + sumw + 0.5*ridgesum + 0.5*lassosum + 0.5*nigmixsum);
    }
  else
    {
    scale(0,0) = rand_invgamma(a_invgamma + 1.5*nrobsmweightzero,
                    b_invgamma + 0.5*sumres + sumw);
    }

  DISTRIBUTION_gaussian::update();

  scale(0,0) *= sigma02;

  }


bool DISTRIBUTION_QUANTREG::posteriormode(void)
  {

  return true;
  }

bool DISTRIBUTION_QUANTREG::posteriormode_converged_fc(const datamatrix & beta,
                                  const datamatrix & beta_mode,
                                  const unsigned & itnr)
  {
  return true;
  }


void DISTRIBUTION_QUANTREG::outresults(void)
  {
  DISTRIBUTION_gaussian::outresults();
  }
#endif

//------------------------------------------------------------------------------
//----------------------- END: DISTRIBUTION_QUANTREG ---------------------------
//------------------------------------------------------------------------------



/*double DISTRIBUTION_gaussian::compute_logmsep(void)
  {
  unsigned i;

  double * worklin = (*linpred_current).getV();
  double * workresp = response.getV();
  double * workweight = weight.getV();
  double * work2 = weight2.getV();

  double sum = 0;
  double help;
  double nr = 0;
  double sigma2_halbe = trmult(0,0) * compute_rss()/(2 * get_nrobs_wpw());

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++,work2++)
    {
    if (*workweight == 0)
      {
      help = exp(trmult(0,0) * *workresp) - exp(trmult(0,0) * *worklin + sigma2_halbe);
      //sum += *workweight*help*help;     // u.U. ändern, so dass zwei Variablen, eine mit
                                          // Gewichten != 0, die andere 0/1 Variable
      sum += help*help * *work2;
      nr += 1;
      }
    }

  return  sum;     // /nr;
  } */


double DISTRIBUTION_gaussian::compute_msep(void)
  {
  unsigned i;

  double * worklin = (*linpred_current).getV();
  double * workresp = response.getV();
  double * workweight = weight.getV();
  double * work2 = weight2.getV();

  double sum = 0;
  double help;
  double nr = 0;

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++,work2++)
    {
    if (*workweight == 0)
      {
      help = *workresp - *worklin;
      //sum += *workweight*help*help;     // u.U. ändern, so dass zwei Variablen, eine mit
                                          // Gewichten != 0, die andere 0/1 Variable
      sum += help*help * *work2;
      nr += 1;
      }
    }
  sum = trmult(0,0)*trmult(0,0)*sum;

  return  sum;     // /nr;
  }


double DISTRIBUTION_gaussian::compute_gcv(const double & df)
  {
  double help2 = 1 - gcvfactor * df/get_nrobs_wpw();
  double gcv = compute_rss() / (get_nrobs_wpw()*help2*help2);
  return gcv;
  }


double DISTRIBUTION_gaussian::compute_aic(const double & df)
  {
  double aic = log(compute_rss() / get_nrobs_wpw()) * get_nrobs_wpw() + 2*df;
  return aic;
  }

double DISTRIBUTION_gaussian::compute_improvedaic(const double & df)
  {
  double impaic = log(compute_rss() / get_nrobs_wpw()) * get_nrobs_wpw()
                  + 2*df + 2*df*(df+1)/(get_nrobs_wpw()-df-1);
  return impaic;
  }

double DISTRIBUTION_gaussian::compute_bic(const double & df)
  {
  double bic = log(compute_rss() / get_nrobs_wpw()) * get_nrobs_wpw()
               + log(static_cast<double>(get_nrobs_wpw()))*df;
  return bic;
  }


void DISTRIBUTION_gaussian::compute_bootstrap_data(const double * linpred,const double * weight,double * wresp)
  {
  double pres = 0;
  if(*weight > 0)
    pres = *linpred + sqrt(scale(0,0) / *weight) * rand_normal();
  *wresp = pres;
  }

// ------------- END: stepwise -------------------------------------------------


bool DISTRIBUTION_gaussian::posteriormode(void)
  {

  unsigned i;

  double * worklin = (*linpred_current).getV();
  double * workresp = response.getV();
  double * workweight = weight.getV();

  double sum = 0;
  double help;
  double sumweight=0;

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - *worklin;
    sum += *workweight*help*help;
    sumweight += *workweight;
    }


  if (constscale==false)
    scale(0,0) = (1.0/sumweight)*sum;


  return true;


  }


// constructor without offset
DISTRIBUTION_gaussian::DISTRIBUTION_gaussian(const double & a,
                                             const double & b,
                                             MCMCoptions * o,
                                             const datamatrix & r,
                                             const ST::string & p,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTRIBUTION(o,r,w,p,ps)

  {
// auskommentiert um nichtinformative Prioris zu erlauben
//  assert (a > 0);
//  assert (b > 0);

  constscale=false;
  uniformprior=false;

/* ---------------------- für sigma2_i bei fMRI !!! ---------------------------
  scale(0,0)=1;
  changingweight=true;
  constscale=true;
  scaleexisting=false;
//----------------------------------------------------------------------------*/

  a_invgamma = a;
  b_invgamma = b;
  family = "Gaussian";

  standardise();

  acceptancescale=100;

  varianceest=false;

  constant_iwlsweights=true;
  iwlsweights_notchanged_df = true;

  }


// constructor with offset
DISTRIBUTION_gaussian::DISTRIBUTION_gaussian(const datamatrix & offset,
                      const double & a, const double & b, MCMCoptions * o,
                      const datamatrix & r,const ST::string & p,
                      const ST::string & ps,const datamatrix & w)
  : DISTRIBUTION(offset,o,r,w,p,ps)

  {

// auskommentiert um nichtinformative Prioris zu erlauben
//  assert (a > 0);
//  assert (b > 0);

  constscale=false;
  uniformprior=false;

  a_invgamma = a;
  b_invgamma = b;
  family = "Gaussian";

  standardise();

  acceptancescale=100;

  varianceest=false;

  constant_iwlsweights=true;
  iwlsweights_notchanged_df = true;

  }


const DISTRIBUTION_gaussian & DISTRIBUTION_gaussian::operator=(
                                      const DISTRIBUTION_gaussian & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));
  a_invgamma = nd.a_invgamma;
  b_invgamma = nd.b_invgamma;
  varianceest = nd.varianceest;
  constscale = nd.constscale;
  uniformprior = nd.uniformprior;
  return *this;
  }

// BEGIN: DSB //
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
double
    DISTRIBUTION_gaussian::sample_from_likelihood(const double weight,
                                                  const double mu) const
    {
        // get the current "scale" parameter of the normal likelihood
        double s =  scale(0,0) * pow(trmult(0,0), 2);
        double sigma2 = (weight != 0.0) ? s / weight : s;

        // and then use the sample function from Random.cpp
        return mu + sqrt(sigma2) * randnumbers::rand_normal();
    }
#endif
// END: DSB //


void DISTRIBUTION_gaussian::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const
  {
  if ((*weight != 0) && (varianceest == false))
    {
    double s = scale(0,0)*pow(trmult(0,0),2);
    double r = *response*trmult(0,0)-*mu;
    *deviance =  (*weight/s)*r*r+log(2*M_PI*s/(*weight));
    *deviancesat = (*weight/s)*r*r;
    }
  else
    {
    *deviance = 0;
    *deviancesat = 0;
    }
  }


double DISTRIBUTION_gaussian::compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col) const
  {
  return *weight;
  }


void DISTRIBUTION_gaussian::compute_mu(const double * linpred,double * mu) const
  {
  *mu = trmult(0,0)* *linpred;
  }


void DISTRIBUTION_gaussian::compute_mu_notransform(
const double * linpred,double * mu) const
  {
  *mu = *linpred;
  }


void  DISTRIBUTION_gaussian::tr_nonlinear(vector<double *> b,
                                          vector<double *> br,
                                          vector<FULLCOND*> & fcp,
                                          unsigned & nr,
                                          unsigned & it,ST::string & trtype)
  {
  if (trtype == "exp")
    DISTRIBUTION::tr_nonlinear(b,br,fcp,nr,it,trtype);
  else if(trtype == "lognormal")
    {
    datamatrix help(1,1);
    Scalesave.readsample2(help,it);
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0)+*b[i]+help(0,0)/2.0);
      }
    }
  else if (trtype == "elasticity")
    {
    if (b.size() == 2)
      {
      *br[1] = *b[1] * fcp[0]->get_data(nr,0) /(interceptsample(it,0)+ *b[0]);
      }
    }
  else if (trtype=="marginal")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = interceptsample(it,0)+ *b[i];
      }
    }
  else if (trtype=="marginalintercept")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = interceptsample(it,0);
      }
    }
  else if (trtype == "lognormalintercept")
    {
    datamatrix help(1,1);
    Scalesave.readsample2(help,it);
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0)+help(0,0)/2.0);
      }
    }

  }


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_gaussian_re -----------------------
//------------------------------------------------------------------------------

void DISTRIBUTION_gaussian_re::set_constscale(double s)
  {
  scale(0,0) = s/(trmult(0,0)*trmult(0,0));
  constscale=true;
  }

void DISTRIBUTION_gaussian_re::undo_constscale(void)
  {
  constscale = false;
  }

void DISTRIBUTION_gaussian_re::set_uniformprior(void)
  {
  uniformprior=true;
  }



void DISTRIBUTION_gaussian_re::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("  Response function: identity\n");
  if(uniformprior)
    {
    optionsp->out("  Uniform prior on sigma\n");
    }
  else
    {
    optionsp->out("  Hyperparameter a: " + ST::doubletostring(a_invgamma,6) + "\n");
    optionsp->out("  Hyperparameter b: " + ST::doubletostring(b_invgamma,6) + "\n");
    }
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTRIBUTION_gaussian_re::update(void)
  {

  register unsigned i;

  double help;

  double * worklin;
  double * workresp;
  double * workweight;


  // scaleparameter

  double sum = 0;

  worklin = (*linpred_current).getV();
  workresp = response.getV();
  workweight = weight.getV();

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {

    help = *workresp - *worklin;
    sum += *workweight*help*help;
    }

  if(uniformprior==true)
    {
    double help = 1000000;
    while (help > 200000)
      help = rand_invgamma(-0.5+0.5*nrobsmweightzero,0.5*sum);
    scale(0,0) = help;
    }
  else
    {
    if(shrinkage)
      {
      scale(0,0) = rand_invgamma(a_invgamma+0.5*nrobsmweightzero + 0.5*nrridge + 0.5*nrlasso,
                     b_invgamma+0.5*sum + 0.5*ridgesum + 0.5*lassosum);
      }
    else
      {
      scale(0,0) = rand_invgamma(a_invgamma+0.5*nrobsmweightzero,
                     b_invgamma+0.5*sum);
      }
    }

  // Prediction

  if (predictresponse==true)
    {
    worklin = (*linpred_current).getV();
    workresp = response.getV();
    workweight = weight.getV();
    double * workpredictind = predictindicator.getV();
    double sscale = sqrt(scale(0,0));

    for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++,workpredictind++)
      {
      if (*workpredictind == 0)
        {
        *workresp = *worklin + (sscale/(*workweight))*rand_normal();
        }

      }
    }


  DISTRIBUTION::update();

  }


double DISTRIBUTION_gaussian_re::loglikelihood(double * res,
                       double * lin,
                       double * w,
                       const int & i) const
  {
  double help = *res-*lin;
  return  - *w * ( help * help )/(2* scale(0,0));
  }


bool DISTRIBUTION_gaussian_re::posteriormode(void)
  {

  unsigned i;

  double * worklin = (*linpred_current).getV();
  double * workresp = response.getV();
  double * workweight = weight.getV();

  double sum = 0;
  double help;
  double sumweight=0;

  for (i=0;i<nrobs;i++,worklin++,workresp++,workweight++)
    {
    help = *workresp - *worklin;
    sum += *workweight*help*help;
    sumweight += *workweight;
    }


  if (constscale==false)
    scale(0,0) = (1.0/sumweight)*sum;


  return true;

  }


// constructor without offset
DISTRIBUTION_gaussian_re::DISTRIBUTION_gaussian_re(const double & a,
                                             const double & b,
                                             MCMCoptions * o,
                                             const datamatrix & r,
                                             const ST::string & p,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTRIBUTION(o,r,w,p,ps)

  {
// auskommentiert um nichtinformative Prioris zu erlauben
//  assert (a > 0);
//  assert (b > 0);

  constscale=false;
  uniformprior=false;

  a_invgamma = a;
  b_invgamma = b;
  family = "Gaussian_re";

  acceptancescale=100;

  constant_iwlsweights=true;
  iwlsweights_notchanged_df = true;

  }


// constructor with offset
DISTRIBUTION_gaussian_re::DISTRIBUTION_gaussian_re(const datamatrix & offset,
                      const double & a, const double & b, MCMCoptions * o,
                      const datamatrix & r,const ST::string & p,
                      const ST::string & ps,const datamatrix & w)
  : DISTRIBUTION(offset,o,r,w,p,ps)

  {

// auskommentiert um nichtinformative Prioris zu erlauben
//  assert (a > 0);
//  assert (b > 0);

  constscale=false;
  uniformprior=false;

  a_invgamma = a;
  b_invgamma = b;
  family = "Gaussian_re";

  acceptancescale=100;

  constant_iwlsweights=true;
  iwlsweights_notchanged_df = true;

  }


const DISTRIBUTION_gaussian_re & DISTRIBUTION_gaussian_re::operator=(
                                      const DISTRIBUTION_gaussian_re & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));
  a_invgamma = nd.a_invgamma;
  b_invgamma = nd.b_invgamma;
  constscale = nd.constscale;
  uniformprior = nd.uniformprior;
  return *this;
  }


// BEGIN: DSB //
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
// this is the same as DISTRIBUTION_gaussian::sample_from_likelihood
// todo: perhaps the other function is not needed because there cannot be any random effects??
double
    DISTRIBUTION_gaussian_re::sample_from_likelihood(const double weight,
                                                     const double mu) const
    {
        // get the current "scale" parameter of the normal likelihood
        double s =  scale(0,0) * pow(trmult(0,0), 2);
        double sigma2 = (weight != 0.0) ? s / weight : s;

        // and then use the sample function from Random.cpp
        return mu + sqrt(sigma2) * randnumbers::rand_normal();
    }
#endif
// END: DSB //

void DISTRIBUTION_gaussian_re::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const
  {
  if ((*weight != 0))
    {
    double s = scale(0,0)*pow(trmult(0,0),2);
    double r = *response*trmult(0,0)-*mu;
    *deviance =  (*weight/s)*r*r+log(2*M_PI*s/(*weight));
    *deviancesat = (*weight/s)*r*r;
    }
  else
    {
    *deviance = 0;
    *deviancesat = 0;
    }
  }


double DISTRIBUTION_gaussian_re::compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col) const
  {
  return *weight;
  }


void DISTRIBUTION_gaussian_re::compute_mu(const double * linpred,double * mu) const
  {
  *mu = trmult(0,0)* *linpred;
  }


void DISTRIBUTION_gaussian_re::compute_mu_notransform(
const double * linpred,double * mu) const
  {
  *mu = *linpred;
  }


void  DISTRIBUTION_gaussian_re::tr_nonlinear(vector<double *> b,
                                          vector<double *> br,
                                          vector<FULLCOND*> & fcp,
                                          unsigned & nr,
                                          unsigned & it,ST::string & trtype)
  {
  if (trtype == "exp")
    DISTRIBUTION::tr_nonlinear(b,br,fcp,nr,it,trtype);
  else if(trtype == "lognormal")
    {
    datamatrix help(1,1);
    Scalesave.readsample2(help,it);
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0)+*b[i]+help(0,0)/2.0);
      }
    }
  else if (trtype == "elasticity")
    {
    if (b.size() == 2)
      {
      *br[1] = *b[1] * fcp[0]->get_data(nr,0) /(interceptsample(it,0)+ *b[0]);
      }
    }
  else if (trtype=="marginal")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = interceptsample(it,0)+ *b[i];
      }
    }
  else if (trtype=="marginalintercept")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = interceptsample(it,0);
      }
    }
  else if (trtype == "lognormalintercept")
    {
    datamatrix help(1,1);
    Scalesave.readsample2(help,it);
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0)+help(0,0)/2.0);
      }
    }

  }


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_lognormal -------------------------
//------------------------------------------------------------------------------


double DISTRIBUTION_lognormal::loglikelihood(double * res,
                       double * lin,
                       double * w,
                       const int & i) const
  {
  double help = *res-*lin;
  return  - *w * ( help * help )/(2* scale(0,0));
  }


// constructor without offset
DISTRIBUTION_lognormal::DISTRIBUTION_lognormal(const double & a,
                                             const double & b,
                                             MCMCoptions * o,
                                             const datamatrix & r,
                                             const ST::string & p,
                                             const ST::string & ps,
                                             const datamatrix & w)
  : DISTRIBUTION_gaussian(a,b,o,r,p,ps,w)

  {
  family = "lognormal";
  }


// constructor with offset
DISTRIBUTION_lognormal::DISTRIBUTION_lognormal(const datamatrix & offset,
                      const double & a, const double & b, MCMCoptions * o,
                      const datamatrix & r,const ST::string & p,
                      const ST::string & ps,const datamatrix & w)
  : DISTRIBUTION_gaussian(offset,a,b,o,r,p,ps,w)

  {
  family = "lognormal";
  }


const DISTRIBUTION_lognormal & DISTRIBUTION_lognormal::operator=(
                                      const DISTRIBUTION_lognormal & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION_gaussian::operator=(DISTRIBUTION_gaussian(nd));
  return *this;
  }


// BEGIN: DSB //
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
double
    DISTRIBUTION_lognormal::sample_from_likelihood(const double weight,
                                                   const double mu) const
    {
        // get the current scale parameter of the normal likelihood
        double s = scale(0,0) * pow(trmult(0,0), 2);

        // and the current mean parameter
        double m = log(mu) - 0.5 * s;
        // todo: why doesn't this depend on sigma2 but on s? It is also like that in the compute_deviance function below!

        // the variance parameter then depends again on the weight
        double sigma2 = (weight != 0.0) ? s / weight : s;

        // finally we can draw from the lognormal distribution
        return exp(m + sqrt(sigma2) * randnumbers::rand_normal());
    }
#endif
// END: DSB //


void DISTRIBUTION_lognormal::compute_deviance(const double * response,
                           const double * weight,
                           const double * mu, double * deviance,
                           double * deviancesat,
                           const datamatrix & scale,const int & i) const
  {
  if ((*weight != 0) && (varianceest == false))
    {
    double s = scale(0,0)*pow(trmult(0,0),2);
    double r = *response*trmult(0,0)-(log(*mu)-0.5*s);
    *deviance = (*weight/s)*r*r+log(2*M_PI*s/(*weight)) + 2.0**response*trmult(0,0);
    *deviancesat = (*weight/s)*r*r;
    }
  else
    {
    *deviance = 0;
    *deviancesat = 0;
    }
  }


void DISTRIBUTION_lognormal::compute_mu(const double * linpred,double * mu) const
  {
  *mu = exp(trmult(0,0)* *linpred + scale(0,0)*trmult(0,0)*trmult(0,0)/2.0);
  }


void DISTRIBUTION_lognormal::compute_mu_notransform(
const double * linpred,double * mu) const
  {
  *mu = exp(*linpred + scale(0,0)*pow(trmult(0,0),2)/2.0);
  }


void DISTRIBUTION_lognormal::tr_nonlinear(vector<double *> b,vector<double *> br,
                    vector<FULLCOND*> & fcp,unsigned & nr,
                    unsigned & it,ST::string & trtype)
  {

  if (trtype == "exp")
    DISTRIBUTION::tr_nonlinear(b,br,fcp,nr,it,trtype);
  else if (trtype == "elasticity")
    {
    if (b.size() == 2)
      {
      *br[1] = *b[1] * fcp[0]->get_data(nr,0);
      }
    }
  else if (trtype=="marginal")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0)+ *b[i]);
      }
    }
  else if (trtype=="marginalintercept")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0));
      }
    }

  }


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_binomial --------------------------
//------------------------------------------------------------------------------


void DISTRIBUTION_binomial::tr_nonlinear(vector<double *> b,vector<double *> br,
                                         vector<FULLCOND*> & fcp,unsigned & nr,
                                         unsigned & it, ST::string & trtype)
  {
  if (trtype == "exp")
    DISTRIBUTION::tr_nonlinear(b,br,fcp,nr,it,trtype);
  else if ((trtype == "logit") || (trtype=="marginal") )
    {
    double h;
    double eh;
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      h = *b[i]+interceptsample(it,0);
      eh = exp(h);
      *br[i] = eh/(1+eh);
      }
    }
  else if (trtype == "oddsratio")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(*b[i]);
      }
    }
  else if ( (trtype == "logitintercept") || (trtype=="marginalintercept")  )
    {
    double eh;
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      eh = exp(interceptsample(it,0));
      *br[i] = eh/(1+eh);
      }
    }

  }


void DISTRIBUTION_binomial::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("  Response function: logistic distribution function\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTRIBUTION_binomial::create(void)
  {

  family = "Binomial (logit link)";
  scale(0,0) = 1;
  scaleexisting = false;

  bool error=false;
  unsigned i=0;
  double * workresp = response.getV();
  double * workweight = weight.getV();
  while ( (i<nrobs) && (error==false) )
    {

    if (*workweight > 0)
      {
      if (*workresp != int(*workresp))
        {
        error=true;
        errors.push_back("ERROR: response cannot be binomial; values must be integer numbers\n");
        }

      if (*workresp < 0)
        {
        error=true;
        errors.push_back("ERROR: response cannot be binomial; some values are negative\n");
        }

      if (*workresp > *workweight)
        {
        error = true;
        errors.push_back("ERROR: response cannot be binomial;\n");
        errors.push_back("       number of successes larger than number of trials for some values\n");
        }

      *workresp = *workresp/(*workweight);
      }

    i++;
    workresp++;
    workweight++;

    }

  }


void DISTRIBUTION_binomial::compute_mu(const double * linpred,
                                       double * mu) const
  {
  double el = exp(*linpred);
  *mu = el/(1+el);
  }


void DISTRIBUTION_binomial::compute_mu_notransform(
const double * linpred,double * mu) const
  {
  double el = exp(*linpred);
  *mu = el/(1+el);
  }


DISTRIBUTION_binomial::DISTRIBUTION_binomial(MCMCoptions * o,
                         const datamatrix & r,const datamatrix & w)
  : DISTRIBUTION(o,r,w)
  {

  create();

  }


DISTRIBUTION_binomial::DISTRIBUTION_binomial(const datamatrix & offset,
                MCMCoptions * o,const datamatrix & r, const datamatrix & w)
  : DISTRIBUTION(offset,o,r,w)
  {

  create();

  }


double DISTRIBUTION_binomial::loglikelihood(double * response,
                                            double * linpred,
                                            double * weight,
                                            const int & i) const
  {

  if (*linpred >= 10)
    return *weight *(*response * *linpred - *linpred);
  else
    return *weight *(*response * *linpred - log(1+exp(*linpred)));

  }

double DISTRIBUTION_binomial::compute_weight(double * linpred, double * weight,
                                             const int & i, const unsigned & col)
                                             const
  {
  double el = exp(*linpred);
  double mu = el/(1+el);
  return mu*(1-mu)* *weight;
  }


void DISTRIBUTION_binomial::compute_IWLS_weight_tildey(double * response,
                      double * linpred, double * weight,
                      const int & i,double * weightiwls,
                      double * tildey, const unsigned & col)
  {

  double el = exp(*linpred);
  double mu = el/(1+el);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;

  double v =mu*(1-mu);

  *weightiwls = v* *weight;

  *tildey = (*response - mu)/v;

  }


double DISTRIBUTION_binomial::compute_IWLS(double * response,double * linpred,
                            double * weight,const int & i,
                            double * weightiwls, double * tildey,
                            bool weightyes,
                            const unsigned & col)
  {
  double el = exp(*linpred);
  double mu = el/(1+el);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;
  double v = mu*(1-mu);

  if (weightyes)
    *weightiwls = *weight * v;

  *tildey = (*response - mu)/v;

  if (*linpred >= 10)
    return *weight *(*response * *linpred - *linpred);
  else
    return *weight *(*response * *linpred - log(1+el));

  }


double DISTRIBUTION_binomial::compute_gmu(double * linpred,
                                          const unsigned & col) const
  {
  double el = exp(*linpred);
  double mu = el/(1+el);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;
  return 1.0/(mu*(1-mu));
  }


void DISTRIBUTION_binomial::update(void)
  {

  DISTRIBUTION::update();
  }


bool DISTRIBUTION_binomial::posteriormode(void)
  {

  return true;

  }


bool DISTRIBUTION_binomial::posteriormode_converged(const unsigned & itnr)
  {

  return true;

  }


// BEGIN: DSB //
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
double
    DISTRIBUTION_binomial::sample_from_likelihood(const double weight,
                                                  const double mu) const
    {
        // determine the size parameter
        double size = (weight != 0.0) ? weight : 1.0;

        // use the sample function from Random.cpp
        return randnumbers::rand_binom(size, mu);
    }
#endif
// END: DSB //


void DISTRIBUTION_binomial::compute_deviance(const double * response,
               const double * weight,const double * mu,double * deviance,
               double * deviancesat,
               const datamatrix & scale,const int & i) const
  {

  if (*weight > 0)
    {
    if (*response==0)
      {
      *deviance = -2* *weight * log(1-*mu);
      *deviancesat = *deviance;
      }
    else if (*response == 1)
      {
      *deviance = -2* *weight*log(*mu);
      *deviancesat = *deviance;
      }
    else
      {
      *deviance = -2* *weight*( *response*log(*mu)+(1-*response)*log(1-*mu) );
      *deviancesat = *deviance +
      2* *weight*( *response*log(*response)+(1-*response)*log(1-*response) );
      }
    }
  else
    {
    *deviance = 0;
    *deviancesat = 0;
    }

  }


double DISTRIBUTION_binomial::compute_auc(void)
  {
  datamatrix linpred_null_hilf = datamatrix(nrobs,1,0);       // Erzeugen einer Matrix die Einträge für Beob.
                                                              // mit weight=0 enthält und darunter Nuller
  datamatrix response_null_hilf = datamatrix(nrobs,1,0);      // Erzeugen einer Matrix die Einträge für Beob.
                                                              // mit weight=0 enthält und darunter Nuller
  double * worklin = linearpred.getV();
  double * lin_null = linpred_null_hilf.getV();
  double * workweight = weight.getV();
  double * workresp = response.getV();
  double * resp_null = response_null_hilf.getV();
  unsigned nrnull = 0;
  unsigned nrdef = 0;
  unsigned i;
  for(i=0;i<response.rows();i++,workresp++,worklin++,workweight++)
    {
    if(*workweight == 0)
      {
      *lin_null = *worklin;
      lin_null++;
      *resp_null = *workresp;
      resp_null++;
      nrnull += 1;
      if(*workresp==1)
         nrdef += 1;
      }
    }

  datamatrix linpred_null = datamatrix(nrnull,1,0);  // Kürzen des "Null"-Prädiktors, d.h. Weglassen der Nullen
  datamatrix response_null = datamatrix(nrnull,1,0); // Kürzen des "Null"-Response, d.h. Weglassen der Nullen
  double * lin_hilf = linpred_null_hilf.getV();
  lin_null = linpred_null.getV();
  double * resp_hilf = response_null_hilf.getV();
  resp_null = response_null.getV();
  for(i=0;i<nrnull;i++,lin_null++,lin_hilf++,resp_null++,resp_hilf++)
     {
     *lin_null = *lin_hilf;
     *resp_null = *resp_hilf;
     }

  statmatrix<int> index_gesamt(linpred_null.rows(),1);      // Sortieren des linearen Prädiktors
  index_gesamt.indexinit();
  linpred_null.indexsort(index_gesamt,0,linpred_null.rows()-1,0,0);
  statmatrix<double> rang_gepoolt(linpred_null.rows(),1);
  linpred_null.rank(rang_gepoolt,index_gesamt,0,linpred_null.rows()-1,0);  // Ränge berechnen

  int* gesamt = index_gesamt.getV();
  double* rang_gep = rang_gepoolt.getV();
  unsigned rang_default = 1;
  double auc = 0;

  for(i=0;i<nrnull;i++,gesamt++,rang_gep++)
    {
    if(response_null.get(*gesamt,0) == 1)
      {
      //auc += i+1 - rang_default;
      auc += *rang_gep - rang_default;
      rang_default += 1;
      }
    }

  return auc/(nrdef*(nrnull-nrdef));
  }


void DISTRIBUTION_binomial::compute_bootstrap_data(const double * linpred,const double * weight,double * wresp)
  {
  double el = exp(*linpred);
  double mu = el/(1+el);
  double pres = 0;
  unsigned i = 1;
  double u;
  if(*weight>0)
    {
    while(i<=*weight)
      {
      u = uniform();
      if(u <= mu)
        pres += 1;
      i += 1;
      }
    pres /= *weight;
    }

  *wresp = pres;
  }

//------------------------------------------------------------------------------
//-------------------- CLASS DISTRIBUTION_binomial_latent ----------------------
//------------------------------------------------------------------------------

void DISTRIBUTION_binomial_latent::tr_nonlinear(vector<double *> b,
                                                vector<double *> br,
                                                vector<FULLCOND*> & fcp,
                                                unsigned & nr,
                                                unsigned & it,
                                                ST::string & trtype)
  {
  if (trtype == "exp")
    DISTRIBUTION::tr_nonlinear(b,br,fcp,nr,it,trtype);
  else if ((trtype== "probit") || (trtype=="marginal"))
    {
    unsigned i;
    double h;
    double eh;
    for (i=0;i<b.size();i++)
      {
      h = *b[i]+interceptsample(it,0);
      if ((tlink) && (nu==8))
        {
        eh = exp(h);
        *br[i]= eh/(1+eh);
        }
      else
        {
        *br[i]=randnumbers::Phi2(h);
        }
      }

    }
  else if ( (trtype== "probitintercept") || (trtype=="marginalintercept"))
    {

    unsigned i;
    double eh;
    for(i=0;i<b.size();i++)
      {
      if ((tlink) && (nu==8))
        {
        eh = exp(interceptsample(it,0));
        *br[i]=eh/(1+eh);
        }
      else
        {
        *br[i] = randnumbers::Phi2(*b[i]);
        }

      }

    }

  }


void DISTRIBUTION_binomial_latent::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  if (!tlink)
    optionsp->out("  Response function: standard normal (probit link)\n");
  else
    {
    optionsp->out("  Response function: t-distribution function\n");
    optionsp->out("  Degrees of freedom: " + ST::inttostring(nu) + "\n");
    }

  optionsp->out("\n");
  optionsp->out("\n");
  }


// BEGIN: DSB //
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
double
    DISTRIBUTION_binomial_latent::sample_from_likelihood(const double weight,
                                                         const double mu) const
    {
        // determine the size parameter
        double size = (weight != 0.0) ? weight : 1.0;

        // use the sample function from Random.cpp
        return randnumbers::rand_binom(size, mu);
    }
#endif
// END: DSB //

void DISTRIBUTION_binomial_latent::compute_deviance(const double * response,
          const double * weight, const double * mu,double * deviance,
          double * deviancesat,
          const datamatrix & scale,const int & i) const
  {
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
  }


void DISTRIBUTION_binomial_latent::create(const bool & tl, const unsigned & n)
  {

  family = "Binomial (probit link)";
  scale(0,0) = 1;
  scaleexisting = false;

  tlink = tl;

  if (tlink==true)
    {
    changingweight = true;
    res = datamatrix(nrobs,1);
    nu = n;
    trmult = datamatrix(1,1,1.577);
    }

  unsigned i = 0;
  double * workresp=response.getV();
  double * workweight = weight.getV();
  bool error = false;

  while ( (i<nrobs) && (error==false )  )
    {

    if ( (*workresp != 0) && (*workresp != 1) )
      {
      error=true;
      errors.push_back("ERROR: response must be either zero or one\n");
      }

    if ( (*workweight != 0) && (*workweight != 1) )
      {
      error=true;
      errors.push_back("ERROR: weights must be either zero or one\n");
      }

    workresp++;
    workweight++;
    i++;
    }


  }


DISTRIBUTION_binomial_latent::DISTRIBUTION_binomial_latent(MCMCoptions * o,
                                     const datamatrix & r,
                                     const datamatrix & w,
                                     const bool & tl,
                                     const unsigned & n)
  : DISTRIBUTION(o,r,w)
  {

  create(tl,n);

  }


DISTRIBUTION_binomial_latent::DISTRIBUTION_binomial_latent(
                             const datamatrix & offset,MCMCoptions * o,
                             const datamatrix & r,const datamatrix & w,
                             const bool & tl,
                             const unsigned & n)
               : DISTRIBUTION(offset,o,r,w)
  {
  create(tl,n);
  }


DISTRIBUTION_binomial_latent::DISTRIBUTION_binomial_latent(
          const DISTRIBUTION_binomial_latent & nd)
   : DISTRIBUTION(DISTRIBUTION(nd))
  {
  res = nd.res;
  nu = nd.nu;
  tlink = nd.tlink;
  }


const DISTRIBUTION_binomial_latent &
DISTRIBUTION_binomial_latent::operator=(const DISTRIBUTION_binomial_latent & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));
  res = nd.res;
  nu = nd.nu;
  tlink = nd.tlink;
  return *this;
  }


double DISTRIBUTION_binomial_latent::loglikelihood(double * resp,
                                                   double * lin,
                                                   double * w,
                                                   const int & i) const
  {
  if (*w!=0)
    {
    double mu = randnumbers::Phi2(*lin);
    if (*resp > 0)
      return log(mu);
    else
      return log(1-mu);
    }
  else
    return 0;
  }



double DISTRIBUTION_binomial_latent::compute_weight(double * linpred,
                                             double * weight, const int & i,
                                             const unsigned & col) const
  {
  double  mu = randnumbers::Phi2(*linpred);
  double g =compute_gmu(linpred);
  g = g*g;
  return *weight/(mu*(1-mu)*g);
  }


double DISTRIBUTION_binomial_latent::compute_gmu(double * linpred,
                                                 const unsigned & col) const
  {
  double h = 0.39894228*exp(-0.5 * *linpred * *linpred);
  return 1.0/h;
  }


void DISTRIBUTION_binomial_latent::compute_mu(const double * linpred,
                                              double * mu) const
  {
  if ((tlink) && (nu==8))
    {
    double el = exp(trmult(0,0)* *linpred);
    *mu = el/(1+el);
    }
  else
    *mu = randnumbers::Phi2(*linpred);
  }


void DISTRIBUTION_binomial_latent::compute_mu_notransform(const double * linpred,
                                              double * mu) const
  {
  if ((tlink) && (nu==8))
    {
    double el = exp(*linpred);
    *mu = el/(1+el);
    }
  else
    *mu = randnumbers::Phi2(*linpred);
  }


void DISTRIBUTION_binomial_latent::outresults(void)
  {
  DISTRIBUTION::outresults();
  }


void DISTRIBUTION_binomial_latent::update(void)
  {

  double * worklin;
  double * workresp;
  double * weightwork;

  register unsigned i;

  if (!tlink)
    {

    worklin = (*linpred_current).getV();
    workresp = response.getV();
    weightwork = weight.getV();
    for(i=0;i<nrobs;i++,worklin++,workresp++,weightwork++)
      {

      if (*weightwork != 0)
        {
        if (*workresp > 0)
          *workresp = trunc_normal2(0,20,*worklin,1);
        else
          *workresp = trunc_normal2(-20,0,*worklin,1);
        }

      }

    }
  else   // tlink
    {

    worklin = (*linpred_current).getV();
    workresp = response.getV();

    res.minus(response,*linpred_current);

    double * reswork = res.getV();
    weightwork = weight.getV();

    double s;

    for(i=0;i<nrobs;i++,worklin++,workresp++,reswork++,
                       weightwork++)
      {

      if (*weightwork != 0)
        {
        *weightwork = 1.0/rand_invgamma(nu/2+0.5,nu/2+(*reswork * *reswork)/2);

        s = sqrt(1.0/(*weightwork));

        if (*workresp>0)
          *workresp =  trunc_normal2(0,20,*worklin,s);
        else
          *workresp =  trunc_normal2(-20,0,*worklin,s);
        }

      }

    }

  DISTRIBUTION::update();

  }




bool DISTRIBUTION_binomial_latent::posteriormode(void)
  {

  return true;

  }


bool DISTRIBUTION_binomial_latent::posteriormode_converged(const unsigned & itnr)
  {
  return true;
  }


void DISTRIBUTION_binomial_latent::compute_bootstrap_data(const double * linpred,const double * weight,double * wresp)
  {
  double mu = randnumbers::Phi2(*linpred);
  double pres = 0;
  unsigned i = 1;
  double u;
  if(*weight>0)
    {
    while(i<=*weight)
      {
      u = uniform();
      if(u <= mu)
        pres += 1;
      i += 1;
      }
    pres /= *weight;
    }

  *wresp = pres;
  }

//------------------------------------------------------------------------------
//----------------- CLASS DISTRIBUTION_binomial_logit_latent -------------------
//------------------------------------------------------------------------------


void DISTRIBUTION_binomial_logit_latent::tr_nonlinear(
                                        vector<double *> b,vector<double *> br,
                                        vector<FULLCOND*> & fcp, unsigned & nr,
                                        unsigned & it, ST::string & trtype)
  {
  if (trtype == "exp")
    DISTRIBUTION::tr_nonlinear(b,br,fcp,nr,it,trtype);
  else if ( (trtype == "logit") || (trtype == "marginal") )
    {
    double h;
    double eh;
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      h = *b[i]+interceptsample(it,0);
      eh = exp(h);
      *br[i] = eh/(1+eh);
      }
    }
  else if ( (trtype == "logitintercept") || (trtype== "marginalintercept") )
    {
    double eh;
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      eh = exp(interceptsample(it,0));
      *br[i] = eh/(1+eh);
      }
    }

  }


void DISTRIBUTION_binomial_logit_latent::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("  Response function: logistic distribution function\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


// BEGIN: DSB //
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
double
    DISTRIBUTION_binomial_logit_latent::sample_from_likelihood(const double weight,
                                                               const double mu) const
    {
        // determine the size parameter
        double size = (weight != 0.0) ? weight : 1.0;

        // use the sample function from Random.cpp
        return randnumbers::rand_binom(size, mu);
    }
#endif
// END: DSB //

void DISTRIBUTION_binomial_logit_latent::compute_deviance(
          const double * response, const double * weight,
          const double * mu,double * deviance, double * deviancesat,
          const datamatrix & scale,const int & i) const
  {
  if (*weight !=  0)
    {
    if (response_org(i,0)==0)
      {
      *deviance = -2*log(1-*mu);
      *deviancesat = *deviance;
      }
    else
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
  }


void DISTRIBUTION_binomial_logit_latent::create(const bool & tl)
  {
  family = "Binomial (logit link)";
  scale(0,0) = 1;
  scaleexisting = false;

  changingweight = true;
  res = datamatrix(nrobs,1);
  trmult = datamatrix(1,1,1);

  unsigned i = 0;
  double * workresp=response.getV();
  double * workweight = weight.getV();
  bool error = false;

  while ( (i<nrobs) && (error==false )  )
    {

    if ( (*workresp != 0) && (*workresp != 1) )
      {
      error=true;
      errors.push_back("ERROR: response must be either zero or one\n");
      }

    if ( (*workweight != 0) && (*workweight != 1) )
      {
      error=true;
      errors.push_back("ERROR: weights must be either zero or one\n");
      }

    workresp++;
    workweight++;
    i++;
    }

  response_org = response;

  }


DISTRIBUTION_binomial_logit_latent::DISTRIBUTION_binomial_logit_latent(
                                     MCMCoptions * o,
                                     const datamatrix & r,
                                     const datamatrix & w,
                                     const bool & tl)
  : DISTRIBUTION(o,r,w)
  {

  create(tl);

  }


DISTRIBUTION_binomial_logit_latent::DISTRIBUTION_binomial_logit_latent(
                             const datamatrix & offset,MCMCoptions * o,
                             const datamatrix & r,const datamatrix & w,
                             const bool & tl)
               : DISTRIBUTION(offset,o,r,w)
  {

  create(tl);

  }


DISTRIBUTION_binomial_logit_latent::DISTRIBUTION_binomial_logit_latent(
          const DISTRIBUTION_binomial_logit_latent & nd)
   : DISTRIBUTION(DISTRIBUTION(nd))
  {
  res = nd.res;
  response_org = nd.response_org;
  }


const DISTRIBUTION_binomial_logit_latent &
DISTRIBUTION_binomial_logit_latent::operator=(
const DISTRIBUTION_binomial_logit_latent & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));
  res = nd.res;
  response_org = nd.response_org;
  return *this;
  }


double DISTRIBUTION_binomial_logit_latent::loglikelihood(double * resp,
                                                   double * lin,
                                                   double * w,
                                                   const int & i) const
  {
  if (*w!=0)
    {

    if (*lin >= 10)
      return response_org(i,0) * *lin - *lin;
    else
      return response_org(i,0) * *lin - log(1+exp(*lin));

    }
  else
    return 0;
  }



double DISTRIBUTION_binomial_logit_latent::compute_weight(double * linpred,
                                             double * weight, const int & i,
                                             const unsigned & col) const
  {
  double el = exp(*linpred);
  double mu = el/(1+el);
  return mu*(1-mu);
  }


double DISTRIBUTION_binomial_logit_latent::compute_gmu(double * linpred,
                                                 const unsigned & col) const
  {
  double el = exp(*linpred);
  double mu = el/(1+el);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;
  return 1.0/(mu*(1-mu));
  }


void DISTRIBUTION_binomial_logit_latent::compute_mu(const double * linpred,
                                                    double * mu) const
  {

  double el = exp(*linpred);
  *mu = el/(1+el);
  }


void DISTRIBUTION_binomial_logit_latent::compute_mu_notransform(
                                                  const double * linpred,
                                                    double * mu) const
  {

  double el = exp(*linpred);
  *mu = el/(1+el);
  }


void DISTRIBUTION_binomial_logit_latent::outresults(void)
  {
  DISTRIBUTION::outresults();
  }


void DISTRIBUTION_binomial_logit_latent::update(void)
  {

  register unsigned i;

  double * worklin = (*linpred_current).getV();
  double * workresp = response.getV();
  double * workresp_org = response_org.getV();
  double * weightwork = weight.getV();

  double help;


  for(i=0;i<nrobs;i++,worklin++,workresp++,//reswork++,
                       weightwork++,workresp_org++)
    {

    if (*workresp_org==1)
      *workresp = randnumbers::trunc_logistic(*worklin, 0);

    else
      *workresp = randnumbers::trunc_logistic(*worklin, 1);


    help = (*workresp-*worklin)*(*workresp-*worklin);
    *weightwork =  1.0/randnumbers::lambda_fc(help);

    }

//  ofstream out("c:\\tmp\\weight.raw");
//  weight.prettyPrint(out);
//  out.close();


/*
  register unsigned i;

  double * worklin = (*linpred_current).getV();
  double * workresp = response.getV();
  double * workresp_org = response_org.getV();
  double * weightwork = weight.getV();
  double * reswork = res.getV();

  double help;
  double old;
  double prop;
  double u;
  double s;
  double alpha;

  for(i=0;i<nrobs;i++,worklin++,workresp++,reswork++,
                       weightwork++)
    {

    if (*weightwork != 0)
      {
      old = 1.0/(*weightwork);
      help = kssample();
      prop = 4*help*help;
      u = uniform();
//      nrtrials++;
      alpha = 0.5*(log(old)-log(prop)) + 0.5*(*reswork * *reswork)*
              (1.0/old-1.0/prop);

      if (log(u) <= alpha)
        {
        *weightwork = 1.0/prop;
        acceptancescale++;
//        (*aswork)++;
        }

      s = sqrt(1.0/(*weightwork));

      if (*workresp>0)
        *workresp =  trunc_normal2(0,20,*worklin,s);
      else
        *workresp =  trunc_normal2(-20,0,*worklin,s);
      }

    }
*/

  DISTRIBUTION::update();

  }




bool DISTRIBUTION_binomial_logit_latent::posteriormode(void)
  {

  return true;

  }


bool DISTRIBUTION_binomial_logit_latent::posteriormode_converged(
     const unsigned & itnr)
  {

  return true;

  }


//------------------------------------------------------------------------------
//----------------------- CLASS DISTRIBUTION_poisson ---------------------------
//------------------------------------------------------------------------------

DISTRIBUTION_poisson::DISTRIBUTION_poisson(MCMCoptions * o,
                                          const datamatrix & r,
                                          const datamatrix & w)

   : DISTRIBUTION(o,r,w)
     {

     family = "Poisson";
     scale(0,0) = 1;
     scaleexisting = false;

    }


DISTRIBUTION_poisson::DISTRIBUTION_poisson(const datamatrix & offset,
                      MCMCoptions * o, const datamatrix & r,
                      const datamatrix & w)
  : DISTRIBUTION(offset,o,r,w)
  {

  family = "Poisson";
  scale(0,0) = 1;
  scaleexisting = false;

  }


double DISTRIBUTION_poisson::compute_weight(double * linpred, double * weight,
                                            const int & i,
                                            const unsigned & col) const
  {

  return *weight * exp(*linpred);

  }

void DISTRIBUTION_poisson::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("  Response function: exponential\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


void DISTRIBUTION_poisson::update(void)
  {
  DISTRIBUTION::update();
  }


bool DISTRIBUTION_poisson::posteriormode(void)
  {
  return true;
  }


bool DISTRIBUTION_poisson::posteriormode_converged(const unsigned & itnr)
  {

  return true;

  }


double DISTRIBUTION_poisson::loglikelihood(double * response,double * linpred,
                                           double * weight,const int & i) const
  {
  return *weight * (*response * *linpred - exp(*linpred));
  }


void DISTRIBUTION_poisson::compute_mu(const double * linpred,double * mu) const
  {
  *mu = exp(*linpred);
  }


void DISTRIBUTION_poisson::compute_mu_notransform(
                       const double * linpred,double * mu) const
  {
  *mu = exp(*linpred);
  }


double DISTRIBUTION_poisson::compute_gmu(double * linpred,
                                         const unsigned & col) const
  {
  return 1.0/exp(*linpred);
  }


double DISTRIBUTION_poisson::compute_IWLS(double * response,double * linpred,
                            double * weight,const int & i,
                            double * weightiwls, double * tildey,
                            bool weightyes,
                            const unsigned & col)
  {
  double mu = exp(*linpred);

  if (weightyes)
    *weightiwls = *weight * mu;

  *tildey = (*response - mu)/mu;

  return *weight * (*response * *linpred - mu);

  }


void DISTRIBUTION_poisson::compute_IWLS_weight_tildey(double * response,
                      double * linpred, double * weight,
                      const int & i,double * weightiwls,
                      double * tildey, const unsigned & col)
  {

  double mu = exp(*linpred);

  *weightiwls = *weight * mu;

  *tildey =  (*response - mu)/mu;

  }

// BEGIN: DSB //
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
double
    DISTRIBUTION_poisson::sample_from_likelihood(const double weight,
                                                 const double mu) const
    {
        // derive the mean for the Poisson distribution
        double m = (weight != 0.0) ? weight * mu : mu;

        // then use the sample function from Random.cpp
        return randnumbers::rand_pois(m);
    }
#endif
// END: DSB //

void DISTRIBUTION_poisson::compute_deviance(const double * response,
                                            const double * weight,
                           const double * mu,double * deviance, double *
                           deviancesat,
                           const datamatrix & scale,const int & i) const
    {
    if (*response==0)
      {
      *deviance = 2* *weight * *mu;
      *deviancesat = *deviance;
      }
    else
      {
      *deviance = -2* *weight*(*response*log(*mu)-*mu);
      *deviancesat = *deviance+2 *
                     *weight*(*response * log(*response) - *response);
      }

    }


void DISTRIBUTION_poisson::tr_nonlinear(vector<double *> b,vector<double *> br,
                    vector<FULLCOND*> & fcp,unsigned & nr,
                    unsigned & it,ST::string & trtype)
  {

  if (trtype == "exp")
    DISTRIBUTION::tr_nonlinear(b,br,fcp,nr,it,trtype);
  else if (trtype == "elasticity")
    {
    if (b.size() == 2)
      {
      *br[1] = *b[1] * fcp[0]->get_data(nr,0);
      }
    }
  else if (trtype=="marginal")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0)+ *b[i]);
      }
    }
  else if (trtype=="marginalintercept")
    {
    unsigned i;
    for (i=0;i<b.size();i++)
      {
      *br[i] = exp(interceptsample(it,0));
      }
    }

  }


void DISTRIBUTION_poisson::compute_bootstrap_data(const double * linpred,const double * weight,double * wresp)
  {
  double mu = exp(*linpred);
  double pres = 0;
  if(*weight>0)
    {
    mu *= *weight;
    double zeit = 0;
    while(zeit<=1)
      {
      zeit += randnumbers::rand_expo(mu);
      pres += 1;
      }
    pres -= 1;
    pres / *weight;
    }
  *wresp = pres;
  }


//------------------------------------------------------------------------------
//------------------------ CLASS DISTRIBUTION_multinom -------------------------
//------------------------------------------------------------------------------

DISTRIBUTION_multinom::DISTRIBUTION_multinom(MCMCoptions * o,
                       const datamatrix & r,
                      const double & refvalue,
                      const datamatrix & w)

  : DISTRIBUTION(o,r,w)
  {
  muhelp = datamatrix(response.cols(),1);
  family = "Multinomial (logit link)";
  scale(0,0) = 1;
  scaleexisting = false;
  reference = ST::doubletostring(refvalue,6);
  }


DISTRIBUTION_multinom::DISTRIBUTION_multinom(const DISTRIBUTION_multinom & nd) :
                         DISTRIBUTION(DISTRIBUTION(nd))
  {
  reference = nd.reference;
  muhelp = nd.muhelp;
  }


const DISTRIBUTION_multinom & DISTRIBUTION_multinom::operator=(
                                         const DISTRIBUTION_multinom & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));
  reference = nd.reference;
  muhelp = nd.muhelp;
  return *this;
  }


void DISTRIBUTION_multinom::compute_mu(const double * linpred,double * mu) const
  {

  unsigned i;
  double sum = 1;

  const double * linpredstart = linpred;

  for(i=0;i<linearpred.cols();i++,linpred++)
    sum+= exp(*linpred);

  for(i=0;i<linearpred.cols();i++,linpredstart++,mu++)
    *mu = exp(*linpredstart)/sum;

  }


void DISTRIBUTION_multinom::compute_mu_notransform(
const double * linpred,double * mu) const
  {

  unsigned i;
  double sum = 1;

  const double * linpredstart = linpred;

  for(i=0;i<linearpred.cols();i++,linpred++)
    sum+= exp(*linpred);

  for(i=0;i<linearpred.cols();i++,linpredstart++,mu++)
    *mu = exp(*linpredstart)/sum;

  }


double DISTRIBUTION_multinom::compute_gmu(double * linpred,
                                          const unsigned & col) const
  {

  double * linpredstart = linpred;

  double sum = 1;

  unsigned i;
  for(i=0;i<linearpred.cols();i++,linpredstart++)
    sum+= exp(*linpredstart);

  linpredstart = linpred+col;

  double el = exp(*linpredstart);
  double mu = el/sum;
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;
  return 1.0/(mu*(1-mu));

  }


double DISTRIBUTION_multinom::compute_weight(double * linpred, double * weight,
                                             const int & i,
                                             const unsigned & col) const
  {

  register unsigned j;
  double expcol=0.0;
  double sumexp=0;
  for (j=0;j<linearpred.cols();j++,linpred++)
    {
    if (j==col)
      {
      expcol=exp(*linpred);
      sumexp+= expcol;
      }
    else
      {
      sumexp+= exp(*linpred);
      }

    }

  double mu = expcol/(1+sumexp);
  return mu*(1-mu);

  }



void DISTRIBUTION_multinom::compute_IWLS_weight_tildey(double * response,
                      double * linpred, double * weight,
                      const int & i,double * weightiwls,
                      double * tildey, const unsigned & col)
  {

  register unsigned j;
  double expcol=0.0;
  double sumexp=0;
  for (j=0;j<linearpred.cols();j++,linpred++)
    {
    if (j==col)
      {
      expcol=exp(*linpred);
      sumexp+= expcol;
      }
    else
      {
      sumexp+= exp(*linpred);
      }

    }

  double mu = expcol/(1+sumexp);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;

  *weightiwls = mu*(1-mu);

  double *respc = response+col;

  *tildey = (*respc - mu)/(*weightiwls);

  }


double DISTRIBUTION_multinom::compute_IWLS(double * resp,double * linpred,
                            double * weight,const int & i,
                            double * weightiwls, double * tildey,
                            bool weightyes,
                            const unsigned & col)
  {

  unsigned dim = response.cols();

  register unsigned j;
  double expcol=0.0;
  double sumexp=0;
  double * worklin = linpred;
  for (j=0;j<dim;j++,worklin++)
    {
    if (j==col)
      {
      expcol=exp(*worklin);
      sumexp+= expcol;
      }
    else
      {
      sumexp+= exp(*worklin);
      }

    }

  double mu = expcol/(1+sumexp);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;

   double v = mu*(1-mu);

  if (weightyes)
    *weightiwls = v;

  double * respc = resp+col;

  *tildey = (*respc - mu)/v;

  worklin = linpred;

  double logl = 0;
  double sum2=0;

  double * workresp = resp;

  for (j=0;j<dim;j++,workresp++,worklin++)
    {
    if ((*workresp) == 1)
      {
      sum2+= *workresp;
      logl +=  *worklin - log(1+sumexp);
      }
    }

  if (1-sum2 > 0)                                       // reference category
    {
    logl -= log(1+sumexp);
    }

  return logl;

  }


void DISTRIBUTION_multinom::compute_deviance(const double * response,
                                     const double * weight,
                                     const double * mu,double * deviance,
                                     double * deviancesat,
                                     const datamatrix & scale,const int & i)
                                     const
  {

  *deviance = 0;
  *deviancesat = 0;
  unsigned j=0;
  double sumy = 0;
  double summu = 0;
  for(j=0;j<linearpred.cols();j++,response++,mu++)
    {
    if (*response == 1)
      {
      sumy += 1;
      *deviance += log(*mu);
      }

    summu += *mu;

    }

  if ( 1-sumy > 0)
    {
    *deviance+= log(1-summu);
    }

  *deviance = -2* *deviance;
  *deviancesat = *deviance;

  }


void DISTRIBUTION_multinom::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("  Response function: logistic distribution function\n");
  optionsp->out("  Reference category: " + reference + "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTRIBUTION_multinom::loglikelihood(double * resp,double * linpred,
                                            double * weight,const int & i) const
  {

  unsigned j;
  unsigned dim = response.cols();

  double * worklin = linpred;

  double sum=0;
  for (j=0;j<dim;j++,worklin++)
    {
    sum+= exp(*worklin);
    }

  worklin = linpred;

  double logl = 0;
  double sum2=0;

  for (j=0;j<dim;j++,resp++,worklin++)
    {
    if ((*resp) == 1)
      {
      sum2+= 1;
      logl += *worklin - log(1+sum);
      }
    }

  if (sum2 == 0)                 // reference category
    {
    logl -= log(1+sum);
    }

  return logl;
  }


void DISTRIBUTION_multinom::update(void)
  {
  DISTRIBUTION::update();
  }


void DISTRIBUTION_multinom::compute_iwls(void)
  {

  register unsigned i,j;
  unsigned dim = response.cols();

  double * worklin = (*linpred_current).getV();

  double * workres = response.getV();
  double * ywork = tildey.getV();
  double * workweightiwls = weightiwls.getV();
  double mu;
  double * muhelpp;

  for (i=0;i<nrobs;i++)
    {
    muhelpp = muhelp.getV();
    compute_mu(worklin,muhelpp);
    for(j=0;j<dim;j++,worklin++,ywork++,workres++,workweightiwls++)
      {
      mu = muhelp(j,0);
      if(mu > 0.999)
        mu = 0.999;
      if(mu < 0.001)
        mu = 0.001;

      *workweightiwls = mu*(1-mu);
      *ywork = *worklin + (*workres - mu)/(*workweightiwls);
      }
    }
  }


bool DISTRIBUTION_multinom::posteriormode(void)
  {

  return true;

  }

bool DISTRIBUTION_multinom::posteriormode_converged(const unsigned & itnr)
  {
  return true;
  }


//------------------------------------------------------------------------------
//------------------------ CLASS DISTRIBUTION_multinom2 ------------------------
//------------------------------------------------------------------------------

void DISTRIBUTION_multinom2::create(void)
  {
  bool error=false;
  unsigned i=0;
  unsigned j;
  double * workresp = response.getV();
  double * workweight = weight.getV();
  while ( (i<nrobs) && (error==false) )
    {
    j=0;
    while ( (j<response.cols()) && (error==false) )
      {
      if (*workweight > 0)
        {
        if (*workresp != int(*workresp))
          {
          error=true;
          errors.push_back("ERROR: response cannot be multinomial; values must be integer numbers\n");
          }

        if (*workresp < 0)
          {
          error=true;
          errors.push_back("ERROR: response cannot be multinomial; some values are negative\n");
          }

        if (*workresp > *workweight)
          {
          error = true;
          errors.push_back("ERROR: response cannot be multinomial;\n");
          errors.push_back("       number of successes larger than number of trials for some values\n");
          }

        *workresp = *workresp/(*workweight);
        }
      workresp++;
      j++;
      }
    i++;
    workweight++;
    }
  }


DISTRIBUTION_multinom2::DISTRIBUTION_multinom2(MCMCoptions * o,
                       const datamatrix & r,
                      const double & refvalue,
                      const datamatrix & w)

  : DISTRIBUTION(o,r,w)
  {
  muhelp = datamatrix(response.cols(),1);
  family = "Multinomial (logit link)";
  scale(0,0) = 1;
  scaleexisting = false;
  reference = ST::doubletostring(refvalue,6);
  if(refvalue == -100)
    create();
  }


DISTRIBUTION_multinom2::DISTRIBUTION_multinom2(const DISTRIBUTION_multinom2 & nd) :
                         DISTRIBUTION(DISTRIBUTION(nd))
  {
  reference = nd.reference;
  muhelp = nd.muhelp;
  }


const DISTRIBUTION_multinom2 & DISTRIBUTION_multinom2::operator=(
                                         const DISTRIBUTION_multinom2 & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));
  reference = nd.reference;
  muhelp = nd.muhelp;
  return *this;
  }


void DISTRIBUTION_multinom2::compute_mu(const double * linpred,double * mu) const
  {

  unsigned i;
  double sum = 1;

  const double * linpredstart = linpred;

  for(i=0;i<linearpred.cols();i++,linpred++)
    sum+= exp(*linpred);

  for(i=0;i<linearpred.cols();i++,linpredstart++,mu++)
    *mu = exp(*linpredstart)/sum;

  }


void DISTRIBUTION_multinom2::compute_mu_notransform(
const double * linpred,double * mu) const
  {

  unsigned i;
  double sum = 1;

  const double * linpredstart = linpred;

  for(i=0;i<linearpred.cols();i++,linpred++)
    sum+= exp(*linpred);

  for(i=0;i<linearpred.cols();i++,linpredstart++,mu++)
    *mu = exp(*linpredstart)/sum;

  }


double DISTRIBUTION_multinom2::compute_gmu(double * linpred,
                                          const unsigned & col) const
  {

  double * linpredstart = linpred;

  double sum = 1;

  unsigned i;
  for(i=0;i<linearpred.cols();i++,linpredstart++)
    sum+= exp(*linpredstart);

  linpredstart = linpred+col;

  double el = exp(*linpredstart);
  double mu = el/sum;
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;
  return 1.0/(mu*(1-mu));

  }


double DISTRIBUTION_multinom2::compute_weight(double * linpred, double * weight,
                                             const int & i,
                                             const unsigned & col) const
  {

  register unsigned j;
  double expcol = 0.0;
  double sumexp=0;
  for (j=0;j<linearpred.cols();j++,linpred++)
    {
    if (j==col)
      {
      expcol=exp(*linpred);
      sumexp+= expcol;
      }
    else
      {
      sumexp+= exp(*linpred);
      }

    }

  double mu = expcol/(1+sumexp);
  return mu*(1-mu)* *weight;

  }



void DISTRIBUTION_multinom2::compute_IWLS_weight_tildey(double * response,
                      double * linpred, double * weight,
                      const int & i,double * weightiwls,
                      double * tildey, const unsigned & col)
  {

  register unsigned j;
  double expcol=0.0;
  double sumexp=0;
  for (j=0;j<linearpred.cols();j++,linpred++)
    {
    if (j==col)
      {
      expcol=exp(*linpred);
      sumexp+= expcol;
      }
    else
      {
      sumexp+= exp(*linpred);
      }

    }

  double mu = expcol/(1+sumexp);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;

  *weightiwls = mu*(1-mu);

  double *respc = response+col;

  *tildey = (*respc - mu)/(*weightiwls);
  *weightiwls = *weightiwls * *weight;
  }


double DISTRIBUTION_multinom2::compute_IWLS(double * resp,double * linpred,
                            double * weight,const int & i,
                            double * weightiwls, double * tildey,
                            bool weightyes,
                            const unsigned & col)
  {

  unsigned dim = response.cols();

  register unsigned j;
  double expcol=0.0;
  double sumexp=0;
  double * worklin = linpred;
  for (j=0;j<dim;j++,worklin++)
    {
    if (j==col)
      {
      expcol=exp(*worklin);
      sumexp+= expcol;
      }
    else
      {
      sumexp+= exp(*worklin);
      }

    }

  double mu = expcol/(1+sumexp);
  if(mu > 0.999)
    mu = 0.999;
  if(mu < 0.001)
    mu = 0.001;

   double v = mu*(1-mu);

  if (weightyes)
    *weightiwls = v * *weight;

  double * respc = resp+col;

  *tildey = (*respc - mu)/v;

  worklin = linpred;

  double logl = 0;
  double sum2=0;

  double * workresp = resp;

  for (j=0;j<dim;j++,workresp++,worklin++)
    {
    if ((*workresp) > 0)
      {
      sum2+= *workresp;
      logl +=  *workresp * *worklin;
      }
    }

  logl -= log(1+sumexp);

  return logl * *weight;
  }


void DISTRIBUTION_multinom2::compute_overall_deviance(double & deviance,double & deviancesat)
  {

  unsigned i;
  double * workresp=response.getV();
  double * worklin = (*linpred_current).getV();
  double * workweight = weight.getV();
  double dev=0;
  double devsat=0;
  //double mu;
  datamatrix mu = datamatrix(linearpred.cols(),1,0);
  double * mum = mu.getV();

  for(i=0;i<nrobs;i++,workresp++,worklin++,workweight++)
    {
    if (*workweight != 0)
      {
      compute_mu(worklin,mum);
      compute_deviance(workresp,workweight,mum,&dev,&devsat,scale,0);
      deviance+=dev;
      deviancesat+=devsat;
      worklin += linearpred.cols() - 1;
      workresp += linearpred.cols() - 1;
      }
    }
  }


void DISTRIBUTION_multinom2::compute_deviance(const double * response,
                                     const double * weight,
                                     const double * mu,double * deviance,
                                     double * deviancesat,
                                     const datamatrix & scale,const int & i)
                                     const
  {
  *deviance = 0;
  *deviancesat = 0;
  if(*weight > 0)
    {
    unsigned j=0;
    double sumy = 0;
    double summu = 0;
    for(j=0;j<linearpred.cols();j++,response++,mu++)
      {
      double muv = *mu;
      if(muv > 0.999)
        muv = 0.999;
      if(muv < 0.001)
        muv = 0.001;

      if (*response > 0)
        {
        sumy += *response;
        *deviance += *response * log(muv);
        *deviancesat += *response * log(*response);
        }
      summu += muv;
      }

    if(summu > 0.999)
      summu = 0.999;

    if ( 1-sumy > 0)
      {
      *deviance += (1-sumy) * log(1-summu);
      *deviancesat += (1-sumy) * log(1-sumy);
      }

    *deviance = -2 * *weight * *deviance;
    *deviancesat = *deviance + 2 * *weight * *deviancesat;
    }
  }


void DISTRIBUTION_multinom2::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("  Response function: logistic distribution function\n");
  optionsp->out("  Reference category: " + reference + "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTRIBUTION_multinom2::loglikelihood(double * resp,double * linpred,
                                            double * weight,const int & i) const
  {

  unsigned j;
  unsigned dim = response.cols();

  double * worklin = linpred;

  double sum=0;
  for (j=0;j<dim;j++,worklin++)
    {
    sum+= exp(*worklin);
    }

  worklin = linpred;

  double logl = 0;

  for (j=0;j<dim;j++,resp++,worklin++)
    {
    if ((*resp) > 1)
      {
      logl += *resp * *worklin;
      }
    }

  logl -= log(1+sum);

  return logl * *weight;
  }


void DISTRIBUTION_multinom2::update(void)
  {
  DISTRIBUTION::update();
  }


void DISTRIBUTION_multinom2::compute_iwls(void)
  {
  iwlsweights_notchanged_df = false;

  register unsigned i,j;
  unsigned dim = response.cols();

  double * worklin = (*linpred_current).getV();

  double * workres = response.getV();
  double * ywork = tildey.getV();
  double * workweightiwls = weightiwls.getV();
  double * wweight = weight.getV();
  double mu;
  double * muhelpp;

  for (i=0;i<nrobs;i++)
    {
    muhelpp = muhelp.getV();
    compute_mu(worklin,muhelpp);
    for(j=0;j<dim;j++,worklin++,ywork++,workres++,workweightiwls++)
      {
      mu = muhelp(j,0);
      if(mu > 0.999)
        mu = 0.999;
      if(mu < 0.001)
        mu = 0.001;

      *workweightiwls = mu*(1-mu);
      *ywork = *worklin + (*workres - mu)/(*workweightiwls);
      *workweightiwls = *workweightiwls * *wweight;
      }
    wweight++;
    }
  }


bool DISTRIBUTION_multinom2::posteriormode(void)
  {

  return true;

  }

bool DISTRIBUTION_multinom2::posteriormode_converged(const unsigned & itnr)
  {
  return true;
  }


void DISTRIBUTION_multinom2::compute_bootstrap_data(const double * linpred,const double * weight,double * wresp)
  {
  datamatrix mu = datamatrix(linearpred.cols(),1,0);
  double * wmu = mu.getV();
  double sum = 1;
  const double * linpredstart = linpred;

  unsigned i,j;
  for(i=0;i<linearpred.cols();i++,linpred++)
    sum+= exp(*linpred);
  linpred = linpredstart;
  for(i=0;i<linearpred.cols();i++,linpred++,wmu++)
    *wmu = exp(*linpred)/sum;

  double * wrespstart = wresp;
  for(i=0;i<linearpred.cols();i++,wresp++)
    *wresp = 0;
  wresp = wrespstart;
  wmu = mu.getV();

  double u;
  i = 1;
  if(*weight>0)
    {
    while(i<=*weight)
      {
      sum = 0;
      u = uniform();
      wresp = wrespstart;
      wmu = mu.getV();
      bool erledigt = false;
      for(j=0;j<linearpred.cols();j++,wresp++,wmu++)
        {
        sum += *wmu;
        if(u <= sum && !erledigt)
          {
          *wresp = *wresp + 1;
          erledigt = true;
          }
        }
      i += 1;
      }
    wresp = wrespstart;
    for(j=0;j<linearpred.cols();j++,wresp++)
      *wresp = *wresp / *weight;
    }
  wresp = wrespstart; // + linearpred.cols() - 1;
  linpred = linpredstart;
  }

//------------------------------------------------------------------------------
//------------------ CLASS DISTRIBUTION_multinomial_latent ---------------------
//------------------------------------------------------------------------------

void DISTRIBUTION_multinomial_latent::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  optionsp->out("  Response function: multivariate normal (independent probit)\n");
  optionsp->out("  Reference category: " + ST::doubletostring(refvalue,6) +
  "\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }


double DISTRIBUTION_multinomial_latent::maxutility(double * r,
const unsigned & cat)
  {
  unsigned j;
  double max = 0;
  double * rwork = r;
  for (j=0;j<nrcat;j++,rwork++)
    {
    if ( (j!=cat) && (*rwork > max) )
      max = *rwork;
    }

  return max;
  }


DISTRIBUTION_multinomial_latent::DISTRIBUTION_multinomial_latent(MCMCoptions * o,
                                     const datamatrix & r,
                                     const double & rv)
  : DISTRIBUTION(o,r)
  {
  family = "Multinomial (probit link)";
  scale(0,0) = 1;
  scaleexisting = false;

  unsigned i,j,l;

  posbeg.push_back(0);
  for(i=1;i<nrobs;i++)
    {
    if ( response(i,0) != response(i-1,0) )
      {
      posbeg.push_back(i);
      posend.push_back(i-1);
      }
    else if (i==nrobs-1)
      posend.push_back(i);
    }

  nrcat = posbeg.size()-1;

  if (nrcat == 0)
    errors.push_back("ERROR: response variable does not vary\n");

  if (nrcat > 10)
    errors.push_back("ERROR: too many values for the response variable\n");

  if (errors.size() == 0)
    {

    refvalue = rv;
    bool catfound = false;
    responsecat = datamatrix(nrcat+1,1);
    for(i=0;i<posbeg.size();i++)
      {
      responsecat(i,0) = response(posbeg[i],0);
      if (response(posbeg[i],0) == refvalue)
        {
        refcat = i;
        catfound = true;
        }
      }

    if (catfound == false)
      {
      refcat = 0;
      refvalue = response(posbeg[0],0);
      }

    response = datamatrix(nrobs,nrcat,0);

    linearpred = datamatrix(nrobs,nrcat,0);
    linearpredprop = linearpred;
    linpred_current = &linearpred;
    linpred_proposed = &linearpredprop;

    unsigned nr =0;

    for(j=0;j<nrcat+1;j++)
      {
      if (j==refcat)
        {
        for(i=posbeg[refcat];i<=posend[refcat];i++)
          {
          for (l=0;l<nrcat;l++)
            response(i,l) = trunc_normal(-20,0,0);
          }
        } // end: if (j==refcat)
      else
        {
        for(i=posbeg[j];i<=posend[j];i++)
          {
          for (l=0;l<nrcat;l++)
            {
            if (l==nr)
              {
              response(i,l) = trunc_normal(maxutility(&response(i,0),l),20,0);
              }
            else
              response(i,l) = trunc_normal(-20,response(i,nr),0);
            }

          }

        nr++;
        } // end: (j!=refcat)

      } // end: for(j=0;j<nrcat+1;j++)

    trmult = datamatrix(nrcat,1,1.0);

    }

  }


DISTRIBUTION_multinomial_latent::DISTRIBUTION_multinomial_latent(
                            const DISTRIBUTION_multinomial_latent & nd)
   : DISTRIBUTION(DISTRIBUTION(nd))
  {
  responsecat = nd.responsecat;
  refcat = nd.refcat;
  refvalue = nd.refvalue;
  nrcat = nd.nrcat;
  posbeg = nd.posbeg;
  posend = nd.posend;
  }



const DISTRIBUTION_multinomial_latent &
   DISTRIBUTION_multinomial_latent::operator=(
   const DISTRIBUTION_multinomial_latent & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));
  responsecat = nd.responsecat;
  refcat = nd.refcat;
  nrcat = nd.nrcat;
  refvalue = nd.refvalue;
  posbeg = nd.posbeg;
  posend = nd.posend;
  return *this;
  }



double DISTRIBUTION_multinomial_latent::loglikelihood(double * resp,
                                                      double * lin,
                                                      double * w,
                                                      const int & i) const
  {
  double help = *resp-*lin;
  return  - 0.5 * help * help;
  }


void DISTRIBUTION_multinomial_latent::compute_mu(const double * linpred,
                                      double * mu) const
    {
    // NOCH FALSCH
    unsigned i;
    for (i=0;i<linearpred.cols();i++,mu++)
      *mu=0;
    }


void DISTRIBUTION_multinomial_latent::compute_mu_notransform(
const double * linpred,double * mu) const
    {
    // NOCH FALSCH
    unsigned i;
    for (i=0;i<linearpred.cols();i++,mu++)
      *mu=0;
    }


void DISTRIBUTION_multinomial_latent::compute_deviance(const double *
                                               response,const double * weight,
                                     const double * mu, double * deviance,
                                     double * deviancesat,
                                     const datamatrix & scale,const int & i) const
    {
    // NOCH FALSCH
    *deviance = 0;
    *deviancesat = *deviance;
    }


double DISTRIBUTION_multinomial_latent::compute_weight(double * linpred,
                                             double * weight, const int & i,
                                             const unsigned & col) const
  {
  return 1.0;
  }


double DISTRIBUTION_multinomial_latent::compute_gmu(double * linpred,
                                                    const unsigned & col) const
  {
  return  1;
  }


void DISTRIBUTION_multinomial_latent::update(void)
  {

  unsigned i,j,l;
  unsigned nr =0;

//  int * workindex = index.getV();

  double * worklin = (*linpred_current).getV();
  double * respwork = response.getV();
  std::vector<unsigned>::iterator begwork = posbeg.begin();
  std::vector<unsigned>::iterator endwork = posend.begin();
  double linnr;

  for(j=0;j<nrcat+1;j++,++begwork,++endwork)
    {
    if (j==refcat)
      {

      for(i=*begwork;i<=*endwork;i++)
        {
        for (l=0;l<nrcat;l++,worklin++,respwork++)
          *respwork = *worklin+truncnormal(-20-*worklin,-*worklin);
        }
      } // end: if (j==refcat)
    else
      {
      for(i=*begwork;i<=*endwork;i++)
        {
        linnr = (*linpred_current)(i,nr);
        response(i,nr) = linnr+ truncnormal(maxutility(respwork,nr)-linnr,
                         20-linnr);

        for (l=0;l<nrcat;l++,worklin++,respwork++)
          {
          if (l!=nr)
            *respwork = *worklin+truncnormal(-20-*worklin,response(i,nr)-*worklin);
          }

        }

      nr++;
      } // end: (j!=refcat)

    } // end: for(j=0;j<nrcat+1;j++)

  DISTRIBUTION::update();

  }

void DISTRIBUTION_multinomial_latent::outresults(void)
  {
  DISTRIBUTION::outresults();
  }



//------------------------------------------------------------------------------
//------------------ CLASS DISTRIBUTION_cumulative_latent3 ---------------------
//------------------------------------------------------------------------------

void DISTRIBUTION_cumulative_latent3::update_utilities(void)
  {

  double * resp = response.getV();
  double * worklin = (*linpred_current).getV();
  std::vector<unsigned>::iterator begwork = posbeg.begin();
  std::vector<unsigned>::iterator endwork = posend.begin();
  double sqrtscale = sqrt(scale(0,0));
  double * workweight = weight.getV();

  unsigned l=2;
  if (posbeg.size()==4)
    l=3;

  unsigned i,j;
  for (i=0;i<=l;i++,++begwork,++endwork)
    {
    if (i==0)
      {
      for (j=*begwork;j<=*endwork;j++,resp++,worklin++,workweight++)
        if (*workweight != 0)
          {
          *resp = trunc_normal2(-20.0,0.0,*worklin,sqrtscale);
          }
      }
    else if (i==2)   // last category
      {
      for (j=*begwork;j<=*endwork;j++,resp++,worklin++,workweight++)
        if (*workweight != 0)
          *resp = trunc_normal2(1.0,20.0,*worklin,sqrtscale);
      }
    else if (i==3)  // missing values
      {
      for (j=*begwork;j<=*endwork;j++,resp++,worklin++,workweight++)
        if (*workweight != 0)
          *resp = *worklin + sqrtscale*rand_normal();
      }
    else
      {
      for (j=*begwork;j<=*endwork;j++,resp++,worklin++,workweight++)
        if (*workweight !=0)
          *resp = trunc_normal2(0.0,1.0,*worklin,sqrtscale);
      }

    } // end: for (i=0;i<=nrcat;i++,++begwork,++endwork)



  }

DISTRIBUTION_cumulative_latent3::DISTRIBUTION_cumulative_latent3(MCMCoptions * o,
                                      const datamatrix & r,const datamatrix & w,
                                      const double & a,const double & b,
                                      const ST::string & p,const ST::string & ps)
  : DISTRIBUTION(o,r,w,p,ps)
  {

  acceptancescale=100;
  family = "Multinomial with ordered categories (probit link)";
  scale(0,0) = 1;
  scaleexisting = true;

  a_invgamma = a;
  b_invgamma = b;

  unsigned i;

  sumweight=0;
  double * workweight = weight.getV();
  for(i=0;i<nrobs;i++,workweight++)
    {
    if ( (*workweight == 1) || (*workweight==0) )
      sumweight += *workweight;
    else
      errors.push_back("ERROR: weights must be either zero or one\n");
    }

  posbeg.push_back(0);
  for(i=1;i<nrobs;i++)
    {
    if ( response(i,0) != response(i-1,0) )
      {
      posbeg.push_back(i);
      posend.push_back(i-1);
      }
    else if (i==nrobs-1)
      posend.push_back(i);
    }

  unsigned nrcat = posbeg.size()-1;

  refvalue = response(nrobs-1,0);

  if (nrcat == 0)
    errors.push_back("ERROR: response variable does not vary\n");

  if (nrcat > 3)
    errors.push_back("ERROR: response variable must be three categorical\n");

  if (posbeg.size() == 4)
    {
    optionsp->out("\n");
    optionsp->out("WARNING: response has 4 categories.\n");
    optionsp->out("         BayesX can only estimate models with 3 categories.\n");
    optionsp->out("         The largest category is assumed to indicate missing response values.\n");
    }

  if (errors.size() == 0)
    {
    unsigned k;
    for (k=0;k<3;k++)
      {
      update_utilities();
      }
    }

  }


void DISTRIBUTION_cumulative_latent3::outoptions(void)
  {
  DISTRIBUTION::outoptions();
  int h = 0;
  if (posbeg.size() == 4)
    {
    h = posend[3]-posbeg[3];
    }
  optionsp->out("  Number of missing observations: " + ST::inttostring(h) + "\n");
  optionsp->out("  Response function: standard normal distribution function\n");
  optionsp->out("\n");
  optionsp->out("\n");
  }



void DISTRIBUTION_cumulative_latent3::set_predict_cum(const ST::string & path,
const ST::string & pathdev,
                       datamatrix * p, vector<ST::string> & Dn)
  {
  Dp = p;
  Dnames = Dn;
  predict = true;
  predictpath = path;
  deviancepath = pathdev;
  linpredmean = datamatrix(nrobs,1,0);
  mumean = datamatrix(nrobs,2,0);
  deviancemean = datamatrix(nrobs,1,0);
  deviancemean_sat = datamatrix(nrobs,1,0);
  deviance = datamatrix(optionsp->compute_samplesize(),2,0);
  }


DISTRIBUTION_cumulative_latent3::DISTRIBUTION_cumulative_latent3(
                            const DISTRIBUTION_cumulative_latent3 & nd)
   : DISTRIBUTION(DISTRIBUTION(nd))
  {
  refvalue = nd.refvalue;
  posbeg = nd.posbeg;
  posend = nd.posend;
  a_invgamma = nd.a_invgamma;
  b_invgamma = nd.b_invgamma;
  sumweight = nd.sumweight;
  }



const DISTRIBUTION_cumulative_latent3 &
   DISTRIBUTION_cumulative_latent3::operator=(
   const DISTRIBUTION_cumulative_latent3 & nd)
  {
  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));
  refvalue = nd.refvalue;
  posbeg = nd.posbeg;
  posend = nd.posend;
  a_invgamma = nd.a_invgamma;
  b_invgamma = nd.b_invgamma;
  sumweight = nd.sumweight;
  return *this;
  }


double DISTRIBUTION_cumulative_latent3::loglikelihood(double * resp,
                                                      double * lin,
                                                      double * w,
                                                      const int & i) const
  {
  // FEHLT
  return 0;
  }


void DISTRIBUTION_cumulative_latent3::compute_mu(const double * linpred,
double * mu) const
  {

  double help;
  *mu = randnumbers::Phi2(-(*linpred)/sqrt(scale(0,0)));
  help = *mu;
  mu++;
  *mu = randnumbers::Phi2((1.0-(*linpred))/sqrt(scale(0,0))) - help;
  }


void DISTRIBUTION_cumulative_latent3::compute_mu_notransform(
const double * linpred, double * mu) const
  {

  double help;
  *mu = randnumbers::Phi2(-*linpred);
  help = *mu;
  mu++;
  double scalemu = Scalesave.get_betamean(0,0);
  *mu = randnumbers::Phi2( 1.0/sqrt(scalemu)-(*linpred) ) - help;
  }


void DISTRIBUTION_cumulative_latent3::compute_deviance(const double * response,
                           const double * weight,const double * mu,
                           double * deviance,double * deviancesat,
                           const datamatrix & scale,const int & i) const
  {

  if (*weight !=0)
    {
    if (*response < 0)
      *deviance = -2*log(*mu);
    else if (*response >=0 && *response <= 1)
      {
      mu++;
      *deviance = -2*log(*mu);
      }
    else
      {
      double h = *mu;
      mu++;
      h+= *mu;
      *deviance = -2*log(1-h);
      }
    }
  else
    *deviance = 0;

  *deviancesat = *deviance;
  }


void DISTRIBUTION_cumulative_latent3::update(void)
  {

  register unsigned i;

  update_utilities();

  // updating the scaleparameter

  double nrobsd = sumweight;

  double sum = 0;
  double help;

  double * worklin = (*linpred_current).getV();
  double * resp = response.getV();
  double * workweight = weight.getV();

  for (i=0;i<nrobs;i++,worklin++,resp++,workweight++)
    {
    if (*workweight !=0)
      {
      help = *resp - *worklin;
      sum += help*help;
      }
    }

  scale(0,0) = rand_invgamma(a_invgamma+0.5*nrobsd, b_invgamma+0.5*sum);

  trmult(0,0) = 1.0/sqrt(scale(0,0));

  DISTRIBUTION::update();

  }




void DISTRIBUTION_cumulative_latent3::update_predict(void)
  {

  if (predict)
    {

    unsigned samplesize = optionsp->get_samplesize();

    if(
      (optionsp->get_nriter() > optionsp->get_burnin())
      &&
      ((optionsp->get_nriter()-optionsp->get_burnin()-1)
       % (optionsp->get_step()) == 0)
      )
      {

      register unsigned i,j;
      double * worklin = (*linpred_current).getV();
      double * workmean = linpredmean.getV();
      double * workmumean = mumean.getV();
      datamatrix muhelp(mumean.cols(),1,0);

      double * workdevmean = deviancemean.getV();
      double * workdevmean_sat = deviancemean_sat.getV();
      double * workresponse = response.getV();
      double * workw = weight.getV();

      unsigned size2 = mumean.cols();

      double * mumeanhelp;
      double reshelp;
      double devhelp;

      double * musavep = musave.getbetapointer();

      if (samplesize==1)
        {
        for (i=0;i<nrobs;i++,worklin++,workmean++,workmumean+=size2
               ,workdevmean++,workdevmean_sat++,workw++,workresponse++)
          {
          mumeanhelp = workmumean;
          compute_mu(worklin,mumeanhelp);
          compute_deviance(workresponse,workw,workmumean,workdevmean,&devhelp,scale,i);
          deviance(0,0) += *workdevmean;
          deviance(0,1) += devhelp;
          *workdevmean_sat = devhelp;

          *workmean = *worklin/sqrt(scale(0,0));

          if ((predictfull) && (i<firstobs))
            {
            mumeanhelp = workmumean;
            for(j=0;j<size2;j++,musavep++,mumeanhelp++)
              {
              *musavep = *mumeanhelp;
              }
            }

          }

        }
        else
          {
          for (i=0;i<nrobs;i++,workmean++,worklin++,workdevmean++,
                           workdevmean_sat++,workw++,workresponse++)
            {

            mumeanhelp = muhelp.getV();
            compute_mu(worklin,mumeanhelp);
            mumeanhelp = muhelp.getV();

            compute_deviance(workresponse,workw,mumeanhelp,&reshelp,&devhelp,scale,i);
            deviance(samplesize-1,0) += reshelp;
            deviance(samplesize-1,1) += devhelp;


            *workmean = (1.0)/samplesize * ( (samplesize-1)* *workmean +
                           (*worklin)/sqrt(scale(0,0)));

            mumeanhelp = muhelp.getV();

            for(j=0;j<size2;j++,workmumean++,mumeanhelp++)
              {
              *workmumean = (1.0)/samplesize * ( (samplesize-1)* *workmumean +
                           *mumeanhelp);
              }

            *workdevmean = (1.0)/samplesize * ( (samplesize-1)* *workdevmean +
                            reshelp);

            *workdevmean_sat = (1.0)/samplesize * ( (samplesize-1)*
                               *workdevmean_sat + devhelp);

            if ((predictfull) && (i<firstobs))
              {
              mumeanhelp = muhelp.getV();
              for(j=0;j<size2;j++,musavep++,mumeanhelp++)
                {
                *musavep = *mumeanhelp;
                }
              }  // end: if (predictfull)

            }

          }

        }

     if (predictfull)
       musave.update();

    } // end: if predict

  }


void DISTRIBUTION_cumulative_latent3::outresults(void)
  {

  DISTRIBUTION::outresults();

  // threshold parameters

  datamatrix scalesample;
  datamatrix thetasample;

  scalesample = datamatrix(optionsp->compute_samplesize(),1);
  Scalesave.readsample(scalesample,0);

  thetasample = datamatrix(optionsp->compute_samplesize(),2);

  double * worksample = interceptsample.getV();
  double * workscale = scalesample.getV();
  unsigned i;
  for(i=0;i<interceptsample.rows();i++,worksample++,workscale++)
    {
    thetasample(i,0) = -*worksample;
    thetasample(i,1) = 1.0/(sqrt(*workscale))-*worksample;
    }


  double lower1 = Scalesave.get_lower1();
  double upper2 = Scalesave.get_upper2();

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string u2 = ST::doubletostring(upper2,4);



  optionsp->out("  Threshold parameters:\n");
  optionsp->out("\n");

  optionsp->out("            mean           Std. Dev.      " +
                  l1 + "% quant.     median         " + u2 + "% quant.\n");


  optionsp->out(ST::outresults(3,"theta_1",thetasample.mean(0),
                               sqrt(thetasample.var(0)),
                               thetasample.quantile(lower1,0),
                               thetasample.quantile(50,0),
                               thetasample.quantile(upper2,0)) + "\n");
  optionsp->out(ST::outresults(3,"theta_2",thetasample.mean(1),
                               sqrt(thetasample.var(1)),
                               thetasample.quantile(lower1,1),
                               thetasample.quantile(50,1),
                               thetasample.quantile(upper2,1)) + "\n");

  optionsp->out("\n");

  }

} // end: namespace MCMC





