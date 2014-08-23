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



#include "first.h"

#include "fullcond_projection.h"

namespace MCMC
{

using std::ios;

  // CONSTRUCTOR 1  (for additive models)

FULLCOND_projection::FULLCOND_projection(MCMCoptions * o,
                      DISTRIBUTION * dp, FULLCOND_const * fcc, const datamatrix & d,
                      const unsigned & nrk, const unsigned & degr, const knotpos & kp,
                      const fieldtype & ft, const ST::string & monotone, const ST::string & ti,
                      const ST::string & fp, const ST::string & pres, const bool & deriv,
                      const int & gs, vector<FULLCOND_projection*> & zeiger, const unsigned & nterms,
                      const unsigned & c)
  : spline_basis(o,dp,fcc,ft,ti,nrk,degr,kp,gs,fp,pres,deriv,c)
  {

  lambda = 100;                                   // noch verbessern!!!

  nrterms = nterms;
  nrvar = d.cols();  // + 1;
//  double ppw = 1/sqrt(double(nrvar));
//  pp_weights = datamatrix(nrvar,1,ppw);
  int ppw = zeiger.size()%nrvar; // + 1;
  pp_weights = datamatrix(nrvar,1,0);
  pp_weights(ppw,0) = 1;

/*  original_data = datamatrix(d.rows(),nrvar,1);
  for(unsigned i=1;i<nrvar;i++)
    {
    original_data.putCol(i, d.getCol(i-1));
    }      */
  original_data = d;
  data_forfixed = datamatrix(d.rows(),1,0);
  pp_pointer = zeiger;
  gesamt = datamatrix(d.rows(),1,0);

  if(monotone == "increasing")
    increasing = true;
  else if(monotone == "decreasing")
    decreasing = true;

  varcoeff = false;

  transform = likep->get_trmult(c);
  splinederivative = datamatrix(data_forfixed.rows(),1,0);

  compute_betaweight();

  compute_Kweights();

  // Penalty Matrix erstellen

  if (type == RW1)
    {
    K = Krw1band(weight);
    Kenv = Krw1env(weight);
    rankK = nrpar-nrparpredictleft-nrparpredictright-1;
    }
  else if (type == RW2)
    {
    K = Krw2band(weight);
    Kenv = Krw2env(weight);
    rankK = nrpar-nrparpredictleft-nrparpredictright-2;
    }

  XX_env = envmatdouble(bandmatdouble(nrpar,degree,0));
  if (type==RW1)
    prec_env = envmatdouble(0.0,nrpar,degree>1?degree:1);
  else if (type == RW2)
    prec_env = envmatdouble(0.0,nrpar,degree>2?degree:2);

  mu = datamatrix(likep->get_nrobs(),1,0);
  muy = datamatrix(nrpar,1,0);
  betahelp = muy;

  identifiable = false;

  likep->compute_iwls();
  compute_linear_combination(true);
  spline = datamatrix(spline.rows(),1,0);
  if(pp_pointer.size() == nrterms-1)     // beim letzten Term wird der lineare Prädiktor wieder auf Null gesetzt.
    likep->substr_linearpred_m(likep->get_linearpred(true),column,true);
  }


  // COPY CONSTRUCTOR

FULLCOND_projection::FULLCOND_projection(const FULLCOND_projection & fc)
  : spline_basis(spline_basis(fc))
  {
  nrterms = fc.nrterms;
  nrvar = fc.nrvar;
  pp_weights = fc.pp_weights;
  original_data = fc.original_data;
  pp_pointer = fc.pp_pointer;
  gesamt = fc.gesamt;
  Bderiv = fc.Bderiv;

  df_lambdaold = fc.df_lambdaold;
  lambdaold = fc.lambdaold;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const FULLCOND_projection & FULLCOND_projection::operator=(
                                            const FULLCOND_projection & fc)
  {
  if (this == &fc)
    return *this;
  spline_basis::operator=(spline_basis(fc));

  nrterms = fc.nrterms;
  nrvar = fc.nrvar;
  pp_weights = fc.pp_weights;
  original_data = fc.original_data;
  pp_pointer = fc.pp_pointer;
  gesamt = fc.gesamt;
  Bderiv = fc.Bderiv;

  df_lambdaold = fc.df_lambdaold;
  lambdaold = fc.lambdaold;

  return *this;
  }


void FULLCOND_projection::compute_linear_combination(bool eins)
  {
  unsigned i;
  //compute_Kweights();                           // wäre nötig für nicht-äquidistante Knoten
  bool abbruch = false;
  datamatrix pp_weights_alt = pp_weights;
datamatrix pp_weights_alt2 = datamatrix(pp_weights.rows(),1,0);
datamatrix pp_weights_alt3 = datamatrix(pp_weights.rows(),1,0);
double krit1 = 0;
double krit2 = 0;

  //datamatrix splinederivative_sort = datamatrix(data_forfixed.rows(),1,0);
  datamatrix splinederiv_square = datamatrix(data_forfixed.rows(),1,0);
  datamatrix ziel = datamatrix(data_forfixed.rows(),1,0);
  datamatrix X = datamatrix(nrvar,nrvar,0);

ST::string wfile = "h:\\simulationen\\sim_ppr\\adaptiv\\" + likep->get_responsename() + "_weights.txt";
outw.open(wfile.strtochar(),ios::app); // append = anfügen
outw << "lambda   w1   w2" << endl;

  unsigned z = 0;
  while(!abbruch)
    {
    if((z==0 && eins==true) || z>0)
      {
      data_forfixed.mult(original_data,pp_weights);   // berechnet Linearkombination

      // 1) Spline + Ableitung berechnen
      freq.erase(freq.begin(),freq.end());
      knot.erase(knot.begin(),knot.end());
      firstnonzero.erase(firstnonzero.begin(),firstnonzero.end());
      lastnonzero.erase(lastnonzero.begin(),lastnonzero.end());
      index = statmatrix<int>(1,1,0);
      begcol.erase(begcol.begin(),begcol.end());
      make_index(data_forfixed);
      make_Bspline(data_forfixed,true);     // eigene Funktion, die Matrix für Spline und Ableitung erstellt
      // index2 initialisieren
      index2.erase(index2.begin(),index2.end());
      index2.push_back(index(0,0));
      for(i=1;i<likep->get_nrobs();i++)
        index2.push_back(index(i,0)-index(i-1,0));
      }

    compute_XWXenv(likep->get_weightiwls(),column);
    prec_env.addto(XX_env,Kenv,1.0,lambda);

    likep->substr_linearpred_m(spline,column,true);
    likep->compute_workingresiduals(column);
    compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);
    prec_env.solve(muy,beta);

    add_linearpred_multBS(true);

    // 2) Gewichte berechnen aus KQ-Schätzung
    double * zielwork = ziel.getV();                                     // enthält Zielvariable für KQ-Schätzung
    double * squarework = splinederiv_square.getV();                     // enthält Gewichte für KQ-Schätzung
    double * rd = data_forfixed.getV();
    double * fdwork = splinederivative.getV();
    likep->compute_workingresiduals(column);                             // berechnet (y-spline)
    double * works = likep->get_workingresiduals().getV();
    double * workweight = likep->get_weightiwls().getV();                // für gewichtete KQ-Schätzung
    for(i=0;i<data_forfixed.rows();i++,zielwork++,squarework++,rd++,fdwork++,works++,workweight++)
      {
      *squarework = *fdwork * *fdwork * *workweight;
      *zielwork = *works / *fdwork + *rd;
      }

    // computes X = (data' W data)
    unsigned p,k;
    squarework = splinederiv_square.getV();
    double * workXp;
    double * workXk;
    for (p=0;p<nrvar;p++)
      for (k=p;k<nrvar;k++)
        {
        X(p,k)=0;
        squarework = splinederiv_square.getV();
        workXp = original_data.getV()+p;
        workXk = original_data.getV()+k;
        for(i=0;i<data_forfixed.rows();i++,squarework++,workXp+=nrvar,workXk+=nrvar)
          X(p,k)+= *squarework  *  *workXp * *workXk;
        X(k,p) = X(p,k);
        }

    X.assign(X.cinverse());

    // berechnet W * ziel
    squarework = splinederiv_square.getV();
    zielwork = ziel.getV();
    datamatrix workingres = datamatrix(data_forfixed.rows(),1,0);
    double * workres = workingres.getV();
    for(i=0;i<data_forfixed.rows();i++,squarework++,zielwork++,workres++)
      *workres = *squarework * *zielwork;

    pp_weights = X*original_data.transposed()*workingres;

/*    // Gewichte so normieren, dass Norm = 1;
    works = pp_weights.getV();
    double sum = 0;
    for(i=0;i<pp_weights.rows();i++,works++)
      sum = sum + *works * *works;
    works = pp_weights.getV();
    sum = sqrt(sum);
    for(i=0;i<pp_weights.rows();i++,works++)
      *works = *works/sum;    */

    // Gewichte so normieren, dass w1 = 1;
    works = pp_weights.getV();
    double sum = pp_weights(0,0);
    for(i=0;i<pp_weights.rows();i++,works++)
      *works = *works/sum;

     // Abbruchkriterium
    double diffmean = norm(pp_weights-pp_weights_alt);
    z += 1;
//double test = pp_weights(0,0);
//test = pp_weights(1,0);

krit1 = krit2;
krit2 = likep->compute_rss();
//test += Kenv.compute_quadform(beta,0);
//test = test;

    //if(diffmean < 0.00001)
    if(diffmean < 0.0001)
      {
      abbruch = true;
      pp_weights = pp_weights_alt;
      outw << ST::doubletostring(lambda,6) << "   "
           << ST::doubletostring(pp_weights_alt(0,0),6) << "   " << ST::doubletostring(pp_weights_alt(1,0),6) << endl;
      outw.close();
      }
    else
      {
/*      if(z==80 || z==160 || z==240)
        {
        int t = z/80;
        int ppw = pp_pointer.size()%nrvar;
        pp_weights = datamatrix(nrvar,1,0);
        if(ppw+t < int(nrvar))
          pp_weights(ppw+1,0) = 1;
        else if(ppw-t >= 0)
          pp_weights(ppw-1,0) = 1;
        else
          pp_weights = datamatrix(nrvar,1,1/sqrt(double(nrvar)));
        }       */
      //if(z==320)
      //  {
      //  abbruch = true;
      //  // Fehler angeben
      //  }
      if(norm(pp_weights-pp_weights_alt2)<0.0001 && norm(pp_weights_alt-pp_weights_alt3)<0.0001 && krit2<=krit1)
        {
        abbruch = true;

      outw << "n " << ST::doubletostring(lambda,6) << "   "
           << ST::doubletostring(pp_weights_alt(0,0),6) << "   " << ST::doubletostring(pp_weights_alt(1,0),6)
           << "   " << ST::doubletostring(pp_weights(0,0),6) << "   " << ST::doubletostring(pp_weights(1,0),6) << endl;
      outw.close();

        pp_weights = pp_weights_alt;
        }
      else
        {
        pp_weights_alt3 = pp_weights_alt2;
        pp_weights_alt2 = pp_weights_alt;
        pp_weights_alt = pp_weights;
        }
      }
    }

  init_fchelp(data_forfixed);
  compute_XWXenv(likep->get_weight());
  }


bool FULLCOND_projection::posteriormode(void)
  {

compute_linear_combination(false);     // *** Versuch ***

  unsigned i;
  transform = likep->get_trmult(column);
  fchelp.set_transform(transform);

/*  likep->substr_linearpred_m(spline,column,true);

  if ( (lambda_prec != lambda) || (likep->iwlsweights_constant() == false) )
    {
    if (likep->iwlsweights_constant() == false)
      {
      compute_XWXenv(likep->get_weightiwls(),column);
      }
    prec_env.addto(XX_env,Kenv,1.0,lambda);
    lambda_prec = lambda;
    }
  likep->compute_workingresiduals(column);
  compute_XWtildey(likep->get_weightiwls(),likep->get_workingresiduals(),1.0,column);

  prec_env.solve(muy,beta);

// monotone Regression!
  if(decreasing)
    {
    bool ok = false;
    while(!ok)
      {
      bool ok2 = true;
      for(unsigned i=1;i<nrpar;i++)
        {
        double diff = beta(i,0)-beta(i-1,0);
        if(diff > 0.0001)
          {
          ok2 = false;
          double mean = 0.5*( beta(i-1,0)+beta(i,0) );
          beta(i-1,0) = mean;
          beta(i,0) = mean;
          }
        if(diff > 0)
          {
          double help = beta(i,0);
          beta(i,0) = beta(i-1,0);
          beta(i-1,0) = help;
          }
        }
      ok = ok2;
      }
    beta.sortcol(0,nrpar-1,0);
    datamatrix bsort = beta;
    for(unsigned j=0;j<nrpar;j++)
      beta(j,0) = bsort(nrpar-1-j,0);
    }

  if(increasing)
    {
    bool ok = false;
    while(!ok)
      {
      bool ok2 = true;
      for(unsigned i=1;i<nrpar;i++)
        {
        double diff = beta(i-1,0)-beta(i,0);
        if(diff > 0.0001)
          {
          ok2 = false;
          double mean = 0.5*(beta(i-1,0)+beta(i,0));
          beta(i-1,0) = mean;
          beta(i,0) = mean;
          }
        if(diff > 0)
          {
          double help = beta(i,0);
          beta(i,0) = beta(i-1,0);
          beta(i-1,0) = help;
          }
        }
      ok = ok2;
      }
    beta.sortcol(0,nrpar-1,0);
    }
// ENDE: monoton

  add_linearpred_multBS(false);      */

  if(center)
    {
    compute_intercept();
    fcconst->posteriormode_intercept(intercept);
    }

  if(center)
    {
    int * workindex = index.getV();
    for(i=0;i<spline.rows();i++,workindex++)
      spline(*workindex,0) -= intercept;
    }
  double * fchelpbetap = fchelp.getbetapointer();

  if(gridsize < 0)
    {
    vector<int>::iterator freqwork = freqoutput.begin();
    int * workindex = index.getV();
    for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
      {
      if(freqwork==freqoutput.begin() || *freqwork!=*(freqwork-1))
        {
        *fchelpbetap = spline(*workindex,0);
        fchelpbetap++;
        }
      }
    }
  else
    {
    multDG(splinehelp,beta);
    for(i=0;i<gridsize;i++,fchelpbetap++)
      *fchelpbetap = splinehelp(i,0) - intercept;
    }

  intercept = 0.0;

  /*write_derivative();
  if(derivative)
    fcderivative.posteriormode();*/

  fchelp.posteriormode();
  return FULLCOND_nonp_basis::posteriormode();
  }


void FULLCOND_projection::create_weight(datamatrix & w)
  {
  vector<int>::iterator freqwork = freqoutput.begin();
  int * workindex = index.getV();
  unsigned i;
  for(i=0;i<likep->get_nrobs();i++,freqwork++,workindex++)
    {
    if(freqwork==freqoutput.begin() || (*freqwork!=*(freqwork-1) && *freqwork == *(freqoutput.end()-1) ))
      w(*workindex,0) = 1;
    }
  }


void FULLCOND_projection::make_Bspline(const datamatrix & md, const bool & deriv)
  {

  unsigned i,j,k;
  double value;
  double * work;
  double * workd;

  vector<int>::iterator freqwork;
  datamatrix help;

  double min = md(index(0,0),0);
  double max = md(index(md.rows()-1,0),0);
  double dist = max-min;

  min -= 0.01*dist;
  max += 0.01*dist;

// Knoten berechnen
  if(knpos == equidistant)
    {
    dist = (max - min)/(nrknots-1);
    knot.push_back(min - degree*dist);
    for(i=1;i<nrknots+2*degree;i++)
      knot.push_back(knot[i-1] + dist);
    }
  else if(knpos == quantiles)
    {
    double distfirst, distlast;
/*
    knot.push_back(min);
    for(i=1;i<nrknots-1;i++)
      knot.push_back(md.quantile((i*100)/double(nrknots-1),0));
    knot.push_back(max);
*/

    knot.push_back(min);
    for(i=1;i<nrknots-1;i++)
      {
      double q1 = md.quantile((i*100)/double(nrknots-1),0);
      double q2 = md.quantile(((i-1)*100)/double(nrknots-1),0);
      if(q1 != q2)
        knot.push_back(q1);
      }
    knot.push_back(max);

    nrknots = knot.size();

    distfirst = knot[1] - knot[0];
    distlast = knot[nrknots-1] - knot[nrknots-2];

    for(i=0;i<degree;i++)
      {
      knot.push_front(min - (i+1)*distfirst);
      knot.push_back(max + (i+1)*distlast);
      }
    }

  if(knpos==quantiles)
    {
    if(nrpar > nrknots-1+degree)
      {
      optionsp->outerror("\n");
      optionsp->outerror("WARNING: Reducing the number of basis functions for term " + title + "\n");
      optionsp->outerror("         due to equal quantiles for the knot positions.\n");
      }
    }

  setbeta(nrknots-1+degree,1,0);

// Designmatrix BS bzw. B (bei VCM), lastnonzero, firstnonzero und Bcolmean berechnen
  help = datamatrix(nrpar,1,0.0);
  Bcolmean = datamatrix(nrpar,1,0.0);

  for(i=0;i<nrpar;i++)
    {
    lastnonzero.push_back(-1);
    firstnonzero.push_back(0);
    }

  if(varcoeff)
    {
    B = datamatrix(*(freq.end()-1)+1,degree+1,0.0);
    work = B.getV();
    if(deriv == true)
      {
      Bderiv = datamatrix(*(freq.end()-1)+1,degree+1,0.0);
      workd = Bderiv.getV();
      }
    }
  else
    {
    BS = datamatrix(nrdiffobs,degree+1,0.0);
    work = BS.getV();
    if(deriv == true)
      {
      Bderiv = datamatrix(nrdiffobs,degree+1,0.0);
      workd = Bderiv.getV();
      }
    }

  freqwork = freq.begin();
  for(i=0;i<md.rows();i++,++freqwork)
//  for(freqwork=freq.begin();freqwork<freq.end();++freqwork)
    {
    value = md(index(i,0),0);
//    value = md(index(*freqwork,0),0);
    if(freqwork == freq.begin() || *freqwork != *(freqwork-1))
      {
      j=0;
      while(knot[degree+j+1] <= value)
        j++;
      begcol.push_back(j);

      // Werte der Designmatrix bestimmen
      help.assign(bspline(value));
      for(k=0;k<degree+1;k++,work++)
        {
        *work = help(k+j,0);
        Bcolmean(k+j,0) += *work;
        }

      if(deriv == true)
        {
        // Werte der Ableitungsmatrix bestimmen
        help.assign(bspline_derivative(value));
        for(k=0;k<degree+1;k++,workd++)
          {
          *workd = help(k+j,0);
          }
        }
      }

    for(k=j;k<nrpar;k++)
      lastnonzero[k] += 1;
    for(k=j+degree+1;k<nrpar;k++)
      firstnonzero[k] += 1;
    }

  for(i=0;i<nrpar;i++)
    Bcolmean(i,0) /= double(nrdiffobs);
  }


datamatrix FULLCOND_projection::bspline_derivative(const double & x)
  {

  datamatrix b(nrpar,1,0.0);
  datamatrix help(nrpar+degree,1,0.0);
  unsigned j;
  double * bwork;
  double * helpwork;

// Grad 0

  for(j=0;j<nrpar;j++)
    if( knot[j]<=x && x<knot[j+1])
      b(j,0) = 1.0;

  for(unsigned l=1;l<=degree-1;l++)
    {
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
      *helpwork = *bwork;
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
      {
      *bwork = (x-knot[j])**helpwork/(knot[j+l]-knot[j])
                  + (knot[j+l+1]-x)**(helpwork+1)/(knot[j+l+1]-knot[j+1]);
      }
    }

// Hämmerlin/Hoffmann Seite 263

  bwork = b.getV();
  helpwork = help.getV();
  for(j=0;j<nrpar;j++,++helpwork,++bwork)
    *helpwork = *bwork;
  bwork = b.getV();
  helpwork = help.getV();
  for(j=0;j<nrpar;j++,++helpwork,++bwork)
    *bwork = degree*( *helpwork/(knot[j+degree]-knot[j]) - *(helpwork+1)/(knot[j+degree+1]-knot[j+1]) );

  return b;
  }


void FULLCOND_projection::add_linearpred_multBS(const bool & deriv, const bool & current)
  {

  double *workBS;
  double * workBd;
  double *workbeta;
  double *lp;

  unsigned j,k;
  unsigned col = degree+1;
  unsigned lpcols = likep->get_linearpred(current).cols();
  int i,stop;

  vector<int>::iterator freqwork = freq.begin();
  vector<int>::iterator workindex2 = index2.begin();
// spline = 0 setzen
  double * workspline = spline.getV();
  double * workderiv = splinederivative.getV();
  for(j=0;j<spline.rows();j++,workspline++,workderiv++)
    {
    *workspline = 0.0;
    *workderiv = 0.0;
    }

  lp = likep->get_linearpred(current).getV() + column + *workindex2*lpcols;
  workspline = spline.getV() + *workindex2;
  if(deriv == true)
    workderiv = splinederivative.getV() + *workindex2;
// Spaltenweise multiplizieren
  i = 0;
  k = 0;
  workBS = BS.getV();
  if(deriv == true)
    workBd = Bderiv.getV();
  while (k<nrpar)
    {
    stop = lastnonzero[k];
//    while (i<lastnonzero[k]+1)
    while (i <= stop)
      {
      workbeta = beta.getV() + k;
      for(j=0;j<col;j++,workBS++,workbeta++,workBd++)
        {
        *lp += *workBS * *workbeta;
        *workspline += *workBS * *workbeta;
        if(deriv == true)
          *workderiv += *workBd * * workbeta;
        }
      if((freqwork+1)!=freq.end() && *freqwork==*(freqwork+1))
        {
        workBS -= col;
        workbeta -= col;
        if(deriv == true)
          workBd -= col;
        }
      i++;
      freqwork++;
      workindex2++;
      workspline += *workindex2;
      if(deriv == true)
        workderiv += *workindex2;
      lp += *workindex2*lpcols;
      }
    k++;
    }

  }


void FULLCOND_projection::reset_effect(const unsigned & pos)
  {
  subtr_spline();
  unsigned i;
  double * work;
  work = spline.getV();
  for(i=0;i<spline.rows();i++,work++)
    *work = 0.0;
  work = beta.getV();
  for(i=0;i<nrpar;i++,work++)
    *work = 0.0;

  intercept = 0.0;
  }


void FULLCOND_projection::update_stepwise(double la)
  {
  if(smoothing == "global")
    {
    lambda=la;

    /*if(likep->iwlsweights_constant() == true)
      {
      bool gefunden = false;
      unsigned i = 0;
      while(i<lambdavec.size() && gefunden == false)
        {
        if(lambda == lambdavec[i])
          gefunden = true;
        i++;
        }
      if(gefunden == true)
        {
        prec_env = all_precenv[i-1];
        lambda_prec = lambda;
        }
      }*/
    }
/*  else
    {
    lambda_prec = -1;
    lambdaold = -1;
    lambda = 1;
    lambdas_local(lambda_nr+nrpar-rankK,0) = 1/la;
    updateK(lambdas_local);
    }  */
  }


double FULLCOND_projection::compute_df(void)
  {
  double df = 0;
  //if(inthemodel == false && fixornot == true)
  //  df = 1;
  if(inthemodel == true)
    {
    if (lambdaold == lambda && likep->get_iwlsweights_notchanged() == true)
      {
      df = df_lambdaold;
      }
    else if (lambdaold != lambda || likep->get_iwlsweights_notchanged() == false)
      {
      if (likep->get_iwlsweights_notchanged() == false)
        compute_XWXenv(likep->get_weightiwls(),column);
      if(lambda != lambda_prec || likep->iwlsweights_constant() == false)
        {
        prec_env.addto(XX_env,Kenv,1.0,lambda);
        lambda_prec = lambda;
        }
      invprec = envmatdouble(0,nrpar,prec_env.getBandwidth());
      prec_env.inverse_envelope(invprec);
      df = df + invprec.traceOfProduct(XX_env);
      if(!identifiable)
        df -= 1;

if(smoothing == "local")
  df = df/rankK;

      df_lambdaold = df;
      lambdaold = lambda;
      }
    }

  return df;
  }


void FULLCOND_projection::hierarchie_rw1(vector<double> & untervector, int dfo)
  {

  unsigned number = untervector.size()-1;

  update_stepwise(untervector[0]);
  double df_max = compute_df();

  update_stepwise(untervector[number]);
  double df_min = compute_df();

  if(df_max > dfo && df_min < dfo)
     {
     bool geordnet = false;
     unsigned stelle_oben = number;
     unsigned stelle_unten = 0;
     while(geordnet==false)
        {
        unsigned stelle = stelle_oben + stelle_unten;
        update_stepwise(untervector[stelle/2]);
        double df_mitteunten = compute_df();
        update_stepwise(untervector[stelle/2 + 1]);
        double df_mitteoben = compute_df();

        if(df_mitteunten > dfo && df_mitteoben > dfo)
          stelle_unten = stelle/2;
        else if(df_mitteunten < dfo && df_mitteoben < dfo)
          stelle_oben = stelle/2 + 1;
        else
          {
          geordnet = true;
          vector<double> hilf;
          unsigned i;
          stelle_unten = stelle/2;
          stelle_oben = stelle/2 + 1;
          for(i=0;i<=stelle_unten;i++)
             hilf.push_back(untervector[i]);
          hilf.push_back(-1);
          for(i=stelle_oben;i<untervector.size();i++)
            hilf.push_back(untervector[i]);
          untervector = hilf;
          }
        }
     }
  else if(df_min >= dfo)
     {
     untervector.push_back(-1);
     }
  else
     {
     vector<double> hilf;
     hilf.push_back(-1);
     unsigned i;
     for(i=0;i<untervector.size();i++)
        hilf.push_back(untervector[i]);
     untervector = hilf;
     }
  }


void FULLCOND_projection::compute_lambdavec(
vector<double> & lvec, int & number)
  {
  if(number>0)
    {
    if (df_equidist==true && spfromdf==true)
       FULLCOND::compute_lambdavec_equi(lvec,number);
    else
       FULLCOND::compute_lambdavec(lvec,number);
    }

  /*lambdavec = lvec;
  compute_XWXenv(likep->get_weightiwls(),column);
  for(unsigned i=0;i<lambdavec.size();i++)
    {
    prec_env.addto(XX_env,Kenv,1.0,lambdavec[i]);
    prec_env.decomp();
    all_precenv.push_back(prec_env);
    }*/

/*    ist der fixe Effekt hier sinnvoll???
  if(type==RW1 && number>0)
    hierarchie_rw1(lvec,1);
  else  // if(type==RW2 || (type==RW1 && number==-1))
    lvec.push_back(-1);
*/

  if(forced_into==false)
     lvec.push_back(0);

  // Startwert für lambda aus df:
  if(spfromdf==true)
    {
    double lambdavorg = 1000;
    if(!varcoeff)
      {
      if(dfstart==1)
        lambdastart = -1;
      else if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    else
      {
      if(dfstart==1 && identifiable)
        lambdastart = -1;
      else if((dfstart==2 && identifiable) || (dfstart==1 && !identifiable))
        lambdastart = -2;
      else if(dfstart==0)
        lambdastart = 0;
      else
        lambdastart = lambda_from_df(dfstart,lambdavorg);
      }
    if(lambdastart==-9 || lambdastart==1000000000)    // falls dfstart nicht erreicht werden kann
      lambdastart = 0;
    }

  }


const datamatrix & FULLCOND_projection::get_data_forfixedeffects(void)
  {
  return data_forfixed;
  }


ST::string FULLCOND_projection::get_effect(void)
  {
  ST::string h;

  h = datanames[0] + "(projection_rw";

  ST::string w = "weights=(";
  double * workg = pp_weights.getV();
  for(unsigned i=0;i<pp_weights.rows()-1;i++,workg++)
    w = w + ST::doubletostring(*workg,4) + ",";
  w = w + ST::doubletostring(*workg,4) + ")";

  if (type== MCMC::RW1)
    h = h + "1," + w + "df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";
  else
    h = h + "2," + w + ",df=" + ST::doubletostring(compute_df(),6) + ",(lambda=" + ST::doubletostring(lambda,6) + "))";

  return h;
  }


void FULLCOND_projection::init_names(const vector<ST::string> & na)
  {

//nrvar -= 1;

  ST::string underscore = "\\_";
  ST::string helpname = na[0];
  unsigned i;
  for(i=1;i<nrvar;i++)
    helpname = helpname + "*" + na[i];

  if(pp_pointer.size()>0)
    helpname = "(" + helpname + ")_" + ST::inttostring(pp_pointer.size()+1);

  datanames = vector<ST::string>(1,helpname);

  helpname = helpname.insert_string_char('_',underscore);
  term_symbolic = "f_{" + helpname + "}(" + helpname + ")";
  priorassumptions.push_back("$" + term_symbolic + "$:");

  if(column > 0)
    //term_symbolic = term_symbolic + " (" + ST::inttostring(column+1) + ". response category)";
    priorassumptions.push_back("$" + term_symbolic
           + " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)$:");

  if(type==MCMC::RW1)
     priorassumptions.push_back("Linear combination using P-spline with first order random walk penalty");
  else if(type==MCMC::RW2)
     priorassumptions.push_back("Linear combination using P-spline with second order random walk penalty");

  ST::string knotstr;
  if (knpos == equidistant)
    knotstr = "equidistant";
  else if (knpos == quantiles)
    knotstr = "quantiles";
  priorassumptions.push_back("Number of knots: " + ST::inttostring(nrknots));
  priorassumptions.push_back("Knot choice: " + knotstr);
  priorassumptions.push_back("Degree of Splines: " + ST::inttostring(degree));

//nrvar += 1;
  }


/*ST::string FULLCOND_projection::get_befehl(void)
  {
  ST::string h;

  if(varcoeff)
    h = datanames[1] + "(psplinerw";
  else
    h = datanames[0] + "(psplinerw";
  if (type== MCMC::RW1)
    h = h + "1";
  else
    h = h + "2";

  h = h + ", lambda=" + ST::doubletostring(lambda,6);
  if(degree!=3)
    h = h + ", degree=" + ST::inttostring(degree);
  if(nrknots!=20)
    h = h + ",nrknots=" + ST::inttostring(nrknots);
  h = h + ")";

  return h;
  }*/


void FULLCOND_projection::outresults(void)
  {
  //spline_basis::outresults();

  optionsp->out("  Results are stored in file\n");
  optionsp->out("  " + pathcurrent + "\n");
  optionsp->out("\n");
  #if defined(BORLAND_OUTPUT_WINDOW)
  optionsp->out("  Results may be visualized using the S-Plus function 'plotnonp'\n");
  ST::string doublebackslash = "\\\\";
  ST::string spluspath = pathcurrent.insert_string_char('\\',doublebackslash);
  optionsp->out("  Type for example:\n");
  optionsp->out("  plotnonp(\"" + spluspath + "\")");
  optionsp->out("\n");
  #elif defined(JAVA_OUTPUT_WINDOW)
  optionsp->out("  Postscript file is stored in file\n");
  ST::string psfile = pathcurrent.substr(0,pathcurrent.length()-4) + ".ps";
  optionsp->out("  " + psfile + "\n");
  optionsp->out("\n");
  optionsp->out("  Results may be visualized using method 'plotnonp'\n");
  optionsp->out("  Type for example: objectname.plotnonp " + ST::inttostring(fcnumber) + "\n");
  #endif
  optionsp->out("\n");

  unsigned i;
  ofstream outres(pathcurrent.strtochar());

  outres << "intnr" << "   ";
  outres << datanames[0] << "   ";
  outres << "pmean   ";

  outres << endl;

  double * workmean = fchelp.get_betameanp();
  double * workxvalues = xvalues.getV();

  for(i=0;i<xvalues.rows();i++,workmean++,workxvalues++)
    {
    outres << (i+1) << "   ";
    outres << *workxvalues << "   ";
    outres << *workmean << "   ";
    outres << endl;
    }

// Ausgabe als mehrdimensionale Funktion

  ST::string pathgesamt = pathcurrent.substr(0,pathcurrent.length()-4)+"_complete.res";
  ofstream outges(pathgesamt.strtochar());
  int j;
  outges << "intnr" << "   ";

  char stern = '*';
  ST::string leer = "   ";
  ST::string names = datanames[0];
  names = names.insert_string_char(stern,leer);
  if( names[0] == '(' )
    names = names.substr(1,names.length()-1);
  if( names[names.length()-3] == ')' )
    names = names.substr(0,names.length()-3);
  outges << names << "   ";

  outges << "pmean" << "   " << endl;

  double * workg = gesamt.getV();
  workmean = spline.getV();

  j=pp_pointer.size()-1;
  bool gefunden = false;
  bool drin,fix;
  double *workges;
  while(gefunden == false && j>=0)
    {
    pp_pointer[j]->get_inthemodel(drin,fix);
    if(drin==true)
      {
      gefunden = true;
      workges = pp_pointer[j]->get_gesamt().getV();
      }
    j -= 1;
    }

  for(i=0;i<spline.rows();i++,workg++,workmean++,workges++)
    {
    if(gefunden)
      *workg = *workges + *workmean;
    else
      *workg = *workmean;
    }

  double * dataori = original_data.getV();
  workg = gesamt.getV();
  double hilf;
  transform = likep->get_trmult(column);
  for(i=0;i<spline.rows();i++,workg++,workmean++)
    {
    outges << (i+1) << "   ";
    for(j=0;j<int(nrvar);j++,dataori++)
      outges << *dataori << "   ";
    hilf = *workg * transform;
    outges << ST::doubletostring(hilf,6) << "   " << endl;
    }
  }

} // end: namespace MCMC


#ifdef __BUILDING_GNU
int main()
{
	return(0);
}
#endif















