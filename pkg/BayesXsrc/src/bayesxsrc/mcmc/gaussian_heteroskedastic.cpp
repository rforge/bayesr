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



#include "gaussian_heteroskedastic.h"

namespace MCMC
{


void DISTRIBUTION_gaussianh::standardise(void)
  {
/*
  double s = sqrt(response.var(0,weight));

  trmult = datamatrix(2,1,s);
  trmult(1,0) = 1.0;

  unsigned i;
  double * workresp = response.getV();
  double * worklin = (*linpred_current).getV();
  for (i=0;i<nrobs;i++,workresp++,worklin++)
   {
   *workresp = *workresp/trmult(0,0);
   workresp++;
   *workresp = *workresp/trmult(0,0);
   *worklin = *worklin/trmult(0,0);
   worklin++;
   *worklin=*worklin-2*log(s);
   }
*/

//  datamatrix tr(1,1,trmult(0,0)*trmult(0,0));
//  Scalesave.set_transformmult(tr);

  }


DISTRIBUTION_gaussianh::DISTRIBUTION_gaussianh(const double & a,
                   const datamatrix & b, MCMCoptions * o, const datamatrix & r,
                   const datamatrix & w)
  : DISTRIBUTION(o,r,w)
  {

  nrcat = response.cols(); //neu

  family = "Gaussian with heteroscedastic errors";

  standardise();

  scale(0,0) = 1;
  scaleexisting = false;

  constant_iwlsweights=true;
  // für den Fall des Mittelwertschätzers wäre dies falsch
  // für den Fall des Varianzschätzers ist dies richtig

  }


DISTRIBUTION_gaussianh::DISTRIBUTION_gaussianh(
const DISTRIBUTION_gaussianh & nd)
: DISTRIBUTION(DISTRIBUTION(nd))
    {
        nrcat = nd.nrcat;//neu
        family = nd.family;
        scale = nd.scale;
        scaleexisting = nd.scaleexisting;
        constant_iwlsweights = nd.constant_iwlsweights;
    }


const DISTRIBUTION_gaussianh & DISTRIBUTION_gaussianh::operator=(
const DISTRIBUTION_gaussianh & nd)
  {

  if (this==&nd)
    return *this;
  DISTRIBUTION::operator=(DISTRIBUTION(nd));

  nrcat = nd.nrcat;//neu
  family = nd.family;
  scale = nd.scale;
  scaleexisting = nd.scaleexisting;
  constant_iwlsweights = nd.constant_iwlsweights;

  return *this;
  }



double DISTRIBUTION_gaussianh::loglikelihood(double * response,
                      double * linpred,
                      double * weight,const int & i) const//für eine Beob.
  {
        double * worklin = linpred;
        double workmu = (*worklin); //erster Prediktor, also der für mu
        worklin++;
        double s = exp((*worklin)); //zweiter Prediktor, also der für die
                                    //Varianz
        double help = (*response) - workmu;
        return -0.5 * (*worklin) - 0.5*(help*help)/s;

        // DISTRIBUTION_gaussian abschauen
        // DISTRIBUTION_multinomial abschauen
        // Die Variablen weight ist hier überflüssig und wird nicht genutzt
  }



void DISTRIBUTION_gaussianh::compute_mu(const double * linpred,double * mu)
                                           const
  {
        const double * worklin = linpred;
        double * workmu = mu;

        //Zunächst den Wert für den Erwartungswertschätzer zuweisen
        (*mu) = trmult(0,0)* (*linpred);

        //Zeiger auf die Einträge des Varianzschätzers richten
        worklin++;
        workmu++;

        //Wert für den Varianzschätzer zuweisen, Transformation nicht notwendig,
        //da in Funktion standardise nicht vorgenommen
        (*workmu) = exp((*worklin));

  // DISTRIBUTION_gaussian abschauen
  // DISTRIBUTION_multinomial abschauen
  // zweite Spalte: exp(linpred für varianz)

  }


void DISTRIBUTION_gaussianh::compute_mu_notransform(const double * linpred,
double * mu) const
  {
        const double * worklin = linpred;
        double * workmu = mu;

        //Zunächst den Wert für den Erwartungswertschätzer zuweisen
        (*mu) = (*linpred);

        //Zeiger auf die Einträge des Varianzschätzers richten
        worklin++;
        workmu++;

        //Wert für den Varianzschätzer zuweisen
        (*workmu) = exp((*worklin));


  // DISTRIBUTION_gaussian abschauen
  // DISTRIBUTION_multinomial abschauen
  // zweite Spalte: exp(linpred für varianz)

  }


void DISTRIBUTION_gaussianh::compute_deviance(const double * response,
                             const double * weight,const double * mu,
                             double * deviance,double * deviancesat,
                             const datamatrix & scale,const int & i) const
  {
      const double  * workmu = mu;
      double r = (*response)* trmult(0,0) - (*mu);
      workmu++;
      double s = (*workmu)*pow(trmult(0,0),2);
      *deviance =  ((double)1.0/s)*r*r+log(2*M_PI*s);
      *deviancesat = ((double)1.0/s)*r*r;

  // DISTRIBUTION_gaussian abschauen
  // DISTRIBUTION_multinomial abschauen
  // zweite Spalte: exp(linpred für varianz)
  // die Variable weight als Argument der Funktion ist wahrscheinlich
  // überflüssig, diese enthält die IWLS-Gewichte nach distribution.h

  }


double DISTRIBUTION_gaussianh::compute_weight(double * linpred, double * weight,
                        const int & i, const unsigned & col) const
  {
  double * worklinpred = linpred;
  worklinpred++; //zeigt jetzt auf eta
  double s = exp((*worklinpred));

  if(col == 0) //Berechnung für den Prädiktor des Mittelwertes
    {
    return   1.0/s;//1/exp(eta)
    }
  else  //Berechnung für den Prädiktor der Varianz if(col == 1)
    {
    return 0.5;
    }

  }


double DISTRIBUTION_gaussianh::compute_IWLS(double * response,double * linpred,
                                           double * weight,
                      const int & i,double * weightiwls,double * tildey,
                      bool weightyes, const unsigned & col)
  {

    double * workweightiwls = weightiwls;
    double * worktildey = tildey;
    double * worklinpred = linpred;
    double workmu = (*worklinpred);
    worklinpred++; //zeigt jetzt auf eta
    double s = exp((*worklinpred));
    double help = (*response)-workmu;


    if(col == 0) //Berechnung für den Prädiktor des Mittelwertes
    {
        if(weightyes)
        {
            (*workweightiwls) = (double)1.0/s;//1/exp(eta)
        }

        (*worktildey) =  (*response); //Für den Erwartungswertschätzer stimmen
                                   //working-obs mit obs überein
    }
    if(col == 1) //Berechnung für den Prädiktor der Varianz
    {
        workweightiwls++;
        worktildey++;

        if(weightyes)
        {
            (*workweightiwls) = 0.5;
//            (*workweightiwls) = 1.0;
        }

//        (*worktildey) = (*worklinpred) + ((help*help)/s) - 1;
        (*worktildey) =  ((help*help)/s) - 1;
    }

    return - 0.5 * (*worklinpred) - 0.5 * (help*help)/s;


  // vgl. distribution h datei
  // Berechnet die IWLS-Gewichte und speichert diese in weightiwls. Als iwls
  // Gewichte werden die negativen Erwartungswerte der jeweiligen Gewichts-
  // matrizen verwendet. Dies entspricht Fisher-Scoring. Inwiefern dies im
  // Rahmen der Programmlogik zulässig ist, ist noch unklar.
  // Berechnet die working-observations ~y und speichert diese in tildey
  // Der i-te Summand der Loglikelihood wird zurückgegeben
  // Die Berechnung erfolgt für die i-te Beobachtung
  // Der Übergabewert weight wird nicht beachtet, da Heteroskedastizität mittels
  // sigma_{i}^{2} modelliert wird.
  // Der Übergabewert col kennzeichnet, ob Berechnungen für den Prädiktor
  // des Mittelwertes (col=0) oder für den der Varianz (col=1) durchgeführt
  // werden.
  // weightyes gibt an, ob die IWLS-Gewichte berechnet werden sollen, d.h.
  // weightyes = true, falls dies erfolgen sollen, weightyes = false, sonst

  }

void DISTRIBUTION_gaussianh::compute_IWLS_weight_tildey(double * response,
                              double * linpred,
                              double * weight,const int & i,
                              double * weightiwls,double * tildey,
                              const unsigned & col)
  {

    double * workweightiwls = weightiwls;
    double * worktildey = tildey;
    double * worklinpred = linpred;
    double workmu = (*worklinpred);
    worklinpred++; //zeigt jetzt auf eta
    double s = exp((*worklinpred));
    double help = (*response)-workmu;


    if(col == 0) //Berechnung für den Prädiktor des Mittelwertes
    {
        (*workweightiwls) = (double)1.0/s;//1/exp(eta)

        (*worktildey) =  (*response); //Für den Erwartungswertschätzer stimmen
                                   //working-obs mit obs überein
    }
    else if(col == 1) //Berechnung für den Prädiktor der Varianz
    {
        workweightiwls++;
        worktildey++;

        (*workweightiwls) = (double)0.5;

//        (*worktildey) = (*worklinpred) + ((help*help)/s) - 1;
        (*worktildey) =  ((help*help)/s) - 1;
    }

  // vgl. distribution h datei
  // im Wesentlichen identisch mit compute_IWLS, allerdings wird nicht
  // abgefragt, ob IWLS-Gewichte berechnet werden sollen und es existiert
  // kein Rückgabewert

  }


void DISTRIBUTION_gaussianh::compute_iwls(void)
  {
    register unsigned i,j;
    unsigned dim = response.cols(); //Der Wert zur Laufzeit müsste 2 sein

    double * worklin = (*linpred_current).getV();

    double * workres = response.getV();
    double * worktildey = tildey.getV();
    double * workweightiwls = weightiwls.getV();
    double help=0.0;


    for (i=0;i<nrobs;i++)
    {
        for(j=0;j<dim;j++,worktildey++,workres++,workweightiwls++)
        {
          //compute_weight(worklin, workweightiwls, &i, &j);
          if(j==0)// Entspricht dem Mittelwertschätzer
          {
            help = (*workres)-(*worklin);
            worklin++;
            double s = exp((*worklin));

            (*workweightiwls) = 1.0/s;//1/exp(eta)
            (*worktildey) = (*workres);
          }
          else if(j==1)//Entspricht dem Varianzschätzer
          {
            (*workweightiwls) = 0.5;
            (*worktildey) = (*worklin) + ((help*help)/exp((*worklin))) - 1;
//            (*worktildey) =  ((help*help)/exp((*worklin))) - 1;
            worklin++;
          }

        }
    }


 }


double DISTRIBUTION_gaussianh::compute_gmu(double * linpred,
const unsigned & col) const
  {

  double * worklinpred = linpred;

  if(col == 0)
  {
    return (double)0;
  }
  else if(col == 1)
  {
    worklinpred++;
    return (double)1.0/exp((*worklinpred));
  }

  // vgl. distribution h datei
  // berechne g'(mu_i) = 1/h'(eta_i)

  return 0; //ist dies sinnvoll?
  }


void DISTRIBUTION_gaussianh::outoptions(void)
  {
  DISTRIBUTION::outoptions();

  optionsp->out("\n");

  }


void DISTRIBUTION_gaussianh::update(void)
  {

    DISTRIBUTION::update();

  }


void DISTRIBUTION_gaussianh::update_predict(void)
  {
    DISTRIBUTION::update_predict();
  }


void DISTRIBUTION_gaussianh::outresults(void)
  {



  }

bool DISTRIBUTION_gaussianh::posteriormode(void)
  {

  return true;

  }

/*
bool DISTRIBUTION_gaussianh::posteriormode_converged_fc(const datamatrix & beta,
                                  const datamatrix & beta_mode,
                                  const unsigned & itnr)
  {
  return true;
  }
*/




} // end: namespace MCMC

