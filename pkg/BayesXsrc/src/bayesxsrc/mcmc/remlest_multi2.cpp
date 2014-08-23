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



#include "remlest_multi2.h"

#if defined(BORLAND_OUTPUT_WINDOW)
#include "StatResults.h"
#include "statwinframe.h"

#endif

using std::ofstream;
using std::flush;

//------------------------------------------------------------------------------
//-------------------------- CLASS: remlest_ordinal ----------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------------- Constructor ------------------------------------
//------------------------------------------------------------------------------

remlest_ordinal::remlest_ordinal(
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
  unsigned i, j, k;

  xcut.push_back(0);
  zcut.push_back(0);

  for(i=0;i<fullcond.size();i++)
    {
    xcut.push_back(xcut[i]+fullcond[i]->get_dimX());
    if (i>0)
      {
      zcut.push_back(zcut[i-1]+fullcond[i]->get_dimZ());
      }
    }

  X = datamatrix(re.rows(),xcut[xcut.size()-1],0);
  Z = datamatrix(re.rows(),zcut[zcut.size()-1],0);

  fullcond[0]->createreml(X,Z,xcut[0],0);

  for(i=1;i<fullcond.size();i++)
    {
    fullcond[i]->createreml(X,Z,xcut[i],zcut[i-1]);
    }

  catspec=false;
  catspecific_fixed = (dynamic_cast <MCMC::FULLCOND_const*> (fullcond[0]))->get_catspecific_fixed();
  for(i=0; i<fullcond.size(); i++)
    {
    catspecific.push_back(fullcond[i]->get_catspecific());
    }
  for(i=1; i<fullcond.size(); i++)
    {
    if(catspecific[i])
      {
      catspec=true;
      }
    }

  for(i=1; i<catspecific_fixed.size(); i++)
    {
    if(catspecific_fixed[i])
      {
      catspec=true;
      }
    }

  if(catspec)
    {
    unsigned nrcatspecific_fixed=0;
    for(i=0; i<catspecific_fixed.size(); i++)
      {
      if(catspecific_fixed[i])
        {
        nrcatspecific_fixed++;
        }
      }

    xcutbeta.push_back(0);
    zcutbeta.push_back(0);
    xcutbeta.push_back(fullcond[0]->get_dimX()+nrcatspecific_fixed*(nrcat2-1));
    for(i=1; i<fullcond.size(); i++)
      {
      if(catspecific[i])
        {
        for(j=0; j<nrcat2; j++)
          {
          xcutbeta.push_back(xcutbeta[xcutbeta.size()-1]+fullcond[i]->get_dimX());
          zcutbeta.push_back(zcutbeta[zcutbeta.size()-1]+fullcond[i]->get_dimZ());
          }
        }
      else
        {
        xcutbeta.push_back(xcutbeta[xcutbeta.size()-1]+fullcond[i]->get_dimX());
        zcutbeta.push_back(zcutbeta[zcutbeta.size()-1]+fullcond[i]->get_dimZ());
        }
      }
    totalnrfixed = xcutbeta[xcutbeta.size()-1];
    totalnrpar = totalnrfixed+zcutbeta[zcutbeta.size()-1];

    beta=statmatrix<double>(totalnrpar,1,0);
    theta=statmatrix<double>(zcutbeta.size()-1,1,0);

    k=0;
    for(i=1; i<fullcond.size(); i++)
      {
      if(catspecific[i])
        {
        for(j=0; j<nrcat2; j++)
          {
          theta(k,0) = fullcond[i]->get_startlambda();
          k++;
          }
        }
      else
        {
        theta(k,0) = fullcond[i]->get_startlambda();
        k++;
        }
      }
    }
  else
    {
    totalnrfixed=X.cols()+nrcat2-1;
    totalnrpar=totalnrfixed+Z.cols();

    beta=statmatrix<double>(totalnrpar,1,0);
    theta=statmatrix<double>(zcut.size()-1,1,0);

    for(i=1; i<fullcond.size(); i++)
      {
      theta(i-1,0) = fullcond[i]->get_startlambda();
      }
    }

// Startwerte für die Schwellenwerte bestimmen
// Berechne Häufigkeiten
  if(respfamily=="cumlogit" || respfamily=="cumprobit")
    {
    datamatrix freq(nrcat2,1,0);
    for(i=0; i<nrobs; i++)
      {
      if(weight(i,0)!=0)
        {
        for(j=0; j<nrcat2; j++)
          {
          if(re(i,0)==cats(j,0))
            {
            freq(j,0) += 1;
            }
          }
        }
      }
// Kumulieren
    for(j=1; j<nrcat2; j++)
      {
      freq(j,0) += freq(j-1,0);
      }
// Kumulierte relative Häufigkeiten bestimmen und Schwellenwerte berechnen
    for(j=0; j<nrcat2; j++)
      {
      freq(j,0) /= nrobspos;
      if(family=="cumlogit")
        {
        beta(j,0) = log(freq(j,0)/(1-freq(j,0)));
        }
     else
        {
        beta(j,0) = randnumbers::invPhi2(freq(j,0));
        }
      }
    }
  }

//------------------------------------------------------------------------------
//----------------------------- Estimation -------------------------------------
//------------------------------------------------------------------------------

bool remlest_ordinal::estimate2(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight)
  {
  unsigned i, j, k, l;

  outoptions();
  out("\n");

  for(i=0;i<fullcond.size();i++)
    fullcond[i]->outoptionsreml();

  //----------------------------------------------------------------------------
  //--- Konstruiere großes X und großes Z
  //----------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------------------
// Zu X: Ausgabe der Testdateien zu X

// Ausgabe von Xalt unter conny1
/*  ofstream out1("c:\\bayesx\\mcmc\\conny1");
        X.prettyPrint(out1);
        out1.close();*/

// Hilfsmatrix zu xcut unter help2
/*  datamatrix help2(xcut.size(),1,0);
  for(i=0; i<xcut.size(); i++)
  {
  help2(i,0) = xcut[i];
  }

// Ausgabe von xcut unter conny2
   ofstream out2("c:\\bayesx\\mcmc\\conny2");
        help2.prettyPrint(out2);
        out2.close();*/


// Hilfsmatrix zu xcutbeta unter help4
/*  datamatrix help4(xcutbeta.size(),1,0);
  for(i=0; i<xcutbeta.size(); i++)
  {
  help4(i,0) = xcutbeta[i];
  }

// Ausgabe von xcutbeta unter conny4
  ofstream out4("c:\\bayesx\\mcmc\\conny4");
        help4.prettyPrint(out4);
        out4.close();*/

//--------------------------------------------------------------------------------------------------------
// Zu Z: Ausgabe der Testdateien zu Z

// Ausgabe von Zalt unter conny5
/*  ofstream out5("c:\\bayesx\\mcmc\\conny5");
        Z.prettyPrint(out5);
        out5.close();*/

// Hilfsmatrix zu zcut unter help6
/*  datamatrix help6(zcut.size(),1,0);
  for(i=0; i<zcut.size(); i++)
  {
  help6(i,0) = zcut[i];
  }

// Ausgabe von zcut unter conny6
  ofstream out6("c:\\bayesx\\mcmc\\conny6");
        help6.prettyPrint(out6);
        out6.close();*/


// Hilfsmatrix zu zcutbeta unter help8
/*  datamatrix help8(zcutbeta.size(),1,0);
  for(i=0; i<zcutbeta.size(); i++)
  {
  help8(i,0) = zcutbeta[i];
  }

// Ausgabe von zcutbeta unter conny8
  ofstream out8("c:\\bayesx\\mcmc\\conny8");
        help8.prettyPrint(out8);
        out8.close();*/

//-----------------------------------------------------------------------------------------------
// Zu catspecific und catspecific_fixed: Ausgabe der Testdateien

// Hilfsmatrix zu catspecific
/*  datamatrix help0(catspecific.size(),1,0);
  for(i=0; i<catspecific.size(); i++)
  {
  help0(i,0) = catspecific[i];
  }

// Ausgabe von catspecific unter catspecific
  ofstream outcatspecific("c:\\bayesx\\mcmc\\catspecific");
        help0.prettyPrint(outcatspecific);
        outcatspecific.close(); */


// Hilfsmatrix zu catspecific_fixed
/*  datamatrix help01(catspecific_fixed.size(),1,0);
  for(i=0; i<catspecific_fixed.size(); i++)
  {
  help01(i,0) = catspecific_fixed[i];
  }

// Ausgabe von catspecific unter conny01
  ofstream out01("c:\\bayesx\\mcmc\\conny01");
        help01.prettyPrint(out01);
        out01.close();*/




//------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------
// Basteln an X

//---------------------------------------------
// Zu X: Hilfsvektor bzw. Hilfsmatrix für die if - Schleifen

// xcutbetalength gibt die Länge der Abschnitte von xcutbeta an
  vector <unsigned> xcutbetalength (xcutbeta.size()-1);
  for (i=0; i<xcutbetalength.size(); i++)
  {
  xcutbetalength[i]=xcutbeta[i+1]-xcutbeta[i];
  }

// Hilfsmatrix zu xcutbetalength namens Xlength
/* datamatrix Xlength(xcutbetalength.size(), 1, 0);
  for (i=0; i<xcutbetalength.size(); i++)
  {
  Xlength(i,0)=xcutbeta[i+1]-xcutbeta[i];
  }

// Ausgabe von Xlength
  ofstream outXlength("c:\\bayesx\\mcmc\\Xlength");
        Xlength.prettyPrint(outXlength);
        outXlength.close();   */

//----------------------------------------------------------------------------------------------------------------------------------
// xcutlength gibt die Länge der Abschnitte von xcut an
  vector <unsigned> xcutlength (xcut.size()-1);
  for (i=0; i<xcutlength.size(); i++)
  {
  xcutlength[i]=xcut[i+1]-xcut[i];
  }

// Hilfsmatrix zu xcutlength namens Xlength
/* datamatrix Xcutlength(xcutlength.size(), 1, 0);
  for (i=0; i<xcutlength.size(); i++)
  {
  Xcutlength(i,0)=xcut[i+1]-xcut[i];
  }

// Ausgabe von Xlength
  ofstream outXcutlength("c:\\bayesx\\mcmc\\Xcutlength");
        Xcutlength.prettyPrint(outXcutlength);
        outXcutlength.close(); */


//----------------------------------------------
// Erzeugt Xneu mit 0
   datamatrix Xneu (nrobs*nrcat2,xcutbeta[xcutbeta.size()-1],0);

//----------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------
// Belegt die ersten Spalten von Xneu
// Durchlaufen von catspecific_fixed

// Zwei neue Variablen, beginnend bei 0
   unsigned XaltSpalte_fixed=0;                                                  // gibt die aktuelle Spalte in Xalt an
   unsigned XneuSpalte_fixed=0;                                                  // gibt die aktuelle Spalte in Xneu an

// j durchläuft den Bereich 0..catspecific_fixed.size()           // catspecific_fixed.size() = xcut[1]

   for (j=0; j<catspecific_fixed.size(); j++)                                    // Durchlauf solange Bedingung vorhanden
   {                                                                             // Anfang zur for-Schleife für X

      //----------------------------------
      // 1. Skalar und nicht catspecific
      // -> nrcat2 x 1 - Matrizen
      // Buchstaben E,F

      if ( !catspecific_fixed[j])
      {
         for(i=0; i<nrobs; i++)                                                  // i durchläuft die Zeilen von Xalt
         {
             for(k=0; k<nrcat2; k++)                                             // k durchläuft nrcat2
             {
             Xneu(i*nrcat2+k, XneuSpalte_fixed)= X (i, XaltSpalte_fixed );
             }
         }

      XaltSpalte_fixed=XaltSpalte_fixed+1;                                       // aktualisiert XaltSpalte_fixed
      XneuSpalte_fixed=XneuSpalte_fixed+1;                                       // aktualisiert XneuSpalte_fixed
      }

      //----------------------------------
      // 2. Skalar und catspecific
      // -> nrcat2 x nrcat2 - Matrizen
      // Buchstaben G, H

      else
      {
         for(i=0; i<nrobs; i++)                                                  // i durchläuft die Zeilen von Xalt
         {
             for(k=0; k<nrcat2; k++)                                             // k durchläuft nrcat2
             {
             Xneu(i*nrcat2+k, XneuSpalte_fixed+k)= X (i, XaltSpalte_fixed );
             }
         }

      XaltSpalte_fixed=XaltSpalte_fixed+1;                                       // aktualisiert XaltSpalte_fixed
      XneuSpalte_fixed=XneuSpalte_fixed+nrcat2;                                  // aktualisiert XneuSpalte_fixed
      }

    }

//---------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------
// Schleife, um Xneu zu erweitern
// Durchlaufen von catspecific

// Zwei neue Variablen, beginnend bei vektor[1], da vorher catspecific_fixed
   unsigned XaltSpalte=xcut[1];                                                  // gibt die aktuelle Spalte in Xalt an
   unsigned XneuSpalte=xcutbeta[1];                                              // gibt die aktuelle Spalte in Xneu an

// j durchläuft den Bereich 1..catspecific.size()

   for (j=1; j<catspecific.size(); j++)                                          // Durchlauf solange Bedingung vorhanden
   {                                                                             // Anfang zur for-Schleife für X

      //----------------------------------
      // 1. Vektor und nicht catspecific
      // -> nrcat2 x 1*xcutbetalength[j] - Matrizen
      // Buchstaben I, K

      if ( !catspecific[j])
      {
            for(i=0; i<nrobs; i++)                                               // i durchläuft die Zeilen von Xalt
            {
               for(k=0; k<nrcat2; k++)                                           // k durchläuft nrcat2
               {
                  for(l=0; l<xcutlength[j]; l++)                                 // l durchläuft xcutbetalength[j]
                  {
                  Xneu(i*nrcat2+k, XneuSpalte+l)= X (i, XaltSpalte+l );
                  }
               }
            }

      XaltSpalte=XaltSpalte+xcutlength[j];                                       // aktualisiert XaltSpalte
      XneuSpalte=XneuSpalte+xcutlength[j];                                       // aktualisiert XneuSpalte
      }

      //----------------------------------
      // 2. Vektor und catspecific
      // -> nrcat2 x nrcat2*xcutbetalength[j] - Matrizen
      // Buchstaben L, M

      else
      {
            for(i=0; i<nrobs; i++)                                               // i durchläuft die Zeilen von Xalt
            {
               for(k=0; k<nrcat2; k++)                                           // k durchläuft nrcat2
               {
                  for(l=0; l<xcutlength[j]; l++)                                 // l durchläuft xcutbetalength[j]
                  {
                  Xneu(i*nrcat2+k, XneuSpalte+(k*xcutlength[j]) +l)= X (i, XaltSpalte+l );
                  }
               }
            }

      XaltSpalte=XaltSpalte+xcutlength[j];                                       // aktualisiert XaltSpalte
      XneuSpalte=XneuSpalte+(xcutlength[j]*nrcat2);                              // aktualisiert XneuSpalte
      }

   }                                                                             // Ende zur for-Schleife für X

   // Ausgabe von Xneu unter Xneu
   /*  ofstream outXneu("c:\\bayesx\\mcmc\\Xneu");
        Xneu.prettyPrint(outXneu);
        outXneu.close();   */


//-------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------
//  Basteln an Z

//---------------------------------------------
// Zu Z: Hilfsvektor bzw. Hilfsmatrix für die if - Schleifen

// zcutbetalength gibt die Länge der Abschnitte von zcutbeta an
  vector <unsigned> zcutbetalength (zcutbeta.size()-1);
  for (i=0; i<zcutbetalength.size(); i++)
  {
  zcutbetalength[i]=zcutbeta[i+1]-zcutbeta[i];
  }

// Hilfsmatrix zu zcutbetalength namens Zlength
/* datamatrix Zlength(zcutbetalength.size(), 1, 0);
  for (i=0; i<zcutbetalength.size(); i++)
  {
  Zlength(i,0)=zcutbeta[i+1]-zcutbeta[i];
  }

// Ausgabe von Zlength
  ofstream outZlength("c:\\bayesx\\mcmc\\Zlength");
        Zlength.prettyPrint(outZlength);
        outZlength.close();  */

//------------------------------------------------------------------
// zcutlength gibt die Länge der Abschnitte von zcut an
  vector <unsigned> zcutlength (zcut.size()-1);
  for (i=0; i<zcutlength.size(); i++)
  {
  zcutlength[i]=zcut[i+1]-zcut[i];
  }

// Hilfsmatrix zu zcutlength namens Zcutlength
/* datamatrix Zcutlength(zcutlength.size(), 1, 0);
  for (i=0; i<zcutlength.size(); i++)
  {
  Zcutlength(i,0)=zcut[i+1]-zcut[i];
  }

// Ausgabe von Zcutlength
  ofstream outZcutlength("c:\\bayesx\\mcmc\\Zcutlength");
        Zcutlength.prettyPrint(outZcutlength);
        outZcutlength.close();     */

//------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------
// Erzeugt Zneu mit Nullen
datamatrix Zneu (nrobs*nrcat2,zcutbeta[zcutbeta.size()-1],0);

// Zwei neue Variablen, beginnend bei 0
   unsigned ZaltSpalte=0;                                                        // gibt die aktuelle Spalte in Zalt an
   unsigned ZneuSpalte=0;                                                        // gibt die aktuelle Spalte in Zneu an

// j durchläuft den Bereich 1..catspecific.size()

   for (j=1; j<catspecific.size(); j++)                                          // Durchlauf solange Bedingung vorhanden
   {                                                                             // Anfang zur for-Schleife für Z

      //----------------------------------
      // 1. Vektor und nicht catspecific
      // -> nrcat2 x 1*xcutbetalength[j-1] - Matrizen
      // Buchstaben I, K

     if (!catspecific[j])
       {
            for(i=0; i<nrobs; i++)                                               // i durchläuft die Zeilen von Zalt
            {
               for(k=0; k<nrcat2; k++)                                           // k durchläuft nrcat2
               {
                  for(l=0; l<zcutlength[j-1]; l++)                               // l durchläuft zcutbetalength[j-1]
                  {
                  Zneu(i*nrcat2+k, ZneuSpalte+l)= Z (i, ZaltSpalte+l );
                  }
               }
        }

      ZaltSpalte=ZaltSpalte+zcutlength[j-1];                                     // aktualisiert ZaltSpalte
      ZneuSpalte=ZneuSpalte+zcutlength[j-1];                                     // aktualisiert ZneuSpalte
      }

      //----------------------------------
      // 2. Vektor und catspecific
      // -> nrcat2 x nrcat2*xcutbetalength[j-1] - Matrizen
      // Buchstaben L, M

      else
      {
            for(i=0; i<nrobs; i++)                                               // i durchläuft die Zeilen von Zalt
            {
               for(k=0; k<nrcat2; k++)                                           // k durchläuft nrcat2
               {
                  for(l=0; l<zcutlength[j-1]; l++)                               // l durchläuft zcutbetalength[j-1]
                  {
                  Zneu(i*nrcat2+k, ZneuSpalte+(k*zcutlength[j-1]) +l)= Z (i, ZaltSpalte+l );
                  }
               }
            }

      ZaltSpalte=ZaltSpalte+zcutlength[j-1];                                     // aktualisiert ZaltSpalte
      ZneuSpalte=ZneuSpalte+(zcutlength[j-1]*nrcat2);                            // aktualisiert ZneuSpalte
      }

   }                                                                             // Ende zur for-Schleife für Z

// Ausgabe von Zneu unter Zneu
/*  ofstream outZneu("c:\\bayesx\\mcmc\\Zneu");
        Zneu.prettyPrint(outZneu);
        outZneu.close();  */

//-----------------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------------

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
  statmatrix<double>Qinv(Zneu.cols(),1,0);

  // Matrices for Fisher scoring (regression parameters)
  statmatrix<double>H(beta.rows(),beta.rows(),0);
  statmatrix<double>H1(beta.rows(),1,0);

  // Matrices for Fisher scoring (variance parameters)
  statmatrix<double>Hinv(beta.rows(),beta.rows(),0);
  statmatrix<double>wresid(1,nrcat2*resp.rows(),0);                  //matrix containing row vectors !!

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

//------------------------------------------------------------------------------
// neue Funktion compute_eta2
//    compute_eta2(eta);
    eta = Xneu*beta.getRowBlock(0,totalnrfixed)+Zneu*beta.getRowBlock(totalnrfixed,totalnrpar);
//------------------------------------------------------------------------------

    compute_weights(mu,workweight,worky,eta,respind,weight);

    stop = check_pause();
    if (stop)
      return true;

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// neue Funktion compute_sscp2 (Conny)

// Ausgabe von workweight
/*  ofstream outWorkweight("c:\\bayesx\\mcmc\\Workweight");
        workweight.prettyPrint(outWorkweight);
        outWorkweight.close();*/


//---------------------------------------------------------------------------------------------------------------------------------------------
// zu Xneu und Zneu: Dimensionen
// von oben:
// datamatrix   Xneu            (   nrobs*nrcat2,    xcutbeta[xcutbeta.size()-1],  0);
// datamatrix   Zneu            (   nrobs*nrcat2,    zcutbeta[zcutbeta.size()-1],  0);
//
// datamatrix   workweight   (   nrobs*nrcat2,     nrcat2,      0);
// datamatrix   worky          (    nrobs*nrcat2,     1,        0);

        // Ausgabe von Worky unter Worky
/*           ofstream outWorky("c:\\bayesx\\mcmc\\Worky");
           worky.prettyPrint(outWorky);
           outWorky.close();*/



//--------------------------------------------------------------------------------------------------------------------------------------------
// Berechnen von Xneu' unter XneuTr

   datamatrix XneuTr (xcutbeta[xcutbeta.size()-1], nrobs*nrcat2, 0);
   XneuTr = Xneu.transposed();

        // Ausgabe von XneuTr unter XneuTr
       /*   ofstream outXneuTr ("c:\\bayesx\\mcmc\\XneuTr");
           XneuTr.prettyPrint(outXneuTr);
           outXneuTr.close();  */

//--------------------------------------------------------------------------------------------------------------------------------------------
// Berechnen von Zneu' unter ZneuTr

   datamatrix ZneuTr (zcutbeta[zcutbeta.size()-1], nrobs*nrcat2, 0);
   ZneuTr = Zneu.transposed();

        // Ausgabe von ZneuTr unter ZneuTr
      /*    ofstream outZneuTr ("c:\\bayesx\\mcmc\\ZneuTr");
           ZneuTr.prettyPrint(outZneuTr);
           outZneuTr.close(); */

//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------
// Berechnen von H
//
//               Xneu' workweight Xneu           Xneu' workweight Zneu              H_00       H_01
//    H  =                                                                   =
//               Zneu' workweight Xneu           Zneu' workweight Zneu              H_10       H_11

//---------------------------------------------------------------------------------------------------------------------------------------
// Berechnen von H1
//
//                 Xneu' workweight worky                 H1_0
//   H1 =                                          =
//                 Zneu' workweight worky                 H1_1


//----------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------
// Berechnen von H_00, H_01, H1_0

//----------------------------------------------------------------------
// Erzeugen von H_00 = Xneu' workweight Xneu
   datamatrix  H_00 (xcutbeta[xcutbeta.size()-1],  xcutbeta[xcutbeta.size()-1],  0);

//----------------------------------------------------------------------
// Erzeugen von H_01 = Xneu' workweight Zneu
   datamatrix  H_01 (xcutbeta[xcutbeta.size()-1],  zcutbeta[zcutbeta.size()-1],  0);

//----------------------------------------------------------------------
// Erzeugen von H1_0 = Xneu' workweight worky
   datamatrix  H1_0 (xcutbeta[xcutbeta.size()-1],  1, 0);


//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Hilfsmatrix zum Berechnen der Zeilen von Xneu' workweight
   datamatrix H_00hilf (1, nrobs*nrcat2,  0);

for (l=0; l<xcutbeta[xcutbeta.size()-1]; l++ )                                   // l durchläuft die Zeilen von Xneu' und die Zeilen von H_00
{
     for  (j=0;  j<nrobs*nrcat2;  j=j+nrcat2)                                                                                                                      // j durchläuft die Spalten von Xneu'
     {
          for (i=0; i<nrcat2; i++)                                                                                                                                       // i durchläuft nrcat2
          {
          H_00hilf (0, j+i) = ((XneuTr.getBlock( l, j, l+1, j+nrcat2)) * (workweight.getBlock(j, i,  j+nrcat2, i+1)))(0,0);
          }
     }

                  // Ausgabe von H_00hilf  unter H_00hilf
/*                     ofstream outH_00hilf  ("c:\\bayesx\\mcmc\\H_00hilf ");
                     H_00hilf .prettyPrint(outH_00hilf );
                     outH_00hilf .close();*/


//----------------------------------------------------------------------------
// Berechnen von H_00 = Xneu' workweight Xneu

     for (k=l; k<xcutbeta[xcutbeta.size()-1]; k++ )                              // k durchläuft die Spalten von Xneu
     {
     H_00(l, k) = ((H_00hilf) * (Xneu.getCol(k)))(0,0);
     H_00(k,l)=H_00(l,k);
     }

                  // Ausgabe von H_00  unter H_00
/*                     ofstream outH_00 ("c:\\bayesx\\mcmc\\H_00 ");
                     H_00.prettyPrint(outH_00 );
                     outH_00.close();*/
//-----------------------------------------------------------------------
// Berechnen von H_01 = Xneu' workweight Zneu
     for (k=0; k<zcutbeta[zcutbeta.size()-1]; k++ )                              // k durchläuft die Spalten von Zneu
     {
     H_01(l, k) = ((H_00hilf) * (Zneu.getCol(k))) (0,0);
     }

//----------------------------------------------------------------------
// Berechnen von H_10 = Zneu' workweight Xneu
// durch Einsetzen von H_01.transposed() in H

//----------------------------------------------------------------------
// Berechnen von H1_0 = Xneu' workweight worky
    H1_0(l, 0) = ((H_00hilf) * (worky)) (0,0);
//----------------------------------------------------------------------

}                                                                                // Ende der l-Schleife

//----------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------
// Berechnen von H_11, H1_1

//----------------------------------------------------------------------------------------------------------------------------------
// Erzeugen von H_11 = Zneu' workweight Zneu
   datamatrix  H_11 (zcutbeta[zcutbeta.size()-1],  zcutbeta[zcutbeta.size()-1],  0);

//----------------------------------------------------------------------------------------------------------------------------------
// Erzeugen von H1_1 = Zneu' workweight worky
   datamatrix  H1_1 (zcutbeta[zcutbeta.size()-1],  1, 0);


//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Hilfsmatrix zum Berechnen der Zeilen von Zneu' workweight
   datamatrix H_11hilf (1, nrobs*nrcat2,  0);

for (l=0; l<zcutbeta[zcutbeta.size()-1]; l++ )                                   // l durchläuft die Zeilen von Zneu' und die Zeilen von H_11
{

     for  (j=0;  j<nrobs*nrcat2; j=j+nrcat2)                                                                                                                      // j durchläuft die Spalten von Zneu'
     {
          for (i=0; i<nrcat2; i++)                                                                                                                                       // i durchläuft nrcat2
          {
          H_11hilf (0, j+i) = ((ZneuTr.getBlock( l, j, l+1, j+nrcat2)) * (workweight.getBlock(j, i,  j+nrcat2, i+1))) (0,0);
          }
     }

                  // Ausgabe von H_11hilf  unter H_11hilf
/*                     ofstream outH_11hilf  ("c:\\bayesx\\mcmc\\H_11hilf ");
                     H_11hilf .prettyPrint(outH_11hilf );
                     outH_11hilf .close();*/

//-----------------------------------------------------------------------------------------------------------------------------------
// Berechnen von H_11 = Zneu' workweight Zneu
     for (k=l; k<zcutbeta[zcutbeta.size()-1]; k++ )                              // k durchläuft die Spalten von Zneu
     {
     H_11(l, k) = ((H_11hilf) * (Zneu.getCol(k))) (0,0);
     H_11(k,l)=H_11(l,k);
     }

                  // Ausgabe von H_11  unter H_11
/*                     ofstream outH_11 ("c:\\bayesx\\mcmc\\H_11 ");
                     H_11.prettyPrint(outH_11 );
                     outH_11 .close();*/

//-----------------------------------------------------------------------------------------------------------------------------------
// Berechnen von H1_1 = Zneu' workweight worky
    H1_1(l, 0) = ((H_11hilf) * (worky)) (0,0);
}                                                                                // Ende der l-Schleife

//-----------------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------------
// Zusammensetzen von H

// datamatrix H (xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1], xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1],  0);

   H.putBlock (H_00 , 0, 0,
                                       xcutbeta[xcutbeta.size()-1],  xcutbeta[xcutbeta.size()-1] );

   H.putBlock (H_01 , 0, xcutbeta[xcutbeta.size()-1],
                                       xcutbeta[xcutbeta.size()-1],  xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1] );

   H.putBlock (H_01.transposed() , xcutbeta[xcutbeta.size()-1],  0,
                                       xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1] , xcutbeta[xcutbeta.size()-1] );

   H.putBlock (H_11 , xcutbeta[xcutbeta.size()-1],  xcutbeta[xcutbeta.size()-1],
                                       xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1] , xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1] );

        // Ausgabe von H unter H
/*          ofstream outH ("c:\\bayesx\\mcmc\\H");
           H.prettyPrint(outH);
           outH.close();  */

//-----------------------------------------------------------------------------------------------------------------------------------
// Zusammensetzen von H1

//   datamatrix H1 (xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1], 1, 0);

   H1.putRowBlock(0, xcutbeta[xcutbeta.size()-1], H1_0 );
   H1.putRowBlock(xcutbeta[xcutbeta.size()-1], xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1], H1_1 );


        // Ausgabe von H1 unter H1
   /*        ofstream outH1 ("c:\\bayesx\\mcmc\\H1");
           H1.prettyPrint(outH1);
           outH1.close();  */

//---------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------


//    compute_sscp2(H,H1,workweight,worky);



//------------------------------------------------------------------------------

    H.addtodiag(Qinv,totalnrfixed,totalnrpar);

    // Fisher-Scoring for beta
    beta=H.solve(H1);

    // update linear predictor and compute residuals

//------------------------------------------------------------------------------
// neue Funktion compute_eta2
//    compute_eta2(eta);
    eta = Xneu*beta.getRowBlock(0,totalnrfixed)+Zneu*beta.getRowBlock(totalnrfixed,totalnrpar);
//------------------------------------------------------------------------------

    worky = worky - eta;
    for(i=0; i<resp.rows(); i++)
      {
      for(j=0; j<nrcat2; j++)
        {
        wresid(0,i*nrcat2+j)=(workweight.getRow(i*nrcat2+j)*worky.getRowBlock(i*nrcat2,(i+1)*nrcat2))(0,0);
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

//------------------------------------------------------------------------------
// anpassen (zcutbeta)

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
/*     for(i=0; i<theta.rows(); i++)
      {
      for(l=zcutbeta[i]; l<zcutbeta[i+1]; l++)
        {
        help=0;
        for(j=0; j<nrcat2; j++)
          {
          help += (wresid.getRow(j)*(Z.getCol(l)))(0,0);
          }
        score(i,0) += 0.5*help*help*thetaold(i,0)*2;
        }
      }*/
     for(i=0; i<theta.rows(); i++)
      {
      for(l=zcutbeta[i]; l<zcutbeta[i+1]; l++)
        {
        help = ( wresid*( Zneu.getCol(l) ) )(0,0);
        score(i,0) += help*help*thetaold(i,0);
        }
      }
//------------------------------------------------------------------------------

/*ofstream out20("c:\\temp\\score.raw");
score.prettyPrint(out20);
out20.close();
ofstream out21("c:\\temp\\Fisher.raw");
Fisher.prettyPrint(out21);
out21.close();*/


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

    //compute norm of eta for the different categories
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
    help=helpmat.max(0);

//------------------------------------------------------------------------------
// anpassen

    // compute norm of the random parts
/*   for(i=0; i<theta.rows(); i++)
     {
     helpmat2=Z.getColBlock(zcut[i],zcut[i+1])*beta.getRowBlock(totalnrfixed+zcut[i],totalnrfixed+zcut[i+1]);
     stopcrit[i]=helpmat2.norm(0)/help;
     if(stopcrit[i]<lowerlim)
       {
       theta(i,0)=thetaold(i,0);
       }
     else
       {
       its[i]=it;
       }
     } */
   k=0;
   for(i=1; i<fullcond.size(); i++)
     {
     if(catspecific[i])
       {
       for(j=0; j<nrcat2; j++)
         {
         helpmat2=Z.getColBlock(zcut[i-1],zcut[i])*beta.getRowBlock(totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[k+1]);
         stopcrit[k]=helpmat2.norm(0)/help;
         k++;
         }
       }
     else
       {
       helpmat2=Z.getColBlock(zcut[i-1],zcut[i])*beta.getRowBlock(totalnrfixed+zcutbeta[k],totalnrfixed+zcutbeta[k+1]);
       stopcrit[k]=helpmat2.norm(0)/help;
       k++;
       }
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

//theta(0,0)=thetaold(0,0);
//theta(1,0)=thetaold(1,0);

//------------------------------------------------------------------------------

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

  for(i=nrcat2; i<beta.rows(); i++)
    {
    beta(i,0) = -beta(i,0);
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

  k=1;
  for(i=1; i<fullcond.size(); i++)
    {
    if(catspecific[i])
      {
      for(j=0; j<nrcat2; j++)
        {
        beta(j,0) -= fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],k-1,false,xcutbeta[k],totalnrfixed+zcutbeta[k-1],cats(j,0),true,k);
        k++;
        }
      }
    else
      {
//      help = fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],k-1,false,xcutbeta[k]+nrcat2-1,totalnrfixed+zcutbeta[k-1],0,false,k);
      help = fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],k-1,false,xcutbeta[k],totalnrfixed+zcutbeta[k-1],0,false,k);
      k++;
      for(j=0; j<nrcat2; j++)
        {
        beta(j,0) -= help;
        }
      }
    }

// -----------------------------------------------------------------------------
// anpassen:
  ( dynamic_cast <MCMC::FULLCOND_const*> (fullcond[0]) )->outresultsreml_ordinal(X,Z,beta,Hinv,nrcat2);
// -----------------------------------------------------------------------------

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

bool remlest_ordinal::estimate(const datamatrix resp, const datamatrix & offset,
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
  statmatrix<double>Qinv(Z.cols(),1,0);

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
      for(k=zcut[i]; k<zcut[i+1]; k++, l++)
        {
        Qinv(l,0)=1/theta(i,0);
        }
      }

    compute_eta2(eta);

    compute_weights(mu,workweight,worky,eta,respind,weight);

    stop = check_pause();
    if (stop)
      return true;

    compute_sscp2(H,H1,workweight,worky);
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
      score(i,0)=-1*(H.getBlock(totalnrfixed+zcut[i],totalnrfixed+zcut[i],totalnrfixed+zcut[i+1],totalnrfixed+zcut[i+1])*thetaold(i,0)).trace()+
                 ((H.getRowBlock(totalnrfixed+zcut[i],totalnrfixed+zcut[i+1]))*Hinv*(H.getColBlock(totalnrfixed+zcut[i],totalnrfixed+zcut[i+1]))*thetaold(i,0)).trace();
      for(k=0; k<theta.rows(); k++)
        {
        Fisher(i,k)=2*(H.getBlock(totalnrfixed+zcut[i],totalnrfixed+zcut[k],totalnrfixed+zcut[i+1],totalnrfixed+zcut[k+1])*H.getBlock(totalnrfixed+zcut[k],totalnrfixed+zcut[i],totalnrfixed+zcut[k+1],totalnrfixed+zcut[i+1])*thetaold(i,0)*thetaold(k,0)).trace()-
                    4*(H.getRowBlock(totalnrfixed+zcut[k],totalnrfixed+zcut[k+1])*Hinv*H.getColBlock(totalnrfixed+zcut[i],totalnrfixed+zcut[i+1])*H.getBlock(totalnrfixed+zcut[i],totalnrfixed+zcut[k],totalnrfixed+zcut[i+1],totalnrfixed+zcut[k+1])*thetaold(i,0)*thetaold(k,0)).trace()+
                    2*(H.getRowBlock(totalnrfixed+zcut[i],totalnrfixed+zcut[i+1])*Hinv*H.getColBlock(totalnrfixed+zcut[k],totalnrfixed+zcut[k+1])*H.getRowBlock(totalnrfixed+zcut[k],totalnrfixed+zcut[k+1])*Hinv*H.getColBlock(totalnrfixed+zcut[i],totalnrfixed+zcut[i+1])*thetaold(i,0)*thetaold(k,0)).trace();
        Fisher(k,i)=Fisher(i,k);
        }
      }

    for(i=0; i<theta.rows(); i++)
      {
      for(l=zcut[i]; l<zcut[i+1]; l++)
        {
        help=0;
        for(j=0; j<nrcat2; j++)
          {
          help += (wresid.getRow(j)*(Z.getCol(l)))(0,0);
          }
        score(i,0) += help*help*thetaold(i,0);
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
    help=helpmat.max(0);
    // compute norm of the random parts
   for(i=0; i<theta.rows(); i++)
     {
     helpmat2=Z.getColBlock(zcut[i],zcut[i+1])*beta.getRowBlock(totalnrfixed+zcut[i],totalnrfixed+zcut[i+1]);
     stopcrit[i]=helpmat2.norm(0)/help;
     if(stopcrit[i]<lowerlim || theta(i,0)>maxvar)
       {
       theta(i,0)=thetaold(i,0);
       }
     else
       {
       its[i]=it;
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
    if(stopcrit[i]<lowerlim)
      {
      thetareml(i,1)=1;
      }
    else if(theta(i,0)>maxvar)
      {
      thetareml(i,1)=-1;
      }
    thetareml(i,2)=its[i];
    thetareml(i,3)=xcut[i+2]-xcut[i+1]+(Hhelp.getBlock(totalnrfixed+zcut[i],totalnrfixed+zcut[i],totalnrfixed+zcut[i+1],totalnrfixed+zcut[i+1])).trace();
    }

  for(i=nrcat2; i<beta.rows(); i++)
    {
    beta(i,0) = -beta(i,0);
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

  for(i=1; i<fullcond.size(); i++)
    {
    help = fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],i-1,false,xcut[i]+nrcat2-1,totalnrfixed+zcut[i-1],0,false,i);
    for(j=0; j<nrcat2; j++)
      {
      beta(j,0) -= help;
      }
    }
  ( dynamic_cast <MCMC::FULLCOND_const*> (fullcond[0]) )->outresultsreml_ordinal(X,Z,beta,Hinv,nrcat2);

// log-likelihood, AIC, BIC, etc.
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

bool remlest_ordinal::estimate_glm(const datamatrix resp,
                  const datamatrix & offset, const datamatrix & weight)
  {
  unsigned i, j;

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

    compute_sscp(H,H1,workweight,worky);

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
  for(i=nrcat2; i<beta.rows(); i++)
    {
    beta(i,0) = -beta(i,0);
    }

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

  ( dynamic_cast <MCMC::FULLCOND_const*> (fullcond[0]) )->outresultsreml_ordinal(X,Z,beta,H,nrcat2);

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

bool remlest_ordinal::estimate_glm2(const datamatrix resp,
                  const datamatrix & offset, const datamatrix & weight)
  {
  unsigned i, j, k, l;                                                           // k und l dazugefügt

  outoptions();
  out("\n");

//-------------------------------------------------------------------------------------------------------------------------------------
// Conny: Xneu ausrechnen
  datamatrix Xneu(nrobs*nrcat2,beta.rows(),0);

//---------------------------------------------------------
// Belegt die ersten Spalten von Xneu
// Durchlaufen von catspecific_fixed

// Zwei neue Variablen, beginnend bei 0
   unsigned XaltSpalte_fixed=0;                                                  // gibt die aktuelle Spalte in Xalt an
   unsigned XneuSpalte_fixed=0;                                                  // gibt die aktuelle Spalte in Xneu an

// j durchläuft den Bereich 0..catspecific_fixed.size()           // catspecific_fixed.size() = xcut[1]

   for (j=0; j<catspecific_fixed.size(); j++)                                    // Durchlauf solange Bedingung vorhanden
   {                                                                             // Anfang zur for-Schleife für X

      //----------------------------------
      // 1. Skalar und nicht catspecific
      // -> nrcat2 x 1 - Matrizen
      // Buchstaben E,F
      // verbesserte Version!!!

      if ( !catspecific_fixed[j])
      {
         for(i=0; i<nrobs; i++)                                                  // i durchläuft die Zeilen von Xalt
         {
             for(k=0; k<nrcat2; k++)                                             // m durchläuft nrcat2
             {
             Xneu(i*nrcat2+k, XneuSpalte_fixed)= X (i, XaltSpalte_fixed );
             }
         }

      XaltSpalte_fixed=XaltSpalte_fixed+1;                                       // aktualisiert XaltSpalte_fixed
      XneuSpalte_fixed=XneuSpalte_fixed+1;                                       // aktualisiert XneuSpalte_fixed
      }

      //----------------------------------
      // 2. Skalar und catspecific
      // -> nrcat2 x nrcat2 - Matrizen
      // Buchstaben G, H

      else
      {
         for(i=0; i<nrobs; i++)                                                  // i durchläuft die Zeilen von Xalt
         {
             for(k=0; k<nrcat2; k++)                                             // k durchläuft nrcat2
             {
             Xneu(i*nrcat2+k, XneuSpalte_fixed+k)= X (i, XaltSpalte_fixed );
             }
         }

      XaltSpalte_fixed=XaltSpalte_fixed+1;                                       // aktualisiert XaltSpalte_fixed
      XneuSpalte_fixed=XneuSpalte_fixed+nrcat2;                                  // aktualisiert XneuSpalte_fixed
      }

    }
//-------------------------------------------------------------------------------------------------------------------------------------


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

//--------------------------------------------------------------------------------------------------------------------------------------------
// Berechnen von Xneu' unter XneuTr

   datamatrix XneuTr (xcutbeta[xcutbeta.size()-1], nrobs*nrcat2, 0);
   XneuTr = Xneu.transposed();

        // Ausgabe von XneuTr unter XneuTr
/*           ofstream outXneuTr ("c:\\bayesx\\mcmc\\XneuTr");
           XneuTr.prettyPrint(outXneuTr);
           outXneuTr.close();*/
//-------------------------------------------------------------------------------


  // Estimation loop
  while(test==true)
    {
    // store current values in betaold
    betaold=beta;

//    compute_eta(eta);
    eta = Xneu*beta;

    compute_weights(mu,workweight,worky,eta,respind,weight);

//-----------------------------------------------------------------------------------------------------------------------------------------------
// Conny: H und H1 nur mit Xneu ausrechnen
//    compute_sscp(H,H1,workweight,worky);


//----------------------------------------------------------------------------------------------------------------------------------
// Berechnen von H, H1

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// Erzeugen von H = Xneu' workweight Xneu
//   datamatrix  H (xcutbeta[xcutbeta.size()-1],  xcutbeta[xcutbeta.size()-1],  0);

//----------------------------------------------------------------------
// Erzeugen von H1 = Xneu' workweight worky
//   datamatrix  H1 (xcutbeta[xcutbeta.size()-1],  1, 0);

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Hilfsmatrix zum Berechnen der Zeilen von Xneu' workweight
   datamatrix H_00hilf (1, nrobs*nrcat2,  0);

for (l=0; l<xcutbeta[xcutbeta.size()-1]; l++ )                                   // l durchläuft die Zeilen von Xneu' und die Zeilen von H_00
{
     for  (j=0;  j<nrobs*nrcat2;  j=j+nrcat2)                                                                                                                      // j durchläuft die Spalten von Xneu'
     {
          for (i=0; i<nrcat2; i++)                                                                                                                                       // i durchläuft nrcat2
          {
          H_00hilf (0, j+i) = ((XneuTr.getBlock( l, j, l+1, j+nrcat2)) * (workweight.getBlock(j, i,  j+nrcat2, i+1)))(0,0);
          }
     }

                  // Ausgabe von H_00hilf  unter H_00hilf
/*                     ofstream outH_00hilf  ("c:\\bayesx\\mcmc\\H_00hilf ");
                     H_00hilf .prettyPrint(outH_00hilf );
                     outH_00hilf .close();*/


//----------------------------------------------------------------------------
// Berechnen von H_00 = Xneu' workweight Xneu

     for (k=l; k<xcutbeta[xcutbeta.size()-1]; k++ )                              // k durchläuft die Spalten von Xneu
     {
     H(l, k) = ((H_00hilf) * (Xneu.getCol(k)))(0,0);
     H(k,l)=H(l,k);
     }

                  // Ausgabe von H_00  unter H_00
/*                     ofstream outH_00 ("c:\\bayesx\\mcmc\\H_00 ");
                     H_00.prettyPrint(outH_00 );
                     outH_00.close();*/

//----------------------------------------------------------------------
// Berechnen von H1 = Xneu' workweight worky
    H1(l, 0) = ((H_00hilf) * (worky)) (0,0);
//-----------------------------------------------------------

}                                                                                                // Ende der l-Schleife


//-----------------------------------------------------------------------------------------------------------------------------------------------







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
  for(i=nrcat2; i<beta.rows(); i++)
    {
    beta(i,0) = -beta(i,0);
    }

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

  ( dynamic_cast <MCMC::FULLCOND_const*> (fullcond[0]) )->outresultsreml_ordinal(X,Z,beta,H,nrcat2);

  loglike=0;
  aic=0;
  bic=0;
  gcv=0;
  df=beta.rows();
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

void remlest_ordinal::compute_respind(const datamatrix & re, datamatrix & respind)
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

void remlest_ordinal::compute_eta(datamatrix & eta)
  {
  unsigned i,j;
  double help;
  for(i=0; i<nrobs; i++)
    {
    help=(X.getBlock(i,1,i+1,X.cols())*beta.getRowBlock(nrcat2,totalnrfixed))(0,0);
    for(j=0; j<nrcat2; j++)
      {
      eta(i*nrcat2+j,0)=beta(j,0)+help;
      }
    }
  }

void remlest_ordinal::compute_eta2(datamatrix & eta)
  {
  unsigned i,j;
  double help;
  if(X.cols()>1)
    {
    for(i=0; i<nrobs; i++)
      {
      help=(X.getBlock(i,1,i+1,X.cols())*beta.getRowBlock(nrcat2,totalnrfixed))(0,0)
           +(Z.getRow(i)*beta.getRowBlock(totalnrfixed,beta.rows()))(0,0);
      for(j=0; j<nrcat2; j++)
        {
        eta(i*nrcat2+j,0)=beta(j,0)+help;
        }
      }
    }
  else
    {
    for(i=0; i<nrobs; i++)
      {
      help=(Z.getRow(i)*beta.getRowBlock(totalnrfixed,beta.rows()))(0,0);
      for(j=0; j<nrcat2; j++)
        {
        eta(i*nrcat2+j,0)=beta(j,0)+help;
        }
      }
    }
  }

void remlest_ordinal::compute_weights(datamatrix & mu, datamatrix & workweights,
                  datamatrix & worky, datamatrix & eta, datamatrix & respind,
                  const datamatrix & weights)
  {
  unsigned i,k,j,l;
  double expo;
  datamatrix expos = datamatrix(nrcat2,1,0);

// Compute mu
  if(respfamily=="cumlogit")
    {
    for(i=0; i<nrobs; i++)
      {
      expo=exp(eta(i*nrcat2,0));
      mu(i*nrcat2,0)=expo/(1+expo);
      for(j=1; j<nrcat2; j++)
        {
        expo=exp(eta(i*nrcat2+j,0));
        mu(i*nrcat2+j,0)=expo/(1+expo);
        for(k=0; k<j; k++)
          {
          mu(i*nrcat2+j,0) -= mu(i*nrcat2+k,0);
          }
        }
      }
    }
  else if(respfamily=="cumprobit")
    {
    for(i=0; i<nrobs; i++)
      {
      mu(i*nrcat2,0)=randnumbers::Phi2(eta(i*nrcat2,0));
      for(j=1; j<nrcat2; j++)
        {
        mu(i*nrcat2+j,0)=randnumbers::Phi2(eta(i*nrcat2+j,0));
        for(k=0; k<j; k++)
          {
          mu(i*nrcat2+j,0) -= mu(i*nrcat2+k,0);
          }
        }
      }
    }
  else if(respfamily=="seqlogit")
    {
    for(i=0; i<nrobs; i++)
      {
      for(j=0; j<nrcat2; j++)
        {
        expos(j,0)=exp(eta(i*nrcat2+j,0));
        }
      mu(i*nrcat2,0)=expos(0,0)/(1+expos(0,0));
      for(j=1; j<nrcat2; j++)
        {
        mu(i*nrcat2+j,0)=expos(j,0)/(1+expos(j,0));
        for(k=0; k<j; k++)
          {
          mu(i*nrcat2+j,0) *= 1/(1+expos(k,0));
          }
        }
      }
    }
  else if(respfamily=="seqprobit")
    {
    for(i=0; i<nrobs; i++)
      {
      mu(i*nrcat2,0)=randnumbers::Phi2(eta(i*nrcat2,0));
      for(j=1; j<nrcat2; j++)
        {
        mu(i*nrcat2+j,0)=randnumbers::Phi2(eta(i*nrcat2+j,0));
        for(k=0; k<j; k++)
          {
          mu(i*nrcat2+j,0) *= 1-randnumbers::Phi2(eta(i*nrcat2+k,0));
          }
        }
      }
    }

// Compute weights
  datamatrix D(nrcat2,nrcat2,0);
  datamatrix S(nrcat2,nrcat2,0);
  for(i=0; i<nrobs; i++)
    {
    if(weights(i,0)>0)
      {
      for(j=0; j<nrcat2; j++)
        {
        S(j,j) = mu(i*nrcat2+j,0)*(1-mu(i*nrcat2+j,0));
        for(k=j+1; k<nrcat2; k++)
          {
          S(j,k) = -mu(i*nrcat2+j,0)*mu(i*nrcat2+k,0);
          S(k,j) = S(j,k);
          }
        }
      if(respfamily=="cumlogit")
        {
        for(j=0; j<nrcat2; j++)
          {
          expo = exp(eta(i*nrcat2+j,0));
          D(j,j) = expo/((1+expo)*(1+expo));
          }
        for(j=0; j<nrcat2-1; j++)
          {
          D(j,j+1) = -D(j,j);
          }
        }
      else if(respfamily=="cumprobit")
        {
        for(j=0; j<nrcat2; j++)
          {
          D(j,j) = randnumbers::phi(eta(i*nrcat2+j,0));
          }
        for(j=0; j<nrcat2-1; j++)
          {
          D(j,j+1) = -D(j,j);
          }
        }
      else if(respfamily=="seqlogit")
        {
        for(j=0; j<nrcat2; j++)
          {
          expos(j,0)=exp(eta(i*nrcat2+j,0));
          }
        for(j=0; j<nrcat2; j++)
          {
          D(j,j)=expos(j,0)/((1+expos(j,0))*(1+expos(j,0)));
          for(k=0; k<j; k++)
            {
            D(j,j) *= 1/(1+expos(k,0));
            }
          }
        for(j=0; j<nrcat2; j++)
          {
          for(k=j+1; k<nrcat2; k++)
            {
            D(j,k) = -expos(k,0)/(1+expos(k,0)) * expos(j,0)/((1+expos(j,0))*(1+expos(j,0)));
            for(l=0; l<j && l<k; l++)
              {
              D(j,k) *= 1/(1+expos(l,0));
              }
            for(l=j+1; l<k; l++)
              {
              D(j,k) *= 1/(1+expos(l,0));
              }
            }
          }
        }
      else if(respfamily=="seqprobit")
        {
        for(j=0; j<nrcat2; j++)
          {
          D(j,j) = randnumbers::phi(eta(i*nrcat2+j,0));
          for(k=0; k<j; k++)
            {
            D(j,j) *= 1-randnumbers::Phi2(eta(i*nrcat2+k,0));
            }
          }
        for(j=0; j<nrcat2; j++)
          {
          for(k=j+1; k<nrcat2; k++)
            {
            D(j,k) = -randnumbers::Phi2(eta(i*nrcat2+k,0))*randnumbers::phi(eta(i*nrcat2+j,0));
            for(l=0; l<j && l<k; l++)
              {
              D(j,k) *= 1-randnumbers::Phi2(eta(i*nrcat2+l,0));
              }
            for(l=j+1; l<k; l++)
              {
              D(j,k) *= 1-randnumbers::Phi2(eta(i*nrcat2+l,0));
              }
            }
          }
        }

      workweights.putBlock(D*S.inverse()*(D.transposed()),i*nrcat2,0,(i+1)*nrcat2,nrcat2);

  // Compute worky;
      worky.putRowBlock(i*nrcat2,(i+1)*nrcat2,eta.getRowBlock(i*nrcat2,(i+1)*nrcat2)+
                        (D.inverse()).transposed()*(respind.getRowBlock(i*nrcat2,(i+1)*nrcat2)-mu.getRowBlock(i*nrcat2,(i+1)*nrcat2)));
      }
    }

  }

void remlest_ordinal::compute_sscp(datamatrix & H, datamatrix & H1,
                      datamatrix & workweight, datamatrix & worky)
  {
  unsigned i;
  H=datamatrix(H.rows(),H.cols(),0);
  H1=datamatrix(H1.rows(),1,0);
  datamatrix weighti=datamatrix(nrcat2,nrcat2,0);
  datamatrix yi=datamatrix(nrcat2,1,0);
  datamatrix xi=datamatrix(X.cols()-1,1,0);
  datamatrix Htemp=datamatrix(H.rows(),H.cols(),0);
  datamatrix H1temp=datamatrix(H1.rows(),1,0);
  datamatrix colsums=datamatrix(nrcat2,1,0);
  double totalsum;
  for(i=0; i<nrobs; i++)
    {
    weighti=workweight.getBlock(i*nrcat2,0,(i+1)*nrcat2,nrcat2);
    colsums=weighti.sum();
    totalsum=colsums.sum(0);

    xi = X.getBlock(i,1,i+1,X.cols()).transposed();

    Htemp.putBlock(weighti,0,0,nrcat2,nrcat2);
    Htemp.putBlock(colsums*xi.transposed(),0,nrcat2,nrcat2,H.cols());
    Htemp.putBlock((Htemp.getBlock(0,nrcat2,nrcat2,H.cols())).transposed(),nrcat2,0,H.cols(),nrcat2);
    Htemp.putBlock(xi*xi.transposed()*totalsum,nrcat2,nrcat2,H.rows(),H.cols());

    H.plus(Htemp);

    yi = worky.getRowBlock(i*nrcat2,(i+1)*nrcat2);
    H1temp.putRowBlock(0,nrcat2,weighti*yi);
    H1temp.putRowBlock(nrcat2,H1temp.rows(),xi*((colsums.transposed()*yi)(0,0)));

    H1.plus(H1temp);
    }
  }

void remlest_ordinal::compute_sscp2(datamatrix & H, datamatrix & H1,
                      datamatrix & workweight, datamatrix & worky)
  {
  unsigned i;
  H=datamatrix(H.rows(),H.cols(),0);
  H1=datamatrix(H1.rows(),1,0);
  datamatrix weighti=datamatrix(nrcat2,nrcat2,0);
  datamatrix yi=datamatrix(nrcat2,1,0);
  datamatrix xi=datamatrix(totalnrpar-nrcat2,1,0);
  datamatrix Htemp=datamatrix(H.rows(),H.cols(),0);
  datamatrix H1temp=datamatrix(H1.rows(),1,0);
  datamatrix colsums=datamatrix(nrcat2,1,0);
  double totalsum;
  for(i=0; i<nrobs; i++)
    {
    weighti=workweight.getBlock(i*nrcat2,0,(i+1)*nrcat2,nrcat2);
    colsums=weighti.sum();
    totalsum=colsums.sum(0);

    if(X.cols()>1)
      {
      xi.putRowBlock(0,X.cols()-1,X.getBlock(i,1,i+1,X.cols()).transposed());
      xi.putRowBlock(X.cols()-1,totalnrpar-nrcat2,Z.getRow(i).transposed());
      }
    else
      {
      xi=Z.getRow(i).transposed();
      }

    Htemp.putBlock(weighti,0,0,nrcat2,nrcat2);
    Htemp.putBlock(colsums*xi.transposed(),0,nrcat2,nrcat2,H.cols());
    Htemp.putBlock((Htemp.getBlock(0,nrcat2,nrcat2,H.cols())).transposed(),nrcat2,0,H.cols(),nrcat2);
    Htemp.putBlock(xi*xi.transposed()*totalsum,nrcat2,nrcat2,H.rows(),H.cols());

    H.plus(Htemp);

    yi = worky.getRowBlock(i*nrcat2,(i+1)*nrcat2);
    H1temp.putRowBlock(0,nrcat2,weighti*yi);
    H1temp.putRowBlock(nrcat2,H1temp.rows(),xi*((colsums.transposed()*yi)(0,0)));

    H1.plus(H1temp);
    }
  }

//------------------------------------------------------------------------------
//----------------------------- Object description -----------------------------
//------------------------------------------------------------------------------

  void remlest_ordinal::outoptions()
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
    ST::string familyname;
    if(respfamily=="cumlogit")
      {
      familyname="cumulative logit (ordered categories)";
      }
    else if(respfamily=="cumprobit")
      {
      familyname="cumulative probit (ordered categories)";
      }
    else if(respfamily=="seqlogit")
      {
      familyname="sequential logit";
      }
    else if(respfamily=="seqprobit")
      {
      familyname="sequential probit";
      }
    out("  Family:                 "+familyname+"\n");
    out("  Number of observations: "+ST::inttostring(X.rows())+"\n");
    out("  Number of observations with positive weight: "+ST::inttostring(nrobspos)+"\n");
    }

//------------------------------------------------------------------------------
//----------------------------- Writing results --------------------------------
//------------------------------------------------------------------------------

void remlest_ordinal::make_plots(ofstream & outtex,ST::string path_batch,
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

    unsigned k;

    // durchlaufen der Fullconditionals
    for(j=0;j<fullcond.size();j++)
      {

      // Pfad der Regr.-Ergebnisse
      pathresult = fullcond[j]->get_pathresult();

      // Plotstyle: noplot, plotnonp, drawmap, drawmapgraph
      plst = fullcond[j]->get_plotstyle();

      if(!catspecific[j])
        {
        k=1;
        }
      else
        {
        k=nrcat2;
        }

      for(i=0; i<k; i++)
        {
        if (plst != MCMC::noplot)
          {
          pathresult = fullcond[j]->get_pathresult();
          if(catspecific[j])
            {
            pathresult = pathresult.insert_after_string(ST::doubletostring(cats(i,0),6)+"_","_f_");
            }

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
            if(catspecific[j])
              {
              outtex << " (Category " << cats(i,0) << ")";
              }
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
              if(catspecific[j])
                {
                outtex << " (Category " << cats(i,0) << ")";
                }
              outtex << ". Shown are the posterior modes.}" << endl
                     << "\\end{figure}" << endl;
              outtex << "\n\\begin{figure}[htb]" << endl
                     << "\\centering" << endl
                     << "\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                     << u_str << "}" << endl
                     << "\\caption{Non--linear Effect of '" << effect << "'";
              if(catspecific[j])
                {
                outtex << " (Category " << cats(i,0) << ")";
                }
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
              if(catspecific[j])
                {
                outtex << " (Category " << cats(i,0) << ")";
                }
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
              if(catspecific[j])
                {
                outtex << " (Category " << cats(i,0) << ")";
                }
              outtex << ". Shown are the posterior modes.}" << endl
                     << "%\\end{figure}" << endl;
              outtex << "\n%\\begin{figure}[htb]" << endl
                     << "%\\centering" << endl
                     << "%\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                     << u_str << "}" << endl
                     << "%\\caption{Non--linear Effect of '" << effect << "'";
              if(catspecific[j])
                {
                outtex << " (Category " << cats(i,0) << ")";
                }
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
              if(catspecific[j])
                {
                outtex << " (Category " << cats(i,0) << ")";
                }
              outtex << ". Posterior probabilities for a nominal level of "
                     << o_str << "\\%." << endl
                     << "%Black denotes regions with strictly negative credible intervals,"
                     << endl
                     << "%white denotes regions with strictly positive credible intervals.}"
                     << endl << "%\\end{figure}" << endl;
              }

            }
          }
        }
      }
    }
  }

void remlest_ordinal::make_model(ofstream & outtex, const ST::string & rname)
  {
    ST::string familyname;
    if(respfamily=="cumlogit")
      {
      familyname="cumulative logit (ordered categories)";
      }
    else if(respfamily=="cumprobit")
      {
      familyname="cumulative probit (ordered categories)";
      }
    else if(respfamily=="seqlogit")
      {
      familyname="sequential logit";
      }
    else if(respfamily=="seqprobit")
      {
      familyname="sequential probit";
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
         << "\n\\noindent {\\bf \\large Predictor:}" << endl << endl;
  }

void remlest_ordinal::make_predictor(ofstream & outtex)
  {

  unsigned j;

  ST::string term2 = fullcond[0]->get_term_symbolic();
  ST::string term = "$\\eta_j = " + term2;    //linearer Prädiktor wird erweitert
  for(j=1;j<fullcond.size();j++)
    {
    term2 = fullcond[j]->get_term_symbolic();
    if(catspecific[j])
      {
      term2 = term2.insert_after_all_string("^{(j)}","f");;
      }
    term = term + " - " + term2;    //linearer Prädiktor wird erweitert
    }
  outtex << term << "$\\\\\n";
  }

void remlest_ordinal::make_prior(ofstream & outtex)
  {
  unsigned i,j;
  outtex << "\n\\noindent {\\bf \\large Priors:}\\\\" << endl << "\\\\" << endl;
  for(j=0;j<fullcond.size();j++)
    {
    vector<ST::string> prior = fullcond[j]->get_priorassumptions();
    if(prior.size() != 0)// nur wenn Priors da sind (d.h. Vektor hat Elemente)
      {
      if(fullcond[j]->get_results_type()!="fixed" && catspecific[j])
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

void remlest_ordinal::make_options(ofstream & outtex)
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

void remlest_ordinal::make_fixed_table(ofstream & outtex)
  {

  // falls andere Quantile gewünscht werden
  double u = fullcond[0]->get_level1();
  ST::string u_str = ST::doubletostring(u,0);

  vector<ST::string> h;

  unsigned j;
  unsigned r;

  r = 2;

  // Tabelle im Tex-File mit fixen Effekten
  outtex << "\n\\newpage \n" << endl << "\n\\noindent {\\bf \\large Fixed Effects:}\\\\"
       << endl << "\\\\" << endl;

  outtex << "\\begin{tabular}{|r|rrrrr|}" << endl << "\\hline" << endl
         << "Variable & Post. Mode & Std. Dev. & p-value & \\multicolumn{2}{r|}{" << u << "\\% confidence interval}\\\\"
         << endl << "\\hline" << endl;

  h = fullcond[0]->get_results_latex();
  for (j=0;j<h.size();j++)
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
             << "\n\\noindent {\\bf \\large Fixed Effects (continued):}\\\\"
             << endl << "\\\\" << endl;

      outtex << "\\begin{tabular}{|r|rrrrr|}" << endl << "\\hline" << endl
             << "Variable & Post. Mode & Std. Dev. & p-value & \\multicolumn{2}{r|}{" << u << "\\% confidence interval}\\\\"
             << endl << "\\hline" << endl;

      outtex << h[j] << endl;
      }
    }
  outtex << "\\hline \n\\end{tabular}" << endl;
  }

void remlest_ordinal::make_graphics(const ST::string & title,
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

bool remlest_ordinal::check_pause()
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

void remlest_ordinal::out(const ST::string & s,bool thick,bool italic,
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


void remlest_ordinal::outerror(const ST::string & s)
  {
  out(s,true,true,12,255,0,0);
  }

