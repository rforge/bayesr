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



#include "remlest_multi3.h"

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

remlest_multinomial_catsp::remlest_multinomial_catsp(
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
  unsigned i,j,k;

  catspecific_fixed = (dynamic_cast <MCMC::FULLCOND_const*> (fullcond[0]))->get_catspecific_fixed();
  unsigned nrcatspecific_fixed=0;
  for(i=0; i<catspecific_fixed.size(); i++)
    {
    if(catspecific_fixed[i])
      nrcatspecific_fixed++;
    }

  for(i=0; i<fullcond.size(); i++)
    {
    catspecific.push_back(fullcond[i]->get_catspecific());
    }

  xcut.push_back(0);
  xcutbeta.push_back(0);

  xcut.push_back(fullcond[0]->get_dimX() + nrcatspecific_fixed*(nrcat2-1));
  xcutbeta.push_back((fullcond[0]->get_dimX()-nrcatspecific_fixed)*nrcat2 + nrcatspecific_fixed);

  zcut.push_back(0);
  zcutbeta.push_back(0);

  for(i=1;i<fullcond.size();i++)
    {
    if(catspecific[i])
      {
      xcut.push_back(xcut[xcut.size()-1]+fullcond[i]->get_dimX()*nrcat2);
      zcut.push_back(zcut[zcut.size()-1]+fullcond[i]->get_dimZ()*nrcat2);

      xcutbeta.push_back(xcutbeta[xcutbeta.size()-1]+fullcond[i]->get_dimX());
      zcutbeta.push_back(zcutbeta[zcutbeta.size()-1]+fullcond[i]->get_dimZ());
      }
    else
      {
      xcut.push_back(xcut[xcut.size()-1]+fullcond[i]->get_dimX());
      zcut.push_back(zcut[zcut.size()-1]+fullcond[i]->get_dimZ());

      for(j=0; j<nrcat2; j++)
        {
        xcutbeta.push_back(xcutbeta[xcutbeta.size()-1]+fullcond[i]->get_dimX());
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

  totalnrfixed = xcutbeta[xcutbeta.size()-1];
  totalnrpar = totalnrfixed + zcutbeta[zcutbeta.size()-1];
  totalvars = zcutbeta.size()-1;

  beta=statmatrix<double>(totalnrpar,1,0);
  for(i=0; i<nrcat2; i++)
    {
    beta(i,0) = randnumbers::uniform();
    }
  theta=statmatrix<double>(totalvars,1,0);

  k=0;
  for(i=1; i<fullcond.size(); i++)
    {
    if(catspecific[i])
      {
      theta(k,0) = fullcond[i]->get_startlambda();
      k++;
      }
    else
      {
      for(j=0; j<nrcat2; j++)
        {
        theta(k,0) = fullcond[i]->get_startlambda();
        k++;
        }
      }
    }
  }

//------------------------------------------------------------------------------
//----------------------------- Estimation -------------------------------------
//------------------------------------------------------------------------------

bool remlest_multinomial_catsp::estimate(const datamatrix resp, const datamatrix & offset,
                const datamatrix & weight, const datamatrix & naindicator)
  {
  unsigned i, j, k, l;

  outoptions();
  out("\n");

  for(i=0;i<fullcond.size();i++)
    fullcond[i]->outoptionsreml();

  vector<int>nasum(nrobs,0);
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      if(naindicator(i,j)==1)
        {
        nasum[i]++;
        }
      }
    }

// ----------------------------------------------------------------------------
// --- Conny: Konstruiere großes X und Z
// ----------------------------------------------------------------------------

 // xcutbetalength gibt die Länge der Abschnitte von xcutbeta an
   vector <unsigned> xcutbetalength (xcutbeta.size()-1);
   for (i=0; i<xcutbetalength.size(); i++)
   {
   xcutbetalength[i]=xcutbeta[i+1]-xcutbeta[i];
   }

 // xcutlength gibt die Länge der Abschnitte von xcut an
   vector <unsigned> xcutlength (xcut.size()-1);
   for (i=0; i<xcutlength.size(); i++)
   {
   xcutlength[i]=xcut[i+1]-xcut[i];
   }

 // Erzeugt Xneu mit 0
    datamatrix Xneu (nrobs*nrcat2,xcutbeta[xcutbeta.size()-1],0);

 // Belegt die ersten Spalten von Xneu
 // Durchlaufen von catspecific_fixed

 // Zwei neue Variablen, beginnend bei 0
    unsigned XaltSpalte_fixed=0;                                                  // gibt die aktuelle Spalte in Xalt an
    unsigned XneuSpalte_fixed=0;                                                  // gibt die aktuelle Spalte in Xneu an

 // j durchläuft den Bereich 0..catspecific_fixed.size()           // catspecific_fixed.size() = xcut[1]

    for (j=0; j<catspecific_fixed.size(); j++)                                    // Durchlauf solange Bedingung vorhanden
    {                                                                             // Anfang zur for-Schleife für X
       // 1. Skalar und catspecific
       // -> nrcat2 x 1 - Matrizen

       if (catspecific_fixed[j])
       {
          for(i=0; i<nrobs; i++)                                                  // i durchläuft die Zeilen von Xalt
          {
              for(k=0; k<nrcat2; k++)                                             // k durchläuft nrcat2
              {
              Xneu(i*nrcat2+k, XneuSpalte_fixed)= X (i, XaltSpalte_fixed + k );
              }
          }

       XaltSpalte_fixed=XaltSpalte_fixed+nrcat2;                                       // aktualisiert XaltSpalte_fixed
       XneuSpalte_fixed=XneuSpalte_fixed+1;                                       // aktualisiert XneuSpalte_fixed
       }

       // 2. Skalar und nicht catspecific
       // -> nrcat2 x nrcat2 - Matrizen

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

 // Schleife, um Xneu zu erweitern
 // Durchlaufen von catspecific

 // Zwei neue Variablen, beginnend bei vektor[1], da vorher catspecific_fixed
    unsigned XaltSpalte=xcut[1];                                                  // gibt die aktuelle Spalte in Xalt an
    unsigned XneuSpalte=xcutbeta[1];                                              // gibt die aktuelle Spalte in Xneu an

 // j durchläuft den Bereich 1..catspecific.size()

    for (j=1; j<catspecific.size(); j++)                                          // Durchlauf solange Bedingung vorhanden
    {                                                                             // Anfang zur for-Schleife für X

       // 1. Vektor und catspecific
       // -> nrcat2 x 1*xcutlength[j] - Matrizen

       if ( catspecific[j])
       {
             for(i=0; i<nrobs; i++)                                               // i durchläuft die Zeilen von Xalt
             {
                for(k=0; k<nrcat2; k++)                                           // k durchläuft nrcat2
                {
                   for(l=0; l<(xcutlength[j]/nrcat2); l++)                                 // l durchläuft xcutbetalength[j]
                    {
                   Xneu(i*nrcat2+k, XneuSpalte+l)= X (i, XaltSpalte+l+(k*xcutlength[j]/nrcat2) );               //+(k*xcutlength[j])
                   }
                }
             }

       XaltSpalte=XaltSpalte+(xcutlength[j]);                                       // aktualisiert XaltSpalte
       XneuSpalte=XneuSpalte+(xcutlength[j]/nrcat2);                                                // aktualisiert XneuSpalte

       }

       // 2. Vektor und nicht catspecific
       // -> nrcat2 x nrcat2*xcutlength[j] - Matrizen

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

 //  Basteln an Z

 // Zu Z: Hilfsvektor bzw. Hilfsmatrix für die if - Schleifen

 // zcutbetalength gibt die Länge der Abschnitte von zcutbeta an
   vector <unsigned> zcutbetalength (zcutbeta.size()-1);
   for (i=0; i<zcutbetalength.size(); i++)
   {
   zcutbetalength[i]=zcutbeta[i+1]-zcutbeta[i];
   }

 // zcutlength gibt die Länge der Abschnitte von zcut an
   vector <unsigned> zcutlength (zcut.size()-1);
   for (i=0; i<zcutlength.size(); i++)
   {
   zcutlength[i]=zcut[i+1]-zcut[i];
   }

 // Erzeugt Zneu mit Nullen
 datamatrix Zneu (nrobs*nrcat2,zcutbeta[zcutbeta.size()-1],0);

 // Zwei neue Variablen, beginnend bei 0
    unsigned ZaltSpalte=0;                                                        // gibt die aktuelle Spalte in Zalt an
    unsigned ZneuSpalte=0;                                                        // gibt die aktuelle Spalte in Zneu an

 // j durchläuft den Bereich 1..catspecific.size()

    for (j=1; j<catspecific.size(); j++)                                          // Durchlauf solange Bedingung vorhanden
    {                                                                             // Anfang zur for-Schleife für Z
       // 1. Vektor und catspecific
       // -> nrcat2 x 1*xcutbetalength[j-1] - Matrizen

      if (catspecific[j])
        {
             for(i=0; i<nrobs; i++)                                               // i durchläuft die Zeilen von Zalt
             {
                for(k=0; k<nrcat2; k++)                                           // k durchläuft nrcat2
                {
                   for(l=0; l< (zcutlength[j-1]/nrcat2); l++)                               // l durchläuft zcutbetalength[j-1]
                   {
                   Zneu(i*nrcat2+k, ZneuSpalte+l)= Z (i, ZaltSpalte+l+(k*zcutlength[j-1]/nrcat2) );
                   }
                }
             }

       ZaltSpalte=ZaltSpalte+zcutlength[j-1];                                     // aktualisiert ZaltSpalte
       ZneuSpalte=ZneuSpalte+(zcutlength[j-1]/nrcat2);                                     // aktualisiert ZneuSpalte
       }

       // 2. Vektor und nicht catspecific                                           sollte ok sein !!!
       // -> nrcat2 x nrcat2*xcutbetalength[j-1] - Matrizen

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

/*
// Berechnen von Xneu' unter XneuTr

   datamatrix XneuTr = Xneu.transposed();

// Berechnen von Zneu' unter ZneuTr

   datamatrix ZneuTr = Zneu.transposed();

// Erzeugen von H_00 = Xneu' workweight Xneu
   datamatrix  H_00 (xcutbeta[xcutbeta.size()-1],  xcutbeta[xcutbeta.size()-1],  0);

// Erzeugen von H_01 = Xneu' workweight Zneu
   datamatrix  H_01 (xcutbeta[xcutbeta.size()-1],  zcutbeta[zcutbeta.size()-1],  0);

// Erzeugen von H1_0 = Xneu' workweight worky
   datamatrix  H1_0 (xcutbeta[xcutbeta.size()-1],  1, 0);

// Hilfsmatrix zum Berechnen der Zeilen von Xneu' workweight
   datamatrix H_00hilf (1, nrobs*nrcat2,  0);

// Erzeugen von H_11 = Zneu' workweight Zneu
   datamatrix  H_11 (zcutbeta[zcutbeta.size()-1],  zcutbeta[zcutbeta.size()-1],  0);

// Erzeugen von H1_1 = Zneu' workweight worky
   datamatrix  H1_1 (zcutbeta[zcutbeta.size()-1],  1, 0);

// Hilfsmatrix zum Berechnen der Zeilen von Zneu' workweight
   datamatrix H_11hilf (1, nrobs*nrcat2,  0);*/

// ----------------------------------------------------------------------------
// --- Ende: Konstruiere großes X und Z
// ----------------------------------------------------------------------------

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
//  statmatrix<double>wresid(nrcat2,resp.rows(),0);                  //matrix containing row vectors !!
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

    eta = offset + Xneu*beta.getRowBlock(0,totalnrfixed)+Zneu*beta.getRowBlock(totalnrfixed,totalnrpar);
    compute_weights(mu,workweight,worky,eta,respind,weight,naindicator,nasum);
    worky = worky - offset;

    stop = check_pause();
    if (stop)
      return true;

/*
// zu Xneu und Zneu: Dimensionen
// von oben:
// datamatrix   Xneu            (   nrobs*nrcat2,    xcutbeta[xcutbeta.size()-1],  0);
// datamatrix   Zneu            (   nrobs*nrcat2,    zcutbeta[zcutbeta.size()-1],  0);
//
// datamatrix   workweight   (   nrobs*nrcat2,     nrcat2,      0);
// datamatrix   worky          (    nrobs*nrcat2,     1,        0);

// Berechnen von H
//
//               Xneu' workweight Xneu           Xneu' workweight Zneu              H_00       H_01
//    H  =                                                                   =
//               Zneu' workweight Xneu           Zneu' workweight Zneu              H_10       H_11

// Berechnen von H1
//
//                 Xneu' workweight worky                 H1_0
//   H1 =                                          =
//                 Zneu' workweight worky                 H1_1

// Berechnen von H_00, H_01, H1_0

for (l=0; l<xcutbeta[xcutbeta.size()-1]; l++ )                                   // l durchläuft die Zeilen von Xneu' und die Zeilen von H_00
{
     for  (j=0;  j<nrobs*nrcat2;  j=j+nrcat2)                                                                                                                      // j durchläuft die Spalten von Xneu'
     {
          for (i=0; i<nrcat2; i++)                                                                                                                                       // i durchläuft nrcat2
          {
          H_00hilf (0, j+i) = ((XneuTr.getBlock( l, j, l+1, j+nrcat2)) * (workweight.getBlock(j, i,  j+nrcat2, i+1)))(0,0);
          }
     }

// Berechnen von H_00 = Xneu' workweight Xneu

     for (k=l; k<xcutbeta[xcutbeta.size()-1]; k++ )                              // k durchläuft die Spalten von Xneu
     {
     H_00(l, k) = ((H_00hilf) * (Xneu.getCol(k)))(0,0);
     H_00(k,l)=H_00(l,k);
     }

// Berechnen von H_01 = Xneu' workweight Zneu
     for (k=0; k<zcutbeta[zcutbeta.size()-1]; k++ )                              // k durchläuft die Spalten von Zneu
     {
     H_01(l, k) = ((H_00hilf) * (Zneu.getCol(k))) (0,0);
     }

// Berechnen von H_10 = Zneu' workweight Xneu
// durch Einsetzen von H_01.transposed() in H

// Berechnen von H1_0 = Xneu' workweight worky
    H1_0(l, 0) = ((H_00hilf) * (worky)) (0,0);

}                                                                                // Ende der l-Schleife

// Berechnen von H_11, H1_1

for (l=0; l<zcutbeta[zcutbeta.size()-1]; l++ )                                   // l durchläuft die Zeilen von Zneu' und die Zeilen von H_11
{

     for  (j=0;  j<nrobs*nrcat2; j=j+nrcat2)                                                                                                                      // j durchläuft die Spalten von Zneu'
     {
          for (i=0; i<nrcat2; i++)                                                                                                                                       // i durchläuft nrcat2
          {
          H_11hilf (0, j+i) = ((ZneuTr.getBlock( l, j, l+1, j+nrcat2)) * (workweight.getBlock(j, i,  j+nrcat2, i+1))) (0,0);
          }
     }

// Berechnen von H_11 = Zneu' workweight Zneu
     for (k=l; k<zcutbeta[zcutbeta.size()-1]; k++ )                              // k durchläuft die Spalten von Zneu
     {
     H_11(l, k) = ((H_11hilf) * (Zneu.getCol(k))) (0,0);
     H_11(k,l)=H_11(l,k);
     }

// Berechnen von H1_1 = Zneu' workweight worky
    H1_1(l, 0) = ((H_11hilf) * (worky)) (0,0);
}                                                                                // Ende der l-Schleife

// Zusammensetzen von H

   H.putBlock (H_00 , 0, 0,
                                       xcutbeta[xcutbeta.size()-1],  xcutbeta[xcutbeta.size()-1] );

   H.putBlock (H_01 , 0, xcutbeta[xcutbeta.size()-1],
                                       xcutbeta[xcutbeta.size()-1],  xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1] );

   H.putBlock (H_01.transposed() , xcutbeta[xcutbeta.size()-1],  0,
                                       xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1] , xcutbeta[xcutbeta.size()-1] );

   H.putBlock (H_11 , xcutbeta[xcutbeta.size()-1],  xcutbeta[xcutbeta.size()-1],
                                       xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1] , xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1] );

// Zusammensetzen von H1

   H1.putRowBlock(0, xcutbeta[xcutbeta.size()-1], H1_0 );
   H1.putRowBlock(xcutbeta[xcutbeta.size()-1], xcutbeta[xcutbeta.size()-1]+zcutbeta[zcutbeta.size()-1], H1_1 );
*/

    compute_sscp_resp2(H1,workweight,worky,Xneu,Zneu);
    compute_sscp2(H,workweight,Xneu,Zneu);

    H.addtodiag(Qinv,totalnrfixed,totalnrpar);

    // Fisher-Scoring for beta
    beta=H.solve(H1);

    // update linear predictor and compute residuals
    eta = Xneu*beta.getRowBlock(0,totalnrfixed)+Zneu*beta.getRowBlock(totalnrfixed,totalnrpar);
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

    for(i=0; i<theta.rows(); i++)
      {
      for(l=zcutbeta[i]; l<zcutbeta[i+1]; l++)
        {
        help = ( wresid*( Zneu.getCol(l) ) )(0,0);
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

    // compute norm of the random parts
    for(i=0; i<theta.rows(); i++)
      {
      helpmat = Zneu.getColBlock(zcutbeta[i],zcutbeta[i+1])*beta.getRowBlock(totalnrfixed+zcutbeta[i],totalnrfixed+zcutbeta[i+1]);
      stopcrit[i] = helpmat.norm(0)/help;
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

// -----------------------------------------------------------------------------
// Thomas: modifizieren
// -----------------------------------------------------------------------------

  // update linear predictor and working weights
  eta = offset + Xneu*beta.getRowBlock(0,totalnrfixed)+Zneu*beta.getRowBlock(totalnrfixed,totalnrpar);
  compute_weights(mu,workweight,worky,eta,respind,weight,naindicator,nasum);

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

  double mean=0;
  k=0;
  for(i=1; i<fullcond.size(); i++)
    {
    if(catspecific[i])
      {
      mean = fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],k,false,xcutbeta[k+1],totalnrfixed+zcutbeta[k],0,false,k+1);
      for(j=0; j<nrcat2; j++)
        {
        beta(j,0) += mean;
        }
      k++;
      }
    else
      {
      out("\n");
      for(j=0; j<nrcat2; j++)
        {
        beta(j,0) += fullcond[i]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[i],zcut[i-1],k,false,xcutbeta[k+1],totalnrfixed+zcutbeta[k],cats(j,0),true,k+1);
        k++;
        }
      }
    }
  beta(0,0) += fullcond[0]->outresultsreml(X,Z,beta,Hinv,thetareml,xcut[0],0,0,false,xcutbeta[0],0,0,false,0);

// compute log-likelihood, AIC, BIC, GCV etc.

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
      if(naindicator(i,j)==0)
        {
        outpredict << eta(i*nrcat2+j,0) << " " << mu(i*nrcat2+j,0) << " ";
        }
      else
        {
        outpredict << "NA" << " " << 0 << " ";
        }
      }
    outpredict << endl;
    }
  outpredict.close();

  return false;
  }

bool remlest_multinomial_catsp::estimate_glm(const datamatrix resp,
                  const datamatrix & offset, const datamatrix & weight,
                  const datamatrix & naindicator)
  {
  unsigned i,j,k,l;

  outoptions();
  out("\n");

  vector<int>nasum(nrobs,0);
  for(i=0; i<nrobs; i++)
    {
    for(j=0; j<nrcat2; j++)
      {
      if(naindicator(i,j)==1)
        {
        nasum[i]++;
        }
      }
    }

  for(i=0;i<fullcond.size();i++)
    fullcond[i]->outoptionsreml();

// ----------------------------------------------------------------------------
// --- Conny: Konstruiere großes X
// ----------------------------------------------------------------------------

// Ausgabedateien gegeben



// Ausgabe von Xalt unter Xalt
/*  ofstream outXalt("c:\\bayesx\\mcmc\\Xalt");
        X.prettyPrint(outXalt);
        outXalt.close();*/

// Hilfsmatrix zu catspecific_fixed
  datamatrix help01(catspecific_fixed.size(),1,0);                          // Anzahl Zeilen= zcutbeta.size(), Anzahl Spalten=1, alle Werte=0
  for(i=0; i<catspecific_fixed.size(); i++)
  {
  help01(i,0) = catspecific_fixed[i];
  }

// Ausgabe von catspecific_fixed unter catspecific_fixed
/*  ofstream outcatspecific_fixed("c:\\bayesx\\mcmc\\catspecific_fixed");
        help01.prettyPrint(outcatspecific_fixed);
        outcatspecific_fixed.close();*/
// ----------------------------------------------------------------------------

// Dimension von Xneu
   unsigned Xneu_rows = nrobs*nrcat2;
   unsigned Xneu_cols = xcutbeta[xcutbeta.size()-1];


/*nicht mehr notwendig
   for(i=0; i<catspecific_fixed.size(); i++)
      if (catspecific_fixed[i])
         {
          Xneu_cols = Xneu_cols + 1
         }
      else
         {
          Xneu_cols = Xneu_cols + nrcat2
         }
*/

   datamatrix Xneu (Xneu_rows,Xneu_cols,0);          // Initialisiert mit 0

// ----------------------------------------------------------------------------

 // Belegt die ersten Spalten von Xneu
// Durchlaufen von catspecific_fixed

// Zwei neue Variablen, beginnend bei 0
   unsigned XaltSpalte_fixed=0;                                                  // gibt die aktuelle Spalte in Xalt an
   unsigned XneuSpalte_fixed=0;                                                  // gibt die aktuelle Spalte in Xneu an

// j durchläuft den Bereich 0..catspecific_fixed.size()           // catspecific_fixed.size() = xcut[1]

   for (j=0; j<catspecific_fixed.size(); j++)                                    // Durchlauf solange Bedingung vorhanden
   {                                                                             // Anfang zur for-Schleife für X

      //----------------------------------
      // 1. Skalar und catspecific
      // -> nrcat2 x 1 - Matrizen

      if (catspecific_fixed[j])
      {
         for(i=0; i<nrobs; i++)                                                  // i durchläuft die Zeilen von Xalt
         {
             for(k=0; k<nrcat2; k++)                                             // k durchläuft nrcat2
             {
             Xneu(i*nrcat2+k, XneuSpalte_fixed)= X (i, XaltSpalte_fixed + k );
             }
         }

      XaltSpalte_fixed=XaltSpalte_fixed+nrcat2;                                       // aktualisiert XaltSpalte_fixed
      XneuSpalte_fixed=XneuSpalte_fixed+1;                                       // aktualisiert XneuSpalte_fixed
      }

      //----------------------------------
      // 2. Skalar und nicht catspecific
      // -> nrcat2 x nrcat2 - Matrizen

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



//-------------------------
// Interessante Matrixoperationen, Matrizen A, B
     // datamatrix A (zeilen, spalten, wert);  z.B. datamatrix A(2, 2, 0);       Anlegen einer Matrix
     // A.diag(dimension, wert);                                                 Diagonalmatrix mit "wert" auf der Diagonalen
     // A.diag(werte);                                                           "werte" als Spaltenvektor mit dem gew³nschten Werten der Diagonalen

// Zugriff auf einzelne Teile
    // A.rows();  Anzahl der Zeilen
    // A.cols();  Anzahl der Spalten
    // A.getBlock(rowfirst, colfirst, rowlast+1, collast+1);                     Blockweiser Elementzugriff

// Besondere Matrizen:
   // A.blockdiag(B);     Blockdiagonalmatrix liefert (A 0
   //                                                  0 B)
   // A.hcat(B);          Horizontale Konkatenation liefert Matrix (A B)           z.B. help2=help2.vcat(help6);
   // A.vcat(B);          Vertikale Konkatenation liefert Matrix (A
   //                                                             B)
   // A.putBlock(B, rowfirst, colfirst, rowlast+1, collast+1);                   ??? siehe Collegeblock
   // A.putColBlock(colfirst, collast+1, B);                                     Vor: B hat genauso viele Zeilen wie A
   // A.putRowBlock(rowfirst, rowlast+1, B);                                     Vor: B hat genauso viele Spalten wie A


   // Matrix X
// von oben:  X = datamatrix(re.rows(),xcut[xcut.size()-1],0);



  // c anfang
  // Schreibt die Matrix X bzw. Z in eine Textdatei
        //ofstream out("c://bayesx/mcmc/conny/test.txt");
        //X.prettyPrint(out);
                                  //oder
                                       //Z.prettyPrint(out);
        //out.close();
  // c ende
  //--------------------------------------------------------------


//  datamatrix Xneu (nrobs*nrcat2,xcutbeta[xcutbeta.size()-1],0);

// ----------------------------------------------------------------------------
// Ausgabe von Xneu unter Xneu
/*     ofstream outXneu("c:\\bayesx\\mcmc\\Xneu");
        Xneu.prettyPrint(outXneu);
        outXneu.close();*/

// ----------------------------------------------------------------------------
// --- Conny: Konstruiere großes X    Ende
// ----------------------------------------------------------------------------

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

/*    ofstream out1("c:\\temp\\offset.raw");
    offset.prettyPrint(out1);
    out1.close();*/

    eta = offset + Xneu*beta;
    compute_weights(mu,workweight,worky,eta,respind,weight,naindicator,nasum);
    worky = worky - offset;

// ----------------------------------------------------------------------------
// --- Conny: Berechne H=Xneu'WXneu und H1=Xneu'Wy
// ----------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------------------------------------------------
// Berechnen von Xneu' unter XneuTr

   datamatrix XneuTr (xcutbeta[xcutbeta.size()-1], nrobs*nrcat2, 0);
   XneuTr = Xneu.transposed();

//----------------------------------------------------------------------
// Erzeugen von H = Xneu' workweight Xneu
   H = datamatrix(xcutbeta[xcutbeta.size()-1],  xcutbeta[xcutbeta.size()-1],  0);

//----------------------------------------------------------------------
// Erzeugen von H1 = Xneu' workweight worky
   H1 = datamatrix(xcutbeta[xcutbeta.size()-1],  1, 0);

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
// Berechnen von H = Xneu' workweight Xneu

     for (k=l; k<xcutbeta[xcutbeta.size()-1]; k++ )                              // k durchläuft die Spalten von Xneu
     {
     H(l, k) = ((H_00hilf) * (Xneu.getCol(k)))(0,0);
     H(k,l)=H(l,k);
     }

//----------------------------------------------------------------------
// Berechnen von H1 = Xneu' workweight worky
    H1(l, 0) = ((H_00hilf) * (worky)) (0,0);
//----------------------------------------------------------------------
}

                  // Ausgabe von H  unter H
/*                     ofstream outH ("c:\\bayesx\\mcmc\\H ");
                     H.prettyPrint(outH );
                     outH.close();*/

                  // Ausgabe von H1  unter H1
/*                     ofstream outH1 ("c:\\bayesx\\mcmc\\H1 ");
                     H1.prettyPrint(outH1 );
                     outH1.close();*/




//    compute_sscp(H,workweight);
//    compute_sscp_resp(H1,workweight,worky);

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

// ----------------------------------------------------------------------------
// --- Thomas: Ergebnisausgabe modifizieren
// ----------------------------------------------------------------------------

  out("\n");
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

  beta(0,0) += fullcond[0]->outresultsreml(X,Z,beta,H,helpmat,xcut[0],0,0,false,xcutbeta[0],0,0,false,0);

  // update eta and working weights
  eta = offset + Xneu*beta;
  compute_weights(mu,workweight,worky,eta,respind,weight,naindicator,nasum);

  loglike=0;
  aic=0;
  bic=0;
  gcv=0;
  df=beta.rows();
  double refprob;
//  unsigned k;

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
      if(naindicator(i,j)==0)
        {
        outpredict << eta(i*nrcat2+j,0) << " " << mu(i*nrcat2+j,0) << " ";
        }
      else
        {
        outpredict << "NA" << " " << 0 << " ";
        }
      }
    outpredict << endl;
    }
  outpredict.close();

  return false;
  }

//------------------------------------------------------------------------------
//------------- Weights, expectation, linear predictor, etc --------------------
//------------------------------------------------------------------------------

void remlest_multinomial_catsp::compute_respind(const datamatrix & re, datamatrix & respind)
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

void remlest_multinomial_catsp::compute_weights(datamatrix & mu, datamatrix & workweights,
                  datamatrix & worky, datamatrix & eta, datamatrix & respind,
                  const datamatrix & weight, const datamatrix & naindicator,
                  const vector<int> & nasum)

  {
  unsigned i,j,k1,k2,l;

  datamatrix workweighthelp;
  datamatrix workweighthelpinverse;

// Compute mu
  datamatrix expos(nrcat2,1,0);
  double exposum;
  for(i=0; i<nrobs; i++)
    {
    expos=datamatrix(nrcat2,1,0);
    exposum=0;
    for(j=0; j<nrcat2; j++)
      {
      if(naindicator(i,j)==0)
        {
        expos(j,0)=exp(eta(i*nrcat2+j,0));
        exposum+=expos(j,0);
        }
      }
    if(naindicator(i,nrcat2)==0)
      {
      exposum+=1;
      }
    for(j=0; j<nrcat2; j++)
      {
      mu(i*nrcat2+j,0)=expos(j,0)/exposum;
      }
    }
/*  datamatrix expos(nrcat2,1,0);
  double exposum;
  for(i=0; i<nrobs; i++)
    {
    expos=datamatrix(nrcat2,1,0);
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
    }*/

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
      if(nasum[i]>0)
        {
//        int test=nasum[i];
        workweighthelp = datamatrix(nrcat2-nasum[i],nrcat2-nasum[i],0);
        workweighthelpinverse = datamatrix(nrcat2,nrcat2,0);
        k1=0;
        for(j=0; j<nrcat2; j++)
          {
          k2=0;
          if(naindicator(i,j)==0)
            {
            for(l=0; l<nrcat2; l++)
              {
              if(naindicator(i,l)==0)
                {
                workweighthelp(k1,k2) = workweights(i*nrcat2+j,l);
                k2++;
                }
              }
            k1++;
            }
          }
        workweighthelp=workweighthelp.inverse();
        k1=0;
        for(j=0; j<nrcat2; j++)
          {
          k2=0;
          if(naindicator(i,j)==0)
            {
            for(l=0; l<nrcat2; l++)
              {
              if(naindicator(i,l)==0)
                {
                workweighthelpinverse(j,l) = workweighthelp(k1,k2);
                k2++;
                }
              }
            k1++;
            }
          }
        worky.putRowBlock(i*nrcat2,(i+1)*nrcat2,eta.getRowBlock(i*nrcat2,(i+1)*nrcat2)+
                        workweighthelpinverse*(respind.getRowBlock(i*nrcat2,(i+1)*nrcat2)-mu.getRowBlock(i*nrcat2,(i+1)*nrcat2)));
        }
      else
        {
        worky.putRowBlock(i*nrcat2,(i+1)*nrcat2,eta.getRowBlock(i*nrcat2,(i+1)*nrcat2)+
                        workweights.getRowBlock(i*nrcat2,(i+1)*nrcat2).inverse()*
                        (respind.getRowBlock(i*nrcat2,(i+1)*nrcat2)-mu.getRowBlock(i*nrcat2,(i+1)*nrcat2)));
        }
      }
    }
/*  for(i=0; i<nrobs; i++)
    {
    if(weight(i,0)>0)
      {
      worky.putRowBlock(i*nrcat2,(i+1)*nrcat2,eta.getRowBlock(i*nrcat2,(i+1)*nrcat2)+
                        workweights.getRowBlock(i*nrcat2,(i+1)*nrcat2).inverse()*
                        (respind.getRowBlock(i*nrcat2,(i+1)*nrcat2)-mu.getRowBlock(i*nrcat2,(i+1)*nrcat2)));
      }
    }*/
  }

void remlest_multinomial_catsp::compute_sscp_resp2(datamatrix & H1, datamatrix & workweight, datamatrix & worky,
                                const datamatrix & Xneu, const datamatrix & Zneu)
  {
  unsigned i,j;
  unsigned xcols = Xneu.cols();
  unsigned zcols = Zneu.cols();
//  H1=datamatrix(H1.rows(),1,0);
  datamatrix weightytemp=datamatrix(worky.rows(),1,0);
  for(i=0; i<nrobs;i++)
    {
    weightytemp.putRowBlock(i*nrcat2,(i+1)*nrcat2,workweight.getRowBlock(i*nrcat2,(i+1)*nrcat2)*worky.getRowBlock(i*nrcat2,(i+1)*nrcat2));
    }
  for(j=0; j<xcols; j++)
    {
    H1(j,0) = (((Xneu.getCol(j)).transposed())*weightytemp)(0,0);
    }
  for(j=0; j<zcols; j++)
    {
    H1(xcols+j,0) = (((Zneu.getCol(j)).transposed())*weightytemp)(0,0);
    }
//  H1.putRowBlock(0,Xneu.cols(),Xneu.transposed()*weightytemp);
//  H1.putRowBlock(Xneu.cols(),Xneu.cols()+Zneu.cols(),Zneu.transposed()*weightytemp);
  }

void remlest_multinomial_catsp::compute_sscp2(datamatrix & H, datamatrix & workweight,
                                const datamatrix & Xneu, const datamatrix & Zneu)
  {
  unsigned i,j,k;
  unsigned xcols = Xneu.cols();
  unsigned zcols = Zneu.cols();
//  H1=datamatrix(H1.rows(),1,0);
  datamatrix temp=datamatrix(Zneu.rows(),1,0);
//  datamatrix temp2;
  for(j=0; j<zcols; j++)
    {
    for(i=0; i<nrobs;i++)
      {
      temp.putRowBlock(i*nrcat2,(i+1)*nrcat2,workweight.getRowBlock(i*nrcat2,(i+1)*nrcat2)*Zneu.getBlock(i*nrcat2,j,(i+1)*nrcat2,j+1));
      }
    for(k=j; k<zcols; k++)
      {
      H(xcols+j,xcols+k) = (((Zneu.getCol(k)).transposed())*temp)(0,0);
      H(xcols+k,xcols+j) = H(xcols+j,xcols+k);
      }
//    temp2 = Zneu.transposed()*temp;
//    H.putBlock(temp2,xcols,xcols+j,xcols+zcols,xcols+j+1);
    }
  for(j=0; j<xcols; j++)
    {
    for(i=0; i<nrobs;i++)
      {
      temp.putRowBlock(i*nrcat2,(i+1)*nrcat2,workweight.getRowBlock(i*nrcat2,(i+1)*nrcat2)*Xneu.getBlock(i*nrcat2,j,(i+1)*nrcat2,j+1));
      }
    for(k=0; k<zcols; k++)
      {
      H(j,xcols+k) = (((Zneu.getCol(k)).transposed())*temp)(0,0);
      H(xcols+k,j) = H(j,xcols+k);
      }
//    temp2 = Zneu.transposed()*temp;
//    H.putBlock(temp2,xcols,j,xcols+zcols,j+1);
//    H.putBlock(temp2.transposed(),j,xcols,j+1,xcols+zcols);

    for(k=j; k<xcols; k++)
      {
      H(j,k) = (((Xneu.getCol(k)).transposed())*temp)(0,0);
      H(k,j) = H(j,k);
      }
//    temp2 = Xneu.transposed()*temp;
//    H.putBlock(temp2,0,j,xcols,j+1);
    }
  }

//------------------------------------------------------------------------------
//----------------------------- Object description -----------------------------
//------------------------------------------------------------------------------

  void remlest_multinomial_catsp::outoptions()
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

void remlest_multinomial_catsp::make_plots(ofstream & outtex,ST::string path_batch,
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

      if(catspecific[j])
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
          if(!catspecific[j])
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
            if(!catspecific[j])
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
              if(!catspecific[j])
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
              if(!catspecific[j])
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
              if(!catspecific[j])
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
              if(!catspecific[j])
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
              if(!catspecific[j])
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
              if(!catspecific[j])
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

void remlest_multinomial_catsp::make_model(ofstream & outtex, const ST::string & rname)
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

void remlest_multinomial_catsp::make_predictor(ofstream & outtex)
  {

  unsigned j;

  ST::string term2 = fullcond[0]->get_term_symbolic();
  ST::string term = "$\\eta^{(j)} = " + term2;    //linearer Prädiktor wird erweitert
  for(j=1;j<fullcond.size();j++)
    {
    term2 = fullcond[j]->get_term_symbolic();
    if(!catspecific[j])
      {
      term2 = term2.insert_after_all_string("^{(j)}","f");;
      }
    term = term + " + " + term2;    //linearer Prädiktor wird erweitert
    }
  outtex << term << "$\\\\\n";
  }

void remlest_multinomial_catsp::make_prior(ofstream & outtex)
  {
  unsigned i,j;
  outtex << "\n\\noindent {\\bf \\large Priors:}\\\\" << endl << "\\\\" << endl;
  for(j=0;j<fullcond.size();j++)
    {
    vector<ST::string> prior = fullcond[j]->get_priorassumptions();
    if(prior.size() != 0)// nur wenn Priors da sind (d.h. Vektor hat Elemente)
      {
      if(fullcond[j]->get_results_type()!="fixed" && !catspecific[j])
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

void remlest_multinomial_catsp::make_options(ofstream & outtex)
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

void remlest_multinomial_catsp::make_fixed_table(ofstream & outtex)
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

void remlest_multinomial_catsp::make_graphics(const ST::string & title,
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

bool remlest_multinomial_catsp::check_pause()
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

void remlest_multinomial_catsp::out(const ST::string & s,bool thick,bool italic,
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


void remlest_multinomial_catsp::outerror(const ST::string & s)
  {
  out(s,true,true,12,255,0,0);
  }







