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



#include "mcmc_ridge.h"

using std::ios;

//------------------------------------------------------------------------------
//------------ CLASS: FULLCOND_ridge implementation of member functions --------
//------------------------------------------------------------------------------



namespace MCMC
{

//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------
//______________________________________________________________________________
//
// CONSTRUCTOR with Parameters
//
// o         : pointer to MCMCoptions object
// dp        : pointer to DISTRIBUTION object
// d         : reference to datamatrix - designmatrix X
// t         : reference to title of the full conditional
//             (for example "fixed effects")
// fs        : reference to file path for storing sampled parameters
// fr        : reference to filename for storing results
// vars      : reference to vector of variances for the ridge penalty
// c         : reference to responsecategory
//             (important in the case of multivariate response)
// d.cols()  : number of rows of the beta matrix (i.e. number of parameters)
// 1         : number of columns of the beta matrix
//             (i.e. number of categories of the response variable)
//______________________________________________________________________________

FULLCOND_ridge::FULLCOND_ridge(MCMCoptions * o, DISTRIBUTION * dp,
                const datamatrix & d, const ST::string & t,
                const ST::string & fs, const ST::string & fr,
                const vector<double> & vars, const unsigned & c)
  : FULLCOND(o,d,t,d.cols(),1,fs)
  {

  // initialize attributes of the class FULLCOND
  //--------------------------------------------
  nrpar = d.cols();
  column = c;
  pathresult = fr;
  pathcurrent = fr;
  setbeta(d.cols(),1,0);                         //initialize beta matrices

  // initialize attributes of the class FULLCOND_ridge
  //--------------------------------------------------
  variances = vars;
  lasso = 0.01;                                                              //SET!
  likep = dp;
  linold = datamatrix(likep->get_nrobs(),1,0);
  mu1 = datamatrix(likep->get_nrobs(),1,0);
  XX = d.transposed()*d;
  X1 = datamatrix(d.cols(),d.cols(),0);
  X2 = datamatrix(d.cols(),likep->get_nrobs(),0);

//-Temorär ---------------------------------------------------------------------
  // Outputmatrix für Varianzschätzer und Lassoschätzer und Hyperparameter
  hypr = 0.01;                                                           //SET!
  hyps = 0.01;                                                           //SET!
  estimVariances = datamatrix(1,d.cols(),0);
  estimLasso = datamatrix(1,1,0);

  //Ausgabe der Startwerte
  ofstream outpHyperPar("c:/bayesx/test/results/HyperPar.txt", ios::out);
  outpHyperPar << hypr << ' ' << hyps << endl;

  ofstream outpInitialVal("c:/bayesx/test/results/InitialVal.txt", ios::out);
  outpInitialVal << lasso << ' ';
  for(unsigned int i=0; i<nrpar; i++)
    {
    outpInitialVal << variances[i] << ' ';
    }
  outpInitialVal << endl;
//-Temorär Ende-----------------------------------------------------------------
  }


//______________________________________________________________________________
//
// COPY CONSTRUCTOR
//______________________________________________________________________________

FULLCOND_ridge::FULLCOND_ridge(const FULLCOND_ridge & m)
  : FULLCOND(FULLCOND(m))
  {
  variances = m.variances;
  lasso = m.lasso;
  likep = m.likep;
  linold = m.linold;
  mu1 = m.mu1;
  XX = m.XX;
  X1 = m.X1;
  X2 = m.X2;
//-Temorär ---------------------------------------------------------------------
  hypr = m.hypr;
  hyps = m.hyps;
  estimVariances = m.estimVariances;
  estimLasso = m.estimLasso;
//-Temorär Ende-----------------------------------------------------------------
  }


//______________________________________________________________________________
//
// OVERLOADED ASSIGNMENT OPERATOR
//______________________________________________________________________________

const FULLCOND_ridge & FULLCOND_ridge::operator=(const FULLCOND_ridge & m)
  {
  if (this == &m)
    return *this;
  FULLCOND::operator=(FULLCOND(m));

  variances = m.variances;
  lasso = m.lasso;
  likep = m.likep;
  linold = m.linold;
  mu1 = m.mu1;
  XX = m.XX;
  X1 = m.X1;
  X2 = m.X2;
//-Temorär ---------------------------------------------------------------------
  hypr = m.hypr;
  hyps = m.hyps;
  estimVariances = m.estimVariances;
  estimLasso = m.estimLasso;
//-Temorär Ende-----------------------------------------------------------------

  return *this;
  }


//-------------------------- UPDATE and related methods-------------------------
//______________________________________________________________________________
//
// FUNCTION: update
// TASK: - stores sampled parameters in file 'samplepath'
//         storing order: first row, second row, ...
// STATUS: In Bearbeitung
//______________________________________________________________________________

void FULLCOND_ridge::update(void)
  {

  // Gibbs-Update of parameter beta with Normaldistribution
  //--------------------------------------------------------

  // stores sampled betaparameters in file 'samplepath'
  // computes samplemean, samplevariance,...
  // storing order: first row, second row, ...
  FULLCOND::update();


  // DO: likep.'linpred_current'[,column]=likep.'linpred_current'[,column]-linold
  // RESULT: likep.'linpred_current'= PredictorpartNotUpdatedHere
  likep->substr_linearpred_m(linold,column);


  // DO: mu1=likep.'response'[,column]-likep.'linpred_current'[,column]
  // RESULT: mu1=response-PredictorpartNotUpdatedHere
  likep->compute_respminuslinpred(mu1,column);


  // Gibbs-Update of parameter beta
  // RESULT: beta=X2*mu1 + sqrt(sigma)*X1*rand_normvek(nrpar)
  for(unsigned int i=0; i<nrpar; i++)
    {
    XX(i,i) = XX(i,i) + likep->get_scale(column)/variances[i];
    }
  X1 = (XX.cinverse()).root();                      // X1=XX^-0.5
  X2 = XX.cinverse()*(data.transposed());           // X2=XX^-1 * X'
  beta.mult(X2,mu1);
  beta+= sqrt(likep->get_scale(column))*X1*rand_normvek(nrpar);


  // UpdatedPredictorpart=linold=X*beta
  linold.mult(data,beta);


  // DO: likep.'linpred_current'[,column]=likep.'linpred_current'[,column]+linold
  // RESULT: likep.'linpred_current'=PredictorpartNotUpdatedHere + UpdatedPredictorpart
  likep->add_linearpred_m(linold,column);


  // number of accepted iterations
  acceptance++;


  // transform: factor with which all beta's will be multiplied before storing
  // trmult: multiplicative constant with which response has been transformed
  transform = likep->get_trmult(column);



  // Gibbs-Update of parameter lasso with Gammadistribution
  //--------------------------------------------------------
  double sumvariances = 0;
  for(unsigned int i=0; i<nrpar; i++)
    {
    sumvariances = sumvariances + variances[i];
    }
  lasso = sqrt(rand_gamma(nrpar + hypr, hyps + 0.5*sumvariances));
//-Temorär----------------------------------------------------------------------
// Update der Outputmatrix für Lassoschätzer
  estimLasso = estimLasso.vcat(datamatrix(1,1,lasso));
//-Temorär Ende-----------------------------------------------------------------


  // Gibbs-Update of parameters 1/tau^2: with Inverse Normaldistribution
  //--------------------------------------------------------------------

  double* workbeta = beta.getV();
  datamatrix matrixvariances = datamatrix(1,nrpar,0);
  for(unsigned int i=0; i<nrpar; i++, workbeta++)
    {
    if (*workbeta>0)
      {
       variances[i] = 1.0/rand_inv_gaussian(lasso/(*workbeta), lasso*lasso);
      }
   else
      {
      variances[i] = 1.0/rand_inv_gaussian(-1.0*lasso/(*workbeta), lasso*lasso);
      }
   matrixvariances(0,i) = variances[i];
    }
//-Temorär----------------------------------------------------------------------
// Update der Outputmatrix für Varianzschätzer und Lassoschätzer
  estimVariances = estimVariances.vcat(matrixvariances);
//-Temorär Ende-----------------------------------------------------------------


  // Bedingugng zur Erzeugung des Ouutputs
  //--------------------------------------
  if  (optionsp->get_nriter() == optionsp->get_iterations())
    {
    FULLCOND::outresults();
    }

  }



//______________________________________________________________________________
//
// FUNCTION: outresults
// TASK: - write results to output window and files
// Status: Noch nicht bearbeitet
//______________________________________________________________________________

void FULLCOND_ridge::outresults(void)
  {
  FULLCOND::outresults();

//-Temorär----------------------------------------------------------------------
  // Output der Matrizen für Varianzschätzer und Lassoschätzer
  ofstream outpLasso("c:/bayesx/test/results/Lasso.txt", ios::out);
  ofstream outpVariances("c:/bayesx/test/results/Variances.txt", ios::out);
  estimLasso.prettyPrint(outpLasso);
  estimVariances.prettyPrint(outpVariances);
//-Temorär Ende-----------------------------------------------------------------


  ofstream outp(pathcurrent.strtochar());

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  vector<ST::string> resultstable(6);

  outp << "paramnr varname pmean pstd pqu" << l1 << " pqu" << l2 <<
            " pmed pqu" << u1 << " pqu" << u2 << " pcat" << level1
            << " pcat" << level2 << endl;

  unsigned i;

  optionsp->out("\n");

  ST::string l;
  int maxvarnamelength = 0;
  int len;

  for(i=0;i<nrpar;i++)
    {
    len = datanames[i].length();
    if (len > maxvarnamelength)
      maxvarnamelength = len;
    }

  if (maxvarnamelength>10)
    l = ST::string(' ',maxvarnamelength-6);
  else
    l = "  ";

    ST::string help =  ST::doubletostring(lower1,4) + "% quant.";
    ST::string levell = help + ST::string(' ',15-help.length());
    help = ST::doubletostring(upper2,4) + "% quant.";
    ST::string levelu = help + ST::string(' ',15-help.length());

    optionsp->out("  Variable" + l +
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

    for (i=0;i<nrpar;i++)
      {

      if (maxvarnamelength  > 10)
        nsp = 2+maxvarnamelength-datanames[i].length();
      else
        nsp = 10-datanames[i].length();

      m= betamean(i,0);

      if (betavar(i,0) == 0)
        stddouble = 0;
      else
        stddouble = sqrt(betavar(i,0));

      outp << (i+1) << "   ";
      outp << datanames[i] << "   ";
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

      optionsp->out(ST::outresults(nsp,datanames[i],m,
                      stddouble,betaqu_l1_lower(i,0),
                      betaqu50(i,0),betaqu_l1_upper(i,0)) + "\n");


      char hchar = '_';
      ST::string hstring = "\\_";
      resultstable[0] = datanames[i].insert_string_char(hchar,hstring);
      resultstable[1] = ST::doubletostring(m,6);
      resultstable[2] = ST::doubletostring(stddouble,6);
      resultstable[3] = ST::doubletostring(betaqu_l1_lower(i,0),6);
      resultstable[4] = ST::doubletostring(betaqu50(i,0),6);
      resultstable[5] = ST::doubletostring(betaqu_l1_upper(i,0),6);

      results_latex.push_back(ST::make_latextable(resultstable));

      }

    optionsp->out("\n");

    optionsp->out("  Results for shrinked effects are also stored in file\n");
    optionsp->out("  " + pathcurrent + "\n");

    optionsp->out("\n");

  }


}

// end: namespace MCMC
