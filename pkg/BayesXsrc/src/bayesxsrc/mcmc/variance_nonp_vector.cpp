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



#include"variance_nonp_vector.h"



//------------------------------------------------------------------------------
//-- CLASS: FULLCOND_variance_nonp_vector implementation of member functions ---
//------------------------------------------------------------------------------


namespace MCMC
{

using randnumbers::rand_inv_gaussian;

//______________________________________________________________________________
//
// CONSTRUCTOR with Parameters
//
// o                 : pointer to MCMCoptions object
// p                 : pointer to FULLCOND_const object
// d                 : pointer to DISTRIBUTION object
// ti                : reference to title of the full conditional(for example "fixed effects")
// fp                : reference to file path for storing sampled parameters
// fr                : reference to filename for storing results
// shrinkagestart    : Starting value for the shrinkageparameter
// ashrinkage        : Hyperparameter for gammaprior of shrinkageparameter
// bshrinkage        : Hyperparameter for gammaprior of shrinkageparameter
// shrinkagefix     : Should the shrinkageparameter be fixed at the value "shrinkagestart"
// is.ridge          : variable that indicates if L2- or L1- penalty for regressioncoefficients is uesd
// ct                : blocks of regression coefficients
// c                 : reference to responsecategory (important in the case of multivariate response)
//______________________________________________________________________________

FULLCOND_variance_nonp_vector::FULLCOND_variance_nonp_vector(MCMCoptions * o,
                         vector<FULLCOND_const*> & p,DISTRIBUTION * d,
                         const ST::string & ti, const ST::string & fp,
                         const ST::string & fr, const vector<double> & shrinkagestart,
                         const vector<double> & ashrinkage, const vector<double> & bshrinkage,
                         const vector<bool> & shrinkagefix, const vector<double> & shrinkageweight,
//                         const datamatrix start_data, 
                         const vector<bool> & shrinkageadaptive,
                         const bool & isridge, const vector<unsigned> & ct, 
                         const unsigned & c)
                         : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
    {

    fctype = MCMC::variance;

    update_sigma2 = true;

    column = c;

    pathresults = fr;

    Cp = p;

    distrp = d;

    cut = ct;
    is_ridge = isridge;

    priorassumptions.push_back("\\\\");

  // Initialize the matrix for the variances
    datamatrix helpvariances;
    helpvariances = datamatrix(cut[cut.size()-1],1,0);
    for(unsigned i=0; i<cut.size()-1; i++)
      {
      helpvariances.putRowBlock(cut[i],cut[i+1],Cp[i]->get_variances());
      }
    setbeta(helpvariances);
    
    // Startwerte setzen aus Option
//    startdata = start_data;      // Kontrolle der uebergebenen Startwerte
    lassosum = 0;
    ridgesum = 0;
    shrinkage_fix = shrinkagefix[0];
    shrinkage_adaptive = shrinkageadaptive[0];
    weight = shrinkageweight;
    a_shrinkage = ashrinkage;
    b_shrinkage = bshrinkage;


    // Set pathes for the shrinkage  fullcond object
    ST::string path = pathresults.substr(0,pathresults.length()-4)+"_shrinkage.raw";

    // Fullcondobject for shrinkageparameter
    fc_shrinkage = FULLCOND(o,datamatrix(cut[cut.size()-1],1),Cp[0]->get_title()+"_shrinkage",cut[cut.size()-1],1,path);
    fc_shrinkage.setflags(MCMC::norelchange);

    // Set the variablenames of the fullcond object
    vector<ST::string> varnam(nrpar);
    vector<ST::string> helpvarnam;
  
    for(unsigned j=0; j<cut.size()-1; j++)
      {
      helpvarnam = Cp[j]->get_datanames();
      unsigned i = 0;
      for(unsigned k=cut[j]; k<cut[j+1]; k++, i++)
        {
        if(shrinkage_adaptive==true)
          {
          varnam[k] = helpvarnam[i];
          }
        if(shrinkage_adaptive==false)
          {
          varnam[k] = "shrinkage";
          }
        }
      }
    fc_shrinkage.init_names(varnam);


    // Getcurrent value of shrinkageparameter lambda
    double * shrinkagep = fc_shrinkage.getbetapointer();
    for(unsigned j=0; j<nrpar; j++, shrinkagep++)
      {
       * shrinkagep = shrinkagestart[j];
      }
//TEMP:BEGIN--------------------------------------------------------------------
// Kontrolle der uebergebenen Werte
  // vector<ST::string> vnam = fc_shrinkage.get_datanames();
  // double * helpshrinkagep = fc_shrinkage.getbetapointer(); 
  
  // ofstream outstartvectors_vnv("c:/bayesx/test/startvectors_vnv.txt", ios::trunc|ios::app);
  // outstartvectors_vnv << "datamatrix with starting values"<< "\n";
  // outstartvectors_vnv << "effect lambda shrinkage weight a b shrinkagefix adaptive\n";    
  // outstartvectors_vnv << "\n";
  // outstartvectors_vnv << "vektoren"<< "\n" << "varnam lambda tau2 shrinkagep shrinkage weight a b shrinkagefix adaptive" << "\n";
  // for(unsigned i=0; i<shrinkageweight.size(); i++, helpshrinkagep++)
    // {
    // outstartvectors_vnv << varnam[i]  << "  " 
                        // << 1/helpvariances(i,0) << "  " 
                        // << helpvariances(i,0) << "  " 
                        // << * helpshrinkagep  << "  " 
                        // << shrinkagestart[i] << "  " 
                        // << shrinkageweight[i] << "  " 
                        // << ashrinkage[i] << "  "
                        // << bshrinkage[i] << "  " 
                        // << shrinkagefix[i] << "  " 
                        // << shrinkageadaptive[i] << "  " 
                        // << "\n";
    // }

//TEMP:END----------------------------------------------------------------------  

  }
    
//______________________________________________________________________________
//
// COPY CONSTRUCTOR
//______________________________________________________________________________

FULLCOND_variance_nonp_vector::FULLCOND_variance_nonp_vector(const FULLCOND_variance_nonp_vector & t)
    : FULLCOND(FULLCOND(t))
  {
  update_sigma2 = t.update_sigma2;
  column = t.column;
  pathresults = t.pathresults;
  Cp = t.Cp;
  distrp = t.distrp;
  fc_shrinkage = t.fc_shrinkage;
  shrinkage_fix = t.shrinkage_fix;
  shrinkage_adaptive = t.shrinkage_adaptive;
  a_shrinkage = t.a_shrinkage;
  b_shrinkage = t.b_shrinkage;
  startdata = t.startdata;
  weight = t.weight;
  lassosum = t.lassosum;
  ridgesum = t.ridgesum;
  cut = t.cut;
  is_ridge = t.is_ridge;
  }


//______________________________________________________________________________
//
// OVERLOADED ASSIGNMENT OPERATOR
//______________________________________________________________________________

const FULLCOND_variance_nonp_vector & FULLCOND_variance_nonp_vector::operator=(
                                        const FULLCOND_variance_nonp_vector & t)
  {
  if (this == &t)
    return *this;
  FULLCOND::operator=(FULLCOND(t));
  update_sigma2 = t.update_sigma2;
  column = t.column;
  pathresults = t.pathresults;
  Cp = t.Cp;
  distrp = t.distrp;
  fc_shrinkage = t.fc_shrinkage;
  shrinkage_fix = t.shrinkage_fix;
  shrinkage_adaptive = t.shrinkage_adaptive;
  a_shrinkage = t.a_shrinkage;
  b_shrinkage = t.b_shrinkage;
  startdata = t.startdata;
  weight = t.weight;
  lassosum = t.lassosum;
  ridgesum = t.ridgesum;
  cut = t.cut;
  is_ridge = t.is_ridge;
  return *this;
  }


// Pointer auf das Shrinkagearameter lambda Fullcond-Objekt
FULLCOND * FULLCOND_variance_nonp_vector::get_shrinkagepointer()
  {
  return &fc_shrinkage;
  }

//______________________________________________________________________________
//
// FUNCTION: update
// TASK: - stores sampled parameters in file 'samplepath'
//         storing order: first row, second row, ...
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::update(void)
  {
  acceptance++;
  unsigned i, j, k;
  
  double rand_invgaussian = 0;
  double helpshrinkage;
    
  // variables for summs
  lassosum = 0;
  ridgesum = 0;
  double sumvariances = 0;
  double sumregcoeff = 0;
  double sumabsregcoeff = 0;

  // get current value of (first)
  int iteration = optionsp->get_nriter();                  // Iteration
  double * shrinkagep = fc_shrinkage.getbetapointer();     // shrinkagearameter
  
//TEMP:BEGIN--------------------------------------------------------------------
//ofstream outputl("c:/bayesx/test/test_lasso_vnv.txt", ios::out|ios::app);
// ofstream outputr("c:/bayesx/test/test_ridge_vnv.txt", ios::out|ios::app);
// outputl << "Penalty  Termnr  Iteration  Beta   Var  Weight Shrp" << "\n";
// outputr << "Penalty  Termnr  Iteration  Beta   Var  Weight Shrp" << "\n";
// int nrridge = distrp->get_ridge();
// int nrlasso = distrp->get_lasso();
//TEMP:END----------------------------------------------------------------------

  
  double * workbeta;                                       // regressioncoefficients
  double * workvariance;                                   // variance
  double sigma2 = distrp->get_scale(column);                // scale parameter


  // Datei-Ausgabe der Startwerte
  //------------------------------
  if(iteration==1)
    {
    get_startvalues(); 
    }




  // Gibbs-Update of varianceparameters 1/tau^2: with Inverse Normaldistribution
  // if L1-penalty is used
  //---------------------------------------------------------------------------
  if (is_ridge == 0)                              // L1-penalty
    {
    //reset shrinkagepointer
    shrinkagep = fc_shrinkage.getbetapointer();
    k = 0;
    for(j=0; j<cut.size()-1; j++)
      {
      workbeta = Cp[j]->getbetapointer();               // current value of first regressionparameter
      for(k=cut[j]; k<cut[j+1]; k++, workbeta++, shrinkagep++)
        {
        if (*workbeta>0 && *shrinkagep>0)
          {
          rand_invgaussian = rand_inv_gaussian((sqrt(sigma2) * *shrinkagep)/(*workbeta), 
                                               (*shrinkagep * *shrinkagep));
          beta(k,0) = 1.0/rand_invgaussian;
          }
        if (*workbeta<0 && *shrinkagep>0)
          {
          rand_invgaussian = rand_inv_gaussian(-1.0 * (sqrt(sigma2) * *shrinkagep)/(*workbeta), 
                                              (*shrinkagep * *shrinkagep));
          beta(k,0) = 1.0/(rand_invgaussian);
          }
        if (*workbeta==0 || *shrinkagep<=0)
          {
          beta(k,0) = 1E-6;
          }
        lassosum = lassosum + ((*workbeta)*(*workbeta))/beta(k,0);              // sum(beta^2/tau^2) for scale update
        sumabsregcoeff = sumabsregcoeff + sqrt((*workbeta)*(*workbeta));              // sum(beta^2/tau^2) for scale update
        sumvariances = sumvariances + beta(k,0);        // sum(tau^2) of current variances
        }
      }
     
    // transfer lassosum to scale update
    distrp->update_lasso(lassosum);
    }

  if (is_ridge == 1)                              // L2-penalty
    {
    //reset shrinkagepointer
    shrinkagep = fc_shrinkage.getbetapointer();
    k = 0;
    for(j=0; j<cut.size()-1; j++)
      {
      workbeta = Cp[j]->getbetapointer();               // current value of first regressionparameter
      for(k=cut[j]; k<cut[j+1]; k++, workbeta++, shrinkagep++)   
        {

         beta(k,0) = 1/(2 * *shrinkagep);
         ridgesum = ridgesum + ((*workbeta)*(*workbeta))/beta(k,0);  // sum(beta^2/tau^2)
         sumregcoeff = sumregcoeff + (*workbeta * *workbeta);
        }
      }
    // transfer ridgesum to scale update
    distrp->update_ridge(ridgesum);
    }
    


  // Gibbs-Update of Shrinkageparameter with Gammadistribution
  //----------------------------------------------------------
  if(shrinkage_fix==false)
    {
    // reset shrinkagepointer
    shrinkagep = fc_shrinkage.getbetapointer();
    
    if(is_ridge == 0 && shrinkage_adaptive == false)            // L1-penalty
      {
      //helpshrinkage = rand_gamma(nrpar + a_shrinkage[0], b_shrinkage[0] + sumabsregcoeff/sqrt(sigma2));//Armagan+Dunson
      helpshrinkage =  sqrt(rand_gamma(nrpar + a_shrinkage[0], b_shrinkage[0] + 0.5*sumvariances));
      for(i=0; i<nrpar; i++, shrinkagep++)
        {
        *shrinkagep = helpshrinkage;
        }
      }
    if(is_ridge == 0 && shrinkage_adaptive == true)            // L1-penalty adaptive
      {
      k = 0;
      for(j=0; j<cut.size()-1; j++)
        {
        // reset betapointer
        workbeta = Cp[j]->getbetapointer();               // current value of first regressionparameter
        for(k=cut[j]; k<cut[j+1]; k++, workbeta++, shrinkagep++)
          {
          //*shrinkagep = rand_gamma(1 + a_shrinkage[k], b_shrinkage[k] + abs(*workbeta)/sqrt(sigma2));//Armagan+Dunson
          *shrinkagep = sqrt(rand_gamma(1 + a_shrinkage[k], b_shrinkage[k] + 0.5*beta(k,0)));
          }
        }
      }
      
    if(is_ridge == 1 && shrinkage_adaptive == false)            // L2-penalty
      {
      helpshrinkage = rand_gamma(0.5*nrpar + a_shrinkage[0], b_shrinkage[0] + sumregcoeff/sigma2);
      for(i=0; i<nrpar; i++, shrinkagep++)
        {
        *shrinkagep = helpshrinkage;
        }
      }
    if(is_ridge == 1 && shrinkage_adaptive == true)            // L2-penalty adaptive
      {  
      k = 0;
      for(j=0; j<cut.size()-1; j++)
        {
        // reset betapointer
        workbeta = Cp[j]->getbetapointer();               // current value of first regressionparameter
        for(k=cut[j]; k<cut[j+1]; k++, workbeta++, shrinkagep++)
          {
          *shrinkagep = rand_gamma(0.5 + a_shrinkage[k], b_shrinkage[k] + ((*workbeta)*(*workbeta))/sigma2);
          }
        }
      }
    }
    
  // Transfer of the updates
  //------------------------
  k = 0;
  for(j=0; j<cut.size()-1; j++)                     // variances
    {
    workvariance = Cp[j]->getvariancespointer();
    for(k=cut[j]; k<cut[j+1]; k++, workvariance++)
      {
      *workvariance = beta(k,0);
      }
    }

  fc_shrinkage.update();                           // shrinkageparameter
  fc_shrinkage.acceptance++;

  FULLCOND::update();
  }



//______________________________________________________________________________
//
// FUNCTION: outresults
// TASK: - write results for varianceparameters to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::outresults(void)
  {
  FULLCOND::outresults();
  
  unsigned int i,j,k;
  vector<ST::string> vnames(nrpar);
  vector<ST::string> helpvarnames;

  for(j=0; j<cut.size()-1; j++)
    {
    helpvarnames = Cp[j]->get_datanames();
    i = 0;
    for(k=cut[j]; k<cut[j+1]; k++, i++)
      {
      vnames[k] = helpvarnames[i];
      }
    }

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  // Dateiausgabe
  ofstream outp(pathresults.strtochar());

  // Kopfzeile Dateiausgabe
  if (pathresults.isvalidfile() != 1)
    outp << "paramnr varname pmean pstd pqu" << l1 << " pqu" << l2 <<
           " pmed pqu" << u1 << " pqu" << u2 << " pcat" << level1 << 
           " pcat" << level2 << endl;

  optionsp->out("\n");

  ST::string l;
  int maxvarnamelength = 0;
  int len;

  for(i=0;i<nrpar;i++)
    {
    len = vnames[i].length();
    if (len > maxvarnamelength)
      maxvarnamelength = len;
    }

  if (maxvarnamelength>10)
    l = ST::string(' ',maxvarnamelength-4);
  else
    l = "  ";

  ST::string help =  ST::doubletostring(lower1,4) + "% quant.";
  ST::string levell = help + ST::string(' ',15-help.length());
  help = ST::doubletostring(upper2,4) + "% quant.";
  ST::string levelu = help + ST::string(' ',15-help.length());

  // Kopfzeile Bildschirmausgabe
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

  for (i=0;i<nrpar;i++)
    {
    if (maxvarnamelength  > 10)
      nsp = 4+maxvarnamelength-vnames[i].length();
    else
      nsp = 10-vnames[i].length();

    m = betamean(i,0);

    if (betavar(i,0) == 0)
      stddouble = 0;
    else
      stddouble = sqrt(betavar(i,0));
     
    // Dateiausgabe
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

      // Bildschirmausgabe
      optionsp->out(ST::outresults(nsp,vnames[i],m,stddouble,
                    betaqu_l1_lower(i,0),betaqu50(i,0),betaqu_l1_upper(i,0)) + "\n");

      }

    }

  optionsp->out("\n");
  optionsp->out("  Results for the variances are also stored in file\n");
  optionsp->out("  " + pathresults + "\n");
  optionsp->out("\n");

  // Ausgabe der restlichen Ergebnisse
  outresults_shrinkage();
  }


//______________________________________________________________________________
//
// FUNCTION: outresults_shrinkage
// TASK: - write results for shrinkaageparameter to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::outresults_shrinkage(void)
  {

  if(shrinkage_fix==false)
  {
  fc_shrinkage.outresults();
  
  ST::string shrinkage_pathresults = pathresults.substr(0,pathresults.length()-7) + "shrinkage.res";

  unsigned nr;
  if(shrinkage_adaptive==true)
    {
    nr = nrpar;
    }
  if(shrinkage_adaptive==false)
    {
    nr =1;
    }
  unsigned i;

  vector<ST::string> vnames = fc_shrinkage.get_datanames();
  
  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  // Dateiausgabe
  ofstream outp(shrinkage_pathresults.strtochar());

  // Kopfzeile Dateiausgabe
  if (shrinkage_pathresults.isvalidfile() != 1)
    outp << "paramnr varname pmean pstd pqu" << l1 << " pqu" << l2 <<
           " pmed pqu" << u1 << " pqu" << u2 << " pcat" << level1 << 
           " pcat" << level2 << endl;

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

  ST::string help =  ST::doubletostring(lower1,4) + "% quant.";
  ST::string levell = help + ST::string(' ',15-help.length());
  help = ST::doubletostring(upper2,4) + "% quant.";
  ST::string levelu = help + ST::string(' ',15-help.length());

  // Kopfzeile Bildschirmausgabe
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

    m = fc_shrinkage.get_betamean(i,0);

    if (fc_shrinkage.get_betavar(i,0) == 0)
      stddouble = 0;
    else
      stddouble = sqrt(fc_shrinkage.get_betavar(i,0));
     
    // Dateiausgabe
    if (shrinkage_pathresults.isvalidfile() != 1)
      {
      outp << (i+1) << "   ";
      outp << vnames[i] << "   ";
      outp << m << "   ";
      outp << stddouble << "   ";

      outp << fc_shrinkage.get_beta_lower1(i,0) << "   ";
      outp << fc_shrinkage.get_beta_lower2(i,0) << "   ";
      outp << fc_shrinkage.get_betaqu50(i,0) << "   ";
      outp << fc_shrinkage.get_beta_upper2(i,0) << "   ";
      outp << fc_shrinkage.get_beta_upper1(i,0) << "   ";
      if (fc_shrinkage.get_beta_lower1(i,0) > 0)
        outp << "1   ";
      else if (fc_shrinkage.get_beta_upper1(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";
      if (fc_shrinkage.get_beta_lower2(i,0) > 0)
        outp << "1   ";
      else if (fc_shrinkage.get_beta_upper2(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";
      outp << endl;

      // Bildschirmausgabe
      optionsp->out(ST::outresults(nsp,vnames[i],m,stddouble,
                    fc_shrinkage.get_beta_lower1(i,0),
                    fc_shrinkage.get_betaqu50(i,0),
                    fc_shrinkage.get_beta_upper1(i,0)) + "\n");

      }

    }

  optionsp->out("\n");
  optionsp->out("  Results for the shrinkage parameter are also stored in file\n");
  optionsp->out("  " + shrinkage_pathresults + "\n");
  optionsp->out("\n");

  }
}

//______________________________________________________________________________
//
// FUNCTION: outoptions
// TASK: - write options to output window
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::outoptions(void)
  {
  
  vector<ST::string> varnames(nrpar);
  vector<ST::string> helpvarnames;

  for(unsigned j=0; j<cut.size()-1; j++)
    {
    helpvarnames = Cp[j]->get_datanames();
    unsigned i = 0;
    for(unsigned k=cut[j]; k<cut[j+1]; k++, i++)
      {
      varnames[k] = helpvarnames[i];
      }
    }

  
  vector <ST::string> helppath;
  for(unsigned j=0; j<cut.size()-1; j++)
    {
    helppath.push_back(Cp[j]->get_title());
    }  
  optionsp->out("  OPTIONS FOR SHRINKAGE EFFECTS: " + helppath[0] + "\n",true);
  if(helppath.size()>1)
    { 
    for(unsigned j=1; j<helppath.size(); j++)
      {
       optionsp->out( ST::string(' ',33) + helppath[j] + "\n",true);
      }
    }
  optionsp->out("\n");
  if(is_ridge == 0)
    {
    optionsp->out("  Priors: LASSO shrinkage priors\n");
    }
  if(is_ridge == 1)
    {
    optionsp->out("  Priors: RIDGE shrinkage priors\n");
    }
    optionsp->out("\n");

  if(shrinkage_adaptive == false && shrinkage_fix==false)
    {
    optionsp->out("  Hyperparameter a for shrinkage: " +
                     ST::doubletostring(a_shrinkage[0]) + "\n" );
    optionsp->out("  Hyperparameter b for shrinkage: " +
                     ST::doubletostring(b_shrinkage[0]) + "\n" );
    }
    
  if(shrinkage_adaptive == false && shrinkage_fix==true)
    {
    optionsp->out("  Shrinkage is fixed at value: " +
                     ST::doubletostring(fc_shrinkage.getbeta(0,0)) + "\n" );
    }
    
  if(shrinkage_adaptive == true && shrinkage_fix==false)
    {
    for(unsigned i=0; i<nrpar; i++)
      {
      optionsp->out("  Hyperparameter a for shrinkage of " +  varnames[i] + ": " +
                       ST::doubletostring(a_shrinkage[i]) + "\n" );
      optionsp->out("  Hyperparameter b for shrinkage of " +  varnames[i] + ": " +
                       ST::doubletostring(b_shrinkage[i]) + "\n" );
      optionsp->out("\n");
      }
    }
    
  if(shrinkage_adaptive == true && shrinkage_fix==true)
    {
    for(unsigned i=0; i<nrpar; i++)
      {
      optionsp->out("  Shrinkage of " + varnames[i] + " is fixed at value: " +
                       ST::doubletostring(fc_shrinkage.getbeta(i,0)) + "\n" );
      }
    }
    
    optionsp->out("\n");
  }
  
//______________________________________________________________________________
//
// FUNCTION: get_samples
// TASK: - write samples to files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::get_samples(const ST::string & filename,
                                                const unsigned & step) const
  {
  FULLCOND::get_samples(filename, step);

  ST::string pathhelp = pathresults.substr(0,pathresults.length()-7)+"shrinkage_sample.raw";

  optionsp->out(pathhelp + "\n");
  optionsp->out("\n");
  fc_shrinkage.get_samples(pathhelp);
  }
  
//______________________________________________________________________________
//
// FUNCTION: get_startvalues
// TASK: - write startvalues to files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector::get_startvalues(void)
  {
  unsigned i, j, k;
    
  // Hilfsvariablen
//  double * helpbeta;
  double * helpvariance;
  double * helpshrinkage = fc_shrinkage.getbetapointer();
  vector<ST::string> helpvarname_variance;
  vector<ST::string> helpvarname_shrinkage = fc_shrinkage.get_datanames();
  
  // Variances
  ST::string variance_pathstartdata = pathresults.substr(0,pathresults.length()-7) + "variance_startdata.raw";
  ofstream variance_outoutstartdata(variance_pathstartdata.strtochar());
  variance_outoutstartdata << "varname startvalue" << endl;
  
  // Shrinkageparameter
  ST::string shrinkage_pathstartdata = pathresults.substr(0,pathresults.length()-7) + "shrinkage_startdata.raw";
  ofstream shrinkage_outoutstartdata(shrinkage_pathstartdata.strtochar());
  shrinkage_outoutstartdata << "varname startvalue" << endl;

  // Hyperparameter
  ST::string hyperpar_pathstartdata = pathresults.substr(0,pathresults.length()-7) + "hyperpar_startdata.raw";
  ofstream hyperpar_outoutstartdata(hyperpar_pathstartdata.strtochar());
  hyperpar_outoutstartdata << "varname weight a b shrinkagefix adptive" << endl;
  
  for(j=0; j<cut.size()-1; j++)
    {
    helpvariance = Cp[j]->getvariancespointer();
    helpvarname_variance = Cp[j]->get_datanames();
    i = 0;
    for(k=cut[j]; k<cut[j+1]; k++, i++, helpvariance++, helpshrinkage++)
      {
      variance_outoutstartdata  << helpvarname_variance[i] << " " << *helpvariance << endl;
      shrinkage_outoutstartdata << helpvarname_variance[i] << " " << *helpshrinkage << "  " << endl;
      hyperpar_outoutstartdata  << helpvarname_variance[i]  << " " 
                                << weight[k]  << " " 
                                << a_shrinkage[k] << " " 
                                << b_shrinkage[k] << " " 
                                << shrinkage_fix << " " 
                                << shrinkage_adaptive << endl;
      }
    }
  }
  
} // end: namespace MCMC





