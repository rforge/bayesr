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




#include"variance_nonp_vector_nigmix.h"



//------------------------------------------------------------------------------
//-- CLASS: FULLCOND_variance_nonp_vector_nigmix implementation of member functions ---
//------------------------------------------------------------------------------

namespace MCMC
{

using randnumbers::rand_inv_gaussian;
using randnumbers::bernoulli;
using randnumbers::rand_beta;

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
// indicators        : Starting value for varianceparametercomponent indicator
// v0                : Point for Indicatorfunction
// v1                : Point for Indicatorfunction
// t2s               : Starting value for varianceparametercomponent t2
// at2               : Hyperparameter for inverse gammaprior of vainaceparametercomponent t2
// bt2               : Hyperparameter for inverse gammaprior of vainaceparametercomponent t2
// aomega;           : Hyperparameter for Prameter w
// bomega;           : Hyperparameter for Prameter w
// omegas            : Starting value for omega
// omegaf            : Should the omega be fixed at the value "omegas"    
// ct                : blocks of regression coefficients
// c                 : reference to responsecategory (important in the case of multivariate response)
//______________________________________________________________________________

FULLCOND_variance_nonp_vector_nigmix::FULLCOND_variance_nonp_vector_nigmix(MCMCoptions * o,
                         vector<FULLCOND_const*> & p,DISTRIBUTION * d,
                         const ST::string & ti, const ST::string & fp,const ST::string & fr, 
                         const vector<unsigned long> & indicators, const vector<double> & v0, const vector<double> & v1,
                         const vector<double> & t2s, const vector<double> & at2, const vector<double> & bt2,
                         const vector<double> & omegas, const vector<double> & aomega, const vector<double> & bomega,
                         const vector<bool> & omegaf,
//                         const datamatrix start_data, 
                         const vector<bool> & omegaad,
                         const vector<unsigned> & ct, const unsigned & c)
                         : FULLCOND(o,datamatrix(1,1),ti,1,1,fp)
  {

  fctype = MCMC::variance;

  update_sigma2 = true;

  column = c;

  pathresults = fr;

  Cp = p;

  distrp = d;

  cut = ct;

  priorassumptions.push_back("\\\\");

  // Initialize the matrix for the variances
  datamatrix helpvariances;
  helpvariances = datamatrix(cut[cut.size()-1],1,0);
  for(unsigned i=0; i<cut.size()-1; i++)
    {
    helpvariances.putRowBlock(cut[i],cut[i+1],Cp[i]->get_variances());
    }
  setbeta(helpvariances);
  
  // Set the starting values 
//  startdata = start_data;      // Kontrolle der uebergebenen Startwerte
  nigmixsum = 0;
  omega_fix = omegaf[0];
  omega_adaptive = omegaad[0];
  v_0 = v0;
  v_1 = v1;
  a_t2 = at2;
  b_t2 = bt2;
  a_omega = aomega;
  b_omega = bomega;
  
    
  // Set pathes for the fullcond objects
  ST::string path1 = pathresults.substr(0,pathresults.length()-4)+"_shrinkage.raw";
  ST::string path2 = pathresults.substr(0,pathresults.length()-4)+"_indicator.raw";
  ST::string path3 = pathresults.substr(0,pathresults.length()-4)+"_t2.raw";

  // Fullcondobject for shrinkageparameter omega
  ST::string helptitle = Cp[0]->get_title();
  ST::string title = helptitle.substr(0,helptitle.length()-1);
  fc_shrinkage = FULLCOND(o,datamatrix(nrpar,1),title+"_shrinkage",nrpar,1,path1);
  fc_shrinkage.setflags(MCMC::norelchange);

  // Fullcondobject for Variancekomponents indicator and t2
  fc_indicator = FULLCOND(o,datamatrix(nrpar,1),title+"_indicator",nrpar,1,path2);
  fc_indicator.setflags(MCMC::norelchange);
  fc_t2 = FULLCOND(o,datamatrix(nrpar,1),title+"_t2",nrpar,1,path3);
  fc_t2.setflags(MCMC::norelchange);
  
  // Set the variablenames of the fullcond objects
  unsigned int i, k  = 0;
  vector<ST::string> varnames1(nrpar);
  vector<ST::string> varnames2(nrpar);
  vector<ST::string> helpvarnames;

  for(unsigned j=0; j<cut.size()-1; j++)
    {
    helpvarnames = Cp[j]->get_datanames(); 
    i = 0;   
    for(k=cut[j]; k<cut[j+1]; k++, i++)
      { 
      varnames1[k] = helpvarnames[i];
      if(omega_adaptive==true)
        {
        varnames2[k] = helpvarnames[i];
        }
      if(omega_adaptive==false)
        {
        varnames2[k] = "w";
        }
      }
    }
    
    fc_t2.init_names(varnames1);
    fc_indicator.init_names(varnames1); 
    fc_shrinkage.init_names(varnames2); 

    // Get current value of shrinkageparameter omega and variancecomponents
    double * workomega = fc_shrinkage.getbetapointer();
    double * workindicator = fc_indicator.getbetapointer();
    double * workt2 = fc_t2.getbetapointer();
    

    for(unsigned int i=0; i<nrpar; i++, workindicator++, workt2++, workomega++)
      {
      *workindicator = indicators[i];
      *workt2 = t2s[i];
      *workomega = omegas[i];
      }
//TEMP:BEGIN--------------------------------------------------------------------
// Kontrolle der uebergebenen Werte
  // ofstream outstartvectors_vnv("c:/bayesx/test/startvectors_nigmix_vnv.txt", ios::trunc|ios::app);
  // outstartvectors_vnv << "datamatrix with starting values"<< "\n";
  // outstartvectors_vnv << "effect I t2 w v0 v1 a b aw bw wfix \n";    
  // outstartvectors_vnv << "\n";
  // outstartvectors_vnv << "vektoren"<< "\n" << "varnam lambda tau2 I t2 w v0 v1 a b aw bw wfix " << "\n";
  // for(i=0; i<indicators.size(); i++)
    // {
    // outstartvectors_vnv << varnames1[i]  << "  " 
                        // << 1/helpvariances(i,0) << "  " 
                        // << helpvariances(i,0) << "  " 
                        // << indicators[i] << "  " 
                        // << t2s[i] << "  " 
                        // << omegas[i] << "  " 
                        // << v0[i] << "  "
                        // << v1[i] << "  " 
                        // << at2[i] << "  " 
                        // << bt2[i] << "  " 
                        // << aomega[i] << "  " 
                        // << bomega[i] << "  " 
                        // << omegaf[i] << "  " 
                        // << "\n";
    // }

//TEMP:END----------------------------------------------------------------------  
    }
//______________________________________________________________________________
//
// COPY CONSTRUCTOR
//______________________________________________________________________________

FULLCOND_variance_nonp_vector_nigmix::FULLCOND_variance_nonp_vector_nigmix(const FULLCOND_variance_nonp_vector_nigmix & t)
    : FULLCOND(FULLCOND(t))
  {
  update_sigma2 = t.update_sigma2;
  column = t.column;
  pathresults = t.pathresults;
  Cp = t.Cp;
  distrp = t.distrp;
  fc_shrinkage = t.fc_shrinkage;
  v_0 = t.v_0;
  v_1 = t.v_1;
  a_t2 = t.a_t2;
  b_t2 = t.b_t2;
  a_omega = t.a_omega;
  b_omega = t.b_omega;
  omega_fix = t.omega_fix;
  omega_adaptive = t.omega_adaptive;
  fc_indicator = t.fc_indicator;
  fc_t2 = t.fc_t2;
  nigmixsum = t.nigmixsum;
  startdata = t.startdata;
  cut = t.cut;

  }


//______________________________________________________________________________
//
// OVERLOADED ASSIGNMENT OPERATOR
//______________________________________________________________________________

const FULLCOND_variance_nonp_vector_nigmix & FULLCOND_variance_nonp_vector_nigmix::operator=(
                                        const FULLCOND_variance_nonp_vector_nigmix & t)
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
  v_0 = t.v_0;
  v_1 = t.v_1;
  a_t2 = t.a_t2;
  b_t2 = t.b_t2;
  a_omega = t.a_omega;
  b_omega = t.b_omega;
  omega_fix = t.omega_fix;    
  omega_adaptive = t.omega_adaptive;
  fc_indicator = t.fc_indicator;
  fc_t2 = t.fc_t2;
  nigmixsum = t.nigmixsum;
  startdata = t.startdata;
  cut = t.cut;

  return *this;
  }


  // Pointer auf das Shrinkagearameter lambda Fullcond-Objekt
FULLCOND * FULLCOND_variance_nonp_vector_nigmix::get_shrinkagepointer()
  {
  return &fc_shrinkage;
  }

//______________________________________________________________________________
//
// FUNCTION: update
// TASK: - stores sampled parameters in file 'samplepath'
//         storing order: first row, second row, ...
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::update(void)
  {
  acceptance++;
  unsigned i, j, k;
  
  double probv1;
  double rand_bernoulli;
  double helpomega;

  // Variables for summs
  nigmixsum = 0;
  unsigned int sumindicatorv0 = 0;
  unsigned int sumindicatorv1 = 0;

  // Get current value of (first)
  int iteration = optionsp->get_nriter();                  // Iteration
  double * workindicator = fc_indicator.getbetapointer();  // indicator
  double * workt2 = fc_t2.getbetapointer();                // t2
  double * workvariance;                                   // variance
  double * workomega = fc_shrinkage.getbetapointer();      // shrinkagearameter omega
  double * workbeta;                                       // regressioncoefficients
  double sigma2 = distrp->get_scale(column);                // scale parameter


  // Datei-Ausgabe der Startwerte
  //------------------------------
  if(iteration == 1)
    {
    get_startvalues();
    }


    
//TEMP:BEGIN--------------------------------------------------------------------
//ofstream output_beta("c:/bayesx/test/nigmix_beta.txt", ios::out|ios::app);
//ofstream output_scale("c:/bayesx/test/nigmix_scale.txt", ios::out|ios::app);
//ofstream output_var("c:/bayesx/test/nigmix_var.txt", ios::out|ios::app);
//ofstream output_w("c:/bayesx/test/w_nigmix.txt", ios::out|ios::app);
//ofstream outputn("c:/bayesx/test/test_nigmix_vnv.txt", ios::out|ios::app);
//outputn << "Termnr  Iteration  Beta Indicator" << "\n";
//outputn << iteration << "\n";
// Sacling of the Kovariates:
//output_beta << distrp->get_trmult(column)<< "\n" ;
//output_var << "\n" ;
//output_scale << sigma2 <<   "\n" ; 
//TEMP:END----------------------------------------------------------------------



  // Gibbs-Update of varianceparameters t2 with Inverse Gammadistribution and
  // indicator with Binomialdistribution
  //---------------------------------------------------------------------------
  workt2 = fc_t2.getbetapointer(); 
  workindicator = fc_indicator.getbetapointer();
  workomega = fc_shrinkage.getbetapointer();
  
  k = 0;
  for(j=0; j<cut.size()-1; j++)
    {
    workbeta = Cp[j]->getbetapointer();
    for(k=cut[j]; k<cut[j+1]; k++, workbeta++, workindicator++, workt2++, workomega++)
      {
      probv1 = 1/(1+((1-(*workomega))/(*workomega)*sqrt(v_1[0]/v_0[0])*exp(-(1/v_0[0]-1/v_1[0])*(*workbeta)*(*workbeta)/(2*sigma2 * (*workt2)))));
      rand_bernoulli = bernoulli(probv1);

      if(rand_bernoulli==0)
        {
        *workindicator = 0.0;
        sumindicatorv0 = sumindicatorv0 + 1;
        *workt2 = rand_invgamma(0.5+a_t2[0],b_t2[0]+(*workbeta)*(*workbeta)/(2* sigma2 * v_0[0]));
        beta(k,0) = v_0[0] * *workt2;
       }
      if(rand_bernoulli==1)
        {
        *workindicator = 1.0;
        sumindicatorv1 = sumindicatorv1 + 1;
        *workt2 = rand_invgamma(0.5+a_t2[0],b_t2[0]+(*workbeta)*(*workbeta)/(2* sigma2 * v_1[0]));
        beta(k,0) = v_1[0] * *workt2;
        }

      nigmixsum = nigmixsum + ((*workbeta)*(*workbeta))/beta(k,0);  // sum(beta^2/tau^2)

      }
    }
  fc_indicator.update();
  fc_t2.update();
  fc_indicator.acceptance++;
  fc_t2.acceptance++;

  // Transfer of the variance updates
  //---------------------------------
  for(j=0; j<cut.size()-1; j++)
    {
    workvariance = Cp[j]->getvariancespointer();
    for(k=cut[j]; k<cut[j+1]; k++, workvariance++)
      {
      *workvariance = beta(k,0);
      }
    }

  // Transfer of nigmixsum to the distribution objekt
  //-------------------------------------------------
  distrp->update_nigmix(nigmixsum);
  
  
  // Gibbs-Update of Shrinkageparameter omega with Betadistribution
  //----------------------------------------------------------------
  if(omega_fix==false && omega_adaptive==false)
    { 
    workomega = fc_shrinkage.getbetapointer();
    helpomega = rand_beta(a_omega[0] + sumindicatorv1, b_omega[0] + sumindicatorv0);
    for(i=0;i<nrpar;i++,workomega++)
      {
      *workomega = helpomega;
      }
    }
  if(omega_fix==false && omega_adaptive==true)
    { 
    workomega = fc_shrinkage.getbetapointer();
    for(i=0;i<nrpar;i++,workomega++)
      {
      *workomega = rand_beta(a_omega[i] + 1, b_omega[i] + 1);
      }
    }

  fc_shrinkage.update();
  fc_shrinkage.acceptance++;


  FULLCOND::update();

  }



//______________________________________________________________________________
//
// FUNCTION: outresults
// TASK: - write results for varianceparameters=idicator*t2 
//         to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::outresults(void)
  {
  FULLCOND::outresults();
  vector<ST::string> vnames = fc_indicator.get_datanames();

  unsigned int i;
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
  optionsp->out("  Results for the variances are stored in file\n");
  optionsp->out("  " + pathresults + "\n");
  optionsp->out("\n");

  // Ausgabe der restlichen Ergebnisse
  outresults_indicator();
  outresults_t2();
  outresults_shrinkage();
  
  }


//______________________________________________________________________________
//
// FUNCTION: outresults_indicator
// TASK: - write results for indicator to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::outresults_indicator(void)
  {

  fc_indicator.outresults();
  ST::string indicator_pathresults = pathresults.substr(0,pathresults.length()-7) + "indicator.res";
  vector<ST::string> vnames = fc_indicator.get_datanames();
  unsigned int i;
  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  // Dateiausgabe
  ofstream outp(indicator_pathresults.strtochar());

  // Kopfzeile Dateiausgabe
  if (indicator_pathresults.isvalidfile() != 1)
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

    m = fc_indicator.get_betamean(i,0);

    if (fc_indicator.get_betavar(i,0) <= 0)
      stddouble = 0;
    else
      stddouble = sqrt(fc_indicator.get_betavar(i,0));
     
    // Dateiausgabe
    if (indicator_pathresults.isvalidfile() != 1)
      {
      outp << (i+1) << "   ";
      outp << vnames[i] << "   ";
      outp << m << "   ";
      outp << stddouble << "   ";

      outp << fc_indicator.get_beta_lower1(i,0) << "   ";
      outp << fc_indicator.get_beta_lower2(i,0) << "   ";
      outp << fc_indicator.get_betaqu50(i,0) << "   ";
      outp << fc_indicator.get_beta_upper2(i,0) << "   ";
      outp << fc_indicator.get_beta_upper1(i,0) << "   ";
      if (fc_indicator.get_beta_lower1(i,0) > 0)
        outp << "1   ";
      else if (fc_indicator.get_beta_upper1(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";
      if (fc_indicator.get_beta_lower2(i,0) > 0)
        outp << "1   ";
      else if (fc_indicator.get_beta_upper2(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";
      outp << endl;

      // Bildschirmausgabe
      optionsp->out(ST::outresults(nsp,vnames[i],m,stddouble,
                    fc_indicator.get_beta_lower1(i,0),
                    fc_indicator.get_betaqu50(i,0),
                    fc_indicator.get_beta_upper1(i,0)) + "\n");

      }

    }

    
  optionsp->out("\n");
  optionsp->out("  Results for the variance component I are also stored in file\n");
  optionsp->out("  " + indicator_pathresults + "\n");
  optionsp->out("\n");
  }


//______________________________________________________________________________
//
// FUNCTION: outresults_shrinkage
// TASK: - write results for t2 to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::outresults_t2(void)
  {
  fc_t2.outresults();
  ST::string t2_pathresults = pathresults.substr(0,pathresults.length()-7) + "t2.res";
  vector<ST::string> vnames = fc_t2.get_datanames();
  unsigned int i;
  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  // Dateiausgabe
  ofstream outp(t2_pathresults.strtochar());

  // Kopfzeile Dateiausgabe
  if (t2_pathresults.isvalidfile() != 1)
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

    m = fc_t2.get_betamean(i,0);

    if (fc_t2.get_betavar(i,0) == 0)
      stddouble = 0;
    else
      stddouble = sqrt(fc_t2.get_betavar(i,0));
     
    // Dateiausgabe
    if (t2_pathresults.isvalidfile() != 1)
      {
      outp << (i+1) << "   ";
      outp << vnames[i] << "   ";
      outp << m << "   ";
      outp << stddouble << "   ";

      outp << fc_t2.get_beta_lower1(i,0) << "   ";
      outp << fc_t2.get_beta_lower2(i,0) << "   ";
      outp << fc_t2.get_betaqu50(i,0) << "   ";
      outp << fc_t2.get_beta_upper2(i,0) << "   ";
      outp << fc_t2.get_beta_upper1(i,0) << "   ";
      if (fc_t2.get_beta_lower1(i,0) > 0)
        outp << "1   ";
      else if (fc_t2.get_beta_upper1(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";
      if (fc_t2.get_beta_lower2(i,0) > 0)
        outp << "1   ";
      else if (fc_t2.get_beta_upper2(i,0) < 0)
        outp << "-1   ";
      else
        outp << "0   ";
      outp << endl;

      // Bildschirmausgabe
      optionsp->out(ST::outresults(nsp,vnames[i],m,stddouble,
                    fc_t2.get_beta_lower1(i,0),
                    fc_t2.get_betaqu50(i,0),
                    fc_t2.get_beta_upper1(i,0)) + "\n");

      }

    }
 
  
  optionsp->out("\n");
  optionsp->out("  Results for the variance component t2 are also stored in file\n");
  optionsp->out("  " + t2_pathresults + "\n");
  optionsp->out("\n");
  }


//______________________________________________________________________________
//
// FUNCTION: outresults_shrinkage
// TASK: - write results for shrinkaageparameter omega to output window and files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::outresults_shrinkage(void)
  {
  if(omega_fix==false)
    {
    unsigned int nr = nrpar;
    if(omega_adaptive==false)
      {
      nr = 1;
      }
   
    fc_shrinkage.outresults();
    ST::string shrinkage_pathresults = pathresults.substr(0,pathresults.length()-7) + "shrinkage.res";
    vector<ST::string> vnames = fc_shrinkage.get_datanames();
    unsigned int i;
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
    optionsp->out("  Results for the shrinkage parameter omega are also stored in file\n");
    optionsp->out("  " + shrinkage_pathresults + "\n");
    optionsp->out("\n");
    }
  }


//______________________________________________________________________________
//
// FUNCTION: outoptions
// TASK: - write options to output window
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::outoptions(void)
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
  optionsp->out("  Priors: NIGMIX shrinkage priors\n");
  optionsp->out("\n");
  
  if(omega_adaptive == false)
    {
    optionsp->out("  Hyperparameter v0 for variance component I: " +
                     ST::doubletostring(v_0[0]) + "\n" );
    optionsp->out("  Hyperparameter v1 for variance component I: " +
                     ST::doubletostring(v_1[0]) + "\n" );
    optionsp->out("\n");
    optionsp->out("  Hyperparameter a for variance component t2: " +
                     ST::doubletostring(a_t2[0]) + "\n" );
    optionsp->out("  Hyperparameter b for variance component t2: " +
                     ST::doubletostring(b_t2[0]) + "\n" );
    optionsp->out("\n");
    if(omega_fix==false)
      {
      optionsp->out("  Hyperparameter aw for w: " +
                       ST::doubletostring(a_omega[0]) + "\n" );
      optionsp->out("  Hyperparameter bw for w: " +
                       ST::doubletostring(b_omega[0]) + "\n" );
      }
    if(omega_fix==true)
      {
      optionsp->out("  Parameter w is fixed at value: " +
                       ST::doubletostring(fc_shrinkage.getbeta(0,0)) + "\n" );
      }
    optionsp->out("\n");
    }
 
  if(omega_adaptive == true)
    {
    for(unsigned i=0; i<nrpar; i++)
      {
      optionsp->out("  Hyperparameter v0 for variance component I of " +  varnames[i] + ": " +
                       ST::doubletostring(v_0[i]) + "\n" );
      optionsp->out("  Hyperparameter v1 for variance component I of " +  varnames[i] + ": " +
                       ST::doubletostring(v_1[i]) + "\n" );
      optionsp->out("  Hyperparameter a for variance component t2 of " +  varnames[i] + ": " +
                       ST::doubletostring(a_t2[i]) + "\n" );
      optionsp->out("  Hyperparameter b for variance component t2 of " +  varnames[i] + ": " +
                       ST::doubletostring(b_t2[i]) + "\n" );
      if(omega_fix==false)
        {
        optionsp->out("  Hyperparameter aw for w of " +  varnames[i] + ": " +
                         ST::doubletostring(a_omega[i]) + "\n" );
        optionsp->out("  Hyperparameter bw for w of " +  varnames[i] + ": " +
                         ST::doubletostring(b_omega[i]) + "\n" );
        }
      if(omega_fix==true)
        {
        optionsp->out("  Parameter w of " +  varnames[i] + " is fixed at value: " +
                         ST::doubletostring(fc_shrinkage.getbeta(i,0)) + "\n" );
        }
      optionsp->out("\n");
      }
    }
    
  optionsp->out("\n");

  }

//______________________________________________________________________________
//
// FUNCTION: get_samples
// TASK: - write samples to files
//______________________________________________________________________________

void FULLCOND_variance_nonp_vector_nigmix::get_samples(const ST::string & filename,
                                                       const unsigned & step) const
  {
  FULLCOND::get_samples(filename, step);

  ST::string pathhelp1 = pathresults.substr(0,pathresults.length()-7)+"indicator_sample.raw";
  ST::string pathhelp2 = pathresults.substr(0,pathresults.length()-7)+"t2_sample.raw";
  ST::string pathhelp3 = pathresults.substr(0,pathresults.length()-7)+"shrinkage_sample.raw";

  optionsp->out(pathhelp1 + "\n");
  optionsp->out("\n");
  fc_indicator.get_samples(pathhelp1);

  optionsp->out(pathhelp2 + "\n");
  optionsp->out("\n");
  fc_t2.get_samples(pathhelp2);

  optionsp->out(pathhelp3 + "\n");
  optionsp->out("\n");
  fc_shrinkage.get_samples(pathhelp3);
  }

//______________________________________________________________________________
//
// FUNCTION: get_startvalues
// TASK: - write startvalues to files
//______________________________________________________________________________

  
void FULLCOND_variance_nonp_vector_nigmix::get_startvalues(void)
  {
  
    unsigned i, j, k;
   
    // Hilfsvariablen
    double * helpvariance;
    double * helpindicator = fc_indicator.getbetapointer();
    double * helpt2 = fc_t2.getbetapointer();
    double * helpshrinkage = fc_shrinkage.getbetapointer();
    vector<ST::string> helpvarname_indicator = fc_indicator.get_datanames();
    vector<ST::string> helpvarname_t2 = fc_t2.get_datanames();
    vector<ST::string> helpvarname_variance;

    // Varianzkomponenten indikator und t2 und Varianzen
    ST::string indicator_pathstartdata = pathresults.substr(0,pathresults.length()-7) + "indicator_startdata.raw";
    ofstream indicator_outoutstartdata(indicator_pathstartdata.strtochar());
    indicator_outoutstartdata << "varname startvalue" << endl;

    ST::string t2_pathstartdata = pathresults.substr(0,pathresults.length()-7) + "t2_startdata.raw";
    ofstream t2_outoutstartdata(t2_pathstartdata.strtochar());
    t2_outoutstartdata << "varname startvalue" << endl;

    ST::string variance_pathstartdata = pathresults.substr(0,pathresults.length()-7) + "variance_startdata.raw";
    ofstream variance_outoutstartdata(variance_pathstartdata.strtochar());
    variance_outoutstartdata << "varname startvalue" << endl;
    
    // Shrinkageparameter omega
    ST::string shrinkage_pathstartdata = pathresults.substr(0,pathresults.length()-7) + "shrinkage_startdata.raw";
    ofstream shrinkage_outoutstartdata(shrinkage_pathstartdata.strtochar());
    shrinkage_outoutstartdata << "varname startvalue" << endl;
    
    // Hyperparameter
    ST::string hyperpar_pathstartdata = pathresults.substr(0,pathresults.length()-7) + "hyperpar_startdata.raw";
    ofstream hyperpar_outoutstartdata(hyperpar_pathstartdata.strtochar());
    hyperpar_outoutstartdata << "varname v0 v1 a b aw bw wfix adaptive" << endl;

    for(j=0; j<cut.size()-1; j++)
    {
    helpvariance = Cp[j]->getvariancespointer();
    helpvarname_variance = Cp[j]->get_datanames();
    i = 0;
    for(k=cut[j]; k<cut[j+1]; k++, i++, helpvariance++, helpindicator++, helpt2++,helpshrinkage++)
      {
      indicator_outoutstartdata << helpvarname_indicator[k] << " " << *helpindicator << endl;
      t2_outoutstartdata << helpvarname_t2[k] << " " << *helpt2 << endl;
      variance_outoutstartdata << helpvarname_variance[i] << " " << *helpvariance << endl;
      shrinkage_outoutstartdata << helpvarname_t2[k] << " " << *helpshrinkage << "  " << endl;
      hyperpar_outoutstartdata << helpvarname_t2[k]  << " " << v_0[0] << " " << v_1[0] << " " 
                               << a_t2[0] << " " << b_t2[0] << " " 
                               << a_omega[0] << " " << b_omega[0] << " " 
                               << omega_fix  << " " << omega_adaptive << endl;
      }
    }

  }
} // end: namespace MCMC






