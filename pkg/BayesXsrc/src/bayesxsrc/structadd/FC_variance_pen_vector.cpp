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



#include"FC_variance_pen_vector.h"



//------------------------------------------------------------------------------
//-- CLASS: FC_variance_pen_vector implementation of member functions ----------
//------------------------------------------------------------------------------


namespace MCMC
{

using randnumbers::rand_inv_gaussian;
using randnumbers::rand_gamma;


void FC_variance_pen_vector::add_variable(datamatrix & x,vector<ST::string> & op,
                                          vector<ST::string> & vn)
  {

  int f;
  double ta;
  double shrink;
  bool sfix;
  bool adsh;
  double a,b;
  double shw;

  // read options
  f = op[33].strtodouble(ta);
  f = op[29].strtodouble(shrink);
  if (op[30]=="true")
    sfix = true;
  else
    sfix = false;
  if (op[32]=="true")
    adsh = true;
  else
    adsh = false;
  f = op[5].strtodouble(a);
  f = op[6].strtodouble(b);
  f = op[31].strtodouble(shw);


  // options from each term
  tau2.push_back(ta);
  shrinkageweight.push_back(shw);

  // options from first term or all terms. Depends on adaptiveshrinkage
  shrinkagefix.push_back(sfix);
  adaptiveshrinkage.push_back(adsh);

  // shrinkagefix and adaptiveshrinkage are set from first term
  is_fix = shrinkagefix[0];
  is_adaptive = adaptiveshrinkage[0];

  if(is_ridge==1)                        // for Ridge
    shrinkagestart.push_back(1/(2*ta));
  if(is_ridge==0)                        // for Lasso
    shrinkagestart.push_back(shrink);

  a_shrinkagegamma.push_back(a);
  b_shrinkagegamma.push_back(b);

  if(is_adaptive==false)
    {
    shrinkagestart[nrpen] = shrinkagestart[0];
    a_shrinkagegamma[nrpen] = a_shrinkagegamma[0];
    b_shrinkagegamma[nrpen] = b_shrinkagegamma[0];
    }

  // optionsp->out("\n");
  // optionsp->out("OPTIONS FROM TERMS\n");
  // optionsp->out("penalty = " + ST::doubletostring(shelp(0,0)) + "\n");
  // optionsp->out("nrpen = " + ST::doubletostring(nrpen) + "\n");
  // optionsp->out("tau2 = " + ST::doubletostring(tau2[nrpen]) + "\n");
  // optionsp->out("weight = " + ST::doubletostring(shrinkageweight[nrpen]) + "\n");
  // optionsp->out("shrinkage = " + ST::doubletostring(shrinkagestart[nrpen]) + "\n");
  // optionsp->out("a_shrinkage = " + ST::doubletostring(a_shrinkagegamma[nrpen]) + "\n");
  // optionsp->out("b_shrinkage = " + ST::doubletostring(b_shrinkagegamma[nrpen]) + "\n");
  // optionsp->out("\n");

  nrpen++;
  Cp->tau2 = datamatrix(nrpen,1,0);
  Cp->tau2oldinv = datamatrix(nrpen,1,0);
  for(int i=0;i<nrpen;i++)
    {
    Cp->tau2(i,0) = tau2[i];
    }
  }


//______________________________________________________________________________
//
// CONSTRUCTOR with Parameters
//
// mp                : pointer to MASTER_object
// o                 : pointer to GENERAL_options object
// p                 : pointer to FULLCOND_lin_pen object
// d                 : pointer to DISTRIBUTION object
// ti                : reference to title of the full conditional(for example "fixed effects")
// fp                : reference to file path for storing sampled parameters
// isr               : variable,indicates if L2- or L1- penalty is uesd
// fr                : reference to filename for storing results
// shrinkage_start   : Starting value for the shrinkageparameter
// a_shrinkage_gamma : Hyperparameter for gammaprior of shrinkageparameter
// b_shrinkagegamma  : Hyperparameter for gammaprior of shrinkageparameter
// shrinkage_fix     : Should the shrinkageparameter be fixed at the value "shrinkage_start"
// ct                : blocks of regression coefficients
// c                 : reference to responsecategory (important in the case of multivariate response)
//______________________________________________________________________________
//==============================================================================
// KOMMENTAR Susanne:
//-------------------
// Derzeit ist die scale-dependent Variante implementiert, d.h. die
// Shrinkageprioris hängen in Gaussfall vom Skalenparameter sigma2 ab.
// Dies hat zur Folge, das zum update von sigma2 der Wert sum(beta^2/tau^2),
// der in lassosum und ridgesum berechnet wird, benötigt wird.
// Für den Fall das kein adaptives shrinkage gewählt wird, werden die
// ersten Einträge aus den Vektoroptionen für das Verfahren verwendet.
// Einentsprechender Hinweis sollte in das Handbuch.
//==============================================================================

FC_variance_pen_vector::FC_variance_pen_vector(MASTER_OBJ * mp,
                        GENERAL_OPTIONS * o, FC_linear_pen * p,
                        DISTR * d,const ST::string & ti,
                        const ST::string & fp, bool isr)
                        :FC(o,ti,1,1,fp)
  {

  is_ridge = isr;

  update_sigma2 = true;

  Cp = p;
  shelp = datamatrix(3,1,0);
  if (is_ridge)
    shelp(0,0) = 2;
  else
    shelp(0,0) = 1;

  distrp = d;

  priorassumptions.push_back("\\\\");

  // vector<ST::string> vnames = FC_shrinkage.get_datanames();
  //Initialisieren der Betamatrizen für die Varianzan + Übergabe der Startwerte
  //is_fix = shrinkagefix[0];
  //is_adaptive = adaptiveshrinkage[0];

  nrpen = 0;
  pensum = 0;

  }
//______________________________________________________________________________
//
// COPY CONSTRUCTOR
//______________________________________________________________________________

FC_variance_pen_vector::FC_variance_pen_vector(const FC_variance_pen_vector & t)
    : FC(FC(t))
  {
  tau2 = t.tau2;
  shelp = t.shelp;
  update_sigma2 = t.update_sigma2;
  Cp = t.Cp;
  distrp = t.distrp;
  FC_shrinkage = t.FC_shrinkage;
  shrinkagefix = t.shrinkagefix;
  adaptiveshrinkage = t.adaptiveshrinkage;
  a_shrinkagegamma = t.a_shrinkagegamma;
  b_shrinkagegamma = t.b_shrinkagegamma;
  shrinkagestart = t.shrinkagestart;
  shrinkageweight = t.shrinkageweight;
  nrpen = t.nrpen;
  pensum = t.pensum;
  is_ridge = t.is_ridge;
  is_fix = t.is_fix;
  is_adaptive = t.is_adaptive;
  }


//______________________________________________________________________________
//
// OVERLOADED ASSIGNMENT OPERATOR
//______________________________________________________________________________

const FC_variance_pen_vector & FC_variance_pen_vector::operator=(
                                        const FC_variance_pen_vector & t)
  {
  if (this == &t)
    return *this;
  FC::operator=(FC(t));
  tau2 = t.tau2;
  shelp = t.shelp;
  update_sigma2 = t.update_sigma2;
  Cp = t.Cp;
  distrp = t.distrp;
  FC_shrinkage = t.FC_shrinkage;
  shrinkagefix = t.shrinkagefix;
  adaptiveshrinkage = t.adaptiveshrinkage;
  a_shrinkagegamma = t.a_shrinkagegamma;
  b_shrinkagegamma = t.b_shrinkagegamma;
  shrinkagestart = t.shrinkagestart;
  shrinkageweight = t.shrinkageweight;
  nrpen = t.nrpen;
  pensum = t.pensum;
  is_ridge = t.is_ridge;
  is_fix = t.is_fix;
  is_adaptive = t.is_adaptive;
  return *this;
  }


bool FC_variance_pen_vector::posteriormode(void)
  {

  shelp(1,0) = double(tau2.size());

  datamatrix tau2inv(tau2.size(),1,0);
  unsigned i;
  for (i=0;i<tau2inv.rows();i++)
    tau2inv(i,0) = 1/tau2[i];

  setbeta(tau2inv);

  return true;
  }


//______________________________________________________________________________
//
// FUNCTION: update
// TASK: - stores sampled parameters in file 'samplepath'
//         storing order: first row, second row, ...
//______________________________________________________________________________

void FC_variance_pen_vector::update(void)
  {
  acceptance++;
  int i;
  double helpshrinkage;


  // reset Variables for summs
//  int nrpen = beta.rows();
  double rand_invgaussian = 0;
  double sumvariances = 0;
  double sumregcoeff = 0;

  pensum = 0;

  if (optionsp->nriter == 1)
    {

    FC_shrinkage = FC(optionsp,"",shrinkagefix.size(),1,samplepath + ".shrinkage");
//    Cp->tau2 = datamatrix(nrpen,1,0);
  // get current value of first shrinkagearameter
    double * shrinkagep = FC_shrinkage.beta.getV();
    for(i=0;i<nrpen;i++,shrinkagep++)
      {
      *shrinkagep = shrinkagestart[i];
      }
    }

  // get current value of first shrinkagearameter
  double * shrinkagep = FC_shrinkage.beta.getV();

  // get current value of first regressionparameter
  double * workbeta = Cp->beta.getV();

  // getcurrent value of sqrt(scale) parameter
  double sigma = sqrt(distrp->get_scale());

//TEMP:BEGIN--------------------------------------------------------------------
  // ofstream outpenreg("c:/bayesx/test/outpenreg.txt", ios::out|ios::app);
  // outpenreg << "beta ";
  // for(unsigned i=0; i<nrpen; i++,workbeta++) {outpenreg << *workbeta  << " ";}
  // outpenreg << "\n";
  // outpenreg << "tau2 ";
  // for(unsigned i=0; i<nrpen; i++) {outpenreg << beta(i,0)  << " ";}
  // outpenreg << "\n";
  // outpenreg << "lambda ";
  // for(unsigned i=0; i<nrpen; i++,shrinkagep++) {outpenreg << *shrinkagep  << " ";}
  // outpenreg << "\n";
  // outpenreg << "sigma "<< sigma <<"\n";

  // shrinkagep = FC_shrinkage.beta.getV();
  // workbeta = Cp->beta.getV();
  // sigma = sqrt(distrp->get_scale());
//TEMP:END----------------------------------------------------------------------



  // Gibbs-Update of varianceparameters tau^2
  //-----------------------------------------
  if (is_ridge == 0 && is_adaptive == false)   // Lasso L1-Penalty
   {
   for(i=0; i<nrpen; i++, workbeta++, shrinkagep++)
    {
    if (*workbeta>0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(sigma*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
      beta(i,0) = 1.0/rand_invgaussian;
      }
    if (*workbeta<0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(-1.0*sigma*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
      beta(i,0) = 1.0/(rand_invgaussian);
      }
    if (*workbeta==0 || *shrinkagep<=0)
      {
      beta(i,0) = 1E-6;
      }
    pensum += ((*workbeta)*(*workbeta))/beta(i,0);                          // sum(beta^2/tau^2)
    sumvariances = sumvariances + beta(i,0);  // sum(tau^2/weights^2) of current variances
   }
  }

  if (is_ridge == 0 && is_adaptive == true)   // Lasso L1-Penalty adaptive
   {
   for(i=0; i<nrpen; i++, workbeta++, shrinkagep++)
    {
    if (*workbeta>0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(sigma*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
      beta(i,0) = 1.0/rand_invgaussian;
      }
    if (*workbeta<0 && *shrinkagep>0)
      {
      rand_invgaussian = rand_inv_gaussian(-1.0*sigma*(*shrinkagep)/(*workbeta), (*shrinkagep)*(*shrinkagep));
      beta(i,0) = 1.0/(rand_invgaussian);
      }
    if (*workbeta==0 || *shrinkagep<=0)
      {
      beta(i,0) = 1E-6;
      }
    pensum +=  ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)
    //sumvariances = sumvariances + beta(i,0);                    // sum(tau^2) of current variances
   }
  }

  if (is_ridge == 1 && is_adaptive == false)   // Ridge L2-penalty
    {
    for(i=0; i<nrpen; i++, workbeta++, shrinkagep++)
      {
      beta(i,0) = 1/(2*(*shrinkagep));
      pensum +=  ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)
      sumregcoeff = sumregcoeff + (*workbeta)*(*workbeta);
      }
    }

  if (is_ridge == 1 && is_adaptive == true)   // Ridge L2-penalty adaptive
    {
    for(i=0; i<nrpen; i++, workbeta++, shrinkagep++)
      {
      beta(i,0) = 1/(2*(*shrinkagep));
      pensum += ((*workbeta)*(*workbeta))/beta(i,0);  // sum(beta^2/tau^2)
      //sumregcoeff = sumregcoeff + (*workbeta)*(*workbeta);
      }
    }

  // Gibbs-Update of Shrinkageparameter with Gammadistribution
  //----------------------------------------------------------
  if(is_fix==false)
    {
    shrinkagep = FC_shrinkage.beta.getV();
    workbeta = Cp->beta.getV();
    if(is_ridge == 0 && is_adaptive == false)            // Lasso L1-penalty
      {
      helpshrinkage = sqrt(rand_gamma(nrpen + a_shrinkagegamma[0], b_shrinkagegamma[0] + 0.5*sumvariances));
      for(i=0; i<nrpen; i++, shrinkagep++)
        {
        *shrinkagep = helpshrinkage;
        }
      }
    if(is_ridge == 0 && is_adaptive == true)            // Lasso L1-penalty adaptive
      {
      for(i=0; i<nrpen; i++, shrinkagep++)
        {
        *shrinkagep = sqrt(rand_gamma(1 + a_shrinkagegamma[i], b_shrinkagegamma[i] + 0.5*beta(i,0)));
        }
      }


    if(is_ridge == 1 && is_adaptive == false)            // Ridge L2-penalty
      {
      helpshrinkage = rand_gamma(0.5*nrpen + a_shrinkagegamma[0], b_shrinkagegamma[0] + sumregcoeff/(sigma*sigma));
      for(i=0; i<nrpen; i++, shrinkagep++)
        {
        *shrinkagep = helpshrinkage;
        }
      }
    if(is_ridge == 1 && is_adaptive == true)            // Ridge L2-penalty adaptive
      {
      for(i=0; i<nrpen; i++, shrinkagep++, workbeta++)
        {
        *shrinkagep = rand_gamma(0.5*1 + a_shrinkagegamma[i], b_shrinkagegamma[i] + ((*workbeta)*(*workbeta))/(sigma*sigma));
        }
      }
    }


  shelp(2,0) = pensum;
  distrp->update_scale_hyperparameters(shelp);

  Cp->tau2.assign(beta);

  // Update Shrinkageparameter
  FC_shrinkage.acceptance++;
  FC_shrinkage.update();

//  acceptance++;
  FC::update();


  }


//______________________________________________________________________________
//
// FUNCTION: get_samples
// TASK: - write samples to files
//______________________________________________________________________________

void FC_variance_pen_vector::get_samples(const ST::string & filename,ofstream & outg) const
  {

  FC::get_samples(filename,outg);

  ST::string pathhelp = filename.substr(0,filename.length()-14)+"shrinkage_sample.raw";

  FC_shrinkage.get_samples(pathhelp,outg);

  }

//______________________________________________________________________________
//
// FUNCTION: outresults
// TASK: - write results for varianceparameters to output window and files
//______________________________________________________________________________

void FC_variance_pen_vector::outresults(ofstream & out_stata, ofstream & out_R,
                             const ST::string & pathresults)
  {

//  int i;
  vector<ST::string> vnames;
  vnames = Cp->datanames;

  FC::outresults(out_stata,out_R,"");
  FC::outresults_help(out_stata,out_R,pathresults,vnames);

  optionsp->out("\n");
  optionsp->out("    Results for variances are also stored in file\n");
  optionsp->out("    " + pathresults + "\n");
  optionsp->out("\n");


  // outresults shrinkage
  ST::string shrinkage_pathresults = pathresults.substr(0,pathresults.length()-7) + "shrinkage.res";
  outresults_shrinkage(shrinkage_pathresults);

  }

//______________________________________________________________________________
//
// FUNCTION: outresults_shrinkage
// TASK: - write results for shrinkageparameter to output window and files
//______________________________________________________________________________

//void FC_variance_pen_vector::outresults_shrinkage(void)
void FC_variance_pen_vector::outresults_shrinkage(const ST::string & pathresults)
  {

  if(is_fix==false)
    {
    optionsp->out("\n");
    optionsp->out("  MAIN_REGRESSOIN: linear effects (shrinkage) \n");
    optionsp->out("\n");
    ST::string shrinkage_pathresults = pathresults;
    vector<ST::string> vnames;
    vnames = Cp->datanames;

    ofstream out1;
    ofstream out2;

    FC_shrinkage.outresults(out1,out2,"");
    FC_shrinkage.outresults_help(out1,out2,shrinkage_pathresults,vnames);

    // unsigned nr;
    // if(is_adaptive==true)
      // {
      // FC_shrinkage.outresults_help(out1,out2,shrinkage_pathresults,vnames);
      // }

    // if(is_adaptive==false)
      // {
      // FC_shrinkage.outresults_help(out1,out2,shrinkage_pathresults,vnames);
      // }

    optionsp->out("\n");
    optionsp->out("    Results for shrinkage parameters are also stored in file\n");
    optionsp->out("    " + shrinkage_pathresults + "\n");
    optionsp->out("\n");
    }
  }


//______________________________________________________________________________
//
// FUNCTION: outoptions
// TASK: - write options to output window
//______________________________________________________________________________

void FC_variance_pen_vector::outoptions(void)
  {

//TEMP:BEGIN--------------------------------------------------------------------
// KOMMENTAR Susanne:  FC_shrinkage kann dann hier wieder rausgenommen weden
//   FC_shrinkage = FC(optionsp,"",shrinkagefix.size(),1,samplepath + ".shrinkage");
//TEMP:END----------------------------------------------------------------------

  int i;
//  int nrpen = beta.rows();
  vector<ST::string> vnames;   // = FC_shrinkage.datanames.getV();
  vnames = Cp->datanames;

  if(is_ridge == 0)
    {
    optionsp->out("  LINEAR EFFECTS WITH LASSO PENALTY: \n");
    }
  if(is_ridge == 1)
    {
    optionsp->out("  LINEAR EFFECTS WITH RIDGE PENALTY: \n");
    }

  if(is_adaptive==false)
    {
    optionsp->out("  Hyperparameter a for shrinkage: " +
                     ST::doubletostring(a_shrinkagegamma[0]) + "\n" );
    optionsp->out("  Hyperparameter b for shrinkage: " +
                     ST::doubletostring(b_shrinkagegamma[0]) + "\n" );
    if(is_fix==true)
      {
      optionsp->out("  Shrinkage is fixed at value = " +
                       ST::doubletostring(shrinkagestart[0]) + "\n" );
      }
    optionsp->out("\n");
    }

  if(is_adaptive==true)
    {
    for(i=0; i<nrpen; i++)
      {
      optionsp->out("  Hyperparameter a for shrinkage of " + vnames[i] + ": " +
                       ST::doubletostring(a_shrinkagegamma[i]) + "\n" );
      optionsp->out("  Hyperparameter b for shrinkage of " + vnames[i] + ": " +
                       ST::doubletostring(b_shrinkagegamma[i]) + "\n" );
      if(is_fix==true)
        {
        optionsp->out("  Shrinkage of " + vnames[i] + " is fixed at value = " +
                         ST::doubletostring(shrinkagestart[i]) + "\n" );
        }
      }
    optionsp->out("\n");

    }

  }


} // end: namespace MCMC





