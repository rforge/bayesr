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



#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#endif

#include"mcmcsimul2.h"
#include<time.h>
#include"clstring.h"
#include <stdlib.h>
#include<math.h>


namespace MCMC
{


STEPWISErun::STEPWISErun(MCMCoptions * go,DISTRIBUTION * dp,
vector<FULLCOND*> & fc)
: MCMCsimulate(go,dp,fc)
  {
  }


STEPWISErun::STEPWISErun(const STEPWISErun & st)
  : MCMCsimulate(MCMCsimulate(st))
  {
  }


const STEPWISErun & STEPWISErun::operator=(
const STEPWISErun & st)
  {
  if (this==&st)
    return *this;
  MCMCsimulate::operator=(MCMCsimulate(st));
  return *this;
  }


// -----------------------------------------------------------------------------
// ------- die grundlegenden Funktionen ----------------------------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::posteriormode(const vector<ST::string> & header,
                                 const bool & presim)
  {
  unsigned i,j;
  //unsigned nrmodels = genoptions_mult.size();
  bool errors=false;
  errors = checkerrors(likep_mult[0],fullcondp,begin[0],end[0]);
  if (errors)
    return true;

  if (header[0] !="")
    {
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out(header[0] + "\n", true,false,16);
    genoptions_mult[0]->out("\n");
    }

  if (!presim)
    {
    if (likepexisting)
      {
      genoptions_mult[0]->out("RESPONSE DISTRIBUTION:\n",true);
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  " + likep_mult[0]->get_family() + "\n");
      genoptions_mult[0]->out("  Number of observations: " +
             ST::inttostring(likep_mult[0]->get_response().rows()) + "\n");
      genoptions_mult[0]->out("\n");
      }

    if (likepexisting)
      {
      set_center(likep_mult[0],fullcondp,begin[0],end[0]);
      }
    } // end: if (!presim)

  unsigned it=1;
  bool converged=false;
  bool allconverged=false;

  bool convergedouterloop=false;

  unsigned k=0;
  while ( (!convergedouterloop) && (k< 100) )
    {
    k++;

    likep_mult[0]->compute_iwls();

    if(likep_mult[0]->iwlsweights_constant() == false)
      {
      for(i=0;i<fullcond_alle.size();i++)
        fullcond_alle[i]->set_calculate_xwx();
      }

    converged=false;
    it=1;
    while ((!converged) && (it <= 100))
      {
      allconverged = true;

      if(likepexisting)
        if(likep_mult[0]->posteriormode() == false)
           allconverged = false;

      if(!isboot)
        {
        for(j=begin[0];j<=end[0];j++)
          {
          if (fullcondp[j]->posteriormode() == false)
            allconverged = false;
          } // end: for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
        }
      else
        {
        double start = begin[0];
        if(start==0)
          {
          if (fullcondp[start]->posteriormode() == false)
            allconverged = false;
          start += 1;
          }
        for(j=start;j<=end[0];j++)
          {
          if(fullcondp[j]->get_lambda()!=0)
            {
            if (fullcondp[j]->posteriormode() == false)
              allconverged = false;
            }
          } // end: for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
        }

      if (allconverged && it>1)         // stellt sicher, dass nicht zufällig das "betaold" vom
        converged = true;               // Auswählen der Lambdas zu einer fälschlichen Konvergenz führt

      #if defined(BORLAND_OUTPUT_WINDOW)
      Application->ProcessMessages();
      if (Frame->stop)
        {
        break;
        }
      if (Frame->pause)
        {
        genoptions_mult[0]->out("\n");
        genoptions_mult[0]->out("SIMULATION PAUSED\n");
        genoptions_mult[0]->out("Click CONTINUE to proceed\n");
        genoptions_mult[0]->out("\n");
        while (Frame->pause)
          {
          Application->ProcessMessages();
          }
        genoptions_mult[0]->out("SIMULATION CONTINUED\n");
        genoptions_mult[0]->out("\n");
        }

      #elif defined(JAVA_OUTPUT_WINDOW)
      bool stop = genoptions_mult[0]->adminb_p->breakcommand();
      if(stop)
        break;
      #endif

      it++;

      } // end:   while ((!converged) && (it <= 100))

//if(it==101 && !converged)
//  genoptions_mult[0]->out("BACKFITTING ALGORITHM DID NOT CONVERGE\n",true,true,12,255,0,0);

    allconverged = true;

    if (likepexisting)
      if (likep_mult[0]->posteriormode_converged(k) == false)
        allconverged = false;

    if(!isboot)
      {
      for(j=begin[0];j<=end[0];j++)
        {
        if (fullcondp[j]->posteriormode_converged(k) == false)
          allconverged = false;
        } // end: for(j=0;j<nrfullcond;j++)
      }
    else
      {
      double start = begin[0];
      if(start==0)
        {
        if (fullcondp[start]->posteriormode_converged(k) == false)
          allconverged = false;
        start += 1;
        }
      for(j=start;j<=end[0];j++)
        {
        if(fullcondp[j]->get_lambda()!=0)
          {
          if (fullcondp[j]->posteriormode_converged(k) == false)
            allconverged = false;
          }
        } // end: for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
      }

    if (allconverged)
      convergedouterloop = true;

    if (likepexisting)
      likep_mult[0]->posteriormode_set_beta_mode();

    for(j=begin[0];j<=end[0];j++)
       fullcondp[j]->posteriormode_set_beta_mode();

    } // end: while ( (!convergedouterloop) && (k< 100) )

//if(k==100 && !convergedouterloop)
//  genoptions_mult[0]->out("LOCAL SCORING PROCEDURE DID NOT CONVERGE\n",true,true,12,255,0,0);

  if (!presim)
    {
    #if defined(BORLAND_OUTPUT_WINDOW)
    if (!Frame->stop)
    #elif defined(JAVA_OUTPUT_WINDOW)
    if (!genoptions_mult[0]->adminb_p->get_stop())
    #endif
      {
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("ESTIMATION RESULTS:\n",true);
      genoptions_mult[0]->out("\n");

      if (!converged)
        genoptions_mult[0]->out("BACKFITTING ALGORITHM DID NOT CONVERGE\n",true,true,12,255,0,0);
      if (!convergedouterloop)
        genoptions_mult[0]->out("LOCAL SCORING PROCEDURE DID NOT CONVERGE\n",true,true,12,255,0,0);
      genoptions_mult[0]->out("\n");

      if (likepexisting)
        {
        likep_mult[0]->outresults();
        }

      for(j=begin[0];j<=end[0];j++)
        fullcondp[j]->outresults();

//        return false;

      } // end: if Frame->stop

    #if defined(BORLAND_OUTPUT_WINDOW)
    else
      {
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out(
      "ESTIMATION TERMINATED BY USER BREAK\n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("Estimation results: none\n");
      genoptions_mult[0]->out("\n");

      if (likepexisting)
        likep_mult[0]->reset();

      for(j=begin[0];j<=end[0];j++)
          fullcondp[j]->reset();

      errors=true;
//        return true;
      }
    #elif defined(JAVA_OUTPUT_WINDOW)
    else
      {
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("Estimation results: none\n");
      genoptions_mult[0]->out("\n");

      if (likepexisting)
        likep_mult[0]->reset();

      for(j=begin[0];j<=end[0];j++)
        fullcondp[j]->reset();

      errors=true;
//        return true;
      }
    #endif
    } // end: if (!presim)

  return errors;
  } // end: function posteriormode


bool STEPWISErun::stepwise(const ST::string & procedure, const ST::string & minimum,
        const ST::string & crit, const int & stp, const ST::string & trac,
        const int & number, const ST::string & stam, const int & inc,
        const int & boot, const bool & uncond,
        const datamatrix & Da, const vector<ST::string> & modelvar,
        const ST::string & name, vector<FULLCOND*> & fullcond_z, ST::string & path,
        const ST::string & CI, bool & hier, const double & prop)
  {

  D = Da;
  modelv = modelvar;
  algorithm = procedure;
  minim = minimum;
  minim2 = minimum;
  criterion = crit;
  increment = inc;
  steps = stp;
  startmodel = stam;
  trace = trac;
  hierarchical = hier;
  bootstrap = boot;
  isboot = false;
  unconditional = uncond;       // falsch rum definiert!

  modell_alt.erase(modell_alt.begin(),modell_alt.end());
  text_alt = " ";
  vector<double> modell_final;
  double kriterium_final;
  ST::string text_final;

  ST::string tr_akt = "trace_on";
  posttitle.push_back("");

  lambdavec.erase(lambdavec.begin(),lambdavec.end());
  names_fixed.erase(names_fixed.begin(),names_fixed.end());
  names_nonp.erase(names_nonp.begin(),names_nonp.end());

  fullcond_alle = fullcondp;

  set_center(likep_mult[0],fullcond_alle,begin[0],end[0]);  // sorgt dafür, daß Funktionen zentriert werden!

  bool gewichte = false;
  if(likep_mult[0]->get_family() != "Gaussian")
    gewichte = true;
  initialise_lambdas(names_nonp,names_fixed,lambdavec,number,gewichte);

  if(criterion == "MSEP" || criterion == "AUC" ||criterion == "CV5" || criterion == "CV10")
    {
    initialise_weights(prop);
    for(unsigned y=0;y<fullcond_alle.size();y++)
      fullcond_alle[y]->set_calculate_xwx();
   }

  unsigned i;

  // Fehlerabfragen
  unsigned j;
  for(j=0;j<fullcond_alle.size();j++)
     {
     if(fullcond_alle[j]->geterrors().size()>0)
      {
      for(i=0;i<fullcond_alle[j]->geterrors().size();i++)
         genoptions_mult[0]->out(fullcond_alle[j]->geterrors()[i]);
      return true;
      }
    }
  // überprüfen, dass Randomslopes nicht auch als fixe Effekte angegeben werden!
  if( vcm_doppelt() == true)      // Für VCM-Modelle!!! --> muß die Funktion raus?
     return true;

  vector<vector<unsigned> > startindex;
  vector<vector<double> > startfix;
  startwerte(startmodel,startindex,startfix);

  options_text(number,startfix,startindex,name);

  ST::string path_tex = path + "_model_summary.tex";
  outtex.open(path_tex.strtochar());
  make_graphics(name,startindex);

  bool first = true;
  bool abbruch = false;

  ST::string pathmodels = path + "_models.raw";
  ST::string pathcrit = path + "_criterium.raw";
  ST::string pathtrace = path + "_criterium_trace.raw";
  outcriterium.open(pathcrit.strtochar());
  outcriterium << "step   var   " << criterion << endl;
  outmodels.open(pathmodels.strtochar());
  outmodels << "step   " << criterion << "   model" << endl << endl;
//  outtrace.open(pathtrace.strtochar());
//  outtrace << "step   " << criterion << endl << endl;

  for(i=0;i<startindex.size();i++)
     {
     abbruch = single_stepwise(startindex[i],startfix[i],true);

     if(abbruch==true)
        return true;

     if(minim=="adaptiv")
       {
       schaetzen(0,kriterium_alt,true,"backfitting");
       outcriterium << "B   0   " << ST::doubletostring(kriterium_alt,8) << endl;
       outmodels << endl << "B   " << ST::doubletostring(kriterium_alt,8) << endl;
//       outtrace << ST::inttostring(steps_aktuell+1) << "  " << ST::doubletostring(kriterium_alt,8) << endl;
       }

     if(first)
       {
       first = false;
       kriterium_final = kriterium_alt;
       modell_final = modell_alt;
       text_final = text_alt;
       }
     else
       {
       if(kriterium_final>kriterium_alt)
         {
         kriterium_final = kriterium_alt;
         modell_final = modell_alt;
         text_final = text_alt;
         }
       }
     }
  ST::string header = "  Final Model:";
  fix_komplett(modell_final);
  fullcond_komplett(modell_final);
  tr_akt = "trace_on";
  if(steps==0)
    maketext(header,modell_final,kriterium_final,text_final,true,tr_akt,true);
  else
    maketext(header,modell_final,kriterium_final,text_final,false,tr_akt,false);
  genoptions_mult[0]->out("\n\n");
  kriterium_tex = kriterium_final;

  if(abbruch==true)
    return true;

  if(criterion == "MSEP" || criterion == "AUC")
    {
    likep_mult[0]->weight_for_all();
    for(unsigned y=0;y<fullcond_alle.size();y++)
      fullcond_alle[y]->set_calculate_xwx();
    }

  // Files werden nicht mehr gebraucht und müssen wieder geschlossen werden!!!
  outcriterium.close();
  outmodels.close();

  if(CI == "none")
    {
    fullcond_z = fullcondp;
    for(i=0;i<fullcond_z.size();i++)
      fullcond_z[i]->set_fcnumber(i);
    posteriormode(posttitle,false); // Problem: linearer Prädiktor bei "true" standardisiert! Hier wird zurückgerechnet!
                                      // danach nicht mehr compute_criterion() aufrufen!!!
    }
  else
    {
    if(criterion == "MSEP" || criterion == "AUC")
      posteriormode(posttitle,true);
    abbruch = confidence_intervals(CI,modell_final,kriterium_final,fullcond_z);
    if(abbruch == true)
      return true;
    }

  make_tex_end(path,modell_final,CI);
  outtex.close();

// FÜR SIMULATIONEN:
/*
  // gibt Lambdas aus, damit man die richtig bestimmten Variablen zählen kann!
  ST::string zaehlername = path + "_lambdas_" + likep_mult[0]->get_responsename() + ".ascii";
  ofstream out(zaehlername.strtochar());
  ST::string beschriftung = "  krit   ";
  ST::string eintrag = "  " + ST::doubletostring(kriterium_final) + "   ";
  for(i=1;i<names_fixed.size();i++)
     beschriftung = beschriftung + names_fixed[i] + "   ";
  for(i=0;i<names_nonp.size();i++)
     beschriftung = beschriftung + names_nonp[i][0] + "   ";
  for(i=0;i<modell_final.size();i++)
     eintrag = eintrag + ST::doubletostring(modell_final[i]) + "   ";
  out << beschriftung << endl;
  out << eintrag << endl;
  out.close();

  // gibt df für nichtlineare Funktionen aus!
  zaehlername = path + "_df_" + likep_mult[0]->get_responsename() + ".ascii";
  ofstream out2(zaehlername.strtochar());
  beschriftung = " ";
  eintrag = " ";
  for(i=0;i<names_nonp.size();i++)
     beschriftung = beschriftung + names_nonp[i][0] + "   ";
  for(i=(names_fixed.size()-1);i<modell_final.size();i++)
     eintrag = eintrag + ST::doubletostring(fullcond_alle[i+1]->compute_df()) + "   ";
  out2 << beschriftung << endl;
  out2 << eintrag << endl;
  out2.close();
// ENDE: FÜR SIMULATIONEN
*/

  return false;
  }


bool STEPWISErun::single_stepwise(const vector<unsigned> & start,
                            const vector<double> & startfix, const bool & tex)
  {
  modell_neu.erase(modell_neu.begin(),modell_neu.end());
  modellematrix.erase(modellematrix.begin(),modellematrix.end());
  steps_aktuell = 0;
  ST::string tr_akt = "trace_on";
  vector<vector<double> > startiteration;

  unsigned i;
  for(i=0;i<names_fixed.size()-1;i++)
     modell_neu.push_back(startfix[i]);
  for(i=1;i<fullcond_alle.size();i++)
     {
     double lambda = lambdavec[i-1][start[i-1]];
     modell_neu.push_back(lambda);
     }

  modell_alt = modell_neu;
  startiteration.push_back(modell_alt);
  modellematrix.push_back(startiteration);
  fix_komplett(modell_alt);
  fullcond_komplett(modell_alt);
  if(likep_mult[0]->get_family() == "Gamma")
    {
    likep_mult[0]->reset();
    for(i=0;i<fullcond_alle.size();i++)
      fullcond_alle[i]->reset();
    fullcond_alle[0]->setbeta(fullcond_alle[0]->get_nrpar(),1,0);
    }

  if(hierarchical == true)
    {
    for(i=fullcond_alle.size()-1;i>=1;i--)   // Abfrage, ob Startmodell hierarchisch ist!
       {
       ST::string possible = "alles";
       fullcond_alle[i]->hierarchical(possible);
       bool falsch = true;

       if((modell_alt[names_fixed.size()-2+i] == -1 || modell_alt[names_fixed.size()-2+i] > 0)
                           && (possible == "raus" || possible == "vraus"))
         falsch = false;
       if(modell_alt[names_fixed.size()-2+i] > 0 && possible == "rfix")
         falsch = false;
       if((modell_alt[names_fixed.size()-2+i] == -1 || modell_alt[names_fixed.size()-2+i] == 0)
                           && (possible == "spline" || possible == "vspline"))
         falsch = false;
       if(modell_alt[names_fixed.size()-2+i] == 0
                           && (possible == "spfix" || possible == "vspfix" || possible == "vfix"))
         falsch = false;

       if(falsch == false)
         {
         genoptions_mult[0]->out("  NOTE: The startmodel is no hierarchical model! Choose another one.");
         return true;
         }
       }
    }

  schaetzen(0,kriterium_alt,true,"backfitting");

  if(tex==true)
     {
     kriterium_tex = kriterium_alt;
     make_predictor();
     }

  kriterium_neu = kriterium_alt;
  outcriterium << steps_aktuell << "   0   " << ST::doubletostring(kriterium_neu,8) << endl;
  outmodels << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << "   ";
//  outtrace << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << endl;
  ST::string header;
  fertig = false;    // überprüft, ob es noch nicht gerechnete Modelle gibt
  ST::string text_neu;

  bool abbruch = false;
  if(algorithm != "coorddescent")
    abbruch = stepfunctions();
  else
    {
    abbruch = koordabstieg();

    /*if(likep_mult[0]->get_family() == "Gamma" && minim=="adaptiv")
      {
      while(modell_neu != modellematrix[0][0])
        {
        fertig = false;
        modell_neu = modell_alt;
        startiteration.erase(startiteration.begin(),startiteration.end());
        startiteration.push_back(modell_alt);
        modellematrix.erase(modellematrix.begin(),modellematrix.end());
        modellematrix.push_back(startiteration);
        schaetzen(0,kriterium_alt,true,"backfitting");
        kriterium_neu = kriterium_alt;
        abbruch = koordabstieg();
        if(abbruch==true)
         return true;
        }
      }*/
    }

  if(abbruch==true)
    return true;

  header = "  Final Model:";
  tr_akt = "trace_on";
  maketext(header,modell_alt,kriterium_alt,text_alt,false,tr_akt,false);
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  Used number of iterations: " + ST::inttostring(steps_aktuell));
  genoptions_mult[0]->out("\n\n");
  genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
  genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");

  return false;
  }


void STEPWISErun::schaetzen(int z, double & kriterium, bool neu, ST::string variante)
  {
  if(variante == "backfitting")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      posteriormode(posttitle,true);
      kriterium = compute_criterion();
      }
    else
      {
      likep_mult[0]->save_weightiwls();
      kriterium = 0;
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        for(unsigned y=0;y<fullcond_alle.size();y++)
          fullcond_alle[y]->set_calculate_xwx();
        posteriormode(posttitle,true);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        for(unsigned y=0;y<fullcond_alle.size();y++)
          fullcond_alle[y]->set_calculate_xwx();
        }
      }
    }

  else if(variante == "nonp" || variante == "fixnonp")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      fullcond_alle[z]->posteriormode();
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[z]->set_calculate_xwx();
        fullcond_alle[z]->posteriormode();
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        fullcond_alle[z]->set_calculate_xwx();
        }
      }
    }

  else if(variante == "nonpnonp")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
//      if(neu)
//        fullcond_alle[z]->const_varcoeff();
      fullcond_alle[z]->posteriormode();
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[0]->set_calculate_xwx();
        fullcond_alle[z]->set_calculate_xwx();
//        if(neu)
//          fullcond_alle[z]->const_varcoeff();
        fullcond_alle[z]->posteriormode();
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        fullcond_alle[0]->set_calculate_xwx();
        fullcond_alle[z]->set_calculate_xwx();
        }
      }
    }

  else if(variante == "factor" || variante == "nonpfix" || variante == "fixfix")
    {
    if(variante == "nonpfix")
      neu = false;
    else if(variante == "fixfix")
      neu = true;
    if(criterion != "CV5" && criterion != "CV10")
      {
      if(!neu)
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),true);
      else
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),false);
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      if(!neu)
        fullcond_alle[0]->include_effect(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects());
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[0]->set_calculate_xwx();
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                               fullcond_alle[z]->get_data_forfixedeffects(),false);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        fullcond_alle[0]->set_calculate_xwx();
        }
      }
    }

  else if(variante == "fix")
    {
    int col = column_for_fix(names_fixed[z]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[z]);

    if(criterion != "CV5" && criterion != "CV10")
      {
      if(neu)
        fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(col)),false);
      else
        fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(col)),true);
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      if(!neu)
        fullcond_alle[0]->include_effect(name_help,datamatrix(D.getCol(col)));
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[0]->set_calculate_xwx();
        fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(col)),false);
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        fullcond_alle[0]->set_calculate_xwx();
        }
      }
    }

  else if(variante == "leer")
    {
    if(criterion != "CV5" && criterion != "CV10")
      {
      fullcond_alle[0]->posteriormode_const();
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[0]->set_calculate_xwx();
        fullcond_alle[0]->posteriormode_const();
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        fullcond_alle[0]->set_calculate_xwx();
        }
      }
    }

  else if(variante == "nonpleer")
    {
    ST::string possible = "alles";
    if(hierarchical == true)
      fullcond_alle[z]->hierarchical(possible);

    if(criterion != "CV5" && criterion != "CV10")
      {
      if((possible != "valles" && possible != "vrfix" && possible != "vraus")
          || fullcond_alle[z]->is_identifiable() == true)
        fullcond_alle[0]->posteriormode_const();
      else // if(possible == "valles")  // bei Rauslassen von VC muß zugehöriger fixer Effekt upgedatet werden!
        {
        vector<ST::string> help;
        help.push_back(fullcond_alle[z]->get_datanames()[1]);
        fullcond_alle[0]->posteriormode_single(help,
                               fullcond_alle[z]->get_data_forfixedeffects(),false);
        }
      kriterium = compute_criterion();
      }
    else
      {
      kriterium = 0;
      likep_mult[0]->save_weightiwls();
      unsigned pcv;
      if(criterion == "CV5")
        pcv = 5;
      else
        pcv = 10;
      for(unsigned c=0;c<pcv;c++)
        {
        likep_mult[0]->compute_cvweights(c);
        fullcond_alle[0]->set_calculate_xwx();
        if((possible != "valles" && possible != "vrfix" && possible != "vraus")
            || fullcond_alle[z]->is_identifiable() == true)
          fullcond_alle[0]->posteriormode_const();
        else // if(possible == "valles")  // bei Rauslassen von VC muß zugehöriger fixer Effekt upgedatet werden!
          {
          vector<ST::string> help;
          help.push_back(fullcond_alle[z]->get_datanames()[1]);
          fullcond_alle[0]->posteriormode_single(help,
                               fullcond_alle[z]->get_data_forfixedeffects(),false);
          }
        kriterium += compute_criterion();
        likep_mult[0]->compute_cvweights(-1);
        fullcond_alle[0]->set_calculate_xwx();
        }
      }
    }

  }


// -----------------------------------------------------------------------------
// ------------------- Funktionen, für Stepwise / Stepmin ----------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::stepfunctions(void)
  {
  ST::string tr_akt = "trace_on";
  ST::string text_neu;
  bool eins = true;
       // Schleife für Minimierung
  while(kriterium_neu<=kriterium_alt && fertig==false && steps_aktuell<steps)
       {
       steps_aktuell = steps_aktuell + 1;
       vector<ST::string> textiteration;
       fertig = true;
       kriterium_alt = kriterium_neu;
       modell_alt = modell_neu;
       ST::string header = "  Startmodel:";
       bool neutext;
       if(eins == true)
         {
         eins = false;
         neutext = true;  // Modell rausschreiben !!
         }
       else
         {
         if(trace == "trace_off")
           tr_akt = "trace_off";
         neutext = false;
         }
       maketext(header,modell_alt,kriterium_alt,text_neu,neutext,
                tr_akt,neutext);
       text_alt = text_neu;

       vector<vector<double> > modeliteration;
       double kriterium;
       vector<double> kriteriumiteration2;
       unsigned j;

       if(algorithm == "stepwise")
         {
         unsigned z = stepwise_fixfactor(kriteriumiteration2,modeliteration,textiteration);
         stepwise_nonp(kriteriumiteration2,modeliteration,textiteration,z);
         }
       else // if(algorithm == "stepmin")
         {
         if(minim == "exact")
           {
           unsigned z = stepwise_fixfactor(kriteriumiteration2,modeliteration,textiteration);
           minexact_nonp(kriteriumiteration2,modeliteration,textiteration,z);
           }
         else
           {
           korrektur();  //fullcondp[0]->posteriormode_const();
           posteriormode(posttitle,true);
           step_minfix(kriteriumiteration2,modeliteration,textiteration);
           unsigned z = step_minfactor(kriteriumiteration2,modeliteration,textiteration);
           stepmin_nonp(kriteriumiteration2,modeliteration,textiteration,z);
           }
         }

       if(fertig==false)
         {
         kriterium = kriteriumiteration2[0];
         for(j=0;j<kriteriumiteration2.size();j++)  //berechnet den besten Wert
            {
            if(kriteriumiteration2[j]<=kriterium)
              {
              kriterium = kriteriumiteration2[j];
              modell_neu = modeliteration[j];
              text_neu = textiteration[j];
              }
            }
         kriterium_neu = kriterium;
         if(isboot==false)
           {
           outcriterium << steps_aktuell << "   0   " << ST::doubletostring(kriterium_neu,8) << endl;
           outmodels << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << "   ";
//           outtrace << steps_aktuell << "   " << ST::doubletostring(kriterium_neu,8) << endl;
           }
         modellematrix.push_back(modeliteration);

         header = "\n\nBest Model of this iteration:";
         fix_komplett(modell_neu);
         fullcond_komplett(modell_neu);
         neutext = false;
         maketext(header,modell_neu,kriterium_neu,text_neu,neutext,
                  trace,true);                    // Modell rausschreiben!!

         if(steps_aktuell==steps)
           {
           if(kriterium_alt>kriterium_neu)
             {
             kriterium_alt = kriterium_neu;
             modell_alt = modell_neu;
             text_alt = text_neu;
             }
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("  Attention: Maximum number " + ST::inttostring(steps) + " of iterations was reached! \n");
           }
         }
       else
         {
         if(trace == "trace_on" || trace == "trace_minim")
           {
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("  There are no new models for this iteration! \n");
           }
         if(isboot==false)
           {
           outcriterium << ST::inttostring(steps_aktuell) << "   0   " << ST::doubletostring(kriterium_neu,8) << endl;
           outmodels << ST::inttostring(steps_aktuell) << "   " << ST::doubletostring(kriterium_neu,8) << "   ";
           }
         }
       if(trace == "trace_on" || trace == "trace_minim")
         {
         genoptions_mult[0]->out("\n\n");
         genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
         genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
         }

       if(make_pause() == true)
          return true;
       }

  if(minim == "apprexact")
    {
    minim = "exact";
    if(fertig==false)
      {
      kriterium_neu = kriterium_alt;
      modell_neu = modell_alt;
      fix_komplett(modell_alt);
      fullcond_komplett(modell_alt);
      }
    //if(trace == "trace_on" || trace == "trace_minim")
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  Beginning of the exact minimization! \n");
    fertig = false;
    stepfunctions();
    }

  return false;
  }


unsigned STEPWISErun::stepwise_fixfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration)
  {
  unsigned i;
  for(i=1;i<names_fixed.size();i++)
    {
    modell_neu = modell_alt;
    if(modell_alt[i-1]==-1)
      modell_neu[i-1]= 0;
    else if(modell_alt[i-1]==0)
      modell_neu[i-1] = -1;
    if(modelcomparison(modell_neu,modellematrix)==false)
      newmodel_fix(modell_neu[i-1],kriteriumiteration2,modeliteration,
                            textiteration,names_fixed[i]);
    }

  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
     {
     modell_neu = modell_alt;
     ST::string possible = "alles";
     // VCM
     if(hierarchical)
       fullcond_alle[z]->hierarchical(possible);

     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false
                                               && possible == "alles")
       modell_neu[z+names_fixed.size()-2]= 0;
     else if(modell_alt[z+names_fixed.size()-2]==0)
       {
       modell_neu[z+names_fixed.size()-2] = -1;
       // VCM
       if(possible == "vfix")
         {
         for(i=0;i<names_nonp[z-1].size();i++)
           reset_fix(names_nonp[z-1][i]);
         kriterium_alt = MAXDOUBLE;
         }
       }
     if(modelcomparison(modell_neu,modellematrix)==false)
       newmodel_factor(modell_neu[z+names_fixed.size()-2],z,kriteriumiteration2,
                 modeliteration,textiteration,names_nonp[z-1]);
     z = z + 1;
     }
  return z;
  }

void STEPWISErun::stepwise_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  for(i=z;i<fullcond_alle.size();i++)
    {

    ST::string possible = "alles";
    if(hierarchical == true)
      fullcond_alle[i]->hierarchical(possible);

    unsigned sch;
    for(sch=1;sch<=unsigned(increment);sch++)
       {
       modell_neu = modell_alt;
       bool lambda_exist;    // zum Überprüfen, ob neues Lambda im Vektor enthalten ist
       unsigned index = search_lambdaindex(modell_alt[names_fixed.size()-2+i],
                                lambdavec[i-1],lambda_exist);
       lambda_exist = false;
       if(index < lambdavec[i-1].size()-sch)
          lambda_exist = true;
       if(lambda_exist==true && hierarchical == true)
         {
         // für 2d-Spline
         if(lambdavec[i-1][index+sch] == 0 && (possible == "spline" || possible == "vspline"
                                               || possible == "spfix" || possible == "vspfix"))
           lambda_exist = false;
         if(lambdavec[i-1][index+sch] == -1 && (possible == "spline" || possible == "raus"))
           lambda_exist = false;
         if(lambdavec[i-1][index+sch] == -2 && (possible == "vspline" || possible == "vraus"))
           lambda_exist = false;
         if(lambdavec[i-1][index+sch] > 0 && (possible == "rfix" || possible == "raus"
                                              || possible == "vrfix" || possible == "vraus"))
           lambda_exist = false;

         // für VCM
         if(lambdavec[i-1][index+sch] == 0 && possible == "vfix")
           lambda_exist = false;
         if(lambdavec[i-1][index+sch] == -1 && possible == "vfix")
           {   // bedeutet, daß fixer Effekt vorher nicht im Modell war, aber durch VC aufgenommen wurde.
           for(unsigned j=0;j<names_nonp[z-1].size();j++)
            reset_fix(names_nonp[z-1][j]);
           }
         }
       if(lambda_exist==true)
          {
          modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][index+sch];
          if(modelcomparison(modell_neu,modellematrix)==false)
            newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
          }

       lambda_exist = false;
       modell_neu = modell_alt;
       if(index >= sch)
          lambda_exist = true;
       if(lambda_exist==true && hierarchical == true)
         {
         // für 2d-Spline
         if(lambdavec[i-1][index-sch] == 0 && (possible == "spline" || possible == "vspline"
                                               || possible == "spfix" || possible == "vspfix"))
           lambda_exist = false;
         if(lambdavec[i-1][index-sch] == -1 && (possible == "spline" || possible == "raus"))
           lambda_exist = false;
         if(lambdavec[i-1][index-sch] == -2 && (possible == "vspline" || possible == "vraus"))
           lambda_exist = false;
         if(lambdavec[i-1][index-sch] > 0 && (possible == "rfix" || possible == "raus"
                                              || possible == "vrfix" || possible == "vraus"))
           lambda_exist = false;

         // für VCM
         if(lambdavec[i-1][index-sch] == 0 && possible == "vfix")
           lambda_exist = false;
         if(lambdavec[i-1][index-sch] == -1 && possible == "vfix")
           {   // bedeutet, daß fixer Effekt vorher nicht im Modell war, aber durch VC aufgenommen wurde.
           for(unsigned j=0;j<names_nonp[z-1].size();j++)
             reset_fix(names_nonp[z-1][j]);
           }
         }
       if(lambda_exist==true)
          {
          modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][index-sch];
          if(modelcomparison(modell_neu,modellematrix)==false)
            newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
          }
       }
    }
  }


void STEPWISErun::stepmin_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  for(i=z;i<fullcond_alle.size();i++)
    {
    modell_neu = modell_alt;
    unsigned lambda_ind;

    unsigned y;
    for(y=1;y<fullcond_alle.size();y++)
      fullcond_alle[y]->set_center(false);    // sorgt dafür, daß Funktionen nicht zentriert werden!

    vector<double> krit_fkt;
    if(modell_alt[i+names_fixed.size()-2]==0)
      stepmin_nonp_leer(i,krit_fkt,kriterium_alt);
    else if(modell_alt[i+names_fixed.size()-2]==-1)
      stepmin_nonp_fix(i,krit_fkt,kriterium_alt);
    else
      stepmin_nonp_nonp(i,krit_fkt,kriterium_alt);

    double kriterium_min = krit_fkt[0];
    unsigned j;
    lambda_ind = 0;
    for(j=1;j<krit_fkt.size();j++)  //berechnet den besten Wert
       {
       if(krit_fkt[j]<=kriterium_min)
         {
         kriterium_min = krit_fkt[j];
         lambda_ind = j;
         }
       }

    for(y=1;y<fullcond_alle.size();y++)
     {
     if(fullcond_alle[y]->is_identifiable() == false)
       fullcond_alle[y]->set_center(true);    // sorgt dafür, daß Funktionen zentriert werden!
     }

    modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][lambda_ind];
    if(modell_neu[names_fixed.size()-2+i] != modell_alt[names_fixed.size()-2+i])
      {
      if(modelcomparison(modell_neu,modellematrix)==false)
        {
        newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
        // Stellt linearen Prädiktor wieder her. Besser wäre, den lin. Prädiktor zu speichern!!!
        korrektur(); //fullcondp[0]->posteriormode_const();
        posteriormode(posttitle,true);
        }
      }
    }
  }

void STEPWISErun::minexact_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  unsigned j;
  for(i=z;i<fullcond_alle.size();i++)
    {
    modell_neu = modell_alt;
    unsigned lambda_ind;
    vector<double> krit_fkt;
    if(modell_alt[i+names_fixed.size()-2]==0)
      minexact_nonp_leer(i,krit_fkt,kriterium_alt);
    else if(modell_alt[i+names_fixed.size()-2]==-1)
      {
      reset_fix(names_nonp[i-1][0]);
      minexact_nonp_fix(i,krit_fkt,kriterium_alt);
      }
    else
      minexact_nonp_nonp(i,krit_fkt,kriterium_alt);

    double kriterium_min = krit_fkt[0];
    lambda_ind = 0;
    for(j=1;j<krit_fkt.size();j++)  //berechnet den besten Wert
       {
       if(krit_fkt[j]<=kriterium_min)
         {
         kriterium_min = krit_fkt[j];
         lambda_ind = j;
         }
       }

    modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][lambda_ind];
    if(modell_neu[names_fixed.size()-2+i] != modell_alt[names_fixed.size()-2+i])
      {
      if(modelcomparison(modell_neu,modellematrix)==false)
        {
        newmodel_nonp(i,kriteriumiteration2,modeliteration,textiteration);
        // Stellt linearen Prädiktor wieder her.
        korrektur();  // fullcondp[0]->posteriormode_const();
        posteriormode(posttitle,true);
        }
      }
    }
  }


// -----------------------------------------------------------------------------
// ------------------ Funktionen für Stepmin -----------------------------------
// -----------------------------------------------------------------------------

void STEPWISErun::step_minfix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration)
  {
  unsigned i;
  for(i=1;i<names_fixed.size();i++)
    {
    if(modell_alt[i-1]==-1)
      stepmin_fix_leer(kriteriumiteration2,modeliteration,textiteration,i);
    else if(modell_alt[i-1]==0)
      stepmin_leer_fix(kriteriumiteration2,modeliteration,textiteration,i);
    }
  }

void STEPWISErun::stepmin_fix_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration, unsigned & i)
  {
  fullcond_alle[0]->safe_const();
  reset_fix(names_fixed[i]);
  schaetzen(0,kriterium_neu,true,"leer");

  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    include_fix(names_fixed[i]);
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    reset_fix(names_fixed[i]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_alt,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
    {
    modell_neu[i-1] = 0;
    if(modelcomparison(modell_neu,modellematrix)==false)
      {
      newmodel(kriteriumiteration2,modeliteration,textiteration);
      include_fix(names_fixed[i]);
      korrektur();  // fullcondp[0]->posteriormode_const();
      posteriormode(posttitle,true);
      }
    else
      {
      int c = column_for_fix(names_fixed[i]);
      vector<ST::string> name_help;
      name_help.push_back(names_fixed[i]);
      fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
      }
    modell_neu[i-1] = -1;
    }
  else
    {
    int c = column_for_fix(names_fixed[i]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[i]);
    fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
    }
  }

void STEPWISErun::stepmin_leer_fix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration, unsigned & i)
  {
  fullcond_alle[0]->safe_const();
  schaetzen(i,kriterium_neu,false,"fix");
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    reset_fix(names_fixed[i]);
    int c = column_for_fix(names_fixed[i]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[i]);
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_alt,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
     {
     modell_neu[i-1] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       reset_fix(names_fixed[i]);
       korrektur();  // fullcondp[0]->posteriormode_const();
       posteriormode(posttitle,true);
       }
     else
       reset_fix(names_fixed[i]);
     modell_neu[i-1] = 0;
     }
  else
     reset_fix(names_fixed[i]);
  }

unsigned STEPWISErun::step_minfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration)
  {
  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
     {
     ST::string possible = "alles";
     // VCM
     if(hierarchical)
       fullcond_alle[z]->hierarchical(possible);

     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false
                                               && possible == "alles")
       stepmin_factor_leer(kriteriumiteration2,modeliteration,textiteration,z);
     else if(modell_alt[z+names_fixed.size()-2]==0)
       {
       // VCM
       if(possible == "vfix")
         {
         for(unsigned i=0;i<names_nonp[z-1].size();i++)
           reset_fix(names_nonp[z-1][i]);
         kriterium_alt = MAXDOUBLE;
         }
       stepmin_leer_factor(kriteriumiteration2,modeliteration,textiteration,z);
       }
     z = z + 1;
     }
  return z;
  }

void STEPWISErun::stepmin_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  fullcond_alle[0]->safe_const();
  for(i=0;i<names_nonp[z-1].size();i++)
    reset_fix(names_nonp[z-1][i]);
  schaetzen(0,kriterium_neu,true,"leer");
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_alt,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
     {
     modell_neu[z+names_fixed.size()-2] = 0;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
       korrektur();  // fullcondp[0]->posteriormode_const();
       posteriormode(posttitle,true);
       }
     else
       fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects(),true);
     modell_neu[z+names_fixed.size()-2] = -1;
     }
  else
     fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                 fullcond_alle[z]->get_data_forfixedeffects(),true);
  }

void STEPWISErun::stepmin_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z)
  {
  unsigned i;
  fullcond_alle[0]->safe_const();
  schaetzen(z,kriterium_neu,false,"factor");
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                 fullcond_alle[z]->get_data_forfixedeffects(),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_alt,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(kriterium_neu < kriterium_alt)
     {
     modell_neu[z+names_fixed.size()-2] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       for(i=0;i<names_nonp[z-1].size();i++)
         reset_fix(names_nonp[z-1][i]);
       korrektur();  // fullcondp[0]->posteriormode_const();
       posteriormode(posttitle,true);
       }
     else
       {
       for(i=0;i<names_nonp[z-1].size();i++)
         reset_fix(names_nonp[z-1][i]);
       }
     modell_neu[z+names_fixed.size()-2] = 0;
     }
  else
     {
     for(i=0;i<names_nonp[z-1].size();i++)
       reset_fix(names_nonp[z-1][i]);
     }
  }


void STEPWISErun::stepmin_nonp_nonp(unsigned & z, vector<double> & krit_fkt,double & kriterium)
  {
  unsigned i;
  ST::string possible = "alles";
  //if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);

  vector<FULLCOND*> fullcond_ori = fullcondp;
  fullcondp = fullcond_ori;

  if(minim == "adaptiv" || minim == "adap_exact" || criterion == "CV5" || criterion == "CV10")
    {
    fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
    schaetzen(z,kriterium,true,"nonpnonp");
    if(possible == "valles" || possible == "vrfix")
      fullcond_alle[0]->posteriormode_const();
    }

  if(hierarchical == false)
    possible == "alles";

  fullcondp[0]->safe_const();
  bool interact = false;
  fullcond_alle[z]->safe_splines(interact);

  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=modell_alt[z+names_fixed.size()-2])
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=-2 && lambdavec[z-1][i]!=-1 && lambdavec[z-1][i]!=0)
        {
        if(possible == "alles" || possible == "valles" || possible == "spline"
           || possible == "vspline" || possible == "spfix" || possible == "vspfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          schaetzen(z,kriterium_versuch,false,"nonpnonp");
          fullcond_alle[0]->set_const_old();
          }
        }
      else if(lambdavec[z-1][i]==-2)
        {
        if(possible == "alles" || possible == "valles" || possible == "spline" || possible == "vrfix"
              || possible == "spfix" || possible == "vspfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          schaetzen(z,kriterium_versuch,false,"nonpnonp");
          fullcond_alle[0]->set_const_old();
          }
        }
      else if(lambdavec[z-1][i]==-1)
        {
        if(possible == "alles" || possible == "valles" || possible == "spfix")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          fullcond_alle[z]->reset_effect(0);
          schaetzen(z,kriterium_versuch,false,"nonpfix");
          fullcond_alle[0]->set_const_old();
          reset_fix(names_nonp[z-1][0]);
          }
        }
      else
        {
        if(possible == "alles" || possible == "valles" || possible == "vrfix")
          {
          fullcond_alle[z]->set_inthemodel(0);
          fullcond_alle[z]->reset_effect(0);
          vector<FULLCOND*> fullcond_start;
          schaetzen(z,kriterium_versuch,false,"nonpleer");
          fullcond_alle[0]->set_const_old();
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      {
      krit_fkt.push_back(kriterium);
      if(interact == true && i+1<lambdavec[z-1].size())
        {
        if(modell_alt[z+names_fixed.size()-2]>0 && lambdavec[z-1][i+1]<=0)
          fullcond_alle[z]->set_splines_old();
        }
      }
    }

  fullcond_alle[z]->set_inthemodel(1);
  fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
  if(possible != "spline" && fullcond_alle[z]->is_identifiable() == false)   // bei Haupteffekten der ANOVA Zerlegung muß center=false bleiben
    fullcond_alle[z]->set_center(true);                                    // sonst wird Interaktion nicht geschätzt.
  fullcond_alle[z]->posteriormode();
  fullcond_alle[0]->update_linold();

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8)
                               + "   " + ST::doubletostring(krit_fkt[i],12) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
    vector<double> kriterium_control;
    genoptions_mult[0]->out("\n\n");
    unsigned y;
    for(y=1;y<fullcond_alle.size();y++)
      {
      if(fullcond_alle[y]->is_identifiable() == false)
        fullcond_alle[y]->set_center(true);
      }
    minexact_nonp_nonp(z,kriterium_control,kriterium);
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx    exact \n");
    for(i=0;i<kriterium_control.size();i++)
      genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],12) + "   "
                        + ST::doubletostring(kriterium_control[i],6) + "\n");
    }
  }


void STEPWISErun::stepmin_nonp_fix(unsigned & z, vector<double> & krit_fkt, double & kriterium)
  {
  unsigned i;

  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);
  if(possible == "valles")
    possible = "alles";

  vector<FULLCOND*> fullcond_ori = fullcondp;
  fullcondp = fullcond_ori;

  if(minim == "adaptiv" || minim == "adap_exact" || criterion == "CV5" || criterion == "CV10")
    {
    schaetzen(z,kriterium,true,"fixfix");
    }

  fullcond_alle[0]->safe_const();
  reset_fix(names_nonp[z-1][0]);
  fullcondp.push_back(fullcond_alle[z]);
  fullcond_alle[z]->set_inthemodel(1);

  unsigned m = 0;
  bool interact = false;
  fullcond_alle[z]->safe_splines(interact);  // für ANOVA-Zerlegung: speichert Haupteffekte

  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=-1)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=0)
        {
        if(possible == "alles" || possible == "spfix" || possible == "vfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          schaetzen(z,kriterium_versuch,false,"fixnonp");
          fullcond_alle[0]->set_const_old();
          }
        }
      else
        {
        if(possible == "alles" || possible == "rfix")
          {
          fullcond_alle[z]->set_inthemodel(0);
          fullcond_alle[z]->reset_effect(0);
          schaetzen(0,kriterium_versuch,true,"leer");
          fullcond_alle[0]->set_const_old();
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      {
      m = i;
      krit_fkt.push_back(kriterium);
      }
    }

  fullcond_alle[z]->set_inthemodel(-1);
  fullcond_alle[z]->reset_effect(0);
  fullcondp.erase(fullcondp.end()-1,fullcondp.end());
  fullcondp[0]->posteriormode_single(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects(),true);
  fullcond_alle[0]->update_linold();
  if(interact && possible=="alles")      // stellt den alten Zustand für ANOVA-Zerlegung wieder her
    {
    krit_fkt[m] = compute_criterion();
    fullcond_alle[z]->set_splines_old();
    fullcondp[0]->posteriormode_single(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects(),false);
    }

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8)
                              + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
    reset_fix(names_nonp[z-1][0]);
    vector<double> kriterium_control;

    unsigned y;
    for(y=1;y<fullcond_alle.size();y++)
      {
      if(fullcond_alle[y]->is_identifiable() == false)
        fullcond_alle[y]->set_center(true);
      }
    minexact_nonp_fix(z,kriterium_control,kriterium);
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx    exact \n");
    for(i=0;i<kriterium_control.size();i++)
      genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "   "
             + ST::doubletostring(kriterium_control[i],6) + "\n");
    }
  }


void STEPWISErun::stepmin_nonp_leer(unsigned & z, vector<double> & krit_fkt, double & kriterium)
  {
  unsigned i=0;

  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);
  if(possible == "valles")
    possible = "alles";

  // VCM
  if(possible == "vfix")
    {
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    kriterium = MAXDOUBLE;
    fullcond_alle[z]->set_inthemodel(0);
    modell_alt[names_fixed.size()-2+z] = -1;
    }

   vector<FULLCOND*> fullcond_ori = fullcondp;
   fullcondp = fullcond_ori;

  // if(minim == "adaptiv" || minim == "adap_exact")
           // ---> hier nicht nötig (siehe koordmin_leer_fix)
  if( ((criterion == "CV5" || criterion == "CV10") && possible != "vfix")
     || ((minim == "adaptiv" || minim == "adap_exact") && likep_mult[0]->get_family()=="Gamma") )
    {
    schaetzen(i,kriterium,true,"leer");
    }

  fullcondp = fullcond_ori;
  fullcond_alle[z]->const_varcoeff();   // Konstante anpassen für varrierenden Koeffizienten
  fullcond_alle[0]->safe_const();
  fullcondp.push_back(fullcond_alle[z]);
  fullcond_alle[z]->set_inthemodel(1);

  unsigned m = 0;
  bool interact = false;
  fullcond_alle[z]->safe_splines(interact);  // für ANOVA-Zerlegung: speichert Haupteffekte

  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=0)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=-1 && lambdavec[z-1][i]!=-2)
        {
        if(possible == "alles")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          schaetzen(z,kriterium_versuch,false,"nonp");
          fullcond_alle[0]->set_const_old();
          }
        }
      else if(lambdavec[z-1][i]==-2)
        {
        if(possible == "alles" || possible == "vrfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          schaetzen(z,kriterium_versuch,false,"nonp");
          fullcond_alle[0]->set_const_old();
          }
        }
      else
        {
        if(possible == "rfix" || possible == "alles" || possible == "vfix")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          fullcond_alle[z]->reset_effect(0);
          schaetzen(z,kriterium_versuch,false,"factor");
          reset_fix(names_nonp[z-1][0]);
          fullcond_alle[0]->set_const_old();
          }
        }
      krit_fkt.push_back(kriterium_versuch); // Länge des Vektors muß zu Anzahl Möglichkeiten passen!
      }
    else
      {
      m = i;
      krit_fkt.push_back(kriterium);
      }
    }

  fullcond_alle[z]->set_inthemodel(0);
  fullcond_alle[z]->reset_effect(0);
  fullcondp.erase(fullcondp.end()-1,fullcondp.end());
  fullcond_alle[0]->posteriormode_const();
  if(interact && possible=="alles")      // stellt den alten Zustand für ANOVA-Zerlegung wieder her
    {
    krit_fkt[m] = compute_criterion();
    fullcond_alle[z]->set_splines_old();
    fullcondp[0]->posteriormode_const();
    }

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8)
                              + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim == "approx_control")
    {
    vector<double> kriterium_control;
    unsigned y;
    for(y=1;y<fullcond_alle.size();y++)
      {
      if(fullcond_alle[y]->is_identifiable() == false)
        fullcond_alle[y]->set_center(true);
      }
    minexact_nonp_leer(z,kriterium_control,kriterium);
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalues: approx    exact \n");
    for(i=0;i<kriterium_control.size();i++)
      genoptions_mult[0]->out("        " + ST::doubletostring(krit_fkt[i],6) + "   "
             + ST::doubletostring(kriterium_control[i],6) + "\n");
    }
  }

//------------------------------------------------------------------------------------

void STEPWISErun::minexact_nonp_nonp(unsigned & z, vector<double> & krit_fkt,
                  double & kriterium)
  {
  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);
  if(possible == "valles")
    possible = "alles";

  unsigned i;
  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=modell_alt[z+names_fixed.size()-2])
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=-2 && lambdavec[z-1][i]!=-1 && lambdavec[z-1][i]!=0)
        {
        if(possible == "alles" || possible == "spline" || possible == "spfix"
           || possible == "vspline" || possible == "vspfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          }
        }
      else if(lambdavec[z-1][i]==-2)
        {
        if(possible == "alles" || possible == "spline" || possible == "spfix"
            || possible == "vrfix" || possible == "vspfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          }
        }
      else if(lambdavec[z-1][i]==-1)
        {
        if(possible == "alles" || possible == "spfix")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          vector<FULLCOND*> fullcond_start = fullcondp;
          vector<double> modell1 = modell_alt;
          modell1[z+names_fixed.size()-2] = -1;
          fullcond_einzeln(modell1,modell_alt,z);  // hier muß der Fullcond-Vekor angepaßt werden!!!
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          fullcondp = fullcond_start;
          end[0] = fullcondp.size()-1;
          reset_fix(names_nonp[z-1][0]);
          }
        }
      else
        {
        if(possible == "alles" || possible == "vrfix")
          {
          fullcond_alle[z]->set_inthemodel(0);
          vector<FULLCOND*> fullcond_start = fullcondp;
          vector<double> modell1 = modell_alt;
          modell1[z+names_fixed.size()-2] = 0;
          fullcond_einzeln(modell1,modell_alt,z);  // hier muß der Fullcond-Vekor angepaßt werden!!!
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          fullcondp = fullcond_start;
          end[0] = fullcondp.size()-1;
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->set_inthemodel(1);
  fullcond_alle[z]->update_stepwise(modell_alt[z+names_fixed.size()-2]);
  korrektur();  // fullcondp[0]->posteriormode_const();
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (exact): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8) + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }


void STEPWISErun::minexact_nonp_fix(unsigned & z, vector<double> & krit_fkt,
          double & kriterium)
  {
  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);
  if(possible == "valles")
    possible = "alles";

  unsigned i;
  vector<FULLCOND*> fullcond_begin = fullcondp;
  vector<double> modell1 = modell_alt;
  modell1[z+names_fixed.size()-2] = 1;
  fullcond_einzeln(modell1,modell_alt,z);  // hier muß der Fullcond-Vekor angepaßt werden!!!
  fullcond_alle[z]->set_inthemodel(1);
  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=-1)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=0)
        {
        if(possible == "alles" || possible == "spfix" || possible == "vfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          }
        }
      else
        {
        if(possible == "alles" || possible == "rfix")
          {
          fullcond_alle[z]->set_inthemodel(0);
          fullcondp = fullcond_begin;
          end[0] = fullcondp.size()-1;
          fullcond_alle[z]->reset_effect(0);
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->set_inthemodel(-1);
  fullcond_alle[z]->reset_effect(0);
  fullcond_alle[0]->include_effect(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects());
  fullcondp = fullcond_begin;
  end[0] = fullcondp.size()-1;
  korrektur();  // fullcondp[0]->posteriormode_const();
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (exact): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8) + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }


void STEPWISErun::minexact_nonp_leer(unsigned & z, vector<double> & krit_fkt,
                  double & kriterium)
  {
  ST::string possible = "alles";
  if(hierarchical == true)
    fullcond_alle[z]->hierarchical(possible);
  if(possible == "valles")
    possible = "alles";

  unsigned i;
  // VCM
  if(possible == "vfix")
    {
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    kriterium = MAXDOUBLE;
    }

  vector<FULLCOND*> fullcond_begin = fullcondp;
  vector<double> modell1 = modell_alt;
  modell1[z+names_fixed.size()-2] = 1;
  fullcond_einzeln(modell1,modell_alt,z);  // hier muß der Fullcond-Vekor angepaßt werden!!!
  fullcond_alle[z]->set_inthemodel(1);

  for(i=0;i<lambdavec[z-1].size();i++)
    {
    if(lambdavec[z-1][i]!=0)
      {
      double kriterium_versuch = MAXDOUBLE;
      if(lambdavec[z-1][i]!=-1 && lambdavec[z-1][i]!=-2)
        {
        if(possible == "alles")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          }
        }
      else if(lambdavec[z-1][i]==-2)
        {
        if(possible == "alles" || possible == "vrfix")
          {
          fullcond_alle[z]->update_stepwise(lambdavec[z-1][i]);
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          }
        }
      else
        {
        if(possible == "rfix" || possible == "alles" || possible == "vfix")
          {
          fullcond_alle[z]->set_inthemodel(-1);
          fullcondp = fullcond_begin;
          end[0] = fullcondp.size()-1;
          fullcond_alle[z]->reset_effect(0);
          fullcondp[0]->include_effect(names_nonp[z-1],
                                fullcond_alle[z]->get_data_forfixedeffects());
          korrektur();  // fullcondp[0]->posteriormode_const();
          schaetzen(z,kriterium_versuch,false,"backfitting");
          reset_fix(names_nonp[z-1][0]);
          }
        }
      krit_fkt.push_back(kriterium_versuch);
      }
    else
      krit_fkt.push_back(kriterium);
    }
  fullcond_alle[z]->set_inthemodel(0);
  fullcond_alle[z]->reset_effect(0);
  fullcondp = fullcond_begin;
  end[0] = fullcondp.size()-1;
  korrektur();  // fullcondp[0]->posteriormode_const();
  posteriormode(posttitle,true);

  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (exact): \n");
     for(i=0;i<krit_fkt.size();i++)
       genoptions_mult[0]->out(" " + ST::doubletostring(lambdavec[z-1][i],6).helpfill(8) + "   " + ST::doubletostring(krit_fkt[i],6) + "\n");
     genoptions_mult[0]->out("\n");
     }
  }

//-----------------------------------------------------------------------------------

double STEPWISErun::criterion_min(const double & df)
  {
  double df1 = df;

  double kriterium;
  if(criterion=="GCV")
    kriterium = likep_mult[0]->compute_gcv(df1);
  else if(criterion=="GCV2")
    kriterium = likep_mult[0]->compute_gcv2(df1);
  else if(criterion=="AIC")
    kriterium = likep_mult[0]->compute_aic(df1);
  else if(criterion=="BIC")
    kriterium = likep_mult[0]->compute_bic(df1);
  else if(criterion=="AIC_imp")
    kriterium = likep_mult[0]->compute_improvedaic(df1);
  else if(criterion=="MSEP" || criterion=="CV5" || criterion=="CV10")
    kriterium = likep_mult[0]->compute_msep();
  else //if(criterion=="AUC")
    kriterium = -1 * likep_mult[0]->compute_auc();

  if(criterion=="CV5" || criterion=="CV10")
    kriterium = kriterium / likep_mult[0]->get_nrobs();
  else if(criterion=="MSEP")
    kriterium = kriterium /(likep_mult[0]->get_nrobs() - likep_mult[0]->get_nrobs_wpw());

  return kriterium;
  }

double STEPWISErun::criterion_min(const double & df, const ST::string & auswahl)
  {
  double df1 = df;

  double kriterium;
  if(auswahl=="GCV")
    kriterium = likep_mult[0]->compute_gcv(df1);
  else if(auswahl=="GCV2")
    kriterium = likep_mult[0]->compute_gcv2(df1);
  else if(auswahl=="AIC")
    kriterium = likep_mult[0]->compute_aic(df1);
  else if(auswahl=="BIC")
    kriterium = likep_mult[0]->compute_bic(df1);
  else if(auswahl=="AIC_imp")
    kriterium = likep_mult[0]->compute_improvedaic(df1);
  else if(auswahl=="MSEP")
    kriterium = likep_mult[0]->compute_msep();
  //else //if(auswahl=="AUC")
  //  kriterium = -1 * likep_mult[0]->compute_auc();

  return kriterium;
  }

// -----------------------------------------------------------------------------
// ------------------ Funktionen für Koordinatenmethode ------------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::koordabstieg(void)
  {
      // Schleife für Minimierung
  ST::string tr_akt = "trace_on";
  ST::string text_neu;
  bool eins = true;
  double kriterium_aktuell;
  while(kriterium_neu <= kriterium_alt && fertig==false && steps_aktuell<steps)
       {
       steps_aktuell = steps_aktuell + 1;
       vector<ST::string> textiteration;
       fertig = true;
       kriterium_aktuell = kriterium_neu;
       kriterium_alt = kriterium_neu;
       modell_alt = modell_neu;

       ST::string header = "  Startmodel:";
       maketext(header,modell_alt,kriterium_alt,text_neu,true,tr_akt,true);
       text_alt = text_neu;

       if(eins == true)
         {
         eins = false;
         if(trace == "trace_off")
           tr_akt = "trace_off";
         }

       vector<vector<double> > modeliteration;
       vector<double> kriteriumiteration2;

       if(minim == "exact")
         {
         unsigned z = koordexact_fixfactor(kriteriumiteration2,modeliteration,
                            textiteration,kriterium_aktuell);
         koordexact_nonp(kriteriumiteration2,modeliteration,textiteration,z,kriterium_aktuell);
         }

       else
         {
         koord_minfix(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell);
         unsigned z = koord_minfactor(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell);
         koord_minnonp(kriteriumiteration2,modeliteration,textiteration,z,kriterium_aktuell);

         if((minim == "adaptiv" || minim == "adap_exact")
                   && likep_mult[0]->get_family() != "Gaussian")
           {
           //likep_mult[0]->compute_iwls();
           likep_mult[0]->posteriormode();
           likep_mult[0]->compute_iwls();
           if(likep_mult[0]->iwlsweights_constant() == false)
             {
             for(unsigned y=0;y<fullcond_alle.size();y++)
               fullcond_alle[y]->set_calculate_xwx();
             }
           }

         if(isboot == false && (minim == "adaptiv" || minim == "adap_exact")
                                     && modellematrix.size()>=3)
           {
           unsigned hilfe = modellematrix.size()-3;
           if(modell_alt == modellematrix[hilfe][modellematrix[hilfe].size()-1])
              fertig = true;
           }
         else if(isboot == true && (minim == "adaptiv" || minim == "adap_exact")
                                     && modellematrix.size()>=2)
           {
           unsigned hilfe = modellematrix.size()-2;
           if(modell_alt == modellematrix[hilfe][modellematrix[hilfe].size()-1])
              fertig = true;
           }
         }

       kriterium_neu = kriterium_aktuell;

       if(fertig==false)
         {
         if(isboot==false)
           {
           outcriterium << ST::inttostring(steps_aktuell) << "   0   " << ST::doubletostring(kriterium_neu,8) << endl;
           outmodels << ST::inttostring(steps_aktuell) << "   " << ST::doubletostring(kriterium_neu,8) << "   ";
           }
         modellematrix.push_back(modeliteration);

         if(steps_aktuell==steps)
           {
           if(kriterium_alt>kriterium_neu)
             {
             kriterium_alt = kriterium_neu;
             modell_alt = modell_neu;
             text_alt = text_neu;
             }
           }
         }
       else
         {
         if(trace == "trace_on" || trace == "trace_minim")
           {
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("\n\n");
           genoptions_mult[0]->out("  There are no new models for this iteration! \n");
           }
         if(isboot==false)
           {
           outcriterium << ST::inttostring(steps_aktuell) << "   0   " << ST::doubletostring(kriterium_neu,8) << endl;
           outmodels << ST::inttostring(steps_aktuell) << "   " << ST::doubletostring(kriterium_neu,8) << "   ";
           }
         }
       if(trace == "trace_on" || trace == "trace_minim")
         {
         genoptions_mult[0]->out("\n\n");
         genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
         genoptions_mult[0]->out("  ------------------------------------------------------------------------ \n");
         }

       if(make_pause() == true)
          return true;
       }

  if(minim == "apprexact" || minim == "adap_exact")
    {
    if(fertig==false)
      {
      modell_neu = modell_alt;
      fix_komplett(modell_alt);
      fullcond_komplett(modell_alt);
      if(minim == "adap_exact")
        {
        schaetzen(0,kriterium_alt,true,"backfitting");
        }
      kriterium_neu = kriterium_alt;
      }
    else if(fertig == true && minim == "adap_exact")
      {
      schaetzen(0,kriterium_alt,true,"backfitting");
      kriterium_neu = kriterium_alt;
      }
    minim = "exact";
    //if(trace == "trace_on" || trace == "trace_minim")
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out("  Beginning of the exact minimization! \n");
    fertig = false;
    modellematrix.erase(modellematrix.begin(),modellematrix.end());
    vector<vector<double> > startiteration;
    startiteration.push_back(modell_alt);
    modellematrix.push_back(startiteration);
    koordabstieg();
    }

  return false;
  }


void STEPWISErun::koord_minfix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell)
  {
  unsigned i;
  for(i=1;i<names_fixed.size();i++)
    {
    if(modell_alt[i-1]==-1)
      koord_fix_leer(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,i);
    else if(modell_alt[i-1]==0)
      koord_leer_fix(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,i);
    modell_alt = modell_neu;
    outcriterium << ST::inttostring(steps_aktuell-1) << "   " << i << "   "
                 << ST::doubletostring(kriterium_aktuell,8) << endl;
    }
  }

void STEPWISErun::koord_fix_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & i)
  {
  double kriterium_adaptiv = kriterium_aktuell;
  if(minim == "adaptiv" || minim == "adap_exact"
                        || criterion == "CV5" || criterion == "CV10")
    {
    schaetzen(i,kriterium_aktuell,true,"fix");
    }

  modell_neu[i-1] = 0;
  fullcond_alle[0]->safe_const();
  reset_fix(names_fixed[i]);
  schaetzen(0,kriterium_neu,false,"leer");
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    include_fix(names_fixed[i]);
    posteriormode(posttitle,true);
    reset_fix(names_fixed[i]);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

 if(minim != "adaptiv" && minim != "adap_exact")
   {
   if(kriterium_neu < kriterium_aktuell)
     {
     kriterium_aktuell = kriterium_adaptiv;
     bool neu = modelcomparison(modell_neu,modellematrix);
     if(neu==false)
       {
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
       }
     else
       kriterium_neu = kriterium_aktuell;

     if(neu==true || kriterium_aktuell < kriterium_neu)
       {
       int c = column_for_fix(names_fixed[i]);
       vector<ST::string> name_help;
       name_help.push_back(names_fixed[i]);
       fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
       modell_neu[i-1] = -1;
       if(kriterium_aktuell < kriterium_neu) // verhindert, daß "approx" schlechter wird!
         {
         posteriormode(posttitle,true);
         if( (trace == "trace_minim" || trace == "trace_on") && neu == false)
             genoptions_mult[0]->out("\n\n  Trial won't become the new model! \n");
         }
       }
     else
       kriterium_aktuell = kriterium_neu;
     }
   else //if(kriterium_neu > kriterium_aktuell)
     {
     kriterium_neu = kriterium_adaptiv;
     kriterium_aktuell = kriterium_adaptiv;
     int c = column_for_fix(names_fixed[i]);
     vector<ST::string> name_help;
     name_help.push_back(names_fixed[i]);
     fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
     modell_neu[i-1] = -1;
     }
   }

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    if(kriterium_aktuell >= kriterium_neu)
      kriterium_aktuell = kriterium_neu;
    else //if(kriterium_neu > kriterium_aktuell)
      {
      kriterium_neu = kriterium_aktuell;
      int c = column_for_fix(names_fixed[i]);
      vector<ST::string> name_help;
      name_help.push_back(names_fixed[i]);
      fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
      modell_neu[i-1] = -1;
      }

    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= std::pow(10,-6.0))
      fertig = false;
    if(modell_alt[i-1] != modell_neu[i-1] && (trace == "trace_on" || trace == "trace_minim"))
      {
      ST::string text;
      maketext("  Trial:",modell_neu,kriterium_aktuell,text,true,trace,false);
      }
    kriterium_alt = kriterium_aktuell;
    modell_alt[i-1] = modell_neu[i-1];
    modeliteration.push_back(modell_alt);
    }
  }

void STEPWISErun::koord_leer_fix(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & i)
  {
  double kriterium_adaptiv = kriterium_aktuell;
  // if(minim == "adaptiv" || minim == "adap_exact")
     // --> hier nicht nötig, weil Intercept immer erneuert wird!
     // --> bei 1. Koeffizienten und "Gamma" doch nötig, weil Phi neu geschätzt wurde!
  if(criterion == "CV5" || criterion == "CV10"
      || (minim == "adaptiv" || minim == "adap_exact") && likep_mult[0]->get_family()=="Gamma")
    {
    schaetzen(i,kriterium_aktuell,true,"leer");
    }

  modell_neu[i-1] = -1;
  fullcond_alle[0]->safe_const();
  schaetzen(i,kriterium_neu,false,"fix");
  fullcond_alle[0]->set_const_old();
  if(minim == "approx_control")
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_fixed[i] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    reset_fix(names_fixed[i]);
    int c = column_for_fix(names_fixed[i]);
    vector<ST::string> name_help;
    name_help.push_back(names_fixed[i]);
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(name_help,datamatrix(D.getCol(c)),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_fixed[i] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim != "adaptiv" && minim != "adap_exact")
    {
    if(kriterium_neu < kriterium_aktuell)
      {
      kriterium_aktuell = kriterium_adaptiv;
      bool neu = modelcomparison(modell_neu,modellematrix);
      if(neu==false)
        {
        newmodel(kriteriumiteration2,modeliteration,textiteration);
        kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
        }
      else
        kriterium_neu = kriterium_aktuell;

      if(neu==true || kriterium_aktuell < kriterium_neu)
        {
        reset_fix(names_fixed[i]);
        modell_neu[i-1] = 0;
        if(kriterium_aktuell < kriterium_neu)
          {
          posteriormode(posttitle,true);
          if( (trace == "trace_minim" || trace == "trace_on") && neu == false)
             genoptions_mult[0]->out("\n\n  Trial won't become the new model! \n");
          }
        }
      else
        kriterium_aktuell = kriterium_neu;
      }
    else //if(kriterium_neu >= kriterium_aktuell)
      {
      reset_fix(names_fixed[i]);
      modell_neu[i-1] = 0;
      }
    }

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    if(kriterium_aktuell >= kriterium_neu)
      kriterium_aktuell = kriterium_neu;
    else  // if(kriterium_neu >= kriterium_aktuell)
      {
      reset_fix(names_fixed[i]);
      modell_neu[i-1] = 0;
      fullcond_alle[0]->posteriormode_const();
      }

    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= std::pow(10,-6.0))
      fertig = false;
    if(modell_alt[i-1] != modell_neu[i-1] && (trace == "trace_on" || trace == "trace_minim"))
      {
      ST::string text;
      maketext("  Trial:",modell_neu,kriterium_aktuell,text,true,trace,false);
      }
    kriterium_alt = kriterium_aktuell;
    modell_alt[i-1] = modell_neu[i-1];
    modeliteration.push_back(modell_alt);
    }
  }

unsigned STEPWISErun::koord_minfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell)
  {
  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
     {
     ST::string possible = "alles";
     // VCM
     if(hierarchical)
       fullcond_alle[z]->hierarchical(possible);

     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false
                                               && possible == "alles")
       koord_factor_leer(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,z);
     else if(modell_alt[z+names_fixed.size()-2]==-1 && (fullcond_alle[z]->get_forced()==true
                                               || possible == "vfix"))
       {
       if(minim == "adaptiv" || minim == "adap_exact")
                 // || criterion == "CV5" || criterion == "CV10")  --> hier nicht nötig!
         {
         kriterium_aktuell = MAXDOUBLE;
         koord_factor_leer(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,z);
         }
       }
     else if(modell_alt[z+names_fixed.size()-2]==0)
       {
       // VCM
       if(possible == "vfix")
         {
         for(unsigned i=0;i<names_nonp[z-1].size();i++)
           reset_fix(names_nonp[z-1][i]);
         kriterium_aktuell = MAXDOUBLE;
         fullcond_alle[z]->set_inthemodel(0);
         }
       koord_leer_factor(kriteriumiteration2,modeliteration,textiteration,kriterium_aktuell,z);
       }
     modell_alt = modell_neu;
     outcriterium << ST::inttostring(steps_aktuell-1) << "   " << ST::inttostring(names_fixed.size()-1+z)
                  << "   " << ST::doubletostring(kriterium_aktuell,8) << endl;
     z = z + 1;
     }
  return z;
  }

void STEPWISErun::koord_factor_leer(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z)
  {
  unsigned i;
  vector<FULLCOND*> fullcond_ori = fullcondp;
  fullcondp = fullcond_ori;

  double kriterium_adaptiv = kriterium_aktuell;
  if(minim == "adaptiv" || minim == "adap_exact"
                        || criterion == "CV5" || criterion == "CV10")
    {
    schaetzen(z,kriterium_aktuell,true,"factor");
    }

  if(kriterium_adaptiv < MAXDOUBLE)
    {
    modell_neu[z+names_fixed.size()-2] = 0;
    fullcond_alle[z]->set_inthemodel(0);
    fullcond_alle[0]->safe_const();
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    schaetzen(0,kriterium_neu,true,"leer");
    fullcond_alle[0]->set_const_old();

    if(minim == "approx_control")
      {
      double kriterium_test;
      schaetzen(-1,kriterium_test,false,"backfitting");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
            + ST::doubletostring(kriterium_neu,6) + " exact = "
            + ST::doubletostring(kriterium_test,6) + "\n");
      fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
      korrektur();  // fullcondp[0]->posteriormode_const();
      posteriormode(posttitle,true);
      for(i=0;i<names_nonp[z-1].size();i++)
        reset_fix(names_nonp[z-1][i]);
      }
    if(trace == "trace_minim" && minim != "approx_control")
      {
      genoptions_mult[0]->out("\n\n");
      genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
      genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
      genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
      genoptions_mult[0]->out("\n");
      }

    if(minim != "adaptiv" && minim != "adap_exact")
      {
      if(kriterium_neu < kriterium_aktuell)
        {
        kriterium_aktuell = kriterium_adaptiv;
        bool neu = modelcomparison(modell_neu,modellematrix);
        if(neu==false)
          {
          newmodel(kriteriumiteration2,modeliteration,textiteration);
          kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
          }
        else
          kriterium_neu = kriterium_aktuell;

        if(neu==true || kriterium_aktuell < kriterium_neu)
          {
          fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                         fullcond_alle[z]->get_data_forfixedeffects(),true);
          modell_neu[z+names_fixed.size()-2] = -1;
          fullcond_alle[z]->set_inthemodel(-1);
          if(kriterium_aktuell < kriterium_neu)
            {
            posteriormode(posttitle,true);
            if( (trace == "trace_minim" || trace == "trace_on") && neu == false)
              genoptions_mult[0]->out("\n\n  Trial won't become the new model! \n");
            }
          }
        else
          kriterium_aktuell = kriterium_neu;
        }
      else //if(kriterium_neu >= kriterium_aktuell)
        {
        kriterium_aktuell = kriterium_adaptiv;
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                     fullcond_alle[z]->get_data_forfixedeffects(),true);
        modell_neu[z+names_fixed.size()-2] = -1;
        fullcond_alle[z]->set_inthemodel(-1);
        }
      }

    if(minim == "adaptiv" || minim == "adap_exact")
      {
      if(kriterium_aktuell >= kriterium_neu)
        kriterium_aktuell = kriterium_neu;
      else //if(kriterium_neu >= kriterium_aktuell)
        {
        fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                     fullcond_alle[z]->get_data_forfixedeffects(),true);
        modell_neu[z+names_fixed.size()-2] = -1;
        fullcond_alle[z]->set_inthemodel(-1);
        }

      if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= std::pow(10,-6.0))
        fertig = false;
      if(modell_alt[z+names_fixed.size()-2] != modell_neu[z+names_fixed.size()-2]
                      && (trace == "trace_on" || trace == "trace_minim"))
        {
        ST::string text;
        maketext("  Trial:",modell_neu,kriterium_aktuell,text,true,trace,false);
        }
      kriterium_alt = kriterium_aktuell;
      modell_alt[z+names_fixed.size()-2] = modell_neu[z+names_fixed.size()-2];
      modeliteration.push_back(modell_alt);
      }
    }
  else
    {
    if(trace == "trace_minim" && minim != "approx_control")
      {
      genoptions_mult[0]->out("\n\n");
      genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
      genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
      genoptions_mult[0]->out("\n");
      }
    }
  }

void STEPWISErun::koord_leer_factor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell, unsigned & z)
  {
  unsigned i;
  vector<FULLCOND*> fullcond_ori = fullcondp;
  fullcondp = fullcond_ori;

  double kriterium_adaptiv = kriterium_aktuell;
  // if(minim == "adaptiv" || minim == "adap_exact")
      // ---> hier überflüssig (siehe oben)!
  if(criterion == "CV5" || criterion == "CV10"
     || (minim == "adaptiv" || minim == "adap_exact") && likep_mult[0]->get_family()=="Gamma")
    {
    schaetzen(z,kriterium_aktuell,true,"leer");
    }

  modell_neu[z+names_fixed.size()-2] = -1;
  fullcond_alle[z]->set_inthemodel(-1);
  fullcond_alle[0]->safe_const();
  schaetzen(z,kriterium_neu,false,"factor");
  fullcond_alle[0]->set_const_old();

  if(minim == "approx_control" && kriterium_aktuell < MAXDOUBLE)
    {
    double kriterium_test;
    schaetzen(-1,kriterium_test,false,"backfitting");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  " + names_nonp[z-1][0] + " Testvalue: approx = "
          + ST::doubletostring(kriterium_neu,6) + " exact = "
          + ST::doubletostring(kriterium_test,6) + "\n");
    for(i=0;i<names_nonp[z-1].size();i++)
      reset_fix(names_nonp[z-1][i]);
    korrektur();  // fullcondp[0]->posteriormode_const();
    posteriormode(posttitle,true);
    fullcond_alle[0]->posteriormode_single(names_nonp[z-1],
                                  fullcond_alle[z]->get_data_forfixedeffects(),true);
    }
  if(trace == "trace_minim" && minim != "approx_control")
     {
     genoptions_mult[0]->out("\n\n");
     genoptions_mult[0]->out("  " + names_nonp[z-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  Lambda   Testvalue (approx): \n");
     genoptions_mult[0]->out(" " + ST::doubletostring(-1).helpfill(8) + "   " + ST::doubletostring(kriterium_neu,6) + "\n");
     genoptions_mult[0]->out(" " + ST::doubletostring(0).helpfill(8) + "   " + ST::doubletostring(kriterium_aktuell,6) + "\n");
     genoptions_mult[0]->out("\n");
     }

  if(minim != "adaptiv" && minim != "adap_exact")
    {
    if(kriterium_neu < kriterium_aktuell)
      {
      kriterium_aktuell = kriterium_adaptiv;
      bool neu = modelcomparison(modell_neu,modellematrix);
      if(neu==false)
        {
        newmodel(kriteriumiteration2,modeliteration,textiteration);
        kriterium_neu = kriteriumiteration2[kriteriumiteration2.size()-1];
        }
      else
        kriterium_neu = kriterium_aktuell;

      if(neu==true || kriterium_aktuell < kriterium_neu)
        {
        for(i=0;i<names_nonp[z-1].size();i++)
          reset_fix(names_nonp[z-1][i]);
        modell_neu[z+names_fixed.size()-2] = 0;
        if(kriterium_aktuell < kriterium_neu)
          {
          posteriormode(posttitle,true);
          if( (trace == "trace_minim" || trace == "trace_on") && neu == false)
             genoptions_mult[0]->out("\n\n  Trial won't become the new model! \n");
          }
        }
      else
        kriterium_aktuell = kriterium_neu;
      }
    else // if(kriterium_neu >= kriterium_aktuell)
      {
      for(i=0;i<names_nonp[z-1].size();i++)
        reset_fix(names_nonp[z-1][i]);
      modell_neu[z+names_fixed.size()-2] = 0;
      fullcond_alle[z]->set_inthemodel(0);
      }
    }

  if(minim == "adaptiv" || minim == "adap_exact")
    {
    if(kriterium_aktuell >= kriterium_neu)
      kriterium_aktuell = kriterium_neu;
    else  // if(kriterium_neu >= kriterium_aktuell)
      {
      for(i=0;i<names_nonp[z-1].size();i++)
        reset_fix(names_nonp[z-1][i]);
      modell_neu[z+names_fixed.size()-2] = 0;
      fullcond_alle[z]->set_inthemodel(0);
      fullcond_alle[0]->posteriormode_const();
      }

    if(fabs((kriterium_adaptiv - kriterium_aktuell)/kriterium_adaptiv) >= std::pow(10,-6.0))
      fertig = false;
    if(modell_alt[z+names_fixed.size()-2] != modell_neu[z+names_fixed.size()-2]
                    && (trace == "trace_on" || trace == "trace_minim"))
      {
      ST::string text;
      maketext("  Trial:",modell_neu,kriterium_aktuell,text,true,trace,false);
      }
    kriterium_alt = kriterium_aktuell;
    modell_alt[z+names_fixed.size()-2] = modell_neu[z+names_fixed.size()-2];
    modeliteration.push_back(modell_alt);
    }
  }

void STEPWISErun::koord_minnonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z, double & kriterium_aktuell)
  {
  unsigned i;
  for(i=z;i<fullcond_alle.size();i++)
    {
    unsigned lambda_ind;
    double kriterium_min;
    double kriterium_test = kriterium_aktuell;

    unsigned y;
    for(y=1;y<fullcond_alle.size();y++)
      fullcond_alle[y]->set_center(false);    // sorgt dafür, daß Funktionen nicht zentriert werden!

    vector<double> krit_fkt;
    if(modell_alt[i+names_fixed.size()-2]==0)
      stepmin_nonp_leer(i,krit_fkt,kriterium_aktuell);
    else if(modell_alt[i+names_fixed.size()-2]==-1)
      stepmin_nonp_fix(i,krit_fkt,kriterium_aktuell);
    else
      stepmin_nonp_nonp(i,krit_fkt,kriterium_aktuell);

    kriterium_min = krit_fkt[0];
    unsigned j;
    lambda_ind = 0;
    for(j=1;j<krit_fkt.size();j++)  //berechnet den besten Wert
       {
       if(krit_fkt[j]<=kriterium_min)
         {
         kriterium_min = krit_fkt[j];
         lambda_ind = j;
         }
       }

    for(y=1;y<fullcond_alle.size();y++)
     {
     if(fullcond_alle[y]->is_identifiable() == false)
       fullcond_alle[y]->set_center(true);    // sorgt dafür, daß Funktionen zentriert werden!
     }

    modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][lambda_ind];
    if(minim != "adaptiv" && minim != "adap_exact")
      {
      kriterium_aktuell = kriterium_test;
      if(modell_neu[names_fixed.size()-2+i] != modell_alt[names_fixed.size()-2+i])
        {
        bool neu = modelcomparison(modell_neu,modellematrix);
        fullcond_einzeln(modell_neu,modell_alt,i);
        if(neu==false)
          {
          korrektur();  // fullcondp[0]->posteriormode_const();
          newmodel(kriteriumiteration2,modeliteration,textiteration);
          kriterium_test = kriteriumiteration2[kriteriumiteration2.size()-1];
          }
        if(kriterium_aktuell > kriterium_test)
          {
          modell_alt = modell_neu;
          kriterium_aktuell = kriterium_test;
          }
        else //if(kriterium_aktuell <= kriterium_test)
          {
          if( (trace == "trace_minim" || trace == "trace_on") && neu == false)
            genoptions_mult[0]->out("\n\n  Trial won't become the new model! \n");
          fullcond_einzeln(modell_alt,modell_neu,i);
          modell_neu = modell_alt;
          posteriormode(posttitle,true);
          }
        }
      /*else
        {
        posteriormode(posttitle,true);  // nur Versuch!!!
        }*/
      }
    else //if(minim == "adaptiv" || minim == "adap_exact")
      {
      if(modell_alt[names_fixed.size()-2+i] != modell_neu[names_fixed.size()-2+i])
        {
        fullcond_einzeln(modell_neu,modell_alt,i);

        vector<FULLCOND*> fullcond_ori = fullcondp;
        fullcondp = fullcond_ori;

        if(modell_neu[names_fixed.size()-2+i] == 0)      // noch mal überprüfen!!!
          {
          //fullcond_alle[i]->reset_effect(0);   // Nicht nötig, wegen "fullcond_einzeln"-> fullcond_alle[i] nicht in fullcondp!!!
          if(modell_alt[names_fixed.size()-2+i] != 0)
            {
            if(hierarchical)                        // neu für VCM!  Versuch!!!
              {
              ST::string possible = "alles";
              fullcond_alle[i]->hierarchical(possible);
              if((possible == "valles" || possible == "vrfix" || possible == "vraus")
                   && fullcond_alle[i]->is_identifiable() == false)
                {
                vector<ST::string> help;
                help.push_back(fullcond_alle[i]->get_datanames()[1]);
                fullcond_alle[i]->posteriormode_single(help,
                                 fullcond_alle[i]->get_data_forfixedeffects(),false);
                }
              }
            }
          fullcond_alle[0]->posteriormode_const();
          }
        else if(modell_neu[names_fixed.size()-2+i] == -1)
          {
          fullcond_alle[i]->reset_effect(0);
          fullcond_alle[0]->posteriormode_single(names_nonp[i-1],
                               fullcond_alle[i]->get_data_forfixedeffects(),false);
          }
        else
          {
          if(modell_alt[names_fixed.size()-2+i] == 0)   // approximiert Intercept bei VCs (g(x)*z)
            fullcond_alle[i]->const_varcoeff();         // durch b0 aus "b0 + b1*z"
          fullcond_alle[i]->update_stepwise(modell_neu[names_fixed.size()-2+i]);
          fullcond_alle[i]->posteriormode();
          fullcond_alle[0]->update_linold();

          ST::string possible = "alles";               // setzt Intercept bei VC so, dass eta_quer = y_quer
          fullcond_alle[i]->hierarchical(possible);
          if(possible == "valles" || possible == "vrfix")
            fullcond_alle[0]->posteriormode_const();
          }

        if(trace == "trace_on" || trace == "trace_minim")
          {
          ST::string text;
          maketext("  Trial:",modell_neu,kriterium_min,text,true,trace,false);
          }

        kriterium_aktuell = kriterium_min;
        }
      //else
      //  maketext(header,modell_neu,kriterium_alt,text_neu,true,tr_akt,true);
      modell_alt = modell_neu;
      kriterium_alt = kriterium_aktuell;
      if(fabs((kriterium_test - kriterium_aktuell)/kriterium_test) >= std::pow(10,-6.0))
        fertig = false;
      modeliteration.push_back(modell_alt);
      }
    outcriterium << ST::inttostring(steps_aktuell-1) << "   " << ST::inttostring(names_fixed.size()-1+i)
                 << "   " << ST::doubletostring(kriterium_aktuell,8) << endl;
    }
  }


//------------------------------------------------------------------------------

unsigned STEPWISErun::koordexact_fixfactor(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      double & kriterium_aktuell)
  {
  bool help_fertig = true;
  unsigned i;
  for(i=1;i<names_fixed.size();i++)
     {
     if(modell_alt[i-1]==-1)
       modell_neu[i-1]= 0;
     else if(modell_alt[i-1]==0)
       modell_neu[i-1] = -1;
     if(modelcomparison(modell_neu,modellematrix)==false)
       {
       if(modell_neu[i-1]==0)
         reset_fix(names_fixed[i]);
       else
         include_fix(names_fixed[i]);
       korrektur();  // fullcondp[0]->posteriormode_const();
       newmodel(kriteriumiteration2,modeliteration,textiteration);
       if(kriteriumiteration2[kriteriumiteration2.size()-1] > kriterium_aktuell)
          {
          if(modell_neu[i-1]==0)
            include_fix(names_fixed[i]);
          else
            reset_fix(names_fixed[i]);
          modell_neu = modell_alt;
          }
       else
          {
          help_fertig = false;
          modell_alt = modell_neu;
          kriterium_aktuell = kriteriumiteration2[kriteriumiteration2.size()-1];
          }
       }
     else
       modell_neu = modell_alt;
     outcriterium << ST::inttostring(steps_aktuell-1) << "   " << i << "   " << ST::doubletostring(kriterium_aktuell,8) << endl;
     }

  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
     {
     ST::string possible = "alles";
     // VCM
     if(hierarchical)
       fullcond_alle[z]->hierarchical(possible);

     if(modell_alt[z+names_fixed.size()-2]==-1 && fullcond_alle[z]->get_forced()==false
                                               && possible == "alles")
       modell_neu[z+names_fixed.size()-2]= 0;
     else if(modell_alt[z+names_fixed.size()-2]==0)
       {
       if(possible == "vfix")         // VCM
         {
         for(i=0;i<names_nonp[z-1].size();i++)
           reset_fix(names_nonp[z-1][i]);
         kriterium_aktuell = MAXDOUBLE;
         }
       modell_neu[z+names_fixed.size()-2] = -1;
       }
     if(modelcomparison(modell_neu,modellematrix)==false)
        {
        if(modell_neu[z+names_fixed.size()-2]==0)
          {
          for(i=0;i<names_nonp[z-1].size();i++)
            reset_fix(names_nonp[z-1][i]);
          }
        else
          fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
        korrektur();  // fullcondp[0]->posteriormode_const();
        newmodel(kriteriumiteration2,modeliteration,textiteration);
        if(kriteriumiteration2[kriteriumiteration2.size()-1] > kriterium_aktuell)
           {
           if(modell_neu[z+names_fixed.size()-2]==0)
             fullcondp[0]->include_effect(names_nonp[z-1],fullcond_alle[z]->get_data_forfixedeffects());
           else
             {
             for(i=0;i<names_nonp[z-1].size();i++)
               reset_fix(names_nonp[z-1][i]);
             }
           modell_neu = modell_alt;
           }
        else
           {
           help_fertig = false;
           modell_alt = modell_neu;
           kriterium_aktuell = kriteriumiteration2[kriteriumiteration2.size()-1];
           }
        }
     else
        modell_neu = modell_alt;
     outcriterium << ST::inttostring(steps_aktuell-1) << "   " << ST::inttostring(names_fixed.size()-1+z) << "   "
                  << ST::doubletostring(kriterium_aktuell,8) << endl;
     z = z + 1;
     }
  fertig = help_fertig;
  return z;
  }

void STEPWISErun::koordexact_nonp(vector<double> & kriteriumiteration2,
      vector<vector<double> > & modeliteration, vector<ST::string> & textiteration,
      unsigned & z, double & kriterium_aktuell)
  {
  unsigned i;
  unsigned j;
  for(i=z;i<fullcond_alle.size();i++)
    {
    modell_neu = modell_alt;
    unsigned lambda_ind;
    vector<double> krit_fkt;
    if(modell_alt[i+names_fixed.size()-2]==0)
      minexact_nonp_leer(i,krit_fkt,kriterium_aktuell);
    else if(modell_alt[i+names_fixed.size()-2]==-1)
      {
      reset_fix(names_nonp[i-1][0]);
      minexact_nonp_fix(i,krit_fkt,kriterium_aktuell);
      }
    else
      minexact_nonp_nonp(i,krit_fkt,kriterium_aktuell);

    double kriterium_min = krit_fkt[0];
    lambda_ind = 0;
    for(j=1;j<krit_fkt.size();j++)  //berechnet den besten Wert
       {
       if(krit_fkt[j]<=kriterium_min)
         {
         kriterium_min = krit_fkt[j];
         lambda_ind = j;
         }
       }

    modell_neu[names_fixed.size()-2+i] = lambdavec[i-1][lambda_ind];
    if(modell_neu[names_fixed.size()-2+i] != modell_alt[names_fixed.size()-2+i])
      {
      if(modelcomparison(modell_neu,modellematrix)==false)
        {
        fullcond_einzeln(modell_neu,modell_alt,i);
        korrektur();  // fullcondp[0]->posteriormode_const();
        newmodel(kriteriumiteration2,modeliteration,textiteration);
        kriterium_aktuell = kriteriumiteration2[kriteriumiteration2.size()-1];
        }
      else
        modell_neu = modell_alt;
      }
    modell_alt = modell_neu;
    outcriterium << ST::inttostring(steps_aktuell-1) << "   " << ST::inttostring(names_fixed.size()-1+i) << "   "
                 << ST::doubletostring(kriterium_aktuell,8) << endl;
    }
  }


// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des Startmodels -----------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::vcm_doppelt(void)
  {
  unsigned i = 1;
  unsigned j;
  unsigned k;
  bool fehler_randomslope = false;
  while(i<fullcond_alle.size() && fehler_randomslope == false)
     {
     j = 1;
     while(j<names_fixed.size() && fehler_randomslope == false)
       {
       if(names_fixed[j]==names_nonp[i-1][0])
         {
         k = j;
         fehler_randomslope = true;
         }
       j++;
       }
     i++;
     }
  if(fehler_randomslope==true)
    {
    genoptions_mult[0]->outerror("\n\n ERROR: You must not put fixed effect " +
       names_fixed[k] + " in the model! \n");
    return true;
    }
  else
    return false;
  }


void STEPWISErun::initialise_lambdas(vector<vector<ST::string> > & namen_nonp,
       vector<ST::string> & namen_fix, vector<vector<double> > & lambdavector,
       const int & number, const bool & gewichte)
  {
  namen_fix = fullcondp[0]->get_datanames();
  unsigned i;

       // berechnet ein Modell, um Gewichte für Augabe usw. zu erhalten!!!
  if(gewichte == true)
    {
    vector<double> modell_init;
    for(i=1;i<namen_fix.size();i++)
       modell_init.push_back(-1);
    for(i=1;i<fullcond_alle.size();i++)
      {
      if(fullcond_alle[i]->get_fctype() != MCMC::factor)
         {
         fullcond_alle[i]->update_stepwise(100);
         fullcond_alle[i]->set_inthemodel(1);
         }
      else if(fullcond_alle[i]->get_fctype() == MCMC::factor)
         {
         fullcondp.erase(fullcondp.begin() + 1);   //löscht das Element an Pos.1 aus fullcondp-Vektor
         fullcondp[0]->include_effect(fullcond_alle[i]->get_datanames(),
                                   fullcond_alle[i]->get_data_forfixedeffects());
         fullcond_alle[i]->set_inthemodel(-1);
         }
      }
    end[0] = fullcondp.size()-1;

//likep_mult[0]->compute_iwls();
    posteriormode(posttitle,true);
    fullcondp = fullcond_alle;
    end[0] = fullcondp.size()-1;
    }
  else
   likep_mult[0]->compute_iwls();

  for(i=1;i<fullcond_alle.size();i++)
     {
     int nummer = number;
     if(fullcond_alle[i]->get_data_forfixedeffects().cols() > 1)  //bei Faktor-Variablen
        namen_nonp.push_back(fullcond_alle[i]->get_datanames());
     else
        {
        vector<ST::string> names_help;
        names_help.push_back(fullcond_alle[i]->get_datanames()[fullcond_alle[i]->get_datanames().size()-1]);
        namen_nonp.push_back(names_help);
        nummer = fullcond_alle[i]->get_number();
        if(nummer==0)
           nummer = number;
        fullcond_alle[i]->set_stepwise_options(fullcond_alle[i]->get_lambdastart(),
                       fullcond_alle[i]->get_lambdamax(),fullcond_alle[i]->get_lambdamin(),
                       fullcond_alle[i]->get_forced(),fullcond_alle[i]->get_df_lambdamax(),
                       fullcond_alle[i]->get_df_lambdamin(),fullcond_alle[i]->get_spfromdf(),
                       nummer,fullcond_alle[i]->get_df_equidist());
        }
     vector<double> untervector;
     fullcond_alle[i]->set_inthemodel(1);
     if(fullcond_alle[i]->get_matrixnumber() == 1)
       fullcond_alle[i]->compute_lambdavec(untervector,nummer);
     else
       untervector = lambdavector[lambdavector.size()-1];
     fullcond_alle[i]->set_inthemodel(0);
     lambdavector.push_back(untervector);
     }
  }


void STEPWISErun::initialise_weights(double prop)
  {
  if(criterion == "MSEP" || criterion == "AUC")
    {
    datamatrix weight = likep_mult[0]->get_weight();
    bool gewicht = false;
    if(weight.min(0) > 0)     // wenn Gewichtsmatrix 0er enthält, werden diese 0er für die Aufteilung verwendet.
      {                       // sonst: bestimmte Werte sollen im Testdatensatz enthalten sein (z.B. Min., Max. bei P-Splines)
      weight = datamatrix(weight.rows(),1,0);
      unsigned i;
      for(i=1;i<fullcond_alle.size();i++)
        fullcond_alle[i]->create_weight(weight);
      gewicht = true;
      }
    likep_mult[0]->create_weight(weight,prop,gewicht,false);
    }
  else
    {
    datamatrix w = datamatrix(1,1,0);
    bool gewicht = false;
    if(criterion== "CV5")
      likep_mult[0]->create_weight(w,5,gewicht,true);
    else
      likep_mult[0]->create_weight(w,10,gewicht,true);
    }
  }


unsigned STEPWISErun::search_lambdaindex(const double & m,
                       const vector<double> lam, bool & b) const
  {
  unsigned index = 0;
  double lambda = m;
  unsigned j=0;
  b = false;
  while(j<lam.size() && b==false)
       {
       if(lambda==lam[j])
         {
         index = j;
         b = true;
         }
       j++;
       }
  return index;
  }


unsigned STEPWISErun::search_lambdastartindex(const double & start,
                       const vector<double> & lambdas) const
  {
  bool gefunden = false;
  unsigned index = search_lambdaindex(start,lambdas,gefunden);

  if(gefunden==false)
    {
    vector<double> diff;
    unsigned i;
    for(i=0;i<lambdas.size();i++)
       {
       if(lambdas[i]!=0 && lambdas[i]!=-1)
         diff.push_back(fabs(lambdas[i]-start));
       else
         diff.push_back(MAXDOUBLE);
       }
    double irgendwas = diff[0];
    for(i=1;i<diff.size();i++)
       {
       if(diff[i]<=irgendwas)
         {
         irgendwas = diff[i];
         index = i;
         }
       }
    }

  return index;
  }


void STEPWISErun::startwerte(const ST::string & startmodel,
          vector<vector<unsigned> > & startindex,
          vector<vector<double> > & startfix)
  {
  unsigned i;

  if(startmodel == "empty" || startmodel == "both" || startmodel == "emplin")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;

    for(i=1;i<names_fixed.size();i++)
       fixhelp.push_back(0);

    for(i=1;i<fullcond_alle.size();i++)
       indexhelp.push_back(lambdavec[i-1].size()-1);

    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }

  if(startmodel == "both" || startmodel == "full")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;

    for(i=1;i<names_fixed.size();i++)
       fixhelp.push_back(-1);
    for(i=1;i<fullcond_alle.size();i++)
       indexhelp.push_back(0);

    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }

  if(startmodel == "emplin")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;

    for(i=1;i<names_fixed.size();i++)
       fixhelp.push_back(-1);
    for(i=1;i<fullcond_alle.size();i++)
       {
       unsigned index = search_lambdastartindex(-1,lambdavec[i-1]);
       indexhelp.push_back(index);
       }

    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }

  if(startmodel == "userdefined")
    {
    vector<unsigned> indexhelp;
    vector<double> fixhelp;

    for(i=1;i<names_fixed.size();i++)
       fixhelp.push_back(-1);         //fehlt: vom Benutzer angeben lassen!!!
    for(i=1;i<fullcond_alle.size();i++)
       {
       double start = fullcond_alle[i]->get_lambdastart();
       unsigned index = search_lambdastartindex(start,lambdavec[i-1]);
       indexhelp.push_back(index);
       }

    startindex.push_back(indexhelp);
    startfix.push_back(fixhelp);
    }
  }


// -----------------------------------------------------------------------------
// ------- Funktionen für die Berechnung neuer Modelle -------------------------
// -----------------------------------------------------------------------------

double STEPWISErun::compute_criterion(void)
  {
  double df = 0;
  if(criterion != "MSEP" && criterion != "AUC" && criterion != "CV5" && criterion != "CV10")
    {
    unsigned i;
    for(i=0;i<fullcond_alle.size();i++)
      df = df + fullcond_alle[i]->compute_df();
    likep_mult[0]->set_iwlsweights_notchanged(true);
    }
  double kriterium;

  if(criterion=="GCV")
    kriterium = likep_mult[0]->compute_gcv(df);
  else if(criterion=="GCV2")
    kriterium = likep_mult[0]->compute_gcv2(df);
  else if(criterion=="AIC")
    kriterium = likep_mult[0]->compute_aic(df);
  else if(criterion=="BIC")
    kriterium = likep_mult[0]->compute_bic(df);
  else if(criterion=="AIC_imp")
    kriterium = likep_mult[0]->compute_improvedaic(df);
  else if(criterion=="MSEP" || criterion=="CV5" || criterion=="CV10")
    kriterium = likep_mult[0]->compute_msep();
  else //if(criterion=="AUC")
    kriterium = -1 * likep_mult[0]->compute_auc();

  if(criterion=="CV5" || criterion=="CV10")
    kriterium = kriterium / likep_mult[0]->get_nrobs();
  else if(criterion=="MSEP")
    kriterium = kriterium /(likep_mult[0]->get_nrobs() - likep_mult[0]->get_nrobs_wpw());

  return kriterium;
  }


void STEPWISErun::newmodel(vector<double> & krit,
  vector<vector<double> > & mi, vector<ST::string> & textit)
  {
  fertig = false;
  mi.push_back(modell_neu);

/*BIC_min += 1;
if(BIC_min == 1500) //1698)
  {
  //ofstream out("c:\\bayesx\\output\\nr.txt");
  //out << modell_neu[0] << "  " << modell_neu[1] << "  " << modell_neu[2] << "  " << modell_neu[3] << "  " << modell_neu[4] << endl;
  ST::string scheisse = "Scheisse";
  } */

  double kriterium;
  schaetzen(0,kriterium,true,"backfitting");
  ST::string header = "  Trial: ";
  ST::string text;
  maketext(header,modell_neu,kriterium,text,true,trace,false);
  textit.push_back(text);
  krit.push_back(kriterium);
  }


void STEPWISErun::newmodel_fix(const double & mo, vector<double> & krit,
   vector<vector<double> > & mi, vector<ST::string> & textit,
   const ST::string & name)
  {
  if(mo==0)
    reset_fix(name);
  else
    include_fix(name);
  korrektur();  // fullcondp[0]->posteriormode_const();
  newmodel(krit,mi,textit);
  if(mo==0)
    include_fix(name);
  else
    reset_fix(name);
  //fullcondp[0]->posteriormode_const();
  }


void STEPWISErun::newmodel_factor(const double & mo, const unsigned & index,
    vector<double> & krit, vector<vector<double> > & mi,
    vector<ST::string> & textit, const vector<ST::string> & name)
  {
  unsigned i;
  if(mo==0)
    {
    for(i=0;i<name.size();i++)
      reset_fix(name[i]);
    }
  else
    fullcondp[0]->include_effect(name,fullcond_alle[index]->get_data_forfixedeffects());
  korrektur();  // fullcondp[0]->posteriormode_const();
  newmodel(krit,mi,textit);
  if(mo==0)
    fullcondp[0]->include_effect(name,fullcond_alle[index]->get_data_forfixedeffects());
  else
    {
    for(i=0;i<name.size();i++)
      reset_fix(name[i]);
    }
  //fullcondp[0]->posteriormode_const();
  }


void STEPWISErun::newmodel_nonp(const unsigned & index,
    vector<double> & krit, vector<vector<double> > & mi,
    vector<ST::string> & textit)
  {
  fullcond_einzeln(modell_neu,modell_alt,index);
  korrektur();  // fullcondp[0]->posteriormode_const();
  newmodel(krit,mi,textit);
  fullcond_einzeln(modell_alt,modell_neu,index);
  //fullcondp[0]->posteriormode_const();
  }


bool STEPWISErun::modelcomparison(const vector<double> & m,
                   const vector<vector<vector<double> > > & mmatrix)
  {
  bool s = false;
  int x = mmatrix.size()-1;
  while(x>=0 && s==false)
       {
       int y = mmatrix[x].size()-1;
       while(y>=0 && s==false)
            {
            if(mmatrix[x][y]==m)
              s = true;
            y = y - 1;
            }
       x = x - 1;
       }
  return s;
  }


// -----------------------------------------------------------------------------
// ------- Funktionen für die Erstellung des fullcondp-Vektors -----------------
// -----------------------------------------------------------------------------

void STEPWISErun::fullcond_einzeln(const vector<double> & modell1,
         const vector<double> & modell2, const unsigned & index)
  {
  vector<FULLCOND*> fullcond_neu;
  unsigned i;
  fullcond_neu.push_back(fullcondp[0]);

  unsigned z = 1;
  while(z<fullcond_alle.size() && fullcond_alle[z]->get_fctype()==factor)
    {
    z = z + 1;
    }

  for(i=z;i<fullcond_alle.size();i++)
     {
     //fullcond_alle[i]->set_inthemodel(modell1[names_fixed.size()-2+i]);
     if(modell2[names_fixed.size()-2+i]==-1 && index==i)
        reset_fix(names_nonp[i-1][0]);
     if(modell1[names_fixed.size()-2+i] != -1 && modell1[names_fixed.size()-2+i] != 0)
        {
        fullcond_neu.push_back(fullcond_alle[i]);
        if(i == index)
          fullcond_alle[i]->update_stepwise(modell1[names_fixed.size()-2+i]);
        }
     else if(modell1[names_fixed.size()-2+i]==0)    // && i == index)
        fullcond_alle[i]->reset_effect(0);
     else if(modell1[names_fixed.size()-2+i] == -1) // && i == index)
        {
        fullcond_alle[i]->reset_effect(0);
        if(i == index)
           fullcond_neu[0]->include_effect(names_nonp[i-1],
                                  fullcond_alle[i]->get_data_forfixedeffects());
        }
     }
  fullcond_alle[index]->set_inthemodel(modell1[names_fixed.size()-2+index]);

  fullcondp = fullcond_neu;
  end[0] = fullcondp.size()-1;
  }


void STEPWISErun::fullcond_komplett(const vector<double> & m)
  {
  vector<FULLCOND*> fullcond_neu;
  unsigned i;

  fullcond_neu.push_back(fullcondp[0]);
  for(i=1;i<fullcond_alle.size();i++)
     {
     fullcond_alle[i]->set_inthemodel(m[names_fixed.size()-2+i]);
     if(m[names_fixed.size()-2+i] != 0 && m[names_fixed.size()-2+i] != -1)
        {
        fullcond_alle[i]->update_stepwise(m[names_fixed.size()-2+i]);
        fullcond_neu.push_back(fullcond_alle[i]);
        }
     else if(m[names_fixed.size()-2+i]==0)
        fullcond_alle[i]->reset_effect(0);
     else if(m[names_fixed.size()-2+i] == -1)
        {
        fullcond_alle[i]->reset_effect(0);
        fullcond_neu[0]->include_effect(names_nonp[i-1],
                                  fullcond_alle[i]->get_data_forfixedeffects());
        }
     }

  fullcondp = fullcond_neu;
  end[0] = fullcondp.size()-1;
  korrektur();  // fullcondp[0]->posteriormode_const();
  }


void STEPWISErun::fix_komplett(const vector<double> &  modell)
  {
  unsigned z;
  for(z=0;z<names_fixed.size()-1;z++)
    {
    if(modell[z]==0)
      reset_fix(names_fixed[z+1]);
    else if(modell[z]==-1)
      {
      unsigned i = 1;
      bool rein = false;
      while(i<fullcondp[0]->get_datanames().size() && rein==false)
          {
          if(fullcondp[0]->get_datanames()[i]==names_fixed[z+1])
            rein = true;
          i = i + 1;
          }
      if(rein==false)
        include_fix(names_fixed[z+1]);
      }
    }

  for(z=names_fixed.size()-1;z<modell.size();z++)
    {
    bool gefunden = false;
    unsigned i = 1;
    while(i<fullcondp[0]->get_datanames().size() && gefunden==false)
      {
      if(fullcondp[0]->get_datanames()[i]==names_nonp[z-names_fixed.size()+1][0])    // ersetzen durch reset_fix?
        {
        gefunden = true;
        fullcondp[0]->reset_effect(i);
        }
      i = i + 1;
      }
    if(gefunden==true && names_nonp[z-names_fixed.size()+1].size()>1)
      {
      unsigned j;
      for(j=1;j<names_nonp[z-names_fixed.size()+1].size();j++)
         reset_fix(names_nonp[z-names_fixed.size()+1][j]);
      }
    }
  }


void STEPWISErun::fix_ganz_komplett(const vector<double> &  modell)
  {
  unsigned z;
  for(z=0;z<names_fixed.size()-1;z++)
     reset_fix(names_fixed[z+1]);

  for(z=0;z<names_fixed.size()-1;z++)
     {
     if(modell[z]==-1)
       {
       include_fix(names_fixed[z+1]);
       }
     }

  for(z=names_fixed.size()-1;z<modell.size();z++)
     {
     bool gefunden = false;
     unsigned i = 1;
     while(i<fullcondp[0]->get_datanames().size() && gefunden==false)
        {
        if(fullcondp[0]->get_datanames()[i]==names_nonp[z-names_fixed.size()+1][0])    // ersetzen durch reset_fix?
           {
           gefunden = true;
           fullcondp[0]->reset_effect(i);
           }
        i = i + 1;
        }
     if(gefunden==true && names_nonp[z-names_fixed.size()+1].size()>1)
        {
        unsigned j;
        for(j=1;j<names_nonp[z-names_fixed.size()+1].size();j++)
           reset_fix(names_nonp[z-names_fixed.size()+1][j]);
        }
     }
  }


void STEPWISErun::reset_fix(const ST::string & name)
  {
  bool raus = false;
  unsigned j = 1;
  while(j<fullcondp[0]->get_datanames().size() && raus==false)
     {
     if(fullcondp[0]->get_datanames()[j]==name)
              // || fullcondp[0]->get_datanames()[j]== name+"_1")    // neu: Zusatz für VCM???
        {
        raus = true;
        fullcondp[0]->reset_effect(j);
        }
     j = j + 1;
     }
  }


void STEPWISErun::include_fix(const ST::string & name)
  {
  int i = column_for_fix(name);
  vector<ST::string> help_name;
  help_name.push_back(name);
  fullcondp[0]->include_effect(help_name,datamatrix(D.getCol(i)));
  }

int STEPWISErun::column_for_fix(const ST::string & name)
  {
  bool gefunden = false;
  unsigned i = 0;
  while(i<modelv.size() && gefunden==false)
     {
     if(name==modelv[i])
        gefunden = true;
     i = i + 1;
     }
  return i-1;
  }


void STEPWISErun::korrektur(void)
  {
  //fullcond_alle[0]->posteriormode_const();
fullcond_alle[0]->posteriormode();
  }

// -----------------------------------------------------------------------------
// ------- Funktionen für die Ausgabe im Output-Fenster ------------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::make_pause(void)
  {
  #if defined(BORLAND_OUTPUT_WINDOW)
  Application->ProcessMessages();

  if (Frame->stop)
    {
    //break;
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("STEPWISE PROCEDURE TERMINATED BY USER BREAK\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("Estimation results: none\n");
    genoptions_mult[0]->out("\n");

    likep_mult[0]->reset();
    for(unsigned j=0;j<fullcond_alle.size();j++)
      fullcond_alle[j]->reset();

    return true;
    }

  if (Frame->pause)
    {
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("STEPWISE PROCEDURE PAUSED\n");
    genoptions_mult[0]->out("Click CONTINUE to proceed\n");
    genoptions_mult[0]->out("\n");

    while (Frame->pause)
      {
      Application->ProcessMessages();
      }

    genoptions_mult[0]->out("STEPWISE PROCEDURE CONTINUED\n");
    genoptions_mult[0]->out("\n");
    }

  #elif defined(JAVA_OUTPUT_WINDOW)
  bool stop = genoptions_mult[0]->adminb_p->breakcommand();
  if(stop)
    {
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("STEPWISE PROCEDURE TERMINATED BY USER BREAK\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("Estimation results: none\n");
    genoptions_mult[0]->out("\n");
    //break;
    likep_mult[0]->reset();
    for(unsigned j=0;j<fullcond_alle.size();j++)
      fullcond_alle[j]->reset();

    return true;
    }
  #endif

  return false;
  }


void STEPWISErun::maketext(const ST::string & h, const vector<double> & m,
                          const double & a, ST::string & text,
                          const bool & neutext, const ST::string & tr,
                          const bool & datei)
  {
if(isboot==false)
  {
  if(tr == "trace_on" || trace == "trace_minim")
    {
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out(h);
    }
  ST::string modeltext;
  if(neutext==true)
    {
    modeltext = "  " + likep_mult[0]->get_responsename() + " = ";
    unsigned i;
    modeltext = modeltext + fullcondp[0]->get_effect();
    for(i=1;i<fullcondp.size();i++)
       modeltext = modeltext + " + " + fullcondp[i]->get_effect();
    text = modeltext;
    }
  else
    {
    modeltext = text;
    }
  if(tr == "trace_on" || trace == "trace_minim")
    {
    genoptions_mult[0]->out("\n\n");
    genoptions_mult[0]->out(modeltext);
    genoptions_mult[0]->out("\n " + criterion + " = " + ST::doubletostring(a,8));
    }
  if(datei==true)
    outmodels << modeltext << endl << endl;
  }
  }


void STEPWISErun::options_text(const int & number,
         const vector<vector<double> > & startfix,
         const vector<vector<unsigned> > & startindex, const ST::string & name)
  {
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("STEPWISE OBJECT " + name + ": stepwise procedure \n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("GENERAL OPTIONS: \n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  Performance criterion: " + criterion + " \n");
  genoptions_mult[0]->out("  Maximum number of iterations: " + ST::inttostring(steps) + "\n");
  //genoptions_mult[0]->out("  Number of different smoothing parameters: "
  //   + ST::inttostring(number) + "\n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  RESPONSE DISTRIBUTION: \n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  Family: " + likep_mult[0]->get_family() + "\n");
  genoptions_mult[0]->out("  Number of observations: "
     + ST::inttostring(likep_mult[0]->get_nrobs()) + "\n");
  //genoptions_mult[0]->out("  Number of observations with positive weights: " + );
  //genoptions_mult[0]->out("  Response function: " + );
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("OPTIONS FOR STEPWISE PROCEDURE: \n");
  unsigned i;
  for(i=1;i<names_fixed.size();i++)
     {
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  OPTIONS FOR LINEAR EFFECTS TERM: "
        + names_fixed[i] + "\n");
     genoptions_mult[0]->out("\n");
          unsigned j;
     for(j=0;j<startfix.size();j++)
        if(startfix[j][i-1] == 0)
           genoptions_mult[0]->out("  Startvalue of the "
              + ST::doubletostring(j+1) + ". startmodel is \"effect excluded\" \n");
        else
           genoptions_mult[0]->out("  Startvalue of the "
              + ST::doubletostring(j+1) + ". startmodel is the fixed effect \n");
     }

  for(i=1;i<fullcondp.size();i++)
     {
     fullcondp[i]->set_inthemodel(1);
     genoptions_mult[0]->out("\n");
     genoptions_mult[0]->out("  OPTIONS FOR NONPARAMETRIC TERM: "
          + names_nonp[i-1][0] + "\n");
     genoptions_mult[0]->out("\n");
     if(fullcondp[i]->get_lambdamin()!=0 && fullcondp[i]->get_lambdamin()!=-1)
          {
          genoptions_mult[0]->out("  Minimum value for the smoothing parameter: "
             + ST::doubletostring(fullcondp[i]->get_lambdamin()) + "\n");
          fullcondp[i]->update_stepwise(fullcondp[i]->get_lambdamin());
          if(fullcondp[i]->get_df_equidist()==false || fullcondp[i]->get_spfromdf()=="direct")
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: "
               + ST::doubletostring(fullcondp[i]->compute_df(),6) + "\n");
          else
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: approximately "
               + ST::doubletostring(fullcondp[i]->get_df_lambdamin()) + ", exact "
               + ST::doubletostring(fullcondp[i]->compute_df(),6) + "\n");
          }
     if(fullcondp[i]->get_lambdamax()!=0 && fullcondp[i]->get_lambdamax()!=-1)
          {
          genoptions_mult[0]->out("  Maximum value for the smoothing parameter: "
             + ST::doubletostring(fullcondp[i]->get_lambdamax(),6) + "\n");
          fullcondp[i]->update_stepwise(fullcondp[i]->get_lambdamax());
          if(fullcondp[i]->get_df_equidist()==false || fullcondp[i]->get_spfromdf()=="direct")
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: "
               + ST::doubletostring(fullcondp[i]->compute_df(),6) + "\n");
          else
            genoptions_mult[0]->out("  This is equivalent to degrees of freedom: approximately "
               + ST::doubletostring(fullcondp[i]->get_df_lambdamax()) + ", exact "
               + ST::doubletostring(fullcondp[i]->compute_df(),6) + "\n");
          }
     if(fullcondp[i]->get_df_equidist()==true && fullcondp[i]->get_spfromdf()!="direct")
         genoptions_mult[0]->out("  Number of different smoothing parameters with equidistant degrees of freedom: "
            + ST::doubletostring(fullcondp[i]->get_number()) + "\n");
     else if(fullcondp[i]->get_fctype()!=MCMC::factor)
         genoptions_mult[0]->out("  Number of different smoothing parameters on a logarithmic scale: "
            + ST::doubletostring(fullcondp[i]->get_number()) + "\n");
     unsigned j;
     for(j=0;j<startindex.size();j++)
          {
          if(lambdavec[i-1][startindex[j][i-1]]==0)
             genoptions_mult[0]->out("  Startvalue of the "
                + ST::doubletostring(j+1) + ". startmodel is \"effect excluded\" \n");
          else if(lambdavec[i-1][startindex[j][i-1]]==-1)
             genoptions_mult[0]->out("  Startvalue of the "
                + ST::doubletostring(j+1) + ". startmodel is the fixed effect \n");
          else
             {
             genoptions_mult[0]->out("  Startvalue of the smoothing parameter for the "
                + ST::doubletostring(j+1) + ". startmodel: "
                + ST::doubletostring(lambdavec[i-1][startindex[j][i-1]],6) + "\n");
             fullcondp[i]->update_stepwise(lambdavec[i-1][startindex[j][i-1]]);
             genoptions_mult[0]->out("  This is equivalent to degrees of freedom: "
                + ST::doubletostring(fullcondp[i]->compute_df(),6) + "\n");
             }
          }
     fullcondp[i]->set_inthemodel(0);
     }
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("STEPWISE PROCEDURE STARTED \n");
  genoptions_mult[0]->out("\n");
  }


// -----------------------------------------------------------------------------
// ------- Funktionen für die Ausgabe im Tex-File ------------------------------
// -----------------------------------------------------------------------------

void STEPWISErun::make_graphics(const ST::string & name,
                  vector<vector<unsigned> > & startindex)
  {
  ST::string title = "STEPWISEREG OBJECT " + name + ": " + algorithm + " procedure";
    //erzeugt den Kopf des Tex-Files
  outtex << "\\documentclass[a4paper, 12pt]{article}" << endl
         << "\n" << "\\usepackage{graphicx}" << endl
         << "\\parindent0em" << endl
         << "\\textheight22cm \\textwidth15cm \\oddsidemargin0.5cm" << endl
         << "\n\\begin{document}" << endl
         << "\\begin{center}" << endl
         << "\\LARGE{\\bf " << title << "}"
         << endl << "\\end{center} \n\\vspace{1cm}" << endl;

  make_model();

  make_options();

  make_prior(startindex);

  outtex << "\n\\noindent {\\bf \\large Start Predictor";
  if(startindex.size()>1)
     outtex << "s";
  outtex << ":}\\\\" << endl;

  }


void STEPWISErun::make_tex_end(ST::string & path, vector<double> & modell,const ST::string & CI)
  {
  ST::string path_batch = path + "_graphics.prg";
  ST::string path_splus = path +  "_r.R";
  //ST::string path_stata = path +  "_stata.do";
  ST::string path_tex = path + "_model_summary.tex";

  outtex << "\n\\noindent {\\bf \\large Final Predictor:}\\\\" << endl;
  make_predictor();

  unsigned j;
  outtex << "\n\\noindent {\\bf \\large Final Properties:}\\\\ \n\\\\" << endl;
  for(j=1;j<fullcond_alle.size();j++)
     {
     if(modell[names_fixed.size()-2+j]!=0 && modell[names_fixed.size()-2+j]!=-1)
         {
         vector<ST::string> prior = fullcond_alle[j]->get_priorassumptions();
         outtex << prior[0] << "\\\\" << endl
                << "smoothing parameter: $\\lambda = "
                << ST::doubletostring(modell[names_fixed.size()-2+j],6)
                << " \\,\\, \\hat{=} \\,\\, \\mbox{df} = "
                << ST::doubletostring(fullcond_alle[j]->compute_df(),6)
                << "$ \\\\ \n\\\\" << endl;
         }
     }

  vector<ST::string> distr_results;
  distr_results = likep_mult[0]->get_results_latex();
  unsigned i;
  for (i=0;i<distr_results.size();i++)
     {
     outtex << distr_results[i] << endl;
     }

  make_fixed_table(CI);

    // Pfade der Files
    //werden im BayesX-Output angegeben
  genoptions_mult[0]->out("  Files of model summary: \n" , true);
  genoptions_mult[0]->out("\n");

  make_plots(path_batch,path_splus);

  genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  Latex file of model summary is stored in file \n");
  genoptions_mult[0]->out("  " + path_tex + "\n");
  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
  genoptions_mult[0]->out("\n");

  outtex << "\\end{document}" << endl;
  }


void STEPWISErun::make_options(void)
  {
  //double l1 = genoptions_mult[0]->get_level1();
  //double l2 = genoptions_mult[0]->get_level2();
  char hcharu = '_';
  ST::string hstringu = "\\_";

  //schreibt Stepwise options ins Tex-File
  outtex << "\n\\noindent {\\bf \\large Stepwise Options:}" << endl
         << "\\begin{tabbing}" << endl
         //<< "Levels for credible intervals: \\= \\\\" << endl
         //<< "Level 1: \\> " << l1 << "\\\\" << endl
         //<< "Level 2: \\> " << l2 << "\\\\" << endl
         << "Maximum number of Iterations: \\= " << steps << " \\\\" << endl
         << "Performance criterion: \\> " << criterion.insert_string_char(hcharu,hstringu) << " \\\\" << endl
         << "Startmodel: \\> " << startmodel << " \\\\" << endl
         //<< "Used number of Iterations: \\> blabla \\\\" << endl
         << "Increment: \\> " << increment << " \\\\" << endl;
  outtex << "\\end{tabbing}\n"  << "\\vspace{0.5cm}" <<  endl;
  }


void STEPWISErun::make_predictor(void)
  {
  unsigned i;
  ST::string term = "$\\eta$ & $=$ & $\\gamma_0";

  char hcharu = '_';
  ST::string hstringu = "\\_";
  for(i=1;i<fullcondp[0]->get_datanames().size();i++)
     term = term + " + " + fullcondp[0]->get_datanames()[i].insert_string_char(hcharu,hstringu);
  for(i=1;i<fullcondp.size();i++)
     term = term + " + " + fullcondp[i]->get_term_symbolic();

  outtex << endl << "\n\\begin{tabular}{ccp{12cm}}\n" << term
         << "$\n\\end{tabular}\n\\\\ \n\\\\" << endl;
  outtex << criterion.insert_string_char(hcharu,hstringu) << " = " << ST::doubletostring(kriterium_tex,6) << " \\\\ \n\\\\" << endl;
  }


void STEPWISErun::make_model(void)
  {
  //Vert-Fam wird übergeben
  ST::string fam = likep_mult[0]->get_family();
  fam = fam.replaceallsigns('_', ' ');

  //Anz. Beob. wird übergeben
  unsigned obs = likep_mult[0]->get_nrobs();

  //Name der Resp.-Var. übergeben
  ST::string resp = likep_mult[0]->get_responsename();
  char hcharu = '_';
  ST::string hstringu = "\\_";
  resp = resp.insert_string_char(hcharu,hstringu);

  //schreibt das Modell und die Prior-annahmen ins Tex-File
  outtex << "\n\\noindent {\\bf \\large Response:}" << endl
         << "\\begin{tabbing}\n"
         << "Number of observations: \\= " << obs << "\\\\" << endl
         << "Response Variable: \\> " << resp << "\\\\" << endl
         << "Family: \\> " << fam << "\\\\" << endl
         << "\\end{tabbing} \n" << endl;
  }


void STEPWISErun::make_prior(vector<vector<unsigned> > & startindex)
  {
  vector<ST::string> names_fixed = fullcond_alle[0]->get_datanames();
  outtex << "\n\\noindent {\\bf \\large Assumptions:}\\\\" << endl << "\\\\" << endl;
  unsigned i;
  if(names_fixed.size()>1)
     {
     outtex << "Linear Effects:\\\\";
     ST::string term = "";
     for(i=1;i<names_fixed.size()-1;i++)
        term = term + "$" + names_fixed[i] + "$, ";
     term = term +  "$" + names_fixed[names_fixed.size()-1];
     outtex << endl << "\\begin{tabular}{p{12cm}}\n" << term
            << "$\n\\end{tabular}\n" << endl;
     if(startmodel == "full" || startmodel == "userdefined")
        outtex << "Startvalue is a linear effect \\\\ \n\\\\" << endl;
     else if(startmodel == "empty")
        outtex << "Startvalue is 'effect excluded' \\\\ \n\\\\" << endl;
     else
        outtex << "1. Startvalue is 'effect excluded' \\\\" << endl
               << "2. Startvalue is a linear effect \\\\ \n\\\\" << endl;
    }

  unsigned j;
  for(i=1;i<fullcond_alle.size();i++)
     {
     j = 0;
     fullcond_alle[i]->set_inthemodel(1);
     vector<ST::string> prior = fullcond_alle[i]->get_priorassumptions();
     for(j=0;j<prior.size()-1;j++)
        outtex << prior[j] << "\\\\" << endl;
     if(fullcond_alle[i]->get_lambdamin()!=0 && fullcond_alle[i]->get_lambdamin()!=-1)
        {
        outtex << "Minimum value for the smoothing parameter: $\\lambda = "
               << ST::doubletostring(fullcond_alle[i]->get_lambdamin())
               << " \\,\\, \\hat{=} \\,\\, \\mbox{df} = ";
        fullcond_alle[i]->update_stepwise(fullcond_alle[i]->get_lambdamin());
        if(fullcond_alle[i]->get_spfromdf()=="direct")
           outtex << ST::doubletostring(fullcond_alle[i]->compute_df(),6);
        else
           outtex << ST::doubletostring(fullcond_alle[i]->compute_df(),6)
                  << " \\approx " << ST::doubletostring(fullcond_alle[i]->get_df_lambdamin());
        outtex << "$ \\\\ \n";

        outtex << "Maximum value for the smoothing parameter: $\\lambda = "
               << ST::doubletostring(fullcondp[i]->get_lambdamax(),6)
               << " \\,\\, \\hat{=} \\,\\, \\mbox{df} = ";
        fullcond_alle[i]->update_stepwise(fullcond_alle[i]->get_lambdamax());
        if(fullcond_alle[i]->get_spfromdf()=="direct")
           outtex << ST::doubletostring(fullcond_alle[i]->compute_df(),6);
        else
           outtex << ST::doubletostring(fullcond_alle[i]->compute_df(),6)
                  << " \\approx " << ST::doubletostring(fullcond_alle[i]->get_df_lambdamax());
        outtex << "$ \\\\ \n";

        if(fullcondp[i]->get_df_equidist()==true && fullcondp[i]->get_spfromdf()!="direct")
            outtex << "Number of different smoothing parameters with equidistant degrees of freedom: "
                   << ST::doubletostring(fullcond_alle[i]->get_number()) << " \\\\ \n";
        else
            outtex << "Number of different smoothing parameters on a logarithmic scale: "
                   << ST::doubletostring(fullcondp[i]->get_number()) << " \\\\ \n";
        }

     if(fullcond_alle[i]->get_forced()==true)
         outtex << "Without the excluded effect" << " \\\\ \n";

     j = 0;
     for(j=0;j<startindex.size();j++)
        {
        if(lambdavec[i-1][startindex[j][i-1]]==0)
           {
           if(startindex.size()>1)
              outtex << ST::doubletostring(j+1) << ". ";
           outtex << "Startvalue is 'effect excluded' \\\\ \n";
           }
        else if(lambdavec[i-1][startindex[j][i-1]]==-1)
           {
           if(startindex.size()>1)
              outtex << ST::doubletostring(j+1) << ". ";
           outtex << "Startvalue is the fixed effect \\\\ \n";
           }
        else
           {
           if(startindex.size()>1)
              outtex << ST::doubletostring(j+1) << ". ";
           outtex << "Startvalue of the smoothing parameter: $\\lambda = "
                  << ST::doubletostring(lambdavec[i-1][startindex[j][i-1]],6)
                  << " \\,\\, \\hat{=} \\,\\, \\mbox{df} = ";
           fullcond_alle[i]->update_stepwise(lambdavec[i-1][startindex[j][i-1]]);
           outtex << ST::doubletostring(fullcond_alle[i]->compute_df(),6) << "$ \\\\ \n";
           }
        }
     outtex << "\\\\" << endl;
     fullcond_alle[i]->set_inthemodel(0);
     }
  }


void STEPWISErun::make_fixed_table(const ST::string & CI)
  {

  // falls andere Quantile gewünscht werden
  double u = fullcondp[begin[0]]->get_level1();
  double o = fullcondp[begin[0]]->get_level2();
  double u1 = fullcondp[begin[0]]->get_lower1();
  double u2 = fullcondp[begin[0]]->get_upper2();
  double o1 = fullcondp[begin[0]]->get_lower2();
  double o2 = fullcondp[begin[0]]->get_upper1();
  ST::string u_str = ST::doubletostring(u,0);
  ST::string o_str = ST::doubletostring(o,0);
  ST::string u1_str = ST::doubletostring(u1,5);
  ST::string u2_str = ST::doubletostring(u2,5);
  ST::string o1_str = ST::doubletostring(o1,5);
  ST::string o2_str = ST::doubletostring(o2,5);

  vector<ST::string> h;
  unsigned j;
  unsigned r = 2;
  outtex << "\n\\newpage \n" << endl << "\n\\noindent {\\bf \\large Linear Effects:}\\\\"
           << endl << "\\\\" << endl;



  outtex << "\\begin{tabular}{|r|r|r|r|r|r|}" << endl << "\\hline" << endl
           << "Variable & Mean & Std & CI" << ST::doubletostring(u)
           << "lower & Median & CI" << ST::doubletostring(u)
           << "upper \\\\" << endl << "\\hline" << endl;

  h = fullcondp[0]->get_results_latex();
  for(j=0;j<h.size();j++)
     {
     r++;
     if (r < 39)
 //       outtex << h[j].substr(0,h[j].length()-17) << "\\\\" << endl;
      outtex << h[j]  << endl;
     else
        {
        r=1;
        outtex << "\\hline \n\\end{tabular}" << endl;
        outtex << "\n\\newpage \n" << endl
               << "\n\\noindent {\\bf \\large Linear Effects (continued):}\\\\"
               << endl << "\\\\" << endl;
        outtex << "\\begin{tabular}{|r|r|}" << endl << "\\hline" << endl
               << "Variable & Mean\\\\" << endl << "\\hline" << endl;
        outtex << h[j] << endl;
        }
     }
  outtex << "\\hline \n\\end{tabular}" << endl;
  }


void STEPWISErun::make_plots(ST::string & path_batch,
                             ST::string & path_splus)      //,ST::string & path_stata)
  {

  char hcharu = '_';
  ST::string hstringu = "\\_";

  unsigned j;

  ST::string pathresult;

  bool stil = false;

  // Schleife überprüft, ob es ein fullcond-Object
  // gibt, bei dem Effekt gezeichnet werden kann
  MCMC::plotstyles plst;
  for(j=0;j<fullcondp.size();j++)
    {
    plst = fullcondp[j]->get_plotstyle();
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
/*    outsplus << "# NOTE: 'directory' must be substituted by the directory"
             << " where the sfunctions are stored \n"
             << endl
    // einlesen der Source-Files für S-Plus
             << "source(\"'directory'\\\\sfunctions\\\\plotsample.s\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\plotnonp.s\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\plotsurf.s\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\drawmap.s\")" << endl
             << "source(\"'directory'\\\\sfunctions\\\\readbndfile.s\")\n" << endl;*/

    genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  Batch file for visualizing effects of nonlinear functions is stored in file \n");
    genoptions_mult[0]->out("  " + path_batch + "\n");
    genoptions_mult[0]->out("\n");

    bool stil2 = true;
    for(j=begin[0];j<=end[0];j++)  //Schleife überprüft, ob es map-Objekt gibt
      {
      plst = fullcondp[j]->get_plotstyle();
      if(plst == MCMC::drawmap || plst == MCMC::drawmapgraph)
        stil2 = false;
      }

    if(stil2 == true)
      {
      genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  Batch file for visualizing effects of nonlinear functions ");
      genoptions_mult[0]->out("  in R is stored in file \n");
      genoptions_mult[0]->out("  " + path_splus + "\n");
      genoptions_mult[0]->out("\n");
      }

    if(stil2 == false)
      {
      genoptions_mult[0]->out("  NOTE: 'input filename' must be substituted by the filename of the boundary-file \n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  Batch file for visualizing effects of nonlinear functions ");
      genoptions_mult[0]->out("  in R is stored in file \n");
      genoptions_mult[0]->out("  " + path_splus + "\n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  NOTE: 'input filename' must be substituted by the filename of the boundary-file \n");
      genoptions_mult[0]->out("\n");
      }

    outbatch << "% usefile " << path_batch << endl;

    // falls andere Quantile gewünscht werden
    double u = fullcondp[begin[0]]->get_level1();
    double o = fullcondp[begin[0]]->get_level2();
    double u1 = fullcondp[begin[0]]->get_lower1();
    double u2 = fullcondp[begin[0]]->get_upper2();
    double o1 = fullcondp[begin[0]]->get_lower2();
    double o2 = fullcondp[begin[0]]->get_upper1();
    ST::string u_str = ST::doubletostring(u,0);
    ST::string o_str = ST::doubletostring(o,0);
    ST::string u1_str = ST::doubletostring(u1,5);
    ST::string u2_str = ST::doubletostring(u2,5);
    ST::string o1_str = ST::doubletostring(o1,5);
    ST::string o2_str = ST::doubletostring(o2,5);

    // durchlaufen der Fullconditionals
    for(j=0;j<fullcondp.size();j++)
      {
      // Pfad der Regr.-Ergebnisse
      pathresult = fullcondp[j]->get_pathresult();

      // Plotstyle: noplot, plotnonp, drawmap, drawmapgraph
      plst = fullcondp[j]->get_plotstyle();

      if (plst != MCMC::noplot)
        {
        // Pfade für ps-, tex-, SPlus-files
        ST::string pathps = pathresult.substr(0, pathresult.length()-4);
        ST::string pathgr = pathps.replaceallsigns('\\', '/');

        char hchar = '\\';
        ST::string hstring = "/";

        ST::string pathps_spl = pathps.insert_string_char(hchar,hstring);
        ST::string pathres_spl = pathresult.insert_string_char(hchar,hstring);
        if(plst == MCMC::plotnonp)
           {
           outbatch << "\n";                // Befehle f. d. batch-file
           outbatch << "dataset _dat" << endl;
           outbatch << "_dat.infile using " << pathresult << endl;
           outbatch << "graph _g" << endl;
           vector<ST::string> varnames = fullcondp[j]->get_datanames();
           ST::string xvar = varnames[0];
           outbatch << "_g.plot " << xvar
                    << " pmean pqu" << u1_str.replaceallsigns('.','p')
                    << " pqu"
                    << o1_str.replaceallsigns('.','p') << " pqu"
                    << o2_str.replaceallsigns('.','p') << " pqu"
                    << u2_str.replaceallsigns('.','p')
                    << ", " << "title = \"Effect of " << xvar << "\" xlab = " << xvar
                    << " ylab = \" \" " << "outfile = " << pathps
                    << ".ps replace using _dat" << endl;
           outbatch << "drop _dat" << endl;
           outbatch << "drop _g" << endl;
            // Plot-Befehle f. d. SPlus-file
           outsplus << "plotnonp(\"" << pathres_spl << "\")" << endl;
            // Plot-Befehle f. d. tex-file
           outtex << "\n\\begin{figure}[h!]" << endl
                  << "\\centering" << endl
                  << "\\includegraphics[scale=0.6]{" << pathgr << "}" << endl
                  << "\\caption{Non--linear Effect of '"
                  << xvar.insert_string_char(hcharu,hstringu) << "'";
           outtex << "." << endl << "Shown are the posterior means.}" << endl
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
           vector<ST::string> varnames = fullcondp[j]->get_datanames();
           ST::string regionvar = varnames[0];
           outbatch << "_g.drawmap " << "pmean" << " " << regionvar
                    << ", map = _map color outfile = " << pathps
                    << "_pmean.ps replace using _dat" << endl;
           /*
           outbatch << "_g.drawmap " << "pcat" << u_str << " " << regionvar
                      << ", map = _map nolegend pcat outfile = " << pathps
                      << "_pcat" << u_str << ".ps replace using _dat" << endl;
           outbatch << "_g.drawmap " << "pcat" << o_str << " " << regionvar
                      << ", map = _map nolegend pcat outfile = " << pathps
                      << "_pcat" << o_str << ".ps replace using _dat" << endl;
           */
           outbatch << "drop _dat" << endl;
           outbatch << "drop _g" << endl;
           outbatch << "drop _map" << endl;
            // Plot-Befehle f. d. SPlus-file
           outsplus << "# NOTE: 'input_filename' must be substituted by the "
                    << "filename of the boundary-file \n"
//                    << "# NOTE: choose a 'name' for the map \n" << endl
                     << "m <- read.bnd(\"'input_filename'\")" << endl
                    << "drawmap(data = \"" << pathres_spl << "\", map = m, plotvar = \"pmean\", regionvar = \""
                    << regionvar << "\")" << endl;
           /*
           outsplus << "drawmap(map = 'name', outfile = \"" << pathps_spl
                    <<"_pcat" << u_str << ".ps\", dfile = \"" << pathres_spl
                    << "\" ,plotvar = \"pcat" << u_str << "\", regionvar = \""
                    << regionvar << "\", legend = F, pcat = T)" << endl;
           outsplus << "drawmap(map = 'name', outfile = \"" << pathps_spl
                    <<"_pcat" << o_str << ".ps\", dfile = \"" << pathres_spl
                    << "\",plotvar = \"pcat" << o_str << "\", regionvar = \""
                    << regionvar << "\", legend = F, pcat = T)" << endl;
           */
            // Plot-Befehle f. d. tex-file
           if(plst == MCMC::drawmap)
             {
             outtex << "\n\\begin{figure}[h!]" << endl
                    << "\\centering" << endl
                    << "\\includegraphics[scale=0.6]{" << pathgr << "_pmean}"
                    << endl
                    << "\\caption{Non--linear Effect of '" <<
                    regionvar.insert_string_char(hcharu,hstringu) << "'";
             outtex << ". Shown are the posterior means.}" << endl
                    << "\\end{figure}" << endl;
             }
           else if(plst == MCMC::drawmapgraph)
             {
             outtex << "\n%\\begin{figure}[h!]" << endl
                    << "%\\centering" << endl
                    << "%\\includegraphics[scale=0.6]{" << pathgr << "_pmean}"
                    << endl
                    << "%\\caption{Non--linear Effect of '" <<
                    regionvar.insert_string_char(hcharu,hstringu) << "'";
             outtex << ". Shown are the posterior means.}" << endl
                    << "%\\end{figure}" << endl;
             }
           /*
           outtex << "\n\\begin{figure}[htb]" << endl
                  << "\\centering" << endl
                  << "\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                  << u_str << "}" << endl
                  << "\\caption{Non--linear Effect of '" << regionvar << "'";
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
                  << "\\caption{Non--linear Effect of '" << regionvar << "'";
           outtex << ". Posterior probabilities for a nominal level of "
                  << o_str << "\\%." << endl
                  << "Black denotes regions with strictly negative credible intervals,"
                  << endl
                  << "white denotes regions with strictly positive credible intervals.}"
                  << endl << "\\end{figure}" << endl;
           */

           } // end: else if (plst == MCMC::drawmap && fullcondp[j]->get_col()==i)
        } // end: if (plst != MCMC::noplot)
      } // end: for(j=begin[0];j<=end[nr];j++)
    }
  }


// -----------------------------------------------------------------------------
// ------------- Model Averaging -----------------------------------------------
// -----------------------------------------------------------------------------

bool STEPWISErun::confidence_intervals(const ST::string & CI, const vector<double> & modell_final,
              const double & kriterium_final, vector<FULLCOND*> & fullcond_z)
  {
  bool abbruch=false;
  if(CI=="MCMCselect")
    {
    abbruch = confidence_MCMCselect(modell_final,kriterium_final,fullcond_z);
    }
  else if(CI=="MCMCvariable")
    {
    }
  else if(CI=="MCMCbootstrap")
    {
    abbruch = confidence_MCMCbootstrap(modell_final,kriterium_final,fullcond_z);
    }
  else if(CI=="bootstrap")
    {
    abbruch = confidence_bootstrap(modell_final,kriterium_final,fullcond_z);
    }

  return abbruch;
  }


bool STEPWISErun::confidence_MCMCselect(const vector<double> & modell_final,
              const double & kriterium_final, vector<FULLCOND*> & fullcond_z)
  {
  bool abbruch;
  unsigned i;
  fullcond_z = fullcondp;
  for(i=0;i<fullcond_z.size();i++)
    fullcond_z[i]->set_fcnumber(i);

  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("CALCULATION OF CONFIDENCE INTERVALS STARTED:\n",true);
  genoptions_mult[0]->out("\n");

  if(likep_mult[0]->get_scaleexisting()==true)
    {
    double scale = likep_mult[0]->get_scale();
    double transform = likep_mult[0]->get_trmult(0);
    scale *= transform*transform;
    likep_mult[0]->set_constscale(scale);
    }

  //bool neuschaetzen = false;
  for(i=1;i<fullcond_alle.size();i++)
    {
    if(modell_final[names_fixed.size()-2+i] == -2)
      {
      fullcond_alle[i]->change_Korder(modell_final[names_fixed.size()-2+i]);
      //neuschaetzen = true;
      }
    else
      fullcond_alle[i]->set_lambdaconst(modell_final[names_fixed.size()-2+i]);
    }

  //if(neuschaetzen == true)
    schaetzen(0,kriterium_alt,true,"backfitting");  // nötig für Startwerte, wenn mind. ein lambda=-2

  unsigned iterations = genoptions_mult[0]->get_iterations();
  unsigned start = 1;
  int seed = likep_mult[0]->get_seed();
  abbruch = simulate(posttitle,seed,start,iterations);
  if(abbruch == true)
    return true;

  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("ESTIMATION RESULTS:\n",true);
  genoptions_mult[0]->out("\n");
  likep_mult[0]->outresults();
  for(unsigned f=0;f<fullcondp.size();f++)
    {
    fullcondp[f]->outresults();
    }
  return abbruch;
  }


bool STEPWISErun::confidence_bootstrap(const vector<double> & modell_final,
              const double & kriterium_final, vector<FULLCOND*> & fullcond_z)
  {
  unsigned samplesize_df = bootstrap + 1;
  unsigned zaehler = 1;
  bool abbruch;
  if(unconditional == true)
   {
   steps = 0;
   }

  isboot = true;
  trace = "trace_off";
  vector<double> modell_boot = modell_final;
  double kriterium_boot = kriterium_final;
  modell_alt = modell_final;
  kriterium_alt = kriterium_final;
  fix_ganz_komplett(modell_alt);
  fullcond_komplett(modell_alt);
  schaetzen(0,kriterium_alt,true,"backfitting");
  update_bootstrap(zaehler);

  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("BEGINNING OF BOOTSTRAP:\n",true);
  genoptions_mult[0]->out("\n");

  while(bootstrap > 0)
    {
    zaehler += 1;
    bootstrap -= 1;
    minim = minim2;
    fertig = false;

    likep_mult[0]->create_bootstrap_weights();

    modellematrix.erase(modellematrix.begin(),modellematrix.end());
    vector<vector<double> > startiteration;
    startiteration.push_back(modell_boot);
    modellematrix.push_back(startiteration);
    fix_ganz_komplett(modell_boot);
    fullcond_komplett(modell_boot);
    modell_alt = modell_boot;
    modell_neu = modell_boot;

    if(criterion == "MSEP" || criterion == "AUC")
      {
      likep_mult[0]->weight_for_all();   // macht hier das Gegenteil und setzt "weight" für MSEP ein!
      for(unsigned y=0;y<fullcond_alle.size();y++)
        fullcond_alle[y]->set_calculate_xwx();
      }
    schaetzen(0,kriterium_alt,true,"backfitting");
    kriterium_neu = kriterium_alt;
    steps_aktuell = 0;

    if(algorithm != "coorddescent")
      abbruch = stepfunctions();
    else
      abbruch = koordabstieg();
    if(abbruch==true)
      return true;

    fix_ganz_komplett(modell_alt);
    fullcond_komplett(modell_alt);
    if(criterion == "MSEP" || criterion == "AUC")
      {
      likep_mult[0]->weight_for_all();
      for(unsigned y=0;y<fullcond_alle.size();y++)
        fullcond_alle[y]->set_calculate_xwx();
      }
    schaetzen(0,kriterium_alt,true,"backfitting");
    update_bootstrap(zaehler);
    }

  modell_alt = modell_boot;
  kriterium_alt = kriterium_boot;

   genoptions_mult[0]->out("\n");
   genoptions_mult[0]->out("ESTIMATION RESULTS:\n",true);
   genoptions_mult[0]->out("\n");

   likep_mult[0]->set_original_response();
   likep_mult[0]->update_bootstrap_betamean();
   likep_mult[0]->outresults();

  if(unconditional)
    {
    fullcond_z = fullcondp;
    for(unsigned f=0;f<fullcondp.size();f++)
      {
      fullcondp[f]->set_fcnumber(f);
      fullcondp[f]->update_bootstrap_betamean();
      fullcondp[f]->outresults_df(samplesize_df);
      fullcondp[f]->outresults();
      }
    }
  else
    {
    fullcond_z = fullcond_alle;
    fullcondp = fullcond_alle;
    for(unsigned f=0;f<fullcond_alle.size();f++)
      {
      fullcond_alle[f]->update_bootstrap_betamean();
      fullcond_alle[f]->outresults_df(samplesize_df);
      fullcond_alle[f]->outresults();
      }
    }
  return abbruch;
  }


bool STEPWISErun::confidence_MCMCbootstrap(const vector<double> & modell_final,
              const double & kriterium_final, vector<FULLCOND*> & fullcond_z)
  {
  unsigned samplesize_df = bootstrap + 1;
  unsigned zaehler = 1;
  bool abbruch;
  unsigned i;
  isboot = true;
  trace = "trace_off";
  vector<double> modell_boot = modell_final;
  double kriterium_boot = kriterium_final;
  modell_alt = modell_final;
  kriterium_alt = kriterium_final;
  fix_ganz_komplett(modell_alt);
  fullcond_komplett(modell_alt);
  schaetzen(0,kriterium_alt,true,"backfitting");
  likep_mult[0]->save_betamean();
  for(i=1;i<fullcond_alle.size();i++)
    {
    fullcond_alle[i]->update_bootstrap_df();
    fullcond_alle[i]->save_betamean();
    }
  fullcond_alle[0]->update_bootstrap_df();
  fullcond_alle[0]->save_betamean();

  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("CALCULATION OF CONFIDENCE BANDS STARTED:\n",true);
  genoptions_mult[0]->out("\n");

  // MCMC-Teil
  int seed = likep_mult[0]->get_seed();
  if(likep_mult[0]->get_scaleexisting()==true)
    {
    double scale = likep_mult[0]->get_scale();
    double transform = likep_mult[0]->get_trmult(0);
    scale *= transform*transform;
    likep_mult[0]->set_constscale(scale);
    }

  vector<double> modell_mcmc = modell_alt;
  for(i=1;i<fullcond_alle.size();i++)
    {
    if((modell_alt[names_fixed.size()-2+i] == -1 || modell_alt[names_fixed.size()-2+i] == -2)
                      && fullcond_alle[i]->get_fctype() != MCMC::factor)
      {
      fullcond_alle[i]->change_Korder(modell_alt[names_fixed.size()-2+i]);
      modell_mcmc[names_fixed.size()-2+i] = 1000000000;
      }
    else if(modell_alt[names_fixed.size()-2+i] == 0 && fullcond_alle[i]->get_fctype() != MCMC::factor)
      modell_mcmc[names_fixed.size()-2+i] = 1000000000;
    else
      fullcond_alle[i]->set_lambdaconst(modell_alt[names_fixed.size()-2+i]);
    }

  fix_ganz_komplett(modell_mcmc);        // stellt "fullcondp" neu auf (für Funktionen mit lambda=-1 -> lambda=10^8)
  fullcond_komplett(modell_mcmc);        // fullcond-Objekte für leere Funktionen sind auch enthalten.

  for(i=1;i<fullcond_alle.size();i++)
    {
    if(modell_alt[names_fixed.size()-2+i] == 0 && fullcond_alle[i]->get_fctype() != MCMC::factor)
      {
      fullcond_alle[i]->set_lambdaconst(modell_alt[names_fixed.size()-2+i]);
      fullcond_alle[i]->set_inthemodel(0);
      }
    else if((modell_alt[names_fixed.size()-2+i] == -1 || modell_alt[names_fixed.size()-2+i] == -2)
                             && fullcond_alle[i]->get_fctype() != MCMC::factor)
        {
        fullcond_alle[i]->set_inthemodel(-1);
        }
    }

  unsigned iterations = genoptions_mult[0]->get_iterations();
  unsigned startit = 1;
  unsigned endit = double(iterations) / double(bootstrap+1);
  unsigned itpmod = endit;
  schaetzen(0,kriterium_alt,true,"backfitting");  // nötig, da fixe Effekte in "fix_ganz_komplett" gleich Null gesetzt werden.
  for(i=1;i<fullcond_alle.size();i++)             // speichert MA-Schätzer für Funktionen
    fullcond_alle[i]->update_beta_average(zaehler);
  fullcond_alle[0]->update_beta_average(zaehler);

  abbruch = simulate(posttitle,seed,startit,endit);
  fullcond_alle[0]->update_linold();
  fullcond_alle[0]->update_linold_vc();
  if(abbruch==true)
    return true;

  for(i=1;i<fullcond_alle.size();i++)
    {
    if((modell_final[names_fixed.size()-2+i] == -1 || modell_final[names_fixed.size()-2+i] == -2)
                   && fullcond_alle[i]->get_fctype() != MCMC::factor)
      {
      fullcond_alle[i]->undo_Korder();
      fullcond_alle[i]->set_lambdaconst(modell_final[names_fixed.size()-2+i]);
      }
    }

  endit += itpmod;
  startit += itpmod;

  while(bootstrap > 0)
    {
    genoptions_mult[0]->out("  BOOTSTRAPSAMPLE " + ST::inttostring(bootstrap) + "\n",false);

    zaehler += 1;
    bootstrap -= 1;
    minim = minim2;
    fertig = false;

    // neues Modell selektieren
    likep_mult[0]->create_bootstrap_weights();

    modellematrix.erase(modellematrix.begin(),modellematrix.end());
    vector<vector<double> > startiteration;
    startiteration.push_back(modell_boot);
    modellematrix.push_back(startiteration);
    fix_ganz_komplett(modell_boot);
    fullcond_komplett(modell_boot);
    modell_alt = modell_boot;
    modell_neu = modell_boot;
    likep_mult[0]->undo_constscale();

    if(criterion == "MSEP" || criterion == "AUC")
      {
      likep_mult[0]->weight_for_all();   // macht hier das Gegenteil und setzt "weight" für MSEP ein!
      for(unsigned y=0;y<fullcond_alle.size();y++)
        fullcond_alle[y]->set_calculate_xwx();
      }
    schaetzen(0,kriterium_alt,true,"backfitting");
    kriterium_neu = kriterium_alt;
    steps_aktuell = 0;

    if(algorithm != "coorddescent")
      abbruch = stepfunctions();
    else
      abbruch = koordabstieg();
    if(abbruch==true)
      return true;

    fix_ganz_komplett(modell_alt);
    fullcond_komplett(modell_alt);

/*for(unsigned c=0;c<modell_alt.size();c++)
  {
  genoptions_mult[0]->out(ST::doubletostring(modell_alt[c],6) + "   ");
  }
genoptions_mult[0]->out("\n");*/

    likep_mult[0]->set_original_response();
    if(criterion == "MSEP" || criterion == "AUC")
      {
      likep_mult[0]->weight_for_all();
      for(unsigned y=0;y<fullcond_alle.size();y++)
        fullcond_alle[y]->set_calculate_xwx();
      }
    schaetzen(0,kriterium_alt,true,"backfitting");
    for(i=0;i<fullcond_alle.size();i++)
      {
      fullcond_alle[i]->update_bootstrap_df();
      }

  // MCMC-Teil
    seed += 1;
    if(likep_mult[0]->get_scaleexisting()==true)
      {
      double scale = likep_mult[0]->get_scale();
      double transform = likep_mult[0]->get_trmult(0);
      scale *= transform*transform;
      likep_mult[0]->set_constscale(scale);
      }

    vector<double> modell_mcmc = modell_alt;

    for(i=1;i<fullcond_alle.size();i++)
      {
      if((modell_alt[names_fixed.size()-2+i] == -1 || modell_alt[names_fixed.size()-2+i] == -2)
                             && fullcond_alle[i]->get_fctype() != MCMC::factor)
        {
        fullcond_alle[i]->change_Korder(modell_alt[names_fixed.size()-2+i]);
        modell_mcmc[names_fixed.size()-2+i] = 1000000000;
        }
      else if(modell_alt[names_fixed.size()-2+i] == 0 && fullcond_alle[i]->get_fctype() != MCMC::factor)
        modell_mcmc[names_fixed.size()-2+i] = 1000000000;
      else
        fullcond_alle[i]->set_lambdaconst(modell_alt[names_fixed.size()-2+i]);
      }

    fix_ganz_komplett(modell_mcmc);
    fullcond_komplett(modell_mcmc);

    for(i=1;i<fullcond_alle.size();i++)
      {
      if(modell_alt[names_fixed.size()-2+i] == 0 && fullcond_alle[i]->get_fctype() != MCMC::factor)
        {
        fullcond_alle[i]->set_lambdaconst(modell_alt[names_fixed.size()-2+i]);
        fullcond_alle[i]->set_inthemodel(0);
        }
      else if((modell_alt[names_fixed.size()-2+i] == -1 || modell_alt[names_fixed.size()-2+i] == -2)
                             && fullcond_alle[i]->get_fctype() != MCMC::factor)
        {
        fullcond_alle[i]->set_inthemodel(-1);
        }
      }

    // neu schätzen, da fixe Effekte in "fix_ganz_komplett" gleich Null gesetzt werden.
    schaetzen(0,kriterium_alt,true,"backfitting");
    for(i=1;i<fullcond_alle.size();i++)   // speichert MA-Schätzer für Funktionen
      fullcond_alle[i]->update_beta_average(zaehler);
    fullcond_alle[0]->update_beta_average(zaehler);

    // MCMC-Simulation
    abbruch = simulate(posttitle,seed,startit,endit);
    if(abbruch==true)
      return true;
    // rechnet die Zentrierungskonstanten zum linearen Prädiktor
    fullcond_alle[0]->update_linold();
    fullcond_alle[0]->update_linold_vc();

    for(i=1;i<fullcond_alle.size();i++)
      {
      if((modell_alt[names_fixed.size()-2+i] == -1 || modell_alt[names_fixed.size()-2+i] == -2)
                   && fullcond_alle[i]->get_fctype() != MCMC::factor)
        {
        fullcond_alle[i]->undo_Korder();
        fullcond_alle[i]->set_lambdaconst(modell_alt[names_fixed.size()-2+i]);
        }
      }

    endit += itpmod;
    startit += itpmod;
    }

  modell_alt = modell_boot;
  kriterium_alt = kriterium_boot;
  fix_ganz_komplett(modell_alt);
  fullcond_komplett(modell_alt);
  likep_mult[0]->set_original_response();
  //schaetzen(0,kriterium_alt,true,"backfitting");

  fullcondp = fullcond_alle;
  for(i=1;i<fullcond_alle.size();i++)
    {
    if(fullcond_alle[i]->get_fctype() == MCMC::factor)
      fullcondp.erase(fullcondp.begin()+1,fullcondp.begin()+2);
    }
  fullcond_z = fullcondp;
  for(i=0;i<fullcond_z.size();i++)
    fullcond_z[i]->set_fcnumber(i);

  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("ESTIMATION RESULTS:\n",true);
  genoptions_mult[0]->out("\n");
  likep_mult[0]->update_bootstrap_betamean();
  likep_mult[0]->outresults();
  for(unsigned f=0;f<fullcond_alle.size();f++)
    {
    fullcond_alle[f]->update_bootstrap_betamean();
    fullcond_alle[f]->outresults_df(samplesize_df);
    fullcond_alle[f]->outresults();
    }

  return abbruch;
  }


void STEPWISErun::update_bootstrap(unsigned & zaehler)
  {
  unsigned i;
  genoptions_mult[0]->update_bootstrap();

  if(likepexisting)
    likep_mult[0]->update_bootstrap();

  if(unconditional)
    {
    for(i=1;i<fullcondp.size();i++)
      fullcondp[i]->update_bootstrap();
    }
  else
    {
    for(i=1;i<fullcond_alle.size();i++)
      {
      fullcond_alle[i]->update_bootstrap();
      fullcond_alle[i]->update_beta_average(zaehler);
      }
    fullcond_alle[0]->update_beta_average(zaehler);
    }
  fullcond_alle[0]->update_bootstrap(unconditional);

  if(likepexisting)
    likep_mult[0]->update_predict_bootstrap(bootstrap);    //likep_mult[0]->update_predict_bootstrap();
  }


bool STEPWISErun::simulate(const vector<ST::string> & header, const int & seed,
                           const unsigned & startit, const unsigned & endit)
  {
  unsigned i,j;
  unsigned nrmodels = genoptions_mult.size();
  bool errors=false;
  i=0;
  while( (i<nrmodels) && (errors==false) )
    {
    errors = checkerrors(likep_mult[i],fullcondp,begin[i],end[i]);
    i++;
    }

  if (errors==false)
  {
  unsigned it;
  //unsigned iterations = genoptions_mult[0]->get_iterations();

  #if defined(MICROSOFT_VISUAL)
    {

    }
  #elif!defined(__BUILDING_GNU)
    {
    srand((unsigned)time(0));
    }
  #else
    {
    srand((unsigned)time(0));
    }
  #endif

  if(seed >= 0)
    srand(seed);

  for (it=startit;it<=endit;it++)
    {

    for(i=0;i<nrmodels;i++)
      {
      genoptions_mult[nrmodels-1-i]->update();

      if (likepexisting)
        likep_mult[nrmodels-1-i]->update();        // DO NOT CHANGE ORDER !!!!

      for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
        {
        fullcondp[j]->update();
        } // end: for(j=0;j<nrfullcond;j++)

      if (likepexisting)
        likep_mult[nrmodels-1-i]->update_predict();
      }


#if defined(BORLAND_OUTPUT_WINDOW)
    Application->ProcessMessages();

    if (Frame->stop)
      {
//      genoptions->out("USER BREAK\n");
      break;
      }

    if (Frame->pause)
      {
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("SIMULATION PAUSED\n");
      genoptions_mult[0]->out("Click CONTINUE to proceed\n");
      genoptions_mult[0]->out("\n");

      while (Frame->pause)
        {
        Application->ProcessMessages();
        }

      genoptions_mult[0]->out("SIMULATION CONTINUED\n");
      genoptions_mult[0]->out("\n");
      }
#elif defined(JAVA_OUTPUT_WINDOW)
      bool stop = genoptions_mult[0]->adminb_p->breakcommand();
      if(stop)
        break;
#endif

    } // end: for (i=1;i<=genoptions->iterations;i++)


#if defined(BORLAND_OUTPUT_WINDOW)
    if (!Frame->stop)
#elif defined(JAVA_OUTPUT_WINDOW)
    if (!genoptions_mult[0]->adminb_p->get_stop())
#endif
      {
      return false;
      } // end: if Frame->stop

#if defined(BORLAND_OUTPUT_WINDOW)
    else
      {
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("SIMULATION TERMINATED BY USER BREAK\n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("Estimation results: none\n");
      genoptions_mult[0]->out("\n");

      for(i=0;i<nrmodels;i++)
        {
        if (likepexisting)
          likep_mult[i]->reset();
        }

      for(j=0;j<fullcondp.size();j++)
        fullcondp[j]->reset();

      return true;
      }
#elif defined(JAVA_OUTPUT_WINDOW)
    else
      {
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("Estimation results: none\n");
      genoptions_mult[0]->out("\n");


      for(i=0;i<nrmodels;i++)
        {
        if (likepexisting)
          likep_mult[i]->reset();
        }

      for(j=0;j<fullcondp.size();j++)
        fullcondp[j]->reset();

      return true;
      }
#endif

  } // end: no errors

  return true;

  }


}


















