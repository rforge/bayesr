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

#include"mcmcsimul.h"
#include<time.h>
#include"clstring.h"
#include <stdlib.h>
#include<math.h>

using std::ifstream;

namespace MCMC
{


MCMCsimulate::MCMCsimulate(MCMCoptions * go,DISTRIBUTION * dp,
                           vector<FULLCOND*> & fc)
  {
  genoptions_mult.reserve(1);
  likep_mult.reserve(1);
  genoptions_mult.push_back(go);
  likep_mult.push_back(dp);
  fullcondp = fc;
  likepexisting = true;
  begin.reserve(1);
  end.reserve(1);
  begin.push_back(0);
  end.push_back(fc.size()-1);
  }


MCMCsimulate::MCMCsimulate(MCMCoptions * go,DISTRIBUTION * dp,FULLCOND* fc)
  {
  genoptions_mult.reserve(1);
  likep_mult.reserve(1);
  genoptions_mult.push_back(go);
  likep_mult.push_back(dp);
  fullcondp.reserve(1);
  fullcondp.push_back(fc);
  likepexisting = true;
  begin.reserve(1);
  end.reserve(1);
  begin.push_back(0);
  end.push_back(0);
  }


MCMCsimulate::MCMCsimulate(MCMCoptions * go,vector<FULLCOND*> & fc)
  {
  genoptions_mult.reserve(1);
  likep_mult.reserve(1);
  genoptions_mult.push_back(go);
  fullcondp = fc;
  likepexisting = false;
  begin.reserve(1);
  end.reserve(1);
  begin.push_back(0);
  end.push_back(fc.size()-1);
  }

MCMCsimulate::MCMCsimulate(vector<MCMCoptions *> go,vector<DISTRIBUTION *> dp,
                           vector<FULLCOND*>  & fc,vector<unsigned> & be,
                           vector<unsigned> & en)
  {
  genoptions_mult = go;
  likep_mult = dp;
  fullcondp = fc;
  begin = be;
  end = en;
  likepexisting = true;
  }


MCMCsimulate::MCMCsimulate(const MCMCsimulate & s)
  {
//  genoptions = s.genoptions;
//  likep = s.likep;

  fullcondp = s.fullcondp;
  genoptions_mult = s.genoptions_mult;
  likep_mult = s.likep_mult;
  begin = s.begin;
  end = s.end;
  likepexisting = s.likepexisting;
  }


const MCMCsimulate & MCMCsimulate::operator=(const MCMCsimulate & s)
  {
  if (this == &s)
    return *this;

//  genoptions = s.genoptions;
//  likep = s.likep;

  fullcondp = s.fullcondp;
  genoptions_mult = s.genoptions_mult;
  likep_mult = s.likep_mult;
  begin = s.begin;
  end = s.end;
  likepexisting = s.likepexisting;

  return *this;
  }



void MCMCsimulate::set_center(DISTRIBUTION * lp,vector<FULLCOND *> fp,
                              unsigned beg, unsigned en)
  {
  unsigned i,j;

  for(j=0;j<lp->get_responsedim();j++)
    {

    for(i=beg;i<=en;i++)
      {

      if ( (fp[i]->get_col() == j) &&
         (fp[i]->is_identifiable() == false) )
        {
        fp[i]->set_center(true);
        }

      }

    }  // end: for(j=0;j<likep->get_responsedim();j++)

  }


bool MCMCsimulate::checkerrors(DISTRIBUTION * lp,vector<FULLCOND *> fp,
                              unsigned beg, unsigned en)
  {
  bool errors = false;
  unsigned i;

  if (likepexisting)
    {
    if (lp->geterrors().size() > 0)
      {
      errors=true;
      lp->outerrors();
      }
    }

  for (i=beg;i<=en;i++)
    {
    if (fp[i]->geterrors().size() > 0)
      {
      errors = true;
      fp[i]->outerrors();
      }

    }

  return errors;
  }


unsigned MCMCsimulate::compute_nrpar(void)
  {
  unsigned j;
  unsigned nrpar=0;

  if (likepexisting)
    {
    for(j=0;j<likep_mult.size();j++)
      nrpar += likep_mult[j]->get_nrpar();
    }

  for(j=0;j<fullcondp.size();j++)
    {
    if (fullcondp[j]->sample_stored())
      nrpar += fullcondp[j]->get_nrpar();

    }
  return nrpar;
  }


void MCMCsimulate::setflags(const bitset<flagnr> & newflags)
  {
  unsigned j;

  for(j=0;j<fullcondp.size();j++)
    {
    fullcondp[j]->setflags(newflags);
    }

  }


bool MCMCsimulate::simulate(const vector<ST::string> & header, const int & seed,
                                 const bool & computemode)
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
  unsigned iterations = genoptions_mult[0]->get_iterations();


  for (i=0;i<nrmodels;i++)
    {


    genoptions_mult[nrmodels-1-i]->out("\n");
    genoptions_mult[nrmodels-1-i]->out("\n");
    genoptions_mult[nrmodels-1-i]->out(header[nrmodels-1-i] + "\n",true,false,16);
    genoptions_mult[nrmodels-1-i]->out("\n");
    if (i==0)
      {
      genoptions_mult[nrmodels-1-i]->outoptions();
      genoptions_mult[nrmodels-1-i]->out("\n");
      }

    if (likepexisting)
      likep_mult[nrmodels-1-i]->outoptions();

    }

    genoptions_mult[0]->out("OPTIONS FOR ESTIMATION:\n",true);
    genoptions_mult[0]->out("\n");
    for(j=0;j<fullcondp.size();j++)
      fullcondp[j]->outoptions();
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("MCMC SIMULATION STARTED\n",true);
    genoptions_mult[0]->out("\n");

    for (i=0;i<nrmodels;i++)
      {
      if (likepexisting)
        set_center(likep_mult[i],fullcondp,begin[i],end[i]);
      }


  //--------------- Compute posterior mode as starting value -------------------

  if (computemode)
    {
    vector<ST::string> he;
    he.push_back("Computing starting values (may take some time)");
    if (nrmodels > 1)
      {
      for(i=1;i<nrmodels;i++)
        he.push_back("");
      }
    posteriormode(he,true);
    }

  //-------------- end: Compute posterior mode as starting value ---------------


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

  clock_t beginsim = clock();
  clock_t it1per;
  clock_t endsim;
  bool runtime=false;

  #if defined(__BUILDING_GNU)
    double clk = (double)CLOCKS_PER_SEC;
  #else
    double clk = (double)CLK_TCK;
  #endif

  for (it=1;it<=iterations;it++)
    {

    if ( (runtime ==false) && (iterations/it == 100) )
      {
      runtime = true;
      it1per = clock();
      double sec = (it1per-beginsim)/clk;
      long timeleft = long(double(iterations-it)*(double(sec)/double(it)));
      long min = timeleft/60;
      sec = timeleft-min*60;
      if (min == 0)
        {
        genoptions_mult[0]->out("\n");
        genoptions_mult[0]->out
         ("  APPROXIMATE RUN TIME: " + ST::inttostring(long(sec)) +
                        " seconds\n");
        genoptions_mult[0]->out("\n");
        }

      else if (min == 1)
        {
        genoptions_mult[0]->out("\n");
        genoptions_mult[0]->out
                 ("  APPROXIMATE RUN TIME: " + ST::inttostring(min) +
                        " minute " + ST::inttostring(long(sec)) + " seconds\n");
        genoptions_mult[0]->out("\n");
        }
      else
        {
        genoptions_mult[0]->out("\n");
        genoptions_mult[0]->out
                     ("  APPROXIMATE RUN TIME: " + ST::inttostring(min)
                          + " minutes "
                          + ST::inttostring(long(sec)) + " seconds\n");
        genoptions_mult[0]->out("\n");
        }
      }  // end: if ( (runtime ==false) && (iterations/it == 100) )


    for(i=0;i<nrmodels;i++)
      {
      genoptions_mult[nrmodels-1-i]->update();

      if (likepexisting)
        likep_mult[nrmodels-1-i]->update();        // DO NOT CHANGE ORDER !!!!

    for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
      {
      fullcondp[j]->update();
      } // end: for(j=0;j<nrfullcond;j++)

 //   for(i=0;i<nrmodels;i++)
//      {
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
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("SIMULATION TERMINATED\n",true);
      genoptions_mult[0]->out("\n");
      endsim = clock();
      long sec = (endsim-beginsim)/clk;
      long min = sec/60;
      sec = sec-min*60;
      if (min == 0)
        {
        genoptions_mult[0]->out("SIMULATION RUN TIME: "
                                + ST::inttostring(sec) +
                                " seconds\n");
        genoptions_mult[0]->out("\n");
        }
      else if (min == 1)
        {
        genoptions_mult[0]->out("SIMULATION RUN TIME: "
                                + ST::inttostring(min) + " minute "
                                +  ST::inttostring(sec) + " seconds\n");
        genoptions_mult[0]->out("\n");
        }
      else
        {
        genoptions_mult[0]->out("SIMULATION RUN TIME: " +
                                ST::inttostring(min) +  " minutes "
        + ST::inttostring(sec) + " seconds\n");
        genoptions_mult[0]->out("\n");
        }

      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("ESTIMATION RESULTS:\n",true);
      genoptions_mult[0]->out("\n");

      for (i=0;i<nrmodels;i++)
        {
        if (nrmodels > 1)
          {
          genoptions_mult[0]->out("\n");
          genoptions_mult[0]->out(header[i],true);
          genoptions_mult[0]->out("\n");
          }

       if (likep_mult[i]->get_family() == "cox" && likep_mult[i]->get_predict() == true)
          {

//          for(j=begin[i];j<=end[i];j++)
//            if(fullcondp[j]->is_baseline() == true)
          j = begin[i];
          while( !fullcondp[j]->is_baseline() )
            j++;
          ( dynamic_cast<pspline_baseline*>(fullcondp[j]) )->compute_int_ti_mean();
          }
       else if (likep_mult[i]->get_family() == "multistate" && likep_mult[i]->get_predict() == true)
          {
          j = begin[i];
          while( !fullcondp[j]->is_baseline() )
            j++;
          ( dynamic_cast<pspline_multibaseline*>(fullcondp[j]) )->compute_int_ti_mean();
          }

        if (likepexisting)
          likep_mult[i]->outresults();

        for(j=begin[i];j<=end[i];j++)
          fullcondp[j]->outresults();
        }

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



bool MCMCsimulate::posteriormode(const vector<ST::string> & header,
                                 const bool & presim)
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

  if (errors)
    return true;

  for (i=0;i<nrmodels;i++)
    {

    if (header[nrmodels-1-i] !="")
      {
      genoptions_mult[nrmodels-1-i]->out("\n");
      genoptions_mult[nrmodels-1-i]->out("\n");
      genoptions_mult[nrmodels-1-i]->out(header[nrmodels-1-i] + "\n",
                                         true,false,16);

      genoptions_mult[nrmodels-1-i]->out("\n");
      }



    if (!presim)
      {

      if (likepexisting)
        {
        genoptions_mult[nrmodels-1-i]->out("RESPONSE DISTRIBUTION:\n",true);
        genoptions_mult[nrmodels-1-i]->out("\n");
        genoptions_mult[nrmodels-1-i]->out("  "
         + likep_mult[nrmodels-1-i]->get_family() + "\n");
        genoptions_mult[nrmodels-1-i]->out("  Number of observations: " +
        ST::inttostring(likep_mult[nrmodels-1-i]->get_response().rows()) +
        "\n");
        genoptions_mult[nrmodels-1-i]->out("\n");
        }


      if (likepexisting)
        {
        set_center(likep_mult[nrmodels-1-i],fullcondp,begin[nrmodels-i-1],
                        end[nrmodels-1-i]);
        }

      } // end: if (!presim)


    unsigned it=1;
    bool converged;
    bool allconverged;
    bool convergedouterloop=false;

    unsigned k=0;
    while ( (!convergedouterloop) && (k< 100) )
      {

      k++;

      likep_mult[nrmodels-1-i]->compute_iwls();

      converged=false;

      it=1;

      while ((!converged) && (it <= 100))
        {

        allconverged = true;

        if (likepexisting)
          if (likep_mult[nrmodels-1-i]->posteriormode() == false)
             allconverged = false;

        for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
          {
          if (fullcondp[j]->posteriormode() == false)
            allconverged = false;
          } // end: for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)


        if (allconverged)
          converged = true;

        #if defined(BORLAND_OUTPUT_WINDOW)

        Application->ProcessMessages();

        if (Frame->stop)
          {
          break;
          }

        if (Frame->pause)
          {
          genoptions_mult[nrmodels-1-i]->out("\n");
          genoptions_mult[nrmodels-1-i]->out("SIMULATION PAUSED\n");
          genoptions_mult[nrmodels-1-i]->out("Click CONTINUE to proceed\n");
          genoptions_mult[nrmodels-1-i]->out("\n");

          while (Frame->pause)
            {
            Application->ProcessMessages();
            }

          genoptions_mult[nrmodels-1-i]->out("SIMULATION CONTINUED\n");
          genoptions_mult[nrmodels-1-i]->out("\n");
          }
        #elif defined(JAVA_OUTPUT_WINDOW)
        bool stop = genoptions_mult[0]->adminb_p->breakcommand();
        if(stop)
          break;
        #endif

        it++;

        } // end:   while ((!converged) && (it <= 100))


      allconverged = true;

      if (likepexisting)
        if (likep_mult[nrmodels-1-i]->posteriormode_converged(k) == false)
          allconverged = false;

      for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
        {
        if (fullcondp[j]->posteriormode_converged(k) == false)
          allconverged = false;
        } // end: for(j=0;j<nrfullcond;j++)

      if (allconverged)
        convergedouterloop = true;


      if (likepexisting)
        likep_mult[nrmodels-1-i]->posteriormode_set_beta_mode();

      for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
         fullcondp[j]->posteriormode_set_beta_mode();

      } // end: while ( (!convergedouterloop) && (k< 100) )

    if (!presim)
      {
      #if defined(BORLAND_OUTPUT_WINDOW)
      if (!Frame->stop)
      #elif defined(JAVA_OUTPUT_WINDOW)
      if (!genoptions_mult[0]->adminb_p->get_stop())
      #endif
        {
        genoptions_mult[nrmodels-1-i]->out("\n");
        genoptions_mult[nrmodels-1-i]->out("ESTIMATION RESULTS:\n",true);
        genoptions_mult[nrmodels-1-i]->out("\n");

        genoptions_mult[nrmodels-1-i]->out("Number of Iterations: "
         + ST::inttostring(k) + "\n");
        if (!convergedouterloop)
          genoptions_mult[nrmodels-1-i]->out("ALGORITHM DID NOT CONVERGE\n",
                                        true,true,12,255,0,0);
        genoptions_mult[nrmodels-1-i]->out("\n");

        if (likepexisting)
          {
          likep_mult[nrmodels-1-i]->outresults();
          }

        for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
          fullcondp[j]->outresults();

//        return false;

        } // end: if Frame->stop
      #if defined(BORLAND_OUTPUT_WINDOW)
      else
        {

        genoptions_mult[nrmodels-1-i]->out("\n");
        genoptions_mult[nrmodels-1-i]->out(
        "ESTIMATION TERMINATED BY USER BREAK\n");
        genoptions_mult[nrmodels-1-i]->out("\n");
        genoptions_mult[nrmodels-1-i]->out("Estimation results: none\n");
        genoptions_mult[nrmodels-1-i]->out("\n");

        if (likepexisting)
          likep_mult[nrmodels-1-i]->reset();

        for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
          fullcondp[j]->reset();

        errors=true;
//        return true;
        }
      #elif defined(JAVA_OUTPUT_WINDOW)
      else
        {

        genoptions_mult[nrmodels-1-i]->out("\n");
        genoptions_mult[nrmodels-1-i]->out("Estimation results: none\n");
        genoptions_mult[nrmodels-1-i]->out("\n");

        if (likepexisting)
          likep_mult[nrmodels-1-i]->reset();

        for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
          fullcondp[j]->reset();

        errors=true;
//        return true;
        }
      #endif

      } // end: if (!presim)

    } // end:   for (i=0;i<nrmodels;i++)

  return errors;

  } // end: function posteriormode


void MCMCsimulate::autocorr(const unsigned & lag,const ST::string & path)
  {

  unsigned nrmodels = likep_mult.size();

  unsigned p = compute_nrpar();
  ofstream out(path.strtochar());
  assert(!out.fail());
  assert(!out.bad());
  if (p > 0)
    {
    ST::string name;
    datamatrix cmat;
    unsigned i,j,k,l,nr,nrm;
    double min,mean,max,n;
    bool miss,misstot;
    misstot = false;
    genoptions_mult[0]->out("Computing autocorrelation functions...\n");
    autocorr(lag,cmat);


    #if defined(BORLAND_OUTPUT_WINDOW)
    if (!Frame->stop)
    #elif defined(JAVA_OUTPUT_WINDOW)
    if (!genoptions_mult[0]->adminb_p->get_stop())
    #endif
      {
      out << "lag ";


      if (likepexisting)
        {

        for (nrm=0;nrm<nrmodels;nrm++)
          {

          if (likep_mult[nrm]->get_scaleexisting())
            {
            for (k=0;k<likep_mult[nrm]->get_scale_rows();k++)
              for (i=0;i<likep_mult[nrm]->get_scale_cols();i++)
                {
                out << "scale_" << (k+1) << "_" << (i+1) << " ";
                }

            out << "scale_min " << "scale_mean " << "scale_max ";
            }  // end: if (likep_mult[nrm]->get_scaleexisting())

          }  // end: for (i=0;i<nrmodels;i++)

        } // end: if (likepexisting)


      for(j=0;j<fullcondp.size();j++)
        {
        if (fullcondp[j]->sample_stored())
          {
          name = fullcondp[j]->get_title();
          for (k=0;k<fullcondp[j]->getbeta().cols();k++)
            for(i=0;i<fullcondp[j]->getbeta().rows();i++)
              {
              if (fullcondp[j]->getbeta().cols() == 1)
                out << name << "_" << (i+1) << " ";
              else
                out << name << (i+1) << "_" << (k+1) << " ";
              }

          out << name << "_min " << name << "_mean " << name << "_max ";
          } // end: if (fullcondp[j]->sample_stored())

        }  // end: for(j=0;j<fullcondp.size();j++)

      out << endl;

      for(l=0;l<lag;l++)
        {
        nr = 0;
        out << (l+1) << " ";

        if (likepexisting)
          {

          for (nrm=0;nrm<nrmodels;nrm++)
            {

            if (likep_mult[nrm]->get_scaleexisting())
              {

              min = 1;
              max = -1;
              mean = 0;
              miss = true;

              for (k=0;k<likep_mult[nrm]->get_scale_rows();k++)
                for (i=0;i<likep_mult[nrm]->get_scale_cols();i++)
                  {
                  if (cmat(l,nr) <= 1)
                    {
                    miss = false;
                    if (cmat(l,nr) > max)
                      max = cmat(l,nr);
                    if (cmat(l,nr) < min)
                      min = cmat(l,nr);
                    mean+=cmat(l,nr);
                    out << cmat(l,nr) << " ";
                    }
                  else
                    {
                    misstot = true;
                    out << "NA ";
                    }

                  nr++;
                  }

              if (miss)
                out << "NA NA NA ";
              else
                {
                n = likep_mult[nrm]->get_scale_rows()*likep_mult[nrm]->get_scale_cols();
                out << min << " " << (mean/n) << " " << max << " ";
                }

              }  // end: if (likep->get_scaleexisting())

            }  // end: for (i=0;i<nrmodels;i++)

          } // end: if (likepexisting)


        for(j=0;j<fullcondp.size();j++)
          {
          if (fullcondp[j]->sample_stored())
            {

            min = 1;
            max = -1;
            mean = 0;
            miss = true;

            for (k=0;k<fullcondp[j]->getbeta().cols();k++)
              for(i=0;i<fullcondp[j]->getbeta().rows();i++)
                {
                if (cmat(l,nr) <= 1)
                  {
                  miss = false;
                  if (cmat(l,nr) > max)
                    max = cmat(l,nr);
                  if (cmat(l,nr) < min)
                    min = cmat(l,nr);
                  mean+=cmat(l,nr);
                  out << cmat(l,nr) << " ";
                  }
                else
                  {
                  misstot = true;
                  out << "NA ";
                  }

                nr++;
                }

            if (miss)
              out << "NA NA NA ";
            else
              {
              n = fullcondp[j]->getbeta().cols()*fullcondp[j]->getbeta().rows();
              out << min << " " << (mean/n) << " " << max << " ";
              }

            } // end: if (fullcondp[j]->sample_stored())

          }  // end: for(j=0;j<fullcondp.size();j++)
        out << endl;
        } // end: for(l=0;l<lag;l++)

      genoptions_mult[0]->out("Autocorrelation functions computed and stored in file\n");
      genoptions_mult[0]->out(path + ".\n");
      genoptions_mult[0]->out("\n");
      if (misstot)
        {
        genoptions_mult[0]->out("WARNING: There were undefined autocorrelations\n",true,true);
        genoptions_mult[0]->out("\n");
        }

    #if !defined(JAVA_OUTPUT_WINDOW)
      genoptions_mult[0]->out("They may be visualized using the R function 'plotautocor'.\n");
      genoptions_mult[0]->out("\n");
    #endif
      } // end: if (!Frame->stop)
    #if defined(BORLAND_OUTPUT_WINDOW)
    else
      {
      genoptions_mult[0]->out("USER BREAK\n");
      genoptions_mult[0]->out("No autocorrelation functions computed\n");
      genoptions_mult[0]->out("\n");
      out.close();
      remove(path.strtochar());
      }
    #elif defined(JAVA_OUTPUT_WINDOW)
    else
      {
//      genoptions->out("SIMULATION TERMINATED BY USER BREAK\n");
      genoptions_mult[0]->out("No autocorrelation functions computed\n");
      genoptions_mult[0]->out("\n");
      out.close();
      remove(path.strtochar());
      }
    #endif

    }  // end: if (p > 0)
  else
    {
    out << "unable to compute autocorrelation functions" << endl;
    genoptions_mult[0]->outerror("ERROR: No access to sampled parameters\n");
    }

  }   // end: autocorr



void MCMCsimulate::autocorr(const unsigned & lag,datamatrix & cmat)
  {


  unsigned p = compute_nrpar();
  if (p > 0)
    cmat = datamatrix(lag,p);

  unsigned j,i,k,nr;
  unsigned col =0;

  unsigned nrmodels = likep_mult.size();


  if (likepexisting)
    {

    for (nr=0;nr<nrmodels;nr++)
      {

      if (likep_mult[nr]->get_scaleexisting())
        {
        for (k=0;k<likep_mult[nr]->get_scale_rows();k++)
          for (i=0;i<likep_mult[nr]->get_scale_cols();i++)
            {
            cmat.putCol(col,likep_mult[nr]->compute_autocor_scale(lag,k,i));
            col++;
            }
        }  // end: if (likep->get_scaleexisting())


      } // end: for (nr=0;nr<nrmodels;nr++)

    } // end: if (likepexisting)


  for(j=0;j<fullcondp.size();j++)
    {
    if (fullcondp[j]->sample_stored())
      {
      for (k=0;k<fullcondp[j]->getbeta().cols();k++)
        for(i=0;i<fullcondp[j]->getbeta().rows();i++)
          {
          cmat.putCol(col,fullcondp[j]->compute_autocorr(lag,i,k));
          col++;

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

          }

      } // end: if (fullcondp[j]->sample_stored())

    #if defined(BORLAND_OUTPUT_WINDOW)
    Application->ProcessMessages();
    if (Frame->stop)
      {
      cmat = datamatrix(1,1);
      break;
      }
    #elif defined(JAVA_OUTPUT_WINDOW)
//    Application->ProcessMessages();
    if (genoptions_mult[0]->adminb_p->get_stop())
      {
      cmat = datamatrix(1,1);
      break;
      }
    #endif

    } // end:  for(j=0;j<fullcondp.size();j++)


  }


void MCMCsimulate::make_options(ofstream & outtex,const unsigned & nr)
  {

  double it = genoptions_mult[nr]->get_iterations();
  double bu = genoptions_mult[nr]->get_burnin();
  double step = genoptions_mult[nr]->get_step();
  double l1 = genoptions_mult[nr]->get_level1();
  double l2 = genoptions_mult[nr]->get_level2();

  //schreibt MCMC options ins Tex-File
  outtex << "\n\\noindent {\\bf \\large MCMC Options:}" << endl
         << "\\begin{tabbing}" << endl
         << "Levels for credible intervals: \\= \\\\" << endl
         << "Level 1: \\> " << l1 << "\\\\" << endl
         << "Level 2: \\> " << l2 << "\\\\" << endl
         << "Number of Iterations: \\> " << it << "\\\\" << endl
         << "Burn in: \\> " << bu << "\\\\" << endl
         << "Thinning Parameter: \\> " << step << endl
         << "\\end{tabbing}\n"  << "\\vspace{0.5cm}" <<  endl;


  }


void MCMCsimulate::make_predictor(ofstream & outtex,const unsigned & nr)
  {

  unsigned i,j;

  //Anz. d. Kategorien
  unsigned lp = (likep_mult[nr]->get_linearpred()).cols();

  for(i=0; i<lp; i++)                //durchläuft alle Kategorien
    {
    ST::string term = "";         //einführen der Formel für den lin. Prädiktor
    if(lp==1)
      term = "$\\eta$ & $=$ & $";
    else
      term = "$\\eta_" + ST::inttostring(i+1) + "$ & $=$ & $";

    //Schleife, die fullcond-Obj. durchläuft
    bool first = true;
    for(j=begin[nr];j<=end[nr];j++)
      {
      ST::string term2 = fullcondp[j]->get_term_symbolic();
      //nur für richtige fullcond-Obj. der Kateg. "i" (d.h. String nicht leer)
      if(term2 != "" && fullcondp[j]->get_col()==i)
        {
        if (first)
          {
          term = term + term2;            //linearer Prädiktor wird erweitert
          first=false;
          }
        else
          term = term + " + " + term2;    //linearer Prädiktor wird erweitert
        }
      }

    outtex << endl << "\n\\begin{tabular}{ccp{12cm}}\n" << term
           << "$\n\\end{tabular}\n\\\\ \n\\\\" << endl;

    } // end: for(i=0; i<lp; i++)

  }


void MCMCsimulate::make_model(ofstream & outtex,const unsigned & nr)
  {

  //Vert-Fam wird übergeben
  ST::string fam = likep_mult[nr]->get_family();

  fam = fam.replaceallsigns('_', ' ');

  //Anz. Beob. wird übergeben
  unsigned obs = likep_mult[nr]->get_nrobs();

  //Name der Resp.-Var. übergeben
  ST::string resp = likep_mult[nr]->get_responsename();
  char charh = '_';
  ST::string stringh = "\\_";
  ST::string helprname = resp.insert_string_char(charh,stringh);

  //schreibt das Modell und die Prior-annahmen ins Tex-File
  outtex << "\n\\noindent {\\bf \\large Response:}" << endl
         << "\\begin{tabbing}\n"
         << "Number of observations: \\= " << obs << "\\\\" << endl
         << "Response Variable: \\> " << helprname << "\\\\" << endl
         << "Family: \\> " << fam << "\\\\" << endl
         << "\\end{tabbing}" << endl
         << "\n\\noindent {\\bf \\large Predictor:}\\\\" << endl;

  }


void MCMCsimulate::make_prior(ofstream & outtex,const unsigned & nr)
  {

  unsigned i,j;

  outtex << "\n\\noindent {\\bf \\large Priors:}\\\\" << endl << "\\\\" << endl;

  for(j=begin[nr];j<=end[nr];j++)
    {
    vector<ST::string> prior = fullcondp[j]->get_priorassumptions();
    if(prior.size() != 0)// nur wenn Priors da sind (d.h. Vektor hat Elemente)
      {
      for(i=0;i<prior.size();i++)
        {
        outtex << prior[i] << "\\\\" << endl;
        }
      }
    }

  }


void MCMCsimulate::make_fixed_table(ofstream & outtex,const unsigned & nr)
  {

  // falls andere Quantile gewünscht werden
  double u = fullcondp[begin[nr]]->get_level1();
  double o = fullcondp[begin[nr]]->get_level2();
  double u1 = fullcondp[begin[nr]]->get_lower1();
  double u2 = fullcondp[begin[nr]]->get_upper2();
  double o1 = fullcondp[begin[nr]]->get_lower2();
  double o2 = fullcondp[begin[nr]]->get_upper1();
  ST::string u_str = ST::doubletostring(u,0);
  ST::string o_str = ST::doubletostring(o,0);
  ST::string u1_str = ST::doubletostring(u1,5);
  ST::string u2_str = ST::doubletostring(u2,5);
  ST::string o1_str = ST::doubletostring(o1,5);
  ST::string o2_str = ST::doubletostring(o2,5);


  vector<ST::string> h;

  unsigned i,j,l;
  unsigned r;

  unsigned lp = (likep_mult[nr]->get_linearpred()).cols();

  for (l=0;l<lp;l++)
    {

    r = 2;

    // Tabelle im Tex-File mit fixen Effekten
    if (l == 0)
      {
      outtex << "\n\\newpage \n" << endl << "\n\\noindent {\\bf \\large Fixed Effects:}\\\\"
           << endl << "\\\\" << endl;
      }
    else
      {
      outtex << "\n\\newpage \n" << endl <<
            "\n\\noindent {\\bf \\large Fixed Effects (" << ST::inttostring(l+1)
             << ". Response Category):}\\\\"
              << endl << "\\\\" << endl;
      }

    outtex << "\\begin{tabular}{|r|rrrrr|}" << endl << "\\hline" << endl
           << "Variable & Mean & STD & " << u1 << "\\%-Quant. & Median & "
           << u2 << "\\%-Quant.\\\\" << endl << "\\hline" << endl;


    for (i=begin[nr];i<=end[nr];i++)
      {
      if ( (fullcondp[i]->get_results_type()=="fixed") &&
           (fullcondp[i]->get_col()==l)
         )
        {
        h = fullcondp[i]->get_results_latex();
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
                   << "Variable & Mean & STD & " << u1
                   << "\\%-Quant. & Median & " << u2 << "\\%-Quant.\\\\" << endl
                   << "\\hline" << endl;

            outtex << h[j] << endl;

            }
          }

        }

      }

    outtex << "\\hline \n\\end{tabular}" << endl;

    } // end: for (l=0;l<lp;l++)

  }



void MCMCsimulate::make_plots(ofstream & outtex,const unsigned nr,
                              const ST::string & path_batch, const ST::string & path_splus,
                              const ST::string & path_stata)
  {

  char hcharu = '_';
  ST::string hstringu = "\\_";


  unsigned i,j;

  ST::string pathresult;

  bool stil = false;
  bool stata = false;

  // Schleife überprüft, ob es ein fullcond-Object
  // gibt, bei dem Effekt gezeichnet werden kann
  MCMC::plotstyles plst;
  for(j=begin[nr];j<=end[nr];j++)
    {
    plst = fullcondp[j]->get_plotstyle();
    if(plst != MCMC::noplot)
      stil = true;
    if(plst == MCMC::plotnonp)
      stata = true;
    }


  if(stil == true)
    {

    //erzeugt File, das Plot-Befehle für Java-Version enthält
    ofstream outbatch(path_batch.strtochar());

    //erzeugt File, das SPlus-Befehle zum Plotten enthält
    ofstream outsplus(path_splus.strtochar());

    //erzeugt File, das Stata-Befehle zum Plotten enthält
    ofstream outstata;
    if(stata == true)
       outstata.open(path_stata.strtochar());

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

    genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  Batch file for visualizing effects of nonlinear functions is stored in file \n");
    genoptions_mult[0]->out("  " + path_batch + "\n");
    genoptions_mult[0]->out("\n");

    bool stil2 = true;
    for(j=begin[nr];j<=end[nr];j++)  //Schleife überprüft, ob es map-Objekt gibt
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
      genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  Batch file for visualizing effects of nonlinear functions ");
      genoptions_mult[0]->out("  in Stata is stored in file \n");
      genoptions_mult[0]->out("  " + path_stata + "\n");
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
//      genoptions_mult[0]->out("  NOTE: 'input filename' must be substituted by the filename of the boundary-file \n");
//      genoptions_mult[0]->out("\n");
      }
    /*
    if(stata == true)
      {
      genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
      genoptions_mult[0]->out("\n");
      genoptions_mult[0]->out("  Batch file for visualizing effects of nonlinear functions ");
      genoptions_mult[0]->out("  in Stata is stored in file \n");
      genoptions_mult[0]->out("  " + path_stata + "\n");
      genoptions_mult[0]->out("\n");
      }
    */

    outbatch << "% usefile " << path_batch << endl;

    unsigned lp = (likep_mult[nr]->get_linearpred()).cols();

    // falls andere Quantile gewünscht werden
    double u = fullcondp[begin[nr]]->get_level1();
    double o = fullcondp[begin[nr]]->get_level2();
    double u1 = fullcondp[begin[nr]]->get_lower1();
    double u2 = fullcondp[begin[nr]]->get_upper2();
    double o1 = fullcondp[begin[nr]]->get_lower2();
    double o2 = fullcondp[begin[nr]]->get_upper1();
    ST::string u_str = ST::doubletostring(u,0);
    ST::string o_str = ST::doubletostring(o,0);
    ST::string u1_str = ST::doubletostring(u1,5);
    ST::string u2_str = ST::doubletostring(u2,5);
    ST::string o1_str = ST::doubletostring(o1,5);
    ST::string o2_str = ST::doubletostring(o2,5);

    // durchlaufen der Fullconditionals
    for(j=begin[nr];j<=end[nr];j++)
      {

      // Pfad der Regr.-Ergebnisse
      pathresult = fullcondp[j]->get_pathresult();

      // Plotstyle: noplot, plotnonp, drawmap
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

        for(i=0;i<lp;i++)
          {             // für nonparametrische Funktionen
          if (plst == MCMC::plotnonp && fullcondp[j]->get_col()==i)
            {
            outbatch << "\n";                // Befehle f. d. batch-file
            outbatch << "dataset _dat" << endl;
            outbatch << "_dat.infile using " << pathresult << endl;
            outbatch << "graph _g" << endl;
            vector<ST::string> varnames = fullcondp[j]->get_datanames();
            ST::string xvar = varnames[0];
            outbatch << "_g.plot " << xvar
                     << " pmean pqu" << u1_str.replaceallsigns('.','p') << " pqu"
                     << o1_str.replaceallsigns('.','p') << " pqu"
                     << o2_str.replaceallsigns('.','p') << " pqu"
                     << u2_str.replaceallsigns('.','p') << ", "
                     << "title = \"Effect of " << xvar << "\" xlab = " << xvar
                     << " ylab = \" \" " << "outfile = " << pathps
                     << ".ps replace using _dat" << endl;
            outbatch << "drop _dat" << endl;
            outbatch << "drop _g" << endl;
            // Plot-Befehle f. d. SPlus-file
            outsplus << "plotnonp(\"" << pathres_spl << "\")" << endl;

            //Plot-Befehle f. d. Stata-File
            ST::string pu1_str = u1_str.replaceallsigns('.','p');
            ST::string pu2_str = u2_str.replaceallsigns('.','p');
            ST::string po1_str = o1_str.replaceallsigns('.','p');
            ST::string po2_str = o2_str.replaceallsigns('.','p');
            ST::string pu_str = u_str.replaceallsigns('.','p');
            ST::string po_str = o_str.replaceallsigns('.','p');

            outstata << "clear" << endl
                     << "infile intnr " << xvar << " pmean pqu" << pu1_str
                     << " pqu" << po1_str << " pmed pqu" << po2_str << " pqu" << pu2_str
                     << " pcat" << pu_str << " pcat" << po_str << " using "
                     << pathresult << endl
                     << "drop in 1" << endl
                     << "graph twoway rarea pqu" << pu1_str << " pqu" << pu2_str
                     << " " << xvar << ", bcolor(gs13) || rarea pqu" << po1_str
                     << " pqu" << po2_str << " " << xvar << " , bcolor(gs10) || /*"
                     << endl << " */ scatter pmean "
                     //pqu" << pu1_str << " pqu" << po1_str << " pqu" << po2_str << " pqu" << pu2_str << " "
                     << xvar << ", c(l) m(i) clpattern(l) clcolor(gs0) /* "
                     << endl << " */ ytitle(\"Effect of "
                     << xvar << "\") xtitle(\"" << xvar
                     << "\") xlab(,grid) ylab(,grid) legend(off)"
                     // legend(order(3 - \"pointwise credible intervals:\" - 2 1) label(3 \"mean\") holes(2) label(2 \""
                     //<< pu_str << " %\") label(1 \"" << po_str << " %\"))"
                     << endl << "graph export " << pathps << ".eps, replace"
                     << endl << endl;

            // Plot-Befehle f. d. tex-file
            outtex << "\n\\begin{figure}[h!]" << endl
                    << "\\centering" << endl
                    << "\\includegraphics[scale=0.6]{" << pathgr << "}" << endl
                    << "\\caption{Non--linear Effect of '" <<
                    xvar.insert_string_char(hcharu,hstringu) << "'";
            if(lp>1)
                outtex << " (" << ST::inttostring(i+1) << ". response category)";
            outtex << "." << endl << "Shown are the posterior means together with "
                   << u_str << "\\% and " << o_str
                   << "\\% pointwise credible intervals.}" << endl
                   << "\\end{figure}" << endl;
            }
          // für map-Funktionen
          else if ((plst == MCMC::drawmap || plst == MCMC::drawmapgraph) && fullcondp[j]->get_col()==i)
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
                     << "drawmap(data = \"" << pathres_spl << "\", map = m, plotvar = \"pmean\", regionvar = \""
                     << regionvar << "\")" << endl;
            outsplus << "drawmap(data = \"" << pathres_spl << "\", map = m, plotvar = \"pcat" << u_str << "\", regionvar = \""
                    << regionvar << "\", legend = F, pcat = T)" << endl;
            outsplus << "drawmap(data = \"" << pathres_spl << "\", map = m, plotvar = \"pcat" << o_str << "\", regionvar = \""
                      << regionvar << "\", legend = F, pcat = T)" << endl;
            // Plot-Befehle f. d. tex-file
            if(plst == MCMC::drawmap)
              {
              outtex << "\n\\begin{figure}[h!]" << endl
                     << "\\centering" << endl
                     << "\\includegraphics[scale=0.6]{" << pathgr << "_pmean}"
                     << endl
                     << "\\caption{Non--linear Effect of '" <<
                     regionvar.insert_string_char(hcharu,hstringu) << "'";
              if(lp>1)
                outtex << " (" << ST::inttostring(i+1) << ". response category)";
              outtex << ". Shown are the posterior means.}" << endl
                     << "\\end{figure}" << endl;
              outtex << "\n\\begin{figure}[htb]" << endl
                     << "\\centering" << endl
                     << "\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                     << u_str << "}" << endl
                     << "\\caption{Non--linear Effect of '" << regionvar << "'";
              if(lp>1)
                outtex << " (" <<  ST::inttostring(i+1) << ". response category)";
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
              if(lp>1)
                outtex << " (" << ST::inttostring(i+1) << ". response category)";
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
                     << "%\\includegraphics[scale=0.6]{" << pathgr << "_pmean}"
                     << endl
                     << "%\\caption{Non--linear Effect of '" <<
                     regionvar.insert_string_char(hcharu,hstringu) << "'";
              if(lp>1)
                outtex << " (" << ST::inttostring(i+1) << ". response category)";
              outtex << ". Shown are the posterior means.}" << endl
                     << "%\\end{figure}" << endl;
              outtex << "\n%\\begin{figure}[htb]" << endl
                     << "%\\centering" << endl
                     << "%\\includegraphics[scale=0.6]{" << pathgr << "_pcat"
                     << u_str << "}" << endl
                     << "%\\caption{Non--linear Effect of '" << regionvar << "'";
              if(lp>1)
                outtex << " (" <<  ST::inttostring(i+1) << ". response category)";
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
                     << "%\\caption{Non--linear Effect of '" << regionvar << "'";
              if(lp>1)
                outtex << " (" << ST::inttostring(i+1) << ". response category)";
              outtex << ". Posterior probabilities for a nominal level of "
                     << o_str << "\\%." << endl
                     << "%Black denotes regions with strictly negative credible intervals,"
                     << endl
                     << "%white denotes regions with strictly positive credible intervals.}"
                     << endl << "%\\end{figure}" << endl;
              }

            } // end: else if (plst == MCMC::drawmap && fullcondp[j]->get_col()==i)


          } // end: for(i=0;i<lp;i++)


        } // end: if (plst != MCMC::noplot)


      } // end: for(j=begin[nr];j<=end[nr];j++)


    }

  }


 void MCMCsimulate::out_effects(const vector<ST::string> & paths)
   {

/*   unsigned i,j,k;

   unsigned nrmodels = genoptions_mult.size();

   datamatrix e;


   for (i=0;i<nrmodels;i++)
     {

     vector<unsigned> be;
     vector<unsigned> en;
     unsigned nr=0;

     ofstream oute(paths[nrmodels-1-i].strtochar());

     for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++)
       {
       be.push_back(nr);
       nr+=fullcondp[j]->get_nreffects(MCMC::mean);
       en.push_back(nr);
       }

     e = datamatrix(likep_mult[nrmodels-1-i]->get_nrobs(),nr);
     vector<ST::string> enames;
//    vector<unsigned>::iterator itbe = be.begin();
//    vector<unsigned>::iterator iten = en.begin();
     unsigned it=0;
     for(j=begin[nrmodels-1-i];j<=end[nrmodels-1-i];j++,it++)
       {
       fullcondp[j]->get_effectmatrix(e,enames,be[it],en[it],MCMC::mean);
       }

     for(j=0;j<enames.size();j++)
       {
       oute << enames[j] << "   ";
       }

     oute << endl;

     for(j=0;j<e.rows();j++)
       {
       for(k=0;k<e.cols();k++)
         oute << e(j,k) << "   ";
       oute << endl;
       }

     }  */

   }


void MCMCsimulate::make_graphics(const vector<ST::string> & title,
const vector<ST::string> & path_batch,
const vector<ST::string> & path_tex,
const vector<ST::string> & path_splus,
const vector<ST::string> & path_stata)
  {

  unsigned i,nr;

  unsigned nrmodels = likep_mult.size();

  ST::string pathresult;                 //Pfad des Ergebnis-Files

  vector<ST::string> distr_results;

  for (nr=0;nr<nrmodels;nr++)
    {

    // erzeugt Tex-File
    ofstream outtex(path_tex[nrmodels-nr-1].strtochar());

    //erzeugt den Kopf des Tex-Files
    outtex << "\\documentclass[a4paper, 12pt]{article}" << endl
           << "\n" << "\\usepackage{graphicx}" << endl
           << "\\parindent0em" << endl
           << "\n\\begin{document}" << endl
           << "\\begin{center}" << endl
           << "\\LARGE{\\bf " << title[nrmodels-1-nr] << "}"
           << endl << "\\end{center} \n\\vspace{1cm}" << endl;

    make_model(outtex,nrmodels-1-nr);

    make_predictor(outtex,nrmodels-1-nr);

    make_prior(outtex,nrmodels-1-nr);

    make_options(outtex,nrmodels-1-nr);

    distr_results = likep_mult[nrmodels-1-nr]->get_results_latex();
    for (i=0;i<distr_results.size();i++)
      {
      outtex << distr_results[i] << endl;
      }

    make_fixed_table(outtex,nrmodels-1-nr);

    // Pfade der Files
    //werden im BayesX-Output angegeben
    genoptions_mult[0]->out("  Files of model summary: \n" , true);
    genoptions_mult[0]->out("\n");


    make_plots(outtex,nrmodels-1-nr,path_batch[nrmodels-1-nr],
           path_splus[nrmodels-1-nr],path_stata[nrmodels-1-nr]);


    genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  Latex file of model summaries is stored in file \n");
    genoptions_mult[0]->out("  " + path_tex[nrmodels-1-nr] + "\n");
    genoptions_mult[0]->out("\n");
    genoptions_mult[0]->out("  --------------------------------------------------------------------------- \n");
    genoptions_mult[0]->out("\n");


    outtex << "\\end{document}" << endl;

    }

  }


void MCMCsimulate::get_samples(
#if defined(JAVA_OUTPUT_WINDOW)
vector<ST::string> & newc,
#endif
const ST::string & path,const unsigned & step)
  {
  unsigned i;
  ST::string filename;
  ST::string help;
  ST::string psname;
  genoptions_mult[0]->out("Storing sampled parameters...\n");
  genoptions_mult[0]->out("Sampled parameters are stored in file(s):\n");
  genoptions_mult[0]->out("\n");

  for(i=0;i<fullcondp.size();i++)
    {
    if (fullcondp[i]->sample_stored())
      {
      filename = path + fullcondp[i]->get_title() + "_sample.raw";
      fullcondp[i]->get_samples(filename,step);
      genoptions_mult[0]->out(filename + "\n");
      #if defined(JAVA_OUTPUT_WINDOW)

      psname = path + fullcondp[i]->get_title() + "_sample.ps";
      newc.push_back("dataset _dat");
      newc.push_back("_dat.infile , nonote using " + filename);
      newc.push_back("graph _g");
      newc.push_back("_g.plotsample , replace outfile=" +
                  psname  + " using _dat");
      genoptions_mult[0]->out(psname + " (graphs)\n");
      newc.push_back("drop _dat _g");

      #endif
      genoptions_mult[0]->out("\n");
      }
    }

  if (likepexisting)
    {
    for(i=0;i<likep_mult.size();i++)
      {

      if (likep_mult[i]->get_predictfull())
        {
        filename = likep_mult[i]->get_mean_sample();
        genoptions_mult[0]->out(filename+"\n");
        genoptions_mult[0]->out("\n");
        #if defined(JAVA_OUTPUT_WINDOW)

        psname = filename.substr(0,filename.length()-4) +   + ".ps";
        newc.push_back("dataset _dat");
        newc.push_back("_dat.infile , nonote using " + filename);
        newc.push_back("graph _g");
        newc.push_back("_g.plotsample , replace outfile=" +
                      psname  + " using _dat");
        genoptions_mult[0]->out(psname + " (graphs)\n");
        newc.push_back("drop _dat _g");

        #endif
        }

      if (likep_mult[i]->get_scaleexisting())
        {
        filename = likep_mult[i]->get_scale_sample();
        genoptions_mult[0]->out(filename+"\n");
        genoptions_mult[0]->out("\n");
        #if defined(JAVA_OUTPUT_WINDOW)

        psname = filename.substr(0,filename.length()-4) +   + ".ps";
        newc.push_back("dataset _dat");
        newc.push_back("_dat.infile , nonote using " + filename);
        newc.push_back("graph _g");
        newc.push_back("_g.plotsample , replace outfile=" +
                      psname  + " using _dat");
        genoptions_mult[0]->out(psname + " (graphs)\n");
        newc.push_back("drop _dat _g");

        #endif

        }

// Begin: DSB
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
      if (likep_mult[i]->get_mscheck())
        {
        likep_mult[i]->get_mssamples();
        }
#endif
// End: DSB
      }
    }

  genoptions_mult[0]->out("\n");
  genoptions_mult[0]->out("Storing completed\n");
  genoptions_mult[0]->out("\n");
  #if defined(BORLAND_OUTPUT_WINDOW)
  genoptions_mult[0]->out(
  "Sampled parameters may be visualized using the R\n");
  genoptions_mult[0]->out("function 'plotsample'.\n");
  genoptions_mult[0]->out("\n");
  #endif

  }


//------------------------------------------------------------------------------
//------------------------- end: class MCMCsimulate ----------------------------
//------------------------------------------------------------------------------


void compare(const vector<ST::string> & files,ostream & out)
  {
  unsigned i,j;
  vector<ST::string> token;
  ST::string li;
  datamatrix m0,m1;
  datamatrix result;

  for(i=0;i<files.size();i++)
    {
    ifstream in(files[i].strtochar());
    ST::getline(in,li);
    token = li.strtoken(" ");
    if (i==0)
      {
      result = datamatrix(files.size(),token.size());
      m0.prettyScan(in);
      m1 = m0;
      }
    else
      {
      m1.prettyScan(in);
      }

    for(j=0;j<m0.cols();j++)
      result(i,j) = norm(m1-m0,j)/norm(m0,j);
    }

  for (j=0;j<token.size();j++)
    out << token[j] << "   ";
  out << endl;
  result.prettyPrint(out);

  }


} // end: namespace MCMC

#if defined(BORLAND_OUTPUT_WINDOW)
//---------------------------------------------------------------------------
#pragma package(smart_init)
#endif






