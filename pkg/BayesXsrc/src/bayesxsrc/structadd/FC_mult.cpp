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



#include "FC_mult.h"


//------------------------------------------------------------------------------
//----------------- CLASS: FC_mult implementation of member functions ----------
//------------------------------------------------------------------------------


namespace MCMC
{

void FC_mult::set_effectp(DESIGN * d,FC_nonp * fp)
  {
  FCnp = fp;
  dp2 = d;
  }

void FC_mult::set_intp(DESIGN * d,FC_nonp * fp)
  {
  dp1 = d;
  FCnp2 = fp;
  }


void FC_mult::set_multeffects(MASTER_OBJ * mp,unsigned & enr, GENERAL_OPTIONS * o,
                              const ST::string & t, const ST::string & fp,
                              bool sm,bool meane, double meanec)
  {

  masterp = mp;
  equationnr = enr;

  unsigned rows = dp1->Zout.rows()*dp2->Zout.rows();

  samplemult = sm;
  compmeaneffect = meane;
  meaneffectconstant = meanec;

  if (samplemult)
    {
    FCmulteffect = FC(o,t,rows,1,fp);
    if (compmeaneffect==true)
      FCmulteffect_mean = FC(o,"",rows,1,fp);
    }

  }


void FC_mult::compute_autocorr_all(const ST::string & path, unsigned lag,
                                    ofstream & outg) const
  {

  }


FC_mult::FC_mult(void)
  {
  nosamples = true;
  samplemult = false;
  compmeaneffect = false;
  meaneffectconstant = 0;
  multexp=false;
  }


FC_mult::FC_mult(bool reu,bool mexp)
     : FC()
  {
  multexp = mexp;
  nosamples = true;
  samplemult = false;
  compmeaneffect = false;
  meaneffectconstant = 0;
  RE_update = reu;
  }


FC_mult::FC_mult(const FC_mult & m)
  : FC(FC(m))
  {
  equationnr = m.equationnr;
  masterp = m.masterp;
  multexp = m.multexp;
  FCmulteffect = m.FCmulteffect;
  FCmulteffect_mean = m.FCmulteffect_mean;
  samplemult = m.samplemult;
  compmeaneffect = m.compmeaneffect;
  meaneffectconstant = m.meaneffectconstant;
  dp1 = m.dp1;
  dp2 = m.dp2;
  FCnp = m.FCnp;
  FCnp2 = m.FCnp2;
  RE_update = m.RE_update;
  effect = m.effect;
  }



const FC_mult & FC_mult::operator=(const FC_mult & m)
  {

  if (this==&m)
	 return *this;
  FC::operator=(FC(m));
  equationnr = m.equationnr;
  masterp = m.masterp;
  multexp = m.multexp;
  FCmulteffect = m.FCmulteffect;
  FCmulteffect_mean = m.FCmulteffect_mean;
  samplemult = m.samplemult;
  compmeaneffect = m.compmeaneffect;
  meaneffectconstant = m.meaneffectconstant;
  dp1 = m.dp1;
  dp2 = m.dp2;
  FCnp = m.FCnp;
  FCnp2 = m.FCnp2;
  RE_update = m.RE_update;
  effect = m.effect;
  return *this;
  }


void FC_mult::update(void)
  {
  double add;

  if (RE_update)
    {
    add=0;
    }
  else
    {
    add=1;
    }

  dp2->compute_effect(effect,FCnp->beta,MCMC::Function);
  dp1->set_intvar(effect,add);

  if ((RE_update==false) && (samplemult))
    {
    update_multeffect();
    FCmulteffect.acceptance++;
    FCmulteffect.update();
    if (compmeaneffect)
      {
      FCmulteffect_mean.acceptance++;
      FCmulteffect_mean.update();
      }
    }

  }


void FC_mult::update_multeffect(void)
  {
  /*
  FCnp RE
  dp2 RE

  FCnp2 nonl
  dp1 nonl
  */

  unsigned i,j;
  double * mebetap = FCmulteffect.beta.getV();
  double * FCnpbetap = FCnp->beta.getV();
  double * FCnp2betap;

  // TEST
  //  ofstream out("c:\\bayesx\\testh\\results\\beta_RE.res");
  //  FCnp->beta.prettyPrint(out);
  // TEST

  // TEST
  // ofstream out2("c:\\bayesx\\testh\\results\\betaspline.res");
  // FCnp2->beta.prettyPrint(out2);
  // TEST

  if (compmeaneffect)
    {
    double meanlin;
    if (meaneffectconstant==0)
      meanlin = masterp->level1_likep[equationnr]->meaneffect-FCnp2->meaneffect;
    else
      meanlin = meaneffectconstant;

    double lin;
    double * FCmp = FCmulteffect_mean.beta.getV();
    for (i=0;i<FCnp->beta.rows();i++,FCnpbetap++)
      {
      FCnp2betap = FCnp2->beta.getV();
      for (j=0;j<FCnp2->beta.rows();j++,mebetap++,FCnp2betap++,FCmp++)
        {
        *mebetap = (*FCnpbetap+1)*(*FCnp2betap);
        lin = meanlin + *mebetap;
        masterp->level1_likep[equationnr]->compute_mu(&lin,FCmp);
        }
      }
    }
  else
    {
    for (i=0;i<FCnp->beta.rows();i++,FCnpbetap++)
      {
      FCnp2betap = FCnp2->beta.getV();
      for (j=0;j<FCnp2->beta.rows();j++,mebetap++,FCnp2betap++)
        *mebetap = (*FCnpbetap+1)*(*FCnp2betap);
      }

    }

  }


bool FC_mult::posteriormode(void)
  {

  if (multexp==false)
    {
    double add;

    if (RE_update)
      {
      add=0;
      }
    else
      {
      add=1;
      }

    dp2->compute_effect(effect,FCnp->beta,MCMC::Function);
    dp1->set_intvar(effect,add);

    if ((RE_update==false) && (samplemult))
      {
      update_multeffect();
      bool h = FCmulteffect.posteriormode();

      if (compmeaneffect)
        {
        h = FCmulteffect_mean.posteriormode();
        }

      }
    }
  else
    {

    /*
    FCnp RE
    dp2 RE

    FCnp2 nonl
    dp1 nonl
    */


    double add=0;


    dp2->compute_effect(effect,FCnp->beta,MCMC::Function);

    if (RE_update)
      {
      }
    else
      {
      unsigned i;
      double * effectp = effect.getV();
      for (i=0;i<effect.rows();i++,effectp++)
        *effectp = exp(*effectp);
      }


    dp1->set_intvar(effect,add);

    if ((RE_update==false) && (samplemult))
      {
      update_multeffect();
      bool h = FCmulteffect.posteriormode();

      if (compmeaneffect)
        {
        h = FCmulteffect_mean.posteriormode();
        }

      }

    }

  return true;
  }


void FC_mult::outgraphs(ofstream & out_stata, ofstream & out_R,const ST::string & path)
  {

  ST::string pathps = path.substr(0,path.length()-4) + "_statagraph";

  double u = FCmulteffect.optionsp->level1;
  double o = FCmulteffect.optionsp->level2;
  double u1 = FCmulteffect.optionsp->lower1;
  double u2 = FCmulteffect.optionsp->upper2;
  double o1 = FCmulteffect.optionsp->lower2;
  double o2 = FCmulteffect.optionsp->upper1;
  ST::string u_str = ST::doubletostring(u,0);
  ST::string o_str = ST::doubletostring(o,0);
  ST::string u1_str = ST::doubletostring(u1,5);
  ST::string u2_str = ST::doubletostring(u2,5);
  ST::string o1_str = ST::doubletostring(o1,5);
  ST::string o2_str = ST::doubletostring(o2,5);

  ST::string pu1_str = u1_str.replaceallsigns('.','p');
  ST::string pu2_str = u2_str.replaceallsigns('.','p');
  ST::string po1_str = o1_str.replaceallsigns('.','p');
  ST::string po2_str = o2_str.replaceallsigns('.','p');
  ST::string pu_str = u_str.replaceallsigns('.','p');
  ST::string po_str = o_str.replaceallsigns('.','p');

  ST::string xvar1 = dp1->datanames[0] + "   ";
  ST::string xvar2 = dp2->datanames[0] + "   ";

  out_stata << "clear" << endl
            << "infile intnr " << xvar1 << xvar2
            << " pmean pqu" << pu1_str
            << " pqu" << po1_str << " pmed pqu" << po2_str << " pqu" << pu2_str
            << " pcat" << pu_str << " pcat" << po_str;



  if (compmeaneffect==true)
    {
    out_stata << " pmean_mu pqu"
              << pu1_str << "_mu"
              << " pqu" << po1_str << "_mu"
              << " pmed_d pqu" << po2_str << "_mu"
              << " pqu" << pu2_str << "_mu";
    }


  out_stata << " using "
            << path << endl
            << "drop in 1" << endl;

  out_stata << "scatter pmean " << xvar2 <<  endl << endl;


  }




void FC_mult::outresults(ofstream & out_stata, ofstream & out_R,
                         const ST::string & pathresults)
  {
  if ((RE_update==false) && (samplemult))
    {
    FCmulteffect.outresults(out_stata,out_R,"");

    if (compmeaneffect)
      FCmulteffect_mean.outresults(out_stata,out_R,"");

    if (pathresults.isvalidfile() != 1)
      {

      outgraphs(out_stata,out_R,pathresults);

      FCmulteffect.optionsp->out("    Results are stored in file\n");
      FCmulteffect.optionsp->out("    " +  pathresults + "\n");
      FCmulteffect.optionsp->out("\n");

      ofstream outres(pathresults.strtochar());

      FCmulteffect.optionsp->out("\n");

      unsigned i,j;

      ST::string l1 = ST::doubletostring(FCmulteffect.optionsp->lower1,4);
      ST::string l2 = ST::doubletostring(FCmulteffect.optionsp->lower2,4);
      ST::string u1 = ST::doubletostring(FCmulteffect.optionsp->upper1,4);
      ST::string u2 = ST::doubletostring(FCmulteffect.optionsp->upper2,4);
      l1 = l1.replaceallsigns('.','p');
      l2 = l2.replaceallsigns('.','p');
      u1 = u1.replaceallsigns('.','p');
      u2 = u2.replaceallsigns('.','p');

      outres << "intnr" << "   ";

      /*
      FCnp RE
      dp2 RE

      FCnp2 nonl
      dp1 nonl
     */


      outres << dp1->datanames[0] << "   ";
      outres << dp2->datanames[0] << "   ";


      outres << "pmean   ";

      if (FCmulteffect.optionsp->samplesize > 1)
        {
        outres << "pqu"  << l1  << "   ";
        outres << "pqu"  << l2  << "   ";
        outres << "pmed   ";
        outres << "pqu"  << u1  << "   ";
        outres << "pqu"  << u2  << "   ";
        outres << "pcat" << FCmulteffect.optionsp->level1 << "   ";
        outres << "pcat" << FCmulteffect.optionsp->level2 << "   ";
        }

      if (compmeaneffect)
        {
        outres << "pmean_mu   ";

        if (FCmulteffect.optionsp->samplesize > 1)
          {
          outres << "pqu"  << l1  << "_mu   ";
          outres << "pqu"  << l2  << "_mu   ";
          outres << "pmed_mu   ";
          outres << "pqu"  << u1  << "_mu   ";
          outres << "pqu"  << u2  << "_mu   ";
          }
        }


      outres << endl;

      double * workmean = FCmulteffect.betamean.getV();
      double * workbetaqu_l1_lower_p = FCmulteffect.betaqu_l1_lower.getV();
      double * workbetaqu_l2_lower_p = FCmulteffect.betaqu_l2_lower.getV();
      double * workbetaqu_l1_upper_p = FCmulteffect.betaqu_l1_upper.getV();
      double * workbetaqu_l2_upper_p = FCmulteffect.betaqu_l2_upper.getV();
      double * workbetaqu50 = FCmulteffect.betaqu50.getV();

      double * workmean_mu=NULL;
      double * workbetaqu_l1_lower_p_mu=NULL;
      double * workbetaqu_l2_lower_p_mu=NULL;
      double * workbetaqu_l1_upper_p_mu=NULL;
      double * workbetaqu_l2_upper_p_mu=NULL;
      double * workbetaqu50_mu=NULL;

      if (compmeaneffect)
        {
        workmean_mu = FCmulteffect_mean.betamean.getV();
        workbetaqu_l1_lower_p_mu = FCmulteffect_mean.betaqu_l1_lower.getV();
        workbetaqu_l2_lower_p_mu = FCmulteffect_mean.betaqu_l2_lower.getV();
        workbetaqu_l1_upper_p_mu = FCmulteffect_mean.betaqu_l1_upper.getV();
        workbetaqu_l2_upper_p_mu = FCmulteffect_mean.betaqu_l2_upper.getV();
        workbetaqu50_mu = FCmulteffect_mean.betaqu50.getV();
        }

      for (i=0;i<FCnp->beta.rows();i++)
        {
        for (j=0;j<FCnp2->beta.rows();j++,workmean++,workbetaqu_l1_lower_p++,
                              workbetaqu_l2_lower_p++,workbetaqu50++,
                              workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
          {
          outres << (i+1) << "   ";
          outres << dp2->effectvalues[i] << "   ";
          outres << dp1->effectvalues[j] << "   ";
          outres << *workmean << "   ";

          if (FCmulteffect.optionsp->samplesize > 1)
            {
            outres << *workbetaqu_l1_lower_p << "   ";
            outres << *workbetaqu_l2_lower_p << "   ";
            outres << *workbetaqu50 << "   ";
            outres << *workbetaqu_l2_upper_p << "   ";
            outres << *workbetaqu_l1_upper_p << "   ";

            if (*workbetaqu_l1_lower_p > 0)
              outres << 1 << "   ";
            else if (*workbetaqu_l1_upper_p < 0)
              outres << -1 << "   ";
            else
              outres << 0 << "   ";

            if (*workbetaqu_l2_lower_p > 0)
              outres << 1 << "   ";
            else if (*workbetaqu_l2_upper_p < 0)
              outres << -1 << "   ";
            else
              outres << 0 << "   ";
            }


          if (compmeaneffect)
            {

            outres << *workmean_mu << "   ";

            if (FCmulteffect.optionsp->samplesize > 1)
              {
              outres << *workbetaqu_l1_lower_p_mu << "   ";
              outres << *workbetaqu_l2_lower_p_mu << "   ";
              outres << *workbetaqu50_mu << "   ";
              outres << *workbetaqu_l2_upper_p_mu << "   ";
              outres << *workbetaqu_l1_upper_p_mu << "   ";

              }

            if ( !( (i == FCnp->beta.rows()-1)  && (j==FCnp2->beta.rows()-1) ) )
              {
              workmean_mu++;
              workbetaqu_l1_lower_p_mu++;
              workbetaqu_l2_lower_p_mu++;
              workbetaqu50_mu++;
              workbetaqu_l1_upper_p_mu++;
              workbetaqu_l2_upper_p_mu++;
              }

            }


          outres << endl;
          }
        }
      }
    }
  }


void FC_mult::reset(void)
  {

  }


} // end: namespace MCMC



