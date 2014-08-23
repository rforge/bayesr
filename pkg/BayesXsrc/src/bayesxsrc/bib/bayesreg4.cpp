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


//---------------------------------------------------------------------------
#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#include<statwin_haupt.h>

#endif

#include<bayesreg.h>
#include<bayesreg4.h>
#include<typeinfo.h>
#include<stddef.h>


bool bayesreg::create_varcoeffmultibaseline(const unsigned & collinpred)
  {

  long h;
  unsigned min,max,degree,nrknots;
  double lambda,a1,b1;
  bool gl;

  int gridsize;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( varcoeffbaseline.checkvector(terms,i) == true )
      {

      // --------------- reading options, term information ---------------------

      MCMC::fieldtype type;
      type = MCMC::RW2;
      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interacting var
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv); // effectmod

      f = (terms[i].options[1]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[2]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[3]).strtolong(h);
      degree = unsigned(h);

      f = (terms[i].options[4]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[5]).strtodouble(lambda);

      f = (terms[i].options[6]).strtolong(h);
      gridsize = unsigned(h);

      f = (terms[i].options[7]).strtodouble(a1);

      f = (terms[i].options[8]).strtodouble(b1);

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant" && terms[i].options[10] == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;


      datamatrix beg;
      if (begin.getvalue() == "")
        beg = datamatrix(1,1);
      else
        beg = D.getCol(begpos);

      datamatrix statemat;
      if (state.getvalue() == "")
        statemat = datamatrix(1,1);
      else
        statemat = D.getCol(statepos);

      gl=true;

      if (f==1)
        return true;


      // -------------end: reading options, term information -------------------

      //--------- creating path for samples and and results, creating title ----
      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_logbaseline.raw","_logbaseline.res","_logbaseline");

      //----- end: creating path for samples and and results, creating title ---


      fcmultibaseline.push_back( pspline_multibaseline(&generaloptions[generaloptions.size()-1],
                                                distr[distr.size()-1],
                                                fcconst_intercept,
                                                D.getCol(j2),
                                                D.getCol(j1),
                                                nrknots,
                                                degree,
                                                po,
                                                lambda,
                                                min,
                                                max,
                                                type,
                                                title,
                                                pathnonp,
                                                pathres,
                                                gridsize,
                                                collinpred,
                                                statemat,
                                                beg,
                                                gl
                                               )
                             );

      if (constlambda.getvalue() == true)
        fcmultibaseline[fcmultibaseline.size()-1].set_lambdaconst(lambda);

      vector<ST::string> na;
      na.push_back(terms[i].varnames[0]);
      na.push_back(terms[i].varnames[1]);
      fcmultibaseline[fcmultibaseline.size()-1].init_names(na);
      fcmultibaseline[fcmultibaseline.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcmultibaseline[fcmultibaseline.size()-1]);

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
        terms[i].varnames[0],
        "_logbaseline_var.raw","_logbaseline_var.res","_logbaseline_variance");

      fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                              &fcmultibaseline[fcmultibaseline.size()-1],
                              distr[distr.size()-1],a1,b1,title,pathnonp,pathres,
                              false,collinpred)
                              );

             if (constlambda.getvalue() == false)
        {
        if(terms[i].options[9]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
        }

      }

    }

  return false;
  }


