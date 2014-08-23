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
#include<statwin_haupt.h>

#endif

#include"bayesreg.h"
#include"bayesreg3.h"

// Vorschlag:
//#include<typeinfo.h>

#include<stddef.h>

using std::ios;

bool bayesreg::create_varcoeffmerror(const unsigned & collinpred)
  {

  ST::string monotone;
  ST::string proposal;

  long h;
  unsigned min,max,degree,nrknots;
  double lambda, mvar1, mvar2, arva, arpa1, arpa2, bmean, bvar;
  double a1,b1;
  int gridsize, contourprob;
  int f;
  bool diagtransform,derivative,center;

  unsigned i,j,k;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpvarcoeffmerror.checkvector(terms,i) == true )
      {

      MCMC::fieldtype type;
      if (terms[i].options[0] == "merrorrw1")
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      // extract interaction variables
      ST::string test = terms[i].varnames[0].substr(0,terms[i].varnames[0].length()-7);
      datamatrix medata = datamatrix(D.rows(),2,0);
      for(k=0; k<2; k++)
        {
        j1 = (test+ST::inttostring(k+1)).isinlist(modelvarnamesv);
        medata.putCol(k, D.getCol(j1));
        }
      terms[i].varnames[0] = test;

      datamatrix meandata = datamatrix(D.rows(),1,0);
      unsigned mecols = medata.cols();
      for(k=0; k<D.rows(); k++)
        {
        for(j=0; j<mecols; j++)
          {
          meandata(k,0) = medata(k,j);
          }
        meandata(k,0) /= mecols;
        }

      // extract effect modifier
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

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

      proposal = terms[i].options[9];
      monotone = terms[i].options[10];

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant" && terms[i].options[19] == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      if (terms[i].options[14] == "false")
        diagtransform = false;
      else
        diagtransform = true;

      if (terms[i].options[15] == "false")
        derivative = false;
      else
        derivative = true;

      f = (terms[i].options[16]).strtolong(h);
      contourprob = unsigned(h);

      if (terms[i].options[20] == "false")
        center = false;
      else
        center = true;

      // SUSI: read new opttion form the options list

      f = (terms[i].options[21]).strtodouble(mvar1);
      f = (terms[i].options[22]).strtodouble(mvar2);
      f = (terms[i].options[23]).strtodouble(arva);
      f = (terms[i].options[24]).strtodouble(arpa1);
      f = (terms[i].options[25]).strtodouble(arpa2);
      f = (terms[i].options[26]).strtodouble(bmean);
      f = (terms[i].options[27]).strtodouble(bvar);

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_pspline.raw","_pspline.res","_pspline");

      if(proposal == "cp")
        {

        fcpspline.push_back( FULLCOND_pspline(&generaloptions[generaloptions.size()-1],
                                              distr[distr.size()-1],
                                              fcconst_intercept,
                                              D.getCol(j2),
                                              medata,
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
                                              derivative,
                                              gridsize,
                                              collinpred
                                             )
                           );

        if (constlambda.getvalue() == true)
          fcpspline[fcpspline.size()-1].set_lambdaconst(lambda);

        vector<ST::string> na;
        na.push_back(terms[i].varnames[1]);
        na.push_back(terms[i].varnames[0]);
        fcpspline[fcpspline.size()-1].init_names(na);
        fcpspline[fcpspline.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcpspline[fcpspline.size()-1]);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                   terms[i].varnames[0],"_pspline_var.raw",
                   "_pspline_var.res","_pspline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                     &fcpspline[fcpspline.size()-1],
                                     distr[distr.size()-1],
                                     a1,
                                     b1,
                                     title,pathnonp,pathres,
                                     false,collinpred)
                           );
         if (constlambda.getvalue() == false)
          {
          if(terms[i].options[17]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }
        }
      else
        {

        bool iwlsmode;
        if(proposal == "iwlsmode")
          iwlsmode = true;
        else
          iwlsmode = false;
        f = (terms[i].options[11]).strtolong(h);
        unsigned updateW;
        updateW = unsigned(h);
        bool updatetau;
        if(terms[i].options[12] == "false" || constlambda.getvalue() == true)
          updatetau = false;
        else
          updatetau = true;
        double fstart;
          f = (terms[i].options[13]).strtodouble(fstart);

        fciwlspspline.push_back( IWLS_pspline(&generaloptions[generaloptions.size()-1],
                                              distr[distr.size()-1],
                                              fcconst_intercept,
                                              D.getCol(j2),
                                              medata,
                                              iwlsmode,
                                              nrknots,
                                              degree,
                                              po,
                                              lambda,
                                              type,
                                              monotone,
                                              updateW,
                                              updatetau,
                                              fstart,
                                              a1,b1,
                                              title,
                                              pathnonp,
                                              pathres,
                                              derivative,
                                              gridsize,
                                              diagtransform,
                                              collinpred
                                             )
                             );

        if (constlambda.getvalue() == true)
          fciwlspspline[fciwlspspline.size()-1].set_lambdaconst(lambda);

        vector<ST::string> na;
        na.push_back(terms[i].varnames[1]);
        na.push_back(terms[i].varnames[0]);
        fciwlspspline[fciwlspspline.size()-1].init_names(na);
        fciwlspspline[fciwlspspline.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fciwlspspline[fciwlspspline.size()-1]);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                     terms[i].varnames[0],"_pspline_var.raw",
                     "_pspline_var.res","_pspline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fciwlspspline[fciwlspspline.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                             );

        if (constlambda.getvalue() == false)
          {
          if(updatetau)
            fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();
          if(terms[i].options[17]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        }

        // Fullcond-Objekt zur Generierung der wahren Kovariablenwerte

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],"_merror.raw",
                 "_merror.res","_merror");

      //SUSI: add value to the merror fullcond
      fcmerror.push_back(fullcond_merror(&generaloptions[generaloptions.size()-1],
                         &fciwlspspline[fciwlspspline.size()-1],
                         distr[distr.size()-1],
                                   medata,
                                   D.getCol(j2),
                                   title,
                                   pathres,
                                   mvar1, mvar2, arva, arpa1, arpa2, bmean, bvar)
                         );

//      fcmerror[fcmerror.size()-1].setmerroroptions(mvar1, mvar2, arva, arpa1, arpa2);
      fcmerror[fcmerror.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcmerror[fcmerror.size()-1]);
      }

    }


  return false;
  }

bool bayesreg::create_varcoeffpspline(const unsigned & collinpred)
  {

  ST::string monotone;
  ST::string proposal;

  long h;
  unsigned min,max,degree,nrknots;
  double lambda;
  double a1,b1;
  int gridsize,contourprob;
  int f;
  bool diagtransform,derivative;
  bool center;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpvarcoeffpspline.checkvector(terms,i) == true )
      {

      MCMC::fieldtype type;
      if (terms[i].options[0] == "varpsplinerw1")
        type = MCMC::RW1;
      else
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

      proposal = terms[i].options[9];
      monotone = terms[i].options[10];

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant" && terms[i].options[19] == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      if (terms[i].options[14] == "false")
        diagtransform = false;
      else
        diagtransform = true;

      if (terms[i].options[15] == "false")
        derivative = false;
      else
        derivative = true;

      f = (terms[i].options[16]).strtolong(h);
      contourprob = unsigned(h);

      if (terms[i].options[20] == "false")
        center = false;
      else
        center = true;


      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_pspline.raw","_pspline.res","_pspline");

      if (check_gaussian(collinpred))
        {

        fcpsplinegaussian.push_back( FULLCOND_pspline_gaussian(&generaloptions[generaloptions.size()-1],
                                              distr[distr.size()-1],
                                              fcconst_intercept,
                                              D.getCol(j2),
                                              D.getCol(j1),
                                              nrknots,
                                              degree,
                                              po,
                                              type,
                                              monotone,
                                              title,
                                              pathnonp,
                                              pathres,
                                              derivative,
                                              lambda,
                                              gridsize,
                                              center,
                                              collinpred
                                             )
                           );

        //
        datamatrix beta_0;
        if(terms[i].options[18]!="")
          {
          dataobject * datap;                           // pointer to datasetobject
          int objpos = findstatobject(*statobj,terms[i].options[18],"dataset");
          if (objpos >= 0)
            {
            statobject * s = statobj->at(objpos);
            datap = dynamic_cast<dataobject*>(s);
            if (datap->obs()==0 || datap->getVarnames().size()==0)
              {
              outerror("ERROR: dataset object " + terms[i].options[18] + " does not contain any data\n");
              return true;
              }
            else if (datap->getVarnames().size()>1)
              {
              outerror("ERROR: dataset object " + terms[i].options[18] + " contains more than one variable\n");
              return true;
              }
            }
          else
            {
            outerror("ERROR: dataset object " + terms[i].options[18] + " is not existing\n");
            return true;
            }
          list<ST::string> names = datap->getVarnames();
          ST::string expr = "";
          datap->makematrix(names,beta_0,expr);
          }
        else
          {
          beta_0 = datamatrix(1,1,0);
          }
        //

        fcpsplinegaussian[fcpsplinegaussian.size()-1].set_contour(contourprob,
            pseudocontourprob.getvalue(),approx.getvalue(),lengthstart.getvalue(),beta_0);

        if (constlambda.getvalue() == true)
          fcpsplinegaussian[fcpsplinegaussian.size()-1].set_lambdaconst(lambda);

        vector<ST::string> na;
        na.push_back(terms[i].varnames[1]);
        na.push_back(terms[i].varnames[0]);
        fcpsplinegaussian[fcpsplinegaussian.size()-1].init_names(na);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                   terms[i].varnames[0],"_pspline_var.raw",
                   "_pspline_var.res","_pspline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fcpsplinegaussian[fcpsplinegaussian.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

        fcpsplinegaussian[fcpsplinegaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcpsplinegaussian[fcpsplinegaussian.size()-1]);

        if (constlambda.getvalue() == false)
          {
          if(terms[i].options[17]=="true")
                fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        }


      else
        {

        if(proposal == "cp")
          {

          fcpspline.push_back( FULLCOND_pspline(&generaloptions[generaloptions.size()-1],
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
                                              derivative,
                                              gridsize,
                                              collinpred
                                             )
                           );

          if (constlambda.getvalue() == true)
            fcpspline[fcpspline.size()-1].set_lambdaconst(lambda);

          vector<ST::string> na;
          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          fcpspline[fcpspline.size()-1].init_names(na);
          fcpspline[fcpspline.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpspline[fcpspline.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                     terms[i].varnames[0],"_pspline_var.raw",
                     "_pspline_var.res","_pspline_variance");


          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fcpspline[fcpspline.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                             );

          if (constlambda.getvalue() == false)
            {
            if(terms[i].options[17]=="true")
                fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          }
        else
          {

          bool iwlsmode;
          if(proposal == "iwlsmode")
            iwlsmode = true;
          else
            iwlsmode = false;
          f = (terms[i].options[11]).strtolong(h);
          unsigned updateW;
          updateW = unsigned(h);
          bool updatetau;
          if(terms[i].options[12] == "false" || constlambda.getvalue() == true)
            updatetau = false;
          else
            updatetau = true;
          double fstart;
            f = (terms[i].options[13]).strtodouble(fstart);

          fciwlspspline.push_back( IWLS_pspline(&generaloptions[generaloptions.size()-1],
                                                distr[distr.size()-1],
                                                fcconst_intercept,
                                                D.getCol(j2),
                                                D.getCol(j1),
                                                iwlsmode,
                                                nrknots,
                                                degree,
                                                po,
                                                lambda,
                                                type,
                                                monotone,
                                                updateW,
                                                updatetau,
                                                fstart,
                                                a1,b1,
                                                title,
                                                pathnonp,
                                                pathres,
                                                derivative,
                                                gridsize,
                                                diagtransform,
                                                collinpred
                                               )
                             );

          if (constlambda.getvalue() == true)
            fciwlspspline[fciwlspspline.size()-1].set_lambdaconst(lambda);

          vector<ST::string> na;
          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          fciwlspspline[fciwlspspline.size()-1].init_names(na);
          fciwlspspline[fciwlspspline.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fciwlspspline[fciwlspspline.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                     terms[i].varnames[0],"_pspline_var.raw",
                     "_pspline_var.res","_pspline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fciwlspspline[fciwlspspline.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                             );

          if (constlambda.getvalue() == false)
            {
            if(updatetau)
              fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();
            if(terms[i].options[17]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          }

        }

      }

    }

  return false;
  }


void mregressrun(bayesreg & b)
  {

  vector<ST::string> header;

  b.resultsyesno = false;
  if (b.modeonly.getvalue() == true)
    b.posteriormode = true;
  else
    b.posteriormode = false;

  b.termsmult = b.modregmult.getterms();

  b.describetext.erase(b.describetext.begin(),b.describetext.end());
  b.describetext.push_back("LAST ESTIMATED MODEL: \n");
  b.describetext.push_back("\n");
  b.describetext.push_back(b.modregmult.getModelText());
  b.describetext.push_back("\n");

  b.clear();

  b.outfiles.push_back(b.outfile.getvalue()+b.add_name);

  bool failure = false;

  if ((b.family.getvalue() != "multgaussian") && (b.family.getvalue() != "multistate"))
    {
    failure = true;
    b.out("ERROR: family " + b.family.getvalue() + " is not allowed for method mregress\n");
    }

  if (!failure)
    failure = b.create_generaloptions();

  if (!failure)
    failure = b.create_distribution();

  unsigned i,j;
  unsigned nrbaseline,nrbaseline_i;
  unsigned nrbaseline_help=0;

  if (!failure)
    {
    for (i=0;i<b.nrcategories;i++)
      {

      b.terms = b.termsmult[i];

      if (!failure)
        failure = b.create_const(i);

      if (!failure)
        failure = b.create_multibaseline(i);

      if (!failure)
        failure = b.create_varcoeffmultibaseline(i);

      if (!failure)
        failure = b.create_nonprw1rw2(i);

      if (!failure)
        failure = b.create_pspline(i);

      if (!failure)
        failure = b.create_spatial(i);

      if (!failure)
        failure = b.create_random(i);

      nrbaseline = b.fcmultibaseline.size();
      nrbaseline_i = nrbaseline - nrbaseline_help;

      if(nrbaseline_i>1)
        {
        vector<MCMC::pspline_multibaseline*> basep;
        for(j=0;j<nrbaseline_i;j++)
          basep.push_back(&b.fcmultibaseline[nrbaseline_help+j]);
        for(j=0;j<nrbaseline_i;j++)
          b.fcmultibaseline[nrbaseline_help+j].set_baselinep(basep);
        }

      nrbaseline_help = nrbaseline;

      }
    }


  if (!failure)
    {

    b.simobj = MCMCsimulate(&b.generaloptions[0],b.distr[0],b.fullcond);

    if (b.modeonly.getvalue())
      {
//      vector<ST::string> header;
      header.push_back("BAYESREG OBJECT " + b.name.to_bstr() +
                                    ": regression procedure");

      failure = b.simobj.posteriormode(header);
      }
    else
      {
//      vector<ST::string> header;
      header.push_back("BAYESREG OBJECT " + b.name.to_bstr() +
                                    ": regression procedure");
      failure = b.simobj.simulate(header,b.setseed.getvalue(),!b.noposteriormode.getvalue());
      }
    }


  if (!failure)
    {
    b.resultsyesno = true;
    }
  else
    {
    b.describetext.erase(b.describetext.begin(),b.describetext.end());
    b.describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
    b.resultsyesno = false;
    }

  if (!failure && (b.missingreg.getvalue()==false))
    {

    vector<ST::string> path;
    vector<ST::string> path2;
    vector<ST::string> path3;
    vector<ST::string> path4;

    for (i=0;i<b.outfiles.size();i++)
      {
      path.push_back(b.outfiles[i] + "_graphics.prg");
      path2.push_back(b.outfiles[i] + "_model_summary.tex");
      path3.push_back(b.outfiles[i] +  "_r.R");
      path4.push_back(b.outfiles[i] +  "_stata.do");
      }



    b.simobj.make_graphics(header,path,path2,path3,path4);

#if defined(JAVA_OUTPUT_WINDOW)
    if(b.nographs.getvalue() == false)
    {
    for(unsigned j=0;j<b.fullcond.size();j++)
       {
       MCMC::plotstyles plst = b.fullcond[j]->get_plotstyle();
       if(plst != MCMC::noplot)
         {
         vector<ST::string> varnames = b.fullcond[j]->get_datanames();
         ST::string xvar = varnames[0];
         ST::string pathresult = b.fullcond[j]->get_pathresult();
         ST::string pathps = pathresult.substr(0, pathresult.length()-4);
         if(plst == MCMC::plotnonp)
                 {
                 b.newcommands.push_back(b.name + ".plotnonp " + ST::inttostring(j)
                 + ", title = \"Effect of " + xvar +"\" xlab = " + xvar
                 + " ylab = \" \" outfile = " + pathps + ".ps replace");
                 }

         if(plst==MCMC::drawmap)  // || plst==MCMC::drawmapgraph)
                 {
                 double u = b.fullcond[j]->get_level1();
                 double o = b.fullcond[j]->get_level2();
                 ST::string u_str = ST::doubletostring(u,0);
                 ST::string o_str = ST::doubletostring(o,0);
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", color outfile = " + pathps + "_pmean.ps replace");
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", plotvar = pcat" + u_str + " nolegend  pcat outfile = " + pathps
                 + "_pcat" + u_str + ".ps replace");
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", plotvar = pcat" + o_str + " nolegend  pcat outfile = " + pathps
                 + "_pcat" + o_str + ".ps replace");

                 }
         }
       }

    b.newcommands.push_back(b.name + ".texsummary");
    }
#endif
    }

  }



void outresultsrun(bayesreg & b)
  {

  if (b.resultsyesno==true)
    {

    unsigned i;
    ST::string suffix;
    ST::string header;
    ST::string trtype=b.transformtype.getvalue();
    ST::string trtypeintercept;

    bool error=false;

    if (b.transformtype.getvalue() == "exp")
      {
      trtypeintercept = "exp";
      suffix = "_exp.res";
      header = "ESTIMATION RESULTS FOR PARAMETER TRANSFORMATION \"exp\"";
      }
    else if (b.transformtype.getvalue() == "oddsratio")
      {
      if ( (b.distrstring[0]=="binomial") ||
           (b.distrstring[0]=="binomlogitlat")
         )
        {
        trtypeintercept = "oddsratiointercept";
        suffix = "_oddsratio.res";
        header = "ESTIMATION RESULTS FOR PARAMETER TRANSFORMATION \"oddsratio\"";
        }
      else
        {
        b.outerror("ERROR: odds ratio transformation only allowed for logit models\n");
        error=true;
        }

      }
    else if (b.transformtype.getvalue() == "marginal")
      {
      trtypeintercept = "marginalintercept";
      suffix = "_marginal.res";
      header = "ESTIMATION RESULTS FOR MARGINAL EFFECTS";

      }
    else if (b.transformtype.getvalue() == "elasticity")
      {
      trtypeintercept = "none";
      suffix = "_elasticity.res";
      header = "ESTIMATION RESULTS FOR ELASTICITIES";
      }
    else if (b.transformtype.getvalue() == "lognormal")
      {
      if (b.distrstring[0]=="gaussian")
        {
        trtypeintercept = "lognormalintercept";
        suffix = "_lognormal.res";
        header = "ESTIMATION RESULTS FOR PARAMETER TRANSFORMATION \"lognormal\"";
        }
      else
        {
        b.outerror("ERROR: lognormal transformation only allowed for Gaussian models\n");
        error=true;
        }
      }
    else if (b.transformtype.getvalue() == "logit")
      {
      if ( (b.distrstring[0]=="binomial") ||
           (b.distrstring[0]=="binomlogitlat")
         )
        {
        trtypeintercept = "logitintercept";
        suffix = "_logit.res";
        header = "ESTIMATION RESULTS FOR PARAMETER TRANSFORMATION \"logit\"";
        }
      else
        {
        b.outerror("ERROR: logit transformation only allowed for logit models\n");
        error=true;
        }
      }
    else if (b.transformtype.getvalue() == "probit")
      {
      if ( b.distrstring[0]=="binomlat")
        {
        trtypeintercept = "probitintercept";
        suffix = "_probit.res";
        header = "ESTIMATION RESULTS FOR PARAMETER TRANSFORMATION \"probit\"";
        }
      else
        {
        b.outerror("ERROR: probit transformation only allowed for probit models\n");
        error=true;
        }
      }


    if (error==false)
      {

      for(i=0;i<b.fcbaseline.size();i++)
        b.fcbaseline[i].set_transform(suffix,trtype);

//      for(i=0;i<b.fcbaselineiwls.size();i++)
//        b.fcbaselineiwls[i].set_transform(suffix,trtype);

      for(i=0;i<b.fcmultibaseline.size();i++)
        b.fcmultibaseline[i].set_transform(suffix,trtype);

      for(i=0;i<b.fcpspline.size();i++)
        b.fcpspline[i].set_transform(suffix,trtype);

      for(i=0;i<b.fciwlspspline.size();i++)
        b.fciwlspspline[i].set_transform(suffix,trtype);

      for(i=0;i<b.fcpsplinegaussian.size();i++)
        b.fcpsplinegaussian[i].set_transform(suffix,trtype);


      if (b.transformtype.getvalue() != "elasticity")
        {
        for(i=0;i<b.fcpsplinesurfgaussian.size();i++)
          b.fcpsplinesurfgaussian[i].set_transform(suffix,trtype);

        for(i=0;i<b.normalconst.size();i++)
          b.normalconst[i].set_transform(suffix,trtype);

        for(i=0;i<b.nongaussianconst.size();i++)
          b.nongaussianconst[i].set_transform(suffix,trtype);

        for(i=0;i<b.fcnonp.size();i++)
          b.fcnonp[i].set_transform(suffix,trtype);

        for(i=0;i<b.fcnonpgaussian.size();i++)
          b.fcnonpgaussian[i].set_transform(suffix,trtype);


        for(i=0;i<b.fcpsplinesurf.size();i++)
          b.fcpsplinesurf[i].set_transform(suffix,trtype);

        for(i=0;i<b.fcrandom.size();i++)
          b.fcrandom[i].set_transform(suffix,trtype);

        for(i=0;i<b.fcrandomgaussian.size();i++)
          b.fcrandomgaussian[i].set_transform(suffix,trtype);

        for(i=0;i<b.fckriging.size();i++)
          b.fckriging[i].set_transform(suffix,trtype);

        }


      b.out("\n");
      b.out(header,true,false);
      b.out("\n");

//      b.distr[0]->outintercept(trtypeintercept,suffix);

      vector<FULLCOND *> fcp(1);

      for (i=0;i<b.fullcond.size();i++)
        {
        if (b.fullcond[i]->transform_yes() == true)
          {
          fcp[0] = b.fullcond[i];
          ST::string helptr = b.transformtype.getvalue();
          b.distr[0]->transform_nonlinear(fcp,helptr);
          b.fullcond[i]->outresults();
          }

        }

#if defined(JAVA_OUTPUT_WINDOW)
      if(b.nographs.getvalue() == false)
      {
      for(unsigned j=0;j<b.fullcond.size();j++)
         {
         MCMC::plotstyles plst = b.fullcond[j]->get_plotstyle();
         if(plst != MCMC::noplot)
           {
           vector<ST::string> varnames = b.fullcond[j]->get_datanames();
           ST::string xvar = varnames[0];
           ST::string pathresult = b.fullcond[j]->get_pathcurrent();
           ST::string pathps = pathresult.substr(0, pathresult.length()-4);
           if(plst == MCMC::plotnonp)
                 {
                 b.newcommands.push_back(b.name + ".plotnonp " + ST::inttostring(j)
                 + ", title = \"Effect of " + xvar +"\" xlab = " + xvar
                 + " ylab = \" \" outfile = " + pathps + ".ps replace");
                 }

           if(plst==MCMC::drawmap)  // || plst==MCMC::drawmapgraph)
                 {
                 double u = b.fullcond[j]->get_level1();
                 double o = b.fullcond[j]->get_level2();
                 ST::string u_str = ST::doubletostring(u,0);
                 ST::string o_str = ST::doubletostring(o,0);
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", color outfile = " + pathps + "_pmean.ps replace");
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", plotvar = pcat" + u_str + " nolegend  pcat outfile = " + pathps
                 + "_pcat" + u_str + ".ps replace");
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", plotvar = pcat" + o_str + " nolegend  pcat outfile = " + pathps
                 + "_pcat" + o_str + ".ps replace");

                 }
           }
         }
      }
#endif

      } // end: if (error==false)

    } // end: if (b.resultsyesno==true)

  }


void autocorrrun(bayesreg & b)
  {

  if (b.resultsyesno==true)
	 {
     if (b.posteriormode == false)
       {
       ST::string path = b.outfile.getvalue()+  "_autocor" + ".raw";
       if (b.generaloptions[0].get_samplesize() < unsigned(b.maxlag.getvalue()*4))
         b.outerror("ERROR: samplesize too small\n");
       else
         b.simobj.autocorr(b.maxlag.getvalue(),path);
       }
     else
       b.outerror("ERROR: no MCMC simulation results\n");
	 }
  else
	 b.outerror("ERROR: no regression results\n");

  }


void drawmaprun(bayesreg & b)
  {

#if defined(BORLAND_OUTPUT_WINDOW)

  b.outerror("ERROR: method drawmap is not available in this version\n");

#elif defined(__BUILDING_GNU)

  b.outerror("ERROR: method drawmap is not available in this version\n");

#elif defined(JAVA_OUTPUT_WINDOW)

  bool error = false;

  vector<ST::string> varnames = b.mdrawmap.getModelVarnamesAsVector();
  if (varnames.size() != 1)
    {
    b.outerror("ERROR: syntax error for method drawmap\n");
    error = true;
    }

  long nr;
  if (varnames[0].strtolong(nr) != 0)
    {
    b.outerror("ERROR: syntax error for method drawmap\n");
    error = true;
    }

  if (nr < 0 || nr >= b.fullcond.size())
    {
    b.outerror("ERROR: syntax error for method drawmap\n");
    error = true;
    }

  if (error == false)
    {
    if (b.fullcond[nr]->get_plotstyle() != MCMC::drawmap
                      && b.fullcond[nr]->get_plotstyle() != MCMC::drawmapgraph)
      {
      error = true;
      b.outerror("ERROR: results cannot be visualized with method drawmap\n");
      }
    else if (b.fullcond[nr]->get_plotstyle() == MCMC::drawmapgraph)
      {
      error = true;
      b.outerror("ERROR: boundaries of the regions are not available from the graph-file \n");
      }
    }

  if (error==false)
    {

    ST::string path = b.fullcond[nr]->get_pathresult();

    vector<ST::string> vnames;
    ifstream in(path.strtochar());
    ST::string h;
    ST::getline(in,10000,h);
    h = h.eatallcarriagereturns();

    vnames = 	h.strtoken(" ");

    ST::string graphname = "_" + b.name + "_graph";
    b.newcommands.push_back("graph " + graphname);

    ST::string datasetname = "_" + b.name + "_r0";
    b.newcommands.push_back("dataset " + datasetname);
    b.newcommands.push_back(datasetname + ".infile , nonote using " + path);

    ST::string plotvar;
      plotvar = b.plotvar.getvalue() + " " + vnames[1] + " ";

    ST::string ot="map=" + b.fullcond[nr]->getinfo() + " ";

    ot = ot + "nrcolors="+b.nrcolors.getValueAsString()+" ";
    ot = ot + "title=\""+b.title2.getvalue() + "\" ";
    if (b.outfile4.getvalue().length() > 0)
      ot = ot + "outfile=\""+b.outfile4.getvalue() + "\" ";
    if (b.nolegend.getvalue() == true)
      ot = ot + "nolegend ";
    if (b.color.getvalue() == true)
      ot = ot + "color ";
    if (b.swapcolors.getvalue() == true)
      ot = ot + "swapcolors ";
    if (b.replace.getvalue() == true)
      ot = ot + "replace ";
    if (b.lowerlimit.changed() == true)
      ot = ot + "lowerlimit="  + b.lowerlimit.getValueAsString() + " ";
    if (b.upperlimit.changed() == true)
      ot = ot + "upperlimit=" + b.upperlimit.getValueAsString() + " ";
    if (b.pcat.getvalue() == true)
      ot = ot + "pcat ";
    if (b.drawnames.getvalue() == true)
      ot = ot + "drawnames ";
    if (b.hclcolors.getvalue() == true)
      ot = ot + "hcl ";
    if (b.fontsize.changed() == true)
      ot = ot + "fontsize=" + b.fontsize.getValueAsString() + " ";
    if (b.titlescale.changed() == true)
      ot = ot + "titlesize=" + b.titlescale.getValueAsString() + " ";

    if (ot.length() == 0)
      b.newcommands.push_back(
      graphname + ".drawmap " + plotvar + " using " + datasetname);
    else
      b.newcommands.push_back(
      graphname + ".drawmap " + plotvar + "," + ot + " using "
      + datasetname);

    b.newcommands.push_back("drop " + graphname + " " + datasetname);
    }

#endif

  }


void plotnonprun(bayesreg & b)
  {



#if defined(BORLAND_OUTPUT_WINDOW)

  b.outerror("ERROR: method plotnonp is not available in this version\n");

#elif defined(__BUILDING_GNU)

  b.outerror("ERROR: method plotnonp is not available in this version\n");

#elif defined(JAVA_OUTPUT_WINDOW)

  bool error = false;

  vector<ST::string> varnames = b.mplotnonp.getModelVarnamesAsVector();
  if (varnames.size() != 1)
    {
    b.outerror("ERROR: syntax error for method plotnonp\n");
    error = true;
    }

  long nr;
  if (varnames[0].strtolong(nr) != 0)
    {
    b.outerror("ERROR: syntax error for method plotnonp\n");
    error = true;
    }

  if (nr < 0 || nr >= b.fullcond.size())
    {
    b.outerror("ERROR: syntax error for method plotnonp\n");
    error = true;
    }

  if (error == false)
    {
    if (b.fullcond[nr]->get_plotstyle() != MCMC::plotnonp)
      {
      error = true;
      b.outerror("ERROR: results cannot be visualized with method plotnonp\n");
      }
    }

  if (error==false)
    {

    ST::string path = b.fullcond[nr]->get_pathcurrent();

    if(path == "")
      path = b.fullcond[nr]->get_pathresult();

    vector<ST::string> vnames;
    ifstream in(path.strtochar());
    ST::string h;
    ST::getline(in,10000,h);
    h = h.eatallcarriagereturns();
    vnames = 	h.strtoken(" ");

    ST::string graphname = "_" + b.name + "_graph";
    b.newcommands.push_back("graph " + graphname);

    ST::string datasetname = "_" + b.name + "_r0";
    b.newcommands.push_back("dataset " + datasetname);
    b.newcommands.push_back(datasetname + ".infile , nonote using " + path);

    ST::string plotvar;
    if (b.levels.getvalue()=="all")
      {
      plotvar = vnames[1] + " ";
      if (b.median.getvalue() == true)
        plotvar = plotvar + vnames[5] + " ";
      else
        plotvar = plotvar + vnames[2] + " ";
      plotvar = plotvar + vnames[3] + " " +
                         vnames[4] + " " +
                         vnames[6] + " " +
                         vnames[7] + " ";
      }
    else if (b.levels.getvalue()=="none")
      {
      if (b.median.getvalue() == true)
        plotvar = vnames[1] + " " + vnames[5];
      else
        plotvar = vnames[1] + " " + vnames[2];

      }
    else if (b.levels.getvalue()=="1")
      {
      if (b.median.getvalue() == true)
        plotvar = vnames[1] + " " + vnames[5] + " " + vnames[3] + " " +
                  vnames[7] + " ";
      else
        plotvar = vnames[1] + " " + vnames[2] + " " + vnames[3] + " " +
                  vnames[7] + " ";
      }
    else
      {
      if (b.median.getvalue() == true)
        plotvar = vnames[1] + " " + vnames[5] + " " + vnames[4] + " " +
                  vnames[6] + " ";
      else
        plotvar = vnames[1] + " " + vnames[2] + " " + vnames[4] + " " +
                  vnames[6] + " ";
      }


    ST::string ot;
    ot = "xlab=\""+b.xlab.getvalue() + "\" ";
    ot = ot + "ylab=\""+b.ylab.getvalue() + "\" ";
    ot = ot + "title=\""+b.title0.getvalue() + "\" ";
    if (b.outfile2.getvalue().length() > 0)
      ot = ot + "outfile=\""+b.outfile2.getvalue() + "\" ";
    ot = ot + "height="+b.height.getValueAsString() + " ";
    ot = ot + "width="+b.width.getValueAsString() + " ";
    if (b.replace2.getvalue() == true)
      ot = ot + " replace ";
    if (b.connect.changed() == true)
      ot = ot + "connect="+b.connect.getvalue() + " ";
    if (b.ylimbottom.changed() == true)
      ot = ot + "ylimbottom="+b.ylimbottom.getValueAsString() + " ";
    if (b.ylimtop.changed() == true)
      ot = ot + "ylimtop="+b.ylimtop.getValueAsString() + " ";
    if (b.ystep.changed() == true)
      ot = ot + "ystep="+b.ystep.getValueAsString() + " ";
    if (b.ystart.changed() == true)
      ot = ot + "ystart="+b.ystart.getValueAsString() + " ";
    if (b.xlimbottom.changed() == true)
      ot = ot + "xlimbottom="+b.xlimbottom.getValueAsString() + " ";
    if (b.xlimtop.changed() == true)
      ot = ot + "xlimtop="+b.xlimtop.getValueAsString() + " ";
    if (b.xstep.changed() == true)
      ot = ot + "xstep="+b.xstep.getValueAsString() + " ";
    if (b.xstart.changed() == true)
      ot = ot + "xstart="+b.xstart.getValueAsString() + " ";
    if (b.linewidth.changed() == true)
      ot = ot + "linewidth="+b.linewidth.getValueAsString() + " ";
    if (b.fontsize.changed() == true)
      ot = ot + "fontsize="+b.fontsize.getValueAsString() + " ";
    if (b.pointsize.changed() == true)
      ot = ot + "pointsize="+b.pointsize.getValueAsString() + " ";
    if (b.linecolor.changed() == true)
      ot = ot + "linecolor="+b.linecolor.getValueAsString() + " ";
    if (b.titlescale.changed() == true)
      ot = ot + "titlesize="+b.titlescale.getValueAsString() + " ";

    if (ot.length() == 0)
      b.newcommands.push_back(graphname + ".plot " + plotvar + " using " + datasetname);
    else
      b.newcommands.push_back(graphname + ".plot " + plotvar + "," + ot + " using "
      + datasetname);

    b.newcommands.push_back("drop " + graphname + " " + datasetname);
    }

#endif

  b.plotnonpoptions.setdefault();

  }


void texsummaryrun(bayesreg & b)
  {

#if defined(BORLAND_OUTPUT_WINDOW)

  b.outerror("ERROR: method texsummary is not available in this version\n");

#elif defined(JAVA_OUTPUT_WINDOW)

  ST::string path = b.outfiles[0];
  ST::string path2 = path;

  int i = path2.length()-1;
  bool gefunden = false;
  while(i>=0 && gefunden == false)
    {
    if(path2[i] == '\\' || path2[i]=='/')
      gefunden = true;
    path2 = path2.substr(0,i);
    i--;
    }

  ST::string helpbat = path2 + "_latexcommands.bat";
  ofstream outbat(helpbat.strtochar());
  outbat << "cd " << path2 << endl;
  outbat << path.substr(0,1) << ":" << endl;
  outbat << "latex " << path << "_model_summary.tex" << endl;
  outbat << "dvips " << path << "_model_summary.dvi" << endl;
  outbat.close();
  system(helpbat.strtochar());
  remove(helpbat.strtochar());
#endif

  }


bool bayesreg::create_nonpseason(const unsigned & collinpred)
  {

  ST::string pathnonpv;
  ST::string pathresv;

  bool iwls;
  unsigned updateW;
  bool updatetau;

  ST::string proposal;

  long h;
  double ftune;
  unsigned min;
  unsigned max;
  unsigned per;
  double lambda,a1,b1;
  int f;

  unsigned i;
  int j1=0,j2=0;
  bool varcoeff;

  for(i=0;i<terms.size();i++)
    {
    if ( nonpseason.checkvector(terms,i) == true )
      {

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);

      if (terms[i].type=="season")
        varcoeff = false;
      else
        {
        varcoeff = true;
        j2 = terms[i].varnames[1].isinlist(modelvarnamesv);
        }

      f = (terms[i].options[1]).strtolong(h);
      per = unsigned(h);
      f = (terms[i].options[2]).strtolong(h);
      min = unsigned(h);
      f = (terms[i].options[3]).strtolong(h);
      max = unsigned(h);
      f = (terms[i].options[4]).strtodouble(lambda);
      f = (terms[i].options[5]).strtodouble(a1);
      f = (terms[i].options[6]).strtodouble(b1);
      f = (terms[i].options[8]).strtolong(h);
      updateW = unsigned(h);
      f = (terms[i].options[10]).strtodouble(ftune);
      proposal = terms[i].options[7];
      if (terms[i].options[9] == "true")
        updatetau=true;
      else
        updatetau=false;

      if (f==1)
        return true;

      ST::string titlev;

      if (varcoeff==false)
        {
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_season.raw","_season.res","_season");

        make_paths(collinpred,pathnonpv,pathresv,titlev,terms[i].varnames[0],"",
                   "_season.raw","_season_var.res","_season_variance");
        }
      else
        {
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                   terms[i].varnames[0],
                   "_season.raw","_season.res","_season");

        make_paths(collinpred,pathnonpv,pathresv,titlev,terms[i].varnames[1],
                   terms[i].varnames[0],
                   "_season.raw","_season_var.res","_season_variance");

        }

      if (proposal != "cp")
        iwls=true;
      else
        iwls=false;

      if ((check_gaussian(collinpred)) || (check_iwls(iwls,collinpred)) )
        {

        if (varcoeff==false)
          {
          fcnonpgaussian.push_back( FULLCOND_nonp_gaussian(
                       &generaloptions[generaloptions.size()-1],
                                         distr[distr.size()-1],
                                         D.getCol(j1),
                                         fcconst_intercept,
                                         unsigned(maxint.getvalue()),
                                         MCMC::seasonal,
                                         title,
                                         pathnonp,
                                         pathres,
                                         collinpred,lambda,per
                                        )
                          );

          fcnonpgaussian[fcnonpgaussian.size()-1].init_name(
                  terms[i].varnames[0]);
          }
        else
          {

          fcnonpgaussian.push_back(
          FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
          distr[distr.size()-1], D.getCol(j2),D.getCol(j1),
          fcconst_intercept,
          unsigned(maxint.getvalue()),MCMC::seasonal,title,pathnonp,
          pathres,collinpred,lambda,per)
                                );

          vector<ST::string> na;
          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);

          }

        if (constlambda.getvalue() == true)
          {
          if (check_nongaussian(collinpred))
            fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(updateW,true);
          fcnonpgaussian[fcnonpgaussian.size()-1].set_lambdaconst(lambda);
          }
        else
          {
          if ( (check_nongaussian(collinpred)) && (proposal == "iwls")
              && (updatetau==false) )
            fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(updateW);

          if ( (check_nongaussian(collinpred)) && (proposal == "iwlsmode")
             && (updatetau==false) )
             fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS(updateW,true);

          if ( (check_nongaussian(collinpred)) && (proposal == "iwls")
              && (updatetau==true) )
            fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS_hyperblock(
                                                                 updateW,a1,b1);

          if ( (check_nongaussian(collinpred)) && (proposal == "iwlsmode")
              && (updatetau==true) )
            fcnonpgaussian[fcnonpgaussian.size()-1].set_IWLS_hyperblock(
                                                    updateW,a1,b1,true);
          }


        fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);


        if (constlambda.getvalue() == false)
          {
          fcvarnonp.push_back(FULLCOND_variance_nonp(
                                         &generaloptions[generaloptions.size()-1],
                                         &fcnonpgaussian[fcnonpgaussian.size()-1],
                                         distr[distr.size()-1],
                                         a1,
                                         b1,
                                         titlev,pathnonpv,pathresv,
                                         false,collinpred)
                             );

          if ( (check_nongaussian(collinpred)) && (proposal == "iwls")
            && (updatetau==true) )
            fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();

          if ( (check_nongaussian(collinpred)) && (proposal == "iwlsmode")
            && (updatetau==true) )
            fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();

          if(terms[i].options[14]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();


          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        }
      else
        {

        if (varcoeff==false)
          {
          Pmatrices.push_back(PenaltyMatrix(D.getCol(j1),
                                            terms[i].varnames[0],
                                            unsigned(maxint.getvalue()),
                                            min,
                                            max,
                                            MCMC::seasonal,
                                            per
                                            )
                              );
           }
         else
           {
           Pmatrices.push_back(PenaltyMatrix(D.getCol(j2),
                                            terms[i].varnames[1],
                                            unsigned(maxint.getvalue()),
                                            min,
                                            max,
                                            MCMC::seasonal,
                                            per
                                            )
                              );

           }

        fcnonp.push_back( FULLCOND_nonp(&generaloptions[generaloptions.size()-1],
                                         distr[distr.size()-1],
                                         &Pmatrices[Pmatrices.size()-1],
                                         fcconst_intercept,
                                         lambda,
                                         pathnonp,
                                         pathres,title," ",
                                         collinpred
                                        )
                        );

        if (varcoeff==false)
          fcnonp[fcnonp.size()-1].init_name(terms[i].varnames[0]);
        else
          {
          vector<ST::string> na;
          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          fcnonp[fcnonp.size()-1].init_names(na);
          }

        fcnonp[fcnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcnonp[fcnonp.size()-1]);

        fcvarnonp.push_back(FULLCOND_variance_nonp(
                                       &generaloptions[generaloptions.size()-1],
                                         &fcnonp[fcnonp.size()-1],
                                         distr[distr.size()-1],
                                         a1,
                                         b1,
                                         titlev,pathnonpv,pathresv,
                                         false,collinpred)
                             );

        if (constlambda.getvalue() == true)
          fcvarnonp[fcvarnonp.size()-1].set_constlambda();

        if(terms[i].options[14]=="true")
          fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);


        }

      }

    }

  return false;
  }


void plotautocorrun(bayesreg & b)
  {

#if defined(BORLAND_OUTPUT_WINDOW)

  b.outerror("ERROR: method plotautocor is not available in this version\n");

#elif defined(__BUILDING_GNU)

  b.outerror("ERROR: method plotautocor is not available in this version\n");

  #elif defined(JAVA_OUTPUT_WINDOW)

  if (b.resultsyesno==true)
	 {
     if (b.posteriormode == false)
       {
       ST::string path = b.outfile.getvalue()+"_autocor" + ".raw";
       if (b.generaloptions[0].get_samplesize() < b.maxlag2.getvalue()*4)
         b.outerror("ERROR: samplesize too small\n");
       else
         {
         b.simobj.autocorr(b.maxlag2.getvalue(),path);

         ST::string graphname = "_" + b.name + "_graph";
         b.newcommands.push_back("graph " + graphname);

         ST::string datasetname = "_" + b.name + "_autocor";
         b.newcommands.push_back("dataset " + datasetname);
         b.newcommands.push_back(datasetname + ".infile using " + path);

         ST::string ot="";
         if (b.outfile3.getvalue().length() > 0)
           ot = ot + "outfile=\""+b.outfile3.getvalue() + "\" ";
         if (b.meanonly.getvalue() == true)
           ot = ot + " mean ";
         if (b.replaceautocor.getvalue() == true)
           ot = ot + " replace ";

//         b.out(ot);

         if (ot.length() == 0)
           b.newcommands.push_back(graphname + ".plotautocor using " + datasetname);
         else
           b.newcommands.push_back(graphname + ".plotautocor , " + ot + " using "
           + datasetname);

         b.newcommands.push_back("drop " + graphname + " " + datasetname);

         }
       }
     else
       b.outerror("ERROR: no MCMC simulation results\n");
	 }
  else
	 b.outerror("ERROR: no regression results\n");

#endif

  }


void getsamplerun(bayesreg & b)
  {
  if (b.resultsyesno == true)
    {
    if (b.posteriormode == false)
      {
      #if defined(JAVA_OUTPUT_WINDOW)

      b.simobj.get_samples(b.newcommands,b.outfile.getvalue() + "_");
      #else
      b.simobj.get_samples(b.outfile.getvalue() + "_");
      #endif
      }
    else
      b.outerror("ERROR: no MCMC simulation results\n");
    }
  else
    b.outerror("ERROR: no regression results\n");

  }


void bayesreg::describe(const optionlist & globaloptions)
  {
  statobject::describe(globaloptions);
  }


bool bayesreg::create_varcoeffbaseline(const unsigned & collinpred)
  {

  long h;
  unsigned min,max,degree,nrknots;
  double lambda,a1,b1;
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

      if (f==1)
        return true;


      // -------------end: reading options, term information -------------------

      //--------- creating path for samples and and results, creating title ----

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_logbaseline.raw","_logbaseline.res","_logbaseline");

      //----- end: creating path for samples and and results, creating title ---


      fcbaseline.push_back( pspline_baseline(&generaloptions[generaloptions.size()-1],
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
                                                beg
                                               )
                             );

      if (constlambda.getvalue() == true)
        fcbaseline[fcbaseline.size()-1].set_lambdaconst(lambda);

      vector<ST::string> na;
      na.push_back(terms[i].varnames[0]);
      na.push_back(terms[i].varnames[1]);
      fcbaseline[fcbaseline.size()-1].init_names(na);
      fcbaseline[fcbaseline.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fcbaseline[fcbaseline.size()-1]);

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_logbaseline_var.raw","_logbaseline_var.res","_logbaseline_variance");

      fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                              &fcbaseline[fcbaseline.size()-1],
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


bool bayesreg::create_multibaseline(const unsigned & collinpred)
  {

  ST::string proposal;

  long h;
  unsigned min,max,degree,nrknots;
  double lambda,a1,b1;
  bool ub, gl;
  int gridsize;
  int f;

  unsigned i;
  int j=0,k=0;
  for(i=0;i<terms.size();i++)
    {
    if ( baseline.checkvector(terms,i) == true )
      {

      // --------------- reading options, term information ---------------------

      MCMC::fieldtype type;
      type = MCMC::RW2;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[2]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[3]).strtolong(h);
      degree = unsigned(h);

      f = (terms[i].options[4]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[5]).strtodouble(lambda);

      f = (terms[i].options[6]).strtodouble(a1);

      f = (terms[i].options[7]).strtodouble(b1);

      if (terms[i].options[8] == "false")
        ub = false;
      else
        ub = true;

      f = (terms[i].options[9]).strtolong(h);
      gridsize = unsigned(h);

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant" && terms[i].options[14] == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      proposal = terms[i].options[11];



      if (f==1)
        return true;

      if (terms[i].options[13] != "" && begin.getvalue() != "")
        outerror("WARNING: begin variable specified twice");

      if (terms[i].options[13] != "")
        gl = false;
      else
        gl = true;

      if (terms[i].options[13] != "")
        k = terms[i].options[13].isinlist(modelvarnamesv);

      datamatrix beg;
      if (begin.getvalue() == "" && terms[i].options[13] == "")
        beg = datamatrix(1,1);
      else if (terms[i].options[13] != "")
        beg = D.getCol(k);
      else
        beg = D.getCol(begpos);

      datamatrix statemat;
      if (state.getvalue() == "")
        statemat = datamatrix(1,1);
      else
        statemat = D.getCol(statepos);

      // -------------end: reading options, term information -------------------

      //--------- creating path for samples and and results, creating title ----

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_logbaseline.raw","_logbaseline.res","_logbaseline");

      //----- end: creating path for samples and and results, creating title ---


      if(proposal == "cp")
        {
        fcmultibaseline.push_back( pspline_multibaseline(&generaloptions[generaloptions.size()-1],
                                                distr[distr.size()-1],
                                                fcconst_intercept,
                                                D.getCol(j),
                                                a1,
                                                b1,
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

        fcmultibaseline[fcmultibaseline.size()-1].init_name(terms[i].varnames[0]);
        fcmultibaseline[fcmultibaseline.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcmultibaseline[fcmultibaseline.size()-1]);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_logbaseline_var.raw","_logbaseline_var.res","_logbaseline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                &fcmultibaseline[fcmultibaseline.size()-1],
                                distr[distr.size()-1],a1,b1,
                                title,pathnonp,pathres,ub,collinpred)
                                );

        }
      else if(proposal == "iwls")
        {

        }

      if (constlambda.getvalue() == false)
        {
        if(terms[i].options[10]=="true")
          fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
        }
      }

    }

  return false;
  }

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

bool bayesreg::create_ridge(const unsigned & collinpred)
  {

  // helpvariables
  double tau2start;
  double shrinkagestart;
  double a_shrinkage;
  double b_shrinkage;
  bool shrinkagefix=false;
  double effectstart;
  double shrinkageweight;
  bool shrinkageadaptive=false;

  // Vektoren mit den Startwerten die uebergeben werden für adaptive shrinkkage  
  vector<double> variances_vec;
  vector<double> shrinkagestart_vec;    // for adaptive shrinkage (individual shrinkageparameters)
  vector<double> a_shrinkage_vec;       // for adaptive shrinkage
  vector<double> b_shrinkage_vec;       // for adaptive shrinkage
  vector<bool> shrinkagefix_vec;
  vector<double> effectstart_vec;
  vector<double> shrinkageweight_vec;   // for weighted shrinkage (weights for variances)
  vector<bool> shrinkageadaptive_vec;
  
  // Matrizen fuer Startwerte 
  datamatrix variances;    //inverse Varianzparameter lambdastart=tau^2
  datamatrix effects;      //Regressionskoeffizienten
  datamatrix startdata;    //alle Startwerte also mit Hyperparametern

  
 // Sollen die Sartwerte der Effekte verwendet werden. Bei use_effectstart="false" 
 // werden die Startwerte in mcmc_const via ::posteriormode() berechnet.  
  bool use_effectstart = true;
  bool external = false;    // stammen Hyperparameter/Startwerte aus externer Datei
  int termnr1 = -1;           // Zaehler fuer die erste Ridgeterm Nummer aus der derzeit der
                            // Startdataset eingelesen wird.
  unsigned ridgecount = 0;  // Zaehler fuer die Anzahl der Ridgeterme
  unsigned readline = 0;    // Zeilennummer die aus Dataset eingelesen wird
  dataobject * datap;       // pointer to datasetobject
  statobject * s;           // pointer to statobject
  
  // names of startdata
  list<ST::string> startnames;
  startnames.push_back("effect");
  startnames.push_back("tau2");
  startnames.push_back("shrinkage");
  startnames.push_back("weight");
  startnames.push_back("a");
  startnames.push_back("b");
  startnames.push_back("shrinkagefix");
  startnames.push_back("adaptive");
  
  unsigned i, j, f;

  vector<ST::string> varnames;

  bool check = false;
  bool isridge = true;

  vector<FULLCOND_const*> fc;


  for(i=0;i<terms.size();i++)
    {
    if (shrinkage.checkvector(terms,i) == true)
      {
      if(terms[i].options[0] == "ridge")
        {
        check = true;
        varnames.push_back(terms[i].varnames[0]);
        if(termnr1==-1) 
          termnr1 = i; 
        ridgecount = ridgecount + 1;                              
        
        // Optionen aus externer Datendatei auslesen
        if (terms[i].options[1] != "")
          {
          external = true;
          }      
        }
      }     
    }

  if(check)
    {
    if(external == false) // Startwerte aus Termen einlesen 
      {
      for(i=0;i<terms.size();i++)
        { 
        if (shrinkage.checkvector(terms,i) == true)
          {
          if(terms[i].options[0] == "ridge")
            {
            // Folgende Werte werden nur aus dem 1. Term gesetzt
            if (terms[termnr1].options[8] == "true")                   //Option aus letztem Term
              {   
              shrinkagefix = true;
              }
            else
              {
              shrinkagefix = false;
              }
            if (terms[termnr1].options[9] == "true")                   //Option aus letztem Term
              {   
              shrinkageadaptive = true;
              }
            else
              {
              shrinkageadaptive = false;
              }
            
            // Folgende Werte werden aus jedem Term gesetzt
            f = terms[i].options[2].strtodouble(effectstart);
            if(effectstart==1E8)
              {
              use_effectstart = false;
              }       
            f = terms[i].options[5].strtodouble(shrinkageweight);
            
            // Folgende Werte werden nur aus jedem Term gesetzt wenn adaptive gewählt
            // ansonsten werden die Werte aus dem 1. Term genommen
            if(shrinkageadaptive==false)
              {
              f = terms[termnr1].options[3].strtodouble(tau2start);    //Option aus 1. Term
              f = terms[termnr1].options[4].strtodouble(shrinkagestart); //Option aus 1. Term
              f = terms[termnr1].options[6].strtodouble(a_shrinkage);    //Option aus 1. Term
              f = terms[termnr1].options[7].strtodouble(b_shrinkage);    //Option aus 1. Term
              }
            if(shrinkageadaptive==true)
              {
              f = terms[i].options[3].strtodouble(tau2start);    //Option aus jedem Term
              f = terms[i].options[4].strtodouble(shrinkagestart); //Option aus jedem Term
              f = terms[i].options[6].strtodouble(a_shrinkage);    //Option aus jedem Term
              f = terms[i].options[7].strtodouble(b_shrinkage);    //Option aus jedem Term
              }
              
            // Vektoren der uebergebenen Optionen
            effectstart_vec.push_back(effectstart);
            variances_vec.push_back(tau2start);
            shrinkagestart_vec.push_back(shrinkagestart);
            shrinkageweight_vec.push_back(shrinkageweight); 
            a_shrinkage_vec.push_back(a_shrinkage);
            b_shrinkage_vec.push_back(b_shrinkage);
            shrinkagefix_vec.push_back(shrinkagefix);
            shrinkageadaptive_vec.push_back(shrinkageadaptive);
            }
          }
        }

      // Startwerte in Datenmatrix schreiben
      startdata = datamatrix(variances_vec.size(),startnames.size(),0);          
      for(readline=0; readline<variances_vec.size(); readline++)
        {
        startdata(readline,0) = effectstart_vec[readline];
        startdata(readline,1) = variances_vec[readline];
        startdata(readline,2) = shrinkagestart_vec[readline];
        startdata(readline,3) = shrinkageweight_vec[readline];
        startdata(readline,4) = a_shrinkage_vec[readline];
        startdata(readline,5) = b_shrinkage_vec[readline];
        startdata(readline,6) = shrinkagefix_vec[readline];
        startdata(readline,7) = shrinkageadaptive_vec[readline];
        }
      }

    
    if(external == true) // Startwerte aus externem Dataset einlesen 
      {
      
      int objpos = findstatobject(*statobj,terms[termnr1].options[1],"dataset");
      if (objpos >= 0)
        {
        s = statobj->at(objpos);
        datap = dynamic_cast<dataobject*>(s);
        if (datap->obs()==0 || datap->getVarnames().size()==0)
          {
          outerror("ERROR: dataset object " + terms[termnr1].options[1] + 
                   " does not contain any data or the length of rows differ\n");
          return true;
          }
        else if (datap->getVarnames().size()>8)
          {
          outerror("ERROR: dataset object " + terms[termnr1].options[1] + 
                   " contains more than seven variables\n");
          return true;
          }
        else if (datap->getVarnames().size()<8)
          {
          outerror("ERROR: dataset object " + terms[termnr1].options[1] + 
                   " contains less than six variables\n");
          return true;
          }
        else if (datap->getVarnames().size()==8)
          {
          outerror("NOTE: dataset " + terms[termnr1].options[1] + 
                   " with starting values for the ridge-variables is assigned in variable " + 
                   terms[termnr1].varnames[0] + "\n");
          }  
        }
      else
        {
        outerror("ERROR: dataset object with the starting values for ridge is not existing.\n"
                 "Check if variable -" + terms[termnr1].varnames[0] + "- contains the startdata option.\n");
        return true;
        }
  
      // create datamatrix
      startnames = datap->getVarnames();
      ST::string expr = "";
      datap->makematrix(startnames,startdata,expr);

      // column arrangement in matrix startdata:
      // effect | tau2 | shrinkage | weight | a | b | shrinkagefix | adaptive
  
      if(ridgecount==startdata.rows())
        {
        for(readline=0;readline<startdata.rows();readline++)
          {
          // Folgende Werte werden nur aus der 1. Zeile gesetzt
          if (startdata.get(0,6)==0.0)           //Option aus letzter Zeile
            {
            shrinkagefix = false;
            }
          if (startdata.get(0,6)==1.0)
            {
            shrinkagefix = true;
            }
          if (startdata.get(0,7)==0.0)           //Option aus letzter Zeile
            {
            shrinkageadaptive = false;
            }
          if (startdata.get(0,7)==1.0)
            {
            shrinkageadaptive = true;
            }
            
          // Folgende Werte werden aus jeder Zeile gesetzt
          effectstart = startdata.get(readline,0);
          if(effectstart==1E8)
            {
            use_effectstart = false;
            }    
          shrinkageweight = startdata.get(readline,3);
          
          // Folgende Werte werden nur aus jeder Zeile gesetzt wenn adaptive gewählt
          // ansonsten werden die Werte aus der 1. zeile genommen
          if(shrinkageadaptive==false)
            {
            tau2start = startdata.get(0,1);
            shrinkagestart = startdata.get(0,2);   //Option aus 1. Zeile
            a_shrinkage = startdata.get(0,4);      //Option aus 1. Zeile
            b_shrinkage = startdata.get(0,5);      //Option aus 1. Zeile
            }
          if(shrinkageadaptive==true)
            {
            tau2start = startdata.get(readline,1);      //Option aus jeder Zeile
            shrinkagestart = startdata.get(readline,2);   //Option aus jeder Zeile
            a_shrinkage = startdata.get(readline,4);      //Option aus jeder Zeile
            b_shrinkage = startdata.get(readline,5);      //Option aus jeder Zeile
            }
          
          // Vektoren der Startwerte die übergeben werden sollen
          effectstart_vec.push_back(effectstart); 
          variances_vec.push_back(tau2start);
          shrinkagestart_vec.push_back(shrinkagestart);
          shrinkageweight_vec.push_back(shrinkageweight); 
          a_shrinkage_vec.push_back(a_shrinkage);
          b_shrinkage_vec.push_back(b_shrinkage);
          shrinkagefix_vec.push_back(shrinkagefix);
          shrinkageadaptive_vec.push_back(shrinkageadaptive);
          }                                              
        }
        
        else
          {
          outerror("ERROR: Number of rows in dataset" + terms[termnr1].options[1] + 
                   " don't coincide with the number of ridge terms in the model formula \n");
          return true;          
          }
      }
    if(use_effectstart == false)
      {
      outerror("NOTE: Starting values of ridge effects are comuted as posteriormode \n"); 
      }
    if(use_effectstart == true)
      {
      outerror("NOTE: Starting values of ridge effects are assigned by the user \n"); 
      }      
      

    // Cut-vektor zum identifizieren der Bloecke
    unsigned nr = varnames.size();
    unsigned bs = blocksize.getvalue();
    unsigned nrblocks = 1;
    vector<unsigned> cut;
    cut.push_back(0);
    i = bs;
    while(i<nr)
      {
      cut.push_back(i);
      i += bs;
      nrblocks++;
      }
    cut.push_back(nr);

    // Varianzparameter und Startwerte fuer Regressionskoeffizienten
    variances = datamatrix(variances_vec.size(),1,0);
    effects = datamatrix(variances_vec.size(),1,0); 

    for(i=0; i<variances_vec.size(); i++)
      {
      variances(i,0) = variances_vec[i];
      effects(i,0) = effectstart_vec[i];
      }
      
    // Daten
    datamatrix data(D.rows(),varnames.size(),0);

    for(i=0; i<varnames.size(); i++)
      {
      j = varnames[i].isinlist(modelvarnamesv);
      data.putCol(i, D.getCol(j));
      }

     // Titel und Pfade zur Datenspeicherung
    ST::string title, titlehelp;
    ST::string pathconst;
    ST::string pathconstres;

    // keine Intercept
    int constpos=-1;

    vector<ST::string> varnameshelp;

    // Case: Gaussian
    if ( check_gaussian(collinpred))
      {
      for(i=0; i<nrblocks; i++)
        {
        // Erstellen der Dateien fuer die Ergebnisse der Koeffizientenschaetzung
        varnameshelp = vector<ST::string>();
        for(j=cut[i]; j<cut[i+1]; j++)
          varnameshelp.push_back(varnames[j]);

        title = "ridge_Effects" + ST::inttostring(i+1);
        #if defined(__BUILDING_LINUX)
        pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #else
        pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #endif
        pathconstres = outfile.getvalue() + add_name + "_" + title + ".res";

        if (pathconst.isvalidfile() == 1)
          {
          errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
          return true;
          }

        // Uebergabe der Optionen an Constuctor FULLCOND_const_gaussian
        normalshrinkage.push_back(FULLCOND_const_gaussian(&generaloptions[generaloptions.size()-1],
                        distr[distr.size()-1], data.getColBlock(cut[i], cut[i+1]), 
                        title, constpos, pathconst, pathconstres, 
                        true, variances.getRowBlock(cut[i], cut[i+1]), 
                        use_effectstart, effects.getRowBlock(cut[i], cut[i+1]),
                        collinpred));

        normalshrinkage[normalshrinkage.size()-1].init_names(varnameshelp);
        normalshrinkage[normalshrinkage.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&normalshrinkage[normalshrinkage.size()-1]);
        fc.push_back(&normalshrinkage[normalshrinkage.size()-1]);
        }

      // Erstellen der Dateien fuer die Ergebnisse der Varianzparameterschaetzung
      title = "ridge";
      make_paths(collinpred,pathnonp,pathres,title,title,"","_var.raw","_var.res","_variance");

      distr[distr.size()-1]->set_ridge(data.cols());
      distr[distr.size()-1]->update_ridge(0.0);


      // Uebergabe der Optionen an Constuctor FULLCOND_variance_nonp_vector
      fcvarnonpvec.push_back(FULLCOND_variance_nonp_vector(&generaloptions[generaloptions.size()-1],
                   fc,distr[distr.size()-1],title,pathnonp,pathres,
                   shrinkagestart_vec,a_shrinkage_vec,b_shrinkage_vec,
                   shrinkagefix_vec,shrinkageweight_vec,
//                   startdata,
                   shrinkageadaptive_vec,
                   isridge,cut,collinpred));

      fullcond.push_back(&fcvarnonpvec[fcvarnonpvec.size()-1]);
      }


    // Case: NonGaussian
    else
      {

      // Erstellen der Dateien fuer die Ergebnisse der Koeffizientenschaetzung
      for(i=0; i<nrblocks; i++)
        {
        varnameshelp = vector<ST::string>();
        for(j=cut[i]; j<cut[i+1]; j++)
          varnameshelp.push_back(varnames[j]);

        title = "ridge_Effects" + ST::inttostring(i+1);
        #if defined(__BUILDING_LINUX)
        pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #else
        pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #endif
        pathconstres = outfile.getvalue() + add_name + "_" + title + ".res";
        if (pathconst.isvalidfile() == 1)
          {
          errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
          return true;
          }

        // Uebergabe der Optionen an Constuctor FULLCOND_const_nongaussian
        nongaussianshrinkage.push_back(FULLCOND_const_nongaussian(&generaloptions[generaloptions.size()-1],
                             distr[distr.size()-1], data.getColBlock(cut[i], cut[i+1]), title, constpos,
                             pathconst, pathconstres, true, variances.getRowBlock(cut[i], cut[i+1]), 
                             use_effectstart, effects.getRowBlock(cut[i], cut[i+1]),
                             collinpred));
                             
        nongaussianshrinkage[nongaussianshrinkage.size()-1].init_names(varnameshelp);
        nongaussianshrinkage[nongaussianshrinkage.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&nongaussianshrinkage[nongaussianshrinkage.size()-1]);
        fc.push_back(&nongaussianshrinkage[nongaussianshrinkage.size()-1]);
        }

      // Erstellen der Dateien fuer die Ergebnisse der Varianzparameterschaetzung
      title = "ridge";
      make_paths(collinpred,pathnonp,pathres,title,title,"","_var.raw","_var.res","_variance");

     // Uebergabe der Optionen an Constuctor FULLCOND_variance_nonp_vector
      fcvarnonpvec.push_back(FULLCOND_variance_nonp_vector(&generaloptions[generaloptions.size()-1],
                   fc,distr[distr.size()-1],title,pathnonp,pathres,
                   shrinkagestart_vec,a_shrinkage_vec,b_shrinkage_vec,
                   shrinkagefix_vec,shrinkageweight_vec,
//                   startdata,
                   shrinkageadaptive_vec,
                   isridge,cut,collinpred));
      fullcond.push_back(&fcvarnonpvec[fcvarnonpvec.size()-1]);
//      fullcond.push_back(fcvarnonpvec[fcvarnonpvec.size()-1].get_shrinkagepointer());
      }
    }
  return false;
  }

bool bayesreg::create_lasso(const unsigned & collinpred)
  {

  // helpvariables
  double tau2start;
  double shrinkagestart;
  double a_shrinkage;
  double b_shrinkage;
  bool shrinkagefix;
  double effectstart;
  double shrinkageweight;
  bool shrinkageadaptive;

  // Vektoren mit den Startwerten die uebergeben werden  
  vector<double> variances_vec;
  vector<double> shrinkagestart_vec;    // for adaptive shrinkage (individual shrinkageparameters)
  vector<double> a_shrinkage_vec;       // for adaptive shrinkage
  vector<double> b_shrinkage_vec;       // for adaptive shrinkage
  vector<bool> shrinkagefix_vec;
  vector<double> effectstart_vec;
  vector<double> shrinkageweight_vec;   // for weighted shrinkage (weights for variances)
  vector<bool> shrinkageadaptive_vec;

  // Matrizen fuer Startwerte 
  datamatrix variances;    //inverse Varianzparameter lambdastart=tau^2
  datamatrix effects;      //Regressionskoeffizienten
  datamatrix startdata;    //alle Startwerte also mit Hyperparametern


  // Sollen die Sartwerte der Effekte verwendet werden. Bei use_effectstart="false"
  // werden die Startwerte in mcmc_const via ::posteriormode() berechnet. 
  bool use_effectstart = true;
  bool external = false;    // stammen Hyperparameter/Startwerte aus externer Datei
  int termnr1 = -1;           // Zaehler fuer die erste Lassoterm Nummer aus der derzeit der
                            // Startdataset üebergeben wird.
  unsigned lassocount = 0;  // Zaehler fuer die Anzahl der Lassoterme
  unsigned readline = 0;    // Zeilennummer die aus Dataset eingelesen wird
  dataobject * datap;       // pointer to datasetobject
  statobject * s;           // pointer to statobject
  
 // names of startdata
  list<ST::string> startnames;
  startnames.push_back("effect");
  startnames.push_back("tau2");
  startnames.push_back("shrinkage");
  startnames.push_back("weight");
  startnames.push_back("a");
  startnames.push_back("b");
  startnames.push_back("shrinkagefix");
  startnames.push_back("adaptive");

  unsigned i, j, f;

  vector<ST::string> varnames;

  bool check = false;
  bool isridge = false;

  vector<FULLCOND_const*> fc;


  for(i=0;i<terms.size();i++)
    {

    if ( shrinkage.checkvector(terms,i) == true )
      {
      if(terms[i].options[0] == "lasso")
        {
        check = true;
        varnames.push_back(terms[i].varnames[0]);
        if(termnr1==-1) 
          termnr1 = i; 
        lassocount = lassocount + 1;                              
        
        // Optionen aus externer Datendatei auslesen
        if (terms[i].options[1] != "")
          {
          external = true;
          }  
        }
      }
    }

  if(check)
    {
    if(external == false) // Startwerte aus Termen einlesen 
      {
      for(i=0;i<terms.size();i++)
        { 
        if (shrinkage.checkvector(terms,i) == true)
          {
          if(terms[i].options[0] == "lasso")
            {
            // Folgende Werte werden nur aus dem 1. Term gesetzt
            if (terms[termnr1].options[8] == "true")                   //Option aus letztem Term
              {   
              shrinkagefix = true;
              }
            else
              {
              shrinkagefix = false;
              }
            if (terms[termnr1].options[9] == "true")                   //Option aus letztem Term
              {   
              shrinkageadaptive = true;
              }
            else
              {
              shrinkageadaptive = false;
              }
            
            // Folgende Werte werden aus jedem Term gesetzt
            f = terms[i].options[2].strtodouble(effectstart);
            if(effectstart==1E8)
              {
              use_effectstart = false;
              }       
            f = terms[i].options[3].strtodouble(tau2start); 
            f = terms[i].options[5].strtodouble(shrinkageweight);
            
            // Folgende Werte werden nur aus jedem Term gesetzt wenn adaptive gewählt
            // ansonsten werden die Werte aus dem 1. Term genommen
            if(shrinkageadaptive==false)
              {
              f = terms[termnr1].options[4].strtodouble(shrinkagestart); //Option aus 1. Term
              f = terms[termnr1].options[6].strtodouble(a_shrinkage);    //Option aus 1. Term
              f = terms[termnr1].options[7].strtodouble(b_shrinkage);    //Option aus 1. Term
              }
            if(shrinkageadaptive==true)
              {
              f = terms[i].options[4].strtodouble(shrinkagestart); //Option aus jedem Term
              f = terms[i].options[6].strtodouble(a_shrinkage);    //Option aus jedem Term
              f = terms[i].options[7].strtodouble(b_shrinkage);    //Option aus jedem Term
              }
              
            // Vektoren der uebergebenen Optionen
            effectstart_vec.push_back(effectstart);
            variances_vec.push_back(tau2start);
            shrinkagestart_vec.push_back(shrinkagestart);
            shrinkageweight_vec.push_back(shrinkageweight); 
            a_shrinkage_vec.push_back(a_shrinkage);
            b_shrinkage_vec.push_back(b_shrinkage);
            shrinkagefix_vec.push_back(shrinkagefix);
            shrinkageadaptive_vec.push_back(shrinkageadaptive);
            }
          }
        }

      // Startwerte in Datenmatrix schreiben
      startdata = datamatrix(variances_vec.size(),startnames.size(),0);          
      for(readline=0; readline<variances_vec.size(); readline++)
        {
        startdata(readline,0) = effectstart_vec[readline];
        startdata(readline,1) = variances_vec[readline];
        startdata(readline,2) = shrinkagestart_vec[readline];
        startdata(readline,3) = shrinkageweight_vec[readline];
        startdata(readline,4) = a_shrinkage_vec[readline];
        startdata(readline,5) = b_shrinkage_vec[readline];
        startdata(readline,6) = shrinkagefix_vec[readline];
        startdata(readline,7) = shrinkageadaptive_vec[readline];
        }
      }

    if(external == true) // Startwerte aus externem Dataset einlesen 
      {
      
      int objpos = findstatobject(*statobj,terms[termnr1].options[1],"dataset");
      if (objpos >= 0)
        {
        s = statobj->at(objpos);
        datap = dynamic_cast<dataobject*>(s);
        if (datap->obs()==0 || datap->getVarnames().size()==0)
          {
          outerror("ERROR: dataset object " + terms[termnr1].options[1] + 
                   " does not contain any data or the length of rows differ\n");
          return true;
          }
        else if (datap->getVarnames().size()>8)
          {
          outerror("ERROR: dataset object " + terms[termnr1].options[1] + 
                   " contains more than seven variables\n");
          return true;
          }
        else if (datap->getVarnames().size()<8)
          {
          outerror("ERROR: dataset object " + terms[termnr1].options[1] + 
                   " contains less than six variables\n");
          return true;
          }
        else if (datap->getVarnames().size()==8)
          {
          outerror("NOTE: dataset " + terms[termnr1].options[1] + 
                   " with starting values for the lasso-variables is assigned in variable " + 
                   terms[termnr1].varnames[0] + "\n");
          }  
        }
      else
        {
        outerror("ERROR: dataset object with the starting values for lasso is not existing.\n"
                 "Check if variable -" + terms[termnr1].varnames[0] + "- contains the startdata option.\n");
        return true;
        }
  
      // create datamatrix
      startnames = datap->getVarnames(); 
      ST::string expr = "";
      datap->makematrix(startnames,startdata,expr);

      // column arrangement in matrix startdata:
      // effect | tau2 | shrinkage | weight | a | b | shrinkagefix | adaptive

      if(lassocount==startdata.rows())
        {
        for(readline=0;readline<startdata.rows();readline++)
          {
          // Folgende Werte werden nur aus der 1. Zeile gesetzt
          if (startdata.get(0,6)==0.0)           //Option aus letzter Zeile
            {
            shrinkagefix = false;
            }
          if (startdata.get(0,6)==1.0)
            {
            shrinkagefix = true;
            }
          if (startdata.get(0,7)==0.0)           //Option aus letzter Zeile
            {
            shrinkageadaptive = false;
            }
          if (startdata.get(0,7)==1.0)
            {
            shrinkageadaptive = true;
            }
            
          // Folgende Werte werden aus jeder Zeile gesetzt
          effectstart = startdata.get(readline,0);
          if(effectstart==1E8)
            {
            use_effectstart = false;
            }    
          tau2start = startdata.get(readline,1);
          shrinkageweight = startdata.get(readline,3);
          
          // Folgende Werte werden nur aus jeder Zeile gesetzt wenn adaptive gewählt
          // ansonsten werden die Werte aus der 1. zeile genommen
          if(shrinkageadaptive==false)
            {
            shrinkagestart = startdata.get(0,2);   //Option aus 1. Zeile
            a_shrinkage = startdata.get(0,4);      //Option aus 1. Zeile
            b_shrinkage = startdata.get(0,5);      //Option aus 1. Zeile
            }
          if(shrinkageadaptive==true)
            {
            shrinkagestart = startdata.get(readline,2);   //Option aus jeder Zeile
            a_shrinkage = startdata.get(readline,4);      //Option aus jeder Zeile
            b_shrinkage = startdata.get(readline,5);      //Option aus jeder Zeile
            }
          
          // Vektoren der Startwerte die übergeben werden sollen
          effectstart_vec.push_back(effectstart); 
          variances_vec.push_back(tau2start);
          shrinkagestart_vec.push_back(shrinkagestart);
          shrinkageweight_vec.push_back(shrinkageweight); 
          a_shrinkage_vec.push_back(a_shrinkage);
          b_shrinkage_vec.push_back(b_shrinkage);
          shrinkagefix_vec.push_back(shrinkagefix);
          shrinkageadaptive_vec.push_back(shrinkageadaptive);   
          }                                              
        }
        
        else
          {
          outerror("ERROR: Number of rows in dataset" + terms[termnr1].options[1] + 
                   " don't coincide with the number of lasso terms in the model formula \n");
          return true;          
          }
      }
    if(use_effectstart == false)
      {
      outerror("NOTE: Starting values of lasso effects are comuted as posteriormode \n"); 
      }
    if(use_effectstart == true)
      {
      outerror("NOTE: Starting values of lasso effects are assigned by the user \n"); 
      }      
      

    // Cut-vektor zum identifizieren der Bloecke
    unsigned nr = varnames.size();
    unsigned bs = blocksize.getvalue();
    unsigned nrblocks = 1;
    vector<unsigned> cut;
    cut.push_back(0);
    i = bs;
    while(i<nr)
      {
      cut.push_back(i);
      i += bs;
      nrblocks++;
      }
    cut.push_back(nr);

    // Varianzparameter und Startwerte fuer Regressionskoeffizienten
    variances = datamatrix(variances_vec.size(),1,0);
    effects = datamatrix(variances_vec.size(),1,0); 

    for(i=0; i<variances_vec.size(); i++)
      {
      variances(i,0) = variances_vec[i];
      effects(i,0) = effectstart_vec[i];
      }

    // Daten
    datamatrix data(D.rows(),varnames.size(),0);

    for(i=0; i<varnames.size(); i++)
      {
      j = varnames[i].isinlist(modelvarnamesv);
      data.putCol(i, D.getCol(j));
      }

     // Titel und Pfade zur Datenspeicherung
    ST::string title, titlehelp;
    ST::string pathconst;
    ST::string pathconstres;

    // keine Intercept
    int constpos=-1;

    vector<ST::string> varnameshelp;

    // Case: Gaussian
    if ( check_gaussian(collinpred))
      {
      for(i=0; i<nrblocks; i++)
        {
        // Erstellen der Dateien fuer die Ergebnisse der Koeffizientenschaetzung
        varnameshelp = vector<ST::string>();
        for(j=cut[i]; j<cut[i+1]; j++)
          varnameshelp.push_back(varnames[j]);

        title = "lasso_Effects" + ST::inttostring(i+1);
        #if defined(__BUILDING_LINUX)
        pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #else
        pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #endif
        pathconstres = outfile.getvalue() + add_name + "_" + title + ".res";

        if (pathconst.isvalidfile() == 1)
          {
          errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
          return true;
          }

        // Uebergabe der Optionen an Constuctor FULLCOND_const_gaussian
        normalshrinkage.push_back(FULLCOND_const_gaussian(&generaloptions[generaloptions.size()-1],
                        distr[distr.size()-1], data.getColBlock(cut[i], cut[i+1]),
                        title, constpos, pathconst, pathconstres, 
                        true, variances.getRowBlock(cut[i], cut[i+1]), 
                        use_effectstart, effects.getRowBlock(cut[i], cut[i+1]),
                        collinpred));

        normalshrinkage[normalshrinkage.size()-1].init_names(varnameshelp);
        normalshrinkage[normalshrinkage.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&normalshrinkage[normalshrinkage.size()-1]);
        fc.push_back(&normalshrinkage[normalshrinkage.size()-1]);
        }

      // Erstellen der Dateien fuer die Ergebnisse der Varianzparameterschaetzung
      title = "lasso";
      make_paths(collinpred,pathnonp,pathres,title,title,"","_var.raw","_var.res","_variance");

      distr[distr.size()-1]->set_lasso(data.cols());
      distr[distr.size()-1]->update_lasso(0.0);

      // Uebergabe der Optionen an Constuctor FULLCOND_variance_nonp_vector
      fcvarnonpvec.push_back(FULLCOND_variance_nonp_vector(&generaloptions[generaloptions.size()-1],
                   fc,distr[distr.size()-1],title,pathnonp,pathres,
                   shrinkagestart_vec,a_shrinkage_vec,b_shrinkage_vec,
                   shrinkagefix_vec,shrinkageweight_vec,
//                   startdata,
                   shrinkageadaptive_vec,
                   isridge,cut,collinpred));

      fullcond.push_back(&fcvarnonpvec[fcvarnonpvec.size()-1]);
      }


    // Case: NonGaussian
    else
      {

      // Erstellen der Dateien fuer die Ergebnisse der Koeffizientenschaetzung
      for(i=0; i<nrblocks; i++)
        {
        varnameshelp = vector<ST::string>();
        for(j=cut[i]; j<cut[i+1]; j++)
          varnameshelp.push_back(varnames[j]);

        title = "lasso_Effects" + ST::inttostring(i+1);
        #if defined(__BUILDING_LINUX)
        pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #else
        pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #endif
        pathconstres = outfile.getvalue() + add_name + "_" + title + ".res";
        if (pathconst.isvalidfile() == 1)
          {
          errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
          return true;
          }

        // Uebergabe der Optionen an Constuctor FULLCOND_const_nongaussian
        nongaussianshrinkage.push_back(FULLCOND_const_nongaussian(&generaloptions[generaloptions.size()-1],
                              distr[distr.size()-1], data.getColBlock(cut[i], cut[i+1]), title, constpos,
                              pathconst, pathconstres, true, variances.getRowBlock(cut[i], cut[i+1]), 
                              use_effectstart, effects.getRowBlock(cut[i], cut[i+1]),
                              collinpred));

        nongaussianshrinkage[nongaussianshrinkage.size()-1].init_names(varnameshelp);
        nongaussianshrinkage[nongaussianshrinkage.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&nongaussianshrinkage[nongaussianshrinkage.size()-1]);
        fc.push_back(&nongaussianshrinkage[nongaussianshrinkage.size()-1]);
        }

      // Erstellen der Dateien fuer die Ergebnisse der Varianzparameterschaetzung
      title = "lasso";
      make_paths(collinpred,pathnonp,pathres,title,title,"","_var.raw","_var.res","_variance");

     // Uebergabe der Optionen an Constuctor FULLCOND_variance_nonp_vector
      fcvarnonpvec.push_back(FULLCOND_variance_nonp_vector(&generaloptions[generaloptions.size()-1],
                   fc,distr[distr.size()-1],title,pathnonp,pathres,
                   shrinkagestart_vec,a_shrinkage_vec,b_shrinkage_vec,
                   shrinkagefix_vec,shrinkageweight_vec,
//                   startdata,
                   shrinkageadaptive_vec,
                   isridge,cut,collinpred));
      fullcond.push_back(&fcvarnonpvec[fcvarnonpvec.size()-1]);
//      fullcond.push_back(fcvarnonpvec[fcvarnonpvec.size()-1].get_shrinkagepointer());
      }
    }

  return false;
  }


bool bayesreg::create_nigmix(const unsigned & collinpred)
  {

  // helpvariables
  double effectstart;
  unsigned long indicator;
  double v0;
  double v1;
  double t2;
  double a_t2;
  double b_t2;
  double a_omega;
  double b_omega;
  double omega;
  bool omegafix=false;
  bool omegaadaptive=false;
  
  double helpindicator;
  double helpvariances;
 
  // Vector-options (ersetzen später die 1-dim optionen von oben) 
  vector<double> variances_vec;
  vector<double> effectstart_vec;
  vector<unsigned long> indicator_vec;
  vector<double> v0_vec;
  vector<double> v1_vec;
  vector<double> t2_vec;
  vector<double> a_t2_vec;
  vector<double> b_t2_vec;
  vector<double> a_omega_vec;
  vector<double> b_omega_vec;
  vector<double> omega_vec;
  vector<bool> omegafix_vec;
  vector<bool> omegaadaptive_vec;



  // Matrizen fuer Startwerte 
  datamatrix variances;    //inverse Varianzparameter tau^2
  datamatrix effects;      //Regressionskoeffizienten
  datamatrix startdata;    //alle Startwerte also mit Hyperparametern
  
  // Sollen die Sartwerte der Effekte verwendet werden. Bei use_effectstart="false"
  // werden die Startwerte in mcmc_const via ::posteriormode() berechnet. 
  bool use_effectstart = true;
  bool external = false;    // stammen Hyperparameter/Startwerte aus externer Datei
  int termnr1 = -1;           // Zaehler fuer die erste NIGMIXterm Nummer aus der derzeit der
                            // Startdataset eingelesen werden.
  unsigned nigmixcount = 0; // Zaehler fuer die Anzahl der Nigmixterme
  unsigned readline = 0;    // Zeilennummer die aus Dataset eingelesen wird
  dataobject * datap;       // pointer to datasetobject
  statobject * s;           // pointer to statobject
  
 // names of startdata
  list<ST::string> startnames;
  startnames.push_back("effect");
  startnames.push_back("I");
  startnames.push_back("t2");
  startnames.push_back("w");
  startnames.push_back("v0");
  startnames.push_back("v1");
  startnames.push_back("a");
  startnames.push_back("b");
  startnames.push_back("aw");
  startnames.push_back("bw");
  startnames.push_back("wfix");
  startnames.push_back("adaptive");

  
  unsigned i, j, f;
  long help;
  
  vector<ST::string> varnames;

  bool check=false;
  
  vector<FULLCOND_const*> fc;


  for(i=0;i<terms.size();i++)
    {

    if ( nigmix.checkvector(terms,i) == true )
      {
      if(terms[i].options[0] == "nigmix")
        {
        check=true;
        varnames.push_back(terms[i].varnames[0]);
        if(termnr1==-1) 
          termnr1 = i; 
        nigmixcount = nigmixcount + 1;
        
        // Optionen aus externer Datendatei auslesen
        if (terms[i].options[1] != "")
          {
          external = true;
          }
       
        }
      }     
    }


  if(check)
    {
    if(external == false) // Startwerte aus Termen einlesen 
      {
      for(i=0;i<terms.size();i++)
        { 
        if (nigmix.checkvector(terms,i) == true)
          {
          if(terms[i].options[0] == "nigmix")
            {
            // Folgende Werte werden nur aus dem 1. Term gesetzt
            if (terms[termnr1].options[12] == "true")                  //Option aus 1. Term
              omegafix = true;
            else
              omegafix = false;
            if (terms[termnr1].options[13] == "true")                  //Option aus 1. Term
              omegaadaptive = true;
            else
              omegaadaptive = false;
              
            // Folgende Werte werden aus jedem Term gesetzt
            f = terms[i].options[2].strtodouble(effectstart);
            if(effectstart==1E8)
              {
              use_effectstart = false;
              }
            f = (terms[i].options[3]).strtolong(help);
            indicator = help;
            f = (terms[i].options[4]).strtodouble(t2);
            
            // Folgende Werte werden nur aus jedem Term gesetzt wenn adaptive gewählt
            // ansonsten werden die Werte aus dem 1. Term genommen
          if(omegaadaptive==false)
            {
            f = (terms[termnr1].options[5]).strtodouble(omega);              //Option aus 1. Term
            f = (terms[termnr1].options[6]).strtodouble(v0);                 //Option aus 1. Term
            f = (terms[termnr1].options[7]).strtodouble(v1);                 //Option aus 1. Term
            f = (terms[termnr1].options[8]).strtodouble(a_t2);               //Option aus 1. Term
            f = (terms[termnr1].options[9]).strtodouble(b_t2);               //Option aus 1. Term
            f = (terms[termnr1].options[10]).strtodouble(a_omega);           //Option aus 1. Term
            f = (terms[termnr1].options[11]).strtodouble(b_omega);           //Option aus 1. Term
            }
          if(omegaadaptive==true)
            {
            f = (terms[i].options[5]).strtodouble(omega);              //Option aus jedem Term
            f = (terms[i].options[6]).strtodouble(v0);                 //Option aus jedem Term
            f = (terms[i].options[7]).strtodouble(v1);                 //Option aus jedem Term
            f = (terms[i].options[8]).strtodouble(a_t2);               //Option aus jedem Term
            f = (terms[i].options[9]).strtodouble(b_t2);               //Option aus jedem Term
            f = (terms[i].options[10]).strtodouble(a_omega);           //Option aus jedem Term
            f = (terms[i].options[11]).strtodouble(b_omega);           //Option aus jedem Term
            }
            // Vektoren der uebergebenen Optionen
            if(indicator==0)
              {
              helpvariances = (v0 * t2);
              }
            if(indicator==1)
              {
              helpvariances = (v1 * t2);
              }
            variances_vec.push_back(helpvariances);
            effectstart_vec.push_back(effectstart);
            indicator_vec.push_back(indicator);
            t2_vec.push_back(t2);
            omega_vec.push_back(omega);
            v0_vec.push_back(v0);
            v1_vec.push_back(v1); 
            a_t2_vec.push_back(a_t2);
            b_t2_vec.push_back(b_t2);
            a_omega_vec.push_back(a_omega);
            b_omega_vec.push_back(b_omega);
            omegafix_vec.push_back(omegafix);  
            omegaadaptive_vec.push_back(omegaadaptive);              
            }
          }
        }


      // Startwerte in Datenmatrix schreiben
      startdata = datamatrix(variances_vec.size(),startnames.size(),0);          
      for(readline=0; readline<variances_vec.size(); readline++)
        {
        startdata(readline,0) = effectstart_vec[readline];
        startdata(readline,1) = indicator_vec[readline];
        startdata(readline,2) = t2_vec[readline];
        startdata(readline,3) = omega_vec[readline];
        startdata(readline,4) = v0_vec[readline];
        startdata(readline,5) = v1_vec[readline];
        startdata(readline,6) = a_t2_vec[readline];
        startdata(readline,7) = b_t2_vec[readline];
        startdata(readline,8) = a_omega_vec[readline];
        startdata(readline,9) = b_omega_vec[readline];
        startdata(readline,10) = omegafix_vec[readline];
        startdata(readline,11) = omegaadaptive_vec[readline];
        }
     }
    
    if(external == true) // Startwerte aus externem Dataset einlesen 
      {

      
      int objpos = findstatobject(*statobj,terms[termnr1].options[1],"dataset");
      if (objpos >= 0)
        {
        s = statobj->at(objpos);
        datap = dynamic_cast<dataobject*>(s);

        if (datap->obs()==0 || datap->getVarnames().size()==0)
          {
          outerror("ERROR: dataset object " + terms[termnr1].options[1] + 
                   " does not contain any data or the length of rows differ\n");
          return true;
          }
        else if (datap->getVarnames().size()>12)
          {
          outerror("ERROR: dataset object " + terms[termnr1].options[1] + 
                   " contains more than eleven variables\n");
          return true;
          }
        else if (datap->getVarnames().size()<12)
          {
          outerror("ERROR: dataset object " + terms[termnr1].options[1] + 
                   " contains less than eleven variables\n");
          return true;
          }
        else if (datap->getVarnames().size()==12)
          {
          outerror("NOTE: dataset " + terms[termnr1].options[1] + 
                   " with starting values for the nigmix-variables is assigned in variable " + 
                   terms[termnr1].varnames[0] + "\n");
          }  
        }
      else
        {
        outerror("ERROR: dataset object with the starting values for nigmix is not existing.\n"
                 "Check if variable -" + terms[termnr1].varnames[0] + "- contains the startdata option.\n");
        return true;
        }
  
      // create datamatrix
      startnames = datap->getVarnames();
      ST::string expr = "";
      datap->makematrix(startnames,startdata,expr);

      // column arrangement in matrix startdata:
      // effect | I | t2 | w | v0 | v1 | a | b | aw | bw| wfix | adaptive
   
      if(nigmixcount==startdata.rows())
        {
        for(readline=0;readline<startdata.rows();readline++)
          {
          // Folgende Werte werden nur aus der 1. Zeile gesetzt
          if (startdata.get(0,10)==0.0)                    //Option 1. Term
            {
            omegafix = false;
            }
          if (startdata.get(0,10)==1.0)
            {
            omegafix = true;
            }
          if (startdata.get(0,11)==0.0)                     //Option 1. Term
            {
            omegaadaptive = false;
            }
          if (startdata.get(0,11)==1.0)
            {
            omegaadaptive = true;
            }
          
          // Folgende Werte werden aus jeder Zeile gesetzt
          effectstart = startdata.get(readline,0);
          if(effectstart==1E8)
            {
            use_effectstart = false;
            }         
          
          helpindicator = startdata.get(readline,1);
          
          if(helpindicator!=1 && helpindicator!=0)
            {
            outerror("ERROR: Indicators need to be 0 or 1 \n"); 
            return true;
            }
          else
            indicator = startdata.get(readline,1);
            
          t2 = startdata.get(readline,2);
          
          // Folgende Werte werden nur aus jeder Zeile gesetzt wenn adaptive gewählt
          // ansonsten werden die Werte aus der 1. zeile genommen
          if(omegaadaptive==false)
            {
            omega = startdata.get(0,3);             //Option aus 1. Term
            v0 = startdata.get(0,4);                //Option aus 1. Term
            v1 = startdata.get(0,5);                //Option aus 1. Term
            a_t2 = startdata.get(0,6);              //Option aus 1. Term
            b_t2 = startdata.get(0,7);              //Option aus 1. Term
            a_omega = startdata.get(0,8);           //Option aus 1. Term
            b_omega = startdata.get(0,9);           //Option aus 1. Term
            }
          if(omegaadaptive==true)
            {
            omega = startdata.get(readline,3);             //Option aus jedem Term
            v0 = startdata.get(readline,4);                //Option aus jedem Term
            v1 = startdata.get(readline,5);                //Option aus jedem Term
            a_t2 = startdata.get(readline,6);              //Option aus jedem Term
            b_t2 = startdata.get(readline,7);              //Option aus jedem Term
            a_omega = startdata.get(readline,8);           //Option aus jedem Term
            b_omega = startdata.get(readline,9);           //Option aus jedem Term
            }

          if(indicator==0)
            {
            helpvariances = (v0 * t2);
            }
          if(indicator==1)
            {
            helpvariances = (v1 * t2);
            }
          // Vektoren denen Startwerte die übergeben werden sollen
          variances_vec.push_back(helpvariances);
          effectstart_vec.push_back(effectstart);
          indicator_vec.push_back(indicator);
          t2_vec.push_back(t2);
          omega_vec.push_back(omega);
          v0_vec.push_back(v0);
          v1_vec.push_back(v1); 
          a_t2_vec.push_back(a_t2);
          b_t2_vec.push_back(b_t2);
          a_omega_vec.push_back(a_omega);
          b_omega_vec.push_back(b_omega);
          omegafix_vec.push_back(omegafix); 
          omegaadaptive_vec.push_back(omegaadaptive);             
          }                                              
        }
        
        else
          {
          outerror("ERROR: Number of rows in dataset" + terms[termnr1].options[1] + 
                   " don't coincide with the number of nigmix terms in the model formula \n");
          return true;          
          }
      }
    if(use_effectstart == false)
      {
      outerror("NOTE: Starting values of nigmix effects are comuted as posteriormode \n"); 
      }
    if(use_effectstart == true)
      {
      outerror("NOTE: Starting values of nigmix effects are assigned by the user \n"); 
      }

    
    // Cut-vektor zum identifizieren der Bloecke
    unsigned nr = varnames.size();
    unsigned bs = blocksize.getvalue();
    unsigned nrblocks = 1;
    vector<unsigned> cut;
    cut.push_back(0);
    i = bs;
    while(i<nr)
      {
      cut.push_back(i);
      i += bs;
      nrblocks++;
      }
    cut.push_back(nr);

    // Varianzparameter und Startwerte fuer Regressionskoeffizienten
    variances = datamatrix(variances_vec.size(),1,0);
    effects = datamatrix(variances_vec.size(),1,0); 

    for(i=0; i<variances_vec.size(); i++)
      {
      variances(i,0) = variances_vec[i];
      effects(i,0) = effectstart_vec[i];
      }

    // Daten
    datamatrix data(D.rows(),varnames.size(),0);

    for(i=0; i<varnames.size(); i++)
      {
      j = varnames[i].isinlist(modelvarnamesv);
      data.putCol(i, D.getCol(j));
      }

     // Titel und Pfade zur Datenspeicherung
    ST::string title, titlehelp;
    ST::string pathconst;
    ST::string pathconstres;

    // keine Intercept
    int constpos=-1;

    vector<ST::string> varnameshelp;

    // Case: Gaussian
    if ( check_gaussian(collinpred))
      {
      for(i=0; i<nrblocks; i++)
        {
        // Erstellen der Dateien fuer die Ergebnisse der Koeffizientenschaetzung
        varnameshelp = vector<ST::string>();
        for(j=cut[i]; j<cut[i+1]; j++)
          varnameshelp.push_back(varnames[j]);

        title = "nigmix_Effects" + ST::inttostring(i+1);
        #if defined(__BUILDING_LINUX)
        pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #else
        pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #endif
        pathconstres = outfile.getvalue() + add_name + "_" + title + ".res";

        if (pathconst.isvalidfile() == 1)
          {
          errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
          return true;
          }

        // Uebergabe der Optionen an Constuctor FULLCOND_const_gaussian
        normalshrinkage.push_back(FULLCOND_const_gaussian(&generaloptions[generaloptions.size()-1],
                        distr[distr.size()-1], data.getColBlock(cut[i], cut[i+1]), 
                        title, constpos, pathconst, pathconstres, 
                        true, variances.getRowBlock(cut[i], cut[i+1]), 
                        use_effectstart, effects.getRowBlock(cut[i], cut[i+1]),
                        collinpred));

        normalshrinkage[normalshrinkage.size()-1].init_names(varnameshelp);
        normalshrinkage[normalshrinkage.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&normalshrinkage[normalshrinkage.size()-1]);
        fc.push_back(&normalshrinkage[normalshrinkage.size()-1]);
        }

      // Erstellen der Dateien fuer die Ergebnisse der Varianzparameterschaetzung
      title = "nigmix";
      make_paths(collinpred,pathnonp,pathres,title,title,"","_var.raw","_var.res","_variance");

      distr[distr.size()-1]->set_nigmix(data.cols());
      distr[distr.size()-1]->update_nigmix(0.0);


      // Uebergabe der Optionen an Constuctor FULLCOND_variance_nonp_vector_nigmix
      fcvarnonpvecnigmix.push_back(FULLCOND_variance_nonp_vector_nigmix(&generaloptions[generaloptions.size()-1],
                         fc,distr[distr.size()-1],title,pathnonp,pathres,
                         indicator_vec,v0_vec,v1_vec,t2_vec,a_t2_vec,b_t2_vec,
                         omega_vec,a_omega_vec,b_omega_vec,omegafix_vec,
//                         startdata,
                         omegaadaptive_vec,
                         cut,collinpred));

      fullcond.push_back(&fcvarnonpvecnigmix[fcvarnonpvecnigmix.size()-1]);
      }


    // Case: NonGaussian
    else
      {

      // Erstellen der Dateien fuer die Ergebnisse der Koeffizientenschaetzung
      for(i=0; i<nrblocks; i++)
        {
        varnameshelp = vector<ST::string>();
        for(j=cut[i]; j<cut[i+1]; j++)
          varnameshelp.push_back(varnames[j]);

        title = "nigmix_Effects" + ST::inttostring(i+1);
        #if defined(__BUILDING_LINUX)
        pathconst = defaultpath.to_bstr() + "/temp/" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #else
        pathconst = defaultpath.to_bstr() + "\\temp\\" + name.to_bstr()
                         + add_name + "_" + title + ".raw";
        #endif
        pathconstres = outfile.getvalue() + add_name + "_" + title + ".res";
        if (pathconst.isvalidfile() == 1)
          {
          errormessages.push_back("ERROR: unable to open file " + pathconst +
                                 " for writing\n");
          return true;
          }

        // Uebergabe der Optionen an Constuctor FULLCOND_const_nongaussian
        nongaussianshrinkage.push_back(FULLCOND_const_nongaussian(&generaloptions[generaloptions.size()-1],
                            distr[distr.size()-1], data.getColBlock(cut[i], cut[i+1]), title, constpos,
                            pathconst, pathconstres, true, variances.getRowBlock(cut[i], cut[i+1]), 
                            use_effectstart, effects.getRowBlock(cut[i], cut[i+1]), 
                            collinpred));
                            
        nongaussianshrinkage[nongaussianshrinkage.size()-1].init_names(varnameshelp);
        nongaussianshrinkage[nongaussianshrinkage.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&nongaussianshrinkage[nongaussianshrinkage.size()-1]);
        fc.push_back(&nongaussianshrinkage[nongaussianshrinkage.size()-1]);
        }

      // Erstellen der Dateien fuer die Ergebnisse der Varianzparameterschaetzung
      title = "nigmix";
      make_paths(collinpred,pathnonp,pathres,title,title,"","_var.raw","_var.res","_variance");

     // Uebergabe der Optionen an Constuctor FULLCOND_variance_nonp_vector_nigmix
      fcvarnonpvecnigmix.push_back(FULLCOND_variance_nonp_vector_nigmix(&generaloptions[generaloptions.size()-1],
                         fc,distr[distr.size()-1],title,pathnonp,pathres,
                         indicator_vec,v0_vec,v1_vec,t2_vec,a_t2_vec,b_t2_vec,
                         omega_vec,a_omega_vec,b_omega_vec,omegafix_vec,
//                         startdata,
                         omegaadaptive_vec,
                         cut,collinpred));
      fullcond.push_back(&fcvarnonpvecnigmix[fcvarnonpvecnigmix.size()-1]);
//      fullcond.push_back(fcvarnonpvec[fcvarnonpvec.size()-1].get_shrinkagepointer());
      }
    }
  return false;
  }


void regressrun(bayesreg & b)
  {

  vector<ST::string> header;

  b.resultsyesno = false;
  if (b.modeonly.getvalue() == true)
    b.posteriormode = true;
  else
    b.posteriormode = false;

  b.terms = b.modreg.getterms();

  b.describetext.erase(b.describetext.begin(),b.describetext.end());
  b.describetext.push_back("LAST ESTIMATED MODEL: \n");
  b.describetext.push_back("\n");
  b.describetext.push_back(b.modreg.getModelText());
  b.describetext.push_back("\n");

  if (b.varianceest==false && b.missingest==false)
    b.clear();


  if (b.family.getvalue()=="vargaussian")
    {
    b.varianceest=true;
    b.add_name="_variancereg";
    }
  else if (b.missingreg.getvalue()==true)
    {

    vector<ST::string> modelvarnamesv = b.modreg.getModelVarnamesAsVector();

    b.missingest=true;
    b.add_name="_missingreg_" + modelvarnamesv[0];
    }
  else
    b.add_name="";

  b.outfiles.push_back(b.outfile.getvalue()+b.add_name);


  bool failure = false;

  if (b.family.getvalue() == "multgaussian")
    {
    failure = true;
    b.out("ERROR: family multivariate gaussian is not allowed for method regress\n");
    }

  if (b.family.getvalue() == "multistate")
    {
    failure = true;
    b.out("ERROR: family multistate is not allowed for method regress\n");
    }

  if (!failure)
    failure = b.create_generaloptions();

  if (!failure)
    failure = b.create_distribution();

// Speicherplatz für normalconst/nongaussianconst/nbinomialconst reservieren
  unsigned nrfcfixed = b.fixedeffects.get_constvariables(b.terms).size()+1;
  unsigned blocksize_fixed = 10;
  unsigned reserved = 20;
  if ( nrfcfixed*b.nrcategories > reserved*blocksize_fixed )
    {
    nrfcfixed = ceil(nrfcfixed/double(blocksize_fixed));
    b.normalconst.reserve(nrfcfixed*b.nrcategories);
    b.nongaussianconst.reserve(nrfcfixed*b.nrcategories);
    b.nbinomialconst.reserve(nrfcfixed*b.nrcategories);
    }

  unsigned i;

  if (!failure)
    {
    for (i=0;i<b.nrcategories;i++)
      {

      if (!failure)
        failure = b.create_const(i);

      if (!failure)
        failure = b.create_baseline(i);

      if (!failure)
        failure = b.create_varcoeffbaseline(i);

      if (!failure)
        failure = b.create_nonprw1rw2(i);

      if (!failure)
        failure = b.create_pspline(i);

      if (!failure)
        failure = b.create_nonpseason(i);

      if (!failure)
        failure = b.create_spatial(i);

      if (!failure)
        failure = b.create_geospline(i);

      if (!failure)
        failure = b.create_varcoeff_geospline(i);

      if (!failure)
        failure = b.create_spatialxy(i);

      if (!failure)
        failure = b.create_varcoeffpspline(i);

      if (!failure)
        failure = b.create_random(i);

      if (!failure)
        failure = b.create_randomslope(i);

      if (!failure)
        failure = b.create_mixture(i);

      if (!failure)
        failure = b.create_interactionspspline(i);

      if (!failure)
        failure = b.create_geokriging(i);

      if(!failure)
        failure = b.create_varcoeffmerror(i);

     if(!failure)
        failure = b.create_ridge(i);

     if(!failure)
        failure = b.create_lasso(i);

     if(!failure)
        failure = b.create_nigmix(i);

      if (!failure)
        failure = b.create_random_rw1rw2(i);

      if (!failure)
        failure = b.create_spatial_rw1rw2(i);

      if (!failure)
        failure = b.create_random_pspline(i);

      } // end: for (i=0;i<b.nrcategories;i++)
    } // end: if (!failure)


  if(b.fcbaseline.size()>1)
    {
    vector<MCMC::pspline_baseline*> basep;
    for(i=0;i<b.fcbaseline.size();i++)
      basep.push_back(&b.fcbaseline[i]);
    for(i=0;i<b.fcbaseline.size();i++)
      b.fcbaseline[i].set_baselinep(basep);
    }

/*  if(b.fcbaselineiwls.size()>1 ||
     (b.fcbaselineiwls.size()>0 && b.fcbaseline.size()>0))
    {
    failure = true;
    b.out("ERROR: geht nicht!\n");
    }*/

  if (!failure                          &&
      (b.family.getvalue() != "vargaussian") &&
      (b.missingreg.getvalue()==false)
     )
    {

    if (b.varianceest==true)
      {

      vector<unsigned> begin;
      vector<unsigned> end;
      begin.push_back(0);
      end.push_back(b.varianceend_fc);
      begin.push_back(b.varianceend_fc+1);
      end.push_back(b.fullcond.size()-1);

      b.generaloptions[0].set_nrout(b.generaloptions[0].get_iterations()+1);

      vector<MCMCoptions*> mo;
      mo.push_back(&(b.generaloptions[0]));  // variance regression
      mo.push_back(&(b.generaloptions[1]));  // regression

      header.push_back("BAYESREG OBJECT " + b.name.to_bstr() +
                       ": variance regression" );
      header.push_back("BAYESREG OBJECT " + b.name.to_bstr() +
                       ": regression procedure" );

      b.simobj = MCMCsimulate(mo,b.distr,b.fullcond,begin,end);

      if (b.modeonly.getvalue())
        {
        failure = b.simobj.posteriormode(header,false);
        }
      else
        failure = b.simobj.simulate(header,b.setseed.getvalue(),!b.noposteriormode.getvalue());

      b.varianceest=false;
      }
    else if (b.missingest==true)
      {

      b.begin_fc.push_back(b.missingend_fc+1);
      b.end_fc.push_back(b.fullcond.size()-1);

      vector<MCMCoptions*> mo;

      unsigned i;
      for (i=0;i<b.generaloptions.size()-1;i++)
        {
        b.generaloptions[i].set_nrout(b.generaloptions[i].get_iterations()+1);
        mo.push_back(&(b.generaloptions[i]));  // missing regressions
        header.push_back("BAYESREG OBJECT " + b.name.to_bstr() +
                       ": missing value regression (" + b.distr[i]->get_responsename() + ")");

        }

      mo.push_back(&(b.generaloptions[b.generaloptions.size()-1]));// regression


      header.push_back("BAYESREG OBJECT " + b.name.to_bstr() +
                       ": regression procedure" );

      unsigned regpos = b.generaloptions.size()-1;

      for (i=0;i<b.generaloptions.size()-1;i++)
        {

#if defined(__BUILDING_LINUX)
        ST::string pathtemp = b.defaultpath + "/temp/" + b.name + "_missingreg_"
                              + b.distr[i]->get_responsename() + ".raw";
#else
        ST::string pathtemp = b.defaultpath + "\\temp\\" + b.name + "_missingreg_"
                              + b.distr[i]->get_responsename() + ".raw";
#endif
        ST::string pathres = b.outfiles[i] + "_missing_" +
                           b.distr[i]->get_responsename() + ".res";

        b.distr[i]->set_missings(b.fullcond,b.begin_fc[regpos],b.end_fc[regpos],
                                 b.mind[i],pathtemp,pathres);
        }


//      b.simobj = MCMCsimulate(mo,b.distr,b.fullcond,b.begin_fc,b.end_fc);


      if (b.modeonly.getvalue())
        {
        failure = b.simobj.posteriormode(header,false);
        }
      else
        {
        if (b.family.getvalue() == "cumprobit")
          {
          failure = b.simobj.simulate(header,b.setseed.getvalue(),false);
          }
        else
          failure = b.simobj.simulate(header,b.setseed.getvalue(),!b.noposteriormode.getvalue());
        }

      b.missingest=false;

      }
    else
      {

      header.push_back("BAYESREG OBJECT " + b.name.to_bstr() +
                       ": regression procedure" );

      b.simobj = MCMCsimulate(&b.generaloptions[0],b.distr[0],b.fullcond);
      if (b.modeonly.getvalue())
        {
        failure = b.simobj.posteriormode(header);
        }
      else
        {
        if (b.nosamples.getvalue() == true)
          b.simobj.setflags(MCMC::nosamples);
        if ( (b.family.getvalue() == "cumprobit") ||
             (b.family.getvalue() == "multinomialprobit") ||
             (b.family.getvalue() == "binomialtlink")
           )
           {
           failure = b.simobj.simulate(header,b.setseed.getvalue(),false);
           }
        else
          {
          failure = b.simobj.simulate(header,b.setseed.getvalue(),!b.noposteriormode.getvalue());
          }

        }

      }

    }


  if (!failure && (b.family.getvalue() != "vargaussian") &&
     (b.missingreg.getvalue()==false))
    {

    vector<ST::string> path;
    vector<ST::string> path2;
    vector<ST::string> path3;
    vector<ST::string> path4;
    vector<ST::string> path5;


    for (i=0;i<b.outfiles.size();i++)
      {
      path.push_back(b.outfiles[i] + "_graphics.prg");
      path2.push_back(b.outfiles[i] + "_model_summary.tex");
      path3.push_back(b.outfiles[i] +  "_r.R");
      path4.push_back(b.outfiles[i] +  "_stata.do");
      path5.push_back(b.outfiles[i] +  "_effects.res");
      }

    b.simobj.out_effects(path5);

    b.simobj.make_graphics(header,path,path2,path3,path4);

#if defined(JAVA_OUTPUT_WINDOW)
    if(b.nographs.getvalue() == false)
    {
    for(unsigned j=0;j<b.fullcond.size();j++)
       {
       MCMC::plotstyles plst = b.fullcond[j]->get_plotstyle();
       if(plst != MCMC::noplot)
         {
         vector<ST::string> varnames = b.fullcond[j]->get_datanames();
         ST::string xvar = varnames[0];
         ST::string pathresult = b.fullcond[j]->get_pathresult();
         ST::string pathps = pathresult.substr(0, pathresult.length()-4);
         if(plst == MCMC::plotnonp)
                 {
                 b.newcommands.push_back(b.name + ".plotnonp " + ST::inttostring(j)
                 + ", title = \"Effect of " + xvar +"\" xlab = " + xvar
                 + " ylab = \" \" outfile = " + pathps + ".ps replace");
                 }

         if(plst==MCMC::drawmap)  // || plst==MCMC::drawmapgraph)
                 {
                 double u = b.fullcond[j]->get_level1();
                 double o = b.fullcond[j]->get_level2();
                 ST::string u_str = ST::doubletostring(u,0);
                 ST::string o_str = ST::doubletostring(o,0);
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", color outfile = " + pathps + "_pmean.ps replace");
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", plotvar = pcat" + u_str + " nolegend  pcat outfile = " + pathps
                 + "_pcat" + u_str + ".ps replace");
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", plotvar = pcat" + o_str + " nolegend  pcat outfile = " + pathps
                 + "_pcat" + o_str + ".ps replace");

                 }
         }
       }

    b.newcommands.push_back(b.name + ".texsummary");
    }
#endif

    }

  if (!failure && (b.family.getvalue() != "vargaussian")
     && (b.missingreg.getvalue()==false))
    {
    b.resultsyesno = true;
    }
  else
    {
    b.describetext.erase(b.describetext.begin(),b.describetext.end());
    b.describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
    b.resultsyesno = false;
    }

  if (b.family.getvalue() == "vargaussian")
    b.varianceend_fc = b.fullcond.size()-1;
  else if (b.missingreg.getvalue()==true)
    {
    if (b.begin_fc.size()==0)
      b.begin_fc.push_back(0);
    else
      b.begin_fc.push_back(b.end_fc[b.end_fc.size()-1]+1);

    b.end_fc.push_back(b.fullcond.size()-1);

    b.missingend_fc = b.fullcond.size()-1;

    }

  }


void hregressrun(bayesreg & b)
  {

  vector<ST::string> header;

  b.resultsyesno = false;
  if (b.modeonly.getvalue() == true)
    b.posteriormode = true;
  else
    b.posteriormode = false;

  b.terms = b.modreg.getterms();

  b.describetext.erase(b.describetext.begin(),b.describetext.end());
  b.describetext.push_back("LAST ESTIMATED MODEL: \n");
  b.describetext.push_back("\n");
  b.describetext.push_back(b.modreg.getModelText());
  b.describetext.push_back("\n");

  if (b.RE_est==false)
    b.clear();

  if (b.family.getvalue()=="gaussian_re")
    {
    b.RE_est=true;
    b.add_name="_gaussian_re";
    }
  else
    b.add_name="";

  b.outfiles.push_back(b.outfile.getvalue()+b.add_name);

  bool failure = false;

  if (b.family.getvalue() == "multgaussian")
    {
    failure = true;
    b.out("ERROR: family multivariate gaussian is not allowed for method hregress\n");
    }

  if (b.family.getvalue() == "multistate")
    {
    failure = true;
    b.out("ERROR: family multistate is not allowed for method hregress\n");
    }

  if (!failure)
    failure = b.create_generaloptions();

  if (!failure)
    failure = b.create_distribution();

// Speicherplatz für normalconst/nongaussianconst/nbinomialconst reservieren
/*
  unsigned nrfcfixed = b.fixedeffects.get_constvariables(b.terms).size()+1;
  unsigned blocksize_fixed = 10;
  unsigned reserved = 20;
  if ( nrfcfixed*b.nrcategories > reserved*blocksize_fixed )
    {
    nrfcfixed = ceil(nrfcfixed/double(blocksize_fixed));
    b.normalconst.reserve(nrfcfixed*b.nrcategories);
    b.nongaussianconst.reserve(nrfcfixed*b.nrcategories);
    b.nbinomialconst.reserve(nrfcfixed*b.nrcategories);
    }
*/

  unsigned i;

  if (!failure)
    {
    for (i=0;i<b.nrcategories;i++)
      {

      if (!failure)
        failure = b.create_const(i);

      if (!failure)
        failure = b.create_baseline(i);

      if (!failure)
        failure = b.create_varcoeffbaseline(i);

      if (!failure)
        failure = b.create_nonprw1rw2(i);

      if (!failure)
        failure = b.create_pspline(i);

      if (!failure)
        failure = b.create_nonpseason(i);

      if (!failure)
        failure = b.create_spatial(i);

      if (!failure)
        failure = b.create_geospline(i);

      if (!failure)
        failure = b.create_varcoeff_geospline(i);

      if (!failure)
        failure = b.create_spatialxy(i);

      if (!failure)
        failure = b.create_varcoeffpspline(i);

      if (!failure)
        failure = b.create_random(i);

      if (!failure)
        failure = b.create_hrandom(i);

      if (!failure)
        failure = b.create_randomslope(i);

      if (!failure)
        failure = b.create_mixture(i);

      if (!failure)
        failure = b.create_interactionspspline(i);

      if (!failure)
        failure = b.create_geokriging(i);

      if(!failure)
        failure = b.create_varcoeffmerror(i);

     if(!failure)
        failure = b.create_ridge(i);

     if(!failure)
        failure = b.create_lasso(i);

     if(!failure)
        failure = b.create_nigmix(i);

      if (!failure)
        failure = b.create_random_rw1rw2(i);

      if (!failure)
        failure = b.create_spatial_rw1rw2(i);

      if (!failure)
        failure = b.create_random_pspline(i);

      } // end: for (i=0;i<b.nrcategories;i++)
    } // end: if (!failure)


  if (!failure    &&
      (b.family.getvalue() != "gaussian_re")
     )
    {

    if (b.RE_est==true)
      {

      vector<unsigned> begin;
      vector<unsigned> end;
      vector<MCMCoptions*> mo;

      begin.push_back(0);


      unsigned j;
      for (j=0;j<b.REest_end_fc.size();j++)
        {
        end.push_back(b.REest_end_fc[j]);
        begin.push_back(b.REest_end_fc[j]+1);

        mo.push_back(&(b.generaloptions[j]));

      header.push_back("BAYESREG OBJECT " + b.name.to_bstr() +
                       ": Random_Effect_" + ST::inttostring(j+1) + " regression");
        }

      end.push_back(b.fullcond.size()-1);

//      b.generaloptions[0].set_nrout(b.generaloptions[0].get_iterations()+1);

      mo.push_back(&(b.generaloptions[b.generaloptions.size()-1]));

      header.push_back("BAYESREG OBJECT " + b.name.to_bstr() +
                       ": regression procedure" );


      for (j=0;j<b.normalconst_re.size();j++)
        b.normalconst_re[j].set_fcintercept(b.fcconst_intercept);


      b.simobj = MCMCsimulate(mo,b.distr,b.fullcond,begin,end);

      if (b.modeonly.getvalue())
        {
        failure = b.simobj.posteriormode(header,false);
        }
      else
        failure = b.simobj.simulate(header,b.setseed.getvalue(),!b.noposteriormode.getvalue());

      b.RE_est=false;
      }
    else
      {

      header.push_back("BAYESREG OBJECT " + b.name.to_bstr() +
                       ": regression procedure" );

      b.simobj = MCMCsimulate(&b.generaloptions[0],b.distr[0],b.fullcond);
      if (b.modeonly.getvalue())
        {
        failure = b.simobj.posteriormode(header);
        }
      else
        {
        if (b.nosamples.getvalue() == true)
          b.simobj.setflags(MCMC::nosamples);
        if ( (b.family.getvalue() == "cumprobit") ||
             (b.family.getvalue() == "multinomialprobit") ||
             (b.family.getvalue() == "binomialtlink")
           )
           {
           failure = b.simobj.simulate(header,b.setseed.getvalue(),false);
           }
        else
          {
          failure = b.simobj.simulate(header,b.setseed.getvalue(),!b.noposteriormode.getvalue());
          }

        }

      }

    }


  if (!failure && (b.family.getvalue() != "gaussian_re"))
    {

    vector<ST::string> path;
    vector<ST::string> path2;
    vector<ST::string> path3;
    vector<ST::string> path4;
    vector<ST::string> path5;


    for (i=0;i<b.outfiles.size();i++)
      {
      path.push_back(b.outfiles[i] + "_graphics.prg");
      path2.push_back(b.outfiles[i] + "_model_summary.tex");
      path3.push_back(b.outfiles[i] +  "_r.R");
      path4.push_back(b.outfiles[i] +  "_stata.do");
      path5.push_back(b.outfiles[i] +  "_effects.res");
      }

    b.simobj.out_effects(path5);

    b.simobj.make_graphics(header,path,path2,path3,path4);

#if defined(JAVA_OUTPUT_WINDOW)
    if(b.nographs.getvalue() == false)
    {
    for(unsigned j=0;j<b.fullcond.size();j++)
       {
       MCMC::plotstyles plst = b.fullcond[j]->get_plotstyle();
       if(plst != MCMC::noplot)
         {
         vector<ST::string> varnames = b.fullcond[j]->get_datanames();
         ST::string xvar = varnames[0];
         ST::string pathresult = b.fullcond[j]->get_pathresult();
         ST::string pathps = pathresult.substr(0, pathresult.length()-4);
         if(plst == MCMC::plotnonp)
                 {
                 b.newcommands.push_back(b.name + ".plotnonp " + ST::inttostring(j)
                 + ", title = \"Effect of " + xvar +"\" xlab = " + xvar
                 + " ylab = \" \" outfile = " + pathps + ".ps replace");
                 }

         if(plst==MCMC::drawmap)  // || plst==MCMC::drawmapgraph)
                 {
                 double u = b.fullcond[j]->get_level1();
                 double o = b.fullcond[j]->get_level2();
                 ST::string u_str = ST::doubletostring(u,0);
                 ST::string o_str = ST::doubletostring(o,0);
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", color outfile = " + pathps + "_pmean.ps replace");
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", plotvar = pcat" + u_str + " nolegend  pcat outfile = " + pathps
                 + "_pcat" + u_str + ".ps replace");
                 b.newcommands.push_back(b.name + ".drawmap " + ST::inttostring(j)
                 + ", plotvar = pcat" + o_str + " nolegend  pcat outfile = " + pathps
                 + "_pcat" + o_str + ".ps replace");

                 }
         }
       }

    b.newcommands.push_back(b.name + ".texsummary");
    }
#endif

    }

  if (!failure && (b.family.getvalue() != "gaussian_re"))
    {
    b.resultsyesno = true;
    }
  else
    {
    b.describetext.erase(b.describetext.begin(),b.describetext.end());
    b.describetext.push_back("CURRENT REGRESSION RESULTS: none\n");
    b.resultsyesno = false;
    }

  if (b.family.getvalue() == "gaussian_re")
    b.REest_end_fc.push_back(b.fullcond.size()-1);

  }







