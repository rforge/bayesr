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
#include"bayesreg2.h"

// Vorschlag:
//#include<typeinfo.h>

#include<stddef.h>

bool bayesreg::create_geokriging(const unsigned & collinpred)
  {

  ST::string pathnonp;
  ST::string pathres;

  long h;
  double nu,maxdist,p,q,lambda,startlambda;
  double a1,b1;
  unsigned nrknots, maxsteps;
  bool full;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial_geokriging.checkvector(terms,i) == true )
      {

      MCMC::fieldtype type;
      if (terms[i].options[0] == "geokriging")
        type = MCMC::kriging;

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[2]).strtodouble(nu);
      if(nu!=0.5 && nu!=1.5&& nu!=2.5 && nu!=3.5)
        {
        outerror("ERROR: Invalid value for nu\n");
        return true;
        }
      f = (terms[i].options[3]).strtodouble(maxdist);
      if(maxdist<=0) // wähle maxdist so, dass Korrelation für Punkte mitmaximalem Abstand = 0.0001
        {
        if(nu==0.5)
          {
          maxdist=9.21034037;//4.605170186;
          }
        else if(nu==1.5)
          {
          maxdist=11.75637122;//6.638352068;
          }
        else if(nu==2.5)
          {
          maxdist=13.53592464;//8.022007057;
          }
        else if(nu==3.5)
          {
          maxdist=15.01510426;//9.158140446;
          }
        }

      if(terms[i].options[4] == "true")
        {
        full=true;
        }
      else
        {
        full=false;
        }

      f = (terms[i].options[6]).strtodouble(p);
      f = (terms[i].options[7]).strtodouble(q);
      f = (terms[i].options[8]).strtolong(h);
      maxsteps = unsigned(h);

      f = (terms[i].options[9]).strtodouble(lambda);
      f = (terms[i].options[10]).strtodouble(startlambda);

      f = (terms[i].options[12]).strtodouble(a1);
      f = (terms[i].options[13]).strtodouble(b1);

      if (f==1)
        return true;

      mapobject * mapp;                           // pointer to mapobject
      int objpos = findstatobject(*statobj,terms[i].options[11],"map");
      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[11] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[11] + " is not a map object\n");
        return true;
        }
      MAP::map m = mapp->getmap();
      if(!m.centroids_existing())
        {
        outerror("ERROR: map object does not contain centroids\n");
        return true;
        }

      datamatrix knotdata;
      if(terms[i].options[5]!="" && !full)
        {
        dataobject * datap;                           // pointer to datasetobject
        int objpos = findstatobject(*statobj,terms[i].options[5],"dataset");
        if (objpos >= 0)
          {
          statobject * s = statobj->at(objpos);
          datap = dynamic_cast<dataobject*>(s);
          if (datap->obs()==0 || datap->getVarnames().size()==0)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " does not contain any data\n");
            return true;
            }
          else if (datap->getVarnames().size()>2)
            {
            outerror("ERROR: dataset object " + terms[i].options[5] + " contains more than two variables\n");
            return true;
            }
          }
        else
          {
          outerror("ERROR: dataset object " + terms[i].options[5] + " is not existing\n");
          return true;
          }
        list<ST::string> knotnames = datap->getVarnames();
        ST::string expr = "";
        datap->makematrix(knotnames,knotdata,expr);
        }
      else
        {
        knotdata = datamatrix(1,1,0);
        }

      ST::string title;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geokriging.raw","_geokriging.res","_geokriging");

      if ( check_gaussian(collinpred) )
        {

        fckriging.push_back(
        FULLCOND_kriging2(&generaloptions[generaloptions.size()-1],
                distr[distr.size()-1],fcconst_intercept,
                D.getCol(j),m,terms[i].options[11],
                knotdata,
                nrknots,nu,maxdist,p,q,maxsteps,full,type,
                title,
                pathnonp,
                pathres,
                lambda,
                startlambda,
                collinpred
                ));

        }
      else if( check_iwls(true,collinpred) )
        {

        ST::string proposal = terms[i].options[14];

        bool iwlsmode;
        if(proposal == "iwlsmode")
          iwlsmode = true;
        else
          iwlsmode = false;

        unsigned updateW;
        f = (terms[i].options[15]).strtolong(h);
        updateW = unsigned(h);

        bool updatetau;
        if(terms[i].options[16] == "false" || constlambda.getvalue() == true)
          updatetau = false;
        else
          updatetau = true;

        double fstart;
          f = (terms[i].options[17]).strtodouble(fstart);

        fckriging.push_back(
        FULLCOND_kriging2(&generaloptions[generaloptions.size()-1],
                distr[distr.size()-1],fcconst_intercept,
                D.getCol(j),m,terms[i].options[11],
                knotdata,
                nrknots,nu,maxdist,p,q,maxsteps,full,type,
                title,
                pathnonp,
                pathres,
                lambda,
                startlambda,
                iwlsmode,
                updateW,
                updatetau,
                fstart,
                collinpred
                ));

        }

      if (constlambda.getvalue() == true)
        fckriging[fckriging.size()-1].set_lambdaconst(lambda);

      vector<ST::string> na;
      na.push_back(terms[i].varnames[0]);

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                     "_geokriging_var.raw","_geokriging_var.res","_geokriging_variance");

      fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                     &fckriging[fckriging.size()-1],
                                     distr[distr.size()-1],
                                     a1,
                                     b1,
                                     title,pathnonp,pathres,
                                     false,collinpred)
                         );

      fckriging[fckriging.size()-1].init_names(na);
      fckriging[fckriging.size()-1].set_fcnumber(fullcond.size());
      fullcond.push_back(&fckriging[fckriging.size()-1]);

     if (constlambda.getvalue() == false)
        {
        if(terms[i].options[18]=="true")
          fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
        }

      }
    }

  return false;
  }


bool bayesreg::create_spatial(const unsigned & collinpred)
  {

  ST::string pathnonpv;
  ST::string pathresv;

  long h;
  double hd;
  unsigned min,max;
  int f;
  double lambda,a1,b1,alpha;
  unsigned i;
  int j1=0,j2=0;
  bool iwls;
  bool varcoeff=false;
  unsigned updateW;
  bool updatetau,Laplace;
  double ftune;
  ST::string proposal;
  vector<ST::string> na;
  bool center;

  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatial.checkvector(terms,i) == true)
      {

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);  // interacting var
      if (terms[i].type == "varcoeffspatial" || terms[i].type == "tvarcoeffspatial")
        {
        varcoeff=true;
        j2 = terms[i].varnames[1].isinlist(modelvarnamesv);  // effect modifier
        }

      mapobject * mapp;                           // pointer to mapobject

      int objpos = findstatobject(*statobj,terms[i].options[1],"map");

      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          {
          if ((terms[i].options[1] == "") || (terms[i].options[1] == " "))
            outerror("ERROR: map object must be specified to estimate a spatial effect\n");
          else
            outerror("ERROR: map object " + terms[i].options[1] + " is not existing\n");
          }
        else
          outerror("ERROR: " + terms[i].options[1] + " is not a map object\n");
        return true;
        }

      MAP::map m = mapp->getmap();
      bool isconnected = m.isconnected();
      if (isconnected==false)
        {
        outerror("ERROR: map is disconnected, spatial effect cannot be estimated\n");
        return true;
        }


      f = (terms[i].options[2]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[3]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[4]).strtodouble(hd);
      lambda = hd;

      f = (terms[i].options[5]).strtodouble(a1);

      f = (terms[i].options[6]).strtodouble(b1);

      proposal = terms[i].options[7];

      f = (terms[i].options[8]).strtolong(h);
      updateW = unsigned(h);

      if (terms[i].options[9] == "true")
        updatetau=true;
      else
        updatetau=false;

      f = (terms[i].options[10]).strtodouble(ftune);

      if (terms[i].options[16] == "true")
        Laplace=true;
      else
        Laplace=false;

      f = (terms[i].options[18]).strtodouble(alpha);

      if (f==1)
        return true;

      if (terms[i].options[20] == "true")
        center=true;
      else
        center=false;


      ST::string titlev;

      if (varcoeff == true)
        {
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                   terms[i].varnames[0],"_spatial.raw","_spatial.res","_spatial");

        make_paths(collinpred,pathnonpv,pathresv,titlev,terms[i].varnames[1],
                   terms[i].varnames[0],"_spatial_var.raw","_spatial_var.res",
                   "_spatial_variance");
        }
      else
        {
        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_spatial.raw","_spatial.res","_spatial");

        make_paths(collinpred,pathnonpv,pathresv,titlev,terms[i].varnames[0],"",
                 "_spatial_var.raw","_spatial_var.res","_spatial_variance");
        }

      if (proposal != "cp")
        iwls=true;
      else
        iwls=false;

      if ( (check_gaussian(collinpred)) || (check_iwls(iwls,collinpred)) )
        {

        if (varcoeff == true)
          {
          fcnonpgaussian.push_back(
          FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
          distr[distr.size()-1],fcconst_intercept,m,terms[i].options[1],
          D.getCol(j2),D.getCol(j1),title,pathnonp,pathres,collinpred,lambda,
          center));

          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);
          }
        else
          {
          fcnonpgaussian.push_back(
          FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
          distr[distr.size()-1],D.getCol(j1),fcconst_intercept,m,
          terms[i].options[1],title,pathnonp,pathres,collinpred,lambda));

          fcnonpgaussian[fcnonpgaussian.size()-1].init_name(terms[i].varnames[0]);
          }

        if(Laplace)
          fcnonpgaussian[fcnonpgaussian.size()-1].set_Laplace();

        if (fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size() > 0)
          {
          unsigned i;
          for(i=0;i<fcnonpgaussian[fcnonpgaussian.size()-1].get_errors().size();i++)
            errormessages.push_back(
            fcnonpgaussian[fcnonpgaussian.size()-1].get_errors()[i]);
          return true;
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

        if (terms[i].options[17] == "true")
          fcnonpgaussian[fcnonpgaussian.size()-1].set_stationary(alpha);

        fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);

        if (constlambda.getvalue() == false)
          {

          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
          &fcnonpgaussian[fcnonpgaussian.size()-1],distr[distr.size()-1],a1,b1,
          titlev,pathnonpv,pathresv,false,collinpred));

          if(terms[i].options[14]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          if ( (check_nongaussian(collinpred)) && (proposal == "iwls")
            && (updatetau==true) )
            fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();

          if ( (check_nongaussian(collinpred)) && (proposal == "iwlsmode")
            && (updatetau==true) )
            fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();

          bool alphafix = false;
          if (terms[i].options[19] == "true")
            alphafix = true;
          if (terms[i].options[17] == "true")
            fcvarnonp[fcvarnonp.size()-1].set_stationary(alpha,alphafix);

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

          if(Laplace)
            fcvarnonp[fcvarnonp.size()-1].set_Laplace();

          }

        if ((terms[i].options[0] == "tspatial") ||
            (terms[i].options[0] == "tvarcoeffspatial")
           )
          {

          if(varcoeff==true)
            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],terms[i].varnames[0],
                   "_spatial_tvar.raw","_spatial_tvar.res","_spatial_tvariance");
          else
            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_spatial_tvar.raw","_spatial_tvar.res","_spatial_tvariance");

          unsigned v = nu.getvalue();
          unsigned nrrows;
          f = (terms[i].options[15]).strtolong(h);
          nrrows = unsigned(h);

          fctvariance2dim.push_back(FULLCOND_tvariance2dim(&generaloptions[generaloptions.size()-1],
          &fcnonpgaussian[fcnonpgaussian.size()-1],v,title,
          pathnonp,pathres,nrrows,false));

          fctvariance2dim[fctvariance2dim.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fctvariance2dim[fctvariance2dim.size()-1]);

          if(Laplace)
            fctvariance2dim[fctvariance2dim.size()-1].set_Laplace();

          }

        }
      else
        {

        if (varcoeff == true)
          {
          Pmatrices.push_back(PenaltyMatrix(D.getCol(j2),terms[i].varnames[1],
          m,min,max));
          }
        else
          Pmatrices.push_back(PenaltyMatrix(D.getCol(j1),terms[i].varnames[0],
          m,min,max));

        if (Pmatrices[Pmatrices.size()-1].get_errormessages().size() > 0)
          {
          outerror(Pmatrices[Pmatrices.size()-1].get_errormessages());
          return true;
          }

        fcnonp.push_back( FULLCOND_nonp(
                                       &generaloptions[generaloptions.size()-1],
                                        distr[distr.size()-1],
                                        &Pmatrices[Pmatrices.size()-1],
                                        fcconst_intercept,
                                        lambda,
                                        pathnonp,
                                        pathres,title,terms[i].options[1],
                                        collinpred
                                        )
                           );

        fcnonp[fcnonp.size()-1].init_name(terms[i].varnames[0]);

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


      }   // end: if ( nonpspatial.checkvector(terms,i) == true )

    } //  end:  for(i=0;i<terms.size();i++)


  return false;
  }


bool bayesreg::create_spatialxy(const unsigned & collinpred)
  {

  ST::string pathmap;
  ST::string mapname;

//  long h;
//  unsigned min,max;
  double a1,b1;
  double lambda;
  double maxdist;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpspatialxy.checkvector(terms,i) == true )
      {

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);


//      f = (terms[i].options[1]).strtolong(h);
//      min = unsigned(h);

//      f = (terms[i].options[2]).strtolong(h);
//      max = unsigned(h);

      f = (terms[i].options[3]).strtodouble(lambda);

      f = (terms[i].options[4]).strtodouble(a1);

      f = (terms[i].options[5]).strtodouble(b1);

      f = (terms[i].options[6]).strtodouble(maxdist);

      if (f==1)
        return true;

      pathmap = outfile.getvalue()+ add_name + "_" + terms[i].varnames[0] + "_" +
                terms[i].varnames[1] + "_dist" + ST::doubletostring(maxdist)
                + ".bnd";

      mapname = terms[i].varnames[0] + "_" + terms[i].varnames[1] + "_dist"
                + ST::doubletostring(maxdist) + add_name;

      ST::string help = terms[i].varnames[0]+"_"+terms[i].varnames[1];

      make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_spatial.raw","_spatial.res","_spatial");

      if (check_gaussian(collinpred))
        {

        fcnonpgaussian.push_back( FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
                                        distr[distr.size()-1],
                                        D.getCol(j1),D.getCol(j2),
                                        fcconst_intercept,
                                        lambda,maxdist,mapname,title,
                                        pathnonp, pathres,pathmap,collinpred
                                        )
                           );

        fcnonpgaussian[fcnonpgaussian.size()-1].init_name("regionnr");

        help = terms[i].varnames[0]+"_"+terms[i].varnames[1];

        make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_spatial_var.raw","_spatial_var.res","_spatial_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fcnonpgaussian[fcnonpgaussian.size()-1],
                                       distr[distr.size()-1],a1,b1,title,
                                       pathnonp,pathres,false,collinpred)
                           );

        if (constlambda.getvalue() == true)
          fcvarnonp[fcvarnonp.size()-1].set_constlambda();



        fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);
        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

        }
      else
        {


        }


      }   // end: if ( nonpspatial.checkvector(terms,i) == true )

    } //  end:  for(i=0;i<terms.size();i++)


  return false;
  }


bool bayesreg::create_randomslope(const unsigned & collinpred)
  {

  ST::string pathnonp2;
  ST::string pathres2;
  ST::string pathresfixed;
  ST::string title2;

  bool iwlsmode=false;
  bool updatetau;

  unsigned i;
  int j1,j2;
  bool inclf;
  double a1,b1,lambda;
  int f;
  for(i=0;i<terms.size();i++)
    {
    if ( randomeffslope.checkvector(terms,i) == true )
      {

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);

      if (terms[i].options[1] == "true")
        inclf = false;
      else
        inclf = true;

      f = (terms[i].options[2]).strtodouble(lambda);

      f = (terms[i].options[3]).strtodouble(a1);

      f = (terms[i].options[4]).strtodouble(b1);

      if (terms[i].options[5] == "iwlsmode")
        iwlsmode = true;

      if (terms[i].options[6] == "true")
        updatetau = true;
      else
        updatetau = false;

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
                 terms[i].varnames[0],
                 "_random.raw","_random.res","_random");

      pathresfixed = pathres.substr(0,pathres.length()-4) + "_fixed.res";

      make_paths(collinpred,pathnonp2,pathres2,title2,terms[i].varnames[1],
                 terms[i].varnames[0],"_random_var.raw",
                 "_random_var.res","_random_variance");






      if (check_gaussian(collinpred))
        {

        fcrandomgaussian.push_back(
        FULLCOND_random_gaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j1),
                                                        D.getCol(j2),
                                                        title,
                                                        pathnonp,
                                                        pathres,pathresfixed,
                                                        lambda,
                                                        inclf,collinpred
                                                        )
                          );
// Begin: DSB
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
        if(mscheck.getvalue()==true)
          fcrandomgaussian[fcrandomgaussian.size()-1].set_mscheck(true);
#endif
// End: DSB
        vector<ST::string> na;
        na.push_back(terms[i].varnames[1]);
        na.push_back(terms[i].varnames[0]);
        fcrandomgaussian[fcrandomgaussian.size()-1].init_names(na);

        fcrandomgaussian[fcrandomgaussian.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true || terms[i].options[8]=="true")
          {
          fcrandomgaussian[fcrandomgaussian.size()-1].set_lambdaconst(lambda);
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);

          }
        else
          {
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);
          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandomgaussian[fcrandomgaussian.size()-1],
                            distr[distr.size()-1],a1,b1,title2,pathnonp2,
                            pathres2,false,collinpred));
          if(terms[i].options[7]=="true")
                fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }


        }

      else
        {

        fcrandom.push_back(
        FULLCOND_random_nongaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j1),
                                                        D.getCol(j2),
                                                        title,
                                                        pathnonp,
                                                        pathres,pathresfixed,
                                                        lambda,iwlsmode,inclf,
                                                        collinpred
                                                        )
                          );
// Begin: DSB
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
        if(mscheck.getvalue()==true)
          fcrandom[fcrandom.size()-1].set_mscheck(true);
#endif
// End: DSB

        vector<ST::string> na;
        na.push_back(terms[i].varnames[1]);
        na.push_back(terms[i].varnames[0]);
        fcrandom[fcrandom.size()-1].init_names(na);

        fcrandom[fcrandom.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true || terms[i].options[8]=="true")
          {
          fcrandom[fcrandom.size()-1].set_lambdaconst(lambda);
          fullcond.push_back(&fcrandom[fcrandom.size()-1]);
          }
        else
          {

          fullcond.push_back(&fcrandom[fcrandom.size()-1]);

          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandom[fcrandom.size()-1],
                            distr[distr.size()-1],a1,b1,
                            title2,pathnonp2,pathres2,false,collinpred)
                            );

          if(updatetau)
            fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();
          if(terms[i].options[7]=="true")
                fcvarnonp[fcvarnonp.size()-1].set_uniformprior();


          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

          }

        }

      }

    }

  return false;

  }


bool bayesreg::create_random(const unsigned & collinpred)
  {

  ST::string pathnonp2;
  ST::string pathres2;
  ST::string title2;
  double a1,b1,lambda;
  int f;
//  long h;
  bool iwlsmode=false;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( randomeff.checkvector(terms,i) == true )
      {

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtodouble(lambda);

      f = (terms[i].options[2]).strtodouble(a1);

      f = (terms[i].options[3]).strtodouble(b1);

      if (terms[i].options[4]=="iwlsmode")
        iwlsmode=true;


      if (f==1)
        return true;


      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_random.raw","_random.res","_random");

      make_paths(collinpred,pathnonp2,pathres2,title2,terms[i].varnames[0],"",
                   "_random_var.raw","_random_var.res","_random_variance");


      if (check_gaussian(collinpred))
        {

        FULLCOND_nonp_gaussian * structuredp = NULL;
        unsigned structured=0;

        unsigned j1;
        for (j1=0;j1<fcnonpgaussian.size();j1++)
          {
          if  ( ((fcnonpgaussian[j1].get_datanames()).size() == 1) &&
                (fcnonpgaussian[j1].get_datanames()[0]==terms[i].varnames[0]) &&
                (fcnonpgaussian[j1].get_col() == collinpred) &&
                (fcnonpgaussian[j1].get_type() == MCMC::mrf)
              )
            {
            structuredp = &fcnonpgaussian[j1];
            structured ++;
            }
          }

        fcrandomgaussian.push_back(
        FULLCOND_random_gaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j),
                                                        title,
                                                        pathnonp,pathres,
                                                        lambda,collinpred
                                                        )
                          );
// Begin: DSB
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
        if(mscheck.getvalue()==true)
          fcrandomgaussian[fcrandomgaussian.size()-1].set_mscheck(true);
#endif
// End: DSB
        if (structured==1)
          {
#if defined(__BUILDING_LINUX)
          ST::string pathnonpt = defaultpath + "/temp/" + name + add_name +
                 terms[i].varnames[0] +
                 "_spatialtotal.raw";
#else
          ST::string pathnonpt = defaultpath + "\\temp\\" + name + add_name +
                 terms[i].varnames[0] +
                 "_spatialtotal.raw";
#endif
          ST::string pathrest = outfile.getvalue() + add_name + "_" + terms[i].varnames[0] +
                                "_spatialtotal.res";

          fcrandomgaussian[fcrandomgaussian.size()-1].init_spatialtotal(
          structuredp,pathnonpt,pathrest);
          }
        else if (structured==0)
          {
          }
        else
          {
          outerror("ERROR: more than one spatial effect specified for variable "
                   + terms[i].varnames[0] + "\n");
          return true;
          }

        fcrandomgaussian[fcrandomgaussian.size()-1].init_name(terms[i].varnames[0]);
        fcrandomgaussian[fcrandomgaussian.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true || terms[i].options[7]=="true")
          {
          fcrandomgaussian[fcrandomgaussian.size()-1].set_lambdaconst(lambda);
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);
          }
        else
          {
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);
          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandomgaussian[fcrandomgaussian.size()-1],
                            distr[distr.size()-1],a1,b1,title2,
                            pathnonp2,pathres2,false,collinpred)
                            );
          if(terms[i].options[6]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

          }

        } // end: gaussian, probit, etc. ...

      else
        {

        FULLCOND_nonp * structuredp1 = NULL;
        unsigned structured1=0;

        unsigned j1;
        for (j1=0;j1<fcnonp.size();j1++)
          {
          if  ( ((fcnonp[j1].get_datanames()).size() == 1) &&
                (fcnonp[j1].get_datanames()[0]==terms[i].varnames[0]) &&
                (fcnonp[j1].get_col() == collinpred) //&&
                // (fcnonp[j1].get_type() == MCMC::mrf)
              )
            {
            structuredp1 = &fcnonp[j1];
            structured1 ++;
            }
          }


        FULLCOND_nonp_gaussian * structuredp2 = NULL;
        unsigned structured2=0;

        for (j1=0;j1<fcnonpgaussian.size();j1++)
          {
          if  ( ((fcnonpgaussian[j1].get_datanames()).size() == 1) &&
                (fcnonpgaussian[j1].get_datanames()[0]==terms[i].varnames[0]) &&
                (fcnonpgaussian[j1].get_col() == collinpred) &&
                (fcnonpgaussian[j1].get_type() == MCMC::mrf)
              )
            {
            structuredp2 = &fcnonpgaussian[j1];
            structured2 ++;
            }
          }



        fcrandom.push_back( FULLCOND_random_nongaussian(
        &generaloptions[generaloptions.size()-1], distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j),
                                                        title,
                                                        pathnonp,
                                                        pathres,
                                                        lambda,
                                                        iwlsmode,
                                                        collinpred
                                                        )
                          );
// Begin: DSB
#if !defined (__BUILDING_THE_DLL) & !defined(__BUILDING_GNU)
        if(mscheck.getvalue()==true)
          fcrandom[fcrandom.size()-1].set_mscheck(true);
#endif
// End: DSB

       if  ( ( (structured1==1) && (structured2==0) ) ||
             ( (structured1==0) && (structured2==1) )
           )
          {
#if defined(__BUILDING_LINUX)
          ST::string pathnonpt = defaultpath + "/temp/" + name + add_name +
                 terms[i].varnames[0] +
                 "_spatialtotal.raw";
#else
          ST::string pathnonpt = defaultpath + "\\temp\\" + name + add_name +
                 terms[i].varnames[0] +
                 "_spatialtotal.raw";
#endif
          ST::string pathrest = outfile.getvalue()
                                + add_name + "_" + terms[i].varnames[0] +
                                "_spatialtotal.res";

          if (structured1==1)
            fcrandom[fcrandom.size()-1].init_spatialtotal(
            structuredp1,pathnonpt,pathrest);
          else
            fcrandom[fcrandom.size()-1].init_spatialtotal(
            structuredp2,pathnonpt,pathrest);
          }
        else if (structured1==0 && structured2==0)
          {
          }
        else
          {
          outerror("ERROR: more than one spatial effect specified for variable "
                   + terms[i].varnames[0] + "\n");
          return true;
          }


        fcrandom[fcrandom.size()-1].init_name(terms[i].varnames[0]);
        fcrandom[fcrandom.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true || terms[i].options[7]=="true")
          {
          fcrandom[fcrandom.size()-1].set_lambdaconst(lambda);
          fullcond.push_back(&fcrandom[fcrandom.size()-1]);
          }
        else
          {

          fullcond.push_back(&fcrandom[fcrandom.size()-1]);

          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandom[fcrandom.size()-1],
                            distr[distr.size()-1],a1,b1,title2,
                            pathnonp2,pathres2,false,collinpred)
                            );
          if(terms[i].options[6]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

          }

        }

      }

    }

  return false;

  }


bool bayesreg::create_hrandom(const unsigned & collinpred)
  {

  ST::string pathnonp2;
  ST::string pathres2;
  ST::string title2;
  double a1,b1,lambda;
  int f;
//  long h;
  bool iwlsmode=false;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( hrandomeff.checkvector(terms,i) == true )
      {

      j = terms[i].varnames[0].isinlist(modelvarnamesv);

      f = (terms[i].options[1]).strtodouble(lambda);

      f = (terms[i].options[2]).strtodouble(a1);

      f = (terms[i].options[3]).strtodouble(b1);

      if (terms[i].options[4]=="iwlsmode")
        iwlsmode=true;


      if (f==1)
        return true;


      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_hrandom.raw","_hrandom.res","_hrandom");

      make_paths(collinpred,pathnonp2,pathres2,title2,terms[i].varnames[0],"",
                   "_hrandom_var.raw","_hrandom_var.res","_hrandom_variance");


      if (check_gaussian(collinpred))
        {

        fchrandom.push_back(
        FULLCOND_hrandom(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j),
                                                        title,
                                                        pathnonp,pathres,
                                                        lambda,collinpred
                                                        )
                          );

        fchrandom[fchrandom.size()-1].init_name(terms[i].varnames[0]);
        fchrandom[fchrandom.size()-1].set_fcnumber(fullcond.size());

        fchrandom[fchrandom.size()-1].set_REdistr(&distr_gaussian_re[0]);

        if (constlambda.getvalue() == true || terms[i].options[7]=="true")
          {
          fchrandom[fchrandom.size()-1].set_lambdaconst(lambda);
          fullcond.push_back(&fchrandom[fchrandom.size()-1]);
          }
        else
          {
          fullcond.push_back(&fchrandom[fchrandom.size()-1]);
          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fchrandom[fchrandom.size()-1],
                            distr[distr.size()-1],a1,b1,title2,
                            pathnonp2,pathres2,false,collinpred)
                            );
          if(terms[i].options[6]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

          }

        } // end: gaussian, probit, etc. ...

      else
        {

        /*

        fcrandom.push_back( FULLCOND_random_nongaussian(
        &generaloptions[generaloptions.size()-1], distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j),
                                                        title,
                                                        pathnonp,
                                                        pathres,
                                                        lambda,
                                                        iwlsmode,
                                                        collinpred
                                                        )
                          );



        fcrandom[fcrandom.size()-1].init_name(terms[i].varnames[0]);
        fcrandom[fcrandom.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true || terms[i].options[7]=="true")
          {
          fcrandom[fcrandom.size()-1].set_lambdaconst(lambda);
          fullcond.push_back(&fcrandom[fcrandom.size()-1]);
          }
        else
          {

          fullcond.push_back(&fcrandom[fcrandom.size()-1]);

          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandom[fcrandom.size()-1],
                            distr[distr.size()-1],a1,b1,title2,
                            pathnonp2,pathres2,false,collinpred)
                            );
          if(terms[i].options[6]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

          }

         */

        }

      }

    }

  return false;

  }



bool bayesreg::create_mixture(const unsigned & collinpred)
  {
  int f;
  unsigned nrcomp,aclag;
  double wprior,mpriorm,mpriorv,vpriora,vpriorb;
  long h;
  bool nosamples=false,vpriorbunif=false,vpriorbgamma=false;
  ST::string ordertype="n";

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( mixtureeff.checkvector(terms,i) == true )
      {

      j = terms[i].varnames[0].isinlist(modelvarnamesv);
      f = (terms[i].options[1]).strtolong(h);
      nrcomp = unsigned(h);
      f = (terms[i].options[2]).strtodouble(wprior);
      f = (terms[i].options[3]).strtodouble(mpriorm);
      f = (terms[i].options[4]).strtodouble(mpriorv);
      f = (terms[i].options[5]).strtodouble(vpriora);
      f = (terms[i].options[6]).strtodouble(vpriorb);
      if(terms[i].options[7] == "true")
        nosamples = true;
      f = (terms[i].options[8]).strtolong(h);
      aclag = unsigned(h);
      if(terms[i].options[9] == "w")
        ordertype = "w";
      if(terms[i].options[10] == "true")
        vpriorbunif = true;
      if(terms[i].options[11] == "true")
        vpriorbgamma = true;

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_mixture.raw","_mixture.res","_mixture");

      if (check_gaussian(collinpred))
        {


           fcmixture.push_back(
           FULLCOND_mixture(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        D.getCol(j),
                                                        title,
                                                        pathnonp,pathres,
                                                        nrcomp,wprior,
                                                        mpriorm,mpriorv,
                                                        vpriora,vpriorb,
                                                        nosamples,aclag,
                                                        ordertype,
                                                        vpriorbunif,vpriorbgamma,
                                                        collinpred
                                                        )
                          );

          fcmixture[fcmixture.size()-1].init_name(terms[i].varnames[0]);
          fcmixture[fcmixture.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcmixture[fcmixture.size()-1]);
        } // end: gaussian,???

      else
        {
        outerror("ERROR: only family=gaussian allowed for variable "
                   + terms[i].varnames[0] + "\n");
        }

      }

    }

  return false;

  }


bool bayesreg::create_interactionspspline(const unsigned & collinpred)
  {

  ST::string proposal;
  bool varcoeff;
  int j1=0,j2=0,j3=0;

  long h;
  double lambda;
//  bool reduced, singleblock;
  bool  singleblock;
  unsigned min,max,nrknots,degree,blocksize;
  bool center;
  double a1,b1;
  int gridsize;
  int f;

  unsigned i;


  for(i=0;i<terms.size();i++)
    {
    if ( nonpinteractpspline.checkvector(terms,i) == true )
      {

      MCMC::fieldtype type;

      if (terms[i].varnames.size()==2)
        varcoeff=false;
      else
        varcoeff=true;

      if ((terms[i].options[0] == "pspline2dimrw1")   ||
          (terms[i].options[0] == "tpspline2dimrw1") ||
          (terms[i].options[0] == "varpspline2dimrw1")
         )
        type = MCMC::mrflinear;
      else if ( (terms[i].options[0] == "pspline2dimrw2") ||
                (terms[i].options[0] == "varpspline2dimrw2")
              )
        type = MCMC::mrfquadratic8;
      else if ((terms[i].options[0] == "pspline2dimband")   ||
              (terms[i].options[0] == "varpspline2dimband") ||
          (terms[i].options[0] == "tpspline2dimband")
         )
        type = MCMC::mrflinearband;
      else if ( (terms[i].options[0] == "psplinekrrw1") ||
                (terms[i].options[0] == "varpsplinekrrw1")
              )
        type = MCMC::mrfkr1;
      else if ( (terms[i].options[0] == "psplinekrrw2") ||
                (terms[i].options[0] == "varpsplinekrrw2")
              )
        type = MCMC::mrfkr2;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);


      if (varcoeff==true)
        j3 = terms[i].varnames[2].isinlist(modelvarnamesv);


      f = (terms[i].options[1]).strtolong(h);
      min = unsigned(h);

      f = (terms[i].options[2]).strtolong(h);
      max = unsigned(h);

      f = (terms[i].options[3]).strtolong(h);
      degree = unsigned(h);

      f = (terms[i].options[4]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[5]).strtodouble(lambda);

//      if (terms[i].options[6] == "false")
//        reduced = false;
//      else
//        reduced = true;

      f = (terms[i].options[7]).strtodouble(a1);

      f = (terms[i].options[8]).strtodouble(b1);

      if (terms[i].options[9] == "false")
        singleblock = false;
      else
        singleblock = true;

      f = (terms[i].options[10]).strtolong(h);
      gridsize = unsigned(h);

      proposal = terms[i].options[11];

      if (terms[i].options[17] == "false")
        center = false;
      else
        center = true;


      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      if (f==1)
        return true;

      ST::string help;

      if (varcoeff==false)
        {
        help  = terms[i].varnames[0] + "_" + terms[i].varnames[1];

        make_paths(collinpred,pathnonp,pathres,title,help,"",
                   "_pspline.raw","_pspline.res","_pspline");
        }
      else
        {
        help  = terms[i].varnames[0] + "_" + terms[i].varnames[1]+
                           "_" + terms[i].varnames[2];

        make_paths(collinpred,pathnonp,pathres,title,help,"",
                   "_pspline.raw","_pspline.res","_pspline");
        }


      if (check_gaussian(collinpred))
        {

        FULLCOND_pspline_gaussian * mainp1=NULL;
        FULLCOND_pspline_gaussian * mainp2=NULL;
        unsigned main1=0;
        unsigned main2=0;

        unsigned j;
        for (j=0;j<fcpsplinegaussian.size();j++)
          {
          if  ( ((fcpsplinegaussian[j].get_datanames()).size() == 1) &&
               (fcpsplinegaussian[j].get_datanames()[0] == terms[i].varnames[0]) &&
                fcpsplinegaussian[j].get_col() == collinpred
              )
            {
            mainp1 = &fcpsplinegaussian[j];
            if(mainp1->get_gridsize() != gridsize)
              {
              outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
              return true;
              }
            if(mainp1->get_nrknots() != nrknots)
              {
              outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
              return true;
              }
            if(mainp1->get_degree() != degree)
              {
              outerror("ERROR: degree for interaction term and main effects must be the same\n");
              return true;
              }
            main1 ++;
            }


          if  ( ((fcpsplinegaussian[j].get_datanames()).size() == 1) &&
              (fcpsplinegaussian[j].get_datanames()[0] == terms[i].varnames[1]) &&
              fcpsplinegaussian[j].get_col() == collinpred
              )
            {
            mainp2 = &fcpsplinegaussian[j];
            if(mainp2->get_gridsize() != gridsize)
              {
              outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
              return true;
              }
            if(mainp2->get_nrknots() != nrknots)
              {
              outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
              return true;
              }
            if(mainp2->get_degree() != degree)
              {
              outerror("ERROR: degree for interaction term and main effects must be the same\n");
              return true;
              }
            main2 ++;
            }

          }

        if (varcoeff==false)
          {
          fcpsplinesurfgaussian.push_back(
          FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],
                                        fcconst_intercept,
                                        D.getCol(j1),
                                        D.getCol(j2),
                                        title,
                                        nrknots,degree,po,
                                        lambda,
                                        gridsize,
                                        type,
                                        pathnonp,
                                        pathres,
                                        outfile.getvalue(),
                                        singleblock,
                                        collinpred
                                        ));
           }
         else
           {
           fcpsplinesurfgaussian.push_back(
           FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],
                                          fcconst_intercept,
                                          D.getCol(j1),
                                          D.getCol(j2),
                                          D.getCol(j3),
                                          title,
                                          nrknots,degree,po,
                                          lambda,
                                          gridsize,
                                          type,
                                          pathnonp,
                                          pathres,
                                          outfile.getvalue(),
                                          singleblock,
                                          center,
                                          collinpred
                                          ));

           }

        if (constlambda.getvalue() == true)
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

        if ( (main1==1) && (main2==1) )
          {

          ST::string pathnonpt;
          ST::string pathrest;
          ST::string titlet;

          make_paths(collinpred,pathnonpt,pathrest,titlet,help,"",
                 "_pspline_total.raw","_pspline_total.res","_pspline_total");

          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].
          init_maineffects(mainp1,mainp2,pathnonpt,pathrest);
          mainp1->set_interaction();
          mainp2->set_interaction();
          }
        else if ( (main1==0) && (main2==0) )
          {
          }
        else
          {
          // FEHLERMELDUNG
          }

        vector<ST::string> na;
        if(varcoeff==false)
          {
          na.push_back(terms[i].varnames[0]);
          na.push_back(terms[i].varnames[1]);
          }
        else if (varcoeff==true)
          {
          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          na.push_back(terms[i].varnames[2]);
          }

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);

        make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline_var.raw","_pspline_var.res","_pspline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

        bool rowwise = false;
        if(terms[i].options[0] == "tpspline2dimband")
          rowwise = true;

        if ((terms[i].options[0] == "tpspline2dimrw1") ||
            (terms[i].options[0] == "tpspline2dimband")
           )
          {

          make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline_tvar.raw","_pspline_tvar.res","_pspline_tvariance");

          unsigned v = nu.getvalue();
          f = (terms[i].options[16]).strtolong(h);
          blocksize = unsigned(h);

          fctvariance2dim.push_back(FULLCOND_tvariance2dim(&generaloptions[generaloptions.size()-1],
          &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],v,title,
          pathnonp,pathres,blocksize,rowwise));

          }

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

        if (constlambda.getvalue() == false)
          {
          if(terms[i].options[15]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        if (terms[i].options[0] == "tpspline2dimrw1" ||
           (terms[i].options[0] == "tpspline2dimband")
           )
          {
          fctvariance2dim[fctvariance2dim.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fctvariance2dim[fctvariance2dim.size()-1]);
          }

        }
      else
        {
        // NONGAUSSIAN CASE
        if(proposal=="iwls" || proposal=="iwlsmode")
          {
          // IWLS

          IWLS_pspline * mainp1=NULL;
          IWLS_pspline * mainp2=NULL;
          unsigned main1=0;
          unsigned main2=0;

          bool iwlsmode;
          if(proposal=="iwlsmode")
            iwlsmode = true;
          else
            iwlsmode = false;
          f = (terms[i].options[12]).strtolong(h);
          unsigned updateW;
          updateW = unsigned(h);
          bool updatetau;
          if(terms[i].options[13] == "false" || constlambda.getvalue() == true)
            updatetau = false;
          else
            updatetau = true;
          double fstart;
            f = (terms[i].options[14]).strtodouble(fstart);

          unsigned j;
          for (j=0;j<fciwlspspline.size();j++)
            {
            if  ( ((fciwlspspline[j].get_datanames()).size() == 1) &&
                   (fciwlspspline[j].get_datanames()[0] == terms[i].varnames[0]) &&
                   fciwlspspline[j].get_col() == collinpred
                )
              {
              mainp1 = &fciwlspspline[j];
              if(mainp1->get_gridsize() != gridsize)
                {
                outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp1->get_nrknots() != nrknots)
                {
                outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp1->get_degree() != degree)
                {
                outerror("ERROR: degree for interaction term and main effects must be the same\n");
                return true;
                }
              main1 ++;
              }


            if  ( ((fciwlspspline[j].get_datanames()).size() == 1) &&
                (fciwlspspline[j].get_datanames()[0] == terms[i].varnames[1]) &&
                fciwlspspline[j].get_col() == collinpred
                )
              {
              mainp2 = &fciwlspspline[j];
              if(mainp2->get_gridsize() != gridsize)
                {
                outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp2->get_nrknots() != nrknots)
                {
                outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp2->get_degree() != degree)
                {
                outerror("ERROR: degree for interaction term and main effects must be the same\n");
                return true;
                }
              main2 ++;
              }

            }


          if (varcoeff==false)
            {
            fcpsplinesurfgaussian.push_back( FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],
                                           distr[distr.size()-1],
                                           fcconst_intercept,
                                           D.getCol(j1),
                                           D.getCol(j2),
                                           iwlsmode,
                                           title,
                                           nrknots,degree,
                                           po,
                                           lambda,
                                           updateW,
                                           updatetau,
                                           fstart,
                                           a1,b1,
                                           gridsize,
                                           type,
                                           pathnonp,
                                           pathres,
                                           outfile.getvalue(),
                                           true,
                                           singleblock,
                                           collinpred
                                          )
                          );
            }
          else
            {

            fcpsplinesurfgaussian.push_back( FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],
                                           distr[distr.size()-1],
                                           fcconst_intercept,
                                           D.getCol(j1),
                                           D.getCol(j2),
                                           D.getCol(j3),
                                           title,
                                           nrknots,degree,
                                           po,
                                           lambda,
                                           updateW,
                                           updatetau,
                                           fstart,
                                           a1,b1,
                                           gridsize,
                                           type,
                                           pathnonp,
                                           pathres,
                                           outfile.getvalue(),
                                           true,
                                           singleblock,
                                           center,
                                           collinpred
                                          )
                          );
            }

          if (constlambda.getvalue() == true)
            fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

          if ( (main1==1) && (main2==1) )
            {

            ST::string pathnonpt;
            ST::string pathrest;
            ST::string titlet;

            make_paths(collinpred,pathnonpt,pathrest,titlet,help,"",
                 "_pspline_total.raw","_pspline_total.res","_pspline_total");

            fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].
            init_maineffects(mainp1,mainp2,pathnonpt,pathrest);
            mainp1->set_interaction();
            mainp2->set_interaction();
            }
          else if ( (main1==0) && (main2==0) )
            {
            }
          else
            {
            // FEHLERMELDUNG
            }


          vector<ST::string> na;
          na.push_back(terms[i].varnames[0]);
          na.push_back(terms[i].varnames[1]);

          if (varcoeff==true)
            na.push_back(terms[i].varnames[1]);


          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);

          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline_var.raw","_pspline_var.res","_pspline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

          bool rowwise = false;
          if(terms[i].options[0] == "tpspline2dimband")
            rowwise = true;

          if ((terms[i].options[0] == "tpspline2dimrw1") ||
              (terms[i].options[0] == "tpspline2dimband")
             )
            {

            make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline_tvar.raw","_pspline_tvar.res","_pspline_tvariance");

            unsigned v = nu.getvalue();
            f = (terms[i].options[16]).strtolong(h);
            blocksize = unsigned(h);

            fctvariance2dim.push_back(FULLCOND_tvariance2dim(&generaloptions[generaloptions.size()-1],
            &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],v,title,
            pathnonp,pathres,blocksize,rowwise));

            }

          if (constlambda.getvalue() == false)
            {
            if(terms[i].options[15]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            if(updatetau)
              fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();
            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          if (terms[i].options[0] == "tpspline2dimrw1")
            {
            fctvariance2dim[fctvariance2dim.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fctvariance2dim[fctvariance2dim.size()-1]);
            }

          }
        else
          {
          // CONDITIONAL PRIOR

          if (terms[i].options[0] == "tpspline2dimrw1")
            {
            outerror("ERROR: '" + terms[i].options[0] + "' not available\n");
            return true;
            }

          FULLCOND_pspline * mainp1=NULL;
          FULLCOND_pspline * mainp2=NULL;
          unsigned main1=0;
          unsigned main2=0;

          unsigned j;
          for (j=0;j<fcpspline.size();j++)
            {
            if  ( ((fcpspline[j].get_datanames()).size() == 1) &&
                 (fcpspline[j].get_datanames()[0] == terms[i].varnames[0]) &&
                  fcpspline[j].get_col() == collinpred
                )
              {
              mainp1 = &fcpspline[j];
              if(mainp1->get_gridsize() != gridsize)
                {
                outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp1->get_nrknots() != nrknots)
                {
                outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp1->get_degree() != degree)
                {
                outerror("ERROR: degree for interaction term and main effects must be the same\n");
                return true;
                }
              main1 ++;
              }


            if  ( ((fcpspline[j].get_datanames()).size() == 1) &&
                (fcpspline[j].get_datanames()[0] == terms[i].varnames[1]) &&
                fcpspline[j].get_col() == collinpred
                )
              {
              mainp2 = &fcpspline[j];
              if(mainp2->get_gridsize() != gridsize)
                {
                outerror("ERROR: gridsize for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp2->get_nrknots() != nrknots)
                {
                outerror("ERROR: number of knots for interaction term and main effects must be the same\n");
                return true;
                }
              if(mainp2->get_degree() != degree)
                {
                outerror("ERROR: degree for interaction term and main effects must be the same\n");
                return true;
                }
              main2 ++;
              }

            }

          fcpsplinesurf.push_back( FULLCOND_pspline_surf(&generaloptions[generaloptions.size()-1],
                                         distr[distr.size()-1],
                                         fcconst_intercept,
                                         D.getCol(j1),
                                         D.getCol(j2),
                                         title,
                                         nrknots,degree,
                                         po,
                                         min,max,
                                         lambda,
                                         gridsize,
                                         type,
                                         pathnonp,
                                         pathres,
                                         outfile.getvalue(),
                                         collinpred
                                        )
                        );

          if (constlambda.getvalue() == true)
            fcpsplinesurf[fcpsplinesurf.size()-1].set_lambdaconst(lambda);

          if ( (main1==1) && (main2==1) )
            {

            ST::string pathnonpt;
            ST::string pathrest;
            ST::string titlet;

            make_paths(collinpred,pathnonpt,pathrest,titlet,help,"",
                   "_pspline_total.raw","_pspline_total.res","_pspline_total");

            fcpsplinesurf[fcpsplinesurf.size()-1].
            init_maineffects(mainp1,mainp2,pathnonpt,pathrest);
            mainp1->set_interaction();
            mainp2->set_interaction();
            }
          else if ( (main1==0) && (main2==0) )
            {
            }
          else
            {
            // FEHLERMELDUNG
            }

          vector<ST::string> na;
          na.push_back(terms[i].varnames[0]);
          na.push_back(terms[i].varnames[1]);

          fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);

          fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,help,"",
                 "_pspline_var.raw","_pspline_var.res","_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                       &fcpsplinesurf[fcpsplinesurf.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );


          if (constlambda.getvalue() == false)
            {
            if(terms[i].options[15]=="true")
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

/*
bool bayesreg::create_varcoeff_interactionspspline(const unsigned & collinpred=0)
  {

  }
*/

bool bayesreg::create_geospline(const unsigned & collinpred)
  {

  ST::string proposal;

  long h;
  double lambda;
  double a1,b1;
//  bool reduced,singleblock;
  bool singleblock;
  unsigned min,max,nrknots,degree;
  int gridsize;
  int f;

  unsigned i;
  int j;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpgeospline.checkvector(terms,i) == true )
      {

      MCMC::fieldtype type;
      if ((terms[i].options[0] == "geospline") || (terms[i].options[0] == "geosplinerw1"))
        type = MCMC::mrflinear;
      else if (terms[i].options[0] == "geosplinerw2")
        type = MCMC::mrfquadratic8;

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

//      if (terms[i].options[6] == "false")
//        reduced = false;
//      else
//        reduced = true;

      mapobject * mapp;                           // pointer to mapobject

      int objpos = findstatobject(*statobj,terms[i].options[7],"map");

      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[7] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[7] + " is not a map object\n");
        return true;
        }

      MAP::map m = mapp->getmap();

      if(!m.centroids_existing())
        {
        outerror("ERROR: map object does not contain centroids\n");
        return true;
        }

      if (terms[i].options[8] == "false")
        singleblock = false;
      else
        singleblock = true;

      f = (terms[i].options[9]).strtodouble(a1);
      f = (terms[i].options[10]).strtodouble(b1);

      proposal = terms[i].options[11];

      gridsize = -1;

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      if (f==1)
        return true;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline.raw","_geospline.res","_geospline");

      if (check_gaussian(collinpred))
        {

        fcpsplinesurfgaussian.push_back(
        FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],
                                      distr[distr.size()-1],fcconst_intercept,
                                      D.getCol(j),m,terms[i].options[7],
                                      title,
                                      nrknots,degree,po,
                                      lambda,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      singleblock,
                                      collinpred
                                      ));

        if (constlambda.getvalue() == true)
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

        vector<ST::string> na;
        na.push_back(terms[i].varnames[0]);

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline_var.raw","_geospline_var.res","_geospline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                      &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

        if (constlambda.getvalue() == false)
          {
          if(terms[i].options[15]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        }
      else
        {
        // NONGAUSSIAN CASE
        if(proposal=="iwls" || proposal=="iwlsmode")
          {
          // IWLS
          bool iwlsmode;
          if(proposal=="iwlsmode")
            iwlsmode = true;
          else
            iwlsmode = false;
          f = (terms[i].options[12]).strtolong(h);
          unsigned updateW;
          updateW = unsigned(h);
          bool updatetau;
          if(terms[i].options[13] == "false" || constlambda.getvalue() == true)
            updatetau = false;
          else
            updatetau = true;
          double fstart;
            f = (terms[i].options[14]).strtodouble(fstart);

          fcpsplinesurfgaussian.push_back(
          FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],
                                      distr[distr.size()-1],fcconst_intercept,
                                      D.getCol(j),m,terms[i].options[7],
                                      iwlsmode,
                                      title,
                                      nrknots,degree,po,
                                      lambda,
                                      updateW,
                                      updatetau,
                                      fstart,
                                      a1,b1,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      true,
                                      singleblock,
                                      collinpred
                                      ));

          if (constlambda.getvalue() == true)
            fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

          vector<ST::string> na;
          na.push_back(terms[i].varnames[0]);
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline_var.raw","_geospline_var.res","_geospline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                        &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
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
            if(terms[i].options[15]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          }
        else
          {
          // CONDITIONAL PRIOR
          fcpsplinesurf.push_back(
          FULLCOND_pspline_surf(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],
                                      fcconst_intercept,
                                      D.getCol(j),m,terms[i].options[7],
                                      title,
                                      nrknots,degree,po,
                                      min,max,
                                      lambda,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      collinpred
                                      ));

          if (constlambda.getvalue() == true)
            fcpsplinesurf[fcpsplinesurf.size()-1].set_lambdaconst(lambda);

          vector<ST::string> na;
          na.push_back(terms[i].varnames[0]);
          fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
          fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_geospline_var.raw","_geospline_var.res","_geospline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                        &fcpsplinesurf[fcpsplinesurf.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

          if (constlambda.getvalue() == false)
            {
            if(terms[i].options[15]=="true")
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


bool bayesreg::create_varcoeff_geospline(const unsigned & collinpred)
  {

  ST::string proposal;

  long h;
  double lambda;
  double a1,b1;
//  bool reduced,singleblock;
  bool singleblock;
  unsigned min,max,nrknots,degree;
  int gridsize;
  bool center;
  int f;

  unsigned i;
  int j1,j2;
  for(i=0;i<terms.size();i++)
    {
    if ( nonpvarcoeffgeospline.checkvector(terms,i) == true )
      {

      MCMC::fieldtype type = MCMC::mrflinear;

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv); // interaction var
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

//      if (terms[i].options[6] == "false")
//        reduced = false;
//      else
//        reduced = true;

      mapobject * mapp;                           // pointer to mapobject

      int objpos = findstatobject(*statobj,terms[i].options[7],"map");

      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          outerror("ERROR: map object " + terms[i].options[7] + " is not existing\n");
        else
          outerror("ERROR: " + terms[i].options[7] + " is not a map object\n");
        return true;
        }

      MAP::map m = mapp->getmap();

      if(!m.centroids_existing())
        {
        outerror("ERROR: map object does not contain centroids\n");
        return true;
        }

      if (terms[i].options[8] == "false")
        singleblock = false;
      else
        singleblock = true;

      if (terms[i].options[16] == "false")
        center = false;
      else
        center = true;

      f = (terms[i].options[9]).strtodouble(a1);
      f = (terms[i].options[10]).strtodouble(b1);

      proposal = terms[i].options[11];

      gridsize = -1;

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      if(f==1)
        return false;

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],terms[i].varnames[0],
                 "_geospline.raw","_geospline.res","_geospline");

      if (check_gaussian(collinpred))
        {

        fcpsplinesurfgaussian.push_back(
        FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],
                                      fcconst_intercept,
                                      D.getCol(j1),D.getCol(j2),m,terms[i].options[7],
                                      title,
                                      nrknots,degree,po,
                                      lambda,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      singleblock,
                                      center,
                                      collinpred
                                      ));

        if(constlambda.getvalue() == true)
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

        vector<ST::string> na;
        na.push_back(terms[i].varnames[1]);
        na.push_back(terms[i].varnames[0]);

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],terms[i].varnames[1],
                   "_geospline_var.raw","_geospline_var.res","_geospline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                      &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );

        fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

        if (constlambda.getvalue() == false)
          {
          if(terms[i].options[15]=="true")
            fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        }
      else
        {
        // NONGAUSSIAN CASE
        if(proposal=="iwls" || proposal=="iwlsmode")
          {
          // IWLS
          bool iwlsmode;
          if(proposal=="iwlsmode")
            iwlsmode = true;
          else
            iwlsmode = false;
          f = (terms[i].options[12]).strtolong(h);
          unsigned updateW;
          updateW = unsigned(h);
          bool updatetau;
          if(terms[i].options[13] == "false")
            updatetau = false;
          else
            updatetau = true;
          double fstart;
            f = (terms[i].options[14]).strtodouble(fstart);

          fcpsplinesurfgaussian.push_back(
          FULLCOND_pspline_surf_gaussian(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],
                                      fcconst_intercept,
                                      D.getCol(j1),D.getCol(j2),m,terms[i].options[7],
                                      title,
                                      nrknots,degree,po,
                                      lambda,
                                      updateW,
                                      updatetau,
                                      fstart,
                                      a1,b1,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      iwlsmode,
                                      singleblock,
                                      center,
                                      collinpred
                                      ));

          if (constlambda.getvalue() == true)
            fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_lambdaconst(lambda);

          vector<ST::string> na;
          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].init_names(na);
          fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline_var.raw","_geospline_var.res","_geospline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                        &fcpsplinesurfgaussian[fcpsplinesurfgaussian.size()-1],
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
            if(terms[i].options[15]=="true")
              fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          }
        else
          {
          // CONDITIONAL PRIOR
          fcpsplinesurf.push_back(
          FULLCOND_pspline_surf(&generaloptions[generaloptions.size()-1],distr[distr.size()-1],
                                      fcconst_intercept,
                                      D.getCol(j1),D.getCol(j2),m,terms[i].options[7],
                                      title,
                                      nrknots,degree,po,
                                      min,max,
                                      lambda,
                                      gridsize,
                                      type,
                                      pathnonp,
                                      pathres,
                                      collinpred
                                      ));

          if (constlambda.getvalue() == true)
            fcpsplinesurf[fcpsplinesurf.size()-1].set_lambdaconst(lambda);

          vector<ST::string> na;
          na.push_back(terms[i].varnames[1]);
          na.push_back(terms[i].varnames[0]);
          fcpsplinesurf[fcpsplinesurf.size()-1].init_names(na);
          fcpsplinesurf[fcpsplinesurf.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpsplinesurf[fcpsplinesurf.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_geospline_var.raw","_geospline_var.res","_geospline_variance");

          fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                        &fcpsplinesurf[fcpsplinesurf.size()-1],
                                       distr[distr.size()-1],
                                       a1,
                                       b1,
                                       title,pathnonp,pathres,
                                       false,collinpred)
                           );


          if (constlambda.getvalue() == false)
            {
            if(terms[i].options[15]=="true")
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



bool bayesreg::create_baseline(const unsigned & collinpred)
  {

  ST::string proposal;

  long h;
  unsigned min,max,degree,nrknots;
  double lambda,a1,b1;
  bool ub, wb;
  int gridsize;
  int f;
  // NEW FOR PARTIALLIKELIHOOD
  bool partlik;

  unsigned i;
  int j;
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

      if (terms[i].options[12] == "false")
        wb = false;
      else
        wb = true;

     // NEW FOR PARTIALLIKELIHOOD
     if (terms[i].options[15] == "false")
       {
       partlik = false;
       }
     else
       {
       partlik = true;
       }
      if (f==1)
        return true;

      datamatrix beg;
      if (begin.getvalue() == "")
        beg = datamatrix(1,1);
      else
        beg = D.getCol(begpos);

      // -------------end: reading options, term information -------------------

      //--------- creating path for samples and and results, creating title ----

      make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_logbaseline.raw","_logbaseline.res","_logbaseline");

      //----- end: creating path for samples and and results, creating title ---


      if(proposal == "cp")
        {
        fcbaseline.push_back( pspline_baseline(&generaloptions[generaloptions.size()-1],
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
                                                beg,
                                                wb,
                                                partlik    // NEW FOR PARTIALLIKELIHOOD
                                               )
                             );

        if (constlambda.getvalue() == true)
          fcbaseline[fcbaseline.size()-1].set_lambdaconst(lambda);

        fcbaseline[fcbaseline.size()-1].init_name(terms[i].varnames[0]);
        fcbaseline[fcbaseline.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcbaseline[fcbaseline.size()-1]);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_logbaseline_var.raw","_logbaseline_var.res","_logbaseline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                &fcbaseline[fcbaseline.size()-1],
                                distr[distr.size()-1],a1,b1,
                                title,pathnonp,pathres,ub,collinpred)
                                );

        }
      else if(proposal == "iwls")
        {

/*        fcbaselineiwls.push_back( IWLS_baseline(&generaloptions[generaloptions.size()-1],
                                                distr[distr.size()-1],
                                                fcconst_intercept,
                                                D.getCol(j),
                                                false,
                                                nrknots,
                                                degree,
                                                po,
                                                lambda,
                                                type,
                                                "unrestricted",
                                                1,
                                                false,
                                                2,
                                                a1,
                                                b1,
                                                title,
                                                pathnonp,
                                                pathres,
                                                false,
                                                gridsize,
                                                false,
                                                collinpred,
                                                beg
                                               )
                             );

        if (constlambda.getvalue() == true)
          fcbaselineiwls[fcbaselineiwls.size()-1].set_lambdaconst(lambda);

        fcbaselineiwls[fcbaselineiwls.size()-1].init_name(terms[i].varnames[0]);
        fcbaselineiwls[fcbaselineiwls.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcbaselineiwls[fcbaselineiwls.size()-1]);

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                   "_logbaseline_var.raw","_logbaseline_var.res","_logbaseline_variance");

        fcvarnonp.push_back(FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                                &fcbaselineiwls[fcbaselineiwls.size()-1],
                                distr[distr.size()-1],a1,b1,
                                title,pathnonp,pathres,ub,collinpred)
                                );   */
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



bool bayesreg::create_random_rw1rw2(const unsigned & collinpred)
  {

  double lambda_r;
  double a_r,b_r;
//  bool updatetau_r;
  ST::string proposal_r;

  long h;
  double lambda,a1,b1,alpha;
//  bool updatetau;
  double ftune;
  unsigned updateW;
  ST::string proposal;
  int f;

  unsigned i;
  int j1,j2;

  for(i=0;i<terms.size();i++)
    {
    if ( randomrw.checkvector(terms,i) == true)
      {

      // -------------- reading options, term information ----------------------

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);  // random effect
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);  // nonlinear function

      MCMC::fieldtype type;

      if (terms[i].options[0] == "random_rw1")
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      f = (terms[i].options[1]).strtodouble(lambda_r);

      f = (terms[i].options[2]).strtodouble(a_r);

      f = (terms[i].options[3]).strtodouble(b_r);

      proposal_r = terms[i].options[4];


      f = (terms[i].options[8]).strtodouble(lambda);

      f = (terms[i].options[9]).strtodouble(a1);

      f = (terms[i].options[10]).strtodouble(b1);

      proposal = terms[i].options[11];

      f = (terms[i].options[12]).strtolong(h);
      updateW = unsigned(h);

 //     if (terms[i].options[13] == "true")
 //       updatetau=true;
 //     else
 //       updatetau=false;

      f = (terms[i].options[14]).strtodouble(ftune);

      f = (terms[i].options[19]).strtodouble(alpha);

      // -------------- reading options, term information ----------------------

      // -------- creating paths for samples and results, titles ---------------


      // -------- end: creating paths for samples and results, titles ----------


      //-------------------- gaussian response, etc. ---------------------------
      if ( (check_gaussian(collinpred)) || (check_iwls(true,collinpred)) )
        {

        // Include nonlinear function f

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
        "f_"+terms[i].varnames[0]+"_random",
                 "_rw_mult1.raw","_rw_mult1.res","_rw_mult1");

        datamatrix eins(D.rows(),1,1);
        datamatrix null(D.rows(),1,0);

        fcnonpgaussian.push_back(
        FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
        distr[distr.size()-1],D.getCol(j2),eins,fcconst_intercept,
        unsigned(maxint.getvalue()),type,title,pathnonp,pathres,collinpred,
        lambda,false));

        fcnonpgaussian[fcnonpgaussian.size()-1].init_name(
        terms[i].varnames[1]);

        fcnonpgaussian[fcnonpgaussian.size()-1].set_changingweight();
//        fcnonpgaussian[fcnonpgaussian.size()-1].set_center(true);


        if (constlambda.getvalue() == true)
          {
          fcnonpgaussian[fcnonpgaussian.size()-1].set_lambdaconst(lambda);
          }


        // Include variance of nonlinear f

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
        "f_"+terms[i].varnames[0]+"_random",
                 "_rw_mult1_var.raw","_rw_mult1_var.res","_rw_mult1_var");

        fcvarnonp.push_back(
        FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
        &fcnonpgaussian[fcnonpgaussian.size()-1],distr[distr.size()-1],
        a1,b1,title,pathnonp,pathres,false,collinpred));

      // Include random effect

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
        "f_"+terms[i].varnames[0]+"_random",
                 "_rw_mult2.raw","_rw_mult2.res","_rw_mult2");

        fcrandomgaussian.push_back(
        FULLCOND_random_gaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        null,
                                                        D.getCol(j1),
                                                        title,
                                                        pathnonp,
                                                        pathres,"",
                                                        lambda_r,
                                                        false,collinpred
                                                        )
                          );

        vector<ST::string> na;
        na.push_back(terms[i].varnames[0]);
        na.push_back(terms[i].varnames[1]);
        fcrandomgaussian[fcrandomgaussian.size()-1].init_names(na);


        fcrandomgaussian[fcrandomgaussian.size()-1].set_notransform();

        fcrandomgaussian[fcrandomgaussian.size()-1].set_changingweight();

       // Include first fcmult

       fcmult.push_back(FULLCOND_mult(&generaloptions[generaloptions.size()-1],
       distr[distr.size()-1],&fcrandomgaussian[fcrandomgaussian.size()-1],
       &fcnonpgaussian[fcnonpgaussian.size()-1],true,"","","",collinpred));

       // Include second fcmult

       fcmult.push_back(FULLCOND_mult(&generaloptions[generaloptions.size()-1],
       distr[distr.size()-1],&fcrandomgaussian[fcrandomgaussian.size()-1],
       &fcnonpgaussian[fcnonpgaussian.size()-1],false,"","","",collinpred));


        // Reinhängen in fullcond

        fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

        fcmult[fcmult.size()-2].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcmult[fcmult.size()-2]);


        fcrandomgaussian[fcrandomgaussian.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true)
          {
          fcrandomgaussian[fcrandomgaussian.size()-1].set_lambdaconst(lambda_r);
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);
          }
        else
          {
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
          "f_"+terms[i].varnames[0]+"_random",
                   "_rw_mult2_var.raw","_rw_mult2_var.res","_rw_mult2_var");

          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandomgaussian[fcrandomgaussian.size()-1],
                            distr[distr.size()-1],a_r,b_r,title,pathnonp,
                            pathres,false,collinpred));

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        fcmult[fcmult.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcmult[fcmult.size()-1]);


          //------------------- end: gaussian response, etc. -------------------
        }


      } // end: if ( nonprw1rw2.checkvector(terms,i) == true )

    }

  return false;

  }


bool bayesreg::create_spatial_rw1rw2(const unsigned & collinpred)
  {

  double lambda_s;
  double a_s,b_s;
  double alpha_s;
  bool center_s;
//  bool updatetau_s;
//  bool Laplace_s;
//  unsigned updateW_s;
  double ftune_s;
  ST::string proposal_s;

//  long h;
  double lambda,a1,b1,alpha;
//  bool updatetau;
  double ftune;
//  unsigned updateW;
  ST::string proposal;
  int f;

  unsigned i;
  int j1,j2;

  for(i=0;i<terms.size();i++)
    {
    if ( spatialrw.checkvector(terms,i) == true)
      {

      // -------------- reading options, term information ----------------------

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);  // spatial effect
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);  // nonlinear function

      MCMC::fieldtype type;

      if (terms[i].options[0] == "spatial_rw1")
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      mapobject * mapp;                           // pointer to mapobject

      int objpos = findstatobject(*statobj,terms[i].options[1],"map");

      if (objpos >= 0)
        {
        statobject * s = statobj->at(objpos);
        mapp = dynamic_cast<mapobject*>(s);
        }
      else
        {
        if (objpos == -1)
          {
          if ((terms[i].options[1] == "") || (terms[i].options[1] == " "))
            outerror("ERROR: map object must be specified to estimate a spatial effect\n");
          else
            outerror("ERROR: map object " + terms[i].options[1] + " is not existing\n");
          }
        else
          outerror("ERROR: " + terms[i].options[1] + " is not a map object\n");
        return true;
        }

      MAP::map m = mapp->getmap();
      bool isconnected = m.isconnected();
      if (isconnected==false)
        {
        outerror("ERROR: map is disconnected, spatial effect cannot be estimated\n");
        return true;
        }

      f = (terms[i].options[2]).strtodouble(lambda_s);

      f = (terms[i].options[3]).strtodouble(a_s);

      f = (terms[i].options[4]).strtodouble(b_s);

      proposal_s = terms[i].options[5];


//      f = (terms[i].options[6]).strtolong(h);
//      updateW_s = unsigned(h);

//      if (terms[i].options[7] == "true")
//        updatetau_s=true;
//      else
//        updatetau_s=false;

      f = (terms[i].options[8]).strtodouble(ftune_s);

//      if (terms[i].options[11] == "true")
//        Laplace_s=true;
//      else
//        Laplace_s=false;

      f = (terms[i].options[13]).strtodouble(alpha_s);

      if (f==1)
        return true;

      if (terms[i].options[15] == "true")
        center_s=true;
      else
        center_s=false;


      f = (terms[i].options[16]).strtodouble(lambda);

      f = (terms[i].options[17]).strtodouble(a1);

      f = (terms[i].options[18]).strtodouble(b1);

      proposal = terms[i].options[19];

 //     f = (terms[i].options[20]).strtolong(h);
 //     updateW = unsigned(h);

//      if (terms[i].options[21] == "true")
//        updatetau=true;
//      else
//        updatetau=false;

      f = (terms[i].options[22]).strtodouble(ftune);

      f = (terms[i].options[27]).strtodouble(alpha);

      // -------------- reading options, term information ----------------------

      // -------- creating paths for samples and results, titles ---------------


      // -------- end: creating paths for samples and results, titles ----------


      //-------------------- gaussian response, etc. ---------------------------
      if ( (check_gaussian(collinpred)) || (check_iwls(true,collinpred)) )
        {

        // Include nonlinear function f

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
        "f_"+terms[i].varnames[0]+"_spatial",
                 "_rw_mult1.raw","_rw_mult1.res","_rw_mult1");

        datamatrix eins(D.rows(),1,1);
        datamatrix null(D.rows(),1,0);

        fcnonpgaussian.push_back(
        FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
        distr[distr.size()-1],D.getCol(j2),eins,fcconst_intercept,
        unsigned(maxint.getvalue()),type,title,pathnonp,pathres,collinpred,
        lambda,false));

        fcnonpgaussian[fcnonpgaussian.size()-1].init_name(
        terms[i].varnames[1]);

        fcnonpgaussian[fcnonpgaussian.size()-1].set_changingweight();
//        fcnonpgaussian[fcnonpgaussian.size()-1].set_center(true);


        if (constlambda.getvalue() == true)
          {
          fcnonpgaussian[fcnonpgaussian.size()-1].set_lambdaconst(lambda);
          }


        // Include variance of nonlinear f

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
        "f_"+terms[i].varnames[0]+"_spatial",
                 "_rw_mult1_var.raw","_rw_mult1_var.res","_rw_mult1_var");

        fcvarnonp.push_back(
        FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
        &fcnonpgaussian[fcnonpgaussian.size()-1],distr[distr.size()-1],
        a1,b1,title,pathnonp,pathres,false,collinpred));

      // Include spatial effect

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
        "f_"+terms[i].varnames[0]+"_spatial",
                 "_rw_mult2.raw","_rw_mult2.res","_rw_mult2");

        fcnonpgaussian.push_back(
        FULLCOND_nonp_gaussian(&generaloptions[generaloptions.size()-1],
        distr[distr.size()-1],fcconst_intercept,m,terms[i].options[1],
        D.getCol(j1),null,title,pathnonp,pathres,collinpred,lambda_s,
        center_s));

        vector<ST::string> na;
        na.push_back(terms[i].varnames[1]);
        na.push_back(terms[i].varnames[0]);
        fcnonpgaussian[fcnonpgaussian.size()-1].init_names(na);

        fcnonpgaussian[fcnonpgaussian.size()-1].set_notransform();

        fcnonpgaussian[fcnonpgaussian.size()-1].set_changingweight();

//        fcnonpgaussian[fcnonpgaussian.size()-1].set_center(true);

       // Include first fcmult

       fcmult.push_back(FULLCOND_mult(&generaloptions[generaloptions.size()-1],
       distr[distr.size()-1],&fcnonpgaussian[fcnonpgaussian.size()-1],
       &fcnonpgaussian[fcnonpgaussian.size()-2],true,"","","",collinpred));

       // Include second fcmult

       fcmult.push_back(FULLCOND_mult(&generaloptions[generaloptions.size()-1],
       distr[distr.size()-1],&fcnonpgaussian[fcnonpgaussian.size()-1],
       &fcnonpgaussian[fcnonpgaussian.size()-2],false,"","","",collinpred));


        // Reinhängen in fullcond

        fcnonpgaussian[fcnonpgaussian.size()-2].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-2]);

        fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);

        fcmult[fcmult.size()-2].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcmult[fcmult.size()-2]);



        fcnonpgaussian[fcnonpgaussian.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true)
          {
          fcnonpgaussian[fcnonpgaussian.size()-1].set_lambdaconst(lambda_s);
          fullcond.push_back(&fcnonpgaussian[fcrandomgaussian.size()-1]);
          }
        else
          {
          fullcond.push_back(&fcnonpgaussian[fcnonpgaussian.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
          "f_"+terms[i].varnames[0]+"_spatial",
                   "_rw_mult2_var.raw","_rw_mult2_var.res","_rw_mult2_var");

          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcnonpgaussian[fcnonpgaussian.size()-1],
                            distr[distr.size()-1],a_s,b_s,title,pathnonp,
                            pathres,false,collinpred));

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        fcmult[fcmult.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcmult[fcmult.size()-1]);

          //------------------- end: gaussian response, etc. -------------------
        }


      } // end: if ( nonprw1rw2.checkvector(terms,i) == true )

    }

  return false;

  }



bool bayesreg::create_random_pspline(const unsigned & collinpred)
  {

  double lambda_r;
  double a_r,b_r;
//  bool updatetau_r;
  ST::string proposal_r;

  ST::string monotone;
  ST::string proposal;

  long h;
  unsigned degree,nrknots;
//  unsigned min, max;
  double lambda,a1,b1,alpha;
  bool ub,derivative,bsplinebasis;
//  bool diagtransform;
  int gridsize,contourprob;
  int f;
  ST::string test ="test";
//  double lowerknot=0;
//  double upperknot=0;

  unsigned i;
  int j1,j2;

  for(i=0;i<terms.size();i++)
    {
    if ( randompspline.checkvector(terms,i) == true )
      {

      j1 = terms[i].varnames[0].isinlist(modelvarnamesv);  // random effect
      j2 = terms[i].varnames[1].isinlist(modelvarnamesv);  // nonlinear function

      // --------------- reading options, term information ---------------------
      MCMC::fieldtype type;
      if ( (terms[i].options[0] == "random_psplinerw1") ||
           (terms[i].options[0] == "random_tpsplinerw1") ||
           (terms[i].options[0] == "random_psplinerw1vrw1") ||
           (terms[i].options[0] == "random_psplinerw1vrw2") )
        type = MCMC::RW1;
      else
        type = MCMC::RW2;

      f = (terms[i].options[1]).strtodouble(lambda_r);

      f = (terms[i].options[2]).strtodouble(a_r);

      f = (terms[i].options[3]).strtodouble(b_r);

      proposal_r = terms[i].options[4];


//      f = (terms[i].options[8]).strtolong(h);
//      min = unsigned(h);

//      f = (terms[i].options[9]).strtolong(h);
//      max = unsigned(h);

      f = (terms[i].options[10]).strtolong(h);
      degree = unsigned(h);

      f = (terms[i].options[11]).strtolong(h);
      nrknots = unsigned(h);

      f = (terms[i].options[12]).strtodouble(lambda);

      f = (terms[i].options[13]).strtodouble(a1);

      f = (terms[i].options[14]).strtodouble(b1);

      if (terms[i].options[15] == "false")
        ub = false;
      else
        ub = true;

      f = (terms[i].options[16]).strtolong(h);
      gridsize = unsigned(h);

      if (f==1)
        return true;

      MCMC::knotpos po;

      if (knots.getvalue() == "equidistant" && terms[i].options[36] == "equidistant")
        po = MCMC::equidistant;
      else
        po = MCMC::quantiles;

      proposal = terms[i].options[20];
      monotone = terms[i].options[21];

//      if (terms[i].options[25] == "false")
//        diagtransform = false;
//      else
//        diagtransform = true;

      if (terms[i].options[26] == "false")
        derivative = false;
      else
        derivative = true;

      if (terms[i].options[27] == "false")
        bsplinebasis = false;
      else
        bsplinebasis = true;

      f = (terms[i].options[28]).strtolong(h);
      contourprob = unsigned(h);

      f = (terms[i].options[34]).strtodouble(alpha);


      // -------------end: reading options, term information -------------------


      // ---------------------- gaussian response, etc. ------------------------
      if (check_gaussian(collinpred))
        {

        // Include nonlinear f

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
        "f_"+terms[i].varnames[0]+"_random",
                 "_pspline_mult1.raw","_pspline_mult1.res","_pspline_mult1");

        datamatrix eins(D.rows(),1,1);
        datamatrix null(D.rows(),1,0);

        fcpsplinegaussian.push_back(
        FULLCOND_pspline_gaussian(&generaloptions[generaloptions.size()-1],
        distr[distr.size()-1],fcconst_intercept,D.getCol(j2),eins,nrknots,degree,po,
        type,monotone,title,pathnonp,pathres,derivative,lambda,gridsize,false,
        collinpred));

        //
        datamatrix beta_0;
        if(terms[i].options[30]!="")
          {
          dataobject * datap;                           // pointer to datasetobject
          int objpos = findstatobject(*statobj,terms[i].options[23],"dataset");
          if (objpos >= 0)
            {
            statobject * s = statobj->at(objpos);
            datap = dynamic_cast<dataobject*>(s);
            if (datap->obs()==0 || datap->getVarnames().size()==0)
              {
              outerror("ERROR: dataset object " + terms[i].options[23] + " does not contain any data\n");
              return true;
              }
            else if (datap->getVarnames().size()>1)
              {
              outerror("ERROR: dataset object " + terms[i].options[23] + " contains more than one variable\n");
              return true;
              }
            }
          else
            {
            outerror("ERROR: dataset object " + terms[i].options[23] + " is not existing\n");
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

        if (bsplinebasis == true)
          fcpsplinegaussian[fcpsplinegaussian.size()-1].set_outbsplines();

        if (terms[i].options[33] == "true")
          fcpsplinegaussian[fcpsplinegaussian.size()-1].set_stationary(alpha);

        fcpsplinegaussian[fcnonpgaussian.size()-1].set_changingweight();

        fcpsplinegaussian[fcpsplinegaussian.size()-1].init_name(terms[i].varnames[0]);
        fcpsplinegaussian[fcpsplinegaussian.size()-1].set_fcnumber(fullcond.size());

        // Include variance of nonlinear f

        if(terms[i].options[0] != "psplinerw1vrw1" && terms[i].options[0]
           != "psplinerw1vrw2" && terms[i].options[0] != "psplinerw2vrw1"
           && terms[i].options[0] != "psplinerw2vrw2")
          {


          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
          "f_"+terms[i].varnames[0]+"_random",
                 "_pspline_mult1_var.raw","_pspline_mult1_var.res",
                 "_pspline_mult1_var");


          fcvarnonp.push_back(FULLCOND_variance_nonp(
          &generaloptions[generaloptions.size()-1],
          &fcpsplinegaussian[fcpsplinegaussian.size()-1],distr[distr.size()-1],
          a1,b1,title,pathnonp,pathres,ub,collinpred));

          if (constlambda.getvalue() == false)
            {

            if(terms[i].options[29]=="true")
               fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
            if(terms[i].options[31]=="true")
              {
              f = terms[i].options[32].strtolong(h);
              fcvarnonp[fcvarnonp.size()-1].set_discrete(unsigned(h));
              }
            bool alphafix = false;
            if (terms[i].options[35] == "true")
              alphafix = true;
            if (terms[i].options[36] == "true")
              fcvarnonp[fcvarnonp.size()-1].set_stationary(alpha,alphafix);

            }

          }



        /*
        if (terms[i].options[0] == "tpsplinerw1" || terms[i].options[0] == "tpsplinerw2")
          {
          fcpsplinegaussian[fcpsplinegaussian.size()-1].set_adaptiv();

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_tvar.raw","_pspline_tvar.res","_pspline_tvariance");

          unsigned v = nu.getvalue();

          fctvariance.push_back(FULLCOND_tvariance(&generaloptions[generaloptions.size()-1],
                                &fcpsplinegaussian[fcpsplinegaussian.size()-1],
                                v,title,pathnonp,pathres));

          fctvariance[fctvariance.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fctvariance[fctvariance.size()-1]);

          } // end: if (terms[i].options[0] == "tpsplinerw1" ...

        if ( terms[i].options[0] == "psplinerw1vrw1" || terms[i].options[0] == "psplinerw1vrw2" ||
             terms[i].options[0] == "psplinerw2vrw1" || terms[i].options[0] == "psplinerw2vrw2" )
          {
          fcpsplinegaussian[fcpsplinegaussian.size()-1].set_adaptiv();

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_vvar.raw","_pspline_vvar.res","_pspline_vvariance");


          f = (terms[i].options[10]).strtolong(h);
          unsigned minvar = unsigned(h);
          f = (terms[i].options[11]).strtolong(h);
          unsigned maxvar = unsigned(h);
          double startv;
          f = (terms[i].options[12]).strtodouble(startv);
          MCMC::fieldtype typevar;
          if (terms[i].options[0] == "psplinerw1vrw1" ||
          terms[i].options[0] == "psplinerw2vrw1")
            typevar = MCMC::RW1;
          else
            typevar = MCMC::RW2;

          fcadaptiv.push_back(FULLCOND_adaptiv(&generaloptions[generaloptions.size()-1],
                              &fcpsplinegaussian[fcpsplinegaussian.size()-1],
                              typevar,title,a1,b1,ub,startv,minvar,maxvar,pathnonp,pathres));

          fcadaptiv[fcadaptiv.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcadaptiv[fcadaptiv.size()-1]);

          } // end: if (terms[i].options[0] == "psplinerw1vrw1" ...
          */


      // Include random effect

        make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
        "f_"+terms[i].varnames[0]+"_random",
                 "_pspline_mult2.raw","_pspline_mult2.res","_pspline_mult2");

        fcrandomgaussian.push_back(
        FULLCOND_random_gaussian(&generaloptions[generaloptions.size()-1],
                                                        distr[distr.size()-1],
                                                        fcconst_intercept,
                                                        null,
                                                        D.getCol(j1),
                                                        title,
                                                        pathnonp,
                                                        pathres,"",
                                                        lambda_r,
                                                        false,collinpred
                                                        )
                          );

        vector<ST::string> na;
        na.push_back(terms[i].varnames[0]);
        na.push_back(terms[i].varnames[1]);
        fcrandomgaussian[fcrandomgaussian.size()-1].init_names(na);

        fcrandomgaussian[fcrandomgaussian.size()-1].set_notransform();

       // Include first fcmult

       fcmult.push_back(FULLCOND_mult(&generaloptions[generaloptions.size()-1],
       distr[distr.size()-1],&fcrandomgaussian[fcrandomgaussian.size()-1],
       &fcpsplinegaussian[fcpsplinegaussian.size()-1],true,"","","",collinpred));

       // Include second fcmult

       fcmult.push_back(FULLCOND_mult(&generaloptions[generaloptions.size()-1],
       distr[distr.size()-1],&fcrandomgaussian[fcrandomgaussian.size()-1],
       &fcpsplinegaussian[fcpsplinegaussian.size()-1],false,"","","",collinpred));


        // Reinhängen in fullcond

        fcpsplinegaussian[fcpsplinegaussian.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcpsplinegaussian[fcpsplinegaussian.size()-1]);

        if (constlambda.getvalue() == false)
          {
          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        fcmult[fcmult.size()-2].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcmult[fcmult.size()-2]);


        fcrandomgaussian[fcrandomgaussian.size()-1].set_fcnumber(fullcond.size());

        if (constlambda.getvalue() == true)
          {
          fcrandomgaussian[fcrandomgaussian.size()-1].set_lambdaconst(lambda_r);
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);
          }
        else
          {
          fullcond.push_back(&fcrandomgaussian[fcrandomgaussian.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[1],
          "f_"+terms[i].varnames[0]+"_random",
                   "_pspline_mult2_var.raw","_pspline_mult2_var.res",
                   "_pspline_mult2_var");

          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
                            &fcrandomgaussian[fcrandomgaussian.size()-1],
                            distr[distr.size()-1],a_r,b_r,title,pathnonp,
                            pathres,false,collinpred));

          fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
          }

        fcmult[fcmult.size()-1].set_fcnumber(fullcond.size());
        fullcond.push_back(&fcmult[fcmult.size()-1]);



        // ---------------------- end: gaussian response, etc. -----------------
        }
      else
        {
        // -------------------- non-gaussian response, etc. --------------------

        /*

        if (proposal == "cp")
          {

          if (terms[i].options[26] == "true")
            {
            outerror("ERROR: option 'alpha' not available\n");
            return true;
            }

          if (terms[i].options[0] == "tpsplinerw1" || terms[i].options[0] == "tpsplinerw2" ||
              terms[i].options[0] == "psplinerw1vrw1" || terms[i].options[0] == "psplinerw1vrw2" ||
              terms[i].options[0] == "psplinerw2vrw1" || terms[i].options[0] == "psplinerw2vrw2" )
            {
            outerror("ERROR: '" + terms[i].options[0] + "' not available\n");
            return true;
            }

          fcpspline.push_back( FULLCOND_pspline(&generaloptions[generaloptions.size()-1],
          distr[distr.size()-1],fcconst_intercept,meandata,nrknots,degree,po,
          lambda,min,max,type,title,pathnonp,pathres,derivative,lowerknot,
          upperknot,gridsize,
          collinpred));

          if (constlambda.getvalue() == true)
            fcpspline[fcpspline.size()-1].set_lambdaconst(lambda);

          fcpspline[fcpspline.size()-1].init_name(terms[i].varnames[0]);
          fcpspline[fcpspline.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fcpspline[fcpspline.size()-1]);

          make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_var.raw","_pspline_var.res","_pspline_variance");


          fcvarnonp.push_back(
          FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
          &fcpspline[fcpspline.size()-1],distr[distr.size()-1],a1,b1,title,
          pathnonp,pathres,ub,collinpred));

          if (constlambda.getvalue() == false)
            {

            if(terms[i].options[22]=="true")
                fcvarnonp[fcvarnonp.size()-1].set_uniformprior();

            fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
            }

          }
        else // iwls
          {

          bool iwlsmode;
          if(proposal == "iwlsmode")
            iwlsmode = true;
          else
            iwlsmode = false;
          f = (terms[i].options[15]).strtolong(h);
          unsigned updateW;
          updateW = unsigned(h);
          bool updatetau;
          if(terms[i].options[16] == "false" || constlambda.getvalue() == true)
            updatetau = false;
          else
            updatetau = true;

          double fstart;
            f = (terms[i].options[17]).strtodouble(fstart);

          fciwlspspline.push_back(
          IWLS_pspline(&generaloptions[generaloptions.size()-1],
          distr[distr.size()-1],fcconst_intercept,meandata,iwlsmode,nrknots,degree,po,
          lambda,type,monotone,updateW,updatetau,fstart,a1,b1,title,pathnonp,
          pathres,derivative,gridsize,diagtransform,lowerknot,upperknot,
          collinpred));

          if (constlambda.getvalue() == true)
            fciwlspspline[fciwlspspline.size()-1].set_lambdaconst(lambda);

          if (terms[i].options[26] == "true")
            fciwlspspline[fciwlspspline.size()-1].set_stationary(alpha);

          fciwlspspline[fciwlspspline.size()-1].init_name(terms[i].varnames[0]);
          fciwlspspline[fciwlspspline.size()-1].set_fcnumber(fullcond.size());
          fullcond.push_back(&fciwlspspline[fciwlspspline.size()-1]);

          if(terms[i].options[0] != "psplinerw1vrw1" && terms[i].options[0] != "psplinerw1vrw2" &&
             terms[i].options[0] != "psplinerw2vrw1" && terms[i].options[0] != "psplinerw2vrw2")
            {

            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_var.raw","_pspline_var.res","_pspline_variance");

            fcvarnonp.push_back(
            FULLCOND_variance_nonp(&generaloptions[generaloptions.size()-1],
            &fciwlspspline[fciwlspspline.size()-1],distr[distr.size()-1],a1,b1,
            title,pathnonp,pathres,ub,collinpred));

            bool alphafix = false;
            if (terms[i].options[28] == "true")
              alphafix = true;
            if (terms[i].options[26] == "true")
              fcvarnonp[fcvarnonp.size()-1].set_stationary(alpha,alphafix);

            if (constlambda.getvalue() == false)
              {

              if(terms[i].options[22]=="true")
                fcvarnonp[fcvarnonp.size()-1].set_uniformprior();
              if(updatetau)
                fcvarnonp[fcvarnonp.size()-1].set_update_sigma2();
              fcvarnonp[fcvarnonp.size()-1].set_fcnumber(fullcond.size());
              fullcond.push_back(&fcvarnonp[fcvarnonp.size()-1]);
              }

            }

          if (terms[i].options[0] == "tpsplinerw1" || terms[i].options[0] == "tpsplinerw2")
            {
            fciwlspspline[fciwlspspline.size()-1].set_adaptiv();

            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_tvar.raw","_pspline_tvar.res","_pspline_tvariance");

            unsigned v = nu.getvalue();

            fctvariance.push_back(FULLCOND_tvariance(&generaloptions[generaloptions.size()-1],
            &fciwlspspline[fciwlspspline.size()-1],v,title,pathnonp,pathres)
                                );

            fctvariance[fctvariance.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fctvariance[fctvariance.size()-1]);

            } // end: if (terms[i].options[0] == "tpsplinerw1" ...

          if (terms[i].options[0] == "psplinerw1vrw1" || terms[i].options[0] == "psplinerw1vrw2" ||
              terms[i].options[0] == "psplinerw2vrw1" || terms[i].options[0] == "psplinerw2vrw2")
            {
            fciwlspspline[fciwlspspline.size()-1].set_adaptiv();

            make_paths(collinpred,pathnonp,pathres,title,terms[i].varnames[0],"",
                 "_pspline_vvar.raw","_pspline_vvar.res","_pspline_vvariance");

            f = (terms[i].options[10]).strtolong(h);
            unsigned minvar = unsigned(h);
            f = (terms[i].options[11]).strtolong(h);
            unsigned maxvar = unsigned(h);
            double startv;
            f = (terms[i].options[12]).strtodouble(startv);
            MCMC::fieldtype typevar;
            if (terms[i].options[0] == "psplinerw1vrw1" || terms[i].options[0] == "psplinerw2vrw1")
              typevar = MCMC::RW1;
            else
              typevar = MCMC::RW2;

            fcadaptiv.push_back(FULLCOND_adaptiv(&generaloptions[generaloptions.size()-1],
                                &fciwlspspline[fciwlspspline.size()-1],
                                typevar,title,a1,b1,ub,startv,minvar,maxvar,pathnonp,pathres));

            fcadaptiv[fcadaptiv.size()-1].set_fcnumber(fullcond.size());
            fullcond.push_back(&fcadaptiv[fcadaptiv.size()-1]);

            } // end: if (terms[i].options[0] == "psplinerw1vrw1" ...

          }

        */

        // ----------------- end:  non-gaussian response, etc. -----------------
        }

      }

    }

  return false;
  }


