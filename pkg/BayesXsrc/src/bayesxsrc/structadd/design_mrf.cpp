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



#include "design_mrf.h"
#include "clstring.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//-------------- CLASS: DESIGN_mrf implementation of member functions ----------
//------------------------------------------------------------------------------


void DESIGN_mrf::read_options(vector<ST::string> & op,vector<ST::string> & vn)
  {

  /*
  1       degree
  2       numberknots
  3       difforder
  4       lambda
  5       a
  6       b
  7       center
  8       map
  9       lambda_re
  10      a_re
  11      b_re
  12      internal_mult
  13      samplemult
  14      constraints
  15      round
  16      centermethod
  */


  if (op[7] == "false")   //nocenter==false, i.e. center
    center = true;
  else
    center = false;

  if (op[16]=="meancoeff" || op[16] == "nullspace")
    centermethod = meancoeff;
  else if (op[16] == "meansimple")
    centermethod = meansimple;
  else if (op[16] == "meaninvvar")
    centermethod = cmeaninvvar;
  else if (op[16] == "meanintegral")
    centermethod = cmeanintegral;
  else if (op[16] == "meanf")
    centermethod = meanf;

  datanames = vn;

  }


DESIGN_mrf::DESIGN_mrf(void) : DESIGN()
  {

  }

  // CONSTRUCTOR 1
  // Spatial covariates

DESIGN_mrf::DESIGN_mrf(const datamatrix & dm,const datamatrix & iv,
                       GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl,
                       const MAP::map & m,vector<ST::string> & op,
                       vector<ST::string> & vn)
                      : DESIGN(o,dp,fcl)
  {

  if (errors==false)
    {

    read_options(op,vn);

    discrete = true;

    ma = m;
    type = Mrf;

    nrpar = ma.get_nrregions();
    consecutive=true;

    Zout = datamatrix(nrpar,1,1);
    index_Zout = statmatrix<int>(Zout.rows(),1);
    index_Zout.indexinit();
    }

  if (errors==false)
    init_data(dm,iv);

  if (errors==false)
    {
    compute_penalty();

    XWX = envmatdouble(0,nrpar);
    XWres = datamatrix(nrpar,1);
    Wsum = datamatrix(nrpar,1,1);

    compute_precision(1.0);

    compute_basisNull();

    identity=true;
    }
  }


  // COPY CONSTRUCTOR

DESIGN_mrf::DESIGN_mrf(const DESIGN_mrf & m)
    : DESIGN(DESIGN(m))
  {
  ma = m.ma;
  }

  // OVERLOADED ASSIGNMENT OPERATOR

const DESIGN_mrf & DESIGN_mrf::operator=(const DESIGN_mrf & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));
  ma = m.ma;
  return *this;
  }


void DESIGN_mrf::init_data(const datamatrix & dm, const datamatrix & iv)
  {

  if (ma.get_bandsize() > 40)
    ma.reorderopt();

  ma.compute_reg(dm,posbeg,posend,effectvalues,index_data);

  vector<ST::string> errormessages = ma.get_errormessages();
  if (errormessages.size() >= 1)
    {

    errors=true;
    unsigned i;
    for (i=0;i<errormessages.size();i++)
      optionsp->out(errormessages[i]);
    }

  if (errors==false)
    {
    vector<int> posbeg_help = posbeg;
    vector<int> posend_help = posend;
    vector<ST::string> effectvalues_help = effectvalues;

    DESIGN::init_data(dm,iv);

    posbeg = posbeg_help;
    posend = posend_help;
    effectvalues = effectvalues_help;

    int k,j;
    int * workindex = index_data.getV();
    ind = statmatrix<unsigned>(dm.rows(),1);
    for (j=0;j<posend.size();j++)
      {

      if (posbeg[j] != -1)
        {
        for (k=posbeg[j];k<=posend[j];k++,workindex++)
          ind(*workindex,0) = j;
        }
      else
        {
        optionsp->out("NOTE: no observations for region " + effectvalues[j] + "\n");
        }

      }


    // TEST
    /*
    int j;
    ofstream out2("c:\\bayesx\\testh\\results\\posend_neu.res");
    for (j=0;j<posend.size();j++)
      out2 << posend[j] << endl;

    ofstream out2a("c:\\bayesx\\testh\\results\\posbeg_neu.res");
    for (j=0;j<posbeg.size();j++)
      out2a << posbeg[j] << endl;

    ofstream out3("c:\\bayesx\\testh\\results\\effectvalues_neu.res");
    for (j=0;j<effectvalues.size();j++)
      out3 << effectvalues[j] << endl;

    ofstream out4("c:\\bayesx\\testh\\results\\index_data_neu.res");
    index_data.prettyPrint(out4);
    */
    // TEST


    meaneffectnr = compute_modecategorie();

    if (ma.get_errormessages().size() > 0)
      {
      //  FEHLT!!
      }

    } // end: if (errors==false)

  }


void DESIGN_mrf::outbasis_R(ofstream & out)
  {
  unsigned i;
  out << "BayesX.design.matrix <- function(x, ...) {" << endl;
  out << "  w <- getOption(\"warn\")" << endl;
  out << "  options(warn = -1)" << endl;
  out << "  tx <- as.integer(as.character(unlist(x)))" << endl;
  out << "  options(\"warn\" = w)" << endl;
  out << "  x <- if(!any(is.na(tx))) tx else as.integer(x)" << endl;
  out << "  levels <- c(";
  for(i = 0; i < effectvalues.size() - 1; i++)
    out << effectvalues[i] << ", ";
  out << effectvalues[effectvalues.size() - 1] << ")" << endl;
  out << "  x <- factor(x, levels = levels)" << endl;
  out << "  X <- diag(length(levels))[x, ]" << endl;
  out << "  attr(X, \"type\") <- \"mrf\"" << endl;
  out << "  return(X)" << endl;
  out << "}" << endl;
  }


void DESIGN_mrf::compute_penalty(void)
  {
  if (type==Mrf)
    K = Kmrfenv(ma);
  rankK = ma.get_nrregions()-1;
  }



void DESIGN_mrf::compute_basisNull(void)
  {
  int i,j;

  basisNull = datamatrix(1,nrpar,1);

  if (centermethod==meanf || centermethod==cmeanintegral)
    {

    unsigned k;
    for (k=0;k<nrpar;k++)
      {
      if (posbeg[k] != -1)
        basisNull(0,k) = posend[k]-posbeg[k]+1;
      else
        basisNull(0,k) = 0;
      }

    }


  position_lin = -1;


  for(i=0;i<basisNull.rows();i++)
    {
    basisNullt.push_back(datamatrix(basisNull.cols(),1));
    for(j=0;j<basisNull.cols();j++)
      basisNullt[i](j,0) = basisNull(i,j);
    }


  // TEST
  /*
    ofstream out("c:\\bayesx\\test\\results\\data.res");
    data.prettyPrint(out);

    ofstream out2("c:\\bayesx\\test\\results\\designlin.res");
    designlinear.prettyPrint(out2);
   */
  // TEST

  }


void DESIGN_mrf::compute_XtransposedWres(datamatrix & partres, double l)
  {
  XWres_p = &partres;
  }


void DESIGN_mrf::compute_precision(double l)
  {

  if (precisiondeclared==false)
    {
    precision = envmatdouble(K.getXenv(),0,nrpar);
    precisiondeclared = true;
    }

  precision.addtodiag(XWX,K,1.0,l);

  /*
  // TEST
  ofstream out2("c:\\bayesx\\test\\results\\XWX.res");
  XWX.print2(out2);

  ofstream out3("c:\\bayesx\\test\\results\\K.res");
  K.print2(out3);


  ofstream out("c:\\bayesx\\test\\results\\precision.res");
  precision.print2(out);
  // TEST
  */

  }


void DESIGN_mrf::outoptions(GENERAL_OPTIONS * op)
  {

  }


/*
void DESIGN::compute_partres(datamatrix & res, datamatrix & f,bool cwsum)
  {

  double * workingresponsep = likep->workingresponse.getV();
  double * worklinp;
  if (likep->linpred_current==1)
    worklinp = likep->linearpred1.getV();
  else
    worklinp = likep->linearpred2.getV();

  double * workingweightp = likep->workingweight.getV();
  unsigned * indp = ind.getV();

  int i,j;
  double * resp = res.getV();
  for (i=0;i<res.rows();i++,resp++)
    *resp =  0;


  if (intvar.rows()==data.rows())   // varying coefficient
    {

    double * workintvar = intvar.getV();

    if ((likep->wtype==wweightsnochange_one) && (cwsum==false) && (changingdesign==false))
      {

      for (i=0;i<ind.rows();i++,workingresponsep++,indp++,worklinp++,workintvar++)
        res(*indp,0) += *workintvar * ( *workingresponsep - *worklinp +
                        (*workintvar) * f(*indp,0));

      }
    else
      {
      if ((likep->wtype==wweightsnochange_constant) && (cwsum==false) && (changingdesign==false))
        {

        for (i=0;i<ind.rows();i++,workingresponsep++,indp++,worklinp++,
                                workingweightp++,workintvar++)
          res(*indp,0) +=  (*workingweightp) *  (*workintvar) *
                           (*workingresponsep - *worklinp +
                           (*workintvar) * f(*indp,0));


        }
      else  // wweightschange
        {

        double * Wsump = Wsum.getV();
        for (i=0;i<Wsum.rows();i++,Wsump++)
          *Wsump =  0;

        double * workintvar2 = intvar2.getV();

        for (i=0;i<ind.rows();i++,workingresponsep++,indp++,worklinp++,
                                workingweightp++,workintvar++,workintvar2++)
          {
          res(*indp,0) +=  (*workingweightp) * (*workintvar) *
                            (*workingresponsep - *worklinp +
                            (*workintvar) * f(*indp,0));

          Wsum(*indp,0) += *workingweightp * (*workintvar2);
          }

        }

      }

    }
  else                              // additive
    {
    if ((likep->wtype==wweightsnochange_one) && (cwsum==false) && (changingdesign==false))
      {

      for (j=0;j<posend.size();j++)
        {
        if (posbeg[j] != -1)
          for (i=posbeg[j];i<=posend[j];i++,workingresponsep++,indp++,worklinp++)
            res(*indp,0) += *workingresponsep - *worklinp + f(*indp,0);
        }

      }
    else
      {
      if ((likep->wtype==wweightsnochange_constant) && (cwsum==false) && (changingdesign==false))
        {

        for (i=0;i<ind.rows();i++,workingresponsep++,indp++,worklinp++,
                                workingweightp++)
          res(*indp,0) +=  (*workingweightp) *
                          (*workingresponsep - *worklinp + f(*indp,0));


        }
      else  // wweightschange
        {

        double * Wsump = Wsum.getV();
        for (i=0;i<Wsum.rows();i++,Wsump++)
          *Wsump =  0;

        for (i=0;i<ind.rows();i++,workingresponsep++,indp++,worklinp++,
                                workingweightp++)
          {
          res(*indp,0) +=  (*workingweightp) *
                            (*workingresponsep - *worklinp + f(*indp,0));

          Wsum(*indp,0) += *workingweightp;
          }

        }

      }

    }

  // TEST

  ofstream out0("c:\\bayesx\\testh\\results\\residuum.res");
  res.prettyPrint(out0);

  // ofstream out("c:\\bayesx\\test\\results\\tildey.res");
  // (likep->workingresponse).prettyPrint(out);
  // TEST

  }
*/


} // end: namespace MCMC



