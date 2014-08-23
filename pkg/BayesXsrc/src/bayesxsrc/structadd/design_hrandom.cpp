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



#include "design_hrandom.h"
#include "clstring.h"

namespace MCMC
{


//------------------------------------------------------------------------------
//------------ CLASS: DESIGN_hrandom implementation of member functions --------
//------------------------------------------------------------------------------


void DESIGN_hrandom::read_options(vector<ST::string> & op,
vector<ST::string> & vn)
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
  */

  datanames = vn;

  if (op[16]=="meancoeff")
    centermethod = meancoeff;
  else if (op[16] == "meansimple")
    centermethod = meansimple;
  else if (op[16] == "integralsimple")
    centermethod = integralsimple;
  else if (op[16] == "nullspace")
    centermethod = nullspace;
  else if (op[16] == "meaninvvar")
    centermethod = cmeaninvvar;
  else if (op[16] == "meanintegral")
    centermethod = cmeanintegral;
  else if (op[16] == "meanf")
    centermethod = meanf;
  else if (op[16] == "meanfd")
    centermethod = meanfd;
  else if (op[16] == "meansum2")
    centermethod = meansum2;


  if (op[46] == "true")
    center = true;
  else
    center = false;

  }


DESIGN_hrandom::DESIGN_hrandom(void)
  {

  }





DESIGN_hrandom::DESIGN_hrandom(const datamatrix & dm, const datamatrix & iv,
                               GENERAL_OPTIONS * o,DISTR * dp, FC_linear * fcl,
                               DISTR * dp_RE,
                               vector<ST::string> & op,vector<ST::string> & vn)
      : DESIGN(o,dp,fcl)
  {

  read_options(op,vn);

  likep_RE = dp_RE;

  discrete = true;

  type = Hrandom;

  init_data(dm,iv);
  nrpar = posbeg.size();

  Zout = datamatrix(posbeg.size(),1,1);
  index_Zout = statmatrix<int>(Zout.rows(),1);
  index_Zout.indexinit();

  consecutive = 1;

  compute_Zout_transposed();

  compute_penalty();

  XWX = envmatdouble(0,nrpar);
  XWres = datamatrix(nrpar,1);
  Wsum = datamatrix(nrpar,1,1);

  compute_precision(1.0);

  centermethod=meansimplevar;

  }


DESIGN_hrandom::DESIGN_hrandom(const DESIGN_hrandom & m)
  : DESIGN(DESIGN(m))
  {
  likep_RE = m.likep_RE;
  }


const DESIGN_hrandom & DESIGN_hrandom::operator=(const DESIGN_hrandom & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));
  likep_RE = m.likep_RE;
  return *this;

  }


void DESIGN_hrandom::init_data(const datamatrix & dm, const datamatrix & iv)
  {


  DESIGN::init_data(dm,iv);

  meaneffectnr = compute_modecategorie();

  }


void DESIGN_hrandom::outbasis_R(ofstream & out)
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
  out << "  attr(X, \"type\") <- \"re\"" << endl;
  out << "  return(X)" << endl;
  out << "}" << endl;
  }


void DESIGN_hrandom::compute_penalty(void)
  {
  K =   envmatrix<double>(1,nrpar);
  rankK = nrpar;
  }


void DESIGN_hrandom::compute_penalty2(const datamatrix & pen)
  {
  if (K.getDim() != pen.rows())
    {
    K =   envmatrix<double>(1,nrpar);
    rankK = nrpar;
    }
  else
    {
    unsigned i;
    double * workpen = pen.getV();
    for (i=0;i<pen.rows();i++,workpen++)
      {
      K.setDiag(i,1/(*workpen));
      }
    }
  }



void DESIGN_hrandom::compute_XtransposedWres(datamatrix & partres, double l)
  {

  double * workXWres = XWres.getV();

  double * linpredREp;
  if (likep_RE->linpred_current==1)
    linpredREp = likep_RE->linearpred1.getV();
  else
    linpredREp = likep_RE->linearpred2.getV();

  double * partresp = partres.getV();

  unsigned i;

  for(i=0;i<nrpar;i++,workXWres++,linpredREp++,partresp++)
    *workXWres =  l*(*linpredREp)+(*partresp);

  XWres_p = &XWres;

  }



void DESIGN_hrandom::compute_basisNull(void)
  {

  if (center==true)
    {
    unsigned i,j;

    if (centermethod==meancoeff || centermethod==meansimple || centermethod==meanfd)
      {
      basisNull = datamatrix(1,nrpar,1);
      position_lin = -1;
      }

    for(i=0;i<basisNull.rows();i++)
      {
      basisNullt.push_back(datamatrix(basisNull.cols(),1));
      for(j=0;j<basisNull.cols();j++)
        basisNullt[i](j,0) = basisNull(i,j);
      }

    }

  }





void DESIGN_hrandom::compute_precision(double l)
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



void DESIGN_hrandom::compute_meaneffect(DISTR * level1_likep,double & meaneffect,
                                datamatrix & beta,datamatrix & meaneffectbeta,
                                bool computemeaneffect,double meaneffectconstant)

  {

  level1_likep->meaneffect -= meaneffect;

  double * linpredREp;
  double linm;
  if (likep_RE->linpred_current==1)
    {
    linpredREp = likep_RE->linearpred1.getV();
    linm = likep_RE->linearpred1(meaneffectnr,0);
    }
  else
    {
    linpredREp = likep_RE->linearpred2.getV();
    linm = likep_RE->linearpred2(meaneffectnr,0);
    }

  meaneffect = meaneffectintvar*(beta(meaneffectnr,0) - linm);

  if (computemeaneffect==true)
    {
    double co;

    if (meaneffectconstant == 0)
      co = level1_likep->meaneffect;
    else
      co = meaneffectconstant;

    unsigned i;
    double * betap = beta.getV();
    double * meffectp = meaneffectbeta.getV();
    double l;
    for(i=0;i<beta.rows();i++,meffectp++,betap++,linpredREp++)
      {
      l=co+meaneffectintvar*((*betap)- (*linpredREp));
      level1_likep->compute_mu(&l,meffectp);
      }
    }

  level1_likep->meaneffect += meaneffect;

  }




void DESIGN_hrandom::outoptions(GENERAL_OPTIONS * op)
  {

  }


} // end: namespace MCMC



