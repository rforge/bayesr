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



#include "design_pspline.h"

using std::deque;

namespace MCMC
{


//------------------------------------------------------------------------------
//-------------- CLASS: DESIGN_mrf implementation of member functions ----------
//------------------------------------------------------------------------------

void DESIGN_pspline::read_options(vector<ST::string> & op,
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
  13      samplemult
  14      constraints
  15      round
  16      centermethod
  17      internal_mult
  18      pvalue
  19      meaneffect
  20      binning
  21      update
  22      nu
  23      maxdist
  24      ccovariate
  */


  int f;

  f = op[1].strtolong(degree);
  f = op[2].strtolong(nrknots);

  f = op[3].strtolong(difforder);
  if (difforder == 1)
    type = Rw1;
 else if (difforder==2)
   type = Rw2;
 else
   type = Rw3;

  if (op[7] == "false")   //nocenter==false, i.e. center
    center = true;
  else
    center = false;

  f = op[15].strtodouble(round);

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


  if (op[12] == "true")
    multeffect=true;
  else
    multeffect=false;

  f = op[20].strtodouble(binning);

  if (op[24] == "false")
    ccov = false;
  else
    ccov = true;

  if (op[26] == "false")
    derivative = false;
  else
    derivative = true;


  datanames = vn;

  }


DESIGN_pspline::DESIGN_pspline(void) : DESIGN()
  {

  }

  // CONSTRUCTOR

DESIGN_pspline::DESIGN_pspline(datamatrix & dm,datamatrix & iv,
                       GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl,
                       vector<ST::string> & op,
                       vector<ST::string> & vn)
                      : DESIGN(o,dp,fcl)
  {

  read_options(op,vn);

  datamatrix dmr = dm;
  if (round != -1)
    {
    unsigned i;
    for(i=0;i<dmr.rows();i++)
      dmr(i,0) =  floor(dm(i,0) * pow( 10, round) + 0.5) * pow(10, -round);

    // TEST
    // ofstream out("c:\\bayesx\\test\\results\\dmr.res");
    // dmr.prettyPrint(out);
    // TEST
    }

  if (binning != -1)
    {

    unsigned i;

    double a = dmr.min(0);
    double b = dmr.max(0);
    double delta = (b-a)/binning;
    double u1 = a+delta/2;
    double h;
    double dmrw;

    for(i=0;i<dmr.rows();i++)
      {
      dmrw = dmr(i,0);
      h = floor((dmrw- a)/delta);
      if (h >= binning)
      h-=1;
      dmr(i,0) = u1+h*delta;
      }

    }

  // TEST
  //  ofstream out("c:\\bayesx\\testh\\results\\dmr.res");
  //  dmr.prettyPrint(out);
  // TEST

  // centering dm
  if (ccov == true)
    {
    double dmr_mean = dmr.mean(0);
    double * workd = dmr.getV();
    unsigned i;
    for (i=0;i<dmr.rows();i++,workd++)
      *workd = *workd-dmr_mean;
    }
  // end centering

  init_data(dmr,iv);

  nrpar = nrknots-1+degree;

  make_Bspline();

  compute_Zout_transposed();

  weightK = vector<double>(nrpar,1);
  compute_penalty();

  XWX = envmatdouble(bandmatdouble(nrpar,degree,0));
  Wsum =datamatrix(posbeg.size(),1,1);
  XWres = datamatrix(nrpar,1);

  compute_precision(1.0);

  compute_basisNull();

  }


  // COPY CONSTRUCTOR

DESIGN_pspline::DESIGN_pspline(const DESIGN_pspline & m)
    : DESIGN(DESIGN(m))
  {
  multeffect=m.multeffect;
  knot = m.knot;
  nrknots = m.nrknots;
  degree = m.degree;
  difforder = m.difforder;
  weightK = m.weightK;
  round = m.round;
  binning = m.binning;
  ccov = m.ccov;
  minBS = m.minBS;
  maxBS = m.maxBS;
  }


  // OVERLOADED ASSIGNMENT OPERATOR

const DESIGN_pspline & DESIGN_pspline::operator=(const DESIGN_pspline & m)
  {
  if (this == &m)
    return *this;
  DESIGN::operator=(DESIGN(m));
  multeffect=m.multeffect;
  knot = m.knot;
  nrknots = m.nrknots;
  degree = m.degree;
  difforder = m.difforder;
  weightK = m.weightK;
  round = m.round;
  binning = m.binning;
  ccov = m.ccov;
  minBS = m.minBS;
  maxBS = m.maxBS;
  return *this;
  }

void DESIGN_pspline::outbasis_R(ofstream & out)
  {
  unsigned i;
  out << "BayesX.design.matrix<-function(x, ...) {" << endl;
  out << "  require(\"splines\")" << endl;
  out << "  x <- unlist(x)" << endl;
  out << "  knots <- c(";
  for(i = 0; i < knot.size() - 1; i++)
    out << knot[i] << ", ";
  out << knot[knot.size() - 1] << ")" << endl;
  out << "  degree <- " << degree << endl;
  out << "  ll <- knots[degree + 1]" << endl;
  out << "  ul <- knots[length(knots) - degree]" << endl;
  out << "  degree <- degree + 1" << endl;
  out << "  n <- length(x)" << endl;
  out << "  ind <- x <= ul & x >= ll" << endl;
  out << "  if(sum(ind) == n) {" << endl;
  out << "    X <- spline.des(knots, x, degree)$design" << endl;
  out << "  } else {" << endl;
  out << "    D <- spline.des(knots, c(ll, ll, ul, ul), degree, c(0, 1, 0, 1))$design" << endl;
  out << "    X <- matrix(0, n, ncol(D))" << endl;
  out << "    X[ind, ] <- spline.des(knots, x[ind], degree)$design" << endl;
  out << "    ind <- x < ll" << endl;
  out << "    if(sum(ind) > 0) " << endl;
  out << "    X[ind, ] <- cbind(1, x[ind] - ll) %*% D[1:2, ]" << endl;
  out << "    ind <- x > ul" << endl;
  out << "    if(sum(ind) > 0)" << endl;
  out << "    X[ind, ] <- cbind(1, x[ind] - ul) %*% D[3:4, ]" << endl;
  out << "  }" << endl;
  out << "  attr(X, \"type\") <- \"ps\"" << endl;
  out << "  return(X)" << endl;
  out << "}" << endl;
  }

void DESIGN_pspline::make_Bspline(void)
  {

  unsigned i,j,k;
  double value;

  datamatrix help;

// berechne x_min, x_max

  minBS = data(0,0);


  maxBS= data(data.rows()-1,0);

  double dist = maxBS-minBS;

  minBS -= 0.01*dist;
  maxBS += 0.01*dist;


// Knoten berechnen

  dist = (maxBS - minBS)/(nrknots-1);
  knot.push_back(minBS - degree*dist);
  for(i=1;i<nrknots+2*degree;i++)
    knot.push_back(knot[i-1] + dist);


  Zout = datamatrix(posbeg.size(),degree+1,0.0);
  index_Zout = statmatrix<int>(posbeg.size(),degree+1,0.0);
  double * work = Zout.getV();
  int * work_index = index_Zout.getV();

  help = datamatrix(nrpar,1,0.0);



  for (i=0;i<posbeg.size();i++)
    {
    value = data(posbeg[i],0);
    j=0;
    while(knot[degree+j+1] <= value)
      j++;
    help.assign(bspline(value));
    for (k=0;k<Zout.cols();k++,work++,work_index++)
      {
      *work = help(k+j,0);
      *work_index = j+k;
      }

    }

  consecutive = 1;


  if (derivative)
    make_Bspline_derivative();

  // TEST
  /*
  bool t = check_Zout_consecutive();

  ofstream out("c:\\bayesx\\test\\results\\Zout.res");
  Zout.prettyPrint(out);

  ofstream out2("c:\\bayesx\\test\\results\\index_Zout.res");
  index_Zout.prettyPrint(out2);

  datamatrix Zoutm(data.rows(),nrpar,0);
  for (i=0;i<posbeg.size();i++)
    {
    for (j=posbeg[i];j<=posend[i];j++)
      {
      for(k=0;k<Zout.cols();k++)
        Zoutm(j,index_Zout(i,k)) = sqrt(likep->workingweight(index_data(j,0),0))*Zout(i,k)*intvar(j,0);
      }

    }

  datamatrix h(data.rows(),1,1);
  datamatrix Zteins =  Zoutm.transposed()*h;

  datamatrix ZtZ = Zoutm.transposed()*Zoutm;

  ofstream out3("c:\\bayesx\\test\\results\\Zteins.res");
  Zteins.prettyPrint(out3);

  ofstream out4("c:\\bayesx\\test\\results\\ZtZ.res");
  ZtZ.prettyPrint(out4);
  */
  // TEST


  }


datamatrix DESIGN_pspline::bspline(const double & x)
  {
// nach Hämmerlin/Hoffmann
  datamatrix b(nrpar,1,0.0);
  datamatrix help(nrpar+degree,1,0.0);
  unsigned j;
  double * bwork;
  double * helpwork;

// Grad 0

  for(j=0;j<nrpar;j++)
    if( knot[j]<=x && x<knot[j+1])
      b(j,0) = 1.0;

  for(int l=1;l<=degree;l++)
    {
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
//      help(j,0) = b(j,0);
      *helpwork = *bwork;
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
      {
//      b(j,0) = (x-knot[j])*help(j,0)/(knot[j+l]-knot[j])
//                  + (knot[j+l+1]-x)*help(j+1,0)/(knot[j+l+1]-knot[j+1]);
      *bwork = (x-knot[j])**helpwork/(knot[j+l]-knot[j])
                  + (knot[j+l+1]-x)**(helpwork+1)/(knot[j+l+1]-knot[j+1]);

      }
    }

  return b;

  }


void DESIGN_pspline::make_Bspline_derivative(void)
  {

  unsigned i,j,k;
  double value;

  datamatrix help;


  Zout_derivative = datamatrix(posbeg.size(),degree+1,0.0);
  double * work = Zout_derivative.getV();

  help = datamatrix(nrpar,1,0.0);

  for (i=0;i<posbeg.size();i++)
    {
    value = data(posbeg[i],0);
    j=0;
    while(knot[degree+j+1] <= value)
      j++;
    help.assign(bspline_derivative(value));

    for (k=0;k<Zout_derivative.cols();k++,work++)
      {
      *work = help(k+j,0);
      }

    }


  // TEST

  // ofstream out("c:\\bayesx\\testh\\results\\Zout.res");
  // Zout_derivative.prettyPrint(out);

  /*
  datamatrix Zoutm(data.rows(),nrpar,0);
  for (i=0;i<posbeg.size();i++)
    {
    for (j=posbeg[i];j<=posend[i];j++)
      {
      for(k=0;k<Zout.cols();k++)
        Zoutm(j,index_Zout(i,k)) = sqrt(likep->workingweight(index_data(j,0),0))*Zout(i,k)*intvar(j,0);
      }

    }

  datamatrix h(data.rows(),1,1);
  datamatrix Zteins =  Zoutm.transposed()*h;

  datamatrix ZtZ = Zoutm.transposed()*Zoutm;

  ofstream out3("c:\\bayesx\\test\\results\\Zteins.res");
  Zteins.prettyPrint(out3);

  ofstream out4("c:\\bayesx\\test\\results\\ZtZ.res");
  ZtZ.prettyPrint(out4);
  */
  // TEST


  }



datamatrix DESIGN_pspline::bspline_derivative(const double & x)
  {


// nach Hämmerlin/Hoffmann
  datamatrix b(nrpar,1,0.0);
  datamatrix help(nrpar+degree,1,0.0);
  unsigned j;
  double * bwork;
  double * helpwork;

// Grad 0

  for(j=0;j<nrpar;j++)
    if( knot[j]<=x && x<knot[j+1])
      b(j,0) = 1.0;

  for(int l=1;l<=degree-1;l++)
    {
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
//      help(j,0) = b(j,0);
      *helpwork = *bwork;
    bwork = b.getV();
    helpwork = help.getV();
    for(j=0;j<nrpar;j++,++helpwork,++bwork)
      {
//      b(j,0) = (x-knot[j])*help(j,0)/(knot[j+l]-knot[j])
//                  + (knot[j+l+1]-x)*help(j+1,0)/(knot[j+l+1]-knot[j+1]);
      *bwork = (x-knot[j])**helpwork/(knot[j+l]-knot[j])
                  + (knot[j+l+1]-x)**(helpwork+1)/(knot[j+l+1]-knot[j+1]);

      }
    }

// Hämmerlin/Hoffmann Seite 263

  bwork = b.getV();
  helpwork = help.getV();
  for(j=0;j<nrpar;j++,++helpwork,++bwork)
    *helpwork = *bwork;
  bwork = b.getV();
  helpwork = help.getV();
  for(j=0;j<nrpar;j++,++helpwork,++bwork)
    *bwork = degree*( *helpwork/(knot[j+degree]-knot[j]) - *(helpwork+1)/(knot[j+degree+1]-knot[j+1]) );

  return b;

  }


void DESIGN_pspline::compute_penalty2(const datamatrix & pen)
  {
  if (type==Rw1)
    {
    K = Krw1env(pen);
    rankK = nrpar-1;
    }
  else if (type==Rw2)
    {
    K = Krw2env(nrpar);
    rankK = nrpar-2;
    }
  else if (type==Rw3)
    {
    K = Krw3env(nrpar);
    rankK = nrpar-3;
    }
  }


void DESIGN_pspline::compute_penalty(void)
  {
  if (type==Rw1)
    {
    K = Krw1env(weightK);
    rankK = nrpar-1;
    }
  else if (type==Rw2)
    {
    K = Krw2env(nrpar);
    rankK = nrpar-2;
    }
  else if (type==Rw3)
    {
    K = Krw3env(nrpar);
    rankK = nrpar-3;
    }
  }


void DESIGN_pspline::compute_betaweight(datamatrix & betaweight)
  {

  unsigned i;

  if(degree==1)
    {
    betaweight = datamatrix(1,nrpar,1);

    betaweight(0,0) = 0.5;
    betaweight(0,nrpar-1) = 0.5;

    for(i=0;i<betaweight.cols();i++)
      betaweight(0,i) /= (nrknots-1);
    }
  else if(degree==2)
    {
    betaweight = datamatrix(1,nrpar,1);

    betaweight(0,0) = 1/6.0;
    betaweight(0,nrpar-1) = 1/6.0;
    betaweight(0,1) = 5/6.0;
    betaweight(0,nrpar-2) = 5/6.0;

    for(i=0;i<betaweight.cols();i++)
      betaweight(0,i) /= (nrknots-1);
    }
  else if(degree==3)
    {
    betaweight = datamatrix(1,nrpar,1);

    betaweight(0,0) = 1/24.0;
    betaweight(0,nrpar-1) = 1/24.0;
    betaweight(0,1) = 12/24.0;
    betaweight(0,nrpar-2) = 12/24.0;
    betaweight(0,2) = 23/24.0;
    betaweight(0,nrpar-3) = 23/24.0;

    for(i=0;i<betaweight.cols();i++)
      betaweight(0,i) /= (nrknots-1);
    }
  else
    {
    betaweight = datamatrix(1,nrpar,1.0/double(nrpar));
    }

  }


void DESIGN_pspline::compute_basisNull(void)
  {

  if (center==true)
  {
  unsigned i,j;

  if (centermethod==meancoeff || centermethod==meansimple)
    {
    basisNull = datamatrix(1,nrpar,1);
    position_lin = -1;
    }
  else if (centermethod==cmeaninvvar)
    {
    compute_precision(10);

    envmatdouble precisioninv;
    precisioninv = envmatdouble(0.0,nrpar,degree>2?degree:2);
    precision.inverse_envelope(precisioninv);

    basisNull = datamatrix(1,nrpar,1);

    unsigned k;
    for (k=0;k<nrpar;k++)
      basisNull(0,k) = 1/precisioninv.getDiag(k);

    position_lin = -1;
    }
  else if ((centermethod==cmeanintegral) || (centermethod==integralsimple))  // integral f = 0
    {
    compute_betaweight(basisNull);
    position_lin = -1;
    }
  else if (centermethod==meanf)            // sum of f's zero (over all observations)
    {

    basisNull = datamatrix(1,nrpar,1);

    unsigned k;
    for (k=0;k<nrpar;k++)
      basisNull(0,k) = compute_sumBk(k);

    // TEST
    // ofstream out("c:\\bayesx\\testh\\results\\basisnull.res");
    // basisNull.prettyPrint(out);
    // ende: TEST

    position_lin = -1;

    }
  else if (centermethod==meanfd)            // sum of f's zero
    {

    basisNull = datamatrix(1,nrpar,1);

    unsigned k;
    for (k=0;k<nrpar;k++)
      basisNull(0,k) = compute_sumBk_different(k);

    // TEST
    // ofstream out("c:\\bayesx\\testh\\results\\basisnull.res");
    // basisNull.prettyPrint(out);
    // ende: TEST

    position_lin = -1;

    }
  else if (centermethod==nullspace)
    {

    if (type==Rw1)
      {
      basisNull = datamatrix(1,nrpar,1);
      position_lin = -1;
      }
    else if (type==Rw2)
      {
      basisNull = datamatrix(2,nrpar,1);
      deque<double>::iterator it = knot.begin();
      for (i=0;i<nrpar;i++,++it)
        basisNull(1,i) = *it;
      }

    }


  for(i=0;i<basisNull.rows();i++)
    {
    basisNullt.push_back(datamatrix(basisNull.cols(),1));
    for(j=0;j<basisNull.cols();j++)
      basisNullt[i](j,0) = basisNull(i,j);
    }


  if (basisNull.rows() > 1)
    {
    designlinear = datamatrix(posbeg.size(),basisNull.rows()-1);

    double * workdl = designlinear.getV();
    double h;
    for(i=0;i<posbeg.size();i++)
      for(j=0;j<designlinear.cols();j++,workdl++)
        {
        h = data(posbeg[i],0);
        *workdl =  pow(static_cast<double>(h),static_cast<double>(j+1));
        }
    }

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




void DESIGN_pspline::compute_precision(double l)
  {
  if (precisiondeclared==false)
    {
    precision = envmatdouble(0.0,nrpar,degree>2?degree:2);
    precisiondeclared = true;
    }

  precision.addto(XWX,K,1.0,l);

  /*
  // TEST
  ofstream out2("c:\\bayesx\\test\\results\\XWX.res");
  XWX.print2(out2);

  ofstream out3("c:\\bayesx\\test\\results\\K.res");
  K.print2(out3);

  ofstream out("c:\\bayesx\\testh\\results\\precision.res");
  precision.print2(out);

  ofstream out2("c:\\bayesx\\testh\\results\\env.res");
  vector<double> e = precision.getEnv();
  unsigned i;
  for (i=0;i!= e.size();i++)
    out2 << e[i] << endl;

  ofstream out3("c:\\bayesx\\testh\\results\\xenv.res");
  vector<unsigned> xe = precision.getXenv();
  for (i=0;i!= xe.size();i++)
    out3 << xe[i] << endl;
  */

  // TEST

  }

void DESIGN_pspline::outoptions(GENERAL_OPTIONS * op)
  {

  ST::string typestr;

 if (type == Rw1)
    typestr = "first order random walk";
  else if (type == Rw2)
    typestr = "second order random walk";
  else if (type == Rw3)
    typestr = "third order random walk";

  ST::string centerm;

  if (center==false)
    centerm = "uncentered sampling";
  else
    centerm = "centered sampling";

  double min;
  min = data(0,0);

  double max;
  max = data(data.rows()-1,0);

  double dist = max-min;

  min -= 0.01*dist;
  max += 0.01*dist;

  dist = (max - min)/(nrknots-1);

  op->out("  Prior: " + typestr + "\n");
  op->out("  Number of knots: " + ST::inttostring(nrknots) + "\n" );
  op->out("  First knot: " + ST::doubletostring(min,8) + "\n");
  op->out("  Last knot: " + ST::doubletostring(max,8) + "\n");
  op->out("  Knot distance: " + ST::doubletostring(dist,8) + "\n");
  op->out("  Degree of Splines: " + ST::inttostring(degree) + "\n" );
  op->out("  " + centerm + "\n" );


  op->out("  B-spline basis may be created e.g. in STATA by: \n");
  op->out("  bspline , xvar(x) generate(Bs) power(" + ST::inttostring(degree) +
         ") knots(" + ST::doubletostring(min,8) + "(" +
          ST::doubletostring(dist,8) + ")" + ST::doubletostring(max,8) + ") \n"
         );
  op->out("  where x contains the values for which the basis functions should be created\n");
  op->out("\n");

  }

} // end: namespace MCMC



