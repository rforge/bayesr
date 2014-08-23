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


#include "fullcond_surf_gaussian.h"

namespace MCMC
{

SparseMatrix FULLCOND_surf_gaussian::Krw1(const vector<double> & weight)
  {

  unsigned S = weight.size();
  datamatrix K(S,S,0);
  unsigned i;
  for (i=1;i<S-1;i++)
    {
    K(i,i) = 1.0/weight[i]+1.0/weight[i+1];
    K(i,i-1) = -1.0/(weight[i]);
    K(i,i+1) = -1.0/(weight[i+1]);
    }
  K(0,0) = 1.0/weight[1];
  K(0,1) = -1.0/weight[1];
  K(S-1,S-1) = 1.0/weight[S-1];
  K(S-1,S-2) = -1.0/weight[S-1];

  return SparseMatrix(K,true);


  }


SparseMatrix FULLCOND_surf_gaussian::Krw2(const vector<double> & weight)
  {
  unsigned i;
  int S = weight.size();

  datamatrix F(S-2,S,0);
  for (i=0;i<F.rows();i++)
    {
    F(i,i)   = weight[2+i]/weight[1+i];
    F(i,i+1) = -(1+weight[2+i]/weight[1+i]);
    F(i,i+2) = 1;
    }
  datamatrix Q(S-2,S-2,0);
  for(i=0;i<Q.rows();i++)
	 Q(i,i) =   weight[2+i]*(1+weight[2+i]/weight[1+i]);

  datamatrix K = F.transposed()*Q.inverse()*F;

  return SparseMatrix(K,true);
  }


SparseMatrix FULLCOND_surf_gaussian::Kmrflinear(const unsigned & nr1,
const unsigned & nr2)
  {

  SparseMatrix K(nr1*nr2,nr1*nr2,4);

  unsigned i,j;
  unsigned row = 0;
  for(i=0;i<nr1;i++)
    {

    for (j=0;j<nr2;j++)
      {

      if ( (i==0) && (j==0) )
        {
        K.put(0,0,2);
        K.put(0,1,-1);
        K.put(0,nr2,-1);
        }
      else if ( (i==0) && (j==nr2-1) )
        {
        K.put(row,row,2);
        K.put(row,row-1,-1);
        K.put(row,row+nr2,-1);
        }
      else if ( (i==nr1-1) && (j==0) )
        {
        K.put(row,row,2);
        K.put(row,row+1,-1);
        K.put(row,row-nr2,-1);
        }
      else if ( (i==nr1-1) && (j==nr2-1) )
        {
        K.put(row,row,2);
        K.put(row,row-1,-1);
        K.put(row,row-nr2,-1);
        }
      else if ( (i==0) && (j>0) && (j < nr2-1) )
        {
        K.put(row,row,3);
        K.put(row,row-1,-1);
        K.put(row,row+1,-1);
        K.put(row,row+nr2,-1);
        }
      else if ( (i==0) && (j>0) && (j < nr2-1) )
        {
        K.put(row,row,3);
        K.put(row,row-1,-1);
        K.put(row,row+1,-1);
        K.put(row,row+nr2,-1);
        }
      else if ( (i==nr1-1) && (j>0) && (j < nr2-1) )
        {
        K.put(row,row,3);
        K.put(row,row-1,-1);
        K.put(row,row+1,-1);
        K.put(row,row-nr2,-1);
        }
      else if ( (j==0) && (i>0) && (i < nr1-1) )
        {
        K.put(row,row,3);
        K.put(row,row+1,-1);
        K.put(row,row-nr2,-1);
        K.put(row,row+nr2,-1);
        }
      else if ( (j==nr2-1) && (i>0) && (i < nr1-1) )
        {
        K.put(row,row,3);
        K.put(row,row-1,-1);
        K.put(row,row-nr2,-1);
        K.put(row,row+nr2,-1);
        }
      else
        {
        K.put(row,row,4);
        K.put(row,row+1,-1);
        K.put(row,row-1,-1);
        K.put(row,row-nr2,-1);
        K.put(row,row+nr2,-1);
        }


      row++;
      }

    }

  return K;

  }


void FULLCOND_surf_gaussian::make_moddata(const datamatrix & moddata1,
                                          const datamatrix & moddata2)
  {

  unsigned i,j;

  datamatrix moddata = datamatrix(moddata1.rows(),1);
  for(i=0;i<moddata.rows();i++)
    moddata(i,0) = (moddata1(i,0)-1) * sizeK2 + moddata2(i,0);

  index = statmatrix<int>(moddata.rows(),1);
  index.indexinit();
  moddata.indexsort(index,0,moddata.rows()-1,0,0);

  posbeg = vector<int>(sizeK,-1);
  posend = vector<int>(sizeK,-1);

  posbeg[unsigned(moddata(index(0,0),0))-1] = 0;
  for(j=1;j<moddata.rows();j++)
    {
    if (moddata(index(j,0),0) != moddata(index(j-1,0),0))
      {
      posbeg[unsigned(moddata(index(j,0),0)) -1] = j;
      posend[unsigned(moddata(index(j-1,0),0))-1] = j-1;
      }

    }

  posend[unsigned(moddata(index(moddata.rows()-1,0),0))-1] = moddata.rows()-1;

  }


datamatrix FULLCOND_surf_gaussian::make_categories(const datamatrix & moddata,
                                                 vector<ST::string> & effvalues)
  {

  posbeg.erase(posbeg.begin(),posbeg.end());
  posend.erase(posend.begin(),posend.end());
  weight.erase(weight.begin(),weight.end());
  effvalues.erase(effvalues.begin(),effvalues.end());

  unsigned j;                              // Loopingvariable

  // Initialization of the index, Indexsort of moddata

  index = statmatrix<int>(moddata.rows(),1);
  index.indexinit();
  moddata.indexsort(index,0,moddata.rows()-1,0,0);

  double diff;                             // difference of two succeeding
                                           // observations

  unsigned nr=1;

  datamatrix result(moddata.rows(),1);

  posbeg.push_back(0);
  weight.push_back(0);
  effvalues.push_back(ST::doubletostring(moddata(index(0,0),0),6));

  result(index(0,0),0) = 1;

  for(j=1;j<moddata.rows();j++)
    {

    diff = moddata(index(j,0),0) - moddata(index(j-1,0),0);

    if (diff > 0)
      {
      weight.push_back(diff);
//      weight.push_back(1.0);  
      effvalues.push_back(ST::doubletostring(moddata(index(j,0),0),6));
      posbeg.push_back(j);
      posend.push_back(j-1);
      nr++;
      }

    result(index(j,0),0) = nr;

    }  // end: for(j=1;j<moddata.rows();j++)

  posend.push_back(moddata.rows()-1);

  return result;
  }  // end: function make_categories


FULLCOND_surf_gaussian::FULLCOND_surf_gaussian(MCMCoptions * o,DISTRIBUTION * dp,
                          const datamatrix & d1, const datamatrix & d2,
                          const double & a, const double & b,
                          const fieldtype & ft, const ST::string & ti,
                          const ST::string & fp, const ST::string & pres,
                          const unsigned & c,
                          bool sb
                        )
                       : FULLCOND_nonp_basis(o,dp,ft,ti,fp,pres,c)

  {

  unsigned i,j,k;

  singleblock = sb;


  centertotal = true;

  datamatrix d1kat;
  datamatrix d2kat;

  vector<ST::string> effvalue1;
  vector<ST::string> effvalue2;
  SparseMatrix K1,K2;


  unsigned bands;
  unsigned bandssingle;

  if (type == mrfkr1)
    {
    d1kat = make_categories(d1,effvalue1);
    K1 = Krw1(weight);


    d2kat = make_categories(d2,effvalue2);
    K2 = Krw1(weight);

    Ksp = K1.kronecker(K2);
    sizeK1 = K1.get_rows();
    sizeK2 = K2.get_rows();

    rankK = (sizeK1-1)*(sizeK2-1);

    bands=1;
    bandssingle = sizeK2+1;

    }
  else if (type == mrfkr2)
    {

    d1kat = make_categories(d1,effvalue1);
    K1 = Krw2(weight);

    d2kat = make_categories(d2,effvalue2);
    K2 = Krw2(weight);

    Ksp = K1.kronecker(K2);

    sizeK1 = K1.get_rows();
    sizeK2 = K2.get_rows();

    rankK = (sizeK1-2)*(sizeK2-2);

    bands=2;
    bandssingle = 2*(sizeK2+1);

    }
  else if (type==mrflinear)
    {

    d1kat = make_categories(d1,effvalue1);

    d2kat = make_categories(d2,effvalue2);

    sizeK1 = effvalue1.size();
    sizeK2 = effvalue2.size();

    Ksp = Kmrflinear(sizeK1,sizeK2);

    rankK = sizeK1*sizeK2-1;

    bands=1;
    bandssingle = sizeK2;

    }


  sizeK = sizeK1*sizeK2;

  sumx1 = datamatrix(sizeK1,1,0);
  sumx2 = datamatrix(sizeK2,1,0); 

  make_moddata(d1kat,d2kat);

  setbeta(Ksp.get_rows(),1,0);

  effectvalues.reserve(nrpar);
  ST::string v;

  for(i=0;i<sizeK1;i++)
    {
    for(j=0;j<sizeK2;j++)
      {
      effectvalues.push_back(effvalue1[i] + " " + effvalue2[j]);
      }
    }

  if (singleblock)
    {

    standnormal = datamatrix(sizeK,1);

    datamatrix xh(sizeK,1);
    datamatrix de(sizeK,1);
    datamatrix ud(sizeK,bandssingle);

    for(i=0;i<sizeK;i++)
      {
      xh(i,0) = 0;
      if (posbeg[i] != -1)
        {

        for (j=posbeg[i];j<=posend[i];j++)
          xh(i,0) += likep->get_weight(index(j,0),0);

        }

      de(i,0) = Ksp(i,i);

      for(j=0;j<ud.cols();j++)
        {
        if (i+j+1 < sizeK)
          ud(i,j) = Ksp(i,i+j+1);
        }

      } // end: for(i=0;i<sizeK;i++)

    XX.reserve(1);
    Kab.reserve(1);

    XX.push_back(bandmatdouble(xh));
    Kab.push_back(bandmatdouble(de,ud));

    prec = bandmatdouble(sizeK,bandssingle);

    muy = datamatrix(sizeK,1);
    betahelp = muy;
    betahelp2 = muy;

    }
  else
    {

    standnormal = datamatrix(sizeK2,1);

    datamatrix xh(sizeK2,1);
    datamatrix de(sizeK2,1);
    datamatrix ud(sizeK2,bands);
    unsigned an=0;
    unsigned en = sizeK2-1;

    Kright.reserve(sizeK1-1);
    Kleft.reserve(sizeK1-1);

    for(i=0;i<sizeK1;i++,an+=sizeK2,en+=sizeK2)
      {
      for (j=an;j<=en;j++)
        {
        xh(j-an,0) = 0;
        if (posbeg[j] != -1)
          {
          for (k=posbeg[j];k<=posend[j];k++)
            xh(j-an,0) += likep->get_weight(index(k,0),0);
          }

        de(j-an,0) = Ksp(j,j);

        for (k=0;k<bands;k++)
          {
          if (j+1+k < sizeK)
            ud(j-an,k) = Ksp(j,j+1+k);
          }

        }

      XX.push_back(bandmatdouble(xh));
      Kab.push_back(bandmatdouble(de,ud));

      if (en < sizeK)   // Kright
        {
        Kright.push_back(Ksp.getBlockasSparse(an,en+1,en+1,sizeK));
        }

      if (an > 0)       // Kleft
        {
        Kleft.push_back(Ksp.getBlockasSparse(an,0,en+1,an));
        }

      } // end: for(i=0;i<sizeK1;i++,an+=sizeK2,en+=sizeK2)


    prec = bandmatdouble(sizeK2,bands);

    muy = datamatrix(sizeK2,1);
    betahelp = muy;
    betahelp2 = muy;

    }  // end: !singleblock


  mu = datamatrix(likep->get_nrobs(),1,0);

  weight = vector<double>(nrpar,1.0/double(nrpar));

  identifiable = false;

  }


void FULLCOND_surf_gaussian::init_maineffects(FULLCOND_nonp_gaussian * mp1,
                                              FULLCOND_nonp_gaussian * mp2,
                                              const ST::string & pnt,
                                              const ST::string & prt)
  {
  mainp1 = mp1;
  mainp2 = mp2;
  centertotal = false;

  fchelprespath = prt;
  datamatrix h(1,1,0);
  fchelp = FULLCOND(optionsp,h,title+"total",nrpar,1,pnt);
  fchelp.setflags(MCMC::norelchange | MCMC::nooutput);
  }



FULLCOND_surf_gaussian::FULLCOND_surf_gaussian(const FULLCOND_surf_gaussian & fc)
  : FULLCOND_nonp_basis(FULLCOND_nonp_basis(fc))
  {
  fchelp = fc.fchelp;
  fchelprespath = fc.fchelprespath;
  singleblock = fc.singleblock;
  centertotal = fc.centertotal;
  mu = fc.mu;
  muy = fc.muy;
  betahelp = fc.betahelp;
  betahelp2 = fc.betahelp2;
  Kab = fc.Kab;
  Kleft = fc.Kleft;
  Kright = fc.Kright;
  sizeK = fc.sizeK;
  sizeK1 = fc.sizeK1;
  sizeK2 = fc.sizeK2;
  standnormal = fc.standnormal;
  XX = fc.XX;
  prec = fc.prec;
  lambda = fc.lambda;
  lambdaold = fc.lambdaold;
  sumx1 = fc.sumx1;
  sumx2 = fc.sumx2;
  }


const FULLCOND_surf_gaussian & FULLCOND_surf_gaussian::operator=(
                                            const FULLCOND_surf_gaussian & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND_nonp_basis::operator=(FULLCOND_nonp_basis(fc));
  fchelp = fc.fchelp;
  fchelprespath = fc.fchelprespath;
  singleblock = fc.singleblock;
  centertotal = fc.centertotal;
  mu = fc.mu;
  muy = fc.muy;
  betahelp = fc.betahelp;
  betahelp2 = fc.betahelp2;
  Kab = fc.Kab;
  Kleft = fc.Kleft;
  Kright = fc.Kright;
  rankK = fc.rankK;
  sizeK = fc.sizeK;
  sizeK1 = fc.sizeK1;
  sizeK2 = fc.sizeK2;
  standnormal = fc.standnormal;
  XX = fc.XX;
  prec = fc.prec;
  lambda = fc.lambda;
  lambdaold = fc.lambdaold;
  sumx1 = fc.sumx1;
  sumx2 = fc.sumx2;
  return *this;
  }

void FULLCOND_surf_gaussian::update(void)
  {

  unsigned i,j,k;
  double * workbetahelp;
  double * workbetahelp2;

  double lambda = sigma2/likep->get_scale();

  double * workbeta = beta.getV();
  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
  for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
    {
    if (*itbeg != -1)
      likep->add_linearpred(-*workbeta,*itbeg,*itend,index,0);
    }

  mu.minus(likep->get_response(),likep->get_linearpred());

  double sigma = sqrt(sigma2);

  double * work = standnormal.getV();

  double * workmuy = muy.getV();

  if (singleblock)
    {

    prec.addto(XX[0],Kab[0],lambda,1);

    for(k=0;k<sizeK;k++,work++)
      *work = sigma*rand_normal();

    prec.solveL(standnormal,betahelp2);

    int * workindex = index.getV();


    for(k=0;k<sizeK;k++,workmuy++)
      {
      *workmuy = 0;
      if (posbeg[k] != -1)
        {
        for(j=posbeg[k];j<=posend[k];j++,workindex++)
          *workmuy+= likep->get_weight(*workindex,0)*mu(*workindex,0);
        *workmuy *= lambda;
        }

      }

    prec.solve(muy,betahelp,0,0);

    workbeta = beta.getV();
    workbetahelp = betahelp.getV();
    workbetahelp2 = betahelp2.getV();
    for (j=0;j<sizeK;j++,workbeta++,workbetahelp++,workbetahelp2++)
      *workbeta = *workbetahelp2 + *workbetahelp;

    } // end: if (singleblock)
  else
    {

    unsigned an =0;
    unsigned en =sizeK2-1;

    for (i=0;i<sizeK1;i++,an+=sizeK2,en+=sizeK2)
      {

      prec.addto(XX[i],Kab[i],1.0/likep->get_scale(),1.0/sigma2);

      work = standnormal.getV();
      for(k=0;k<sizeK2;k++,work++)
        *work = rand_normal();

      prec.solveL(standnormal,betahelp2);

//      int * workindex = index.getV()+posbeg[an];
      workmuy = muy.getV();


      for(k=0;k<sizeK2;k++,workmuy++)
        {
        *workmuy = 0;
        if (posbeg[an+k] != -1)
          {
          for(j=posbeg[an+k];j<=posend[an+k];j++)
            *workmuy+= likep->get_weight(index(j,0),0)*mu(index(j,0),0);
          }
        *workmuy /= likep->get_scale();
        }

      datamatrix muhelp(sizeK2,1,0);

      if (en < sizeK)
        Kright[i].substr_mult(beta,en+1,0,muhelp);
      if (an > 0)
        Kleft[i-1].substr_mult(beta,0,0,muhelp);

      workmuy = muhelp.getV();
      for (k=0;k<sizeK2;k++,workmuy++)
        *workmuy /=sigma2;

      muy.plus(muy,muhelp);

      prec.solve(muy,betahelp,0,0);

      workbeta = beta.getV()+an;
      workbetahelp = betahelp.getV();
      workbetahelp2 = betahelp2.getV();
      for (j=0;j<sizeK2;j++,workbeta++,workbetahelp++,workbetahelp2++)
        *workbeta = *workbetahelp2 + *workbetahelp;

      }  // end: for (i=0;i<sizeK1;i++,an+=sizeK2,en+=sizeK2)

    } // end: !singleblock


  workbeta = beta.getV();
  itbeg = posbeg.begin();
  itend = posend.begin();
  for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
    {
    if (*itbeg != -1)
      likep->add_linearpred(*workbeta,*itbeg,*itend,index,0);
    }


  if (center)
    {
    if (centertotal)
      {
      double m = centerbeta();
      likep->update_intercept(m,column);
      }
    else
      {
      double m = centerbeta2(sumx1,sumx2);
      likep->update_intercept(-m,column);
      mainp1->changebeta(sumx1);
      mainp2->changebeta(sumx2);

      double * workbetahelp = fchelp.getbetapointer();
      double * workbeta = beta.getV();
      double * workbeta1 = mainp1->getbetapointer();
      double * workbeta2;

      unsigned i,j;

      for (i=0;i<sizeK1;i++,workbeta1++)
        {
        workbeta2 = mainp2->getbetapointer();
        for(j=0;j<sizeK2;j++,workbetahelp++,workbeta++,workbeta2++)
          *workbetahelp = *workbeta + *workbeta1 + *workbeta2;
        }

      fchelp.update();  

      }

    }


  acceptance++;

  FULLCOND::update();
  }


bool FULLCOND_surf_gaussian::posteriormode(void)
  {
/*
  unsigned i,j,k;

  double * workbeta = beta.getV();
  vector<int>::iterator itbeg = posbeg.begin();
  vector<int>::iterator itend = posend.begin();
  for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
    {
    if (*itbeg != -1)
      likep->add_linearpred(-*workbeta,*itbeg,*itend,index,0);
    }

  mu.minus(likep->get_response(),likep->get_linearpred());

  double lambda = likep->get_scale()/sigma2;

  prec.addto(XX[0],Kab[0],1,lambda);


//  int * workindex = index.getV();
  double * workmuy = muy.getV();
  for(k=0;k<sizeK;k++,workmuy++)
    {
    *workmuy = 0;
    if (posbeg[k] != -1)
      {
      for(j=posbeg[k];j<=posend[k];j++)
        *workmuy+= mu(index(j,0),0);
      }

    }

  prec.solve(muy,beta,0,0);

  if (center)
    centerbeta(weight);

  workbeta = beta.getV();
  itbeg = posbeg.begin();
  itend = posend.begin();
  for (i=0;i<nrpar;i++,workbeta++,++itbeg,++itend)
    {
    if (*itbeg != -1)
      likep->add_linearpred(*workbeta,*itbeg,*itend,index,0);
    }

  double a = a_invgamma+0.5*rankK;
  double b = b_invgamma+0.5*K.compute_quadform(beta,0);

  sigma2 = b/(a+1);         // mode

  return FULLCOND::posteriormode();

  */
  return true;
  }


void FULLCOND_surf_gaussian::outoptions(void)
  {
  FULLCOND_nonp_basis::outoptions();
  ST::string h;
  if (centertotal)
    h = "total";
  else
    h = "row- and columnwise";

  optionsp->out("  Centermethod: " + h + "\n");
  }


void FULLCOND_surf_gaussian::outresults(void)
  {

  FULLCOND::outresults();

  optionsp->out("  Results are stored in file " + pathresults + "\n");
  optionsp->out("\n");
  if (centertotal == false)
    {
    optionsp->out("  Results for total effect are stored in file " +
                  fchelprespath + "\n");
    optionsp->out("\n");
    }

  unsigned i;

  ofstream outres(pathresults.strtochar());

  ST::string name = title;

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');

  outres << "intnr" << "   ";
  for (i=0;i<datanames.size();i++)
    outres << datanames[i] << "   ";
  outres << name << "mean   ";
  outres << name << "qu" << l1 << "   ";
  outres << name << "qu" << l2 << "   ";
  outres << name << "med   ";
  outres << name << "qu" << u1 << "   ";
  outres << name << "qu" << u2 << "   ";
  outres << endl;

  double * workmean = betamean.getV();
  double * workbetaqu_l1 = betaqu_l1_lower.getV();
  double * workbetaqu_l2 = betaqu_l2_lower.getV();
  double * workbetaqu50 = betaqu50.getV();
  double * workbetaqu_u1 = betaqu_l2_upper.getV();
  double * workbetaqu_u2 = betaqu_l1_upper.getV();
  vector<ST::string>::iterator effit = effectvalues.begin();

   for(i=0;i<nrpar;i++,workmean++,workbetaqu_l1++,workbetaqu_l2++,
          workbetaqu50++,workbetaqu_u1++,workbetaqu_u2++,++effit)
     {
     outres << (i+1) << "   ";
     outres << *effit << "   ";
     outres << *workmean << "   ";
     outres << *workbetaqu_l1 << "   ";
     outres << *workbetaqu_l2 << "   ";
     outres << *workbetaqu50 << "   ";
     outres << *workbetaqu_u1 << "   ";
     outres << *workbetaqu_u2 << "   ";
     outres << endl;
     }


  if (centertotal == false)
    {

    fchelp.outresults();

    ofstream outrestotal(fchelprespath.strtochar());

    optionsp->out("  Results for total effect are stored in file " +
                  fchelprespath + "\n");
    optionsp->out("\n");


    outrestotal << "intnr" << "   ";
    for (i=0;i<datanames.size();i++)
      outrestotal << datanames[i] << "   ";

    outrestotal << name << "mean   ";
    outrestotal << name << "qu" << l1 << "   ";
    outrestotal << name << "qu" << l2 << "   ";
    outrestotal << name << "med   ";
    outrestotal << name << "qu" << u1 << "   ";
    outrestotal << name << "qu" << u2 << "   ";
    outrestotal << endl;

    workmean = fchelp.get_betameanp();
    workbetaqu_l1 = fchelp.get_beta_lower1_p();
    workbetaqu_l2 = fchelp.get_beta_lower2_p();
    workbetaqu50 = fchelp.get_betaqu50p();
    workbetaqu_u1 = fchelp.get_beta_upper2_p();
    workbetaqu_u2 = fchelp.get_beta_upper1_p();
    effit = effectvalues.begin();

   for(i=0;i<nrpar;i++,workmean++,workbetaqu_l1++,
           workbetaqu_l2++,workbetaqu50++,workbetaqu_u1++,workbetaqu_u2++,
           ++effit)
     {
     outrestotal << (i+1) << "   ";
     outrestotal << *effit << "   ";
     outrestotal << *workmean << "   ";
     outrestotal << *workbetaqu_l1 << "   ";
     outrestotal << *workbetaqu_l2 << "   ";
     outrestotal << *workbetaqu50 << "   ";
     outrestotal << *workbetaqu_u1 << "   ";
     outrestotal << *workbetaqu_u2 << "   ";
     outrestotal << endl;
     }
   }

  } // end: outresults




} // end: namespace MCMC
