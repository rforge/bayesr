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



#include "mcmc_nonp.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//----------- class PenaltyMatrix: implementation of member functions ----------
//------------------------------------------------------------------------------


void PenaltyMatrix::make_categories(const datamatrix & moddata,
                                    const unsigned & maxnrint)
  {

  unsigned j;                              // Loopingvariable

  // Initialization of the index, Indexsort of moddata

  index = statmatrix<int>(moddata.rows(),1);
  index.indexinit();
  moddata.indexsort(index,0,moddata.rows()-1,0,0);

  double diff;                             // difference of two succeeding
                                           // observations
  double maxdiff = 1.0/double(maxnrint);
  double beg;                              // first observation in the
                                           // current intervall
  double last;                             // (j-1). observation
  double intbeg;

  double distance =
  moddata(index(moddata.rows()-1,0),0) - moddata(index(0,0),0);

  if (distance == 0)
    errormessages.push_back("ERROR: not enough different covariate values (" + modname +
                     ")\n");
  else
    {

    beg = moddata(index(0,0),0);
    last = moddata(index(0,0),0);
    intbeg = moddata(index(0,0),0);

    posbeg.push_back(0);

    for(j=1;j<moddata.rows();j++)
      {


      diff = (moddata(index(j,0),0)-beg)/distance;

      if (diff > maxdiff)
        {
        weight.push_back(last - intbeg);
        effectvalues.push_back(ST::doubletostring(last));
        effectvdouble.push_back(last);
        posbeg.push_back(j);
        posend.push_back(j-1);

        beg = moddata(index(j,0),0);
        intbeg = last;
        }

      if (j == moddata.rows()-1)
        {
        weight.push_back(moddata(index(j,0),0)-intbeg);
        effectvalues.push_back((ST::doubletostring(moddata(index(j,0),0))));
        effectvdouble.push_back(moddata(index(j,0),0));
        }

      last = moddata(index(j,0),0);


      }  // end: for(j=1;j<moddata.rows();j++)

    posend.push_back(moddata.rows()-1);

    double sum=0;
    for(j=1;j<weight.size();j++)
      sum+= weight[j];

    if (type==RW1)
      sum = double(weight.size()-1)/sum;
    else
      sum = 0.5*double(weight.size()-1)/sum;

    for(j=1;j<weight.size();j++)
      weight[j] *= sum;

    if (posbeg.size() < 6)
      errormessages.push_back("ERROR: not enough different covariate values (" +
                          modname + ")\n");
    }

  }  // end: function make_categories


statmatrix<unsigned> PenaltyMatrix::make_categories2(const datamatrix & moddata,
                                     const unsigned & maxnrint,unsigned & size,
                                     vector<ST::string> & values)
  {

  unsigned j;                              // Loopingvariable

  statmatrix<unsigned> result(moddata.rows(),1,1);

  // Initialization of the index, Indexsort of moddata

  statmatrix<int> index(moddata.rows(),1);
  index.indexinit();
  moddata.indexsort(index,0,moddata.rows()-1,0,0);

  double diff;                             // difference of two succeeding
                                           // observations
  double maxdiff = 1.0/double(maxnrint);
  double beg;                              // first observation in the
                                           // current intervall
  double last;                             // (j-1). observation

  unsigned nr = 1;

  double distance =
  moddata(index(moddata.rows()-1,0),0) - moddata(index(0,0),0);

  beg = moddata(index(0,0),0);
  last = moddata(index(0,0),0);

  for(j=1;j<moddata.rows();j++)
    {

    diff = (moddata(index(j,0),0)-beg)/distance;

    if (diff > maxdiff)
      {
      nr++;
      beg = moddata(index(j,0),0);
      values.push_back((ST::doubletostring(last)));
      }


    if (j == moddata.rows()-1)
      values.push_back((ST::doubletostring(moddata(index(j,0),0))));

    result(index(j,0),0) = nr;
    last = moddata(index(j,0),0);

    }  // end: for(j=1;j<moddata.rows();j++)


  size = result(index(moddata.rows()-1,0),0);
  return result;

  }  // end: function make_categories2


void PenaltyMatrix::make_Kab_list(void)
  {

  unsigned i,j;
  unsigned nrupdate;
  unsigned a,b;
  datamatrix kab;
  datamatrix help;
  datamatrix null(1,1,0);

//  ofstream o("c:\\daten\\kab.raw");

  unsigned sum=0;
  for(i=min;i<=max;i++)
    {
    nrupdate = K.get_rows()/i;
    if ((nrupdate*i) < K.get_rows())
      nrupdate++;
    sum+= nrupdate;
    }


  KAB.reserve(sum);
  KABroot.reserve(sum);
  KABr_sp.reserve(sum);
  KABl_sp.reserve(sum);

  for (i=min;i<=max;i++)
    {

    nrupdate = K.get_rows()/i;
    if ((nrupdate*i) < K.get_rows())
      nrupdate++;

    begin.push_back(KAB.size());
    matquant.push_back(nrupdate);

    for(j=1;j <= nrupdate;j++)
      {

      a = 1+(j-1)*i;
      if (j == nrupdate)
        b = K.get_rows();
      else
        b = j*i;

      kab = (K.getBlock(a-1,a-1,b,b)).inverse();
      KAB.push_back(kab);
      KABroot.push_back(kab.root());

      if (b != K.get_cols())
        {
        help = kab*K.getBlock(a-1,b,b,K.get_cols());
        KABr_sp.push_back(SparseMatrix(help,true));
        }
      else
        KABr_sp.push_back(SparseMatrix());

      if (a != 1)
        {
        help = kab*K.getBlock(a-1,0,b,a-1);
        KABl_sp.push_back(SparseMatrix(help,true));
        }
      else
        KABl_sp.push_back(SparseMatrix());

      }  // end: for(j=1;j <= nrupdate;j++)

    } // end: for (i=min;i<=max;i++)

  }  // end: function make_Kab_list



PenaltyMatrix::PenaltyMatrix(const datamatrix & md, const ST::string & na,
                             const unsigned & nrint, const unsigned & minb,
                             const unsigned & maxb, const fieldtype & ft,const
                             unsigned & per)
  {

  unsigned i;
  type = ft;

  min = minb;
  max = maxb;

  fc_random.reserve(maxb);
  randnorm.reserve(maxb);

  for(i=0;i<maxb;i++)
    {
    fc_random.push_back(datamatrix(i+1,1,0));
    randnorm.push_back(datamatrix(i+1,1,0));
    }

  modname = na;
  make_categories(md,nrint);

  if (errormessages.empty())
    {
    if (type == RW1)
      {
      K = Krw1(weight);
      sizeK=K.get_rows();
      rankK = sizeK-1;
      }
    else if (type == RW2)
      {
      K = Krw2(weight);
      Kband = Krw2band(weight);
      sizeK=K.get_rows();
      rankK = sizeK-2;
      }
    else if (type == seasonal)
      {
      K = Kseason(per,weight.size());
      sizeK = K.get_rows();
      rankK = sizeK - per+1;
      period = per;
      }

    make_Kab_list();

    weight = vector<double>(sizeK,1.0/double(sizeK));

    randnormal = datamatrix(sizeK,1);
    }

  }


PenaltyMatrix::PenaltyMatrix(const datamatrix & md,const ST::string & na,
                             MAP::map & m,const unsigned & minb,
                             const unsigned & maxb)
  {
  type = mrf;

  min = minb;
  max = maxb;

  modname = na;

  unsigned i;

  if (m.polygones_existing())
    polex = true;
  else
    polex = false;

  m.compute_reg(md,posbeg,posend,effectvalues,index);

  if (m.get_errormessages().size() > 0)
    errormessages = m.get_errormessages();
  else
    {

    for(i=0;i<maxb;i++)
      {
      fc_random.push_back(datamatrix(i+1,1,0));
      randnorm.push_back(datamatrix(i+1,1,0));
      }

    K = Kmrf(m);

    sizeK = K.get_rows();
    rankK = sizeK-1;

    weight = vector<double>(sizeK,1.0/double(sizeK));

    make_Kab_list();

    randnormal = datamatrix(sizeK,1);
    }

  }


void PenaltyMatrix::make_moddata(const PenaltyMatrix & p1,
                                 const PenaltyMatrix & p2,
                                 const datamatrix & moddata1,
                                 const datamatrix & moddata2)
  {

  unsigned i,j;


  datamatrix moddata = datamatrix(p1.index.rows(),1);
  for(i=0;i<moddata.rows();i++)
    {
    moddata(i,0) = (  p1.get_category(moddata1(i,0)) )*p2.sizeK +
                   (p2.get_category(moddata2(i,0)) +1) ;
    }

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


void PenaltyMatrix::make_moddata2(const statmatrix<unsigned> & moddata1,
                                  const unsigned & size1,
                                  const statmatrix<unsigned> & moddata2,
                                  const unsigned & size2)
  {

  unsigned i,j;

  datamatrix moddata = datamatrix(moddata1.rows(),1);

  for(i=0;i<moddata.rows();i++)
    {
    moddata(i,0) = (moddata1(i,0)-1) * size2 + moddata2(i,0);
    }

  index = statmatrix<int>(moddata.rows(),1);
  index.indexinit();
  moddata.indexsort(index,0,moddata.rows()-1,0,0);

  posbeg = vector<int>(size1*size2,-1);
  posend = vector<int>(size1*size2,-1);

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


PenaltyMatrix::PenaltyMatrix(const PenaltyMatrix & p1,const PenaltyMatrix & p2,
                             const datamatrix & moddata1,
                             const datamatrix & moddata2,
                             const unsigned & minb, const unsigned & maxb)
  {

  type = mrfkronecker;

  unsigned i,j;

  min = minb;
  max = maxb;

  sizeK = p1.sizeK*p2.sizeK;
  sizeK1 = p1.sizeK;
  sizeK2 = p2.sizeK;
  rankK = p1.rankK*p2.rankK;

  modname = p1.modname + p2.modname;

  make_moddata(p1,p2,moddata1,moddata2);

  K = p1.K.kronecker(p2.K);

//  ofstream out("c:\\daten\\Ktest.raw");

//  K.print(out);


  make_Kab_list();

  randnormal = datamatrix(sizeK,1);

  for(i=0;i<maxb;i++)
    {
    fc_random.push_back(datamatrix(i+1,1,0));
    randnorm.push_back(datamatrix(i+1,1,0));
    }

  effectvalues = vector<ST::string>(sizeK);
  weight = vector<double>(sizeK,1.0/double(sizeK));

  unsigned k=0;
  for(i=0;i<p1.sizeK;i++)
    for(j=0;j<p2.sizeK;j++)
      {
      effectvalues[k] = p1.effectvalues[i] + " " + p2.effectvalues[j];
      k++;
      }

  }



PenaltyMatrix::PenaltyMatrix(const datamatrix & moddata1,
                             const datamatrix & moddata2,
                             const ST::string & na1,const ST::string & na2,
                             const unsigned & minb, const unsigned & maxb,
                             const unsigned & maxnrint, const fieldtype & ft)
  {

  type = ft;
  min = minb;
  max = maxb;

  unsigned i,j;

  vector<ST::string> values1;
  vector<ST::string> values2;

  statmatrix<unsigned> m1 = make_categories2(moddata1,maxnrint,sizeK1,values1);
  statmatrix<unsigned> m2 = make_categories2(moddata2,maxnrint,sizeK2,values2);

  sizeK = sizeK1*sizeK2;
  rankK = sizeK-1;

  make_moddata2(m1,sizeK1,m2,sizeK2);

  K = Kmrflinear(sizeK1,sizeK2);

  make_Kab_list();

  randnormal = datamatrix(sizeK,1);

  for(i=0;i<maxb;i++)
    {
    fc_random.push_back(datamatrix(i+1,1,0));
    randnorm.push_back(datamatrix(i+1,1,0));
    }

  effectvalues = vector<ST::string>(sizeK);
  weight = vector<double>(sizeK,1.0/double(sizeK));


  modname = na1 + na2;

  unsigned k=0;
  for(i=0;i<sizeK1;i++)
    for(j=0;j<sizeK2;j++)
      {
      effectvalues[k] = values1[i] + " " + values2[j];
      k++;
      }

  }



PenaltyMatrix::PenaltyMatrix(const PenaltyMatrix & pm)
  {

  type = pm.type;
  period = pm.period;
  modname = pm.modname;
  index = pm.index;
  posend = pm.posend;
  posbeg = pm.posbeg;
  fc_random = pm.fc_random;
  randnorm = pm.randnorm;
  randnormal = pm.randnormal;


  weight = pm.weight;
  effectvalues = pm.effectvalues;
  effectvdouble = pm.effectvdouble;

  K = pm.K;
  Kband = pm.Kband;
  rankK = pm.rankK;
  sizeK = pm.sizeK;
  sizeK1 = pm.sizeK1;
  sizeK2 = pm.sizeK2;

  KAB = pm.KAB;
  KABl_sp = pm.KABl_sp;
  KABr_sp = pm.KABr_sp;
  KABroot = pm.KABroot;
  min = pm.min;
  max = pm.max;
  begin = pm.begin;
  matquant = pm.matquant;

  errormessages = pm.errormessages;

  polex = pm.polex;

  }


const PenaltyMatrix & PenaltyMatrix::operator=(const PenaltyMatrix & pm)
  {
  if (this == &pm)
    return *this;
  type = pm.type;
  period = pm.period;
  modname = pm.modname;
  index = pm.index;
  posbeg = pm.posbeg;
  posend = pm.posend;

  fc_random = pm.fc_random;
  randnorm = pm.randnorm;
  randnormal = pm.randnormal;

  weight = pm.weight;
  effectvalues = pm.effectvalues;
  effectvdouble = pm.effectvdouble;

  K = pm.K;
  Kband = pm.Kband;

  rankK = pm.rankK;
  sizeK = pm.sizeK;
  sizeK1 = pm.sizeK1;
  sizeK2 = pm.sizeK2;

  KAB = pm.KAB;
  KABl_sp = pm.KABl_sp;
  KABr_sp = pm.KABr_sp;
  KABroot = pm.KABroot;
  min = pm.min;
  max = pm.max;
  begin = pm.begin;
  matquant = pm.matquant;

  errormessages = pm.errormessages;

  polex = pm.polex;

  return *this;
  }


double PenaltyMatrix::compute_Kab_quadform(const datamatrix & beta,
                                           const datamatrix & vec,const unsigned & start,
                                         const unsigned & a,const unsigned & b,
                                         const unsigned & bs)
  {

  datamatrix help(b-a+1,1);

  datamatrix mu(b-a+1,1,0);
  unsigned matnr = begin[bs-min]+ ((a-1)/bs);
  if (a==1)
    KABr_sp[matnr].substr_mult(beta,b,0,mu);
  else if (b == sizeK)
    KABl_sp[matnr].substr_mult(beta,0,0,mu);
  else
    {
    KABr_sp[matnr].substr_mult(beta,b,0,mu);
    KABl_sp[matnr].substr_mult(beta,0,0,mu);
    }

  unsigned i;
  for(i=0;i<b-a+1;i++)
    help(i,0) = vec(start+i,0)-mu(i,0);

  return Kband.compute_quadformblock(help,0,a-1,b-1);;

  }


void PenaltyMatrix::compute_proposal(const datamatrix & beta, const unsigned & bs,
                               const unsigned & a,const unsigned & b,
                               const double & Q,const unsigned & v)
  {
  unsigned matnr = begin[bs-min]+ ((a-1)/bs );
  unsigned l = b-a+1;
  register unsigned i,j;

  double * fc_randwork = fc_random[b-a].getV();
  double * KABrootwork = KABroot[matnr].getV();

  double * randnormwork = randnorm[b-a].getV();

  for(i=0;i<l;i++,randnormwork++)
    *randnormwork = rand_normal();

  for (i=0;i<l;i++,fc_randwork++)
    {
    *fc_randwork = 0;
    randnormwork = randnorm[b-a].getV();
    for(j=0;j<l;j++,KABrootwork++,randnormwork++)
      *fc_randwork += *KABrootwork * *randnormwork;

    *fc_randwork *= Q;

    *fc_randwork += beta(a-1+i,v);

    }

  }


double PenaltyMatrix::compute_quadform_prec(const datamatrix & beta,
                                            const datamatrix & prop,
                                            const bandmatdouble & prec,
                                            const unsigned & a,
                             const unsigned & b,const unsigned & bs)
  {

  unsigned matnr = begin[bs-min]+ ((a-1)/bs);
  unsigned l = b-a+1;

  unsigned i;
  double * fc_randwork = fc_random[b-a].getV();
  double * propwork = prop.getV()+a-1;

  for(i=0;i<l;i++,fc_randwork++,propwork++)
    {
    *fc_randwork = -(*propwork);
    }


  if (a==1)
    KABr_sp[matnr].substr_mult(beta,b,0,fc_random[b-a]);
  else if (b == sizeK)
    KABl_sp[matnr].substr_mult(beta,0,0,fc_random[b-a]);
  else
    {
    KABr_sp[matnr].substr_mult(beta,b,0,fc_random[b-a]);
    KABl_sp[matnr].substr_mult(beta,0,0,fc_random[b-a]);
    }

  return prec.compute_quadform(fc_random[b-a],0);

  }


unsigned PenaltyMatrix::get_category(const double & v) const
  {
  unsigned cat;
  unsigned i=0;
  double value;

  if (type != mrf)
    {

    while (v > effectvdouble[i])
      i++;

    cat = i;

    }
  else if (type==mrf)
    {

    effectvalues[i].strtodouble(value);

    while (v != value)
      {
      i++;
      effectvalues[i].strtodouble(value);
      }

    cat = i;

    }

  return cat;
  }



//------------------------------------------------------------------------------
//---------- CLASS: FULLCOND_nonp (implementation of member fucntions) ---------
//------------------------------------------------------------------------------

void FULLCOND_nonp::init_priorassumptions(const ST::string & na)
    {
    if(Pmatrix->get_type()==MCMC::RW1)
       priorassumptions.push_back("first order random walk");
    else if(Pmatrix->get_type()==MCMC::RW2)
       priorassumptions.push_back("second order random walk");
    else if(Pmatrix->get_type()==MCMC::mrf)
       priorassumptions.push_back("Markov random field");
    else if(Pmatrix->get_type()==MCMC::seasonal)
       priorassumptions.push_back("seasonal component");
    }


void FULLCOND_nonp::init_name(const ST::string & na)
    {
    FULLCOND::init_name(na);
    ST::string temp = "\\_";
    ST::string helpname = na.insert_string_char('_', temp);
    if (Pmatrix->get_type()==MCMC::seasonal)
      term_symbolic = "f^{Season}_{" +  helpname + "}("+helpname+")";
    else
      term_symbolic = "f_{" +  helpname + "}("+helpname+")";

    if (column > 0)
      priorassumptions.push_back("$" + term_symbolic + "$" +
       " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");
    init_priorassumptions(na);
    }


void FULLCOND_nonp::init_names(const vector<ST::string> & na)
    {
    FULLCOND::init_names(na);
    if (na.size()==1)
      {
      ST::string temp = "\\_";
      ST::string helpname = na[0].insert_string_char('_', temp);
      if (Pmatrix->get_type()==MCMC::seasonal)
        term_symbolic = "f^{Season}_{" +  helpname + "}("+helpname+")";
      else
        term_symbolic = "f_{" +  helpname + "}("+helpname+")";
      }
    else
      {
      ST::string temp = "\\_";
      ST::string helpname1 = na[0].insert_string_char('_', temp);
      ST::string helpname2 = na[1].insert_string_char('_', temp);
      if (Pmatrix->get_type()==MCMC::seasonal)
        term_symbolic = "f^{Season}_{" +  helpname1 + "}("+helpname1+") \\cdot "
                        + helpname2;
      else
        term_symbolic = "f_{" +  helpname1 + "}("+helpname1+") \\cdot "
                        + helpname2;
      }

    if (column > 0)
      priorassumptions.push_back("$" + term_symbolic + "$" +
       " (" + ST::inttostring(column+1) + ". \\mbox{ } response \\mbox{ } category)");
    else
      priorassumptions.push_back("$" + term_symbolic + "$");
    priorassumptions.push_back(term_symbolic);
    init_priorassumptions(na[0]);
    }


/*
void FULLCOND_nonp::init_name(const ST::string & na)
    {
    FULLCOND::init_name(na);
    ST::string helpname = na.insert_string_char('_', "\\_");
    term_symbolic = "f_{" +  helpname + "}("+helpname+")";
    priorassumptions.push_back("$" + term_symbolic + "$:");
    init_priorassumptions(na);
    }

void FULLCOND_nonp::init_names(const vector<ST::string> & na)
    {
    FULLCOND::init_names(na);
    if (na.size()==1)
      {
      ST::string helpname = na[0].insert_string_char('_', "\\_");
      term_symbolic = "f_{" +  helpname + "}("+helpname+")";
      }
    else
      {
      ST::string helpname1 = na[0].insert_string_char('_', "\\_");
      ST::string helpname2 = na[1].insert_string_char('_', "\\_");
      term_symbolic = "f_{" +  helpname1 + "}("+helpname1+") \\cdot " + helpname2;
      }
    priorassumptions.push_back("$" + term_symbolic + "$:");
    init_priorassumptions(na[0]);
    }
*/


double FULLCOND_nonp::compute_quadform(void)
  {
  return Pmatrix->compute_quadform(beta,0);
  }


FULLCOND_nonp::FULLCOND_nonp(MCMCoptions * o,DISTRIBUTION * dp,PenaltyMatrix* K,
                             FULLCOND_const * fcco,
                             const double & l,
                             const ST::string & fp, const ST::string & pres,
                             const ST::string & t, const ST::string & mn,
                             const unsigned & c)
  : FULLCOND(o,datamatrix(1,2),K->get_modname(),K->get_sizeK(),1,fp)
  {

  fcconst = fcco;
  title = t;
  lambda = l;
  sigma2=10;
  if (K->get_type() != mrf)
    {
    plotstyle = plotnonp;
    }
  else if(K->get_polex()==false)
    plotstyle = drawmapgraph;
  else
    plotstyle = drawmap;

  mapname = mn;
  pathresult = pres;
  pathcurrent = pres;
  likep = dp;
  column = c;
  centerm = total;
  Pmatrix = K;
  if (Pmatrix->get_errormessages().empty())
    {
    if (K->get_type() != seasonal)
      {
      identifiable = false;
      }
    else
      {
      identifiable = true;
      }

    pathresult = pres;
    pathcurrent = pres;
    accept = datamatrix(beta.rows(),beta.cols(),0);
    effectmod = datamatrix(1,1);
    varcoeff=false;

    effectvalues = Pmatrix->get_values();
    effectvdouble = Pmatrix->get_valuesd();

    unsigned i;
    for (i=0;i<nrpar;i++)
      {
      if (Pmatrix->get_posbeg(i) == -1)
        optionsp->out("NOTE: no observations for covariate value " +
                      effectvalues[i] + "\n");

      }

    weight = vector<double>(nrpar,1.0/double(nrpar));
    }
  else
    {
    errors = Pmatrix->get_errormessages();
    }

  }


// varying coefficients
FULLCOND_nonp::FULLCOND_nonp(MCMCoptions * o,DISTRIBUTION * dp,PenaltyMatrix* K,
                             FULLCOND_const * fcco,
                             const double & l,
                             const datamatrix & effmod, const ST::string & ti,
                             const ST::string & fp, const ST::string & pres,
                              const ST::string & mn,
                             const unsigned & c)
  : FULLCOND(o,datamatrix(1,2),K->get_modname(),K->get_sizeK(),1,fp)
  {

  fcconst = fcco;
  lambda = l;
  sigma2=10;
  if (K->get_type() != mrf)
    {
    plotstyle = plotnonp;
    }
  else if(K->get_polex()==false)
    plotstyle = drawmapgraph;
  else
    plotstyle = drawmap;

  mapname = mn;
  pathresult = pres;
  pathcurrent = pres;
  likep = dp;
  column = c;
  centerm = total;
  Pmatrix = K;
  if (Pmatrix->get_errormessages().empty())
    {
    identifiable = true;

    pathresult = pres;
    pathcurrent = pres;
    accept = datamatrix(beta.rows(),beta.cols(),0);

    title = ti;
    effectmod = effmod;
    varcoeff = true;

    effectvalues = Pmatrix->get_values();
    effectvdouble = Pmatrix->get_valuesd();
    }
  else
    {
    errors = Pmatrix->get_errormessages();
    }

  }


FULLCOND_nonp::FULLCOND_nonp(const FULLCOND_nonp & fc)
  : FULLCOND(FULLCOND(fc))
  {
  fcconst = fc.fcconst;
  likep = fc.likep;
  centerm = fc.centerm;
  Pmatrix = fc.Pmatrix;
  sigma2 = fc.sigma2;
  lambda = fc.lambda;
  accept = fc.accept;
  effectmod = fc.effectmod;
  varcoeff = fc.varcoeff;
  effectvalues = fc.effectvalues;
  effectvdouble = fc.effectvdouble;
  mapname = fc.mapname;
  }


const FULLCOND_nonp & FULLCOND_nonp::operator=(const FULLCOND_nonp & fc)
  {
  if (this == &fc)
    return *this;
  FULLCOND::operator=(FULLCOND(fc));
  fcconst = fc.fcconst;
  likep = fc.likep;
  centerm = fc.centerm;
  Pmatrix = fc.Pmatrix;
  sigma2 = fc.sigma2;
  lambda = fc.lambda;
  accept = fc.accept;
  effectmod = fc.effectmod;
  varcoeff = fc.varcoeff;
  effectvalues = fc.effectvalues;
  effectvdouble = fc.effectvdouble;
  mapname = fc.mapname;
  return *this;
  }


double FULLCOND_nonp::centerbeta(vector<double> & weight)
  {
  unsigned i;
  double sum=0;
  double * workbeta = beta.getV();

#if defined(MICROSOFT_VISUAL)
  for (i=0;i<nrpar;i++,workbeta++)
    sum+= weight[i] * *workbeta;
#else
  vector<double>::iterator weightwork = weight.begin();

  for (i=0;i<nrpar;i++,workbeta++,++weightwork)
    sum+= *weightwork * *workbeta;
#endif

  workbeta = beta.getV();

  for (i=0;i<nrpar;i++,workbeta++)
    *workbeta-= sum;

  return sum;
  }


void FULLCOND_nonp::centerbeta2(void)
  {
  unsigned i,j,k;
  double betaold;

  double * workbeta = beta.getV();
  unsigned rK1 = Pmatrix->get_sizeK1();
  unsigned rK2 = Pmatrix->get_sizeK2();
  datamatrix sumk1(rK1,1,0);
  datamatrix sumk2(rK2,1,0);
  double sumtotal=0;

  for (j=0;j<rK1;j++)
    {
    for (i=0;i<rK2;i++,workbeta++)
      {
      sumtotal+= *workbeta;
      sumk1(j,0) += *workbeta;
      sumk2(i,0) += (1.0/double(rK1))* *workbeta;
      }

    sumk1(j,0) = sumk1(j,0)/double(rK2);

    }

  sumtotal/= double(rK1*rK2);

  workbeta = beta.getV();

  k=0;
  for (j=0;j<rK1;j++)
    {
    for (i=0;i<rK2;i++,workbeta++)
      {
      betaold = *workbeta;
      *workbeta-= sumk1(j,0)+sumk2(i,0)-sumtotal;
      if (Pmatrix->get_posbeg(k) != -1)
        likep->add_linearpred(*workbeta-betaold,Pmatrix->get_posbeg(k),
                               Pmatrix->get_posend(k),Pmatrix->get_index(),column);
      k++;
      }

    }

  }



void FULLCOND_nonp::outoptions(void)
  {

  if (Pmatrix->get_type()==seasonal)
    optionsp->out("  OPTIONS FOR FLEXIBLE SEASONAL COMPONENT: " + title + "\n",true);
  else
    optionsp->out("  OPTIONS FOR NONPARAMETRIC TERM: " + title + "\n",true);
  optionsp->out("\n");
  optionsp->out("  Prior: " + Pmatrix->get_typeasstring() + "\n");
  optionsp->out("  Minimum blocksize: " + ST::inttostring(Pmatrix->get_minsize())
                  + "\n");
  optionsp->out("  Maximum blocksize: " + ST::inttostring(Pmatrix->get_maxsize())
                  + "\n");
  if (Pmatrix->get_type()==seasonal)
  optionsp->out("  Period of seasonal effect: " +
                ST::inttostring(Pmatrix->get_period()) + "\n");

  optionsp->out("\n");

  }



void FULLCOND_nonp::update(void)
  {

  unsigned blocksize = int(Pmatrix->get_minsize() +
	  uniform()*(Pmatrix->get_maxsize()-Pmatrix->get_minsize()+1));

  double u;
  unsigned an = 1;
  unsigned en = blocksize;
  int beg;
  int end;

  double logold;
  double logprop;
  double * workbeta;
  double value;
  unsigned j,k,l;
  int ind;

  if (!varcoeff)
    {
    for(j=0;j<Pmatrix->get_nrblocks(blocksize);j++)
      {

      nrtrials++;

      Pmatrix->compute_fc(beta,blocksize,an,en,sqrt(sigma2));

      logold = 0;
      logprop = 0;
      for (k=an;k<=en;k++)
        {
        beg = Pmatrix->get_posbeg(k-1);
        end = Pmatrix->get_posend(k-1);

        if (beg != -1)
          logold += likep->loglikelihood(beg,end,Pmatrix->get_index());

        if (beg != -1)
          likep->add_linearpred(
          Pmatrix->get_fc_random()[en-an](k-an,0)-beta(k-1,0),beg,end,
                                  Pmatrix->get_index(),column);


        if (beg != -1)
          logprop += likep->loglikelihood(beg,end,Pmatrix->get_index());

        }  // end: for (k=an;k<=en;k++)

      u = log(uniform());

      if (u <= (logprop-logold))    // accept
        {
        workbeta = beta.getV()+an-1;
        for(k=an-1;k<en;k++,workbeta++)
          {
          *workbeta = Pmatrix->get_fc_random()[en-an](k-an+1,0);
          accept(k,0)++;
          }

	    acceptance++;
        }
      else
        {
        for (k=an;k<=en;k++)
         if (Pmatrix->get_posbeg(k-1) != -1)
           {
           likep->add_linearpred(
           beta(k-1,0)-Pmatrix->get_fc_random()[en-an](k-an,0),
           Pmatrix->get_posbeg(k-1),
           Pmatrix->get_posend(k-1),Pmatrix->get_index(),column);
           }
        }

      an+=blocksize;
      if (j ==  Pmatrix->get_nrblocks(blocksize)-2)
        en = Pmatrix->get_sizeK();
      else
        en+=blocksize;

      }    // end: for(j=0;j<Pmatrix->getnrblock(blocksize);j++)

    } // end:   if (!varcoeff)
  else   // varcoeff == true
    {
    for(j=0;j<Pmatrix->get_nrblocks(blocksize);j++)
      {

      nrtrials++;

      Pmatrix->compute_fc(beta,blocksize,an,en,sqrt(sigma2));

      logold = 0;
      logprop = 0;
      for (k=an;k<=en;k++)
        {
        beg = Pmatrix->get_posbeg(k-1);
        end = Pmatrix->get_posend(k-1);

        if (beg != -1)
          logold += likep->loglikelihood(beg,end,Pmatrix->get_index());

        if (beg != -1)
          {
          value = Pmatrix->get_fc_random()[en-an](k-an,0)-beta(k-1,0);
          for (l=beg;l<=end;l++)
            {
            ind = Pmatrix->get_index()(l,0);
            likep->add_linearpred(value*effectmod(ind,0),ind,column);
            }
          }


        if (beg != -1)
          logprop += likep->loglikelihood(beg,end,Pmatrix->get_index());

        }  // end: for (k=an;k<=en;k++)

      u = log(uniform());

      if (u <= (logprop-logold))    // accept
        {
        workbeta = beta.getV()+an-1;
        for(k=an-1;k<en;k++,workbeta++)
          {
          *workbeta = Pmatrix->get_fc_random()[en-an](k-an+1,0);
          accept(k,0)++;
          }

	    acceptance++;
        }
      else
        {
        for (k=an;k<=en;k++)
          if (Pmatrix->get_posbeg(k-1) != -1)
            {
            value = beta(k-1,0)-Pmatrix->get_fc_random()[en-an](k-an,0);
            beg = Pmatrix->get_posbeg(k-1);
            end = Pmatrix->get_posend(k-1);
            for (l=beg;l<=end;l++)
              {
              ind = Pmatrix->get_index()(l,0);
              likep->add_linearpred(value*effectmod(ind,0),ind,column);
              }
            }

        }

      an+=blocksize;
      if (j ==  Pmatrix->get_nrblocks(blocksize)-2)
        en = Pmatrix->get_sizeK();
      else
        en+=blocksize;

      }    // end: for(j=0;j<Pmatrix->getnrblock(blocksize);j++)

    } // end varcoeff


  if (center)
    {
    double m = FULLCOND::centerbeta();
    fcconst->update_intercept(m);
    }

  FULLCOND::update();
  } // end: function update



void FULLCOND_nonp::outresults(void)
  {
  FULLCOND::outresults();

  optionsp->out("  Results are stored in file\n");
  optionsp->out("  " + pathcurrent + "\n");
  optionsp->out("\n");
  if (Pmatrix->get_type() == MCMC::mrf)
    {
    #if defined(JAVA_OUTPUT_WINDOW)

    if (Pmatrix->get_polex() == true)
      {
      optionsp->out("  Postscript file is stored in file\n");
      ST::string psfile = pathcurrent.substr(0,pathcurrent.length()-4) + ".ps";
      optionsp->out("  " + psfile + "\n");
      optionsp->out("\n");
      }

    optionsp->out("  Results may be visualized using method 'drawmap'\n");
    optionsp->out("  Type for example: objectname.drawmap " +
    ST::inttostring(fcnumber) + "\n");
    optionsp->out("\n");
    #else
    optionsp->out("  Results may be visualized using the R function");
    optionsp->out(" 'drawmap'\n");
    optionsp->out("\n");
    #endif
    }
  else
    {
    #if defined(JAVA_OUTPUT_WINDOW)
    optionsp->out("  Postscript file is stored in file\n");
    ST::string psfile = pathcurrent.substr(0,pathcurrent.length()-4) + ".ps";
    optionsp->out("  " + psfile + "\n");
    optionsp->out("\n");
    optionsp->out("  Results may be visualized in BayesX using method 'plotnonp'\n");
    optionsp->out("  Type for example: objectname.plotnonp " +
    ST::inttostring(fcnumber) + "\n");
    optionsp->out("\n");
    #else
    ST::string doublebackslash = "/";
    ST::string pathresultsplus = pathcurrent.insert_string_char('\\',doublebackslash);
    optionsp->out("  Results may be visualized using the R function 'plotnonp'");
    optionsp->out("\n");
    optionsp->out("  Type for example:\n");
    optionsp->out("\n");
    optionsp->out("  plotnonp(\""+ pathresultsplus + "\")\n");
    optionsp->out("\n");
    #endif
    }
  optionsp->out("\n");

  unsigned i;

  ofstream outres(pathcurrent.strtochar());

  ST::string name = Pmatrix->get_modname();

  ST::string l1 = ST::doubletostring(lower1,4);
  ST::string l2 = ST::doubletostring(lower2,4);
  ST::string u1 = ST::doubletostring(upper1,4);
  ST::string u2 = ST::doubletostring(upper2,4);
  l1 = l1.replaceallsigns('.','p');
  l2 = l2.replaceallsigns('.','p');
  u1 = u1.replaceallsigns('.','p');
  u2 = u2.replaceallsigns('.','p');


  outres << "intnr" << "   ";
  outres << datanames[0] << "   ";
  outres << "pmean   ";
  outres << "pqu" << l1 << "   ";
  outres << "pqu" << l2 << "   ";
  outres << "pmed   ";
  outres << "pqu" << u1 << "   ";
  outres << "pqu" << u2 << "   ";
  outres << "pcat" << level1 << "   ";
  outres << "pcat" << level2 << "   ";
  outres << endl;

  double * workmean = betamean.getV();
  double * workbetaqu_l1_lower_p = betaqu_l1_lower.getV();
  double * workbetaqu_l2_lower_p = betaqu_l2_lower.getV();
  double * workbetaqu_l1_upper_p = betaqu_l1_upper.getV();
  double * workbetaqu_l2_upper_p = betaqu_l2_upper.getV();
  double * workbetaqu50 = betaqu50.getV();

  vector<ST::string>::iterator effit = effectvalues.begin();

  for(i=0;i<nrpar;i++,++effit,workmean++,workbetaqu_l1_lower_p++,
                           workbetaqu_l2_lower_p++,workbetaqu50++,
                           workbetaqu_l1_upper_p++,workbetaqu_l2_upper_p++)
    {
    outres << (i+1) << "   ";
    outres << *effit << "   ";
    outres << *workmean << "   ";
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

    outres << endl;
    }

  }

void FULLCOND_nonp::get_effectmatrix(datamatrix & e,
                                        vector<ST::string> & enames,unsigned be,
                                        unsigned en, effecttype t)
  {

  statmatrix<int> ind = Pmatrix->get_index();

  unsigned i;
  unsigned k=0;
  unsigned anf, end, j;
  for(i=0; i<beta.rows(); i++)
    {
    anf = Pmatrix->get_posbeg(i);
    end = Pmatrix->get_posend(i);
    for(j=anf; j<=end; j++)
      {
      e(j,0) = beta(ind(i,0),0);
      k++;
      }
    }

  }

void FULLCOND_nonp::get_acceptance(const ST::string & filename)
  {

  ofstream out(filename.strtochar());
  assert( !out.fail() );
  get_acceptance(out);

  }


void FULLCOND_nonp::get_acceptance(ostream & out)
  {

  out << "parnr ";
  unsigned i,j;
  for (i=0;i<beta.cols();i++)
    out << "rate_" << (i+1) << " ";

  out << endl;

  for(j=0;j<beta.rows();j++)
    {
    out << (j+1) << " ";
    for(i=0;i<beta.cols();i++)
      out << accept(j,i) << " ";
    out << endl;
    }

  }


void PenaltyMatrix::compute_fc(const datamatrix & beta, const unsigned & bs,
                               const unsigned & a,const unsigned & b,
                               const double & Q,const unsigned & v)
  {
  unsigned matnr = begin[bs-min]+ ((a-1)/bs );
  unsigned l = b-a+1;
  register unsigned i,j;

  double * fc_randwork = fc_random[b-a].getV();
  double * KABrootwork = KABroot[matnr].getV();

  double * randnormwork = randnorm[b-a].getV();

  for(i=0;i<l;i++,randnormwork++)
    *randnormwork = rand_normal();

  for (i=0;i<l;i++,fc_randwork++)
    {
    *fc_randwork = 0;
    randnormwork = randnorm[b-a].getV();
    for(j=0;j<l;j++,KABrootwork++,randnormwork++)
      *fc_randwork += *KABrootwork * *randnormwork;

    *fc_randwork *= Q;
    }

  compute_mu(beta,bs,a,b,v);

  }



void PenaltyMatrix::compute_mu2(const datamatrix & beta,datamatrix & res,
      const unsigned & resstart,
      const unsigned & bs,const unsigned & a,const unsigned & b,
      const unsigned & v)
  {
  unsigned matnr = begin[bs-min]+ ((a-1)/bs);

  if (a==1)
    KABr_sp[matnr].substr_mult(beta,b,v,res,resstart);
  else if (b == sizeK)
    KABl_sp[matnr].substr_mult(beta,0,v,res,resstart);
  else
    {
    KABr_sp[matnr].substr_mult(beta,b,v,res,resstart);
    KABl_sp[matnr].substr_mult(beta,0,v,res,resstart);
    }
  }


void PenaltyMatrix::compute_mu(const datamatrix & beta,
      const unsigned & bs,const unsigned & a,const unsigned & b,
      const unsigned & v)
  {
  unsigned matnr = begin[bs-min]+ ((a-1)/bs);

  if (a==1)
    KABr_sp[matnr].substr_mult(beta,b,v,fc_random[b-a]);
  else if (b == sizeK)
    KABl_sp[matnr].substr_mult(beta,0,v,fc_random[b-a]);
  else
    {
    KABr_sp[matnr].substr_mult(beta,b,v,fc_random[b-a]);
    KABl_sp[matnr].substr_mult(beta,0,v,fc_random[b-a]);
    }
  }



void PenaltyMatrix::compute_fc2(const datamatrix & beta,  datamatrix & res,
                                const unsigned & a,const unsigned & b,
                                const unsigned & bs, const double & Q,
                                const unsigned & v)
  {

  register unsigned i,j;

  unsigned matnr = begin[bs-min]+ ((a-1)/bs );
  unsigned l = b-a+1;

  double * fc_randwork = res.getV()+a-1;
  double * KABrootwork = KABroot[matnr].getV();

  double * randnormwork = randnorm[b-a].getV();

  for(i=0;i<l;i++,randnormwork++)
    *randnormwork = rand_normal();

  for (i=0;i<l;i++,fc_randwork++)
    {
    *fc_randwork = 0;
    randnormwork = randnorm[b-a].getV();
    for(j=0;j<l;j++,KABrootwork++,randnormwork++)
      *fc_randwork += *KABrootwork * *randnormwork;

    *fc_randwork *= Q;
    }

  compute_mu2(beta,res,a-1,bs,a,b,v);

  }


void PenaltyMatrix::compute_proposal2(bandmatdouble & prec,const datamatrix & muab,
                            const datamatrix & beta, datamatrix & res,
                            const unsigned & bs,
                            const unsigned & a,const unsigned & b,
                            const unsigned & v)
  {

  double * randnormwork = randnorm[b-a].getV();

  unsigned i;
  unsigned l = b-a+1;
  for(i=0;i<l;i++,randnormwork++)
    *randnormwork = rand_normal();

  prec.solveL(randnorm[b-a],fc_random[b-a]);

   for (i=0;i<l;i++)
     res(a-1+i,0) = fc_random[b-a](i,0) + muab(i,0);

//  compute_mu2(beta,res,a-1,bs,a,b,v);

  }


} // end: namespace MCMC


