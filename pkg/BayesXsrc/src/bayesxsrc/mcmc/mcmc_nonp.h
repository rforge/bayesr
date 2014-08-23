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



#if !defined (MCMCnonp_INCLUDED)

#define MCMCnonp_INCLUDED

#include"../export_type.h"
#include"mcmc.h"
#include"fullcond.h"
#include"distribution.h"
#if !defined (MAP_INCLUDED)
#include"map.h"
#endif
#include"bandmat.h"
#include"bandmat_penalty.h"
#include"mcmc_nonpbasis.h"

namespace MCMC
{

//------------------------------------------------------------------------------
//-------------------------- CLASS: PenaltyMatrix ------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE PenaltyMatrix
  {

  protected:

  fieldtype type;                    // type of the markov random field
                                     // (first or second order random walk,
                                     //  seasonal, Markov random field or
                                     //  Markov random field based on the
                                     //  kronecker product of two penalty
                                     //  matrices)

  unsigned period;                   // period of seasonal effect

  ST::string modname;                // Name of the effectmodifier
                                     // (usefull for estimation output)

  statmatrix<int> index;             // rank vector of moddata
                                     // access to the i th order observation:
                                     // moddata(index(i,0),0)


  vector<int>     posbeg;            // Let x_1 < x_2 < ... < x_m be the
  vector<int>     posend;            // m different observations in moddata
                                     // Then posbeg[i] - posend[i] indicates
                                     // the number and position of observed
                                     // x_i
                                     // access to all observations with value
                                     // x_i:
                                     // for (j=posbeg[i];j<=posend[i];j++)
                                     //   {
                                     //   ...
                                     //   moddata(index(j,0),0)
                                     //   ...
                                     //   }

  vector<double> weight;             // Length of the intervals between
                                     // successive different observations
  vector<ST::string> effectvalues;       // vector of different values of moddata

  vector<double> effectvdouble;

  SparseMatrix K;                    // Penalty Matrix
  bandmatdouble Kband;
  unsigned rankK;                    // rank of the penalty matrix
  unsigned sizeK;                    // number of rows/columns of K;

  unsigned sizeK1;
  unsigned sizeK2;

  unsigned min;                      // Minimum Blocksize
  unsigned max;                      // Maximum Blocksize

  vector<datamatrix> KAB;            // vector of all possible K_ab^-1
  vector<datamatrix> KABroot;        // vector of all possible sqrt(K_ab^-1)

  vector<SparseMatrix> KABr_sp;      // vector of all possible K_left matrices
                                     // (sparse matrix version)
  vector<SparseMatrix> KABl_sp;      // vector of all possible K_right matrices
                                     // (sparse matrix version)

  vector<unsigned> begin;            // begin[size-min] gives the Position of
                                     // the first KAB,KABl,KABr,KABroot Matrices
                                     // for blocksize 'size'
  vector<unsigned> matquant;         // matquant[size-min] gives the number of
                                     // blocks for blocksize 'size'


  vector<datamatrix> fc_random;
  vector<datamatrix> randnorm;

  datamatrix randnormal;


  vector<ST::string> errormessages;

  bool polex;

  // FUNCTION: make_categories
  // TASK: devides the data in moddata into categories, maximal number of
  //       categories is 'maxnrint'
  // ADDITIONAL NOTE: Implementation is independent of the type of MRF

  void make_categories(const datamatrix & moddata,const unsigned & maxintnr);

  statmatrix<unsigned> make_categories2(const datamatrix & moddata,
                                        const unsigned & maxintnr,
                                        unsigned & size,
                                        vector<ST::string> & values);


  // FUNCTION: make_moddata
  // TASK: computes moddata, if the penalty matrix is the kronecker product
  //       of p1 and p2

  void make_moddata(const PenaltyMatrix & p1,const PenaltyMatrix & p2,
                    const datamatrix & moddata1, const datamatrix & moddata2);


  void make_moddata2(const statmatrix<unsigned> & moddata1,
                                  const unsigned & size1,
                                  const statmatrix<unsigned> & moddata2,
                                  const unsigned & size2);


  // FUNCTION: make_Kab_list
  // TASK: creates KAB, KAbroot,KABr,KABl matrices
  // NOTE: functions depends on:
  //       - Penalty matrix 'K'
  //       - Minimum and Maximum Block size 'min' and 'max'
  //       - Type of MRF 'type'

  void make_Kab_list(void);

  public:


  // DEFAULT CONSTRUCTOR

  PenaltyMatrix(void)
    {
    }

  // CONSTRUCTOR 1 (for RW1,RW2,seasonal component)
  // TASK: md                     Effectmodifier
  //       na      = name         Name of the Effectmodifier
  //       nrint                  Maximum number of intervals
  //       minb    = min          Minumum block size
  //       maxb    = max          Maximum block size
  //       ft      = type         Type of MRF (should be RW1,RW2 or seasonal)
  //       per                    period of seasonal effect, if type = seasonal
  //      function devides moddata into intervals (maximum number 'maxnrint',
  //      computes the Penalty Matrix 'K', KAB,KABroot,KABr and KABl matrices

  PenaltyMatrix(const datamatrix & md,const ST::string & na,const unsigned & nrint,
                const unsigned & minb,const unsigned & maxb,
                const  fieldtype & ft,const unsigned & per = 12);

  // CONSTRUCTOR 2  (for geographical maps)
  // TASK: md                      covariate vector
  //       na                      Name of the covariate
  //       m                       map object
  //       minb = min              Minimum block size
  //       maxb = max              Maximum block size

  PenaltyMatrix(const datamatrix & md,const ST::string & na, MAP::map & m,
                const unsigned & minb, const unsigned & maxb);

  // CONSTRUCTOR 3 (for 2 dimensional surfaces)
  // TASK: computes the kroneckerproduct of the two penalty matrices p1 and p2
  //       minb = min                Minimum block size
  //       maxb = max                Maximum block size

  PenaltyMatrix(const PenaltyMatrix & p1,const PenaltyMatrix & p2,
                const datamatrix & moddata1,const datamatrix & moddata2,
                const unsigned & minb, const unsigned & maxb);

  // CONSTRUCTOR 4 (for 2 dimensional surfaces)
  // TASK: If t = mrflinear, t = mrfquadratic8, t =
  //       mrfquadratic12 the penalty matrix corresponds to a local linear
  //       polynomial, a local quadratic with the 8 nearest neigbors or
  //       a local quadractic ploynomial with the 12 nearest neighbors
  //       minb = min                Minimum block size
  //       maxb = max                Maximum block size
  //       t = type                  Type of the mrf
  //       (allowed: mrflinear, mrfquadratic8, mrfquadratic12)

  PenaltyMatrix(const datamatrix & moddata1,const datamatrix & moddata2,
                const ST::string & na1,const ST::string & na2,
                const unsigned & minb, const unsigned & maxb,
                const unsigned & maxnrint, const fieldtype & ft);

  // COPY CONSTRUCTOR

  PenaltyMatrix(const PenaltyMatrix & pm);


  // OVERLOADED ASSIGNMENT OPERATOR

  const PenaltyMatrix & operator=(const PenaltyMatrix & pm);

  // FUNCTION: compute mu
  // TASK: computes the mean of the conditional prior f[a,b]
  //       'beta' is the current state of the Markov chain
  //        bs = blocksize
  //        v = column of beta

  void compute_mu(const datamatrix & beta, const unsigned & bs,
                  const unsigned & a,const unsigned & b,
                  const unsigned & v);

  void compute_mu2(const datamatrix & beta, datamatrix & res,
                   const unsigned & resstart,const unsigned & bs,
                  const unsigned & a,const unsigned & b,
                  const unsigned & v);

  // FUNCTION: compute_fc
  // TASK:
  // ADDITIONAL NOTE: implementation is independent of the type of MRF

  void compute_fc(const datamatrix & beta, const unsigned & bs,
                  const unsigned & a,const unsigned & b,
                  const double & Q,const unsigned & v=0);

  void compute_fc2(const datamatrix & beta,datamatrix & res, const unsigned & bs,
                   const unsigned & a,const unsigned & b,
                  const double & Q,const unsigned & v=0);


double compute_Kab_quadform(const datamatrix & beta,
                            const datamatrix & vec,const unsigned & start,
                            const unsigned & a,const unsigned & b,
                            const unsigned & bs);

double compute_quadform_prec(const datamatrix & beta,const datamatrix & prop,
                             const bandmatdouble & prec,const unsigned & a,
                             const unsigned & b,const unsigned & bs);


  // FUNCTION: compute_proposal
  // TASK:
  // ADDITIONAL NOTE: implementation is independent of the type of MRF

  void compute_proposal(const datamatrix & beta, const unsigned & bs,
                  const unsigned & a,const unsigned & b,
                  const double & Q,const unsigned & v=0);


  void compute_proposal2(bandmatdouble & prec, const datamatrix & muab,
                         const datamatrix & beta, datamatrix & res ,const unsigned & bs,
                  const unsigned & a,const unsigned & b,const unsigned & v=0);



   // FUNCTION: compute_quadform
   // TASK: returns beta(.,v)' K beta(.,v) where K is the penalty matrix

  double compute_quadform(const datamatrix & beta,const unsigned & v=0)
    {
    return K.compute_quadform(beta,v);
    }

  // FUNCTION: rank_K
  // TASK: returns the rank of K

  const unsigned & rank_K(void) const
    {
    return rankK;
    }


  const unsigned & get_bandsize(void) const
    {
    return Kband.bandsize();
    }

  const bandmatdouble & get_Kband(void) const
    {
    return Kband;
    }

   // FUNCTION: size_K
   // TASK: returns the number of columns/rows of K

  const unsigned & get_sizeK(void) const
    {
    return sizeK;
    }

  // FUNCTION: get_sizeK1
  // TASK: returns sizeK1 (useful only, if K is the kronecker product
  //       of two Penalty Matrizes K1 and K2

  const unsigned & get_sizeK1(void) const
    {
    return sizeK1;
    }


  // FUNCTION: get_sizeK2
  // TASK: returns sizeK2 (useful only, if K is the kronecker product
  //       of two Penalty Matrizes K1 and K2

  const unsigned & get_sizeK2(void) const
    {
    return sizeK2;
    }


  const int & get_posbeg(const unsigned & i) const
    {
    return posbeg[i];
    }


  const int & get_posend(const unsigned & i) const
    {
    return posend[i];
    }


  const unsigned & get_nrblocks(const unsigned & blocksize)
    {
    return matquant[blocksize-min];
    }

 const unsigned & get_minsize(void)
   {
   return min;
   }

 const unsigned & get_maxsize(void)
   {
   return max;
   }

 const fieldtype & get_type(void)
   {
   return type;
   }

 const unsigned & get_period(void)
   {
   return period;
   }

 ST::string get_typeasstring(void)
   {
   if (type==RW1)
     return "first order random walk";
   else if (type==RW2)
     return "second order random walk";
   else if (type==mrf)
     return "spatial Markov random field";
   else if (type==seasonal)
     return ("seasonal component");
   else if (type==mrflinear)
     return "2 dimensional first order random walk";
   else if (type==mrfkronecker)
     return "Kronecker product interaction";
   else return "";
   }


 ST::string get_typeshort(void)
   {
   if (type==RW1)
     return "rw1";
   else if (type==RW2)
     return "rw2";
   else if (type==mrf)
     return "spatial";
   else if (type==seasonal)
     return ("season");
   else if (type==mrflinear)
     return "2dimrw1";
   else if (type==mrfkronecker)
     return "kronecker";
   else return "";
   }

  bool get_polex(void)
    {
    return polex;
    }

 const statmatrix<int> & get_index(void)
   {
   return index;
   }

 const vector<double> & get_weight(void)
   {
   return weight;
   }

 const vector<ST::string> & get_values(void)
   {
   return effectvalues;
   }

 const vector<double> & get_valuesd(void)
   {
   return effectvdouble;
   }

 const ST::string & get_modname(void)
   {
   return modname;
   }

 const vector<datamatrix> & get_fc_random(void)
   {
   return fc_random;
   }

  unsigned get_category(const double & v) const;

  const vector<ST::string> & get_errormessages(void)
    {
    return errormessages;
    }

  // DESTRUCTOR

  ~PenaltyMatrix() {}

  };


//------------------------------------------------------------------------------
//---------------------------- class: FULLCOND_nonp -----------------------------
//------------------------------------------------------------------------------


enum centermethod {total,rowwise};


class __EXPORT_TYPE FULLCOND_nonp : public FULLCOND
  {


  protected:

  FULLCOND_const * fcconst;

  PenaltyMatrix * Pmatrix;           // Pointer to PenaltyMatrix

  DISTRIBUTION * likep;

  double sigma2;                     // Varianze parameter /tau^2 in the paper
  double lambda;

  datamatrix effectmod;             // Effectmodifier for varying coefficients
                                    // models
  bool varcoeff;                    // true for varying coefficient model,
                                    // otherwise false


  centermethod centerm;

  datamatrix accept;

  vector<ST::string> effectvalues;
  vector<double> effectvdouble;

  ST::string mapname;

  // FUNCTION: centerbeta
  // TASK: centers the current parametermatrix about the mean

  double centerbeta (vector<double> & weight);

  // FUNCTION: centerbeta2

  void centerbeta2 (void);

  public:

  // DEFAULT CONSTRUCTOR

  FULLCOND_nonp(void) : FULLCOND()
    {
    }

  // CONSTRUCTOR 1  (for additive models)
  // o    : pointer to MCMCoptions object
  // dp   : pointer to distribution object
  // K    : pointer to PenaltyMatrix object
  // minb : Minimum blocksize (minblock)
  // maxb : Maximum blocksize (maxblock)
  // fp   : file where sampled parameters are stored
  // pres : file where results are stored

  FULLCOND_nonp(MCMCoptions * o,DISTRIBUTION * dp,PenaltyMatrix * K,
                FULLCOND_const * fcco,
                const double & l,
                const ST::string & fp, const ST::string & pres,
                const ST::string & t,const ST::string & mn,
                const unsigned & c=0);


  // CONSTRUCTOR 2 (for varying coefficients models)
  // o      : pointer to MCMCoptions object
  // dp     : pointer to distribution object
  // K      : pointer to PenaltyMatrix object
  // effmod : Effectmodifier
  // minb   : Minimum blocksize (minblock)
  // maxb   : Maximum blocksize (maxblock)
  // fp     : file where sampled parameters are stored
  // pres   : file where results are stored

  FULLCOND_nonp(MCMCoptions * o,DISTRIBUTION * dp,PenaltyMatrix * K,
                FULLCOND_const * fcco,
                const double & l,
                const datamatrix & effmod,
                const ST::string & ti, const ST::string & fp,
                const ST::string & pres, const ST::string & mn,
                const unsigned & c=0);

  // COPY CONSTRUCTOR

  FULLCOND_nonp(const FULLCOND_nonp & fc);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_nonp & operator=(const FULLCOND_nonp & fc);

  void update(void);

  bool posteriormode(void)
    {
    return true;
    }

  void outresults(void);

  void get_effectmatrix(datamatrix & e,vector<ST::string> & enames,
                        unsigned be, unsigned en,effecttype t);

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND::reset();
    accept = datamatrix(beta.rows(),beta.cols(),0);
    }


  // FUNCTION: get_acceptance
  // TASK: stores the acceptance rates in file 'filename'

  void get_acceptance(const ST::string & filename);

  // FUNCTION: get_acceptance
  // TASK: writes the acceptance rates to filestream 'out'

  void get_acceptance(ostream & out);


  // FUNCTION: set_center

  void set_center(const centermethod & m)
    {
    centerm = m;
    }

  void set_effectmod(datamatrix v)
    {
    effectmod = v;
    }

  const vector<ST::string> & get_effectvalues(void)
    {
    return effectvalues;
    }

  // FUNCTION: update_sigma2
  // TASK: updates sigma2

  void update_sigma2(const double & s)
    {
    sigma2 = s;
    }

  // FUNCTION: compute_quadform
  // TASK: computes beta' K beta

  double compute_quadform(void);

  // FUNCTION: get_rankK
  // TASK: returns the rank of the penalty matrix K

  const unsigned & get_rankK(void)
    {
    return Pmatrix->rank_K();
    }

  const double & getlambda(void)
    {
    return lambda;
    }

  ST::string getinfo(void)
    {
    if (Pmatrix->get_type() == mrf)
      return mapname;
    else
      return title;
    }


  void init_name(const ST::string & na);

  void init_names(const vector<ST::string> & na);

  void init_priorassumptions(const ST::string & na);


  // DESTRUCTOR

  ~FULLCOND_nonp() {}

  };


}   // end: namespace MCMC

#endif
