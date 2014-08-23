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



#if !defined (DISTRcategorical_mult_INCLUDED)
#define DISTRcategorical_mult_INCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"Random.h"
#include"GENERAL_OPTIONS.h"
#include"FC.h"
#include"distr.h"


namespace MCMC
{




//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_multinomprobit ----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_multinomprobit : public DISTR
  {

  protected:

  unsigned catnr;                              // category number

  bool master;                                 // master equation that updates
                                               // utilities yes/no

  vector<DISTR*> othercat;
  unsigned nrcat;                              // total number of categories
                                               // including master
  unsigned nrothercat;                         // number of servant categories

  double maxutility(vector<datamatrix*>
                    responsep,
                    const unsigned & i,
                    const unsigned & cat);


  public:

  datamatrix responsecat;

  void create_responsecat(void);

  void assign_othercat(DISTR* o);
  void assign_distributions(vector<DISTR * > dp);


   // DEFAULT CONSTRUCTOR

   DISTR_multinomprobit(void) : DISTR()
     {
     }

   // CONSTRUCTOR

   DISTR_multinomprobit(GENERAL_OPTIONS * o, unsigned & cnr, bool mast,
                        const datamatrix & r,
                        const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_multinomprobit(const DISTR_multinomprobit & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_multinomprobit & operator=(const DISTR_multinomprobit & nd);

   // DESTRUCTOR

   ~DISTR_multinomprobit() {}

  void compute_mu(const double * linpred,double * mu);

  void compute_deviance(const double * response, const double * weight,
                        const double * mu,double * deviance, double * scale)
                        const;

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void outoptions(void);

  void update(void);

  };


//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_multgaussian ------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_multgaussian : public DISTR_multinomprobit
  {

  protected:

  double devianceconst;


  void outresults_help(ST::string t,datamatrix & r);
  void initpointer(unsigned j,double* & worklin, double* & workresp);

  ST::string pathscale;
  FC FC_scale;
  FC FC_corr;

  double A;                    // Hyperparameter for Inverse Wishart
  datamatrix B;                // Hyperparameter for Inverse Wishart
  datamatrix sumB;
  datamatrix diff;

  datamatrix SIGMA_mr;         // Sigma_-r
  datamatrix SIGMA_rmr;        // Sigma_r|-r (row vector)

                               // helpmat1 stores SIGMA_rmr * SIGMA_mr

  datamatrix offset;           // stores the offsets in a nrobs x nrcat matrix
                                // i -th row corresponds to i-th observation
                                // r - th column stores o_r|-r 's


   // FUNCTION: compute_IWproduct
   // TASK: computes 0.5*SUM_{i=1}^{n} (y_i-eta_i)(y_i-eta_i)' and stores
   //       the result in sumB
   //       y_i = (y_i1,...y_ik)'

   void compute_IWproduct(void);

   // FUNCTION: compute_SIGMA_mr
   // TASK: computes the submatrix of Sigma where the elements corresponding to
   //       the r-th component of Sigma are removed
   //       stores the INVERSE of SIGMA_mr

   void compute_SIGMA_mr(unsigned r);

   // FUNCTION: compute_SIGMA_rmr
   // TASK: computes SIGMA_r|-r as a row vector

   void compute_SIGMA_rmr(unsigned r);

  // FUNCTION: compute_sigmarmr
  // TASK: computes for r=1,...,k: sigma2_r|-r =
  //       sigma2_r - SIGMA_r|-r * (SIGMA_-r)^-1 * (SIGMA_r|-r)'
  //       stores the results in the column vector sigma_rmr

  void compute_sigmarmr(void);

  // FUNCTION: compute_offset
  // TASK computes the offset and stores it in the matrix 'offset'

  void compute_offset(void);


  public:

   // DEFAULT CONSTRUCTOR

   DISTR_multgaussian(void) : DISTR_multinomprobit()
     {
     }

   // CONSTRUCTOR

   DISTR_multgaussian(const double & a,const double & b, unsigned & cnr,
                      GENERAL_OPTIONS * o, bool mast,
                      const datamatrix & r,const ST::string & ps,
                      const datamatrix & w=datamatrix());

   // COPY CONSTRUCTOR

   DISTR_multgaussian(const DISTR_multgaussian & nd);

   // OVERLOADED ASSIGNMENT OPERATOR

   const DISTR_multgaussian & operator=(const DISTR_multgaussian & nd);

   // DESTRUCTOR

   ~DISTR_multgaussian() {}

  void compute_mu(const double * linpred,double * mu);

  void compute_mu_mult(vector<double *> linpred,double * mu);


  void compute_deviance_mult(vector<double *> response,
                             vector<double *> weight,
                             vector<double *> linpred,
                             double * deviance,
                             vector<datamatrix *> aux);

  datamatrix * get_auxiliary_parameter(auxiliarytype t = auxcurrent);

  double loglikelihood(double * response, double * linpred,
                       double * weight);

  double loglikelihood_weightsone(double * response, double * linpred);

  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void outoptions(void);

  bool posteriormode(void);

  void update(void);

  void outresults(ofstream & out_stata,ofstream & out_R,ST::string pathresults="");

  void get_samples(const ST::string & filename,ofstream & outg) const;

  };




//------------------------------------------------------------------------------
//-------------------- CLASS: DISTRIBUTION_multinomlogit -----------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DISTR_multinomlogit : public DISTR_multinomprobit
{
 protected:

	int H;
  datamatrix SQ;
  datamatrix weights_mixed;


 public:

 	// DEFAULT CONSTRUCTOR

 	DISTR_multinomlogit(void) : DISTR_multinomprobit()
 		{
 		}

 	// CONSTRUCTOR1
 	DISTR_multinomlogit(GENERAL_OPTIONS * o,unsigned & cnr,
                        bool mast,
                        const datamatrix & r,
                        const datamatrix & w=datamatrix());


 	// COPY CONSTRUCTOR
 	DISTR_multinomlogit(const DISTR_multinomlogit & nd);


 	// OVERLOADED ASSIGNMENT OPERATOR
 	const DISTR_multinomlogit & operator=(const DISTR_multinomlogit & nd);


 	// DESTRUCTOR
 	~DISTR_multinomlogit()
 	{
 	}
////////////////////

/*
 	double compute_MSE();
*/

//  basis class implementation
// 	void compute_mu(const double * linpred,double * mu);

// basis class implementation
//  void compute_deviance(const double * response, const double * weight,
//                        const double * mu,double * deviance,
//                        double * deviancesat, double * scale) const;

// basis class implementation
//  double loglikelihood(double * response, double * linpred,
//                       double * weight) const;


// basis class implementation
//  double loglikelihood_weightsone(double * response, double * linpred) const;


// basis class implementation
/*
  double compute_iwls(double * response, double * linpred,
                      double * weight, double * workingweight,
                      double * workingresponse, const bool & like);

  void compute_iwls_wweightschange_weightsone(
                                         double * response, double * linpred,
                                         double * workingweight,
                                         double * workingresponse,double & like,
                                         const bool & compute_like);

  void compute_iwls_wweightsnochange_constant(double * response,
                                              double * linpred,
                                              double * workingweight,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

  void compute_iwls_wweightsnochange_one(double * response,
                                              double * linpred,
                                              double * workingresponse,
                                              double & like,
                                              const bool & compute_like);

*/


	void outoptions();

	// FUNCTION: update
	// TASK: uptdates the scale parameter

	void update(void);



  // no results
  // void outresults();

  // not required
	// double get_scalemean(void);

//  basis class implementation
// 	void sample_responses(unsigned i,datamatrix & sr);

//  basis class implementation
//	void sample_responses_cv();

//  basis class implementation
//	void outresults_predictive_check();

};




} // end: namespace MCMC


#endif
