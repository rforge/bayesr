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



#if !defined (TVARIANCE2DIM_INCLUDED)
#define TVARIANCE2DIM_INCLUDED

#include"../export_type.h"
#include "fullcond_pspline_surf_gaussian.h"
#include "fullcond_nonp_gaussian.h"
#include "distribution.h"

#include<cmath>
#include<limits>

namespace MCMC
{

class __EXPORT_TYPE FULLCOND_tvariance2dim : public FULLCOND
  {


  protected:

  DISTRIBUTION_gamma dgam;

  FULLCOND_pspline_surf_gaussian * Kp;  // pointer to psplines full conditional
  FULLCOND_nonp_gaussian * Kp_spat;  // pointer to psplines full conditional
  bool spatial;

  bool Laplace;
  statmatrix<int> indexmat;

  double * Kmatdiag;                    //
  double * Kmatupper;

  ST::string pathresults;               // path for results

  unsigned nu;                          // hyperparameter nu

  unsigned m;                           // number of parameters per row
                                        // (or column) in the psplines fc

  bool rowwise;
  datamatrix u;

  envmatdouble K11;
  double detalt;
  double detneu;

  unsigned nrrows;

  vector<double> deltapropvec;
  vector<double> rowvec;
  vector<double> colvec;
  vector<double> betakvec;


  public:


  // DEFAULT CONSTRUCTOR

  FULLCOND_tvariance2dim(void) : FULLCOND()
    {
    }


  // CONSTRUCTOR

  FULLCOND_tvariance2dim(MCMCoptions * o,FULLCOND_pspline_surf_gaussian * p,
                     unsigned & v,const ST::string & ti, const ST::string & fp,
                     const ST::string & pres,const unsigned & bs,const bool & rw = false);

  FULLCOND_tvariance2dim(MCMCoptions * o,FULLCOND_nonp_gaussian * p,
                     unsigned & v,const ST::string & ti, const ST::string & fp,
                     const ST::string & pres,const unsigned & bs,const bool & rw = false);


  // COPY CONSTRUCTOR

  FULLCOND_tvariance2dim(const FULLCOND_tvariance2dim & t);

  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_tvariance2dim & operator=(const FULLCOND_tvariance2dim & t);

  void update(void);

  void update_2dim(void);

  void update_spat(void);

  void update_spat_laplace(void);

  bool posteriormode(void)
    {
    return true;
    }

  void set_Laplace(void)
    {
    Kp_spat->set_deltadim(nrpar);
    Laplace = true;
    }

  void outresults(void);

  void outoptions(void);

  // FUNCTION: reset
  // TASK: resets all parameters

  void reset(void)
    {
    FULLCOND::reset();
    setbeta(nrpar,1,1);
    }

  // FUNCTION: predict (virtual)
  // TASK: predicts the mean for a new observation Xnew

  void predict(const datamatrix & newX, datamatrix & linpred)
    {

    }

  // DESTRUCTOR

  ~FULLCOND_tvariance2dim() {}

  }; // end: class FULLCOND_tvariance

template <class T>
class NRVec {
private:
	int nn;	// size of array. upper index is nn-1
	T *v;
public:
	NRVec();
	explicit NRVec(int n);		// Zero-based array
	NRVec(const T &a, int n);	//initialize to constant value
	NRVec(const T *a, int n);	// Initialize to array
	NRVec(const NRVec &rhs);	// Copy constructor
	NRVec & operator=(const NRVec &rhs);	//assignment
	NRVec & operator=(const T &a);	//assign a to every element
	inline T & operator[](const int i);	//i'th element
	inline const T & operator[](const int i) const;
	inline int size() const;
	~NRVec();
};

template <class T>
NRVec<T>::NRVec() : nn(0), v(0) {}

template <class T>
NRVec<T>::NRVec(int n) : nn(n), v(new T[n]) {}

template <class T>
NRVec<T>::NRVec(const T& a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = a;
}

template <class T>
NRVec<T>::NRVec(const T *a, int n) : nn(n), v(new T[n])
{
	for(int i=0; i<n; i++)
		v[i] = *a++;
}

template <class T>
NRVec<T>::NRVec(const NRVec<T> &rhs) : nn(rhs.nn), v(new T[nn])
{
	for(int i=0; i<nn; i++)
		v[i] = rhs[i];
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const NRVec<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
	if (this != &rhs)
	{
		if (nn != rhs.nn) {
			if (v != 0) delete [] (v);
			nn=rhs.nn;
			v= new T[nn];
		}
		for (int i=0; i<nn; i++)
			v[i]=rhs[i];
	}
	return *this;
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const T &a)	//assign a to every element
{
	for (int i=0; i<nn; i++)
		v[i]=a;
	return *this;
}

template <class T>
inline T & NRVec<T>::operator[](const int i)	//subscripting
{
	return v[i];
}

template <class T>
inline const T & NRVec<T>::operator[](const int i) const	//subscripting
{
	return v[i];
}

template <class T>
inline int NRVec<T>::size() const
{
	return nn;
}

template <class T>
NRVec<T>::~NRVec()
{
	if (v != 0)
		delete[] (v);
}

typedef NRVec<double> Vec_DP;
typedef const NRVec<double> Vec_I_DP;

double besselK(const double x, const double xnu);
void bessik(const double x, const double xnu, double &ri, double &rk, double &rip, double &rkp);
void beschb(const double x, double &gam1, double &gam2, double &gampl, double &gammi);
double chebev(const double a, const double b, Vec_I_DP &c, const int m, const double x);

double log_besselK(const double x, const double xnu);

} // end: namespace MCMC

#endif




