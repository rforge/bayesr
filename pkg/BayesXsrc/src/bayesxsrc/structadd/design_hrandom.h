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



#if !defined (DESIGNhrandomINCLUDED)

#define DESIGNhrandomINCLUDED

#include"../export_type.h"
#include"statmat.h"
#include"sparsemat.h"

#include"Random.h"
#include"envmatrix_penalty.h"
#include"../values.h"
#include<fstream>
#include<vector>
#include<bitset>
#include"GENERAL_OPTIONS.h"
#include"clstring.h"
#include"distr.h"
#include"design.h"
#include<cmath>


namespace MCMC
{


//------------------------------------------------------------------------------
//--------------------------- CLASS: DESIGN_hrandom ----------------------------
//------------------------------------------------------------------------------

class __EXPORT_TYPE DESIGN_hrandom : public DESIGN
  {

  // beta contains linpredRE+random effect


  protected:

  DISTR * likep_RE;

  public:


//----------------------- CONSTRUCTORS, DESTRUCTOR -----------------------------

  // DEFAULT CONSTRUCTOR

  DESIGN_hrandom(void);

  // CONSTRUCTOR 1
  // Spatial covariates

  DESIGN_hrandom(const datamatrix & dm, const datamatrix & iv,
             GENERAL_OPTIONS * o,DISTR * dp,FC_linear * fcl, DISTR * dp_RE,
             vector<ST::string> & op,
             vector<ST::string> & vn);

  // COPY CONSTRUCTOR

  DESIGN_hrandom(const DESIGN_hrandom & m);

  // OVERLOADED ASSIGNMENT OPERATOR

  const DESIGN_hrandom & operator=(const DESIGN_hrandom & m);

  // virtual functions

  void init_data(const datamatrix & dm,const datamatrix & iv);

  void compute_penalty(void);

  void compute_penalty2(const datamatrix & pen);

  void compute_XtransposedWres(datamatrix & partres, double l);

  void compute_precision(double l);

  void compute_meaneffect(DISTR * level1_likep,double & meaneffect,
                          datamatrix & beta,datamatrix & meaneffectbeta,
                          bool computemeaneffect, double meaneffectconstant);

  void read_options(vector<ST::string> & op,vector<ST::string> & vn);

  void compute_basisNull(void);

  void outoptions(GENERAL_OPTIONS * op);

  void outbasis_R(ofstream & out);

  // DESTRUCTOR

  ~DESIGN_hrandom() {}

  };



} // end: namespace MCMC

#endif


