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





#include "fullcond_dag_d.h"


namespace MCMC
{

	// constructor 1

	FULLCOND_dag_d::FULLCOND_dag_d (double s_i, unsigned int num, 
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp)
							: FULLCOND_dag(s_i, num, o,d,t,r,c,fp)
	{
		ini_dag_d();
	}



	// constructor 2
	FULLCOND_dag_d::FULLCOND_dag_d (double v_a, double v_b, ST::string prio_s, 
							bool d_all, const datamatrix & res, 
							double s_i, unsigned int num, 
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp)
							: FULLCOND_dag(v_a, v_b, prio_s, d_all, res, s_i,num, 
											o,d,t,r,c,fp)
	{
		ini_dag_d();
	}

	


  ST::string FULLCOND_dag_d::get_family(void)
  {
	  return "binary without interactions";
  }



	void FULLCOND_dag_d::ini_dag_d(void)
	{
		unsigned i;

		double * worky_true;
		double * workdata;

		y_true = datamatrix(nobs,1,0);
		
		worky_true = y_true.getV();
		workdata = data.getV() + self;

		for(i=0; i<nobs; i++, worky_true++)
		{
			*worky_true = *workdata;
			workdata = workdata + nvar; 
		}
	}


void FULLCOND_dag_d::draw_utilities(void)
{
	unsigned i;

	 double * worky_true = y_true.getV(); 
	 double * worky = y.getV();
	 double * worklinp = lin.getV(); 

	 for(i=0; i<nobs; i++, worky_true++, worky++, worklinp++)
	 {
        if(*worky_true==1)
			*worky = *worklinp + truncnormal(-*worklinp, 20- *worklinp);
		else if(*worky_true==0)
			*worky = *worklinp + truncnormal(-20-*worklinp, - *worklinp);
	 }
}


	// FUNCTION update
  // TASK: updates parameters (i.e. matrix beta and beta_help)
  void FULLCOND_dag_d::update(void)
  {
	  draw_utilities ();

	  calc_beta_mean();

	  FULLCOND_dag::update();
  }





	} // namespace MCMC
