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





#include "fullcond_dag.h"

#include <iterator>


namespace MCMC
{

	FULLCOND_dag::FULLCOND_dag (double s_i, unsigned int num,
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp)
							: FULLCOND(o,d,t,r,c,fp)
	{

		 nvar = d.cols();
		 nobs = d.rows();

		 adcol= Matrix<unsigned>(nvar,1,0);

		 ncoef_m = num_pa()+1;
		 ncoef_ia = 0;
		 ncoef = ncoef_m+ncoef_ia;
		 self=num;

		 lin = datamatrix (nobs, 1);
		 lin_prop = datamatrix (nobs, 1);

		 beta_help = datamatrix(ncoef,1,1);
		 x = datamatrix(nobs,ncoef,1); //matrix of covariates; first column is 1
		 xx = datamatrix(ncoef,ncoef,1);
		 Sigma = datamatrix(ncoef,ncoef,0);
		 alpha=nvar;

		 create_matrices();

		 priori_sigma = "non_inf";
		 priori_beta="inf";

		 if(priori_sigma=="non_inf")
		 {
			 a_invg = 1;
			 b_invg = 0.005;
		 }
		 else if (priori_sigma=="inf")
		 {
			a_invg = nvar/2;  //(nvar+1)/2;
			b_invg = 0.5;
		 }
		 sigma_i = s_i;
		 sigma_prop =1;
		 
		 print_dags = false;
 }


	// constructor 2
	FULLCOND_dag::FULLCOND_dag (double value_a, double value_b, ST::string prio_sig, 
							bool dags_all, const datamatrix & res, 
							double s_i, unsigned int num,
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp)
							: FULLCOND(o,d,t,r,c,fp)
	{
		y = res;

		adcol=Matrix<unsigned>(r,1,0);

		 nvar = d.cols();
		 nobs = d.rows();
		 ncoef_m = num_pa()+1;
		 ncoef_ia = 0;
		 ncoef = ncoef_m+ncoef_ia;
		 self=num;

		 lin = datamatrix (nobs, 1,0);
		 lin_prop = datamatrix (nobs, 1,0);

		 x = datamatrix(nobs,ncoef,0); //matrix of covariates; first column is 1
		 xx = datamatrix(ncoef,ncoef,0);
         beta_help = datamatrix(ncoef,1,1);
		 Sigma = datamatrix(ncoef,ncoef,0);
		 
		 alpha=nvar;

		 create_matrices();

		 priori_sigma= prio_sig;
		 priori_beta= "non_inf";
		 
		 if(priori_sigma=="non_inf")
		 {
			 a_invg = 1;
			 b_invg = 0.005;
		 }
		 else if (priori_sigma=="inf")
		 {
			a_invg = (nvar+1)/2;
			b_invg = 0.5;
		 }
		 else
		 {
			 a_invg = value_a;
			 b_invg = value_b;
		 }

		 sigma_i = s_i;

		 sigma_prop=1;

		 print_dags = dags_all;

		 if(print_dags==true)
				flags[0]=0;


 }




  ST::string FULLCOND_dag::get_family(void)
  {
	  return "Gaussian";
  }




void FULLCOND_dag::calc_xx(void)
{
	unsigned i,j,k;

	double value;

	double * workxx;
	double * workx1;
	double * workx2; 

	workxx = xx.getV();

	for(i=0; i<ncoef; i++)
	{
		for(j=0; j<ncoef; j++, workxx++)
		{
			workx1 = x.getV()+i;
			workx2 = x.getV()+j;
			value =0;

			for(k=0; k<nobs; k++)
			{
				value = value + (*workx1) * (*workx2);
				workx1 = workx1 + ncoef;
				workx2 = workx2 + ncoef;
			}
			*workxx = value;
		}	
	}
}




void FULLCOND_dag::calc_Sigma (void)
{
	

	if(priori_beta == "inf")
	{
	   unsigned int i;
	   Sigma.assign(xx);

	  for(i=0; i<ncoef; i++)
		  Sigma(i,i) = Sigma(i,i)+alpha;

	  Sigma.assign(Sigma.cinverse());
	  Sigma.assign(Sigma * sigma_i );
	}
	else
		Sigma.assign(xx.cinverse());
	
}







	// FUNCTION initialize
	// TASK: initializes x and xx  and y for pred_mod[i] (regression model i)
	void FULLCOND_dag::initialize (const adja & zeta, unsigned int j)
	{
		unsigned k, i;

		y = data.getCol(j);




		for(k=0; k<nvar; k++)	
			adcol(k,0) = unsigned(zeta(k,j));
		
		ncoef_m=1;
		for(i=0; i<nvar; i++)
		{
			if(zeta(i,j) ==1)
				ncoef_m++;
		}

		ncoef_ia = 0;

		initialize_ia(zeta,j); 
		
		ncoef=ncoef_ia+ncoef_m;


		//setbeta
		x = datamatrix (nobs, ncoef,1);
		xx = datamatrix(ncoef, ncoef);
		Sigma = datamatrix(ncoef,ncoef,0);
		beta_mean = datamatrix(ncoef,1,1);
		
		write_to_x(zeta);	// writes the data of the main effects to x

		write_ia_to_x();

		calc_xx();

		beta_help = datamatrix(ncoef,1,1);

		// initializing of beta_help
		unsigned m;
		for(m=0; m<ncoef; m++)
			beta_help(m,0) = rand_normal();


		// corresponding "initializing" of beta
		write_to_beta();

		calc_lin();

		SQT_x = calc_SQT_x();
		SQT_b = calc_SQT_b();
		SQT_x_n = 0;
		SQT_b_n = 0;
		calc_Sigma ();
		calc_beta_mean ();

		create_matrices();
	}



  


  void FULLCOND_dag::create_matrices (void)
  {
	  // parameters for a birth-step (one dimension more)
	  // if(ncoef-ncoef_ia<nvar)
	  {
		  b_new_b = datamatrix(ncoef+1,1,0);
		  x_new_b  = datamatrix(nobs,ncoef+1,0);
		  xx_new_b = datamatrix(ncoef+1,ncoef+1,0);
	  }

	  // parameters for a death-step (one dimension less)
	  if(ncoef>1)
	  {
		  b_new_d  = datamatrix(ncoef-1,1,0);
		  x_new_d  = datamatrix(nobs,ncoef-1,0);
		  xx_new_d = datamatrix(ncoef-1,ncoef-1,0);
	  }
  }



  void FULLCOND_dag::create_matrices (ST::string step, unsigned ncoef_new)
  {
	  if(step=="b")
	  {
		  b_new_b = datamatrix(ncoef_new,1,0);
		  x_new_b  = datamatrix(nobs,ncoef_new,0);
		  xx_new_b = datamatrix(ncoef_new,ncoef_new,0);
	  }
	  // parameters for a death-step (one dimension less)
	  else if(step=="d")
	  {
		  b_new_d  = datamatrix(ncoef_new,1,0);
		  x_new_d  = datamatrix(nobs,ncoef_new,0);
		  xx_new_d = datamatrix(ncoef_new,ncoef_new,0);
	  }
  }



  // FUNCTION calc_kq_est
  // TASK: sets beta_n equal to the kq-estimator for given x, xx  
  // where beta_n denotes ALL  regression coefficions of the models (old ones and added)
	void FULLCOND_dag::calc_kq_est (datamatrix & beta_n, const datamatrix & x_new, 
														 const datamatrix & xx_new)
	{
		unsigned coef_n=xx_new.cols();
		datamatrix sigma (coef_n, coef_n);
        
		sigma.assign(xx_new.cinverse());

		unsigned i,j,k;

		double sum1, sum2;

		double * worky; 
		double * workx;
		double * worksigma;
		double * workbeta_n;

		workbeta_n = beta_n.getV();

		for(k=0; k<coef_n; k++, workbeta_n++)
		{
			worksigma = sigma.getV() + k;
			
			sum1 = 0;
			for(i=0; i<coef_n; i++)
			{
				workx = x_new.getV()+i;
				sum2 = 0;
				worky = y.getV();

				for(j=0; j<nobs; j++, worky++)
				{
					sum2 = sum2 + (*workx) * (*worky);
					workx = workx + xx_new.cols();
				}

				sum1 = sum1 + (*worksigma) * sum2;
				worksigma = worksigma + coef_n;
			}
			
			*workbeta_n = sum1;
		}
	}



  // FUNCTION calc_kq_est
  // TASK: sets beta_n equal to the kq-estimator for given x_ia, xx_ia  
  // where beta_n denotes ONLY the regression coefficions of the added covariables
	void FULLCOND_dag::calc_kq_est (datamatrix & beta_n, const datamatrix & x_ia,
									const datamatrix & xx_ia, const datamatrix & y_ia)
	{
		unsigned i,j,k;

		unsigned coef_n = xx_ia.cols();

		datamatrix sigma (coef_n, coef_n);




		/*******************************
		for(i=0;i<coef_n; i++)
		{
			for(j=0; j<coef_n; j++)
				cout<<xx_ia(i,j)<<" ";

			cout<<endl;
		}
		cout<<endl<<endl;
		*******************************/

		
		sigma.assign(xx_ia.cinverse());
		
		double sum1, sum2;

		double * worky; 
		double * workx;
		double * worksigma;
		double * workbeta_n;

		workbeta_n = beta_n.getV();

		for(k=0; k<coef_n; k++, workbeta_n++)
		{
			worksigma = sigma.getV() + k;
			
			sum1 = 0;
			for(i=0; i<coef_n; i++)
			{
				workx = x_ia.getV()+i;
				sum2 = 0;
				worky = y_ia.getV();

				for(j=0; j<nobs; j++, worky++)
				{
					sum2 = sum2 + (*workx) * (*worky);
					workx = workx + coef_n;
				}

				sum1 = sum1 + (*worksigma) * sum2;
				worksigma = worksigma + coef_n;
			}
			
			*workbeta_n = sum1;
		}
	}



	void FULLCOND_dag::calc_beta_mean (void)
	{
		unsigned i,j,k;

		double sum1, sum2;

		double * worky; 
		double * workx;
		double * worksigma;
		double * workbeta_mean;
		
		workbeta_mean = beta_mean.getV();

		if(priori_beta == "non_inf")
		{
			// beta_mean = SIGMA * X'y; SIGMA = inv(X'X)

			for(k=0; k<ncoef; k++, workbeta_mean++)
			{
				worksigma = Sigma.getV() + k;
				
				sum1 = 0;
				for(i=0; i<ncoef; i++)
				{
					workx = x.getV()+i;
					sum2 = 0;
					worky = y.getV();
					for(j=0; j<nobs; j++, worky++)
					{
						sum2 = sum2 + (*workx) * (*worky);
						workx = workx + ncoef;
					}
					sum1 = sum1 + (*worksigma) * sum2;
					worksigma = worksigma + ncoef;
				}
				
				*workbeta_mean = sum1;
			}		
		}
		else if (priori_beta == "inf")
		{
			worksigma = Sigma.getV();
			for(k=0; k<ncoef; k++, workbeta_mean++)
			{
				sum1 = 0;
				for(i=0; i<ncoef; i++, worksigma++)
				{
					worky = y.getV();
					workx = x.getV()+i;
					sum2 = 0;	

					for(j=0; j<nobs; j++, worky++)
					{
						sum2 = sum2 + (*workx) * (*worky);  // + sigma_i/alpha * I * b_i
						workx = workx + ncoef;
					}

					sum1 = sum1 + (*worksigma) * sum2;
				}
				*workbeta_mean = sum1/sigma_i;
			} 
		}
	}





	double FULLCOND_dag::log_p_x(void) 
  {
	  return -1/(2*sigma_i) * get_SQT_x();
  }
	



double FULLCOND_dag::calc_yXb( const datamatrix & yy, const datamatrix & XX,
							  const datamatrix & bb) 
{
		// (yy-Xb)' (yy-Xb);


		double value, help;
		unsigned i, j;
		unsigned cols = XX.cols();
		unsigned rows = XX.rows();

		double * workX = XX.getV();
		double * workyy = yy.getV();
		double * workb;

		value=0;

		for(i=0; i<rows; i++, workyy++)
		{
			help=0;
			workb = bb.getV();
			for(j=0; j<cols; j++, workX++, workb++)
			{
				help = help + (*workX) * (*workb);
			}
			help = help - (*workyy);
			value = value + help*help;
		}

		return value;
}




double FULLCOND_dag::calc_yXb( const datamatrix & XX, const datamatrix & bb) 
{
		// (y-Xb)' (y-Xb);


		double value, help;
		unsigned i, j;
		unsigned cols = XX.cols();
		unsigned rows = XX.rows();

		double * workX = XX.getV();
		double * worky = y.getV();
		double * workb;

		value=0;

		for(i=0; i<rows; i++, worky++)
		{
			help=0;
			workb = bb.getV();
			for(j=0; j<cols; j++, workX++, workb++)
			{
				help = help + (*workX) * (*workb);
			}
			help = help - (*worky);
			value = value + help*help;
		}

		return value;
}




		
  double FULLCOND_dag::log_p_x(const datamatrix & b_new, const datamatrix & x_new) 
  {
	  // (y - x_new*b_new)' (y-x_new*b_new);

		unsigned i;
		double value, help;
		
		double * worky = y.getV();
		double * worklin = lin_prop.getV();

		value=0;

		for(i=0; i<nobs; i++, worky++, worklin++)
		{
			help = *worklin - (*worky);
			value = value + help*help;
		}

		return -1/(2*sigma_i) * value; 
  }





  double FULLCOND_dag::log_p_x(const datamatrix & b_new, const datamatrix & x_new, double sig_i_new) 
  {
	  // -(y-x_new*b_new)' (y-x_new*b_new) / 2sig;

		unsigned i;
		double value, help;
		
		double * worky = y.getV();
		double * worklin = lin_prop.getV();

		value=0;

		for(i=0; i<nobs; i++, worky++, worklin++)
		{
			help = *worklin - (*worky);
			value = value + help*help;
		}

	  return -1/(2*sig_i_new) * value; 
  }







  double FULLCOND_dag::log_p_b1(void)
  {
	  return -1/(2*sigma_i) *  get_SQT_b();
  }







double FULLCOND_dag::log_p_b1(const datamatrix & b_n)
{
	unsigned i;
	double value, help;
	double * workb_n;

	value=0;
	workb_n = b_n.getV();

	for(i=0; i<ncoef; i++, workb_n++)
	{
		help = *workb_n;
		value = value + help*help;
	}

	return -0.5 * value / sigma_i;
}







  double FULLCOND_dag::log_p_b2(void)
  {
	  return  -0.5*(ncoef*log(sigma_i) + get_SQT_b()/sigma_i);
  }

	  



   double FULLCOND_dag::log_p_b(const datamatrix & b_new)
   {
	   // value = b_new' * b_new
	   unsigned i;

	   double * workb_new;
	   workb_new = b_new.getV();

	   double value = 0;
	   double help; 

	   for(i=0; i<b_new.rows(); i++, workb_new++)
	   {
		   help = *workb_new;
		   value = value + help*help;
	   }

	   return -0.5* value / sigma_i ;
   }





	// FUNCTION: p_prop()
	// TARGET: returns the density of the proposal u, which is normaldistributed

	double FULLCOND_dag::p_prop(double prop)
	{
		double pi = 3.14159;

        double help =   2*pi*sigma_prop;

		return -1/2 * log(help) - prop*prop / (2*sigma_prop);
	}


   

	
	// FUNCTION: p_prop()
	// TARGET: returns the log-density of the proposal u, which is normaldistributed
	double FULLCOND_dag::p_prop(double prop, double mu, double sigma)
	{
		double pi = 3.14159;
		double help1, help2a, help2b, help3;

		help1 = (prop-mu)*(prop-mu);
		help2a = 2*pi*sigma;
        help2b = -1/2 * log(help2a);
        help3 = help2b - help1 / (2*sigma) ;

		return help3;
	}


	
	// FUNCTION: p_prop()
	// TARGET: returns the log-density of the proposal u, which is normaldistributed
	double FULLCOND_dag::p_prop(const datamatrix & prop, const datamatrix & mu, double sigma)
	{
		unsigned t; 
		double sum;
		

		double * workprop;
		double * workmu;

		sum=0;
		workprop = prop.getV();
		workmu = mu.getV();

		for(t=0; t<prop.rows(); t++, workprop++, workmu++)
		{
			sum = sum + p_prop(*workprop, *workmu, sigma);
		}
		return sum;
	}



	// FUNCTION: p_prop()
	// TARGET: returns the log-density of the proposal u, which is normaldistributed
	double FULLCOND_dag::p_prop(const datamatrix & prop, const datamatrix & mu, const datamatrix & concentration)
	{
		unsigned n = mu.rows();


		datamatrix vec (n,1);
		datamatrix help2 (1,n);
		datamatrix help3 (1,1);

		vec.minus(prop,mu);
		help2.mult(vec.transposed(), concentration);
		help3.mult(help2,vec);

		return help3(0,0);	
	}



	// FUNCTION: p_prop()
	// TARGET: returns the log-density of the proposal u, which is normaldistributed
	double FULLCOND_dag::b_distr(void)
	{


		datamatrix vec (ncoef,1);
		datamatrix help2 (1,ncoef);
		datamatrix help3 (1,1);

		vec.minus(beta_mean,beta_help);
		help2.mult(vec.transposed(), xx);
		help3.mult(help2,vec);

		return -0.5* ( Sigma.det() + help3(0,0));
	}








    // FUNCTION: p_prop()
	// TARGET: returns the log-density of beta_help, which is normaldistributed 
	// with mean beta_mean and variance sigma_i
	double FULLCOND_dag::p_prop(void)
	{
		return p_prop(beta_help, beta_mean, sigma_i);
	}







double FULLCOND_dag::calc_SQT_b(const datamatrix & b_n)
	{
		unsigned i;
		double value, help;
		double * workb_n;

		value=0;
		workb_n = b_n.getV();


		for(i=0; i<b_n.rows(); i++, workb_n++)
		{
			help = *workb_n;
			value = value + help*help;
		}

		return value;
	}












    double FULLCOND_dag::log_u(void)
	{
		unsigned i,j;
		double diff1, diff2, sum1, sum2, help1a, help1, help2, help2a;

		double * workxx; 
		double * workbeta_h1;
		double * workbeta_m1;


		workxx = xx.getV();

        help1a =   2*3.14159 ;
		help1 = -0.5*ncoef*log(help1a);


        help2a =  sigma_i *Sigma.det() ;
		help2 = - (ncoef * log(help2a));

		sum1 = 0;
		workbeta_h1 = beta_help.getV();
		workbeta_m1 = beta_mean.getV();

		for(i=0; i<ncoef; i++, workbeta_h1++, workbeta_m1++)
		{
			diff1 = (*workbeta_h1) - (*workbeta_m1);
			sum2 = 0;

			for(j=0; j<ncoef; j++, workxx++)
			{
				diff2 = beta_help(j,0)-beta_mean(j,0);
				sum2 = sum2 + diff2* (*workxx);
			}
			sum1 = sum1 + diff1*sum2;
		}

		return help1+help2-0.5*sigma_i*sum1;
	}










	 double FULLCOND_dag::log_u(const datamatrix & mean, const datamatrix & beta, 
								const datamatrix & Sigma, unsigned int ncoef)
	{
        double help1a = 2*3.14159;
		double help1 = -0.5*ncoef*log(help1a);
		double help2 = - (ncoef * log(Sigma.det()));
        
		datamatrix help3 (1,1,1);

		datamatrix help4 = beta-mean;
		
	//	help3.mult(help4.transposed(), Sigma.inverse()*help4));

		return help1+help2-0.5*help3;
	}








   void FULLCOND_dag::change(unsigned i, const datamatrix & beta_help_new, const datamatrix & x_new,
								const datamatrix & xx_new, unsigned int ncoef_new)
   {

	   lin.assign(lin_prop);
	   beta_help = beta_help_new;

	   x = x_new;	
	   xx = xx_new;

	   if(ncoef_new>ncoef)
			ncoef_m++;
	   else
		   ncoef_m--;




       //ncoef_ia is changed in change_current

	   ncoef = ncoef_new;
	   SQT_x = calc_SQT_x();					
	   SQT_b = calc_SQT_b();

	   Sigma = datamatrix (ncoef, ncoef,0); 
	   calc_Sigma ();
							
	   beta_mean = datamatrix (ncoef,1);
	   calc_beta_mean ();	

	   create_matrices();  //has to stand AFTER the changing of ncoef!!!
   }










	double FULLCOND_dag::new_comp_xx(unsigned int m, unsigned int n, unsigned int i, 
										unsigned int t)
	{
		assert( (t==m) || (t==n));
		assert(i<data.cols());

		double * workdata;
		double * workx;

		double sum = 0;
		unsigned k;

		if( m < t)
		{
			workdata = data.getV()+i;	
			workx = x.getV()+m;

			for(k=0; k<nobs; k++)
			{
				sum = sum + (*workx) * (*workdata);	//x(k,m)*data(k,i);
				workdata = workdata + nvar;
				workx = workx + ncoef;
			}
		}

		else if( m > t)
		{
			workdata = data.getV()+i;	
			workx = x.getV()+m-1;

			for(k=0; k<nobs; k++)
			{
				sum = sum +  (*workx) * (*workdata);	// x(k,m-1)*data(k,i);
				workdata = workdata + nvar;
				workx = workx + ncoef;
			}
		}

		else if( n < t) 
		{
			workdata = data.getV()+i;	
			workx = x.getV()+n;

			for(k=0; k<nobs; k++)
			{
				sum = sum + (*workx) * (*workdata);		//data(k,i)*x(k,n);
				workdata = workdata + nvar;
				workx = workx + ncoef;
			}
		}

		else if( n > t) 
		{ 
			workdata = data.getV()+i;	
			workx = x.getV()+n-1;

			for(k=0; k<nobs; k++)
			{
				sum = sum + (*workx) * (*workdata);		// data(k,i)*x(k,n-1);
				workdata = workdata + nvar;
				workx = workx + ncoef;
			}
		}

		else //m=n=t
		{		
			workdata = data.getV()+i;	
			for(k=0; k<nobs; k++)
			{
				sum = sum + (*workdata) * (*workdata);		// data(k,i)*data(k,i); 
				workdata = workdata + nvar;
			}
		}

		return sum;
	}










void FULLCOND_dag::calc_lin_prop(const datamatrix & x_prop, const datamatrix & b_prop)
	{
		unsigned i,j;

		double value; 
		double * workx = x_prop.getV();
		double * worklin = lin_prop.getV();
		double * workbeta_help;
		double intercept = *(b_prop.getV());

		if(b_prop.rows()==1) // so linpred=intercept 
		{
			for(i=0; i<nobs; i++, worklin++)
            {
				*worklin = intercept;
            }
		}
		else
		{
			for(i=0; i<nobs; i++, worklin++)  
			{
				value = intercept;				// intercept is everywhere the same 
				workbeta_help = b_prop.getV();
				workbeta_help++;
				workx++;

				for(j=1; j<x_prop.cols(); j++, workx++, workbeta_help++)
                {
					value = value + (*workx) * (*workbeta_help);
                }
				*worklin = value;
			}
		}
	}




	



	void FULLCOND_dag::calc_lin(void)
	{
		unsigned i,j;

		double value; 
		double * workx = x.getV();
		double * worklinp = lin.getV();
		double * workbeta_help;


		double intercept= beta_help(0,0);

		if(beta_help.rows()==1)  // so linpred=intercept 
		{
			for(i=0; i<nobs; i++, worklinp++)
				*worklinp = intercept;
		}
		else
		{
			for(i=0; i<nobs; i++, worklinp++)
			{
				value = intercept; // intercept is everywhere the same 
				workbeta_help = beta_help.getV();
				workbeta_help++;
				workx++;

				for(j=1; j<ncoef; j++, workx++, workbeta_help++)
					value = value + (*workx) * (*workbeta_help);

				*worklinp = value;
			}
		}

	}

		









	double FULLCOND_dag::calc_SQT_x(void)
	{	
		// (Xb-y)' (Xb-y)

		unsigned i;
		double value, help;
		
		double * worky = y.getV();
		double * worklinp = lin.getV();

		value=0;

		for(i=0; i<nobs; i++, worky++, worklinp++)
		{
			help = *worklinp - (*worky);
			value = value + help*help;
		}
		
		return value;
	}




	double FULLCOND_dag::calc_SQT_x(const datamatrix & X_n, const datamatrix & b_n)
	{	
		// (X_n*b_n -y)' (X_n*b_n-y)  calculates the SQT_x for the proposed X and b

		unsigned i;
		double value, help;

		//calculate the linear predictor for the proposed values
		calc_lin_prop(X_n, b_n);
		
		double * worky = y.getV();
		double * worklinp = lin_prop.getV();

		value=0;
        
		for(i=0; i<nobs; i++, worky++, worklinp++)
		{
            help = *worklinp - (*worky);
			value = value + help*help;
		}
		return value;
	}





	double FULLCOND_dag::calc_SQT_b(void)
	{
		unsigned i;
		double value, help;
		double * workbeta_help;

		value=0;
		workbeta_help = beta_help.getV();


		for(i=0; i<ncoef; i++, workbeta_help++)
		{
			help = *workbeta_help;
			value = value + help*help;
		}

		return value;
	}





	






	void FULLCOND_dag::update_sigma_i(void)
	{
		double a = a_invg + 0.5 * (nobs+ncoef);
		double b = b_invg + 0.5* (SQT_x + SQT_b);
		sigma_i = rand_invgamma(a, b);
	}







	void FULLCOND_dag::update(void)
	{

		beta_help.mult(Sigma.root(),rand_normvek(ncoef));
		beta_help.plus(beta_help, beta_mean);

		// corresponding updating of beta
		write_to_beta();
		
		calc_lin();

		SQT_x = calc_SQT_x();
		SQT_b = calc_SQT_b();

		update_sigma_i();


		
		if(print_dags==true)
			  FULLCOND::update();		// this command should be included (at the end of the
								// function) to update automatically the curent
								// mean, variance,minimum,maximum etc.

    }








	void FULLCOND_dag::write_to_beta(void)
	{
		if((optionsp->get_nriter() > optionsp->get_burnin()) &&
		   (optionsp->get_nriter() % (optionsp->get_step()) == 0))
		{
			unsigned m,l;

			double * workbeta = beta.getV(); // +1;
			double * workbeta_help = beta_help.getV(); //+1;

			*workbeta = *workbeta_help;

			workbeta++;
			workbeta_help++;

			l=1;
			for(m=0; m<nvar; m++, workbeta++)
			{
				if(adcol(m,0)==1)
				{
					*workbeta = *workbeta_help;
					 //beta(m,0) = beta_help(l,0); //
					workbeta_help ++;
					l++ ;
				}
				else 
				{
					if (m==self)
					{
						workbeta --;
						//m++;
					}
					else
					{
						*workbeta=0;
						//beta(m,0)=0;
					}
				}
			}
		}
	}







	void FULLCOND_dag::write_to_x(const adja & zeta)
	{
		unsigned m,k,l;

		double *workdata;
		double *workx;

		
		
		l=1; 
		for(m=0; m<nvar; m++)
		{
			if(zeta(m,self) ==1)
			{
				workdata=data.getV()+m;
				workx = x.getV()+l;

				for(k=0;k<nobs; k++)  //rows xx
				{
					*workx= *workdata; // AENDERUNG  x(k,l)= data(k,m);
					workdata = workdata+nvar;
					workx = workx+ncoef; 
				}

				l++;  // increase column of x
			}
		}
	}





	void FULLCOND_dag::write_to_y(const datamatrix & uti)
	{
		double * workuti =  uti.getV() + self;
		double * worky = y.getV();
		unsigned i; 

		for(i=0; i<nobs; i++, worky++)
		{
			*worky = *workuti;
			workuti = workuti + nvar;
		}
	}
			






	void FULLCOND_dag::write_to_x(const adja & zeta, const datamatrix & uti)
	{
		unsigned l=1; 
		unsigned m,k ;
		double * workuti;
		double * workx; 

		unsigned *workzeta = zeta.getV() + self;

		for(m=0; m<nvar; m++)
		{
			if(*workzeta==1)//if(zeta(m,self)==1)
			{
				workuti = uti.getV() + m;
				workx = x.getV() + l;
				for(k=0;k<nobs; k++)  //rows xx
				{
					*workx = *workuti;   //x(k,l) =uti(k,m)
					workuti = workuti + nvar;
					workx = workx + ncoef;
				}

				 l++; // increase column of x
			}
			workzeta = workzeta + nvar;
		}
	}










	void FULLCOND_dag::outresults(void)
    {
		
		if(print_dags==true)
		{
		FULLCOND::outresults(); // this command should be included at the beginning
                            // of the function, because the outresults function
                            // of the base class automatically computes
                            // estimated 10, 50, 90 percent quantiles of the
                            // parameters. In addition the acceptance rate will
                            // be computed and writen to (*optionsp->logout)
                            // as well as the Title/Name of the full conditional

    // It is very likely that the following variables from the base class will
    // be needed:

    // betamean             estimated posteriori mean
    // betavar              posteriori variance
    // betaqu10             posteriori 10% quantile
    // betaqu50             posteriori 50% quantile
    // betaqu90             posteriori 90% quantile
    // betamin              posteriori Minimum
    // betamax              posteriori Maximum


    // Write the posterior characteristics of beta either to a
    // file or to (*optionsp->logout)


		optionsp->out(" **********  CHARACTERISTICS OF REGRESSION MODEL "
		+ ST::inttostring(self)
		+" *************\n"),
		optionsp->out("\n");
		optionsp->out("\n");

		for(unsigned r=0; r<nvar; r++)
		{
			
			/*
			optionsp->out(ST::doubletostring(flags[0]) + " "
                        + ST::doubletostring(flags[1]) + " "
                        + ST::doubletostring(flags[2])+ "\n");
			*/


			if(r==0)
			{
				if(flags[0] == 1)
				{
					optionsp->out("Intercept: "
						+  ST::doubletostring(betamean(r,0),5) + "\n");
				}

				else if(flags[0] == 0)
				{
					optionsp->out("Intercept: \n");
					optionsp->out("\n");

                    ST::string l1 = ST::doubletostring(lower1,4);
                    ST::string l2 = ST::doubletostring(lower2,4);
                    ST::string u1 = ST::doubletostring(upper1,4);
                    ST::string u2 = ST::doubletostring(upper2,4);


					optionsp->out("mean: "		+  ST::doubletostring(betamean(r,0),5) + "\n");

					optionsp->out(l1 + "% quantile: " 
									 + ST::doubletostring(betaqu_l1_lower(r,0),5) + "\n");
					optionsp->out(l2 + "% quantile: " 
									 + ST::doubletostring(betaqu_l2_lower(r,0),5) + "\n");
					optionsp->out("50% quantile: " +  ST::doubletostring(betaqu50(r,0),5) + "\n");
					optionsp->out(u1 + "% quantile: " 
									 + ST::doubletostring(betaqu_l2_upper(r,0),5) + "\n");
					optionsp->out(u2 + "% quantile: " 
									 + ST::doubletostring(betaqu_l1_upper(r,0),5) + "\n");
					optionsp->out("\n");
				}
				optionsp->out("\n");
			}
			else if(r<self+1)
			{
				if(flags[0] == 1)
				{
					optionsp->out("regression coefficient of variable "
									+  ST::inttostring(r-1)
									+ " : "
									+ ST::doubletostring(betamean(r,0),5) + "\n");
				}

				else if(flags[0] == 0)
				{
					optionsp->out("regression coefficient of variable:"
								+  ST::inttostring(r-1)
								+ " : \n");
					optionsp->out("\n");

                    ST::string l1 = ST::doubletostring(lower1,4);
                    ST::string l2 = ST::doubletostring(lower2,4);
                    ST::string u1 = ST::doubletostring(upper1,4);
                    ST::string u2 = ST::doubletostring(upper2,4);


					optionsp->out("mean: "	   +  ST::doubletostring(betamean(r,0),5) + "\n");
					optionsp->out(l1 + "% quantile: " +
                    ST::doubletostring(betaqu_l1_lower(r,0),5) + "\n");
					optionsp->out(l2 + "% quantile: " +
                    ST::doubletostring(betaqu_l2_lower(r,0),5) + "\n");
					optionsp->out("50% quantile: " +  ST::doubletostring(betaqu50(r,0),5) + "\n");
					optionsp->out(u1 + "% quantile: " +
                    ST::doubletostring(betaqu_l2_upper(r,0),5) + "\n");
					optionsp->out(u2 + "% quantile: " +
                    ST::doubletostring(betaqu_l1_upper(r,0),5) + "\n");
					optionsp->out("\n");
				}
				optionsp->out("\n");
			}
			else if(r>=self+1)
			{
				if(flags[0] == 1)
				{

					optionsp->out("regression coefficient of variable:"
									+  ST::inttostring(r)
									+ " : "
									+  ST::doubletostring(betamean(r,0),5) + "\n");
				}

				if(flags[0] == 0)
				{
					optionsp->out("regression coefficient of variable "
								+  ST::inttostring(r)
								+ " : \n");
					optionsp->out("\n");

                    ST::string l1 = ST::doubletostring(lower1,4);
                    ST::string l2 = ST::doubletostring(lower2,4);
                    ST::string u1 = ST::doubletostring(upper1,4);
                    ST::string u2 = ST::doubletostring(upper2,4);

					optionsp->out("mean: "	   +  ST::doubletostring(betamean(r,0),5) + "\n");
					optionsp->out(l1 + "% quantile: " +
                    ST::doubletostring(betaqu_l1_lower(r,0),5) + "\n");
					optionsp->out(l2 + "% quantile: " +
                    ST::doubletostring(betaqu_l2_lower(r,0),5) + "\n");
					optionsp->out("50% quantile: " +  ST::doubletostring(betaqu50(r,0),5) + "\n");
					optionsp->out(u1 + "% quantile: " +
                    ST::doubletostring(betaqu_l2_upper(r,0),5) + "\n");
					optionsp->out(u2 + "% quantile: " +
                    ST::doubletostring(betaqu_l1_upper(r,0),5) + "\n");
					optionsp->out("\n");

				}
				optionsp->out("\n");
			}
		}
		optionsp->out("\n");
		optionsp->out("\n");

		}    // if (print_dags==true)




    }
	
	
	// FUNCTION: outoptions
    // TASK: writes estimation options (hyperparameters, etc.) to outputstream
	  void FULLCOND_dag::outoptions(void)
	  {
			if(self==0)
			{
				if(print_dags==true)
				  {
					  optionsp->out("Estimations of regression coefficents are listed. \n");
				}
		  }
	  }



	} // namespace MCMC
