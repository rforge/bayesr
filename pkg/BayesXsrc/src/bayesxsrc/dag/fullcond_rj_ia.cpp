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





#include "fullcond_rj_ia.h"
#include <set>
#include <algorithm>



namespace MCMC
{


FULLCOND_rj_ia::FULLCOND_rj_ia (vector < FULLCOND_dag_d_ia* > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp)
				: FULLCOND_rj( o,d,t,r,c,fp)
{
	
	//nvar=c;
	//nobs=d.rows();
	//assert(c==r);

	change_preg_mods(dagp);

	ini_structure();
	//ini_ratio();
	//ini_hyperpar();

	//zetamean = datamatrix(nvar,nvar,0);
    //nrpar = nvar*nvar;
    //setzeta(zeta);

	//set_options();*/


}  //constructor_1








FULLCOND_rj_ia::FULLCOND_rj_ia (unsigned int l, double a , ST::string swi, ST::string pm, unsigned & ty,
						  vector < FULLCOND_dag * > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp)
				: FULLCOND_rj(l,a, swi, pm, ty, dagp, o,d,t,r,c,fp)
{
//	preg_mods = dagp;
	ini_structure(type);
	

}  //constructor_2






 // COPY CONSTRUCTOR

  FULLCOND_rj_ia::FULLCOND_rj_ia(const FULLCOND_rj_ia & fc) : FULLCOND_rj(FULLCOND_rj(fc))
  {
  }



  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_rj_ia& FULLCOND_rj_ia::operator=(const FULLCOND_rj_ia & fc)
  {
	  if (this==&fc)
		  return *this;

     FULLCOND_rj::operator=(FULLCOND_rj(fc));
	  
	  return *this;
  }






  // FUNCTION: change_preg_mods
  // TASK: changes preg_mods
  void FULLCOND_rj_ia::change_preg_mods( vector <FULLCOND_dag_d_ia* > dagp)
  {

	  for(unsigned i=0; i<nvar; i++)
	  {
		 preg_mods.push_back(dagp[i]);
	  }
  }



// FUNCTION: death_step
// makes death step, tries to delete edge i->j

void FULLCOND_rj_ia::death_step(unsigned int i, unsigned int j)
{

	// 1. Bestimme Anzahl der zu entferneneden Regressionskoeffizienten

	unsigned ncoef_old = preg_mods[j]->get_ncoef();
	unsigned num_ia_del = preg_mods[j]->ia_of_i(i);		 //number of to deletig interactions

	unsigned ncoef_new = ncoef_old - 1 - num_ia_del;
	
	if(num_ia_del==0)
	{
		// instead of: datamatrix b_new (ncoef_new,1);
		datamatrix & b_new = preg_mods[j]->get_b_new_d(); 	
		// instead of: datamatrix x_new (nobs,ncoef_new);
		datamatrix & x_new = preg_mods[j]->get_x_new_d();
		// instead of: datamatrix xx_new (ncoef_new,ncoef_new);
		datamatrix & xx_new = preg_mods[j]->get_xx_new_d();

		double beta_old; //coefficient which will vanish
		
		// computing of the new values
		make_new_d("d", i,j,xx_new, beta_old, b_new, x_new);
		


		// calculate ratio

		double ratio;
		double log_num;
		double log_denom;

		log_num = preg_mods[j]->calc_SQT_x(x_new, b_new) + preg_mods[j]->calc_SQT_b(b_new);
		log_denom = preg_mods[j]->calc_SQT_x() + preg_mods[j]->calc_SQT_b();

		ratio = -1/(2*preg_mods[j]->get_sigma_i()) 
				*(log_num - log_denom) + p_prop(beta_old) ;




		// accept proposal and change the corresponding values

		if(func::accept(ratio) == true)
		{
			zeta(i,j)=0;	
			zeta.change_list(i,j,1);	
			preg_mods[j]->change_adcol(i,0);
			preg_mods[j]->change(b_new, x_new, xx_new, ncoef_new);

			acceptance_d ++;
			nrtrials_d ++;
			zeta.edge_minus();

			cout<<"successful death"<<endl;
		}

		nrtrials_d ++;
		step_aborted = false;
	}
	else
	{
		datamatrix b_new (ncoef_new,1);
		datamatrix x_new (nobs,ncoef_new);
		datamatrix xx_new (ncoef_new,ncoef_new);

		// vector of interactions that will vanish
		vector <vector <unsigned > > ias_del;
		preg_mods[j]->ia_of_i( i, ias_del);


		// vector of coefficients which will vanish
		datamatrix beta_old (1+num_ia_del,1);
		
		// current_ia of model j, after successful death-step
		vector <vector <unsigned > > current_ia_n;

	
		// computing of the new values
		make_new_d_ia("d", i, j, num_ia_del, beta_old, current_ia_n, xx_new, b_new, x_new);




	// calculate ratio

  
     double ratio;
	 double log_num;
	double log_denom;

	log_num = preg_mods[j]->calc_SQT_x(x_new, b_new) + preg_mods[j]->calc_SQT_b(b_new);
	log_denom = preg_mods[j]->calc_SQT_x() + preg_mods[j]->calc_SQT_b();

	double help1 = -1/(2*preg_mods[j]->get_sigma_i()) ;
	double help2 = log_num;
	double help3 = log_denom;
	double help4 = p_prop(beta_old) ;


	ratio = -1/(2*preg_mods[j]->get_sigma_i()) 
			*(log_num - log_denom) + p_prop(beta_old) ;




	// accept proposal and change the corresponding values

	if(func::accept(ratio) == true)
	{
		zeta(i,j)=0;	
		zeta.change_list(i,j,1);	
		preg_mods[j]->change_adcol(i,0);
		preg_mods[j]->change(b_new, x_new, xx_new, ncoef_new);

		preg_mods[j]->change_occur('d', ias_del);
		preg_mods[j]->change_current('d', ias_del);


		acceptance_d ++;
		nrtrials_d ++;
		zeta.edge_minus();

		cout<<"successful death_ia"<<endl;
	}

	nrtrials_d ++;
	step_aborted = false;
	//*************/

			}
}







/*

// FUNCTION: switch_step
// makes switch step from j->i to i->j

void FULLCOND_rj::switch_step(unsigned int i, unsigned int j)
{

	zeta(j,i)=0;

	zeta.change_list(j,i,1);

	if(zeta.azy_test(i,j)==true)
	{
		zeta(j,i)=1;
		zeta.change_list(j,i,0); 


		//************ DEATH-step for edge j->i  ******************

		unsigned int ncoef_old_i = preg_mods[i]->get_ncoef();
		unsigned int ncoef_new_i = ncoef_old_i - 1;
	


		// instead of: datamatrix b_new_i (ncoef_new_i,1);
		datamatrix & b_new_i = preg_mods[i]->get_b_new_d();
	
		// instead of: datamatrix x_new_i (nobs,ncoef_new_i);
		datamatrix & x_new_i = preg_mods[i]->get_x_new_d();

		// instead of: datamatrix xx_new_i (ncoef_new_i,ncoef_new_i);
		datamatrix & xx_new_i = preg_mods[i]->get_xx_new_d();

	
		
		

		double beta_old_ij; //coefficient which will vanish
		double sigma_new_i; //coefficient which will vanish
		double sigma_new_j; //coefficient which will vanish
		
		// computing of the new values
		make_new_d("s", j,i,xx_new_i, beta_old_ij, b_new_i, x_new_i);

		datamatrix help_i (ncoef_new_i,ncoef_new_i);

		help_i.assign(xx_new_i.inverse());

		datamatrix mean_i (ncoef_new_i,1);
		mean_i.mult(help_i, x_new_i.transposed()*preg_mods[i]->get_y());

		sigma_new_i = sample_sigma('i', i, ncoef_new_i, mean_i, x_new_i);

		b_new_i.mult(help_i.root(),rand_normvek(ncoef_new_i));
		b_new_i.plus(sqrt(sigma_new_i)*b_new_i,mean_i);

		preg_mods[i]->calc_lin_prop(x_new_i, b_new_i);


	//************ BIRTH-step for edge i->j ******************

		unsigned int ncoef_old_j = preg_mods[j]->get_ncoef();
		unsigned int ncoef_new_j = ncoef_old_j + 1;

		// instead of: datamatrix b_new_j (ncoef_new_j,1);
		datamatrix & b_new_j = preg_mods[j]->get_b_new_b();
		
		// instead of: datamatrix x_new_j (nobs,ncoef_new_j);
		datamatrix & x_new_j = preg_mods[j]->get_x_new_b() ;

		// instead of: datamatrix xx_new_j (ncoef_new_j,ncoef_new_j);
		datamatrix & xx_new_j = preg_mods[j]->get_xx_new_b() ;

	
		
		

		double beta_new_ji=0;

	
		// computing of the new values
		make_new_b("s", i,j, beta_new_ji, xx_new_j,b_new_j, x_new_j);

		datamatrix help_j (ncoef_new_j,ncoef_new_j);
		help_j.assign(xx_new_j.inverse());

		datamatrix mean_j (ncoef_new_j,1);
		mean_j.mult(help_j, x_new_j.transposed()*preg_mods[j]->get_y());

		sigma_new_j = sample_sigma('j', j, ncoef_new_j, mean_j, x_new_j);

		b_new_j.mult(help_j.root(),rand_normvek(ncoef_new_j));
		b_new_j.plus(sqrt(sigma_new_j)*b_new_j,mean_j);

		preg_mods[j]->calc_lin_prop(x_new_j, b_new_j);


		// calculate ratio
		double ratio = 	ratio_s( i, j, b_new_i, b_new_j, x_new_i, x_new_j, mean_i, mean_j,
									 sigma_new_i, sigma_new_j);
			
		//ratio_s(i, j, ncoef_new_i, ncoef_new_j, b_new_i, b_new_j, x_new_i, x_new_j, beta_old_i, beta_new_j, mean_prop_i, mean_prop_j);

		if(func::accept(ratio) == true)
		{
			zeta(j,i)=0;
			zeta(i,j)=1;
			zeta.change_list(i,j,2);  //last argument=2, because switch-step

			preg_mods[i]->change_adcol(j,0);
			preg_mods[j]->change_adcol(i,1);

			preg_mods[i]->change(b_new_i, x_new_i, xx_new_i, ncoef_new_i);
			preg_mods[j]->change(b_new_j, x_new_j, xx_new_j, ncoef_new_j);

			acceptance_s ++;
		}

		step_aborted = false;
		
	} // end of "if(azy_test(i,j)==true)"
	
	else // azy_test(i,j)==false"
	{
		zeta(j,i)=1;
		zeta.change_list(j,i,0);
	}

	nrtrials_s ++;	
}







	// FUNCTION: make_new_b 
	// TASK: computes the new values for a birth-step
	void FULLCOND_rj::make_new_b (ST::string step, unsigned int i, unsigned int j, double beta_new, 
				datamatrix & xx_new, datamatrix & b_new, datamatrix & x_new)
	{
		unsigned int ncoef_new = preg_mods[j]->get_ncoef() +1;
		unsigned int k, kk, l;
		

		// Definition of t:
		// the variable i is the t-th regression coefficient of j
		// without counting the intercept as coefficient
		unsigned t=1;
		for(k=0; k<nvar; k++)
		{
			if( (k<i) && (zeta(k,j)==1))
				t++;
		}

		

	   //compute x_new
		double * workx_new;
		double * workx;
		double * workdata;

		workx_new = x_new.getV();
		workx = preg_mods[j]->getV_x();
		workdata = data.getV()+i;

		for(k=0; k<nobs; k++)
		{
			for(l=0; l<ncoef_new; l++, workx_new++, workx++)
			{
				if(l != t)
					*workx_new = *workx;
				else 
				{
					*workx_new = *workdata;
					workdata = workdata + nvar;
					workx--;
				}
			}
		}
		
		//compute xx_new

		double * workxx_new;
		double * workxx;
		double * workhelp;
		double * workx1;
		double * workx2; 

		double value;

		workxx_new = xx_new.getV();
		workxx = preg_mods[j]->getV_xx();

		for(k=0; k<ncoef_new; k++)
		{
			if(k != t)
			{
				for(l=0; l<ncoef_new; l++, workxx++, workxx_new++)
				{
					if(l != t)
						*workxx_new = *workxx;
					else
					{
						workx1 = x_new.getV() + k;
						workx2 = x_new.getV() +l;
						value =0;

						for(kk=0; kk<nobs; kk++)
						{
							value = value + (*workx1) * (*workx2);
							workx1 = workx1 + ncoef_new;
							workx2 = workx2 + ncoef_new;
						}
									
						*workxx_new = value;  

						workxx--;
					}
				}
			}
			else
			{
				for(l=0; l<ncoef_new; l++,workxx_new++)
				{
					workx1 = x_new.getV()+k;
					workx2 = x_new.getV()+l;
					value =0;

					for(kk=0; kk<nobs; kk++)
					{
						value = value + (*workx1) * (*workx2);
						workx1 = workx1 + ncoef_new;
						workx2 = workx2 + ncoef_new;
					}

					*workxx_new = value;
				}
			}
		}


		if(step != "s") //if the birth-step is NOT part of a switch step
		{
			//compute proposed beta
			double * workb_new;
			double * workbeta;

			workb_new = b_new.getV(); //b_new.getV();
			workbeta = preg_mods[j]->getV_beta_help();

			for(k=0; k<ncoef_new ;	workb_new++, workbeta++, k++)
			{
				if(k!= t)
					*workb_new = *workbeta;
				else 
				{
					*workb_new = beta_new;
					workbeta--;
				}
			}

			//compute proposed linear predictor
			preg_mods[j]->calc_lin_prop( x_new, b_new); 
		}
	}
*/





	
   // FUNCTION: make_new_d_ia
   // TASK: computes the new values for a death-step
   void FULLCOND_rj_ia::make_new_d_ia (ST::string step, unsigned i, unsigned j, 
						unsigned num_ia_del, datamatrix & beta_old, 
						vector <vector <unsigned > > & current_ia_n, 
						datamatrix & xx_new, datamatrix & b_new, datamatrix & x_new)
	{
	    unsigned  k,l;
		unsigned  ncoef=preg_mods[j]->get_ncoef();

		//void ia_of_i( vector <vector <unsigned > > & v);


		// Definition of ia
		// contains interactions of j-th regression model that contain i
		vector< vector< unsigned> > ia;
		preg_mods[j]->ia_of_i(i,ia);	


		

		// Definition of pos_of_del
		// contains positions of the main effect and the corresponding 
		// interactions in j-th regression model that are going to be deleted
		// (position of intercept is 0)
		vector <unsigned> pos_of_del;	
		preg_mods[j]->get_pos(i,pos_of_del);

	





/*****************
		unsigned t=0;
		cout<<"Haupteffekte: ";
		for(k=0; k<nvar; k++)
		{
			if( (preg_mods[j]->get_adcol())(k,0)==1)
				cout<<k<<"  ";	
		}

		cout<<endl<<"Interaktionen: ";

		for(k=0; k<preg_mods[j]->get_ncoef_ia(); k++)
		{
			cout<<preg_mods[j]->get_current_ia(k)[0]
				<<preg_mods[j]->get_current_ia(k)[1]<<"  ";	
		}

		cout<<endl; 

		

		cout<<"Position of del for " <<i<<": ";
		for(k=0; k<pos_of_del.size(); k++)
		{
			cout<<pos_of_del[k]<<" ";
		}
		cout<<endl;
		cout<<endl;
**************/


		//compute x_new
		double * workx_new;
		double * workx;

		workx_new = x_new.getV();
		workx = preg_mods[j]->getV_x();

		unsigned pos1, pos2;
		unsigned t1, t2;
		unsigned sign=1;

		

		for(k=0; k<nobs; k++)
		{
			t1=0;
			pos1 = pos_of_del[t1];

			for(l=0; l<ncoef; l++, workx_new++, workx++)
			{
				if(l != pos1)
					*workx_new = *workx;
				else 
				{
					workx_new--;
					t1++;
					pos1 = pos_of_del[t1];
				}
			}
		}




		//compute xx_new
		double * workxx_new;
		double * workxx;

		workxx_new = xx_new.getV();
		workxx = preg_mods[j]->getV_xx();

		
		t1=0;
		pos1 = pos_of_del[t1];

		for(k=0; k<ncoef; k++)
		{
			t2=0;
			pos2 = pos_of_del[t2];

			if(k != pos1)
			{
				for(l=0; l<ncoef; l++, workxx++, workxx_new++)
				{
					if(l != pos2)
						*workxx_new = *workxx;
					else
					{
						workxx_new--;
						t2++;
						pos2 = pos_of_del[t2];
					}
				}
			}
			else
			{
				workxx = workxx + ncoef;
				t1++;
				pos1 = pos_of_del[t1];
			}
		}



		//save the coefficients which are going to be deleted
		t1 = 0;
		pos1 = pos_of_del[t1];
		double help;

		for(k=0; k<ncoef; k++)
		{
			if(k==pos1)
			{
				help = preg_mods[j]->get_beta_help(k,0);
				beta_old(t1,0) = help;
				t1++;
				pos1 = pos_of_del[t1];	
			}
		}

		if(step != "s") //if the death-step is NOT part of a switch step
		{
			//compute proposed beta
			double * workb_new;
			double * workbeta;

			workb_new = b_new.getV(); 
			workbeta = preg_mods[j]->getV_beta_help();

			t1 = 0;
			pos1 = pos_of_del[t1];

			for(k=0; k<ncoef; workb_new++, workbeta++, k++)
			{
				if(k!= pos1)
					*workb_new = *workbeta;
				else 
				{
					workb_new--;
					t1++;
					pos1 = pos_of_del[t1];
				}
			}


			unsigned stop =1;

			preg_mods[j]->calc_lin_prop(x_new, b_new);


		}



	}




/*
	// FUNCTION: sample_sigma
	// TARGET: samples the new variance of the regression model i in the switch step
	
	double FULLCOND_rj::sample_sigma(char vertex, unsigned int i, unsigned int ncoef_new_i, 
							const datamatrix & mean_i, const datamatrix & x_new_i)
	{
		double value; 

		if(vertex== 'i')
		{
			value = preg_mods[i]->calc_yXb(preg_mods[i]->get_y(),x_new_i, mean_i);

			alpha_sig_i = 0.5* (nobs - ncoef_new_i);
			beta_sig_i = 0.5 * value;

			
			return rand_invgamma(alpha_sig_i ,beta_sig_i);
		}
		else if(vertex=='j')
		{
			value = preg_mods[i]->calc_yXb(preg_mods[i]->get_y(),x_new_i, mean_i);

			alpha_sig_j = 0.5* (nobs - ncoef_new_i);
			beta_sig_j = 0.5*value;

			return rand_invgamma(alpha_sig_j ,beta_sig_j);
		}
		else
		{
			double a, b;

			value = preg_mods[i]->calc_yXb(preg_mods[i]->get_y(),x_new_i, mean_i);

			a = 0.5* (nobs - ncoef_new_i);
			b= 0.5 * value; 

			return rand_invgamma(a, b);
		}
	}

*/



/**** gibt´s nicht mehr**** 

// FUNCTION: ratio_b()
// TARGET: calculates the ratio of a birth-step

double FULLCOND_rj::ratio_b(unsigned int j, double u,
						const datamatrix & b_new, const datamatrix & x_new)
{
	double log_num = preg_mods[j]->log_p_x(b_new, x_new) 
							+ preg_mods[j]->log_p_b(b_new);
		
	double log_denom = preg_mods[j]->log_p_x()
						+ preg_mods[j]->log_p_b1();

	return log_num - log_denom - p_prop(u) ;
}





// FUNCTION: ratio_d()
// TARGET: calculates the ratio of a death-step

double FULLCOND_rj::ratio_d(unsigned int j, double u,
							const datamatrix & b_new, const datamatrix & x_new)
{
	double log_num = preg_mods[j]->log_p_x(b_new, x_new) 
						+ preg_mods[j]->log_p_b(b_new);

	double log_denom = preg_mods[j]->log_p_x()
						+ preg_mods[j]->log_p_b1();

	return log_num - log_denom  + p_prop(u);   
}


**** gibt´s nicht mehr**** */









/*
 // FUNCTION: ratio_s
 // TARGET: computes acceptance ratio in the switch step
 double FULLCOND_rj::ratio_s(unsigned int i,unsigned int j, 
							const datamatrix & b_new_i, const datamatrix & b_new_j, 
							const datamatrix & x_new_i, const datamatrix & x_new_j,
							const datamatrix & mean_i, const datamatrix & mean_j,
							double sigma_new_i, double sigma_new_j)
{
	double sigma_old_i;
	double sigma_old_j;
	unsigned ncoef_i;
	unsigned ncoef_j;

	sigma_old_i = preg_mods[i]->get_sigma_i();
	sigma_old_j = preg_mods[j]->get_sigma_i();
	ncoef_i = preg_mods[i]->get_ncoef();
	ncoef_j = preg_mods[j]->get_ncoef(); 


	double log_num1 = -0.5* (nobs*log(sigma_new_j) 
						+ preg_mods[j]->calc_yXb(x_new_j, b_new_j) / sigma_new_j);
	
	double log_den1 = -0.5* (nobs*log(sigma_old_j) 
						+ preg_mods[j]->get_SQT_x()/ sigma_old_j);

	double log_num2 = -0.5* (nobs*log(sigma_new_i) 
						+ preg_mods[i]->calc_yXb(x_new_i, b_new_i) / sigma_new_i);

	double log_den2 = -0.5* (nobs*log(sigma_old_i) 
						+ preg_mods[i]->get_SQT_x()/ sigma_old_i);
	
	double log_num3 = -0.5*(b_new_j.cols()*log(sigma_new_j) 
						+ preg_mods[j]->calc_SQT_b(b_new_j)/sigma_new_j);

	double log_den3 = -0.5*(ncoef_j*log(sigma_old_j) 
						+ preg_mods[j]->get_SQT_b()/sigma_old_j);


	double log_num4 = -0.5*(b_new_i.cols()*log(sigma_new_i) 
						+ preg_mods[i]->calc_SQT_b(b_new_i)/sigma_new_i);

	double log_den4 = -0.5*(ncoef_i*log(sigma_old_i) 
						+ preg_mods[i]->get_SQT_b()/sigma_old_i);



	double log_prop_b = 0;
	double log_prop_s = 0;

	double help_i_old, help_i_new, help_j_old, help_j_new;
	double a_sig_i, b_sig_i, a_sig_j, b_sig_j; 


	a_sig_i = 0.5* (nobs - preg_mods[i]->get_ncoef());
	a_sig_j = 0.5* (nobs - preg_mods[j]->get_ncoef());

	double value;
	
	value = preg_mods[i]->get_SQT_x();	
	b_sig_i =0.5*value;

	value = preg_mods[j]->get_SQT_x();
	b_sig_j =0.5*value;

	help_i_old = a_sig_i*log(b_sig_i) - log_gamma1(a_sig_i)
					- (a_sig_i+1)*log (sigma_old_i) 
					- b_sig_i/sigma_old_i;

	help_i_new =  alpha_sig_i*log(beta_sig_i) - log_gamma1(alpha_sig_i)
					- (alpha_sig_i+1)*log (sigma_new_i) - beta_sig_i/sigma_new_i;

	help_j_old = a_sig_j*log(b_sig_j) - log_gamma1(a_sig_j)
					- (a_sig_j+1)*log (sigma_old_j) 
					- b_sig_j/sigma_old_j;
	
	help_j_new =  alpha_sig_j*log(beta_sig_j) - log_gamma1(alpha_sig_j)
					- (alpha_sig_j+1)*log (sigma_new_j) - beta_sig_j/sigma_new_j;

	log_prop_s = help_i_old - help_i_new + help_j_old- help_j_new ;

	double log_ratio = log_num1 - log_den1 
						+ log_num2 - log_den2
						+ log_num3 - log_den3
						+ log_num4 - log_den4
						+ log_prop_b
						+ log_prop_s;

	return log_ratio; 	
 }












// FUNCTION: log_gamma
// TASK: returns the logarithm of the gammafunction when value=0.5*(unsigned)

double FULLCOND_rj::log_gamma(double value) const
{
	assert( (int(value) % 1==0) || (int(value) % 1==0.5));

    if (value==0 || value ==1) 
		return 0;
	else if(value==0.5) 
		return log(3.14159);
    else 
		return log(value) + log_gamma(value-1);
}







// FUNCTION: log_gamma1
// TASK: returns the logarithm of the gammafunction when value=0.5*(nobs-k)

double FULLCOND_rj::log_gamma1(double x) const
{
	assert( (int(x) % 1==0) || (int(x) % 1==0.5));
	assert( x<nobs );

	double k;
	double sum=0;

	double start; 

	if(nobs % 2 ==0)
	{
		if( (nobs-int(2*x)) % 2==0)
			start=0.5*nobs-1;
		else
			start=0.5*(nobs+1)-1;
	}
	else
	{
		if( (nobs-int(2*x)) % 2==0)
			start=0.5*nobs-1;
		else
			start=0.5*(nobs-1);
	}

	for(k=start; k>x-1; k--)
	{
		sum= sum+log(k);
	}

	if( (nobs-int(2*x)) % 2==0)
		return gamma_nobs1 - sum;
	else
		return gamma_nobs2 - sum;
}




// FUNCTION: p_prop()
// TARGET: returns the density of the proposal u, which is normaldistributed

double FULLCOND_rj::p_prop(double prop)
{
	double pi = 3.14159;
    double help = 2*pi*sigma_prop;
	return -1/2 * log(help) - prop*prop / (2*sigma_prop);
}





// FUNCTION: accept
// TASK: returns true with probability ratio
//
//bool FULLCOND_rj::accept (double ratio)
//{
//	double u = uniform();
//		
//	if( log(u) > ratio )
//		return false;
//	else 
//		return true ;
//}
	







  // FUNCTION: outresults
  // TASK: writes estimation results to logout or into a file (after estimation)

	void FULLCOND_rj::outresults(void)
    {
		int i,j;
        int end;
		unsigned number_total=0;

//		FULLCOND::outresults();

		 #if defined(MICROSOFT_VISUAL)
		{
			number = __min(number, freq.size());
		}
		#else
		{
			number = min(number, freq.size());
		}
		#endif


		std::sort(freq.begin(), freq.end());

		for(i=freq.size()-1; i>=0; i--)
		{
			number_total = number_total+freq[i].freq;
		}



//		ST::string path;
//		path= "c:\\eva\\mod_freq.txt";
//		ofstream fout(path.strtochar());


		optionsp->out("******** DIFFERENT MODELS sorted by frequencies  ********\n");
						
		double rel_freq;
		double prob_cum = 0;

		double help1 = 0;
		double hilfe;
		


		if(print_models=="all")
		{
			for(j=freq.size()-1; j>=0; j--)
			{
				hilfe = double(freq[freq.size()-1-j].freq);
				help1 = help1+double(freq[freq.size()-1-j].freq)/double(number_total); 
			}
		}
		if(number==freq.size() && print_models!="all")
		{
			for(j=freq.size()-1; j>=0; j--)
			{
                hilfe = double(freq[freq.size()-1-j].freq);
				help1 = help1+double(freq[freq.size()-1-j].freq)/double(number_total); 
			}
		}
		else
		{
            end=  freq.size()-number ;
			for(j=freq.size()-1; j>=end; j--)
			{
				hilfe = double(freq[freq.size()-1-j].freq);
				help1 = help1+double(freq[freq.size()-1-j].freq)/double(number_total); 
			}
		}




		if(	(help1<1-alpha && print_models=="normal")
			|| print_models=="limit")
		{
			optionsp->out("******** first " + ST::inttostring(number) + " of the most important models ********\n");
			optionsp->out("\n");

            end = freq.size()-number;

			for(i=freq.size()-1; i>=end; i--)
			{
				optionsp->out(freq[i].model
					+ "	"
					+ ST::inttostring(freq[i].nedges)
					+ "	"
					+ ST::inttostring(freq[i].freq)
					+ "	"
					+ ST::doubletostring(double(freq[i].freq)/number_total,3)
					+ "\n");
			}
		}
		else if( print_models=="all")
		{
			optionsp->out("******** all models ********\n");
			optionsp->out("\n");

			for(i=freq.size()-1; i>=0; i--)
			{
				optionsp->out(freq[i].model
					+ "	"
					+ ST::inttostring(freq[i].nedges)
					+ "	"
					+ ST::inttostring(freq[i].freq)
					+ "	"
					+ ST::doubletostring(double(freq[i].freq)/number_total,3)
					+ "\n");
			}
		}


		else if(print_models=="normal"|| print_models=="prob")
		{

			optionsp->out(" ****  alltogether " + ST::doubletostring(100-alpha*100) + " % of posterior probability ********\n");
			optionsp->out("\n");

			for(i=freq.size()-1; i>=0; i--)
			{
				if(prob_cum < 1-alpha)
				{
					rel_freq = double(freq[i].freq)/double(number_total);
					optionsp->out(freq[i].model
						+ "	"
						+ ST::inttostring(freq[i].nedges)
						+ "	"
						+ ST::inttostring(freq[i].freq)
						+ "	"
						+ ST::doubletostring(rel_freq,3)
						+ "\n");
					prob_cum = prob_cum + rel_freq;
				}
			}
		}

        optionsp->out("\n");
        optionsp->out("\n");
        optionsp->out("\n");
		optionsp->out("******** averaged adjacency matrix: ********\n");
		optionsp->out("\n");

        ST::string help;

		for(i=0; i<nvar; i++)
		{
        help="";
			for(j=0; j<nvar; j++)
				help = help + ST::doubletostring(zetamean(i,j),2) + "	";

            optionsp->out(help + "\n");
			optionsp->out("\n");
		}


		optionsp->out("\n");
        optionsp->out("\n");
        optionsp->out("\n");
		optionsp->out("********  mean of skelets:  ********\n");
		optionsp->out("\n");

		for(i=0; i<nvar; i++)
		{
        help="";
			for(j=0; j<nvar; j++)
			{
				if(i<j)
					help = help + ST::doubletostring(zetamean(i,j)+zetamean(j,i),2) + "	";
				else
					help = help + "*	";
			}
			optionsp->out(help + "\n");
			optionsp->out("\n");
		}




		optionsp->out("\n");
		optionsp->out("\n");
		optionsp->out("\n");
		optionsp->out("******** correlation ********\n");
		optionsp->out("\n");


		datamatrix correlation (nvar,nvar);
		datamatrix partial_corr (nvar,nvar);

		correlation.assign(data.corr());
		partial_corr.assign(data.partial_var());

		
		for(i=0; i<nvar; i++)
		{
        help="";
			for(j=0; j<nvar; j++)
				help = help + ST::doubletostring(correlation(i,j),2) + "	";

            optionsp->out(help + "\n");
			optionsp->out("\n");
		}
		
		
		optionsp->out("\n");
		optionsp->out("\n");
		optionsp->out("******** partial correlation ********\n");
		optionsp->out("\n");



		for(i=0; i<nvar; i++)
		{
        help="";
			for(j=0; j<nvar; j++)
				help = help + ST::doubletostring(partial_corr(i,j),2) + "	";

            optionsp->out(help + "\n");
			optionsp->out("\n");
		}



		optionsp->out("\n");


		
		optionsp->out("acceptance ratio for a birth-step: " 
			+ ST::doubletostring(double(acceptance_b)/nrtrials,3) + "\n");
        optionsp->out("acceptance ratio for a death-step: "  
			+ ST::doubletostring(double(acceptance_d)/nrtrials,3) + "\n");
        optionsp->out("acceptance ratio for a switch-step: "
			+ ST::doubletostring(double(acceptance_s)/nrtrials,3) + "\n");

		optionsp->out("\n");
		
		optionsp->out("acceptance ratio if a birth-step has already been proposed: "
			+ ST::doubletostring(double(acceptance_b)/nrtrials_b,3) + "\n");
        optionsp->out("acceptance ratio if a death-step has already been proposed: "  
			+ ST::doubletostring(double(acceptance_d)/nrtrials_d,3) + "\n");
        optionsp->out("acceptance ratio if a switch-step has already been proposed: "
			+ ST::doubletostring(double(acceptance_s)/nrtrials_s,3) + "\n");

		optionsp->out("\n");

		optionsp->out("A birth-step has been proposed: "
			+ ST::inttostring(nrtrials_b) + " times.\n");
        optionsp->out("A death-step has been proposed: "
			+ ST::inttostring(nrtrials_d) + " times.\n");
        optionsp->out("A switch-step has been proposed: "
			+ ST::inttostring(nrtrials_s) + " times.\n");

		optionsp->out("\n");
		optionsp->out("\n");
}





	void FULLCOND_rj::update(void)
    {
		  rj_step ();

		  nrtrials++;


		  update_zeta();
		  
		  
		  //std::set < ST::string, unsigned int> freq;

		  if((optionsp->get_nriter() > optionsp->get_burnin()) &&
			  (optionsp->get_nriter() % (optionsp->get_step()) == 0))
		  {
			  store_model();
		  }

		 /* {
			  unsigned int r,c, count;
			  ST::string s;
			  bool found = false;

			  for(r=0; r<nvar; r++)
				  for(c=0, count=0; c<nvar; c++, count++)
				  {
					  if(((count+1)%nvar==0) && (count!=nvar*nvar))
						  s = s + ST::inttostring(int(zeta(r,c)))+ " ";
					  else
						  s = s + ST::inttostring(int(zeta(r,c)));
				  }
				  
			  if(freq.size()==0)
			  {
				  freq.push_back( modfreq(s,zeta.get_nedge(), 1));
			  }
			  else
			  {
				  for(unsigned int i=0; i<freq.size(); i++)
				  {
					  if(s==freq[i].model)
					  {
						  freq[i].freq++;
						  found=true;
						  i=freq.size();
					  }
				  }
				  if(found==false)
					  freq.push_back(modfreq(s,zeta.get_nedge(),1));
			  }
		  }  */

		  	 //	  FULLCOND::update();    // this command should be included (at the end of the
                           // function) to update automatically the curent
                           // mean, variance,minimum,maximum etc
/*	}










 // FUNCTION: store_model
 // TASK: stores the models in a vector
 
 void FULLCOND_rj::store_model(void)
{
	 unsigned int r,c, i, count;

	 int * zetap;
	 zetap = zeta.getV();

	  ST::string s, s_actual;
	  bool found = false;

	  for(r=0; r<nvar; r++)
	  {
		  for(c=0, count=0; c<nvar; c++, count++, zetap++)
		  {
			  if(((count+1)%nvar==0) && (count!=nvar*nvar))
			  {
				  // s = s + ST::inttostring(int(zeta(r,c)))+ " ";
				  s = s + ST::inttostring(*zetap)+ " ";

			  }
			  else
			  {
				  // s = s + ST::inttostring(int(zeta(r,c)));
				  s = s + ST::inttostring(*zetap);
			  }
		  }
	  }
		  
	  if(freq.size()==0)
	  {
		  freq.push_back( modfreq(s,zeta.get_nedge(), 1));
	  }
	  else
	  {
		 unsigned end=freq.size() ; 
		  
		 for(i=0; i<end; i++)
		  {
			  if(s==freq[i].model)
			  {
				  freq[i].freq++;
				  found=true;
				  i=end;
			  }
		  }
		  if(found==false)
			  freq.push_back(modfreq(s,zeta.get_nedge(),1)); 
		
/*
		  i=0;
		  std::vector < modfreq>::iterator pos_i;
		  pos_i = freq.begin();

		  i=0;
		  unsigned size=freq.size(); 

	//	  cout<<s<<endl; 

		  while(found==false)
		  {
			  s_actual = (*pos_i).model;
			  if(s > s_actual)
			  {
				  if(i+1<size)
				  {
					i++;
					++pos_i;
				  }
				  else
				  {
					  freq.push_back(modfreq(s,zeta.get_nedge(),1));
					  found=true;
				  }
			  }
			  else if(s==s_actual)
			  {
				  (*pos_i).freq++;
				  found=true;
			  }
			  else
			  {
				  //freq.push_back(modfreq(s,zeta.get_nedge(),1));
				  modfreq new_elem = modfreq(s,zeta.get_nedge(),1); 
				  freq.insert(pos_i,new_elem);
				  found=true;
			  }
		  } */


	//	  for(i=0; i<freq.size(); i++)
	//		  cout<<freq[i].model<<"		"<<freq[i].freq<<endl; 


	//	  cout<<endl<<endl<<endl; 

  
/*	  }
  } 




 // FUNCTION: update_zeta
 // TASK: updates zetamean and the auxiliary variables like zeta_max etc.
 
 void FULLCOND_rj::update_zeta(void)
 {

	  if((optionsp->get_nriter() > optionsp->get_burnin()) &&
		 (optionsp->get_nriter() % (optionsp->get_step()) == 0))
	  {
		  register unsigned i;

		  int * workzeta = zeta.getV();
		  double * workzetamean = zetamean.getV();

		  unsigned samplesize = optionsp->get_samplesize();

		  for(i=0;i<nrpar;i++,workzeta++,workzetamean++)
		  {
			  // updating betamean
              if (samplesize==1)
			    *workzetamean= double(*workzeta);
			  else
			    *workzetamean = (1.0/(samplesize))*
					 ((samplesize-1)*(*workzetamean) + double(*workzeta));
		  }  // end: for i=0; ...
	  } // end: if ( (nriter > optionsp->burnin) && (nriter % optionsp->step == 0) )
  }




  // FUNCTION: setzeta
  // TASK: initializes zeatmean etc.

  void FULLCOND_rj::setzeta(const adja & zetanew)
  {
	  zetamean = zeta;
	  zetameanold = zeta;

   /* zetas2 = beta;
	  betavar = beta;
	  betamin = beta;
	  betamax = beta;
	  betaqu10 = beta;
	  betaqu50 = beta;
	  betaqu90 = beta; 
	  betavarold = beta;
	  betaminold = beta;
	  betamaxold = beta;  */
//  }








} //namespace




