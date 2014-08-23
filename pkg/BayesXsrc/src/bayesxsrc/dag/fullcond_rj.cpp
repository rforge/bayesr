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





#include "fullcond_rj.h"
#include <set>
#include <algorithm>
#include <math.h>

using std::ios;

namespace MCMC
{



  FULLCOND_rj::FULLCOND_rj (MCMCoptions * o, const datamatrix & d,
						  const ST::string & t,const unsigned & r,
						  const unsigned & c, const ST::string & fp)
				: FULLCOND(o,d,t,r,c,fp)
{
	nvar=c;
	nobs=d.rows();
	assert(c==r);

//	std::vector < MCMC::FULLCOND_dag *> test;
//	preg_mods = test;

	//ini_structure();
	ini_ratio();
	ini_hyperpar();

	zetamean = datamatrix(nvar,nvar,0);
    nrpar = nvar*nvar;
//	setzeta(zeta);

	set_options();

	file_of_results = true;
	path_res = "c:\\results.res";


}  //constructor_1






FULLCOND_rj::FULLCOND_rj (vector < FULLCOND_dag * > dagp, MCMCoptions * o, const datamatrix & d,
						  const ST::string & t,const unsigned & r,
						  const unsigned & c, const ST::string & fp)
				: FULLCOND(o,d,t,r,c,fp)
{
	nvar=c;
	nobs=d.rows();
	assert(c==r);

	preg_mods= dagp ;
	ini_structure();

	ini_ratio();
	ini_hyperpar();

	zetamean = datamatrix(nvar,nvar,0);
    nrpar = nvar*nvar;
//	setzeta(zeta);

	set_options();

	mixed_case=false;

	file_of_results = false;
	path_res ="c:\\results.res";


}  //constructor_1







// CONSTRUCTOR_2

FULLCOND_rj::FULLCOND_rj (ST::string fix, const ST::string & rp, unsigned int lim,
						  double alph, ST::string switch_t, ST::string print_mod,
						  unsigned & type, vector < FULLCOND_dag * > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp)
				: FULLCOND( o,d,t,r,c,fp)
{
	setbeta(1,1,0);
	nvar=c;
	nobs=d.rows();

	preg_mods = dagp;

	set_options( lim, alph, switch_t,print_mod, fix);

    ini_structure(type);

	family = preg_mods[0]->get_family();

	ini_ratio();
	ini_hyperpar();

	zetamean = datamatrix(nvar,nvar,0);
    nrpar = nvar*nvar;

//	setzeta(zeta);

//	set_options( lim, alph, switch_t,print_mod, fix);

	mixed_case=false;

	file_of_results = true;
	path_res =rp;

}  //constructor_2






 // COPY CONSTRUCTOR

  FULLCOND_rj::FULLCOND_rj(const FULLCOND_rj & fc) : FULLCOND(FULLCOND(fc))
  {
	  nvar = fc.nvar;
	  nobs = fc.nobs;
	  type = fc.type;

	  zeta=fc.zeta;
	  zeta_fix=fc.zeta_fix;

	  preg_mods = fc.preg_mods;

	  freq = fc.freq;

	  sigma_prop = fc.sigma_prop;

	  acceptance_b = fc.acceptance_b;
      acceptance_d = fc.acceptance_d;
      acceptance_s = fc.acceptance_s;

	  nrtrials_b = fc.nrtrials_b;
	  nrtrials_d = fc.nrtrials_d;
	  nrtrials_s = fc.nrtrials_s;

	  alpha_sig_i = fc.alpha_sig_i;
	  beta_sig_i  = fc.beta_sig_i;
	  alpha_sig_j = fc.alpha_sig_j;
	  beta_sig_j  = fc.beta_sig_j;

	  gamma_nobs1= fc.gamma_nobs1;
	  gamma_nobs2= fc.gamma_nobs2;

	  step_aborted = fc.step_aborted;

	  zetamean = fc.zetamean;
	  zetameanold = fc.zetameanold;

	  limit_number = fc.limit_number;
      alpha = fc.alpha;

	  switch_type = fc.switch_type;

	  print_models = fc.print_models;
	  mixed_case = fc.mixed_case;
	  file_of_results = fc.file_of_results;
	  conditions = fc.conditions;

	  path_res = fc.path_res;
      family = fc.family;
  }




  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_rj & FULLCOND_rj::operator=(const FULLCOND_rj & fc)
  {
	  if (this==&fc)
		  return *this;

     FULLCOND::operator=(FULLCOND(fc));

	  nvar = fc.nvar;
	  nobs = fc.nobs;
	  type = fc.type;

	  preg_mods = fc.preg_mods;

	  zeta=fc.zeta;
	  zeta_fix=fc.zeta_fix;

	  freq = fc.freq;

	  sigma_prop = fc.sigma_prop;

	  acceptance_b = fc.acceptance_b;
      acceptance_d = fc.acceptance_d;
      acceptance_s = fc.acceptance_s;

	  nrtrials_b = fc.nrtrials_b;
	  nrtrials_d = fc.nrtrials_d;
	  nrtrials_s = fc.nrtrials_s;

	  alpha_sig_i = fc.alpha_sig_i;
	  beta_sig_i  = fc.beta_sig_i;
	  alpha_sig_j = fc.alpha_sig_j;
	  beta_sig_j  = fc.beta_sig_j;

	  gamma_nobs1= fc.gamma_nobs1;
	  gamma_nobs2= fc.gamma_nobs2;

	  step_aborted = fc.step_aborted;

	  zetamean = fc.zetamean;
	  zetameanold = fc.zetameanold;

	  limit_number = fc.limit_number;
      alpha = fc.alpha;

	  switch_type = fc.switch_type;

	  print_models = fc.print_models;
	  mixed_case = fc.mixed_case;
	  file_of_results = fc.file_of_results;
	  conditions = fc.conditions;

	  path_res = fc.path_res;
      family = fc.family;

	  return *this;
  }




// FUNCTION: rj_step
// chooses randomly new edge, decides which kind of step and makes it

void FULLCOND_rj::rj_step(void)
{
	step_aborted = true;

	unsigned int vertex_i;
	unsigned int vertex_j;

	while (step_aborted==true)
	{
		// two variables are randomly chosen
		vertex_i = rand() % nvar;
		vertex_j = vertex_i;

		while(vertex_i==vertex_j)
			vertex_j= rand() % nvar;





		// ************** decide which kind of step (birth, death, switch) *************************
		if(conditions == false)
		{
			if (zeta(vertex_i,vertex_j) == 1)
			{
				death_step(vertex_i,vertex_j);
			}
			else if (zeta(vertex_j,vertex_i)==1)
			{
				switch_step(vertex_i,vertex_j);
			}
			else
			{
				birth_step(vertex_i,vertex_j);
			}
		}
		else	// i.e. conditions == true
		{
			if (zeta(vertex_i,vertex_j) == 1 )
			{
				if(conditions_okay_d(vertex_i,vertex_j) == true)
				{
					//cout<<"try_d: "<<vertex_i<<" -> "<<vertex_j<<endl;
					death_step(vertex_i,vertex_j);
				}
			}
			else if (zeta(vertex_j,vertex_i) == 1)
			{
				if(conditions_okay_s(vertex_i,vertex_j) == true)
					switch_step(vertex_i,vertex_j);
			}
			else
			{
				if(conditions_okay_b(vertex_i,vertex_j) == true)
				{
					//cout<<"try_b: "<<vertex_i<<" -> "<<vertex_j<<endl;
					birth_step(vertex_i,vertex_j);
				}
			}
		}
	}

	step_aborted=false;
}




// FUNCTION: birth_step
// makes Birth step

void FULLCOND_rj::birth_step(unsigned int v_i, unsigned int v_j)
{
	if(zeta.azy_test(v_i,v_j)== true)
	{

		unsigned int ncoef_old = preg_mods[v_j]->get_ncoef();
		unsigned int ncoef_new = ncoef_old + 1;

		if(mixed_case==true)
				preg_mods[v_j]->create_matrices("b", ncoef_new);

		// instead of: datamatrix b_new (ncoef_new,1);
// Vorschlag:
//		datamatrix & b_new = preg_mods[v_j]->get_b_new_b();
		datamatrix b_new = preg_mods[v_j]->get_b_new_b();
		// instead of: datamatrix x_new (nobs,ncoef_new);
// Vorschlag:
//		datamatrix & x_new = preg_mods[v_j]->get_x_new_b() ;
		datamatrix x_new = preg_mods[v_j]->get_x_new_b() ;
		// instead of: datamatrix xx_new (ncoef_new,ncoef_new);
// Vorschlag:
//		datamatrix & xx_new = preg_mods[v_j]->get_xx_new_b() ;
		datamatrix xx_new = preg_mods[v_j]->get_xx_new_b() ;

		assert(ncoef_new==b_new.rows());

		double beta_new = rand_normal(); //coefficient which will be added

		// computing of the new values
		make_new_b("b", v_i,v_j,beta_new, xx_new,b_new, x_new);


		// calculate ratio

		double log_denom;
		double log_num;
		double ratio;

		log_num = preg_mods[v_j]->calc_SQT_x(x_new, b_new) + preg_mods[v_j]->calc_SQT_b(b_new);
		log_denom = preg_mods[v_j]->calc_SQT_x() +  preg_mods[v_j]->calc_SQT_b();

		ratio = -1/(2*preg_mods[v_j]->get_sigma_i())
				*(log_num - log_denom) - p_prop(beta_new) ;


		//accept proposal and change the corresponding values

		if(func::accept(ratio) == true)
		{
			preg_mods[v_j]->change_adcol(v_i,1);
			preg_mods[v_j]->change(v_i, b_new, x_new, xx_new, ncoef_new);

			acceptance_b ++;
			zeta.edge_plus();

			zeta(v_i,v_j)=1;
			zeta.change_list(v_i,v_j,0);

			/*********************
				unsigned k,l;
				for(k=0; k<nvar; k++)
				{
					for(l=0; l<nvar; l++)
						cout<<zeta(k,l)<<" ";
					cout<<endl;
				}
				cout<<endl
					<<endl;;
			*********************/
		}

		nrtrials_b ++;
		step_aborted = false;

	} // end of "if(azy_test(i,j)== true) ...."
}




void FULLCOND_rj::death_step(unsigned int v_i, unsigned int v_j)
{
	unsigned int ncoef_old = preg_mods[v_j]->get_ncoef();
	unsigned int ncoef_new = ncoef_old - 1;

	if(mixed_case==true)
				preg_mods[v_j]->create_matrices("d", ncoef_new);

	// instead of: datamatrix b_new (ncoef_new,1);
	datamatrix & b_new = preg_mods[v_j]->get_b_new_d();
	// instead of: datamatrix x_new (nobs,ncoef_new);
	datamatrix & x_new = preg_mods[v_j]->get_x_new_d();
	// instead of: datamatrix xx_new (ncoef_new,ncoef_new);
	datamatrix & xx_new = preg_mods[v_j]->get_xx_new_d();


	if(b_new.rows() != ncoef_new)
	{
		cout<<"main_effects: "<<preg_mods[v_j]->get_ncoef_m()<<endl;
		cout<<"interactions: "<<preg_mods[v_j]->get_ncoef_ia()<<endl;
		cout<<"all: "<<preg_mods[v_j]->get_ncoef()<<endl;
	}


	double beta_old; //coefficient which will vanish

	// computing of the new values
	make_new_d("d", v_i,v_j,xx_new, beta_old, b_new, x_new);



	// calculate ratio

    double ratio;
	double log_num;
	double log_denom;

	log_num = preg_mods[v_j]->calc_SQT_x(x_new, b_new) + preg_mods[v_j]->calc_SQT_b(b_new);
	log_denom = preg_mods[v_j]->calc_SQT_x() + preg_mods[v_j]->calc_SQT_b();

	ratio = -1/(2*preg_mods[v_j]->get_sigma_i())
			*(log_num - log_denom) + p_prop(beta_old) ;




	// accept proposal and change the corresponding values

	if(func::accept(ratio) == true)
	{
        zeta(v_i,v_j)=0;
		zeta.change_list(v_i,v_j,1);
		preg_mods[v_j]->change_adcol(v_i,0);
		preg_mods[v_j]->change(v_i, b_new, x_new, xx_new, ncoef_new);

		acceptance_d ++;
		nrtrials_d ++;
		zeta.edge_minus();

		/*********************
			unsigned k,l;
			for(k=0; k<nvar; k++)
			{
				for(l=0; l<nvar; l++)
					cout<<zeta(k,l)<<" ";
				cout<<endl;
			}
			cout<<endl
				<<endl;;
		*********************/
	}

	nrtrials_d ++;
	step_aborted = false;
}





// FUNCTION: switch_step
// makes switch step from j->i to i->j

void FULLCOND_rj::switch_step(unsigned int v_i, unsigned int v_j)
{

	zeta(v_j,v_i)=0;
	zeta.change_list(v_j,v_i,1);

	if(zeta.azy_test(v_i,v_j)==true)
	{

		zeta(v_j,v_i)=1;
		zeta.change_list(v_j,v_i,0);

		if(switch_type=="equi"
			&& zeta.equi_test(v_j,v_i)==true)
		{
			FULLCOND_rj::switch_version_2(v_i,v_j);
		}
		else if (switch_type=="normal")
		{
			FULLCOND_rj::switch_version_1(v_i,v_j);
		}
		else if (switch_type=="mix")
		{
			if(zeta.equi_test(v_j,v_i)==true)
			{
				FULLCOND_rj::switch_version_2(v_i,v_j);
			}
			else
			{
				FULLCOND_rj::switch_version_1(v_i,v_j);
			}
		}


	} // end of "if(azy_test(i,j)==true)"

	else // azy_test(i,j)==false"
	{
		zeta(v_j,v_i)=1;
		zeta.change_list(v_j,v_i,0);
	}

	nrtrials_s ++;
}












void FULLCOND_rj::switch_version_1(unsigned i, unsigned j)
{


	//************ DEATH-step for edge j->i  ******************

	unsigned int ncoef_old_i = preg_mods[i]->get_ncoef();
	unsigned int ncoef_new_i = ncoef_old_i - 1;

	if(mixed_case==true)
			preg_mods[i]->create_matrices("d", ncoef_new_i);

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

	datamatrix sig_mean_i (ncoef_new_i,ncoef_new_i);
	sig_mean_i.assign(xx_new_i.inverse());

	datamatrix mean_i (ncoef_new_i,1);

	mean_i.mult(sig_mean_i, x_new_i.transposed()*preg_mods[i]->get_y());

	sigma_new_i = sample_sigma('i', i, ncoef_new_i, mean_i, x_new_i);

	b_new_i.mult(sig_mean_i.root(),rand_normvek(ncoef_new_i));
	b_new_i.plus(b_new_i,mean_i);


	preg_mods[i]->calc_lin_prop(x_new_i, b_new_i);








	//************ BIRTH-step for edge i->j ******************

	unsigned int ncoef_old_j = preg_mods[j]->get_ncoef();
	unsigned int ncoef_new_j = ncoef_old_j + 1;

	if(mixed_case==true)
				preg_mods[j]->create_matrices("b", ncoef_new_j);

	// instead of: datamatrix b_new_j (ncoef_new_j,1);
// Vorschlag:
//  datamatrix & b_new_j = preg_mods[j]->get_b_new_b();
	datamatrix b_new_j = preg_mods[j]->get_b_new_b();
	// instead of: datamatrix x_new_j (nobs,ncoef_new_j);
// Vorschlag:
// 	datamatrix & x_new_j = preg_mods[j]->get_x_new_b() ;
	datamatrix x_new_j = preg_mods[j]->get_x_new_b() ;
	// instead of: datamatrix xx_new_j (ncoef_new_j,ncoef_new_j);
// Vorschlag:
//	datamatrix & xx_new_j = preg_mods[j]->get_xx_new_b() ;
	datamatrix xx_new_j = preg_mods[j]->get_xx_new_b() ;

	double beta_new_ji=0;

	// computing of the new values
	make_new_b("s", i,j, beta_new_ji, xx_new_j,b_new_j, x_new_j);


	datamatrix sig_mean_j (ncoef_new_j,ncoef_new_j);
	sig_mean_j.assign(xx_new_j.inverse());

	datamatrix mean_j (ncoef_new_j,1);
	mean_j.mult(sig_mean_j, x_new_j.transposed()*preg_mods[j]->get_y());

	sigma_new_j = sample_sigma('j', j, ncoef_new_j, mean_j, x_new_j);

	b_new_j.mult(sig_mean_j.root(),rand_normvek(ncoef_new_j));
	b_new_j.plus(b_new_j,mean_j);

	preg_mods[j]->calc_lin_prop(x_new_j, b_new_j);









	// **************  calculate ratio and accept or reject *************************************

	double ratio = 	ratio_s( i, j, b_new_i, b_new_j, x_new_i, x_new_j, mean_i, mean_j,
								     xx_new_i, xx_new_j, sig_mean_i, sig_mean_j,
									 sigma_new_i, sigma_new_j);

	//ratio_s(i, j, ncoef_new_i, ncoef_new_j, b_new_i, b_new_j, x_new_i, x_new_j, beta_old_i, beta_new_j, mean_prop_i, mean_prop_j);

	if(func::accept(ratio) == true)
	{
        zeta(j,i)=0;
		zeta(i,j)=1;
		zeta.change_list(i,j,2);  //last argument=2, because switch-step

		preg_mods[i]->change_adcol(j,0);
		preg_mods[j]->change_adcol(i,1);

		preg_mods[i]->change(j, b_new_i, x_new_i, xx_new_i, ncoef_new_i);
		preg_mods[j]->change(i, b_new_j, x_new_j, xx_new_j, ncoef_new_j);

		acceptance_s ++;
	}

	step_aborted = false;

	}












void FULLCOND_rj::switch_version_2(unsigned v_i, unsigned v_j)
{
	if(randnumbers::uniform()<0.5)
	{


		//************ DEATH-step for edge j->i  ******************

		unsigned int ncoef_old_i = preg_mods[v_i]->get_ncoef();
		unsigned int ncoef_new_i = ncoef_old_i - 1;

		if(mixed_case==true)
				preg_mods[v_i]->create_matrices("d", ncoef_new_i);

		// instead of: datamatrix b_new_i (ncoef_new_i,1);
		datamatrix & b_new_i = preg_mods[v_i]->get_b_new_d();
		// instead of: datamatrix x_new_i (nobs,ncoef_new_i);
		datamatrix & x_new_i = preg_mods[v_i]->get_x_new_d();
		// instead of: datamatrix xx_new_i (ncoef_new_i,ncoef_new_i);
		datamatrix & xx_new_i = preg_mods[v_i]->get_xx_new_d();

		double beta_old_ij; //coefficient which will vanish
		//double sigma_new_i; //coefficient which will vanish
		//double sigma_new_j; //coefficient which will vanish

		// computing of the new values
		make_new_d("d", v_j,v_i,xx_new_i, beta_old_ij, b_new_i, x_new_i);





		//************ BIRTH-step for edge i->j ******************

		unsigned int ncoef_old_j = preg_mods[v_j]->get_ncoef();
		unsigned int ncoef_new_j = ncoef_old_j + 1;

		if(mixed_case==true)
				preg_mods[v_j]->create_matrices("b", ncoef_new_j);

		// instead of: datamatrix b_new_j (ncoef_new_j,1);
// Vorschlag:
//		datamatrix & b_new_j = preg_mods[v_j]->get_b_new_b();
		datamatrix b_new_j = preg_mods[v_j]->get_b_new_b();
		// instead of: datamatrix x_new_j (nobs,ncoef_new_j);
// Vorschlag:
//		datamatrix & x_new_j = preg_mods[v_j]->get_x_new_b() ;
		datamatrix x_new_j = preg_mods[v_j]->get_x_new_b() ;
		// instead of: datamatrix xx_new_j (ncoef_new_j,ncoef_new_j);
// Vorschlag:
//		datamatrix & xx_new_j = preg_mods[v_j]->get_xx_new_b() ;
		datamatrix xx_new_j = preg_mods[v_j]->get_xx_new_b() ;

		double beta_new_ji=0;


		// computing of the new values
		make_new_b("b", v_i,v_j, beta_new_ji, xx_new_j,b_new_j, x_new_j);

		zeta(v_j,v_i)=0;
		zeta(v_i,v_j)=1;
		zeta.change_list(v_i,v_j,2);  //last argument=2, because switch-step

		preg_mods[v_i]->change_adcol(v_j,0);
		preg_mods[v_j]->change_adcol(v_i,1);

		preg_mods[v_i]->change(v_j, b_new_i, x_new_i, xx_new_i, ncoef_new_i);
		preg_mods[v_j]->change(v_i, b_new_j, x_new_j, xx_new_j, ncoef_new_j);

		acceptance_s ++;

		/*********************
		unsigned k,l;
		for(k=0; k<nvar; k++)
		{
			for(l=0; l<nvar; l++)
				cout<<zeta(k,l)<<" ";
			cout<<endl;
		}
		cout<<endl
			<<endl;;
		*********************/
	}

	step_aborted = false;
}






  // FUNCTION: ratio_s
 // TARGET: computes acceptance ratio in the switch step
 double FULLCOND_rj::ratio_s(unsigned int i,unsigned int j,
							const datamatrix & b_new_i, const datamatrix & b_new_j,
							const datamatrix & x_new_i, const datamatrix & x_new_j,
							const datamatrix & mean_i, const datamatrix & mean_j,
							const datamatrix & xx_new_i, const datamatrix & xx_new_j,
							const datamatrix & sig_mean_i, const datamatrix & sig_mean_j,
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

	help_i_old = preg_mods[i]->p_prop(preg_mods[i]->get_beta_mean(), preg_mods[i]->get_beta_help(), preg_mods[i]->get_sigma_i());



	double help1 = preg_mods[i]->b_distr();
    double help2 = preg_mods[j]->b_distr();
    double help3 = -0.5* ( sig_mean_i.det() + preg_mods[i]->p_prop(b_new_i, mean_i, xx_new_i));
    double help4 = -0.5* ( sig_mean_j.det() + preg_mods[i]->p_prop(b_new_j, mean_j, xx_new_j));


     log_prop_b = help1+help2-help3-help4;


	double log_ratio = log_num1 - log_den1
						+ log_num2 - log_den2
						+ log_num3 - log_den3
						+ log_num4 - log_den4
						+ log_prop_b
						+ log_prop_s;

	return log_ratio;
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



	   //*********************** compute x_new ************************************
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






		//*********************** compute xx_new ************************************

		double * workxx_new;
		double * workxx;
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







   // FUNCTION: make_new_d
   // TASK: computes the new values for a death-step

   void FULLCOND_rj::make_new_d (ST::string step, unsigned int i, unsigned int j, datamatrix & xx_new,
							double & beta_old, datamatrix & b_new, datamatrix & x_new)
	{
		unsigned int k,l;
		unsigned int ncoef=preg_mods[j]->get_ncoef();

		unsigned * workzeta = zeta.getV() + j;





		// Definition of t:
		// the variable i is the t-th regression coefficient of j (without intercept)
		unsigned t=0;
		for(k=0; k<nvar; k++)
		{
			if( (k<=i) && (*workzeta==1))
				t++;
			workzeta = workzeta+nvar;
		}





		//************************** compute x_new ************************************
		double * workx_new;
		double * workx;

		workx_new = x_new.getV();
		workx = preg_mods[j]->getV_x();

		for(k=0; k<nobs; k++)
		{
			for(l=0; l<ncoef; l++, workx_new++, workx++)
			{
				if(l != t)
					*workx_new = *workx;
				else
					workx_new--;
			}
		}




		//************************* compute xx_new ************************************
		double * workxx_new;
		double * workxx;

		workxx_new = xx_new.getV();
		workxx = preg_mods[j]->getV_xx();

		for(k=0; k<ncoef; k++)
		{
			if(k != t)
			{
				for(l=0; l<ncoef; l++, workxx++, workxx_new++)
				{
					if(l != t)
						*workxx_new = *workxx;
					else
						workxx_new--;
				}
			}
			else
				workxx = workxx + ncoef;
		}









		//save the t-th coefficient (which is going to be deleted)
		beta_old = preg_mods[j]->get_beta_help(t,0);



		if(step != "s") //if the death-step is NOT part of a switch step
		{
			//compute proposed beta
			double * workb_new;
			double * workbeta;

			workb_new = b_new.getV();
			workbeta = preg_mods[j]->getV_beta_help();

			for(k=0; k<ncoef; workb_new++, workbeta++, k++)
			{
				if(k!= t)
					*workb_new = *workbeta;
				else
					workb_new--;

			}

			//compute proposed linear predictor
			preg_mods[j]->calc_lin_prop(x_new, b_new);
		}
	}





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
	// TARGET: returns the log-density of the proposal u, which is normaldistributed

	double FULLCOND_rj::p_prop(double prop)
	{
		double pi = 3.14159;
		double help =  2*pi*sigma_prop ;
		return -1/2 * log(help) - prop*prop / (2*sigma_prop);
	}








	// FUNCTION: p_prop()
	// TARGET: returns the log-density of the proposal u, which is normaldistributed

	double FULLCOND_rj::p_prop(const datamatrix & prop)
	{
		unsigned t;
		double sum;
		sum=0;

		for(t=0; t<prop.rows(); t++)
			sum = sum + p_prop(prop(t,0));

		return sum;
	}









	// FUNCTION: outresults
	// TASK: writes estimation results to logout or into a file (after estimation)

	void FULLCOND_rj::outresults(void)
    {

		unsigned k, l;


//		FULLCOND::outresults();

		std::sort(freq.begin(), freq.end());

		make_list_essential();
		std::sort(list_ess.begin(), list_ess.end());

		outres_dags();
		outres_essentials();




        optionsp->out("\n");
        optionsp->out("\n");
        optionsp->out("\n");
		optionsp->out("******** averaged adjacency matrix: ********\n");
		optionsp->out("\n");

        ST::string help;

		for(k=0; k<nvar; k++)
		{
        help="";
			for(l=0; l<nvar; l++)
				help = help + ST::doubletostring(zetamean(k,l),2) + "	";

            optionsp->out(help + "\n");
			optionsp->out("\n");
		}


		optionsp->out("\n");
        optionsp->out("\n");
        optionsp->out("\n");
		optionsp->out("********  mean of sceletons:  ********\n");
		optionsp->out("\n");

		for(k=0; k<nvar; k++)
		{
        help="";
			for(l=0; l<nvar; l++)
			{
				if(k<l)
					help = help + ST::doubletostring(zetamean(k,l)+zetamean(l,k),2) + "	";
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


		for(k=0; k<nvar; k++)
		{
        help="";
			for(l=0; l<nvar; l++)
				help = help + ST::doubletostring(correlation(k,l),2) + "	";

            optionsp->out(help + "\n");
			optionsp->out("\n");
		}


		optionsp->out("\n");
		optionsp->out("\n");
		optionsp->out("******** partial correlation ********\n");
		optionsp->out("\n");



		for(k=0; k<nvar; k++)
		{
        help="";
			for(l=0; l<nvar; l++)
				help = help + ST::doubletostring(partial_corr(k,l),2) + "	";

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

		if( nrtrials_b !=0)
            optionsp->out("acceptance ratio if a birth-step has already been proposed: "
			                          + ST::doubletostring(double(acceptance_b)/nrtrials_b,3) + "\n");
        if( nrtrials_d !=0)
            optionsp->out("acceptance ratio if a death-step has already been proposed: "
			                          + ST::doubletostring(double(acceptance_d)/nrtrials_d,3) + "\n");
        if( nrtrials_s !=0)
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













	// FUNCTION: outres_dags
	// TASK: writes out the different dags
	void FULLCOND_rj::outres_dags(void)
	{
		int k;
        unsigned freq_size = freq.size();	// size of vector freq

		unsigned total_number =0;		// sample size
		double prob_of_limit = 0;		// probability of first n dags, where n=limit_number

		ST::string print_mod;			// local variable for print_models






		//******************** calculate necessary numbers ********************************************

		for(k=freq_size-1; k>=0; k--)
		{
			total_number = total_number + freq[k].freq;
		}

		if(limit_number<freq_size)
		{
			for(k=freq_size-1; k>=freq_size-limit_number; k--)
			{
				prob_of_limit = prob_of_limit + double(freq[k].freq)/total_number;
			}
		}
		else
		{
			for(k=freq_size-1; k>=0; k--)
			{
				prob_of_limit = prob_of_limit + double(freq[k].freq)/total_number;
			}
		}






		//******************** check which "print_model"-option can be realized ***********************************

		if(		(print_models=="all")
			||	((print_models=="limit") && (limit_number >= freq_size))
			||  ((print_models=="normal") && (1-alpha >= prob_of_limit) && (limit_number >= freq_size)))
		{
			print_mod="all";
		}
		else if ((print_models=="limit")
			||   ((print_models=="normal") && (1-alpha >= prob_of_limit)))
		{
			print_mod="limit";
		}
		else if ((print_models=="prob")
			||   (print_models=="normal") && (1-alpha < prob_of_limit))
		{
			print_mod="prob";
		}
		else
			cout<<"strange...."<<endl;





		//********************* start to write output ************************************

		optionsp->out("Number of different dags visited by the algorithm: "
						+ ST::inttostring(freq.size())
						+ "\n"
						+ "\n");

		optionsp->out("******** DIFFERENT MODELS sorted by frequencies  ********\n");








		//************** different outputs with respect to the option chosen above ************

		if(print_mod=="all")
		{
			optionsp->out("********************** all models **********************\n");
			optionsp->out("\n");

			for(k=freq_size-1; k>=0; k--)
			{
				optionsp->out(freq[k].model
					+ "	"
					+ ST::inttostring(freq[k].nedges)
					+ "	"
					+ ST::inttostring(freq[k].freq)
					+ "	"
					+ ST::doubletostring(double(freq[k].freq)/total_number,3)
					+ "\n");
			}
		}
		else if (print_mod=="limit")
		{
			optionsp->out("******** first "	+ ST::inttostring(limit_number)
											+ " of the most important models ********\n");
			optionsp->out("\n");

			for(k=freq_size-1; k>=freq_size-limit_number; k--)
			{
				optionsp->out(freq[k].model
					+ "	"
					+ ST::inttostring(freq[k].nedges)
					+ "	"
					+ ST::inttostring(freq[k].freq)
					+ "	"
					+ ST::doubletostring(double(freq[k].freq)/total_number,3)
					+ "\n");
			}

		}
		else if (print_mod=="prob")
		{
			double rel_freq;
			double prob_cum = 0;

			optionsp->out("******* at least " + ST::doubletostring(100-alpha*100)
											  + " % of posterior probability ******** \n");
			optionsp->out("\n");

			k=freq_size-1;

			while(prob_cum<1-alpha && k>=0)
			{
				rel_freq = double(freq[k].freq)/double(total_number);
					optionsp->out(freq[k].model
						+ "	"
						+ ST::inttostring(freq[k].nedges)
						+ "	"
						+ ST::inttostring(freq[k].freq)
						+ "	"
						+ ST::doubletostring(rel_freq,3)
						+ "\n");

				prob_cum = prob_cum + rel_freq;
				k--;
			}
		}
	}














  // FUNCTION: outres_essentials
  // TASK: writes out the different essential graphs
  void FULLCOND_rj::outres_essentials(void)
  {


	  unsigned total_number_ess;	// sample size
	  unsigned list_size;
      int k;


	  double prob_of_limit;		// probability of first n dags, where n=limit_number
	  double rel_freq;
	  double prob_cum = 0;

	  ST::string print_mod;		// local variable for print_models







	  //************** calculate necessary numbers *************************************

	  total_number_ess = 0;
      list_size = list_ess.size();

	  for(k=list_size-1; k>=0; k--)
	  {
			total_number_ess = total_number_ess +list_ess[k].freq;
	  }

	  prob_of_limit = 0;
	  if(limit_number<list_size)
		{
			for(k=list_size-1; k>=list_size-limit_number; k--)
			{
				prob_of_limit = prob_of_limit + double(list_ess[k].freq)/total_number_ess;
			}
		}
		else
		{
			for(k=list_size-1; k>=0; k--)
			{
				prob_of_limit = prob_of_limit + double(list_ess[k].freq)/total_number_ess;
			}
		}






		//******************** check which "print_model"-option can be realized ***********************************

		if(		(print_models=="all")
			||	((print_models=="limit") && (limit_number >= list_size))
			||  ((print_models=="normal") && (1-alpha >= prob_of_limit) && (limit_number >= list_size)))
		{
			print_mod="all";
		}
		else if ((print_models=="limit")
			||   ((print_models=="normal") && (1-alpha >= prob_of_limit)))
		{
			print_mod="limit";
		}
		else if ((print_models=="prob")
			||   (print_models=="normal") && (1-alpha < prob_of_limit))
		{
			print_mod="prob";
		}
		else
			cout<<"strange...."<<endl;







		//********************* start to write output ************************************

		optionsp->out("\n");
		optionsp->out("\n");
		optionsp->out("\n");
		optionsp->out("\n");
		optionsp->out("\n");
		optionsp->out("\n");
		optionsp->out("Number of different equivalent classes visited by the algorithm: "
						+ ST::inttostring(list_size)
						+ "\n"
						+ "\n"
						+ "\n");
		optionsp->out("************* DIFFERENT EQUIVALENCE CLASSES sorted by frequencies  *************\n");









		//************** different outputs with respect to the option chosen above ************

		if(print_mod=="all")
		{
			optionsp->out("********************************** all models **********************************\n");
			optionsp->out("\n");

			for(k=list_size-1; k>=0; k--)
				write_out_essential(list_ess[k], total_number_ess);
		}
		else if (print_mod=="limit")
		{
			optionsp->out("************** first " + ST::inttostring(limit_number)
				+ " of the most important equivalent classes **************\n");
			optionsp->out("\n");

            for(k=list_size-1; k>=list_size-limit_number; k--)
				write_out_essential(list_ess[k], total_number_ess);
		}
		else if (print_mod=="prob")
		{
			optionsp->out("******************* at least "
							+ ST::doubletostring(100-alpha*100)
							+ " % of posterior probability ********************* \n");
			optionsp->out("\n");
			optionsp->out("\n");

			k=list_size-1;
			while(prob_cum<1-alpha && k>=0)
			{
				rel_freq = double(list_ess[k].freq)/double(total_number_ess);

				write_out_essential(list_ess[k], total_number_ess);

				prob_cum = prob_cum + rel_freq;
				k--;
			}
		}
		if( file_of_results == true)
			write_out_resfile();
  }











  // FUNCTION: write_out_essential
  // TASK: writes out the essential graph m and its properties
  void FULLCOND_rj::write_out_essential(essfreq & ess, unsigned total_number_ess)
  {
	  unsigned c,r,k, size;

	  Matrix<unsigned> m;

	  m = ess.sceleton;
      ST::string s;

	  optionsp->out("Sceleton: ");

	  for(c=0; c<nvar; c++)
	  {

		  for(r=0; r<nvar; r++)
		  {
			  s = s + ST::inttostring(m(c,r));
		  }

          s = s +  " ";
	  }

		  optionsp->out(s + "\n");

		  size= (ess.immoral).size();

		  optionsp->out("\n");

		  if(size>0)
		  {
			  optionsp->out("Immoralities: ");
			  for(k=0; k<size; k++)
			  {
				  optionsp->out("("
					+ ST::inttostring(ess.immoral[k][0])
					+ ";"
					+ ST::inttostring(ess.immoral[k][1])
					+","
					+ ST::inttostring(ess.immoral[k][2])
					+")"
					+ " ");
				}
		  }
		  else
			  optionsp->out("No immoralities.");

		optionsp->out("\n");
		optionsp->out("Number of edges: "
				+ ST::inttostring(ess.nedges)
				+ "\n"
				+ "Abs.freq.: "
				+ ST::inttostring(ess.freq)
				+ "\n"
				+ "Rel.freq.: "
				+ ST::doubletostring(double(ess.freq)/total_number_ess,3)
				+ "\n");

		optionsp->out("\n");
  }









  // FUNCTION: write_out_resfile
  // TASK: writes out results (= 10 most important essential graphs and the adjacency matrix)
  // into a separate file which has to be named before; good for simulation studies

	void FULLCOND_rj::write_out_resfile(void)
	{
		int i, loop_end;
	    unsigned size_list;
		unsigned total_number_ess=0;
		adja adja_help (nvar);

        size_list = list_ess.size();

        std::sort(list_ess.begin(), list_ess.end());

		for(i=size_list-1; i>=0; i--)
		{
			total_number_ess = total_number_ess+list_ess[i].freq;
		}

	  	#if defined(MICROSOFT_VISUAL)
		{
			loop_end = __min(size_list, 10);
		}
		#else
		{
// Vorschlag:
//			loop_end = min(size_list, unsigned (10)) ;
			loop_end = std::min(size_list, unsigned (10)) ;
		}
		#endif


        loop_end = size_list - loop_end;

        for(i=size_list-1; i>=loop_end; i--)
			adja_help.write_out_ess_short(list_ess[i], path_res, total_number_ess);

		ofstream fout(path_res.strtochar(),ios::app);

		for(i=0; i<nvar; i++)
		{
			for(unsigned j=0; j<nvar; j++)
			{
				fout<<zetamean(i,j)<<" ";

			}
		}

		fout<<endl;
		fout.close();


	}





	void FULLCOND_rj::update(void)
    {

		/*********************
		unsigned k,l;
		for(k=0; k<nvar; k++)
		{
			for(l=0; l<nvar; l++)
				cout<<zeta(k,l)<<" ";
			cout<<endl;
		}
		cout<<endl;
		/*********************/

		  rj_step ();



		/*********************
		for(k=0; k<nvar; k++)
		{
			for(l=0; l<nvar; l++)
				cout<<zeta(k,l)<<" ";
			cout<<endl;
		}

		cout<<endl;
		cout<<endl;
		cout<<endl;
		*********************/

		  nrtrials++;

		  update_zeta();

		  //std::set < ST::string, unsigned int> freq;

		  if((optionsp->get_nriter() > optionsp->get_burnin()) &&
			  (optionsp->get_nriter() % (optionsp->get_step()) == 0))
		  {
			  store_model();
		  }

          /*

          // Falls nach jeder 5000. Beobachtung rausgeschrieben werden soll:

           if( (optionsp->get_nriter() - optionsp->get_burnin()) % 5000 == 0
               &&
               optionsp->get_nriter() > optionsp->get_burnin())
		  {
             outresults();
             list_ess = vector < essfreq > (0);
		  }

          */



		 // if(print_dags==true)
			  //FULLCOND::update();    // this command should be included (at the end of the
									// function) to update automatically the curent
									// mean, variance,minimum,maximum etc
	}






 // FUNCTION: make_list_essential
 // TASK: computes the list of the essential graphs from the list of all dags
 void FULLCOND_rj::make_list_essential(void)
 {
	 unsigned k;
	 unsigned num_edge, frequency;

	Matrix <unsigned> scel;

	 adja adja_help(nvar);
	 modfreq model_help;

	 vector <Matrix <unsigned> > essential(2);

	 for(k=0; k<freq.size(); k++)
	 {

		 model_help = freq[k];

		 num_edge = model_help.nedges;
		 frequency = model_help.freq;

		 scel = Matrix <unsigned> (nvar, nvar, 0);
		 vector< vector <unsigned> > imm;

		 adja_help.string_to_adja (model_help.model);

		 adja_help.adja_to_ess(scel, imm);

		 essfreq ess_new;
		 ess_new = essfreq(scel, imm, num_edge, frequency);

		 adja_help.add_ess_to_list (list_ess, ess_new);
	 }
 }




 // FUNCTION: store_model
 // TASK: stores the models in a vector

 void FULLCOND_rj::store_model(void)
{
	 unsigned int r,c, i, count;

	 unsigned * zetap;
	 zetap = zeta.getV();

	  ST::string s, s_actual;
	  bool found = false;

	  for(r=0; r<nvar; r++)
	  {
		  for(c=0, count=0; c<nvar; c++, count++, zetap++)
		  {
			  s = s + ST::inttostring(*zetap);

			  if(((count+1)%nvar==0) && (count!=nvar*nvar))
				  s = s +  " ";
		  }
	  }

	  if(freq.size()==0)
	  {
		  freq.push_back( modfreq(s,zeta.get_nedge(), 1));
	  }
	  else
	  {
  	 unsigned end=freq.size();

         std::vector < modfreq>::iterator pos_i;
         pos_i = freq.begin();

		 for(i=0; i<end; i++,++pos_i)
		  {
			  if(s==(*pos_i).model) //if(s==freq[i].model)
			  {
				  (*pos_i).freq++;
                  // freq[i].freq++;
				  found=true;
				  i=end;
			  }
		  }
		  if(found==false)
			  freq.push_back(modfreq(s,zeta.get_nedge(),1));

	  }
  }




 // FUNCTION: update_zeta
 // TASK: updates zetamean and the auxiliary variables like zeta_max etc.

 void FULLCOND_rj::update_zeta(void)
 {

	  if((optionsp->get_nriter() > optionsp->get_burnin()) &&
		 (optionsp->get_nriter() % (optionsp->get_step()) == 0))
	  {
		  register unsigned i;

		  unsigned * workzeta = zeta.getV();
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



   // FUNCTION: ini_structure
   // TASK: initializes structure
	void FULLCOND_rj::ini_structure(void)
	{

		unsigned int i;
		unsigned type = 0; 			// start with independent model

        zeta = adja(nvar,type);		// start with independent model

		for(i=0; i<nvar; i++)
		{
			preg_mods[i]->initialize(zeta,i);
		}
	}




	// FUNCTION: ini_structure
   // TASK: initializes structure
	void  FULLCOND_rj::ini_structure(unsigned t)
	{

		unsigned int k,l;
		bool okay=true;

		type = t;
		zeta = adja(nvar, t);

		if(conditions==true)
		{
			for(k=0; k<nvar; k++)
			{
				for(l=0; l<nvar; l++)
				{
					if(conditions_okay(k,l)==false)
						okay=false;
				}
			}
		}

		if(okay==false)
		{
			optionsp->out("IMPROPER CONDITIONS on the adjacency matrix!\n");
			optionsp->out("\n");

			optionsp->out("The adjacency matrix to start with \n"
						  "and the conditions onto the graph \n"
						  "do not fit together.\n");

			optionsp->out("\n");

			optionsp->out("The STARTING GRAPH is therefore fixed \n"
						  "to contain (only)m those EDGES \n"
						  "that are GIVEN by the CONDITIONS.\n");

			optionsp->out("\n");
			optionsp->out("\n");
			optionsp->out("\n");





			//*********** adjust zeta to conditions ************
 			zeta = adja(nvar, 0);
			for(k=0; k<nvar; k++)
			{
				for(l=0;l<nvar;l++)
				{
					if(zeta_fix(k,l)==1)
						zeta(k,l)=1;
				}
			}
		}


		for(k=0; k<nvar; k++)
		{
			preg_mods[k]->initialize(zeta,k);
		}
	}





  // FUNCTION: ini_ratio
  // TASK: initializes ratios.
  void FULLCOND_rj::ini_ratio(void)
  {
		acceptance_b = 0;
		acceptance_d = 0;
		acceptance_s = 0;

		nrtrials_b = 0;
		nrtrials_d = 0;
		nrtrials_s = 0;

		step_aborted = true;
  }



  // FUNCTION: ini_hyperpar
  // TASK: initializes hyperparamaters etc
  void FULLCOND_rj::ini_hyperpar(void)
  {
	  sigma_prop=1;

	  if(nobs%2==0)
	  {
		  gamma_nobs1 = log_gamma(0.5*nobs);
		  gamma_nobs2 = log_gamma(0.5*(nobs+1));
	  }
	  else
	  {
		  gamma_nobs2 = log_gamma(0.5*(nobs));
		  gamma_nobs1 = log_gamma(0.5*(nobs+1));
	  }
  }




  // FUNCTION: set_options
  // TASK: sets the options
  void FULLCOND_rj::set_options(void)
  {
	  switch_type = "mix";
      family=" ";
	  limit_number = 10;
	  alpha = 0.05;

	  print_models="normal";

	  conditions=false;
  }






	// FUNCTION: conditions_okay
	// TASK: returns true if conditions are fullfilled
	bool FULLCOND_rj::conditions_okay (unsigned int i, unsigned int j)
	{
		if(		(zeta_fix(i,j)==0 && zeta(i,j)==1)
			||	(zeta_fix(i,j)==1 && zeta(i,j)==0))
			return false;
		else
			return true;
	}





	// FUNCTION: conditions_okay_d
	// TASK: returns true if conditions are fullfilled
	bool FULLCOND_rj::conditions_okay_d (unsigned int i, unsigned int j)
	{
		assert(i!=j);

		if(zeta_fix(i,j) != 1)
			return true;
		else
			return false;
	}


	// FUNCTION: conditions_okay_b
	// TASK: returns true if conditions are fullfilled
	bool FULLCOND_rj::conditions_okay_b (unsigned int i, unsigned int j)
	{
		assert(i!=j);

		if(zeta_fix(i,j)!=0)
			return true;
		else
			return false;
	}


	// FUNCTION: conditions_okay_s
	// TASK: returns true if conditions are fullfilled
	bool FULLCOND_rj::conditions_okay_s (unsigned int i, unsigned int j)
	{
		assert(i!=j);

		if(zeta_fix(i,j)!=0 && zeta_fix(j,i)!=1)
			return true;
		else
			return false;
	}




  // FUNCTION: set_options
  // TASK: sets the options
  void FULLCOND_rj::set_options( unsigned lim, double alph, ST::string switch_t,
									ST::string print_mod, ST::string fix_path)
  {
	  unsigned k,l;

	  switch_type = switch_t;
      print_models = print_mod;
	  limit_number = lim;
	  alpha = alph;

	  conditions=false;

	  if(fix_path != "")
	  {

		  ifstream fin (fix_path.strtochar());
		  zeta_fix.prettyScan(fin);
		  fin.close();

		  if(zeta_fix.cols()==nvar && zeta_fix.rows()==nvar)
		  {
			  for(k=0;k<nvar;k++)
			  {
				  for(l=0; l<nvar;l++)
				  {
					  if(zeta_fix(k,l)==2 && k!=l)
						  conditions=true;
				  }
			  }
		  }

		  if(conditions ==false)
		  {
			  optionsp->out("Improper conditions on the adjacency matrix!");
			  optionsp->out("\n");
			  optionsp->out("Simulation runs without conditions.");
			  optionsp->out("\n");
		  }
	  }
  }








  // FUNCTION: setzeta
  // TASK: initializes zeatmean etc.

  void FULLCOND_rj::setzeta(const adja & zetanew)
  {
	  zetamean = Matrix<double>(zeta);
	  zetameanold = Matrix<double>(zeta);

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
  }


	void FULLCOND_rj::outoptions(void)
	{
		unsigned k,l, value;

		optionsp->out("Type of switch-step: " + switch_type + "\n");
		optionsp->out("Type of starting dag: " + ST::inttostring(type) + "\n");
        optionsp->out("Distribution family: " + preg_mods[0]->get_family()   + "\n");


		if(conditions == true)
		{
			optionsp->out("Conditions are given by: \n");
			optionsp->out("\n");

			for(k=0; k<nvar; k++)
			{
				ST::string s;
				optionsp->out("   ");

				for (l=0; l<nvar; l++)
				{
					value = zeta_fix(k,l);

					if(value==0 || value==1)
						s = s + ST::inttostring(value) + "  ";
					else
						s = s + "*  ";
				}
				optionsp->out(s + "\n");
			}
			optionsp->out("\n");
			optionsp->out("\n");
			optionsp->out("\n");
		}


    }






} //namespace





