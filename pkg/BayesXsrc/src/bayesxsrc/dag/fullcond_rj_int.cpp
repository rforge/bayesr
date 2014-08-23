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





#include "fullcond_rj_int.h"
#include <set>
#include <algorithm>


namespace MCMC
{


FULLCOND_rj_int::FULLCOND_rj_int (vector < FULLCOND_dag_ia* > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp)
				: FULLCOND_rj(o,d,t,r,c,fp)
{

	unsigned i;

	change_preg_mods(dagp);
	ini_structure();

	mixed_case=false;

	for(i=0; i<nvar; i++)
	{
		if(dagp[i]->tell_var_type()=='c')
			mixed_case=true;
	}


}  //constructor_1



FULLCOND_rj_int::FULLCOND_rj_int (vector < FULLCOND_dag_ia_mixed * > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp)
				: FULLCOND_rj(o,d,t,r,c,fp)
{

	unsigned k;

	change_preg_mods(dagp);
	ini_structure();

	mixed_case=false;

	for(k=0; k<nvar; k++)
	{
		if(dagp[k]->tell_var_type()=='c')
			mixed_case=true;
	}



}  //constructor_1a




//constructor_2
FULLCOND_rj_int::FULLCOND_rj_int (ST::string fix, const ST::string & rp, unsigned int lim, double alph ,
								  ST::string switch_t, ST::string print_mod, 
								  unsigned & ty, vector < FULLCOND_dag_ia* > dagp,
								  MCMCoptions * o, const datamatrix & d, const ST::string & t,
								  const unsigned & r, const unsigned & c, const ST::string & fp)
									: FULLCOND_rj(o,d,t,r,c,fp)
{
	//		comment: the following should be done by 
	//		: FULLCOND_rj(lim, alpha, switch_t, print_mod, type,dagp, o,d,t,r,c,fp)
	//		but this is not possible as dagp is of type FULLCOND_dag_d_ia* whereas 
	//		FULLCOND_dag* would be demanded		


	unsigned i;

	setbeta(1,1,0);
	nvar=c;
	nobs=d.rows();

	ini_ratio();
	ini_hyperpar();

	zetamean = datamatrix(nvar,nvar,0);
    nrpar = nvar*nvar;

	set_options( lim, alph, switch_t,print_mod, fix);

    // change dagp from initialized type vector<FULLCOND_dag*> to a vector
	// that contains the FULLCOND_dag_d_ia*
	change_preg_mods(dagp);

	ini_structure(ty);

	mixed_case=false;

	for(i=0; i<nvar; i++)
	{
		if(dagp[i]->tell_var_type()=='c')
			mixed_case=true;
	}

    path_res = rp;


}  //constructor_2




//constructor_2a
FULLCOND_rj_int::FULLCOND_rj_int (ST::string fix, const ST::string & rp, unsigned int lim, double alph,
								  ST::string switch_t,  ST::string print_mod, 
								  unsigned & ty, vector < FULLCOND_dag_ia_mixed* > dagp,
								  MCMCoptions * o, const datamatrix & d, const ST::string & t,
								  const unsigned & r, const unsigned & c, const ST::string & fp)
									: FULLCOND_rj(o,d,t,r,c,fp)
{
	//		comment: the following should be done by 
	//		: FULLCOND_rj(lim, alpha, switch_t, print_mod, type,dagp, o,d,t,r,c,fp)
	//		but this is not possible as dagp is of type FULLCOND_dag_d_ia* whereas 
	//		FULLCOND_dag* would be demanded		


	unsigned i;

	setbeta(1,1,0);
	nvar=c;
	nobs=d.rows();

	ini_ratio();
	ini_hyperpar();

	zetamean = datamatrix(nvar,nvar,0);
    nrpar = nvar*nvar;


	ifstream fin (fix.strtochar());
	zeta_fix.prettyScan(fin);
	fin.close();


	set_options( lim, alph, switch_t, print_mod,fix);

    // change dagp from initialized type vector<FULLCOND_dag*> to a vector
	// that contains the FULLCOND_dag_d_ia*
	change_preg_mods(dagp);

	ini_structure(ty);

	mixed_case=false;

	for(i=0; i<nvar; i++)
	{
		if(dagp[i]->tell_var_type()=='c')
			mixed_case=true;
	}

    path_res = rp;


}  //constructor_2a





 // COPY CONSTRUCTOR

  FULLCOND_rj_int::FULLCOND_rj_int(const FULLCOND_rj_int & fc) : FULLCOND_rj(FULLCOND_rj(fc))
  {
	  
  }



  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_rj_int& FULLCOND_rj_int::operator=(const FULLCOND_rj_int & fc)
  {
	  if (this==&fc)
		  return *this;

     FULLCOND_rj::operator=(FULLCOND_rj(fc));
	  
	  return *this;
  }




  



  // FUNCTION: change_preg_mods
  // TASK: changes preg_mods
  void FULLCOND_rj_int::change_preg_mods( vector <FULLCOND_dag_ia* > dagp)
  {
	  for(unsigned i=0; i<nvar; i++)
		 preg_mods.push_back(dagp[i]);
  }



  // FUNCTION: change_preg_mods
  // TASK: changes preg_mods
  void FULLCOND_rj_int::change_preg_mods( vector <FULLCOND_dag_ia_mixed* > dagp)
  {
	  for(unsigned i=0; i<nvar; i++)
		 preg_mods.push_back(dagp[i]);
  }



// FUNCTION: death_step
// makes death step, tries to delete edge i->j

void FULLCOND_rj_int::death_step(unsigned int i, unsigned int j)
{

	// 1. Bestimme Anzahl der zu entferneneden Regressionskoeffizienten

	unsigned ncoef_old = preg_mods[j]->get_ncoef();
	unsigned num_ia_del = preg_mods[j]->ia_of_i(i);		 //number of to deletig interactions
	unsigned ncoef_new = ncoef_old - 1 - num_ia_del;
	
	if(num_ia_del==0)
	{
		FULLCOND_rj::death_step(i,j);
	}
	else 
	{
		if(mixed_case==true)
				preg_mods[j]->create_matrices("d", ncoef_new);

		// instead of: datamatrix b_new (ncoef_new,1);
		datamatrix & b_new = preg_mods[j]->get_b_new_d();
		// instead of: datamatrix x_new (nobs,ncoef_new);
		datamatrix & x_new = preg_mods[j]->get_x_new_d();
		// instead of: datamatrix xx_new (ncoef_new,ncoef_new);
		datamatrix & xx_new = preg_mods[j]->get_xx_new_d();
		

		// vector of interactions that will vanish
		vector <vector <unsigned > > ias_del;
		preg_mods[j]->ia_of_i( i, ias_del);


		// vector of coefficients which will vanish
		datamatrix beta_old (1+num_ia_del,1);
		
		// current_ia of model j, after successful death-step
		vector <vector <unsigned > > current_ia_n;


		// computing of the new values
		preg_mods[j]->make_new_d_int("d", i, j, num_ia_del, beta_old, 
			current_ia_n, xx_new, b_new, x_new);



		// calculate ratio

		double ratio;
		double log_num;
		double log_denom;

		log_num = preg_mods[j]->calc_SQT_x(x_new, b_new) + preg_mods[j]->calc_SQT_b(b_new);
		log_denom = preg_mods[j]->calc_SQT_x() + preg_mods[j]->calc_SQT_b();

		double prop = preg_mods[j]->get_ln_prop_beta();


		ratio = -1/(2*preg_mods[j]->get_sigma_i()) 
				*(log_num - log_denom) +prop ;

                
		// accept proposal and change the corresponding values

		if(func::accept(ratio) == true)
		{
			zeta(i,j)=0;
			zeta.change_list(i,j,1);
			preg_mods[j]->change_adcol(i,0);

			preg_mods[j]->change(i, b_new, x_new, xx_new, ncoef_new);
			preg_mods[j]->change_occur('d', ias_del);
			preg_mods[j]->change_current('d', ias_del);

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
}







// FUNCTION: birth_step
// makes Birth step

void FULLCOND_rj_int::birth_step(unsigned int v_i, unsigned int v_j)
{
	if(preg_mods[v_j]->get_ncoef()<2)
	{
		FULLCOND_rj::birth_step(v_i,v_j);
	}
	else
	{	
		if(zeta.azy_test(v_i,v_j)== true)
		{
			unsigned k;

			unsigned int ncoef_old = preg_mods[v_j]->get_ncoef();
			unsigned int num_ia_new = preg_mods[v_j]->num_ia_new(v_i);
			unsigned int ncoef_new = ncoef_old + num_ia_new +1; 
			
			/*
			datamatrix b_new (ncoef_new,1);
			datamatrix x_new (nobs,ncoef_new);
			datamatrix xx_new (ncoef_new,ncoef_new);
			*/

			if(mixed_case==true)
				preg_mods[v_j]->create_matrices("b", ncoef_new);
			
			// instead of: datamatrix b_new (ncoef_new,1);
// Vorschlag:
//			datamatrix & b_new = preg_mods[v_j]->get_b_new_b();
			datamatrix b_new = preg_mods[v_j]->get_b_new_b();
			// instead of: datamatrix x_new (nobs,ncoef_new);
// Vorschlag:
//			datamatrix & x_new = preg_mods[v_j]->get_x_new_b() ;
			datamatrix x_new = preg_mods[v_j]->get_x_new_b() ;
			// instead of: datamatrix xx_new (ncoef_new,ncoef_new);
// Vorschlag:
//			datamatrix & xx_new = preg_mods[v_j]->get_xx_new_b() ;
			datamatrix xx_new = preg_mods[v_j]->get_xx_new_b() ;

			datamatrix beta_new (num_ia_new+1,1);
				
			for(k=0; k<num_ia_new+1; k++)
				beta_new(k,0)= 3*rand_normal(); //coefficient which will be added

			
			vector <vector <unsigned > > ia_new; 
			preg_mods[v_j]->new_ia_of_i(v_i, ia_new);


			preg_mods[v_j]->make_new_b_int ("b", v_i, ia_new, beta_new, xx_new, b_new, x_new);


			// calculate ratio

			double log_denom;
			double log_num;
			double prop; 
			double ratio; 

			log_num = preg_mods[v_j]->calc_SQT_x(x_new, b_new)
						+ preg_mods[v_j]->calc_SQT_b(b_new);
			log_denom = preg_mods[v_j]->calc_SQT_x()
						+  preg_mods[v_j]->calc_SQT_b();
			prop = preg_mods[v_j]->get_ln_prop_beta();


			ratio = -1/(2*preg_mods[v_j]->get_sigma_i())
					*(log_num - log_denom)
					- prop;


			//accept proposal and change the corresponding values

			if(func::accept(ratio) == true)
			{
                zeta(v_i,v_j)=1;
				zeta.change_list(v_i,v_j,0);
				zeta.edge_plus();
				preg_mods[v_j]->change_adcol(v_i,1);

				preg_mods[v_j]->change(v_i, b_new, x_new, xx_new, ncoef_new);
				preg_mods[v_j]->change_occur('b', ia_new);
				preg_mods[v_j]->change_current('b', ia_new);

				acceptance_b ++;


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
				/*********************/
			}

			nrtrials_b ++;
			step_aborted = false;

		} // end of "if(azy_test(i,j)== true) ...."
	
    }

} 








void FULLCOND_rj_int::switch_version_1(unsigned v_i, unsigned j)
{
	// ***** properties for DEATH-STEP j->i **********
	unsigned ncoef_old_i = preg_mods[v_i]->get_ncoef();
	unsigned num_ia_del_i = preg_mods[v_i]->ia_of_i(j);		 //number of to deletig interactions
	unsigned ncoef_new_i = ncoef_old_i - 1 - num_ia_del_i;

	// ***** properties for BIRTH-STEP i->j **********
	unsigned int ncoef_old_j = preg_mods[j]->get_ncoef();
	unsigned int num_ia_new_j = preg_mods[j]->num_ia_new(v_i);
	unsigned int ncoef_new_j = ncoef_old_j +1+num_ia_new_j;


	if(num_ia_del_i==0  && preg_mods[j]->get_ncoef()<2)
	{
		FULLCOND_rj::switch_version_1(v_i,j) ;
	}

	else //if(num_ia_del_i == num_ia_new_j)
	{


		//************ DEATH-step for edge j->i  ******************

		if(mixed_case==true)
				preg_mods[v_i]->create_matrices("d", ncoef_new_i);

		// instead of: datamatrix b_new (ncoef_new,1);
		datamatrix & b_new_i = preg_mods[v_i]->get_b_new_d();
		// instead of: datamatrix x_new (nobs,ncoef_new);
		datamatrix & x_new_i = preg_mods[v_i]->get_x_new_d();
		// instead of: datamatrix xx_new (ncoef_new,ncoef_new);
		datamatrix & xx_new_i = preg_mods[v_i]->get_xx_new_d();


		// vector of interactions that will vanish
		vector <vector <unsigned > > ias_del;
		preg_mods[v_i]->ia_of_i( j, ias_del);


		// vector of coefficients which will vanish
		datamatrix beta_old_i (1+num_ia_del_i,1);

		// current_ia of model i, after successful death-step
			vector <vector <unsigned > > current_ia_n_i;

		// computing of the new values
		preg_mods[v_i]->make_new_d_int("s", j, v_i, num_ia_del_i,
			beta_old_i, current_ia_n_i, xx_new_i, b_new_i, x_new_i);


		//************ERGAENZUNG ANFANG*************
		datamatrix sig_mean_i (ncoef_new_i,ncoef_new_i);
		sig_mean_i.assign(xx_new_i.inverse());

		datamatrix mean_i (ncoef_new_i,1);
		mean_i.mult(sig_mean_i, x_new_i.transposed()*preg_mods[v_i]->get_y());

        double sigma_new_i = sample_sigma('i', v_i, ncoef_new_i, mean_i, x_new_i);
		b_new_i.mult(sig_mean_i.root(),rand_normvek(ncoef_new_i));
		b_new_i.plus(b_new_i,mean_i);

		preg_mods[v_i]->calc_lin_prop(x_new_i, b_new_i);
		//********ERGAENZUNG ENDE *******************







		//************ BIRTH-step for edge i->j  ******************


		//datamatrix b_new (ncoef_new,1);
		//datamatrix x_new (nobs,ncoef_new);
		//datamatrix xx_new (ncoef_new,ncoef_new);
		

		if(mixed_case==true)
				preg_mods[j]->create_matrices("b", ncoef_new_j);

		// instead of: datamatrix b_new (ncoef_new,1);
// Vorschlag:
//		datamatrix & b_new_j = preg_mods[j]->get_b_new_b();
		datamatrix b_new_j = preg_mods[j]->get_b_new_b();
		// instead of: datamatrix x_new (nobs,ncoef_new);
// Vorschlag:
//		datamatrix & x_new_j = preg_mods[j]->get_x_new_b() ;
		datamatrix x_new_j = preg_mods[j]->get_x_new_b() ;
		// instead of: datamatrix xx_new (ncoef_new,ncoef_new);
// Vorschlag:
//		datamatrix & xx_new_j = preg_mods[j]->get_xx_new_b() ;
		datamatrix xx_new_j = preg_mods[j]->get_xx_new_b() ;

		datamatrix beta_new_j (num_ia_new_j+1,1);

		vector <vector <unsigned > > ia_new_j;
		preg_mods[j]->new_ia_of_i(v_i, ia_new_j);

		preg_mods[j]->make_new_b_int ("s", v_i, ia_new_j, beta_new_j, xx_new_j, b_new_j, x_new_j);


		//************ERGAENZUNG ANFANG*************
		datamatrix sig_mean_j (ncoef_new_j,ncoef_new_j);
		sig_mean_j.assign(xx_new_j.inverse());
		datamatrix mean_j (ncoef_new_j,1);

		mean_j.mult(sig_mean_j, x_new_j.transposed()*preg_mods[j]->get_y());

        double sigma_new_j = sample_sigma('j', j, ncoef_new_j, mean_j, x_new_j);
		b_new_j.mult(sig_mean_j.root(),rand_normvek(ncoef_new_j));
		b_new_j.plus(b_new_j,mean_j);

		preg_mods[j]->calc_lin_prop(x_new_j, b_new_j);
        //************ERGAENZUNG ENDE **************



		// calculate ratio
		double ratio =  ratio_s( v_i, j, b_new_i, b_new_j, x_new_i, x_new_j, mean_i, mean_j,
								     xx_new_i, xx_new_j, sig_mean_i, sig_mean_j,
									 sigma_new_i, sigma_new_j);

        //ratio_s( i, j, b_new_i, b_new_j, x_new_i, x_new_j, mean_i, mean_j, sigma_new_i, sigma_new_j);

		if(func::accept(ratio) == true)
		{
			zeta(v_i,j)=1;
			zeta(j,v_i)=0;
			zeta.change_list(v_i,j,2);

			preg_mods[j]->change_adcol(v_i,1);
			preg_mods[v_i]->change_adcol(j,0);

			preg_mods[j]->change(v_i, b_new_j, x_new_j, xx_new_j, ncoef_new_j);
			preg_mods[v_i]->change(j, b_new_i, x_new_i, xx_new_i, ncoef_new_i);

			preg_mods[j]->change_occur('b', ia_new_j);
			preg_mods[v_i]->change_occur('d', ias_del);

			preg_mods[j]->change_current('b', ia_new_j);
			preg_mods[v_i]->change_current('d', ias_del);

			acceptance_s ++;
		}
	}
}






 void FULLCOND_rj_int::switch_version_2(unsigned v_i, unsigned v_j)
{
	if(randnumbers::uniform()<0.5)
	{

		if(preg_mods[v_i]->tell_var_type()=='c')
		{
			assert(preg_mods[v_i]->ia_of_i(v_j)==0);
			assert(preg_mods[v_j]->tell_var_type()=='d');
		}

		if(preg_mods[v_j]->tell_var_type()=='c')
		{
			assert(preg_mods[v_j]->num_ia_new(v_i)==0);
			assert(preg_mods[v_i]->tell_var_type()=='d');
		}


			

        // ***** properties for DEATH-STEP j->i **********
	    unsigned ncoef_old_i = preg_mods[v_i]->get_ncoef();
	    unsigned num_ia_del_i = preg_mods[v_i]->ia_of_i(v_j);		 //number of to deleting  interactions
	    unsigned ncoef_new_i = ncoef_old_i - 1 - num_ia_del_i;

	    // ***** properties for BIRTH-STEP i->j **********
	    unsigned int ncoef_old_j = preg_mods[v_j]->get_ncoef();
	    unsigned int num_ia_new_j = preg_mods[v_j]->num_ia_new(v_i);
		unsigned int ncoef_new_j = ncoef_old_j + num_ia_new_j +1;


	    if(preg_mods[v_i]->get_ncoef()<3  && preg_mods[v_j]->get_ncoef()<2)
	    {
	 	    FULLCOND_rj::switch_version_2(v_i,v_j) ;
        }

		
		// eigentlich else ohne bedingung!!!

     	else
	    {

		    //************ DEATH-step for edge j->i  ******************

			if(mixed_case==true)
				preg_mods[v_i]->create_matrices("d", ncoef_new_i);

            // instead of: datamatrix b_new (ncoef_new,1);
		    datamatrix & b_new_i = preg_mods[v_i]->get_b_new_d();
		    // instead of: datamatrix x_new (nobs,ncoef_new);
		    datamatrix & x_new_i = preg_mods[v_i]->get_x_new_d();
		    // instead of: datamatrix xx_new (ncoef_new,ncoef_new);
		    datamatrix & xx_new_i = preg_mods[v_i]->get_xx_new_d();

		    // vector of interactions that will vanish
		    vector <vector <unsigned > > ias_del;
		    preg_mods[v_i]->ia_of_i( v_j, ias_del);

		    // vector of coefficients which will vanish
		    datamatrix beta_old_i (1+num_ia_del_i,1);

		    // current_ia of model i, after successful death-step
			vector <vector <unsigned > > current_ia_n_i;

		    // computing of the new values
		    preg_mods[v_i]->make_new_d_int("s", v_j, v_i, num_ia_del_i,
                beta_old_i, current_ia_n_i, xx_new_i, b_new_i, x_new_i);



			//************ BIRTH-step for edge i->j ******************

			if(mixed_case==true)
					preg_mods[v_j]->create_matrices("b", ncoef_new_j);

			// instead of: datamatrix b_new (ncoef_new,1);
// Vorschlag:
//			datamatrix & b_new_j = preg_mods[v_j]->get_b_new_b();
			datamatrix b_new_j = preg_mods[v_j]->get_b_new_b();
			// instead of: datamatrix x_new (nobs,ncoef_new);
// Vorschlag:
//			datamatrix & x_new_j = preg_mods[v_j]->get_x_new_b() ;
			datamatrix x_new_j = preg_mods[v_j]->get_x_new_b() ;
			// instead of: datamatrix xx_new (ncoef_new,ncoef_new);
// Vorschlag:
//			datamatrix & xx_new_j = preg_mods[v_j]->get_xx_new_b() ;
			datamatrix xx_new_j = preg_mods[v_j]->get_xx_new_b() ;

			datamatrix beta_new_j (num_ia_new_j+1,1);

			vector <vector <unsigned > > ia_new_j;

			preg_mods[v_j]->new_ia_of_i(v_i, ia_new_j);

			preg_mods[v_j]->make_new_b_int ("s", v_i, ia_new_j, beta_new_j, xx_new_j, b_new_j, x_new_j);


			zeta(v_i,v_j)=1;
			zeta(v_j,v_i)=0;
			zeta.change_list(v_i,v_j,2);

			preg_mods[v_j]->change_adcol(v_i,1);
			preg_mods[v_i]->change_adcol(v_j,0);

			preg_mods[v_j]->change(v_i, b_new_j, x_new_j, xx_new_j, ncoef_new_j);
			preg_mods[v_i]->change(v_j, b_new_i, x_new_i, xx_new_i, ncoef_new_i);

			preg_mods[v_j]->change_occur('b', ia_new_j);
			preg_mods[v_i]->change_occur('d', ias_del);

			preg_mods[v_j]->change_current('b', ia_new_j);
			preg_mods[v_i]->change_current('d', ias_del);

			acceptance_s ++;
		}
    }
	
	step_aborted = false;
}





// FUNCTION: ratio_s
 // TARGET: computes acceptance ratio in the switch step
 double FULLCOND_rj_int::ratio_s_int(unsigned int i,unsigned int j, 
							const datamatrix & b_new_i, const datamatrix & b_new_j, 
							const datamatrix & x_new_i, const datamatrix & x_new_j)
 {
	 double log_num_j, log_num_i, log_denom_j, log_denom_i;
	 double	prop_j, prop_i;


	 log_num_j = preg_mods[j]->calc_SQT_x(x_new_j, b_new_j) 
					+ preg_mods[j]->calc_SQT_b(b_new_j);

	 log_denom_j = preg_mods[j]->calc_SQT_x()
					+  preg_mods[j]->calc_SQT_b();

	 log_num_i = preg_mods[i]->calc_SQT_x(x_new_i, b_new_i) 
				+ preg_mods[i]->calc_SQT_b(b_new_i);

	 log_denom_i = preg_mods[i]->calc_SQT_x() 
				+ preg_mods[i]->calc_SQT_b();

	prop_j = preg_mods[j]->get_ln_prop_beta();
	prop_i = preg_mods[i]->get_ln_prop_beta();



	double log_ratio = -1/(2*preg_mods[i]->get_sigma_i()) 
						*(log_num_i - log_denom_i )
					   - 1/(2*preg_mods[j]->get_sigma_i()) 
					    *(log_num_j - log_denom_j)
						+ prop_i - prop_j ;

	return log_ratio;
 }



 // FUNCTION: make_new_b 
  // TASK: computes the new values for a birth-step
  void FULLCOND_rj_int::make_new_b (ST::string step, unsigned int i, unsigned int j, unsigned ia_new,
						datamatrix & beta_new, datamatrix & xx_new, datamatrix & b_new, 
						datamatrix & x_new)
  {
	  void make_new_b(vector<unsigned> ia_new, double beta_new, 
				datamatrix & xx_new, datamatrix & b_new, datamatrix & x_new);
  }



  // FUNCTION: switch_step
// makes switch step from j->i to i->j
void FULLCOND_rj_int::switch_step(unsigned int v_i, unsigned int v_j)
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
			FULLCOND_rj_int::switch_version_2(v_i,v_j);
		}
		else if (switch_type=="normal")
		{
			FULLCOND_rj_int::switch_version_1(v_i,v_j);
		}
		else if (switch_type=="mix")
		{
			if(zeta.equi_test(v_j,v_i)==true)
			{
				FULLCOND_rj_int::switch_version_2(v_i,v_j);
			}
			else
			{
				FULLCOND_rj_int::switch_version_1(v_i,v_j);
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


} //namespace








