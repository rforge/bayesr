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





#include "fullcond_rj_mix.h"
#include <set>
#include <algorithm>


namespace MCMC
{


FULLCOND_rj_mix::FULLCOND_rj_mix (vector < FULLCOND_dag_ia * > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp)
				: FULLCOND_rj_int(dagp, o,d,t,r,c,fp)
{

}  //constructor_1



FULLCOND_rj_mix::FULLCOND_rj_mix (vector < FULLCOND_dag_ia_mixed * > dagp,
				MCMCoptions * o, const datamatrix & d, const ST::string & t,
				const unsigned & r, const unsigned & c, const ST::string & fp)
				: FULLCOND_rj_int(dagp, o,d,t,r,c,fp)
{

}  //constructor_1



FULLCOND_rj_mix::FULLCOND_rj_mix (ST::string fix, const ST::string & rp, unsigned int lim,
								  double alph , ST::string switch_t, ST::string print_mod,
								  unsigned & ty, vector < FULLCOND_dag_ia* > dagp,
								  MCMCoptions * o, const datamatrix & d, const ST::string & t,
								  const unsigned & r, const unsigned & c, const ST::string & fp)
									: FULLCOND_rj_int(fix, rp, lim, alph, switch_t,print_mod,ty,
																		dagp,o,d,t,r,c,fp)
{
}  //constructor_2


FULLCOND_rj_mix::FULLCOND_rj_mix (ST::string fix, const ST::string & rp, unsigned int lim,
								  double alph , ST::string switch_t, ST::string print_mod, 
								  unsigned & ty, vector < FULLCOND_dag_ia_mixed * > dagp,
								  MCMCoptions * o, const datamatrix & d, const ST::string & t,
								  const unsigned & r, const unsigned & c, const ST::string & fp)
									: FULLCOND_rj_int(fix, rp, lim, alph, switch_t,print_mod,ty,
																		dagp,o,d,t,r,c,fp)
{
}  //constructor_2



 // COPY CONSTRUCTOR

  FULLCOND_rj_mix::FULLCOND_rj_mix(const FULLCOND_rj_mix & fc) 
									: FULLCOND_rj_int(FULLCOND_rj_int(fc))
  {
  }
  


  // OVERLOADED ASSIGNMENT OPERATOR

  const FULLCOND_rj_mix& FULLCOND_rj_mix::operator=(const FULLCOND_rj_mix & fc)
  {
	  if (this==&fc)
		  return *this;

     FULLCOND_rj_int::operator=(FULLCOND_rj_int(fc));
	  
	  return *this;
  }




  	void FULLCOND_rj_mix::update(void)
    {
		  step_aborted = true;


		/*********	in principle like FULLCOND_rj::rj_step, but now 
					depending on the var_type functions of 
					FULLCOND_rj or FULLCOND_rj_int are used **********************/
		
		  while (step_aborted==true)
		  {
				// two variables are randomly chosen
				unsigned int i,j;

				i= rand() % nvar;
				j=i;

				while(i==j)
				j= rand() % nvar;


				// decide which kind of step (birth, death, switch)
				if (zeta(i,j) == 1)
				{
					 if(preg_mods[j]->tell_var_type()=='d')
					 {
						 if(conditions == false)
							FULLCOND_rj_int::death_step(i,j);
						 else
							if(conditions_okay_d(i,j) ==true)
								FULLCOND_rj_int::death_step(i,j);
							

					 }
					 else
					 {
						 if(conditions == false)
							FULLCOND_rj::death_step(i,j);
						 else
							if(conditions_okay_d(i,j) ==true)
							
								FULLCOND_rj::death_step(i,j);

					 }
				}
				else if (zeta(j,i)==1)	
				{
					if(preg_mods[j]->tell_var_type()=='d' ||
						preg_mods[i]->tell_var_type()=='d')
					{
						if(conditions == false)
							FULLCOND_rj_int::switch_step(i,j);
						
						 else
							if(conditions_okay_s(i,j) ==true)
							{
								FULLCOND_rj_int::switch_step(i,j);
								//cout<<"s_int: "<<i<<" "<<j<<endl;
							}
						 
					}
					else  if(preg_mods[j]->tell_var_type()=='c' && 
								preg_mods[i]->tell_var_type()=='c')
					{
                          if(conditions == false)
							FULLCOND_rj::switch_step(i,j);
						  else
							if(conditions_okay_s(i,j) ==true)
								FULLCOND_rj::switch_step(i,j);
					}
				}
				else
				{
					if(preg_mods[j]->tell_var_type()=='d')
					{
						 if(conditions == false)
							FULLCOND_rj_int::birth_step(i,j);
						  else
							if(conditions_okay_b(i,j) ==true)
								FULLCOND_rj_int::birth_step(i,j);
					}
					else
					{
						if(conditions == false)
							FULLCOND_rj::birth_step(i,j);
						  else
							if(conditions_okay_b(i,j) ==true)
								FULLCOND_rj::birth_step(i,j);
					}
				}
		}



		  /************* from this point on like in Fullcond_rj::update ***************/
		  nrtrials++;
		  update_zeta();


		  
		  //std::set < ST::string, unsigned int> freq;

		  if((optionsp->get_nriter() > optionsp->get_burnin()) &&
			  (optionsp->get_nriter() % (optionsp->get_step()) == 0))
		  {
			  store_model();
		  }

	}





} //namespace








