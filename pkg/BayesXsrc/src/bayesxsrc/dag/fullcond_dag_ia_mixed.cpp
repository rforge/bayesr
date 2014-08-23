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





#include "fullcond_dag_ia_mixed.h"
#include <algorithm>
#include <iterator>


namespace MCMC                                              
{

	


// constructor 3;  used in mixed case, corresponds to constructor1 in discrete case

	FULLCOND_dag_ia_mixed::FULLCOND_dag_ia_mixed (IA * iap, double s_i, unsigned int num,
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp)
							: FULLCOND_dag_ia(iap, s_i, num, o,d,t,r,c,fp)
	{
		unsigned i;


		
		// *********** check, if which ias are allowed **************
		for(i=0; i<occurrence.size(); i++)
		{
			if(pia->ia_okay(i)==false)
				occurrence[i] =-1;	
		}


		num_continuous_pa = 0;			
	    num_discrete_pa = 0;

		
		pia->give_var_kind (adcol,num_continuous_pa, num_discrete_pa);
  
	}



	// constructor 4;  used in mixed case, corresponds to constructor2 in discrete case

	FULLCOND_dag_ia_mixed::FULLCOND_dag_ia_mixed ( bool detail_ia, IA * iap, double v_a, double v_b,
							ST::string prio_s, 
							bool d_all, const datamatrix & res, double s_i, unsigned int num, 
							MCMCoptions * o, const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c, const ST::string & fp)
							: FULLCOND_dag_ia(detail_ia, iap, v_a, v_b, prio_s, d_all, res, s_i,num,
											o,d,t,r,c,fp)
	{
        unsigned i;
        
		// *********** check, if which ias are allowed **************
		for(i=0; i<occurrence.size(); i++)
		{
			if(pia->ia_okay(i)==false)
				occurrence[i] = -1;	
		}

		num_continuous_pa = 0;			
	    num_discrete_pa = 0;


		//******* hier falsch, da adcol nicht bekannt
		pia->give_var_kind (adcol,num_continuous_pa, num_discrete_pa);

	}



	ST::string FULLCOND_dag_ia_mixed::get_family(void)
	{
		return "mixed";
	}




	// FUNCTION initialize
	// TASK: initializes x and xx  and y for pred_mod[i] (regression model i)
	void FULLCOND_dag_ia_mixed::initialize (const adja & zeta, unsigned int j)
	{
		//pia->give_var_kind (adcol,num_continuous_pa, num_discrete_pa);

		FULLCOND_dag::initialize (zeta,j);
	}






	// FUNCTION: ini_ia
	// TASK: creates and adds the interactions of the given main effects
	void FULLCOND_dag_ia_mixed::initialize_ia (const adja & zeta, unsigned int j)
	{
		unsigned k,l;
		vector <unsigned> new_ia(2);


		pia->give_var_kind (adcol,num_continuous_pa, num_discrete_pa);

		unsigned not_allowed = num_continuous_pa *(num_continuous_pa-1)/2;

		ncoef_ia = (ncoef_m-1)*(ncoef_m-2)/2 - not_allowed;
		ncoef=ncoef_m+ncoef_ia;

		for(k=0; k<nvar; k++)
		{	
			if(zeta(k,j)==1)
			{
				for(l=k+1; l<nvar; l++)
				{
					if(zeta(l,j)==1)
					{
						new_ia[0]=k;
						new_ia[1]=l;

						if(pia->ia_okay(new_ia)==true) 	//interaction is an allowed one
						{

							if(pia->already_there(new_ia)==false)
								pia->make_ia(new_ia);		// new ia_term is created and added

							current_ia.push_back(new_ia);
							change_occur('b', new_ia);
						}
					}
				}
			}
		}
	}






	 



	
	
	// FUNCTION: ia_of_i
	// TASK: counts the number of interactions containing variable i
	unsigned FULLCOND_dag_ia_mixed::ia_of_i(unsigned i)
	{
		unsigned count;
		count=0;
       /*
        if( var_type=='d')
        {
            std::vector< vector< unsigned> > :: iterator it_cur;
            std::vector<int> :: iterator it_occ;

            it_cur = current_ia.begin();
            it_occ = occurrence.begin();

            for(k=0; k<ncoef_ia; k++, ++it_cur, ++it_occ)
            {
			         if(
				     ((*it_cur)[0]==i ||	(*it_cur)[1]==i )
                     &&
                     *it_occ != -1)
				     count++;
            }

            if(pia->tell_var_type(i)=='d')
               assert(count==ncoef_m-2)    ;
            else
               assert(count==num_discrete_pa);
        }
        */

        if(var_type=='d')
        {
           if(pia->tell_var_type(i)=='d' )
               count = ncoef_m-2;
           else
               count= num_discrete_pa;
        }
        else
            count=0;
            

		return count; 
	}










	// FUNCTION: new_ia_of_i
	// TASK:adds ALLOWED interactions containing variable i to v
	// when i IS NOT already a main effect
	void FULLCOND_dag_ia_mixed::new_ia_of_i( unsigned i, vector <vector <unsigned > > & v)
	{
		assert(v.size()==0);
        assert(i !=self);
		unsigned k,l;

		l=0;
		if(var_type == 'd')
		{
			for(k=0; k<nvar, l<ncoef_m-1; k++)
			{		
				if(adcol(k,0)==1)
				{
					vector <unsigned> ia;
					if(k<i)
					{
						ia.push_back(k);
						ia.push_back(i);
					}
					else
					{
						ia.push_back(i);
						ia.push_back(k);
					}
					
					if (pia->ia_okay(ia) == true)
						v.push_back(ia);
				
					l++;
				}
			}
		}
	}
 




  // FUNCTION: num_ia_of_i
  // TASK: returns the number of allowed interactions of the existing main effect i
  unsigned FULLCOND_dag_ia_mixed::num_ia_of_i(unsigned i)
  { 

	  if(var_type == 'd')
	  {
		  if(pia->tell_var_type(i) == 'd')
			  return ncoef_m-2;
		  else
			  return num_discrete_pa;

	  }
	  else
		  return 0;
  }


   // FUNCTION: num_ia_new
   // TASK: returns the number of allowed new interactions of the new main effect i
   unsigned FULLCOND_dag_ia_mixed::num_ia_new(unsigned i)
   {
	  
	   if(var_type=='d')
	   {
		   if(pia->tell_var_type(i) == 'd')
			  return ncoef_m-1;
		   else
			  return num_discrete_pa;
	   }
	   else
		   return 0;
   }





   void FULLCOND_dag_ia_mixed::change(unsigned i, const datamatrix & beta_help_new, const datamatrix & x_new, 
								const datamatrix & xx_new, unsigned int ncoef_new)
   {

	   
	   if(ncoef_new<ncoef) //death
	   {
		   if(pia->tell_var_type(i) == 'd')
			   num_discrete_pa--;
		   else
			   num_continuous_pa--;
	   }
	   else
	   {
		   if(pia->tell_var_type(i) == 'd')
			   num_discrete_pa++;
		   else
			   num_continuous_pa++;
	   }
		  


	   FULLCOND_dag::change(i, beta_help_new, x_new, xx_new, ncoef_new);
   }

   




	} // namespace MCMC



