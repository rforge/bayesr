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





#include "fullcond_dag_ia.h"
#include <algorithm>
#include <iterator>



namespace MCMC                                              
{

	// constructor 1;  used in discrete case

	FULLCOND_dag_ia::FULLCOND_dag_ia (IA * iap, double s_i, unsigned int num, 
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp)
							: FULLCOND_dag_d(s_i, num, o,d,t,r,c,fp)
	{

		
		unsigned i,j;
		detail=false;

		pia = iap;

		var_type = pia->tell_var_type(self);		 
		mixed_case = pia->tell_mixed_case();  
		max_ia_order = pia->get_max_order();

		all_possible_ia = nvar*(nvar-1)/2;

	
		if(detail==true)
		{
			nrpar = nvar+all_possible_ia; 
		}


		
		occurrence = vector <int> (all_possible_ia,0);
		occmean = datamatrix(all_possible_ia,1,0);

		y_ia = datamatrix(nobs,1);
		x_ia_b = datamatrix(1, 1);
		xx_ia_b= datamatrix(1, 1);
		x_ia_d = datamatrix(1, 1);
		xx_ia_d= datamatrix(1, 1);

		ia_b_there =false;
		ia_d_there = false;

		proposal_version=2;

		// *********** compute all_ia ********************************************
		vector <unsigned> help (2);
		for(i=0; i<nvar; i++)
		{
			for(j=i+1; j<nvar; j++)
			{
				help[0] = i;
				help[1] = j;

				all_ia.push_back(help);
			}
		}
	}


	// constructor 2;  used in dicrete case
	FULLCOND_dag_ia::FULLCOND_dag_ia (bool detail_ia, IA * iap, double v_a, double v_b, ST::string prio_s,
							bool d_all, const datamatrix & res, double s_i, unsigned int num, 
							MCMCoptions * o, const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c, const ST::string & fp)
							: FULLCOND_dag_d(v_a, v_b, prio_s, d_all, res, s_i,num, 
											o,d,t,r,c,fp)
	{
		unsigned i,j;
		   

		detail=detail_ia;
		pia = iap;

		var_type = pia->tell_var_type(self);		 
		mixed_case = pia->tell_mixed_case();  
		max_ia_order = pia->get_max_order();

		all_possible_ia = nvar*(nvar-1)/2;

		
		if(detail==true)
		{
			nrpar = nvar+all_possible_ia; 
			setbeta(nvar+all_possible_ia,1,1);
		}


		
		occurrence = vector <int> (all_possible_ia,0);
		occmean = datamatrix(all_possible_ia,1,0);

		
		y_ia = datamatrix(nobs,1);
		x_ia_b = datamatrix(1, 1);
		xx_ia_b= datamatrix(1, 1);
		x_ia_d = datamatrix(1, 1);
		xx_ia_d= datamatrix(1, 1);

		ia_b_there =false;
		ia_d_there = false;

        proposal_version=2; 

		// *********** compute all_ia ******************************************************
		vector <unsigned> help (2);
		for(i=0; i<nvar-1; i++)
		{
			for(j=i+1; j<nvar; j++)
			{
				help[0] = i;
				help[1] = j;

				all_ia.push_back(help);
			}
		}
	}




// constructor 3;  used in mixed case, corresponds to constructor1 in discrete case

	FULLCOND_dag_ia::FULLCOND_dag_ia (bool detail_ia, char typ, IA * iap, double s_i, unsigned int num,
							MCMCoptions * o,
							const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c,
							const ST::string & fp)
							: FULLCOND_dag_d(s_i, num, o,d,t,r,c,fp)
	{
		var_type = typ;		// discrete variables
		mixed_case = true;  // variables can be dicrete OR continuous

		unsigned i,j;
		detail=true;

		pia = iap;
		max_ia_order = pia->get_max_order();

		all_possible_ia = nvar*(nvar-1)/2;
		
		if(detail==true)
		{
			nrpar = nvar+all_possible_ia; 
			setbeta(nvar+all_possible_ia,1,0);
		}


		
		occurrence = vector <int> (all_possible_ia,0);
		occmean = datamatrix(all_possible_ia,1,0);

		y_ia = datamatrix(nobs,1);
		x_ia_b = datamatrix(1, 1);
		xx_ia_b= datamatrix(1, 1);
		x_ia_d = datamatrix(1, 1);
		xx_ia_d= datamatrix(1, 1);

		ia_b_there =false;
		ia_d_there = false;

		proposal_version=2;

		// *********** compute all_ia ********************************************
		vector <unsigned> help (2);
		for(i=0; i<nvar; i++)
		{
			for(j=i+1; j<nvar; j++)
			{
				help[0] = i;
				help[1] = j;

				all_ia.push_back(help);
			}
		}
	}



	// constructor 4;  used in mixed case, corresponds to constructor2 in discrete case

	FULLCOND_dag_ia::FULLCOND_dag_ia (bool detail_ia, char typ, IA * iap, double v_a, double v_b, ST::string prio_s,
							bool d_all, const datamatrix & res, double s_i, unsigned int num,
							MCMCoptions * o, const datamatrix & d, const ST::string & t,
							const unsigned & r, const unsigned & c, const ST::string & fp)
							: FULLCOND_dag_d(v_a, v_b, prio_s, d_all, res, s_i,num,
											o,d,t,r,c,fp)
	{
		var_type = typ;
		mixed_case = true;  // variables can be dicrete OR continuous

		unsigned i,j;

		detail=false;
		pia = iap;
		max_ia_order = pia->get_max_order();

		all_possible_ia = nvar*(nvar-1)/2;

		
		if(detail==true)
		{
			nrpar = nvar+all_possible_ia;
			setbeta(nvar+all_possible_ia,1,1);
		}



		occurrence = vector <int> (all_possible_ia,0);
		occmean = datamatrix(all_possible_ia,1,0);


		y_ia = datamatrix(nobs,1);
		x_ia_b = datamatrix(1, 1);
		xx_ia_b= datamatrix(1, 1);
		x_ia_d = datamatrix(1, 1);
		xx_ia_d= datamatrix(1, 1);

		ia_b_there =false;
		ia_d_there = false;

        proposal_version=2;

		// *********** compute all_ia ******************************************************
		vector <unsigned> help (2);
		for(i=0; i<nvar-1; i++)
		{
			for(j=i+1; j<nvar; j++)
			{
				help[0] = i;
				help[1] = j;

				all_ia.push_back(help);
			}
		}
	}




	  ST::string FULLCOND_dag_ia::get_family(void)
	  {
		  return "binary with interactions";
	  }






	// FUNCTION: ini_ia
	// TASK: creates and adds the interactions of the given main effects
	void FULLCOND_dag_ia::initialize_ia (const adja & zeta, unsigned int j)
	{
		unsigned k,l;
		vector <unsigned> new_ia(2);

		ncoef_ia = (ncoef_m-1)*(ncoef_m-2)/2;
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

						if(pia->already_there(new_ia)==false)
							pia->make_ia(new_ia);		// new ia_term is created and added

						current_ia.push_back(new_ia);
						change_occur('b', new_ia);
					}
				}
			}
		}
	}



	// FUNCTION: write_ia_to_x
	// TASK: adds the interaction to x when x is initialized
	// is called in FULLCOND_dag::initialize
	void FULLCOND_dag_ia::write_ia_to_x(void)
	{

		unsigned k,l;
		double * workx;
		double * workia;

		vector<unsigned> ia_new (2);

		for(k=0; k<ncoef_ia; k++)
		{
			ia_new = current_ia[k];




			/**************************
			if((current_ia[k]).size()==0)
			{
				unsigned stop=0;
			}

			cout<<k<<": "
				<<ia_new[0]
				<<ia_new[1]<<endl;
			****************************/

			workx = getV_x()+ncoef_m+k;
			workia = pia->get_ia(ia_new);

			for(l=0; l<nobs; l++, workia++)
			{
				*workx = *workia;
				workx = workx + ncoef;
			}
		}
	}


	void FULLCOND_dag_ia::create_matrices (void)
	{
		if(mixed_case==false)
		{
		
				unsigned num;
			// parameters for a birth-step (one dimension more)
			// if(ncoef-ncoef_ia<nvar)
			{
				b_new_b = datamatrix(ncoef+ncoef_m,1,0);
				x_new_b  = datamatrix(nobs,ncoef+ncoef_m,0);
				xx_new_b = datamatrix(ncoef+ncoef_m,ncoef+ncoef_m,0);
			}

			// parameters for a death-step (one dimension less)
			if(ncoef>1)
			{
				num = ncoef_m-1+(ncoef_m-2)*(ncoef_m-3) /2;

				b_new_d  = datamatrix(num,1,0);
				x_new_d  = datamatrix(nobs,num,0);
				xx_new_d = datamatrix(num,num,0);
			}
		}
	}
	 


	
	
	// FUNCTION: ia_of_i
	// TASK: counts the number of interactions containing variable i
	unsigned FULLCOND_dag_ia::ia_of_i(unsigned i)
	{
		unsigned count,k;
		count=0;

		std::vector< vector< unsigned> > :: iterator it_cur;
		it_cur = current_ia.begin();

		for(k=0; k<ncoef_ia; k++, ++it_cur)
		{
			if(		(*it_cur)[0]==i 
				||	(*it_cur)[1]==i)
				count++;
		}

		return count; 
	}


	// FUNCTION: ia_of_i
	// TASK:adds interactions containing variable i to v
	void FULLCOND_dag_ia::ia_of_i( unsigned i, vector <vector <unsigned > > & v)
	{
		assert(v.size()==0);
		unsigned k;

		std::vector< vector< unsigned> > :: iterator it_cur;

		it_cur = current_ia.begin();

		for(k=0; k<ncoef_ia; k++, ++it_cur)
		{
			if(		(*it_cur)[0]==i 
				||	(*it_cur)[1]==i)
				v.push_back(*it_cur);
		}
	}


	// FUNCTION: ia_of_i
	// TASK:adds interactions containing variable i to v
	// when i IS NOT already a main effect
	void FULLCOND_dag_ia::new_ia_of_i( unsigned i, vector <vector <unsigned > > & v)
	{
		assert(v.size()==0);
        assert(i !=self);
		unsigned k,l;

		l=0;
		for(k=0; k<nvar, l<ncoef_m-1; k++)
		{
			if(adcol(k,0)==1)
			{
				assert(k!=i);

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
				
				v.push_back(ia);
				l++;
			}
		}
	}
 




  // FUNCTION: ia_of_i
  // TASK: returns the number of interactions of the main effect i
  // if an interaction between all pairs of main effects is assumed
  unsigned FULLCOND_dag_ia::ia_of_i(void)
  {
	  return ncoef_m-2;
  }





  








	void FULLCOND_dag_ia::update(void)
	{
		if(var_type=='d')
		{
			FULLCOND_dag_d::update();	
		}
		else
			FULLCOND_dag::update();	


		if(detail==true)
			write_to_beta_ia();
		update_occ();	
	}



	unsigned FULLCOND_dag_ia::get_nvar(void)
	{
		return nvar;
	}






	void FULLCOND_dag_ia::birth_step (vector<unsigned> ia_new)
	{
		unsigned int ncoef_old = ncoef;
		unsigned int ncoef_new = ncoef_old + 1;


		datamatrix b_new;
		datamatrix x_new;
		datamatrix xx_new;

		if(mixed_case==false)
		{
			// instead of: datamatrix b_new (ncoef_new,1);
			datamatrix  b_new = get_b_new_b();
			// instead of: datamatrix x_new (nobs,ncoef_new);
			datamatrix  x_new = get_x_new_b() ;
			// instead of: datamatrix xx_new (ncoef_new,ncoef_new);
			datamatrix  xx_new = get_xx_new_b() ;
		}

		double beta_new = rand_normal(); //coefficient which will be added


		// ***** computing of the new values *****************************************
		make_new_b (ia_new, beta_new, xx_new,  b_new, x_new);


		
		// *****  calculate ratio ****************************************************
		double log_denom;
		double log_num;
		double ratio; 

		log_num = calc_SQT_x(x_new, b_new) + calc_SQT_b(b_new);
		log_denom = calc_SQT_x() +  calc_SQT_b();

		ratio = -1/(2* sigma_i) 
			*(log_num - log_denom)  - p_prop(beta_new) ;


		// *****  accept proposal and change the corresponding values ****************
		if(func::accept(ratio) == true)
		{
			change_current('b', ia_new);
			change_occur('b', ia_new);
			change(99, b_new, x_new, xx_new, ncoef_new);  
		}
 

// ************************************************************************************


	}




	// FUNCTION: make_new_b
	// TASK: computes the new values for a birth-step
	void FULLCOND_dag_ia::make_new_b (vector<unsigned> ia_new, double beta_new, 
				datamatrix & xx_new, datamatrix & b_new, datamatrix & x_new)
	{
		unsigned ncoef_new = ncoef + 1;
		unsigned k, kk, l;


		// Definition of ttt:
		// interaction ia_new would the ttt-th interaction in pia 
		// starting with 0
		//unsigned ttt= pia->get_pos(ia_new);
	
		
		// Definition of tt:
		// interaction ia_new would the tt-th interaction in ia_current 
		// starting with 0
		unsigned tt;
		tt=get_pos_cur(ia_new);
		
		// Definition of t:
		// interaction ia_new is the t-th regression coefficient 
		// without counting the intercept as coefficient
		unsigned t;
		t = ncoef-ncoef_ia+tt;

	

	   // ********** compute x_new *******************************************************
		double * workx_new;
		double * workx;
		double * workia;

		workx_new = x_new.getV();
		workx = getV_x();
		workia = pia->get_ia(ia_new);

		for(k=0; k<nobs; k++)
		{
			for(l=0; l<ncoef_new; l++, workx_new++, workx++)
			{
				if(l != t)
				{
					*workx_new = *workx;
				}
				else 
				{
					*workx_new = *workia;
					workia++;	
					workx--;
				}
			}
		}
		

		// *********** compute xx_new ********************************************

		double * workxx_new;
		double * workxx;
		double * workx1;
		double * workx2; 

		double value;

		workxx_new = xx_new.getV();
		workxx = getV_xx();

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


		// *********** compute proposed beta ********************************************

		double * workb_new;			// b_new: new regression coefficient
		double * workbeta;			// beta: old regression coefficient

		workb_new = b_new.getV(); //b_new.getV();
		workbeta = getV_beta_help();

		for(k=0; k<ncoef_new;	workb_new++, workbeta++, k++)
		{
			if(k!= t)
				*workb_new = *workbeta;
			else 
			{
				beta_new = rand_normal();		// the new coefficient (=u)
				*workb_new = beta_new;
				workbeta--;
			}
		}


		//*********** compute proposed linear predictor **********************************
		calc_lin_prop( x_new, b_new); 
		
	}
















	// FUNCTION: make_new_d_int
   // TASK: computes the new values for a death-step
   void FULLCOND_dag_ia::make_new_d_int (ST::string step, unsigned i, unsigned j, 
						unsigned num_ia_del, datamatrix & beta_old, 
						vector <vector <unsigned > > & current_ia_n, 
						datamatrix & xx_new, datamatrix & b_new, datamatrix & x_new)
	{
	    unsigned  k,l;
		unsigned  ncoef_new = ncoef-num_ia_del-1;

		if(proposal_version==2 && ia_d_there==false && step!="s")
		{
			x_ia_d = datamatrix(nobs, num_ia_del+1);
			xx_ia_d = datamatrix(num_ia_del+1, num_ia_del+1);

			if(mixed_case==false)
				ia_d_there=true;
		}

		double * workx_ia;
		double * workxx_ia;

		workx_ia = x_ia_d.getV();
		workxx_ia = xx_ia_d.getV();


		// Definition of ia
		// contains interactions of j-th regression model that contain i
		// which are those which are going to be deleted
		vector< vector< unsigned> > ia;
		ia_of_i(i,ia);


		

		// Definition of pos_of_del
		// contains positions of the main effect and the corresponding 
		// interactions in j-th regression model that are going to be deleted
		// (position of intercept is 0)
		vector <unsigned> pos_of_del;	
		get_pos(i,pos_of_del);

		//compute x_new
		double * workx_new;
		double * workx;

		workx_new = x_new.getV();
		workx = x.getV();


        unsigned t, pos;

		for(k=0; k<nobs; k++)
		{
            t=0;
            pos =  pos_of_del[t];

			for(l=0; l<ncoef; l++, workx_new++, workx++)
			{
				if ( l!=pos )
					*workx_new = *workx;
				else 
				{
                     if(proposal_version==2 && step!="s")
                     {
					    *workx_ia= *workx;
						workx_ia++;
                     }
					workx_new--;
                    t++;
                    if(t<num_ia_del+1)
                         pos =  pos_of_del[t];
				}
			}
		}





		//compute xx_new
		double * workxx_new;
		double * workxx;

		workxx_new = xx_new.getV();
		workxx = xx.getV();

		unsigned t1, pos1;
        unsigned t2, pos2;
        unsigned t3, pos3;

        t1=0;
        pos1 = pos_of_del[t1];

		for(k=0; k<ncoef; k++)
		{
            t2=0;
            pos2 = pos_of_del[t2];

			if( k!=pos1)
			{
				for(l=0; l<ncoef; l++, workxx++, workxx_new++)
				{
					if(l!=pos2)
						*workxx_new = *workxx;
					else
					{
						workxx_new--;

                        t2++;
                        if(t2<num_ia_del+1)
                            pos2 = pos_of_del[t2];
					}
				}
			}
			else
			{
                t3=0;
                pos3 = pos_of_del[t3];

				for(l=0; l<ncoef; l++, workxx++)
				{
					if(l==pos3)
					{
						if(proposal_version==2 && step!="s")
						{
							*workxx_ia = *workxx;
							workxx_ia++;
						}

                        t3++;
                        if(t3<num_ia_del+1)
                           pos3=pos_of_del[t3];
					}
                }

                t1++;
                if(t1<num_ia_del+1)
                   pos1 = pos_of_del[t1];
			}
		}






		if(step != "s") //if the death-step is NOT part of a switch step
		{
			make_prop_beta ('d', b_new, beta_old, x_new, xx_new, ncoef_new,
							pos_of_del);

			calc_lin_prop(x_new, b_new);
		}
	}













	// FUNCTION: make_new_b_int
	// TASK: computes the new values for a birth-step when one new main effect
	// and ALL corresponding interactions are added
	// is called from birth_step in rj_int

	void FULLCOND_dag_ia::make_new_b_int (ST::string step, unsigned i, vector <vector <unsigned> > ia_new,
											datamatrix & beta_new, datamatrix & xx_new,
											datamatrix & b_new, datamatrix & x_new)
	{
		unsigned num_ia_new = ia_new.size();
		unsigned ncoef_new = ncoef + num_ia_new +1;
		unsigned k, kk, l, m;

		// Definition of vec_t:
		// the t-th interaction of ia_new is the vec_t[t]-th regression coefficient
		// without counting the intercept as coefficient; t>0
		// t=0 corresponds to the main effect
		vector<unsigned> vec_t;

		// Definition of vec_t[0]:
		// the added main effect i is the vec_t[0]-th regression coefficient of j
		// without counting the intercept as coefficient
		unsigned t;
		t=1;
		for(k=0; k<nvar, k<i; k++)
		{
			if(adcol(k,0)==1)
				t++;
		}
		vec_t.push_back(t);


		for(k=0; k<num_ia_new; k++)
		{
			// Definition of tt:
			// interaction ia_new would the tt-th interaction in ia_current
			// starting with 0
			unsigned tt;
			tt=get_pos_cur(ia_new[k])+k;


			// Definition of t:
			// interaction ia_new is the t-th regression coefficient
			// without counting the intercept as coefficient
			// +1 stands for the added main effect
			t = ncoef-ncoef_ia+tt+1;
			vec_t.push_back(t);
		}


		
		


	   // ********** compute x_new *******************************************************

		double * workx_new;
		double * workx;


		// vec_workia contains pointers to the main effect resp. the interactions which are added 
		double * workia;
		vector <double *> vec_workia;

		//first element is pointer to the first observation of the added main effect i
		workia = data.getV()+i;
		vec_workia.push_back(workia);



		for(k=0; k<num_ia_new; k++)
		{
			workia = pia->get_ia(ia_new[k]);
			vec_workia.push_back(workia);
		}


// *********  x_ia and xx_ia are only needed when proposal_version==2  **********
		if(proposal_version==2 && ia_b_there==false && step!="s")
		{
			x_ia_b = datamatrix(nobs, num_ia_new+1);
			xx_ia_b = datamatrix(num_ia_new+1, num_ia_new+1);

			if(mixed_case==false)
				ia_b_there=true;
		}

		double * workx_ia;
		double * workxx_ia;

		workx_ia = x_ia_b.getV();
		workxx_ia = xx_ia_b.getV();
// ***************************************************************************



		workx_new = x_new.getV();
		workx = x.getV();

        unsigned position=0;


		for(k=0; k<nobs; k++)
		{	
			m = 0;
			for(l=0; l<ncoef_new; l++, workx_new++, workx++)
			{
				if(m<vec_t.size())
					position = vec_t[m];

				if(l != position)          //column does not belong to an added coefficient
				{
                    *workx_new = *workx;    //x_new(k,l)= x(k,l-m);
                    if(mixed_case ==false)
						assert(*workx_new==0 ||*workx_new==1);
				}
				else
				{
                   *workx_new = *(vec_workia[m]);   //	x_new(k,l) = *(vec_workia[m]) ;

					if(proposal_version==2 && step!="s")
					{
						if(mixed_case ==false)
							assert(*workx_new==0 ||*workx_new==1);
						*workx_ia = *workx_new;
						workx_ia ++;
					}

					if(m==0)	
						vec_workia[m]= vec_workia[m]+nvar;
					else
						vec_workia[m]++;

					workx--;
					m++;
				}
			}
		}


		// *********** compute xx_new ********************************************

		double * workxx_new;
		double * workxx;
		double * workx1;
		double * workx2;

		double value;
		unsigned m1, m2, m3, t1, t2, t3;

		m1 = 0;
		m2 = 0;

		t1 = vec_t[m1];
		t2 = vec_t[m2];

		workxx_new = xx_new.getV();
		workxx = getV_xx();

		for(k=0; k<ncoef_new; k++)
		{
			if(k != t1)
			{
				m2 = 0;
                t2 =vec_t[m2];

				for(l=0; l<ncoef_new; l++, workxx++, workxx_new++)
				{
					if(l != t2)
						*workxx_new = *workxx;
					else
					{
						workx1 = x_new.getV() + k;
						workx2 = x_new.getV() + l;
						value =0;

						for(kk=0; kk<nobs; kk++)
						{
							value = value + (*workx1) * (*workx2);
							workx1 = workx1 + ncoef_new;
							workx2 = workx2 + ncoef_new;
						}

						*workxx_new = value;

						workxx--;

                        m2++;
						if(m2<vec_t.size())
						    t2 = vec_t[m2];
					}
				}
			}
			else
			{
				m3 = 0;
				t3 = vec_t[m3];

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

					if(l == t3)
					{
						if(proposal_version==2 && step!="s")
						{
                            *workxx_ia = value;
                            workxx_ia++;
						}

                        m3++;
                        if(m3<vec_t.size())
                            t3 = vec_t[m3];
					}
				}

                m1++;
                if(m1<vec_t.size())
				    t1 = vec_t[m1];
			}
		}

 

		if(step !="s")
		{
			// *********** compute proposed beta ********************************************
			//make_prop_beta ('b', b_new, beta_new);
			make_prop_beta ('b', b_new, beta_new, x_new, xx_new, ncoef_new, vec_t);

			//*********** compute proposed linear predictor **********************************
			calc_lin_prop( x_new, b_new);
		}
	}












	// FUNCTION: get_pos_cur
	// TASK: gives back the position of ia in the vector current_ia
	//		 starting with 0; independent if it is already there or not
	unsigned FULLCOND_dag_ia::get_pos_cur(vector <unsigned> ia)
	{
		unsigned t;
		t=0;

		if(ncoef_ia>0)
		{
			assert(current_ia.size()==ncoef_ia);

			if(ia > current_ia[ncoef_ia-1])
				t=ncoef_ia; 
			else
			{
				std::vector< vector< unsigned> > :: iterator it_k;
				it_k = current_ia.begin();

				while( ia > *it_k )
				{
					t++;
					++it_k;
				}
			}
		}	
		return t;
	}




	// FUNCTION: make_prop_beta
	// TASK: computes the proposed regression coefficient beta_new
	// stores its new components in u
	void FULLCOND_dag_ia::make_prop_beta (char step,  datamatrix & beta_new,
			datamatrix & u, const datamatrix & x_new, const datamatrix & xx_new, 
			unsigned ncoef_new, const vector <unsigned> vec_t)
	{
		unsigned k, t, m;
		double sum;


		if(step=='d')
		{
			datamatrix mu (ncoef-ncoef_new,1);

			//compute proposed beta
			double * workb_new;
			double * workbeta;
			double * worku;

			workb_new = beta_new.getV();
			workbeta = beta_help.getV();
			worku = u.getV();

			m = 0;
			t = vec_t[m];

			sum=0;

			//calculate the conditioned mu

			for(k=0; k<ncoef; workb_new++, workbeta++, k++)
			{
				if(k!= t)
					*workb_new = *workbeta;
				else 
				{
					*worku= *workbeta;
					worku++;
					workb_new--;
					m++;

                    if(m<vec_t.size())
                       t = vec_t[m];
				}
			}

			calc_lin_prop(x_new, beta_new);

			y_ia.minus(y,lin_prop);

			calc_kq_est (mu, x_ia_d, xx_ia_d, y_ia);

			ln_prop_beta = p_prop(u, mu,1);

            //ln_prop_beta = sum;
		}
		else if(step == 'b')
		{
			if(proposal_version==0) //just proposal with N(0,sigma)
			{
			}
			else if (proposal_version==1) // all components are new 
										  // and correspond to the ls-estimator
			{
				calc_kq_est (beta_new, x_new, xx_new);

				double mu,  std, help;
				double sigma = 0.1;
				std = sqrt(sigma);

				double * workbeta_new;
				workbeta_new = beta_new.getV();
                
				sum=0;
				for(k=0; k<ncoef_new; k++, workbeta_new++)
				{
                    mu = *workbeta_new;
					help = mu  + std*rand_normal();
					*workbeta_new = help;
					sum = sum + p_prop(help, mu,sigma);
				}

				ln_prop_beta = sum; 

			}
			else if (proposal_version==2)	// only the added components are sampled
											// and correspond to ls-estimator giving the others
			{
				datamatrix y_ia (nobs,1);
				y_ia.minus(y,lin);
				calc_kq_est (u, x_ia_b, xx_ia_b, y_ia);
			}

			if(proposal_version==0 || proposal_version==2 )
			{
				double * pbeta;
				double * pbeta_new;
				double * pu;                  //mean of proposal

				pbeta = beta_help.getV();
				pbeta_new = beta_new.getV();
				pu = u.getV();


                m=0;
				t=vec_t[m];

                sum=0;
				for(k=0; k<ncoef_new; k++, pbeta_new++)
				{
					if(t !=k)
					{
						*pbeta_new = *pbeta ;
						pbeta++;
					}
					else
					{
						*pbeta_new = *pu + rand_normal();     // variance=1
						 sum = sum + p_prop(*pbeta_new, *pu, 1);
                         *pu = *pbeta_new;

						pu++;
						m++;
                        if(m<vec_t.size())
						    t=vec_t[m];
					}
				}
				ln_prop_beta = sum; 
			}
		}
	}



	// FUNCTION: get_pos
	// TASK: gives back the position of main effect i and the corresponding 
	// ia in the regression model starting with 0 for the intercept
	void FULLCOND_dag_ia::get_pos(unsigned i, vector <unsigned> & pos)
	{
		unsigned k;
		unsigned position;

		assert(pos.size()==0);

		unsigned t=1;
		for(k=0; k<nvar, k<i; k++)
		{
			if(adcol(k,0)==1)
				t++;	
		}
		pos.push_back(t);

		std::vector <vector <unsigned > > ::iterator it_cur;
		it_cur = current_ia.begin();

		for(k=0; k<ncoef_ia; k++, ++it_cur)
		{
			if(	(*it_cur)[0]==i  ||	(*it_cur)[1]==i)
			{
				position = ncoef-ncoef_ia+k;
				pos.push_back(position);
			}
		}	
	}








	void FULLCOND_dag_ia::write_to_beta_ia(void)
	{
		if((optionsp->get_nriter() > optionsp->get_burnin()) &&
		   (optionsp->get_nriter() % (optionsp->get_step()) == 0))
		{

			/*
			unsigned i,k;

			double * pbeta_help;
			double * pbeta;
			std::vector< unsigned> :: iterator  it_occ;
			
			k = ncoef-ncoef_ia;  // beginning of the interactions in the vector beta_help

			it_occ = occurrence.begin();
			pbeta = beta.getV() + k;
			pbeta_help = beta_help.getV();

			for(i=nvar; i<max_ncoef; pbeta++, ++it_occ, i++)
			{
				if( *it_occ == 1 )
				{
					*pbeta = *pbeta_help;
					pbeta_help++;
				}
				else
					*pbeta  = 0;	
			}	*/


			
			unsigned i,k;

			double * pbeta_help;
			double * pbeta;
			std::vector< int> :: iterator  it_occ;
			
			k = ncoef-ncoef_ia;  // beginning of the interactions in the vector beta_help

			it_occ = occurrence.begin();
			pbeta = beta.getV() + nvar;
			pbeta_help = beta_help.getV() + k;

			for(i=nvar; i<nvar+nvar*(nvar-1)/2;  pbeta++, ++it_occ, i++)   
			{
				if( *it_occ == 1 )
				{
					*pbeta = *pbeta_help;
					pbeta_help++;
				}
				else
					*pbeta  = 0;	
			}	





		}
	}





	// FUNCTION: change_occur()
	// TARGET: adds or delets ia in occur after birth or death step
	void FULLCOND_dag_ia::change_occur(char step, vector<unsigned> ia_new)
	{
        assert(step=='b' || step=='d');
		unsigned pos = pia->get_pos (ia_new);

		if(step == 'b')
			occurrence[pos] = 1;
		else
			occurrence[pos] = 0; 
	}




	// FUNCTION: change_occur()
	// TARGET: adds or delets ia in occur after birth or death step
	void FULLCOND_dag_ia::change_occur(char step, vector <vector <unsigned > > ia_vec)
	{
        assert(step=='b' || step=='d');

		unsigned t;
		unsigned pos; 
		vector<unsigned> ia_new;

		// strage but true: iterator decreases speed!!!

		//std::vector<unsigned> :: iterator it_occ ; 
		//it_occ = occurrence.begin();

		//k=0;
		for(t=0; t<ia_vec.size(); t++)
		{
			ia_new=ia_vec[t];
			pos = pia->get_pos (ia_new);

			if(step == 'b')
				occurrence[pos]=1;
				//*it_occ=1;	
			else
				occurrence[pos]=0;
				//*it_occ=0;
		}
	}





	// FUNCTION: change_current()
	// TARGET: adds or delets ia.term to/from current_ia
	void FULLCOND_dag_ia::change_current(char step, vector<unsigned> term)
	{
        assert(step=='d' || step=='b');

		if(step=='b')
		{
			if(ncoef_ia==0 || current_ia[ncoef_ia-1]<term)
				current_ia.push_back(term);
			else
			{
				std::vector < vector <unsigned> > :: iterator pos; 
				pos = current_ia.begin();

				while (*pos <term)
					pos++;

				current_ia.insert(pos, term);
			}
			ncoef_ia++;					
		}
		else
		{
			std::vector < vector <unsigned> > :: iterator pos; 
			pos = current_ia.begin() + get_pos_cur(term);
			current_ia.erase(pos);
			ncoef_ia--;
		}

	}






	// FUNCTION: change_current()
	// TARGET: adds or delets ia.term to/from current_ia
	void FULLCOND_dag_ia::change_current(char step, vector <vector <unsigned > > term)
	{
         assert(step=='b' ||step=='d');
         
		unsigned t, size;
		vector<unsigned> ia;

		size = term.size();
		if(step=='b')
		{
			for(t=0; t<size; t++)
			{
				ia = term[t];
				if(ncoef_ia==0 || current_ia[ncoef_ia-1]<ia)
				{
					current_ia.push_back(ia);
					ncoef_ia++;
				}
				else
				{
					std::vector < vector <unsigned> > :: iterator pos;
					pos = current_ia.begin();
					
					while(*pos <ia)
						++pos;

					current_ia.insert(pos, ia);
					ncoef_ia++;
				}
			}
		}
		else
		{
			std::vector < vector <unsigned> > :: iterator pos;
			for(t=0; t<size; t++)
			{	
				ia = term[t];
				
				pos = current_ia.begin() + get_pos_cur(ia);

				current_ia.erase(pos);
				ncoef_ia--;
			}	
		}

		ia_d_there=false;
		ia_b_there=false;
	}







	void FULLCOND_dag_ia::death_step (vector<unsigned> old_ia)
	{
		unsigned int ncoef_new = ncoef - 1;

		// instead of: datamatrix b_new (ncoef_new,1);
		datamatrix & b_new = get_b_new_d(); 	
		// instead of: datamatrix x_new (nobs,ncoef_new);
		datamatrix & x_new = get_x_new_d();
		// instead of: datamatrix xx_new (ncoef_new,ncoef_new);
		datamatrix & xx_new = get_xx_new_d();


		double beta_old; //coefficient which will vanish
			
		// computing of the new values
		make_new_d(old_ia, xx_new, beta_old, b_new, x_new);
		


		// calculate ratio

		double ratio;
		double log_num;
		double log_denom;

		log_num = calc_SQT_x(x_new, b_new) + calc_SQT_b(b_new);
		log_denom = calc_SQT_x() + calc_SQT_b();

		ratio = -1/(2*sigma_i) 
				*(log_num - log_denom) + p_prop(beta_old) ;




		// accept proposal and change the corresponding values

		if(func::accept(ratio) == true)
		{
			change_current('d', old_ia);
			change_occur('d', old_ia);
			change(99, b_new, x_new, xx_new, ncoef_new);

		}

	}





	 // FUNCTION: make_new_d
   // TASK: computes the new values for a death-step

   void FULLCOND_dag_ia::make_new_d ( vector<unsigned> ia_old, datamatrix & xx_new, 
							double & beta_old, datamatrix & b_new, datamatrix & x_new)
   {	
	    unsigned int k,l;

		

		// Definition of t:
		// the interaction ai_old i is the t-th regression coefficient 
		// of j (without intercept)

		unsigned t;
		t = ncoef-ncoef_ia+get_pos_cur(ia_old);


		//compute x_new
		double * workx_new;
		double * workx;

		workx_new = x_new.getV();
		workx = getV_x();

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

		//compute xx_new
		double * workxx_new;
		double * workxx;

		workxx_new = xx_new.getV();
		workxx = getV_xx();

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
		beta_old = beta_help(t,0);	

		//compute proposed beta
		double * workb_new;
		double * workbeta;

		workb_new = b_new.getV(); 
		workbeta = getV_beta_help();

		for(k=0; k<ncoef; workb_new++, workbeta++, k++)
		{
			if(k!= t)
				*workb_new = *workbeta;
			else 
				workb_new--;

		}

		//compute proposed linear predictor
		calc_lin_prop(x_new, b_new);

	}






 // FUNCTION: update_occ
 // TASK: updates occurrence and the auxiliary variables like oc_old, etc.
 
 void FULLCOND_dag_ia::update_occ(void)
 {

	  if((optionsp->get_nriter() > optionsp->get_burnin()) &&
		 (optionsp->get_nriter() % (optionsp->get_step()) == 0))
	  {
		  register unsigned i;

		  std::vector<int>::iterator workocc;
		  workocc = occurrence.begin();
		  double * workoccmean = occmean.getV();

		  unsigned samplesize = optionsp->get_samplesize();

		  for(i=0;i<(nvar)*(nvar-1)/2;i++,workocc++,workoccmean++)
		  {
			  // updating betamean
              if (samplesize==1)
			    *workoccmean= double(*workocc);
			  else
			    *workoccmean = (1.0/(samplesize))*
					 ((samplesize-1)*(*workoccmean) + double(*workocc));
		  }  // end: for i=0; ...
	  } // end: if ( (nriter > optionsp->burnin) && (nriter % optionsp->step == 0) )
  }












void FULLCOND_dag_ia::outresults(void)
    {		
	
	FULLCOND_dag::outresults(); // this command should be included at the beginning
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

	 optionsp->out(" **********  INTERACTIONS OF REGRESSION MODEL "
		+ ST::inttostring(self)
		+" *************\n"),
		optionsp->out("\n");
		optionsp->out("\n");

		unsigned r;
		bool any_ia_there = false;

	 if(detail==false)
	 {

		 for(r=0; r < nvar*(nvar-1)/2;r++)
		 {
			 unsigned first = (all_ia[r])[0];
			 unsigned second = (all_ia[r])[1];

		     if(first!=self && second !=self)
			 {
				 if(occmean(r,0)>0)
				 {
					 optionsp->out("rel. frequency of ia " 
									+  ST::inttostring(first) 
									+  ST::inttostring(second) 
									+ ": " 
									+  ST::doubletostring(occmean(r,0),5) + "\n");
					 optionsp->out("\n");

					 any_ia_there = true;
				 }
			}
		 }
	 }
	 else 
	 {

		 for(r=nvar; r<nvar+nvar*(nvar-1)/2; r++)				
		{
			unsigned first = (all_ia[r-nvar])[0];    // all_ia.size = nvar(nvar-1)/2
			unsigned second = (all_ia[r-nvar])[1];	

			if(first!=self && second !=self)
			{

				if(occmean(r-nvar,0)>0)    //occmean.size()=nvar*(nvar-1)/2
				{
					
					optionsp->out("Interaction  " 
									+  ST::inttostring(first) 
									+  ST::inttostring(second) 
									+ " : \n");
					optionsp->out("\n");
 

					optionsp->out("mean: "	   +  ST::doubletostring(betamean(r,0),5) + "\n");
					if(flags[0] == 0)
					{

                    ST::string l1 = ST::doubletostring(lower1,4);
                    ST::string l2 = ST::doubletostring(lower2,4);
                    ST::string u1 = ST::doubletostring(upper1,4);
                    ST::string u2 = ST::doubletostring(upper2,4);

                    optionsp->out(l1 + "% quantile: " 
									 + ST::doubletostring(betaqu_l1_lower(r,0),5) + "\n");

                    optionsp->out(l2 + "% quantile: " 
									 + ST::doubletostring(betaqu_l2_lower(r,0),5) + "\n");

                    optionsp->out("50% quantile: " +  ST::doubletostring(betaqu50(r,0),5) + "\n");

                    optionsp->out(u1 + "% quantile: " 
									 + ST::doubletostring(betaqu_l2_upper(r,0),5) + "\n");

                    optionsp->out(u2 + "% quantile: " 
									 + ST::doubletostring(betaqu_l1_upper(r,0),5) + "\n");


					}

					optionsp->out("\n");
					optionsp->out("\n");

					any_ia_there = true;
				} //if
			} // if...
		} //for
	 } // else .... detail=true

	 if(any_ia_there == false)
		 optionsp->out("No interactions observed.");

	optionsp->out("\n");
	optionsp->out("\n");
}



	void FULLCOND_dag_ia::outoptions(void)
	{
		FULLCOND_dag::outoptions();

		if(self==0)
		{
			if(detail==true)
			{
			  if(print_dags==false)
			  {
				  optionsp->out("NOTE: The option 'detail_ia' presupposes the option 'print_dags'.\n");
			  }
			  else
				   optionsp->out("Estimations of interactions are listed in detail. \n");

			}
		}

		if(detail==true)
		{
			if(print_dags==false)
				detail=false;
		}
    }






	// FUNCTION: num_ia_new
   // TASK: returns the number of allowed new interactions of the new main effect i
   unsigned FULLCOND_dag_ia::num_ia_new(unsigned i)
   {
	   return ncoef_m-1;
   }




	} // namespace MCMC




