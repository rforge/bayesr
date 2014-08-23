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



#include"../export_type.h"

#if !defined (IA_INCLUDED)

#define IA_INCLUDED

#include"clstring.h"
#include"statmat.h"


namespace MCMC
{


struct interact
{

	vector<unsigned> ia_term;
	datamatrix ia_dat;


	interact(void)
	{
		vector<unsigned> help;
		ia_term = help ;
	}

	interact(vector<unsigned> & inter)
	{
		ia_term = inter;
	}

	interact(vector<unsigned> &  inter, datamatrix   ia)
	{
		ia_term = inter;
		ia_dat = ia;
	}

//	interact(const ST::string &  inter, datamatrix * pointer_ia)
//	{
//		interaction= inter ;
//		ia = ia;
//	}

	interact(const interact & ia)
	{
		ia_term = ia.ia_term ;
		ia_dat = ia.ia_dat;
	}


    const interact & operator=(const interact & ia)
	{
		if (this == &ia)
			return *this;

		ia_term= ia.ia_term ;
		ia_dat = ia.ia_dat;

		return *this;
	}



	int operator<(const interact & ia) const
	{

/*		// the following holds true for interactions with order >2

		unsigned i;

		if(interaction.size()<(ia.interaction).size())
			return 1;
		if(interaction.size()>(ia.interaction).size())
			return 0;
		else
		{
			for(i=0; i<interaction.size(); )
			{
				if(interaction[i]<ia.interaction[i])
					return 1;
				else if (interaction[i]==ia.interaction[i])
					i++;
				else
					return 0;
			}
			return 0;
		}
*/


		// the following holds true for interactions with order equal to 2
		assert(ia_term != ia.ia_term);

		if(ia_term[0]>(ia.ia_term)[0])
			return 0;
		else if(ia_term[0]<(ia.ia_term)[0])
			return 1;
		else
		{
			if(ia_term[1]>(ia.ia_term)[1])
				return 0;
			else if(ia_term[1]<(ia.ia_term)[1])
				return 1;
      else
        return 0;
		}
	}




	friend int operator>(const interact & ia1, const interact & ia2)
	{
		return ia2 < ia1;
	}

};   //struct interact




class __EXPORT_TYPE IA

{
	protected:

		unsigned nobs;
		unsigned nvar;

		unsigned max_ia_order;		// maximal order of interaction terms
		unsigned max_num;			// maximal number of interaction terms


		unsigned number_of_continuous;
		unsigned number_of_discrete;


		datamatrix data;

		vector <interact> ia_var;	// contains the vectors of ALL possible interactions
									// but the corresponding variable ONLY
									// if it has already occured

		vector<int> occurred;	// contains 1 if interactionvariable of this position
									// has already been constructed, otherwise 0

		vector <char> var_type;

		bool mixed_case;

		/*
		unsigned max_nia;			// maximal number of interaction
									// terms in one regression model

		double prob_ia;				// probability of one interaction term
		vector <interact>  ia;		// vector of structures interact
	   */


	public:

		// DEFAULT CONSTRUCTOR:
		IA(void);


		// CONSTRUCTOR_1
		IA(const datamatrix & d);



		// CONSTRUCTOR_2
		// for interactions of order>2 (some day in future....)
		IA(unsigned order, const datamatrix & d);

		// COPY CONSTRUCTOR
		IA (const IA & a);

		 // OVERLOADED ASSIGNMENT OPERATOR
		const IA & operator=(const IA & a);


		// DESTRUCTOR
		~IA() {}


		// FUNCTION: make_list
		// TASK: represents adja-matrix m as a list
		// vector < list <unsigned int> > make_list (const adja & m) const;


		// FUNCTION: make_ia_term
		// TASK: creates a new interaction term
		void make_ia (vector<unsigned> terms);


	  // FUNCTION: choose_ia_term
      // TASK: chooses a new interaction term of order ord
      // which is NOT already in current_ia
		vector<unsigned> choose_ia (const Matrix<unsigned> & col,
									 vector <vector <unsigned> > & current_ia);

	  // FUNCTION: choose_ia_term
      // TASK: chooses a new interaction term of order ord
      // regardless if it is already in current_ia or not
		vector<unsigned> choose_ia ( const  Matrix<unsigned>  & col);



	  // FUNCTION: choose_ia_term
      // TASK: chooses a new interaction term of an orderorder ord
      // which is NOT already in current_ia
	//	vector<unsigned> IA::choose_ia ( unsigned ord, const datamatrix & col,
	//								 const vector <vector <unsigned> > & current_ia);


	   // FUNCTION: already_there (vec_ia, current_ia)
       // TASK: returns true if the interaction vec_ia is already in the current model
		bool already_there ( const vector<unsigned> & vec_ia,
						vector <vector <unsigned> > & current_ia);


		// FUNCTION: already_there (vec_ia)
       // TASK: returns true if the interaction vec_ia is already in ia_var
		bool already_there ( const vector<unsigned> & vec_ia);




		// FUNCTION: string vec_to_str (terms)
		// TASK: changes terms (vector of numbers) into a string
		ST::string vec_to_str (vector<unsigned> terms);



		// FUNCTION: add_ia
		// TASK: adds datamatix ia.ia_dat to the corresponding entry of ia_var
		void add_ia(interact ia) ;




		// FUNCTION: add_ia
		// TASK: adds datamatix ia.ia_dat to interaction at ia_var[pos]
		void add_ia(datamatrix & data, unsigned pos);




		// FUNCTION: get_interaction
		// TASK: returns the vector of the k-th interaction of ia_var
		vector<unsigned> get_interaction(unsigned k)
		{
			return ia_var[k].ia_term;
		}



		// FUNCTION: get_ia
		// TASK: returns pointer to the first element of the matrix of the k--th interaction
		double * get_ia(unsigned k)
		{
			return (ia_var[k].ia_dat).getV();
		}


		// FUNCTION: get_ia
		// TASK: returns pointer to the first element of the matrix of interaction ia
		// regardless if it has already existed before or not
		double * get_ia(vector<unsigned> ia);






		// FUNCTION: get_max_order
		// TASK: returns the maximal allowed order of the interaction terms
		unsigned get_max_order(void)
		{
			return max_ia_order;;
		}



		// FUNCTION: get_pos
		// TASK: gives position of ia if all possible interactions of order 2
		// of nvar variables are stored in an ordered vector
		unsigned get_pos(vector<unsigned> ia);


		// FUNCTION: ia_okay
		// TASK: returns true if interaction is allowed
		bool ia_okay (vector<unsigned> & terms);

		// FUNCTION: ia_okay
		// TASK: returns true if interaction is allowed
		bool ia_okay (unsigned i);


		// FUNCTION: tell_var_type
		// TASK: gives back the var_kind of the i--th variable
		char tell_var_type (unsigned i);


		// FUNCTION: mixed_case
		// TASK: says if mixed_case or not
		bool tell_mixed_case (void)
		{
			return mixed_case;
		}


		// FUNCTION: give_var_kind
		// TASK: gives back the right values for num_cont and num_disc
		void give_var_kind (const Matrix<unsigned> & adc,
						unsigned & num_cont, unsigned & num_disc);



		// FUNCTION: update
		// TASK: updates regression model by reversible jump mcmc
		//void update(void);








};    // class

}  //namespace

#endif



