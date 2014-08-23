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


#if !defined (IA_MIXED_INCLUDED)

#define IA_MIXED_INCLUDED

#include "clstring.h"
#include "statmat.h"
#include "ia.h"

namespace MCMC
{



class IA_MIXED  : public IA 
{

	protected:

		


	public:



		// DEFAULT CONSTRUCTOR:
		IA_MIXED(void);

		// CONSTRUCTOR_1
		IA_MIXED(const datamatrix & d);

		// CONSTRUCTOR_2
		// for interactions of order>2 (some day in future....) 
		IA_MIXED(unsigned order, const datamatrix & d);

		// COPY CONSTRUCTOR
		IA_MIXED (const IA_MIXED & a);

		 // OVERLOADED ASSIGNMENT OPERATOR
		const IA_MIXED & operator=(const IA_MIXED & a);


		// DESTRUCTOR
		~IA_MIXED() {}  


/*		// FUNCTION: make_list
		// TASK: represents adja-matrix m as a list 
		// vector < list <unsigned int> > make_list (const adja & m) const;


		// FUNCTION: make_ia_term
		// TASK: creates a new interaction term
		void make_ia (vector<unsigned> terms);

		
	  // FUNCTION: choose_ia_term
      // TASK: chooses a new interaction term of order ord 
      // which is NOT already in current_ia
		vector<unsigned> IA::choose_ia (const Matrix<unsigned> & col, 
									 vector <vector <unsigned> > & current_ia);

	  // FUNCTION: choose_ia_term
      // TASK: chooses a new interaction term of order ord 
      // regardless if it is already in current_ia or not
		vector<unsigned> IA::choose_ia ( const  Matrix<unsigned>  & col);




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


*/


	

  

};    // class

}  //namespace

#endif



