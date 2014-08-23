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


#if defined(BORLAND_OUTPUT_WINDOW)
using namespace std;
#endif

#if !defined (adja_INCLUDED)

#define adja_INCLUDED

#include "tmatrix.h"
#include <vector>
#include "clstring.h"

#include "statmat.h"
#include <list>

using std::vector;
using std::list;
using std::ifstream;

struct essfreq
{
	Matrix <unsigned >  sceleton;
	vector< vector <unsigned> > immoral;
	unsigned nedges;
	double freq;

	essfreq(void);

	essfreq(Matrix <unsigned> scel, vector< vector <unsigned> > & imm, unsigned num_edges, double frequency ) ;

	essfreq(unsigned dim);

	essfreq(const essfreq & ess) ;

    const essfreq & operator=(const essfreq & ess)
	{
		if (this == &ess)
			return *this;

		sceleton = ess.sceleton;
		immoral = ess.immoral;
		nedges = ess.nedges;
		freq = ess.freq;

		return *this;
	}


	// "<" (smaller) in the sense of smaller frequency
	int operator<(const essfreq & m1) const
	{
		if (freq < m1.freq)
			return 1;
		return 0;
	}


	friend int operator>(const essfreq & m1, const essfreq & m2)
	{
		return m2 < m1;
	}


	// "<<" (smaller) in the sense of number of edges
	int operator<<(const essfreq & m1) const
	{
		if (nedges < m1.nedges )
			return 1;
		return 0;
	}


	friend int operator>>(const essfreq & m1, const essfreq & m2)
	{
		return m2 << m1;
	}


};   //struct interact





//------------------------------------------------------------------------------
//---------------------------- CLASS: adjacency -------------------------------
//------------------------------------------------------------------------------




class adja : public Matrix<unsigned>

{

	public:

		//datamatrix adja;					// adjacency matrix
		vector <list <unsigned> > ladja;	// representastion of adja by a list


		unsigned nvar;		 // number of variables (=ncols=nrows)
		unsigned nedge;		 // number of edges
		unsigned start_type; // start_type=0: independent
							 // start_type=1: dependent1
							 // start_type=2: dependent2
							 // start_type=3: chain




		// DEFAULT CONSTRUCTOR:
		adja(void)
		{
			nvar=1;
		}

		// CONSTRUCTOR_1
		adja(unsigned n);

		// CONSTRUCTOR_2
		adja( unsigned n, unsigned type);

		// COPY CONSTRUCTOR
		adja(const adja & a);

		 // OVERLOADED ASSIGNMENT OPERATOR
	//	const adja & operator=(const adja & a);

		 // OVERLOADED ASSIGNMENT OPERATOR
	//	const datamatrix  & operator= (const adja & a);


		//const SparseMatrix & operator=(const SparseMatrix & m);


		// DESTRUCTOR
		~adja() {}

		// FUNCTION: make_list
		// TASK: represents adja-matrix m as a list
		 vector < list <unsigned int> > make_list (const adja & m) const;

		 // FUNCTION: make_list
		// TASK: represents Matrix<int> m as a list
		 vector < list <unsigned int> > make_list (const Matrix<unsigned> & m) const;

		 // FUNCTION: make_list
		// TASK: represents calling matrix as a list
		 vector < list <unsigned int> > make_list (void) const;


		// FUNCTION: change_list
		// TASK: changes the existing list ladja
		// when there is a birth(0), death(1) or switch(2) step
		void change_list (unsigned int i, unsigned int j, unsigned int step);


		// FUNCTION: compare
		// TASK: "multiplies" the lists l1 and l2
		// returns true, if the "product is unequal to zero
		bool compare(list <unsigned int> & l1, list <unsigned int>  & l2) const;


		// FUNCTION: azy_test
		// TASK: returns true, if no cycles when adding (ii,jj)
		bool azy_test(unsigned int ii, unsigned int jj);


		// FUNCTION: azy_test
		// TASK: returns true, if no cycles in the calling matrix
		bool azy_test(void);


		// FUNCTION: equi_test
		// TASK: tests if i->j is covered,
		// e.g. if by changing i->j into j->i an equivalent graph is got
		bool equi_test(unsigned i, unsigned j);


		// FUNCTION: edge_plus
		// TASK: increases number of edges by one
		void edge_plus (void)
		{
			nedge ++;
		}

		// FUNCTION: edge_minus
		// TASK: increases number of edges by one
		void edge_minus (void)
		{
			nedge --;
		}

		// FUNCTION: get_nedge
		// TASK: returns number of edges
		unsigned get_nedge (void)
		{
			return nedge;
		}

		// FUNCTION: get_nvar
		// TASK: returns number of variables
		unsigned get_nvar (void)
		{
			return nvar;
		}



        // FUNCTION: num_pa
		// TASK: returns the number of parents of variable i
		unsigned int  num_pa(unsigned i) const
		{
			unsigned int num =0;
			for(unsigned int k=0; k<nvar; k++)
			{
			  num = num + get(k,i);
			}
			return num;
		 }



		// FUNCTION: is_pa
		// TASK: returns true if i is parent of j
		bool  is_pa(unsigned i, unsigned j) const
		{
			if(get(i,j)==1)
				return true;
			else
				return false;
		 }



		// FUNCTION: string_to_adja
		// TASK: changes calling adjacency matix into adjacency
		// matrix that corresponds to list of modfreq
		void string_to_adja (ST::string model);



		// FUNCTION: string_to_adja
		// TASK: changes calling adjacency matix into adjacency
		// matrix that corresponds to list of modfreq
		// and sets num_edges equal to the number of edges
		void string_to_adja (ST::string model, unsigned & num_edges);



		// FUNCTION: adja_to_ess
		// TASK: changes calling adjacency matix into
		// sceleton and vector of immoralities
		void adja_to_ess ( Matrix <unsigned> & scel, vector< vector <unsigned> > & imm)  ;



		// FUNCTION: assign
		// TASK: assigns the elements of A to the elements of the calling matrix
		//       faster than B = A or B.putBlock(A,0,0,B.rows,B.cols())
		// void assign(const statmatrix & A);



		//	Datenzeiger zugreifbar
		unsigned *getV() const { return m_v; }




		// FUNCTION: adja_to_ess
		// TASK: changes calling adjacency matix into
		// sceleton and vector of immoralities
		void read_ess ( vector< essfreq> & list_ess, datamatrix & mean_all, datamatrix & sum_square,
						ifstream & fin, unsigned number);


		// FUNCTION: add_ess_to_list
		// TASK: adds ess_new to ther lists and updates frequencies
		void add_ess_to_list(vector <essfreq> & list_ess, essfreq & ess_new);



		// FUNCTION: add_to_mean
		// TASK: takes matrix_new into account for the matrix_mean
		void add_to_mean(const datamatrix & matrix_new, datamatrix & matrix_mean,  unsigned n);


		// FUNCTION: add_to_mean
		// TASK: takes matrix_new into account for the matrix_mean
		void add_to_square(const datamatrix & matrix_new, datamatrix & sum_square, unsigned n);


		// FUNCTION: write_out_ess_short
		// TASK: writes out the essential graph ess into a separate file
		void write_out_ess_short(essfreq & ess, ST::string path_res, unsigned n);













};


//} //namespace std


#endif
