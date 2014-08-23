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





#include "adjacency.h"

#include <iostream>
#include "clstring.h"
#include <algorithm>


#include <stdlib.h>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;

essfreq::essfreq(void)
	{
		Matrix <unsigned> help1 ;
		sceleton = help1;

		vector< vector <unsigned> > help2;
		immoral = help2;

		nedges = 0;
		freq=0;
	}


essfreq::essfreq(Matrix <unsigned> scel, vector< vector <unsigned> >  & imm, unsigned num_edges, double frequency )
	{

		sceleton = scel;

        if(imm.size() ==0)
        {
            vector< vector <unsigned> > help;
		   immoral = help;
        }
        else
		   immoral = imm;

		nedges =  num_edges;
		freq = frequency;
	}

 essfreq::	essfreq(unsigned dim)
	{
		sceleton = Matrix <unsigned> (dim,dim, 0);
		vector< vector <unsigned> > help2;
		immoral = help2;
		nedges = 0;
		freq=0;
	}

  essfreq::	essfreq(const essfreq & ess)
	{
		sceleton = ess.sceleton;
		immoral = ess.immoral;
		nedges = ess.nedges;
		freq = ess.freq;
	}












//namespace std
//{


	// CONSTRUCTOR_1
	adja::adja(unsigned n): Matrix<unsigned > (n, n, 0)
	{
		nvar = n;
		start_type = 0; //independence
		nedge = 0;
		ladja = make_list();
	}


	// CONSTRUCTOR_2
	adja::adja( unsigned n, unsigned type): Matrix<unsigned> (n, n, 0)
	{
		unsigned i,j;
		start_type = type;
		nvar = n;
		nedge = 0;

		assert(type<5);

		if(start_type == 1)  // complete model; edges always to "higher" variables
		{
			for(i=0; i<nvar; i++)
			{
				for(j=0; j<nvar; j++)
				{
					if(i<j)
						put(i,j,1);
				}
			}
		}
		else if (start_type == 2)  // complete model; edges always to "lower" variables
		{
			for(i=0; i<nvar; i++)
			{
				for(j=0; j<nvar; j++)
				{
					if(i>j)
						put(i,j,1);
				}
			}
		}
		else if (start_type == 3) // "chain"; edges always to "higher" variables
		{
			for(i=0; i<nvar; i++)
			{
				for(j=0; j<nvar; j++)
				{
					if(i+1==j)
						put(i,j,1);
				}
			}
		}
		else if (start_type == 4) // "chain"; edges always to "lower" variables
		{
			for(i=0; i<nvar; i++)
			{
				for(j=0; j<nvar; j++)
				{
					if(i==j+1)
						put(i,j,1);
				}
			}
		}

		if(start_type !=0)
		{
			assert(nvar >0);

			for(i=0; i<nvar; i++)
			{
				for(j=0; j<nvar; j++)
				{
					if(get(i,j)==1)
						nedge++;
				}
			}
		}

		ladja = make_list();
	}



	// COPY CONSTRUCTOR
	adja::adja(const adja & a)
	{
		nvar = a.nvar;
		nedge = a.nedge;
		start_type = a.start_type;
		ladja = a.ladja;
	}


	// OVERLOADED ASSIGNMENT OPERATOR
/*	const adja & adja::operator=(const adja & a)
    {
		if (this==&a)
			return *this;

		unsigned i,j,k;

	k=a.rows();
	k=rows();

		for( i=0; i<a.rows(); i++)
			for(j=0; j<a.cols(); j++)
				put(i,j,a(i,j));*/

/*		for( i=0; i<a.rows(); i++)
			for(j=0; j<a.cols(); j++)
				cout<<get(i,j)<<endl;;

*/

/*		nvar = a.nvar;
		nedge = a.nedge;
		start_type = a.start_type;
		ladja = a.ladja;

		return *this;
    }

*/



/*	// OVERLOADED ASSIGNMENT OPERATOR
	const datamatrix  & adja::operator= (const adja& a)
	{
		datamatrix res(a.rows(),a.cols());

		double * work = res.getV();
		register unsigned i,j;
		for(i=0;i<res.rows();i++)
			for(j=0;j<res.cols();j++,work++)
				 *work = a(i,j);
		return res;
	}

*/


	vector < list <unsigned > > adja::make_list (const adja & m) const
	{
		unsigned  i,j;
		vector< list <unsigned > >  help (nvar);

		//i is going through nvar components (which are lists) of the vector ladja

		for(i=0; i<nvar; i++)
		{
			list <unsigned > l;
			for(j=0; j<nvar; j++)
			{
				if(m.get(i,j)==1)
				{
					l.push_back(j);
				}
			}
			help[i] = l;
		}

		return help;
	}



	vector < list <unsigned int> > adja::make_list (const Matrix<unsigned> & m) const
	{
		unsigned int i,j;
		vector< list <unsigned int> >  help (nvar);

		//i is going through nvar components (which are lists) of the vector ladja

		for(i=0; i<nvar; i++)
		{
			list <unsigned int> l;
			for(j=0; j<nvar; j++)
			{
				if(m.get(i,j)==1)
				{
					l.push_back(j);
				}
			}
			help[i] = l;
		}

		return help;
	}




	vector < list <unsigned int> > adja::make_list (void) const
	{
		unsigned int i,j;
		vector< list <unsigned int> >  help (nvar);

		//i is going through nvar components (which are lists) of the vector ladja

		for(i=0; i<nvar; i++)
		{
			list <unsigned int> l;
			for(j=0; j<nvar; j++)
			{
				if(get(i,j)==1)
				{
					l.push_back(j);
				}
			}
			help[i] = l;
		}

		return help;
	}







	// FUNCTION: change_list
		// TASK: changes the existing list ladja_new
		// when there is a birth(0), death(1) or switch(2) step

	void adja::change_list (unsigned int i, unsigned int j, unsigned int step)
	{
		unsigned int k;
		std::list < unsigned int>::iterator pos_j;
		unsigned int t=0;


		if(step==0) // change of list after a birth-step (i has become parent of j)
		{
			 // there are t variables "in front of" i which are parents of j
			// so if i is a parent of j, it should be the (t+1)-th place of the list ladja[j]
			for(k=0; k<i; k++)
			{
				if((k,j)==1)
					t++;
			}

			for(k=0, pos_j = ladja[j].begin(); k<t+1; k++, ++pos_j)
			{
					if(k==t)
					 ladja[j].insert(pos_j,i);
			}
		}

		else if (step==1) // change list after a death-step (i is not parent of j anymore)
		{
			ladja[j].remove(i);
		}

		else // change list after a switch-step (now i is parent of j and not viceversa)
		{
			ladja[i].remove(j);

			 // there are t variables "in front of" i which are parents of j
			 // so if i is a parent of j, it should be the (t+1)-th place of the list ladja[j]

			for(k=0; k<i; k++)
			{
				 if((k,j)==1)
					 t++;
			}

			 for(k=0, pos_j = ladja[j].begin(); k<t+1; k++, ++pos_j)
			 {
				 if(k==t)
					 ladja[j].insert(pos_j,i);
			 }
		 } // else
	 }





	// FUNCTION: compare
	// TASK: "multiplies" the lists l1 and l2
	// returns true, if the "product is unequal to zero
	bool adja::compare (list <unsigned int> & l1, list <unsigned int>  & l2) const
	{
		bool equal = false;

		std::list < unsigned int>::iterator pos1,pos2;

		pos1 = l1.begin();
		pos2 = l2.begin();

		if(l1.empty() || l2.empty())
			return equal;

		else
		{
			// idea: if the x-th element of the i-th row of the proposed adja
			// and the x-th element of the j-th column of l2 of the transposed new adja
			// are both unequal to zero
			// then the ij-th element of the product is different from zero, too.
			// here the i-th row and the j-th colum are represented by the lists l1 and l2

			while(pos1 != l1.end() && pos2 != l2.end())
			{

				if( *pos1 < *pos2)
					while( (*pos1 < *pos2) && (pos1 != l1.end()) )
					{
						++pos1;
					}

				else if( *pos1 > *pos2)
					while( (*pos1 > *pos2) && (pos2 != l2.end()) )
					{
						++pos2;
					}

				else if(*pos1 == *pos2)
				{
					equal = true;
					pos2 = l2.end();
				}
			}

			return equal;
		}
	}




	// FUNCTION: azy_test
	// TASK: returns true, if no cycles when adding (ii,jj)

		bool adja::azy_test(unsigned int ii, unsigned int jj)
	{
		bool azy = true;
		unsigned int i,j,l, limit;

		// change of ladja as if ii->jj would be added
		// last argument=0, because of birth-step
		change_list(ii,jj,0);

		// change of beta as if ii_>jj would be added
		 put(ii,jj,1);

		vector <list <unsigned > > ladja_pot, ladja_t;
		ladja_pot = make_list();

		ladja_t = make_list(transposed());

		std::vector <list <unsigned > > :: iterator pos_i; // like i-th row of adja_pot
		std::vector <list <unsigned > > :: iterator pos_j; // like j-th column of adja

		list <unsigned int> new_pot_i;

		#if defined(MICROSOFT_VISUAL)
		{
			limit = __min(nedge+1, nvar);
		}
		#else
		{
// Vorschlag:
//			limit = min(nedge+1, nvar);
			limit = std::min(nedge+1, nvar);
		}
		#endif

		for(l=0; l<limit; l++)
		{
			for(pos_i = ladja_pot.begin(),i=0; i<nvar-l; ++pos_i,i++)
			{
				for(pos_j = ladja_t.begin(),j=0; j<nvar; ++pos_j,j++)
				{
					if ( compare(*pos_i,*pos_j) == true)
					{
						new_pot_i.push_back(j);
						if(i==j)
						{
							// cout<<"Zyklus von "<< i+1<<" nach "<<i+1<<" in "<<l+2<<" Schritten !!!"<<endl;
							azy=false;
							i=nvar;
							j=nvar;
							l=nvar;
						}
					}
				}
				if(azy==true)
				{
					*pos_i = new_pot_i;
					new_pot_i.clear();
				}
			}
		}

		put(ii,jj,0);			// beta is set back to the old value
		change_list(ii,jj,1);   // the birth-step is set back

		return azy;
	}




	// FUNCTION: azy_test
	// TASK: returns true, if no cycles in the calling matrix
	bool adja::azy_test(void)
	{
		bool azy = true;
		unsigned int i,j,l, limit;

		vector <list <unsigned> > ladja_pot, ladja_t;
		ladja_pot = make_list();


	//	datamatrix help (nvar, nvar);
	//	help.assign(transposed());

		ladja_t = make_list(transposed());

		std::vector <list <unsigned int> > :: iterator pos_i; // like i-th row of adja_pot
		std::vector <list <unsigned int> > :: iterator pos_j; // like j-th column of adja

		list <unsigned int> new_pot_i;

		#if defined(MICROSOFT_VISUAL)
		{
			limit = __min(nedge+1, nvar);
		}
		#else
		{
// Vorschlag:
//			limit = min(nedge+1, nvar);
			limit = std::min(nedge+1, nvar);
		}
		#endif

		for(l=0; l<limit; l++)
		{
			for(pos_i = ladja_pot.begin(),i=0; i<nvar-l; ++pos_i,i++)
			{
				for(pos_j = ladja_t.begin(),j=0; j<nvar; ++pos_j,j++)
				{
					if ( compare(*pos_i,*pos_j) == true)
					{
						new_pot_i.push_back(j);
						if(i==j)
						{
							// cout<<"Zyklus von "<< i+1<<" nach "<<i+1<<" in "<<l+2<<" Schritten !!!"<<endl;
							azy=false;
							i=nvar;
							j=nvar;
							l=nvar;
						}
					}
				}
				if(azy==true)
				{
					*pos_i = new_pot_i;
					new_pot_i.clear();
				}
			}
		}

		return azy;
	}


	// FUNCTION: equi_test
	// TASK: tests if i->j is covered,
	// e.g. if by changing i->j into j->i an equivalent graph is got
	bool adja::equi_test(unsigned i, unsigned j)
	{
		unsigned k;
		bool equi=true;

		unsigned pa_j = num_pa(j); //number of parents of j


		if( num_pa(i)+1  == pa_j)  //
		{
			unsigned * work_i;
			unsigned * work_j;

			work_i = getV()+i;
			work_j = getV()+j;

			for(k=0; k<nvar; k++)
			{
				if(*work_i == *work_j || k == i )
				{
					work_i= work_i+nvar;
					work_j= work_j+nvar;
				}
				else
				{
					k=nvar;
					equi=false;
				}
			}
		}
		else
			equi=false;



		return equi;
	}



	// FUNCTION: string_to_adja
	// TASK: changes calling adjacency matix into adjacency
	// matrix that corresponds to list of modfreq
	void adja::string_to_adja (ST::string model)
	{
		unsigned k,l,pos,space,count;

		space=0;
		count=0;

		for(k=0; k<nvar; k++)
		{
			for(l=0; l<nvar; count++, l++)
			{
				pos = count+space;
				if(model[pos]=='1')
					put(k,l,1);
				else
					put(k,l,0);
			}
			space++;
		}
	 }



	// FUNCTION: string_to_adja
	// TASK: changes calling adjacency matix into adjacency
	// matrix that corresponds to list of modfreq
	// and sets num_edges equal to the number of edges
	void adja::string_to_adja (ST::string model, unsigned & num_edges)
	{
		unsigned k,l,count;

		count=0;
		num_edges=0;

		for(k=0; k<nvar; k++)
		{
			for(l=0; l<nvar; count++, l++)
			{
				if(model[count]=='1')
				{
					put(k,l,1);
					num_edges++;
				}
				else
				{
					put(k,l,0);
				}
			}
		}
	 }


	// FUNCTION: adja_to_ess
	// TASK: changes calling adjacency matix into adjacency
	// matrix that corresponds to list of modfreq
	void adja::adja_to_ess ( Matrix <unsigned> & scel, vector< vector <unsigned> > & immoralities)
	{
		unsigned k,l,m;

		assert(nvar==scel.cols());
		assert(immoralities.size()==0);


		for(k=0; k<nvar; k++)
		{
			for(l=0; l<nvar; l++)
			{
				if(get(k,l)==1)
				{
					if(k<l)
					{
						scel(k,l)=1;
					}
					else
					{
						scel(l,k)=1;
					}
				}
			}
		}

		for(l=0; l<nvar; l++)
		{
			if(num_pa(l)>1)
			{
				for(k=0; k<nvar; k++)
				{
					if(get(k,l)==1)
					{
						for(m=k+1; m<nvar; m++)
						{
							if(get(m,l)==1 && is_pa(m,k)==false &&is_pa(k,m)==false)
							{
								vector<unsigned> imm;

								imm.push_back(l);
								imm.push_back(k);
								imm.push_back(m);

								immoralities.push_back(imm);
							}
						}
					}
				}
			}
		}
		std::sort(immoralities.begin(), immoralities.end());
	}








	void adja::read_ess ( vector< essfreq> & list_ess, datamatrix & mean_all,
							datamatrix & sum_square, ifstream & fin, unsigned number)
	{



		if(fin.is_open()==false)
			cout<<"fin is not open!!!"<<endl; // fin.open();

		unsigned k,l;
		essfreq ess_new (nvar);

		ST::string s;
		//char c;

		adja help_adja (nvar);
		adja scel(nvar);

		bool already_there;
		already_there=false;

		int elem_imm;
		unsigned num_edge;
		unsigned num_imm;

		vector <unsigned> help_c, help_int ;

		for(l=0; l<10;l++)
		{
			vector< vector <unsigned> > immo_all;
			fin>>s;

			if(s.length()>1)
			{
				scel.string_to_adja(s, num_edge);

				fin>>num_imm;

				unsigned i=0;
				while(i<num_imm)
				{
					vector <unsigned> immo;

					for(k=0; k<3; k++)
					{
						fin>>elem_imm;
						immo.push_back(elem_imm);
					}

					immo_all.push_back(immo);

					i++;
				}

				double freq;
				fin>>freq;

				ess_new = essfreq(scel, immo_all, num_edge, freq);

				help_adja.add_ess_to_list(list_ess, ess_new);
			}
			else
			{
				l=10;
				already_there=true;
			}
		}





		double element;
		datamatrix adja_mean(nvar, nvar);
		for(k=0; k<nvar; k++)
		{
			for(l=0; l<nvar; l++)
			{
				if(already_there==false)
					fin>>element;
				else
				{
					element = 0;
					already_there=false;
				}

				adja_mean(k,l) = element;
			}
		}

		help_adja.add_to_mean(adja_mean, mean_all, number);
		help_adja.add_to_square(adja_mean, sum_square, number);
	}






	// FUNCTION: add_ess_to_list
	// TASK: adds ess_new to ther lists and updates frequencies
	void adja::add_ess_to_list(vector <essfreq> & list_ess, essfreq & ess_new)
	{
		unsigned l,m;
		bool ess_is_there=false;

		if(list_ess.size()==0)
		 {
			 list_ess.push_back(ess_new);
			 ess_is_there=true;
		 }
		 else
		 {
			 l=0;
			 ess_is_there=false;

			 while(ess_is_there==false && l<list_ess.size())
			 {
				 if(ess_new.nedges==list_ess[l].nedges)
				 {
					 if(ess_new.sceleton==list_ess[l].sceleton)
					 {
						 if((ess_new.immoral).size() == (list_ess[l].immoral).size())
						 {
							 ess_is_there=true;
							 for(m=0; m<(ess_new.immoral).size(); m++)
							 {
								 if((ess_new.immoral)[m] != (list_ess[l].immoral)[m])
									  ess_is_there=false;
							 }

							 if(ess_is_there==true)
								list_ess[l].freq = list_ess[l].freq + ess_new.freq;
						  }
					 }
				 }
				 l++;

			 }
			 if( ess_is_there==false)
					 list_ess.push_back(ess_new);

		 }
	}



	// FUNCTION: add_to_mean
	// TASK: takes matrix_new into account for the matrix_mean
	void adja::add_to_mean(const datamatrix & matrix_new, datamatrix & matrix_mean, unsigned n)
	{
		assert(matrix_new.cols() == matrix_mean.cols());
		assert(matrix_new.rows() == matrix_mean.rows());
		assert(matrix_new.rows() == nvar);
		assert(matrix_new.cols() == nvar);

		unsigned i,j;

		for(i=0; i<nvar; i++)
		{
			for(j=0; j<nvar; j++)
			{
				if(n==0)
					matrix_mean(i,j) = matrix_new(i,j);
				else
					matrix_mean(i,j) = (n)* matrix_mean(i,j) /(n+1)  + matrix_new(i,j)/(n+1);
			}
		}
	}


	// FUNCTION: add_to_mean
	// TASK: takes matrix_new into account for the matrix_mean
	void adja::add_to_square(const datamatrix & matrix_new, datamatrix & sum_square, unsigned n)
	{
		unsigned i,j;

		for(i=0; i<nvar; i++)
		{
			for(j=0; j<nvar; j++)
			{
				if(n==0)
					sum_square(i,j) = matrix_new(i,j);
				else
					sum_square(i,j)  = sum_square(i,j)  + matrix_new(i,j)*matrix_new(i,j);
			}
		}
	}



	// FUNCTION: write_out_ess_short
	// TASK: writes out the essential graph ess into a separate file
	void adja::write_out_ess_short(essfreq & ess, ST::string path_res, unsigned n)
	{
//		int i,j;
		int num_imm;

		ofstream fout(path_res.strtochar(),ios::app);

		for(unsigned i=0; i<nvar; i++)
		{
			for(unsigned j=0; j<nvar; j++)
			{
				fout<<ess.sceleton(i,j);
			}
		}

		fout<<endl;

        num_imm = (ess.immoral).size();
		fout<<num_imm<<"\t";

		if(num_imm>0)
		{
			for(int i=num_imm-1; i>=0; i--)
			{
				fout<<(ess.immoral)[i][0]<<" ";
				fout<<(ess.immoral)[i][1]<<" ";
				fout<<(ess.immoral)[i][2]<<"\t";

			}
		}
		fout<<endl;

		double rel_freq;
		rel_freq = ess.freq/n;

		fout<<rel_freq<<endl;
		fout.close();
	}




//}  // namespace std











