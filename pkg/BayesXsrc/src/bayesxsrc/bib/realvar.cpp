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





#include "realvar.h"

#include <list>
#include <iterator>

using namespace std;

namespace realob
{

//------------------------------------------------------------------------------
//------------- CLASS realvar: implementation of member functions --------------
//------------------------------------------------------------------------------


realobs realvar::min()
  {
  realobs m = NA;
  if (! empty())
	 {
	 realvar::iterator it;
	 for (it=begin();it != end();++it)
		if ((*it) < m)
		  m = (*it);
	 }
  return m;
  }


realobs realvar::sum()
  {
  bool allNA = true;
  realobs s = 0;
  if (! empty())
	 {
	 realvar::iterator it;
	 for (it=begin();it != end();++it)
		if (*it != NA)
		  {
		  allNA = false;
		  s = s + (*it);
		  }
	 }
  if (allNA == true)
	 return NA;
  else
	 return s;
  }



realvar lag(realvar & v)
  {
  realvar help;
  help.push_back(NA);
  if (v.size() > 1)
	 help.insert(help.end(),v.begin(),v.end()-1);
  return help;
  }




realvar lagrealvar(realvar & v,vector<int> & index)
  {
  realvar help(v.size());
  help[index[0]] = NA;
  unsigned i;
  for (i=1;i<v.size();i++)
	 help[index[i]] = v[index[i-1]];
  return help;
  }



realvar power(const realvar & v1,const realvar & v2)
  {
  realvar h(v1.size());
  for (unsigned i=0;i<v1.size();i++)
	 h[i] = pow(v1[i],v2[i]);
  return h;
  }



/*
realvar realvar::applied(realobs (*func)(realobs &))
  {
  realvar h(size());
  int i;
  for (i=0;i<size();i++)
	 h[i] = func(operator[](i));
  return h;
  }
*/

realvar sqrt(const realvar & v)
  {
  return v.applied(realob::sqrt);
  }

realvar abs(const realvar & v)
  {
  return v.applied(realob::abs);
  }

realvar exp(const realvar & v)
  {
  return v.applied(realob::exp);
  }

realvar cos(const realvar & v)
  {
  return v.applied(realob::cos);
  }

realvar sin(const realvar & v)
  {
  return v.applied(realob::sin);
  }

realvar log(const realvar & v)
  {
  return v.applied(realob::log);
  }

realvar log10(const realvar & v)
  {
  return v.applied(realob::log10);
  }

realvar floor(const realvar & v)
  {
  return v.applied(realob::floor);
  }

realvar cumul(realvar & v,vector<int> & index)
  {
  int j;
  double value;

  vector<int> helpindex = index;

  realvar h(v.size());
  v.sort(helpindex,0,v.size()-1);

  j = v.size()-1;
  int n = v.size();

  while (j >= 0)
	 {
	 if (v[helpindex[j]] == NA)
		{
		n--;
		h[helpindex[j]] == NA;
		}
	 else
		break;
	 j--;
	 }

  if (n > 0)
	 {

	 realobs pred = v[helpindex[0]];
	 int nr = 0;

	 for (int i=0;i<n;i++)
		{
		if (i == n-1)
		  {
		  if (v[helpindex[i]] == pred)
			 value = 1;
		  else
			 value = i/double(n);
		  for (j=i-nr;j<i;j++)
			 h[helpindex[j]] = value;
		  h[helpindex[i]] = 1;
		  }
		else
		  {
		  if (v[helpindex[i]] == pred)
			 nr++;
		  else
			 {
			 value = i/double(n);
			 for (j=i-nr;j<i;j++)
				h[helpindex[j]] = value;
			 nr = 1;
			 pred = v[helpindex[i]];
			 }
		  }
		}
	 return h;

	 }

  else
	 return realvar(v.size(),NA);

  }



realvar uniform(unsigned obs)
  {
  srand((unsigned)time(0));
  realvar h(obs);
  register unsigned i;
  realvar::iterator it = h.begin();
  for (i=0;i<obs;i++,++it)
	 *it = randnumbers::uniform();
  return h;
  }


realvar cumulnorm(realvar & v)
  {
  realvar h(v.size());
  register unsigned i;
  realvar::iterator it = h.begin();
  realvar::iterator vit = v.begin();
  for (i=0;i<v.size();i++,++it,++vit)
	 *it = randnumbers::Phi((*vit).getvalue());
  return h;

  }

realvar normal(unsigned obs)
  {
  srand((unsigned)time(0));
  realvar h(obs);
  register unsigned i;
  realvar::iterator it = h.begin();
  for (i=0;i<obs;i++,++it)
	 *it = randnumbers::rand_normal();
  return h;
  }



realvar exponential(unsigned obs,realobs lambda)
  {
  srand((unsigned)time(0));
  realvar h(obs);

  for (unsigned i=0;i<obs;i++)
	 h[i] = _exponential(lambda);
  return h;
  }



realvar exponential(realvar & lambda)
  {
  srand((unsigned)time(0));
  realvar h(lambda.size());

  realvar::iterator lambdait = lambda.begin();
  realvar::iterator hit = h.begin();


  for (unsigned i=0;i<lambda.size();i++,++hit,++lambdait)
	 {
	 if ( (*lambdait <= 0.0) || (*lambdait==NA) )
		*hit = NA;
	 else
		*hit = _exponential(*lambdait);
	 }
  return h;
  }



realvar bernoulli(realvar & p)
  {
  srand((unsigned)time(0));
  realobs u;
  realvar h(p.size());

  realvar::iterator pit = p.begin();
  realvar::iterator hit = h.begin();

  for(unsigned i=0;i<p.size();i++,++pit,++hit)
	 {
	 if ((*pit > 1.0) || (*pit < 0.0) || (*pit == NA) )
		{
		*hit = NA;
		}
	 else
		{
		u = _uniform();
		if (u <= *pit)
		  *hit = 1;
		else
		  *hit = 0;
		}
	 }  // end: for(int i=0;i<p.size();i++)
  return h;
  }


realvar binomial(realvar & n,realvar & p)
  {
  srand((unsigned)time(0));
  unsigned i,j;
  realobs u;
  realvar h(n.size());

  realvar::iterator pit = p.begin();
  realvar::iterator nit = n.begin();
  realvar::iterator hit = h.begin();

  for (i=0;i<p.size();i++,++hit,++nit,++pit)
	 {

	 if ( (*nit < 1.0) || (*pit > 1.0) || (*pit < 0.0)
          || *nit == NA || *pit == NA)
		{
		*hit = NA;
		}
	 else
		{
		*hit = 0;
		for (j=1;j<=(*nit);j++)
		  {
		  u = _uniform();
		  if (u <= (*pit))
//			 (*hit)++;
			 ++(*hit);
             }
		}
	 }
  return h;
  }


realvar gamma(realvar & mu, realvar & nu)
  {
  srand((unsigned)time(0));
  register unsigned  i;
  realobs u;
  realvar h(mu.size());

  realvar::iterator muit = mu.begin();
  realvar::iterator nuit = nu.begin();
  realvar::iterator hit = h.begin();

  for (i=0;i<mu.size();i++,++hit,++nuit,++muit)
    {

    if ( (*muit <= 0.0) || (*nuit <= 0.0) || (*muit == NA) || (*nuit == NA) )
      {
      *hit = NA;
      }
    else
      {

      *hit = randnumbers::rand_gamma((*nuit).getvalue(),
                                    (*nuit).getvalue() / (*muit).getvalue() );

      }

    }

  return h;

  }


} // end: namespace realobs

