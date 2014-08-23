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
#include <vcl.h>
#pragma hdrstop

#include<StatwinFrame.h>
#include<statwin_haupt.h>
#endif

#include "graph.h"

using std::endl;
using std::ifstream;
using std::ofstream;

bool graph::isedgeof(unsigned i,unsigned j)
  {

  unsigned k=0;
  bool found = false;
  while (k < edges[j].size() && (found==false) )
    {
    if (edges[j][k] == i)
      found = true;

    k++;

    }

  return found;

  }


  // TASK: tests wheater node i is an edge of node j
bool graph::issymmetric(unsigned & node,unsigned & edge)
  {
  unsigned i,j;

  bool correct = true;

  i=0;

  while ( (i<nodes.size()) && (correct == true) )
    {
    j=0;
    while ( (j< edges[i].size()) && (correct == true) )
      {
      if (isedgeof(i,edges[i][j]) == false)
        {
        correct = false;
        node = i;
        edge = edges[i][j];
        }
      j++;
      }

    i++;
    }


  return correct;


  }


void graph::take(vector<unsigned> & alreadyvis,unsigned node)
  {
  if (alreadyvis[node] == 0)
    {
    alreadyvis[node] = 1;

    unsigned i;

    for (i=0;i<edges[node].size();i++)
      {
      mask[edges[node][i]] = 1;
      connnodes[edges[node][i]] = nrgraphs;
      take(alreadyvis,edges[node][i]);
      }

    }

  }

int graph::findfirstzero(void)
  {
  unsigned i=0;
  int f = -1;
  while ( (i<mask.size()) && (f<0) )
    {
    if (mask[i] == 0)
      f = i;
    i++;
    }

  return f;
  }


void graph::checkconnectivity(void)
  {

  nrgraphs = 1;

  mask[0] = 1;

  vector<unsigned> h(nodes.size(),0);

  connnodes[0] = 1;

  take(h,0);

  int firstzero = findfirstzero();

  while(firstzero > 0)
    {
    nrgraphs++;
    mask[firstzero] = 1;
    connnodes[firstzero] = nrgraphs;
    take(h,firstzero);
    firstzero = findfirstzero();
    }

  unsigned i;
  for (i=0;i<mask.size();i++)
    mask[i] = 0;

/*
  ofstream out("c:\\tmp\\connodes.raw");
  for (i=0;i<nodes.size();i++)
    out << nodes[i] << "   " << connnodes[i] << endl;
*/
  }


// CONSTRUCTOR 1

graph::graph(const ST::string & path)
  {
  ifstream in(path.strtochar());
  assert(!in.fail());

  bool weightex = false;

  linenr=1;
  ST::string help, help2;
  vector<ST::string> token;

  long helplong;
  double helpdouble;
  long nrnodes,edgesi;
  nrnodes=0;

  ST::getline(in,help);
  help = help.eatallcarriagereturns();

  error = help.strtolong(nrnodes);

  if (error == false)
    {

    nodes = vector<ST::string>(nrnodes);
    edges = vector< vector<unsigned> >(nrnodes);
    weights = vector< vector<double> >(nrnodes);

    unsigned i,j,k;

    unsigned weightcheck=0;

    i=0;
    while ( (i<nrnodes) && (error != true) )
      {

      linenr++;
      ST::getline(in,help2);          // reading node i
      help2 = help2.eatallcarriagereturns();
      nodes[i] = help2;

      linenr++;
      ST::getline(in,help);              // reading number of edges
      help = help.eatallcarriagereturns();
      error = help.strtolong(edgesi);
      if (error==false)
        {
        edges[i] = vector<unsigned>(edgesi);
        weights[i] = vector<double>(edgesi);
        linenr++;
        ST::getline(in,help);
        help = help.eatallcarriagereturns();
        token = help.strtoken(" ");

        if( (i==weightcheck) && (token.size()==2*edgesi) )
          {
          if( edgesi > 0 )
            weightex = true;
          else
            weightcheck++;
          }
//        if ( (i==0) && (token.size() == 2*edgesi))
//          weightex=true;

        if (weightex == false)
          {

          if (token.size() == edgesi)
            {

            j=0;
            while ( (j<edgesi) && (error==false) )
              {
              error = token[j].strtolong(helplong);
              if ( (error==false) && (helplong >= 0) && (helplong < nrnodes) )
                {
                edges[i][j] = helplong;
                weights[i][j] = 1.0;
                j++;
                }
              else
                error = true;

              }

            if (error==false)
              i++;

            } // end: if (token.size() == edgesi)
          else
            error = true;

          } // end: if (weightex == false)
        else
          {

          if (token.size() == edgesi*2)
            {

            j=0;
            while ( (j<edgesi) && (error==false) )
              {
              error = token[j].strtolong(helplong);
              if ( (error==false) && (helplong >= 0) && (helplong < nrnodes) )
                {
                edges[i][j] = helplong;
                j++;
                }
              else
                error = true;

              } // end: while ...

            k=0;
            while ( (j<edgesi*2) && (error==false) )
              {
              error = token[j].strtodouble(helpdouble);
              if ( (error==false) && (helpdouble > 0) )
                {
                weights[i][k] = helpdouble;
                j++;
                k++;
                }
              else
                error = true;

              }


            if (error==false)
              i++;

            } // end: if (token.size() == edgesi)
          else
            error = true;

          } // weightex = true


        } // end: if (error == false)

      }

    mask = vector<unsigned>(nrnodes,0);

    }

  if (error == false)
    {

    unsigned n;
    unsigned e;

    bool sym = issymmetric(n,e);

    if (sym == false)
      {
      errormessage = "ERROR: graph is not symmetric. node " +
       ST::inttostring(e) + " is an edge of node " + ST::inttostring(n) +
       " but not vice versa\n";
      error = true;
      }
    else
      {
      connnodes = vector<unsigned>(nodes.size(),0);
      checkconnectivity();
      }


    }
  else
    {
    errormessage =  "ERROR: graph file invalid in line " +
      ST::inttostring(linenr) + "\n";
    }


  if (error == true)
    {
    nodes = vector<ST::string>();
    edges = vector< vector<unsigned> >();
    }

  }


graph::graph(const vector<ST::string> & no,
             const vector< vector<unsigned> > & ed)
  {
  nodes = no;
  edges = ed;
  unsigned i;
  weights = vector< vector<double> >(nodes.size());
  for(i=0;i<nodes.size();i++)
    weights[i] = vector<double>(edges[i].size(),1.0);
  mask = vector<unsigned>(nodes.size(),0);

  unsigned n;
  unsigned e;
  bool error = false;

  bool sym = issymmetric(n,e);

  if (sym == false)
    {
    errormessage = "ERROR: graph is not symmetric. node " +
    ST::inttostring(e) + " is an edge of node " + ST::inttostring(n) +
    " but not vice versa\n";
    error = true;
    }
  else
    {
    connnodes = vector<unsigned>(nodes.size(),0);
    checkconnectivity();
    }

  if (error == true)
    {
    nodes = vector<ST::string>();
    edges = vector< vector<unsigned> >();
    }

  }


graph::graph(const vector<ST::string> & no,
             const vector< vector<unsigned> > & ed,
             const vector< vector<double> > & we)
  {
  nodes = no;
  edges = ed;
  weights = we;
  mask = vector<unsigned>(nodes.size(),0);

  unsigned n;
  unsigned e;
  bool error = false;

  bool sym = issymmetric(n,e);

  if (sym == false)
    {
    errormessage = "ERROR: graph is not symmetric. node " +
    ST::inttostring(e) + " is an edge of node " + ST::inttostring(n) +
    " but not vice versa\n";
    error = true;
    }
  else
    {
    connnodes = vector<unsigned>(nodes.size(),0);
    checkconnectivity();
    }

  if (error == true)
    {
    nodes = vector<ST::string>();
    edges = vector< vector<unsigned> >();
    }

  }



graph::graph(const graph & g)
  {
  nodes = g.nodes;
  edges = g.edges;
  weights = g.weights;
  mask = g.mask;
  error = g.error;
  linenr = g.linenr;
  errormessage = g.errormessage;
  nrgraphs = g.nrgraphs;
  connnodes = g.connnodes;
  }

const graph & graph::operator=(const graph & g)
  {
  if (this == &g)
	 return *this;
  nodes = g.nodes;
  edges = g.edges;
  weights = g.weights;
  mask = g.mask;
  error = g.error;
  linenr = g.linenr;
  errormessage = g.errormessage;
  nrgraphs = g.nrgraphs;
  connnodes = g.connnodes;
  return *this;
  }


ostream & operator<<(ostream & c, const graph & g)
  {
  c << "Number of nodes: " << g.nodes.size() << endl;

  return c;
  }


void graph::outgraph(const ST::string & path,const bool & we) const
  {
  ofstream out(path.strtochar());
  unsigned i,j;
  out << nodes.size() << endl;
  for(i=0;i<nodes.size();i++)
    {
    out << nodes[i] << endl;
    out << edges[i].size() << endl;
    for (j=0;j<edges[i].size();j++)
      out << edges[i][j] << "  ";
    if (we == true)
      {
      for (j=0;j<weights[i].size();j++)
        out << weights[i][j] << "  ";
      }

    out << endl;
    }

  }


void graph::outgraph2(const ST::string & path,const bool & we) const
  {
  ofstream out(path.strtochar());
  unsigned i,j;
  out << nodes.size() << endl;
  for(i=0;i<nodes.size();i++)
    {
    out << nodes[i] << endl;
    out << edges[i].size() << endl;
    for (j=0;j<edges[i].size();j++)
      out << nodes[edges[i][j]] << "  ";
    if (we == true)
      {
      for (j=0;j<weights[i].size();j++)
        out << weights[i][j] << "  ";
      }

    out << endl << endl;
    }
  }


vector<unsigned> graph::CM(const unsigned & start)
  {

  vector<unsigned> perm(nodes.size());
  vector<unsigned> permrev(nodes.size());
  unsigned i,j;
  unsigned nr;

  perm[0] = start;
  mask[perm[0]] = 1;

  nr=1;

  for (i=0;i<nodes.size();i++)
    {

    if (nr < nodes.size())
      {

      for (j=0;j<edges[perm[i]].size();j++)
        {
        if (mask[edges[perm[i]][j]] == 0)
          {
          perm[nr] = edges[perm[i]][j];
          mask[edges[perm[i]][j]] = 1;
          nr++;
          }

        }

      }  // end: if (nr < nodes.size())

    }

  for(i=0;i<nodes.size();i++)
    mask[i] = 0;

  for (i=0;i<nodes.size();i++)
    permrev[i] = perm[nodes.size()-1-i];
  return permrev;


//  return perm;
  }


vector<unsigned> graph::CMopt(void)
  {
  graph g = graph(nodes,edges);
  unsigned i;
//  unsigned minsize=nodes.size();
  unsigned minsize = nodes.size()*nodes.size();
  unsigned minind;
  unsigned bs;

  for(i=0;i<nodes.size();i++)
    {
    g = graph(nodes,edges);
    g.reorder(i);
//    bs = g.maxbeta();
    bs = g.sumbeta();
    if (bs < minsize)
      {
      minind = i;
      minsize = bs;
      }
    }


  return CM(minind);
  }


void graph::outindizes(const ST::string & path) const
  {
  ofstream out(path.strtochar());
  unsigned i,j;
  for (i=0;i<nodes.size();i++)
    for(j=0;j<edges[i].size();j++)
      out << i << "  " << edges[i][j] << endl;
  }



void graph::reorder(const unsigned & start)
  {
  vector<unsigned> perm;
  vector<unsigned> invp(nodes.size());

  perm = CM(start);
  unsigned i,j;
  for(i=0;i<nodes.size();i++)
    invp[perm[i]] = i;

  vector<ST::string> nodesnew(nodes.size());
  for(i=0;i<nodes.size();i++)
    nodesnew[i] = nodes[perm[i]];



  vector< vector<unsigned> > edgesnew(nodes.size());


  for(i=0;i<nodes.size();i++)
    {
    edgesnew[i] = vector<unsigned>(edges[perm[i]].size());
    for(j=0;j<edges[perm[i]].size();j++)
      edgesnew[i][j] = invp[edges[perm[i]][j]];
    }


  vector< vector<double> > weightsnew(nodes.size());

  for(i=0;i<nodes.size();i++)
    {
    weightsnew[i] = vector<double>(weights[perm[i]].size());
    for(j=0;j<weights[perm[i]].size();j++)
      weightsnew[i][j] = invp[weights[perm[i]][j]];
    }

  nodes = nodesnew;
  edges = edgesnew;
  weights = weightsnew;

  }


unsigned graph::f(const unsigned & i) const
  {
  unsigned f = i;
  unsigned j;
  for (j=0;j<edges[i].size();j++)
    if ( (edges[i][j] < i) && (edges[i][j] < f) )
      f = edges[i][j];
  return f;
  }


unsigned graph::beta(const unsigned & i) const
  {
  return i-f(i);
  }


unsigned graph::maxbeta(void) const
  {
  unsigned maxbeta = beta(0);
  unsigned i;
  unsigned betacurrent;
  for(i=1;i<nodes.size();i++)
    {
    betacurrent = beta(i);
    if (betacurrent > maxbeta)
      maxbeta=betacurrent;
    }

  return maxbeta;
  }


unsigned graph::sumbeta(void) const
  {
  unsigned sumbeta=0;
  unsigned i;
  for(i=0;i<nodes.size();i++)
    {
    sumbeta += beta(i);
    }

  return sumbeta;
  }


