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



#if !defined (GRAPH_included)

#define GRAPH_included

#include"../export_type.h"
#include"clstring.h"
#include<vector>


class __EXPORT_TYPE graph

  {

  protected:

  vector<ST::string> nodes;                   // nodes of the Graph
  vector< vector<unsigned> > edges;           // edges of the Graph
  vector< vector<double> > weights;           // weights of the edges
  vector <unsigned> mask;                     // helping vector

  unsigned nrgraphs;
  vector<unsigned> connnodes;

  bool error;
  unsigned linenr;
  ST::string errormessage;


  // FUNCTION: isedgeof
  // TASK: tests wheather node i is an edge of node j
  //       returns true if correct

  bool isedgeof(unsigned i,unsigned j);

  // FUNCTION: issymmetric
  // TASK: tests if the graph is symmetric, function stops and returns false
  //       if the first incorrect pair of nodes is detected
  //       in that case 'node' has an edge 'edge' but 'edge' has not an edge
  //       'node'

  bool issymmetric(unsigned & node,unsigned & edge);

  void take(vector<unsigned> & alreadyvis,unsigned  node);

  int findfirstzero(void);

  void checkconnectivity(void);

  public:



  // DEFAULT CONSTRUCTOR

  graph(void)
    {
    }

  // CONSTRUCTOR 1
  // TASK: reads in the Graph stored in file path

  graph(const ST::string & path);

  // CONSTRUCTOR 2
  // TASK: constructs a graph out of the nodes 'no' and the edges 'ed'
  //       weights are set equal to 1

  graph(const vector<ST::string> & no,
        const vector< vector<unsigned> > & ed);

  // CONSTRUCTOR 3
  // TASK: constructs a graph out of the nodes 'no' and the edges 'ed'

  graph(const vector<ST::string> & no,
        const vector< vector<unsigned> > & ed,
        const vector< vector<double> > & we);

  // COPY CONSTRUCTOR

  graph(const graph & g);

  // OVERLOADED ASSIGNMENT CONSTRUCTOR

  const graph & operator=(const graph & g);

  // OVERLOADED << OPERATOR
  // TASK: writes charakteristics of the graph  'g' to the stream 'c'

  friend ostream & operator<<(ostream & c, const graph & g);

  // DESTRUCTOR

  ~graph() {}


  // FUNCTION: geterror
  // TASK: returns variable error
  //

  const bool & geterror(void)
    {
    return error;
    }

  // FUNCTION: getlinenr
  // TASK: returns the line number where an error occured

  const unsigned & getlinenr(void)
    {
    return linenr;
    }

  const ST::string & geterrormessage(void)
    {
    return errormessage;
    }


  // FUNCTION: outgraph
  // TASK: writes the graph to file 'path'
  // FORMAT:
  //  number of nodes
  //  node 1
  //  edges of node 1 (indizes)
  //  .
  //  .
  //  .

  void outgraph(const ST::string & path,const bool & we=false) const;

  // FUNCTION: outgraph2
  // TASK: writes the graph to file 'path'
  // FUNCTION: outgraph
  // TASK: writes the graph to file 'path'
  // FORMAT:
  //  number of nodes
  //  node 1
  //  edges of node 1 (NAMES)
  //  .
  //  .
  //  .

  void outgraph2(const ST::string & path,const bool & we=false) const;

  // FUNCTION: CM
  // TASK: finds a new order of the graph using the Cuthill Mc Kee algorithm
  //       starting node is 'start'
  //       returns the permutation vector of the new order

  vector<unsigned> CM(const unsigned & start);

  // FUNCTION: CMopt
  // TASK: finds a new order of the graph using the Cuthill Mc Kee algorithm
  //       with optimal starting node
  //       returns the permutation vector of the new order

  vector<unsigned> CMopt(void);

  // FUNCTION: reorder
  // TASK: reorder the graph using the Cuthill Mc Kee algorithm
  //       starting node is 'start'

  void reorder(const unsigned & start);

  unsigned f(const unsigned & i) const;

  // FUNCTION: beta
  // TASK: returns the bandsize of row i

  unsigned beta(const unsigned & i) const;

  // FUNCTION: maxbeta
  // TASK: returns the maximum bandsize of the correponding matrix

  unsigned maxbeta(void) const;

  unsigned sumbeta(void) const;

  void outindizes(const ST::string & path) const;


  //----------------------------------------------------------------------------
  //------------------------ Access to nodes and edges -------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: get_nrnodes
  // TASK: returns the number of nodes

  unsigned get_nrnodes(void) const
    {
    return nodes.size();
    }

  // FUNCTION: get_node
  // TASK: returns node 'i'

  const ST::string & get_node(const unsigned & i) const
    {
    return nodes[i];
    }

  // FUNCTION: get_edges
  // TASK: returns the edges of the graph

  const vector < vector<unsigned> > & get_edges(void) const
    {
    return edges;
    }


  // FUNCTION: get_nrgraphs
  // TASK: returns the number of separate graphs

  const unsigned & get_nrgraphs(void)
    {
    return nrgraphs;
    }

  //----------------------------------------------------------------------------
  //--------------------- END: Access to nodes and edges -----------------------
  //----------------------------------------------------------------------------

  // FUNCTION: get_weights
  // TASK: returns weights

  const vector< vector<double> > & get_weights(void) const
    {
    return weights;
    }

  };



#endif
