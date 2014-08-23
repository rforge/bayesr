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



#if !defined (MAP_INCLUDED)

#define MAP_INCLUDED

#include"../export_type.h"
#if defined(MICROSOFT_VISUAL)
#include<limits>
#else
#include "../values.h"
#endif
#include<fstream>
#include "clstring.h"
#include<assert.h>
#include<vector>
#include<algorithm>
#include "matrix.h"
#include "statmat.h"
#include "graph.h"

#if defined(JAVA_OUTPUT_WINDOW)
#include"adminparse_basic.h"
#endif


using std::endl;

namespace MAP
{

//#if defined(JAVA_OUTPUT_WINDOW)
 #if defined(MICROSOFT_VISUAL)
 const double NA = DBL_MAX;
 #else
 const double NA = MAXDOUBLE;
 #endif
//#endif
//#if defined(__BUILDING_GNU)
// const double NA = MAXDOUBLE;
//#endif


enum order {names, xvalues, yvalues};
enum metric {adjacent, centroid, combnd};

//------------------------------------------------------------------------------
//------------------------------- struct: line ---------------------------------
//------------------------------------------------------------------------------


// struct 'line' repräsentiert eine Gerade in der 2-Punkt Darstellung

struct line
  {

  double x1;                                     // first  X-Coordinate
  double y1;                                     // first  Y-Coordinate
  double x2;                                     // second X-Coordinate
  double y2;                                     // second Y-Coordinate
  double slope;                                  // slope of the line

  // CONSTRUCTOR
  // TASK: initializes a line with end points (nx1,ny1) and (nx2,ny2)
  //       computes the slope of the line
  //       slope = NA if nx1 = nx2, i.e. if the line is parallel to the
  //       y-axis

  line(double nx1, double ny1, double nx2, double ny2);

  // DEFAULTCONSTRUCTOR

  line(void)
    {
    x1 = 0;
    y1=  0;
    x2 = 1;
    y2 = 1;
    slope = 1;
    }


  // COPY CONSTRUCTOR

  line(const line & l);


  // OVERLOADED ASSIGNMENT OPERATOR

  const line & operator=(const line & l);


  // DESTRUCTOR

  ~line() {}

  // FUNCTION: isparallelto
  // TASK: returns true, if line 'l' is parallel to the calling object

  bool isparallelto(const line & l) const
    {
    if (slope == l.slope)
      return true;
    else
      return false;
    }

   // FUNCTION: isonline
   // TASK: returns true, if x,y is a point of the calling line

  bool isonline(const double & x,const double & y) const;

  // FUNCTION: isinsideline
  // TASK: returns true, if x,y is a inner point of the line

  bool isinsideline(const double x,const double y) const;

  // FUNCTION: isconnected
  // TASK: returns true, if both lines are parallel and in addition have
  //       common pieces

  bool isconnected(const line & l) const;

  // FUNCTION: commonlength
  // TASK: returns length of the common segment of the two lines

  double commonlength(const line & l) const;

  // FUNCTION: print
  // TASK: mainly for testing
  //       writes endpoints and slope to stream object 'o'

  void print(ostream & o) const;


  };


//------------------------------------------------------------------------------
//---------------------------- class: polygone --------------------------------
//------------------------------------------------------------------------------

// class polygone repräsentiert einen Polygonzug

class polygone
  {


  private:


  unsigned nrlines;                     // number of lines of the polygone
  vector<line> lines;                   // vector of lines

  double xmin;                     // minimum x coordinate
  double xmax;                     // maximum x coordinate
  double ymin;                     // minimum y coordinate
  double ymax;                     // maximum y coordinate

  void compute_min_max(void);

  public:


  //----------------------- CONSTRUCTORS AND DESTRUCTOR ------------------------


  // DEFAULT CONSTRUCTOR

  polygone(void);

  // CONSTRUCTOR

  polygone(const vector<line> & l);

  // COPY CONSTRUCTOR

  polygone(const polygone & p);

  // OVERLOADED ASSIGNMENT OPERATOR

  const polygone & operator=(const polygone & p);

  // DESTRUCTOR

  ~polygone() {}


//--------------------- ACCESS TO POLYGONE CHARAKTERISTICS ---------------------

  // FUNCTION: get_nrlines
  // TASK: returns the number of lines of the polygone

  const unsigned & get_nrlines(void) const
    {
    return nrlines;
    }

  // FUNCTION: get_line
  // TASK: returns line 'i'

  const line & get_line(const unsigned & i) const
    {
    assert (lines.size() > 0);
    assert (i < lines.size());
    return lines[i];
    }

  // FUNCTION: add_line
  // TASK: adds a line to the current polygone

  void add_line(const line & l);

  // FUNCTION: commonborderlength
  // TASK: returns the length of the common border segment of the two polygones

  double commonborderlength(const polygone & p) const;

  // FUNCTION: get_xmin
  // TASK: returns xmin

  double get_xmin(void)
    {
    return xmin;
    }

  // FUNCTION: get_xmax
  // TASK: returns xmax

  double get_xmax(void)
    {
    return xmax;
    }

  // FUNCTION: get_ymin
  // TASK: returns ymin

  double get_ymin(void)
    {
    return ymin;
    }
  // FUNCTION: get_ymax
  // TASK: returns ymax

  double get_ymax(void)
    {
    return ymax;
    }

  };


//------------------------------------------------------------------------------
//------------------------------ class: region ---------------------------------
//------------------------------------------------------------------------------


class region
  {


  private:

  order orderrelation;             // order relation according to which
                                   // the sorting is done

  double xcenter;                  // x-coordinate of the centroid
  double ycenter;                  // y-coordinate of the centroid


  ST::string name;                     // name of the region
  vector<polygone> polygones;      // vector of polygones the region consists of
  unsigned nrpoly;                 // number of polygones

  ST::string isin;                     // isin != "" indicates, that the region is
                                   // completely surrounded by a larger
                                   // region; isin stores the name of that
                                   // region

  double xmin;                     // minimum x coordinate
  double xmax;                     // maximum x coordinate
  double ymin;                     // minimum y coordinate
  double ymax;                     // maximum y coordinate

  void compute_min_max(void);

  public:


  // FUNCTION: x_center
  // TASK: computes x-coordinate of the centroid and sets xcenter

  void x_center(void);

  // FUNCTION: y_center
  // TASK: computes y-coordinate of the centroid and sets ycenter

  void y_center(void);

  void set_center(const double & x, const double & y);

  const double & get_xcenter(void) const
    {
    return xcenter;
    }


  const double & get_ycenter(void) const
    {
    return ycenter;
    }


  // FUNCTION: set_order
  // TASK: sets the order relation according to which the regions are sorted

/*  void set_order(const order & or)
    {
    orderrelation = or;
    }*/

  // FUNCTION: order
  // TASK: returns the order relation according to which the regions are sorted

  const order & get_order(void) const
    {
    return orderrelation;
    }


  //--------------------- CONSTRUCTORS AND DESTRUCTOR --------------------------

  // DEFAULT CONSTRUCTOR

  region(void)
    {
    nrpoly=0;
    orderrelation = names;
    }

  // CONSTRUCTOR 0
  // TASK: initializes a region with no polygones

  region(const ST::string & n)
    {
    nrpoly=0;
    name=n;
    orderrelation = names;
    }

  // CONSTRUCTOR 1
  // TASK: initializes a region that consists of just one polygone 'p'

  region(const ST::string & n,const polygone & p)
    {
    name = n;
    polygones.push_back(p);
    nrpoly = 1;
    orderrelation = names;
    compute_min_max();
    }

  // CONSTRUCTOR 2
  // TASK: initializes a region that consists of the polygones stored in vector
  //       'p'

  region(const ST::string & n,const vector<polygone> & p)
    {
    name = n;
    polygones = p;
    nrpoly = p.size();
    orderrelation = names;
    compute_min_max();
    }


  // COPY CONSTRUCTOR

  region(const region & r);

  // OVERLOADED ASSIGNMENT OPERATOR

  const region & operator=(const region & r);


  // DESTRUCTOR

  ~region() {}


  // OVERLOADED << OPERATOR

  friend ostream & operator<<(ostream & c, const region & r)
    {
    c << "region name:                             " << r.get_name() << endl;
    c << "order definition:                        " << r.orderrelation << endl;
    c << "nrpoly:                                  " << r.nrpoly << endl;
    c << "xcenter:                                 " << r.xcenter << endl;
    c << "ycenter:                                 " << r.ycenter << endl << endl;
    return c;
    }


  // OVERLOADED < OPERATOR

  int operator<(const region & p2) const
    {
    if(orderrelation == names)
      return name < p2.name;
    else if(orderrelation == xvalues)
      return xcenter < p2.xcenter;
    else if(orderrelation == yvalues)
      return ycenter < p2.ycenter;
    else return 0;
    }


//------------------ GETTING ACCESS TO REGION CHARACTERISTICS ------------------

  // FUNCTION: get_name
  // TASK: returns the name of the region

  const ST::string & get_name(void) const
    {
    return name;
    }

  // FUNCTION: set_name
  // TASK: sets the name of the region to 'n'

  void set_name(const ST::string & n)
    {
    name = n;
    }

  // FUNCTION: get_polygones
  // TASK: returns the vector of polygones the region consists of

  const vector<polygone> & get_polygones(void) const
    {
    return polygones;
    }

  // FUNCTION: get_polygone
  // TASK: returns polygone 'nr'

  const polygone & get_polygone(const unsigned & nr) const
    {
    assert(polygones.size() > 0);
    assert(nr < polygones.size());
    return polygones[nr];
    }

  // FUNCTION: get_nrpoly
  // TASK: returns the number of polygones the region consists of

  const unsigned & get_nrpoly(void) const
    {
    return nrpoly;
    }

  // FUNCTION: add_polygone
  // TASK: adds a polygone to the region

  void add_polygone(const polygone & p)
    {
    polygones.push_back(p);
    nrpoly++;
    compute_min_max();
    }

  // FUNCTION: set_isin
  // TASK: initializes string 'isin' (see above)

  void set_isin(const ST::string & is)
    {
    isin = is;
    }

  const ST::string get_isin(void) const
    {
    return isin;
    }

  // FUNCTION: compare
  // TASK: returns true, if the two regions are direct neigbors, i.e. have a
  //       common boundary

  bool compare(const region & r) const;


  // FUNCTION: distance
  // TASK: computes the distance between two regions according to metric m

  double distance(const region & r, const metric & m) const;

  // FUNCTION: get_xmin
  // TASK: returns xmin

  double get_xmin(void)
    {
    return xmin;
    }

  // FUNCTION: get_xmax
  // TASK: returns xmax

  double get_xmax(void)
    {
    return xmax;
    }

  // FUNCTION: get_ymin
  // TASK: returns ymin

  double get_ymin(void)
    {
    return ymin;
    }
  // FUNCTION: get_ymax
  // TASK: returns ymax

  double get_ymax(void)
    {
    return ymax;
    }

  };


//------------------------------------------------------------------------------
//---------------------------- CLASS: map --------------------------------------
//------------------------------------------------------------------------------


class __EXPORT_TYPE map
  {


  private:

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adminb_p;
  #endif

  metric weightmode;                       // metric according to which the
                                           // weights are computed

  vector<ST::string> errormessages;        // errormessages for errors occuring
                                           // while reading boundary information

  vector<region> regions;                  // vector of regions
  bool nopolygones;
  bool nocentroids;
  unsigned nrregions;                      // number of regions


  double minX;                           // minimum X coordinate
  double maxX;                           // maximum X coordinate
  double minY;                           // minimum Y coordinate
  double maxY;                           // maximum Y coordinate


  vector < vector<unsigned> > neighbors;   // neighbors of the regions
  unsigned minn;                           // Minimum number of neighbors
  unsigned maxn;                           // Maximum number of neighbors

  unsigned bandsize;


  vector < vector<double> > weights;       // weights of the neigbors
  double mindistance;
  double maxdistance;

  // FUNCTION: identify_regions
  // TASK: drops regions, that appear more than once in the map. boundaries of
  //       such regions will be stored in the remaining region.

  void identify_regions(void);


  // FUNCTION: infile
  // TASK: reads in map information

  void infile(const ST::string & path);

  // FUNCTION: infile_neighbors
  // TASK: reads in neighborhood information
  //       neighborhood information must be given in index form

  void infile_neighbors(const ST::string & path);

  // FUNCTION: compute_minmaxn
  // TASK: computes minimum and maximum number of neighbors
  //       (stored in 'minn' and 'maxn')

  void compute_minmaxn(void);

  void reset(void);


  public:


  //----------------------------------------------------------------------------
  //------------ CONSTRUCTORS, DESTRUCTOR AND OVERLOADED OPERATORS -------------
  //----------------------------------------------------------------------------

  // DEFAULT CONSTRUCTOR

  map(void)
    {
    nrregions =  0;
    minn = 0;
    maxn = 0;
    nopolygones = true;
    nocentroids = true;
    }

  // CONSTRUCTOR 1
  // TASK: reads the polygones stored in file 'path'
  //       computes neigbors (first order) and
  //       weights = number of neighbors

  map(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * abp,
  #endif
  const ST::string & path,const metric & m /*= adjacent*/);

  // CONSTRUCTOR 2
  // TASK: reads the polygones stored in file 'path'
  //       reads the neighbors stored in file 'npath'
  //       weights = number of neighbors

  map(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * abp,
  #endif
  const ST::string & bpath,const ST::string & npath,
      const metric & m = adjacent);

  // CONSTRUCTOR 3
  // TASK: reads in the values of two variables x and y stored in path
  //

  map(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * abp,
  #endif
  const ST::string & path);

  // CONSTRUCTOR 4
  // TASK: reads in the values of two variables x and y stored in path
  //       creates a map, neighbors of a particular point are all adjacent
  //       points within a distance 'md' from the point
  // IMPORTANT: metric=combnd is invalid

  map(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * abp,
  #endif
  const datamatrix & x,const double & md, const metric & m);

  // CONSTRUCTOR 5
  // TASK: reads in the stored in graph 'g'

  map(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * abp,
  #endif
  const graph & g);

  // COPY CONSTRUCTOR

 map(const map & m);

 // OVERLOADED ASSIGNMENT OPERATOR

 const map & operator=(const map & m);

  // DESTRUCTOR

  ~map(void) {}

  // OVERLOADED << OPERATOR
  // TASK: writes charakteristics of the map 'm' to the stream 'c'

  friend ostream & operator<<(ostream & c, const map & m);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //------------------------ SORTING AND REORDERING ----------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: sortmap
  // sorts the regions according to their names

  void sortmap(void)
    {
    std::sort(regions.begin(),regions.end());
    }

  // FUNCTION: reorder
  // reorders the regions of the map according to the ordering in 'path'
  // the file must contain a list of numbers where the i'th number 'inr'
  // in the list indicates that the new ith region is the inr th region in the
  // old order

  void reorder(const ST::string & path);

  // FUNCTION: reorderopt
  // TASK: aims at finding a optimal ordering with minimal bandwidth of
  //       the resulting adjacency matrix

  void reorderopt(void);

  // FUNCTION: isconnected
  // TASK: returns true, if map is connected, false else

  bool isconnected(void);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //--- functions for computing neighbors and getting access to neighbor info --
  //----------------------------------------------------------------------------

  // FUNCTION: outneighbors
  // TASK: writes neighbors to file 'path'
  //       first column represents the region name followed by the neighbors of
  //       that regions
  //       if names = true, neighbors will be written as regions names
  //       if names = false, neighbors will be written as indices of regions
  //                         (index starts with 0)

  void outneighbors(const ST::string & path,const bool & names) const;

  void outneighbors2(const ST::string & path) const;

  // FUNCTION: getneighbors
  // TASK: für Guenter

  matrix getneighbors(void) const;

  // FUNCTION: computeneighbors
  // TASK: computes the neighbors and stores the result in 'neighbors'
  //       computes weights

  void computeneighbors(void);

  // FUNCTION: get_neighbors
  // TASK: returns the neighbors

  const vector< vector<unsigned> > & get_neighbors(void) const
    {
    return neighbors;
    }

  // FUNCTION: get_minn
  // TASK: returns the minimum number of neighbors

  const unsigned & get_minn(void) const
    {
    return minn;
    }

  // FUNCTION: get_maxn
  // TASK: returns the maximum number of neighbors

  const unsigned get_maxn(void) const
    {
    return maxn;
    }

  //----------------------------------------------------------------------------
  //------------------ end: functions for computing neigbors -------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //-------------- functions for computing and accessing weights ---------------
  //----------------------------------------------------------------------------

  // FUNCTION: compute_weights
  // TASK: computes weights for neighbors
  //       weights = 1/number of neighbors

  void compute_weights(const metric);

  // FUNCTION: out_weights
  // TASK: TASK: writes weights to file 'path'

  void out_weights(const ST::string & path) const;

  // FUNCTION: get_weights
  // TASK: returns weights

  const vector< vector<double> > & get_weights(void) const
    {
    return weights;
    }

  // FUNCTION: get_mindistance
  // TASK: returns

  const double & get_mindistance(void) const
    {
    return mindistance;
    }

  const double & get_maxdistance(void) const
    {
    return maxdistance;
    }

  // FUNCTION: get_weightssum
  // TASK: returns the sum of weights for the 'i' th region

  double get_weightssum(const unsigned & i) const;

  //----------------------------------------------------------------------------
  //------------------- end: functions for computing weights -------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //----------------------- access to map characteristics ----------------------
  //----------------------------------------------------------------------------

  bool polygones_existing(void) const
    {
    return !nopolygones;
    }

  bool centroids_existing(void) const
    {
    return !nocentroids;
    }

  // FUNCTION: outmap
  // TASK: writes the map to file 'path'

  void outmap(const ST::string & path) const;

  // FUNCTION: outgraph
  // TASK: writes the graph assiciated with the map to file 'path'

  void outgraph(const ST::string & path,const bool & we=false) const;

  // FUNCTION get_graph
  // TASK: returns the graph associated with the map

  graph get_graph(void) const;


  void outindices (const ST::string & path) const;


  // FUNCTION: get_nrregions
  // TASK: returns the number of polygones currently in memory

  const unsigned & get_nrregions(void) const
    {
    return nrregions;
    }


  // FUNCTION: drop_region
  // TASK: drops region 'nr'
  //       recomputes neighbors and weights

  void drop_region(const unsigned & nr);


  // FUNCTION: getname
  // TASK: returns the name of the 'nr'

  const ST::string & getname(const unsigned & nr) const
    {
    assert (regions.size() > 0);
    assert (nr < regions.size());
    return regions[nr].get_name();
    }


  // FUNCTION: get_region
  // TASK: returns region 'nr'

  const region & get_region(const unsigned & nr) const
    {
    assert(regions.size() > 0);
    assert(nr < regions.size());
    return regions[nr];
    }


  // FUNCTION: get minX
  // TASK: returns minimum X-Coordinate

  const double & get_minX(void) const
    {
    return minX;
    }

  // FUNCTION: get maxX
  // TASK: returns maximum X-Coordinate

  const double & get_maxX(void) const
    {
    return maxX;
    }

  // FUNCTION: get minY
  // TASK: returns minimum Y-Coordinate

  const double & get_minY(void) const
    {
    return minY;
    }

  // FUNCTION: get maxY
  // TASK: returns maximum Y-Coordinate

  const double & get_maxY(void) const
    {
    return maxY;
    }

  // FUNCTION: getnr
  // TASK: returns number of the region 'name' if existing, and -1 otherwise

  const int getnr(const ST::string & name) const
    {
    int i = 0;
    while(i<int(regions.size()) && !(regions[i].get_name()==name))
      i++;
    return i;
    }

  // FUNCTION: get_bandsize
  // TASK: returns the maximum bandisze of the map

  const unsigned & get_bandsize(void) const
    {
    return bandsize;
    }

  // FUNCTION: outcentroids
  // TASK: returns the centroid of each region

  void outcentroids(const ST::string & path);


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //-------------- COMPUTING INDIZES, ETC. FOR REGRESSION MODELS ---------------
  //----------------------------------------------------------------------------

  // FUNCTION: compute_reg
  // TASK:
  // effectvalues = vector of regions in the map in the order of appearence
  // posbeg(j) - posend(j) gives the "position" of the j-th value in effectvalue
  // in the data
  // posbeg(j) = -1 denotes value not in d
  // i: posbeg(j) - posend(j): d(index(i,0)) gives the values in d with
  // effectvalue(j)

  void compute_reg(const datamatrix & d,vector<int> & posbeg,
                   vector<int> & posend,vector<ST::string> & effectvalues,
                   statmatrix<int> & index);

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  //----------------------------------------------------------------------------
  //------------------------- ACCESS TO ERRORMESSAGES --------------------------
  //----------------------------------------------------------------------------

  // FUNCTION: get_errormessages

  const vector<ST::string> & get_errormessages(void)  const
    {
    return errormessages;
    }

  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------


  };


} // end: namespace map


#endif


