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

#include"map.h"

using std::ifstream;
using std::ofstream;

namespace MAP
{

//------------------------------------------------------------------------------
//---------------- struct line: implementation of member functions -------------
//------------------------------------------------------------------------------

line::line(double nx1, double ny1, double nx2, double ny2)
  {

  if (nx1 < nx2)
	 {
	 x1 = nx1;
	 y1 = ny1;
	 x2 = nx2;
	 y2 = ny2;
	 }
  else if (nx1 == nx2)
    {
    x1 = nx1;
    x2 = x1;
    if (ny2 < ny1)
      {
      y1 = ny2;
      y2 = ny1;
      }
    else
      {
      y1 = ny1;
      y2 = ny2;
      }

    }
  else
	 {
	 x1 = nx2;
	 y1 = ny2;
	 x2 = nx1;
	 y2 = ny1;
	 }

  if (x1 == x2)
	 slope = NA;
  else
	 slope = (y2-y1)/(x2-x1);

  }


line::line(const line & l)
  {
  x1 = l.x1;
  y1 = l.y1;
  x2 = l.x2;
  y2 = l.y2;
  slope = l.slope;
  }


const line & line::operator=(const line & l)
  {
  if (this == &l)
    return *this;
  x1 = l.x1;
  y1 = l.y1;
  x2 = l.x2;
  y2 = l.y2;
  slope = l.slope;
  return *this;
  }


bool line::isonline(const double & x,const double & y) const
  {
  bool ison=false;

  if ( (x >= x1) && (x <= x2) )
    {
    if (slope==NA)
      {
      if (y1 <= y2)
        {
        if ( (y >= y1) && (y <= y2) )
          ison = true;
        }
      else
        {
        if ( (y >= y2) && (y <= y1) )
          ison = true;
        }

      }
    else
      {
      if (y == (y1+slope*(x-x1)))
        ison = true;
      }

    } // end: if ( (x >= x1) && (x <= x2) )

  return ison;
  }


bool line::isinsideline(const double x,const double y) const
  {
  bool isin = false;
  if (isonline(x,y) == true)
    {
    if ( ((x1 != x) || (y1 != y)) && ((x2 != x) || (y2 != y)) )
      isin = true;
    }

  return isin;
  }


bool line::isconnected(const line & l) const
  {
  bool c = false;
  if (isparallelto(l) == true)
  {
  if ( (x1 == l.x1) && (y1 == l.y1) && (x2==l.x2) && (y2==l.y2) )
    c = true;
  else if ( (isinsideline(l.x1,l.y1) == true) || (isinsideline(l.x2,l.y2) == true) )
    c = true;
  else if ( (l.isinsideline(x1,y1) == true) || (l.isinsideline(x2,y2) == true) )
    c = true;
  }

  return c;
  }


double line::commonlength(const line & l) const
  {
  double length = 0.0;
  vector<double> segment;
  if (isconnected(l) == true)
    {
    if (isonline(l.x1,l.y1) == true)
      {
      segment.push_back(l.x1);
      segment.push_back(l.y1);
      if (isonline(l.x2,l.y2) == true)
        {
        segment.push_back(l.x2);
        segment.push_back(l.y2);
        }
      else
        {
        segment.push_back(x2);
        segment.push_back(y2);
        }
      length = sqrt((segment[3] - segment[1])*(segment[3] - segment[1])
             + (segment[2] - segment[0])*(segment[2] - segment[0]));
      }
    else if (l.isinsideline(x1,y1) == true)
      {
      segment.push_back(x1);
      segment.push_back(y1);
      if (isonline(l.x2,l.y2) == true)
        {
        segment.push_back(l.x2);
        segment.push_back(l.y2);
        }
      else
        {
        segment.push_back(x2);
        segment.push_back(y2);
        }
      length = sqrt((segment[3] - segment[1])*(segment[3] - segment[1])
             + (segment[2] - segment[0])*(segment[2] - segment[0]));
      }
    }
  return length;
  }


void line::print(ostream & o) const
  {
  o << "X1= " << x1 << endl;
  o << "Y1= " << y1 << endl;
  o << "X2= " << x2 << endl;
  o << "y2= " << y2 << endl;
  o << "slope= " << slope << endl;
  o << endl;
  }


//------------------------------------------------------------------------------
//-------------- class polygone: implementation of member functions ------------
//------------------------------------------------------------------------------

polygone::polygone(void)
  {
  nrlines = 0;
  }


void polygone::add_line(const line & l)
  {
  lines.push_back(l);
  nrlines++;
  if (nrlines == 1)
    {

    if (lines[0].x1 <= lines[0].x2)
      {
      xmin = lines[0].x1;
      xmax = lines[0].x2;
      }
    else
      {
      xmin = lines[0].x2;
      xmax = lines[0].x1;
      }

    if (lines[0].y1 <= lines[0].y2)
      {
      ymin = lines[0].y1;
      ymax = lines[0].y2;
      }
    else
      {
      ymin = lines[0].y2;
      ymax = lines[0].y1;
      }

    }
  else
    {

    unsigned i = nrlines-1;

    if (lines[i].x1 < xmin)
      xmin = lines[i].x1;
    if (lines[i].x2 < xmin)
      xmin = lines[i].x2;

    if (lines[i].x1 > xmax)
      xmax = lines[i].x1;
    if (lines[i].x2 > xmax)
      xmax = lines[i].x2;

    if (lines[i].y1 < ymin)
      ymin = lines[i].y1;
    if (lines[i].y2 < ymin)
      ymin = lines[i].y2;

    if (lines[i].y1 > ymax)
      ymax = lines[i].y1;
    if (lines[i].y2 > ymax)
      ymax = lines[i].y2;

    }

  }


void polygone::compute_min_max(void)
  {
  if (lines[0].x1 <= lines[0].x2)
    {
    xmin = lines[0].x1;
    xmax = lines[0].x2;
    }
  else
    {
    xmin = lines[0].x2;
    xmax = lines[0].x1;
    }

  if (lines[0].y1 <= lines[0].y2)
    {
    ymin = lines[0].y1;
    ymax = lines[0].y2;
    }
  else
    {
    ymin = lines[0].y2;
    ymax = lines[0].y1;
    }

  unsigned i;
  for(i=1;i<lines.size();i++)
    {
    if (lines[i].x1 < xmin)
      xmin = lines[i].x1;
    if (lines[i].x2 < xmin)
      xmin = lines[i].x2;

    if (lines[i].x1 > xmax)
      xmax = lines[i].x1;
    if (lines[i].x2 > xmax)
      xmax = lines[i].x2;

    if (lines[i].y1 < ymin)
      ymin = lines[i].y1;
    if (lines[i].y2 < ymin)
      ymin = lines[i].y2;

    if (lines[i].y1 > ymax)
      ymax = lines[i].y1;
    if (lines[i].y2 > ymax)
      ymax = lines[i].y2;

    }

  }


polygone::polygone(const vector<line> & l)
  {
  lines = l;
  nrlines = lines.size();
  compute_min_max();
  }


polygone::polygone(const polygone & p)
  {
  lines = p.lines;
  nrlines = p.nrlines;
  xmin = p.xmin;
  xmax = p.xmax;
  ymin = p.ymin;
  ymax = p.ymax;
  }


const polygone & polygone::operator=(const polygone & p)
  {
  if (this == &p)
    return *this;
  lines = p.lines;
  nrlines = p.nrlines;
  xmin = p.xmin;
  xmax = p.xmax;
  ymin = p.ymin;
  ymax = p.ymax;
  return *this;
  }


double polygone::commonborderlength(const polygone & p) const
  {
  double length = 0.0;
  for(unsigned i = 0; i < nrlines; i++)
    {
    for(unsigned j = 0; j < p.nrlines; j++)
      {
      if (lines[i].isconnected(p.lines[j]) == true)
        length += lines[i].commonlength(p.lines[j]);
      }
    }
  return length;
  }


//------------------------------------------------------------------------------
//---------------- class region: implementation of member functions ------------
//------------------------------------------------------------------------------


void region::compute_min_max(void)
  {

  xmin = polygones[0].get_xmin();
  xmax = polygones[0].get_xmax();
  ymin = polygones[0].get_ymin();
  ymax = polygones[0].get_ymax();

  unsigned i;
  for(i=1;i<polygones.size();i++)
    {
    if (polygones[i].get_xmin() < xmin)
      xmin = polygones[i].get_xmin();

    if (polygones[i].get_xmax() > xmax)
      xmax = polygones[i].get_xmax();

    if (polygones[i].get_ymin() < ymin)
      ymin = polygones[i].get_ymin();

    if (polygones[i].get_ymax() > ymax)
      ymax = polygones[i].get_ymax();

    }

  }


region::region(const region & r)
  {
  orderrelation = r.orderrelation;
  xcenter = r.xcenter;
  ycenter = r.ycenter;
  name = r.name;
  polygones = r.polygones;
  nrpoly = r.nrpoly;
  isin = r.isin;
  xmin = r.xmin;
  xmax = r.xmax;
  ymin = r.ymin;
  ymax = r.ymax;
  }



const region & region::operator=(const region & r)
  {
  if (this == &r)
    return *this;
  orderrelation = r.orderrelation;
  xcenter = r.xcenter;
  ycenter = r.ycenter;
  name = r.name;
  polygones = r.polygones;
  nrpoly = r.nrpoly;
  isin = r.isin;
  xmin = r.xmin;
  xmax = r.xmax;
  ymin = r.ymin;
  ymax = r.ymax;
  return *this;
  }



void region::x_center(void)
  {
  double center = 0.0;
  int norm = 0;
  for(unsigned i=0;i<nrpoly;i++)
    {
    norm += get_polygone(i).get_nrlines();
    for(unsigned j=0;j<get_polygone(i).get_nrlines();j++)
      center += get_polygone(i).get_line(j).x1 + get_polygone(i).get_line(j).x2;
    }
  if(norm > 0)
    xcenter = center/(2*double(norm));
  }


void region::y_center(void)
  {
  double center = 0.0;
  int norm = 0;
  for(unsigned i=0;i<nrpoly;i++)
    {
    norm += get_polygone(i).get_nrlines();
    for(unsigned j=0;j<get_polygone(i).get_nrlines();j++)
      center += get_polygone(i).get_line(j).y1 + get_polygone(i).get_line(j).y2;
    }
  if(norm > 0)
    ycenter = center/(2*double(norm));
  }


void region::set_center(const double & x, const double & y)
  {
  xcenter = x;
  ycenter = y;
  }


bool region::compare(const region & r) const
  {

  if ( (isin == r.name) || (r.isin == name) )
    return true;

  if ((ymax < r.ymin) || (r.ymax < ymin))
    return false;

  if ((xmax < r.xmin) || (r.xmax < xmin))
    return false;

  bool neighbors = false;
  unsigned i,j,k,l;
  i=0;
  while ( (i<nrpoly) && (neighbors==false) )
    {
    j=0;
    while ( (j<r.get_nrpoly()) && (neighbors==false) )
      {
      l=0;
      while ( (l<polygones[i].get_nrlines()) && (neighbors==false) )
        {
        k=0;
        while ( (k<r.polygones[j].get_nrlines()) && (neighbors==false) )
          {
          if (polygones[i].get_line(l).isconnected(r.polygones[j].get_line(k)))
            neighbors = true;
          k++;
          }

        l++;
        }

      j++;
      }

    i++;
    }

  return neighbors;

  }


double region::distance(const region & r,const metric & m) const
  {

  double dist = 0.0;
  switch(m)
    {
    case adjacent:
      if(compare(r) == true)
        dist = 1.0;
      else
        dist = 0.0;
    break;
    case centroid:
      dist = sqrt( (r.xcenter - xcenter)*(r.xcenter - xcenter)+ (r.ycenter - ycenter)*(r.ycenter - ycenter) );
    break;
    case combnd:
      for(unsigned i=0;i<nrpoly;i++)
        for(unsigned j=0;j<r.nrpoly;j++)
          {
/*          if (isin == r.name)
            dist += polygones[i].commonborderlength(polygones[i]);
          else if (r.isin == name)
            dist += r.polygones[j].commonborderlength(r.polygones[j]);
          else */
            dist += polygones[i].commonborderlength(r.polygones[j]);
          }

    break;
    }
  return dist;
  }


//------------------------------------------------------------------------------
//-------------- class map: implementation of member functions -----------------
//------------------------------------------------------------------------------


map::map(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adbp,
#endif
const ST::string & path,const metric & m)
  {
  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = adbp,
  #endif

  nopolygones = false;
  weightmode = m;
  infile(path);
  if (errormessages.empty() )
    {
    for(unsigned i=0;i<nrregions;i++)
      {
      regions[i].x_center();
      regions[i].y_center();
      }
    nocentroids = false;
    computeneighbors();
    }

  }


map::map(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adbp,
#endif
const ST::string & bpath,const ST::string & npath,const metric & m)
  {
  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = adbp,
  #endif
  nopolygones = false;
  weightmode = m;
  infile(bpath);
  infile_neighbors(npath);
  for(unsigned i=0;i<nrregions;i++)
    {
    regions[i].x_center();
    regions[i].y_center();
    }
  nocentroids = false;
  }


map::map(
  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_basic * adbp,
  #endif
  const ST::string & path)
  {

  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = adbp,
  #endif

  nopolygones = true;

//  weightmode = adjacent;

  ifstream in(path.strtochar());
  double x,y;
  ST::string name;

  in >> nrregions;
  for(unsigned i=0;i<nrregions;i++)
    {
    in >> name;
    in >> x;
    in >> y;

    regions.push_back(region(name));
    regions[i].set_center(x,y);
    }

  nocentroids = false;

  if (errormessages.empty() )
    {
    computeneighbors();
    }


  }


map::map(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adbp,
#endif
const graph & g)
  {

  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = adbp;
  #endif


  unsigned i;

  nocentroids = true;
  nrregions = g.get_nrnodes();
  regions = vector<region>(g.get_nrnodes());
  for(i=0;i<nrregions;i++)
     regions[i].set_name(g.get_node(i));
  nopolygones = true;

  neighbors = g.get_edges();
  compute_minmaxn();
  bandsize = g.maxbeta();

  //--------------------------- computing weights ------------------------------

  weightmode = adjacent;

  weights = g.get_weights();

  }


map::map(
#if defined(JAVA_OUTPUT_WINDOW)
administrator_basic * adbp,
#endif
const datamatrix & xo,const double & md, const metric & m)
  {

  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = adbp;
  #endif

  datamatrix x = xo;

  assert(m!=combnd);
  assert(x.cols() == 2);

  nocentroids = true;
  nopolygones = false;
  int i,j;

  //----------- sort x according to col 1, then according to col 2 -------------

  x.sort(0,x.rows()-1,0);

  unsigned beg = 0;
  unsigned end;
  for (i=1;i<x.rows();i++)
    {
    if ( x(i,0) != x(i-1,0) )
      {
      end = i-1;
      if (end-beg > 0)
        x.sort(beg,end,1);

      beg = i;
      }
    else if (i==x.rows()-1)
      {
      end = i;
      if (end-beg > 0)
        x.sort(beg,end,1);
      }

    }

  //--------- end: sort x according to col 1, then according to col 2 ----------


  //---------------- computing region names, region centroids ------------------

  ST::string name;

  region r;
//  name = ST::doubletostring(x(0,0),6) + "_" + ST::doubletostring(x(0,1),6);
  name = "1";
  r = region(name);
  r.set_center(x(0,0),x(0,1));
  regions.push_back(r);

  unsigned n=2;
  for(i=1;i<x.rows();i++)
    {
    if ( ( x(i,0) != x(i-1,0) ) || ( x(i,1) != x(i-1,1) ) )
      {
//      name = ST::doubletostring(x(i,0),6) + "_" + ST::doubletostring(x(i,1),6);
      name = ST::inttostring(n);
      n++;
      r = region(name);
      r.set_center(x(i,0),x(i,1));
      regions.push_back(r);
      }

    }

  nrregions = regions.size();

  //--------------- end: computing region names, region centroids --------------


  //---------------------- computing neighbors, bandsize -----------------------

  neighbors = vector< vector<unsigned> >(nrregions);


  double distance;
  double h1,h2;
  double xi,yi;
  double xcomp=0.0,ycomp=0.0;

  mindistance=NA;
  maxdistance=0;
  bandsize = 0;

  for (i=0;i<nrregions;i++)
    {

    xi = regions[i].get_xcenter();
    yi = regions[i].get_ycenter();

    j = i-1;
    if (j >= 0)
      {
      xcomp = regions[j].get_xcenter();
      ycomp = regions[j].get_ycenter();
      }
    while ((j >= 0) && (xi-xcomp <= md))
      {
      h1 = xi-xcomp;
      h2 = yi-ycomp;
      distance = sqrt(h1*h1+h2*h2);
      if (distance <= md )
        {
        if (distance < mindistance)
          mindistance = distance;
        if (distance > maxdistance)
          maxdistance = distance;
        neighbors[i].push_back(j);
        if (abs(i-j) > bandsize)
          {
          bandsize = abs(i-j);
          }

        }
      j--;
      if (j >= 0)
        {
        xcomp = regions[j].get_xcenter();
        ycomp = regions[j].get_ycenter();
        }
      }

    j = i+1;
    if (j < nrregions)
      {
      xcomp = regions[j].get_xcenter();
      ycomp = regions[j].get_ycenter();
      }
    while ((j < nrregions) && (xcomp-xi <= md))
      {
      h1 = xi-xcomp;
      h2 = yi-ycomp;
      distance = sqrt(h1*h1+h2*h2);
      if (distance <= md )
        {
        if (distance < mindistance)
          mindistance = distance;
        if (distance > maxdistance)
          maxdistance = distance;
        neighbors[i].push_back(j);
        if (abs(i-j) > bandsize)
          {
          bandsize = abs(i-j);
          }
        }
      j++;
      if (j < nrregions)
        {
        xcomp = regions[j].get_xcenter();
        ycomp = regions[j].get_ycenter();
        }
      }

    } // end: for (i=0;i<nrregions;i++)

  compute_minmaxn();

  //------------------------ end: computing neighbors --------------------------

  //--------------------------- computing weights ------------------------------

  weights = vector< vector<double> > (nrregions);

  if (m== adjacent)
    {
    for (i=0;i<nrregions;i++)
      {
      weights[i] = vector<double>(neighbors[i].size());
      for(j=0;j<neighbors[i].size();j++)
        weights[i][j] = 1.0;
      }
    }
  else // centroid
    {

    double dist;
    double sum=0;
    double nrneighbors=0;
    for (i=0;i<nrregions;i++)
      {
      for(j=0;j<neighbors[i].size();j++)
        {
        nrneighbors++;
        dist = regions[i].distance(regions[neighbors[i][j]], centroid);
        sum+= exp(-dist);
        }
      }

    double normconst = nrneighbors/sum;
    for (i=0;i<nrregions;i++)
      {
      weights[i] = vector<double>(neighbors[i].size());
      for(j=0;j<neighbors[i].size();j++)
        {
        dist = regions[i].distance(regions[neighbors[i][j]], centroid);
        weights[i][j] = normconst*exp(-dist);
        }
      }

    }  // end: centroid

  //------------------------ end: computing weights ----------------------------

  //-------------------- computing polygones of regions ------------------------

  #if defined(MICROSOFT_VISUAL)
    {
    minX = DBL_MAX;
    maxX = -DBL_MAX;
    minY = DBL_MAX;
    maxY = -DBL_MAX;
    }
  #else
    {
    minX = MAXDOUBLE;
    maxX = -MAXDOUBLE;
    minY = MAXDOUBLE;
    maxY = -MAXDOUBLE;
    }
  #endif

  polygone help;

  double s = mindistance/2;
  double xc,yc;
  for (i=0;i<nrregions;i++)
    {
    help = polygone();
    xc = regions[i].get_xcenter();
    yc = regions[i].get_ycenter();
    help.add_line(line(xc-s,yc-s,xc-s,yc+s));
    help.add_line(line(xc-s,yc+s,xc+s,yc+s));
    help.add_line(line(xc+s,yc+s,xc+s,yc-s));
    help.add_line(line(xc+s,yc-s,xc-s,yc-s));
    if (xc-s < minX)
      minX = xc-s;
    if (xc+s > maxX)
      maxX = xc+s;
    if (yc-s < minY)
      minY = yc-s;
    if (yc+s > maxY)
      maxY = yc+s;

    regions[i].add_polygone(help);
    }

  //------------------ end: computing polygones of regions ---------------------

  }


map::map(const map & m)
  {
  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = m.adminb_p,
  #endif
  weightmode = m.weightmode;
  errormessages = m.errormessages;
  regions = m.regions;
  nopolygones = m.nopolygones;
  nocentroids = m.nocentroids;
  nrregions = m.nrregions;
  minX = m.minX;
  maxX = m.maxX;
  minY = m.minY;
  maxY = m.maxY;
  neighbors = m.neighbors;
  minn = m.minn;
  maxn = m.maxn;
  bandsize = m.bandsize;
  weights = m.weights;
  mindistance = m.mindistance;
  maxdistance = m.maxdistance;
  }

const map & map::operator=(const map & m)
  {
  if (this == &m)
	 return *this;
  #if defined(JAVA_OUTPUT_WINDOW)
  adminb_p = m.adminb_p,
  #endif
  weightmode = m.weightmode;
  errormessages = m.errormessages;
  regions = m.regions;
  nopolygones = m.nopolygones;
  nocentroids = m.nocentroids;
  nrregions = m.nrregions;
  minX = m.minX;
  maxX = m.maxX;
  minY = m.minY;
  maxY = m.maxY;
  neighbors = m.neighbors;
  minn = m.minn;
  maxn = m.maxn;
  bandsize = m.bandsize;
  weights = m.weights;
  mindistance = m.mindistance;
  maxdistance = m.maxdistance;
  return *this;
  }



ostream & operator<<(ostream & c, const map & m)
  {
  if (m.errormessages.empty())
    {
    c << "Number of regions:                       " << m.nrregions << endl;
    if (!m.nopolygones)
      {
      c << "Minimum X-coordinate:                    " << m.minX << endl;
      c << "Maximum X-coordinate:                    " << m.maxX << endl;
      c << "Minimum Y-coordinate:                    " << m.minY << endl;
      c << "Maximum Y-coordinate:                    " << m.maxY << endl;
      }
    c << "Neighborhood definition:                 " << "first order" << endl;
    c << "Minimum number of neighbors:             " << m.minn << endl;
    c << "Maximum number of neighbors:             " << m.maxn << endl;
    if (m.weightmode == adjacent)
      c << "Weight definition:                       " << "adjacent" << endl;
    if (m.weightmode == centroid)
      c << "Weight definition:                       " << "centroid" << endl;
    if (m.weightmode == combnd)
      c << "Weight definition:                       " << "combnd" << endl;
    c << "Bandsize:                                " << m.bandsize << endl;
    return c;
    }
  else
    {
    unsigned i;
    for(i=0;i<m.errormessages.size();i++)
      c << m.errormessages[i];
    return c;
    }

  }


void map::outmap(const ST::string & path) const
  {
  assert(!nopolygones);
  ofstream out(path.strtochar());
  unsigned i,j,k;
  line l;
  line lalt;
  for(i=0;i<nrregions;i++)
    {
    for(j=0;j<regions[i].get_nrpoly();j++)
      {

      out << "\"" << regions[i].get_name() << "\","
          << ((regions[i].get_polygone(j)).get_nrlines()+1) << endl;
      if (regions[i].get_isin().length() > 0)
        out << "is.in,\"" << regions[i].get_isin() << "\"\n";
      for (k=0;k<(regions[i].get_polygone(j)).get_nrlines();k++)
        {
        if (k > 0)
          lalt = l;
        l = (regions[i].get_polygone(j)).get_line(k);
        if (k == 0)
          {
          out << l.x1 << "," << l.y1 << endl;
          out << l.x2 << "," << l.y2 << endl;
          }
        else
          {
          if ( ((lalt.x1 == l.x1) && (lalt.y1 == l.y1)) ||
              ((lalt.x2 == l.x1) && (lalt.y2 == l.y1)) )
            out << l.x2 << "," << l.y2 << endl;
          else
            out << l.x1 << "," << l.y1 << endl;
          }

        }  // end: for k=0;...

      }

    }


  }


void map::outindices (const ST::string & path) const
  {

  ofstream out(path.strtochar());
  unsigned i,j;

  for (i=0;i<nrregions;i++)
    for (j=0;j<neighbors[i].size();j++)
      out << i << "  " << neighbors[i][j] << endl;

  }



void map::outcentroids(const ST::string & path)
  {

  ofstream out(path.strtochar());
  unsigned i;

  out << nrregions << endl;
  for (i=0;i<nrregions;i++)
    {
    out << regions[i].get_name() << "  ";
    out << regions[i].get_xcenter() << "  ";
    out << regions[i].get_ycenter() << endl;
    }


  }


void map::identify_regions(void)
  {
  unsigned j,i;
  statmatrix<ST::string> names(regions.size(),1);
  for(j=0;j<regions.size();j++)
    names(j,0) = regions[j].get_name();

  statmatrix<int> index(regions.size(),1);
  index.indexinit();
  names.indexsort(index,0,names.rows()-1,0,0);

  statmatrix<unsigned> mask(regions.size(),1,0);

  j=0;
  for(i=1;i<regions.size();i++)
    {
    if (regions[index(i,0)].get_name() == regions[index(j,0)].get_name())
      {
      regions[index(j,0)].add_polygone(regions[index(i,0)].get_polygone(0));
      mask(index(i,0),0) = 1;
      }
    else
      j=i;

    }


  vector<region> newregions;
  newregions.reserve(regions.size());
  for(i=0;i<regions.size();i++)
    if (mask(i,0) == 0)
      newregions.push_back(regions[i]);

  regions = newregions;

/*
  unsigned i=1;
  while (i<regions.size())
    {
    if (regions[i].get_name() == regions[i-1].get_name())
      {
      regions[i-1].add_polygone(regions[i].get_polygone(0));
      regions.erase(regions.begin()+i,regions.begin()+i+1);
      }
    else
      i++;
    } // end: while (i<regions.size())
*/


  nrregions = regions.size();
  }


void map::outneighbors(const ST::string & path,const bool & names) const
  {
  ofstream out(path.strtochar());
  assert(!out.fail());
  unsigned i,j;
  for(i=0;i<nrregions;i++)
    {
    out << regions[i].get_name() << " ";
    for (j=0;j<neighbors[i].size();j++)
      {
      if (names == true)
        out << regions[neighbors[i][j]].get_name() << " ";
      else
        out << neighbors[i][j] << " ";
      }
    out << endl;
    }

  }


matrix map::getneighbors(void) const
  {

  matrix help(nrregions,maxn+1,0);

  unsigned i,j;
  for (i=0;i<nrregions;i++)
    {
    help(i,0) = neighbors[i].size();
    for (j=0;j<neighbors[i].size();j++)
      help(i,j+1) = neighbors[i][j]+1;
    }

  return help;

  }


void map::outneighbors2(const ST::string & path) const
  {

  ofstream out(path.strtochar());
  assert(!out.fail());
  unsigned i,j;
  out << nrregions << endl;
  for(i=0;i<nrregions;i++)
    {
    out << regions[i].get_name() << endl;
    out << neighbors[i].size() << endl;
    for (j=0;j<neighbors[i].size();j++)
      out << neighbors[i][j] << " ";

    out << endl;
    }

  }


void map::outgraph(const ST::string & path,const bool & we) const
  {

  ofstream out(path.strtochar());
  assert(!out.fail());
  unsigned i,j;

  out << nrregions;
  out << endl;

  for(i=0;i<nrregions;i++)
    {
    out << regions[i].get_name() << endl;
    out << neighbors[i].size() << endl;
    for (j=0;j<neighbors[i].size();j++)
      {
      out << neighbors[i][j] << " ";
      }


    if (we==true)
      {

      for (j=0;j<weights[i].size();j++)
        {
        out << weights[i][j] << " ";
        }

      }


    out << endl;
    }

  }


graph map::get_graph(void) const
  {
  vector<ST::string> regnames(regions.size());
  unsigned i;
  for(i=0;i<regions.size();i++)
    regnames[i] = regions[i].get_name();

  return graph(regnames,neighbors,weights);
  }


void map::infile(const ST::string & path)
  {
  nrregions = 0;

  bool stop = false;
  ifstream fin(path.strtochar());
  assert(!fin.fail());

  ST::string li;               // stores lines read successively from file path
  ST::string help;
  long n=0.0;                      // number of points that build a polygone
  double n1=0.0,n2=0.0;                // x,y coordinates of a point
  double begx=0.0,begy=0.0;
  int nrlines=0;               // number of lines read (so far) from file 'path'
  int linenr=0;
  polygone helppoly;
  vector<ST::string> token;

  line lhelp;
  bool first=false;

#if defined(MICROSOFT_VISUAL)
  {
  minX = DBL_MAX;
  maxX = -DBL_MAX;
  minY = DBL_MAX;
  maxY = -DBL_MAX;
  }
#else
  {
  minX = MAXDOUBLE;
  maxX = -MAXDOUBLE;
  minY = MAXDOUBLE;
  maxY = -MAXDOUBLE;
  }
#endif

  while ( (!fin.eof()) && (stop == false) )
	 {
	 ST::getline(fin,li);
	 nrlines++;
	 li = li.eatallwhitespace();
     li = li.eatallcarriagereturns();
	 if (li.length() > 0)
		{
		token = li.strtoken(",");
		if (token.size() == 3)
		  {
		  if (token[0][0] == '"')                        // new region
			 {
			 if (token[2].strtolong(n) == 1)
				{
				stop = true;
				errormessages.push_back(
                "ERROR: " + token[2] +
                " cannot be read as a number in line " +
                (ST::inttostring(nrlines)) + "\n");
				}

			 if (token[0][token[0].length()-1] != '"')
				{
				stop = true;
				errormessages.push_back(
                "ERROR: missing \" in line " +
                (ST::inttostring(nrlines)) + "\n");
				}

			 if (stop == false)
				{
                nrregions++;
				help = token[0].deleteallsigns('"');
				regions.push_back(region(help));
                helppoly = polygone();
				linenr=0;
				first=true;
				}

			 } // end: if (token[0][0] == '"')
          else if (token[0] == "is.in")
            {
            ST::string h = token[2].deleteallsigns('"');
            regions[regions.size()-1].set_isin(h);
            }
		  else
			 {

			 if (token[0].strtodouble(n1) == 1)
				{
				stop = true;
				errormessages.push_back(
                "ERROR: " + token[0] +
                " cannot be read as a number in line " +
                (ST::inttostring(nrlines)) + "\n");
				}

			 if (token[2].strtodouble(n2) == 1)
				{
				stop = true;
				errormessages.push_back(
                "ERROR: " + token[2] +
                " cannot be read as a number in line " +
                (ST::inttostring(nrlines)) + "\n");
				}

			 if (linenr >= n-1)
				{
				stop = true;
				errormessages.push_back("ERROR: too many lines for region "
						 + regions[nrregions-1].get_name() + "\n");
				}

			 if (stop == false)
				{

                if (n1 < minX)
                  minX = n1;
                if (n1 > maxX)
                  maxX = n1;
                if (n2 < minY)
                  minY = n2;
                if (n2 > maxY)
                  maxY = n2;

				if (first == true)
				  {
				  begx = n1;
				  begy = n2;
				  first=false;
				  }
				else
				  {
				  lhelp = line(begx,begy,n1,n2);
				  begx = n1;
				  begy = n2;
                  helppoly.add_line(lhelp);
				  linenr++;
                  if (linenr == n-1)
                    regions[regions.size()-1].add_polygone(helppoly);
				  }


				} // end: if (stop == false)

			 }

		  } // end: if (token.size() == 3;
		else // token.size() != 3
		  {
		  stop = true;
		  errormessages.push_back("ERROR: " + li + " invalid\n");
		  }
        }

     if (stop==false)
       {
       #if defined(BORLAND_OUTPUT_WINDOW)
       stop = hauptformular->breakcommand();
       #elif defined(JAVA_OUTPUT_WINDOW)
       stop = adminb_p->breakcommand();
       #endif
       if (stop)
         errormessages.push_back("ERROR: reading map info not completed due to user break\n");
       }

	 }  // end:  while (!fin.eof())


  if (stop == false)
    {
//    sortmap();
    identify_regions();
    }
  else
    {
    regions.erase(regions.begin(),regions.end());
    nrregions = 0;
    minn = 0;
    maxn = 0;

   #if defined(MICROSOFT_VISUAL)
	{
    minX = DBL_MAX;
    maxX = -DBL_MAX;
    minY = DBL_MAX;
    maxY = -DBL_MAX;
	}
  #else
	{
    minX = MAXDOUBLE;
    maxX = -MAXDOUBLE;
    minY = MAXDOUBLE;
    maxY = -MAXDOUBLE;
	}
  #endif
   }

  }


void map::reset(void)
  {
  regions.erase(regions.begin(),regions.end());
  nrregions = 0;
  minn = 0;
  maxn = 0;

#if defined(MICROSOFT_VISUAL)
  {
  minX = DBL_MAX;
  maxX = -DBL_MAX;
  minY = DBL_MAX;
  maxY = -DBL_MAX;
  }
#else
  {
  minX = MAXDOUBLE;
  maxX = -MAXDOUBLE;
  minY = MAXDOUBLE;
  maxY = -MAXDOUBLE;
  }
#endif


  bandsize = 0;

  neighbors.erase(neighbors.begin(),neighbors.end());
  weights.erase(weights.begin(),weights.end());

  }


void map::infile_neighbors(const ST::string & path)
  {
  ifstream fin(path.strtochar());
  assert(!fin.fail());
  neighbors = vector< vector<unsigned> >(nrregions,vector<unsigned>());
  unsigned i,j;
  ST::string li;
  vector<ST::string> token;
  long v;
  int test;
  for(i=0;i<nrregions;i++)
    {
    ST::getline(fin,li);
    li = li.eatallcarriagereturns();
    assert (!fin.fail());
    token = li.strtoken(" ");
    assert(token[0] == regions[i].get_name());
    for(j=1;j<token.size();j++)
      {
      test = token[j].strtolong(v);
      assert(test==0);
      neighbors[i].push_back(unsigned(v));
      }

    } // end:   for(i=0;i<nrregions;i++)

  compute_minmaxn();

  compute_weights(weightmode);

  } // end: function


void map::computeneighbors(void)
  {
  neighbors.erase(neighbors.begin(),neighbors.end());
  neighbors = vector< vector<unsigned> >(nrregions,vector<unsigned>());
  int i,j;
  bool stop=false;

  bandsize = 0;

  unsigned bshelp;

  for(i=0;i<nrregions;i++)
    {

    bshelp = 0;

    for(j=0;j<nrregions;j++)
      {
      if ( (i!=j) && (regions[i].compare(regions[j])) )
        {
        neighbors[i].push_back(j);
        if (abs(i-j) > bshelp)
          {
          bshelp = abs(i-j);
          }
        }

      #if defined(BORLAND_OUTPUT_WINDOW)
      stop = hauptformular->breakcommand();
      #elif defined(JAVA_OUTPUT_WINDOW)
      stop = adminb_p->breakcommand();
      #endif
      if (stop)
        errormessages.push_back("ERROR: reading map info not completed due to user break\n");
      if (stop)
        break;

      }

    if (bshelp > bandsize)
      bandsize = bshelp;

    if (stop)
      break;

    }

  if (!stop)
    {
    compute_minmaxn();
    compute_weights(weightmode);
    }
  else
    {
    reset();
    }

  }


void map::reorder(const ST::string & path)
  {

  vector<region> reghelp(nrregions);

  ifstream in(path.strtochar());

  vector<unsigned> perm(nrregions);
  vector<unsigned> invp(nrregions);

  unsigned i,j;

  for(i=0;i<nrregions;i++)
    {
    in >> perm[i];
    invp[perm[i]] = i;
    }

  for(i=0;i<nrregions;i++)
    reghelp[i] = regions[perm[i]];

  for(i=0;i<nrregions;i++)
    regions[i] = reghelp[i];

  for(i=0;i<nrregions;i++)
    {
    regions[i].x_center();
    regions[i].y_center();
    }


  vector< vector<unsigned> > neighborsnew;

  bandsize = 0;
  for (i=0;i<nrregions;i++)
    for(j=0;j<neighbors[perm[i]].size();j++)
      {
      neighborsnew[i][j] = invp[neighbors[perm[i]][j]];
      if (abs(static_cast<int>(i)-static_cast<int>(neighborsnew[i][j])) > bandsize)
        bandsize = abs(static_cast<int>(i)-static_cast<int>(neighborsnew[i][j]));
      }

  neighbors = neighborsnew;

  vector< vector<double> > weightsnew(nrregions);

  for (i=0;i<nrregions;i++)
    for(j=0;j<weights[perm[i]].size();j++)
      weightsnew[i][j] = invp[weights[perm[i]][j]];

   weights = weightsnew;

  }


void map::reorderopt(void)
  {

  errormessages.erase(errormessages.begin(),errormessages.end());

  unsigned i,j;

  vector<unsigned> perm(nrregions);
  vector<unsigned> invp(nrregions);

  graph g = get_graph();

  if (g.get_nrgraphs()== 1)
    {

    perm = g.CMopt();

    for(i=0;i<nrregions;i++)
      invp[perm[i]] = i;

    vector<region> reghelp(nrregions);

    for(i=0;i<nrregions;i++)
      reghelp[i] = regions[perm[i]];

    for(i=0;i<nrregions;i++)
      regions[i] = reghelp[i];

    for(i=0;i<nrregions;i++)
      {
      regions[i].x_center();
      regions[i].y_center();
      }

    vector< vector<unsigned> > neighborsnew(nrregions);

    for (i=0;i<nrregions;i++)
      {
      neighborsnew[i] = vector<unsigned>(neighbors[perm[i]].size());
      for(j=0;j<neighbors[perm[i]].size();j++)
        {
      neighborsnew[i][j] = invp[neighbors[perm[i]][j]];
        }
      }


    neighbors = neighborsnew;

    bandsize = 0;
    for (i=0;i<nrregions;i++)
      {
      for(j=0;j<neighbors[i].size();j++)
        if (abs(static_cast<int>(i)-static_cast<int>(neighbors[i][j])) > bandsize)
          bandsize = abs(static_cast<int>(i)-static_cast<int>(neighbors[i][j]));

      }


    vector< vector<double> > weightsnew(nrregions);

    for (i=0;i<nrregions;i++)
      {
      weightsnew[i] = vector<double>(weights[perm[i]].size());
      for(j=0;j<weights[perm[i]].size();j++)
        weightsnew[i][j] = weights[perm[i]][j];

      }

    weights = weightsnew;

//    out_weights("d:\\daten\\angela\\weights.raw");

    } // end: if (g.nr_graphs() == 1)
  else
    errormessages.push_back("ERROR: Reordering is not possible, map is disconnected.\n");


  }


bool map::isconnected(void)
  {

  graph g = get_graph();

  if (g.get_nrgraphs() == 1)
    return true;
  else
    return false;
  }



void map::compute_minmaxn(void)
  {
  maxn = 0;
  minn = nrregions-1;
  unsigned i;
  for(i=0;i<nrregions;i++)
    {
    if (neighbors[i].size() > maxn)
      maxn = neighbors[i].size();
    if (neighbors[i].size() < minn)
      minn = neighbors[i].size();
    }
  }


void map::compute_weights(const metric m)
  {
  weightmode = m;
  unsigned i,j;
  weights = vector< vector<double> >(nrregions,vector<double>());
  switch(m)
    {
    case adjacent:
      for(i=0;i<nrregions;i++)
        {
        weights[i] = vector<double>(neighbors[i].size());
        for(j=0;j<neighbors[i].size();j++)
          weights[i][j]=1.0;
        }
    break;

    case centroid:
      {
      double dist;
      double sum=0;
      double nrneighbors=0;
      for (i=0;i<nrregions;i++)
        {
        for(j=0;j<neighbors[i].size();j++)
          {
          nrneighbors++;
          dist = regions[i].distance(regions[neighbors[i][j]], centroid);
          sum+= exp(-dist);
          }
        }

      double normconst = nrneighbors/sum;
      for (i=0;i<nrregions;i++)
        {
        weights[i] = vector<double>(neighbors[i].size());
        for(j=0;j<neighbors[i].size();j++)
          {
          dist = regions[i].distance(regions[neighbors[i][j]], centroid);
          weights[i][j] = normconst*exp(-dist);
          }
        }

      }
    break;

    case combnd:
      {
      double dist;
      double sum=0;
      double nrneighbors=0;
      for (i=0;i<nrregions;i++)
        {
        for(j=0;j<neighbors[i].size();j++)
          {
          nrneighbors++;
          dist = regions[i].distance(regions[neighbors[i][j]], combnd);
          sum+= dist;
          }
        }

      double normconst = nrneighbors/sum;

      for(i=0;i<nrregions;i++)
        {
        weights[i] = vector<double>(neighbors[i].size());
        for(j=0;j<neighbors[i].size();j++)
          weights[i][j]= normconst*
                         regions[neighbors[i][j]].distance(regions[i], combnd);
        }
      }
    break;

    }
  }


void map::out_weights(const ST::string & path) const
  {
  ofstream out(path.strtochar());
  assert(!out.fail());

  unsigned i,j;

  for(i=0;i<nrregions;i++)
    {
    out << regions[i].get_name() << " ";
    for (j=0;j<weights[i].size();j++)
      {
      out << weights[i][j] << " ";
      }
    out << endl;
    }

  }


double map::get_weightssum(const unsigned & i) const
  {
  unsigned j;
  double sum=0;
  for (j=0;j<weights[i].size();j++)
    sum+= weights[i][j];
  return sum;
  }


void map::drop_region(const unsigned & nr)
  {
  assert(regions.size() > 0);
  assert(nr < regions.size());
  nrregions--;
  regions.erase(regions.begin()+nr,regions.begin()+nr+1);
  computeneighbors();
  }



void map::compute_reg(const datamatrix & d,vector<int> & posbeg,
                      vector<int> & posend,vector<ST::string> & effectvalues,
                      statmatrix<int> & index)
  {

  errormessages.erase(errormessages.begin(),errormessages.end());

  bool error = false;

  unsigned j,i,nr;

  // Initialization of the index, Indexsort of moddata
  index = statmatrix<int>(d.rows(),1);
  index.indexinit();
  d.indexsort(index,0,d.rows()-1,0,0);

  if (nrregions <= 3)
    {
    error=true;
    errormessages.push_back("ERROR: not enough regions to estimate spatial effect\n");
    }


  vector<int> posbeghelp;
  vector<int> posendhelp;
  datamatrix effvalueshelp;
  statmatrix<int> effindex;


  if (error==false)
    {

    posbeg = vector<int>(nrregions,-1);
    posend = vector<int>(nrregions,-1);

    posbeghelp = vector<int>(nrregions,-1);
    posendhelp = vector<int>(nrregions,-1);

    effectvalues = vector<ST::string>(nrregions);

    effvalueshelp = datamatrix(nrregions,1);

    effindex = statmatrix<int>(nrregions,1);

    for(j=0;j<effectvalues.size();j++)
      {
      effectvalues[j] = getname(j);
      effectvalues[j].strtodouble(effvalueshelp(j,0));
      }

    effindex.indexinit();
    effvalueshelp.indexsort(effindex,0,effvalueshelp.rows()-1,0,0);

    nr =  0;
    for (j=0;j<effvalueshelp.rows();j++)
      {
      if ( nr < d.rows() )
        {

        if ( d(index(nr,0),0) == effvalueshelp(effindex(j,0),0) )
          {
          posbeghelp[j] = nr;
          while  ( (nr < d.rows()) && ( d(index(nr,0),0)
                  == effvalueshelp(effindex(j,0),0)) )
            nr++;
          posendhelp[j] = nr-1;
          }

        }  // end: if ( nr <.d.rows() )


      } // end:  for (j=0;j<effvalueshelp.rows();j++)

    if ( nr < d.rows() )
      {

      error=true;
      errormessages.push_back("ERROR: region " +
                        ST::doubletostring(d(index(nr,0),0)) +
                         " in the dataset is not contained in the map object\n");
      }


    } // end: if (error == false)


  if (error == false)
    {
    // index, posbeg, posend ändern

    statmatrix<int> inveffindex(effindex.rows(),1);
    for(j=0;j<posbeg.size();j++)
      inveffindex(effindex(j,0),0) = j;

    // inveffindex(i,0) = rank of the i-th effectvalue

    // TEST
    /*
    ofstream out("c:\\bayesx\\testh\\results\\inveffindex.res");
    inveffindex.prettyPrint(out);

    ofstream out2("c:\\bayesx\\testh\\results\\effindex.res");
    effindex.prettyPrint(out2);

    ofstream out3("c:\\bayesx\\testh\\results\\effvalues.res");
    effvalueshelp.prettyPrint(out3);
    */
    // TEST

    statmatrix<int> indexneu(index.rows(),1);

    unsigned k=0;
    for(j=0;j<posbeg.size();j++)
      {
      if (posbeghelp[inveffindex(j,0)] != -1)
        {
        for(i=posbeghelp[inveffindex(j,0)];i<=posendhelp[inveffindex(j,0)];i++)
          {
          indexneu(k,0) = index(i,0);
          k++;
          }
        }

      }

    k=0;
    for(j=0;j<posbeg.size();j++)
      {
      if (posbeghelp[inveffindex(j,0)] != -1)
        {
        posbeg[j] = k;
        posend[j] = k+posendhelp[inveffindex(j,0)] -posbeghelp[inveffindex(j,0)];
        k = posend[j]+1;
        }

      }

    index = indexneu;

    }


    // TEST
    /*
    ofstream out("c:\\bayesx\\testh\\results\\index.res");
    index.prettyPrint(out);

    ofstream out2("c:\\bayesx\\testh\\results\\posend.res");
    for (j=0;j<posend.size();j++)
      out2 << posend[j] << endl;

    ofstream out2a("c:\\bayesx\\testh\\results\\posbeg.res");
    for (j=0;j<posbeg.size();j++)
      out2a << posbeg[j] << endl;

    ofstream out3("c:\\bayesx\\testh\\results\\effvalues.res");
    effvalueshelp.prettyPrint(out3);

    ofstream out4("c:\\bayesx\\testh\\results\\index_data.res");
    index.prettyPrint(out4);
    */

    // TEST



  }



} // end: namespace map


#if defined(BORLAND_OUTPUT_WINDOW)
//---------------------------------------------------------------------------
#pragma package(smart_init)
#endif














