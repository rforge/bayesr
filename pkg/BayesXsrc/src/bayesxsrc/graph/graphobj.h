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


#if !defined (GRAPHOBJ_INCLUDED)

#define GRAPHOBJ_INCLUDED

#include"statobj.h"
#include"map.h"
#include"mapobject.h"


#if defined(JAVA_OUTPUT_WINDOW)
#include"adminparse_pointers.h"
#endif


class __EXPORT_TYPE graphobj : public statobject
  {

  private:

  #if defined(JAVA_OUTPUT_WINDOW)
  administrator_pointer * adminp_p;
  #endif


  MAP::map mapinfo;
  datamatrix D;
  vector<ST::string> varnames;

  vector<statobject*> * statobj;

  //----------------------- PRIVATE VARIABLES ----------------------------------

  typedef void (* runpointer)(graphobj &);

  runpointer functions[8];

  modelStandard m;

  expression e;

  use u;

  // for method 'drawmap'

  stroption mapname;
  stroption psname;
  stroption title;
  intoption nrcolors;
  doubleoption upperlimit;
  doubleoption lowerlimit;
  simpleoption color;
  simpleoption nolegend;
  simpleoption swapcolors;
  simpleoption drawnames;
  simpleoption hclcolors;
  simpleoption replace;
  simpleoption pcat;

  optionlist drawmapoptions;

  friend void drawmaprun(graphobj & o);

  // for method 'plotnonp'

  stroption xlab;
  stroption ylab;
  stroption connect;
  intoption height;
  intoption width;
  doubleoption xlimtop;
  doubleoption xlimbottom;
  doubleoption ylimtop;
  doubleoption ylimbottom;
  doubleoption xstep;
  doubleoption xstart;
  doubleoption ystep;
  doubleoption ystart;
  intoption year;
  intoption month;
  intoption linewidth;
  intoption fontsize;
  intoption pointsize;
  stroption linecolor;
  doubleoption titlescale;

  optionlist plotnonpoptions;

  friend void plotnonprun(graphobj & o);

  // for method 'plotsample'

  optionlist plotsampleoptions;

  friend void plotsamplerun(graphobj & o);

  // for method 'plotautocor'

  simpleoption mean;

  optionlist plotautocoroptions;

  friend void plotautocorrun(graphobj & o);

  // for method 'plotsurf'

  stroption zlab;
  doubleoption xrot;
  doubleoption yrot;
  doubleoption zrot;
  doubleoption zmin;
  doubleoption zmax;
  doubleoption zstart;
  doubleoption zstep;
  doubleoption zlimtop;
  doubleoption zlimbottom;
  intoption gridsize;

  optionlist plotsurfoptions;

  friend void plotsurfrun(graphobj & o);

  //------------------------ PRIVATE FUNCTIONS ---------------------------------

  void create(void);

  void changedescription(void);

  public:


  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR

  graphobj(void) : statobject()
	 {
	 type = "graph";
	 }

  // CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - name = n

  #if defined(JAVA_OUTPUT_WINDOW)
  graphobj(administrator_basic * adb, administrator_pointer * adp,
           const ST::string & n,ofstream * lo,istream * in,
           vector<statobject*> * st);
  #else
  graphobj(const ST::string & n,ofstream * lo,istream * in,
           vector<statobject*> * st);
  #endif

  // COPY CONSTRUCTOR

  graphobj(const graphobj & o);

  // OVERLOADED ASSIGNMENT OPERATOR

  const graphobj & operator=(const graphobj & o);

  // FUNCTION: parse
  // TASK: parses command c

  int parse(const ST::string & c);

  void describe(const optionlist & globaloptions = optionlist());


  };


void plotnonprun(graphobj & o);
void drawmaprun(graphobj & o);
void plotsamplerun(graphobj & o);
void plotsurfrun(graphobj & o);
void plotautocorrun(graphobj & o);

#endif
