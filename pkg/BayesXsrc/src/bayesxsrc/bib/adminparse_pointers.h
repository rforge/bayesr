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

#ifndef ADMINPARSEPOINTER
#define ADMINPARSEPOINTER

#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include "StatReview.h"
#include<StatwinFrame.h>
#include<statwin_haupt.h>
#endif

#include<fstream>
#include<string.h>
#include"data.h"
#include"map.h"
//#include<dir.h>
//------------------------------------------------------------------------------


class __EXPORT_TYPE administrator_pointer
  {

  private:

  //------------------------ PRIVATE VARIABLES ---------------------------------

  dataset * datap;
  MAP::map * mapinfo;
  datamatrix * Dp;
  vector<ST::string> * varnamesp;


  public:


  // CONSTRUCTOR

  administrator_pointer(void);

  // DESTRUCTOR

  ~administrator_pointer() {}

  dataset * get_datap(void)
    {
    return datap;
    }

  void set_datap(dataset *d)
    {
    datap = d;
    }

  MAP::map * get_mapinfo(void)
    {
    return mapinfo;
    }

  void set_mapinfo(MAP::map *m)
    {
    mapinfo = m;
    }

  datamatrix * get_Dp(void)
    {
    return Dp;
    }

  void set_Dp(datamatrix *D)
    {
    Dp = D;
    }

  vector<ST::string> * get_varnamesp(void)
    {
    return varnamesp;
    }

  void set_varnamesp(vector<ST::string> *varn)
    {
    varnamesp = varn;
    }


  };


#endif
