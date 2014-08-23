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

#ifndef ADMINPARSEBASIC
#define ADMINPARSEBASIC

#if defined(BORLAND_OUTPUT_WINDOW)
#include <vcl.h>
#pragma hdrstop

#include "StatReview.h"
#include<StatwinFrame.h>
#include<statwin_haupt.h>
#elif defined(JAVA_OUTPUT_WINDOW)
#include<jni.h>
#endif


//------------------------------------------------------------------------------


class __EXPORT_TYPE administrator_basic
  {

  private:

  //------------------------ PRIVATE VARIABLES ---------------------------------


  public:

  bool pause;
  bool stop;
  bool processrunning;
  bool suppressoutput;


#if defined(JAVA_OUTPUT_WINDOW)
  JNIEnv* Java;
  jclass BayesX_cls;
  jobject BayesX_obj;
  jmethodID javaoutput;
#endif

  // CONSTRUCTOR

  administrator_basic(void);

  // DESTRUCTOR

  ~administrator_basic() {}

  bool breakcommand(void);

  bool get_pause(void)
    {
    return pause;
    }

  bool get_stop(void)
    {
    return stop;
    }

  bool get_processrunning(void)
    {
    return processrunning;
    }

  bool get_suppressoutput(void)
    {
    return suppressoutput;
    }

  void set_pause(bool b)
    {
    pause = b;
    }

  void set_stop(bool b)
    {
    stop = b;
    }

  void set_processrunning(bool b)
    {
    processrunning = b;
    }

  void set_suppressoutput(bool b)
    {
    suppressoutput = b;
    }



  };


#endif
