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


//---------------------------------------------------------------------------
#include <vcl.h>
#pragma hdrstop

#include "describe_map.h"
#include "StatObjects.h"
#include "statwin_haupt.h"
#include "mapobject.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
Tmapform *mapform;
//---------------------------------------------------------------------------
__fastcall Tmapform::Tmapform(TComponent* Owner)
    : TForm(Owner)
{
}
//---------------------------------------------------------------------------

void Tmapform::drawmap(void)
{

  Caption = ("map " + mapname).strtochar();

  double maxX = mapinfo->get_maxX();
  double maxY = mapinfo->get_maxY();
  double minX = mapinfo->get_minX();
  double minY = mapinfo->get_minY();

  ClientHeight = ClientWidth*(maxY-minY)/(maxX-minX);

  for(int i=0;i<mapinfo->get_nrregions();i++)
    {
    for(int j=0;j<mapinfo->get_region(i).get_nrpoly();j++)
      {
      for(int k=0;k<mapinfo->get_region(i).get_polygone(j).get_nrlines();k++)
        {
        Canvas->MoveTo((mapinfo->get_region(i).get_polygone(j).get_line(k).x1-minX)
                        *0.9*ClientWidth/(maxX-minX)+0.05*ClientWidth,
                       -(mapinfo->get_region(i).get_polygone(j).get_line(k).y1-minY)
                        *0.9*ClientHeight/(maxY-minY)+0.95*ClientHeight);
        Canvas->LineTo((mapinfo->get_region(i).get_polygone(j).get_line(k).x2-minX)
                        *0.9*ClientWidth/(maxX-minX)+0.05*ClientWidth,
                       -(mapinfo->get_region(i).get_polygone(j).get_line(k).y2-minY)
                        *0.9*ClientHeight/(maxY-minY)+0.95*ClientHeight);
        }
      }
    }

}
//---------------------------------------------------------------------------


void __fastcall Tmapform::FormPaint(TObject *Sender)
{
  drawmap();
}
//---------------------------------------------------------------------------


void __fastcall Tmapform::FormResize(TObject *Sender)
{
  Canvas->FillRect(ClientRect);
  drawmap();
}
//---------------------------------------------------------------------------

