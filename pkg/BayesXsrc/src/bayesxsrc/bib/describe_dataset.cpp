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

#include "describe_dataset.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
Tdatasetform *datasetform;
//---------------------------------------------------------------------------
__fastcall Tdatasetform::Tdatasetform(TComponent* Owner)
    : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void Tdatasetform::fillGrid()
{

  Caption = ("dataset " + dataname).strtochar();

  list<ST::string> liste = datap->getVarnames();
  list<ST::string>::iterator it;

  StringGrid->Top = 0;
  StringGrid->Left = 0;
  StringGrid->Height = ClientHeight - 17;
  StringGrid->Width = ClientWidth - 17;

  HorzScrollBar->Width = ClientWidth - 17;
  HorzScrollBar->Top = ClientHeight - 17;
  VertScrollBar->Left = ClientWidth - 17;
  VertScrollBar->Height = ClientHeight - 17;

  int datarows = datap->obs();
  int datacols = datap->getVarnames().size();
  StringGrid->RowCount = ClientHeight/(StringGrid->DefaultRowHeight+1)+1;
//  StringGrid->ColCount = ClientWidth/(StringGrid->DefaultColWidth+1)+1;
  StringGrid->ColCount = datacols+1+ClientWidth/(StringGrid->DefaultColWidth+1)+1;

  VertScrollBar->Min = 1;
  VertScrollBar->Max = datap->obs();
  VertScrollBar->LargeChange = ClientHeight/(StringGrid->DefaultRowHeight+1)-2;
  HorzScrollBar->Min = 1;
  HorzScrollBar->Max = datap->getVarnames().size();
  HorzScrollBar->LargeChange = ClientWidth/(StringGrid->DefaultColWidth+1)-1;

  int i, j, beginrow, endrow, begincol, endcol;

  beginrow = VertScrollBar->Position;
  endrow = VertScrollBar->Position + ClientHeight/(StringGrid->DefaultRowHeight+1);
  begincol = HorzScrollBar->Position;
//  endcol = HorzScrollBar->Position + ClientWidth/(StringGrid->DefaultColWidth+1);

  it = liste.begin();
  for(i=1;i<begincol;i++)
    it++;
//  for(i=begincol;i<endcol+1;i++,it++)
  for(i=begincol;i<StringGrid->ColCount+1;i++,it++)
    if(i<=datacols)
      StringGrid->Cells[i-begincol+1][0] = (*it).strtochar();
    else
      StringGrid->Cells[i-begincol+1][0] = "";

  for(j=beginrow;j<endrow+1;j++)
    if(j<=datarows)
      StringGrid->Cells[0][j-beginrow+1] = j;
    else
      StringGrid->Cells[0][j-beginrow+1] = "";

//  for (i=begincol; i<endcol+1; i++)
  for(i=begincol;i<StringGrid->ColCount+1;i++,it++)
    {
    datap->set_iterator(i);
    for (j=beginrow; j<endrow+1; j++)
      {
      if(j <= datarows && i <= datacols)
        {
        if(datap->getvalue(j-1)==MAXDOUBLE)
          StringGrid->Cells[i-begincol+1][j-beginrow+1] = ".";
        else
          StringGrid->Cells[i-begincol+1][j-beginrow+1] =
            ST::doubletostring(datap->getvalue(j-1),8).strtochar();
        }
      else
        StringGrid->Cells[i-begincol+1][j-beginrow+1] = "";
      }
    }
}
//---------------------------------------------------------------------------
void __fastcall Tdatasetform::FormActivate(TObject *Sender)
{
  fillGrid();
}
//---------------------------------------------------------------------------

void __fastcall Tdatasetform::VertScrollBarChange(TObject *Sender)
{
  fillGrid();
}
//---------------------------------------------------------------------------

void __fastcall Tdatasetform::FormResize(TObject *Sender)
{
  fillGrid();
}
//---------------------------------------------------------------------------


void __fastcall Tdatasetform::FormClose(TObject *Sender,
      TCloseAction &Action)
{
  VertScrollBar->Position = VertScrollBar->Min;
  HorzScrollBar->Position = HorzScrollBar->Min;
}
//---------------------------------------------------------------------------


void __fastcall Tdatasetform::HorzScrollBarChange(TObject *Sender)
{
  fillGrid();
}
//---------------------------------------------------------------------------


