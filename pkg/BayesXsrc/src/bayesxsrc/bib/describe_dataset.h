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
#ifndef describe_datasetH
#define describe_datasetH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <DBGrids.hpp>
#include <Grids.hpp>
#include <data.h>
//---------------------------------------------------------------------------
class Tdatasetform : public TForm
{
__published:	// Komponenten, die von der IDE verwaltet werden
    TStringGrid *StringGrid;
    TScrollBar *VertScrollBar;
    TScrollBar *HorzScrollBar;
    void __fastcall FormActivate(TObject *Sender);
    
    
    
    void __fastcall VertScrollBarChange(TObject *Sender);
    void __fastcall FormResize(TObject *Sender);
    
    void __fastcall FormClose(TObject *Sender, TCloseAction &Action);
    
        void __fastcall HorzScrollBarChange(TObject *Sender);
    
private:	// Benutzerdeklarationen
public:		// Benutzerdeklarationen
dataset * datap;
ST::string dataname;
    __fastcall Tdatasetform(TComponent* Owner);
void fillGrid(void);    
};
//---------------------------------------------------------------------------
extern PACKAGE Tdatasetform *datasetform;
//---------------------------------------------------------------------------
#endif
