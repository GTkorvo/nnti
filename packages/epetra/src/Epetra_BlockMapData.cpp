
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER

#include "Epetra_BlockMapData.h"
#include "Epetra_HashTable.h"
#include "Epetra_Directory.h"
#include "Epetra_Comm.h"
//#include "Epetra_ConfigDefs.h" //DATA_DEBUG

//=============================================================================
Epetra_BlockMapData::Epetra_BlockMapData(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm & Comm) 
	: Comm_(Comm.Clone()),
		Directory_(0),
		NumGlobalElements_(NumGlobalElements),
		IndexBase_(IndexBase),
		ElementSize_(ElementSize),
		LastContiguousGIDLoc_(0),
		LIDHash_(0)
{
	//cout << "--BMD created, addr: " << this << endl; //DATA_DEBUG
}

//=============================================================================
Epetra_BlockMapData::~Epetra_BlockMapData() {
  if(LIDHash_ != 0) {
    delete LIDHash_;
    LIDHash_ = 0;
  }
  if(Directory_ != 0) {
    delete Directory_;
    Directory_ = 0;
  }
	if(Comm_ != 0) {
		delete Comm_;
		Comm_ = 0;
	}
	//cout << "--BMD destroyed, addr: " << this << endl; //DATA_DEBUG
}
