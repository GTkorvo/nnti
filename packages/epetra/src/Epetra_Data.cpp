
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

#include "Epetra_Data.h"

//=============================================================================
Epetra_Data::Epetra_Data()
	: ReferenceCount_(1) {}

//=============================================================================
Epetra_Data::Epetra_Data(const Epetra_Data & Data)
	: ReferenceCount_(1) {}

//=============================================================================
Epetra_Data::~Epetra_Data() {}

//=============================================================================
void Epetra_Data::IncrementReferenceCount() {
	ReferenceCount_++;
}

//=============================================================================
void Epetra_Data::DecrementReferenceCount() {
	ReferenceCount_--;
}

//=============================================================================
int Epetra_Data::ReferenceCount() const {
		return(ReferenceCount_);
}
