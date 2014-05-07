//@HEADER
// ***********************************************************************
//
//                     Aristos Optimization Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "GLdistApp_SchurOp.hpp"

namespace GLdistApp {

int SchurOp::Apply(const Epetra_MultiVector& x, Epetra_MultiVector& y) const
{
  Epetra_MultiVector statetmp1( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector statetmp2( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector nn( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector an( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector na( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector aa( (Epetra_BlockMap&)A_->DomainMap(), 1);
  Epetra_MultiVector ctrldomain( (Epetra_BlockMap&)B_->DomainMap(), 1);
  Epetra_MultiVector bb( (Epetra_BlockMap&)B_->RangeMap(), 1);

  N_->Multiply(true, x, statetmp1); // N'*x
  A_->Multiply(true, x, statetmp2); // A'*x

  N_->Multiply(false, statetmp1, nn); //N*N'*x
  A_->Multiply(false, statetmp1, an); //A*N'*x
  N_->Multiply(false, statetmp2, na); //N*A'*x
  A_->Multiply(false, statetmp2, aa); // A*A'*x
  
  B_->Multiply(true, x, ctrldomain); // B'*x
  B_->Multiply(false, ctrldomain, bb); // B*B'*x

  y.Update(1.0, nn, 1.0, an, 0.0);
  y.Update(1.0, na, 1.0, aa, 1.0);
  y.Update(1.0, bb, 1.0);

  return 0;
}

}  // End namespace GLdistApp.
