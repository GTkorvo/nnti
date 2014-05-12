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

#include "SimpleApp_SimpleEpetraFeasStep.hpp"
#include "Epetra_Map.h"

namespace SimpleApp {

SimpleEpetraFeasStep::SimpleEpetraFeasStep( Teuchos::RefCountPtr<SimpleEpetraDataPool> dat )
  : dat_(dat)
{}


void SimpleEpetraFeasStep::getValue( const Aristos::Vector &x, const Aristos::Vector &c,
                                 Aristos::Vector &fs, double &tol) const
{

  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_Vector> ex =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(x))).getVector();
  Teuchos::RefCountPtr<const Epetra_Vector> ec =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(c))).getVector();
  Teuchos::RefCountPtr<Epetra_Vector> efs =
      Teuchos::rcp_const_cast<Epetra_Vector>((Teuchos::dyn_cast<Aristos::EpetraVector>(fs)).getVector());

  // Now do arithmetic in Epetra-land.

  Epetra_Vector rhs((dat_->getAugmat())->RowMap());
  
  for (int i=0; i < 5; i++)
    rhs[i] = 0.0;
  for (int i=5; i < 8; i++)
    rhs[i] = (*ec)[i-5];

  Epetra_Vector result((dat_->getAugmat())->RowMap());

  dat_->solveAugsys(rhs, result);

  Epetra_Vector tmp(ec->Map(), false);

  for (int i=0; i < 3; i++)
    tmp[i] = result[i+5];

  Epetra_CrsMatrix* myjac = dat_->getJacobian();

  myjac->Multiply(true, tmp, *efs);

}

} // namespace SimpleApp
