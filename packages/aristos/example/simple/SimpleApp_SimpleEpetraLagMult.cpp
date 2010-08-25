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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "SimpleApp_SimpleEpetraLagMult.hpp"
#include "Epetra_Map.h"

namespace SimpleApp {

SimpleEpetraLagMult::SimpleEpetraLagMult( Teuchos::RefCountPtr<SimpleEpetraDataPool> dat )
  : dat_(dat)
{}


void SimpleEpetraLagMult::getValue( const Aristos::Vector &x, const Aristos::Vector &g,
                                 Aristos::Vector &l, double &tol) const
{

  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_Vector> ex =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(x))).getVector();
  Teuchos::RefCountPtr<const Epetra_Vector> eg =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(g))).getVector();
  Teuchos::RefCountPtr<Epetra_Vector> el =
      Teuchos::rcp_const_cast<Epetra_Vector>((Teuchos::dyn_cast<Aristos::EpetraVector>(l)).getVector());

  // Now do arithmetic in Epetra-land.

  Epetra_Vector rhs((dat_->getAugmat())->RowMap());
  
  for (int i=0; i < 5; i++)
    rhs[i] = (*eg)[i];
  for (int i=5; i < 8; i++)
    rhs[i] = 0.0;

  Epetra_Vector result((dat_->getAugmat())->RowMap());

  dat_->solveAugsys(rhs, result);

  for (int i=0; i < 3; i++)
    (*el)[i] = result[i+5];

}

} // namespace SimpleApp
