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

#include "SimpleApp_SimpleEpetraConstraints.hpp"
#include "Aristos_EpetraVector.hpp"
#include "Epetra_Map.h"

namespace SimpleApp {

SimpleEpetraConstraints::SimpleEpetraConstraints( Teuchos::RefCountPtr<SimpleEpetraDataPool> dat )
  : dat_(dat)
{}

void SimpleEpetraConstraints::getValue( const Aristos::Vector &x, Aristos::Vector &c ) const
{
  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_Vector> ex =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(x))).getVector();
  Teuchos::RefCountPtr<Epetra_Vector> ec =
      Teuchos::rcp_const_cast<Epetra_Vector>((Teuchos::dyn_cast<Aristos::EpetraVector>(c)).getVector());

  // Now do arithmetic in Epetra-land.

  double x0 = (*ex)[0];
  double x1 = (*ex)[1];
  double x2 = (*ex)[2];
  double x3 = (*ex)[3];
  double x4 = (*ex)[4];

  (*ec)[0] = pow(x0,2) + pow(x1,2) + pow(x2,2) + pow(x3,2) + pow(x4,2) - 10;
  (*ec)[1] = x1*x2 - 5.0*x3*x4;
  (*ec)[2] = pow(x0,3) + pow(x1,3) + 1;

}


void SimpleEpetraConstraints::applyJacobian( bool Trans, const Aristos::Vector &x, const Aristos::Vector &v,
                                         Aristos::Vector &Jv) const
{

  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_Vector> ex =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(x))).getVector();
  Teuchos::RefCountPtr<const Epetra_Vector> ev =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(v))).getVector();
  Teuchos::RefCountPtr<Epetra_Vector> eJv =
      Teuchos::rcp_const_cast<Epetra_Vector>((Teuchos::dyn_cast<Aristos::EpetraVector>(Jv)).getVector());

  // Now do arithmetic in Epetra-land.
 
  Epetra_CrsMatrix* myjac = dat_->getJacobian();

  if (Trans)
    myjac->Multiply(true, *ev, *eJv);
  else
    myjac->Multiply(false, *ev, *eJv);

}


void SimpleEpetraConstraints::applyNullSp( bool Trans, const Aristos::Vector &x, const Aristos::Vector &v,
                                       Aristos::Vector &Wv, double * tol) const
{

  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_Vector> ex =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(x))).getVector();
  Teuchos::RefCountPtr<const Epetra_Vector> ev =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(v))).getVector();
  Teuchos::RefCountPtr<Epetra_Vector> eWv =
      Teuchos::rcp_const_cast<Epetra_Vector>((Teuchos::dyn_cast<Aristos::EpetraVector>(Wv)).getVector());

  // Now do arithmetic in Epetra-land.

  Epetra_Vector rhs((dat_->getAugmat())->RowMap());
  
  for (int i=0; i < 5; i++)
    rhs[i] = (*ev)[i];
  for (int i=5; i < 8; i++)
    rhs[i] = 0.0;

  Epetra_Vector res((dat_->getAugmat())->RowMap());

  dat_->solveAugsys(rhs, res);

  for (int i=0; i < 5; i++)
    (*eWv)[i] = res[i];

}

} // namespace SimpleApp
