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

#include "Aristos_EpetraVector.hpp"

namespace Aristos {

EpetraVector::EpetraVector( const Teuchos::RefCountPtr<Epetra_Vector> &epetra_vec )
  :epetra_vec_(epetra_vec)
    {}

// Overridden from Vector

double EpetraVector::innerProd( const Vector &x ) const
{
  double dot[1];
  EpetraVector &ex = Teuchos::dyn_cast<EpetraVector>(const_cast <Vector&>(x));
  //EpetraVector &ex = Teuchos::dyn_cast<EpetraVector>(x);
  epetra_vec_->Dot( *ex.epetra_vec_, dot );
  return dot[0];
}

void EpetraVector::linComb( const double &alpha, const Vector &x, const double &beta )
{
  EpetraVector &ex = Teuchos::dyn_cast<EpetraVector>(const_cast <Vector&>(x));
  //EpetraVector &ex = Teuchos::dyn_cast<EpetraVector>(x);
  epetra_vec_->Update( alpha, *ex.epetra_vec_, beta );
}

void EpetraVector::Scale( const double &alpha )
{
  epetra_vec_->Scale( alpha );
}

void EpetraVector::Set( const double &alpha )
{
  epetra_vec_->PutScalar( alpha );
}

void EpetraVector::Set( const double &alpha, const Vector &x )
{
  EpetraVector &ex = Teuchos::dyn_cast<EpetraVector>(const_cast <Vector&>(x));
  epetra_vec_->Scale( alpha, *ex.epetra_vec_ );
}

Teuchos::RefCountPtr<Vector> EpetraVector::createVector() const
{
  return Teuchos::rcp( new EpetraVector( Teuchos::rcp(new Epetra_Vector(epetra_vec_->Map(),false)) ));
}

Teuchos::RefCountPtr<const Epetra_Vector> EpetraVector::getVector() const
{
  return epetra_vec_;
}


} // namespace Aristos
