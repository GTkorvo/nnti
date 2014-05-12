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

#include "Aristos_YUEpetraVector.hpp"

namespace Aristos {

YUEpetraVector::YUEpetraVector( const Teuchos::RefCountPtr<Epetra_MultiVector> &y_epetra_vec,
                                const Teuchos::RefCountPtr<Epetra_MultiVector> &u_epetra_vec )
  :y_epetra_vec_(y_epetra_vec), u_epetra_vec_(u_epetra_vec)
    {}

// Overridden from Vector

double YUEpetraVector::innerProd( const Vector &x ) const
{
  double ydot[1];
  double udot[1];
  YUEpetraVector &ex = Teuchos::dyn_cast<YUEpetraVector>(const_cast <Vector&>(x));
  y_epetra_vec_->Dot( *ex.y_epetra_vec_, ydot );
  if (u_epetra_vec_.get() == 0)
    udot[0] = 0.0;
  else
    u_epetra_vec_->Dot( *ex.u_epetra_vec_, udot );
  return (ydot[0] + udot[0]);
}

void YUEpetraVector::linComb( const double &alpha, const Vector &x, const double &beta )
{
  YUEpetraVector &ex = Teuchos::dyn_cast<YUEpetraVector>(const_cast <Vector&>(x));
  y_epetra_vec_->Update( alpha, *ex.y_epetra_vec_, beta );
  if (u_epetra_vec_.get() != 0)
    u_epetra_vec_->Update( alpha, *ex.u_epetra_vec_, beta );
}

void YUEpetraVector::Scale( const double &alpha )
{
  y_epetra_vec_->Scale( alpha );
  if (u_epetra_vec_.get() != 0)
    u_epetra_vec_->Scale( alpha );
}

void YUEpetraVector::Set( const double &alpha )
{
  y_epetra_vec_->PutScalar( alpha );
  if (u_epetra_vec_.get() != 0)
    u_epetra_vec_->PutScalar( alpha );
}

void YUEpetraVector::Set( const double &alpha, const Vector &x )
{
  YUEpetraVector &ex = Teuchos::dyn_cast<YUEpetraVector>(const_cast <Vector&>(x));
  y_epetra_vec_->Scale( alpha, *ex.y_epetra_vec_ );
  if (u_epetra_vec_.get() != 0)
    u_epetra_vec_->Scale( alpha, *ex.u_epetra_vec_ );
}

Teuchos::RefCountPtr<Vector> YUEpetraVector::createVector() const
{
  Teuchos::RefCountPtr<Epetra_MultiVector> yptr =
      Teuchos::rcp(new Epetra_MultiVector(y_epetra_vec_->Map(),1,false));
  Teuchos::RefCountPtr<Epetra_MultiVector> uptr = Teuchos::null;
  if (u_epetra_vec_.get() != 0)
    uptr = Teuchos::rcp(new Epetra_MultiVector(u_epetra_vec_->Map(),1,false));
  
  return Teuchos::rcp( new YUEpetraVector( yptr, uptr ));
}

Teuchos::RefCountPtr<const Epetra_MultiVector> YUEpetraVector::getYVector() const
{
  return y_epetra_vec_;
}

Teuchos::RefCountPtr<const Epetra_MultiVector> YUEpetraVector::getUVector() const
{
  return u_epetra_vec_;
}

} // namespace Aristos
