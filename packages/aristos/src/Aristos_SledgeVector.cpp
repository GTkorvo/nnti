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

#include "Aristos_SledgeVector.hpp"

namespace Aristos {

SledgeVector::SledgeVector( const Teuchos::RefCountPtr<USEME::umDVec> &sledge_vec )
  :sledge_vec_(sledge_vec)
    {}

// Overridden from Vector

double SledgeVector::innerProd( const Vector &x ) const
{
  SledgeVector &ex = Teuchos::dyn_cast<SledgeVector>(const_cast <Vector&>(x));
  return sledge_vec_->d_dot(*ex.sledge_vec_);
}

void SledgeVector::linComb( const double &alpha, const Vector &x, const double &beta )
{
  SledgeVector &ex = Teuchos::dyn_cast<SledgeVector>(const_cast <Vector&>(x));

  (*sledge_vec_) = alpha*(*ex.sledge_vec_) + beta*(*sledge_vec_);
}

void SledgeVector::Scale( const double &alpha )
{
  (*sledge_vec_) *= alpha;
}

void SledgeVector::Set( const double &alpha )
{
  (*sledge_vec_) = alpha;
}

void SledgeVector::Set( const double &alpha, const Vector &x )
{
  SledgeVector &ex = Teuchos::dyn_cast<SledgeVector>(const_cast <Vector&>(x));
  (*sledge_vec_) = alpha * (*ex.sledge_vec_);
}

Teuchos::RefCountPtr<Vector> SledgeVector::createVector() const
{
  return Teuchos::rcp( new SledgeVector( Teuchos::rcp(new USEME::umDVec(sledge_vec_->size(), "vec")) ));
}

Teuchos::RefCountPtr<const USEME::umDVec> SledgeVector::getVector() const
{
  return sledge_vec_;
}

} // namespace Aristos
