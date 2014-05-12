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

#include "GLdistApp_GLdistYUEpetraLagMult.hpp"

namespace GLdistApp {

GLdistYUEpetraLagMult::GLdistYUEpetraLagMult(
    Teuchos::RefCountPtr<GLdistYUEpetraDataPool> dat )
  :dat_(dat)
{}


void GLdistYUEpetraLagMult::getValue( const Aristos::Vector &x,
    const Aristos::Vector &g, Aristos::Vector &l, double &tol) const
{

  // Dynamic cast back to Epetra_MultiVectors here.
  Teuchos::RefCountPtr<const Epetra_MultiVector> egy =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(g))).getYVector();
  Teuchos::RefCountPtr<const Epetra_MultiVector> egu =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(g))).getUVector();
  Teuchos::RefCountPtr<Epetra_MultiVector> el =
      Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(l)).getYVector());
  
  Teuchos::RefCountPtr<Epetra_MultiVector> egp = Teuchos::rcp(new Epetra_MultiVector(egy->Map(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> ey = Teuchos::rcp(new Epetra_MultiVector(egy->Map(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> eu = Teuchos::rcp(new Epetra_MultiVector(egu->Map(),1));

  dat_->solveAugsysDyn( egy, egu, egp, ey, eu, el, &tol );

}

} // namespace GLdistApp
