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

#include "GLdistApp_GLdistYUEpetraFeasStep.hpp"

namespace GLdistApp {

GLdistYUEpetraFeasStep::GLdistYUEpetraFeasStep(
    Teuchos::RefCountPtr<GLdistYUEpetraDataPool> dat )
  :dat_(dat)
{}


void GLdistYUEpetraFeasStep::getValue( const Aristos::Vector &x,
    const Aristos::Vector &c, Aristos::Vector &fs, double &tol) const
{

  // Dynamic cast back to Sledge vectors here.
  Teuchos::RefCountPtr<const Epetra_MultiVector> ecy =
    (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(c))).getYVector();
  Teuchos::RefCountPtr<Epetra_MultiVector> efsy =
    Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(fs)).getYVector());
  Teuchos::RefCountPtr<Epetra_MultiVector> efsu =
    Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(fs)).getUVector());

  Teuchos::RefCountPtr<Epetra_MultiVector> ey = Teuchos::rcp(new Epetra_MultiVector(efsy->Map(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> eu = Teuchos::rcp(new Epetra_MultiVector(efsu->Map(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> eyres = Teuchos::rcp(new Epetra_MultiVector(efsy->Map(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> eures = Teuchos::rcp(new Epetra_MultiVector(efsu->Map(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> epres = Teuchos::rcp(new Epetra_MultiVector(efsy->Map(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> efsy1 = Teuchos::rcp(new Epetra_MultiVector(efsy->Map(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> efsy2 = Teuchos::rcp(new Epetra_MultiVector(efsy->Map(),1));
      

  dat_->solveAugsysDyn( ey, eu, ecy, eyres, eures, epres, &tol );


  (dat_->getA())->Multiply(true, *epres, *efsy1);
  (dat_->getNpy())->Multiply(true, *epres, *efsy2);
  efsy->Update(1.0, *efsy1, 1.0, *efsy2, 0.0);

  (dat_->getB())->Multiply(true, *epres, *efsu);
  
}

} // namespace GLsApp
