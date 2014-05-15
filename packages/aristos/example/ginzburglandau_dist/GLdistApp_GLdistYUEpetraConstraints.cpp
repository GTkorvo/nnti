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

#include "GLdistApp_GLdistYUEpetraConstraints.hpp"

namespace GLdistApp {

GLdistYUEpetraConstraints::GLdistYUEpetraConstraints(
    Teuchos::RefCountPtr<GLdistYUEpetraDataPool> dat )
  :dat_(dat)
{}

void GLdistYUEpetraConstraints::getValue( const Aristos::Vector &x,
    Aristos::Vector &c ) const
{
  // Dynamic cast back to Epetra_MultiVectors here.
  Teuchos::RefCountPtr<const Epetra_MultiVector> exy =
    (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(x))).getYVector();
  Teuchos::RefCountPtr<const Epetra_MultiVector> exu =
    (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(x))).getUVector();
  Teuchos::RefCountPtr<Epetra_MultiVector> ec =
    Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(c)).getYVector());

  Teuchos::RefCountPtr<Epetra_MultiVector> ay =
      Teuchos::rcp(new Epetra_MultiVector((dat_->getA())->RangeMap(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> bu =
      Teuchos::rcp(new Epetra_MultiVector((dat_->getB())->RangeMap(),1));
  
  (dat_->getA())->Multiply(false, *exy, *ay);
  (dat_->getB())->Multiply(false, *exu, *bu);

  ec->Update(1.0, *ay,  1.0, *(dat_->getNy()), 0.0);
  ec->Update(1.0, *bu, -1.0, *(dat_->getb()),  1.0);
}


void GLdistYUEpetraConstraints::applyJacobian( bool Trans, const Aristos::Vector &x,
    const Aristos::Vector &v, Aristos::Vector &Jv) const
{
  if (Trans) {
    // Dynamic cast back to YUEpetra vectors here.
    Teuchos::RefCountPtr<const Epetra_MultiVector> evy =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(v))).getYVector();
    Teuchos::RefCountPtr<Epetra_MultiVector> eJvy =
      Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(Jv)).getYVector());
    Teuchos::RefCountPtr<Epetra_MultiVector> eJvu =
      Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(Jv)).getUVector());
  
    Teuchos::RefCountPtr<Epetra_MultiVector> atvy =
        Teuchos::rcp(new Epetra_MultiVector((dat_->getA())->DomainMap(),1));
    Teuchos::RefCountPtr<Epetra_MultiVector> ntvy =
        Teuchos::rcp(new Epetra_MultiVector((dat_->getNpy())->DomainMap(),1));
    Teuchos::RefCountPtr<Epetra_MultiVector> btvy =
        Teuchos::rcp(new Epetra_MultiVector((dat_->getB())->DomainMap(),1));
      
    (dat_->getA())->Multiply(true, *evy, *atvy);
    (dat_->getNpy())->Multiply(true, *evy, *ntvy);
    (dat_->getB())->Multiply(true, *evy, *btvy);
    
    eJvy->Update(1.0, *atvy, 1.0, *ntvy, 0.0);
    eJvu->Update(1.0, *btvy, 0.0);
  }
  else {
    // Dynamic cast back to YUEpetra vectors here.
    Teuchos::RefCountPtr<const Epetra_MultiVector> evy =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(v))).getYVector();
    Teuchos::RefCountPtr<const Epetra_MultiVector> evu =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(v))).getUVector();
    Teuchos::RefCountPtr<Epetra_MultiVector> eJvy =
      Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(Jv)).getYVector());
  
    Teuchos::RefCountPtr<Epetra_MultiVector> avy =
        Teuchos::rcp(new Epetra_MultiVector((dat_->getA())->RangeMap(),1));
    Teuchos::RefCountPtr<Epetra_MultiVector> nvy =
        Teuchos::rcp(new Epetra_MultiVector((dat_->getNpy())->RangeMap(),1));
    Teuchos::RefCountPtr<Epetra_MultiVector> bvu =
        Teuchos::rcp(new Epetra_MultiVector((dat_->getB())->RangeMap(),1));
      
    (dat_->getA())->Multiply(false, *evy, *avy);
    (dat_->getNpy())->Multiply(false, *evy, *nvy);
    (dat_->getB())->Multiply(false, *evu, *bvu);
    
    eJvy->Update(1.0, *avy, 1.0, *nvy, 0.0);
    eJvy->Update(1.0, *bvu, 1.0);
  }
}


void GLdistYUEpetraConstraints::applyNullSp( bool Trans, const Aristos::Vector &x,
    const Aristos::Vector &v, Aristos::Vector &Wv, double * tol) const
{
  assert(Trans==false);

  // Dynamic cast back to YUEpetra vectors here.
  Teuchos::RefCountPtr<const Epetra_MultiVector> evy =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(v))).getYVector();
  Teuchos::RefCountPtr<const Epetra_MultiVector> evu =
      (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(v))).getUVector();
                      
  Teuchos::RefCountPtr<Epetra_MultiVector> eWvy =
      Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(Wv)).getYVector());
  Teuchos::RefCountPtr<Epetra_MultiVector> eWvu =
      Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(Wv)).getUVector());

  Teuchos::RefCountPtr<Epetra_MultiVector> ep1 =
      Teuchos::rcp(new Epetra_MultiVector((dat_->getA())->DomainMap(),1));
  Teuchos::RefCountPtr<Epetra_MultiVector> ep2 =
      Teuchos::rcp(new Epetra_MultiVector((dat_->getA())->RangeMap(),1));

  dat_->solveAugsysDyn( evy, evu, ep2, eWvy, eWvu, ep1, tol );
}

} // namespace GLdistApp
