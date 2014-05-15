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

#include "GLdistApp_GLdistYUEpetraHessVec.hpp"

namespace GLdistApp {

GLdistYUEpetraHessVec::GLdistYUEpetraHessVec( Teuchos::RefCountPtr<GLdistApp::GLdistYUEpetraDataPool> dat )
    : dat_(dat)
{}


void GLdistYUEpetraHessVec::getValue( const Aristos::Vector &x, const Aristos::Vector &l,
    const Aristos::Vector &s, Aristos::Vector &Hs) const
{

  // Dynamic cast back to Epetra_MultiVectors here.
  Teuchos::RefCountPtr<const Epetra_MultiVector> ey =
    (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(x))).getYVector();
  Teuchos::RefCountPtr<const Epetra_MultiVector> eu =
    (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(x))).getUVector();
  Teuchos::RefCountPtr<const Epetra_MultiVector> ely =
    (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(l))).getYVector();
  Teuchos::RefCountPtr<const Epetra_MultiVector> esy =
    (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(s))).getYVector();
  Teuchos::RefCountPtr<const Epetra_MultiVector> esu =
    (Teuchos::dyn_cast<Aristos::YUEpetraVector>(const_cast<Aristos::Vector&>(s))).getUVector();
  // Dynamic cast back to FE_Vectors here.
  Teuchos::RefCountPtr<Epetra_MultiVector> eHsy =
    Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(Hs)).getYVector());
  Teuchos::RefCountPtr<Epetra_MultiVector> eHsu =
    Teuchos::rcp_const_cast<Epetra_MultiVector>((Teuchos::dyn_cast<Aristos::YUEpetraVector>(Hs)).getUVector());

  double beta = dat_->getbeta();
  
  // Get volume and face mass matrices.
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> H = dat_->getH();
  Teuchos::RefCountPtr<Epetra_FECrsMatrix> R = dat_->getR();
  
  // Create standard and overlap map.
  Epetra_Map overlapmap(-1, (dat_->getpindx())->M(), (int*)(dat_->getpindx())->A(), 1, *(dat_->getCommPtr()));
  Epetra_Map standardmap(-1, (dat_->getipindx())->M(), (int*)(dat_->getipindx())->A(), 1, *(dat_->getCommPtr()));
  
  // Create overlapped versions of ey, ely, and esy, necessary for the call to nonlinhessvec.
  // Use the Epetra_Import utility.
  Epetra_Import Importer(overlapmap, standardmap);
  Teuchos::RefCountPtr<Epetra_MultiVector> impy  = Teuchos::rcp(new Epetra_MultiVector(overlapmap, 1));
  Teuchos::RefCountPtr<Epetra_MultiVector> impsy = Teuchos::rcp(new Epetra_MultiVector(overlapmap, 1));
  Teuchos::RefCountPtr<Epetra_MultiVector> imply = Teuchos::rcp(new Epetra_MultiVector(overlapmap, 1));
  impy->Import(*ey, Importer, Insert);
  impsy->Import(*esy, Importer, Insert);
  imply->Import(*ely, Importer, Insert);

  // Create temporary FEVector.
  Teuchos::RefCountPtr<Epetra_FEVector> tmpHsy = Teuchos::rcp(new Epetra_FEVector(standardmap));

  nonlinhessvec(*(dat_->getCommPtr()), *(dat_->getipindx()), *(dat_->getipcoords()), *(dat_->getpindx()),
                *(dat_->getpcoords()), *(dat_->gett()), impy, impsy, imply, tmpHsy);

  // Copy into the MultiVector.
  *(eHsy) = *(tmpHsy);

  Teuchos::RefCountPtr<Epetra_MultiVector> mass_sy  = Teuchos::rcp(new Epetra_MultiVector(standardmap, 1));
  H->Multiply(false, *esy, *mass_sy);
  eHsy->Update(1.0, *mass_sy, -1.0);
  
  R->Multiply(false, *esu, *eHsu);
  eHsu->Scale(beta);
}

}  // End namespace GLdistApp.
