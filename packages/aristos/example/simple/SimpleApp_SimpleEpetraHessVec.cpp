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

#include "SimpleApp_SimpleEpetraHessVec.hpp"
#include "Epetra_Map.h"
#include <math.h>

namespace SimpleApp {

SimpleEpetraHessVec::SimpleEpetraHessVec(){}


void SimpleEpetraHessVec::getValue( const Aristos::Vector &x, const Aristos::Vector &l, const Aristos::Vector &s,
                                Aristos::Vector &Hs) const
{

  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_Vector> ex =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(x))).getVector();
  Teuchos::RefCountPtr<const Epetra_Vector> el =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(l))).getVector();
  Teuchos::RefCountPtr<const Epetra_Vector> es =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(s))).getVector();
  Teuchos::RefCountPtr<Epetra_Vector> eHs =
      Teuchos::rcp_const_cast<Epetra_Vector>((Teuchos::dyn_cast<Aristos::EpetraVector>(Hs)).getVector());

  // Now do arithmetic in Epetra-land.

  double x1 = (*ex)[0];
  double x2 = (*ex)[1];
  double x3 = (*ex)[2];
  double x4 = (*ex)[3];
  double x5 = (*ex)[4];

  double l1 = (*el)[0];
  double l2 = (*el)[1];
  double l3 = (*el)[2];

  double s1 = (*es)[0];
  double s2 = (*es)[1];
  double s3 = (*es)[2];
  double s4 = (*es)[3];
  double s5 = (*es)[4];

  double expxi = exp(x1*x2*x3*x4*x5);

  (*eHs)[0] = ( pow(x2,2)*pow(x3,2)*pow(x4,2)*pow(x5,2)*expxi-9.0*pow(x1,4)-6.0*(pow(x1,3)+pow(x2,3)+1.0)*x1-2.0*l1-6.0*l3*x1 ) * s1  +
              ( x3*x4*x5*expxi+x2*pow(x3,2)*pow(x4,2)*pow(x5,2)*x1*expxi-9.0*pow(x2,2)*pow(x1,2) ) * s2 +
              ( x2*x4*x5*expxi+pow(x2,2)*x3*pow(x4,2)*pow(x5,2)*x1*expxi ) * s3 +
              ( x2*x3*x5*expxi+pow(x2,2)*pow(x3,2)*x4*pow(x5,2)*x1*expxi ) * s4 +
              ( x2*x3*x4*expxi+pow(x2,2)*pow(x3,2)*pow(x4,2)*x5*x1*expxi ) * s5;

  (*eHs)[1] = ( x3*x4*x5*expxi+x2*pow(x3,2)*pow(x4,2)*pow(x5,2)*x1*expxi-9.0*pow(x2,2)*pow(x1,2) ) * s1  +
              ( pow(x1,2)*pow(x3,2)*pow(x4,2)*pow(x5,2)*expxi-9.0*pow(x2,4)-6.0*(pow(x1,3)+pow(x2,3)+1.0)*x2-2.0*l1-6.0*l3*x2 ) * s2  +
              ( x1*x4*x5*expxi+pow(x1,2)*x3*pow(x4,2)*pow(x5,2)*x2*expxi-l2 ) * s3  +
              ( x1*x3*x5*expxi+pow(x1,2)*pow(x3,2)*x4*pow(x5,2)*x2*expxi ) * s4  +
              ( x1*x3*x4*expxi+pow(x1,2)*pow(x3,2)*pow(x4,2)*x5*x2*expxi ) * s5;

  (*eHs)[2] = ( x2*x4*x5*expxi+pow(x2,2)*x3*pow(x4,2)*pow(x5,2)*x1*expxi ) * s1  +
              ( x1*x4*x5*expxi+pow(x1,2)*x3*pow(x4,2)*pow(x5,2)*x2*expxi-l2 ) * s2  +
              ( pow(x1,2)*pow(x2,2)*pow(x4,2)*pow(x5,2)*expxi-2*l1 ) * s3  +
              ( x1*x2*x5*expxi+pow(x1,2)*pow(x2,2)*x4*pow(x5,2)*x3*expxi ) * s4  +
              ( x1*x2*x4*expxi+pow(x1,2)*pow(x2,2)*pow(x4,2)*x5*x3*expxi ) * s5;

  (*eHs)[3] = ( x2*x3*x5*expxi+pow(x2,2)*pow(x3,2)*x4*pow(x5,2)*x1*expxi ) * s1  +
              ( x1*x3*x5*expxi+pow(x1,2)*pow(x3,2)*x4*pow(x5,2)*x2*expxi ) * s2  +
              ( x1*x2*x5*expxi+pow(x1,2)*pow(x2,2)*x4*pow(x5,2)*x3*expxi ) * s3  +
              ( pow(x1,2)*pow(x2,2)*pow(x3,2)*pow(x5,2)*expxi-2*l1 ) * s4  +
              ( x1*x2*x3*expxi+pow(x1,2)*pow(x2,2)*pow(x3,2)*x5*x4*expxi+5*l2 ) * s5;

  (*eHs)[4] = ( x2*x3*x4*expxi+pow(x2,2)*pow(x3,2)*pow(x4,2)*x5*x1*expxi ) * s1  +
              ( x1*x3*x4*expxi+pow(x1,2)*pow(x3,2)*pow(x4,2)*x5*x2*expxi ) * s2  +
              ( x1*x2*x4*expxi+pow(x1,2)*pow(x2,2)*pow(x4,2)*x5*x3*expxi ) * s3  +
              ( x1*x2*x3*expxi+pow(x1,2)*pow(x2,2)*pow(x3,2)*x5*x4*expxi+5*l2 ) * s4  +
              ( pow(x1,2)*pow(x2,2)*pow(x3,2)*pow(x4,2)*expxi-2*l1 ) * s5;

}

} // namespace SimpleApp
