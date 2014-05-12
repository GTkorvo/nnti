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

#include "SimpleApp_SimpleEpetraObjective.hpp"
#include <math.h>

namespace SimpleApp {

SimpleEpetraObjective::SimpleEpetraObjective(){}

double SimpleEpetraObjective::getValue( const Aristos::Vector &x ) const
{
  double tmp = 0.0;
  Teuchos::RefCountPtr<const Epetra_Vector> ex = (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(x))).getVector();

  tmp = exp( (*ex)[0] * (*ex)[1] * (*ex)[2] * (*ex)[3] * (*ex)[4] ) -
        0.5 * pow( (pow((*ex)[0],3) + pow((*ex)[1],3) + 1), 2);

  return tmp;
}

void SimpleEpetraObjective::getGradient( const Aristos::Vector &x, Aristos::Vector &g ) const
{
  double sc(3.0);

  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_Vector> ex = (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(x))).getVector();
  Teuchos::RefCountPtr<Epetra_Vector> eg = Teuchos::rcp_const_cast<Epetra_Vector>((Teuchos::dyn_cast<Aristos::EpetraVector>(g)).getVector());

  // Now do arithmetic in Epetra-land.

  double x0 = (*ex)[0];
  double x1 = (*ex)[1];
  double x2 = (*ex)[2];
  double x3 = (*ex)[3];
  double x4 = (*ex)[4];

  (*eg)[0] = x1*x2*x3*x4 * exp( x0*x1*x2*x3*x4 ) -
             3*pow(x0,2) * (pow(x0,3) + pow(x1,3) + 1);
  (*eg)[1] = x0*x2*x3*x4 * exp( x0*x1*x2*x3*x4 ) -
             3*pow(x1,2) * (pow(x0,3) + pow(x1,3) + 1);
  (*eg)[2] = x0*x1*x3*x4 * exp( x0*x1*x2*x3*x4 );
  (*eg)[3] = x0*x1*x2*x4 * exp( x0*x1*x2*x3*x4 );
  (*eg)[4] = x0*x1*x2*x3 * exp( x0*x1*x2*x3*x4 );

}

} // namespace SimpleApp
