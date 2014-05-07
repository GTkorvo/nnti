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

#ifndef SIMPLEAPP_DATAPOOL_H
#define SIMPLEAPP_DATAPOOL_H

#include "Aristos_DataPool.hpp"
#include "Aristos_EpetraVector.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

namespace SimpleApp {

class SimpleEpetraDataPool : public Aristos::DataPool {
private:

  Epetra_CrsMatrix jacobian_;
  Epetra_CrsMatrix augmat_;
  Epetra_LinearProblem augsys_;
  
public:

  SimpleEpetraDataPool( Epetra_CrsMatrix &jacobian, Epetra_CrsMatrix &augmat );
  
  void computeAll( const Aristos::Vector &x );

  void solveAugsys( const Epetra_Vector &rhs, Epetra_Vector &result );

  Epetra_CrsMatrix* getJacobian();
  Epetra_CrsMatrix* getAugmat();

};

} // namespace SimpleApp

#endif
