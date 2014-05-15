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

#include "SimpleApp_SimpleEpetraDataPool.hpp"
#include "Amesos.h"

namespace SimpleApp {

SimpleEpetraDataPool::SimpleEpetraDataPool( Epetra_CrsMatrix &jacobian, Epetra_CrsMatrix &augmat )
  : jacobian_(jacobian), augmat_(augmat), augsys_()
{}

void SimpleEpetraDataPool::computeAll( const Aristos::Vector &x )
{
  // Dynamic cast back to Epetra vectors here.
  Teuchos::RefCountPtr<const Epetra_Vector> ex =
      (Teuchos::dyn_cast<Aristos::EpetraVector>(const_cast<Aristos::Vector&>(x))).getVector();

  // Now do arithmetic in Epetra-land.

  double x0 = (*ex)[0];
  double x1 = (*ex)[1];
  double x2 = (*ex)[2];
  double x3 = (*ex)[3];
  double x4 = (*ex)[4];

  int indices0_4[5] = {0, 1, 2, 3, 4};
  double r0[5];
    r0[0] = 2*x0; r0[1] = 2*x1; r0[2] = 2*x2; r0[3] = 2*x3;  r0[4] = 2*x4;
  double r1[5]; 
    r1[0] = 0;    r1[1] = x2;   r1[2] = x1;   r1[3] = -5*x4; r1[4] = -5*x3;
  double r2[5];
    r2[0] = 3*pow(x0,2); r2[1] = 3*pow(x1,2); r2[2] = 0; r2[3] = 0; r2[4] = 0;

  jacobian_.ReplaceGlobalValues(0, 5, r0, indices0_4);
  jacobian_.ReplaceGlobalValues(1, 5, r1, indices0_4);
  jacobian_.ReplaceGlobalValues(2, 5, r2, indices0_4);

  augmat_.ReplaceGlobalValues(5, 5, r0, indices0_4);
  augmat_.ReplaceGlobalValues(6, 5, r1, indices0_4);
  augmat_.ReplaceGlobalValues(7, 5, r2, indices0_4);

  int indices5_7[3] = {5, 6, 7};
  double ar0[3]; ar0[0] = r0[0]; ar0[1] = r1[0]; ar0[2] = r2[0];
  double ar1[3]; ar1[0] = r0[1]; ar1[1] = r1[1]; ar1[2] = r2[1];
  double ar2[3]; ar2[0] = r0[2]; ar2[1] = r1[2]; ar2[2] = r2[2];
  double ar3[3]; ar3[0] = r0[3]; ar3[1] = r1[3]; ar3[2] = r2[3];
  double ar4[3]; ar4[0] = r0[4]; ar4[1] = r1[4]; ar4[2] = r2[4];

  augmat_.ReplaceGlobalValues(0, 3, ar0, indices5_7);
  augmat_.ReplaceGlobalValues(1, 3, ar1, indices5_7);
  augmat_.ReplaceGlobalValues(2, 3, ar2, indices5_7);
  augmat_.ReplaceGlobalValues(3, 3, ar3, indices5_7);
  augmat_.ReplaceGlobalValues(4, 3, ar4, indices5_7);

  augsys_.SetOperator(&augmat_);

}


void SimpleEpetraDataPool::solveAugsys(const Epetra_Vector &rhs, Epetra_Vector &result)
{
  augsys_.SetRHS( const_cast<Epetra_Vector *>( &rhs ) );
  augsys_.SetLHS( &result );

  Amesos_BaseSolver* Solver;
  // Amesos_Factory is the function class used to create the solver.
  // This class contains no data.
  Amesos Amesos_Factory;

  // empty parameter list
  Teuchos::ParameterList List;
  
  // may also try: "Amesos_Umfpack", "Amesos_Lapack", ...
  std::string SolverType = "Amesos_Umfpack";
  
  Solver = Amesos_Factory.Create(SolverType, augsys_);
  // Amesos_Factory returns 0 is the selected solver is not
  // available
  if (!Solver) {
    SolverType = "Amesos_Klu";
    Solver = Amesos_Factory.Create(SolverType, augsys_);
  }
  assert (Solver);

  // start solving
  Solver->SymbolicFactorization();
  Solver->NumericFactorization();
  Solver->Solve();

  // delete Solver
  delete Solver;
    
}


Epetra_CrsMatrix* SimpleEpetraDataPool::getJacobian()
{
  return &jacobian_;
}


Epetra_CrsMatrix* SimpleEpetraDataPool::getAugmat()
{
  return &augmat_;
}



} //namespace SimpleApp
