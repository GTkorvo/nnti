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

#include "Aristos_EpetraVector.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Util.h"
#include "SimpleApp_SimpleEpetraObjective.hpp"
#include "SimpleApp_SimpleEpetraConstraints.hpp"
#include "SimpleApp_SimpleEpetraFeasStep.hpp"
#include "SimpleApp_SimpleEpetraLagMult.hpp"
#include "SimpleApp_SimpleEpetraHessVec.hpp"
#include "SimpleApp_SimpleEpetraDataPool.hpp"
#include "Aristos_SQPAlgo.hpp"

#include "Teuchos_GlobalMPISession.hpp"

int main(int argc, char *argv[])
{

  // This is a standard communicator declaration.
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
#ifdef HAVE_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // create maps for 3, 5, and 8 dimensions.
  Epetra_Map map3(3, 0, Comm);
  Epetra_Map map5(5, 0, Comm);
  Epetra_Map map8(8, 0, Comm);


  // create vectors
  
  Teuchos::RefCountPtr<Epetra_Vector> eresvec3 = Teuchos::rcp( new Epetra_Vector(map3) );
  eresvec3->PutScalar( 0.0 );
  Teuchos::RefCountPtr<Aristos::Vector> resvec3 = 
    Teuchos::rcp( new Aristos::EpetraVector( eresvec3 ) );

  Teuchos::RefCountPtr<Epetra_Vector> eresvec3s = Teuchos::rcp( new Epetra_Vector(map3) );
  eresvec3s->PutScalar( 0.0 );
  Teuchos::RefCountPtr<Aristos::Vector> resvec3s = 
    Teuchos::rcp( new Aristos::EpetraVector( eresvec3s ) );

  Teuchos::RefCountPtr<Epetra_Vector> eresvec5 = Teuchos::rcp( new Epetra_Vector(map5) );
  eresvec5->PutScalar( 0.0 );
  Teuchos::RefCountPtr<Aristos::Vector> resvec5 = 
    Teuchos::rcp( new Aristos::EpetraVector( eresvec5 ) );

  Teuchos::RefCountPtr<Epetra_Vector> eresvec5s = Teuchos::rcp( new Epetra_Vector(map5) );
  eresvec5s->PutScalar( 0.0 );
  Teuchos::RefCountPtr<Aristos::Vector> resvec5s = 
    Teuchos::rcp( new Aristos::EpetraVector( eresvec5s ) );


///////////////////////////////////////////////////////////////////

  // problem description
  std::cout << "\nThe below pertains to problem 18.2 p. 574 in \"Numerical Optimization\" by Nocedal and Wright:\n\n";
  std::cout << "    min   exp(x1*x2*x3*x4*x5) - 0.5*(x1^3+x2^3+1)^2\n";
  std::cout << "    s.t.  x1^2 + x2^2 + x3^2 + x4^2 + x5^2 - 10 = 0\n";
  std::cout << "                                x2*x3 - 5*x4*x5 = 0\n";
  std::cout << "                                x1^3 + x2^3 + 1 = 0\n\n";

///////////////////////////////////////////////////////////////////


  // construct Objective object
  SimpleApp::SimpleEpetraObjective obj;

///////////////////////////////////////////////////////////////////

  // For calls to the Constraints object, we'll first need to construct a 'workspace' structure,
  // which will contain the Jacobian matrix as well as the augmented system matrix.


  // In the following lines, we create a 5x3 workspace Jacobian and fill it with 1's.
  // Note that the graph (i.e. sparsity pattern) is unique and will never be changed.
  Epetra_CrsMatrix jac(Copy, map3, map5, 0);

  int indices[5] = {0, 1, 2, 3, 4};
  double r0[5];
    r0[0] = 2; r0[1] = 2; r0[2] = 2; r0[3] = 2;  r0[4] = 2;

  jac.InsertGlobalValues(0, 5, r0, indices);
  jac.InsertGlobalValues(1, 5, r0, indices);
  jac.InsertGlobalValues(2, 5, r0, indices);

  jac.PutScalar(1.0);

  jac.FillComplete(map5, map3);

  // In the following lines, we create a 8x8 workspace augmented matrix  and fill it with 1's.
  // Note that the graph (i.e. sparsity pattern) is unique and will never be changed.
  Epetra_CrsMatrix aug(Copy, map8, map8, 0);

  int indices0[4]   = {0, 5, 6, 7};
  int indices1[4]   = {1, 5, 6, 7};
  int indices2[4]   = {2, 5, 6, 7};
  int indices3[4]   = {3, 5, 6, 7};
  int indices4[4]   = {4, 5, 6, 7};
  double r0_4[4]    = {1.0, 1.0, 1.0, 1.0};
  int indices5_7[5] = {0, 1, 2, 3, 4};
  double r5_7[5]    = {1.0, 1.0, 1.0, 1.0, 1.0};

  aug.InsertGlobalValues(0, 4, r0_4, indices0);
  aug.InsertGlobalValues(1, 4, r0_4, indices1);
  aug.InsertGlobalValues(2, 4, r0_4, indices2);
  aug.InsertGlobalValues(3, 4, r0_4, indices3);
  aug.InsertGlobalValues(4, 4, r0_4, indices4);
  aug.InsertGlobalValues(5, 5, r5_7, indices5_7);
  aug.InsertGlobalValues(6, 5, r5_7, indices5_7);
  aug.InsertGlobalValues(7, 5, r5_7, indices5_7);

  aug.FillComplete(map8, map8);



  // Now we can build the DataPool object ...
  SimpleApp::SimpleEpetraDataPool dat(jac, aug);

  Teuchos::RefCountPtr<SimpleApp::SimpleEpetraDataPool> rcpdat = Teuchos::rcp( &dat, false );


///////////////////////////////////////////////////////////////////

  // construct Constraints object
  SimpleApp::SimpleEpetraConstraints constr(rcpdat);

//////////////////////////////////////////////////////////////////

  // construct feasibility step object
  SimpleApp::SimpleEpetraFeasStep feasstep(rcpdat);
  
//////////////////////////////////////////////////////////////////

  // construct Lagrange multiplier object
  SimpleApp::SimpleEpetraLagMult lagmult(rcpdat);
  
//////////////////////////////////////////////////////////////////

  // construct Hessian of Lagrangian times vector object
  SimpleApp::SimpleEpetraHessVec hessvec;
  
//////////////////////////////////////////////////////////////////

  // construct SQP algorithm
  Aristos::SQPAlgo sqp(dat, obj, constr, hessvec, lagmult, feasstep);


//////////////////////////////////////////////////////////////////

  // run derivative checks
  eresvec5->Random();
  eresvec5s->Random();
  eresvec3->Random();
  
  sqp.runDerivativeCheck(*resvec5, *resvec3, *resvec5s);  


//////////////////////////////////////////////////////////////////

  Teuchos::ParameterList parlist;
  parlist.set("Max Number of SQP Iterations", 50);

  (*eresvec5)[0]=3.0; (*eresvec5)[1]=2.0; (*eresvec5)[2]=2.0; (*eresvec5)[3]=1.0; (*eresvec5)[4]=1.0;
  eresvec3->PutScalar(0.0);
  eresvec3s->PutScalar(0.0);

  int iter=0;
  int iflag=0;

  // run SQP algorithm
  sqp.run(*resvec5, *resvec3, *resvec3s, iter, iflag, parlist);
  std::cout << endl << *eresvec5 << endl;

  std::cout << "End Result: TEST PASSED\n";

  return 0;

}
