// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
//  This test is for the BlockDavidson solver using the Thyra interface
//  The Thyra objects will be extracted from Epetra objects using the
//  Epetra-Thyra interface.
//  Therefore, this test should yield identical results compared against
//  the Epetra-only BlockDavidson solver test.
//
#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBlockDavidson.hpp"
#include "AnasaziEpetraAdapter.hpp"
#include "AnasaziBasicSort.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

#ifdef HAVE_EPETRA_THYRA
#include "AnasaziThyraAdapter.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#endif

#include "ModeLaplace1DQ1.h"

int main(int argc, char *argv[]) 
{
  int i;
  int info = 0;

#ifdef EPETRA_MPI

  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

  Epetra_SerialComm Comm;

#endif

  int MyPID = Comm.MyPID();

  bool testFailed = false;
  bool verbose = 0;
  std::string which("SM");
  if (argc>1) {
    if (argv[1][0]=='-' && argv[1][1]=='v') {
      verbose = true;
    }
    else {
      which = argv[1];
    }
  }
  if (argc>2) {
    if (argv[2][0]=='-' && argv[2][1]=='v') {
      verbose = true;
    }
    else {
      which = argv[2];
    }
  }

  if (verbose && MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }

#ifndef HAVE_EPETRA_THYRA
  if (verbose && MyPid == 0) {
      cout << "Please configure Anasazi with:" << endl;
      cout << "--enable-epetra-thyra" << endl;
      cout << "--enable-anasazy-thyra" << endl;
  }
  return 0;
#endif

  Anasazi::ReturnType returnCode = Anasazi::Ok;  

  typedef Thyra::MultiVectorBase<double> MV;
  typedef Thyra::LinearOpBase<double>    OP;
  typedef Anasazi::MultiVecTraits<double,MV>    MVT;
  typedef Anasazi::OperatorTraits<double,MV,OP> OPT;

  //  Problem information
  int space_dim = 1;
  std::vector<double> brick_dim( space_dim );
  brick_dim[0] = 1.0;
  std::vector<int> elements( space_dim );
  elements[0] = 100;

  // Eigensolver parameters
  int nev = 4;
  int blockSize = 5;
  int maxBlocks = 8;
  int maxIters = 500;
  double tol = 1.0e-6;

  // Create problem
  Teuchos::RefCountPtr<ModalProblem> testCase = 
    Teuchos::rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );

  // Get the stiffness and mass matrices
  Teuchos::RefCountPtr<Epetra_Operator> K = 
    Teuchos::rcp( const_cast<Epetra_Operator *>(testCase->getStiffness()), false );
  Teuchos::RefCountPtr<Epetra_Operator> M = 
    Teuchos::rcp( const_cast<Epetra_Operator *>(testCase->getMass()), false );

  // Get a pointer to the Epetra_Map
  Teuchos::RefCountPtr<const Epetra_Map> Map =  
    Teuchos::rcp( &K->OperatorDomainMap(), false );

  // create an epetra multivector
  Teuchos::RefCountPtr<Epetra_MultiVector> ivec = 
      Teuchos::rcp( new Epetra_MultiVector(K->OperatorDomainMap(), blockSize) );
  ivec->Random();

  // create a Thyra::VectorSpaceBase
  Teuchos::RefCountPtr<const Thyra::MPIVectorSpaceBase<double> > epetra_vs = 
    Thyra::create_MPIVectorSpaceBase(Map);

  // create a ScalarProdVectorSpaceBase
  Teuchos::RefCountPtr<const Thyra::ScalarProdVectorSpaceBase<double> > sp_domain = 
    Teuchos::rcp_dynamic_cast<const Thyra::ScalarProdVectorSpaceBase<double> >(epetra_vs,true);

  // create a MultiVectorBase (from the Epetra_MultiVector)
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<double> > thyra_ivec = 
    Thyra::create_MPIMultiVectorBase(Teuchos::rcp_implicit_cast<Epetra_MultiVector>(ivec), 
                                     epetra_vs,sp_domain);

  // Create Thyra LinearOpBase objects from the Epetra_Operator objects
  Teuchos::RefCountPtr<Thyra::LinearOpBase<double> > thyra_K = 
    Teuchos::rcp( new Thyra::EpetraLinearOp(K) );
  Teuchos::RefCountPtr<Thyra::LinearOpBase<double> > thyra_M = 
    Teuchos::rcp( new Thyra::EpetraLinearOp(M) );

  // Create parameter list to pass into solver
  Teuchos::ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Max Iters", maxIters );
  MyPL.set( "Max Blocks", maxBlocks );
  MyPL.set( "Tol", tol );

  // Create default output manager 
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );

  // Set verbosity level
  if (verbose) {
    MyOM->SetVerbosity( Anasazi::FinalSummary + Anasazi::TimingDetails 
                      + Anasazi::IterationDetails + Anasazi::OrthoDetails );
  }

  // Create the sort manager
  Teuchos::RefCountPtr<Anasazi::BasicSort<double, MV, OP> > MySM = 
     Teuchos::rcp( new Anasazi::BasicSort<double, MV, OP>(which) );

  // Create eigenproblem

  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > MyProblem =
    Teuchos::rcp( 
      new Anasazi::BasicEigenproblem<double, MV, OP>(thyra_K, thyra_M, thyra_ivec) 
    );
  //  MyProblem->SetPrec( Teuchos::rcp( const_cast<Epetra_Operator *>(opStiffness->getPreconditioner()), false ) );

  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->SetSymmetric(true);

  // Set the number of eigenvalues requested and the blocksize the solver should use
  MyProblem->SetNEV( nev );

  // Inform the eigenproblem that you are finishing passing it information
  info = MyProblem->SetProblem();
  if (info)
    cout << "Anasazi::BasicEigenproblem::SetProblem() returned with code : "<< info << endl;

  // Create the eigensolver 
  Anasazi::BlockDavidson<double, MV, OP> MySolver(MyProblem, MySM, MyOM, MyPL);



  // Solve the problem to the specified tolerances or length
  returnCode = MySolver.solve();
  if (returnCode != Anasazi::Ok) {
    testFailed = true;
  }

  if (!testFailed) {
    // Get the eigenvalues and eigenvectors from the eigenproblem
    Teuchos::RefCountPtr<std::vector<double> > evals = MyProblem->GetEvals();
    Teuchos::RefCountPtr<MV> evecs = MyProblem->GetEvecs();
    
    // eigenCheck() doesn't know Thyra::MultiVecBase objects
    //if (verbose) 
    //  info = testCase->eigenCheck( *evecs, &(*evals)[0], 0 );
  
    // Compute the direct residual
    Teuchos::RefCountPtr<MV> Kvec, Mvec; 
    int numVecs = MVT::GetNumberVecs(*evecs);
    std::vector<double> normV( numVecs );
    Teuchos::SerialDenseMatrix<int,double> T(numVecs,numVecs);
    Kvec = MVT::Clone( *evecs, numVecs );
    Mvec = MVT::Clone( *evecs, numVecs );
  
    // Put eigenvalues on the diagonal of T
    for (i=0; i<numVecs; i++) {
      T(i,i) = (*evals)[i];
    }
    // Compute K*evecs
    returnCode = OPT::Apply( *thyra_K, *evecs, *Kvec );
    assert( returnCode==Anasazi::Ok );
    // Compute M*evecs
    returnCode = OPT::Apply( *thyra_M, *evecs, *Mvec );
    assert( returnCode==Anasazi::Ok );
    // Compute residuals: K*evecs - M*evecs*T
    MVT::MvTimesMatAddMv( -1.0, *Mvec, T, 1.0, *Kvec );
    // Compute norms of residuals
    MVT::MvNorm( *Kvec, &normV );
    
    for (i=0; i<nev; i++ ) {
      if ( Teuchos::ScalarTraits<double>::magnitude(normV[i]/(*evals)[i]) > 5.0e-5) {
        testFailed = true;
      }
    }
  
  }

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  if (testFailed) {
    if (verbose && MyPID==0)
      cout << "End Result: TEST FAILED" << endl;	
    return -1;
  }
  //
  // Default return value
  //
  if (verbose && MyPID==0)
    cout << "End Result: TEST PASSED" << endl;
  return 0;

}	
