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
#include "AnasaziBlockDavidson.hpp"
#include "AnasaziBasicSort.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"

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
  MPI_Comm mpiComm(Comm.Comm());
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

  bool testFailed = false;
  bool verbose = true;
  std::string which("LM");
  int numel;
  if (argc > 1) {
    numel = atoi(argv[1]);  
  }
  else {
    numel = 100;
  }
  

  if (verbose && MyPID == 0) {
    cout << Anasazi::Anasazi_Version() << endl;
#ifdef EPETRA_MPI
    cout << "MPI" << endl;
#else
    cout << "SERIAL" << endl;
#endif
    cout << "BlockDavidson, Thyra-Epetra" << endl;
    cout << "Num elements: " << numel << endl << endl;
  }

#ifndef HAVE_EPETRA_THYRA
  if (verbose && MyPid == 0) {
      cout << "Please configure Anasazi with:" << endl;
      cout << "--enable-epetra-thyra" << endl;
      cout << "--enable-anasazi-thyra" << endl;
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
  elements[0] = numel;

  // Eigensolver parameters
  int nev = 4;
  int blockSize = 5;
  int maxBlocks = 5;
  int maxIters = 50;
  double tol = 1.0e-10;

  // Create problem
  Teuchos::RefCountPtr<ModalProblem> testCase = 
    Teuchos::rcp( new ModeLaplace1DQ1(Comm, brick_dim[0], elements[0]) );

  // Get the stiffness and mass matrices
  Teuchos::RefCountPtr<Epetra_Operator> e_K = 
    Teuchos::rcp( const_cast<Epetra_Operator *>(testCase->getStiffness()), false );
  Teuchos::RefCountPtr<Epetra_Operator> e_M = 
    Teuchos::rcp( const_cast<Epetra_Operator *>(testCase->getMass()), false );

  // Get a pointer to the Epetra_Map
  Teuchos::RefCountPtr<const Epetra_Map> Map =  
    Teuchos::rcp( &e_K->OperatorDomainMap(), false );


  /********************* Create wrapped Epetra multivector ********************/
  // create an epetra multivector
  Teuchos::RefCountPtr<Epetra_MultiVector> e_ivec = 
      Teuchos::rcp( new Epetra_MultiVector(e_K->OperatorDomainMap(), blockSize) );
  e_ivec->Random();

  // create a thyra vectorspacebase from the epetra map
  Teuchos::RefCountPtr<const Thyra::MPIVectorSpaceBase<double> > e_vs = 
    Thyra::create_MPIVectorSpaceBase(Map);

  // then, a ScalarProdVectorSpaceBase
  Teuchos::RefCountPtr<const Thyra::ScalarProdVectorSpaceBase<double> > sp_domain = 
    Teuchos::rcp_dynamic_cast<const Thyra::ScalarProdVectorSpaceBase<double> >(
      e_vs->smallVecSpcFcty()->createVecSpc(e_ivec->NumVectors())
    );

  // create a MultiVectorBase (from the Epetra_MultiVector)
  Teuchos::RefCountPtr<Thyra::MultiVectorBase<double> > te_ivec = 
    Thyra::create_MPIMultiVectorBase(Teuchos::rcp_implicit_cast<Epetra_MultiVector>(e_ivec), 
                                     e_vs,sp_domain);


  // Create Thyra LinearOpBase objects from the Epetra_Operator objects
  Teuchos::RefCountPtr<Thyra::LinearOpBase<double> > te_K = 
    Teuchos::rcp( new Thyra::EpetraLinearOp(e_K) );
  Teuchos::RefCountPtr<Thyra::LinearOpBase<double> > te_M = 
    Teuchos::rcp( new Thyra::EpetraLinearOp(e_M) );

  // Create parameter list to pass into solver
  Teuchos::ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Max Iters", maxIters );
  MyPL.set( "Max Blocks", maxBlocks );
  MyPL.set( "Tol", tol );
  MyPL.set( "Restart Timers", true);

  // Create default output manager 
  Teuchos::RefCountPtr<Anasazi::OutputManager<double> > MyOM = Teuchos::rcp( new Anasazi::OutputManager<double>( MyPID ) );

  // Set verbosity level
  if (verbose) {
    MyOM->SetVerbosity( Anasazi::FinalSummary + Anasazi::TimingDetails );
  }

  // Create the sort manager
  Teuchos::RefCountPtr<Anasazi::BasicSort<double, MV, OP> > MySM = 
     Teuchos::rcp( new Anasazi::BasicSort<double, MV, OP>(which) );

  Teuchos::RefCountPtr<Anasazi::BasicEigenproblem<double, MV, OP> > te_Problem =
    Teuchos::rcp( 
      new Anasazi::BasicEigenproblem<double, MV, OP>(te_K, te_M, te_ivec) 
    );

  // Inform the eigenproblem that the operator A is symmetric
  te_Problem->SetSymmetric(true);

  // Set the number of eigenvalues requested and the blocksize the solver should use
  te_Problem->SetNEV( nev );

  // Inform the eigenproblem that you are finishing passing it information
  info = te_Problem->SetProblem();
  if (info)
    cout << "Anasazi::BasicEigenproblem::SetProblem() returned with code : "<< info << endl;

  // Create the eigensolvers 
  Anasazi::BlockDavidson<double, MV, OP> te_Solver(te_Problem, MySM, MyOM, MyPL);

  returnCode = te_Solver.solve();
  if (returnCode != Anasazi::Ok) {
    testFailed = true;
  }
  if (!testFailed) {
    // Get the eigenvalues and eigenvectors from the eigenproblem
    Teuchos::RefCountPtr<std::vector<double> > evals = te_Problem->GetEvals();
    Teuchos::RefCountPtr<MV> evecs = te_Problem->GetEvecs();
    
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
    returnCode = OPT::Apply( *te_K, *evecs, *Kvec );
    assert( returnCode==Anasazi::Ok );
    // Compute M*evecs
    returnCode = OPT::Apply( *te_M, *evecs, *Mvec );
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
