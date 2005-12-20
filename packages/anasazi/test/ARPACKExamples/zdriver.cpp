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
// This test compares the Anasazi solvers against ARPACK. The eigenproblems
// used are from the ARPACK examples: SYM, NONSYM, and COMPLEX
// See ARPACK_Operators.hpp and examlpesdesc for more information.

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCG.hpp"
#include "AnasaziBlockDavidson.hpp"
// #include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBasicSort.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "AnasaziMVOPTester.hpp"

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif

// templated multivector 
#include "MyMultiVec.hpp"
// ARPACK test problems
#include "ARPACK_Operators.hpp"

using namespace Teuchos;

int main(int argc, char *argv[]) 
{
  int info = 0;

#ifdef EPETRA_MPI
  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();

  bool testFailed;
  bool verbose = 0;
  std::string which("auto");
  int nx = 10;
  std::string problem("SDRV1");
  bool isherm;
  std::string solver("auto");

  CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("sort",&which,"Targetted eigenvalues (auto, SR, LR, SI, LI, SM, LM).");
  cmdp.setOption("nx",&nx,"Number of interior elements.");
  cmdp.setOption("problem",&problem,"Problem to solve.");
  cmdp.setOption("solver",&solver,"Eigensolver to use (LOBPCG, BKS, BD, auto)");
  if (cmdp.parse(argc,argv) != CommandLineProcessor::PARSE_SUCCESSFUL) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }

#ifdef HAVE_COMPLEX
  typedef std::complex<double> ST;
#elif HAVE_COMPLEX_H
  typedef ::complex<double> ST;
#else
  typedef double ST;
  // no complex. quit with failure.
  if (verbose && MyPID == 0) {
    cout << "Not compiled with complex support." << endl;
    cout << "End Result: TEST FAILED" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
#endif
  typedef ScalarTraits<ST>                   SCT;
  typedef SCT::magnitudeType                  MT;
  typedef Anasazi::MultiVec<ST>               MV;
  typedef Anasazi::Operator<ST>               OP;
  typedef Anasazi::MultiVecTraits<ST,MV>     MVT;
  typedef Anasazi::OperatorTraits<ST,MV,OP>  OPT;
  ST ONE  = SCT::one();
  ST ZERO = SCT::zero();


  // Create default output manager 
  RefCountPtr<Anasazi::OutputManager<ST> > MyOM 
    = rcp( new Anasazi::OutputManager<ST>( MyPID ) );
  // Set verbosity level
  if (verbose) {
    MyOM->SetVerbosity( Anasazi::Warning + Anasazi::FinalSummary );
  }

  if (MyOM->isVerbosityAndPrint(Anasazi::Warning)) {
    cout << Anasazi::Anasazi_Version() << endl << endl;
  }

  Anasazi::ReturnType returnCode = Anasazi::Ok;  

  // Eigensolver parameters
  int dim = nx*nx;
  int nev = 4;
  int blockSize = 4;
  int maxIters = 500;
  int maxBlocks = 5;
  MT tol = 1.0e-10;

  // Create parameter list to pass into solver
  ParameterList MyPL;
  MyPL.set( "Block Size", blockSize );
  MyPL.set( "Max Iters", maxIters );
  MyPL.set( "Tol", tol );
  MyPL.set( "Max Blocks", maxBlocks );
  MyPL.set( "Max Restarts", maxIters );

  // Create initial vectors
  RefCountPtr<MV> ivec = rcp( new MyMultiVec<ST>(dim,blockSize) );
  ivec->MvRandom();

  // Create matrices
  RefCountPtr< ARPACK_Example<ST> > prob;
  RefCountPtr<OP> A, M, Op;

  prob = GetARPACKExample<ST>(problem,dim);
  if (!prob.get()) {
    if ( MyOM->isVerbosityAndPrint(Anasazi::Warning) ) {
      cout << "Invalid driver name. Try something like ""ndrv3"" or ""sdrv2""." << endl;
      cout << "End Result: TEST FAILED" << endl;	
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  A = prob->getA();
  M = prob->getM();
  Op = prob->getOp();
  isherm = prob->isHerm();

  // determine solver
  if (solver == "auto") {
    if (isherm) {
      solver = "LOBPCG";
    }
    else {
      solver = "BKS";
    }
  }

  // determine sort
  if (which == "auto") {
    if (solver == "LOBPCG" || solver == "BD") {
      which = "SM";
    }
    else {
      which = prob->getSort();
    }
  }

  // test multivector and operators
  int ierr;
  cout << "Testing MultiVector" << endl;
  ierr = Anasazi::TestMultiVecTraits<ST,MV>(MyOM,ivec);
  if (ierr != Anasazi::Ok) {
    cout << "MultiVec failed TestMultiVecTraits()" << endl;
  }
  cout << "Testing OP" << endl;
  ierr = Anasazi::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,A);
  if (ierr != Anasazi::Ok) {
    cout << "OP failed TestOperatorTraits()" << endl;
  }
  cout << "Testing M" << endl;
  ierr = Anasazi::TestOperatorTraits<ST,MV,OP>(MyOM,ivec,M);
  if (ierr != Anasazi::Ok) {
    cout << "M failed TestOperatorTraits()" << endl;
  }

  // Create the sort manager
  RefCountPtr<Anasazi::BasicSort<ST,MV,OP> > MySM = 
     rcp( new Anasazi::BasicSort<ST,MV,OP>(which) );

  // Create eigenproblem
  RefCountPtr<Anasazi::BasicEigenproblem<ST,MV,OP> > MyProblem;
  if (solver == "LOBPCG" || solver == "BD") {
    MyProblem = rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(A, M, ivec) );
  }
  else {
    MyProblem = rcp( new Anasazi::BasicEigenproblem<ST,MV,OP>(Op, M, ivec) );
  }
  // Inform the eigenproblem that the operator A is symmetric
  MyProblem->SetSymmetric(isherm);

  // Set the number of eigenvalues requested and the blocksize the solver should use
  MyProblem->SetNEV( nev );

  // Inform the eigenproblem that you are done passing it information
  info = MyProblem->SetProblem();
  if (info) {
    if (MyOM->isVerbosityAndPrint(Anasazi::Warning)) {
      cout << "Anasazi::BasicEigenproblem::SetProblem() returned with code : "<< info << endl;
      cout << "End Result: TEST FAILED" << endl;	
    }
#ifdef EPETRA_MPI
    MPI_Finalize() ;
#endif
    return -1;
  }

  // Create the eigensolver 
  RefCountPtr< Anasazi::Eigensolver<ST,MV,OP> > MySolver;
  if (solver == "LOBPCG") {
    MySolver = rcp( new Anasazi::LOBPCG<ST,MV,OP>(MyProblem, MySM, MyOM, MyPL));
  }
  /*
  else if (solver == "BKS") {
    MyPL.set( "Block Size", 1 );
    MyPL.set( "Max Blocks", 20 );
    MySolver = rcp( new Anasazi::BlockKrylovSchur<ST,MV,OP>(MyProblem, MySM, MyOM, MyPL));
  }
  */
  else if (solver == "BD") {
    MySolver = rcp( new Anasazi::BlockDavidson<ST,MV,OP>(MyProblem, MySM, MyOM, MyPL));
  }
  else {
    if (MyOM->isVerbosityAndPrint(Anasazi::Warning)) {
      cout << "Invalid solver: " << solver << endl;
      cout << "End Result: TEST FAILED" << endl;	
    }
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return -1;
  }
  if ( MyOM->isVerbosityAndPrint(Anasazi::Warning) ) {
    cout << "Using solver: " << solver << endl;
  }

  // Solve the problem to the specified tolerances or length
  returnCode = MySolver->solve();
  testFailed = false;
  if (returnCode != Anasazi::Ok) {
    testFailed = true; 
  }

  // Get the eigenvalues and eigenvectors from the eigenproblem
  RefCountPtr<std::vector<ST> > evals = MyProblem->GetEvals();
  RefCountPtr<MV > evecs = MyProblem->GetEvecs();
  int nevecs = MVT::GetNumberVecs(*evecs);

  // Perform spectral transform on eigenvalues, if we used the 
  // spectral xformed operator (Op)
  if (solver == "BKS") {
    prob->xformeval(*evals);
  }

  // Compute the direct residual
  std::vector<MT> normV( nevecs );
  SerialDenseMatrix<int,ST> L(nevecs,nevecs);
  for (int i=0; i<nevecs; i++) {
    L(i,i) = (*evals)[i];
  }
  RefCountPtr<MV > Avecs = MVT::Clone( *evecs, nevecs );
  RefCountPtr<MV > Mvecs = MVT::Clone( *evecs, nevecs );

  OPT::Apply( *A, *evecs, *Avecs );
  OPT::Apply( *M, *evecs, *Mvecs );
  // Compute A*evecs - M*evecs*L
  MVT::MvTimesMatAddMv( -ONE, *Mvecs, L, ONE, *Avecs );
  MVT::MvNorm( *Avecs, &normV );

  // check residuals
  for (int i=0; i<nevecs; i++) {
    if ( (*evals)[i] != SCT::zero() ) {
      normV[i] = SCT::magnitude(normV[i]/(*evals)[i]);
    }
    if ( normV[i] > ((MT)2.0)*tol ) {
      testFailed = true;
    }
  }

  if ( MyOM->isVerbosityAndPrint(Anasazi::FinalSummary) ) {
    //      28,5,22
    cout << "Back transformed eigenvalues     Relative Residual Norm" << endl
         << "-------------------------------------------------------" << endl;
    for (int i=0; i<nevecs; i++) {
      cout.setf(ios::scientific, ios::floatfield);  
      cout.precision(10);
      cout << std::setw(28) << std::right << (*evals)[i] 
           << "     "
           << std::setw(22) << std::right << normV[i] 
           << endl;
    }
  }

  // Exit
#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

  if (testFailed) {
    if (MyOM->isVerbosityAndPrint(Anasazi::Warning)) {
      cout << "End Result: TEST FAILED" << endl;	
    }
    return -1;
  }
  //
  // Default return value
  //
  if (MyOM->isVerbosityAndPrint(Anasazi::Warning)) {
    cout << "End Result: TEST PASSED" << endl;
  }
  return 0;

}	
