#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"


#ifdef PACKAGE
#undef PACKAGE
#endif

#ifdef PACKAGE_NAME
#undef PACKAGE_NAME
#endif

#ifdef PACKAGE_BUGREPORT
#undef PACKAGE_BUGREPORT
#endif

#ifdef PACKAGE_STRING
#undef PACKAGE_STRING
#endif

#ifdef PACKAGE_TARNAME
#undef PACKAGE_TARNAME
#endif

#ifdef PACKAGE_VERSION
#undef PACKAGE_VERSION
#endif

#ifdef VERSION
#undef VERSION
#endif

#include "Teuchos_ParameterList.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "AztecOO.h"

#include "Trilinos_Util_CrsMatrixGallery.h"

static bool verbose = false;
static bool ml_verbose = false;

void PrintLine() 
{
  cout << endl;
  for( int i=0 ; i<80 ; ++i )
    cout << "=";
  cout << endl;
  cout << endl;
  
  return;
}

int TestMultiLevelPreconditioner(char ProblemType[],
				 Teuchos::ParameterList & MLList,
				 Epetra_LinearProblem & Problem, double & TotalErrorResidual,
				 double & TotalErrorExactSol)
{
  
  Epetra_MultiVector   * lhs     = Problem.GetLHS();
  Epetra_MultiVector   * rhs     = Problem.GetRHS();
  Epetra_RowMatrix     * A       = Problem.GetMatrix();
  
  // ======================================== //
  // create a rhs corresponding to lhs or 1's //
  // ======================================== //
  
  lhs->PutScalar(1.0);
  A->Multiply(false,*lhs,*rhs);

  lhs->PutScalar(0.0);
  
  Epetra_Time Time(A->Comm());
  
  // =================== //
  // call ML and AztecOO //
  // =================== //
  
  AztecOO solver(Problem);
  
  if (ml_verbose)
    MLList.set("output", 10);
  else
    MLList.set("output", 0);

  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);
  
  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  
  solver.SetAztecOption(AZ_solver, AZ_gmres);
  solver.SetAztecOption(AZ_output, 32);
  
  solver.SetAztecOption(AZ_kspace, 160);
  
  solver.Iterate(1550, 1e-12);
  
  delete MLPrec;
  
  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //
  
  double d = 0.0, d_tot = 0.0;
  
  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    d += ((*lhs)[0][i] - 1.0) * ((*lhs)[0][i] - 1.0);
  
  A->Comm().SumAll(&d,&d_tot,1);
  
  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  
  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A->Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);
  
  string msg = ProblemType;
  
  if(verbose) {
    cout << msg << "......Using " << A->Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }
  
  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;
  
  return( solver.NumIters() );
  
}

using namespace Trilinos_Util;

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  if (Comm.MyPID() == 0)
    verbose = true;

  if (argc >= 2)
    if (strcmp(argv[1], "-v") == 0) 
      ml_verbose = true;

  // initialize the random number generator

  srandom((unsigned int)1);

  int NumProcs = Comm.NumProc();
  
  // ===================== //
  // create linear problem //
  // ===================== //
  
  CrsMatrixGallery Gallery("recirc_2d", Comm);
  Gallery.Set("problem_size", 10000);
  
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  int TotalFailed = 0;
  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;

  int iters;
  
  // ================== //
  // no-default options //
  // ================== //

  if (true) {
    
    if( Comm.MyPID() == 0 ) PrintLine();
    
    Teuchos::ParameterList MLList;
    iters = TestMultiLevelPreconditioner("no defaults", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );
    
    if (verbose) {
#ifdef HAVE_ML_AMESOS
      // expected iterations
      switch( NumProcs ) {
      case 1:
        if( iters != 23 ) {
          ++TotalFailed;
          cout << endl << "### TEST FAILED : expecting 23 iterations, got " 
               << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      case 4:
        if( iters != 29 ) {
          ++TotalFailed;
          cout << endl << "### TEST FAILED : expecting 29 iterations, got " 
               << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      }
#else
      cout << "### Cannot check the number of iterations (no `enable-amesos')" << endl;
      cout << "### Checking linear system residual only (at the bottom of the file)" << endl;
#endif
    }
  }

  // ====================== //
  // default options for DD //
  // ====================== //

  if (true) {
    
    if( Comm.MyPID() == 0 ) PrintLine();
    
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("DD",MLList);
    iters = TestMultiLevelPreconditioner("DD", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    if (verbose) {
#ifdef HAVE_ML_AMESOS
      // expected iterations
      switch( NumProcs ) {
      case 1:
        if( iters != 70 ) {
          ++TotalFailed;
          cerr << endl << "### TEST FAILED : expecting 70 iterations, got "
            << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      case 4:
        if( iters != 82 ) {
          ++TotalFailed;
          cerr << endl << "### TEST FAILED : expecting 82 iterations, got "
            << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      }
#else
      cout << "### Cannot check the number of iterations (no `enable-amesos')" << endl;
      cout << "### Checking linear system residual only (at the bottom of the file)" << endl;
#endif
    }
  }

  // ========================================== //
  // default options for DD -- 16 aggr per proc //
  // ========================================== //

  if (true) {

    if( Comm.MyPID() == 0 ) PrintLine();

    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("DD",MLList);
    MLList.set("aggregation: local aggregates", 16);
    iters = TestMultiLevelPreconditioner("DD", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    if (verbose) {
#if defined(HAVE_ML_AMESOS) && defined(HAVE_ML_METIS)
      // expected iterations
      switch( NumProcs ) {
      case 1:
        if( iters != 65 ) {
          ++TotalFailed;
          cerr << endl << "### TEST FAILED : expecting 65 iterations, got "
            << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      case 4:
        if( iters != 66 ) {
          ++TotalFailed;
          cerr << endl << "### TEST FAILED : expecting 66 iterations, got "
            << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      }
#else
      cout << "### Cannot check the number of iterations" << endl;
      cout << "### (no `enable-amesos --enable-metis')" << endl;
      cout << "### Checking linear system residual only (at the bottom of the file)" << endl;
#endif
    }
  }

  // ========================= //
  // default options for DD-ML //
  // ========================= //
  
  if (true) {

    if( Comm.MyPID() == 0 ) PrintLine();

    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("DD-ML",MLList);
    iters = TestMultiLevelPreconditioner("DD-ML", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    if (verbose) {
#if defined(HAVE_ML_AMESOS) && defined(HAVE_ML_METIS)
      // can check iteration number only with Amesos coarse solver,
      // and METIS installed
      // expected iterations
      switch( NumProcs ) {
      case 1:
        if( iters != 64 ) {
          ++TotalFailed;
          cerr << endl << "### TEST FAILED : expecting 64 iterations, got "
            << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      case 4:
        if( iters != 77 ) {
          ++TotalFailed;
          cerr << endl << "### TEST FAILED : expecting 77 iterations, got "
            << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      }
#else
      cout << "### Cannot check the number of iterations" << endl;
      cout << "### (no `enable-amesos --enable-metis')" << endl;
      cout << "### Checking linear system residual only (at the bottom of the file)" << endl;
#endif
    }
  }

  // ========================= //
  // default options for DD-ML //
  // ========================= //

  if (true) {
    
    if( Comm.MyPID() == 0 ) PrintLine();
    
    Teuchos::ParameterList MLList;
    ML_Epetra::SetDefaults("DD-ML",MLList);
    MLList.set("aggregation: nodes per aggregate (level 0)", 64);
    MLList.set("aggregation: nodes per aggregate (level 1)", 27);
    iters = TestMultiLevelPreconditioner("DD-ML", MLList, *Problem, TotalErrorResidual, TotalErrorExactSol );

    if (verbose) {
#if defined(HAVE_ML_AMESOS) && defined(HAVE_ML_METIS) && defined(HAVE_ML_PARMETIS_3x)
      // can check iteration number only with Amesos coarse solver,
      // and METIS, ParMETIS installed
      // expected iterations
      switch( NumProcs ) {
      case 1:
        if( iters != 48 ) {
          ++TotalFailed;
          cerr << endl << "### TEST FAILED : expecting 48 iterations, got "
               << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      case 4:
        // it seems that parmetis has some random stuff goingon
        // inside. At least on stratus and on s850675 I
        // got different iteration counts, 57 vs 59.
        if( iters > 100 ) {
          ++TotalFailed;
            cerr << endl << "### TEST FAILED : expecting less than 100 iterations, got "
                 << iters << endl << endl;
        }
        else
          cout << endl << "### TEST PASSED" << endl << endl;
        break;
      }
#else
      cout << "### Cannot check the number of iterations" << endl;
      cout << "### (no `enable-amesos --enable-parmetis3x')" << endl;
      cout << "### Checking linear system residual only (at the bottom of the file)" << endl;
#endif
    }
  }

  // ===================== //
  // print out total error //
  // ===================== //

  if( Comm.MyPID() == 0 ) {
    cout << endl;
    cout << "......Total error for residual        = " << TotalErrorResidual << endl;
    cout << "......Total error for exact solution  = " << TotalErrorExactSol << endl;
    cout << "......Total # of failed tests         = " << TotalFailed << endl;
    cout << endl;
  }

  if( TotalErrorResidual > 1e-8 ) return( EXIT_FAILURE );
#ifdef HAVE_ML_AMESOS
  if( TotalFailed ) return( EXIT_FAILURE );
#endif

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (verbose)
    cout << endl << "### ALL TESTS PASSED" << endl;

  return( EXIT_SUCCESS );

}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{

  // still need to deal with MPI, some architecture don't like
  // an exit(0) without MPI_Finalize()
#ifdef ML_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");

#ifdef ML_MPI
  MPI_Finalize();
#endif

  return(0);
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) HAVE_ML_AZTECOO */
