#include "Amesos_ConfigDefs.h"
#ifdef HAVE_AMESOS_SUPERLUDIST

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_Util.h"
#include "Amesos_Superludist.h"
#include "Amesos_TestRowMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include <vector>
using namespace Trilinos_Util;

//=============================================================================
bool CheckError(const Epetra_RowMatrix& A,
		const Epetra_MultiVector& x,
		const Epetra_MultiVector& b,
		const Epetra_MultiVector& x_exact)
{
  vector<double> Norm;
  int NumVectors = x.NumVectors();
  Norm.resize(NumVectors);
  Epetra_MultiVector Ax(x);
  A.Multiply(false,x,Ax);
  Ax.Update(1.0,b,-1.0);
  Ax.Norm2(&Norm[0]);
  bool TestPassed = false;
  double TotalNorm = 0.0;
  for (int i = 0 ; i < NumVectors ; ++i) {
    TotalNorm += Norm[i];
  }
  if (A.Comm().MyPID() == 0)
    cout << "||Ax - b||  = " << TotalNorm << endl;
  if (TotalNorm < 1e-5 )
    TestPassed = true;
  else
    TestPassed = false;

  Ax.Update (1.0,x,-1.0,x_exact,0.0);
  Ax.Norm2(&Norm[0]);
  for (int i = 0 ; i < NumVectors ; ++i) {
    TotalNorm += Norm[i];
  }
  if (A.Comm().MyPID() == 0)
    cout << "||Ax - b||  = " << TotalNorm << endl;
  if (TotalNorm < 1e-5 )
    TestPassed = true;
  else
    TestPassed = false;

  return(TestPassed);
}

int sub_main( Epetra_Comm &Comm ) { 
  //  Allow destruction of the Amesos class(es) before the
  //  call to MPI_Finalize()

  int NumGlobalElements = 100;   // kludge  see bug #1142
  int NumVectors = 7;
  
  // =================== //
  // create a random map //
  // =================== //
 
  int* part = new int[NumGlobalElements];

  if (Comm.MyPID() == 0) {
    Epetra_Util Util;

    for( int i=0 ; i<NumGlobalElements ; ++i ) {
      unsigned int r = Util.RandomInt();	
      part[i] = r%(Comm.NumProc());
    }
  }

  Comm.Broadcast(part,NumGlobalElements,0);

  // count the elements assigned to this proc
  int NumMyElements = 0;
  for (int i = 0 ; i < NumGlobalElements ; ++i) {
    if (part[i] == Comm.MyPID()) 
      NumMyElements++;
  }

  // get the loc2global list
  int* MyGlobalElements = new int[NumMyElements];
  int count = 0;
  for (int i = 0 ; i < NumGlobalElements ; ++i) {
    if (part[i] == Comm.MyPID() ) 
      MyGlobalElements[count++] = i;
  }

  Epetra_Map Map(NumGlobalElements,NumMyElements,MyGlobalElements,
		 0,Comm);

  delete [] part;

  // ===================== //
  // Create a dense matrix //
  // ===================== //
 
  Epetra_CrsMatrix Matrix(Copy,Map,NumGlobalElements);

  int* Indices = new int[NumGlobalElements];
  double* Values = new double[NumGlobalElements];

  for (int i = 0 ; i < NumGlobalElements ; ++i) 
    Indices[i] = i;

  for (int i = 0 ; i < NumMyElements ; ++i) {
    int iGlobal = MyGlobalElements[i];
    for (int jGlobal = 0 ; jGlobal < NumGlobalElements ; ++jGlobal) {
      if (iGlobal == jGlobal) 
	Values[jGlobal] = 1.0 * (NumGlobalElements + 1 ) *
	  (NumGlobalElements + 1);
      else if (iGlobal > jGlobal)
	Values[jGlobal] = -1.0*(jGlobal+1);
      else
	Values[jGlobal] = 1.0*(iGlobal+1);
    }
    Matrix.InsertGlobalValues(MyGlobalElements[i],
                              NumGlobalElements, Values, Indices);

  }

  Matrix.FillComplete();

  delete [] MyGlobalElements;
  delete [] Indices;
  delete [] Values;
 
  // ======================== //
  // other data for this test //
  // ======================== //

  Amesos_TestRowMatrix A(&Matrix);
  Epetra_MultiVector x(Map,NumVectors);
  Epetra_MultiVector x_exact(Map,NumVectors);
  Epetra_MultiVector b(Map,NumVectors);
  x_exact.Random();
  A.Multiply(false,x_exact,b);

  // =========== //
  // AMESOS PART //
  // =========== //

  Epetra_LinearProblem Problem;
  Amesos_Superludist Solver(Problem);

  Problem.SetOperator(&A);
  Problem.SetLHS(&x);
  Problem.SetRHS(&b);

  AMESOS_CHK_ERR(Solver.SymbolicFactorization());
  AMESOS_CHK_ERR(Solver.NumericFactorization());
#if 1
  AMESOS_CHK_ERR(Solver.Solve());
#endif 

  bool TestPassed = true;

  TestPassed = TestPassed &&
    CheckError(A,x,b,x_exact);

  if (TestPassed) {
    if (Comm.MyPID() == 0)
      cout << endl << "TEST PASSED" << endl << endl;
  }
  else {
    if (Comm.MyPID() == 0)
      cout << endl << "TEST FAILED" << endl << endl;
  }

  AMESOS_CHK_ERR( ! TestPassed ) ; 

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);
}

//=============================================================================
int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int retvalue = sub_main(Comm) ; 

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return( retvalue ) ;   
}

#else

// Triutils is not available. Sorry, we have to give up.

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#else
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  puts("Please configure AMESOS with --enable-amesos-superludist");
  puts("to run this example");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(0);
}

#endif


