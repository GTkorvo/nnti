#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <string>
#ifdef EPETRA_MPI
#include "mpi.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif

// Trilinos Objects
#include <Epetra_Comm.h>
#include <Epetra_SerialComm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Time.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

// NLSPack Objects
#include <NLS_Interface.h>
#include <NLS_PetraGroup.H>
#include <NLS_MethodManager.H>

// User specific files 
#include <fill.h>

using namespace std;

int main(int argc, char *argv[])
{
  int ierr = 0, i, j;
  bool debug = false;

  // Initialize MPI
#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
  int size = 1; // Serial case (not using MPI)
  int rank = 0;
#endif

#ifdef EPETRA_MPI
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  // Get the number of local equations from the command line
  if (argc!=2) { 
    cout << "Usage: " << argv[0] << " number_of_elements" << endl;
    exit(1);
  }
  int NumGlobalElements = atoi(argv[1]) + 1;
  int IndexBase = 0;

  if (NumGlobalElements < NumProc) {
    cout << "numGlobalBlocks = " << NumGlobalElements 
	 << " cannot be < number of processors = " << NumProc << endl;
    exit(1);
  }

  // Construct a Source Map that puts approximately the same 
  // Number of equations on each processor in uniform global ordering

  Epetra_Map StandardMap(NumGlobalElements, 0, Comm);
  int NumMyElements = StandardMap.NumMyElements();
  int StandardMyGlobalElements[NumMyElements];
  StandardMap.MyGlobalElements(StandardMyGlobalElements);

  int OverlapNumMyElements;
  int OverlapMinMyGID;

  OverlapNumMyElements = NumMyElements + 2;
  if ((MyPID == 0) || (MyPID == NumProc - 1)) 
    OverlapNumMyElements --;

  if (MyPID==0) 
    OverlapMinMyGID = StandardMap.MinMyGID();
  else 
    OverlapMinMyGID = StandardMap.MinMyGID() - 1;

  int OverlapMyGlobalElements[OverlapNumMyElements];

  for (i = 0; i < OverlapNumMyElements; i ++) 
    OverlapMyGlobalElements[i] = OverlapMinMyGID + i;

  Epetra_Map* tmpMap;
  if (size == 1) 
    tmpMap = new Epetra_Map(StandardMap);  
  else 
    tmpMap = new Epetra_Map(-1, OverlapNumMyElements, 
			    OverlapMyGlobalElements, 0, Comm);
  
  Epetra_Map& OverlapMap(*tmpMap);

  // Construct Linear Objects  
  Epetra_Import Importer(StandardMap, OverlapMap);
  Epetra_Vector soln(StandardMap);
  Epetra_CrsMatrix AA(Copy, StandardMap, 3);

  // Construct the Matrix
  Fill LO;  
  LO.registerFillObjects(StandardMap, OverlapMap, Importer, Comm);

  // Fix things so we can do multiple imports on objects (make graph static)
  // Get Matrix structure
  LO.fillMatrix(&soln, NULL, &AA);
  // Create a graph
  Epetra_CrsMatrix A(Copy, AA.Graph());
  A.TransformToLocal();

  // Re-register so that static matrix A is used
  LO.registerFillObjects(StandardMap, OverlapMap, Importer, Comm);

  // Initialize Solution
  i = soln.PutScalar(1.0);
  
  // Begin Nonlinear Solver ************************************

  // Create and set the parameters

  NLS_ParameterList nlParams;
  nlParams.setParameter("Nonlinear Solver", "Newton");   
  nlParams.setParameter("Output Level", 4);
  nlParams.setParameter("MyPID", MyPID); 

  // Sublist for convergence tests
  NLS_ParameterList& convParams = nlParams.sublist("Convergence Tests");
  convParams.setParameter("Max Iterations", 10);  
  convParams.setParameter("Absolute Residual Tolerance", 1.0e-6); 
  convParams.setParameter("Relative Residual Tolerance", 1.0e-4); 

  
  // Sublist for linear solver
  NLS_ParameterList& lsParams = nlParams.sublist("Linear Solver Parameters");
  lsParams.setParameter("Max Iterations", 400);  
  lsParams.setParameter("Tolerance", 1e-4); 

  // Sublist for line search
  NLS_ParameterList& searchParams = nlParams.sublist("Line Search Parameters");
  searchParams.setParameter("Method", "Interval Halving");
  searchParams.setParameter("Default Step", 10.0);

  // Create the interface between the test problem and the nonlinear solver
  // This is created by the user using inheritance of abstract base class:
  // NLS_PetraGroupInterface
  NLS_Interface Interface;
  Interface.registerFill(&LO);

  // Crate the shared Jacobian
  NLS_PetraSharedJacobian SharedA(A);

  // Create the Groups 
  NLS_PetraGroup initialGuess(soln, SharedA, Interface);
  NLS_PetraGroup solutionSpace(soln, SharedA, Interface);  

  // Create the method
  NLS_MethodManager Newton(initialGuess, solutionSpace, nlParams);
  ierr = Newton.solve();

  // If converged, get the group with the final solution
  if (ierr >= 0) {
    NLS_Group& finalSolution(Newton.getSolutionGroup());
    // Copy the solution into a petra vector for printing
    soln = dynamic_cast<const NLS_PetraVector&> 
      (finalSolution.getX()).getPetraVector();
  }  

  // End Nonlinear Solver **************************************

  // Print solution
  char file_name[25];
  FILE *ifp;
  (void) sprintf(file_name, "output.%d",MyPID);
  ifp = fopen(file_name, "w");
  for (i=0; i<NumMyElements; i++)
    fprintf(ifp, "%d  %E\n", StandardMap.MinMyGID()+i, soln[i]);
  fclose(ifp);

  //cout << "Final Solution" << (&soln)[0] <<endl;

  delete tmpMap;


  //delete &StandardGraph;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
