#ifdef HAVE_MPI
/*--------------------------------------------------------------------*/
#include <time.h>
#include "mpi.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"

//This program was contributed by a user (Eric Marttila) to
//demonstrate a performance problem when filling an entire
//FECrsMatrix on a single processor and then letting
//GlobalAssemble distribute the data according to the Epetra_Map
//that the matrix was constructed with.
//
//Using profile data generated by running this program, we
//made changes to FECrsMatrix and improved the performance of
//inserting non-local data.

int main(int argCount, char **argValue)
{
  int ierr;
  MPI_Init(&argCount,&argValue);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
  const int rank = Comm.MyPID();

  // Construct a Map 
  int nGlobalElements = 10000;//10,000 is deliberately small for nightly testing purposes.
                             // Set to 1 million for performance testing.

  Epetra_Map Map(nGlobalElements, 0, Comm);

  // Create a matrix
  Epetra_FECrsMatrix A(Copy, Map, 1);

  time_t startTime = 0;
  if (rank == 0) {
    startTime = time(0);
  }

  // Fill matrix on the master process
  if (rank == 0) {
    double values[1];
    int    indices[1];
    const int numEntries = 1;

    for (int globalRowIdx=0; globalRowIdx<nGlobalElements; ++globalRowIdx) {
      indices[0] = globalRowIdx;
      values[0] = 3.2 + globalRowIdx*0.01;

      if (globalRowIdx % 5000 == 0) {
  cerr << "About to insert row " << globalRowIdx << "\n";
      }

      ierr = A.InsertGlobalValues( globalRowIdx, numEntries,
           (const double *)&values[0],
           (const int *)&indices[0] );

      assert(ierr==0);
    }
  }

  double insertionTime = 0;
  if (rank == 0) {
    time_t endTime = time(0);
    insertionTime = difftime(endTime, startTime);
  }

  // Finish up
  ierr = A.GlobalAssemble();
  assert(ierr==0);

  if (rank == 0) {
    cerr << "insertion time = " << insertionTime << " (seconds)\n";
  }


  MPI_Finalize();

  return 0;
}
/*--------------------------------------------------------------------*/
#else
int main(int,char**) {
  return 0;
}
#endif

