// CrsMatrix_MapColoring Test routine

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_IntVector.h"
#include "Epetra_MapColoring.h"
#include "EDT_CrsGraph_MapColoring.h"
#include "EDT_CrsGraph_MapColoringIndex.h"
#include "../epetra_test_err.h"
#include <stdio.h>

#include <vector>

int main(int argc, char *argv[]) {

  int i, ierr=0, returnierr=0;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init( &argc, &argv );
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;

#endif

  bool verbose = false;

  if( argc < 2 || argc > 4 )
  {
    cout << "Usage: " << argv[0] << " [-v] nx [ny]" << endl;
    exit(1);
  }

  int loc = 1;
  // Check if we should print results to standard out
  if(argv[loc][0]=='-' && argv[loc][1]=='v')
  { verbose = true; ++loc; }

  int nx = atoi( argv[loc++] );
  int ny = 1;
  if( argc == (loc+1) ) ny = atoi( argv[loc] );

#ifdef EPETRA_MPI
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();

  if(verbose) cout << Comm << endl << flush;
  Comm.Barrier();

  bool verbose1 = false;
  if(verbose) verbose1 = (MyPID==0);

  int NumGlobalElements = nx * ny;
  if( NumGlobalElements < NumProc )
  {
    cout << "NumGlobalElements = " << NumGlobalElements <<
            " cannot be < number of processors = " << NumProc;
    exit(1);
  } 
	
  int IndexBase = 0;
  Epetra_Map Map( NumGlobalElements, IndexBase, Comm );

  // Extract the global indices of the elements local to this processor
  int NumMyElements = Map.NumMyElements();
  std::vector<int> MyGlobalElements( NumMyElements );
  Map.MyGlobalElements( &MyGlobalElements[0] );
  if( verbose ) cout << Map;

  // Create the number of non-zeros for a tridiagonal (1D problem) or banded
  // (2D problem) matrix
  std::vector<int> NumNz( NumMyElements, 5 );
  int global_i;
  int global_j;
  for (int i = 0; i < NumMyElements; ++i)
  {
    global_j = MyGlobalElements[i] / nx;
    global_i = MyGlobalElements[i] - global_j * nx;
    if (global_i == 0)    NumNz[i] -= 1;  // By having separate statements,
    if (global_i == nx-1) NumNz[i] -= 1;  // this works for 2D as well as 1D
    if (global_j == 0)    NumNz[i] -= 1;  // systems (i.e. nx x 1 or 1 x ny)
    if (global_j == ny-1) NumNz[i] -= 1;  // or even a 1 x 1 system
  }
  if(verbose)
  { 
    cout << endl << "NumNz: ";
    for (int i = 0; i < NumMyElements; i++) cout << NumNz[i] << " ";
    cout << endl;
  } // end if
  
  // Create the Epetra Compressed Row Sparse Graph
  Epetra_CrsGraph A( Copy, Map, &NumNz[0] );
  
  std::vector<int> Indices(5);
  int NumEntries;
  
  for (int i = 0; i < NumMyElements; ++i )
  {
    global_j = MyGlobalElements[i] / nx;
    global_i = MyGlobalElements[i] - global_j * nx;
    NumEntries = 0;
    // (i,j-1) entry
    if (global_j > 0 && ny > 1)
      Indices[NumEntries++] = global_i   + (global_j-1)*nx;
    // (i-1,j) entry
    if (global_i > 0)
      Indices[NumEntries++] = global_i-1 +  global_j   *nx;
    // (i,j) entry
    Indices[NumEntries++] = MyGlobalElements[i];
    // (i+1,j) entry
    if (global_i < nx-1)
      Indices[NumEntries++] = global_i+1 +  global_j   *nx;
    // (i,j+1) entry
    if (global_j < ny-1 && ny > 1)
      Indices[NumEntries++] = global_i   + (global_j+1)*nx;

    // Insert the global indices
    assert( A.InsertGlobalIndices( MyGlobalElements[i], NumEntries, &Indices[0] )  == 0 );
  } // end i loop

  // Finish up graph construction
  assert(A.TransformToLocal() == 0);

  EpetraExt::CrsGraph_MapColoring MapColoringTransform(verbose);
  Epetra_MapColoring & ColorMap = MapColoringTransform( A );

  int NumColors = ColorMap.NumColors();
  int * ListOfColors = ColorMap.ListOfColors();

  EpetraExt::CrsGraph_MapColoringIndex MapColoringIndexTransform( ColorMap );
  vector<Epetra_IntVector> & ColIndices = MapColoringIndexTransform( A );

  if( verbose )
  {

    cout << endl;
    cout << "***************************************\n";
    cout << "Column Indexing by Color:\n";
    cout << "***************************************\n";
    cout << endl;
    for( int i = 0; i < NumColors; ++i )
    {
      cout << "COLOR: " << ListOfColors[i] << endl;
      cout << ColIndices[i];
    }
    cout << endl;
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return returnierr;
}

