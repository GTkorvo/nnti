// Epetra_Comm Test routine

#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include <mpi.h>
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "checkmap.h"

int main(int argc, char *argv[]) {


#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;


  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << Comm << endl;

  bool verbose1 = verbose;
  if (verbose) verbose = (MyPID==0);

  int NumMyElements = 10000;
  int NumMyElements1 = NumMyElements; // Used for local map
  int NumGlobalElements = NumMyElements*NumProc+EPETRA_MIN(NumProc,3);
  if (MyPID < 3) NumMyElements++;
  int IndexBase = 0;
  bool DistributedGlobal = (NumGlobalElements>NumMyElements);
  
  Epetra_Map * Map;

  // Test exceptions

  if (verbose)
    cout << "*******************************************************************************************" << endl
	 << "        Testing Exceptions (Expect error messages if EPETRA_NO_ERROR_REPORTS is not defined" << endl
	 << "*******************************************************************************************" << endl
	 << endl << endl;

  try {
    if (verbose) cout << "Checking Epetra_Map(-2, IndexBase, Comm)" << endl;
    Map = new Epetra_Map(-2, IndexBase, Comm);
  }
  catch (int Error) {
    if (Error!=-1) {
      if (verbose) cout << "Error code = " << Error << "Should be -1" << endl;
      return 1;
    }
    if (verbose) cout << "Checked OK\n\n" << endl;
  }

  try {
    if (verbose) cout << "Checking Epetra_Map(2, 3, IndexBase, Comm)" << endl;
    Map = new Epetra_Map(2, 3, IndexBase, Comm);
  }
  catch (int Error) {
    if (Error!=-4) {
      if (verbose) {
	cout << "Error code = " << Error << "Should be -4" << endl;
      }
      return 1;
    }
    if (verbose) {
      cout << "Checked OK\n\n" << endl;
    }
  }

  if (verbose) cerr << flush;
  if (verbose) cout << flush;
  Comm.Barrier();
  if (verbose)
    cout << endl << endl
      << "*******************************************************************************************" << endl
      << "        Testing valid constructor now......................................................" << endl
      << "*******************************************************************************************" << endl
      << endl << endl;

  // Test Epetra-defined uniform linear distribution constructor
  Map = new Epetra_Map(NumGlobalElements, IndexBase, Comm);
  if (verbose) cout << "Checking Epetra_Map(NumGlobalElements, IndexBase, Comm)" << endl;
  int ierr = checkmap(*Map, NumGlobalElements, NumMyElements, 0, 
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  delete Map;

  // Test User-defined linear distribution constructor
  Map = new Epetra_Map(NumGlobalElements, NumMyElements, IndexBase, Comm);

  if (verbose) cout << "Checking Epetra_Map(NumGlobalElements, NumMyElements, IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, 0, 
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  delete Map;

  // Test User-defined arbitrary distribution constructor
  // Generate Global Element List.  Do in reverse for fun!

  int * MyGlobalElements = new int[NumMyElements];
  int MaxMyGID = (Comm.MyPID()+1)*NumMyElements-1+IndexBase;
  if (Comm.MyPID()>2) MaxMyGID+=3;
  for (int i = 0; i<NumMyElements; i++) MyGlobalElements[i] = MaxMyGID-i;

  Map = new Epetra_Map(NumGlobalElements, NumMyElements, MyGlobalElements, 
		      IndexBase, Comm);
  if (verbose) cout << "Checking Epetra_Map(NumGlobalElements, NumMyElements, MyGlobalElements,  IndexBase, Comm)" << endl;
  ierr = checkmap(*Map, NumGlobalElements, NumMyElements, MyGlobalElements, 
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  // Test Copy constructor
  Epetra_Map * Map1 = new Epetra_Map(*Map);

  if (verbose) cout << "Checking Epetra_Map(*Map)" << endl;
  ierr = checkmap(*Map1, NumGlobalElements, NumMyElements, MyGlobalElements, 
		  IndexBase, Comm, DistributedGlobal);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  Epetra_Map * SmallMap = 0;
  if (verbose1) {
    // Build a small map for test cout.  Use 10 elements from current map
    int * MyEls = Map->MyGlobalElements();
    int IndBase = Map->IndexBase();
    int MyLen = EPETRA_MIN(10+Comm.MyPID(),Map->NumMyElements());
    SmallMap = new Epetra_Map(-1, MyLen, MyEls, IndBase, Comm);
  }

  delete [] MyGlobalElements;
  delete Map;
  delete Map1;


  // Test LocalMap constructor
  Epetra_LocalMap *LocalMap = new Epetra_LocalMap(NumMyElements1, IndexBase, 
						Comm);
  if (verbose) cout << "Checking Epetra_LocalMap(NumMyElements1, IndexBase, Comm)" << endl;
  ierr = checkmap(*LocalMap, NumMyElements1, NumMyElements1, 0, 
		  IndexBase, Comm, false);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  // Test Copy constructor
  Epetra_LocalMap * LocalMap1 = new Epetra_LocalMap(*LocalMap);
  if (verbose) cout << "Checking Epetra_LocalMap(*LocalMap)" << endl;
  ierr = checkmap(*LocalMap1, NumMyElements1, NumMyElements1, 
		  0, IndexBase, Comm, false);

  if (verbose)
    if (ierr==0) cout << "Checked OK\n\n" <<endl;
    else cout << "Error code: "<< ierr << endl;

  assert(ierr==0);

  delete LocalMap1;
  delete LocalMap;
  if (verbose1) {
    if (verbose) cout << "Test ostream << operator" << endl << flush;
    cout << *SmallMap;
    delete SmallMap;
  }

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  return 0;
}

