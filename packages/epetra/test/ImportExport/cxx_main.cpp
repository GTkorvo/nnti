#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#ifndef __cplusplus
#define __cplusplus
#endif
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
 
int main(int argc, char *argv[])
{
  int ierr = 0, i, j;
  bool debug = false;

#ifdef EPETRA_MPI

  // Initialize MPI

  MPI_Init(&argc,&argv);
  int size, rank; // Number of MPI processes, My process ID

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

#else

  int size = 1; // Serial case (not using MPI)
  int rank = 0;
  Epetra_SerialComm Comm;

#endif

  bool verbose = false;

  // Check if we should print results to standard out
  if (argc>1) if (argv[1][0]=='-' && argv[1][1]=='v') verbose = true;




  //char tmp;
  //if (rank==0) cout << "Press any key to continue..."<< endl;
  //if (rank==0) cin >> tmp;
  //Comm.Barrier();

  int MyPID = Comm.MyPID();
  int NumProc = Comm.NumProc();
  if (verbose) cout << "Processor "<<MyPID<<" of "<< NumProc
              << " is alive."<<endl;

  bool verbose1 = verbose;

  // Redefine verbose to only print on PE 0
  if (verbose && rank!=0) verbose = false;

  int NumMyEquations = 10000;
  int NumGlobalEquations = NumMyEquations*NumProc+EPETRA_MIN(NumProc,3);
  if (MyPID < 3) NumMyEquations++;
  int IndexBase = 0;
  bool DistributedGlobal = (NumGlobalEquations>NumMyEquations);

  // Construct a Source Map that puts approximately the same Number of equations on each processor in 
  // uniform global ordering

  Epetra_Map& SourceMap = *new Epetra_Map(NumGlobalEquations, NumMyEquations, 0, Comm);
  
  // Get update list and number of local equations from newly created Map
  int NumMyElements = SourceMap.NumMyElements();
  int * SourceMyGlobalElements = new int[NumMyElements];
  SourceMap.MyGlobalElements(SourceMyGlobalElements);


  // Construct a Target Map that will contain:
  //  some unchanged elements (relative to the soure map),
  //  some permuted elements
  //  some off-processor elements
  Epetra_Vector & RandVec = *new Epetra_Vector(SourceMap);
  RandVec.Random(); // This creates a vector of random numbers between negative one and one.

  int *TargetMyGlobalElements = new int[NumMyElements];

  for (i=0; i< NumMyEquations/2; i++) TargetMyGlobalElements[i] = i; // Half will be the same...
  for (i=NumMyEquations/2; i<NumMyEquations; i++) {
    int index = abs((int)(((double) (NumGlobalEquations-1) ) * RandVec[i]));
    TargetMyGlobalElements[i] = EPETRA_MIN(NumGlobalEquations-1,EPETRA_MAX(0,index));
  }

  int NumSameIDs = 0;
  int NumPermutedIDs = 0;
  int NumRemoteIDs = 0;
  bool StillContiguous = true;
  for (i=0; i < NumMyEquations; i++) {
    if (SourceMyGlobalElements[i]==TargetMyGlobalElements[i] && StillContiguous)
      NumSameIDs++;
    else if (SourceMap.MyGID(TargetMyGlobalElements[i])) {
      StillContiguous = false;
      NumPermutedIDs++;
    }
    else {
      StillContiguous = false;
      NumRemoteIDs++;
    }
  }
  assert(NumMyEquations==NumSameIDs+NumPermutedIDs+NumRemoteIDs);

  Epetra_Map & TargetMap = *new Epetra_Map(-1, NumMyElements, TargetMyGlobalElements, 0, Comm);

  // Create a multivector whose elements are GlobalID * (column number +1)

  int NumVectors = 3;
  Epetra_MultiVector & SourceMultiVector = *new Epetra_MultiVector(SourceMap, NumVectors);
  for (j=0; j < NumVectors; j++)
    for (i=0; i < NumMyElements; i++)
      SourceMultiVector[j][i] = (double) SourceMyGlobalElements[i]*(j+1);

  // Create a target multivector that we will fill using an Import

  Epetra_MultiVector & TargetMultiVector = *new Epetra_MultiVector(TargetMap, NumVectors);

  Epetra_Import & Importer = *new Epetra_Import(TargetMap, SourceMap);

  assert(TargetMultiVector.Import(SourceMultiVector, Importer, Insert)==0);

  // Test Target against expected values

  for (j=0; j < NumVectors; j++)
    for (i=0; i < NumMyElements; i++) {
      if (TargetMultiVector[j][i]!= (double) TargetMyGlobalElements[i]*(j+1))
	cout << "TargetMultiVector["<<i<<"]["<<j<<"] = " << TargetMultiVector[j][i] 
	     <<  "  TargetMyGlobalElements[i]*(j+1) = " <<  TargetMyGlobalElements[i]*(j+1) << endl;
      assert(TargetMultiVector[j][i]== (double) TargetMyGlobalElements[i]*(j+1));
    }

  if (verbose) cout << "MultiVector Import using Importer Check OK" << endl << endl;


  //////////////////////////////////////////////////////////////////////////////

  // Now use Importer to do an export

  Epetra_Vector & TargetVector = *new  Epetra_Vector(SourceMap);
  Epetra_Vector & ExpectedTarget = *new  Epetra_Vector(SourceMap);
  Epetra_Vector & SourceVector = *new  Epetra_Vector(TargetMap);

  NumSameIDs = Importer.NumSameIDs();
  int NumPermuteIDs = Importer.NumPermuteIDs();
  int NumExportIDs = Importer.NumExportIDs();
  int *PermuteToLIDs = Importer.PermuteToLIDs();
  int *PermuteFromLIDs = Importer.PermuteFromLIDs();
  int *ExportLIDs = Importer.ExportLIDs();
  int *ExportPIDs = Importer.ExportPIDs();

  for (i=0; i < NumSameIDs; i++) ExpectedTarget[i] = (double) (MyPID+1);
  for (i=0; i < NumPermuteIDs; i++) ExpectedTarget[PermuteFromLIDs[i]] = 
				      (double) (MyPID+1);
  for (i=0; i < NumExportIDs; i++) ExpectedTarget[ExportLIDs[i]] += 
				     (double) (ExportPIDs[i]+1);

  for (i=0; i < NumMyElements; i++) SourceVector[i] =  (double) (MyPID+1);

  assert(TargetVector.Export(SourceVector, Importer, Add)==0);

    for (i=0; i < NumMyElements; i++) {
      if (TargetVector[i]!= ExpectedTarget[i])
	cout <<  "     TargetVector["<<i<<"] = " << TargetVector[i] 
	     <<  "   ExpectedTarget["<<i<<"] = " <<  ExpectedTarget[i] << " on PE " << MyPID << endl;
      assert(TargetVector[i]== ExpectedTarget[i]);
    }

  if (verbose) cout << "Vector Export using Importer Check OK" << endl << endl;



  //////////////////////////////////////////////////////////////////////////////////////////
  //  Build a tridiagonal system two ways:
  //  1) From "standard" matrix view where equations are uniquely owned.
  //  2) From 1D PDE view where nodes (equations) between processors are shared and partial contributions are done
  //     in parallel, then merged together at the end of the construction process.
  //
  //////////////////////////////////////////////////////////////////////////////////////////



  // Construct a Standard Map that puts approximately the same number of equations on each processor in 
  // uniform global ordering

  Epetra_Map& StandardMap = *new Epetra_Map(NumGlobalEquations, NumMyEquations, 0, Comm);
  
  // Get update list and number of local equations from newly created Map
  NumMyElements = StandardMap.NumMyElements();
  int * StandardMyGlobalElements = new int[NumMyElements];
  StandardMap.MyGlobalElements(StandardMyGlobalElements);


  // Create a standard Epetra_CrsGraph

  Epetra_CrsGraph& StandardGraph = *new Epetra_CrsGraph(Copy, StandardMap, 3);
  assert(!StandardGraph.IndicesAreGlobal());
  assert(!StandardGraph.IndicesAreLocal());
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  int *Indices = new int[2];
  int NumEntries;
  
  for (i=0; i<NumMyEquations; i++)
    {
    if (StandardMyGlobalElements[i]==0)
      {
	Indices[0] = 1;
	NumEntries = 1;
      }
    else if (StandardMyGlobalElements[i] == NumGlobalEquations-1)
      {
	Indices[0] = NumGlobalEquations-2;
	NumEntries = 1;
      }
    else
      {
	Indices[0] = StandardMyGlobalElements[i]-1;
	Indices[1] = StandardMyGlobalElements[i]+1;
	NumEntries = 2;
      }
     assert(StandardGraph.InsertGlobalIndices(StandardMyGlobalElements[i], NumEntries, Indices)==0);
     assert(StandardGraph.InsertGlobalIndices(StandardMyGlobalElements[i], 1, StandardMyGlobalElements+i)==0); // Put in the diagonal entry
    }
  
  // Finish up
  assert(StandardGraph.IndicesAreGlobal());
  assert(StandardGraph.TransformToLocal()==0);
  assert(StandardGraph.IndicesAreLocal());
  assert(!StandardGraph.StorageOptimized());
  StandardGraph.OptimizeStorage();
  assert(StandardGraph.StorageOptimized());
  assert(!StandardGraph.UpperTriangular());
  assert(!StandardGraph.LowerTriangular());


  // Create Epetra_CrsMatrix using the just-built graph

  Epetra_CrsMatrix& StandardMatrix = *new Epetra_CrsMatrix(Copy, StandardGraph);
  assert(!StandardMatrix.IndicesAreGlobal());
  assert(StandardMatrix.IndicesAreLocal());
  
  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1


  double *Values = new double[2];
  Values[0] = -1.0; Values[1] = -1.0;
  double two = 2.0;
  
  for (i=0; i<NumMyEquations; i++)
    {
    if (StandardMyGlobalElements[i]==0)
      {
	Indices[0] = 1;
	NumEntries = 1;
      }
    else if (StandardMyGlobalElements[i] == NumGlobalEquations-1)
      {
	Indices[0] = NumGlobalEquations-2;
	NumEntries = 1;
      }
    else
      {
	Indices[0] = StandardMyGlobalElements[i]-1;
	Indices[1] = StandardMyGlobalElements[i]+1;
	NumEntries = 2;
      }
     assert(StandardMatrix.ReplaceGlobalValues(StandardMyGlobalElements[i], NumEntries, Values, Indices)==0);
     assert(StandardMatrix.ReplaceGlobalValues(StandardMyGlobalElements[i], 1, &two, StandardMyGlobalElements+i)==0); // Put in the diagonal entry
    }
  
  // Finish up
  assert(StandardMatrix.IndicesAreLocal());
  assert(StandardMatrix.TransformToLocal()==0);
  assert(StandardMatrix.IndicesAreLocal());
  assert(StandardMatrix.StorageOptimized());
  assert(!StandardMatrix.UpperTriangular());
  assert(!StandardMatrix.LowerTriangular());

  // Construct an Overlapped Map of StandardMap that include the endpoints from two neighboring processors.

  int OverlapNumMyElements;
  int OverlapMinMyGID;

  OverlapNumMyElements = NumMyElements + 1;
  if (MyPID==0) OverlapNumMyElements--;

  if (MyPID==0) OverlapMinMyGID = StandardMap.MinMyGID();
  else OverlapMinMyGID = StandardMap.MinMyGID()-1;

  int * OverlapMyGlobalElements = new int[OverlapNumMyElements];

  for (i=0; i< OverlapNumMyElements; i++) OverlapMyGlobalElements[i] = OverlapMinMyGID + i;

  Epetra_Map& OverlapMap = *new Epetra_Map(-1, OverlapNumMyElements, OverlapMyGlobalElements, 0, Comm);

  // Create the Overlap Epetra_Matrix

  Epetra_CrsMatrix& OverlapMatrix = *new Epetra_CrsMatrix(Copy, OverlapMap, 4);
  assert(!OverlapMatrix.IndicesAreGlobal());
  assert(!OverlapMatrix.IndicesAreLocal());
  
  // Add  matrix element one cell at a time.
  // Each cell does an incoming and outgoing flux calculation


  double pos_one = 1.0;
  double neg_one = -1.0;

  for (i=0; i<OverlapNumMyElements; i++)
    {
      int node_left = OverlapMyGlobalElements[i]-1;
      int node_center = node_left + 1;
      int node_right = node_left + 2;
      if (i>0) {
	if (node_left>-1)
	  assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &neg_one, &node_left)==0);
	assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0);
      }
      if (i<OverlapNumMyElements-1) {
	assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0);
	if (node_right<NumGlobalEquations) 
	  assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &neg_one, &node_right)==0);
      }
    }

  // Handle endpoints
  if (MyPID==0) {
    int node_center = 0;
    assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0);
  }
  if (MyPID==NumProc-1) {
    int node_center = OverlapMyGlobalElements[OverlapNumMyElements-1];
    assert(OverlapMatrix.InsertGlobalValues(node_center, 1, &pos_one, &node_center)==0);
  }
    
  assert(OverlapMatrix.TransformToLocal()==0);

  // Make a gathered matrix from OverlapMatrix.  It should be identical to StandardMatrix

  Epetra_CrsMatrix& GatheredMatrix = *new Epetra_CrsMatrix(Copy, StandardGraph);
  Epetra_Export & Exporter = *new Epetra_Export(OverlapMap, StandardMap);
  assert(GatheredMatrix.Export(OverlapMatrix, Exporter, Add)==0);
  assert(GatheredMatrix.TransformToLocal()==0);

  // Check if entries of StandardMatrix and GatheredMatrix are identical

  int StandardNumEntries, GatheredNumEntries;
  int * StandardIndices, * GatheredIndices;
  double * StandardValues, * GatheredValues;

  int StandardNumMyNonzeros = StandardMatrix.NumMyNonzeros();
  int GatheredNumMyNonzeros = GatheredMatrix.NumMyNonzeros();
  assert(StandardNumMyNonzeros==GatheredNumMyNonzeros);

  int StandardNumMyRows = StandardMatrix.NumMyRows();
  int GatheredNumMyRows = GatheredMatrix.NumMyRows();
  assert(StandardNumMyRows==GatheredNumMyRows);

  for (i=0; i< StandardNumMyRows; i++)
    {
      assert(StandardMatrix.ExtractMyRowView(i, StandardNumEntries, StandardValues, StandardIndices)==0);
      assert(GatheredMatrix.ExtractMyRowView(i, GatheredNumEntries, GatheredValues, GatheredIndices)==0);
      assert(StandardNumEntries==GatheredNumEntries);
      for (j=0; j < StandardNumEntries; j++) {
	//if (StandardIndices[j]!=GatheredIndices[j])
	// cout << "MyPID = " << MyPID << " i = " << i << "   StandardIndices[" << j << "] = " << StandardIndices[j] 
	//      << "   GatheredIndices[" << j << "] = " << GatheredIndices[j] << endl;
	//if (StandardValues[j]!=GatheredValues[j])
	//cout << "MyPID = " << MyPID << " i = " << i << "    StandardValues[" << j << "] = " <<  StandardValues[j] 
	//     << "    GatheredValues[" << j << "] = " <<  GatheredValues[j] << endl;
	assert(StandardIndices[j]==GatheredIndices[j]);
	assert(StandardValues[j]==GatheredValues[j]);
      }
    }

  if (verbose) cout << "Matrix Export Check OK" << endl;
  // Release all objects

  delete &SourceVector;
  delete &TargetVector;
  delete &ExpectedTarget;


  delete &Importer;
  delete &SourceMap;
  delete &TargetMap;

  delete [] SourceMyGlobalElements;
  delete [] TargetMyGlobalElements;

  delete &SourceMultiVector;
  delete &TargetMultiVector;
  delete &RandVec;

  delete &Exporter;
  delete &GatheredMatrix;
  delete &OverlapMatrix;
  delete &OverlapMap;
  delete [] OverlapMyGlobalElements;

  delete &StandardMatrix;
  delete &StandardGraph;
  delete &StandardMap;
  delete [] StandardMyGlobalElements;

  delete [] Values;
  delete [] Indices;

#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

/* end main
*/
return ierr ;
}
