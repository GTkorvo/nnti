#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Map.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include <vector>

//============================================================================
void Ifpack_BreakForDebugger(Epetra_Comm& Comm)
{
  char hostname[80];
  char buf[80];
  if (Comm.MyPID()  == 0) cout << "Host and Process Ids for tasks" << endl;
  for (int i = 0; i <Comm.NumProc() ; i++) {
    if (i == Comm.MyPID() ) {
      gethostname(hostname, sizeof(hostname));
      sprintf(buf, "Host: %s\tComm.MyPID(): %d\tPID: %d",
	      hostname, Comm.MyPID(), getpid());
      printf("%s\n",buf);
      fflush(stdout);
      sleep(1);
    }
  }
  if(Comm.MyPID() == 0) {
    printf("\n");
    printf("** Pausing to attach debugger...\n");
    printf("** You may now attach debugger to the processes listed above.\n");
    printf( "**\n");
    printf( "** Enter a character to continue > "); fflush(stdout);
    char go;
    scanf("%c",&go);
  }

  Comm.Barrier();

}

//============================================================================
Epetra_CrsMatrix* Ifpack_CreateOverlappingCrsMatrix(Epetra_RowMatrix* Matrix,
						    const int OverlappingLevel)
{

  if (OverlappingLevel == 0) 
    return(0); // All done
  if (Matrix->Comm().NumProc() == 1) 
    return(0); // All done

  Epetra_CrsMatrix* OverlappingMatrix;
  OverlappingMatrix = 0;
  Epetra_Map* OverlappingMap;
  OverlappingMap = (Epetra_Map*)&(Matrix->RowMatrixRowMap());

  Epetra_RowMatrix* OldMatrix;
  const Epetra_Map* DomainMap = &(Matrix->OperatorDomainMap());
  const Epetra_Map* RangeMap = &(Matrix->OperatorRangeMap());

  for (int level = 1; level <= OverlappingLevel ; ++level) {

    if (OverlappingMatrix)
      OldMatrix = OverlappingMatrix;
    else
      OldMatrix = Matrix;

    Epetra_Import* OverlappingImporter;
    OverlappingImporter = (Epetra_Import*)OldMatrix->RowMatrixImporter();
    int NumMyElements = OverlappingImporter->TargetMap().NumMyElements();
    int* MyGlobalElements = OverlappingImporter->TargetMap().MyGlobalElements();

    // need to build an Epetra_Map in this way because Epetra_CrsMatrix
    // requires Epetra_Map and not Epetra_BlockMap

    OverlappingMap = new Epetra_Map(-1,NumMyElements,MyGlobalElements,
				    0, Matrix->Comm());

    if (level < OverlappingLevel)
      OverlappingMatrix = new Epetra_CrsMatrix(Copy, *OverlappingMap, 0);
    else
      // On last iteration, we want to filter out all columns except 
      // those that correspond
      // to rows in the graph.  This assures that our matrix is square
      OverlappingMatrix = new Epetra_CrsMatrix(Copy, *OverlappingMap, 
					       *OverlappingMap, 0);

    OverlappingMatrix->Import(*OldMatrix, *OverlappingImporter, Insert);
    if (level < OverlappingLevel) {
      OverlappingMatrix->TransformToLocal(DomainMap, RangeMap);
    }
    else {
      OverlappingMatrix->TransformToLocal(DomainMap, RangeMap);
    }

    delete OverlappingMap;

    if (level > 1) {
      delete OldMatrix;
    }
    OverlappingMatrix->FillComplete();

  }

  return(OverlappingMatrix);
}

//============================================================================
Epetra_CrsGraph* Ifpack_CreateOverlappingCrsMatrix(Epetra_CrsGraph* Graph,
						   const int OverlappingLevel)
{

  if (OverlappingLevel == 0) 
    return(0); // All done
  if (Graph->Comm().NumProc() == 1) 
    return(0); // All done

  Epetra_CrsGraph* OverlappingGraph;
  Epetra_BlockMap* OverlappingMap;
  OverlappingGraph = const_cast<Epetra_CrsGraph*>(Graph);
  OverlappingMap = const_cast<Epetra_BlockMap*>(&(Graph->RowMap()));

  Epetra_CrsGraph* OldGraph;
  Epetra_BlockMap* OldMap;
  const Epetra_BlockMap* DomainMap = &(Graph->DomainMap());
  const Epetra_BlockMap* RangeMap = &(Graph->RangeMap());

  for (int level = 1; level <= OverlappingLevel ; ++level) {

    OldGraph = OverlappingGraph;
    OldMap = OverlappingMap;

    Epetra_Import* OverlappingImporter;
    OverlappingImporter = const_cast<Epetra_Import*>(OldGraph->Importer());
    OverlappingMap = new Epetra_BlockMap(OverlappingImporter->TargetMap());

    if (level < OverlappingLevel)
      OverlappingGraph = new Epetra_CrsGraph(Copy, *OverlappingMap, 0);
    else
      // On last iteration, we want to filter out all columns except 
      // those that correspond
      // to rows in the graph.  This assures that our matrix is square
      OverlappingGraph = new Epetra_CrsGraph(Copy, *OverlappingMap, 
					  *OverlappingMap, 0);

    OverlappingGraph->Import(*OldGraph, *OverlappingImporter, Insert);
    if (level < OverlappingLevel) 
      OverlappingGraph->TransformToLocal(DomainMap, RangeMap);
    else {
      // Copy last OverlapImporter because we will use it later
      OverlappingImporter = new Epetra_Import(*OverlappingMap, *DomainMap);
      OverlappingGraph->TransformToLocal(DomainMap, RangeMap);
    }

    if (level > 1) {
      delete OldGraph;
      delete OldMap;
    }

    delete OverlappingMap;
    OverlappingGraph->FillComplete();

  }

  return(OverlappingGraph);
}

//============================================================================
string Ifpack_toString(const int& x)
{
  char s[100];
  sprintf(s, "%d", x);
  return string(s);
}

//============================================================================
string Ifpack_toString(const double& x)
{
  char s[100];
  sprintf(s, "%g", x);
  return string(s);
}

//============================================================================
int Ifpack_PrintResidual(char* Label, const Epetra_RowMatrix& A,
                         const Epetra_MultiVector& X, const Epetra_MultiVector&Y)
{
  if (X.Comm().MyPID() == 0) {
    cout << "***** " << Label << endl;
  }
  Ifpack_PrintResidual(0,A,X,Y);

  return(0);
}

//============================================================================
int Ifpack_PrintResidual(const int iter, const Epetra_RowMatrix& A,
                         const Epetra_MultiVector& X, const Epetra_MultiVector&Y)
{
  Epetra_MultiVector RHS(X);
  std::vector<double> Norm2;
  Norm2.resize(X.NumVectors());

  IFPACK_CHK_ERR(A.Multiply(false,X,RHS));
  RHS.Update(1.0, Y, -1.0);

  RHS.Norm2(&Norm2[0]);

  if (X.Comm().MyPID() == 0) {
    cout << "***** iter: " << iter << ":  ||Ax - b||_2 = " 
         << Norm2[0] << endl;
  }

  return(0);
}

//============================================================================
void Ifpack_PrintSparsity(Epetra_RowMatrix& A)
{
  int MaxEntries = A.MaxNumEntries();
  vector<int> Indices(MaxEntries);
  vector<double> Values(MaxEntries);
  vector<bool> FullRow(A.NumMyRows());

  cout << "+-";
  for (int j = 0 ; j < A.NumMyRows() ; ++j)
    cout << '-';
  cout << "-+" << endl;

  for (int i = 0 ; i < A.NumMyRows() ; ++i) {

    int Length;
    A.ExtractMyRowCopy(i,MaxEntries,Length,
                       &Values[0], &Indices[0]);

    for (int j = 0 ; j < A.NumMyRows() ; ++j)
      FullRow[j] = false;

    for (int j = 0 ; j < Length ; ++j) {
      FullRow[Indices[j]] = true;
    }

    cout << "| ";
    for (int j = 0 ; j < A.NumMyRows() ; ++j) {
      if (FullRow[j])
        cout << '*';
      else
        cout << ' ';
    }
    cout << " |" << endl;
  }

  cout << "+-";
  for (int j = 0 ; j < A.NumMyRows() ; ++j)
    cout << '-';
  cout << "-+" << endl << endl;

}


