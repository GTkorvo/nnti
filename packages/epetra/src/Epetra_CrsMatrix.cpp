
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int *NumEntriesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyRows_(RowMap.NumMyPoints()),
    CV_(CV)
{
  Graph_ = new Epetra_CrsGraph(CV, RowMap, NumEntriesPerRow);
  InitializeDefaults();
  assert(Allocate()==0);
}

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, int NumEntriesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyRows_(RowMap.NumMyPoints()),
    CV_(CV)
{
  Graph_ = new Epetra_CrsGraph(CV, RowMap, NumEntriesPerRow);
  InitializeDefaults();
  assert(Allocate()==0);
}
//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, 
				   const Epetra_Map& ColMap, int *NumEntriesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyRows_(RowMap.NumMyPoints()),
    CV_(CV)
{
  Graph_ = new Epetra_CrsGraph(CV, RowMap, ColMap, NumEntriesPerRow);
  InitializeDefaults();
  assert(Allocate()==0);
}

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, 
				   const Epetra_Map& ColMap, int NumEntriesPerRow) 
  : Epetra_DistObject(RowMap, "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(false),
    StaticGraph_(false),
    NumMyRows_(RowMap.NumMyPoints()),
    CV_(CV)
{
  Graph_ = new Epetra_CrsGraph(CV, RowMap, ColMap,  NumEntriesPerRow);
  InitializeDefaults();
  assert(Allocate()==0);
}
//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_CrsGraph & Graph) 
  : Epetra_DistObject(Graph.Map(), "Epetra::CrsMatrix"),
    Epetra_CompObject(),
    Epetra_BLAS(),
    Graph_((Epetra_CrsGraph*) &Graph),
    Allocated_(false),
    StaticGraph_(true),
    NumMyRows_(Graph.NumMyRows_),
    CV_(CV)
{
  InitializeDefaults();
  assert(Allocate()==0);
}

//==============================================================================
Epetra_CrsMatrix::Epetra_CrsMatrix(const Epetra_CrsMatrix & Matrix) 
  : Epetra_DistObject(Matrix),
    Epetra_CompObject(Matrix),
    Epetra_BLAS(),
    Graph_(0),
    Allocated_(Matrix.Allocated_),
    StaticGraph_(false),
    UseTranspose_(Matrix.UseTranspose_),
    Values_(0),
    All_Values_(0),
    NormInf_(-1.0),
    NormOne_(-1.0),
    NumMyRows_(Matrix.NumMyRows_),
    ImportVector_(0),
    CV_(Copy)
{
  Graph_ = new Epetra_CrsGraph(Matrix.Graph());
  assert(Allocate()==0);
  for (int i=0; i<NumMyRows_; i++) {
    int NumEntries = NumEntriesPerRow_[i];
    for (int j=0; j< NumEntries; j++) Values_[i][j] = Matrix.Values_[i][j];
  }
}

//==============================================================================
void Epetra_CrsMatrix::InitializeDefaults() { // Initialize all attributes that have trivial default values

  UseTranspose_ = false;
  Values_ = 0;
  All_Values_ = 0;
  NormInf_ = -1.0;
  NormOne_ = -1.0;
  ImportVector_ = 0;

  NumEntriesPerRow_  = 0;
  NumAllocatedEntriesPerRow_ = 0;
  Indices_ = 0;

  return;
}

//==============================================================================
int Epetra_CrsMatrix::Allocate() {

  int i, j;
  
  // Set direct access pointers to graph info (needed for speed)
  NumEntriesPerRow_ = Graph_->NumIndicesPerRow();
  NumAllocatedEntriesPerRow_ = Graph_->NumAllocatedIndicesPerRow();
  Indices_ = Graph_->Indices();

  // Allocate Values array
  Values_ = new double*[NumMyRows_];

  // Allocate and initialize entries if we are copying data
  if (CV_==Copy) {
    for (i=0; i<NumMyRows_; i++) {
      int NumAllocatedEntries = NumAllocatedEntriesPerRow_[i];

      if (NumAllocatedEntries > 0) Values_[i] = new double[NumAllocatedEntries];
      else Values_[i] = 0;

      for (j=0; j< NumAllocatedEntries; j++) Values_[i][j] = 0.0; // Fill values with zero
    }
  }	 
  else {
    for (i=0; i<NumMyRows_; i++) {
      Values_[i] = 0;
    }
  }
    SetAllocated(true);
    return(0);
}
//==============================================================================
Epetra_CrsMatrix::~Epetra_CrsMatrix(){

  int i;

  if (CV_==Copy) {
    if (All_Values_!=0) delete [] All_Values_;
    else for (i=0; i<NumMyRows_; i++) if (NumAllocatedEntriesPerRow_[i] >0) delete [] Values_[i];
  }

  if (ImportVector_!=0) delete ImportVector_;
  ImportVector_=0;
    
    
  delete [] Values_;
  if (!StaticGraph()) delete Graph_; // We created the graph, so must delete it.

  NumMyRows_ = 0;
  
  Allocated_ = false;
}

//==============================================================================
int Epetra_CrsMatrix::PutScalar(double ScalarConstant) 
{
  for (int i=0; i<NumMyRows_; i++) {
    int NumEntries = NumEntriesPerRow_[i];
    for (int j=0; j< NumEntries; j++) Values_[i][j] = ScalarConstant;
  }
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::InsertGlobalValues(int Row, int NumEntries, double * Values, int *Indices) {

  if (IndicesAreLocal()) EPETRA_CHK_ERR(-2); // Cannot insert global values into local graph
  if (IndicesAreContiguous()) EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and newed
  Graph_->SetIndicesAreGlobal(true);
  Row = Graph_->LRID(Row); // Find local row number for this global row index

  EPETRA_CHK_ERR( InsertValues(Row, NumEntries, Values, Indices) );

//  if( !StaticGraph() ) EPETRA_CHK_ERR( Graph_->InsertGlobalIndices( Row, NumEntries, Indices ) );

  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::InsertMyValues(int Row, int NumEntries, double * Values, int *Indices) {

  if (IndicesAreGlobal()) EPETRA_CHK_ERR(-2); // Cannot insert global values into filled graph
  if (IndicesAreContiguous()) EPETRA_CHK_ERR(-3); // Indices cannot be individually deleted and new
  Graph_->SetIndicesAreLocal(true);

  EPETRA_CHK_ERR( InsertValues(Row, NumEntries, Values, Indices) );

//  if( !StaticGraph() ) EPETRA_CHK_ERR( Graph_->InsertMyIndices( Row, NumEntries, Indices ) );

  return(0);

}

//==========================================================================
int Epetra_CrsMatrix::InsertValues(int Row, int NumEntries, double * Values, int *Indices) {

  int j;
  double * tmp_Values;
  int ierr = 0;

  if (Row < 0 || Row >= NumMyRows_) EPETRA_CHK_ERR(-1); // Not in Row range
    
  if (CV_==View) {

    //test indices in static graph
    if( StaticGraph() )
    {
      int testNumEntries;
      int * testIndices;
      int testRow = Row;
      if( IndicesAreGlobal() ) testRow = Graph_->LRID( Row );
      EPETRA_CHK_ERR( Graph_->ExtractMyRowView( testRow, testNumEntries, testIndices ) );

      bool match = true;
      if( NumEntries != testNumEntries ) match = false;
      for( int i = 0; i < NumEntries; ++i ) match = match && (Indices[i]==testIndices[i]);

      if( !match ) ierr = -3;
    }

    if (Values_[Row]!=0) ierr = 2; // This row has be defined already.  Issue warning.
    Values_[Row] = Values;

  }
  else {
    
    if (StaticGraph()) EPETRA_CHK_ERR(-2); // If the matrix graph is fully constructed, we cannot insert new values

    int start = NumEntriesPerRow_[Row];
    int stop = start + NumEntries;
    int NumAllocatedEntries = NumAllocatedEntriesPerRow_[Row];
    if (stop > NumAllocatedEntries){
      if (NumAllocatedEntries==0) Values_[Row] = new double[NumEntries]; // Row was never allocated, so do it
      else {
	ierr = 1; // Out of room.  Must delete and allocate more space...
	tmp_Values = new double[stop];
	for (j=0; j< start; j++) tmp_Values[j] = Values_[Row][j]; // Copy existing entries
	delete [] Values_[Row]; // Delete old storage
	Values_[Row] = tmp_Values; // Set pointer to new storage
      }
    }
        
    for (j=start; j<stop; j++) Values_[Row][j] = Values[j-start];

  }

  if( !StaticGraph() ) Graph_->InsertIndices( Row, NumEntries, Indices );


  EPETRA_CHK_ERR(ierr);
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  return(0);

}

//==========================================================================
int Epetra_CrsMatrix::ReplaceGlobalValues(int Row, int NumEntries, double * Values, int *Indices) {

  int j;
  int ierr = 0;
  int Loc;

  if (CV_==View) EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

   Row = Graph_->LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyRows_) EPETRA_CHK_ERR(-1); // Not in Row range
    
  for (j=0; j<NumEntries; j++) {
    int Index = Indices[j];
    if (Graph_->FindGlobalIndexLoc(Row,Index,j,Loc)) Values_[Row][Loc] = Values[j];
    else EPETRA_CHK_ERR(-2); // Value not found
  }

  EPETRA_CHK_ERR(ierr);
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::ReplaceMyValues(int Row, int NumEntries, double * Values, int *Indices) {

  if (!IndicesAreLocal()) EPETRA_CHK_ERR(-4); // Indices must be local.

  int j;
  int ierr = 0;
  int Loc;

  if (CV_==View) EPETRA_CHK_ERR(-3); // This is a view only.  Cannot remove entries.

    
  if (Row < 0 || Row >= NumMyRows_) EPETRA_CHK_ERR(-1); // Not in Row range
    
  for (j=0; j<NumEntries; j++) {
    int Index = Indices[j];
    if (Graph_->FindMyIndexLoc(Row,Index,j,Loc)) Values_[Row][Loc] = Values[j];
    else EPETRA_CHK_ERR(-2); // Value not found
  }

  EPETRA_CHK_ERR(ierr);
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::SumIntoGlobalValues(int Row, int NumEntries, double * Values, int *Indices) {

  int j;
  int ierr = 0;
  int Loc;

  Row = Graph_->LRID(Row); // Normalize row range
    
  if (Row < 0 || Row >= NumMyRows_) EPETRA_CHK_ERR(-1); // Not in Row range
    
  for (j=0; j<NumEntries; j++) {
    int Index = Indices[j];
    if (Graph_->FindGlobalIndexLoc(Row,Index,j,Loc)) Values_[Row][Loc] += Values[j];
    else EPETRA_CHK_ERR(-2); // Value not found
  }

  EPETRA_CHK_ERR(ierr);
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::SumIntoMyValues(int Row, int NumEntries, double * Values, int *Indices) {

  if (!IndicesAreLocal()) EPETRA_CHK_ERR(-4); // Indices must be local.

  int j;
  int ierr = 0;
  int Loc;

  if (Row < 0 || Row >= NumMyRows_) EPETRA_CHK_ERR(-1); // Not in Row range
    
  for (j=0; j<NumEntries; j++) {
    int Index = Indices[j];
    if (Graph_->FindMyIndexLoc(Row,Index,j,Loc)) Values_[Row][Loc] += Values[j];
    else EPETRA_CHK_ERR(-2); // Value not found
  }

  EPETRA_CHK_ERR(ierr);
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::TransformToLocal() {
  EPETRA_CHK_ERR(TransformToLocal((Epetra_Map *) (&RowMap()),(Epetra_Map *) (&RowMap())));
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::TransformToLocal(Epetra_Map *DomainMap, Epetra_Map *RangeMap) {
  
  if (!StaticGraph()) EPETRA_CHK_ERR(Graph_->MakeIndicesLocal(*DomainMap, *RangeMap));
  SortEntries();  // Sort column entries from smallest to largest
  MergeRedundantEntries(); // Get rid of any redundant index values
  if (!StaticGraph()) EPETRA_CHK_ERR(Graph_->TransformToLocal(DomainMap, RangeMap));


  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::SortEntries() {

  if (!IndicesAreLocal()) EPETRA_CHK_ERR(-1);
  if (Sorted()) return(0);

  // For each row, sort column entries from smallest to largest.
  // Use shell sort. Stable sort so it is fast if indices are already sorted.

  
  for (int i=0; i<NumMyRows_; i++){

    double * Values = Values_[i];
    int NumEntries = NumEntriesPerRow_[i];
    int * Indices = Indices_[i];

    int n = NumEntries;
    int m = n/2;
    
    while (m > 0) {
      int max = n - m;
      for (int j=0; j<max; j++)
        {
	  for (int k=j; k>=0; k-=m)
            {
	      if (Indices[k+m] >= Indices[k])
		break;
	      double dtemp = Values[k+m];
	      Values[k+m] = Values[k];
	      Values[k] = dtemp;
	      int itemp = Indices[k+m];
	      Indices[k+m] = Indices[k];
	      Indices[k] = itemp;
            }
        }
      m = m/2;
    }
  }
  Graph_->SetSorted(true); // This also sorted the graph
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::MergeRedundantEntries() {

  int i;

  if (NoRedundancies()) return(0);
  if (!Sorted()) EPETRA_CHK_ERR(-1);  // Must have sorted entries

  // For each row, remove column indices that are repeated.
  // Also, determine if matrix is upper or lower triangular or has no diagonal (Done in graph)
  // Note:  This function assumes that SortEntries was already called.

  for (i=0; i<NumMyRows_; i++) {
    int NumEntries = NumEntriesPerRow_[i];
    if (NumEntries>1) {
      double * const Values = Values_[i];
      int * const Indices = Indices_[i];
			
			int curEntry =0;
			double curValue = Values[0];
			for (int k=1; k<NumEntries; k++) {
				if (Indices[k]==Indices[k-1]) curValue += Values[k];
				else {
					Values[curEntry++] = curValue;
					curValue = Values[k];
				}
      }
			Values[curEntry] = curValue;

    }
  }

  EPETRA_CHK_ERR(Graph_->RemoveRedundantIndices()); // Remove redundant indices and then return
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::OptimizeStorage() {

  int i, j;

  if (StorageOptimized()) return(0); // Have we been here before?

  bool Contiguous = true; // Assume contiguous is true
  for (i=1; i<NumMyRows_; i++){
    int NumEntries = NumEntriesPerRow_[i];
    int NumAllocatedEntries = NumAllocatedEntriesPerRow_[i];
      
    // Check if NumEntries is same as NumAllocatedEntries and 
    // check if end of beginning of current row starts immediately after end of previous row.
    if ((NumEntries!=NumAllocatedEntries) || (Values_[i]!=Values_[i-1]+NumEntries)) {
      Contiguous = false;
      break;
    }
  }

  // NOTE:  At the end of the above loop set, there is a possibility that NumEntries and NumAllocatedEntries
  //        for the last row could be different, but I don't think it matters.


  if ((CV_==View) && !Contiguous) EPETRA_CHK_ERR(-1);  // This is user data, it's not contiguous and we can't make it so.

  int ierr = Graph_->OptimizeStorage(); // Make sure graph has optimized storage
  if (ierr) EPETRA_CHK_ERR(ierr);

  if (Contiguous) return(0); // Everything is done.  Return

 // Compute Number of Nonzero entries (Done in FillComplete, but we may not have been there yet.)
  int NumMyNonzeros = Graph_->NumMyNonzeros();

  // Allocate one big integer array for all index values
  All_Values_ = new double[NumMyNonzeros];
  
  // Set Entries_ to point into All_Entries_
  
  double * tmp = All_Values_;
  for (i=0; i<NumMyRows_; i++) {
    int NumEntries = NumEntriesPerRow_[i];
    for (j=0; j<NumEntries; j++) tmp[j] = Values_[i][j];
    if (Values_[i] !=0) delete [] Values_[i];
    Values_[i] = tmp;
    tmp += NumEntries;
  }
  
  
    return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractGlobalRowCopy(int Row, int Length, int & NumEntries, double * Values,
					 int * Indices) const 
{

  int ierr = Graph_->ExtractGlobalRowCopy(Row, Length, NumEntries, Indices);
  if (ierr) EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractGlobalRowCopy(Row, Length, NumEntries, Values));
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, double * Values,
					 int * Indices) const 
{

  int ierr = Graph_->ExtractMyRowCopy(Row, Length, NumEntries, Indices);
  if (ierr) EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractMyRowCopy(Row, Length, NumEntries, Values));
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::NumMyRowEntries(int Row, int & NumEntries) const 
{

  if (!MyLRID(Row)) EPETRA_CHK_ERR(-1); // Not in the range of local rows
  NumEntries = NumMyEntries(Row);
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractGlobalRowCopy(int Row, int Length, int & NumEntries, double * Values) const 
{

  int Row0 = Graph_->RowMap().LID(Row); // Normalize row range

  EPETRA_CHK_ERR(ExtractMyRowCopy(Row0, Length, NumEntries, Values));
  return(0);
}


//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowCopy(int Row, int Length, int & NumEntries, double * Values) const 
{
  int j;

  if (Row < 0 || Row >= NumMyRows_) EPETRA_CHK_ERR(-1); // Not in Row range

  NumEntries = NumEntriesPerRow_[Row];
  if (Length < NumEntries) EPETRA_CHK_ERR(-2); // Not enough space for copy. Needed size is passed back in NumEntries


  for(j=0; j<NumEntries; j++)Values[j] = Values_[Row][j];
  
  return(0);
}


//==============================================================================
int Epetra_CrsMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const {
	
  if (!Filled()) EPETRA_CHK_ERR(-1); // Can't get diagonal unless matrix is filled
  if (!RowMap().SameAs(Diagonal.Map())) EPETRA_CHK_ERR(-2); // Maps must be the same

  int Base = IndexBase();
  for(int i=0; i<NumMyRows_; i++){
    int Row = i + Base;
    int NumEntries = NumEntriesPerRow_[i];
    int * Indices = Indices_[i];
    Diagonal[i] = 0.0;
    for (int j=0; j<NumEntries; j++) {
      int Col = Indices[j];
      if (Row==Col) {
	Diagonal[i] = Values_[i][j];
	break;
      }
    }
  }
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractGlobalRowView(int Row, int & NumEntries, double *& Values, int *& Indices) const 
{

  int ierr = Graph_->ExtractGlobalRowView(Row, NumEntries, Indices);
  if (ierr) EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractGlobalRowView(Row, NumEntries, Values));
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowView(int Row, int & NumEntries, double *& Values, int *& Indices) const 
{

  int ierr = Graph_->ExtractMyRowView(Row, NumEntries, Indices);
  if (ierr) EPETRA_CHK_ERR(ierr);

  EPETRA_CHK_ERR(ExtractMyRowView(Row, NumEntries, Values));
  return(0);
}
//==========================================================================
int Epetra_CrsMatrix::ExtractGlobalRowView(int Row, int & NumEntries, double *& Values) const 
{

  int Row0 = Graph_->RowMap().LID(Row); // Normalize row range

  EPETRA_CHK_ERR(ExtractMyRowView(Row0, NumEntries, Values));
  return(0);
}

//==========================================================================
int Epetra_CrsMatrix::ExtractMyRowView(int Row, int & NumEntries, double *& Values) const 
{

  if (Row < 0 || Row >= NumMyRows_) EPETRA_CHK_ERR(-1); // Not in Row range

  NumEntries = NumEntriesPerRow_[Row];

  Values = Values_[Row];
  
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Multiply(bool TransA, const Epetra_Vector& x, Epetra_Vector& y) const {
//
// This function forms the product y = A * x or y = A' * x
//

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;
  double * xp = (double*)x.Values();
  double *yp = (double*)y.Values();
  int NumMyCols_ = NumMyCols();


  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()!=1) { delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),1); // Create import vector if needed
      ImportVector_->Import(x, *Importer(), Insert);
      xp = (double*)ImportVector_->Values();
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()!=1) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
      yp = (double*)ExportVector_->Values();
    }

    // Do actual computation

    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      double sum = 0.0;
      for (j=0; j < NumEntries; j++) sum += RowValues[j] * xp[RowIndices[j]];

      yp[i] = sum;

    }
    if (Exporter()!=0) y.Export(*ExportVector_, *Exporter(), Add); // Fill y with Values from export vector
  }

  else { // Transpose operation


    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors

    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()!=1) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),1); // Create Export vector if needed
      ExportVector_->Import(x, *Exporter(), Insert);
      xp = (double*)ExportVector_->Values();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()!=1) { delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),1); // Create import vector if needed
      yp = (double*)ImportVector_->Values();
    }

    // Do actual computation

    for (i=0; i < NumMyCols_; i++) yp[i] = 0.0; // Initialize y for transpose multiply
        
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      for (j=0; j < NumEntries; j++) yp[RowIndices[j]] += RowValues[j] * xp[i];
    }
    if (Importer()!=0) y.Export(*ImportVector_, *Importer(), Add); // Fill y with Values from export vector
  }

    UpdateFlops(2*NumGlobalNonzeros());
    return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
//
// This function forms the product Y = A * Y or Y = A' * X
//
  if (X.NumVectors()==1 && Y.NumVectors()==1) {
    double * xp = (double *) X[0];
    double * yp = (double *) Y[0];
    Epetra_Vector x(View, X.Map(), xp);
    Epetra_Vector y(View, Y.Map(), yp);
    return(Multiply(TransA, x, y));
  }
  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  int i, j, k;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;

  double **Xp = (double**)X.Pointers();
  double **Yp = (double**)Y.Pointers();

  int NumVectors = X.NumVectors();
  int NumMyCols_ = NumMyCols();


  // Need to better manage the Import and Export vectors:
  // - Need accessor functions
  // - Need to make the NumVector match (use a View to do this)
  // - Need to look at RightScale and ColSum routines too.

  if (!TransA) {

    // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()!=NumVectors) { delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors); // Create import vector if needed
      ImportVector_->Import(X, *Importer(), Insert);
      Xp = (double**)ImportVector_->Pointers();
    }

    // If we have a non-trivial exporter, we must export elements that are permuted or belong to other processors
    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()!=NumVectors) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
      Yp = (double**)ExportVector_->Pointers();
    }

    // Do actual computation

    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      for (k=0; k<NumVectors; k++) {
	double sum = 0.0;
	for (j=0; j < NumEntries; j++) sum += RowValues[j] * Xp[k][RowIndices[j]];
	Yp[k][i] = sum;
      }
    }
    if (Exporter()!=0) Y.Export(*ExportVector_, *Exporter(), Add); // Fill Y with Values from export vector
  }
  else { // Transpose operation


    // If we have a non-trivial exporter, we must import elements that are permuted or are on other processors

    if (Exporter()!=0) {
      if (ExportVector_!=0) {
	if (ExportVector_->NumVectors()!=NumVectors) { delete ExportVector_; ExportVector_= 0;}
      }
      if (ExportVector_==0) ExportVector_ = new Epetra_MultiVector(RowMap(),NumVectors); // Create Export vector if needed
      ExportVector_->Import(X, *Exporter(), Insert);
      Xp = (double**)ExportVector_->Pointers();
    }

    // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
    if (Importer()!=0) {
      if (ImportVector_!=0) {
	if (ImportVector_->NumVectors()!=NumVectors) { delete ImportVector_; ImportVector_= 0;}
      }
      if (ImportVector_==0) ImportVector_ = new Epetra_MultiVector(ColMap(),NumVectors); // Create import vector if needed
      Yp = (double**)ImportVector_->Pointers();
    }

    // Do actual computation



        for (k=0; k<NumVectors; k++) 
	  for (i=0; i < NumMyCols_; i++) Yp[k][i] = 0.0; // Initialize y for transpose multiply
    
    for (i=0; i < NumMyRows_; i++) {
      int      NumEntries = *NumEntriesPerRow++;
      int *    RowIndices = *Indices++;
      double * RowValues  = *Values++;
      for (k=0; k<NumVectors; k++) {
	for (j=0; j < NumEntries; j++) Yp[k][RowIndices[j]] += RowValues[j] * Xp[k][i];
      }
    }
    if (Importer()!=0) Y.Export(*ImportVector_, *Importer(), Add); // Fill Y with Values from export vector
  }

  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_Vector& x, Epetra_Vector& y) const {
//
// This function find y such that Ly = x or Uy = x or the transpose cases.
//

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  if ((Upper) && (!UpperTriangular())) EPETRA_CHK_ERR(-2);
  if ((!Upper) && (!LowerTriangular())) EPETRA_CHK_ERR(-3);
  if ((!UnitDiagonal) && (NoDiagonal())) EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_)) EPETRA_CHK_ERR(-5); // Need each row to have a diagonal
      

  int i, j, j0;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;
  int NumMyCols_ = NumMyCols();

  // If upper, point to last row
  if ((Upper && !Trans) || (!Upper && Trans)) {
    NumEntriesPerRow += NumMyRows_-1;
    Indices += NumMyRows_-1;
    Values += NumMyRows_-1;
  }
    
    double *xp = (double*)x.Values();
    double *yp = (double*)y.Values();

  if (!Trans) {

    if (Upper) {

      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=NumMyRows_-1; i >=0; i--) {
	int      NumEntries = *NumEntriesPerRow--;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	double sum = 0.0;
	for (j=j0; j < NumEntries; j++) sum += RowValues[j] * yp[RowIndices[j]];
	
	if (UnitDiagonal) yp[i] = xp[i] - sum;
	else yp[i] = (xp[i] - sum)/RowValues[0];

      }
    }
    else {
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++ - j0;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	double sum = 0.0;
	for (j=0; j < NumEntries; j++) sum += RowValues[j] * yp[RowIndices[j]];
	
	if (UnitDiagonal) yp[i] = xp[i] - sum;
	else yp[i] = (xp[i] - sum)/RowValues[NumEntries];

      }
    }
  }

  // ***********  Transpose case *******************************

  else {

    if (xp!=yp) for (i=0; i < NumMyCols_; i++) yp[i] = xp[i]; // Initialize y for transpose solve
    
    if (Upper) {

      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
    
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	if (!UnitDiagonal) yp[i] = yp[i]/RowValues[0];
	for (j=j0; j < NumEntries; j++) yp[RowIndices[j]] -= RowValues[j] * yp[i];
      }
    }
    else {

      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
    
      for (i=NumMyRows_-1; i >= 0; i--) {
	int      NumEntries = *NumEntriesPerRow-- - j0;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	if (!UnitDiagonal)  yp[i] = yp[i]/RowValues[NumEntries];
	for (j=0; j < NumEntries; j++) yp[RowIndices[j]] -= RowValues[j] * yp[i];
      }
    }

  }
    UpdateFlops(2*NumGlobalNonzeros());
    return(0);
}

//=============================================================================
int Epetra_CrsMatrix::Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
//
// This function find Y such that LY = X or UY = X or the transpose cases.
//
  if (X.NumVectors()==1 && Y.NumVectors()==1) {
    double * xp = (double *) X[0];
    double * yp = (double *) Y[0];
    Epetra_Vector x(View, X.Map(), xp);
    Epetra_Vector y(View, Y.Map(), yp);
    return(Solve(Upper, Trans, UnitDiagonal, x, y));
  }
  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  if ((Upper) && (!UpperTriangular())) EPETRA_CHK_ERR(-2);
  if ((!Upper) && (!LowerTriangular())) EPETRA_CHK_ERR(-3);
  if ((!UnitDiagonal) && (NoDiagonal())) EPETRA_CHK_ERR(-4); // If matrix has no diagonal, we must use UnitDiagonal
  if ((!UnitDiagonal) && (NumMyDiagonals()<NumMyRows_)) EPETRA_CHK_ERR(-5); // Need each row to have a diagonal

  int i, j, j0, k;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;
  double diag;

  // If upper, point to last row
  if ((Upper && !Trans) || (!Upper && Trans)) {
    NumEntriesPerRow += NumMyRows_-1;
    Indices += NumMyRows_-1;
    Values += NumMyRows_-1;
  }

  double **Xp = (double**)X.Pointers();
  double **Yp = (double**)Y.Pointers();

  int NumVectors = X.NumVectors();

  if (!Trans) {
    

    if (Upper) {
      
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=NumMyRows_-1; i >=0; i--) {
	int      NumEntries = *NumEntriesPerRow--;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	if (!UnitDiagonal) diag = 1.0/RowValues[0]; // Take inverse of diagonal once for later use
	for (k=0; k<NumVectors; k++) {
	  double sum = 0.0;
	  for (j=j0; j < NumEntries; j++) sum += RowValues[j] * Yp[k][RowIndices[j]];
	  
	  if (UnitDiagonal) Yp[k][i] = Xp[k][i] - sum;
	  else Yp[k][i] = (Xp[k][i] - sum)*diag;
	}
      }
    }
    else {
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++ - j0;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	if (!UnitDiagonal) diag = 1.0/RowValues[NumEntries]; // Take inverse of diagonal once for later use
	for (k=0; k<NumVectors; k++) {
	  double sum = 0.0;
	  for (j=0; j < NumEntries; j++) sum += RowValues[j] * Yp[k][RowIndices[j]];
	  
	  if (UnitDiagonal) Yp[k][i] = Xp[k][i] - sum;
	  else Yp[k][i] = (Xp[k][i] - sum)*diag;
	}
      }
    }
  }
  // ***********  Transpose case *******************************

  else {

    for (k=0; k<NumVectors; k++) 
      if (Yp[k]!=Xp[k]) for (i=0; i < NumMyRows_; i++) Yp[k][i] = Xp[k][i]; // Initialize y for transpose multiply
    
    if (Upper) {
      
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      
      for (i=0; i < NumMyRows_; i++) {
	int      NumEntries = *NumEntriesPerRow++;
	int *    RowIndices = *Indices++;
	double * RowValues  = *Values++;
	if (!UnitDiagonal) diag = 1.0/RowValues[j0]; // Take inverse of diagonal once for later use
	for (k=0; k<NumVectors; k++) {
	  if (!UnitDiagonal) Yp[k][i] = Yp[k][i]*diag;
	  for (j=j0; j < NumEntries; j++) Yp[k][RowIndices[j]] -= RowValues[j] * Yp[k][i];
	}
      }
    }
    else {
      
      j0 = 1;
      if (NoDiagonal()) j0--; // Include first term if no diagonal
      
      for (i=NumMyRows_-1; i>=0; i--) {
	int      NumEntries = *NumEntriesPerRow-- - j0;
	int *    RowIndices = *Indices--;
	double * RowValues  = *Values--;
	for (k=0; k<NumVectors; k++) {
	   if (!UnitDiagonal)  Yp[k][i] = Yp[k][i]/Xp[k][i];
	   for (j=0; j < NumEntries; j++) Yp[k][RowIndices[j]] -= RowValues[j] * Yp[k][i];
        }
      }
    }
  }
  
  UpdateFlops(2*NumVectors*NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::InvRowSums(Epetra_Vector& x) const {
//
// Put inverse of the sum of absolute values of the ith row of A in x[i].
//

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  if (!Graph().RowMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the range of A
  int ierr = 0;
  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  double ** Values = Values_;
  double * xp = (double*)x.Values();


  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++;
    double * RowValues  = *Values++;
    double scale = 0.0;
    for (j=0; j < NumEntries; j++) scale += fabs(RowValues[j]);
    if (scale<Epetra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      xp[i] = Epetra_MaxDouble;
    }
    else
      xp[i] = 1.0/scale;
  }
  UpdateFlops(NumGlobalNonzeros());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
//=============================================================================
int Epetra_CrsMatrix::InvColSums(Epetra_Vector& x) const {
//
// Put inverse of the sum of absolute values of the jth column of A in x[j].
//

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  if (!Graph().DomainMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A
  
  double * xp = (double*)x.Values();
  Epetra_Vector * x_tmp = 0;
  int NumMyCols_ = NumMyCols();
  

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if (Importer()!=0) {
    x_tmp = new Epetra_Vector(ColMap()); // Create import vector if needed
    xp = (double*)x_tmp->Values();
  }
  int ierr = 0;
  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;

  for (i=0; i < NumMyCols_; i++) xp[i] = 0.0;

  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++;
    int *    ColIndices = *Indices++;
    double * RowValues  = *Values++;
    for (j=0; j < NumEntries; j++) xp[ColIndices[j]] += fabs(RowValues[j]);
  }

  if (Importer()!=0){
    x.Export(*x_tmp, *Importer(), Add); // Fill x with Values from import vector
    delete x_tmp;
    xp = (double*) x.Values();
  }
  // Invert values, don't allow them to get too large
  for (i=0; i < NumMyRows_; i++) {
    double scale = xp[i];
    if (scale<Epetra_MinDouble) {
      if (scale==0.0) ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
      else if (ierr!=1) ierr = 2;
      xp[i] = Epetra_MaxDouble;
    }
    else
      xp[i] = 1.0/scale;
  }
  UpdateFlops(NumGlobalNonzeros());
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::LeftScale(const Epetra_Vector& x) {
//
// This function scales the ith row of A by x[i].
//

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  if (!Graph().RangeMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the range of A

  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  double ** Values = Values_;
  double * xp = (double*)x.Values();


  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++;
    double * RowValues  = *Values++;
    double scale = xp[i];
    for (j=0; j < NumEntries; j++)  RowValues[j] *= scale;
  }
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  UpdateFlops(NumGlobalNonzeros());
  return(0);
}

//=============================================================================
int Epetra_CrsMatrix::RightScale(const Epetra_Vector& x) {
//
// This function scales the jth row of A by x[j].
//

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.
  if (!Graph().DomainMap().SameAs(x.Map())) EPETRA_CHK_ERR(-2); // x must have the same distribution as the domain of A

  double *xp = (double*)x.Values();
  Epetra_MultiVector * x_tmp = 0;

  // If we have a non-trivial importer, we must import elements that are permuted or are on other processors
  if (Importer()!=0) {
    x_tmp = new Epetra_Vector(ColMap()); // Create import vector if needed
    x_tmp->Import(x,*Importer(), Insert); // x_tmp will have all the values we need
    xp = (double*)x_tmp->Values();
  }

  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;

  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++;
    int *    ColIndices = *Indices++;
    double * RowValues  = *Values++;
    for (j=0; j < NumEntries; j++)  RowValues[j] *=  xp[ColIndices[j]];
  }
  if (x_tmp!=0) delete x_tmp;
  NormOne_ = -1.0; // Reset Norm so it will be recomputed.
  NormInf_ = -1.0; // Reset Norm so it will be recomputed.
  UpdateFlops(NumGlobalNonzeros());
  return(0);
}

//=============================================================================
double Epetra_CrsMatrix::NormInf() const {

  if (NormInf_>-1.0) return(NormInf_);

  int * NumEntriesPerRow = NumEntriesPerRow_;
  double ** Values = Values_;
  double Local_NormInf = 0.0;
  for (int i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++ ;
    double * RowValues  = *Values++;
    double sum = 0.0;
    for (int j=0; j < NumEntries; j++) sum += fabs(RowValues[j]);
    
    Local_NormInf = EPETRA_MAX(Local_NormInf, sum);
  }
  Comm().MaxAll(&Local_NormInf, &NormInf_, 1);
  UpdateFlops(NumGlobalNonzeros());
  return(NormInf_);
}
//=============================================================================
double Epetra_CrsMatrix::NormOne() const {

  if (NormOne_>-1.0) return(NormOne_);

  if (!Filled()) EPETRA_CHK_ERR(-1); // Matrix must be filled.

  Epetra_Vector * x = new Epetra_Vector(RowMap()); // Need temp vector for column sums
  
  double * xp = (double*)x->Values();
  Epetra_MultiVector * x_tmp = 0;
  int NumMyCols_ = NumMyCols();
  

  // If we have a non-trivial importer, we must export elements that are permuted or belong to other processors
  if (Importer()!=0) {
    x_tmp = new Epetra_Vector(ColMap()); // Create temporary import vector if needed
    xp = (double*)x_tmp->Values();
  }
  int i, j;
  int * NumEntriesPerRow = NumEntriesPerRow_;
  int ** Indices = Indices_;
  double ** Values = Values_;

  for (i=0; i < NumMyCols_; i++) xp[i] = 0.0;

  for (i=0; i < NumMyRows_; i++) {
    int      NumEntries = *NumEntriesPerRow++;
    int *    ColIndices = *Indices++;
    double * RowValues  = *Values++;
    for (j=0; j < NumEntries; j++) xp[ColIndices[j]] += fabs(RowValues[j]);
  }
  if (Importer()!=0) x->Export(*x_tmp, *Importer(), Add); // Fill x with Values from temp vector
  x->MaxValue(&NormOne_); // Find max
  if (x_tmp!=0) delete x_tmp;
  delete x;
  UpdateFlops(NumGlobalNonzeros());
  return(NormOne_);
}
//=========================================================================
int Epetra_CrsMatrix::CopyAndPermute(const Epetra_DistObject & Source, int NumSameIDs, 
				     int NumPermuteIDs, int * PermuteToLIDs,
				     int *PermuteFromLIDs) {
  
  const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Source);

  int i;
  
  int Row, NumEntries;
  int * Indices;
  double * Values;
  int FromRow, ToRow;
  
  // Do copy first
  if (NumSameIDs>0) {
    if (A.IndicesAreLocal()) {
      int MaxNumEntries = A.MaxNumEntries();
      Indices = new int[MaxNumEntries];  // Need some temporary space
      Values = new double[MaxNumEntries];  // Need some temporary space
      
      for (i=0; i<NumSameIDs; i++) {
				Row = GRID(i);
				assert(A.ExtractGlobalRowCopy(Row, MaxNumEntries, NumEntries, Values, Indices)==0); // Set pointers
				// Place into target matrix.  Depends on Epetra_DataAccess copy/view and static/dynamic graph.
				if (StaticGraph() || IndicesAreLocal())
					assert(ReplaceGlobalValues(Row, NumEntries, Values, Indices)==0);
				else
					assert(InsertGlobalValues(Row, NumEntries, Values, Indices)==0); 
      }
      delete [] Values;
      delete [] Indices;
    }
    else { // A.IndiceAreGlobal()
      for (i=0; i<NumSameIDs; i++) {
				Row = GRID(i);
				assert(A.ExtractGlobalRowView(Row, NumEntries, Values, Indices)==0); // Set pointers
				// Place into target matrix.  Depends on Epetra_DataAccess copy/view and static/dynamic graph.
				if (StaticGraph() || IndicesAreLocal())
					assert(ReplaceGlobalValues(Row, NumEntries, Values, Indices)==0); 
				else
					assert(InsertGlobalValues(Row, NumEntries, Values, Indices)==0); 
      }
    }	
  }
	
  // Do local permutation next
  if (NumPermuteIDs>0) {
    if (A.IndicesAreLocal()) {
      int MaxNumEntries = A.MaxNumEntries();
      Indices = new int[MaxNumEntries];  // Need some temporary space
      Values = new double[MaxNumEntries];  // Need some temporary space
      
      for (i=0; i<NumPermuteIDs; i++) {
				FromRow = A.GRID(PermuteFromLIDs[i]);
				ToRow = GRID(PermuteToLIDs[i]);
				assert(A.ExtractGlobalRowCopy(FromRow, MaxNumEntries, NumEntries, Values, Indices)==0); // Set pointers
				// Place into target matrix.  Depends on Epetra_DataAccess copy/view and static/dynamic graph.
				if (StaticGraph() || IndicesAreLocal())
					assert(ReplaceGlobalValues(ToRow, NumEntries, Values, Indices)==0);
				else
					assert(InsertGlobalValues(ToRow, NumEntries, Values, Indices)==0); 
      }
      delete [] Values;
      delete [] Indices;
    }
    else { // A.IndiceAreGlobal()
      for (i=0; i<NumPermuteIDs; i++) {
				FromRow = A.GRID(PermuteFromLIDs[i]);
				ToRow = GRID(PermuteToLIDs[i]);
				assert(A.ExtractGlobalRowView(FromRow, NumEntries, Values, Indices)==0); // Set pointers
				// Place into target matrix.  Depends on Epetra_DataAccess copy/view and static/dynamic graph.
				if (StaticGraph() || IndicesAreLocal())
					assert(ReplaceGlobalValues(ToRow, NumEntries, Values, Indices)==0); 
				else
					assert(InsertGlobalValues(ToRow, NumEntries, Values, Indices)==0); 
      }
    }
  }	
	
  return(0);
}

//=========================================================================
int Epetra_CrsMatrix::PackAndPrepare(const Epetra_DistObject & Source, 
				     int NumExportIDs, int * ExportLIDs,
				     int Nsend, int Nrecv,
				     int & LenExports, char * & Exports, int & LenImports, 
				     char * & Imports, 
				     int & SizeOfPacket, Epetra_Distributor & Distor) {
  

  const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Source);

  double * DoubleExports = 0;
  double * DoubleImports = 0;
  int GlobalMaxNumEntries = A.GlobalMaxNumEntries();
  // Will have GlobalMaxNumEntries doubles, GlobalMaxNumEntries +2 ints, pack them interleaved
  int DoublePacketSize = GlobalMaxNumEntries +  
                        (((GlobalMaxNumEntries+2)+sizeof(int)-1)*sizeof(int))/sizeof(double);
  SizeOfPacket = DoublePacketSize * sizeof(double); 




  if (DoublePacketSize*Nsend>LenExports) {
    if (LenExports>0) delete [] Exports;
    LenExports = DoublePacketSize*Nsend;
    DoubleExports = new double[LenExports];
    Exports = (char *) DoubleExports;
  }

  if (DoublePacketSize*Nrecv>LenImports) {
    if (LenImports>0) delete [] Imports;
    LenImports = DoublePacketSize*Nrecv;
    DoubleImports = new double[LenImports];
    Imports = (char *) DoubleImports;
  }

  if (NumExportIDs<=0) return(0); // All done if nothing to pack

  int NumEntries;
  int * Indices;
  double * Values;
  int FromRow;
  double * valptr, * dintptr;
  int * intptr;
  

  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  NumEntries, Number of indices in row.
  // next NumEntries: The actual indices for the row.
  // Any remaining space (of length GlobalMaxNumEntries - NumEntries ints) will be wasted but we need fixed
  //   sized segments for current communication routines.

  valptr = (double *) Exports;
  dintptr = valptr + GlobalMaxNumEntries;
  intptr = (int *) dintptr;
  for (int i=0; i<NumExportIDs; i++) {
    FromRow = A.GRID(ExportLIDs[i]);
    *intptr = FromRow;
    Values = valptr; 
    Indices = intptr + 2;
    EPETRA_CHK_ERR(A.ExtractGlobalRowCopy(FromRow, GlobalMaxNumEntries, NumEntries, Values, Indices));
    intptr[1] = NumEntries; // Load second slot of segment
    valptr += DoublePacketSize; // Point to next segment
    dintptr = valptr + GlobalMaxNumEntries;
    intptr = (int *) dintptr;
  }
    
  return(0);
}

//=========================================================================
int Epetra_CrsMatrix::UnpackAndCombine(const Epetra_DistObject & Source, 
				       int NumImportIDs, int * ImportLIDs, 
				       char * Imports, int & SizeOfPacket, 
				       Epetra_Distributor & Distor, 
				       Epetra_CombineMode CombineMode) {

  if (NumImportIDs<=0) return(0);

  if (   CombineMode != Add
      && CombineMode != Insert
      && CombineMode != Zero )
    EPETRA_CHK_ERR(-1); //Unsupported CombineMode, defaults to Zero


  const Epetra_CrsMatrix & A = dynamic_cast<const Epetra_CrsMatrix &>(Source);
  int NumEntries;
  int * Indices;
  double * Values;
  int ToRow;
  int i;
  
  double * valptr, *dintptr;
  int * intptr;
  int GlobalMaxNumEntries = A.GlobalMaxNumEntries();
  // Will have GlobalMaxNumEntries doubles, GlobalMaxNumEntries +2 ints, pack them interleaved
  int DoublePacketSize = GlobalMaxNumEntries +  
                        (((GlobalMaxNumEntries+2)+sizeof(int)-1)*sizeof(int))/sizeof(double);
  // Unpack it...


  // Each segment of Exports will be filled by a packed row of information for each row as follows:
  // 1st int: GRID of row where GRID is the global row ID for the source matrix
  // next int:  NumEntries, Number of indices in row.
  // next NumEntries: The actual indices for the row.
  // Any remaining space (of length GlobalMaxNumEntries - NumEntries ints) will be 
  //  wasted but we need fixed sized segments for current communication routines.

  valptr = (double *) Imports;
  dintptr = valptr + GlobalMaxNumEntries;
  intptr = (int *) dintptr;
    
  for (i=0; i<NumImportIDs; i++) {
    ToRow = GRID(ImportLIDs[i]);
    assert((intptr[0])==ToRow); // Sanity check
    NumEntries = intptr[1];
    Values = valptr; 
    Indices = intptr + 2; 
    if (CombineMode==Add) {
      if (StaticGraph() || IndicesAreLocal())
				// Add to any current values
				assert(SumIntoGlobalValues(ToRow, NumEntries, Values, Indices)==0);
      else
				// Insert values
				assert(InsertGlobalValues(ToRow, NumEntries, Values, Indices)>=0);
    }
    else if (CombineMode==Insert) {
      if (StaticGraph() || IndicesAreLocal())
				// Replace any current values
				assert(ReplaceGlobalValues(ToRow, NumEntries, Values, Indices)==0);
      else
				// Insert values
				assert(InsertGlobalValues(ToRow, NumEntries, Values, Indices)>=0);
    }
		
    valptr += DoublePacketSize; // Point to next segment
    dintptr = valptr + GlobalMaxNumEntries;
    intptr = (int *) dintptr;
  }
  return(0);
}

//=========================================================================

void Epetra_CrsMatrix::Print(ostream& os) const {
  int MyPID = RowMap().Comm().MyPID();
  int NumProc = RowMap().Comm().NumProc();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      const Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
      const Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
      const int             oldp = os.precision(12);
      if (MyPID==0) {
	os <<  "\nNumber of Global Rows        = "; os << NumGlobalRows(); os << endl;
	os <<    "Number of Global Cols        = "; os << NumGlobalCols(); os << endl;
	os <<    "Number of Global Diagonals   = "; os << NumGlobalDiagonals(); os << endl;
	os <<    "Number of Global Nonzeros    = "; os << NumGlobalNonzeros(); os << endl;
	os <<    "Global Maximum Num Entries   = "; os << GlobalMaxNumEntries(); os << endl;
	if (LowerTriangular()) os <<    " ** Matrix is Lower Triangular **"; os << endl;
	if (UpperTriangular()) os <<    " ** Matrix is Upper Triangular **"; os << endl;
	if (NoDiagonal())      os <<    " ** Matrix has no diagonal     **"; os << endl; os << endl;
      }

      os <<  "\nNumber of My Rows        = "; os << NumMyRows(); os << endl;
      os <<    "Number of My Cols        = "; os << NumMyCols(); os << endl;
      os <<    "Number of My Diagonals   = "; os << NumMyDiagonals(); os << endl;
      os <<    "Number of My Nonzeros    = "; os << NumMyNonzeros(); os << endl;
      os <<    "My Maximum Num Entries   = "; os << MaxNumEntries(); os << endl; os << endl;

      os << flush;
      
      // Reset os flags
      
      os.setf(olda,ios::adjustfield);
      os.setf(oldf,ios::floatfield);
      os.precision(oldp);
    }
    // Do a few global ops to give I/O a chance to complete
    Comm().Barrier();
    Comm().Barrier();
    Comm().Barrier();
  }

  {for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumMyRows1 = NumMyRows();
      int MaxNumIndices = MaxNumEntries();
      int * Indices  = new int[MaxNumIndices];
      double * Values  = new double[MaxNumIndices];
      int NumIndices;
      int i, j;

      if (MyPID==0) {
	os.width(8);
	os <<  "   Processor ";
	os.width(10);
	os <<  "   Row Index ";
	os.width(10);
	os <<  "   Col Index ";
	os.width(20);
	os <<  "   Value     ";
	os << endl;
      }
      for (i=0; i<NumMyRows1; i++) {
	int Row = GRID(i); // Get global row number
	ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, Values, Indices);
	
	for (j = 0; j < NumIndices ; j++) {   
	  os.width(8);
	  os <<  MyPID ; os << "    ";	
	  os.width(10);
	  os <<  Row ; os << "    ";	
	  os.width(10);
	  os <<  Indices[j]; os << "    ";
	  os.width(20);
	  os <<  Values[j]; os << "    ";
	  os << endl;
	}
      }

      delete [] Indices;
      delete [] Values;
      
      os << flush;
      
    }
    // Do a few global ops to give I/O a chance to complete
    RowMap().Comm().Barrier();
    RowMap().Comm().Barrier();
    RowMap().Comm().Barrier();
  }}

  return;
}
