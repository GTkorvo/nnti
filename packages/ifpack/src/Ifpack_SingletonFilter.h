#ifndef IFPACK_SINGLETONFILTER_H
#define IFPACK_SINGLETONFILTER_H

#include "Ifpack_ConfigDefs.h"
#include "Epetra_RowMatrix.h"
class Epetra_Comm;
class Epetra_Map;
class Epetra_MultiVector;
class Epetra_Import;
class Epetra_BlockMap;

//! Ifpack_SingletonFilter: Filter based on matrix entries.
//
class Ifpack_SingletonFilter : public virtual Epetra_RowMatrix {

public:
  //! Constructor.
  Ifpack_SingletonFilter(Epetra_RowMatrix* Matrix);

  //! Destructor.
  virtual ~Ifpack_SingletonFilter();

  //! Returns the number of entries in MyRow.
  virtual inline int NumMyRowEntries(int MyRow, int & NumEntries) const
  {
    return(NumEntries_[MyRow]);
  }

  //! Returns the maximum number of entries.
  virtual int MaxNumEntries() const
  {
    return(MaxNumEntries_);
  }

  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;

  virtual int Multiply(bool TransA, const Epetra_MultiVector& X, 
		       Epetra_MultiVector& Y) const;

  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, 
		    const Epetra_MultiVector& X,
		    Epetra_MultiVector& Y) const;

  virtual int Apply(const Epetra_MultiVector& X,
		    Epetra_MultiVector& Y) const;

  virtual int ApplyInverse(const Epetra_MultiVector& X,
			   Epetra_MultiVector& Y) const;

  virtual int InvRowSums(Epetra_Vector& x) const
  {
    return(-98); // NOT IMPLEMENTED
  }

  virtual int LeftScale(const Epetra_Vector& x)
  {
    return(-98); // NOT IMPLEMENTED
  }

  virtual int InvColSums(Epetra_Vector& x) const
  {
    return(-98); // NOT IMPLEMENTED
  }

  virtual int RightScale(const Epetra_Vector& x) 
  {
    return(-98); // NOT IMPLEMENTED
  }

  virtual bool Filled() const
  {
    return(A_.Filled());
  }

  virtual double NormInf() const
  {
    return(-1.0);
  }

  virtual double NormOne() const
  {
    return(-1.0);
  }

  virtual int NumGlobalNonzeros() const
  {
    return(NumNonzeros_);
  }

  virtual int NumGlobalRows() const
  {
    return(NumRows_);
  }

  virtual int NumGlobalCols() const
  {
    return(NumRows_);
  }

  virtual int NumGlobalDiagonals() const
  {
    return(NumRows_);
  }

  virtual int NumMyNonzeros() const
  {
    return(NumNonzeros_);
  }

  virtual int NumMyRows() const
  {
    return(NumRows_);
  }

  virtual int NumMyCols() const
  {
    return(NumRows_);
  }

  virtual int NumMyDiagonals() const
  {
    return(NumRows_);
  }

  virtual bool LowerTriangular() const
  {
    return(false);
  }

  virtual bool UpperTriangular() const
  {
    return(false);
  }

  virtual const Epetra_Map & RowMatrixRowMap() const
  {
    return(*Map_);
  }

  virtual const Epetra_Map & RowMatrixColMap() const
  {
    return(*Map_);
  }

  virtual const Epetra_Import * RowMatrixImporter() const
  {
    return(A_.RowMatrixImporter());
  }

  int SetUseTranspose(bool UseTranspose)
  {
    return(A_.SetUseTranspose(UseTranspose));
  }

  bool UseTranspose() const 
  {
    return(A_.UseTranspose());
  }

  bool HasNormInf() const
  {
    return(false);
  }

  const Epetra_Comm & Comm() const
  {
    return(A_.Comm());
  }

  const Epetra_Map & OperatorDomainMap() const 
  {
    return(*Map_);
  }

  const Epetra_Map & OperatorRangeMap() const 
  {
    return(*Map_);
  }

  const Epetra_BlockMap& Map() const 
  {
    return(*(const Epetra_BlockMap*)Map_);
  }

  const char* Label() const{
    return(Label_);
  }

  int SolveSingletons(const Epetra_MultiVector& RHS, 
		      Epetra_MultiVector& LHS);

  int CreateReducedRHS(const Epetra_MultiVector& LHS,
		       const Epetra_MultiVector& RHS, 
		       Epetra_MultiVector& ReducedRHS);

  int UpdateLHS(const Epetra_MultiVector& ReducedLHS,
		Epetra_MultiVector& LHS);

private:

  //! Pointer to the matrix to be preconditioned.
  Epetra_RowMatrix& A_;

  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable vector<int> Indices_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable vector<double> Values_;
  //! Label for \c this object.
  char Label_[80];
  
  int NumSingletons_;
  vector<int> SingletonIndex_;

  vector<int> Reorder_;
  vector<int> InvReorder_;

  vector<int> NumEntries_;

  int NumRows_;
  int NumRowsA_;
  int MaxNumEntries_;
  int MaxNumEntriesA_;
  int NumNonzeros_;
  Epetra_Map* Map_;
  
  Epetra_Vector* Diagonal_;

};

#endif /* IFPACK_SINGLETONFILTER_H */
