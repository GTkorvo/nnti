#ifndef MLAPI_SERIALMATRIX_H
#define MLAPI_SERIALMATRIX_H

#include "ml_common.h"
#ifdef HAVE_ML_MLAPI

#include "ml_include.h"
//#include "ml_lapack.h"
#include "ml_comm.h"
#include "MLAPI_Error.h"
#include "MLAPI_BaseObject.h"
#include "MLAPI_Space.h"
#include "Epetra_RowMatrix.h"
#include <iomanip>

namespace MLAPI {

class SerialMatrix : public Epetra_RowMatrix, public BaseObject {

public:

  SerialMatrix(const Space& RowSpace, const Space& ColSpace)
  {
    RowSpace_ = RowSpace;
    ColSpace_ = ColSpace;

    NumMyRows_ = RowSpace_.GetNumMyElements();
    NumMyCols_ = ColSpace_.GetNumMyElements();
    
    NumMyNonzeros_ = 0;
    NumMyDiagonals_ = 0;

    if (GetNumProcs() != 1)
      ML_THROW("Class SerialMatrix can only be used for serial computations.", -1);

    RowMap_ = new Epetra_Map(NumMyRows_,0,GetEpetra_Comm());
    ColMap_ = new Epetra_Map(NumMyCols_,0,GetEpetra_Comm());

    ptr_.resize(NumMyRows_);
  }

  virtual int NumMyRowEntries(int MyRow, int & NumEntries) const
  {
#ifdef MLAPI_CHECK
    if (MyRow < 0 || MyRow >= NumMyRows())
      ML_THROW("Requested not valid row (" + GetString(MyRow) +").", -1);
#endif
    NumEntries = ptr_[MyRow].size();

    return(0);
  }

  virtual int MaxNumEntries() const
  {
    int res = 0, res_i = 0;

    for (int i = 0 ; i < NumMyRows() ; ++i) {
      NumMyRowEntries(i, res_i);
      if (res_i > res)
        res = res_i;
    }

    return(res);
  }

  virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
                               double *Values, int * Indices) const
  {
    NumMyRowEntries(MyRow, NumEntries);
    if (Length < NumEntries) ML_CHK_ERR(-1);
    if (MyRow < 0 || MyRow >= NumMyRows()) 
      ML_CHK_ERR(-2);

    int count = 0;
    for (where_ = ptr_[MyRow].begin() ; where_ != ptr_[MyRow].end() ; ++where_) {
      Indices[count] = where_->first;
      Values[count] = where_->second;
      ++count;
    }
    return(0);
  }

  virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
  {
#ifdef MLAPI_CHECK
    if (!Diagonal.Map().SameAs(RowMatrixRowMap()))
      ML_CHK_ERR(-1);
#endif

    Diagonal.PutScalar(0.0);
                       
    for (int i = 0 ; i < NumMyRows() ; ++i) {
      for (where_ = ptr_[i].begin() ; where_ != ptr_[i].end() ; ++where_) {
        if (where_->first == i) {
          Diagonal[i] = where_->second;
          break;
        }
      }
    }
    return(0);
  }

  virtual int Multiply(bool TransA, const Epetra_MultiVector& X, 
                       Epetra_MultiVector& Y) const
  {

    Y.PutScalar(0.0);

    if (!TransA) {
      for (int v = 0 ; v < X.NumVectors() ; ++v) {
        for (int i = 0 ; i < NumMyRows() ; ++i) {
          for (where_ = ptr_[i].begin() ; where_ != ptr_[i].end() ; ++where_) {
            Y[v][i] += (where_->second) * X[v][where_->first];
          }
        }
      }
    }
    else {
      for (int v = 0 ; v < X.NumVectors() ; ++v) {
        for (int i = 0 ; i < NumMyRows() ; ++i) {
          for (where_ = ptr_[i].begin() ; where_ != ptr_[i].end() ; ++where_) {
            Y[v][where_->first] += (where_->second) * X[v][i];
          }
        }
      }
    }
    
    return(0);
  }

  virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
                    Epetra_MultiVector& Y) const
  {
    ML_CHK_ERR(-1);
  }

  virtual int InvRowSums(Epetra_Vector& x) const
  {
    ML_CHK_ERR(-1);
  }

  virtual int LeftScale(const Epetra_Vector& x)
  {
    ML_CHK_ERR(-1);
  }

  virtual int InvColSums(Epetra_Vector& x) const
  {
    ML_CHK_ERR(-1);
  }

  virtual int RightScale(const Epetra_Vector& x)
  {
    ML_CHK_ERR(-1);
  }

  virtual bool Filled() const
  {
    return(true);
  }

  virtual double NormInf() const 
  {
    ML_CHK_ERR(-1);
  }

  virtual double NormOne() const
  {
    ML_CHK_ERR(-1);
  }

  virtual int NumGlobalNonzeros() const
  {
    return(NumMyNonzeros_);
  }

  virtual int NumGlobalRows() const
  {
     return(NumMyRows_);
  }

  virtual int NumGlobalCols() const
  {
    return(NumMyCols_);
  }

  virtual int NumGlobalDiagonals() const
  {
    return(NumMyDiagonals_);
  }

  virtual int NumMyNonzeros() const
  {
    return(NumMyNonzeros_);
  }

  virtual int NumMyRows() const
  {
    return(NumMyRows_);
  }

  virtual int NumMyCols() const
  {
    return(NumMyCols_);
  }

  virtual int NumMyDiagonals() const
  {
    return(NumMyDiagonals_);
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
    return(*RowMap_);
  }

  virtual const Epetra_Map & RowMatrixColMap() const
  {
    return(*ColMap_);
  }

  virtual const Epetra_Import * RowMatrixImporter() const
  {
    return(0);
  }

  virtual const Epetra_Map& OperatorDomainMap() const
  {
    return(*ColMap_);
  }

  virtual const Epetra_Map& OperatorRangeMap() const
  {
    return(*RowMap_);
  }

  virtual const Epetra_Map& Map() const
  {
    return(*ColMap_);
  }
    
  //@}

  virtual int SetUseTranspose(bool)
  {
    ML_CHK_ERR(-1);
  }
  
  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
  {
    return(Multiply(false, X, Y));
  }

  virtual int ApplyInverse(const Epetra_MultiVector& X,
                           Epetra_MultiVector& Y) const
  {
    ML_CHK_ERR(-1);
  }

  virtual const char* Label() const
  {
    return(GetLabel().c_str());
  }

  virtual bool UseTranspose() const
  {
    return(false);
  }

  virtual bool HasNormInf() const
  {
    return(false);
  }

  virtual const Epetra_Comm& Comm() const
  {
    return(GetEpetra_Comm());
  }

  inline double& operator()(const int row, const int col)
  {
#ifdef MLAPI_CHECK
    if (row < 0 || row >= NumMyRows())
      ML_THROW("Requested not valid row (" + GetString(row) +").", -1);
    if (col < 0 || row >= NumMyCols())
      ML_THROW("Requested not valid column (" + GetString(col) +").", -1);
#endif
    where_ = ptr_[row].find(col);

    if (where_ != ptr_[row].end())
      // return a reference to this guy
      return(where_->second);
    else {
      double& val = ptr_[row][col];
      val = 0.0;
      // track number of stored elements
      ++NumMyNonzeros_;
      // track number of diagonals 
      if (row == col)
        ++NumMyDiagonals_;
      // return a reference to this guy
      return(val);
    }
  }
           
  std::ostream& Print(std::ostream& os, const bool verbose = true) const
  {
    int Length = MaxNumEntries();
    vector<double> Values(Length);
    vector<int>    Indices(Length);

    os << endl;
    os << "*** MLAPI::SerialMatrix ***" << endl;
    os << "Label = " << GetLabel() << endl;
    os << "Number of rows = " << NumMyRows() << endl;
    os << "Number of columns = " << NumMyCols() << endl;
    os << endl;
    os.width(10); os << "row ID";
    os.width(10); os << "col ID";
    os.width(30); os << "value";
    os << endl;
    os << endl;

    for (int i = 0 ; i < NumMyRows() ; ++i) {
      int NnzRow = 0;
      ExtractMyRowCopy(i, Length, NnzRow, &Values[0], &Indices[0]);
      for (int j = 0 ; j < NnzRow ; ++j) {
        os.width(10); os << i;
        os.width(10); os << Indices[j];
        os.width(30); os << Values[j];
        os << endl;
      }
    }
    return(os);
  }

  Space RowSpace() const
  {
    return(RowSpace_);
  }

  Space ColSpace() const
  {
    return(ColSpace_);
  }


private:

  SerialMatrix(const SerialMatrix& rhs)
  {
  }

  SerialMatrix& operator=(const SerialMatrix& rhs)
  {
    return(*this);
  }

  int NumMyRows_;
  int NumMyCols_;
  int NumMyDiagonals_;
  int NumMyNonzeros_;

  mutable map<int,double>::iterator where_;
  mutable vector<map<int,double> > ptr_;

  Space ColSpace_;
  Space RowSpace_;

  Epetra_Map* RowMap_;
  Epetra_Map* ColMap_;

}; // class SerialMatrix

} // namespace MLAPI

#endif // ifdef HAVE_ML_MLAPI
#endif // ifndef MLAPI_SERIALMATRIX_H
