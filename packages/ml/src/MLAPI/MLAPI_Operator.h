#ifndef ML_OPERATOR_H
#define ML_OPERATOR_H

#include "ml_include.h"
#include <iostream>
#include "ml_operator.h"
#include "ml_epetra_utils.h"
#include "ml_RowMatrix.h"
#include "Teuchos_RefCountPtr.hpp"
#include "MLAPI_Error.h"
#include "MLAPI_Space.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Operator_Box.h"

using namespace std;

namespace MLAPI {

/*!
 * \class Operator
 *
 * \brief Operator: basic class to define operators within MLAPI.
 *
 * \author Marzio Sala, SNL 9214
 *
 * \date Last updated on Feb-05.
 */

class Operator : public BaseOperator, public CompObject, public TimeObject {

public:
  //@{ \name Constructors and destructors.

  //! Default constructor.
  Operator() 
  {
    RCPOperatorBox_ = Teuchos::null;
  }

  //! Constructor with given already computed ML_Operator pointer.
  Operator(const Space& DomainSpace, const Space& RangeSpace,
           ML_Operator* Op, bool Ownership = true,
           Teuchos::RefCountPtr<ML_Operator_Box> AuxOp = Teuchos::null)
  {
    Reshape(DomainSpace, RangeSpace, Op, Ownership, AuxOp);
  }

  //! Constructor with given already FillComplete()'d object.
  Operator(const Space& DomainSpace, const Space& RangeSpace,
           Epetra_RowMatrix* Matrix, bool Ownership = true,
           Teuchos::RefCountPtr<ML_Operator_Box> AuxOp = Teuchos::null)
  {
    Reshape(DomainSpace, RangeSpace, Matrix, Ownership, AuxOp);
  }

  //! Copy constructor.
  Operator(const Operator& RHS) 
  {
    DomainSpace_       = RHS.GetDomainSpace();
    RangeSpace_        = RHS.GetRangeSpace();
    ColumnSpace_       = RHS.GetColumnSpace();
    RCPOperatorBox_    = RHS.GetRCPOperatorBox();
    RCPAuxOperatorBox_ = RHS.GetRCPAuxOperatorBox();
    RCPRowMatrix_      = RHS.GetRCPRowMatrix();

    SetLabel(RHS.GetLabel());
  }

  //! Destructor.
  ~Operator()
  {
    Destroy();
  }

  // @}
  // @{ \name Reshape methods

  //! Resets \c this object.
  void Reshape()
  {
    Destroy();
  }

  //! Reshape with given already computed ML_Operator pointer.
  void Reshape(const Space& DomainSpace, const Space& RangeSpace,
               ML_Operator* Op, bool Ownership = true,
               Teuchos::RefCountPtr<ML_Operator_Box> AuxOp = Teuchos::null)
  {
    RangeSpace_        = RangeSpace;
    DomainSpace_       = DomainSpace;

    RCPOperatorBox_    = Teuchos::rcp(new ML_Operator_Box(Op,Ownership));
    RCPRowMatrix_      = Teuchos::rcp(new ML_Epetra::RowMatrix(Op,&(GetEpetra_Comm())),Ownership);
    RCPAuxOperatorBox_ = AuxOp;

    BuildColumnSpace();
  }

  //! Reshape with given already FillComplete()'d object.
  void Reshape(const Space& DomainSpace, const Space& RangeSpace,
               Epetra_RowMatrix* Matrix, bool Ownership = true,
               Teuchos::RefCountPtr<ML_Operator_Box> AuxOp = Teuchos::null)
  {
    RangeSpace_ = RangeSpace;
    DomainSpace_ = DomainSpace;


    ML_Operator* Op = ML_Operator_Create(MLAPI::GetML_Comm());
    RCPOperatorBox_ = Teuchos::rcp(new ML_Operator_Box(Op,true));
    RCPAuxOperatorBox_ = AuxOp;

    RCPRowMatrix_ = Teuchos::rcp(Matrix,Ownership);
    Epetra2MLMatrix(RCPRowMatrix_.get(), GetML_Operator());

    BuildColumnSpace();

  }

  // @}
  // @{ \name Overloaded operators

  //! Makes \c this object equivalent to \c RHS.
  Operator& operator=(const Operator& RHS) 
  {
    Destroy();

    DomainSpace_    = RHS.GetDomainSpace();
    RangeSpace_     = RHS.GetRangeSpace();
    ColumnSpace_    = RHS.GetColumnSpace();
    RCPOperatorBox_ = RHS.GetRCPOperatorBox();
    RCPRowMatrix_   = RHS.GetRCPRowMatrix();
    
    SetLabel(RHS.GetLabel());
    return(*this);
  }

  //! Sets the label of \c this object.
  inline Operator& operator=(const string& Label)
  {
    SetLabel(Label);
    return(*this);
  }

  // @}
  // @{ \name Get and Set methods
  
  //! Returns a reference to the internally stored domain space.
  const Space GetOperatorDomainSpace() const {
    return(DomainSpace_);
  }

  //! Returns a reference to the internally stored range space.
  const Space GetOperatorRangeSpace() const {
    return(RangeSpace_);
  }

  //! Returns a reference to the internally stored domain space.
  inline const Space GetDomainSpace() const {
    return(DomainSpace_);
  }

  //! Returns a reference to the internally stored range space.
  inline const Space GetRangeSpace() const {
    return(RangeSpace_);
  }

  //! Returns a reference to the internally stored column space.
  inline const Space GetColumnSpace() const 
  {
    return(ColumnSpace_);
  }

  inline int GetNumGlobalRows() const 
  {
    return(GetRangeSpace().GetNumGlobalElements());
  }

  inline int GetNumMyRows() const 
  {
    return(GetRangeSpace().GetNumMyElements());
  }

  inline int GetNumGlobalCols() const 
  {
    return(GetDomainSpace().GetNumGlobalElements());
  }

  inline int GetNumMyCols() const 
  {
    return(GetDomainSpace().GetNumMyElements());
  }

  inline int GetNumGlobalNonzeros() const 
  {
    return(GetRowMatrix()->NumGlobalNonzeros());
  }

  inline int GetNumMyNonzeros() const 
  {
    return(GetRowMatrix()->NumMyNonzeros());
  }

  //! Returns the RefCountPtr of OperatorBox_.
  inline const Epetra_RowMatrix* GetRowMatrix() const
  {
    return(RCPRowMatrix_.get());
  }
  
  //! Returns the RefCountPtr of OperatorBox_.
  inline ML_Operator* GetML_Operator() const
  {
    return(GetRCPOperatorBox()->GetData());
  }

  //! Returns the RefCountPtr of OperatorBox_.
  inline const Teuchos::RefCountPtr<ML_Operator_Box>& GetRCPOperatorBox() const
  {
    return(RCPOperatorBox_);
  }

  //! Returns the RefCountPtr of AuxOperatorBox_.
  inline const Teuchos::RefCountPtr<ML_Operator_Box>& GetRCPAuxOperatorBox() const
  {
    return(RCPAuxOperatorBox_);
  }
  //! Returns the RefCountPtr of RowMatrix_
  inline const Teuchos::RefCountPtr<Epetra_RowMatrix>& GetRCPRowMatrix() const
  {
    return(RCPRowMatrix_);
  }

  // @}
  // @{ \name Mathematical methods.
  
  //! Applies \c this operator to LHS, returns the result in \c RHS.
  int Apply(const MultiVector& X, MultiVector& Y) const
  {
    ResetTimer();

    if (GetDomainSpace() != X.GetVectorSpace())
      ML_THROW("Domain spaces differ", -1);
    if (GetRangeSpace() != Y.GetVectorSpace())
      ML_THROW("Range spaces differ", -1);
    if (X.GetNumVectors() != Y.GetNumVectors())
      ML_THROW("Number of vectors differ", -1);
    if (GetML_Operator() == 0)
      ML_THROW("Operator not set", -1);
      
    int (*func)(ML_Operator*,int,double*,int,double*) = 
      GetML_Operator()->matvec->func_ptr;

    for (int v = 0 ; v < X.GetNumVectors() ; ++v) {
      double* x_ptr = (double*)&X(0) + v * X.GetMyLength();
      double* y_ptr = (double*)&Y(0) + v * Y.GetMyLength();
      (*func)(GetML_Operator(),X.GetMyLength(),x_ptr,
              Y.GetMyLength(), y_ptr);
    }

    // FIXME-RST: is there a way to get the flop count from
    // ML for this operator? The following is not always correct...

    UpdateFlops(2.0 * GetNumGlobalNonzeros());
    UpdateTime();

    return(0);
  }

  // @}
  // @{ \name Miscellaneous methods
  
  //! Prints basic information about \c this object.
  ostream& Print(std::ostream& os, const bool verbose = true) const
  {
    if (GetRCPOperatorBox().get() == 0) {
      if (GetMyPID() == 0) {
        os << endl;
        os << "*** MLAPI::Operator ***" << endl;
        os << "Label  = " << GetLabel() << endl;
        os << "Status = empty" << endl;
        os << endl;
      }
      return(os);
    }

    int    *bindx;
    double *val;
    int    allocated, row_length;
    ML_Operator* matrix = GetML_Operator();

    if (matrix->getrow == NULL) 
      ML_THROW("getrow not set", -1);

    if (GetMyPID() == 0) {
      os << endl;
      os << "*** MLAPI::Operator ***" << endl;
      os << "Label             = " << GetLabel() << endl;
      os << "Number of rows    = " << GetRangeSpace().GetNumGlobalElements() << endl;
      os << "Number of columns = " << GetDomainSpace().GetNumGlobalElements() << endl;
      os << "Flop count        = " << GetFlops() << endl;
      os << "Cumulative time   = " << GetTime() << endl;
      os << "MFlops rate       = " << 1.0e-6 * GetFlops() / GetTime() << endl;
      os << endl;
    }

    if (!verbose) 
      return(os);

    allocated = 100;
    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
    val   = (double *)  ML_allocate(allocated*sizeof(double));

    for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

      if (GetMyPID() == 0) {
        os.width(10);
        os << "ProcID";
        os.width(20);
        os << "Global Row";
        os.width(20);
        os << "Global Col";
        os.width(20);
        os << "Value" << endl;
        os << endl;
      }

      if (GetMyPID() == iproc) {

        for (int i = 0 ; i < matrix->getrow->Nrows; i++) {
          ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                            &row_length, 0);
          for  (int j = 0; j < row_length; j++) {
            int GlobalRow = GetRangeSpace()(i);
            int GlobalCol = GetColumnSpace()(bindx[j]);
            os.width(10);
            os << iproc;
            os.width(20);
            os << GlobalRow;
            os.width(20);
            os << GlobalCol;
            os.width(20);
            os << val[j] << endl;
          }
        }
      }
      Barrier();
    }

    if (GetMyPID() == 0)
      os << endl;

    ML_free(val);
    ML_free(bindx);
    return (os);
  }

  // @}
  
private:
  
  //! Build the column space, by computing the GID of all local columns.
  void BuildColumnSpace()
  {

    vector<double> dtemp;
    vector<int> GlobalElements;

    int Nrows = GetML_Operator()->getrow->Nrows;
    int Nghosts;
    if (GetML_Operator()->getrow->pre_comm == NULL) Nghosts = 0;
    else {
      if (GetML_Operator()->getrow->pre_comm->total_rcv_length <= 0)
        ML_CommInfoOP_Compute_TotalRcvLength(GetML_Operator()->getrow->pre_comm);
      Nghosts = GetML_Operator()->getrow->pre_comm->total_rcv_length;
    }

    dtemp.resize(Nrows + Nghosts);

    for (int i = 0 ; i < Nrows ; ++i) 
      dtemp[i] = 1.0 * GetDomainSpace()(i);
    for (int i = 0 ; i < Nghosts; ++i) 
      dtemp[i + Nrows] = -1;

    ML_exchange_bdry(&dtemp[0],GetML_Operator()->getrow->pre_comm,
                     GetML_Operator()->outvec_leng,
                     GetML_Comm(), ML_OVERWRITE,NULL);

    GlobalElements.resize(Nrows + Nghosts);

    for (int i = 0 ; i < Nrows + Nghosts ; ++i)
      GlobalElements[i] = (int)dtemp[i];

    ColumnSpace_.Reshape(-1, Nrows + Nghosts, &GlobalElements[0]);

    return;
  }

  //! Destroys all internal data and resets \c this object.
  void Destroy() 
  { 
    RangeSpace_.Reshape();
    DomainSpace_.Reshape();
    RCPOperatorBox_    = Teuchos::null;
    RCPRowMatrix_      = Teuchos::null;
    RCPAuxOperatorBox_ = Teuchos::null;
  }

  //! Domain space.
  Space DomainSpace_;
  //! Range space.
  Space RangeSpace_;
  //! Column space.
  Space ColumnSpace_;
  //! Container for the underlying ML_Operator pointer.
  Teuchos::RefCountPtr<ML_Operator_Box> RCPOperatorBox_;
  //! Container for the underlying ML_Operator pointer.
  Teuchos::RefCountPtr<ML_Operator_Box> RCPAuxOperatorBox_;
  //! Container for the underlying Epetra_RowMatrix pointer
  Teuchos::RefCountPtr<Epetra_RowMatrix> RCPRowMatrix_;

}; // Operator

} // namespace MLAPI

#endif // ML_OPERATOR_H
