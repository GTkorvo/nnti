#ifndef IFPACK_BLOCKPRECONDITIONER_H
#define IFPACK_BLOCKPRECONDITIONER_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_Preconditioner.h" 
#include "Ifpack_Partitioner.h"
#include "Ifpack_LinearPartitioner.h"
#include "Ifpack_GreedyPartitioner.h"
#include "Ifpack_METISPartitioner.h"
#include "Ifpack_EquationPartitioner.h"
#include "Ifpack_UserPartitioner.h"
#include "Ifpack_Graph_Epetra_RowMatrix.h"
#include "Ifpack_DenseContainer.h" 
#include "Ifpack_Utils.h" 
#include "Teuchos_ParameterList.hpp"
#include "Epetra_RowMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_Import.h"

static const int IFPACK_JACOBI = 0;
static const int IFPACK_GS = 1;
static const int IFPACK_SGS = 2;


//! Ifpack_BlockRelaxation: a class to define block relaxation preconditioners of Epetra_RowMatrix's.

/*! The Ifpack_BlockRelaxation class enables the construction of
  block relaxation
  preconditioners of an Epetra_RowMatrix. Ifpack_PointRelaxation 
  is derived from 
  the Ifpack_Preconditioner class, which is derived from Epetra_Operator.
  Therefore this object can be used as preconditioner everywhere an
  ApplyInverse() method is required in the preconditioning step.
 
  The class currently support:
  - block Jacobi;
  - block Gauss-Seidel;
  - symmetric block Gauss-Seidel.
  
  The idea of block relaxation method is to extend their point relaxation
  counterpart (implemented in Ifpack_PointRelaxation), by working on a
  group of equation simulteneously. Generally, larger blocks result
  in better convergence and increased robusteness.

  The user can decide:
  - the number of blocks (say, NumBlocks). If NumBlocks is equal to the
    number of rows, then the resulting scheme is equivalent to
    a point relaxation scheme;
  - how to apply the inverse of each diagonal block, by choosing a dense
    container or a sparse container. The implementation of
    block relaxation schemes requires the application of the
    inverse of each diagonal block. This can be done using LAPACK (dense 
    container), or any Ifpack_Preconditioner derived class (sparse
    container);
  - blocks can be defined using a linear decomposition, by a simple greedy
    algorithm, or by resorting to METIS.

The following is an example of usage of this preconditioner with dense
containers. First, we include the header files:
\code
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_BlockPreconditioner.h"
#include "Ifpack_DenseContainer.h"
\endcode

Then, we declare the preconditioner. Note that this is done through
the class Ifpack_AdditiveSchwarz (see note below in this section).
\code
// A is an Epetra_RowMatrix
// List is a Teuchos::ParameterList
Ifpack_AdditiveSchwarz<Ifpack_BlockRelaxation<Ifpack_DenseContainer> > > Prec(A);
IFPACK_CHK_ERR(Prec.SetParameters(List));
IFPACK_CHK_ERR(Prec.Initialize());
IFPACK_CHK_ERR(Prec.Compute());

// action of the preconditioner is given by ApplyInverse()
// Now use it in AztecOO, solver is an AztecOO object
solver.SetPrecOperator(&Prec);
\endcode

<P>The complete list of supported parameters is reported in page \ref ifp_params. For a presentation of basic relaxation schemes, please refer to page
\ref Ifpack_PointRelaxation.

\note The ApplyInverse() implementation of this class is \e not AztecOO
complaint, as it does assume that the two input vectors X and Y actually
refer to two different memory location. In fact, this case is handled
by class Ifpack_AdditiveSchwarz, which takes care of calling methods ApplyInverse()
of Ifpack_BlockRelaxation with two vectors pointing to different memory
locations.

\author Marzio Sala, SNL 9214.

\date Last modified: Nov-04.
  
*/
template<typename T>
class Ifpack_BlockRelaxation : public Ifpack_Preconditioner {

public:

  //@{ \name Constructors/Destructors
  //! Ifpack_BlockRelaxation constructor with given Epetra_RowMatrix.
  /*! Creates an Ifpack_Preconditioner preconditioner. 
   *
   * \param In
   * Matrix - Pointer to matrix to be preconditioned.
   */
  Ifpack_BlockRelaxation(const Epetra_RowMatrix* Matrix);

  virtual ~Ifpack_BlockRelaxation();

  //@}

  //@{ \name Mathematical functions.

  //! Applies the matrix to an Epetra_MultiVector.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
    */
  virtual int Apply(const Epetra_MultiVector& X, 
		    Epetra_MultiVector& Y) const;

  //! Applies the block Jacobi preconditioner to X, returns the result in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method 
    must support the case where X and Y are the same object.
    */
  virtual int ApplyInverse(const Epetra_MultiVector& X, 
			   Epetra_MultiVector& Y) const;

  //! Returns the infinity norm of the global matrix (not implemented)
  virtual double NormInf() const
  {
    return(-1.0);
  }
  //@}

  //@{ \name Atribute access functions

  virtual int SetUseTranspose(bool UseTranspose)
  {
    IFPACK_CHK_ERR(-1); // FIXME: can I work with the transpose?
  }

  virtual char * Label() const;
 
  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const
  {
    return(false);
  }

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const
  {
    return(false);
  }

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const;

  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map & OperatorDomainMap() const;

  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map & OperatorRangeMap() const;
  //@}

  //! Returns the number local blocks.
  int NumLocalBlocks() const 
  {
    return(NumLocalBlocks_);
  }

  //! Returns \c true if the preconditioner has been successfully computed.
  virtual bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Returns \c true if the preconditioner has been successfully computed.
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Sets all the parameters for the preconditioner.
  virtual int SetParameters(Teuchos::ParameterList& List);

  //! Initializes the preconditioner.
  virtual int Initialize();

  //! Computes the preconditioner.
  virtual int Compute();

  virtual const Epetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                         const int MaxIters = 1550,
                         const double Tol = 1e-9,
			 Epetra_RowMatrix* Matrix = 0)
  {
    return(-1.0);
  }

  virtual double Condest() const
  {
    return(-1.0);
  }

  std::ostream& Print(std::ostream& os) const;

  //! Returns the number of calls to Initialize().
  virtual int NumInitialize() const
  {
    return(NumInitialize_);
  }

  //! Returns the number of calls to Compute().
  virtual int NumCompute() const
  {
    return(NumCompute_);
  }

  //! Returns the number of calls to ApplyInverse().
  virtual int NumApplyInverse() const
  {
    return(NumApplyInverse_);
  }

  //! Returns the time spent in Initialize().
  virtual double InitializeTime() const
  {
    return(InitializeTime_);
  }

  //! Returns the time spent in Compute().
  virtual double ComputeTime() const
  {
    return(ComputeTime_);
  }

  //! Returns the time spent in ApplyInverse().
  virtual double ApplyInverseTime() const
  {
    return(ApplyInverseTime_);
  }

  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

private:

  //! Copy constructor (PRIVATE, should not be used).
  Ifpack_BlockRelaxation(const Ifpack_BlockRelaxation& rhs);

  //! operator= (PRIVATE, should not be used).
  Ifpack_BlockRelaxation & operator=(const Ifpack_BlockRelaxation& rhs)
  {
    return(*this);
  }

  virtual int ApplyInverseJacobi(const Epetra_MultiVector& X, 
                                 Epetra_MultiVector& Y) const;

  virtual int DoJacobi(const Epetra_MultiVector& X, 
                                  Epetra_MultiVector& Y) const;

  virtual int ApplyInverseGS(const Epetra_MultiVector& X, 
                             Epetra_MultiVector& Y) const;

  virtual int DoGaussSeidel(Epetra_MultiVector& X, 
                            Epetra_MultiVector& Y) const;

  virtual int ApplyInverseSGS(const Epetra_MultiVector& X, 
                              Epetra_MultiVector& Y) const;

  virtual int DoSGS(const Epetra_MultiVector& X,
                    Epetra_MultiVector& Xtmp,
                    Epetra_MultiVector& Y) const;

  int ExtractSubmatrices();

  // @{ Initializations, timing and flops

  //! If true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to ApplyInverse().
  mutable int NumApplyInverse_;
  //! Contains the time for all successful calls to Initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to ApplyInverse().
  mutable double ApplyInverseTime_;
  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable double ApplyInverseFlops_;
  // @}

  // @{ Settings
  //! Number of preconditioning sweeps.
  int NumSweeps_;
  //! Damping parameter.
  double DampingFactor_;
  //! Number of local blocks
  int NumLocalBlocks_;
  //! Parameters list to be used to solve on each subblock
  Teuchos::ParameterList List_;
  // @}

  // @{ Other data
  //! Containers_[i] contains all the necessary information to solve on each subblock.
  //! Pointers to the matrix to be preconditioned.
  const Epetra_RowMatrix* Matrix_;
  mutable vector<T*> Containers_;
  //! Contains information about non-overlapping partitions.
  Ifpack_Partitioner* Partitioner_;
  string PartitionerType_;
  int PrecType_;
  //! Label for \c this object
  string Label_;
  //! If \c true, starting solution is the zero vector.
  bool ZeroStartingSolution_;
  Ifpack_Graph* Graph_;
  //! Weights for overlapping Jacobi only.
  Epetra_Vector* W_;
  // Level of overlap among blocks (for Jacobi only).
  int OverlapLevel_;
  Epetra_Time* Time_;
  bool IsParallel_;
  Epetra_Import* Importer_;
  // @}

}; // class Ifpack_BlockRelaxation

//==============================================================================
template<typename T>
Ifpack_BlockRelaxation<T>::
Ifpack_BlockRelaxation(const Epetra_RowMatrix* Matrix) :
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0),
  ApplyInverseFlops_(0),
  NumSweeps_(1),
  DampingFactor_(1.0),
  NumLocalBlocks_(1),
  Matrix_(Matrix),
  Containers_(0),
  Partitioner_(0),
  PartitionerType_("greedy"),
  PrecType_(IFPACK_JACOBI),
  ZeroStartingSolution_(true),
  Graph_(0),
  W_(0),
  OverlapLevel_(0),
  Time_(0),
  IsParallel_(false),
  Importer_(0)
{
  Time_ = new Epetra_Time(Matrix->Comm());
  if (Matrix->Comm().NumProc() != 1)
    IsParallel_ = true;
}

//==============================================================================
template<typename T>
Ifpack_BlockRelaxation<T>::~Ifpack_BlockRelaxation()
{
  for (int i = 0 ; i < NumLocalBlocks() ; ++i)
    delete Containers_[i];
  if (Partitioner_)
    delete Partitioner_;
  if (Graph_)
    delete Graph_;
  if (W_)
    delete W_;
  if (Time_)
    delete Time_;
  if (Importer_)
    delete Importer_;
}

//==============================================================================
template<typename T>
char* Ifpack_BlockRelaxation<T>::Label() const
{
  return(const_cast<char*>(Label_.c_str()));
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  IFPACK_RETURN(Matrix().Apply(X,Y));
}

//==============================================================================
template<typename T>
const Epetra_Comm& Ifpack_BlockRelaxation<T>::
Comm() const
{
  return(Matrix().Comm());
}

//==============================================================================
template<typename T>
const Epetra_Map& Ifpack_BlockRelaxation<T>::
OperatorDomainMap() const
{
  return(Matrix().OperatorDomainMap());
}

//==============================================================================
template<typename T>
const Epetra_Map& Ifpack_BlockRelaxation<T>::
OperatorRangeMap() const
{
  return(Matrix().OperatorRangeMap());
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::ExtractSubmatrices()
{

  if (Partitioner_ == 0)
    IFPACK_CHK_ERR(-1);

  NumLocalBlocks_ = Partitioner_->NumLocalParts();

  Containers_.resize(NumLocalBlocks());

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    int rows = Partitioner_->NumRowsInPart(i);
    Containers_[i] = new T(rows);
    
    Ifpack_DenseContainer* DC = 0;
    DC = dynamic_cast<Ifpack_DenseContainer*>(Containers_[i]);

    if (Containers_[i] == 0)
      IFPACK_CHK_ERR(-10);
    
    IFPACK_CHK_ERR(Containers_[i]->SetParameters(List_));
    IFPACK_CHK_ERR(Containers_[i]->Initialize());

    // set "global" ID of each partitioner row
    for (int j = 0 ; j < rows ; ++j) {
      int LRID = (*Partitioner_)(i,j);
      Containers_[i]->ID(j) = LRID;
    }

    IFPACK_CHK_ERR(Containers_[i]->Compute(*Matrix_));
    ComputeFlops_ += Containers_[i]->ComputeFlops();

  }

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::Compute()
{

  Time_->ResetStartTime();

  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());
  IsComputed_ = false;

  if (Matrix().NumGlobalRows() != Matrix().NumGlobalCols())
    IFPACK_CHK_ERR(-1); // only square matrices

  IFPACK_CHK_ERR(ExtractSubmatrices());
  
  if (IsParallel_ && PrecType_ != IFPACK_JACOBI) {
    // not needed by Jacobi (done by matvec)
    Importer_ = new Epetra_Import(Matrix().RowMatrixColMap(),
                                  Matrix().RowMatrixRowMap());
    assert (Importer_ != 0);
  }
  IsComputed_ = true;
  ComputeTime_ += Time_->ElapsedTime();
  ++NumCompute_;

  return(0);

}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  Time_->ResetStartTime();

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-4);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-3);

  switch (PrecType_) {
  case IFPACK_JACOBI:
    IFPACK_CHK_ERR(ApplyInverseJacobi(X,Y));
    break;
  case IFPACK_GS:
    IFPACK_CHK_ERR(ApplyInverseGS(X,Y));
    break;
  case IFPACK_SGS:
    IFPACK_CHK_ERR(ApplyInverseSGS(X,Y));
    break;
  }

  ApplyInverseTime_ += Time_->ElapsedTime();
  ++NumApplyInverse_;

  return(0);
}

//==============================================================================
// This method in general will not work with AztecOO if used
// outside Ifpack_AdditiveSchwarz and OverlapLevel_ != 0
//
template<typename T>
int Ifpack_BlockRelaxation<T>::
ApplyInverseJacobi(const Epetra_MultiVector& X, 
                   Epetra_MultiVector& Y) const
{

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  // do not compute the residual in this case
  if (NumSweeps_ == 1 && ZeroStartingSolution_) {
    IFPACK_RETURN(DoJacobi(X,Y));
  }

  Epetra_MultiVector AX(Y);

  for (int j = 0; j < NumSweeps_ ; j++) {
    IFPACK_CHK_ERR(Apply(Y,AX));
    IFPACK_CHK_ERR(AX.Update(1.0,X,-1.0));
    IFPACK_CHK_ERR(DoJacobi(AX,Y));
  }

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
DoJacobi(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (OverlapLevel_ == 0) {

    for (int i = 0 ; i < NumLocalBlocks() ; ++i) {
     
      // may happen that a partition is empty
      if (Containers_[i]->NumRows() == 0) 
        continue;

      int LID;

      // extract RHS from X
      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
        LID = Containers_[i]->ID(j);
        for (int k = 0 ; k < X.NumVectors() ; ++k) {
          Containers_[i]->RHS(j,k) = X[k][LID];
        }
      }

      // apply the inverse of each block
      IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());

      // copy back into solution vector Y
      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
        LID = Containers_[i]->ID(j);
        for (int k = 0 ; k < Y.NumVectors() ; ++k) {
          Y[k][LID] += DampingFactor_ * Containers_[i]->LHS(j,k);
        }
      }
    }
  }
  else {

    for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

      // may happen that a partition is empty
      if (Containers_[i]->NumRows() == 0) 
        continue;

      int LID;

      // extract RHS from X
      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
        LID = Containers_[i]->ID(j);
        for (int k = 0 ; k < X.NumVectors() ; ++k) {
          Containers_[i]->RHS(j,k) = (*W_)[LID] * X[k][LID];
        }
      }

      // apply the inverse of each block
      IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());

      // copy back into solution vector Y
      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
        LID = Containers_[i]->ID(j);
        for (int k = 0 ; k < Y.NumVectors() ; ++k) {
          Y[k][LID] += DampingFactor_ * (*W_)[LID] * Containers_[i]->LHS(j,k);
        }
      }
    }
  }
  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
ApplyInverseGS(const Epetra_MultiVector& X, 
               Epetra_MultiVector& Y) const
{

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  Epetra_MultiVector Xcopy(X);
  for (int j = 0; j < NumSweeps_ ; j++) {
    IFPACK_CHK_ERR(DoGaussSeidel(Xcopy,Y));
    if (j != NumSweeps_ - 1)
      Xcopy = X;
  }

  return(0);

}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
DoGaussSeidel(Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  // cycle over all local subdomains

  int Length = Matrix().MaxNumEntries();
  vector<int> Indices(Length);
  vector<double> Values(Length);

  int NumMyRows = Matrix().NumMyRows();
  int NumVectors = X.NumVectors();

  // an additonal vector is needed by parallel computations
  // (note that applications through Ifpack_AdditiveSchwarz
  // are always seen are serial)
  Epetra_MultiVector* Y2;
  if (IsParallel_)
    Y2 = new Epetra_MultiVector(Importer_->TargetMap(), NumVectors);
  else
    Y2 = &Y;

  double** y_ptr;
  double** y2_ptr;
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);

  // data exchange is here, once per sweep
  if (IsParallel_)
    IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    // may happen that a partition is empty
    if (Containers_[i]->NumRows() == 0) 
      continue;

    int LID;

    // update from previous block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);

      int NumEntries;
      IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(LID, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int k = 0 ; k < NumEntries ; ++k) {
        int col = Indices[k];

          for (int kk = 0 ; kk < NumVectors ; ++kk) {
            X[kk][LID] -= Values[k] * y2_ptr[kk][col];
          }
      }
    }

    // solve with this block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        Containers_[i]->RHS(j,k) = X[k][LID];
      }
    }

    IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());
    ApplyInverseFlops_ += Containers_[i]->ApplyInverseFlops();

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        y2_ptr[k][LID] += DampingFactor_ * Containers_[i]->LHS(j,k);
      }
    }
  }

  // Attention: this is delicate... Not all combinations
  // of Y2 and Y will always work (tough for ML it should be ok)
  if (IsParallel_)
    for (int m = 0 ; m < NumVectors ; ++m) 
      for (int i = 0 ; i < NumMyRows ; ++i)
        y_ptr[m][i] = y2_ptr[m][i];

  if (IsParallel_)
    delete Y2;

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
ApplyInverseSGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  Epetra_MultiVector Xcopy(X);
  for (int j = 0; j < NumSweeps_ ; j++) {
    IFPACK_CHK_ERR(DoSGS(X,Xcopy,Y));
    if (j != NumSweeps_ - 1)
      Xcopy = X;
  }
  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::
DoSGS(const Epetra_MultiVector& X, Epetra_MultiVector& Xcopy, 
      Epetra_MultiVector& Y) const
{

  int NumMyRows = Matrix().NumMyRows();
  int NumVectors = X.NumVectors();

  int Length = Matrix().MaxNumEntries();
  vector<int> Indices;
  vector<double> Values;
  Indices.resize(Length);
  Values.resize(Length);

  // an additonal vector is needed by parallel computations
  // (note that applications through Ifpack_AdditiveSchwarz
  // are always seen are serial)
  Epetra_MultiVector* Y2;
  if (IsParallel_)
    Y2 = new Epetra_MultiVector(Importer_->TargetMap(), NumVectors);
  else
    Y2 = &Y;

  double** y_ptr;
  double** y2_ptr;
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);

  // data exchange is here, once per sweep
  if (IsParallel_)
    IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    // may happen that a partition is empty
    if (Containers_[i]->NumRows() == 0) 
      continue;

    int LID;

    // update from previous block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);

      int NumEntries;
      IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(LID, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int k = 0 ; k < NumEntries ; ++k) {
        int col = Indices[k];

        for (int kk = 0 ; kk < NumVectors ; ++kk) {
          Xcopy[kk][LID] -= Values[k] * y2_ptr[kk][col];
        }
      }
    }

    // solve with this block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        Containers_[i]->RHS(j,k) = Xcopy[k][LID];
      }
    }

    IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());
    ApplyInverseFlops_ += Containers_[i]->ApplyInverseFlops();

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        y2_ptr[k][LID] += DampingFactor_ * Containers_[i]->LHS(j,k);
      }
    }
  }

  Xcopy = X;

  for (int i = NumLocalBlocks() - 1; i >=0 ; --i) {

    if (Containers_[i]->NumRows() == 0) 
      continue;

    int LID;

    // update from previous block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);

      int NumEntries;
      IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(LID, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int k = 0 ; k < NumEntries ; ++k) {
        int col = Indices[k];

          for (int kk = 0 ; kk < NumVectors ; ++kk) {
            Xcopy[kk][LID] -= Values[k] * y2_ptr[kk][col];
          }
      }
    }

    // solve with this block

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        Containers_[i]->RHS(j,k) = Xcopy[k][LID];
      }
    }

    IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());
    ApplyInverseFlops_ += Containers_[i]->ApplyInverseFlops();

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      LID = Containers_[i]->ID(j);
      for (int k = 0 ; k < NumVectors ; ++k) {
        y2_ptr[k][LID] += DampingFactor_ * Containers_[i]->LHS(j,k);
      }
    }
  }

  // Attention: this is delicate... Not all combinations
  // of Y2 and Y will always work (tough for ML it should be ok)
  if (IsParallel_)
    for (int m = 0 ; m < NumVectors ; ++m) 
      for (int i = 0 ; i < NumMyRows ; ++i)
        y_ptr[m][i] = y2_ptr[m][i];

  if (IsParallel_)
    delete Y2;

  return(0);
}

//==============================================================================
template<typename T>
ostream& Ifpack_BlockRelaxation<T>::Print(ostream & os) const
{

  string PT;
  if (PrecType_ == IFPACK_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == IFPACK_GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == IFPACK_SGS)
    PT = "symmetric Gauss-Seidel";

  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_BlockRelaxation, " << PT << endl;
    os << "Sweeps = " << NumSweeps_ << endl;
    os << "Damping factor = " << DampingFactor_;
    if (ZeroStartingSolution_) 
      cout << ", using zero starting solution" << endl;
    else
      cout << ", using input starting solution" << endl;
    os << "Number of local blocks = " << Partitioner_->NumLocalParts() << endl;
    //os << "Condition number estimate = " << Condest_ << endl;
    os << "Global number of rows            = " << Matrix_->NumGlobalRows() << endl;
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize_ 
       << "  " << std::setw(15) << InitializeTime_ 
       << "                -              -" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute_ 
       << "  " << std::setw(15) << ComputeTime_
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_ 
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_ / ComputeTime_ << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse_ 
       << "  " << std::setw(15) << ApplyInverseTime_
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_ 
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_ / ApplyInverseTime_ << endl;
    os << "================================================================================" << endl;
    os << endl;
  }

  return(os);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::SetParameters(Teuchos::ParameterList& List)
{

  string PT;
  if (PrecType_ == IFPACK_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == IFPACK_GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == IFPACK_SGS)
    PT = "symmetric Gauss-Seidel";

  PT = List.get("relaxation: type", PT);

  if (PT == "Jacobi") {
    PrecType_ = IFPACK_JACOBI;
  }
  else if (PT == "Gauss-Seidel") {
    PrecType_ = IFPACK_GS;
  }
  else if (PT == "symmetric Gauss-Seidel") {
    PrecType_ = IFPACK_SGS;
  } else {
    cerr << "Option `relaxation: type' has an incorrect value ("
      << PT << ")'" << endl;
    cerr << "(file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
    exit(EXIT_FAILURE);
  }

  NumSweeps_ = List.get("relaxation: sweeps",NumSweeps_);
  DampingFactor_ = List.get("relaxation: damping factor", DampingFactor_);
  ZeroStartingSolution_ = List.get("relaxation: zero starting solution", 
                                   ZeroStartingSolution_);
  PartitionerType_ = List.get("partitioner: type", PartitionerType_);
  NumLocalBlocks_ = List.get("partitioner: local parts", NumLocalBlocks_);
  // only Jacobi can work with overlap among local domains,
  OverlapLevel_ = List.get("partitioner: overlap", OverlapLevel_);

  // check parameters
  if (PrecType_ != IFPACK_JACOBI)
    OverlapLevel_ = 0;
  if (NumLocalBlocks_ < 0)
    NumLocalBlocks_ = Matrix().NumMyRows() / (-NumLocalBlocks_);
  // other checks are performed in Partitioner_
  
  // copy the list as each subblock's constructor will
  // require it later
  List_ = List;

  // set the label
  string PT2;
  if (PrecType_ == IFPACK_JACOBI)
    PT2 = "BJ";
  else if (PrecType_ == IFPACK_GS)
    PT2 = "BGS";
  else if (PrecType_ == IFPACK_SGS)
    PT2 = "BSGS";
  Label_ = "IFPACK (" + PT2 + ", sweeps=" 
    + Ifpack_toString(NumSweeps_) + ", damping="
    + Ifpack_toString(DampingFactor_) + ", blocks="
    + Ifpack_toString(NumLocalBlocks()) + ")";

  return(0);
}

//==============================================================================
template<typename T>
int Ifpack_BlockRelaxation<T>::Initialize()
{
  IsInitialized_ = false;

  if (Partitioner_)
    delete Partitioner_;
  if (Graph_)
    delete Graph_;

  Graph_ = new Ifpack_Graph_Epetra_RowMatrix(&Matrix());
  assert (Graph_ != 0);

  if (PartitionerType_ == "linear")
    Partitioner_ = new Ifpack_LinearPartitioner(Graph_);
  else if (PartitionerType_ == "greedy")
    Partitioner_ = new Ifpack_GreedyPartitioner(Graph_);
  else if (PartitionerType_ == "metis")
    Partitioner_ = new Ifpack_METISPartitioner(Graph_);
  else if (PartitionerType_ == "equation")
    Partitioner_ = new Ifpack_EquationPartitioner(Graph_);
  else if (PartitionerType_ == "user")
    Partitioner_ = new Ifpack_UserPartitioner(Graph_);
  else
    IFPACK_CHK_ERR(-1);

  assert (Partitioner_ != 0);

  // need to partition the graph of A
  IFPACK_CHK_ERR(Partitioner_->SetParameters(List_));
  IFPACK_CHK_ERR(Partitioner_->Compute());

  // get actual number of partitions
  NumLocalBlocks_ = Partitioner_->NumLocalParts();

  // weight of each vertex
  if (W_)
    delete W_;
  W_ = new Epetra_Vector(Matrix().RowMatrixRowMap());
  W_->PutScalar(0.0);

  for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

    for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
      int LID = (*Partitioner_)(i,j);
      (*W_)[LID]++;
    }
  }
  W_->Reciprocal(*W_);

  IsInitialized_ = true;
  ++NumInitialize_;

  return(0);
}

//==============================================================================
#endif // HAVE_IFPACK_TEUCHOS
#endif // IFPACK_BLOCKPRECONDITIONER_H
