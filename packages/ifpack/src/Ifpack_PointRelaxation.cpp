#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include <iomanip>
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_Import.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Condest.h"

static const int IFPACK_JACOBI = 0;
static const int IFPACK_GS = 1;
static const int IFPACK_SGS = 2;

//==============================================================================
Ifpack_PointRelaxation::
Ifpack_PointRelaxation(const Epetra_RowMatrix* Matrix) :
  IsInitialized_(false),
  IsComputed_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  NumSweeps_(1),
  DampingFactor_(1.0),
  UseTranspose_(false),
  Condest_(-1.0),
  ComputeCondest_(false),
  PrecType_(IFPACK_JACOBI),
  MinDiagonalValue_(0.0),
  NumMyRows_(0),
  NumMyNonzeros_(0),
  NumGlobalRows_(0),
  NumGlobalNonzeros_(0),
  Matrix_(Matrix),
  Importer_(0),
  Diagonal_(0),
  Time_(0),
  IsParallel_(false),
  ZeroStartingSolution_(true)
{
}

//==============================================================================
Ifpack_PointRelaxation::~Ifpack_PointRelaxation()
{
  if (Diagonal_)
    delete Diagonal_;
  if (Time_)
    delete Time_;
  if (Importer_)
    delete Importer_;
}

//==============================================================================
int Ifpack_PointRelaxation::SetParameters(Teuchos::ParameterList& List)
{

  string PT;
  if (PrecType_ == IFPACK_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == IFPACK_GS)
    PT = "Gauss-Seidel";
  else if (PrecType_ == IFPACK_SGS)
    PT = "symmetric Gauss-Seidel";

  PT = List.get("relaxation: type", PT);

  if (PT == "Jacobi")
    PrecType_ = IFPACK_JACOBI;
  else if (PT == "Gauss-Seidel")
    PrecType_ = IFPACK_GS;
  else if (PT == "symmetric Gauss-Seidel")
    PrecType_ = IFPACK_SGS;
  else {
    IFPACK_CHK_ERR(-2);
  }
  
  NumSweeps_            = List.get("relaxation: sweeps",NumSweeps_);
  DampingFactor_        = List.get("relaxation: damping factor", 
                                   DampingFactor_);
  MinDiagonalValue_     = List.get("relaxation: min diagonal value", 
                                   MinDiagonalValue_);
  ZeroStartingSolution_ = List.get("relaxation: zero starting solution", 
                                   ZeroStartingSolution_);

  SetLabel();

  return(0);
}

//==============================================================================
const Epetra_Comm& Ifpack_PointRelaxation::Comm() const
{
  return(Matrix_->Comm());
}

//==============================================================================
const Epetra_Map& Ifpack_PointRelaxation::OperatorDomainMap() const
{
  return(Matrix_->OperatorDomainMap());
}

//==============================================================================
const Epetra_Map& Ifpack_PointRelaxation::OperatorRangeMap() const
{
  return(Matrix_->OperatorRangeMap());
}

//==============================================================================
int Ifpack_PointRelaxation::
Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  if (IsComputed() == false)
    IFPACK_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  IFPACK_CHK_ERR(Matrix_->Multiply(UseTranspose(),X,Y));
  return(0);
}

//==============================================================================
int Ifpack_PointRelaxation::Initialize()
{
  IsInitialized_ = false;

  if (Matrix_ == 0)
    IFPACK_CHK_ERR(-2);

  if (Time_ == 0)
    Time_ = new Epetra_Time(Comm());

  if (Matrix().NumGlobalRows() != Matrix().NumGlobalCols())
    IFPACK_CHK_ERR(-2); // only square matrices

  NumMyRows_ = Matrix_->NumMyRows();
  NumMyNonzeros_ = Matrix_->NumMyNonzeros();
  NumGlobalRows_ = Matrix_->NumGlobalRows();
  NumGlobalNonzeros_ = Matrix_->NumGlobalNonzeros();

  if (Comm().NumProc() != 1)
    IsParallel_ = true;
  else
    IsParallel_ = false;

  ++NumInitialize_;
  InitializeTime_ += Time_->ElapsedTime();
  IsInitialized_ = true;
  return(0);
}

//==============================================================================
int Ifpack_PointRelaxation::Compute()
{
  if (!IsInitialized())
    IFPACK_CHK_ERR(Initialize());

  Time_->ResetStartTime();

  // reset values
  IsComputed_ = false;
  Condest_ = -1.0;

  if (NumSweeps_ <= 0)
    IFPACK_CHK_ERR(-2); // at least one application
  
  if (Diagonal_)
    delete Diagonal_;
  Diagonal_ = new Epetra_Vector(Matrix().RowMatrixRowMap());

  if (Diagonal_ == 0)
    IFPACK_CHK_ERR(-5);

  IFPACK_CHK_ERR(Matrix().ExtractDiagonalCopy(*Diagonal_));

  // check diagonal elements, replace zeros with 1.0
  for (int i = 0 ; i < NumMyRows_ ; ++i) {
    double diag = (*Diagonal_)[i];
    if (IFPACK_ABS(diag) < MinDiagonalValue_)
      (*Diagonal_)[i] = MinDiagonalValue_;
  }

  // some methods require the inverse of the diagonal, compute it
  // now and store it in `Diagonal_'
  if ((PrecType_ == IFPACK_JACOBI) || (PrecType_ == IFPACK_GS)) {
    Diagonal_->Reciprocal(*Diagonal_);
    // update flops
    ComputeFlops_ += NumMyRows_;
  }

  // We need to import data from external processors. Here I create an
  // Epetra_Import object because I cannot assume that Matrix_ has one.
  // This is a bit of waste of resources (but the code is more robust).
  // Note that I am doing some strange stuff to set the components of Y
  // from Y2 (to save some time).
  //
  if (IsParallel_ && ((PrecType_ == IFPACK_GS) || (PrecType_ == IFPACK_SGS))) {
    Importer_ = new Epetra_Import(Matrix().RowMatrixColMap(),
                                  Matrix().RowMatrixRowMap());
    if (Importer_ == 0) IFPACK_CHK_ERR(-5);
  }

  ++NumCompute_;
  ComputeTime_ += Time_->ElapsedTime();
  IsComputed_ = true;

  return(0);
}

//==============================================================================
ostream& Ifpack_PointRelaxation::Print(ostream & os) const
{

  double MyMinVal, MyMaxVal;
  double MinVal, MaxVal;

  if (IsComputed_) {
    Diagonal_->MinValue(&MyMinVal);
    Diagonal_->MaxValue(&MyMaxVal);
    Comm().MinAll(&MyMinVal,&MinVal,1);
    Comm().MinAll(&MyMaxVal,&MaxVal,1);
  }

  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_PointRelaxation" << endl;
    os << "Sweeps         = " << NumSweeps_ << endl;
    os << "damping factor = " << DampingFactor_;
    if (PrecType_ == IFPACK_JACOBI)
      os << "Type           = Jacobi" << endl;
    else if (PrecType_ == IFPACK_GS)
      os << "Type           = Gauss-Seidel" << endl;
    else if (PrecType_ == IFPACK_SGS)
      os << "Type           = symmetric Gauss-Seidel" << endl;
    if (ZeroStartingSolution_) 
      os << "Using zero starting solution" << endl;
    else
      os << "Using input starting solution" << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << "Global number of rows            = " << Matrix_->NumGlobalRows() << endl;
    if (IsComputed_) {
      os << "Minimum value on stored diagonal = " << MinVal << endl;
      os << "Maximum value on stored diagonal = " << MaxVal << endl;
    }
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize_ 
       << "  " << std::setw(15) << InitializeTime_ 
       << "              0.0              0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute_ 
       << "  " << std::setw(15) << ComputeTime_
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_;
    if (ComputeTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_ / ComputeTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse_ 
       << "  " << std::setw(15) << ApplyInverseTime_
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_;
    if (ApplyInverseTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_ / ApplyInverseTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }

  return(os);
}

//==============================================================================
double Ifpack_PointRelaxation::
Condest(const Ifpack_CondestType CT, 
        const int MaxIters, const double Tol,
	Epetra_RowMatrix* Matrix)
{
  if (!IsComputed()) // cannot compute right now
    return(-1.0);

  // always computes it. Call Condest() with no parameters to get
  // the previous estimate.
  Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix);

  return(Condest_);
}

//==============================================================================
void Ifpack_PointRelaxation::SetLabel()
{
  string PT;
  if (PrecType_ == IFPACK_JACOBI)
    PT = "Jacobi";
  else if (PrecType_ == IFPACK_GS)
    PT = "GS";
  else if (PrecType_ == IFPACK_SGS)
    PT = "SGS";

  Label_ = "IFPACK (" + PT + ", sweeps=" + Ifpack_toString(NumSweeps_)
    + ", damping=" + Ifpack_toString(DampingFactor_) + ")";
}

//==============================================================================
// Note that Ifpack_PointRelaxation and Jacobi is much faster than
// Ifpack_AdditiveSchwarz<Ifpack_PointRelaxation> (because of the
// way the matrix-vector produce is performed).
//
// Another ML-related observation is that the starting solution (in Y)
// is NOT supposed to be zero. This may slow down the application of just
// one sweep of Jacobi.
//
int Ifpack_PointRelaxation::
ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (!IsComputed())
    IFPACK_CHK_ERR(-3);

  if (X.NumVectors() != Y.NumVectors())
    IFPACK_CHK_ERR(-2);

  Time_->ResetStartTime();

  // AztecOO gives X and Y pointing to the same memory location,
  // need to create an auxiliary vector, Xcopy
  const Epetra_MultiVector* Xcopy;
  if (X.Pointers()[0] == Y.Pointers()[0])
    Xcopy = new Epetra_MultiVector(X);
  else
    Xcopy = &X;
    
  if (ZeroStartingSolution_)
    Y.PutScalar(0.0);

  // Flops are updated in each of the following. 
  switch (PrecType_) {
  case IFPACK_JACOBI:
    IFPACK_CHK_ERR(ApplyInverseJacobi(*Xcopy,Y));
    break;
  case IFPACK_GS:
    IFPACK_CHK_ERR(ApplyInverseGS(*Xcopy,Y));
    break;
  case IFPACK_SGS:
    IFPACK_CHK_ERR(ApplyInverseSGS(*Xcopy,Y));
    break;
  default:
    IFPACK_CHK_ERR(-1); // something wrong
  }

  if (Xcopy != &X)
    delete Xcopy;

  ++NumApplyInverse_;
  ApplyInverseTime_ += Time_->ElapsedTime();
  return(0);
}

//==============================================================================
// This preconditioner can be much slower than AztecOO and ML versions
// if the matrix-vector product requires all ExtractMyRowCopy() 
// (as done through Ifpack_AdditiveSchwarz).
int Ifpack_PointRelaxation::
ApplyInverseJacobi(const Epetra_MultiVector& RHS, Epetra_MultiVector& LHS) const
{

  int NumVectors = LHS.NumVectors();
  Epetra_MultiVector* A_times_LHS;
  A_times_LHS = new Epetra_MultiVector(LHS.Map(),NumVectors);
  if (A_times_LHS == 0) IFPACK_CHK_ERR(-5);

  for (int j = 0; j < NumSweeps_ ; j++) {

    IFPACK_CHK_ERR(Apply(LHS,*A_times_LHS));
    IFPACK_CHK_ERR(A_times_LHS->Update(1.0,RHS,-1.0));
    IFPACK_CHK_ERR(LHS.Multiply(DampingFactor_, *A_times_LHS, *Diagonal_, 1.0));

  }
  delete A_times_LHS;

  // Flops:
  // - matrix vector              (2 * NumGlobalNonzeros_)
  // - update                     (2 * NumGlobalRows_)
  // - Multiply:
  //   - DampingFactor            (NumGlobalRows_)
  //   - Diagonal                 (NumGlobalRows_)
  //   - A + B                    (NumGlobalRows_)
  //   - 1.0                      (NumGlobalRows_)
  ApplyInverseFlops_ += NumVectors * (6 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);

  return(0);

}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int NumVectors = X.NumVectors();

  int Length = Matrix().MaxNumEntries();
  vector<int> Indices(Length);
  vector<double> Values(Length);

  Epetra_MultiVector* Y2;
  if (IsParallel_)
    Y2 = new Epetra_MultiVector(Importer_->TargetMap(), NumVectors);
  else
    Y2 = &Y;
                        
  // extract views (for nicer and faster code)
  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);

  for (int j = 0; j < NumSweeps_ ; j++) {

    // data exchange is here, once per sweep
    if (IsParallel_)
      IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    // FIXME: do I really need this code below?
    if (NumVectors == 1) {

      double* y0_ptr = y_ptr[0];
      double* y20_ptr = y2_ptr[0];
      double* x0_ptr = x_ptr[0];

      for (int i = 0 ; i < NumMyRows_ ; ++i) {

        int NumEntries;
        int col;
        IFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                                 &Values[0], &Indices[0]));

        double dtemp = 0.0;
        for (int k = 0 ; k < NumEntries ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y20_ptr[col];
        }

        y20_ptr[i] += DampingFactor_ * d_ptr[i] * (x0_ptr[i] - dtemp);
      }
      // using Export() sounded quite expensive
      if (IsParallel_)
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y0_ptr[i] = y20_ptr[i];

    }
    else {

      for (int i = 0 ; i < NumMyRows_ ; ++i) {

        int NumEntries;
        int col;
        IFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                                 &Values[0], &Indices[0]));

        for (int m = 0 ; m < NumVectors ; ++m) {

          double dtemp = 0.0;
          for (int k = 0 ; k < NumEntries ; ++k) {

            col = Indices[k];
            dtemp += Values[k] * y2_ptr[m][col];
          }

          y2_ptr[m][i] += DampingFactor_ * d_ptr[i] * (x_ptr[m][i] - dtemp);
        }
        // using Export() sounded quite expensive
      }

      if (IsParallel_)
        for (int m = 0 ; m < NumVectors ; ++m) 
          for (int i = 0 ; i < NumMyRows_ ; ++i)
            y_ptr[m][i] = y2_ptr[m][i];

    }
  }

  if (IsParallel_)
    delete Y2;

  ApplyInverseFlops_ += NumVectors * (4 * NumGlobalRows_ + 2 * NumGlobalNonzeros_);

  return(0);
}

//==============================================================================
int Ifpack_PointRelaxation::
ApplyInverseSGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  int NumVectors = X.NumVectors();
  int Length = Matrix().MaxNumEntries();
  vector<int> Indices(Length);
  vector<double> Values(Length);

  Epetra_MultiVector* Y2;
  if (IsParallel_) {
    Y2 = new Epetra_MultiVector(Importer_->TargetMap(), NumVectors);
  }
  else
    Y2 = &Y;

  double** y_ptr, ** y2_ptr, ** x_ptr, *d_ptr;
  X.ExtractView(&x_ptr);
  Y.ExtractView(&y_ptr);
  Y2->ExtractView(&y2_ptr);
  Diagonal_->ExtractView(&d_ptr);
  
  for (int iter = 0 ; iter < NumSweeps_ ; ++iter) {
    
    // only one data exchange per sweep
    if (IsParallel_)
      IFPACK_CHK_ERR(Y2->Import(Y,*Importer_,Insert));

    for (int i = 0 ; i < NumMyRows_ ; ++i) {

      int NumEntries;
      int col;
      double diag = d_ptr[i];
      double dtemp = 0.0;

      IFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

        for (int k = 0 ; k < NumEntries ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) / diag;
      }
    }

    for (int i = NumMyRows_  - 1 ; i > -1 ; --i) {

      int NumEntries;
      int col;
      double diag = d_ptr[i];
      double dtemp = 0.0;

      IFPACK_CHK_ERR(Matrix_->ExtractMyRowCopy(i, Length,NumEntries,
                                               &Values[0], &Indices[0]));

      for (int m = 0 ; m < NumVectors ; ++m) {

        for (int k = 0 ; k < NumEntries ; ++k) {

          col = Indices[k];
          dtemp += Values[k] * y2_ptr[m][col];
        }

        y2_ptr[m][i] += DampingFactor_ * (x_ptr[m][i] - dtemp) / diag;

      }
    }

    if (IsParallel_)
      for (int m = 0 ; m < NumVectors ; ++m) 
        for (int i = 0 ; i < NumMyRows_ ; ++i)
          y_ptr[m][i] = y2_ptr[m][i];
  }

  if (IsParallel_)
    delete Y2;

  ApplyInverseFlops_ += NumVectors * (8 * NumGlobalRows_ + 4 * NumGlobalNonzeros_);
  return(0);
}

#endif
