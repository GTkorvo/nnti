// @HEADER
// ***********************************************************************
// 
//                Amesos: Direct Sparse Solver Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Amesos_Umfpack.h"
extern "C" {
#include "umfpack.h"
}
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"

//=============================================================================
Amesos_Umfpack::Amesos_Umfpack(const Epetra_LinearProblem &prob ) :
  Symbolic(0),
  Numeric(0),
  SerialMatrix_(0), 
  UseTranspose_(false),
  Problem_(&prob), 
  Rcond_(0.0), 
  RcondValidOnAllProcs_(true)
{
  
  // MS // move declaration of Problem_ above because I need it
  // MS // set up before calling Comm()
  Teuchos::ParameterList ParamList ;
  SetParameters( ParamList ) ; 
}

//=============================================================================
Amesos_Umfpack::~Amesos_Umfpack(void) 
{
  if (Symbolic) umfpack_di_free_symbolic (&Symbolic);
  if (Numeric) umfpack_di_free_numeric (&Numeric);

  // print out some information if required by the user
  if ((verbose_ && PrintTiming_) || verbose_ == 2) PrintTiming();
  if ((verbose_ && PrintStatus_) || verbose_ == 2) PrintStatus();
}

//=============================================================================
// If FirstTime is true, then build SerialMap and ImportToSerial,
// otherwise simply re-ship the matrix, so that the numerical values
// are updated.
int Amesos_Umfpack::ConvertToSerial(const bool FirstTime) 
{ 
  ResetTime(0);
  ResetTime(1);
  
  const Epetra_Map &OriginalMap = Matrix()->RowMatrixRowMap() ; 

  NumGlobalElements_ = Matrix()->NumGlobalRows();
  numentries_ = Matrix()->NumGlobalNonzeros();
  assert (NumGlobalElements_ == Matrix()->NumGlobalCols());

  int NumMyElements_ = 0 ;
  if (MyPID_ == 0) NumMyElements_ = NumGlobalElements_;

  IsLocal_ = ( OriginalMap.NumMyElements() == 
	       OriginalMap.NumGlobalElements() )?1:0;
  Comm().Broadcast( &IsLocal_, 1, 0 ) ; 

  //  Convert Original Matrix to Serial (if it is not already) 
  //
  if (IsLocal_== 1) {
     SerialMatrix_ = Matrix();
  } 
  else 
  {
    if (FirstTime)
    {
      SerialMap_ = rcp(new Epetra_Map(NumGlobalElements_,NumMyElements_,
                                      0,Comm()));

      if (SerialMap_.get() == 0)
        AMESOS_CHK_ERR(-1);

      ImportToSerial_ = rcp(new Epetra_Import (SerialMap(),OriginalMap));

      if (ImportToSerial_.get() == 0)
        AMESOS_CHK_ERR(-1);
    }

    SerialCrsMatrixA_ = rcp(new Epetra_CrsMatrix(Copy,SerialMap(),0));

    if (SerialCrsMatrixA_.get() == 0)
      AMESOS_CHK_ERR(-1);

    SerialCrsMatrix().Import(*Matrix(), Importer(),Insert); 
    
    SerialCrsMatrix().FillComplete(); 
    SerialMatrix_ = &SerialCrsMatrix();
  }

  AddTime("matrix redistribution", 0);
  AddTime("overhead", 1);
  
  return(0);
} 

//=============================================================================
int Amesos_Umfpack::ConvertToUmfpackCRS()
{
  ResetTime(0);
  ResetTime(1);
  
  // Convert matrix to the form that Umfpack expects (Ap, Ai, Aval),
  // only on processor 0. The matrix has already been assembled in
  // SerialMatrix_; if only one processor is used, then SerialMatrix_
  // points to the problem's matrix.

  if (MyPID_ == 0) 
  {
    Ap.resize( NumGlobalElements_+1 );
    Ai.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ; 
    Aval.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ; 

    int NumEntries = SerialMatrix_->MaxNumEntries();

    int NumEntriesThisRow;
    int Ai_index = 0 ; 
    int MyRow;
    for (MyRow = 0 ; MyRow < NumGlobalElements_; MyRow++) 
    {
      int ierr;
      Ap[MyRow] = Ai_index ; 
      ierr = SerialMatrix_->ExtractMyRowCopy(MyRow, NumEntries, 
					     NumEntriesThisRow, 
					     &Aval[Ai_index], &Ai[Ai_index]);
      if (ierr)
	AMESOS_CHK_ERR(-1);

      Ai_index += NumEntriesThisRow;
    }

    Ap[MyRow] = Ai_index ; 
  }

  AddTime("conversion", 0);
  AddTime("overhead", 1);
  
  return 0;
}   

//=============================================================================
int Amesos_Umfpack::SetParameters( const Teuchos::ParameterList &ParameterList ) 
{
  // ========================================= //
  // retrive UMFPACK's parameters from list.   //
  // default values defined in the constructor //
  // ========================================= //

  // retrive general parameters
  SetStatusParameters( ParameterList ) ;
  SetControlParameters( ParameterList ) ;

  return 0;
}

//=============================================================================
int Amesos_Umfpack::PerformSymbolicFactorization() 
{
  // MS // no overhead time in this method
  ResetTime(0);  
  
  double *Control = (double *) NULL, *Info = (double *) NULL;
  
  if (Symbolic) 
    umfpack_di_free_symbolic (&Symbolic) ;
  if (MyPID_== 0) {
    (void) umfpack_di_symbolic (NumGlobalElements_, NumGlobalElements_, &Ap[0], 
				&Ai[0], &Aval[0], 
				&Symbolic, Control, Info) ;
  }

  AddTime("symbolic", 0);

  return 0;
}

//=============================================================================
int Amesos_Umfpack::PerformNumericFactorization( ) 
{
  // MS // no overhead time in this method
  ResetTime(0);

  RcondValidOnAllProcs_ = false ; 
  if (MyPID_ == 0) {
    vector<double> Control(UMFPACK_CONTROL);
    vector<double> Info(UMFPACK_INFO);
    umfpack_di_defaults( &Control[0] ) ; 
    if (Numeric) umfpack_di_free_numeric (&Numeric) ;
    int status = umfpack_di_numeric (&Ap[0], 
				     &Ai[0], 
				     &Aval[0], 
				     Symbolic, 
				     &Numeric, 
				     &Control[0], 
				     &Info[0]) ;
    Rcond_ = Info[UMFPACK_RCOND]; 

#if NOT_DEF
    cout << " Rcond_ = " << Rcond_ << endl ; 

    int lnz1 = 1000 ;
    int unz1 = 1000 ;
    int n = 4;
    int * Lp = (int *) malloc ((n+1) * sizeof (int)) ;
    int * Lj = (int *) malloc (lnz1 * sizeof (int)) ;
    double * Lx = (double *) malloc (lnz1 * sizeof (double)) ;
    int * Up = (int *) malloc ((n+1) * sizeof (int)) ;
    int * Ui = (int *) malloc (unz1 * sizeof (int)) ;
    double * Ux = (double *) malloc (unz1 * sizeof (double)) ;
    int * P = (int *) malloc (n * sizeof (int)) ;
    int * Q = (int *) malloc (n * sizeof (int)) ;
    double * Dx = (double *) NULL ;	/* D vector not requested */
    double * Rs  = (double *) malloc (n * sizeof (double)) ;
    if (!Lp || !Lj || !Lx || !Up || !Ui || !Ux || !P || !Q || !Rs)
    {
      assert( false ) ; 
    }
    int do_recip;
    status = umfpack_di_get_numeric (Lp, Lj, Lx, Up, Ui, Ux,
	P, Q, Dx, &do_recip, Rs, Numeric) ;
    if (status < 0)
    {
      assert( false ) ; 
    }

    printf ("\nL (lower triangular factor of C): ") ;
    (void) umfpack_di_report_matrix (n, n, Lp, Lj, Lx, 0, &Control[0]) ;
    printf ("\nU (upper triangular factor of C): ") ;
    (void) umfpack_di_report_matrix (n, n, Up, Ui, Ux, 1, &Control[0]) ;
    printf ("\nP: ") ;
    (void) umfpack_di_report_perm (n, P, &Control[0]) ;
    printf ("\nQ: ") ;
    (void) umfpack_di_report_perm (n, Q, &Control[0]) ;
    printf ("\nScale factors: row i of A is to be ") ;

#endif

    assert( status == 0 ) ; 
  }
  
  AddTime("numeric", 0);

  return 0;
}

//=============================================================================
double Amesos_Umfpack::GetRcond() const 
{
  if ( !RcondValidOnAllProcs_ ) {
    Comm().Broadcast( &Rcond_, 1, 0 ) ; 
    RcondValidOnAllProcs_ = true; 
  }
  return(Rcond_);
}; 

//=============================================================================
bool Amesos_Umfpack::MatrixShapeOK() const 
{ 
  bool OK = true;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


//=============================================================================
int Amesos_Umfpack::SymbolicFactorization() 
{
  // MS // NOTE: If you change this method, also change
  // MS // NumericFactorization(), because it performs part of the actions
  // MS // of this method. This is to avoid to ship the matrix twice
  // MS // (once for the symbolic factorization, and once for the numeric
  // MS // factorization) when it is not necessary.
  
  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;

  InitTime(Comm(), 2);
  // MS // Initialize two timers:
  // MS // timer 1: this will track all time spent in Amesos most
  // MS //          important functions, *including* UMFPACK functions
  // MS // timer 2: this will track all time spent in this function
  // MS //          that is not due to UMFPACK calls, and therefore it
  // MS //          tracks down how much Amesos costs. The timer starts
  // MS //          and ends in *each* method, unless the method does not
  // MS //          perform any real operation. If a method calls another
  // MS //          method, the timer will be stopped before the called
  // MS //          method, then restared.
  // MS //          All the time of this timer goes into "overhead"

  MyPID_    = Comm().MyPID();
  NumProcs_ = Comm().NumProc();

  AMESOS_CHK_ERR(ConvertToSerial(true)); 
  AMESOS_CHK_ERR(ConvertToUmfpackCRS());
  
  AMESOS_CHK_ERR(PerformSymbolicFactorization());

  IsSymbolicFactorizationOK_ = true; 

  NumSymbolicFact_++;  

  return 0;
}

//=============================================================================
int Amesos_Umfpack::NumericFactorization() 
{
  IsNumericFactorizationOK_ = false;

  if (IsSymbolicFactorizationOK_ == false)
  {
    // Call here what is needed, to avoid double shipping of the matrix
    InitTime(Comm(), 2);

    MyPID_    = Comm().MyPID();
    NumProcs_ = Comm().NumProc();

    AMESOS_CHK_ERR(ConvertToSerial(true)); 
    AMESOS_CHK_ERR(ConvertToUmfpackCRS());
  
    AMESOS_CHK_ERR(PerformSymbolicFactorization());

    IsSymbolicFactorizationOK_ = true; 

    NumSymbolicFact_++;  

    AMESOS_CHK_ERR(PerformNumericFactorization());
  }
  else
  {
    // need to reshuffle and reconvert because entry values may have changed
    AMESOS_CHK_ERR(ConvertToSerial(false));
    AMESOS_CHK_ERR(ConvertToUmfpackCRS());

    AMESOS_CHK_ERR(PerformNumericFactorization());
  }

  NumNumericFact_++;  

  IsNumericFactorizationOK_ = true;

  return 0;
}

//=============================================================================
int Amesos_Umfpack::Solve() 
{ 
  // if necessary, perform numeric factorization. 
  // This may call SymbolicFactorization() as well.
  if (!IsNumericFactorizationOK_)
    AMESOS_CHK_ERR(NumericFactorization()); 

  ResetTime(1);

  Epetra_MultiVector* vecX = Problem_->GetLHS(); 
  Epetra_MultiVector* vecB = Problem_->GetRHS(); 

  if ((vecX == 0) || (vecB == 0))
    AMESOS_CHK_ERR(-1);

  int NumVectors = vecX->NumVectors(); 
  if (NumVectors != vecB->NumVectors())
    AMESOS_CHK_ERR(-1);

  Epetra_MultiVector *SerialB, *SerialX; 

  //  Extract Serial versions of X and B 
  //
  double *SerialXvalues ;
  double *SerialBvalues ;

  Epetra_MultiVector* SerialXextract = 0;
  Epetra_MultiVector* SerialBextract = 0;
    
  //  Copy B to the serial version of B
  //
  ResetTime(0);
  
  if (IsLocal_ == 1) { 
    SerialB = vecB ; 
    SerialX = vecX ; 
  } else { 
    assert (IsLocal_ == 0);
    SerialXextract = new Epetra_MultiVector(SerialMap(),NumVectors); 
    SerialBextract = new Epetra_MultiVector(SerialMap(),NumVectors); 

    SerialBextract->Import(*vecB,Importer(),Insert);
    SerialB = SerialBextract; 
    SerialX = SerialXextract; 
  } 

  AddTime("vector redistribution", 0);
  
  //  Call UMFPACK to perform the solve
  //  Note:  UMFPACK uses a Compressed Column Storage instead of compressed row storage, 
  //  Hence to compute A X = B, we ask UMFPACK to perform A^T X = B and vice versa

  AddTime("overhead", 1);

  ResetTime(0);

  int SerialBlda, SerialXlda ; 
  int UmfpackRequest = UseTranspose()?UMFPACK_A:UMFPACK_At ;
  int status = 0;

  if ( MyPID_ == 0 ) {
    int ierr;
    ierr = SerialB->ExtractView(&SerialBvalues, &SerialBlda);
    assert (ierr == 0);
    ierr = SerialX->ExtractView(&SerialXvalues, &SerialXlda);
    assert (ierr == 0);
    assert( SerialBlda == NumGlobalElements_ ) ; 
    assert( SerialXlda == NumGlobalElements_ ) ; 
    
    for ( int j =0 ; j < NumVectors; j++ ) { 
      double *Control = (double *) NULL, *Info = (double *) NULL ;

      status = umfpack_di_solve (UmfpackRequest, &Ap[0], 
				     &Ai[0], &Aval[0], 
				     &SerialXvalues[j*SerialXlda], 
				     &SerialBvalues[j*SerialBlda], 
				     Numeric, Control, Info) ;
    }
  }
    
  if (status) AMESOS_CHK_ERR(status);

  AddTime("solve", 0);
  
  //  Copy X back to the original vector
  
  ResetTime(0);
  ResetTime(1);

  if ( IsLocal_ == 0 ) {
    vecX->Export(*SerialX, Importer(), Insert ) ;
    if (SerialBextract) delete SerialBextract ;
    if (SerialXextract) delete SerialXextract ;
  }

  AddTime("vector redistribution", 0);

  if (ComputeTrueResidual_)
  {
    Epetra_RowMatrix* Matrix = 
      dynamic_cast<Epetra_RowMatrix*>(Problem_->GetOperator());
    ComputeTrueResidual(*Matrix, *vecX, *vecB, UseTranspose(), "Amesos_Umfpack");
  }

  if (ComputeVectorNorms_) {
    ComputeVectorNorms(*vecX, *vecB, "Amesos_Umfpack");
  }

  NumSolve_++;

  AddTime("overhead", 1); // Amesos overhead

  return(0);
}

// ====================================================================== 
void Amesos_Umfpack::PrintStatus() const
{
  if (MyPID_ != 0) return;

  PrintLine();

  cout << "Amesos_Umfpack : Matrix has " << NumGlobalElements_ << " rows"
       << " and " << numentries_ << " nonzeros" << endl;
  cout << "Amesos_Umfpack : Nonzero elements per row = "
       << 1.0*numentries_/NumGlobalElements_ << endl;
  cout << "Amesos_Umfpack : Percentage of nonzero elements = "
       << 100.0*numentries_/(pow(NumGlobalElements_,2.0)) << endl;
  cout << "Amesos_Umfpack : Use transpose = " << UseTranspose_ << endl;

  PrintLine();

  return;
}

// ====================================================================== 
void Amesos_Umfpack::PrintTiming() const
{
  if (Problem_->GetOperator() == 0 || MyPID_ != 0)
    return;

  double ConTime = GetTime("conversion");
  double MatTime = GetTime("matrix redistribution");
  double VecTime = GetTime("vector redistribution");
  double SymTime = GetTime("symbolic");
  double NumTime = GetTime("numeric");
  double SolTime = GetTime("solve");
  double OveTime = GetTime("overhead");

  if (NumSymbolicFact_)
    SymTime /= NumSymbolicFact_;

  if (NumNumericFact_)
    NumTime /= NumNumericFact_;

  if (NumSolve_)
    SolTime /= NumSolve_;

  string p = "Amesos_Umfpack : ";
  PrintLine();

  cout << p << "Time to convert matrix to Umfpack format = "
       << ConTime << " (s)" << endl;
  cout << p << "Time to redistribute matrix = "
       << MatTime << " (s)" << endl;
  cout << p << "Time to redistribute vectors = "
       << VecTime << " (s)" << endl;
  cout << p << "Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << p << "Time for sym fact = "
       << SymTime * NumSymbolicFact_ << " (s), avg = " 
       << SymTime << " (s)" << endl;
  cout << p << "Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << p << "Time for num fact = "
       << NumTime * NumNumericFact_ << " (s), avg = " 
       << NumTime << " (s)" << endl;
  cout << p << "Number of solve phases = "
       << NumSolve_ << endl;
  cout << p << "Time for solve = "
       << SolTime * NumSolve_ << " (s), avg = " << SolTime << " (s)" << endl;
  double tt = SymTime * NumSymbolicFact_ + NumTime * NumNumericFact_ + SolTime * NumSolve_;
  if (tt != 0)
  {
    cout << p << "Total time spent in Amesos = " << tt << " (s) " << endl;
    cout << p << "Total time spent in the Amesos interface = " << OveTime << " (s)" << endl;
    cout << p << "(the above time does not include UMFPACK time)" << endl;
    cout << p << "Amesos interface time / total time = " << OveTime / tt << endl;
  }

  PrintLine();

  return;
}
