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

#include "Amesos_Klu.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "CrsMatrixTranspose.h"
extern "C" {
  // #include "amd.h"
#include "klu_btf.h"
#if 0
#include "klu_dump.h"
#endif
}

class Amesos_Klu_Pimpl {
public:
   klu_symbolic *Symbolic_ ;
   klu_numeric *Numeric_ ;

  Amesos_Klu_Pimpl():
    Symbolic_(0),
    Numeric_(0)
  {}

  ~Amesos_Klu_Pimpl(void){

    if ( Symbolic_ ) klu_btf_free_symbolic (&Symbolic_) ;
    if ( Numeric_ ) klu_btf_free_numeric (&Numeric_) ;
  }

} ;


//=============================================================================
Amesos_Klu::Amesos_Klu(const Epetra_LinearProblem &prob ) :
  PrivateKluData_( new Amesos_Klu_Pimpl() ),
  SerialMap_(0),
  SerialCrsMatrixA_(0),
  SerialMatrix_(0),
  Matrix_(0),
  UseTranspose_(false),
  Problem_(&prob),
  refactorize_(false),
  rcond_threshold_(1e-12),
  ScaleMethod_(1),
  ImportToSerial_(0)
{
  // MS // move declaration of Problem_ above because I need it
  // MS // set up before calling Comm()
  Teuchos::ParameterList ParamList ;
  SetParameters( ParamList ) ;

}

//=============================================================================
Amesos_Klu::~Amesos_Klu(void) {

  if (SerialMap_) 
    delete SerialMap_;
  if (SerialCrsMatrixA_) 
    delete SerialCrsMatrixA_;

  delete PrivateKluData_;

  if (ImportToSerial_) 
    delete ImportToSerial_;

  // print out some information if required by the user
  if( (verbose_ && PrintTiming_) || verbose_ == 2 ) PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2 ) PrintStatus();


}

//=============================================================================
int Amesos_Klu::ConvertToSerial() 
{

  ResetTime();

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());

  iam = Comm().MyPID() ;

  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ;

  NumGlobalElements_ = RowMatrixA->NumGlobalRows();
  numentries_ = RowMatrixA->NumGlobalNonzeros();
  assert( NumGlobalElements_ == RowMatrixA->NumGlobalCols() );

  //
  //  Create a serial matrix
  //
  assert( NumGlobalElements_ == OriginalMap.NumGlobalElements() ) ;
  int NumMyElements_ = 0 ;
  if (iam==0) NumMyElements_ = NumGlobalElements_;


  IsLocal_ = ( OriginalMap.NumMyElements() ==
	       OriginalMap.NumGlobalElements() )?1:0;
  Comm().Broadcast( &IsLocal_, 1, 0 ) ;

  //
  //  KEN:  Consider giving Epetra_RowMatrix_Transposer a shot
  //  I am not confident that  Epetra_RowMatrix_Transposer works,
  //  but it is worth a try.
  //
  //  Convert Original Matrix to Serial (if it is not already)
  //
  if (SerialMap_) { 
    delete SerialMap_ ; SerialMap_ = 0 ;
  }
  if (SerialCrsMatrixA_) { 
    delete SerialCrsMatrixA_ ; SerialCrsMatrixA_ = 0;
  }
  if (IsLocal_ == 1) {
     SerialMatrix_ = RowMatrixA;
  } else {
    if ( SerialMap_ ) delete SerialMap_;
    SerialMap_ = new Epetra_Map(NumGlobalElements_, NumMyElements_, 0, Comm());

    // check whether the stored ImportToSerial_ (if allocated) 
    // is still valid or not.
    if (ImportToSerial_ != 0) {
      if (!(OriginalMap.SameAs(ImportToSerial_->TargetMap()))) {
	delete SerialMap_;
	AMESOS_CHK_ERR(CreateSerialMap());
	delete ImportToSerial_;
	ImportToSerial_ = 0;
      }
    }

    if (ImportToSerial_ == 0) {
      ImportToSerial_ = new Epetra_Import(*SerialMap_,OriginalMap);
      assert (ImportToSerial_ != 0);
    }

    if (SerialCrsMatrixA_) 
      delete SerialCrsMatrixA_ ;
    SerialCrsMatrixA_ = new Epetra_CrsMatrix(Copy, *SerialMap_, 0);
    AMESOS_CHK_ERR(SerialCrsMatrixA_->Import(*RowMatrixA, 
					     *ImportToSerial_, Add));

    AMESOS_CHK_ERR(SerialCrsMatrixA_->FillComplete());
    SerialMatrix_ = SerialCrsMatrixA_ ;
  }

  AddTime("matrix redistribution");

  return 0;
}

//=============================================================================
int Amesos_Klu::CreateSerialMap()
{
  Epetra_RowMatrix *RowMatrixA;
  RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  if (RowMatrixA == 0)
    AMESOS_CHK_ERR(-1);

  iam = Comm().MyPID();
  NumGlobalElements_ = RowMatrixA->NumGlobalRows();
  numentries_ = RowMatrixA->NumGlobalNonzeros();
  int NumMyElements = 0;
  if (iam == 0) 
    NumMyElements = NumGlobalElements_;

  SerialMap_ = new Epetra_Map(NumGlobalElements_,NumMyElements,0,Comm());
  if (SerialMap_ == 0)
    AMESOS_CHK_ERR(-1);
  
  return(0);
}

//=============================================================================
//
//  See also pre and post conditions in Amesos_Klu.h
//  Preconditions:
//    firsttime specifies that this is the first time that 
//    ConertToKluCrs has been called, i.e. in symbolic factorization.  
//    No data allocation should happen unless firsttime=true.
//    SerialMatrix_ points to the matrix to be factored and solved
//    NumGlobalElements_ has been set to the dimension of the matrix
//    numentries_ has been set to the number of non-zeros in the matrix
//      (i.e. ConvertToSerial() has been callded)
//
//  Postconditions:
//    Ap, Ai, Aval contain the matrix as Klu needs it
//
//
int Amesos_Klu::ConvertToKluCRS(bool firsttime)
{

  ResetTime();

  Matrix_ = SerialMatrix_ ;
  //
  //  Convert matrix to the form that Klu expects (Ap, Ai, Aval)
  //

  if (iam==0) {
    assert( NumGlobalElements_ == Matrix_->NumGlobalRows());
    assert( NumGlobalElements_ == Matrix_->NumGlobalCols());
    assert( numentries_ == Matrix_->NumGlobalNonzeros());
    if ( firsttime ) { 
      Ap.resize( NumGlobalElements_+1 );
      Ai.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ;
      Aval.resize( EPETRA_MAX( NumGlobalElements_, numentries_) ) ;
    }

    int NumEntriesThisRow;
    int Ai_index = 0 ;
    int MyRow;
    Epetra_CrsMatrix *CrsMatrix = dynamic_cast<Epetra_CrsMatrix *>(Matrix_);
    //    Epetra_CrsMatrix *CrsMatrix = 0 ;  //  Uncomment this and comment the above line to test how we do with Row Matrices

    int MaxNumEntries_ = Matrix_->MaxNumEntries();
    if ( firsttime && CrsMatrix == 0 ) {
      ColIndicesV_.resize(MaxNumEntries_);
      RowValuesV_.resize(MaxNumEntries_);
    }
    double *RowValues;
    int *ColIndices;

    for ( MyRow = 0; MyRow <NumGlobalElements_; MyRow++ ) {
      if ( CrsMatrix != 0 ) {
	EPETRA_CHK_ERR( CrsMatrix->
			ExtractMyRowView( MyRow, NumEntriesThisRow, RowValues,
					  ColIndices ) != 0 ) ;
      } else {
	EPETRA_CHK_ERR( Matrix_->
			ExtractMyRowCopy( MyRow, MaxNumEntries_,
					  NumEntriesThisRow, &RowValuesV_[0],
					  &ColIndicesV_[0] ) != 0 ) ;
	RowValues =  &RowValuesV_[0];
	ColIndices = &ColIndicesV_[0];
      }

      if ( firsttime ) {
	Ap[MyRow] = Ai_index ;
	for ( int j = 0; j < NumEntriesThisRow; j++ ) {
	  Ai[Ai_index] = ColIndices[j] ;
	  Ai_index++;
	}
      } else { 
	for ( int j = 0; j < NumEntriesThisRow; j++ ) {
	  Aval[Ai_index] = RowValues[j] ;     
          if (ColIndices[j] == MyRow) {
            Aval[Ai_index] += AddToDiag_;     // Bug #1405   - this fails if the matrix is missing diagonal entries 
	  }
	  Ai_index++;
	}
      }
    }
    Ap[MyRow] = Ai_index ;
  }

  AddTime("matrix conversion");

  return 0;
}

//=============================================================================
int Amesos_Klu::SetParameters( Teuchos::ParameterList &ParameterList ) {

  // ========================================= //
  // retrive KLU's parameters from list.       //
  // default values defined in the constructor //
  // ========================================= //

  // retrive general parameters

  // solve problem with transpose
  if( ParameterList.isParameter("UseTranspose") )
    SetUseTranspose(ParameterList.get("UseTranspose",false));

  // print some timing information (on process 0)
  if( ParameterList.isParameter("PrintTiming") )
    PrintTiming_ = ParameterList.get("PrintTiming", false);

  // print some statistics (on process 0). Do not include timing
  if( ParameterList.isParameter("PrintStatus") )
    PrintStatus_ = ParameterList.get("PrintStatus", false);

  // add this value to diagonal
  if( ParameterList.isParameter("AddToDiag") )
    AddToDiag_ = ParameterList.get("AddToDiag", 0.0);

  // compute norms of some vectors
  if( ParameterList.isParameter("ComputeVectorNorms") )
    ComputeVectorNorms_ = ParameterList.get("ComputeVectorNorms",false);

  // compute the true residual Ax-b after solution
  if( ParameterList.isParameter("ComputeTrueResidual") )
    ComputeTrueResidual_ = ParameterList.get("ComputeTrueResidual",false);

  // some verbose output:
  // 0 - no output at all
  // 1 - output as specified by other parameters
  // 2 - all possible output
  if( ParameterList.isParameter("OutputLevel") )
    verbose_ = ParameterList.get("OutputLevel",1);

  // refactorize
  if( ParameterList.isParameter("Refactorize") )
    refactorize_ = ParameterList.get("Refactorize", false);

  // threshold for determining if refactorize worked OK
  // UNUSED at present - KSS June 2004
  if( ParameterList.isParameter("RcondThreshold") )
    rcond_threshold_ = ParameterList.get("RcondThreshold", 1e-12);

  // scaling method: 0: none, 1: use method's default, 2: use
  // the method's 1st alternative, 3: etc.
  if( ParameterList.isParameter("ScaleMethod") )
    ScaleMethod_ = ParameterList.get("ScaleMethod", 1);

  // MS // now comment it out, if we have parameters for KLU sublist
  // MS // uncomment it
  /*
  if (ParameterList.isSublist("Klu") ) {
    Teuchos::ParameterList KluParams = ParameterList.sublist("Klu") ;
  }
  */

  return 0;
}


//=============================================================================
int Amesos_Klu::PerformSymbolicFactorization() 
{
  ResetTime();

  if (iam == 0) {
    if (PrivateKluData_->Symbolic_) {
	klu_btf_free_symbolic (&(PrivateKluData_->Symbolic_)) ;
    }

    PrivateKluData_->Symbolic_ =
	klu_btf_analyze (NumGlobalElements_, &Ap[0], &Ai[0], (klu_control *) 0);
    if ( PrivateKluData_->Symbolic_ == 0 ) EPETRA_CHK_ERR( 1 ) ;
  }

  AddTime("symbolic");

  return 0;
}

//=============================================================================
int Amesos_Klu::PerformNumericFactorization( ) 
{
  ResetTime();

  if (iam == 0) {

    bool factor_with_pivoting = true ;

    // set the default parameters
    klu_control control ;
    klu_btf_defaults (&control) ;
    control.scale = ScaleMethod_ ;

    // see if we can "refactorize"
    if ( refactorize_ && PrivateKluData_->Numeric_ ) {

	// refactorize using the existing Symbolic and Numeric objects, and
	// using the identical pivot ordering as the prior klu_btf_factor.
	// No partial pivoting is done.
	int result = klu_btf_refactor (&Ap[0], &Ai[0], &Aval[0],
		    PrivateKluData_->Symbolic_, &control,
		    PrivateKluData_->Numeric_) ;

	// Did it work?
	if ( result == KLU_OK) {

	    // Get the largest and smallest entry on the diagonal of U
	    double umin = (PrivateKluData_->Numeric_)->umin ;
	    double umax = (PrivateKluData_->Numeric_)->umax ;

	    // compute a crude estimate of the reciprocal of
	    // the condition number
	    double rcond = umin / umax ;

	    if ( rcond > rcond_threshold_ ) {
		// factorizing without pivot worked fine.  We are done.
		factor_with_pivoting = false ;
	    }

	}
    }

    if ( factor_with_pivoting ) {

	// factor with partial pivoting:
	// either this is the first time we are factoring the matrix, or the
	// refactorize parameter is false, or we tried to refactorize and
	// found it to be too inaccurate.

	// destroy the existing Numeric object, if it exists
	if ( PrivateKluData_->Numeric_ ) {
	    klu_btf_free_numeric (&(PrivateKluData_->Numeric_)) ;
	}

	// factor the matrix using partial pivoting
	PrivateKluData_->Numeric_ =
	    klu_btf_factor (&Ap[0], &Ai[0], &Aval[0],
		    PrivateKluData_->Symbolic_, &control) ;
	if ( PrivateKluData_->Numeric_ == 0 ) EPETRA_CHK_ERR( 2 ) ;
    }

  }

  AddTime("numeric");

  return 0;
}

//=============================================================================
bool Amesos_Klu::MatrixShapeOK() const {
  bool OK = true;

  // Comment by Tim:  The following code seems suspect.  The variable "OK"
  // is not set if the condition is true.
  // Does the variable "OK" default to true?
  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() !=
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) {
    OK = false;
  }
  return OK;
}

//=============================================================================
int Amesos_Klu::SymbolicFactorization() 
{

  IsSymbolicFactorizationOK_ = false;
  IsNumericFactorizationOK_ = false;
  
  InitTime(Comm());

  NumSymbolicFact_++;

  ConvertToSerial() ;

  ConvertToKluCRS(true);

  PerformSymbolicFactorization();

  IsSymbolicFactorizationOK_ = true;
  
  return 0;
}

//=============================================================================
int Amesos_Klu::NumericFactorization() 
{
 
  IsNumericFactorizationOK_ = false;
  if (IsSymbolicFactorizationOK_ == false)
    AMESOS_CHK_ERR(SymbolicFactorization());

  NumNumericFact_++;

  ConvertToSerial() ;
  ConvertToKluCRS(false);

  PerformNumericFactorization();

  IsNumericFactorizationOK_ = true;
  
  return 0;
}

//=============================================================================
int Amesos_Klu::Solve() 
{

  if (IsNumericFactorizationOK_ == false)
    AMESOS_CHK_ERR(NumericFactorization());
  
  ++NumSolve_;

  Epetra_MultiVector* vecX = Problem_->GetLHS() ;
  Epetra_MultiVector* vecB = Problem_->GetRHS() ;

  if ((vecX == 0) || (vecB == 0))
    AMESOS_CHK_ERR(-1); // something wrong in input
  
  int NumVectors = vecX->NumVectors();
  if (NumVectors != vecB->NumVectors())
    AMESOS_CHK_ERR(-1); // something wrong in input

  // vectors with SerialMap_
  Epetra_MultiVector* SerialB = 0;
  Epetra_MultiVector* SerialX = 0;

  //  Extract Serial versions of X and B
  double *SerialXvalues ;

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  Epetra_MultiVector* SerialXextract = 0;
  Epetra_MultiVector* SerialBextract = 0;

  ResetTime();

  //  Copy B to the serial version of B
  //
  if (IsLocal_ == 1) {
    SerialB = vecB;
    SerialX = vecX;
  } else {
    assert (IsLocal_ == 0);
    const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap();

    // check whether the stored ImportToSerial_ (if allocated) 
    // is still valid or not.
    if (ImportToSerial_ != 0) {
      if (!(OriginalMap.SameAs(ImportToSerial_->TargetMap()))) {
	delete SerialMap_;
	AMESOS_CHK_ERR(CreateSerialMap());
	delete ImportToSerial_;
	ImportToSerial_ = 0;
      }
    }

    if (ImportToSerial_ == 0) {
      ImportToSerial_ = new Epetra_Import(*SerialMap_,OriginalMap);
      assert (ImportToSerial_ != 0);
    }

    assert ( SerialBextract == 0 ) ; 
    assert ( SerialXextract == 0 ) ; 

    SerialXextract = new Epetra_MultiVector(*SerialMap_,NumVectors);
    SerialBextract = new Epetra_MultiVector(*SerialMap_,NumVectors);
    
    SerialBextract->Import(*vecB,*ImportToSerial_,Insert);
    SerialB = SerialBextract ;
    SerialX = SerialXextract ;
  }

  AddTime("vector redistribution");

  //  Call KLU to perform the solve

  SerialX->Scale(1.0, *SerialB) ;

  ResetTime();

  int SerialXlda ;
  if (iam == 0) {
    AMESOS_CHK_ERR(SerialX->ExtractView(&SerialXvalues,&SerialXlda ));

    if (SerialXlda != NumGlobalElements_)
      AMESOS_CHK_ERR(-1);

    if (UseTranspose()) {
      klu_btf_solve( PrivateKluData_->Symbolic_, PrivateKluData_->Numeric_,
		     SerialXlda, NumVectors, &SerialXvalues[0] );
    } else {
      klu_btf_tsolve( PrivateKluData_->Symbolic_, PrivateKluData_->Numeric_,
		      SerialXlda, NumVectors, &SerialXvalues[0] );
    }
  }

  AddTime("solve");

  //  Copy X back to the original vector

  ResetTime();

  if (IsLocal_ == 0) {
    vecX->Export( *SerialX, *ImportToSerial_, Insert ) ;
    delete SerialBextract ;
    delete SerialXextract ;

  } // otherwise we are already in place.

  AddTime("vector redistribution");

#if 0
  //
  //  ComputeTrueResidual causes TestOptions to fail on my linux box 
  //  Bug #1147
  if (ComputeTrueResidual_)
    ComputeTrueResidual(*Matrix_, *vecX, *vecB, UseTranspose(), "Amesos_Klu");
#endif

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*vecX, *vecB, "Amesos_Klu");

  return(0) ;
}

// ================================================ ====== ==== ==== == =

void Amesos_Klu::PrintStatus() const
{

  if( iam != 0  ) return;

  PrintLine();

  cout << "Amesos_Klu : Matrix has " << NumGlobalElements_ << " rows"
       << " and " << numentries_ << " nonzeros" << endl;
  cout << "Amesos_Klu : Nonzero elements per row = "
       << 1.0*numentries_/NumGlobalElements_ << endl;
  cout << "Amesos_Klu : Percentage of nonzero elements = "
       << 100.0*numentries_/(pow(NumGlobalElements_,2.0)) << endl;
  cout << "Amesos_Klu : Use transpose = " << UseTranspose_ << endl;

  PrintLine();

  return;

}

// ================================================ ====== ==== ==== == =

void Amesos_Klu::PrintTiming() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  double ConTime = GetTime("matrix conversion");
  double MatTime = GetTime("matrix redistribution");
  double VecTime = GetTime("vector redistribution");
  double SymTime = GetTime("symbolic");
  double NumTime = GetTime("numeric");
  double SolTime = GetTime("solve");

  if (NumSymbolicFact_)
    SymTime /= NumSymbolicFact_;

  if (NumNumericFact_)
    NumTime /= NumNumericFact_;

  if (NumSolve_)
    SolTime /= NumSolve_;

  string p = "Amesos_Klu : ";
  PrintLine();

  cout << p << "Time to convert matrix to Klu format = "
       << ConTime << " (s)" << endl;
  cout << p << "Time to redistribute matrix = "
       << MatTime << " (s)" << endl;
  cout << p << "Time to redistribute vectors = "
       << VecTime << " (s)" << endl;
  cout << p << "Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << p << "Time for sym fact = "
       << SymTime << " (s), avg = " << SymTime << " (s)" << endl;
  cout << p << "Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << p << "Time for num fact = "
       << NumTime << " (s), avg = " << NumTime << " (s)" << endl;
  cout << p << "Number of solve phases = "
       << NumSolve_ << endl;
  cout << p << "Time for solve = "
       << SolTime << " (s), avg = " << SolTime << " (s)" << endl;

  PrintLine();

  return;
}
