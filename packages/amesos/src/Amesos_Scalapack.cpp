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

#include "Amesos_Scalapack.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "CrsMatrixTranspose.h"


//=============================================================================
Amesos_Scalapack::Amesos_Scalapack(const Epetra_LinearProblem &prob ):
  ictxt_(-1313),
  ScaLAPACK1DMap_(0), 
  ScaLAPACK1DMatrix_(0), 
  SerialMap_(0),
  UseTranspose_(false),
  Problem_(&prob), 
  PrintTiming_(false),
  PrintStatus_(false),
  ComputeVectorNorms_(false),
  ComputeTrueResidual_(false),
  verbose_(1),
  debug_(0),
  ConTime_(0.0),
  SymTime_(0.0),
  NumTime_(0.0),
  SolTime_(0.0),
  VecTime_(0.0),
  MatTime_(0.0),
  NumSymbolicFact_(0),
  NumNumericFact_(0),
  NumSolve_(0),
  Time(Comm())
{
  // MS // Ken, I moved Problem_ above because I need it to be set
  // MS // before calling Comm()
    
    Teuchos::ParameterList ParamList ;
    SetParameters( ParamList ) ; 
}

//=============================================================================
Amesos_Scalapack::~Amesos_Scalapack(void) {

  if ( ScaLAPACK1DMap_ ) delete ScaLAPACK1DMap_ ; 
  if ( ScaLAPACK1DMatrix_ ) delete ScaLAPACK1DMatrix_ ; 

  // print out some information if required by the user
  if( (verbose_ && PrintTiming_) || verbose_ == 2 ) PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2 ) PrintStatus();
  
}
//  See  pre and post conditions in Amesos_Scalapack.h

//
//  Distribution of the matrix:
//    The current code uses a data distribution which both ScaLAPACK
//    and Trilinos support, though ScaLAPACK sees the transpose.
//    Trilinos produces a blocked distributed matrix with an even 
//    number of rows owned by each process.  
//    ScaLAPACK factors the transpose of this matrix, hence NPROW=1,
//    NPCOL = NumberOfProcesses and NB (the block size) is set equal 
//    to the number of number of rows owned by each process.
//
//    This data distribution should be reasonably efficient for the matrices
//    that we expect to deal with, i.e. n on the order of 4,000 
//    with up to 16 processes.  For n >= 10,000 we should switch to 
//    a 2D block-cyclis data distribution.  
//
//    The number of messages in pdgetrf is O(n lg(nprow) + n/nb lg(npcol)).
//    Hence with nprow=1, the number of messages is O(n/nb) 
//    For small problems message latency is a significnat cost and hence this
//    is a reasonably efficient data layout.
//    The total volume of message traffic in pdgetrf is O(n^2/nprow + n^2/npcol)
//    The load imbalance is O(n^2 nb/nprow + n^2 nb/npcol) 
//    For medium sized problems, message throughput and load imbalance are 
//    important and hence this is a less efficient data layout for medium 
//    sized problems.  
//    
//  Distribution of the vector(s)
//    ScaLAPACK requires that the vectors be distributed in the same manner
//    as the matrices.  Since PDGETRF factors the transpose of the matrix, 
//    using NPROW=1, the vectors must also be distributed with NPROW=1, i.e.
//    each vector fully owned by a particular process.  And, the first nb
//    vectors belong on the first process. 
//
//    The easiest way to deal with this is to limit ourselves to nb right hand 
//    sides at a time (this is not a significant limitation as nb >= n/p ) 
//    and using our basic heuristic for the number of processes to use 
//    nb >= min(200,n) as well. 
//   
//    Limiting the number of right hand sides to <= nb means that all right hand 
//    sides are stored on process 0.  Hence, they can be treated as a serial 
//    matrix of vectors.
//
//
//  Note:
//    An efficient ScaLAPACK solution would differ from this in two ways:
//    1)  We would first redistribute the data in a ScaLAPACK form, but 
//    on all processors.
//    2)  We would then use the ScaLAPACK redistribution routine 
//    to redistribute the data to a 2D block cyclic data distribution
//    3)  We would move the vector(s) to 
//


//
int Amesos_Scalapack::RedistributeA( ) {

  if( debug_ == 1 ) cout << "Entering `RedistributeA()'" << endl;

  Time.ResetStartTime();
  
  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  EPETRA_CHK_ERR( RowMatrixA == 0 ) ; 

  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ; 
  int NumberOfProcesses = Comm().NumProc() ; 

  //
  //  Compute a uniform distribution as ScaLAPACK would want it
  //    MyFirstElement - The first element which this processor would have
  //    NumExpectedElemetns - The number of elements which this processor would have
  //

  int NumRows_ = RowMatrixA->NumGlobalRows() ; 
  if ( MaxProcesses_ > 0 ) {
    NumberOfProcesses = EPETRA_MIN( NumberOfProcesses, MaxProcesses_ ) ; 
  }
  else {
    int ProcessNumHeuristic = (1+NumRows_/200)*(1+NumRows_/200);
    NumberOfProcesses = EPETRA_MIN( NumberOfProcesses,  ProcessNumHeuristic );
  }

  //
  //   Compute nprow_ and npcol_ 
  //
  nprow_ = 1 ;
  npcol_ = NumberOfProcesses / nprow_ ;
  assert ( nprow_ * npcol_ == NumberOfProcesses ) ; 

  m_per_p_ = ( NumRows_ + NumberOfProcesses - 1 ) / NumberOfProcesses ;
  int MyFirstElement = EPETRA_MIN( iam_ * m_per_p_, NumRows_ ) ;
  int MyFirstNonElement = EPETRA_MIN( (iam_+1) * m_per_p_, NumRows_ ) ;
  int NumExpectedElements = MyFirstNonElement - MyFirstElement ; 

  assert( NumRows_ ==  RowMatrixA->NumGlobalRows() ) ; 
  if ( ScaLAPACK1DMap_ ) delete( ScaLAPACK1DMap_ ) ; 
  ScaLAPACK1DMap_ = new Epetra_Map( NumRows_, NumExpectedElements, 0, Comm() );
  if ( ScaLAPACK1DMatrix_ ) delete( ScaLAPACK1DMatrix_ ) ; 
  ScaLAPACK1DMatrix_ = new Epetra_CrsMatrix(Copy, *ScaLAPACK1DMap_, 0);
  Epetra_Export ExportToScaLAPACK1D_( OriginalMap, *ScaLAPACK1DMap_);

  ScaLAPACK1DMatrix_->Export( *RowMatrixA, ExportToScaLAPACK1D_, Add ); 
  
  ScaLAPACK1DMatrix_->FillComplete() ; 

  int info; 
  const int zero = 0 ; 
  if ( ictxt_ == -1313 ) {
    ictxt_ = 0 ; 
    SL_INIT_F77(&ictxt_, &nprow_, &npcol_) ; 
  }
  int nprow;
  int npcol;
  int myrow;
  int mycol;
  BLACS_GRIDINFO_F77(&ictxt_, &nprow, &npcol, &myrow, &mycol) ; 
  if ( iam_ < nprow_ * npcol_ ) { 
    assert( nprow == nprow_ ) ; 
    assert( npcol == npcol_ ) ; 
    assert( myrow == 0 ) ; 
    assert( mycol == iam_ ) ; 
    DESCINIT_F77(DescA_, 
		 &NumGlobalElements_, 
		 &NumGlobalElements_, 
		 &m_per_p_,            //  Irrelevant as nprow_ = 1, but may affect blocking
		 &m_per_p_,
		 &zero,
		 &zero,
		 &ictxt_,
		 &NumGlobalElements_,
		 &info) ;
    assert( info == 0 ) ; 
  } else {
    DescA_[0] = -13;
    assert( nprow == -1 ) ; 
  }

  MatTime_ += Time.ElapsedTime();
  
  return 0;
}


int Amesos_Scalapack::ConvertToScalapack(){
  
  //
  //  Convert matrix and vector to the form that Scalapack expects
  //  ScaLAPACK accepts the matrix to be in any 2D block-cyclic form
  //  We have chosen the simplest 2D block-cyclic form, a 1D blocked (not-cyclic)
  //  data distribution, for the matrix A.
  //  We use the same distribution for the multivectors X and B.  However, 
  //  except for very large numbers of right hand sides, this places all of X and B
  //  on process 0, making it effectively a serial matrix.  
  //  
  //  For now, we simply treat X and B as serial matrices (as viewed from epetra)
  //  though ScaLAPACK treats them as distributed matrices. 
  //
  //  int MyActualFirstElement = ScaLAPACK1DMatrix_->RowMatrixRowMap().MinMyGID() ;

  if( debug_ == 1 ) cout << "Entering `ConvertToScalapack()'" << endl;

  Time.ResetStartTime();
  
  int NumMyElements = ScaLAPACK1DMatrix_->NumMyRows() ; 
  Comm().Barrier();
  if ( iam_ < nprow_ * npcol_ ) { 
    
    //  int nnz_loc = ScaLAPACK1DMatrix_->NumMyNonzeros() ;
    assert( NumGlobalElements_ ==ScaLAPACK1DMatrix_->NumGlobalRows());
    assert( NumGlobalElements_ ==ScaLAPACK1DMatrix_->NumGlobalCols());
    DenseA_.resize( NumGlobalElements_ * NumMyElements ) ;
    for ( int i = 0 ; i < (int)DenseA_.size() ; i++ ) DenseA_[i] = 0 ; 
    
    int NzThisRow ;
    int MyRow;
    int num_my_cols = ScaLAPACK1DMatrix_->NumMyCols() ; 
    
    double *RowValues;
    int *ColIndices;
    int MaxNumEntries = ScaLAPACK1DMatrix_->MaxNumEntries();
    
    assert( DescA_[8] == NumGlobalElements_ ) ; //  Double check Lda
    vector<int>ColIndicesV(MaxNumEntries);
    vector<double>RowValuesV(MaxNumEntries);
    
    Epetra_CrsMatrix *SuperluCrs = dynamic_cast<Epetra_CrsMatrix *>(ScaLAPACK1DMatrix_);
    for ( MyRow = 0; MyRow < NumMyElements ; MyRow++ ) {
      if ( SuperluCrs != 0 ) {
	EPETRA_CHK_ERR( SuperluCrs->
			ExtractMyRowView( MyRow, NzThisRow, RowValues, 
					  ColIndices ) != 0 ) ;
      }
      else {
	//      ColIndicesV.resize(MaxNumEntries);
	//      RowValuesV.resize(MaxNumEntries);
	EPETRA_CHK_ERR( ScaLAPACK1DMatrix_->
			ExtractMyRowCopy( MyRow, MaxNumEntries, 
					  NzThisRow, &RowValuesV[0], 
					  &ColIndicesV[0] ) != 0 );
	RowValues =  &RowValuesV[0];
	ColIndices = &ColIndicesV[0];
      }
      for ( int j = 0; j < NzThisRow; j++ ) { 
	DenseA_[ ( ScaLAPACK1DMatrix_->RowMatrixColMap().GID( ColIndices[j] ) ) 
		 + MyRow * NumGlobalElements_ ] = RowValues[j] ; 
      }
    }
  }
  
  //
  //  Create a map to allow us to redistribute the vectors X and B 
  //
  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ; 
  assert( NumGlobalElements_ == OriginalMap.NumGlobalElements() ) ;
  int NumMyElements_ = 0 ;
  if (iam_==0) NumMyElements_ = NumGlobalElements_;
  
  if (SerialMap_) { delete SerialMap_ ; SerialMap_ = 0 ; } 
  SerialMap_ = new Epetra_Map( NumGlobalElements_, NumMyElements_, 0, Comm() );

  ConTime_ += Time.ElapsedTime();
  
  return 0;
}   


int Amesos_Scalapack::SetParameters( Teuchos::ParameterList &ParameterList ) {

  if( debug_ == 1 ) cout << "Entering `SetParameters()'" << endl;

  //
  //  We have to set these to their defaults here because user codes 
  //  are not guaranteed to have a "Scalapack" parameter list.
  //
  MaxProcesses_ = - 1; 

  //  Some compilers reject the following cast:
  //  if( &ParameterList == 0 ) return 0;

  // MS // now comment it out (only because the list if empty).
  // MS // When we will have parameters for KLU sublist
  // MS // uncomment it

  // define on how many processes matrix should be converted into ScaLAPCK
  // format. (this value must be less than available procs)
  if( ParameterList.isParameter("MaxProcs") )
    MaxProcesses_ = ParameterList.get("MaxProcs",-1);
  
  // ========================================= //
  // retrive ScaLAPACK's parameters from list. //
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

  // possible debug statements
  // 0 - no debug
  // 1 - debug
  if( ParameterList.isParameter("DebugLevel") )
    debug_ = ParameterList.get("DebugLevel",0);
  
  // MS // now comment it out, if we have parameters for KLU sublist
  // MS // uncomment it
  /*
  if (ParameterList.isSublist("Scalapack") ) {
    Teuchos::ParameterList ScalapackParams = ParameterList.sublist("Scalapack") ;
    MaxProcesses_ = ScalapackParams.get("MaxProcesses",MaxProcesses_);
  }  
  */
  
  return 0;
}

int Amesos_Scalapack::PerformNumericFactorization( ) {

  if( debug_ == 1 ) cout << "Entering `PerformNumericFactorization()'" << endl;
  
  Time.ResetStartTime();  

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ; 
  Ipiv_.resize(NumGlobalElements_) ;

  int Ierr[1] ; 
  Ierr[0] = 0 ; 
  const int one = 1 ; 
  if ( iam_ < nprow_ * npcol_ ) {
    PDGETRF_F77(&NumGlobalElements_,  
		&NumGlobalElements_, 
		&DenseA_[0],
		&one,
		&one, 
		DescA_,
		&Ipiv_[0],
		Ierr) ;
  }

  //  All processes should return the same error code
  if ( nprow_ * npcol_ < Comm().NumProc() ) 
    Comm().Broadcast( Ierr, 1, 0 ) ; 

  NumTime_ += Time.ElapsedTime();

  return Ierr[0];
}




bool Amesos_Scalapack::MatrixShapeOK() const { 
  bool OK ;

  if ( GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
       GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints() ) OK = false;
  return OK; 
}


int Amesos_Scalapack::SymbolicFactorization() {

  if( debug_ == 1 ) cout << "Entering `PerformSymbolicFactorization()'" << endl;

  NumSymbolicFact_++;
  
  return 0;
}

int Amesos_Scalapack::NumericFactorization() {

  if( debug_ == 1 ) cout << "Entering `NumericFactorization()'" << endl;

  NumNumericFact_++;

  iam_ = Comm().MyPID();
  
  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  const Epetra_Map &OriginalMap = RowMatrixA->RowMatrixRowMap() ; 
  NumGlobalElements_ = OriginalMap.NumGlobalElements();

  RedistributeA();
  ConvertToScalapack();

  return PerformNumericFactorization( );
}


int Amesos_Scalapack::Solve() { 

  if( debug_ == 1 ) cout << "Entering `Solve()'" << endl;

  NumSolve_++;

  Epetra_MultiVector   *vecX = Problem_->GetLHS() ; 
  Epetra_MultiVector   *vecB = Problem_->GetRHS() ; 

  //
  //  Compute the number of right hands sides 
  //  (and check that X and B have the same shape) 
  //
  int nrhs; 
  if ( vecX == 0 ) { 
    nrhs = 0 ;
    EPETRA_CHK_ERR( vecB != 0 ) ; 
  } else { 
    nrhs = vecX->NumVectors() ; 
    EPETRA_CHK_ERR( vecB->NumVectors() != nrhs ) ; 
  }

  Epetra_MultiVector *SerialB =0;
  Epetra_MultiVector *SerialX =0;
  //
  //  Extract Serial versions of X and B 
  //
  double *SerialXvalues ;

  Epetra_RowMatrix *RowMatrixA = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  Epetra_CrsMatrix *CastCrsMatrixA = dynamic_cast<Epetra_CrsMatrix*>(RowMatrixA) ; 
  Time.ResetStartTime(); // track time to broadcast vectors
  //
  //  Copy B to the serial version of B
  //
  const Epetra_Map &OriginalMap = CastCrsMatrixA->RowMap();
  Epetra_MultiVector *SerialXextract = new Epetra_MultiVector( *SerialMap_, nrhs ) ; 
  Epetra_MultiVector *SerialBextract = new Epetra_MultiVector( *SerialMap_, nrhs ) ; 
  
  Epetra_Import ImportToSerial( *SerialMap_, OriginalMap );
  SerialBextract->Import( *vecB, ImportToSerial, Insert ) ;
  SerialB = SerialBextract ; 
  SerialX = SerialXextract ; 

  VecTime_ += Time.ElapsedTime();

  //
  //  Call SCALAPACKs PDGETRS to perform the solve
  //

  int DescX[10];  
  
  SerialX->Scale(1.0, *SerialB) ;  

  int SerialXlda ; 

  Time.ResetStartTime(); // tract time to solve

  //
  //  Setup DescX 
  //
  assert( nrhs <= m_per_p_ ) ;   

  int Ierr[1] ; 
  Ierr[0] = 0 ; 
  const int zero = 0 ; 
  const int one = 1 ; 
  if ( iam_ < nprow_ * npcol_ ) {
    assert( SerialX->ExtractView( &SerialXvalues, &SerialXlda ) == 0 ) ; 

    assert( iam_ >0 || SerialXlda == NumGlobalElements_ ) ; 
    
    DESCINIT_F77(DescX, 
		 &NumGlobalElements_, 
		 &nrhs, 
		 &m_per_p_,            //  Irrelevant as nprow_ = 1
		 &m_per_p_,            //  nrhs would be a better choice but not legal
		 &zero,
		 &zero,
		 &ictxt_,
		 &NumGlobalElements_,
		 Ierr ) ;
    assert( Ierr[0] == 0 ) ; 
		
    //
    //  In ScaLAPACK we factor the transposed matrix, hence
    //  we must invert the sense of the transposition
    //
    char trans = 'T';
    if ( UseTranspose() ) trans = 'N' ;

    PDGETRS_F77(&trans,
		&NumGlobalElements_,  
		&nrhs, 
		&DenseA_[0],
		&one,
		&one, 
		DescA_,
		&Ipiv_[0],
		SerialXvalues,
		&one,
		&one, 
		DescX,
		Ierr ) ;
  }

  SolTime_ += Time.ElapsedTime();

  Time.ResetStartTime();  // track time to broadcast vectors
  //
  //  Copy X back to the original vector
  // 
  Epetra_Import ImportFromSerial( OriginalMap, *SerialMap_ );
  vecX->Import( *SerialX, ImportFromSerial, Insert ) ;
  delete SerialBextract ;
  delete SerialXextract ;

  VecTime_ += Time.ElapsedTime();

  //  All processes should return the same error code
  if ( nprow_ * npcol_ < Comm().NumProc() ) 
    Comm().Broadcast( Ierr, 1, 0 ) ; 

  // MS // compute vector norms
  if( ComputeVectorNorms_ == true || verbose_ == 2 ) {
    double NormLHS, NormRHS;
    for( int i=0 ; i<nrhs ; ++i ) {
      assert((*vecX)(i)->Norm2(&NormLHS)==0);
      assert((*vecB)(i)->Norm2(&NormRHS)==0);
      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Scalapack : vector " << i << ", ||x|| = " << NormLHS
	     << ", ||b|| = " << NormRHS << endl;
      }
    }
  }
  
  // MS // compute true residual
  if( ComputeTrueResidual_ == true || verbose_ == 2  ) {
    double Norm;
    Epetra_MultiVector Ax(vecB->Map(),nrhs);
    for( int i=0 ; i<nrhs ; ++i ) {
      (Problem_->GetMatrix()->Multiply(UseTranspose(), *((*vecX)(i)), Ax));
      (Ax.Update(1.0, *((*vecB)(i)), -1.0));
      (Ax.Norm2(&Norm));
      
      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Scalapack : vector " << i << ", ||Ax - b|| = " << Norm << endl;
      }
    }
  }
  
  return Ierr[0];

}


// ================================================ ====== ==== ==== == =

void Amesos_Scalapack::PrintStatus() 
{

  if( Comm().MyPID() != 0  ) return;

  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Scalapack : Matrix has " << Problem_->GetMatrix()->NumGlobalRows() << " rows"
       << " and " << Problem_->GetMatrix()->NumGlobalNonzeros() << " nonzeros" << endl;
  cout << "Amesos_Scalapack : Nonzero elements per row = "
       << 1.0*Problem_->GetMatrix()->NumGlobalNonzeros()/Problem_->GetMatrix()->NumGlobalRows() << endl;
  cout << "Amesos_Scalapack : Percentage of nonzero elements = "
       << 100.0*Problem_->GetMatrix()->NumGlobalNonzeros()/
    (pow(Problem_->GetMatrix()->NumGlobalRows(),2.0)) << endl;
  cout << "Amesos_Scalapack : Use transpose = " << UseTranspose_ << endl;
  cout << "----------------------------------------------------------------------------" << endl;

  return;
  
}

// ================================================ ====== ==== ==== == =

void Amesos_Scalapack::PrintTiming()
{
  if( Comm().MyPID() ) return;
  
  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Scalapack : Time to convert matrix to ScaLAPACK format = "
       << ConTime_ << " (s)" << endl;
  cout << "Amesos_Scalapack : Time to redistribute matrix = "
       << MatTime_ << " (s)" << endl;
  cout << "Amesos_Scalapack : Time to redistribute vectors = "
       << VecTime_ << " (s)" << endl;
  cout << "Amesos_Scalapack : Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << "Amesos_Scalapack : Time for sym fact = "
       << SymTime_ << " (s), avg = " << SymTime_/NumSymbolicFact_
       << " (s)" << endl;
  cout << "Amesos_Scalapack : Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << "Amesos_Scalapack : Time for num fact = "
       << NumTime_ << " (s), avg = " << NumTime_/NumNumericFact_
       << " (s)" << endl;
  cout << "Amesos_Scalapack : Number of solve phases = "
       << NumSolve_ << endl;
  cout << "Amesos_Scalapack : Time for solve = "
       << SolTime_ << " (s), avg = " << SolTime_/NumSolve_
       << " (s)" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
   
  return;
}
