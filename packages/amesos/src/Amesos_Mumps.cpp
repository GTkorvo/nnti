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
 
#include "Amesos_Mumps.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Epetra_Time.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Vector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Amesos_EpetraBaseSolver.h"
#include "EpetraExt_Redistor.h"
#include "Epetra_Util.h"

#define ICNTL(I) icntl[(I)-1]
#define CNTL(I)  cntl[(I)-1]
#define INFOG(I) infog[(I)-1]
#define INFO(I) info[(I)-1]
#define RINFOG(I) rinfog[(I)-1]

const int DEF_VALUE_INT = -123456789;
const double DEF_VALUE_DOUBLE = -123456.789;
  
//=============================================================================

Amesos_Mumps::Amesos_Mumps(const Epetra_LinearProblem &prob ) :
  Amesos_EpetraBaseSolver(prob),
  SymbolicFactorizationOK_(false), 
  NumericFactorizationOK_(false),
  KeepMatrixDistributed_(true),
  MaxProcs_(-1),
  MaxProcsInputMatrix_(-1),
  UseMpiCommSelf_(false),
  Map_(0),
  NumMUMPSNonzeros_(0),
  NumMyMUMPSNonzeros_(0),
  Col(0),
  Row(0),
  Val(0),
  NumSchurComplementRows_(-1),
  SchurComplementRows_(0),
  CrsSchurComplement_(0),
  DenseSchurComplement_(0),
  IsComputeSchurComplementOK_(false),
  RowSca_(0),
  ColSca_(0),
  PermIn_(0),
  Maxis_(DEF_VALUE_INT),
  Maxs_(DEF_VALUE_INT),
  OldMatrix_(0),
  Redistor_(0),
  TargetVector_(0),
  verbose_(1),
  debug_(0),
  AddToDiag_(0.0),
  AddDiagElement_(false),
  PrintTiming_(false),
  PrintStatus_(false),
  Threshold_(0.0),
  MUMPSComm_(0),
  UseTranspose_(false),
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
  // -777 is for me. It means : never called MUMPS, so
  // SymbolicFactorization will not call Destroy();
  MDS.job = -777;
  
  // set to -1 icntl_ and cntl_. The use can override default values by using
  // SetICNTL(pos,value) and SetCNTL(pos,value).
  for( int i=0 ; i<40 ; ++i ) icntl_[i] = DEF_VALUE_INT;
  for( int i=0 ; i<5 ; ++i ) cntl_[i] = DEF_VALUE_DOUBLE;
}

//=============================================================================

void Amesos_Mumps::Destroy()
{
  if( debug_ == 1 ) cout << "Entering `Destroy()'..." << endl;
  
  // destroy instance of the package
  MDS.job = -2;

  if( Comm().MyPID() < MaxProcs_ ) dmumps_c(&MDS);
  
  if( Redistor_ ) delete Redistor_;
  if( TargetVector_ ) delete TargetVector_;
  
  if( Row ) { delete Row; Row = 0; }
  if( Col ) { delete Col; Col = 0; }
  if( Val ) { delete Val; Val = 0; }

  if( IsComputeSchurComplementOK_ && Comm().MyPID() == 0 ) {
    delete [] MDS.schur;
  }

  if( MUMPSComm_ ) MPI_Comm_free( &MUMPSComm_ );
  MUMPSComm_ = 0;

  if( (verbose_ && PrintTiming_) || verbose_ == 2 ) PrintTiming();
  if( (verbose_ && PrintStatus_) || verbose_ == 2 ) PrintStatus();

  return;
}

//=============================================================================

Amesos_Mumps::~Amesos_Mumps(void)
{

  Destroy();
 
}

//=============================================================================

int Amesos_Mumps::ConvertToTriplet()
{

  if( debug_ == 1 ) cout << "Entering `ConvertToTriplet()'" << endl;
  
  Time.ResetStartTime();

  // MS // convert to MUMPS format, keeping in distributed form.
  // MS // This doesn't require the matrix to be shipped to proc 0,
  // MS // then decomposed and reshipped by MUMPS

  // MS // be sure that those vectors are empty. Allocate them

  if( Row ) delete Row;
  if( Col ) delete Col;
  if( Val ) delete Val;

  Row = new Epetra_IntSerialDenseVector( NumMyNonzeros()+1);
  Col = new Epetra_IntSerialDenseVector( NumMyNonzeros()+1 );
  Val = new Epetra_SerialDenseVector( NumMyNonzeros()+1 );

  NumMyMUMPSNonzeros_ = 0;
  
  // MS // retrive already allocated pointers for rows, cols, and vals
  // MS // Those arrays are allocated and freed by the Amesos_EpetraBaseSolver
  // MS // class. Note that, if a View method is implemented, those pointers
  // MS // will change to reflect the action of the View. Otherwise,
  // MS // the allocated vectors will be used to copy data  

  int NumIndices;
  
  int * RowIndices = GetRowIndices();
  int * ColIndices = GetColIndices();
  double * MatrixValues = GetValues();
  
  // MS // for CrsA and VbrA, handle the case of IndexBase different from 0
  // MS // diff is the difference between what we have and what we want 
  // MS // (1, FORTRAN style)
  int diff = 1-IndexBase();
  
  // MS // do all in terms of BlockRows. This way we can support VBR matrices

  for( int LocalBlockRow=0; LocalBlockRow<NumMyBlockRows() ; ++LocalBlockRow ) {

    bool FoundDiagonal = false;
    int Diag;
    
    GetRow(LocalBlockRow,NumIndices,RowIndices,ColIndices,MatrixValues);

    for( int i=0 ; i<NumIndices ; ++i ) {
      if( abs(MatrixValues[i]) >= Threshold_ || true ) {
	(*Row)[NumMyMUMPSNonzeros_] = RowIndices[i]+diff;
	(*Col)[NumMyMUMPSNonzeros_] = ColIndices[i]+diff;
	(*Val)[NumMyMUMPSNonzeros_] = MatrixValues[i];
	if( RowIndices[i] == ColIndices[i] ) {
	  Diag = RowIndices[i];
	  FoundDiagonal = true;
	  (*Val)[NumMyMUMPSNonzeros_] += AddToDiag_;
	}
	NumMyMUMPSNonzeros_++;
      }
      // add diagonal element if not found
      if( AddDiagElement_ && FoundDiagonal == false ) {
	(*Row)[NumMyMUMPSNonzeros_] = Diag;
	(*Col)[NumMyMUMPSNonzeros_] = Diag;
	(*Val)[NumMyMUMPSNonzeros_] = AddToDiag_;
	NumMyMUMPSNonzeros_++;
      }
    }
  }

  ConTime_ == Time.ElapsedTime(); 
  
  assert(NumMyMUMPSNonzeros_<=NumMyNonzeros());

  // MS // bring matrix to proc zero if required. Note that I first
  // MS // convert to COO format (locally), then I redistribute COO vectors.

  if( MaxProcsInputMatrix_ != Comm().NumProc() ) {

    RedistributeMatrix(MaxProcsInputMatrix_);
    
  } else {

    NumMUMPSNonzeros_ = NumMyMUMPSNonzeros_;

  }

  return 0;

}

//=============================================================================

void Amesos_Mumps::RedistributeMatrix(const int NumProcs)
{

  if( debug_ == 1 ) cout << "Entering `RedistributeMatrix()' ..." << endl;

  Time.ResetStartTime();
  
  Epetra_IntSerialDenseVector * OldRow = Row;
  Epetra_IntSerialDenseVector * OldCol = Col;
  Epetra_SerialDenseVector * OldVal = Val;

  Comm().SumAll(&NumMyMUMPSNonzeros_,&NumMUMPSNonzeros_,1);

  int local = NumMUMPSNonzeros_ / NumProcs;
  if( Comm().MyPID() == 0 ) local += NumMUMPSNonzeros_%NumProcs;

  if( Comm().MyPID() < NumProcs ) NumMUMPSNonzeros_ = local;
  else                            NumMUMPSNonzeros_ = 0;

  
  Row = new Epetra_IntSerialDenseVector(NumMUMPSNonzeros_+1);
  Col = new Epetra_IntSerialDenseVector(NumMUMPSNonzeros_+1);
  Val = new Epetra_SerialDenseVector(NumMUMPSNonzeros_+1);

  Epetra_Map OldNnz(-1,NumMyMUMPSNonzeros_,0,Comm());
    
  Epetra_Map NewNnz(-1,NumMUMPSNonzeros_,0,Comm());
  Epetra_Import OldToNew(NewNnz,OldNnz);

  Epetra_IntVector GOldRow(View,OldNnz,OldRow->Values());
  Epetra_IntVector GOldCol(View,OldNnz,OldCol->Values());
  Epetra_Vector GOldVal(View,OldNnz,OldVal->Values());
  
  Epetra_IntVector GRow(View,NewNnz,Row->Values());
  Epetra_IntVector GCol(View,NewNnz,Col->Values());
  Epetra_Vector GVal(View,NewNnz,Val->Values());

  GRow.Import(GOldRow,OldToNew,Insert);
  GCol.Import(GOldCol,OldToNew,Insert);
  GVal.Import(GOldVal,OldToNew,Insert);

  if( OldRow ) { delete OldRow; OldRow = 0; }
  if( OldCol ) { delete OldCol; OldCol = 0; }
  if( OldVal ) { delete OldVal; OldVal = 0; }

  MatTime_ += Time.ElapsedTime();
  
}

//=============================================================================

void Amesos_Mumps::RedistributeMatrixValues(const int NumProcs)
{

  if( debug_ == 1 ) cout << "Entering `RedistributeMatrixValues()' ..." << endl;

  Time.ResetStartTime();
  
  // I suppose that NumMyMUMPSNonzeros_ has been computed
  // before calling this method.

  Comm().SumAll(&NumMyMUMPSNonzeros_,&NumMUMPSNonzeros_,1);    

  Epetra_SerialDenseVector * OldVal = Val;
  
  int local = NumMUMPSNonzeros_ / NumProcs;
  if( Comm().MyPID() == 0 ) local += NumMUMPSNonzeros_%NumProcs;
  
  if( Comm().MyPID() < NumProcs ) NumMUMPSNonzeros_ = local;
  else                            NumMUMPSNonzeros_ = 0;

  Val = new Epetra_SerialDenseVector(NumMUMPSNonzeros_);
  
  Epetra_Map OldNnz(-1,NumMyMUMPSNonzeros_,0,Comm());
  
  Epetra_Map NewNnz(-1,NumMUMPSNonzeros_,0,Comm());
  Epetra_Import OldToNew(NewNnz,OldNnz);
  
  Epetra_Vector GOldVal(View,OldNnz,OldVal->Values());
  
  Epetra_Vector GVal(View,NewNnz,Val->Values());
  
  GVal.Import(GOldVal,OldToNew,Insert);
  
  if( OldVal ) { delete OldVal; OldVal = 0; }

  MatTime_ += Time.ElapsedTime();
  
}

//=============================================================================

int Amesos_Mumps::ConvertToTripletValues()
{

  if( debug_ == 1 ) cout << "Entering `ConvertToTripletValues()'" << endl;

  Time.ResetStartTime();  

  // MS // convert to MUMPS format, keeping in distributed form.
  // MS // This doesn't require the matrix to be shipped to proc 0,
  // MS // then decomposed and reshipped by MUMPS

  // MS // be sure that those vectors are empty. Allocate them
  
  if( Val ) delete Val;

  Val = new Epetra_SerialDenseVector( NumMyNonzeros() ) ;  

  int NumMUMPSNonzerosValues = 0;
  
  // MS // retrive already allocated pointers for rows, cols, and vals
  // MS // Those arrays are allocated and freed by the Amesos_EpetraBaseSolver
  // MS // class. Note that, if a View method is implemented, those pointers
  // MS // will change to reflect the action of the View. Otherwise,
  // MS // the allocated vectors will be used to copy data  

  int NumIndices;
  
  int * RowIndices = GetRowIndices();
  int * ColIndices = GetColIndices();
  double * MatrixValues = GetValues();
  
  // MS // do all in terms of BlockRows. This way we can support VBR matrices
  
  for( int LocalBlockRow=0; LocalBlockRow<NumMyBlockRows() ; ++LocalBlockRow ) {

    bool FoundDiagonal = false;

    GetRow(LocalBlockRow,NumIndices,RowIndices,ColIndices,MatrixValues);
    /*    
    for( int i=0 ; i<NumIndices ; ++i ) {
      if( abs(MatrixValues[i]) >= Threshold_ ) {
	
	(*Val)[NumMUMPSNonzerosValues] = MatrixValues[i];
	NumMUMPSNonzerosValues++;
      }
    }
    */

    for( int i=0 ; i<NumIndices ; ++i ) {
      if( abs(MatrixValues[i]) >= Threshold_ || true ) {
	(*Val)[NumMUMPSNonzerosValues] = MatrixValues[i];
	if( RowIndices[i] == ColIndices[i] ) {
	  FoundDiagonal = true;
	  (*Val)[NumMUMPSNonzerosValues] += AddToDiag_;
	}
	NumMUMPSNonzerosValues++;
      }

      // add diagonal element if not found
      if( AddDiagElement_ && FoundDiagonal == false ) {
	(*Val)[NumMUMPSNonzerosValues++] = AddToDiag_;
      }
    }
    
  }

  if( NumMUMPSNonzerosValues != NumMyMUMPSNonzeros_ ) {
     cerr << "Amesos_Mumps ERROR : it appears that matrix has changed\n"
          << "Amesos_Mumps ERROR : its structure since SymbolicFactorization() has\n"
	  << "Amesos_Mumps ERROR : called. # nonzeros (sym) = " << NumMyMUMPSNonzeros_ << endl
	  << "Amesos_Mumps ERROR : called. # nonzeros (num) = " << NumMUMPSNonzerosValues << endl;
     return -2;
  }  

  ConTime_ += Time.ElapsedTime();
  
  // MS // bring matrix to proc zero if required
  
  if( MaxProcsInputMatrix_ != Comm().NumProc() ) {

    RedistributeMatrixValues(MaxProcsInputMatrix_);

  } else {
    
    NumMUMPSNonzeros_ = NumMyMUMPSNonzeros_;

  }

  return 0;

}

//=============================================================================

int Amesos_Mumps::SetICNTL(int pos, int value)
{
  // NOTE: suppose first position is 1 (as in FORTRAN)
  
  if( pos>0 && pos<41 ) {
    icntl_[pos-1] = value;
    return 0;
  } else {
    return -1;
  }
}

//=============================================================================

int Amesos_Mumps::SetCNTL(int pos, double value)
{
  // NOTE: suppose first position is 1 (as in FORTRAN)
  if( pos>0 && pos<6 ) {
    cntl_[pos-1] = value;
    return 0;
  } else {
    return -1;
  }
}

//=============================================================================

void Amesos_Mumps::SetICNTLandCNTL()
{
  
  MDS.ICNTL(1)  = -1;  // Turn off error messages
  MDS.ICNTL(2)  = -1;  // Turn off diagnostic printing
  MDS.ICNTL(3)  = -1;  // Turn off global information messages
  MDS.ICNTL(4)  = -1;
  MDS.ICNTL(5)  = 0;   // Matrix is given in elemental (i.e. triplet) from
  MDS.ICNTL(6)  = 7;   // Choose column permutation automatically
  MDS.ICNTL(7)  = 7;   // Choose ordering method automatically
  MDS.ICNTL(8)  = 7;   // Choose scaling automatically
  MDS.ICNTL(9)  = 1;   // Compute A x = b 
  MDS.ICNTL(10) = 0;   // Maximum steps of iterative refinement
  MDS.ICNTL(11) = 0;   // Do not collect statistics
  MDS.ICNTL(12) = 0;   // Use Node level parallelism
  MDS.ICNTL(13) = 0;   // Use ScaLAPACK for root node 
  MDS.ICNTL(14) = 20;  // Increase memory allocation 20% at a time 
  MDS.ICNTL(15) = 0;   // Minimize memory use (not flop count)
  MDS.ICNTL(16) = 0;   // Do not perform null space detection
  MDS.ICNTL(17) = 0;   // Unused (null space dimension)
  
  if( KeepMatrixDistributed_ ) MDS.ICNTL(18)= 3;
  else                         MDS.ICNTL(18)= 0;
  
  if ( UseTranspose() )  MDS.ICNTL(9) = 0 ; 
  else                   MDS.ICNTL(9) = 1 ;
  
  if( IsComputeSchurComplementOK_ ) MDS.ICNTL(19) = 1;
  else                              MDS.ICNTL(19) = 0;

  // something to do if the Schur complement is required.
  if( IsComputeSchurComplementOK_ && Comm().MyPID() == 0 ) {
    MDS.size_schur = NumSchurComplementRows_;
    MDS.listvar_schur = SchurComplementRows_;
    MDS.schur = new double[NumSchurComplementRows_*NumSchurComplementRows_];
  }
  
  // retrive user's specified options
  for( int i=0 ; i<40 ; ++i ) {
    if( icntl_[i] != DEF_VALUE_INT )
      MDS.ICNTL(i+1) = icntl_[i];
  }
  for( int i=0 ; i<5 ; ++i ) {
    if( icntl_[i] != DEF_VALUE_DOUBLE )
      MDS.CNTL(i+1) = cntl_[i];
  }
  
  // take care that required options are not overwritten
  assert(MDS.ICNTL(5)== 0); 
  
  return;
  
}

//=============================================================================

int Amesos_Mumps::SetParameters( Teuchos::ParameterList & ParameterList)
{

  if( debug_ == 1 ) cout << "Entering `SetParameters()' ..." << endl;
  
  // ========================================= //
  // retrive MUMPS' parameters from list.      //
  // default values defined in the constructor //
  // ========================================= //
  
  // retrive general parameters

  // solve problem with transpose
  if( ParameterList.isParameter("UseTranspose") )
    SetUseTranspose(ParameterList.get("UseTranspose",false));

  // ignore all elements below given tolerance
  if( ParameterList.isParameter("Threshold") )
    Threshold_ = ParameterList.get("Threshold", 0.0);

  // add zero to diagonal if diagonal element is not present
  if( ParameterList.isParameter("AddZeroToDiag") )
    AddDiagElement_ = ParameterList.get("AddZeroToDiag", false);

  // add this value to diagonal
  if( ParameterList.isParameter("AddToDiag") )
    AddToDiag_ = ParameterList.get("AddToDiag", 0.0);

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

  // on how many processes distribute the matrix
  if( ParameterList.isParameter("MaxProcsMatrix") )
    MaxProcsInputMatrix_ = ParameterList.get("MaxProcsMatrix",-1);

  // define on how many processes matrix should be converted into MUMPS
  // format. (this value must be less than available procs)
  if( ParameterList.isParameter("MaxProcs") )
    MaxProcs_ = ParameterList.get("MaxProcs",-1);

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
  
  // retrive MUMPS' specific parameters
  
  if (ParameterList.isSublist("mumps") ) {
    Teuchos::ParameterList MumpsParams = ParameterList.sublist("mumps") ;
    // integer array of parameters
    if( MumpsParams.isParameter("ICNTL") ) {
      int * ICNTL = 0;
      ICNTL = MumpsParams.get("ICNTL", ICNTL);
      if( ICNTL ) SetICNTL(ICNTL);
    }
    // double array of parameters
    if( MumpsParams.isParameter("CNTL") ) {
      double * CNTL = 0;
      CNTL = MumpsParams.get("CNTL", CNTL);
      if( CNTL ) SetCNTL(CNTL);
    }
    // ordering
     if( MumpsParams.isParameter("PermIn") ) {
      int * PermIn = 0;
      PermIn = MumpsParams.get("PermIn", PermIn);
      if( PermIn ) SetOrdering(PermIn);
    }
     // Maxis
     if( MumpsParams.isParameter("Maxis") ) {
       int Maxis = 0;
       Maxis = MumpsParams.get("Maxis", Maxis);
       SetMaxis(Maxis);
     }
     // Maxs
     if( MumpsParams.isParameter("Maxs") ) {
       int Maxs = 0;
       Maxs = MumpsParams.get("Maxs", Maxs);
       SetMaxs(Maxs);
     }
     // Col scaling
     if( MumpsParams.isParameter("ColScaling") ) {
       double * ColSca = 0;
       ColSca = MumpsParams.get("ColScaling", ColSca);
       if( ColSca ) SetColScaling(ColSca);
     }
     // Row scaling
     if( MumpsParams.isParameter("RowScaling") ) {
       double * RowSca = 0;
       RowSca = MumpsParams.get("RowScaling", RowSca);
       if( RowSca ) SetRowScaling(RowSca);
     }
     // that's all folks
  }  

  return 0;
}

//=============================================================================

void Amesos_Mumps::CheckParameters() 
{

  // FIXME ???
#ifndef HAVE_AMESOS_MPI_C2F
  MaxProcs_ = -3;
  MaxProcsInputMatrix_ = -3;
#endif

  // check parameters and fix values of MaxProcs_ and MaxProcsInputMatrix_

  int NumGlobalNonzeros=0, NumRows=0;
  
  NumGlobalNonzeros = GetMatrix()->NumGlobalNonzeros(); 
  NumRows == GetMatrix()->NumGlobalRows(); 

  // optimal value for MaxProcs == -1
  
  int OptNumProcs1 = 1+EPETRA_MAX( NumRows/10000, NumGlobalNonzeros/1000000 );
  OptNumProcs1 = EPETRA_MIN(Comm().NumProc(),OptNumProcs1 );

  // optimal value for MaxProcs == -2

  int OptNumProcs2 = (int)sqrt(1.0*Comm().NumProc());
  if( OptNumProcs2 < 1 ) OptNumProcs2 = 1;

  // fix the value of MaxProcs

  switch( MaxProcs_ ) {
  case -1:
    MaxProcs_ = OptNumProcs1;
    break;
  case -2:
    MaxProcs_ = OptNumProcs2;
    break;
  case -3:
    MaxProcs_ = Comm().NumProc();
    break;
  }

  // fix the value of MaxProcsInputMatrix
  
  switch( MaxProcsInputMatrix_ ) {
  case -1:
    MaxProcsInputMatrix_ = OptNumProcs1;
    break;
  case -2:
    MaxProcsInputMatrix_ = OptNumProcs2;
    break;
  case -3:
    MaxProcsInputMatrix_ = MaxProcs_;
    break;  
  }
  
  // check available processes; -1 means use all available processes
  if( MaxProcs_ < -3 || MaxProcs_ > Comm().NumProc() ) MaxProcs_ = Comm().NumProc();

  // check available processes; -1 means use all available processes
  if( MaxProcsInputMatrix_ < -3 ||  MaxProcsInputMatrix_ > MaxProcs_ )
    MaxProcsInputMatrix_ = MaxProcs_;

  // cannot distribute input matrix to this number,
  // use all the processes instead
  if( MaxProcsInputMatrix_ > MaxProcs_ ) {
    MaxProcsInputMatrix_ = MaxProcs_;
  }
  
  if( Comm().NumProc() == 1 || MaxProcsInputMatrix_ == 1 ) KeepMatrixDistributed_ = false;

  return;
  
}

//=============================================================================

int Amesos_Mumps::PerformSymbolicFactorization()
{

  if( debug_ == 1 ) cout << "Entering `PerformSymbolicFactorization()'" << endl;

  // erase data if present
  if( SymbolicFactorizationOK_ && MDS.job != -777 ) {
   Destroy();
  }

  if( IsLocal() || UseMpiCommSelf_ ) {
    // use --with-amesos-mpi-c2f to create sub-communicators
    // C to FORTRAN communicator is not standard
#if defined(EPETRA_MPI) && defined(HAVE_AMESOS_MPI_C2F)
    MPI_Comm MPIC = MPI_COMM_SELF ;
    MDS.comm_fortran = (F_INT) MPI_Comm_c2f( MPIC ) ;   // Compiled on cygwin but not on Atlantis
#else
    MDS.comm_fortran = -987654 ;  // Needed by MUMPS 4.3 
#endif
  } else {
    
#if defined(EPETRA_MPI) && defined(HAVE_AMESOS_MPI_C2F)
    if( MaxProcs_ != Comm().NumProc() ) {

      if( debug_ == 1 ) cout << "Creating MPI Communicator with "
			     << MaxProcs_ << " procs" << endl;
      
      if( MUMPSComm_ ) MPI_Comm_free( &MUMPSComm_ );
      
      int * ProcsInGroup = new int[MaxProcs_];
      for( int i=0 ; i<MaxProcs_ ; ++i ) ProcsInGroup[i] = i;
      MPI_Group OrigGroup, MumpsGroup;
      MPI_Comm_group(MPI_COMM_WORLD, &OrigGroup);

      MPI_Group_incl(OrigGroup, MaxProcs_, ProcsInGroup, &MumpsGroup);
      
      MPI_Comm_create(MPI_COMM_WORLD, MumpsGroup, &MUMPSComm_);

      MDS.comm_fortran = (F_INT) MPI_Comm_c2f( MUMPSComm_);

      delete [] ProcsInGroup;

    } else {
      // use MPI_COMM_WORLD
      // FIXME: should use the comm in Epetra_MpiComm()
      MDS.comm_fortran = -987654;  
    }
#else
    MDS.comm_fortran = -987654;
#endif
  }
  
  //  MDS.comm_fortran = (F_INT) MPIR_FromPointer( MPIC ) ;  // Other recommendation from the MUMPS manual - did not compile on Atlantis either
  
  MDS.job = -1  ;     //  Initialization
  MDS.par = 1 ;       //  Host IS involved in computations
  MDS.sym = MatrixProperty();

  if( Comm().MyPID() < MaxProcs_ ) dmumps_c( &MDS ) ;   //  Initialize MUMPS 
  CheckError();
  
  // check the error flag. MUMPS is successful only if
  //  infog(1) is zereo
  MDS.n = NumGlobalRows() ;

  // fix pointers for nonzero pattern of A. Numerical values
  // will be entered in PerformNumericalFactorization()
  if( KeepMatrixDistributed_ ) {

    MDS.nz_loc = NumMUMPSNonzeros_;
    
    if( Comm().MyPID() < MaxProcs_ ) {
      MDS.irn_loc = Row->Values(); 
      MDS.jcn_loc = Col->Values();

    }
  } else {
    if( Comm().MyPID() == 0 ) {
      MDS.nz = NumMUMPSNonzeros_;
      MDS.irn = Row->Values(); 
      MDS.jcn = Col->Values(); 
    }
  }

  // scaling if provided by the user
  if( RowSca_ != 0 ) {
    MDS.rowsca = RowSca_;
    MDS.colsca = ColSca_;
  }

  // given ordering if provided by the user
  if( PermIn_ != 0 ) {
    MDS.perm_in = PermIn_;
  }
  
  //  if( Maxis_ != DEF_VALUE_INT ) MDS.maxis = Maxis_;
  //  if( Maxs_ != DEF_VALUE_INT ) MDS.maxs = Maxs_;
  
  MDS.job = 1  ;     // Request symbolic factorization

  SetICNTLandCNTL(); // initialize icntl and cntl. NOTE: I initialize those vectors
                     // here, and I don't change them anymore. This is not that much
                     // a limitation, though

  // Perform symbolic factorization

  Time.ResetStartTime();  
  if( Comm().MyPID() < MaxProcs_ ) dmumps_c( &MDS ) ;
  SymTime_ += Time.ElapsedTime();

  CheckError();

  SymbolicFactorizationOK_ = true ;

  return 0;

}

//=============================================================================

int Amesos_Mumps::PerformNumericFactorization( )
{

  if( debug_ == 1 ) cout << "Entering `PerformNumericFactorization()'" << endl;

  // set vector for matrix entries. User may have changed it
  // since PerformSymbolicFactorization() has been called

  if( KeepMatrixDistributed_ ) {
    if( Comm().MyPID() < MaxProcs_ ) MDS.a_loc = Val->Values();
  } else {
    if( Comm().MyPID() == 0 ) {
      MDS.a = Val->Values();
    }
  }

  MDS.job = 2  ;     // Request numeric factorization
  // Perform numeric factorization
  Time.ResetStartTime();
  if( Comm().MyPID() < MaxProcs_ ) dmumps_c( &MDS ) ;
  NumTime_ += Time.ElapsedTime();
  
  CheckError();

  NumericFactorizationOK_ = true;

  return 0;
}

//=============================================================================

int Amesos_Mumps::SymbolicFactorization()
{

  NumSymbolicFact_++;  
  
  if( debug_ == 1 ) cout << "Entering `SymbolicFactorization()'" << endl;
  
  CheckParameters();
  
  // first, gather matrix property using base class
  // Amesos_EpetraBaseSolver

  Epetra_RowMatrix * Matrix = GetMatrix();
  assert( Matrix != NULL );
  if( Matrix != OldMatrix_ ) {
    // initialize redistor, target map is on proc 0.
    // Redistor will be used to bring rhs to proc 0, and solution to
    // all procs
    EPETRA_CHK_ERR( SetInterface(Matrix) );
    OldMatrix_ = Matrix;
    SymbolicFactorizationOK_ = false;
    NumericFactorizationOK_ = false;
  }

  // now, create a map, in which all nodes are stored on processor 0.
  // Also creates import objects. This is done using the base class
  // Amesos_EpetraRedistributor. NOTE: the following classes return
  // if they have already been called.  

  EPETRA_CHK_ERR(ConvertToTriplet());
  EPETRA_CHK_ERR(PerformSymbolicFactorization());

  return 0;
}

//=============================================================================

int Amesos_Mumps::NumericFactorization()
{

  NumNumericFact_++;  
  
  if( debug_ == 1 ) cout << "Entering `NumericFactorization()'" << endl;

  // MS // As in SymbolicFactorization. All those functions
  // MS // returns if they have already been called.
  
  Epetra_RowMatrix * Matrix = GetMatrix();
  assert( Matrix != NULL );  
  if( Matrix != OldMatrix_ ) {
    EPETRA_CHK_ERR( SetInterface(Matrix) );
    OldMatrix_ = Matrix;
    SymbolicFactorizationOK_ = false;
    NumericFactorizationOK_ = false;
  }

  // MS // reship the matrix. User can update numerical values
  // MS // Should the values remain the same, it is better to use
  // MS // Solve() and NOT SymbolicFactorization() + NumericFactorization()
  // MS // + Solve()
  
  if( SymbolicFactorizationOK_ == false ) {
    EPETRA_CHK_ERR( ConvertToTriplet() );
    EPETRA_CHK_ERR( PerformSymbolicFactorization() );
  } else {
    EPETRA_CHK_ERR( ConvertToTripletValues() );
  }  
  
  EPETRA_CHK_ERR( PerformNumericFactorization() );
  
  return 0;
}
 
//=============================================================================

int Amesos_Mumps::Solve()
{ 

  NumSolve_++;  
  
  if( debug_ == 1 ) cout << "Entering `Solve()'" << endl;

  Epetra_RowMatrix * Matrix = GetMatrix();
  assert( Matrix != NULL );
  if( Matrix != OldMatrix_ ) {
    EPETRA_CHK_ERR( SetInterface(Matrix) );
    OldMatrix_ = Matrix;
    SymbolicFactorizationOK_ = false;
    NumericFactorizationOK_ = false;
  }
  
  if( SymbolicFactorizationOK_ == false && NumericFactorizationOK_ == false ) {

    // MS // it is the first time the user call something. The
    // MS // matrix is neither sym fact not num fact
    // MS // We convert it to COO format, then sym+fact
    EPETRA_CHK_ERR( ConvertToTriplet() );
    EPETRA_CHK_ERR( PerformSymbolicFactorization() );
    EPETRA_CHK_ERR( PerformNumericFactorization() );

  }  else if( NumericFactorizationOK_  == false )  {

    // MS // User has already called SymFact, but maybe he changed
    // MS // the matrix entries. Ww reship the matrix
    EPETRA_CHK_ERR( ConvertToTripletValues() );
    EPETRA_CHK_ERR( PerformNumericFactorization() );

  } else {

    // MS // just for comment it out: This case, the user has already
    // MS // numerically factorized the matrix. We can solve without
    // MS // shipping the matrix
    
  }

  Epetra_MultiVector * vecX = GetLHS() ; 
  Epetra_MultiVector * vecB = GetRHS() ;


  EPETRA_CHK_ERR( UpdateLHS() );

  // create redistor the first time enters here. I suppose that LHS and RHS
  // (all of them) have the same map.
  if( IsLocal() == false && UseMpiCommSelf_ == false && Redistor_ == 0 ) {
    Redistor_ = new EpetraExt_Redistor(vecX->Map(),1);
    assert( Redistor_ != 0 );
  }
  
  //  Compute the number of right hands sides (and check that X and B have the same shape) 

  if( vecB == NULL || vecX == NULL ) {
    if( Comm().MyPID() == 0 ) cout << "Amesos_Mumps ERROR : LHS or RHS not set" << endl;
    return -1;
  }

  int nrhs = vecX->NumVectors() ;
  if( nrhs != vecB->NumVectors() ) {
    if( Comm().MyPID() == 0 ) cout << "Amesos_Mumps ERROR : multivectors LHS or RHS have different sizes" << endl;
    return -2;
  }

  // will be used only if IsLocal == false
  if( IsLocal() == false && UseMpiCommSelf_ == false && TargetVector_ == 0 ) {
    TargetVector_ = new Epetra_MultiVector(*(Redistor_->TargetMap()),nrhs);
    assert( TargetVector_ != 0 );
  }
  
  // MS // Extract Target versions of X and B. 
  // MS // NOTE: both Target and distributed version of MUMPS require
  // MS // rhs to be entirely stored on processor 0. This is because we
  // MS // still need TargetMap_ also for distributed input

  if( IsLocal() || UseMpiCommSelf_ ) {

    if( Comm().MyPID() == 0 ) {
      for ( int j =0 ; j < nrhs; j++ ) {
	for( int i=0 ; i<NumMyRows() ; ++i ) (*vecX)[j][i] = (*vecB)[j][i];
	MDS.rhs = (*vecX)[j];
	MDS.job = 3  ;     // Request solve
	Time.ResetStartTime();
	if( Comm().MyPID() < MaxProcs_ ) dmumps_c( &MDS ) ;  // Perform solve
	SolTime_ += Time.ElapsedTime();
	
	CheckError();
      }
    }
  } else {

    // bring rhs to process 0 and take timing for this phase
    Time.ResetStartTime();
    Redistor_->TargetImport( *vecB, *TargetVector_ ) ;
    VecTime_ += Time.ElapsedTime();
    
    for ( int j =0 ; j < nrhs; j++ ) { 
      if ( Comm().MyPID() == 0 ) {
	MDS.rhs = (*TargetVector_)[j];
      }
      // solve the linear system and take time
      MDS.job = 3  ;     
      Time.ResetStartTime();
      if( Comm().MyPID() < MaxProcs_ )dmumps_c( &MDS ) ;  // Perform solve
      SolTime_ += Time.ElapsedTime();
      
      CheckError();

    }

    // ship solution back and take timing
    Time.ResetStartTime();
    Redistor_->SourceImport( *vecX, *TargetVector_ ) ;
    VecTime_ += Time.ElapsedTime();

  }

  EPETRA_CHK_ERR( UpdateLHS() );

  // compute vector norms
  if( ComputeVectorNorms_ == true || verbose_ == 2 ) {
    double NormLHS, NormRHS;
    for( int i=0 ; i<nrhs ; ++i ) {
      assert((*vecX)(i)->Norm2(&NormLHS)==0);
      assert((*vecB)(i)->Norm2(&NormRHS)==0);
      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Mumps : vector " << i << ", ||x|| = " << NormLHS
	     << ", ||b|| = " << NormRHS << endl;
      }
    }
  }
  
  // compute true residual
  if( ComputeTrueResidual_ == true || verbose_ == 2  ) {
    double Norm;
    Epetra_MultiVector Ax(vecB->Map(),nrhs);
    for( int i=0 ; i<nrhs ; ++i ) {
      (GetMatrix()->Multiply(UseTranspose(), *((*vecX)(i)), Ax));
      (Ax.Update(1.0, *((*vecB)(i)), -1.0));
      (Ax.Norm2(&Norm));
      
      if( verbose_ && Comm().MyPID() == 0 ) {
	cout << "Amesos_Mumps : vector " << i << ", ||Ax - b|| = " << Norm << endl;
      }
    }
  }

  return(0) ; 
}

// ================================================ ====== ==== ==== == =

Epetra_CrsMatrix * Amesos_Mumps::GetCrsSchurComplement() 
{

  if( IsComputeSchurComplementOK_ ) {
    
    if( Comm().MyPID() != 0 ) NumSchurComplementRows_ = 0;
    
    Epetra_Map SCMap(-1,NumSchurComplementRows_, 0, Comm());

    CrsSchurComplement_ = new Epetra_CrsMatrix(Copy,SCMap,NumSchurComplementRows_);

    if( Comm().MyPID() == 0 )
      for( int i=0 ; i<NumSchurComplementRows_ ; ++i ) {
	for( int j=0 ; j<NumSchurComplementRows_ ; ++j ) {
	  int pos = i+ j *NumSchurComplementRows_;
	  CrsSchurComplement_->InsertGlobalValues(i,1,&(MDS.schur[pos]),&j);
	}
      }
    
    CrsSchurComplement_->FillComplete();

    return CrsSchurComplement_;
  }

  return 0;
  
}

// ================================================ ====== ==== ==== == =

Epetra_SerialDenseMatrix * Amesos_Mumps::GetDenseSchurComplement() 
{
  
  if( IsComputeSchurComplementOK_ ) {
    
    if( Comm().MyPID() != 0 ) return 0;
    
    DenseSchurComplement_ = new Epetra_SerialDenseMatrix(NumSchurComplementRows_,
							NumSchurComplementRows_);
    
    for( int i=0 ; i<NumSchurComplementRows_ ; ++i ) {
      for( int j=0 ; j<NumSchurComplementRows_ ; ++j ) {
	int pos = i+ j *NumSchurComplementRows_;
	(*DenseSchurComplement_)(i,j) = MDS.schur[pos];
      }
    }
    
    return DenseSchurComplement_;
    
  }
  
  return 0;
  
}

// ================================================ ====== ==== ==== == =

int Amesos_Mumps::ComputeSchurComplement(bool flag, int NumSchurComplementRows,
					 int * SchurComplementRows)
{
  NumSchurComplementRows_ = NumSchurComplementRows;
  
  SchurComplementRows_ = SchurComplementRows;

  // modify because MUMPS is fortran-driven
  if( Comm().MyPID() == 0 )
    for( int i=0 ; i<NumSchurComplementRows ; ++i ) SchurComplementRows_[i]++;
  
  IsComputeSchurComplementOK_ = flag;

  return 0;
}

// ================================================ ====== ==== ==== == =

void Amesos_Mumps::PrintStatus() 
{

  if( Comm().MyPID() != 0  ) return;

  cout << "----------------------------------------------------------------------------" << endl;
  cout << "Amesos_Mumps : Matrix has " << NumGlobalRows() << " rows"
       << " and " << NumGlobalNonzeros() << " nonzeros" << endl;
  cout << "Amesos_Mumps : Nonzero elements per row = "
       << 1.0*NumGlobalNonzeros()/NumGlobalRows() << endl;
  cout << "Amesos_Mumps : Percentage of nonzero elements = "
       << 100.0*NumGlobalNonzeros()/(pow(NumGlobalRows(),2.0)) << endl;
  cout << "Amesos_Mumps : Use transpose = " << UseTranspose_ << endl;
  cout << "Amesos_Mumps : Available process(es) = " << Comm().NumProc() << endl;
  cout << "Amesos_Mumps : Using " << MaxProcs_ << " process(es)" << endl;
  cout << "Amesos_Mumps : Input matrix distributed over " << MaxProcsInputMatrix_ << " process(es)" << endl;
  
  cout << "Amesos_Mumps : Estimated FLOPS for elimination = "
       << MDS.RINFOG(1) << endl;
  cout << "Amesos_Mumps : Total FLOPS for assembly = "
       << MDS.RINFOG(2) << endl;
  cout << "Amesos_Mumps : Total FLOPS for elimination = "
       << MDS.RINFOG(3) << endl;
  
  cout << "Amesos_Mumps : Total real space to store the LU factors = "
       << MDS.INFOG(9) << endl;
  cout << "Amesos_Mumps : Total integer space to store the LU factors = "
       << MDS.INFOG(10) << endl;
  cout << "Amesos_Mumps : Total number of iterative steps refinement = "
       << MDS.INFOG(15) << endl;
  cout << "Amesos_Mumps : Estimated size of MUMPS internal data\n"
       << "Amesos_Mumps : for running factorization = "
       << MDS.INFOG(16) << " Mbytes" << endl;
  cout << "Amesos_Mumps : for running factorization = "
       << MDS.INFOG(17) << " Mbytes" << endl;
  cout << "Amesos_Mumps : Allocated during factorization = "
       << MDS.INFOG(19) << " Mbytes" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
 
  return;
  
}

// ================================================ ====== ==== ==== == =

int Amesos_Mumps::SetICNTL(int * icntl)
{
  for( int i=0 ; i<40 ; ++i ) icntl_[i] = icntl[i];
  return 0;

}

// ================================================ ====== ==== ==== == =

int Amesos_Mumps::SetCNTL(double * cntl)  
{
  for( int i=0 ; i<5 ; ++i ) cntl_[i] = cntl[i];
  return 0;
}
 
// ================================================ ====== ==== ==== == =

void Amesos_Mumps::CheckError() 
{
  
  bool Wrong = MDS.INFOG(1) != 0 && Comm().MyPID() < MaxProcs_;
  
  // an error occurred in MUMPS. Print out information and quit.

  if( Comm().MyPID() == 0 && Wrong ) {
    cerr << "Amesos_Mumps : ERROR" << endl;
    cerr << "Amesos_Mumps : INFOG(1) = " << MDS.INFOG(1) << endl;
    cerr << "Amesos_Mumps : INFOG(2) = " << MDS.INFOG(2) << endl;
  }
  
  if( MDS.INFO(1) != 0 && Wrong ) {
    cerr << "Amesos_Mumps : On process " << Comm().MyPID()
	 << ", INFO(1) = " << MDS.INFO(1) << endl;
    cerr << "Amesos_Mumps : On process " << Comm().MyPID()
	 << ", INFO(2) = " << MDS.INFO(2) << endl;
  }

  if( Wrong ) exit( EXIT_FAILURE );
  
}

// ================================================ ====== ==== ==== == =

void Amesos_Mumps::PrintTiming()
{
  if( Comm().MyPID() ) return;
  
  cout << "-------------- -------------------------------------------------------------" << endl;
  cout << "Amesos_Mumps : Time to convert matrix to MUMPS format = "
       << ConTime_ << " (s)" << endl;
  if( MaxProcsInputMatrix_ != Comm().NumProc() ) 
    cout << "Amesos_Mumps : Time to redistribute matrix = "
	 << MatTime_ << " (s)" << endl;
  if( ! (IsLocal() || UseMpiCommSelf_) )
    cout << "Amesos_Mumps : Time to redistribute vectors = "
       << VecTime_ << " (s)" << endl;
  cout << "Amesos_Mumps : Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << "Amesos_Mumps : Time for sym fact = "
       << SymTime_ << " (s), avg = " << SymTime_/NumSymbolicFact_
       << " (s)" << endl;
  cout << "Amesos_Mumps : Number of numeric factorizations = "
       << NumNumericFact_ << endl;
  cout << "Amesos_Mumps : Time for num fact = "
       << NumTime_ << " (s), avg = " << NumTime_/NumNumericFact_
       << " (s)" << endl;
  cout << "Amesos_Mumps : Number of solve phases = "
       << NumSolve_ << endl;
  cout << "Amesos_Mumps : Time for solve = "
       << SolTime_ << " (s), avg = " << SolTime_/NumSolve_
       << " (s)" << endl;
  cout << "----------------------------------------------------------------------------" << endl;
   
  return;
}
