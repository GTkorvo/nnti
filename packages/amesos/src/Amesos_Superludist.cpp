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

// TO DO: use Stat structure to print statistics???
// allow users to specify usermap ???

#include "Amesos_Superludist.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Util.h"
// #include "CrsMatrixTranspose.h"

// ====================================================================== 
int Superludist_NumProcRows( int NumProcs ) {
#ifdef TFLOP
  //  Else, parameter 6 of DTRSV CTNLU is incorrect 
  return 1;
#else
  int i;
  int NumProcRows ;
  for ( i = 1; i*i <= NumProcs; i++ ) 
    ;
  bool done = false ;
  for ( NumProcRows = i-1 ; done == false ; ) {
    int NumCols = NumProcs / NumProcRows ; 
    if ( NumProcRows * NumCols == NumProcs ) 
      done = true; 
    else 
      NumProcRows-- ; 
  }
  return NumProcRows;
#endif
}

// ====================================================================== 
int SetNPRowAndCol(const int MaxProcesses, int& nprow, int& npcol)
{
  nprow = Superludist_NumProcRows(MaxProcesses);
  if (nprow < 1 ) nprow = 1;
  npcol = MaxProcesses / nprow;
  
  if( nprow <=0 || npcol <= 0 || MaxProcesses <=0 ) {
    cerr << "Amesos_Superludist: wrong value for MaxProcs ("
	 << MaxProcesses << "), or nprow (" << nprow 
	 << ") or npcol (" << npcol << ")" << endl;
    AMESOS_CHK_ERR(-1);
  }
}

//=============================================================================
Amesos_Superludist::Amesos_Superludist(const Epetra_LinearProblem &prob) :
  Problem_(&prob),
  GridCreated_(0), 
  FactorizationDone_(0), 
  NumGlobalRows_(0), 
  RowMatrixA_(0), 
  nprow_(0),               // Overwritten by call to SetParameters below
  npcol_(0),               // Overwritten by call to SetParameters below
  PrintNonzeros_(false),
  ColPerm_("NOT SET"),     // Overwritten by call to SetParameters below
  RowPerm_("NOT SET"),     // Overwritten by call to SetParameters below
  perm_c_(0),              // Overwritten by call to SetParameters below
  perm_r_(0),              // Overwritten by call to SetParameters below
  IterRefine_("NOT SET"),  // Overwritten by call to SetParameters below
  ReplaceTinyPivot_(true), // Overwritten by call to SetParameters below - WAS false
  FactorizationOK_(false),    
  NumNumericFact_(0),
  NumSolve_(0)
{
  Redistribute_ = true;
  AddZeroToDiag_ = false;
  FactOption_ = SamePattern_SameRowPerm ;
  ReuseSymbolic_ = false ; 

  MaxProcesses_ = - 1; 
  nprow_ = 0;
  npcol_ = 0;
  Equil_ = true;
  ColPerm_ = "MMD_AT_PLUS_A";
  perm_c_ = 0;
  RowPerm_ = "LargeDiag";
  perm_r_ = 0;
  IterRefine_ = "DOUBLE";
  ReplaceTinyPivot_ = true;
  
  PrintNonzeros_ = false;

  ComputeTrueResidual_ = false;
  ComputeVectorNorms_ = false;
  PrintStatus_ = false ; 
  PrintTiming_ = false ; 
  
  Teuchos::ParameterList ParamList;
  SetParameters(ParamList); 
}

//=============================================================================
Amesos_Superludist::~Amesos_Superludist(void) 
{
  if (PrintTiming_) PrintTiming();
  if (PrintStatus_) PrintStatus();

  if ( FactorizationDone_ ) {
    SUPERLU_FREE( SuperluA_.Store );
    ScalePermstructFree(&ScalePermstruct_);
    Destroy_LU(NumGlobalRows_, &grid_, &LUstruct_);
    LUstructFree(&LUstruct_);
    if ( options_.SolveInitialized ) {
      dSolveFinalize(&options_, &SOLVEstruct_ ) ; 
    }
  }
  if ( GridCreated_ ) {
    superlu_gridexit(&grid_);
  }
}

// ====================================================================== 
int Amesos_Superludist::SetParameters( Teuchos::ParameterList &ParameterList ) 
{
  if (ParameterList.isParameter("Redistribute"))
    Redistribute_ = ParameterList.get("Redistribute",Redistribute_);  

  if (ParameterList.isParameter("AddZeroToDiag"))
    AddZeroToDiag_ = ParameterList.get("AddZeroToDiag",AddZeroToDiag_); 

  if (ParameterList.isParameter("AddToDiag"))
    AddToDiag_ = ParameterList.get("AddToDiag", AddToDiag_);

  // print some statistics (on process 0). Do not include timing
  if (ParameterList.isParameter("PrintStatus"))
    PrintStatus_ = ParameterList.get("PrintStatus", PrintStatus_);

  // print some statistics (on process 0). Do not include timing
  if (ParameterList.isParameter("PrintTiming"))
    PrintTiming_ = ParameterList.get("PrintTiming", PrintTiming_);

  if (ParameterList.isParameter("MaxProcs")) 
    MaxProcesses_ = ParameterList.get("MaxProcs",MaxProcesses_);

  if (ParameterList.isParameter("ComputeTrueResidual"))
    ComputeTrueResidual_ = ParameterList.get("ComputeTrueResidual",ComputeTrueResidual_);

  if (ParameterList.isParameter("ComputeVectorNorms") )
    ComputeVectorNorms_ = ParameterList.get("ComputeVectorNorms",ComputeVectorNorms_);

  // parameters for Superludist only

  if (ParameterList.isSublist("Superludist") ) 
  {
    Teuchos::ParameterList& SuperludistParams = 
      ParameterList.sublist("Superludist") ;

    if( SuperludistParams.isParameter("ReuseSymbolic") )
      ReuseSymbolic_ = SuperludistParams.get("ReuseSymbolic",ReuseSymbolic_);
    string FactOption = "NotSet";

    if( SuperludistParams.isParameter("Fact") )
      FactOption = SuperludistParams.get("Fact", FactOption);

    if( FactOption == "SamePattern_SameRowPerm" ) FactOption_ = SamePattern_SameRowPerm;
    else if( FactOption == "SamePattern" ) FactOption_ = SamePattern;
    else if ( FactOption != "NotSet" ) 
      AMESOS_CHK_ERR(-2); // input not valid

    if (SuperludistParams.isParameter("Equil"))
      Equil_ = SuperludistParams.get("Equil",Equil_);

    if (SuperludistParams.isParameter("ColPerm"))
      ColPerm_ = SuperludistParams.get("ColPerm",ColPerm_);

    if (ColPerm_ == "MY_PERMC")
      if( SuperludistParams.isParameter("perm_c"))
        perm_c_ = SuperludistParams.get("perm_c",perm_c_);

    if (SuperludistParams.isParameter("RowPerm"))
      RowPerm_ = SuperludistParams.get("RowPerm",RowPerm_);
    if( RowPerm_ == "MY_PERMR" ) {
      if (SuperludistParams.isParameter("perm_r"))
        perm_r_ = SuperludistParams.get("perm_r",perm_r_);
    }

    if (SuperludistParams.isParameter("IterRefine"))
      IterRefine_ = SuperludistParams.get("IterRefine",IterRefine_);

    if (SuperludistParams.isParameter("ReplaceTinyPivot"))
      ReplaceTinyPivot_ = SuperludistParams.get("ReplaceTinyPivot",ReplaceTinyPivot_);

    if (SuperludistParams.isParameter("PrintNonzeros"))
      PrintNonzeros_ = SuperludistParams.get("PrintNonzeros",PrintNonzeros_);
  }

  return(0);
}

// ====================================================================== 
// Tasks of this method:
// 1) To set the required number of processes;
// 2) To set nprow_ and npcol_
// 3) To create a linear distribution (map) with elements only in the
//    active processes
// 4) To redistribute the matrix from the original to the linear map.
// ====================================================================== 
int Amesos_Superludist::RedistributeA()
{
  ResetTime();
  
  if (NumGlobalRows_ != RowMatrixA_->NumGlobalRows())
    AMESOS_CHK_ERR(-1); // something has changed

  int iam = Comm().MyPID();
  int NumberOfProcesses = Comm().NumProc();
  if (MaxProcesses_ > NumberOfProcesses)
    MaxProcesses_ = NumberOfProcesses;

  SetMaxProcesses(MaxProcesses_, *RowMatrixA_);
  SetNPRowAndCol(MaxProcesses_, nprow_, npcol_);

  int m_per_p = NumGlobalRows_ / MaxProcesses_ ;
  int remainder = NumGlobalRows_ - ( m_per_p * MaxProcesses_ );
  int MyFirstElement = iam * m_per_p + EPETRA_MIN( iam, remainder );
  int MyFirstNonElement = (iam+1) * m_per_p + EPETRA_MIN( iam+1, remainder );
  int NumExpectedElements = MyFirstNonElement - MyFirstElement ; 

  if (iam >= MaxProcesses_)
    NumExpectedElements = 0; 

  if (NumGlobalRows_ !=  RowMatrixA_->NumGlobalRows())
    AMESOS_CHK_ERR(-1);

  const Epetra_Map &OriginalMap = RowMatrixA_->RowMatrixRowMap();

  UniformMap_       = rcp(new Epetra_Map(NumGlobalRows_, NumExpectedElements, 
                                         0, Comm()));
  CrsUniformMatrix_ = rcp(new Epetra_CrsMatrix(Copy, *UniformMap_, 0));
  UniformMatrix_    = rcp(&CrsUniformMatrix(), false);
  Importer_         = rcp(new Epetra_Import(*UniformMap_, OriginalMap));

  CrsUniformMatrix_->Import(*RowMatrixA_, Importer(), Add); 
  CrsUniformMatrix_->FillComplete(); 

  AddTime("matrix conversion");
  
  return(0);
}

// ====================================================================== 
int Amesos_Superludist::Factor()
{
  // FIXME????
  //  For now, if you change the shape of a matrix, you need to 
  //  create a new Amesos instance.
  //  
  //
  if (NumGlobalRows_ != 0 && NumGlobalRows_ != RowMatrixA_->NumGlobalRows())
    AMESOS_CHK_ERR(-5);

  NumGlobalRows_ = RowMatrixA_->NumGlobalRows() ; 

  if (Comm().NumProc() == 1)
    Redistribute_ = false;

  // Set the matrix and grid shapes. Time is tracked within
  // the RedistributeA() function

  if (Redistribute_)
    RedistributeA() ; 
  else 
  {
    if (Comm().NumProc() == 1)
    {
      nprow_ = 1;
      npcol_ = 1;
    }
    else 
    {
      if (!(RowMatrixA_->RowMatrixRowMap().LinearMap())) 
        AMESOS_CHK_ERR(-2);
      SetNPRowAndCol(Comm().NumProc(), nprow_, npcol_);
    }

    UniformMatrix_ = rcp(RowMatrixA_, false);
  }

  //  Extract Ai_, Ap_ and Aval_ from UniformMatrix_

  ResetTime();

  int MyActualFirstElement = UniformMatrix().RowMatrixRowMap().MinMyGID() ; 
  int NumMyElements = UniformMatrix().NumMyRows() ; 
  int nnz_loc = UniformMatrix().NumMyNonzeros() ;
  Ap_.resize( NumMyElements+1 );
  Ai_.resize( EPETRA_MAX( NumMyElements, nnz_loc) ) ; 
  Aval_.resize( EPETRA_MAX( NumMyElements, nnz_loc) ) ; 
  
  int NzThisRow ;
  int Ai_index = 0 ; 
  int MyRow;
  int num_my_cols = UniformMatrix().NumMyCols() ; 
  double *RowValues;
  int *ColIndices;
  int MaxNumEntries_ = UniformMatrix().MaxNumEntries();
  vector<double> RowValuesV_(MaxNumEntries_);
  vector<int>    ColIndicesV_(MaxNumEntries_);

  Global_Columns_ = UniformMatrix().RowMatrixColMap().MyGlobalElements();

  const Epetra_CrsMatrix *SuperluCrs = dynamic_cast<const Epetra_CrsMatrix *>(&UniformMatrix());

  int ierr;

  for (MyRow = 0; MyRow < NumMyElements ; MyRow++) 
  {
    if (SuperluCrs != 0) 
    {
      ierr = SuperluCrs->ExtractMyRowView(MyRow, NzThisRow, RowValues, 
                                      ColIndices);

    }
    else {
      ierr = UniformMatrix().ExtractMyRowCopy(MyRow, MaxNumEntries_, NzThisRow,
                                              &RowValuesV_[0], &ColIndicesV_[0]);
      RowValues =  &RowValuesV_[0];
      ColIndices = &ColIndicesV_[0];
    }

    AMESOS_CHK_ERR(ierr);

    // MS // Added on 15-Mar-05
    if (AddToDiag_ != 0.0 || AddZeroToDiag_) {
      for (int i = 0 ; i < NzThisRow ; ++i) {
        if (ColIndices[i] == MyRow) {
          RowValues[i] += AddToDiag_;
          break;
        }
      }
    }

    Ap_[MyRow] = Ai_index ; 
    for ( int j = 0; j < NzThisRow; j++ ) { 
      Ai_[Ai_index] = Global_Columns_[ColIndices[j]] ; 
      Aval_[Ai_index] = RowValues[j] ;
      Ai_index++;
    }
  }
  assert( NumMyElements == MyRow );
  Ap_[ NumMyElements ] = Ai_index ; 

  //
  //  Setup Superlu's grid 
  //
  const Epetra_MpiComm & comm1 = dynamic_cast<const Epetra_MpiComm &> (Comm());

  if ( ! GridCreated_ ) {
    // NOTE: nprow_ and npcol_ cannot be changed by the user
    GridCreated_ = true;
    superlu_gridinit(comm1.Comm(), nprow_, npcol_, &grid_);
  }

  if ( FactorizationDone_ ) {
    SUPERLU_FREE( SuperluA_.Store );
    ScalePermstructFree(&ScalePermstruct_);
    Destroy_LU(NumGlobalRows_, &grid_, &LUstruct_);
    LUstructFree(&LUstruct_);
    if ( options_.SolveInitialized ) {
      dSolveFinalize(&options_, &SOLVEstruct_ ) ; 
    }
  }

  AddTime("matrix conversion");
  ResetTime();

  //
  //  Only those processes in the grid participate from here on
  //
  if (Comm().MyPID() < nprow_ * npcol_) {
    //
    //  Set up Superlu's data structures
    //
    set_default_options_dist(&options_);

    dCreate_CompRowLoc_Matrix_dist( &SuperluA_, NumGlobalRows_, NumGlobalRows_, 
				    nnz_loc, NumMyElements, MyActualFirstElement,
				    &Aval_[0], &Ai_[0], &Ap_[0], 
				    SLU_NR_loc, SLU_D, SLU_GE );

    FactorizationDone_ = true;   // i.e. clean up Superlu data structures in the destructor

    ScalePermstructInit(NumGlobalRows_, NumGlobalRows_, &ScalePermstruct_);
    LUstructInit(NumGlobalRows_, NumGlobalRows_, &LUstruct_);

    // stick options from ParameterList to options_ structure
    // Here we follow the same order of the SuperLU_dist 2.0 manual (pag 55/56)
    
    assert( options_.Fact == DOFACT );  
    options_.Fact = DOFACT ;       

    if( Equil_ ) options_.Equil = (yes_no_t)YES;
    else         options_.Equil = NO;

    if( ColPerm_ == "NATURAL" ) options_.ColPerm = NATURAL;
    else if( ColPerm_ == "MMD_AT_PLUS_A" ) options_.ColPerm = MMD_AT_PLUS_A;
    else if( ColPerm_ == "MMD_ATA" ) options_.ColPerm = MMD_ATA;
    else if( ColPerm_ == "COLAMD" ) options_.ColPerm = COLAMD;
    else if( ColPerm_ == "MY_PERMC" ) {
      options_.ColPerm = MY_PERMC;
      ScalePermstruct_.perm_c = perm_c_;
    }

    if( RowPerm_ == "NATURAL" ) options_.RowPerm = (rowperm_t)NATURAL;
    if( RowPerm_ == "LargeDiag" ) options_.RowPerm = LargeDiag;
    else if( ColPerm_ == "MY_PERMR" ) {
      options_.RowPerm = MY_PERMR;
      ScalePermstruct_.perm_r = perm_r_;
    }

    if( ReplaceTinyPivot_ ) options_.ReplaceTinyPivot = (yes_no_t)YES;
    else                    options_.ReplaceTinyPivot = (yes_no_t)NO;

    if( IterRefine_ == "NO" ) options_.IterRefine = (IterRefine_t)NO;
    else if( IterRefine_ == "DOUBLE" ) options_.IterRefine = DOUBLE;
    else if( IterRefine_ == "EXTRA" ) options_.IterRefine = EXTRA;

    //  Without the following two lines, SuperLU_DIST cannot be made
    //  quiet.
    if (PrintNonzeros_) options_.PrintStat = (yes_no_t)YES;
    else                options_.PrintStat = (yes_no_t)NO;
    
    SuperLUStat_t stat;
    PStatInit(&stat);    /* Initialize the statistics variables. */

    //
    //  Factor A using Superludsit (via a call to pdgssvx)
    //
    int info ;
    double berr ;    //  Should be untouched
    double xValues;  //  Should be untouched
    int nrhs = 0 ;   //  Prevents forward and back solves
    int ldx = NumGlobalRows_;     //  Should be untouched

    pdgssvx(&options_, &SuperluA_, &ScalePermstruct_, &xValues, ldx, nrhs, &grid_,
	    &LUstruct_, &SOLVEstruct_, &berr, &stat, &info);

    if ( options_.SolveInitialized ) {
      dSolveFinalize(&options_, &SOLVEstruct_ ) ; 
    }
    AMESOS_CHK_ERR(info);

    PStatFree(&stat);
  }

  AddTime("numeric");

  return 0;
}

// ====================================================================== 
//   Refactor - Refactor the matrix 
//
//     Preconditions:
//       The non-zero pattern of the matrix must not have changed since the 
//         previous call to Factor().  Refactor ensures that each process owns 
//         the same number of columns that it did on the previous call to Factor()
//         and returns -4 if a discrepancy is found.  However, that check does not
//         guarantee that no change was made to the non-zero structure of the matrix.
//       No call to SetParameters should be made between the call to Factor()
//         and the call to Refactor().  If the user does not call SetParameters, 
//         as they need never do, they are safe on this.
//
//     Postconditions:
//       The matrix specified by Problem_->Operator() will have been redistributed,
//         converted to the form needed by Superludist and factored.
//       Ai_, Aval_ 
//       SuperluA_
//       SuperLU internal data structures reflecting the LU factorization
//         ScalePermstruct_
//         LUstructInit_
//       
//     Performance notes:
//       Refactor does not allocate or de-allocate memory.
//         
//     Return codes:
//       -4 if we detect a change to the non-zero structure of the matrix.
//
int Amesos_Superludist::ReFactor( ) 
{
  ResetTime();

  //
  //  Update Ai_ and Aval_ (while double checking Ap_)
  //
  if (Redistribute_)  
    if(CrsUniformMatrix().Import(*RowMatrixA_, Importer(), Insert)) 
      AMESOS_CHK_ERR(-4); 

  AddTime("matrix redistribution");
  ResetTime();

  const Epetra_CrsMatrix *SuperluCrs = dynamic_cast<const Epetra_CrsMatrix *>(&UniformMatrix());

  double *RowValues;
  int *ColIndices;
  int MaxNumEntries_ = UniformMatrix().MaxNumEntries();
  int NumMyElements  = UniformMatrix().NumMyRows() ; 
  vector<int> ColIndicesV_(MaxNumEntries_);
  vector<double> RowValuesV_(MaxNumEntries_);

  int NzThisRow ;
  int Ai_index = 0 ; 
  int MyRow;
  int ierr;

  for ( MyRow = 0; MyRow < NumMyElements ; MyRow++ ) {
    if ( SuperluCrs != 0 ) {
      ierr = SuperluCrs->ExtractMyRowView(MyRow, NzThisRow, RowValues, 
                                          ColIndices);
    }
    else {
      ierr = UniformMatrix().ExtractMyRowCopy( MyRow, MaxNumEntries_,
                                              NzThisRow, &RowValuesV_[0], 
                                              &ColIndicesV_[0]);
      RowValues =  &RowValuesV_[0];
      ColIndices = &ColIndicesV_[0];
    }

    AMESOS_CHK_ERR(ierr);

    if ( Ap_[MyRow] != Ai_index ) AMESOS_CHK_ERR(-4);
    for ( int j = 0; j < NzThisRow; j++ ) { 
      //  pdgssvx alters Ai_, so we have to set it again.
      Ai_[Ai_index] = Global_Columns_[ColIndices[j]];
      Aval_[Ai_index] = RowValues[j] ;  
      Ai_index++;
    }
  }
  if( Ap_[ NumMyElements ] != Ai_index ) AMESOS_CHK_ERR(-4); 

  AddTime("matrix conversion");
  ResetTime();

  if (Comm().MyPID() < nprow_ * npcol_) {

    set_default_options_dist(&options_);

    options_.Fact = FactOption_;
    SuperLUStat_t stat;
    PStatInit(&stat);    /* Initialize the statistics variables. */
    int info ;
    double berr ;    //  Should be untouched
    double xValues;  //  Should be untouched
    int nrhs = 0 ;   //  Prevents forward and back solves
    int ldx = NumGlobalRows_;     //  Should be untouched
    pdgssvx(&options_, &SuperluA_, &ScalePermstruct_, &xValues, ldx, nrhs, &grid_,
            &LUstruct_, &SOLVEstruct_, &berr, &stat, &info);
    PStatFree(&stat);
    AMESOS_CHK_ERR( info ) ;
  } 

  AddTime("numeric");

  return 0;
}

// ====================================================================== 
bool Amesos_Superludist::MatrixShapeOK() const 
{ 
  if (GetProblem()->GetOperator()->OperatorRangeMap().NumGlobalPoints() != 
      GetProblem()->GetOperator()->OperatorDomainMap().NumGlobalPoints()) 
    return(false);
  else
    return(true);
}

// ====================================================================== 
int Amesos_Superludist::SymbolicFactorization() 
{
  FactorizationOK_ = false ; 

  return(0);
}

// ====================================================================== 
int Amesos_Superludist::NumericFactorization() 
{
  RowMatrixA_ = dynamic_cast<Epetra_RowMatrix *>(Problem_->GetOperator());
  if (RowMatrixA_ == 0) 
    AMESOS_CHK_ERR(-1); // Linear problem does not contain Epetra_RowMatrix

  // reset factorization
  FactorizationOK_ = false; 

  if (!MatrixShapeOK())
    AMESOS_CHK_ERR(-1); // matrix not square

  InitTime(Comm());

  if (FactorizationOK_ && ReuseSymbolic_)
    ReFactor();
  else
    Factor();

  FactorizationOK_ = true;   
  NumNumericFact_++;
  
  return(0);
}

// ====================================================================== 
// We proceed as follows:
// - perform numeric factorization if not done;
// - extract solution and right-hand side;
// - redistribute them if necessary by creating additional
//   working vectors, or use vecX and vecB otherwise;     
// - call SuperLU_DIST's solve routine;                    
// - re-ship data if necessary.                             
// ====================================================================== 
int Amesos_Superludist::Solve() 
{
  if (!FactorizationOK_)
    AMESOS_CHK_ERR(NumericFactorization());

  Epetra_MultiVector* vecX = Problem_->GetLHS(); 
  Epetra_MultiVector* vecB = Problem_->GetRHS(); 

  if (vecX == 0 || vecB == 0)
    AMESOS_CHK_ERR(-1);

  int nrhs = vecX->NumVectors() ; 
  if (vecB->NumVectors() != nrhs)
    AMESOS_CHK_ERR(-1);

  double *values;
  int ldx;

  RefCountPtr<Epetra_MultiVector> vec_uni;

  if (Redistribute_) 
  {
    vec_uni = Teuchos::rcp(new Epetra_MultiVector(*UniformMap_, nrhs)); 
    ResetTime();
    vec_uni->Import(*vecB, Importer(), Insert);
    AddTime("vector redistribution");
  } 
  else 
  {
    vecX->Update(1.0, *vecB, 0.0);
    vec_uni = Teuchos::rcp(vecX, false); 
  }

  int NumMyElements = vec_uni->MyLength(); 
  AMESOS_CHK_ERR(vec_uni->ExtractView(&values, &ldx)); 

  ResetTime();
  
  /* Bail out if I do not belong in the grid. */
  if (Comm().MyPID() < nprow_ * npcol_) 
  {
    int info ;
    vector<double>berr(nrhs);
    SuperLUStat_t stat;
    PStatInit(&stat);    /* Initialize the statistics variables. */
    
    if (!GridCreated_ || !FactorizationDone_)
      AMESOS_CHK_ERR(-1); // internal error
    options_.Fact = FACTORED ;       

    bool BlockSolve = true ; 
    pdgssvx(&options_, &SuperluA_, &ScalePermstruct_, &values[0], ldx, 
            nrhs, &grid_, &LUstruct_, &SOLVEstruct_, &berr[0], 
            &stat, &info);
    AMESOS_CHK_ERR(info);
    
    PStatFree(&stat);
  }

  AddTime("solve");
  
  if (Redistribute_) 
  {
    ResetTime();
    vecX->Export(*vec_uni, Importer(), Insert);
    AddTime("vector redistribution");
  }

  if (ComputeTrueResidual_)
    ComputeTrueResidual(*RowMatrixA_, *vecX, *vecB, UseTranspose(), 
                        "Amesos_Superludist");

  if (ComputeVectorNorms_)
    ComputeVectorNorms(*vecX, *vecB, "Amesos_Superludist");

  NumSolve_++;
  
  return(0);
}

// ====================================================================== 
void Amesos_Superludist::PrintStatus() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;
  
  string p = "Amesos_Superludist : ";
  int NNZ = RowMatrixA_->NumGlobalNonzeros();

  PrintLine();

  cout << p << "Matrix has " << NumGlobalRows_ << " rows"
       << " and " << NNZ << " nonzeros" << endl;
  cout << p << "Nonzero elements per row = "
       << 1.0 * NNZ / NumGlobalRows_ << endl;
  cout << p << "Percentage of nonzero elements = "
       << 100.0 * NNZ /pow(NumGlobalRows_, 2.0) << endl;
  cout << p << "Use transpose = " << UseTranspose() << endl;
  cout << p << "Redistribute = " << Redistribute_ << endl;
  cout << p << "# available processes = " << Comm().NumProc() << endl;
  cout << p << "# processes used in computation = " << nprow_ * npcol_
       << " ( = " << nprow_ << "x" << npcol_ << ")" << endl;
  cout << p << "Equil = " << Equil_ << endl;
  cout << p << "ColPerm = " << ColPerm_ << endl;
  cout << p << "RowPerm = " << RowPerm_ << endl;
  cout << p << "IterRefine = " << IterRefine_ << endl;
  cout << p << "ReplaceTinyPivot = " << ReplaceTinyPivot_ << endl;
  cout << p << "AddZeroToDiag = " << AddZeroToDiag_ << endl;
  cout << p << "Redistribute = " << Redistribute_ << endl;
  
  PrintLine();

  return;
}

// ====================================================================== 
void Amesos_Superludist::PrintTiming() const
{
  if (Problem_->GetOperator() == 0 || Comm().MyPID() != 0)
    return;

  double ConTime = GetTime("matrix conversion");
  double MatTime = GetTime("matrix redistribution");
  double VecTime = GetTime("vector redistribution");
  double NumTime = GetTime("numeric");
  double SolTime = GetTime("solve");

  if (NumNumericFact_)
    NumTime /= NumNumericFact_;

  if (NumSolve_)
    SolTime /= NumSolve_;

  string p = "Amesos_Superludist : ";
  PrintLine();

  cout << p << "Time to convert matrix to Klu format = "
       << ConTime << " (s)" << endl;
  cout << p << "Time to redistribute matrix = "
       << MatTime << " (s)" << endl;
  cout << p << "Time to redistribute vectors = "
       << VecTime << " (s)" << endl;
  cout << p << "Number of symbolic factorizations = "
       << NumSymbolicFact_ << endl;
  cout << p << "Time for sym fact = 0.0 (s), avg = 0.0 (s)" << endl;
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
