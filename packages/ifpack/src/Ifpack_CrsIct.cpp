/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#include "Ifpack_CrsIct.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Util.h"
#include "icrout_cholesky_mex.h"

//==============================================================================
Ifpack_CrsIct::Ifpack_CrsIct(const Epetra_CrsMatrix & A, double Droptol, int Lfil) 
  : A_(A),
    Comm_(A.Comm()),
    Allocated_(false),
    ValuesInitialized_(false),
    Factored_(false),
    Condest_(-1.0),
    Athresh_(0.0),
    Rthresh_(1.0),
    Droptol_(Droptol),
    Lfil_(Lfil),
    OverlapX_(0),
    OverlapY_(0),
    LevelOverlap_(0),
    OverlapMode_(Zero),
    Aict_(0),
    Lict_(0),
    Ldiag_(0)
{
  int ierr = Allocate();
}

//==============================================================================
Ifpack_CrsIct::Ifpack_CrsIct(const Ifpack_CrsIct & FactoredMatrix) 
  : A_(FactoredMatrix.A_),
    Comm_(FactoredMatrix.Comm_),
    Allocated_(FactoredMatrix.Allocated_),
    ValuesInitialized_(FactoredMatrix.ValuesInitialized_),
    Factored_(FactoredMatrix.Factored_),
    Condest_(FactoredMatrix.Condest_),
    Athresh_(FactoredMatrix.Athresh_),
    Rthresh_(FactoredMatrix.Rthresh_),
    Droptol_(FactoredMatrix.Droptol_),
    Lfil_(FactoredMatrix.Lfil_),
    OverlapX_(0),
    OverlapY_(0),
    LevelOverlap_(FactoredMatrix.LevelOverlap_),
    OverlapMode_(FactoredMatrix.OverlapMode_),
    Aict_(0),
    Lict_(0),
    Ldiag_(0)
  
{
  U_ = new Epetra_CrsMatrix(FactoredMatrix.U());
  D_ = new Epetra_Vector(FactoredMatrix.D());
  
}

//==============================================================================
int Ifpack_CrsIct::Allocate() {

  // Allocate Epetra_CrsMatrix using ILUK graphs
  if (LevelOverlap_==0) {
    U_ = new Epetra_CrsMatrix(Copy, A_.RowMatrixRowMap(), A_.RowMatrixRowMap(), 0);
    D_ = new Epetra_Vector(A_.RowMatrixRowMap());
  }
  else {
    EPETRA_CHK_ERR(-1); // LevelOverlap > 0 not implemented yet
    //    U_ = new Epetra_CrsMatrix(Copy, OverlapRowMap());
    //    D_ = new Epetra_Vector(OverlapRowMap());
  }
    
  
  
    SetAllocated(true);
    return(0);
}
//==============================================================================
Ifpack_CrsIct::~Ifpack_CrsIct(){


  delete U_;
  delete D_; // Diagonal is stored separately.  We store the inverse.

  if (OverlapX_!=0) delete OverlapX_;
  if (OverlapY_!=0) delete OverlapY_;

  if (Lict_!=0) {
    Matrix * Lict = (Matrix *) Lict_;
    free(Lict->ptr);
    free(Lict->col);
    free(Lict->val);
    delete Lict;
  }
  if (Aict_!=0) {
    Matrix * Aict = (Matrix *) Aict_;
    delete Aict;
  }
  if (Ldiag_!=0) free(Ldiag_);

  ValuesInitialized_ = false;
  Factored_ = false;
  Allocated_ = false;
}

//==========================================================================
int Ifpack_CrsIct::InitValues(const Epetra_CrsMatrix & A) {

  int ierr = 0;
  int i, j;
  int * InI=0, * LI=0, * UI = 0;
  double * InV=0, * LV=0, * UV = 0;
  int NumIn, NumL, NumU;
  bool DiagFound;
  int NumNonzeroDiags = 0;

  Epetra_CrsMatrix * OverlapA = (Epetra_CrsMatrix *) &A_;

  if (LevelOverlap_>0) {
    EPETRA_CHK_ERR(-1); // Not implemented yet
    //OverlapA = new Epetra_CrsMatrix(Copy, *Graph_.OverlapGraph());
    //EPETRA_CHK_ERR(OverlapA->Import(A, *Graph_.OverlapImporter(), Insert));
    //EPETRA_CHK_ERR(OverlapA->TransformToLocal());
  }
  // Get Maximun Row length
  int MaxNumEntries = OverlapA->MaxNumEntries();

  InI = new int[MaxNumEntries]; // Allocate temp space
  UI = new int[MaxNumEntries];
  InV = new double[MaxNumEntries];
  UV = new double[MaxNumEntries];

  double *DV;
  ierr = D_->ExtractView(&DV); // Get view of diagonal
    

  // First we copy the user's matrix into diagonal vector and U, regardless of fill level

  int NumRows = OverlapA->NumMyRows();

  for (i=0; i< NumRows; i++) {

    OverlapA->ExtractMyRowCopy(i, MaxNumEntries, NumIn, InV, InI); // Get Values and Indices
    
    // Split into L and U (we don't assume that indices are ordered).
    
    NumL = 0; 
    NumU = 0; 
    DiagFound = false;
    
    for (j=0; j< NumIn; j++) {
      int k = InI[j];

      if (k==i) {
	DiagFound = true;
	DV[i] += Rthresh_ * InV[j] + EPETRA_SGN(InV[j]) * Athresh_; // Store perturbed diagonal in Epetra_Vector D_
      }

      else if (k < 0) return(-1); // Out of range
      else if (i<k && k<NumRows) {
	UI[NumU] = k;
	UV[NumU] = InV[j];
	NumU++;
      }
    }
    
    // Check in things for this row of L and U

    if (DiagFound) NumNonzeroDiags++;
    if (NumU) U_->InsertMyValues(i, NumU, UV, UI);
    
  }

  delete [] UI;
  delete [] UV;
  delete [] InI;
  delete [] InV;

  if (LevelOverlap_>0 && U().DistributedGlobal()) delete OverlapA;


  U_->TransformToLocal(&A_.OperatorDomainMap(), &A_.OperatorRangeMap());
  SetValuesInitialized(true);
  SetFactored(false);

  int ierr1 = 0;
  if (NumNonzeroDiags<U_->NumMyRows()) ierr1 = 1;
  A_.Comm().MaxAll(&ierr1, &ierr, 1);
  EPETRA_CHK_ERR(ierr);
  return(0);
}

//==========================================================================
int Ifpack_CrsIct::Factor() {

  // if (!Allocated()) return(-1); // This test is not needed at this time.  All constructors allocate.
  if (!ValuesInitialized_) EPETRA_CHK_ERR(-2); // Must have values initialized.
  if (Factored()) EPETRA_CHK_ERR(-3); // Can't have already computed factors.

  SetValuesInitialized(false);

  int i;

  int m, n, nz, Nrhs, ldrhs, ldlhs;
  int * ptr=0, * ind;
  double * val, * rhs, * lhs;

  int ierr = Epetra_Util_ExtractHbData(U_, 0, 0, m, n, nz, ptr, ind,
			    val, Nrhs, rhs, ldrhs, lhs, ldlhs);
  if (ierr<0) EPETRA_CHK_ERR(ierr);

  Matrix * Aict;
  if (Aict_==0) {
    Aict = new Matrix;
    Aict_ = (void *) Aict;
  }
  else Aict = (Matrix *) Aict_;
  Matrix * Lict;
  if (Lict_==0) {
    Lict = new Matrix;
    Lict_ = (void *) Lict;
  }
  else Lict = (Matrix *) Lict_;
  Aict->val = val;
  Aict->col = ind;
  Aict->ptr = ptr;
  double *DV;
  EPETRA_CHK_ERR(D_->ExtractView(&DV)); // Get view of diagonal
    
  crout_ict(m, Aict, DV, Droptol_, Lfil_, Lict, &Ldiag_);

  // Get rid of unnecessary data
  delete [] ptr;
  delete U_;
  delete D_;

  // Create Epetra View of L from crout_ict

  if (LevelOverlap_==0) {
    U_ = new Epetra_CrsMatrix(View, A_.RowMatrixRowMap(), A_.RowMatrixRowMap(),0);
    D_ = new Epetra_Vector(View, A_.RowMatrixRowMap(), Ldiag_);
  }
  else {
    EPETRA_CHK_ERR(-1); // LevelOverlap > 0 not implemented yet
    //    U_ = new Epetra_CrsMatrix(Copy, OverlapRowMap());
    //    D_ = new Epetra_Vector(OverlapRowMap());
  }

  ptr = Lict->ptr;
  ind = Lict->col;
  val = Lict->val;
    
  for (i=0; i< m; i++) {
    int NumEntries = ptr[i+1]-ptr[i];
    int * Indices = ind+ptr[i];
    double * Values = val+ptr[i];
    U_->InsertMyValues(i, NumEntries, Values, Indices);
  }

  U_->TransformToLocal(&(A_.OperatorDomainMap()), &(A_.OperatorRangeMap()));
  
  D_->Reciprocal(*D_); // Put reciprocal of diagonal in this vector
  // Add up flops
 
  double current_flops = 2 * nz; // Just an estimate
  double total_flops = 0;
    
  A_.Comm().SumAll(&current_flops, &total_flops, 1); // Get total madds across all PEs

  // Now count the rest
  total_flops += (double) U_->NumGlobalNonzeros(); // Accounts for multiplier above
  total_flops += (double) D_->GlobalLength(); // Accounts for reciprocal of diagonal

  UpdateFlops(total_flops); // Update flop count

  SetFactored(true);

  return(0);

}

//=============================================================================
int Ifpack_CrsIct::Solve(bool Trans, const Epetra_MultiVector& X, 
				Epetra_MultiVector& Y) const {
//
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//

  if (X.NumVectors()!=Y.NumVectors()) EPETRA_CHK_ERR(-1); // Return error: X and Y not the same size

  bool Upper = true;
  bool UnitDiagonal = true;

  Epetra_MultiVector * X1 = (Epetra_MultiVector *) &X;
  Epetra_MultiVector * Y1 = (Epetra_MultiVector *) &Y;


  U_->Solve(Upper, true, UnitDiagonal, *X1, *Y1);
  Y1->Multiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
  U_->Solve(Upper, false, UnitDiagonal, *Y1, *Y1); // Solve Uy = y
  return(0);
}
//=============================================================================
int Ifpack_CrsIct::Multiply(bool Trans, const Epetra_MultiVector& X, 
				Epetra_MultiVector& Y) const {
//
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
//

  if (X.NumVectors()!=Y.NumVectors()) EPETRA_CHK_ERR(-1); // Return error: X and Y not the same size

  bool Upper = true;
  bool Lower = false;
  bool UnitDiagonal = true;

  Epetra_MultiVector * X1 = (Epetra_MultiVector *) &X;
  Epetra_MultiVector * Y1 = (Epetra_MultiVector *) &Y;

  U_->Multiply(false, *X1, *Y1);
  Y1->Update(1.0, *X1, 1.0); // Y1 = Y1 + X1 (account for implicit unit diagonal)
  Y1->ReciprocalMultiply(1.0, *D_, *Y1, 0.0); // y = D*y (D_ has inverse of diagonal)
  Epetra_MultiVector Y1temp(*Y1); // Need a temp copy of Y1
  U_->Multiply(true, Y1temp, *Y1);
  Y1->Update(1.0, Y1temp, 1.0); // (account for implicit unit diagonal)
  return(0);
}
//=============================================================================
int Ifpack_CrsIct::Condest(bool Trans, double & ConditionNumberEstimate) const {

  if (Condest_>=0.0) {
    ConditionNumberEstimate = Condest_;
    return(0);
  }
  // Create a vector with all values equal to one
  Epetra_Vector Ones(A_.RowMap());
  Epetra_Vector OnesResult(Ones);
  Ones.PutScalar(1.0);

  EPETRA_CHK_ERR(Solve(Trans, Ones, OnesResult)); // Compute the effect of the solve on the vector of ones
  EPETRA_CHK_ERR(OnesResult.Abs(OnesResult)); // Make all values non-negative
  EPETRA_CHK_ERR(OnesResult.MaxValue(&ConditionNumberEstimate)); // Get the maximum value across all processors
  Condest_ = ConditionNumberEstimate; // Save value for possible later calls
  return(0);
}
//=============================================================================
// Non-member functions

ostream& operator << (ostream& os, const Ifpack_CrsIct& A)
{
  // int LevelOverlap = A.LevelOverlap();
  Epetra_CrsMatrix & U = (Epetra_CrsMatrix &) A.U();
  Epetra_Vector & D = (Epetra_Vector &) A.D();

  //os.width(14);
  //os <<  "     Level of Overlap = "; os << LevelOverlap;
  //os << endl;

  os.width(14);
  os <<  "     Inverse of Diagonal = ";
  os << endl;
  os << D; // Let Epetra_Vector handle the rest.
  os << endl;

  os.width(14);
  os <<  "     Upper Triangle = ";
  os << endl;
  os << U; // Let Epetra_CrsMatrix handle the rest.
  os << endl;
 
  // Reset os flags

  return os;
}
