
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_SerialSymDenseMatrix.h"
#include "Epetra_Util.h"

//=============================================================================
Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(bool set_object_label)
  : Epetra_CompObject(),
    Epetra_Object(-1, false),
    M_(0),
    N_(0),
    A_Copied_(false),
    CV_(Copy),
    LDA_(0),
    A_(0),
    UseTranspose_(false)
{
  if (set_object_label) {
    SetLabel("Epetra::SerialDenseMatrix");
  }
}

//=============================================================================
Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(int NumRows, int NumCols,
                                                   bool set_object_label)
  : Epetra_CompObject(),
    Epetra_Object(-1, false),
    M_(0),
    N_(0),
    A_Copied_(false),
    CV_(Copy),
    LDA_(0),
    A_(0),
    UseTranspose_(false)
{
  if (set_object_label) {
    SetLabel("Epetra::SerialDenseMatrix");
  }
  if(NumRows < 0)
	throw ReportError("NumRows = " + toString(NumRows) + ". Should be >= 0", -1);
  if(NumCols < 0)
	throw ReportError("NumCols = " + toString(NumCols) + ". Should be >= 0", -1);

  int errorcode = Shape(NumRows, NumCols);
  if(errorcode != 0)
    throw ReportError("Shape returned non-zero value", errorcode);
}

//=============================================================================
Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(Epetra_DataAccess CV, double* A, int LDA,
                                                   int NumRows, int NumCols,
                                                   bool set_object_label)
  : Epetra_CompObject(),
    Epetra_Object(-1, false),
    M_(NumRows),
    N_(NumCols),
    A_Copied_(false),    
    CV_(CV),
    LDA_(LDA),
    A_(A),
    UseTranspose_(false)
{
  if (set_object_label) {
    SetLabel("Epetra::SerialDenseMatrix");
  }
  if(A == 0)
	throw ReportError("Null pointer passed as A parameter.", -3);
  if(NumRows < 0)
    throw ReportError("NumRows = " + toString(NumRows) + ". Should be >= 0", -1);
  if(NumCols < 0)
	throw ReportError("NumCols = " + toString(NumCols) + ". Should be >= 0", -1);
  if(LDA < 0)
	throw ReportError("LDA = " + toString(LDA) + ". Should be >= 0", -1);

  if (CV == Copy) {
    LDA_ = M_;
    const int newsize = LDA_ * N_;
    if (newsize > 0) {
      A_ = new double[newsize];
      CopyMat(A, LDA, M_, N_, A_, LDA_);
      A_Copied_ = true;
    }
    else {
      A_ = 0;
    }
  }
}
//=============================================================================
Epetra_SerialDenseMatrix::Epetra_SerialDenseMatrix(const Epetra_SerialDenseMatrix& Source)
  : Epetra_CompObject(Source),  
    M_(Source.M_),
    N_(Source.N_),
    A_Copied_(false),
    CV_(Source.CV_),
    LDA_(Source.LDA_),
    A_(Source.A_),
    UseTranspose_(false)
{
	SetLabel(Source.Label());
	if(CV_ == Copy) {
		LDA_ = M_;
		const int newsize = LDA_ * N_;
		if(newsize > 0) {
			A_ = new double[newsize];
			CopyMat(Source.A_, Source.LDA_, M_, N_, A_, LDA_);
			A_Copied_ = true;
		}
		else {
			A_ = 0;
		}
	}
}

//=============================================================================
int Epetra_SerialDenseMatrix::Reshape(int NumRows, int NumCols) {
	if(NumRows < 0 || NumCols < 0)
		return(-1);

	double* A_tmp = 0;
	const int newsize = NumRows * NumCols;

	if(newsize > 0) {
		// Allocate space for new matrix
		A_tmp = new double[newsize];
		for (int k = 0; k < newsize; k++) 
			A_tmp[k] = 0.0; // Zero out values
		int M_tmp = EPETRA_MIN(M_, NumRows);
		int N_tmp = EPETRA_MIN(N_, NumCols);
		if (A_ != 0) 
			CopyMat(A_, LDA_, M_tmp, N_tmp, A_tmp, NumRows); // Copy principal submatrix of A to new A
  }
  CleanupData(); // Get rid of anything that might be already allocated  
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
	if(newsize > 0) {
		A_ = A_tmp; // Set pointer to new A
		A_Copied_ = true;
	}

  return(0);
}
//=============================================================================
int Epetra_SerialDenseMatrix::Shape(int NumRows, int NumCols) {
	if(NumRows < 0 || NumCols < 0)
		return(-1);

  CleanupData(); // Get rid of anything that might be already allocated
  M_ = NumRows;
  N_ = NumCols;
  LDA_ = M_;
	const int newsize = LDA_ * N_;
	if(newsize > 0) {
		A_ = new double[newsize];
		for (int k = 0; k < newsize; k++) 
			A_[k] = 0.0; // Zero out values
		A_Copied_ = true;
	}

  return(0);
}
//=============================================================================
Epetra_SerialDenseMatrix::~Epetra_SerialDenseMatrix()
{
  CleanupData();
}
//=============================================================================
void Epetra_SerialDenseMatrix::CleanupData()
{
  if (A_Copied_)
    delete [] A_;
	A_ = 0;
	A_Copied_ = false;
	M_ = 0;
	N_ = 0;
	LDA_ = 0;
}
//=============================================================================
Epetra_SerialDenseMatrix& Epetra_SerialDenseMatrix::operator = (const Epetra_SerialDenseMatrix& Source) {
  if(this == &Source)
		return(*this); // Special case of source same as target
	if((CV_ == View) && (Source.CV_ == View) && (A_ == Source.A_))
		return(*this); // Special case of both are views to same data.

	if(strcmp(Label(), Source.Label()) != 0)
		throw ReportError("operator= type mismatch (lhs = " + string(Label()) + 
											", rhs = " + string(Source.Label()) + ").", -5);
	
	if(Source.CV_ == View) {
		if(CV_ == Copy) { // C->V only
			CleanupData();
			CV_ = View;
		}
		M_ = Source.M_; // C->V and V->V
		N_ = Source.N_;
		LDA_ = Source.LDA_;
		A_ = Source.A_;
	}
	else {
		if(CV_ == View) { // V->C
			CV_ = Copy;
			M_ = Source.M_;
			N_ = Source.N_;
			LDA_ = Source.M_;
			const int newsize = LDA_ * N_;
			if(newsize > 0) {
				A_ = new double[newsize];
				A_Copied_ = true;
			}
			else {
				A_ = 0;
				A_Copied_ = false;
			}
		}
		else { // C->C
			if((Source.M_ <= LDA_) && (Source.N_ == N_)) { // we don't need to reallocate
				M_ = Source.M_;
				N_ = Source.N_;
			}
			else { // we need to allocate more space (or less space)
				CleanupData();
				M_ = Source.M_;
				N_ = Source.N_;
				LDA_ = Source.M_;
				const int newsize = LDA_ * N_;
				if(newsize > 0) {
					A_ = new double[newsize];
					A_Copied_ = true;
				}
			}
		}
		CopyMat(Source.A_, Source.LDA_, M_, N_, A_, LDA_); // V->C and C->C
	}
	
  return(*this);
}

//=============================================================================
Epetra_SerialDenseMatrix& Epetra_SerialDenseMatrix::operator+= ( const Epetra_SerialDenseMatrix & Source) {
  if (M() != Source.M()) 
		throw ReportError("Row dimension of source = " + toString(Source.M()) +
											" is different than  row dimension of target = " + toString(LDA()), -1);
  if (N() != Source.N()) 
		throw ReportError("Column dimension of source = " + toString(Source.N()) +
											" is different than column dimension of target = " + toString(N()), -2);

  CopyMat(Source.A(), Source.LDA(), Source.M(), Source.N(), A(), LDA(), true);
  return(*this);
}
//=============================================================================
void Epetra_SerialDenseMatrix::CopyMat(double* Source,
                                       int Source_LDA,
                                       int NumRows,
                                       int NumCols, 
                                       double* Target,
                                       int Target_LDA,
                                       bool add)
{
  int i;
  double* tptr = Target;
  double* sptr = Source;
  if (add) {
    for(int j=0; j<NumCols; ++j) {
      for(i=0; i<NumRows; ++i) {
        tptr[i] += sptr[i];
      }

      tptr += Target_LDA;
      sptr += Source_LDA;
    }
  }
  else {
    for(int j=0; j<NumCols; ++j) {
      for(i=0; i<NumRows; ++i) {
        tptr[i] = sptr[i];
      }

      tptr += Target_LDA;
      sptr += Source_LDA;
    }
  }
/*
  double* targetPtr = Target;
  double* sourcePtr = Source;
  double* targetEnd = 0;
 
  for (j = 0; j < NumCols; j++) {
    targetPtr = Target + j * Target_LDA;
    targetEnd = targetPtr+NumRows;
    sourcePtr = Source + j * Source_LDA;
    if (add) {
      while(targetPtr != targetEnd) 
		*targetPtr++ += *sourcePtr++;
    }
    else {
      while(targetPtr != targetEnd) 
		*targetPtr++ = *sourcePtr++;
    }
  }
*/
  return;
}
//=============================================================================
double Epetra_SerialDenseMatrix::NormOne() const {

  int i, j;

    double anorm = 0.0;
    double * ptr;
    for (j=0; j<N_; j++) {
      double sum=0.0;
      ptr = A_ + j*LDA_;
      for (i=0; i<M_; i++) sum += fabs(*ptr++);
      anorm = EPETRA_MAX(anorm, sum);
    }
    UpdateFlops(N_*N_);
    return(anorm);
}
//=============================================================================
double Epetra_SerialDenseMatrix::NormInf() const {

  int i, j;

    double anorm = 0.0;
    double * ptr;

    // Loop across columns in inner loop.  Most expensive memory access, but 
    // requires no extra storage.
    for (i=0; i<M_; i++) {
      double sum=0.0;
      ptr = A_ + i;
      for (j=0; j<N_; j++) {
	sum += fabs(*ptr);
	ptr += LDA_;
      }
      anorm = EPETRA_MAX(anorm, sum);
    }
    UpdateFlops(N_*N_);
    return(anorm);
}
//=============================================================================
int Epetra_SerialDenseMatrix::Scale(double ScalarA) {

  int i, j;
  
  double * ptr;
  for (j=0; j<N_; j++) {
    ptr = A_ + j*LDA_;
    for (i=0; i<M_; i++) { *ptr = ScalarA * (*ptr); ptr++; }
  }
  UpdateFlops(N_*N_);
  return(0);
  
}
//=========================================================================
int  Epetra_SerialDenseMatrix::Multiply (char TransA, char TransB, double ScalarAB, 
				      const Epetra_SerialDenseMatrix& A, 
				      const Epetra_SerialDenseMatrix& B,
				      double ScalarThis ) {
  // Check for compatible dimensions

  if (TransA!='T' && TransA!='N') EPETRA_CHK_ERR(-2); // Return error
  if (TransB!='T' && TransB!='N') EPETRA_CHK_ERR(-3);
  
  int A_nrows = (TransA=='T') ? A.N() : A.M();
  int A_ncols = (TransA=='T') ? A.M() : A.N();
  int B_nrows = (TransB=='T') ? B.N() : B.M();
  int B_ncols = (TransB=='T') ? B.M() : B.N();
  
  if (M_        != A_nrows     ||
      A_ncols   != B_nrows     ||
      N_        != B_ncols  ) EPETRA_CHK_ERR(-1); // Return error

    
  // Call GEMM function
  GEMM(TransA, TransB, M_, N_, A_ncols, ScalarAB, A.A(), A.LDA(), 
       B.A(), B.LDA(), ScalarThis, A_, LDA_);
  long int nflops = 2*M_;
  nflops *= N_;
  nflops *= A_ncols;
  if (ScalarAB != 1.0) nflops += M_*N_;
  if (ScalarThis != 0.0) nflops += M_*N_;
  UpdateFlops(nflops);

  return(0);
}
//=========================================================================
int  Epetra_SerialDenseMatrix::Multiply (char SideA, double ScalarAB, 
				      const Epetra_SerialSymDenseMatrix& A, 
				      const Epetra_SerialDenseMatrix& B,
				      double ScalarThis ) {
  // Check for compatible dimensions
  
  if (SideA=='R') {
    if (M_ != B.M() || 
	N_ != A.N() ||
	B.N() != A.M() ) EPETRA_CHK_ERR(-1); // Return error
  }
  else if (SideA=='L') {
    if (M_ != A.M() || 
	N_ != B.N() ||
	A.N() != B.M() ) EPETRA_CHK_ERR(-1); // Return error
  }
  else {
    EPETRA_CHK_ERR(-2); // Return error, incorrect value for SideA
  }
  
  // Call SYMM function
  SYMM(SideA, A.UPLO(), M_, N_, ScalarAB, A.A(), A.LDA(), 
       B.A(), B.LDA(), ScalarThis, A_, LDA_);
  long int nflops = 2*M_;
  nflops *= N_;
  nflops *= A.N();
  if (ScalarAB != 1.0) nflops += M_*N_;
  if (ScalarThis != 0.0) nflops += M_*N_;
  UpdateFlops(nflops);

  return(0);
}
//=========================================================================
void Epetra_SerialDenseMatrix::Print(ostream& os) const {
	os << endl;
	if(CV_ == Copy)
		os << "Data access mode: Copy" << endl;
	else
		os << "Data access mode: View" << endl;
	if(A_Copied_)
		os << "A_Copied: yes" << endl;
	else
		os << "A_Copied: no" << endl;
	os << "Rows(M): " << M_ << endl;
	os << "Columns(N): " << N_ << endl;
	os << "LDA: " << LDA_ << endl;
	if(M_ == 0 || N_ == 0)
		os << "(matrix is empty, no values to display)" << endl;
	else
		for(int i = 0; i < M_; i++) {
			for(int j = 0; j < N_; j++){
				os << (*this)(i,j) << " ";
			}
			os << endl;
		}
}
//=========================================================================
int Epetra_SerialDenseMatrix::Random() {

	Epetra_Util util;

	for(int j = 0; j < N_; j++) {
		double* arrayPtr = A_ + (j * LDA_);
		for(int i = 0; i < M_; i++) {
			*arrayPtr++ = util.RandomDouble();
		}
	}
	
  return(0);
}
  
