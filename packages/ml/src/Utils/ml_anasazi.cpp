#include "ml_common.h"

#if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_ANASAZI)

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_Time.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Teuchos_ParameterList.hpp"

#include "AnasaziPetraInterface.hpp"
#include "AnasaziBlockArnoldi.hpp"
#include "AnasaziConfigDefs.hpp"
#include "AnasaziReturnType.hpp"
#include "AnasaziMultiVec.hpp"
#include "AnasaziMatrix.hpp"
#include "AnasaziOperator.hpp"

#include "ml_anasazi.h"

using namespace Anasazi;

// ================================================ ====== ==== ==== == =

enum MLMatOp { A_MATRIX, I_MINUS_A_MATRIX, A_PLUS_AT_MATRIX, A_MINUS_AT_MATRIX };

// ================================================ ====== ==== ==== == =

template <class TYPE> 
class MLMat : public virtual Matrix<TYPE> {
public:
  MLMat(const Epetra_RowMatrix &, const MLMatOp MatOp, const bool UseScaling );
  ~MLMat();
  ReturnType ApplyMatrix ( const MultiVec<TYPE>& x, 
			   MultiVec<TYPE>& y ) const;
  void SetMatOp(const MLMatOp MatOp) 
  {
    MatOp_ = MatOp;
  }

  void SetScaling(const double Scale) 
  {
    Scale_ = Scale;
  }
  
private:
  const Epetra_RowMatrix & Mat_;
  Epetra_Vector * Diagonal_;
  MLMatOp MatOp_;
  const bool UseScaling_;
  double Scale_;
  
};

// ================================================ ====== ==== ==== == =

template <class TYPE>
MLMat<TYPE>::MLMat(const Epetra_RowMatrix & Matrix,
		   const MLMatOp MatOp, const bool UseScaling ) 
  : Mat_(Matrix),
    MatOp_(MatOp),
    UseScaling_(UseScaling),
    Scale_(1.0)
{
  if( UseScaling_ ) {
    
    Epetra_Vector * Diagonal_ = new Epetra_Vector(Matrix.RowMatrixRowMap());
    Matrix.ExtractDiagonalCopy(*Diagonal_);
    
    int NumMyElements = Matrix.RowMatrixRowMap().NumMyElements();
    for( int i=0 ; i<NumMyElements ; ++i ) {
      if( (*Diagonal_)[i] != 0.0 ) (*Diagonal_)[i] = 1.0/(*Diagonal_)[i];
      else  (*Diagonal_)[i] = 0.0;
    }
    
  }
  
}

// ================================================ ====== ==== ==== == =

template <class TYPE>
MLMat<TYPE>::~MLMat() 
{
  if( UseScaling_ ) delete Diagonal_;
}

// ================================================ ====== ==== ==== == =

template <class TYPE>
ReturnType MLMat<TYPE>::ApplyMatrix ( const MultiVec<TYPE>& x, 
				      MultiVec<TYPE>& y ) const 
{

  int info = 0;
  MultiVec<TYPE> & temp_x = const_cast<MultiVec<TYPE> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);

  if( vec_x==0 || vec_y==0 ) return Failed;

  info = const_cast<Epetra_RowMatrix &>(Mat_).Apply( *vec_x, *vec_y );
  vec_y->Scale(Scale_);

  if( info ) return Failed;

  // FIXME: put me in constructor....
  Epetra_MultiVector * tmp ;
  if( MatOp_ == A_PLUS_AT_MATRIX || MatOp_ == A_MINUS_AT_MATRIX ) {
    tmp = new Epetra_MultiVector(*vec_x);
    const Epetra_RowMatrix & RowMatrix = dynamic_cast<const Epetra_RowMatrix &>(Mat_);
    info = RowMatrix.Multiply(true,*vec_x,*tmp);
    tmp->Scale(Scale_);
    if( info ) return Failed;
  }
  if( MatOp_ == A_PLUS_AT_MATRIX ) vec_y->Update(1.0,*tmp,1.0);
  else if( MatOp_ == A_MINUS_AT_MATRIX ) vec_y->Update(-1.0,*tmp,1.0);

  if(  MatOp_ == A_PLUS_AT_MATRIX || MatOp_ == A_MINUS_AT_MATRIX ) delete tmp;
  
  // diagonal scaling
  if( UseScaling_ == true ) {
    for( int j=0 ; j<vec_y->NumVectors() ; ++j ) 
      for( int i=0 ; i<vec_y->Map().NumMyElements() ; ++i ) (*vec_y)[j][i] *= (*Diagonal_)[i];
  }
  
  if( MatOp_ == I_MINUS_A_MATRIX ) {
    vec_y->Update (1.0,*vec_x,-1.0);
  }

  return Ok; 
}

// ================================================ ====== ==== ==== == =

int ML_Anasazi_Interface(const Epetra_RowMatrix * RowMatrix, Epetra_MultiVector & EigenVectors,
			 double RealEigenvalues[], double ImagEigenvalues[],
			 Teuchos::ParameterList & List) 
{

  int MyPID = RowMatrix->Comm().MyPID();
  
  /* ********************************************************************** */
  /* Retrive parameters' choices                                            */
  /* ********************************************************************** */

  MLMatOp MatOp;
  string MatOpStr = "A";
  MatOpStr = List.get("eigen-analysis: matrix operation", MatOpStr);
  if( MatOpStr == "A" ) MatOp = A_MATRIX;
  else if( MatOpStr == "I-A" ) MatOp = I_MINUS_A_MATRIX;
  else if( MatOpStr == "A+A^T" ) MatOp = A_PLUS_AT_MATRIX;
  else if( MatOpStr == "A-A^T" ) MatOp = A_MINUS_AT_MATRIX;

  bool UseScaling = true;
  UseScaling = List.get("eigen-analysis: use scaling", UseScaling);

  int length = List.get("eigen-analysis: length", 20);
  double tol = List.get("eigen-analysis: tolerance", 1.0e-5);
  string which = List.get("eigen-analysis: action", "LM");
  int restarts = List.get("eigen-analysis: restart", 100);
  bool isSymmetric = List.get("eigen-analysis: symmetric problem", false);

  int output = List.get("eigen-analysis: output", 5);
  
  double Scaling = List.get("eigen-analysis: scaling", 1.0);
  
  if( output > 5 && MyPID == 0 ) {
    if( MatOp == A_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of A" << endl;
    if( MatOp == I_MINUS_A_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of I - A" << endl;
    if( MatOp == A_PLUS_AT_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of A + A^T" << endl;
    if( MatOp == A_MINUS_AT_MATRIX ) cout << "ML_Anasazi : Computing eigenvalues of A - A^T" << endl;
    if( UseScaling ) cout << "ML_Anasazi : where A is scaled by D^{-1}" << endl;
    if( isSymmetric ) cout << "ML_Anasazi : Problem is symmetric" << endl;
    cout << "ML_Anasazi : Tolerance = " << tol << endl;
    cout << "ML_Anasazi : Scaling = " << Scaling << endl;
    cout << "ML_Anasazi : Required Action = " << which << endl;
	
  }

  bool PrintCurrentStatus =  List.get("eigen-analysis: print current status", false);

  /* ********************************************************************** */
  /* view mode for PetraVectors                                             */
  /* ********************************************************************** */
  
  int NumBlocks = EigenVectors.NumVectors();

  int indices[NumBlocks];
  for( int i=0 ; i<NumBlocks ; ++i ) indices[i] = i;
  
  Anasazi::PetraVec<double> Vectors(View,EigenVectors,indices,NumBlocks);
  
  /* ********************************************************************** */
  /* Perform required action                                                */
  /* ********************************************************************** */

  int step = restarts*length*NumBlocks;
  
  // Call the ctor that calls the petra ctor for a matrix
  MLMat<double> Amat(*RowMatrix,MatOp,UseScaling);
  Amat.SetScaling(Scaling);
  
  Anasazi::Eigenproblem<double> MyProblem(&Amat, &Vectors);

  // Initialize the Block Arnoldi solver
  Anasazi::BlockArnoldi<double> MyBlockArnoldi1(MyProblem, tol, NumBlocks, length, NumBlocks, 
						which, step, restarts);
	
  // Inform the solver that the problem is symmetric
  if( isSymmetric ) MyBlockArnoldi1.setSymmetric(true);
  MyBlockArnoldi1.setDebugLevel(0);

  //  MyBlockArnoldi.iterate(5);
  
  // Solve the problem to the specified tolerances or length
  MyBlockArnoldi1.solve();
  
  // Obtain results directly

  double * evalr = MyBlockArnoldi1.getEvals(); 
  double * evali = MyBlockArnoldi1.getiEvals();

  for( int i=0 ; i<NumBlocks ; ++i ) {
    RealEigenvalues[i] = evalr[i];
    ImagEigenvalues[i] = evali[i];
  }
  
  // Output results to screen
  if( PrintCurrentStatus > 5 && MyPID == 0 ) MyBlockArnoldi1.currentStatus();

  MyBlockArnoldi1.getEvecs(Vectors);
  
  /* ********************************************************************** */
  /* Scale vectors (real and imag part) so that max is 1 for each of them   */
  /* ********************************************************************** */

  if( List.get("eigen-analysis: normalize eigenvectors",false) ) {
    
    double MaxVals[NumBlocks];
    double MinVals[NumBlocks];
    double Vals[NumBlocks];
    EigenVectors.MaxValue(MaxVals);
    EigenVectors.MinValue(MinVals);

    for( int i=0 ; i<NumBlocks ; ++i ) {
      Vals[i] = EPETRA_MAX(abs(MaxVals[i]),abs(MinVals[i]));
      if( Vals[i] == abs(MaxVals[i]) && MaxVals[i]<0 ) Vals[i] *= -1.0;
      if( Vals[i] == abs(MinVals[i]) && MinVals[i]<0 ) Vals[i] *= -1.0;
      Vals[i] = 1.0/Vals[i];
      EigenVectors(i)->Scale(Vals[i]);

    }
  }

  return 0;
  
}

// ================================================ ====== ==== ==== == =

int ML_Anasazi_Get_FiledOfValuesBox(const Epetra_RowMatrix * RowMatrix, 
				    double & MaxReal, double & MaxImag,
				    Teuchos::ParameterList & List ) 
{

  int MyPID = RowMatrix->Comm().MyPID();

  bool UseScaling = true;
  UseScaling = List.get("field-of-values: use scaling", UseScaling);
  
  int length = List.get("field-of-values: length", 20);
  double tol = List.get("field-of-values: tolerance", 1.0e-5);
  int restarts = List.get("field-of-values: restart", 100);

  int output = List.get("output", 5);
  
  if( output > 5 && MyPID == 0 ) {
    cout << "ML_Anasazi : Estimate box containing the field of values" << endl;
    cout << "ML_Anasazi : Computing eigenvalues of A + A^T" << endl;
    if( UseScaling ) cout << "ML_Anasazi : where A is scaled by D^{-1}" << endl;
    cout << "ML_Anasazi : Tolerance = " << tol << endl;	
  }

  bool PrintCurrentStatus =  List.get("field-of-values: print current status", false);
  
  /* ********************************************************************** */
  /* First compute A + A^T to get the real bound                            */
  /* ********************************************************************** */
  
  Anasazi::PetraVec<double> Vectors(RowMatrix->RowMatrixRowMap(), 1);
  Vectors.MvRandom();
  
  int step = restarts*length*1;
  
  // Call the ctor that calls the petra ctor for a matrix
  MLMat<double> Amat(*RowMatrix,A_PLUS_AT_MATRIX,UseScaling);
  
  Anasazi::Eigenproblem<double> MyProblem(&Amat, &Vectors);

  // Initialize the Block Arnoldi solver
  Anasazi::BlockArnoldi<double> MyBlockArnoldi1(MyProblem, tol, 1, length, 1, 
						"LM", step, restarts);

  // Inform the solver that the problem is symmetric
  MyBlockArnoldi1.setSymmetric(true);
  MyBlockArnoldi1.setDebugLevel(0);

  //  MyBlockArnoldi.iterate(5);
  
  // Solve the problem to the specified tolerances or length
  MyBlockArnoldi1.solve();
  
  // Obtain results directly

  double * evalr = MyBlockArnoldi1.getEvals(); 

  MaxReal = evalr[0] / 2;

  if( PrintCurrentStatus && MyPID == 0 ) MyBlockArnoldi1.currentStatus();

  /* ********************************************************************** */
  /* First compute A - A^T to get the real bound                            */
  /* ********************************************************************** */
  
  if( output > 5 && MyPID == 0 ) {
    cout << "ML_Anasazi : Computing eigenvalues of A - A^T" << endl;
  }
  
  Vectors.MvRandom();
  
  // Call the ctor that calls the petra ctor for a matrix
  MLMat<double> Amat2(*RowMatrix,A_MINUS_AT_MATRIX,UseScaling);
  
  Anasazi::Eigenproblem<double> MyProblem2(&Amat2, &Vectors);

  // Initialize the Block Arnoldi solver
  Anasazi::BlockArnoldi<double> MyBlockArnoldi2(MyProblem2, tol, 1, length, 1, 
						"LM", step, restarts);

  // Inform the solver that the problem is symmetric
  MyBlockArnoldi2.setSymmetric(false);
  MyBlockArnoldi2.setDebugLevel(0);

  MyBlockArnoldi2.iterate(5);
  
  // Solve the problem to the specified tolerances or length
  MyBlockArnoldi2.solve();
  
  // Obtain results directly

  double * evali = MyBlockArnoldi2.getiEvals();

  MaxImag = evali[0] / 2;

  if( PrintCurrentStatus && MyPID == 0 ) MyBlockArnoldi2.currentStatus();
  
  return 0;
  
}

#endif
