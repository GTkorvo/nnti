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

/*!
 * \file Amesos_Klu.h
 *
 * \class Amesos_Klu
 *
 * \brief Interface to KLU internal solver.Interface to KLU internal solver.
 *
 * \date Last updated on 24-May-05.
 */

#ifndef AMESOS_KLU_H
#define AMESOS_KLU_H

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_NoCopiable.h"
#include "Amesos_Utils.h"
#include "Amesos_Time.h"
#include "Amesos_Status.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Epetra_Import.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"


//! Amesos_Klu:  A serial, unblocked code ideal for getting started and for very sparse matrices, such as circuit matrces.

/*! 

Class Amesos_Klu is an object-oriented wrapper for KLU. KLU, whose sources
are distributed
within Amesos, is a serial solver for sparse matrices. KLU will solve a 
linear system of equations: \f$A X = B\f$, where
<TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
Epetra_MultiVector objects.

Amesos_Klu computes \f$A^T X = B\f$ 
more efficiently than \f$>A X = B\f$.  The
latter requires a matrix transpose -- which costs both time and space.

KLU is Tim Davis' implementation of Gilbert-Peierl's left-looking
sparse partial pivoting algorithm, with Eisenstat & Liu's symmetric
pruning.  Gilbert's version appears as \c [L,U,P]=lu(A) in MATLAB.
It doesn't exploit dense matrix kernels, but it is the only sparse
LU factorization algorithm known to be asymptotically optimal,
in the sense that it takes time proportional to the number of
floating-point operations.  It is the precursor to SuperLU,
thus the name ("clark Kent LU").  For very sparse matrices that
do not suffer much fill-in (such as most circuit matrices when
permuted properly) dense matrix kernels do not help, and the
asymptotic run-time is of practical importance.

The \c klu_btf code first permutes the matrix to upper block
triangular form (using two algorithms by Duff and Reid,
MC13 and MC21, in the ACM Collected Algorithms).  It then permutes
each block via a symmetric minimum degree ordering (AMD, by Amestoy,
Davis, and Duff).  This ordering phase can be done just once
for a sequence of matrices.  Next, it factorizes each reordered
block via the klu routine, which also attempts to preserve
diagonal pivoting, but allows for partial pivoting if the diagonal
is to small.    

*/

// Amesos_Klu_Pimpl contains a pointer to two structures defined in 
// klu.h:  klu_symbolic and klu_numeric.  This prevents Amesos_Klu.h 
// from having to include klu.h.
//
//  Doxygen does not handle forward class references well.
#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Amesos_Klu_Pimpl ; 
#endif

class Amesos_Klu: public Amesos_BaseSolver,  
                  private Amesos_Time, 
                  private Amesos_NoCopiable, 
                  private Amesos_Utils, 
                  private Amesos_Status { 

public: 

  //@{ \name Constructors and Destructors
  //! Amesos_Klu Constructor.
  /*! Creates an Amesos_Klu instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Klu(const Epetra_LinearProblem& LinearProblem );

  //! Amesos_Klu Destructor.
  ~Amesos_Klu(void);
  
  //@}
  //@{ \name Mathematical functions.

  int SymbolicFactorization() ;

  int NumericFactorization() ;

  int Solve();

  //@}
  //@{ \name 

  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Returns true if KLU can handle this matrix shape 
  /*! Returns true if the matrix shape is one that KLU can
    handle. KLU only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  //! SetUseTranpose(true) is more efficient in Amesos_Klu
  /*! 
    If SetUseTranspose() is set to true, 
    \f$A^T X = B\f$ is computed.
  */  
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

  bool UseTranspose() const {return(UseTranspose_);};

  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};

  int SetParameters( Teuchos::ParameterList &ParameterList );

  //! Prints timing information
  void PrintTiming() const;
  
  //! Prints information about the factorization and solution phases.
  void PrintStatus() const;
  
private:  
  
  //@}
  //@{ \name Utility methods

  /*
  CreateLocalMatrixAndExporters - Prepare to convert matrix and vectors to serial 
    Preconditions:
      Problem_ must be set 
      SerialMap and SerialCrsMatrix must either be 0 or be pointers to 
        appropriatly allocate objects.  If they are non-zero, those objects
	will be deleted (and possibly recreated).  
	
    Postconditions:
      UseDataInPlace_ is set to 1 if the input matrix can be used in place, i.e.
        1)  is entirely stored on process 0
        2)  is storage optimized
        3)  AddToDiag_ is not set
      SerialMap points to a serial map if UseDataInPlace_==1
      SerialCrsMatrix contains a serial version of the matrix A if UseDataInPlace_==1
      SerialMatrix points to a serial copy of the matrix
      NumGlobalElements_   is set to the number of rows in the matrix
      numentries_ is set to the number of non-zeroes in the matrix 
   */
  int CreateLocalMatrixAndExporters() ;
  /*
    ExportToSerial
    Preconditions:
       UseDataInPlace_ must be set
       ImportToSerial and SerialCrsMatrixA_ must be set if UseDataInPlace_ != 1
       AddToDiag_ 
    Postconditions
       SerialMatrix_ points to a serial version of the matrix
         With AddToDiag_ added to the diagonal if AddToDiag_ is non-zero 
   */
  int ExportToSerial() ;
  /*
    ConvertToKluCRS - Convert matrix to form expected by Klu: Ai, Ap, Aval
    Preconditions:
      numentries_, NumGloalElements_ and SerialMatrix_ must be set.
    Postconditions:
      Ai, Ap, and Aval are resized and populated with a compresses row storage 
      version of the input matrix A.
  */
  int ConvertToKluCRS(bool firsttime);     

  /*
    PerformSymbolicFactorization - Call Klu to perform symbolic factorization
    Preconditions:
      UseDataInPlace_ must be set to 1 if the input matrix is entirely stored on process 0
      Ap, Ai and Aval are a compressed row storage version of the input matrix A.
    Postconditions:
      Symbolic points to an KLU internal opaque object containing the
        symbolic factorization and accompanying information.  
      SymbolicFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
      
  int PerformSymbolicFactorization(); 

  /*
    PerformNumericFactorization - Call Klu to perform numeric factorization
    Preconditions:
      UseDataInPlace_ must be set 
      Ap, Ai and Aval are a compressed row storage version of the input matrix A.
      Symbolic must be set
    Postconditions:
      Numeric points to an KLU internal opaque object containing the
        numeric factorization and accompanying information.  
      NumericFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
  int PerformNumericFactorization(); 

  // @}
  
  //! Creates SerialMap_
  int CreateSerialMap();

  int *Lp, *Li, *Up, *Ui, *P ;	
  double *Lx, *Ux ;
  Amesos_Klu_Pimpl *PrivateKluData_; 

  //! Ap, Ai, Aval form the compressed row storage used by Klu
  vector <int> Ap;
  vector <int> Ai;
  vector <double> Aval;

  //! Process number (i.e. Comm().MyPID()).
  int iam;
  //! Number of processes in computation.
  int NumProcs_;
  //! 1 if Problem_->GetOperator() is stored entirely on process 0
  int UseDataInPlace_;
  //! Number of non-zero entries in Problem_->GetOperator()
  int numentries_;
  //! Number of rows and columns in the Problem_->GetOperator()
  int NumGlobalElements_;

  //! Operator converted to a RowMatrix
  Epetra_RowMatrix *RowMatrixA_;
  //! Points to a Serial Map (unused if UseDataInPlace_ == 1 )
  Epetra_Map *SerialMap_;
  //! Points to a Serial Copy of A (unused if UseDataInPlace_==1)
  Epetra_CrsMatrix *SerialCrsMatrixA_;
  //! Points to a Serial Copy of A 
  Epetra_RowMatrix *SerialMatrix_ ; 
  //! Points to the matrix that is used to compute
  Epetra_RowMatrix *Matrix_;

  //! If \c true, the transpose of A is used.
  bool UseTranspose_;
  //! Pointer to the linear system problem.
  const Epetra_LinearProblem * Problem_;

  //! Only used for RowMatrices to extract copies.
  vector<int>ColIndicesV_;
  //! Only used for RowMatrices to extract copies.
  vector<double>RowValuesV_;

  bool refactorize_;	    // if true, and if the Symbolic and Numeric
			    // objects have already been created, then
			    // attempt to "refactorize" (factor the matrix
			    // with no changes to the pivot order since the
			    // last call the klu_btf_factor).

  double rcond_threshold_;  // if we refactorize, the factorization may suffer
			    // in numeric quality.  We compute rcond =
			    // min (abs (diag (U))) / max (abs (diag (U))).
			    // If this ratio is <= rcond_threshold_, then
			    // the "refactorization" is scrapped, and we factor
			    // with full partial pivoting instead.

  int ScaleMethod_;	    // most methods (KLU, UMFPACK, Mumps, ...) can scale
			    // the input matrix prior to factorization.  This can
			    // improve pivoting, reduce fill-in, and lead to a
			    // better quality factorization.  The options are:
			    // 0: no scaling
			    // 1: use the default method for the specific package
			    // 2: use the method's 1st alternative (if it has one)
			    // 3: use the method's 2nd alternative, and so on.

  //! Importer to process 0.
  Epetra_Import * ImportToSerial_;
  
};  // class Amesos_Klu  

#endif /* AMESOS_KLU_H */
