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

#ifndef AMESOS_LAPACK_H
#define AMESOS_LAPACK_H

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_NoCopiable.h"
#include "Amesos_Utils.h"
#include "Amesos_Time.h"
#include "Amesos_Status.h"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_SerialDenseSolver.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
class Epetra_RowMatrix;
class Epetra_LinearProblem;
#include "Teuchos_RefCountPtr.hpp"
using namespace Teuchos;

//! Amesos_Lapack: an interface to LAPACK.
/*!
Class Amesos_Lapack enables the solution of the distributed linear
system, defined by an Epetra_LinearProblem, using LAPACK.

Amesos_Lapack stores the lineaar system matrix as an
Epetra_SerialDensMatrix. The linear problem is an Epetra_SerialDenseProblem.
Amesos_Lapack factorizes the matrix using DGETRF().

\date Last updated on 16-Mar-05.

\author Marzio Sala, 9214.

*/
class Amesos_Lapack: public Amesos_BaseSolver,
                     private Amesos_Time,
                     private Amesos_NoCopiable,
                     private Amesos_Utils,
                     private Amesos_Status 
{
public: 

  //@{ \name Constructor methods
  //! Amesos_Lapack Constructor.
  /*! Creates an Amesos_Lapack instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Lapack(const Epetra_LinearProblem& LinearProblem );

  //! Amesos_Lapack Destructor.
  /*! Completely deletes an Amesos_Lapack object.  
  */
  ~Amesos_Lapack(void);
  
  //@}
  //@{ \name Mathematical functions.

  int SymbolicFactorization() ;

  int NumericFactorization() ;

  int Solve();

  //@}
  //@{ \name Additional methods required to support the Epetra_Operator interface.

  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  bool MatrixShapeOK() const;

  int SetUseTranspose(bool UseTranspose) {
    UseTranspose_ = UseTranspose; 
    return(0);
  }

  bool UseTranspose() const {return(UseTranspose_);};

  const Epetra_Comm & Comm() const {
    return(GetProblem()->GetOperator()->Comm());
  }

  int SetParameters( Teuchos::ParameterList &ParameterList )  ;

  //! Computes the eigenvalues of the linear system matrix using DGEEV.
  /*!
    \param Er - (Out) On processor zero only, it will contain the 
                      real component of the eigenvalues.

    \param Ei - (Out) On processor zero only, it will contain the 
                      imaginary component of the eigenvalues.

    \note Er and Ei must have been allocated so that the local
    length on processor 0 equals the global size of the matrix.
    */
  int GEEV(Epetra_Vector& Er, Epetra_Vector& Ei);

  //! Print timing information
  void PrintTiming();
  
  //! Print information about the factorization and solution phases.
  void PrintStatus();
  
  //@}

protected:

  //! Returns true if SymbolicFactorization() has been successfully called.
  bool IsSymbolicFactorizationOK()
  {
    return(IsSymbolicFactorizationOK_);
  }

  //! Returns true if SymbolicFactorization() has been successfully called.
  bool IsNumericFactorizationOK()
  {
    return(IsNumericFactorizationOK_);
  }
  
  //! Returns a pointer to the linear system matrix.
  const Epetra_RowMatrix* Matrix() const
  {
    return(Problem_->GetMatrix());
  }

  //! Returns the number of global rows, or -1 if Matrix() returns 0.
  inline int NumGlobalRows() const
  {
    return(Matrix()->NumGlobalRows());
  }

  //! Returns the number of local rows, or -1 if Matrix() returns 0.
  inline int NumMyRows() const
  {
    return(Matrix()->NumMyRows());
  }

  //! Returns a reference to serial map (that with all elements on process 0). Builds SerialMap_ if necessary or required.
  inline const Epetra_Map& SerialMap()
  {
    return(*(SerialMap_.get()));
  }

  //! Returns a reference to serial matrix (that with all rows on process 0). Builds SerialMap_ if necessary or required.
  inline Epetra_RowMatrix& SerialMatrix()
  {
    return(*(SerialMatrix_.get()));
  }

  inline Epetra_CrsMatrix& SerialCrsMatrix()
  {
    return(*(SerialCrsMatrix_.get()));
  }

  //! Returns a reference to the importer map. Builds SerialMap_ if necessary or required.
  const Epetra_Import& Importer()
  {
    return(*(Importer_.get()));
  }

  //! Solves the linear system, when only one process is used.
  int SolveSerial(Epetra_MultiVector& X,
		  const Epetra_MultiVector& B);
  
  //! Solves the linear system, when more than one process is used.
  int SolveDistributed(Epetra_MultiVector& X,
		       const Epetra_MultiVector& B);
  
  //! Converts a distributed matrix to serial matrix.
  int DistributedToSerial();

  //! Converts a serial matrix to dense format.
  int SerialToDense();

  //! Factors the matrix using LAPACK.
  int DenseToFactored();

  RefCountPtr<Epetra_RowMatrix> SerialMatrix_;
  RefCountPtr<Epetra_CrsMatrix> SerialCrsMatrix_;
  RefCountPtr<Epetra_Map> SerialMap_;
  RefCountPtr<Epetra_Import> Importer_;
  //! Dense matrix.
  Epetra_SerialDenseMatrix DenseMatrix_;
  //! Dense LHS.
  Epetra_SerialDenseMatrix DenseLHS_;
  //! Dense RHS.
  Epetra_SerialDenseMatrix DenseRHS_;
  //! Linear problem for dense matrix and vectors.
  Epetra_SerialDenseSolver DenseSolver_;

  //! If \c true, the linear system with the transpose will be solved.
  bool UseTranspose_;
  //! Pointer to the linear problem.
  const Epetra_LinearProblem* Problem_;

  //! Number of calls to SymbolicFactorization().
  int NumSymbolicFact_;
  //! Number of calls to NumericFactorization().
  int NumNumericFact_;
  //! Number of calls to Solver().
  int NumSolve_;

};  // End of  class Amesos_Lapack
#endif /* AMESOS_LAPACK_H */
