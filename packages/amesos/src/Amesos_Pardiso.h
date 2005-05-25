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

#ifndef AMESOS_PARDISO_H
#define AMESOS_PARDISO_H

#include <vector>
#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_NoCopiable.h"
#include "Amesos_Utils.h"
#include "Amesos_Time.h"
#include "Amesos_Status.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Comm.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

using namespace Teuchos;

//! Amesos_Pardiso: Interface to the PARDISO package.

/*!
  \author Marzio Sala, SNL 9214

  \date Last updated on June 2005
*/

class Amesos_Pardiso: public Amesos_BaseSolver, 
                      private Amesos_Time, 
                      private Amesos_NoCopiable, 
                      private Amesos_Utils, 
                      private Amesos_Status { 

public: 

  //@{ \name Constructor methods
  //! Constructor.
  Amesos_Pardiso(const Epetra_LinearProblem& LinearProblem );

  //! Destructor.
  ~Amesos_Pardiso(void);
  //@}

  //@{ \name Mathematical functions.

  //! Performs SymbolicFactorization on the matrix A.
  int SymbolicFactorization() ;

  //! Performs NumericFactorization on the matrix A.
  int NumericFactorization() ;

  //! Solves A X = B (or A<SUP>T</SUP> X = B) 
  int Solve();
  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

  //! Get a pointer to the Problem.
  const Epetra_LinearProblem* GetProblem() const { return(Problem_); };

  //! Returns true if KLU can handle this matrix shape 
  /*! Returns true if the matrix shape is one that KLU can
    handle. KLU only works with square matrices.  
  */
  bool MatrixShapeOK() const;

  //! SetUseTranpose(true) is more efficient in Amesos_Klu
  /*! 
    If SetUseTranspose() is set to true, 
    \f$A^T X = B\f$ is computed.
  */  
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm& Comm() const {return(GetProblem()->GetOperator()->Comm());};

  //! Set parameters from the input parameters list, returns 0 if successful.
  int SetParameters(Teuchos::ParameterList &ParameterList);

  //! Prints timing information
  void PrintTiming() const;
  
  //! Prints information about the factorization and solution phases.
  void PrintStatus() const;
  
  //@}

private:  
  
  inline const Epetra_Map& Map() const
  {
    return(Matrix_->RowMatrixRowMap());
  }
  
  inline const Epetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  inline Epetra_Map& SerialMap() 
  {
    return(*(SerialMap_.get()));
  }
  
  inline Epetra_RowMatrix& SerialMatrix()
  {
    return(*(SerialMatrix_.get()));
  }

  inline Epetra_CrsMatrix& SerialCrsMatrix()
  {
    return(*(SerialCrsMatrix_.get()));
  }

  inline Epetra_Import& Importer()
  {
    return(*(Importer_.get()));
  }
  
  int ConvertToSerial();
  int ConvertToPardiso();
  int PerformSymbolicFactorization();
  int PerformNumericFactorization(); 

  RefCountPtr<Epetra_Map> SerialMap_;
  RefCountPtr<Epetra_CrsMatrix> SerialCrsMatrix_;
  RefCountPtr<Epetra_RowMatrix> SerialMatrix_;
  RefCountPtr<Epetra_Import> Importer_;

  const Epetra_Map* Map_;
  const Epetra_RowMatrix* Matrix_;

  //! If \c true, the transpose of A is used.
  bool UseTranspose_;
  //! Pointer to the linear system problem.
  const Epetra_LinearProblem* Problem_;

  // Data for PARDISO
  vector<double> aa_;
  vector<int>    ia_;
  vector<int>    ja_;

  int mtype_;
  void* pt_[64];

  int iparm_[64];
  int maxfct_, mnum_, msglvl_;
  int nrhs_;

};  // class Amesos_Pardiso  
#endif
