/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER
*/

#ifndef IFPACK_EUCLID_H
#define IFPACK_EUCLID_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_EUCLID

#include "Ifpack_Condest.h"
#include "Ifpack_ScalingType.h"
#include "Epetra_CompObject.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"

#include "Mem_dh.h"
#include "io_dh.h"
#include "TimeLog_dh.h"
#include "Parser_dh.h"
#include "Euclid_dh.h"

namespace Teuchos {
  class ParameterList;
}

//! Ifpack_Euclid: A class for constructing and using an ILU factorization of a given Epetra_CrsMatrix, using the Euclid library by Argonne National Laboratories.


class Ifpack_Euclid: public Epetra_Object, public Epetra_CompObject, public virtual Epetra_Operator {
      
  friend ostream& operator << (ostream& os, const Ifpack_Euclid& A);

public:
  // @{ Constructors and destructors.
  //! Constructor
  Ifpack_Euclid(Epetra_CrsMatrix* A);
  
  //! Destructor
  ~Ifpack_Euclid()
  {
    Destroy();
  }

  // @}
  // @{ Construction methods
  
  //! Initialize the preconditioner, does not touch matrix values.
  int Initialize();
  
  //! Returns \c true if the preconditioner has been successfully initialized.
  bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  //! Compute ILU factors L and U using the specified graph, diagonal perturbation thresholds and relaxation parameters.
  /*! This function computes the ILU(k) factors.
   */
  int Compute();

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool IsComputed() const 
  {
    return(IsComputed_);
  }


  //! Set parameters using a Teuchos::ParameterList object.
  /*! This method is only available if the Teuchos package is enabled.
  \param ParameterList (In) - The Parameter list. Options are:
       SetLevel (int)
       SetBJ (int)
       SetStats (int)
       SetMem (int)
       SetSparse (double)
       SetRowScale (int)
       SetILUT (double)

  \return Integer error code, set to 0 if successful.
  */
  int SetParameters(Teuchos::ParameterList& parameterlist);

  //! Set a parameter that takes a single int.
  /*!
  \param name (In) -The parameter that is getting set.
  \param Value (In) -An integer value corresponding to the parameter.

  \return Integer error code, set to 0 if successful.
 */
  int SetParameter(string name, int Value);

  //! Set a parameter that takes a single double.
  /*!
  \param name (In) -The parameter that is getting set.
  \param Value (In) -A double value corresponding to the parameter.

  \return Integer error code, set to 0 if successful.
  */
  int SetParameter(string name, double Value);

  int SetUseTranspose(bool UseTranspose_in) {UseTranspose_ = UseTranspose_in; return(0);};
  // @}

  // @{ Mathematical functions.
  // Applies the matrix to X, returns the result in Y.
  int Apply(const Epetra_MultiVector& X, 
	       Epetra_MultiVector& Y) const
  {
    return(Multiply(false,X,Y));
  }

  int Multiply(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{ return A_->Multiply(Trans, X, Y); }

  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! In this implementation, we use several existing attributes to determine how virtual
      method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(), 
      the Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
      if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
      be accounted for when doing a triangular solve.

    \param 
	   X - (In) A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y - (Out) A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful. -1 if compute() hasn't been called. -2 if the multivectors have differing numbers of vectors.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Computes the estimated condition number and returns the value.
  double Condest(const Ifpack_CondestType CT = Ifpack_Cheap, 
                 const int MaxIters = 1550,
                 const double Tol = 1e-9,
		 Epetra_RowMatrix* Matrix_in = 0);

  //! Returns the computed estimated condition number, or -1.0 if not computed.
  double Condest() const
  {
    return(Condest_);
  }

  // @}
  // @{ Query methods
  
  //! Returns a character string describing the operator
  const char* Label() const {return(Label_);}

  //! Sets label for \c this object.
  void SetLabel(const char* Label_in)
  {
    strcpy(Label_,Label_in);
  }
  
  //! Returns the domain map from the creating matrix.
  const Epetra_Map &OperatorDomainMap() const{return A_->DomainMap();}

  //! Returns the range map from the creating matrix.
  const Epetra_Map &OperatorRangeMap() const{return A_->RangeMap();}

  //! Returns 0.0 because this class cannot compute Inf-norm.
  double NormInf() const {return(0.0);};

  //! Returns false because this class cannot compute an Inf-norm.
  bool HasNormInf() const {return(false);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
  const Epetra_Comm & Comm() const{return(A_->Comm());};

  //! Returns a reference to the matrix to be preconditioned.
  const Epetra_CrsMatrix& Matrix() const
  { 
    return(*A_);
  }

  //! Returns the number of calls to Initialize().
  virtual int NumInitialize() const
  {
    return(NumInitialize_);
  }

  //! Returns the number of calls to Compute().
  virtual int NumCompute() const
  {
    return(NumCompute_);
  }

  //! Returns the number of calls to ApplyInverse().
  virtual int NumApplyInverse() const
  {
    return(NumApplyInverse_);
  }

  //! Returns the time spent in Initialize().
  virtual double InitializeTime() const
  {
    return(InitializeTime_);
  }

  //! Returns the time spent in Compute().
  virtual double ComputeTime() const
  {
    return(ComputeTime_);
  }

  //! Returns the time spent in ApplyInverse().
  virtual double ApplyInverseTime() const
  {
    return(ApplyInverseTime_);
  }

  //! Returns the number of flops in the initialization phase.
  virtual double InitializeFlops() const
  {
    return(0.0);
  }

  //! Returns the number of flops in the compute phase.
  virtual double ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  //! Returns the number of flops in the applyinverse phase.
  virtual double ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }

private:

  // @}
  // @{ Private methods

  //! Copy constructor (should never be used)
  Ifpack_Euclid(const Ifpack_Euclid& RHS) :
    Time_(RHS.Comm())
  {}

  //! operator= (should never be used)
  Ifpack_Euclid& operator=(const Ifpack_Euclid& RHS)
  {
    return(*this);
  }

  //! Destroys all internal data
  void Destroy();

  //! Returns the MPI comm used in the matrix that created the preconditioner.
  MPI_Comm GetMpiComm() const{
    return (dynamic_cast<const Epetra_MpiComm*>(&A_->Comm()))->GetMpiComm();
  }

  //! Internal method to call the euclid solve method.
  int CallEuclid(double *x, double *y) const;

  //! Returns the result of a Ifpack_ILU forward/back solve on a Epetra_MultiVector X in Y.
  /*! 
    \param In
    Trans -If true, solve transpose problem.
    \param 
    X - (In) A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y - (Out) A Epetra_MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;


  //! Returns the number of global matrix rows.
  int NumGlobalRows() const {return(A_->NumGlobalRows());};
  
  //! Returns the number of global matrix columns.
  int NumGlobalCols() const {return(A_->NumGlobalCols());};
  
  //! Returns the number of local matrix rows.
  int NumMyRows() const {return(A_->NumMyRows());};
  
  //! Returns the number of local matrix columns.
  int NumMyCols() const {return(A_->NumMyCols());};
  
  // @}
  // @{ Internal data
  
  //! Pointer to the Epetra_CrsMatrix to factorize
  Teuchos::RefCountPtr<Epetra_CrsMatrix> A_;
  //! This objects copy of the ParamterList.
  Teuchos::ParameterList List_;
  //! If true, use transpose operator operations.
  bool UseTranspose_;
  //! The condition estimate for this preconditioner, will be -1 for now. 
  double Condest_;
  //! If \c true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If \c true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! Label of \c this object.
  char Label_[160];
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to ApplyInverse().
  mutable int NumApplyInverse_;
  //! Contains the time for all successful calls to Initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to ApplyInverse().
  mutable double ApplyInverseTime_;
  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable double ApplyInverseFlops_;
  //! Used for timing issues.
  mutable Epetra_Time Time_;
  //! This is the Euclid solver.
  Euclid_dh eu;
  //! true if the row map of provided matrix is in form that Euclid likes
  bool NiceRowMap_;
  //! Set livel k for ILU(k) factorization
  int SetLevel_;
  //! block-jacobi solver
  int SetBJ_;
  //! print stats
  int SetStats_;
  //! print memory usage information
  int SetMem_;
  //! define drop-tolerance
  double SetSparse_;
  //! scale values prior to factorization
  int SetRowScale_;
  //! drop tolerance relative to the absolute value of any entry in the row being factored
  double SetILUT_;
};

//! This is the print function.
ostream& operator << (ostream& os, const Ifpack_Euclid& A);

#endif // HAVE_EUCLID
#endif /* IFPACK_EUCLID_H */
