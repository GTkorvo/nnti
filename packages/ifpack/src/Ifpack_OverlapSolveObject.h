/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * October 20, 2002, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef IFPACK_SOLVEOBJECT_H
#define IFPACK_SOLVEOBJECT_H

#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_CombineMode.h"
class Epetra_Flops;

//! Ifpack_OverlapSolveObject: Provides Overlapped Forward/back solve services for Ifpack.

class Ifpack_OverlapSolveObject: public virtual Epetra_Operator {

 public:
  //@{ \name Constructors/Destructor
  //! Constructor.
  Ifpack_OverlapSolveObject(char * Label, bool IsOverlapped, const Epetra_Comm & Comm);

  //! Copy constructor.
  Ifpack_OverlapSolveObject(const Ifpack_OverlapSolveObject & Source);

  //! Ifpack_OverlapSolveObject Destructor
  virtual ~Ifpack_OverlapSolveObject();
  //@}

  //@{ \name Initialization methods.

  //! Generate Ifpack_OverlapGraph object using current settings.
  void SetOverlapMode( Epetra_CombineMode OverlapMode) {OverlapMode_ = OverlapMode; return;}

  //! Define the operator to be used for the lower triangle
  int SetLowerOperator (Epetra_CrsMatrix * L, bool UseLtrans) {L_ = L; return(0);};

  //! Define the vector to be used for the diagonal
  int SetDiagonal (Epetra_Vector * D) {D_ = D; return(0);};

  //! Define the operator to be used for the upper triangle
  int SetUpperOperator (Epetra_CrsMatrix * U, bool UseUtrans) {U_ = U; return(0);};
  //@}

  //@{ \name Mathematical functions.
  
  
  //! Returns the result of a Ifpack_CrsIlut forward/back solve on a Epetra_MultiVector X in Y (works for Epetra_Vectors also).
  /*! 
    \param In
    Trans -If true, solve transpose problem.
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of multiplying U, D and L in that order on an Epetra_MultiVector X in Y.
  /*! 
    \param In
    Trans -If true, multiply by L^T, D and U^T in that order.
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectorscontaining result.
    
    \return Integer error code, set to 0 if successful.
  */
  int Multiply(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the maximum over all the condition number estimate for each local ILU set of factors.
  /*! This functions computes a local condition number estimate on each processor and return the
      maximum over all processor of the estimate.
   \param In
    Trans -If true, solve transpose problem.
    \param Out
    ConditionNumberEstimate - The maximum across all processors of 
    the infinity-norm estimate of the condition number of the inverse of LDU.
  */
  int Condest(bool Trans, double & ConditionNumberEstimate) const;
  //@}

  //@{ \name Attribute access methods.

  //! Returns the overlap mode used to combine terms that are redundantly computed.
  /*! Since rows of the graph, and any related matrices are multiply owned, some values
      in the subdomain solves will be computed on multiple processors.  The overlap mode
      is used to determine how the redundant values that come in from other processors will
      be handled.
  */
  Epetra_CombineMode OverlapMode() const {return(OverlapMode_);}
  
  //! Returns the number of nonzero entries in the global graph.
  int NumGlobalNonzeros() const {return(L().NumGlobalNonzeros()+U().NumGlobalNonzeros());};
  
  //! Returns the number of nonzero entries in the local graph.
  int NumMyNonzeros() const {return(L().NumMyNonzeros()+U().NumMyNonzeros());};
  
  //! Returns the address of the L factor associated with this factored matrix.
  const Epetra_CrsMatrix & L() const {return(*L_);};
    
  //! Returns the address of the D factor associated with this factored matrix.
  const Epetra_Vector & D() const {return(*D_);};
    
  //! Returns the address of the L factor associated with this factored matrix.
  const Epetra_CrsMatrix & U() const {return(*U_);};
  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

    //! Returns a character string describing the operator
    char * Label() const {return(Label_);};
    
    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param In
	   UseTranspose -If true, multiply by the transpose of operator, otherwise just use operator.

    \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

    //! Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*! Note that this implementation of Apply does NOT perform a forward back solve with
        the LDU factorization.  Instead it applies these operators via multiplication with 
	U, D and L respectively.  The ApplyInverse() method performs a solve.

    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Multiply(Ifpack_OverlapSolveObject::UseTranspose(), X, Y));};

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*! In this implementation, we use several existing attributes to determine how virtual
        method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(), 
	the Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
	if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
	be accounted for when doing a triangular solve.

    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(Solve(Ifpack_OverlapSolveObject::UseTranspose(), X, Y));};

    //! Returns 0.0 because this class cannot compute Inf-norm.
    double NormInf() const {return(0.0);};

    //! Returns false because this class cannot compute an Inf-norm.
    bool HasNormInf() const {return(false);};

    //! Returns the current UseTranspose setting.
    bool UseTranspose() const {return(UseTranspose_);};

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map & OperatorDomainMap() const {return(U_->OperatorDomainMap());};

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map & OperatorRangeMap() const{return(L_->OperatorRangeMap());};

    //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
    const Epetra_Comm & Comm() const{return(Comm_);};
  //@}

 protected:
    virtual int SetupXY(bool Trans, 
			const Epetra_MultiVector& Xin, const Epetra_MultiVector& Yin,
			Epetra_MultiVector * & Xout, Epetra_MultiVector * & Yout) const=0;
 private:
  char * Label_;
  Epetra_CrsMatrix * L_;
  bool UseLTrans_;
  Epetra_Vector * D_;
  Epetra_CrsMatrix * U_;
  bool UseUTrans_;
  bool UseTranspose_;
  const Epetra_Comm & Comm_;
  mutable double Condest_;
  Epetra_Flops * Counter_;
  bool IsOverlapped_;
  Epetra_CombineMode OverlapMode_;

};
#endif // IFPACK_SOLVEOBJECT_H
