#ifndef IFPACK_AMESOS_H
#define IFPACK_AMESOS_H

#include "Ifpack_ConfigDefs.h"
#if defined(HAVE_IFPACK_AMESOS) && defined(HAVE_IFPACK_TEUCHOS)

#include "Ifpack_Preconditioner.h"
#include "Epetra_Operator.h"
#include "Teuchos_ParameterList.hpp"
class Epetra_Map;
class Epetra_Comm;
class Amesos_BaseSolver;
class Epetra_LinearProblem;
class Epetra_RowMatrix;

//! Ifpack_Amesos: a class to use Amesos's LU solvers as preconditioners.
/*!
Class Ifpack_Amesos enables the use of Amesos's LU factorizations as preconditioners.
Like most other Ifpack_Preconditioner objects, this class requires the
pointer of the precondioned matrix in the constructor. Parameters are
set using \c SetParameters(), and the preconditioner is computed using
Compute().

It is important to note that this class solves the linear system defined by
the matrix, with \e all the processes contained in the \c Comm object of this
matrix. As such, this preconditioner is seldomly used alone; moreover, it
is used as local solver for Ifpack_AdditiveSchwarz or Ifpack_CrsAdditiveSchwarz preconditioners. (In these cases, the distributed matrix is localized using
Ifpack_LocalRowMatrix, than this local matrix is given to Ifpack_Amesos).

Ifpack_Amesos is just a bare-bone wrap to Amesos. The most
important parameter
required by this class is \c "amesos: solver type" (defaulted to \c
"Amesos_Klu"), which defined the Amesos solver. However, the \e same
list is used to set Amesos parameters.

\author Marzio Sala, SNL 9214

\date First development Sep-04, last update Oct-04

*/
class Ifpack_Amesos : public Ifpack_Preconditioner {
      
public:

  Ifpack_Amesos(Epetra_RowMatrix* Matrix);

  //@{ \name Destructor.
  //! Destructor
  virtual ~Ifpack_Amesos();
  //@}

  //@{ \name Atribute set methods.

    //! If set true, transpose of this operator will be applied (not implemented).
    /*! This flag allows the transpose of the given operator to be used 
     * implicitly.  
      
    \param 
	   UseTranspose - (In) If true, multiply by the transpose of operator, 
	   otherwise just use operator.

    \return Integer error code, set to 0 if successful.  Set to -1 if this implementation does not support transpose.
  */
    virtual int SetUseTranspose(bool UseTranspose);
  //@}
  
  //@{ \name Mathematical functions.

    //! Applies the matrix to an Epetra_MultiVector.
  /*! 
    \param
    X - (In) A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param 
    Y - (Out) A Epetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
    */
    virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Applies the preconditioner to X, returns the result in Y.
  /*! 
    \param 
    X - (In) A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param 
    Y - (Out) A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method 
    must support the case where X and Y are the same object.
    */
    virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the infinity norm of the global matrix (not implemented)
    virtual double NormInf() const;
  //@}
  
  //@{ \name Atribute access functions

    //! Returns a character string describing the operator
    virtual char * Label() const;

    //! Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;

    //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const;

    //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const Epetra_Comm & Comm() const;

    //! Returns the Epetra_Map object associated with the domain of this operator.
    virtual const Epetra_Map & OperatorDomainMap() const;

    //! Returns the Epetra_Map object associated with the range of this operator.
    virtual const Epetra_Map & OperatorRangeMap() const;
  //@}

  virtual bool IsInitialized() const
  {
    return(IsInitialized_);
  }

  virtual int Initialize();

  //! Returns \c true if the preconditioner has been successfully computed.
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Computes the preconditioners.
  /*! Computes the preconditioners: 
   *
   * \return
   * 0 if successful, 1 if problems occurred.
   */
  virtual int Compute();

  //! Sets all the parameters for the preconditioner.
  /*! Only two parameters are recognized by Ifpack_Amesos:
   * - \c "amesos: solver type" : Specifies the solver type
   *   for Amesos. Default: \c Amesos_Klu.
   * - \c "amesos: compute condest" : if \c true,
   *   computes the estimated condition number. Default: \c false.
   */   
  virtual int SetParameters(Teuchos::ParameterList& List);

  //! Sets integer parameters.
  virtual int SetParameter(const string Name, const int value)
  {
    return(0);
  }

  //! Sets double parameters.
  virtual int SetParameter(const string Name, const double value)
  {
    return(0);
  }
  //! Returns a reference to the internally stored matrix.
  virtual const Epetra_RowMatrix& Matrix() const
  {
    return(*Matrix_);
  }

  //! Returns a reference to the internally stored matrix.
  virtual Epetra_RowMatrix& Matrix()
  {
    return(*Matrix_);
  }

  //! Returns the estimated condition number (if computed).
  virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
			 Epetra_RowMatrix* Matrix= 0)
  {
    // FIXME
    return(Condest_);
  }
  
  virtual double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
			 Epetra_RowMatrix* Matrix= 0) const
  {
    return(Condest_);
  }

  virtual std::ostream& Print(std::ostream& os) const;

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

  virtual long int ComputeFlops() const
  {
    return(ComputeFlops_);
  }

  virtual long int ApplyInverseFlops() const
  {
    return(ApplyInverseFlops_);
  }


private:

  //! Pointers to the matrix to be preconditioned.
  Epetra_RowMatrix* Matrix_;
  //! Amesos solver, use to apply the inverse of the local matrix.
  Amesos_BaseSolver* Solver_;
  //! Linear problem required by Solver_.
  Epetra_LinearProblem* Problem_;
  //! Contains the label of \c this object.
  string Label_;
  //! If true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! Contains the estimated condition number.
  double Condest_;
  //! If true, compute the estimated condition number.
  bool ComputeCondest_;
  //! Contains a copy of the input parameter list.
  Teuchos::ParameterList List_;
 
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to ApplyInverse().
  int NumApplyInverse_;

  //! Contains the time for all successful calls to Initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to ApplyInverse().
  double ApplyInverseTime_;

  //! Contains the number of flops for Compute().
  long int ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  long int ApplyInverseFlops_;

};

#endif // HAVE_IFPACK_AMESOS && HAVE_IFPAC_TEUCHOS
#endif // IFPACK_AMESOS_H
