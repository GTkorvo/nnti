
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
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

#ifndef _AZTECOO_H_
#define _AZTECOO_H_

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "Epetra_Vector.h"
#include "az_aztec.h"


//! AztecOO:  An object-oriented wrapper for Aztec.
/*! Currently it accepts a Petra matrix, initial guess and RHS as
  separate arguments, or alternatively, accepts a Epetra_LinearProblem.
  If constructed using a Epetra_LinearProblem, AztecOO will infer some
  solver/preconditioner, etc., options and parameters. Users may override
  these choices and manually choose from among the full set of Aztec options
  using the SetAztecOption() and SetAztecParam() functions.

  AztecOO will solve a linear systems of equations: \f$ AX=B \f$, using Epetra
  objects and the Aztec solver library, where \f$A\f$ is an Epetra_Operator or Epetra_RowMatrix (note
  that the Epetra_Operator class is a base class for Epetra_RowMatrix so that Epetra_RowMatrix \e isa
  Epetra_Operator.) \f$X\f$ and \f$B\f$ are Epetra_MultiVector objects.

  \warning AztecOO does not presently support solution of more than one simultaneous right-hand-side.
*/

class AztecOO {
    
  public:
  //@{ \name Constructors/destructors.
  //!  AztecOO Constructor.
  /*! Creates a AztecOO instance, passing in already-defined objects for the linear operator
      (as an Epetra_Operator),
      left-hand-side and right-hand-side. 

      Note: Use of this constructor may prohibit use of native AztecOO preconditioners, since
      an Epetra_Operator is not necessarily an Epetra_RowMatrix and all AztecOO incomplete
      factorization preconditioners are based on having explicit access to matrix coefficients.
      Polynomial preconditioners are available if the Epetra_Operator passed in here has a
      non-trivial definition of the NormInf() method and HasNormInf() returns true.
  */
  AztecOO(Epetra_Operator * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
  //!  AztecOO Constructor.
  /*! Creates a AztecOO instance, passing in already-defined objects for the linear operator
      (as an Epetra_RowMatrix),
      left-hand-side and right-hand-side. 

      Note: Use of this constructor allows full access to native AztecOO preconditioners, using
      the Epetra_RowMatrix A passed in here as the basis for computing the preconditioner.  
      All AztecOO incomplete
      factorization preconditioners are based on having explicit access to matrix coefficients.
      Polynomial preconditioners are also available.  It is possible to change the matrix used for
      computing incomplete factorization by calling the SetPrecMatrix() method.  It is 
      also possible to provide a user-supplied preconditioner by call SetPrecOperator().
  */
  AztecOO(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B);

  //! AztecOO Constructor.
  /*! Creates a AztecOO instance, using a Epetra_LinearProblem, 
      passing in an already-defined Epetra_LinearProblem object. The Epetra_LinearProblem class
      is the preferred method for passing in the linear problem to AztecOO because this class
      provides scaling capabilities and self-consistency checks that are not available when
      using other constructors.

      Note: If the Epetra_LinearProblem passed in here has a non-trivial pointer to an Epetra_Matrix
      then use of this constructor allows full access to native AztecOO preconditioners, using
      the Epetra_RowMatrix A passed in here as the basis for computing the preconditioner.  
      All AztecOO incomplete
      factorization preconditioners are based on having explicit access to matrix coefficients.
      Polynomial preconditioners are also available.  It is possible to change the matrix used for
      computing incomplete factorization by calling the SetPrecMatrix() method.  It is 
      also possible to provide a user-supplied preconditioner by call SetPrecOperator().

      If the Epetra_LinearProblems passed in here has only an Epetra_Operator, then use 
      of this constructor may prohibit use of native AztecOO preconditioners, since
      an Epetra_Operator is not necessarily an Epetra_RowMatrix and all AztecOO incomplete
      factorization preconditioners are based on having explicit access to matrix coefficients.
      Polynomial preconditioners are available if the Epetra_Operator passed in here has a
      non-trivial definition of the NormInf() method and HasNormInf() returns true.
  */
  AztecOO(const Epetra_LinearProblem& LinearProblem);

  //! AztecOO Default constructor.
  AztecOO();

  //! AztecOO Copy Constructor.
  /*! Makes copy of an existing AztecOO instance.
  */
  AztecOO(const AztecOO& Solver);

  //! AztecOO Destructor.
  /*! Completely deletes a AztecOO object.  
  */
  virtual ~AztecOO(void);
  //@}
  
  //@{ \name Post-construction setup methods.

  //! AztecOO Epetra_LinearProblem Set
  /*! Associates an already defined Epetra_LinearProblem as the problem that will be solved during
      iterations.  This method allows the user to change which problem is being solved by an existing
      AztecOO object.
      \warning If a preconditioner has been pre-built and associated with this AztecOO object, the 
      Epetra_LinearProblem being passed in to this method \e must have compatible domain and range maps.
   */
  int SetProblem(const Epetra_LinearProblem& prob);

  //! AztecOO User Operator Set
  /*! Associates an already defined Epetra_Operator as the linear operator for the linear system
      system that will be solved during
      iterations.
      This set method allows the user to pass any type of linear operator to AztecOO, as long
      as the operator implements the Epetra_Operator pure virtual class, and has proper
      domain and range map dimensions. Epetra_CrsMatrix and Epetra_VbrMatrix objects can be passed in through 
      this method.
   */
  int SetUserOperator(Epetra_Operator * UserOperator);

  //! AztecOO User Matrix Set
  /*! Associates an already defined Epetra_Matrix as the matrix that will be used by
      AztecOO as the linear operator when solving the linear system.  
      Epetra_CrsMatrix and Epetra_VbrMatrix objects can be passed in through 
      this method.  This method also sets the preconditioner matrix to the matrix passed in here.
   */
  int SetUserMatrix(Epetra_RowMatrix * UserMatrix);

  //! AztecOO LHS Set
  /*! Associates an already defined Epetra_MultiVector (or Epetra_Vector) as the initial guess
      and location where the solution will be return.
   */
  int SetLHS(Epetra_MultiVector * X);


  //! AztecOO RHS Set
  /*! Associates an already defined Epetra_MultiVector (or Epetra_Vector) as the right-hand-side of 
      the linear system.
   */
  int SetRHS(Epetra_MultiVector * B);

  //! AztecOO Preconditioner Matrix Set
  /*! Associates an already defined Epetra_Matrix as the matrix that will be used by
      AztecOO when constructing native AztecOO preconditioners.  By default,
      if AztecOO native preconditioners are used, the original operator matrix will be used as
      the source for deriving the preconditioner.  However, there are instances where a user would like
      to have the preconditioner be defined using a different matrix than the original operator matrix.
      Another common situation is where the user may not have the operator in matrix form but has a matrix
      that approximates the operator and can be used as the basis for an incomplete factorization.
      This set method allows the user to pass any Epetra_RowMatrix to AztecOO for use in constructing an AztecOO
      native preconditioner, as long
      as the matrix implements the Epetra_RowMatrix pure virtual class, and has proper
      domain and range map dimensions.  Epetra_CrsMatrix and Epetra_VbrMatrix objects can be passed in through 
      this method.
   */
  int SetPrecMatrix(Epetra_RowMatrix * PrecMatrix);

  //! AztecOO External Preconditioner Set
  /*! Associates an already defined Epetra_Operator as the preconditioner that will be called during
      iterations.
      This set method allows the user to pass any type of preconditioner to AztecOO, as long
      as the preconditioner implements the Epetra_Operator pure virtual class, and has proper
      domain and range map dimensions.  Ifpack preconditioners can be passed in through this method.
   */
  int SetPrecOperator(Epetra_Operator * PrecOperator);
  //@}
  
  //@{ \name Post-construction setup methods (classic approach: ONLY EXPERTS SHOULD USE THESE METHODS).
  //! AztecOO External Preconditioner Set (object)
  /*! Associates an already defined Aztec preconditioner with this solve.
   */
  int SetPreconditioner(AZ_PRECOND * Prec) {Prec_ = Prec; return(0);};


  //! AztecOO External Preconditioner Set (function and data)
  /*! Associates an external function and data pointer with preconditioner
   */
  int SetPreconditioner(void  (*prec_function)(double *, int *, int *, double *,
					   struct AZ_MATRIX_STRUCT  *, 
					   struct AZ_PREC_STRUCT *),
		    void *prec_data);

  //! AztecOO External Scaling Set
  /*! Associates an already defined Aztec scaling object with this solve.
   */
  int SetScaling(struct AZ_SCALING * Scaling) {Scaling_ = Scaling; return(0);};


  //! AztecOO Label Matrix for Aztec
  /*! This is used to label individual matrices within Aztec. This might
      be useful if several Aztec invokations are involved corresponding
      to different matrices.
   */
  int  SetMatrixName(int label);
  //@}
  
  //@{ \name Explicit preconditioner construction/assessment/destruction methods.
  //! Forces explicit construction and retention of an AztecOO native preconditioner.
  /*! AztecOO typically constructs the preconditioner on the first call to the solve function.
      However, there are situations where we would like to compute the preconditioner ahead
      of time.  One particular case is when we want to confirm that the preconditioner 
      well-conditioned.  This method allows us to precompute the preconditioner.  It also
      provides a estimate of the condition number of the preconditioner.  If \it condest is
      large, e.g., > 1.0e+14, it is likely the preconditioner will fail.  In this case, using 
      threshold values (available in the incomplete factorizations) can be used to reduce
      the condition number.

      Note: This method does not work for user-defined preconditioners (defined via calls to 
      SetPrecOperator().  It will return with an error code  of -1 for this case.
  */
  int ConstructPreconditioner(double & condest);

  //! Destroys a preconditioner computed using ConstructPreconditioner().
  /*! The ConstructPreconditioner() method creates a persistent preconditioner.
      In other words the preconditioner will be used by all calls to the Iterate() 
      method.  DestroyPreconditioner() deletes the current preconditioner and restores
      AztecOO to a state where the preconditioner will computed on first use of the
      preconditioner solve.
  */
  int DestroyPreconditioner();

  //! Returns the condition number estimate for the current, if one exists, returns -1.0 if no estimate
  double Condest() const {return(condest_);};
  //@}

  //@{ \name Check/Attribute Access Methods.
    
  //! Prints a summary of solver parameters, performs simple sanity checks.
  int CheckInput() const {
    return(AZ_check_input(Amat_->data_org, options_, params_, proc_config_));};

  //! Get a pointer to the user operator A.
  Epetra_Operator * GetUserOperator() const {return(UserOperatorData_->A);};
  //! Get a pointer to the user matrix A.
  Epetra_RowMatrix * GetUserMatrix() const {return(UserMatrixData_->A);};
  //! Get a pointer to the preconditioner operator.
  Epetra_Operator * GetPrecOperator() const {return(PrecOperatorData_->A);};
  //! Get a pointer to the matrix used to construct the preconditioner.
  Epetra_RowMatrix * GetPrecMatrix() const {return(PrecMatrixData_->A);};
  //! Get a pointer to the left-hand-side X.
  Epetra_MultiVector * GetLHS() const {return(X_);};
  //! Get a pointer to the right-hand-side B.
  Epetra_MultiVector * GetRHS() const {return(B_);};
  //@}

  //@{ \name Standard AztecOO option and parameter setting methods.
 
  //! AztecOO function to restore default options/parameter settings.
  /*! This function is called automatically within AztecOO's constructor,
   but if constructed using a Epetra_LinearProblem object, some options are
   reset based on the ProblemDifficultyLevel associated with the 
   Epetra_LinearProblem. 
   
   See the Aztec 2.1 User Guide for a complete list of these options.

  */
    int SetAztecDefaults();

  //! AztecOO option setting function.
  /*! Set a specific Aztec parameter value.
      Example: problem.SetAztecOption(AZ_precond, AZ_Jacobi)
      
      See the Aztec 2.1 User Guide for a complete list of these options.
   */
  int SetAztecOption(int option, int value)
    {options_[option] = value; return(0);};

  //! AztecOO param setting function.
  /*! Set a specific Aztec parameter value.
      Example: problem.SetAztecParam(AZ_drop, 1.0E-6)
      
      See the Aztec 2.1 User Guide for a complete list of these parameters.
   */
    int SetAztecParam(int param, double value)
    {params_[param] = value; return(0);};

  //! AztecOO option setting function.
  /*! Set all Aztec option values using an existing Aztec options array.
   */
    int SetAllAztecOptions(int * options)
      {for (int i=0; i<AZ_OPTIONS_SIZE; i++) options_[i] = options[i]; return(0);};

  //! AztecOO param setting function.
  /*! Set all Aztec parameter values using an existing Aztec params array. 
   */
    int SetAllAztecParams(double * params)
      {for (int i=0; i<AZ_PARAMS_SIZE; i++) params_[i] = params[i]; return(0);};
  //@}

  //@{ \name Standard AztecOO solve methods.
  //! AztecOO iteration function.
  /*! Iterates on the current problem until MaxIters or Tolerance is reached.
  */
  int Iterate(int MaxIters, double Tolerance);

  //! AztecOO iteration function.
  /*! Iterates on the specified matrix and vectors until MaxIters or Tolerance
      is reached..
  */
  int Iterate(Epetra_RowMatrix * A,
                      Epetra_MultiVector * X,
                      Epetra_MultiVector * B,
                      int MaxIters, double Tolerance);


  //@}

  //@{ \name Specialist AztecOO solve method.
  //! AztecOO iteration functions.
  /*! Iterates on the current problem until MaxIters or Tolerance is reached..
      This one should be suitable for recursive invocations of Aztec.
  */
  int recursiveIterate(int MaxIters, double Tolerance);

  //@}

  //@{ \name Adaptive Solve methods.

    //! Force the AdaptiveIterate() method to use default adaptive strategy.
    int SetUseAdaptiveDefaultsTrue(){useAdaptiveDefaults_ = true;return(0);};

  //! Set the parameter that control the AdaptiveIterate() method.
  /*! The AdaptiveIterate() method attempts to solve a given problem using multiple preconditioner
      and iterative method tuning parameters.  There are defaults that are coded into AdaptiveIterate()
      method, but the defaults can be over-ridden by the use of the SetAdaptiveParams() method. Details of 
      condition number management follow:
\verbinclude Managing_conditioning_howto.txt

    \param NumTrials In
           The number of Athresh and Rthresh pairs that should be tried when attempting to stabilize 
	   the preconditioner.
    \param athresholds In
           The list of absolute threshold values that should be tried when attempting to stabilize 
	   the preconditioner.
    \param rthresholds In
           The list of relative threshold values that should be tried when attempting to stabilize 
	   the preconditioner.
    \param condestThreshold In
           If the condition number estimate of the preconditioner is above this number, no attempt will be
	   made to try iterations.  Instead a new preconditioner will be computed using the next threshold 
	   pair.
    \param maxFill In
           In addition to managing the condest, the AdaptiveIterate() method will also try to increase
	   the preconditioner fill if it is determined that this might help.  maxFill specifies the
	   maximum fill allowed.
    \param maxKspace In
           In addition to managing the condest, the AdaptiveIterate() method will also try to increase
	   the Krylov subspace size if GMRES is being used and it is determined that this might help.
	   maxKspace specifies the maximum Krylov subspace allowed.
	   
  */
    int SetAdaptiveParams(int NumTrials, double * athresholds, double * rthresholds,
			  double condestThreshold, double maxFill, int maxKspace);

    //! Attempts to solve the given linear problem using an adaptive strategy.
    int AdaptiveIterate(int MaxIters, int MaxSolveAttempts, double Tolerance);
  //@}

  //@{ \name Post-solve access functions

    //! Returns the total number of iterations performed on this problem.
    int NumIters() const {return((int) status_[AZ_its]);};
    
    //! Returns the true unscaled residual for this problem.
    double TrueResidual() const {return(status_[AZ_r]);};
    
    //! Returns the true scaled residual for this problem.
    double ScaledResidual() const {return(status_[AZ_scaled_r]);};
    
    //! AztecOO status extraction function.
    /*! Extract Aztec status array into user-provided array.  The array must be of 
        length AZ_STATUS_SIZE as defined in the az_aztec.h header file.
     */
    int GetAllAztecStatus(double * status)
      {for (int i=0; i<AZ_STATUS_SIZE; i++) status[i] = status_[i]; return(0);};
  //@}


	struct MatrixData {
		Epetra_RowMatrix * A;
		Epetra_Vector * X;
		Epetra_Vector * Y;
		Epetra_Vector * SourceVec;
		Epetra_Vector * TargetVec;

		MatrixData(Epetra_RowMatrix * inA = 0, Epetra_Vector * inX = 0, 
							 Epetra_Vector * inY = 0, Epetra_Vector * inSourceVec = 0,
							 Epetra_Vector * inTargetVec = 0)
			: A(inA), X(inX), Y(inY), SourceVec(inSourceVec), TargetVec(inTargetVec){};

		~MatrixData(void) {
			if (X!=0) delete X;
			if (Y!=0) delete Y;
			if (SourceVec!=0) delete SourceVec;
			if (TargetVec!=0) delete TargetVec;
		};
	};

	struct OperatorData {
		Epetra_Operator * A;
		Epetra_Vector * X;
		Epetra_Vector * Y;
		OperatorData(Epetra_Operator * inA = 0, Epetra_Vector * inX = 0, 
							 Epetra_Vector * inY = 0)
			: A(inA), X(inX), Y(inY) {};
		~OperatorData(void) {
			if (X!=0) delete X;
			if (Y!=0) delete Y;
		};
	};

 protected:

  int AllocAzArrays();
  void DeleteAzArrays();
  int SetAztecVariables();
  int SetProblemOptions(ProblemDifficultyLevel PDL,
                            bool ProblemSymmetric);
  int SetProcConfig(const Epetra_Comm & Comm);

  void DeleteMemory();
  

  Epetra_MultiVector * X_;
  Epetra_MultiVector * B_;

  int N_local_;
  int x_LDA_;
  double *x_;
  int b_LDA_;
  double *b_;
  int * proc_config_;
  int * options_;
  double * params_;
  double * status_;
  AZ_MATRIX *Amat_;
  AZ_MATRIX *Pmat_;
  AZ_PRECOND *Prec_;
  struct AZ_SCALING * Scaling_;

  double condest_;
  bool useAdaptiveDefaults_;
  int NumTrials_;
  double maxFill_;
  int maxKspace_;
  double * athresholds_;
  double * rthresholds_;
  double condestThreshold_;
  bool inConstructor_; // Shuts down zero pointer error reporting while in a constructor
  bool procConfigSet_;

	MatrixData * UserMatrixData_;
	MatrixData * PrecMatrixData_;
	OperatorData * UserOperatorData_;
	OperatorData * PrecOperatorData_;

};

// External prototypes
extern "C" void Epetra_Aztec_matvec(double x[], double y[], AZ_MATRIX *Amat, int proc_config[]);
extern "C" void Epetra_Aztec_operatorvec(double x[], double y[], AZ_MATRIX *Amat, int proc_config[]);
extern "C" void Epetra_Aztec_precond(double x[], int input_options[],
				     int proc_config[], double input_params[], AZ_MATRIX *Amat,
				     AZ_PRECOND *prec);
extern "C" int Epetra_Aztec_getrow(int columns[], double values[], int row_lengths[], 
				   AZ_MATRIX *Amat, int N_requested_rows,
				   int requested_rows[], int allocated_space);
extern "C" int Epetra_Aztec_comm_wrapper(double vec[], AZ_MATRIX *Amat);
		
#endif /* _AZTECOO_H_ */

