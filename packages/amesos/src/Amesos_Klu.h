
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

#ifndef _EPETRA_KLU_H_
#define _EPETRA_KLU_H_

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_CrsGraph.h"

extern "C" {
  // #include "amd.h"
#include "klu_btf.h"
}

//! Amesos_Klu:  An object-oriented wrapper for Klu.
/*!  Amesos_Klu will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the Klu solver library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.

  Klu execution can be tuned through a variety of parameters.
  Amesos_Klu.h allows control of these parameters through the
  following named parameters, ignoring parameters with names that it
  does not recognize.  Where possible, the parameters are common to
  all direct solvers (although some may ignore them).  However, some
  parameters, in particular tuning parameters, are unique to each
  solver.
    
*/
class Amesos_Klu: public Amesos_BaseSolver { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Klu Constructor.
  /*! Creates an Amesos_Klu instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Klu(const Epetra_LinearProblem& LinearProblem, const AMESOS::Parameter::List &ParameterList );

  //! Amesos_Klu Destructor.
  /*! Completely deletes an Amesos_Klu object.  
  */
  ~Amesos_Klu(void);
  //@}

  //@{ \name Mathematical functions.

    //! Performs SymbolicFactorization on the matrix A.
    /*! 
      In addition to performing symbolic factorization on the matrix A, 
      the call to SymbolicFactorization() implies that no change will
      be made to the non-zero structure of the underlying matrix without 
      a subsequent call to SymbolicFactorization().
      
      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      </ul>

      postconditions:<ul>
      <li>Symbolic Factorization will be performed (or marked to be performed) 
      allowing NumericFactorization() and Solve() to be called.
      </ul>

    \return Integer error code, set to 0 if successful.
  */
    int SymbolicFactorization() ;

    //! Performs NumericFactorization on the matrix A.
    /*!  In addition to performing numeric factorization (and symbolic
      factorization if necessary) on the matrix A, the call to
      NumericFactorization() implies that no change will be made to
      the underlying matrix without a subsequent call to
      NumericFactorization().  

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().  
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      </ul>

      postconditions:<ul>
      <li>Numeric Factorization will be performed (or marked to be performed) 
      allowing Solve() to be performed correctly despite a potential change in 
      in the matrix values (though not in the non-zero structure).
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int NumericFactorization() ;

    //! Solves A X = B (or A<SUP>T</SUP> x = B) 
    /*! 

      preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>GetProblem()->CheckInput (see Epetra_LinearProblem::CheckInput() for return values)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      <li>The matrix should not have changed
          since the last call to NumericFactorization().
      </ul>

      postconditions:<ul> 
      <li>X will be set such that A X = B (or
      A<SUP>T</SUP> X = B), within the limits of the accuracy of the
      underlying solver.  
      </ul>

     \return Integer error code, set to 0 if successful.
  */
    int Solve();

  //@}
  
  //@{ \name Additional methods required to support the Epetra_Operator interface.

#if 0
  //! Returns a character string describing the operator
  char * Label() const {return(Epetra_Object::Label());};
#endif
    
  //! Get a pointer to the Problem.
  const Epetra_LinearProblem *GetProblem() const { return(Problem_); };

  //! Get a pointer to the ParameterList.
  const AMESOS::Parameter::List *GetParameterList() const { return(ParameterList_); };

  //! Returns true if KLU can handle this matrix shape 
  /*! Returns true if the matrix shape is one that KLU can
    handle. KLU only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};
  //@}

 private:  
  /*
  ConvertToSerial - Convert matrix to a serial Epetra_CrsMatrix
    Preconditions:
      Problem_ must be set 
      SerialMap and SerialCrsMatrix must either be 0 or be pointers to 
        appropriatly allocate objects.  If they are non-zero, those objects
	will be deleted (and possibly recreated).  
	
    Postconditions:
      IsLocal is set to 1 if the input matrix is entirely stored on process 0
      SerialMap points to a serial map if IsLocal==1
      SerialCrsMatrix contains a serial version of the matrix A if IsLocal==1
      SerialMatrix points to a serial copy of the matrix
      NumGlobalElements_   is set to the number of rows in the matrix
      numentries_ is set to the number of non-zeroes in the matrix 
   */
  int ConvertToSerial();

  /*
    ConvertToKluCRS - Convert matirx to form expected by Klu: Ai, Ap, Aval
    Preconditions:
      numentries_, NumGloalElements_ and SerialMatrix_ must be set.
    Postconditions:
      Ai, Ap, and Aval are resized and populated with a compresses row storage 
      version of the input matrix A.
  */
  int ConvertToKluCRS();     

  /*
    PerformSymbolicFactorization - Call Klu to perform symbolic factorization
    Preconditions:
      IsLocal must be set to 1 if the input matrix is entirely stored on process 0
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
      IsLocal must be set 
      Ap, Ai and Aval are a compressed row storage version of the input matrix A.
      Symbolic must be set
    Postconditions:
      Numeric points to an KLU internal opaque object containing the
        numeric factorization and accompanying information.  
      NumericFactorizationOK_ = true; 
    Note:  All action is performed on process 0
  */
  int PerformNumericFactorization(); 

 protected:

    int *Lp, *Li, *Up, *Ui, *P ;	
    double *Lx, *Ux ;
    klu_symbolic *Symbolic_ ;
    klu_numeric *Numeric_ ;

    

  //
  //  Ap, Ai, Aval form the compressed row storage used by Klu
  //
  vector <int> Ap;
  vector <int> Ai;
  vector <double> Aval;

  int iam;                 //  Process number (i.e. Comm().MyPID() 
  
  int IsLocal_;            //  1 if Problem_->GetOperator() is stored entirely on process 0
                           //  Note:  Local Problems do not require redistribution of
                           //  the matrix A or vectors X and B.
  int numentries_;         //  Number of non-zero entries in Problem_->GetOperator()
  int NumGlobalElements_;  //  Number of rows and columns in the Problem_->GetOperator()

  Epetra_Map *SerialMap_ ;               //  Points to a Serial Map (unused if IsLocal == 1 ) 
  Epetra_CrsMatrix *SerialCrsMatrixA_ ;  //  Points to a Serial Copy of A (unused if IsLocal==1)
  Epetra_CrsMatrix *SerialMatrix_ ;      //  Points to a Serial Copy of A 
                                         //  IsLocal==1 - Points to the original matix 
                                         //  IsLocal==0 - Points to SerialCrsMatrixA
  Epetra_CrsMatrix *TransposeMatrix_ ;   //  Points to a Serial Transposed Copy of A 
  Epetra_CrsMatrix *Matrix_ ;            //  Points to the matrix that is used to compute
                                         //  the values passed to Klu
                                     

  bool UseTranspose_;      //  Set by 
  const Epetra_LinearProblem * Problem_;
  const AMESOS::Parameter::List * ParameterList_ ; 

};  // End of  class Amesos_Klu  
#endif /* _EPETRA_KLU_H_ */
