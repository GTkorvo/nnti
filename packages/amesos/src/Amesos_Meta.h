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

#ifndef _AMESOS_KLU_H_
#define _AMESOS_KLU_H_

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_CrsGraph.h"


//! Amesos_Klu:  A serial, unblocked code ideal for getting started and for very sparse matrices, such as circuit matrces.  AmesosKlu computes <p class="code">A<sup>T</sup> X = B</p> more efficiently than <p class="code">A X = B</p>.
/*!  Amesos_Klu, an object-oriented wrapper for Klu, will solve a linear systems of equations: <TT>A X = B</TT>
   using Epetra objects and the Klu solver library, where
  <TT>A</TT> is an Epetra_RowMatrix and <TT>X</TT> and <TT>B</TT> are 
  Epetra_MultiVector objects.

<br /><br /><p>AmesosKlu computes <p class="code">A<sup>T</sup> X = B</p> 
more efficiently than <p class="code">A X = B</p>.  The
latter requires a matrix transpose - which costs both time and space.

<br /><br /><p>klu is Davis' implementation of Gilbert-Peierl's left-looking
sparse partial pivoting algorithm, with Eisenstat & Liu's symmetric
pruning.  Gilbert's version appears as [L,U,P]=lu(A) in MATLAB.
It doesn't exploit dense matrix kernels, but it is the only sparse
LU factorization algorithm known to be asymptotically optimal,
in the sense that it takes time proportional to the number of
floating-point operations.  It is the precursor to SuperLU,
thus the name ("clark Kent LU").  For very sparse matrices that
do not suffer much fill-in (such as most circuit matrices when
permuted properly) dense matrix kernels do not help, and the
asymptotic run-time is of practical importance.

<br /><br /><p>The klu_btf code first permutes the matrix to upper block
triangular form (using two algorithms by Duff and Reid,
MC13 and MC21, in the ACM Collected Algorithms).  It then permutes
each block via a symmetric minimum degree ordering (AMD, by Amestoy,
Davis, and Duff).  This ordering phase can be done just once
for a sequence of matrices.  Next, it factorizes each reordered
block via the klu routine, which also attempts to preserve
diagonal pivoting, but allows for partial pivoting if the diagonal
is to small.    
<!-- A fast-factorization version of klu is
available, which does not do any partial pivoting at all.  -->


<br /><br /><p>  Klu execution can be tuned through a variety of parameters.
  Amesos_Klu.h allows control of these parameters through the
  following named parameters, ignoring parameters with names that it
  does not recognize.  Where possible, the parameters are common to
  all direct solvers (although some may ignore them).  However, some
  parameters, in particular tuning parameters, are unique to each
  solver.
    
*/

// Amesos_Klu_Pimpl contains a pointer to two sructures defined in 
// klu.h:  klu_symbolic and klu_numeric.  This prevents Amesos_Klu.h 
// from having to include klu.h.
//
//  Doxygen does not handle forward class references well.
#ifndef DOXYGEN_SHOULD_SKIP_THIS
class Amesos_Klu_Pimpl ; 
#endif

class Amesos_Klu: public Amesos_BaseSolver { 

public: 

  //@{ \name Constructor methods
  //! Amesos_Klu Constructor.
  /*! Creates an Amesos_Klu instance, using an Epetra_LinearProblem,
      passing in an already-defined Epetra_LinearProblem object. 

      Note: The operator in LinearProblem must be an
      Epetra_RowMatrix.

  */
  Amesos_Klu(const Epetra_LinearProblem& LinearProblem );

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

    //! Solves A X = B (or A<SUP>T</SUP> X = B) 
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

  //! Returns true if KLU can handle this matrix shape 
  /*! Returns true if the matrix shape is one that KLU can
    handle. KLU only works with square matrices.  
  */
  bool MatrixShapeOK() const ;

  //! SetUseTranpose(true) is more efficient in Amesos_Klu
  /*! 
<ul>
  <li>If SetUseTranspose() is set to true, 
    <ul>
       <li><p class="code">A<sup>T</sup> X = B</p> is computed</li>
       <li>(This is the more efficient operation)</li>
    </ul></li>
  <li>else
    <ul>
       <li><p class="code">A X = B</p> is computed</li>
       <li>(This requires a matrix transpose)</li>
    </ul></li>
</ul>
  */  
  int SetUseTranspose(bool UseTranspose) {UseTranspose_ = UseTranspose; return(0);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns a pointer to the Epetra_Comm communicator associated with this matrix.
  const Epetra_Comm & Comm() const {return(GetProblem()->GetOperator()->Comm());};

  //!  Updates internal variables. 
  /*!  
      <br \>Preconditions:<ul>
      <li>None.</li>
      </ul>

      <br \>Postconditions:<ul> 
      <li>Internal variables controlling the factorization and solve will
      be updated and take effect on all subsequent calls to NumericFactorization() 
      and Solve().</li>
      <li>All parameters whose value are to differ from the default values must 
be included in ParameterList.  Parameters not specified in ParameterList 
revert to their default values.
      </ul>

    \return Integer error code, set to 0 if successful. 
   */
  int SetParameters( Teuchos::ParameterList &ParameterList )  ;

  //! Print timing information
  void PrintTiming();
  
  //! Print information about the factorization and solution phases.
  void PrintStatus();
  
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
  int ConvertToKluCRS(bool firsttime);     

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
    Amesos_Klu_Pimpl *PrivateKluData_; 
    

  //
  //  Ap, Ai, Aval form the compressed row storage used by Klu
  //
  vector <int> Ap;
  vector <int> Ai;
  vector <double> Aval;

  int iam;                 //  Process number (i.e. Comm().MyPID() 
  int NumProcs_;
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
  //
  //  This USE_VIEW define shows that we can could easily support the 
  //  Epetra_RowMatrix interface.  
  //
#define USE_VIEW
#ifdef USE_VIEW
  Epetra_CrsMatrix *Matrix_ ;            //  Points to the matrix that is used to compute
#else
  Epetra_RowMatrix *Matrix_ ;            //  Points to the matrix that is used to compute
#endif                                         //  the values passed to Klu
                                     

  bool UseTranspose_;
  const Epetra_LinearProblem * Problem_;

  bool PrintTiming_;
  bool PrintStatus_;
  bool ComputeVectorNorms_;
  bool ComputeTrueResidual_;
  
  int verbose_;
  int debug_;

  // some timing internal, copied from MUMPS
  double ConTime_;                        // time to convert to KLU format
  double SymTime_;                        // time for symbolic factorization
  double NumTime_;                        // time for numeric factorization
  double SolTime_;                        // time for solution
  double VecTime_;                        // time to redistribute vectors
  double MatTime_;                        // time to redistribute matrix
  
  int NumSymbolicFact_;
  int NumNumericFact_;
  int NumSolve_;  

  Epetra_Time * Time_;
  
};  // End of  class Amesos_Klu  
#endif /* _AMESOS_KLU_H_ */
