
/*
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
*/

/*!
 * \file Amesos_BaseSolver.h
 *
 * \class Amesos_BaseSolver
 *
 * \brief Pure virtual class for all Amesos concrete implementions
 *
 * \date Last updated on 24-May-05.
 */

#ifndef _AMESOS_BASESOLVER_H_
#define _AMESOS_BASESOLVER_H_

#include "Teuchos_ParameterList.hpp"
#include "Epetra_LinearProblem.h"
class Epetra_LinearProblem;
class Epetra_MultiVector;
class Epetra_Map;
class Epetra_Comm;

//! Amesos_BaseSolver: A pure virtual class for direct solution of real-valued double-precision operators.
/*! 

The Amesos_BaseSolver class is a pure virtual class (that is, it specifies
interface only) that enables the use of real-valued
double-precision direct sparse solvers.  Every Amesos class named 
Amesos_<i>SolverName</i> derives from Amesos_BaseSolver.

<H1>Usage Examples</H1> 

<H2>Basic calling sequence</H2> The basic
calling sequence solves A x = b or A<SUP>T</SUP> x = b without
specifying how A has changed between each call to Solve().

\code
    Epetra_LinearProblem Problem(A,X,B);
    Amesos_SolverName Solver(Problem);
    Problem.SetUseTranspose(false); 
    while( ... ) { 

      <Here code which may change A or B>

      Solver.SymbolicFactorization(); 
      Solver.NumericFactorization(); 
      Solver.Solve();
    }
\endcode   
    

<H2>Re-using the symbolic factorization</H2> The following 
calling sequence performs multiple solves of A x = b or A<SUP>T</SUP> x = b 
in cases where the non-zero structure of A remains unchanged
between each call to Solve().

\code
    Epetra_LinearProblem Problem(A,X,B);
    Amesos_SolverName Solver(Problem);
    Problem.SetUseTranspose(false); 
    Solver.SymbolicFactorization(); 
    while( ... ) { 

      <Here code which may change B or the non-zero values
      (but not the non-zero structure) of A>

      Solver.NumericFactorization() ; 
      Solver.Solve() ; 
    }
\endcode
    

<H2>Pivot-less refactorization</H2> The following 
calling sequence performs multiple solves of A x = b or A<SUP>T</SUP> x = b 
provided that the change to A is small. 

\code
    Epetra_LinearProblem Problem(A,X,B);
    Amesos_SolverName Solver(Problem);
    Problem.SetTranspose( false ); 
    Solver.NumericFactorization(); 
    while( ... ) { 

      <Here code which may change A slightly>

      Solver.PivotlessRefactorization(); 
      Solver.Solve(); 
    }
\endcode


<H2>Re-using the numeric factorization</H2> The following 
calling sequence performs multiple solves of A x = b or A<SUP>T</SUP> x = b 
provided that A remains unchanged between each call to Solve().

\code
    Epetra_LinearProblem Problem(A,X,B);
    Amesos_SolverName Solver(Problem);
    Problem.SetTranspose( false ); 
    Solver.NumericFactorization(); 
    while( ... ) { 

      <Here code which may change B but not A>

      Solver.Solve();
    }
\endcode
    

<H2>Constructor requirements</H2>
 Every Amesos_SolverName class should accept an Epetra_LinearProblem 


<H2>Mathematical methods</H2> Four mathematical methods are
defined in the base class Amesos_BaseSolver:
SymbolicFactorization(), NumericFactorization(),
PivotelssRefactorization() and Solve().  


<H2>Switching concrete classes</H2>
Different concrete
classes, each based on a different third party solver, will have
different performance characteristics and will accept different
parameters.


<H2>Changing the underlying matrix operator.</H2> 
In the basic
calling sequence (no calls to SymbolicFactorization() or
NumericFactorization()), the underlying matrix can be modified
between any two calls to Solve() - any class implementing this
interface must perform the solve based on the values in the matrix
at the time that Solve() is called.  

<p>Once SymbolicFactorization() has been called, classes implementing
this interface may assume that any change made to the non-zero
structure of the underlying matrix will be accompanied by a call
to SymbolicFactorization() prior to a subsequent call to
NumericFactorization or Solve().


<H2>Named Parameters</H2>
    
Parameters can be changed or added at any time by calling 
SetParameters(ParamList) with the new parameters specified in ParamList.


It is left to the user to be sure that changes made to the
parameters are appropriate for the concrete class that they are
using.

Examples of appropriate changes in parameters include:
<ul>
<li><b>Changing iterative refinement rules between calls to Solve()</b>
<li><b>Changing drop tolerance rules between calls to NumericFactorization()</b>
</ul>

Examples of inappropriate changes in parameters include: <ul>
<li><b>Changing drop tolerance rules between solve steps.</b>
<pre> Solver.NumericFactorization();
Solver.getList()->set("DropTolerance",.001);
Solver.Solve(); </pre> </ul> Results of making inappropriate
changes in parameters is unpredictable and could include an error
return, a bogus result or ignoring the parameter change.
    
<H2>Transpose solve</H2> Any class implementing Amesos_BaseSolver
should handle calls to SetTranspose() at any point.  Some third
party libraries are able to solve A<SUP>T</SUP> x = b and Ax = b
using the same factorization.  Others will require a new
factorization anytime that a call to SetTranspose() changes the
intended solve from A<SUP>T</SUP> x = b to Ax = b or vice-versa.


<H1>Performance expectations</H1>

The following is a list of performance guidelines that classes 
which implement the Amesos_BaseSolver class are expected to maintain.


<H1>Memory usage:</H1>
For serial codes, no more than one extra copy of the original matrix 
should be required.  Except that some codes require matrix transpostion 
which requires additional copies of the input matrix.  

For distributed memory codes, no serial copies of the original matrix 
should be required.  


<H1>Robustness requirements</H1>

<p>Failures should be caught by AMESOS_CHK_ERR().  
   The following error codes should be used:
     - 1: Singular matrix
     - 2: Non-symmetric matrix
     - 3: Matrix is not positive definite
     - 4: Insufficient memory

<p>Because we do not check to see if a matrix has changed
between the call to SymbolicFactorization() and the call to
NumericFactorization(), it is possible that a change to the 
matrix will cause a potentially catastrophic error.  
*/    

class Amesos_BaseSolver {
      
 public:

  //@{ \name Destructor
  //! Destructor
  virtual ~Amesos_BaseSolver() {};
  
  //@}
  //@{ \name Mathematical functions

  //! Performs SymbolicFactorization on the matrix A.
  /*! 
      In addition to performing symbolic factorization on the matrix A, 
      the call to SymbolicFactorization() implies that no change will
      be made to the non-zero structure of the underlying matrix without 
      a subsequent call to SymbolicFactorization().
      
      <br \>Preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      </ul>

      <br \>Postconditions:<ul> 
      <li>Symbolic Factorization will be performed (or marked to be performed) 
      allowing NumericFactorization() and Solve() to be called.
      </ul>

    \return Integer error code, set to 0 if successful.
  */
  virtual int SymbolicFactorization() = 0;

  //! Performs NumericFactorization on the matrix A.
  /*!  In addition to performing numeric factorization on the 
      matrix A, the call to
      NumericFactorization() implies that no change will be made to
      the underlying matrix without a subsequent call to
      NumericFactorization().  

      <br \>Preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().  
          (return -2 if the number of non-zeros changes) 
          Other changes can have arbitrary consequences.
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      <li>The matrix should be indexed from 0 to n-1, unless the parameter "Reindex" 
          was set to "true" prior to the call to SymbolicFactorization().  
          (return -3 - if caught)
      <li>The paremeter "Reindex" should not be set to "true" except on CrsMatrices.
          (return -4) 
      <li>The paremeter "Reindex" should not be set to "true" unless Amesos was built
          with EpetraExt, i.e. with --enable-epetraext on the configure line.
          (return -4) 
      <li>Internal errors retur -5.
      </ul>

      <br \>Postconditions:<ul> 
      <li>Numeric Factorization will be performed (or marked to be performed) 
      allowing Solve() to be performed correctly despite a potential change in 
      in the matrix values (though not in the non-zero structure).
      </ul>

     \return Integer error code, set to 0 if successful.
  */
  virtual int NumericFactorization() = 0;


    //! Performs Pivotless Numeric Factorization on the matrix A.
    /*! PivotlessReFactorization() performs the same functionality 
as NumericFactorization() except that typically 
pivoting is performed.

<p>PivotlessFactorization() should compare the ratio of the largest pivot 
to the smallest pivot against the ratio after the most recent call to 
NumericFactorization(), if pivot growth in pivot-less factorization 
is too large, NumericFactorization() should be called by 
PivotlessFactorization().  Such calls to NumericFactorization()
are controlled by two parameters, either of which can be set to -1 to indicate no maximum.  
<ul>
  <li>"Maximum Relative Pivot Growth", default=10 
  <li>"Maximum Pivot Growth", default=-1 
</ul>

<p>In solvers which do not offer pivot-less factorization, 
PivotlessFactorization() defaults to a call to NumericFactorization().

<p>If no previous cal to NumericFactorization() has been made, 
PivotlessFactorization() defaults to a call to NumericFactorization().

<p>Ideally, a call to PivotlessFactorization() should not involve 
any memory allocation or re-allocation.


      <br \>Preconditions:<ul>
      <li>GetProblem().GetOperator() != 0 (return -1)
      <li>MatrixShapeOk(GetProblem().GetOperator()) == true (return -6)
      <li>The non-zero structure of the matrix should not have changed
          since the last call to SymbolicFactorization().  
      <li>The distribution of the matrix should not have changed 
          since the last call to SymbolicFactorization()
      <li>The change in the values in the matrix since the last call
          to NumericFactorization() should be small enough that 
	  pivot-less factorization is sufficiently accurate. 
          (Or, expected to be sufficiently accurate.)
      </ul>

      <br \>Postconditions:<ul> 
      <li>Numeric Factorization will be performed
      allowing Solve() to be performed correctly despite a 
      potential change in 
      in the matrix values (though not in the non-zero structure).
      </ul>

     \return Integer error code, set to 0 if successful.

    virtual int PivotlessFactorization() = 0;
  */

    //! Solves A X = B (or A<SUP>T</SUP> x = B) 
    /*! 

      <br \>Preconditions:<ul>
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

      <br \>Postconditions:<ul> 
      <li>X will be set such that A X = B (or
      A<SUP>T</SUP> X = B), within the limits of the accuracy of the
      underlying solver.  
      </ul>

     \return Integer error code, set to 0 if successful.
  */
  virtual int Solve() = 0;

  //@}
  //@{ \name Atribute set methods

  //! If set true, X will be set to the solution of A<SUP>T</SUP> X = B (not A X = B) 
  /*! If the implementation of this interface does not support
           transpose use, this method should return a value of -1.
      
      <br \>Preconditions:<ul>
      <li>SetTranspose() should be called prior to the call to SymbolicFactorization() If NumericFactorization() or Solve() is called after SetTranspose() without an intervening call to SymbolicFactorization() the result is implementation dependent. </li>
      </ul>

      <br \>Postconditions:<ul> 
      <li>The next factorization and solve will be performed with
      the new value of UseTranspose. </li>
      </ul>

    \param
	   UseTranspose -- (In) If true, solve A<SUP>T</SUP> X = B, otherwise
	   solve A X = B.  

    \return Integer error code, set to 0 if successful.  Set to -1 if
           this implementation does not support transpose.
  */
  virtual int SetUseTranspose(bool UseTranspose) = 0;
  
  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const = 0;

  //!  Updates internal variables. 
  /*!  
      <br \>Preconditions:<ul>
      <li>None.</li>
      </ul>

      <br \>Postconditions:<ul> 
      <li>Internal variables controlling the factorization and solve will
      be updated and take effect on all subseuent calls to NumericFactorization() 
      and Solve().</li>
      <li>All parameters whose value are to differ from the default values must 
be included in ParameterList.  Parameters not specified in ParameterList 
revert to their default values.
      </ul>

    \return Integer error code, set to 0 if successful. 
  */
  virtual int SetParameters( Teuchos::ParameterList &ParameterList ) = 0 ;

  //! Returns the Epetra_LinearProblem
  virtual const Epetra_LinearProblem* GetProblem() const = 0;

  //! Returns true if the solver can handle this matrix shape 
  /*! Returns true if the matrix shape is one that the underlying
    sparse direct solver can handle. Classes that work only on square
    matrices should return false for rectangular matrices.  Classes
    that work only on symmetric matrices whould return false for
    non-symmetric matrices.
  */
  virtual bool MatrixShapeOK() const = 0;

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const = 0;

  //@}

};

#endif /* _AMESOS_BASESOLVER_H_ */
