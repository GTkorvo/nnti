//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef TPETRA_MATRIXMATRIX_DECL_HPP
#define TPETRA_MATRIXMATRIX_DECL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultSparseMultiply.hpp>
#include <Kokkos_DefaultSparseSolve.hpp>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"

/*! \file Tpetra_MMMultiply_decl.hpp 

    The declarations for the class Tpetra::MMMultiMultiply and related non-member constructors.
 */

namespace Tpetra {
#ifndef DOXYGEN_SHOULD_SKIP_THIS  
	// forward declaration
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class SpMatVec, class SpMatSlv>
class CrsMatrix;
#endif
  /** Collection of matrix-matrix operations. This class basically
      functions as a namespace, containing only static methods.
      See the program epetraext/test/MatrixMatrix/cxx_main.cpp for
      a usage example.
   */
template <class Scalar, 
	class LocalOrdinal=int, 
	class GlobalOrdinal=LocalOrdinal, 
	class Node=Kokkos::DefaultNode::DefaultNodeType, 
	class SpMatVec=Kokkos::DefaultSparseMultiply<Scalar, LocalOrdinal, Node>, 
	class SpMatSlv=Kokkos::DefaultSparseSolve<Scalar, LocalOrdinal, Node> >
class MatrixMatrix {
  typedef CrsMatrix<LocalOrdinal,
  					GlobalOrdinal,
					Node,
					SpMatVec,
					SpMatSlv> CrsMatixType;
  public:
    /** destructor */
    virtual ~MatrixMatrix(){}

    /** Given CrsMatrix objects A, B and C, form the product C = A*B.
	In a parallel setting, A and B need not have matching distributions,
	but C needs to have the same row-map as A.

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param C Result. On entry to this method, it doesn't matter whether
             FillComplete() has already been called on C or not. If it has,
	     then C's graph must already contain all nonzero locations that
	     will be produced when forming the product A*B. On exit,
	     C.FillComplete() will have been called, unless the last argument
             to this function is specified to be false.
    @param call_FillComplete_on_result Optional argument, defaults to true.
           Power users may specify this argument to be false if they *DON'T*
           want this function to call C.FillComplete. (It is often useful
           to allow this function to call C.FillComplete, in cases where
           one or both of the input matrices are rectangular and it is not
           trivial to know which maps to use for the domain- and range-maps.)

    @return error-code, 0 if successful. non-zero returns may result if A or
             B are not already Filled, or if errors occur in putting values
             into C, etc.
     */
    static int Multiply(const CrsMatrixType& A,
			bool transposeA,
			const CrsMatrixType& B,
			bool transposeB,
			CrsMatrixType& C,
                        bool call_FillComplete_on_result=true);

    /** Given CrsMatrix objects A and B, form the sum B = a*A + b*B

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param scalarA Input, scalar multiplier for matrix A.
    @param B Result. On entry to this method, it doesn't matter whether
             FillComplete() has already been called on B or not. If it has,
	     then B's graph must already contain all nonzero locations that
	     will be produced when forming the sum.
    @param scalarB Input, scalar multiplier for matrix B.

    @return error-code, 0 if successful. non-zero returns may result if A is
             not already Filled, or if errors occur in putting values
             into B, etc.
     */
    static int Add(const CrsMatrixType& A,
                   bool transposeA,
                   Scalar scalarA,
                   CrsMatrixType& B,
                   Scalar scalarB);

    /** Given CrsMatrix objects A and B, form the sum C = a*A + b*B

    @param A Input, must already have had 'FillComplete()' called.
    @param transposeA Input, whether to use transpose of matrix A.
    @param scalarA Input, scalar multiplier for matrix A.
    @param B Input, must already have had 'FillComplete()' called.
    @param transposeB Input, whether to use transpose of matrix B.
    @param scalarB Input, scalar multiplier for matrix B.
    @param C Result. On entry to this method, C can be NULL or a pointer
             to an unfilled or filled matrix. If C is NULL then a new
             object is allocated and must be deleted by the user.
             If C is not NULL and FillComplete has already
             been called then the sparsity pattern is assumed to be fixed
             and compatible  with the sparsity of A+B. If FillComplete has
             not been called then the sum is completed and the function
             returns without calling FillComplete on C.

    @return error-code, 0 if successful. non-zero returns may result if A or is
             not already Filled, or if errors occur in putting values
             into C, etc.
     */
    static int Add(const CrsMatrixType& A,
                   bool transposeA,
                   Scalar scalarA,
                   const CrsMatrixType& B,
                   bool transposeB,
                   Scalar scalarB,
                   CrsMatrixType * & C);

};//class MatrixMatrix


/**
 *Method for internal use... sparsedot forms a dot-product between two
 *sparsely-populated 'vectors'.
 *Important assumption: assumes the indices in u_ind and v_ind are sorted.
 */
 //double sparsedot(double* u, int* u_ind, int u_len,
//		  double* v, int* v_ind, int v_len);
 template<Scalar, LocalOrdinal>
 Scalar sparsedot(Teuchos::ArrayRCP<Scalar> u, Teuchos::ArrayRCP<LocalOrdinal> u_ind, 
		  Teuchos::ArrayRCP<Scalar> v, Teuchos::ArrayRCP<LocalOrdinal> v_ind);


}
#endif // TPETRA_MATRIXMATRIX_DECL_HPP

