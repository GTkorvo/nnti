// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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

#ifndef ANASAZI_MULTI_VEC_HPP
#define ANASAZI_MULTI_VEC_HPP

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziConfigDefs.hpp"

/*! 	\class Anasazi::MultiVec

	\brief Anasazi's templated pure virtual class for constructing multivectors that 
	are used by the eigensolver.

	A concrete implementation of this class is necessary.  The user can create
	their own implementation if those supplied are not suitable for their needs.

	\author Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {

template <class TYPE>
class MultiVec {
public:
	//@{ \name Constructor/Destructor.
	//! %Anasazi::MultiVec constructor.
	MultiVec() {};

	//! %Anasazi::MultiVec destructor.
	virtual ~MultiVec () {};

	//@}
	//@{ \name Creation methods for new multivectors.

	/*! \brief Creates a new empty %Anasazi::MultiVec containing \c numvecs columns.

	    \return Pointer to the new multivector	
	*/

	virtual MultiVec<TYPE> * Clone ( const int numvecs ) const = 0;

	/*! \brief Creates a new %Anasazi::MultiVec and copies contents of \c *this into
	    the new vector (deep copy).
	
	    \return Pointer to the new multivector	
	*/
	
	virtual MultiVec<TYPE> * CloneCopy () const = 0;
	
	/*! \brief Creates a new %Anasazi::MultiVec and copies the selected contents of \c *this 
	    into the new vector (deep copy).  The number (\c numvecs) of copied 
	    vectors from \c *this are indicated by the indices in \c index.

	    \return Pointer to the new multivector	
	*/

	virtual MultiVec<TYPE> * CloneCopy ( int index[], int numvecs ) const = 0;
	
	/*! \brief Creates a new %Anasazi::MultiVec that shares the selected contents of \c *this.
	    The index of the \c numvecs vectors copied from \c *this are indicated by the
	    indices given in \c index.

	    \return Pointer to the new multivector	
	*/

	virtual MultiVec<TYPE> * CloneView ( int index[], int numvecs ) = 0;
	//@}

	//@{ \name Dimension information methods.	
	//! Obtain the vector length of *this multivector block.

	virtual int GetVecLength () const = 0;

	//! Obtain the number of vectors in *this multivector block.

	virtual int GetNumberVecs () const = 0;

	//@}
	//@{ \name Update methods.
	/*! \brief Update \c *this with \c alpha * \c A * \c B + \c beta * (\c *this).
	*/

	virtual void MvTimesMatAddMv ( const TYPE alpha, const MultiVec<TYPE>& A, 
		const Teuchos::SerialDenseMatrix<int,TYPE>& B, const TYPE beta ) = 0;

	/*! \brief Replace \c *this with \c alpha * \c A + \c beta * \c B.
	*/

	virtual void MvAddMv ( const TYPE alpha, const MultiVec<TYPE>& A, const TYPE beta, const MultiVec<TYPE>& B ) = 0;

	/*! \brief Compute a dense matrix \c B through the matrix-matrix multiply 
	   \c alpha * \c A^T * (\c *this).
	*/

	virtual void MvTransMv ( const TYPE alpha, const MultiVec<TYPE>& A, Teuchos::SerialDenseMatrix<int,TYPE>& B) const = 0;

	/*! \brief Compute a vector \c b where the components are the individual dot-products, i.e.\c b[i] = \c A[i]^T*\c this[i] where \c A[i] is the i-th column of A.
	*/

	virtual void MvDot ( const MultiVec<TYPE>& A, TYPE b[] ) const = 0;

	//@}
	//@{ \name Norm method.

	/*! \brief Compute the 2-norm of each individual vector of \c *this.  
	   Upon return, \c normvec[i] holds the 2-norm of the \c i-th vector of \c *this
	*/

	virtual void MvNorm ( TYPE* normvec ) const = 0;

	//@}
	//@{ \name Initialization methods.
	/*! \brief Copy the vectors in \c A to a set of vectors in \c *this.  The \c 
  	    numvecs vectors in \c A are copied to a subset of vectors in \c *this
	    indicated by the indices given in \c index.
	*/

	virtual void SetBlock ( const MultiVec<TYPE>& A, int index[], int numvecs ) = 0;
	
	/*! \brief Replace the vectors in \c *this with random vectors.
	*/

	virtual void MvRandom () = 0;

	/*! \brief Replace each element of the vectors in \c *this with \c alpha.
	*/

	virtual void MvInit ( const TYPE alpha ) = 0;

	//@}
	//@{ \name Print method.
	/*! \brief Print the \c *this multivector.
	*/
	virtual void MvPrint () const = 0;
	//@}
};


  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Anasazi::MultiVecTraits for Anasazi::MultiVec.
  //
  ////////////////////////////////////////////////////////////////////


  template<class TYPE>
  class MultiVecTraits<TYPE,MultiVec<TYPE> >
  {
  public:
    ///
    static Teuchos::RefCountPtr<MultiVec<TYPE> > Clone( const MultiVec<TYPE>& mv, const int numvecs )
    { return Teuchos::rcp( const_cast<MultiVec<TYPE>&>(mv).Clone(numvecs) ); }
    ///
    static Teuchos::RefCountPtr<MultiVec<TYPE> > CloneCopy( const MultiVec<TYPE>& mv )
    { return Teuchos::rcp( const_cast<MultiVec<TYPE>&>(mv).CloneCopy() ); }
    ///
    static Teuchos::RefCountPtr<MultiVec<TYPE> > CloneCopy( const MultiVec<TYPE>& mv, int index[], int numvecs )
    { return Teuchos::rcp( const_cast<MultiVec<TYPE>&>(mv).CloneCopy(index,numvecs) ); }
    ///
    static Teuchos::RefCountPtr<MultiVec<TYPE> > CloneView( MultiVec<TYPE>& mv, int index[], int numvecs )
    { return Teuchos::rcp( mv.CloneView(index,numvecs) ); }
    ///
    static Teuchos::RefCountPtr<const MultiVec<TYPE> > CloneView( const MultiVec<TYPE>& mv, int index[], int numvecs )
    { return Teuchos::rcp( const_cast<MultiVec<TYPE>&>(mv).CloneView(index,numvecs) ); }
    ///
    static int GetVecLength( const MultiVec<TYPE>& mv )
    { return mv.GetVecLength(); }
    ///
    static int GetNumberVecs( const MultiVec<TYPE>& mv )
    { return mv.GetNumberVecs(); }
    ///
    static void MvTimesMatAddMv( TYPE alpha, const MultiVec<TYPE>& A, 
				 const Teuchos::SerialDenseMatrix<int,TYPE>& B, 
				 TYPE beta, MultiVec<TYPE>& mv )
    { mv.MvTimesMatAddMv(alpha, A, B, beta); }
    ///
    static void MvAddMv( TYPE alpha, const MultiVec<TYPE>& A, TYPE beta, const MultiVec<TYPE>& B, MultiVec<TYPE>& mv )
    { mv.MvAddMv(alpha, A, beta, B); }
    ///
    static void MvTransMv( TYPE alpha, const MultiVec<TYPE>& A, const MultiVec<TYPE>& mv, Teuchos::SerialDenseMatrix<int,TYPE>& B )
    { mv.MvTransMv(alpha, A, B); }
    ///
    static void MvDot( const MultiVec<TYPE>& mv, const MultiVec<TYPE>& A, TYPE b[] )
    { mv.MvDot( A, b ); }
    ///
    static void MvNorm( const MultiVec<TYPE>& mv, TYPE *normvec )
    { mv.MvNorm(normvec); }
    ///
    static void SetBlock( const MultiVec<TYPE>& A, int index[], int numvecs, MultiVec<TYPE>& mv )
    { mv.SetBlock(A, index, numvecs); }
    ///
    static void MvRandom( MultiVec<TYPE>& mv )
    { mv.MvRandom(); }
    ///
    static void MvInit( MultiVec<TYPE>& mv, TYPE alpha = Teuchos::ScalarTraits<TYPE>::zero() )
    { mv.MvInit(alpha); }
    ///
    static void MvPrint( const MultiVec<TYPE>& mv, ostream& os )
    { mv.MvPrint(os); }
    
  };


} // namespace Anasazi

#endif

// end of file AnasaziMultiVec.hpp
