// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
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

#ifndef TSFCORE_SERIAL_MULTI_VECTOR_STD_DECL_HPP
#define TSFCORE_SERIAL_MULTI_VECTOR_STD_DECL_HPP

#include "TSFCoreSerialMultiVectorBaseDecl.hpp"

namespace TSFCore {

///
/** \brief Highly optimized concrete implementation subclass for
 * serial shared-memory multi-vectors.
 *
 * This subclass provides a very efficient and very general concrete
 * implementation of a <tt>TSFCore::MultiVector</tt> object.
 *
 * Objects of this type generally should not be constructed directly
 * by a client but instead by using the concrete vector space subclass
 * <tt>TSFCore::SerialVectorSpace</tt> and using the function
 * <tt>TSFCore::SerialVectorSpace::createMembers()</tt>.
 *
 * The storage type can be anything since a
 * <tt>Teuchos::RefCountPtr</tt> is used to pass in the values pointer
 * into the constructor and <tt>initialize()</tt>.
 *
 * \ingroup TSFCore_adapters_serial_concrete_std_grp
 */
template<class Scalar>
class SerialMultiVectorStd : virtual public SerialMultiVectorBase<Scalar> {
public:

  ///
  using SerialMultiVectorBase<Scalar>::subView;
  ///
  using SerialMultiVectorBase<Scalar>::col;

  /** @name Constructors/initializers/accessors */
  //@{

	/// Construct to uninitialized
	SerialMultiVectorStd();

  /// Calls <tt>initialize()</tt>
	SerialMultiVectorStd(
    const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
    ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
    );

  /// Calls <tt>initialize()</tt>
	SerialMultiVectorStd(
    const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
    ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
    ,const Teuchos::RefCountPtr<Scalar>                                    &values
    ,const Index                                                           leadingDim
    );

  ///
  /** Initialize.
   *
   * @param  range     [in] Smart pointer to <tt>SerialVector</tt> object
   *                   that defines the data distribution for <tt>mpiSpace()</tt> and <tt>range()</tt>.
   * @param  domain    [in] Smart pointer to <tt>ScalarProdVectorSpaceBase</tt> object
   *                   that defines <tt>domain()</tt> space.
   *
   * Preconditions:<ul>
   * <li><tt>range.get()!=NULL</tt>
   * <li><tt>domain.get()!=NULL</tt>
   * </ul>
	 *
	 * This function simply calls <tt>initialize(range,domain,...)</tt>
	 * and passes in dynamically allocated data for the values.
   */
	void initialize(
    const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
    ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
    );

  ///
  /** Initialize.
   *
   * @param  range     [in] Smart pointer to <tt>SerialVector</tt> object
   *                   that defines the data distribution for <tt>mpiSpace()</tt> and <tt>range()</tt>.
   * @param  domain    [in] Smart pointer to <tt>ScalarProdVectorSpaceBase</tt> object
   *                   that defines <tt>domain()</tt> space.
   * @param  values    [in] Smart pointer to beginning of Fortran-style column-major
   *                   array that defines the local values in the multi-vector.
   *                   This array must be at least of dimension <tt>range->leadingDim()*domain->dim()</tt>
   *                   and <tt>(&*values)[ (i-1) + (j-1)*leadingDim ]</tt> gives the local value
   *                   of the one-based <tt>(i,j)</tt> entry where <tt>i=1...mpiSpace()->localSubDim()</tt>
   *                   and <tt>j=1...domain->dim()</tt>.
	 * @param  leadingDim
	 *                   [in] The leading dimension of the multi-vector.
   *
   * Preconditions:<ul>
   * <li><tt>range.get()!=NULL</tt>
   * <li><tt>domain.get()!=NULL</tt>
   * <li><tt>values.get()!=NULL</tt>
   * <li><tt>leadingDim >= range->localSubDim()</tt>
   * </ul>
   */
	void initialize(
    const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >   &range
    ,const Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >  &domain
    ,const Teuchos::RefCountPtr<Scalar>                                    &values
    ,const Index                                                           leadingDim
    );

  ///
  /** Set to an uninitialized state.
   *
   * Postconditions:<ul>
   * <li><tt>this->mpiSpace().get() == NULL</tt>.
   */
	void uninitialize(
    Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >         *range         = NULL
    ,Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >        *domain        = NULL
    ,Teuchos::RefCountPtr<Scalar>                                          *values        = NULL
    ,Index                                                                 *leadingDim    = NULL
    );

  //@}

	/** @name Overridden from EuclideanLinearOpBase */
	//@{
	///
	Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const;
	///
	Teuchos::RefCountPtr< const ScalarProdVectorSpaceBase<Scalar> > domainScalarProdVecSpc() const;
  //@}

	/** @name Overridden from MultiVector */
	//@{
	///
	Teuchos::RefCountPtr<Vector<Scalar> > col(Index j);
	///
	Teuchos::RefCountPtr<MultiVector<Scalar> > subView( const Range1D& col_rng );
	//@}

	/** @name Overridden from SerialMultiVectorBase */
	//@{

	///
	void getData( const Scalar **values, Index *leadingDim ) const;
	///
	void freeData( const Scalar *values ) const;
	///
	void getData( Scalar **values, Index *leadingDim );
	///
	void commitData( Scalar *values );
	//@}
	
private:
	
	// ///////////////////////////////////////
	// Private data members

  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >     range_;
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> >     domain_;
  Teuchos::RefCountPtr<Scalar>                                       values_;
  Index                                                              leadingDim_;
	Index                                                              numRows_; // Cached
	Index                                                              numCols_; // Cached
	
}; // end class SerialMultiVectorStd

} // end namespace TSFCore

#endif // TSFCORE_SERIAL_MULTI_VECTOR_STD_DECL_HPP
