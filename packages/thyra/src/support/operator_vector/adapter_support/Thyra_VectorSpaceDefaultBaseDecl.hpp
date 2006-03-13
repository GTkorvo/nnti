// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_VECTOR_SPACE_DEFAULT_BASE_DECL_HPP
#define THYRA_VECTOR_SPACE_DEFAULT_BASE_DECL_HPP

#include "Thyra_OperatorVectorAdapterSupportTypes.hpp"
#include "Thyra_VectorSpaceBaseDecl.hpp"

namespace Thyra {

/** \brief Node <tt>VectorSpaceBase</tt> subclass that provides default
 * implementations for many functions using a default multi-vectors
 * implementation.
 *
 * <b>Notes to Subclass Developers</b>
 *
 * Because of the default multi-vector implementation used in this node
 * subclasses, a concrete subclass is only required to override three methods:
 * <tt>dim()</tt>, <tt>isCompatible()</tt> and <tt>createMember()</tt>.  Note
 * that implementing the <tt>createMember()</tt> method also entails defining
 * a concrete <tt>VectorBase</tt> subclass.
 *
 * If a subclass can support specialized multi-vectors, then the
 * <tt>createMembers()</tt> function should be overridden as well.  Note that
 * implementing <tt>createMembers()</tt> also entails defining a concrete
 * <tt>MultiVectorBase</tt> subclass.  For some types of concrete
 * <tt>MultiVectorBase</tt> subclass implementations (e.g. serial
 * multi-vectors), the same default <tt>VectorSpaceFactoryBase</tt> typed
 * object returned from the <tt>smallVecSpcFcty()</tt> override provided in
 * this node subclass can be reused in the concrete subclass.  However, more
 * specialized <tt>MultiVectorBase</tt> subclasses (e.g. distributed memory
 * parallel) may require an override of the <tt>smallVecSpcFcty()</tt> method
 * to return a specialized type of <tt>VectorSpaceFactoryBase</tt> object.
 *
 * \ingroup Thyra_Op_Vec_Adapters_grp
 */
template<class Scalar>
class VectorSpaceDefaultBase : virtual public VectorSpaceBase<Scalar> {
public:

  /** @name Public functions overridden from VectorSpaceBase */ 
  //@{

  /** \brief .
   *
   * The default implementation returns <tt>dynamic_cast<
   * SerialVectorSpaceFactory<Scalar> >(return.get())!=NULL</tt>.  Note that
   * if a subclass overrides <tt>createMembers()</tt> then it may also need to
   * override this method as well.
   */
  Teuchos::RefCountPtr< const VectorSpaceFactoryBase<Scalar> > smallVecSpcFcty() const;

  //@}

protected:

  /** @name Protected functions overridden from VectorSpaceBase */
  //@{

  /** \brief .
   *
   * The default implementation returns <tt>dynamic_cast<
   * DefaultColumnwiseMultiVector<Scalar> >(return.get())!=NULL</tt>.
   */
  Teuchos::RefCountPtr< MultiVectorBase<Scalar> > createMembers(int numMembers) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>this->createMember()</tt> to create a vector then uses the explicit
   * element access functions to set the elements and then only when the
   * vector is destroyed is the data copied out of the vector and back into
   * the elements pointed to by <tt>raw_v.values()</tt>.
   */
  Teuchos::RefCountPtr<VectorBase<Scalar> > createMemberView( const RTOpPack::SubVectorView<Scalar> &raw_v ) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>this->createMember()</tt> to create a vector then uses the explicit
   * element access functions to set the elements from
   * <tt>raw_v.values()</tt>.
   */
  Teuchos::RefCountPtr<const VectorBase<Scalar> > createMemberView( const RTOpPack::ConstSubVectorView<Scalar> &raw_v ) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>this->createMembers(raw_mv.numSubCols())</tt> to create a
   * multi-vector then uses the explicit element access functions to set the
   * elements and then only when the multi-vector is destroyed is the data
   * copied out of the multi-vector and back into the elements pointed to by
   * <tt>raw_mv.values()</tt>.
   */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > createMembersView( const RTOpPack::SubMultiVectorView<Scalar> &raw_mv ) const;

  /** \brief .
   *
   * The default implementation of this function simply calls
   * <tt>this->createMembers()</tt> to create a multi-vector then uses the
   * explicit element access functions to set the elements.
   */
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > createMembersView( const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv ) const;

  //@}

}; // end class VectorSpaceBase

} // end namespace Thyra

#endif  // THYRA_VECTOR_SPACE_DEFAULT_BASE_DECL_HPP
