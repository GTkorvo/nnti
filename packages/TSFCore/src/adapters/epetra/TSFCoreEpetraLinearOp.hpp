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

// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraLinearOp.hpp

#ifndef TSFCORE_EPETRA_LINEAR_OP_HPP
#define TSFCORE_EPETRA_LINEAR_OP_HPP

#include "TSFCoreLinearOp.hpp"
#include "ReleaseResource.hpp"

class Epetra_Operator;

namespace TSFCore {

class EpetraVectorSpace;

///
/** Implementation of <tt>LinearOp</tt> using an <tt>Epetra_Operator</tt> object.
 *
 * This subclass can be used to represent the non-transposed operator
 * or transposed operator defined by an <tt>Epetra_Operator</tt>
 * object.  This class assumes that both the non-transposed and
 * transposed applications of the operator can be performed.
 */
class EpetraLinearOp : public LinearOp<RTOp_value_type> {
public:

	/** @name Public types */
	//@{

	///
	typedef RTOp_value_type Scalar;

	//@}

	/** @name Constructors / initializers / accessors */
	//@{

	///
	/** Construct to uninitialized.
	 *
	 * See the postconditions for <tt>setUninitialized()</tt>
	 */
	EpetraLinearOp();

	/// Calls <tt>initialize()</tt>.
	EpetraLinearOp(
		const Teuchos::RefCountPtr<Epetra_Operator>   &op
		,ETransp                                      opTrans = NOTRANS
		);

	///
	/** Initialize
	 *
	 * @param  op       [in] The <tt>Epetra_Operator</tt> this <tt>*this</tt> will wrap.
	 * @param  opTrans  [in] If <tt>opTrans==NOTRANS</tt> then <tt>op</tt> will be viewed as <tt>op</tt>
	 *                  and if <tt>opTrans==TRANS</tt> then <tt>op</tt> will be viewed as its transpose
	 *                  <tt>op'</tt> for the behavior of <tt>apply()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>op.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->domain().get() != NULL</tt>
	 * <li> <tt>this->range().get() != NULL</tt>
	 * </ul>
	 */
	void initialize(
		const Teuchos::RefCountPtr<Epetra_Operator>   &op
		,ETransp                                      opTrans = NOTRANS
		);
	
	///
	/** Set to uninitialized and return the current state.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->domain().get() == NULL</tt>
	 * <li> <tt>this->range().get() == NULL</tt>
	 * </ul>
	 */
	void setUninitialized(
		Teuchos::RefCountPtr<Epetra_Operator>    *op      = NULL
		,ETransp                                 *opTrans = NULL
		);

  ///
  /** Return a smart pointer to the EpetraVectorSpace object for range.
   *
	 * Postconditions:<ul>
	 * <li> [<tt>this->range().get() != NULL</tt>] <tt>return.get() != NULL</tt>
	 * <li> [<tt>this->range().get() == NULL</tt>] <tt>return.get() == NULL</tt>
	 */
	Teuchos::RefCountPtr<const EpetraVectorSpace> epetraRange() const;

  ///
  /** Return a smart pointer to the EpetraVectorSpace object for domain.
   *
	 * Postconditions:<ul>
	 * <li> [<tt>this->domain().get() != NULL</tt>] <tt>return.get() != NULL</tt>
	 * <li> [<tt>this->domain().get() == NULL</tt>] <tt>return.get() == NULL</tt>
	 */
	Teuchos::RefCountPtr<const EpetraVectorSpace> epetraDomain() const;

  ///
  /** Return a smart pointer to the Epetra_Operator object. Added by KL
   * to support operator loading in TSFExtended */
	Teuchos::RefCountPtr<Epetra_Operator> epetra_op();

  /** Return a smart pointer to the Epetra_Operator object. Added by KL
   * to support operator loading in TSFExtended */
	const Teuchos::RefCountPtr<Epetra_Operator>& epetra_op() const ;
	//@}
	
	/** @name Overridden from OpBase */
	//@{
	
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > range() const;
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > domain() const;
	
	//@}
	
	/** @name Overridden from LinearOp */
	//@{
	
	///
	Teuchos::RefCountPtr<const LinearOp<Scalar> > clone() const;
	///
	void apply(
		const ETransp            M_trans
		,const Vector<Scalar>    &x
		,Vector<Scalar>          *y
		,const Scalar            alpha
		,const Scalar            beta
		) const;
	///
	void apply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha
		,const Scalar                 beta
		) const;

	//@}

  
protected:
  // ////////////////////////////////////
  // protected methods


  /** \name Allocators for domain and range spaces */
  //@{
  /** Allocate the domain space of the operator. Purpose: In
   * TSFExtended, both EpetraLinearOp and EpetraVectorSpace are
   * extended from the TSFCore versions by inheritance, and the
   * TSFExtended operator subclasses expect to work with an extended
   * vector space subclass. Thus, it is necessary for the base
   * operator class to never directly allocate vector space objects,
   * and allocation is delegated to a virtual allocator function. 
   * KRL and RAB, 2/18/04. */
  virtual Teuchos::RefCountPtr<const EpetraVectorSpace> 
  allocateDomain(const Teuchos::RefCountPtr<Epetra_Operator>  &op 
                 ,ETransp  op_trans 
                 )  const ; 
  
  /** Allocate the range space of the operator.Purpose: In
   * TSFExtended, both EpetraLinearOp and EpetraVectorSpace are
   * extended from the TSFCore versions by inheritance, and the
   * TSFExtended operator subclasses expect to work with an extended
   * vector space subclass. Thus, it is necessary for the base
   * operator class to never directly allocate vector space objects,
   * and allocation is delegated to a virtual allocator function. 
   * KRL and RAB, 2/18/04. */
  virtual Teuchos::RefCountPtr<const EpetraVectorSpace> allocateRange( 
    const Teuchos::RefCountPtr<Epetra_Operator>  &op 
    ,ETransp  op_trans 
    )  const ; 
  //@}
private:

	// ////////////////////////////////////
	// Private data members

	Teuchos::RefCountPtr<Epetra_Operator>          op_;
	ETransp                                        opTrans_;
	Teuchos::RefCountPtr<const EpetraVectorSpace>  domain_;
	Teuchos::RefCountPtr<const EpetraVectorSpace>  range_;

};	// end class EpetraLinearOp

}	// end namespace TSFCore

#endif	// TSFCORE_EPETRA_LINEAR_OP_HPP
