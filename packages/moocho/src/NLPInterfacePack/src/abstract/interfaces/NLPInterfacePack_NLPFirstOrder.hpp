// ///////////////////////////////////////////////////////////////////////
// NLPFirstOrder.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef NLP_FIRST_ORDER_INFO_H
#define NLP_FIRST_ORDER_INFO_H

#include "NLPObjGrad.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace NLPInterfacePack {
///
/** NLP first order information interface class {abstract}.
 *
 * <b>Overview:</b>
 *
 * This class adds Jacobian information for the constraints.  This augments the
 * information provided by the \c NLP and \c NLPObjGrad interfaces.  This interface
 * includes access to matrix space objects which must be used by the client
 * to create the matrix objects that are used with this interface.  This
 * totally decouples the client from the implementation of these matrix
 * objects.
 *
 * <b>Client Usage:</b>
 *
 * As with the <tt>NLP</tt> base interface, the <tt>initialize()</tt> method must be called before
 * the %NLP object can be used.
 *
 * The matrix space object returned from \c factory_Gc()must be used to create
 * the matrix objects for \c Gc used with this interface.  Note that the matrix objects
 * returned from this matrix space object can not be expected to be usable until they are
 * passed to the calculation routines.
 *
 * The method \c set_Gc() si used to set a pointer to matrix object to be updated
 * when \c Gc is computed using \c calc_Gc().
 *
 * The number of evaluations of \c Gc using calc_Gc() and calc_Gh() is returned by
 * \c num_Gc_evals().
 * 
 * <b>Subclass developer's notes:</b>
 *
 * In addition to the methods that must be overridden by the \c NLPObjGrad interface
 * (<A HREF="classNLPInterfacePack_1_1NLPObjGradient.html#must_override">see</A>) the following methods
 * must also be overridden: \c factory_Gc(), \c imp_calc_Gc().
 *
 * In addition to the methods that should be overridden from <tt>%NLPObjGrad</tt> by most subclasses
 * (<A HREF="classNLPInterfacePack_1_1NLPObjGradient.html#should_override">see</A>), the following
 * additional methods should be overridden: \c initialize().
 *
 * The following methods should never have to be overridden by most subclasses except in some very
 * specialized situations: \c set_Gc(), \c get_Gc(), \c Gc(), \c num_Gc_evals().
 */
class NLPFirstOrder : virtual public NLPObjGrad {
public:

	///
	typedef Teuchos::RefCountPtr<
		const Teuchos::AbstractFactory<MatrixOp> >        mat_fcty_ptr_t;
	///
	typedef Teuchos::RefCountPtr<BasisSystem>           basis_sys_ptr_t;

	/** @name Constructors */
	//@{

	/// Initialize to no reference set to calculation quanities
	NLPFirstOrder();

	//@}

	/** @name NLP initialization */
	//@{

	///
	/** Initialize the NLP for its first use.
	 *
	 * This function implementation should be called by subclass implementations
	 * in order to reset counts for \c f(x), \c c(x), \c h(x), \c Gf(x), \c Gc(x)
	 * and \c Gh(x) evaluations.  This implementation calls
	 * <tt>this->NLPObjGrad::initialize()</tt>
	 *
	 * Postconditions:<ul>
	 * <li> See <tt>NLPObjGrad::initialize()</tt>
	 * <li> <tt>this->num_Gc_evals() == 0</tt>
	 * <li> <tt>this->num_Gh_evals() == 0</tt>
	 * </ul>
	 */
	void initialize(bool test_setup);

	//@}

	/** @name Matrix factory objects */
	//@{
	
	///
	/** Return a matrix factory object for creating <tt>Gc</tt>.
	 *
	 * This method may return <tt>return.get() == NULL</tt> if <tt>m() == 0</tt>.
	 * Otherwise, it must must return a valid matrix factory object.
	 */
	virtual const mat_fcty_ptr_t factory_Gc() const = 0;

	//@}

	/** @name BasisSystem */
	//@{

	///
	/** Return a <tt>BasisSystem</tt> object compatible with <tt>Gc</tt> and <tt>Gh</tt>.
	 *
	 * Note that multiple calls to this method may return the same <tt>return.get()</tt>
	 * value so the client must not assume that they are unique.
	 *
	 * The default implementation returns <tt>return.get() == NULL</tt>.
	 */
	virtual const basis_sys_ptr_t basis_sys() const;

	//@}

	/** @name <<std aggr>> members for the gradient of the objective function Gc(x) */
	//@{

	///
	/** Set a pointer to a matrix object to be updated when <tt>this->calc_Gc()</tt> is called.
	 *
	 * @param  Gc  [in] Pointer to matrix of gradients.  May be \c NULL.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> [<tt>Gc != NULL</tt>] <tt>Gc->space().is_compatible(*this->factory_Gc(),no_trans) == true</tt>
	 *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_Gc() == Gc</tt>
	 * </ul>
	 */
	virtual void set_Gc(MatrixOp* Gc);
	///
	/** Return pointer passed to <tt>this->set_Gc()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual MatrixOp* get_Gc();
	///
	/** Returns non-<tt>const</tt> <tt>*this->get_Gc()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_Gc() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual MatrixOp& Gc();
	///
	/** Returns <tt>const</tt> <tt>*this->get_Gc()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_Gc() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual const MatrixOp& Gc() const;

	//@}

	/** @name Unset calculation quantities */
	//@{
	
	///
	/** Call to unset all storage quantities (both in this class and all subclasses).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> See <tt>NLPObjGrad::unset_quantities()</tt>
	 * <li> <tt>this->get_Gc() == NULL</tt>
	 * </ul>
	 *
	 * This method must be called by all subclasses that override it.
	 */
	void unset_quantities();

	//@}

	/** @name Calculation Members */
	//@{

	///
	/** Update the matrix for \c Gc at the point \c x and put it in the stored reference.
	 *
	 * @param  x     [in] Point at which to calculate the matrix of gradients <tt>Gc(x)</tt>.
	 * @param  newx  [in] (default \c true) If \c false, the values in \c x are assumed to be the same as
	 *               the last call to a <tt>this->imp_calc_*(x,newx)</tt> member.
	 *               If \c true, the values in \c x are assumed to not be the same as the last call to a
	 *               <tt>this->imp_calc_*(x,newx)</tt> member.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>this->get_Gc() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->Gc()</tt> is updated to \c Gc(x)
	 * </ul>
	 *
	 * The storage reference for <tt>f</tt> and/or <tt>c</tt> may also be updated at this point
	 * (if <tt>get_f() != NULL</tt> and/or <tt>get_c() != NULL</tt>)
	 * but is not guarentied to be.  But no other quanities from possible subclasses are allowed
	 * to be updated as a side effect.
	 */ 
	virtual void calc_Gc(const Vector& x, bool newx = true) const;

	//@}

	/** @name Function evaluation counts */
	//@{

	///
	/** Gradient of constraints matrix \c Gc evaluations count.
	 *
	 * This function can be called to find out how many evaluations
	 * \c this->calc_Gc() the client requested since \c this->initialize() was called.
	 */
	virtual size_type num_Gc_evals() const;

	//@}

protected:

	///
	/** Struct for zero and first order quantities (pointers)
	 */
	struct FirstOrderInfo {
		///
		FirstOrderInfo()
			: Gc(NULL), Gf(NULL), f(NULL), c(NULL)
		{}
		///
		FirstOrderInfo( MatrixOp* Gc_in, const ObjGradInfo& obj_grad )
			: Gc(Gc_in), Gf(obj_grad.Gf), f(obj_grad.f), c(obj_grad.c)
		{}
		/// Pointer to Jacobian of equality constraints <tt>Gc</tt> (may be NULL if not set)
		MatrixOp*           Gc;
		/// Pointer to gradient of objective function <tt>Gf</tt> (may be NULL if not set)
		VectorMutable*      Gf;
		/// Pointer to objective function <tt>f</tt> (may be NULL if not set)
		value_type*         f;
		/// Pointer to equality constraints residule <tt>c</tt> (may be NULL if not set)
		VectorMutable*      c;
	}; // end struct FirstOrderInfo

	/// Return objective gradient and zero order information.
	const FirstOrderInfo first_order_info() const;

	/** @name Protected methods to be overridden by subclasses */
	//@{

	///
	/** Overridden to compute \a Gc(x) and perhaps \a Gf(x), \a f(x) and \a c(x).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
	 * <li> <tt>obj_grad_info.Gc != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*obj_grad_info.Gc</tt> is updated to \a Gc(x).
	 * </ul>
	 *
	 * @param x       [in]  Unknown vector (size n).
	 * @param  newx   [in] (default \c true) If \c false, the values in \c x are assumed to be the same as
	 *                the last call to a <tt>this->imp_calc_*(x,newx)</tt> member.
	 *                If \c true, the values in \c x are assumed to not be the same as the last call to a
	 *                <tt>this->imp_calc_*(x,newx)</tt> member.
	 * @param obj_grad_info
	 *                [out] Pointers to \c f, \c c, \c Gf and \c Gc
	 *                On output <tt>*obj_grad_info.Gc</tt> is updated to \a Gc(x).
	 *                Any of the other objects pointed to in \c obj_grad_info may
	 *                also be set but are now guaranteed to be.
	 */
	virtual void imp_calc_Gc(const Vector& x, bool newx, const FirstOrderInfo& first_order_info) const = 0;

	//@}

private:

#ifdef DOXYGEN_COMPILE
	AbstractLinAlgPack::BasisSystem                            *basis_sys;
	Teuchos::AbstractFactory<AbstractLinAlgPack::MatrixOp>  *factory_Gc;
#endif
	mutable MatrixOp      *Gc_;
	mutable size_type     num_Gc_evals_;

};	// end class NLPFirstOrder

// /////////////////////
// Inline members

inline
const NLPFirstOrder::FirstOrderInfo NLPFirstOrder::first_order_info() const
{
	return FirstOrderInfo(Gc_,obj_grad_info());
}

}	// end namespace NLPInterfacePack 

#endif // NLP_FIRST_ORDER_INFO_H
