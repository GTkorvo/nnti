// ////////////////////////////////////////////////////////////////////
// MatrixOpNonsingThyra.hpp
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
//

#ifndef ALAP_MATRIX_OP_NONSING_Thyra_HPP
#define ALAP_MATRIX_OP_NONSING_Thyra_HPP

#include "AbstractLinAlgPack/src/abstract/interfaces/MatrixOpNonsing.hpp"
#include "MatrixOpThyra.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"

namespace AbstractLinAlgPack {

///
/** <tt>MatrixOpNonsing</tt> adapter subclass for <tt>Thyra::Nonlin::LinearOpWithSolve</tt>.
 */
class MatrixOpNonsingThyra
	:virtual public MatrixOpNonsing
	,virtual public MatrixOpThyra
{
public:

	/** @name Constructors / Initializers */
	//@{

	///
	/** Construct to uninitialized.
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_vec().get() == NULL</tt>
	 * <li>See the post conditions for <tt>MatrixOpThyra::MatrixOpThyra()</tt>
	 * </ul>
	 */
	MatrixOpNonsingThyra();
	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	MatrixOpNonsingThyra(
		const Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<value_type> >  &thyra_linear_op_ns
		,BLAS_Cpp::Transp                                                            thyra_linear_op_trans = BLAS_Cpp::no_trans
		);
	///
	/** Initalize given a smart pointer to a <tt>Thyra::LinearOpWithSolveBase</tt> object.
	 *
	 * @param  thyra_linear_op_ns  [in] Smart pointer to Thyra vector <tt>this</tt> will adapt.
	 *
	 * Preconditioins:<ul>
	 * <li><tt>thyra_linear_op_ns.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li><tt>thyra_linear_op_ns->opSupported(Thyra::NOTRANS) && thyra_linear_op_ns->opSupported(Thyra::TRANS)
	 *     (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_linear_op_ns().get() == thyra_linear_op_ns.get()</tt>
	 * <li>See the post conditions for <tt>MatrixOpThyra::initialize()</tt>
	 * </ul>
	 */
	void initialize(
		const Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<value_type> >  &thyra_linear_op_ns
		,BLAS_Cpp::Transp                                                            thyra_linear_op_trans = BLAS_Cpp::no_trans
		);
	///
	/** Set to uninitialized and return smart pointer to the internal <tt>Thyra::LinearOpWithSolveBase</tt> object.
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_linear_op_ns().get() == NULL</tt>
	 * </ul>
	 *
	 * Note that his nonvirtual function hides the nonvirtual function
	 * <tt>MatrixOpThyra::set_uninitialized()</tt>.
	 */
	Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<value_type> > set_uninitialized();
	///
	/** Return a smart pointer to the <tt>Thyra::LinearOpWithSolveBase</tt> object.
	 */
	Teuchos::RefCountPtr<const Thyra::LinearOpWithSolveBase<value_type> > thyra_linear_op_ns() const;

	//@}
  
	/** @name Overridden from MatrixOp (needed to remove ambiguities) */
	//@{

	/// Overridden to call <tt>MatrixOpThyra::clone()</tt>
	mat_mut_ptr_t clone();

	//@}

	/** @name Overridden from MatrixNonsing */
	//@{

	///
	void V_InvMtV(
		VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
		,const Vector& v_rhs2
		) const;
	///
	void M_StInvMtM(
		MatrixOp* m_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
		) const;

	//@}

	/** @name Overridden from MatrixOpNonsing */
	//@{

	///
	mat_mwons_ptr_t clone_mwons() const;

	//@}

};	// end class MatrixOpNonsingThyra

} // end namespace AbstractLinAlgPack

#endif	// ALAP_MATRIX_OP_NONSING_Thyra_HPP
