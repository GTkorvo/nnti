// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef MATRIX_SYM_HESSIAN_RELAX_NON_SING_H
#define MATRIX_SYM_HESSIAN_RELAX_NON_SING_H

#include "ConstrainedOptPack_Types.hpp"
#include "AbstractLinAlgPack_MatrixSymDiagStd.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace ConstrainedOptPack {

/** \brief Matrix class for non-singular Hessian matrix augmented with a terms for
 * "Big M" relaxation variables.
 *
 * The matrix that is formed is:
 \begin{verbatim}

 H = [ G    ]
     [    M ]
 \end{verbatim}
 * where <tt>M</tt> is a diagonal matrix made up of entries <tt>M_diag</tt>.
 */
class MatrixSymHessianRelaxNonSing
  : public AbstractLinAlgPack::MatrixSymOpNonsing
{
public:

  /** \brief . */
  typedef Teuchos::RefCountPtr<const MatrixSymOpNonsing>  G_ptr_t;
  /** \brief . */
  typedef Teuchos::RefCountPtr<VectorMutable>               vec_mut_ptr_t;
  /** \brief . */
  typedef Teuchos::RefCountPtr<const VectorSpace>                 space_ptr_t;
  
  /** @name Constructors/initializers */
  //@{

  /** \brief Construct uninitialized.
   *
   * Postconditions:
   * <ul>
   * <li> <tt>this->rows() == 0</tt>
   * <li> <tt>&this->G_ptr().get() == NULL</tt>
   * <li> <tt>&this->M_diag_ptr().get() == NULL</tt>
   * </ul>
   */
  MatrixSymHessianRelaxNonSing();

  /** \brief Constructor (calls \c initialize()).
   */
  MatrixSymHessianRelaxNonSing(
    const G_ptr_t         &G_ptr
    ,const vec_mut_ptr_t  &M_diag_ptr
    ,const space_ptr_t    &space = Teuchos::null
    );

  /** \brief Initialize the Hessian and the relaxation terms.
   *
   * Preconditions:
   * <ul>
   * <li> <tt>G_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>G_ptr->rows() > 0</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>M_diag_ptr.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * <li> <tt>M_diag_ptr->dim() > 0</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   *
   * Postconditions:
   * <ul>
   * <li> <tt>this->rows() == G_ptr->rows() + M_diag->dim()</tt>
   * <li> <tt>&this->G() == G_ptr.get()</tt>
   * <li> <tt>&this->M().diag() == M_diag_ptr.get()</tt>
   * </ul>
   *
   * @param  G_ptr   [in] Smart pointer to matrix that this object will represent.  The underlying
   *                 matrix object <tt>*G_ptr.get()</tt> should not be modified without calling \c this->initialize()
   *                 again.
   * @param  M_diag_ptr
   *                 [in] Smart pointer to the diagonal for <tt>M</tt>.  All of the elements in this vector must
   *                 be nonzero!  This vector is used to initalize the diagonal matrix <tt>M</tt> and this vector should
   *                 not be modified externally without calling \c this->initialize() again.
   * @param  space   [in] Smart pointer to the space used for <tt>this->space_cols()</tt> and <tt>this->space_rows()</tt>.
   *                 If <tt>space.get() == NULL</tt> then <tt>VectorSpaceCompoisteStd</tt> will be used which is
   *                 constructed from the spaces <tt>G_ptr->space_cols()</tt> and <tt>M_diag_ptr->space()</tt>.
   */
  void initialize(
    const G_ptr_t         &G_ptr
    ,const vec_mut_ptr_t  &M_diag_ptr
    ,const space_ptr_t    &space = Teuchos::null
    );

  /** \brief Set uninitalized.
   *
   */
  void set_uninitialized();
  
  /** \brief . */
  const G_ptr_t& G_ptr() const;

  /** \brief . */
  const vec_mut_ptr_t& M_diag_ptr() const;

  /** \brief . */
  const MatrixSymOpNonsing& G() const;

  /** \brief . */
  const AbstractLinAlgPack::MatrixSymDiagStd& M() const;

  //@}
  
  /** @name Overridden from MatrixOp */
  //@{

  /** \brief . */
  const VectorSpace& space_cols() const;
  /** \brief . */
  bool Mp_StM(
    MatrixOp* mwo_lhs, value_type alpha
    , BLAS_Cpp::Transp trans_rhs) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StMtV(
    VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
    ,const SpVectorSlice& sv_rhs2, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(
    VectorMutable* v_lhs, value_type alpha
    ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    ,BLAS_Cpp::Transp M_rhs2_trans
    ,const Vector& v_rhs3, value_type beta) const;
  /** \brief . */
  void Vp_StPtMtV(
    VectorMutable* v_lhs, value_type alpha
    ,const GenPermMatrixSlice& P_rhs1, BLAS_Cpp::Transp P_rhs1_trans
    ,BLAS_Cpp::Transp M_rhs2_trans
    ,const SpVectorSlice& sv_rhs3, value_type beta) const;

  //@}

  /** @name Overridden form MatrixSymOp */
  //@{

  /** \brief . */
  void Mp_StPtMtP(
    MatrixSymOp* sym_lhs, value_type alpha
    ,EMatRhsPlaceHolder dummy_place_holder
    ,const GenPermMatrixSlice& gpms_rhs, BLAS_Cpp::Transp gpms_rhs_trans
    ,value_type beta
    ) const;

  //@}

  /** @name Overridden from MatrixNonsing */
  //@{

  /** \brief . */
  void V_InvMtV(
    VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const Vector& v_rhs2) const;
  /** \brief . */
  void V_InvMtV(
    VectorMutable* v_lhs, BLAS_Cpp::Transp trans_rhs1
    ,const SpVectorSlice& sv_rhs2) const;

  //@}

private:
  
  // ///////////////////////////////
  // Private data members

  space_ptr_t              vec_space_;
  G_ptr_t                  G_ptr_;
  MatrixSymDiagStd     M_;

  // ///////////////////////////////
  // Private member functions

  void assert_initialized() const;

};

// ////////////////////////////////////
// Inline members

inline
const MatrixSymHessianRelaxNonSing::G_ptr_t&
MatrixSymHessianRelaxNonSing::G_ptr() const
{
  return G_ptr_;
}

inline
const MatrixSymHessianRelaxNonSing::vec_mut_ptr_t&
MatrixSymHessianRelaxNonSing::M_diag_ptr() const
{
  return M_.diag_ptr();
}

inline
const MatrixSymOpNonsing&
MatrixSymHessianRelaxNonSing::G() const
{
  assert_initialized();
  return *G_ptr_;
}

inline
const MatrixSymDiagStd&
MatrixSymHessianRelaxNonSing::M() const
{
  assert_initialized();
  return M_;
}
  
} // end namespace ConstrainedOptPack

#endif // MATRIX_SYM_HESSIAN_RELAX_NON_SING_H
