// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
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

#ifndef RTOPPACK_ROP_GET_SUB_VECTOR_DECL_HPP
#define RTOPPACK_ROP_GET_SUB_VECTOR_DECL_HPP


#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Reduction operator that extracts a sub-vector in the range of
 * global zero-based indexes [l,u].
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class ROpGetSubVector : public RTOpT<Scalar> {
public:

  /** \brief . */
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;

  /** \brief . */
  ROpGetSubVector( const index_type l = 0, const index_type u = 0 );

  /** \brief Set the range of global indexes to extract elements for. */
  void set_range( const index_type l, const index_type u );

  /** \brief Extract the subvector after all of the reductions are
   * completed.
   */
  const ConstSubVectorView<Scalar> operator()(
    const ReductTarget& reduct_obj ) const;

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void get_reduct_type_num_entries_impl(
    const Ptr<int> &num_values,
    const Ptr<int> &num_indexes,
    const Ptr<int> &num_chars
    ) const;
  /** \brief . */
  Teuchos::RCP<ReductTarget> reduct_obj_create_impl() const;
  /** \brief . */
  void reduce_reduct_objs_impl( const ReductTarget &in_reduct_obj,
    const Ptr<ReductTarget> &inout_reduct_obj) const;
  /** \brief . */
  void reduct_obj_reinit_impl( const Ptr<ReductTarget> &reduct_obj ) const;
  /** \brief . */
  void extract_reduct_obj_state_impl(
    const ReductTarget &reduct_obj,
    const ArrayView<primitive_value_type> &value_data,
    const ArrayView<index_type> &index_data,
    const ArrayView<char_type> &char_data
    ) const;
  /** \brief . */
  void load_reduct_obj_state_impl(
    const ArrayView<const primitive_value_type> &value_data,
    const ArrayView<const index_type> &index_data,
    const ArrayView<const char_type> &char_data,
    const Ptr<ReductTarget> &reduct_obj
    ) const;
  /** \brief . */
  bool coord_invariant_impl() const;
  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj
    ) const;

  //@}

private:

  index_type l_;
  index_type u_;

}; // class ROpGetSubVector


} // namespace RTOpPack


#endif // RTOPPACK_ROP_GET_SUB_VECTOR_DECL_HPP
