// @HEADER
// ***********************************************************************
// 
//      TSFCoreUtils: Trilinos Solver Framework Utilities Package 
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

// //////////////////////////////////////////////
// RTOpPack_MPI_apply_op_decl.hpp

#ifndef RTOPPACK_MPI_APPLY_OP_DECL_HPP
#define RTOPPACK_MPI_APPLY_OP_DECL_HPP

#include "RTOpPack_RTOpT.hpp"

namespace RTOpPack {

///
/** Initialize MPI compatible type signature arrays for
 * reduction/transformation operator object instance data and
 * reduction target object data.
 *
 * @param num_values    [in] Number of primitive_value members
 * @param num_indexes   [in] Number of index members
 * @param num_chars     [in] Number of character members
 * @param num_entries   [out] Number of entries in output arrays set
 * @param block_lengths [out] array (length<tt> >= RTOp_NUM_DATA_TYPES</tt>)
 * @param displacements [out] array (length<tt> >= RTOp_NUM_DATA_TYPES</tt>)
 * @param datatypes     [out] array (length<tt> >= RTOp_NUM_DATA_TYPES</tt>)
 *
 * See the MPI function <tt>MPI_Type_struct(...)</tt> for a discription of these arrays.
 */
template<class primitive_value_type>
void MPI_type_signature(
	const int num_values
	,const int num_indexes
	,const int num_chars
	,int* num_entries
	,int block_lengths[]
	,MPI_Aint displacements[]
	,MPI_Datatype datatypes[]
	);

/// Return the size in bytes of an external representation of <tt>reduct_obj</tt>
template<class primitive_value_type>
int reduct_obj_ext_size(
	int   num_values
	,int  num_indexes
	,int  num_chars
	)
{
  return (3 + num_values) * sizeof(primitive_value_type)
    + num_indexes         * sizeof(index_type)
    + num_chars           * sizeof(char_type);
}

///
template<class Scalar>
void extract_reduct_obj_ext_state(
  const RTOpT<Scalar>    &op
	,const ReductTarget    &reduct_obj
	,int                   num_values
	,int                   num_indexes
	,int                   num_chars
	,void                  *reduct_obj_ext
	);

///
template<class Scalar>
void load_reduct_obj_ext_state(
  const RTOpT<Scalar>    &op
	,const void            *reduct_obj_ext
	,ReductTarget          *reduct_obj
	);

///
template<class Scalar>
void MPI_apply_op(
	MPI_Comm                                      comm
  ,const RTOpT<Scalar>                          &op
  ,const int                                    root_rank
	,const int                                    num_vecs
  ,const RTOpPack::SubVectorT<Scalar>           sub_vecs[]
	,const int                                    num_targ_vecs
  ,const RTOpPack::MutableSubVectorT<Scalar>    targ_sub_vecs[]
	,ReductTarget                                 *reduct_obj
	);

///
template<class Scalar>
void MPI_apply_op(
	MPI_Comm                                           comm
  ,const RTOpT<Scalar>                               &op
  ,const int                                         root_rank
	,const int                                         num_cols
	,const int                                         num_multi_vecs
  ,const RTOpPack::SubMultiVectorT<Scalar>           sub_multi_vecs[]
	,const int                                         num_targ_multi_vecs
  ,const RTOpPack::MutableSubMultiVectorT<Scalar>    targ_sub_multi_vecs[]
	,ReductTarget*                                     reduct_objs[]
	);

///
template<class Scalar>
void  MPI_apply_op(
	MPI_Comm                                  comm
  ,const RTOpT<Scalar>                      &op
  ,int                                      root_rank
	,const int                                num_cols
	,const int                                num_vecs
  ,const SubVectorT<Scalar>                 sub_vecs[]
	,const int                                num_targ_vecs
  ,const MutableSubVectorT<Scalar>          sub_targ_vecs[]
	,ReductTarget*                            reduct_objs[]
  );

} // end namespace RTOpPack

#endif // RTOPPACK_MPI_APPLY_OP_DECL_HPP
