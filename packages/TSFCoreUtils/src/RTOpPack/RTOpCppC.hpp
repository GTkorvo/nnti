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

#ifndef RTOP_CPP_C_H
#define RTOP_CPP_C_H

#include "RTOpCpp.hpp"

namespace RTOpPack {

///
/** C++ subclass for RTOp_RTOp objects.
 *
 * ToDo: Finish Documentation
 */
class RTOpC : public RTOp {
public:

	///
	RTOpC( const RTOp_RTOp_vtbl_t* = NULL );
	///
	~RTOpC();
	///
	RTOp_RTOp& op();
	///
	const RTOp_RTOp& op() const;

	/** @name Overridden from RTOp */
	//@{

	///
	const char* op_name() const;
	///
	op_create_func_t get_op_create_func() const;
	///
	op_free_func_t get_op_free_func() const;
	///
	void get_op_type_num_entries(
		int*  num_values
		,int* num_indexes
		,int* num_chars
		) const;
	///
	void extract_op_state(
		int               num_values
		,RTOp_value_type  value_data[]
		,int              num_indexes
		,RTOp_index_type  index_data[]
		,int              num_chars
		,RTOp_char_type   char_data[]
		) const;
	///
	void load_op_state(
		int                       num_values
		,const RTOp_value_type    value_data[]
		,int                      num_indexes
		,const RTOp_index_type    index_data[]
		,int                      num_chars
		,const RTOp_char_type     char_data[]
		);
	///
	void get_reduct_type_num_entries(
		int*   num_values
		,int*  num_indexes
		,int*  num_chars
		) const;
	///
	void reduct_obj_create_raw( RTOp_ReductTarget* reduct_obj ) const;
	///
	void reduct_obj_reinit( RTOp_ReductTarget reduct_obj ) const;
	///
	void reduct_obj_free( RTOp_ReductTarget* reduct_obj ) const;
	///
	void extract_reduct_obj_state(
		const RTOp_ReductTarget   reduct_obj
		,int                      num_values
		,RTOp_value_type          value_data[]
		,int                      num_indexes
		,RTOp_index_type          index_data[]
		,int                      num_chars
		,RTOp_char_type           char_data[]
		) const;
	///
	void load_reduct_obj_state(
		int                      num_values
		,const RTOp_value_type   value_data[]
		,int                     num_indexes
		,const RTOp_index_type   index_data[]
		,int                     num_chars
		,const RTOp_char_type    char_data[]
		,RTOp_ReductTarget       reduct_obj
		) const;
	///
	void apply_op(
		const int   num_vecs,       const SubVector         sub_vecs[]
		,const int  num_targ_vecs,  const MutableSubVector  targ_sub_vecs[]
		,RTOp_ReductTarget reduct_obj
		) const;
	///
	void reduce_reduct_objs(
		RTOp_ReductTarget in_reduct_obj, RTOp_ReductTarget inout_reduct_obj
		) const;
	///
    void get_reduct_op( RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr ) const;

	//@}

private:

	const RTOp_RTOp_vtbl_t  *vtbl_;
	RTOp_RTOp               op_;

}; // end class RTOpC

// //////////////////////////////////////////
// Inline member functions

// Inline members for RTOpC

inline
RTOp_RTOp& RTOpC::op()
{
	return op_;
}

inline
const RTOp_RTOp& RTOpC::op() const
{
	return op_;
}

} // end namespace RTOpPack

#endif // RTOP_CPP_C_H
