/*
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
*/

#ifndef RTOP_REDUCT_SUM_VALUES_H
#define RTOP_REDUCT_SUM_VALUES_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @name Definitions of reduction functions for summing a list of scalar values.
 *
 * These functions perform a simple sum of a list of scalar objects
 * as defined by the virtual function table \Ref{RTOp_obj_values_vtbl}.
 */
/*@{ */

/* */
/** Use this function for <tt>reduce_reduct_objs</tt> in the RTOp_RTOp_vtbl_t virtual
 * function table.
 */
int RTOp_reduct_sum_values(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
	, RTOp_ReductTarget in_targ_obj, RTOp_ReductTarget inout_targ_obj );

/* */
/** Use this function for <tt>get_reduct_op</tt> in the RTOp_RTOp_vtbl_t virtual
 * function table.
 */
int RTOp_get_reduct_sum_values_op(
	const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
	, RTOp_reduct_op_func_ptr_t* reduct_op_func_ptr );

/*@} */

#ifdef __cplusplus
}
#endif

#endif /* RTOP_REDUCT_SUM_VALUES_H */
