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

#include <assert.h>
#include <string.h>
#include <malloc.h>

#include "RTOp_ROp_get_ele.h"
#include "RTOp_obj_index_vtbl.h"
#include "RTOp_obj_value_vtbl.h"
#include "RTOp_reduct_sum_value.h"

/* Implementation functions */

static int RTOp_ROp_get_ele_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  RTOp_index_type        i;
  RTOp_index_type        global_offset;
  RTOp_index_type        sub_dim;
  const RTOp_value_type  *v0_val;
  ptrdiff_t              v0_val_s;

  /* */
  /* Validate the input */
  /* */
  if( num_vecs != 1 )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 0 )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  assert(targ_obj);
  assert(vecs);

  /* */
  /* Get pointers to data */
  /* */

  /* i */
  i = *((RTOp_index_type*)obj_data);

  /* v0 */
  global_offset  = vecs[0].global_offset;
  sub_dim        = vecs[0].sub_dim;
  v0_val         = vecs[0].values;
  v0_val_s       = vecs[0].values_stride;

  /* */
  /* Lookup the element! */
  /* */

  if( i < global_offset + 1 || global_offset + sub_dim < i )
    return 0; /* The element we are looking for is not here. */

  *((RTOp_value_type*)targ_obj) = *(v0_val+v0_val_s*(i-global_offset-1));

  return 0; /* success? */
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_ROp_get_ele_vtbl =
{
  &RTOp_obj_index_vtbl  /* use simple scalar index type for object instance data */
  ,&RTOp_obj_value_vtbl /* use simple scalar value type for target object */
  ,"ROp_get_ele"
  ,NULL
  ,RTOp_ROp_get_ele_apply_op
  ,RTOp_reduct_sum_value
  ,RTOp_get_reduct_sum_value_op
};

/* Class specific functions */

int RTOp_ROp_get_ele_construct( RTOp_index_type i, struct RTOp_RTOp* op )
{
  op->vtbl = &RTOp_ROp_get_ele_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  *((RTOp_index_type*)op->obj_data) = i;
  return 0;
}

int RTOp_ROp_get_ele_set_i( RTOp_index_type i, struct RTOp_RTOp* op )
{
  *((RTOp_index_type*)op->obj_data) = i;
  return 0;
}

int RTOp_ROp_get_ele_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->vtbl     = NULL;
  return 0;
}

RTOp_value_type RTOp_ROp_get_ele_val(RTOp_ReductTarget targ_obj)
{
  return *((RTOp_value_type*)targ_obj);
}
