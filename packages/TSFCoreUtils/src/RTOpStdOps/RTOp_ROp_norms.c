// /////////////////////////////////////////////
// RTOp_ROp_norms.c
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

#include <assert.h>
#include <math.h>

#include "RTOp_ROp_norms.h"
#include "RTOp_obj_null_vtbl.h"
#include "RTOp_obj_value_vtbl.h"
#include "RTOp_reduct_sum_value.h"
#include "RTOp_reduct_max_value.h"

#define MY_MAX(a,b) a > b ? a : b

// One norm reduction operator class.

static int RTOp_ROp_norms_apply_op_norm_1(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  RTOp_index_type        sub_dim;
  const RTOp_value_type  *v0_val;
  ptrdiff_t              v0_val_s;
  register RTOp_index_type k;
  RTOp_value_type        *norm = NULL;

  //
  // Validate the input
  //
  if( num_vecs != 1 )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 0 )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  assert(targ_obj);
  assert(vecs);

  //
  // Get pointers to data
  //

  // v0
  sub_dim        = vecs[0].sub_dim;
  v0_val         = vecs[0].values;
  v0_val_s       = vecs[0].values_stride;

  //
  // Perform the reduction
  //
  norm = (RTOp_value_type*)targ_obj;
  for( k = 0; k < sub_dim; ++k, v0_val += v0_val_s ) {
    *norm += fabs(*v0_val);  // ||v[0]||_1
  }

  return 0; // success?
}

const char RTOp_ROp_norm_1_name[] = "RTOp_ROp_norm_1";

const struct RTOp_RTOp_vtbl_t RTOp_ROp_norm_1_vtbl =
{
  &RTOp_obj_null_vtbl   // use null type for instance data
  ,&RTOp_obj_value_vtbl // use simple scalar type for target object
  ,NULL
  ,RTOp_ROp_norms_apply_op_norm_1
  ,RTOp_reduct_sum_value
  ,RTOp_get_reduct_sum_value_op
};

int RTOp_ROp_norm_1_construct( struct RTOp_RTOp* op )
{
  op->vtbl     = &RTOp_ROp_norm_1_vtbl;
  op->obj_data = NULL;
  return 0;
}

RTOp_value_type RTOp_ROp_norm_1_val(RTOp_ReductTarget targ_obj)
{
#ifdef RTOp_DEBUG
  assert(targ_obj != RTOp_REDUCT_OBJ_NULL );
#endif
  return *(RTOp_value_type*)targ_obj;
}

// Two norm reduction operator class.

static int RTOp_ROp_norms_apply_op_norm_2(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  RTOp_index_type        sub_dim;
  const RTOp_value_type  *v0_val;
  ptrdiff_t              v0_val_s;
  register RTOp_index_type k;
  RTOp_value_type        *norm = NULL;

  //
  // Validate the input
  //
  if( num_vecs != 1 )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 0 )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  assert(targ_obj);
  assert(vecs);

  //
  // Get pointers to data
  //

  // v0
  sub_dim        = vecs[0].sub_dim;
  v0_val         = vecs[0].values;
  v0_val_s       = vecs[0].values_stride;

  //
  // Perform the reduction
  //
  norm = (RTOp_value_type*)targ_obj;
  for( k = 0; k < sub_dim; ++k, v0_val += v0_val_s )
    *norm += (*v0_val)*(*v0_val);  // (||v[0]||_2)^2

  return 0; // success?
}

const char RTOp_ROp_norm_2_name[] = "RTOp_ROp_norm_2";

const struct RTOp_RTOp_vtbl_t RTOp_ROp_norm_2_vtbl =
{
  &RTOp_obj_null_vtbl   // use null type for instance data
  ,&RTOp_obj_value_vtbl // use simple scalar type for target object
  ,NULL
  ,RTOp_ROp_norms_apply_op_norm_2
  ,RTOp_reduct_sum_value
  ,RTOp_get_reduct_sum_value_op
};

int RTOp_ROp_norm_2_construct( struct RTOp_RTOp* op )
{
  op->vtbl     = &RTOp_ROp_norm_2_vtbl;
  op->obj_data = NULL;
  return 0;
}

RTOp_value_type RTOp_ROp_norm_2_val(RTOp_ReductTarget targ_obj)
{
#ifdef RTOp_DEBUG
  assert(targ_obj != RTOp_REDUCT_OBJ_NULL );
#endif
  return sqrt(*(RTOp_value_type*)targ_obj);
}

// Infinity norm reduction operator class.

static int RTOp_ROp_norms_apply_op_norm_inf(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  RTOp_index_type        sub_dim;
  const RTOp_value_type  *v0_val;
  ptrdiff_t              v0_val_s;
  register RTOp_index_type k;
  RTOp_value_type        *norm = NULL;

  //
  // Validate the input
  //
  if( num_vecs != 1 )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 0 )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;
  assert(targ_obj);
  assert(vecs);

  //
  // Get pointers to data
  //

  // v0
  sub_dim        = vecs[0].sub_dim;
  v0_val         = vecs[0].values;
  v0_val_s       = vecs[0].values_stride;

  //
  // Perform the reduction
  //
  norm = (RTOp_value_type*)targ_obj;
  for( k = 0; k < sub_dim; ++k, v0_val += v0_val_s )
    *norm = MY_MAX( fabs(*v0_val), (*(RTOp_value_type*)targ_obj) );  // ||v[0]||_inf

  return 0; // success?
}

const char RTOp_ROp_norm_inf_name[] = "RTOp_ROp_norm_inf";

const struct RTOp_RTOp_vtbl_t RTOp_ROp_norm_inf_vtbl =
{
   &RTOp_obj_null_vtbl  // use null type for instance data
  ,&RTOp_obj_value_vtbl // use simple scalar type for target object
  ,NULL
  ,RTOp_ROp_norms_apply_op_norm_inf
  ,RTOp_reduct_max_value
  ,RTOp_get_reduct_max_value_op
};

int RTOp_ROp_norm_inf_construct( struct RTOp_RTOp* op )
{
  op->vtbl     = &RTOp_ROp_norm_inf_vtbl;
  op->obj_data = NULL;
  return 0;
}

RTOp_value_type RTOp_ROp_norm_inf_val(RTOp_ReductTarget targ_obj)
{
#ifdef RTOp_DEBUG
  assert(targ_obj != RTOp_REDUCT_OBJ_NULL );
#endif
  return *(RTOp_value_type*)targ_obj;
}

// Common functions

int RTOp_ROp_norm_destroy( struct RTOp_RTOp* op )
{
  op->vtbl = NULL;
  return 0;
}
