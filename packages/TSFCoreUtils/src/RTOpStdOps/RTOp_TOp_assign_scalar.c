/* ///////////////////////////////////////////// */
/* RTOp_TOp_assign_scalar.c */
/* */
/* Copyright (C) 2001 Roscoe Ainsworth Bartlett */
/* */
/* This is free software; you can redistribute it and/or modify it */
/* under the terms of the "Artistic License" (see the web site */
/*   http://www.opensource.org/licenses/artistic-license.html). */
/* This license is spelled out in the file COPYING. */
/* */
/* This software is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* above mentioned "Artistic License" for more details. */

#include <assert.h>
#include <malloc.h>

#include "RTOp_TOp_assign_scalar.h"
#include "RTOp_obj_value_vtbl.h"
#include "RTOp_obj_null_vtbl.h"

/* Implementation functions for RTOp_RTOp */

static int RTOp_TOp_assign_scalar_apply_op(
  const struct RTOp_RTOp_vtbl_t* vtbl, const void* obj_data
  , const int num_vecs, const struct RTOp_SubVector vecs[]
  , const int num_targ_vecs, const struct RTOp_MutableSubVector targ_vecs[]
  , RTOp_ReductTarget targ_obj )
{
  /* */
  /* Get pointers to data */
  /* */

  /* alpha */
  const RTOp_value_type alpha = *((RTOp_value_type*)obj_data);

  /* z */
  RTOp_index_type        z_sub_dim     = targ_vecs[0].sub_dim;
  RTOp_value_type        *z_val        = targ_vecs[0].values;
  ptrdiff_t              z_val_s       = targ_vecs[0].values_stride;

  register RTOp_index_type k;

  /* */
  /* Validate the input */
  /* */

  if( num_vecs != 0 || vecs != NULL )
    return RTOp_ERR_INVALID_NUM_VECS;
  if( num_targ_vecs != 1 || targ_vecs == NULL )
    return RTOp_ERR_INVALID_NUM_TARG_VECS;

  /* */
  /* Assign the elements */
  /* */

  for( k = 0; k < z_sub_dim; ++k, z_val += z_val_s )
    *z_val = alpha;

  return 0; /* success? */
}

/* Virtual function table */
const struct RTOp_RTOp_vtbl_t RTOp_TOp_assign_scalar_vtbl =
{
  &RTOp_obj_value_vtbl /* use simple scalar value type for object instance data */
  ,&RTOp_obj_null_vtbl /* use null type for target object */
  ,"TOp_assign_scalar"
  ,NULL /* use default from reduct_vtbl */
  ,RTOp_TOp_assign_scalar_apply_op
  ,NULL
  ,NULL
};

/* Class specific functions */

int RTOp_TOp_assign_scalar_construct( RTOp_value_type alpha, struct RTOp_RTOp* op )
{
  op->vtbl = &RTOp_TOp_assign_scalar_vtbl;
  op->vtbl->obj_data_vtbl->obj_create(NULL,NULL,&op->obj_data);
  *((RTOp_value_type*)op->obj_data) = alpha;
  return 0; /* success? */
}

int RTOp_TOp_assign_scalar_destroy( struct RTOp_RTOp* op )
{
  op->vtbl->obj_data_vtbl->obj_free(NULL,NULL,&op->obj_data);
  op->vtbl = NULL;
  return 0; /* success? */
}

int RTOp_TOp_assign_scalar_set_alpha( RTOp_value_type alpha, struct RTOp_RTOp* op )
{
  *((RTOp_value_type*)op->obj_data) = alpha;
  return 0; /* success? */
}
