// ////////////////////////////////////////////////////////////////////
// RTOp_ROp_max_step.h
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

#ifndef RTOP_ROP_MAX_STEP_H
#define RTOP_ROP_MAX_STEP_H

#include "RTOp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** @name Reduction operator for finding the maximum step for feasibility.
 *
 * <tt>targ_obj <- { max alpha | v[0] + alpha * v[1] >= beta }</tt>
 *
 * This is a specialized reduction operation that is used in many optimization
 * methods to find the maximum step length \c alpha such that the
 * iterates remain feasible.  It is assumed that <tt>v[0] > beta</tt>
 * so that <tt>alpha == 0.0</tt> is a valid return.  If the step is
 * unrestricted, the function \c RTOp_ROp_max_step_val() will return
 * \c RTOp_ROp_max_step_inf.
 *
 * This operator is defined to allow exactly two vector arguments
 * (<tt>num_vecs == 2</tt>) \c v[0], \c v[1], and can only handle
 * dense vectors.
 */
//@{

/// Name of this reduction operator class
extern const char RTOp_ROp_max_step_name[];

/// Virtual function table
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_max_step_vtbl;

/// Constructor
int RTOp_ROp_max_step_construct( RTOp_value_type beta, struct RTOp_RTOp* op );

/// Destructor
int RTOp_ROp_max_step_destroy( struct RTOp_RTOp* op );

/// Reset beta
int RTOp_ROp_max_step_set_beta( RTOp_value_type beta, struct RTOp_RTOp* op );

/// Value of max step if step is unrestricted.
extern RTOp_value_type  RTOp_ROp_max_step_inf;

/// Extract the concrete value from a reduction target object.
RTOp_value_type  RTOp_ROp_max_step_val(RTOp_ReductTarget reduct_obj);

//@}

#ifdef __cplusplus
}
#endif

#endif  // RTOP_ROP_MAX_STEP_H
