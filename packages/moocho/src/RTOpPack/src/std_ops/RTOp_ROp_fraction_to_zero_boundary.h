/*
// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/

/* ///////////////////////////////////////////// */
/* RTOp_ROp_fraction_to_zero_boundary.h */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 7/8/2002 at 19:19 */
/* */

#ifndef RTOp_ROp_fraction_to_zero_boundary_H
#define RTOp_ROp_fraction_to_zero_boundary_H

#include "RTOp.h"
#include "RTOp_obj_value_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_fraction_to_zero_boundary.h
 *
 \verbatim

element-wise reduction:
    alpha_max_ith = (v1 > 0) ? 1.0 : tau*v0/(-v1+1e-50);


 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 *
 * ToDo: Write the documentation for this class!
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_fraction_to_zero_boundary_vtbl;

/* Constructor */
int RTOp_ROp_fraction_to_zero_boundary_construct( RTOp_value_type tau,  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_fraction_to_zero_boundary_destroy( struct RTOp_RTOp* op );

/* Initialize the state of the operator object */
int RTOp_ROp_fraction_to_zero_boundary_init( RTOp_value_type tau, struct RTOp_RTOp* op );

/* Extract the value of the reduction object */
RTOp_value_type RTOp_ROp_fraction_to_zero_boundary_val(RTOp_ReductTarget reduct_obj);

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_ROp_fraction_to_zero_boundary_H */
