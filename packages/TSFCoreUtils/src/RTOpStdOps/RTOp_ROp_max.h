/* ///////////////////////////////////////////// */
/* RTOp_ROp_max.h */

/* */
/* Note: This file was created automatically by 'new_rtop.pl' */
/*       on 7/15/2002 at 16:59 */
/* */

#ifndef RTOp_ROp_max_H
#define RTOp_ROp_max_H

#include "RTOp.h"
#include "RTOp_obj_value_vtbl.h"  /* vtbl for reduction object data */

#ifdef __cplusplus
extern "C" {
#endif

/** \file RTOp_ROp_max.h
 *
 \verbatim

 max_ele = max{ v0(i), i = 1...n }
 \endverbatim
 *
 * This operator class implementation was created
 * automatically by 'new_rtop.pl'.
 */
/*@{ */

/* Virtual function table */
extern const struct RTOp_RTOp_vtbl_t RTOp_ROp_max_vtbl;

/* Constructor */
int RTOp_ROp_max_construct(  struct RTOp_RTOp* op );

/* Destructor */
int RTOp_ROp_max_destroy( struct RTOp_RTOp* op );


/* Extract the value of the reduction object */
RTOp_value_type RTOp_ROp_max_val(RTOp_ReductTarget reduct_obj);

/*@} */

#ifdef __cplusplus
}
#endif

#endif  /* RTOp_ROp_max_H */
