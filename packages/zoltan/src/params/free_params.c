/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "lb_const.h"


void Zoltan_Free_Params(
ZOLTAN_PARAM **params)				/* load balance structure */
{
/*
 * Free the list of new parameter values.
 */
    ZOLTAN_PARAM *ptr, *next;

    if (params == NULL) return;

    ptr = *params;
    while (ptr != NULL) {
	next = ptr->next;
	ZOLTAN_FREE(&(ptr->name));
	ZOLTAN_FREE(&(ptr->new_val));
	ZOLTAN_FREE(&ptr);
	ptr = next;
    }

    *params = NULL;
}
