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

#include "lb_const.h"
#include "ha_const.h"

int LB_Get_Processor_Name(
   LB *lb,             /* The load-balancing structure.                */
   char *name          /* A string uniquely identifying the processor. 
                          We assume that at least MAX_PROC_NAME_LEN 
                          characters have been allocated.              */
)
{
/* This routine gets the name of the physical processor
   that an MPI process is running on.
 */

  int ierr = ZOLTAN_OK;
  int length;

  if (lb->Get_Processor_Name != NULL) {
    /* Use application-registered function */
    lb->Get_Processor_Name(lb->Get_Processor_Name_Data,
            name, &length, &ierr);
  }
  else {
    /* Use MPI_Get_processor_name by default */
    ierr = MPI_Get_processor_name(name, &length);
  }

  /* Add a trailing \0 to mark end of string */
  name[length] = '\0';

  return ierr;
}

