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
#include "par_const.h"


#define PRINT_SYNC 5000   /* definition needed for print sync */

/************ R O U T I N E S   I N   T H I S   F I L E  **********************

       NAME                             TYPE
----------------------------------------------------------------------
	LB_Print_Sync_Start		void
	LB_Print_Sync_End		void

******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_Print_Sync_Start(MPI_Comm communicator, int do_print_line)
{
/* 
 * Routine to allow I/O between LB_Print_Sync_Start and LB_Print_Sync_End to be 
 * printed by each processor in the communicator entirely before the next
 * processor begins its I/O.  The printing sequence is from proc = 0 to the
 * last processor, where the last processor is num_proc - 1.
 *
 * The do_print_line argument is a boolean variable.  If true, a line of # 
 * is printed to indicate the start of a Print_Sync I/O block.
 *
 * NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.
 *
 * Author: John Shadid (9221, SNL)
 */

int        flag = 1, from, type;
static int offset = 0;
MPI_Status st;
char *yo = "LB_Print_Sync_Start";
int proc;

  MPI_Comm_rank(communicator, &proc);

  offset = (offset + 1)%100;
  type   = PRINT_SYNC + offset;

  if (proc != 0) {
    from = proc -1;
    if (MPI_Recv((void *) &flag, 1, MPI_INT, from, type, communicator, &st)
        != 0) {
      fprintf(stderr, "%s: ERROR on processor %d\n", yo, proc);
      fprintf(stderr, "MPI_Recv failed, message type %d\n", type);
      exit (-1);
    }
  }
  else {
    if (do_print_line) {
      printf("\n");
      for (flag = 0; flag < 37; flag++) printf("#");
      printf(" PRINT_SYNC_START ");
      for (flag = 0; flag < 25; flag++) printf("#");
      printf("\n");
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_Print_Sync_End(MPI_Comm communicator, int do_print_line)
{
/*
 * Routine to allow I/O between LB_Print_Sync_Start and LB_Print_Sync_End to be 
 * printed by each processor in the communicator entirely before the next
 * processor begins its I/O.  The printing sequence is from proc = 0 to the
 * last processor, where the last processor is num_proc - 1.
 *
 * The do_print_line argument is a boolean variable.  If true, a line of # 
 * is printed to indicate the start of a Print_Sync I/O block.
 *
 * NOTE: THERE CAN BE NO COMMUNICATON BETWEEN THESE CALLS.
 *
 * Author: John Shadid (9221, SNL)
 */

int         flag = 1, from, type, to;
static int  offset = 0;
MPI_Status  st;
int proc, num_proc;
char *yo = "LB_Print_Sync_End";

  MPI_Comm_rank(communicator, &proc);
  MPI_Comm_size(communicator, &num_proc);

  fflush(stdout);

  offset = (offset + 1)%100;
  type   = PRINT_SYNC + offset;

  if (proc < num_proc -1)
    to = proc + 1;
  else {
    to = 0;
    if (do_print_line) {
      printf("\n");
      for (flag = 0; flag < 37; flag++) printf("#");
      printf(" PRINT_SYNC_END__ ");
      for (flag = 0; flag < 25; flag++) printf("#");
      printf("\n\n");
    }
  }

  if (MPI_Send((void *) &flag, 1, MPI_INT, to, type, communicator) != 0 ) {
    fprintf(stderr, "%s: ERROR on node %d\n", yo, proc);
    fprintf(stderr, "MPI_Send failed, message type %d\n", type);
    exit (-1);
  }
  if (proc == 0) {
    from = num_proc -1;
    if (MPI_Recv((void *) &flag, 1, MPI_INT, from, type, communicator, &st)
        != 0) {
      fprintf(stderr, "%s: ERROR on node %d\n", yo, proc);
      fprintf(stderr, "MPI_Recv failed, message type %d/n", type);
      exit (-1);
    }
  }

  /*
   * Do a final sync among all the processors, so that all of the other
   * processors must wait for Proc 0 to receive the final message from Proc
   * (num_proc-1)
   */

  MPI_Barrier(communicator);
}
