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

#ifndef __SHARED_CONST_H
#define __SHARED_CONST_H

/* Definitions shared by parallel RCB and RIB */
#define DEFAULT_CHECK_GEOM 1
#define TINY 1.0e-6

/* Data structures shared by parallel RCB and RIB */
/* dot to balance on for RCB and RIB */ 
struct Dot_Struct {	        /* dot = point in 3-space */
  double    X[3];		/* location of dot */
  double    Weight;             /* weight of dot - if used must be > 0 */
  int Proc;                     /* Processor ID for processor owning a dot.
                                   For now, we'll keep it with a dot, even 
                                   though the global and local ids are now
                                   stored separately.                       */
};

extern int LB_RB_Build_Structure(LB *, LB_ID_PTR *, LB_ID_PTR *, 
  struct Dot_Struct **, int *, int *, int *, int, int);

extern void LB_RB_Print_All(LB *, LB_ID_PTR , struct Dot_Struct *,
  int , int , int , LB_ID_PTR , int *);

extern int LB_RB_Send_Outgoing(LB *, LB_ID_PTR *, LB_ID_PTR *,
  struct Dot_Struct **, int *, int *, int *, int *, int, int *, double, int,
  int *, int, MPI_Comm, int, int);

extern int LB_RB_Send_Dots(LB *, LB_ID_PTR *, LB_ID_PTR *,
  struct Dot_Struct **, int *, int *, int, int *, int *, int, int *, double,
  int, int *, int, MPI_Comm);

extern int LB_RB_Return_Arguments(LB *, LB_ID_PTR, LB_ID_PTR, 
  struct Dot_Struct *, int, LB_ID_PTR *, LB_ID_PTR *, int **, int);

extern int LB_RB_check_geom_input(LB *, struct Dot_Struct *, int);

extern int LB_RB_check_geom_output(LB *, struct Dot_Struct *, 
  int, int, void *);

extern void LB_RB_stats(LB *, double, struct Dot_Struct *, int , double *, 
  int *, int, int *, void *, int);

extern int LB_Use_IDs(LB *);

#endif
