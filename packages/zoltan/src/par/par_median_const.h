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

#ifndef __PAR_UTIL_CONST_H
#define __PAR_UTIL_CONST_H

#include <mpi.h>

extern int Zoltan_RB_find_median(int, double *, double *, int *,
  int, int, double, MPI_Comm,
  double *, int, int *, int, int,
  int, int, double, double, double,
  double *, double *, int *, int);

#endif
