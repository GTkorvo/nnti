/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */

/* ******************************************************************** */
/* Functions for the BICGSTABL Krylov solver                            */
/* ******************************************************************** */
/* Author        : Charles Tong (LLNL)                                  */
/* Date          : December, 1999                                       */
/* ******************************************************************** */

#include "ml_krylov.h"

#ifndef __MLCGSTABL__
#define __MLCGSTABL__

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

extern int ML_BICGSTABL_Solve(ML_Krylov *,int,double *rhs,double *sol);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif

