/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * Zoltan is distributed under the GNU Lesser General Public License 2.1.    * 
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __ALL_ALLO_H
#define __ALL_ALLO_H

#include "lb_const.h"

#ifndef HAVE_PROTOTYPES
#   if defined(__STDC__) || defined(__GNUC__) || defined(__cplusplus) || defined(c_plusplus)
#       define	HAVE_PROTOTYPES
#   endif
#endif

#undef PROTO
#ifdef HAVE_PROTOTYPES
#   define	PROTO(x)	x
#else
#   define	PROTO(x)	()
#endif

#define LB_MALLOC(a) LB_Malloc((a), __FILE__, __LINE__)
#define LB_REALLOC(a, b) LB_Realloc((a), (b), __FILE__, __LINE__)
#define LB_FREE(a) LB_Free((void **) (a), __FILE__, __LINE__)

/* function declarations for dynamic array allocation */

#ifdef __STDC__
extern double *LB_Array_Alloc(char *file, int lineno, int numdim, ...);
#else
extern double *LB_Array_Alloc();
#endif

extern int LB_Set_Malloc_Param(char *, char *);
extern void LB_Free(void **ptr, char *file, int lineno);
extern double *LB_Malloc(int n, char *file, int lineno);
extern double *LB_Realloc(void *ptr, int n, char *filename, int lineno);
extern void LB_Memory_Stats();
extern int LB_Memory_Num();
extern void LB_Free_Structure(LB *);

#ifdef __STDC__
extern void LB_Multifree(int n, ...);
#else
extern void LB_Multifree();
#endif

/* function prototypes for Fortran allocation functions */

#ifdef PGI
typedef void LB_FORT_MALLOC_INT_FN(int *arg, int *size, int **ret, int *hidden);
typedef void LB_FORT_MALLOC_GID_FN(LB_GID *arg, int *size, int **ret, int *hidden);
typedef void LB_FORT_MALLOC_LID_FN(LB_LID *arg, int *size, int **ret, int *hidden);
typedef void LB_FORT_FREE_INT_FN(int *arg, int *hidden);
typedef void LB_FORT_FREE_GID_FN(LB_GID *arg, int *hidden);
typedef void LB_FORT_FREE_LID_FN(LB_LID *arg, int *hidden);
#else
#ifdef FUJITSU
typedef void LB_FORT_MALLOC_INT_FN(int *arg, int *size, int **ret, int *hidden, int hidden2, int hidden3);
typedef void LB_FORT_MALLOC_GID_FN(LB_GID *arg, int *size, int **ret, int *hidden, int hidden2, int hidden3);
typedef void LB_FORT_MALLOC_LID_FN(LB_LID *arg, int *size, int **ret, int *hidden, int hidden2, int hidden3);
typedef void LB_FORT_FREE_INT_FN(int *arg, int *hidden);
typedef void LB_FORT_FREE_GID_FN(LB_GID *arg, int *hidden);
typedef void LB_FORT_FREE_LID_FN(LB_LID *arg, int *hidden);
#else
typedef void LB_FORT_MALLOC_INT_FN(int *arg, int *size, int **ret);
typedef void LB_FORT_MALLOC_GID_FN(LB_GID *arg, int *size, int **ret);
typedef void LB_FORT_MALLOC_LID_FN(LB_LID *arg, int *size, int **ret);
typedef void LB_FORT_FREE_INT_FN(int *arg);
typedef void LB_FORT_FREE_GID_FN(LB_GID *arg);
typedef void LB_FORT_FREE_LID_FN(LB_LID *arg);
#endif
#endif

/* type selector for LB_Special_Malloc */

enum LB_Special_Malloc_Type {
  LB_SPECIAL_MALLOC_INT,
  LB_SPECIAL_MALLOC_GID,
  LB_SPECIAL_MALLOC_LID
};

typedef enum LB_Special_Malloc_Type LB_SPECIAL_MALLOC_TYPE;

/* function declarations for special malloc */

extern int LB_Special_Malloc(struct LB_Struct *lb, void **array, int size,
                      LB_SPECIAL_MALLOC_TYPE type);
extern int LB_Special_Free(struct LB_Struct *lb, void **array,
                      LB_SPECIAL_MALLOC_TYPE type);
extern void LB_Register_Fort_Malloc(LB_FORT_MALLOC_INT_FN *fort_malloc_int,
                             LB_FORT_MALLOC_GID_FN *fort_malloc_GID,
                             LB_FORT_MALLOC_LID_FN *fort_malloc_LID,
                             LB_FORT_FREE_INT_FN *fort_free_int,
                             LB_FORT_FREE_GID_FN *fort_free_GID,
                             LB_FORT_FREE_LID_FN *fort_free_LID);

#endif
