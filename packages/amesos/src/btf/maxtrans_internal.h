#ifndef _MAXTRANS_INTERNAL_H
#define _MAXTRANS_INTERNAL_H

/* ========================================================================== */
/* === maxtrans_internal include file ======================================= */
/* ========================================================================== */

/* Not to be included in any user program. */

/* ========================================================================== */
/* make sure debugging and printing is turned off */
#ifndef NDEBUG
#define NDEBUG
#endif
#ifndef NPRINT
#define NPRINT
#endif
/* To enable debugging and assertions, uncomment this line: 
#undef NDEBUG
*/
/* To enable diagnostic printing, uncomment this line: 
#undef NPRINT
*/
/* ========================================================================== */

#ifdef MATLAB_MEX_FILE
#include "matrix.h"
#include "mex.h"
#define ASSERT(a) mxAssert(a, "")
#else
#include <stdio.h>
#include <assert.h>
#define ASSERT(a) assert(a)
#endif

#undef TRUE
#undef FALSE
#undef PRINTF
#undef MIN

#ifndef NPRINT
#define PRINTF(s) { printf s ; } ;
#else
#define PRINTF(s)
#endif

#define TRUE 1
#define FALSE 0
#define EMPTY (-1)
#define MIN(a,b) (((a) < (b)) ?  (a) : (b))

#endif
