/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 *====================================================================*/

#ifndef PARAMS_CONST_H
#define PARAMS_CONST_H

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

/*
 * Type used to pass list of available parameters to lower routines.
 */

typedef struct Param_Vars {
  char *name;			/* Parameter variable name (all CAPS) */
  void *ptr;			/* Pointer to parameter variable */
  char *type;			/* type of parameter: */
				/* INT, DOUBLE, LONG, STRING, or CHAR */
} PARAM_VARS;

/* string length limit for param val. Allocate space of this + 1 */
#define MAX_PARAM_STRING_LEN 50

typedef struct Param_Utype {
  int ival;
  double dval;
  long lval;
  char sval[MAX_PARAM_STRING_LEN + 1];
  char cval;
} PARAM_UTYPE;

/* API for general parameter setting functions */
typedef int LB_SET_PARAM_FN(char *, char *); 

/* function declarations for parameter modification routines */

extern int LB_Assign_Param_Vals(LB_PARAM *, PARAM_VARS *, int, int);
extern int LB_Bind_Param(PARAM_VARS *, char *, void *);
extern int LB_Set_Param(LB *, char *, char *);
extern int LB_Set_Key_Param(LB *, char *, char *);
extern void LB_Print_Key_Params(LB *);
extern int LB_Check_Param(char *, char *, PARAM_VARS *,
    PARAM_UTYPE *, int *);
extern void LB_Free_Params(LB_PARAM **);

#endif
