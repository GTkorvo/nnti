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
 * $Name$
 *====================================================================*/
#ifndef lint
static char *cvs_all_allo_c =
  "$Id$";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "mpi.h"
#include "lb_const.h"
#include "all_allo_const.h"
#include "params_const.h"

#ifdef __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif

static int DEBUG_MEMORY = 0;	/* Flag for detecting memory leaks */
static int bytes_used = 0;	/* Sum of active allocations */
static int bytes_max = 0;	/* Largest total of active allocations */

static int nmalloc = 0;         /* number of calls to malloc */
static int nfree = 0;           /* number of calls to free */

static struct malloc_debug_data {
  int       order;			/* which malloc call is it? */
  int       size;			/* size of malloc invocation */
  double   *ptr;			/* memory location returned */
  struct malloc_debug_data *next;	/* pointer to next element */
} *top = NULL;



int LB_Set_Malloc_Param(
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    int status;
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */
    PARAM_VARS malloc_params[] = {
	{ "DEBUG_MEMORY", &DEBUG_MEMORY, "INT" },
	{ NULL, NULL, NULL }
    };

    status = LB_Check_Param(name, val, malloc_params, &result, &index);
    if (status == 0 && index == 0) {
	DEBUG_MEMORY = result.ival;
	status = 3;
    }

    return(status);
}


/******************************************************************************
 *
 *                    Dynamic Allocation of Multidimensional Arrays
 *-----------------------------------------------------------------------------
 *
 * Example Usage:
 *
 *     typedef  struct
 *       {      int     bus1;
 *              int     bus2;
 *              int     dest;
 *      }       POINT;
 *
 *      POINT    **points, corner;
 *
 *      points = (POINT **)LB_Array_Alloc(file, lineno, 2, x, y, sizeof(POINT));
 *                                        ^     ^       ^  ^  ^
 *                                        |     |       |  |  |
 *                 name of calling file---*     |       |  |  |
 *                                              |       |  |  |
 *                  line number of call---------*       |  |  |
 *                                                      |  |  |
 *                 number of dimensions-----------------+  |  |
 *                                                         |  |
 *                  first dimension max--------------------+  |
 *                                                            |
 *                 second dimension max-----------------------+
 *
 *         (points may be now be used as if it were declared
 *          POINT points[x][y])
 *
 *          This particular version is limited to dimensions of 4 or less.
 *
 *      corner = points[2][3]; (refer to the structure as you would any array)
 *
 *      free (points); (frees the entire structure in one fell swoop)
 *
 *****************************************************************************/
/******************************************************************************
*       The following section is a commented section containing
*       an example main code:
*******************************************************************************
*double *LB_Array_Alloc();
*main()
*{
*  int ***temp;
*   int *temp2;
*   int i, j, k;
*   int il, jl, kl;
*
*   malloc_debug(2);
*   il = 2;
*   jl = 3;
*   kl = 3;
*   temp = (int ***) Array_Alloc(__FILE__,__LINE__,3,il,jl,kl,sizeof(int));
*   for (i=0; i<il; i++) {
*      for (j=0; j<jl; j++) {
*         for (k=0; k<kl; k++) temp[i][j][k] = 1;
*      }
*   }
*
*   temp2 = (int *) malloc(10*sizeof(int));
*   for (i=0; i<10; i++) temp2[i] = 0;
*
*   for (i=0; i<il; i++) {
*      for (j=0; j<jl; j++) {
*         for (k=0; k<kl; k++) (void) printf(" %d\n", temp[i][j][k]);
*      }
*   }
*   malloc_verify();
*}
******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef __STDC__

double *LB_Array_Alloc(char *file, int lineno, int numdim, ...)

#else

double *LB_Array_Alloc(va_alist)
va_dcl

#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

{
  char *yo = "LB_Array_Alloc";
  int i, j;
  struct dimension {
    long index;  /* Number of elements in the dimension  */
    long total;  /* Total number of elements             */
    long size;   /* Size of a single element in bytes    */
    long off;    /* offset from beginning of array       */
  } dim[4];      /* Info about each dimension            */

#ifndef __STDC__
  char *file;           /* Filename of source code from call.   */
  int lineno;           /* Line number of call.                 */
  int numdim;           /* Number of dimensions                 */
#endif

  long    total;        /* Total size of the array              */
  double *dfield;       /* ptr to avoid lint complaints         */
  char   *field;        /* The multi-dimensional array          */
  char  **ptr;          /* Pointer offset                       */
  char   *data;         /* Data offset                          */
  va_list va;           /* Current pointer in the argument list */

#ifdef __STDC__
  va_start(va, numdim);
#else
  va_start(va);
  file = va_arg(va, char *);
  lineno = va_arg(va, int);
  numdim = va_arg(va, int);
#endif

  if (numdim <= 0) {
    fprintf(stderr, "%s (%s: %d) ERROR: number of dimensions, %d, is <=0\n",
            yo, file, lineno, numdim);
    return((double *) NULL);
  }
  else if (numdim > 4) {
    fprintf(stderr, "%s (%s: %d) ERROR: number of dimensions, %d, is > 4\n",
            yo, file, lineno, numdim);
    return((double *) NULL);
  }

  dim[0].index = va_arg(va, int);

  if (dim[0].index <= 0) {
#ifdef DEBUG
    fprintf(stderr, "WARNING, %s (%s: %d) called with first "
            "dimension <= 0, %d; will return NULL\n", 
            yo, file, lineno, dim[0].index);
#endif
    return((double *) NULL);
  }

  dim[0].total = dim[0].index;
  dim[0].size  = sizeof(void *);
  dim[0].off   = 0;
  for (i = 1; i < numdim; i++) {
    dim[i].index = va_arg(va, int);
    if (dim[i].index <= 0) {
      fprintf(stderr, "WARNING: %s (%s: %d) called with dimension %d <= 0, "
              "%d; will return NULL\n",
              yo, file, lineno, i+1, dim[i].index);
      return((double *) NULL);
    }
    dim[i].total = dim[i-1].total * dim[i].index;
    dim[i].size  = sizeof(void *);
    dim[i].off   = dim[i-1].off + dim[i-1].total * dim[i-1].size;
  }

  dim[numdim-1].size = va_arg(va, int);
  va_end(va);

  /* Round up the last offset value so data is properly aligned. */

  dim[numdim-1].off = dim[numdim-1].size *
    ((dim[numdim-1].off+dim[numdim-1].size-1)/dim[numdim-1].size);

  total = dim[numdim-1].off + dim[numdim-1].total * dim[numdim-1].size;

  dfield = (double *) LB_Malloc((int) total, file, lineno); 

  if (dfield != NULL) {
    field  = (char *) dfield;
    for (i = 0; i < numdim - 1; i++) {
      ptr  = (char **) (field + dim[i].off);
      data = (char *) (field + dim[i+1].off);
      for (j = 0; j < dim[i].total; j++) {
        ptr[j] = data + j * dim[i+1].size * dim[i+1].index;
      }
    }
  }

  return dfield;

} /* LB_Array_Alloc */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/


/* Safe version of malloc.  Does not initialize memory .*/

double *LB_Malloc(int n, char *filename, int lineno)

{
  char *yo = "LB_Malloc";
  struct malloc_debug_data *new;     /* data structure for malloc data */
  int       proc;             /* processor ID for debugging msg */
  double *pntr;           /* return value */

  if (n > 0) {
    nmalloc++;
    pntr = (double *) malloc((unsigned) n);
    if (pntr == NULL) {
      MPI_Comm_rank(MPI_COMM_WORLD, &proc);
      fprintf(stderr, "%s (from %s,%d) No space on proc %d - number of bytes "
              "requested = %d\n", yo, filename, lineno, proc, n);
      return ((double *) NULL);
    }
  }
  else if (n == 0)
    pntr = NULL;
  else {		/* n < 0 */
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    fprintf(stderr, "%s (from %s,%d) ERROR on proc %d: "
	    "Negative malloc argument. (%d)\n", yo, filename, lineno, proc, n);
    return ((double *) NULL);
  }

  if (DEBUG_MEMORY > 1 && pntr != NULL) {
    /* Add allocation to linked list of active allocations. */
    new = (struct malloc_debug_data *)
      malloc(sizeof(struct malloc_debug_data));

    if (new == NULL) {
      MPI_Comm_rank(MPI_COMM_WORLD, &proc);
      fprintf(stderr, "WARNING: No space on proc %d for malloc_debug %d.\n",
	proc, n);
      return (pntr);
    }

    new->order = nmalloc;
    new->size = n;
    new->ptr = pntr;
    new->next = top;
    top = new;
    bytes_used += n;
    if (bytes_used > bytes_max) {
      bytes_max = bytes_used;
    }
  }

  if (DEBUG_MEMORY > 2) {
    /* Print out details of allocation. */
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    fprintf(stderr, "Proc %d: order=%d, size=%d, location=0x%lx\n",
      proc, nmalloc, n, (long) pntr);
  }


  return pntr;

} /* LB_Malloc */

/* Safe version of realloc. Does not initialize memory. */

double *LB_Realloc(void *ptr, int n, char *filename, int lineno)
{
  char *yo = "LB_Realloc";
  struct malloc_debug_data *dbptr;   /* loops through debug list */
  int       proc;             /* processor ID */
  double   *p;                /* returned pointer */

  if (ptr == NULL) {	/* Previous allocation not important */
    if (n == 0) {
      p = NULL;
    }
    else {
      p = LB_Malloc(n, filename, lineno);
    }
  }
  else {
    if (n == 0) {
      LB_Free((void **) &ptr);
      p = NULL;
    }
    else {
      p = (double *) realloc((char *) ptr, n);

      if (DEBUG_MEMORY > 1) {
        /* Need to replace item in allocation list */
        for (dbptr = top; dbptr != NULL && (char *) dbptr->ptr != ptr;
	   dbptr = dbptr->next);
	if (dbptr == NULL) {	/* previous allocation not found in list. */
	   MPI_Comm_rank(MPI_COMM_WORLD, &proc);
	   fprintf(stderr, "Proc %d: Memory error: "
	     "In realloc, address not found in debug list (0x%lx)\n",
	     proc, (long) ptr);
	}
	else {	/* Update entry in allocation list */
	  bytes_used += n - dbptr->size;
	  dbptr->size = n;
	  dbptr->ptr = p;
	  if (bytes_used > bytes_max) {
	    bytes_max = bytes_used;
	  }
	}
      }

      if (p == NULL) {
        MPI_Comm_rank(MPI_COMM_WORLD, &proc);
        fprintf(stderr, "%s (from %s,%d) No space on proc %d - "
		"number of bytes requested = %d\n",
		yo, filename, lineno, proc, n);
      }
    }
  }

  return (p);
} /* LB_Realloc */


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_Free (void **ptr)
{
  struct malloc_debug_data *dbptr;   /* loops through debug list */
  struct malloc_debug_data **prev;   /* holds previous pointer */
  int       proc;             /* processor ID */

/*
 *  This version of free calls the system's free function
 *  with maximum error checking. It also doesn't call free if ptr is
 *  the NULL pointer.
 */
 
  if (ptr == NULL || *ptr == NULL) 
    return;

  nfree++;

  if (DEBUG_MEMORY > 1) {
    /* Remove allocation of list of active allocations */
    prev = &top;
    for (dbptr = top; dbptr != NULL && (char *) dbptr->ptr != *ptr;
      dbptr = dbptr->next) {
      prev = &(dbptr->next);
    }
    if (dbptr == NULL) {
      MPI_Comm_rank(MPI_COMM_WORLD, &proc);
      fprintf(stderr, "Proc %d: Memory error: In free, address (0x%lx) "
	"not found in debug list\n", proc, (long) *ptr);
   }
   else {
       *prev = dbptr->next;
       bytes_used = -dbptr->size;
       free((char *) dbptr);
       }
   }
 
  free(*ptr);
 
  /* Set value of ptr to NULL, to flag further references to it. */
  *ptr = NULL;

}  /* LB_Free */


/* Print out status of malloc/free calls.  Flag any memory leaks. */

void      LB_Memory_Stats()
{
    struct malloc_debug_data *dbptr;	/* loops through debug list */
    int       proc;		/* processor ID */

    MPI_Comm_rank(MPI_COMM_WORLD, &proc);

    if (DEBUG_MEMORY == 1) {
	fprintf(stderr, "Proc %d: Calls to malloc = %d,  Calls to free = %d\n",
		proc, nmalloc, nfree);
    }
    if (DEBUG_MEMORY > 1) {
	fprintf(stderr, "Proc %d: Calls to malloc = %d,  Calls to free = %d, maximum bytes = %d\n",
		proc, nmalloc, nfree, bytes_max);
	if (top != NULL) {
	    fprintf(stderr, "Proc %d: Remaining allocations:\n", proc);
	    for (dbptr = top; dbptr != NULL; dbptr = dbptr->next) {
		fprintf(stderr, " order=%d, size=%d, location=0x%lx\n", dbptr->order,
			dbptr->size, (long) dbptr->ptr);
	    }
	}
    }
} /* LB_Memory_Stats */


/* Return number associated with the most recent malloc call. */
int       LB_Malloc_Num()
{
  return (nmalloc);
} /* LB_Memory_Num */


/*****************************************************************************/
/*                      END of all_allo.c                                     */
/*****************************************************************************/








