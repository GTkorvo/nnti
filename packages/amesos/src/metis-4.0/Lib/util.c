/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * util.c
 *
 * This function contains various utility routines
 *
 * Started 9/28/95
 * George
 *
 * $Id$
 */

#include <metis.h>

#define Memory_MAXBLOCKS 10000
int GKmalloc_count = 0, GKmalloc_highwater = 0 ;
int Memory_nblocks = 0 ;
void *Memory_pointer [Memory_MAXBLOCKS] ;
int Memory_size [Memory_MAXBLOCKS] ;
int Memory_highwater = 0 ;
int Memory_nmalloc = 0 ;

void *(*my_malloc) (size_t) ;		/* pointer to malloc */
void (*my_free) (void *) ;		/* pointer to free */

void Memory_free (void *ptr)
{
    /* find a block of memory in the list and set its size to zero.
     * Do so in reverse order to save time (assume memory is typically
     * malloc'ed and free'd in a stack order) */
    int i, found = -1 ;
    for (i = 0 ; i < Memory_nblocks ; i++)
    {
	if (ptr == Memory_pointer [i])
	{
	    Memory_size [i] = 0 ;
	    found = i ;
	    break ;
	}
    }
    /* printf ("free: %u ", ptr) ; */
    if (found == -1) errexit ("(not found)") ;
    /* printf ("\n") ; */
    /* now reduce the number of memory blocks if possible by assuming that
     * all 'd and free'd blocks are returned to the memory pool */
    Memory_nblocks-- ;
    Memory_pointer [i] = Memory_pointer [Memory_nblocks] ;
    Memory_size    [i] = Memory_size    [Memory_nblocks] ;
}

void Memory_picture (void)
{
    /* print a picture of the memory map */
    int i ;
    printf ("\n--------------------------------------------------\n") ;
    for (i = 0 ; i < Memory_nblocks ; i++)
    {
	printf ("%15u: ", Memory_pointer [i]) ;
	if (Memory_size [i] > 0)
	{
	    printf ("#") ;
	}
	else
	{
	    printf (".") ;
	}
	/* if (i % 50 == 49) */ printf ("\n") ;
    }
    printf ("\n") ;
}

/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void errexit(char *f_str,...)
{
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

  sprintf(out2, "Error! %s", out1);

  fprintf(stdout, out2);
  fflush(stdout);

  abort();
}



#ifndef DMALLOC
/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
int *imalloc(int n, char *msg)
{
  if (n == 0)
    return NULL;

  return (int *)GKmalloc(sizeof(int)*n, msg);
}


/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
idxtype *idxmalloc(int n, char *msg)
{
  if (n == 0)
    return NULL;

  return (idxtype *)GKmalloc(sizeof(idxtype)*n, msg);
}


/*************************************************************************
* The following function allocates an array of float 
**************************************************************************/
float *fmalloc(int n, char *msg)
{
  if (n == 0)
    return NULL;

  return (float *)GKmalloc(sizeof(float)*n, msg);
}


/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
int *ismalloc(int n, int ival, char *msg)
{
  if (n == 0)
    return NULL;

  return iset(n, ival, (int *)GKmalloc(sizeof(int)*n, msg));
}



/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
idxtype *idxsmalloc(int n, idxtype ival, char *msg)
{
  if (n == 0)
    return NULL;

  return idxset(n, ival, (idxtype *)GKmalloc(sizeof(idxtype)*n, msg));
}


/*************************************************************************
* This function is my wrapper around malloc
**************************************************************************/
void *GKmalloc(int nbytes, char *msg)
{
  void *ptr;

  int i, usage ;

  if (nbytes < 0)
  {
      errexit ("Hey, METIS is asking for negative block size %d\n", nbytes) ;
  }

  if (nbytes == 0)
    return NULL;

  /* printf ("%d: %s: ", nbytes, msg) ; */
  /*
  ptr = (void *)my_malloc(nbytes);
  */
  ptr = (void *)malloc(nbytes);

  /* printf ("%u block %d\n", ptr, Memory_nblocks) ;*/
  if (ptr == NULL) 
    errexit("***Memory allocation failed for %s. Requested size: %d bytes", msg, nbytes);


  GKmalloc_count++ ;
/*  printf ("GKmalloc %15d (%d) %s  malloc_count %d\n", ptr, nbytes, msg, GKmalloc_count) ;*/
  if (GKmalloc_count > GKmalloc_highwater)
  {
      GKmalloc_highwater = GKmalloc_count ;
  }

  /* put the pointer into the memory table */
  if (Memory_nblocks > Memory_MAXBLOCKS)
  {
      errexit ("Out of table space!\n") ;
  }
  Memory_pointer [Memory_nblocks] = ptr ;
  Memory_size    [Memory_nblocks] = nbytes ;
  Memory_nblocks++ ;

  Memory_nmalloc++ ;

  usage = 0 ;
  for (i = 0 ; i < Memory_nblocks ; i++)
  {
      usage += Memory_size [i] ;
  }

  if (usage > Memory_highwater)
  {
      Memory_highwater = usage ;
  }

  /* Memory_picture ( ) ; */

  return ptr;
}
#endif



/*************************************************************************
* This function is my wrapper around free, allows multiple pointers    
**************************************************************************/
void GKfree(void **ptr1,...)
{
  va_list plist;
  void **ptr;

  if (*ptr1 != NULL)
  {
    GKmalloc_count-- ;
    /* printf ("GKfree   %15d malloc count %d\n", *ptr1, GKmalloc_count) ; */
    Memory_free (*ptr1) ;
    /*
    my_free(*ptr1);
    */
    free(*ptr1);
  }
  *ptr1 = NULL;

  va_start(plist, ptr1);

  /* while ((int)(ptr = va_arg(plist, void **)) != -1) { */
  while ((ptr = va_arg(plist, void **)) != LTERM) {
    if (*ptr != NULL)
    {
      GKmalloc_count-- ;
      /* printf ("GKfree   %15d multiple, malloc_count %d \n", *ptr, GKmalloc_count) ; */
      Memory_free (*ptr) ;
      /*
      my_free(*ptr);
      */
      free(*ptr);
    }
    *ptr = NULL;
  }

  /* Memory_picture ( ) ; */

  va_end(plist);
}            


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
int *iset(int n, int val, int *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
idxtype *idxset(int n, idxtype val, idxtype *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
float *sset(int n, float val, float *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}



/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int iamax(int n, int *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}


/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int idxamax(int n, idxtype *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}

/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int idxamax_strd(int n, idxtype *x, int incx)
{
  int i, max=0;

  n *= incx;
  for (i=incx; i<n; i+=incx)
    max = (x[i] > x[max] ? i : max);

  return max/incx;
}



/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int samax(int n, float *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}

/*************************************************************************
* These functions return the index of the almost maximum element in a vector
**************************************************************************/
int samax2(int n, float *x)
{
  int i, max1, max2;

  if (x[0] > x[1]) {
    max1 = 0;
    max2 = 1;
  }
  else {
    max1 = 1;
    max2 = 0;
  }

  for (i=2; i<n; i++) {
    if (x[i] > x[max1]) {
      max2 = max1;
      max1 = i;
    }
    else if (x[i] > x[max2])
      max2 = i;
  }

  return max2;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int idxamin(int n, idxtype *x)
{
  int i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int samin(int n, float *x)
{
  int i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int idxsum(int n, idxtype *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int idxsum_strd(int n, idxtype *x, int incx)
{
  int i, sum = 0;

  for (i=0; i<n; i++, x+=incx) {
    sum += *x;
  }

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
void idxadd(int n, idxtype *x, idxtype *y)
{
  for (n--; n>=0; n--)
    y[n] += x[n];
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int charsum(int n, char *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int isum(int n, int *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
float ssum(int n, float *x)
{
  int i;
  float sum = 0.0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
float ssum_strd(int n, float *x, int incx)
{
  int i;
  float sum = 0.0;

  for (i=0; i<n; i++, x+=incx)
    sum += *x;

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
void sscale(int n, float alpha, float *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] *= alpha;
}


/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
float snorm2(int n, float *v)
{
  int i;
  float partial = 0;
 
  for (i = 0; i<n; i++)
    partial += v[i] * v[i];

  return sqrt(partial);
}



/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
float sdot(int n, float *x, float *y)
{
  int i;
  float partial = 0;
 
  for (i = 0; i<n; i++)
    partial += x[i] * y[i];

  return partial;
}


/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
void saxpy(int n, float alpha, float *x, int incx, float *y, int incy)
{
  int i;
 
  for (i=0; i<n; i++, x+=incx, y+=incy) 
    *y += alpha*(*x);
}




/*************************************************************************
* This file randomly permutes the contents of an array.
* flag == 0, don't initialize perm
* flag == 1, set p[i] = i 
**************************************************************************/
void RandomPermute(int n, idxtype *p, int flag)
{
  int i, u, v;
  idxtype tmp;

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  if (n <= 4)
    return;

  for (i=0; i<n; i+=16) {
    u = RandomInRangeFast(n-4);
    v = RandomInRangeFast(n-4);
    SWAP(p[v], p[u], tmp);
    SWAP(p[v+1], p[u+1], tmp);
    SWAP(p[v+2], p[u+2], tmp);
    SWAP(p[v+3], p[u+3], tmp);
  }
}



/*************************************************************************
* This function returns true if the a is a power of 2
**************************************************************************/
int ispow2(int a)
{
  for (; a%2 != 1; a = a>>1);
  return (a > 1 ? 0 : 1);
}


/*************************************************************************
* This function initializes the random number generator
**************************************************************************/
void InitRandom(int seed)
{
  if (seed == -1) {
#ifndef __VC__
    srand48(7654321L);  
#endif
    srand(4321);  
  }
  else {
#ifndef __VC__
    srand48(seed);  
#endif
    srand(seed);  
  }
}

/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
int METIS_log2(int a)
{
  int i;

  for (i=1; a > 1; i++, a = a>>1);
  return i-1;
}

