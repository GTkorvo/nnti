/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_rand.h"


/****************************************************************************/
/* Random Number generator due to Knuth found in Numerical Recipes in C
 * (2nd edition) by Press, Vetterling, Teukolsky, Flannery (Page 284.)
 * Needed because different random number implementations on different 
 * machines produced different answers!  This generator provides a portable, 
 * fast, algorithm with adequate random number generation. 
 * This generator was designed for 32 bit ints but works for 64, too.
 */

static unsigned int zidum = ZOLTAN_RAND_INIT;

unsigned int Zoltan_Seed()
{
/* Function that returns the current value of the Zoltan seed. */
  return zidum;
}


unsigned Zoltan_Rand(unsigned int *myidum) {
/* 
 * If myidum is non-NULL, use *myidum as the generator value.  This feature
 * allows synchronization of the RNG across processors.
 * If myidum is NULL, use zidum.
 */
unsigned int *idum;

  if (myidum) 
    idum = myidum;
  else
    idum = &zidum;
  *idum = ((1664525U * *idum) + 1013904223U) % ZOLTAN_RAND_MAX; /* mod is slow! */
  return (*idum);
}



void Zoltan_Srand (unsigned int seed, unsigned int *myidum) {
/* 
 * If myidum is non-NULL, set *myidum to the seed.  
 * If myidum is NULL, set zidum.
 */
unsigned int *idum;

  if (myidum) 
    idum = myidum;
  else
    idum = &zidum;
  *idum = seed;
}



void Zoltan_Srand_Sync(
  unsigned int seed, 
  unsigned int *myidum,
  MPI_Comm comm
)
{
/* Synchronize the random number seed across processor within a communicator.
 * Proc 0's seed is the broadcast seed used by all procs. 
 * If myidum is non-NULL, set *myidum to the seed.  
 * If myidum is NULL, set zidum.
 */

unsigned int *idum;

  if (myidum) 
    idum = myidum;
  else
    idum = &zidum;
  *idum = seed;
  MPI_Bcast(idum, 1, MPI_UNSIGNED, 0, comm);
}


unsigned int Zoltan_Rand_InRange (unsigned int *myidum, unsigned int n)
{
  double denom = ZOLTAN_RAND_MAX + 1.0;

  return (int) ((double) n * (double) Zoltan_Rand(myidum) / denom);
}
    
/* Randomly permute an array of ints. */
void Zoltan_Rand_Perm_Int (int *data, int n, unsigned int *myidum)
{
int i, number, temp;
double denom = ZOLTAN_RAND_MAX + 1.0;
/* Scaling of random number to appropriate range is done as recommended
 * in Numerical Recipes in C.
 */

  for (i = n; i > 0; i--) {
    number       = (int) ((double) i * (double) Zoltan_Rand(myidum) / denom);
    temp         = data[number];
    data[number] = data[i-1];
    data[i-1]    = temp;
  }
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
