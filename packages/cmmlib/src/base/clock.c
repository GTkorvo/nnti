/*! \file clock.c
    \brief Clock and date functions.
    \ingroup base

    These functions are obviously very central, and unfortunately
    differ extensively among os.

*/

#include "cmm.h"

#ifndef DOXYGEN_SKIP
#include "time.h"
#endif

/*!
    \brief Return elapsed CPU time.
*/
#if defined(USE_DCLOCK)
double dclock(void);
r8 cmm_clock(void)
   {
    return dclock();
   }

#elif defined(USE_DCLOCK_)
double dclock_();
r8 cmm_clock(void)
   {
    return dclock_();
   }
#elif defined(USE_MPI)
double mpi_wtime();
r8 cmm_clock(void)
   {
    return mpi_wtime();
   }
#else
// Use a standard POSIX routine.

#ifndef DOXYGEN_SKIP
#ifdef linux
#include  <linux/param.h>
#endif
#endif

r8 cmm_clock(void)
{
  double ct;
  clock_t t0;

  t0 = clock();
  ct = (double) t0;
  ct = ct / (double) CLOCKS_PER_SEC;

  return ct;
}
#endif

#ifndef DOXYGEN_SKIP
#include <sys/time.h>
#include <unistd.h>
#endif

/*!
   \brief Return current date and time.
*/
char *cmm_date(void)
{
  static char datebuf[256];
  char *s;
  time_t the_time;
  time(&the_time);
  s = asctime(localtime(&the_time));
  cmm_str_copyn(s, datebuf, 256);
  cmm_str_trim_eol(datebuf);
  return datebuf;
}
