/*! \file cmm.h
    \brief Main include file for CMMLIB
    \ingroup base
     Author:  Mark P. Sears

     Copyright: See file COPYRIGHT in the CMMLIB source distribution.

*/

#ifndef CMM_H
#define CMM_H 1

// DEPENDENT INCLUDE FILES:
#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
extern "C"
{
#endif                          /* __cplusplus */

// CONTENTS:

// Basic types
typedef int Boolean;

#define True  1
#define False 0

/* Real type */
typedef double r8;

/* Complex type */
typedef struct
  {
   r8 r, i;
  } c8;

/*! \def CMM_BYTE
    type identifier for unsigned char
*/
#define CMM_BYTE                1

/*! \def CMM_CHAR
 type identifier for signed char
*/
#define CMM_CHAR                2

/*! \def CMM_UINT
 type identifier for unsigned integer
*/
#define CMM_UINT                3

/*! \def CMM_INTEGER
 type identifier for signed integer
*/
#define CMM_INTEGER             4

/*! \def CMM_SINGLE_REAL
 type identifier for single precision real
*/
#define CMM_SINGLE_REAL         5

/*! \def CMM_REAL
 type identifier for double precision real
*/
#define CMM_REAL                6

/*! \def CMM_SINGLE_COMPLEX
 type identifier for single precision complex
*/
#define CMM_SINGLE_COMPLEX      7     // double precision complex (8 bytes normally)

/*! \def CMM_COMPLEX
 type identifier for double precision complex
*/
#define CMM_COMPLEX             8     // double precision complex (16 bytes normally)

int   cmm_typesize(int type);

char *cmm_typename(int type);

// Basic macros (inline functions)

// There are still problems with these macros,
// for example x could evaluate a function with a side effect.
// Inline functions would be better, in fact I could use
// an ifdef.

#ifndef ABS
/*! \def ABS(x)
 Compute absolute value of x.
*/
#define ABS(x)   ((x) >= 0. ? (x) : (-(x)))
#endif

#ifndef MAX
/*! \def MAX(x)
 Compute maximum of x, y.
*/
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#endif

#ifndef MIN
/*! \def MIN(x)
 Compute minimum of x, y.
*/
#define MIN(x,y) ((x) <= (y) ? (x) : (y))
#endif

/*! \def DBLE(x)
  Convert x to double precision real.
*/
#define DBLE(x)   ((r8)(x))

/*! \def INT(x)
  Convert x to integer.
*/
#define INT(x)    ((int)(x))

/*! \def IFLOOR(x)
  Compute floor(x) as integer.
*/
#define IFLOOR(x) ((int)(floor(x)))

/*! \def ICEIL(x)
  Compute ceil(x) as integer.
*/
#define ICEIL(x)  ((int)(ceil(x)))

/*! \def FRAC(x)
  Compute fractional part of x.
*/
#define FRAC(x)   ((x) - floor(x))

/*! \def RND(x)
  Round x to nearest integer, return as integer.
*/
#define RND(x)    ((int)((x) + .5))

// Math constants

// ref Abramowitz and Stegun
/*! \def CMM_SQRT2 
  \f$ \sqrt{2} \f$.
*/
#define CMM_SQRT2         (1.414213562373095145)

/*! \def CMM_PI
  \f$ \pi \f$.
*/
#define CMM_PI            (3.141592653589793238462643)

/*! \def CMM_PI_OVER_2 
  \f$ {\pi}\over{2} \f$.
*/
#define CMM_PI_OVER_2     (.5 * CMM_PI)

/*! \def CMM_2PI 
  \f$ 2\pi \f$.
*/
#define CMM_2PI           (2. * CMM_PI)

/*! \def CMM_4PI 
  \f$ 4\pi \f$.
*/
#define CMM_4PI           (4. * CMM_PI)

/*! \def CMM_SQRTPI 
  \f$ \sqrt{\pi} \f$.
*/
#define CMM_SQRTPI        (1.7724538509055160272)

/*! \def CMM_E
  \f$ e \f$, base of natural logarithm.
*/
#define CMM_E             (2.7182818284590452353)

/*! \def CMM_EULER 
  \f$ \gamma \f$, Euler's constant.
*/
#define CMM_EULER         (.577215664901532860606)

/*! \def CMM_LN2 
  \f$ \ln{2} \f$.
*/
#define CMM_LN2           (.693147180559945309417)

/*! \def CMM_LN10 
  \f$ \ln(10) \f$.
*/
#define CMM_LN10          (2.3025850929940456840) 

// Machine constants. These should work with 64 bit IEEE double precision and 32 bit integers
/*! \def CMM_EPS
  \f$ \epsilon \f$, smallest double precision number such that \f$ 1+\epsilon \ne 1 \f$.
*/
#define CMM_EPS           (1.110223e-16)        /* smallest number x such that 1+x != 1 */

/*! \def CMM_MAXREAL
  Largest possible (magnitude) double precision number.
*/
#define CMM_MAXREAL       (1.e255)              // should really be 1.e305 or so

/*! \def CMM_MINREAL
  Smallest possible (magnitude) double precision number.
*/
#define CMM_MINREAL       (1.e-255)             // should really be 1.e-320 or so

/*! \def CMM_MAXEXPONENT
  Largest possible double precision exponent.
*/
#define CMM_MAXEXPONENT   ( 255)

/*! \def CMM_MINEXPONENT
  Smallest possible double precision exponent.
*/
#define CMM_MINEXPONENT   (-255)

/*! \def CMM_MAXINT
  Largest possible positive signed integer.
*/
#define CMM_MAXINT        (2147483647)

/*! \def CMM_MININT
  Largest possible negative signed integer.
*/
#define CMM_MININT        (-2147483647)

// Physics constants and conversion factors

#include "phys.h"

// Conversion and string routines
r8    cmm_str_cvt_r8(char *s);
int   cmm_str_cvt_int(char *s);
int   cmm_str_cvt_boolean(char *s);
int   cmm_str_cvt_choice(char *str, int nchoices, char **choices);
char *cmm_str_clone(char *s);
int   cmm_str_copy(char *src, char *dst);
int   cmm_str_copyn(char *src, char *dst, int n);
int   cmm_str_copyd(char *src, char *dst, char *delims);
int   cmm_str_len(char *src);
int   cmm_str_trim_eol(char *src);
Boolean cmm_str_compare(char *s, char *t);
Boolean cmm_str_compareto(char *s, char *t, int maxlen);
void  cmm_str_toupper(char *s);
void  cmm_str_tolower(char *s);

// Basic parsing routines
int   cmm_str_arglist(char *str, int maxargs, int maxbuf, char **args, char *buf);
void  cmm_str_encode_qs(char *src, char *dst, int maxlen);
int   cmm_str_parse_file(char *name, char *path, char *file, char *ext);

// Memory copy
void  cmm_mem_copy(int l, void *src, void *dst);

// Warning/error/exit 
void  cmm_exit(int status);
void  cmm_print(char *fmt, ...);
void  cmm_warn(char *fmt, ...);
void  cmm_fatal(char *fmt, ...);
void  cmm_pause(char *fmt, ...);

// cpu clock, date
r8    cmm_clock(void);
char *cmm_date(void);

// memory allocation
void *cmm_alloc(long l);
r8   *cmm_r8_alloc(int n);
c8   *cmm_c8_alloc(int n);
int  *cmm_int_alloc(int n);
void  cmm_free(void *p);
void *cmm_realloc(void *oldp, long l);

void  cmm_preallocate(int blocks);
void  cmm_limitheap(int nlim);
void  cmm_memstat(void *base);
int   cmm_memcheck(void *base);

// miscellaneous math
Boolean cmm_nbrhood(r8 x, r8 y, r8 eps);  // return true if   |x-y| <= eps.
r8    cmm_sq(r8 x);                      // return x*x
r8    cmm_ipow(r8 x, int n);             // return x^n
r8    cmm_ipow10(int n);                 // return 10.^n
r8    cmm_factorial(int n);              // return n!
r8    cmm_combination(int n, int m);     // compute combination n!/(m!(n-m)!)

// END CONTENTS

#ifdef __cplusplus
}
#endif                          /* __cplusplus */

#endif
