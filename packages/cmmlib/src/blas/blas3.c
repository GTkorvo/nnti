/*! \file blas3.c
    \brief  Level 3 BLAS interface.
    \ingroup blas

    These routines define a C interface to level 3 BLAS routines. If
    possible, we use an external BLAS library, otherwise we use internal
    routines.

    The BLAS standard is unfortunately a Fortran standard, which means
    it is not well defined from our point of view. The purpose of these
    routines is to cover up this ugliness. That is why these routines
    themselves have so much conditional compilation garbage (which I
    abhor). We use the CBLAS interface if available, a Fortran interface
    if we can. The backup is to use an internal reference code.

    The actual part of the BLAS that is used is very narrow. Essentially
    we only call dgemm and zgemm.



   Supported forms:

     C  = fc * C + fa * A*B       -- c8_gemm_XX  (a near-traditional GEMM interface)
     C  = A*B                     -- c8_mm_XX
     C += A*B                     -- c8_mm_a_XX
     C -= A*B                     -- c8_mm_s_XX
     C  = A*C                     -- c8_tm
     C  = C*A                     -- c8_mt

   where XX is nn, hn, or nh (this list could be extended).

   Current implementation is that most of these map back
   to c8_gemm_XX. Note that c8_mt requires allocation
   of storage; this will come off the stack if possible
   and is allocated/freed if not.

   Also included:

     C  = f*A                     -- c8_m_scale
     C += f*A                     -- c8_m_axpy
     C  = fa*A + fb*B             -- c8_m_axpy2
     C  = fa*A + fc*C             -- c8_m_axpby

   A and B are not allowed to overlap C.

*/

#define CMM_SHORTHAND
#include "cmm.h"
#include "cmmblas.h"

/*  Only one of the following should be defined.

    USE_CBLAS
    USE_FBLAS
    USE_MKL_CBLAS
    USE_MYZGEMM

*/

#ifdef USE_CBLAS
#include "cblas.h"
#define CMM_HAVE_BLAS
#endif

#ifdef USE_MKL_CBLAS
#include "mkl_cblas.h"
#define USE_CBLAS
#define CMM_HAVE_BLAS
#endif

#ifdef USE_FBLAS

/* 
  If CBLAS is not defined, then define a more or less 
  standard interface to call the Fortran BLAS routines. This is 
  used by CXML BLAS and GOTO BLAS, which do not have a CBLAS 
  interface.  

  If none of these work, complain! (NOT TO ME!) complain to your
  Fortran/BLAS vendor. CBLAS is the standard interface.
*/

#ifdef FORTRAN_MANGLE_NONE
#define DGEMM_FTN dgemm
#define ZGEMM_FTN zgemm
#endif

#ifdef FORTRAN_MANGLE_ 
#define DGEMM_FTN dgemm_
#define ZGEMM_FTN zgemm_
#endif

#ifdef FORTRAN_MANGLE_UPPER
#define DGEMM_FTN DGEMM
#define ZGEMM_FTN ZGEMM
#endif

#ifdef FORTRAN_MANGLE_UPPER_
#define DGEMM_FTN DGEMM_
#define ZGEMM_FTN ZGEMM_
#endif

// default
#ifndef DGEMM_FTN
#define DGEMM_FTN dgemm_
#define ZGEMM_FTN zgemm_
#endif

#define DGEMM_NORMAL              "N"
#define DGEMM_TRANSPOSE           "T"

   void DGEMM_FTN(char *ta, char *tb, int *m, int *n, int *k,
         r8 *alpha, r8 *a, int *lda,
                    r8 *b, int *ldb,
         r8 *beta,  r8 *c, int *ldc);


#define ZGEMM_NORMAL              "N"
#define ZGEMM_TRANSPOSE           "T"
#define ZGEMM_CONJUGATE_TRANSPOSE "C"

   void ZGEMM_FTN(char *ta, char *tb, int *m, int *n, int *k,
         c8 *alpha, c8 *a, int *lda,
                    c8 *b, int *ldb,
         c8 *beta,  c8 *c, int *ldc);

#define CMM_HAVE_BLAS
#endif

/*
    If we DONT have any blas, use a default.
*/

#ifndef CMM_HAVE_BLAS
#include "myzgemm.h"
#define CMM_HAVE_BLAS
#endif

static r8 d_one  =  1.;
static r8 d_mone = -1.;
static r8 d_zero =  0.;

static c8 z_one  = { 1., 0.};
static c8 z_mone = {-1., 0.};
static c8 z_zero = { 0., 0.};

// C = fc*C + fa*A*B
void cmm_c8_gemm_nn(c8 fc, c8 fa,
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
#ifdef USE_CBLAS
//   cmm_warn("calling cblas_zgemm %d %d %d\n",cm,cn,an);
   cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
        cm, cn, an,
        &fa, a, as, b, bs, &fc, c, cs);
#elif USE_MYZGEMM
   myzgemm_nn(fc, fa,
           a, as, am, an,
           b, bs, bm, bn,
           c, cs, cm, cn);
#else // Fortran BLAS
   ZGEMM_FTN(ZGEMM_NORMAL, ZGEMM_NORMAL,
        &cm, &cn, &an,
        &fa, a, &as, b, &bs, &fc, c, &cs);
#endif
  }

// C = fc*C + fa*A'*B
void cmm_c8_gemm_hn(c8 fc, c8 fa,
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
#ifdef USE_CBLAS
   cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans,
        cm, cn, am, &fa, a, as, b, bs, &fc, c, cs);
#elif USE_MYZGEMM
   myzgemm_hn(fc, fa,
           a, as, am, an,
           b, bs, bm, bn,
           c, cs, cm, cn);
#else // Fortran BLAS
   ZGEMM_FTN(ZGEMM_CONJUGATE_TRANSPOSE, ZGEMM_NORMAL,
        &cm, &cn, &am,
        &fa, a, &as, b, &bs, &fc, c, &cs);
#endif
  }

// C = fc*C + fa*A*B'
void cmm_c8_gemm_nh(c8 fc, c8 fa,
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
#ifdef USE_CBLAS
   cblas_zgemm(CblasColMajor, CblasNoTrans, CblasConjTrans,
        cm, cn, an, &fa, a, as, b, bs, &fc, c, cs);
#elif USE_MYZGEMM
   myzgemm_nh(fc, fa,
           a, as, am, an,
           b, bs, bm, bn,
           c, cs, cm, cn);
#else // Fortran BLAS
   ZGEMM_FTN(ZGEMM_NORMAL, ZGEMM_CONJUGATE_TRANSPOSE,
        &cm, &cn, &an,
        &fa, a, &as, b, &bs, &fc, c, &cs);
#endif
  }

/*

   That should be it for the conditional code.

*/


// C = A*B
void cmm_c8_mm_nn(
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
   cmm_c8_gemm_nn(z_zero, z_one, a, as, am, an, b, bs, bm, bn, c, cs, cm, cn);
  }

// C = A'*B
void cmm_c8_mm_hn(
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
   cmm_c8_gemm_hn(z_zero, z_one, a, as, am, an, b, bs, bm, bn, c, cs, cm, cn);
  }

//HERE

// C = A*B'
void cmm_c8_mm_nh(
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
   cmm_c8_gemm_nh(z_zero, z_one, a, as, am, an, b, bs, bm, bn, c, cs, cm, cn);
  }

// C = C + A*B
void cmm_c8_mm_a_nn(
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
   cmm_c8_gemm_nn(z_one, z_one, a, as, am, an, b, bs, bm, bn, c, cs, cm, cn);
  }

// C = C + A'*B
void cmm_c8_mm_a_hn(
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
   cmm_c8_gemm_hn(z_one, z_one, a, as, am, an, b, bs, bm, bn, c, cs, cm, cn);
  }

// C = C + A*B'
void cmm_c8_mm_a_nh(
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
   cmm_c8_gemm_nh(z_one, z_one, a, as, am, an, b, bs, bm, bn, c, cs, cm, cn);
  }

// C = C - A*B
void cmm_c8_mm_s_nn(
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
   cmm_c8_gemm_nn(z_one, z_mone, a, as, am, an, b, bs, bm, bn, c, cs, cm, cn);
  }

// C = C - A'*B
void cmm_c8_mm_s_hn(
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
   cmm_c8_gemm_hn(z_one, z_mone, a, as, am, an, b, bs, bm, bn, c, cs, cm, cn);
  }

// C = C - A*B'
void cmm_c8_mm_s_nh(
        c8 *a, int as, int am, int an,
        c8 *b, int bs, int bm, int bn,
        c8 *c, int cs, int cm, int cn)
  {
   cmm_c8_gemm_nh(z_one, z_mone, a, as, am, an, b, bs, bm, bn, c, cs, cm, cn);
  }

// C = C * A
#define TMP_SIZE 4096

void cmm_c8_mt(c8 *cp, int cs, int cm, int cn, c8 *ap, int as, int am, int an)
  {
   c8 tmp[TMP_SIZE];
   int msize;
   c8 *ct;

   msize = cn*cm;
   if(msize > TMP_SIZE)
      ct = (c8 *) cmm_alloc(sizeof(c8)*msize);
   else
      ct = tmp;

   cmm_c8_m_copy(cp, cs, cm, cn, ct, cm);

   cmm_c8_mm_nn(ct, cm, cm, cn,
            ap, as, am, an,
            cp, cs, cm, cn);

   if(msize > TMP_SIZE)
     cmm_free(ct);

  }

void cmm_c8_m_copy(c8 *src, int ss, int sm, int sn, c8 *dst, int ds)
  {
   int j;

   for(j=0;j<sn;j++)
     cmm_c8_copy(sm, src + j*ss, 1, dst + j*ds, 1);
  }

c8 *cmm_c8_m_alloc(int am, int an)
  {
   return cmm_alloc(sizeof(c8)*am*an);
  }

void cmm_c8_m_clear(c8 *a, int as, int am, int an)
  {
   int i;
   for(i=0;i<an;i++)
     cmm_c8_zero(am, a+as*i);
  }

void cmm_c8_m_apply(c8 *a, int as, int am, int an, AppFunc func, void *p)
  {
   int i, j;
   for(j=0;j<an;j++)
     for(i=0;i<am;i++)
        func(i, j, a + i + j*as, p);
  }

void cmm_c8_m_print(c8 *a, int as, int am, int an, char *fmt)
  {
   int i, j;

   if(fmt == NULL)
     fmt = "%8.3lf";

   for(i=0;i<am;i++)
    {
     for(j=0;j<an;j++)
      {
       c8 *ap =  a + i + j*as;
       char rs[32],is[32];
       sprintf(rs, fmt, ap->r);
       sprintf(is, fmt, ap->i);
       printf("(%s %s) ",rs,is);
      }
     printf("\n");
    }
  }

// C = C + fa*A
void  cmm_c8_m_axpy(c8 fa, c8 *a, int astride, int am, int an, c8 *c, int cstride)
  {
   int i, j;

   for(j=0;j<an;j++)
     {
      c8 *ap = a + j*astride;
      c8 *cp = c + j*cstride;
      for(i=0;i<am;i++)
        {
         cp[i].r += fa.r * ap[i].r - fa.i * ap[i].i;
         cp[i].i += fa.r * ap[i].i + fa.i * ap[i].r;
        }
     }
  }

static void _func_identity(int i, int j, c8 *e, void *p)
  {
   if(i == j)
     {
      e->r = 1.;
      e->i = 0.;
     }
   else
     {
      e->r = 0.;
      e->i = 0.;
     }
  }

void cmm_c8_m_identity(c8 *a, int as, int n)
  {
   cmm_c8_m_apply(a, as, n, n, _func_identity, NULL);
  }

static void _func_norms(int i, int j, c8 *e, void *p)
  {
   r8 *norms = (r8 *) p;
   r8 av;

   av = ABS(e->r) + ABS(e->i);

   if(i == j)
     {
      norms[0] += e->r;
      norms[1] += e->i;
     }

   norms[2] += e->r*e->r + e->i*e->i;
   norms[3] = MAX(norms[3], av);
  }

void cmm_c8_m_norms(c8 *a, int as, int am, int an, r8 *norms)
  {
   r8_zero(4,norms);
   cmm_c8_m_apply(a, as, am, an, _func_norms, norms);
  }

static void _func_clear_upper(int i, int j, c8 *e, void *p)
  {
   if(i < j)
     {
      e->r = 0.;
      e->i = 0.;
     }
  }

void cmm_c8_m_clear_upper(c8 *a, int as, int am, int an)
  {
   cmm_c8_m_apply(a, as, am, an, _func_clear_upper, NULL);
  }

static void _func_clear_lower(int i, int j, c8 *e, void *p)
  {
   if(i > j)
     {
      e->r = 0.;
      e->i = 0.;
     }
  }

void cmm_c8_m_clear_lower(c8 *a, int as, int am, int an)
  {
   cmm_c8_m_apply(a, as, am, an, _func_clear_lower, NULL);
  }


