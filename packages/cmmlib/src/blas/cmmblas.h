/*! \file cmmblas.h
    \brief Interface for BLAS routines.
   \ingroup blas
*/
#ifndef CMMBLAS_H
#define CMMBLAS_H 1

#ifdef __cplusplus
extern "C"
{
#endif                          /* __cplusplus */

// Level 1 BLAS and miscellaneous routines

// add some routines to get machine precision info, byte size, byte ordering
// add some routines to copy blocks, gather/scatter blocks

void cmm_i_zero(int n, int *iv);
void cmm_i_set(int n, int c, int *iv);
void cmm_i_copy(int n, int *src, int src_stride, int *dst, int dst_stride);
void cmm_i_gather(int n, int *list, int *iu, r8 *iv);
void cmm_i_scatter(int n, int *list, int *iu, r8 *iv);
void cmm_i_range(int n, int *iv, int *ivmin, int *imin, int *ivmax, int *imax); 
void cmm_i_sort(int n, int *iv, int *list);

void cmm_r8_copy(int n, r8 *src, int src_stride, r8 *dst, int dst_stride);
void cmm_r8_gather(int n, int *list, r8 *u, r8 *v);
void cmm_r8_scatter(int n, int *list, r8 *u, r8 *v);
void cmm_r8_zero(int n, r8 *v); 
void cmm_r8_set(int n, r8 c, r8 *v); 
void cmm_r8_scale(int n, r8 c, r8 *v); 
r8   cmm_r8_sum(int n, r8 *v);
r8   cmm_r8_dot(int n, r8 *v, r8 *w);
r8   cmm_r8_norm(int n, r8 *v);
r8   cmm_r8_normsq(int n, r8 *v);
r8   cmm_r8_norm1(int n, r8 *v);
r8   cmm_r8_normi(int n, r8 *v);
void cmm_r8_range(int n, r8 *v, r8 *vmin, int *imin, r8 *vmax, int *imax);
void cmm_r8_axpy(int n, r8 a, r8 *x, r8 *y);
void cmm_r8_sort(int n, r8 *v, int *list);

void cmm_c8_copy(int n, c8 *src, int src_stride, c8 *dst, int dst_stride);
void cmm_c8_gather(int n, int *list, c8 *u, c8 *v);
void cmm_c8_scatter(int n, int *list, c8 *u, c8 *v);
void cmm_c8_zero(int n, c8 *v); 
void cmm_c8_set(int n, c8 c, c8 *v); 
void cmm_c8_scale(int n, c8 c, c8 *v); 
c8   cmm_c8_sum(int n, c8 *v);
c8   cmm_c8_dot(int n, c8 *v, c8 *w);
c8   cmm_c8_hdot(int n, c8 *v, c8 *w);
r8   cmm_c8_norm(int n, c8 *v);
r8   cmm_c8_normsq(int n, c8 *v);
void cmm_c8_axpy(int n, c8 a, c8 *x, c8 *y);

// Might want some routines that rethink a complex array as a real array
// and vice versa:
//
//   r8 *cmm_c8_cast_real(c8 *v)
//   c8 *cmm_r8_cast_complex(r8 *v)
//
//   void cmm_c8_conj(int n, c8 *v);
//   void cmm_c8_setreal(int n, r8 *src, c8 *dst);
//   void clbas_c8_getreal(int n, c8 *src, r8 *dst);
//   void cmm_c8_setimag(int n, r8 *src, c8 *dst);
//   void clbas_c8_getimag(int n, c8 *src, r8 *dst);
//

// r8 *cmm_r8_allocate(int n);
// c8 *cmm_c8_allocate(int n);

// void cmm_r8_free(r8 *);
// void cmm_c8_free(c8 *);

// Level 2 BLAS routines.

// Level 3 BLAS routines

// C = fc*C + fa*A*B
void cmm_c8_gemm_nn(c8 fc, c8 fa,
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = fc*C + fa*A'*B
void cmm_c8_gemm_hn(c8 fc, c8 fa,
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = fc*C + fa*A*B'
void cmm_c8_gemm_nh(c8 fc, c8 fa,
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = A*B
void cmm_c8_mm_nn(
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = A'*B
void cmm_c8_mm_hn(
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = A*B'
void cmm_c8_mm_nh(
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = C + A*B
void cmm_c8_mm_a_nn(
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = C + A'*B
void cmm_c8_mm_a_hn(
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = C + A*B'
void cmm_c8_mm_a_nh(
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = C - A*B
void cmm_c8_mm_s_nn(
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = C - A'*B
void cmm_c8_mm_s_hn(
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = C - A*B'
void cmm_c8_mm_s_nh(
        c8 *a, int astride, int am, int an,
        c8 *b, int bstride, int bm, int bn,
        c8 *c, int cstride, int cm, int cn);

// C = C * A
void cmm_c8_mt(c8 *cp, int cs, int cm, int cn, c8 *ap, int as, int am, int an);

// C = A * C
void cmm_c8_tm(c8 *cp, int cs, int cm, int cn, c8 *ap, int as, int am, int an);

// C = fa*A
void cmm_c8_m_scale(c8 fa, c8 *a, int astride, int am, int an, c8 *c, int cstride);

// C += fa*A
void  cmm_c8_m_axpy(c8 fa, c8 *a, int astride, int am, int an, c8 *c, int cstride);

// C += fa*A + fb*B
void cmm_c8_m_axpy2(c8 fa, c8 *a, int astride, int am, int an, c8 fb, c8 *b, int bstride, c8 *c, int cstride);

// C = fa*A + fc*C
void  cmm_c8_m_axpby(c8 fa, c8 *a, int astride, int am, int an, c8 fc, c8 *c, int cstride);

void cmm_c8_m_copy(c8 *src, int ss, int sm, int sn, c8 *dst, int ds);

// Utility functions
typedef void (*AppFunc)(int i, int j, c8 *e, void *p);

c8  *cmm_c8_m_alloc(int am, int an);
void cmm_c8_m_clear(c8 *a, int as, int am, int an);
void cmm_c8_m_apply(c8 *a, int as, int am, int an, AppFunc func, void *p);
void cmm_c8_m_print(c8 *a, int as, int am, int an, char *fmt);

// Useful macro
#define C8_M_ELEMENT(a, as, i, j) (((c8*)(a)) + (i) + (j)*(as))

// Create identity matrix
void cmm_c8_m_identity(c8 *a, int as, int an);

// Compute trace (2), sum of squares, max absolute value
void cmm_c8_m_norms(c8 *a, int as, int am, int an, r8 *norms);

// Clear upper half of matrix
void cmm_c8_m_clear_upper(c8 *a, int as, int am, int an);

 // Clear lower half of matrix
void cmm_c8_m_clear_lower(c8 *a, int as, int am, int an);


#ifdef CMM_SHORTHAND

#define int_zero(n, ix)          cmm_i_zero(n, ix)
#define int_set(n, c, ix)        cmm_i_set(n, c, ix)
#define r8_copy(n, x, ss, y, ds) cmm_r8_copy(n, x, ss, y, ds)
#define r8_zero(n, x)            cmm_r8_zero(n, x)
#define r8_axpy(n, a, x, y)      cmm_r8_axpy(n, a, x, y)
#define r8_dot(n, x, y)          cmm_r8_dot(n, x, y)
#define r8_norm(n, x)            cmm_r8_norm(n, x)
#define r8_scale(n, s, x)        cmm_r8_scale(n, s, x)
#define r8_sum(n, x)             cmm_r8_sum(n, x)

#define c8_copy(n, x, ss, y, ds) cmm_c8_copy(n, x, ss, y, ds)
#define c8_zero(n, x)            cmm_c8_zero(n, x)
#define c8_axpy(n, a, x, y)      cmm_c8_axpy(n, a, x, y)
#define c8_dot(n, x, y)          cmm_c8_dot(n, x, y)
#define c8_hdot(n, x, y)         cmm_c8_hdot(n, x, y)  // (X) dot (Y*)
#define c8_norm(n, v)            cmm_c8_norm(n, v)

#endif

#ifdef __cplusplus
}
#endif                          /* __cplusplus */

#endif
