/*! \file cmmatrix.h
    \brief Interface for matrix routines.
   \ingroup matrix
*/
#ifndef CMMMATRIX_H
#define CMMMATRIX_H 1

#ifdef __cplusplus
extern "C"
{
#endif                          /* __cplusplus */

// Compute cholesky factorization
void cmm_c8_cholesky_unblocked(c8 *a, int as, int an, int m);
void cmm_c8_cholesky(c8 *a, int as, int an, int m);

// Compute inv(L)*A
void cmm_c8_tsolve_la_unblocked(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn);
void cmm_c8_tsolve_la(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn);

// Compute A*inv(L')
void cmm_c8_tsolve_alt_unblocked(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn);
void cmm_c8_tsolve_alt(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn);

void cmm_c8_tsolve_lta_unblocked(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn);

void cmm_c8_tsolve_lta(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn);

// Build block of Householder vectors
void c8_build_house(c8 *a, int as, int am, int an, c8 *b, int bs);

// Build WY transform from Householder vectors
void c8_build_wy(c8 *y, int ys, int ym, int yn, c8 *w, int ws);

// Update A' = (I+WY')*A*(I+WY')', accessing lower half only
void c8_update_wy(
             c8 *y, int ys, int ym, int yn,
             c8 *w, int ws, int wn, int wm,
             c8 *a, int as, int am, int an,
             c8 *tmp, int ntmp
#ifdef DO_CLOCKS
 ,int *clocks
#endif
);

// Hermitian matrix-vector product, accessing lower half only
void c8_mvl(c8 *a, int as, int am, int an,
            c8 *x, int xs, int xm, int xn,
            c8 *y, int ys, int ym, int yn);

// Hermitian matrix update, accessing lower half only
void   c8_updl(
           c8 *a, int as,  int am, int an,
           c8 *w, int ws,  int wn, int wm,
           c8 *v, int vs,  int vn, int vm);

// Condense matrix to banded form.
void c8_condense(c8 *a, int as, int an, int b, c8 *aband);

// Back-transform eigenvectors 
void c8_house_back(c8 *a, int as, int am, int an, int bw, c8 *v, int vs, int vn, int vm);

// Compute eigenvalues and eigenvectors of ordinary Hermitian system
void c8_eigen(c8 *a, int as, int an, c8 *v, int vs, r8 *eigenvals);

// Compute eigenvalues and eigenvectors of generalized Hermitian system
void c8_geigen(c8 *a, int as, int an, c8 *b, int bs, c8 *v, int vs, r8 *eigenvals);

#ifdef __cplusplus
}
#endif                          /* __cplusplus */

#endif
