/*! \file c8chol.c
    \brief Cholesky factorization.

   Cholesky factorization.

   Let A be a Hermitian matrix. This routine
   computes the factorization

      A = C * C' where C is lower triangular.

   Only the lower part of A is accessed. A is
   replaced with the factor C.


   Notes:

     For the blocked form I should try going
   backwards in the update loop. This means
   updating the next block at the end of the
   loop, so that it may still be in cache
   for the next step.

     It is probably useful to block the "unblocked"
   routine by hand.

*/

#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmatrix.h"

void cmm_c8_cholesky_unblocked(c8 *a, int as, int m, int n)
  {
   int i,j,k;

   for(j=0;j<m;j++)
     {
      c8 *lp = a + j + j*as;
      int mj = m-j;
      r8 adiag;

      if(j > 0)
        c8_zero(j, lp-j);

      adiag = 1./sqrt(lp->r);
      for(i=0;i<mj;i++)
        {
         lp[i].r *= adiag;
         lp[i].i *= adiag;
        }

      for(k=1; k<mj;k++)  //
        {
         c8 *lk = a + j + (j+k)*as;
         for(i=k;i<mj;i++)
           {
            lk[i].r -= lp[i].r * lp[k].r + lp[i].i * lp[k].i;
            lk[i].i -= lp[i].i * lp[k].r - lp[i].r * lp[k].i;
           }
        }
     }
  }

/*! \brief Cholesky factorization.
   \param a Hermitian matrix to factor (base)
   \param as Stride.
   \param m Number of rows.
   \param n Number of columns.

   The complex Hermitian matrix \f$ A \f$ is factored:

  \f$ A = L^{\dagger} L \f$.

   The factorization is in-place, the matrix \f$ A \f$ is overwritten
   with the lower triangular factor \f$ L \f$.

   This function is blocked and does use BLAS. For small matrices
   it will call the unblocked version.

*/
void cmm_c8_cholesky(c8 *a, int as, int ni, int nj)
  {
   int j;
   int b;

// Block size
   b = 16;

   for(j=0;j<nj;j+=b)
     {
      int nr = nj-j;
      int nr1 = nj-(j+b);
      int k;
      c8 *a00 = C8_M_ELEMENT(a, as, j,   j);
      c8 *ab0 = C8_M_ELEMENT(a, as, j+b, j);

// Handle remainder case using the unblocked algorithm
      if(nr <= b)
        {
         cmm_c8_cholesky_unblocked(a00, as, nr, nr);
         break;
        }

// Build the diagonal block
      cmm_c8_cholesky_unblocked(a00, as, b, b);

// Build the subdiagonal blocks
      cmm_c8_tsolve_alt_unblocked(a00, as, b,   b,
                              ab0, as, nr1, b,
                              ab0, as, nr1, b);

// Transform the remainder of the matrix
      for(k=j+b;k<nj;k+=b)
        {
         int nk = nj-k;
         int rk = MIN(b, nk);
         int nk1 = nj-(k+b);
         int k1;

       // Triangular part of the update
         for(k1=0;k1<rk;k1++) // this block
           {
            c8 *lk = C8_M_ELEMENT(a, as, k, k+k1);
            c8 *lj = C8_M_ELEMENT(a, as, k, j);
            int i,j1;

            for(j1=0;j1<b;j1++)
              {
               for(i=k1;i<rk;i++)
                 {
                  lk[i].r -= lj[i].r * lj[k1].r + lj[i].i * lj[k1].i;
                  lk[i].i -= lj[i].i * lj[k1].r - lj[i].r * lj[k1].i;
                 }
               lj += as;
              }
           }

         // Rectangular part of the update
         // A(k+bn:, k:k+b-1) = A(k+b:n, k:k+b-1) - A() * A()'

         if(nk1>0)
         cmm_c8_mm_s_nh(C8_M_ELEMENT(a, as, k+b, j), as, nk1, b,
                    C8_M_ELEMENT(a, as, k, j),   as, b,   b,
                    C8_M_ELEMENT(a, as, k+b, k), as, nk1, b);

        }
    }
  }

