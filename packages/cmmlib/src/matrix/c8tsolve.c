/*! \file c8tsolve.c
    \brief Triangular solve routines.

    Let A be a lower triangular matrix. Then

    tsolve_la (L, A, Y)   Y = inv(L)*A
    tsolve_lta(L, A, Y)   Y = inv(L')*A
    tsolve_al (L, A, Y)   Y = A*inv(L)
    tsolve_alt(L, A, Y)   Y = A*inv(L')

    Comments: 

     1) We allow A to be overwritten with
        the solution Y. This is very useful
        in practice.

     2) L is square, but A, Y need not be
        square.

     3) Make notes on parallel implementation.

     4) The main routines are blocked and recursive.
        They call the unblocked routines to implement
        the triangular solve for the diagonal blocks
        and then use standard matrix multiply routines
        to update the remainder of the solution.

    To do:

       tsolve_al

       upper triangular versions

       versions which assume that A is triangular.

       versions which assume that A is Hermitian and
       only access the lower part.

       need to clean up the unblocked code for
       efficiency.

*/

#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmatrix.h"

// unblocked code for tsolve_la (L, A, Y),   Y = inv(L)*A
void cmm_c8_tsolve_la_unblocked(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn)
  {
   int i, j, k;

   if(y != a)
     cmm_c8_m_copy(a, as, am, an, y, ys);   // First step is to copy A to Y if not the same.

   for(i=0;i<ym;i++)
     {
      c8 *lp = l + i*ls;
      c8 t; r8 tn;

      t = lp[i];
      tn = 1./(t.r*t.r + t.i*t.i);

      for(j=0;j<yn;j++)
        {
         c8 *yp = y + j*ys;
         c8 y0, y1;

         y0 = yp[i];
         y1.r = (y0.r * t.r + y0.i * t.i)*tn;
         y1.i = (y0.i * t.r - y0.r * t.i)*tn;
         yp[i] = y1;

         for(k=i+1;k<ym;k++)
           {
            yp[k].r -= lp[k].r * y1.r - lp[k].i * y1.i;
            yp[k].i -= lp[k].r * y1.i + lp[k].i * y1.r;
           }
        } // end j loop
     } // end i loop

  }

/*! \brief Triangular solve tsolve_la (L, A, Y)   Y = inv(L)*A
   \param l  Triangular factor \f$ L \f$ (base).
   \param ls Stride for \f$ L \f$.
   \param lm Number of rows for \f$ L \f$.
   \param ln Number of columns for \f$ L \f$.
   \param a  Full matrix \f$ A \f$ (base).
   \param as Stride for \f$ A \f$.
   \param am Number of rows for \f$ A \f$.
   \param an Number of columns for \f$ A \f$.
   \param y  Output matrix \f$ Y \f$ (base).
   \param ys Stride for \f$ Y \f$.
   \param ym Number of rows for \f$ Y \f$.
   \param yn Number of columns for \f$ Y \f$.

   Given the lower triangular matrix \f$ L \f$ and the full matrix
   \f$ A \f$, we compute  \f$ Y = L^{-1} A \f$.

   No storage is allocated. If memory for \f$ A \f$ and
   \f$ Y \f$ are the same then the computation is done in-place
   and a copy is avoided.

   Note that \f$ L \f$ must be square and dimensions of
  \f$ A, Y \f$ must match.

*/
void cmm_c8_tsolve_la(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn)
  {
   int i;
   int b;

   b = 16;

   if(y != a)
     cmm_c8_m_copy(a, as, am, an, y, ys);   // First step is to copy A to Y if not the same.

   for(i=0;i<lm;i += b)
     {
      int rb = MIN(b, lm-i);
      int ni = lm - (i+rb);
      c8 *lii = l + i + i * ls;
      c8 *yi = y + i;

      cmm_c8_tsolve_la_unblocked(lii, ls, rb, rb,
                              yi, ys, rb, yn,
                              yi, ys, rb, yn);

      cmm_c8_mm_s_nn(lii + rb, ls, ni, rb,
                 yi,       ys, rb, yn,
                 y + i+rb, ys, ni, yn);

     }

  }

void cmm_c8_tsolve_la_old(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn)
  {
   int i, j;
   int b;

   b = 16;

   if(y != a)
     cmm_c8_m_copy(a, as, am, an, y, ys);   // First step is to copy A to Y if not the same.

   for(j=0;j<yn;j+=b)  // blocked j-loop
     {
      int rb = MIN(b, yn-j);
      for(i=0;i<ym;i+=b) // blocked i-loop
        {
         int nr = ym-i;
         int r = MIN(b, nr);
         c8 *lp0 = l  + i + i*ls;
         c8 *yp0 = y  + i + j*ys;
         c8 *lp = lp0;
         int i1,j1,k1;

         for(i1=0;i1<r;i1++)
           {
            c8 l0 = lp[i1];
            r8 tn = 1./(l0.r*l0.r + l0.i*l0.i);
            c8 *yp = yp0;

            for(j1=0;j1<rb;j1++)
              {
               c8 yt1 = yp[i1];
               c8 yt2;

               yt2.r = tn*(yt1.r*l0.r + yt1.i*l0.i);
               yt2.i = tn*(yt1.i*l0.r - yt1.r*l0.i);
               yp[i1] = yt2;

               for(k1=i1+1;k1<r;k1++)
                 {
                  yp[k1].r -= lp[k1].r * yt2.r - lp[k1].i * yt2.i;
                  yp[k1].i -= lp[k1].i * yt2.r + lp[k1].r * yt2.i;
                 }
               yp += ys;
              }
            lp += ls;
           }

         cmm_c8_mm_s_nn(
                    C8_M_ELEMENT(l, ls, i+r, i), ls, ym-(i+r), r, 
                    C8_M_ELEMENT(y, ys, i,   j), ys, r,        rb,
                    C8_M_ELEMENT(y, ys, i+r, j), ys, ym-(i+r), rb);

        }  // end loop over i-block
     } // end loop over j-block
  }

void cmm_c8_tsolve_alt_unblocked(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn)
  {
   int i, j, k;
   int n = lm;

   if(y != a)
     cmm_c8_m_copy(a, as, am, an, y, ys);

   for(j=0;j<n;j++)
     {
      c8 *lp = C8_M_ELEMENT(l, ls, 0, j);
      c8 *yp = C8_M_ELEMENT(y, ys, 0, j);
      c8 ljj = lp[j];
      c8 *yk;

     {
       r8 t =  1./(ljj.r*ljj.r + ljj.i*ljj.i);
       ljj.r = t*ljj.r;
       ljj.i = t*ljj.i;
     }

      for(i=0;i<ym;i++)
        {
         c8 cij = yp[i];
         yp[i].r = cij.r * ljj.r - cij.i * ljj.i;
         yp[i].i = cij.i * ljj.r + cij.r * ljj.i;
        }

      yk = yp;
      for(k=j+1;k<n;k++)
        {
         c8 lkj = lp[k];
         yk += ys;
         for(i=0;i<ym;i++)
           {
            yk[i].r -= yp[i].r * lkj.r + yp[i].i * lkj.i;
            yk[i].i -= yp[i].i * lkj.r - yp[i].r * lkj.i;
           }
        }
     }
  }

/*! \brief Triangular solve tsolve_alt (L, A, Y)   Y = A*inv(L')
   \param l  Triangular factor \f$ L \f$ (base).
   \param ls Stride for \f$ L \f$.
   \param lm Number of rows for \f$ L \f$.
   \param ln Number of columns for \f$ L \f$.
   \param a  Full matrix \f$ A \f$ (base).
   \param as Stride for \f$ A \f$.
   \param am Number of rows for \f$ A \f$.
   \param an Number of columns for \f$ A \f$.
   \param y  Output matrix \f$ Y \f$ (base).
   \param ys Stride for \f$ Y \f$.
   \param ym Number of rows for \f$ Y \f$.
   \param yn Number of columns for \f$ Y \f$.

   Given the lower triangular matrix \f$ L \f$ and the full matrix
   \f$ A \f$, we compute  \f$ Y = A L^{-\dagger} \f$.

   No storage is allocated. If memory for \f$ A \f$ and
   \f$ Y \f$ are the same then the computation is done in-place
   and a copy is avoided.

   Note that \f$ L \f$ must be square and dimensions of
  \f$ A, Y \f$ must match.

*/
void cmm_c8_tsolve_alt(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn)
  {
   int i;
   int b;

   b = 16;

   if(y != a)
     cmm_c8_m_copy(a, as, am, an, y, ys);

   for(i=0;i<lm;i += b)
     {
      int rb;
      int ni;

      rb = MIN(b, lm-i);
      cmm_c8_tsolve_alt_unblocked(l + i + i *ls, ls, rb, rb,
                              y + i*ys,      ys, ym, rb,
                              y + i*ys,      ys, ym, rb);

      ni = yn-(i+rb);
      if(ni > 0)
      cmm_c8_mm_s_nh(y + (i)*ys,      ys,        ym,        rb, 
                 l + (i+rb)+i*ls, ls,        ni,        rb, 
                 y + (i+rb)*ys,   ys,        ym,        ni);

     }
  }

void cmm_c8_tsolve_alt_old(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn)
  {
   int i, j;
   int b, bi;

/*
   If ym gets large (>128 for laptop) then need
   to block over i, this is not implemented yet.
*/

   b = 16;
   bi = 256;

   if(y != a)
     cmm_c8_m_copy(a, as, am, an, y, ys);

   for(j=0;j<lm;j += b)
     {
      int nr = lm-j;
      int r = MIN(b, nr);
      int j1,k1;
      c8 *lp = l + j*ls;
      c8 *yp = y + j*ys;

      for(j1=0;j1<r;j1++)
        {
         c8 t = lp[j+j1];
         r8 tn = 1./(t.r*t.r + t.i*t.i);
         c8 *yk;

         t.r *= tn;
         t.i *= tn;

         for(i=0;i<ym;i++)
           {
            c8 yt;
            yt.r = yp[i].r*t.r - yp[i].i*t.i;
            yt.i = yp[i].r*t.i + yp[i].i*t.r;
            yp[i] = yt;
           }

         yk = yp;

         for(k1=j1+1;k1<r;k1++)
           {
            c8 t1 = lp[j+k1];
            yk += ys;
            for(i=0;i<ym;i++)
              {
               yk[i].r -= yp[i].r * t1.r + yp[i].i * t1.i;
               yk[i].i -= yp[i].i * t1.r - yp[i].r * t1.i;
              }
           }

        yp += ys;
        lp += ls;
        }

      cmm_c8_mm_s_nh(
            C8_M_ELEMENT(y, ys,   0,   j), ys,   ym,    r,
            C8_M_ELEMENT(l, ls, j+r,   j), ls, nr-r,    r,
            C8_M_ELEMENT(y, ys,   0, j+r), ys,   ym, nr-r);

     }
  }

void cmm_c8_tsolve_lta_unblocked(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn)
  {
   int i, j, k;

   if(y != a)
     cmm_c8_m_copy(a, as, am, an, y, ys);

   for(i=lm-1;i>=0;i--)
     {
      c8 *lp = l + i*ls;
      c8 t = lp[i];
      r8 tn = 1./(t.r*t.r + t.i*t.i);

      t.r *= tn;
      t.i *= tn;
         
      for(j=0;j<yn;j++)
        {
         c8 *yp = y + j*ys;

         c8 y0;

         y0.r = yp[i].r * t.r - yp[i].i * t.i;
         y0.i = yp[i].r * t.i + yp[i].i * t.r;
         yp[i] = y0;

         for(k=0;k<i;k++)
           {
            c8 t1 = l[i+k*ls];
            yp[k].r -= t1.r * y0.r + t1.i * y0.i;
            yp[k].i -= t1.r * y0.i - t1.i * y0.r;
           }
        }
     }

/*
   for(j=0;j<yn;j++)
     {
      c8 *yp = y + j*ys;

      for(i=lm-1;i>=0;i--)
        {
         c8 *lp = l + i*ls;
         c8 t = lp[i];
         r8 tn = 1./(t.r*t.r + t.i*t.i);
         c8 y0;

         t.r *= tn;
         t.i *= tn;

         y0.r = yp[i].r * t.r - yp[i].i * t.i;
         y0.i = yp[i].r * t.i + yp[i].i * t.r;
         yp[i] = y0;

         for(k=0;k<i;k++)
           {
            c8 t1 = l[i+k*ls];
            yp[k].r -= t1.r * y0.r + t1.i * y0.i;
            yp[k].i -= t1.r * y0.i - t1.i * y0.r;
           }
        }
     }
*/

  }

void cmm_c8_tsolve_lta(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn)
  {
   int i, r;
   int b;

   b = 16;

   if(y != a)
     cmm_c8_m_copy(a, as, am, an, y, ys);

   r = lm % b;
   if(r == 0)
     r = b;

   for( i=lm-r ; i >= 0; i -= b)
     {
      cmm_c8_tsolve_lta_unblocked(l + i + i*ls, ls, r,  r,
                              y + i,        ys, r, yn,
                              y + i,        ys, r, yn);
 
      if(i>0)
      cmm_c8_mm_s_hn(l+i,           ls,  r, i,
                 y+i,           ys,  r, yn,
                 y,             ys,  i, yn);

      r = b;
     }
   
  }

void cmm_c8_tsolve_lta_old(
             c8 *l, int ls, int lm, int ln,
             c8 *a, int as, int am, int an,
             c8 *y, int ys, int ym, int yn)
  {
   int j, k;
   int b;

   b = 16;

   if(y != a)
     cmm_c8_m_copy(a, as, am, an, y, ys);

   for(j=0;j<yn;j+=b)
     {
      int rb = MIN(b, yn-j);
      c8 *yp0 = y + j*ys;
      int j1,i0;

      for(i0=ym;i0>0;i0 -= b)
        {
         int i1 = i0-1;
         int ii;
         int  r = MIN(b, i0);
         c8 *lp = l + i1*ls;
         int i2;

         for(i2=0;i2<r;i2++)
           {
            c8 *yp = yp0;
            c8 t = lp[i1-i2];
            r8 tn = 1./(t.r*t.r + t.i*t.i);

            t.r *= tn;
            t.i *= tn;
            for(j1=0;j1<rb;j1++)
              {
               c8 y0, y1;
               y0 = yp[i1-i2];
               y1.r = y0.r * t.r - y0.i * t.i;
               y1.i = y0.r * t.i + y0.i * t.r;
               yp[i1-i2] = y1;

               for(k=i0-r;k<i1-i2;k++)
                 {
                  c8 t1 = l[i1-i2+k*ls];
                  yp[k].r -= t1.r * y1.r + t1.i * y1.i;
                  yp[k].i -= t1.r * y1.i - t1.i * y1.r;
                 }
               yp += ys;
              }
            lp -= ls;
           }

         ii = i0-r;
         if(ii > 0)
         cmm_c8_mm_s_hn(
                     C8_M_ELEMENT(l, ls,    ii, 0), ls,  r,  ii,
                     C8_M_ELEMENT(y, ys,    ii, j), ys,  r,  rb,
                     C8_M_ELEMENT(y, ys,     0, j), ys, ii,  rb);
        } // end of i-block
     } // end of j-block
  }
