
#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmband.h"

/*

    Compute trace and sum of squares of a Hermitian banded matrix

    The banded matrix has 2*b+1 diagonals. We only access the
    lower diagonals.

*/

void c8band_trace(int n, int b, r8 *A, r8 *tr, r8 *e2)
    {
     r8 e0,e1;
     int i,j;

     e0 = 0.;
     e1 = 0.;

     for(i=0;i<n;i++)
       {
        r8 z0, z1;
        r8 *trow = A + 2*(2*b+1)*i;

        z0 = trow[2*b];
        z1 = 0.;

        for(j=0;j<=b;j++)  // assume zero padding
          {
           r8 tr, ti;
           tr = trow[2*j]; ti = trow[2*j+1];
           if(j == b)
              z1 = z1 + (tr*tr+ti*ti);
           else
              z1 = z1 + 2.*(tr*tr+ti*ti);
          }

        e0 += z0;
        e1 += z1;

       }

     *tr = e0;
     *e2 = e1;
    }

/*
   Create a banded complex matrix
*/
r8 *c8band_new(int n, int bl, int bu, void *buf)
    {
     int rs;
     r8 *a;

     rs = 2*(bl+bu+1);

     if(buf)
       a = buf;
     else
       a = cmm_alloc(sizeof(r8)*n*rs);

     r8_zero(n*rs, a);
     return a;

    }

/*

   Create a Hermitian test matrix

*/
void c8band_test1(int n, int b, r8 *A)
    {
     int i, j;
     int rs;

     rs = (2*b+1)*2;

     for(i=0;i<n;i++)
       {
        r8 *tr = A + rs*i;

        for(j=0;j<2*b+1;j++)
          {
           int jj = i+j-b;
           if(jj < 0 || jj >= n)
             {
              tr[2*j] = 0.;
              tr[2*j+1] = 0.;
              continue;
             }

           if(jj == i)
             {
              tr[2*j] = DBLE(5.*i);
              tr[2*j+1] = 0.;
             }
           else
             {
              tr[2*j]   = DBLE(i+jj);
              tr[2*j+1] = DBLE(i-jj);
             }

          }
       }
    }

/*

  Print a banded matrix

*/
void c8band_print(int n, int bl, int bu, r8 *A)
    {
     int i,j,rs;
     rs = 2*(bl+bu+1);
     for(i=0;i<n;i++)
       {
        r8 *trow = A + rs * i;

        for(j=0;j<i-bl;j++)
          printf("(%lf %lf) ",0.,0.);

        for(j=i-bl;j<=i+bu;j++)
          {
          int jj = j-i+bl;

          if(j < 0 || j >= n)
            continue;
          printf("(%8.3lf %8.3lf) ",trow[2*jj],trow[2*jj+1]);
          }

        for(j=i+bu+1;j<n;j++)
          printf("(%lf %lf) ",0.,0.);

        printf("\n");
       }

    }
/*
void c8band_matvec(int n, int bl, int bu, r8 *A, r8 *x, r8 *y)
    {
     int m = bl + bu + 1;
     int i, j;
     c8 *xp;
     c8 *yp;
     c8 *ap;

     ap = A;
     xp = x;
     yp = y;

     for(i=0;i<n;i++)
       {
        int j0;
        int j1;

        j0 = i-bl;
        j1 = i+bu;

        if(j0 < 0)  j0 = 0;
        if(j1 >= n) j1 = n-1;

        yp[i] = c8_dot(j1-j0+1, ap+j0-i+bl, xp+j0);

        ap += m;
       }

    }
*/

void c8band_matvec(int n, int bl, int bu, r8 *a, r8 *x, r8 *y)
    {
     int i;
     int m = bl + bu + 1;
     c8 *xp;
     c8 *yp;
     c8 *ap;

     ap = (c8 *) a;
     xp = (c8 *) x;
     yp = (c8 *) y;

     r8_zero(2*n, y);
     for(i=0;i<n;i++)
       {
        int i0, i1;  // true extent of column i of a
        int l, off;  // length of extent and offset into packed storage

        i0  = MAX(0, i-bu);
        i1  = MIN(i+bl, n-1);
        l   = i1-i0+1;
        off = i0 - (i-bu);

        c8_axpy(l, xp[i], ap+off+m*i, yp+i0);
       }

    }

void c8band_matvec1(int n, int bl, int bu, r8 *a, r8 *x, r8 *y)
    {
     int i;
     int m = bl + bu + 1;
     c8 *xp;
     c8 *yp;
     c8 *ap;

     ap = (c8 *) a;
     xp = (c8 *) x;
     yp = (c8 *) y;

     r8_zero(2*n, y);
     for(i=0;i<n;i++)
       {
        int i0, i1;  // true extent of column i of a
        int l, off;  // length of extent and offset into packed storage

        i0  = MAX(0, i-bu);
        i1  = MIN(i+bl, n-1);
        l   = i1-i0+1;
        off = i0 - (i-bu);

        c8_axpy(l, xp[i], ap+off+m*i, yp+i0);

        {int k;
         c8 *v1 = ap + off + m*i;
         for(k=0;k<l;k++)
           {
           printf("       (%8.3lf %8.3lf) * (%8.3lf %8.3lf) = (%8.3lf %8.3lf)\n",
              xp[i].r, xp[i].i,v1->r,v1->i,xp[i].r*v1->r - xp[i].i*v1->i,  xp[i].i*v1->r + xp[i].r*v1->i);
           v1++;
           }
         }
       }

    }

r8 *c8band_element(int n, int bl, int bu, r8 *A, int i, int j)
    {
     int b1 = bu+bl+1;
    
     return A + 2*(bl+i-j + b1*j);
    }

c8 *cband_element(int n, int bl, int bu, c8 *A, int i, int j)
    {
     int b1 = bu+bl+1;
    
     return A + (bl+i-j + b1*j);
    }

