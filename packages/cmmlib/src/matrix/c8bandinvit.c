/*

   Routines for inverse iteration on banded complex matrices.

*/

#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmband.h"

void c8_householder(int n, r8 *a, r8 *v);
void c8_apply_householder(int n, r8 *v, r8 *x);

#ifdef DO_CLOCKS
int pcc_delta();
#endif

/*

    Compute the iteration (RR') x = y 

*/
static void c8band_invit_general(int n, int bl, int bu, r8 *r, r8 *y, r8 *x);
static void c8band_invit_2(int n, r8 *r, r8 *y, r8 *x);

void c8band_invit(int n, int bl, int bu, r8 *r, r8 *y, r8 *x)
    {
     if((n < 8) || (bl != 2) || (bu != 2))
       c8band_invit_general(n, bl, bu, r, y, x);
     else
       c8band_invit_2(n, r, y, x); 
    }

static void c8band_invit_general(int n, int bl, int bu, r8 *r, r8 *y, r8 *x)
    {
     int i,mq,mr,bmu;
     r8 *y1, *y2, *v;
     int iclock;

     mq = bl+1;
     bmu = bl + bu;
     mr = bmu+1;

#ifdef DO_CLOCKS
//     iclock = pcc_delta();
#endif
     y1 = cmm_alloc(sizeof(r8)*n*2);
     y2 = cmm_alloc(sizeof(r8)*n*2);
#ifdef DO_CLOCKS
//     iclock = pcc_delta();
//     printf("malloc clocks %d\n",iclock);
#endif

     // Solve R*y1 = y. Use y2 as temp storage

     v = y2;
     c8_zero(n, (c8 *) v);

#ifdef DO_CLOCKS
//     iclock = pcc_delta();
#endif
     for(i=n-1;i>0;i--)
       {
        int i0, i1;
        int l, off;
        r8 *rp;
        r8 zr, zi, tr, ti, rn;
        c8 zx;

        // solve for x[i]

        rp = r + 2*(bmu + i*mr);

        zr = rp[0]; // diagonal element
        zi = rp[1];

        rn = zr*zr + zi*zi;

        tr = y[2*i]   - v[2*i];
        ti = y[2*i+1] - v[2*i+1];

        zx.r = (tr*zr + ti*zi)/rn;
        zx.i = (ti*zr - tr*zi)/rn;
        y1[2*i]   = zx.r;
        y1[2*i+1] = zx.i;

        // update v(0:i-1) = v(0:i-1) + x(i) * r(0:i-1, i);
        i0  = MAX(0, i-bmu);
        i1  = MIN(i-1, n-1);  // only update to column i-1.
        l   = i1-i0+1;
        off = i0 - (i-bmu);

//        printf("forward i %d l %d i0 %d i1 %d off %d\n",i,l,i0,i1,off);

        if(l == 4) // common case
          {
           register r8 *xp = (r8 *) (r+2*(off+mr*i));
           register r8 *yp = (r8 *) (v+2*i0);
           register r8 tr = zx.r;
           register r8 ti = zx.i;

           yp[0] += tr * xp[0] - ti * xp[1];
           yp[1] += ti * xp[0] + tr * xp[1];
           yp[2] += tr * xp[2] - ti * xp[3];
           yp[3] += ti * xp[2] + tr * xp[3];
           yp[4] += tr * xp[4] - ti * xp[5];
           yp[5] += ti * xp[4] + tr * xp[5];
           yp[6] += tr * xp[6] - ti * xp[7];
           yp[7] += ti * xp[6] + tr * xp[7];

          }
        else
           c8_axpy(l, zx, (c8 *) (r+2*(off+mr*i)), (c8 *) (v+2*i0));

       }
#ifdef DO_CLOCKS
//     iclock = pcc_delta();
//     printf("forward clocks %d\n",iclock);
#endif

       {
        r8 *rp;
        r8 zr, zi, tr, ti, rn;
        c8 zx;

        // solve for x[0]

        rp = r + 2*bmu;

        zr = rp[0];
        zi = rp[1];

        rn = zr*zr + zi*zi;

        tr = y[0] - v[0];
        ti = y[1] - v[1];

        zx.r = (tr*zr + ti*zi)/rn;
        zx.i = (ti*zr - tr*zi)/rn;
        y1[0] = zx.r;
        y1[1] = zx.i;
       }

/*
     // check the solve.

      c8band_matvec(n, 0, bl+bu, r, y1, y2);

      for(i=0;i<n;i++)
        {
         r8 er, ei;
         er = y2[2*i]   - y[2*i];
         ei = y2[2*i+1] - y[2*i+1];
         printf("i %d y - R*x  %lf %lf\n",i,er,ei);
        }

     m_pause("check y-R*x");
*/

     // Solve R'*x = y1.  We don't actually modify R to do this. No temp storage needed.

     for(i=0;i<n;i++)
       {
        int i0, i1;
        int l, off;
        r8 *rp;
        r8 zr, zi, tr, ti, rn;
        c8 zx;

        // compute x(1:i-1)*conj(r(1:i-1,i)

        i0  = MAX(0, i-bmu);
        i1  = MIN(i-1, n-1);  // only update to row i-1.
        l   = i1-i0+1;
        off = i0 - (i-bmu);

//        printf("bakward i %d l %d i0 %d i1 %d off %d\n",i,l,i0,i1,off);

        rp = r + 2*(off+mr*i);

        if(l == 4) // common case
          {
           register r8 tr,ti;
           register r8 *xp = (r8 *) (x + 2*i0);
           register r8 *yp = (r8 *) (r + 2*(off+mr*i));

           tr  = xp[0] * yp[0] + xp[1] * yp[1];
           ti  = xp[1] * yp[0] - xp[0] * yp[1];
           tr += xp[2] * yp[2] + xp[3] * yp[3];
           ti += xp[3] * yp[2] - xp[2] * yp[3];
           tr += xp[4] * yp[4] + xp[5] * yp[5];
           ti += xp[5] * yp[4] - xp[4] * yp[5];
           tr += xp[6] * yp[6] + xp[7] * yp[7];
           ti += xp[7] * yp[6] - xp[6] * yp[7];

           zx.r = tr;
           zx.i = ti;
          }

        else if(l > 0)
          {
           zx = c8_hdot(l, (c8 *) (x + 2*i0), (c8 *) (r + 2*(off+mr*i)));
          }

        else
          {
           zx.r = 0.;
           zx.i = 0.;
          }

        // solve for y[i]

        rp = r + 2*(bmu + i*mr);

        zr = rp[0];
        zi = rp[1];
        rn = zr*zr + zi*zi;

        tr = y1[2*i]   - zx.r;
        ti = y1[2*i+1] - zx.i;

        x[2*i]   = (tr*zr - ti*zi)/rn;
        x[2*i+1] = (tr*zi + ti*zr)/rn;

       }

#ifdef DO_CLOCKS
//     iclock = pcc_delta();
#endif
     cmm_free(y1);
     cmm_free(y2);
#ifdef DO_CLOCKS
//     iclock = pcc_delta();
//     printf("free clocks %d\n",iclock);
#endif
    }

// Special case for bl = bu = 2
static void c8band_invit_2(int n, r8 *r, r8 *y, r8 *x)
    {
     int bl = 2;
     int bu = 2;
     int i,mq,mr,bmu;
     r8 *y1, *y2, *v;
     int iclock;

     mq  = 3;
     bmu = 4;
     mr  = 5;

     y1 = cmm_alloc(sizeof(r8)*n*2);
     y2 = cmm_alloc(sizeof(r8)*n*2);

     // Solve R*y1 = y. Use y2 as temp storage

     v = y2;
     c8_zero(n, (c8 *) v);

// Break forward loop into two peices
    {
     r8 *rp = r + 2*(bmu + (n-1)*mr);
     r8 *vp = v + 2*(n-1);
     r8 *yp  = y + 2*(n-1);
     r8 *y1p = y1 + 2*(n-1);
     register int ii;

     for(ii=n-4;ii>0;ii--) // for(i=n-1;i>3;i--)
       {
        r8 zr, zi, tr, ti, rn;
        r8 yr, yi;

        // solve for x[i]

        zr = rp[0]; // diagonal element
        zi = rp[1];

        rn = 1./(zr*zr + zi*zi);

        tr = yp[0] - vp[0];
        ti = yp[1] - vp[1];

        yr = rn*(tr*zr + ti*zi);
        yi = rn*(ti*zr - tr*zi);
        y1p[0] = yr;
        y1p[1] = yi;

// 14 + 32 flops per iteration = 46 flops/it
        vp[-8+0] += yr * rp[-8+0] - yi * rp[-8+1];
        vp[-8+1] += yi * rp[-8+0] + yr * rp[-8+1];
        vp[-8+2] += yr * rp[-8+2] - yi * rp[-8+3];
        vp[-8+3] += yi * rp[-8+2] + yr * rp[-8+3];
        vp[-8+4] += yr * rp[-8+4] - yi * rp[-8+5];
        vp[-8+5] += yi * rp[-8+4] + yr * rp[-8+5];
        vp[-8+6] += yr * rp[-8+6] - yi * rp[-8+7];
        vp[-8+7] += yi * rp[-8+6] + yr * rp[-8+7];

        rp -= 10;
        vp -= 2;
        yp -= 2;
        y1p -= 2;
       }

     }

     for(i=3;i>0;i--)
       {
        int i0, i1;
        int l, off;
        r8 *rp;
        r8 zr, zi, tr, ti, rn;
        c8 zx;

        // solve for x[i]

        rp =  r + 2*(bmu + i*mr);
        zr = rp[0]; // diagonal element
        zi = rp[1];

        rn = zr*zr + zi*zi;

        tr = y[2*i]   - v[2*i];
        ti = y[2*i+1] - v[2*i+1];

        zx.r = (tr*zr + ti*zi)/rn;
        zx.i = (ti*zr - tr*zi)/rn;
        y1[2*i]   = zx.r;
        y1[2*i+1] = zx.i;

        // update v(0:i-1) = v(0:i-1) + x(i) * r(0:i-1, i);
        i0  = MAX(0, i-bmu);
        i1  = MIN(i-1, n-1);  // only update to column i-1.
        l   = i1-i0+1;
        off = i0 - (i-bmu);

//        printf("forward i %d l %d i0 %d i1 %d off %d\n",i,l,i0,i1,off);

        if(l == 4) // common case
          {
           register r8 *xp = (r8 *) (r+2*(off+mr*i));
           register r8 *yp = (r8 *) (v+2*i0);
           register r8 tr = zx.r;
           register r8 ti = zx.i;

           yp[0] += tr * xp[0] - ti * xp[1];
           yp[1] += ti * xp[0] + tr * xp[1];
           yp[2] += tr * xp[2] - ti * xp[3];
           yp[3] += ti * xp[2] + tr * xp[3];
           yp[4] += tr * xp[4] - ti * xp[5];
           yp[5] += ti * xp[4] + tr * xp[5];
           yp[6] += tr * xp[6] - ti * xp[7];
           yp[7] += ti * xp[6] + tr * xp[7];

          }
        else
           c8_axpy(l, zx, (c8 *) (r+2*(off+mr*i)), (c8 *) (v+2*i0));

       }

       {
        r8 *rp;
        r8 zr, zi, tr, ti, rn;
        c8 zx;

        // solve for x[0]

        rp = r + 2*bmu;

        zr = rp[0];
        zi = rp[1];

        rn = zr*zr + zi*zi;

        tr = y[0] - v[0];
        ti = y[1] - v[1];

        zx.r = (tr*zr + ti*zi)/rn;
        zx.i = (ti*zr - tr*zi)/rn;
        y1[0] = zx.r;
        y1[1] = zx.i;
       }

/*
     // check the solve.

      c8band_matvec(n, 0, bl+bu, r, y1, y2);

      for(i=0;i<n;i++)
        {
         r8 er, ei;
         er = y2[2*i]   - y[2*i];
         ei = y2[2*i+1] - y[2*i+1];
         printf("i %d y - R*x  %lf %lf\n",i,er,ei);
        }

     m_pause("check y-R*x");
*/

     // Solve R'*x = y1.  We don't actually modify R to do this. No temp storage needed.

// Break reverse loop into two peices

     for(i=0;i<4;i++)
       {
        int i0, i1;
        int l, off;
        r8 *rp;
        r8 zr, zi, tr, ti, rn;
        c8 zx;

        // compute x(1:i-1)*conj(r(1:i-1,i)

        i0  = MAX(0, i-bmu);
        i1  = MIN(i-1, n-1);  // only update to row i-1.
        l   = i1-i0+1;
        off = i0 - (i-bmu);

//        printf("bakward i %d l %d i0 %d i1 %d off %d\n",i,l,i0,i1,off);

        rp = r + 2*(off+mr*i);

        if(l == 4) // common case
          {
           register r8 tr,ti;
           register r8 *xp = (r8 *) (x + 2*i0);
           register r8 *yp = (r8 *) (r + 2*(off+mr*i));

           tr  = xp[0] * yp[0] + xp[1] * yp[1];
           ti  = xp[1] * yp[0] - xp[0] * yp[1];
           tr += xp[2] * yp[2] + xp[3] * yp[3];
           ti += xp[3] * yp[2] - xp[2] * yp[3];
           tr += xp[4] * yp[4] + xp[5] * yp[5];
           ti += xp[5] * yp[4] - xp[4] * yp[5];
           tr += xp[6] * yp[6] + xp[7] * yp[7];
           ti += xp[7] * yp[6] - xp[6] * yp[7];

           zx.r = tr;
           zx.i = ti;
          }

        else if(l > 0)
          {
           zx = c8_hdot(l, (c8 *) (x + 2*i0), (c8 *) (r + 2*(off+mr*i)));
          }

        else
          {
           zx.r = 0.;
           zx.i = 0.;
          }

        // solve for y[i]

        rp = r + 2*(bmu + i*mr);

        zr = rp[0];
        zi = rp[1];
        rn = zr*zr + zi*zi;

        tr = y1[2*i]   - zx.r;
        ti = y1[2*i+1] - zx.i;

        x[2*i]   = (tr*zr - ti*zi)/rn;
        x[2*i+1] = (tr*zi + ti*zr)/rn;

       }

   {
    r8 *xp = x;
    r8 *rp = r;
    r8 *y1p = y1 + 8;
    register int ii;

//     for(i=4;i<n;i++)
    for(ii=n-4;ii>0;ii--)
       {
        r8 zr, zi, tr, ti, rn;
        r8 zxr, zxi;

        zxr  = xp[0] * rp[40 + 0] + xp[1] * rp[40 + 1];
        zxi  = xp[1] * rp[40 + 0] - xp[0] * rp[40 + 1];
        zxr += xp[2] * rp[40 + 2] + xp[3] * rp[40 + 3];
        zxi += xp[3] * rp[40 + 2] - xp[2] * rp[40 + 3];
        zxr += xp[4] * rp[40 + 4] + xp[5] * rp[40 + 5];
        zxi += xp[5] * rp[40 + 4] - xp[4] * rp[40 + 5];
        zxr += xp[6] * rp[40 + 6] + xp[7] * rp[40 + 7];
        zxi += xp[7] * rp[40 + 6] - xp[6] * rp[40 + 7];

        zr = rp[48+0];
        zi = rp[48+1];
        rn = zr*zr + zi*zi;

        tr = y1p[0] - zxr;
        ti = y1p[1] - zxi;

        xp[8+0]   = (tr*zr - ti*zi)/rn;
        xp[8+1]   = (tr*zi + ti*zr)/rn;

        xp += 2;
        y1p += 2;
        rp += 10;
       }
      
     }
#ifdef DO_CLOCKS
//     iclock = pcc_delta();
#endif
     cmm_free(y1);
     cmm_free(y2);
#ifdef DO_CLOCKS
//     iclock = pcc_delta();
//     printf("free clocks %d\n",iclock);
#endif
    }

/*

   Apply a unitary matrix in packed (banded) Householder form to a vector

      x <- Q x

*/

void c8band_q_apply(int n, int bl, r8 *Q, r8 *x)
    {
     int i;
     int mq = bl + 1;

     for(i=n-1;i>=0;i--)
       {
        int i0, i1;  // true extent of column i of q
        int l, off;  // length of extent and offset into packed storage

        i0  = MAX(0, i);
        i1  = MIN(i+bl, n-1);
        l   = i1-i0+1;
        off = i0 - i;

        c8_apply_householder(l, Q+2*i*mq, x+2*i0);
       }

    }

/*

   Apply a unitary matrix in packed (banded) Householder form to a vector

      x <- Q' x

*/

void c8band_qh_apply(int n, int bl, r8 *Q, r8 *x)
    {
     int i;
     int mq = bl + 1;

     for(i=0;i<n-1;i++)
       {
        int i0, i1;  // true extent of column i of q
        int l, off;  // length of extent and offset into packed storage

        i0  = MAX(0, i);
        i1  = MIN(i+bl, n-1);
        l   = i1-i0+1;
        off = i0 - i;

        c8_apply_householder(l, Q+2*i*mq, x+2*i0);
       }
     }
