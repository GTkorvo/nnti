/*

   Routines for qr on banded complex matrices.

*/

#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmband.h"

//#include "chouse.h"

void c8_householder(int n, r8 *a, r8 *v);
void c8_apply_householder(int n, r8 *v, r8 *x);

/*

   Compute the QR factorization of a banded matrix. See band.h for the storage
   arrangement. This routine does an implicit complex shift of the input matrix,
   which is left unmodified.

   If the input matrix A has bandwidths [bl, bu] (zero for a diagonal matrix), then
   the Q and R factors will have bandwidths [bl+1,0] and [0,bl+bu+1]. Note that
   the Q factor is stored as the Householder vectors, not as a matrix (which would
   be almost full). Thus the allocation of Q is (bl+1)*n and R is (bl+bu+1)*n (complex).

     n      -- order of matrix A
     bl     -- lower bandwidth
     bu     -- upper bandwidth
     A      -- pointer to storage of A
     z      -- complex shift
     Q      -- pointer to storage of Q
     R      -- pointer to storage of R

*/
static void c8band_qr_shifted_general(int n, int bl, int bu, r8 *A, c8 z, r8 *Q, r8 *R);
static void c8band_qr_shifted_2(int n, r8 *A, c8 z, r8 *Q, r8 *R);

#define MAX_TMP 256

void c8band_qr_shifted(int n, int bl, int bu, r8 *A, c8 z, r8 *Q, r8 *R)
    {

     if((n < 8) || (bl != 2) || (bu != 2))
       c8band_qr_shifted_general(n, bl, bu, A, z, Q, R);
     else
       c8band_qr_shifted_2(n, A, z, Q, R);

    }

#ifdef DO_CLOCKS
int pcc_delta();
#endif

static void c8band_qr_shifted_general(int n, int bl, int bu, r8 *A, c8 z, r8 *Q, r8 *R)
    {
     int i,j,k;
     int m = bl + bu + 1;
     int mq = bl + 1;
     int mr = m;
     int mt = bl*2 + bu + 1;
     r8 *t;
     int kmax;
     r8 tmp[MAX_TMP];
     int ntmp;

// window size
     kmax = MIN(n,m); // but it will change.

// allocate the window.
     ntmp = mt*kmax*4;
     if(ntmp > MAX_TMP)
       t = cmm_alloc(sizeof(r8)*ntmp);
     else
       t = tmp;
//     t = cmm_alloc(sizeof(r8)*mt*kmax*2);
     r8_zero(ntmp, t);

// initialize the window
     for(k=0;k<kmax;k++)
       {
        int i0,i1,l,off;
        r8 *tp;

        i0 = MAX(0, k - bu);
        i1 = MIN(k+bl, n-1);
        l = i1-i0+1;
//        off = MAX(0, bu-k);   // buffer offset
        off = i0 - (k-bu);

        tp = t + 2*(bl+off+k*mt);

         {
          r8 *ap = A + 2*(off+m*k);
          for(i=0;i<l;i++)
            {
             tp[0] = ap[0];
             tp[1] = ap[1];
             tp += 2;
             ap += 2;
            }
         }

     // do the shift on the diagonal
        tp = t + 2*(bl+bu+k*mt);
        tp[0] += z.r;
        tp[1] += z.i;
       }

// note that the process includes j=n, although this is trivial
     for(j=0;j<n;j++)
       {
        int lh;
        r8 *vh; // pointer to householder vector

         // build the householder vector from t0 and copy to Q.
         lh = MIN(bu+1, n-j);
         vh = Q + 2*j*mq;

         c8_householder(lh, t+2*(bl+bu), vh);

         // apply to t0 and copy upper part to R
         c8_apply_householder(lh, vh, t+2*(bl+bu));

         c8_copy(bl+bu+1, (c8 *) t, 1, (c8 *) (R + 2*j*mr), 1);

         // apply to remainder, compressing the window.
         for(k=1;k<kmax;k++)
           {
            c8_apply_householder(lh, vh, t+2*(bl+bu-k + k*mt));
            c8_copy(mt, (c8 *) (t+2*k*mt), 1, (c8 *)(t+2*(k-1)*mt), 1);
           }

         kmax = kmax - 1;

         // get another column of a and append to the window
         if(kmax + j + 1< n)
           {
            int i0,i1,l,off;
            r8 *tp;

            k = kmax + j + 1;

            i0 = MAX(0, k - bu);
            i1 = MIN(k+bl, n-1);
            l = i1-i0+1;
//             off = MAX(0, bu-k);
            off = i0 - (k-bu);

            r8_zero(2*mt, (t + 2*kmax*mt));

            tp = t + 2*(bl+off + kmax*mt);
              {
               r8 *ap = A + 2*(off+m*k);
               for(i=0;i<l;i++)
                 {
                  tp[0] = ap[0];
                  tp[1] = ap[1];
                  tp += 2;
                  ap += 2;
                 }
              }

     // do the shift
            tp = t + 2*(bl+bu+kmax*mt);
            tp[0] += z.r;
            tp[1] += z.i;

            kmax = kmax + 1;

           }
         // (if we cant then reduce kmax)

       }

     if(ntmp > MAX_TMP)
       cmm_free(t);

    }

/*

   For this routine we are guaranteed that n >= 8 and bl = bu = 2
*/
static void c8band_qr_shifted_2(int n, r8 *A, c8 z, r8 *Q, r8 *R)
    {
     int bl = 2;
     int bu = 2;
     int i,j,k;
     int m  = 5;
     int mq = 3;
     int mr = 5;
     int mt = 7;
     r8 *t;
     int kmax;
     r8 tmp[MAX_TMP];
     int ntmp;
     int iclock, jclock[2];

// window size
     kmax = 5;

// allocate the window.
     ntmp = mt*kmax*2*2;
     if(ntmp > MAX_TMP)
       t = cmm_alloc(sizeof(r8)*ntmp);
     else
       t = tmp;
//     t = cmm_alloc(sizeof(r8)*mt*kmax*2);
     r8_zero(ntmp, t);

// initialize the window
     for(k=0;k<kmax;k++)
       {
        int i0,i1,l,off;
        r8 *tp;

        i0 = MAX(0, k - bu);
        i1 = MIN(k+bl, n-1);
        l = i1-i0+1;

        off = i0 - (k-bu);

        tp = t + 2*(bl+off+k*mt);

         {
          r8 *ap = A + 2*(off+m*k);
          for(i=0;i<l;i++)
            {
             tp[0] = ap[0];
             tp[1] = ap[1];
             tp += 2;
             ap += 2;
            }
         }

     // do the shift on the diagonal
        tp = t + 2*(bl+bu+k*mt);
        tp[0] += z.r;
        tp[1] += z.i;
       }

// note that the process includes j=n, although this is trivial

     for(j=0;j<2;j++)
       {
        int lh;
        r8 *vh; // pointer to householder vector

         // build the householder vector from t0 and copy to Q.
         lh = MIN(bu+1, n-j);
         vh = Q + 2*j*mq;

         c8_householder(lh, t+8, vh);

         // apply to t0 and copy upper part to R
         c8_apply_householder(lh, vh, t+8);

         c8_copy(bl+bu+1, (c8 *) t, 1, (c8 *) (R + 2*j*mr), 1);

         // apply to remainder, compressing the window.
         for(k=1;k<kmax;k++)
           {
            c8_apply_householder(lh, vh, t+2*(bl+bu-k + k*mt));
            c8_copy(mt, (c8 *) (t+2*k*mt), 1, (c8 *)(t+2*(k-1)*mt), 1);
           }

         kmax = kmax - 1;

         // get another column of a and append to the window
         if(kmax + j + 1< n) // j < n-5
           {
            int i0,i1,l,off;
            r8 *tp;

            k = kmax + j + 1;

            i0 = MAX(0, k - bu);
            i1 = MIN(k+bl, n-1);
            l = i1-i0+1;
//             off = MAX(0, bu-k);
            off = i0 - (k-bu);

            r8_zero(2*mt, (t + 2*kmax*mt));

            tp = t + 2*(bl+off + kmax*mt);
              {
               r8 *ap = A + 2*(off+m*k);
               for(i=0;i<l;i++)
                 {
                  tp[0] = ap[0];
                  tp[1] = ap[1];
                  tp += 2;
                  ap += 2;
                 }
              }

     // do the shift
            tp = t + 2*(bl+bu+kmax*mt);
            tp[0] += z.r;
            tp[1] += z.i;

            kmax = kmax + 1;

           }
         // (if we cant then reduce kmax)

       }

/*

   In this loop everything should be constant, no special cases,
   so we can unroll everything.

     int bl = 2;
     int bu = 2;
     int m  = 5;
     int mq = 3;
     int mr = 5;
     int mt = 7;

*/
#ifdef DO_CLOCKS
     iclock = pcc_delta();
     jclock[0] = 0;
     jclock[1] = 0;
#endif

    {
     r8 *vh = Q + 12; // pointer to householder vector
     r8 *rp = R + 20; // pointer to R matrix

     for(j=2;j<n-5;j++)
       {

      // build the householder vector from t0 and copy to Q.
#ifdef DO_CLOCKS
//         iclock = pcc_delta();
#endif

// TOFF = 8
//         c8_householder(3, t+8, vh);
#define TOFF 8
         {
          r8 vr, vi, vn, v1, zn, tn, ts;

          vr = t[TOFF + 0];
          vi = t[TOFF + 1];
          vn = (vr*vr + vi*vi);
          zn = t[TOFF+2]*t[TOFF+2] + t[TOFF+3]*t[TOFF+3]
             + t[TOFF+4]*t[TOFF+4] + t[TOFF+5]*t[TOFF+5];

          tn = sqrt(vn + zn);

          if(tn == 0.) // The Householder vector vanishes.
            {
             vr = 0.;
             vi = 0.;
             t[TOFF+0] = 0.;
             t[TOFF+1] = 0.;
             vh[0] = 0.;
             vh[1] = 0.;
             vh[2] = 0.;
             vh[3] = 0.;
             vh[4] = 0.;
             vh[5] = 0.;
            }
          else
            {
             v1 = sqrt(vn);

             if(v1 != 0.)
              {
               vr = vr/v1;
               vi = vi/v1;
              }
             else
              {
               vr = 1.;
               vi = 0.;
              }

         // compute modified first element
             t[TOFF+0] = -tn*vr;
             t[TOFF+1] = -tn*vi;
  
             tn += v1;

       // save the scaled Householder vector
             ts = 1./sqrt(tn*tn + zn);
             vh[0] = ts*tn*vr;
             vh[1] = ts*tn*vi;
             vh[2] = ts*t[TOFF+2]; 
             vh[3] = ts*t[TOFF+3];  
             vh[4] = ts*t[TOFF+4];  
             vh[5] = ts*t[TOFF+5];  
            }

         }

#ifdef DO_CLOCKS
//         jclock[0] += pcc_delta();
#endif

      // apply to t0
/* skip this, already applied (52 flops!)
         {register r8 tr, ti;
          tr  = t[8]*vh[0]  + t[9]*vh[1];
          ti  = t[9]*vh[0]  - t[8]*vh[1];
          tr += t[10]*vh[2] + t[11]*vh[3];
          ti += t[11]*vh[2] - t[10]*vh[3];
          tr += t[12]*vh[4] + t[13]*vh[5];
          ti += t[13]*vh[4] - t[12]*vh[5];
          tr *= -2.;
          ti *= -2.;
          t[8]  += vh[0]*tr - vh[1]*ti;
          t[9]  += vh[1]*tr + vh[0]*ti;
          t[10] += vh[2]*tr - vh[3]*ti;
          t[11] += vh[3]*tr + vh[2]*ti;
          t[12] += vh[4]*tr - vh[5]*ti;
          t[13] += vh[5]*tr + vh[4]*ti;
         }
*/
     //  copy nonzero part of t to R
         rp[0] = t[0];
         rp[1] = t[1];
         rp[2] = t[2];
         rp[3] = t[3];
         rp[4] = t[4];
         rp[5] = t[5];
         rp[6] = t[6];
         rp[7] = t[7];
         rp[8] = t[8];
         rp[9] = t[9];

      // apply to remainder, compressing the window.

         {register r8 tr, ti;                                
          tr  = t[20+0]*vh[0] + t[20+1]*vh[1];    
          ti  = t[20+1]*vh[0] - t[20+0]*vh[1];    
          tr += t[20+2]*vh[2] + t[20+3]*vh[3];       
          ti += t[20+3]*vh[2] - t[20+2]*vh[3];       
          tr += t[20+4]*vh[4] + t[20+5]*vh[5];       
          ti += t[20+5]*vh[4] - t[20+4]*vh[5];       
          tr *= -2.;       
          ti *= -2.;       
          t[20-20] = t[20-6];       
          t[20-19] = t[20-5];       
          t[20-18] = t[20-4];       
          t[20-17] = t[20-3];       
          t[20-16] = t[20-2];       
          t[20-15] = t[20-1];       
          t[20-14] = t[20+0] + vh[0]*tr - vh[1]*ti;       
          t[20-13] = t[20+1] + vh[1]*tr + vh[0]*ti;       
          t[20-12] = t[20+2] + vh[2]*tr - vh[3]*ti;       
          t[20-11] = t[20+3] + vh[3]*tr + vh[2]*ti;       
          t[20-10] = t[20+4] + vh[4]*tr - vh[5]*ti;       
          t[20- 9] = t[20+5] + vh[5]*tr + vh[4]*ti;       
          t[20- 8] = t[20+6]; 
          t[20- 7] = t[20+7]; 
         }

         {register r8 tr, ti;                                
          tr  = t[32+0]*vh[0] + t[32+1]*vh[1];    
          ti  = t[32+1]*vh[0] - t[32+0]*vh[1];    
          tr += t[32+2]*vh[2] + t[32+3]*vh[3];       
          ti += t[32+3]*vh[2] - t[32+2]*vh[3];       
          tr += t[32+4]*vh[4] + t[32+5]*vh[5];       
          ti += t[32+5]*vh[4] - t[32+4]*vh[5];       
          tr *= -2.;       
          ti *= -2.;       
          t[32-18] = t[32-4];       
          t[32-17] = t[32-3];       
          t[32-16] = t[32-2];       
          t[32-15] = t[32-1];       
          t[32-14] = t[32+0] + vh[0]*tr - vh[1]*ti;       
          t[32-13] = t[32+1] + vh[1]*tr + vh[0]*ti;       
          t[32-12] = t[32+2] + vh[2]*tr - vh[3]*ti;       
          t[32-11] = t[32+3] + vh[3]*tr + vh[2]*ti;       
          t[32-10] = t[32+4] + vh[4]*tr - vh[5]*ti;       
          t[32- 9] = t[32+5] + vh[5]*tr + vh[4]*ti;       
          t[32- 8] = t[32+6]; 
          t[32- 7] = t[32+7]; 
          t[32- 6] = t[32+8];       
          t[32- 5] = t[32+9];       
         }
         {register r8 tr, ti;                                
          tr  = t[44+0]*vh[0] + t[44+1]*vh[1];    
          ti  = t[44+1]*vh[0] - t[44+0]*vh[1];    
          tr += t[44+2]*vh[2] + t[44+3]*vh[3];       
          ti += t[44+3]*vh[2] - t[44+2]*vh[3];       
          tr += t[44+4]*vh[4] + t[44+5]*vh[5];       
          ti += t[44+5]*vh[4] - t[44+4]*vh[5];       
          tr *= -2.;       
          ti *= -2.;       
          t[44-16] = t[44-2];       
          t[44-15] = t[44-1];       
          t[44-14] = t[44+0] + vh[0]*tr - vh[1]*ti;       
          t[44-13] = t[44+1] + vh[1]*tr + vh[0]*ti;       
          t[44-12] = t[44+2] + vh[2]*tr - vh[3]*ti;       
          t[44-11] = t[44+3] + vh[3]*tr + vh[2]*ti;       
          t[44-10] = t[44+4] + vh[4]*tr - vh[5]*ti;       
          t[44- 9] = t[44+5] + vh[5]*tr + vh[4]*ti;       
          t[44- 8] = t[44+6]; 
          t[44- 7] = t[44+7]; 
          t[44- 6] = t[44+8];       
          t[44- 5] = t[44+9];       
          t[44- 4] = t[44+10];       
          t[44- 3] = t[44+11];       
         }
         {register r8 tr, ti;                                
          tr  = t[56+0]*vh[0] + t[56+1]*vh[1];    
          ti  = t[56+1]*vh[0] - t[56+0]*vh[1];    
          tr += t[56+2]*vh[2] + t[56+3]*vh[3];       
          ti += t[56+3]*vh[2] - t[56+2]*vh[3];       
          tr += t[56+4]*vh[4] + t[56+5]*vh[5];       
          ti += t[56+5]*vh[4] - t[56+4]*vh[5];       
          tr *= -2.;       
          ti *= -2.;       
          t[56-14] = t[56+0] + vh[0]*tr - vh[1]*ti;       
          t[56-13] = t[56+1] + vh[1]*tr + vh[0]*ti;       
          t[56-12] = t[56+2] + vh[2]*tr - vh[3]*ti;       
          t[56-11] = t[56+3] + vh[3]*tr + vh[2]*ti;       
          t[56-10] = t[56+4] + vh[4]*tr - vh[5]*ti;       
          t[56- 9] = t[56+5] + vh[5]*tr + vh[4]*ti;       
          t[56- 8] = t[56+6]; 
          t[56- 7] = t[56+7]; 
          t[56- 6] = t[56+8];       
          t[56- 5] = t[56+9];       
          t[56- 4] = t[56+10];       
          t[56- 3] = t[56+11];       
          t[56- 2] = t[56+12];       
          t[56- 1] = t[56+13];       
         }

         {  // get another column of a and append to the window
          r8 *ap = A + 10*(j+5);

          t[56] = 0.;
          t[57] = 0.;
          t[58] = 0.;
          t[59] = 0.;
          t[60] = ap[0];
          t[61] = ap[1];
          t[62] = ap[2];
          t[63] = ap[3];
          t[64] = ap[4]+z.r;
          t[65] = ap[5]+z.i;
          t[66] = ap[6];
          t[67] = ap[7];
          t[68] = ap[8];
          t[69] = ap[9];
         }

        vh += 6;
        rp += 10;

#ifdef DO_CLOCKS
//         jclock[1] += pcc_delta();
#endif

       } // end inner loop over j

     }

#ifdef DO_CLOCKS
//     iclock = pcc_delta();
//     printf("inner clocks %d %d\n",jclock[0],jclock[1]);
#endif

     for(j=n-5;j<n;j++)
       {
        int lh;
        r8 *vh; // pointer to householder vector

         // build the householder vector from t0 and copy to Q.
         lh = MIN(bu+1, n-j);
         vh = Q + 2*j*mq;

         c8_householder(lh, t+8, vh);

         // apply to t0 and copy upper part to R
         c8_apply_householder(lh, vh, t+8);

         c8_copy(bl+bu+1, (c8 *) t, 1, (c8 *) (R + 2*j*mr), 1);

         // apply to remainder, compressing the window.
         for(k=1;k<kmax;k++)
           {
            c8_apply_householder(lh, vh, t+2*(bl+bu-k + k*mt));
            c8_copy(mt, (c8 *) (t+2*k*mt), 1, (c8 *)(t+2*(k-1)*mt), 1);
           }

         kmax = kmax - 1;

         // get another column of a and append to the window
         if(kmax + j + 1< n) // j < n-5
           {
            int i0,i1,l,off;
            r8 *tp;

            k = kmax + j + 1;

            i0 = MAX(0, k - bu);
            i1 = MIN(k+bl, n-1);
            l = i1-i0+1;
//             off = MAX(0, bu-k);
            off = i0 - (k-bu);

            r8_zero(2*mt, (t + 2*kmax*mt));

            tp = t + 2*(bl+off + kmax*mt);
              {
               r8 *ap = A + 2*(off+m*k);
               for(i=0;i<l;i++)
                 {
                  tp[0] = ap[0];
                  tp[1] = ap[1];
                  tp += 2;
                  ap += 2;
                 }
              }

     // do the shift
            tp = t + 2*(bl+bu+kmax*mt);
            tp[0] += z.r;
            tp[1] += z.i;

            kmax = kmax + 1;

           }
         // (if we cant then reduce kmax)

       }

     if(ntmp > MAX_TMP)
       cmm_free(t);
    }

/*
void c8_householder(int n, r8 *a, r8 *v)
    {
     r8 v1r,v1i,v1n;
     c8 zn; r8 vn;
     r8 h0, h1, t;
     int i;

     if(n == 1)
       {
        v[0] = 0.;
        v[1] = 0.;
        return;
       }

     if(a != v)
       c8_copy(n, (c8 *) a, 1, (c8 *) v, 1);

     v1r = v[0];
     v1i = v[1];
     v1n = sqrt(v1r*v1r + v1i*v1i);
     zn = c8_hdot(n-1, (c8 *) (v+2) ,(c8 *) (v+2));
     vn = zn.r;

     t = sqrt(v1n*v1n + vn);

     if(t == 0.)
       {
        v[0] = 0.;
        v[1] = 0.;
        return;
       }

     if(v1n != 0.)
       {
        v1r = v1r/v1n;
        v1i = v1i/v1n;
       }
     else
       {
        v1r = 1.;
        v1i = 0.;
       }

     t = v1n+t;
     v[0] = v1r*t;
     v[1] = v1i*t;

     t = sqrt(t*t + vn);

     r8_scale(2*n, 1./t, v);

    }

void c8_apply_householder(int n, r8 *v, r8 *x)
    {
     c8 t;
     int i;

     if(n == 0)
      return;

     t = c8_hdot(n, (c8 *) x, (c8 *) v);
     t.r = -2.*t.r;
     t.i = -2.*t.i;

     c8_axpy(n, t, (c8 *)v, (c8 *)x);

    }
*/

/*

    The construction requires a complex dot product and a real 
  scaling of the complex vector, so E(n) is 10n flops.  

    For a matrix of size n we build n vectors of average length 
  n/2. For n=1024 we have E = 10n^2/2 = 50MFlops. Laptop (800Mhz)
  does this in .1 sec, so 500MFlops, not too shabby.
 */

void c8_householder(int n, r8 *a, r8 *v)
   {
     r8   v1r,v1i,v1n;
     c8 zn;
     r8 vn;
     r8 t;

     if(n <= 0)
      return;

#ifdef HOUSE_IGNORE_1
     if(n == 1)
       {
        v[0] = 0.;
        v[1] = 0.;
        return;
       }
#endif

     if(a != v)
       c8_copy(n, (c8 *) a, 1, (c8 *) v, 1);

     v1r = v[0];
     v1i = v[1];
     v1n = sqrt(v1r*v1r + v1i*v1i);
     zn = c8_hdot(n-1, (c8 *) (v+2) ,(c8 *) (v+2));
     vn = zn.r;

     // printf("v1r %le v1i %le vn %le\n",v1r, v1i, vn);

     t = sqrt(v1n*v1n + vn);

     if(t == 0.)
       {
        v[0] = 0.;
        v[1] = 0.;
        return;
       }

     if(v1n != 0.)
       {
        v1r = v1r/v1n;
        v1i = v1i/v1n;
       }
     else
       {
        v1r = 1.;
        v1i = 0.;
       }

     t = v1n+t;
     v[0] = v1r*t;
     v[1] = v1i*t;

     t = sqrt(t*t + vn);

     r8_scale(2*n, 1./t, v);

    }

void c8_apply_householder(int n, r8 *v, r8 *x)
    {
     c8 t;

     if(n == 0)
      return;

     t = c8_hdot(n, (c8 *) x, (c8 *) v);  // this is x.v'

     t.r = -2.*t.r;
     t.i = -2.*t.i;

     c8_axpy(n, t, (c8 *)v, (c8 *)x);

    }

