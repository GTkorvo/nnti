#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmband.h"

/*

   Compute the number of eigenvalues of a banded Hermitian matrix A greater than x.

   We do this by computing the number of eigenvalues of A-xI greater than zero.

*/

int c8band_oracle(int n, int b, r8 *A, r8 x, r8 *tmp);

int oracle_count;

#define DO_PHASE 1

static void shiftrow(int b, r8 *t, int s)  // shift complex vector t s places left and pad with zeros
    {
     int i;
     int l1,l2;
     r8 *v;

     l2 = (2*b+1);
     l1 = l2-s;      // if s == b then this is b+1

     v = t + 2*s;    // t + 2*b
     for(i=0;i<l1;i++)
       {
        t[0] = v[0];
        t[1] = v[1];
        t += 2;
        v += 2;
       }

     for(i=l1;i<l2;i++)
       {
        t[0] = 0.;
        t[1] = 0.;
        t += 2;
       }
    }

static c8 compute_phase(r8 zr, r8 zi)  // compute phase of (zr, zi)
    {
     c8 phase;
#ifdef DO_PHASE
     r8 av = zr*zr+zi*zi;

     if(av == 0.)
       {
        phase.r = 1.;
        phase.i = 0.;
       }
     else
       {
        av = 1./sqrt(av);
        phase.r = zr*av;
        phase.i = zi*av;
       }

#else
     phase.r = zr;
     phase.i = zi;
#endif
     return phase;
    }

static c8 c8_mult(c8 za, c8 zb)
    {
     c8 z;
     z.r = za.r*zb.r - za.i*zb.i;
     z.i = za.r*zb.i + za.i*zb.r;
     return z;
    }


static int c8band_oracle_general(int n, int b, r8 *A, r8 x, r8 *tmp);
static int c8band_oracle_1(int n, int b, r8 *A, r8 x, r8 *tmp);
static int c8band_oracle_2(int n, int b, r8 *A, r8 x, r8 *tmp);
static int c8band_oracle_2x(int n, int b, r8 *A, r8 x, r8 *tmp);
static int c8band_oracle_4(int n, int b, r8 *A, r8 x, r8 *tmp);

int c8band_oracle(int n, int b, r8 *A, r8 x, r8 *tmp)
    {
     int k;

     oracle_count++;

     switch(b)
       {
      case 1:
         k =  c8band_oracle_1(n, b, A, x, tmp);
         break;
      case 2:
         k =  c8band_oracle_2x(n, b, A, x, tmp);
//         k = c8band_oracle_general(n, b, A, x, tmp);
         break;
      case 4:
         k =  c8band_oracle_4(n, b, A, x, tmp);
         break;
      default:
         k =  c8band_oracle_general(n, b, A, x, tmp);
         break;
       }

//     printf("oracle %lf count %d\n", x, k);
     return k;

    }

#define B_LOOP1 1

static int c8band_oracle_1(int n, int bb, r8 *A, r8 x, r8 *tmp)
    {
     int ks;
     r8 ws;
     int i,j,k;
     int rowsize;
     r8 *t;
     c8 zphase; // current phase product
     c8 sphase; // saved phase product

     rowsize = 6;

     r8_copy(rowsize*n, A, 1, tmp, 1);  // Make a temp copy

     t = tmp;
     for(k=0;k<n;k++)  // shift by -xI
       {
        t[2*B_LOOP1] -= x;
        t += rowsize;
       }

     shiftrow(B_LOOP1, tmp, B_LOOP1);

     ks = 1;

     if(tmp[0] < 0.)
        ks = 0;

     zphase = compute_phase(tmp[0], 0.);

     ws = 1.;  // accumulate sign due to swaps in each k step

     for(k=1;k<n;k++)
       {
        int l0;
        r8 *tk = tmp + rowsize*k;

        l0 = k-B_LOOP1;
        if(l0 < 0)
          {
           shiftrow(B_LOOP1, tk, -l0);
           l0 = 0;
          }

        for(i=l0;i<k;i++)
          {
           r8 *ti = tmp + rowsize*i;
           r8 ar, ai,ak;
           c8 tv;

           ai = ti[0]*ti[0] + ti[1]*ti[1];
           ak = tk[0]*tk[0] + tk[1]*tk[1];

           if(ak > ai)
             {
              for(j=0;j<rowsize;j++)
                {
                 r8 tx;
                 tx = ti[j]; ti[j] = tk[j]; tk[j] = tx;
                }
              ws = -ws;
             }

           ar = ti[0];
           ai = ti[1];
           ak = -1./(ar*ar + ai*ai);

           tv.r = ak*(tk[0]*ar + tk[1]*ai);
           tv.i = ak*(tk[1]*ar - tk[0]*ai);

           {
            int jj;
            r8 *xp,*yp;

            yp = tk;
            xp = ti+2;

            for(jj=1;jj<2*B_LOOP1+1;jj++)
              {
               yp[0] = yp[2] + tv.r*xp[0] - tv.i*xp[1];
               yp[1] = yp[3] + tv.i*xp[0] + tv.r*xp[1];
               yp += 2;
               xp += 2;
              }

            yp[0] = 0.;
            yp[1] = 0.;
           }

          }

      {
         c8 yphase;
         int k1;

         k1 = k-B_LOOP1;
         if(k1 <= 0)
           {
            k1 = 0;
            sphase.r = 1.;  // dont use saved phase yet.
            sphase.i = 0.;
           }

         {  // Saving the initial phase product converts the effort from O(n) to O(b).
          c8 tphase;
          r8 *tp = tmp + rowsize*k1;
          tphase = compute_phase(tp[0],tp[1]);
          sphase = c8_mult(sphase, tphase);
         }

         yphase = sphase;

         for(j=k1+1;j<=k;j++)
           {
            c8 tphase;
            r8 *tp = tmp + rowsize*j;

            tphase.r = tp[0];
            tphase.i = tp[1];
            yphase = c8_mult(yphase, tphase);            
           }

         yphase.r *= ws;
         yphase.i *= ws;

// count agreements in sign
         if(     (yphase.r >= 0. && zphase.r >= 0.)
              || (yphase.r < 0. && zphase.r < 0.))
           ks++;

// save computed phase for next comparison
         zphase = yphase;

         }

       }

     return ks;
    }

#define B_LOOP2 2

static int clocks[8];
static void c8band_oracle_2_init_clocks()
    {
     int i;
     for(i=0;i<8;i++)
       clocks[i] = 0;
    }
static int *c8band_oracle_2_clocks()
    {
     return clocks;
    }

static int c8band_oracle_2(int n, int bb, r8 *A, r8 x, r8 *tmp)
    {
     int ks;
     r8 ws;
     int i,j,k;
     int iclock, iflops;
     int rowsize;
     r8 *t;
     c8 zphase; // current phase product
     c8 sphase; // saved phase product

     c8band_oracle_2_init_clocks();
#ifdef DO_CLOCKS
     i = pcc_delta();
#endif
     iflops = 0;

     rowsize = 10; // 2*(2*b+1)

     r8_copy(10*n, A, 1, tmp, 1);  // Make a temp copy of A

     t = tmp;
     for(k=0;k<n;k++)  // shift by -xI
       {
        t[2*B_LOOP2] -= x;
        t += 10;
       }

     shiftrow(B_LOOP2, tmp, B_LOOP2);

     ks = 1;

     if(tmp[0] < 0.)
        ks = 0;

     zphase = compute_phase(tmp[0], 0.);

     ws = 1.;  // accumulate sign due to swaps in each k step

#ifdef DO_CLOCKS
     clocks[0] = pcc_delta();
#endif
     for(k=1;k<n;k++)
       {
        int l0;
        r8 *tk = tmp + 10*k;
        r8 *ti;

        l0 = k-B_LOOP2;
        if(l0 < 0)
          {
           shiftrow(B_LOOP2, tk, -l0);
           l0 = 0;
          }

#ifdef DO_CLOCKS
        clocks[1] += pcc_delta();
#endif
        ti = tmp + 10*l0;
        for(i=l0;i<k;i++)
          {
           r8 ar, ai,ak;

           ai = ti[0]*ti[0] + ti[1]*ti[1];
           ak = tk[0]*tk[0] + tk[1]*tk[1];

           if(ak > ai)
             {
              r8 txr, txi;

              txr = ti[0];
              txi = ti[1];
              ti[0] = tk[0];
              ti[1] = tk[1];
              tk[0] = txr;
              tk[1] = txi;

              txr = ti[2+0];
              txi = ti[2+1];
              ti[2+0] = tk[2+0];
              ti[2+1] = tk[2+1];
              tk[2+0] = txr;
              tk[2+1] = txi;

              txr = ti[4+0];
              txi = ti[4+1];
              ti[4+0] = tk[4+0];
              ti[4+1] = tk[4+1];
              tk[4+0] = txr;
              tk[4+1] = txi;

              txr = ti[6+0];
              txi = ti[6+1];
              ti[6+0] = tk[6+0];
              ti[6+1] = tk[6+1];
              tk[6+0] = txr;
              tk[6+1] = txi;

              txr = ti[8+0];
              txi = ti[8+1];
              ti[8+0] = tk[8+0];
              ti[8+1] = tk[8+1];
              tk[8+0] = txr;
              tk[8+1] = txi;

/*
              for(j=0;j<10;j++)
                {
                 r8 tx;
                 tx = ti[j]; ti[j] = tk[j]; tk[j] = tx;
                }
*/
              ws = -ws;
             }

           ar = ti[0];
           ai = ti[1];

           ak = 1./(ar*ar + ai*ai);
           {
            r8 xr, xi;
            xr = ak*(tk[0]*ar + tk[1]*ai);
            xi = ak*(tk[1]*ar - tk[0]*ai);

//            printf("k %d i %d xr %lf xi %lf\n",k,k-i,xr,xi);

            tk[0] = tk[2] - xr*ti[2+0] + xi*ti[2+1];
            tk[1] = tk[3] - xi*ti[2+0] - xr*ti[2+1];

            tk[2] = tk[4] - xr*ti[2+2] + xi*ti[2+3];
            tk[3] = tk[5] - xi*ti[2+2] - xr*ti[2+3];

            tk[4] = tk[6] - xr*ti[2+4] + xi*ti[2+5];
            tk[5] = tk[7] - xi*ti[2+4] - xr*ti[2+5];

            tk[6] = tk[8] - xr*ti[2+6] + xi*ti[2+7];
            tk[7] = tk[9] - xi*ti[2+6] - xr*ti[2+7];

            tk[8] = 0.;
            tk[9] = 0.;
          }

           iflops += 8 + 6 + 6 + 32;

          ti += 10;
          } // end loop over i

#ifdef DO_CLOCKS
        clocks[2] += pcc_delta();
#endif

      {
         c8 yphase;
         int k1;

         k1 = k-B_LOOP2;
         if(k1 <= 0)
           {
            k1 = 0;
            sphase.r = 1.;  // dont use saved phase yet.
            sphase.i = 0.;
           }

/*
         {  // Saving the initial phase product converts the effort from O(n) to O(b).
          c8 tphase;
          r8 *tp = tmp + 10*k1;
          tphase = compute_phase(tp[0],tp[1]);
          sphase = c8_mult(sphase, tphase);
         }
*/

         {
          r8 tr, ti, tt;
          r8 sr, si;
          r8 *tp = tmp + 10*k1;
          tr = tp[0];
          ti = tp[1];
          tt = tr*tr + ti*ti;
          tt = sqrt(tt);
          if(tt > 0.)
            {
             tt = 1./tt;
             sr = sphase.r;
             si = sphase.i;
             sphase.r = tt*(sr*tr - si*ti);
             sphase.i = tt*(sr*ti + si*tr);
            }
         }

//         printf("k %d sphase %lf %lf\n",k,sphase.r, sphase.i);

         yphase = sphase;


        {
         r8 *tp = tmp + 10*(k1+1);
         r8 xr, xi;
         int nj = k-k1;
         int jj;
         xr = yphase.r;
         xi = yphase.i;
//         for(j=k1+1;j<=k;j++)
         for(jj=nj;jj>0;jj--)
           {
            r8 tr, ti, xt;
            tr = tp[0];
            ti = tp[1]; tp += 10;
            xt = xr;
            xr = xt * tr - xi * ti;
            xi = xt * ti + xi * tr;
           }
         yphase.r = xr*ws;
         yphase.i = xi*ws;
        }

//       printf("k %d phase %lf %lf ks %d ws %lf\n",k,yphase.r,yphase.i,ks, ws);

// count agreements in sign
         if(     (yphase.r >= 0. && zphase.r >= 0.)
              || (yphase.r < 0. && zphase.r < 0.))
           ks++;

// save computed phase for next comparison
         zphase = yphase;

         }

#ifdef DO_CLOCKS
        clocks[3] += pcc_delta();
#endif
       }

/*
     printf("oracle %d c0 %d c1 %d c2 %d c3 %d flops %d\n",n, clocks[0], clocks[1], clocks[2],clocks[3],iflops);

*/

//     printf("oracle x %lf ks %d ws %lf\n",x,ks,ws);

     return ks;
    }

static int c8band_oracle_2_saved(int n, int bb, r8 *A, r8 x, r8 *tmp)
    {
     int ks;
     r8 ws;
     int i,j,k;
     int rowsize;
     r8 *t;
     c8 zphase; // current phase product
     c8 sphase; // saved phase product

     rowsize = 10; // 2*(2*b+1)

     r8_copy(10*n, A, 1, tmp, 1);  // Make a temp copy

     t = tmp;
     for(k=0;k<n;k++)  // shift by -xI
       {
        t[2*B_LOOP2] -= x;
        t += 10;
       }

     shiftrow(B_LOOP2, tmp, B_LOOP2);

     ks = 1;

     if(tmp[0] < 0.)
        ks = 0;

     zphase = compute_phase(tmp[0], 0.);

     ws = 1.;  // accumulate sign due to swaps in each k step

     for(k=1;k<n;k++)
       {
        int l0;
        r8 *tk = tmp + 10*k;

        l0 = k-B_LOOP2;
        if(l0 < 0)
          {
           shiftrow(B_LOOP2, tk, -l0);
           l0 = 0;
          }

        for(i=l0;i<k;i++)
          {
           r8 *ti = tmp + 10*i;
           r8 ar, ai,ak;
           c8 tv;

           ai = ti[0]*ti[0] + ti[1]*ti[1];
           ak = tk[0]*tk[0] + tk[1]*tk[1];

           if(ak > ai)
             {
              for(j=0;j<10;j++)
                {
                 r8 tx;
                 tx = ti[j]; ti[j] = tk[j]; tk[j] = tx;
                }
              ws = -ws;
             }

           ar = ti[0];
           ai = ti[1];

           ak = 1./(ar*ar + ai*ai);
           tv.r = ak*(tk[0]*ar + tk[1]*ai);
           tv.i = ak*(tk[1]*ar - tk[0]*ai);

/*
           {
            int jj;
            r8 *xp,*yp;

            yp = tk;
            xp = ti+2;

            for(jj=1;jj<2*B_LOOP2+1;jj++)
              {
               r8 xr,xi;

               yp[0] = yp[2] - tv.r*xp[0] + tv.i*xp[1];
               yp[1] = yp[3] - tv.i*xp[0] - tv.r*xp[1];
               yp += 2;
               xp += 2;
              }

            yp[0] = 0.;
            yp[1] = 0.;
           }
*/
           {
            register r8 *xp;
            register r8 *yp;
            register r8 xr, xi;

            xr = tv.r;
            xi = tv.i;

            yp = tk;
            xp = ti+2;

            yp[0] = yp[2] - xr*xp[0] + xi*xp[1];
            yp[1] = yp[3] - xi*xp[0] - xr*xp[1];

            yp[2] = yp[4] - xr*xp[2] + xi*xp[3];
            yp[3] = yp[5] - xi*xp[2] - xr*xp[3];

            yp[4] = yp[6] - xr*xp[4] + xi*xp[5];
            yp[5] = yp[7] - xi*xp[4] - xr*xp[5];

            yp[6] = yp[8] - xr*xp[6] + xi*xp[7];
            yp[7] = yp[9] - xi*xp[6] - xr*xp[7];

            yp[8] = 0.;
            yp[9] = 0.;
          }

          } // end loop over i


      {
         c8 yphase;
         int k1;

         k1 = k-B_LOOP2;
         if(k1 <= 0)
           {
            k1 = 0;
            sphase.r = 1.;  // dont use saved phase yet.
            sphase.i = 0.;
           }

/*
         {  // Saving the initial phase product converts the effort from O(n) to O(b).
          c8 tphase;
          r8 *tp = tmp + 10*k1;
          tphase = compute_phase(tp[0],tp[1]);
          sphase = c8_mult(sphase, tphase);
         }
*/

         {
          r8 tr, ti, tt;
          r8 sr, si;
          r8 *tp = tmp + 10*k1;
          tr = tp[0];
          ti = tp[1];
          tt = tr*tr + ti*ti;
          tt = sqrt(tt);
          if(tt > 0.)
            {
             tt = 1./tt;
             sr = sphase.r;
             si = sphase.i;
             sphase.r = tt*(sr*tr - si*ti);
             sphase.i = tt*(sr*ti + si*tr);
            }
         }

         yphase = sphase;
/*
         for(j=k1+1;j<=k;j++)
           {
            c8 tphase;
            r8 *tp = tmp + 10*j;

            tphase.r = tp[0];
            tphase.i = tp[1];
            yphase = c8_mult(yphase, tphase);            
           }

         yphase.r *= ws;
         yphase.i *= ws;

*/

        {
         r8 *tp = tmp + 10*(k1+1);
         r8 xr, xi;
         int nj = k-k1;
         int jj;
         xr = yphase.r;
         xi = yphase.i;
//         for(j=k1+1;j<=k;j++)
         for(jj=nj;jj>0;jj--)
           {
            r8 tr, ti, xt;
            tr = tp[0];
            ti = tp[1]; tp += 10;
            xt = xr;
            xr = xt * tr - xi * ti;
            xi = xt * ti + xi * tr;
           }
         yphase.r = xr*ws;
         yphase.i = xi*ws;
        }

// count agreements in sign
         if(     (yphase.r >= 0. && zphase.r >= 0.)
              || (yphase.r < 0. && zphase.r < 0.))
           ks++;

// save computed phase for next comparison
         zphase = yphase;

         }

       }

     return ks;
    }

#define B_LOOP4 4
#define B_ROWSIZE4 (2*(2*B_LOOP4+1))

static int c8band_oracle_4(int n, int bb, r8 *A, r8 x, r8 *tmp)
    {
     int ks;
     r8 ws;
     int i,j,k;
//     int rowsize;
     r8 *t;
     c8 zphase; // current phase product
     c8 sphase; // saved phase product

//     rowsize = 2*(2*b+1);

     r8_copy(B_ROWSIZE4*n, A, 1, tmp, 1);  // Make a temp copy

     t = tmp;
     for(k=0;k<n;k++)  // shift by -xI
       {
        t[2*B_LOOP4] -= x;
        t += B_ROWSIZE4;
       }

     shiftrow(B_LOOP4, tmp, B_LOOP4);

     ks = 1;

     if(tmp[0] < 0.)
        ks = 0;

     zphase = compute_phase(tmp[0], 0.);

     ws = 1.;  // accumulate sign due to swaps in each k step

     for(k=1;k<n;k++)
       {
        int l0;
        r8 *tk = tmp + B_ROWSIZE4*k;

        l0 = k-B_LOOP4;
        if(l0 < 0)
          {
           shiftrow(B_LOOP4, tk, -l0);
           l0 = 0;
          }

        for(i=l0;i<k;i++)
          {
           r8 *ti = tmp + B_ROWSIZE4*i;
           r8 ar, ai,ak;
           c8 tv;

           ai = ti[0]*ti[0] + ti[1]*ti[1];
           ak = tk[0]*tk[0] + tk[1]*tk[1];

           if(ak > ai)
             {
              r8 zr,zi;
              // swap ti, tk

              zr = ti[0]; zi = ti[1]; ti[0] = tk[0]; ti[1] = tk[1]; tk[0] = zr; tk[1] = zi;
              zr = ti[2]; zi = ti[3]; ti[2] = tk[2]; ti[3] = tk[3]; tk[2] = zr; tk[3] = zi;
              zr = ti[4]; zi = ti[5]; ti[4] = tk[4]; ti[5] = tk[5]; tk[4] = zr; tk[5] = zi;
              zr = ti[6]; zi = ti[7]; ti[6] = tk[6]; ti[7] = tk[7]; tk[6] = zr; tk[7] = zi;
              zr = ti[8]; zi = ti[9]; ti[8] = tk[8]; ti[9] = tk[9]; tk[8] = zr; tk[9] = zi;
              zr = ti[10]; zi = ti[11]; ti[10] = tk[10]; ti[11] = tk[11]; tk[10] = zr; tk[11] = zi;
              zr = ti[12]; zi = ti[13]; ti[12] = tk[12]; ti[13] = tk[13]; tk[12] = zr; tk[13] = zi;
              zr = ti[14]; zi = ti[15]; ti[14] = tk[14]; ti[15] = tk[15]; tk[14] = zr; tk[15] = zi;
              zr = ti[16]; zi = ti[17]; ti[16] = tk[16]; ti[17] = tk[17]; tk[16] = zr; tk[17] = zi;

/*
              for(j=0;j<B_ROWSIZE4;j++)
                {
                 r8 tx;
                 tx = ti[j]; ti[j] = tk[j]; tk[j] = tx;
                }

*/
              ws = -ws;
             }

           ar = ti[0];
           ai = ti[1];
           ak = -1./(ar*ar + ai*ai);

           tv.r = ak*(tk[0]*ar + tk[1]*ai);
           tv.i = ak*(tk[1]*ar - tk[0]*ai);

           {
            c8 *xp, *yp;
            yp = (c8 *) tk;
            xp = (c8 *) ti;
            yp[0].r = yp[1].r + tv.r*xp[1].r - tv.i*xp[1].i;
            yp[0].i = yp[1].i + tv.i*xp[1].r + tv.r*xp[1].i;
            yp[1].r = yp[2].r + tv.r*xp[2].r - tv.i*xp[2].i;
            yp[1].i = yp[2].i + tv.i*xp[2].r + tv.r*xp[2].i;
            yp[2].r = yp[3].r + tv.r*xp[3].r - tv.i*xp[3].i;
            yp[2].i = yp[3].i + tv.i*xp[3].r + tv.r*xp[3].i;
            yp[3].r = yp[4].r + tv.r*xp[4].r - tv.i*xp[4].i;
            yp[3].i = yp[4].i + tv.i*xp[4].r + tv.r*xp[4].i;
            yp[4].r = yp[5].r + tv.r*xp[5].r - tv.i*xp[5].i;
            yp[4].i = yp[5].i + tv.i*xp[5].r + tv.r*xp[5].i;
            yp[5].r = yp[6].r + tv.r*xp[6].r - tv.i*xp[6].i;
            yp[5].i = yp[6].i + tv.i*xp[6].r + tv.r*xp[6].i;
            yp[6].r = yp[7].r + tv.r*xp[7].r - tv.i*xp[7].i;
            yp[6].i = yp[7].i + tv.i*xp[7].r + tv.r*xp[7].i;
            yp[7].r = yp[8].r + tv.r*xp[8].r - tv.i*xp[8].i;
            yp[7].i = yp[8].i + tv.i*xp[8].r + tv.r*xp[8].i;
            yp[8].r = 0.;
            yp[8].i = 0.;
            
           }


          }

      {
         c8 yphase;
         int k1;

         k1 = k-B_LOOP4;
         if(k1 <= 0)
           {
            k1 = 0;
            sphase.r = 1.;  // dont use saved phase yet.
            sphase.i = 0.;
           }

         {  // Saving the initial phase product converts the effort from O(n) to O(b).
          r8 tx,ta,zr,zi;
          r8 *tp = tmp + B_ROWSIZE4*k1;

//          tphase = compute_phase(tp[0],tp[1]);
//          sphase = c8_mult(sphase, tphase);

          zr = tp[0]; zi = tp[1];
          ta = 1./sqrt(zr*zr + zi*zi);
//          ta = 1.;

          tx = sphase.r;
          sphase.r = ta*(tx*zr - sphase.i*zi);
          sphase.i = ta*(tx*zi + sphase.i*zr);

         }

         yphase = sphase;

         for(j=k1+1;j<=k;j++)
           {
            c8 tphase; r8 tx;
            r8 *tp = tmp + B_ROWSIZE4*j;

            tphase.r = tp[0];
            tphase.i = tp[1];
            tx = yphase.r;
            yphase.r = tx*tphase.r - yphase.i*tphase.i;
            yphase.i = tx*tphase.i + yphase.i*tphase.r;
//            yphase = c8_mult(yphase, tphase);            
           }

         yphase.r *= ws;
         yphase.i *= ws;

// count agreements in sign
         if(     (yphase.r >= 0. && zphase.r >= 0.)
              || (yphase.r < 0. && zphase.r < 0.))
           ks++;

// save computed phase for next comparison
         zphase = yphase;

         }

       }

     return ks;
    }

static int c8band_oracle_general(int n, int b, r8 *A, r8 x, r8 *tmp)
    {
     int ks;
     r8 ws;
     int i,j,k;
     int rowsize;
     r8 *t;
     c8 zphase; // current phase product
     c8 sphase; // saved phase product
#ifdef DO_SLR
     r8 slr, sli;
#endif
//     int nflops;

     rowsize = 2*(2*b+1);

//     printf("copy %d\n",rowsize*n*sizeof(r8));
//     memcpy(tmp, A, rowsize*n*sizeof(r8));

     r8_copy(rowsize*n, A, 1, tmp, 1);  // Make a temp copy

//     nflops = 0;

     t = tmp;
     for(k=0;k<n;k++)  // shift by -xI
       {
        t[2*b] -= x;
        t += rowsize;
       }
//     nflops += n;

     shiftrow(b, tmp, b);

     ks = 1;

     if(tmp[0] < 0.)
        ks = 0;

#ifdef DO_SLR
     slr = tmp[0];
     sli = 0.;
     ks1 = ks;
#endif

     //printf("ks start %d\n",ks);

     zphase = compute_phase(tmp[0], 0.);

     ws = 1.;  // accumulate sign due to swaps in each k step

     for(k=1;k<n;k++)
       {
        int l0;
        r8 *tk = tmp + rowsize*k;

        l0 = k-b;
        if(l0 < 0)
          {
           shiftrow(b, tk, -l0);
           l0 = 0;
          }

        for(i=l0;i<k;i++)
          {
           r8 *ti = tmp + rowsize*i;
           r8 ar, ai,ak;
           c8 tv;

           ai = ti[0]*ti[0] + ti[1]*ti[1];
           ak = tk[0]*tk[0] + tk[1]*tk[1];

//           nflops += 6;

           if(ak > ai)
             {
              for(j=0;j<rowsize;j++)
                {
                 r8 tx;
                 tx = ti[j]; ti[j] = tk[j]; tk[j] = tx;
                }
              ws = -ws;
             }

           ar = ti[0];
           ai = ti[1];
           ak = -1./(ar*ar + ai*ai);

           tv.r = ak*(tk[0]*ar + tk[1]*ai);
           tv.i = ak*(tk[1]*ar - tk[0]*ai);

//           nflops += 4 + 5 + 6;

//           c8_axpy(2*b+1, tv, ti, tk);
//           shiftrow(b, tk, 1);

           {
            int jj,b21;
            r8 *xp,*yp;

            b21 = 2*b+1;
            yp = tk;
            xp = ti+2;

// This loop implies about 16 n b^2 flops per oracle
//
// To get all the eigenvalues for n=128, b=8 we have about 500MFlops
// which is 240 * 128^3, which is way too slow. If I had b = 2 then
// should get factor of 16 speedup. Might also hardwire b into the loop.
// In fact I only get a factor of 5, not 16

// (2*8+1)*8 = 140
// (2*2+1)*2 =  10
// this implies 7x, pretty close to 5x.

            for(jj=1;jj<b21;jj++)
              {
//               nflops += 8;
               yp[0] = yp[2] + tv.r*xp[0] - tv.i*xp[1];
               yp[1] = yp[3] + tv.i*xp[0] + tv.r*xp[1];
               yp += 2;
               xp += 2;
              }

            yp[0] = 0.;
            yp[1] = 0.;
           }

          }

#ifdef DO_SLR

      {
         r8 yr,yi,av;
         c8 yv;
         c8 tphase;
         int k1;

         k1 = k-b;
         if(k1 < 0)
            k1 = 0;

         yv.r = 1.;
         yv.i = 0.;

         yr = 1.;
         yi = 0.;

         for(j=0;j<=k;j++)
           {
            c8 cv,tv;

            t = tmp + rowsize*j;

            tv.r = t[0];
            tv.i = t[1];
            yv = c8_mult(yv, tv);

            tx = yr;
            yr = (tx*t[0] - yi*t[1]);
            yi = (tx*t[1] + yi*t[0]);

            if(j <= k1)
              {
               c8 tphase, yphase;
               yphase = compute_phase(yv.r,yv.i);
               tphase = compute_phase(t[0],t[1]);
               //printf("t phase %lf %lf y %lf %lf j %d\n",tphase.r, tphase.i, yphase.r, yphase.i, j);
              }

           }

         yr *= ws;
         yi *= ws;

// count agreements in sign


         if(yr >= 0. && slr >= 0.)
           ks1++;
         else if(yr < 0. && slr < 0.)
           ks1++;


         slr = yr;
         sli = yi;

        tphase = compute_phase(slr,sli);
        //printf("ks1 %d (%lf %lf) \n",ks1,tphase.r,tphase.i);

       }
#endif

/*
   The stabilized method computes the phase of the principal minor,
   which is all we really need. This should be numerically very
   stable since we are multiplying complex numbers of modulus 1.

   This version does O(n) work. We only need to do O(b) work, by
   saving the product up to k-b, since it cannot change.
*/
      {
         c8 yphase;
         int k1;

         k1 = k-b;
         if(k1 <= 0)
           {
            k1 = 0;
            sphase.r = 1.;  // dont use saved phase yet.
            sphase.i = 0.;
           }

//         nflops += 2;

         {  // Saving the initial phase product converts the effort from O(n) to O(b).
          c8 tphase;
          r8 *tp = tmp + rowsize*k1;
          tphase = compute_phase(tp[0],tp[1]);

//         nflops += 20;

          //printf("s phase %lf %lf\n",sphase.r,sphase.i);
          sphase = c8_mult(sphase, tphase);
          //printf("t  phase %lf %lf y phase %lf %lf\n",tphase.r,tphase.i,sphase.r,sphase.i);
         }

         yphase = sphase;
         //printf("c initial phase %lf %lf k1 %d\n",yphase.r,yphase.i,k1);

         for(j=k1+1;j<=k;j++)
           {
            c8 tphase;
            r8 *tp = tmp + rowsize*j;
//            tphase = compute_phase(tp[0],tp[1]);
            tphase.r = tp[0];
            tphase.i = tp[1];
            yphase = c8_mult(yphase, tphase);            

//            nflops += 8;
           }

         yphase.r *= ws;
         yphase.i *= ws;

//         nflops += 2;

// count agreements in sign
         if(     (yphase.r >= 0. && zphase.r >= 0.)
              || (yphase.r < 0. && zphase.r < 0.))
           ks++;

//         nflops += 4;

         //printf("ks %d (%lf %lf) (%lf %lf)\n\n",ks,yphase.r,yphase.i,zphase.r,zphase.i);

// save computed phase for next comparison
         zphase = yphase;

         }

       }

//     printf("nflops %d\n",nflops);

//     printf("x %lf ks %ld\n",x,ks);
     return ks;
    }

/*

   Fast version for b=2

   In this we use a window rather than copying A to
   a temporary shifted array, and everything is unrolled.

   This works and is about as fast as it will go. The
   end cases are moved out of the main loop so there
   are no extra tests, and all array references are
   simplified as much as possible.

*/
static int c8band_oracle_2x(int n, int bb, r8 *A, r8 x, r8 *tmp)
    {
     r8 ws;
     int ks;
     int k;
     r8 *t0, *t1, *t2, *t3;
     r8 tmp0[4*10];
     int iclocks;
     c8 sphase, zphase;

     t0 = tmp0;
     t1 = tmp0 + 10;
     t2 = tmp0 + 20;
     t3 = tmp0 + 30;

// Move first column of A into t0 shifted and padded appropriately

     t1[0] = A[4] - x;
     t1[1] = 0.;
     t1[2] = A[6];
     t1[3] = A[7];
     t1[4] = A[8];
     t1[5] = A[9];
     t1[6] = 0.;
     t1[7] = 0.;
     t1[8] = 0.;
     t1[9] = 0.;

     A += 10;

     ws = 1.;
#ifdef DO_CLOCKS
//     iclocks = pcc_delta();
#endif
     if(t1[0] < 0.)
       {
        zphase.r = -1.;
        zphase.i = 0.;
        ks = 0;
       }
     else
       {
        zphase.r = 1.;
        zphase.i = 0.;
        ks = 1;
       }
 

// k = 1 case, moved outside k loop

       {
        r8 ar, ai, ak;
        k = 1;

// Move second column of A into t1 shifted and padded appropriately

        t2[0] = A[2];
        t2[1] = A[3];
        t2[2] = A[4] - x;
        t2[3] = 0.;
        t2[4] = A[6];
        t2[5] = A[7];
        t2[6] = A[8];
        t2[7] = A[9];
        t2[8] = 0.;
        t2[9] = 0.;

        A += 10;

        ai = t1[0]*t1[0] + t1[1]*t1[1];
        ak = t2[0]*t2[0] + t2[1]*t2[1];

        if(ak > ai)
          { r8 *t = t1; t1 = t2; t2 = t; ws = -ws;} // swap pointers

        ar = t1[0];
        ai = t1[1];

        {
         r8 xr, xi, ak;
         ak = 1./(ar*ar + ai*ai);
         xr = ak*(t2[0]*ar + t2[1]*ai);
         xi = ak*(t2[1]*ar - t2[0]*ai);

         t3[0] = t2[2] - xr*t1[2+0] + xi*t1[2+1];
         t3[1] = t2[3] - xi*t1[2+0] - xr*t1[2+1];

         t3[2] = t2[4] - xr*t1[2+2] + xi*t1[2+3];
         t3[3] = t2[5] - xi*t1[2+2] - xr*t1[2+3];

         t3[4] = t2[6] - xr*t1[2+4] + xi*t1[2+5];
         t3[5] = t2[7] - xi*t1[2+4] - xr*t1[2+5];

         t3[6] = t2[8] - xr*t1[2+6] + xi*t1[2+7];
         t3[7] = t2[9] - xi*t1[2+6] - xr*t1[2+7];

         t3[8] = 0.;
         t3[9] = 0.;
        }

      { r8 *t = t3; t3 = t2; t2 = t;} // swap ptrs to t2, t3

// k == 1 case
      {
         c8 yphase;
            sphase.r = 1.;  // dont use saved phase yet.
            sphase.i = 0.;

         {
          r8 tr, ti, tt;
          r8 sr, si;
          tr = t1[0];
          ti = t1[1];
          tt = tr*tr + ti*ti;
          tt = sqrt(tt);
          if(tt > 0.)
            {
             tt = 1./tt;
             sr = sphase.r;
             si = sphase.i;
             sphase.r = tt*(sr*tr - si*ti);
             sphase.i = tt*(sr*ti + si*tr);
            }
         }

//         printf("k %d sphase %lf %lf\n",1,sphase.r, sphase.i);
         yphase = sphase;

        {
         r8 xr, xi;
         xr = yphase.r;
         xi = yphase.i;
           {
            r8 tr, ti, xt;
            tr = t2[0];
            ti = t2[1];
            xt = xr;
            xr = xt * tr - xi * ti;
            xi = xt * ti + xi * tr;
           }
         yphase.r = xr*ws;
         yphase.i = xi*ws;
        }

//       printf("k = 1 phase %lf %lf\n",yphase.r, yphase.i);

// count agreements in sign
         if(     (yphase.r >= 0. && zphase.r >= 0.)
              || (yphase.r  < 0. && zphase.r  < 0.))
           ks++;

// save computed phase for next comparison
         zphase = yphase;

         }

     // Cycle pointers
       {r8 *t = t0;
        t0 = t1; t1 = t2; t2 = t3; t3 = t;
       }

      }

     sphase.r = 1.;
     sphase.i = 0.;

     for(k=2;k<n;k++)
       {
        r8 ar, ai, ak;

    // Move next column of A into t2 shifted and padded appropriately

        t2[0] = A[0];
        t2[1] = A[1];
        t2[2] = A[2];
        t2[3] = A[3];
        t2[4] = A[4] - x;
        t2[5] = 0.;
        t2[6] = A[6];
        t2[7] = A[7];
        t2[8] = A[8];
        t2[9] = A[9];

        A += 10;

        ai = t0[0]*t0[0] + t0[1]*t0[1];
        ak = t2[0]*t2[0] + t2[1]*t2[1];

        if(ak > ai)
          { r8 *t = t0; t0 = t2; t2 = t; ws = -ws;} // swap pointers

        ar = t0[0];
        ai = t0[1];

        {
         r8 xr, xi, ak;
         ak = 1./(ar*ar + ai*ai);
         xr = ak*(t2[0]*ar + t2[1]*ai);
         xi = ak*(t2[1]*ar - t2[0]*ai);

//         printf("%d %d %lf %lf\n",k, 1, xr, xi);

         t3[0] = t2[2] - xr*t0[2+0] + xi*t0[2+1];
         t3[1] = t2[3] - xi*t0[2+0] - xr*t0[2+1];

         t3[2] = t2[4] - xr*t0[2+2] + xi*t0[2+3];
         t3[3] = t2[5] - xi*t0[2+2] - xr*t0[2+3];

         t3[4] = t2[6] - xr*t0[2+4] + xi*t0[2+5];
         t3[5] = t2[7] - xi*t0[2+4] - xr*t0[2+5];

         t3[6] = t2[8] - xr*t0[2+6] + xi*t0[2+7];
         t3[7] = t2[9] - xi*t0[2+6] - xr*t0[2+7];

         t3[8] = 0.;
         t3[9] = 0.;
        }

      { r8 *t = t3; t3 = t2; t2 = t;} // swap ptrs to t2, t3

        ai = t1[0]*t1[0] + t1[1]*t1[1];
        ak = t2[0]*t2[0] + t2[1]*t2[1];

        if(ak > ai)
          { r8 *t = t1; t1 = t2; t2 = t; ws = -ws;} // swap pointers

        ar = t1[0];
        ai = t1[1];

        {
         r8 xr, xi, ak;
         ak = 1./(ar*ar + ai*ai);
         xr = ak*(t2[0]*ar + t2[1]*ai);
         xi = ak*(t2[1]*ar - t2[0]*ai);

//         printf("%d %d %lf %lf\n",k, 2, xr, xi);

         t3[0] = t2[2] - xr*t1[2+0] + xi*t1[2+1];
         t3[1] = t2[3] - xi*t1[2+0] - xr*t1[2+1];

         t3[2] = t2[4] - xr*t1[2+2] + xi*t1[2+3];
         t3[3] = t2[5] - xi*t1[2+2] - xr*t1[2+3];

         t3[4] = t2[6] - xr*t1[2+4] + xi*t1[2+5];
         t3[5] = t2[7] - xi*t1[2+4] - xr*t1[2+5];

         t3[6] = t2[8] - xr*t1[2+6] + xi*t1[2+7];
         t3[7] = t2[9] - xi*t1[2+6] - xr*t1[2+7];

         t3[8] = 0.;
         t3[9] = 0.;
        }

      { r8 *t = t3; t3 = t2; t2 = t;} // swap ptrs to t2, t3

// k >= 2 case
      {
         c8 yphase;

         {
          r8 tr, ti, tt;
          r8 sr, si;
          tr = t0[0];
          ti = t0[1];
          tt = tr*tr + ti*ti;

          if(tt > 0.)
            {
/*
             tt = sqrt(tt);   // this is for scaling purposes only, to avoid overflow
             tt = 1./tt;
*/

// Well, I dont really like to use obscure extensions, but here goes. We
// should just be able to approximate the square root. ilogb is the binary
// exponent, scalbn scales 1 by minus half the exponent. Note that CPLANT
// (alpha?) does not have scalbn/ilogb.

#ifdef CPLANT
             tt = 1./sqrt(tt);
#else
             tt = scalbn(1., -(ilogb(tt)>>1));  // approximation to the square
#endif

             sr = sphase.r;
             si = sphase.i;
             sphase.r = tt*(sr*tr - si*ti);
             sphase.i = tt*(sr*ti + si*tr);
            }
         }

         yphase = sphase;
//         printf("k %d sphase %lf %lf\n",k,sphase.r, sphase.i);

        {
         r8 *tp = t1;
         r8 xr, xi;
         xr = yphase.r;
         xi = yphase.i;
           {
            r8 tr, ti, xt;
            tr = t1[0];
            ti = t1[1];
            xt = xr;
            xr = xt * tr - xi * ti;
            xi = xt * ti + xi * tr;
           }
           {
            r8 tr, ti, xt;
            tr = t2[0];
            ti = t2[1];
            xt = xr;
            xr = xt * tr - xi * ti;
            xi = xt * ti + xi * tr;
           }
         yphase.r = xr*ws;
         yphase.i = xi*ws;
        }

// count agreements in sign

//       printf("k %d phase %lf %lf ks %d ws %lf\n",k,yphase.r,yphase.i,ks,ws);

         if(     (yphase.r >= 0. && zphase.r >= 0.)
              || (yphase.r < 0. && zphase.r < 0.))
           ks++;

// save computed phase for next comparison
         zphase = yphase;

         }

     // Cycle pointers
       {r8 *t = t0;
        t0 = t1; t1 = t2; t2 = t3; t3 = t;
       }

      } // end k loop

#ifdef DO_CLOCKS
//     iclocks = pcc_delta();
//     printf("clocks %d\n",iclocks);
#endif
//     printf("oracle x %lf ks %d ws %lf\n",x,ks,ws);

     return ks;
    }
