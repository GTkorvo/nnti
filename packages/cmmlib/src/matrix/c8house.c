#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmatrix.h"

#define BMAX 16

void c8_mvl(c8 *a, int as, int am, int an,
            c8 *x, int xs, int xm, int xn,
            c8 *y, int ys, int ym, int yn);

void   c8_updl(
           c8 *a, int as,  int am, int an,
           c8 *w, int ws,  int wm, int wn,
           c8 *v, int vs,  int vm, int vn);

void c8_mvl_2(c8 *a, int as, int am, int an,
            c8 *x, int xs, int xm, int xn,
            c8 *y, int ys, int ym, int yn);

void   c8_updl_2(
           c8 *a, int as,  int am, int an,
           c8 *w, int ws,  int wm, int wn,
           c8 *v, int vs,  int vm, int vn);

/*

    This routine builds a block of Householder vectors from a 
  block matrix A which is n by m. A is overwritten with the 
  Householder vectors, and the matrix B is filled with the 
  modified nonzero elements of A. B is thus upper triangular m 
  by m, and the Householder vectors are lower triangular
  n by m.

    This routine is not intended for the case where m is very
  large. In that case the Householder vectors can be built in
  blocks and then applied as a WY transform to the remainder
  of the matrix.

*/
void c8_build_house(c8 *a, int as, int m, int n, c8 *b, int bs)
  {
   int j, k;

   cmm_c8_m_clear(b, bs, n, n);

   for(j=0;j<n;j++)
     {
      c8 *hv = a + j*as;
      int mj = m-j; // length of this Householder vector
      c8 a0;

     // Build the Householder vector in place

     // Put away modified elements into upper half of B
      cmm_c8_copy(j, hv, 1, b + bs*j, 1);

     // clear initial elements
      c8_zero(j, hv);

      if(mj <= 0)
        break;

      {
       r8   v1r,v1i,v1n;
       c8 zn;
       r8 vm;
       r8 t;
       c8 h0;

       a0 = hv[j];
       h0 = a0;

       v1r = h0.r;
       v1i = h0.i;
       v1n = sqrt(v1r*v1r + v1i*v1i);
       zn = c8_hdot(mj-1, hv+j+1 ,hv+j+1);
       vm = zn.r;

       t = sqrt(v1n*v1n + vm);

       if(t == 0.)
         {
          h0.r = 0.;
          h0.i = 0.;
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

       // compute the modified first element of a
       a0.r = -t*v1r;
       a0.i = -t*v1i;

       // scale Householder vector
       t = v1n+t;
       h0.r = v1r*t;
       h0.i = v1i*t;
       hv[j] = h0;

       t = sqrt(t*t + vm);

       if(t > 0.)
         r8_scale(2*mj, 1./t, (r8 *) (hv+j));

      }

      C8_M_ELEMENT(b, bs, j, j)[0] = a0;  // put away modified element

     // Apply to remainder of the block

      for(k=j+1;k<n;k++)
        {
         c8 h0;
         c8 *ap = a + k*as;

         h0 = c8_hdot(mj, ap+j, hv+j);
         h0.r *= -2.;
         h0.i *= -2.;
         c8_axpy(mj, h0, hv+j, ap+j);
        }
     }
  }

/*

    Build the WY transform from a block of Householder vectors.

*/
void c8_build_wy(c8 *y, int ys, int m, int n, c8 *w, int ws)
  {
   int i, j, k;

   c8 ztmp[BMAX*BMAX];
   c8 *zp;

   if(n<=BMAX)
     zp = ztmp;
   else
     zp = cmm_alloc(sizeof(c8)*m*m);

   for(j=0;j<n;j++)
    for(k=0;k<j;k++)
      {
       c8 z0;
       z0 = c8_hdot(m, y+k*ys, y+j*ys);
       z0.r *= -2.;
       z0.i *=  2.;
       zp[k+j*n] = z0;
      }

   for(j=0;j<n;j++)
    {
     c8 *wp = w + j*ws;
     c8 *yp = y + j*ys;
     for(i=0;i<m;i++)
       {
        c8 t = yp[i];
        t.r *= -2.;
        t.i *= -2.;
        wp[i] = t;
       }
     for(k=0;k<j;k++)
       c8_axpy(m, zp[k+j*n], w+k*ws, wp);
    }

   if(n > BMAX)
    cmm_free(zp);
  }

/*

    Update A' = (I+WY')'*A*(I+WY')

*/
#ifdef DO_CLOCKS
int pcc_delta();
#endif
void c8_update_wy(
             c8 *y, int ys, int ym, int yn,
             c8 *w, int ws, int wm, int wn,
             c8 *a, int as, int am, int an,
             c8 *t, int nt
#ifdef DO_CLOCKS
, int *clks
#endif
)
  {
   c8 *x;   // X = A*Y
   c8 *u;
   c8 *v;
   c8 z[BMAX*BMAX];
   int i, j, k;
   int m = ym;
   int b = yn;

// Allocate temp storage (3*n*b complex)
   x = t;
   u = x + m*b;
   v = u + m*b;

// Compute X = A*W

//   printf("called update_wy %d\n",b);

   if(b == 2)
     c8_mvl_2(a, as, am, an,
          w, ws, m, b,
          x, m, m, b);
   else
     c8_mvl(a, as, am, an,
          w, ws, m, b,
          x, m, m, b);

#ifdef DO_CLOCKS
   clks[0] = pcc_delta();
#endif

// Compute Z = .5*Y'*X
   for(j=0;j<b;j++)
   for(i=0;i<b;i++)
     {
       c8 sum;

       sum.r = 0.;
       sum.i = 0.;

       for(k=0;k<m;k++)
         {
          c8 we = w[k + ws*i];
          c8 xe = x[k + m*j];
          sum.r += we.r * xe.r + we.i * xe.i;
          sum.i += we.r * xe.i - we.i * xe.r;
         }

       sum.r *= .5;
       sum.i *= .5;
       z[i+b*j] = sum;
     }

// Compute  V = X + Y*Z

   for(i=0;i<m;i++)
    for(j=0;j<b;j++)
      {
       c8 sum = x[i+m*j];

       for(k=0;k<b;k++)
         {
          c8 ye = y[i + ys*k];
          c8 ze = z[k + b*j];
          sum.r += ye.r * ze.r - ye.i * ze.i;
          sum.i += ye.r * ze.i + ye.i * ze.r;
         }

       v[i+m*j] = sum;
      }

// Update A = A + V*Y' + Y*V'
   if(b == 2)
     c8_updl_2(a, as, m, m,
           y, ys, m, b,
           v, m,  m, b);
   else
     c8_updl(a, as, m, m,
           y, ys, m, b,
           v, m,  m, b);

#ifdef DO_CLOCKS
   clks[1] = pcc_delta();
#endif

  }

/*
     Compute Y = A*X using only lower half of A (assumed Hermitian)
     This is a general version.
*/
void c8_mvl(c8 *a, int as, int am, int an,
            c8 *x, int xs, int xm, int xn,
            c8 *y, int ys, int ym, int yn)
  {
   int i, j, k;

   for(k=0;k<xn;k++)
     c8_zero(ym, y+k*ys);

// This is 16 flops per inner step
// times b * n^2 / 2. For b=2 we have
// 16 n^2 flops

   for(j=0;j<am;j++)
     {
      for(k=0;k<xn;k++)
        {
         c8 t;
         c8 *ap0 = a + as*j;
         c8 *xp0 = x + xs*k;
         c8 *yp0 = y + ys*k;
         c8 ae   = ap0[j];
         c8 xe   = xp0[j];

         t.r = yp0[j].r + ae.r * xe.r - ae.i * xe.i;
         t.i = yp0[j].i + ae.r * xe.i + ae.i * xe.r;

         for(i=j+1;i<am;i++)
           {
            ae = ap0[i];

            yp0[i].r += ae.r * xe.r     - ae.i * xe.i;
            yp0[i].i += ae.r * xe.i     + ae.i * xe.r;
                 t.r += ae.r * xp0[i].r + ae.i * xp0[i].i;
                 t.i += ae.r * xp0[i].i - ae.i * xp0[i].r;

           }

         yp0[j] = t;
        }
     }
  }

// mvl unrolled b loop
void c8_mvl_2(c8 *a, int as, int am, int an,
            c8 *x, int xs, int xm, int xn,
            c8 *y, int ys, int ym, int yn)
  {
   int i, j, k;
   c8 *y0 = y;
   c8 *y1 = y + ys;
   c8 *x0 = x;
   c8 *x1 = x + xs;

   c8_zero(ym, y0);
   c8_zero(ym, y1);

   for(j=0;j<am;j++)
     {
      c8 *ap = a + as*j;
      c8 ae = ap[j];
      c8 xe0 = x0[j];
      c8 xe1 = x1[j];
      c8 t0, t1;

      t0.r = y0[j].r + ae.r * xe0.r - ae.i * xe0.i;
      t0.i = y0[j].i + ae.r * xe0.i + ae.i * xe0.r;
      t1.r = y1[j].r + ae.r * xe1.r - ae.i * xe1.i;
      t1.i = y1[j].i + ae.r * xe1.i + ae.i * xe1.r;

      for(i=j+1;i<am;i++)
        {
         c8 ae = ap[i];

         y0[i].r += ae.r * xe0.r     - ae.i * xe0.i;
         y0[i].i += ae.r * xe0.i     + ae.i * xe0.r;
         y1[i].r += ae.r * xe1.r     - ae.i * xe1.i;
         y1[i].i += ae.r * xe1.i     + ae.i * xe1.r;

         t0.r += ae.r * x0[i].r + ae.i * x0[i].i;
         t0.i += ae.r * x0[i].i - ae.i * x0[i].r;
         t1.r += ae.r * x1[i].r + ae.i * x1[i].i;
         t1.i += ae.r * x1[i].i - ae.i * x1[i].r;

        }
      y0[j] = t0;
      y1[j] = t1;
     }
  }

/*
    Compute the update

     A = A + W*V' + V*W'

     where W and V are thin.

     A(i,j) = A(i,j) + sum(k) (W(i,k)*V'(k,j) + V(i,k)*W'(k,j))
            = A(i,j) + sum(k) (W(i,k)*Vc(j,k) + V(i,k)*Wc(j,k))
*/
void   c8_updl(
           c8 *a, int as,  int am, int an,
           c8 *w, int ws,  int wm, int wn,
           c8 *v, int vs,  int vm, int vn)
  {
   int i, j, k;
   int m, b;

   b = wn;
   m = am;

   for(j=0;j<m;j++)
    for(i=j;i<m;i++)
      {
       c8 sum = a[i + as*j];
       c8 *wp0  = w + i;
       c8 *whp0 = w + j;
       c8 *vp0  = v + i;
       c8 *vhp0 = v + j;

       for(k=0;k<b;k++)
         {
          c8 we, whe;
          c8 ve, vhe;

          we  = *wp0;  wp0  += ws;
          whe = *whp0; whp0 += ws;
          ve  = *vp0;  vp0  += vs;
          vhe = *vhp0; vhp0 += vs;
          sum.r += (we.r * vhe.r + we.i * vhe.i + ve.r * whe.r + ve.i * whe.i);
          sum.i += (we.i * vhe.r - we.r * vhe.i + ve.i * whe.r - ve.r * whe.i);
         }
       a[i+as*j] = sum;
      }
  }

// sane but unrolled wn=2
void   c8_updl_2(
           c8 *a, int as,  int am, int an,
           c8 *w, int ws,  int wm, int wn,
           c8 *v, int vs,  int vm, int vn)
  {
   int i, j, k;
   int m, b;
   c8 *ap = a;

   b = wn;
   m = am;

   for(j=0;j<m;j++)
    {
    for(i=j;i<m;i++)
      {
       c8 sum = ap[i];
       c8 *wp0  = w + i;
       c8 *whp0 = w + j;
       c8 *vp0  = v + i;
       c8 *vhp0 = v + j;
       c8 we, whe;
       c8 ve, vhe;

       we  = *wp0;  wp0  += ws;
       whe = *whp0; whp0 += ws;
       ve  = *vp0;  vp0  += vs;
       vhe = *vhp0; vhp0 += vs;
       sum.r += (we.r * vhe.r + we.i * vhe.i + ve.r * whe.r + ve.i * whe.i);
       sum.i += (we.i * vhe.r - we.r * vhe.i + ve.i * whe.r - ve.r * whe.i);
       we  = *wp0;  
       whe = *whp0; 
       ve  = *vp0;  
       vhe = *vhp0; 
       sum.r += (we.r * vhe.r + we.i * vhe.i + ve.r * whe.r + ve.i * whe.i);
       sum.i += (we.i * vhe.r - we.r * vhe.i + ve.i * whe.r - ve.r * whe.i);

       ap[i] = sum;
      }
    ap += as;
    }
  }

/*

   Apply a Householder matrix to back-transform the matrix V.

*/
void c8_house_back_unblocked(c8 *a, int as, int am, int an, int bw, c8 *v, int vs, int vm, int vn)
  {
   int i,j;

   for(j=am-bw-1;j>=0;j--) // Transforms must be applied in reverse order
     {
      int icol = j+bw;    // starting location for Householder vector
      int mj = am-icol;  // length of the Householder vector
      c8 *hv = a + as*j;

      if(mj < 1)
        continue;

      for(i=0;i<vm;i++)
        {
         c8 *vp = v + i*vs;
         c8 t;

         t = c8_hdot(mj, vp+icol, hv+icol);
         t.r *= -2.;
         t.i *= -2.;
         c8_axpy(mj, t, hv+icol, vp+icol);
        }
     }
  }

/*

   This routine works quite well.

   Notes:

     The end case loop is not quite right, and the routine
   should work for the case where a is rectangular and bw=0,
   which could come from a QR decomposition.

*/
void c8_house_back(c8 *a, int as, int am, int an, int bw, c8 *v, int vs, int vm, int vn)
  {
   int i,j,jj;
   int b, mb;

   b = 16;

   mb = (am-bw)/b;

// take care of end cases
   for(j=am-bw;j>= mb*b;j--) // Transforms must be applied in reverse order
     {
      int icol = j+bw;    // starting location for Householder vector
      int nj = am-icol;  // length of the Householder vector
      c8 *hv = a + as*j;

      if(nj < 1)
        continue;

      for(i=0;i<vm;i++)
        {
         c8 *vp = v + i*vs;
         c8 t;

         t = c8_hdot(nj, vp+icol, hv+icol);
         t.r *= -2.;
         t.i *= -2.;
         c8_axpy(nj, t, hv+icol, vp+icol);
        }
     }

/*
   for(jj=mb;jj>0;jj--)
     {
      int j1;

      for(j1=0;j1<b;j1++)
        {
         int j = jj*b-j1-1;
         int icol = j+bw;   // starting col index for Householder vector
         int nj = am-icol;  // length of the Householder vector
         c8 *hv = a + as*j; // pointer to column with Householder vector

         if(nj < 1)
           continue;

         for(i=0;i<vm;i++)
           {
            c8 *vp = v + i*vs;
            c8 t;

            t = c8_hdot(nj, vp+icol, hv+icol);
            t.r *= -2.;
            t.i *= -2.;
            c8_axpy(nj, t, hv+icol, vp+icol);
           }
        }
     }
*/

// remaining block transforms

   {
    c8 *w;
    c8 *t;
    c8 *tmp;

// temp storage
// w is am*b
// t is b*vm

    tmp = cmm_alloc(sizeof(c8)*(b*am + b*vn));
    w = tmp;
    t = w + b*am; 

   for(jj=mb;jj>0;jj--)
     {
      int j1;
      int icol = ((jj-1)*b + bw);
      int nj = am - icol;
      int mj = b;
      c8 *y  = a + icol + as*(jj-1)*b;
      c8 *v1 = v + icol;

    // clear upper part of Y just in case
      {
       c8 *yp = y;
       for(j1=0;j1<b;j1++)
        {
         for(i=0;i<j1;i++)
           {
            yp[i].r = 0.;
            yp[i].i = 0.;
           }
         yp += as;
        }
      }

    // Build the WY transform
      c8_build_wy(y, as, nj, b, w, nj);
 
    // compute T = Y'*V
      cmm_c8_mm_hn(y,  as, nj, b,
               v1, vs, nj, vn,
               t,  b,   b, vn);

    // compute V = V + W*T
      cmm_c8_mm_a_nn(w, nj, nj, b,
                 t, b,   b, vn,
                v1, vs, nj, vn);

     } // end loop over block transforms

    cmm_free(tmp);
   }
  }
