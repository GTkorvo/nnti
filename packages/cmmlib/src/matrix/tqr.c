#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmatrix.h"

//void  cmm_c8_m_norms(c8 *a, int as, int an, int am, r8 *norms);

static void usage()
  {
   cmm_warn("Usage: tqr m [n]  where m >= n. \n");
   cmm_exit(0);
  }

static c8 z_mone = {-1., 0.};

void init_func(int i, int j, c8 *e, void *p)
  {
   c8 *tpar = (c8 *) p;
   r8 di, dj;

   di = DBLE(i);
   dj = DBLE(j);

   if(i == j)
     {
      e->r = 10. + (di*di + dj*dj);
      e->i = 0.;
     }
   else
     {
      e->r = 1./(di*di + dj*dj);
      e->i = e->r * (di - dj);
     }
  }

int main(int argc, char **argv)
  {
   int m,n;
   static c8 tpar[6] =
               {{.2, .43}, {.52, -.46},
               {1.389, .73},  {-1.01, .779},
               {1.987, .211}, {-.735, .015}
               }; 

   c8 *a;
   c8 *q; int qs;
   c8 *r; int rs;

   c8 *acopy;
   c8 *atmp;
   c8 *ytmp;
   c8 *itmp;
   r8 tn[4];
   r8 t0, mflops;
   r8 t1, mf1;
   int as;

   cmm_preallocate(80*1024);

   if(argc < 2)
     usage();

   m = cmm_str_cvt_int(argv[1]);
   if(argc == 3)
     n = cmm_str_cvt_int(argv[2]);
   else
     n = m;

   printf("tqr m %d n %d\n",m,n);

   as = m;
//   if ((as % 64) == 0) as += 5;
   a     =  cmm_c8_m_alloc(as, n);
   acopy =  cmm_c8_m_alloc(as, n);

   qs = m;
   q     =  cmm_c8_m_alloc(qs, n);

   rs = n;
   r     =  cmm_c8_m_alloc(rs, n);

    cmm_c8_m_apply(a, as, m, n, init_func, tpar);

    cmm_c8_m_copy(a, as, m, n, acopy, m);

   printf("Check QR factorization m %d n %d\n",m,n);
   t0 = cmm_clock();

   cmm_c8_qr(a, as, m, n,
             q, qs, m, n, 0,  // use default block size
             r, rs, n, n);

   t0 = cmm_clock() - t0;
   t0 = MAX(t0, 1.e-6);

// not sure about correct prefactor

   mflops = 5.3333e-6 * DBLE(m)* DBLE(n)* DBLE(n);  // MFlops

   printf("QR time %lf MFlops %lf  MFlop/sec %lf\n",t0,mflops,mflops/t0);

if(m < 8)
  {
   cmm_c8_m_print(q, qs, m, n, NULL);

   cmm_c8_m_print(r, rs, n, n, NULL);
  }


/*
   Test the result
*/
// for testing
   atmp  =  cmm_c8_m_alloc(as, n);
   ytmp  =  cmm_c8_m_alloc(as, n);
   itmp  =  cmm_c8_m_alloc(as, n);

/*
   printf("A'*A\n");
   cmm_c8_mm_hn(a, as, m, n, a, as, m, n, ytmp, as, n, n);
   cmm_c8_m_print(ytmp, as, n, n, NULL);

   printf("R'*R\n");
   cmm_c8_mm_hn(r, rs, n, n, r, rs, n, n, ytmp, as, n, n);
   cmm_c8_m_print(ytmp, as, n, n, NULL);
*/

   return 0;
  }

