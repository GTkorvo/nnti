#include "cmm.h"
#include "cmmblas.h"

/*

    Test matrix multiply C = C + A*B

       C is N by M
       A is N by K
       B is K by M

   We assume that the multiply works correctly, see
   the tcheck program. A, B, and C are initialized
   to random matrices. Then we run a performance
   test.

*/

static void tperf_nn(int n, int m, int k);
static void tperf_hn(int n, int m, int k);

static void usage()
  {
   cmm_warn("Usage: tperf n [m k]\n");
   cmm_exit(0);
  }

int main(int argc, char **argv)
  {
   int n, m, k;

   cmm_preallocate(80*1024);

   if(argc < 2)
     usage();

   n = cmm_str_cvt_r8(argv[1]);
   m = n;
   k = n;

   if(argc > 2)
     m = cmm_str_cvt_r8(argv[2]);

   if(argc > 3)
     k = cmm_str_cvt_r8(argv[3]);

   printf("Performance complex matrix multiply A*B  N %d M %d K %d\n", n, m, k);
   tperf_nn(n, m, k);

   printf("Performance complex matrix multiply A'*B N %d M %d K %d\n", n, m, k);
   tperf_hn(n, m, k);

   return 0;
  }

static void init_func(int i, int j, c8 *e, void *p)
  {
   c8 *tp = (c8 *) p;

   e->r = DBLE(i)*tp[0].r + DBLE(j)*tp[1].r;
   e->i = DBLE(i)*tp[0].i + DBLE(j)*tp[1].i;
  }

static c8 *init_test(int n, int m, c8 *tpar)
  {
   c8 *ap = cmm_c8_m_alloc(n, m);

   cmm_c8_m_apply(ap, n, n, m, init_func, (void *) tpar);
   return ap;
  }

static void tperf_nn(int n, int m, int k)
  {
   r8 mflops;
   int ntest;
   int itest;
   r8 t0;
   c8 p[2] = {{1.23, -.46},{.343, .233}};
   c8 *a, *b, *c;

   mflops = 8.*DBLE(n)*DBLE(m)*DBLE(k);

   a = init_test(n, k, p);
   b = init_test(k, m, p);
   c = init_test(n, m, p);

   ntest = INT(1.e9 / mflops);

   if(ntest < 1)
    ntest = 1;

//   ntest *= 10;

   printf("ntest %d\n",ntest);
   fflush(stdout);

   t0 = cmm_clock();
   for(itest = 0;itest < ntest; itest++)
     {
      cmm_c8_mm_a_nn(a, n, n, k, b, k, k, m, c, n, n, m);
     }
   t0 = cmm_clock() - t0;

   mflops = mflops * 1.e-6 *DBLE(ntest);

   printf("time %lf mflops %lf mflops/sec %lf\n",t0, mflops, mflops/t0);
   fflush(stdout);

  }

static void tperf_hn(int n, int m, int k)
  {
   r8 mflops;
   int ntest;
   int itest;
   r8 t0;
   c8 p[2] = {{1.23, -.46},{.343, .233}};
   c8 *a, *b, *c;

   mflops = 8.*DBLE(n)*DBLE(m)*DBLE(k);

   a = init_test(k, n, p);
   b = init_test(k, m, p);
   c = init_test(n, m, p);

   ntest = INT(1.e9 / mflops);

   if(ntest < 1)
    ntest = 1;

   printf("ntest %d\n",ntest);

   t0 = cmm_clock();
   for(itest = 0;itest < ntest; itest++)
     {
      cmm_c8_mm_a_hn(a, k, k, n, b, k, k, m, c, n, n, m);
     }
   t0 = cmm_clock() - t0;

   mflops = mflops * 1.e-6 *DBLE(ntest);

   printf("time %lf mflops %lf mflops/sec %lf\n",t0, mflops, mflops/t0);

  }
