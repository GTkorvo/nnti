#include "cmm.h"
#include "cmmblas.h"

/*

    Test matrix multiply C = C + A*B

       C is N by M
       A is N by K
       B is K by M

   Initially

     A(i,j) = a1*i+a2*j;
     B(i,j) = b1*i+b2*j;
     C(i,j) = c1*i+c2*j;

   so after the operation

     C(i,j) = c1*i + c2*j + a1*i*(b1*S(1,K) + b2*K*j) + a2*(b1*S(2,K) + b2*S(1,K)*j)

            = a2*b1*S(2,K) + i*(c1 + a1*b1*S(1,K)) + j*(c2 + a2*b2*S(1,K))
                  + i*j*a1*b2*K

   where

     S(1,K) = sum(k<K) k = K*(K-1)/2
     S(2,K) = sum(k<K) k^2 = K*(K-1)*(2*K-1)/6

   We generate some random parameters for a1, etc

   Have similar results for A'*B

*/

static void testmm_nn(int n, int m, int k, c8 *tpar);
static void testmm_hn(int n, int m, int k, c8 *tpar);
static void testmm_nh(int n, int m, int k, c8 *tpar);

static void usage()
  {
   cmm_warn("Usage: tcheck n [m k]\n");
   cmm_exit(0);
  }

int main(int argc, char **argv)
  {
   int n, m, k;
   static c8 tpar[6] =
               {{.2, .43}, {.52, -.46},
               {1.389, .73},  {-1.01, .779},
               {1.987, .211}, {-.735, .015}
               }; 

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

   printf("Check complex matrix multiply A*B N %d M %d K %d\n", n, m, k);

   testmm_nn(n, m, k, tpar);

   printf("Check complex matrix multiply A'*B N %d M %d K %d\n", n, m, k);

   testmm_hn(n, m, k, tpar);

   printf("Check complex matrix multiply A*B' N %d M %d K %d\n", n, m, k);

   testmm_nh(n, m, k, tpar);

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

static int test_status;
static int K;

static void test_func_nn(int i, int j, c8 *e, void *p)
  {
   c8 *tp = (c8 *) p;
   c8 vcheck;
   c8 a1 = tp[0];
   c8 a2 = tp[1];
   c8 b1 = tp[2];
   c8 b2 = tp[3];
   c8 c1 = tp[4];
   c8 c2 = tp[5];
   r8 dk = DBLE(K);
   r8 sk = dk*(dk-1.)/2.;
   r8 sk2 = dk*(dk-1.)*(2.*dk-1.)/6.;
   r8 er, ei, err;

/*
            = a2*b1*S(2,K) + i*(c1 + a1*b1*S(1,K)) + j*(c2 + a2*b2*S(1,K))
                  + i*j*a1*b2*K
*/

   vcheck.r = (a2.r*b1.r - a2.i*b1.i)*sk2;
   vcheck.i = (a2.r*b1.i + a2.i*b1.r)*sk2;

   vcheck.r += DBLE(i)*(c1.r + (a1.r*b1.r - a1.i*b1.i)*sk);
   vcheck.i += DBLE(i)*(c1.i + (a1.r*b1.i + a1.i*b1.r)*sk);

   vcheck.r += DBLE(j)*(c2.r + (a2.r*b2.r - a2.i*b2.i)*sk);
   vcheck.i += DBLE(j)*(c2.i + (a2.r*b2.i + a2.i*b2.r)*sk);

   vcheck.r += DBLE(i)*DBLE(j)*dk*(a1.r*b2.r - a1.i*b2.i);
   vcheck.i += DBLE(i)*DBLE(j)*dk*(a1.r*b2.i + a1.i*b2.r);

   er = vcheck.r - e->r;
   ei = vcheck.i - e->i;

   err = er*er + ei*ei;

   if(err > 1.e-6)
     {
      test_status = 1;
      printf("error for %d %d e %le %le vcheck %le %le\n",i,j,e->r,e->i,vcheck.r,vcheck.i);
     }

  }

static void check_test_nn(c8 *c, int cs, int cn, int cm, c8 *tpar)
  {
   test_status = 0;

   cmm_c8_m_apply(c, cs, cn, cm, test_func_nn, (void *) tpar);

   if(test_status)
     printf("check_test failed\n");
   else
     printf("check_test successful\n");
  }

static void testmm_nn(int n, int m, int k, c8 *tpar)
  {
   c8 *a, *b, *c;

   a = init_test(n, k, tpar);
   b = init_test(k, m, tpar+2);
   c = init_test(n, m, tpar+4);

   if(n < 10)
     {
      printf("am\n");
      cmm_c8_m_print(a, n, n, k, NULL);
      printf("bm\n");
      cmm_c8_m_print(b, k, k, m, NULL);
      printf("cm\n");
      cmm_c8_m_print(c, n, n, m, NULL);
     }

// C = C + A*B

   cmm_c8_mm_a_nn(a, n, n, k, b, k, k, m, c, n, n, m);

   if(n < 10)
     {
      printf("cm after\n");
      cmm_c8_m_print(c, n, n, m, NULL);
     }

   K = k;
   check_test_nn(c, n, n, m, tpar);

  }

static void test_func_hn(int i, int j, c8 *e, void *p)
  {
   c8 *tp = (c8 *) p;
   c8 vcheck;
   c8 a1 = tp[0];
   c8 a2 = tp[1];
   c8 b1 = tp[2];
   c8 b2 = tp[3];
   c8 c1 = tp[4];
   c8 c2 = tp[5];
   r8 dk = DBLE(K);
   r8 sk = dk*(dk-1.)/2.;
   r8 sk2 = dk*(dk-1.)*(2.*dk-1.)/6.;
   r8 er, ei, err;

   vcheck.r = (a1.r*b1.r + a1.i*b1.i)*sk2;
   vcheck.i = (a1.r*b1.i - a1.i*b1.r)*sk2;

   vcheck.r += DBLE(i)*(c1.r + (a2.r*b1.r + a2.i*b1.i)*sk);
   vcheck.i += DBLE(i)*(c1.i + (a2.r*b1.i - a2.i*b1.r)*sk);

   vcheck.r += DBLE(j)*(c2.r + (a1.r*b2.r + a1.i*b2.i)*sk);
   vcheck.i += DBLE(j)*(c2.i + (a1.r*b2.i - a1.i*b2.r)*sk);

   vcheck.r += DBLE(i)*DBLE(j)*dk*(a2.r*b2.r + a2.i*b2.i);
   vcheck.i += DBLE(i)*DBLE(j)*dk*(a2.r*b2.i - a2.i*b2.r);

   er = vcheck.r - e->r;
   ei = vcheck.i - e->i;

   err = er*er + ei*ei;

   if(err > 1.e-6)
     {
      test_status = 1;
      printf("error for %d %d e %le %le vcheck %le %le\n",i,j,e->r,e->i,vcheck.r,vcheck.i);
     }

  }

static void check_test_hn(c8 *c, int cs, int cn, int cm, c8 *tpar)
  {
   test_status = 0;

   cmm_c8_m_apply(c, cs, cn, cm, test_func_hn, (void *) tpar);

   if(test_status)
     printf("check_test failed\n");
   else
     printf("check_test successful\n");
  }

static void testmm_hn(int n, int m, int k, c8 *tpar)
  {
   c8 *a, *b, *c;

   a = init_test(k, n, tpar);
   b = init_test(k, m, tpar+2);
   c = init_test(n, m, tpar+4);

   if(n < 10)
     {
      printf("am\n");
      cmm_c8_m_print(a, k, k, n, NULL);
      printf("bm\n");
      cmm_c8_m_print(b, k, k, m, NULL);
      printf("cm\n");
      cmm_c8_m_print(c, n, n, m, NULL);
     }

// C = C + A*B

   cmm_c8_mm_a_hn(a, k, k, n, b, k, k, m, c, n, n, m);

   if(n < 10)
     {
      printf("cm after\n");
      cmm_c8_m_print(c, n, n, m, NULL);
     }

   K = k;
   check_test_hn(c, n, n, m, tpar);

  }

static void test_func_nh(int i, int j, c8 *e, void *p)
  {
   c8 *tp = (c8 *) p;
   c8 vcheck;
   c8 a1 = tp[0];
   c8 a2 = tp[1];
   c8 b1 = tp[2];
   c8 b2 = tp[3];
   c8 c1 = tp[4];
   c8 c2 = tp[5];
   r8 dk = DBLE(K);
   r8 sk = dk*(dk-1.)/2.;
   r8 sk2 = dk*(dk-1.)*(2.*dk-1.)/6.;
   r8 er, ei, err;

   vcheck.r = (a2.r*b2.r + a2.i*b2.i)*sk2;
   vcheck.i = (a2.i*b2.r - a2.r*b2.i)*sk2;

   vcheck.r += DBLE(i)*(c1.r + (a1.r*b2.r + a1.i*b2.i)*sk);
   vcheck.i += DBLE(i)*(c1.i + (a1.i*b2.r - a1.r*b2.i)*sk);

   vcheck.r += DBLE(j)*(c2.r + (a2.r*b1.r + a2.i*b1.i)*sk);
   vcheck.i += DBLE(j)*(c2.i + (a2.i*b1.r - a2.r*b1.i)*sk);

   vcheck.r += DBLE(i)*DBLE(j)*dk*(a1.r*b1.r + a1.i*b1.i);
   vcheck.i += DBLE(i)*DBLE(j)*dk*(a1.i*b1.r - a1.r*b1.i);

   er = vcheck.r - e->r;
   ei = vcheck.i - e->i;

   err = er*er + ei*ei;

   if(err > 1.e-6)
     {
      test_status = 1;
      printf("error for %d %d e %le %le vcheck %le %le\n",i,j,e->r,e->i,vcheck.r,vcheck.i);
     }

  }

static void check_test_nh(c8 *c, int cs, int cn, int cm, c8 *tpar)
  {
   test_status = 0;

   cmm_c8_m_apply(c, cs, cn, cm, test_func_nh, (void *) tpar);

   if(test_status)
     printf("check_test failed\n");
   else
     printf("check_test successful\n");
  }


static void testmm_nh(int n, int m, int k, c8 *tpar)
  {
   c8 *a, *b, *c;

   a = init_test(n, k, tpar);
   b = init_test(m, k, tpar+2);
   c = init_test(n, m, tpar+4);

   if(n < 10)
     {
      printf("am\n");
      cmm_c8_m_print(a, n, n, k, NULL);
      printf("bm\n");
      cmm_c8_m_print(b, m, m, k, NULL);
      printf("cm\n");
      cmm_c8_m_print(c, n, n, m, NULL);
     }

// C = C + A*B

   cmm_c8_mm_a_nh(a, n, n, k, b, m, m, k, c, n, n, m);

   if(n < 10)
     {
      printf("cm after\n");
      cmm_c8_m_print(c, n, n, m, NULL);
     }

   K = k;
   check_test_nh(c, n, n, m, tpar);

  }
