#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmatrix.h"

//void  cmm_c8_m_norms(c8 *a, int as, int an, int am, r8 *norms);

static void usage()
  {
   cmm_warn("Usage: tchol n\n");
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
   int n;
   static c8 tpar[6] =
               {{.2, .43}, {.52, -.46},
               {1.389, .73},  {-1.01, .779},
               {1.987, .211}, {-.735, .015}
               }; 

   c8 *a;
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

   n = cmm_str_cvt_r8(argv[1]);
   as = n;
   if (as %64 == 0) as += 5;
   a     =  cmm_c8_m_alloc(as, n);
   acopy =  cmm_c8_m_alloc(as, n);
   atmp  =  cmm_c8_m_alloc(as, n);
   ytmp  =  cmm_c8_m_alloc(as, n);
   itmp  =  cmm_c8_m_alloc(as, n);

    cmm_c8_m_apply(a, as, n, n, init_func, tpar);

    cmm_c8_m_copy(a, as, n, n, acopy, n);

   printf("Check cholesky factorization n %d\n",n);
   t0 = cmm_clock();

//   cmm_c8_cholesky_unblocked(a, as, n, n);
   cmm_c8_cholesky(a, as, n, n);

   t0 = cmm_clock() - t0;
   t0 = MAX(t0, 1.e-6);

   mflops = 1.3333e-6 * DBLE(n)* DBLE(n)* DBLE(n);  // MFlops

   printf("cholesky time %lf MFlops %lf  MFlop/sec %lf\n",t0,mflops,mflops/t0);

    cmm_c8_m_clear_upper(a, as, n, n);
   if(n <= 8)
    {
     printf("C\n");
      cmm_c8_m_print(a, as, n, n, NULL);
    }

   t0 = cmm_clock();
    cmm_c8_mm_nh(a, as, n, n, a, as, n, n, atmp, as, n, n); // Atmp = C*C'
   t0 = cmm_clock() - t0;
   t0 = MAX(t0, 1.e-6);
   mflops = 8.e-6 * DBLE(n)* DBLE(n)* DBLE(n);  // MFlops
   printf("mm_nh time %lf MFlops %lf  MFlop/sec %lf\n",t0,mflops,mflops/t0);

   if(n <= 8)
    {
     printf("acopy\n");
      cmm_c8_m_print(acopy, as, n, n, NULL);
    }

   if(n <= 8)
    {
     printf("C * C'\n");
      cmm_c8_m_print(atmp, as, n, n, NULL);
    }

    cmm_c8_m_axpy(z_mone, acopy, as, n, n, atmp, as);      // Atmp = C*C' - A, should be zero 

   if(n <= 8)
    {
     printf("C * C'\n");
      cmm_c8_m_print(atmp, as, n, n, NULL);
    }

    cmm_c8_m_norms(atmp, as, n, n, tn);
   printf("|A - C*C'| %le\n",tn[3]);

// Test tsolve_la

    cmm_c8_m_identity(itmp, as, n);
   t0 = cmm_clock();
   c8_tsolve_la(a, as, n, n,
             itmp, as, n, n, 
             ytmp, as, n, n);
   t0 = cmm_clock() - t0;
   t0 = MAX(t0, 1.e-6);
   mflops = 4.e-6 * DBLE(n)* DBLE(n)* DBLE(n);  // MFlops
   printf("tsolve_la time %lf MFlops %lf  MFlop/sec %lf\n",t0,mflops,mflops/t0);

// Check tsolve_la

   t0 = cmm_clock();
    cmm_c8_mm_nn(a, as, n, n, ytmp, as, n, n, atmp, as, n, n); // Atmp = L*Y
   t0 = cmm_clock() - t0;
   t0 = MAX(t0, 1.e-6);
   mflops = 8.e-6 * DBLE(n)* DBLE(n)* DBLE(n);  // MFlops

    cmm_c8_m_axpy(z_mone, itmp, as, n, n, atmp, as);      // Atmp = Itmp - Y, should be zero 
   cmm_c8_m_norms(atmp, as, n, n, tn);
   printf("|I-L*Y| %le\n",tn[3]);

// Test tsolve_alt

    cmm_c8_m_identity(itmp, as, n);
   t0 = cmm_clock();
   c8_tsolve_alt(a, as, n, n,
              itmp, as, n, n, 
              ytmp, as, n, n);
   t0 = cmm_clock() - t0;
   t0 = MAX(t0, 1.e-6);
   mflops = 4.e-6 * DBLE(n)* DBLE(n)* DBLE(n);  // MFlops
   printf("tsolve_alt time %lf MFlops %lf  MFlop/sec %lf\n",t0,mflops,mflops/t0);

// Check tsolve_alt

   t0 = cmm_clock();
    cmm_c8_mm_nh(ytmp, as, n, n, a, as, n, n, atmp, as, n, n); // Atmp = Y*L'
   t0 = cmm_clock() - t0;
   t0 = MAX(t0, 1.e-6);
   mflops = 8.e-6 * DBLE(n)* DBLE(n)* DBLE(n);  // MFlops

    cmm_c8_m_axpy(z_mone, itmp, as, n, n, atmp, as);      // Atmp = Itmp - Y, should be zero 
    cmm_c8_m_norms(atmp, as, n, n, tn);
   printf("|I-Y*L'| %le\n",tn[3]);

// Test tsolve_lta

    cmm_c8_m_identity(itmp, as, n);
   t0 = cmm_clock();
   c8_tsolve_lta(a, as, n, n,
              itmp, as, n, n, 
              ytmp, as, n, n);
   t0 = cmm_clock() - t0;
   t0 = MAX(t0, 1.e-6);
   mflops = 4.e-6 * DBLE(n)* DBLE(n)* DBLE(n);  // MFlops
   printf("tsolve_lta time %lf MFlops %lf  MFlop/sec %lf\n",t0,mflops,mflops/t0);

// Check tsolve_lta

   t0 = cmm_clock();
    cmm_c8_mm_hn(a, as, n, n, ytmp, as, n, n, atmp, as, n, n); // Atmp = Y*L'
   t0 = cmm_clock() - t0;
   t0 = MAX(t0, 1.e-6);
   mflops = 8.e-6 * DBLE(n)* DBLE(n)* DBLE(n);  // MFlops

    cmm_c8_m_axpy(z_mone, itmp, as, n, n, atmp, as);      // Atmp = Itmp - Y, should be zero 
    cmm_c8_m_norms(atmp, as, n, n, tn);
   printf("|I-L'*Y| %le\n",tn[3]);

   return 0;
  }

