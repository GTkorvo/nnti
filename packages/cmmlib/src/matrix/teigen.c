#include "cmm.h"
#include "cmmblas.h"
#include "cmmatrix.h"

#include "cmmrand.h"

static void usage()
  {
   cmm_warn("Usage: teigen n\n");
   cmm_exit(0);
  }

static c8 z_mone = {-1., 0.};

static void init_func_1(int i, int j, c8 *e, void *p)
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
   else  // if(i < 11 && j < 11)
     {
      e->r = 1./(di*di + dj*dj);
      e->i = e->r * (di - dj);
     }
/*
   else
     {
      e->r = 0.;
      e->i = 0.;
     }
*/
  }

/* This is Wilkinson's W21 matrix */
static void init_func_2(int i, int j, c8 *e, void *p)
  {
   c8 *tpar = (c8 *) p;
   r8 di, dj;
   r8 small = 0.;
   
   di = DBLE(i);
   dj = DBLE(j);

   if(i == j)
     {
      if(i < 10)
        e->r = DBLE(10-i);
      else
        e->r = DBLE(i-10);

      e->i = 0.;
     }
   else if(i == (j-1))
     {
      e->r = 1.;
      e->i = 0.;
     }
   else if(i == (j+1))
     {
      e->r = 1.;
      e->i = 0.;
     }
   else
     {
      e->r = small;
      e->i = 0.;
     }
  }

static void init_func_3(int i, int j, c8 *e, void *p)
  {
   c8 *tpar = (c8 *) p;
   r8 di, dj;

   di = DBLE(i);
   dj = DBLE(j);

   if(i == j)
     {
      e->r = -2.;
      e->i = 0.;
     }
   else if(i == j+1)
     {
      e->r = 1.;
      e->i = 0.;
     }
   else if(i == j-1)
     {
      e->r = 1.;
      e->i = 0.;
     }
   else
     {
      e->r = 0.;
      e->i = 0.;
     }
  }


static void scale_eigenvalues(int i, int j, c8 *e, void *p)
  {
   r8 *evals = (r8 *) p;

   e->r *= evals[j];
   e->i *= evals[j];

  }

int main(int argc, char **argv)
  {
   int n, b;
   static c8 tpar[6] =
               {{.2, .43}, {.52, -.46},
               {1.389, .73},  {-1.01, .779},
               {1.987, .211}, {-.735, .015}
               }; 

   c8 *a;
   c8 *acopy;
   c8 *v;
   c8 *itmp;
   c8 *ytmp;
   c8 *vtmp;

   r8 *evals;
   r8 tn[4];
   r8 t0, mflops;
   r8 t1, mf1;

   cmm_preallocate(80*1024);
   cmm_rand_init_seed(3498739, 838973);

   if(argc < 2)
     usage();

   n = cmm_str_cvt_r8(argv[1]);
   b = 2;

//   n = 21;

   a     = cmm_c8_m_alloc(n, n);
   acopy = cmm_c8_m_alloc(n, n);
   v     = cmm_c8_m_alloc(n, n);
   itmp  = cmm_c8_m_alloc(n, n);
   ytmp  = cmm_c8_m_alloc(n, n);
   vtmp  = cmm_c8_m_alloc(n, n);
   evals = cmm_alloc(sizeof(r8) * n);

   cmm_c8_m_apply(a, n, n, n, init_func_3, tpar);
   if(n < 16)
      cmm_c8_m_print(a, n, n, n, NULL);

   cmm_c8_m_copy(a, n, n, n, acopy, n);

   cmm_c8_m_norms(a, n, n, n, tn);

   printf("Initial tr(A) %23.18lf %23.18lf, e2(A)  %23.18lf\n",tn[0],tn[1],tn[2]);

   t0 = cmm_clock();
   c8_eigen(a, n, n, v, n, evals);
   t0 = cmm_clock() - t0;

   mflops = (5.333 + 8.)* 1.e-6 * DBLE(n) * DBLE(n) * DBLE(n);
   printf("Eigensolve time %lf MFlops %lf MFlops/sec %lf\n",t0,mflops,mflops/t0);

   if(n < 16)
     {
      printf("V\n");
      cmm_c8_m_print(v, n, n, n, NULL);
     }

// Now test the results
   {
    int i;
    r8 esum, e2;
    esum = 0.;
    e2 = 0.;
    for(i=0;i<n;i++)
      {
       if(n < 40)
         printf("ev[%d] %18.14lf\n",i,evals[i]);
       esum += evals[i];
       e2 += evals[i]*evals[i];
      }
    printf("Eigenvals tr(A) %23.18lf, e2(A) %23.18lf\n",esum, e2);
   }

// Test orthogonality, |V'*V - I| should be small compared to 1.
   cmm_c8_m_identity(itmp, n, n);
   t0 = cmm_clock();
   cmm_c8_mm_hn(v, n, n, n, v, n, n, n, ytmp, n, n, n); // Ytmp = V' * V
   t0 = cmm_clock() - t0;
   mflops = 8.e-6 *  DBLE(n) * DBLE(n) * DBLE(n);
   printf("mult hn time %lf mflops %lf MFlop/sec %lf\n",t0,mflops,mflops/t0);

   cmm_c8_m_axpy(z_mone, itmp, n, n, n, ytmp, n);      // Ytmp = Ytmp - I, should be zero 
   cmm_c8_m_norms(ytmp, n, n, n, tn);
   printf("|I-V'*V| %le\n",tn[3]);

// Test eigenpairs |A*V - V*E| should be small compared to |A|.
   t0 = cmm_clock();
   cmm_c8_mm_nn(acopy, n, n, n, v, n, n, n, ytmp, n, n, n);     // Ytmp = A * V
   t0 = cmm_clock() - t0;
   mflops = 8.e-6 *  DBLE(n) * DBLE(n) * DBLE(n);

   if( n < 32)
     {
      printf("A copy:\n");
      cmm_c8_m_print(acopy, n, n, n, NULL);     
      printf("A*V:\n");
      cmm_c8_m_print(ytmp, n, n, n, NULL);
     }

   printf("mult nn time %lf mflops %lf MFlop/sec %lf\n",t0,mflops,mflops/t0);

   cmm_c8_m_copy(v, n, n, n, vtmp, n);                      // Vtmp = V
   cmm_c8_m_apply(vtmp, n, n, n, scale_eigenvalues, evals); // Vtmp = V*E

   if(n < 32)
     cmm_c8_m_print(vtmp, n, n, n, NULL);
   cmm_c8_m_axpy(z_mone, vtmp, n, n, n, ytmp, n);           // Ytmp = A*V-V*E
   if(n < 32)
     cmm_c8_m_print(ytmp, n, n, n, NULL);
   cmm_c8_m_norms(ytmp, n, n, n, tn);
   printf("|A*V - V*E| %le\n",tn[3]);

   return 0;
  }
