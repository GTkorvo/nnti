#include "cmm.h"
#include "cmmblas.h"
#include "cmmatrix.h"
#include "cmmband.h"

/*

   Compute eigenvalues and eigenvectors of Hermitian
   matrix

   Ordinary form:

     A*V = V*E

   The matrix A is destroyed.

*/

void c8_eigen(c8 *a, int as, int an, c8 *v, int vs, r8 *eigenvals)
  {
   c8 *aband;
   int bw;
   r8 t0;
   r8 mflops;

// Check memory
   printf("Memory summary before eigensolve\n");
   cmm_memstat(NULL);

// Reduce to condensed form
   bw = 2;
   aband = cmm_alloc(sizeof(c8)*an*(2*bw+1));  // allocate temp
   t0 = cmm_clock();
   c8_condense(a, as, an, bw, aband);
   t0 = cmm_clock() - t0;
   mflops = 5.3333e-6 * DBLE(an)* DBLE(an)* DBLE(an);
   printf("condense time %lf mflops %lf MFlop/sec %lf\n",t0,mflops,mflops/t0);

//   c8band_print(an, bw, bw, aband);

// Compute eigenpairs of condensed form
   c8band_eigen(an, bw, aband, 0, an, eigenvals, v, vs);
   printf("c8band time %lf mflops %lf mflop/sec %lf\n",
        c8band_eigen_time(), c8band_eigen_mflops(),
        c8band_eigen_mflops()/(1.e-12 + c8band_eigen_time()));

   cmm_free(aband);  // free temp storage

// Transform eigenvectors back to dense system
   if(an < 10)
     {
      printf("V before back transform\n");
      cmm_c8_m_print(v, vs, an, an, NULL);
     }

   t0 = cmm_clock();
   c8_house_back(a, as, an, an, bw, v, vs, an, an);
   t0 = cmm_clock() - t0;
   mflops = 8.e-6 * DBLE(an)* DBLE(an)* DBLE(an);
   printf("back tr time %lf mflops %lf MFlop/sec %lf\n",t0,mflops,mflops/t0);

// Check memory again
   printf("Memory summary after eigensolve\n");
   cmm_memstat(NULL);

  }

/*

   Compute eigenvalues and eigenvectors of Hermitian
   system

   Generalized form:

     A*V = B*V*E

   The matrices A and B are destroyed.

*/

void c8_geigen(c8 *a, int as, int an,
               c8 *b, int bs,
               c8 *v, int vs, r8 *eigenvals)
  {
   c8 *aband;
   int bw;
   r8 t0;
   r8 mflops;

// Check memory
   printf("Memory summary before eigensolve\n");
   cmm_memstat(NULL);

// Cholesky factorization of B
   t0 = cmm_clock();
   cmm_c8_cholesky(b, bs, an, an);
   t0 = cmm_clock() - t0;
   mflops = 1.3333e-6 * DBLE(an)* DBLE(an)* DBLE(an);
   printf("chol fact time %lf mflops %lf MFlop/sec %lf\n",t0, mflops, mflops/t0);


// Reduce A to ordinary form, in place
   t0 = cmm_clock();
   cmm_c8_tsolve_la( b, bs, an, an,
                 a, as, an, an,
                 a, as, an, an);
   cmm_c8_tsolve_alt(b, bs, an, an,
                 a, as, an, an,
                 a, as, an, an);
   t0 = cmm_clock() - t0;
   mflops = 8.e-6 * DBLE(an)* DBLE(an)* DBLE(an);
   printf("triangle time %lf mflops %lf MFlop/sec %lf\n",t0, mflops, mflops/t0);


// Reduce to condensed form
   bw = 2;
   aband = cmm_alloc(sizeof(c8)*an*(2*bw+1));
   t0 = cmm_clock();
   c8_condense(a, as, an, bw, aband);
   t0 = cmm_clock() - t0;
   mflops = 5.3333e-6 * DBLE(an)* DBLE(an)* DBLE(an);
   printf("condense time %lf mflops %lf MFlop/sec %lf\n",t0, mflops, mflops/t0);


// Compute eigenpairs of condensed form
   c8band_eigen(an, bw, aband, 0, an, eigenvals, v, vs);
   printf("c8band time %lf mflops %lf mflop/sec %lf\n",
        c8band_eigen_time(), c8band_eigen_mflops(),
        c8band_eigen_mflops()/(1.e-12 + c8band_eigen_time()));

// Free temp storage
   cmm_free(aband);

// Transform eigenvectors back to dense system
   t0 = cmm_clock();
   c8_house_back(a, as, an, an, bw, v, vs, an, an);
   t0 = cmm_clock() - t0;
   mflops = 8.e-6 * DBLE(an)* DBLE(an)* DBLE(an);
   printf("hback time %lf mflops %lf MFlop/sec %lf\n",t0, mflops, mflops/t0);

// Transform eigenvectors to generalized system, in place
   t0 = cmm_clock();
   cmm_c8_tsolve_lta(b, bs, an, an,
                 v, vs, an, an,
                 v, vs, an, an);
   t0 = cmm_clock() - t0;
   mflops = 4.e-6 * DBLE(an)* DBLE(an)* DBLE(an);
   printf("back tri time %lf mflops %lf MFlop/sec %lf\n",t0, mflops, mflops/t0);


// Check memory
   printf("Memory summary after eigensolve\n");
   cmm_memstat(NULL);

  }
