
#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmband.h"
#include "cmmrand.h"

/*

    This routine computes a contiguous subset of all the eigenpairs
    of a banded Hermitian matrix A.

    The approach is to use bisection to isolate eigenvalues
    or clusters and then to use inverse iteration based on QR
    factorization with a complex shift to compute eigenvectors
    and refine the eigenvalues.

    On output the eigenvectors should be orthogonal to "working
    accuracy", or something like it.

    n      -- matrix order
    b      -- bandwidth (note that a typical row contains 2*b+1 nonzero elements
    A      -- pointer to band matrix data
    kstart -- starting eigenpair (0 based index)
    kcount -- number of eigenpairs to compute
    e      -- pointer to output list of eigenvalues
    v      -- pointer to output eigenvectors

    Band matrix data is stored row-wise. Each row is of length 2*b+1 and there
    are n rows. The band matrix data is therefore dimensioned (2*b+1)*n complex
    numbers.

    The output eigenvectors are stored column wise as a complex n by n array.
    Column i is the eigenvector associated with eigenvalue e[i]. Eigenvalues
    are sorted smallest first.

TO DO:

    add code for case where cluster size is > 1

*/
void c8band_eigen(int n, int b, r8 *A, int kstart, int kcount, r8 *ev, r8 *v, int vstride);

// Accumulated stats:
static int zeigen_n;
static int zeigen_count_oracle;
static int zeigen_count_qr;
static int zeigen_count_invit;
static int zeigen_count_q_apply;
static int zeigen_count_matvec;
static r8  zeigen_time;
static r8  zeigen_tr, zeigen_e2;
static r8  zeigen_tol;
extern int oracle_count;

// Return MFlops in most recent call
r8 c8band_eigen_mflops()
  {
   r8 mflops = 0;

  // a divide counts 3 and sqrt counts 10 (somewhat arbitrary)

   mflops  = DBLE(zeigen_n*zeigen_count_oracle) * DBLE(129 + 6 + 10);
   mflops += DBLE(zeigen_n*zeigen_count_qr)     * DBLE(220 + 3*3 + 10*3);
   mflops += DBLE(zeigen_n*zeigen_count_q_apply)* DBLE(90);
   mflops += DBLE(zeigen_n*zeigen_count_invit)  * DBLE( 45+3);
   mflops += DBLE(zeigen_n*zeigen_count_matvec) * DBLE(40);

   return mflops * 1.e-6;  
  }

// Return time used in most recent call
r8 c8band_eigen_time()
  {
   return zeigen_time;
  }

// A cluster is defined as a group of eigenvalues in
// an interval [e0,e1] with gaps to nearest clusters
// of g0, g1. The isolation routine builds the clusters
// using the oracle function.

typedef struct
  {
    r8 e0, e1;    // lower, upper interval limits
    r8 g0, g1;    // lower, upper gaps
    r8 q;         // quality
    int k0, k1;   // starting eigenvalue, ending eigenvalue
    int state;
    int indx;
    void *prev;
    void *next;
  } BCluster;

// build the clusters using the banded Wilkinson algorithm for the oracle
static int bisect_isolate(int n, int b, r8 *A, r8 emin, r8 emax, int k1, int kc,
      int kmax, r8 sigma, BCluster *clusters);

// build the eigenvectors and eigenvalues for a cluster using the banded complex shift QR method
static void invit_qr(int n, int b, r8 *A, BCluster *cluster, r8 *v, r8 *ev,
     r8 *q, r8 *r, r8 *xt);

void c8band_eigen(int n, int b, r8 *A, int kstart, int kcount, r8 *ev, r8 *v, int vstride)
  {
   r8 emin, emax;
   BCluster *clusters;
   int nclusters;
   int i;
   int kmax;
   r8 sigma;
   r8 t0;  // timer
   r8 *q, *r, *xt;    // temp working storage for invit_qr
   r8 *etmp;          // temp storage for cluster eigenvalues
   r8 *vtmp;          // temp storage for cluster eigenvectors
   int invit_count;
   int neig;

   // Initialize stats
   zeigen_n                = n;
   zeigen_count_oracle     = 0;
   zeigen_count_qr         = 0;
   zeigen_count_invit      = 0;
   zeigen_count_q_apply    = 0;
   zeigen_time  = cmm_clock();      // total time

   t0 = cmm_clock();

/* 

   Using bisection, isolate clusters of eigenvalues. The parameter
   kmax specifies the maximum number of eigenvalues in a cluster.
   The parameter sigma defines the target quality factor for 
   each cluster, which is the ratio of the cluster interval to 
   the minimum of the upper and lower gaps on each side. The result is
   stored in the vectors el,eu which define the interval and ke,
   which defines the number of eigenvalues in the interval after
   isolation. The tolerance zeigen_tol is used to converge eigenvectors.

*/

   zeigen_tol = 1.e-15;
//   kmax = 2;
     kmax = 1;
//     sigma = .003;
     sigma = .003;
//     sigma = 1.e-5;

  // Allocate O(n) temp storage. Note that we must allocate
  // a temp block of size n*kmax to build eigenvectors for
  // a cluster. The reason for this is that the requested
  // block of eigenvectors (kstart to kstart + kcount-1) may
  // be part of a cluster. This also means that construction
  // of the eigenpairs must be determinate so that two calls
  // which split a cluster will give the same results for
  // eigenpairs even on different processors.

   clusters = cmm_alloc(sizeof(BCluster)*(n+2));
   q   = c8band_new(n, b, 0, NULL);
   r   = c8band_new(n, 0, 2*b+1, NULL);

   xt  = cmm_alloc(sizeof(r8)*n*2);   // ?

   etmp = cmm_alloc(sizeof(r8)*kmax);
   vtmp = cmm_alloc(sizeof(c8)*n*kmax);

  // Compute trace, e2 from the initial matrix
   c8band_trace(n, b, A, &zeigen_tr, &zeigen_e2);
//   printf("band eigensolver begin, tr %23.18lf e2 %23.18lf (direct)\n",
//      zeigen_tr, zeigen_e2);

  // Compute range of eigenvalues using Gershgorin
   c8band_gershgorin(n, b, A, &emin, &emax);  // compute range of eigenvalues

  // Enlarge the Gershgorin interval a little bit
   {r8 de = emax-emin;
    emin = emin - .01*(de);
    emax = emax + .01*(de);
   }

  // Build the clusters
   oracle_count = 0;
   nclusters = bisect_isolate(n, b, A, emin, emax, kstart, kcount, kmax, sigma, clusters);
   zeigen_count_oracle = oracle_count;

   t0 = cmm_clock() - t0;
//   printf("bisect time %lf  ni %d kstart %d kcount %d\n",t0,nclusters,kstart,kcount);

   t0 = cmm_clock();

//   Now we use inverse iteration to compute the eigenvectors in
//   each cluster.

   for(i=0;i<nclusters;i++)
     {
      int j,k;
      BCluster *c = clusters+i;

      if(c->k1 < kstart || c->k0 >= kstart+kcount)
        continue;

    // number of eigenvalues in the cluster
      k = c->k1-c->k0 + 1;

    // build cluster of eigenvectors using complex shift qr inverse iteration.
      invit_qr(n, b, A, c, vtmp, etmp, q, r, xt);

    // copy eigenpairs to output
      for(j=0;j<k;j++)
        {
         int ieig;

         ieig = c->k0 + j;
         if(ieig < kstart || ieig >= kstart + kcount)
           continue;

         r8_copy(2*n, vtmp + 2*j*n, 1, v+2*vstride*(ieig-kstart), 1);
         ev[ieig-kstart] = etmp[j];

        }

     }

// Free temp storage
   cmm_free(q);
   cmm_free(r);
   cmm_free(xt);
   cmm_free(vtmp);
   cmm_free(etmp);
   cmm_free(clusters);

// Final stats
   zeigen_tr = 0.;
   zeigen_e2 = 0.;
   for(i=0;i<kcount;i++)
     {
      zeigen_tr += ev[i];
      zeigen_e2 += ev[i]*ev[i];
     }

   t0 = cmm_clock() - t0;

   zeigen_time  = cmm_clock() - zeigen_time;

  }

/*

   Isolate a group of kc eigenvalues starting at k1.

   Return value:

          ni -- number of intervals isolated

   Return intervals:

          el -- lower limit on each interval
          eu -- upper limit on each interval
          ke -- number of eigenvalues in each interval
          kl -- number of eigenvalues less than interval

*/
static int bisect_isolate(int n, int b, r8 *A, r8 emin, r8 emax, int kstart,
      int kc, int kmax, r8 sigma, BCluster *clusters)
    {
     int i,ni;
     int *ku;
     r8 *tmp;
     BCluster *c;

     if(kc == 0)
       return 0;

     c = clusters;

     tmp = cmm_alloc(sizeof(r8)*2*(2*b+1)*n*4);  // I don't actually need this

// Set up the initial interval
     ni = 1;

     c->e0 = emin;
     c->e1 = emax;
     c->g0 = 1.e10 * (emax - emin);
     c->g1 = 1.e10 * (emax - emin);
     c->k0 = c8band_oracle(n, b, A, emin, tmp);
     c->k1 = c8band_oracle(n, b, A, emax, tmp);
     c->state = 0;

     for(i=0;i<n;i++)
       {
        clusters[i].indx = i;
       }

// Loop until all the clusters are smaller than kmax with quality less than sigma (small is good)
     while(1)
       {
        r8 qmax;
        int klg,knew;
        r8 ex;

      // Loop over clusters and get the worst one. Skip clusters that are finished.

      // This is a stupid implementation. A smart one would use a priority queue.

        c = NULL;
        qmax = 0.;
        klg  = kmax;
        for(i=0;i<ni;i++)
          {
           BCluster *c1;
           int k;
           r8 q;

           c1 = clusters + i;

//           printf("look at cluster %d of %d state %d  [%lf %lf] [%d %d]\n",
//             c1->indx,ni,c1->state,c1->e0,c1->e1,n-c1->k0,n-c1->k1);

           if(c1->state) // skip good cluster
             continue;

           if(n - c1->k1 < kstart - kmax*2)  // skip irrelevant cluster
             continue;

           if(n - c1->k0 >= kstart + kc + kmax*2) // skip irrelevant cluster
             continue;

           k = c1->k0 - c1->k1; // number of eigenvalues in the cluster

           if(c1->g0 == 0. || c1->g1 == 0.)
              q = 1.e18;
           else
              q = (c1->e1 - c1->e0)/MIN(c1->g0, c1->g1);

//           printf("k %d q %le\n",k,q);

           if(k <= kmax && q <= sigma)  // skip good cluster
             {
              c1->state = 1;
              continue;
             }

           if(k > klg)  // split the largest cluster
             {
              klg = k;
              c = c1;
             }

           else if(q > qmax)  // or the one with worst quality
             {
              qmax = q;
              c = c1;
             }

          }

   //      if(c == NULL)
   //        printf("all clusters ok qmax %le klg %d\n",qmax,klg);

         if(c == NULL)  // clusters are all ok, we are done.
           break;

   //      printf("bad cluster %d of %d   [%lf %lf] [%d %d]\n",
   //        c->indx,ni,c->e0,c->e1,n-c->k0,n-c->k1);


       // At this point c is the worst cluster. Use the oracle function to either
       // split it or improve the quality.

         ex = .5*(c->e0 + c->e1);
         knew = c8band_oracle(n, b, A, ex, tmp); // number of eigenvalues > ex

         if(knew == c->k0)  // there are no eigenvalues between e0 and ex
           {
            c->e0 = ex;
           }
         else if(knew == c->k1) // there are no eigenvalues between ex and e1
           {
            c->e1 = ex;
           }
         else // split the cluster
           {
            BCluster *c1;

            for(i=ni;i>c->indx;i--)
              {
               clusters[i] = clusters[i-1];  // make room
               clusters[i].indx = i;
              }

            ni = ni + 1;

            c->e1 = ex;
            c->k1 = knew;  // c->g0 is correct. c->g1 will be updated below

            c1 = c + 1;
            c1->e0 = ex;
            c1->k0 = knew;
            c1->g0 = 0.;   // c1->g1 is correct.
           }

         // Recompute gaps for c and its neighbors.

         {
          BCluster *cprev;
          BCluster *cnext;

          cprev = (c->indx>0)?(c-1):NULL;
          if(cprev)
           {
            cprev->g1 = c->e0 - cprev->e1;
            c->g0 = cprev->g1;
           }

          cnext = (c->indx<ni-1)?(c+1):NULL;
          if(cnext)
           {
            cnext->g0 = cnext->e0 - c->e1;
            c->g1 = cnext->g0;
           }

         }

        }  // end iteration loop

   // update clusters
      for(i=0;i<ni;i++)
        {
         c = clusters + i;
         c->k0 = n-c->k0;
         c->k1 = n-c->k1-1;
        }

// debug, print clusters

/*
     {
      r8 tr,e2;
      printf("isolated clusters %d\n",ni);
      tr = 0.;
      e2 = 0.;
      for(i=0;i<ni;i++)
        {
         BCluster *c = clusters + i;
         r8 ev = .5*(c->e0 + c->e1);
         int nc = c->k1 - c->k0 + 1;

         printf("cluster %d size %d [%lf %lf] [%d %d]\n",c->indx,nc,c->e0,c->e1,c->k0,c->k1);

         tr += ev*DBLE(nc);
         e2 += ev*ev*DBLE(nc);
        }

     }
*/
     cmm_free(tmp);

   // Return number of clusters built
      return ni;
     }

// Handle cluster with single eigenvalue
static void  invit_qr_1(int n, int b, r8 *A, BCluster *c, r8 *vtmp, r8 *etmp, r8 *q, r8 *r, r8 *xt);

// Handle cluster with multiple eigenvalues
static void  invit_qr_general(int n, int b, r8 *A, BCluster *c, r8 *vtmp, r8 *etmp, r8 *q, r8 *r, r8 *xt)
  {
   cmm_warn("Cannot deal with multiple eigenvalues");
   cmm_exit(0);
  }

// Driver for inverse iteration
static void  invit_qr(int n, int b, r8 *A, BCluster *c, r8 *vtmp, r8 *etmp, r8 *q, r8 *r, r8 *xt)
  {
   int i,k;
   c8 z;

// Build shifted qr factorization
   z.r = -.5*(c->e0 + c->e1);
   z.i =   .5*(c->e1 - c->e0);

   c8band_qr_shifted(n, b, b, A, z, q, r);
   zeigen_count_qr++;

// inverse iteration in factored space
   k = c->k1 - c->k0 + 1;

   // printf("invit %d k %d e_est %lf de %lf\n",
   //    c->k0, k, z.r, z.i);

   switch(k)
    {
    case 1:  // cluster with only one eigenvalue
     invit_qr_1(n, b, A, c, vtmp, etmp, q, r, xt);
     break;

    default: // cluster with multiple eigenvalues
     invit_qr_general(n, b, A, c, vtmp, etmp, q, r, xt);
     break;
    }

// Transform eigenvectors back to band basis and build eigenvalue estimates.
   for(i=0;i<k;i++)
     {
      r8 *vp = vtmp + 2*n*i;

     // Apply Q to get eigenvector back in band space

      c8band_q_apply(n, b, q, vp);
      zeigen_count_q_apply++;

     // Compute the eigenvalue from the Ritz estimate

      r8_zero(2*n,xt);  // this is important, should it be?
      c8band_matvec(n, b, b, A, vp, xt);
      zeigen_count_matvec++;

      z = c8_hdot(n, (c8 *) vp, (c8 *) xt);
   //   printf("eigenvalue from Ritz estimate %le %le\n",z.r,z.i);
      etmp[i] = z.r;

      {
       int j;
       r8 e_est = z.r;
       r8 t;

       t = 0.;
       for(j=0;j<n;j++)
         {
          r8 xr, xi;

          xr = xt[2*j]   - e_est*vp[2*j];
          xi = xt[2*j+1] - e_est*vp[2*j+1];

          t += xr*xr + xi*xi;
          // printf("j %d %18.16lf %18.16lf\n",j,xr, xi);
         }

   //    printf("|A*v - eig*v| %le eig = %20.18lf\n",t,e_est);
   //    for(j=0;j<n;j++)
   //       printf("j %d %18.16lf %18.16lf\n",j,xt[2*j]-e_est*vp[2*j],xt[2*j+1]-e_est*vp[2*j+1]);
       }

     }

  }

// This is the most common case in practice
static void  invit_qr_1(int n, int b, r8 *A, BCluster *c, r8 *vtmp, r8 *etmp, r8 *q, r8 *r, r8 *xt)
  {
   int j,it;
   r8 tol = zeigen_tol;
   r8 xn;

  // Create the initial vector

#ifdef RAND_INITIALIZE
      cmm_r8_randv(2*n, vtmp);
#else
   for(j=0;j<n;j++)
     {
      vtmp[2*j]   = DBLE(j*c->k1 + 1.);
      vtmp[2*j+1] = 0.;
     }
#endif
  // normalize it

   xn = c8_norm(n, (c8 *) vtmp);
   r8_scale(2*n, 1./xn, vtmp);

  // inverse iteration
   for(it=0;it<15;it++)
      {
       c8 zdot;

       c8band_invit(n, b, b, r, vtmp, xt);
       zeigen_count_invit++;

       xn = c8_norm(n, (c8 *) xt);
       r8_scale(2*n, 1./xn, xt);

// I am not sure this is the correct way to do this. Perhaps
// |vnew - vold| would be better, although then I would have
// to normalize the phase also.

       zdot = c8_hdot(n, (c8 *) xt, (c8 *) vtmp);

       // printf("xn %le zdot %le %le\n",xn,zdot.r-1., zdot.i);

       xn = zdot.r-1.;

       xn = ABS(xn) + ABS(zdot.i);

       c8_copy(n, (c8 *) xt, 1, (c8 *) vtmp, 1);

       if(xn < tol)
          break;

      }

   // one extra iteration
   c8band_invit(n, b, b, r, vtmp, xt);
   c8_copy(n, (c8 *) xt, 1, (c8 *) vtmp, 1);
   zeigen_count_invit++;

   // final normalization
   xn = c8_norm(n, (c8 *) vtmp);
   r8_scale(2*n, 1./xn, vtmp);

  }

