#define CMM_SHORTHAND

#include "cmm.h"
#include "cmmblas.h"
#include "cmmatrix.h"

/*

   Condense a Hermitian matrix to banded form.

   The matrix A is overwritten with the Householder
   vectors used to build the banded matrix.

   The bandwidth b must be fairly small (b <= 4)

*/

#ifdef DO_CLOCKS
#define MHZ  698.e6
void pcc_init();
int pcc_delta();
#endif

void c8_condense(c8 *a, int as, int n, int b, c8 *aband)
  {
   c8 *tmp;
   c8 *tmp1;
   int ntmp;
   c8 *w;
   c8 atmp[4*4];
   int i, j, k;
   int jlast;
   r8 t0;

#ifdef DO_CLOCKS
   int *clocks;
   int *clkp;
   int nstep;
#endif

   t0 = cmm_clock();

#ifdef DO_CLOCKS
#define NUM_CLOCKS 8
   nstep = 0;
   clocks = cmm_alloc(sizeof(int)*n*NUM_CLOCKS);
   clkp = clocks;
   pcc_init();
#endif

   // Allocate temp storage
   ntmp = 4*n*b;
   tmp = cmm_alloc(sizeof(c8) * ntmp);

   tmp1 = tmp;
#define TMP_ALLOC(ptr, l) ptr = tmp1; tmp1 += l; ntmp -= l;

   TMP_ALLOC(w, n*b);   

   jlast = 0;
   for(j=0;j<n-b-1;j += b)  // Loop over blocks
     {
      int icol = j+b;
      c8 *a0 = C8_M_ELEMENT(a, as, icol, j);
      int ni = n-icol;
      int r;

#ifdef DO_CLOCKS
      clkp[0] = pcc_delta();
#endif

      r = MIN(ni, b);
      jlast = j + r;

//      cmm_warn("j %d icol %d r %d jlast %d\n",j,icol,r,jlast);

#define C8_BAND_OFFSET(bw, icol, irow) \
  ( (bw) + (icol) - (irow) + (irow)*(2*bw+1))

   // Put away unmodified diagonal/subdiagonal elements of A
      for(k=0;k<r;k++)
        for(i=k;i<b;i++)
          {
           c8 aik = C8_M_ELEMENT(a, as, j+i,j+k)[0];
           int offik = C8_BAND_OFFSET(b, j+i,j+k);
           int offki = C8_BAND_OFFSET(b, j+k,j+i); 

           aband[offik] = aik;

           aik.i = -aik.i;
           if(i != k)
             aband[offki] = aik;
          }

   // Build Householder vectors for this block
      c8_build_house(a0, as, ni, r, atmp, r);

   // Put away modified elements of A (now stored in atmp)
      for(k=0;k<r;k++)
        for(i=0;i<k+1;i++)
          {
           c8 aik = atmp[i+k*r];
//           c8 aik = C8_M_ELEMENT(a, as, j+b+i,j+k)[0]; // for debugging
           int offik = C8_BAND_OFFSET(b, j+b+i,j+k);
           int offki = C8_BAND_OFFSET(b, j+k,j+b+i); 

           aband[offik] = aik;

           aik.i = -aik.i;
           if(i+b != k)
             aband[offki] = aik;
          }

      if(ni == 0)
        break;  // we are done

#ifdef DO_CLOCKS
      clkp[1] = pcc_delta();
#endif

   // Build WY from Householder vectors
      c8_build_wy(a0, as, ni, r,    // Y
                   w, ni);          // W

#ifdef DO_CLOCKS
      clkp[2] = pcc_delta();
#endif

   // Update remainder of A
      c8_update_wy(a0,  as, ni, r,  // Y
                    w,  ni, ni, r,  // W
                    C8_M_ELEMENT(a, as, j+b, j+b), as, ni, ni,  // A submatrix
                    tmp1, ntmp
#ifdef DO_CLOCKS
, clkp + 3
#endif
                  );

#ifdef DO_CLOCKS
      clkp[5] = pcc_delta();
#endif

#ifdef DO_CLOCKS
      nstep++;
      clkp += NUM_CLOCKS;
#endif

     }

   // Put away final elements of A
    {
     int r;

     j = jlast;
     r = n-j;

//     cmm_warn("Put away r %d\n",r);

     for(k=0;k<r;k++)
       for(i=k;i<r;i++)
        {
         c8 aik = C8_M_ELEMENT(a, as, j+i,j+k)[0];
         int offik = C8_BAND_OFFSET(b, j+i,j+k);
         int offki = C8_BAND_OFFSET(b, j+k,j+i); 

         aband[offik] = aik;

         aik.i = -aik.i;
         if(j+i != j+k)
           aband[offki] = aik;
        }

      for(k=0;k<r;k++)
        c8_zero(r-k+1, C8_M_ELEMENT(a, as, jlast+k-1, jlast+k));
     }

   // Free temp storage
   cmm_free(tmp);

   t0 = cmm_clock() - t0;

#ifdef DO_CLOCKS
   {
    r8 sum;
    sum = 0;
    clkp = clocks;
    for(i=0;i<nstep;i++)
      {
       printf(" %d    %d    %d  %d  %d  %d\n",clkp[0],clkp[1],clkp[2],clkp[3],clkp[4],clkp[5]);
       sum += DBLE(clkp[3]);
       clkp += NUM_CLOCKS;
      }
    printf("time in update_wy %lf total %lf\n",sum/MHZ,t0);
    cmm_free(clocks);
   }
#endif
  }
