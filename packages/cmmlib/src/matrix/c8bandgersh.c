#include "cmm.h"
#include "cmmblas.h"
#include "cmmband.h"

/*
    Compute estimates for range of eigenvalues of a banded
    matrix using Gershgorin discs.
 */
void c8band_gershgorin(int n, int b, r8 *A, r8 *emin, r8 *emax)
    {
     int i,j;
     r8 e0,e1;

     e0 =  1.e99;
     e1 = -1.e99;

     for(i=0;i<n;i++)
       {
        r8 z0, z1;
        r8 *trow = A + 2*(2*b+1)*i;

        z0 = trow[2*b];  // diagonal element is the center of the disc
        z1 = 0.;

        for(j=0;j<2*b+1;j++)  // assume zero padding
          {
           r8 tr, ti;

           if(j == b)
             continue; // skip diagonal element

           if(i+j-b >= n)
             continue;
           if(i+j-b < 0)
             continue;

           tr = trow[2*j];
           ti = trow[2*j+1];
           z1 = z1 + sqrt(tr*tr+ti*ti);
          }

        e0 = MIN(e0, z0-z1);
        e1 = MAX(e1, z0+z1);

//        printf("gersh i %d z0 %lf z1 %lf e0 %lf e1 %lf\n",i,z0,z1,e0,e1);
           
       }


     *emin = e0;
     *emax = e1;
    }
