#include "cmm.h"
#include "cmmrand.h"

// Need to add some mean and correlation tests

int main(int argc, char **argv)
    {
     int i,n;
     r8 v;
     r8 v_1;
     r8 v_2;
     r8 v_3;
     r8 v_4;

     cmm_rand_init_seed(55978,0);

     v_1 = 0.;
     v_2 = 0.;
     v_3 = 0.;
     v_4 = 0.;
     n = 1000000;

     for(i=0;i<n;i++)
       {
        v = cmm_rand();
        v_1 = v_1 + v;
        v_2 = v_2 + v*v;
        v_3 = v_3 + v*v*v;
        v_4 = v_4 + v*v*v*v;
//        printf("%d %le\n",i,v);
       }

     printf("moments 1 %le 2 %le 3 %le 4 %le\n",
        v_1/DBLE(n),v_2/DBLE(n),v_3/DBLE(n),v_4/DBLE(n));

     return 0;
    }
