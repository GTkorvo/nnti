#include "cmm.h"
#include "cmmrand.h"

// seeds:

static unsigned long seed0;
static unsigned long seed1;

// recommended by Paul Coddington, paulc@npac.syr.edu
#define A1 40014
#define C1 0
#define M1 2147483563

#define A2 40692
#define C2 0
#define M2 2147483399

// on a pIII, 800Mhz, this routine can compute 10^6 random numbers in .3 sec

r8 cmm_rand()
  {
   r8 xr;

   if(seed1 == 0)
     cmm_rand_init_seed(398745, 0);

   seed0 = (A1*seed0 + C1)%M1;
   seed1 = (A2*seed1 + C2)%M2;

   xr = DBLE(seed0) + DBLE(seed1);
   if(xr >= DBLE(M1)) // because M1 = MAX(M1, M2)
     xr -= DBLE(M1);
   xr = xr/DBLE(M1);
   return xr;
  }

void cmm_randv(int n, r8 *v)
  {
   int i;

   if(seed1 == 0)
     cmm_rand_init_seed(398745, 0);

   for(i=0;i<n;i++)
     {
      r8 xr;

      seed0 = (A1*seed0 + C1)%M1;
      seed1 = (A2*seed1 + C2)%M2;

      xr = DBLE(seed0) + DBLE(seed1);
      if(xr >= DBLE(M1)) // because M1 = MAX(M1, M2)
        xr -= DBLE(M1);
      xr = xr/DBLE(M1);

      v[i] = xr;
     }
  }

void cmm_rand_init_seed(unsigned long seed0v, unsigned long seed1v)
  {
   seed0 = seed0v;
   if(seed0 == 0)
     seed0 = 954837833;

   seed1 = seed1v;
   if(seed1 == 0)
     seed1 = 3987583;
  }

void  cmm_rand_get_seed(unsigned long *seed0p, unsigned long *seed1p)
  {
   *seed0p = seed0;
   *seed1p = seed1;
  }

//  Park and Miller RNG
r8 r8_rand_park(int *seed)
  {
   static r8 aa=16807.;
   static r8 mm=2147483647.;
   static r8 t;
   int j;

   t = DBLE(*seed);
   t = t*aa;
   j = INT(t/mm);

   t -= j*mm;

   *seed = INT(t);
   return t/mm;
  }
