#include "cmm.h"

// Test cmm_clock() and cmm_date()
static r8 tx(int j)
  {
   r8 x;

   x = DBLE(j)/1000000;
   return sin(x)*exp(-x);
  }

static void t1(void)
  {
   r8 t0;
   int i,j;
   r8 s;

   cmm_print("Current date and time %s\n",cmm_date());
   cmm_print("doing a few seconds of cpu work\n");
   t0 = cmm_clock();

   for(j=0;j<10;j++)
    {
     s = 0.;
     for(i=0;i<1000000;i++)
       s += tx(i);
    }

   cmm_print("Elapsed time %lf\n",cmm_clock() - t0);

  }

int main(int argc, char **argv)
  {
   t1();
   return 0;
  }
