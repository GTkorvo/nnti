#include "cmm.h"

// Test type info routines
static void t1()
  {
   printf("sizeof r8 %d =? %d\n",sizeof(r8),cmm_typesize(CMM_REAL));
   printf("sizeof c8 %d =? %d\n",sizeof(c8),cmm_typesize(CMM_COMPLEX));
  }

// Test cmm_print, cmm_warn, cmm_pause, 
static void t2()
  {
   cmm_print("test cmm_print: Hello, world.\n");
   cmm_pause("test cmm_pause");
   cmm_warn("test cmm_warn");
  }

int main(int argc, char **argv)
  {
   t1();
   t2();

   cmm_fatal("test cmm_fatal");

   printf("cmm_fatal: Error: Got here, and should not.\n");
   return 0;
  }
