/*

   Test memory allocation/free

*/

#include "cmm.h"

int main(int argc, char **argv)
  {
   void *p;

   cmm_preallocate(10*1024);

   p = cmm_alloc(12345);

   cmm_free(p);

   cmm_memstat(NULL);

   return 0;
  }
