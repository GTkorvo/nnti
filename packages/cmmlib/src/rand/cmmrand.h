#ifndef CMMM_RAND_H
#define CMMM_RAND_H 1

#ifdef __cplusplus
extern "C"
{
#endif                          /* __cplusplus */

r8   cmm_rand();
void cmm_randv(int n, r8 *v);
void cmm_rand_init_seed(unsigned long seed0, unsigned long seed1);
void cmm_rand_get_seed(unsigned long *seed0p,  unsigned long *seed1p);

#ifdef __cplusplus
}
#endif                          /* __cplusplus */

#endif
