/*! \file math.c
    \brief Simple math functions.
    \ingroup base
*/

#include "cmm.h"

/*!
  \brief test if value is in a neighborhood.
  \param x Test value.
  \param y Center of neighborhood.
  \param epsilon Size of neighborhood.

  \return True if \f$ \left| x - y \right| \le \epsilon \f$, otherwise False.
*/
Boolean cmm_nbrhood(r8 x, r8 y, r8 eps)
    {
     r8 t = x-y;
     if(t < -eps)
       return False;
     if(t > eps)
       return False;
     return True;
    }

/*!
  \brief Compute square of value.
  \param x value.
  \return \f$ x^2 \f$.
*/
r8 cmm_sq(r8 x)
    {
     return x*x;
    }

/*!
  \brief Compute integer power of value.
  \param x value.
  \param n power.
  \return \f$ x^n \f$.

  Works for zero and negative integer powers.

*/
r8 cmm_ipow(r8 x, int n)
    {
     int nn;
     r8 t;

     nn = n;
     if(nn<0) nn = -nn;
     t = 1.;
     while(nn)
       {
        t *= x;
        nn--;
       }
     if(n >= 0)
      return t;
     else
      return 1./t;
    }

/*!
  \brief Compute integer power of 10.
  \param n power.
  \return \f$ 10^n \f$.

  Works for zero and negative integer powers.

*/
r8 cmm_ipow10(int n)
    {
     r8 t;

     if(n == 0)
       return 1.;

     if(n < 0)
       return 1./cmm_ipow10(-n);

     t = 1.;
     while(n>0)
       {
        t *= 10.;
        n--;
       }

     return t;
    }

// Using 64 bit IEEE double precision, factorials
// are exact up to about 16! to 17!, assuming
// 48 bits of fraction available. 70! is about 10^100,
// so we don't go beyond that.

#define MAX_FACTORIAL 70
static int nfactorials;
static r8 factorials[MAX_FACTORIAL+1];

static void setup_factorials()
    {
     int i;
     r8 t;

     t = 1.;
     for(i=0;i<=MAX_FACTORIAL;i++)
       {
        factorials[i] = t;
        t *= DBLE(i+1);
       }
     nfactorials = MAX_FACTORIAL;
    }

/*!
   \brief Compute factorials.
   \param n value.
   \return \f$ n! \f$.

*/
r8 cmm_factorial(int n)
    {
     if(nfactorials == 0)
       setup_factorials();

     if(n < 0 || n > MAX_FACTORIAL)
       {
        cmm_fatal("Attempt to compute factorial of %d\n",n);
       }

     return factorials[n];
    }

/*!
   \brief Compute combination.
   \param n value.
   \return \f$ {n!}\over{m! (n-m)!} \f$.

*/
r8 cmm_combination(int n, int m)
    {
     if(nfactorials == 0)
       setup_factorials();

     return factorials[n]/(factorials[m]*factorials[n-m]);
    }

/*

   Set up permutation lists

*/
static int *p;
static int *plist[10];

void gen_permutations(int n)
    {
     int i,j,k;
     int ifact;
     int nfact;

     if(p != NULL)
       return;

// allocate lists
     ifact = 0;
     for(i=0;i<=n;i++)
       {
        ifact = ifact + i*cmm_factorial(i);
       }

     p = cmm_alloc(sizeof(int)*ifact);

     ifact = 0;
     nfact = 1;
     for(i=1;i<=n;i++)
       {
        plist[i] = p + ifact;
        nfact *= i;
        ifact = ifact + nfact*i;
       }

     p[0] = 0;

// build list i from list i-1
     ifact = 1;
     for(i=1;i<=n;i++)
       {
        int *pout = plist[i];
        int *pin  = plist[i-1];

        for(k=0;k<ifact;k++)
          {
           for(j=0;j<i;j++)          // generate all perms where i-1 is in location j.
              {
               int ii;
               for(ii=0;ii<i;ii++)
                 {
                  int kk;
                  if(ii < j)
                    kk = pin[ii];
                  else if(ii == j)
                    kk = i-1;
                  else
                    kk = pin[ii-1];
                  pout[ii] = kk;
                 }
                pout = pout + i;
              }
           pin = pin + i-1;
          }
        ifact = ifact * i;
       }

    }

//?
int *get_permutation(int n, int k)
    {
     return plist[n] + k*n;
    }

