/*! \file myzgemm.c
    \brief reference implementation.

     This implementation is intended to be clear and correct,
     although it is probably not very efficient.

*/

#include "cmm.h"

void myzgemm_nn(c8 fc, c8 fa,
               c8 *a, int as, int am, int an,
               c8 *b, int bs, int bm, int bn,
               c8 *c, int cs, int cm, int cn
              )
  {
   int i, j, k;

   if( (an != bm) || (am != cm) || (bn != cn))
     {
      if(an!=bm)
        printf("an!=bm\n");

      if(am!=cm)
        printf("am!=cm\n");

      if(bn!=cn)
        printf("bn!=cn\n");

      cmm_warn("error in myzgemm_nn\n %d %d %d\n %d %d %d\n %d %d %d\n",
        as,am,an,bs,bm,bn,cs,cm,cn);
      cmm_exit(0);
     }

   for(j=0;j<cn;j++)
     {
      c8 *cp = c + j*cs;
      for(i=0;i<cm;i++)
        {
         c8 *ap = a + i;
         c8 *bp = b + j*bs;
         r8 tr, ti, cr, ci;
         tr = 0.;
         ti = 0.;
         for(k=0;k<an;k++)
          {
           tr += ap->r * bp->r - ap->i * bp->i;
           ti += ap->r * bp->i + ap->i * bp->r;
           ap += as;
           bp ++;
          }
        cr = cp[i].r;
        ci = cp[i].i;

        cp[i].r = cr*fc.r - ci*fc.i + tr * fa.r - ti * fa.i;
        cp[i].i = ci*fc.r + cr*fc.i + tr * fa.i + ti * fa.r;

        }
     }
  }

void myzgemm_hn(c8 fc, c8 fa,
               c8 *a, int as, int am, int an,
               c8 *b, int bs, int bm, int bn,
               c8 *c, int cs, int cm, int cn
              )
  {
   int i, j, k;

   if( (am != bm) || (an != cm) || (bn != cn))
     {
      cmm_warn("error in myzgemm_hn\n %d %d %d\n %d %d %d\n %d %d %d\n",
        as,am,an,bs,bm,bn,cs,cm,cn);
      cmm_exit(0);
     }

    for(j=0;j<cn;j++)
     {
      c8 *cp = c + j*cs;
      for(i=0;i<cm;i++)
       {
        r8 tr, ti, cr, ci;
        c8 *ap = a + i*as;
        c8 *bp = b + j*bs;

        tr = 0.;
        ti = 0.;
        for(k=0;k<am;k++)
          {
           tr += ap->r * bp->r + ap->i * bp->i;
           ti += ap->r * bp->i - ap->i * bp->r;
           ap++;
           bp++;
          }

        cr = cp[i].r;
        ci = cp[i].i;

        cp[i].r = cr*fc.r - ci*fc.i + tr * fa.r - ti * fa.i;
        cp[i].i = ci*fc.r + cr*fc.i + tr * fa.i + ti * fa.r;
       }
     }
  }

void myzgemm_nh(c8 fc, c8 fa,
               c8 *a, int as, int am, int an,
               c8 *b, int bs, int bm, int bn,
               c8 *c, int cs, int cm, int cn
              )
  {
   int i, j, k;

   if( (an != bn) || (am != cm) || (bm != cn))
     {
      cmm_warn("error in myzgemm_nh\n %d %d %d\n %d %d %d\n %d %d %d\n",
        as,am,an,bs,bm,bn,cs,cm,cn);
      cmm_exit(0);
     }

   for(j=0;j<cn;j++)
    {
     c8 *cp = c + j*cs;
     for(i=0;i<cm;i++)
       {
        r8 tr, ti, cr, ci;
        c8 *ap = a + i;
        c8 *bp = b + j;

        tr = 0.;
        ti = 0.;
        for(k=0;k<an;k++)
          {
           tr += ap->r * bp->r + ap->i * bp->i;
           ti += ap->i * bp->r - ap->r * bp->i;
           ap += as;
           bp += bs;
          }

        cr = cp[i].r;
        ci = cp[i].i;

        cp[i].r = cr*fc.r - ci*fc.i + tr * fa.r - ti * fa.i;
        cp[i].i = ci*fc.r + cr*fc.i + tr * fa.i + ti * fa.r;
       }
     }

  }
