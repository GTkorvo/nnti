/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

/*  This file is included only in hg_test, not in Zoltan. */

#include "hypergraph.h"
#include "dr_hg_readfile.h"

#define BUF_LEN 1000000


int hg_readfile (ZZ *zz, HGraph *hg, char *hgraphfile, int *base)
{
FILE *f;
int err = ZOLTAN_OK;
char errstr[200];
char *yo = "hg_readfile";

   ZOLTAN_TRACE_ENTER(zz, yo) ;

   Zoltan_HG_HGraph_Init(hg);

   f = fopen (hgraphfile, "r");
   if (!f) {
      sprintf(errstr, "ERROR...not able to open file %s!\n", hgraphfile);
      ZOLTAN_PRINT_ERROR (zz->Proc, yo, errstr);
      err = ZOLTAN_FATAL;
      goto End;
      }

   err = Zoltan_HG_Readfile (0, f, &hg->nVtx, &hg->nEdge, &hg->nInput,
    &hg->hindex, &hg->hvertex, &hg->VertexWeightDim, &hg->vwgt,
    &hg->EdgeWeightDim, &hg->ewgt, base);
   if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      fclose(f);
      goto End;
      }

   if (*base > 0) {
      /* Convert to zero-based vertex numbers */
      int i;
      for (i = 0; i < hg->nInput; i++)
         hg->hvertex[i] -= *base;
      }

   err = Zoltan_HG_Create_Mirror (zz, hg);
   if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      fclose(f);
      goto End;
      }

   if (fclose(f))
      err = ZOLTAN_WARN;

End:
   ZOLTAN_TRACE_EXIT(zz, yo);
   return err;
   }


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
