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

#define BUF_LEN 1000000


int hg_readfile (ZZ *zz, HGraph *hg, char *hgraphfile, int *base)
{
FILE *f;
int err;
char errstr[200];
char *yo = "hg_readfile";

   Zoltan_HG_HGraph_Init(hg);

   f = fopen (hgraphfile, "r");
   if (!f) {
      sprintf(errstr, "ERROR...not able to open file %s!\n", hgraphfile);
      ZOLTAN_PRINT_ERROR (zz->Proc, yo, errstr);
      return ZOLTAN_FATAL;
      }

   err = Zoltan_HG_Readfile (0, f, &hg->nVtx, &hg->nEdge, &hg->nInput,
    &hg->hindex, &hg->hvertex, &hg->VertexWeightDim, &hg->vwgt,
    &hg->EdgeWeightDim, &hg->ewgt, base);
   if (err != ZOLTAN_OK && err != ZOLTAN_WARN) {
      fclose(f);
      return err;
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
      return err;
      }

   if (fclose(f))
      return ZOLTAN_WARN;
   return ZOLTAN_OK;
   }


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
