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

#include "hypergraph.h"



/*****************************************************************************/



void Zoltan_HG_HGraph_Init(
  HGraph *hgraph
)
{
  hgraph->info  = 0;
  hgraph->nVtx  = 0;
  hgraph->nEdge = 0;
  hgraph->nInput= 0;
  hgraph->nDim  = 0;
  hgraph->EdgeWeightDim   = 0;
  hgraph->VertexWeightDim = 0;

  hgraph->coor    = NULL;
  hgraph->vwgt    = NULL;
  hgraph->ewgt    = NULL;
  hgraph->hindex  = NULL;
  hgraph->hvertex = NULL;
  hgraph->vindex  = NULL;
  hgraph->vedge   = NULL;
}

/*****************************************************************************/



void Zoltan_HG_Graph_Init(
  Graph *graph
)
{
  graph->info  = 0;
  graph->nVtx  = 0;
  graph->nEdge = 0;
  graph->nDim  = 0;
  graph->EdgeWeightDim   = 0;
  graph->VertexWeightDim = 0;

  graph->vtxdist = NULL;

  graph->coor    = NULL;
  graph->vwgt    = NULL;
  graph->ewgt    = NULL;
  graph->nindex  = NULL;
  graph->neigh   = NULL;
}

/****************************************************************************/



int Zoltan_HG_HGraph_Free(
  HGraph *hg
)
{
/* Frees all memory associated with a hypergraph;
 * does not free the hypergraph itself.
 */
  if (hg) {
     ZOLTAN_FREE ((void**) &hg->coor);
     ZOLTAN_FREE ((void**) &hg->vwgt);
     ZOLTAN_FREE ((void**) &hg->ewgt);
     ZOLTAN_FREE ((void**) &hg->hindex);
     ZOLTAN_FREE ((void**) &hg->hvertex);
     ZOLTAN_FREE ((void**) &hg->vindex);
     ZOLTAN_FREE ((void**) &hg->vedge);
     }
  return ZOLTAN_OK;
}

/****************************************************************************/



/* Frees all memory associated with a graph; does not free the graph itself. */
int Zoltan_HG_Graph_Free(
  Graph *g
)
{
  if (g) {
     ZOLTAN_FREE ((void**) &g->coor);
     ZOLTAN_FREE ((void**) &g->vwgt);
     ZOLTAN_FREE ((void**) &g->ewgt);
     ZOLTAN_FREE ((void**) &g->vtxdist);
     ZOLTAN_FREE ((void**) &g->nindex);
     ZOLTAN_FREE ((void**) &g->neigh);
     }
  return ZOLTAN_OK;
}

/****************************************************************************/



int Zoltan_HG_Info (
  ZZ *zz,
  HGraph *hg
)
{
  int i, size, size_min, size_max;
  float wgt_min, wgt_max, wgt_tot;

  if (zz->Debug_Level < ZOLTAN_DEBUG_LIST)
     return ZOLTAN_OK;

  printf("--------- HGraph Information (min/ave/max/tot) ------------------\n");
  printf("info:%d |V|=%d |E|=%d |P|=%d \n", hg->info, hg->nVtx, hg->nEdge,
   hg->nInput);

  /* print weights */
  if (hg->vwgt) {
     wgt_tot = 0.0;
     wgt_min = FLT_MAX;
     wgt_max = FLT_MIN;
     for (i = 0; i < hg->nVtx; i++) {
        wgt_tot += hg->vwgt[i];
        wgt_min = MIN(wgt_min, hg->vwgt[i]);
        wgt_max = MAX(wgt_max, hg->vwgt[i]);
        }
     printf("Vertex weights   :    %9.2f %9.2f %9.2f %12.2f\n", wgt_min,
      wgt_tot/hg->nVtx, wgt_max, wgt_tot);
     }
  if (hg->ewgt) {
     wgt_tot = 0.0;
     wgt_min = FLT_MAX;
     wgt_max = FLT_MIN;
     for (i = 0; i < hg->nEdge; i++) {
        wgt_tot += hg->ewgt[i];
        wgt_min = MIN(wgt_min, hg->ewgt[i]);
        wgt_max = MAX(wgt_max, hg->ewgt[i]);
        }
     printf("HEdge weights    :    %9.2f %9.2f %9.2f %12.2f\n", wgt_min,
      wgt_tot/hg->nEdge, wgt_max, wgt_tot);
     }

  /* print sizes */
  if (hg->hindex) {
     size_min = INT_MAX;
     size_max = INT_MIN;
     for (i = 0; i < hg->nEdge; i++) {
        size     = hg->hindex[i+1] - hg->hindex[i];
        size_min = MIN(size_min, size);
        size_max = MAX(size_max, size);
        }
     printf("Edge sizes       :    %6d    %9.2f %6d    %9d\n", size_min,
      (float)(hg->nInput)/hg->nEdge,size_max,hg->nInput);
     }
  if (hg->vindex) {
     size_min = INT_MAX;
     size_max = INT_MIN;
     for (i = 0; i < hg->nVtx; i++) {
        size     = hg->vindex[i+1] - hg->vindex[i];
        size_min = MIN(size_min,size);
        size_max = MAX(size_max,size);
        }
     printf("Vertex sizes     :    %6d    %9.2f %6d    %9d\n",size_min,
      (float)(hg->nInput)/hg->nEdge,size_max,hg->nInput);
     }

  printf("-----------------------------------------------------------------\n");
  return ZOLTAN_OK;
}



/****************************************************************************/
/**  Given arrays to lookup a vertex for a given hyperedge (hindex, hvertex)
 **  or arrays to lookup a hyperedge given a vertex (vindex, vedge) create the
 **  other set of arrays (mirror).
 **  The mirror (to be filled in) should have NULL pointers in hg on entry!
 **  The calling program is responsible for freeing the mirror's memory  ***/

int Zoltan_HG_Create_Mirror (
  ZZ *zz,
  HGraph *hg
)
{
   int i, j;                  /* loop counters */
   int inlength, outlength;   /* input/output array lengths */
   int *index, *data;         /* pointers to input information */
   int *outindex, *outdata;
   char *yo = "Zoltan_HG_Create_Mirror";

   ZOLTAN_TRACE_ENTER(zz, yo);

   /* determine which data to "mirror" and set corresponding data pointers. */
   if (hg &&  hg->hindex && (hg->nInput == 0 || hg->hvertex)
    && !(hg->vindex) && !(hg->vedge)) {
      ZOLTAN_TRACE_DETAIL(zz, yo, "Have hindex; building vindex.");

      inlength  = hg->nEdge;
      outlength = hg->nVtx;
      index     = hg->hindex;
      data      = hg->hvertex;
      outindex  = hg->vindex = (int*)ZOLTAN_MALLOC((hg->nVtx+1) * sizeof(int));
      outdata   = hg->vedge  = (int*)ZOLTAN_MALLOC (hg->nInput  * sizeof(int));
      if (outindex == NULL || (hg->nInput > 0 && outdata == NULL)) {
         ZOLTAN_FREE ((void **) &(hg->vindex));
         ZOLTAN_FREE ((void **) &(hg->vedge));
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         ZOLTAN_TRACE_EXIT(zz, yo);
         return ZOLTAN_MEMERR;
         }
      }
   else if (hg && hg->vindex && (hg->nInput == 0 || hg->vedge)
    && !(hg->hindex) && !(hg->hvertex)) {
      ZOLTAN_TRACE_DETAIL(zz, yo, "Have vindex; building hindex.");

      inlength  = hg->nVtx;
      outlength = hg->nEdge;
      index     = hg->vindex;
      data      = hg->vedge;
      outindex  = hg->hindex  = (int*)ZOLTAN_MALLOC((hg->nEdge+1)* sizeof(int));
      outdata   = hg->hvertex = (int*)ZOLTAN_MALLOC(hg->nInput    *sizeof(int));
      if (outindex == NULL || (hg->nInput > 0 && outdata == NULL)) {
         ZOLTAN_FREE ((void **) &(hg->hindex));
         ZOLTAN_FREE ((void **) &(hg->hvertex));
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         ZOLTAN_TRACE_EXIT(zz, yo);
         return ZOLTAN_MEMERR;
         }
      }
   else {
      ZOLTAN_TRACE_EXIT(zz, yo);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input error.");
      return ZOLTAN_FATAL;  /* unable to proceed */
      }

   /* count number of data objects in out index space */
   for (i = 0; i < outlength+1; i++)
      outindex[i] = 0;
   for (i = 0; i < inlength; i++)
      for (j = index[i];  j < index[i+1]; j++)
          (outindex[data[j]+1])++;

   /* compute partial sum */
   for (i = 1; i < outlength; i++)
      outindex [i] += outindex[i-1];

   /* store data in out index location and increment index */
   for (i = 0; i < inlength; i++)
      for (j = index[i]; j < index[i+1]; j++)
         outdata [(outindex[data[j]])++ ] = i;

   /* out index now points past end of sublist. Shift array one position up */
   for (i = outlength; i > 0; i--)
      outindex[i] = outindex[i-1];
   outindex[0] = 0;

   ZOLTAN_TRACE_EXIT(zz, yo);
   return ZOLTAN_OK;
}



/****************************************************************************/
/* check that (hindex, hvertex) and (vindex, vedge) are consistant mappings */
/* returns ZOLTAN_WARN if not consistant, returns ZOLTAN_OK otherwise */
int Zoltan_HG_Check (
  ZZ *zz,
  HGraph *hg
)
{
  int i;
  int iedge;               /* iedge denotes an hyperedge */
  int j;                   /* j is the index of a vertex */
  int k;                   /* k is an index of hyperedge */
  int *check;
  int ierr = ZOLTAN_OK;
  char str[256];
  char *yo = "Zoltan_HG_Check";

  if (!(hg->hindex) || !(hg->hvertex) || !(hg->vindex) || !(hg->vedge) )
     return ZOLTAN_WARN;

  check = (int*) ZOLTAN_MALLOC(((hg->nEdge > hg->nVtx) ? hg->nEdge : hg->nVtx)
   * sizeof(int));
  if (check == NULL) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Unable to allocate memory.");
     return ZOLTAN_MEMERR;
     }

  for (i = 0; i < hg->nEdge; i++)
     check[i] = -1;
  for (i = 0; i < hg->nVtx; i++)
     for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)
        if (check[hg->vedge[j]] < i)
           check[hg->vedge[j]] = i;
        else {
           ZOLTAN_PRINT_WARN(zz->Proc, yo,
            "Found multiple hedge for same vertex.");
           ierr = ZOLTAN_WARN;
           }
  for (i = 0; i < hg->nVtx; i++)
     check[i] = -1;
  for (i = 0; i < hg->nEdge; i++)
     for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
        if (check[hg->hvertex[j]] < i)
           check[hg->hvertex[j]] = i;
        else {
           ZOLTAN_PRINT_WARN(zz->Proc, yo,
            "Found multiple vertex for same hedge.");
           ierr =  ZOLTAN_WARN;
           }
  ZOLTAN_FREE ((void**) &check);

  for (i = 0; i < hg->VertexWeightDim * hg->nVtx; i += hg->VertexWeightDim)
     if (hg->vwgt[i] < 0.0) {
        ZOLTAN_PRINT_WARN(zz->Proc, yo, "Found negative vertex weight.");
        ierr = ZOLTAN_WARN;
        }

  for (i = 0; i < hg->EdgeWeightDim * hg->nEdge; i += hg->EdgeWeightDim)
     if (hg->ewgt[i] < 0.0) {
        ZOLTAN_PRINT_WARN(zz->Proc, yo, "Found negative edge weight.");
        ierr = ZOLTAN_WARN;
        }

  for (i = 0; i < hg->nEdge; i++)
     if ((hg->hindex[i+1] - hg->hindex[i]) < 2) {
        sprintf (str, "Found hedge with less than two vertices: "
         "edge %d has %d vtxs\n", i, (hg->hindex[i+1] - hg->hindex[i]));
        ZOLTAN_PRINT_WARN(zz->Proc, yo, str);
        ierr = ZOLTAN_WARN;
        }

  for (i = 0; i < hg->nEdge; i++)
     for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
        if (hg->hvertex[j] < 0  ||  hg->hvertex[j] >= hg->nVtx) {
           ZOLTAN_PRINT_WARN(zz->Proc, yo,
            "Found vertex out of range in hvertex.");
           ierr = ZOLTAN_WARN;
           }

  for (i = 0; i < hg->nVtx; i++)
     for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)
        if (hg->vedge[j] < 0  ||  hg->vedge[j] >= hg->nEdge) {
           ZOLTAN_PRINT_WARN(zz->Proc, yo, "Found edge out of range in vedge.");
           ierr = ZOLTAN_WARN;
           }

  /* starting from (hindex,hvertex), for each edge determine each associated
   * vertex. Then for each vertex lookup associated edges using (vindex, vedge) */
  for (iedge = 0; iedge < hg->nEdge; iedge++)
     for (j = hg->hindex[iedge]; j < hg->hindex[iedge+1]; j++)
        /* for each hyperedge get index to vertices */
        {
        for (k=hg->vindex[hg->hvertex[j]]; k<hg->vindex[hg->hvertex[j]+1]; k++)
           /* for each vertex of current hyperedge get index to hyperedges */
           if (hg->vedge[k] == iedge)     /* does it match with original edge? */
              break;
        if (k == hg->vindex[hg->hvertex[j]+1])   /* if no match was found then */
           return ZOLTAN_WARN;                   /* failure, else keep on */
        }

  return ierr;
}

/****************************************************************************/



int Zoltan_HG_HGraph_to_Graph(
  ZZ *zz,
  HGraph *hg,
  Graph *g
)
{
  int 	i, j, k, e, roughly_e, pins, *_neigh, *degrees, current_degree,
        vertex1, vertex2, start, end, empty;
  float	*w, *_ewgt, weight;
  char  *yo = "Zoltan_HG_HGraph_to_Graph";

  Zoltan_HG_Graph_Init(g);
  g->info = hg->info;
  g->nVtx = hg->nVtx;
  g->nDim = hg->nDim;
  g->VertexWeightDim = hg->VertexWeightDim;
  g->EdgeWeightDim = hg->EdgeWeightDim;
  g->vtxdist = NULL;

   /* copy coordinates */
   /* not necessary at the moment, because matching does not use them */

  /* copy vertex weights */
  if (hg->vwgt) {
     if (!(g->vwgt = (float*) ZOLTAN_MALLOC (g->nVtx * sizeof(float)))) {
        ZOLTAN_FREE ((void**) &g->vwgt);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        return ZOLTAN_MEMERR;
        }
     for (i=0; i<g->nVtx; i++)
        g->vwgt[i] = hg->vwgt[i];
     }

  /* calculate the roughly |E| and each degree */
  roughly_e = 0;
  if (!(degrees = (int*) ZOLTAN_CALLOC (hg->nVtx, sizeof(int))))
     {
     ZOLTAN_FREE ((void**) &degrees);
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i = 0; i < hg->nEdge; i++) {
     /* if (hg->hindex[i+1] - hg->hindex[i] <= 10) */
     pins = hg->hindex[i+1] - hg->hindex[i];
     roughly_e += pins*(pins-1);
     pins--;
     for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
        degrees[hg->hvertex[j]] += pins;
     }

  /* Calculate the initial nindex */
  if (!(g->nindex = (int*) ZOLTAN_MALLOC ((hg->nVtx + 1) * sizeof(int)))) {
     ZOLTAN_FREE ((void**) &g->nindex);
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  e = 0;
  for (i = 0; i < hg->nVtx; i++) {
     current_degree = degrees[i];
     degrees[i]     = g->nindex[i] = e;
     e += current_degree;
     }
  g->nindex[hg->nVtx] = e;

  /* Calculate the initial neigh and ewgt */
  if (!(_neigh = (int*)   ZOLTAN_MALLOC (roughly_e * sizeof(int)))
   || !(_ewgt  = (float*) ZOLTAN_MALLOC (roughly_e * sizeof(float))) ) {
      ZOLTAN_FREE ((void**) &_neigh);
      ZOLTAN_FREE ((void**) &_ewgt);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
  for (i = 0; i < hg->nEdge; i++) {
     /* if (hg->hindex[i+1]-hg->hindex[i] <= 10) */
     pins = hg->hindex[i+1] - hg->hindex[i];
     /* weight = 1.0 / (pins-1); */
     weight = 2.0 / ((pins-1) * pins);
     if (hg->ewgt)
        weight *= hg->ewgt[i];
     for (j  = hg->hindex[i]; j < hg->hindex[i+1]; j++) {
        vertex1 = hg->hvertex[j];
        for (k = j + 1; k < hg->hindex[i+1]; k++) {
           vertex2                   = hg->hvertex[k];
           _neigh[degrees[vertex1]]  = vertex2;
           _neigh[degrees[vertex2]]  = vertex1;
           _ewgt[degrees[vertex1]++] = _ewgt[degrees[vertex2]++] = weight;
           }
        }
     }
  ZOLTAN_FREE ((void **) &degrees);

  /* Compact identical incident edges and their weight */
  if (!(w = (float*) ZOLTAN_CALLOC (hg->nVtx, sizeof(float)))) {
     ZOLTAN_FREE ((void**) &w);
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  empty = 0;
  for (i = 0; i < hg->nVtx; i++)
     {
     start = g->nindex[i];
     end   = g->nindex[i+1];
     g->nindex[i] -= empty;
     for (j = start; j < end; j++)
        if (_ewgt[j] > 0.0)
           w[_neigh[j]] += _ewgt[j];
     for (j = start; j < end; j++)
        if (_ewgt[j] > 0.0  &&  w[_neigh[j]] > 0.0)
           {
           _neigh[j-empty] = _neigh[j];
           _ewgt[j-empty]  = w[_neigh[j]];
           w[_neigh[j]]    = 0.0;
           }
        else
           empty++;
  }
  g->nindex[hg->nVtx] -= empty;
  ZOLTAN_FREE ((void**) &w);

/* Reallocate to the exact size */
  g->nEdge = g->nindex[g->nVtx];

  if (!(g->neigh = (int*) ZOLTAN_MALLOC (g->nindex[hg->nVtx] * sizeof(int)))) {
     ZOLTAN_FREE ((void**) &g->neigh);
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  memcpy(g->neigh, _neigh, g->nEdge * sizeof(int));
  ZOLTAN_FREE ((void**) &_neigh);

  if (!(g->ewgt = (float*)ZOLTAN_MALLOC(g->nindex[hg->nVtx] * sizeof(float)))) {
     ZOLTAN_FREE ((void**) &g->ewgt);
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  memcpy(g->ewgt,_ewgt,g->nEdge*sizeof(float));
  ZOLTAN_FREE ((void **) &_ewgt);

  return ZOLTAN_OK;
}

/****************************************************************************/



int Zoltan_HG_Graph_to_HGraph(
  ZZ *zz,
  Graph *g,        /* Input graph */
  HGraph *hg       /* Ouput hypergraph */
)
{
/*
 *  Converts a graph g into a hypergraph hg.
 *  One hyperedge is created for each vertex of g.
 *  Hyperedge i consists of vertex i + all edge neighbors of i in graph g.
 */
char *yo = "Zoltan_HG_Graph_to_HGraph";
int i, j;
int *hindex = NULL, *hvertex = NULL;  /* temporary array pointers */
int cnt;
int ierr = ZOLTAN_OK;

  ZOLTAN_TRACE_ENTER(zz, yo);

  Zoltan_HG_HGraph_Init(hg);
  hg->info  = g->info;
  hg->nVtx  = g->nVtx;
  hg->nEdge = g->nVtx;
  hg->nInput  = g->nVtx + g->nEdge;
  hg->nDim  = g->nDim;
  hg->VertexWeightDim = g->VertexWeightDim;
  hg->EdgeWeightDim = g->EdgeWeightDim;

  /* Copy coordinates from graph to hypergraph */
  if (g->coor) {
     ZOLTAN_TRACE_DETAIL(zz, yo, "Copying coordinates");
     cnt = hg->nVtx * g->nDim;
     if (!(hg->coor = (double *) ZOLTAN_MALLOC(cnt * sizeof(double)))) {
        ierr = ZOLTAN_MEMERR;
        goto End;
        }
    memcpy(hg->coor, g->coor, cnt * sizeof(double));
    }

  /* Copy vertex weights from graph to hypergraph */
  if (g->vwgt) {
     ZOLTAN_TRACE_DETAIL(zz, yo, "Copying vertex weights");
     cnt = hg->nVtx * hg->VertexWeightDim;
     if (!(hg->vwgt = (float *) ZOLTAN_MALLOC(cnt * sizeof(float)))) {
        ierr = ZOLTAN_MEMERR;
        goto End;
        }
     memcpy(hg->vwgt, g->vwgt, cnt * sizeof(float));
     }

  /* Allocate memory for hypergraph edges */
  ZOLTAN_TRACE_DETAIL(zz, yo, "Allocating edge arrays");
  if (hg->nEdge > 0 && hg->nInput > 0)
    if (!(hindex  = hg->hindex =(int*)ZOLTAN_MALLOC((hg->nEdge+1)*sizeof(int)))
     || !(hvertex = hg->hvertex=(int*)ZOLTAN_MALLOC(hg->nInput*sizeof(int))) ) {
        ierr = ZOLTAN_MEMERR;
        goto End;
        }

  /* at the moment, the hyperedge weight will be constant !
  if (g->ewgt) {
     ZOLTAN_TRACE_DETAIL(zz, yo, "Allocating edge weights");
     cnt = hg->nEdge * hg->EdgeWeightDim;
     if (!(hg->ewgt = (float *) ZOLTAN_MALLOC(cnt * sizeof(float)))) {
        ierr = ZOLTAN_MEMERR;
        goto End;
        }
     }
   */

  /* KDD  The following section will have to be reworked for parallel HGs. */
  ZOLTAN_TRACE_DETAIL(zz, yo, "Filling edge arrays");
  cnt = 0;
  for (i = 0; i < hg->nEdge; i++) {
     hindex[i] = cnt;
     hvertex[cnt++] = i;
     for (j = g->nindex[i]; j < g->nindex[i+1]; j++) {
        hvertex[cnt++] = g->neigh[j];
        /* at the moment, the hyperedge weight will be constant
        for (k = 0; k < hg->EdgeWeightDim; k++)
           hg->ewgt[i*(hg->EdgeWeightDim)+k] += g->ewgt[j*(hg->EdgeWeightDim)+k];
        */
        }
     }
  hindex[hg->nEdge] = cnt;

  /* Sanity check */
  if (hg->nInput != cnt) {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Sanity check failed: nInput != cnt.");
     ierr = ZOLTAN_FATAL;
     goto End;
     }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Creating mirror");
  ierr = Zoltan_HG_Create_Mirror(zz, hg);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error in building mirror.");


End:
  if (ierr == ZOLTAN_MEMERR)
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory.");
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
     ZOLTAN_FREE ((void**) &(hg->coor));
     ZOLTAN_FREE ((void**) &(hg->vwgt));
     ZOLTAN_FREE ((void**) &(hg->ewgt));
     ZOLTAN_FREE ((void**) &(hg->hindex));
     ZOLTAN_FREE ((void**) &(hg->hvertex));
  }
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/****************************************************************************/



/* Routine to print hypergraph weights and edges. Assumes serial execution;
 * put inside Zoltan_Print_Sync_Start/Zoltan_Print_Sync_End for parallel
 * programs. */
void Zoltan_HG_Print(
  ZZ *zz,
  HGraph *hg
)
{

int i, j;
int num_vwgt;
int num_ewgt;

  if (hg == NULL)
     return;

  num_vwgt = hg->VertexWeightDim;
  num_ewgt = hg->EdgeWeightDim;
  if (hg->vwgt != NULL) {
     printf("Vertices: [weights])\n");
     for (i = 0; i < hg->nVtx; i++) {
        printf("%d:  [", i);
        for (j = 0; j < num_vwgt; j++)
           printf("%f ", hg->vwgt[i*num_vwgt + j]);
        printf("])\n");
        }
     }

  /* Print Hyperedge Info */
  printf("Hyperedges:  (vertices)\n");
  for (i = 0; i < hg->nEdge; i++) {
     printf("%d:  ", i);
     printf("(");
     for (j = hg->hindex[i]; j < hg->hindex[i+1]; j++)
       printf("%d ", hg->hvertex[j]);
     printf(")\n");
     }

  if (hg->ewgt != NULL) {
     printf("Hyperedge Weights:  [weights]\n");
     for (i = 0; i < hg->nEdge; i++) {
        printf("%d:  ", i);
        printf("[");
        for (j = 0; j < num_ewgt; j++)
           printf("%f ", hg->ewgt[i*num_ewgt + j]);
        printf("])\n");
        }
     }
}

/****************************************************************************/



/*  Random Number generator due to Knuth found in Numerical Recipes in C
   (2nd edition) by Press, Vetterling, Teukolsky, Flannery (Page 284.)
   Currently, for use in the hypergraph library only. Needed because different
   random number implimentations on different machines produced different
   answers! Hence this is a single portable, fast, algorithm with adequate
   random number generation. */


static unsigned long idum = 123456789;

unsigned long Zoltan_HG_Rand (void) {
   return idum = 1664525L * idum + 1013904223L;
   }

void Zoltan_HG_Srand (unsigned long seed) {
   idum = seed;
   }



/* Randomly permute an array of ints. */
void Zoltan_HG_Rand_Perm_Int (int *data, int n)
{
  int i, number, temp;

  for (i = n; i > 0; i--) {
     number       = Zoltan_HG_Rand() % i;
     temp         = data[number];
     data[number] = data[i-1];
     data[i-1]    = temp;
     }
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
