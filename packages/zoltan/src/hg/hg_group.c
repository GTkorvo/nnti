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

static ZOLTAN_HG_GROUPING_FN grouping_mxg;  /* maximal grouping */
static ZOLTAN_HG_GROUPING_FN grouping_reg;  /* random edge grouping */
static ZOLTAN_HG_GROUPING_FN grouping_heg;  /* heavy edge grouping */
static ZOLTAN_HG_GROUPING_FN grouping_grg;  /* greedy grouping */
static ZOLTAN_HG_GROUPING_FN grouping_deg;  /* decreasing grouping */
static ZOLTAN_HG_GROUPING_FN grouping_rrg;  /* random vertex, random edge grouping */


/****************************************************************************/

ZOLTAN_HG_GROUPING_FN *Zoltan_HG_Set_Grouping_Fn(char *str)
{
  static int srand_set ;
  if (srand_set == 0)
     {
     srand_set = 1 ;
     srand ((unsigned long) RANDOM_SEED) ;
     }

  if      (strcasecmp(str, "mxg") == 0)  return grouping_mxg;
  else if (strcasecmp(str, "reg") == 0)  return grouping_reg;
  else if (strcasecmp(str, "heg") == 0)  return grouping_heg;
  else if (strcasecmp(str, "grg") == 0)  return grouping_grg;
  else if (strcasecmp(str, "deg") == 0)  return grouping_deg;
  else if (strcasecmp(str, "rrg") == 0)  return grouping_rrg;
  else                                   return NULL;
}


/****************************************************************************/

int Zoltan_HG_Grouping (ZZ *zz, HGraph *hg, Packing pack, HGParams *hgp)
{ int   i, j;
  int   limit=0 ;   /* reserved for future use */
  int   ierr = ZOLTAN_OK;
  float *old_ewgt=NULL, weight, sum1, sum2;
  char  *yo = "Zoltan_HG_Grouping";

  old_ewgt = hg->ewgt;
  hg->ewgt = (float *) ZOLTAN_MALLOC (sizeof (float) * hg->nEdge);
  if (hg->ewgt == NULL)
     {
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }

  for (i=0; i<hg->nEdge; i++)
  { sum1 = sum2 = 0.0;
    if (hg->vwgt)
    { for (j=hg->hindex[i]; j<hg->hindex[i+1]; j++)
      { weight = hg->vwgt[hg->hvertex[j]];
        sum1 += weight;
        sum2 += weight*weight;
    } }
    else
      sum1 = sum2 = (float)(hg->hindex[i+1]-hg->hindex[i]);
    sum1 = (sum1*sum1-sum2)/2.0;
    if (sum1 == 0.0)
      hg->ewgt[i] = FLT_MAX;
    else
      hg->ewgt[i] = (old_ewgt?old_ewgt[i]:1.0)/sum1;
  }

  ierr = hgp->grouping(zz,hg,pack,limit);
  ZOLTAN_FREE ((void **) &hg->ewgt);
  hg->ewgt = old_ewgt;
  return ierr;
}


/****************************************************************************/

static int grouping_mxg (ZZ *zz, HGraph *hg, Packing pack, int limit)
   {                                    /* limit is defined for future use */
   int i, j, vertex, first_vertex ;
   char *yo = "grouping_mxg" ;

   for (i = 0 ; i < hg->nVtx ; i++)
      pack[i] = i ;

   for (i = 0 ; i < hg->nEdge ; i++)
      for (j = hg->hindex[i] ; j < hg->hindex[i+1] ; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[i+1] ; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j])
                  vertex = pack[vertex] = hg->hvertex[j] ;
            pack[vertex] = first_vertex ;
            break ;       /* not required, might improve speed */
            }
   return ZOLTAN_OK ;
   }


/****************************************************************************/

static int grouping_reg (ZZ *zz, HGraph *hg, Packing pack, int limit)
   {                               /* limit is defined for future use */
   int i, j, *edges=NULL, edge, random, vertex, first_vertex ;
   char *yo = "grouping_reg" ;

   edges = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge) ;
   if (edges == NULL)
      {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
   for (i = 0 ; i < hg->nEdge ; i++)
      edges[i] = i ;
   for (i = 0 ; i < hg->nVtx ;  i++)
      pack[i]  = i ;

   for (i = hg->nEdge ; i > 0 ; i--)
      {
      random = rand() % i ;
      edge = edges[random] ;
      edges[random] = edges[i-1] ;

      for (j = hg->hindex[edge] ; j < hg->hindex[edge+1] ; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[edge+1] ; j++)
                if (pack[hg->hvertex[j]] == hg->hvertex[j])
                   vertex = pack[vertex] = hg->hvertex[j] ;
            pack[vertex] = first_vertex ;
            break ;
            }
      }
   ZOLTAN_FREE ((void **) &edges) ;
   return ZOLTAN_OK ;
   }


/****************************************************************************/

static int grouping_rrg (ZZ *zz, HGraph *hg, Packing pack, int limit)
   {                                     /* limit is defined for future use */
   int i, j, edge, random, *vertices=NULL, vertex, first_vertex, count ;
   char *yo = "grouping_rrg" ;

   vertices  = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx) ;
   if (vertices == NULL)
      {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
   for (i = 0 ; i < hg->nVtx ;  i++)
      vertices[i] = pack[i] = i ;

   for (i = hg->nVtx ; i > 0 ; i--)
      {
      random = rand() % i ;
      vertex = vertices[random] ;
      vertices[random] = vertices[i-1] ;
      if (pack[vertex] != vertex)
         continue ;          /* vertex already packed, move on */

      count = hg->vindex[vertex+1] - hg->vindex[vertex] ;  /* count edges */
      if (count == 0)
         continue ;

      edge = hg->vedge[hg->vindex[vertex] + (rand() % count)] ;  /* random edge */
      for (j = hg->hindex[edge] ; j < hg->hindex[edge+1] ; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[edge+1] ; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j])
                  vertex = pack[vertex] = hg->hvertex[j] ;
            pack[vertex] = first_vertex ;
            break ;
            }
      }
   ZOLTAN_FREE ((void **) &vertices) ;
   return ZOLTAN_OK ;
   }


/****************************************************************************/

static int grouping_heg (ZZ *zz, HGraph *hg, Packing pack, int limit)
{                                /* limit is defined for future use */
   int   i, j, *vertices=NULL, *del_edges=NULL, vertex, first_vertex, edge,
         number, best_edge, best_size;
   float best_ewgt;
   char  *yo = "grouping_heg" ;

   vertices  = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx) ;
   del_edges = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge) ;
   if (vertices == NULL  ||  del_edges == NULL)
      {
      ZOLTAN_FREE ((void **) &vertices) ;
      ZOLTAN_FREE ((void **) &del_edges) ;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
      }
   for (i = 0 ; i < hg->nVtx ; i++)
      pack[i] = vertices[i] = i ;
   for (i = 0 ; i < hg->nEdge ; i++)
      del_edges[i] = 0;

   for (i = hg->nVtx ; i > 0 ; i--)
      {
      number = rand() % i ;
      vertex = vertices[number] ;
      vertices[number] = vertices[i-1] ;
      if (pack[vertex] != vertex)
         continue ;            /* vertex is already matched, move on */

      best_edge = best_size = -1;
      best_ewgt = -1.0 ;
      for (j = hg->vindex[vertex] ; j < hg->vindex[vertex+1] ; j++)
         {
         int size ;
         edge = hg->vedge[j];
         size = hg->hindex[edge+1] - hg->hindex[edge] ;

         if (del_edges[edge]==0 && ((hg->ewgt == NULL  && size < best_size)
          || (hg->ewgt && (hg->ewgt[edge] >  best_ewgt
                       || (hg->ewgt[edge] == best_ewgt && size < best_size)))))
              {
              best_edge = edge ;
              best_ewgt = hg->ewgt[edge] ;
              best_size = hg->hindex[edge+1]-hg->hindex[edge];
              }
         }
      if (best_edge == -1)
         continue ;                       /* no suitable edge found */

      for (j = hg->hindex[best_edge] ; j < hg->hindex[best_edge+1] ; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[best_edge+1] ; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j])
                  vertex = pack[vertex] = hg->hvertex[j] ;
            pack[vertex] = first_vertex ;
            break ;
            }
      del_edges[best_edge] = 1 ;
      }
   ZOLTAN_FREE ((void **) &vertices) ;
   ZOLTAN_FREE ((void **) &del_edges) ;
   return ZOLTAN_OK ;
}


/****************************************************************************/

static void quickpart_dec_float_int (float *val1, int *val2, int *sorted,
int start, int end, int *equal, int *smaller)
{ int	i, next, key2, key2_next;
  float key1, key1_next;

  i = (end+start)/2;
  key1 = val1?val1[sorted[i]]:1.0;
  key2 = val2?val2[sorted[i]]:1;

  (*equal) = (*smaller) = start;
  for (i=start; i<=end; i++)
  { next = sorted[i];
    key1_next = val1?val1[next]:1.0;
    key2_next = val2?val2[next]:1;
    if (key1_next>key1 || (key1_next==key1 && key2_next>key2))
    { sorted[i] = sorted[(*smaller)];
      sorted[(*smaller)++] = sorted[(*equal)];
      sorted[(*equal)++] = next;
    }
    else if (key1_next==key1 && key2_next==key2)
    { sorted[i] = sorted[(*smaller)];
      sorted[(*smaller)++] = next;
  } }
}


static void quicksort_dec_float_int (float* val1, int *val2, int* sorted,
int start, int end)
{ int  equal, smaller;

  if (start < end)
  { quickpart_dec_float_int (val1,val2,sorted,start,end,&equal,&smaller);
    quicksort_dec_float_int (val1,val2,sorted,start,equal-1);
    quicksort_dec_float_int (val1,val2,sorted,smaller,end);
  }
}


/****************************************************************************/

static int grouping_grg (ZZ *zz, HGraph *hg, Packing pack, int limit)
{                                      /* limit is defined for future use */
  int   i, j, *size=NULL, *sorted=NULL, first_vertex, vertex ;
  char *yo = "grouping_grg" ;

  for (i=0; i<hg->nVtx; i++)
    pack[i] = i;

/* Sort the hyperedges according to their weight and size */
  size   = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge);
  sorted = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge);
  if (size == NULL || sorted == NULL)
     {
     ZOLTAN_FREE ((void **) &size) ;
     ZOLTAN_FREE ((void **) &sorted) ;
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i=0; i<hg->nEdge; i++)
     size[i] = -(hg->hindex[i+1]-hg->hindex[i]);
  for (i=0; i<hg->nEdge; i++)
     sorted[i] = i;
  quicksort_dec_float_int(hg->ewgt,size,sorted,0,hg->nEdge-1);
  ZOLTAN_FREE ((void **) &size);

  /* Match hyperedges along decreasing weight */
  for (i=0; i<hg->nEdge; i++)
      for (j = hg->hindex[sorted[i]] ; j < hg->hindex[sorted[i]+1] ; j++)
         if (pack[hg->hvertex[j]] == hg->hvertex[j])
            {
            first_vertex = vertex = hg->hvertex[j] ;
            for (j++ ; j < hg->hindex[sorted[i]+1] ; j++)
               if (pack[hg->hvertex[j]] == hg->hvertex[j])
                  vertex = pack[vertex] = hg->hvertex[j] ;
            pack[vertex] = first_vertex ;
            break ;
            }
  ZOLTAN_FREE ((void **) &sorted);
  return ZOLTAN_OK;
  }


/****************************************************************************/

static int grouping_deg (ZZ *zz, HGraph *hg, Packing pack, int limit)
{                                        /* limit is defined for future use */
  int   i, j, *size=NULL, *sorted=NULL, vertex, first_vertex;
  char *yo = "grouping_deg" ;

  for (i=0; i<hg->nVtx; i++)
    pack[i] = i;

/* Sort the hyperedges according to their weight and size */
  size   = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge);
  sorted = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nEdge);
  if (size == NULL || sorted == NULL)
     {
     ZOLTAN_FREE ((void **) &size) ;
     ZOLTAN_FREE ((void **) &sorted) ;
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     return ZOLTAN_MEMERR;
     }
  for (i=0; i<hg->nEdge; i++)
    size[i] = -(hg->hindex[i+1]-hg->hindex[i]);
  for (i=0; i<hg->nEdge; i++)
    sorted[i] = i;
  quicksort_dec_float_int(hg->ewgt,size,sorted,0,hg->nEdge-1);
  ZOLTAN_FREE ((void **) &size);

/* Pack hyperedges along decreasing weight */
  for (i=0; i<hg->nEdge; i++)
  { j = hg->hindex[sorted[i]];
    while (j<hg->hindex[sorted[i]+1] && pack[hg->hvertex[j]]!=hg->hvertex[j])
      j++;
    if (j<hg->hindex[sorted[i]+1])
    { first_vertex = vertex = hg->hvertex[j];
      j++;
      while (j<hg->hindex[sorted[i]+1])
      { if (pack[hg->hvertex[j]] == hg->hvertex[j])
        { pack[vertex] = hg->hvertex[j];
          vertex = hg->hvertex[j];
        }
        j++;
      }
      pack[vertex] = first_vertex;
  } }
  ZOLTAN_FREE ((void **) &sorted);
  return ZOLTAN_OK;
}


/****************************************************************************/

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
