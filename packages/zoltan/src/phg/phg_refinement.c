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

#include "phg.h"

static ZOLTAN_PHG_REFINEMENT_FN refine_no;
static ZOLTAN_PHG_REFINEMENT_FN refine_fm2;
static ZOLTAN_PHG_REFINEMENT_FN refine_grkway;

/****************************************************************************/

ZOLTAN_PHG_REFINEMENT_FN *Zoltan_PHG_Set_Refinement_Fn(char *str)
{
  
  if      (!strcasecmp(str, "grkway"))         return refine_grkway;
  if      (!strcasecmp(str, "fm2"))            return refine_fm2;  
  else if (!strcasecmp(str, "no"))             return refine_no;
  else                                         return NULL;
}



/****************************************************************************/
int Zoltan_PHG_Refinement (ZZ *zz, PHGraph *hg, int p, Partition part,
PHGPartParams *hgp)
{
  return hgp->Refinement(zz, hg, p, part, hgp, hgp->bal_tol);
}



/****************************************************************************/
static int refine_no (
  ZZ *zz,     /* Zoltan data structure */
  PHGraph *hg,
  int p,
  Partition part,
  PHGPartParams *hgp,
  float bal_tol
)
{
  return ZOLTAN_OK;
}



/*****************************************************************************/
/* 2-way Parallel FM refinement                                              */
static char *uMe(PHGComm *hgc)
{
    static char msg[1024];

    sprintf(msg, "<%d/%d>: (%d,%d)/[%d,%d] ->", hgc->Proc, hgc->Num_Proc, hgc->myProc_x, hgc->myProc_y, hgc->nProc_x, hgc->nProc_y);
    return msg;
}

static int refine_fm2 (
  ZZ *zz,
  PHGraph *hg,
  int p,
  Partition part,
  PHGPartParams *hgp,
  float bal_tol
)
{
    int    i, j, vertex, edge, *pins[2], *tpins[2], *moves;
    double total_weight, weights[2], tweights[2], max_weight[2];
    double cutsize, best_cutsize, *gain = 0;
    HEAP   heap[2];
    char   *yo="local_fm2";
    double ratio = hg->ratio, error, best_error;
    int    best_imbalance, imbalance;
    PHGComm *hgc=&hgp->comm;

#ifdef _DEBUG    
    FILE *fp;
    char fname[512];
    int i;

    sprintf(fname, "uvc-phg-debug.%d.txt", zz->Proc);    
    fp = fopen(fname, "w");
    if (!fp) {
        printf ("failed to open debug file '%s'\n", fname);
        return ZOLTAN_FATAL;
    }
    fprintf(fp, "Process (%d, %d) of [%d, %d]\n", hg->myProc_x, hg->myProc_y, hg->nProc_x, hg->nProc_y);
    fprintf(fp, "H(%d, %d, %d)\n", hg->nVtx, hg->nEdge, hg->nNonZero);
    fprintf(fp, "p=%d  bal_tol = %.3f\n PartVec:\n", p, bal_tol);    
    for (i = 0; i < hg->nVtx; ++i)
        fprintf(fp, "%d ", part[i]);
    fprintf(fp, "\n\n");
    Zoltan_PHG_Print(zz, hg, fp);
    fclose(fp);
#endif

    if (p != 2) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "p!=2 not allowed for phg_fm2.");
        return ZOLTAN_FATAL;
    }
    
    if (hg->nEdge == 0)
        return ZOLTAN_OK;
    
    /* Calculate the weights in each partition and total, then maxima */
    weights[0] = weights[1] = 0.0;
    tweights[0] = tweights[1] = 0.0;
    if (hg->vwgt)  
        for (i = 0; i < hg->nVtx; i++)
            tweights[part[i]] += hg->vwgt[i];
    else  
        for (i = 0; i < hg->nVtx; i++)
            tweights[part[i]] += 1.0;

    MPI_Allreduce(tweights, weights, 2, MPI_DOUBLE, MPI_SUM, hgc->row_comm);
    total_weight = weights[0] + weights[1];
    max_weight[0] = total_weight * bal_tol *      ratio;
    max_weight[1] = total_weight * bal_tol * (1 - ratio);

    printf("%s W:(%.2lf, %.2lf) MaxW:(%.2lf, %.2lf)\n", uMe(hgc), weights[0], weights[1], max_weight[0], max_weight[1]);
    
    if (!(pins[0]    = (int*)   ZOLTAN_CALLOC(2 * hg->nEdge, sizeof(int)))
        || !(tpins[0] = (int*)   ZOLTAN_CALLOC(2 * hg->nEdge, sizeof(int)))
        || !(moves  = (int*)   ZOLTAN_CALLOC    (hg->nVtx,  sizeof(int)))
        || !(gain   = (double*)ZOLTAN_CALLOC    (hg->nVtx,  sizeof(double))) ) {
         Zoltan_Multifree(__FILE__,__LINE__, 4, &pins[0], &tpins[0], &moves, &gain);
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
    }
    pins[1] = &(pins[0][hg->nEdge]);
    tpins[1] = &(tpins[0][hg->nEdge]);

    /* Initial calculation of the pins distribution  */
    for (i = 0; i < hg->nEdge; ++i)
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j)
            ++(tpins[part[hg->hvertex[j]]][i]);

    MPI_Allreduce(tpins[0], pins[0], 2*hg->nEdge, MPI_INT, MPI_SUM, hgc->row_comm); 
    cutsize = 0.0;
    for (i=0; i < hg->nEdge; ++i) {
        if (pins[0][i] && pins[1][i])
            cutsize += (hg->ewgt ? hg->ewgt[i] : 1.0);
    }

    /* Just for debugging */
    best_cutsize = Zoltan_PHG_hcut_size_total(hgc, hg, part, p);
    error = Zoltan_PHG_hcut_size_links(hgc, hg, part, p);
    printf("%s: Initial cutsize=%.2lf Verify: total=%.2lf links=%.2lf\n", uMe(hgc), cutsize,
           best_cutsize,
           error);
    if (error!=cutsize || best_cutsize!=cutsize) {
        Zoltan_Multifree(__FILE__, __LINE__, 4, &pins[0], &tpins[0], &moves, &gain);
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid cut!!!");
    }
    /* debuggging code ends here */
        

    
    Zoltan_Multifree(__FILE__, __LINE__, 4, &pins[0], &tpins[0], &moves, &gain);
    return ZOLTAN_OK;
}




/*****************************************************************************/
/* This algorithm is loosely based on "A Coarse-Grained Parallel Formulation */
/* of Multilevel k-way Graph Partitioning Algorithm", Karypis & Kumar, 1997. */
/* It is implemented in serial as a testbed for future parallel development  */

typedef struct
   {
   double weight;
   float gain;
   int   vertex;
   int   source;
   int   destination;
   } Vdata;


static int comparison_gain_weight_no  (const void*, const void*);
static int comparison_weight (const void*, const void*);

static int refine_grkway (
  ZZ *zz,
  PHGraph *hg,
  int p,
  Partition part,
  PHGPartParams *hgp,
  float bal_tol
)
{
  const int MAX_LOOP = 7;
  int     i, j, loop, vertex, edge, ipart;  /* loop counters */
  double *part_weight, total_weight, max_weight;
  double *gain, tgain;
  int   **cuts, *store1, *listend, *movect;
  Vdata **lists, *store2;
  int     bestpart;
  char   *yo= "refine_grkway";

  double smallest;
  int found, smallpart;
  int **up, *upstore;

  /* Allocate necessary storage for heaps and weight calculation */
  if  (!(part_weight = (double*) ZOLTAN_CALLOC (p,         sizeof (double)))
   ||  !(gain        = (double*) ZOLTAN_CALLOC (p,         sizeof (double)))
   ||  !(cuts        = (int**)   ZOLTAN_CALLOC (hg->nEdge, sizeof (int)))
   ||  !(store1      = (int*)    ZOLTAN_CALLOC (hg->nEdge * p, sizeof (int)))
   ||  !(lists       = (Vdata**) ZOLTAN_CALLOC (p,         sizeof (Vdata)))
   ||  !(store2      = (Vdata*)  ZOLTAN_CALLOC (hg->nVtx * p, sizeof (Vdata)))
   ||  !(listend     = (int*)    ZOLTAN_CALLOC (p,         sizeof (int)))
   ||  !(movect      = (int*)    ZOLTAN_CALLOC (hg->nVtx,  sizeof (int))))  {
    Zoltan_Multifree(__FILE__,__LINE__, 8, &part_weight, &gain, &cuts,
      &store1, &lists, &store2, &listend, &movect);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
     
upstore = (int*)  ZOLTAN_CALLOC (p*p, sizeof (int));
up      = (int**) ZOLTAN_CALLOC (p,   sizeof (int));

for (i = 0; i < p; i++)
   up[i] = upstore + i * p;
   
for (i = 0; i < p; i++)
   for (j = 0; j < p; j++)
      if ((j > i && j <= i + p/2) || (j < i && j <= i - p/2))
         up[i] [j] = 1;

  /* simulate 2 dimensional arrays whose dimensions are known at run time */
  for (edge = 0; edge < hg->nEdge; edge++)
    cuts[edge]   = store1 + edge  * p;
  for (ipart = 0; ipart < p; ipart++)
    lists[ipart] = store2 + ipart * hg->nVtx;     
  for (i = 0; i < hg->nVtx; i++)
    movect [i] = -1;
      
  /* Calculate the total weights (local vertices weight) */
  /* most general case, need MPI global communication for total weight */
  total_weight = 0.0;
  memset (part_weight, 0, p * sizeof (double));
  if (hg->vwgt)
    for (vertex = 0; vertex < hg->nVtx; vertex++)  {
      total_weight              += hg->vwgt[vertex];
      part_weight[part[vertex]] += hg->vwgt[vertex];
    }
  else {
    total_weight = hg->nVtx;
    for (vertex = 0; vertex < hg->nVtx; vertex++)
      part_weight[part[vertex]] += 1.0;
  }
  max_weight = bal_tol * total_weight / p;
  
  /* determine if there are any overfilled partitions */
  smallest = part_weight[0];
  smallpart = 0;
  found = 0;
  memset (listend, 0, p * sizeof (int));
  for (ipart = 0; ipart < p; ipart++)  {
    if (part_weight[ipart] < smallest)  {
      smallest = part_weight[ipart];
      smallpart = ipart;
    }
    if (part_weight[ipart] > max_weight)
      found = 1;
  }

  /* overfilled partitions found, move vertices out of it */     
  if (found)  {
    for (vertex = 0; vertex < hg->nVtx; vertex++)  {
      ipart = part[vertex];

      lists[ipart][listend[ipart]].weight      = hg->vwgt[vertex];
      lists[ipart][listend[ipart]].vertex      = vertex;
      lists[ipart][listend[ipart]].source      = ipart;
      lists[ipart][listend[ipart]].destination = -1;
      lists[ipart][listend[ipart]].gain        = 0.0;
      ++listend[ipart];
    }
    for (ipart = 0; ipart < p; ipart++)
      qsort (lists[ipart], listend[ipart], sizeof (Vdata), comparison_weight);

    for (ipart = 0; ipart < p; ipart++)
      for (i = 0;  (part_weight[ipart] > max_weight) && (i < listend[ipart]); i++)  {
        part_weight[ipart] -= lists[ipart][i].weight;
        part[lists[ipart][i].vertex] = smallpart;
        }
  }

  /* algorithm loops to create interprocessor communication sections */
  for (loop = 0; loop < MAX_LOOP; loop++)   {
    int oddloop = loop & 1;         /* determines direction of legal moves */
    memset (listend, 0, p * sizeof (int));

    /* Calculate the total weights (local vertices weight) */
    /* most general case, need MPI global communication for total weight */
    total_weight = 0.0;
    memset (part_weight, 0, p * sizeof (double));
    if (hg->vwgt)
      for (vertex = 0; vertex < hg->nVtx; vertex++) {
        total_weight              += hg->vwgt[vertex];
        part_weight[part[vertex]] += hg->vwgt[vertex];
      }
    else  {
      total_weight = hg->nVtx;
      for (vertex = 0; vertex < hg->nVtx; vertex++)
        part_weight[part[vertex]] += 1.0;
    }
    max_weight = bal_tol * total_weight / p;

    /* Initial calculation of the cut distribution */
    memset (store1, 0, hg->nEdge * p * sizeof (int));
    for (edge = 0; edge < hg->nEdge; edge++)
      for (i = hg->hindex[edge]; i < hg->hindex[edge+1]; i++)
        ++((cuts[edge])[part[hg->hvertex[i]]]);

    for (vertex = 0; vertex < hg->nVtx; vertex++)  {
      /* calculate gains */
      for (i = 0; i < p; i++)
        gain [i] = 0.0;
      for (i = hg->vindex[vertex]; i < hg->vindex[vertex+1]; i++)  {
        edge = hg->vedge[i];
        for (ipart = 0; ipart < p; ipart++)  {
          if (ipart == part[vertex])
            continue;
          if ( (cuts[edge]) [ipart] != 0 && (cuts[edge]) [part[vertex]] == 1)
            gain[ipart] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
          if ( (cuts[edge]) [ipart] == 0)
            gain[ipart] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
        }
      }

      /* save best move, if any, for each vertex */
      /* oddloop, isabove, isbelow control move direction each pass */
      bestpart = -1;                            /* arbitrary illegal value */
      tgain    = -1.0;                          /* arbitrary small value */
      for (ipart = 0; ipart < p; ipart++)
        if (ipart != part[vertex] && gain[ipart] >= tgain
 /*      && !( (loop & 1)  ^  up[part[vertex]] [ipart]) */   )   {
          bestpart = ipart;
          tgain    = gain[ipart];
        }

      /* fill lists with the best, legal gain value per vertex */
      if (bestpart != -1 && tgain >= 0.0)  {
        (lists[bestpart]) [listend[bestpart]].weight      = hg->vwgt[vertex];
        (lists[bestpart]) [listend[bestpart]].vertex      = vertex;
        (lists[bestpart]) [listend[bestpart]].source      = part[vertex];
        (lists[bestpart]) [listend[bestpart]].destination = bestpart;
        (lists[bestpart]) [listend[bestpart]].gain        = tgain;
        (listend[bestpart])++;
      }
    } /* end of loop over all vertices */

    for (ipart = 0; ipart < p; ipart++)
      qsort (lists[ipart], listend[ipart], sizeof (Vdata), comparison_gain_weight_no);

    /* make moves while balance is OK */
    for (ipart = 0; ipart < p; ipart++)
      for (i = 0; i < listend[ipart]; i++)  {
        vertex = (lists[ipart]) [i].vertex;
        if (((hg->vwgt[vertex] + part_weight[ipart]) < max_weight)
         && (movect[vertex] < loop)) {
          part[vertex]        = ipart;
          part_weight[ipart] += hg->vwgt[vertex];
          movect[vertex] = loop+1;   /* sets when vertex may move again */
        }
      }
               
    /* communicate back moves made/rejected, update info (faked above) */
    /* update "ghost" vertices by either the previous comm is all to all, */
    /* or by a third comm by vertex owners to ghost owners */
    /* NOTE: this too is implicit in this serial version */

  }   /* end of loop over loop */
  
  Zoltan_Multifree(__FILE__,__LINE__, 8, &part_weight, &gain, &cuts, &store1,
   &store2, &lists, &listend, &movect);
  Zoltan_Multifree(__FILE__,__LINE__, 2, &upstore, &up);

  return ZOLTAN_OK;
}

   

/* note: sort is normally asending order, need the opposite!  */
static int comparison_gain_weight_no (const void *a, const void *b)
{
  if (((Vdata*)a)->gain < ((Vdata*)b)->gain)       return  1;
  if (((Vdata*)a)->gain > ((Vdata*)b)->gain)       return -1;

  if (((Vdata*)a)->weight < ((Vdata*)b)->weight)   return  1;
  if (((Vdata*)a)->weight > ((Vdata*)b)->weight)   return -1;

  if (((Vdata*)a)->vertex > ((Vdata*)b)->vertex)   return  1;
  if (((Vdata*)a)->vertex < ((Vdata*)b)->vertex)   return -1;

  return 0;
}



static int comparison_weight (const void *a, const void *b)
{
  if (((Vdata*)a)->weight > ((Vdata*)b)->weight)   return  1;
  if (((Vdata*)a)->weight < ((Vdata*)b)->weight)   return -1;

  return 0;
}



/******************************************************************************/

int Zoltan_PHG_move_vertex (PHGraph *hg, int vertex, int sour, int dest,
 int *part, int **cut, double *gain, HEAP *heap)
{
  int i, j, edge, v;

  gain[vertex] = 0.0;
  part[vertex] = dest;

  for (i = hg->vindex[vertex]; i < hg->vindex[vertex+1]; i++) {
    edge = hg->vedge[i];
    if (cut[sour][edge] == 1) {
      for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
        v = hg->hvertex[j];
        gain[v] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
        if (heap)
          Zoltan_heap_change_value(&heap[part[v]], v, gain[v]);
      }
    }
    else if (cut[sour][edge] == 2) {
      for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
        v = hg->hvertex[j];
        if (part[v] == sour) {
          gain[v] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
          if (heap)
            Zoltan_heap_change_value(&heap[part[v]], v, gain[v]);
          break;
        }
      }
    }

    if (cut[dest][edge] == 0) {
      for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
        v = hg->hvertex[j];
        gain[v] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
        if (heap)
          Zoltan_heap_change_value(&heap[part[v]], v, gain[v]);
      }
    }
    else if (cut[dest][edge] == 1) {
      for (j = hg->hindex[edge]; j < hg->hindex[edge+1]; j++) {
        v = hg->hvertex[j];
        if (v != vertex && part[v] == dest) {
          gain[v] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
          if (heap)
            Zoltan_heap_change_value(&heap[part[v]], v, gain[v]);
          break;
        }
      }
    }
    cut[sour][edge]--;
    cut[dest][edge]++;
  }
  return ZOLTAN_OK;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
