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

#include <stdarg.h>
#include "phg.h"


#define _DEBUG
#define _DEBUG2
    
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

/*************************************************************************
* -------------------------- Error Exit ----------------------------------
**************************************************************************/
void errexit(char *f_str,...)
{
va_list argp;

fflush(stdout);
fflush(stderr);
fprintf(stderr, "\n****** Error:\n");
va_start(argp, f_str);
vfprintf(stderr, f_str, argp);
va_end(argp);

fprintf(stderr," ******\n");
fflush(stderr);
exit(1);
}


typedef int  (*SelectFunc)(HEAP heap[2], double *weights, double *max_weight);

static int fm2_select(HEAP heap[2], double *weights, double *max_weight)
{
    int from;
    /* select a vertex with max gain; if possible */
    if (Zoltan_heap_not_empty(&heap[0]) && Zoltan_heap_not_empty(&heap[1])) {
        if (Zoltan_heap_max_value(&heap[0])==Zoltan_heap_max_value(&heap[1]))
            from = (weights[0] < max_weight[0]) ? 1 : 0;
        else
            from = (Zoltan_heap_max_value(&heap[0])>Zoltan_heap_max_value(&heap[1])) ? 0 : 1;
    } else if (Zoltan_heap_empty(&heap[0])) {
        if (Zoltan_heap_empty(&heap[1])) /* too bad both are empty */
            return -1; /* nothing to select */
        else
            from = 1;
    } else
        from = 0;
    return Zoltan_heap_extract_max(&heap[from]);    
}


static void fm2_move_vertex(int v, PHGraph *hg, Partition part, float *gain, HEAP *heap, int *pins[2], int *lpins[2], double *weights, double *lweights, int *mark, int *adj)
{
    float oldgain=gain[v];
    int   pno=part[v], vto=1-pno, adjsz=0, j, i;
    
    mark[v] = 1;
    part[v] = vto;
    weights[pno] -= (hg->vwgt ? hg->vwgt[v] : 1.0);
    weights[vto] += (hg->vwgt ? hg->vwgt[v] : 1.0);
    lweights[pno] -= (hg->vwgt ? hg->vwgt[v] : 1.0);
    lweights[vto] += (hg->vwgt ? hg->vwgt[v] : 1.0);

    for (j = hg->vindex[v]; j < hg->vindex[v+1]; j++) {
        int n = hg->vedge[j];
        int w = hg->ewgt ? hg->ewgt[n] : 1.0;
    
        --pins[pno][n];
        --lpins[pno][n];
#ifdef _DEBUG2
        if (pins[pno][n] < 0)
            errexit("move of %d makes pin[%d][%d]=%d", v, pno, n, pins[pno][n]);
#endif

        if (!pins[pno][n]) {  /* no pin in source part */ 
            for (i = hg->hindex[n]; i < hg->hindex[n+1]; ++i) {
                int u = hg->hvertex[i]; 
                gain[u] -= w;
                if (!mark[u]) {
                    adj[adjsz++] = u;
                    mark[u] = -1;
		}
	    }
        } else if (pins[pno][n]==1) {
            for (i = hg->hindex[n]; i < hg->hindex[n+1]; ++i) {
                int u = hg->hvertex[i]; 
                if (part[u]==pno) {
                    gain[u] -= w;
                    if (!mark[u]) {
                        adj[adjsz++] = u;
                        mark[u] = -1;
                    }
                }
            }
        }

        ++pins[vto][n];
        ++lpins[vto][n];
        if (pins[vto][n]==1) { /* now there is at least one pin here */
            for (i = hg->hindex[n]; i < hg->hindex[n+1]; ++i) {
                int u = hg->hvertex[i];
                if (part[u]==pno) { 
                    gain[u] += w;
                    if (!mark[u]) {
                        adj[adjsz++] = u;
                        mark[u] = -1;
                    }
                }
            }
        } else if (pins[vto][n]==2)
            for (i = hg->hindex[n]; i < hg->hindex[n+1]; ++i) {
                int u = hg->hvertex[i];
                if (part[u]==vto) { 
                    gain[u] -= w;
                    if (!mark[u]) {
                        adj[adjsz++] = u;
                        mark[u] = -1;
                    }
                }
            }
    }

    gain[v] = -oldgain;    
    for (i=0; i<adjsz; i++) {
        int u=adj[i], p=part[u];
        
        mark[u] = 0;
        if (Zoltan_heap_has_elem(&heap[p], u))
            Zoltan_heap_change_value(&heap[p], u, gain[u]);
    }
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
    int    i, j,  *pins[2], *lpins[2], *moves=0, *mark=0, *adj=0, round=0, best_cutsizeat, cont;
    double total_weight, weights[2], lweights[2], max_weight[2];
    double cutsize, best_cutsize, ratio = hg->ratio, best_imbalance, imbalance;
    float  *gain=0, *lgain=0;
    HEAP   heap[2];
    char   *yo="local_fm2";
    PHGComm *hgc=&hgp->comm;

    SelectFunc select_func = fm2_select;
        
#ifdef _DEBUG    
    FILE *fp;
    char fname[512];

    sprintf(fname, "uvc-phg-debug.%d.txt", zz->Proc);    
    fp = fopen(fname, "w");
    if (!fp) {
        printf ("failed to open debug file '%s'\n", fname);
        return ZOLTAN_FATAL;
    }
    fprintf(fp, "%s\n", uMe(hgc));
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
    lweights[0] = lweights[1] = 0.0;
    if (hg->vwgt)  
        for (i = 0; i < hg->nVtx; i++)
            lweights[part[i]] += hg->vwgt[i];
    else  
        for (i = 0; i < hg->nVtx; i++)
            lweights[part[i]] += 1.0;

    MPI_Allreduce(lweights, weights, 2, MPI_DOUBLE, MPI_SUM, hgc->row_comm);
    total_weight = weights[0] + weights[1];
    max_weight[0] = total_weight * bal_tol *      ratio;
    max_weight[1] = total_weight * bal_tol * (1 - ratio);

    printf("%s W:(%.2lf, %.2lf) tol: %.3f MaxW:(%.2lf, %.2lf)\n", uMe(hgc), weights[0], weights[1], bal_tol, max_weight[0], max_weight[1]);
    
    if (!(pins[0]     = (int*) ZOLTAN_CALLOC(2 * hg->nEdge, sizeof(int)))
        || !(lpins[0] = (int*) ZOLTAN_CALLOC(2 * hg->nEdge, sizeof(int)))
        || !(mark     = (int*)   ZOLTAN_CALLOC(hg->nVtx,  sizeof(int)))
        || !(adj      = (int*)   ZOLTAN_CALLOC(hg->nVtx,  sizeof(int)))                
        || !(moves    = (int*)   ZOLTAN_CALLOC(hg->nVtx,  sizeof(int)))
        || !(gain     = (float*) ZOLTAN_CALLOC(hg->nVtx,  sizeof(float)))
        || !(lgain    = (float*) ZOLTAN_CALLOC(hg->nVtx,  sizeof(float))) ) {
         Zoltan_Multifree(__FILE__,__LINE__, 7, &pins[0], &lpins[0], &mark, &adj, &moves, &gain, &lgain);
         ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
         return ZOLTAN_MEMERR;
    }
    pins[1] = &(pins[0][hg->nEdge]);
    lpins[1] = &(lpins[0][hg->nEdge]);

    /* Initial calculation of the pins distribution  */
    for (i = 0; i < hg->nEdge; ++i)
        for (j = hg->hindex[i]; j < hg->hindex[i+1]; ++j)
            ++(lpins[part[hg->hvertex[j]]][i]);

    do {
        int v=1, movecnt=0, neggaincnt=0;
        int maxneggain = (hgp->fm_max_neg_move < 0) ? hg->nVtx : hgp->fm_max_neg_move;

        MPI_Allreduce(lpins[0], pins[0], 2*hg->nEdge, MPI_INT, MPI_SUM, hgc->row_comm);
        best_cutsizeat=0;
        cutsize = 0.0;
        for (i=0; i < hg->nEdge; ++i) {
            if (pins[0][i] && pins[1][i])
                cutsize += (hg->ewgt ? hg->ewgt[i] : 1.0);
        }
        MPI_Allreduce(&cutsize, &best_cutsize, 1, MPI_DOUBLE, MPI_SUM, hgc->col_comm);
        cutsize = best_cutsize;
        imbalance = MAX (weights[0]-max_weight[0], weights[1]-max_weight[1]);
    
#ifdef _DEBUG
    /* Just for debugging */
        best_cutsize = Zoltan_PHG_hcut_size_total(hgc, hg, part, p);
        printf("%s: Initial cutsize=%.2lf Verify: total=%.2lf\n", uMe(hgc), cutsize,
               best_cutsize);
        if (best_cutsize!=cutsize) {
            ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Invalid cut!!!");
        }
        /* debuggging code ends here */
#endif
        
        for (i = 0; i < hg->nVtx; i++)
            for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++) {
                int edge = hg->vedge[j];
                if (pins[part[i]][edge] == 1)
                    lgain[i] += (hg->ewgt ? hg->ewgt[edge] : 1.0);
                else if (pins[1-part[i]][edge] == 0)
                    lgain[i] -= (hg->ewgt ? hg->ewgt[edge] : 1.0);
            }
        /* now sum up all gains on only first row of processors */
        MPI_Reduce(lgain, gain, hg->nVtx, MPI_FLOAT, MPI_SUM, 0, hgc->col_comm);

        best_cutsize = cutsize;
        best_imbalance = imbalance;
        printf("%s: FM PassCnt= %d  Cut=%.2lf  Imbal= %.2lf\n", uMe(hgc), round+1, cutsize, imbalance);

        if (!hgc->myProc_y) { /* those are the lucky ones; each proc in column-group
                                 could have compute the same moves concurrently; but for this
                                 version we'll do in in the root procs and broadcast */
            /* Initialize the heaps and fill them with the gain values */
            Zoltan_heap_init(zz, &heap[0], hg->nVtx);
            Zoltan_heap_init(zz, &heap[1], hg->nVtx);  
            for (i = 0; i < hg->nVtx; ++i)
                Zoltan_heap_input(&heap[part[i]], i, gain[i]);
            Zoltan_heap_make(&heap[0]);
            Zoltan_heap_make(&heap[1]);
            
            while ((v>=0) && (neggaincnt < maxneggain)) {
                int from, to;

                v = (*select_func)(heap, weights, max_weight);
                if (v<0) /* there was nothing to select */
                    break;
                from = part[v];
                to = 1-from;

                mark[v]++;
                if (weights[to]+((hg->vwgt)? hg->vwgt[v] : 1.0) > max_weight[to]) {
                    printf("%s %4d: %6d (g: %5.1lf), p:%2d [%4.0lf, %4.0lf] NF\n", uMe(hgc), movecnt, v, gain[v], from, weights[0], weights[1]);
                    moves[movecnt++] = -(v+1);
                    continue;
                } 
	
                moves[movecnt] = (v+1);
                ++neggaincnt;
                cutsize -= gain[v];

                fm2_move_vertex(v, hg, part, gain, heap, pins, lpins, weights, lweights, mark, adj);
                imbalance = MAX (weights[0]-max_weight[0], weights[1]-max_weight[1]);
                if ((cutsize<best_cutsize) || (cutsize==best_cutsize && imbalance < best_imbalance)) {
                    printf("%s %4d: %6d (g: %5.1lf), p:%2d [%4.0lf, %4.0lf] %.1lf <-- Best\n", uMe(hgc), movecnt, v, -gain[v], from, weights[0], weights[1], cutsize); /* after move gain is -oldgain */
                    
                    best_cutsize = cutsize;
                    best_cutsizeat = movecnt+1;
                    best_imbalance = imbalance;
                    neggaincnt = 0;
                } else
                    printf("%s %4d: %6d (g: %5.1lf), p:%2d [%4.0lf, %4.0lf] %.1lf\n", uMe(hgc), movecnt, v, -gain[v], from, weights[0], weights[1], cutsize);
                ++movecnt;
            }

            printf("%s END of Pass bestcut= %.1lf @ %d curcut: %.1lf\n", uMe(hgc), best_cutsize, best_cutsizeat, cutsize); 
            /* roll back the moves without any improvement */
            for (i=movecnt-1; i>=best_cutsizeat; --i) {
                int v = moves[i];
                if (v>0)
                    fm2_move_vertex(v-1, hg, part, gain, heap, pins, lpins, weights, lweights, mark, adj);
            }
            printf("%s Roll back completed!\n", uMe(hgc));
        }
        
        printf("%s before broadcast best_cutsizeat=%d\n", uMe(hgc), best_cutsizeat);
        /* now root bcast moves to column procs */
        MPI_Bcast(&best_cutsizeat, 1, MPI_INT, 0, hgc->col_comm);
        MPI_Bcast(moves, best_cutsizeat, MPI_INT, 0, hgc->col_comm);
        printf("%s after bcast best_cutsizeat=%d\n", uMe(hgc), best_cutsizeat);
        if (hgc->myProc_y) { /* now non-root does move simulation */
            for (i=0; i<best_cutsizeat; ++i) {
                int v = moves[i];
                if (v>0)
                    fm2_move_vertex(v-1, hg, part, gain, heap, pins, lpins, weights, lweights, mark, adj);
            }
        }
        MPI_Allreduce(lweights, weights, 2, MPI_DOUBLE, MPI_SUM, hgc->row_comm);
        cont = 0;
        MPI_Allreduce(&best_cutsizeat, &cont, 1, MPI_INT, MPI_LOR, hgc->row_comm);
        printf("%s after allreduce best_cutsizeat=%d,  cont=%d,  round=%d  limit=%d\n", uMe(hgc), best_cutsizeat, cont, round, hgp->fm_loop_limit);
    } while (cont &&  (round++ < hgp->fm_loop_limit));

    Zoltan_Multifree(__FILE__, __LINE__, 7, &pins[0], &lpins[0], &mark, &adj, &moves, &gain, &lgain);
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
      /* int oddloop = loop & 1; */        /* determines direction of legal moves */
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
