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

static ZOLTAN_HG_GLOBAL_PART_FN global_ran;
static ZOLTAN_HG_GLOBAL_PART_FN global_lin;
static ZOLTAN_HG_GLOBAL_PART_FN global_bfs;
static ZOLTAN_HG_GLOBAL_PART_FN global_bfsr;

/****************************************************************************/

ZOLTAN_HG_GLOBAL_PART_FN *Zoltan_HG_Set_Global_Part_Fn(char *str)
{
  if      (!strcasecmp(str, "ran")) return global_ran;
  else if (!strcasecmp(str, "lin")) return global_lin;
  else if (!strcasecmp(str, "bfs")) return global_bfs;
  else if (!strcasecmp(str, "bfsr")) return global_bfsr;
  else                              return NULL;
}

/****************************************************************************/

/* Zoltan_HG_Global computes a global partitioning of a hypergraph.
 * Typically, this routine is called at the bottom level in a 
 * multilevel scheme (V-cycle).
 */ 

int Zoltan_HG_Global (ZZ *zz, HGraph *hg, int p, Partition part, HGPartParams *hgp)
{
  return hgp->global_part(zz,hg,p,part);
}

/****************************************************************************/

/* Sequence partitioning on the vertices of a hypergraph
   in a given order. Currently, even partition sizes
   are assumed. Multi-weights are not yet supported. 

   This function is called by global_lin and global_ran. 

   EB: This is a quick heuristic. We could alternatively use 
   a more expensive but optimal algorithm, see e.g. Ali Pinar's thesis. */

static int seq_part (ZZ *zz, HGraph *hg, int *order, int p, Partition part)
{ 
  int i, j, number;
  float weight_sum = 0.0, part_sum = 0.0, old_sum, cutoff;

  /* First sum up all the weights. */
  if (hg->vwgt){
    for (i=0; i<hg->nVtx; i++)
      weight_sum += hg->vwgt[i];
  }
  else
    weight_sum = (float)hg->nVtx;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("GLOBAL_PART weight_sum=%f, avg=%f\n",weight_sum, weight_sum/hg->nVtx);
  number = 0; /* Assign next vertex to partition no. number */
  cutoff = weight_sum/p;  /* Cutoff for current partition */
  for (i=0; i<hg->nVtx; i++)
  { 
    /* If order==NULL, then use linear order. */
    j = order ? order[i] : i;
    part[j] = number;
    old_sum = part_sum;
    part_sum += hg->vwgt?hg->vwgt[j]:1.0;
    /* Check if we passed the cutoff and should start a new partition */
    if ((number+1)<p && part_sum > cutoff){
      number++;
      /* Decide if current vertex should be moved to the next partition */
      if (part_sum-cutoff > cutoff-old_sum){
        part[j]++;
        part_sum = old_sum; 
      }
      weight_sum -= part_sum;
      cutoff = weight_sum/(p-number);
      if (part[j] == number)
        part_sum = hg->vwgt?hg->vwgt[j]:1.0;
      else
        part_sum = 0.0;
    }
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("GLOBAL_PART i=%2d, part[%2d] = %2d, part_sum=%f\n",i,j,part[j],part_sum);
  }

  return ZOLTAN_OK;
}

/****************************************************************************/

/* Linear partitioning. Sequence partitioning with vertices in linear order. */

static int global_lin (ZZ *zz, HGraph *hg, int p, Partition part)
{ 
  /* Call sequence partitioning with no order array. */
  return seq_part( zz, hg, NULL, p, part);
}

/****************************************************************************/

/* Random partitioning. Sequence partitioning with vertices in random order. */
static int global_ran (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int i, ierr, number, temp, *order=NULL;
  char *yo = "global_ran" ;

  if (!(order  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))
  { ZOLTAN_FREE ((void **) &order) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }
  for (i=0; i<hg->nVtx; i++)
    order[i] = i;

  /* Randomly permute order array */
  for (i=hg->nVtx; i>0; i--)
  { number=Zoltan_HG_Rand()%i;
    temp = order[number];
    order[number] = order[i-1];
    order[i-1] = temp;
  }
 
  /* Call sequence partitioning with random order array. */
  ierr = seq_part( zz, hg, order, p, part);

  ZOLTAN_FREE ((void **) &order);
  return (ierr);
}


/****************************************************************************/

/* BFS partitioning. Sequence partitioning with vertices 
   in breadth-first search order. */

static int global_bfs (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int i, ierr, start, *order=NULL;
  static int bfs_order();
  char *yo = "global_bfs" ;

  if (!(order  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))
  { ZOLTAN_FREE ((void **) &order) ;
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    return ZOLTAN_MEMERR;
  }

  /* Find pseudo-peripheral start vertex */
  start = Zoltan_HG_Rand()%(hg->nVtx);
  for (i=0; i<2; i++){
    ierr = bfs_order(zz, hg, order, &start, NULL, NULL, NULL);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      goto error;
    start = order[hg->nVtx -1];
  }


  /* Compute BFS order */
  ierr = bfs_order(zz, hg, order, &start, NULL, NULL, NULL);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
    goto error;

  /* Call sequence partitioning with BFS order array. */
  ierr = seq_part( zz, hg, order, p, part);

error:
  ZOLTAN_FREE ((void **) &order);
  return (ierr);
}
 
/****************************************************************************/

/* Compute BFS order on a hypergraph. 
   order[0] is the first vertex, order[1] the next, etc. */

static int bfs_order (
  ZZ *zz, 
  HGraph *hg,		/* Hypergraph. */
  int *order,		/* Order array. May be partially filled on entry. */
 			/*           On exit, order[i] is the i'th vertex. */
  int *start_vtx,	/* Start the BFS from this vertex. */
                        /*           On exit, an unmarked vertex (or -1). */
  int *start_index,	/* Optional: Number vertices starting at this value. */
  			/*           On exit, highest ordered vertex + 1 */
  int *rank,            /* Optional: Rank of vertices numbered so far (array). */
  float *cutoff		/* Optional: Number vertices until the sum > cutoff. */
			/*           On exit, actual size of labeled piece. */
)
{
  int i, j, vtx, edge, number, nbor, unmarked, free_rank; 
  int first, last;
  int ierr=ZOLTAN_OK;
  float weight_sum = 0.0;
  char msg[128];
  static char *yo = "bfs_order";

/*
    Breadth first search algorithm:
    --------------------------------
    unmark all vertices 
    num = 0
    order[start_vtx] = num
    queue Q = start_vtx
    while Q nonempty and size < cutoff
      choose a vertex v from front of queue
      order[v] = num++
      remove v from front of queue 
      for each hyperedge v shares
        for each unmarked neighbor w
          add w to end of queue
*/

  number = start_index ? (*start_index) : 0; /* Assign next vertex this number */

  /* If rank array is given, use it, otherwise allocate it here. */
  if (rank){
    free_rank = 0;
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
      /* Verify input data. */
      for (i=0; i<number; i++)
        if (rank[order[i]] != i){
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Input arrays order and rank are inconsistent.");
          ierr = ZOLTAN_FATAL;
          goto error;
        }
    }
  }
  else { /* !rank */
    free_rank = 1;
    if (!(rank  = (int *)   ZOLTAN_MALLOC (sizeof (int) * hg->nVtx))) {
      ZOLTAN_FREE ((void **) &rank) ;
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
    }
    for (i=0; i<hg->nVtx; i++)
      rank[i] = -1;  /* -1 means this vtx has not yet been numbered */
    /* Fill in rank information from partial order array. */
    for (i=0; i<number; i++)
      rank[order[i]] = i;
  }

  /* Use order array as a queue. Put start_vtx in queue. */
  first = number;
  last = first+1;
  if (start_vtx)
    order[first] = *start_vtx;
  else
    ZOLTAN_PRINT_ERROR(0, yo, "Input start_vtx is NULL.");
 
  /* printf("Starting new BFS at vertex %d\n", *start); */
  unmarked = 0; /* Next unmarked vertex */

  while (number < hg->nVtx && !(cutoff && (weight_sum >= *cutoff))) {
    /* Is queue empty? */
    if (last == first){
      ZOLTAN_PRINT_WARN(0, yo, "Queue is empty; hypergraph must be disconnected");
      /* printf("debug: first = last = %d, value is %d\n", first, order[first]);*/

      /* Find an unmarked vertex to put in queue */
      while (unmarked<hg->nVtx && rank[unmarked]>=0) unmarked++;
      if (unmarked==hg->nVtx){
        ZOLTAN_PRINT_ERROR(0, yo, "All vertices looks to be visited, but that cant be!");
        ierr = ZOLTAN_FATAL;
        goto error;
      }
      order[last++] = unmarked;
    }
    /* Get next vertex from queue */
    vtx = order[first++];
    if (rank[vtx]<0)
      rank[vtx] = number++;
    else{
      sprintf(msg, "Vertex %d in queue already labeled", vtx);
      ZOLTAN_PRINT_ERROR(0, yo, msg);
      ierr = ZOLTAN_FATAL;
      goto error;
    }
    weight_sum += hg->vwgt?hg->vwgt[vtx]:1.0;
     
    /* Add nbors to queue. */
    /* Possible variation: pick heaviest edge first. */
    for (j=hg->vindex[vtx]; j<hg->vindex[vtx+1]; j++){
      edge = hg->vedge[j];
      for (i=hg->hindex[edge]; i<hg->hindex[edge+1]; i++){
        nbor = hg->hvertex[i];
        /* printf("debug: vtx=%2d, nbor=%2d, rank[nbor] = %2d\n", vtx, nbor, rank[nbor]); */
        if (rank[nbor] == -1){
          order[last++] = nbor;
          if (last > hg->nVtx) {
            ZOLTAN_PRINT_ERROR(0, yo, "Queue is full");
            ierr = ZOLTAN_FATAL;
            goto error;
          }
          rank[nbor] = -2; /* nbor is now in queue */
        }
      }
    }
  }

  /* Order should be the inverse permutation of rank. */
  for (i=0; i<hg->nVtx; i++){
    if (rank[i]>=0) 
      if (order[rank[i]] != i)
         ZOLTAN_PRINT_WARN(0, yo, "Arrays order and rank are inconsistent.");
    /* Clean out queue for possible future bfs. */
    if (rank[i] == -2) rank[i] = -1; 
  }

/*
  printf("BFS order = ");
  for (i=0; i<hg->nVtx; i++) printf("%3d ", order[i]);
  printf("\n");
*/

  /* Update return arguments. */
  if (start_vtx)    *start_vtx = (first<hg->nVtx ? order[first] : -1);
  if (start_index)  *start_index = number;
  if (cutoff)       *cutoff = weight_sum;

error:
  if (free_rank)
    ZOLTAN_FREE ((void **) &rank);
  return ierr;
}

/****************************************************************************/

/* BFS partitioning with restart. Sequence partitioning with vertices 
   in breadth-first search order, breaking of pieces as we go; 
   that is, the BFS is restarted for each partition. */

static int global_bfsr (
  ZZ *zz, 
  HGraph *hg,
  int p,
  Partition part
)
{ 
  int i, j, number, start, index, old_index, *rank, *order;
  int ierr = ZOLTAN_OK;
  float weight_sum = 0.0, cutoff;
  char *yo = "global_bfsr" ;

  if ((!(order  = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx))) ||
      (!(rank   = (int *) ZOLTAN_MALLOC (sizeof (int) * hg->nVtx)))) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
    ierr = ZOLTAN_MEMERR;
    goto error;
  }

  /* Compute weight sum. */
  if (hg->vwgt)
  { for (i=0; i<hg->nVtx; i++)
      weight_sum += hg->vwgt[i];
  }
  else
    weight_sum = (float)hg->nVtx;

  if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
    printf("GLOBAL_PART weight_sum:%f\n",weight_sum);

  /* Find pseudo-peripheral start vertex */
  start = Zoltan_HG_Rand()%(hg->nVtx);
  for (j=0; j<2; j++){
    ierr = bfs_order(zz, hg, order, &start, NULL, NULL, NULL);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN)
      goto error;
    start = order[hg->nVtx -1];
  }

  /* Unmark all vertices */
  for (i=0; i<hg->nVtx; i++)
    rank[i] = -1;  /* -1 means this vtx has not yet been assigned to a part */

  /* Do BFS until right size, then break off a partition. Repeat. */
  index = 0;
  for (number=0; number<p-1; number++){
    cutoff = weight_sum/((float)(p-number));
    old_index = index;
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("GLOBAL_PART before bfs: number=%2d, start = %2d, index = %2d, cutoff = %f\n", number, start, index, cutoff);
    ierr = bfs_order(zz, hg, order, &start, &index, rank, &cutoff);
    if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
      printf("GLOBAL_PART after bfs: number=%2d, start = %2d, index = %2d, cutoff = %f\n", number, start, index, cutoff);
    weight_sum -= cutoff;
    /* Assign partition number to newly ordered vertices */
    for (i=old_index; i<index; i++){
      part[order[i]] = number;
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL)
        printf("GLOBAL_PART i=%2d, part[%2d]= %2d\n",i,order[i],part[order[i]]);
    }
  }
  /* Put remaining unmarked vertices in the last partition. */
  for (i=0; i<hg->nVtx; i++)
    if (rank[i]<0) part[i] = p-1;

error:
  /* Free data and return. */
  ZOLTAN_FREE ((void **) &rank) ;
  ZOLTAN_FREE ((void **) &order) ;
  return ierr;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
