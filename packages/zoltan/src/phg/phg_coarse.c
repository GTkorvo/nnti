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


/* Procedure to coarsen a hypergraph based on a matching. All vertices of one
   match are clustered to one vertex. Identical hyperedges are collapsed to a
   single hyperedge with combined weight. The array LevelMap is the mapping of
   the old vertices to the new vertices. It will be used to pass a partition
   of the coarse graph back to the original graph.                         */
   
int Zoltan_PHG_Coarsening
( ZZ       *zz,         /* the Zoltan data structure */
  HGraph  *hg,         /* information about hypergraph, weights, etc. */
  int      *match,      /* Matching, Packing or Grouping array */
  HGraph  *c_hg,       /* points to a working copy of hg structure */
  int      *LevelMap)   /* information to reverse coarsenings later */
{
  int i, j, vertex, edge, *ip, me, size, count;
  int *cmatch=NULL, *list=NULL, *used_edges=NULL, *c_vindex=NULL, *c_vedge=NULL;
  float *c_ewgt=NULL, *pwgt;
  char *buffer=NULL, *rbuffer=NULL;
  int *displs=NULL, *each_size=NULL;
  PHGComm *hgc = hg->comm;
  char *yo = "Zoltan_PHG_Coarsening";

  ZOLTAN_TRACE_ENTER (zz, yo);
  
  Zoltan_HG_HGraph_Init (c_hg);   /* inits working copy of hypergraph info */
  c_hg->info  = hg->info + 1;      /* for debugging */
  c_hg->ratio = hg->ratio;         /* for "global" recursive bisectioning */
  c_hg->redl  = hg->redl;          /* to stop coarsening near desired count */
    
  if (!(cmatch    = (int*) ZOLTAN_MALLOC (hg->nVtx     * sizeof(int)))
   || !(list      = (int*) ZOLTAN_MALLOC (hg->nVtx     * sizeof(int)))
   || !(displs    = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int)))
   || !(each_size = (int*) ZOLTAN_MALLOC (hgc->nProc_x * sizeof(int))))  {
     Zoltan_Multifree (__FILE__, __LINE__, 4, &cmatch, &displs, &list, 
      &each_size);     
     ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
     ZOLTAN_TRACE_EXIT (zz, yo);
     return ZOLTAN_MEMERR;
  }     
  for (i = 0; i < hg->nVtx; i++)
     cmatch[i] = match[i];         /* working copy of match array */

  /* Assume all rows in a column have the entire (column's) matching info */
  /* Calculate the number of coarse vertices. */
  c_hg->nVtx = 0;                 /* counts number of new (coarsened) vertices */
  me = hgc->myProc_x;             /* short name, convenience variable */
  size = 0;                       /* size (in ints) to communicate */
  count = 0;                      /* number of messages to communicate */
  for (i = 0; i < hg->nVtx; i++)  {    /* loop over every local vertice */
    if (match[i] < 0)  {               /* external processor match */
      int proc, gx = -match[i]-1;
      
      /* rule to determine "ownership" of coarsened vertices between procs */
      proc = ((gx + VTX_LNO_TO_GNO (hg,i)) & 1) ? MIN(gx, me) : MAX(gx, me);
      if (proc != me)   {             /* another processor owns this vertex */
        size += hg->vindex[i+1] - hg->vindex[i];  /* send buffer sizing */ 
        list[count++] = i;                        /* list of vtx's to send */
        }
      else 
        c_hg->nVtx++;         /* myProc owns the matching across processors */ 
    }
    else if (cmatch[i] >= 0)  { /* allow (local only) packing and groupings */    
      c_hg->nVtx++;
      vertex = i;
      while (cmatch[vertex] >= 0)  {
        j              =  cmatch[vertex];
        cmatch[vertex] = -cmatch[vertex] - 1;  /* flag this as done already */     
        vertex         =  j;
      }
    }
  }

  /* size and allocate the send buffer */
  /* fix for size/count = 0 */
  size += (4*count);          /* add gno's, weight, edge counts to message size */
  if (size > 0)
    if (!(buffer = (char*) ZOLTAN_MALLOC (size * sizeof(int)))  
     || !(c_hg->vwgt = (float*) ZOLTAN_CALLOC (c_hg->nVtx, sizeof(float))))  {
       Zoltan_Multifree (__FILE__, __LINE__, 2, c_hg->vwgt, buffer);
       ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
       ZOLTAN_TRACE_EXIT (zz, yo);
       return ZOLTAN_MEMERR;
    }
        
  /* Message is list of <gno, gno's edge count, list of edge gno's> */
  ip = (int*) buffer;
  for (i = 0; i < count; i++)  {
    *ip++ = 98765;
    *ip++ = VTX_LNO_TO_GNO (hg, list[i]);           /* destination vertex gno */ 
    pwgt = (float*) ip++;                                    /* vertex weight */
    *pwgt = (hg->vwgt == NULL) ? 1 : hg->vwgt[list[i]] ;
    *ip++ = hg->vindex[list[i]+1] - hg->vindex[list[i]];             /* count */                                                             /* weights??? */
    for (j = hg->vindex[list[i]]; j < hg->vindex[list[i]+1]; j++)
      *ip++ = EDGE_LNO_TO_GNO (hg, hg->vedge[j]);                    /* edges */
  }    
    
  MPI_Allgather (&size, 1, MPI_INT, each_size, 1, MPI_INT, hgc->row_comm);  

  size = 0;
  for (i = 0; i < hgc->nProc_x; i++)
    size += each_size[i];
    
  if (size > 0)  {  
    if (!(rbuffer = (char*) ZOLTAN_MALLOC (size * sizeof(int))))   {
      ZOLTAN_TRACE_EXIT (zz, yo);
      ZOLTAN_PRINT_ERROR (zz->Proc, yo, "Insufficient memory.");
      return ZOLTAN_MEMERR;
    }
  }  

  displs[0] = 0;
  for (i = 1; i < hgc->nProc_x; i++)
     displs[i] = displs[i-1] + each_size[i-1];
                
  MPI_Allgatherv (buffer, count, MPI_INT, rbuffer, each_size, displs, MPI_INT,
   hgc->row_comm);            

  /* index all received data for rapid lookup */
  ip = (int*) rbuffer;
  for (i = 0; i < size; i++)  {
    if (ip[i++] != 98765)
       uprintf (hgc, "RTHRTH alignment error, i %d of %d\n", i, size);
       
    if (VTX_TO_PROC_X (hg, ip[i]) == me)
      cmatch[VTX_GNO_TO_LNO(hg,ip[i])] = i;
    i++;                /* destination gno */

    i++;                /* vertex weight */


    i += (ip[i]+1);     /* count of hyperedges + hyperedges */
  }
     
  if (!(used_edges = (int*)   ZOLTAN_CALLOC (hg->nEdge,       sizeof(int)))
   || !(c_hg->vwgt = (float*) ZOLTAN_CALLOC (hg->nVtx,        sizeof(float)))
   || (hg->ewgt && !(c_ewgt = (float*) ZOLTAN_MALLOC (hg->nEdge * sizeof(float))))
   || !(c_vindex   = (int*)   ZOLTAN_MALLOC ((hg->nVtx+1)   * sizeof(int)))
   || !(c_vedge    = (int*)   ZOLTAN_MALLOC (hg->nPins   * sizeof(int)))) {
      Zoltan_Multifree (__FILE__, __LINE__, 7, c_hg->vwgt,
       buffer, rbuffer, &used_edges, &c_ewgt, &c_vindex, &c_vedge);
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
  }

  if (c_ewgt)   /* TODO: if we collapse edges as well we need to compute actual weights*/
      for (i=0; i<hg->nEdge; ++i)
          c_ewgt[i] = hg->ewgt[i];


  /* Construct the LevelMap; match[vertex] is changed back to original value */
  /* Coarsen vertices (create vindex, vedge), sum up coarsened vertex weights */
  c_hg->nPins = 0;                       /* count of coarsened pins */
  c_hg->nVtx     = 0;                       /* count of coarsened vertices */
  for (i = 0; i < hg->nVtx; i++)  {
    if (match[i] < 0 && cmatch[i] < 0)      /* match to external vertex */                  
       LevelMap [i] = match[i];             /* negative value => external vtx */         
    else if (match[i] < 0) {                /* match from external vertex */
       LevelMap[i] = c_hg->nVtx;
      
       c_vindex[c_hg->nVtx] = c_hg->nPins;
       ip =  ((int*) rbuffer) + cmatch[i];      /* point to receieved data */ 
       ip++;                                    /* skip over vtx gno */
       pwgt = (float*) &ip[i++];
       c_hg->vwgt[c_hg->nVtx] = *pwgt;  
       count = *ip++;           /* extract edge count, advance to first edge */
       while (count--)  {
          edge = EDGE_GNO_TO_LNO (hg, *ip++);
          used_edges [edge]         = i+1;
          c_vedge[c_hg->nPins++] = edge;
       }
       
       c_hg->vwgt[c_hg->nVtx] += hg->vwgt ? hg->vwgt[vertex] : 1.0;
       for (j = hg->vindex[i]; j < hg->vindex[i+1]; j++)  {
         if (used_edges [hg->vedge[j]] <= i)  {
           used_edges [hg->vedge[j]]     = i+1;          
           c_vedge[c_hg->nPins++] = hg->vedge[j];
         }      
       }        
       c_hg->nVtx++;          
    }
    else if (match[i] >= 0 && cmatch[i] < 0)  { /* match/pack/group local vtx's */
      c_vindex[c_hg->nVtx] = c_hg->nPins;
      vertex = i;
      while (cmatch[vertex] < 0)  {    
        LevelMap[vertex] = c_hg->nVtx;    
        c_hg->vwgt[c_hg->nVtx] += hg->vwgt ? hg->vwgt[vertex] : 1.0;

        for (j = hg->vindex[vertex]; j < hg->vindex[vertex+1]; j++)  {
          if (used_edges [hg->vedge[j]] <= i)  {
            used_edges [hg->vedge[j]]     = i+1;          
            c_vedge[c_hg->nPins++] = hg->vedge[j];
          }              
        }                       
        cmatch[vertex] = -cmatch[vertex] - 1;
        vertex         =  cmatch[vertex];
      }
      c_hg->nVtx++;
    }
  }
  c_vindex[c_hg->nVtx] = c_hg->nPins;
      
  ZOLTAN_FREE ((void**) &used_edges);
  
  /* Update the dist_x array after coarsening: allocate and fill */
  if (!(c_hg->dist_x = (int*)ZOLTAN_MALLOC ((hgc->nProc_x+1) * sizeof(int)))) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
  }  
  size = c_hg->nVtx;
  MPI_Allgather (&size, 1, MPI_INT, each_size, 1, MPI_INT, hgc->row_comm);  
  c_hg->dist_x[0] = 0;
  for (i = 1; i < hgc->nProc_x; i++)
     c_hg->dist_x[i] = c_hg->dist_x[i-1] + each_size[i-1];
  size = 0;
  for (i = 0; i < hgc->nProc_x; i++)
     size += each_size[i];
  c_hg->dist_x[hgc->nProc_x] = size;

  /* Done if there are no remaining vertices */
  if (c_hg->nVtx == 0)  {
    if (!(c_hg->vindex = (int*) ZOLTAN_CALLOC (1, sizeof(int))))  {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT (zz, yo);
      return ZOLTAN_MEMERR;
    };
  }
  else  {
    c_hg->ewgt   = c_ewgt; 
    c_hg->vindex = c_vindex;
    c_hg->vedge  = c_vedge;
    c_hg->nEdge  = hg->nEdge;
    c_hg->comm   = hg->comm;
  }
  
  Zoltan_Multifree (__FILE__, __LINE__, 6, &buffer, &rbuffer, &list, &cmatch,
   &displs, &each_size);  

  ZOLTAN_TRACE_EXIT (zz, yo);
  return Zoltan_HG_Create_Mirror(zz, c_hg);
}



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
