/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef __PHG_HYPERGRAPH_H
#define __PHG_HYPERGRAPH_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


typedef struct {
   int info;    /* primarily for debugging recursive algorithms;initially 0 */
   int nVtx;    /* number of vertices, |V| */
   int nEdge;   /* Size of neigh array; 2|E| */
   int nDim;    /* Number of dimensions for a vertex's coordinate */
   int VertexWeightDim;  /* number of weight dimensions for a vertex */
   int EdgeWeightDim;    /* number of weight dimensions for an edge */
   int redl;             /* Working Reduction limit. */

   int *vtxdist;  /* distributions of vertices to processors, as in ParMETIS.
                     Vertices vtxdist[n] to vtxdist[n+1]-1 are stored on
                     processor n.   KDD:  temporary; may change later. */

   /* physical coordinates of each vertex, optional */
   double *coor; /*  |V| long by CoordinateDim */

   /* arrays with vertex and edge weights */
   float *vwgt;  /* weights of vertices, |V| long by VtxWeightDim */
   float *ewgt;  /* weights of hypergraph edges, 2|E| long by EdgeWeightDim */

   /* arrays to look up the neighbors of a vertex */
   int *nindex;  /* length |V|+1 index to neigh, last is 2|E| */
   int *neigh;   /* length 2|E|, list of neighbors for each vertex */
} PGraph;
  
typedef struct {
  int info;       /* primarily for debugging recursive algorithms;initially 0 */
  int nProc_x;    /* number of processors in x-direction of 2D data distrib.  */
  int nProc_y;    /* number of processors in y-direction of 2D data distrib.  */
                  /* nProc_x * nProc_y should equal number of processors!     */
  int myProc_x;   /* my processor's row block number in [0,nProc_x-1] */
  int myProc_y;   /* my processor's column block number in [0,nProc_y-1] */
  int *dist_x;    /* distributions of vertices to processor columns. Vertices
                   * dist_x[n] to dist_x[n+1]-1 are stored in col block n */
  int *dist_y;    /* distribution of hyperedges to processor rows as above */                  
  int nVtx;             /* number of vertices on this processor */
  int nEdge;            /* number of hyperedges on this processor */
  int nNonZero;         /* number of nonzeros (pins) on this processor */
  
  int VertexWeightDim;  /* number of weight dimensions for a vertex */
  int EdgeWeightDim;    /* number of weight dimensions for a hyperedge */

  int redl;             /* working reduction limit */

  /* physical coordinates of each vertex, optional */
  int nDim;         /* number of coordinate dimensions for a vertex */
  double *coor;     /* |V| long by CoordinateDim */

  /* arrays with vertex and edge weights */
  float *vwgt;      /* weights of vertices, |V| long by VtxWeightDim */
  float *ewgt;      /* weights of hypergraph edges, |E| long by EdgeWeightDim */

  /* arrays to look up vertices given a hyperedge */
  int *hindex;      /* length |E|+1 index into hvertex, last is |P| */
  int *hvertex;     /* length |P| array containing associated vertices */

  /* arrays to look up hyperedges given a vertex */
  int *vindex;      /* length |V|+1 index into vedge, last is |P| */
  int *vedge;       /* length |P| array containing associated hyperedges */
  
  int *vmap;        /* used when recursively dividing for p > 2 */
  double ratio;     /* split when recursively dividing for p > 2 */
} PHGraph;


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* ZOLTAN_HYPERGRAPH_H */
