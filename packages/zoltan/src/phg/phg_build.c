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
#include "phg_hypergraph.h"
#include "zz_const.h"
#include "parmetis_jostle.h"
#include "zz_util_const.h"

#define MEMORY_ERROR { \
  ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Memory error."); \
  ierr = ZOLTAN_MEMERR; \
  goto End; \
}

#define TWOD_GNO_TO_LNO(gno, dist, myblock) (gno) - (dist)[(myblock)]

#define ZOLTAN_PHG_EDGE_TO_ROW(gno,dist_y,nProc_y) \
        Zoltan_PHG_Gno_To_Proc_Block(gno, dist_y, nProc_y);
#define ZOLTAN_PHG_VTX_TO_COL(gno,dist_x,nProc_x) \
        Zoltan_PHG_Gno_To_Proc_Block(gno, dist_x, nProc_x);

int Zoltan_PHG_Gno_To_Proc_Block(
  int gno,
  int *dist_dim,
  int nProc_dim
)
{
/* Function that locates a given global number gno within a distribution 
 * vector dist.
 * Works for both vtx and edges.
 * Takes an initial guess based on equal distribution of gno's.
 * Modifies guess based on actual distribution.
 */

int idx;
int maxgno = dist_dim[nProc_dim];

  idx = gno * nProc_dim / maxgno;

  while (gno < dist_dim[idx]) idx--;
  while (gno >= dist_dim[idx+1]) idx++;

  return idx;
}

/*****************************************************************************/
/* Function prototypes */

static int hash_lookup (ZZ*, ZOLTAN_ID_PTR, int, struct Hash_Node**);
static int Zoltan_PHG_Fill_Hypergraph (ZZ*, ZHG*, PHGPartParams*);

/*****************************************************************************/

int Zoltan_PHG_Build_Hypergraph(
  ZZ *zz,                            /* Zoltan data structure */
  ZHG **zoltan_hg,                   /* Hypergraph to be allocated and built.*/
  PHGPartParams *hgp                 /* Parameters for HG partitioning.*/
)
{
/* allocates and builds hypergraph data structure using callback routines */
ZHG *zhg;                     /* Temporary pointer to Zoltan_PHGraph. */
PHGraph *phgraph;             /* Temporary pointer to HG field */
int ierr = ZOLTAN_OK;
char *yo = "Zoltan_PHG_Build_Hypergraph";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Allocate a Zoltan hypergraph.  */
  zhg = *zoltan_hg = (ZHG*) ZOLTAN_MALLOC (sizeof(ZHG));
  if (zhg == NULL) MEMORY_ERROR;

  /* Initialize the Zoltan hypergraph data fields. */
  zhg->Global_IDs = NULL;
  zhg->Local_IDs = NULL;
  zhg->Parts = NULL;

  phgraph = &(zhg->HG);
  Zoltan_PHG_HGraph_Init(phgraph);

  /* Use callback functions to build the hypergraph. */
  if (zz->Get_Num_HG_Edges && zz->Get_HG_Edge_List && zz->Get_Num_HG_Pins){
    /* 
     * Hypergraph callback functions exist; 
     * call them and build the HG directly.
     */
    ZOLTAN_TRACE_DETAIL(zz, yo, "Using Hypergraph Callbacks.");

printf("%d KDDKDD NVTXS = %d\n", zz->Proc, phgraph->nVtx);
    ierr = Zoltan_PHG_Fill_Hypergraph(zz, zhg, hgp);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building hypergraph");
      goto End;
    }
  }

  else if ((zz->Get_Num_Edges != NULL || zz->Get_Num_Edges_Multi != NULL) &&
           (zz->Get_Edge_List != NULL || zz->Get_Edge_List_Multi != NULL)) {
    /* 
     * Hypergraph callback functions don't exist, but graph functions do;
     * call the graph callback, build a graph, and convert it to a hypergraph. 
     */
    PGraph graph;             /* Temporary graph. */

    ZOLTAN_TRACE_DETAIL(zz, yo, "Using Graph Callbacks.");

    ZOLTAN_PRINT_WARN(zz->Proc, yo, "GRAPH TO HGRAPH CONVERSION MAY NOT BE "
                      "CORRECT IN PARALLEL YET -- KDD KDDKDD");

    Zoltan_PHG_Graph_Init(&graph);
    ierr = Zoltan_Get_Obj_List(zz, &(graph.nVtx), &(zhg->Global_IDs),
     &(zhg->Local_IDs), zz->Obj_Weight_Dim, &(graph.vwgt), &(zhg->Parts));
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error getting object data");
      Zoltan_PHG_Graph_Free(&graph);
      goto End;
    }

    ierr = Zoltan_Build_Graph(zz, 1, hgp->check_graph, graph.nVtx,
     zhg->Global_IDs, zhg->Local_IDs, zz->Obj_Weight_Dim, zz->Edge_Weight_Dim,
     &(graph.vtxdist), &(graph.nindex), &(graph.neigh), &(graph.ewgt));
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error building graph");
      Zoltan_PHG_Graph_Free(&graph);
      goto End;
    }

    graph.nEdge = graph.nindex[graph.nVtx];
    ierr = Zoltan_PHG_Graph_to_HGraph(zz, &graph, phgraph);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error converting graph to hypergraph");
      Zoltan_PHG_Graph_Free(&graph);
      goto End;
    }

    Zoltan_PHG_Graph_Free(&graph);
  }

  if (zz->Get_Num_Geom != NULL && 
      (zz->Get_Geom != NULL || zz->Get_Geom_Multi != NULL)) {
     /* Geometric callbacks are registered;       */
     /* get coordinates for hypergraph objects.   */
     ZOLTAN_TRACE_DETAIL(zz, yo, "Getting Coordinates.");
     ierr = Zoltan_Get_Coordinates(zz, phgraph->nVtx, zhg->Global_IDs,
      zhg->Local_IDs, &(phgraph->nDim), &(phgraph->coor));
  }

  if (hgp->check_graph) {
    ierr = Zoltan_PHG_Check(zz, phgraph);
    if (ierr == ZOLTAN_WARN) {
      ZOLTAN_PRINT_WARN(zz->Proc, yo, "Warning returned from Zoltan_PHG_Check");
    }
    else if (ierr != ZOLTAN_OK) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Zoltan_PHG_Check");
      goto End;     
    }
  }

  if (hgp->output_level >= PHG_DEBUG_PRINT)
    Zoltan_PHG_HGraph_Print(zz, zhg, &(zhg->HG), stdout);

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    /* Return NULL zhg */
    Zoltan_PHG_HGraph_Free(&(zhg->HG));
    Zoltan_Multifree(__FILE__, __LINE__, 4, &(zhg->Global_IDs),
     &(zhg->Local_IDs), &(zhg->Parts), zoltan_hg);
  }
    
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*****************************************************************************/


static int Zoltan_PHG_Fill_Hypergraph(
  ZZ *zz,
  ZHG *zhg,     /* Description of hypergraph provided by the application. */
  PHGPartParams *hgp
)
{
/* Routine to call HG query function and build HG data structure. 
 * Assumes (for now) that input is given as follows from application:
 *  - Application already provided list of objects (vertices) assigned
 *    to proc.  (KDD -- perhaps move that call here).
 *  - Application gives hyperedges it owns; it knows GIDs and processor owner
 *    of all vertices of each owned hyperedge.
 *  - Application gives each hyperedge only once (KDDKDD -- MAY REMOVE THIS
 *    RESTRICTION LATER. KDDKDD)
 * Output is a fully functioning parallel hypergraph with 2D distribution of
 * pins (non-zeros).
 */

char *yo = "Zoltan_PHG_Fill_Hypergraph";
struct application_input {     /* Data provided by application callbacks. */
  int nVtx;                         /* # objects (vertices) on proc. */
  int nEdge;                        /* # hyperedges on proc. */
  int nPins;                        /* # pins (nonzeros) on proc. */
  int GnVtx;                        /* Total nVtx across all procs. */
  int GnEdge;                       /* Total nEdge across all procs. */
  int *edge_sizes;                  /* # of GIDs in each hyperedge. */
  ZOLTAN_ID_PTR pins;               /* Object GIDs (vertices) belonging to 
                                       hyperedges.  */
  int *pin_procs;                   /* Processor owning each pin vertex of 
                                       hyperedges. */
  int *pin_gno;                     /* Global numbers in range [0,GnVtx-1]
                                       for pins. */
  int *vtxdist;                     /* Distribution of vertices
                                       across original owning processors;
                                       # vtxs on proc i == 
                                       vtxdist[i+1] - vtxdist[i]. */
  int *edgedist;                    /* Distribution of edges
                                       across original owning processors;
                                       # edges on proc i == 
                                       edgedist[i+1] - edgedist[i]. */
  float *vwgt;                      /* Vertex weights. */
  float *ewgt;                      /* Edge weights. */
} app;

struct Hash_Node *hash_nodes = NULL;  /* Hash table variables for mapping   */
struct Hash_Node **hash_tab = NULL;   /* GIDs to global numbering system.   */
ZOLTAN_COMM_OBJ *plan;

int i, j, cnt;
int msg_tag = 30000;
int numwgts = 0;
int ierr = ZOLTAN_OK;
int nProc = zz->Num_Proc;
int nRequests;
ZOLTAN_ID_PTR pin_requests;
int *request_gno = NULL;
int edge_gno, edge_Proc_y;
int vtx_gno, vtx_Proc_x;
int *proclist = NULL;
int *sendbuf = NULL;
int nnz, idx;
int *nonzeros = NULL;
int *tmp = NULL;
int *hindex = NULL, *hvertex = NULL;
int *dist_x = NULL, *dist_y = NULL;
int nEdge, nVtx;

ZOLTAN_ID_PTR global_ids = zhg->Global_IDs;  
PHGraph *phg = &(zhg->HG);
int num_gid_entries = zz->Num_GID;
int nProc_x = hgp->nProc_x;
int nProc_y = hgp->nProc_y;
int myProc_x = hgp->myProc_x;
int myProc_y = hgp->myProc_y;
float frac_x, frac_y;

  ZOLTAN_TRACE_ENTER(zz, yo);

  /**************************************************/
  /* Obtain vertex information from the application */
  /**************************************************/

  ierr = Zoltan_Get_Obj_List(zz, &(zhg->nObj), &(zhg->Global_IDs),
                             &(zhg->Local_IDs), 
                             zz->Obj_Weight_Dim, &app.vwgt, 
                             &(zhg->Parts));
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error getting object data");
    goto End;
  }

  app.nVtx = zhg->nObj;
  app.edge_sizes = NULL;
  app.pins = NULL;
  app.pin_procs = NULL;
  app.pin_gno = NULL;
  app.vtxdist = NULL;
  app.edgedist = NULL;
  app.vwgt = NULL;
  app.ewgt = NULL;

  /* Build app.vtxdist as in Zoltan_Build_Graph. */

  app.vtxdist = (int *)ZOLTAN_MALLOC(2 * (nProc+1) * sizeof(int));
  if (!(app.vtxdist)) MEMORY_ERROR;

  app.edgedist = app.vtxdist + nProc + 1;

  /* Construct app.vtxdist[i] = the number of vertices on all procs < i. */
  /* Scan to compute partial sums of the number of objs */

  MPI_Scan (&app.nVtx, app.vtxdist, 1, MPI_INT, MPI_SUM, zz->Communicator);

  /* Gather data from all procs */

  MPI_Allgather (&(app.vtxdist[0]), 1, MPI_INT,
                 &(app.vtxdist[1]), 1, MPI_INT, zz->Communicator);
  app.vtxdist[0] = 0;
  app.GnVtx = app.vtxdist[nProc];

  /* 
   * Correlate GIDs in edge_verts with local indexing in zhg to build the
   * input HG.
   * Use hash table to map global IDs to local position in zhg->Global_IDs.
   * Based on hashing code in Zoltan_Build_Graph.
   * KDD -- This approach is serial for now; look more closely at 
   * KDD -- Zoltan_Build_Graph when move to parallel.
   */

  /* Construct local hash table mapping GIDs to global number (gno) */
  if (app.nVtx) {
    hash_nodes = (struct Hash_Node *) ZOLTAN_MALLOC(app.nVtx * 
                                                    sizeof(struct Hash_Node));
    hash_tab = (struct Hash_Node **) ZOLTAN_MALLOC(app.nVtx *
                                                   sizeof(struct Hash_Node *));
    if (!hash_nodes || !hash_tab) MEMORY_ERROR;

    /* Assign consecutive numbers based on the order of the ids */
    for (i=0; i< app.nVtx; i++) {
      hash_tab[i] = NULL;
      hash_nodes[i].gid = &(global_ids[i*num_gid_entries]);
      hash_nodes[i].gno = app.vtxdist[zz->Proc]+i;
    }

    for (i=0; i< app.nVtx; i++){
      /* insert hashed elements into hash table */
      j = Zoltan_Hash(&(global_ids[i*num_gid_entries]), num_gid_entries,
                      (unsigned int) app.nVtx);
      hash_nodes[i].next = hash_tab[j];
      hash_tab[j] = &hash_nodes[i];
    }
  }

  /***********************************************************************/
  /* Get hyperedge information from application through query functions. */
  /***********************************************************************/

  app.nEdge = zz->Get_Num_HG_Edges(zz->Get_Num_HG_Edges_Data, &ierr);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_Num_HG_Edges");
    goto End;
  }
  
  /* KDD:  question:  How do we compute size to malloc array for HG Edges? 
   * KDD:  We can't have a size function unless we assume the application
   * KDD:  can "name" the hyperedges.
   * KDD:  For now, assume application can return number of pins.
   */

  app.nPins = zz->Get_Num_HG_Pins(zz->Get_Num_HG_Pins_Data, &ierr);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "Error returned from Get_Max_HG_Edge_Size");
    goto End;
  }

  if (app.nEdge > 0) {
    numwgts = app.nEdge * zz->Edge_Weight_Dim;
    app.pins = ZOLTAN_MALLOC_GID_ARRAY(zz, app.nPins);
    app.edge_sizes = (int *) ZOLTAN_MALLOC(app.nEdge * sizeof(int));
    app.pin_procs = (int *) ZOLTAN_MALLOC(app.nPins * sizeof(int));
    app.pin_gno = (int *) ZOLTAN_MALLOC(app.nPins * sizeof(int));
    if (numwgts) 
      app.ewgt = (float *) ZOLTAN_MALLOC(numwgts * sizeof(float));
    if (!app.pins || !app.edge_sizes || !app.pin_procs || !app.pin_gno ||
        (numwgts && !app.ewgt)) MEMORY_ERROR;

    ierr = zz->Get_HG_Edge_List(zz->Get_HG_Edge_List_Data, num_gid_entries,
                                zz->Edge_Weight_Dim, app.nEdge, app.nPins,
                                app.edge_sizes, app.pins, 
                                app.pin_procs, app.ewgt);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,"Error returned from Get_HG_Edge_List");
      goto End;
    }
  }
  
  /* 
   * KDDKDD -- Assuming hyperedges are given to Zoltan by one processor only.
   * KDDKDD -- Eventually, will loosen that constraint and remove duplicates.
   * KDDKDD -- Or the code might work (although with extra communication)
   * KDDKDD -- with the duplicates.
   * KDDKDD -- Might be easier once we have edge GIDs.
   */
  /* Impose a global hyperedge numbering */
  /* Construct app.edgedist[i] = the number of edges on all procs < i. */
  /* Scan to compute partial sums of the number of edges */

  MPI_Scan (&app.nEdge, app.edgedist, 1, MPI_INT, MPI_SUM, zz->Communicator);

  /* Gather data from all procs */

  MPI_Allgather (&(app.edgedist[0]), 1, MPI_INT,
                 &(app.edgedist[1]), 1, MPI_INT, zz->Communicator);
  app.edgedist[0] = 0;
  app.GnEdge = app.edgedist[nProc];

  /* 
   * Obtain the global num in range 0 .. (total_num_vtx-1) 
   * for each vertex pin.
   * For each pin GID in app.pins, request the global number (gno) from the
   * input processor.
   * Fill requests (using hash table) for GIDs local to this processor.
   * Upon completion, app.pin_gno will contain the global nums.
   */

  ierr = Zoltan_Comm_Create(&plan, app.nPins, app.pin_procs, zz->Communicator,
                            msg_tag, &nRequests);

  if (nRequests) {
    pin_requests = ZOLTAN_MALLOC_GID_ARRAY(zz, nRequests);
    request_gno = (int *) ZOLTAN_MALLOC(nRequests * sizeof(int));
    if (!pin_requests || !request_gno)
      MEMORY_ERROR;
  }

  msg_tag--;
  ierr = Zoltan_Comm_Do(plan, msg_tag, (char *) app.pins, 
                        sizeof(ZOLTAN_ID_TYPE) * num_gid_entries,
                        (char *) pin_requests);

  for (i = 0; i < nRequests; i++)
    request_gno[i] = hash_lookup(zz, &(pin_requests[i*num_gid_entries]),
                                 app.nVtx, hash_tab);

  msg_tag--;
  Zoltan_Comm_Do_Reverse(plan, msg_tag, (char *) request_gno, sizeof(int), NULL,
                         (char *) app.pin_gno);

  Zoltan_Comm_Destroy(&plan);

  ZOLTAN_FREE(&pin_requests);
  ZOLTAN_FREE(&request_gno);

  /* 
   * Compute the distribution of vertices and edges to the 2D data
   * distribution's processor columns and rows.
   * For now, these distributions are described by arrays dist_x
   * and dist_y; in the future, we may prefer a hashing function
   * mapping GIDs to processor columns and rows. KDDKDD
   */

  dist_x = (int *) ZOLTAN_CALLOC((nProc_x+1), sizeof(int));
  dist_y = (int *) ZOLTAN_CALLOC((nProc_y+1), sizeof(int));

  if (!dist_x || !dist_y) MEMORY_ERROR;

  frac_x = (float) app.GnVtx / (float) nProc_x;
  for (i = 1; i < nProc_x; i++)
    dist_x[i] = (int) (i * frac_x);
  dist_x[nProc_x] = app.GnVtx;
  
  frac_y = (float) app.GnEdge / (float) nProc_y;
  for (i = 1; i < nProc_y; i++)
    dist_y[i] = (int) (i * frac_y);
  dist_y[nProc_y] = app.GnEdge;
  
  nEdge = dist_y[myProc_y+1] - dist_y[myProc_y];
  nVtx = dist_x[myProc_x+1] - dist_x[myProc_x];

  /*
   * Build comm plan for sending non-zeros to their target processors in
   * 2D data distribution. 
   */

  cnt = 0; 
  for (i = 0; i < app.nEdge; i++) {
    /* processor row for the edge */
    edge_gno = app.edgedist[zz->Proc] + i;
    edge_Proc_y = ZOLTAN_PHG_EDGE_TO_ROW(edge_gno, dist_y, myProc_y);

    for (j = 0; j < app.edge_sizes[i]; j++) {
      /* processor column for the vertex */
      vtx_gno = app.pin_gno[cnt];
      vtx_Proc_x = ZOLTAN_PHG_VTX_TO_COL(vtx_gno, dist_x, myProc_x);

      proclist[cnt] = edge_Proc_y * nProc_x + vtx_Proc_x;
      sendbuf[2*cnt] = edge_gno;
      sendbuf[2*cnt+1] = vtx_gno;
      cnt++;
    } 
  }

  Zoltan_Multifree(__FILE__, __LINE__, 8, &app.pins, 
                                          &app.edge_sizes, 
                                          &app.pin_procs, 
                                          &app.pin_gno, 
                                          &app.vwgt,
                                          &app.ewgt,
                                          &hash_nodes,
                                          &hash_tab);
  /*
   * Send pins to their target processors.
   * They become non-zeros in the 2D data distribution.
   */

  msg_tag--;
  ierr = Zoltan_Comm_Create(&plan, cnt, proclist, zz->Communicator,
                     msg_tag, &nnz);

  msg_tag--;
  Zoltan_Comm_Do(plan, msg_tag, (char *) sendbuf, 2*sizeof(int),
                 (char *) nonzeros);

  Zoltan_Comm_Destroy(&plan);

  ZOLTAN_FREE(&proclist);
  ZOLTAN_FREE(&sendbuf);

  /* Unpack the non-zeros received. */

  tmp = (int *) ZOLTAN_CALLOC(nEdge + 1, sizeof(int));
  hindex = (int *) ZOLTAN_CALLOC(nEdge + 1, sizeof(int));
  if (!tmp || !hindex) MEMORY_ERROR;

  /* Count the number of nonzeros per hyperedge */
  for (i = 0; i < nnz; i++) {
    idx = TWOD_GNO_TO_LNO(nonzeros[2*i], dist_y, myProc_y); 
    tmp[idx]++;
  }

  /* Compute prefix sum to represent hindex correctly. */
  for (i = 0; i < nEdge; i++)  {
    hindex[i+1] = hindex[i] + tmp[i];
    tmp[i] = 0;
  }
       
  for (i = 0; i < nnz; i++) {
    idx = TWOD_GNO_TO_LNO(nonzeros[2*i], dist_y, myProc_y); 
    hvertex[hindex[idx] + tmp[idx]] = 
            TWOD_GNO_TO_LNO(nonzeros[2*i+1], dist_x, myProc_x);
    tmp[idx]++;
  }

  ZOLTAN_FREE(&tmp);

  phg->nVtx = nVtx;
  phg->nEdge = nEdge;
  phg->nNonZero = nnz;
  phg->hindex = hindex;
  phg->hvertex = hvertex;
  phg->dist_x = dist_x;
  phg->dist_y = dist_y;

  /*
   * KDDKDD -- Not yet handling vertex and edge weights.
   * KDDKDD -- What should be done with them?
   */

  ierr = Zoltan_PHG_Create_Mirror(zz, phg);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error from Zoltan_PHG_Create_Mirror");
    goto End;
  }

End:
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    Zoltan_PHG_HGraph_Free(phg);
  }
  
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/*****************************************************************************/

static int hash_lookup(
  ZZ *zz,
  ZOLTAN_ID_PTR key,
  int nVtx,
  struct Hash_Node **hash_tab
)
{
/* Looks up a key GID in the hash table; returns its gno. */
/* Based on hash_lookup in build_graph.c. */

  int i;
  struct Hash_Node *ptr;

  i = Zoltan_Hash(key, zz->Num_GID, (unsigned int) nVtx);
  for (ptr=hash_tab[i]; ptr != NULL; ptr = ptr->next){
    if (ZOLTAN_EQ_GID(zz, ptr->gid, key))
      return (ptr->gno);
  }
  /* Key not in hash table */
  return -1;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
