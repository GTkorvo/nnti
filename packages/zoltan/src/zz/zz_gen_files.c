/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002,2003, Sandia National Laboratories.          *
 * For more info, see the README file in the top-level Zoltan directory.     *   *****************************************************************************/
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


#include "zz_const.h"
#include "parmetis_jostle.h"

#define ZOLTAN_PRINT_VTX_NUM  0  /* print vertex number at beginning of line? */


/* Temporary prototypes. These functions are HG routines
   currently not compiled into Zoltan, but duplicated in this file. */
static int Zoltan_HG_Get_Hedges(ZZ *zz, int **p_hindex, 
           ZOLTAN_ID_PTR *p_edge_verts, int **p_edge_procs, 
           float **p_edge_wgts, int *glob_hedges, int *glob_pins);
static int Zoltan_HG_Print_Hedges(ZZ *zz, FILE *fp, 
           int *hindex, ZOLTAN_ID_PTR hevtxs, float *hewgts);

/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines for generating output files
 *  from Zoltan that describe the data given to Zoltan. These 
 *  functions are all callable by the application, and can be
 *  invoked independent of load-balancing (Zoltan_LB_Balance).
 */
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Generate_Files(ZZ *zz, char *fname, int base_index,
int gen_geom, int gen_graph, int gen_hg)
{
/*
 *  Generate up to four output files:
 *   a) Current assignment of objects to partitions (and procs?)
 *   b) Geometry of the objects, if geometry query functions are registered.
 *   c) Graph if graph query functions are available.
 *   d) Hypergraph if hypergraph query functions are available.
 *
 *  Input parameters:
 *   zz,         pointer to Zoltan struct
 *   fname,      the basename for the output files
 *   base_index, first vertex (object) is labelled 0 or 1?
 *   gen_geom,   generate geometry file? 
 *   gen_graph,  generate graph file? 
 *   gen_hg,     generate hypergraph file? 
 *
 *  The output is serialized, such that each processor
 *  will open and close each output file in turn.
 */

  int error=ZOLTAN_OK;
  ZOLTAN_ID_PTR local_ids = NULL;
  ZOLTAN_ID_PTR global_ids = NULL;
  FILE *fp;
  char full_fname[256];
  int *vtxdist, *xadj, *adjncy, *part;
  int *heprocs, *hindex;
  ZOLTAN_ID_PTR hevtxs;
  float *float_vwgt, *ewgts, *hewgts;
  double *xyz;
  int i, j, k, num_obj, num_geom, num_edges;
  int glob_nvtxs, glob_edges, glob_hedges, glob_pins;
  int print_vtx_num = ZOLTAN_PRINT_VTX_NUM;
  char *yo = "Zoltan_Generate_Files";

  ZOLTAN_TRACE_ENTER(zz, yo);

  /* Initialize all local pointers to NULL. This is necessary
   * because we free all non-NULL pointers upon errors.
   */
  vtxdist = xadj = adjncy = part = NULL;
  float_vwgt = ewgts = hewgts = NULL;
  xyz = NULL;
  heprocs = hindex = NULL;
  hevtxs = NULL;

  /* Assign default file name if none was given. */
  if (fname==NULL) fname = "noname";

  /* Zoltan_Get_Obj_List allocates memory for all return lists. */
  error = Zoltan_Get_Obj_List(zz, &num_obj, &global_ids, &local_ids,
                              zz->Obj_Weight_Dim, &float_vwgt, &part);
  if (error != ZOLTAN_OK){
    /* Return error */
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Get_Obj_List returned error.");
    error = ZOLTAN_FATAL;
    goto End;
  }

  if (gen_graph){
    /* Build (ParMetis) graph data structures. */
    error = Zoltan_Build_Graph(zz, 1, 1, num_obj,
           global_ids, local_ids, zz->Obj_Weight_Dim, zz->Edge_Weight_Dim,
           &vtxdist, &xadj, &adjncy, &ewgts);
    if (error != ZOLTAN_OK && error != ZOLTAN_WARN){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Zoltan_Build_Graph returned error.");
      goto End;
    }
    glob_nvtxs = vtxdist[zz->Num_Proc];
  }
  else{
    /* Compute global number of vertices. */
    MPI_Reduce(&num_obj, &glob_nvtxs, 1, MPI_INT, MPI_SUM, 0, 
        zz->Communicator);  
  }

  /* Local number of edges. */
  if (xadj==NULL || xadj[num_obj]==0)
    num_edges = 0;
  else
    num_edges = xadj[num_obj];
  /* Compute global number of edges. */
  MPI_Reduce(&num_edges, &glob_edges, 1, MPI_INT, MPI_SUM, 0, 
      zz->Communicator);  
  /* Assume no self-edges! */
  glob_edges /= 2;

  /* Build hypergraph, or get hypergraph data. */
  if (gen_hg){
    if (zz->Get_Num_HG_Edges == NULL || zz->Get_HG_Edge_List == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Hypergraph output requested, but no corresponding query function was found.\n");
      error = ZOLTAN_FATAL;
      goto End;
    }
    /* error = Zoltan_HG_Build_Hypergraph(zz, &zhg, NULL); */
    /* Get data in parallel. Zoltan_HG_Build_Hypergraph
       currently only works in serial. */
    error = Zoltan_HG_Get_Hedges(zz, &hindex, &hevtxs, &heprocs, &hewgts,
            &glob_hedges, &glob_pins);
  }

  /**********************************************************/
  /* Write to files, serialized.                            */
  /* Note: This will be slow (not scalable) for many procs. */
  /**********************************************************/
  Zoltan_Print_Sync_Start(zz->Communicator, 0); 

  /* Write object assignments to file. */
  /* For now, only write partition number. */
  if (num_obj > 0){
    sprintf(full_fname, "%s.assign", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_FATAL;
      goto End;
    }
    for (i=0; i<num_obj; i++)
      fprintf(fp, "%d\n", part[i]);
      /* fprintf(fp, "%d\n", zz->Proc); */
    fclose(fp);
  }

  /* Write geometry to file, if applicable. */
  if (gen_geom){
    if (zz->Get_Num_Geom == NULL ||
     (zz->Get_Geom == NULL && zz->Get_Geom_Multi == NULL)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Geometry output requested, but no corresponding query function was found.\n");
      error = ZOLTAN_FATAL;
      goto End;
    }
    error = Zoltan_Get_Coordinates(zz, num_obj, global_ids, local_ids,
                                   &num_geom, &xyz);

    sprintf(full_fname, "%s.coords", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_FATAL;
      goto End;
    }
    for (i=0; i<num_obj; i++){
      for (j=0; j<num_geom; j++)
        fprintf(fp, "%f ", xyz[i*num_geom + j]);
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  /* Write graph to file, if applicable. */
  /* Also create a minimal .graph file with no edges for geometric methods. */
  if (gen_geom || gen_graph){
    sprintf(full_fname, "%s.graph", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_FATAL;
      goto End;
    }

    /* If proc 0, write first line. */
    if (zz->Proc == 0){
      fprintf(fp, "%% First line: #vertices #edges weight_flag\n");
      fprintf(fp, "%d %d %1d%1d%1d", glob_nvtxs, glob_edges, 
        print_vtx_num, (zz->Obj_Weight_Dim>0), (zz->Edge_Weight_Dim>0));
      if (zz->Obj_Weight_Dim>1 || zz->Edge_Weight_Dim>1)
        fprintf(fp, " %d %d", zz->Obj_Weight_Dim, zz->Edge_Weight_Dim);
      fprintf(fp, "\n");
    }


    /* Print edge list for each node (object). */
    for (i=0; i<num_obj; i++){
      /* Print vertex number at beginning of line? */
      if (print_vtx_num){
        fprintf(fp, "%d ", vtxdist[zz->Proc]+base_index+i);
      }
      /* First print object (vertex) weight, if any. */
      for (k=0; k<zz->Obj_Weight_Dim; k++)
        fprintf(fp, "%f ", float_vwgt[i*(zz->Obj_Weight_Dim)+k]);
      if (gen_graph){
        /* If graph, then print neighbor list */
        for (j=xadj[i]; j<xadj[i+1]; j++){
          fprintf(fp, "%d ", adjncy[j]+base_index);
          /* Also print edge weight, if any. */
          for (k=0; k<zz->Edge_Weight_Dim; k++)
            fprintf(fp, "%f ", ewgts[j*(zz->Edge_Weight_Dim)+k]);
        }
      }
      fprintf(fp, "\n");
    }
    fclose(fp);
  }

  Zoltan_Print_Sync_End(zz->Communicator, 0); 

  /* Separate synchronization for hypergraphs; this could be merged
     into the previous synchronization. */

  Zoltan_Print_Sync_Start(zz->Communicator, 0); 

  /* Write hypergraph to file, if applicable. */
  if (gen_hg){
    if (zz->Get_Num_HG_Edges == NULL || zz->Get_HG_Edge_List == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Hypergraph output requested, but no corresponding query function was found.\n");
      error = ZOLTAN_FATAL;
      goto End;
    }
    sprintf(full_fname, "%s.hg", fname);
    if (zz->Proc == 0)
      fp = fopen(full_fname, "w");
    else
      fp = fopen(full_fname, "a");
    if (fp==NULL){
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Could not open file for writing.\n");
      error = ZOLTAN_FATAL;
      goto End;
    }

    /* If proc 0, write first line. */
    if (zz->Proc == 0){
      fprintf(fp, "%% First line: #vertices #hyperedges #pins weight_flag\n");
      fprintf(fp, "%d %d %d %1d%1d%1d", glob_nvtxs, glob_hedges, 
        glob_pins, ZOLTAN_PRINT_VTX_NUM, 
        (zz->Obj_Weight_Dim>0), (zz->Edge_Weight_Dim>0));
      if (zz->Obj_Weight_Dim>1 || zz->Edge_Weight_Dim>1)
        fprintf(fp, " %d %d", zz->Obj_Weight_Dim, zz->Edge_Weight_Dim);
      fprintf(fp, "\n");
    }

    /* Each proc prints its part of the hgraph. */
    Zoltan_HG_Print_Hedges(zz, fp, hindex, hevtxs, hewgts);

    fclose(fp);
  }
  
  Zoltan_Print_Sync_End(zz->Communicator, 0); 

End:
  ZOLTAN_FREE(&xyz);
  ZOLTAN_FREE(&vtxdist);
  ZOLTAN_FREE(&xadj);
  ZOLTAN_FREE(&adjncy);
  ZOLTAN_FREE(&global_ids);
  ZOLTAN_FREE(&local_ids);
  ZOLTAN_FREE(&float_vwgt);
  ZOLTAN_FREE(&ewgts);
  ZOLTAN_FREE(&part);
  if ( zz->Get_Num_HG_Edges != NULL && zz->Get_HG_Edge_List != NULL) {
    ZOLTAN_FREE(&hindex);
    ZOLTAN_FREE(&hevtxs);
    ZOLTAN_FREE(&heprocs);
    ZOLTAN_FREE(&hewgts);
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return error;
}

/************************************************************************ 
 * EBEB: The routines below were copied from the hg directory because   *
 * Zoltan is distributed without hg. Duplicate functions should be      *
 * removed later.                                                       *
 ************************************************************************/

/* Each proc gets hypergraph data from query functions.
 * Compute some global values. 
 * All procs must participate due to collective communication.
 */
static int Zoltan_HG_Get_Hedges(ZZ *zz, int **p_hindex, 
           ZOLTAN_ID_PTR *p_edge_verts, int **p_edge_procs, 
           float **p_edge_wgts, int *glob_hedges, int *glob_pins)
{
  int i, ierr, cnt, j, nEdge, npins, numwgts, minproc;
  int loc_hedges, loc_pins;
  int *hindex = NULL, *edge_procs = NULL;
  ZOLTAN_ID_PTR edge_verts = NULL;
  ZOLTAN_ID_PTR edge_gids = NULL, edge_lids = NULL;
  float *edge_wgts = NULL;
  static char *yo = "Zoltan_HG_Get_Hedges";

  ZOLTAN_TRACE_ENTER(zz, yo);
  ierr = ZOLTAN_OK;

  /* Get hyperedge information from application through query functions. */

  nEdge = zz->Get_Num_HG_Edges(zz->Get_Num_HG_Edges_Data, &ierr);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from Get_Num_HG_Edges");
    goto End;
  }
  if (nEdge > 0) {
    edge_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, nEdge);
    edge_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, nEdge);
    (*p_hindex) = (int *) ZOLTAN_MALLOC((nEdge+1) * sizeof(int));
    hindex = *p_hindex;
    numwgts = nEdge * zz->Edge_Weight_Dim;
    if (numwgts){
       (*p_edge_wgts) = (float *) ZOLTAN_MALLOC(numwgts * sizeof(float));
       edge_wgts = *p_edge_wgts;
    }
    if (!edge_gids || (zz->Num_LID && !edge_lids) || 
        !hindex || (numwgts && !edge_wgts)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    ierr = zz->Get_HG_Edge_Info(zz->Get_HG_Edge_Info_Data, 
                                zz->Num_GID, zz->Num_LID, nEdge,
                                zz->Edge_Weight_Dim,
                                edge_gids, edge_lids, 
                                hindex, edge_wgts);
    npins = 0;
    for (i = 0; i < nEdge; i++) npins += hindex[i];

    (*p_edge_verts) = ZOLTAN_MALLOC_GID_ARRAY(zz, npins);
    edge_verts = *p_edge_verts;
    (*p_edge_procs) = (int *) ZOLTAN_MALLOC(npins * sizeof(int));
    edge_procs = *p_edge_procs;
    if (npins && (!edge_verts || !edge_procs)) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient Memory");
      ierr = ZOLTAN_MEMERR;
      goto End;
    }
    ierr = zz->Get_HG_Edge_List(zz->Get_HG_Edge_List_Data, 
                                zz->Num_GID, zz->Num_LID, 
                                nEdge, edge_gids, edge_lids,
                                hindex, edge_verts, edge_procs);
    if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo,"Error returned from Get_HG_Edge_List");
      goto End;
    }
  }
  /* Build hindex from edge_sizes */
  cnt = 0;
  for (i = 0; i < nEdge; i++) {
    j = hindex[i];
    hindex[i] = cnt;
    cnt += j;
  }
  hindex[nEdge] = cnt;

  /* Compute local hg statistics. */
  /* Make hindex negative if hedge is owned by another proc. */
  loc_hedges = loc_pins = 0;
  for (i = 0; i < nEdge; i++) {
    minproc = zz->Num_Proc;
    for (j=hindex[i]; j<hindex[i+1]; j++)
      if (edge_procs[j]<minproc)
        minproc = edge_procs[j];
    if (minproc == zz->Proc){  /* my hyperedge */
      loc_hedges++;
      loc_pins += (hindex[i+1] - hindex[i]); /* edge_size[i] */
    }
    else  /* lowest proc owns hyperedge, not me */
      hindex[i] = -hindex[i];
  }

  /* Sanity check */
  if (cnt != npins) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, 
                       "Input error:  Number of pins != sum of edge sizes");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  /* Compute global #pins and #hyperedges, no duplicates */
  MPI_Allreduce(&loc_hedges, glob_hedges, 1, MPI_INT, MPI_SUM,
      zz->Communicator);
  MPI_Allreduce(&loc_pins, glob_pins, 1, MPI_INT, MPI_SUM,
      zz->Communicator);

End:
  /* Memory will be freed in calling function. */

  ZOLTAN_TRACE_EXIT(zz, yo);

  return ierr;
}

#define ABS(x) ((x)<0 ? -(x) : (x))

/* Each processor prints its hyperedges to file. */
static int Zoltan_HG_Print_Hedges(ZZ *zz, FILE *fp, 
           int *hindex, ZOLTAN_ID_PTR hevtxs, float *hewgts)
{
  int i, j, ierr, num_edges;
  char *yo = "Zoltan_HG_Print_Hedges";

  ZOLTAN_TRACE_ENTER(zz, yo);

  ierr = ZOLTAN_OK;
  num_edges = zz->Get_Num_HG_Edges(zz->Get_Num_HG_Edges_Data, &ierr);

  for (i=0; i<num_edges; i++){
    if (hindex[i]>=0){
      /* Only print hyperedges owned by me (to avoid duplicate hedges) */
      for (j=hindex[i]; j<ABS(hindex[i+1]); j++){ 
        /* EBEB - Print global ids as integers. */
        fprintf(fp, "%d ", (int) hevtxs[j]);
      }
    }
    fprintf(fp, "\n");
  }

  /* Print weights. EBEB - Not yet impl. */

  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}


#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

