/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/
#ifndef lint
static char *cvs_loadbal = "$Id$";
#endif

/*--------------------------------------------------------------------------*/
/* Purpose: Call Zoltan to determine a new load balance.                    */
/*          Contains all of the callback functions that Zoltan needs        */
/*          for the load balancing.                                         */
/*--------------------------------------------------------------------------*/
/* Author(s):  Matthew M. St.John (9226)                                    */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/* Revision History:                                                        */
/*    10 May 1999:       Date of creation.                                  */
/*--------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>

#include "dr_const.h"
#include "dr_err_const.h"
#include "dr_loadbal_const.h"
#include "dr_eval_const.h"
#include "all_allo_const.h"

/*
 *  PROTOTYPES for load-balancer interface functions.
 */

LB_NUM_OBJ_FN get_num_elements;
/* not used right now --->
LB_OBJ_LIST_FN get_elements;
*/
LB_FIRST_OBJ_FN get_first_element;
LB_NEXT_OBJ_FN get_next_element;

LB_NUM_GEOM_FN get_num_geom;
LB_GEOM_FN get_geom;

LB_NUM_EDGES_FN get_num_edges;
LB_EDGE_LIST_FN get_edge_list;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int run_zoltan(int Proc, PROB_INFO_PTR prob, ELEM_INFO *elements[])
{
/* Local declarations. */
  char *yo = "run_zoltan";
  struct LB_Struct *lb_obj;

  /* Variables returned by the load balancer */
  LB_GID *import_gids = NULL;    /* Global node nums of nodes to be imported */
  LB_LID *import_lids = NULL;    /* Pointers to nodes to be imported         */
  int   *import_procs = NULL;    /* Proc IDs of procs owning nodes to be
                                    imported.                                */
  LB_GID *export_gids = NULL;    /* Global node nums of nodes to be exported */
  LB_LID *export_lids = NULL;    /* Pointers to nodes to be exported         */
  int   *export_procs = NULL;    /* Proc IDs of destination procs for nodes
                                    to be exported.                          */
  int num_imported;              /* Number of nodes to be imported.          */
  int num_exported;              /* Number of nodes to be exported.          */
  int new_decomp;                /* Flag indicating whether the decomposition
                                    has changed                              */
  int i;                         /* Loop index                               */

/***************************** BEGIN EXECUTION ******************************/

  DEBUG_TRACE_START(Proc, yo);

  LB_Set_Param(NULL, "DEBUG_MEMORY", "1");

  /*
   *  Create a load-balancing object.
   */
  if ((lb_obj = LB_Create_Object(MPI_COMM_WORLD)) == NULL) {
    Gen_Error(0, "fatal:  NULL object returned from LB_Create_Object()\n");
    return 0;
  }

  /* set the user-specified parameters */
  for (i = 0; i < prob->num_params; i++) {
    LB_Set_Param(lb_obj, prob->params[i][0], prob->params[i][1]);
  }

  /* set the method */
  if (LB_Set_Method(lb_obj, prob->method) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Method()\n");
    return 0;
  }

  /*
   * set the callback functions
   */

  if (LB_Set_Fn(lb_obj, LB_NUM_OBJ_FN_TYPE, (void *) get_num_elements,
                NULL) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Fn()\n");
    return 0;
  }

  if (LB_Set_Fn(lb_obj, LB_FIRST_OBJ_FN_TYPE, (void *) get_first_element,
                (void *) *elements) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Fn()\n");
    return 0;
  }

  if (LB_Set_Fn(lb_obj, LB_NEXT_OBJ_FN_TYPE, (void *) get_next_element,
                (void *) *elements) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Fn()\n");
    return 0;
  }

  /* functions for geometry based algorithms */
  if (LB_Set_Fn(lb_obj, LB_NUM_GEOM_FN_TYPE, (void *) get_num_geom,
                NULL) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Fn()\n");
    return 0;
  }

  if (LB_Set_Fn(lb_obj, LB_GEOM_FN_TYPE, (void *) get_geom,
                (void *) *elements) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Fn()\n");
    return 0;
  }

  /* functions for geometry based algorithms */
  if (LB_Set_Fn(lb_obj, LB_NUM_EDGES_FN_TYPE, (void *) get_num_edges,
                (void *) *elements) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Fn()\n");
    return 0;
  }

  if (LB_Set_Fn(lb_obj, LB_EDGE_LIST_FN_TYPE, (void *) get_edge_list,
                (void *) *elements) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Set_Fn()\n");
    return 0;
  }

  /* Evaluate the old balance */
  if (Debug_Driver > 0) {
    if (lb_obj->Proc == 0) printf("\nBEFORE load balancing\n");
    driver_eval();
    LB_Eval(lb_obj, 2, 0, 0, NULL, NULL, &i);
  }

  /*
   * call the load balancer
   */
  if (LB_Balance(lb_obj, &new_decomp, &num_imported, &import_gids,
                 &import_lids, &import_procs, &num_exported, &export_gids,
                 &export_lids, &export_procs) == LB_FATAL) {
    Gen_Error(0, "fatal:  error returned from LB_Balance()\n");
    return 0;
  }

  /*
   * call another routine to perform the migration
   */
  if (new_decomp) {
    if (!migrate_elements(Proc, elements, lb_obj, num_imported, import_gids,
                          import_lids, import_procs, num_exported, export_gids,
                          export_lids, export_procs)) {
      Gen_Error(0, "fatal:  error returned from migrate_elements()\n");
      return 0;
    }
  }

  /* Evaluate the new balance */
  if (Debug_Driver > 0) {
    if (lb_obj->Proc == 0) printf("\nAFTER load balancing\n");
    driver_eval();
    LB_Eval(lb_obj, 2, 0, 0, NULL, NULL, &i);
  }

  /* Clean up */
  (void) LB_Free_Data(&import_gids, &import_lids, &import_procs,
                      &export_gids, &export_lids, &export_procs);

  LB_Destroy_Object(&lb_obj);

  LB_Memory_Stats();

  DEBUG_TRACE_END(Proc, yo);
  return 1;

}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_num_elements(void *data, int *ierr)
{
  *ierr = LB_OK; /* set error code */

  return(Mesh.num_elems);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_first_element(void *data, LB_GID *global_id, LB_LID *local_id,
                      int wdim, float *wgt, int *ierr)
{
  ELEM_INFO *elem;

  if (data == NULL) {
    *ierr = LB_FATAL;
    return 0;
  }
  
  elem = (ELEM_INFO *) data;

  *local_id = 0;
  *global_id = elem[*local_id].globalID;

  if (wdim>0)
    *wgt = elem[*local_id].cpu_wgt;

  if (wdim>1)
    *ierr = LB_WARN; /* we didn't expect multidimensional weights */
  else
    *ierr = LB_OK; 

  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_next_element(void *data, LB_GID global_id, LB_LID local_id,
                     LB_GID *next_global_id, LB_LID *next_local_id, 
                     int wdim, float *next_wgt, int *ierr)
{
  int found = 0;
  ELEM_INFO *elem;

  if (data == NULL) {
    *ierr = LB_FATAL;
    return 0;
  }
  
  elem = (ELEM_INFO *) data;

  if (local_id+1 < Mesh.num_elems) { 
    found = 1;
    *next_local_id = local_id + 1;
    *next_global_id = elem[*next_local_id].globalID;

    if (wdim>0)
      *next_wgt = elem[*next_local_id].cpu_wgt;

    if (wdim>1)
      *ierr = LB_WARN; /* we didn't expect multidimensional weights */
    else
      *ierr = LB_OK; 
  }

  return(found);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_num_geom(void *data, int *ierr)
{

  *ierr = LB_OK; /* set error flag */

  return(Mesh.num_dims);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void get_geom(void *data, LB_GID global_id, LB_LID local_id,
              double *coor, int *ierr)
{
  ELEM_INFO *elem;
  int i, j;
  double tmp;

  if (data == NULL) {
    *ierr = LB_FATAL;
    return;
  }

  elem = (ELEM_INFO *) data;

  if (Mesh.eb_nnodes[elem[local_id].elem_blk] == 0) {
    /* No geometry info was read. */
    *ierr = LB_FATAL;
    return;
  }
  
  /*
   * calculate the geometry of the element by averaging
   * the coordinates of the nodes in its connect table
   */
  for (i = 0; i < Mesh.num_dims; i++) {
    tmp = 0.0;
    for (j = 0; j < Mesh.eb_nnodes[elem[local_id].elem_blk]; j++)
      tmp += elem[local_id].coord[j][i];

    coor[i] = tmp / Mesh.eb_nnodes[elem[local_id].elem_blk];
  }

  *ierr = LB_OK;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int get_num_edges(void *data, LB_GID global_id, LB_LID local_id, int *ierr)
{
  ELEM_INFO *elem;

  if (data == NULL) {
    *ierr = LB_FATAL;
    return 0;
  }

  elem = (ELEM_INFO *) data;

  *ierr = LB_OK;

  return(elem[local_id].nadj);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void get_edge_list (void *data, LB_GID global_id, LB_LID local_id,
                   LB_GID *nbor_global_id, int *nbor_procs,
                   int get_ewgts, int *nbor_ewgts, int *ierr)
{
  ELEM_INFO *elem;
  int i, j, proc, local_elem;

  if (data == NULL) {
    *ierr = LB_FATAL;
    return;
  }

  elem = (ELEM_INFO *) data;

  /* get the processor number */
  MPI_Comm_rank(MPI_COMM_WORLD, &proc);

  j = 0;
  for (i = 0; i < elem[local_id].adj_len; i++) {

    /* Skip NULL adjacencies (sides that are not adjacent to another elem). */
    if (elem[local_id].adj[i] == -1) continue;

    if (elem[local_id].adj_proc[i] == proc) {
      local_elem = elem[local_id].adj[i];
      nbor_global_id[j] = elem[local_elem].globalID;
    }
    else { /* adjacent element on another processor */
      nbor_global_id[j] = elem[local_id].adj[i];
    }
    nbor_procs[j] = elem[local_id].adj_proc[i];

    if (get_ewgts) {
      if (elem[local_id].edge_wgt == NULL)
        nbor_ewgts[j] = 1; /* uniform weights is default */
      else
        nbor_ewgts[j] = (int) elem[local_id].edge_wgt[i];
    }
    j++;
  }

  *ierr = LB_OK;
}
