/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <stdio.h>
#include "lb_const.h"
#include "reftree_const.h"
#include "params_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Prototypes for functions internal to this file */

static void LB_Reftree_Free_Subtree(LB *lb, LB_REFTREE *subroot);
static int order_tri_bisect(LB *lb, int *vert1, int *order, int *vertices,
                     int *in_vertex, int *out_vertex, LB_REFTREE *subroot);
static int order_quad_quad(LB *lb, int *vert1, int *order, int *vertices,
                     int *in_vertex, int *out_vertex, LB_REFTREE *subroot);
static int order_other_ref(LB *lb, LB_REFTREE *parent, int num_child, 
                    int *num_vert,
                    int *vert1, int *vertices, int *order, int *in_vertex,
                    int *out_vertex);
static void order_other_ref_recur(int new_entry, int level, int *order, 
                          int *on_path,
                          int num_child, int *has_out, int **share_vert,
                          int max_share, int *solved);
static int find_inout(int level, int num_child, int *num_vert, int *vert1,
               int *vertices, int *in_vertex, int *out_vertex, int *order);
static int LB_Reftree_Reinit_Coarse(LB *lb);
static int LB_Reftree_Build_Recursive(LB *lb,LB_REFTREE *subroot);
static int alloc_reftree_nodes(LB *lb, LB_REFTREE **node, int num_node,
                               int *num_vert);
void free_reftree_nodes(LB_REFTREE **node);

static ZOLTAN_ID_PTR slocal_gids;  /* coarse element Global IDs from user */
static ZOLTAN_ID_PTR slocal_lids;  /* coarse element Local IDs from user */
static int *sassigned;         /* 1 if the element is assigned to this proc */
static int *snum_vert;         /* number of vertices for each coarse element */
static int *svertices;         /* vertices for the coarse elements */
static int *sin_vertex;        /* "in" vertex for each coarse element */
static int *sout_vertex;       /* "out" vertex for each coarse element */
static int *svert1;        /* array containing the first vertex for each child*/
static int *sorder;            /* order of the children */
static int ssize;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*  Parameters structure for reftree methods */
static PARAM_VARS REFTREE_params[] = {
        { "REFTREE_HASH_SIZE", NULL, "INT" },
        { NULL, NULL, NULL } };

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Set_Reftree_Param(
char *name,                     /* name of variable */
char *val)                      /* value of variable */
{
    int status;
    PARAM_UTYPE result;         /* value returned from Check_Param */
    int index;                  /* index returned from Check_Param */

    status = LB_Check_Param(name, val, REFTREE_params, &result, &index);

    return(status);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Reftree_Init(LB *lb)

{
/*
 *  Function to initialize a refinement tree.  This creates the root and
 *  the first level of the tree, which corresponds to the initial coarse grid
 */
char *yo = "LB_Reftree_Init";
char msg[256];
struct LB_reftree_data_struct *reftree_data; /* data pointed to by lb */
LB_REFTREE *root;          /* Root of the refinement tree */
struct LB_reftree_hash_node **hashtab; /* hash table */
int nproc;                 /* number of processors */
ZOLTAN_ID_PTR local_gids;      /* coarse element Global IDs from user */
ZOLTAN_ID_PTR local_lids;      /* coarse element Local IDs from user */
ZOLTAN_ID_PTR lid, prev_lid;   /* temporary coarse element Local ID; used to pass
                              NULL to query functions when NUM_LID_ENTRIES=0 */
ZOLTAN_ID_PTR all_gids;        /* coarse element Global IDs from all procs */
int *assigned;             /* 1 if the element is assigned to this proc */
int *num_vert;             /* number of vertices for each coarse element */
int *reorder_nvert;        /* num_vert reordered by permutation "order" */
int root_vert[1];          /* fake number of vertices for the root */
int *vertices;             /* vertices for the coarse elements */
int *in_vertex;            /* "in" vertex for each coarse element */
int *out_vertex;           /* "out" vertex for each coarse element */
int in_order;              /* 1 if user is supplying order of the elements */
int num_obj;               /* number of coarse objects known to this proc */
int *num_obj_all;          /* num_obj from each processor */
int *displs;               /* running sum of num_obj_all */
int sum_num_obj;           /* full sum of num_obj_all */
int total_num_obj;         /* number of objects in the whole coarse grid */
int ierr;                  /* error flag from calls */
int final_ierr;            /* error flag returned by this routine */
int wdim;                  /* dimension of object weights */
int count;                 /* counter for number of objects */
int sum_vert;              /* summation of number of vertices of objects */
int *order;                /* permutation array for ordering coarse elements */
int found;                 /* flag for terminating first/next query loop */
int hashsize;              /* size of the hash table */
int i, j;                  /* loop counters */
int num_gid_entries = lb->Num_GID;  /* number of array entries in a global ID */
int num_lid_entries = lb->Num_LID;  /* number of array entries in a local ID */

  ZOLTAN_LB_TRACE_ENTER(lb, yo);

  ssize = 0;
  final_ierr = ZOLTAN_OK;

  if (lb->Obj_Weight_Dim == 0) {
    wdim = 1;
  } else {
    wdim = lb->Obj_Weight_Dim;
  }

  nproc = lb->Num_Proc;

  /*
   * Allocate the root of the refinement tree for this load balancing structure.
   * If a tree already exists, destroy it first.
   */

  if (lb->Data_Structure != NULL) LB_Reftree_Free_Structure(lb);

  root_vert[0] = 1;
  ierr = alloc_reftree_nodes(lb, &root, 1, root_vert);
  if (ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Error returned by alloc_reftree_nodes.");
    return(ierr);
  }

  /*
   * Initialize the root
   */

  root->children       = (LB_REFTREE *) NULL;
  root->num_child      = 0;
  root->num_vertex     = 0;
  root->in_vertex      = (int) NULL;
  root->out_vertex     = (int) NULL;
  root->assigned_to_me = 0;
  root->partition      = 0;

  for (i=0; i<wdim; i++) {
    root->weight[i] = 0.0;
    root->summed_weight[i] = 0.0;
    root->my_sum_weight[i] = 0.0;
  }

  /*
   * Allocate and initialize the hash table.
   */

  LB_Bind_Param(REFTREE_params, "REFTREE_HASH_SIZE", (void *) &hashsize);
  hashsize = DEFAULT_HASH_TABLE_SIZE;
  LB_Assign_Param_Vals(lb->Params, REFTREE_params, lb->Debug_Level, lb->Proc,
                       lb->Debug_Proc);

  hashtab = (struct LB_reftree_hash_node **)
            ZOLTAN_MALLOC(sizeof(struct LB_reftree_hash_node *)*hashsize);
  if (hashtab == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    LB_Reftree_Free_Structure(lb);
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return(ZOLTAN_MEMERR);
  }
  for (i=0; i<hashsize; i++)
    hashtab[i] = (struct LB_reftree_hash_node *)NULL;

  /*
   * set the lb pointer for later access to the refinement tree and hash table
   */

  reftree_data = (struct LB_reftree_data_struct *)
                 ZOLTAN_MALLOC(sizeof(struct LB_reftree_data_struct));
  reftree_data->reftree_root = root;
  reftree_data->hash_table = hashtab;
  reftree_data->hash_table_size = hashsize;
  lb->Data_Structure = (void *) reftree_data;

  /*
   * Get the list of initial elements known to this processor
   */


  /*
   * Get the number of objects
   */

  if (lb->Get_Num_Coarse_Obj == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Must register ZOLTAN_NUM_COARSE_OBJ_FN.");
    LB_Reftree_Free_Structure(lb);
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return(ZOLTAN_FATAL);
  }

  num_obj = lb->Get_Num_Coarse_Obj(lb->Get_Num_Coarse_Obj_Data, &ierr);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned from user function Get_Num_Coarse_Obj.");
    LB_Reftree_Free_Structure(lb);
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return(ierr);
  }

  /*
   * Get the objects, if the number is not 0
   */

  if (num_obj > 0) {

    num_obj += 1; /* allocate one extra spot for the last call to NEXT_OBJ */
    local_gids = ZOLTAN_ZOLTAN_MALLOC_GID_ARRAY(lb, num_obj);
    local_lids = ZOLTAN_ZOLTAN_MALLOC_LID_ARRAY(lb, num_obj);
    assigned   = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
    num_vert   = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
    vertices   = (int *) ZOLTAN_MALLOC(MAXVERT*num_obj*sizeof(int));
    in_vertex  = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
    out_vertex = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
    num_obj -= 1;

    if (local_gids == NULL || (num_lid_entries > 0 && local_lids == NULL) || 
        assigned   == NULL ||
        num_vert   == NULL || vertices   == NULL || in_vertex == NULL ||
        out_vertex == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&local_gids);
      ZOLTAN_FREE(&local_lids);
      ZOLTAN_FREE(&assigned);
      ZOLTAN_FREE(&num_vert);
      ZOLTAN_FREE(&vertices);
      ZOLTAN_FREE(&in_vertex);
      ZOLTAN_FREE(&out_vertex);
      LB_Reftree_Free_Structure(lb);
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      return(ZOLTAN_MEMERR);
    }

    if (lb->Get_Coarse_Obj_List != NULL) {

  /*
   * Get objects via list
   */

      lb->Get_Coarse_Obj_List(lb->Get_Coarse_Obj_List_Data, 
                              num_gid_entries, num_lid_entries,
                              local_gids, local_lids, 
                              assigned, num_vert, vertices,
                              &in_order, in_vertex, out_vertex, &ierr);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                      "Error returned from user function Get_Coarse_Obj_List.");
        ZOLTAN_FREE(&local_gids);
        ZOLTAN_FREE(&local_lids);
        ZOLTAN_FREE(&assigned);
        ZOLTAN_FREE(&num_vert);
        ZOLTAN_FREE(&vertices);
        ZOLTAN_FREE(&in_vertex);
        ZOLTAN_FREE(&out_vertex);
        LB_Reftree_Free_Structure(lb);
        ZOLTAN_LB_TRACE_EXIT(lb, yo);
        return(ierr);
      }

    }

    else if (lb->Get_First_Coarse_Obj != NULL &&
             lb->Get_Next_Coarse_Obj  != NULL) {

  /*
   * Get objects via first/next
   */

      sum_vert = 0;
      count = 0;
      lid = (num_lid_entries ? &(local_lids[count*num_lid_entries]) : NULL);
      found = lb->Get_First_Coarse_Obj(lb->Get_First_Coarse_Obj_Data,
                                       num_gid_entries, num_lid_entries,
                                       &(local_gids[count*num_gid_entries]), 
                                       lid,
                                       &assigned[count],
                                       &num_vert[count], &vertices[sum_vert],
                                       &in_order,
                                       &in_vertex[count], &out_vertex[count],
                                       &ierr);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                     "Error returned from user function Get_First_Coarse_Obj.");
        ZOLTAN_FREE(&local_gids);
        ZOLTAN_FREE(&local_lids);
        ZOLTAN_FREE(&assigned);
        ZOLTAN_FREE(&num_vert);
        ZOLTAN_FREE(&vertices);
        ZOLTAN_FREE(&in_vertex);
        ZOLTAN_FREE(&out_vertex);
        LB_Reftree_Free_Structure(lb);
        ZOLTAN_LB_TRACE_EXIT(lb, yo);
        return(ierr);
      }

      while (found && count <= num_obj) {
        sum_vert += num_vert[count];
        count += 1;
        prev_lid = (num_lid_entries ? &(local_lids[(count-1)*num_lid_entries]) 
                                    : NULL);
        lid = (num_lid_entries ? &(local_lids[count*num_lid_entries]) : NULL);
        found = lb->Get_Next_Coarse_Obj(lb->Get_Next_Coarse_Obj_Data,
                                      num_gid_entries, num_lid_entries,
                                      &(local_gids[(count-1)*num_gid_entries]), 
                                      prev_lid,
                                      &(local_gids[count*num_gid_entries]), 
                                      lid,
                                      &assigned[count],
                                      &num_vert[count], &vertices[sum_vert],
                                      &in_vertex[count], &out_vertex[count],
                                      &ierr);
        if (ierr) {
          ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                      "Error returned from user function Get_Next_Coarse_Obj.");
          ZOLTAN_FREE(&local_gids);
          ZOLTAN_FREE(&local_lids);
          ZOLTAN_FREE(&assigned);
          ZOLTAN_FREE(&num_vert);
          ZOLTAN_FREE(&vertices);
          ZOLTAN_FREE(&in_vertex);
          ZOLTAN_FREE(&out_vertex);
          LB_Reftree_Free_Structure(lb);
          ZOLTAN_LB_TRACE_EXIT(lb, yo);
          return(ierr);
        }
      }
      if (count != num_obj) {
        sprintf(msg, "Number of objects returned by "
                     "First/Next_Coarse_Obj = %d is not equal to the "
                     "number returned by Num_Coarse_Obj = %d\n",
                     count, num_obj);
        ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
        final_ierr = ZOLTAN_WARN;
      }
    }

    else {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Must define and register either "
        "ZOLTAN_COARSE_OBJ_LIST_FN or "
        "ZOLTAN_FIRST_COARSE_OBJ_FN/ZOLTAN_NEXT_COARSE_OBJ_FN pair.");
      ZOLTAN_FREE(&local_gids);
      ZOLTAN_FREE(&local_lids);
      ZOLTAN_FREE(&assigned);
      ZOLTAN_FREE(&num_vert);
      ZOLTAN_FREE(&vertices);
      ZOLTAN_FREE(&in_vertex);
      ZOLTAN_FREE(&out_vertex);
      LB_Reftree_Free_Structure(lb);
      ZOLTAN_LB_TRACE_EXIT(lb, yo);
      return(ZOLTAN_FATAL);
    }
  } /* endif (num_obj > 0) */

  /*
   * Communicate to get coarse grid objects unknown to this processor.
   */

  /*
   * First determine how many coarse objects are on each processor
   */

  num_obj_all = (int *)ZOLTAN_MALLOC(nproc*sizeof(int));
  displs = (int *)ZOLTAN_MALLOC(nproc*sizeof(int));
  if (num_obj_all == NULL || displs == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&local_gids);
    ZOLTAN_FREE(&local_lids);
    ZOLTAN_FREE(&assigned);
    ZOLTAN_FREE(&num_vert);
    ZOLTAN_FREE(&vertices);
    ZOLTAN_FREE(&in_vertex);
    ZOLTAN_FREE(&out_vertex);
    ZOLTAN_FREE(&num_obj_all);
    ZOLTAN_FREE(&displs);
    LB_Reftree_Free_Structure(lb);
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return(ZOLTAN_MEMERR);
  }

  MPI_Allgather((void *)&num_obj,1,MPI_INT,(void *)num_obj_all,1,MPI_INT,
                lb->Communicator);
  displs[0] = 0;
  for (i=1; i<nproc; i++) displs[i] = displs[i-1]+num_obj_all[i-1];
  sum_num_obj = displs[nproc-1] + num_obj_all[nproc-1];

  /*
   * Then get the coarse objects from all processors
   */

  all_gids = ZOLTAN_ZOLTAN_MALLOC_GID_ARRAY(lb, sum_num_obj);
  if (all_gids == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&local_gids);
    ZOLTAN_FREE(&local_lids);
    ZOLTAN_FREE(&assigned);
    ZOLTAN_FREE(&num_vert);
    ZOLTAN_FREE(&vertices);
    ZOLTAN_FREE(&in_vertex);
    ZOLTAN_FREE(&out_vertex);
    ZOLTAN_FREE(&num_obj_all);
    ZOLTAN_FREE(&displs);
    ZOLTAN_FREE(&all_gids);
    LB_Reftree_Free_Structure(lb);
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return(ZOLTAN_MEMERR);
  }

  /* KDDKDD Changed MPI_BYTE to ZOLTAN_ID_MPI_TYPE  */

  /* Account for number of array entries in an ID. */
  for (i=0; i<nproc; i++) {
    num_obj_all[i] = num_obj_all[i]*num_gid_entries;
    displs[i] = displs[i]*num_gid_entries;
  }

  MPI_Allgatherv((void *)local_gids,num_obj*num_gid_entries,ZOLTAN_ID_MPI_TYPE,
                 (void *)all_gids,num_obj_all,displs,ZOLTAN_ID_MPI_TYPE,
                 lb->Communicator);

  ZOLTAN_FREE(&displs);
  ZOLTAN_FREE(&num_obj_all);

  /*
   * Finally, build a list with each coarse grid element, beginning with
   * the ones this processor knows.  Also set the default order of the
   * elements as given by the user, with processor rank resolving duplicates
   */

  local_gids = ZOLTAN_LB_REALLOC_GID_ARRAY(lb, local_gids, sum_num_obj);
  order = (int *) ZOLTAN_MALLOC(sum_num_obj*sizeof(int));
  if (local_gids == NULL || order == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&local_gids);
    ZOLTAN_FREE(&local_lids);
    ZOLTAN_FREE(&assigned);
    ZOLTAN_FREE(&num_vert);
    ZOLTAN_FREE(&vertices);
    ZOLTAN_FREE(&in_vertex);
    ZOLTAN_FREE(&out_vertex);
    ZOLTAN_FREE(&all_gids);
    ZOLTAN_FREE(&order);
    LB_Reftree_Free_Structure(lb);
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return(ZOLTAN_MEMERR);
  }

/* TEMP this is terribly inefficient.  Probably better to sort all_gids to
        identify the duplicates.  Of course, it's not bad if the
        initial grid is really coarse */

  total_num_obj = num_obj;
  count = 0;
  for (i=0; i<sum_num_obj; i++) order[i] = -1;

  for (i=0; i<sum_num_obj; i++) {
    found = 0;
    for (j=0; j<total_num_obj && !found; j++) {
      if (ZOLTAN_LB_EQ_GID(lb, &(all_gids[i*num_gid_entries]),
                    &(local_gids[j*num_gid_entries]))) 
        found = 1;
    }
    if (found) {
      if (order[j-1] == -1) {
        order[j-1] = count;
        count += 1;
      }
    }
    else {
      ZOLTAN_LB_SET_GID(lb, &(local_gids[total_num_obj*num_gid_entries]), 
                     &(all_gids[i*num_gid_entries]));
      order[total_num_obj] = count;
      count += 1;
      total_num_obj += 1;
    }
  }

  if (count != total_num_obj) {
    sprintf(msg, "Number of objects counted while "
                 "setting default order = %d is not equal to the "
                 "number counted while getting objects from other procs "
                 "= %d.", count, total_num_obj);
    ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
    final_ierr = ZOLTAN_WARN;
  }

  ZOLTAN_FREE(&all_gids);

  num_vert = (int *) ZOLTAN_REALLOC(num_vert,total_num_obj*sizeof(int));
  if (num_vert == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&local_gids);
    ZOLTAN_FREE(&local_lids);
    ZOLTAN_FREE(&assigned);
    ZOLTAN_FREE(&num_vert);
    ZOLTAN_FREE(&vertices);
    ZOLTAN_FREE(&in_vertex);
    ZOLTAN_FREE(&out_vertex);
    ZOLTAN_FREE(&order);
    LB_Reftree_Free_Structure(lb);
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return(ZOLTAN_MEMERR);
  }

  for (i=num_obj; i<total_num_obj; i++) num_vert[i] = -1;

  /*
   * Determine the order of the coarse grid elements.
   * If the user supplies the order, it was set above.
   */

  if (!in_order) {

  /*
   * TEMP For now, require that the user provide the order.
   */

    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Currently not supporting automatic "
                    "determination of the order of the coarse grid objects.  "
                    "Using the order in which they were provided.");
    final_ierr = ZOLTAN_WARN;

  }

  /*
   * Copy the elements into the child list of the root
   */

  /*
   * Allocate the children of the root
   */

  reorder_nvert = (int *) ZOLTAN_MALLOC(total_num_obj*sizeof(int));
  if (reorder_nvert == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&local_gids);
    ZOLTAN_FREE(&local_lids);
    ZOLTAN_FREE(&assigned);
    ZOLTAN_FREE(&num_vert);
    ZOLTAN_FREE(&vertices);
    ZOLTAN_FREE(&in_vertex);
    ZOLTAN_FREE(&out_vertex);
    ZOLTAN_FREE(&order);
    LB_Reftree_Free_Structure(lb);
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return(ZOLTAN_MEMERR);
  }
/* use MAXVERT for coarse grid objects to avoid complicated reallocation
   during Reinit_Coarse */
  for (i=0; i<total_num_obj; i++) {
/*    reorder_nvert[order[i]] = num_vert[i]; */
    reorder_nvert[i] = MAXVERT;
  }

  ierr = alloc_reftree_nodes(lb, &(root->children), total_num_obj,
                             reorder_nvert);
  ZOLTAN_FREE(&reorder_nvert);
  if (ierr == ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Error returned by alloc_reftree_nodes.");
    ZOLTAN_FREE(&local_gids);
    ZOLTAN_FREE(&local_lids);
    ZOLTAN_FREE(&assigned);
    ZOLTAN_FREE(&num_vert);
    ZOLTAN_FREE(&vertices);
    ZOLTAN_FREE(&in_vertex);
    ZOLTAN_FREE(&out_vertex);
    ZOLTAN_FREE(&order);
    LB_Reftree_Free_Structure(lb);
    return(ierr);
  }

  root->num_child = total_num_obj;

  /*
   * Make sure the weights have been provided, if needed
   */

  if (lb->Obj_Weight_Dim != 0 && lb->Get_Child_Weight == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Must register ZOLTAN_CHILD_WEIGHT_FN.");
    ZOLTAN_FREE(&local_gids);
    ZOLTAN_FREE(&local_lids);
    ZOLTAN_FREE(&assigned);
    ZOLTAN_FREE(&num_vert);
    ZOLTAN_FREE(&vertices);
    ZOLTAN_FREE(&in_vertex);
    ZOLTAN_FREE(&out_vertex);
    ZOLTAN_FREE(&order);
    LB_Reftree_Free_Structure(lb);
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return(ZOLTAN_FATAL);
  }

  /*
   * For each coarse grid object ...
   */

  sum_vert = 0;
  for (i=0; i<total_num_obj; i++) {

  /*
   * Get the weight
   */

    if (lb->Obj_Weight_Dim == 0) {
  /* if an initial element is a leaf, the weight gets set to 1 later */
       *(root->children[order[i]].weight) = 0.0;
    }
    else if (num_vert[i] == -1) {
  /* if the element is not known to this processor, the weight is 0 */
       *(root->children[order[i]].weight) = 0.0;
    }
    else {
      lid = (num_lid_entries ? &(local_lids[i*num_lid_entries]) : NULL);
      lb->Get_Child_Weight(lb->Get_Child_Weight_Data,
                           num_gid_entries, num_lid_entries,
                           &(local_gids[i*num_gid_entries]),
                           lid, lb->Obj_Weight_Dim, 
                           root->children[order[i]].weight, &ierr);
    }
    for (j=0; j<wdim; j++) {
      root->children[order[i]].summed_weight[j] = 0.0;
      root->children[order[i]].my_sum_weight[j] = 0.0;
    }

  /*
   * Copy the vertices
   */

    for (j=0; j<num_vert[i]; j++) 
      root->children[order[i]].vertices[j] = vertices[sum_vert+j];
    if (num_vert[i] > 0) sum_vert += num_vert[i];

  /*
   * Copy from temporary arrays and set empty defaults
   */

    if (num_vert[i] == -1) {
  /* elements not known to this processor have more empty entries */
      ZOLTAN_LB_SET_GID(lb, root->children[order[i]].global_id,
                 &(local_gids[i*num_gid_entries]));
      ZOLTAN_LB_INIT_LID(lb, root->children[order[i]].local_id);
      root->children[order[i]].children       = (LB_REFTREE *) NULL;
      root->children[order[i]].num_child      = 0;
      root->children[order[i]].num_vertex     = num_vert[i];
      root->children[order[i]].in_vertex      = 0;
      root->children[order[i]].out_vertex     = 0;
      root->children[order[i]].assigned_to_me = 0;
      root->children[order[i]].partition      = 0;
    }
    else {
      ZOLTAN_LB_SET_GID(lb, root->children[order[i]].global_id,
                 &(local_gids[i*num_gid_entries]));
      ZOLTAN_LB_SET_LID(lb, root->children[order[i]].local_id,
                 &(local_lids[i*num_lid_entries]));
      root->children[order[i]].children       = (LB_REFTREE *) NULL;
      root->children[order[i]].num_child      = 0;
      root->children[order[i]].num_vertex     = num_vert[i];
      root->children[order[i]].in_vertex      = in_vertex[i];
      root->children[order[i]].out_vertex     = out_vertex[i];
      root->children[order[i]].assigned_to_me = assigned[i];
      root->children[order[i]].partition      = 0;
    }

  /*
   * Add it to the hash table
   */

    LB_Reftree_Hash_Insert(lb, &(root->children[order[i]]),hashtab,hashsize);

  }

  /*
   * clean up and return error code
   */

  ZOLTAN_FREE(&local_gids);
  ZOLTAN_FREE(&local_lids);
  ZOLTAN_FREE(&assigned);
  ZOLTAN_FREE(&num_vert);
  ZOLTAN_FREE(&vertices);
  ZOLTAN_FREE(&in_vertex);
  ZOLTAN_FREE(&out_vertex);
  ZOLTAN_FREE(&order);
  ZOLTAN_LB_TRACE_EXIT(lb, yo);
  return(final_ierr);
}


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Reftree_Build(LB *lb)

{
/*
 * Function to build a refinement tree
 */
char *yo = "LB_Reftree_Build";
LB_REFTREE *root;          /* Root of the refinement tree */
int ierr;                  /* Error code returned by called functions */
int i;                     /* loop counter */

  /*
   * Initialize the tree, if not already there, and set the root.  If already
   * there, reinitialize coarse grid.
   */

  if (lb->Data_Structure == NULL) {
    ierr = LB_Reftree_Init(lb);
    if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Error returned from LB_Reftree_Init.");
      return(ierr);
    }
  }
  else {
    LB_Reftree_Reinit_Coarse(lb);
  }
  root = ((struct LB_reftree_data_struct *)lb->Data_Structure)->reftree_root;

  /*
   * Verify the required child query functions are registered
   */

  if (lb->Get_Num_Child == NULL || lb->Get_Child_List == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Must register ZOLTAN_NUM_CHILD_FN"
            " and ZOLTAN_CHILD_LIST_FN.");
    LB_Reftree_Free_Structure(lb);
    return(ZOLTAN_FATAL);
  }

  /*
   * For each child of the root, build the subtree rooted at that child
   * and, if the weights are not provided, set its weight if it is a leaf.
   * Skip elements not known to this processor.
   */

  for (i=0; i<root->num_child; i++) {
    if ( (root->children[i]).num_vertex != -1 ) {
      ierr = LB_Reftree_Build_Recursive(lb,&(root->children[i]));
      if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
        ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                       "Error returned from LB_Reftree_Build_Recursive.");
        return(ierr);
      }
    }
  }

  return(ZOLTAN_OK);
}

static int LB_Reftree_Build_Recursive(LB *lb,LB_REFTREE *subroot)

{
/*
 * Recursive function to traverse a tree while building it
 */
static int TEMP_first_warning = 1; /* TEMP until ref_type is fully supported */
char *yo = "LB_Reftree_Build_Recursive";
char msg[256];
int ierr;                  /* error code called routines */
int final_ierr;            /* error code returned by this routine */
int num_obj;               /* number of children returned by user */
ZOLTAN_ID_PTR lid;             /* temporary coarse element Local ID; used to pass
                              NULL to query functions when NUM_LID_ENTRIES=0 */
int *reorder_nvert;        /* num_vert reordered by permutation "order" */
ZOLTAN_REF_TYPE ref_type;  /* type of refinement that creates children */
int wdim;                  /* dimension for weights */
int i, j;                  /* loop counters */
int sum_vert;              /* running sum of the number of vertices */
struct LB_reftree_hash_node **hashtab; /* hash tree */
int hashsize;              /* size of the hash table */
int num_gid_entries = lb->Num_GID;  /* number of array entries in a global ID */
int num_lid_entries = lb->Num_LID;  /* number of array entries in a local ID */
int children_agree;        /* flag, true if all children of a node in the
                              refinement tree agree with data from GET_CHILD */
int existing;              /* existing child that agrees with GET_CHILD data */

  final_ierr = ZOLTAN_OK;
  if (lb->Obj_Weight_Dim == 0) {
    wdim = 1;
  } else {
    wdim = lb->Obj_Weight_Dim;
  }

  /*
   * Print a warning if a nonexistent subroot is passed in
   */

  if (subroot == NULL) {
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Called with nonexistent subroot.");
    return(ZOLTAN_WARN);
  }

  /*
   * Get the number of children of this node
   */

  num_obj = lb->Get_Num_Child(lb->Get_Num_Child_Data, 
                              num_gid_entries, num_lid_entries,
                              subroot->global_id, subroot->local_id, &ierr);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned from user function Get_Num_Child.");
    LB_Reftree_Free_Structure(lb);
    return(ierr);
  }

  /*
   * If there are no children, make sure the tree has no children, and
   * set the weight if it is not user provided,
   * and return.  The default is to use 1.0 for leaves and 0.0 for others.
   */

  if (num_obj == 0) {
    if (subroot->num_child != 0) {
      LB_Reftree_Free_Subtree(lb, subroot);
    }
    if (lb->Obj_Weight_Dim == 0) *(subroot->weight) = 1.0;
    return(ZOLTAN_OK);
  }

  /*
   * Get the children
   */

  if (num_obj > ssize) {
    if (ssize > 0) {
      ZOLTAN_FREE(&slocal_gids);
      ZOLTAN_FREE(&slocal_lids);
      ZOLTAN_FREE(&sassigned);
      ZOLTAN_FREE(&snum_vert);
      ZOLTAN_FREE(&svertices);
      ZOLTAN_FREE(&sin_vertex);
      ZOLTAN_FREE(&sout_vertex);
      ZOLTAN_FREE(&svert1);
    }
    slocal_gids = ZOLTAN_ZOLTAN_MALLOC_GID_ARRAY(lb, num_obj);
    slocal_lids = ZOLTAN_ZOLTAN_MALLOC_LID_ARRAY(lb, num_obj);
    sassigned   = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
    snum_vert   = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
    svertices   = (int *) ZOLTAN_MALLOC(MAXVERT*num_obj*sizeof(int));
    sin_vertex  = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
    sout_vertex = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
    svert1      = (int *) ZOLTAN_MALLOC((num_obj+1)*sizeof(int));
    sorder      = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
    ssize = num_obj;
    if (slocal_gids == NULL || (num_lid_entries > 0 && slocal_lids == NULL) || 
        sassigned   == NULL ||
        snum_vert   == NULL || svertices   == NULL || sin_vertex == NULL ||
        sout_vertex == NULL || svert1      == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&slocal_gids);
      ZOLTAN_FREE(&slocal_lids);
      ZOLTAN_FREE(&sassigned);
      ZOLTAN_FREE(&snum_vert);
      ZOLTAN_FREE(&svertices);
      ZOLTAN_FREE(&sin_vertex);
      ZOLTAN_FREE(&sout_vertex);
      ZOLTAN_FREE(&svert1);
      ZOLTAN_FREE(&sorder);
      ssize = 0;
      LB_Reftree_Free_Structure(lb);
      return(ZOLTAN_MEMERR);
    }
  }
  lb->Get_Child_List(lb->Get_Child_List_Data, 
                     num_gid_entries, num_lid_entries,
                     subroot->global_id, subroot->local_id, 
                     slocal_gids, slocal_lids, sassigned,
                     snum_vert, svertices, &ref_type, sin_vertex, sout_vertex,
                     &ierr);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned from user function Get_Child_List.");
    LB_Reftree_Free_Structure(lb);
    return(ierr);
  }

  /*
   * Set the start of the list of vertices for each child
   */

  svert1[0] = 0;
  for (i=0; i<num_obj; i++) svert1[i+1] = svert1[i] + snum_vert[i];

  /*
   * If there already exist children, then make sure they are exactly the
   * same.  Otherwise, delete them and rebuild the tree from here.
   */

  /*
   * check that the number of children agree and that each child agrees
   * with an existing child in GID, LID and vertices
   */

  children_agree = 1;
  if (subroot->num_child == 0) {
    children_agree = 0;
  } else {
    if (subroot->num_child != num_obj) {
      children_agree = 0;
    } else {
      for (i=0; i<num_obj && children_agree; i++) {
        existing = -1;
        for (j=0; j<subroot->num_child && existing==-1; j++) {
          if (ZOLTAN_LB_EQ_GID(lb, subroot->children[j].global_id,
                        &(slocal_gids[i*num_gid_entries]))) {
            existing = j;
          }
        }
        if (existing == -1) {
          children_agree = 0;
        } else {
          for (j=0; j<num_lid_entries; j++) {
            if (subroot->children[existing].local_id[j] !=
                slocal_lids[i*num_lid_entries+j]) {
              children_agree = 0;
            }
          }
          if (subroot->children[existing].num_vertex != snum_vert[i]) {
            children_agree = 0;
          } else {
            if (snum_vert[i] != 0) {
              for (j=0; j<snum_vert[i] && children_agree; j++) {
                if (subroot->children[existing].vertices[j] != svertices[svert1[i]+j]) {
                  children_agree = 0;
                }
              }
            }
          }
/* set new value of assigned just in case we keep the children */
          subroot->children[existing].assigned_to_me = sassigned[i];
        }
      }
    }
  }

  /*
   * If the children do not agree, then get rid of them and rebuild
   */

  if (!children_agree) LB_Reftree_Free_Subtree(lb, subroot);

  if (subroot->num_child != 0) {

  /*
   * If the children were kept, set the current weights
   */

    for (i=0; i<subroot->num_child; i++) {
      if (lb->Obj_Weight_Dim == 0) {
         *(subroot->children[i].weight) = 0.0;
      }
      else {
        lid = (num_lid_entries ? subroot->children[i].local_id : NULL);
        lb->Get_Child_Weight(lb->Get_Child_Weight_Data,
                             num_gid_entries, num_lid_entries,
                             subroot->children[i].global_id,
                             lid, 
                             lb->Obj_Weight_Dim,
                             subroot->children[i].weight, &ierr);
      }
      for (j=0; j<wdim; j++) {
        subroot->children[i].summed_weight[j] = 0.0;
        subroot->children[i].my_sum_weight[j] = 0.0;
      }
    }

  } else {

  /*
   * If the children did not already exist or were removed, add them
   */

  /*
   * Determine the order of the children
   */

  /*
   * TEMP until code is supplied to support these refinement types
   */

    switch (ref_type) {
    case ZOLTAN_HEX3D_OCT:
      if (TEMP_first_warning) {
        ZOLTAN_PRINT_WARN(lb->Proc, yo, "Currently not supporting "
                        "automatic ordering of elements for refinement type "
                        "ZOLTAN_HEX3D_OCT.  Using ZOLTAN_OTHER_REF.");
        TEMP_first_warning = 0;
        final_ierr = ZOLTAN_WARN;
      }
      ref_type = ZOLTAN_OTHER_REF;
      break;
    default:
      break;
    }

  /* end TEMP */

  /*
   *  Refinement type dependent determination of the order of the children
   *  and the in/out vertices
   */

    switch (ref_type) {

    case ZOLTAN_IN_ORDER:
      for (i=0; i<num_obj; i++) sorder[i] = i;
      break;
    case ZOLTAN_TRI_BISECT:
      ierr = order_tri_bisect(lb,svert1,sorder,svertices,sin_vertex,sout_vertex,
                              subroot);
      break;
    case ZOLTAN_QUAD_QUAD:
      ierr = order_quad_quad(lb,svert1,sorder,svertices,sin_vertex,sout_vertex,
                             subroot);
      break;
    case ZOLTAN_HEX3D_OCT:
    /* TEMP */
      ZOLTAN_PRINT_WARN(lb->Proc, yo, "Oops, still got into case for HEX3D_OCT.");
      for (i=0; i<num_obj; i++) sorder[i] = i;
      break;
    case ZOLTAN_OTHER_REF:
      ierr = order_other_ref(lb, subroot, num_obj, snum_vert, svert1, svertices,
                             sorder, sin_vertex, sout_vertex);
      break;

  /*
   * Default case if a bad value gets returned; use them in order.
   */
    default:
      sprintf(msg, "Unknown value returned for ref_type"
              " = %d.  Using children in order provided.",ref_type);
      ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
      for (i=0; i<num_obj; i++) sorder[i] = i;
      final_ierr = ZOLTAN_WARN;
    }

  /*
   * Copy the children into the child list of the subroot
   */

  /*
   * Allocate the children
   */

    if (subroot->children != NULL) {
      ZOLTAN_PRINT_WARN(lb->Proc, yo, "children already existed; memory"
                      " leak potential.");
      final_ierr = ZOLTAN_WARN;
    }

    reorder_nvert = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
    if (reorder_nvert == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_Reftree_Free_Structure(lb);
      return(ZOLTAN_MEMERR);
    }
    for (i=0; i<num_obj; i++) {
      reorder_nvert[sorder[i]] = snum_vert[i];
    }

    ierr = alloc_reftree_nodes(lb, &(subroot->children), num_obj, reorder_nvert);

    ZOLTAN_FREE(&reorder_nvert);

    subroot->num_child = num_obj;

    hashtab  = ((struct LB_reftree_data_struct *)lb->Data_Structure)->hash_table;
    hashsize = ((struct LB_reftree_data_struct *)lb->Data_Structure)->hash_table_size;

  /*
   * For each child ...
   */

    sum_vert = 0;
    for (i=0; i<num_obj; i++) {

  /*
   * Get the weight
   */

      if (lb->Obj_Weight_Dim == 0) {
         *(subroot->children[sorder[i]].weight) = 0.0;
      }
      else {
        lid = (num_lid_entries ? &(slocal_lids[i*num_lid_entries]) : NULL);
        lb->Get_Child_Weight(lb->Get_Child_Weight_Data,
                             num_gid_entries, num_lid_entries,
                             &(slocal_gids[i*num_gid_entries]),
                             lid, 
                             lb->Obj_Weight_Dim,
                             subroot->children[sorder[i]].weight, &ierr);
      }
      for (j=0; j<wdim; j++) {
        subroot->children[sorder[i]].summed_weight[j] = 0.0;
        subroot->children[sorder[i]].my_sum_weight[j] = 0.0;
      }

  /*
   * Copy the vertices
   */

      for (j=0; j<snum_vert[i]; j++)
        subroot->children[sorder[i]].vertices[j] = svertices[sum_vert+j];
      if (snum_vert[i] > 0) sum_vert += snum_vert[i];

  /*
   * Copy from temporary arrays and set empty defaults
   */

      ZOLTAN_LB_SET_GID(lb, subroot->children[sorder[i]].global_id,
                 &(slocal_gids[i*num_gid_entries]));
      ZOLTAN_LB_SET_LID(lb, subroot->children[sorder[i]].local_id,
                 &(slocal_lids[i*num_lid_entries]));
      subroot->children[sorder[i]].children       = (LB_REFTREE *) NULL;
      subroot->children[sorder[i]].num_child      = 0;
      subroot->children[sorder[i]].num_vertex     = snum_vert[i];
      subroot->children[sorder[i]].in_vertex      = sin_vertex[i];
      subroot->children[sorder[i]].out_vertex     = sout_vertex[i];
      subroot->children[sorder[i]].assigned_to_me = sassigned[i];
      subroot->children[sorder[i]].partition      = 0;

  /*
   * Add it to the hash table
   */

      LB_Reftree_Hash_Insert(lb, &(subroot->children[sorder[i]]),hashtab,hashsize);

    }
  }

  /*
   * recursively do the children
   */

  for (i=0; i<subroot->num_child; i++) {
    ierr = LB_Reftree_Build_Recursive(lb,&(subroot->children[i]));
    if (ierr) final_ierr = ierr;
  }

  return(final_ierr);

}

/*****************************************************************************/

static int order_tri_bisect(LB *lb, int *vert1, int *order, int *vertices,
                     int *in_vertex, int *out_vertex, LB_REFTREE *subroot)
{
/*
 * Function to determine the order of the children and in/out vertices
 * when refinement is done by bisecting triangles.  Determine which of
 * the first and second child has the in_vertex and out_vertex, and find the
 * common vertex to go between them.
 */
char *yo = "order_tri_bisect";
int i, j;                  /* loop indices */
int parents_vert[6];       /* cross index between children and parent */
int parent_in;             /* index of the parent in vertex */
int parent_out;            /* index of the parent out vertex */
int parent_third;          /* index of the parent triangle non in/out vert */
int not_parent[2];         /* index of a vertex the parent does not have */
int has_in[2];             /* flag for an element having the parent in vert */
int has_out[2];            /* flag for an element having the parent out vert */
int has_third[2];          /* flag for a triangle having the parent non in/out*/
int bad_case;              /* flag for failing to identify order */

  /* verify that 3 vertices were given for each triangle; if not, punt */
  if (vert1[1] != 3 || vert1[2] != 6) {
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Incorrect number of vertices "
                                "given for bisected triangles.");
    order[0] = 0;
    order[1] = 1;
    in_vertex[0] = vertices[vert1[0]];
    out_vertex[0] = vertices[vert1[0]+1];
    in_vertex[1] = vertices[vert1[1]];
    out_vertex[1] = vertices[vert1[1]+1];
    return(ZOLTAN_WARN);
  }

  /* determine the relationship between the parent's vertices and
     the children's vertices */
  for (i=0; i<6; i++) {
    parents_vert[i] = -1;
    for (j=0; j<3; j++) {
      if (vertices[i] == subroot->vertices[j]) parents_vert[i] = j;
    }
  }

  /* determine the location of the parents in and out vertices */
  parent_in = -1; parent_out = -1; parent_third = -1;
  for (i=0; i<3; i++) {
    if (subroot->vertices[i] == subroot->in_vertex) {
      parent_in = i;
    }
    else if (subroot->vertices[i] == subroot->out_vertex) {
      parent_out = i;
    }
    else {
    parent_third = i;
    }
  }
  if (parent_in == -1 || parent_out == -1 || parent_third == -1) {
    /* failed to locate one of them */
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Could not locate in and out "
                                "vertices in the parent.");
    order[0] = 0;
    order[1] = 1;
    in_vertex[0] = vertices[vert1[0]];
    out_vertex[0] = vertices[vert1[0]+1];
    in_vertex[1] = vertices[vert1[1]];
    out_vertex[1] = vertices[vert1[1]+1];
    return(ZOLTAN_WARN);
  }

  /* find the vertex that the parent doesn't have */
  if (parents_vert[0] == -1) not_parent[0] = 0;
  if (parents_vert[1] == -1) not_parent[0] = 1;
  if (parents_vert[2] == -1) not_parent[0] = 2;
  if (parents_vert[3] == -1) not_parent[1] = 3;
  if (parents_vert[4] == -1) not_parent[1] = 4;
  if (parents_vert[5] == -1) not_parent[1] = 5;

  /* see which children have which special vertices */
  if (parents_vert[0] == parent_in || parents_vert[1] == parent_in ||
      parents_vert[2] == parent_in) has_in[0] = 1;
  else has_in[0] = 0;
  if (parents_vert[0] == parent_out || parents_vert[1] == parent_out ||
      parents_vert[2] == parent_out) has_out[0] = 1;
  else has_out[0] = 0;
  if (parents_vert[0] == parent_third || parents_vert[1] == parent_third ||
      parents_vert[2] == parent_third) has_third[0] = 1;
  else has_third[0] = 0;
  if (parents_vert[3] == parent_in || parents_vert[4] == parent_in ||
      parents_vert[5] == parent_in) has_in[1] = 1;
  else has_in[1] = 0;
  if (parents_vert[3] == parent_out || parents_vert[4] == parent_out ||
      parents_vert[5] == parent_out) has_out[1] = 1;
  else has_out[1] = 0;
  if (parents_vert[3] == parent_third || parents_vert[4] == parent_third ||
      parents_vert[5] == parent_third) has_third[1] = 1;
  else has_third[1] = 0;

  /* look for the case for this refinement */
  bad_case = 0;
  if (has_in[0]) {
    if (has_out[1]) {
      order[0] = 0; order[1] = 1;
      in_vertex[0] = subroot->vertices[parent_in];
      out_vertex[1] = subroot->vertices[parent_out];
      if (has_third[0] && has_third[1]) {
        out_vertex[0] = subroot->vertices[parent_third];
        in_vertex[1] = subroot->vertices[parent_third];
      }else{
        out_vertex[0] = vertices[not_parent[0]];
        in_vertex[1] = vertices[not_parent[1]];
      }
    }
    else if (has_in[1]) {
      if (has_out[0]) {
        order[0] = 1; order[1] = 0;
        in_vertex[1] = subroot->vertices[parent_in];
        out_vertex[0] = subroot->vertices[parent_out];
        if (has_third[0] && has_third[1]) {
          out_vertex[1] = subroot->vertices[parent_third];
          in_vertex[0] = subroot->vertices[parent_third];
        }else{
          out_vertex[1] = vertices[not_parent[1]];
          in_vertex[0] = vertices[not_parent[0]];
        }
      }else{ /* impossible case, no one has the out vertex */
        bad_case = 1;
        order[0] = 0; order[1] = 1;
        in_vertex[0] = subroot->vertices[parent_in];
        out_vertex[0] = subroot->vertices[parent_third];
        in_vertex[1] = subroot->vertices[parent_third];
        out_vertex[1] = subroot->vertices[parent_in];
      }
    }else{ /* impossible case, second child has neither in nor out */
      bad_case = 1;
      order[0] = 0; order[1] = 1;
      in_vertex[0] = subroot->vertices[parent_in];
      out_vertex[0] = subroot->vertices[parent_third];
      in_vertex[1] = vertices[3];
      out_vertex[1] = vertices[4];
    }
  }
  else if (has_out[0]) {
    if (has_in[1]) {
      order[0] = 1; order[1] = 0;
      in_vertex[1] = subroot->vertices[parent_in];
      out_vertex[0] = subroot->vertices[parent_out];
      if (has_third[0] && has_third[1]) {
        out_vertex[1] = subroot->vertices[parent_third];
        in_vertex[0] = subroot->vertices[parent_third];
      }else{
        out_vertex[1] = vertices[not_parent[1]];
        in_vertex[0] = vertices[not_parent[0]];
      }
    }else{ /* impossible case, no one has the in vertex */
      bad_case = 1;
      order[0] = 0; order[1] = 1;
      in_vertex[0] = subroot->vertices[parent_out];
      out_vertex[0] = subroot->vertices[parent_third];
      in_vertex[1] = subroot->vertices[parent_third];
      out_vertex[1] = subroot->vertices[parent_out];
    }
  }else{ /* impossible case, first child has neither in nor out */
    bad_case = 1;
    order[0] = 0; order[1] = 1;
    in_vertex[0] = vertices[0];
    out_vertex[0] = vertices[1];
    in_vertex[1] = vertices[3];
    out_vertex[1] = vertices[4];
  }
  if (bad_case) {
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Vertices of children did not "
                    "match the in and out vertices of parent.");
    return(ZOLTAN_WARN);
  }
  else {
    return(ZOLTAN_OK);
  }
}

/*****************************************************************************/

static int order_quad_quad(LB *lb, int *vert1, int *order, int *vertices,
                     int *in_vertex, int *out_vertex, LB_REFTREE *subroot)
{
/*
 * Function to determine the order of the children and in/out vertices
 * when refinement is done by quadrasecting quadrilaterals.
 */

int i,j,k,found,shared[3],ord[4];
char *yo = "order_quad_quad";

  /* verify that 4 vertices were given for each quadrilateral; if not, punt */
  if (vert1[1] != 4 || vert1[2] != 8 || vert1[3] != 12) {
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Incorrect number of vertices "
                                "given for quadrasected quadrilaterals.");
    for (i=0; i<4; i++) {
      order[i] = i;
      in_vertex[i] = vertices[vert1[i]];
      out_vertex[i] = vertices[vert1[i]+3];
    }
    return(ZOLTAN_WARN);
  }

  /* find the child that contains the in_vertex and make it first */

  found = 0;
  for (i=0; i<4 && !found; i++) {
    for (j=0; j<4 && !found; j++) {
      if (vertices[4*i+j] == subroot->in_vertex) {
         ord[0] = i;
         found = 1;
      }
    }
  }
  if (!found) {
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Couldn't find in_vertex in children");
    for (i=0; i<4; i++) {
      order[i] = i;
      in_vertex[i] = vertices[vert1[i]];
      out_vertex[i] = vertices[vert1[i]+3];
    }
    return(ZOLTAN_WARN);
  }

  /* find the child that contains the out_vertex and make it last */

  found = 0;
  for (i=0; i<4 && !found; i++) {
    for (j=0; j<4 && !found; j++) {
      if (vertices[4*i+j] == subroot->out_vertex) {
         ord[3] = i;
         found = 1;
      }
    }
  }
  if (!found) {
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Couldn't find out_vertex in children");
    for (i=0; i<4; i++) {
      order[i] = i;
      in_vertex[i] = vertices[vert1[i]];
      out_vertex[i] = vertices[vert1[i]+3];
    }
    return(ZOLTAN_WARN);
  }

  /* find a child that shares two vertices with the first child and is
     not the last child, and make it second */

  found = 0;
  for (k=0; k<4 && found!=2; k++) {
    if (k != ord[0] && k != ord[3]) {
      found = 0;
      for (j=0; j<4 && found!=2; j++) {
        for (i=0; i<4 && found!=2; i++) {
          if (vertices[4*k+j] == vertices[4*ord[0]+i]) {
            shared[found] = vertices[4*k+j];
            found = found + 1;
          }
        }
      }
    }
    if (found == 2) {
      ord[1] = k;
    }
  }
  if (found != 2) {
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Couldn't find second child of quadrasection");
    for (i=0; i<4; i++) {
      order[i] = i;
      in_vertex[i] = vertices[vert1[i]];
      out_vertex[i] = vertices[vert1[i]+3];
    }
    return(ZOLTAN_WARN);
  }

  /* the remaining child is given by six minus the sum of the others */

  ord[2] = 6 - ord[0] - ord[1] - ord[3];

  /* determine which vertex shared by the first and second children is also
     shared by the last child, and hence is the middle of the parent, and
     place that one in shared[0] */

  found = 0;
  for (j=0; j<4 && !found; j++) {
    if (shared[0] == vertices[4*ord[3]+j]) {
      found = 1;
    }
    if (shared[1] == vertices[4*ord[3]+j]) {
      shared[1] = shared[0];
      shared[0] = vertices[4*ord[3]+j];
      found = 1;
    }
  }
  if (!found) {
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Couldn't find central node of quadrasection");
    for (i=0; i<4; i++) {
      order[i] = i;
      in_vertex[i] = vertices[vert1[i]];
      out_vertex[i] = vertices[vert1[i]+3];
    }
    return(ZOLTAN_WARN);
  }

  /* find the other vertex shared by the third and fourth children */

  found = 0;
  for (j=0; j<4 && !found; j++) {
    if (vertices[4*ord[2]+j] != shared[0]) {
      for (i=0; i<4 && !found; i++) {
        if (vertices[4*ord[2]+j] == vertices[4*ord[3]+i]) {
          shared[2] = vertices[4*ord[2]+j];
          found = 1;
        }
      }
    }
  }
  if (!found) {
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Couldn't find shared vertex of 3rd and 4th child");
    for (i=0; i<4; i++) {
      order[i] = i;
      in_vertex[i] = vertices[vert1[i]];
      out_vertex[i] = vertices[vert1[i]+3];
    }
    return(ZOLTAN_WARN);
  }

  /* invert the permutation matrix */

  for (i=0; i<4; i++) {
    order[ord[i]] = i;
  }

  /* set the in/out vertices */

   in_vertex[ord[0]] = subroot->in_vertex;
  out_vertex[ord[0]] = shared[1];
   in_vertex[ord[1]] = shared[1];
  out_vertex[ord[1]] = shared[0];
   in_vertex[ord[2]] = shared[0];
  out_vertex[ord[2]] = shared[2];
   in_vertex[ord[3]] = shared[2];
  out_vertex[ord[3]] = subroot->out_vertex;

  return(ZOLTAN_OK);
}

/*****************************************************************************/

static int order_other_ref(LB *lb, LB_REFTREE *parent, int num_child, 
                    int *num_vert,
                    int *vert1, int *vertices, int *order, int *in_vertex,
                    int *out_vertex)
{
/*
 * Function to determine the order of the children for an undetermined
 * refinement scheme.  This is expensive as it performs a tree search
 * to solve this NP hard problem, but it should work for any refinement.
 */

char *yo = "order_other_ref";
int i, j, vi, vj;   /* loop counters */
int *has_in;        /* flag for children having in vertex */
int *has_out;       /* flag for children having out vertex */
int **share_vert;   /* number of vertices shared by two elements */
int max_share;      /* maximum number of vertices shared by two elements */
int solved;         /* flag for having found the solution */
int final_ierr;     /* error code returned */
int *on_path;       /* flag for already placed element on path */

  final_ierr = ZOLTAN_OK;

  /*
   * Determine which elements contain the in and out vertices of the parent
   */

  has_in = (int *) ZOLTAN_MALLOC(num_child*sizeof(int));
  has_out = (int *) ZOLTAN_MALLOC(num_child*sizeof(int));
  if (has_in == NULL || has_out == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&has_in);
    ZOLTAN_FREE(&has_out);
    return(ZOLTAN_MEMERR);
  }

  for (i=0; i<num_child; i++) {
    has_in[i] = 0;
    has_out[i] = 0;
    for (j=0; j<num_vert[i] && !has_in[i]; j++)
      if (vertices[vert1[i]+j] == parent->in_vertex) has_in[i] = 1;
    for (j=0; j<num_vert[i] && !has_out[i]; j++)
      if (vertices[vert1[i]+j] == parent->out_vertex) has_out[i] = 1;
  }

  /*
   * Determine which elements share vertices other than the in/out vertices
   */

  share_vert = (int **) ZOLTAN_MALLOC(num_child*sizeof(int *));
  if (share_vert == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&share_vert);
    ZOLTAN_FREE(&has_in);
    ZOLTAN_FREE(&has_out);
    return(ZOLTAN_MEMERR);
  }
  for (i=0; i<num_child; i++) {
    share_vert[i] = (int *) ZOLTAN_MALLOC(num_child*sizeof(int));
    if (share_vert[i] == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      for (j=0; j<=i; j++) ZOLTAN_FREE(&(share_vert[j]));
      ZOLTAN_FREE(&share_vert);
      ZOLTAN_FREE(&has_in);
      ZOLTAN_FREE(&has_out);
      return(ZOLTAN_MEMERR);
    }
  }

  max_share = 0;
  for (i=0; i<num_child; i++) {
    share_vert[i][i] = 1;
    for (j=i+1; j<num_child; j++) {
      share_vert[i][j] = 0;
      share_vert[j][i] = 0;
      for (vi=0; vi<num_vert[i]; vi++) {
        for (vj=0; vj<num_vert[j]; vj++) {
          if (vertices[vert1[i]+vi] == vertices[vert1[j]+vj]) {
            if (vertices[vert1[i]+vi] != parent->in_vertex &&
                vertices[vert1[i]+vi] != parent->out_vertex) {
              share_vert[i][j] = share_vert[i][j] + 1;
              share_vert[j][i] = share_vert[i][j];
            }
          }
        }
      }
      if (share_vert[i][j] > max_share) {
        max_share = share_vert[i][j];
      }
    }
  }

  /*
   * Perform tree search to find solution
   */

  solved = 0;
  on_path = (int *) ZOLTAN_MALLOC(num_child*sizeof(int));
  if (on_path == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    for (j=0; j<=i; j++) ZOLTAN_FREE(&(share_vert[j]));
    ZOLTAN_FREE(&on_path);
    ZOLTAN_FREE(&share_vert);
    ZOLTAN_FREE(&has_in);
    ZOLTAN_FREE(&has_out);
    return(ZOLTAN_MEMERR);
  }
  for (i=0; i<num_child; i++) on_path[i]=0;

  /*
   * Try each element with the in vertex to start the path
   */

  for (i=0; i<num_child && !solved; i++) {
    if (has_in[i]) {
      order_other_ref_recur(i,0,order,on_path,num_child,
                            has_out,share_vert,max_share,&solved);
    }
  }

  /*
   * This should have found a solution, but if not then use given order
   */

  if (!solved) {
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Couldn't find path through children."
                                "  Using given order.");
    for (i=0; i<num_child; i++) order[i] = i;
    final_ierr = ZOLTAN_WARN;
  }

  ZOLTAN_FREE(&on_path);
  ZOLTAN_FREE(&share_vert);
  ZOLTAN_FREE(&has_out);

  /*
   * Finally, determine the in and out vertices of each child
   */

  in_vertex[order[0]] = parent->in_vertex;
  out_vertex[order[num_child-1]] = parent->out_vertex;
  solved = find_inout(0, num_child, num_vert, vert1, vertices, in_vertex,
                      out_vertex, order);
  if (!solved) {
    ZOLTAN_PRINT_WARN(lb->Proc, yo, "Couldn't find good set of in/out"
                    " vertices.  Using first and second.\n");
    for (i=0; i<num_child; i++) {
      in_vertex[i]  = vertices[vert1[i]];
      out_vertex[i] = vertices[vert1[i]+1];
    }
    final_ierr = ZOLTAN_WARN;
  }

  /*
   * Invert the permutation matrix (order) to agree with it's usage in
   * LB_Reftree_Build_Recursive, using has_in as workspace
   */

  for (i=0; i<num_child; i++) {
    has_in[order[i]] = i;
  }
  for (i=0; i<num_child; i++) {
    order[i] = has_in[i];
  }
  ZOLTAN_FREE(&has_in);

  return(final_ierr);
}

/*****************************************************************************/

static void order_other_ref_recur(int new_entry, int level, int *order, 
                          int *on_path,
                          int num_child, int *has_out, int **share_vert,
                          int max_share, int *solved)
{
/*
 * Recursive routine to search the solution space tree
 */
int i, nshare;

  if (level == num_child-1) {

  /*
   * End of a path, success if this element has the out vertex
   */
    if (has_out[new_entry]) {
      order[level] = new_entry;
      *solved = 1;
    }
    else {
      *solved = 0;
    }
  }

  else {

  /*
   * Add this element to the proposed path
   */

    order[level] = new_entry;
    on_path[new_entry] = 1;

  /*
   * Try each element that is not already on the path and shares a vertex
   * with the current new entry, starting first with those that share the
   * most vertices to give a preference to going through faces
   */

    for (nshare = max_share; nshare>0 && !(*solved); nshare--) {
      for (i=0; i<num_child && !(*solved); i++) {
        if (!on_path[i] && share_vert[new_entry][i]==nshare) {
          order_other_ref_recur(i, level+1, order, on_path, num_child, has_out,
                                share_vert, max_share, solved);
        }
      }
    }

    on_path[new_entry] = 0;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int find_inout(int level, int num_child, int *num_vert, int *vert1,
               int *vertices, int *in_vertex, int *out_vertex, int *order)
{
/*
 * Function to find in and out vertices.
 * On first call, the first in_vertex and last out_vertex should already be set.
 * level should be 0 in the first call.
 */
int i, j;                       /* loop counters */
int solved;                     /* found a solution */

  if (level == num_child-1) {

  /*
   * Last element.  Success if the in vertex is not the last out
   */

    if (in_vertex[order[level]] == out_vertex[order[level]])
      solved = 0;
    else
      solved = 1;

  }
  else {

  /*
   * Not last element.
   * Try each vertex that is not the in vertex, and if the next element in
   * the path shares that vertex, move on to the next element
   */

    solved = 0;
    for (i=0; i<num_vert[order[level]] && !solved; i++) {
      if (vertices[vert1[order[level]]+i] != in_vertex[order[level]]) {
        for (j=0; j<num_vert[order[level+1]] && !solved; j++) {
          if (vertices[vert1[order[level+1]]+j] == vertices[vert1[order[level]]+i]) {
            out_vertex[order[level]]  = vertices[vert1[order[level]]+i];
            in_vertex[order[level+1]] = vertices[vert1[order[level]]+i];
            solved = find_inout(level+1, num_child, num_vert, vert1, vertices,
                                in_vertex, out_vertex, order);
          }
        }
      }
    }

  }

  return(solved);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int alloc_reftree_nodes(LB *lb, LB_REFTREE **node, int num_node,
                               int *num_vert)

{
/*
 *  Function to allocate num_node refinement tree nodes
 */

/*
 *  A pointer to the first allocated node is returned in node.
 *  num_vert is input to indicate the number of vertices to allocate for
 *  the element corresponding to each node.
 */

ZOLTAN_ID_PTR gids;     /* pointer to memory for GIDs */
ZOLTAN_ID_PTR lids;     /* pointer to memory for LIDs */
float *float_mem;   /* pointer to memory for floats */
int *int_mem;       /* pointer to memory for ints */
int sum_vert;       /* sum of num_vert */
int wdim;           /* dimension of object weights */
int i;              /* loop counter */

char *yo = "alloc_reftree_nodes";

  if (lb->Obj_Weight_Dim == 0) {
    wdim = 1;
  } else {
    wdim = lb->Obj_Weight_Dim;
  }

/* compute sum of num_vert */

  sum_vert = 0;
  for (i=0; i<num_node; i++) sum_vert = sum_vert + num_vert[i];

/* allocate the structures themselves */

  *node = (LB_REFTREE *) ZOLTAN_MALLOC(num_node*sizeof(LB_REFTREE));

/* allocate memory to be used within the structures */

  gids = ZOLTAN_ZOLTAN_MALLOC_GID_ARRAY(lb, num_node);
  lids = ZOLTAN_ZOLTAN_MALLOC_LID_ARRAY(lb, num_node);
  float_mem = (float *) ZOLTAN_MALLOC(3*wdim*num_node*sizeof(float));
  int_mem   = (int   *) ZOLTAN_MALLOC(sum_vert*sizeof(int));

  if (node == NULL || gids == NULL || lids == NULL || float_mem == NULL ||
      int_mem == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&gids);
    ZOLTAN_FREE(&lids);
    ZOLTAN_FREE(&float_mem);
    ZOLTAN_FREE(&int_mem);
    ZOLTAN_FREE(&node);
    ZOLTAN_LB_TRACE_EXIT(lb, yo);
    return(ZOLTAN_MEMERR);
  }

/* divide the memory up among the nodes */

  for (i=0; i<num_node; i++) {
    (*node)[i].global_id = gids;
    gids += lb->Num_GID;
    (*node)[i].local_id = lids;
    lids += lb->Num_LID;
    (*node)[i].weight = float_mem;
    (*node)[i].summed_weight = float_mem+wdim;
    (*node)[i].my_sum_weight = float_mem+2*wdim;
    float_mem += 3*wdim;
    (*node)[i].vertices = int_mem;
    int_mem += num_vert[i];
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/

void free_reftree_nodes(LB_REFTREE **node)

{
/*
 *  Function to free memory of one or more refinement tree nodes.
 *  node should be a pointer returned by alloc_reftree_nodes; all nodes
 *  allocated by that call are freed.
 */

  ZOLTAN_FREE(&((*node)->global_id));
  ZOLTAN_FREE(&((*node)->local_id));
  ZOLTAN_FREE(&((*node)->weight));
  ZOLTAN_FREE(&((*node)->vertices));
  ZOLTAN_FREE(node);

}

/*****************************************************************************/

void LB_Reftree_Free_Structure(LB *lb)

{
/*
 *  Function to free all the memory of a refinement tree
 */
struct LB_reftree_data_struct *reftree_data; /* data structure from lb */
LB_REFTREE *root;                            /* Root of the refinement tree */
struct LB_reftree_hash_node **hashtab;       /* hash table */
int hashsize;                                /* dimension of hash table */
int i;                                       /* loop counter */

  reftree_data = (struct LB_reftree_data_struct *)lb->Data_Structure;

  root = reftree_data->reftree_root;

  if (root != NULL) {

  /*
   * Make all the children of the root be leaves, recursively
   */

    if (root->children != NULL) {
      for (i=0; i<root->num_child; i++)
        LB_Reftree_Free_Subtree(lb, &(root->children[i]));
    }

  /*
   * Free the memory used by the children, making root a leaf
   */

      free_reftree_nodes(&(root->children));

  /*
   * Free the root
   */

    free_reftree_nodes(&root);
  }

  /*
   * Free the memory in the hash table
   */

  hashtab  = reftree_data->hash_table;
  hashsize = reftree_data->hash_table_size;

  if (hashtab != NULL) {
    LB_Reftree_Clear_Hash_Table(hashtab,hashsize);
    ZOLTAN_FREE(&hashtab);
  }

  ZOLTAN_FREE(&(lb->Data_Structure));

}

static void LB_Reftree_Free_Subtree(LB *lb, LB_REFTREE *subroot)

{
/*
 *  Function to free the memory of a subtree.  Upon return, subroot is a leaf.
 */
int i;   /* loop counter */
struct LB_reftree_data_struct *reftree_data; /* data structure from lb */

  if (subroot != NULL) {

    reftree_data = (struct LB_reftree_data_struct *)lb->Data_Structure;

  /*
   * Turn all the children into leaves and remove them from the hash table
   */

    if (subroot->children != NULL) {
      for (i=0; i<subroot->num_child; i++) {
        LB_Reftree_Free_Subtree(lb,&(subroot->children[i]));
        LB_Reftree_Hash_Remove(lb,&(subroot->children[i]),
                               reftree_data->hash_table,
                               reftree_data->hash_table_size);
      }

  /*
   * Free the memory used by the children, making subroot a leaf
   */

      free_reftree_nodes(&(subroot->children));
      subroot->num_child = 0;
    }

  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int LB_Reftree_Reinit_Coarse(LB *lb)

{
/*
 *  Function to reestablish which coarse grid elements are known to this proc
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

char *yo = "LB_Reftree_Reinit_Coarse";
LB_REFTREE *root;     /* Root of the refinement tree */
struct LB_reftree_hash_node **hashtab; /* hash table */
int hashsize;         /* dimension of hash table */
int i, j;             /* loop counter */
ZOLTAN_ID_PTR local_gids; /* coarse element Global IDs from user */
ZOLTAN_ID_PTR local_lids; /* coarse element Local IDs from user */
ZOLTAN_ID_PTR lid;        /* temporary coarse element Local ID; used to pass
                         NULL to query functions when NUM_LID_ENTRIES=0 */
int *assigned;        /* 1 if the element is assigned to this proc */
int *num_vert;        /* number of vertices for each coarse element */
int *vertices;        /* vertices for the coarse elements */
int *in_vertex;       /* "in" vertex for each coarse element */
int *out_vertex;      /* "out" vertex for each coarse element */
ZOLTAN_ID_PTR slocal_gids;/* coarse element Global IDs from user */
ZOLTAN_ID_PTR slocal_lids;/* coarse element Local IDs from user */
ZOLTAN_ID_PTR plocal_gids;/* previous coarse element Global IDs from user */
ZOLTAN_ID_PTR plocal_lids;/* previous coarse element Local IDs from user */
int sassigned;        /* 1 if the element is assigned to this proc */
int snum_vert;        /* number of vertices for a coarse element */
int sin_vertex;       /* "in" vertex for a coarse element */
int sout_vertex;      /* "out" vertex for a coarse element */
int in_order;         /* 1 if user is supplying order of the elements */
int num_obj;          /* number of coarse objects known to this proc */
int ierr;             /* error flag */
LB_REFTREE *tree_node;/* pointer to an initial grid element in the tree */
int final_ierr;       /* error code returned */
int sum_vert;         /* running total of number of vertices */
int found;            /* flag for another coarse grid element */
int num_gid_entries = lb->Num_GID;  /* number of array entries in a global ID */
int num_lid_entries = lb->Num_LID;  /* number of array entries in a local ID */

  root = ((struct LB_reftree_data_struct *)lb->Data_Structure)->reftree_root;
  hashtab  = ((struct LB_reftree_data_struct *)lb->Data_Structure)->hash_table;
  hashsize = ((struct LB_reftree_data_struct *)lb->Data_Structure)->hash_table_size;
  final_ierr = ZOLTAN_OK;

  /*
   * Mark all coarse elements as unknown
   */

  for (i=0; i<root->num_child; i++) {
    ((root->children)[i]).num_vertex = -1;
  }

  /*
   * Get the coarse grid objects and update the vertices, whether the element
   * is assigned to this processor, weight and, if not already set, in/out vert
   */

  if (lb->Get_Coarse_Obj_List != NULL) {

  /*
   * Get objects via list
   */

    num_obj = lb->Get_Num_Coarse_Obj(lb->Get_Num_Coarse_Obj_Data, &ierr);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                     "Error returned from user function Get_Num_Coarse_Obj.");
      return(ierr);
    }

    if (num_obj > 0) {
      local_gids = ZOLTAN_ZOLTAN_MALLOC_GID_ARRAY(lb, num_obj);
      local_lids = ZOLTAN_ZOLTAN_MALLOC_LID_ARRAY(lb, num_obj);
      assigned   = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
      num_vert   = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
      vertices   = (int *) ZOLTAN_MALLOC(MAXVERT*num_obj*sizeof(int));
      in_vertex  = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));
      out_vertex = (int *) ZOLTAN_MALLOC(num_obj*sizeof(int));

      if (local_gids == NULL || (num_lid_entries > 0 && local_lids == NULL) ||
          assigned   == NULL ||
          num_vert   == NULL || vertices   == NULL || in_vertex == NULL ||
          out_vertex == NULL) {
        ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
        ZOLTAN_FREE(&local_gids);
        ZOLTAN_FREE(&local_lids);
        ZOLTAN_FREE(&assigned);
        ZOLTAN_FREE(&num_vert);
        ZOLTAN_FREE(&vertices);
        ZOLTAN_FREE(&in_vertex);
        ZOLTAN_FREE(&out_vertex);
        return(ZOLTAN_MEMERR);
      }

      lb->Get_Coarse_Obj_List(lb->Get_Coarse_Obj_List_Data, 
                              num_gid_entries, num_lid_entries,
                              local_gids, local_lids, 
                              assigned, num_vert, vertices,
                              &in_order, in_vertex, out_vertex, &ierr);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                      "Error returned from user function Get_Coarse_Obj_List.");
        ZOLTAN_FREE(&local_gids);
        ZOLTAN_FREE(&local_lids);
        ZOLTAN_FREE(&assigned);
        ZOLTAN_FREE(&num_vert);
        ZOLTAN_FREE(&vertices);
        ZOLTAN_FREE(&in_vertex);
        ZOLTAN_FREE(&out_vertex);
        return(ierr);
      }

      sum_vert = 0;
      for (i=0; i<num_obj; i++) {

        tree_node = LB_Reftree_hash_lookup(lb, hashtab,
                                           &(local_gids[i*num_gid_entries]),
                                           hashsize);
        if (tree_node == NULL) {
          ZOLTAN_PRINT_WARN(lb->Proc, yo, "coarse grid element not"
                                      " previously seen.");
          final_ierr = ZOLTAN_WARN;
        }
        else {
/* can just reassign num_vertex instead of doing a complicated reallocation
   because coarse grid objects are allocated with MAXVERT */
          tree_node->num_vertex = num_vert[i];
          for (j=0; j<num_vert[i]; j++)
            tree_node->vertices[j] = vertices[sum_vert+j];
          if (num_vert[i] > 0) sum_vert += num_vert[i];

          tree_node->assigned_to_me = assigned[i];
/* TEMP if not provided in_order, then in/out are not returned and must be
        determined */
          if (tree_node->in_vertex == 0) tree_node->in_vertex = in_vertex[i];
          if (tree_node->out_vertex == 0) tree_node->out_vertex = out_vertex[i];
          lid = (num_lid_entries ? &(local_lids[i*num_lid_entries]) : NULL);
          lb->Get_Child_Weight(lb->Get_Child_Weight_Data, 
                               num_gid_entries, num_lid_entries,
                               &(local_gids[i*num_gid_entries]),
                               lid,
                               lb->Obj_Weight_Dim,
                               tree_node->weight, &ierr);
        }
      }
      ZOLTAN_FREE(&local_gids);
      ZOLTAN_FREE(&local_lids);
      ZOLTAN_FREE(&assigned);
      ZOLTAN_FREE(&num_vert);
      ZOLTAN_FREE(&vertices);
      ZOLTAN_FREE(&in_vertex);
      ZOLTAN_FREE(&out_vertex);
    }

  }
  else {

  /*
   * Get objects via first/next
   */

    slocal_gids = ZOLTAN_ZOLTAN_MALLOC_GID(lb);
    slocal_lids = ZOLTAN_ZOLTAN_MALLOC_LID(lb);
    plocal_gids = ZOLTAN_ZOLTAN_MALLOC_GID(lb);
    plocal_lids = ZOLTAN_ZOLTAN_MALLOC_LID(lb);
    vertices = (int *) ZOLTAN_MALLOC(MAXVERT*sizeof(int));
    if (slocal_gids == NULL || (num_lid_entries > 0 && slocal_lids == NULL) || 
        plocal_gids == NULL || (num_lid_entries > 0 && plocal_lids == NULL) || 
        vertices == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&slocal_gids);
      ZOLTAN_FREE(&slocal_lids);
      ZOLTAN_FREE(&vertices);
      return(ZOLTAN_MEMERR);
    }

    found = lb->Get_First_Coarse_Obj(lb->Get_First_Coarse_Obj_Data,
                                     num_gid_entries, num_lid_entries,
                                     slocal_gids, slocal_lids, &sassigned,
                                     &snum_vert, vertices, &in_order,
                                     &sin_vertex, &sout_vertex, &ierr);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                     "Error returned from user function Get_First_Coarse_Obj.");
      ZOLTAN_FREE(&slocal_gids);
      ZOLTAN_FREE(&slocal_lids);
      ZOLTAN_FREE(&vertices);
      return(ierr);
    }
    while (found) {
      tree_node = LB_Reftree_hash_lookup(lb, hashtab,slocal_gids,hashsize);
      if (tree_node == NULL) {
        ZOLTAN_PRINT_WARN(lb->Proc, yo, "coarse grid element not"
                                    " previously seen.");
        final_ierr = ZOLTAN_WARN;
      }
      else {
/* can just reassign num_vertex instead of doing a complicated reallocation
   because coarse grid objects are allocated with MAXVERT */
        tree_node->num_vertex = snum_vert;
        for (j=0; j<snum_vert; j++)
          tree_node->vertices[j] = vertices[j];
  
        tree_node->assigned_to_me = sassigned;
/* TEMP if not provided in_order, then in/out are not returned and must be
        determined */
        if (tree_node->in_vertex == 0) tree_node->in_vertex = sin_vertex;
        if (tree_node->out_vertex == 0) tree_node->out_vertex = sout_vertex;
        lb->Get_Child_Weight(lb->Get_Child_Weight_Data, 
                             num_gid_entries, num_lid_entries,
                             slocal_gids, slocal_lids, lb->Obj_Weight_Dim,
                             tree_node->weight, &ierr);
      }

      ZOLTAN_LB_SET_GID(lb, plocal_gids, slocal_gids);
      ZOLTAN_LB_SET_LID(lb, plocal_lids, slocal_lids);
      found = lb->Get_Next_Coarse_Obj(lb->Get_Next_Coarse_Obj_Data,
                                      num_gid_entries, num_lid_entries,
                                      plocal_gids, plocal_lids,
                                      slocal_gids, slocal_lids, &sassigned,
                                      &snum_vert, vertices,
                                      &sin_vertex, &sout_vertex, &ierr);
    }
    ZOLTAN_FREE(&slocal_gids);
    ZOLTAN_FREE(&slocal_lids);
    ZOLTAN_FREE(&plocal_gids);
    ZOLTAN_FREE(&plocal_lids);
  }
  return(final_ierr);
}

void LB_Reftree_Print(LB *lb, LB_REFTREE *subroot, int level)
{
/*
 * Print the refinement tree, for debugging
 */

  int i, me;

  if (subroot == NULL) return;

  me = lb->Proc;
  printf("\n");
  printf("[%d] refinement tree node with local id ", me);
  ZOLTAN_LB_PRINT_LID(lb, subroot->local_id);
  printf(" on level %d\n", level);
  printf("[%d]   Global ID ",me);
  ZOLTAN_LB_PRINT_GID(lb, subroot->global_id);
  printf("\n");
  printf("[%d]   first weight %f\n",me,subroot->weight[0]);
  printf("[%d]   first summed weight %f\n",me,subroot->summed_weight[0]);
  printf("[%d]   first my_sum weight %f\n",me,subroot->my_sum_weight[0]);
  printf("[%d]   number of vertices %d\n",me,subroot->num_vertex);
  printf("[%d]   vertices ",me);
  for (i=0; i<subroot->num_vertex; i++) printf("%d ",subroot->vertices[i]);
  printf("\n");
  printf("[%d]   in and out vertices %d %d\n",me,subroot->in_vertex,subroot->out_vertex);
  printf("[%d]   assigned_to_me %d\n",me,subroot->assigned_to_me);
  printf("[%d]   partition %d\n",me,subroot->partition);
  printf("[%d]   number of children %d \n",me,subroot->num_child);
  printf("[%d]   children follow.\n",me);
  for (i=0; i<subroot->num_child; i++)
    LB_Reftree_Print(lb,&(subroot->children[i]),level+1);
}

/* TEMP child_order */

static void get_child_order_recur(LB *lb, LB_REFTREE *subroot, int *isub, int *order)
{

  /*
   * adds the children to order and recursively continues down the tree
   */

int i;

  /*
   * if no children, done with this branch
   */

  if (subroot->num_child == 0) {
    return;
  }

  /*
   * add the subroot and children to order
   */

  order[*isub] = subroot->local_id[0];
  for (i=0; i<subroot->num_child; i++) {
    order[*isub+i+1] = (subroot->children[i]).local_id[0];
  }
  *isub = *isub + subroot->num_child + 1;

  /*
   * traverse the children
   */

  for (i=0; i<subroot->num_child; i++) {
    get_child_order_recur(lb, &(subroot->children[i]), isub, order);
  }
}

void LB_Get_Child_Order(LB *lb, int *order, int *ierr)
{
/*
 * Return the order of the children in the refinement tree.
 * Upon return, order contains LIDs assumed to be an integer.  It contains
 * sets of entries consisting of the LID of an element followed by the
 * LIDs of the children in the order determined by the reftree code.
 * order should be allocated to the correct size by the caller.
 * This is a hack, will be removed in the future, and should not be publicized.
 */

char *yo = "LB_Get_Child_Order";
int isub;
LB_REFTREE *root;

  *ierr = ZOLTAN_OK;

  /*
   * initialize the tree, if not already done
   */

  if (lb->Data_Structure == NULL) {
    *ierr = LB_Reftree_Init(lb);
    if (*ierr==ZOLTAN_FATAL || *ierr==ZOLTAN_MEMERR) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo,
                     "Error returned by LB_Reftree_Init.");
      return;
    }
  }

  /*
   * build the refinement tree
   */

  *ierr = LB_Reftree_Build(lb);
  if (*ierr==ZOLTAN_FATAL || *ierr==ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo,
                   "Error returned by LB_Reftree_Build.");
    return;
  }

  /*
   * traverse the tree to find the child order
   */

  root = ((struct LB_reftree_data_struct *)lb->Data_Structure)->reftree_root;
  isub = 0;
  get_child_order_recur(lb,root,&isub,order);

  /*
   * delete the tree, except for the first level (initial coarse grid)
   */

}

/* end TEMP child_order */
