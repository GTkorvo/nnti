/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#include <stdio.h>
#include "lb_const.h"
#include "reftree.h"
#include "all_allo_const.h"
#include "params_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Prototypes for functions internal to this file */
static int Zoltan_Reftree_Sum_Weights(LB *lb);

int Zoltan_Reftree_Sum_Weights_gather(LB *lb);
int Zoltan_Reftree_Sum_Weights_pairs(LB *lb);
int Zoltan_Reftree_Sum_Weights_bcast(LB *lb);

static void Zoltan_Reftree_Sum_My_Weights(LB *lb, ZOLTAN_REFTREE *subroot, 
       int *count, int wdim);
static void Zoltan_Reftree_Sum_All_Weights(LB *lb, ZOLTAN_REFTREE *subroot, int wdim);
static void Zoltan_Reftree_List_Other_Leaves(LB *lb, ZOLTAN_REFTREE *subroot, 
       ZOLTAN_ID_PTR list, int *count);
static int Zoltan_Reftree_Partition(LB *lb, int *num_export, 
       ZOLTAN_ID_PTR *export_global_ids, ZOLTAN_ID_PTR *export_local_ids, 
       int **export_procs);
static void Zoltan_Reftree_Part_Recursive(LB *lb, ZOLTAN_REFTREE *subroot, int *part,
       float *current_size, int *num_exp, float *cutoff,
       int num_part, float partition_size, float eps);
static void Zoltan_Reftree_Mark_and_Count(ZOLTAN_REFTREE *subroot, int part, 
       int *num_exp);
static void Zoltan_Reftree_Export_Lists(LB *lb, ZOLTAN_REFTREE *subroot, 
       int *num_export, ZOLTAN_ID_PTR *export_global_ids,
       ZOLTAN_ID_PTR *export_local_ids, int **export_procs);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Reftree_Part(

  LB *lb,                       /* The load-balancing structure */
  int *num_import,              /* Not computed, set to -1 */
  ZOLTAN_ID_PTR *import_global_ids, /* Not computed */
  ZOLTAN_ID_PTR *import_local_ids,  /* Not computed */
  int **import_procs,           /* Not computed */
  int *num_export,              /* Number of objects to be exported */
  ZOLTAN_ID_PTR *export_global_ids, /* global ids of objects to be exported */
  ZOLTAN_ID_PTR *export_local_ids,  /* local  ids of objects to be exported */
  int **export_procs            /* list of processors to export to */
)
{
char *yo = "Zoltan_Reftree_Part";
int ierr;       /* error code returned by called routines */
int final_ierr; /* error code returned by this routine */
double time0, time1, time2, time3, time4;

  /* Initializations in case of early exit. */
  *num_export = -1;
  *num_import = -1;
  final_ierr = ZOLTAN_OK;

  /*
   * initialize the tree (first call only)
   */

  if (lb->Data_Structure == NULL) {
    if (lb->Debug_Level >= ZOLTAN_DEBUG_ATIME) time0 = Zoltan_Time(lb->Timer);
    ierr = Zoltan_Reftree_Init(lb);
    if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                     "Error returned by Zoltan_Reftree_Init.");
      return(ierr);
    }
    if (ierr==ZOLTAN_WARN) final_ierr = ZOLTAN_WARN;
    if (lb->Debug_Level >= ZOLTAN_DEBUG_ATIME) time1 = Zoltan_Time(lb->Timer);
  } else {
    if (lb->Debug_Level >= ZOLTAN_DEBUG_ATIME) {
      time1 = Zoltan_Time(lb->Timer);
      time0 = time1 + 1.0;
    }
  }

  /*
   * build the refinement tree
   */

  ierr = Zoltan_Reftree_Build(lb);
  if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned by Zoltan_Reftree_Build.");
    return(ierr);
  }
  if (ierr==ZOLTAN_WARN) final_ierr = ZOLTAN_WARN;
  if (lb->Debug_Level >= ZOLTAN_DEBUG_ATIME) time2 = Zoltan_Time(lb->Timer);

  /*
   * sum the weights in the tree
   */

  ierr = Zoltan_Reftree_Sum_Weights(lb);
  if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned by Zoltan_Reftree_Sum_Weights.");
    return(ierr);
  }
  if (ierr==ZOLTAN_WARN) final_ierr = ZOLTAN_WARN;
  if (lb->Debug_Level >= ZOLTAN_DEBUG_ATIME) time3 = Zoltan_Time(lb->Timer);

  /*
   * determine the new partition
   */

  ierr = Zoltan_Reftree_Partition(lb, num_export, export_global_ids,
                              export_local_ids, export_procs);
  if (ierr==ZOLTAN_FATAL || ierr==ZOLTAN_MEMERR) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned by Zoltan_Reftree_Partition.");
    return(ierr);
  }
  if (ierr==ZOLTAN_WARN) final_ierr = ZOLTAN_WARN;
  if (lb->Debug_Level >= ZOLTAN_DEBUG_ATIME) time4 = Zoltan_Time(lb->Timer);

  if (lb->Debug_Level >= ZOLTAN_DEBUG_ATIME) {
    if (time0 <= time1) {
      Zoltan_Print_Stats(lb->Communicator, lb->Debug_Proc, time1-time0,
                     "REFTREE Time to initialize :");
    }
    Zoltan_Print_Stats(lb->Communicator, lb->Debug_Proc, time2-time1, 
                   "REFTREE Time to build tree :");
    Zoltan_Print_Stats(lb->Communicator, lb->Debug_Proc, time3-time2,
                   "REFTREE Time to sum weights:");
    Zoltan_Print_Stats(lb->Communicator, lb->Debug_Proc, time4-time3,
                   "REFTREE Time to partition  :");
  }

  return(final_ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int Zoltan_Reftree_Sum_Weights(LB *lb)

{
/*
 * Wrapper function for summing the weights.  This just calls one of
 * three functions for summing the weights, one of which uses MPI collective
 * communication routines, another of which uses broadcasts
 * from each processor to all other processors, and the other of which
 * uses O(log p) pairwise communication steps.
 * TEMP eventually determine what is the best approach and get rid of the
 * others, or maybe keep two of them with a switch depending on the number
 * of processors or something.
 */

/* currently using the collective communication version */
  return(Zoltan_Reftree_Sum_Weights_gather(lb));

/* examples of calling the others
   return(Zoltan_Reftree_Sum_Weights_pairs(lb));
   return(Zoltan_Reftree_Sum_Weights_bcast(lb));
*/
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Reftree_Sum_Weights_gather(LB *lb)

{
/*
 * Function to sum the weights in the refinement tree.  On input the
 * refinement tree should be valid and have weight set.  On output the
 * values in summed_weight at each node is the sum of the weights in the
 * subtree with that node as root.
 * This function also sets assigned_to_me for interior nodes to be
 * 1 if the entire subtree is assigned to this processor
 * 0 if none of the subtree is assigned to this processor
 * -1 if some of the subtree is assigned to this processor
 *
 * This version uses MPI collective communication routines.
 */
char *yo = "Zoltan_Reftree_Sum_Weights_gather";
ZOLTAN_REFTREE *root;         /* Root of the refinement tree */
int wdim;                 /* Dimension of the weight array */
int i,j;                  /* loop counters */
int count;                /* counter */
ZOLTAN_ID_PTR leaf_list = NULL;      
                          /* leaves for which some proc requests weight */
ZOLTAN_ID_PTR all_leaflist = NULL;   
                          /* leaf_list from all processors */
int reqsize;              /* length of leaf_list */
int *reqsize_all;         /* reqsize from all processors */
int sum_reqsize;          /* sum of all reqsize */
int *displs;              /* running sum of all reqsize */
int my_start;             /* position in leaf_list of this proc's list */
int nproc;                /* number of processors */
ZOLTAN_REFTREE *node;         /* a node in the refinement tree */
struct Zoltan_Reftree_hash_node **hashtab; /* hash table */
int hashsize;             /* dimension of hash table */
float *send_float;        /* sending message of floats */
float *req_weights;       /* the requested weights */
int num_gid_entries = lb->Num_GID; /* Number of array entries in a global ID */

  /*
   * set the root and hash table
   */

  root = ((struct Zoltan_Reftree_data_struct *)lb->Data_Structure)->reftree_root;
  if (root == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Refinement tree not defined.");
    return(ZOLTAN_FATAL);
  }
  hashtab  = ((struct Zoltan_Reftree_data_struct *)lb->Data_Structure)->hash_table;
  hashsize = ((struct Zoltan_Reftree_data_struct *)lb->Data_Structure)->hash_table_size;

  /*
   * Determine the dimension of the weight array
   */

  if (lb->Obj_Weight_Dim == 0) {
    wdim = 1;
  } else {
    wdim = lb->Obj_Weight_Dim;
  }

  /*
   * In the first pass, sum the weights of the nodes that are assigned to
   * this processor, and count the leaves that are not.
   */

  count = 0;
  for (i=0; i<root->num_child; i++) {
    Zoltan_Reftree_Sum_My_Weights(lb,&(root->children[i]),&count,wdim);
  }
  root->assigned_to_me = -1;

  /*
   * Make a list of the leaves that are not assigned to this processor
   */

  if (count == 0)
    leaf_list = ZOLTAN_MALLOC_GID(lb);
  else
    leaf_list = ZOLTAN_MALLOC_GID_ARRAY(lb, count);
  if (leaf_list == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }

  count = 0;
  Zoltan_Reftree_List_Other_Leaves(lb, root,leaf_list,&count);

  /*
   * Get the unknown leaf weights from other processors.
   */

  nproc = lb->Num_Proc;
  reqsize = count;

  /*
   * Build a list of all processor's request list by concatinating them in
   * the order of the processor ranks
   */

  /*
   * Determine the request size of all processors
   */

  reqsize_all = (int *)ZOLTAN_MALLOC(nproc*sizeof(int));
  displs = (int *)ZOLTAN_MALLOC(nproc*sizeof(int));
  if (reqsize_all == NULL || displs == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&displs);
    ZOLTAN_FREE(&reqsize_all);
    ZOLTAN_FREE(&leaf_list);
    return(ZOLTAN_MEMERR);
  }

  MPI_Allgather((void *)&reqsize,1,MPI_INT,(void *)reqsize_all,1,MPI_INT,
                lb->Communicator);
  displs[0] = 0;
  for (i=1; i<nproc; i++) displs[i] = displs[i-1]+reqsize_all[i-1];
  sum_reqsize = displs[nproc-1] + reqsize_all[nproc-1];
  my_start = displs[lb->Proc];

  /*
   * If sum_reqsize is 0, nothing needs to be communciated
   */

  if (sum_reqsize == 0) {
    ZOLTAN_FREE(&displs);
    ZOLTAN_FREE(&reqsize_all);
    ZOLTAN_FREE(&leaf_list);
  }
  else {

  /*
   * Gather the request list from all processors
   */

    all_leaflist = ZOLTAN_MALLOC_GID_ARRAY(lb, sum_reqsize);
    if (all_leaflist == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&all_leaflist);
      ZOLTAN_FREE(&displs);
      ZOLTAN_FREE(&reqsize_all);
      ZOLTAN_FREE(&leaf_list);
      return(ZOLTAN_MEMERR);
    }

    /* KDDKDD Changed MPI_BYTE to ZOLTAN_ID_MPI_TYPE  */

    /* Account for number of array entries in an ID. */
    for (i=0; i<nproc; i++) {
      reqsize_all[i] = reqsize_all[i]*num_gid_entries;
      displs[i] = displs[i]*num_gid_entries;
    }

    MPI_Allgatherv((void *)leaf_list,reqsize*num_gid_entries,ZOLTAN_ID_MPI_TYPE,
                   (void *)all_leaflist,reqsize_all,displs,ZOLTAN_ID_MPI_TYPE,
                   lb->Communicator);

    ZOLTAN_FREE(&displs);
    ZOLTAN_FREE(&leaf_list);

    for (i=0; i<nproc; i++) reqsize_all[i] = reqsize_all[i]/num_gid_entries;

  /* 
   * Create a list with the partial sums this processor has
   */

    send_float = (float *) ZOLTAN_MALLOC(sizeof(float)*wdim*sum_reqsize);
    if (send_float == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&send_float);
      ZOLTAN_FREE(&all_leaflist);
      ZOLTAN_FREE(&reqsize_all);
      return(ZOLTAN_MEMERR);
    }

    for (i=0; i<sum_reqsize; i++) {
      node = Zoltan_Reftree_hash_lookup(lb, hashtab,
                                    &(all_leaflist[i*num_gid_entries]),
                                    hashsize);
      if (node == NULL)
         for (j=0; j<wdim; j++) send_float[i*wdim+j] = 0.0;
      else
         for (j=0; j<wdim; j++) send_float[i*wdim+j] = node->my_sum_weight[j];
    }

  /*
   * Sum the weights over all the processors
   */

    if (reqsize == 0)
      req_weights = (float *) ZOLTAN_MALLOC(sizeof(float)*wdim);
    else
      req_weights = (float *) ZOLTAN_MALLOC(sizeof(float)*wdim*reqsize);
    if (req_weights == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&req_weights);
      ZOLTAN_FREE(&send_float);
      ZOLTAN_FREE(&all_leaflist);
      ZOLTAN_FREE(&reqsize_all);
      return(ZOLTAN_MEMERR);
    }

    MPI_Reduce_scatter((void *)send_float, (void *)req_weights, reqsize_all,
                       MPI_FLOAT, MPI_SUM, lb->Communicator);

    ZOLTAN_FREE(&send_float);
    ZOLTAN_FREE(&reqsize_all);

  /*
   * Set the weights this processor requested
   */

    for (i=0; i<count; i++) {
      node = Zoltan_Reftree_hash_lookup(lb, hashtab,
                                  &(all_leaflist[(i+my_start)*num_gid_entries]),
                                  hashsize);
      for (j=0; j<wdim; j++) node->summed_weight[j] = req_weights[i*wdim+j];
    }

    ZOLTAN_FREE(&req_weights);
    ZOLTAN_FREE(&all_leaflist);
  }

  /*
   * All the leaves now have summed_weight set.
   * Sum the weights throughout the tree.
   */

  Zoltan_Reftree_Sum_All_Weights(lb,root,wdim);

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Reftree_Sum_Weights_pairs(LB *lb)

{
/*
 * Function to sum the weights in the refinement tree.  On input the
 * refinement tree should be valid and have weight set.  On output the
 * values in summed_weight at each node is the sum of the weights in the
 * subtree with that node as root.
 * This function also sets assigned_to_me for interior nodes to be
 * 1 if the entire subtree is assigned to this processor
 * 0 if none of the subtree is assigned to this processor
 * -1 if some of the subtree is assigned to this processor
 *
 * This version uses O(log p) pairwise communication steps.
 * TEMP Restricted to power-of-2 number of processors.
 */
char *yo = "Zoltan_Reftree_Sum_Weights_pairs";
ZOLTAN_REFTREE *root;         /* Root of the refinement tree */
int wdim;                 /* Dimension of the weight array */
int i,j,comm_loop;        /* loop counters */
int count;                /* counter */
ZOLTAN_ID_PTR leaf_list = NULL; 
                          /* leaves for which some proc requests weight */
int reqsize;              /* length of leaf_list */
ZOLTAN_ID_PTR newlist = NULL; /* building concatinated leaf_list of all procs */
int newsize;              /* new length of leaf_list */
int my_start;             /* position in leaf_list of this proc's list */
int nproc;                /* number of processors */
int log_nproc;            /* base 2 log of number of processors */
int myproc;               /* this processor's processor number */
int other_proc;           /* processor being communicated with */
MPI_Status mess_status;   /* status of an MPI message */
ZOLTAN_REFTREE *node;         /* a node in the refinement tree */
struct Zoltan_Reftree_hash_node **hashtab; /* hash table */
int hashsize;             /* dimension of hash table */
float *send_float;        /* sending message of floats */
float *req_weights;       /* the requested weights */
int gid_off;              /* offset into array of global IDs. */
int num_gid_entries = lb->Num_GID;  /* Number of array entries in a global ID */

  /*
   * set the root and hash table
   */

  root = ((struct Zoltan_Reftree_data_struct *)lb->Data_Structure)->reftree_root;
  if (root == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Refinement tree not defined.");
    return(ZOLTAN_FATAL);
  }
  hashtab  = ((struct Zoltan_Reftree_data_struct *)lb->Data_Structure)->hash_table;
  hashsize = ((struct Zoltan_Reftree_data_struct *)lb->Data_Structure)->hash_table_size;

  /*
   * Determine the dimension of the weight array
   */

  if (lb->Obj_Weight_Dim == 0) {
    wdim = 1;
  } else {
    wdim = lb->Obj_Weight_Dim;
  }

  /*
   * In the first pass, sum the weights of the nodes that are assigned to
   * this processor, and count the leaves that are not.
   */

  count = 0;
  for (i=0; i<root->num_child; i++) {
    Zoltan_Reftree_Sum_My_Weights(lb,&(root->children[i]),&count,wdim);
  }
  root->assigned_to_me = -1;

  /*
   * Make a list of the leaves that are not assigned to this processor
   */

  if (count == 0)
    leaf_list = ZOLTAN_MALLOC_GID(lb);
  else
    leaf_list = ZOLTAN_MALLOC_GID_ARRAY(lb, count);
  if (leaf_list == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }

  count = 0;
  Zoltan_Reftree_List_Other_Leaves(lb, root,leaf_list,&count);

  /*
   * Get the unknown leaf weights from other processors.
   */

  nproc = lb->Num_Proc;
  myproc = lb->Proc;
  log_nproc = -1;
  i = nproc;
  while(i) {log_nproc++; i=i>>1;};

  /*
   * Build a list of all processor's request list by concatinating them in
   * the order of the processor ranks
   */

  reqsize = count;
  my_start = 0;

  for (comm_loop=0; comm_loop<log_nproc; comm_loop++) {

  /*
   * Exchange the leaf_list with the next processor
   */

    other_proc = myproc^(1<<comm_loop);

  /*
   * First exchange the size of the lists
   */

    MPI_Send((void *)&reqsize,1,MPI_INT,other_proc,100+comm_loop,
             lb->Communicator);
    MPI_Recv((void *)&newsize,1,MPI_INT,other_proc,100+comm_loop,
             lb->Communicator, &mess_status);

  /*
   * Allocate space for the concatinated list, copy the current list to it,
   * and free the current list
   */

    if (reqsize+newsize == 0)
      newlist = ZOLTAN_MALLOC_GID(lb);
    else
      newlist = ZOLTAN_MALLOC_GID_ARRAY(lb, (reqsize+newsize));
    if (newlist == NULL) {
      ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&leaf_list);
      return(ZOLTAN_MEMERR);
    }
    if (myproc < other_proc) {
      for (i=0; i<reqsize; i++) {
        gid_off = i * num_gid_entries;
        ZOLTAN_SET_GID(lb, &(newlist[gid_off]),&(leaf_list[gid_off]));
      }
    }
    else {
      for (i=0; i<reqsize; i++) {
        ZOLTAN_SET_GID(lb, &(newlist[(i+newsize)*num_gid_entries]),
                       &(leaf_list[i*num_gid_entries]));
      }
      my_start += newsize;
    }
    ZOLTAN_FREE(&leaf_list);
    leaf_list = newlist;

  /*
   * Finally, exchange the lists
   */
    /* KDDKDD Changed MPI_BYTE to ZOLTAN_ID_MPI_TYPE  */

    if (myproc < other_proc) {
      MPI_Send((void *)leaf_list,reqsize*num_gid_entries,ZOLTAN_ID_MPI_TYPE,
               other_proc, 200+comm_loop,lb->Communicator);
      MPI_Recv((void *)&(leaf_list[reqsize*num_gid_entries]),
               newsize*num_gid_entries,ZOLTAN_ID_MPI_TYPE,
               other_proc,200+comm_loop,lb->Communicator,&mess_status);
    }
    else {
      MPI_Send((void *)&(leaf_list[newsize*num_gid_entries]),
                reqsize*num_gid_entries,ZOLTAN_ID_MPI_TYPE,
               other_proc,200+comm_loop,lb->Communicator);
      MPI_Recv((void *)leaf_list,newsize*num_gid_entries,ZOLTAN_ID_MPI_TYPE,
               other_proc,200+comm_loop,lb->Communicator,&mess_status);
    }

    reqsize = reqsize + newsize;
  }

  /* 
   * Create a list with the partial sums this processor has
   */

  send_float = (float *) ZOLTAN_MALLOC(sizeof(float)*wdim*reqsize);
  if (send_float == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&leaf_list);
    return(ZOLTAN_MEMERR);
  }

  for (i=0; i<reqsize; i++) {
    node = Zoltan_Reftree_hash_lookup(lb, hashtab, &(leaf_list[i*num_gid_entries]),
                                  hashsize);
    if (node == NULL)
       for (j=0; j<wdim; j++) send_float[i*wdim+j] = 0.0;
    else
       for (j=0; j<wdim; j++) send_float[i*wdim+j] = node->my_sum_weight[j];
  }

  /*
   * Sum the weights over all the processors
   */

  req_weights = (float *) ZOLTAN_MALLOC(sizeof(float)*wdim*reqsize);
  if (req_weights == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&leaf_list);
    ZOLTAN_FREE(&send_float);
    return(ZOLTAN_MEMERR);
  }

/* TEMP perhaps it would be better to use MPI_Reduce_Scatter */
  MPI_Allreduce((void *)send_float, (void *)req_weights, reqsize, MPI_FLOAT,
                MPI_SUM, lb->Communicator);

  ZOLTAN_FREE(&send_float);

  /*
   * Set the weights this processor requested
   */

  for (i=0; i<count; i++) {
    node = Zoltan_Reftree_hash_lookup(lb, hashtab,
              &(leaf_list[(i+my_start)*num_gid_entries]),hashsize);
    for (j=0; j<wdim; j++) node->summed_weight[j] = req_weights[(i+my_start)*wdim+j];
  }

  ZOLTAN_FREE(&req_weights);
  ZOLTAN_FREE(&leaf_list);

  /*
   * All the leaves now have summed_weight set.
   * Sum the weights throughout the tree.
   */

  Zoltan_Reftree_Sum_All_Weights(lb,root,wdim);

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Reftree_Sum_Weights_bcast(LB *lb)

{
/*
 * Function to sum the weights in the refinement tree.  On input the
 * refinement tree should be valid and have weight set.  On output the
 * values in summed_weight at each node is the sum of the weights in the
 * subtree with that node as root.
 * This function also sets assigned_to_me for interior nodes to be
 * 1 if the entire subtree is assigned to this processor
 * 0 if none of the subtree is assigned to this processor
 * -1 if some of the subtree is assigned to this processor
 *
 * This version uses all-to-all communciation.
 */
char *yo = "Zoltan_Reftree_Sum_Weights_bcast";
ZOLTAN_REFTREE *root;         /* Root of the refinement tree */
int wdim;                 /* Dimension of the weight array */
int i,j,comm_loop,rproc;  /* loop counters */
int count;                /* counter */
ZOLTAN_ID_PTR leaf_list = NULL; 
                          /* list of the leaves not assigned to this proc */
int nproc;                /* number of processors */
int myproc;               /* this processor's processor number */
MPI_Status mess_status;   /* status of an MPI message */
int reqsize;              /* number of leaves for which weights are requested */
ZOLTAN_ID_PTR recv_gid = NULL;/* received message of global ids */
ZOLTAN_REFTREE *node;         /* a node in the refinement tree */
struct Zoltan_Reftree_hash_node **hashtab; /* hash table */
int hashsize;             /* dimension of hash table */
float *send_float;        /* sending message of floats */
float *recv_float;        /* received message of floats */
int num_gid_entries = lb->Num_GID; /* Number of array entries in a global ID. */

  /*
   * set the root and hash table
   */

  root = ((struct Zoltan_Reftree_data_struct *)lb->Data_Structure)->reftree_root;
  if (root == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Refinement tree not defined.");
    return(ZOLTAN_FATAL);
  }
  hashtab  = ((struct Zoltan_Reftree_data_struct *)lb->Data_Structure)->hash_table;
  hashsize = ((struct Zoltan_Reftree_data_struct *)lb->Data_Structure)->hash_table_size;

  /*
   * Determine the dimension of the weight array
   */

  if (lb->Obj_Weight_Dim == 0) {
    wdim = 1;
  } else {
    wdim = lb->Obj_Weight_Dim;
  }

  /*
   * In the first pass, sum the weights of the nodes that are assigned to
   * this processor, and count the leaves that are not.
   */

  count = 0;
  for (i=0; i<root->num_child; i++) {
    Zoltan_Reftree_Sum_My_Weights(lb,&(root->children[i]),&count,wdim);
  }
  root->assigned_to_me = -1;

  /*
   * Make a list of the leaves that are not assigned to this processor
   */

  if (count == 0)
    leaf_list = ZOLTAN_MALLOC_GID(lb);
  else
    leaf_list = ZOLTAN_MALLOC_GID_ARRAY(lb, count);
  if (leaf_list == NULL) {
    ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }

  count = 0;
  Zoltan_Reftree_List_Other_Leaves(lb, root,leaf_list,&count);

  /*
   * Get the unknown leaf weights from other processors.
   */

  nproc = lb->Num_Proc;
  myproc = lb->Proc;

  /*
   * Communication loop
   */

  for (comm_loop=0; comm_loop<nproc; comm_loop++) {

  /*
   * Broadcast the size of the request
   */

    reqsize = count;
    MPI_Bcast((void *)&reqsize, 1, MPI_INT, comm_loop, lb->Communicator);

  /*
   * If the request size is 0, skip this communication
   */

    if (reqsize != 0) {

  /*
   * Allocate space to receive the request
   */

      if (myproc != comm_loop) {
        recv_gid = ZOLTAN_MALLOC_GID_ARRAY(lb, reqsize);
        if (recv_gid == NULL) {
          ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
          ZOLTAN_FREE(&leaf_list);
          return(ZOLTAN_MEMERR);
        }
      }

  /*
   * broadcast the list of leaves processor comm_loop needs to know
   */

      /* KDDKDD Changed MPI_BYTE to ZOLTAN_ID_MPI_TYPE  */

      if (myproc == comm_loop)
        MPI_Bcast((void *)leaf_list,reqsize*num_gid_entries,ZOLTAN_ID_MPI_TYPE,
                   comm_loop, lb->Communicator);
      else
        MPI_Bcast((void *)recv_gid,reqsize*num_gid_entries,ZOLTAN_ID_MPI_TYPE,
                   comm_loop, lb->Communicator);

  /*
   *  Reply with any weights I have, and 0. where I don't have it
   */

      if (myproc != comm_loop) {

        send_float = (float *) ZOLTAN_MALLOC(sizeof(float)*wdim*reqsize);
        if (send_float == NULL) {
          ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
          ZOLTAN_FREE(&recv_gid);
          ZOLTAN_FREE(&leaf_list);
          return(ZOLTAN_MEMERR);
        }

        for (i=0; i<reqsize; i++) {
          node = Zoltan_Reftree_hash_lookup(lb, hashtab,
                                        &(recv_gid[i*num_gid_entries]),
                                        hashsize);
          if (node == NULL)
             for (j=0; j<wdim; j++) send_float[i*wdim+j] = 0.0;
          else
             for (j=0; j<wdim; j++) send_float[i*wdim+j] = node->my_sum_weight[j];
        }

        MPI_Send((void *)send_float,wdim*reqsize,MPI_FLOAT,comm_loop,
                 100+comm_loop,lb->Communicator);
        ZOLTAN_FREE(&send_float);
        ZOLTAN_FREE(&recv_gid);
      }

  /*
   * Receive the weights and add them into the tree
   */

      if (myproc == comm_loop) {
        recv_float = (float *)ZOLTAN_MALLOC(sizeof(float)*wdim*reqsize);
        if (recv_float == NULL) {
          ZOLTAN_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
          ZOLTAN_FREE(&leaf_list);
          return(ZOLTAN_MEMERR);
        }

        for (rproc=1; rproc<nproc; rproc++) {
          MPI_Recv((void *)recv_float,reqsize*wdim,MPI_FLOAT,MPI_ANY_SOURCE,
                   100+comm_loop, lb->Communicator,&mess_status);

          for (i=0; i<reqsize; i++) {
            node = Zoltan_Reftree_hash_lookup(lb, hashtab,
                                          &(leaf_list[i*num_gid_entries]),
                                          hashsize);
            for (j=0; j<wdim; j++) node->summed_weight[j] += recv_float[i];
          }
        }
        ZOLTAN_FREE(&recv_float);
      }
    }

  } /* end comm_loop */

  ZOLTAN_FREE(&leaf_list);

  /*
   * All the leaves now have summed_weight set.
   * Sum the weights throughout the tree.
   */

  Zoltan_Reftree_Sum_All_Weights(lb,root,wdim);

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Reftree_Sum_My_Weights(LB *lb, ZOLTAN_REFTREE *subroot, 
       int *count, int wdim)

{
/*
 * Function to recursively sum the weights of the nodes assigned
 * to this processor, and set assigned_to_me for interior nodes
 */
int i, j;          /* loop counter */
int none_assigned; /* flag for no children assigned to this proc */
int all_assigned;  /* flag for all children assigned to this proc */

  if (subroot->num_child == 0) {

  /*
   * If there are no children, then the sum is the weight of this node if
   * it is assigned to this processor, or 0 if it is not.
   */

    if (subroot->assigned_to_me) {
      for (i=0; i<wdim; i++) subroot->my_sum_weight[i] = subroot->weight[i];
      for (i=0; i<wdim; i++) subroot->summed_weight[i] = subroot->weight[i];
    }
    else {
      for (i=0; i<wdim; i++) subroot->my_sum_weight[i] = 0.0;
      for (i=0; i<wdim; i++) subroot->summed_weight[i] = 0.0;
      *count += 1;
    }

  }
  else {

  /*
   * If there are children, sum the weights of the children with the
   * node's weight and set assigned to me.
   */

    if (subroot->assigned_to_me) {
      for (i=0; i<wdim; i++) subroot->my_sum_weight[i] = subroot->weight[i];
    }
    else {
      for (i=0; i<wdim; i++) subroot->my_sum_weight[i] = 0.0;
    }
    none_assigned = 1;
    all_assigned = 1;

    for (j=0; j<subroot->num_child; j++) {
      Zoltan_Reftree_Sum_My_Weights(lb,&(subroot->children[j]),count,wdim);
      for (i=0; i<wdim; i++)
        subroot->my_sum_weight[i] += (subroot->children[j]).my_sum_weight[i];
      if ((subroot->children[j]).assigned_to_me == 1) none_assigned = 0;
      if ((subroot->children[j]).assigned_to_me == 0) all_assigned = 0;
      if ((subroot->children[j]).assigned_to_me == -1) {
        none_assigned = 0;
        all_assigned = 0;
      }
    }
    if (none_assigned)
      subroot->assigned_to_me = 0;
    else if (all_assigned)
      subroot->assigned_to_me = 1;
    else
      subroot->assigned_to_me = -1;
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Reftree_Sum_All_Weights(LB *lb, ZOLTAN_REFTREE *subroot, int wdim)

{
/*
 * Function to recursively sum the weights of all the nodes in the tree
 * assuming that summed_weight contains partial sums in the leaves
 */
int i, j;   /* loop counter */

  /*
   * If there are no children, do nothing
   */

  if (subroot->num_child != 0) {

  /*
   * If there are children, sum the weights of the children with the
   * weight of this node.
   */

    for (i=0; i<wdim; i++) subroot->summed_weight[i] = subroot->weight[i];

    for (j=0; j<subroot->num_child; j++) {
      Zoltan_Reftree_Sum_All_Weights(lb,&(subroot->children[j]),wdim);
      for (i=0; i<wdim; i++)
        subroot->summed_weight[i] += (subroot->children[j]).summed_weight[i];
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Reftree_List_Other_Leaves(LB *lb, ZOLTAN_REFTREE *subroot, 
       ZOLTAN_ID_PTR list, int *count)

{
/*
 * Function to make a list of the leaves not assigned to this processor
 */
int j;   /* loop counter */

  if (subroot->num_child == 0) {

  /*
   * If there are no children, then add it to the list if it is not
   * assigned to this processor
   */

    if (!subroot->assigned_to_me) {
      ZOLTAN_SET_GID(lb, &(list[(*count)*lb->Num_GID]),subroot->global_id);
      *count += 1;
    }

  }
  else {

  /*
   * If there are children, traverse the subtrees
   */

    for (j=0; j<subroot->num_child; j++) {
      Zoltan_Reftree_List_Other_Leaves(lb, &(subroot->children[j]),list,count);
    }
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int Zoltan_Reftree_Partition(LB *lb, int *num_export, 
       ZOLTAN_ID_PTR *export_global_ids, ZOLTAN_ID_PTR *export_local_ids, 
       int **export_procs)

{
/*
 * Function to partition the grid and fill the export arrays.
 */

char *yo = "Zoltan_Reftree_Partition";
char msg[256];
int num_exp;          /* count the number of export objects */
ZOLTAN_REFTREE *root;     /* root of the tree */
float partition_size; /* amount of weight for each partition */
float cutoff;         /* the value at which the current partition is full */
int part;             /* partition under construction */
float current_size;   /* amount of weight consumed so far */
float eps;            /* allowed deviation from average partition size */
int num_part;         /* number of partitions */

  root = ((struct Zoltan_Reftree_data_struct *)lb->Data_Structure)->reftree_root;

  /*
   * determine the size of the partitions and tolerance interval
   */

/* TEMP equal sizes does not support heterogeneous computer.  Can at least
        make the partition sizes proportional to the computational power
        of the processors.  Probably should set up an array of cutoff
        points up front, creating intervals whose sizes reflect the power
        of the processors.
        Don't know what to do about unequal communication */

/* TEMP using Num_Proc for number of partitions */

  num_part = lb->Num_Proc;
  partition_size = root->summed_weight[0]/num_part;
  eps = (lb->Imbalance_Tol - 1.0)*partition_size/2.0;

  /*
   * traverse the tree to define the partition and count the number of exports
   */

  num_exp = 0;
  part = 0;
  current_size = 0.0;
  cutoff = partition_size;
  Zoltan_Reftree_Part_Recursive(lb,root,&part,&current_size,&num_exp,&cutoff,
                            num_part,partition_size,eps);

  /*
   * if no exports, we're done
   */

  if (lb->LB_Return_Lists == ZOLTAN_LB_NO_LISTS) {
    return(ZOLTAN_OK);
  }
  else if (num_exp == 0) {
    *num_export = 0;
    return(ZOLTAN_OK);
  }

  /*
   * allocate space for the export lists
   */

  if (!Zoltan_Special_Malloc(lb,(void **)export_global_ids,num_exp,
                         ZOLTAN_SPECIAL_MALLOC_GID))
    return ZOLTAN_MEMERR;
  if (!Zoltan_Special_Malloc(lb,(void **)export_local_ids,num_exp,
                         ZOLTAN_SPECIAL_MALLOC_LID)) {
    Zoltan_Special_Free(lb,(void **)export_global_ids,ZOLTAN_SPECIAL_MALLOC_GID);
    return ZOLTAN_MEMERR;
  }
  if (!Zoltan_Special_Malloc(lb,(void **)export_procs,num_exp,
                         ZOLTAN_SPECIAL_MALLOC_INT)) {
    Zoltan_Special_Free(lb,(void **)export_global_ids,ZOLTAN_SPECIAL_MALLOC_GID);
    Zoltan_Special_Free(lb,(void **)export_local_ids,ZOLTAN_SPECIAL_MALLOC_LID);
    return ZOLTAN_MEMERR;
  }

  /*
   * traverse the tree to make the export lists
   */

  *num_export = 0;
  Zoltan_Reftree_Export_Lists(lb,root,num_export,export_global_ids,
                          export_local_ids,export_procs);

  if (num_exp != *num_export) {
    sprintf(msg, "num_exp = %d not equal to num_export = %d.",
            num_exp,*num_export);
    ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
    return(ZOLTAN_WARN);
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Reftree_Part_Recursive(LB *lb, ZOLTAN_REFTREE *subroot, int *part,
                              float *current_size, int *num_exp, float *cutoff,
                              int num_part, float partition_size, float eps)

{
/*
 * function to recursively define the partition and count the number of exports
 */
int i;         /* loop counter */
float newsize; /* size of partition if this subroot gets added to it */

/* TEMP don't know what to do with multicomponent weights.  Just using
        the first component */
  newsize = *current_size + subroot->summed_weight[0];

  if (newsize <= *cutoff + eps) {

  /*
   * This subtree fits in the current partition
   */

    subroot->partition = *part;
    *current_size = newsize;

  /*
   * If this is this processor's partition, there are no exports below this
   * node, so we don't have to traverse the subtree.
   * If there are no leaves of this subtree assigned to this processor, there
   * are no exports below this node.
   * Otherwise, traverse the subtree setting partition and counting exports
   */

    if (*part != lb->Proc && subroot->assigned_to_me) {
      if (subroot->num_child == 0)
        *num_exp += 1;
      else
        for (i=0; i<subroot->num_child; i++)
          Zoltan_Reftree_Mark_and_Count(&(subroot->children[i]), *part, num_exp);
    }

  /*
   * See if it is close enough to filling the partition
   */

    if (*current_size >= *cutoff - eps && *part < num_part-1) {
      *part += 1;
      *cutoff = (*part + 1)*partition_size;
    }
  }
  else {

  /*
   * This subtree is too big for the current partition
   */

  /*
   * If it has children, traverse them.
   */

    if (subroot->num_child != 0) {
      subroot->partition = -1;
      for (i=0; i<subroot->num_child; i++)
        Zoltan_Reftree_Part_Recursive(lb, &(subroot->children[i]),part,
               current_size, num_exp, cutoff, num_part, partition_size, eps);
    }
    else {

  /*
   * If there are no children, move on to the next partition
   */

      while (newsize > *cutoff+eps && *part < num_part-1) {
        *part += 1;
        *cutoff = (*part + 1)*partition_size;
      }
      subroot->partition = *part;
      *current_size = newsize;
      if (*part != lb->Proc && subroot->assigned_to_me) *num_exp += 1;
    }
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Reftree_Mark_and_Count(ZOLTAN_REFTREE *subroot, int part, 
       int *num_exp)

{
/*
 * Function to set the partition and count exports
 */
int i;

  subroot->partition = part;
  if (subroot->num_child == 0) {
    if (subroot->assigned_to_me) *num_exp += 1;
  } else {
    for (i=0; i<subroot->num_child; i++)
      Zoltan_Reftree_Mark_and_Count(&(subroot->children[i]), part, num_exp);
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void Zoltan_Reftree_Export_Lists(LB *lb, ZOLTAN_REFTREE *subroot, 
                            int *num_export, ZOLTAN_ID_PTR *export_global_ids,
                            ZOLTAN_ID_PTR *export_local_ids, int **export_procs)
{
/*
 * Function to build the export lists
 */
int i;

/*
 * if this subtree has no leaves assigned to this processor, or if the
 * partition of this subtree is this processor, then there can be no
 * exports below it
 */

  if (!subroot->assigned_to_me || subroot->partition == lb->Proc) return;

  if (subroot->num_child == 0) {

/*
 * if this is a leaf, put it on the export lists
 */

    ZOLTAN_SET_GID(lb, &((*export_global_ids)[(*num_export)*lb->Num_GID]),
               subroot->global_id);
    ZOLTAN_SET_LID(lb, &((*export_local_ids)[(*num_export)*lb->Num_LID]),
               subroot->local_id);
    (*export_procs)[*num_export] = subroot->partition;
    *num_export += 1;

  }
  else {

/*
 * if it is not a leaf, traverse the subtree
 */

    for (i=0; i<subroot->num_child; i++)
      Zoltan_Reftree_Export_Lists(lb, &(subroot->children[i]), num_export,
                              export_global_ids,export_local_ids,export_procs);
  }
}
