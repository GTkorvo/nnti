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

#include "lbi_const.h"
#include "lb_const.h"
#include "lb_util_const.h"
#include "rcb_const.h"
#include "rib_const.h"
#include "all_allo_const.h"
#include "comm_const.h"
#include "create_proc_list_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* PROTOTYPES */

static int initialize_dot(LB *, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, struct Dot_Struct *,
                          int, int, float *);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_RB_Build_Structure(
  LB *lb,                       /* load-balancing structure */
  ZOLTAN_ID_PTR *global_ids,        /* pointer to array of global IDs; allocated
                                   in this function.  */
  ZOLTAN_ID_PTR *local_ids,         /* pointer to array of local IDs; allocated
                                   in this function.  */
  struct Dot_Struct **dots,     /* pointer to array of Dots; allocated in this
                                   function. */
  int *num_obj,                 /* number of objects on this processor. */
  int *max_obj,                 /* number of Dots for which storage is 
                                   allocated on this processor. */
  int *num_geom,                /* # values per object used to describe
                                   the geometry.                       */
  int wgtflag,                  /* true is dot weights are to be used. */
  int use_ids                   /* true if global and local IDs are to be
                                   kept for RCB or RIB.  In all cases, the 
                                   IDs are allocated and used in Zoltan 
                                   callback functions.  Then, if use_ids is
                                   false, the IDs are deallocated before
                                   this function exits.                */
)
{
/*
 *  Function to build the geometry-based data structures for 
 *  RCB and RIB.
 */
char *yo = "LB_RB_Build_Structure";
char msg[256];
float *objs_wgt = NULL;               /* Array of object weights returned by 
                                         the application.                    */
int i, ierr = 0;

  /*
   * Compute the number of geometry fields per object.  This
   * value should be one, two or three, describing the x-, y-, and z-coords.
   */

  *num_geom = lb->Get_Num_Geom(lb->Get_Num_Geom_Data, &ierr);
  if (*num_geom > 3 || *num_geom < 1) {
    sprintf(msg, "Number of geometry fields %d is "
                  "invalid; valid range is 1-3\n",
                  *num_geom);
    LB_PRINT_ERROR(lb->Proc, yo, msg);
    return(ZOLTAN_FATAL);
  }
  if (ierr) {
    LB_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned from user function Get_Num_Geom.");
    return(ierr);
  }

  /*
   * Allocate space for objects.  Allow extra space
   * for objects that are imported to the processor.
   */

  *num_obj = lb->Get_Num_Obj(lb->Get_Num_Obj_Data, &ierr);
  if (ierr) {
    LB_PRINT_ERROR(lb->Proc, yo, 
                   "Error returned from user function Get_Num_Obj.");
    return(ierr);
  }

  *max_obj = (int)(1.5 * *num_obj) + 1;
  *global_ids = LB_MALLOC_GID_ARRAY(lb, (*max_obj));
  *local_ids  = LB_MALLOC_LID_ARRAY(lb, (*max_obj));
  *dots = (struct Dot_Struct *)LB_MALLOC((*max_obj)*sizeof(struct Dot_Struct));

  if (!(*global_ids) || (lb->Num_LID && !(*local_ids)) || !(*dots)) {
    LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
    return(ZOLTAN_MEMERR);
  }

  if (*num_obj > 0) {

    if (wgtflag) {

      /* 
       *  Allocate space for object weights.
       */

      objs_wgt    = (float *) LB_MALLOC((*num_obj)*sizeof(float));
      if (!objs_wgt) {
        LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
        return(ZOLTAN_MEMERR);
      }
      for (i = 0; i < *num_obj; i++) objs_wgt[i] = 0.;
    }

    /*
     *  Get list of objects' IDs and weights.
     */

    LB_Get_Obj_List(lb, *global_ids, *local_ids, wgtflag, objs_wgt, &ierr);
    if (ierr) {
      LB_PRINT_ERROR(lb->Proc, yo, 
                     "Error returned from user function LB_Get_Obj_List.");
      LB_FREE(&objs_wgt);
      return(ierr);
    }

    ierr = initialize_dot(lb, *global_ids, *local_ids, *dots,
                          *num_obj, wgtflag, objs_wgt);
    if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) {
      LB_PRINT_ERROR(lb->Proc, yo, 
                     "Error returned from user function initialize_dot.");
      LB_FREE(&objs_wgt);
      return(ierr);
    }

    LB_FREE(&objs_wgt);
  }

  if (!use_ids) {
    /* 
     * Do not need IDs in RCB or RIB.  Don't store, manipulate or
     * communicate any IDs.  
     */
    LB_FREE(global_ids);
    LB_FREE(local_ids);
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int initialize_dot(LB *lb, ZOLTAN_ID_PTR gid, ZOLTAN_ID_PTR lid,
                          struct Dot_Struct *dots, int num_obj,
                          int wgtflag, float *wgt)
{
/*
 *  Function that initializes the dot data structure for RCB and RIB. 
 *  It uses the global ID, coordinates and weight provided by the application.  
 */
int ierr = ZOLTAN_OK;
int i;
int num_gid_entries = lb->Num_GID;
int num_lid_entries = lb->Num_LID;
struct Dot_Struct *dot;
char *yo = "initialize_dot";

  for (i = 0; i < num_obj; i++) {
    dot = &(dots[i]);
    dot->Proc = lb->Proc;
    dot->X[0] = dot->X[1] = dot->X[2] = 0.0;
    lb->Get_Geom(lb->Get_Geom_Data, num_gid_entries, num_lid_entries,
                 &(gid[i*num_gid_entries]), &(lid[i*num_lid_entries]),
                 dot->X, &ierr);
    if (ierr == ZOLTAN_FATAL || ierr == ZOLTAN_MEMERR) {
      LB_PRINT_ERROR(lb->Proc, yo, 
                     "Error returned from user defined Get_Geom function.");
      return(ierr);
    }
    if (wgtflag)
       dot->Weight = wgt[i];
  }
  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_RB_Send_Outgoing(
  LB *lb,                           /* Load-balancing structure. */
  ZOLTAN_ID_PTR *gidpt,                 /* pointer to Global_IDs array. */
  ZOLTAN_ID_PTR *lidpt,                 /* pointer to Local_IDs array.  */
  struct Dot_Struct **dotpt,        /* pointer to Dots array. */
  int *dotmark,                     /* which side of median for each dot */
  int *dottop,                      /* dots >= this index are new */
  int *dotnum,                      /* number of dots */
  int *dotmax,                      /* max # of dots arrays can hold */
  int  set,                         /* which part processor is in = 0/1 */
  int *allocflag,                   /* have to re-allocate space */
  double overalloc,                 /* amount to overallocate by when realloc
                                       of dot array must be done.
                                       1.0 = no extra; 1.5 = 50% extra; etc. */
  int stats,                        /* Print timing & count summary? */
  int counters[],                   /* diagnostic counts
                                       0 = # of median iterations
                                       1 = # of dots sent
                                       2 = # of dots received
                                       3 = most dots this proc ever owns
                                       4 = most dot memory this proc ever allocs
                                       5 = # of times a previous cut is re-used
                                       6 = # of reallocs of dot array */
  int use_ids,                      /* true if global and local IDs are to be
                                       kept for RCB or RIB.  The IDs must be
                                       communicated if use_ids is true.  */
  MPI_Comm local_comm,
  int proclower,                    /* smallest processor for Tflops_Special */
  int numprocs                      /* number of processors for Tflops_Special*/
)
{
/* Routine to determine new processors for outgoing dots. */

  char *yo = "LB_RB_Send_Outgoing";
  int keep, outgoing;               /* message exchange counters */
  int *proc_list = NULL;            /* list of processors to send dots to */
  int i, ierr;

  /* outgoing = number of dots to ship to partner */
  /* dottop = number of dots that have never migrated */

  for (i = 0, keep = 0, outgoing = 0; i < *dotnum; i++) {
    if (dotmark[i] != set)
      outgoing++;
    else if (i < *dottop)
      keep++;
  }
  *dottop = keep;

  if (outgoing)
    if ((proc_list = (int *) LB_MALLOC(outgoing*sizeof(int))) == NULL) {
      LB_TRACE_EXIT(lb, yo);
      return ZOLTAN_MEMERR;
    }

  ierr = LB_Create_Proc_List(lb, set, *dotnum, outgoing, proc_list, local_comm,
                             proclower, numprocs);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    LB_PRINT_ERROR(lb->Proc, yo, "Error returned from LB_Create_Proc_List.");
    LB_FREE(&proc_list);
    LB_TRACE_EXIT(lb, yo);
    return (ierr);
  }

  ierr = LB_RB_Send_Dots(lb, gidpt, lidpt, dotpt, dotmark, proc_list, outgoing,
                         dotnum, dotmax, set, allocflag, overalloc, stats,
                         counters, use_ids, local_comm);

  if (outgoing) LB_FREE(&proc_list);

  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_RB_Send_Dots(
  LB *lb,                           /* Load-balancing structure. */
  ZOLTAN_ID_PTR *gidpt,                 /* pointer to Global_IDs array. */
  ZOLTAN_ID_PTR *lidpt,                 /* pointer to Local_IDs array.  */
  struct Dot_Struct **dotpt,        /* pointer to Dots array. */
  int *dotmark,                     /* which side of median for each dot */
  int *proc_list,                   /* list of processors to send dots to */
  int outgoing,                     /* message exchange counters */
  int *dotnum,                      /* number of dots */
  int *dotmax,                      /* max # of dots arrays can hold */
  int  set,                         /* which part processor is in = 0/1 */
  int *allocflag,                   /* have to re-allocate space */
  double overalloc,                 /* amount to overallocate by when realloc
                                       of dot array must be done.
                                       1.0 = no extra; 1.5 = 50% extra; etc. */
  int stats,                        /* Print timing & count summary? */
  int counters[],                   /* diagnostic counts
                                       0 = # of median iterations
                                       1 = # of dots sent
                                       2 = # of dots received
                                       3 = most dots this proc ever owns
                                       4 = most dot memory this proc ever allocs                                       5 = # of times a previous cut is re-used
                                       6 = # of reallocs of dot array */
  int use_ids,                      /* true if global and local IDs are to be
                                       kept for RCB or RIB (for Return_Lists 
                                       or high Debug_Levels).  The IDs must be
                                       communicated if use_ids is true.  */
  MPI_Comm local_comm
)
{
/* Routine to send outgoing dots to their new processors. */

  char *yo = "LB_RB_Send_Outgoing";
  int dotnew;                       /* # of new dots after send/recv */
  int keep, incoming;               /* message exchange counters */
  ZOLTAN_ID_PTR gidbuf = NULL;          /* communication buffer for global IDs. */
  ZOLTAN_ID_PTR lidbuf = NULL;          /* communication buffer for local IDs.  */
  struct Dot_Struct *dotbuf = NULL; /* communication buffer for dots. */
  struct Comm_Obj *cobj = NULL;     /* pointer for communication object */
  int message_tag = 32760;          /* message tag */
  int num_gid_entries = lb->Num_GID;
  int num_lid_entries = lb->Num_LID;
  int i, ierr;

  incoming = 0;
  ierr = LB_Comm_Create(&cobj, outgoing, proc_list, local_comm, message_tag,
                        &incoming);
  if (ierr != COMM_OK && ierr != COMM_WARN) {
    LB_PRINT_ERROR(lb->Proc, yo, "Error returned from LB_Comm_Create.");
    LB_FREE(&proc_list);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == COMM_MEMERR ? ZOLTAN_MEMERR : ZOLTAN_FATAL);
  }

  /* check if need to malloc more space */

  dotnew = *dotnum - outgoing + incoming;

  if (dotnew > *dotmax) {
    *allocflag = 1;
    *dotmax = (int) (overalloc * dotnew);
    if (*dotmax < dotnew) *dotmax = dotnew;
    if (use_ids) {
      *gidpt = LB_REALLOC_GID_ARRAY(lb, *gidpt, *dotmax);
      *lidpt = LB_REALLOC_LID_ARRAY(lb, *lidpt, *dotmax);
      if (!*gidpt || (num_lid_entries && !*lidpt)) {
        LB_TRACE_EXIT(lb, yo);
        return ZOLTAN_MEMERR;
      }
    }
    *dotpt = (struct Dot_Struct *) 
             LB_REALLOC(*dotpt,(unsigned) *dotmax * sizeof(struct Dot_Struct));
    if (!*dotpt) {
      LB_TRACE_EXIT(lb, yo);
      return ZOLTAN_MEMERR;
    }
    if (stats) counters[6]++;
  }

  if (stats) {
    counters[1] += outgoing;
    counters[2] += incoming;
    if (dotnew > counters[3])  counters[3] = dotnew;
    if (*dotmax > counters[4]) counters[4] = *dotmax;
  }
    
  /* malloc comm send buffer */

  if (outgoing > 0) {
    if (use_ids) {
      gidbuf = LB_MALLOC_GID_ARRAY(lb, outgoing);
      lidbuf = LB_MALLOC_LID_ARRAY(lb, outgoing);
      if (!gidbuf || (num_lid_entries && !lidbuf)) {
        LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
        LB_FREE(&gidbuf);
        LB_FREE(&lidbuf);
        LB_TRACE_EXIT(lb, yo);
        return ZOLTAN_MEMERR;
      }
    }
    dotbuf = (struct Dot_Struct *)
              LB_MALLOC(outgoing * sizeof(struct Dot_Struct));
    if (!dotbuf) {
      LB_PRINT_ERROR(lb->Proc, yo, "Insufficient memory.");
      LB_FREE(&gidbuf);
      LB_FREE(&lidbuf);
      LB_FREE(&dotbuf);
      LB_TRACE_EXIT(lb, yo);
      return ZOLTAN_MEMERR;
    }
  }

  /* fill buffer with dots that are marked for sending */
  /* pack down the unmarked ones */
    
  keep = outgoing = 0;
  for (i = 0; i < *dotnum; i++) {
    if (dotmark[i] != set) {
      if (use_ids) {
        LB_SET_GID(lb, &(gidbuf[outgoing*num_gid_entries]), 
                       &((*gidpt)[i*num_gid_entries]));
        LB_SET_LID(lb, &(lidbuf[outgoing*num_lid_entries]), 
                       &((*lidpt)[i*num_lid_entries]));
      }
      memcpy((char *) &dotbuf[outgoing], (char *) &((*dotpt)[i]), 
             sizeof(struct Dot_Struct));
      outgoing++;
    }
    else {
      if (use_ids) {
        LB_SET_GID(lb, &((*gidpt)[keep*num_gid_entries]), 
                       &((*gidpt)[i*num_gid_entries]));
        LB_SET_LID(lb, &((*lidpt)[keep*num_lid_entries]), 
                       &((*lidpt)[i*num_lid_entries]));
      }
      memcpy((char *) &((*dotpt)[keep]), (char *) &((*dotpt)[i]), 
             sizeof(struct Dot_Struct));
      keep++;
    }
  }

  if (use_ids) {
    /* Communicate Global IDs */
    ierr = LB_Comm_Do(cobj, message_tag, (char *) gidbuf, 
                      sizeof(ZOLTAN_ID_TYPE)*num_gid_entries,
                      (char *) &((*gidpt)[keep*num_gid_entries]));
    if (ierr != COMM_OK && ierr != COMM_WARN) {
      LB_PRINT_ERROR(lb->Proc, yo, "Error returned from LB_Comm_Do.");
      LB_FREE(&gidbuf);
      LB_FREE(&lidbuf);
      LB_FREE(&dotbuf);
      LB_TRACE_EXIT(lb, yo);
      return (ierr == COMM_MEMERR ? ZOLTAN_MEMERR : ZOLTAN_FATAL);
    }

    /* Communicate Local IDs, if any */
    if (num_lid_entries) {
      message_tag--;
      ierr = LB_Comm_Do(cobj, message_tag, (char *) lidbuf, 
                        sizeof(ZOLTAN_ID_TYPE)*num_lid_entries,
                        (char *) &((*lidpt)[keep*num_lid_entries]));
      if (ierr != COMM_OK && ierr != COMM_WARN) {
        LB_PRINT_ERROR(lb->Proc, yo, "Error returned from LB_Comm_Do.");
        LB_FREE(&gidbuf);
        LB_FREE(&lidbuf);
        LB_FREE(&dotbuf);
        LB_TRACE_EXIT(lb, yo);
        return (ierr == COMM_MEMERR ? ZOLTAN_MEMERR : ZOLTAN_FATAL);
      }
    }
  }

  /* Communicate Dots */
  message_tag--;
  ierr = LB_Comm_Do(cobj, message_tag, (char *) dotbuf, 
                    sizeof(struct Dot_Struct), (char *) &((*dotpt)[keep]));
  if (ierr != COMM_OK && ierr != COMM_WARN) {
    LB_PRINT_ERROR(lb->Proc, yo, "Error returned from LB_Comm_Do.");
    LB_FREE(&gidbuf);
    LB_FREE(&lidbuf);
    LB_FREE(&dotbuf);
    LB_TRACE_EXIT(lb, yo);
    return (ierr == COMM_MEMERR ? ZOLTAN_MEMERR : ZOLTAN_FATAL);
  }

  ierr = LB_Comm_Destroy(&cobj);

  if (use_ids) {
    LB_FREE(&gidbuf);
    LB_FREE(&lidbuf);
  }
  LB_FREE(&dotbuf);
    
  *dotnum = dotnew;

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_RB_Print_All(LB *lb, ZOLTAN_ID_PTR global_ids, struct Dot_Struct *dots,
                     int pdotnum, int pdottop, 
                     int num_import, ZOLTAN_ID_PTR import_global_ids, 
                     int *import_procs)
{
/*
 * Routine to print debugging information.  This routine runs serially
 * over all processors (due to LB_Print_Sync_*) and thus should be used
 * only for debugging.
 */
int kk;
int num_gid_entries = lb->Num_GID;

  LB_Print_Sync_Start(lb->Communicator, TRUE);
  printf("ZOLTAN Proc %d Num_Obj=%d Num_Keep=%d Num_Non_Local=%d\n",
         lb->Proc, pdotnum, pdottop, num_import);
  printf("  Assigned objects:\n");
  for (kk = 0; kk < pdotnum; kk++) {
     printf("    Obj:  ");
     LB_PRINT_GID(lb, &(global_ids[kk*num_gid_entries]));
     printf("  Orig: %4d\n", dots[kk].Proc);
  }
  printf("  Non_locals:\n");
  for (kk = 0; kk < num_import; kk++) {
     printf("    Obj:  ");
     LB_PRINT_GID(lb, &(import_global_ids[kk*num_gid_entries]));
     printf("     Orig: %4d\n", import_procs[kk]);
  }
  LB_Print_Sync_End(lb->Communicator, TRUE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_RB_Return_Arguments(
  LB *lb,                  /* Load-balancing structure */
  ZOLTAN_ID_PTR gidpt,         /* pointer to array of global IDs. */
  ZOLTAN_ID_PTR lidpt,         /* pointer to array of local IDs. */
  struct Dot_Struct *dotpt,/* pointer to array of Dots. */
  int num_import,          /* number of objects to be imported. */
  ZOLTAN_ID_PTR *import_global_ids,   /* global IDs of objects to be imported. */
  ZOLTAN_ID_PTR *import_local_ids,    /* local IDs of objects to be imported. */
  int **import_procs,             /* processors from which objects will be 
                                     imported. */
  int dottop               /* index of first dot to import on this processor. */
)
{
/*
 * Function to build the return arguments expected by Zoltan.
 * Allocates, fills and returns import_global_ids, import_local_ids, and
 * import_procs.
 */
char *yo = "LB_RB_Return_Arguments";
int i, ii;
int num_gid_entries = lb->Num_GID;
int num_lid_entries = lb->Num_LID;

  if (!LB_Special_Malloc(lb,(void **)import_global_ids,num_import,
                         LB_SPECIAL_MALLOC_GID)) {
    LB_TRACE_EXIT(lb, yo);
    return ZOLTAN_MEMERR;
  }
  if (!LB_Special_Malloc(lb,(void **)import_local_ids,num_import,
                         LB_SPECIAL_MALLOC_LID)) {
    LB_Special_Free(lb,(void **)import_global_ids,LB_SPECIAL_MALLOC_GID);
    LB_TRACE_EXIT(lb, yo);
    return ZOLTAN_MEMERR;
  }
  if (!LB_Special_Malloc(lb,(void **)import_procs,num_import,
                         LB_SPECIAL_MALLOC_INT)) {
    LB_Special_Free(lb,(void **)import_global_ids,LB_SPECIAL_MALLOC_GID);
    LB_Special_Free(lb,(void **)import_local_ids,LB_SPECIAL_MALLOC_LID);
    LB_TRACE_EXIT(lb, yo);
    return ZOLTAN_MEMERR;
  }

  for (i = 0; i < num_import; i++) {
    ii = i + dottop;
    LB_SET_GID(lb, &((*import_global_ids)[i*num_gid_entries]),
               &(gidpt[ii*num_gid_entries]));
    LB_SET_LID(lb, &((*import_local_ids)[i*num_lid_entries]),
               &(lidpt[ii*num_lid_entries]));
    (*import_procs)[i] = dotpt[ii].Proc;
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_Use_IDs(LB *lb)
{
/* Function that returns a flag indicating whether or not global and
 * local IDs should be stored, manipulated, and communicated within the
 * RCB and RIB algorithms.  
 * IDs are used for two purposes in RCB and RIB:
 * 1.  To build import lists.  When Return_Lists == NONE, import lists are not
 * built, so there is no need to carry along the IDs,
 * 2.  To be printed as debugging information.  Only for high debug levels are
 * the IDs printed; for lower levels, the IDs need not be carried along.
 */
  return (lb->Return_Lists || (lb->Debug_Level >= LB_DEBUG_ALL));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int LB_RB_check_geom_input(
  LB *lb,
  struct Dot_Struct *dotpt,
  int dotnum
)
{
/* Routine to check input to geometric methods for consistency. */
  char *yo = "LB_RB_check_geom_input";
  int i, j, k;
  char msg[256];
  int proc = lb->Proc;
  int ierr = ZOLTAN_OK;

  /* Error check the weights. */
  for (j = i = 0; i < dotnum; i++) if (dotpt[i].Weight == 0.0) j++;
  MPI_Allreduce(&j,&k,1,MPI_INT,MPI_SUM,lb->Communicator);
  if (k > 0 && proc == 0) {
     sprintf(msg, "%d dot weights are equal to 0.", k);
     LB_PRINT_WARN(proc, yo, msg);
     ierr = ZOLTAN_WARN;
  }

  for (j = i = 0; i < dotnum; i++) if (dotpt[i].Weight < 0.0) j++;
  MPI_Allreduce(&j,&k,1,MPI_INT,MPI_SUM,lb->Communicator);
  if (k > 0) {
    if (proc == 0) {
      sprintf(msg, "%d dot weights are < 0.",k);
      LB_PRINT_ERROR(proc, yo, msg);
    }
    ierr = ZOLTAN_FATAL;
  }
  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int LB_RB_check_geom_output(
  LB *lb, 
  struct Dot_Struct *dotpt,
  int dotnum,
  int dotorig,
  void *rcbbox_arg)
{
/* Routine to check output of geometric methods for consistency. */

  char *yo = "LB_RB_check_geom_output";
  char msg[256];
  int i,iflag,proc,nprocs,total1,total2;
  double weight,wtmax,wtmin,wtone,tolerance;
  struct rcb_box *rcbbox = (struct rcb_box *) rcbbox_arg;
  int ierr = ZOLTAN_OK;

  MPI_Comm_rank(lb->Communicator,&proc);
  MPI_Comm_size(lb->Communicator,&nprocs);

  /* check that total # of dots remained the same */

  MPI_Allreduce(&dotorig,&total1,1,MPI_INT,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&dotnum,&total2,1,MPI_INT,MPI_SUM,lb->Communicator);
  if (total1 != total2) {
    if (proc == 0) {
      sprintf(msg, "Points before partitioning = %d, "
                   "Points after partitioning = %d.",
                    total1,total2);
      LB_PRINT_WARN(proc, yo, msg);
      ierr = ZOLTAN_WARN;
    }
  }
  
  /* check that result is load-balanced within log2(P)*max-wt */

  weight = wtone = 0.0;
  for (i = 0; i < dotnum; i++) {
    weight += dotpt[i].Weight;
    if (dotpt[i].Weight > wtone) wtone = dotpt[i].Weight;
  }

  MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  MPI_Allreduce(&wtone,&tolerance,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);

  /* i = smallest power-of-2 >= nprocs */
  /* tolerance = largest-single-weight*log2(nprocs) */

  for (i = 0; (nprocs >> i) != 0; i++);
  tolerance = tolerance * i * (1.0 + TINY);

  if (wtmax - wtmin > tolerance) {
    if (proc == 0) {
      sprintf(msg, "Load-imbalance > tolerance of %g.",
              tolerance);
      LB_PRINT_WARN(proc, yo, msg);
      ierr = ZOLTAN_WARN;
    }
    MPI_Barrier(lb->Communicator);
    if (weight == wtmin) {
      sprintf(msg, "  Proc %d has weight = %g.",proc,weight);
      LB_PRINT_WARN(proc, yo, msg);
      ierr = ZOLTAN_WARN;
    }
    if (weight == wtmax) {
      sprintf(msg, "  Proc %d has weight = %g.",proc,weight);
      LB_PRINT_WARN(proc, yo, msg);
      ierr = ZOLTAN_WARN;
    }
  }
  
  MPI_Barrier(lb->Communicator);
  
  if (lb->Method == RCB) {

    /* check that final set of points is inside RCB box of each proc */
  
    iflag = 0;
    for (i = 0; i < dotnum; i++) {
      if (dotpt[i].X[0] < rcbbox->lo[0] || dotpt[i].X[0] > rcbbox->hi[0] ||
          dotpt[i].X[1] < rcbbox->lo[1] || dotpt[i].X[1] > rcbbox->hi[1] ||
  	  dotpt[i].X[2] < rcbbox->lo[2] || dotpt[i].X[2] > rcbbox->hi[2])
      iflag++;
    }
    if (iflag > 0) {
      sprintf(msg, "%d points are out-of-box on proc %d.", iflag, proc);
      LB_PRINT_ERROR(proc, yo, msg);
      ierr = ZOLTAN_FATAL;
    }
  
    MPI_Barrier(lb->Communicator);
  }
  return(ierr);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void LB_RB_stats(LB *lb, double timetotal, struct Dot_Struct *dotpt,
                 int dotnum, double *timers, int *counters, int stats,
                 int *reuse_count, void *rcbbox_arg, int reuse)

{
  int i,proc,nprocs,sum,min,max,print_proc;
  double ave,rsum,rmin,rmax;
  double weight,wttot,wtmin,wtmax;
  struct rcb_box *rcbbox = (struct rcb_box *) rcbbox_arg;

  MPI_Comm_rank(lb->Communicator,&proc);
  MPI_Comm_size(lb->Communicator,&nprocs);
  print_proc = lb->Debug_Proc;
  
  if (proc == print_proc) 
    printf("Partitioning total time: %g (secs)\n", timetotal);

  if (stats) {
    if (proc == print_proc) printf("Partitioning Statistics:\n");

    MPI_Barrier(lb->Communicator);

    /* distribution info */
  
    for (i = 0, weight = 0.0; i < dotnum; i++) weight += dotpt[i].Weight;
    MPI_Allreduce(&weight,&wttot,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
    MPI_Allreduce(&weight,&wtmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
    MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);

    if (proc == print_proc) {
      printf(" Total weight of dots = %g\n",wttot);
      printf(" Weight on each proc: ave = %g, max = %g, min = %g\n",
	     wttot/nprocs,wtmax,wtmin);
    }

    MPI_Barrier(lb->Communicator);
    if (stats > 1)  
      printf("    Proc %d has weight = %g\n",proc,weight);

    for (i = 0, weight = 0.0; i < dotnum; i++) 
      if (dotpt[i].Weight > weight) weight = dotpt[i].Weight;
    MPI_Allreduce(&weight,&wtmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  
    if (proc == print_proc) 
      printf(" Maximum weight of single dot = %g\n",wtmax);

    MPI_Barrier(lb->Communicator);
    if (stats > 1)  
      printf("    Proc %d max weight = %g\n",proc,weight);
  
    /* counter info */
  
    MPI_Allreduce(&counters[0],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
    MPI_Allreduce(&counters[0],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
    MPI_Allreduce(&counters[0],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == print_proc) 
      printf(" Median iter: ave = %g, min = %d, max = %d\n",ave,min,max);
    MPI_Barrier(lb->Communicator);
    if (stats > 1)  
      printf("    Proc %d median count = %d\n",proc,counters[0]);
  
    MPI_Allreduce(&counters[1],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
    MPI_Allreduce(&counters[1],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
    MPI_Allreduce(&counters[1],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == print_proc) 
      printf(" Send count: ave = %g, min = %d, max = %d\n",ave,min,max);
    MPI_Barrier(lb->Communicator);
    if (stats > 1) 
      printf("    Proc %d send count = %d\n",proc,counters[1]);
    
    MPI_Allreduce(&counters[2],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
    MPI_Allreduce(&counters[2],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
    MPI_Allreduce(&counters[2],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == print_proc) 
      printf(" Recv count: ave = %g, min = %d, max = %d\n",ave,min,max);
    MPI_Barrier(lb->Communicator);
    if (stats > 1) 
      printf("    Proc %d recv count = %d\n",proc,counters[2]);
  
    if (reuse) {
      MPI_Allreduce(&reuse_count[1],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
      MPI_Allreduce(&reuse_count[1],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
      MPI_Allreduce(&reuse_count[1],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
      ave = ((double) sum)/nprocs;
      if (proc == print_proc) 
        printf(" Presend count: ave = %g, min = %d, max = %d\n",ave,min,max);
      MPI_Barrier(lb->Communicator);
      if (stats > 1) 
        printf("    Proc %d presend count = %d\n",proc,reuse_count[1]);
    
      MPI_Allreduce(&reuse_count[2],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
      MPI_Allreduce(&reuse_count[2],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
      MPI_Allreduce(&reuse_count[2],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
      ave = ((double) sum)/nprocs;
      if (proc == print_proc) 
        printf(" Prerecv count: ave = %g, min = %d, max = %d\n",ave,min,max);
      MPI_Barrier(lb->Communicator);
      if (stats > 1) 
        printf("    Proc %d prerecv count = %d\n",proc,reuse_count[2]);
    }
  
    MPI_Allreduce(&counters[3],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
    MPI_Allreduce(&counters[3],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
    MPI_Allreduce(&counters[3],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == print_proc) 
      printf(" Max dots: ave = %g, min = %d, max = %d\n",ave,min,max);
    MPI_Barrier(lb->Communicator);
    if (stats > 1) 
      printf("    Proc %d max dots = %d\n",proc,counters[3]);
  
    MPI_Allreduce(&counters[4],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
    MPI_Allreduce(&counters[4],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
    MPI_Allreduce(&counters[4],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == print_proc) 
      printf(" Max memory: ave = %g, min = %d, max = %d\n",ave,min,max);
    MPI_Barrier(lb->Communicator);
    if (stats > 1) 
      printf("    Proc %d max memory = %d\n",proc,counters[4]);
    
    if (reuse) {
      MPI_Allreduce(&counters[5],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
      MPI_Allreduce(&counters[5],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
      MPI_Allreduce(&counters[5],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
      ave = ((double) sum)/nprocs;
      if (proc == print_proc) 
        printf(" # of Reuse: ave = %g, min = %d, max = %d\n",ave,min,max);
      MPI_Barrier(lb->Communicator);
      if (stats > 1) 
        printf("    Proc %d # of Reuse = %d\n",proc,counters[5]);
    }
  
    MPI_Allreduce(&counters[6],&sum,1,MPI_INT,MPI_SUM,lb->Communicator);
    MPI_Allreduce(&counters[6],&min,1,MPI_INT,MPI_MIN,lb->Communicator);
    MPI_Allreduce(&counters[6],&max,1,MPI_INT,MPI_MAX,lb->Communicator);
    ave = ((double) sum)/nprocs;
    if (proc == print_proc) 
      printf(" # of OverAlloc: ave = %g, min = %d, max = %d\n",ave,min,max);
    MPI_Barrier(lb->Communicator);
    if (stats > 1) 
      printf("    Proc %d # of OverAlloc = %d\n",proc,counters[6]);
  }

  /* timer info */
  
  MPI_Allreduce(&timers[0],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&timers[0],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&timers[0],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  ave = rsum/nprocs;
  if (proc == print_proc) 
    printf(" Start-up time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(lb->Communicator);
  if (stats > 1) 
    printf("    Proc %d start-up time = %g\n",proc,timers[0]);
  
  MPI_Allreduce(&timers[1],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&timers[1],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&timers[1],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  ave = rsum/nprocs;
  if (proc == print_proc) 
    printf(" Pre-median time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(lb->Communicator);
  if (stats > 1) 
    printf("    Proc %d pre-median time = %g\n",proc,timers[1]);
  
  MPI_Allreduce(&timers[2],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&timers[2],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&timers[2],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  ave = rsum/nprocs;
  if (proc == print_proc) 
    printf(" Median time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(lb->Communicator);
  if (stats > 1) 
    printf("    Proc %d median time = %g\n",proc,timers[2]);
  
  MPI_Allreduce(&timers[3],&rsum,1,MPI_DOUBLE,MPI_SUM,lb->Communicator);
  MPI_Allreduce(&timers[3],&rmin,1,MPI_DOUBLE,MPI_MIN,lb->Communicator);
  MPI_Allreduce(&timers[3],&rmax,1,MPI_DOUBLE,MPI_MAX,lb->Communicator);
  ave = rsum/nprocs;
  if (proc == print_proc) 
    printf(" Comm time %%: ave = %g, min = %g, max = %g\n",
	   ave/timetotal*100.0,rmin/timetotal*100.0,rmax/timetotal*100.0);
  MPI_Barrier(lb->Communicator);
  if (stats > 1) 
    printf("    Proc %d comm time = %g\n",proc,timers[3]);
  
  
  if (lb->Method == RCB && stats > 1)  {
    /* RCB boxes for each proc */
    if (proc == print_proc) printf(" RCB sub-domain boxes:\n");
    for (i = 0; i < 3; i++) {
      MPI_Barrier(lb->Communicator);
      if (proc == print_proc) printf("    Dimension %d\n",i+1);
      MPI_Barrier(lb->Communicator);
      printf("      Proc = %d: Box = %g %g\n",
	     proc,rcbbox->lo[i],rcbbox->hi[i]);
    }
  }

  MPI_Barrier(lb->Communicator);
}
