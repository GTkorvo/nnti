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


#include "zz_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  This file contains routines implementing the migration-help tools.
 *  These functions are all callable by the application.  
 */
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static int check_id_lengths(ZZ *);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_Migrate(
  ZZ *zz,                      /* Zoltan structure.                  */
  int num_import,              /* Number of non-local objects assigned to the 
                                  processor in the new decomposition.        */
  ZOLTAN_ID_PTR import_global_ids, /* Array of global IDs for non-local objects 
                                  assigned to this processor in the new
                                  decomposition; this field can be NULL if 
                                  the application doesn't provide import IDs.*/
  ZOLTAN_ID_PTR import_local_ids,  /* Array of local IDs for non-local objects
                                  assigned to the processor in the new
                                  decomposition; this field can be NULL if the 
                                  application does not provide import IDs.   */
  int *import_procs,           /* Array of processor IDs of processors owning
                                  the non-local objects that are assigned to
                                  this processor in the new decomposition; this
                                  field can be NULL if the application does
                                  not provide import IDs.                    */
  int *import_to_part,         /* Array of partition numbers to which imported
                                  objects should be assigned.                */
  int num_export,              /* Number of objs to be exported
                                  to other processors to establish the new
                                  decomposition.                             */
  ZOLTAN_ID_PTR export_global_ids, /* Array of global IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  ZOLTAN_ID_PTR export_local_ids,  /* Array of local IDs of
                                  objects to be exported to other processors
                                  to establish the new decomposition.        */
  int *export_procs,           /* Array of processor IDs
                                  to which objects will be exported 
                                  to establish the new decomposition.        */
  int *export_to_part          /* Array of partition numbers to which exported
                                  objects should be assigned.                */
)
{
/*
 *  Routine to help perform migration.  If migration pre-processing routine
 *  (ZOLTAN_PRE_MIGRATE_FN) is specified, this routine first calls that fn.
 *  It then calls a function to obtain the size of the migrating objects
 *  (ZOLTAN_OBJ_SIZE_FN).  The routine next calls an application-specified
 *  object packing routine (ZOLTAN_PACK_OBJ_FN) for each object
 *  to be exported.  It develops the needed communication map to move the
 *  objects to other processors.  It performs the communication according
 *  to the map, and then calls an application-specified object unpacking
 *  routine (ZOLTAN_UNPACK_OBJ_FN) for each object imported.
 */

char *yo = "Zoltan_Migrate";
char msg[256];
int num_gid_entries, num_lid_entries;  /* lengths of global & local ids */
int *sizes = NULL;       /* sizes (in bytes) of the object data for export. */
int id_size;             /* size (in bytes) of ZOLTAN_GID + padding for 
                            alignment                                       */
int tag_size;            /* size (in bytes) of ZOLTAN_GID + one int 
                            (for message size) */
char *export_buf = NULL; /* buffer for packing export data.                 */
char *import_buf = NULL; /* buffer for receiving imported data.             */
char *tmp;               /* temporary pointer into buffers.                 */
int i, j;                /* loop counter.                                   */
int actual_import;       /* number of objects to be imported.               */
int tmp_size;            /* size of a single object's data.                 */
int *idx = NULL;         /* index used for multi-fn packs and unpacks.      */
int idx_cnt = 0;         /* index counter for idx array.                    */
ZOLTAN_ID_PTR tmp_id = NULL; /* pointer to storage for a global ID in comm  
                                buf  */
ZOLTAN_ID_PTR lid;       /* temporary pointer to a local ID; used to pass
                            NULL to query functions when NUM_LID_ENTRIES=0. */
ZOLTAN_COMM_OBJ *comm_plan;     /* Object returned by communication routines*/
int msgtag, msgtag2;     /* Tags for communication routines                 */
int total_send_size;     /* Total size of outcoming message (in #items)     */
int total_recv_size;     /* Total size of incoming message (in #items)      */
int aligned_int;         /* size of an int padded for alignment             */
int dest;                /* temporary destination partition.                */
int ierr = 0;
int actual_num_exp = 0;
int actual_allocated = 0;
ZOLTAN_ID_PTR actual_exp_gids = NULL;    /* Arrays containing only objs to  */
ZOLTAN_ID_PTR actual_exp_lids = NULL;    /* actually be packed.  Objs that  */
int *actual_exp_procs = NULL;            /* are changing partition but not  */
int *actual_exp_to_part = NULL;          /* processor may not be included.  */

  ZOLTAN_TRACE_ENTER(zz, yo);

  /*
   *  Return if this processor is not in the Zoltan structure's
   *  communicator.
   */

  if (ZOLTAN_PROC_NOT_IN_COMMUNICATOR(zz)) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_OK);
  }

  /*
   *  Return if export lists are not provided (through faulty
   *  value of RETURN_LISTS parameter).
   */

  if (num_export == -1) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Export lists must be "
                                 "provided; change RETURN_LISTS parameter.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return(ZOLTAN_FATAL);
  }

  /*
   *  Check that all procs use the same id types.
   */

  ierr = check_id_lengths(zz);
  if (ierr != ZOLTAN_OK) {
    ZOLTAN_TRACE_EXIT(zz, yo);
    return ierr;
  }
  num_gid_entries = zz->Num_GID;
  num_lid_entries = zz->Num_LID;

  /*
   *  Check that all necessary query functions are available.
   */

  if (zz->Get_Obj_Size == NULL && zz->Get_Obj_Size_Multi == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register a "
           "ZOLTAN_OBJ_SIZE_FN or ZOLTAN_OBJ_SIZE_MULTI_FN function "
           "to use the migration-help tools.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_FATAL);
  }

  if (zz->Pack_Obj == NULL && zz->Pack_Obj_Multi == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register a "
           "ZOLTAN_PACK_OBJ_FN or ZOLTAN_PACK_OBJ_MULTI_FN function "
           "to use the migration-help tools.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_FATAL);
  }

  if (zz->Unpack_Obj == NULL && zz->Unpack_Obj_Multi == NULL) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Must register a "
         "ZOLTAN_UNPACK_OBJ_FN or ZOLTAN_UNPACK_MULTI_OBJ_FN function "
         "to use the migration-help tools.");
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ZOLTAN_FATAL);
  }

  /*
   *  Test whether to pack objects that have changed partition
   *  but not changed processor.  
   *  If packing them, the actual exports == exports passed to this function.
   *  If not packing them, build arrays with them stripped out.
   */

  if (!(zz->Migrate.Only_Proc_Changes)) {
    /* Pack all exports, even if they are not changing processor. */
    actual_num_exp = num_export;
    actual_exp_gids = export_global_ids;
    actual_exp_lids = export_local_ids;
    actual_exp_procs = export_procs;
    actual_exp_to_part = export_to_part;
  }
  else {  /* zz->Migrate.Only_Proc_Changes */
    /* Pack only exports that are actually changing processor. */
    actual_num_exp = 0;
    for (i = 0; i < num_export; i++) 
      if (export_procs[i] != zz->Proc)
        actual_num_exp++;

    if (actual_num_exp == num_export) {
      /*  Number of actual exports == number of exports in input arrays. */
      /*  No stripping needed. */
      actual_exp_gids = export_global_ids;
      actual_exp_lids = export_local_ids;
      actual_exp_procs = export_procs;
      actual_exp_to_part = export_to_part;
    }
    else if (actual_num_exp != num_export && actual_num_exp > 0) {
      /*  Number of actual exports < num_exports.  Build arrays  */
      /*  containing only actual exports. */
      actual_allocated = 1;
      actual_exp_gids = ZOLTAN_MALLOC_GID_ARRAY(zz, actual_num_exp);
      actual_exp_lids = ZOLTAN_MALLOC_LID_ARRAY(zz, actual_num_exp);
      actual_exp_procs = (int *) ZOLTAN_MALLOC(sizeof(int) * actual_num_exp);
      if (export_to_part != NULL)
        actual_exp_to_part = (int *) ZOLTAN_MALLOC(sizeof(int)*actual_num_exp);
      if (actual_exp_gids == NULL || actual_exp_lids == NULL ||
          actual_exp_procs == NULL || 
          (export_to_part != NULL && actual_exp_to_part == NULL)) {
        Zoltan_Multifree(__FILE__, __LINE__, 4, 
                         &actual_exp_gids, &actual_exp_lids, 
                         &actual_exp_procs, &actual_exp_to_part);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return (ZOLTAN_MEMERR);
      }
     
      for (j = 0, i = 0; i < num_export; i++) {
        if (export_procs[i] != zz->Proc) {
          ZOLTAN_SET_GID(zz, 
                        &(actual_exp_gids[j*num_gid_entries]),
                        &(export_global_ids[i*num_gid_entries]));
          if (num_lid_entries)
            ZOLTAN_SET_LID(zz, 
                          &(actual_exp_lids[j*num_lid_entries]),
                          &(export_local_ids[i*num_lid_entries]));
          actual_exp_procs[j] = export_procs[i];
          if (export_to_part) actual_exp_to_part[j] = export_to_part[i];
          j++;
        }
      }
    }
  }

  if (zz->Migrate.Pre_Migrate != NULL) {
    zz->Migrate.Pre_Migrate(zz->Migrate.Pre_Migrate_Data,
                            num_gid_entries, num_lid_entries,
                            num_import, import_global_ids,
                            import_local_ids, import_procs, import_to_part,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, export_to_part,
                            &ierr);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                      "Migrate.Pre_Migrate function.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done pre-migration processing");

  id_size = Zoltan_Align(num_gid_entries * sizeof(ZOLTAN_ID_TYPE));
  /* Note that alignment is not strictly necessary 
     when ZOLTAN_ID_TYPE is int or unsigned int. */
  aligned_int = Zoltan_Align(sizeof(int));
  tag_size = id_size + aligned_int;

  /*
   * For each object, allow space for its global ID and its data plus 
   * one int (for the object data size).
   * Zoltan will pack the global IDs; the application must pack the data
   * through the pack routine.  Zoltan needs the global IDs for unpacking,
   * as the order of the data received during communication is not 
   * necessarily the same order as import_global_ids[].
   * Zoltan also needs to communicate the sizes of the objects because
   * only the sender knows the size of each object.
   */
  if (actual_num_exp > 0) {
    sizes = (int *) ZOLTAN_MALLOC(actual_num_exp * sizeof(int));
    if (!sizes) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }

    if (zz->Get_Obj_Size_Multi != NULL) {
      zz->Get_Obj_Size_Multi(zz->Get_Obj_Size_Multi_Data, 
                             num_gid_entries, num_lid_entries, actual_num_exp,
                             actual_exp_gids, actual_exp_lids, sizes, &ierr);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                        "ZOLTAN_OBJ_SIZE_MULTI function.");
        ZOLTAN_TRACE_EXIT(zz, yo);
        return (ZOLTAN_FATAL);
      }
    }
    else {
      for (i = 0; i < actual_num_exp; i++){
        lid = (num_lid_entries ? &(actual_exp_lids[i*num_lid_entries]) : NULL);
        sizes[i] = zz->Get_Obj_Size(zz->Get_Obj_Size_Data, 
                       num_gid_entries, num_lid_entries,
                       &(actual_exp_gids[i*num_gid_entries]), 
                       lid, &ierr);
        if (ierr) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from "
                          "ZOLTAN_OBJ_SIZE function.");
          ZOLTAN_TRACE_EXIT(zz, yo);
          return (ZOLTAN_FATAL);
        }
      }
    }

    total_send_size = 0;

    for (i = 0; i < actual_num_exp; i++) {
      sizes[i] = Zoltan_Align(sizes[i]);
      total_send_size += sizes[i] + tag_size;
    }
    export_buf = (char *) ZOLTAN_MALLOC(total_send_size);
    if (!export_buf) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&sizes);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }

    if (zz->Pack_Obj_Multi != NULL) {
      /* Allocate an index array for ZOLTAN_PACK_OBJ_MULTI_FN. */
      idx = (int *) ZOLTAN_MALLOC(actual_num_exp * sizeof(int));
      if (!idx) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        ZOLTAN_FREE(&export_buf);
        ZOLTAN_FREE(&sizes);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return (ZOLTAN_FATAL);
      }
    }

    /*
     *  Pack the objects for export.
     */
  
    idx_cnt = 0;
    tmp = export_buf;
    for (i = 0; i < actual_num_exp; i++) {

      /* Pack the object's global ID */
      tmp_id = (ZOLTAN_ID_PTR) tmp;
      ZOLTAN_SET_GID(zz, tmp_id, &(actual_exp_gids[i*num_gid_entries]));
      tmp += id_size;
    
      /* Pack the object's size */
      *((int *)tmp) = sizes[i];
      tmp += aligned_int;

      /* If using ZOLTAN_PACK_OBJ_MULTI_FN, build the index array. */
      idx_cnt += tag_size;
      if (idx != NULL) {
        idx[i] = idx_cnt;
      }
      tmp += sizes[i];
      idx_cnt += sizes[i];
    }

    if (zz->Pack_Obj_Multi != NULL) {
      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] DEBUG in %s: Packing objects with multi-pack\n", 
               zz->Proc, yo);
      }
      zz->Pack_Obj_Multi(zz->Pack_Obj_Multi_Data,
                         num_gid_entries, num_lid_entries, actual_num_exp,
                         actual_exp_gids, actual_exp_lids, 
                         (actual_exp_to_part!=NULL ? actual_exp_to_part : actual_exp_procs),
                         sizes, idx, export_buf, &ierr);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                        "Pack_Obj function.");
        ZOLTAN_FREE(&sizes);
        ZOLTAN_FREE(&export_buf);
        ZOLTAN_FREE(&idx);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return (ZOLTAN_FATAL);
      }
    }
    else {
      tmp = export_buf + tag_size;
      for (i = 0; i < actual_num_exp; i++) {
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
          printf("[%1d] DEBUG in %s: Packing object with gid ", zz->Proc, yo);
          ZOLTAN_PRINT_GID(zz, &(actual_exp_gids[i*num_gid_entries]));
          printf("size = %d bytes\n", sizes[i]); 
        }

        /* Pack the object's data */
        lid = (num_lid_entries ? &(actual_exp_lids[i*num_lid_entries]) : NULL);
        dest = (actual_exp_to_part != NULL ? actual_exp_to_part[i] : actual_exp_procs[i]);
        zz->Pack_Obj(zz->Pack_Obj_Data, 
                           num_gid_entries, num_lid_entries,
                           &(actual_exp_gids[i*num_gid_entries]), lid, dest,
                           sizes[i], tmp, &ierr);
        if (ierr) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                          "Pack_Obj function.");
          ZOLTAN_FREE(&sizes);
          ZOLTAN_FREE(&export_buf);
          ZOLTAN_TRACE_EXIT(zz, yo);
          return (ZOLTAN_FATAL);
        }
        tmp += sizes[i] + tag_size;
      }
    }
    ZOLTAN_FREE(&idx);
    tmp_id = NULL;
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done packing objects");

  /*
   *  Compute communication map and actual_import, the number of objs this
   *  processor has to import to establish the new decomposition.
   */

  msgtag = 32767;
  ierr = Zoltan_Comm_Create(&comm_plan, actual_num_exp, actual_exp_procs, zz->Communicator,
                        msgtag, &actual_import);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    sprintf(msg, "Error %s returned from Zoltan_Comm_Create.", 
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_FREE(&sizes);
    ZOLTAN_FREE(&export_buf);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }

  /* Modify sizes[] to contain message sizes, not object sizes */
  for (i=0; i<actual_num_exp; i++)
    sizes[i] += tag_size;
  msgtag--;
  ierr = Zoltan_Comm_Resize(comm_plan, sizes, msgtag, &total_recv_size);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    sprintf(msg, "Error %s returned from Zoltan_Comm_Create.", 
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_FREE(&sizes);
    ZOLTAN_FREE(&export_buf);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }

  if (actual_import > 0) {
    import_buf = (char *) ZOLTAN_MALLOC(total_recv_size);
    if (!import_buf) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      ZOLTAN_FREE(&sizes);
      ZOLTAN_FREE(&export_buf);
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }
  }

  /*
   *  Send the export data using the communication plan.
   */

  msgtag2 = 32765;
  ierr = Zoltan_Comm_Do(comm_plan, msgtag2, export_buf, 1, import_buf);
  if (ierr != ZOLTAN_OK && ierr != ZOLTAN_WARN) {
    sprintf(msg, "Error %s returned from Zoltan_Comm_Do.", 
            (ierr == ZOLTAN_MEMERR ? "ZOLTAN_MEMERR" : "ZOLTAN_FATAL"));
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
    ZOLTAN_FREE(&sizes);
    ZOLTAN_FREE(&export_buf);
    ZOLTAN_FREE(&import_buf);
    Zoltan_Comm_Destroy(&comm_plan);
    ZOLTAN_TRACE_EXIT(zz, yo);
    return (ierr);
  }

  /*
   *  Free whatever memory we can.
   */

  Zoltan_Comm_Destroy(&comm_plan);
  ZOLTAN_FREE(&export_buf);
  ZOLTAN_FREE(&sizes);

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done communication");

  /* 
   *  Perform application-specified processing before unpacking the data.
   */
  if (zz->Migrate.Mid_Migrate != NULL) {
    zz->Migrate.Mid_Migrate(zz->Migrate.Mid_Migrate_Data,
                            num_gid_entries, num_lid_entries,
                            num_import, import_global_ids,
                            import_local_ids, import_procs, import_to_part,
                            num_export, export_global_ids,
                            export_local_ids, export_procs, export_to_part,
                            &ierr);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                      "Migrate.Mid_Migrate function.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }

    ZOLTAN_TRACE_DETAIL(zz, yo, "Done mid-migration processing");
  }

  /*
   *  Unpack the object data.
   */

  if (actual_import > 0) {

    if (zz->Unpack_Obj_Multi != NULL) {

      /* Allocate and fill input arrays for Unpack_Obj_Multi. */
      sizes = (int *) ZOLTAN_MALLOC(actual_import * sizeof(int));
      tmp_id = (ZOLTAN_ID_PTR) ZOLTAN_MALLOC_GID_ARRAY(zz, actual_import);
      idx = (int *) ZOLTAN_MALLOC(actual_import * sizeof(int));
      if (!sizes || !tmp_id || !idx) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
        ZOLTAN_FREE(&import_buf);
        ZOLTAN_FREE(&sizes);
        ZOLTAN_FREE(&tmp_id);
        ZOLTAN_TRACE_EXIT(zz, yo);
        return (ZOLTAN_FATAL);
      }

      tmp = import_buf;
      idx_cnt = 0;
      for (i = 0; i < actual_import; i++) {

        /* Unpack the object's global ID */
        ZOLTAN_SET_GID(zz, &(tmp_id[i*num_gid_entries]), (ZOLTAN_ID_PTR) tmp);
        tmp += id_size;

        /* Unpack the object's size */
        sizes[i] = *((int *)tmp);
        tmp += aligned_int;

        /* If using ZOLTAN_UNPACK_OBJ_MULTI_FN, build the index array. */
        idx_cnt += tag_size;
        if (idx != NULL) {
          idx[i] = idx_cnt;
        }

        tmp += sizes[i];
        idx_cnt += sizes[i];
      }

      if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
        printf("[%1d] DEBUG in %s: Unpacking objects with multi-fn\n",
               zz->Proc,yo);
      }
      zz->Unpack_Obj_Multi(zz->Unpack_Obj_Multi_Data, num_gid_entries,
                          actual_import, tmp_id, sizes, idx, import_buf, &ierr);
      ZOLTAN_FREE(&import_buf);
      ZOLTAN_FREE(&sizes);
      ZOLTAN_FREE(&tmp_id);
      ZOLTAN_FREE(&idx);
      if (ierr) {
        ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                        "ZOLTAN_UNPACK_OBJ_MULTI_FN.");
        ZOLTAN_TRACE_EXIT(zz, yo);
        return (ZOLTAN_FATAL);
      }
    }
    else {
      tmp = import_buf;
      for (i = 0; i < actual_import; i++) {
        tmp_size = *((int *)(tmp + id_size));
        if (zz->Debug_Level >= ZOLTAN_DEBUG_ALL){
          printf("[%1d] DEBUG in %s: Unpacking object with gid ", zz->Proc, yo);
          ZOLTAN_PRINT_GID(zz, (ZOLTAN_ID_PTR)tmp);
          printf("size = %d bytes\n", tmp_size);
        }

        /* Unpack the object's data */
       
        zz->Unpack_Obj(zz->Unpack_Obj_Data, num_gid_entries,
                       (ZOLTAN_ID_PTR) tmp, tmp_size,
                       tmp + tag_size, &ierr);
        if (ierr) {
          ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                          "ZOLTAN_UNPACK_OBJ_FN.");
          ZOLTAN_FREE(&import_buf);
          ZOLTAN_TRACE_EXIT(zz, yo);
          return (ZOLTAN_FATAL);
        }
        tmp += (tmp_size + tag_size);
      }
      ZOLTAN_FREE(&import_buf);
    }
  }

  ZOLTAN_TRACE_DETAIL(zz, yo, "Done unpacking objects");

  if (zz->Migrate.Post_Migrate != NULL) {
    zz->Migrate.Post_Migrate(zz->Migrate.Post_Migrate_Data,
                             num_gid_entries, num_lid_entries,
                             num_import, import_global_ids,
                             import_local_ids, import_procs, import_to_part,
                             num_export, export_global_ids,
                             export_local_ids, export_procs, export_to_part,
                             &ierr);
    if (ierr) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Error returned from user defined "
                      "Migrate.Post_Migrate function.");
      ZOLTAN_TRACE_EXIT(zz, yo);
      return (ZOLTAN_FATAL);
    }

    ZOLTAN_TRACE_DETAIL(zz, yo, "Done post-migration processing");
  }

  if (actual_allocated) {
    Zoltan_Multifree(__FILE__, __LINE__, 4, 
                     &actual_exp_gids, &actual_exp_lids, 
                     &actual_exp_procs, &actual_exp_to_part);
  }

  ZOLTAN_TRACE_EXIT(zz, yo);
  return (ZOLTAN_OK);
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static int check_id_lengths(
  ZZ *zz 
)
{
/* 
 * Routine to ensure that all processors have the same values of 
 * zz->Num_GID and zz->Num_LID.
 * All processors return the same error code.
 */
char *yo = "check_id_lengths";
char msg[256];
int loc_tmp[2];
int glob_min[2] = {0,0};
int glob_max[2] = {0,0};
int ierr = ZOLTAN_OK;

  loc_tmp[0] = zz->Num_GID;
  loc_tmp[1] = zz->Num_LID;

  /* 
   * Check both max and min values of IDs so that all processors can 
   * return the same error code. 
   */

  MPI_Allreduce(loc_tmp, glob_min, 2,
                MPI_INT, MPI_MIN, zz->Communicator);
  MPI_Allreduce(loc_tmp, glob_max, 2,
                MPI_INT, MPI_MAX, zz->Communicator);

  if ((glob_min[0] != glob_max[0]) ||
      (glob_min[1] != glob_max[1]))
    ierr = ZOLTAN_FATAL;

  if (zz->Num_GID != glob_max[0]){
    sprintf(msg, "Inconsistent global id sizes: Num_GID = %d "
      "but global max is %d\n", zz->Num_GID, glob_max[0]);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
  }

  if (zz->Num_LID != glob_max[1]){
    sprintf(msg, "Inconsistent local id sizes: Num_LID = %d "
      "but global max is %d\n", zz->Num_LID, glob_max[1]);
    ZOLTAN_PRINT_ERROR(zz->Proc, yo, msg);
  }

  return ierr;
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int Zoltan_Help_Migrate(
  ZZ *zz,
  int num_import,
  ZOLTAN_ID_PTR import_global_ids,
  ZOLTAN_ID_PTR import_local_ids,
  int *import_procs,
  int num_export,
  ZOLTAN_ID_PTR export_global_ids,
  ZOLTAN_ID_PTR export_local_ids,
  int *export_procs
)
{
/*
 *  Wrapper around Zoltan_Migrate with NULL pointers for partition arrays.
 *  Maintained for backward compatibility.
 *  Arguments are same as for Zoltan_Migrate.
 */

char *yo = "Zoltan_Help_Migrate";
int ierr;

  ZOLTAN_TRACE_ENTER(zz, yo);

  if (zz->LB.Num_Global_Parts != zz->Num_Proc) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
      "Number of partitions != Number of processors; use Zoltan_Migrate.");
    ierr = ZOLTAN_FATAL;
    goto End;
  }

  /*
   * Wrapper (for backward compatilibity) around Zoltan_Migrate.
   * Passes NULL for partition assignment arrays.
   */
  ierr = Zoltan_Migrate(zz, num_import, import_global_ids, import_local_ids,
                        import_procs, NULL,
                        num_export, export_global_ids, export_local_ids,
                        export_procs, NULL);

End:
  ZOLTAN_TRACE_EXIT(zz, yo);
  return ierr;
}

/****************************************************************************/
#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

