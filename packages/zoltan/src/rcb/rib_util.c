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
#include "rib.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int Zoltan_RIB_Build_Structure(ZZ *zz, int *num_obj, int *max_obj, int wgtflag,
                           int use_ids)
{
/* Function to build the geometry-based data structures for RIB method. */
char           *yo = "Zoltan_RIB_Build_Structure";
RIB_STRUCT     *rib;                  /* Data structure for RIB.             */
struct rib_tree *treeptr;
int            i, ierr = 0;

  /* Allocate an RIB data structure for this Zoltan structure.
     If the previous data structure is still there, free the Dots and IDs first;
     the other fields can be reused. */

  if (zz->LB.Data_Structure == NULL) {
    rib = (RIB_STRUCT *) ZOLTAN_MALLOC(sizeof(RIB_STRUCT));
    if (rib == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      return(ZOLTAN_MEMERR);
    }
    zz->LB.Data_Structure = (void *) rib;
    rib->Tree_Ptr = NULL;
    rib->Global_IDs = NULL;
    rib->Local_IDs = NULL;
    rib->Dots = NULL;

    rib->Tree_Ptr = (struct rib_tree *)
              ZOLTAN_MALLOC(zz->LB.Num_Global_Parts* sizeof(struct rib_tree));
    if (rib->Tree_Ptr == NULL) {
      ZOLTAN_PRINT_ERROR(zz->Proc, yo, "Insufficient memory.");
      Zoltan_RIB_Free_Structure(zz);
      return(ZOLTAN_MEMERR);
    }
    /* initialize Tree_Ptr */
    for (i = 0; i < zz->LB.Num_Global_Parts; i++) {
      treeptr = &(rib->Tree_Ptr[i]);
      treeptr->cm[0] = treeptr->cm[1] = treeptr->cm[2] = 0.0;
      treeptr->ev[0] = treeptr->ev[1] = treeptr->ev[2] = 0.0;
      treeptr->cut = 0.0;
      treeptr->parent = treeptr->left_leaf = treeptr->right_leaf = 0;
    }
  }
  else {
    rib = (RIB_STRUCT *) zz->LB.Data_Structure;
    ZOLTAN_FREE(&(rib->Global_IDs));
    ZOLTAN_FREE(&(rib->Local_IDs));
    ZOLTAN_FREE(&(rib->Dots));
  }

  ierr = Zoltan_RB_Build_Structure(zz, &(rib->Global_IDs), &(rib->Local_IDs),
                               &(rib->Dots), num_obj, max_obj, &(rib->Num_Geom),
                               wgtflag, use_ids);
  if (ierr) {
    ZOLTAN_PRINT_ERROR(zz->Proc, yo,
                       "Error returned from Zoltan_RB_Build_Structure.");
    Zoltan_RIB_Free_Structure(zz);
    return(ierr);
  }

  return(ZOLTAN_OK);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_RIB_Free_Structure(ZZ *zz)
{
/* Deallocate the persistent RIB data structures in zz->Structure.  */
RIB_STRUCT    *rib;                   /* Data structure for RIB. */

  rib = (RIB_STRUCT *) (zz->LB.Data_Structure);

  if (rib != NULL) {
    ZOLTAN_FREE(&(rib->Tree_Ptr));
    ZOLTAN_FREE(&(rib->Global_IDs));
    ZOLTAN_FREE(&(rib->Local_IDs));
    ZOLTAN_FREE(&(rib->Dots));
    ZOLTAN_FREE(&(zz->LB.Data_Structure));
  }
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#define COPY_BUFFER(buf, type, allocLen, copyLen) \
  if (from->buf) { \
    to->buf = (type *)ZOLTAN_MALLOC((allocLen) * sizeof(type)); \
    if (!to->buf) { \
      Zoltan_RIB_Free_Structure(toZZ); \
      ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory."); \
      return(ZOLTAN_MEMERR); \
    } \
    if (copyLen > 0){ \
      memcpy(to->buf, from->buf, (copyLen) * sizeof(type)); \
    } else { \
      memset(to->buf, 0, allocLen); \
    } \
  } \
  else { \
    to->buf = NULL; \
  }

int Zoltan_RIB_Copy_Structure(ZZ *toZZ, ZZ *fromZZ)
{
  char *yo = "Zoltan_RIB_Copy_Structure";
  int num_obj, max_obj, rc;
  RIB_STRUCT *to, *from;
  ZOLTAN_ID_PTR gids, lids;

  from = (RIB_STRUCT *)fromZZ->LB.Data_Structure;
  Zoltan_RIB_Free_Structure(toZZ);

  if (!from){
    return ZOLTAN_OK;
  }

  rc = Zoltan_Copy_Obj_List(fromZZ, from->Global_IDs, &gids,
    from->Local_IDs, &lids,  0, NULL, NULL, NULL, NULL, &num_obj);

  if ((rc != ZOLTAN_OK) && (rc != ZOLTAN_WARN)){
    return rc;
  }

  to = (RIB_STRUCT *)ZOLTAN_MALLOC(sizeof(RIB_STRUCT));
  if (to == NULL) {
    ZOLTAN_PRINT_ERROR(fromZZ->Proc, yo, "Insufficient memory.");
    ZOLTAN_FREE(&gids);
    ZOLTAN_FREE(&lids);
    return(ZOLTAN_MEMERR);
  }

  memset(to, 0, sizeof(RIB_STRUCT));
  toZZ->LB.Data_Structure = (void *)to;

  max_obj = (int)(1.5 * num_obj) + 1;
  
  to->Global_IDs = ZOLTAN_REALLOC_GID_ARRAY(fromZZ, gids, max_obj);
  to->Local_IDs = ZOLTAN_REALLOC_GID_ARRAY(fromZZ, lids, max_obj);

  COPY_BUFFER(Dots, struct Dot_Struct, max_obj, num_obj);

  if (fromZZ->LB.Num_Global_Parts){
    COPY_BUFFER(Tree_Ptr, struct rib_tree,
              fromZZ->LB.Num_Global_Parts, fromZZ->LB.Num_Global_Parts);
  }

  to->Num_Geom = from->Num_Geom;

  return ZOLTAN_OK;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
