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
#include "lb_init_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void Zoltan_Migrate_Init(struct Zoltan_Migrate_Struct *mig)
{
  mig->Auto_Migrate = ZOLTAN_AUTO_MIGRATE_DEF;
  mig->Only_Proc_Changes = ZOLTAN_MIGRATE_ONLY_PROC_CHANGES_DEF;
  mig->Pre_Migrate = NULL;
  mig->Mid_Migrate = NULL;
  mig->Post_Migrate = NULL;
  mig->Pre_Migrate_Fort = NULL;
  mig->Mid_Migrate_Fort = NULL;
  mig->Post_Migrate_Fort = NULL;
}

void Zoltan_LB_Init(struct Zoltan_LB_Struct *lb, int num_proc)
{
  lb->Num_Global_Parts = num_proc;
  lb->Num_Local_Parts = 1;
  lb->PartDist = NULL;
  lb->Method = RCB;
  lb->LB_Fn = Zoltan_RCB;
  lb->Return_Lists = ZOLTAN_LB_RETURN_LISTS_DEF;
  lb->Imbalance_Tol = ZOLTAN_LB_IMBALANCE_TOL_DEF;
  lb->Data_Structure = NULL;
  lb->Free_Structure = Zoltan_RCB_Free_Structure;
  lb->Point_Assign = Zoltan_RB_Point_Assign;
  lb->Box_Assign = Zoltan_RB_Box_Assign;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
