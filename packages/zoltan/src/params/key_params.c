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

#include <stdio.h>
#include "lb_const.h"
#include "params_const.h"

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static PARAM_VARS Key_params[] = {
	{ "IMBALANCE_TOL", NULL, "DOUBLE" },
	{ "AUTO_MIGRATE", NULL, "INT" },
	{ "OBJ_WEIGHT_DIM", NULL, "INT" },
	{ "COMM_WEIGHT_DIM", NULL, "INT" },
	{ "DEBUG_LEVEL", NULL, "INT" },
	{ "DEBUG_PROCESSOR", NULL, "INT" },
	{ "DETERMINISTIC", NULL, "INT" },
	{ "TIMER", NULL, "STRING" },
	{ "NUM_GID_ENTRIES", NULL, "INT" },
	{ "NUM_LID_ENTRIES", NULL, "INT" },
	{ "RETURN_LISTS", NULL, "STRING" },
        { "TFLOPS_SPECIAL", NULL, "INT" },
	{ NULL, NULL, NULL } };
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* 
 * Handle parameter changes for variables stored in LB structure.
 */

int Zoltan_Set_Key_Param(
LB *lb,                         /* load balance structure */
char *name,			/* name of variable */
char *val)			/* value of variable */
{
    char *yo = "Zoltan_Set_Key_Param";
    char msg[256];
    int status;			/* return code */
    PARAM_UTYPE result;		/* value returned from Check_Param */
    int index;			/* index returned from Check_Param */
    int tmp;

    status = Zoltan_Check_Param(name, val, Key_params, &result, &index);

    if (status == 0) {

      switch (index) {

      case 0:  		/* Imbalance_Tol */
        if (result.def) 
            result.dval = ZOLTAN_LB_IMBALANCE_TOL_DEF;
	if (result.dval < 1.0) {
	    sprintf(msg, "Invalid Imbalance_Tol value (%g) "
		"being set to %g.", result.dval, ZOLTAN_LB_IMBALANCE_TOL_DEF);
            ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
	    result.dval = ZOLTAN_LB_IMBALANCE_TOL_DEF;
	}
	lb->Imbalance_Tol = result.dval;
	status = 3;		/* Don't add to Params field of LB */
        break;

      case 1:		/* Help_Migrate */
        if (result.def)
            result.ival = ZOLTAN_AUTO_MIGRATE_DEF;
	lb->Migrate.Auto_Migrate = result.ival;
	status = 3;		/* Don't add to Params field of LB */
        break;

      case 2:		/* Object weight dim.  */
        if (result.def)
            result.ival = ZOLTAN_OBJ_WEIGHT_DEF;
	if (result.ival < 0) {
	    sprintf(msg, "Invalid Obj_Weight_Dim value (%d) "
		"being set to %d.", result.ival, ZOLTAN_OBJ_WEIGHT_DEF);
            ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
	    result.ival = ZOLTAN_OBJ_WEIGHT_DEF;
	}
	lb->Obj_Weight_Dim = result.ival;
	status = 3;		/* Don't add to Params field of LB */
        break;

      case 3:		/* Communication weight dim.  */
        if (result.def)
            result.ival = ZOLTAN_COMM_WEIGHT_DEF;
	if (result.ival < 0) {
	    sprintf(msg, "Invalid Comm_Weight_Dim value (%d) "
		"being set to %d.", result.ival, ZOLTAN_COMM_WEIGHT_DEF);
            ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
	    result.ival = ZOLTAN_COMM_WEIGHT_DEF;
	}
	lb->Comm_Weight_Dim = result.ival;
	status = 3;		/* Don't add to Params field of LB */
        break;

      case 4: 		/* Debug level  */
        if (result.def)
            result.ival = ZOLTAN_DEBUG_LEVEL_DEF;
	if (result.ival < 0) {
	    sprintf(msg, "Invalid Debug_Level value (%d) "
		"being set to %d.", result.ival, ZOLTAN_DEBUG_LEVEL_DEF);
            ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
	    result.ival = ZOLTAN_DEBUG_LEVEL_DEF;
	}
	lb->Debug_Level = result.ival;
	status = 3;		/* Don't add to Params field of LB */
        break;

      case 5: 		/* Debug processor  */
        if (result.def)
            result.ival = ZOLTAN_DEBUG_PROC_DEF;
	if (result.ival < 0 || result.ival > lb->Num_Proc) {
	    sprintf(msg, "Invalid Debug_Processor value (%d) "
		"being set to %d.", result.ival, ZOLTAN_DEBUG_PROC_DEF);
            ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
	    result.ival = ZOLTAN_DEBUG_PROC_DEF;
	}
	lb->Debug_Proc = result.ival;
	status = 3;		/* Don't add to Params field of LB */
        break;
       
      case 6: 		/* Deterministic flag */
        if (result.def)
            result.ival = ZOLTAN_DETERMINISTIC_DEF;
	if (result.ival < 0) {
	    sprintf(msg, "Invalid Deterministic value (%d) "
		"being set to %d.", result.ival, ZOLTAN_DETERMINISTIC_DEF);
            ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
	    result.ival = ZOLTAN_DETERMINISTIC_DEF;
	}
	lb->Deterministic = result.ival;
	status = 3;		/* Don't add to Params field of LB */
        break;

      case 7: 		/* Timer */
	status = Zoltan_Set_Timer_Param(name, val, &tmp);
        lb->Timer = tmp;

	if (status==0) status = 3;	/* Don't add to Params field of LB */
        break;

      case 8:           /* Num_GID_Entries */
        if (result.def)
            result.ival = ZOLTAN_NUM_ID_ENTRIES_DEF;
        if (result.ival < 1) {
	    sprintf(msg, "Invalid Num_GID_Entries value (%d); "
		"being set to %d.", result.ival, ZOLTAN_NUM_ID_ENTRIES_DEF);
            ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
            result.ival = ZOLTAN_NUM_ID_ENTRIES_DEF;
        }
        lb->Num_GID = result.ival;
        status = 3;
        break;

      case 9:           /* Num_LID_Entries */
        if (result.def)
            result.ival = ZOLTAN_NUM_ID_ENTRIES_DEF;
        if (result.ival < 0) {
	    sprintf(msg, "Invalid Num_LID_Entries value (%d); "
		"being set to %d.", result.ival, ZOLTAN_NUM_ID_ENTRIES_DEF);
            ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
            result.ival = ZOLTAN_NUM_ID_ENTRIES_DEF;
        }
        lb->Num_LID = result.ival;
        status = 3;
        break;

      case 10:          /* LB_Return_Lists */
        if (strcmp(result.sval, "ALL") == 0) {
          tmp = ZOLTAN_LB_ALL_LISTS;
          status = 3;
        }
        else if (strcmp(result.sval, "IMPORT")==0) {
          tmp = ZOLTAN_LB_IMPORT_LISTS;
          status = 3;
        }
        else if (strcmp(result.sval, "EXPORT")==0) {
          tmp = ZOLTAN_LB_EXPORT_LISTS;
          status = 3;
        }
        else if (strcmp(result.sval, "NONE")==0) {
          tmp = ZOLTAN_LB_NO_LISTS;
          status = 3;
        }
        else{
          tmp = ZOLTAN_LB_RETURN_LISTS_DEF;
          sprintf(msg, "Unknown return_lists option %s.", result.sval);
          ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
          status = 2; /* Illegal parameter */
        }
	lb->LB_Return_Lists = tmp;
        break;

      case 11: 		/* Tflops Special flag */
        if (result.def)
            result.ival = ZOLTAN_TFLOPS_SPECIAL_DEF;
	if (result.ival < 0) {
	    sprintf(msg, "Invalid Tflops Special value (%d) "
		"being set to %d.", result.ival, ZOLTAN_TFLOPS_SPECIAL_DEF);
            ZOLTAN_PRINT_WARN(lb->Proc, yo, msg);
	    result.ival = ZOLTAN_TFLOPS_SPECIAL_DEF;
	}
	lb->Tflops_Special = result.ival;
	status = 3;		/* Don't add to Params field of LB */
        break;

      }  /* end switch (index) */
    }

    return(status);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*
 *  Print key parameters.
 */
void Zoltan_Print_Key_Params(LB *lb)
{
  printf("ZOLTAN Parameter %s = %f\n", Key_params[0].name, 
         lb->Imbalance_Tol);
  printf("ZOLTAN Parameter %s = %s\n", Key_params[1].name, 
         (lb->Migrate.Auto_Migrate ? "TRUE" : "FALSE"));
  printf("ZOLTAN Parameter %s = %d\n", Key_params[2].name, 
         lb->Obj_Weight_Dim);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[3].name, 
         lb->Comm_Weight_Dim);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[4].name, 
         lb->Debug_Level);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[5].name, 
         lb->Debug_Proc);
  printf("ZOLTAN Parameter %s = %s\n", Key_params[6].name, 
         (lb->Deterministic ? "TRUE" : "FALSE"));
  printf("ZOLTAN Parameter %s = %d ", Key_params[7].name, lb->Timer);
  if (lb->Timer == ZOLTAN_TIME_WALL)
     printf("(wall)");
  else if (lb->Timer == ZOLTAN_TIME_CPU)
     printf("(cpu)");
  printf("\n");
  printf("ZOLTAN Parameter %s = %d\n", Key_params[8].name, 
         lb->Num_GID);
  printf("ZOLTAN Parameter %s = %d\n", Key_params[9].name, 
         lb->Num_LID);
  printf("ZOLTAN Parameter %s = ", Key_params[10].name);
  switch (lb->LB_Return_Lists) {
  case ZOLTAN_LB_ALL_LISTS:
    printf("ALL\n");
    break;
  case ZOLTAN_LB_IMPORT_LISTS:
    printf("IMPORT\n");
    break;
  case ZOLTAN_LB_EXPORT_LISTS:
    printf("EXPORT\n");
    break;
  case ZOLTAN_LB_NO_LISTS:
    printf("NONE\n");
    break;
  }
  if (lb->Tflops_Special)   /* print only if set */
     printf("ZOLTAN Parameter %s = %s\n", Key_params[11].name, "TRUE");
}
