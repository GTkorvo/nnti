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

#ifndef __ZOLTAN_UTIL_CONST_H
#define __ZOLTAN_UTIL_CONST_H

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

extern void Zoltan_Get_Obj_List(LB *, ZOLTAN_ID_PTR, ZOLTAN_ID_PTR, int, float *, int *);
extern int LB_pad_for_alignment(int);
extern unsigned int Zoltan_Hash(ZOLTAN_ID_PTR, int, unsigned int);
extern int Zoltan_Clean_String(char *, char **);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#endif
