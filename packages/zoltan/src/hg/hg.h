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

#ifndef __HG_H
#define __HG_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "zz_const.h"
#include "hg_const.h"
#include "hypergraph.h"



/*****************************************************************************/
/* Data structure for Zoltan's base hypergraph.
 * Includes Zoltan IDs corresponding to local objects (vertices) and
 * a HGraph as used by the algorithms. */

struct Zoltan_HGraph {
  ZOLTAN_ID_PTR Global_IDs; /* Global IDs for on-processor objects.  */
  ZOLTAN_ID_PTR Local_IDs;  /* Local IDs for on-processor objects.   */
  Partition Parts;          /* Initial partition #s for on-processor objects */
                            /* KDD In parallel version Part may be part of HG */
  HGraph HG;                /* Hypergraph for initial objects.       */
};
typedef struct Zoltan_HGraph ZHG;

/*****************************************************************************/



/* Prototypes */
extern int Zoltan_HG_Build_Hypergraph (ZZ*, ZHG**, HGPartParams*);
extern void Zoltan_HG_HGraph_Print(ZZ*, ZHG*,  HGraph*, FILE*);



#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif   /* __HG_H */
