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

#ifndef __PARMETIS_JOSTLE_H
#define __PARMETIS_JOSTLE_H

#include <limits.h>
#include "comm_const.h"
#include "parmetis_jostle_const.h"

/* Do we have the multiconstraint beta version of ParMetis? */
/* #define BETA_PARMETIS */

/* Include ParMetis and/or Jostle header files if necessary. 
 * These include files must be available in the include path set in the 
 * Zoltan configuration file.
 */
#ifdef ZOLTAN_PARMETIS
#include "parmetis.h"
#else 
typedef int idxtype; 
#endif

#ifdef ZOLTAN_JOSTLE
#include "jostle.h"
#endif

/* ParMetis option defs. These must be identical to the defs
 * in defs.h in the version of ParMetis you are using!
 */
#define OPTION_IPART            1
#define OPTION_FOLDF            2
#define OPTION_DBGLVL           3
#define MAX_OPTIONS             4  /* Total number of options +1 */

/* Misc. defs to be used with MPI */
#define TAG1  32001
#define TAG2  32002
#define TAG3  32003
#define TAG4  32004
#define TAG5  32005

/* Misc. local constants */
#define CHUNKSIZE 20  /* Number of nodes to allocate in one chunk. */


/* ParMETIS data types and definitions. */

/* Undefine the following #define in order to use short as the idxtype.
 * Make sure these defs are consistent with those in your 
 * ParMetis installation ! It is strongly recommended to use 
 * integers, not shorts, if you load balance with weights.
*/

#ifdef IDXTYPE_IS_SHORT
/* typedef short idxtype; This should have been done in parmetis.h */
#define IDX_DATATYPE    MPI_SHORT
#define MAX_WGT_SUM (SHRT_MAX/8)
#else /* the default for idxtype is int; this is recommended */
/* typedef int idxtype; This should have been done in parmetis.h */
#define IDX_DATATYPE    MPI_INT
#define MAX_WGT_SUM (INT_MAX/8)
#endif

extern int Zoltan_Verify_Graph(MPI_Comm comm, idxtype *vtxdist, idxtype *xadj,
              idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt,
              int vwgt_dim, int ewgt_dim, int check_graph, int debug_level);
extern int Zoltan_Scatter_Graph(int have_graph, idxtype **vtxdist, idxtype **xadj, idxtype **adjncy,
              idxtype **vwgt, idxtype **adjwgt, float   **xyz, int     ndims,
              LB      *lb, ZOLTAN_COMM_OBJ **plan);

#endif
