/*****************************************************************************
 * Zoltan Dynamic Load-Balancing Library for Parallel Applications           *
 * Copyright (c) 2000, Sandia National Laboratories.                         *
 * Zoltan is distributed under the GNU Lesser General Public License 2.1.    * 
 * For more info, see the README file in the top-level Zoltan directory.     *  
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

#ifndef _DR_OUTPUT_CONST_H_
#define _DR_OUTPUT_CONST_H_

#include "dr_input_const.h"

extern void print_distributed_mesh(
  int Proc,
  int Num_Proc,
  MESH_INFO_PTR mesh
);

extern int output_results(
  char *cmd_file,
  int Proc,
  int Num_Proc,
  PROB_INFO_PTR prob,
  PARIO_INFO_PTR pio_info,
  MESH_INFO_PTR mesh);

#endif /* _DR_OUTPUT_CONST_H_ */
