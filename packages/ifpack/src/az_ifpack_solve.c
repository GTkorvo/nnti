/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "az_aztec.h"
#include "az_ifpack.h"

void AZ_ifpack_solve(double x[], double b[], int options[], double params[],
              int indx[], int bindx[], int rpntr[], int cpntr[], int bpntr[],
              double val[], int data_org[], double status[], int proc_config[])
{
  AZ_MATRIX *Amat;

   Amat    = AZ_matrix_create(data_org[AZ_N_internal]+data_org[AZ_N_border]);

   options[AZ_output] = 1;
   if (data_org[AZ_matrix_type] == AZ_MSR_MATRIX) 
      AZ_set_MSR(Amat, bindx, val, data_org, 0, NULL, AZ_LOCAL);
   else if (data_org[AZ_matrix_type] == AZ_VBR_MATRIX)
      AZ_set_VBR(Amat, rpntr, cpntr, bpntr, indx, bindx, val,
                 data_org, 0, NULL, AZ_LOCAL);
   else {
      fprintf(stderr,"Unknown matrix type (%d)\n",data_org[AZ_matrix_type]);
      fprintf(stderr,"Matrix-free is now available via AZ_iterate()\n");
      exit(1);
   }

  AZ_ifpack_iterate(x, b, options, params, status, proc_config, Amat);

  AZ_matrix_destroy(&Amat);

}
/* AZ_ifpack_solve*/
