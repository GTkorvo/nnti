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

void  scscmv (int isym, int m, int n, 
	      double *val, int *indx, int *pntr, 
	      double *x, double *y)
{
    int i, j, jbgn, jend;


/*     Performs the matrix-vector operation

                      y = A*x 

       where x and y are vectors and A is a sparse matrix stored
       in (Harwell-Boeing) compress column format. */

/*     -------------------------- 
       First executable statement 
       -------------------------- */

/* .....initialize soln */

    for (i = 0; i < m; i++)
	y[i] = 0.0;

/* .....do a series of SPAXPYs (sparse saxpys) */

    for (j = 0; j <n ; j++) 
      {
	jbgn = pntr[j];
	jend = pntr[j + 1];

	for (i = jbgn; i < jend; i++)
	  {
	    y[indx[i]] += val[i] * x[j];
	    if (indx[i] != j && isym) y[j] += val[i]*x[indx[i]];
	  }
      }

    return;
} /* scscmv */

