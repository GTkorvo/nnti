
// @HEADER
// ***********************************************************************
// 
//            Trilinos: An Object-Oriented Solver Framework
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// Example of Epetra_MultiVector; use of ExtractView on MultiVector.

#include "Didasko_ConfigDefs.h"
#if defined(HAVE_DIDASKO_EPETRA)

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif  

  // Total number of elements in the vector
  int NumElements = 10;

  // Construct a Map with NumElements and index base of 0
  Epetra_Map Map(NumElements, 0, Comm);

  // Create x as a 2-component multi-vector
  Epetra_MultiVector x(Map,2);
 
  // get the local size of the vector
  int MyLength = x.MyLength();

  /* First way to define the vector:   */
  /* use the [] operator on the object */

  for( int c=0 ; c<x.NumVectors() ; ++c ) 
    for( int i=0 ; i<MyLength ; ++i ) x[c][i] = 1.0*i+1000*c;

  // need a double pointer because this works with multi-vectors
  double ** pointer;
  
  x.ExtractView( &pointer );

  for( int c=0 ; c<x.NumVectors() ;++c ) 
    for( int i=0 ; i<MyLength ; ++i )
      cout << "on proc " << Comm.MyPID() << ", x["
	   << i << "] = " << pointer[c][i] << endl;

  // now modify the values
  for( int c=0 ; c<x.NumVectors() ;++c ) 
    for( int i=0 ; i<MyLength ; ++i )
      pointer[c][i] *= 10;

  cout << x;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(0);
  
} /* main */

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure Didasko with:\n"
       "--enable-epetra");

  return 0;
}

#endif
