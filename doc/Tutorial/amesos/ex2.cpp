
//@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2001) Sandia Corporation
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

// Trilinos Tutorial
// -----------------
// Solve a linear system with Amesos
//
// Marzio Sala, SNL, 9214, 19-Nov-2003

#include "Epetra_config.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_IntVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "Amesos_Umfpack.h"
#include "Amesos_Parameter_List.h"
#include "Amesos_Superludist.h"

void  get_neighbours( const int i, const int nx, const int ny,
		      int & left, int & right, 
		      int & lower, int & upper);

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int nx = 5;
  int ny = 6;
  int NumGlobalElements = nx * ny;

  // create a map
  Epetra_Map Map(NumGlobalElements,0,Comm);
  // local number of rows
  int NumMyElements = Map.NumMyElements();
  // get update list
  int * MyGlobalElements = new int [NumMyElements];
  Map.MyGlobalElements( MyGlobalElements );

  // Create an integer vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor

  int * NumNz = new int[NumMyElements];
  
  double off_left  = -1.0;
  double off_right = -1.0;
  double off_lower = -1.0;
  double off_upper = -1.0;
  double diag      =  4.0;
  int left, right, lower, upper;
  
  for ( int i=0; i<NumMyElements; i++) {
    NumNz[i] = 1;
    get_neighbours( MyGlobalElements[i], nx, ny, 
		    left, right, lower, upper); 
    if( left  != -1 ) ++NumNz[i];
    if( right != -1 ) ++NumNz[i];
    if( lower != -1 ) ++NumNz[i];
    if( upper != -1 ) ++NumNz[i];
  }
  
  // Create a Epetra_Matrix
  // create a CRS matrix

  Epetra_CrsMatrix A(Copy,Map,NumNz);

  // Add  rows one-at-a-time
  // Need some vectors to help
  // Off diagonal Values will always be -1

  double *Values = new double[4];
  int *Indices = new int[4];
  int NumEntries;

  for( int i=0 ; i<NumMyElements; ++i ) {
    int NumEntries=0;
    get_neighbours(  MyGlobalElements[i], nx, ny, 
		     left, right, lower, upper);
    if( left != -1 ) {
	Indices[NumEntries] = left;
	Values[NumEntries] = off_left;
	++NumEntries;
    }
    if( right != -1 ) {
      Indices[NumEntries] = right;
      Values[NumEntries] = off_right;
      ++NumEntries;
    }
    if( lower != -1 ) {
      Indices[NumEntries] = lower;
      Values[NumEntries] = off_lower;
      ++NumEntries;
    }
    if( upper != -1 ) {
      Indices[NumEntries] = upper;
      Values[NumEntries] = off_upper;
      ++NumEntries;
    }
    // put the off-diagonal entries
    assert(A.InsertGlobalValues(MyGlobalElements[i], NumEntries, 
				Values, Indices)==0);
    // Put in the diagonal entry
    assert(A.InsertGlobalValues(MyGlobalElements[i], 1, 
				&diag, MyGlobalElements+i)==0);
  }
  
  // Finish up
  assert(A.TransformToLocal()==0);

  // create x and b vectors
  Epetra_Vector x(Map);
  Epetra_Vector b(Map);
  x.PutScalar(1.0);

  A.Multiply(false,x,b);

  x.PutScalar(0.0);
  
  // create linear problem
  Epetra_LinearProblem Problem(&A,&x,&b);

  AMESOS::Parameter::List params;
  
  Amesos_Superludist * SuperludistProblem =
    new Amesos_Superludist(Problem,params);
   
  SuperludistProblem->SymbolicFactorization();
  SuperludistProblem->NumericFactorization();
  SuperludistProblem->Solve();
  cout << x;

  delete SuperludistProblem;
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return( EXIT_SUCCESS );

}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void  get_neighbours( const int i, const int nx, const int ny,
		      int & left, int & right, 
		      int & lower, int & upper) 
{

  int ix, iy;
  ix = i%nx;
  iy = (i - ix)/nx;

  if( ix == 0 ) 
    left = -1;
  else 
    left = i-1;
  if( ix == nx-1 ) 
    right = -1;
  else
    right = i+1;
  if( iy == 0 ) 
    lower = -1;
  else
    lower = i-nx;
  if( iy == ny-1 ) 
    upper = -1;
  else
    upper = i+nx;

  return;

}
