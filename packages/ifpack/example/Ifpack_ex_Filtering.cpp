// @HEADER
// ***********************************************************************
// 
//                IFPACK
//                 Copyright (2004) Sandia Corporation
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

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Ifpack_DropFilter.h"
#include "Ifpack_SparsityFilter.h"
#include "Ifpack_SingletonFilter.h"
#include "Ifpack_Utils.h"

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  if (Comm.NumProc() != 1) {
    cerr << "This example must be run with one process only." << endl;
    // exit with success not to break the test harness
    exit(EXIT_SUCCESS);
  }
  
  int NumPoints = 5;
  Epetra_Map Map(NumPoints,0,Comm);

  Epetra_CrsMatrix Matrix(Copy,Map,0);

  vector<int> Indices(NumPoints);
  vector<double> Values(NumPoints);
  double Diag = 0.0;

  for (int i = 0 ; i < NumPoints ; ++i) {
    // add a diagonal
    Matrix.InsertGlobalValues(i,1,&Diag,&i);

    // add off-diagonals
    int NumEntries = 0;
    for (int j = i + 1 ; j < NumPoints ; ++j) {
      Indices[NumEntries] = j;
      Values[NumEntries] = 1.0 * (j - i);
      ++NumEntries;
    }
    Matrix.InsertGlobalValues(i,NumEntries,&Values[0],&Indices[0]);
  }
  Matrix.FillComplete();

  // ================================= //
  // print sparsity of original matrix //
  // ================================= //
 
  cout << "Sparsity, non-dropped matrix" << endl;
  Ifpack_PrintSparsity_Simple(Matrix);

  // ====================================== //
  // create a new matrix, dropping by value //
  // ====================================== //
  //
  // drop all elements below 4.0. Only the upper-right element
  // is maintained, plus all the diagonals that are not
  // considering in dropping.
  Ifpack_DropFilter DropA(&Matrix,4.0);
  assert (DropA.MaxNumEntries() == 2);

  cout << "Sparsity, dropping by value" << endl;
  Ifpack_PrintSparsity_Simple(DropA);

  // ========================================= //
  // create a new matrix, dropping by sparsity //
  // ========================================= //
  //
  // Mantain 2 off-diagonal elements.
  Ifpack_SparsityFilter SparsityA(&Matrix,2);

  cout << "Sparsity, dropping by sparsity" << endl;
  Ifpack_PrintSparsity_Simple(SparsityA);
  assert (SparsityA.MaxNumEntries() == 3);

  // ======================================== //
  // create new matrices, dropping singletons //
  // ======================================== //
  //
  // If we apply this filter NumPoints - 1 times, 
  // we end up with a one-row matrix
  Ifpack_SingletonFilter Filter1(&Matrix);
  Ifpack_SingletonFilter Filter2(&Filter1);
  Ifpack_SingletonFilter Filter3(&Filter2);
  Ifpack_SingletonFilter Filter4(&Filter3);

  cout << "Sparsity, dropping singletons 4 times" << endl;
  Ifpack_PrintSparsity_Simple(Filter4);
  assert (Filter4.NumMyRows() == 1);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  return(EXIT_SUCCESS);
}

