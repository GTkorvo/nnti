//@HEADER
// ************************************************************************
//
//               Isorropia: Partitioning and Load Balancing Package
//                 Copyright (2006) Sandia Corporation
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
// Questions? Contact Alan Williams (william@sandia.gov)
//                 or Erik Boman    (egboman@sandia.gov)
//
// ************************************************************************
//@HEADER

#include <Isorropia_configdefs.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Rebalance.hpp>

#include <Teuchos_CommandLineProcessor.hpp>

#include <ispatest_utils.hpp>
#include <ispatest_epetra_utils.hpp>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef HAVE_EPETRA
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#endif

#ifdef HAVE_EPETRA
Epetra_CrsGraph* create_epetra_test_graph_1(int numProcs,
                                            int localProc,
                                            bool verbose);

Epetra_CrsMatrix* create_epetra_test_matrix_1(int numProcs,
                                              int localProc,
                                              bool verbose);

bool test_rebalance_epetra_matrix(int numProcs, int localProc, bool verbose);
bool test_rebalance_epetra_graph(int numProcs, int localProc, bool verbose);
#endif

int main(int argc, char** argv) {
#ifdef HAVE_MPI

  bool verbose = false;
  int numProcs = 1;
  int localProc = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  try {
    verbose = ispatest::set_verbose(localProc, argc, argv);
  }
  catch(std::exception& exc) {
    std::cout << "err, setting verbosity: " << exc.what() << std::endl;
    std::cout << "End Result: TEST FAILED" << std::endl;
    MPI_Finalize();
    return(-1);
  }

  // Only print on proc 0
  if (localProc != 0) {
    verbose = false;
  }

  bool test1_passed = false, test2_passed = false;
#ifdef HAVE_EPETRA
  test1_passed = test_rebalance_epetra_graph(numProcs, localProc, verbose);
  test2_passed = test_rebalance_epetra_matrix(numProcs, localProc, verbose);
#else
  std::cout << "rebalance_1d_default main: currently can only test "
         << "rebalancing with Epetra enabled." << std::endl;
#endif

  if (test1_passed && test2_passed && verbose) {
    std::cout << "rebalance_1d_default main: tests passed."<<std::endl;
  }

#else
  std::cout << "rebalance_1d_default: don't have MPI, can't run test."
            << std::endl;
#endif

  MPI_Finalize();

  return(0);
}

#ifdef HAVE_EPETRA
//-------------------------------------------------------------------
bool test_rebalance_epetra_matrix(int numProcs, int localProc, bool verbose)
{
  bool test_passed = false;
#ifndef HAVE_MPI
  return(test_passed);
#endif

  Epetra_CrsMatrix* input_matrix =
    create_epetra_test_matrix_1(numProcs, localProc, verbose);

  //Call the Isorropia rebalance method allowing default behavior. (I.e.,
  //we won't specify any parameters, weights, etc.) Default behavior should
  //be to balance the matrix so that the number of nonzeros on each processor
  //is roughly equal. i.e., by default, weights for each row are assumed to
  //be the number of nonzeros in that row.

  Epetra_CrsMatrix* balanced_matrix = 0;
  try {
    if (verbose) {
      std::cout << " calling Isorropia::create_balanced_copy(Epetra matrix)..."
                << std::endl;
    }

    balanced_matrix = Isorropia::create_balanced_copy(*input_matrix);
  }
  catch(std::exception& exc) {
    std::cout << "caught exception: " << exc.what() << std::endl;
    return(false);
  }

  //Now check the result matrix and make sure that the number of nonzeros
  //is indeed equal on each processor. (We constructed the input matrix
  //so that a correct rebalancing would result in the same number of
  //nonzeros being on each processor.)
  const Epetra_Map& bal_rowmap = balanced_matrix->RowMap();
  int bal_local_num_rows = bal_rowmap.NumMyElements();

  //count the local nonzeros.
  if (verbose) {
    std::cout << " counting local num-nonzeros for balanced matrix..."
              << std::endl;
  }

  int num_nonzeros = 0;
  for(int i=0; i<bal_local_num_rows; ++i) {
    num_nonzeros += balanced_matrix->NumMyEntries(i);
  }

  const Epetra_Comm& comm = input_matrix->Comm();

  int global_num_nonzeros;
  comm.SumAll(&num_nonzeros, &global_num_nonzeros, 1);

  int avg_nnz_per_proc = global_num_nonzeros/numProcs;

  if (verbose) {
    std::cout << " making sure local num-nonzeros ("<<num_nonzeros
             <<") == global_nnz/numProcs ("<<avg_nnz_per_proc
              << ") on every proc...\n" << std::endl;
  }

  if (num_nonzeros == avg_nnz_per_proc) test_passed = true;

  int local_int_result = test_passed ? 1 : 0;
  int global_int_result;
  comm.MinAll(&local_int_result, &global_int_result, 1);

  delete balanced_matrix;
  delete input_matrix;

  test_passed = global_int_result==1 ? true : false;

  if (!test_passed && verbose) {
    std::cout << "test FAILED!" << std::endl;
  }

  return(test_passed);
}

//-------------------------------------------------------------------
bool test_rebalance_epetra_graph(int numProcs, int localProc, bool verbose)
{
  bool test_passed = false;
#ifndef HAVE_MPI
  return(test_passed);
#endif

  Epetra_CrsGraph* input_graph =
    create_epetra_test_graph_1(numProcs, localProc, verbose);

  //Call the Isorropia rebalance method allowing default behavior. (I.e.,
  //we won't specify any parameters, weights, etc.) Default behavior should
  //be to balance the graph so that the number of nonzeros on each processor
  //is roughly equal. i.e., by default, weights for each row are assumed to
  //be the number of nonzeros in that row.

  Epetra_CrsGraph* balanced_graph = 0;
  try {
    if (verbose) {
      std::cout << " calling Isorropia::create_balanced_copy(Epetra graph)..."
                << std::endl;
    }

    balanced_graph = Isorropia::create_balanced_copy(*input_graph);
  }
  catch(std::exception& exc) {
    std::cout << "caught exception: " << exc.what() << std::endl;
    return(false);
  }

  //Now check the result graph and make sure that the number of nonzeros
  //is indeed equal on each processor. (We constructed the input graph
  //so that a correct rebalancing would result in the same number of
  //nonzeros being on each processor.)
  const Epetra_BlockMap& bal_rowmap = balanced_graph->RowMap();
  int bal_local_num_rows = bal_rowmap.NumMyElements();

  //count the local nonzeros.
  if (verbose) {
    std::cout << " counting local num-nonzeros for balanced graph..."
              << std::endl;
  }

  int num_nonzeros = 0;
  for(int i=0; i<bal_local_num_rows; ++i) {
    num_nonzeros += balanced_graph->NumMyIndices(i);
  }

  const Epetra_Comm& comm = input_graph->Comm();

  int global_num_nonzeros;
  comm.SumAll(&num_nonzeros, &global_num_nonzeros, 1);

  int avg_nnz_per_proc = global_num_nonzeros/numProcs;

  if (verbose) {
    std::cout << " making sure local num-nonzeros ("<<num_nonzeros
             <<") == global_nnz/numProcs ("<<avg_nnz_per_proc
              << ") on every proc...\n" << std::endl;
  }

  if (num_nonzeros == avg_nnz_per_proc) test_passed = true;

  int local_int_result = test_passed ? 1 : 0;
  int global_int_result;
  comm.MinAll(&local_int_result, &global_int_result, 1);

  delete balanced_graph;
  delete input_graph;

  test_passed = global_int_result==1 ? true : false;

  if (!test_passed && verbose) {
    std::cout << "test FAILED!" << std::endl;
  }

  return(test_passed);
}

Epetra_CrsMatrix* create_epetra_test_matrix_1(int numProcs,
                                              int localProc,
                                              bool verbose)
{
  if (verbose) {
    std::cout << " creating Epetra_CrsMatrix with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsMatrix with rows spread un-evenly over
  //processors.
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int local_num_rows = 800;
  int nnz_per_row = local_num_rows/4+1;
  int global_num_rows = numProcs*local_num_rows;

  int mid_proc = numProcs/2;
  bool num_procs_even = numProcs%2==0 ? true : false;

  int adjustment = local_num_rows/2;

  //adjust local_num_rows so that it's not equal on all procs.
  if (localProc < mid_proc) {
    local_num_rows -= adjustment;
  }
  else {
    local_num_rows += adjustment;
  }

  //if numProcs is not an even number, undo the local_num_rows adjustment
  //on one proc so that the total will still be correct.
  if (localProc == numProcs-1) {
    if (num_procs_even == false) {
      local_num_rows -= adjustment;
    }
  }

  //now we're ready to create a row-map.
  Epetra_Map rowmap(global_num_rows, local_num_rows, 0, comm);

  //create a matrix
  Epetra_CrsMatrix* input_matrix =
    new Epetra_CrsMatrix(Copy, rowmap, nnz_per_row);

  int err = ispatest::fill_matrix(*input_matrix, nnz_per_row, verbose);
  if (err != 0) {
    delete input_matrix;
    return(0);
  }

  return(input_matrix);
}

Epetra_CrsGraph* create_epetra_test_graph_1(int numProcs,
                                              int localProc,
                                              bool verbose)
{
  if (verbose) {
    std::cout << " creating Epetra_CrsGraph with un-even distribution..."
            << std::endl;
  }

  //create an Epetra_CrsGraph with rows spread un-evenly over
  //processors.
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int local_num_rows = 800;
  int nnz_per_row = local_num_rows/4+1;
  int global_num_rows = numProcs*local_num_rows;

  int mid_proc = numProcs/2;
  bool num_procs_even = numProcs%2==0 ? true : false;

  int adjustment = local_num_rows/2;

  //adjust local_num_rows so that it's not equal on all procs.
  if (localProc < mid_proc) {
    local_num_rows -= adjustment;
  }
  else {
    local_num_rows += adjustment;
  }

  //if numProcs is not an even number, undo the local_num_rows adjustment
  //on one proc so that the total will still be correct.
  if (localProc == numProcs-1) {
    if (num_procs_even == false) {
      local_num_rows -= adjustment;
    }
  }

  //now we're ready to create a row-map.
  Epetra_Map rowmap(global_num_rows, local_num_rows, 0, comm);

  //create a graph
  Epetra_CrsGraph* input_graph =
    new Epetra_CrsGraph(Copy, rowmap, nnz_per_row);

  int err = ispatest::fill_graph(*input_graph, nnz_per_row, verbose);
  if (err != 0) {
    delete input_graph;
    return(0);
  }

  return(input_graph);
}

#endif //HAVE_EPETRA

