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
#if defined(HAVE_IFPACK_AZTECOO) && defined(HAVE_IFPACK_AMESOS) && defined(HAVE_IFPACK_TEUCHOS)
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_BlockRelaxation.h"
#include "Ifpack_SparseContainer.h"
#include "Ifpack_Amesos.h"

using namespace Trilinos_Util;

int main(int argc, char *argv[])
{

  // initialize MPI and Epetra communicator
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  // size of the global matrix (must be a square number)
  const int NumPoints = 100;

  // build the matrix corresponding to a 2D Laplacian on a
  // structured grid.
  CrsMatrixGallery Gallery("laplace_2d_bc", Comm);
  Gallery.Set("problem_size", NumPoints);
  // for simplicity, linear map.
  Gallery.Set("map_type", "linear");

  // get the pointer to the linear system matrix
  Epetra_RowMatrix* A = Gallery.GetMatrix();

// cout << *(dynamic_cast<Epetra_CrsMatrix*>(A));
//  exit(0);
  // =============================================================== //
  // B E G I N N I N G   O F   I F P A C K   C O N S T R U C T I O N //
  // =============================================================== //

  Teuchos::ParameterList List;

  // builds an Ifpack_AdditiveSchwarz. This is templated with
  // the local solvers, in this case Ifpack_BlockRelaxation.
  // Ifpack_BlockRelaxation requires as a templated a container
  // class. A container defines
  // how to store the diagonal blocks. Two choices are avaiable:
  // Ifpack_DenseContainer (to store them as dense block,
  // than use LAPACK' factorization to apply the inverse of
  // each block), of Ifpack_SparseContainer (to store
  // the diagonal block as Epetra_CrsMatrix's). 
  // 
  // Here, we use Ifpack_SparseContainer, which in turn is
  // templated with the class to use to apply the inverse
  // of each block. For example, we can use Ifpack_Amesos.
 
  // We still have to decide the overlap among the processes,
  // and the overlap among the blocks. The two values
  // can be different.
  int OverlapProcs = 2;
  int OverlapBlocks = 0;

  Ifpack_AdditiveSchwarz<Ifpack_BlockRelaxation<Ifpack_SparseContainer<Ifpack_Amesos> > > Prec(A, OverlapProcs);

  // other options are "Jacobi", "symmetric Gauss-Seidel"
  List.set("block: type", "Gauss-Seidel");
  List.set("block: overlap", OverlapBlocks);
  // use METIS to create the blocks. This requires --enable-ifpack-metis.
  // If METIS is not installed, the user may select "linear". A
  // simple greedy algorithm can be enabled using option "greedy
  List.set("partitioner: type", "metis");
  List.set("partitioner: use symmetric graph", true);
  // defines here the number of local blocks. If 1,
  // and only one process is used in the computation, then
  // the preconditioner must converge in one iteration.
  List.set("partitioner: local parts", 4);

  // sets the parameters
  IFPACK_CHK_ERR(Prec.SetParameters(List));

  // initialize the preconditioner. 
  IFPACK_CHK_ERR(Prec.Initialize());

  // Builds the preconditioners.
  IFPACK_CHK_ERR(Prec.Compute());

  // =================================================== //
  // E N D   O F   I F P A C K   C O N S T R U C T I O N //
  // =================================================== //

  // At this point, we need some additional objects
  // to define and solve the linear system.

  // defines LHS and RHS
  Epetra_Vector LHS(A->OperatorDomainMap());
  Epetra_Vector RHS(A->OperatorDomainMap());

  LHS.PutScalar(0.0);
  RHS.Random();

  // need an Epetra_LinearProblem to define AztecOO solver
  Epetra_LinearProblem Problem(A,&LHS,&RHS);

  // now we can allocate the AztecOO solver
  AztecOO Solver(Problem);

  // specify solver
  Solver.SetAztecOption(AZ_solver,AZ_gmres);
  Solver.SetAztecOption(AZ_output,32);

  // HERE WE SET THE IFPACK PRECONDITIONER
  Solver.SetPrecOperator(&Prec);

  // .. and here we solve
  // NOTE: with one process, the solver must converge in
  // one iteration.
  Solver.Iterate(1550,1e-8);

#ifdef HAVE_MPI
  MPI_Finalize() ; 
#endif

    return(EXIT_SUCCESS);
}

#else

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm Comm;
#endif

  puts("please configure IFPACK with --eanble-aztecoo --enable-teuchos");
  puts("--enable-amesos to run this test");

#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif
  return(EXIT_SUCCESS);
}

#endif
