
//@HEADER
// ************************************************************************
// 
//               ML: A Multilevel Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"
#include "ml_include.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Space.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_Preconditioner.h"
#include "MLAPI_EpetraPreconditioner.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"

using namespace Teuchos;
using namespace Trilinos_Util;
using namespace MLAPI;

// ======================================= //
// 2-level additive Schwarz preconditioner //
// ======================================= //

class TwoLevelDDAdditive : public Preconditioner {

public:
  // Constructor assumes that all operators and inverse operators are already
  // filled.
  TwoLevelDDAdditive(const Operator FineMatrix, 
                     const InverseOperator FineSolver, 
                     const InverseOperator CoarseSolver, 
                     const Operator R, 
                     const Operator P) :
    FineMatrix_(FineMatrix),
    R_(R),
    P_(P),
    FineSolver_(FineSolver),
    CoarseSolver_(CoarseSolver)
  {}
      
  int Solve(const DoubleVector& r_f, DoubleVector& x_f) const
  {
    
    DoubleVector r_c(FineSolver_.DomainSpace());

    // apply fine level preconditioner
    x_f = FineSolver_ / r_f;
    // restrict to coarse
    r_c = R_ * r_f;
    // solve coarse problem
    r_c = CoarseSolver_ / r_c;
    // prolongate back and add to solution
    x_f = x_f + P_ * r_c;

    return(0);
  }

  const Space DomainSpace() const {
    return(FineMatrix_.DomainSpace());
  }

  const Space RangeSpace() const {
    return(FineMatrix_.RangeSpace());
  }

private:
  const Operator FineMatrix_;
  const Operator R_;
  const Operator P_;
  const InverseOperator FineSolver_;
  const InverseOperator CoarseSolver_;

}; // TwoLevelDDAdditive

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_Comm Comm(MPI_COMM_WORLD);
#endif

  int NumGlobalElements = 10000;
  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", NumGlobalElements);
  Epetra_RowMatrix& A = *(Gallery.GetMatrix());
  int NumMyElements = A.NumMyRows();

  try {

    // Initialize the workspace and set the output level
    Init();
    SetPrintLevel(10);

    Teuchos::ParameterList MLList;

    // define the space for fine level vectors and operators.
    Space FineSpace(NumMyElements);
    Operator FineMatrix(FineSpace,FineSpace,A);

    // CoarseMatrix will contain the coarse level matrix,
    // while FineSolver and CoarseSolver will contain
    // the fine level smoother and the coarse level solver,
    // respectively. We will use symmetric Gauss-Seidel
    // for the fine level, and Amesos (LU) for the coarse level.
    Operator CoarseMatrix;
    InverseOperator FineSolver, CoarseSolver;

    // Now we define restriction (R) and prolongator (P) from the fine space
    // to the coarse space using non-smoothed aggregation.
    // The coarse-level matrix will be defined via a triple
    // matrix-matrix product.
    Operator R, P;

#ifdef HAVE_ML_METIS
    MLList.set("aggregation: type","METIS");
    MLList.set("aggregation: nodes per aggregate",64);
#else
    MLList.set("aggregation: type","Uncoupled");
#endif

    P = BuildP(FineMatrix,MLList);
    R = Transpose(P);

    CoarseMatrix = RAP(R,FineMatrix,P);
    FineSolver.Reshape(FineMatrix,"SGS",MLList);
    CoarseSolver.Reshape(CoarseMatrix,"Amesos",MLList);

    // We can now construct a Preconditioner-derived object, that
    // implements the 2-level hybrid domain decomposition preconditioner.
    // Preconditioner `TwoLevelDDHybrid' can be replaced by
    // `TwoLevelDDAdditive' to define an purely additive preconditioner.
    TwoLevelDDAdditive MLAPIPrec(FineMatrix,FineSolver,CoarseSolver,R,P);

    // Define an AztecOO solver, wrap MLAPIPrec as an Epetra_Operator,
    // then solve the problem
    AztecOO solver(*Gallery.GetLinearProblem());

    EpetraPreconditioner EpetraPrec(A.RowMatrixRowMap(),MLAPIPrec);

    solver.SetPrecOperator(&EpetraPrec);

    solver.SetAztecOption(AZ_solver, AZ_gmres);
    solver.SetAztecOption(AZ_output, 16);

    // solve with 500 iterations and 1e-12 tolerance  
    // The problem should converge as follows:
    solver.Iterate(500, 1e-5);

  }
  catch (const char e[]) {
    cerr << "Caught exception: " << e << endl;
  }
  catch (...) {
    cerr << "Caught exception..." << endl;
  }

#ifdef HAVE_MPI
  // finalize the MLAPI workspace
  Finalize();
  MPI_Finalize() ;
#endif

  return 0 ;
  
}

#else

#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-triutils");
  
  return 0;
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) */
