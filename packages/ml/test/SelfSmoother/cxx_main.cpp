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

#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_IFPACK) && defined(HAVE_MPI) && defined(FIXME

#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"
#include "Trilinos_Util_CrsMatrixGallery.h"
#include "AztecOO.h"
#include "Ifpack_Amesos.h"
#include "ml_Ifpack_ML.h"
#include "ml_MultiLevelPreconditioner.h"

using namespace Teuchos;
using namespace Trilinos_Util;

// ====================================================================== 

template<class T>
int TestAdditiveSchwarz()
{
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 900);
  Epetra_RowMatrix* A = Gallery.GetMatrix();

  Epetra_Vector LHS(A->OperatorDomainMap());
  Epetra_Vector RHS(A->OperatorRangeMap());

  LHS.PutScalar(0.0);
  RHS.Random();

  Epetra_LinearProblem problem(A, &LHS, &RHS);

  AztecOO solver(problem);

  Teuchos::ParameterList List;
  ML_Epetra::SetDefaults("SA", List);
  List.set("output", 0);
  List.set("schwarz: combine mode", Add);
  List.set("cycle applications", 10);

  int Overlap = 2;
  Ifpack_AdditiveSchwarz<T> Prec(A, Overlap);
  Prec.SetParameters(List);
  ML_CHK_ERR(Prec.Initialize());
  ML_CHK_ERR(Prec.Compute());

  solver.SetPrecOperator(&Prec);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

  return(solver.NumIters());
}

// ====================================================================== 
void TestML() 
{

  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("problem_size", 900);

  Epetra_RowMatrix * A = Gallery.GetMatrix();
  Epetra_LinearProblem * Problem = Gallery.GetLinearProblem();

  AztecOO solver(*Problem);

  ParameterList MLList;
  ML_Epetra::SetDefaults("SA",MLList);

  MLList.set("max levels",5);
  MLList.set("increasing or decreasing","decreasing");
  MLList.set("aggregation: type", "Uncoupled");
  MLList.set("smoother: type","self");
  MLList.set("output", 0);

  Teuchos::ParameterList& SelfList = MLList.sublist("smoother: self list");
  ML_Epetra::SetDefaults("DD-ML", SelfList);
  SelfList.set("output", 0);
  SelfList.set("coarse: max size", 128);
  SelfList.set("cycle applications", 1);
  SelfList.set("aggregation: damping factor", 0.0);
  // uncomment the following to use Amesos locally
  //MLList.set("smoother: type", "IFPACK");
  //MLList.set("smoother: ifpack type", "Amesos");
  //MLList.set("smoother: ifpack overlap", 0);

  MLList.set("smoother: pre or post", "both");

#ifdef HAVE_ML_AMESOS
  MLList.set("coarse: type","Amesos-KLU");
#else
  MLList.set("aggregation: type", "MIS");
  MLList.set("smoother: type","Jacobi");
  MLList.set("coarse: type","Jacobi");
#endif

  ML_Epetra::MultiLevelPreconditioner* MLPrec = 
    new ML_Epetra::MultiLevelPreconditioner(*A, MLList);

  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  solver.SetAztecOption(AZ_output, 32);
  solver.Iterate(500, 1e-12);

  delete MLPrec;

  double residual, diff;
  Gallery.ComputeResidual(&residual);
  Gallery.ComputeDiffBetweenStartingAndExactSolutions(&diff);

  if (Comm.MyPID() == 0) {
    cout << "||b-Ax||_2 = " << residual << endl;
    cout << "||x_exact - x||_2 = " << diff << endl;
  }

  if (residual > 1e-5)
    exit(EXIT_FAILURE);
}

// ============== //
// example driver //
// ============== //

int main(int argc, char *argv[])
{
  
  MPI_Init(&argc,&argv);
  int MyPID;
  MPI_Comm_rank(MPI_COMM_WORLD, &MyPID);

  int ItersML = TestAdditiveSchwarz<Ifpack_Amesos>();
  int ItersAm = TestAdditiveSchwarz<ML_Epetra::Ifpack_ML>();

  if (MyPID == 0)
    cout << "ML iterations = " << ItersML
         << ", Amesos iterations = " << ItersAm << endl;

  int diff = ItersML - ItersAm;
  if (diff < 0) diff = -diff;
  if (diff > 5)
    exit(EXIT_FAILURE);
    
  TestML();
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  exit(EXIT_SUCCESS);
 
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with:");
  puts("--enable-epetra");
  puts("--enable-teuchos");
  puts("--enable-aztecoo");
  puts("--enable-triutils");
  puts("--enable-ifpack");
  puts("--enable-mpi");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  exit(EXIT_SUCCESS);
}

#endif /* #if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_TRIUTILS) && defined(HAVE_ML_AZTECOO) */
