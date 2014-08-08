
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "Epetra_LinearProblem.h"
#include "Trilinos_Util_CrsMatrixGallery.h"

using namespace Trilinos_Util;

int main(int argc, char *argv[]) 
{
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  int nx = 4000;
  int ny = 1000 * Comm.NumProc();

  CrsMatrixGallery Gallery("laplace_2d", Comm);
  Gallery.Set("ny", ny);
  Gallery.Set("nx", nx);
  Gallery.Set("problem_size",nx*ny);
  Gallery.Set("map_type", "linear");

  Epetra_LinearProblem* Problem = Gallery.GetLinearProblem();
  assert (Problem != 0);
  // retrive pointers to solution (lhs), right-hand side (rhs)
  // and matrix itself (A)
  Epetra_MultiVector* lhs = Problem->GetLHS();
  Epetra_MultiVector* rhs = Problem->GetRHS();
  Epetra_RowMatrix* A = Problem->GetMatrix();
    
  Epetra_Time Time(Comm);

  for (int i = 0 ; i < 10 ; ++i)
    A->Multiply(false, *lhs, *rhs);

  cout << Time.ElapsedTime() << endl;

  MPI_Finalize();

  return(EXIT_SUCCESS);
} // end of main()
