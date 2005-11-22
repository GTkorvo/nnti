#include "TSFAmesosSolver.hpp"
#include "TSFEpetraVector.hpp"
#include "TSFEpetraMatrix.hpp"

#include "Amesos.h"
#include "Amesos_BaseSolver.h"

using namespace TSFExtended;
using namespace Teuchos;


AmesosSolver::AmesosSolver()
  : LinearSolverBase<double>(ParameterList())
{}



SolverState<double> AmesosSolver::solve(const LinearOperator<double>& op, 
                                       const Vector<double>& rhs, 
                                       Vector<double>& soln) const
{
	TSFExtended::Vector<double> bCopy = rhs.copy();
	TSFExtended::Vector<double> xCopy = rhs.copy();

  Epetra_Vector* b = EpetraVector::getConcretePtr(bCopy);
  Epetra_Vector* x = EpetraVector::getConcretePtr(xCopy);

	Epetra_CrsMatrix& A = EpetraMatrix::getConcrete(op);

  Epetra_LinearProblem prob(&A, x, b);

  Amesos amFactory;
  RefCountPtr<Amesos_BaseSolver> solver = rcp(amFactory.Create("Amesos_Umfpack", prob));
  TEST_FOR_EXCEPTION(solver.get()==0, runtime_error, 
                     "AmesosSolver::solve() failed to instantiate Umfpack "
                     "solver");

  int ierr = solver->Solve();
  
  soln = xCopy;

  SolverStatusCode state;
  string msg;

  switch(ierr)
    {
    case 0:
      state = SolveConverged;
      msg = "converged";
      break;
    default:
      state = SolveCrashed;
      msg = "amesos failed: ierr=" + Teuchos::toString(ierr);
    }

  SolverState<double> rtn(state, "Amesos solver " + msg, 0, 0);
  return rtn;
}

