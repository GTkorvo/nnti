#include "TSFBelosSolver.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosThyraAdapter.hpp"
#include "TSFPreconditioner.hpp"
#include "TSFPreconditionerFactory.hpp"
#include "TSFParameterListPreconditionerFactory.hpp"


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANTIATION
#include "TSFVectorImpl.hpp"
#include "TSFLinearOperatorImpl.hpp"
#include "TSFLinearSolverImpl.hpp"
#endif

using namespace TSFExtended;
using namespace Teuchos;


BelosSolver::BelosSolver(const ParameterList& params)
  : LinearSolverBase<double>(params), pf_()
{
  setName("BelosSolver");
  if (params.isSublist("Preconditioner"))
  {
    ParameterList precParams = params.sublist("Preconditioner");
    pf_ = new ParameterListPreconditionerFactory(precParams);
  }
}



SolverState<double> BelosSolver::solve(const LinearOperator<double>& A, 
  const Vector<double>& rhs, 
  Vector<double>& soln) const
{
  typedef Thyra::MultiVectorBase<double>         MV;
  typedef Thyra::LinearOpBase<double>            OP;
  typedef Belos::LinearProblem<double, MV, OP>   LP;

  /* get Thyra objects */
  RCP<OP> APtr = A.ptr();
  RCP<MV> bPtr = rhs.ptr(); 
  RCP<MV> ansPtr = soln.ptr();
  
  RCP<LP> prob = rcp(new LP(APtr, ansPtr, bPtr));

  TEST_FOR_EXCEPT(!prob->setProblem());

  
  if (pf_.ptr().get())
  {
    Preconditioner<double> P = pf_.createPreconditioner(A);
    if (P.hasLeft())
    {
      prob->setLeftPrec(P.left().ptr());
    }
  
    if (P.hasRight())
    {
      prob->setRightPrec(P.right().ptr());
    }
  }

  ParameterList plist = parameters();
  Out::os() << "Belos parameters = " << plist << endl;
  RCP<ParameterList> belosList = rcp(&plist, false);
  RCP<Belos::SolverManager<double, MV, OP> > solver ;
  string solverType = parameters().get<string>("Method");
  
  if (solverType=="GMRES")
  {
    solver=rcp(new Belos::BlockGmresSolMgr<double, MV, OP>(prob, belosList));
  }
  else if (solverType=="CG")
  {
    solver=rcp(new Belos::BlockCGSolMgr<double, MV, OP>(prob, belosList));
  }
  else
  {
    TEST_FOR_EXCEPT(!(solverType=="GMRES" || solverType=="CG"));
  }

  Out::os() << "Belos thinks that parameters are " << *(solver->getCurrentParameters()) << endl;

  Belos::ReturnType rtn = solver->solve();

  int numIters = solver->getNumIters();
  cerr << "num iters=" << numIters << endl;
  double resid = -1.0;
  
  SolverStatusCode code = SolveFailedToConverge;
  if (rtn==Belos::Converged) code = SolveConverged;
  else cerr << "solver failed to converge!" << endl;
  SolverState<double> state(code, "Belos solver completed", numIters, resid);
  
  return state;
}



