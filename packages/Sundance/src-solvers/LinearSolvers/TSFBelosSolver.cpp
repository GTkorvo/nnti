#include "TSFBelosSolver.hpp"
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
  : LinearSolverBase<double>(params), pf_(), hasSolver_(false)
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

  TEST_FOR_EXCEPT(!A.ptr().get());
  TEST_FOR_EXCEPT(!rhs.ptr().get());

  /* get Thyra objects */
  RCP<OP> APtr = A.ptr();
  RCP<MV> bPtr = rhs.ptr(); 

  if (!soln.ptr().get()) soln = rhs.copy();

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

  if (!hasSolver_)
    {

      ParameterList plist = parameters();

      RCP<ParameterList> belosList = rcp(&plist, false);

      std::string solverType = parameters().get<string>("Method");
      
      if (solverType=="GMRES")
	{
	  solver_=rcp(new Belos::BlockGmresSolMgr<double, MV, OP>(prob, belosList));
	}
      else if (solverType=="CG")
	{
	  solver_=rcp(new Belos::BlockCGSolMgr<double, MV, OP>(prob, belosList));
	}
      else if (solverType=="TFQMR")
	{
	  solver_=rcp(new Belos::TFQMRSolMgr<double, MV, OP>(prob, belosList));
	}
      else if (solverType=="GCRODR")
	{
	  solver_=rcp(new Belos::GCRODRSolMgr<double, MV, OP>(prob, belosList));
	  hasSolver_ = true; // only cache recycling solvers
	}
      else if (solverType=="RCG")
	{
	  solver_=rcp(new Belos::RCGSolMgr<double, MV, OP>(prob, belosList));
	  hasSolver_ = true; // only cache recycling solvers
	}
      else
	{
	  TEST_FOR_EXCEPT(!(solverType=="GMRES" || solverType=="CG"));
	}
    }
  else // reset problem
    {
      solver_->setProblem( prob );
    }
  
  Belos::ReturnType rtn = solver_->solve();

  int numIters = solver_->getNumIters();
  double resid = -1.0;
  
  SolverStatusCode code = SolveFailedToConverge;
  if (rtn==Belos::Converged) code = SolveConverged;
  SolverState<double> state(code, "Belos solver completed", numIters, resid);
  
  return state;
}



