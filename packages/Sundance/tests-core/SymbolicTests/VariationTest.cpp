#include "SundanceExpr.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceDiscreteFuncElement.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceProductTransformation.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceParameter.hpp"
#include "SundanceUnknownParameter.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceRegionQuadCombo.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceEvalVector.hpp"
#include "SundanceSymbPreprocessor.hpp"
#include "SundanceStringEvalMediator.hpp"

using namespace SundanceUtils;
using SundanceCore::List;
using namespace SundanceCore;
using namespace Teuchos;




static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}

static Time& doitTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("doit"); 
  return *rtn;
}



void doVariations(const Expr& e, 
  const Expr& vars,
  const Expr& varEvalPt,
  const Expr& unks,
  const Expr& unkEvalPt, 
  const Expr& unkParams,
  const Expr& unkParamEvalPts,
  const Expr& fixed,
  const Expr& fixedEvalPt, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts, 
  const EvalContext& region)
{
  TimeMonitor t0(doitTimer());
  EvalManager mgr;
  mgr.setRegion(region);

  static RefCountPtr<AbstractEvalMediator> mediator 
    = rcp(new StringEvalMediator());

  mgr.setMediator(mediator);

  const EvaluatableExpr* ev 
    = dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get());

  DerivSet d = SymbPreprocessor::setupVariations(e[0], 
    vars,
    varEvalPt,
    unks,
    unkEvalPt,
    unkParams,
    unkParamEvalPts,
    fixed,
    fixedEvalPt,
    fixedParams,
    fixedParamEvalPts,
    region,
    MatrixAndVector);

  Tabs tab;

  Array<double> constantResults;
  Array<RefCountPtr<EvalVector> > vectorResults;

  ev->evaluate(mgr, constantResults, vectorResults);

  ev->sparsitySuperset(region)->print(cerr, vectorResults, constantResults);
}



void doGradient(const Expr& e, 
  const Expr& vars,
  const Expr& varEvalPt,
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixed,
  const Expr& fixedEvalPts, 
  const EvalContext& region)
{
  TimeMonitor t0(doitTimer());
  EvalManager mgr;
  mgr.setRegion(region);

  static RefCountPtr<AbstractEvalMediator> mediator 
    = rcp(new StringEvalMediator());

  mgr.setMediator(mediator);

  const EvaluatableExpr* ev 
    = dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get());

  DerivSet d = SymbPreprocessor::setupGradient(e[0], 
    vars,
    varEvalPt,
    fixedParams,
    fixedParamEvalPts,
    fixed,
    fixedEvalPts,
    region,
    FunctionalAndGradient);

  Tabs tab;
  //  cerr << tab << *ev->sparsitySuperset(region) << endl;
  //  ev->showSparsity(cerr, region);

  // RefCountPtr<EvalVectorArray> results;

  Array<double> constantResults;
  Array<RefCountPtr<EvalVector> > vectorResults;

  ev->evaluate(mgr, constantResults, vectorResults);

  ev->sparsitySuperset(region)->print(cerr, vectorResults, constantResults);

  
  // results->print(cerr, ev->sparsitySuperset(region).get());
}



void doFunctional(const Expr& e, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixed,
  const Expr& fixedEvalPt, 
  const EvalContext& region)
{
  TimeMonitor t0(doitTimer());
  EvalManager mgr;
  mgr.setRegion(region);

  static RefCountPtr<AbstractEvalMediator> mediator 
    = rcp(new StringEvalMediator());

  mgr.setMediator(mediator);

  const EvaluatableExpr* ev 
    = dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get());

  DerivSet d = SymbPreprocessor::setupFunctional(e[0], 
    fixedParams,
    fixedParamEvalPts,
    fixed,
    fixedEvalPt,
    region,
    FunctionalOnly);

  Tabs tab;
  //  cerr << tab << *ev->sparsitySuperset(region) << endl;
  //  ev->showSparsity(cerr, region);

  // RefCountPtr<EvalVectorArray> results;

  Array<double> constantResults;
  Array<RefCountPtr<EvalVector> > vectorResults;

  ev->evaluate(mgr, constantResults, vectorResults);

  ev->sparsitySuperset(region)->print(cerr, vectorResults, constantResults);

  
  // results->print(cerr, ev->sparsitySuperset(region).get());
}



void testVariations(const Expr& e,  
  const Expr& vars,
  const Expr& varEvalPt,
  const Expr& unks,
  const Expr& unkEvalPt, 
  const Expr& unkParams,
  const Expr& unkParamEvalPts,
  const Expr& fixed,
  const Expr& fixedEvalPt, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,  
  const EvalContext& region)
{
  cerr << endl 
       << "------------------------------------------------------------- " << endl;
  cerr  << "-------- testing " << e.toString() << " -------- " << endl;
  cerr << endl 
       << "------------------------------------------------------------- " << endl;

  try
  {
    doVariations(e, vars, varEvalPt, 
      unks, unkEvalPt, 
      unkParams,
      unkParamEvalPts,
      fixed, fixedEvalPt, 
      fixedParams,
      fixedParamEvalPts,
      region);
  }
  catch(exception& ex)
  {
    cerr << "EXCEPTION DETECTED!" << endl;
    cerr << ex.what() << endl;
    // cerr << "repeating with increased verbosity..." << endl;
    //       cerr << "-------- testing " << e.toString() << " -------- " << endl;
    //       Evaluator::verbosity() = 2;
    //       EvalVector::verbosity() = 2;
    //       EvaluatableExpr::verbosity() = 2;
    //       Expr::showAllParens() = true;
    //       doit(e, region);
    exit(1);
  }
}

void testGradient(const Expr& e,  
  const Expr& vars,
  const Expr& varEvalPt,
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts,
  const Expr& fixed,
  const Expr& fixedEvalPt,  
  const EvalContext& region)
{
  cerr << endl 
       << "------------------------------------------------------------- " << endl;
  cerr  << "-------- testing " << e.toString() << " -------- " << endl;
  cerr << endl 
       << "------------------------------------------------------------- " << endl;

  try
  {
    doGradient(e, vars, varEvalPt,  
      fixedParams,
      fixedParamEvalPts,
      fixed, 
      fixedEvalPt, 
      region);
  }
  catch(exception& ex)
  {
    cerr << "EXCEPTION DETECTED!" << endl;
    cerr << ex.what() << endl;
    // cerr << "repeating with increased verbosity..." << endl;
    //       cerr << "-------- testing " << e.toString() << " -------- " << endl;
    //       Evaluator::verbosity() = 2;
    //       EvalVector::verbosity() = 2;
    //       EvaluatableExpr::verbosity() = 2;
    //       Expr::showAllParens() = true;
    //       doit(e, region);
    exit(1);
  }
}


void testFunctional(const Expr& e, 
  const Expr& fixedParams,
  const Expr& fixedParamEvalPts, 
  const Expr& fixed,
  const Expr& fixedEvalPt,  
  const EvalContext& region)
{
  cerr << endl 
       << "------------------------------------------------------------- " << endl;
  cerr  << "-------- testing " << e.toString() << " -------- " << endl;
  cerr << endl 
       << "------------------------------------------------------------- " << endl;

  try
  {
    doFunctional(e,   
      fixedParams,
      fixedParamEvalPts,fixed, fixedEvalPt, region);
  }
  catch(exception& ex)
  {
    cerr << "EXCEPTION DETECTED!" << endl;
    cerr << ex.what() << endl;
    // cerr << "repeating with increased verbosity..." << endl;
    //       cerr << "-------- testing " << e.toString() << " -------- " << endl;
    //       Evaluator::verbosity() = 2;
    //       EvalVector::verbosity() = 2;
    //       EvaluatableExpr::verbosity() = 2;
    //       Expr::showAllParens() = true;
    //       doit(e, region);
    exit(1);
  }
}


int main(int argc, char** argv)
{
  
  try
  {
    GlobalMPISession session(&argc, &argv);

    TimeMonitor t(totalTimer());

    int maxDiffOrder = 2;

    verbosity<SymbolicTransformation>() = VerbSilent;
    verbosity<Evaluator>() = VerbSilent;
    verbosity<EvalVector>() = VerbSilent;
    verbosity<EvaluatableExpr>() = VerbSilent;
    Expr::showAllParens() = true;
    ProductTransformation::optimizeFunctionDiffOps() = false;

    EvalVector::shadowOps() = true;

    Expr dx = new Derivative(0);
    Expr dy = new Derivative(1);
    Expr grad = List(dx, dy);

    Expr u = new UnknownFunctionStub("u");
    Expr lambda_u = new UnknownFunctionStub("lambda_u");
    Expr T = new UnknownFunctionStub("T");
    Expr lambda_T = new UnknownFunctionStub("lambda_T");
    Expr alpha = new UnknownParameter("alpha");

    Expr u0 = new DiscreteFunctionStub("u0");
    Expr lambda_u0 = new DiscreteFunctionStub("lambda_u0");
    Expr T0 = new DiscreteFunctionStub("T0");
    Expr lambda_T0 = new DiscreteFunctionStub("lambda_T0");
    Expr zero = new ZeroExpr();
    Expr alpha0 = new Parameter(3.14, "alpha0");

    Expr x = new CoordExpr(0);
    Expr y = new CoordExpr(1);

    Expr empty;

    Array<Expr> tests;

//#define BLAHBLAH 1
#ifdef BLAHBLAH
    verbosity<Evaluator>() = VerbExtreme;
    verbosity<SparsitySuperset>() = VerbExtreme;
    verbosity<EvaluatableExpr>() = VerbExtreme;
#endif

    tests.append( 0.5*(u-1.404)*(u-1.404)
      + 0.5*(grad*u)*(grad*u)
      + 0.5*alpha*alpha
      + (grad*lambda_u)*(grad*u)
      + lambda_u*alpha);


    cerr << endl << "============== u STATE EQUATIONS =================" << endl;

    for (int i=0; i<tests.length(); i++)
    {
      RegionQuadCombo rqc(rcp(new CellFilterStub()), 
        rcp(new QuadratureFamilyStub(1)));
      EvalContext context(rqc, makeSet(1,2), EvalContext::nextID());
      testVariations(tests[i], 
        List(lambda_u),
        List(zero),
        List(u),
        List(u0),
        empty,
        empty,
        empty,
        empty,
        alpha,
        alpha0,
        context);
    }

    cerr << endl << "=============== u ADJOINT EQUATIONS =================" << endl;
    for (int i=0; i<tests.length(); i++)
    {
      RegionQuadCombo rqc(rcp(new CellFilterStub()), 
        rcp(new QuadratureFamilyStub(1)));
      EvalContext context(rqc, makeSet(1,2), EvalContext::nextID());
      testVariations(tests[i], 
        List(u),
        List(u0),
        List(lambda_u),
        List(zero),
        empty,
        empty,
        empty,
        empty,
        alpha,
        alpha0,
        context);
    }


    cerr << endl << "================ REDUCED GRADIENT ====================" << endl;
    for (int i=0; i<tests.length(); i++)
    {
      RegionQuadCombo rqc(rcp(new CellFilterStub()), 
        rcp(new QuadratureFamilyStub(1)));
      EvalContext context(rqc, makeSet(0,1), EvalContext::nextID());
      testGradient(tests[i], 
        alpha, 
        alpha0,
        empty,
        empty,
        List(u, lambda_u),
        List(u0, zero),
        context);
    }

    cerr << endl << "=================== FUNCTIONAL ====================" << endl;
    for (int i=0; i<tests.length(); i++)
    {
      RegionQuadCombo rqc(rcp(new CellFilterStub()), 
        rcp(new QuadratureFamilyStub(1)));
      EvalContext context(rqc, makeSet(0), EvalContext::nextID());
      testFunctional(tests[i], 
        alpha,
        alpha0,
        List(u, lambda_u),
        List(u0, zero),
        context);
    }
    TimeMonitor::summarize();
  }
	catch(exception& e)
  {
    Out::println(e.what());
  }


  
}
