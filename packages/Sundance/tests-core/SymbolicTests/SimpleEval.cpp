#include "SundanceExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionBase.hpp"
#include "SundanceTestFunctionBase.hpp"
#include "SundanceDiscreteFunctionBase.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceParameter.hpp"
#include "SundanceOut.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "SundanceDerivSet.hpp"
#include "SundanceEvalRegion.hpp"
#include "SundanceEvalManager.hpp"
#include "SundanceBruteForceEvaluator.hpp"
#include "SundanceEvalVectorArray.hpp"
#include "SundanceSymbPreprocessor.hpp"


using namespace SundanceUtils;
using namespace SundanceCore;
using namespace SundanceCore::Internal;
using namespace Teuchos;

static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}
static Time& ioTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("result output"); 
  return *rtn;
}
static Time& doitTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("doit"); 
  return *rtn;
}



void doit(const Expr& e, 
          const Expr& tests,
          const Expr& unks,
          const Expr& u0, 
          const EvalRegion& region)
{
  TimeMonitor t0(doitTimer());
  EvalManager mgr;
  mgr.setRegion(region);

  RefCountPtr<EvaluatorFactory> factory 
    = rcp(new BruteForceEvaluatorFactory());

  const EvaluatableExpr* ev 
    = dynamic_cast<const EvaluatableExpr*>(e[0].ptr().get());

  DerivSet d = SymbPreprocessor::setupExpr(e[0], 
                                           tests,
                                           unks,
                                           u0,
                                           region, factory.get());

  RefCountPtr<EvalVectorArray> results;

  ev->evaluate(mgr, results);

  results->print(cerr, d);
}

void testExpr(const Expr& e,  
              const Expr& tests,
              const Expr& unks,
              const Expr& u0, 
              const EvalRegion& region)
{
  cerr << endl 
       << "-------- testing " << e.toString() << " -------- " << endl;

  try
    {
      doit(e, tests, unks, u0, region);
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

int main(int argc, void** argv)
{
  
  try
		{
      MPISession::init(&argc, &argv);

      TimeMonitor t(totalTimer());
      SymbolicTransformation::verbosity() = 0;
      Evaluator::verbosity() = 0;
      EvalVector::verbosity() = 0;
      EvaluatableExpr::verbosity() = 0;
      Expr::showAllParens() = true;

      Expr dx = new Derivative(0);
      Expr dy = new Derivative(1);

			Expr u = new UnknownFunctionBase("u");
			Expr w = new UnknownFunctionBase("w");
			Expr v = new TestFunctionBase("v");
			Expr s = new TestFunctionBase("s");

      cerr << "u=" << u << endl;
      cerr << "v=" << v << endl;

      Expr x = new CoordExpr(0);
      Expr y = new CoordExpr(1);

      Expr u0 = new DiscreteFunctionBase("u0");
      Expr w0 = new DiscreteFunctionBase("w0");

      Array<Expr> tests;

      tests.append(v);

      tests.append(s + v);

      tests.append(u*v);

      tests.append(u*v + s);

      tests.append(s + u*v);

      tests.append(s + v+u);

      tests.append(v + v*u*u);

      tests.append(v*u*u + v);

      tests.append(v*w + v*u*u);

      tests.append(v*u*u + v*w);

      tests.append(v*(u+w));

      tests.append((u+w)*v);

      tests.append((v+s)*(u+w));

      tests.append(dx*v);

      tests.append(dx*v + dx*s);

      tests.append((dx*u)*(dx*v));

      tests.append(u*(dx*v));

      tests.append((dx*v)*u);

      tests.append((dx*u)*v);

      tests.append(v*(dx*u));

      tests.append(v*u*dx*u + v*w*dy*u);


      for (int i=0; i<tests.length(); i++)
        {
          testExpr(tests[i], 
                   SundanceCore::List(v, s),
                   SundanceCore::List(u, w),
                   SundanceCore::List(u0, w0),
                   EvalRegion(Teuchos::toString(i)));
        }

      

    }
	catch(exception& e)
		{
			Out::println(e.what());
		}
  TimeMonitor::summarize();

  MPISession::finalize();
}
