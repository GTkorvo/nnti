#include "SundanceExpr.hpp"
#include "SundanceStdMathOps.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceSymbolicTransformation.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceParameter.hpp"
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
#include "SundanceEvaluationTester.hpp"

using namespace SundanceUtils;
using namespace SundanceTesting;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace Teuchos;


static Time& totalTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("total"); 
  return *rtn;
}


#define TESTER(func, equiv)\
{\
Tabs tabs1;\
cerr << tabs1 << endl << tabs1\
<< "------------- Testing " #func " -----------"\
<< endl << tabs1 << endl;\
bool thisTestIsOK = true;\
aSave=a;\
  bSave=b;\
  cSave=c;\
  dSave=d;\
func;\
error = fabs(a - (equiv));\
cerr << tabs1 << "vec op value = " << a << " check=" << equiv \
<< " |value-check|=" << fabs(a-equiv) << endl;\
if (error > tol)\
{\
thisTestIsOK=false;\
cerr << "value computation FAILED" << endl;\
isOK = false;\
}\
if (!thisTestIsOK)\
{\
failures.append(#func);\
}\
}

int main(int argc, char** argv)
{
  int stat = 0;

  try
		{
      GlobalMPISession session(&argc, &argv);
      Tabs tabs;
      TimeMonitor timer(totalTimer());

      

      EvalVector::shadowOps() = true;
      bool isOK = true;
      Array<string> failures;
      

      TempStack s(1);

      typedef RefCountPtr<EvalVector> Vec;

      Vec A = s.popVector();
      Vec B = s.popVector();
      Vec C = s.popVector();
      Vec D = s.popVector();

      double& a = *(A->start());
      double& b = *(B->start());
      double& c = *(C->start());
      double& d = *(D->start());

      a = 1.1;
      b = 0.8;
      c = 0.5;
      d = 1.2;

      A->setString("A");
      B->setString("B");
      C->setString("C");
      D->setString("D");

      cerr << "A = " << *A << endl;

      cerr << "B = " << *B << endl;

      cerr << "C = " << *C << endl;

      cerr << "D = " << *D << endl;

      double aSave;
      double bSave;
      double cSave;
      double dSave;
      double error;
      double tol = 1.0e-13;

      
      
      TESTER((A->add_SV(1.23, B.get())), (aSave + 1.23*b));

      TESTER((A->add_SVV(4.56, B.get(), C.get())), (aSave + 4.56*b*c));

      TESTER((A->add_V(B.get())), (aSave + b));

      TESTER((A->add_S(7.89)), (aSave + 7.89));

      TESTER((A->add_VV(B.get(), C.get())), (aSave + b*c));

      TESTER((A->multiply_S_add_SV(1.23, 4.56, D.get())), (aSave*1.23 + 4.56*d));

      TESTER((A->multiply_S_add_S(0.0123, 0.0456)), (aSave*0.0123 + 0.0456));

      TESTER((A->multiply_V_add_VVV(B.get(), B.get(), C.get(), D.get())), (aSave*b + b*c*d));

      TESTER((A->multiply_V_add_SVV(B.get(), 1.23, C.get(), D.get())), (aSave*b + 1.23*c*d));

      TESTER((A->multiply_V_add_SV(B.get(), 1.23, C.get())), (aSave*b + 1.23*c));

      TESTER((A->multiply_VV(B.get(), C.get())), (aSave*b*c));

      TESTER((A->multiply_SV(4.56, C.get())), (aSave*4.56*c));

      TESTER((A->multiply_V(C.get())), (aSave*c));

      TESTER((A->multiply_S(1.2)), (aSave*1.2));

      
      TESTER((A->setTo_S_add_SVV(4.56, 1.23, C.get(), D.get())), (4.56 + 1.23*c*d));

      TESTER((A->setTo_S_add_VV(4.56, C.get(), D.get())), (4.56 + c*d));

      TESTER((A->setTo_S_add_SV(4.56, 1.23, D.get())), (4.56 + 1.23*d));

      TESTER((A->setTo_S_add_V(4.56, C.get())), (4.56 + c));

      TESTER((A->setTo_V(B.get())), (b));

      TESTER((A->setTo_VV(B.get(), C.get())), (b*c));

      TESTER((A->setTo_SV(1.23, C.get())), (1.23*c));

      TESTER((A->setTo_SVV(1.23, C.get(), D.get())), (1.23*c*d));


      if (isOK)
        {
          cerr << "all tests PASSED!" << endl;
        }
      else
        {
          stat = -1;
          cerr << "test FAILED!" << endl;
        }
      TimeMonitor::summarize();
    }
	catch(exception& e)
		{
      stat = -1;
      cerr << "test FAILED!" << endl;
      cerr << "detected exception: " << e.what() << endl;
		}

  return stat;
  
}
