#include "SundanceDefs.hpp"
#include "SundanceFunctionIdentifier.hpp"
#include "SundanceExceptions.hpp"
#include "SundanceDeriv.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceSymbolicFuncElement.hpp"
#include "SundanceSymbolicFunc.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "SundanceOut.hpp"

using namespace SundanceCore;
using namespace std;

#define TEST_THROW(code, passFail) \
  TEUCHOS_TEST_THROW( code, std::exception, Out::os(), passFail)

#define TEST_NOTHROW(code, passFail) \
  TEUCHOS_TEST_NOTHROW( code, Out::os(), passFail)

Deriv funcDeriv(const Expr& u, const MultiIndex& mi = MultiIndex(0,0,0))
{
  if (u.size()==1U)
  {
    const SymbolicFuncElement* sfe 
      = dynamic_cast<const SymbolicFuncElement*>(u[0].ptr().get());
    TEST_FOR_EXCEPT(sfe==0);
    return SundanceCore::funcDeriv(sfe, mi);
  }
  else
  {
    const SymbolicFunc* sf 
      = dynamic_cast<const SymbolicFunc*>(u.ptr().get());
    TEST_FOR_EXCEPT(sf==0);
    return SundanceCore::funcDeriv(sf);
  }
}

Deriv normalDeriv(const Expr& u)
{
  const SymbolicFuncElement* sfe 
    = dynamic_cast<const SymbolicFuncElement*>(u[0].ptr().get());
  TEST_FOR_EXCEPT(sfe==0);
  return normalDeriv(sfe);
}

Deriv divergenceDeriv(const Expr& u)
{
  const SymbolicFunc* sf 
    = dynamic_cast<const SymbolicFunc*>(u.ptr().get());
  TEST_FOR_EXCEPT(sf==0);
  return divergenceDeriv(sf);
}

FunctionIdentifier fid(const Expr& u)
{
  const SymbolicFuncElement* sfe 
    = dynamic_cast<const SymbolicFuncElement*>(u[0].ptr().get());
  const SymbolicFunc* sf 
    = dynamic_cast<const SymbolicFunc*>(u.ptr().get());
  TEST_FOR_EXCEPT(sf==0 && sfe==0);
  if (sf) return sf->fid();
  else return sfe->fid();

}


bool testVecFunction()
{
  bool pass = true;
  Tabs tab;
  int verb=1;

  /* make a vector function */
  unsigned int dim = 3;
  Expr u = new UnknownFunctionStub("u", 1, 3);
  Expr phi = new UnknownFunctionStub("u", 0, 1);
  Expr v = new TestFunctionStub("v", 1, 3);

  TEUCHOS_TEST_EQUALITY(u.size(), dim, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(v.size(), dim, Out::os(), pass);
  
  /* Verify that D_mi u fails for nonzero multiindex.  */
//  MultiIndex m1(0,1,0);
//  TEST_THROW(
//    funcDeriv(u, m1), pass
//    );
  
  
  Deriv d_du = funcDeriv(u);
  Deriv d_dphi = funcDeriv(phi);
  Deriv d_dphi_x = d_dphi.derivWrtMultiIndex(MultiIndex(1,0,0));
  Deriv divU = divergenceDeriv(u);
  Deriv divV = divergenceDeriv(v);

  TEUCHOS_TEST_EQUALITY(divU.fid(), fid(u), Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking algSpec for divergence");
  TEUCHOS_TEST_EQUALITY(divU.algSpec().isScalar(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.algSpec().isVector(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.algSpec().isNormal(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.algSpec().isCoordinateComponent(), false, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking algSpec for vector operative func");
  TEUCHOS_TEST_EQUALITY(d_du.algSpec().isScalar(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.algSpec().isVector(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.algSpec().isNormal(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.algSpec().isCoordinateComponent(), false, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking algSpec for scalar operative func");
  TEUCHOS_TEST_EQUALITY(d_dphi.algSpec().isScalar(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.algSpec().isVector(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.algSpec().isNormal(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.algSpec().isCoordinateComponent(), false, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking algSpec for d/d(Dx(phi))");
  TEUCHOS_TEST_EQUALITY(d_dphi_x.algSpec().isScalar(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.algSpec().isVector(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.algSpec().isNormal(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.algSpec().isCoordinateComponent(), true, Out::os(), pass);


  SUNDANCE_BANNER1(verb, tab, "checking distinction between functional and spatial derivs");
  TEUCHOS_TEST_EQUALITY(divU.isFunctionalDeriv(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.isCoordDeriv(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.isFunctionalDeriv(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.isCoordDeriv(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.isFunctionalDeriv(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.isCoordDeriv(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.isFunctionalDeriv(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.isCoordDeriv(), false, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking identification of differentiation order");
  TEUCHOS_TEST_EQUALITY(divU.opOnFunc().derivOrder(), 1, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.opOnFunc().derivOrder(), 0, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.opOnFunc().derivOrder(), 0, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.opOnFunc().derivOrder(), 1, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking identification of test and unk funcs");
  TEUCHOS_TEST_EQUALITY(d_du.isTestFunction(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_du.isUnknownFunction(), true, Out::os(), pass);

  TEUCHOS_TEST_EQUALITY(d_dphi.isTestFunction(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi.isUnknownFunction(), true, Out::os(), pass);

  TEUCHOS_TEST_EQUALITY(d_dphi_x.isTestFunction(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(d_dphi_x.isUnknownFunction(), true, Out::os(), pass);

  TEUCHOS_TEST_EQUALITY(divU.isTestFunction(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.isUnknownFunction(), true, Out::os(), pass);

  TEUCHOS_TEST_EQUALITY(divV.isTestFunction(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divV.isUnknownFunction(), false, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking identification of spatial operators acting on operative functions");

  TEUCHOS_TEST_EQUALITY(divU.opOnFunc().isDivergence(), true, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.opOnFunc().isNormal(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.opOnFunc().isPartial(), false, Out::os(), pass);
  TEUCHOS_TEST_EQUALITY(divU.opOnFunc().isIdentity(), false, Out::os(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking that asking for component information from a divergence throws an exception");

  TEST_THROW(divU.opOnFunc().normalDerivOrder(), pass);
  TEST_THROW(divU.opOnFunc().mi(), pass);

  SUNDANCE_BANNER1(verb, tab, "checking that applying a partial derivative to a divergence throws an exception");
  TEST_THROW(divU.opOnFunc().derivWrtMultiIndex(MultiIndex(1,0,0)), pass);
  TEST_THROW(divU.opOnFunc().derivWrtMultiIndex(MultiIndex(0,1,0)), pass);
  TEST_THROW(divU.opOnFunc().derivWrtMultiIndex(MultiIndex(0,0,1)), pass);
  SUNDANCE_BANNER1(verb, tab, "checking that applying a zero-order partial derivative to a divergence does not throw an exception");
  TEST_NOTHROW(divU.opOnFunc().derivWrtMultiIndex(MultiIndex(0,0,0)), pass);

  SUNDANCE_BANNER1(verb, tab, "checking that applying a partial derivative to a divergence throws an exception");
  TEST_THROW(divU.derivWrtMultiIndex(MultiIndex(1,0,0)), pass);
  TEST_THROW(divU.derivWrtMultiIndex(MultiIndex(0,1,0)), pass);
  TEST_THROW(divU.derivWrtMultiIndex(MultiIndex(0,0,1)), pass);
  SUNDANCE_BANNER1(verb, tab, "checking that applying a zero-order partial derivative to a divergence does not throw an exception");
  TEST_NOTHROW(divU.derivWrtMultiIndex(MultiIndex(0,0,0)), pass);


  SUNDANCE_BANNER1(verb, tab, "all done!");

  return pass;
}

int main(int argc, char** argv)
{
  bool pass = true;
  try
    {
      pass = pass && testVecFunction();
    }
  catch(exception& e)
    {
      pass = false;
      Out::os() << "unexpected exception: " << e.what() << endl;
    }
  if (pass) 
  {
    cerr << "test PASSED" << endl;
    return 0;
  }
  else 
  {
    cerr << "test FAILED" << endl;
    return -1;
  }
}

