#include "SundanceEquationSet.hpp"
#include "SundanceExpr.hpp"
#include "SundanceUnknownFunctionStub.hpp"
#include "SundanceUnknownFuncElement.hpp"
#include "SundanceUnknownParameter.hpp"
#include "SundanceTestFunctionStub.hpp"
#include "SundanceDiscreteFunctionStub.hpp"
#include "SundanceCoordExpr.hpp"
#include "SundanceZeroExpr.hpp"
#include "SundanceDerivative.hpp"
#include "SundanceQuadratureFamilyStub.hpp"
#include "SundanceCellFilterStub.hpp"
#include "SundanceIntegral.hpp"



using namespace SundanceUtils;
using namespace SundanceCore;
using namespace SundanceCore;
using namespace Teuchos;

using SundanceCore::List;

int main()
{
  Expr dx = new Derivative(0);
  Expr dy = new Derivative(1);

  Expr u = new UnknownFunctionStub("u");
  Expr alpha = new UnknownParameter("alpha");
  Expr alpha0 = new Parameter(3.14, "alpha0");
  Expr beta = new UnknownParameter("beta");
  Expr beta0 = new Parameter(2.72, "beta0");
  Expr v = new TestFunctionStub("v");

  Out::os() << "u=" << u << endl;
  Out::os() << "v=" << v << endl;
  Out::os() << "alpha=" << alpha << endl;

  Expr x = new CoordExpr(0);
  Expr y = new CoordExpr(1);

  Expr u0 = new DiscreteFunctionStub("u0");
  Expr zero = new ZeroExpr();

  RCP<CellFilterStub> cells = rcp(new CellFilterStub());
  RCP<QuadratureFamilyStub> quad = rcp(new QuadratureFamilyStub(1));

  WatchFlag watchMe("watch eqn");
  watchMe.setParam("symbolic preprocessing", 6);
  

  Expr w = Integral(cells, v*u*alpha, quad, watchMe);
  Expr dum;
  Array<Expr> dum2;

  EquationSet eqn(w, dum, tuple(v), tuple(u), tuple(u0), 
    dum, dum, alpha, alpha0,  dum2, dum2);

  Out::os() << "num unk params=" << eqn.numUnkParams() << endl;
  Out::os() << "num fixed params=" << eqn.numFixedParams() << endl;

  for (unsigned int r=0; r<eqn.numRegions(); r++)
  {
    const RegionQuadCombo& rqc = eqn.regionQuadCombos()[r];
    const DerivSet& derivs = eqn.nonzeroFunctionalDerivs(Sensitivities, rqc);

    for (DerivSet::const_iterator i=derivs.begin(); i!=derivs.end(); i++)
    {
      const MultipleDeriv& d = *i;
      Out::os() << "d=" << d << endl;
      for (MultipleDeriv::const_iterator j=d.begin(); j!=d.end(); j++)
      {
        Out::os() << "j=" << *j;
        if (j->isParameter()) Out::os() << " (parameter) ";
        Out::os() << endl;
      }
    }
  }
}


