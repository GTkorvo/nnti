
#include "Thyra_DelayedLinearOpWithSolve.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_LinearOpWithSolveTester.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "Thyra_UnitTestHelpers.hpp"


namespace Thyra {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::inOutArg;


//
// Unit Tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DelayedLinearOpWithSolve, defaultConstruct,
  Scalar )
{
  const RCP<DelayedLinearOpWithSolve<Scalar> > dlows =
    delayedLinearOpWithSolve<Scalar>();
  TEST_ASSERT(nonnull(dlows));
  TEST_EQUALITY_CONST(dlows->range(), null);
  TEST_EQUALITY_CONST(dlows->domain(), null);
  out << "dlows = " << *dlows;
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DelayedLinearOpWithSolve,
  defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DelayedLinearOpWithSolve, basic,
  Scalar )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

  const Ordinal dim = 4;

  const RCP<const VectorSpaceBase<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(dim);

  const RCP<MultiVectorBase<Scalar> > M = createMembers(vs, dim);
  randomize<Scalar>(-ST::one(), ST::one(), M.ptr());

  const RCP<DelayedLinearOpWithSolve<Scalar> > dlows =
    delayedLinearOpWithSolve<Scalar>(
      defaultLinearOpSource<Scalar>(M),
      defaultSerialDenseLinearOpWithSolveFactory<Scalar>()
      );

  TEST_ASSERT(nonnull(dlows));
  TEST_ASSERT(dlows->range()->isCompatible(*vs));
  TEST_ASSERT(dlows->domain()->isCompatible(*vs));
  out << "dlows = " << *dlows;

  Thyra::LinearOpTester<Scalar> linearOpTester;
  TEST_ASSERT(linearOpTester.check(*dlows, inOutArg(out)));

  Thyra::LinearOpWithSolveTester<Scalar> linearOpWithSolveTester;
  TEST_ASSERT(linearOpWithSolveTester.check(*dlows, &out));

}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DelayedLinearOpWithSolve,
  basic )


} // namespace Thyra
