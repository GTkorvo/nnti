
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
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


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultMultiVectorLinearOpWithSolve, defaultConstruct,
  Scalar )
{
  const RCP<DefaultMultiVectorLinearOpWithSolve<Scalar> > dmvlows =
    multiVectorLinearOpWithSolve<Scalar>();
  TEST_ASSERT(nonnull(dmvlows));
  TEST_EQUALITY_CONST(dmvlows->range(), null);
  TEST_EQUALITY_CONST(dmvlows->domain(), null);
  out << "dmvlows = " << *dmvlows;
}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultMultiVectorLinearOpWithSolve,
  defaultConstruct )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultMultiVectorLinearOpWithSolve, basic,
  Scalar )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

  const Ordinal dim = 4;
  const int numBlocks = 3;

  const RCP<const VectorSpaceBase<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(dim);

  const RCP<MultiVectorBase<Scalar> > M = createMembers(vs, dim);
  randomize<Scalar>(-ST::one(), ST::one(), M.ptr());

  const RCP<const LinearOpWithSolveFactoryBase<Scalar> > lowsf = 
    defaultSerialDenseLinearOpWithSolveFactory<Scalar>();

  const RCP<LinearOpWithSolveBase<Scalar> > Minv = 
    linearOpWithSolve<Scalar>(*lowsf, M);
      
  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > dmvpvs =
    multiVectorProductVectorSpace<Scalar>(vs, numBlocks);

  const RCP<DefaultMultiVectorLinearOpWithSolve<Scalar> > dmvlows =
    multiVectorLinearOpWithSolve<Scalar>(Minv, dmvpvs, dmvpvs);

  TEST_ASSERT(nonnull(dmvlows));
  TEST_EQUALITY(dmvlows->range(), dmvpvs);
  TEST_EQUALITY(dmvlows->domain(), dmvpvs);
  out << "dmvlows = " << *dmvlows;

  Thyra::LinearOpTester<Scalar> linearOpTester;
  TEST_ASSERT(linearOpTester.check(*dmvlows, inOutArg(out)));

  Thyra::LinearOpWithSolveTester<Scalar> linearOpWithSolveTester;
  TEST_ASSERT(linearOpWithSolveTester.check(*dmvlows, &out));

}
THYRA_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultMultiVectorLinearOpWithSolve,
  basic )


} // namespace Thyra
