
#include "RTOpPack_ROpNormInf.hpp"
#include "opsUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNormInf, unitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, v);
  RTOpPack::ROpNormInf<Scalar> normInfOp;
  RCP<RTOpPack::ReductTarget> normInf = normInfOp.reduct_obj_create();
  normInfOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    normInf.ptr()
    );
  TEST_FLOATING_EQUALITY( normInfOp(*normInf), ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNormInf, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNormInf, nonunitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(n, 3, v);
  RTOpPack::ROpNormInf<Scalar> normInfOp;
  RCP<RTOpPack::ReductTarget> normInf = normInfOp.reduct_obj_create();
  normInfOp.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    normInf.ptr()
    );
  TEST_FLOATING_EQUALITY( normInfOp(*normInf), ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNormInf, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNormInf, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::ReductTargetScalar;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;

  const ScalarMag three = as<ScalarMag>(3.0);
  const ScalarMag four = as<ScalarMag>(4.0);
  const ScalarMag two = as<ScalarMag>(2.0);

  RTOpPack::ROpNormInf<Scalar> normInfOp;

  RCP<RTOpPack::ReductTarget> reduct1 = normInfOp.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = normInfOp.reduct_obj_create();

  ReductTargetScalar<ScalarMag> &scalarReduct1 =
    dyn_cast<ReductTargetScalar<ScalarMag> >(*reduct1); 
  ReductTargetScalar<ScalarMag> &scalarReduct2 =
    dyn_cast<ReductTargetScalar<ScalarMag> >(*reduct2); 

  scalarReduct1.set(three);
  scalarReduct2.set(four);
  normInfOp.reduct_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(two);
  normInfOp.reduct_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_EQUALITY( normInfOp(*reduct2), four );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNormInf, reduct )


} // namespace
