//@HEADER
// ***********************************************************************
//
//                           Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Rythmos_Types.hpp"
#include "Rythmos_UnitTestHelpers.hpp"

#include "Rythmos_ExplicitRKStepper.hpp"

#include "../SinCos/SinCosModel.hpp"

#include "Rythmos_DefaultIntegrator.hpp"
#include "Rythmos_InterpolationBuffer.hpp"
#include "Rythmos_SimpleIntegrationControlStrategy.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Thyra_DetachedVectorView.hpp"

namespace Rythmos {

using Thyra::VectorBase;
using Thyra::VectorSpaceBase;
using Teuchos::is_null;

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, create ) {
  RCP<SinCosModel> model = sinCosModel(false);
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>(model);
  TEST_EQUALITY_CONST( is_null(stepper), false );
}

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, setgetRKButcherTableau ) {
  RCP<SinCosModel> model = sinCosModel(false);
  RKButcherTableau<double> rkbt = createExplicit4Stage4thOrder_RKBT<double>();
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>(model,rkbt);
  TEST_EQUALITY_CONST( is_null(stepper), false );
  RKButcherTableau<double> rkbt_out = stepper->getRKButcherTableau();
  TEST_EQUALITY_CONST( rkbt == rkbt_out, true );
}

TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, getTimeRange ) {
  RCP<SinCosModel> model = sinCosModel(false);
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>(model);
  TimeRange<double> tr = stepper->getTimeRange();
  TEST_EQUALITY_CONST( tr.isValid(), true );
  TEST_EQUALITY_CONST( tr.lower(), 0.0 );
  TEST_EQUALITY_CONST( tr.upper(), 0.0 );
} 

// 12/17/08 tscoffe:  I need a model evaluator _without_ a nominal values to
// test the initialization behavior of the ERK stepper (and the ImplicitBDF stepper).

// Test the ERK stepper through the integrator
TEUCHOS_UNIT_TEST( Rythmos_ExplicitRKStepper, defaultIntegrator ) {
  // Integrator
  RCP<DefaultIntegrator<double> > integrator = defaultIntegrator<double>();

  // Stepper
  double finalTime = 1.0;
  RCP<SinCosModel> model = sinCosModel(false);
  DefaultRKButcherTableauFactory<double> rkbtFactory;
  RCP<ParameterList> rkbtPL = Teuchos::parameterList();
  rkbtPL->set("Method by name","Forward Euler");
  RKButcherTableau<double> rkbt = rkbtFactory.create(*rkbtPL);
  RCP<ExplicitRKStepper<double> > stepper = explicitRKStepper<double>(model,rkbt);
  integrator->setStepper(stepper, finalTime);

  // IntegrationControlStrategy to specify fixed steps for this stepper
  double dt = 0.1;
  RCP<ParameterList> pl = Teuchos::parameterList();
  pl->set("Take Variable Steps",false);
  pl->set("Fixed dt", dt);
  RCP<SimpleIntegrationControlStrategy<double> > intCont = simpleIntegrationControlStrategy<double>(pl);
  TimeRange<double> tr(0.0,finalTime);
  intCont->resetIntegrationControlStrategy(tr); // ??? Why do I need to do this?
  integrator->setIntegrationControlStrategy(intCont);

  // Ask integrator for points
  Array<double> time_vec;
  Array<RCP<const VectorBase<double> > > x_vec;
  double N = 10;
  for (int i=0 ; i<=N ; ++i) {
    double t = 0.0 + i*finalTime/N;
    time_vec.push_back(t);
  }
  integrator->getFwdPoints(time_vec,&x_vec,NULL,NULL);

  // Verify that these points are accurate
  // Since we're using Forward Euler, we can write down the exact numerical solution
  double tol = 1.0e-10;
  double exact_x0 = 0.0; // nominal values on SinCosModel
  double exact_x1 = 1.0;
  for (int i=0 ; i<=N ; ++i) {
    {
      Thyra::ConstDetachedVectorView<double> x_vec_view( *(x_vec[i]) );
      TEST_FLOATING_EQUALITY( exact_x0, x_vec_view[0], tol );
      TEST_FLOATING_EQUALITY( exact_x1, x_vec_view[1], tol );
    }
    double x0 = exact_x0;
    double x1 = exact_x1;
    exact_x0 += dt*x1;
    exact_x1 -= dt*x0;
  }
}

} // namespace Rythmos

