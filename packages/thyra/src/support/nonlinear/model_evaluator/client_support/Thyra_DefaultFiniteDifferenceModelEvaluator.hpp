// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_HPP
#define THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_DirectionalFiniteDiffCalculator.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_Time.hpp"

namespace Thyra {

/** \brief This class wraps any ModelEvaluator object and computes certain
 * derivatives using finite differences.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class DefaultFiniteDifferenceModelEvaluator
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief Utility object that computes directional finite differences */
  STANDARD_COMPOSITION_MEMBERS( Thyra::DirectionalFiniteDiffCalculator<Scalar>, direcFiniteDiffCalculator );

  /** \brief . */
  DefaultFiniteDifferenceModelEvaluator();

  /** \brief . */
  DefaultFiniteDifferenceModelEvaluator(
    const Teuchos::RCP<ModelEvaluator<Scalar> >                          &thyraModel
    ,const Teuchos::RCP<Thyra::DirectionalFiniteDiffCalculator<Scalar> > &direcFiniteDiffCalculator
    );

  /** \brief . */
  void initialize(
    const Teuchos::RCP<ModelEvaluator<Scalar> >                          &thyraModel
    ,const Teuchos::RCP<Thyra::DirectionalFiniteDiffCalculator<Scalar> > &direcFiniteDiffCalculator
    );

  /** \brief . */
  void uninitialize(
    Teuchos::RCP<ModelEvaluator<Scalar> >                          *thyraModel
    ,Teuchos::RCP<Thyra::DirectionalFiniteDiffCalculator<Scalar> > *direcFiniteDiffCalculator
    );

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{
  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgs() const;
  /** \brief . */
  void evalModel(
    const ModelEvaluatorBase::InArgs<Scalar>    &inArgs
    ,const ModelEvaluatorBase::OutArgs<Scalar>  &outArgs
    ) const;

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}
  
};

// /////////////////////////////////
// Implementations

// Constructors/initializers/accessors/utilities

template<class Scalar>
DefaultFiniteDifferenceModelEvaluator<Scalar>::DefaultFiniteDifferenceModelEvaluator()
{}

template<class Scalar>
DefaultFiniteDifferenceModelEvaluator<Scalar>::DefaultFiniteDifferenceModelEvaluator(
  const Teuchos::RCP<ModelEvaluator<Scalar> >                          &thyraModel
  ,const Teuchos::RCP<Thyra::DirectionalFiniteDiffCalculator<Scalar> > &direcFiniteDiffCalculator
  )
{
  initialize(thyraModel,direcFiniteDiffCalculator);
}

template<class Scalar>
void DefaultFiniteDifferenceModelEvaluator<Scalar>::initialize(
  const Teuchos::RCP<ModelEvaluator<Scalar> >                          &thyraModel
  ,const Teuchos::RCP<Thyra::DirectionalFiniteDiffCalculator<Scalar> > &direcFiniteDiffCalculator
  )
{
  this->ModelEvaluatorDelegatorBase<Scalar>::initialize(thyraModel);
  direcFiniteDiffCalculator_ = direcFiniteDiffCalculator;
}

template<class Scalar>
void DefaultFiniteDifferenceModelEvaluator<Scalar>::uninitialize(
  Teuchos::RCP<ModelEvaluator<Scalar> >                          *thyraModel
  ,Teuchos::RCP<Thyra::DirectionalFiniteDiffCalculator<Scalar> > *direcFiniteDiffCalculator
  )
{
  if(thyraModel) *thyraModel = this->getUnderlyingModel();
  this->ModelEvaluatorDelegatorBase<Scalar>::uninitialize();
  if(direcFiniteDiffCalculator) *direcFiniteDiffCalculator = direcFiniteDiffCalculator_;
  direcFiniteDiffCalculator_ = Teuchos::null;
}

// Overridden from ModelEvaulator.

template<class Scalar>
ModelEvaluatorBase::OutArgs<Scalar>
DefaultFiniteDifferenceModelEvaluator<Scalar>::createOutArgs() const
{
  typedef ModelEvaluatorBase MEB;
  const Teuchos::RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  const MEB::OutArgs<Scalar> wrappedOutArgs = thyraModel->createOutArgs();
  const int Np = wrappedOutArgs.Np(), Ng = wrappedOutArgs.Ng();
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(Np,Ng);
  outArgs.setSupports(wrappedOutArgs);
  // Just support derivatives of DgDp for now!
  for( int j = 0; j < Ng; ++j ) {
    for( int l = 0; l < Np; ++l ) {
      outArgs.setSupports(MEB::OUT_ARG_DgDp,j,l,MEB::DERIV_TRANS_MV_BY_ROW);
    }
  }
  // ToDo: Add support for more derivatives as needed!
  return outArgs;
}

template<class Scalar>
void DefaultFiniteDifferenceModelEvaluator<Scalar>::evalModel(
  const ModelEvaluatorBase::InArgs<Scalar>     &inArgs
  ,const ModelEvaluatorBase::OutArgs<Scalar>   &outArgs
  ) const
{
  typedef ModelEvaluatorBase MEB;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  typedef typename ST::magnitudeType ScalarMag;

  typedef RCP<VectorBase<Scalar> >         V_ptr;
  typedef RCP<const VectorBase<Scalar> >   CV_ptr;
  typedef RCP<MultiVectorBase<Scalar> >    MV_ptr;

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_BEGIN(
    "Thyra::DefaultFiniteDifferenceModelEvaluator",inArgs,outArgs
    );

  //
  // Just do the g_0(p_0) case for now!
  //

  const RCP<const VectorSpaceBase<Scalar> >
    p_space = thyraModel->get_p_space(0),
    g_space = thyraModel->get_g_space(0);

  //
  // Compute the base point
  //

  if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
    *out << "\nComputing the base point ...\n";

  const int Np = outArgs.Np();
  const int Ng = outArgs.Ng();
  MEB::InArgs<Scalar>  wrappedInArgs = inArgs;
  MEB::OutArgs<Scalar> baseFunc = thyraModel->createOutArgs();
  if( outArgs.supports(MEB::OUT_ARG_f) && outArgs.get_f().get() )
    baseFunc.set_f(outArgs.get_f());
  for( int j = 0; j < Ng; ++j ) {
    V_ptr g_j;
    if( (g_j=outArgs.get_g(j)).get() )
      baseFunc.set_g(j,g_j);
  }
  thyraModel->evalModel(wrappedInArgs,baseFunc);

  bool failed = baseFunc.isFailed();
  
  if(!failed) {
    //
    // Compute the derivatives
    //
    MEB::OutArgs<Scalar> deriv = thyraModel->createOutArgs();
    for( int l = 0; l < Np; ++l ) {
      if( outArgs.supports(MEB::OUT_ARG_DfDp,l).none()==false
          && outArgs.get_DfDp(l).isEmpty()==false )
      {
        deriv.set_DfDp(l,outArgs.get_DfDp(l));
      }
      for( int j = 0; j < Ng; ++j ) {
        if(
          outArgs.supports(MEB::OUT_ARG_DgDp,j,l).none()==false
          &&
          outArgs.get_DgDp(j,l).isEmpty()==false
          )
        {
          deriv.set_DgDp(j,l,outArgs.get_DgDp(j,l));
        }
      }
    }
    direcFiniteDiffCalculator_->calcDerivatives(
      *thyraModel,inArgs,baseFunc,deriv
      );
  }

  if(failed) {
    if(out.get() && static_cast<int>(verbLevel) >= static_cast<int>(Teuchos::VERB_LOW))
      *out
        << "\nEvaluation failed, returning NaNs ...\n";
    outArgs.setFailed();
  }

  THYRA_MODEL_EVALUATOR_DECORATOR_EVAL_MODEL_END();
  
}

// Public functions overridden from Teuchos::Describable

template<class Scalar>
std::string DefaultFiniteDifferenceModelEvaluator<Scalar>::description() const
{
  const Teuchos::RCP<const ModelEvaluator<Scalar> >
    thyraModel = this->getUnderlyingModel();
  std::ostringstream oss;
  oss << "Thyra::DefaultFiniteDifferenceModelEvaluator{";
  oss << "thyraModel=";
  if(thyraModel.get())
    oss << "\'"<<thyraModel->description()<<"\'";
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}

} // namespace Thyra

#endif // THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_HPP
