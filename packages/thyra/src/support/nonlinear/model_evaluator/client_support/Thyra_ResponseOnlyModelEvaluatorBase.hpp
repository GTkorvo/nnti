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

#ifndef THYRA_RESPONSE_ONLY_MODEL_EVALUATOR_BASE_HPP
#define THYRA_RESPONSE_ONLY_MODEL_EVALUATOR_BASE_HPP


#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_Assert.hpp"


namespace Thyra {

/** \brief This base class defines default function implementations
 * appropritate for a response-only model evaluator <tt>(p) -> g(j)</tt>,
 * for <tt>j=0...Ng-1</tt>.
 *
 * The minimum that a subclass must to is to define implementations for
 * <tt>Np()</tt>, <tt>Ng()</tt>, <tt>get_p_space()</tt>,
 * <tt>get_g_space()</tt>, <tt>createInArgs()</tt>, <tt>createOutArgs</tt>,
 * and <tt>evalModel()</tt>.
 */
template<class Scalar>
class ResponseOnlyModelEvaluatorBase : virtual public ModelEvaluator<Scalar> {
public:

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief Throws exception. */
  Teuchos::RCP<const VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief Returns null. */
  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \brief Throws exception. */
  Teuchos::RCP<const VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief Returns this->createInArgs(). */
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief Returns this->createInArgs(). */
  ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief Returns this->createInArgs(). */
  ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;
  /** \brief Thorws exception. */
  Teuchos::RCP<LinearOpWithSolveBase<Scalar> > create_W() const;
  /** \brief Returns Teuchos::null. */
  Teuchos::RCP<LinearOpBase<Scalar> > create_W_op() const;
  /** \brief Throws exception. */
  Teuchos::RCP<LinearOpBase<Scalar> > create_DfDp_op(int l) const;
  /** \brief Throws exception. */
  Teuchos::RCP<LinearOpBase<Scalar> > create_DgDx_op(int j) const;
  /** \brief Throws exception. */
  Teuchos::RCP<LinearOpBase<Scalar> > create_DgDp_op( int j, int l ) const;
  /** \brief Does nothing and ignores input. */
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
    const bool wasSolved
    );

  //@}
  
};


// /////////////////////////////////
// Implementations


// Public functions overridden from ModelEvaulator


template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::get_x_space() const
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error,
    "ModelEvaluator<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::get_x_space(): "
    "Error, this function was not overridden in *this = \'"<<this->description()<<"\'!"
    );
  return Teuchos::null;
}


template<class Scalar>
Teuchos::RCP<const Teuchos::Array<std::string> >
ResponseOnlyModelEvaluatorBase<Scalar>::get_p_names(int l) const
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE( l, 0, this->Np() );
#endif
  return Teuchos::null;
}

template<class Scalar>
Teuchos::RCP<const VectorSpaceBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::get_f_space() const
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error,
    "ModelEvaluator<"<<Teuchos::ScalarTraits<Scalar>::name()<<">::get_f_space(): "
    " Error, this function was not overridden in \'"
    <<this->description()<<"\'!"
    );
  return Teuchos::null;
}


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ResponseOnlyModelEvaluatorBase<Scalar>::getNominalValues() const
{ return this->createInArgs(); }


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ResponseOnlyModelEvaluatorBase<Scalar>::getLowerBounds() const
{ return this->createInArgs(); }


template<class Scalar>
ModelEvaluatorBase::InArgs<Scalar>
ResponseOnlyModelEvaluatorBase<Scalar>::getUpperBounds() const
{ return this->createInArgs(); }


template<class Scalar>
Teuchos::RCP<LinearOpWithSolveBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::create_W() const
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"Error, if \'W\' is supported by the ModelEvaluator subclass then"
    " this function create_W() must be overridden by the subclass to return"
    " a non-null object!"
    );
  return Teuchos::null; // Should never be called!
}


template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::create_W_op() const
{
  TEST_FOR_EXCEPTION(
    true, std::logic_error
    ,"Error, if \'W\' is supported by the ModelEvaluator subclass then"
    " this function create_W() must be overridden by the subclass "
    <<this->description()<<" to return a non-null object!"
    );
  return Teuchos::null; // Should never be called!
}


template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::create_DfDp_op(int l) const
{
  typedef ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = this->createOutArgs();
  TEST_FOR_EXCEPTION(
    outArgs.supports(MEB::OUT_ARG_DfDp,l).supports(MEB::DERIV_LINEAR_OP), std::logic_error
    ,"Error, The ModelEvaluator subclass "<<this->description()<<" says that it"
    " supports the LinearOpBase form of DfDp("<<l<<") (as determined from its OutArgs object created by createOutArgs())"
    " but this function create_DfDp_op(...) has not been overriden to create such an object!"
    );
  return Teuchos::null;
}


template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::create_DgDx_op(int j) const
{
  typedef ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = this->createOutArgs();
  TEST_FOR_EXCEPTION(
    outArgs.supports(MEB::OUT_ARG_DgDx,j).supports(MEB::DERIV_LINEAR_OP), std::logic_error
    ,"Error, The ModelEvaluator subclass "<<this->description()<<" says that it"
    " supports the LinearOpBase form of DgDx("<<j<<") (as determined from its OutArgs object created by createOutArgs())"
    " but this function create_DgDx_op(...) has not been overriden to create such an object!"
    );
  return Teuchos::null;
}


template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
ResponseOnlyModelEvaluatorBase<Scalar>::create_DgDp_op( int j, int l ) const
{
  typedef ModelEvaluatorBase MEB;
  MEB::OutArgs<Scalar> outArgs = this->createOutArgs();
  TEST_FOR_EXCEPTION(
    outArgs.supports(MEB::OUT_ARG_DgDp,j,l).supports(MEB::DERIV_LINEAR_OP), std::logic_error
    ,"Error, The ModelEvaluator subclass "<<this->description()<<" says that it"
    " supports the LinearOpBase form of DgDp("<<j<<","<<l<<") (as determined from its OutArgs object created by createOutArgs())"
    " but this function create_DgDp_op(...) has not been overriden to create such an object!"
    );
  return Teuchos::null;
}


template<class Scalar>
void ResponseOnlyModelEvaluatorBase<Scalar>::reportFinalPoint(
  const ModelEvaluatorBase::InArgs<Scalar>      &finalPoint
  ,const bool                                   wasSolved
  )
{
  // This final point is just ignored by default!
}


} // namespace Thyra


#endif // THYRA_RESPONSE_ONLY_MODEL_EVALUATOR_BASE_HPP
