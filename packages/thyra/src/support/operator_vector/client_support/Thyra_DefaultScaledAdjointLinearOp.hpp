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

#ifndef THYRA_SCALED_ADJOINT_LINEAR_OP_HPP
#define THYRA_SCALED_ADJOINT_LINEAR_OP_HPP

#include "Thyra_DefaultScaledAdjointLinearOpDecl.hpp"
#include "Thyra_ScaledAdjointLinearOpBase.hpp"

namespace Thyra {

//Constructors/initializers/accessors

template<class Scalar>
void DefaultScaledAdjointLinearOp<Scalar>::initialize(
  const Scalar                                               &scalar
  ,const ETransp                                             &transp
  ,const Teuchos::RefCountPtr<LinearOpBase<Scalar> >         &Op
  )
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
#ifdef _DEBUG
  TEST_FOR_EXCEPTION(
    Op.get()==NULL, std::invalid_argument
    ,"DefaultScaledAdjointLinearOp<"<<ST::name()<<">::initialize(scalar,transp,Op): Error!, Op.get()==NULL is not allowed!"
    );
#endif // _DEBUG
  Teuchos::RefCountPtr<DefaultScaledAdjointLinearOp<Scalar> >
    saOp = Teuchos::rcp_dynamic_cast<DefaultScaledAdjointLinearOp<Scalar> >(Op);
  if(saOp.get()) {
    origOp_            = saOp->origOp_;
    overallScalar_     = saOp->overallScalar_*scalar;
    overallTransp_     = trans_trans(saOp->overallTransp_,transp) ;
    my_index_          = saOp->my_index_ + 1;
    allScalarETransp_  = saOp->allScalarETransp_;
  }
  else {
    origOp_            = Op;
    overallScalar_     = scalar;
    overallTransp_     = transp;
    my_index_          = 0;
    allScalarETransp_  = Teuchos::rcp(new allScalarETransp_t());
  }
  allScalarETransp_->push_back(ScalarETransp<Scalar>(scalar,transp));
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultScaledAdjointLinearOp<Scalar>::getOp()
{
  assertInitialized();
  if( my_index_ > 0 ) {
    const ScalarETransp<Scalar> &scalar_transp = allScalarETransp_->at(my_index_);
    Teuchos::RefCountPtr<DefaultScaledAdjointLinearOp<Scalar> >
      Op = Teuchos::rcp(new DefaultScaledAdjointLinearOp<Scalar>());
    Op->origOp_            = origOp_;
    Op->overallScalar_     = overallScalar_/scalar_transp.scalar;
    Op->overallTransp_     = trans_trans(overallTransp_,scalar_transp.transp);
    Op->my_index_          = my_index_ - 1;
    Op->allScalarETransp_  = allScalarETransp_;
    return Op;
  }
  return origOp_;
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
DefaultScaledAdjointLinearOp<Scalar>::getOp() const
{
  return const_cast<DefaultScaledAdjointLinearOp<Scalar>*>(this)->getOp();
}

template<class Scalar>
void DefaultScaledAdjointLinearOp<Scalar>::uninitialize(
  Scalar                                              *scalar
  ,ETransp                                            *transp
  ,Teuchos::RefCountPtr<LinearOpBase<Scalar> >        *Op
  )
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;

  assertInitialized();

  const ScalarETransp<Scalar> &scalar_transp = allScalarETransp_->at(my_index_);
  if(scalar) *scalar = scalar_transp.scalar;
  if(transp) *transp = scalar_transp.transp;
  if(Op)     *Op     = getOp();

  origOp_           = Teuchos::null;
  overallScalar_    = ST::zero();
  overallTransp_    = NOTRANS;
  allScalarETransp_ = Teuchos::null;

}

// Overridden from Teuchos::Describable
                                                
template<class Scalar>
std::string DefaultScaledAdjointLinearOp<Scalar>::description() const
{
  assertInitialized();
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  std::ostringstream oss;
  oss << "DefaultScaledAdjointLinearOp<" << ST::name() << ">{"
      << overallScalar() << ","<<toString(overallTransp())<<","
      << origOp_->description() << "}";
  return oss.str();
}

template<class Scalar>
void DefaultScaledAdjointLinearOp<Scalar>::describe(
  Teuchos::FancyOStream                &out_arg
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::RefCountPtr;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  assertInitialized();
  RefCountPtr<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  switch(verbLevel) {
    case Teuchos::VERB_DEFAULT:
    case Teuchos::VERB_LOW:
      *out << this->description() << std::endl;
      break;
    case Teuchos::VERB_MEDIUM:
    case Teuchos::VERB_HIGH:
    case Teuchos::VERB_EXTREME:
    {
      *out
        << "type = \'DefaultScaledAdjointLinearOp<" << ST::name() << ">\', "
        << "rangeDim = " << this->range()->dim() << ", domainDim = " << this->domain()->dim() << std::endl;
      OSTab tab(out);
      *out
        <<  "overallScalar="<< overallScalar() << std::endl
        <<  "overallTransp="<<toString(overallTransp()) << std::endl
        <<  "Constituent transformations:\n";
      for( int i = 0; i <= my_index_; ++i ) {
        const ScalarETransp<Scalar> &scalar_transp = (*allScalarETransp_)[my_index_-i];
        OSTab tab(out,i+1);
        if(scalar_transp.scalar != ST::one() && scalar_transp.transp != NOTRANS)
          *out << "scalar="<<scalar_transp.scalar<<",transp="<<toString(scalar_transp.transp)<<std::endl;
        else if(scalar_transp.scalar != ST::one())
          *out << "scalar="<<scalar_transp.scalar<<std::endl;
        else if( scalar_transp.transp != NOTRANS )
          *out << "transp="<<toString(scalar_transp.transp)<<std::endl;
        else
          *out << "no-transformation\n";
      }
      tab.incrTab(my_index_+2);
      *out << "origOp =\n" << Teuchos::describe(*origOp_,verbLevel);
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Should never be called!
  }
}

// Overridden from OpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultScaledAdjointLinearOp<Scalar>::range() const
{
  assertInitialized();
  return ( real_trans(this->overallTransp())==NOTRANS ? this->getOrigOp()->range() : this->getOrigOp()->domain() );
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> >
DefaultScaledAdjointLinearOp<Scalar>::domain() const
{
  assertInitialized();
  return ( real_trans(this->overallTransp())==NOTRANS ? this->getOrigOp()->domain() : this->getOrigOp()->range() );
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
DefaultScaledAdjointLinearOp<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be
}

template<class Scalar>
bool DefaultScaledAdjointLinearOp<Scalar>::opSupported(ETransp M_trans) const
{
  assertInitialized();
  return Thyra::opSupported(*this->getOrigOp(),trans_trans(this->overallTransp(),M_trans));
}

// Overridden from LinearOpBase

template<class Scalar>
void DefaultScaledAdjointLinearOp<Scalar>::apply(
  const ETransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  assertInitialized();
  Thyra::apply(*this->getOrigOp(),trans_trans(M_trans,this->overallTransp()),X,Y,Scalar(this->overallScalar()*alpha),beta);
}

// Overridden from ScaledAdointLinearOpBase

template<class Scalar>
Scalar DefaultScaledAdjointLinearOp<Scalar>::overallScalar() const
{
  return overallScalar_;
}

template<class Scalar>
ETransp DefaultScaledAdjointLinearOp<Scalar>::overallTransp() const
{
  return overallTransp_;
}

template<class Scalar>
Teuchos::RefCountPtr<LinearOpBase<Scalar> >
DefaultScaledAdjointLinearOp<Scalar>::getOrigOp()
{
  return origOp_;
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
DefaultScaledAdjointLinearOp<Scalar>::getOrigOp() const
{
  return origOp_;
}

} // namespace Thyra

#endif	// THYRA_SCALED_ADJOINT_LINEAR_OP_HPP
