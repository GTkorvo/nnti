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

#ifndef THYRA_DEFAULT_ADDED_LINEAR_OP_HPP
#define THYRA_DEFAULT_ADDED_LINEAR_OP_HPP

#include "Thyra_DefaultAddedLinearOpDecl.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_Utils.hpp"

namespace Thyra {

// Constructors/initializers/accessors

template<class Scalar>
DefaultAddedLinearOp<Scalar>::DefaultAddedLinearOp()
{}

template<class Scalar>
DefaultAddedLinearOp<Scalar>::DefaultAddedLinearOp(
  const int                                                   numOps
  ,const Teuchos::RCP<LinearOpBase<Scalar> >          Ops[]
  )
{
  initialize(numOps,Ops);
}

template<class Scalar>
DefaultAddedLinearOp<Scalar>::DefaultAddedLinearOp(
  const int                                                   numOps
  ,const Teuchos::RCP<const LinearOpBase<Scalar> >    Ops[]
  )
{
  initialize(numOps,Ops);
}

template<class Scalar>
void DefaultAddedLinearOp<Scalar>::initialize(
  const int                                                   numOps
  ,const Teuchos::RCP<LinearOpBase<Scalar> >          Ops[]
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( numOps <= 0 || Ops == NULL );
#endif
  Ops_.resize(numOps);
  for( int k = 0; k < numOps; ++k )
    Ops_[k].initialize(Ops[k]);
  validateOps();
  setupDefaultObjectLabel();
}

template<class Scalar>
void DefaultAddedLinearOp<Scalar>::initialize(
  const int                                                   numOps
  ,const Teuchos::RCP<const LinearOpBase<Scalar> >    Ops[]
  )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( numOps <= 0 || Ops == NULL );
#endif
  Ops_.resize(numOps);
  for( int k = 0; k < numOps; ++k )
    Ops_[k].initialize(Ops[k]);
  validateOps();
  setupDefaultObjectLabel();
}

template<class Scalar>
void DefaultAddedLinearOp<Scalar>::uninitialize()
{
  Ops_.resize(0);
  setupDefaultObjectLabel();
}

// Overridden form AddedLinearOpBase

template<class Scalar>
int DefaultAddedLinearOp<Scalar>::numOps() const
{
  return Ops_.size();
}

template<class Scalar>
bool DefaultAddedLinearOp<Scalar>::opIsConst(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k].isConst();
}

template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
DefaultAddedLinearOp<Scalar>::getNonconstOp(const int k)
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k].getNonconstObj();
}

template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultAddedLinearOp<Scalar>::getOp(const int k) const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( !( 0 <= k && k < numOps() ) );
#endif
  return Ops_[k].getConstObj();
}

// Overridden from LinearOpBase

template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultAddedLinearOp<Scalar>::range() const
{
  assertInitialized();
  return getOp(0)->range();
}

template<class Scalar>
Teuchos::RCP< const VectorSpaceBase<Scalar> >
DefaultAddedLinearOp<Scalar>::domain() const
{
  assertInitialized();
  return getOp(numOps()-1)->domain();
}

template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
DefaultAddedLinearOp<Scalar>::clone() const
{
  return Teuchos::null; // Not supported yet but could be!
}

// Overridden from Teuchos::Describable
                                                
template<class Scalar>
std::string DefaultAddedLinearOp<Scalar>::description() const
{
  assertInitialized();
  std::ostringstream oss;
  oss << Teuchos::Describable::description() << "{numOps = "<<numOps()<<"}";
  return oss.str();
}

template<class Scalar>
void DefaultAddedLinearOp<Scalar>::describe(
  Teuchos::FancyOStream                &out_arg
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar>  ST;
  using Teuchos::RCP;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  assertInitialized();
  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  const int numOps = Ops_.size();
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
        << Teuchos::Describable::description() << "{"
        << "rangeDim=" << this->range()->dim()
        << ",domainDim="<< this->domain()->dim() << "}\n";
      OSTab tab(out);
      *out
        <<  "numOps="<< numOps << std::endl
        <<  "Constituent LinearOpBase objects for M = Op[0]*...*Op[numOps-1]:\n";
      tab.incrTab();
      for( int k = 0; k < numOps; ++k ) {
        *out << "Op["<<k<<"] = " << Teuchos::describe(*getOp(k),verbLevel);
      }
      break;
    }
    default:
      TEST_FOR_EXCEPT(true); // Should never get here!
  }
}

// protected

// Overridden from SingleScalarLinearOpBase

template<class Scalar>
bool DefaultAddedLinearOp<Scalar>::opSupported(ETransp M_trans) const
{
  bool opSupported = true;
  for( int k = 0; k < static_cast<int>(Ops_.size()); ++k )
    if(!Thyra::opSupported(*getOp(k),M_trans)) opSupported = false;
  return opSupported;
  // ToDo: Cache these?
}

template<class Scalar>
void DefaultAddedLinearOp<Scalar>::apply(
  const ETransp                     M_trans
  ,const MultiVectorBase<Scalar>    &X
  ,MultiVectorBase<Scalar>          *Y
  ,const Scalar                     alpha
  ,const Scalar                     beta
  ) const
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(
    "DefaultAddedLinearOp<Scalar>::apply(...)",*this,M_trans,X,Y
    );
#endif // TEUCHOS_DEBUG  
  //
  // Y = alpha * op(M) * X + beta*Y
  //
  // =>
  //
  // Y = beta*Y + sum(alpha*op(Op[j])*X),j=0...numOps-1)
  //
  const int numOps = Ops_.size();
  for( int j = 0; j < numOps; ++j )
    Thyra::apply(*getOp(j),M_trans,X,Y,alpha,j==0?beta:ST::one());
}


// private


template<class Scalar>
void DefaultAddedLinearOp<Scalar>::validateOps()
{
  typedef std::string s;
  using Teuchos::toString;
#ifdef TEUCHOS_DEBUG
  try {
    const int numOps = Ops_.size();
    for( int k = 0; k < numOps; ++k ) {
      TEST_FOR_EXCEPT( Ops_[k]().get() == NULL );
      if( k > 0 ) {
        THYRA_ASSERT_LINEAR_OP_PLUS_LINEAR_OP_SPACES_NAMES(
          "DefaultMultipliedLinearOp<Scalar>::initialize(...)"
          ,*Ops_[0].getConstObj(),NOTRANS,("Ops[0]")
          ,*Ops_[k].getConstObj(),NOTRANS,("Ops["+toString(k)+"]")
          );
      }
    }
  }
  catch(...) {
    uninitialize();
    throw;
  }
#endif
}


template<class Scalar>
void DefaultAddedLinearOp<Scalar>::setupDefaultObjectLabel()
{
  std::ostringstream label;
  const int numOps = Ops_.size();
  for( int k = 0; k < numOps; ++k ) {
    std::string Op_k_label = Ops_[k].getConstObj()->getObjectLabel();
    if (Op_k_label.length() == 0)
      Op_k_label = "ANYM";
    if (k > 0)
      label << "+";
    label << "("<<Op_k_label<<")";
  }
  this->setObjectLabel(label.str());
  validateOps();
}


}	// end namespace Thyra


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
Thyra::nonconstAdd(
  const Teuchos::RCP<LinearOpBase<Scalar> >    &A
  ,const Teuchos::RCP<LinearOpBase<Scalar> >   &B
  )
{
  using Teuchos::arrayArg;
  using Teuchos::RCP;
  return Teuchos::rcp(
    new DefaultAddedLinearOp<Scalar>(
      2
      ,arrayArg<RCP<LinearOpBase<Scalar> > >(A,B)()
      )
    );
}

template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpBase<Scalar> >
Thyra::add(
  const Teuchos::RCP<const LinearOpBase<Scalar> >    &A
  ,const Teuchos::RCP<const LinearOpBase<Scalar> >   &B
  )
{
  using Teuchos::arrayArg;
  using Teuchos::RCP;
  return Teuchos::rcp(
    new DefaultAddedLinearOp<Scalar>(
      2
      ,arrayArg<RCP<const LinearOpBase<Scalar> > >(A,B)()
      )
    );
}

#endif	// THYRA_DEFAULT_ADDED_LINEAR_OP_HPP
