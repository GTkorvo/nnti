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

#ifndef THYRA_LINEAR_OP_DEFAULT_BASE_HPP
#define THYRA_LINEAR_OP_DEFAULT_BASE_HPP

#include "Thyra_LinearOpDefaultBaseDecl.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_describeLinearOp.hpp"

namespace Thyra {


// Overridden from Teuchos::Describable


template<class RangeScalar, class DomainScalar>
std::string LinearOpDefaultBase<RangeScalar,DomainScalar>::description() const
{
  std::ostringstream oss;
  const Teuchos::RCP<const VectorSpaceBase<RangeScalar> >
    l_range = this->range();
  const Teuchos::RCP<const VectorSpaceBase<DomainScalar> >
    l_domain = this->domain();
  oss << Teuchos::Describable::description();
  if(!l_range.get()) {
    oss << "{range=NULL,domain=NULL}"; 
  }
  else {
    const Index dimDomain = l_domain->dim(), dimRange = l_range->dim();
    oss
      << "{rangeDim=" << dimRange
      << ",domainDim=" << dimDomain << "}";
  }
  return oss.str();
}


template<class RangeScalar, class DomainScalar>
void LinearOpDefaultBase<RangeScalar,DomainScalar>::describe(
  Teuchos::FancyOStream                &out
  ,const Teuchos::EVerbosityLevel      verbLevel
  ) const
{
  describeLinearOp(*this,out,verbLevel);
}


}	// end namespace Thyra

#endif // THYRA_LINEAR_OP_DEFAULT_BASE_HPP
