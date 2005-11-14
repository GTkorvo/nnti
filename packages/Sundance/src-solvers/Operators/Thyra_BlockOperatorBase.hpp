/* @HEADER@ */
/* ***********************************************************************
//
//           TSFExtended: Trilinos Solver Framework Extended
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
// **********************************************************************/
 /* @HEADER@ */

#ifndef THYRA_BLOCKOPERATORBASE_HPP
#define THYRA_BLOCKOPERATORBASE_HPP

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"





namespace Thyra
{
  /**
   * Class BlockOperatorBase provides an abstract interface for 
   * block operators.
   *
   * @author Paul T Boggs (ptboggs@sandia.gov)
   * @author Kevin Long (ptboggs@sandia.gov)
   * @author Ross Bartlett (rabartl@sandia.gov)
   */
  template <class RangeScalar, class DomainScalar>
  class BlockOperatorBase 
  {
  public:

    /** Returns a product vector view of the domain space */
    virtual Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<DomainScalar> > domainProductSpace() const = 0 ;
    

    /** Returns a product vector view of the range space */
    virtual Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<RangeScalar> > rangeProductSpace() const = 0 ;

    /** Returns the number of block rows */
    virtual int numBlockRows() const = 0 ;


    /** Returns  the number of block columns */
    virtual int numBlockCols() const = 0 ;

    /** Returns the (i,j)-th block */
    Teuchos::RefCountPtr<LinearOpBase<RangeScalar, DomainScalar> > 
    getBlock(int i, int j) const ;

  }; 
}

#endif
