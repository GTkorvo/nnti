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
   * 
   *
   * @author Paul T Boggs (ptboggs@sandia.gov)
   * @author Kevin Long (ptboggs@sandia.gov)
   * @author Ross Bartlett (rabartl@sandia.gov)
   */
  template <class RangeScalar, class DomainScalar>
  class DefaultBlockOperator : public BlockOperatorBase<RangeScalar, DomainScalar>,
    public MutableBlockOperatorBase<RangeScalar, DomainScalar>
  {
  public:

    /** Returns a product vector view of the domain space */
    RefCountPtr<const Thyra::ProductVectorSpaceBase<DomainScalar> > domainProductSpace() const ;
    

    /** Returns a product vector view of the range space */
    RefCountPtr<const Thyra::ProductVectorSpaceBase<RangeScalar> > rangeProductSpace() const ;


    /** Returns the domain space */
    RefCountPtr<const Thyra::VectorSpaceBase<DomainScalar> > domain() const ;
    

    /** Returns the range space */
    RefCountPtr<const Thyra::VectorSpaceBase<RangeScalar> > range() const ;

    /** Returns the number of block rows */
    int numBlockRows() const ;


    /** Returns  the number of block columns */
    int numBlockCols() const ;

    /** Returns the (i,j)-th block */
    RefCountPtr<LinearOpBase<RangeScalar, DomainScalar> > 
    getBlock(int i, int j) const ;

    /** Sets the (i,j)-th block */
    void setBlock(int i, int j, const RefCountPtr<LinearOpBase<RangeScalar, DomainScalar> >& subBlock); 


    /** */
    void apply(
    const EConj                             conj
    ,const MultiVectorBase<DomainScalar>    &X
    ,MultiVectorBase<RangeScalar>           *Y
    ,const RangeScalar                      alpha = ScalarTraits<RangeScalar>::one()
    ,const RangeScalar                      beta  = ScalarTraits<RangeScalar>::zero()
    ) const ;


    /** */
    bool applySupports( const EConj conj ) const ;

    /** */
    bool applyTransposeSupports( const EConj conj ) const;

    /** */
    applyTranspose(
    const EConj                            conj
    ,const MultiVectorBase<RangeScalar>    &X
    ,MultiVectorBase<DomainScalar>         *Y
    ,const DomainScalar                     alpha = ScalarTraits<DomainScalar>::one()
    ,const DomainScalar                     beta  = ScalarTraits<DomainScalar>::zero()
    ) const;

    /** Do a block-by-block clone of this operator */
    RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> > clone() const;
    
  private:
    Array<Array<RefCountPtr<LinearOpBase<RangeScalar, DomainScalar> > > > blocks_;
    
  }; 
}

#endif
