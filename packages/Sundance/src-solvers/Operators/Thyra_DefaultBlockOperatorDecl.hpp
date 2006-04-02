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

#ifndef THYRA_DEFAULTBLOCKOPERATORDECL_HPP
#define THYRA_DEFAULTBLOCKOPERATORDECL_HPP


#include "Thyra_LinearOpBase.hpp"
#include "Thyra_BlockOperatorBase.hpp"
#include "Thyra_SetableBlockOperatorBase.hpp"
#include "Teuchos_Array.hpp"

#ifndef TRILINOS_6
#include "Thyra_ProductMultiVectorBase.hpp"
#endif


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
  class DefaultBlockOperator : public Thyra::LinearOpBase<RangeScalar, DomainScalar>,
                               public BlockOperatorBase<RangeScalar, DomainScalar>,
                               public SetableBlockOperatorBase<RangeScalar, DomainScalar>,
    
                               public virtual Teuchos::Describable
  {
  public:

    /** */
    DefaultBlockOperator(const Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<DomainScalar> >& domain,
                         const Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<RangeScalar> >& range);

    /** Returns a product vector view of the domain space */
    Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<DomainScalar> > domainProductSpace() const ;
    

    /** Returns a product vector view of the range space */
    Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<RangeScalar> > rangeProductSpace() const ;


   /** \brief Return a smart pointer for the range space for <tt>this</tt> operator.
   */
  Teuchos::RefCountPtr< const Thyra::VectorSpaceBase<RangeScalar> > range() const ;

  /** \brief Return a smart pointer for the domain space for <tt>this</tt> operator.
   */
  Teuchos::RefCountPtr< const Thyra::VectorSpaceBase<DomainScalar> > domain() const ;

    /** Returns the number of block rows */
    int numBlockRows() const ;


    /** Returns  the number of block columns */
    int numBlockCols() const ;

    /** Returns the (i,j)-th block */
    Teuchos::RefCountPtr<LinearOpBase<RangeScalar, DomainScalar> > 
    getBlock(int i, int j) const ;

    /** Sets the (i,j)-th block */
    void setBlock(int i, int j, const Teuchos::RefCountPtr<LinearOpBase<RangeScalar, DomainScalar> >& subBlock); 


    /** */
    void apply(
    const EConj                             conj
    ,const MultiVectorBase<DomainScalar>    &X
    ,MultiVectorBase<RangeScalar>           *Y
    ,const RangeScalar                      alpha = Teuchos::ScalarTraits<RangeScalar>::one()
    ,const RangeScalar                      beta  = Teuchos::ScalarTraits<RangeScalar>::zero()
    ) const ;


    /** */
    bool applySupports( const EConj conj ) const ;

    /** */
    bool applyTransposeSupports( const EConj conj ) const;

    /** */
    void applyTranspose(
    const EConj                            conj
    ,const MultiVectorBase<RangeScalar>    &X
    ,MultiVectorBase<DomainScalar>         *Y
    ,const DomainScalar                     alpha = ScalarTraits<DomainScalar>::one()
    ,const DomainScalar                     beta  = ScalarTraits<DomainScalar>::zero()
    ) const;


    /** Return a short description of the operator */
    string description() const ;

    /** */
    std::ostream& describe(
    std::ostream                         &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ,const std::string                   leadingIndent
    ,const std::string                   indentSpacer
    ) const;
    

  protected:
    /** */
    void applyToVector(
    const EConj                             conj
    ,const ProductVectorBase<DomainScalar>    &X
    ,ProductVectorBase<RangeScalar>           *Y
    ,const RangeScalar                      alpha
    ,const RangeScalar                      beta
    ) const ;
    /** */
    void applyTransposeToVector(
    const EConj                             conj
    ,const ProductVectorBase<DomainScalar>    &X
    ,ProductVectorBase<RangeScalar>           *Y
    ,const RangeScalar                      alpha
    ,const RangeScalar                      beta
    ) const ;

    /** */
    void applyToMultiVector(
    const EConj                             conj
    ,const ProductMultiVectorBase<DomainScalar>    &X
    ,ProductMultiVectorBase<RangeScalar>           *Y
    ,const RangeScalar                      alpha
    ,const RangeScalar                      beta
    ) const ;
    /** */
    void applyTransposeToMultiVector(
    const EConj                             conj
    ,const ProductMultiVectorBase<DomainScalar>    &X
    ,ProductMultiVectorBase<RangeScalar>           *Y
    ,const RangeScalar                      alpha
    ,const RangeScalar                      beta
    ) const ;

    

    
  private:
    Teuchos::Array<Teuchos::Array<Teuchos::RefCountPtr<LinearOpBase<RangeScalar, DomainScalar> > > > blocks_;

    Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<DomainScalar> > domain_;

    Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<DomainScalar> > range_;
    
  }; 
}

#endif
