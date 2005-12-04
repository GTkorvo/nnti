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

#ifndef THYRA_DEFAULTBLOCKOPERATORIMPL_HPP
#define THYRA_DEFAULTBLOCKOPERATORIMPL_HPP


#include "Thyra_DefaultBlockOperatorDecl.hpp"
#include "Thyra_SetableBlockOperatorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Thyra
{

  template <class RangeScalar, class DomainScalar>
  inline DefaultBlockOperator<RangeScalar, DomainScalar>::DefaultBlockOperator(
                                                                               const Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<DomainScalar> >& domain,
                                                                               const Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<RangeScalar> >& range)
    : blocks_(range->numBlocks(), domain->numBlocks()),
      domain_(domain),
      range_(range)
  {
    /* we leave all blocks uninitialized, their pointers set to null. All
     * uninitialized blocks will be ignored in operations */
  }


  template <class RangeScalar, class DomainScalar>
  inline Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<DomainScalar> > 
  DefaultBlockOperator<RangeScalar, DomainScalar>
  ::domainProductSpace() const 
  {
    return domain_;
  }

  template <class RangeScalar, class DomainScalar>
  inline Teuchos::RefCountPtr<const Thyra::ProductVectorSpaceBase<RangeScalar> > 
  DefaultBlockOperator<RangeScalar, DomainScalar>
  ::rangeProductSpace() const 
  {
    return range_;
  }




  template <class RangeScalar, class DomainScalar>
  inline Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<DomainScalar> > 
  DefaultBlockOperator<RangeScalar, DomainScalar>
  ::domain() const 
  {
    return rcp_dynamic_cast<const Thyra::VectorSpaceBase<DomainScalar> >(domain_);
  }

  template <class RangeScalar, class DomainScalar>
  inline Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<RangeScalar> > 
  DefaultBlockOperator<RangeScalar, DomainScalar>
  ::range() const 
  {
    return rcp_dynamic_cast<const Thyra::VectorSpaceBase<RangeScalar> >(range_);
  }


  template <class RangeScalar, class DomainScalar>
  inline int DefaultBlockOperator<RangeScalar, DomainScalar>::numBlockRows() const
  {
    return range_->numBlocks();
  }


  template <class RangeScalar, class DomainScalar>
  inline int DefaultBlockOperator<RangeScalar, DomainScalar>::numBlockCols() const
  {
    return domain_->numBlocks();
  }


  template <class RangeScalar, class DomainScalar>
  inline Teuchos::RefCountPtr<LinearOpBase<RangeScalar, DomainScalar> > 
  DefaultBlockOperator<RangeScalar, DomainScalar>::getBlock(int i, int j) const 
  {
    return blocks_[i][j];
  }


  template <class RangeScalar, class DomainScalar>
  inline void DefaultBlockOperator<RangeScalar, DomainScalar>
  ::setBlock(int i, int j, 
             const Teuchos::RefCountPtr<LinearOpBase<RangeScalar, DomainScalar> >& subBlock)
  {
    blocks_[i][j] = subBlock;
  }

  template <class RangeScalar, class DomainScalar>
  inline void DefaultBlockOperator<RangeScalar, DomainScalar>
  ::applyToVector(
    const EConj                             conj
    ,const ProductVectorBase<DomainScalar>    &X
    ,ProductVectorBase<RangeScalar>           *Y
    ,const RangeScalar                      alpha
    ,const RangeScalar                      beta
    ) const 
  {
    for (int i=0; i<numBlockRows(); i++)
      {
        Teuchos::RefCountPtr<Thyra::VectorBase<RangeScalar> > tmp
          = createMember(range_->getBlock(i));
        Thyra::assign(tmp.get(), ScalarTraits<RangeScalar>::zero());
        bool rowHasNonzeroBlock = false;
        
        /* multiply the (i,j)-th block by the i-th row of the input vector,
         * and add into the result vector. Note that with beta=one() on
         * the inner calls, we are adding into the result vector */
        for (int j=0; j<numBlockCols(); j++)
          {
            if (blocks_[i][j].get() == 0) continue;
            rowHasNonzeroBlock = true;
            blocks_[i][j]->apply(conj, *(X.getBlock(j)),
                                 tmp.get(), ScalarTraits<RangeScalar>::one(),
                                 ScalarTraits<RangeScalar>::one());
          }
        if (rowHasNonzeroBlock)
          {
            Thyra::scale(alpha, tmp.get());
            if (beta != ScalarTraits<RangeScalar>::zero())
              {
                Thyra::update(beta, *(Y->getBlock(i)), tmp.get());
              }
            Teuchos::RefCountPtr<Thyra::VectorBase<RangeScalar> > outBlock 
              = Y->getBlock(i);
            Thyra::assign(outBlock.get(), *tmp);
          }
        else
          {
            Thyra::assign(Y->getBlock(i).get(), ScalarTraits<RangeScalar>::one());
          }
      }
  }

  template <class RangeScalar, class DomainScalar>
  inline void DefaultBlockOperator<RangeScalar, DomainScalar>
  ::applyTransposeToVector(
    const EConj                             conj
    ,const ProductVectorBase<DomainScalar>    &X
    ,ProductVectorBase<RangeScalar>           *Y
    ,const RangeScalar                      alpha
    ,const RangeScalar                      beta
    ) const
  {
    for (int i=0; i<numBlockCols(); i++)
      {
        Teuchos::RefCountPtr<Thyra::VectorBase<DomainScalar> > tmp
          = createMember(domain_->getBlock(i));
        Thyra::assign(tmp.get(), ScalarTraits<RangeScalar>::zero());
        bool colHasNonzeroBlock = false;
      
        for (int j=0; j<numBlockRows(); j++)
          {
            if (blocks_[j][i].get() == 0) continue;
            colHasNonzeroBlock = true;
            blocks_[j][i]->applyTranspose(conj, *(X.getBlock(j).get()),
                                          tmp.get(), 
                                          ScalarTraits<DomainScalar>::one(),
                                          ScalarTraits<DomainScalar>::one());
          }
        if (colHasNonzeroBlock)
          {
            Thyra::scale(alpha, tmp.get());
            if (beta != ScalarTraits<DomainScalar>::zero())
              {
                Thyra::update(beta, *(Y->getBlock(i)), tmp.get());
              }
            Teuchos::RefCountPtr<Thyra::MultiVectorBase<DomainScalar> > outBlock 
              = Y->getBlock(i);
            Thyra::assign(outBlock.get(), *tmp);
          }
        else
          {
            Thyra::assign(Y->getBlock(i).get(), ScalarTraits<DomainScalar>::zero());
          }
      }
  }


  template <class RangeScalar, class DomainScalar>
  inline void DefaultBlockOperator<RangeScalar, DomainScalar>
  ::applyToMultiVector(
    const EConj                             conj
    ,const ProductMultiVectorBase<DomainScalar>    &X
    ,ProductMultiVectorBase<RangeScalar>           *Y
    ,const RangeScalar                      alpha
    ,const RangeScalar                      beta
    ) const
  {
    Index numMVCols = X.domain()->dim();
    for (int i=0; i<numBlockRows(); i++)
      {
        bool rowHasNonzeroBlock = false;
        Teuchos::RefCountPtr<Thyra::MultiVectorBase<RangeScalar> > tmp
          = createMembers(range_->getBlock(i), numMVCols);
        Thyra::assign(tmp.get(), ScalarTraits<RangeScalar>::zero());
        
        /* multiply the (i,j)-th block by the i-th row of the input vector,
         * and add into the result vector. Note that with beta=one() on
         * the inner calls, we are adding into the result vector */
        for (int j=0; j<numBlockCols(); j++)
          {
            if (blocks_[i][j].get() == 0) continue;
            rowHasNonzeroBlock = true;
            blocks_[i][j]->apply(conj, *(X.getBlock(j).get()),
                                 tmp.get(), ScalarTraits<RangeScalar>::one(),
                                 ScalarTraits<RangeScalar>::one());
          }
        if (rowHasNonzeroBlock)
          {
            Thyra::scale(alpha, tmp.get());
            if (beta != ScalarTraits<RangeScalar>::zero())
              {
                Thyra::update(beta, *(Y->getBlock(i)), tmp.get());
              }
            Teuchos::RefCountPtr<Thyra::MultiVectorBase<RangeScalar> > outBlock 
              = Y->getBlock(i);
            Thyra::assign(outBlock.get(), *tmp);
          }
        else
          {
            Thyra::assign(Y->getBlock(i).get(), ScalarTraits<RangeScalar>::zero());
          }
      }
  }
  
  template <class RangeScalar, class DomainScalar>
  inline void DefaultBlockOperator<RangeScalar, DomainScalar>
  ::applyTransposeToMultiVector(
    const EConj                             conj
    ,const ProductMultiVectorBase<DomainScalar>    &X
    ,ProductMultiVectorBase<RangeScalar>           *Y
    ,const RangeScalar                      alpha
    ,const RangeScalar                      beta
    ) const
  {
    Index numMVCols = X.domain()->dim();

    for (int i=0; i<numBlockCols(); i++)
      {
        bool colHasNonzeroBlock = false;
        Teuchos::RefCountPtr<Thyra::MultiVectorBase<DomainScalar> > tmp
          = createMembers(domain_->getBlock(i), numMVCols);
        Thyra::assign(tmp.get(), ScalarTraits<RangeScalar>::zero());
      
        for (int j=0; j<numBlockRows(); j++)
          {
            if (blocks_[j][i].get() == 0) continue;
            colHasNonzeroBlock = true;
            blocks_[j][i]->applyTranspose(conj, *(X.getBlock(j).get()),
                                          tmp.get(), 
                                          ScalarTraits<DomainScalar>::one(),
                                          ScalarTraits<DomainScalar>::one());
          }
        if (colHasNonzeroBlock)
          {
            Thyra::scale(alpha, tmp.get());
            if (beta != ScalarTraits<DomainScalar>::zero())
              {
                Thyra::update(beta, *(Y->getBlock(i)), tmp.get());
              }
            Teuchos::RefCountPtr<Thyra::MultiVectorBase<DomainScalar> > outBlock 
              = Y->getBlock(i);
            Thyra::assign(outBlock.get(), *tmp);
          }
        else
          {
            Thyra::assign(Y->getBlock(i).get(), 
                          ScalarTraits<DomainScalar>::zero());
          }
      }
  }

  template <class RangeScalar, class DomainScalar>
  inline void DefaultBlockOperator<RangeScalar, DomainScalar>
  ::apply(
          const EConj                             conj
          ,const MultiVectorBase<DomainScalar>    &X
          ,MultiVectorBase<RangeScalar>           *Y
          ,const RangeScalar                      alpha 
          ,const RangeScalar                      beta 
          ) const 
  {
    const Thyra::ProductMultiVectorBase<DomainScalar>* xPMV = 
      dynamic_cast<const Thyra::ProductMultiVectorBase<DomainScalar>* > (&X);

    Thyra::ProductMultiVectorBase<RangeScalar>* yPMV = 
      dynamic_cast< Thyra::ProductMultiVectorBase<RangeScalar>* > (Y);

    const Thyra::ProductVectorBase<DomainScalar>* xPV = 
      dynamic_cast<const Thyra::ProductVectorBase<DomainScalar>* > (&X);

    Thyra::ProductVectorBase<RangeScalar>* yPV = 
      dynamic_cast< Thyra::ProductVectorBase<RangeScalar>* > (Y);

    if (xPMV != 0 && yPMV != 0)
      {
        applyToMultiVector(conj, *xPMV, yPMV, alpha, beta);
      }
    else if (xPV != 0 && yPV != 0)
      {
        applyToVector(conj, *xPV, yPV, alpha, beta);
      }
    else
      {
        TEST_FOR_EXCEPTION(true, runtime_error, 
                           "input to DefaultBlockOperator::apply() is neither "
                           "a ProductVector nor a ProductMultiVector.");
      }
  }


  template <class RangeScalar, class DomainScalar>
  inline void DefaultBlockOperator<RangeScalar, DomainScalar>
  ::applyTranspose(
                   const EConj                            conj
                   ,const MultiVectorBase<RangeScalar>    &X
                   ,MultiVectorBase<DomainScalar>         *Y
                   ,const DomainScalar                     alpha
                   ,const DomainScalar                     beta
                   ) const
  {
    const Thyra::ProductMultiVectorBase<DomainScalar>* xPMV = 
      dynamic_cast<const Thyra::ProductMultiVectorBase<DomainScalar>* > (&X);

    Thyra::ProductMultiVectorBase<RangeScalar>* yPMV = 
      dynamic_cast< Thyra::ProductMultiVectorBase<RangeScalar>* > (Y);

    const Thyra::ProductVectorBase<DomainScalar>* xPV = 
      dynamic_cast<const Thyra::ProductVectorBase<DomainScalar>* > (&X);

    Thyra::ProductVectorBase<RangeScalar>* yPV = 
      dynamic_cast< Thyra::ProductVectorBase<RangeScalar>* > (Y);

    if (xPMV != 0 && yPMV != 0)
      {
        applyTransposeToMultiVector(conj, *xPMV, yPMV, alpha, beta);
      }
    else if (xPV != 0 && yPV != 0)
      {
        applyTransposeToVector(conj, *xPV, yPV, alpha, beta);
      }
    else
      {
        TEST_FOR_EXCEPTION(true, runtime_error, 
                           "input to DefaultBlockOperator::applyTranspose() "
                           "is neither "
                           "a ProductVector nor a ProductMultiVector.");
      }
  }



  template <class RangeScalar, class DomainScalar>
  inline bool DefaultBlockOperator<RangeScalar, DomainScalar>
  ::applySupports( const EConj conj ) const
  {
    for (unsigned int i=0; i<blocks_.size(); i++)
      {
        for (unsigned int j=0; j<blocks_[i].size(); j++)
          {
            if (blocks_[i][j].get() != 0 && !blocks_[i][j]->applySupports(conj))
              return false;
          }
      }
    return true;
  }

  template <class RangeScalar, class DomainScalar>
  inline bool DefaultBlockOperator<RangeScalar, DomainScalar>
  ::applyTransposeSupports( const EConj conj ) const
  {
    for (unsigned int i=0; i<blocks_.size(); i++)
      {
        for (unsigned int j=0; j<blocks_[i].size(); j++)
          {
            if (blocks_[i][j].get() != 0 
                && !blocks_[i][j]->applyTransposeSupports(conj))
              return false;
          }
      }
    return true;
  }

  


  template <class RangeScalar, class DomainScalar>
  inline std::string
  DefaultBlockOperator<RangeScalar, DomainScalar>
  ::description() const
  {
    TeuchosOStringStream ss;
    ss << "DefaultBlockOperator[" << numBlockRows() << " by " 
       << numBlockCols() << "]";
    return ss.str();
  }

  template <class RangeScalar, class DomainScalar>
  inline std::ostream&
  DefaultBlockOperator<RangeScalar, DomainScalar>
  ::describe(std::ostream &   	 os,
             const EVerbosityLevel  	verbLevel,
             const std::string  	leadingIndent,
             const std::string  	indentSpacer
             ) const
  {
    if (verbLevel == VERB_DEFAULT)
      {
        os << leadingIndent<< description() << endl;
      }
    else
      {
        os << leadingIndent << "<DefaultBlockOperator nRows=\"" << numBlockRows()
            << "\" nCols=\"" << numBlockCols() << "\">" << std::endl;
        for (int i=0; i<numBlockRows(); i++)
          {
            os << leadingIndent << indentSpacer 
               << "<BlockRow i=\"" << i << "\">" << std::endl;
            for (int j=0; j<numBlockCols(); j++)
              {
                os << leadingIndent << indentSpacer 
                   << indentSpacer 
                   << "<BlockCol j=\"" << j << "\">" << std::endl;
                os << leadingIndent << indentSpacer 
                   << indentSpacer ;
                if (getBlock(i,j).get() != 0)
                  {
                    os << getBlock(i,i)->description() << std::endl;
                  }
                else
                  {
                    os << "zero block" << std::endl;
                  }
                os << leadingIndent << indentSpacer 
                   << indentSpacer 
                   << "</BlockCol>" << std::endl;
              }
            os << leadingIndent << indentSpacer 
               << "</BlockRow>" << std::endl;
          }
        os << leadingIndent << "</BlockOperator>" << std::endl;
      }
    return os;
  }


}

#endif
