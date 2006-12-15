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

#ifndef TSFVECTORSPACEIMPL_HPP
#define TSFVECTORSPACEIMPL_HPP


#include "Thyra_ProductVectorSpaceBase.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "TSFDescribable.hpp"

using namespace TSFExtended;
using namespace Teuchos;
using std::ostream;

static inline Time& createVecTimer() 
{
  static RefCountPtr<Time> rtn 
    = TimeMonitor::getNewTimer("vector allocation"); 
  return *rtn;
}

 
//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::operator==(const VectorSpace<Scalar>& other) const 
{
  return isCompatible(other);  
}


//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::operator!=(const VectorSpace<Scalar>& other) const 
{
  return !(operator==(other));
}
    


//========================================================================
template <class Scalar>
Vector<Scalar> VectorSpace<Scalar>::createMember() const 
{
  TimeMonitor timer(createVecTimer());
  return Thyra::createMember(this->ptr());
}
    


//========================================================================
template <class Scalar>
int VectorSpace<Scalar>::lowestLocallyOwnedIndex() const
{
  const Thyra::SpmdVectorSpaceBase<Scalar>* mpiSpace 
    = dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>*>(this->ptr().get());
  if (mpiSpace != 0)
    {
      return mpiSpace->localOffset();
    }
  const Thyra::SpmdVectorSpaceBase<Scalar>* serialSpace 
    = dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>*>(this->ptr().get());
   if (serialSpace != 0)
     {
       return 0;
     }
   TEST_FOR_EXCEPTION(mpiSpace == 0 && serialSpace==0, runtime_error,
		      "don't know how to compute lowest local index for "
		      "a vector space that is neither MPI nor serial");
   return 0;
}

//========================================================================
template <class Scalar>
int VectorSpace<Scalar>::numLocalElements() const
{
  const Thyra::SpmdVectorSpaceBase<Scalar>* mpiSpace 
    = dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>*>(this->ptr().get());
  if (mpiSpace != 0)
    {
      return mpiSpace->localSubDim();
    }
   const Thyra::SpmdVectorSpaceBase<Scalar>* serialSpace 
    = dynamic_cast<const Thyra::SpmdVectorSpaceBase<Scalar>*>(this->ptr().get());
   if (serialSpace != 0)
     {
       return dim();
     }
   TEST_FOR_EXCEPTION(mpiSpace == 0 && serialSpace==0, runtime_error,
		      "don't know how to compute number of local elements for "
		      "a vector space that is neither MPI nor serial");
   return 0;
}
    



//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::isCompatible(const VectorSpace<Scalar>& vecSpc) const 
{
  TEST_FOR_EXCEPTION(vecSpc.ptr().get() == 0, runtime_error,
                     "null argument in VectorSpace<Scalar>::isCompatible()");
  return this->ptr().get()->isCompatible(*(vecSpc.ptr().get()));
}





//========================================================================
template <class Scalar>
bool VectorSpace<Scalar>::contains(const Vector<Scalar> &vec) const
{
  return (operator==(vec.space()));
}


//========================================================================
template <class Scalar>
int VectorSpace<Scalar>::numBlocks() const
{
  const Thyra::ProductVectorSpaceBase<Scalar>* pvs = 
    dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar>* > (this->ptr().get());
  if (pvs != 0)
    {
      return pvs->numBlocks();
    }
  return 1;
}



//========================================================================
template <class Scalar>
VectorSpace<Scalar> VectorSpace<Scalar>::getBlock(const int i) const
{
  const Thyra::ProductVectorSpaceBase<Scalar>* pvs = 
    dynamic_cast<const Thyra::ProductVectorSpaceBase<Scalar>* > (this->ptr().get());
  TEST_FOR_EXCEPTION(pvs == 0 && numBlocks()!=1, runtime_error,
		     "Space not a ProductVectorSpace" << endl);
  if (pvs != 0)
    {
      return pvs->getBlock(i);
    }
  return *this;
}


// //========================================================================
// template <class Scalar>
// void VectorSpace<Scalar>::setBlock(int i, 
// 				   const VectorSpace<Scalar>& space)
// {
//   const Thyra::ProductVectorSpace<Scalar>*  pvs = 
//     dynamic_cast<const Thyra::ProductVectorSpace<Scalar>* >  (this->ptr().get());

//   TEST_FOR_EXCEPTION(pvs == 0, runtime_error,
// 		     "Can't set block of vector space that is " <<
// 		     "not a ProductVectorSpace.");

//   Thyra::ProductVectorSpace<Scalar>* pvsc = const_cast<ProductVectorSpace<Scalar>*> (pvs);
//   pvsc->setBlock(i, space);
// }







#endif
