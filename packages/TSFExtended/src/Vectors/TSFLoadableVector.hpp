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

#ifndef TSFLOADABLEVECTOR_HPP
#define TSFLOADABLEVECTOR_HPP

#include "TSFConfigDefs.hpp"
#include "TSFCoreTypes.hpp"

#ifndef DOXYGEN_DEVELOPER_ONLY

namespace TSFExtended
{
  using TSFCore::Index;
  /**
   * LoadableVector defines an interface through which elements can 
   * be loaded into a vector. Element loading is used extensively
   * by application codes in creating vectors, 
   * but should never be used by high-performance solver codes; this 
   * capability is therefore in TSFExtended rather than TSFCore.
   *
   * A TSFExtended vector type that will be
   * used in a context where loading is required should multiply inherit
   * from both TSFCore::Vector and TSFExtended::LoadableVector.
   * 
   * Elements can by loaded one at a time
   * or in batches. The methods to load single elements arew pure virtual
   * and thus must be defined by derived classes. 
   * Loading in batches will usually be more efficient
   * provided the underlying vector implementation supports it. 
   * For those types not supporting batch loading, LoadableVector provides
   * default batch loading functions which delegate to single-element loading.
   *
   * Elements can by loaded either by setting a value, or adding to an 
   * existing value. The latter will typically by used in finite-element
   * codes.
   *
   * @author Kevin Long (krlong@sandia.gov)
   */
  template <class Scalar>
  class LoadableVector 
    {
    public:
      /** virtual dtor */
      virtual ~LoadableVector() {;}

      /** set a single element at the given global index */
      virtual void setElement(Index globalIndex, const Scalar& value) = 0 ;

      /** add to the existing value of 
       * a single element at the given global index */
      virtual void addToElement(Index globalIndex, const Scalar& value) = 0 ;

      /** set a group of elements */
      virtual void setElements(size_t numElems, const Index* globalIndices, 
                               const Scalar* values) ;

      /** add to a group of elements */
      virtual void addToElements(size_t numElems, const Index* globalIndices, 
                         const Scalar* values);

      /** Do whatever finalization steps are needed by the implementation,
       for instance, synchronizing border elements. The default implementation
      * is a no-op. */
      virtual void finalizeAssembly() {;}
    };

  /* Default implementation of setElements makes multiple calls to
   * setElement(). If at all possible, this should be overridden
   * with a method specialized to the underlying type.  */
  template <class Scalar> 
  inline void LoadableVector<Scalar>::setElements(size_t numElems, 
                                                  const Index* globalIndices, 
                                                  const Scalar* values)
  {
    for (int i=0; i<numElems; i++)
      {
        setElement(globalIndices[i], values[i]);
      }
  }

  /* Default implementation of addToElements makes multiple calls to
   * addToElement(). If at all possible, this should be overridden
   * with a method specialized to the underlying type.  */
  template <class Scalar> 
  inline void LoadableVector<Scalar>::addToElements(size_t numElems, 
                                                    const Index* globalIndices, 
                                                    const Scalar* values)
  {
    for (int i=0; i<numElems; i++)
      {
        addToElement(globalIndices[i], values[i]);
      }
  }

  
  
}

#endif /* DOXYGEN_DEVELOPER_ONLY */
#endif
