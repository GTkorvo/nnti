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

#ifndef TSFVECTORTYPEEXTENSIONS_HPP
#define TSFVECTORTYPEEXTENSIONS_HPP

#include "TSFHandle.hpp"
#include "TSFCoreVectorSpaceFactory.hpp"
#include "TSFVectorSpaceDecl.hpp"
#include "TSFLinearOperatorDecl.hpp" //Decl added by ptb
#include "TSFGhostImporter.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /**
   * VectorTypeExtensions provides extensions to 
   * the TSFCore::VectorSpaceFactory
   * object, appropriate to the interface with applications codes.
   * TSFCore::VectorSpaceFactory has a method to create a small serial
   * vector space for use in multivector operators, however, that is 
   * insufficient for the needs of applications.
   *
   * <h4> Notes for subclass developers </h4>
   *
   * Subclasses should also derive from some subclass of
   * TSFCore::VectorSpaceFactory in order to use the createReplicatedSpace()
   * method. 
   *
   * Because applications may use the clean syntax
   * \code
   * VectorType vt = new MyVectorType();
   * \endcode
   * subclasses should also derive from
   * Handleable<VectorTypeBase<Scalar> and implement the getRcp() method.
   *
   * Subclasses might optionally implement the Printable and Describable
   * interfaces
   *
   */
  template <class Scalar>
  class VectorTypeExtensions 
  {
  public:

    /** create a distributed vector space.
     * @param dimension the dimension of the space 
     * @param nLocal number of indices owned by the local processor
     * @param locallyOwnedIndices array of indices owned by this processor  
     */
    virtual RefCountPtr<const TSFCore::VectorSpace<Scalar> >
    createSpace(int dimension, 
                int nLocal,
                const int* locallyOwnedIndices) const = 0 ;

    /**  
     * Create an importer for accessing ghost elements.
     * @param space the distributed vector space on which ghost elements
     * are to be shared
     * @param nGhost number of ghost elements needed by this processor
     * @param ghostIndices read-only C array of off-processor indices needed
     * by this processor.
     * @return A RCP to a GhostImporter object.
     */
    virtual RefCountPtr<GhostImporter<Scalar> > 
    createGhostImporter(const VectorSpace<Scalar>& space,
                        int nGhost,
                        const int* ghostIndices) const = 0 ;

    
    /**
     * Create an empty matrix of type compatible with this vector type,
     * sized according to the given domain and range spaces.
     */
    virtual LinearOperator<Scalar>
    createMatrix(const VectorSpace<Scalar>& domain,
                 const VectorSpace<Scalar>& range) const = 0 ;
  };
  
}

#endif
