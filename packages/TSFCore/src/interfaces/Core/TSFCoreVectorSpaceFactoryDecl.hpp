// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
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

// ///////////////////////////////////////////////////////////////
// TSFCoreVectorSpaceFactoryDecl.hpp

#ifndef TSFCORE_VECTOR_SPACE_FACTORY_DECL_HPP
#define TSFCORE_VECTOR_SPACE_FACTORY_DECL_HPP

#include "TSFCoreTypes.hpp"

namespace TSFCore {

///
/** Abstract interface for objects that can create vector spaces of a specified dimension.
 *
 * The primary role that a <tt>%VectorSpaceFactory</tt> object takes
 * is defined in the documentation for the class <tt>VectorSpace</tt>
 * and is related to the domain space of <tt>MultiVector</tt> objects.
 *
 * \ingroup TSFCore_fundamental_interfaces_code_grp
 */
template<class Scalar>
class VectorSpaceFactory {
public:

  ///
  virtual ~VectorSpaceFactory() {}

	/** @name Pure virtual functions that must be overridden */
	//@{

	///
	/** Create a vector space of the given dimension.
	 *
	 * @param  dim  [in] The dimension of the vector space to create.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>dim > 0</tt> (throw <tt>std::invalid_argument</tt>).
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get() != NULL</tt>
	 * <li> <tt>return->dim() == dim</tt>
	 * </ul>
	 *
	 * @return  Returns a smart reference counted pointer to a dynamically
	 * allocated vector space object that can be used to create vector.
	 */
	virtual Teuchos::RefCountPtr< const VectorSpace<Scalar> > createVecSpc(int dim) const = 0;

	//@}

}; // end class VectorSpaceFactory

} // end namespace TSFCore

#endif  // TSFCORE_VECTOR_SPACE_FACTORY_DECL_HPP
