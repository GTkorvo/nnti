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
// TSFCoreMPIVectorSpaceFactoryStdDecl.hpp

#ifndef TSFCORE_MPI_VECTOR_SPACE_FACTORY_STD_DECL_HPP
#define TSFCORE_MPI_VECTOR_SPACE_FACTORY_STD_DECL_HPP

#include "TSFCoreVectorSpaceFactory.hpp"

namespace TSFCore {

///
/** Implementation of a vector-space factory for a <tt>MPIVectorSpaceStd</tt> objects.
 *
 * This will create either serial (<tt>mpiComm==MPI_COMM_NULL</tt>) or
 * locally replicated (<tt>mpiComm!=MPI_COMM_NULL</tt>) vector space
 * objects (see <tt>createVecSpc()</tt>).  The primary motivation for
 * this subclass is to create locally replicated vector spaces for the
 * domain space of <tt>MPIMultiVectorStd</tt>.  In addition, an object
 * of this type is also returned from
 * <tt>MPIVectorSpaceBase::smallVecSpcFtcy()</tt>.
 *
 * Note that the default constructor is not allowed to avoid mistakes in using
 * this class.
 */
template<class Scalar>
class MPIVectorSpaceFactoryStd : public VectorSpaceFactory<Scalar> {
public:

	///
	/** Construct with an <tt>MPI_Comm</tt> object.
	 *
	 * @param  mpiComm      [in] The MPI communicator.  This object must be maintained
	 *                      by the client the entire time that <tt>this</tt> is in use.
	 * 
	 * Postconditions:<ul>
	 * <li><tt>this->mpiComm() == mpiComm</tt>
	 * </ul>
	 */
	MPIVectorSpaceFactoryStd( MPI_Comm  mpiComm );

	/// Return the MPI communicator.
	MPI_Comm mpiComm() const;

	/** @name Overridden from VectorSpaceFactory */
	//@{

	///
	/** Create a new locally-replicated <tt>MPIVectorSpaceStd</tt> object given its dimension!
	 *
	 * @param  dim  [in] The dimension of the (locally replicated) vector space to create.
	 *
	 * This function returns:
	 *
	 \code
	 return Teuchos::rcp(new MPIVectorSpaceStd(this->mpiComm(),dim,dim))</tt>
	 \endcode
	 *
	 * and therefore <tt>return->dim()==dim</tt> and this implementation
	 * fully satisfies the specification of
	 * <tt>VectorSpaceFactory::createVecSpc()</tt>.
	 */
 	Teuchos::RefCountPtr<const VectorSpace<Scalar> > createVecSpc(int dim) const;

	//@}

private:

	MPI_Comm  mpiComm_;

  MPIVectorSpaceFactoryStd(); // Not defined and not to be called!
  	
}; // end class MPIVectorSpaceFactoryStd

// ///////////////////////////
// Inline members

template<class Scalar>
inline
MPI_Comm MPIVectorSpaceFactoryStd<Scalar>::mpiComm() const
{
	return mpiComm_;
}

} // end namespace TSFCore

#endif  // TSFCORE_MPI_VECTOR_SPACE_FACTORY_STD_DECL_HPP
