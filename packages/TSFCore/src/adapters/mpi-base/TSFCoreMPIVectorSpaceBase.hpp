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

// ////////////////////////////////////////////////////////////////////////////
// TSFCoreMPIVectorSpaceBase.hpp

#ifndef TSFCORE_MPI_VECTOR_SPACE_BASE_HPP
#define TSFCORE_MPI_VECTOR_SPACE_BASE_HPP

#include "TSFCoreMPIVectorSpaceBaseDecl.hpp"
#ifdef RTOp_USE_MPI
#  include "Teuchos_RawMPITraits.hpp"
#endif

namespace TSFCore {

///
template<class Scalar>
MPIVectorSpaceBase<Scalar>::MPIVectorSpaceBase()
	:mapCode_(-1),isInCore_(false),defaultLocalOffset_(-1)
{}

// Virtual methods with default implementations

template<class Scalar>
Index MPIVectorSpaceBase<Scalar>::localOffset() const
{
	return defaultLocalOffset_;
}

template<class Scalar>
Index MPIVectorSpaceBase<Scalar>::mapCode() const
{
	return mapCode_;
}

// Overridden from VectorSpace

template<class Scalar>
bool MPIVectorSpaceBase<Scalar>::isInCore() const
{
	return isInCore_;
}

template<class Scalar>
bool MPIVectorSpaceBase<Scalar>::isCompatible( const VectorSpace<Scalar>& vecSpc ) const
{
	if( isInCore() && vecSpc.isInCore() )
		return this->dim() == vecSpc.dim();
	const MPIVectorSpaceBase<Scalar>
		*mpiVecSpc = dynamic_cast<const MPIVectorSpaceBase<Scalar>*>(&vecSpc);
	if(mpiVecSpc)
		return mapCode() == mpiVecSpc->mapCode();
	return false;
}

// protected

template<class Scalar>
void MPIVectorSpaceBase<Scalar>::updateState()
{
	const Index localSubDim = this->localSubDim(); 
	if( localSubDim > 0 ) {
		const MPI_Comm mpiComm = this->mpiComm();
		const Index globalDim = this->dim(); 
		int numProc = 1;
		int procRank = 0;
#ifdef RTOp_USE_MPI
		if( mpiComm != MPI_COMM_NULL ) {
			MPI_Comm_size( mpiComm, &numProc );
			MPI_Comm_rank( mpiComm, &procRank );
		}
		if( numProc > 1 && localSubDim < globalDim ) {
			//
			// Here we will make a map code out of just the local
			// sub-dimension on each processor.  If each processor
			// has the same number of local elements, then the maps
			// will be the same and this is all you need for
			// RTOp compatibility unless the operations are not
			// coordinate invariant.  I will work on this issue
			// if it becomes a problem.
			//
			Index localCode = localSubDim % (procRank+1) + localSubDim;
			MPI_Allreduce(
				&localCode                              // sendbuf
				,&mapCode_                              // recvbuf
				,1                                      // count
				,Teuchos::RawMPITraits<Index>::type()   // datatype
				,MPI_SUM                                // op
				,mpiComm                                // comm
				);
			// Set the default localOffset automatically
			Index localOffset = localSubDim;
			MPI_Scan(
				&localOffset                            // sendbuf
				,&defaultLocalOffset_                   // recvbuf
				,1                                      // count
				,Teuchos::RawMPITraits<Index>::type()   // datatype
				,MPI_SUM                                // op
				,mpiComm                                // comm
				);
			defaultLocalOffset_ -= localSubDim;
			//int procRank; MPI_Comm_rank( mpiComm, &procRank );
			//std::cout << "\nMPIVectorSpaceBase<Scalar>::updateState(): procRank = " << procRank << ", defaultLocalOffset = " << defaultLocalOffset_ << std::endl;
			// 
			isInCore_ = false;  // This is not an inCore vector
#ifdef _DEBUG
			// Check that the vector does not have any ghost elements
			Index computedGlobalDim = 0;
			MPI_Allreduce(
				(void*)&localSubDim                     // sendbuf
				,&computedGlobalDim                     // recvbuf
				,1                                      // count
				,Teuchos::RawMPITraits<Index>::type()   // datatype
				,MPI_SUM                                // op
				,mpiComm                                // comm
				);
			TEST_FOR_EXCEPTION(
				computedGlobalDim != globalDim, std::logic_error
				,"MPIVectorSpaceBase<Scalar>::updateState(): Error, the computed "
				"global dimension of computedGlobalDim = " << computedGlobalDim
				<< " is not equal to the reported global dimension of globalDim = this->dim() = "
				<< globalDim << "!"
				);
#endif // _DEBUG
		}
		else {
#endif // RTOp_USE_MPI
			// This is a serial or a locally-replicated parallel
			// vector space.
			mapCode_ = localSubDim;
			isInCore_ = true;
			defaultLocalOffset_ = 0;
#ifdef RTOp_USE_MPI
		}
#endif
	}
    else {
		mapCode_  = -1;     // Uninitialized!
		isInCore_ = false;
		defaultLocalOffset_ = -1;
	}
}
	
} // end namespace TSFCoreo

#endif // TSFCORE_MPI_VECTOR_SPACE_BASE_HPP
