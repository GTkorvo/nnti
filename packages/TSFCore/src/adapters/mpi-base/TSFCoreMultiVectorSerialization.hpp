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

// /////////////////////////////////////////////////////////////////////////////
// TSFCoreMultiVectorSerialization.hpp

#ifndef TSFCORE_MULTI_VECTOR_SERIALIZATION_HPP
#define TSFCORE_MULTI_VECTOR_SERIALIZATION_HPP

#include "TSFCoreMultiVectorSerializationDecl.hpp"
#include "TSFCoreMPIVectorSpaceBase.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreExplicitMultiVectorView.hpp"

namespace TSFCore {

template<class Scalar>
MultiVectorSerialization<Scalar>::MultiVectorSerialization(
	const bool  binaryMode
	)
	:binaryMode_(binaryMode)
{}

template<class Scalar>
void MultiVectorSerialization<Scalar>::serialize( const MultiVector<Scalar>& mv, std::ostream& out ) const
{
	Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >
		mpi_vec_spc = Teuchos::rcp_dynamic_cast<const MPIVectorSpaceBase<Scalar> >(mv.range());
	out.precision(std::numeric_limits<Scalar>::digits10+4);
	if( mpi_vec_spc.get() ) {
		// This is a mpi-based vector space so let's just write the local
		// multi-vector elements (row-by-row).
		const Index
			localOffset = mpi_vec_spc->localOffset(),
			localSubDim = mpi_vec_spc->localSubDim();
		const Range1D localRng( localOffset+1, localOffset+localSubDim ); 
		ExplicitMultiVectorView<Scalar> local_mv(mv,localRng,Range1D());
		out << localSubDim << " " << local_mv.numSubCols() << std::endl;
		if( binaryMode() ) {
			// Write column-wise for better cache performance
			for( Index j = 1; j <= local_mv.numSubCols(); ++j )
				out.write( reinterpret_cast<const char*>(&local_mv(1,j)), sizeof(Scalar)*localSubDim );
		}
		else {
			// Write row-wise for better readability
			for( Index i = 1; i <= localSubDim; ++i ) {
				out << " " << i;
				for( Index j = 1; j <= local_mv.numSubCols(); ++j ) {
					out << " " << local_mv(i,j);
				}
				out << std::endl;
			}
		}
	}
	else {
		//  This is a serial (or locally replicated) vector space so
		// just write all of the multi-vector elements here.
		TEST_FOR_EXCEPTION( true, std::logic_error, "Does not handle non-MPI spaces yet" );
	}
}

template<class Scalar>
void MultiVectorSerialization<Scalar>::unserialize( std::istream& in, MultiVector<Scalar>* mv ) const
{
	Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >
		mpi_vec_spc = Teuchos::rcp_dynamic_cast<const MPIVectorSpaceBase<Scalar> >(mv->range());
	if( mpi_vec_spc.get() ) {
		// This is a mpi-based vector space so let's just read the local
		// multi-vector elements (row-by-row).
		const Index
			localOffset = mpi_vec_spc->localOffset(),
			localSubDim = mpi_vec_spc->localSubDim();
		const Range1D localRng( localOffset+1, localOffset+localSubDim ); 
		ExplicitMutableMultiVectorView<Scalar> local_mv(*mv,localRng,Range1D());
#ifdef _DEBUG
		TEST_FOR_EXCEPTION( !in, std::logic_error, "Error, premature end of input!"	);
#endif
		Index localSubDim_in;
		in >> localSubDim_in;
#ifdef _DEBUG
		TEST_FOR_EXCEPTION(
			localSubDim != localSubDim_in, std::logic_error
			, "Error, localSubDim = "<<localSubDim<<" does not match the readin value of "
			"localSubDim_in = "<<localSubDim_in<<"!"
			);
#endif
		Index numSubCols_in;
		in >> numSubCols_in;
#ifdef _DEBUG
		TEST_FOR_EXCEPTION(
			local_mv.numSubCols() != numSubCols_in, std::logic_error
			, "Error, numSubCols = "<<local_mv.numSubCols()<<" does not match the readin value of "
			"numSubCols_in = "<<numSubCols_in<<"!"
			);
#endif
		// Get rid of extra newline after first line
		in >> std::ws;
		// Get the elements
		if( binaryMode() ) {
			// Column-wise
			for( Index j = 1; j <= local_mv.numSubCols(); ++j )
				in.read( reinterpret_cast<char*>(&local_mv(1,j)), sizeof(Scalar)*localSubDim );
		}
		else {
			// Row-wise
			for( Index i = 1; i <= localSubDim; ++i ) {
#ifdef _DEBUG
				TEST_FOR_EXCEPTION( !in, std::logic_error, "Error, premature end of input!"	);
#endif
				Index i_in;
				in >> i_in;
#ifdef _DEBUG
				TEST_FOR_EXCEPTION(
					i != i_in, std::logic_error
					, "Error, i = "<<i<<" does not match the readin value of "
					"i_in = "<<i_in<<"!"
					);
#endif
				for( Index j = 1; j <= local_mv.numSubCols(); ++j ) {
#ifdef _DEBUG
					TEST_FOR_EXCEPTION( !in, std::logic_error, "Error, premature end of input!"	);
#endif
					in >> local_mv(i,j);
				}
			}
		}
	}
	else {
		//  This is a serial (or locally replicated) vector space so
		// just read all of the multi-vector elements here.
		TEST_FOR_EXCEPTION( true, std::logic_error, "Does not handle non-MPI spaces yet" );
	}
}

} // end namespace TSFCore

#endif // TSFCORE_MULTI_VECTOR_SERIALIZATION_HPP
