// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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

#ifndef TPETRA_OMNIPLATFORMDATA_HPP
#define TPETRA_OMNIPLATFORMDATA_HPP

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Object.hpp"
//#include "Tpetra_OmniPlatform.hpp"
#ifdef TPETRA_MPI
#include <mpi.h>
#endif
#ifdef TPETRA_THREADED_MPI
// -- includes for a threaded MPI would go here --
#endif

namespace Tpetra {

	class OmniPlatformData : public Object {
		friend class OmniPlatform;
	public:
		// serial constructor
		OmniPlatformData() 
			: Object("Tpetra::OmniPlatformData")
			, serialEnabled_(true)
			, mpiEnabled_(false)
			, threadedMpiEnabled_(false)
		{}

#ifdef TPETRA_MPI
		// MPI constructor
		OmniPlatformData(MPI_Comm comm)
			: Object("Tpetra::OmniPlatformData")
			, serialEnabled_(true)
			, mpiEnabled_(true)
			, threadedMpiEnabled_(false)
			, MpiComm_(comm)
		{}
#endif

#ifdef TPETRA_THREADED_MPI
		// -- A constructor for a threaded MPI would go here --
#endif

		// destructor. no heap-data, so no need to override
		~OmniPlatformData() {}

	protected:

		bool serialEnabled_;
		bool mpiEnabled_;
		bool threadedMpiEnabled_;
		
#ifdef TPETRA_MPI
		MPI_Comm MpiComm_;
#endif

#ifdef TPETRA_THREADED_MPI
		// -- Data members for a threaded MPI would go here --
#endif

	private:

		//! Copy constructor (declared but not defined, do not use)
		OmniPlatformData(OmniPlatformData const& Source);

		//! Assignment operator (declared but not defined, do not use)
		OmniPlatformData& operator = (OmniPlatformData const& Source);

	}; // class OmniPlatformData

} // namespace Tpetra

#endif // TPETRA_OMNIPLATFORMDATA_HPP
