//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2010 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

///
/// \file Epetra_TsqrMessenger.hpp
/// 
/// Method for fetching TSQR::MessengerBase instance for use in the
/// Epetra_MultiVector adaptor for TSQR.
///
/// \warning Since TSQR lives in Kokkos and depends on Teuchos, this
///   file should only be included in the Epetra build when the Kokkos
///   and Teuchos packages are enabled.
///

#ifndef __Epetra_TsqrMessenger_hpp
#define __Epetra_TsqrMessenger_hpp

#include <Kokkos_ConfigDefs.hpp> // HAVE_KOKKOS_TSQR

#ifdef HAVE_KOKKOS_TSQR
#  include <Epetra_ConfigDefs.h> // EPETRA_MPI
#  include <Epetra_Comm.h>
#  ifdef EPETRA_MPI
#    include <Epetra_MpiComm.h>
#    include <Tsqr_MpiMessenger.hpp>
#  endif // EPETRA_MPI
#  include <Tsqr_TrivialMessenger.hpp>
#  include <algorithm>
#  include <utility> // std::pair

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Epetra { 

#ifdef EPETRA_MPI
    /// If the input Epetra_Comm is really an Epetra_MpiComm, return its
    /// raw MPI_Comm MPI communicator object.  Otherwise, return
    /// MPI_COMM_NULL.  (The Epetra_Comm interface doesn't define sends
    /// and receives, which TSQR needs.  That's why TSQR wants the raw
    /// MPI_COMM object.)
    ///
    /// \return (The MPI_Comm object, and whether it's valid)
    std::pair< MPI_Comm, bool >
    extractRawMpiComm (const Teuchos::RCP< const Epetra_MpiComm >& pComm)
    {
      MPI_Comm rawMpiComm = MPI_COMM_NULL;
      bool haveMpiComm = false;

      Teuchos::RCP< const Epetra_MpiComm > pMpiComm = 
	Teuchos::rcp_dynamic_cast (pComm, false);
      if (pMpiComm.get() == NULL)
	{
	  // See if the input Epetra_Comm is really an
	  // Epetra_MpiSmpComm.  If so, pull out its raw MPI_Comm
	  // object and use that.
	  Teuchos::RCP< const Epetra_MpiSmpComm > pMpiSmpComm = 
	    Teuchos::rcp_dynamic_cast (pComm, false);
	  if (pMpiSmpComm.get() != NULL)
	    {
	      rawMpiComm = pMpiSmpComm->Comm();
	      haveMpiComm = true;
	    }
	}
      else
	{
	  rawMpiComm = pMpiComm->Comm();
	  haveMpiComm = true;
	}
      return std::make_pair (rawMpiComm, haveMpiComm);
    }
#endif // EPETRA_MPI


    /// Given a pointer to an Epetra_Comm object, return a pointer to
    /// a TSQR::MessengerBase< Datum > (actually, a pointer to an
    /// appropriate subclass thereof, depending on whether the given
    /// Epetra_Comm uses MPI).
    template< class Datum >
    Teuchos::RCP< MessengerBase< Datum > >
    makeTsqrMessenger (const Teuchos::RCP< const Epetra_Comm >& pComm)
    {
      using Teuchos::RCP;
      using Teuchos::rcp_implicit_cast;
      typedef TSQR::MessengerBase< Datum > base_mess_type;

#ifdef EPETRA_MPI
      // If the Epetra_Comm wraps MPI, then extract the raw MPI_Comm
      // object and use it to construct a TSQR messenger object that
      // wraps MPI.  Otherwise, return a TSQR messenger object that
      // wraps trivial communication.
      std::pair< MPI_Comm, bool > results = extractRawMpiComm (pComm);
      if (results.second == true)
	{
	  typedef TSQR::MPI::MpiMessenger< Datum > mess_type;

	  RCP< mess_type > pMess (new mess_type (results.first));
	  RCP< base_mess_type > pMessBase = rcp_implicit_cast (pMess);
	  return pMessBase;
	}
      else
#endif // EPETRA_MPI
	{
	  typedef TSQR::TrivialMessenger< Datum > mess_type;
	  
	  RCP< mess_type > pMess (new mess_type);
	  RCP< base_mess_type > pMessBase = rcp_implicit_cast (pMess);
	  return pMessBase;
	}
    }

  } // namespace Epetra
} // namespace TSQR

#endif // HAVE_KOKKOS_TSQR
#endif // __Epetra_TsqrMessenger_hpp

