//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#ifndef __TSQR_TrivialMessenger_hpp
#define __TSQR_TrivialMessenger_hpp

#include <Tsqr_MessengerBase.hpp>

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <vector>


namespace TSQR { 

  /// \class TrivialMessenger
  /// \brief Noncommunicating "communication" object for TSQR.
  ///
  /// The internode parallel part of TSQR communicates via a
  /// MessengerBase<Datum> interface.  TrivialMessenger<Datum>
  /// implements that interface by acting as if running on MPI with
  /// only one rank, though it doesn't require MPI support to build.
  ///
  /// \tparam Datum A class with value-type semantics, whose instances
  ///   are less-than comparable.
  template<class Datum>
  class TrivialMessenger : public MessengerBase<Datum> {
  public:
    //! Trivial / default constructor, since no member data.
    TrivialMessenger () {}

    //! Virtual destructor for memory safety of derived classes.
    virtual ~TrivialMessenger() {}

    /// \brief Send sendData[0:sendCount-1] to process destProc.
    ///
    /// \param sendData [in] Array of value-type elements to send
    /// \param sendCount [in] Number of elements in the array
    /// \param destProc [in] Rank of destination process
    /// \param tag [in] MPI tag (ignored)
    void 
    send (const Datum sendData[], 
	  const int sendCount, 
	  const int destProc, 
	  const int tag) 
    {}

    /// \brief Receive recvData[0:recvCount-1] from process srcProc.
    ///
    /// \param recvData [out] Array of value-type elements to receive
    /// \param recvCount [in] Number of elements to receive in the array
    /// \param srcProc [in] Rank of sending process
    /// \param tag [in] MPI tag (ignored)
    void 
    recv (Datum recvData[], 
	  const int recvCount, 
	  const int srcProc, 
	  const int tag) 
    {}

    /// \brief Exchange data between processors.
    ///
    /// Exchange sencRecvCount elements of sendData with processor
    /// destProc, receiving the result into recvData.  Assume that
    /// sendData and recvData do not alias one another.
    ///
    /// \param sendData [in] Array of value-type elements to send
    /// \param recvData [out] Array of value-type elements to
    ///   receive.  Caller is responsible for making sure that
    ///   recvData does not alias sendData.
    /// \param sendRecvCount [in] Number of elements to send and
    ///   receive in the array
    /// \param destProc [in] The "other" process' rank (to which
    ///   this process is sending data, and from which this process is
    ///   receiving data)
    /// \param tag [in] MPI tag (ignored)
    void 
    swapData (const Datum sendData[], 
	      Datum recvData[], 
	      const int sendRecvCount, 
	      const int destProc, 
	      const int tag)
    {
      if (destProc != rank())
	{
	  std::ostringstream os;
	  os << "Destination rank " << destProc << " is invalid.  The only "
	     << "valid rank for TSQR::TrivialMessenger is 0 (zero).";
	    throw std::invalid_argument (os.str());
	}
      else if (sendRecvCount < 0)
	{
	  std::ostringstream os;
	  os << "sendRecvCount = " << sendRecvCount << " is invalid: "
	     << "only nonnegative values are allowed.";
	  throw std::invalid_argument (os.str());
	}
      else if (sendRecvCount == 0)
	return; // No data to exchange
      else 
	safeCopy (sendData, recvData, sendRecvCount);
    }

    //! Sum inDatum on all processors, and return the result.
    Datum 
    globalSum (const Datum& inDatum) 
    {
      Datum outDatum (inDatum);
      return outDatum;
    }

    /// \brief Compute the global minimum over all processors.
    ///
    /// Assumes that Datum objects are less-than comparable.
    Datum 
    globalMin (const Datum& inDatum)
    {
      Datum outDatum (inDatum);
      return outDatum;
    }

    /// \brief Compute the global maximum over all processors.
    ///
    /// Assumes that Datum objects are less-than comparable.
    Datum 
    globalMax (const Datum& inDatum)
    {
      Datum outDatum (inDatum);
      return outDatum;
    }

    //! Sum inData[0:count-1] over all processors into outData.
    void
    globalVectorSum (const Datum inData[], 
		     Datum outData[], 
		     const int count) 
    {
      safeCopy (inData, outData, count);
    }

    //! Broadcast data[0:count-1] from root to all processors.
    void
    broadcast (Datum data[], 
	       const int count,
	       const int root)
    {}

    //! Return this process' rank.
    int rank () const { return 0; }

    //! Return the total number of ranks in the communicator.
    int size () const { return 1; }

    //! Execute a barrier over the communicator.
    void barrier () const { }

  private:

    /// \brief Copy count elements of inData into outData.
    /// 
    /// Attempt to detect aliasing, and use a method appropriate for
    /// either the nonaliased or the aliased case.
    void
    safeCopy (const Datum inData[],
	      Datum outData[],
	      const int count) 
    {
      // Check for nonaliasing of inData and outData.
      if (&inData[count-1] < &outData[0] ||
	  &outData[count-1] < &inData[0])
	// The arrays don't overlap, so we can call std::copy.
	// std::copy assumes that the third argument does not
	// point to an element in the range of the first two
	// arguments.
	std::copy (inData, inData+count, outData);
      else
	{
	  // If inData and outData do alias one another, use
	  // the buffer as intermediate scratch space.
	  buf_.resize (count);
	  std::copy (inData, inData+count, buf_.begin());
	  std::copy (buf_.begin(), buf_.end(), outData);
	}
    }

    /// Buffer to guard against incorrect behavior for aliased arrays.
    /// Space won't be allocated unless needed.
    std::vector<Datum> buf_;
  };
} // namespace TSQR

#endif // __TSQR_TrivialMessenger_hpp

