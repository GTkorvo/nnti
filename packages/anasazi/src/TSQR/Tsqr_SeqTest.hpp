// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2010) Sandia Corporation
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

#ifndef __TSQR_Test_SeqTest_hpp
#define __TSQR_Test_SeqTest_hpp

#include <Tsqr_Config.hpp>

#include <cstring> // size_t definition
#include <string>
#include <iostream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace Test {

    /// Test the accuracy of sequential TSQR on an nrows by ncols
    /// matrix (using the given cache block size (in bytes)), and
    /// print the results to stdout.
    void
    verifySeqTsqr (std::ostream& out,
		   const int nrows, 
		   const int ncols, 
		   const size_t cache_block_size,
		   const bool test_complex_arithmetic,
		   const bool save_matrices,
		   const bool contiguous_cache_blocks,
		   const bool printFieldNames,
		   const bool human_readable = false,
		   const bool b_debug = false);

    /// Test the accuracy of LAPACK's QR factorization on an nrows by
    /// ncols matrix, and print the results to stdout.
    void
    verifyLapack (std::ostream& out,
		  const int nrows, 
		  const int ncols, 
		  const bool test_complex_arithmetic,
		  const std::string& additionalFieldNames,
		  const std::string& additionalData,
		  const bool printFieldNames,
		  const bool human_readable,
		  const bool b_debug = false);
        
    /// Test the runtime (over ntrials trials) of sequential TSQR, on
    /// an nrows by ncols matrix (using the given cache block size (in
    /// bytes)), and print the results to stdout.
    ///
    /// \param human_readable [in] If true, print the benchmark
    /// results to stdout in human-readable format.  Otherwise, print
    /// them as two rows of comma-delimited ASCII, in an abbreviated
    /// format suitable for automatic processing.
    void
    benchmarkSeqTsqr (std::ostream& out,
		      const int numRows,
		      const int numCols,
		      const int numTrials,
		      const size_t cacheBlockSize,
		      const bool contiguousCacheBlocks,
		      const bool testComplex,
		      const bool printFieldNames,
		      const bool humanReadable);

    /// Test the runtime (over numTrials trials) of LAPACK QR (_GEQRF
    /// + _ORGQR), on a numRows by numCols matrix, and print the
    /// results to out.
    ///
    /// \param humanReadable [in] If true, print the benchmark results
    ///   to out in human-readable format.  Otherwise, print them as
    ///   two rows of comma-delimited ASCII, in an abbreviated format
    ///   suitable for automatic processing.
    void
    benchmarkLapack (std::ostream& out,
		     const int numRows,
		     const int numCols,
		     const int numTrials,
		     const bool testComplex,
		     const std::string& additionalFieldNames,
		     const std::string& additionalData,
		     const bool printFieldNames,
		     const bool humanReadable);

  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_SeqTest_hpp
