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

#ifndef __TSQR_Test_MgsTest_hpp
#define __TSQR_Test_MgsTest_hpp

#include <Tsqr_ConfigDefs.hpp>

#include <Tsqr_Mgs.hpp>
#ifdef HAVE_TSQR_INTEL_TBB
#  include <TbbTsqr_TbbMgs.hpp>
#endif // HAVE_TSQR_INTEL_TBB
#include <Tsqr_TestSetup.hpp>
#include <Tsqr_GlobalVerify.hpp>
#include <Tsqr_printGlobalMatrix.hpp>
#include <Tsqr_verifyTimerConcept.hpp>

#include <Teuchos_RCP.hpp>

#include <algorithm>
#include <sstream>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <utility>


namespace TSQR {
  namespace Test {

    static std::string
    mgs_human_readable_name (const std::string& which)
    {
      if (which == "MpiSeqMGS")
	return std::string ("MPI parallel / sequential MGS");
      else if (which == "MpiTbbMGS")
	{
#ifdef HAVE_TSQR_INTEL_TBB
	  return std::string ("MPI parallel / TBB parallel MGS");
#else
	  throw std::logic_error("MGS not built with Intel TBB support");
#endif // HAVE_TSQR_INTEL_TBB
	}
      else 
	throw std::logic_error("Unknown MGS implementation type \"" + which + "\"");
    }

    template< class MgsType >
    class MgsVerifier {
    public:
      typedef MgsType mgs_type;
      typedef typename MgsType::ordinal_type ordinal_type;
      typedef typename MgsType::scalar_type scalar_type;
      typedef Matrix< ordinal_type, scalar_type > matrix_type;
      typedef MessengerBase< scalar_type > messenger_type;
      typedef Teuchos::RCP< messenger_type > messenger_ptr;

      static void
      verify (mgs_type& orthogonalizer,
	      const messenger_ptr& messenger,
	      matrix_type& Q_local,
	      matrix_type& R,
	      const bool b_debug = false)
      {
	using std::cerr;
	using std::endl;

	// Factor the (copy of the) matrix.  On output, the explicit Q
	// factor (of A_local) is in Q_local and the R factor is in R.
	orthogonalizer.mgs (Q_local.nrows(), Q_local.ncols(), 
			    Q_local.get(), Q_local.lda(),
			    R.get(), R.lda());
	if (b_debug)
	  {
	    messenger->barrier();
	    if (messenger->rank() == 0)
	      cerr << "-- Finished MGS::mgs" << endl;
	  }
      }
    };

    template< class Ordinal, class Scalar, class Generator >
    void
    verifyMgs (const std::string& which,
	       Generator& generator,
	       const Ordinal nrows_global,
	       const Ordinal ncols,
	       const Teuchos::RCP< MessengerBase< Ordinal > >& ordinalComm,
	       const Teuchos::RCP< MessengerBase< Scalar > >& scalarComm,
	       const int num_cores,
	       const bool human_readable,
	       const bool b_debug)
    {
      typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
      using std::cerr;
      using std::cout;
      using std::endl;

      const bool b_extra_debug = false;
      const int nprocs = scalarComm->size();
      const int my_rank = scalarComm->rank();
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "mgs_verify:" << endl;
	  scalarComm->barrier();
	}
      const Ordinal nrows_local = numLocalRows (nrows_global, my_rank, nprocs);

      // Set up storage for the test problem
      Matrix< Ordinal, Scalar > A_local (nrows_local, ncols);
      if (std::numeric_limits< Scalar >::has_quiet_NaN)
	A_local.fill (std::numeric_limits< Scalar >::quiet_NaN());
      Matrix< Ordinal, Scalar > R (ncols, ncols, Scalar(0));

      // Generate the test problem.
      distributedTestProblem (generator, A_local, ordinalComm.get(), scalarComm.get());
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Generated test problem." << endl;
	}

      // Make sure that the test problem (the matrix to factor) was
      // distributed correctly.
      if (b_extra_debug && b_debug)
	{
	  if (my_rank == 0)
	    cerr << "Test matrix A:" << endl;
	  scalarComm->barrier();
	  printGlobalMatrix (cerr, A_local, scalarComm.get(), ordinalComm.get());
	  scalarComm->barrier();
	}

      // Factoring the matrix stored in A_local overwrites it, so we
      // copy A_local into Q_local.  MGS orthogonalization does not
      // support contiguously stored cache blocks, unlike TSQR, so we
      // don't have to consider whether or not to rearrange cache
      // blocks here (unlike with TSQR).
      Matrix< Ordinal, Scalar > Q_local (A_local);

      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Starting verification" << endl;
	}

      if (which == "MpiTbbMGS")
	{
#ifdef HAVE_TSQR_INTEL_TBB
	  typedef TSQR::TBB::TbbMgs< Ordinal, Scalar > mgs_type;
	  mgs_type mgser (scalarComm);
	  MgsVerifier< mgs_type >::verify (mgser, scalarComm, Q_local, R, b_debug);
#else
	  throw std::logic_error("MGS not built with Intel TBB support");
#endif // HAVE_TSQR_INTEL_TBB
	}
      else if (which == "MpiSeqMGS")
	{
	  typedef MGS< Ordinal, Scalar > mgs_type;
	  mgs_type mgser (scalarComm);
	  MgsVerifier< mgs_type >::verify (mgser, scalarComm, Q_local, R, b_debug);
	}
      else
	throw std::logic_error ("Invalid MGS implementation type \"" + which + "\"");

      // Print out the Q and R factors
      if (b_extra_debug && b_debug)
	{
	  if (my_rank == 0)
	    cerr << endl << "Q factor:" << endl;
	  scalarComm->barrier ();
	  printGlobalMatrix (cerr, A_local, scalarComm.get(), ordinalComm.get());
	  scalarComm->barrier ();
	  if (my_rank == 0)
	    {
	      cerr << endl << "R factor:" << endl;
	      print_local_matrix (cerr, ncols, ncols, R.get(), R.lda());
	      cerr << endl;
	    }
	  scalarComm->barrier ();
	}

      // Test accuracy of the resulting factorization
      std::vector< magnitude_type > results = 
	global_verify (nrows_local, ncols, A_local.get(), A_local.lda(),
		       Q_local.get(), Q_local.lda(), R.get(), R.lda(), 
		       scalarComm.get());
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Finished global_verify" << endl;
	  scalarComm->barrier();
	}

      // Print the results on Proc 0.
      if (my_rank == 0)
	{
	  if (human_readable)
	    {
	      cout << mgs_human_readable_name(which) << endl
		   << "# rows = " << nrows_global << endl
		   << "# columns = " << ncols << endl
		   << "# MPI processes = " << nprocs << endl;
	      if (which == "MpiTbbTSQR")
		cout << "# cores per process = " << num_cores << endl;
	      cout << "Absolute residual $\\|A - Q*R\\|_2: "
		   << results[0] << endl
		   << "Absolute orthogonality $\\|I - Q^T*Q\\|_2$: " 
		   << results[1] << endl
		   << "Test matrix norm $\\| A \\|_F$: "
		   << results[2] << endl
		   << endl;
	    }
	  else
	    {
	      cout << which
		   << "," << nrows_global
		   << "," << ncols
		   << "," << nprocs;
	      if (which == "MpiTbbTSQR")
		cout << "," << num_cores << endl;
	      cout << "," << results[0] 
		   << "," << results[1]
		   << "," << results[2]
		   << endl;
	    }
	}
    }


    template< class MgsBase, class TimerType >
    static double // returns timing in s
    do_mgs_benchmark (MgsBase& orthogonalizer,
		      Matrix< typename MgsBase::ordinal_type, typename MgsBase::scalar_type >& Q_local,
		      Matrix< typename MgsBase::ordinal_type, typename MgsBase::scalar_type >& R,
		      const int num_trials,
		      const bool human_readable)
    {
      typedef typename MgsBase::ordinal_type ordinal_type;
      using std::cout;

      TSQR::Test::verifyTimerConcept< TimerType >();

      const ordinal_type nrows_local = Q_local.nrows();
      const ordinal_type ncols = Q_local.ncols();

      // Benchmark MGS for ntrials trials.  The answer (the numerical
      // results of the factorization) is only valid if ntrials == 1,
      // but this is a benchmark and not a verification routine.  Call
      // mgs_verify() if you want to determine whether MGS computes
      // the right answer.
      //
      // Name of timer doesn't matter here; we only need the timing.
      TimerType timer("MGS"); 
      timer.start();
      for (int trial_num = 0; trial_num < num_trials; ++trial_num)
	{
	  // Orthogonalize the columns of A using MGS.  Don't worry about
	  // the fact that we're overwriting the input; this is a
	  // benchmark, not a numerical verification test.  (We have the
	  // latter implemented as mgs_verify() in this file.)
	  orthogonalizer.mgs (nrows_local, ncols, Q_local.get(),
			      Q_local.lda(), R.get(), R.lda());
	  // Timings in debug mode likely won't make sense, because
	  // Proc 0 is outputting the debug messages to cerr.
	  // Nevertheless, we don't put any "if(b_debug)" calls in the
	  // timing loop.
	}
      // Compute the resulting total time (in seconds) to execute
      // num_trials runs of :mgs().  The time may differ on different
      // MPI processes.
      const double mgs_timing = timer.stop();
      return mgs_timing;
    }

    template< class Ordinal, class Scalar, class Generator, class TimerType >
    void
    benchmarkMgs (const std::string& which,
		  Generator& generator,
		  const int ntrials,
		  const Ordinal nrows_global,
		  const Ordinal ncols,
		  const Teuchos::RCP< MessengerBase< Ordinal > >& ordinalComm,
		  const Teuchos::RCP< MessengerBase< Scalar > >& scalarComm,
		  const int num_cores,
		  const bool human_readable,
		  const bool b_debug)
    {
      typedef typename Teuchos::ScalarTraits< Scalar >::magnitudeType magnitude_type;
      using std::cerr;
      using std::cout;
      using std::endl;

      TSQR::Test::verifyTimerConcept< TimerType >();

      const bool b_extra_debug = false;
      const int nprocs = scalarComm->size();
      const int my_rank = scalarComm->rank();
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "mgs_benchmark:" << endl;
	  scalarComm->barrier();
	}
      const Ordinal nrows_local = numLocalRows (nrows_global, my_rank, nprocs);
   
      // Set up storage for the test problem.
      Matrix<Ordinal, Scalar> A_local (nrows_local, ncols);
      if (std::numeric_limits< Scalar >::has_quiet_NaN)
	A_local.fill (std::numeric_limits< Scalar >::quiet_NaN());
      Matrix<Ordinal, Scalar> R (ncols, ncols, Scalar(0));

      // Generate the test problem.
      distributedTestProblem (generator, A_local, ordinalComm.get(), scalarComm.get());
      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Generated test problem." << endl;
	}

      // Make sure that the test problem (the matrix to factor) was
      // distributed correctly.
      if (b_extra_debug && b_debug)
	{
	  if (my_rank == 0)
	    cerr << "Test matrix A:" << endl;
	  scalarComm->barrier ();
	  printGlobalMatrix (cerr, A_local, scalarComm.get(), ordinalComm.get());
	  scalarComm->barrier ();
	}

      // Factoring the matrix stored in A_local overwrites it, so we
      // make a copy of A_local.  MGS orthogonalization does not
      // support contiguously stored cache blocks, unlike TSQR, so we
      // don't have to consider whether or not to rearrange cache
      // blocks here (unlike with TSQR).
      Matrix< Ordinal, Scalar > Q_local (A_local);

      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Starting timing loop" << endl;
	}

      // Set up MGS and run the benchmark.
      double mgs_timing; // Total run time in seconds of all ntrials trials
      if (which == "MpiTbbMGS")
	{
#ifdef HAVE_TSQR_INTEL_TBB
	  typedef TSQR::TBB::TbbMgs< Ordinal, Scalar > mgs_type;
	  mgs_type mgser (scalarComm);
	  mgs_timing = do_mgs_benchmark< mgs_type, TimerType > (mgser, Q_local, R, 
								ntrials, human_readable);
#else
	  throw std::logic_error("MGS not built with Intel TBB support");
#endif // HAVE_TSQR_INTEL_TBB
	}
      else if (which == "MpiSeqMGS")
	{
	  typedef MGS< Ordinal, Scalar > mgs_type;
	  mgs_type mgser (scalarComm);
	  mgs_timing = do_mgs_benchmark< mgs_type, TimerType > (mgser, Q_local, R, 
								ntrials, human_readable);
	}
      else
	throw std::logic_error ("Invalid MGS implementation type \"" + which + "\"");

      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Finished timing loop" << endl;
	}

      // Find the min and max MGS timing on all processors.
      const double min_mgs_timing = scalarComm->globalMin (mgs_timing);
      const double max_mgs_timing = scalarComm->globalMax (mgs_timing);

      if (b_debug)
	{
	  scalarComm->barrier();
	  if (my_rank == 0)
	    cerr << "-- Computed min and max timings" << endl;
	}

      // Print the results on Proc 0.
      if (my_rank == 0)
	{
	  if (human_readable)
	    {
	      cout << mgs_human_readable_name(which) << ":" << endl
		   << "# rows = " << nrows_global << endl
		   << "# columns = " << ncols << endl
		   << "# MPI processes = " << nprocs << endl;
	      if (which == "MpiTbbTSQR")
		cout << "# cores per process = " << num_cores << endl;
	      cout << "# trials = " << ntrials << endl
		   << "Min total time (s) over all MPI processes = " 
		   << min_mgs_timing << endl
		   << "Max total time (s) over all MPI processes = " 
		   << max_mgs_timing << endl
		   << endl;
	    }
	  else
	    {
	      cout << which
		   << "," << nrows_global
		   << "," << ncols 
		   << "," << nprocs;
	      if (which == "MpiTbbTSQR")
		cout << "," << num_cores << endl;
	      cout << "," << ntrials 
		   << "," << min_mgs_timing 
		   << "," << max_mgs_timing 
		   << endl;
	    }
	}
    }
    

  } // namespace Test
} // namespace TSQR

#endif // __TSQR_Test_MgsTest_hpp
