/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_OPENMP_HPP
#define KOKKOS_OPENMP_HPP

#include <Kokkos_Macros.hpp>

#if defined( KOKKOS_HAVE_OPENMP ) && defined( _OPENMP )

#include <omp.h>

#include <cstddef>
#include <iosfwd>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_Parallel.hpp>
#include <Kokkos_Layout.hpp>
#include <impl/Kokkos_Tags.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
class OpenMPexec ;
} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/// \class OpenMP
/// \brief Kokkos device for multicore processors in the host memory space.
class OpenMP {
public:
  //------------------------------------
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  //! The tag (what type of kokkos_object is this).
  typedef Impl::DeviceTag       kokkos_tag ;
  typedef OpenMP                device_type ;
  typedef HostSpace::size_type  size_type ;
  typedef HostSpace             memory_space ;
  typedef OpenMP                scratch_memory_space ;
  typedef LayoutRight           array_layout ;
  typedef OpenMP                host_mirror_device_type ;

  //@}
  //------------------------------------
  //! \name Functions that all Kokkos devices must implement.
  //@{

  inline static bool in_parallel() { return omp_in_parallel(); }

  /** \brief  Set the device in a "sleep" state. A noop for OpenMP.  */
  static bool sleep();

  /** \brief Wake the device from the 'sleep' state. A noop for OpenMP. */
  static bool wake();

  /** \brief Wait until all dispatched functors complete. A noop for OpenMP. */
  static void fence() {}

  /// \brief Print configuration information to the given output stream.
  static void print_configuration( std::ostream & , const bool detail = false );

  /// \brief Free any resources being consumed by the device.
  static void finalize();

  /** \brief  Initialize the device.
   *
   *  1) If the hardware locality library is enabled and OpenMP has not
   *     already bound threads then bind OpenMP threads to maximize
   *     core utilization and group for memory hierarchy locality.
   *
   *  2) Allocate a HostThread for each OpenMP thread to hold its
   *     topology and fan in/out data.
   */
  static void initialize( unsigned thread_count = 0 ,
                          unsigned use_numa_count = 0 ,
                          unsigned use_cores_per_numa = 0 );

  static int is_initialized();

  KOKKOS_FUNCTION static unsigned team_max();
  KOKKOS_FUNCTION static unsigned team_recommended();
  KOKKOS_INLINE_FUNCTION static unsigned hardware_thread_id();
  KOKKOS_INLINE_FUNCTION static unsigned max_hardware_threads();
  //@}
  //------------------------------------
  //! \name Function for the functor device interface */

  KOKKOS_INLINE_FUNCTION void * get_shmem( const int size ) const ;

  explicit KOKKOS_INLINE_FUNCTION OpenMP( Impl::OpenMPexec & );

  //------------------------------------

private:

  Impl::OpenMPexec & m_exec ;

};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

#include <OpenMP/Kokkos_OpenMPexec.hpp>
#include <OpenMP/Kokkos_OpenMP_Parallel.hpp>

/*--------------------------------------------------------------------------*/

#endif /* #if defined( KOKKOS_HAVE_OPENMP ) && defined( _OPENMP ) */
#endif /* #ifndef KOKKOS_OPENMP_HPP */


