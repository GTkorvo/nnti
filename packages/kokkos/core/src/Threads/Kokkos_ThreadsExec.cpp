/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <limits>
#include <utility>
#include <iostream>
#include <Kokkos_Threads.hpp>
#include <Kokkos_hwloc.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
namespace {

ThreadsExec                  s_threads_process ;
ThreadsExec                * s_threads_exec[  ThreadsExec::MAX_THREAD_COUNT ];
std::pair<unsigned,unsigned> s_threads_coord[ ThreadsExec::MAX_THREAD_COUNT ];
std::string                  s_exception_msg ;

unsigned s_threads_count       = 0 ;
unsigned s_threads_per_numa    = 0 ;
unsigned s_threads_per_core    = 0 ;
unsigned s_threads_per_team_default = 0 ;

unsigned s_current_reduce_size = 0 ;
unsigned s_current_shared_size = 0 ;
unsigned s_current_team_alloc  = 0 ;
unsigned s_current_team_size   = 0 ;
unsigned s_current_league_size = 0 ;

void (* volatile s_current_function)( ThreadsExec & , const void * );
const void * volatile s_current_function_arg = 0 ;

struct Sentinel {
  Sentinel()
  {
    HostSpace::register_in_parallel( ThreadsExec::in_parallel );
  }

  ~Sentinel()
  {
    if ( s_threads_count ||
         s_threads_per_numa ||
         s_threads_per_core ||
         s_current_reduce_size ||
         s_current_shared_size ||
         s_current_function ||
         s_current_function_arg ||
         s_threads_exec[0] ) {
      std::cerr << "ERROR : Process exiting without calling Kokkos::Threads::terminate()" << std::endl ;
    }
  }
};

inline
unsigned fan_size( const unsigned rank , const unsigned size )
{
  const unsigned rank_rev = size - ( rank + 1 );
  unsigned count = 0 ;
  for ( unsigned n = 1 ; ( rank_rev + n < size ) && ! ( rank_rev & n ) ; n <<= 1 ) { ++count ; }
  return count ;
}

} // namespace
} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

inline
void ThreadsExec::set_team_relations()
{
  m_team_base        = 0 ;
  m_team_shared      = 0 ;
  m_team_shared_end  = 0 ;
  m_team_size        = 0 ;
  m_team_rank        = 0 ;
  m_team_fan_size    = 0 ;
  m_league_size      = 0 ;
  m_league_rank      = 0 ;
  m_league_end       = 0 ;

  const size_t league_size = s_current_league_size ;

  if ( league_size ) {
    // Execution is using device-team interface:

    const unsigned team_alloc    = s_current_team_alloc ;
    const unsigned team_size     = s_current_team_size ;
    const unsigned pool_rank_rev = m_pool_size - ( m_pool_rank + 1 );
    const unsigned team_rank_rev = pool_rank_rev % team_alloc ;

    // May be using fewer threads per team than a multiple of threads per core,
    // some threads will idle.

    if ( team_rank_rev < team_size ) {
      const size_t pool_league_size     = m_pool_size   / team_alloc ;
      const size_t pool_league_rank_rev = pool_rank_rev / team_alloc ;
      const size_t pool_league_rank     = pool_league_size - ( pool_league_rank_rev + 1 );

      m_team_base        = m_pool_base + team_alloc * pool_league_rank_rev ;
      m_team_shared      = (*m_team_base)->m_alloc_shared ;
      m_team_shared_end  = s_current_shared_size ;
      m_team_size        = team_size ;
      m_team_rank        = team_size - ( team_rank_rev + 1 );
      m_team_fan_size    = fan_size( m_team_rank , team_size );
      m_league_size      = league_size ;
      m_league_rank      = ( league_size *  pool_league_rank    ) / pool_league_size ;
      m_league_end       = ( league_size * (pool_league_rank+1) ) / pool_league_size ;
    }
  }
}

void ThreadsExec::driver(void)
{
  // If hardware locality library unavailable then pass in the rank.

  size_t thread_rank = (size_t) s_current_function_arg ;

  if ( s_threads_count <= thread_rank ) {
    thread_rank = Kokkos::hwloc::bind_this_thread( s_threads_count , s_threads_coord );
  }

  if ( s_threads_count <= thread_rank || 0 != ((ThreadsExec * volatile *)s_threads_exec)[ thread_rank ] ) {

    // An error occured. Inform process that thread is terminating
    s_threads_process.m_pool_state = ThreadsExec::Terminating ;

    return ;
  }

  {
    ThreadsExec this_thread ;

    this_thread.m_pool_state = ThreadsExec::Active ;

    // Try to protect against cache coherency failure by casting to volatile.
    ((ThreadsExec * volatile *)s_threads_exec)[ thread_rank ] = & this_thread ;
    // Really need a memory fence here.

    // Inform spawning process that the threads_exec entry has been set.
    s_threads_process.m_pool_state = ThreadsExec::Active ;

    while ( ThreadsExec::Active == this_thread.m_pool_state ) {

#if 0
      try {
        // Call work function
        (*s_current_function)( this_thread , s_current_function_arg );
      }
      catch( const std::exception & x ) {
        std::ostringstream msg ;
        msg << "Kokkos::Threads[" << thread_rank << "] Uncaught exeception : " << x.what() << std::endl ;
        s_exception_msg.append( msg.str() );
      }
      catch( ... ) {
        std::ostringstream msg ;
        msg << "Kokkos::Threads[" << thread_rank << "] Uncaught exeception"  << std::endl ;
        s_exception_msg.append( msg.str() );
      }
#else
        this_thread.set_team_relations();

        (*s_current_function)( this_thread , s_current_function_arg );
#endif

      // Deactivate thread and wait for reactivation
      this_thread.m_pool_state = ThreadsExec::Inactive ;
      wait_yield( this_thread.m_pool_state , ThreadsExec::Inactive );
    }

    s_threads_process.m_pool_state = ThreadsExec::Terminating ;

    ((ThreadsExec * volatile * )s_threads_exec)[ thread_rank ] = 0 ;
  }
}

void execute_function_noop( ThreadsExec & , const void * ) {}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

ThreadsExec::~ThreadsExec()
{
  m_pool_base   = 0 ;
  m_team_base   = 0 ;

  m_alloc_reduce     = 0 ;
  m_alloc_shared     = 0 ;
  m_team_shared      = 0 ;
  m_team_shared_end  = 0 ;
  m_team_shared_iter = 0 ;

  m_pool_rank     = 0 ;
  m_pool_size     = 0 ;
  m_pool_fan_size = 0 ;

  m_team_rank     = 0 ;
  m_team_size     = 0 ;
  m_team_fan_size = 0 ;

  m_league_rank = 0 ;
  m_league_end  = 0 ;
  m_league_size = 0 ;

  m_pool_state  = ThreadsExec::Terminating ;
  m_team_state  = ThreadsExec::Inactive ;
}

ThreadsExec::ThreadsExec()
  : m_pool_base(0)
  , m_team_base(0)
  , m_alloc_reduce(0)
  , m_alloc_shared(0)
  , m_team_shared(0)
  , m_team_shared_end(0)
  , m_team_shared_iter(0)

  , m_pool_rank(0)
  , m_pool_size(0)
  , m_pool_fan_size(0)

  , m_team_rank(0)
  , m_team_size(0)
  , m_team_fan_size(0)

  , m_league_rank(0)
  , m_league_end(0)
  , m_league_size(0)

  , m_pool_state( ThreadsExec::Terminating )
  , m_team_state( ThreadsExec::Inactive )
{
  if ( & s_threads_process == this ) {
    // Enables 'parallel_for' to execute on unitialized Threads device
    m_pool_state = ThreadsExec::Inactive ;
    m_pool_rank  = 0 ;
    m_pool_size  = 1 ;
  }
}

int ThreadsExec::get_thread_count()
{
  return s_threads_count ;
}

ThreadsExec * ThreadsExec::get_thread( const int init_thread_rank )
{
  ThreadsExec * const th =
    unsigned(init_thread_rank) < s_threads_count
    ? s_threads_exec[ s_threads_count - ( init_thread_rank + 1 ) ] : 0 ;

  if ( 0 == th || th->m_pool_rank != init_thread_rank ) {
    std::ostringstream msg ;
    msg << "Kokkos::Impl::ThreadsExec::get_thread ERROR : "
        << "thread " << init_thread_rank << " of " << s_threads_count ;
    if ( 0 == th ) {
      msg << " does not exist" ;
    }
    else {
      msg << " has wrong thread_rank " << th->m_pool_rank ;
    }
    Kokkos::Impl::throw_runtime_exception( msg.str() );
  }

  return th ;
}

//----------------------------------------------------------------------------

void ThreadsExec::execute_get_binding( ThreadsExec & exec , const void * )
{
  s_threads_coord[ exec.m_pool_rank ] = Kokkos::hwloc::get_this_thread_coordinate();
}

void ThreadsExec::execute_sleep( ThreadsExec & exec , const void * )
{
  ThreadsExec::global_lock();
  ThreadsExec::global_unlock();

  const int n = exec.m_pool_fan_size ;
  const int rank_rev = exec.m_pool_size - ( exec.m_pool_rank + 1 );

  for ( int i = 0 ; i < n ; ++i ) {
    Impl::spinwait( exec.m_pool_base[ rank_rev + (1<<i) ]->m_pool_state , ThreadsExec::Active );
  }

  exec.m_pool_state = ThreadsExec::Inactive ;
}

void ThreadsExec::execute_reduce_resize( ThreadsExec & exec , const void * )
{
  if ( exec.m_alloc_reduce ) {
    HostSpace::decrement( exec.m_alloc_reduce );
    exec.m_alloc_reduce = 0 ;
  }

  if ( s_current_reduce_size ) {

    exec.m_alloc_reduce =
      HostSpace::allocate( "reduce_scratch_space" , typeid(unsigned char) , 1 , s_current_reduce_size );

    // Guaranteed multiple of 'unsigned'

    unsigned * ptr = (unsigned *)( exec.m_alloc_reduce );
    unsigned * const end = ptr + s_current_reduce_size / sizeof(unsigned);

    // touch on this thread
    while ( ptr < end ) *ptr++ = 0 ;
  }
}

void ThreadsExec::execute_shared_resize( ThreadsExec & exec , const void * )
{
  // First thread pinned to a core allocates shared memory
  const int rank_rev = exec.m_pool_size - ( exec.m_pool_rank + 1 );

  if ( ! ( rank_rev % s_threads_per_core ) ) {

    if ( exec.m_alloc_shared ) {
      HostSpace::decrement( exec.m_alloc_shared );
      exec.m_alloc_shared = 0 ;
    }

    if ( s_current_shared_size ) {
      exec.m_alloc_shared =
        HostSpace::allocate( "shared_scratch_space" , typeid(unsigned char) , 1 , s_current_shared_size );

      // Guaranteed multiple of 'unsigned'

      unsigned * ptr = (unsigned *)( exec.m_alloc_shared );
      unsigned * const end = ptr + s_current_shared_size / sizeof(unsigned);

      // touch on this thread
      while ( ptr < end ) *ptr++ = 0 ;
    }
  }
  else {
    exec.m_alloc_shared = 0 ;
  }
}

void * ThreadsExec::get_shmem( const int size )
{
  // m_team_shared_iter is in bytes, convert to integer offsets
  const int offset = m_team_shared_iter >> power_of_two<sizeof(int)>::value ;

  m_team_shared_iter += size ;

  if ( m_team_shared_end < m_team_shared_iter ) {
    Kokkos::Impl::throw_runtime_exception( std::string("ThreadsExec::get_shmem FAILED : exceeded shared memory size" ) );
  }

  return ((int*)m_team_shared) + offset ;
}

}
}

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void ThreadsExec::verify_is_process( const std::string & name , const bool initialized )
{
  if ( ! is_process() ) {
    std::string msg( name );
    msg.append( " FAILED : Called by a worker thread, can only be called by the master process." );
    Kokkos::Impl::throw_runtime_exception( msg );
  }

  if ( initialized && 0 == s_threads_count ) {
    std::string msg( name );
    msg.append( " FAILED : Threads not initialized." );
    Kokkos::Impl::throw_runtime_exception( msg );
  }
}

int ThreadsExec::in_parallel()
{
  // A thread function is in execution and
  // the function argument is not the special threads process argument and
  // the master process is a worker or is not the master process.
  return s_current_function &&
         ( & s_threads_process != s_current_function_arg ) &&
         ( s_threads_process.m_pool_base || ! is_process() );
}

// Wait for root thread to become inactive
void ThreadsExec::fence()
{
  if ( s_threads_count ) {
    // Wait for the root thread to complete:
    Impl::spinwait( s_threads_exec[0]->m_pool_state , ThreadsExec::Active );

    if ( s_exception_msg.size() ) {
      Kokkos::Impl::throw_runtime_exception( s_exception_msg );
    }
  }

  s_current_function     = 0 ;
  s_current_function_arg = 0 ;
  s_current_team_size    = 0 ;
  s_current_team_alloc   = 0 ;
  s_current_league_size  = 0 ;
}

/** \brief  Begin execution of the asynchronous functor */
void ThreadsExec::start( void (*func)( ThreadsExec & , const void * ) , const void * arg ,
                         int work_league_size ,
                         int work_team_size )
{
  verify_is_process("ThreadsExec::start" , false );

  if ( s_current_function || s_current_function_arg ) {
    Kokkos::Impl::throw_runtime_exception( std::string( "ThreadsExec::start() FAILED : already executing" ) );
  }

  s_exception_msg.clear();

  s_current_function     = func ;
  s_current_function_arg = arg ;
  s_current_team_size    = std::min( s_threads_per_numa ,
                                     unsigned( work_team_size ? work_team_size : s_threads_per_team_default ) );
  s_current_team_alloc   = s_threads_per_core * ( ( s_current_team_size + s_threads_per_core - 1 ) / s_threads_per_core );
  s_current_league_size  = work_league_size ;

  // Activate threads:
  for ( int i = s_threads_count ; 0 < i-- ; ) {
    s_threads_exec[i]->m_pool_state = ThreadsExec::Active ;
  }

  if ( s_threads_process.m_pool_size ) {
    // Master process is the root thread, run it:
    s_threads_process.set_team_relations();
    (*func)( s_threads_process , arg );
    s_threads_process.m_pool_state = ThreadsExec::Inactive ;
  }
}

//----------------------------------------------------------------------------

bool ThreadsExec::sleep()
{
  verify_is_process("ThreadsExec::sleep", true );

  if ( & execute_sleep == s_current_function ) return false ;

  fence();

  ThreadsExec::global_lock();

  s_exception_msg.clear();

  s_current_function = & execute_sleep ;

  // Activate threads:
  for ( unsigned i = s_threads_count ; 0 < i ; ) {
    s_threads_exec[--i]->m_pool_state = ThreadsExec::Active ;
  }

  return true ;
}

bool ThreadsExec::wake()
{
  verify_is_process("ThreadsExec::wake", true );

  if ( & execute_sleep != s_current_function ) return false ;

  ThreadsExec::global_unlock();

  if ( s_threads_process.m_pool_base ) {
    execute_sleep( s_threads_process , 0 );
    s_threads_process.m_pool_state = ThreadsExec::Inactive ;
  }

  fence();

  return true ;
}

//----------------------------------------------------------------------------

void ThreadsExec::execute_serial( void (*func)( ThreadsExec & , const void * ) )
{
  s_exception_msg.clear();

  s_current_function = func ;
  s_current_function_arg = & s_threads_process ;

  const unsigned begin = s_threads_process.m_pool_base ? 1 : 0 ;

  for ( unsigned i = s_threads_count ; begin < i ; ) {
    ThreadsExec & th = * s_threads_exec[ --i ];

    th.m_pool_state = ThreadsExec::Active ;

    wait_yield( th.m_pool_state , ThreadsExec::Active );
  }

  if ( s_threads_process.m_pool_base ) {
    s_threads_process.m_pool_state = ThreadsExec::Active ;
    (*func)( s_threads_process , 0 );
    s_threads_process.m_pool_state = ThreadsExec::Inactive ;
  }

  s_current_function_arg = 0 ;
  s_current_function = 0 ;
}

//----------------------------------------------------------------------------

void * ThreadsExec::root_reduce_scratch()
{
  return s_threads_process.reduce_base();
}

void ThreadsExec::resize_reduce_scratch( size_t size )
{
  fence();

  if ( size ) { size += REDUCE_TEAM_BASE ; }

  const size_t rem = size % Kokkos::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Kokkos::Impl::MEMORY_ALIGNMENT - rem ;

  if ( ( s_current_reduce_size < size ) ||
       ( 0 == size && s_current_reduce_size ) ) {

    verify_is_process( "ThreadsExec::resize_reduce_scratch" , true );

    s_current_reduce_size = size ;

    execute_serial( & execute_reduce_resize );

    s_threads_process.m_alloc_reduce = s_threads_exec[0]->m_alloc_reduce ;
  }
}

void ThreadsExec::resize_shared_scratch( size_t size )
{
  fence();

  const size_t rem = size % Kokkos::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Kokkos::Impl::MEMORY_ALIGNMENT - rem ;

  if ( s_current_shared_size < size || ( 0 == size && s_current_shared_size ) ) {

    verify_is_process( "ThreadsExec::resize_shared_scratch" , true );

    s_current_shared_size = size ;

    execute_serial( & execute_shared_resize );
  }
}

//----------------------------------------------------------------------------

void ThreadsExec::print_configuration( std::ostream & s , const bool detail )
{
  verify_is_process("ThreadsExec::print_configuration",false);

  fence();

  const unsigned numa_count       = Kokkos::hwloc::get_available_numa_count();
  const unsigned cores_per_numa   = Kokkos::hwloc::get_available_cores_per_numa();
  const unsigned threads_per_core = Kokkos::hwloc::get_available_threads_per_core();

#if defined( KOKKOS_HAVE_HWLOC )
  s << "macro  KOKKOS_HAVE_HWLOC   : defined" << std::endl ;
#endif
#if defined( KOKKOS_HAVE_PTHREAD )
  s << "macro  KOKKOS_HAVE_PTHREAD : defined" << std::endl ;
#endif

  s << "Kokkos::Threads hwloc[" << numa_count << "x" << cores_per_numa << "x" << threads_per_core << "]" ;

  if ( s_threads_count ) {
    s << " threads[" << s_threads_count << "]"
      << " threads_per_numa[" << s_threads_per_numa << "]"
      << " threads_per_core[" << s_threads_per_core << "]"
      ;
    if ( 0 == s_threads_process.m_pool_base ) { s << " Asynchronous" ; }
    s << " ReduceScratch[" << s_current_reduce_size << "]"
      << " SharedScratch[" << s_current_shared_size << "]" ;
    s << std::endl ;

    if ( detail || 1 ) {

      execute_serial( & execute_get_binding );

      for ( unsigned i = 0 ; i < s_threads_count ; ++i ) {
        ThreadsExec * const th = s_threads_exec[i] ;
        s << "  Thread hwloc("
          << s_threads_coord[i].first << ","
          << s_threads_coord[i].second << ")" ;

        s_threads_coord[i].first  = ~0u ;
        s_threads_coord[i].second = ~0u ;

        if ( th ) {
          const int rank_rev = th->m_pool_size - ( th->m_pool_rank + 1 );

          s << " rank(" << th->m_pool_rank << ")" ;

          if ( th->m_pool_fan_size ) {
            s << " Fan{" ;
            for ( int j = 0 ; j < th->m_pool_fan_size ; ++j ) {
              s << " " << th->m_pool_base[rank_rev+(1<<j)]->m_pool_rank ;
            }
            s << " }" ;
          }

          if ( th->m_team_base && th->m_team_size ) {
            s << " Team[ " << th->m_team_base[0]->m_pool_rank
              << " .. " << th->m_team_base[ th->m_team_size - 1 ]->m_pool_rank
              << " ]" ;
          }

          if ( th == & s_threads_process ) {
            s << " is_process" ;
          }
        }
        s << std::endl ;
      }
    }
  }
  else {
    s << " not initialized" << std::endl ;
  }
}

//----------------------------------------------------------------------------

int ThreadsExec::league_max()
{ return std::numeric_limits<int>::max(); }

int ThreadsExec::team_max()
{ return s_threads_per_numa ; }

//----------------------------------------------------------------------------

int ThreadsExec::is_initialized()
{ return 0 != s_threads_exec[0] ; }

void ThreadsExec::initialize( 
  const std::pair<unsigned,unsigned> team_topology ,
        std::pair<unsigned,unsigned> use_core_topology )
{
  static const Sentinel sentinel ;

  verify_is_process("ThreadsExec::initialize",false);

  std::ostringstream msg ;

  msg << "Kokkos::Threads::initialize("
      << " team_topology(" << team_topology.first << "," << team_topology.second << ")"
      << ", use_core_topology(" << use_core_topology.first << "," << use_core_topology.second << ")"
      << " )" ;

  if ( s_threads_count ) {
    msg << " FAILED : Already initialized" ;
    Kokkos::Impl::throw_runtime_exception( msg.str() );
  }

  const unsigned thread_count = team_topology.first * team_topology.second ;

  if ( 0 == thread_count ) {
    msg << " FAILED : zero thread count" ;
    Kokkos::Impl::throw_runtime_exception( msg.str() );
  }
  //------------------------------------
  // Query hardware topology and capacity, if available.

  const bool hwloc_avail = Kokkos::hwloc::available();

  const std::pair<unsigned,unsigned>
    hwloc_core_topo( Kokkos::hwloc::get_available_numa_count() ,
                     Kokkos::hwloc::get_available_cores_per_numa() );

  std::pair<unsigned,unsigned> master_coord = Kokkos::hwloc::get_this_thread_coordinate();
  bool                         asynchronous = false ;

  if ( hwloc_avail && 1 < thread_count ) {

    if ( 0 == use_core_topology.first || 0 == use_core_topology.second ) {
      use_core_topology = Kokkos::hwloc::use_core_topology( thread_count );
    }

    if ( use_core_topology.first < hwloc_core_topo.first ) {
      // Can omit a (NUMA) group of cores and execute work asynchronously
      // on the other groups.

      Kokkos::hwloc::thread_mapping( team_topology , use_core_topology , hwloc_core_topo , s_threads_coord );

      // Don't use master thread's first core coordinate (NUMA region).
      // Originally mapped:
      //   begin = hwloc_core_topo.first - use_core_topology.first ;
      //   end   = hwloc_core_topo.first ;
      // So can decrement.

      for ( unsigned i = 0 ; i < thread_count ; ++i ) {
        if ( s_threads_coord[i].first <= master_coord.first ) {
          --( s_threads_coord[i].first );
        }
      }

      asynchronous = true ;
    }
    else if ( use_core_topology.second < hwloc_core_topo.second ) {
      // Can omit a core from each group and execute work asynchronously

      Kokkos::hwloc::thread_mapping( team_topology , use_core_topology , hwloc_core_topo , s_threads_coord );

      // Threads' coordinates are in the range
      //   0 <= numa_begin = hwloc_core_topo.first - use_core_topology.first
      //   1 <= numa_end   = hwloc_core_topo.first
      //   1 <= core_begin = hwloc_core_topo.second - use_core_topology.second
      //   1 <= core_end   = hwloc_core_topo.second
      //
      //   range: ( [numa_begin,numa_end) , [core_begin,core_end) )
      //
      // Force master thread onto the highest rank unused core of its current numa region.
      //
      master_coord.second = ( hwloc_core_topo.second - use_core_topology.second ) - 1 ;

      asynchronous = true ;
    }
    else {
      // Spawn threads with root thread on the master process' core

      Kokkos::hwloc::thread_mapping( team_topology , use_core_topology , hwloc_core_topo , master_coord , s_threads_coord );

      s_threads_coord[0] = std::pair<unsigned,unsigned>( ~0u , ~0u );
    }
  }
  else {
    // Not applicable so ignore it.
    use_core_topology.first  = 0 ;
    use_core_topology.second = 0 ;
  }

  //------------------------------------
  // Spawn threads

  {
    const unsigned thread_spawn_begin  = asynchronous ? 0 : 1 ;
    unsigned       thread_spawn_failed = 0 ;

    s_threads_count    = thread_count ;
    s_current_function = & execute_function_noop ; // Initialization work function

    // If not fully utilizing the capacity then spawn threads for asynchronous execution.

    for ( unsigned i = thread_spawn_begin ; i < thread_count ; ++i ) {

      s_threads_process.m_pool_state = ThreadsExec::Inactive ;

      // If hwloc available then spawned thread will choose its own rank,
      // otherwise specify the rank.
      s_current_function_arg = (void*)( hwloc_avail ? ~0u : i );

      // Spawn thread executing the 'driver()' function.
      // Wait until spawned thread has attempted to initialize.
      // If spawning and initialization is successfull then
      // an entry in 's_threads_exec' will be assigned.
      if ( ThreadsExec::spawn() ) {
        wait_yield( s_threads_process.m_pool_state , ThreadsExec::Inactive );
      }
    }

    // Wait for all spawned threads to deactivate before zeroing the function.

    for ( unsigned i = thread_spawn_begin ; i < thread_count ; ++i ) {
      // Try to protect against cache coherency failure by casting to volatile.
      ThreadsExec * const th = ((ThreadsExec * volatile *)s_threads_exec)[i] ;
      if ( th ) {
        wait_yield( th->m_pool_state , ThreadsExec::Active );
      }
      else {
        ++thread_spawn_failed ;
      }
    }

    s_current_function     = 0 ;
    s_current_function_arg = 0 ;

    if ( thread_spawn_failed ) {

      s_threads_count = 0 ;

      msg << " FAILED " << thread_spawn_failed << " attempts to spawn threads" ;

      Kokkos::Impl::throw_runtime_exception( msg.str() );
    }

    if ( 1 < thread_count ) { Kokkos::hwloc::bind_this_thread( master_coord ); }

    // Clear master thread data.
    // The master thread will be unused or initialized
    // as part of the thread pool.

    s_threads_process.m_pool_base = 0 ;
    s_threads_process.m_pool_rank = 0 ;
    s_threads_process.m_pool_size = 0 ;
    s_threads_process.m_pool_fan_size = 0 ;

    s_threads_process.m_team_base = 0 ;
    s_threads_process.m_team_rank = 0 ;
    s_threads_process.m_team_size = 0 ;

    s_threads_process.m_league_rank = 0 ;
    s_threads_process.m_league_end  = 0 ;
    s_threads_process.m_league_size = 0 ;

    s_threads_process.m_pool_state = ThreadsExec::Inactive ;

    if ( thread_spawn_begin ) {
      s_threads_exec[0] = & s_threads_process ; // Include the master thread in pool.
    }
  }

  //------------------------------------
  // Initialize team topology and fan-in/out relationships:

  s_threads_per_numa = use_core_topology.first ? thread_count / use_core_topology.first : thread_count ;
  s_threads_per_core = use_core_topology.first ? thread_count / ( use_core_topology.first * use_core_topology.second ) : 1 ;
  s_threads_per_team_default = team_topology.second ;

  // Set threads' pool size, rank, and fan-in counts.
  for ( unsigned i = 0 ; i < thread_count ; ++i ) {
    ThreadsExec & th = * s_threads_exec[i] ;

    th.m_pool_base     = s_threads_exec ;
    th.m_pool_rank     = thread_count - ( i + 1 ); // Reversed for scan-compatible reductions
    th.m_pool_size     = thread_count ;
    th.m_pool_fan_size = fan_size( th.m_pool_rank , th.m_pool_size );
  }

  // Initial allocations:
  ThreadsExec::resize_reduce_scratch( 4096 - REDUCE_TEAM_BASE );
  ThreadsExec::resize_shared_scratch( 4096 );
}

//----------------------------------------------------------------------------

void ThreadsExec::finalize()
{
  verify_is_process("ThreadsExec::finalize",false);

  fence();

  resize_reduce_scratch(0);
  resize_shared_scratch(0);

  const unsigned begin = s_threads_process.m_pool_base ? 1 : 0 ;

  for ( unsigned i = s_threads_count ; begin < i-- ; ) {

    if ( s_threads_exec[i] ) {

      s_threads_exec[i]->m_pool_state = ThreadsExec::Terminating ;

      wait_yield( s_threads_process.m_pool_state , ThreadsExec::Inactive );

      s_threads_process.m_pool_state = ThreadsExec::Inactive ;
    }
  }

  if ( s_threads_process.m_pool_base ) {
    ( & s_threads_process )->~ThreadsExec();
    s_threads_exec[0] = 0 ;
  }

  Kokkos::hwloc::unbind_this_thread();

  s_threads_count    = 0 ;
  s_threads_per_numa = 0 ;
  s_threads_per_core = 0 ;
  s_threads_per_team_default = 0 ;

  // Reset master thread to run solo.
  s_threads_process.m_pool_base = 0 ;
  s_threads_process.m_pool_rank = 0 ;
  s_threads_process.m_pool_size = 1 ;
  s_threads_process.m_pool_fan_size = 0 ;
  s_threads_process.m_pool_state = ThreadsExec::Inactive ;
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */


