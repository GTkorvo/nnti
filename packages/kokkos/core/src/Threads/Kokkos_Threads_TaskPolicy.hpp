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

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_THREADS_TASKPOLICY_HPP
#define KOKKOS_THREADS_TASKPOLICY_HPP


#include <Kokkos_Threads.hpp>
#include <Kokkos_TaskPolicy.hpp>

#if defined( KOKKOS_HAVE_PTHREAD )

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
namespace Impl {

/** \brief  Base class for all Kokkos::Threads tasks */
template<>
class TaskMember< Kokkos::Threads , void , void > {
public:

  typedef void         (* function_dealloc_type)( TaskMember * );
  typedef TaskMember * (* function_verify_type) ( TaskMember * );
  typedef void         (* function_serial_type) ( TaskMember * );
  typedef void         (* function_team_type)   ( TaskMember * , Kokkos::Impl::ThreadsExecTeamMember & );

private:

  function_dealloc_type  m_dealloc ;      ///< Deallocation
  function_verify_type   m_verify ;       ///< Result type verification
  function_team_type     m_team ;         ///< Apply function
  function_serial_type   m_serial ;       ///< Apply function
  TaskMember **          m_dep ;          ///< Dependences
  TaskMember *           m_wait ;         ///< Linked list of tasks waiting on this task
  TaskMember *           m_next ;         ///< Linked list of tasks waiting on a different task
  int                    m_dep_capacity ; ///< Capacity of dependences
  int                    m_dep_size ;     ///< Actual count of dependences
  int                    m_shmem_size ;
  int                    m_league_size ;
  int                    m_league_rank ;
  int                    m_league_end ;
  int                    m_ref_count ;    ///< Reference count
  int                    m_state ;        ///< State of the task

  // 7 pointers + 8 integers

  TaskMember( const TaskMember & ) = delete ;
  TaskMember & operator = ( const TaskMember & ) = delete ;

  static void * allocate( const unsigned arg_size );
  static void deallocate( void * );

  template< class DerivedTaskType >
  static
  void deallocate( TaskMember * t )
    {
      DerivedTaskType * ptr = static_cast< DerivedTaskType * >(t);
      ptr->~DerivedTaskType();
      deallocate( (void*) ptr );
    }

  static void execute_serial( TaskMember * const );
  static bool execute_ready_team_task( Kokkos::Impl::ThreadsExecTeamMember & );
  static void execute_ready_tasks_driver( Kokkos::Impl::ThreadsExec & , const void * );
  static void throw_error_verify_type();
  static int team_fixed_size();

protected:

  TaskMember()
    : m_dealloc(0)
    , m_verify(0)
    , m_team(0)
    , m_serial(0)
    , m_dep(0)
    , m_wait(0)
    , m_next(0)
    , m_dep_capacity(0)
    , m_dep_size(0)
    , m_shmem_size(0)
    , m_league_size(0)
    , m_league_rank(0)
    , m_league_end(0)
    , m_ref_count(0)
    , m_state(0)
    {}

public:

  ~TaskMember();

  template< typename ResultType >
  KOKKOS_FUNCTION static
  TaskMember * verify_type( TaskMember * t )
    {
      enum { check_type = ! Kokkos::Impl::is_same< ResultType , void >::value };

      if ( check_type && t != 0 ) {

        // Verify that t->m_verify is this function
        const function_verify_type self = & TaskMember::template verify_type< ResultType > ;

        if ( t->m_verify != self ) {
          t = 0 ;
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
          throw_error_verify_type();
#endif
        }
      }
      return t ;
    }

  

  template< class DerivedTaskType , class Tag >
  KOKKOS_FUNCTION static
  void serial_function( typename Kokkos::Impl::enable_if< ! Kokkos::Impl::is_same< typename DerivedTaskType::result_type , void >::value
                                                , TaskMember * >::type t )
    {
      typedef typename DerivedTaskType::functor_type  functor_type ;
      typedef typename DerivedTaskType::result_type   result_type ;

      DerivedTaskType & self = * static_cast< DerivedTaskType * >(t);

      Kokkos::Impl::FunctorApply< functor_type , Tag , result_type & >
        ::apply( (functor_type &) self , & self.m_result );
    }

  template< class DerivedTaskType , class Tag >
  KOKKOS_FUNCTION static
  void serial_function( typename Kokkos::Impl::enable_if< Kokkos::Impl::is_same< typename DerivedTaskType::result_type , void >::value
                                                , TaskMember * >::type t )
    {
      typedef typename DerivedTaskType::functor_type  functor_type ;

      DerivedTaskType & self = * static_cast< DerivedTaskType * >(t);

      Kokkos::Impl::FunctorApply< functor_type , Tag , void >::apply( (functor_type &) self );
    }

  template< class DerivedTaskType , class Tag >
  KOKKOS_FUNCTION static
  void team_function( typename Kokkos::Impl::enable_if<( Kokkos::Impl::is_same<Tag,void>::value ), TaskMember * >::type t
                    , Kokkos::Impl::ThreadsExecTeamMember & member
                    )
    {
      typedef typename DerivedTaskType::functor_type  functor_type ;
      DerivedTaskType & self = * static_cast< DerivedTaskType * >(t);
      self.functor_type::operator()( member );
    }

  //----------------------------------------
  /*  Inheritence Requirements on task types:
   *
   *    class DerivedTaskType
   *      : public TaskMember< Threads , DerivedType::value_type , FunctorType >
   *      { ... };
   *
   *    class TaskMember< Threads , DerivedType::value_type , FunctorType >
   *      : public TaskMember< Threads , DerivedType::value_type , void >
   *      , public Functor
   *      { ... };
   *
   *  If value_type != void
   *    class TaskMember< Threads , value_type , void >
   *      : public TaskMember< Threads , void , void >
   *
   *  Allocate space for DerivedTaskType followed by TaskMember*[ dependence_capacity ]
   *
   */

  /** \brief  Allocate and construct a task */
  template< class DerivedTaskType , class Tag >
  static
  TaskMember * create( const typename DerivedTaskType::functor_type &  arg_functor
                     , const function_team_type                        arg_team_function
                     , const function_serial_type                      arg_serial_function
                     , const unsigned                                  arg_dependence_capacity
                     , const unsigned                                  arg_league_size = 0
                     )
    {
      enum { padding_size = sizeof(DerivedTaskType) % sizeof(TaskMember*)
                          ? sizeof(TaskMember*) - sizeof(DerivedTaskType) % sizeof(TaskMember*) : 0 };
      enum { derived_size = sizeof(DerivedTaskType) + padding_size };

      DerivedTaskType * const task =
        new( allocate( derived_size + sizeof(TaskMember*) * arg_dependence_capacity ) )
          DerivedTaskType( arg_functor );

      task->TaskMember::m_dealloc      = & TaskMember::template deallocate< DerivedTaskType > ;
      task->TaskMember::m_verify       = & TaskMember::template verify_type< typename DerivedTaskType::value_type > ;
      task->TaskMember::m_team         = arg_team_function ;
      task->TaskMember::m_serial       = arg_serial_function ;
      task->TaskMember::m_dep          = (TaskMember**)( ((unsigned char *)task) + derived_size );
      task->TaskMember::m_dep_capacity = arg_dependence_capacity ;
      task->TaskMember::m_shmem_size   = Kokkos::Impl::FunctorTeamShmemSize< typename DerivedTaskType::functor_type >
                                           ::value( arg_functor , team_fixed_size() );
      task->TaskMember::m_league_rank  = arg_league_size ;
      task->TaskMember::m_league_rank  = arg_league_size - 1 ;
      task->TaskMember::m_league_end   = arg_league_size - 1 ;
      task->TaskMember::m_state        = TASK_STATE_CONSTRUCTING ;

      for ( unsigned i = 0 ; i < arg_dependence_capacity ; ++i ) task->TaskMember::m_dep[i] = 0 ;

      return static_cast< TaskMember * >( task );
    }

  void reschedule();
  void schedule();
  static void execute_ready_tasks();
  static void wait( const Future< void , Kokkos::Threads > & );

  //----------------------------------------

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  static
  void assign( TaskMember ** const lhs , TaskMember * const rhs );
#else
  KOKKOS_INLINE_FUNCTION static
  void assign( TaskMember ** const lhs , TaskMember * const rhs ) {}
#endif

  TaskMember * get_dependence( int i ) const ;

  KOKKOS_INLINE_FUNCTION
  int get_dependence() const
    { return m_dep_size ; }

  void clear_dependence();
  void add_dependence( TaskMember * before );

  //----------------------------------------

  typedef FutureValueTypeIsVoidError get_result_type ;

  KOKKOS_INLINE_FUNCTION
  get_result_type get() const { return get_result_type() ; }

  KOKKOS_INLINE_FUNCTION
  Kokkos::Experimental::TaskState get_state() const { return Kokkos::Experimental::TaskState( m_state ); }

};

/** \brief  A Future< Kokkos::Threads , ResultType > will cast
 *          from  TaskMember< Kokkos::Threads , void , void >
 *          to    TaskMember< Kokkos::Threads , ResultType , void >
 *          to query the result.
 */
template< class ResultType >
class TaskMember< Kokkos::Threads , ResultType , void >
  : public TaskMember< Kokkos::Threads , void , void >
{
public:

  typedef ResultType result_type ;

  result_type  m_result ;

  typedef const result_type & get_result_type ;

  KOKKOS_INLINE_FUNCTION
  get_result_type get() const { return m_result ; }

  inline
  TaskMember() : TaskMember< Kokkos::Threads , void , void >(), m_result() {}

  TaskMember( const TaskMember & ) = delete ;
  TaskMember & operator = ( const TaskMember & ) = delete ;
};

/** \brief  Callback functions will cast
 *          from  TaskMember< Kokkos::Threads , void , void >
 *          to    TaskMember< Kokkos::Threads , ResultType , FunctorType >
 *          to execute work functions.
 */
template< class ResultType , class FunctorType >
class TaskMember< Kokkos::Threads , ResultType , FunctorType >
  : public TaskMember< Kokkos::Threads , ResultType , void >
  , public FunctorType
{
public:
  typedef ResultType   result_type ;
  typedef FunctorType  functor_type ;

  inline
  TaskMember( const functor_type & arg_functor )
    : TaskMember< Kokkos::Threads , ResultType , void >()
    , functor_type( arg_functor )
    {}
};

} /* namespace Impl */
} /* namespace Experimental */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {

template<>
class TaskPolicy< Kokkos::Threads >
{
public:

  typedef Kokkos::Threads                      execution_space ;
  typedef TaskPolicy                           execution_policy ;
  typedef Kokkos::Impl::ThreadsExecTeamMember  member_type ;

private:

  typedef Impl::TaskMember< Kokkos::Threads , void , void >  task_root_type ;

  int m_default_dependence_capacity ;

  TaskPolicy & operator = ( const TaskPolicy & ) = delete ;

  template< class FunctorType >
  static inline
  const task_root_type * get_task_root( const FunctorType * f )
    {
      typedef Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType > task_type ;
      return static_cast< const task_root_type * >( static_cast< const task_type * >(f) );
    }

  template< class FunctorType >
  static inline
  task_root_type * get_task_root( FunctorType * f )
    {
      typedef Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType > task_type ;
      return static_cast< task_root_type * >( static_cast< task_type * >(f) );
    }

public:

  TaskPolicy() : m_default_dependence_capacity(4) {}

  TaskPolicy( const TaskPolicy & rhs ) : m_default_dependence_capacity( rhs.m_default_dependence_capacity ) {}
  
  KOKKOS_INLINE_FUNCTION
  explicit
  TaskPolicy( const unsigned arg_default_dependence_capacity )
    : m_default_dependence_capacity( arg_default_dependence_capacity ) {}

  KOKKOS_INLINE_FUNCTION
  TaskPolicy( const TaskPolicy &
            , const unsigned arg_default_dependence_capacity )
    : m_default_dependence_capacity( arg_default_dependence_capacity ) {}

  // Create serial-thread task

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  Future< typename FunctorType::value_type , execution_space >
  create( const FunctorType & functor
        , const unsigned dependence_capacity = ~0u ) const
    {
      typedef typename FunctorType::value_type  value_type ;
      typedef Impl::TaskMember< execution_space , value_type , FunctorType >  task_type ;

      return Future< value_type , execution_space >(
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        task_root_type::create< task_type , void >
          ( functor
          , task_root_type::function_team_type(0)
          , & task_root_type::template serial_function< task_type , void >
          , ( ~0u == dependence_capacity ? m_default_dependence_capacity : dependence_capacity )
          , 0
          )
#endif
        );
    }

#if 0

  // Create parallel foreach task

  template< class PolicyType , class FunctorType >
  KOKKOS_INLINE_FUNCTION
  Future< typename FunctorType::value_type , execution_space >
  create_foreach( const PolicyType  & policy
                , const FunctorType & functor
                , const unsigned      dependence_capacity = ~0u ) const
    {
      typedef typename FunctorType::value_type  value_type ;
      typedef typename PolicyType::work_tag     work_tag ;

      typedef Impl::TaskMember< execution_space , value_type , FunctorType >  task_type ;

      return Future< value_type , execution_space >(
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        task_root_type::create< task_type , work_tag >
          ( functor
          , & task_root_type::template team_function<   task_type , work_tag >
          , & task_root_type::template serial_function< task_type , work_tag > 
          , ( ~0u == dependence_capacity ? m_default_dependence_capacity : dependence_capacity )
          , policy.league_size()
          )
#endif
       );
    }

  // Create parallel reduce task

  template< class PolicyType , class FunctorType >
  KOKKOS_INLINE_FUNCTION
  Future< typename FunctorType::value_type , execution_space >
  create_reduce( const PolicyType  & policy
               , const FunctorType & functor
               , const unsigned      dependence_capacity = ~0u ) const
    {
      typedef typename FunctorType::value_type  value_type ;
      typedef typename PolicyType::work_tag     work_tag ;

      typedef Impl::TaskMember< execution_space , value_type , FunctorType >  task_type ;

      return Future< value_type , execution_space >(
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        task_root_type::create< task_type , work_tag >
          ( functor
          , & task_root_type::template team_function<   task_type , work_tag >
          , & task_root_type::template serial_function< task_type , work_tag > 
          , ( ~0u == dependence_capacity ? m_default_dependence_capacity : dependence_capacity )
          , policy.league_size()
          )
#endif
        );
    }

#endif


  template< class A1 , class A2 , class A3 , class A4 >
  KOKKOS_INLINE_FUNCTION
  void add_dependence( const Future<A1,A2> & after
                     , const Future<A3,A4> & before
                     , typename Kokkos::Impl::enable_if
                        < Kokkos::Impl::is_same< typename Future<A1,A2>::execution_space , execution_space >::value
                          &&
                          Kokkos::Impl::is_same< typename Future<A3,A4>::execution_space , execution_space >::value
                        >::type * = 0
                      ) const
    {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      after.m_task->add_dependence( before.m_task );
#endif
    }

  template< class FunctorType , class A3 , class A4 >
  KOKKOS_INLINE_FUNCTION
  void add_dependence( FunctorType * task_functor
                     , const Future<A3,A4> & before
                     , typename Kokkos::Impl::enable_if
                        < Kokkos::Impl::is_same< typename Future<A3,A4>::execution_space , execution_space >::value
                        >::type * = 0
                      ) const
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { get_task_root(task_functor)->add_dependence( before.m_task ); }
#else
    {}
#endif


  template< class ValueType >
  const Future< ValueType , execution_space > &
    spawn( const Future< ValueType , execution_space > & f ) const
      {
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        f.m_task->schedule();
#endif
        return f ;
      }

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  void respawn( FunctorType * task_functor ) const
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { get_task_root(task_functor)->reschedule(); }
#else
    {}
#endif

  //----------------------------------------
  // Functions for an executing task functor to query dependences,
  // set new dependences, and respawn itself.

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  Future< void , execution_space >
  get_dependence( const FunctorType * task_functor , int i ) const
    {
      return Future<void,execution_space>(
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
        get_task_root(task_functor)->get_dependence(i)
#endif
        );
    }

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  int get_dependence( const FunctorType * task_functor ) const
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { return get_task_root(task_functor)->get_dependence(); }
#else
    { return 0 ; }
#endif

  template< class FunctorType >
  KOKKOS_INLINE_FUNCTION
  void clear_dependence( FunctorType * task_functor ) const
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { get_task_root(task_functor)->clear_dependence(); }
#else
    {}
#endif

};

inline
void wait( TaskPolicy< Kokkos::Threads > & )
{
  Impl::TaskMember< Kokkos::Threads , void , void >::execute_ready_tasks();
}

inline
void wait( const Future< void , Kokkos::Threads > & future )
{ Impl::TaskMember< Kokkos::Threads , void , void >::wait( future ); }


} /* namespace Experimental */
} /* namespace Kokkos */

#endif /* #if defined( KOKKOS_HAVE_PTHREAD ) */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_THREADS_TASKPOLICY_HPP */


