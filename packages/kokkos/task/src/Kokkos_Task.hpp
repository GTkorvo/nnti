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

#ifndef KOKKOS_TASK_HPP
#define KOKKOS_TASK_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_StaticAssert.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {

template< class DeviceType , class ResultType >
class task_serial ;

template< class DeviceType , class WorkType = size_t >
class task_for ;

template< class DeviceType , class ResultType , class WorkType = size_t >
class task_reduce ;

template< class DeviceType , class ResultType , class WorkType = size_t >
class task_scan ;

} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

struct TaskResultTypeIsVoid {};

template< class PatternType
        , class Op1Type = void 
        , class Op2Type = void 
        , class Op3Type = void 
        , class Op4Type = void
        >
class FunctorTraits ;

//----------------------------------------------------------------------------
// Functor must have a pattern_type and must not be the pattern_type

template< class FunctorType >
class FunctorTraits< FunctorType , void , void , void , void > {
public:
  typedef FunctorType                         functor_type ;
  typedef typename FunctorType::pattern_type  pattern_type ;
  typedef typename pattern_type::device_type  device_type ;
  typedef typename pattern_type::result_type  result_type ;
  typedef typename pattern_type::future_type  future_type ;
  typedef typename pattern_type::work_type    work_type ;
};

//----------------------------------------------------------------------------

template< class DeviceType , class Op1Type >
class FunctorTraits< task_serial< DeviceType , void > , Op1Type , void , void , void > {
public:
  typedef task_serial< DeviceType , void >    pattern_type ;
  typedef typename pattern_type::device_type  device_type ;
  typedef typename pattern_type::result_type  result_type ;
  typedef typename pattern_type::future_type  future_type ;
  typedef typename pattern_type::work_type    work_type ;

  class functor_type : public pattern_type {
  private:
    Op1Type m_op1 ;
  public:

    KOKKOS_INLINE_FUNCTION
    void apply() { m_op1(); }

    KOKKOS_INLINE_FUNCTION
    functor_type( const Op1Type & arg_op1 )
      : task_serial< DeviceType , void >()
      , m_op1( arg_op1 )
      {}
  };
};

template< class DeviceType , class ResultType , class Op1Type >
class FunctorTraits< task_serial< DeviceType , ResultType > , Op1Type , void , void , void > {
public:
  typedef task_serial< DeviceType , ResultType >  pattern_type ;
  typedef typename pattern_type::device_type      device_type ;
  typedef typename pattern_type::result_type      result_type ;
  typedef typename pattern_type::future_type      future_type ;
  typedef typename pattern_type::work_type        work_type ;

  class functor_type : public pattern_type {
  private:
    Op1Type m_op1 ;
  public:

    KOKKOS_INLINE_FUNCTION
    void apply( ResultType & result ) { m_op1(result); }

    KOKKOS_INLINE_FUNCTION
    functor_type( const Op1Type & arg_op1 )
      : task_serial< DeviceType , ResultType >()
      , m_op1( arg_op1 )
      {}
  };
};

//----------------------------------------------------------------------------

template< class DeviceType , class ResultType , class WorkType , class Op1Type >
class FunctorTraits< task_reduce< DeviceType , ResultType , WorkType >
                   , Op1Type , void , void , void >
{
  typedef task_reduce< DeviceType , ResultType , WorkType >  pattern_type ;
  typedef typename pattern_type::device_type  device_type ;
  typedef typename pattern_type::result_type  result_type ;
  typedef typename pattern_type::future_type  future_type ;
  typedef typename pattern_type::work_type    work_type ;

  class functor_type : public pattern_type {
  private:
    Op1Type  m_op1 ;
  public:

    KOKKOS_INLINE_FUNCTION
    void operator()( const work_type & work , result_type & update ) const
      { m_op1( work , update ); }

    KOKKOS_INLINE_FUNCTION
    functor_type( const WorkType & arg_work , const Op1Type & arg_op1 )
      : pattern_type( arg_work )
      , m_op1( arg_op1 )
      {}
  };
};

//----------------------------------------------------------------------------

} /* namespace Imple */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class DeviceType , class ResultType = void >
class Task ;

} /* namespace Imple */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {

template< class DeviceType , class ResultType = void >
class Future {
public:

  KOKKOS_INLINE_FUNCTION
  Future() : m_task(0) {}

  KOKKOS_INLINE_FUNCTION
  Future( const Future & rhs )
    : m_task( rhs.m_task )
    { Impl::Task<DeviceType>::increment( m_task ); }

  KOKKOS_INLINE_FUNCTION
  Future & operator = ( const Future & rhs )
    {
      Impl::Task<DeviceType>::decrement( m_task );
      m_task = rhs.m_task ;
      Impl::Task<DeviceType>::increment( m_task );
      return *this ;
    }

  template< class T >
  KOKKOS_INLINE_FUNCTION
  Future( const Future<DeviceType,T> & rhs )
    : m_task( Impl::Task<DeviceType,ResultType>::verify_type( rhs.m_task ) )
    { Impl::Task<DeviceType>::increment( m_task ); }

  template< class T >
  KOKKOS_INLINE_FUNCTION
  Future & operator = ( const Future<DeviceType,T> & rhs )
    {
      Impl::Task<DeviceType>::decrement( m_task );
      m_task = Impl::Task<DeviceType,ResultType>::verify_type( rhs.m_task );
      Impl::Task<DeviceType>::increment( m_task );
      return *this ;
    }

  typedef typename Impl::Task<DeviceType,ResultType>::get_result_type get_result_type ;

  KOKKOS_INLINE_FUNCTION
  get_result_type get() const 
    { return static_cast< Impl::Task<DeviceType,ResultType> *>( m_task )->get(); }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION explicit
  Future( Impl::Task<DeviceType> * t )
    : m_task( Impl::Task<DeviceType,ResultType>::verify_type(t) )
    { Impl::Task<DeviceType>::increment( m_task ); }

  KOKKOS_INLINE_FUNCTION
  Future( const Impl::Task<DeviceType> * const t , int i )
    : m_task( t ? Impl::Task<DeviceType,ResultType>::verify_type( t->task_dependence(i) )
                : (Impl::Task<DeviceType,ResultType>*) 0 )
    { Impl::Task<DeviceType>::increment( m_task ); }

  static inline
  void wait( Future const & f ) { Impl::Task<DeviceType>::task_wait( f.m_task ); }

private:

  Impl::Task<DeviceType> * m_task ;

  template< class , class > friend class Future ;
  template< class , class > friend class Impl::Task ;
};

} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {

template< class DeviceType , class ResultType >
inline
void wait( Future<DeviceType,ResultType> const & f )
{ Future<DeviceType,ResultType>::wait( f ); }

//----------------------------------------------------------------------------
// One argument spawn:

/**\brief  Spawn a Functor derived from a pattern and return a Future */
template< class FunctorType >
KOKKOS_INLINE_FUNCTION
Future< typename Impl::FunctorTraits< FunctorType >::device_type >
spawn( FunctorType const & functor )
{
  typedef Impl::FunctorTraits< FunctorType > traits ;
  typedef typename traits::device_type  device_type ;
  typedef Impl::Task< device_type >     task_type ;

  return Future< device_type >( task_type::template task_spawn< traits >( functor ) );
}

/**\brief  Spawn a task_serial defined as a composition and return a Future */
template< class PatternType , class Closure1Type >
KOKKOS_INLINE_FUNCTION
Future< typename Impl::FunctorTraits< PatternType , Closure1Type >::device_type >
spawn( Closure1Type const & closure1 )
{
  typedef Impl::FunctorTraits< PatternType , Closure1Type > traits ;
  typedef typename traits::device_type  device_type ;
  typedef Impl::Task< device_type >     task_type ;

  return Future< device_type >( task_type::template task_spawn< traits >( closure1 ) );
}

//----------------------------------------------------------------------------
// Two argument spawn:

/**\brief  Spawn a data parallel task defined from a composition and return a Future */
template< class PatternType , class Closure1Type >
KOKKOS_INLINE_FUNCTION
Future< typename Impl::FunctorTraits< PatternType , Closure1Type >::device_type >
spawn( typename PatternType::work_type const & work
     , Closure1Type const & closure1
     )
{
  typedef Impl::FunctorTraits< PatternType , Closure1Type > traits ;
  typedef typename traits::device_type  device_type ;
  typedef Impl::Task< device_type >     task_type ;

  return Future< device_type >( task_type::template task_spawn< traits >( work , closure1 ) );
}

//----------------------------------------------------------------------------
// Three argument spawn:

/**\brief  Spawn a task with dependences on other tasks and return a Future */
template< class FunctorType , class FutureType >
KOKKOS_INLINE_FUNCTION
Future< typename Impl::FunctorTraits< FunctorType >::device_type >
spawn( FunctorType const & functor 
     , FutureType const * const ibegin
     , FutureType const * const iend
     )
{
  typedef Impl::FunctorTraits< FunctorType > traits ;
  typedef typename traits::device_type  device_type ;
  typedef Impl::Task< device_type >     task_type ;

  return Future< device_type >( task_type::template task_spawn< traits >( functor , ibegin , iend ) );
}

/**\brief  Spawn a data parallel task defined as a composition and return a Future */
template< class PatternType , class Closure1Type , class Closure2Type >
KOKKOS_INLINE_FUNCTION
Future< typename Impl::FunctorTraits< PatternType , Closure1Type , Closure2Type >::device_type >
spawn( typename PatternType::work_type const & work 
     , Closure1Type const & closure1
     , Closure2Type const & closure2
     )
{
  typedef Impl::FunctorTraits< PatternType , Closure1Type , Closure2Type > traits ;
  typedef typename traits::device_type  device_type ;
  typedef Impl::Task< device_type >     task_type ;

  return Future< device_type >( task_type::template task_spawn< traits >( work , closure1 , closure2 ) );
}

/**\brief  Spawn a serial task defined as a composition and return a Future */
template< class PatternType , class Closure1Type , class FutureType >
KOKKOS_INLINE_FUNCTION
Future< typename Impl::FunctorTraits< PatternType , Closure1Type >::device_type >
spawn( Closure1Type const & closure1
     , FutureType const * const ibegin
     , FutureType const * const iend
     )
{
  typedef Impl::FunctorTraits< PatternType , Closure1Type > traits ;
  typedef typename traits::device_type  device_type ;
  typedef Impl::Task< device_type >     task_type ;

  return Future< device_type >( task_type::template task_spawn< traits >( closure1 , ibegin , iend ) );
}

//----------------------------------------------------------------------------
// Four argument spawn:

/**\brief  Spawn a task defined as a composition and return a Future */
template< class PatternType , class Closure1Type , class Closure2Type , class Closure3Type >
KOKKOS_INLINE_FUNCTION
Future< typename Impl::FunctorTraits< PatternType , Closure1Type , Closure2Type , Closure3Type >::device_type >
spawn( typename PatternType::work_type const & work
     , Closure1Type const & closure1
     , Closure2Type const & closure2
     , Closure3Type const & closure3
     )
{
  typedef Impl::FunctorTraits< PatternType , Closure1Type , Closure2Type > traits ;
  typedef typename traits::device_type  device_type ;
  typedef Impl::Task< device_type >     task_type ;

  return Future< device_type >( task_type::template task_spawn< traits >( work , closure1 , closure2 , closure3 ) );
}

/**\brief  Spawn a task defined as a composition and return a Future */
template< class PatternType , class Closure1Type , class FutureType >
KOKKOS_INLINE_FUNCTION
Future< typename Impl::FunctorTraits< PatternType , Closure1Type >::device_type >
spawn( typename PatternType::work_type const & work 
     , Closure1Type const & closure1
     , FutureType const * const ibegin
     , FutureType const * const iend
     )
{
  typedef Impl::FunctorTraits< PatternType , Closure1Type > traits ;
  typedef typename traits::device_type  device_type ;
  typedef Impl::Task< device_type >     task_type ;

  return Future< device_type >( task_type::template task_spawn< traits >( work , closure1 , ibegin , iend ) );
}

//----------------------------------------------------------------------------

/**\brief  Respawn an executing task with dependences on other tasks */
template< class FunctorType , class FutureType >
KOKKOS_INLINE_FUNCTION
void respawn( FunctorType * const functor 
            , FutureType const * const ibegin
            , FutureType const * const iend
            )
{
  functor->Impl::Task< typename Impl::FunctorTraits< FunctorType >::device_type >
    ::task_respawn( ibegin , iend - ibegin );
}

//----------------------------------------------------------------------------

/**\brief Query the i^th dependence of this task */
template< class FunctorType >
KOKKOS_INLINE_FUNCTION
Future< typename Impl::FunctorTraits< FunctorType >::device_type >
task_dependence( FunctorType const * const functor , const int i )
{
  return Future< typename Impl::FunctorTraits< FunctorType >::device_type >( functor , i );
}

} /* namespace Kokkos */

//----------------------------------------------------------------------------

#endif /* #define KOKKOS_TASK_HPP */

