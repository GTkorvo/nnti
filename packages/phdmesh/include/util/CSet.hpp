/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   November 2006
 */

#ifndef util_CSet_hpp
#define util_CSet_hpp

#include <typeinfo>
#include <vector>

namespace phdmesh {

//----------------------------------------------------------------------
/**
 * \class CSet
 * \brief Set of entities of arbitrary types.
 *
 *  Example usage of the three methods:
 *
 * <PRE>
 *  class A { ... };
 *  class B { ... };
 *
 *  CSet cset ;
 *
 *  // Insert pointers to objects:
 *
 *  cset.insert<A>( new A ); // Do not delete on destruction
 *  cset.insert<B>( new B , true ); // Delete on destruction
 *
 *  // Query the collection of objects of a given type:
 *
 *  const A * sa = cset.get<A>();
 *  const B * sb = cset.get<B>();
 *
 *  // Remove a member:
 *
 *  {
 *    B * b = ... ;
 *    cset.remove<B>( b ); // Remove never deletes
 *    delete b ;
 *  }
 * </PRE>
 */
class CSet {
public:

  /** Get member conforming to the given type. */
  template<class T> const T * get() const ;

  /** Insert and optionally request deletion upon destruction.
   *  If already exists then return existing member, insert fails.
   */
  template<class T> const T * insert( const T * , bool = false );

  /** Erase a member without deleting.
   *  Return if the remove operation was successful.
   */
  template<class T> bool remove( const T * );

  //--------------------------------

  ~CSet();
  CSet();

private:

  typedef void (*DeleteFunction)(void *);

  typedef std::pair< const std::type_info * , DeleteFunction > Manager ;

  const void * p_get( const std::type_info & ) const ;

  const void * p_insert( const Manager & , const void * );

  bool p_remove( const std::type_info & , const void * );

  std::vector< Manager > m_manager ;
  std::vector< const void * > m_value ;

  CSet( const CSet & );
  CSet & operator = ( const CSet & );
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#ifndef DOXYGEN_COMPILE

// Inlined template methods have casting.

namespace phdmesh {

namespace {
template<class T>
void cset_member_delete( void * v ) { delete reinterpret_cast<T*>( v ); }
}

template<class T>
inline
const T * CSet::get() const
{ return (const T*) p_get( typeid(T) ); }

template<class T>
inline
const T * CSet::insert( const T * arg_value , bool arg_delete )
{
  Manager m ;
  m.first = & typeid(T);
  m.second = NULL ;

  if ( arg_delete ) { m.second = & cset_member_delete<T> ; }

  return (const T *) p_insert( m , arg_value );
}

template<class T>
inline
bool CSet::remove( const T * arg_value )
{ return p_remove( typeid(T) , arg_value ); }

} // namespace phdmesh

#endif /* DOXYGEN_COMPILE */

#endif /* util_CSet_hpp */


