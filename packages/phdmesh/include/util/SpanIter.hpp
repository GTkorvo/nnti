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
 * @date   November 2007
 */

#ifndef util_SpanIter_hpp
#define util_SpanIter_hpp

#include <iterator>

namespace phdmesh {

/** @class SpanIter
 *  Iterate a span of a container defined by begin and end iterators.
 */

template< class IterType ,
          class IterCategory =
             typename std::iterator_traits< IterType >::iterator_category >
class SpanIter ;

//----------------------------------------------------------------------
// Only defined for random access iterators

template< class IterType >
class SpanIter< IterType , std::random_access_iterator_tag > {
public:
  typedef IterType iterator ;

private:
  typedef SpanIter< iterator , std::random_access_iterator_tag > Self ;
  typedef std::iterator_traits< iterator > Traits ;

  iterator m_end ;
  iterator m_iter ;
public:

  //--------------------------------
  // Forward iterator requirements

  typedef typename Traits::value_type      value_type ;
  typedef typename Traits::difference_type difference_type ;
  typedef typename Traits::pointer         pointer ;
  typedef typename Traits::reference       reference ;
  typedef std::forward_iterator_tag        iterator_category ;

  ~SpanIter() {}

  SpanIter() : m_end() { m_iter = m_end ; }

  SpanIter( const Self & rhs ) : m_end( rhs.m_end ) , m_iter( rhs.m_iter ) {}

  Self & operator = ( const Self & rhs )
    { m_end = rhs.m_end ; m_iter = rhs.m_iter ; return *this ; }

  bool operator == ( const Self & rhs ) const
    { return m_end == rhs.m_end && m_iter == rhs.m_iter ; }

  bool operator != ( const Self & rhs ) const
    { return m_end != rhs.m_end || m_iter != rhs.m_iter ; }

  Self & operator ++ () { ++m_iter ; return *this ; }

  Self operator ++ (int) { Self tmp(*this); ++m_iter ; return tmp ; }

  reference operator * ()  const { return *m_iter ; }
  pointer   operator -> () const { return & *m_iter ; }

  //--------------------------------
  // Additional 'span' functionality

  iterator begin() const { return m_iter ; }
  iterator end()   const { return m_end ; }

  SpanIter( iterator i , iterator e ) : m_end(e), m_iter(i) {}

  template<class Container>
  explicit
  SpanIter( const Container & c ) : m_end( c.end() ), m_iter( c.begin() ) {}

  template<class Container>
  explicit
  SpanIter( Container & c ) : m_end( c.end() ), m_iter( c.begin() ) {}

  bool empty () const { return ! ( m_iter < m_end ) ; }

  operator bool () const { return m_iter < m_end ; }

  //--------------------------------
  // Additonal random access functionality

  reference operator [] ( difference_type n ) const { return m_iter[n] ; }

  difference_type size() const { return std::distance( m_iter , m_end ); }
};

} // namespace phdmesh

#endif

