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
 * @author H. Carter Edwards
 */

#ifndef phdmesh_Kernel_hpp
#define phdmesh_Kernel_hpp

//----------------------------------------------------------------------

#include <iosfwd>

#include <util/Setv.hpp>
#include <mesh/Types.hpp>
#include <mesh/Part.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

std::ostream & operator << ( std::ostream & , const Kernel & );

/** A given PartSet, i.e. the intersection of a set of parts,
 *  has one or more associated Kernels.  Each kernel is identified
 *  by its PartSet and ordinal.
 */
struct KernelLess {
  bool operator()( const unsigned * lhs , const unsigned * rhs ) const ;
};

/** Mesh kernel - a homogeneous collection of mesh entities and
 *  their field values.
 *  Homogeneous in that all entities are of the same entity type
 *  and are members of the same set of parts.
 */
class Kernel : private SetvMember<const unsigned * const> {
private:
  struct DataMap {
    const unsigned * m_stride ;
    unsigned         m_base ;
    unsigned         m_size ;
  };

  Mesh      & m_mesh ;        // Mesh in which this kernel resides
  EntityType  m_entity_type ; // Type of mesh entities
  unsigned    m_size ;        // Number of entities
  unsigned    m_capacity ;    // Capacity for entities
  DataMap   * m_field_map ;   // Field value data map, shared
  Entity   ** m_entities ;    // Array of entity pointers,
                              // begining of field value memory.

public:

  Mesh & mesh() const { return m_mesh ; }

  EntityType entity_type() const { return m_entity_type ; }

  /** This kernel's supersets */
  void supersets( PartSet & ) const ;

  /** This kernel's supersets' schema ordinals */
  void supersets( std::vector<unsigned> & ) const ;

  /** The given part is a superset of this kernel */
  bool has_superset( const Part & ) const ;

  /** The given ordered set of parts are all supersets of this kernel */
  bool has_superset( const PartSet & ) const ;

  //--------------------------------
  /** Number of entities */
  unsigned size() const { return m_size ; }

  /** Capacity of entities */
  unsigned capacity() const { return m_capacity ; }

  /** Access the i^th entity */
  Entity * operator[] ( unsigned i ) const { return m_entities[i] ; }

  typedef Entity * const * iterator ;
  iterator begin() const { return m_entities ; }
  iterator end()   const { return m_entities + m_size ; }

  //--------------------------------

  std::ostream & print( std::ostream & , const std::string & ) const ;

  ~Kernel();

private:

  Kernel();
  Kernel( const Kernel & );
  Kernel & operator = ( const Kernel & );

  Kernel( Mesh & , EntityType , const unsigned * );

  void update_state();

  friend class Mesh ;
  friend class Setv<Kernel,KernelLess> ;
  friend class SetvIter<Kernel,true> ;
  friend class SetvIter<Kernel,false> ;
  friend class SetvIter<const Kernel,true> ;
  friend class SetvIter<const Kernel,false> ;

  static void copy_fields( Kernel & k_dst , unsigned i_dst ,
                           Kernel & k_src , unsigned i_src );

  static void zero_fields( Kernel & k_dst , unsigned i_dst );

  template< class field_type >
  friend
  typename field_type::BlockDimension
  field_dimension( const field_type & f , const Kernel & k );

  template< class field_type >
  friend
  typename field_type::Dimension
  field_dimension( const field_type & f , const Entity & e );

  template< class field_type >
  friend
  typename field_type::data_type *
  field_data( const field_type & f , const Kernel & k );

  template< class field_type >
  friend
  typename field_type::data_type *
  field_data( const field_type & f , const Entity & e );

  friend
  unsigned field_data_size( const FieldBase & f , const Kernel & k );
};

/** The set of mesh kernels is dynamic and potentially large. */

typedef Setv<Kernel,KernelLess> KernelSet ;

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

