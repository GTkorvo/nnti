/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_FIXTURES_QUAD_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_QUAD_MESH_FIXTURE_HPP

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/DataTraits.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

class QuadFixture {

public:
  typedef int Scalar ;
  typedef Field<Scalar, Cartesian>    CoordFieldType;
  typedef Field<Scalar*,ElementNode>  CoordGatherFieldType;

  QuadFixture( stk::ParallelMachine pm, unsigned nx , unsigned ny );

  ~QuadFixture();

  stk::mesh::MetaData    meta_data ;
  stk::mesh::BulkData    bulk_data ;
  stk::mesh::Part      & quad_part ;
  CoordFieldType       & coord_field ;
  CoordGatherFieldType & coord_gather_field ;
  const unsigned         NX ;
  const unsigned         NY ;

  stk::mesh::EntityId node_id( unsigned ix , unsigned iy ) const
    { return 1 + ix + ( NX + 1 ) * iy ; }

  stk::mesh::EntityId elem_id( unsigned ix , unsigned iy ) const
    { return 1 + ix + NX * iy ; }

  stk::mesh::Entity * node( unsigned ix , unsigned iy ) const
    { return bulk_data.get_entity( stk::mesh::Node , node_id(ix,iy) ); }

  stk::mesh::Entity * elem( unsigned ix , unsigned iy ) const
    { return bulk_data.get_entity( stk::mesh::Element , elem_id(ix,iy) ); }

  void generate_mesh( std::vector<EntityId> & element_ids_on_this_processor );

  void generate_mesh();

private:

  QuadFixture();
  QuadFixture( const QuadFixture & );
  QuadFixture & operator = ( const QuadFixture & );
};

} // fixtures
} // mesh
} // stk
#endif
