/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 */

#include <stk_percept/fixtures/BeamFixture.hpp>
//
//----------------------------------------------------------------------

#include <iostream>

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/base/TopologyDimensions.hpp>

//----------------------------------------------------------------------

namespace stk{
  namespace percept {

    typedef shards::Beam<2>               Beam2;

    BeamFixture::BeamFixture( stk::ParallelMachine comm, bool doCommit ) :
      m_spatial_dimension(3)
      , m_metaData(m_spatial_dimension )
      , m_bulkData(m_metaData, comm )
      , m_block_beam( m_metaData.declare_part< Beam2 >( "block_2" ) )
      , m_elem_rank( stk::mesh::MetaData::ELEMENT_RANK )
      , m_coordinates_field( m_metaData.declare_field< VectorFieldType >( "coordinates" ))
      , m_centroid_field(    m_metaData.declare_field< VectorFieldType >( "centroid" ))
      , m_temperature_field( m_metaData.declare_field< ScalarFieldType >( "temperature" ))
      , m_volume_field( m_metaData.declare_field< ScalarFieldType >( "volume" ))
    {
      // Define where fields exist on the mesh:
      stk::mesh::Part & universal = m_metaData.universal_part();

      put_field( m_coordinates_field , stk::mesh::MetaData::NODE_RANK , universal );
      put_field( m_centroid_field , m_elem_rank , universal );
      put_field( m_temperature_field, stk::mesh::MetaData::NODE_RANK, universal );
      put_field( m_volume_field, m_elem_rank, m_block_beam );

      if (doCommit)
        m_metaData.commit();
    }

    BeamFixture::~BeamFixture()
    { }

    //------------------------------------------------------------------------------
    // Use case specific mesh generation data:

    enum { SpatialDim = 3 };
    enum { node_count = 4 };
    enum { number_beam = 2 };

    namespace {

      // Hard coded node coordinate data for all the nodes in the entire mesh
      static const double node_coord_data[ node_count ][ SpatialDim ] = {
        { 10 , 10 , 10 } , { 11 , 12 , 13 } , { 24 , 24 , 24 } , { 28 , 29 , 30 }
      };

      // Hard coded beam node ids for all the beam nodes in the entire mesh
      static const stk::mesh::EntityId beam_node_ids[number_beam][ Beam2::node_count ] = {
        { 1,2 } ,
        { 3,4 }
      };

    }

    //------------------------------------------------------------------------------

    void BeamFixture::populate()
    {
      // Populate mesh with all node types

      m_bulkData.modification_begin();

      if (m_bulkData.parallel_rank() == 0)
        {
          stk::mesh::EntityId curr_elem_id = 1;

          // For each element topology declare elements
          for ( unsigned i = 0 ; i < number_beam ; ++i , ++curr_elem_id ) {
            stk::mesh::declare_element( m_bulkData, m_block_beam, curr_elem_id, beam_node_ids[i] );
          }

          // For all nodes assign nodal coordinates
          for ( unsigned i = 0 ; i < node_count ; ++i ) {
            stk::mesh::Entity const node = m_bulkData.get_entity( stk::mesh::MetaData::NODE_RANK , i + 1 );
            double * const coord = m_bulkData.field_data( m_coordinates_field , node );
            coord[0] = node_coord_data[i][0] ;
            coord[1] = node_coord_data[i][1] ;
            coord[2] = node_coord_data[i][2] ;
          }

        }
      m_bulkData.modification_end();

    }

    // Verify mesh for 6 different parts
    bool verifyMesh( const BeamFixture & mesh )
    {
      bool result = true;

      return result;
    }

  } //namespace percept
} //namespace stk
