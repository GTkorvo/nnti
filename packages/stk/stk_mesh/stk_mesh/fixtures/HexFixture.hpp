/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_FIXTURES_HEX_MESH_FIXTURE_HPP
#define STK_MESH_FIXTURES_HEX_MESH_FIXTURE_HPP

#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/DataTraits.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/TopologyDimensions.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

namespace {

const double PI     = 3.14159265358979;
const double TWO_PI = 2 * PI;

} // namespace


/**
 * A mapping for the coordinates as a function of the node indices
 */
class CoordinateMapping
{
public:
  typedef double Scalar;
  CoordinateMapping() {}
  virtual void getNodeCoordinates(Scalar * field, const unsigned nx, const unsigned ny, const unsigned nz) const = 0;
  virtual ~CoordinateMapping() {};
};


/**
 * Standard Cartesian X-Y-Z coordinate mapping
 */
class CartesianCoordinateMapping : public CoordinateMapping
{
public:
  CartesianCoordinateMapping() : CoordinateMapping() {}
  virtual void getNodeCoordinates(Scalar * field, const unsigned nx, const unsigned ny, const unsigned nz) const
  {
    field[0] = nx;
    field[1] = ny;
    field[2] = nz;
  }
};


/**
 * Cylindrical coordinate mapping where the standard X-Y-Z coordinates
 * are mapped onto a wedge slice specified by the inner radius and the
 * maximum angle theta.
 */
class CylindricalCoordinateMapping : public CoordinateMapping
{
public:
  CylindricalCoordinateMapping(Scalar radius, Scalar theta, unsigned numTheta)
      : CoordinateMapping(), m_radius(radius), m_theta(theta), m_numTheta(numTheta)
  { }
  virtual void getNodeCoordinates(Scalar * field, const unsigned nx, const unsigned ny, const unsigned nz) const
  {
    Scalar fracTheta = nx/(m_numTheta - 1);

    //we want the angle to go from pi/2 to pi/2 - theta so we do not
    //invert any elements
    Scalar angle = PI/2.0 + m_theta*fracTheta;
    field[0] = (m_radius + ny)*std::cos(angle);
    field[1] = (m_radius + ny)*std::sin(angle);
    field[2] = nz;
  }
private:
  Scalar m_radius;
  Scalar m_theta;
  unsigned m_numTheta;


};

/**
 * A 3-dimensional X*Y*Z hex fixture.
 *
 * A coordinate field will be added to all nodes, a coordinate-gather field
 * will be added to all elements.
 */
class HexFixture
{
 public:
  typedef double                     Scalar ;
  typedef Field<Scalar, Cartesian>   CoordFieldType;

  /**
   * Set up meta data to support this fixture. Meta data is left uncommitted
   * to allow additional modifications by the client.
   */
  HexFixture(   stk::ParallelMachine pm
              , unsigned nx
              , unsigned ny
              , unsigned nz
              , ConnectivityMap const* connectivity_map = NULL
            );

  const int                     m_spatial_dimension;
  const unsigned                m_nx;
  const unsigned                m_ny;
  const unsigned                m_nz;
  MetaData                      m_meta;
  BulkData                      m_bulk_data;
  PartVector                    m_elem_parts;
  PartVector                    m_node_parts;
  CoordFieldType &              m_coord_field ;


  /**
   * Thinking in terms of a 3D grid of nodes, get the id of the node in
   * the (x, y, z) position.
   */
  EntityId node_id( unsigned x , unsigned y , unsigned z ) const  {
    return 1 + x + ( m_nx + 1 ) * ( y + ( m_ny + 1 ) * z );
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the id of the
   * element in the (x, y, z) position.
   */
  EntityId elem_id( unsigned x , unsigned y , unsigned z ) const  {
    return 1 + x + m_nx * ( y + m_ny * z );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, get the node in the (x, y, z)
   * position. Return NULL if this process doesn't know about this node.
   */
  Entity node( unsigned x , unsigned y , unsigned z ) const {
    return m_bulk_data.get_entity( MetaData::NODE_RANK , node_id(x, y, z) );
  }

  /**
   * Thinking in terms of a 3D grid of elements, get the elements in the
   * (x, y, z) position. Return NULL if this process doesn't know about this
   * element.
   */
  Entity elem( unsigned x , unsigned y , unsigned z ) const {
    return m_bulk_data.get_entity( MetaData::ELEMENT_RANK, elem_id(x, y, z) );
  }

  /**
   * Thinking in terms of a 3D grid of nodes, compute the (x, y, z) position
   * of a node given it's id.
   */
  void node_x_y_z( EntityId entity_id, unsigned &x , unsigned &y , unsigned &z ) const;

  /**
   * Thinking in terms of a 3D grid of elements, compute the (x, y, z) position
   * of an element given it's id.
   */
  void elem_x_y_z( EntityId entity_id, unsigned &x , unsigned &y , unsigned &z ) const;

  /**
   * Create the mesh (into m_bulk_data).
   */
  void generate_mesh(const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  void generate_mesh( std::vector<EntityId> & element_ids_on_this_processor, const CoordinateMapping & coordMap = CartesianCoordinateMapping());

  // When creating entities, you can tell HexFixture what parts to add
  // elements and nodes.

  template <typename Iterator>
  void add_elem_parts(Iterator itr, size_t num)
  {
    ThrowRequire(!m_meta.is_commit());
    m_elem_parts.insert(m_elem_parts.end(), itr, itr + num);
  }

  template <typename Iterator>
  void add_node_parts(Iterator itr, size_t num)
  {
    ThrowRequire(!m_meta.is_commit());
    m_node_parts.insert(m_node_parts.end(), itr, itr + num);
  }

 private:

  HexFixture();
  HexFixture( const HexFixture &);
  HexFixture & operator = (const HexFixture &);
};


} // fixtures
} // mesh
} // stk
#endif
