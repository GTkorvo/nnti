/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <Shards_BasicTopologies.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <stk_mesh/base/GetEntities.hpp>

#include <stk_linsys/DofMapper.hpp>
#include <stk_linsys/ImplDetails.hpp>

#include <unit_tests/UnitTest_helpers.hpp>

static const stk::mesh::EntityRank NODE_RANK = stk::mesh::MetaData::NODE_RANK;

namespace stk_linsys_unit_tests {

static const size_t spatial_dimension = 3;

//------------- here is the DofMapper unit-test... -----------------------

void testDofMapper( MPI_Comm comm )
{
  //First create and fill MetaData and BulkData objects:

  const unsigned bucket_size = 100; //for a real application mesh, bucket_size would be much bigger...

  stk::mesh::MetaData fem_meta;
  fem_meta.initialize(spatial_dimension);

  const stk::mesh::EntityRank element_rank = stk::mesh::MetaData::ELEMENT_RANK;

  stk::mesh::BulkData bulk_data( fem_meta, comm, bucket_size );

  fill_utest_mesh_meta_data( fem_meta );
  fill_utest_mesh_bulk_data( bulk_data );

  stk::mesh::Selector selector = fem_meta.locally_owned_part() | fem_meta.globally_shared_part() ;
  std::vector<unsigned> count;
  stk::mesh::count_entities(selector, bulk_data, count);

  STKUNIT_ASSERT_EQUAL( count[element_rank], (unsigned)4 );
  STKUNIT_ASSERT_EQUAL( count[NODE_RANK],    (unsigned)20 );

  std::vector<stk::mesh::Entity> nodes;
  stk::mesh::get_entities(bulk_data, NODE_RANK, nodes);

  ScalarField* temperature_field = fem_meta.get_field<ScalarField>("temperature");

  //Now we're ready to test the DofMapper:

  stk::linsys::DofMapper dof_mapper(comm);

  const stk::mesh::Selector select_used = fem_meta.locally_owned_part() | fem_meta.globally_shared_part();

  dof_mapper.add_dof_mappings(bulk_data, select_used,
                              NODE_RANK, *temperature_field);

  stk::mesh::EntityRank ent_type=0;
  stk::mesh::EntityId ent_id=0;
  const stk::mesh::FieldBase* field = NULL;
  int offset_into_field=0;
  int index = 0;
  //DofMapper::get_dof can't be called until after DofMapper::finalize() has
  //been called.
  //We'll call it now to verify that an exception is thrown:
  std::cout << "Testing error condition: " << std::endl;
  STKUNIT_ASSERT_THROW(dof_mapper.get_dof(index, ent_type, ent_id, field, offset_into_field), std::runtime_error );
  std::cout << "...Completed testing error condition." << std::endl;

  dof_mapper.finalize();

  //find a node that is in the locally-used part:
  size_t i_node = 0;
  while(! select_used( bulk_data.bucket(nodes[i_node]) ) && i_node<nodes.size()) {
    ++i_node;
  }

  //test the get_global_index function:
  stk::mesh::EntityId node_id = bulk_data.identifier(nodes[i_node]);
  index = dof_mapper.get_global_index(NODE_RANK, node_id, *temperature_field);
  STKUNIT_ASSERT_EQUAL( index, (int)(node_id-1) );

  std::cout << "Testing error condition: " << std::endl;
  //call DofMapper::get_global_index with a non-existent ID and verify that an
  //exception is thrown:
  STKUNIT_ASSERT_THROW(dof_mapper.get_global_index(NODE_RANK, (stk::mesh::EntityId)999999, *temperature_field), std::runtime_error);
  std::cout << "...Completed testing error condition." << std::endl;

  int numProcs = 1;
  numProcs = stk::parallel_machine_size( MPI_COMM_WORLD );

  fei::SharedPtr<fei::VectorSpace> fei_vspace = dof_mapper.get_fei_VectorSpace();
  int numIndices = fei_vspace->getGlobalNumIndices();
  STKUNIT_ASSERT_EQUAL( numIndices, (int)(numProcs*20 - (numProcs-1)*4) );

  dof_mapper.get_dof(index, ent_type, ent_id, field, offset_into_field);

  STKUNIT_ASSERT_EQUAL( ent_type, bulk_data.entity_rank(nodes[i_node]) );
  STKUNIT_ASSERT_EQUAL( ent_id,   bulk_data.identifier(nodes[i_node]) );
  STKUNIT_ASSERT_EQUAL( field->name() == temperature_field->name(), true );
  STKUNIT_ASSERT_EQUAL( offset_into_field, (int)0 );
}

} // namespace stk_linsys_unit_tests

STKUNIT_UNIT_TEST(UnitTestingOfDofMapper, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_linsys_unit_tests::testDofMapper ( MPI_COMM_WORLD );
}

