/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <mesh/UseCase_Skinning.hpp>

#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>

#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

namespace {

static const size_t NODE_RANK = stk::mesh::MetaData::NODE_RANK;

unsigned count_skin_entities( stk::mesh::BulkData & mesh, stk::mesh::Part & skin_part, stk::mesh::EntityRank skin_rank )
{
  const stk::mesh::MetaData & meta = stk::mesh::MetaData::get(mesh);

  stk::mesh::Selector select_skin = skin_part & meta.locally_owned_part()  ;

  const std::vector<stk::mesh::Bucket*>& buckets = mesh.buckets( skin_rank );

  return count_selected_entities( select_skin, buckets);
}

} // empty namespace

bool skinning_use_case_1(stk::ParallelMachine pm)
{
  bool passed = true;
  {
    //setup the mesh
    stk::mesh::fixtures::HexFixture fixture(pm,3,3,3);

    stk::mesh::MetaData & fem_meta = fixture.m_meta;
    stk::mesh::BulkData & mesh = fixture.m_bulk_data;
    const stk::mesh::EntityRank element_rank = stk::mesh::MetaData::ELEMENT_RANK;
    const stk::mesh::EntityRank side_rank    = fem_meta.side_rank();

    stk::mesh::Part & skin_part = fem_meta.declare_part("skin_part");
    fem_meta.commit();

    fixture.generate_mesh();

    skin_mesh(mesh, element_rank, &skin_part);

    std::vector< stk::mesh::EntityId > elements_to_separate;

    //separate out the middle element
    elements_to_separate.push_back(fixture.elem_id(1,1,1));

    separate_and_skin_mesh(
        fem_meta,
        mesh,
        skin_part,
        elements_to_separate,
        element_rank
        );

    // pointer to middle_element after mesh modification.
    stk::mesh::Entity middle_element = mesh.get_entity(element_rank, fixture.elem_id(1,1,1));

    unsigned num_skin_entities = count_skin_entities(mesh, skin_part, side_rank);

    stk::all_reduce(pm, stk::ReduceSum<1>(&num_skin_entities));

    //there should be 66 faces in the skin part
    //54 on the outside
    //6 on the inside attached to the entire mesh
    //6 on the inside attected to the element that was detached
    bool correct_skin = ( num_skin_entities == 66 );
    bool correct_relations = true;
    bool correct_comm = true;

    //all nodes connected to the single element that has been broken off
    //should have relations.size() == 4 and comm.size() == 0
    if (mesh.is_valid(middle_element) && mesh.parallel_owner_rank(middle_element) == mesh.parallel_rank()) {

      stk::mesh::Entity const *rel_nodes_iter = mesh.begin_nodes(middle_element);
      stk::mesh::Entity const *rel_nodes_end = mesh.end_nodes(middle_element);

      for (; rel_nodes_iter != rel_nodes_end; ++rel_nodes_iter) {
        stk::mesh::Entity current_node = *rel_nodes_iter;
        //each node should be attached to only 1 element and 3 faces
        correct_relations &= ( mesh.count_relations(current_node) == 4 );
        //the entire closure of the element should exist on a single process
        correct_comm      &= ( mesh.entity_comm(mesh.entity_key(current_node)).size() == 0 );
      }
    }
    passed &= (correct_skin && correct_relations && correct_comm);
  }

  //seperate the entire middle layer of the mesh
  {
    //setup the mesh
    stk::mesh::fixtures::HexFixture fixture(pm,3,3,3);

    stk::mesh::MetaData & fem_meta = fixture.m_meta;
    stk::mesh::BulkData & mesh = fixture.m_bulk_data;
    const stk::mesh::EntityRank element_rank = stk::mesh::MetaData::ELEMENT_RANK;
    const stk::mesh::EntityRank side_rank    = fem_meta.side_rank();

    stk::mesh::Part & skin_part = fem_meta.declare_part("skin_part");
    fem_meta.commit();

    fixture.generate_mesh();

    skin_mesh(mesh, element_rank, &skin_part);

    std::vector< stk::mesh::EntityId > elements_to_separate;

    //separate out the middle level
    elements_to_separate.push_back(fixture.elem_id(1,0,0));
    elements_to_separate.push_back(fixture.elem_id(1,0,1));
    elements_to_separate.push_back(fixture.elem_id(1,0,2));
    elements_to_separate.push_back(fixture.elem_id(1,1,0));
    elements_to_separate.push_back(fixture.elem_id(1,1,1));
    elements_to_separate.push_back(fixture.elem_id(1,1,2));
    elements_to_separate.push_back(fixture.elem_id(1,2,0));
    elements_to_separate.push_back(fixture.elem_id(1,2,1));
    elements_to_separate.push_back(fixture.elem_id(1,2,2));

    separate_and_skin_mesh(
        fem_meta,
        mesh,
        skin_part,
        elements_to_separate,
        element_rank
        );

    // pointer to middle_element after mesh modification.
    unsigned num_skin_entities = count_skin_entities(mesh, skin_part, side_rank);

    stk::all_reduce(pm, stk::ReduceSum<1>(&num_skin_entities));

    //there should be 90 faces in the skin part
    //30 attached to each level of the mesh
    bool correct_skin = ( num_skin_entities == 90 );

    passed &= correct_skin;
  }

  return passed;
}
