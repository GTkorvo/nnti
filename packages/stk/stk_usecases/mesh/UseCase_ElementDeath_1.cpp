// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_mesh/fixtures/GridFixture.hpp>

#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <mesh/UseCase_ElementDeath_1_validation_helpers.hpp>

/*
The grid fixture creates the mesh below and skins it
1-16 Quadrilateral<4>
17-41 Nodes
skin ids are generated by the distributed index

Note:  "=" and "||" represent side entities.

17===18===19===20===21
|| 1 |  2 |  3 |  4 ||
22---23---24---25---26
|| 5 |  6 |  7 |  8 ||
27---28---29---30---31
|| 9 | 10 | 11 | 12 ||
32---33---34---35---36
|| 13| 14 | 15 | 16 ||
37===38===39===40===41

This use case will iteratively erode the mesh.

Each iteration will move a selection of faces to the 'dead_part"
Create boundaries between live and dead faces
Destroy nodes and sides that are no longer attached to a live face

0:  Init the mesh

17===18===19===20===21
|| 1 |  2 |  3 |  4 ||
22---23---24---25---26
|| 5 |  6 |  7 |  8 ||
27---28---29---30---31
|| 9 | 10 | 11 | 12 ||
32---33---34---35---36
|| 13| 14 | 15 | 16 ||
37===38===39===40===41

1: Move 4, 9 and 10 to the dead part


17===18===19===20
|| 1 |  2 |  3 ||
22---23---24---25===26
|| 5 |  6 |  7 |  8 ||
27===28===29---30---31
          || 11| 12 ||
32===33===34---35---36
|| 13| 14 | 15 | 16 ||
37===38===39===40===41


2: Move faces 2 and 3 to the dead part

17===18
|| 1 ||
22---23===24===25===26
|| 5 |  6 |  7 |  8 ||
27===28===29---30---31
          || 11| 12 ||
32===33===34---35---36
|| 13| 14 | 15 | 16 ||
37===38===39===40===41

3: Move faces 1 and 11 to the dead part

22===23===24===25===26
|| 5 |  6 |  7 |  8 ||
27===28===29===30---31
               || 12||
32===33===34===35---36
|| 13| 14 | 15 | 16 ||
37===38===39===40===41

4: Move faces 6 and 7 to the dead part

22===23        25===26
|| 5 ||        ||  8||
27===28        30---31
               || 12||
32===33===34===35---36
|| 13| 14 | 15 | 16 ||
37===38===39===40===41

5: Move faces 5 and 16 to the dead part

               25===26
               ||  8||
               30---31
               || 12||
32===33===34===35===36
|| 13| 14 | 15 ||
37===38===39===40

6: Move the remaining faces to the dead part


(this space intentionally left blank)
  (nothing to see here)

*/

const int NUM_ITERATIONS = 7;

namespace {

//Finds the sides that need to be created between the live and dead entities
void find_sides_to_be_created(
    const stk::mesh::BulkData &mesh,
    const stk::mesh::EntitySideVector & boundary,
    const stk::mesh::Selector & select,
    std::vector<stk::mesh::EntitySideComponent> & sides,
    stk::mesh::EntityRank side_rank
    );

//Finds entities from the closure of the entities_to_be_killed
//that only have relations with dead entities
void find_lower_rank_entities_to_kill(
    const stk::mesh::BulkData &mesh,
    const stk::mesh::EntityVector & entities_closure,
    stk::mesh::EntityRank closure_rank,
    stk::mesh::EntityRank entity_rank,
    unsigned spatial_dim,
    const stk::mesh::Selector & select_owned,
    const stk::mesh::Selector & select_live,
    stk::mesh::EntityVector & kill_list
    );

}

bool element_death_use_case_1(stk::ParallelMachine pm)
{
  //set up the mesh
  stk::mesh::fixtures::GridFixture fixture(pm);

  stk::mesh::BulkData& mesh = fixture.bulk_data();
  stk::mesh::MetaData& fem_meta = fixture.fem_meta();
  const stk::mesh::EntityRank element_rank = stk::topology::ELEMENT_RANK;
  const stk::mesh::EntityRank side_rank = fem_meta.side_rank();

  fem_meta.commit();

  mesh.modification_begin();
  fixture.generate_grid();
  mesh.modification_end();

  stk::mesh::skin_mesh(mesh);

  // Nothing happens on iteration #0,
  // so the initial mesh should pass this validation.

  if ( ! validate_iteration( pm, fixture, 0) ) { return false ; }

  stk::mesh::Part & dead_part = *fixture.dead_part();

  stk::mesh::PartVector dead_parts;
  dead_parts.push_back( & dead_part);

  bool passed = true;

  stk::mesh::EntityRank mesh_rank = element_rank;
  stk::mesh::Part &line2_part = fem_meta.get_topology_root_part(stk::topology::LINE_2);
  stk::mesh::PartVector add_line2_parts, empty_parts;
  add_line2_parts.push_back(&line2_part);

  for (int iteration = 0; iteration <NUM_ITERATIONS; ++iteration) {
    //find the entities to kill in this iteration
    stk::mesh::EntityVector entities_to_kill = entities_to_be_killed(mesh, iteration, element_rank);

    // find the parallel-consistent closure of the entities to be killed
    // The closure of an entity includes the entity and any lower ranked
    // entities which are reachable through relations.  For example, the
    // closure of an element consist of the element and the faces, edges,
    // and nodes that are attached to the element through relations.
    //
    // The find closure function will return a sorted parallel consistent vector
    // which contains all the entities that make up the closure of the input
    // vector.
    stk::mesh::EntityVector entities_closure;
    stk::mesh::find_closure(mesh,
        entities_to_kill,
        entities_closure);

    // find the boundary of the entities we're killing
    stk::mesh::EntitySideVector boundary;
    stk::mesh::boundary_analysis(mesh,
        entities_closure,
        mesh_rank,
        boundary);

    // Find the sides that need to be created.
    // Sides need to be created when the outside
    // of the boundary is both live and owned and
    // a side separating the live and dead doesn't
    // already exist.
    stk::mesh::Selector select_owned = fem_meta.locally_owned_part();
    stk::mesh::Selector select_live = ! dead_part ;
    stk::mesh::Selector select_live_and_owned = select_live & select_owned;

    std::vector<stk::mesh::EntitySideComponent> skin;
    find_sides_to_be_created(mesh, boundary, select_live_and_owned, skin, side_rank);

    mesh.modification_begin();

    // Kill entities by moving them to the dead part.
    for (stk::mesh::EntityVector::iterator itr = entities_to_kill.begin();
        itr != entities_to_kill.end(); ++itr) {
      mesh.change_entity_parts(*itr, dead_parts);
    }

    // Ask for new entities to represent the sides between the live and dead entities
    //
    std::vector<size_t> requests(fem_meta.entity_rank_count(), 0);
    requests[side_rank] = skin.size();

    // generate_new_entities creates new blank entities of the requested ranks
    stk::mesh::EntityVector requested_entities;
    mesh.generate_new_entities(requests, requested_entities);

    // Create boundaries between live and dead entities
    // by creating a relation between the new entities and the live entities
    for ( size_t i = 0; i < skin.size(); ++i) {
      stk::mesh::Entity entity = skin[i].entity;
      const unsigned side_ordinal  = skin[i].side_ordinal;
      stk::mesh::Entity side   = requested_entities[i];
      mesh.change_entity_parts(side, add_line2_parts, empty_parts);

      stk::mesh::declare_element_side(mesh, entity, side, side_ordinal);
    }

    mesh.modification_end();
    //the modification_end() will communicate which entities have been changed
    //to other processes.

    //find lower ranked entity that are only related to the dead entities
    //and kill them
    for (stk::mesh::EntityRank irank = stk::topology::BEGIN_RANK; irank <= side_rank; ++irank) {
      stk::mesh::EntityRank rank = static_cast<stk::mesh::EntityRank>(side_rank - irank);
      stk::mesh::EntityVector kill_list;
      find_lower_rank_entities_to_kill(
          mesh,
          entities_closure,
          mesh_rank,
          rank,
          fem_meta.spatial_dimension(),
          select_owned,
          select_live,
          kill_list
          );

      //need to communicate killing the higher ranking entities among
      //processors before killing the lower.
      mesh.modification_begin();
      for (stk::mesh::EntityVector::iterator itr = kill_list.begin();
          itr != kill_list.end(); ++itr) {
        mesh.change_entity_parts(*itr, dead_parts);
      }
      mesh.modification_end();
    }

    passed &= validate_iteration( pm, fixture, iteration);
  }

  return passed;
}

//----------------------------------------------------------------------------------
namespace {

//----------------------------------------------------------------------------------
void find_sides_to_be_created(
    const stk::mesh::BulkData &mesh,
    const stk::mesh::EntitySideVector & boundary,
    const stk::mesh::Selector & select,
    std::vector<stk::mesh::EntitySideComponent> & sides,
    stk::mesh::EntityRank side_rank
    )
{
  //look at the outside of the boundary since the inside will be kill this
  //iteration

  for (stk::mesh::EntitySideVector::const_iterator itr = boundary.begin();
      itr != boundary.end(); ++itr) {

    const stk::mesh::EntitySideComponent & outside = itr->outside;

    // examine the boundary of the outside of the closure.
    if ( mesh.is_valid(outside.entity) && select(mesh.bucket(outside.entity)) ) {

      //make sure the side does not already exist
      const int side_ordinal = outside.side_ordinal;
      const stk::mesh::Entity entity = outside.entity;
      stk::mesh::Entity const *existing_sides = mesh.begin(entity, side_rank);
      int num_sides = mesh.num_connectivity(entity, side_rank);
      stk::mesh::ConnectivityOrdinal const *existing_side_ordinals = mesh.begin_ordinals(entity, side_rank);

      int side_num = 0;
      for (;
           (side_num < num_sides)
               && (!existing_sides[side_num].is_local_offset_valid()
                   || (existing_side_ordinals[side_num] != static_cast<unsigned>(side_ordinal))) ;
           ++side_num);

      //reached the end -- a new side needs to be created
      if (side_num == num_sides) {
        sides.push_back(outside);
      }
    }
  }
}

//----------------------------------------------------------------------------------
void find_lower_rank_entities_to_kill(
    const stk::mesh::BulkData &mesh,
    const stk::mesh::EntityVector & entities_closure,
    stk::mesh::EntityRank mesh_rank,
    stk::mesh::EntityRank entity_rank,
    unsigned spatial_dim,
    const stk::mesh::Selector & select_owned,
    const stk::mesh::Selector & select_live,
    stk::mesh::EntityVector & kill_list
    )
{

  kill_list.clear();

  stk::mesh::EntityLess lesser(mesh);

  //find the first entity in the closure
  stk::mesh::EntityVector::const_iterator itr = std::lower_bound(entities_closure.begin(),
      entities_closure.end(),
      stk::mesh::EntityKey(entity_rank, 0),
      lesser);

  const stk::mesh::EntityVector::const_iterator end =
      std::lower_bound(entities_closure.begin(),
                       entities_closure.end(),
                       stk::mesh::EntityKey(static_cast<stk::mesh::EntityRank>(entity_rank+1), 0),
                       lesser);

  for (; itr != end; ++itr) {
    stk::mesh::Entity entity = *itr;

    stk::mesh::Bucket *bucket_ptr = mesh.bucket_ptr(entity);
    if (!bucket_ptr)
      continue;

    if (select_owned(*bucket_ptr)) {
      bool found_live = false;

      for(stk::mesh::EntityRank rank = static_cast<stk::mesh::EntityRank>(entity_rank + 1); rank<=mesh_rank && !found_live; ++rank) {
        if (spatial_dim == 2 && rank == stk::topology::FACE_RANK) {
          continue;
        }

        stk::mesh::Entity const * relations_iter = mesh.begin(entity, rank);
        stk::mesh::Entity const * relations_end = mesh.end(entity, rank);

        for (; relations_iter != relations_end && !found_live; ++relations_iter) {

          if( select_live(mesh.bucket(*relations_iter)) ) {
            found_live = true;
          }
        }
      }

      if (!found_live) {
        kill_list.push_back(entity);
      }
    }
  }

}

}
