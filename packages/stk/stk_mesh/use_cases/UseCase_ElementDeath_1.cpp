/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <use_cases/GridFixture.hpp>
#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

/*
The grid fixture creates the mesh below and skins it
1-16 Quadrilateral<4>
17-41 Nodes
skin ids are generated by the distributed index

17===18===19===20===21
|| 1 |  2 |  3 |  4 ||
22---23---24---25---26
|| 5 |  6 |  7 |  8 ||
27---28---29---30---31
|| 9 | 10 | 11 | 12 ||
32---33---34---35---36
|| 13| 14 | 15 | 16 ||
37===38===39===40===41

This use case will erode the mesh with  iterations

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


  (nothing to see here)

*/

const int NUM_ITERATIONS = 7;
const int NUM_RANK = 3;

const int global_num_dead[NUM_ITERATIONS][NUM_RANK] =
{ //nodes  edges  faces
   {0,     0,     0 }, //0
   {1,     3,     3 }, //1
   {3,     6,     5 }, //2
   {5,     10,    7 }, //3
   {7,     14,    9 }, //4
   {12,    20,    11}, //5
   {25,    34,    16}  //6
};

const int global_num_live[NUM_ITERATIONS][NUM_RANK] =
{ //nodes  edges  faces
   {25,    16,    16}, //0
   {24,    20,    13}, //1
   {22,    20,    11}, //2
   {20,    20,    9 }, //3
   {18,    18,    7 }, //4
   {13,    14,    5 }, //5
   {0,     0,     0 }  //6
};


typedef std::vector<stk::mesh::Entity *> EntityVector;

namespace {

//Generates a vector a entities to be killed in this iteration
EntityVector entities_to_be_kill( const stk::mesh::BulkData & mesh, int pass);

//Validates that the correct entites were killed in this iteration
bool validate_iteration( stk::ParallelMachine pm, GridFixture & fixture, int pass);

//Finds the sides that need to be created between the live and dead entities
void find_sides_to_be_created(
    const stk::mesh::EntitySideVector & boundary,
    const stk::mesh::Selector & select,
    std::vector<stk::mesh::EntitySideComponent> & sides
    );

//Finds entities from the closure of the entities_to_be_killed
//that only have relations with dead entities
void find_entities_to_kill(
    const EntityVector & entities_closure,
    unsigned closure_rank,
    unsigned entity_rank,
    const stk::mesh::Selector & select_owned,
    const stk::mesh::Selector & select_live,
    EntityVector & kill_list
    );

}

bool element_death_use_case_1(stk::ParallelMachine pm)
{
  //setup the mesh
  GridFixture fixture(pm);
  stk::mesh::BulkData& mesh = fixture.bulk_data();
  stk::mesh::MetaData& meta_data = fixture.meta_data();

  stk::mesh::Part & dead_part = *fixture.dead_part();

  stk::mesh::PartVector dead_parts;
  dead_parts.push_back( & dead_part);

  bool passed = true;

  unsigned volume_rank = stk::mesh::Face;

  for (int iteration = 0; iteration <NUM_ITERATIONS; ++iteration) {
    //find the entities to kill in this iteration
    EntityVector entities_to_kill = entities_to_be_kill(mesh, iteration);

    // find the parallel-consistent closure of the entities to be killed
    EntityVector entities_closure;
    stk::mesh::find_closure(mesh,
        entities_to_kill,
        entities_closure);


    // find the boundary of the entities we're killing
    stk::mesh::EntitySideVector boundary;
    stk::mesh::boundary_analysis(mesh,
        entities_closure,
        volume_rank,
        boundary);


    //find the sides that need to be created
    //sides need to be created when the outside
    //of the boundary is both live and owned and
    //a side seperating the live and dead doesn't
    //already exist
    stk::mesh::Selector select_owned = meta_data.locally_owned_part();
    stk::mesh::Selector select_live = ! dead_part ;
    stk::mesh::Selector select_live_and_owned = select_live & select_owned;

    std::vector<stk::mesh::EntitySideComponent> skin;
    find_sides_to_be_created( boundary, select_live_and_owned, skin);


    mesh.modification_begin();
    // Kill entities by moving them to the dead part.
    for (EntityVector::iterator itr = entities_to_kill.begin();
        itr != entities_to_kill.end(); ++itr) {
      mesh.change_entity_parts(**itr, dead_parts);
    }


    // Ask for new entites to represent the sides between the live and dead entities
    std::vector<size_t> requests(meta_data.entity_type_count(), 0);
    EntityVector requested_entities;
    requests[volume_rank-1] = skin.size();
    mesh.generate_new_entities(requests, requested_entities);

    // Create boundaries between live and dead entities
    for ( size_t i = 0; i < skin.size(); ++i) {
      stk::mesh::Entity & entity = *(skin[i].entity);
      const unsigned side_id  = skin[i].side_id;
      stk::mesh::Entity & side   = * (requested_entities[i]);

      stk::mesh::declare_element_side(entity, side, side_id);
    }
    mesh.modification_end();
    //the modification_end() will communicate which entities have been changed
    //to other processes.

    //find lower ranked entity that are only related to the dead entities
    //and kill them
    for (int rank = volume_rank -1; rank >= 0; --rank) {
      EntityVector kill_list;
      find_entities_to_kill( entities_closure, volume_rank, rank, select_owned, select_live, kill_list);

      //need to communicate killing the higher ranking entities among processors before killing the
      //lower.
      mesh.modification_begin();
      for (EntityVector::iterator itr = kill_list.begin();
          itr != kill_list.end(); ++itr) {
        mesh.change_entity_parts(**itr, dead_parts);
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
    const stk::mesh::EntitySideVector & boundary,
    const stk::mesh::Selector & select,
    std::vector<stk::mesh::EntitySideComponent> & sides
    )
{
  //look at the outside of the boundary since the inside will be kill this
  //iteration

  for (stk::mesh::EntitySideVector::const_iterator itr = boundary.begin();
      itr != boundary.end(); ++itr) {

    const stk::mesh::EntitySideComponent & outside = itr->outside;


    // examine the boundary of the outside of the closure.
    /// \TODO conside making a selector operate on entities as well as buckets
    if ( outside.entity != NULL && select(*(outside.entity)) ) {

      //make sure the side does not already exist
      const unsigned side_id = outside.side_id;
      const stk::mesh::Entity & entity = * outside.entity;
      stk::mesh::PairIterRelation existing_sides = entity.relations(entity.entity_rank()-1);

      for (; existing_sides.first != existing_sides.second &&
          existing_sides.first->identifier() != side_id ;
          ++existing_sides.first);

      //reached the end -- a new side needs to be created
      if (existing_sides.first == existing_sides.second) {
        sides.push_back(outside);
      }
    }
  }
}

//----------------------------------------------------------------------------------
void find_entities_to_kill(
    const EntityVector & entities_closure,
    unsigned volume_rank,
    unsigned entity_rank,
    const stk::mesh::Selector & select_owned,
    const stk::mesh::Selector & select_live,
    EntityVector & kill_list
    )
{

  kill_list.clear();

  //find the first entity in the closure
  EntityVector::const_iterator itr = std::lower_bound(entities_closure.begin(),
      entities_closure.end(),
      stk::mesh::EntityKey(entity_rank, 0),
      stk::mesh::EntityLess());

  const EntityVector::const_iterator end = std::lower_bound(entities_closure.begin(),
      entities_closure.end(),
      stk::mesh::EntityKey(entity_rank+1, 0),
      stk::mesh::EntityLess());

  for (; itr != end; ++itr) {
    stk::mesh::Entity & entity = **itr;

    if (select_owned(entity.bucket())) {
      bool found_live = false;

      for(unsigned rank = entity_rank + 1; rank<=volume_rank && !found_live; ++rank) {

        stk::mesh::PairIterRelation relations_pair = entity.relations(rank);

        for (; relations_pair.first != relations_pair.second && !found_live; ++relations_pair.first) {

          if( select_live(*relations_pair.first->entity())) {
            found_live = true;
          }
        }
      }

      if (!found_live) {
        kill_list.push_back(&entity);
      }
    }
  }

}


//----------------------------------------------------------------------------------

EntityVector entities_to_be_kill( const stk::mesh::BulkData & mesh, int pass) {

  std::vector<unsigned> entity_ids_to_kill;
  switch(pass) {
    case 0:
      break;
    case 1:
      entity_ids_to_kill.push_back(4);
      entity_ids_to_kill.push_back(9);
      entity_ids_to_kill.push_back(10);
      break;
    case 2:
      entity_ids_to_kill.push_back(2);
      entity_ids_to_kill.push_back(3);
      break;
    case 3:
      entity_ids_to_kill.push_back(1);
      entity_ids_to_kill.push_back(11);
      break;
    case 4:
      entity_ids_to_kill.push_back(6);
      entity_ids_to_kill.push_back(7);
      break;
    case 5:
      entity_ids_to_kill.push_back(5);
      entity_ids_to_kill.push_back(16);
      break;
    case 6:
      entity_ids_to_kill.push_back(8);
      entity_ids_to_kill.push_back(12);
      entity_ids_to_kill.push_back(13);
      entity_ids_to_kill.push_back(14);
      entity_ids_to_kill.push_back(15);
      break;
    default:
      break;
  }

  EntityVector entities_to_kill;
  for (std::vector<unsigned>::const_iterator itr = entity_ids_to_kill.begin();
      itr != entity_ids_to_kill.end(); ++itr) {
    stk::mesh::Entity * temp = mesh.get_entity(stk::mesh::Face, *itr);
    //select the entity only if the current process in the owner
    if (temp != NULL && temp->owner_rank() == mesh.parallel_rank()) {
      entities_to_kill.push_back(temp);
    }
  }
  return entities_to_kill;
}


//----------------------------------------------------------------------------------

bool validate_iteration( stk::ParallelMachine pm, GridFixture & fixture, int pass) {

  if (pass >= NUM_ITERATIONS || pass < 0) {
    return false;
  }

  stk::mesh::BulkData& mesh = fixture.bulk_data();
  stk::mesh::MetaData& meta_data = fixture.meta_data();

  stk::mesh::Part & dead_part = *fixture.dead_part();

  stk::mesh::Selector select_dead = dead_part & meta_data.locally_owned_part();
  stk::mesh::Selector select_live = !dead_part & meta_data.locally_owned_part();

  int num_dead[NUM_RANK] = {0, 0, 0};
  int num_live[NUM_RANK] = {0, 0, 0};

  for ( int i = 0; i < NUM_RANK ; ++i) {
    const std::vector<stk::mesh::Bucket*>& buckets = mesh.buckets( stk::mesh::fem_entity_type(i));
    num_dead[i] = count_selected_entities( select_dead, buckets);
    num_live[i] = count_selected_entities( select_live, buckets);
  }

  stk::all_reduce(pm, stk::ReduceSum<3>(num_dead) & stk::ReduceSum<3>(num_live));

  bool correct_dead = true;
  bool correct_live = true;

  for (int i=0; i<NUM_RANK; ++i) {
    correct_dead &= global_num_dead[pass][i] == num_dead[i];
    correct_live &= global_num_live[pass][i] == num_live[i];
  }

  return correct_dead && correct_live;
}

}

