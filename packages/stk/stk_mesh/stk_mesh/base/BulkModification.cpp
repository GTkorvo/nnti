#include <set>
#include <stdexcept>
#include <sstream>

#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_util/parallel/ParallelComm.hpp>


namespace stk {
namespace mesh {

typedef std::vector<Entity *> EntityVector;
typedef std::set<Entity *, EntityLess> EntitySet;
typedef std::set<EntityProc , EntityLess> EntityProcSet;

namespace {

void construct_transitive_closure( EntitySet & closure , Entity & entry )
{

  std::pair< EntitySet::const_iterator , bool >
    result = closure.insert( & entry );

  // A new insertion, must also insert the closure
  if ( result.second ) {

    const unsigned etype = entry.entity_type();
    PairIterRelation irel  = entry.relations();

    for ( ; irel.first != irel.second ; ++irel.first ) {
      // insert entities with relations of lower rank into the closure
      if ( irel.first->entity_rank() < etype ) {
        Entity * tmp = irel.first->entity();
        construct_transitive_closure( closure , *tmp );
      }
    }
  }
}

void find_local_closure ( EntitySet & closure, const EntityVector & entities)
{
  for (EntityVector::const_iterator i = entities.begin();
      i != entities.end(); ++i)
  {
    construct_transitive_closure(closure, **i);
  }
}

void construct_communication_set( const BulkData & bulk, const EntitySet & closure, EntityProcSet & communication_set)
{
  if (bulk.parallel_size() < 2) return;

  for ( EntitySet::const_iterator i = closure.begin();
      i != closure.end(); ++i)
  {
    Entity & entity = **i;

    //add sharing processes to communication_set
    PairIterEntityProc sharing = entity.sharing();
    communication_set.insert(sharing.first, sharing.second);

    const std::vector<Ghosting *> & all_ghosting = bulk.ghostings();

    //add Ghosting to the communication_set
    for( std::vector<Ghosting *>::const_iterator g = all_ghosting.begin();
        g != all_ghosting.end(); ++g)
    {
      const Ghosting & ghost = **g;
      const std::vector< EntityProc > & ghost_send = ghost.send();

      //ghost_send is sorted and unique
      //find where this entity is in the Ghosting list
      std::vector< EntityProc >::const_iterator j_begin =
        std::lower_bound(ghost_send.begin(),ghost_send.end(), entity, EntityLess());

      //find the end of this entity in the Ghosting list
      std::vector< EntityProc >::const_iterator j_end = j_begin;
      for(; j_end != ghost_send.end() && j_end->first == &entity; ++j_end);

      //insert this entity Ghosting into the communication_set
      communication_set.insert(j_begin,j_end);
    }

  }
}

size_t count_non_used_entities( const BulkData & bulk, const EntityVector & entities)
{

  size_t non_used_entities = 0;

  const Part & locally_used_part = bulk.mesh_meta_data().locally_used_part();

  //Check that entities are only in the locally_used part
  for (EntityVector::const_iterator i = entities.begin();
      i != entities.end(); ++i)
  {
    const Bucket & b = (**i).bucket();
    if ( ! has_superset(b, locally_used_part)) {
      ++non_used_entities;
    }
  }

  return non_used_entities;

}

}



void find_closure( const BulkData & bulk,
    const std::vector< Entity *> & entities,
    std::vector< Entity *> & entities_closure)
{

  entities_closure.clear();


  EntityProcSet send_list;
  EntitySet     temp_entities_closure;

  const bool bulk_not_synchronized = bulk.synchronized_state() != BulkData::SYNCHRONIZED;
  const size_t non_used_entities = bulk_not_synchronized ? 0 : count_non_used_entities(bulk, entities);

  const bool local_bad_input = bulk_not_synchronized || (0 < non_used_entities);

  //Can skip if error on input
  if ( !local_bad_input) {

    find_local_closure(temp_entities_closure, entities);

    construct_communication_set(bulk, temp_entities_closure, send_list);
  }


  CommAll all( bulk.parallel() );

  //pack send_list for sizing
  for ( EntityProcSet::const_iterator
      ep = send_list.begin() ; ep != send_list.end() ; ++ep ) {
    all.send_buffer( ep->second).pack<EntityKey>(ep->first->key());
  }


  const bool global_bad_input = all.allocate_buffers( bulk.parallel_size() / 4 , false, local_bad_input );

  if (global_bad_input) {

    std::ostringstream msg;
    //parallel consisent throw
    if (bulk_not_synchronized) {
      msg << "stk::mesh::find_closure( const BulkData & bulk, ... ) bulk is not synchronized";
    }
    else if ( 0 < non_used_entities) {
      msg << "stk::mesh::find_closure( const BulkData & bulk, std::vector<Entity *> entities, ... ) \n"
          << "entities contains " << non_used_entities << " non locally used entities \n";
    }

    throw std::runtime_error(msg.str());
  }


  //pack send_list
  for ( EntityProcSet::const_iterator
      ep = send_list.begin() ; ep != send_list.end() ; ++ep ) {
    all.send_buffer( ep->second).pack<EntityKey>(ep->first->key());
  }


  all.communicate();

  //unpack the send_list into the temp entities closure set
  for ( unsigned p = 0 ; p < bulk.parallel_size() ; ++p ) {
    CommBuffer & buf = all.recv_buffer( p );
    EntityKey k ;
    while ( buf.remaining() ) {
      buf.unpack<EntityKey>( k );
      Entity * e = bulk.get_entity(k);
      temp_entities_closure.insert(e);
    }
  }

  //copy the set into the entities_closure vector
  entities_closure.assign(temp_entities_closure.begin(), temp_entities_closure.end());
}

} // namespace mesh
} // namespace stk

